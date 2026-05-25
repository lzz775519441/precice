#include "CoRTCommunication.hpp"
#include <algorithm>
#include <boost/io/ios_state.hpp>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "com/Extra.hpp"
#include "com/Request.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "cort/DoubleDataExchange.hpp"
#include "cort/NodeTopology.hpp"
#include "cort/Transport.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/algorithm.hpp"
#include "utils/assertion.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif

using precice::profiling::Event;

namespace precice::m2n {

namespace {

constexpr int DefaultValueDimension = 1;

class ComRequestAdapter : public cort::Request {
public:
  explicit ComRequestAdapter(com::PtrRequest request)
      : _request(std::move(request))
  {
  }

  void wait() override
  {
    _request->wait();
  }

private:
  com::PtrRequest _request;
};

class ComTransportAdapter : public cort::Transport {
public:
  explicit ComTransportAdapter(com::PtrCommunication communication)
      : _communication(std::move(communication))
  {
  }

  cort::PtrRequest asyncSend(precice::span<double const> itemsToSend, int remoteRank) override
  {
    return std::make_shared<ComRequestAdapter>(_communication->aSend(itemsToSend, remoteRank));
  }

  cort::PtrRequest asyncReceive(precice::span<double> itemsToReceive, int remoteRank) override
  {
    return std::make_shared<ComRequestAdapter>(_communication->aReceive(itemsToReceive, remoteRank));
  }

private:
  com::PtrCommunication _communication;
};

} // namespace

namespace impl {
void send(mesh::Mesh::VertexDistribution const &m,
          int                                   rankReceiver,
          const com::PtrCommunication          &communication)
{
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->send(rank, rankReceiver);
    communication->sendRange(indices, rankReceiver);
  }
}

void receive(mesh::Mesh::VertexDistribution &m,
             int                             rankSender,
             const com::PtrCommunication    &communication)
{
  m.clear();
  int size = 0;
  communication->receive(size, rankSender);

  while (size--) {
    Rank rank = -1;
    communication->receive(rank, rankSender);
    m[rank] = communication->receiveRange(rankSender, com::asVector<int>);
  }
}

void broadcastSend(mesh::Mesh::VertexDistribution const &m,
                   const com::PtrCommunication          &communication = utils::IntraComm::getCommunication())
{
  communication->broadcast(static_cast<int>(m.size()));

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->broadcast(rank);
    communication->broadcast(indices);
  }
}

void broadcastReceive(mesh::Mesh::VertexDistribution &m,
                      int                             rankBroadcaster,
                      const com::PtrCommunication    &communication = utils::IntraComm::getCommunication())
{
  m.clear();
  int size = 0;
  communication->broadcast(size, rankBroadcaster);

  while (size--) {
    Rank rank = -1;
    communication->broadcast(rank, rankBroadcaster);
    communication->broadcast(m[rank], rankBroadcaster);
  }
}

void broadcast(mesh::Mesh::VertexDistribution &m)
{
  if (utils::IntraComm::isPrimary()) {
    m2n::impl::broadcastSend(m);
  } else if (utils::IntraComm::isSecondary()) {
    m2n::impl::broadcastReceive(m, 0);
  }
}

void print(std::map<int, std::vector<int>> const &m)
{
  std::ostringstream oss;
  oss << "rank: " << utils::IntraComm::getRank() << "\n";
  for (auto &i : m) {
    for (auto &j : i.second) {
      oss << i.first << ":" << j << '\n';
    }
  }

  if (utils::IntraComm::isSecondary()) {
    utils::IntraComm::getCommunication()->send(oss.str(), 0);
  } else {
    std::string s;
    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(s, rank);
      oss << s;
    }
    std::cout << oss.str();
  }
}

void printCommunicationPartnerCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = m.size();
  if (utils::IntraComm::isPrimary()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));
      count++;
    }

    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(size, rank);
      total += size;
      if (size) {
        maximum = std::max(maximum, static_cast<size_t>(size));
        minimum = std::min(minimum, static_cast<size_t>(size));
        count++;
      }
    }

    if (minimum > maximum)
      minimum = maximum;

    auto average = static_cast<double>(total);
    if (count != 0)
      average /= count;

    boost::io::ios_all_saver ias{std::cout};
    std::cout << std::fixed << std::setprecision(3)
              << "Number of Communication Partners per Interface Process:\n"
              << "  Total:   " << total << "\n"
              << "  Maximum: " << maximum << "\n"
              << "  Minimum: " << minimum << "\n"
              << "  Average: " << average << "\n"
              << "Number of Interface Processes: " << count << "\n\n";
  } else {
    PRECICE_ASSERT(utils::IntraComm::isSecondary());
    utils::IntraComm::getCommunication()->send(size, 0);
  }
}

void printLocalIndexCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = 0;
  for (auto &i : m)
    size += i.second.size();

  if (utils::IntraComm::isPrimary()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));
      count++;
    }

    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(size, rank);
      total += size;
      if (size) {
        maximum = std::max(maximum, static_cast<size_t>(size));
        minimum = std::min(minimum, static_cast<size_t>(size));
        count++;
      }
    }

    if (minimum > maximum)
      minimum = maximum;

    auto average = static_cast<double>(total);
    if (count != 0)
      average /= count;

    boost::io::ios_all_saver ias{std::cout};
    std::cout << std::fixed << std::setprecision(3)
              << "Number of LVDIs per Interface Process:\n"
              << "  Total:   " << total << '\n'
              << "  Maximum: " << maximum << '\n'
              << "  Minimum: " << minimum << '\n'
              << "  Average: " << average << '\n'
              << "Number of Interface Processes: " << count << "\n\n";
  } else {
    PRECICE_ASSERT(utils::IntraComm::isSecondary());
    utils::IntraComm::getCommunication()->send(size, 0);
  }
}

std::map<int, std::vector<int>> buildCommunicationMap(
    mesh::Mesh::VertexDistribution const &thisVertexDistribution,
    mesh::Mesh::VertexDistribution const &otherVertexDistribution,
    int                                   thisRank = utils::IntraComm::getRank())
{
  auto iterator = thisVertexDistribution.find(thisRank);
  if (iterator == thisVertexDistribution.end()) {
    return {};
  }

  std::map<int, std::vector<int>> communicationMap;
  PRECICE_ASSERT(std::is_sorted(iterator->second.begin(), iterator->second.end()));

  for (const auto &[rank, vertices] : otherVertexDistribution) {
    PRECICE_ASSERT(std::is_sorted(vertices.begin(), vertices.end()));

    if (iterator->second.empty() || vertices.empty() || (vertices.back() < iterator->second.at(0)) || (vertices.at(0) > iterator->second.back())) {
      continue;
    }
    std::vector<int> inters;
    precice::utils::set_intersection_indices(iterator->second.begin(), iterator->second.begin(), iterator->second.end(),
                                             vertices.begin(), vertices.end(),
                                             std::back_inserter(inters));
    if (!inters.empty()) {
      communicationMap.insert({rank, std::move(inters)});
    }
  }
  return communicationMap;
}

} // namespace impl

using namespace impl;

CoRTCommunication::CoRTCommunication(com::PtrCommunicationFactory communicationFactory,
                                     mesh::PtrMesh                mesh)
    : DistributedCommunication(mesh),
      _communicationFactory(std::move(communicationFactory)),
      _nodeTopology(std::make_unique<cort::NodeTopology>()),
      _dataExchange(std::make_unique<cort::DoubleDataExchange>(*_nodeTopology))
{
}

CoRTCommunication::~CoRTCommunication()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

bool CoRTCommunication::isConnected() const
{
  return _isConnected;
}

void CoRTCommunication::setupLocalTopology()
{
  _nodeTopology->initialize();
}

std::set<int> CoRTCommunication::gatherUniqueRemoteProxies() const
{
#ifdef PRECICE_NO_MPI
  throw std::runtime_error{"CoRT communication requires preCICE to be compiled with MPI support."};
#else
  std::set<int> localProxies;
  for (const auto &mapping : _mappings) {
    auto proxy = _remoteRankToProxy.find(mapping.remoteRank);
    if (!mapping.indices.empty() && proxy != _remoteRankToProxy.end()) {
      localProxies.insert(proxy->second);
    }
  }

  std::vector<int> localProxyVector(localProxies.begin(), localProxies.end());
  int              localCount = static_cast<int>(localProxyVector.size());
  std::vector<int> counts(_nodeTopology->localSize());
  MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, _nodeTopology->localCommunicator());

  std::vector<int> displs(_nodeTopology->localSize() + 1, 0);
  for (int i = 0; i < _nodeTopology->localSize(); ++i) {
    displs[i + 1] = displs[i] + counts[i];
  }

  std::vector<int> gathered(displs.back());
  MPI_Allgatherv(localProxyVector.data(), localCount, MPI_INT,
                 gathered.data(), counts.data(), displs.data(), MPI_INT,
                 _nodeTopology->localCommunicator());

  return std::set<int>(gathered.begin(), gathered.end());
#endif
}

void CoRTCommunication::configureDataExchange()
{
  std::vector<cort::RankMapping> rankMappings;
  rankMappings.reserve(_mappings.size());
  for (const auto &mapping : _mappings) {
    rankMappings.push_back({mapping.remoteRank, mapping.indices});
  }
  _dataExchange->configure(std::move(rankMappings), _remoteRankToProxy, _transport);
}

void CoRTCommunication::initializeDefaultDataPatterns()
{
  _dataExchange->initializePatterns(DefaultValueDimension);
}

void CoRTCommunication::fillRemoteProxyMap(const std::vector<int> &remoteProxyMap)
{
  _remoteRankToProxy.clear();
  for (int rank = 0; rank < static_cast<int>(remoteProxyMap.size()); ++rank) {
    _remoteRankToProxy[rank] = remoteProxyMap[rank];
  }
}

void CoRTCommunication::exchangeProxyMaps(std::string const &acceptorName,
                                          std::string const &requesterName,
                                          bool               isAcceptor)
{
  setupLocalTopology();

  const std::vector<int> &localProxyMap = _nodeTopology->proxyRankByRank();
  std::vector<int>        remoteProxyMap;

  if (not utils::IntraComm::isSecondary()) {
    Event             e("m2n.exchangeProxyMaps");
    auto              communication = _communicationFactory->newCommunication();
    const std::string tag           = "TMP-CORT-PROXYMAP-" + _mesh->getName();

    if (isAcceptor) {
      communication->acceptConnection(acceptorName, requesterName, tag, utils::IntraComm::getRank());

      const int localSize = static_cast<int>(localProxyMap.size());
      communication->send(localSize, 0);
      communication->send(precice::span<const int>(localProxyMap.data(), localSize), 0);

      int remoteSize = 0;
      communication->receive(remoteSize, 0);
      remoteProxyMap.resize(remoteSize);
      communication->receive(precice::span<int>(remoteProxyMap.data(), remoteSize), 0);
    } else {
      communication->requestConnection(acceptorName, requesterName, tag, 0, 1);

      int remoteSize = 0;
      communication->receive(remoteSize, 0);
      remoteProxyMap.resize(remoteSize);
      communication->receive(precice::span<int>(remoteProxyMap.data(), remoteSize), 0);

      const int localSize = static_cast<int>(localProxyMap.size());
      communication->send(localSize, 0);
      communication->send(precice::span<const int>(localProxyMap.data(), localSize), 0);
    }

    communication->closeConnection();
  }

  _nodeTopology->broadcastFromPrimary(remoteProxyMap);
  fillRemoteProxyMap(remoteProxyMap);
}

void CoRTCommunication::storeConnectionData(const std::vector<int> &connectedRanks)
{
  _connectionDataVector.clear();
  _connectionDataVector.reserve(connectedRanks.size());
  for (int connectedRank : connectedRanks) {
    _connectionDataVector.push_back({connectedRank});
  }
}

void CoRTCommunication::createDirectControlCommunication(std::string const &acceptorName,
                                                         std::string const &requesterName,
                                                         bool               isAcceptor)
{
  const std::vector<int> connectedRanks = _mesh->getConnectedRanks();
  storeConnectionData(connectedRanks);

  if (connectedRanks.empty()) {
    return;
  }

  _communication = _communicationFactory->newCommunication();
  if (isAcceptor) {
    _communication->acceptConnectionAsServer(acceptorName,
                                             requesterName,
                                             _mesh->getName(),
                                             utils::IntraComm::getRank(),
                                             connectedRanks.size());
  } else {
    std::set<int> acceptingRanks(connectedRanks.begin(), connectedRanks.end());
    _communication->requestConnectionAsClient(acceptorName,
                                              requesterName,
                                              _mesh->getName(),
                                              acceptingRanks,
                                              utils::IntraComm::getRank());
  }
}

void CoRTCommunication::closeDirectControlCommunication()
{
  if (_communication) {
    _communication->closeConnection();
    _communication.reset();
  }
  _transport.reset();
}

void CoRTCommunication::createProxyDataCommunication(std::string const &acceptorName,
                                                     std::string const &requesterName,
                                                     bool               isAcceptor)
{
  std::set<int> uniqueRemoteProxies = gatherUniqueRemoteProxies();

  if (!_nodeTopology->isProxy() || uniqueRemoteProxies.empty()) {
    return;
  }

  const std::string tag = _usesPreConnection ? "CORT-DATA-" + _mesh->getName() : _mesh->getName();

  _communication = _communicationFactory->newCommunication();
  if (isAcceptor) {
    _communication->acceptConnectionAsServer(acceptorName,
                                             requesterName,
                                             tag,
                                             utils::IntraComm::getRank(),
                                             uniqueRemoteProxies.size());
  } else {
    _communication->requestConnectionAsClient(acceptorName,
                                              requesterName,
                                              tag,
                                              uniqueRemoteProxies,
                                              utils::IntraComm::getRank());
  }
  _transport = std::make_shared<ComTransportAdapter>(_communication);
}

void CoRTCommunication::acceptConnection(std::string const &acceptorName,
                                         std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  setupLocalTopology();
  mesh::Mesh::VertexDistribution vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution requesterVertexDistribution;

  const std::vector<int> &localProxyMap = _nodeTopology->proxyRankByRank();
  std::vector<int>        remoteProxyMap;
  PRECICE_DEBUG("Exchange vertex distribution");
  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Exchange vertex distribution & Topology");
    Event e0("m2n.exchangeVertexDistribution");
    auto  c = _communicationFactory->newCommunication();

    c->acceptConnection(acceptorName, requesterName, "TMP-PRIMARYCOM-" + _mesh->getName(), utils::IntraComm::getRank());

    m2n::send(vertexDistribution, 0, c);

    int lSize = localProxyMap.size();
    c->send(lSize, 0);
    c->send(precice::span<const int>(localProxyMap.data(), lSize), 0);

    m2n::receive(requesterVertexDistribution, 0, c);

    int rSize = 0;
    c->receive(rSize, 0);
    remoteProxyMap.resize(rSize);
    c->receive(precice::span<int>(remoteProxyMap.data(), rSize), 0);
  }

  PRECICE_DEBUG("broadcastVertexDistributions");
  Event e1("m2n.broadcastVertexDistributions", profiling::Synchronize);
  m2n::broadcast(vertexDistribution);
  if (utils::IntraComm::isSecondary()) {
    _mesh->setVertexDistribution(vertexDistribution);
  }
  m2n::broadcast(requesterVertexDistribution);

  _nodeTopology->broadcastFromPrimary(remoteProxyMap);
  fillRemoteProxyMap(remoteProxyMap);
  e1.stop();

  PRECICE_DEBUG("buildCommunicationMap");
  Event                           e2("m2n.buildCommunicationMap", profiling::Synchronize);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      vertexDistribution, requesterVertexDistribution);
  e2.stop();

  Event e4("m2n.createCommunications");
  e4.addData("Connections", communicationMap.size());

  for (auto const &comMap : communicationMap) {
    int  globalRequesterRank = comMap.first;
    auto indices             = comMap.second;
    _mappings.push_back({globalRequesterRank, std::move(indices)});
    _connectionDataVector.push_back({globalRequesterRank});
  }

  createProxyDataCommunication(acceptorName, requesterName, true);
  e4.stop();
  _isConnected = true;
  configureDataExchange();
  initializeDefaultDataPatterns();
}

void CoRTCommunication::requestConnection(std::string const &acceptorName,
                                          std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  setupLocalTopology();

  mesh::Mesh::VertexDistribution vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution acceptorVertexDistribution;

  const std::vector<int> &localProxyMap = _nodeTopology->proxyRankByRank();
  std::vector<int>        remoteProxyMap;

  if (not utils::IntraComm::isSecondary()) {
    Event e0("m2n.exchangeVertexDistribution");
    auto  c = _communicationFactory->newCommunication();
    c->requestConnection(acceptorName, requesterName, "TMP-PRIMARYCOM-" + _mesh->getName(), 0, 1);

    m2n::receive(acceptorVertexDistribution, 0, c);

    int rSize = 0;
    c->receive(rSize, 0);
    remoteProxyMap.resize(rSize);
    c->receive(precice::span<int>(remoteProxyMap.data(), rSize), 0);

    m2n::send(vertexDistribution, 0, c);
    int lSize = localProxyMap.size();
    c->send(lSize, 0);
    c->send(precice::span<const int>(localProxyMap.data(), lSize), 0);
  }

  Event e1("m2n.broadcastVertexDistributions", profiling::Synchronize);
  m2n::broadcast(vertexDistribution);
  if (utils::IntraComm::isSecondary()) {
    _mesh->setVertexDistribution(vertexDistribution);
  }
  m2n::broadcast(acceptorVertexDistribution);

  _nodeTopology->broadcastFromPrimary(remoteProxyMap);
  fillRemoteProxyMap(remoteProxyMap);
  e1.stop();

  Event                           e2("m2n.buildCommunicationMap", profiling::Synchronize);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      vertexDistribution, acceptorVertexDistribution);
  e2.stop();

  Event e4("m2n.createCommunications");
  e4.addData("Connections", communicationMap.size());

  for (auto const &comMap : communicationMap) {
    auto globalAcceptorRank = comMap.first;
    auto indices            = comMap.second;
    _mappings.push_back({globalAcceptorRank, std::move(indices)});
    _connectionDataVector.push_back({globalAcceptorRank});
  }

  createProxyDataCommunication(acceptorName, requesterName, false);
  e4.stop();
  _isConnected = true;
  configureDataExchange();
  initializeDefaultDataPatterns();
}

void CoRTCommunication::acceptPreConnection(std::string const &acceptorName,
                                            std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  _connectionRole    = ConnectionRole::Acceptor;
  _acceptorName      = acceptorName;
  _requesterName     = requesterName;
  _usesPreConnection = true;

  exchangeProxyMaps(acceptorName, requesterName, true);
  createDirectControlCommunication(acceptorName, requesterName, true);
  _isConnected = true;
}

void CoRTCommunication::requestPreConnection(std::string const &acceptorName,
                                             std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  _connectionRole    = ConnectionRole::Requester;
  _acceptorName      = acceptorName;
  _requesterName     = requesterName;
  _usesPreConnection = true;

  exchangeProxyMaps(acceptorName, requesterName, false);
  createDirectControlCommunication(acceptorName, requesterName, false);
  _isConnected = true;
}

void CoRTCommunication::completeSecondaryRanksConnection()
{
  PRECICE_ASSERT(_usesPreConnection, "CoRT can only complete a connection created by acceptPreConnection() or requestPreConnection().");

  mesh::Mesh::CommunicationMap localCommunicationMap = _mesh->getCommunicationMap();
  _mappings.clear();
  for (auto &i : _connectionDataVector) {
    _mappings.push_back({i.remoteRank, std::move(localCommunicationMap[i.remoteRank])});
  }

  closeDirectControlCommunication();
  createProxyDataCommunication(_acceptorName, _requesterName, _connectionRole == ConnectionRole::Acceptor);
  _usesPreConnection = false;

  configureDataExchange();
  initializeDefaultDataPatterns();
}

void CoRTCommunication::closeConnection()
{
  PRECICE_TRACE();
  if (not isConnected())
    return;

  _dataExchange->reset();

  if (_communication) {
    _communication->closeConnection();
  }
  _transport.reset();
  _communication.reset();
  _mappings.clear();
  _connectionDataVector.clear();
  _remoteRankToProxy.clear();
  _usesPreConnection = false;
  _connectionRole    = ConnectionRole::Undefined;
  _acceptorName.clear();
  _requesterName.clear();
  _isConnected = false;
  _nodeTopology->reset();
}

void CoRTCommunication::send(precice::span<double const> itemsToSend, int valueDimension)
{
  _dataExchange->send(itemsToSend, valueDimension);
}

void CoRTCommunication::receive(precice::span<double> itemsToReceive, int valueDimension)
{
  _dataExchange->receive(itemsToReceive, valueDimension);
}

void CoRTCommunication::broadcastSend(int itemToSend)
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      _communication->send(itemToSend, connectionData.remoteRank);
    }
    return;
  }

  if (!_nodeTopology->isProxy() || !_communication)
    return;

  std::set<int> sentProxies;
  for (auto &connectionData : _connectionDataVector) {
    if (_remoteRankToProxy.count(connectionData.remoteRank)) {
      int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
      if (sentProxies.find(targetProxy) == sentProxies.end()) {
        _communication->send(itemToSend, targetProxy);
        sentProxies.insert(targetProxy);
      }
    }
  }
}

void CoRTCommunication::broadcastReceiveAll(std::vector<int> &itemToReceive)
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      int data = 0;
      _communication->receive(data, connectionData.remoteRank);
      itemToReceive.push_back(data);
    }
    return;
  }

  std::map<int, int> proxyData;

  if (_nodeTopology->isProxy() && _communication) {
    std::set<int> recvProxies;
    for (auto &connectionData : _connectionDataVector) {
      if (_remoteRankToProxy.count(connectionData.remoteRank)) {
        int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
        if (recvProxies.find(targetProxy) == recvProxies.end()) {
          int val = 0;
          _communication->receive(val, targetProxy);
          proxyData[targetProxy] = val;
          recvProxies.insert(targetProxy);
        }
      }
    }
  }

  for (auto &connectionData : _connectionDataVector) {
    int data = 0;
    if (_nodeTopology->isProxy() && _remoteRankToProxy.count(connectionData.remoteRank)) {
      data = proxyData[_remoteRankToProxy[connectionData.remoteRank]];
    }
    _nodeTopology->broadcastFromLocalProxy(data);
    itemToReceive.push_back(data);
  }
}

void CoRTCommunication::broadcastSendMesh()
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      com::sendMesh(*_communication, connectionData.remoteRank, *_mesh);
    }
    return;
  }

  if (!_nodeTopology->isProxy() || !_communication)
    return;

  std::set<int> sentProxies;
  for (auto &connectionData : _connectionDataVector) {
    if (_remoteRankToProxy.count(connectionData.remoteRank)) {
      int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
      if (sentProxies.find(targetProxy) == sentProxies.end()) {
        com::sendMesh(*_communication, targetProxy, *_mesh);
        sentProxies.insert(targetProxy);
      }
    }
  }
}

void CoRTCommunication::broadcastReceiveAllMesh()
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      com::receiveMesh(*_communication, connectionData.remoteRank, *_mesh);
    }
    return;
  }

  if (_nodeTopology->isProxy() && _communication) {
    std::set<int> recvProxies;
    for (auto &connectionData : _connectionDataVector) {
      if (_remoteRankToProxy.count(connectionData.remoteRank)) {
        int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
        if (recvProxies.find(targetProxy) == recvProxies.end()) {
          com::receiveMesh(*_communication, targetProxy, *_mesh);
          recvProxies.insert(targetProxy);
        }
      }
    }
  }
}

void CoRTCommunication::scatterAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      _communication->sendRange(localCommunicationMap[connectionData.remoteRank], connectionData.remoteRank);
    }
    return;
  }

  if (!_nodeTopology->isProxy() || !_communication)
    return;

  for (auto &connectionData : _connectionDataVector) {
    if (_remoteRankToProxy.count(connectionData.remoteRank)) {
      int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
      _communication->sendRange(localCommunicationMap[connectionData.remoteRank], targetProxy);
    }
  }
}

void CoRTCommunication::gatherAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  if (_usesPreConnection) {
    if (!_communication)
      return;

    for (auto &connectionData : _connectionDataVector) {
      localCommunicationMap[connectionData.remoteRank] = _communication->receiveRange(connectionData.remoteRank, com::asVector<int>);
    }
    return;
  }

  for (auto &connectionData : _connectionDataVector) {
    std::vector<int> recvVec;

    if (_nodeTopology->isProxy() && _communication) {
      if (_remoteRankToProxy.count(connectionData.remoteRank)) {
        int targetProxy = _remoteRankToProxy[connectionData.remoteRank];
        recvVec         = _communication->receiveRange(targetProxy, com::asVector<int>);
      }
    }

    _nodeTopology->broadcastFromLocalProxy(recvVec);

    localCommunicationMap[connectionData.remoteRank] = recvVec;
  }
}
} // namespace precice::m2n
