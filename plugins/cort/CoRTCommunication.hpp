#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "m2n/DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace m2n {

namespace cort {
class DoubleDataExchange;
class NodeTopology;
class Transport;
} // namespace cort

/// CoRT communication implementation of DistributedCommunication.
class CoRTCommunication : public DistributedCommunication {
public:
  CoRTCommunication(com::PtrCommunicationFactory communicationFactory,
                    mesh::PtrMesh                mesh);

  ~CoRTCommunication() override;

  /// Returns true, if a connection to a remote participant has been established.
  bool isConnected() const override;

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param[in] acceptorName  Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptConnection(std::string const &acceptorName,
                        std::string const &requesterName) override;

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestConnection(std::string const &acceptorName,
                         std::string const &requesterName) override;

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestPreConnection().
   *        Only initial connection is created.
   *
   * @param[in] acceptorName  Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptPreConnection(std::string const &acceptorName,
                           std::string const &requesterName) override;

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *        Only initial connection is created.
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestPreConnection(std::string const &acceptorName,
                            std::string const &requesterName) override;

  /// Completes the secondary connections for both acceptor and requester by updating the vertex list in _mappings
  void completeSecondaryRanksConnection() override;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection() override;

  /**
   * @brief Sends a subset of local double values corresponding to local indices
   *        deduced from the current and remote vertex distributions.
   */
  void send(precice::span<double const> itemsToSend, int valueDimension = 1) override;

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  void receive(precice::span<double> itemsToReceive, int valueDimension = 1) override;

  /// Broadcasts an int to connected ranks on remote participant
  void broadcastSend(int itemToSend) override;

  /**
   * @brief Receives an int per connected rank on remote participant
   * @para[out] itemToReceive received ints from remote ranks are stored with the sender rank order
   */
  void broadcastReceiveAll(std::vector<int> &itemToReceive) override;

  /// Broadcasts a mesh to connected ranks on remote participant
  void broadcastSendMesh() override;

  /// Receive mesh partitions per connected rank on remote participant
  void broadcastReceiveAllMesh() override;

  /// Scatters a communication map over connected ranks on remote participant
  void scatterAllCommunicationMap(CommunicationMap &localCommunicationMap) override;

  /// Gathers a communication maps from connected ranks on remote participant
  void gatherAllCommunicationMap(CommunicationMap &localCommunicationMap) override;

private:
  logging::Logger _log{"m2n::CoRTCommunication"};

  com::PtrCommunicationFactory _communicationFactory;

  /// Communication used for the active direct-control or proxy-data channel.
  com::PtrCommunication _communication;

  /**
   * @brief Defines mapping between:
   *        1. global remote process rank;
   *        2. local data indices, which define a subset of local (for process
   *           rank in the current participant) data to be communicated between
   *           the current process rank and the remote process rank.
   */
  struct Mapping {
    int              remoteRank;
    std::vector<int> indices;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        mappings (one to service each point-to-point connection).
   */
  std::vector<Mapping> _mappings;

  /**
   * @brief this data structure is used to store m2n communication information for the 1 step of
   *        bounding box initialization. It stores:
   *        1. global remote process rank;
   */
  struct ConnectionData {
    int remoteRank;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        ConnectionData (one to service each point-to-point connection).
   */
  std::vector<ConnectionData> _connectionDataVector;

  bool _isConnected       = false;
  bool _usesPreConnection = false;

  enum class ConnectionRole {
    Undefined,
    Acceptor,
    Requester
  };

  ConnectionRole _connectionRole = ConnectionRole::Undefined;
  std::string    _acceptorName;
  std::string    _requesterName;

  /// Routing table: remote rank -> remote proxy rank.
  std::map<int, int> _remoteRankToProxy;

  std::unique_ptr<cort::NodeTopology>       _nodeTopology;
  std::unique_ptr<cort::DoubleDataExchange> _dataExchange;
  std::shared_ptr<cort::Transport>          _transport;

  void          setupLocalTopology();
  void          configureDataExchange();
  void          initializeDefaultDataPatterns();
  void          exchangeProxyMaps(std::string const &acceptorName,
                                  std::string const &requesterName,
                                  bool               isAcceptor);
  void          fillRemoteProxyMap(const std::vector<int> &remoteProxyMap);
  void          storeConnectionData(const std::vector<int> &connectedRanks);
  void          createDirectControlCommunication(std::string const &acceptorName,
                                                 std::string const &requesterName,
                                                 bool               isAcceptor);
  void          closeDirectControlCommunication();
  void          createProxyDataCommunication(std::string const &acceptorName,
                                             std::string const &requesterName,
                                             bool               isAcceptor);
  std::set<int> gatherUniqueRemoteProxies() const;
};
} // namespace m2n
} // namespace precice
