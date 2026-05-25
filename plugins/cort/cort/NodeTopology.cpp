#include "cort/NodeTopology.hpp"

#include <stdexcept>
#include <utility>

#include "logging/LogMacros.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif

namespace precice::m2n::cort {

NodeTopology::~NodeTopology()
{
  reset();
}

NodeTopology::NodeTopology(NodeTopology &&other) noexcept
{
  *this = std::move(other);
}

NodeTopology &NodeTopology::operator=(NodeTopology &&other) noexcept
{
  if (this == &other) {
    return *this;
  }

  reset();

  _initialized     = other._initialized;
  _isProxy         = other._isProxy;
  _localRank       = other._localRank;
  _localSize       = other._localSize;
  _proxyRank       = other._proxyRank;
  _proxyRankByRank = std::move(other._proxyRankByRank);

#ifndef PRECICE_NO_MPI
  _participantComm       = other._participantComm;
  _localComm             = other._localComm;
  other._participantComm = MPI_COMM_NULL;
  other._localComm       = MPI_COMM_NULL;
#endif

  other._initialized = false;
  other._isProxy     = false;
  other._localRank   = 0;
  other._localSize   = 1;
  other._proxyRank   = 0;

  return *this;
}

void NodeTopology::initialize()
{
  if (_initialized) {
    return;
  }

#ifdef PRECICE_NO_MPI
  throw std::runtime_error{"CoRT communication requires preCICE to be compiled with MPI support."};
#else
  _participantComm = utils::Parallel::current()->comm;
  PRECICE_ASSERT(_participantComm != MPI_COMM_NULL, "CoRT communication needs a valid participant MPI communicator.");

  int ret = MPI_Comm_split_type(_participantComm, MPI_COMM_TYPE_SHARED, utils::IntraComm::getRank(), MPI_INFO_NULL, &_localComm);
  PRECICE_ASSERT(ret == MPI_SUCCESS, "MPI_Comm_split_type failed for CoRT communication.");

  MPI_Comm_rank(_localComm, &_localRank);
  MPI_Comm_size(_localComm, &_localSize);

  _isProxy   = (_localRank == 0);
  _proxyRank = utils::IntraComm::getRank();
  MPI_Bcast(&_proxyRank, 1, MPI_INT, 0, _localComm);

  _proxyRankByRank.resize(utils::IntraComm::getSize());
  MPI_Allgather(&_proxyRank, 1, MPI_INT, _proxyRankByRank.data(), 1, MPI_INT, _participantComm);

  _initialized = true;
#endif
}

void NodeTopology::reset()
{
#ifndef PRECICE_NO_MPI
  if (_localComm != MPI_COMM_NULL) {
    MPI_Comm_free(&_localComm);
    _localComm = MPI_COMM_NULL;
  }
  _participantComm = MPI_COMM_NULL;
#endif

  _initialized = false;
  _isProxy     = false;
  _localRank   = 0;
  _localSize   = 1;
  _proxyRank   = 0;
  _proxyRankByRank.clear();
}

bool NodeTopology::isInitialized() const
{
  return _initialized;
}

bool NodeTopology::isProxy() const
{
  return _isProxy;
}

int NodeTopology::localRank() const
{
  return _localRank;
}

int NodeTopology::localSize() const
{
  return _localSize;
}

int NodeTopology::proxyRank() const
{
  return _proxyRank;
}

const std::vector<int> &NodeTopology::proxyRankByRank() const
{
  return _proxyRankByRank;
}

void NodeTopology::broadcastFromPrimary(std::vector<int> &values) const
{
#ifndef PRECICE_NO_MPI
  int size = static_cast<int>(values.size());
  MPI_Bcast(&size, 1, MPI_INT, 0, _participantComm);
  values.resize(size);
  if (size > 0) {
    MPI_Bcast(values.data(), size, MPI_INT, 0, _participantComm);
  }
#else
  (void) values;
#endif
}

void NodeTopology::broadcastFromLocalProxy(int &value) const
{
#ifndef PRECICE_NO_MPI
  MPI_Bcast(&value, 1, MPI_INT, 0, _localComm);
#else
  (void) value;
#endif
}

void NodeTopology::broadcastFromLocalProxy(std::vector<int> &values) const
{
#ifndef PRECICE_NO_MPI
  int size = static_cast<int>(values.size());
  MPI_Bcast(&size, 1, MPI_INT, 0, _localComm);
  values.resize(size);
  if (size > 0) {
    MPI_Bcast(values.data(), size, MPI_INT, 0, _localComm);
  }
#else
  (void) values;
#endif
}

void NodeTopology::localBarrier() const
{
#ifndef PRECICE_NO_MPI
  MPI_Barrier(_localComm);
#endif
}

#ifndef PRECICE_NO_MPI
utils::Parallel::Communicator NodeTopology::participantCommunicator() const
{
  return _participantComm;
}

utils::Parallel::Communicator NodeTopology::localCommunicator() const
{
  return _localComm;
}
#endif

} // namespace precice::m2n::cort
