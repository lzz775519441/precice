#pragma once

#include <vector>

#include "utils/Parallel.hpp"

namespace precice::m2n::cort {

class NodeTopology {
public:
  NodeTopology() = default;
  ~NodeTopology();

  NodeTopology(const NodeTopology &)            = delete;
  NodeTopology &operator=(const NodeTopology &) = delete;

  NodeTopology(NodeTopology &&other) noexcept;
  NodeTopology &operator=(NodeTopology &&other) noexcept;

  void initialize();
  void reset();

  bool isInitialized() const;
  bool isProxy() const;
  int  localRank() const;
  int  localSize() const;
  int  proxyRank() const;

  const std::vector<int> &proxyRankByRank() const;

  void broadcastFromPrimary(std::vector<int> &values) const;
  void broadcastFromLocalProxy(int &value) const;
  void broadcastFromLocalProxy(std::vector<int> &values) const;
  void localBarrier() const;

#ifndef PRECICE_NO_MPI
  utils::Parallel::Communicator participantCommunicator() const;
  utils::Parallel::Communicator localCommunicator() const;
#endif

private:
  bool _initialized = false;
  bool _isProxy     = false;
  int  _localRank   = 0;
  int  _localSize   = 1;
  int  _proxyRank   = 0;

  std::vector<int> _proxyRankByRank;

#ifndef PRECICE_NO_MPI
  utils::Parallel::Communicator _participantComm = MPI_COMM_NULL;
  utils::Parallel::Communicator _localComm       = MPI_COMM_NULL;
#endif
};

} // namespace precice::m2n::cort
