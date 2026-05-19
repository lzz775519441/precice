#include "m2n/cort/DoubleDataExchange.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>
#include <utility>

#include "logging/LogMacros.hpp"
#include "m2n/cort/NodeTopology.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

#ifndef PRECICE_NO_MPI
#include <mpi.h>
#endif

namespace precice::m2n::cort {

struct DoubleDataExchange::Impl {
  explicit Impl(const NodeTopology &topology)
      : topology(topology)
  {
  }

  const NodeTopology &topology;

  std::vector<RankMapping> mappings;
  std::map<int, int>       remoteRankToProxy;
  PtrTransport             transport;

  std::vector<PtrRequest> ongoingRequests;

  struct ProxyTask {
    int  targetRank   = -1;
    long shmOffset    = 0;
    long totalDoubles = 0;
  };

  struct WorkerTask {
    double    *shmPtr     = nullptr;
    const int *indicesPtr = nullptr;
    size_t     count      = 0;
  };

  int cachedSendDim = 0;
  int cachedRecvDim = 0;

#ifndef PRECICE_NO_MPI
  MPI_Win winSend = MPI_WIN_NULL;
  MPI_Win winRecv = MPI_WIN_NULL;
#endif

  char *sendBasePtr = nullptr;
  char *recvBasePtr = nullptr;

  std::vector<ProxyTask>  proxySendTasks;
  std::vector<WorkerTask> workerSendTasks;
  std::vector<ProxyTask>  proxyRecvTasks;
  std::vector<WorkerTask> workerRecvTasks;

  void configure(std::vector<RankMapping> newMappings,
                 std::map<int, int>       newRemoteRankToProxy,
                 PtrTransport             newTransport)
  {
    wait();
    freeSendWindow();
    freeRecvWindow();

    mappings          = std::move(newMappings);
    remoteRankToProxy = std::move(newRemoteRankToProxy);
    transport         = std::move(newTransport);

    cachedSendDim = 0;
    cachedRecvDim = 0;
  }

  void reset()
  {
    wait();
    freeSendWindow();
    freeRecvWindow();
    mappings.clear();
    remoteRankToProxy.clear();
    transport.reset();
  }

  void wait()
  {
    for (auto &request : ongoingRequests) {
      request->wait();
    }
    ongoingRequests.clear();
  }

  void initializePatterns(int valueDimension)
  {
    PRECICE_ASSERT(valueDimension > 0, "CoRT communication requires a positive value dimension.");

    if (topology.isProxy()) {
      wait();
    }
    topology.localBarrier();

    if (valueDimension != cachedSendDim) {
      initializeSendPattern(valueDimension);
    }
    if (valueDimension != cachedRecvDim) {
      initializeRecvPattern(valueDimension);
    }

    topology.localBarrier();
  }

  void send(precice::span<double const> itemsToSend, int valueDimension)
  {
    PRECICE_ASSERT(valueDimension > 0, "CoRT communication requires a positive value dimension.");

    if (topology.isProxy()) {
      wait();
    }
    topology.localBarrier();

    if (valueDimension != cachedSendDim) {
      initializeSendPattern(valueDimension);
    }

    for (const auto &task : workerSendTasks) {
      for (size_t i = 0; i < task.count; ++i) {
        const int vertexIndex = task.indicesPtr[i];
        for (int d = 0; d < valueDimension; ++d) {
          task.shmPtr[i * valueDimension + d] = itemsToSend[vertexIndex * valueDimension + d];
        }
      }
    }

    topology.localBarrier();

    if (topology.isProxy() && transport) {
      for (const auto &task : proxySendTasks) {
        double *sendBuf = reinterpret_cast<double *>(sendBasePtr + task.shmOffset);
        ongoingRequests.push_back(transport->asyncSend(precice::span<double const>(sendBuf, task.totalDoubles), task.targetRank));
      }
    }
  }

  void receive(precice::span<double> itemsToReceive, int valueDimension)
  {
    PRECICE_ASSERT(valueDimension > 0, "CoRT communication requires a positive value dimension.");

    std::fill(itemsToReceive.begin(), itemsToReceive.end(), 0.0);

    if (topology.isProxy()) {
      wait();
    }
    topology.localBarrier();

    if (valueDimension != cachedRecvDim) {
      initializeRecvPattern(valueDimension);
    }

    if (topology.isProxy() && transport) {
      std::vector<PtrRequest> currentRecvRequests;
      for (const auto &task : proxyRecvTasks) {
        double *recvBuf = reinterpret_cast<double *>(recvBasePtr + task.shmOffset);
        currentRecvRequests.push_back(transport->asyncReceive(precice::span<double>(recvBuf, task.totalDoubles), task.targetRank));
      }
      for (auto &request : currentRecvRequests) {
        request->wait();
      }
    }

    topology.localBarrier();

    for (const auto &task : workerRecvTasks) {
      for (size_t i = 0; i < task.count; ++i) {
        const int vertexIndex = task.indicesPtr[i];
        for (int d = 0; d < valueDimension; ++d) {
          itemsToReceive[vertexIndex * valueDimension + d] += task.shmPtr[i * valueDimension + d];
        }
      }
    }
  }

  void initializeSendPattern(int valueDimension)
  {
#ifdef PRECICE_NO_MPI
    (void) valueDimension;
    throw std::runtime_error{"CoRT communication requires preCICE to be compiled with MPI support."};
#else
    freeSendWindow();

    std::vector<long> localMeta;
    const int         myRank = utils::IntraComm::getRank();
    for (const auto &mapping : mappings) {
      auto proxy = remoteRankToProxy.find(mapping.remoteRank);
      if (proxy == remoteRankToProxy.end()) {
        continue;
      }

      const long bytes = static_cast<long>(mapping.indices.size() * valueDimension * sizeof(double));
      if (bytes > 0) {
        localMeta.push_back(proxy->second);
        localMeta.push_back(mapping.remoteRank);
        localMeta.push_back(myRank);
        localMeta.push_back(bytes);
      }
    }

    std::vector<long> globalMeta = allGatherLocalLongs(localMeta);

    struct Entry {
      long target;
      long remote;
      long source;
      long bytes;
      long offset;
    };

    std::vector<Entry> requests;
    requests.reserve(globalMeta.size() / 4);
    for (size_t i = 0; i < globalMeta.size(); i += 4) {
      requests.push_back({globalMeta[i], globalMeta[i + 1], globalMeta[i + 2], globalMeta[i + 3], 0});
    }

    std::sort(requests.begin(), requests.end(), [](const Entry &a, const Entry &b) {
      if (a.target != b.target) {
        return a.target < b.target;
      }
      if (a.remote != b.remote) {
        return a.remote < b.remote;
      }
      return a.source < b.source;
    });

    long currentOffset = 0;
    for (auto &request : requests) {
      request.offset = currentOffset;
      if (proxySendTasks.empty() || proxySendTasks.back().targetRank != request.target) {
        proxySendTasks.push_back({static_cast<int>(request.target), currentOffset, 0});
      }
      proxySendTasks.back().totalDoubles += request.bytes / static_cast<long>(sizeof(double));
      currentOffset += request.bytes;
    }

    allocateSharedWindow(currentOffset, &sendBasePtr, &winSend);

    for (const auto &request : requests) {
      if (request.source != myRank) {
        continue;
      }

      auto mapping = std::find_if(mappings.begin(), mappings.end(), [&](const RankMapping &candidate) {
        return static_cast<long>(candidate.remoteRank) == request.remote;
      });
      if (mapping != mappings.end()) {
        workerSendTasks.push_back({reinterpret_cast<double *>(sendBasePtr + request.offset),
                                   mapping->indices.data(),
                                   mapping->indices.size()});
      }
    }

    cachedSendDim = valueDimension;
#endif
  }

  void initializeRecvPattern(int valueDimension)
  {
#ifdef PRECICE_NO_MPI
    (void) valueDimension;
    throw std::runtime_error{"CoRT communication requires preCICE to be compiled with MPI support."};
#else
    freeRecvWindow();

    std::vector<long> localMeta;
    const int         myRank = utils::IntraComm::getRank();
    for (const auto &mapping : mappings) {
      auto proxy = remoteRankToProxy.find(mapping.remoteRank);
      if (proxy == remoteRankToProxy.end()) {
        continue;
      }

      const long bytes = static_cast<long>(mapping.indices.size() * valueDimension * sizeof(double));
      if (bytes > 0) {
        localMeta.push_back(proxy->second);
        localMeta.push_back(myRank);
        localMeta.push_back(mapping.remoteRank);
        localMeta.push_back(bytes);
      }
    }

    std::vector<long> globalMeta = allGatherLocalLongs(localMeta);

    struct Entry {
      long remoteProxy;
      long myRank;
      long remoteRank;
      long bytes;
      long offset;
    };

    std::vector<Entry> requests;
    requests.reserve(globalMeta.size() / 4);
    for (size_t i = 0; i < globalMeta.size(); i += 4) {
      requests.push_back({globalMeta[i], globalMeta[i + 1], globalMeta[i + 2], globalMeta[i + 3], 0});
    }

    std::sort(requests.begin(), requests.end(), [](const Entry &a, const Entry &b) {
      if (a.remoteProxy != b.remoteProxy) {
        return a.remoteProxy < b.remoteProxy;
      }
      if (a.myRank != b.myRank) {
        return a.myRank < b.myRank;
      }
      return a.remoteRank < b.remoteRank;
    });

    long currentOffset = 0;
    for (auto &request : requests) {
      request.offset = currentOffset;
      if (proxyRecvTasks.empty() || proxyRecvTasks.back().targetRank != request.remoteProxy) {
        proxyRecvTasks.push_back({static_cast<int>(request.remoteProxy), currentOffset, 0});
      }
      proxyRecvTasks.back().totalDoubles += request.bytes / static_cast<long>(sizeof(double));
      currentOffset += request.bytes;
    }

    allocateSharedWindow(currentOffset, &recvBasePtr, &winRecv);

    for (const auto &request : requests) {
      if (request.myRank != myRank) {
        continue;
      }

      auto mapping = std::find_if(mappings.begin(), mappings.end(), [&](const RankMapping &candidate) {
        return static_cast<long>(candidate.remoteRank) == request.remoteRank;
      });
      if (mapping != mappings.end()) {
        workerRecvTasks.push_back({reinterpret_cast<double *>(recvBasePtr + request.offset),
                                   mapping->indices.data(),
                                   mapping->indices.size()});
      }
    }

    cachedRecvDim = valueDimension;
#endif
  }

#ifndef PRECICE_NO_MPI
  std::vector<long> allGatherLocalLongs(const std::vector<long> &localValues) const
  {
    std::vector<int> counts(topology.localSize());
    const int        localCount = static_cast<int>(localValues.size());
    MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, topology.localCommunicator());

    std::vector<int> displs(topology.localSize() + 1, 0);
    for (int i = 0; i < topology.localSize(); ++i) {
      displs[i + 1] = displs[i] + counts[i];
    }

    std::vector<long> globalValues(displs.back());
    MPI_Allgatherv(localValues.data(), localCount, MPI_LONG,
                   globalValues.data(), counts.data(), displs.data(), MPI_LONG,
                   topology.localCommunicator());
    return globalValues;
  }

  void allocateSharedWindow(long requiredBytes, char **basePtr, MPI_Win *window) const
  {
    const long totalSize = (requiredBytes + 7) & ~7;

    MPI_Info winInfo;
    MPI_Info_create(&winInfo);
    MPI_Info_set(winInfo, "alloc_shared_noncontig", "true");

    const MPI_Aint localAllocationSize = topology.isProxy() ? static_cast<MPI_Aint>(totalSize) : 0;
    int            ret                 = MPI_Win_allocate_shared(localAllocationSize, sizeof(char), winInfo, topology.localCommunicator(), basePtr, window);
    PRECICE_ASSERT(ret == MPI_SUCCESS, "MPI_Win_allocate_shared failed for CoRT communication.");
    MPI_Info_free(&winInfo);

    if (!topology.isProxy()) {
      MPI_Aint sharedSize = 0;
      int      dispUnit   = 0;
      MPI_Win_shared_query(*window, 0, &sharedSize, &dispUnit, basePtr);
    }
  }
#endif

  void freeSendWindow()
  {
#ifndef PRECICE_NO_MPI
    if (winSend != MPI_WIN_NULL) {
      MPI_Win_free(&winSend);
      winSend = MPI_WIN_NULL;
    }
#endif
    sendBasePtr = nullptr;
    proxySendTasks.clear();
    workerSendTasks.clear();
  }

  void freeRecvWindow()
  {
#ifndef PRECICE_NO_MPI
    if (winRecv != MPI_WIN_NULL) {
      MPI_Win_free(&winRecv);
      winRecv = MPI_WIN_NULL;
    }
#endif
    recvBasePtr = nullptr;
    proxyRecvTasks.clear();
    workerRecvTasks.clear();
  }
};

DoubleDataExchange::DoubleDataExchange(const NodeTopology &topology)
    : _impl(std::make_unique<Impl>(topology))
{
}

DoubleDataExchange::~DoubleDataExchange()
{
  reset();
}

void DoubleDataExchange::configure(std::vector<RankMapping> mappings,
                                   std::map<int, int>       remoteRankToProxy,
                                   PtrTransport             transport)
{
  _impl->configure(std::move(mappings), std::move(remoteRankToProxy), std::move(transport));
}

void DoubleDataExchange::initializePatterns(int valueDimension)
{
  _impl->initializePatterns(valueDimension);
}

void DoubleDataExchange::reset()
{
  _impl->reset();
}

void DoubleDataExchange::send(precice::span<double const> itemsToSend, int valueDimension)
{
  _impl->send(itemsToSend, valueDimension);
}

void DoubleDataExchange::receive(precice::span<double> itemsToReceive, int valueDimension)
{
  _impl->receive(itemsToReceive, valueDimension);
}

void DoubleDataExchange::wait()
{
  _impl->wait();
}

} // namespace precice::m2n::cort
