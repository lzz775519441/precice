#pragma once

#include <memory>
#include <vector>

#include "precice/span.hpp"

namespace precice::m2n::cort {

struct RankMapping {
  int              remoteRank = -1;
  std::vector<int> indices;
};

class Request {
public:
  virtual ~Request() = default;

  virtual void wait() = 0;
};

using PtrRequest = std::shared_ptr<Request>;

class Transport {
public:
  virtual ~Transport() = default;

  virtual PtrRequest asyncSend(precice::span<double const> itemsToSend, int remoteRank) = 0;

  virtual PtrRequest asyncReceive(precice::span<double> itemsToReceive, int remoteRank) = 0;
};

using PtrTransport = std::shared_ptr<Transport>;

} // namespace precice::m2n::cort
