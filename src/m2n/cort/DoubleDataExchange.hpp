#pragma once

#include <map>
#include <memory>
#include <vector>

#include "m2n/cort/Transport.hpp"
#include "precice/span.hpp"

namespace precice::m2n::cort {

class NodeTopology;

class DoubleDataExchange {
public:
  explicit DoubleDataExchange(const NodeTopology &topology);
  ~DoubleDataExchange();

  DoubleDataExchange(const DoubleDataExchange &)            = delete;
  DoubleDataExchange &operator=(const DoubleDataExchange &) = delete;

  void configure(std::vector<RankMapping> mappings,
                 std::map<int, int>       remoteRankToProxy,
                 PtrTransport             transport);

  void initializePatterns(int valueDimension);

  void reset();

  void send(precice::span<double const> itemsToSend, int valueDimension);
  void receive(precice::span<double> itemsToReceive, int valueDimension);

  void wait();

private:
  struct Impl;
  std::unique_ptr<Impl> _impl;
};

} // namespace precice::m2n::cort
