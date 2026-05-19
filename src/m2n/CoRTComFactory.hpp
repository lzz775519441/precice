#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::m2n {

class CoRTComFactory : public DistributedComFactory {

public:
  explicit CoRTComFactory(com::PtrCommunicationFactory communicationFactory);

  DistributedCommunication::SharedPointer newDistributedCommunication(
      mesh::PtrMesh mesh) override;

private:
  com::PtrCommunicationFactory _comFactory;
};

} // namespace precice::m2n