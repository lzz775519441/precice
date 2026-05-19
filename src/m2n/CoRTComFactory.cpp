#include "CoRTComFactory.hpp"

#include "CoRTCommunication.hpp"

#include <utility>

#include "com/SharedPointer.hpp"

namespace precice::m2n {

CoRTComFactory::CoRTComFactory(com::PtrCommunicationFactory comFactory)
    : _comFactory(std::move(comFactory)) {}

DistributedCommunication::SharedPointer
CoRTComFactory::newDistributedCommunication(mesh::PtrMesh mesh)
{
  return DistributedCommunication::SharedPointer(new CoRTCommunication(_comFactory, mesh));
}

} // namespace precice::m2n