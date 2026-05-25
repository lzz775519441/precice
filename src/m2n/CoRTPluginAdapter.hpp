#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"

namespace precice::m2n {

DistributedComFactory::SharedPointer
createCoRTPluginComFactory(com::PtrCommunicationFactory communicationFactory);

} // namespace precice::m2n
