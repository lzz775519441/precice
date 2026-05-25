#include "CoRTComFactory.hpp"

#include "CoRTCommunication.hpp"
#include "m2n/CoRTPluginSymbols.hpp"

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

#if defined(_WIN32)
#define PRECICE_CORT_PLUGIN_EXPORT __declspec(dllexport)
#else
#define PRECICE_CORT_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

extern "C" {

PRECICE_CORT_PLUGIN_EXPORT int precice_cort_plugin_abi_version()
{
  return precice::m2n::cortplugin::ABI_VERSION;
}

PRECICE_CORT_PLUGIN_EXPORT precice::m2n::DistributedComFactory *
precice_cort_create_distributed_com_factory(
    const precice::com::PtrCommunicationFactory *communicationFactory)
{
  if (communicationFactory == nullptr) {
    return nullptr;
  }
  return new precice::m2n::CoRTComFactory(*communicationFactory);
}

PRECICE_CORT_PLUGIN_EXPORT void
precice_cort_destroy_distributed_com_factory(precice::m2n::DistributedComFactory *factory)
{
  delete factory;
}

} // extern "C"

#undef PRECICE_CORT_PLUGIN_EXPORT
