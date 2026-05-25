#include "com/SharedPointer.hpp"
#include "m2n/CoRTPluginSymbols.hpp"
#include "m2n/DistributedComFactory.hpp"

#if defined(_WIN32)
#define PRECICE_CORT_TEST_PLUGIN_EXPORT __declspec(dllexport)
#else
#define PRECICE_CORT_TEST_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

extern "C" {

PRECICE_CORT_TEST_PLUGIN_EXPORT int precice_cort_plugin_abi_version()
{
  return precice::m2n::cortplugin::ABI_VERSION + 1;
}

PRECICE_CORT_TEST_PLUGIN_EXPORT precice::m2n::DistributedComFactory *
precice_cort_create_distributed_com_factory(
    const precice::com::PtrCommunicationFactory *)
{
  return nullptr;
}

PRECICE_CORT_TEST_PLUGIN_EXPORT void
precice_cort_destroy_distributed_com_factory(precice::m2n::DistributedComFactory *)
{
}

} // extern "C"

#undef PRECICE_CORT_TEST_PLUGIN_EXPORT
