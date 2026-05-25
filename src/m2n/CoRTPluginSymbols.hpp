#pragma once

namespace precice::m2n::cortplugin {

constexpr int ABI_VERSION = 1;

constexpr const char *ABI_VERSION_FUNCTION = "precice_cort_plugin_abi_version";
constexpr const char *CREATE_FACTORY_FUNCTION = "precice_cort_create_distributed_com_factory";
constexpr const char *DESTROY_FACTORY_FUNCTION = "precice_cort_destroy_distributed_com_factory";

} // namespace precice::m2n::cortplugin
