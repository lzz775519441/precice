#include "CoRTPluginAdapter.hpp"

#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "CoRTPluginSymbols.hpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace precice::m2n {
namespace {

constexpr const char *EnvironmentVariable = "PRECICE_CORT_PLUGIN";

#ifndef PRECICE_CORT_PLUGIN_DEFAULT_PATH
#define PRECICE_CORT_PLUGIN_DEFAULT_PATH "libprecice-cort-plugin.so"
#endif

using ABIversionFunction = int (*)();
using CreateFactoryFunction = DistributedComFactory *(*)(const com::PtrCommunicationFactory *);
using DestroyFactoryFunction = void (*)(DistributedComFactory *);

struct DynamicLibrary {
#ifdef _WIN32
  using Handle = HMODULE;
#else
  using Handle = void *;
#endif

  Handle                 handle = nullptr;
  ABIversionFunction     abiVersion = nullptr;
  CreateFactoryFunction  createFactory = nullptr;
  DestroyFactoryFunction destroyFactory = nullptr;
  std::string            path;

  explicit DynamicLibrary(std::string libraryPath)
      : path(std::move(libraryPath))
  {
#ifdef _WIN32
    handle = LoadLibraryA(path.c_str());
    if (handle == nullptr) {
      throw std::runtime_error{"Could not load CoRT plugin \"" + path + "\"."};
    }
#else
    handle = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (handle == nullptr) {
      const char *error = dlerror();
      throw std::runtime_error{"Could not load CoRT plugin \"" + path + "\": " + (error == nullptr ? "unknown error" : error)};
    }
#endif

    try {
      abiVersion     = loadSymbol<ABIversionFunction>(cortplugin::ABI_VERSION_FUNCTION);
      createFactory  = loadSymbol<CreateFactoryFunction>(cortplugin::CREATE_FACTORY_FUNCTION);
      destroyFactory = loadSymbol<DestroyFactoryFunction>(cortplugin::DESTROY_FACTORY_FUNCTION);

      const int pluginABI = abiVersion();
      if (pluginABI != cortplugin::ABI_VERSION) {
        throw std::runtime_error{"The CoRT plugin ABI version is " + std::to_string(pluginABI) +
                                 ", but this preCICE build requires ABI version " + std::to_string(cortplugin::ABI_VERSION) + "."};
      }
    } catch (...) {
      close();
      throw;
    }
  }

  ~DynamicLibrary()
  {
    close();
  }

  DynamicLibrary(const DynamicLibrary &) = delete;
  DynamicLibrary &operator=(const DynamicLibrary &) = delete;

private:
  void close()
  {
    if (handle == nullptr) {
      return;
    }
#ifdef _WIN32
    FreeLibrary(handle);
#else
    dlclose(handle);
#endif
    handle = nullptr;
  }

  template <typename Function>
  Function loadSymbol(const char *name)
  {
#ifdef _WIN32
    auto symbol = reinterpret_cast<Function>(GetProcAddress(handle, name));
#else
    dlerror();
    auto symbol = reinterpret_cast<Function>(dlsym(handle, name));
    const char *error = dlerror();
    if (error != nullptr) {
      throw std::runtime_error{"Could not load symbol \"" + std::string{name} + "\" from CoRT plugin \"" + path + "\": " + error};
    }
#endif
    if (symbol == nullptr) {
      throw std::runtime_error{"Could not load symbol \"" + std::string{name} + "\" from CoRT plugin \"" + path + "\"."};
    }
    return symbol;
  }
};

std::vector<std::string> candidatePluginPaths()
{
  if (const char *overridePath = std::getenv(EnvironmentVariable)) {
    if (overridePath[0] != '\0') {
      return {overridePath};
    }
  }

  std::vector<std::string> paths;
#ifdef PRECICE_CORT_PLUGIN_BUILD_PATH
  paths.emplace_back(PRECICE_CORT_PLUGIN_BUILD_PATH);
#endif
  paths.emplace_back(PRECICE_CORT_PLUGIN_DEFAULT_PATH);
  paths.emplace_back("./libprecice-cort-plugin.so");
  paths.emplace_back("libprecice-cort-plugin.so");
  return paths;
}

std::shared_ptr<DynamicLibrary> loadCoRTPlugin()
{
  const auto paths = candidatePluginPaths();
  std::ostringstream errors;

  for (const auto &path : paths) {
    try {
      return std::make_shared<DynamicLibrary>(path);
    } catch (const std::exception &e) {
      errors << "\n  - " << e.what();
    }
  }

  throw std::runtime_error{"Could not load the CoRT m2n communication plugin. "
                           "Set PRECICE_CORT_PLUGIN to the plugin shared library path or install "
                           "libprecice-cort-plugin.so to the preCICE plugin directory. Tried:" +
                           errors.str()};
}

class CoRTPluginComFactory : public DistributedComFactory {
public:
  explicit CoRTPluginComFactory(com::PtrCommunicationFactory communicationFactory)
      : _library(loadCoRTPlugin()),
        _factory(_library->createFactory(&communicationFactory))
  {
    if (_factory == nullptr) {
      throw std::runtime_error{"The CoRT plugin returned a null distributed communication factory."};
    }
  }

  ~CoRTPluginComFactory() override
  {
    if (_factory != nullptr) {
      _library->destroyFactory(_factory);
      _factory = nullptr;
    }
  }

  DistributedCommunication::SharedPointer newDistributedCommunication(mesh::PtrMesh mesh) override
  {
    return _factory->newDistributedCommunication(std::move(mesh));
  }

private:
  std::shared_ptr<DynamicLibrary> _library;
  DistributedComFactory          *_factory = nullptr;
};

} // namespace

DistributedComFactory::SharedPointer
createCoRTPluginComFactory(com::PtrCommunicationFactory communicationFactory)
{
  return std::make_shared<CoRTPluginComFactory>(std::move(communicationFactory));
}

} // namespace precice::m2n
