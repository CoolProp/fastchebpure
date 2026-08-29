cmake_minimum_required(VERSION 3.14)
project(fastcheb)
set(CMAKE_CXX_STANDARD 20)

# CMake 4 removed compatibility with cmake_minimum_required(VERSION < 3.5), and
# some transitively fetched dependencies (nlohmann_json via teqp) still declare
# one.  Without this, configuring with CMake >= 4 aborts before any of our own
# code is reached.  Harmless on CMake 3.x, which ignores the variable.
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

include(FetchContent)
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.4.2
)
FetchContent_MakeAvailable(CLI11)

set(TEQP_DISABLED_FACTORIES "CPA,SAFTVRMIE,GENERICSAFT,SQUAREWELL,EXP6,2CLJF,MIE,MULTIFLUIDACTIVITY,PCSAFT,GERG200X,MULTIFLUIDASSOCIATION,LKP,ECSHUBERELY1994,AMMONIAWATERTILLNERROTH,  ADVANCEDMIXCUBIC,SOFTSAFT" CACHE STRING "Disable some unused models")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp")

set(CHEBTOOLS_STATIC_LIBRARY true CACHE INTERNAL "Build ChebTools as a static library")
add_subdirectory ("${CMAKE_CURRENT_SOURCE_DIR}/externals/ChebTools")

# externals/teqp and externals/ChebTools each force CMAKE_OSX_ARCHITECTURES via
# `set(... CACHE INTERNAL ...)` -- INTERNAL implies FORCE -- to "arm64" and
# "x86_64;arm64" respectively.  Whichever add_subdirectory() runs last silently
# wins, so the executables below can be asked to build slices that teqp never
# produced, and the link fails with "symbol(s) not found for architecture
# x86_64".  Pin it explicitly here, after both subdirectories, so the target
# architecture is ours to choose rather than an artifact of ordering.
if(APPLE AND NOT FASTCHEB_OSX_ARCHITECTURES_PINNED)
  set(CMAKE_OSX_ARCHITECTURES "${CMAKE_HOST_SYSTEM_PROCESSOR}" CACHE INTERNAL "archs for osx")
  set(FASTCHEB_OSX_ARCHITECTURES_PINNED TRUE CACHE INTERNAL "guard")
endif()

add_executable(fitcheb fastcheb.cpp main.cpp)
target_link_libraries(fitcheb teqpinterface ChebTools teqpcpp CLI11::CLI11)
target_include_directories(fitcheb
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/autodiff"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals"
)
target_compile_definitions(fitcheb PUBLIC -DTEQP_MULTIPRECISION_ENABLED)
target_compile_definitions(fitcheb PUBLIC -DTEQP_MULTICOMPLEX_ENABLED)

add_executable(bench bench.cpp)
target_link_libraries(bench teqpinterface ChebTools teqpcpp Catch2WithMain)
target_include_directories(bench
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/autodiff"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals"
)
target_compile_definitions(bench PUBLIC -DTEQP_MULTIPRECISION_ENABLED)
target_compile_definitions(bench PUBLIC -DTEQP_MULTICOMPLEX_ENABLED)

add_executable(test_eos_hash test_eos_hash.cpp)
target_link_libraries(test_eos_hash Catch2WithMain)
target_include_directories(test_eos_hash
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals/teqp/externals/nlohmann_json"
)
