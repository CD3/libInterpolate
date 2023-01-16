cmake_minimum_required(VERSION 3.16)
project(libInterpolate)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CMakePrintHelpers)

option(BUILD_TESTS "Build unit tests" ON)

if(${BUILD_TESTS})
  enable_testing()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(macro-git_version)
git_version(libInterpolate)
message(STATUS "libInterpolate version: ${libInterpolate_VERSION}")

find_package(Boost REQUIRED)
find_package(Eigen3 3.3.7 REQUIRED) # v3.3.7 adds support for cmake targets

# need to run configure_file twice
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/libInterpolate_version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/libInterpolate/libInterpolate_version.h")
configure_file(
  "${CMAKE_CURRENT_BINARY_DIR}/include/libInterpolate/libInterpolate_version.h"
  "${CMAKE_CURRENT_BINARY_DIR}/include/libInterpolate/libInterpolate_version.h")

string(REGEX REPLACE "^lib" "" LIB_NAME libInterpolate)
add_library(Interpolate INTERFACE)
add_library(libInterpolate::Interpolate ALIAS Interpolate)
target_include_directories(
  Interpolate
  INTERFACE
    $<BUILD_INTERFACE:${libInterpolate_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${libInterpolate_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(Interpolate
                           INTERFACE $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>)
set_target_properties(Interpolate PROPERTIES CXX_STANDARD 17)
target_sources(
  Interpolate
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interp.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolate.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/AnyInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/InterpolatorBase.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/ThinPlateSplineInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/BilinearInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/AnyInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/BicubicInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_1D/InterpolatorBase.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_1D/CubicSplineInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_1D/AnyInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_1D/MonotonicInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_1D/LinearInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Utils/Indexing.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Utils/ReadFunction.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Utils/Concepts.hpp>
)

target_link_libraries(Interpolate INTERFACE Boost::boost Eigen3::Eigen)

if(${BUILD_TESTS})
  add_subdirectory(testing)
endif()

install(
  TARGETS Interpolate
  EXPORT libInterpolateTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include)
install(
  DIRECTORY src/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  EXPORT libInterpolateTargets
  FILE libInterpolateTargets.cmake
  NAMESPACE libInterpolate::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/libInterpolateConfig.cmake
  "include(CMakeFindDependencyMacro)
find_dependency(Boost)
find_dependency(Eigen3 3.3.7 REQUIRED)
include(\${CMAKE_CURRENT_LIST_DIR}/libInterpolateTargets.cmake)
")
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libInterpolateConfigVersion.cmake
  VERSION ${libInterpolate_VERSION}
  COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/libInterpolateConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libInterpolateConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "A C++ library for numerical interpolation.")
set(CPACK_PACKAGE_VENDOR "C.D. Clark III")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
if(libInterpolate_VERSION_MAJOR)
  set(CPACK_PACKAGE_VERSION_MAJOR ${libInterpolate_VERSION_MAJOR})
endif()
if(libInterpolate_VERSION_MINOR)
  set(CPACK_PACKAGE_VERSION_MINOR ${libInterpolate_VERSION_MINOR})
endif()
if(libInterpolate_VERSION_PATCH)
  set(CPACK_PACKAGE_VERSION_PATCH ${libInterpolate_VERSION_PATCH})
endif()
include(CPack)

# vvvvDEPRECATEDvvvv

# the code below is provided for backwards compatibility for when the library
# was named libInterp instead of libInterpolate

add_library(Interp INTERFACE)
add_library(libInterp::Interp ALIAS Interp)
target_include_directories(
  Interp
  INTERFACE
    $<BUILD_INTERFACE:${libInterpolate_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${libInterpolate_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${libInterpolate_SOURCE_DIR}/src/libInterpolate>
    $<BUILD_INTERFACE:${libInterpolate_BINARY_DIR}/include/libInterpolate>
    $<INSTALL_INTERFACE:include/libInterpolate>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(Interp
                           INTERFACE $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>)
target_link_libraries(Interp INTERFACE Boost::boost Eigen3::Eigen)

install(TARGETS Interp EXPORT libInterpTargets)
install(
  EXPORT libInterpTargets
  FILE libInterpTargets.cmake
  NAMESPACE libInterp::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/libInterpConfig.cmake
  "message(WARNING \"libInterp has been renamed libInterpolate. Please use find_package(libInterpolate) and the libInterpolate::Interpolate target.\")
include(CMakeFindDependencyMacro)
find_dependency(Boost)
find_dependency(Eigen3 3.3.7 REQUIRED)
include(\${CMAKE_CURRENT_LIST_DIR}/libInterpTargets.cmake)
")
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/libInterpConfigVersion.cmake
  VERSION ${libInterpolate_VERSION}
  COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/libInterpConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libInterpConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake)
