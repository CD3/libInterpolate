cmake_minimum_required(VERSION 3.16)

option(BUILD_TESTS "Build unit tests" ON)


include("version.cmake")
if(NOT DEFINED libInterpolate_FULL_VERSION)
  message(STATUS "Version not set, attempting to detected with git")
  execute_process(
    COMMAND git describe --tags --always --dirty
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE GIT_DESCRIBE_ERROR
    ERROR_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE GIT_DESCRIBE_RESULT
  )
  if("${GIT_DESCRIBE_RESULT}" STREQUAL "0")
    message(STATUS "Version detected with git")
    set(libInterpolate_FULL_VERSION "${GIT_DESCRIBE_VERSION}" CACHE STRING "Full build version number")
  else()
    message(WARNING "Could NOT detect version number with git")
    message(WARNING "STDERR: ${GIT_DESCRIBE_RESULT} ${GIT_DESCRIBE_ERROR}")
    set(libInterpolate_FULL_VERSION "0.0.0-UNKNOWN" CACHE STRING "WARNING!! Full build version number could not be determined.")
  endif()
endif()

string(REGEX MATCH "^[0-9]+(\\.[0-9]+)*" libInterpolate_PROJECT_VERSION_NUMERIC "${libInterpolate_FULL_VERSION}")

message(STATUS "FULL VERSION: ${libInterpolate_FULL_VERSION}")  
message(STATUS "NUMERIC VERSION NUMBER: ${libInterpolate_PROJECT_VERSION_NUMERIC}")  
# NOTE: this will create a variables named libInterpolate_VERSION
project(libInterpolate VERSION ${libInterpolate_PROJECT_VERSION_NUMERIC})

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(CMakePrintHelpers)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Boost REQUIRED)
find_package(Eigen3 3.3.7 REQUIRED) # v3.3.7 adds support for cmake targets

add_library(Interpolate INTERFACE)
add_library(libInterpolate::Interpolate ALIAS Interpolate)
target_include_directories(
  Interpolate
  INTERFACE
    $<BUILD_INTERFACE:${libInterpolate_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_compile_definitions(Interpolate
                           INTERFACE $<$<BOOL:${MSVC}>:_USE_MATH_DEFINES>)
target_compile_features(Interpolate INTERFACE cxx_std_17)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/libInterpolate/version.h"
)
target_sources(
  Interpolate
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolate.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/AnyInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/InterpolatorBase.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/ThinPlateSplineInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/BilinearInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/AnyInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/BicubicInterpolator.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/libInterpolate/Interpolators/_2D/NearestNeighborInterpolator.hpp>
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
  enable_testing()
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
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libInterpolate)
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
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libInterpolate)

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
