cmake_minimum_required( VERSION 3.1 )

if( TARGET libInterp )
  return()
endif( TARGET libInterp )


enable_testing()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GitVersion)
GIT_VERSION(libInterp)
project(libInterp VERSION ${libInterp_VERSION})

# need to run configure_file twice
configure_file( "${PROJECT_SOURCE_DIR}/src/version.h.in"
                "${PROJECT_SOURCE_DIR}/src/version.h")
configure_file( "${PROJECT_SOURCE_DIR}/src/version.h"
                "${PROJECT_SOURCE_DIR}/src/version.h")


set( cxx_needed_features
     cxx_auto_type
   )

find_package( Boost REQUIRED )
find_package( Eigen3 REQUIRED )
add_library( Eigen3::eigen3 INTERFACE IMPORTED )
set_property( TARGET Eigen3::eigen3 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIRS})

string( REGEX REPLACE "^lib" "" LIB_NAME ${PROJECT_NAME} )
add_library( ${LIB_NAME} INTERFACE )
target_include_directories( ${LIB_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/libInterp>
)
target_sources( ${LIB_NAME} INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interp.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/LinearInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/InterpolatorBase.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/MonotonicInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/CubicSplineInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/BicubicInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/ThinPlateSplineInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/BilinearInterpolator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/InterpolatorBase.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/Indexing.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/ReadFunction.hpp
)
  target_link_libraries(${LIB_NAME} INTERFACE Boost::boost)
  target_link_libraries(${LIB_NAME} INTERFACE Eigen3::eigen3)

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
add_subdirectory( testing )
endif()
