cmake_minimum_required( VERSION 3.1 )

if( TARGET libInterp )
  return()
endif( TARGET libInterp )

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  set(STANDALONE ON)
endif( )

enable_testing()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(libInterp)
if( STANDALONE )
  include(GitVersion)
  GIT_VERSION(${PROJECT_NAME})
else()
  set(${PROJECT_NAME}_VERSION "X.X-subproject")
endif()

# need to run configure_file twice
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h.in"
                "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h")
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h"
                "${CMAKE_CURRENT_SOURCE_DIR}/src/version.h")


set( cxx_needed_features
     cxx_auto_type
   )

find_package( Boost REQUIRED )
find_package( Eigen3 REQUIRED )
add_library( Eigen3::eigen3 INTERFACE IMPORTED )
set_property( TARGET Eigen3::eigen3 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})

string( REGEX REPLACE "^lib" "" LIB_NAME ${PROJECT_NAME} )
add_library( ${LIB_NAME} INTERFACE )
target_include_directories( ${LIB_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h
)

target_link_libraries(${LIB_NAME} INTERFACE Boost::boost)
target_link_libraries(${LIB_NAME} INTERFACE Eigen3::eigen3)

if( STANDALONE )
add_subdirectory( testing )
endif()
