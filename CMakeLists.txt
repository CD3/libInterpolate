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
target_sources(
  ${LIB_NAME}
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interp.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/version.h>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/ReadFunction.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Utils/Indexing.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/InterpolatorBase.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/ThinPlateSplineInterpolator.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/BicubicInterpolator.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_2D/BilinearInterpolator.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/LinearInterpolator.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/CubicSplineInterpolator.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/InterpolatorBase.hpp>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolators/_1D/MonotonicInterpolator.hpp>
)

target_link_libraries(${LIB_NAME} INTERFACE Boost::boost)
target_link_libraries(${LIB_NAME} INTERFACE Eigen3::eigen3)

install( TARGETS ${LIB_NAME} EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/${PROJECT_NAME}
)
install(DIRECTORY src/
  DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING PATTERN "*.hpp"
)
install(EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION lib/cmake/${PROJECT_NAME}
)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
"include(CMakeFindDependencyMacro)
find_dependency(Boost)
find_dependency( Eigen3 )
add_library( Eigen3::eigen3 INTERFACE IMPORTED )
set_property( TARGET Eigen3::eigen3 PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})
include(\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake)
"
  )
include(CMakePackageConfigHelpers)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION lib/cmake/${PROJECT_NAME}
  )

if( STANDALONE )
add_subdirectory( testing )
endif()
