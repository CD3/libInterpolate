#          _                                      _        
# ___  ___| |_ _   _ _ __     ___ _ __ ___   __ _| | _____ 
#/ __|/ _ \ __| | | | '_ \   / __| '_ ` _ \ / _` | |/ / _ \
#\__ \  __/ |_| |_| | |_) | | (__| | | | | | (_| |   <  __/
#|___/\___|\__|\__,_| .__/   \___|_| |_| |_|\__,_|_|\_\___|
#                   |_|                                    
cmake_minimum_required( VERSION 3.1 ) # version 3.1 introduced the 'target_comple_features' command

include(CheckIncludeFiles) 
enable_testing()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
# allow compile commands to be exported for YCM
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



#                 _           _                      __ _                       _   _             
# _ __  _ __ ___ (_) ___  ___| |_    ___ ___  _ __  / _(_) __ _ _   _ _ __ __ _| |_(_) ___  _ __  
#| '_ \| '__/ _ \| |/ _ \/ __| __|  / __/ _ \| '_ \| |_| |/ _` | | | | '__/ _` | __| |/ _ \| '_ \ 
#| |_) | | | (_) | |  __/ (__| |_  | (_| (_) | | | |  _| | (_| | |_| | | | (_| | |_| | (_) | | | |
#| .__/|_|  \___// |\___|\___|\__|  \___\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|\__|_|\___/|_| |_|
#|_|           |__/                                       |___/                                   

project(libInterp) # replace <Project Name> with your project's name

# Sometimes you need/want to know if the project is being build as a standalone project (the user
# is calling your CMakeLists.txt file directory) or as a sub-project (the user is calling
# a CMakeLists.txt file for a project that is calling your CMakeLists.txt file).
# This code detects if we are being called from another project by
# checking if the cmake source directory is different than the cmake current source directory.
set( STANDALONE_CONFIG ON )
if( NOT (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR}) )
    message( STATUS "${PROJECT_NAME} being configured as subdirectory" )
    set( STANDALONE_CONFIG OFF )
else()
    message( STATUS "${PROJECT_NAME} being configured as stand-alone" )
endif()

# we want to turn on C++11 support

# You can specify a list of features that you need the compler to support here (see CMake documentation for a list of features).
# The cxx_needed_features variable is passed to `tartget_compile_features` command later one
# when the targets are defined.
set( cxx_needed_features
   )


# OPTIONS
# We only want to setup options if we are in stand-alone mode because
# the calling script may want to set these, and we shouldn't overwrite them.
if (STANDALONE_CONFIG)
  option( BUILD_SWIG_WRAPPERS    "Build interface wrapper libraries for other languages with SWIG" OFF )
  option( GENERATE_SWIG_WRAPPERS "Generate new SWIG wrappers to build. The
                                  differentce between this and BUILD_SWIG_WRAPPERS is that the latter
                                  option will build wrapper libraries from
                                  existing SWIG generated files. If you want
                                  to generate new files to build with,
                                  you should turn both of these options on" OFF )
endif()

# only include CPack if we are being configured as a standalone package
if (STANDALONE_CONFIG)
  include(CPack)
endif()








#                    _                   _      _            _   _             
#__   _____ _ __ ___(_) ___  _ __     __| | ___| |_ ___  ___| |_(_) ___  _ __  
#\ \ / / _ \ '__/ __| |/ _ \| '_ \   / _` |/ _ \ __/ _ \/ __| __| |/ _ \| '_ \ 
# \ V /  __/ |  \__ \ | (_) | | | | | (_| |  __/ ||  __/ (__| |_| | (_) | | | |
#  \_/ \___|_|  |___/_|\___/|_| |_|  \__,_|\___|\__\___|\___|\__|_|\___/|_| |_|
                                                                              


if( NOT ${PROJECT_NAME}_SKIP_VERSION_DETECT )
include(GitVersion)
GIT_VERSION()
endif( NOT ${PROJECT_NAME}_SKIP_VERSION_DETECT )

set( CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR} )
set( CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR} )
set( CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH} )


# rather than call include_directories multiple times
# we will collect all the include directories in a list named _${PROJECT_NAME}_INCLUDE_DIRS.
# then we can add them all at once with a single include_directories call.
# we will also "export" the list by setting a CACHE variable named ${PROJECT_NAME}_INCLUDE_DIRS.
# this will allow other proejcts to access the list if they need.
# we'll do the same for definitions and anything else that might be useful for external projects
# to have access to.
set( _${PROJECT_NAME}_INCLUDE_DIRS "" )
set( _${PROJECT_NAME}_DEFINITIONS "" )



#     _                           _                              _               _        
#  __| | ___ _ __   ___ _ __   __| | ___ _ __   ___ _   _    ___| |__   ___  ___| | _____ 
# / _` |/ _ \ '_ \ / _ \ '_ \ / _` |/ _ \ '_ \ / __| | | |  / __| '_ \ / _ \/ __| |/ / __|
#| (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (__| |_| | | (__| | | |  __/ (__|   <\__ \
# \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___|\__, |  \___|_| |_|\___|\___|_|\_\___/
#           |_|                                     |___/                                 

# Look for project dependencies.
# The code below looks for common dependencies.
# If you have more dependencies, add checks for them here.

# BOOST
# static link to boost
set(Boost_USE_STATIC_LIBS ON )
set(Boost_USE_DYN_LIBS OFF )
# dynamic link to boost
#list( APPEND _${PROJECT_NAME}_DEFINITIONS -DBOOST_ALL_DYN_LINK )

# If you need a specific boost component, check for it here.
MESSAGE( "NOTE: Checking for various boost components. If a specific
component is not found, even if that component is optional, CMake will print a
message saying that Boost was NOT found. This does not mean that no boost
installation was found." )
# The boost command line argument parser.
find_package( Boost COMPONENTS program_options)
list( APPEND _${PROJECT_NAME}_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )
set( Boost_LIBS ${Boost_LIBS} ${Boost_LIBRARIES} )

# The boost header libraries
#find_package( Boost REQUIRED ) # if boost is a required dependency, use this
find_package( Boost )
list( APPEND _${PROJECT_NAME}_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} )

# The eigen linear algebra library.
find_package( Eigen3 REQUIRED )
list( APPEND _${PROJECT_NAME}_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} )

# _           _ _     _ 
#| |__  _   _(_) | __| |
#| '_ \| | | | | |/ _` |
#| |_) | |_| | | | (_| |
#|_.__/ \__,_|_|_|\__,_|
                       

# Add any specific comiler defines that you need. These should be
# defines that your source code uses. Don't add things like --std=g++11 for example,
# this will be taken care of by target_compile_features command.
list( APPEND _${PROJECT_NAME}_DEFINITIONS "" )

# The configure_file command is used to create a source file from a template.
# For example, we can insert version information into our source code by creating
# a version.h.in file and letting CMake create a version.h file with the version
# information detected above.
configure_file( "${PROJECT_SOURCE_DIR}/src/version.h.in"
                "${PROJECT_SOURCE_DIR}/src/version.h")
configure_file( "${PROJECT_SOURCE_DIR}/src/version.h"
                "${PROJECT_SOURCE_DIR}/src/version.h")


# add the source directory to the includes list
list( APPEND _${PROJECT_NAME}_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")

# Add the source code for your project here.
# This command will find all *.cpp, *.hpp, and *.h files in the source directory.
# While not required to build the project, you should
# add .h and .hpp files here so that they will be displayed
# in IDE's like Visual Studio.
file( GLOB_RECURSE ${PROJECT_NAME}_DETECTED_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "src/*.cpp" )
message(STATUS "Detected Source Files:")
foreach(detectedSrc ${${PROJECT_NAME}_DETECTED_SOURCES})
  message(STATUS "  ${detectedSrc}" )
endforeach()
# add additional source code files here.
# Note: you should prefix filenames with ${PROJECT_SOURCE_DIR} so that relative path
#       names don't cause any problems
# Add any additional sources here. For example, if you want to compile source from
# an included external dependency, you could add the sources here.
set( ${PROJECT_NAME}_SOURCES ${${PROJECT_NAME}_DETECTED_SOURCES}
   # ${PROJECT_SOURCE_DIR}/src/SourceCode.cpp
   )


 # add the includes, definitions, etc
 include_directories( ${_${PROJECT_NAME}_INCLUDE_DIRS} )
 add_definitions( ${_${PROJECT_NAME}_DEFINITIONS} )



# If your project is a library, use add_library. The file that is created will be named lib<TARGET>, so
# if your project name is libMPE for example, you would want to use MPE for <TARGET>.
string( REGEX REPLACE "^lib" "" LIB_NAME ${PROJECT_NAME} )
message( STATUS "library target name: ${LIB_NAME}" )
#add_library( ${LIB_NAME} STATIC ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_EXTERNAL_SOURCES})
#set_target_properties( ${LIB_NAME} PROPERTIES VERSION ${${PROJECT_NAME}_VERSION} )
#target_link_libraries( ${LIB_NAME} ${Boost_LIBS} ) # add any link dependencies
#target_compile_features( ${LIB_NAME} PRIVATE ${cxx_needed_features})

# Build applications

if( Boost_PROGRAM_OPTIONS_FOUND )
  add_executable( interp-cli applications/interp-cli.cpp )
  target_link_libraries( interp-cli ${Boost_PROGRAM_OPTIONS_LIBRARY} )
  set_target_properties( interp-cli PROPERTIES CXX_STANDARD 11 )
endif()




#           _         _ _               _             _           
# ___ _   _| |__   __| (_)_ __ ___  ___| |_ ___  _ __(_) ___  ___ 
#/ __| | | | '_ \ / _` | | '__/ _ \/ __| __/ _ \| '__| |/ _ \/ __|
#\__ \ |_| | |_) | (_| | | | |  __/ (__| || (_) | |  | |  __/\__ \
#|___/\__,_|_.__/ \__,_|_|_|  \___|\___|\__\___/|_|  |_|\___||___/
                                                                 


# These directories should contain their own CMakeLists.txt files.
add_subdirectory( testing )
add_subdirectory( doc )
if( BUILD_SWIG_WRAPPERS )
add_subdirectory( swig )
endif( BUILD_SWIG_WRAPPERS )

#                _          
#  ___ __ _  ___| |__   ___ 
# / __/ _` |/ __| '_ \ / _ \
#| (_| (_| | (__| | | |  __/
# \___\__,_|\___|_| |_|\___|


# add variables to the cache so other projects can see them
set( ${PROJECT_NAME}_INCLUDE_DIRS ${_${PROJECT_NAME}_INCLUDE_DIRS} CACHE INTERNAL "${PROJECT_NAME} Include Dirs:" FORCE)
set( ${PROJECT_NAME}_DEFINITIONS  ${_${PROJECT_NAME}_DEFINITIONS}  CACHE INTERNAL "${PROJECT_NAME} Compiler Definitions:" FORCE)
set( ${PROJECT_NAME}_CONFIGURED   1                                CACHE INTERNAL "${PROJECT_NAME} Configured:" FORCE) # that indicates this project has been configured
