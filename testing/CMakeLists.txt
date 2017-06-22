cmake_minimum_required( VERSION 3.1 )
message( "Configuring unit test build." )

OPTION( BUILD_UNIT_TESTS "Build unit tests for the library" ON )
OPTION( BUILD_NONIUS_BENCHMARKS "Build Nonius-based benchmarks." OFF )
OPTION( BUILD_CELERO_BENCHMARKS "Build Celero-based benchmarks. Note that Nonius and Celero benchmarks cannot be enabled at the same time currently." OFF )

if( BUILD_NONIUS_BENCHMARKS AND BUILD_CELERO_BENCHMARKS )
  message(SEND_ERROR "Nonius- and Celero-based benchmarks cannot be built at
  the some time. This will cause compile errors.")
endif()

#  ____      _       _       _   _       _ _     _____         _       
# / ___|__ _| |_ ___| |__   | | | |_ __ (_) |_  |_   _|__  ___| |_ ___ 
#| |   / _` | __/ __| '_ \  | | | | '_ \| | __|   | |/ _ \/ __| __/ __|
#| |__| (_| | || (__| | | | | |_| | | | | | |_    | |  __/\__ \ |_\__ \
# \____\__,_|\__\___|_| |_|  \___/|_| |_|_|\__|   |_|\___||___/\__|___/
                                                                      



include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/include" )
include_directories( "${PROJECT_SOURCE_DIR}/src" )

if(BUILD_UNIT_TESTS)

set( binDir ${CMAKE_BINARY_DIR}/testBin )
set( buildDir ${CMAKE_BINARY_DIR} )

# Unit tests written with the Catch framework are
# automatically detected. Just put the source files
# in the `CatchTests` directory.

set( unitTest_DIR ./CatchTests)
# Search for sources.
file( GLOB_RECURSE unitTest_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "${unitTest_DIR}/*.cpp" )
message(STATUS "Detected Catch-based Unit Test Sources:")
foreach(testSrc ${unitTest_SOURCES})
  message(STATUS "  ${testSrc}" )
endforeach()


# Add unit test executable
set( testName CatchTests )
add_executable(${testName} ${unitTest_SOURCES})
target_compile_definitions( ${testName} PRIVATE -DTESTING )

# Link dependencies
target_link_libraries(${testName} ) # If your project builds a library, link to it here

# Make sure that testing binaries are created in the binDir directory
set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${binDir} CXX_STANDARD 11)

# Now add the unit test binary to the CMake test target
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME ${testName} 
         WORKING_DIRECTORY ${buildDir}
         COMMAND ${binDir}/${testName} )


# Copy the data directory to the build directory
file(COPY data DESTINATION ${buildDir} )

endif()

# ____                  _                          _        
#| __ )  ___ _ __   ___| |__  _ __ ___   __ _ _ __| | _____ 
#|  _ \ / _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ / __|
#| |_) |  __/ | | | (__| | | | | | | | | (_| | |  |   <\__ \
#|____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\___/
                                                           


set( binDir ${CMAKE_BINARY_DIR}/benchmarkBin )

# _   _             _             ____                  _                          _        
#| \ | | ___  _ __ (_)_   _ ___  | __ )  ___ _ __   ___| |__  _ __ ___   __ _ _ __| | _____ 
#|  \| |/ _ \| '_ \| | | | / __| |  _ \ / _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ / __|
#| |\  | (_) | | | | | |_| \__ \ | |_) |  __/ | | | (__| | | | | | | | | (_| | |  |   <\__ \
#|_| \_|\___/|_| |_|_|\__,_|___/ |____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\___/
                                                                                           


# Benchmarks written with the Nonius framework are
# automatically detected. Just put the source files
# in the `NoniusBenchmarks` directory.

if(BUILD_NONIUS_BENCHMARKS)
set( benchmark_DIR ./NoniusBenchmarks)
# Search for sources.
file( GLOB_RECURSE benchmark_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "${benchmark_DIR}/*.cpp" )
    message(STATUS "Detected Nonius-based Benchmark Sources:")
foreach(benchmarkSrc ${benchmark_SOURCES})
  message(STATUS "  ${benchmarkSrc}" )
endforeach()

# Nonius requires a threading library
find_package(Threads)

if( NOT Threads_FOUND )
  message(SEND_ERROR "Nonius requires a threading library, which was not found. Either install a Threads library, or disable Nonius benchmarks." )
endif()

if(Threads_FOUND)

foreach(benchmarkSrc ${benchmark_SOURCES})

  get_filename_component(benchmarkName ${benchmarkSrc} NAME_WE)

  # Add benchmark executable
  add_executable(${benchmarkName} ${benchmarkSrc})

  # Link dependencies
  target_link_libraries(${benchmarkName} Threads::Threads ) # If your project builds a library, link to it here

  # Make sure that testing binaries are created in the binDir directory and use C++11
  set_target_properties(${benchmarkName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${binDir} CXX_STANDARD 11)

endforeach()
endif()
endif()


#  ____     _                  ____                  _                          _        
# / ___|___| | ___ _ __ ___   | __ )  ___ _ __   ___| |__  _ __ ___   __ _ _ __| | _____ 
#| |   / _ \ |/ _ \ '__/ _ \  |  _ \ / _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ / __|
#| |__|  __/ |  __/ | | (_) | | |_) |  __/ | | | (__| | | | | | | | | (_| | |  |   <\__ \
# \____\___|_|\___|_|  \___/  |____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\___/

if(BUILD_CELERO_BENCHMARKS)

set( benchmark_DIR ${CMAKE_CURRENT_SOURCE_DIR}/CeleroBenchmarks)

if( EXISTS ${benchmark_DIR} )

set(CELERO_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/externals/Celero CACHE PATH "")
set(CELERO_INCLUDE_DIR ${CELERO_DIRECTORY}/include CACHE PATH "Celero's include directory.")
set(CELERO_LIBRARY celero CACHE STRING "")

if( NOT EXISTS ${CELERO_DIRECTORY} )
  message(SEND_ERROR "Celero directory not found. Please clone https://github.com/DigitalInBlue/Celero into ${CELERO_DIRECTORY} or disable Celero benchmarks." )
endif()


#If celero is "out of tree", we have to add the subdirectory differently.
add_subdirectory(${CELERO_DIRECTORY} ${CMAKE_CURRENT_BINARY_DIR}/celero)
include_directories(${CELERO_INCLUDE_DIR})


# Benchmarks written with the Celero framework are
# automatically detected. Just put the source files
# in the `CeleroBenchmarks` directory.

# Search for sources.
file( GLOB_RECURSE benchmark_SOURCES
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "${benchmark_DIR}/*.cpp" )
    message(STATUS "Detected Celero-based Benchmark Sources:")
foreach(benchmarkSrc ${benchmark_SOURCES})
  message(STATUS "  ${benchmarkSrc}" )
endforeach()


foreach(benchmarkSrc ${benchmark_SOURCES})

get_filename_component(benchmarkName ${benchmarkSrc} NAME_WE)

# Add benchmark executable
add_executable(${benchmarkName} ${benchmarkSrc})

# Link dependencies
target_link_libraries(${benchmarkName} celero ${CELERO_EXTRA_LIBS} ) # If your project builds a library, link to it here

# Make sure that testing binaries are created in the binDir directory and use C++11
set_target_properties(${benchmarkName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${binDir} CXX_STANDARD 11)

endforeach()
endif()

endif()
