cmake_minimum_required( VERSION 3.1 )
message( "Configuring unit test build." )



#  ____      _       _       _   _       _ _     _____         _       
# / ___|__ _| |_ ___| |__   | | | |_ __ (_) |_  |_   _|__  ___| |_ ___ 
#| |   / _` | __/ __| '_ \  | | | | '_ \| | __|   | |/ _ \/ __| __/ __|
#| |__| (_| | || (__| | | | | |_| | | | | | |_    | |  __/\__ \ |_\__ \
# \____\__,_|\__\___|_| |_|  \___/|_| |_|_|\__|   |_|\___||___/\__|___/
                                                                      



include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/include" )
include_directories( "${PROJECT_SOURCE_DIR}/src" )

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

# _   _             _             ____                  _                          _        
#| \ | | ___  _ __ (_)_   _ ___  | __ )  ___ _ __   ___| |__  _ __ ___   __ _ _ __| | _____ 
#|  \| |/ _ \| '_ \| | | | / __| |  _ \ / _ \ '_ \ / __| '_ \| '_ ` _ \ / _` | '__| |/ / __|
#| |\  | (_) | | | | | |_| \__ \ | |_) |  __/ | | | (__| | | | | | | | | (_| | |  |   <\__ \
#|_| \_|\___/|_| |_|_|\__,_|___/ |____/ \___|_| |_|\___|_| |_|_| |_| |_|\__,_|_|  |_|\_\___/
                                                                                           

set( binDir ${CMAKE_BINARY_DIR}/benchmarkBin )

# Benchmarks written with the Nonius framework are
# automatically detected. Just put the source files
# in the `NoniusBenchmarks` directory.

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

if(Threads_FOUND)

# Add benchmark executable
set( benchmarkName NoniusBenchmarks )
add_executable(${benchmarkName} ${benchmark_SOURCES})
#target_compile_definitions( ${benchmarkName} )

# Link dependencies
target_link_libraries(${benchmarkName} Threads::Threads ) # If your project builds a library, link to it here

# Make sure that testing binaries are created in the binDir directory and use C++11
set_target_properties(${benchmarkName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${binDir} CXX_STANDARD 11)


endif()

