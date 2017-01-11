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
target_compile_features( ${testName} PRIVATE ${cxx_needed_features})
target_compile_definitions( ${testName} PRIVATE -DTESTING )

# Link dependencies
target_link_libraries(${testName} ) # If your project builds a library, link to it here

# Make sure that testing binaries are created in the binDir directory
set_target_properties(${testName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY  ${binDir} )

# Now add the unit test binary to the CMake test target
#Notice the WORKING_DIRECTORY and COMMAND
add_test(NAME ${testName} 
         WORKING_DIRECTORY ${buildDir}
         COMMAND ${binDir}/${testName} )


# Copy the data directory to the build directory
file(COPY data DESTINATION ${buildDir} )


#Copy config files to build directory
add_custom_command( OUTPUT do_copyPythonTests
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    "${CMAKE_SOURCE_DIR}/testing/PythonTests"
                    "${CMAKE_BINARY_DIR}"
                  )

add_custom_target(pythonTests ALL DEPENDS do_copyPythonTests)
