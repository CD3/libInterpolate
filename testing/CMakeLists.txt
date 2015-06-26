cmake_minimum_required( VERSION 2.8 )



find_package( Boost COMPONENTS unit_test_framework )
add_definitions( -DBOOST_TEST_DYN_LINK -DBOOSTTEST_MAIN )
if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  include_directories( ${Boost_INCLUDE_DIRS} )
  set( LIBS ${LIBS} ${Boost_LIBRARIES} )
  set( HAVE_UNITTESTPP 1 )  # compatible with autoconf

  include_directories( "${PROJECT_SOURCE_DIR}/interpLib" )
  include_directories( "${PROJECT_SOURCE_DIR}/interpLib/fileReader" )
     link_directories( "${PROJECT_SOURCE_DIR}/src/build" )

  file( GLOB_RECURSE unitTest_SOURCES
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        *.cpp )
  message(STATUS "Detected Unit Test Sources:")
  foreach(testSrc ${unitTest_SOURCES})
    message(STATUS "  ${testSrc}" )
  endforeach()
  
  set( testDir ${CMAKE_BINARY_DIR}/testBin )
  foreach(testSrc ${unitTest_SOURCES})
    #Extract the filename without an extension (NAME_WE)
    get_filename_component(testName ${testSrc} NAME_WE)

    #Add compile target
    add_executable(${testName} ${testSrc})

    #Link to Boost libraries AND your targets and dependencies
    target_link_libraries(${testName} ${EXTRA_LIBS} ${LIBS} interpLib boost_system boost_filesystem)

    #Move testing binaries into a testDIr directory
    set_target_properties(${testName} PROPERTIES 
        RUNTIME_OUTPUT_DIRECTORY  ${testDir} )

    #Finally add it to test execution - 
    #Notice the WORKING_DIRECTORY and COMMAND
    add_test(NAME ${testName} 
             WORKING_DIRECTORY ${testDir}
             COMMAND ${testDir}/${testName} )
  endforeach(testSrc)
 

endif(Boost_UNIT_TEST_FRAMEWORK_FOUND)
