cmake_minimum_required(VERSION 3.16)

find_package(GSL)
find_package(OpenMP REQUIRED)
find_package(Catch2 REQUIRED)
# Find all sources.
file(
  GLOB_RECURSE unitTest_SOURCES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "./CatchTests/*.cpp")
message(STATUS "Detected Catch-based Unit Test Sources:")
foreach(testSrc ${unitTest_SOURCES})
  message(STATUS "  ${testSrc}")
endforeach()

set(TEST_WITH_WERROR ON CACHE BOOL "Enable -Werror for tests")

set(testName ${PROJECT_NAME}_CatchTests)
add_executable(${testName} ${unitTest_SOURCES})
target_compile_definitions(
  ${testName} PRIVATE TESTING $<$<TARGET_EXISTS:GSL::gsl>:HAVE_GSL>)
target_link_libraries(${testName} Interpolate
  $<$<TARGET_EXISTS:GSL::gsl>:GSL::gsl> OpenMP::OpenMP_CXX Catch2::Catch2WithMain)
set_target_properties(${testName} PROPERTIES CXX_STANDARD 17)
if(MSVC)
  target_compile_options(${testName} PRIVATE /W4 -wd4996)
else()
  target_compile_options(
    ${testName}
    PRIVATE
      -Wall
      -Wextra
      -pedantic
      -Wno-unknown-pragmas
      -Wno-unused-parameter
      -Wno-gnu-zero-variadic-macro-arguments)
  if(TEST_WITH_WERROR)
    target_compile_options(${testName} PRIVATE -Werror)
  endif()
endif()

add_test(
  NAME ${testName}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${testName})

file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
