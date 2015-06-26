cmake_minimum_required (VERSION 2.6)
project(int)


#not happy that I had to include yaml-cpp with a compiler flag like this.
#it would be nice if cmake would find yaml-cpp, find the version installed, and throw an error if the version is under 0.5
set(CMAKE_CXX_FLAGS "-std=c++0x ")


include_directories("${PROJECT_SOURCE_DIR}/interpLib")
include_directories("${PROJECT_SOURCE_DIR}/interpLib/fileReader")


set( EXTERNAL_LIBS boost_system boost_filesystem )
set( SOURCES
            main.cpp
            interpLib/interpLib.cpp)


add_library(interpLib interpLib/interpLib.cpp)
#add_library(beerLib src/lib/beerLib.cpp)


add_executable(int ${SOURCES})
target_link_libraries(int ${EXTERNAL_LIBS})


#find_package(yaml-cpp)
add_subdirectory( testing )
