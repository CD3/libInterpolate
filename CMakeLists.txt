cmake_minimum_required (VERSION 2.6)
project(interpolationExample)


set(CMAKE_CXX_FLAGS "-std=c++0x ")


include_directories("${PROJECT_SOURCE_DIR}/interpLib")


set( EXTERNAL_LIBS boost_system boost_filesystem )
set( SOURCES
            example/main.cpp )


add_executable(interpolateExample ${SOURCES})
target_link_libraries(interpolateExample ${EXTERNAL_LIBS})


add_subdirectory( testing )
