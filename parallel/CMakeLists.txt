cmake_minimum_required(VERSION 3.5)
project(parallel)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS_DEBUG "-fopenmp -g3 -O0 -Wall -Werror")

set(SOURCE_FILES test.cpp)

add_executable(parallel ${SOURCE_FILES})
target_link_libraries(parallel -lpthread)