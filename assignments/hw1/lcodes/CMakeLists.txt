cmake_minimum_required(VERSION 3.26)
project(lcodes)

set(CMAKE_CXX_STANDARD 20)

add_library(lcodes lcodes.h lcodes.cpp)
target_include_directories(lcodes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(lcodes-demo main.cpp)
target_link_libraries(lcodes-demo PRIVATE lcodes gsl gslcblas)
