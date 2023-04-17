cmake_minimum_required(VERSION 3.26)
project(lcodes)

set(CMAKE_CXX_STANDARD 23)

add_library(lcodes lcodes.cpp channel.h channel.cpp)
target_include_directories(lcodes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(lcodes-demo main.cpp channel.cpp)
target_link_libraries(lcodes-demo PRIVATE lcodes)
