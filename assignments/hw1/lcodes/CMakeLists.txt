cmake_minimum_required(VERSION 3.22)
project(lcodes)

set(CMAKE_CXX_STANDARD 20)

# Enable conan.
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

option(LCODES_BUILD_DEMO "Whether the additional demonstration program should be built.")

add_library(lcodes lcodes.h lcodes.cpp gsl_wrapper.h gsl_wrapper.cpp channel.cpp channel.h)
target_include_directories(lcodes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(lcodes PRIVATE
		CONAN_PKG::gsl
		CONAN_PKG::ms-gsl
		CONAN_PKG::fmt)

#if("${LCODES_BUILD_DEMO}" STREQUAL "yes")
	add_executable(lcodes-demo main.cpp)
	target_link_libraries(lcodes-demo PRIVATE lcodes)
#endif()
