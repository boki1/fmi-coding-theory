cmake_minimum_required(VERSION 3.22)
project(gf-polynomials)

set(CMAKE_CXX_STANDARD 20)

# Enable conan.
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

option(GFP_BUILD_TESTS "Whether to build unit tests.")

add_library(gf-polynomials INTERFACE)
target_include_directories(gf-polynomials PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gh-polynomials PRIVATE
		CONAN_PKG::fmt)

if("${GFP_BUILD_TESTS}" STREQUAL "yes")
	message("-- gf-polynomials: Adding tests to build targets.")
	set(TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")
	include(GoogleTest)
	enable_testing()

	function(add_unit_test name source_file)
		message("-- gf-polynomials: Add ${name} unit test (from ${TESTS_DIR}/${source_file})")
		set(unit_test "test_${name}")
		add_executable(${unit_test} ${TESTS_DIR}/${source_file})
		set_target_properties(${unit_test} PROPERTIES LINKER_LANGUAGE CXX)
		target_link_libraries(${unit_test} PRIVATE
			CONAN_PKG::gtest
			gf-polynomials)
		gtest_discover_tests(${unit_test})
		add_test(NAME ${name}
			COMMAND ${CMAKE_CURRENT_BINARY_DIR}/bin/${unit_test})
	endfunction()

	add_unit_test(operations test_operations.cpp)
endif()
