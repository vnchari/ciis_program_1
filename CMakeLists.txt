cmake_minimum_required(VERSION 3.22)
project(programming_assigment_one)

set(CMAKE_CXX_STANDARD 20)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set( CMAKE_CXX_FLAGS  "-Og -Wall -Wextra")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set( CMAKE_CXX_FLAGS  "-march=native -O3 -Wall -Wextra")
endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable(programming_assigment_one main.cpp)
target_link_libraries (programming_assigment_one Eigen3::Eigen)
if (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_custom_command(TARGET programming_assigment_one POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:programming_assigment_one> ${CMAKE_CURRENT_BINARY_DIR}/..)
endif()
