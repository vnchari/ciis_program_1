cmake_minimum_required(VERSION 3.22)
project(programming_assigment_two)

set(CMAKE_CXX_STANDARD 20)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set( CMAKE_CXX_FLAGS  "-O0 -Wall -Wextra")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set( CMAKE_CXX_FLAGS  "-O3 -Wall -Wextra")
endif()
find_package (Eigen3 3.4 REQUIRED NO_MODULE)
add_executable(programming_assigment_two main.cpp)
target_link_libraries (programming_assigment_two Eigen3::Eigen)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_custom_command(TARGET programming_assigment_two POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:programming_assigment_two> ${CMAKE_CURRENT_BINARY_DIR}/..)
endif()
