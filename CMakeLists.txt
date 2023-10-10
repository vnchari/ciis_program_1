cmake_minimum_required(VERSION 3.22)
project(programming_assigment_one)

set(CMAKE_CXX_STANDARD 20)
set( CMAKE_CXX_FLAGS  "-Wall -Wextra")
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(programming_assigment_one main.cpp
        geometry/frame_lib.cpp)
target_link_libraries (programming_assigment_one Eigen3::Eigen)
