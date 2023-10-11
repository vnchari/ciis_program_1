cmake_minimum_required(VERSION 3.22)
project(programming_assigment_one)

set(CMAKE_CXX_STANDARD 20)
set( CMAKE_CXX_FLAGS  "-Wall -Wextra")

add_executable(programming_assigment_one main.cpp
        geometry/cal_body_parser.h
        geometry/cal_readings_parser.h
        geometry/em_pivot_parser.h
        geometry/opt_pivot_parser.h
        geometry/output_file_creator.cpp
        geometry/output_file_creator.h
        geometry/pivot_calibration.cpp
        geometry/pivot_calibration.h
)
