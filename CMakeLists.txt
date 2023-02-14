cmake_minimum_required(VERSION 3.5)

project(mandelbrot LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

# -fopenmp flag (enables OpenMP)
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -O3 -funroll-loops -fopenmp")

add_executable(mandelbrot main.cpp
    complex.h complex.cpp
    colour.h colour.cpp)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(mandelbrot PUBLIC MPI::MPI_CXX)
