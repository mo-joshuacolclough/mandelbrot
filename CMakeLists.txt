cmake_minimum_required(VERSION 3.5)

project(mandelbrot LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_executable(mandelbrot main.cpp
    complex.h complex.cpp
    colour.h colour.cpp)

target_link_libraries(mandelbrot Threads::Threads)
