cmake_minimum_required(VERSION 3.8)
project(jacobi)

set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Ofast -fopenmp")

add_compile_options(-DPARALLEL_REDUCE)

include_directories("~/fastflow")

set(SOURCE_FILES src/main.cpp src/jacobi.hpp src/barrier.hpp src/jacobi_fastflow.hpp
        src/jacobi_thread.hpp src/jacobi_thread.cpp src/utils.hpp src/jacobi_omp.hpp)

add_executable(jacobi ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(
        jacobi
        pthread
)