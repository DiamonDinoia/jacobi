cmake_minimum_required(VERSION 3.8)
project(jacobi)

unset(DEBUG)

set(CMAKE_CXX_COMPILER "g++")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-O3 -static")

include_directories("~/fastflow")

set(SOURCE_FILES src/main.cpp src/jacobi.hpp src/barrier.hpp src/jacobi_map.cpp src/jacobi_map.hpp src/jacobi_par_for.hpp
        src/jacobi_thread.hpp src/jacobi_thread.cpp src/utils.hpp)

add_executable(jacobi ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(
        jacobi
        pthread
)