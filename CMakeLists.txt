cmake_minimum_required(VERSION 3.0)
project(jacobi)

set(CMAKE_CXX_STANDARD 11)

#Please set fastflow directory
set(ff "~/fastflow" CACHE STRING "fastflow directory")

set(CMAKE_CXX_FLAGS "-O3")

#ON enables OpenMP
option(omp "enable OpenMP" OFF)
#ON enable xeon phi flags
option(mic "Xeon phi compilation" OFF)
#ON enable intel compiler
option(intel "enable intel compiler" OFF)

if(intel)
    set(CMAKE_CXX_COMPILER "icc")
else()
    set(CMAKE_CXX_COMPILER "g++")
endif()

if(omp)
    add_compile_options("-DWITHOMP")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -openmp-link=static")
endif()

if(mic)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic -DNO_DEFAULT_MAPPING")
endif()

include_directories(${ff})

set(SOURCE_FILES src/main.cpp src/jacobi.hpp src/barrier.hpp src/jacobi_fastflow.hpp
        src/jacobi_thread.hpp src/jacobi_thread.cpp src/utils.hpp src/jacobi_omp.hpp)

add_executable(jacobi ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(
        jacobi
        pthread
)
