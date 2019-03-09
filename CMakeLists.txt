cmake_minimum_required(VERSION 3.5)
project(matrix_numerov)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")

set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

set(SPECTRA_DIR "/home/artfin/Desktop/repos/spectra-0.7.0/include/")
include_directories(${SPECTRA_DIR})

add_executable(matrix_numerov
        main.cpp
        constants.h
        filereader.cpp
        filereader.h
        generalizedmatrixnumerov.cpp
        generalizedmatrixnumerov.h
        matrixnumerov.cpp
        matrixnumerov.h
        matrixreader.cpp
        matrixreader.h
        parameters.cpp
        parameters.h
        eigenvalue.cpp eigenvalue.h)

find_package(GSL REQUIRED)
target_link_libraries(matrix_numerov GSL::gsl GSL::gslcblas)
