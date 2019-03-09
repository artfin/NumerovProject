cmake_minimum_required(VERSION 3.13)
project(spectra_example)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")

set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

set(SPECTRA_DIR "/home/artfin/Desktop/repos/spectra-0.7.0/include/")
include_directories(${SPECTRA_DIR})

add_executable(spectra_example main.cpp)