cmake_minimum_required(VERSION 3.5)

project(dunker_richardson)

# Eigen library
set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

# ARPACK-C++ library (SPECTRA)
set(SPECTRA_DIR "/home/artfin/Desktop/repos/spectra-0.7.0/include/")
include_directories(${SPECTRA_DIR})

# GSL library
find_package(GSL REQUIRED)

add_executable(dunker_richardson main.cpp dunker.cpp dunker.h dunker.cpp dunker.h)
target_link_libraries(dunker_richardson ${PROJECT_SOURCE_DIR}/../lib/libGMN.a)
target_link_libraries(dunker_richardson GSL::gsl GSL::gslcblas)
