cmake_minimum_required(VERSION 3.5)

project(HO_richardson_extrapolation)

# Eigen library
set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

# ARPACK-C++ library (SPECTRA)
set(SPECTRA_DIR "/home/artfin/Desktop/repos/spectra-0.7.0/include/")
include_directories(${SPECTRA_DIR})

# GSL library
find_package(GSL REQUIRED)

add_executable(HO_richardson_extrapolation HO_richardson_extrapolation.cpp)
target_link_libraries(HO_richardson_extrapolation ${PROJECT_SOURCE_DIR}/../lib/libGMN.a)
target_link_libraries(HO_richardson_extrapolation GSL::gsl GSL::gslcblas)