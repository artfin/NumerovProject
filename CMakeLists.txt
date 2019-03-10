cmake_minimum_required(VERSION 3.5)
project(generalized_matrix_numerov)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")

# Eigen library
set(EIGEN_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN_DIR})

# ARPACK-C++ library (SPECTRA)
set(SPECTRA_DIR "/home/artfin/Desktop/repos/spectra-0.7.0/include/")
include_directories(${SPECTRA_DIR})

# GSL library
find_package(GSL REQUIRED)

# GMN := generalized matrix numerov
file( GLOB LIB_SOURCES lib/*.cpp )
file( GLOB LIB_HEADERS lib/*.h )
add_library( GMN ${LIB_SOURCES} ${LIB_HEADERS} )
install(TARGETS GMN
        ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/lib
        LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/lib)
