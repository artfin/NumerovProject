cmake_minimum_required(VERSION 3.13)
project(brent)

set(CMAKE_CXX_STANDARD 11)

add_executable(brent main.cpp)

find_package(GSL REQUIRED)
target_link_libraries(brent GSL::gsl GSL::gslcblas)

