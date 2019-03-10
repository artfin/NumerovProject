cmake_minimum_required(VERSION 3.13)
project(morse_classical_pf)

set(CMAKE_CXX_STANDARD 11)

set(HEP_DIR "/home/artfin/Desktop/libs/hep-mc-0.7/include")
include_directories(${HEP_DIR})

find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(morse_classical_pf
        main.cpp
        morse.cpp
        morse.h)

target_link_libraries(morse_classical_pf ${MPI_C_LIBRARIES})

