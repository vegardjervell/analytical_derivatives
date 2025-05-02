cmake_minimum_required(VERSION 3.24)
project(polynomials LANGUAGES CXX Fortran)

# set(CMAKE_FORTRAN_COMPILER /opt/homebrew/bin/gfortran)
# set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "The OSX architecture")

set(CMAKE_CXX_FLAGS "-O3 -Wfatal-errors")

add_executable(fpoly polynomials.f90 polynomial_demo.f90)
add_executable(cpoly polynomials.cpp factorial.cpp polynomials_demo.cpp)
target_include_directories(cpoly PUBLIC ${CMAKE_CURRENT_LIST_DIR})