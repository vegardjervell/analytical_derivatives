cmake_minimum_required(VERSION 3.24)
project(polynomials LANGUAGES CXX Fortran)

# set(CMAKE_FORTRAN_COMPILER /opt/homebrew/bin/gfortran)
# set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "The OSX architecture")

set(CMAKE_CXX_FLAGS "-O3 -Wfatal-errors")
set(CMAKE_Fortran_FLAGS "-O3 -fdefault-real-8")

add_executable(fpoly polynomials.f90 polynomial_demo.f90)
add_executable(cpoly polynomials.cpp factorial.cpp polynomials_demo.cpp)
target_include_directories(cpoly PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)
pybind11_add_module(libpyanalytical polynomials.cpp factorial.cpp polynomial_pybind.cpp)
target_include_directories(libpyanalytical PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS libpyanalytical DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})