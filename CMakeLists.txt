cmake_minimum_required(VERSION 3.13)

project(Iterative_Voltera_Propagator)
enable_language(Fortran)
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/Mod )
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-cpp -g -fbacktrace -fcheck=all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -Ofast -march=native")
    set(lapackblas_lib "-L/usr/lib/ -llapack -lblas")
endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    #find_package(mkl REQUIRED PATHS /opt/intel/oneapi/mkl/2021.3.0/latest/lib)
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-cpp -g -traceback -check all -debug all")
    set(CMAKE_Fortran_FLAGS_RELEASE "-cpp -Ofast")
    set(lapackblas_lib "-qmkl")
endif()


add_subdirectory(Library)
add_subdirectory(Source)
