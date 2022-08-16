
add_executable(harmonic_oscil.exe 
                harmonic_oscillator.f90
                banded_matrices.f90
                fconfig.f90
                grid.f90
                integral_method.f90
                parameter_read.f90
                parameters.f90
                potential.f90
                propagator.f90
                pulse.f90           )

target_link_libraries(harmonic_oscil.exe library ${lapackblas_lib})                

target_include_directories(harmonic_oscil.exe PUBLIC $<TARGET_PROPERTY:library,Fortran_MODULE_DIRECTORY>)
