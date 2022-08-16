
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

file(COPY ${CMAKE_BINARY_DIR}/../Input/harmonic_oscillator_input.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Input/harmonic_oscillator_input.in )

#  SET( legacy_link   ./Input/harmonic_oscillator_input.in)
#  SET( legacy_target ${CMAKE_BUILD_DIRECTORY}/../Input/harmonic_oscillator_input.in)
#  add_custom_target(symlink ALL
#     ${CMAKE_COMMAND} -E create_symlink ${legacy_target} ${legacy_link}
#                      )
