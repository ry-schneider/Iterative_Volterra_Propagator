
add_library (src_files SHARED banded_matrices.f90
                fconfig.f90
                grid.f90
                integral_method.f90
                parameter_read.f90
                parameters.f90
                potential.f90
                propagator.f90
                pulse_module.f90 )

add_executable(harmonic_oscil.exe 
                harmonic_oscillator.f90)

target_link_libraries(harmonic_oscil.exe src_files library ${lapackblas_lib})                

target_include_directories(harmonic_oscil.exe PUBLIC $<TARGET_PROPERTY:library,Fortran_MODULE_DIRECTORY>)

add_executable(model_ode.exe 
                model_ode.f90)

target_link_libraries(model_ode.exe src_files library ${lapackblas_lib})                

target_include_directories(model_ode.exe PUBLIC $<TARGET_PROPERTY:library,Fortran_MODULE_DIRECTORY>)

add_executable(two_channel.exe 
                two_channel.f90)

target_link_libraries(two_channel.exe src_files library ${lapackblas_lib})                

target_include_directories(two_channel.exe PUBLIC $<TARGET_PROPERTY:library,Fortran_MODULE_DIRECTORY>)

add_executable(two_level.exe 
                two_level_atom.f90)

target_link_libraries(two_level.exe src_files library ${lapackblas_lib})                

target_include_directories(two_level.exe PUBLIC $<TARGET_PROPERTY:library,Fortran_MODULE_DIRECTORY>)


# file(COPY ${CMAKE_BINARY_DIR}/../Input/harmonic_oscillator_input.in
#        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Input/harmonic_oscillator_input.in )

#  SET( legacy_link   ./Input/harmonic_oscillator_input.in)
#  SET( legacy_target ${CMAKE_BUILD_DIRECTORY}/../Input/harmonic_oscillator_input.in)
#  add_custom_target(symlink ALL
#     ${CMAKE_COMMAND} -E create_symlink ${legacy_target} ${legacy_link}
#                      )
