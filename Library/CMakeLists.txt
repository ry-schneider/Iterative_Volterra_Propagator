
add_library(library STATIC
            class.f
            convpt.f
            gamfun.f
            Gauss_Quadrature.f90
            gaussq.f
            gbtql2.f
            gbslve.f
            general_utility.f90
            hermite.f90
            Lagrange_weights.f90
            Lagrange.f90
            precin.f90
            timing.f90)

set_target_properties(library PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
