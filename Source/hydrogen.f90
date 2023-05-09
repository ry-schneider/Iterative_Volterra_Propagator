program hydrogen
  use parameters
  use parameter_read
  use general_utility
  use banded_matrices
  use propagator
  use pulse_module
  use integral_method
  use timing
  use clebsch_gordan
  implicit none

  type(complex_sb_mat)       :: h_zero, v
  complex(8), allocatable    :: psi(:)
  real(8), allocatable       :: r_grid(:)
  real(8)                    :: t
  integer                    :: i,d
  integer                    :: max_iter

  procedure(pulse_at_t_func), pointer  :: pulse

  call hydrogen_read
  call initial_print
  max_iter = 0

  pulse => select_pulse_type(pulse_name)
  d = r_size

  ! allocate arrays
  allocate(r_grid(1:d), psi(1:l_max*d))

  ! construct the radial grid (equally spaced on [dr, r_max])
  do i = 1,d
     r_grid(i) = i*dr
  end do

  call construct_hzero(h_zero, r_grid)
  call construct_v(v, r_grid)

  ! construction initial wavefunction (starts in 1s)
  psi(:) = 0
  do i = 1,d
     psi(i) = 2d0*exp(-r_grid(i))
  end do

  ! start propagation
  call begin_timing
  t = 0

  do while (t < t_intv)
     if (soln_method /= 'it_volt') then
        stop

     ! call general iteration routines   
     else if (it_type == 'gmres') then
        call linear_solve(h_zero, t, psi, v, max_iter)
     else
        call iterative_loop(h_zero, t, psi, v, max_iter)
     end if

     ! move to next step
     t = t + dt
     print *, 'Time: ', t
  end do

  call stop_timing
  print *, 'Propogation finished -- computing ionization probabilities!'
  print *, '************************************'

  call compute_ionization_probabilities(psi, r_grid)

  print *, 'Propagation Time:'
  call print_timing
  

contains

  
  subroutine initial_print
    
    print *, '*************************************'
    print *, 'Problem: 3D Hydrogen'
    print *, 'Pulse: ', pulse_name
    print *, 'Potential: ', potential_type
    print *, '*************************************'
    print *, 'Expansion paramters:'
    print *, 'Radial grid max:', r_max
    print *, 'Radial grid step size: ', dr
    print *, 'Number of angular term:', l_max
    print *, '*************************************'
    print *, 'Solution details:'
    if (soln_method == 'it_volt') then
       if (it_type == 'short_time') then
          print *, 'Solution method: Short Time Approximation'
       else
          print *, 'Solution method: Iterative Volterra Propagator (ITVOLT)'
          print *, 'Iteration type: ', it_type
          print *, 'Quadrature type: ', quad_type
          print *, 'Number of quadrature points:', quad_pt
       end if
       print *, 'Method for computing exponentials: ', prop_method
    else
       print *, 'Solution method: ', soln_method
    end if
    print *, 'Final propagation time:', t_intv
    print *, 'Propagation step size:', dt
    print *, '************************************'
    
  end subroutine initial_print
  

  subroutine construct_hzero(h_zero, r_grid)
    type(complex_sb_mat)       :: h_zero, h_grid
    real(8)                    :: r_grid(:)
    complex(8), allocatable    :: h_diagonal(:), h_offdiagonal(:,:)
    integer                    :: j, k, p

    allocate(h_diagonal(1:l_max*r_size), h_offdiagonal(1:l_max*r_size-1,1:band_num_sym_mat))

    call h_grid%make_banded_matrix_on_grid(r_grid, dr, band_num_sym_mat) ! construct finite difference approx and potential on the grid
    h_offdiagonal(:,:) = 0
    do j = 1,l_max
       do k = 1,r_size
          h_diagonal((j-1)*r_size+k) = h_grid%diagonal(k) + (j-1)*j/(2d0*r_grid(k)*r_grid(k)) - 1d0/r_grid(k)
       end do
       do p = 1,band_num_sym_mat
          h_offdiagonal((j-1)*r_size+1:j*r_size-p,p) = h_grid%offdiagonal(1:r_size-p,p)
       end do
    end do
    
    call h_zero%initialize(l_max*r_size, band_num_sym_mat, h_diagonal, h_offdiagonal)
    
  end subroutine construct_hzero


  subroutine construct_v(v, r_grid)
    type(complex_sb_mat)       :: v
    real(8)                    :: r_grid(:)
    complex(8), allocatable    :: v_diagonal(:), v_offdiagonal(:,:)
    integer                    :: i

    allocate(v_diagonal(1:l_max*r_size), v_offdiagonal(1:l_max*r_size-1, 1))

    v_diagonal(:) = 0
    v_offdiagonal(:,:) = 0
    do i = 1,l_max-1
       v_offdiagonal((i-1)*r_size+1:i*r_size,1) = ((1d0/3d0)*sqrt(dble((2*i-1)*(2*i+1)))*F_3J(i-1,0,i,0,1,0,.true.)**2)*r_grid(:) ! F_3J is a Clebsch-Grodran coeff
    end do
    
    call v%initialize(l_max*r_size, 1, v_diagonal, v_offdiagonal)
    
  end subroutine construct_v


  subroutine compute_ionization_probabilities(psi, r_grid)
    complex(8)                 :: psi(:)
    real(8)                    :: r_grid(:)
    real(8), allocatable       :: momenta(:), spectra(:)
    real(8), allocatable       :: coulomb_waves(:,:), GC(:), FC(:), FCP(:), GCP(:)
    complex(8), allocatable    :: z_el(:)
    integer                    :: i, j, k, p, ifail

    allocate(momenta(1:coul_num), coulomb_waves(1:r_size,1:l_max), GC(1:l_max), FCP(1:l_max), &
         GCP(1:l_max), z_el(1:l_max), spectra(1:coul_num))

    open(unit=87, file='hydrogen_spectra')
    write(87,*) 'Energy', 'Angle-integrated Probability'

    do i = 1,coul_num
       momenta(i) = sqrt(2*(0.0001d0 + (i-1)*0.002d0)) ! momenta match energies output by Nico's code
    end do

    ! compute photoelectron spectra by projecting on Coulomb waves
    spectra(:) = 0
    do i = 1,coul_num
       ! evaluate Coulomb functions on the grid
       do j = 1,r_size
          call coul90(momenta(i)*r_grid(j), -1d0/momenta(i), 0d0, l_max-1, coulomb_waves(j,:), GC, FCP, GCP, 0, ifail)
       end do

       ! integrate produce of psi and Coulomb functions via the composite trapezoidal rule (over radial grid) and sum
       z_el(:) = 0
       do k = 1,l_max
          do p = 1,r_size-1
             z_el(k) = z_el(k) + sqrt(2/pi/momenta(i))*dr*coulomb_waves(p,k)*psi((k-1)*r_size+p)
          end do
          z_el(k) = z_el(k) + 0.5d0*sqrt(2/pi/momenta(i))*dr*coulomb_waves(r_size,k)*psi(k*r_size)
          spectra(i) = spectra(i) + abs(z_el(k))**2
       end do

       write(87,*) 0.5d0*momenta(i)**2, spectra(i)
    end do

    close(87)

    print *, 'Total Ionization Probability: ', sum(spectra)
    
  end subroutine compute_ionization_probabilities

end program hydrogen
