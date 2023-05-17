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

  type(complex_sb_mat)       :: h_zero, v, h_grid, h_mid
  type(banded_sym_mat)       :: angular_v, v_mid
  complex(8), allocatable    :: psi(:), alternate_psi(:,:)
  real(8), allocatable       :: r_grid(:)
  real(8)                    :: t
  integer                    :: i, j, r_reach, time_step
  integer                    :: max_iter

  procedure(pulse_at_t_func), pointer  :: pulse

  call hydrogen_read
  call initial_print
  max_iter = 0

  pulse => select_pulse_type(pulse_name)

  ! allocate arrays
  allocate(r_grid(1:r_size), alternate_psi(1:r_size,1:l_max))

  ! construct the radial grid (equally spaced on [dr, r_max])
  do i = 1,r_size
     r_grid(i) = i*dr
  end do

  ! construction initial wavefunction (starts in 1s)
  alternate_psi(:,:) = 0d0
  do i = 1,r_size
     alternate_psi(i,1) = 2d0*r_grid(i)*exp(-r_grid(i))
  end do

  ! set initial dynamic grid size
  do i = r_size,1,-2
     if (abs(alternate_psi(i,1)) > 1.d-9) then
        r_reach = i
        exit
     end if
  end do
  r_reach = min(r_reach+100, r_size)
  alternate_psi(r_reach+1:r_size,1) = 0d0

  ! construct radial and angular operators
  call h_grid%make_banded_matrix_on_grid(r_grid, dr, band_num_sym_mat)
  do i = 1,r_size
     h_grid%diagonal(i) = h_grid%diagonal(i) - 1d0/r_grid(i)
  end do
  call construct_angularv(angular_v)
  call v_mid%initialize(l_max, 1, angular_v%diagonal, angular_v%offdiagonal)
  call h_mid%initialize(r_size, band_num_sym_mat, h_grid%diagonal, h_grid%offdiagonal)
  
  ! start propagation
  call begin_timing
  time_step = 1
  t = 0

  do while (t < t_intv)
     ! adjust grid size based on current wavefunction
     if (time_step > 1) then
        do i = 1,l_max
           if (abs(alternate_psi(r_reach-100,i)) > grid_tol) then
              r_reach = r_reach + 100
              exit
           end if
        end do
        r_reach = min(r_reach, r_size)
     end if

     ! adjust midpoint operators
     h_mid%mat_size = r_reach
     h_mid%diagonal = h_grid%diagonal(1:r_reach)
     h_mid%offdiagonal = h_grid%offdiagonal(1:r_reach-1,1:h_grid%bsz)
     v_mid%offdiagonal =  pulse(t+0.5d0*dt)*angular_v%offdiagonal
     
     ! call split operator function (in propagator module)
     if (soln_method == 'split_operator') then
        alternate_psi(1:r_reach,:) = split_operator_cn(h_mid, v_mid, r_grid(1:r_reach), dt, alternate_psi(1:r_reach,:), 1)

     ! call general iteration routines (in integral_method module)   
     ! else if (it_type == 'gmres') then
     !    call linear_solve(h_zero, t, psi, v, max_iter)
     ! else
     !    call iterative_loop(h_zero, t, psi, v, max_iter)
        
     end if

     ! move to next step
     time_step = time_step + 1   
     t = t + dt
     print *, 'Time: ', t, 'Number of radial points: ', r_reach
  end do

  call stop_timing
  print *, 'Propogation finished -- computing ionization probabilities!'
  print *, '************************************'

  ! write final wave function in data file
  open(unit=88, file='final_wave_function.dat')
  do i = 1,r_reach
     write(88,*) alternate_psi(i,1), r_grid(i), 1
  end do
  close(88)

  ! record observables
  allocate(psi(1:r_reach*l_max))
  do j = 1,l_max
    psi(r_reach*(j-1)+1:j*r_reach) = alternate_psi(1:r_reach,j)
  end do

  call compute_ionization_probabilities(psi, r_grid(1:r_reach))

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


  subroutine construct_angularv(angular_v)
    type(banded_sym_mat)       :: angular_v
    real(8), allocatable       :: v_diagonal(:), v_offdiagonal(:,:)
    integer                    :: i

    allocate(v_diagonal(1:l_max), v_offdiagonal(1:l_max-1,1))

    v_diagonal(:) = 0
    do i = 1,l_max-1
       v_offdiagonal(i,1) = (1d0/3d0)*sqrt(dble((2*i-1)*(2*i+1)))*F_3J(i-1,0,i,0,1,0,.true.)**2
    end do

    call angular_v%initialize(l_max, 1, v_diagonal, v_offdiagonal)
    
  end subroutine construct_angularv


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
    real(8)                    :: total_prob
    integer                    :: i, j, k, p, d, ifail

    d = size(r_grid)

    allocate(momenta(1:coul_num), coulomb_waves(1:d,1:l_max), GC(1:l_max), FCP(1:l_max), &
         GCP(1:l_max), z_el(1:l_max), spectra(1:coul_num))

    open(unit=87, file='hydrogen_spectra')
    write(87,*) '         Energy      ', '     Angle-integrated Probability'

    do i = 1,coul_num
       momenta(i) = sqrt(2*(e_min + (i-1)*dE)) ! momenta match energies output by Nico's code
    end do

    ! compute photoelectron spectra by projecting on Coulomb waves
    spectra(:) = 0
    do i = 1,coul_num
       ! evaluate regular Coulomb functions on the grid
       do j = 1,d
          call coul90(momenta(i)*r_grid(j), -1d0/momenta(i), 0d0, l_max-1, coulomb_waves(j,:), GC, FCP, GCP, 0, ifail)
       end do

       ! integrate produce of psi and Coulomb functions via Simpson's rule (over radial grid) and sum
       z_el(:) = 0
       do k = 1,l_max
          if (mod(d,2) == 0) then
             do p = 1,d/2
                z_el(k) = z_el(k) + 4d0*coulomb_waves(2*p-1,k)*psi((k-1)*d+2*p-1)
             end do

             do p = 1,d/2-1
                z_el(k) = z_el(k) + 2d0*coulomb_waves(2*p,k)*psi((k-1)*d+2*p)
             end do

             z_el(k) = z_el(k) + coulomb_waves(d,k)*psi(k*d)
             z_el(k) = (1d0/3d0)*dr*sqrt(2d0/(pi*momenta(i)))*z_el(k)
             spectra(i) = spectra(i) + abs(z_el(k))**2
          else
             do p = 1,(d-1)/2
                z_el(k) = z_el(k) + 4d0*coulomb_waves(2*p-1,k)*psi((k-1)*d+2*p-1)
             end do

             do p = 1,(d-1)/2-1
                 z_el(k) = z_el(k) + 2d0*coulomb_waves(2*p,k)*psi((k-1)*d+2*p)
             end do

             z_el(k) = z_el(k) + coulomb_waves(d-1,k)*psi(k*d-1)
             z_el(k) = (1d0/3d0)*z_el(k)
             z_el(k) = z_el(k) + 0.5d0*coulomb_waves(d-1,k)*psi(k*d-1)
             z_el(k) = z_el(k) + 0.5d0*coulomb_waves(d,k)*psi(k*d)
             z_el(k) = dr*sqrt(2d0/(pi*momenta(i)))*z_el(k)
             spectra(i) = spectra(i) + abs(z_el(k))**2
             
          end if
          
          ! composite trapezoidal rule
          ! do p = 1,d-1
          !    z_el(k) = z_el(k) + sqrt(2d0/(pi*momenta(i)))*dr*coulomb_waves(p,k)*psi((k-1)*d+p)
          ! end do
          ! z_el(k) = z_el(k) + 0.5d0*sqrt(2d0/(pi*momenta(i)))*dr*coulomb_waves(d,k)*psi(k*d)
          ! spectra(i) = spectra(i) + abs(z_el(k))**2
          
       end do

       write(87,*) e_min+(i-1)*dE, spectra(i)
    end do

    close(87)

    ! integrate total ionization probability with Simpson's rule
    if (mod(coul_num,2) == 1) then
       total_prob = spectra(1)
       do i = 1,(coul_num-1)/2
          total_prob = total_prob + 4d0*spectra(2*i)
       end do

       do i = 1,(coul_num-1)/2-1
          total_prob = total_prob + 2d0*spectra(2*i+1)
       end do
       total_prob = total_prob + spectra(coul_num)
       total_prob = (1d0/3d0)*0.002*total_prob
       
    else
       total_prob = spectra(1)
       do i = 1,(coul_num-2)/2
          total_prob = total_prob + 4d0*spectra(2*i)
       end do

       do i = 1,(coul_num-2)/2-1
          total_prob = total_prob + 2d0*spectra(2*i+1)
       end do
       
       total_prob = total_prob + spectra(coul_num-1)
       total_prob = (1d0/3d0)*total_prob
       total_prob = total_prob + 0.5d0*spectra(coul_num-1)
       total_prob = total_prob + 0.5d0*spectra(coul_num)
       total_prob = 0.002*total_prob
  
    end if

    ! trapezoidal rule
    ! total_prob = 0.5d0*0.002d0*spectra(1)
    ! do i = 2,coul_num-1
    !    total_prob = total_prob + 0.002d0*spectra(i)
    ! end do
    ! total_prob = total_prob + 0.5d0*0.002d0*spectra(coul_num)
    
    print *, 'Total Ionization Probability: ', total_prob
    
  end subroutine compute_ionization_probabilities

end program hydrogen
