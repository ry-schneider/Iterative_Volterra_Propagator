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

  type(complex_sb_mat)       :: h_grid, h_zero, v
  complex(8), allocatable    :: psi(:)
  complex(8), allocatable    :: h_diagonal(:), h_offdiagonal(:,:)
  complex(8), allocatable    :: v_diagonal(:), v_offdiagonal(:,:)
  complex(8), allocatable    :: z_el(:)
  real(8), allocatable       :: prj_energies(:)
  real(8), allocatable       :: coulomb_waves(:,:), GC(:), FC(:), FCP(:), GCP(:)
  real(8), allocatable       :: r_grid(:), spectra(:)
  real(8)                    :: t
  integer                    :: i,j,k,d,p,q
  integer                    :: max_iter, ifail

  procedure(pulse_at_t_func), pointer  :: pulse
  
  call hydrogen_read

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

  max_iter = 0

  pulse => select_pulse_type(pulse_name)
  d = ceiling(r_max/dr) ! size of radial grid

  ! allocate arrays
  allocate(r_grid(1:d), psi(1:l_max*d), h_diagonal(1:l_max*d), h_offdiagonal(1:l_max*d-1,1:d), v_diagonal(1:l_max*d), &
       v_offdiagonal(1:l_max*d-1, 1:d), prj_energies(1:coul_num), coulomb_waves(1:d,1:l_max), GC(1:l_max), &
       FCP(1:l_max), GCP(1:l_max), z_el(1:l_max), spectra(1:coul_num))

  ! construct the radial grid (equally spaced on [dr, r_max])
  do i = 1,d
     r_grid(i) = i*dr
  end do

  ! construct time independent Hamiltonian
  call h_grid%make_banded_matrix_on_grid(r_grid, dr, band_num_sym_mat) ! construct finite difference approx and potential on the grid
  h_offdiagonal(:,:) = 0
  do j = 1,l_max
     do k = 1,d
        h_diagonal((j-1)*d+k) = h_grid%diagonal(k) + (l_max-1)*l_max/(2d0*r_grid(k)*r_grid(k)) - 1d0/r_grid(k)
     end do
     do p = 1,band_num_sym_mat
        h_offdiagonal((j-1)*d+1:j*d-p,p) = h_grid%offdiagonal(:,p)
     end do
  end do
  call h_zero%initialize(l_max*d, d, h_diagonal, h_offdiagonal)

  ! consruct V (the fixed matrix multiplying E(t) to give the time dependent piece)
  v_diagonal(:) = 0
  v_offdiagonal(:,:) = 0
  do i = 1,l_max-1
     v_offdiagonal((i-1)*d+1:i*d,d) = ((1d0/3d0)*sqrt(dble((2*i+1)*(2*i+3)))*F_3J(i,0,i+1,0,1,0,.true.)**2)*r_grid(:) ! F_3J is a Clebsch-Grodran coeff
  end do
  call v%initialize(l_max*d, d, v_diagonal, v_offdiagonal)

  ! construction initial wavefunction (starts in 1s)
  psi(:) = 0
  do i = 1,d
     psi(i) = 2d0*exp(-r_grid(i))
  end do

  open(unit=87, file=datafilename)

  ! start propagation
  call begin_timing
  t = 0

  do while (t < t_intv)
     if (soln_method /= 'itvolt') then
        stop

     ! call general iteration routines   
     else if (it_type == 'gmres') then
        call linear_solve(h_zero, t, psi, v, max_iter)
     else
        call iterative_loop(h_zero, t, psi, v, max_iter)
     end if

     ! error analysis at time t

     ! move to next step
     t = t + dt
  end do
  
  ! choose continuum energies to project  (equal to k^2/2)
  do i = 1,coul_num
     prj_energies(i) = 0d0 ! FIGURE OUT HOW TO CHOOSE/EVALUATE THESE
  end do

  ! compute photoelectron spectra by projecting on Coulomb waves
  spectra(:) = 0
  do i = 1,coul_num
     ! evaluate Coulomb functions on the grid
     do j = 1,d
        call coul90(r_grid(j), prj_energies(i), 0d0, l_max-1, coulomb_waves(j,:), GC, FCP, GCP, 0, ifail)
     end do

     ! integrate produce of psi and Coulomb functions via the composite trapezoidal rule (over radial grid) and sum
     z_el(:) = 0
     do k = 1,l_max
        do p = 1,d-1
           z_el(k) = z_el(k) + dr*coulomb_waves(p,k)*psi((k-1)*d+p)
        end do
        z_el(k) = z_el(k) + 0.5d0*dr*coulomb_waves(d,k)*psi(k*d)
        spectra(i) = spectra(i) + abs(z_el(k))**2
     end do
  end do

  
  close(87)
  

end program hydrogen
