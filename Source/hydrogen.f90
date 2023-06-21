program hydrogen
  use parameters
  use parameter_read
  use general_utility
  use banded_matrices
  use propagator
  use pulse_module
  use potential
  use integral_method
  use timing
  use clebsch_gordan
  use Lagrange_weights
  use omp_lib
  implicit none

  type(complex_sb_mat)       :: h_grid, h_mid
  type(banded_sym_mat)       :: angular_v, v_mid, h_real, h_real2
  complex(8), allocatable    :: psi(:,:)
  real(8), allocatable       :: r_grid(:)
  real(8)                    :: t
  integer                    :: i, j, r_reach, time_step, info
  integer                    :: max_iter

  procedure(pulse_at_t_func), pointer  :: pulse

  call hydrogen_read
  call initial_print
  max_iter = 0

  pulse => select_pulse_type(pulse_name)

  allocate(r_grid(1:r_size), psi(1:r_size,1:l_max))
  call initialize_hydrogen(psi, r_grid, r_reach)
  call construct_operators(r_grid, r_reach, h_grid, angular_v, v_mid, h_mid, h_real, h_real2)

  if (soln_method == 'split_operator') then
     call short_time_half_step(h_mid, dt, r_grid(1:r_reach), psi(1:r_reach,:))
  end if

  call begin_timing
  t = 0
  max_iter = 0
  do time_step = 1,int(t_intv/dt)
     call adjust_dynamic_grid(psi, h_grid, h_real,  r_reach, h_mid, h_real2)
     
     if (soln_method == 'split_operator') then
        v_mid%offdiagonal =  pulse(t+0.5d0*dt)*angular_v%offdiagonal
        psi(1:r_reach,:) = split_operator_cn(time_step, h_mid, v_mid, r_grid(1:r_reach), dt, psi(1:r_reach,:), 1)
  
     else if (soln_method == 'itvolt') then
        if (itvolt_version == 'radial_lanczos') then
           call hydrogen_itvolt1(t, h_real2, angular_v, r_grid(1:r_reach), psi(1:r_reach,:), max_iter)
           
        else if (itvolt_version == 'radial_arnoldi') then
           call hydrogen_itvolt2(t, h_mid, angular_v, r_grid(1:r_reach), psi(1:r_reach,:), max_iter)
           
        else if (itvolt_version == 'split_operator_cn') then
           call hydrogen_itvolt3(t, h_mid, angular_v, r_grid(1:r_reach), psi(1:r_reach,:), max_iter)
           
        else if (itvolt_version == 'full_lanczos') then
           call hydrogen_itvolt4(t, h_real2, angular_v, r_grid(1:r_reach), psi(1:r_reach,:), max_iter)

        else
           print *, 'itvolt version not programmed'
           stop
           
        end if
        
     else
        print *, 'solution method not programmed'
        stop
        
     end if

     t = t + dt
     if (mod(time_step,10) == 0) then
        print *, 'Time: ', t, 'Number of radial points: ', r_reach
     end if

  end do 

  print *, 'Final time: ', t
  print *, 'Maximum number of iterations used by ITVOLT: ', max_iter

  if (soln_method == 'split_operator') then
     call short_time_half_step(h_mid, dt, r_grid(1:r_reach), psi(1:r_reach,:))
  end if

  call stop_timing

  call compute_ionization_probabilities(psi(1:r_reach,:), r_grid(1:r_reach), 'simpson')

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
    if (soln_method == 'split_operator') then
       print *, 'Solution method: Short Time Split Operator'
    else if (soln_method == 'itvolt') then
       print *, 'Solution method: Iterative Volterra Propagator (ITVOLT)'
       print *, 'Iteration type: ', it_type
       print *, 'Quadrature type: ', quad_type
       print *, 'Number of quadrature points:', quad_pt
       ! print *, 'Method for computing exponentials: ', prop_method
    else
       print *, 'Solution method: ', soln_method
    end if
    print *, 'Final propagation time:', t_intv
    print *, 'Propagation step size:', dt
    print *, '************************************'
    
  end subroutine initial_print


  subroutine initialize_hydrogen(psi, r_grid, r_reach)
    complex(8)                 :: psi(:,:)
    real(8)                    :: r_grid(:)
    integer                    :: r_reach
    integer                    :: i
    real(8)                    :: junk, wave_value

    ! construct the radial grid (equally spaced on [dr, r_max]
    do i = 1,r_size
       r_grid(i) = i*dr
    end do

    ! construct initial wavefunction (starts in 1s)
    do i = 1,r_size
       psi(i,1) = 2d0*r_grid(i)*exp(-r_grid(i))
    end do
    ! open(unit=65, file='h-0.01.wfn', status='old')
    ! do i = 1,1999
    !    read(65,*) junk, junk, wave_value
    !    psi(i,1) = dcmplx(wave_value,0d0)
    ! end do

    ! set initial dynamic grid and truncate
    do i = r_size,1,-1
       if (abs(psi(i,1)) > 1.d-9) then
          r_reach = i
          exit
       end if
    end do
    r_reach = min(r_reach+100, r_size)
    psi(r_reach+1:r_size,1) = 0d0

    print *, 'Number of steps: ', int(t_intv/dt)
    print *, 'Time: ', 0d0, 'Number of radial points: ', r_reach
    
  end subroutine initialize_hydrogen


  subroutine construct_operators(r_grid, r_reach, h_grid, angular_v, v_mid, h_mid, h_real, h_real2)
    real(8), intent(in)        :: r_grid(:)
    integer, intent(in)        :: r_reach
    type(complex_sb_mat)       :: h_grid, h_mid
    type(banded_sym_mat)       :: angular_v, v_mid, h_real, h_real2
    real(8), allocatable       :: diag(:), offd(:,:)
    integer                    :: i

    ! construct angular operators
    call construct_angularv(angular_v)
    v_mid = angular_v

    ! construct complex radial operators
    call h_grid%make_banded_matrix_on_grid(r_grid, dr, band_num_sym_mat)
    do i = 1,r_size
       h_grid%diagonal(i) = h_grid%diagonal(i) - 1d0/r_grid(i)
    end do
    call h_mid%initialize(r_reach, h_grid%bsz, h_grid%diagonal(1:r_reach), h_grid%offdiagonal(1:r_reach-1,:))

    ! construct real radial operators
    allocate(diag(1:r_size), offd(1:r_size-1,1:1))
    diag(:) = 1d0/dr**2
    offd(:,1) = -0.5d0/dr**2
    do i = 1,r_size
       diag(i) = diag(i) - 1d0/r_grid(i)
    end do
    call h_real%initialize(r_size, 1, diag, offd)
    call h_real2%initialize(r_reach, 1, h_real%diagonal(1:r_reach), h_real%offdiagonal(1:r_reach-1,:))
    
  end subroutine construct_operators

  
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


  subroutine adjust_dynamic_grid(psi, h_grid, h_real, r_reach, h_mid, h_real2)
    complex(8), intent(in)            :: psi(:,:)
    type(complex_sb_mat), intent(in)  :: h_grid
    type(banded_sym_mat), intent(in)  :: h_real
    integer                           :: r_reach
    type(complex_sb_mat)              :: h_mid
    type(banded_sym_mat)              :: h_real2
    integer                           :: i

    ! adjust grid size
    do i = 1,l_max
       if (abs(psi(r_reach-100,i)) > grid_tol) then
          r_reach = r_reach + 100
          exit
       end if
    end do
    r_reach = min(r_reach, r_size)

    ! update operators
    h_mid%mat_size = r_reach
    h_mid%diagonal = h_grid%diagonal(1:r_reach)
    h_mid%offdiagonal = h_grid%offdiagonal(1:r_reach-1,:)
    h_real2%mat_size = r_reach
    h_real2%diagonal = h_real%diagonal(1:r_reach)
    h_real2%offdiagonal = h_real%offdiagonal(1:r_reach-1,:)
    
  end subroutine adjust_dynamic_grid


  subroutine compute_ionization_probabilities(psi, r_grid, quadrature)
    complex(8), intent(in)        :: psi(:,:)
    real(8), intent(in)           :: r_grid(:)
    character(len=*), intent(in)  :: quadrature
    real(8), allocatable          :: momenta(:), spectra(:)
    real(8), allocatable          :: coulomb_waves(:,:), GC(:), FC(:), FCP(:), GCP(:)
    complex(8), allocatable       :: z_el(:)
    real(8)                       :: total_prob
    integer                       :: i, j, k, p, d, ifail

    print *, '************************************'
    print *, 'Computing ionization probabilities!'
    print *, '************************************'

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

       ! integrate produce of psi and Coulomb functions and sum
       z_el(:) = 0
       do k = 1,l_max
          ! Simpson's rule
          if (quadrature == 'simpson') then
             if (mod(d,2) == 0) then
                do p = 1,d/2
                   z_el(k) = z_el(k) + 4d0*coulomb_waves(2*p-1,k)*psi(2*p-1,k)
                end do

                do p = 1,d/2-1
                   z_el(k) = z_el(k) + 2d0*coulomb_waves(2*p,k)*psi(2*p,k)
                end do

                z_el(k) = z_el(k) + coulomb_waves(d,k)*psi(d,k)
                z_el(k) = (1d0/3d0)*dr*sqrt(2d0/(pi*momenta(i)))*z_el(k)
                spectra(i) = spectra(i) + abs(z_el(k))**2
             else
                do p = 1,(d-1)/2
                   z_el(k) = z_el(k) + 4d0*coulomb_waves(2*p-1,k)*psi(2*p-1,k)
                end do

                do p = 1,(d-1)/2-1
                   z_el(k) = z_el(k) + 2d0*coulomb_waves(2*p,k)*psi(2*p,k)
                end do

                z_el(k) = z_el(k) + coulomb_waves(d-1,k)*psi(d-1,k)
                z_el(k) = (1d0/3d0)*z_el(k)
                z_el(k) = z_el(k) + 0.5d0*coulomb_waves(d-1,k)*psi(d-1,k)
                z_el(k) = z_el(k) + 0.5d0*coulomb_waves(d,k)*psi(d,k)
                z_el(k) = dr*sqrt(2d0/(pi*momenta(i)))*z_el(k)
                spectra(i) = spectra(i) + abs(z_el(k))**2
             
             end if

          ! composite trapezoidal rule   
          else if (quadrature == 'trapezoid') then   
             do p = 1,d-1
                z_el(k) = z_el(k) + sqrt(2d0/(pi*momenta(i)))*dr*coulomb_waves(p,k)*psi(p,k)
             end do
             z_el(k) = z_el(k) + 0.5d0*sqrt(2d0/(pi*momenta(i)))*dr*coulomb_waves(d,k)*psi(d,k)
             spectra(i) = spectra(i) + abs(z_el(k))**2

          else
             print *, 'quadrature for ionization probabilities not programmed'
             stop
             
          end if
          
       end do

       write(87,*) e_min+(i-1)*dE, spectra(i)
    end do

    close(87)

    ! integrate total ionization probability
    ! Simpson's rule
    if (quadrature == 'simpson') then
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

    ! composite trapezoidal rule   
    else if (quadrature == 'trapezoid') then  
       total_prob = 0.5d0*0.002d0*spectra(1)
       do i = 2,coul_num-1
          total_prob = total_prob + 0.002d0*spectra(i)
       end do
       total_prob = total_prob + 0.5d0*0.002d0*spectra(coul_num)
       
    end if
    
    print *, 'Total Ionization Probability: ', total_prob

    ! write final wavefunction in data file
    open(unit=88, file='final_wave_function.dat')
    do i = 1,d
       write(88,*) r_grid(i), real(psi(i,:)), aimag(psi(i,:))
    end do
    close(88)
    
  end subroutine compute_ionization_probabilities


  ! e^{-ih_mid(dt/2)}*psi done by crank-nicolson
  subroutine short_time_half_step(h_mid, local_dt, r_grid, psi)
    type(complex_sb_mat), intent(in)    :: h_mid
    real(8), intent(in)                 :: local_dt
    real(8), intent(in)                 :: r_grid(:)
    complex(8)                          :: psi(:,:)
    complex(8), allocatable             :: h_diagonal(:), r_rhs(:)
    complex(8), allocatable             :: r_offdiag(:), r_offdiag2(:), r_diagonal(:)
    integer                             :: d, i, j, info

    d = size(r_grid(:))
    allocate(h_diagonal(1:d), r_rhs(1:d), r_offdiag(1:d-1), r_offdiag2(1:d-1), r_diagonal(1:d))

    !$OMP parallel do private(i,j,h_diagonal,r_rhs,r_diagonal,r_offdiag,r_offdiag2,info)
     do i = 1,l_max
        ! update h_grid with l
        do j = 1,d
           h_diagonal(j) = h_mid%diagonal(j) + dble((i-1)*i)/(2d0*r_grid(j)**2)
        end do

        ! construct RHS vector and coefficient diagonal/off diagonal
        r_rhs(:) = psi(1:d,i) - ii*(local_dt/4d0)*h_diagonal(:)*psi(1:d,i)
        r_rhs(1:d-1) = r_rhs(1:d-1) - ii*(local_dt/4d0)*h_mid%offdiagonal(:,1)*psi(2:d,i)
        r_rhs(2:d) = r_rhs(2:d) - ii*(local_dt/4d0)*h_mid%offdiagonal(:,1)*psi(1:d-1,i)
        do j = 1,d
           r_diagonal(j) = 1d0 + ii*(local_dt/4d0)*h_diagonal(j)
        end do
        r_offdiag(:) = ii*(local_dt/4d0)*h_mid%offdiagonal(:,1)
        r_offdiag2(:) = r_offdiag(:)

        ! call tridiagonal system solve
        call zgtsv(d, 1, r_offdiag, r_diagonal, r_offdiag2, r_rhs, r_reach, info)
        if (info /= 0) then
           print *, 'system solve in split operator failed'
           stop
        end if

        psi(1:d,i) = r_rhs(:)
     end do
     !$OMP end parallel do
    
  end subroutine short_time_half_step
   

  ! real radial exponentials done via lanczos
  subroutine hydrogen_itvolt1(time, h_real, v, r_grid, psi, max_iters)
    real(8), intent(in)                 :: time
    type(banded_sym_mat), intent(in)    :: h_real, v
    real(8), intent(in)                 :: r_grid(:)
    complex(8)                          :: psi(:,:)
    integer                             :: max_iters
    type(banded_sym_mat)                :: mat
    real(8), allocatable                :: pt(:), wt(:,:), comp_wt(:,:), pot(:)
    complex(8), allocatable             :: inhomogeneity(:,:,:), iterative_ans(:,:,:)
    complex(8), allocatable             :: phi(:,:,:), v_psi(:,:,:), b(:,:)
    logical                             :: converged
    integer                             :: d, it_num, n, i, j, k, l, p

    procedure(pulse_at_t_func), pointer  :: pulse
    procedure(potential_func), pointer   :: potential

    ! set up parameters
    pulse => select_pulse_type(pulse_name)
    potential => select_potential_type(potential_type)
    d = size(r_grid)
    it_num = 0
    n = quad_pt

    allocate(pt(1:n), wt(1:n,1:n-1), comp_wt(1:n,1:n-1), inhomogeneity(1:d,1:l_max,1:n), iterative_ans(1:d,1:l_max,1:n), &
         phi(1:d,1:l_max,1:n), v_psi(1:d,1:l_max,1:n), b(1:d,1:l_max), pot(1:d))

    ! compute quadrature points and weights
    call lgr_weights(time, time+dt, pt, wt, n-2, quad_type)

    comp_wt(:,1) = wt(:,1)
    do i = 2,n-1
       comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
    end do

    call mat%initialize(h_real%mat_size, h_real%bsz, h_real%diagonal, h_real%offdiagonal)
    pot(:) = potential(r_grid)
    
    ! evaluate inhomogeneity at the quadrature points
    inhomogeneity(:,:,1) = psi(:,:)
    do j = 2,n
       do i = 1,l_max
          do k = 1,d
             mat%diagonal(k) = h_real%diagonal(k) + dble(i*(i-1))/(2d0*r_grid(k)**2)
          end do
          inhomogeneity(:,i,j) = lanczos_prop(mat, pt(j)-pt(1), psi(:,i))
       end do
    end do

    converged = .FALSE.
    iterative_ans(:,:,:) = inhomogeneity(:,:,:)

    ! iterate until converged
    do while (.not. (converged .or. it_num >= it_cap))
      phi(:,:,:) = iterative_ans(:,:,:)

      if (it_type == 'jacobi') then
         do i = 1,n
            do j = 1,d
               v_psi(j,:,i) = r_grid(j)*pulse(pt(i))*sb_matvec_mul(v,iterative_ans(j,:,i)) - ii*pot(j)*iterative_ans(j,:,i)
            end do
         end do

         do k = 2,n
            b(:,:) = inhomogeneity(:,:,k)

            do j = 1,k-1
               do l = 1,l_max
                  do p = 1,d
                     mat%diagonal(p) = h_real%diagonal(p) + dble(l*(l-1))/(2d0*r_grid(p)**2)
                  end do
                  b(:,l) = b(:,l) - ii*comp_wt(j,k-1)*lanczos_prop(mat, pt(k)-pt(j), v_psi(:,l,j))
               end do
            end do

            b(:,:) = b(:,:) - ii*comp_wt(k,k-1)*v_psi(:,:,k)

            do j = k+1,n
               do l = 1,l_max
                  do p = 1,d
                     mat%diagonal(p) = h_real%diagonal(p) + dble(l*(l-1))/(2d0*r_grid(p)**2)
                  end do
                   b(:,l) = b(:,l) - ii*comp_wt(j,k-1)*lanczos_prop(mat, pt(k)-pt(j), v_psi(:,l,j))
               end do
            end do

            iterative_ans(:,:,k) = b(:,:)
         end do

      else
         print *, 'ITVOLT iteration type not programmed'
         stop
         
      end if

      it_num = it_num + 1

      ! test for convergence
      converged = convergence_check(phi, iterative_ans, it_tolerance)

    end do

    psi(:,:) = iterative_ans(:,:,n)

    if (max_iters < it_num) then
       max_iters = it_num
    end if

  end subroutine hydrogen_itvolt1

  ! complex radial exponentials done via arnoldi
  subroutine hydrogen_itvolt2(time, h_mid, v, r_grid, psi, max_iters)
    real(8), intent(in)                 :: time
    type(complex_sb_mat), intent(in)    :: h_mid
    type(banded_sym_mat), intent(in)    :: v
    real(8), intent(in)                 :: r_grid(:)
    type(complex_sb_mat)                :: mat
    complex(8)                          :: psi(:,:)
    integer                             :: max_iters
    real(8)                             :: E_diff
    type(banded_sym_mat)                :: v_mid
    real(8), allocatable                :: pt(:), wt(:,:), comp_wt(:,:)
    complex(8), allocatable             :: inhomogeneity(:,:,:), iterative_ans(:,:,:), phi(:,:,:)
    complex(8), allocatable             :: v_psi(:,:,:), b(:,:)
    logical                             :: converged
    integer                             :: d, it_num, n, i, j, k, l, p

    procedure(pulse_at_t_func), pointer  :: pulse

    ! set up parameters
    pulse => select_pulse_type(pulse_name)
    d = size(r_grid)
    it_num = 0
    n = quad_pt

    allocate(pt(1:n), wt(1:n,1:n-1), comp_wt(1:n,1:n-1), inhomogeneity(1:d,1:l_max,1:n), iterative_ans(1:d,1:l_max,1:n), &
         phi(1:d,1:l_max,1:n), v_psi(1:d,1:l_max,1:n), b(1:d,1:l_max))

    ! compute quadrature points and weights
    call lgr_weights(time, time+dt, pt, wt, n-2, quad_type)

    comp_wt(:,1) = wt(:,1)
    do i = 2,n-1
       comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
    end do

    call mat%initialize(h_mid%mat_size, h_mid%bsz, h_mid%diagonal, h_mid%offdiagonal)
       
    ! evaluate inhomogeneity at the quadrature points
    inhomogeneity(:,:,1) = psi(:,:)
    do j = 2,n
       do k = 1,l_max
          do p = 1,d
             mat%diagonal(p) = h_mid%diagonal(p) + dble(k*(k-1))/(2d0*r_grid(p)**2)
          end do
          inhomogeneity(:,k,j) = arnoldi_prop(mat, pt(j)-pt(1), psi(:,k))
       end do
    end do

    converged = .FALSE.
    iterative_ans(:,:,:) = inhomogeneity(:,:,:)

    ! iterate until converged
    do while (.not. (converged .or. it_num >= it_cap))
      phi(:,:,:) = iterative_ans(:,:,:)

      if (it_type == 'jacobi') then
         do i = 1,n
            do j = 1,d
               v_psi(j,:,i) = r_grid(j)*pulse(pt(i))*sb_matvec_mul(v,iterative_ans(j,:,i))
             end do
         end do

         do k = 2,n
            b(:,:) = inhomogeneity(:,:,k)

            do j = 1,k-1
               do l = 1,l_max
                  do p = 1,d
                     mat%diagonal(p) = h_mid%diagonal(p) + dble(l*(l-1))/(2d0*r_grid(p)**2)
                  end do
                  b(:,l) = b(:,l) - ii*comp_wt(j,k-1)*arnoldi_prop(mat, pt(k)-pt(j), v_psi(:,l,j))
               end do
            end do

            b(:,:) = b(:,:) - ii*comp_wt(k,k-1)*v_psi(:,:,k)

            do j = k+1,n
               do l = 1,l_max
                  do p = 1,d
                     mat%diagonal(p) = h_mid%diagonal(p) + dble(l*(l-1))/(2d0*r_grid(p)**2)
                  end do
                   b(:,l) = b(:,l) - ii*comp_wt(j,k-1)*arnoldi_prop(mat, pt(k)-pt(j), v_psi(:,l,j))
               end do
            end do

            iterative_ans(:,:,k) = b(:,:)
         end do

       else
         print *, 'ITVOLT iteration type not programmed'
         stop
       end if

       it_num = it_num + 1

       ! test for convergence
       converged = convergence_check(phi, iterative_ans, it_tolerance)
          
    end do

    psi(:,:) = iterative_ans(:,:,n)

    if (max_iters < it_num) then
       max_iters = it_num
    end if
    
  end subroutine hydrogen_itvolt2


  ! complex midpoint exponentials done by split operator + crank-nicolson
  subroutine hydrogen_itvolt3(time, h_mid, v, r_grid, psi, max_iters)
    real(8), intent(in)                 :: time
    type(complex_sb_mat), intent(in)    :: h_mid
    type(banded_sym_mat), intent(in)    :: v
    real(8), intent(in)                 :: r_grid(:)
    type(complex_sb_mat)                :: mat
    complex(8)                          :: psi(:,:)
    integer                             :: max_iters
    real(8)                             :: E_diff
    type(banded_sym_mat)                :: v_mid
    real(8), allocatable                :: pt(:), wt(:,:), comp_wt(:,:)
    complex(8), allocatable             :: inhomogeneity(:,:,:), iterative_ans(:,:,:), phi(:,:,:)
    complex(8), allocatable             :: v_psi(:,:,:), b(:,:)
    logical                             :: converged
    integer                             :: d, it_num, n, i, j, k, l, p

    procedure(pulse_at_t_func), pointer  :: pulse

    ! set up parameters
    pulse => select_pulse_type(pulse_name)
    d = size(r_grid)
    it_num = 0
    n = quad_pt

    allocate(pt(1:n), wt(1:n,1:n-1), comp_wt(1:n,1:n-1), inhomogeneity(1:d,1:l_max,1:n), iterative_ans(1:d,1:l_max,1:n), &
         phi(1:d,1:l_max,1:n), v_psi(1:d,1:l_max,1:n), b(1:d,1:l_max))

    ! compute quadrature points and weights
    call lgr_weights(time, time+dt, pt, wt, n-2, quad_type)

    comp_wt(:,1) = wt(:,1)
    do i = 2,n-1
       comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
    end do

    ! construct midpoint v
    call v_mid%initialize(v%mat_size, v%bsz, v%diagonal, pulse(time+0.5d0*dt)*v%offdiagonal)

    ! evaluate inhomogeneity at the quadrature points
    inhomogeneity(:,:,1) = psi(:,:)
    do j = 2,n
       inhomogeneity(:,:,j) = split_operator_cn(1, h_mid, v_mid, r_grid, pt(j)-pt(1), psi, 2)
    end do

    converged = .FALSE.
    iterative_ans(:,:,:) = inhomogeneity(:,:,:)

    ! iterate until converged
    do while (.not. (converged .or. it_num >= it_cap))
       phi(:,:,:) = iterative_ans(:,:,:)

       if (it_type == 'jacobi') then
          do i = 1,n
             E_diff = pulse(pt(i)) - pulse(time+0.5d0*dt)
             do j = 1,d
                v_psi(j,:,i) = r_grid(j)*E_diff*sb_matvec_mul(v,iterative_ans(j,:,i))
             end do
          end do

          do k = 2,n
             b(:,:) = inhomogeneity(:,:,k)

             do j = 1,k-1
                b(:,:) = b(:,:) - ii*comp_wt(j,k-1)*split_operator_cn(1, h_mid, v_mid, r_grid, pt(k)-pt(j), v_psi(:,:,j), 2)
             end do

             b(:,:) = b(:,:) - ii*comp_wt(k,k-1)*v_psi(:,:,k)

             do j = k+1,n
                b(:,:) = b(:,:) - ii*comp_wt(j,k-1)*split_operator_cn(1, h_mid, v_mid, r_grid, pt(k)-pt(j), v_psi(:,:,j), 2)
             end do

             iterative_ans(:,:,k) = b(:,:)
          end do
          
       else
          print *, 'ITVOLT iteration type not programmed'
          stop
          
       end if
          
       it_num = it_num + 1

       ! test for convergence
       converged = convergence_check(phi, iterative_ans, it_tolerance)
       
    end do

    psi(:,:) = iterative_ans(:,:,n)

    if (max_iters < it_num) then
       max_iters = it_num
    end if
    
  end subroutine hydrogen_itvolt3


  ! real midpoint exponentials done by arnodi (re-using vectors)
  subroutine hydrogen_itvolt4(time, h_real, v, r_grid, psi, max_iters)
    real(8), intent(in)                 :: time
    type(banded_sym_mat), intent(in)    :: h_real
    type(banded_sym_mat), intent(in)    :: v
    real(8), intent(in)                 :: r_grid(:)
    type(complex_sb_mat)                :: mat
    complex(8)                          :: psi(:,:)
    integer                             :: max_iters
    real(8)                             :: E_diff
    type(banded_sym_mat)                :: v_mid
    real(8), allocatable                :: pt(:), wt(:,:), comp_wt(:,:)
    complex(8), allocatable             :: inhomogeneity(:,:,:), iterative_ans(:,:,:), phi(:,:,:)
    complex(8), allocatable             :: v_psi(:,:,:), b(:,:), Q(:,:)
    real(8), allocatable                :: eigenvalues(:,:), eigenvectors(:,:,:), pot(:)
    logical                             :: converged
    integer                             :: d, it_num, n, i, j, k, l, p, lancz_num

    procedure(pulse_at_t_func), pointer  :: pulse
    procedure(potential_func), pointer   :: potential

    ! set up parameters
    pulse => select_pulse_type(pulse_name)
    potential => select_potential_type(potential_type)
    d = size(r_grid)
    it_num = 0
    n = quad_pt

    allocate(pt(1:n), wt(1:n,1:n-1), comp_wt(1:n,1:n-1), inhomogeneity(1:d,1:l_max,1:n), iterative_ans(1:d,1:l_max,1:n), &
         phi(1:d,1:l_max,1:n), v_psi(1:d,1:l_max,1:n), b(1:d,1:l_max))

    ! compute quadrature points and weights
    call lgr_weights(time, time+dt, pt, wt, n-2, quad_type)

    comp_wt(:,1) = wt(:,1)
    do i = 2,n-1
       comp_wt(:,i) = comp_wt(:,i-1) + wt(:,i)
    end do
    
    ! construct midpoint v
    call v_mid%initialize(v%mat_size, v%bsz, v%diagonal, pulse(time+0.5d0*dt)*v%offdiagonal)

    ! get Lanczos info
    lancz_num = 100 ! make sure this is a multiple of 10
    k = lancz_num/10
    allocate(Q(1:d*l_max,1:lancz_num), eigenvalues(1:lancz_num,1:k), eigenvectors(1:lancz_num,1:lancz_num,1:k), pot(1:d))
    pot(:) = potential(r_grid)
    call get_lanczos_vectors(h_real, v_mid, r_grid, psi, Q, eigenvectors, eigenvalues)

    ! evaluate inhomogeneity at the quadrature points
    inhomogeneity(:,:,1) = psi(:,:)
    do j = 2,n
       inhomogeneity(:,:,j) = hydrogen_lanczos(Q, eigenvalues, eigenvectors, pt(j)-pt(1), psi)
    end do

    converged = .FALSE.
    iterative_ans(:,:,:) = inhomogeneity(:,:,:)

    ! iterate until converged
    do while (.not. (converged .or. it_num >= it_cap))
       phi(:,:,:) = iterative_ans(:,:,:)

       if (it_type == 'jacobi') then
          do i = 1,n
             E_diff = pulse(pt(i)) - pulse(time+0.5d0*dt)
             do j = 1,d
                v_psi(j,:,i) = r_grid(j)*E_diff*sb_matvec_mul(v,iterative_ans(j,:,i)) - ii*pot(j)*iterative_ans(j,:,i)
             end do
          end do

          do k = 2,n
             b(:,:) = inhomogeneity(:,:,k)

             do j = 1,k-1
                b(:,:) = b(:,:) - ii*comp_wt(j,k-1)*hydrogen_lanczos(Q, eigenvalues, eigenvectors, pt(k)-pt(j), v_psi(:,:,j))
             end do

             b(:,:) = b(:,:) - ii*comp_wt(k,k-1)*v_psi(:,:,k)

             do j = k+1,n
                b(:,:) = b(:,:) - ii*comp_wt(j,k-1)*hydrogen_lanczos(Q, eigenvalues, eigenvectors, pt(k)-pt(j), v_psi(:,:,j))
             end do

             iterative_ans(:,:,k) = b(:,:)
          end do
          
       else
          print *, 'ITVOLT iteration type not programmed'
          stop
          
       end if
          
       it_num = it_num + 1

       ! test for convergence
       converged = convergence_check(phi, iterative_ans, it_tolerance)
       
    end do

    psi(:,:) = iterative_ans(:,:,n)

    if (max_iters < it_num) then
       max_iters = it_num
    end if
    
  end subroutine hydrogen_itvolt4


  function convergence_check(wave1, wave2, tolerance)  result(ans)
    complex(8), intent(in)        :: wave1(:,:,:), wave2(:,:,:)
    real(8), intent(in)           :: tolerance
    complex(8), allocatable       :: v1(:), v2(:)
    real(8)                       :: diff, max_diff
    logical                       :: ans
    integer                       :: d, n, i, j

    d = size(wave1(:,1,1))
    n = size(wave1(1,1,:))
    max_diff = 0d0

    allocate(v1(1:d*l_max), v2(1:d*l_max))

    do i = 1,n
       do j = 1,l_max
          v1((j-1)*d+1:j*d) = wave1(:,j,i)
          v2((j-1)*d+1:j*d) = wave2(:,j,i)
       end do

       diff = sqrt(dot_product(v1-v2, v1-v2))
       if (max_diff < diff) then
          max_diff = diff
       end if
    end do

    if (max_diff <= tolerance) then
       ans = .TRUE.
    else
       ans = .FALSE.
    end if    

  end function convergence_check


  subroutine get_lanczos_vectors(h_mid, v_mid, r_grid, psi, Q, eigenvectors, eigenvalues)
    type(banded_sym_mat), intent(in)              :: h_mid, v_mid
    real(8), intent(in)                           :: r_grid(:)
    complex(8), intent(in)                        :: psi(:,:)
    complex(8)                                    :: Q(:,:)
    real(8), allocatable                          :: alpha(:), beta(:)
    real(8)                                       :: eigenvalues(:,:), eigenvectors(:,:,:)
    real(8), allocatable                          :: off_diagonal(:), work(:)
    complex(8), allocatable                       :: w(:)
    integer                                       :: d, i, m, p, k, info

    d = h_mid%mat_size
    k = size(Q(1,:))
    p = size(eigenvalues(1,:))
    allocate(w(1:d*l_max), off_diagonal(1:k-1), work(1:2*k-2), alpha(1:k), beta(1:k-1))

    ! first vector is normalized psi
    do i = 1,l_max
       Q((i-1)*d+1:i*d,1) = psi(:,i)
    end do
    Q(:,1) = Q(:,1)/sqrt(dot_product(Q(:,1),Q(:,1)))

    ! compute second vector
    w(:) = hydrogen_matvec_mul(h_mid, v_mid, r_grid, Q(:,1))
    alpha(1) = dot_product(Q(:,1),w)
    w(:) = w(:) - alpha(1)*Q(:,1)
    beta(1) = sqrt(dot_product(w,w))
    Q(:,2) = w(:)/beta(1)

    ! re-orthogonalize against first
    call zschmab(Q(:,1), Q(:,2), 1.d-10, d*l_max, 1, 1, m)

    ! compute remaining lanczos vectors
    do i = 2,k
       w(:) = hydrogen_matvec_mul(h_mid, v_mid, r_grid, Q(:,i))
       alpha(i) = dot_product(Q(:,i),w)
       if (i /= k) then
          w(:) = w(:) - beta(i-1)*Q(:,i-1) - alpha(i)*Q(:,i)
          beta(i) = sqrt(dot_product(w,w))
          Q(:,i+1) = w(:)/beta(i)

          call zschmab(Q(:,1:i), Q(:,i+1), 1.d-10, d*l_max, i, 1, m)
       end if
    end do

    ! get eigenvalues/eigenvectors
    do i = 1,p
       eigenvalues(1:10*i,i) = alpha(1:10*i)
       off_diagonal(1:10*i-1) = beta(1:10*i-1)
       call dstev('V', 10*i, eigenvalues(1:10*i,i), off_diagonal(1:10*i-1), eigenvectors(1:10*i,1:10*i,i), &
            10*i, work(1:20*i-2), info)

       if (info /= 0) then
          print *, 'eigenvsolve in lanczos failed'
          stop
       end if
    end do
    
  end subroutine get_lanczos_vectors


  function hydrogen_lanczos(Q, eigenvalues, eigenvectors, local_dt, psi) result(ans)
    complex(8), intent(in)                     :: Q(:,:)
    real(8), intent(in)                        :: eigenvalues(:,:), eigenvectors(:,:,:)
    real(8), intent(in)                        :: local_dt
    complex(8), intent(in)                     :: psi(:,:)
    complex(8), allocatable                    :: ans(:,:)
    complex(8), allocatable                    :: psi_vec(:), dummy(:), ans_vec(:), phi(:)
    real(8)                                    :: error
    logical                                    :: converged
    integer                                    :: d, k, i, j

    
    ! set parameters
    d = size(psi(:,1))
    k = size(eigenvalues(1,:))
    allocate(ans(1:d,1:l_max), psi_vec(1:d*l_max), dummy(1:10*k), phi(1:d*l_max), ans_vec(1:d*l_max))
    do i = 1,l_max
       psi_vec((i-1)*d+1:i*d) = psi(:,i)
    end do

    ! compute first approximation
    dummy(1:10) = matmul(transpose(eigenvectors(1:10,1:10,1)), matmul(transpose(conjg(Q(:,1:10))), psi_vec))
    do i = 1,10
       dummy(i) = cmplx(cos(-local_dt*eigenvalues(i,1)), sin(-local_dt*eigenvalues(i,1)), 8)*dummy(i)
    end do
    ans_vec(:) = matmul(Q(:,1:10), matmul(eigenvectors(1:10,1:10,1), dummy(1:10)))

    ! iterate (adding 10 more vectors each time)
    converged = .FALSE.
    phi(:) = ans_vec(:)
    i = 2
    do while (.not. (converged .or. i > k))
       dummy(1:10*i) = matmul(transpose(eigenvectors(1:10*i,1:10*i,i)), matmul(transpose(conjg(Q(:,1:10*i))), psi_vec))
       do j = 1,10*i
          dummy(j) = cmplx(cos(-local_dt*eigenvalues(j,i)), sin(-local_dt*eigenvalues(j,i)), 8)*dummy(j)
       end do
       ans_vec(:) = matmul(Q(:,1:10*i), matmul(eigenvectors(1:10*i,1:10*i,i), dummy(1:10*i)))

       error = sqrt(dot_product(ans_vec-phi, ans_vec-phi))
       
       if (error < 1.d-3) then
          converged = .TRUE.
       else
          phi(:) = ans_vec(:)
          i = i+1
       end if
       
    end do

    if (.not. converged) then
       print *, 'error: more lanczos vectors needed'
       print *, error
       stop
    end if

    ! reorganize answer as an array
    do i = 1,l_max
       ans(:,i) = ans_vec((i-1)*d+1:i*d)
    end do
    
  end function hydrogen_lanczos


  function hydrogen_matvec_mul(h_mid, v_mid, r_grid, vec) result(ans)
    type(banded_sym_mat), intent(in)           :: h_mid, v_mid
    complex(8), intent(in)                     :: vec(:)
    real(8), intent(in)                        :: r_grid(:)
    complex(8), allocatable                    :: ans(:)
    type(banded_sym_mat)                       :: h_copy
    complex(8), allocatable                    :: vec_array(:,:), ans_array(:,:)
    integer                                    :: d, i, j

    d = h_mid%mat_size
    allocate(ans(1:d*l_max), vec_array(1:d,1:l_max), ans_array(1:d,1:l_max))
    call h_copy%initialize(d, h_mid%bsz, h_mid%diagonal, h_mid%offdiagonal)

    ! reorganize vec as an array
    do i = 1,l_max
       vec_array(:,i) = vec((i-1)*d+1:i*d)
    end do

    ! multiply by h_mid + v_mid
    do i = 1,d
       ans_array(i,:) = r_grid(i)*sb_matvec_mul(v_mid, vec_array(i,:))
    end do

    do i = 1,l_max
       do j = 1,d
          h_copy%diagonal(j) = h_mid%diagonal(j) + dble(i*(i-1))/(2d0*r_grid(j)**2)
       end do

       ans_array(:,i) = ans_array(:,i) + sb_matvec_mul(h_copy, vec_array(:,i))
    end do

    ! reorganize answer itno a vector
    do i = 1,l_max
       ans((i-1)*d+1:i*d) = ans_array(:,i)
    end do
    
  end function hydrogen_matvec_mul
  
  
end program hydrogen
