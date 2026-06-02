!# 1 "xc_rVV10.f90"
!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE rVV10 
  !--------------------------------------------------------------------------
  !! This module is modeled after the vdW-DF implementation in
  !! 'Modules/xc_vdW_DF.f90'. See that file for references, explanations, and
  !! many useful comments.
!# 15 "xc_rVV10.f90"
  USE kinds,             ONLY : dp
  USE constants,         ONLY : pi
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE io_global,         ONLY : ionode, stdout
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft 
  USE control_flags,     ONLY : gamma_only, iverbosity
!# 24 "xc_rVV10.f90"
  IMPLICIT NONE
  SAVE
!# 27 "xc_rVV10.f90"
  real(dp), parameter :: epsr      = 1.d-12
  real(dp), parameter :: epsg      = 1.D-10
  integer,  parameter :: Nr_points = 1024
  real(dp), parameter :: r_max     = 100.0D0
  real(dp), parameter :: dr        = r_max/Nr_points
  real(dp), parameter :: dk        = 2.0D0*pi/r_max
  real(dp), parameter :: q_min     = 1.0D-4
  real(dp), parameter :: q_cut     = 0.5D0
  integer,  parameter :: Nqs       = 20
  real(dp), parameter, dimension(Nqs):: q_mesh= (/ q_min, 3.0D-4, 5.893850845618885D-4, 1.008103720396345D-3, &
            1.613958359589310D-3, 2.490584839564653D-3, 3.758997979748929D-3, 5.594297198907115D-3, &
            8.249838297569416D-3, 1.209220822453922D-2, 1.765183095571029D-2, 2.569619042667097D-2, &
            3.733577865542191D-2, 5.417739477463518D-2, 7.854595729872216D-2, 0.113805449932145D0,  &
            0.164823306218807D0 , 0.238642339497217D0 , 0.345452975434964D0 , q_cut /)
!# 42 "xc_rVV10.f90"
  real(dp) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )
!# 44 "xc_rVV10.f90"
  real(dp) :: b_value = 6.3_DP
  real(dp) :: C_value = 0.0093 
!# 47 "xc_rVV10.f90"
  private  
  public :: xc_rVV10,  &
            interpolate_kernel, &
            initialize_spline_interpolation, &
            rVV10_stress, b_value, &
            q_mesh, Nr_points, r_max, q_min, q_cut, Nqs
!# 54 "xc_rVV10.f90"
CONTAINS
!# 56 "xc_rVV10.f90"
! #################################################################################################
!                                       |             |
!                                       |  xc_rVV10   |
!                                       |_____________|
!# 61 "xc_rVV10.f90"
  SUBROUTINE xc_rVV10(rho_valence, rho_core, nspin, etxc, vtxc, v, b_value_)
  
    !! Calculate exchange-correlation energy and potential for rVV10.
!# 65 "xc_rVV10.f90"
    ! Modules to include
    ! -------------------------------------------------------------------------
    
    use gvect,           ONLY : ngm, g
    USE fft_base,        ONLY : dfftp
    USE cell_base,       ONLY : omega, tpiba
    ! -------------------------------------------------------------------------
    
    real(dp), intent(IN) :: rho_valence(:)
    !! valence charge density
    real(dp), intent(IN) :: rho_core(:)
    !! core charge density
    INTEGER,  INTENT(IN) :: nspin
    !! number of spin components
    real(dp), intent(inout) :: etxc
    !! total XC energy
    real(dp), intent(inout) :: vtxc
    !! total XC potential
    real(dp), intent(inout) :: v(:,:)
    !! XC potential on rho grid
    real(DP),optional,intent(in) :: b_value_
    
    !
    ! Local variables
    ! ----------------------------------------------------------------------------------
    !   
    
    integer :: i_grid, theta_i, i_proc, I      
    real(dp) :: grid_cell_volume                
    
    real(dp), allocatable :: total_rho(:)
    real(dp), allocatable :: gradient_rho(:,:)  
!# 98 "xc_rVV10.f90"
    real(dp), allocatable :: q0(:)                                        
    real(dp), allocatable :: dq0_drho(:) 
    real(dp), allocatable :: dq0_dgradrho(:) 
    complex(dp), allocatable :: thetas(:,:)    
                                                
    real(dp) :: Ec_nl                          
    real(dp), allocatable :: potential(:)                                                 
!# 106 "xc_rVV10.f90"
    logical, save :: first_iteration = .true.  
!# 108 "xc_rVV10.f90"
    real(dp) ::  beta
 
 
    ! ---------------------------------------------------------------------------------------------
    !   Begin calculations
  
    !call errore('xc_rVV10','rVV10 functional not implemented for spin polarized runs', size(rho_valence,2)-1)
    if (nspin>2) call errore('xc_vdW_DF','vdW functional not implemented for nspin > 2', nspin)
!# 117 "xc_rVV10.f90"
    if(present(b_value_)) b_value = b_value_
    
    ! --------------------------------------------------------------------------------------------------------
!# 121 "xc_rVV10.f90"
    call start_clock( 'rVV10' )
!# 123 "xc_rVV10.f90"
    beta = 0.0625d0 * (3.0d0 / (b_value**2.0D0) )**(0.75d0)
!# 126 "xc_rVV10.f90"
    ! Write parameters during the first iteratio
    !
    if (first_iteration) then
!# 130 "xc_rVV10.f90"
       first_iteration = .false.
      
       CALL generate_kernel
!# 134 "xc_rVV10.f90"
       if (ionode .and. iverbosity > -1 ) then
!# 136 "xc_rVV10.f90"
          WRITE(stdout,'(/ /A )') "---------------------------------------------------------------------------------"
          WRITE(stdout,'(A)') "Carrying out rVV10 run using the following parameters:"
          WRITE(stdout,'(A,I6,A,I6,A,F8.3)') "Nqs =  ",Nqs, "    Nr_points =  ", Nr_points,"   r_max =  ",r_max
          WRITE(stdout, '(A, F8.5, A, F8.5 )') "b_value = ", b_value, "    beta = ", beta        
          WRITE(stdout,'(5X,"q_mesh =",4F12.8)') (q_mesh(I), I=1, 4)
          WRITE(stdout,'(13X,4F12.8)') (q_mesh(I), I=5, Nqs)
                 
          WRITE(stdout,'(/ A )') "Gradients computed in Reciprocal space"
          WRITE(stdout,'(/ A / /)') "---------------------------------------------------------------------------------"
!# 146 "xc_rVV10.f90"
          
       end if
       
    end if
!# 151 "xc_rVV10.f90"
    ! --------------------------------------------------------------------------------------------------
    ! Allocate arrays.   
    ! ---------------------------------------------------------------------------------------
!# 155 "xc_rVV10.f90"
    allocate( q0(dfftp%nnr) )
    allocate( gradient_rho(3,dfftp%nnr) )
    allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
    allocate( total_rho(dfftp%nnr) )
   
 
    ! ---------------------------------------------------------------------------------------
    ! Add together the valence and core charge densities to get the total charge density    
    !
    total_rho = rho_valence(:) + rho_core(:)
!# 166 "xc_rVV10.f90"
    ! -------------------------------------------------------------------------
    ! Here we calculate the gradient in reciprocal space using FFT
    ! -------------------------------------------------------------------------
    call fft_gradient_r2r( dfftp, total_rho, g, gradient_rho)
!# 171 "xc_rVV10.f90"
    ! -------------------------------------------------------------------------
    ! Get Q and all the derivatives
    ! -------------------------------------------------------------------------
    CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
!# 176 "xc_rVV10.f90"
    ! ---------------------------------------------------------------------------------    
!# 178 "xc_rVV10.f90"
    allocate( thetas(dfftp%nnr, Nqs) )
    CALL get_thetas_on_grid(total_rho, q0, thetas)
    
    call start_clock( 'rVV10_energy')
!# 183 "xc_rVV10.f90"
    call vdW_energy(thetas, Ec_nl)
   
    Ec_nl = Ec_nl + beta * SUM(total_rho) * (omega/(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
    etxc = etxc + Ec_nl
!# 188 "xc_rVV10.f90"
    call stop_clock( 'rVV10_energy')
!# 190 "xc_rVV10.f90"
    ! Print stuff if verbose run
    !
    if (iverbosity > 0) then
!# 194 "xc_rVV10.f90"
       call mp_sum(Ec_nl,intra_bgrp_comm)
       if (ionode) write(*,'(/ / A /)') "     ----------------------------------------------------------------"
       if (ionode) write(*,'(A, F22.15 /)') "     Non-local correlation energy =         ", Ec_nl
       if (ionode) write(*,'(A /)') "     ----------------------------------------------------------------"
       
    end if
!# 201 "xc_rVV10.f90"
    ! ----------------------------------------------------------------------------------------
    ! Inverse Fourier transform the u_i(k) to get the u_i(r) 
    !---------------------------------------------------------------------------------------
!# 205 "xc_rVV10.f90"
    call start_clock( 'rVV10_ffts')
    
    do theta_i = 1, Nqs
       CALL invfft('Rho', thetas(:,theta_i), dfftp) 
    end do
!# 211 "xc_rVV10.f90"
    call stop_clock( 'rVV10_ffts')
!# 213 "xc_rVV10.f90"
    ! -------------------------------------------------------------------------
!# 215 "xc_rVV10.f90"
    call start_clock( 'rVV10_v' )
!# 217 "xc_rVV10.f90"
    allocate( potential(dfftp%nnr) )
    call get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, gradient_rho, thetas, potential)
    
    ! -------------------------------------------------------------------------
    ! Add beta
    ! ------------------------------------------------------------------------- 
    potential = potential + beta
    
    v(:,1) = v(:,1) + potential(:)
    if (nspin==2) v(:,2) = v(:,2) + potential(:) 
!# 228 "xc_rVV10.f90"
    call stop_clock( 'rVV10_v' )
!# 230 "xc_rVV10.f90"
    ! -----------------------------------------------------------------------
    ! The integral of rho(r)*potential(r) for the vtxc output variable
    ! --------------------------------------------------------------------
!# 234 "xc_rVV10.f90"
    grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)  
 
    do i_grid = 1, dfftp%nnr
       vtxc = vtxc + grid_cell_volume*rho_valence(i_grid)*potential(i_grid)
    end do
!# 240 "xc_rVV10.f90"
    deallocate(potential)  
!# 242 "xc_rVV10.f90"
    ! ----------------------------------------------------------------------
    
    ! Deallocate all arrays.
    deallocate(q0, gradient_rho, dq0_drho, dq0_dgradrho, total_rho, thetas)  
!# 247 "xc_rVV10.f90"
    call stop_clock('rVV10')
!# 249 "xc_rVV10.f90"
  END SUBROUTINE xc_rVV10 
!# 252 "xc_rVV10.f90"
  ! #################################################################################################
  !                   |                 |
  !                   |  rVV10_STRESS   |
  !                   |_________________|
!# 257 "xc_rVV10.f90"
  SUBROUTINE rVV10_stress (rho_valence, rho_core, nspin, sigma)
!# 259 "xc_rVV10.f90"
      !! Calculate the stress tensor for rVV10.
!# 261 "xc_rVV10.f90"
      USE fft_base,        ONLY : dfftp
      use gvect,           ONLY : ngm, g
      USE cell_base,       ONLY : tpiba
!# 265 "xc_rVV10.f90"
      implicit none
!# 267 "xc_rVV10.f90"
      real(dp), intent(IN) :: rho_valence(:)
      !! valence charge density
      real(dp), intent(IN) :: rho_core(:)
      !! core charge density
      INTEGER,  INTENT(IN) :: nspin
      !! number of spin components
      real(dp), intent(inout) :: sigma(3,3)
      !! stress tensor
!# 276 "xc_rVV10.f90"
      real(dp), allocatable :: gradient_rho(:,:)         !
      real(dp), allocatable :: total_rho(:)              ! Rho values
!# 279 "xc_rVV10.f90"
      real(dp), allocatable :: q0(:)                     !
      real(dp), allocatable :: dq0_drho(:)               ! Q-values
      real(dp), allocatable :: dq0_dgradrho(:)           !
!# 283 "xc_rVV10.f90"
      complex(dp), allocatable :: thetas(:,:)            ! Thetas
      integer :: i_proc, theta_i, l, m
!# 286 "xc_rVV10.f90"
      real(dp)  :: sigma_grad(3,3)
      real(dp)  :: sigma_ker(3,3)
!# 289 "xc_rVV10.f90"
      ! ---------------------------------------------------------------------------------------------
      !   Tests
      ! --------------------------------------------------------------------------------------------------------
!# 293 "xc_rVV10.f90"
      if (nspin>2) call errore('rV10_stress',' rVV10 stress not implemented for nspin > 2', nspin)
!# 295 "xc_rVV10.f90"
      sigma(:,:) = 0.0_DP
      sigma_grad(:,:) = 0.0_DP
      sigma_ker(:,:) = 0.0_DP
!# 299 "xc_rVV10.f90"
      ! ---------------------------------------------------------------------------------------
      ! Allocations
      ! ---------------------------------------------------------------------------------------
!# 303 "xc_rVV10.f90"
      allocate( gradient_rho(3,dfftp%nnr) )
      allocate( total_rho(dfftp%nnr) )
      allocate( q0(dfftp%nnr) )
      allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
      allocate( thetas(dfftp%nnr, Nqs) )
 
      ! ---------------------------------------------------------------------------------------
      ! Charge
      ! ---------------------------------------------------------------------------------------
      total_rho = rho_valence(:) + rho_core(:)
!# 314 "xc_rVV10.f90"
      ! -------------------------------------------------------------------------
      ! Here we calculate the gradient in reciprocal space using FFT
      ! -------------------------------------------------------------------------
      call fft_gradient_r2r( dfftp, total_rho, g, gradient_rho)
      
      ! -------------------------------------------------------------------------------------------------------------
      ! Get q0.
      ! ---------------------------------------------------------------------------------
!# 323 "xc_rVV10.f90"
      CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
!# 325 "xc_rVV10.f90"
      ! ---------------------------------------------------------------------------------
      ! Get thetas in reciprocal space.
      ! ---------------------------------------------------------------------------------
!# 329 "xc_rVV10.f90"
      CALL get_thetas_on_grid(total_rho, q0, thetas)
!# 331 "xc_rVV10.f90"
      ! ---------------------------------------------------------------------------------------
      ! Stress
      ! ---------------------------------------------------------------------------------------
      CALL rVV10_stress_gradient(total_rho, gradient_rho, q0, dq0_drho, &
                                  dq0_dgradrho, thetas, sigma_grad)
!# 337 "xc_rVV10.f90"
      CALL rVV10_stress_kernel(total_rho, q0, thetas, sigma_ker)
!# 339 "xc_rVV10.f90"
      sigma = - (sigma_grad + sigma_ker) 
!# 341 "xc_rVV10.f90"
      do l = 1, 3
         do m = 1, l - 1
            sigma (m, l) = sigma (l, m)
         enddo
      enddo
!# 347 "xc_rVV10.f90"
      deallocate( gradient_rho, total_rho, q0, dq0_drho, dq0_dgradrho, thetas )
 
   END SUBROUTINE rVV10_stress
!# 351 "xc_rVV10.f90"
   ! ###############################################################################################################
   !                             |                          |
   !                             |  rVV10_stress_gradient   |
!# 355 "xc_rVV10.f90"
   SUBROUTINE rVV10_stress_gradient (total_rho, gradient_rho, q0, dq0_drho, &
                                      dq0_dgradrho, thetas, sigma)
!# 358 "xc_rVV10.f90"
      !! Calculate rVV10 stress with gradient correction.
                                      
      !-----------------------------------------------------------------------------------
      ! Modules to include
      ! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, g, gg, igtongl, &
                                        gl, ngl, gstart
      USE fft_base,              ONLY : dfftp
      USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2
!# 368 "xc_rVV10.f90"
      ! ----------------------------------------------------------------------------------
!# 370 "xc_rVV10.f90"
      implicit none
!# 372 "xc_rVV10.f90"
      real(dp), intent(IN) :: total_rho(:)               !
      real(dp), intent(IN) :: gradient_rho(:, :)         ! Input variables
      real(dp), intent(inout) :: sigma(:,:)              !  
      real(dp), intent(IN) :: q0(:)                      !
      real(dp), intent(IN) :: dq0_drho(:)                ! 
      real(dp), intent(IN) :: dq0_dgradrho(:)            !
      complex(dp), intent(IN) :: thetas(:,:)             !
!# 380 "xc_rVV10.f90"
      complex(dp), allocatable :: u_vdW(:,:)             !
!# 382 "xc_rVV10.f90"
      real(dp), allocatable    :: d2y_dx2(:,:)           !
      real(dp) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation
      real(dp) :: dq                                     !
!# 386 "xc_rVV10.f90"
      integer  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points
!# 388 "xc_rVV10.f90"
      integer  :: l, m
      real(dp) :: prefactor                              ! Final summation of sigma
!# 391 "xc_rVV10.f90"
      integer  :: i_proc, theta_i, i_grid, q_i, &        !
                  ix, iy, iz                             ! Iterators
      
      character(LEN=1) :: intvar
      real(dp) :: const
!# 397 "xc_rVV10.f90"
      !real(dp)       :: at_inverse(3,3)
!# 399 "xc_rVV10.f90"
      allocate( d2y_dx2(Nqs, Nqs) ) 
      allocate( u_vdW(dfftp%nnr, Nqs) )
!# 402 "xc_rVV10.f90"
      const = 1.0D0 / (3.0D0 * b_value**(3.0D0/2.0D0) * pi**(5.0D0/4.0D0) )
      sigma(:,:) = 0.0_DP
      prefactor = 0.0_DP
      
      ! --------------------------------------------------------------------------------------------------
      ! Get u in k-space.
      ! ---------------------------------------------------------------------------------------------------
!# 410 "xc_rVV10.f90"
      call thetas_to_uk(thetas, u_vdW)
!# 412 "xc_rVV10.f90"
      ! --------------------------------------------------------------------------------------------------
      ! Get u in real space.
      ! ---------------------------------------------------------------------------------------------------
!# 416 "xc_rVV10.f90"
      call start_clock( 'rVV10_ffts')
!# 418 "xc_rVV10.f90"
      do theta_i = 1, Nqs
         CALL invfft('Rho', u_vdW(:,theta_i), dfftp) 
      end do
!# 422 "xc_rVV10.f90"
      call stop_clock( 'rVV10_ffts')
!# 424 "xc_rVV10.f90"
      ! --------------------------------------------------------------------------------------------------
      ! Get the second derivatives for interpolating the P_i
      ! ---------------------------------------------------------------------------------------------------
!# 428 "xc_rVV10.f90"
      call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
!# 430 "xc_rVV10.f90"
      ! ---------------------------------------------------------------------------------------------
!# 432 "xc_rVV10.f90"
      i_grid = 0
!# 434 "xc_rVV10.f90"
      ! ----------------------------------------------------------------------------------------------------
      ! Do the real space integration to obtain the stress component
      ! ----------------------------------------------------------------------------------------------------
!# 438 "xc_rVV10.f90"
      do i_grid = 1, dfftp%nnr
!# 440 "xc_rVV10.f90"
                  q_low = 1
                  q_hi = Nqs 
!# 443 "xc_rVV10.f90"
                  !
                  ! Figure out which bin our value of q0 is in in the q_mesh
                  !
                  do while ( (q_hi - q_low) > 1)
!# 448 "xc_rVV10.f90"
                      q = int((q_hi + q_low)/2)
!# 450 "xc_rVV10.f90"
                      if (q_mesh(q) > q0(i_grid)) then
                          q_hi = q
                      else 
                          q_low = q
                      end if
!# 456 "xc_rVV10.f90"
                  end do
!# 458 "xc_rVV10.f90"
                  if (q_hi == q_low) call errore('stress_vdW_gradient','qhi == qlow',1)
!# 460 "xc_rVV10.f90"
                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!# 462 "xc_rVV10.f90"
                  dq = q_mesh(q_hi) - q_mesh(q_low)
!# 464 "xc_rVV10.f90"
                  a = (q_mesh(q_hi) - q0(i_grid))/dq
                  b = (q0(i_grid) - q_mesh(q_low))/dq
                  c = (a**3 - a)*dq**2/6.0D0
                  d = (b**3 - b)*dq**2/6.0D0
                  e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
                  f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
                 
                  do q_i = 1, Nqs
!# 473 "xc_rVV10.f90"
                      y(:) = 0.0D0
                      y(q_i) = 1.0D0
!# 476 "xc_rVV10.f90"
                      dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)
!# 478 "xc_rVV10.f90"
                      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                      
                      if (total_rho(i_grid) > epsr) then
!# 482 "xc_rVV10.f90"
                        prefactor = u_vdW(i_grid,q_i) * const * total_rho(i_grid)**(3.0D0/4.0D0) * dP_dq0 * dq0_dgradrho(i_grid) 
                      
                        do l = 1, 3
                          do m = 1, l
                                        
                              sigma (l, m) = sigma (l, m) -  prefactor * &
                                             (gradient_rho(l,i_grid) * gradient_rho(m,i_grid))
                           enddo
                        enddo
                     endif
!# 493 "xc_rVV10.f90"
                     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!# 495 "xc_rVV10.f90"
                 end do
           
      end do
!# 499 "xc_rVV10.f90"
      call mp_sum(  sigma, intra_bgrp_comm )
!# 501 "xc_rVV10.f90"
      call dscal (9, 1.d0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)
!# 503 "xc_rVV10.f90"
      deallocate( d2y_dx2, u_vdW )
!# 505 "xc_rVV10.f90"
   END SUBROUTINE rVV10_stress_gradient
!# 509 "xc_rVV10.f90"
   ! ###############################################################################################################
   !                      |                        |
   !                      |  rVV10_stress_kernel   |
   !                      |                        |
!# 514 "xc_rVV10.f90"
   SUBROUTINE rVV10_stress_kernel (total_rho, q0, thetas, sigma)
!# 516 "xc_rVV10.f90"
      ! Modules to include
      ! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart 
      USE fft_base,              ONLY : dfftp
      USE cell_base,             ONLY : omega, tpiba, tpiba2
      USE constants, ONLY: pi
!# 523 "xc_rVV10.f90"
      implicit none
     
      real(dp), intent(IN) :: q0(:) 
      real(dp), intent(IN) :: total_rho(:)
      real(dp), intent(inout) :: sigma(3,3)                     !  
      complex(dp), intent(IN) :: thetas(:,:) 
!# 530 "xc_rVV10.f90"
      real(dp), allocatable :: dkernel_of_dk(:,:)               !
      
      integer               :: l, m, q1_i, q2_i , g_i           !
      real(dp)              :: g2, ngmod2, g_kernel, G_multiplier             ! 
      integer               :: last_g, theta_i
!# 536 "xc_rVV10.f90"
      allocate( dkernel_of_dk(Nqs, Nqs) )
!# 538 "xc_rVV10.f90"
      sigma(:,:) = 0.0_DP
!# 540 "xc_rVV10.f90"
      ! --------------------------------------------------------------------------------------------------
      ! Integration in g-space
      ! ---------------------------------------------------------------------------------------------------
!# 544 "xc_rVV10.f90"
      last_g = -1
!# 546 "xc_rVV10.f90"
      G_multiplier = 1.0D0
      if (gamma_only) G_multiplier = 2.0D0
!# 549 "xc_rVV10.f90"
      do g_i = gstart, ngm
!# 551 "xc_rVV10.f90"
          g2 = gg (g_i) * tpiba2
          g_kernel = sqrt(g2)
!# 554 "xc_rVV10.f90"
          if ( igtongl(g_i) .ne. last_g) then
!# 556 "xc_rVV10.f90"
             call interpolate_Dkernel_Dk(g_kernel, dkernel_of_dk)  ! Gets the derivatives
             last_g = igtongl(g_i)
!# 559 "xc_rVV10.f90"
          end if
          
          do q2_i = 1, Nqs
             do q1_i = 1, Nqs
                 do l = 1, 3
                     do m = 1, l
!# 566 "xc_rVV10.f90"
                     sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * &
                                     thetas(dfftp%nl(g_i),q1_i)*dkernel_of_dk(q1_i,q2_i)*conjg(thetas(dfftp%nl(g_i),q2_i))* &
                                     (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel 
                     end do
                 end do 
             enddo
         end do      
!# 574 "xc_rVV10.f90"
         if (g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier
         
      enddo
!# 578 "xc_rVV10.f90"
      call mp_sum(  sigma, intra_bgrp_comm )
      
      deallocate( dkernel_of_dk )
      
   END SUBROUTINE rVV10_stress_kernel
!# 585 "xc_rVV10.f90"
  ! ###############################################################################################################
  !                                    |                  |
  !                                    |  GET_Q0_ON_GRID  |
  !                                    |__________________|
!# 590 "xc_rVV10.f90"
  SUBROUTINE get_q0_on_grid (total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
    
    USE fft_base,        ONLY : dfftp
    
    real(dp),  intent(IN)    :: total_rho(:), gradient_rho(:,:)   
    real(dp),  intent(OUT) :: q0(:), dq0_drho(:), dq0_dgradrho(:)
!# 597 "xc_rVV10.f90"
    integer,    parameter      :: m_cut = 12
!# 599 "xc_rVV10.f90"
    real(dp) :: dw0_dn, dk_dn, gmod2
    real(dp) :: mod_grad, wp2, wg2, w0, k   
    real(dp) :: q, exponent, dq0_dq                
    integer  :: i_grid, index, count=0 
 
    ! initialize q0-related arrays ... 
    q0(:) = q_cut
    dq0_drho(:) = 0.0_DP
    dq0_dgradrho(:) = 0.0_DP
  
    do i_grid = 1, dfftp%nnr
          
     if (total_rho(i_grid) > epsr) then
!# 613 "xc_rVV10.f90"
      gmod2 = gradient_rho(1,i_grid)**2 + &
              gradient_rho(2,i_grid)**2 + &
              gradient_rho(3,i_grid)**2
 
       ! Calculate some intermediate values needed to find q
       ! ------------------------------------------------------------------------------------
       mod_grad = sqrt(gmod2)
!# 621 "xc_rVV10.f90"
       wp2= 16.0_dp*pi*total_rho(i_grid)
       wg2 = 4.0_dp*C_value * (mod_grad/total_rho(i_grid))**4
!# 624 "xc_rVV10.f90"
       k = b_value*3.0_dp*pi* ((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
       w0 = sqrt( wg2 + wp2/3.0_dp  )
!# 627 "xc_rVV10.f90"
       q = w0 / k 
     
       ! Here, we calculate q0 by saturating q according 
       ! ---------------------------------------------------------------------------------------
!# 632 "xc_rVV10.f90"
       exponent = 0.0_dp
       dq0_dq = 0.0_dp
     
       do index = 1, m_cut
        
          exponent = exponent + ( (q/q_cut)**index)/index
          dq0_dq = dq0_dq + ( (q/q_cut)**(index-1))
        
       end do
!# 642 "xc_rVV10.f90"
       q0(i_grid) = q_cut*(1.0_dp - exp(-exponent))
       dq0_dq = dq0_dq * exp(-exponent)
!# 645 "xc_rVV10.f90"
       ! ---------------------------------------------------------------------------------------
!# 647 "xc_rVV10.f90"
       if (q0(i_grid) < q_min) then
         q0(i_grid) = q_min
       end if
!# 651 "xc_rVV10.f90"
       !---------------------------------Final values---------------------------------  
!# 653 "xc_rVV10.f90"
       dw0_dn = 1.0_dp/(2.0_dp*w0) * (16.0_dp/3.0_dp*pi - 4.0_dp*wg2 / total_rho(i_grid) )                
       dk_dn = k / ( 6.0_dp * total_rho(i_grid) )
!# 656 "xc_rVV10.f90"
       dq0_drho(i_grid) = dq0_dq / (k**2) * (dw0_dn * k - dk_dn * w0 )   
       IF ( gmod2 > epsr) THEN
          dq0_dgradrho(i_grid) = dq0_dq  / ( 2.0_dp*k*w0 ) * 4.0_dp*wg2 / (mod_grad**2)
       ELSE
          dq0_dgradrho(i_grid) = 0.0_dp
       ENDIF
!# 663 "xc_rVV10.f90"
     endif  
 
    end do
!# 667 "xc_rVV10.f90"
  end SUBROUTINE get_q0_on_grid
!# 670 "xc_rVV10.f90"
! ###############################################################################################################
!                                      |                      |
!                                      |  GET_THETAS_ON_GRID  |
!# 675 "xc_rVV10.f90"
  SUBROUTINE get_thetas_on_grid (total_rho, q0_on_grid, thetas)
!# 677 "xc_rVV10.f90"
    real(dp), intent(in) :: total_rho(:), q0_on_grid(:)  
!# 679 "xc_rVV10.f90"
    complex(dp), intent(inout):: thetas(:,:)         
!# 681 "xc_rVV10.f90"
    integer :: i_grid, Ngrid_points                
    integer :: theta_i                        
  
    Ngrid_points = size(q0_on_grid)
  
    ! Interpolate the P_i polynomials 
    CALL spline_interpolation(q_mesh, q0_on_grid, thetas)
  
    ! Form the thetas where theta is defined as rho*p_i(q0)
    ! ------------------------------------------------------------------------------------
!# 692 "xc_rVV10.f90"
    do i_grid = 1, Ngrid_points
  
      if (total_rho(i_grid) > epsr ) then 
       thetas(i_grid,:) = thetas(i_grid,:) * (1.0 / (3.0 * sqrt(pi) &
         * ( b_value**(3.0/2.0) ) ) ) * (total_rho(i_grid) / pi)**(3.0/4.0)
      else
       thetas(i_grid,:) = 0.0d0
      endif
     
    end do
!# 703 "xc_rVV10.f90"
    ! ------------------------------------------------------------------------------------
  
    call start_clock( 'rVV10_ffts')
!# 707 "xc_rVV10.f90"
    do theta_i = 1, Nqs
!# 709 "xc_rVV10.f90"
     CALL fwfft ('Rho', thetas(:,theta_i), dfftp)
    end do
!# 712 "xc_rVV10.f90"
    call stop_clock( 'rVV10_ffts')
  
  END SUBROUTINE get_thetas_on_grid
!# 717 "xc_rVV10.f90"
! ###############################################################################################################
!                                     |                        | 
!                                     |  SPLINE_INTERPOLATION  |
!                                     |________________________|
!# 723 "xc_rVV10.f90"
SUBROUTINE spline_interpolation (x, evaluation_points, values)
  
  real(dp), intent(in) :: x(:)
  !! The x values used to form the interpolation
  real(dp), intent(in) :: evaluation_points(:)
  !! (q_mesh in this case) and the values of q0 for which we are 
  !! interpolating the function.
  complex(dp), intent(inout) :: values(:,:)
  !! An output array (allocated outside this routine) that stores the
  !! interpolated values of the P_i (SOLER equation 3) polynomials. The
  !! format is values(grid_point, P_i).
  !
  integer :: Ngrid_points, Nx                                 ! Total number of grid points to evaluate and input x points
  
  real(dp), allocatable, save :: d2y_dx2(:,:)                 ! The second derivatives required to do the interpolation
  
  integer :: i_grid, lower_bound, upper_bound, index, P_i     ! Some indexing variables
  
  real(dp), allocatable :: y(:)                               ! Temporary variables needed for the interpolation
  real(dp) :: a, b, c, d, dx                                  !
 
 
  Nx = size(x)
  Ngrid_points = size(evaluation_points)
  
!# 749 "xc_rVV10.f90"
  ! Allocate the temporary array
  allocate( y(Nx) )
!# 752 "xc_rVV10.f90"
  ! If this is the first time this routine has been called we need to get the second
  ! derivatives (d2y_dx2) required to perform the interpolations.  So we allocate the
  ! array and call initialize_spline_interpolation to get d2y_dx2.
  ! ------------------------------------------------------------------------------------
  if (.not. allocated(d2y_dx2) ) then
!# 758 "xc_rVV10.f90"
     allocate( d2y_dx2(Nx,Nx) )
     call initialize_spline_interpolation(x, d2y_dx2)
     
  end if
!# 763 "xc_rVV10.f90"
  ! ------------------------------------------------------------------------------------
  
  
  do i_grid=1, Ngrid_points
     
     lower_bound = 1
     upper_bound = Nx
     
     do while ( (upper_bound - lower_bound) > 1 )
        
        index = (upper_bound+lower_bound)/2
        
        if ( evaluation_points(i_grid) > x(index) ) then
           lower_bound = index 
        else
           upper_bound = index
        end if
        
     end do
     
     dx = x(upper_bound)-x(lower_bound)
     
     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0
!# 790 "xc_rVV10.f90"
     
     do P_i = 1, Nx
        
        y = 0
        y(P_i) = 1
        
        values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
        
     end do
     
  end do
!# 803 "xc_rVV10.f90"
  deallocate( y )
!# 805 "xc_rVV10.f90"
END SUBROUTINE spline_interpolation
!# 807 "xc_rVV10.f90"
  
! ###############################################################################################################
!                                |                                   |
!                                |  INITIALIZE_SPLINE_INTERPOLATION  |
!                                |___________________________________|
!# 814 "xc_rVV10.f90"
SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)
!# 816 "xc_rVV10.f90"
  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, pages 96-97.  It was adapted for Fortran and for the problem at hand.
!# 819 "xc_rVV10.f90"
  real(dp), intent(in)  :: x(:)
  !! The input abscissa values 
  real(dp), intent(inout) :: d2y_dx2(:,:)
  !! The output array (allocated outside this routine) that holds the second derivatives
  !! required for interpolating the function.
!# 825 "xc_rVV10.f90"
  integer :: Nx, P_i, index                        ! The total number of x points and some indexing
  !                                                ! variables
!# 828 "xc_rVV10.f90"
  real(dp), allocatable :: temp_array(:), y(:)     ! Some temporary arrays required.  y is the array
  !                                                ! that holds the funcion values (all either 0 or 1 here).
!# 831 "xc_rVV10.f90"
  real(dp) :: temp1, temp2                         ! Some temporary variables required
!# 834 "xc_rVV10.f90"
  
  Nx = size(x)
  
  allocate( temp_array(Nx), y(Nx) )
!# 839 "xc_rVV10.f90"
  do P_i=1, Nx
     
!# 842 "xc_rVV10.f90"
     ! In the Soler method, the polynomicals that are interpolated are Kroneker delta funcions
     ! at a particular q point.  So, we set all y values to 0 except the one corresponding to 
     ! the particular function P_i.
     ! ----------------------------------------------------------------------------------------
!# 847 "xc_rVV10.f90"
     y = 0.0D0
     y(P_i) = 1.0D0
!# 850 "xc_rVV10.f90"
     ! ----------------------------------------------------------------------------------------
     
     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0
     
     do index = 2, Nx-1
        
        temp1 = (x(index)-x(index-1))/(x(index+1)-x(index-1))
        temp2 = temp1 * d2y_dx2(P_i,index-1) + 2.0D0
        d2y_dx2(P_i,index) = (temp1-1.0D0)/temp2
        temp_array(index) = (y(index+1)-y(index))/(x(index+1)-x(index)) &
             - (y(index)-y(index-1))/(x(index)-x(index-1))
        temp_array(index) = (6.0D0*temp_array(index)/(x(index+1)-x(index-1)) &
             - temp1*temp_array(index-1))/temp2
        
     end do
     
     d2y_dx2(P_i,Nx) = 0.0D0
     
     do index=Nx-1, 1, -1
        
        d2y_dx2(P_i,index) = d2y_dx2(P_i,index) * d2y_dx2(P_i,index+1) + temp_array(index)
        
     end do
!# 875 "xc_rVV10.f90"
  end do
!# 877 "xc_rVV10.f90"
  deallocate( temp_array, y)
!# 879 "xc_rVV10.f90"
end SUBROUTINE initialize_spline_interpolation
!# 882 "xc_rVV10.f90"
! ###############################################################################################################
!                                         |                    |
!                                         | INTERPOLATE_KERNEL |
!                                         |____________________|
!# 888 "xc_rVV10.f90"
subroutine interpolate_kernel(k, kernel_of_k)
!# 890 "xc_rVV10.f90"
  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, page 97.  Adapted for Fortran and the problem at hand.  This function is used to 
  !! find the Phi_alpha_beta needed for equations 11 and 14 of SOLER.
!# 894 "xc_rVV10.f90"
  real(dp), intent(in) :: k
  !! Input value, the magnitude of the g-vector for the current point.
  
  real(dp), intent(inout) :: kernel_of_k(:,:)
  !! An output array (allocated outside this routine) that holds the interpolated value of
  !! the kernel for each pair of q points (i.e. the phi_alpha_beta of the Soler method.
!# 901 "xc_rVV10.f90"
  integer :: q1_i, q2_i, k_i                    ! Indexing variables
 
  real(dp) :: A, B, C, D                        ! Intermediate values for the interpolation
  
!# 906 "xc_rVV10.f90"
  ! Check to make sure that the kernel table we have is capable of dealing with this
  ! value of k.  If k is larger than Nr_points*2*pi/r_max then we can't perform the 
  ! interpolation.  In that case, a kernel file should be generated with a larger number
  ! of radial points.
  ! -------------------------------------------------------------------------------------
!# 912 "xc_rVV10.f90"
  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if
!# 919 "xc_rVV10.f90"
  ! -------------------------------------------------------------------------------------
  
  kernel_of_k = 0.0D0
  
  ! This integer division figures out which bin k is in since the kernel
  ! is set on a uniform grid.
  k_i = int(k/dk)
  
  ! Test to see if we are trying to interpolate a k that is one of the actual
  ! function points we have.  The value is just the value of the function in that
  ! case.
  ! ----------------------------------------------------------------------------------------
!# 932 "xc_rVV10.f90"
  if (mod(k,dk) == 0) then
     
     do q1_i = 1, Nqs
        do q2_i = 1, q1_i
           
           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)
           
        end do
     end do
     
     return
     
  end if
!# 947 "xc_rVV10.f90"
  ! ----------------------------------------------------------------------------------------
!# 950 "xc_rVV10.f90"
  ! If we are not on a function point then we carry out the interpolation
  ! ----------------------------------------------------------------------------------------
  
  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
  C = (A**3-A)*dk**2/6.0D0
  D = (B**3-B)*dk**2/6.0D0
  
  do q1_i = 1, Nqs
     do q2_i = 1, q1_i
        
        kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
             +(C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))
        
        kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)
        
     end do
  end do
!# 969 "xc_rVV10.f90"
  ! ----------------------------------------------------------------------------------------
!# 971 "xc_rVV10.f90"
  
end subroutine interpolate_kernel
!# 975 "xc_rVV10.f90"
! ###############################################################################################################
!                                         |                        |
!                                         | INTERPOLATE_DKERNEL_DK |
!                                         |________________________|
!# 982 "xc_rVV10.f90"
subroutine interpolate_Dkernel_Dk(k, dkernel_of_dk)
  
  implicit none 
!# 986 "xc_rVV10.f90"
  real(dp), intent(in) :: k        
  real(dp), intent(inout) :: dkernel_of_dk(Nqs,Nqs)
!# 989 "xc_rVV10.f90"
  integer :: q1_i, q2_i, k_i  
  real(dp) :: A, B, dAdk, dBdk, dCdk, dDdk    
  
!# 993 "xc_rVV10.f90"
  ! -------------------------------------------------------------------------------------
!# 995 "xc_rVV10.f90"
  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if
!# 1002 "xc_rVV10.f90"
  ! -------------------------------------------------------------------------------------
!# 1004 "xc_rVV10.f90"
  dkernel_of_dk = 0.0D0
!# 1006 "xc_rVV10.f90"
  k_i = int(k/dk)
!# 1008 "xc_rVV10.f90"
  ! ----------------------------------------------------------------------------------------
!# 1010 "xc_rVV10.f90"
  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
!# 1013 "xc_rVV10.f90"
  dAdk = -1.0D0/dk
  dBdk = 1.0D0/dk
  dCdk = -((3*A**2 -1.0D0)/6.0D0)*dk
  dDdk = ((3*B**2 -1.0D0)/6.0D0)*dk
!# 1018 "xc_rVV10.f90"
  do q1_i = 1, Nqs
     do q2_i = 1, q1_i
!# 1021 "xc_rVV10.f90"
        dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
             + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)
!# 1024 "xc_rVV10.f90"
        dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)
!# 1026 "xc_rVV10.f90"
     end do
  end do
!# 1029 "xc_rVV10.f90"
  ! ----------------------------------------------------------------------------------------
!# 1031 "xc_rVV10.f90"
  
end subroutine interpolate_Dkernel_Dk 
!# 1034 "xc_rVV10.f90"
! #################################################################################################
!                                          |              |
!                                          | thetas_to_uk |
!                                          |______________|
!# 1040 "xc_rVV10.f90"
subroutine thetas_to_uk(thetas, u_vdW)
  
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega
!# 1046 "xc_rVV10.f90"
  complex(dp), intent(in) :: thetas(:,:)  
  complex(dp), intent(out) :: u_vdW(:,:)
!# 1049 "xc_rVV10.f90"
  real(dp), allocatable :: kernel_of_k(:,:)    
!# 1051 "xc_rVV10.f90"
  real(dp) :: g
  integer :: last_g, g_i, q1_i, q2_i, count, i_grid   
!# 1054 "xc_rVV10.f90"
  complex(dp) :: theta(Nqs)        
  
  ! -------------------------------------------------------------------------------------------------
  
  allocate( kernel_of_k(Nqs, Nqs) )
!# 1060 "xc_rVV10.f90"
  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP, kind=dp)
  
  last_g = -1 
!# 1064 "xc_rVV10.f90"
  do g_i = 1, ngm
    
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(dfftp%nl(g_i),:)
     
     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i) = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
     end do
!# 1082 "xc_rVV10.f90"
  end do
!# 1084 "xc_rVV10.f90"
  if (gamma_only) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))
  
  deallocate( kernel_of_k )
     
  ! -----------------------------------------------------------------------------------------------
  
end subroutine thetas_to_uk
!# 1092 "xc_rVV10.f90"
! #################################################################################################
!                                              |             |
!                                              | VDW_ENERGY  |
!                                              |_____________|
!# 1097 "xc_rVV10.f90"
subroutine vdW_energy(thetas, vdW_xc_energy)
  
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega
!# 1103 "xc_rVV10.f90"
  complex(dp), intent(inout) :: thetas(:,:)  
  real(dp), intent(out) :: vdW_xc_energy     
!# 1106 "xc_rVV10.f90"
  real(dp), allocatable :: kernel_of_k(:,:)   
!# 1108 "xc_rVV10.f90"
  real(dp) :: g                      
  integer  :: last_g               
  
  integer :: g_i, q1_i, q2_i, count, i_grid  
!# 1113 "xc_rVV10.f90"
  complex(dp) :: theta(Nqs), thetam(Nqs), theta_g(Nqs)  
  real(dp)    :: G0_term, G_multiplier 
!# 1116 "xc_rVV10.f90"
  complex(dp), allocatable :: u_vdw(:,:)   
!# 1118 "xc_rVV10.f90"
  vdW_xc_energy = 0.0D0
 
  allocate (u_vdW(dfftp%nnr,Nqs))
  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP, kind=dp)
!# 1123 "xc_rVV10.f90"
  allocate( kernel_of_k(Nqs, Nqs) )
  
  !
  ! Here we should use gstart,ngm but all the cases are handeld by conditionals inside the loop
  !
  G_multiplier = 1.0D0
  if (gamma_only) G_multiplier = 2.0D0
!# 1131 "xc_rVV10.f90"
  last_g = -1 
!# 1133 "xc_rVV10.f90"
  do g_i = 1, ngm
     
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(dfftp%nl(g_i),:)
!# 1145 "xc_rVV10.f90"
     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i)  = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(dfftp%nl(g_i),q2_i)*conjg(theta(q2_i)))
     end do
     
     if (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier
!# 1154 "xc_rVV10.f90"
  end do
!# 1156 "xc_rVV10.f90"
  if (gamma_only) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))
!# 1158 "xc_rVV10.f90"
  ! Final value 
  vdW_xc_energy = 0.5D0 * omega * vdW_xc_energy   
  
  deallocate( kernel_of_k )
  thetas(:,:) = u_vdW(:,:)
  deallocate (u_vdW)
  ! ---------------------------------------------------------------------------------------------------
  
end subroutine vdW_energy
!# 1169 "xc_rVV10.f90"
! ###############################################################################################################
!                                             |                 |
!                                             |  GET_POTENTIAL  |
!                                             |_________________|
!# 1174 "xc_rVV10.f90"
  subroutine get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, gradient_rho, u_vdW, potential)
!# 1176 "xc_rVV10.f90"
    use gvect,               ONLY : g
    USE fft_base,            ONLY : dfftp
    USE cell_base,           ONLY : alat, tpiba
!# 1180 "xc_rVV10.f90"
    real(dp), intent(in) ::  q0(:), gradient_rho(:,:)   
    real(dp), intent(in) :: dq0_drho(:), dq0_dgradrho(:)
    real(dp), intent(in) :: total_rho(:)
    complex(dp), intent(in) :: u_vdW(:,:)    
    real(dp), intent(inout) :: potential(:)
  
    real(dp), allocatable, save :: d2y_dx2(:,:)         
!# 1188 "xc_rVV10.f90"
    integer :: i_grid, P_i,icar                         
    integer :: q_low, q_hi, q                           
    real(dp) :: dq, a, b, c, d, e, f                    
    real(dp) :: y(Nqs), dP_dq0, P                       
    !                                                   
!# 1194 "xc_rVV10.f90"
    real(dp), allocatable ::h_prefactor(:)
    complex(dp), allocatable ::h(:)
    real(dp) :: dtheta_dn, dtheta_dgradn
!# 1198 "xc_rVV10.f90"
    real(dp) :: const 
!# 1201 "xc_rVV10.f90"
    allocate (h_prefactor(dfftp%nnr),h(dfftp%nnr))
!# 1203 "xc_rVV10.f90"
    const = 1.0D0 / (3.0D0 * b_value**(3.0D0/2.0D0) * pi**(5.0D0/4.0D0) )
    potential = 0.0D0
    h_prefactor   = 0.0D0
!# 1207 "xc_rVV10.f90"
    ! -------------------------------------------------------------------------------------------
    ! Get the second derivatives of the P_i functions for interpolation
    ! ---------------------------------------------------------------------------------------------
!# 1211 "xc_rVV10.f90"
    if (.not. allocated( d2y_dx2) ) then
     
     allocate( d2y_dx2(Nqs, Nqs) )
     call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
     
    end if
  
    ! ---------------------------------------------------------------------------------------------
  
!# 1221 "xc_rVV10.f90"
    do i_grid = 1,dfftp%nnr
           
     q_low = 1
     q_hi = Nqs 
!# 1226 "xc_rVV10.f90"
     ! Figure out which bin our value of q0 is in in the q_mesh
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!# 1229 "xc_rVV10.f90"
     do while ( (q_hi - q_low) > 1)
              
        q = int((q_hi + q_low)/2)
              
        if (q_mesh(q) > q0(i_grid)) then
           q_hi = q
        else 
           q_low = q
        end if
              
     end do
           
     if (q_hi == q_low) call errore('get_potential','qhi == qlow',1)
           
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!# 1245 "xc_rVV10.f90"
     dq = q_mesh(q_hi) - q_mesh(q_low)
           
     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
           
     do P_i = 1, Nqs
        y = 0.0D0
        y(P_i) = 1.0D0
              
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
        P = a*y(q_low) + b*y(q_hi) + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
 
!# 1262 "xc_rVV10.f90"
        ! IF THE CHARGE DENSITY IS NEGATIVE WE PUT POTENTIAL = 0, OUTSIDE THE SUBROUTINE WE ADD BETA. 
        if (total_rho(i_grid) > epsr) then
!# 1265 "xc_rVV10.f90"
          dtheta_dn = const * (3.0D0/4.0D0) / (total_rho(i_grid)**(1.0D0/4.0D0))  * P + &
                      const * total_rho(i_grid)**(3.0D0/4.0D0) *  dP_dq0 * dq0_drho(i_grid)
          dtheta_dgradn = const * total_rho(i_grid)**(3.0D0/4.0D0) * dP_dq0 * dq0_dgradrho(i_grid)
!# 1269 "xc_rVV10.f90"
          potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* dtheta_dn 
!# 1271 "xc_rVV10.f90"
          if (q0(i_grid) .ne. q_mesh(Nqs)) then
            h_prefactor(i_grid) = h_prefactor(i_grid) +  u_vdW(i_grid,P_i)* dtheta_dgradn
          end if
!# 1275 "xc_rVV10.f90"
        end if
!# 1277 "xc_rVV10.f90"
     end do
    end do
!# 1280 "xc_rVV10.f90"
    do icar = 1,3
      h(:) = CMPLX( h_prefactor(:)*gradient_rho(icar,:), 0.0_DP, kind=dp)
      CALL fwfft ('Rho', h, dfftp) 
      h(dfftp%nl(:)) = CMPLX(0.0_DP,1.0_DP,kind=dp)*tpiba*g(icar,:)*h(dfftp%nl(:))
      if (gamma_only) h(dfftp%nlm(:)) = CONJG(h(dfftp%nl(:)))
      CALL invfft ('Rho', h, dfftp) 
      potential(:) = potential(:) - REAL(h(:))
    end do
!# 1289 "xc_rVV10.f90"
    ! ------------------------------------------------------------------------------------------------------------------------
    deallocate (h_prefactor,h)
!# 1292 "xc_rVV10.f90"
  end subroutine get_potential
!# 1295 "xc_rVV10.f90"
! ###############################################################################################################
!                                                 |                   |
!                                                 |  generate_kernel  |
!                                                 |___________________|
SUBROUTINE generate_kernel
!# 1301 "xc_rVV10.f90"
  implicit none
!# 1303 "xc_rVV10.f90"
  integer  :: q1_i, q2_i, r_i                   ! Indexing variables
  real(dp) :: d1, d2                            ! Intermediate values
!# 1307 "xc_rVV10.f90"
  kernel    = 0.0D0
  d2phi_dk2 = 0.0D0
!# 1310 "xc_rVV10.f90"
  do q1_i = 1, Nqs
     do q2_i = 1, q1_i
!# 1313 "xc_rVV10.f90"
        do r_i = 1, Nr_points
           d1 = q_mesh(q1_i) * (dr * r_i)**2    ! Different definition of d1 and d2 for vv10
           d2 = q_mesh(q2_i) * (dr * r_i)**2    ! Different definition of d1 and d2 for vv10
           kernel(r_i, q1_i, q2_i) = -24.0D0 / ( ( d1+1.0 ) * ( d2+1.0 ) * ( d1+d2+2.0 ) )
        end do
!# 1319 "xc_rVV10.f90"
        call radial_fft( kernel(:, q1_i, q2_i) )
        call set_up_splines( kernel(:, q1_i, q2_i), d2phi_dk2(:, q1_i, q2_i) )
!# 1322 "xc_rVV10.f90"
        kernel    (:, q2_i, q1_i) = kernel   (:, q1_i, q2_i)
        d2phi_dk2 (:, q2_i, q1_i) = d2phi_dk2(:, q1_i, q2_i)
!# 1325 "xc_rVV10.f90"
     end do
  end do
!# 1328 "xc_rVV10.f90"
END SUBROUTINE generate_kernel
!# 1331 "xc_rVV10.f90"
! ###############################################################################################################
!                                                 |              |
!                                                 |  radial_fft  |
!                                                 |______________|
SUBROUTINE radial_fft(phi)
!# 1337 "xc_rVV10.f90"
  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  REAL(DP)                :: phi_k(0:Nr_points)
  INTEGER                 :: k_i, r_i
  REAL(DP)                :: r, k
  
  
  phi_k = 0.0D0
  
  DO r_i = 1, Nr_points
     r        = r_i * dr
     phi_k(0) = phi_k(0) + phi(r_i)*r**2
  END DO
  
  phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)
  
  DO k_i = 1, Nr_points
     k = k_i * dk
     DO r_i = 1, Nr_points
        r          = r_i * dr
        phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * SIN(k*r) / k
     END DO
     phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r * SIN(k*r) / k
  END DO
  
  phi = 4.0D0 * pi * phi_k * dr
!# 1363 "xc_rVV10.f90"
END SUBROUTINE radial_fft
!# 1366 "xc_rVV10.f90"
! ###############################################################################################################
!                                              |                  |
!                                              |  set_up_splines  |
!                                              |__________________|
SUBROUTINE set_up_splines(phi, D2)
!# 1372 "xc_rVV10.f90"
  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  REAL(DP), ALLOCATABLE   :: temp_array(:)
  REAL(DP)                :: temp_1, temp_2
  INTEGER                 :: r_i
  
  
  ALLOCATE( temp_array(0:Nr_points) )
  
  D2         = 0
  temp_array = 0
  
  DO r_i = 1, Nr_points - 1
     temp_1  = DBLE(r_i - (r_i - 1))/DBLE( (r_i + 1) - (r_i - 1) )
     temp_2  = temp_1 * D2(r_i-1) + 2.0D0
     D2(r_i) = (temp_1 - 1.0D0)/temp_2
     temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/DBLE( dk*((r_i+1) - r_i) ) - &
          ( phi(r_i) - phi(r_i-1))/DBLE( dk*(r_i - (r_i-1)) )
     temp_array(r_i) = (6.0D0*temp_array(r_i)/DBLE( dk*((r_i+1) - (r_i-1)) )-&
          temp_1*temp_array(r_i-1))/temp_2
  END DO
  
  D2(Nr_points) = 0.0D0
  DO  r_i = Nr_points-1, 0, -1
     D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)
  END DO
  
  DEALLOCATE( temp_array )
!# 1401 "xc_rVV10.f90"
END SUBROUTINE set_up_splines
!# 1404 "xc_rVV10.f90"
END MODULE rVV10
