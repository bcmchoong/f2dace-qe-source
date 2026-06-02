!# 1 "xc_vdW_DF.f90"
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2019 Brian Kolb, Timo Thonhauser
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! ----------------------------------------------------------------------
!# 10 "xc_vdW_DF.f90"
MODULE vdW_DF
!# 12 "xc_vdW_DF.f90"
!! This module calculates the non-local correlation contribution to the
!! energy and potential according to:
!!
!! * M. Dion et al., Phys. Rev. Lett. 92, 246401 (2004),
!!   <https://doi.org/10.1103/PhysRevLett.92.246401>
!!
!! henceforth referred to as DION. Further information about the
!! functional and its corresponding potential can be found in:
!!
!! * T. Thonhauser et al., Phys. Rev. B 76, 125112 (2007),
!!   <https://doi.org/10.1103/PhysRevB.76.125112>
!!
!! The proper spin extension of vdW-DF, i.e. svdW-DF, is derived in:
!!
!! * T. Thonhauser et al., Phys. Rev. Lett. 115, 136402 (2015),
!!   <https://doi.org/10.1103/PhysRevLett.115.136402>
!!
!! henceforth referred to as THONHAUSER. Two review articles show many
!! of the vdW-DF applications:
!!
!! * D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009),
!!   <https://10.1088/0953-8984/21/8/084203>
!!
!! * K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015),
!!   <https://10.1088/0034-4885/78/6/066501>
!!
!! The method implemented is based on the method of G. Roman-Perez and
!! J.M. Soler described in:
!!
!! * G. Roman-Perez and J.M. Soler, Phys. Rev. Lett. 103, 096102 (2009),
!!   <https://doi.org/10.1103/PhysRevLett.103.096102>
!!
!! henceforth referred to as SOLER.
!!
!!
!! Some of the algorithms used in this module are somewhat modified
!! versions of those found in:
!!
!! * Numerical Recipes in C; William H. Press, Brian P. Flannery, Saul
!!   A. Teukolsky, and William T. Vetterling. Cambridge University
!!   Press (1988).
!!
!! hereafter referred to as NUMERICAL_RECIPES. The routines were
!! translated to Fortran, of course, and variable names are generally
!! different.
!!
!!
!! `xc_vdW_DF` and `xc_vdW_DF_spin` are the driver routines for vdW-DF
!! calculations and are called from `Modules/funct.f90.` The routines in
!! this module set up the parallel run (if any) and carry out the calls
!! necessary to calculate the non-local correlation contributions to the
!! energy and potential.
!
!
!  Other files relevant for vdW-DF are:
!
!  * Modules/funct.f90: Definition of functional names
!
!  * XClib/qe_drivers_gga.f90: Driver routines for XC functionals
!
!  * XClib/qe_funct_exch_gga.f90: Code for the exchange
!
!  * XClib/qe_dft_refs.f90: References for all functionals
!
!  * XClib/qe_dft_list.f90: List of all functionals
!# 79 "xc_vdW_DF.f90"
USE kinds,             ONLY : dp
USE constants,         ONLY : pi, fpi, e2
USE mp,                ONLY : mp_sum, mp_barrier, mp_get, mp_size, mp_rank, mp_bcast
USE mp_images,         ONLY : intra_image_comm
USE mp_bands,          ONLY : intra_bgrp_comm
USE io_global,         ONLY : stdout, ionode
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fwfft, invfft
USE control_flags,     ONLY : iverbosity, gamma_only
USE corr_lda,          ONLY : pw, pw_spin
!# 91 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! No implicit variables
!# 94 "xc_vdW_DF.f90"
IMPLICIT NONE
!# 97 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! By default everything is private
!# 100 "xc_vdW_DF.f90"
PRIVATE
!# 103 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! Save all objects in this module
!# 106 "xc_vdW_DF.f90"
SAVE
!# 109 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! Public functions
!# 112 "xc_vdW_DF.f90"
PUBLIC  :: xc_vdW_DF, xc_vdW_DF_spin, vdW_DF_stress,                   &
           vdW_DF_energy, vdW_DF_potential,                            &
           generate_kernel, interpolate_kernel,                        &
           initialize_spline_interpolation, spline_interpolation
!# 118 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! Public variables
!# 121 "xc_vdW_DF.f90"
PUBLIC  :: inlc, vdW_DF_analysis, Nr_points, r_max, q_min, q_cut, Nqs, q_mesh
!# 124 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! General variables
!# 127 "xc_vdW_DF.f90"
INTEGER                  :: inlc            = 1
! The non-local correlation
!# 130 "xc_vdW_DF.f90"
INTEGER                  :: vdW_DF_analysis = 0
! vdW-DF analysis tool as described in PRB 97, 085115 (2018)
!# 133 "xc_vdW_DF.f90"
REAL(DP), PARAMETER      :: epsr            = 1.0D-12
! A small number to cut off densities
!# 136 "xc_vdW_DF.f90"
INTEGER                  :: idx
! Indexing variable
!# 140 "xc_vdW_DF.f90"
! ----------------------------------------------------------------------
! Kernel specific parameters and variables
!# 143 "xc_vdW_DF.f90"
INTEGER, PARAMETER       :: Nr_points = 1024
! The number of radial points (also the number of k points) used in the
! formation of the kernel functions for each pair of q values.
! Increasing this value will help in case you get a run-time error
! saying that you are trying to use a k value that is larger than the
! largest tabulated k point since the largest k point will be 2*pi/r_max
! * Nr_points. Memory usage of the vdW_DF piece of PWSCF will increase
! roughly linearly with this variable.
!# 152 "xc_vdW_DF.f90"
REAL(DP), PARAMETER      :: r_max     = 100.0D0
! The value of the maximum radius to use for the real-space kernel
! functions for each pair of q values. The larger this value is the
! smaller the smallest k value will be since the smallest k point value
! is 2*pi/r_max. Be careful though, since this will also decrease the
! maximum k point value and the vdW_DF code will crash if it encounters
! a g-vector with a magnitude greater than 2*pi/r_max *Nr_points.
!# 160 "xc_vdW_DF.f90"
REAL(DP), PARAMETER      :: dr = r_max/Nr_points, dk = 2.0D0*pi/r_max
! Real space and k-space spacing of grid points.
!# 163 "xc_vdW_DF.f90"
REAL(DP), PARAMETER      :: q_min = 1.0D-5, q_cut = 5.0D0
! The maximum and minimum values of q. During a vdW run, values of q0
! found larger than q_cut will be saturated (SOLER equation 5) to q_cut.
!# 167 "xc_vdW_DF.f90"
INTEGER,  PARAMETER                 :: Nqs    = 20
REAL(DP), PARAMETER, DIMENSION(Nqs) :: q_mesh = (/  &
   q_min              , 0.0449420825586261D0, 0.0975593700991365D0, 0.159162633466142D0, &
   0.231286496836006D0, 0.315727667369529D0 , 0.414589693721418D0 , 0.530335368404141D0, &
   0.665848079422965D0, 0.824503639537924D0 , 1.010254382520950D0 , 1.227727621364570D0, &
   1.482340921174910D0, 1.780437058359530D0 , 2.129442028133640D0 , 2.538050036534580D0, &
   3.016440085356680D0, 3.576529545442460D0 , 4.232271035198720D0 , q_cut /)
!# 175 "xc_vdW_DF.f90"
! The above two parameters define the q mesh to be used in the vdW_DF
! code. These are perhaps the most important to have set correctly.
! Increasing the number of q points will DRAMATICALLY increase the
! memory usage of the vdW_DF code because the memory consumption depends
! quadratically on the number of q points in the mesh. Increasing the
! number of q points may increase accuracy of the vdW_DF code, although,
! in testing it was found to have little effect. The largest value of
! the q mesh is q_cut. All values of q0 (DION equation 11) larger than
! this value during a run will be saturated to this value using equation
! 5 of SOLER. In testing, increasing the value of q_cut was found to
! have little impact on the results, although it is possible that in
! some systems it may be more important. Always make sure that the
! variable Nqs is consistent with the number of q points that are
! actually in the variable q_mesh. Also, do not set any q value to 0.
! This will cause an infinity in the Fourier transform.
!# 191 "xc_vdW_DF.f90"
INTEGER,  PARAMETER      :: Nintegration_points = 256
! Number of integration points for real-space kernel generation (see
! DION equation 14). This is how many a's and b's there will be.
!# 195 "xc_vdW_DF.f90"
REAL(DP), PARAMETER      :: a_min = 0.0D0, a_max = 64.0D0
! Min/max values for the a and b integration in DION equation 14.
!# 198 "xc_vdW_DF.f90"
REAL(DP) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )
! Matrices holding the Fourier transformed kernel function  and its
! second derivative for each pair of q values. The ordering is
! kernel(k_point, q1_value, q2_value).
!# 203 "xc_vdW_DF.f90"
REAL(DP) :: W_ab( Nintegration_points, Nintegration_points )
! Defined in DION equation 16.
!# 206 "xc_vdW_DF.f90"
REAL(DP) :: a_points( Nintegration_points ), a_points2( Nintegration_points )
! The values of the "a" points (DION equation 14) and their squares.
!# 209 "xc_vdW_DF.f90"
CONTAINS
!# 218 "xc_vdW_DF.f90"
  ! ####################################################################
  !                           |             |
  !                           |  functions  |
  !                           |_____________|
  !
  ! Functions to be used in get_q0_on_grid, get_q0_on_grid_spin, and
  ! phi_value.
!# 227 "xc_vdW_DF.f90"
  FUNCTION Fs(s)
!# 229 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP) :: s, Fs, Z_ab = 0.0D0
!# 232 "xc_vdW_DF.f90"
     IF ( inlc == 1 .OR. inlc == 3 ) THEN
        Z_ab = -0.8491D0
     ELSE IF ( inlc == 2 .OR. inlc == 4 .OR. inlc == 5 .OR. inlc == 6 ) THEN
        Z_ab = -1.887D0
     END IF
!# 238 "xc_vdW_DF.f90"
     Fs = 1.0D0 - Z_ab * s * s / 9.0D0
!# 240 "xc_vdW_DF.f90"
  END FUNCTION Fs
!# 245 "xc_vdW_DF.f90"
  FUNCTION dFs_ds(s)
!# 247 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP)             :: s, dFs_ds, Z_ab = 0.0D0
     REAL(DP), PARAMETER  :: prefac = -2.0D0/9.0D0
!# 251 "xc_vdW_DF.f90"
     IF ( inlc == 1 .OR. inlc == 3 ) THEN
        Z_ab = -0.8491D0
     ELSE IF ( inlc == 2 .OR. inlc == 4 .OR. inlc == 5 .OR. inlc == 6 ) THEN
        Z_ab = -1.887D0
     END IF
!# 257 "xc_vdW_DF.f90"
     dFs_ds =  prefac * s * Z_ab
!# 259 "xc_vdW_DF.f90"
  END FUNCTION dFs_ds
!# 264 "xc_vdW_DF.f90"
  FUNCTION kF(rho)
!# 266 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP)             :: rho, kF
     REAL(DP), PARAMETER  :: ex = 1.0D0/3.0D0
!# 270 "xc_vdW_DF.f90"
     kF = ( 3.0D0 * pi * pi * rho )**ex
!# 272 "xc_vdW_DF.f90"
  END FUNCTION kF
!# 277 "xc_vdW_DF.f90"
  FUNCTION dkF_drho(rho)
!# 279 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP)             :: rho, dkF_drho
     REAL(DP), PARAMETER  :: prefac = 1.0D0/3.0D0
!# 283 "xc_vdW_DF.f90"
     dkF_drho = prefac * kF(rho) / rho
!# 285 "xc_vdW_DF.f90"
  END FUNCTION dkF_drho
!# 290 "xc_vdW_DF.f90"
  FUNCTION ds_drho(rho, s)
!# 292 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP) :: rho, s, ds_drho
!# 295 "xc_vdW_DF.f90"
     ds_drho = -s * ( dkF_drho(rho) / kF(rho) + 1.0D0 / rho )
!# 297 "xc_vdW_DF.f90"
  END FUNCTION ds_drho
!# 302 "xc_vdW_DF.f90"
  FUNCTION ds_dgradrho(rho)
!# 304 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP) :: rho, ds_dgradrho
!# 307 "xc_vdW_DF.f90"
     ds_dgradrho = 0.5D0 / (kF(rho) * rho)
!# 309 "xc_vdW_DF.f90"
  END FUNCTION ds_dgradrho
!# 314 "xc_vdW_DF.f90"
  FUNCTION dqx_drho(rho, s)
!# 316 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP) :: rho, s, dqx_drho
!# 319 "xc_vdW_DF.f90"
     dqx_drho = dkF_drho(rho) * Fs(s) + kF(rho) * dFs_ds(s) * ds_drho(rho, s)
!# 321 "xc_vdW_DF.f90"
  END FUNCTION dqx_drho
!# 326 "xc_vdW_DF.f90"
  FUNCTION h_function(y)
!# 328 "xc_vdW_DF.f90"
     IMPLICIT NONE
     REAL(DP)             :: y, y2, y4, h_function
     REAL(DP), PARAMETER  :: g1 = fpi/9.0D0                                     ! vdW-DF1/2
     REAL(DP), PARAMETER  :: a3 = 0.94950D0, g3 = 1.12D0, g32 = g3*g3           ! vdW-DF3-opt1
     REAL(DP), PARAMETER  :: a4 = 0.28248D0, g4 = 1.29D0, g42 = g4*g4           ! vdW-DF3-opt2
     REAL(DP), PARAMETER  :: a5 = 2.01059D0, b5 = 8.17471D0, g5 = 1.84981D0, &  ! vdW-DF-C6
                             AA = ( b5 + a5*(a5/2.0D0-g5) ) / ( 1.0D0+g5-a5 )   !
     REAL(DP), PARAMETER  :: a6 = 0.0532D0,  g6 = 1.42D0, g62 = g6*g6           ! vdW-DF3-mc
!# 338 "xc_vdW_DF.f90"
     y2 = y*y
!# 340 "xc_vdW_DF.f90"
     IF ( inlc == 1 .OR. inlc == 2 ) THEN
!# 342 "xc_vdW_DF.f90"
        h_function = 1.0D0 - EXP( -g1*y2 )
!# 344 "xc_vdW_DF.f90"
     ELSE IF ( inlc == 3 ) THEN
!# 346 "xc_vdW_DF.f90"
        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g3*y2 + g32*y4 + a3*y4*y4 )
!# 349 "xc_vdW_DF.f90"
     ELSE IF ( inlc == 4 ) THEN
!# 351 "xc_vdW_DF.f90"
        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g4*y2 + g42*y4 + a4*y4*y4 )
!# 354 "xc_vdW_DF.f90"
     ELSE IF ( inlc == 5 ) THEN
!# 356 "xc_vdW_DF.f90"
        y4 = y2*y2
        h_function = 1.0D0 - ( 1.0D0 + ( (a5-g5)*y2 + AA*y4 ) / ( 1.0D0+AA*y2 ) ) * EXP( -a5*y2 )
!# 359 "xc_vdW_DF.f90"
     ELSE IF ( inlc == 6 ) THEN
!# 361 "xc_vdW_DF.f90"
        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g6*y2 + g62*y4 + a6*y4*y4 )
!# 364 "xc_vdW_DF.f90"
     END IF
!# 366 "xc_vdW_DF.f90"
  END FUNCTION
!# 375 "xc_vdW_DF.f90"
  ! ####################################################################
  !                           |             |
  !                           |  XC_VDW_DF  |
  !                           |_____________|
!# 380 "xc_vdW_DF.f90"
  SUBROUTINE xc_vdW_DF (rho_valence, rho_core, etxc, vtxc, v)
!# 382 "xc_vdW_DF.f90"
  !! Driver routine for vdW-DF calculations, called from `Modules/funct.f90`.
  !! The routine here sets up the parallel run (if any) and carres out the
  !! calls necessary to calculate the non-local correlation contributions
  !! to the energy and potential. This routine handles the `nspin=1` case,
  !! the `nspin=2` case is handled by `xc_vdW_DF_spin`.
!# 389 "xc_vdW_DF.f90"
  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
!# 392 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 394 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Local variables
!# 397 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: rho_valence(:,:)    !
  REAL(DP), INTENT(IN)    :: rho_core(:)         !  PWSCF input variables
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:)  !
!# 401 "xc_vdW_DF.f90"
  INTEGER :: i_grid, theta_i, i_proc             ! Indexing variables over grid points,
                                                 ! theta functions, and processors.
!# 404 "xc_vdW_DF.f90"
  REAL(DP) :: grid_cell_volume                   ! The volume of the unit cell per G-grid point.
!# 406 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE ::  q0(:)                ! The saturated value of q (equations 11 and 12
                                                 ! of DION). This saturation is that of
                                                 ! equation 5 in SOLER.
!# 410 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: grad_rho(:,:)         ! The gradient of the charge density. The
                                                 ! format is as follows:
                                                 ! grad_rho(cartesian_component,grid_point).
!# 414 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: potential(:)          ! The vdW contribution to the potential.
!# 416 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: dq0_drho(:)           ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the charge density (see
                                                 ! get_q0_on_grid subroutine for details).
!# 421 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: dq0_dgradrho(:)       ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the gradient of the charge density
                                                 ! (again, see get_q0_on_grid subroutine).
!# 426 "xc_vdW_DF.f90"
  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)        ! These are the functions of equation 8 of
                                                 ! SOLER. They will be forward Fourier transformed
                                                 ! in place to get theta(k) and worked on in
                                                 ! place to get the u_alpha(r) of equation 11
                                                 ! in SOLER. They are formatted as follows:
                                                 ! thetas(grid_point, theta_i).
!# 433 "xc_vdW_DF.f90"
  REAL(DP) :: Ec_nl                              ! The non-local vdW contribution to the energy.
!# 435 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: total_rho(:)          ! This is the sum of the valence and core
                                                 ! charge. This just holds the piece assigned
                                                 ! to this processor.
!# 439 "xc_vdW_DF.f90"
  LOGICAL, SAVE :: first_iteration = .TRUE.      ! Whether this is the first time this
                                                 ! routine has been called.
!# 445 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Write out the vdW-DF information and initialize the calculation.
!# 448 "xc_vdW_DF.f90"
  IF ( first_iteration ) THEN
     IF ( inlc > 6 ) CALL errore( 'xc_vdW_DF', 'inlc not implemented', 1 )
     CALL generate_kernel
     IF ( ionode ) CALL vdW_info
     first_iteration = .FALSE.
  END IF
!# 456 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.
!# 460 "xc_vdW_DF.f90"
  allocate( total_rho(dfftp%nnr), grad_rho(3,dfftp%nnr),                &
            potential(dfftp%nnr), thetas(dfftp%nnr, Nqs),               &
            q0(dfftp%nnr), dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
!# 465 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.
!# 471 "xc_vdW_DF.f90"
  total_rho = rho_valence(:,1) + rho_core(:)
!# 474 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.
!# 477 "xc_vdW_DF.f90"
  CALL fft_gradient_r2r (dfftp, total_rho, g, grad_rho)
!# 480 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. This routine also returns the
  ! derivatives of the q0s with respect to the charge-density and the
  ! gradient of the charge-density. These are needed for the potential
  ! calculated below. This routine also calculates the thetas.
!# 488 "xc_vdW_DF.f90"
  CALL get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)
!# 491 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the theta arrays into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.
!# 497 "xc_vdW_DF.f90"
  CALL vdW_DF_energy (thetas, Ec_nl)
  etxc = etxc + Ec_nl
!# 500 "xc_vdW_DF.f90"
  IF ( iverbosity > 0 ) THEN
     CALL mp_sum(Ec_nl, intra_bgrp_comm)
     IF ( ionode ) THEN
        WRITE(stdout,'(/ / A)')       "     -----------------------------------------------"
        WRITE(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        WRITE(stdout,'(A /)')         "     -----------------------------------------------"
     END IF
  END IF
!# 510 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.
!# 518 "xc_vdW_DF.f90"
  DO theta_i = 1, Nqs
     CALL invfft('Rho', thetas(:,theta_i), dfftp)
  END DO
!# 522 "xc_vdW_DF.f90"
  CALL vdW_DF_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, thetas, potential)
  v(:,1) = v(:,1) + e2 * potential(:)
!# 526 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable.
!# 529 "xc_vdW_DF.f90"
  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
!# 531 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * rho_valence(i_grid,1) * potential(i_grid)
  END DO
!# 535 "xc_vdW_DF.f90"
  DEALLOCATE ( total_rho, grad_rho, potential, thetas, q0, dq0_drho, dq0_dgradrho )
!# 537 "xc_vdW_DF.f90"
  END SUBROUTINE xc_vdW_DF
!# 546 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |                  |
  !                          |  XC_VDW_DF_spin  |
  !                          |__________________|
!# 551 "xc_vdW_DF.f90"
  SUBROUTINE xc_vdW_DF_spin (rho_valence, rho_core, etxc, vtxc, v)
!# 553 "xc_vdW_DF.f90"
  !! This subroutine is as similar to `xc_vdW_DF` as possible,
  !! but handles the collinear `nspin=2` case.
!# 557 "xc_vdW_DF.f90"
  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
!# 560 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 562 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Local variables
!# 565 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN) :: rho_valence(:,:)      !
  REAL(DP), INTENT(IN) :: rho_core(:)           ! PWSCF input variables
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:) !
!# 570 "xc_vdW_DF.f90"
  INTEGER :: i_grid, theta_i, i_proc            ! Indexing variables over grid points,
                                                ! theta functions, and processors, and a
                                                ! generic index.
!# 574 "xc_vdW_DF.f90"
  REAL(DP) :: grid_cell_volume                  ! The volume of the unit cell per G-grid point.
!# 576 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: q0(:)                ! The saturated value of q (equations 11 and 12
                                                ! of DION). This saturation is that of
                                                ! equation 5 in SOLER.
!# 580 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: grad_rho(:,:)        ! The gradient of the charge density. The
                                                ! format is as follows:
                                                ! grad_rho(cartesian_component,grid_point).
  REAL(DP), ALLOCATABLE :: grad_rho_up(:,:)     ! The gradient of the up charge density.
                                                ! Same format as grad_rho.
  REAL(DP), ALLOCATABLE :: grad_rho_down(:,:)   ! The gradient of the down charge density.
                                                ! Same format as grad_rho.
!# 588 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: potential_up(:)      ! The vdW contribution to the potential.
  REAL(DP), ALLOCATABLE :: potential_down(:)    ! The vdW contribution to the potential.
!# 591 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: dq0_drho_up(:)       ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_drho_down(:)     ! (equation 5 of SOLER) with respect
                                                ! to the charge density (see
                                                ! get_q0_on_grid subroutine for details).
!# 596 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: dq0_dgradrho_up(:)   ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_dgradrho_down(:) ! (equation 5 of SOLER) with respect
                                                ! to the gradient of the charge density
                                                ! (again, see get_q0_on_grid subroutine).
!# 601 "xc_vdW_DF.f90"
  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)       ! These are the functions of equation 8 of
                                                ! SOLER. They will be forward Fourier transformed
                                                ! in place to get theta(k) and worked on in
                                                ! place to get the u_alpha(r) of equation 11
                                                ! in SOLER. They are formatted as follows:
                                                ! thetas(grid_point, theta_i).
!# 608 "xc_vdW_DF.f90"
  REAL(DP) :: Ec_nl                             ! The non-local vdW contribution to the energy.
!# 610 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: total_rho(:)         ! This is the sum of the valence (up and down)
                                                ! and core charge. This just holds the piece
                                                ! assigned to this processor.
  REAL(DP), ALLOCATABLE :: rho_up(:)            ! This is the just the up valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.
  REAL(DP), ALLOCATABLE :: rho_down(:)          ! This is the just the down valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.
!# 620 "xc_vdW_DF.f90"
  LOGICAL, SAVE :: first_iteration = .TRUE.     ! Whether this is the first time this
                                                ! routine has been called.
!# 626 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Write out the vdW-DF information and initialize the calculation.
!# 629 "xc_vdW_DF.f90"
  IF ( first_iteration ) THEN
     IF ( inlc > 6 ) CALL errore( 'xc_vdW_DF_spin', 'inlc not implemented', 1 )
     CALL generate_kernel
     IF ( ionode ) CALL vdW_info
     first_iteration = .FALSE.
  END IF
!# 637 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.
!# 641 "xc_vdW_DF.f90"
  ALLOCATE( total_rho(dfftp%nnr), rho_up(dfftp%nnr), rho_down(dfftp%nnr),         &
     grad_rho(3,dfftp%nnr), grad_rho_up(3,dfftp%nnr), grad_rho_down(3,dfftp%nnr), &
     potential_up(dfftp%nnr), potential_down(dfftp%nnr), thetas(dfftp%nnr, Nqs),  &
     q0(dfftp%nnr), dq0_drho_up(dfftp%nnr), dq0_dgradrho_up(dfftp%nnr),           &
     dq0_drho_down(dfftp%nnr), dq0_dgradrho_down(dfftp%nnr) )
!# 648 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.
!# 654 "xc_vdW_DF.f90"
  rho_up    = ( rho_valence(:,1) + rho_valence(:,2) + rho_core(:) )*0.5D0
  rho_down  = ( rho_valence(:,1) - rho_valence(:,2) + rho_core(:) )*0.5D0
  total_rho = rho_up + rho_down
!# 659 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.
!# 662 "xc_vdW_DF.f90"
  CALL fft_gradient_r2r (dfftp, total_rho, g, grad_rho)
  CALL fft_gradient_r2r (dfftp, rho_up,    g, grad_rho_up)
  CALL fft_gradient_r2r (dfftp, rho_down,  g, grad_rho_down)
!# 667 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. In the spin case, q0 is defined by
  ! equation 8 (and text above that equation) of THONHAUSER. This
  ! routine also returns the derivatives of the q0s with respect to the
  ! charge-density and the gradient of the charge-density. These are
  ! needed for the potential calculated below.
!# 676 "xc_vdW_DF.f90"
  CALL get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
       grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
       dq0_dgradrho_up, dq0_dgradrho_down, thetas)
!# 681 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.
!# 687 "xc_vdW_DF.f90"
  CALL vdW_DF_energy(thetas, Ec_nl)
  etxc = etxc + Ec_nl
!# 690 "xc_vdW_DF.f90"
  IF ( iverbosity > 0 ) THEN
     CALL mp_sum(Ec_nl, intra_bgrp_comm)
     IF (ionode) THEN
        WRITE(stdout,'(/ / A)')       "     -----------------------------------------------"
        WRITE(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        WRITE(stdout,'(A /)')         "     -----------------------------------------------"
     END IF
  END IF
!# 700 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.
!# 708 "xc_vdW_DF.f90"
  DO theta_i = 1, Nqs
     CALL invfft('Rho', thetas(:,theta_i), dfftp)
  END DO
!# 712 "xc_vdW_DF.f90"
  CALL vdW_DF_potential (q0, dq0_drho_up  , dq0_dgradrho_up  , grad_rho_up  , thetas, potential_up  )
  CALL vdW_DF_potential (q0, dq0_drho_down, dq0_dgradrho_down, grad_rho_down, thetas, potential_down)
!# 715 "xc_vdW_DF.f90"
  v(:,1) = v(:,1) + e2 * potential_up  (:)
  v(:,2) = v(:,2) + e2 * potential_down(:)
!# 719 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable
!# 722 "xc_vdW_DF.f90"
  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
!# 724 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * (rho_valence(i_grid,1) +   &
            rho_valence(i_grid,2)) * 0.5_dp * potential_up  (i_grid)  &
                 + e2 * grid_cell_volume * (rho_valence(i_grid,1) -   &
            rho_valence(i_grid,2)) * 0.5_dp * potential_down(i_grid)
  END DO
!# 731 "xc_vdW_DF.f90"
  DEALLOCATE( total_rho, rho_up, rho_down, grad_rho, grad_rho_up, grad_rho_down, &
              potential_up, potential_down, thetas,                              &
              q0, dq0_drho_up, dq0_dgradrho_up, dq0_drho_down, dq0_dgradrho_down )
!# 735 "xc_vdW_DF.f90"
  END SUBROUTINE xc_vdW_DF_spin
!# 744 "xc_vdW_DF.f90"
  ! ####################################################################
  !                       |                  |
  !                       |  GET_Q0_ON_GRID  |
  !                       |__________________|
!# 749 "xc_vdW_DF.f90"
  SUBROUTINE get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)
!# 751 "xc_vdW_DF.f90"
  !! This routine first calculates the `q` value defined in (DION equations
  !! 11 and 12), then saturates it according to (SOLER equation 5). More
  !! specifically it calculates the following:
  !!
  !! * `q0(ir)` = saturated value of q
  !! * `dq0_drho(ir) = total_rho * d q0/d rho`
  !! * `dq0_dgradrho = total_rho/|grad_rho| * d q0/d |grad_rho|`
!# 760 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 762 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)      :: total_rho(:), grad_rho(:,:)         ! Input variables needed
!# 764 "xc_vdW_DF.f90"
  REAL(DP), INTENT(OUT)     :: q0(:), dq0_drho(:), dq0_dgradrho(:) ! Output variables that have been allocated
                                                                   ! outside this routine but will be set here.
  COMPLEX(DP), INTENT(INOUT):: thetas(:,:)                         ! The thetas from SOLER.
!# 768 "xc_vdW_DF.f90"
  INTEGER, PARAMETER        :: m_cut = 12                          ! How many terms to include in the sum
                                                                   ! of SOLER equation 5.
!# 771 "xc_vdW_DF.f90"
  REAL(DP)                  :: rho                                 ! Local variable for the density
  REAL(DP)                  :: r_s                                 ! Wigner-Seitz radius
  REAL(DP)                  :: s                                   ! Reduced gradient
  REAL(DP)                  :: q
  REAL(DP)                  :: ec
  REAL(DP)                  :: dq0_dq                              ! The derivative of the saturated
                                                                   ! q0 with respect to q.
!# 779 "xc_vdW_DF.f90"
  INTEGER                   :: i_grid                              ! Indexing variable
!# 784 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays
!# 787 "xc_vdW_DF.f90"
  q0(:)           = q_cut
  dq0_drho(:)     = 0.0D0
  dq0_dgradrho(:) = 0.0D0
!# 792 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
!# 794 "xc_vdW_DF.f90"
     rho = total_rho(i_grid)
!# 797 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! This prevents numerical problems. If the charge density is
     ! negative (an unphysical situation), we simply treat it as very
     ! small. In that case, q0 will be very large and will be saturated.
     ! For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
     ! q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
     ! point.
!# 805 "xc_vdW_DF.f90"
     IF ( rho < epsr ) CYCLE
!# 808 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Calculate some intermediate values needed to find q.
!# 811 "xc_vdW_DF.f90"
     r_s = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)
!# 813 "xc_vdW_DF.f90"
     s   = SQRT( grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2 ) / &
           (2.0D0 * kF(rho) * rho )
!# 817 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION.
     ! Use pw() from XClib/qe_funct_corr_lda_lsda.f90 to get
     ! qc = kf/eps_x * eps_c.
!# 822 "xc_vdW_DF.f90"
     CALL pw(r_s, 1, ec, dq0_drho(i_grid))
     q = -4.0D0*pi/3.0D0 * ec + kF(rho) * Fs(s)
!# 826 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.
!# 829 "xc_vdW_DF.f90"
     CALL saturate_q ( q, q_cut, q0(i_grid), dq0_dq )
     IF (q0(i_grid) < q_min) q0(i_grid) = q_min
!# 833 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Here we find derivatives. These are actually the density times
     ! the derivative of q0 with respect to rho and grad_rho. The
     ! density factor comes in since we are really differentiating
     ! theta = (rho)*P(q0) with respect to density (or its gradient)
     ! which will be
     !
     !    dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]
     !
     ! and
     !
     !    dtheta_dgrad_rho = dP_dq0 * [rho * dq0_dq * dq_dgrad_rho]
     !
     ! The parts in square brackets are what is calculated here. The
     ! dP_dq0 term will be interpolated later.
!# 849 "xc_vdW_DF.f90"
     dq0_drho(i_grid)     = dq0_dq * rho * ( -4.0D0*pi/3.0D0 * &
                            (dq0_drho(i_grid) - ec)/rho + dqx_drho(rho, s) )
     dq0_dgradrho(i_grid) = dq0_dq * rho * kF(rho) * dFs_ds(s) * ds_dgradrho(rho)
!# 853 "xc_vdW_DF.f90"
  END DO
!# 856 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the theta functions of SOLER equation 8. These are
  ! defined as
  !
  !    rho * P_i(q0(rho, grad_rho))
  !
  ! where P_i is a polynomial that interpolates a Kroneker delta
  ! function at the point q_i (taken from the q_mesh) and q0 is the
  ! saturated version of q. q is defined in equations 11 and 12 of DION
  ! and the saturation proceedure is defined in equation 5 of SOLER.
  ! This is the biggest memory consumer in the method since the thetas
  ! array is (total # of FFT points)*Nqs complex numbers. In a parallel
  ! run, each processor will hold the values of all the theta functions
  ! on just the points assigned to it. thetas are stored in reciprocal
  ! space as theta_i(k) because this is the way they are used later for
  ! the convolution (equation 8 of SOLER). Start by interpolating the
  ! P_i polynomials defined in equation 3 in SOLER for the particular q0
  ! values we have.
!# 875 "xc_vdW_DF.f90"
  CALL spline_interpolation (q_mesh, q0, thetas)
!# 877 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  END DO
!# 881 "xc_vdW_DF.f90"
  DO idx = 1, Nqs
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  END DO
!# 885 "xc_vdW_DF.f90"
  END SUBROUTINE get_q0_on_grid
!# 894 "xc_vdW_DF.f90"
  ! ####################################################################
  !                       |                       |
  !                       |  GET_Q0_ON_GRID_spin  |
  !                       |_______________________|
!# 899 "xc_vdW_DF.f90"
  SUBROUTINE get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
             grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
             dq0_dgradrho_up, dq0_dgradrho_down, thetas)
!# 903 "xc_vdW_DF.f90"
  !! Find the value of `q0` for all assigned grid points. `q` is defined in
  !! equations 11 and 12 of DION and `q0` is the saturated version of `q`
  !! defined in equation 5 of SOLER. In the spin case, `q0` is defined by
  !! equation 8 (and text above that equation) of THONHAUSER. This
  !! routine also returns the derivatives of the `q0`s with respect to the
  !! charge-density and the gradient of the charge-density. These are
  !! needed for the potential.
!# 912 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 914 "xc_vdW_DF.f90"
  REAL(DP),  INTENT(IN)      :: total_rho(:), grad_rho(:,:)              ! Input variables
  REAL(DP),  INTENT(IN)      :: rho_up(:), grad_rho_up(:,:)              ! Input variables
  REAL(DP),  INTENT(IN)      :: rho_down(:), grad_rho_down(:,:)          ! Input variables
!# 918 "xc_vdW_DF.f90"
  REAL(DP),  INTENT(OUT)     :: q0(:), dq0_drho_up(:), dq0_drho_down(:)  ! Output variables
  REAL(DP),  INTENT(OUT)     :: dq0_dgradrho_up(:), dq0_dgradrho_down(:) ! Output variables
  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)                              ! Thetas from SOLER
!# 922 "xc_vdW_DF.f90"
  REAL(DP)                   :: rho, up, down                            ! Local copy of densities
  REAL(DP)                   :: zeta                                     ! Spin polarization
  REAL(DP)                   :: r_s                                      ! Wigner-Seitz radius
  REAL(DP)                   :: q, qc, qx, qx_up, qx_down                ! q for exchange and correlation
  REAL(DP)                   :: q0x_up, q0x_down                         ! Saturated q values
  REAL(DP)                   :: fac
  REAL(DP)                   :: ec, vc(2)
  REAL(DP)                   :: dq0_dq, dq0x_up_dq, dq0x_down_dq         ! Derivative of q0 w.r.t q
  REAL(DP)                   :: dqc_drho_up, dqc_drho_down               ! Intermediate values
  REAL(DP)                   :: dqx_drho_up, dqx_drho_down               ! Intermediate values
  REAL(DP)                   :: s_up, s_down                             ! Reduced gradients
  INTEGER                    :: i_grid                                   ! Indexing variable
  LOGICAL                    :: calc_qx_up, calc_qx_down
!# 939 "xc_vdW_DF.f90"
  fac = 2.0D0**(-1.0D0/3.0D0)
!# 942 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays
!# 945 "xc_vdW_DF.f90"
  q0(:)                = q_cut
  dq0_drho_up(:)       = 0.0D0
  dq0_drho_down(:)     = 0.0D0
  dq0_dgradrho_up(:)   = 0.0D0
  dq0_dgradrho_down(:) = 0.0D0
!# 952 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
!# 954 "xc_vdW_DF.f90"
     rho  = total_rho(i_grid)
     up   = rho_up(i_grid)
     down = rho_down(i_grid)
!# 959 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! This prevents numerical problems. If the charge density is
     ! negative (an unphysical situation), we simply treat it as very
     ! small. In that case, q0 will be very large and will be saturated.
     ! For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
     ! q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
     ! point.
!# 967 "xc_vdW_DF.f90"
     IF ( rho < epsr ) CYCLE
!# 969 "xc_vdW_DF.f90"
     calc_qx_up   = .TRUE.
     calc_qx_down = .TRUE.
!# 972 "xc_vdW_DF.f90"
     IF ( up   < epsr/2.0D0 ) calc_qx_up   = .FALSE.
     IF ( down < epsr/2.0D0 ) calc_qx_down = .FALSE.
!# 976 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! The spin case is numerically even more tricky and we have to
     ! saturate each spin channel separately. Note that we are
     ! saturating at a higher value here, so that very large q values
     ! get saturated to exactly q_cut in the second, overall saturation.
!# 982 "xc_vdW_DF.f90"
     q0x_up        = 0.0D0
     q0x_down      = 0.0D0
     dqx_drho_up   = 0.0D0
     dqx_drho_down = 0.0D0
!# 988 "xc_vdW_DF.f90"
     IF ( calc_qx_up ) THEN
        s_up    = SQRT( grad_rho_up(1,i_grid)**2 + grad_rho_up(2,i_grid)**2 + &
                  grad_rho_up(3,i_grid)**2 ) / (2.0D0 * kF(up) * up)
        qx_up   = kF(2.0D0*up) * Fs(fac*s_up)
        CALL saturate_q (qx_up, 4.0D0*q_cut, q0x_up, dq0x_up_dq)
     END IF
!# 995 "xc_vdW_DF.f90"
     IF ( calc_qx_down ) THEN
        s_down  = SQRT( grad_rho_down(1,i_grid)**2 + grad_rho_down(2,i_grid)**2 + &
                  grad_rho_down(3,i_grid)**2) / (2.0D0 * kF(down) * down)
        qx_down = kF(2.0D0*down) * Fs(fac*s_down)
        CALL saturate_q (qx_down, 4.0D0*q_cut, q0x_down, dq0x_down_dq)
     END IF
!# 1003 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION and
     ! equation 8 of THONHAUSER (also see text above that equation).
!# 1007 "xc_vdW_DF.f90"
     r_s  = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)
     zeta = (up - down) / rho
     IF ( ABS(zeta) > 1.0D0 ) zeta = SIGN(1.0D0, zeta)
     call pw_spin( r_s, zeta, ec, vc(1), vc(2) )
     dqc_drho_up   = vc(1)
     dqc_drho_down = vc(2)
!# 1014 "xc_vdW_DF.f90"
     qx = ( up * q0x_up + down * q0x_down ) / rho
     qc = -4.0D0*pi/3.0D0 * ec
     q  = qx + qc
!# 1019 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.
!# 1022 "xc_vdW_DF.f90"
     CALL saturate_q (q, q_cut, q0(i_grid), dq0_dq)
     IF (q0(i_grid) < q_min) q0(i_grid) = q_min
!# 1026 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Here we find derivatives. These are actually the density times
     ! the derivative of q0 with respect to rho and grad_rho. The
     ! density factor comes in since we are really differentiating
     ! theta = (rho)*P(q0) with respect to density (or its gradient)
     ! which will be
     !
     !    dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]
     !
     ! and
     !
     !    dtheta_dgrad_rho = dP_dq0 * [rho * dq0_dq * dq_dgrad_rho]
     !
     ! The parts in square brackets are what is calculated here. The
     ! dP_dq0 term will be interpolated later.
!# 1042 "xc_vdW_DF.f90"
     IF ( calc_qx_up ) THEN
        dqx_drho_up   = 2.0D0*dq0x_up_dq*up*dqx_drho(2.0D0*up, fac*s_up) + q0x_up*down/rho
        dq0_dgradrho_up (i_grid) = 2.0D0 * dq0_dq * dq0x_up_dq * up * kF(2.0D0*up) * &
                        dFs_ds(fac*s_up) * ds_dgradrho(2.0D0*up)
     END IF
!# 1048 "xc_vdW_DF.f90"
     IF ( calc_qx_down ) THEN
        dqx_drho_down = 2.0D0*dq0x_down_dq*down*dqx_drho(2.0D0*down, fac*s_down) + q0x_down*up/rho
        dq0_dgradrho_down(i_grid) = 2.0D0 * dq0_dq * dq0x_down_dq * down * kF(2.0D0*down) * &
                        dFs_ds(fac*s_down) * ds_dgradrho(2.0D0*down)
     END IF
!# 1054 "xc_vdW_DF.f90"
     IF ( calc_qx_down ) dqx_drho_up   = dqx_drho_up   - q0x_down*down/rho
     IF ( calc_qx_up )   dqx_drho_down = dqx_drho_down - q0x_up  *up  /rho
!# 1057 "xc_vdW_DF.f90"
     dqc_drho_up   = -4.0D0*pi/3.0D0 * (dqc_drho_up   - ec)
     dqc_drho_down = -4.0D0*pi/3.0D0 * (dqc_drho_down - ec)
!# 1060 "xc_vdW_DF.f90"
     dq0_drho_up  (i_grid) = dq0_dq * (dqc_drho_up   + dqx_drho_up  )
     dq0_drho_down(i_grid) = dq0_dq * (dqc_drho_down + dqx_drho_down)
!# 1063 "xc_vdW_DF.f90"
  END DO
!# 1066 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the theta functions of SOLER equation 8. These are
  ! defined as
  !
  !    rho * P_i(q0(rho, grad_rho))
  !
  ! where P_i is a polynomial that interpolates a Kroneker delta
  ! function at the point q_i (taken from the q_mesh) and q0 is the
  ! saturated version of q. q is defined in equations 11 and 12 of DION
  ! and the saturation proceedure is defined in equation 5 of SOLER.
  ! This is the biggest memory consumer in the method since the thetas
  ! array is (total # of FFT points)*Nqs complex numbers. In a parallel
  ! run, each processor will hold the values of all the theta functions
  ! on just the points assigned to it. thetas are stored in reciprocal
  ! space as theta_i(k) because this is the way they are used later for
  ! the convolution (equation 8 of SOLER). Start by interpolating the
  ! P_i polynomials defined in equation 3 in SOLER for the particular q0
  ! values we have.
!# 1085 "xc_vdW_DF.f90"
  CALL spline_interpolation (q_mesh, q0, thetas)
!# 1087 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  END DO
!# 1091 "xc_vdW_DF.f90"
  DO idx = 1, Nqs
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  END DO
!# 1095 "xc_vdW_DF.f90"
  END SUBROUTINE get_q0_on_grid_spin
!# 1104 "xc_vdW_DF.f90"
  ! ####################################################################
  !                            |              |
  !                            |  saturate_q  |
  !                            |______________|
!# 1109 "xc_vdW_DF.f90"
  SUBROUTINE saturate_q (q, q_cutoff, q0, dq0_dq)
!# 1111 "xc_vdW_DF.f90"
  !! Here, we calculate `q0` by saturating `q` according to equation 5
  !! of SOLER. Also, we find the derivative `dq0/dq` needed for the
  !! derivatives `dq0/drho` and `dq0/dgradrho`.
!# 1116 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1118 "xc_vdW_DF.f90"
  REAL(DP),  INTENT(IN)      :: q             ! Input q
  REAL(DP),  INTENT(IN)      :: q_cutoff      ! Cutoff q
  REAL(DP),  INTENT(OUT)     :: q0            ! Output saturated q
  REAL(DP),  INTENT(OUT)     :: dq0_dq        ! Derivative of dq0/dq
!# 1123 "xc_vdW_DF.f90"
  REAL(DP)                   :: e_exp         ! Exponent
  INTEGER,   PARAMETER       :: m_cut = 12    ! How many terms to include in
                                              ! the sum of SOLER equation 5.
!# 1128 "xc_vdW_DF.f90"
  e_exp  = 0.0D0
  dq0_dq = 0.0D0
!# 1131 "xc_vdW_DF.f90"
  DO idx = 1, m_cut
     e_exp  = e_exp + (q/q_cutoff)**idx/idx
     dq0_dq = dq0_dq + (q/q_cutoff)**(idx-1)
  END Do
!# 1136 "xc_vdW_DF.f90"
  q0     = q_cutoff*(1.0D0 - EXP(-e_exp))
  dq0_dq = dq0_dq * EXP(-e_exp)
!# 1139 "xc_vdW_DF.f90"
  END SUBROUTINE saturate_q
!# 1148 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |               |
  !                          | vdW_DF_energy |
  !                          |_______________|
!# 1153 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_energy (thetas, vdW_xc_energy)
!# 1155 "xc_vdW_DF.f90"
  !! This routine carries out the integration of equation 8 of SOLER. It
  !! returns the non-local exchange-correlation energy and the
  !! `u_alpha(k)` arrays used to find the `u_alpha(r)` arrays via
  !! equations 11 and 12 in SOLER.
!# 1161 "xc_vdW_DF.f90"
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega
!# 1164 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1166 "xc_vdW_DF.f90"
  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)            ! On input this variable holds the theta
                                                       ! functions (equation 8, SOLER) in the format
                                                       ! thetas(grid_point, theta_i). On output
                                                       ! this array holds u_alpha(k) =
                                                       ! Sum_j[theta_beta(k)phi_alpha_beta(k)].
!# 1172 "xc_vdW_DF.f90"
  REAL(DP), INTENT(OUT) :: vdW_xc_energy               ! The non-local correlation energy
!# 1174 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: kernel_of_k(:,:)            ! This array will hold the interpolated kernel
                                                       ! values for each pair of q values in the q_mesh.
!# 1177 "xc_vdW_DF.f90"
  REAL(DP)    :: g                                     ! The magnitude of the current g vector
  INTEGER     :: last_g                                ! The shell number of the last g vector
!# 1181 "xc_vdW_DF.f90"
  INTEGER     :: g_i, q1_i, q2_i, i_grid               ! Index variables
!# 1183 "xc_vdW_DF.f90"
  COMPLEX(DP) :: theta(Nqs), thetam(Nqs), theta_g(Nqs) ! Temporary storage arrays used since we
                                                       ! are overwriting the thetas array here.
  REAL(DP)    :: G0_term, G_multiplier
!# 1187 "xc_vdW_DF.f90"
  COMPLEX(DP), ALLOCATABLE :: u_vdw(:,:)               ! Temporary array holding u_alpha(k)
!# 1192 "xc_vdW_DF.f90"
  ALLOCATE ( u_vdW(dfftp%nnr,Nqs), kernel_of_k(Nqs, Nqs) )
  vdW_xc_energy = 0.0D0
  u_vdW(:,:)    = CMPLX(0.0_DP, 0.0_DP, kind=dp)
!# 1197 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Loop over PWSCF's array of magnitude-sorted g-vector shells. For
  ! each shell, interpolate the kernel at this magnitude of g, then find
  ! all points on the shell and carry out the integration over those
  ! points. The PWSCF variables used here are ngm = number of g-vectors
  ! on this processor, nl = an array that gives the indices into the FFT
  ! grid for a particular g vector, igtongl = an array that gives the
  ! index of which shell a particular g vector is in, gl = an array that
  ! gives the magnitude of the g vectors for each shell. In essence, we
  ! are forming the reciprocal-space u(k) functions of SOLER equation
  ! 11. These are kept in thetas array. Here we should use gstart,ngm
  ! but all the cases are handled by conditionals inside the loop
!# 1210 "xc_vdW_DF.f90"
  G_multiplier = 1.0D0
  IF ( gamma_only ) G_multiplier = 2.0D0
!# 1213 "xc_vdW_DF.f90"
  last_g = -1
!# 1215 "xc_vdW_DF.f90"
  DO g_i = 1, ngm
!# 1217 "xc_vdW_DF.f90"
     IF ( igtongl(g_i) .NE. last_g) THEN
        g = SQRT(gl(igtongl(g_i))) * tpiba
        CALL interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
     END IF
!# 1223 "xc_vdW_DF.f90"
     theta = thetas(dfftp%nl(g_i),:)
!# 1225 "xc_vdW_DF.f90"
     DO q2_i = 1, Nqs
        DO q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i)  = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        END DO
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(dfftp%nl(g_i),q2_i)*CONJG(theta(q2_i)))
     END DO
!# 1232 "xc_vdW_DF.f90"
     IF (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier
!# 1234 "xc_vdW_DF.f90"
  END DO
!# 1236 "xc_vdW_DF.f90"
  IF ( gamma_only ) u_vdW(dfftp%nlm(:),:) = CONJG( u_vdW(dfftp%nl(:),:) )
!# 1239 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Apply scaling factors. The e2 comes from PWSCF's choice of units.
  ! This should be 0.5 * e2 * vdW_xc_energy * (2pi)^3/omega * (omega)^2,
  ! with the (2pi)^3/omega being the volume element for the integral
  ! (the volume of the reciprocal unit cell) and the 2 factors of omega
  ! being used to cancel the factor of 1/omega PWSCF puts on forward
  ! FFTs of the 2 theta factors. 1 omega cancels and the (2pi)^3
  ! cancels because there should be a factor of 1/(2pi)^3 on the radial
  ! Fourier transform of phi that was left out to cancel with this
  ! factor.
!# 1250 "xc_vdW_DF.f90"
  vdW_xc_energy = 0.5D0 * e2 * omega * vdW_xc_energy
!# 1252 "xc_vdW_DF.f90"
  thetas(:,:) = u_vdW(:,:)
  DEALLOCATE ( u_vdW, kernel_of_k )
!# 1255 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_energy
!# 1264 "xc_vdW_DF.f90"
  ! ####################################################################
  !                        |                   |
  !                        |  vdW_DF_potential |
  !                        |___________________|
!# 1269 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, u_vdW, potential)
!# 1271 "xc_vdW_DF.f90"
  !! This routine finds the non-local correlation contribution to the
  !! potential (i.e. the derivative of the non-local piece of the energy
  !! with respect to density) given in SOLER equation 10. The
  !! `u_alpha(k)` functions were found while calculating the energy.
  !! They are passed in as the matrix `u_vdW`. Most of the required
  !! derivatives were calculated in the `get_q0_on_grid` routine, but
  !! the derivative of the interpolation polynomials, `P_alpha(q)`
  !! (SOLER equation 3), with respect to `q` is interpolated here, along
  !! with the polynomials themselves.
!# 1282 "xc_vdW_DF.f90"
  USE gvect,               ONLY : g
  USE cell_base,           ONLY : alat, tpiba
!# 1285 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1287 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN) ::  q0(:), grad_rho(:,:)       ! Input arrays holding the value of q0 for
                                                      ! all points assigned to this processor and
                                                      ! the gradient of the charge density for
                                                      ! points assigned to this processor.
!# 1292 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN) :: dq0_drho(:), dq0_dgradrho(:)! The derivative of q0 with respect to the
                                                      ! charge density and gradient of the charge
                                                      ! density (almost). See comments in the
                                                      ! get_q0_on_grid subroutine above.
!# 1297 "xc_vdW_DF.f90"
  COMPLEX(DP), INTENT(IN) :: u_vdW(:,:)               ! The functions u_alpha(r) obtained by
                                                      ! inverse transforming the functions
                                                      ! u_alph(k). See equations 11 and 12 in SOLER.
!# 1301 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: potential(:)             ! The non-local correlation potential for
                                                      ! points on the grid over the whole cell (not
                                                      ! just those assigned to this processor).
!# 1305 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE, SAVE :: d2y_dx2(:,:)         ! Second derivatives of P_alpha polynomials
                                                      ! for interpolation.
!# 1308 "xc_vdW_DF.f90"
  INTEGER :: i_grid, P_i,icar                         ! Index variables
!# 1310 "xc_vdW_DF.f90"
  INTEGER :: q_low, q_hi, q                           ! Variables to find the bin in the q_mesh that
                                                      ! a particular q0 belongs to (for interpolation).
  REAL(DP) :: dq, a, b, c, d, e, f                    ! Intermediate variables used in the
                                                      ! interpolation of the polynomials.
!# 1315 "xc_vdW_DF.f90"
  REAL(DP) :: y(Nqs), dP_dq0, P                       ! The y values for a given polynomial (all 0
                                                      ! exept for element i of P_i) The derivative
                                                      ! of P at a given q0 and the value of P at a
                                                      ! given q0. Both of these are interpolated
                                                      ! below.
!# 1321 "xc_vdW_DF.f90"
  REAL(DP) :: gradient2                               ! Squared gradient
!# 1323 "xc_vdW_DF.f90"
  REAL(DP)   , ALLOCATABLE ::h_prefactor(:)
  COMPLEX(DP), ALLOCATABLE ::h(:)
!# 1329 "xc_vdW_DF.f90"
  ALLOCATE ( h_prefactor(dfftp%nnr), h(dfftp%nnr) )
!# 1331 "xc_vdW_DF.f90"
  potential     = 0.0D0
  h_prefactor   = 0.0D0
!# 1335 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get the second derivatives of the P_i functions for interpolation.
  ! We have already calculated this once but it is very fast and it's
  ! just as easy to calculate it again.
!# 1340 "xc_vdW_DF.f90"
  IF (.NOT. ALLOCATED( d2y_dx2) ) THEN
!# 1342 "xc_vdW_DF.f90"
     ALLOCATE( d2y_dx2(Nqs, Nqs) )
     CALL initialize_spline_interpolation ( q_mesh, d2y_dx2(:,:) )
!# 1345 "xc_vdW_DF.f90"
  end if
!# 1348 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
!# 1350 "xc_vdW_DF.f90"
     q_low = 1
     q_hi = Nqs
!# 1354 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in in the q_mesh.
!# 1357 "xc_vdW_DF.f90"
     DO WHILE ( (q_hi - q_low) > 1)
!# 1359 "xc_vdW_DF.f90"
        q = INT((q_hi + q_low)/2)
!# 1361 "xc_vdW_DF.f90"
        IF (q_mesh(q) > q0(i_grid)) THEN
           q_hi = q
        ELSE
           q_low = q
        END IF
!# 1367 "xc_vdW_DF.f90"
     END DO
!# 1369 "xc_vdW_DF.f90"
     IF ( q_hi == q_low ) CALL errore('vdW_DF_potential','qhi == qlow',1)
!# 1371 "xc_vdW_DF.f90"
     dq = q_mesh(q_hi) - q_mesh(q_low)
!# 1373 "xc_vdW_DF.f90"
     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
!# 1380 "xc_vdW_DF.f90"
     DO P_i = 1, Nqs
        y = 0.0D0
        y(P_i) = 1.0D0
!# 1384 "xc_vdW_DF.f90"
        P      = a*y(q_low) + b*y(q_hi)  + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
!# 1388 "xc_vdW_DF.f90"
        ! --------------------------------------------------------------
        ! The first term in equation 10 of SOLER.
!# 1391 "xc_vdW_DF.f90"
        potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))
        IF (q0(i_grid) .NE. q_mesh(Nqs)) THEN
           h_prefactor(i_grid) = h_prefactor(i_grid) + u_vdW(i_grid,P_i)*dP_dq0*dq0_dgradrho(i_grid)
        END IF
!# 1396 "xc_vdW_DF.f90"
     END DO
!# 1398 "xc_vdW_DF.f90"
  END DO
!# 1400 "xc_vdW_DF.f90"
  DO icar = 1,3
!# 1402 "xc_vdW_DF.f90"
     h(:) = CMPLX( h_prefactor(:) * grad_rho(icar,:), 0.0_DP, kind=dp )
!# 1404 "xc_vdW_DF.f90"
     DO i_grid = 1, dfftp%nnr
        gradient2 = grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2
        IF ( gradient2 > 0.0D0 ) h(i_grid) = h(i_grid) / SQRT( gradient2 )
     END DO
!# 1409 "xc_vdW_DF.f90"
     CALL fwfft ('Rho', h, dfftp)
     h(dfftp%nl(:)) = CMPLX(0.0_DP,1.0_DP, kind=dp) * tpiba * g(icar,:) * h(dfftp%nl(:))
     IF (gamma_only) h(dfftp%nlm(:)) = CONJG(h(dfftp%nl(:)))
     CALL invfft ('Rho', h, dfftp)
     potential(:) = potential(:) - REAL(h(:))
!# 1415 "xc_vdW_DF.f90"
  END Do
!# 1417 "xc_vdW_DF.f90"
  DEALLOCATE ( h_prefactor, h )
!# 1419 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_potential
!# 1428 "xc_vdW_DF.f90"
  ! ####################################################################
  !                       |                        |
  !                       |  SPLINE_INTERPOLATION  |
  !                       |________________________|
!# 1433 "xc_vdW_DF.f90"
  SUBROUTINE spline_interpolation (x, evaluation_points, values)
!# 1435 "xc_vdW_DF.f90"
  !! This routine is modeled after an algorithm from NUMERICAL_RECIPES.
  !! It was adapted for Fortran, of course, and for the problem at hand,
  !! in that it finds the bin a particular `x` value is in and then loops
  !! over all the `P_i` functions so we only have to find the bin once.
!# 1441 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1443 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN) :: x(:), evaluation_points(:)     ! Input variables. The x values used to
                                                         ! form the interpolation (q_mesh in this
                                                         ! case) and the values of q0 for which we
                                                         ! are interpolating the function.
!# 1448 "xc_vdW_DF.f90"
  COMPLEX(DP), INTENT(INOUT) :: values(:,:)              ! An output array (allocated outside this
                                                         ! routine) that stores the interpolated
                                                         ! values of the P_i (SOLER equation 3)
                                                         ! polynomials. The format is
                                                         ! values(grid_point, P_i).
!# 1454 "xc_vdW_DF.f90"
  INTEGER :: Ngrid_points, Nx                            ! Total number of grid points to evaluate
                                                         ! and input x points.
!# 1457 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE, SAVE :: d2y_dx2(:,:)            ! The second derivatives required to do
                                                         ! the interpolation.
!# 1460 "xc_vdW_DF.f90"
  INTEGER :: i_grid, lower_bound, upper_bound, P_i       ! Some indexing variables.
!# 1462 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: y(:)                          ! Temporary variables needed for the
  REAL(DP) :: a, b, c, d, dx                             ! interpolation.
!# 1468 "xc_vdW_DF.f90"
  Nx = size(x)
  Ngrid_points = size(evaluation_points)
!# 1472 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Allocate the temporary array
!# 1475 "xc_vdW_DF.f90"
  ALLOCATE( y(Nx) )
!# 1478 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! If this is the first time this routine has been called we need to
  ! get the second derivatives (d2y_dx2) required to perform the
  ! interpolations. So we allocate the array and call
  ! initialize_spline_interpolation to get d2y_dx2.
!# 1484 "xc_vdW_DF.f90"
  IF (.NOT. ALLOCATED(d2y_dx2) ) THEN
!# 1486 "xc_vdW_DF.f90"
     ALLOCATE( d2y_dx2(Nx,Nx) )
     CALL initialize_spline_interpolation(x, d2y_dx2)
!# 1489 "xc_vdW_DF.f90"
  END IF
!# 1491 "xc_vdW_DF.f90"
  DO i_grid=1, Ngrid_points
!# 1493 "xc_vdW_DF.f90"
     lower_bound = 1
     upper_bound = Nx
!# 1496 "xc_vdW_DF.f90"
     DO WHILE ( (upper_bound - lower_bound) > 1 )
!# 1498 "xc_vdW_DF.f90"
        idx = (upper_bound+lower_bound) / 2
!# 1500 "xc_vdW_DF.f90"
        IF ( evaluation_points(i_grid) > x(idx) ) THEN
           lower_bound = idx
        ELSE
           upper_bound = idx
        END IF
!# 1506 "xc_vdW_DF.f90"
     END DO
!# 1508 "xc_vdW_DF.f90"
     dx = x(upper_bound)-x(lower_bound)
!# 1510 "xc_vdW_DF.f90"
     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0
!# 1515 "xc_vdW_DF.f90"
     DO P_i = 1, Nx
!# 1517 "xc_vdW_DF.f90"
        y = 0
        y(P_i) = 1
!# 1520 "xc_vdW_DF.f90"
        values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
     END DO
!# 1524 "xc_vdW_DF.f90"
  END DO
!# 1526 "xc_vdW_DF.f90"
  DEALLOCATE( y )
!# 1528 "xc_vdW_DF.f90"
  END SUBROUTINE spline_interpolation
!# 1537 "xc_vdW_DF.f90"
  ! ####################################################################
  !                  |                                   |
  !                  |  INITIALIZE_SPLINE_INTERPOLATION  |
  !                  |___________________________________|
!# 1542 "xc_vdW_DF.f90"
  SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)
!# 1544 "xc_vdW_DF.f90"
  !! This routine is modeled after an algorithm from NUMERICAL_RECIPES.
  !! It was adapted for Fortran and for the problem at hand.
!# 1548 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1550 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: x(:)               ! The input abscissa values.
  REAL(DP), INTENT(INOUT) :: d2y_dx2(:,:)       ! The output array (allocated outside this routine)
                                                ! that holds the second derivatives required for
                                                ! interpolating the function.
!# 1555 "xc_vdW_DF.f90"
  INTEGER :: Nx, P_i                            ! The total number of x points and some indexing
                                                ! variables.
!# 1558 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: temp_array(:), y(:)  ! Some temporary arrays required. y is the array
                                                ! that holds the funcion values (all either 0 or
                                                ! 1 here).
!# 1562 "xc_vdW_DF.f90"
  REAL(DP) :: temp1, temp2                      ! Some temporary variables required.
!# 1567 "xc_vdW_DF.f90"
  Nx = SIZE(x)
!# 1569 "xc_vdW_DF.f90"
  ALLOCATE( temp_array(Nx), y(Nx) )
!# 1571 "xc_vdW_DF.f90"
  DO P_i=1, Nx
!# 1573 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! In the Soler method, the polynomicals that are interpolated are Kroneker
     ! delta funcions at a particular q point. So, we set all y values to 0
     ! except the one corresponding to the particular function P_i.
!# 1578 "xc_vdW_DF.f90"
     y = 0.0D0
     y(P_i) = 1.0D0
!# 1581 "xc_vdW_DF.f90"
     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0
!# 1584 "xc_vdW_DF.f90"
     DO idx = 2, Nx-1
!# 1586 "xc_vdW_DF.f90"
        temp1 = (x(idx)-x(idx-1))/(x(idx+1)-x(idx-1))
        temp2 = temp1 * d2y_dx2(P_i,idx-1) + 2.0D0
        d2y_dx2(P_i,idx) = (temp1-1.0D0)/temp2
!# 1590 "xc_vdW_DF.f90"
        temp_array(idx) = (y(idx+1)-y(idx))/(x(idx+1)-x(idx)) &
             - (y(idx)-y(idx-1))/(x(idx)-x(idx-1))
        temp_array(idx) = (6.0D0*temp_array(idx)/(x(idx+1)-x(idx-1)) &
             - temp1*temp_array(idx-1))/temp2
!# 1595 "xc_vdW_DF.f90"
     END DO
!# 1597 "xc_vdW_DF.f90"
     d2y_dx2(P_i,Nx) = 0.0D0
!# 1599 "xc_vdW_DF.f90"
     DO idx=Nx-1, 1, -1
!# 1601 "xc_vdW_DF.f90"
        d2y_dx2(P_i,idx) = d2y_dx2(P_i,idx) * d2y_dx2(P_i,idx+1) + temp_array(idx)
!# 1603 "xc_vdW_DF.f90"
     END DO
!# 1605 "xc_vdW_DF.f90"
  END DO
!# 1607 "xc_vdW_DF.f90"
  DEALLOCATE( temp_array, y )
!# 1609 "xc_vdW_DF.f90"
  END SUBROUTINE initialize_spline_interpolation
!# 1618 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |                    |
  !                          | INTERPOLATE_KERNEL |
  !                          |____________________|
!# 1623 "xc_vdW_DF.f90"
  SUBROUTINE interpolate_kernel (k, kernel_of_k)
!# 1625 "xc_vdW_DF.f90"
  !! This routine is modeled after an algorithm from NUMERICAL_RECIPES,
  !! adapted for Fortran and the problem at hand. This function is used
  !! to find the Phi-alpha-beta needed for equations 8 and 11 of SOLER.
!# 1630 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1632 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: k                 ! Input value, the magnitude of the g-vector
                                               ! for the current point.
!# 1635 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: kernel_of_k(:,:)  ! An output array (allocated outside this routine)
                                               ! that holds the interpolated value of the kernel
                                               ! for each pair of q points (i.e. the phi_alpha_beta
                                               ! of the Soler method.
!# 1640 "xc_vdW_DF.f90"
  INTEGER  :: q1_i, q2_i, k_i                  ! Indexing variables
!# 1642 "xc_vdW_DF.f90"
  REAL(DP) :: A, B, C, D                       ! Intermediate values for the interpolation.
!# 1647 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Check to make sure that the kernel table we have is capable of
  ! dealing with this value of k. If k is larger than
  ! Nr_points*2*pi/r_max then we can't perform the interpolation. In
  ! that case, a kernel file should be generated with a larger number of
  ! radial points.
!# 1654 "xc_vdW_DF.f90"
  IF ( k >= Nr_points*dk ) THEN
!# 1656 "xc_vdW_DF.f90"
     WRITE(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ", Nr_points*dk
     CALL errore('interpolate kernel', 'k value requested is out of range',1)
!# 1659 "xc_vdW_DF.f90"
  END IF
!# 1661 "xc_vdW_DF.f90"
  kernel_of_k = 0.0D0
!# 1664 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! This integer division figures out which bin k is in since the kernel
  ! is set on a uniform grid.
!# 1668 "xc_vdW_DF.f90"
  k_i = INT(k/dk)
!# 1671 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Test to see if we are trying to interpolate a k that is one of the
  ! actual function points we have. The value is just the value of the
  ! function in that case.
!# 1676 "xc_vdW_DF.f90"
  IF ( MOD(k,dk) == 0 ) THEN
!# 1678 "xc_vdW_DF.f90"
     DO q1_i = 1, Nqs
        DO q2_i = 1, q1_i
!# 1681 "xc_vdW_DF.f90"
           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)
!# 1684 "xc_vdW_DF.f90"
        END DO
     END DO
!# 1687 "xc_vdW_DF.f90"
     RETURN
!# 1689 "xc_vdW_DF.f90"
  END IF
!# 1692 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! If we are not on a function point then we carry out the
  ! interpolation.
!# 1696 "xc_vdW_DF.f90"
  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
  C = (A**3-A)*dk**2/6.0D0
  D = (B**3-B)*dk**2/6.0D0
!# 1701 "xc_vdW_DF.f90"
  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i
!# 1704 "xc_vdW_DF.f90"
        kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
                        + (C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))
!# 1707 "xc_vdW_DF.f90"
        kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)
!# 1709 "xc_vdW_DF.f90"
     END DO
  END DO
!# 1712 "xc_vdW_DF.f90"
  END SUBROUTINE interpolate_kernel
!# 1721 "xc_vdW_DF.f90"
  ! ####################################################################
  !                         |                 |
  !                         |  VDW_DF_STRESS  |
  !                         |_________________|
!# 1726 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_stress (rho_valence, rho_core, nspin, sigma)
!# 1728 "xc_vdW_DF.f90"
  !! This routine calculates the vdW-DF contribution to the stress tensor.
!# 1731 "xc_vdW_DF.f90"
  use gvect,           ONLY : ngm, g
  USE cell_base,       ONLY : tpiba
!# 1734 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1736 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)     :: rho_valence(:,:)       !
  REAL(dp), INTENT(IN)     :: rho_core(:)            ! Input variables
  INTEGER,  INTENT(IN)     :: nspin                  !
  REAL(dp), INTENT(INOUT)  :: sigma(3,3)             !
!# 1741 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE    :: grad_rho(:,:)          !
  REAL(DP), ALLOCATABLE    :: total_rho(:)           ! Rho values
!# 1744 "xc_vdW_DF.f90"
  real(dp), allocatable    :: rho_up(:)              !
  real(dp), allocatable    :: rho_down(:)            !
!# 1747 "xc_vdW_DF.f90"
  real(dp), allocatable    :: grad_rho_up(:,:)       ! Gradient of the up charge density
                                                     ! Same format as grad_rho
  real(dp), allocatable    :: grad_rho_down(:,:)     ! Gradient of the down charge density
                                                     ! Same format as grad_rho
!# 1752 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE    :: q0(:)                  !
  REAL(DP), ALLOCATABLE    :: dq0_drho(:)            ! q-values
  REAL(DP), ALLOCATABLE    :: dq0_dgradrho(:)        !
  real(dp), allocatable    :: dq0_drho_up(:)         ! Derivative of the saturated q0
  real(dp), allocatable    :: dq0_drho_down(:)       ! with respect to the spin charge density
  real(dp), allocatable    :: dq0_dgradrho_up(:)     ! Derivative of the saturated q0 with respect
  real(dp), allocatable    :: dq0_dgradrho_down(:)   ! to the gradient of the spin charge sensity
!# 1760 "xc_vdW_DF.f90"
  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)            ! Thetas
  INTEGER                  :: i_proc, theta_i, l, m
!# 1763 "xc_vdW_DF.f90"
  REAL(DP)                 :: sigma_grad(3,3)
  REAL(DP)                 :: sigma_ker(3,3)
!# 1769 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Tests
!# 1772 "xc_vdW_DF.f90"
  IF ( inlc > 6 ) CALL errore( 'xc_vdW_DF', 'inlc not implemented', 1 )
!# 1774 "xc_vdW_DF.f90"
  IF ( nspin > 2 ) THEN
     CALL errore ('vdW_DF_stress', 'vdW stress not implemented for nspin > 2', 1)
  END IF
!# 1779 "xc_vdW_DF.f90"
  sigma(:,:)      = 0.0_DP
  sigma_grad(:,:) = 0.0_DP
  sigma_ker(:,:)  = 0.0_DP
!# 1784 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Allocations
!# 1787 "xc_vdW_DF.f90"
  ALLOCATE( total_rho(dfftp%nnr), grad_rho(3,dfftp%nnr), thetas(dfftp%nnr, Nqs), q0(dfftp%nnr) )
  ALLOCATE( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
!# 1790 "xc_vdW_DF.f90"
  IF (nspin==2) THEN
     ALLOCATE( rho_up(dfftp%nnr), rho_down(dfftp%nnr) )
     ALLOCATE( grad_rho_up(3,dfftp%nnr), grad_rho_down(3,dfftp%nnr) )
     ALLOCATE( dq0_drho_up (dfftp%nnr), dq0_dgradrho_up  (dfftp%nnr) )
     ALLOCATE( dq0_drho_down(dfftp%nnr), dq0_dgradrho_down(dfftp%nnr) )
  ENDIF
!# 1798 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Charge
!# 1801 "xc_vdW_DF.f90"
  total_rho = rho_valence(:,1) + rho_core(:)
!# 1803 "xc_vdW_DF.f90"
  IF (nspin == 2) THEN
    rho_up    = ( rho_valence(:,1) + rho_valence(:,2) + rho_core(:) )*0.5D0
    rho_down  = ( rho_valence(:,1) - rho_valence(:,2) + rho_core(:) )*0.5D0
  ENDIF
!# 1809 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.
!# 1812 "xc_vdW_DF.f90"
  CALL fft_gradient_r2r (dfftp, total_rho,  g, grad_rho)
!# 1814 "xc_vdW_DF.f90"
  IF (nspin == 2) THEN
     CALL fft_gradient_r2r (dfftp, rho_up,    g, grad_rho_up)
     CALL fft_gradient_r2r (dfftp, rho_down,  g, grad_rho_down)
  ENDIF
!# 1820 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get q0.
!# 1823 "xc_vdW_DF.f90"
  IF (nspin == 1) THEN
     CALL get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)
  ELSEIF (nspin == 2) THEN
     CALL get_q0_on_grid_spin ( total_rho, rho_up, rho_down, grad_rho, grad_rho_up, grad_rho_down, &
          q0, dq0_drho_up, dq0_drho_down, dq0_dgradrho_up, dq0_dgradrho_down, thetas)
  ENDIF
!# 1831 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Stress
!# 1834 "xc_vdW_DF.f90"
  IF (nspin == 1) THEN
     CALL vdW_DF_stress_gradient (total_rho, grad_rho, q0, dq0_drho, &
               dq0_dgradrho, thetas, sigma_grad)
  ELSEIF (nspin == 2) THEN
     CALL vdW_DF_stress_gradient_spin (total_rho, grad_rho_up, grad_rho_down, &
               q0, dq0_dgradrho_up, dq0_dgradrho_down, thetas, sigma_grad)
  ENDIF
!# 1843 "xc_vdW_DF.f90"
  CALL vdW_DF_stress_kernel (total_rho, q0, thetas, sigma_ker)
!# 1845 "xc_vdW_DF.f90"
  sigma = - (sigma_grad + sigma_ker)
!# 1847 "xc_vdW_DF.f90"
  DO l = 1, 3
     DO m = 1, l - 1
        sigma (m, l) = sigma (l, m)
     END DO
  END DO
!# 1853 "xc_vdW_DF.f90"
  DEALLOCATE( total_rho, grad_rho, thetas, q0 )
  DEALLOCATE( dq0_drho, dq0_dgradrho )
!# 1856 "xc_vdW_DF.f90"
  IF (nspin == 2) THEN
     DEALLOCATE( rho_up, rho_down )
     DEALLOCATE( grad_rho_up, grad_rho_down )
     DEALLOCATE( dq0_drho_up, dq0_drho_down, dq0_dgradrho_up, dq0_dgradrho_down )
  ENDIF
!# 1862 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_stress
!# 1871 "xc_vdW_DF.f90"
  ! ####################################################################
  !                     |                          |
  !                     |  VDW_DF_STRESS_GRADIENT  |
  !                     |__________________________|
!# 1876 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_stress_gradient (total_rho, grad_rho, q0, &
             dq0_drho, dq0_dgradrho, thetas, sigma)
!# 1880 "xc_vdW_DF.f90"
  USE gvect,                 ONLY : ngm, g, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2
!# 1883 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 1885 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)     :: total_rho(:)           !
  REAL(DP), INTENT(IN)     :: grad_rho(:, :)         !
  REAL(DP), INTENT(INOUT)  :: sigma(:,:)             !
  REAL(DP), INTENT(IN)     :: q0(:)                  ! Input variables
  REAL(DP), INTENT(IN)     :: dq0_drho(:)            !
  REAL(DP), INTENT(IN)     :: dq0_dgradrho(:)        !
  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)            !
!# 1893 "xc_vdW_DF.f90"
  COMPLEX(DP), ALLOCATABLE :: u_vdW(:,:)             !
!# 1895 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE    :: d2y_dx2(:,:)           !
  REAL(DP) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation
  REAL(DP) :: dq                                     !
!# 1899 "xc_vdW_DF.f90"
  INTEGER  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points
!# 1901 "xc_vdW_DF.f90"
  INTEGER  :: l, m
  REAL(DP) :: prefactor                              ! Final summation of sigma
  REAL(DP) :: grad2                                  ! magnitude of density gradient
!# 1905 "xc_vdW_DF.f90"
  INTEGER  :: i_proc, theta_i, i_grid, q_i, &        ! Iterators
              ix, iy, iz                             !
!# 1908 "xc_vdW_DF.f90"
  CHARACTER(LEN=1) :: intvar
!# 1913 "xc_vdW_DF.f90"
  ALLOCATE( d2y_dx2(Nqs, Nqs) )
  ALLOCATE( u_vdW(dfftp%nnr, Nqs) )
!# 1916 "xc_vdW_DF.f90"
  sigma(:,:) = 0.0_DP
  prefactor  = 0.0_DP
!# 1920 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get u in k-space.
!# 1923 "xc_vdW_DF.f90"
  CALL thetas_to_uk(thetas, u_vdW)
!# 1926 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get u in real space.
!# 1929 "xc_vdW_DF.f90"
  DO theta_i = 1, Nqs
     CALL invfft('Rho', u_vdW(:,theta_i), dfftp)
  END DO
!# 1934 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get the second derivatives for interpolating the P_i.
!# 1937 "xc_vdW_DF.f90"
  CALL initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
!# 1940 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Do the real space integration to obtain the stress component.
!# 1943 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
!# 1945 "xc_vdW_DF.f90"
     IF ( total_rho(i_grid) < epsr ) CYCLE
!# 1947 "xc_vdW_DF.f90"
     q_low = 1
     q_hi  = Nqs
     grad2 = sqrt( grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2 )
!# 1951 "xc_vdW_DF.f90"
     IF ( grad2 == 0.0_dp ) CYCLE
!# 1954 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in the q_mesh.
!# 1957 "xc_vdW_DF.f90"
     DO WHILE ( (q_hi - q_low) > 1)
!# 1959 "xc_vdW_DF.f90"
        q = INT((q_hi + q_low)/2)
!# 1961 "xc_vdW_DF.f90"
        IF (q_mesh(q) > q0(i_grid)) THEN
            q_hi = q
        ELSE
            q_low = q
        END IF
!# 1967 "xc_vdW_DF.f90"
     END DO
!# 1969 "xc_vdW_DF.f90"
     IF (q_hi == q_low) CALL errore('stress_vdW_gradient','qhi == qlow', 1)
!# 1971 "xc_vdW_DF.f90"
     dq = q_mesh(q_hi) - q_mesh(q_low)
     a  = (q_mesh(q_hi) - q0(i_grid))/dq
     b  = (q0(i_grid) - q_mesh(q_low))/dq
     c  = (a**3 - a)*dq**2/6.0D0
     d  = (b**3 - b)*dq**2/6.0D0
     e  = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f  = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
!# 1979 "xc_vdW_DF.f90"
     DO q_i = 1, Nqs
!# 1981 "xc_vdW_DF.f90"
        y(:)   = 0.0D0
        y(q_i) = 1.0D0
!# 1984 "xc_vdW_DF.f90"
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)
!# 1986 "xc_vdW_DF.f90"
        prefactor = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho(i_grid) / grad2
!# 1988 "xc_vdW_DF.f90"
        DO l = 1, 3
        DO m = 1, l
!# 1991 "xc_vdW_DF.f90"
            sigma (l, m) = sigma (l, m) -  e2 * prefactor * &
                           (grad_rho(l,i_grid) * grad_rho(m,i_grid))
        END DO
        END DO
!# 1996 "xc_vdW_DF.f90"
     END DO
!# 1998 "xc_vdW_DF.f90"
  END DO
!# 2000 "xc_vdW_DF.f90"
  CALL mp_sum(  sigma, intra_bgrp_comm )
!# 2002 "xc_vdW_DF.f90"
  CALL dscal (9, 1.0D0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)
!# 2004 "xc_vdW_DF.f90"
  DEALLOCATE( d2y_dx2, u_vdW )
!# 2006 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_stress_gradient
!# 2015 "xc_vdW_DF.f90"
  ! ####################################################################
  !                     |                               |
  !                     |  VDW_DF_STRESS_GRADIENT_SPIN  |
  !                     |_______________________________|
!# 2020 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_stress_gradient_spin (total_rho, grad_rho_up, grad_rho_down, &
                        q0, dq0_dgradrho_up, dq0_dgradrho_down, thetas, sigma)
!# 2023 "xc_vdW_DF.f90"
  !! This routine was implemented Per Hyldgaard (2019, GPL). No Waranties.
  !! Adapted from the original `nspin=1` code by Thonhauser and coauthors.
!# 2027 "xc_vdW_DF.f90"
  USE gvect,                 ONLY : ngm, g, gg, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2
!# 2030 "xc_vdW_DF.f90"
  implicit none
!# 2032 "xc_vdW_DF.f90"
  real(dp), intent(IN)     :: total_rho(:)           !
  real(dp), intent(IN)     :: grad_rho_up (:, :)     ! Input variables
  real(dp), intent(IN)     :: grad_rho_down(:, :)    !
  real(dp), intent(inout)  :: sigma(:,:)             !
  real(dp), intent(IN)     :: q0(:)                  !
  real(dp), intent(IN)     :: dq0_dgradrho_up(:)     !
  real(dp), intent(IN)     :: dq0_dgradrho_down(:)   !
  complex(dp), intent(IN)  :: thetas(:,:)            !
!# 2041 "xc_vdW_DF.f90"
  complex(dp), allocatable :: u_vdW(:,:)             !
!# 2043 "xc_vdW_DF.f90"
  real(dp), allocatable    :: d2y_dx2(:,:)           !
  real(dp) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation
  real(dp) :: dq                                     !
!# 2047 "xc_vdW_DF.f90"
  integer  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points
!# 2049 "xc_vdW_DF.f90"
  integer  :: l, m
  real(dp) :: prefactor_up, prefactor_down           ! Final summation of sigma
  real(dp) :: grad2_up, grad2_down                   ! Magnitude of density gradient
!# 2054 "xc_vdW_DF.f90"
  integer  :: i_proc, theta_i, i_grid, q_i, &        !
              ix, iy, iz                             ! Iterators
!# 2057 "xc_vdW_DF.f90"
  character(LEN=1) :: intvar
!# 2062 "xc_vdW_DF.f90"
  ALLOCATE( d2y_dx2(Nqs, Nqs) )
  ALLOCATE( u_vdW(dfftp%nnr, Nqs) )
!# 2065 "xc_vdW_DF.f90"
  sigma(:,:)      = 0.0_DP
  prefactor_up    = 0.0_DP
  prefactor_down  = 0.0_DP
!# 2070 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get u in k-space
!# 2073 "xc_vdW_DF.f90"
  call thetas_to_uk(thetas, u_vdW)
!# 2076 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get u in real space
!# 2079 "xc_vdW_DF.f90"
  DO theta_i = 1, Nqs
     CALL invfft('Rho', u_vdW(:,theta_i), dfftp)
  END DO
!# 2084 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get the second derivatives for interpolating the P_i.
!# 2087 "xc_vdW_DF.f90"
  CALL initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
!# 2090 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Do the real space integration to obtain the stress component.
!# 2093 "xc_vdW_DF.f90"
  DO i_grid = 1, dfftp%nnr
!# 2095 "xc_vdW_DF.f90"
     IF ( total_rho(i_grid) < epsr ) CYCLE
!# 2097 "xc_vdW_DF.f90"
     q_low = 1
     q_hi  = Nqs
     grad2_up = sqrt( grad_rho_up(1,i_grid)**2 &
                + grad_rho_up(2,i_grid)**2 + grad_rho_up(3,i_grid)**2 )
     grad2_down = sqrt( grad_rho_down(1,i_grid)**2 &
                  + grad_rho_down(2,i_grid)**2 + grad_rho_down(3,i_grid)**2 )
!# 2104 "xc_vdW_DF.f90"
     IF ( grad2_up == 0.0_dp .OR. grad2_down == 0.0_dp) CYCLE
!# 2106 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in the q_mesh.
!# 2109 "xc_vdW_DF.f90"
     DO WHILE ( (q_hi - q_low) > 1)
!# 2111 "xc_vdW_DF.f90"
        q = int((q_hi + q_low)/2)
!# 2113 "xc_vdW_DF.f90"
        IF (q_mesh(q) > q0(i_grid)) THEN
            q_hi  = q
        ELSE
            q_low = q
        END IF
!# 2119 "xc_vdW_DF.f90"
     END DO
!# 2121 "xc_vdW_DF.f90"
     IF (q_hi == q_low) CALL errore('vdW_DF_stress_gradient_spin','qhi == qlow',1)
!# 2123 "xc_vdW_DF.f90"
     dq = q_mesh(q_hi) - q_mesh(q_low)
!# 2125 "xc_vdW_DF.f90"
     a  = (q_mesh(q_hi) - q0(i_grid))/dq
     b  = (q0(i_grid) - q_mesh(q_low))/dq
     c  = (a**3 - a)*dq**2/6.0D0
     d  = (b**3 - b)*dq**2/6.0D0
     e  = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f  = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
!# 2132 "xc_vdW_DF.f90"
     DO q_i = 1, Nqs
!# 2134 "xc_vdW_DF.f90"
        y(:)   = 0.0D0
        y(q_i) = 1.0D0
!# 2137 "xc_vdW_DF.f90"
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)
!# 2139 "xc_vdW_DF.f90"
        prefactor_up = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho_up(i_grid) / grad2_up
        prefactor_down = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho_down(i_grid) / grad2_down
!# 2142 "xc_vdW_DF.f90"
        DO l = 1, 3
            DO m = 1, l
!# 2145 "xc_vdW_DF.f90"
                sigma (l, m) = sigma (l, m) -  e2 * prefactor_up * &
                               (grad_rho_up(l,i_grid) * grad_rho_up(m,i_grid))
                sigma (l, m) = sigma (l, m) -  e2 * prefactor_down * &
                               (grad_rho_down(l,i_grid) * grad_rho_down(m,i_grid))
            END DO
        END DO
!# 2152 "xc_vdW_DF.f90"
     END DO
!# 2154 "xc_vdW_DF.f90"
  END DO
!# 2156 "xc_vdW_DF.f90"
  CALL mp_sum(  sigma, intra_bgrp_comm )
!# 2158 "xc_vdW_DF.f90"
  CALL dscal (9, 1.d0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)
!# 2160 "xc_vdW_DF.f90"
  DEALLOCATE( d2y_dx2, u_vdW )
!# 2162 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_stress_gradient_spin
!# 2171 "xc_vdW_DF.f90"
  ! ####################################################################
  !                      |                        |
  !                      |  VDW_DF_STRESS_KERNEL  |
  !                      |________________________|
!# 2176 "xc_vdW_DF.f90"
  SUBROUTINE vdW_DF_stress_kernel (total_rho, q0, thetas, sigma)
!# 2179 "xc_vdW_DF.f90"
  USE gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, tpiba2
!# 2182 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 2184 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: q0(:)
  REAL(DP), INTENT(IN)    :: total_rho(:)
  REAL(DP), INTENT(INOUT) :: sigma(3,3)
  COMPLEX(DP), INTENT(IN) :: thetas(:,:)
!# 2189 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE   :: dkernel_of_dk(:,:)
!# 2191 "xc_vdW_DF.f90"
  INTEGER                 :: l, m, q1_i, q2_i , g_i
  INTEGER                 :: last_g, theta_i
  REAL(DP)                :: g2, ngmod2, g_kernel, G_multiplier
!# 2198 "xc_vdW_DF.f90"
  ALLOCATE( dkernel_of_dk(Nqs, Nqs) )
!# 2200 "xc_vdW_DF.f90"
  sigma(:,:) = 0.0_DP
!# 2202 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Integration in g-space
!# 2205 "xc_vdW_DF.f90"
  last_g = -1
!# 2207 "xc_vdW_DF.f90"
  G_multiplier = 1.0D0
!# 2209 "xc_vdW_DF.f90"
  IF ( gamma_only ) G_multiplier = 2.0D0
!# 2211 "xc_vdW_DF.f90"
  DO g_i = gstart, ngm
!# 2213 "xc_vdW_DF.f90"
     g2 = gg (g_i) * tpiba2
     g_kernel = SQRT(g2)
!# 2216 "xc_vdW_DF.f90"
     IF ( igtongl(g_i) .NE. last_g) THEN
!# 2218 "xc_vdW_DF.f90"
        CALL interpolate_Dkernel_Dk(g_kernel, dkernel_of_dk)  ! Gets the derivatives
        last_g = igtongl(g_i)
!# 2221 "xc_vdW_DF.f90"
     END IF
!# 2223 "xc_vdW_DF.f90"
     DO q2_i = 1, Nqs
     DO q1_i = 1, Nqs
        DO l = 1, 3
        DO m = 1, l
!# 2228 "xc_vdW_DF.f90"
           sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * e2 * thetas(dfftp%nl(g_i),q1_i) * &
                          dkernel_of_dk(q1_i,q2_i)*conjg(thetas(dfftp%nl(g_i),q2_i))* &
                          (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel
        END DO
        END DO
     END DO
     END DO
!# 2236 "xc_vdW_DF.f90"
     IF ( g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier
!# 2238 "xc_vdW_DF.f90"
  END DO
!# 2240 "xc_vdW_DF.f90"
  CALL mp_sum( sigma, intra_bgrp_comm )
!# 2242 "xc_vdW_DF.f90"
  DEALLOCATE( dkernel_of_dk )
!# 2244 "xc_vdW_DF.f90"
  END SUBROUTINE vdW_DF_stress_kernel
!# 2253 "xc_vdW_DF.f90"
  ! ####################################################################
  !                        |                        |
  !                        | INTERPOLATE_DKERNEL_DK |
  !                        |________________________|
!# 2258 "xc_vdW_DF.f90"
  SUBROUTINE interpolate_Dkernel_Dk (k, dkernel_of_dk)
!# 2261 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 2263 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: k                      ! Input value, the magnitude of the g-vector
                                                    ! for the current point.
!# 2266 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: dkernel_of_dk(Nqs,Nqs) ! An output array (allocated outside this
                                                    ! routine) that holds the interpolated value of
                                                    ! the kernel for each pair of q points (i.e. the
                                                    ! phi_alpha_beta of the Soler method.
!# 2271 "xc_vdW_DF.f90"
  INTEGER :: q1_i, q2_i, k_i                        ! Indexing variables
!# 2273 "xc_vdW_DF.f90"
  REAL(DP) :: A, B, dAdk, dBdk, dCdk, dDdk          ! Intermediate values for the interpolation.
!# 2278 "xc_vdW_DF.f90"
  IF ( k >= Nr_points*dk ) THEN
!# 2280 "xc_vdW_DF.f90"
     WRITE(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ", Nr_points*dk
     CALL errore('interpolate kernel', 'k value requested is out of range',1)
!# 2283 "xc_vdW_DF.f90"
  END IF
!# 2285 "xc_vdW_DF.f90"
  dkernel_of_dk = 0.0D0
!# 2287 "xc_vdW_DF.f90"
  k_i  = INT(k/dk)
!# 2289 "xc_vdW_DF.f90"
  A    = (dk*(k_i+1.0D0) - k)/dk
  B    = (k - dk*k_i)/dk
!# 2292 "xc_vdW_DF.f90"
  dAdk = -1.0D0/dk
  dBdk = 1.0D0/dk
  dCdk = -((3*A**2 -1.0D0)/6.0D0)*dk
  dDdk = ((3*B**2 -1.0D0)/6.0D0)*dk
!# 2297 "xc_vdW_DF.f90"
  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i
!# 2300 "xc_vdW_DF.f90"
        dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
                            + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)
!# 2303 "xc_vdW_DF.f90"
        dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)
!# 2305 "xc_vdW_DF.f90"
     END DO
  END DO
!# 2308 "xc_vdW_DF.f90"
  END SUBROUTINE interpolate_Dkernel_Dk
!# 2317 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |              |
  !                          | thetas_to_uk |
  !                          |______________|
!# 2323 "xc_vdW_DF.f90"
  SUBROUTINE thetas_to_uk (thetas, u_vdW)
!# 2326 "xc_vdW_DF.f90"
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega
!# 2329 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 2331 "xc_vdW_DF.f90"
  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)       ! On input this variable holds the theta functions
                                                ! (equation 8, SOLER) in the format
                                                ! thetas(grid_point, theta_i).
  COMPLEX(DP), INTENT(OUT) :: u_vdW(:,:)        ! On output this array holds u_alpha(k) =
                                                ! Sum_j[theta_beta(k)phi_alpha_beta(k)].
!# 2337 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE    :: kernel_of_k(:,:)  ! This array will hold the interpolated kernel
                                                ! values for each pair of q values in the q_mesh.
!# 2340 "xc_vdW_DF.f90"
  REAL(DP) :: g
  INTEGER  :: last_g, g_i, q1_i, q2_i, i_grid   ! Index variables.
!# 2343 "xc_vdW_DF.f90"
  COMPLEX(DP) :: theta(Nqs)                     ! Temporary storage vector used since we are
                                                ! overwriting the thetas array here.
!# 2349 "xc_vdW_DF.f90"
  ALLOCATE( kernel_of_k(Nqs, Nqs) )
!# 2351 "xc_vdW_DF.f90"
  u_vdW(:,:) = CMPLX(0.0_DP, 0.0_DP, kind=dp)
!# 2353 "xc_vdW_DF.f90"
  last_g = -1
!# 2355 "xc_vdW_DF.f90"
  DO g_i = 1, ngm
!# 2357 "xc_vdW_DF.f90"
     IF ( igtongl(g_i) .ne. last_g) THEN
!# 2359 "xc_vdW_DF.f90"
        g = SQRT(gl(igtongl(g_i))) * tpiba
        CALL interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
!# 2363 "xc_vdW_DF.f90"
     END IF
!# 2365 "xc_vdW_DF.f90"
     theta = thetas(dfftp%nl(g_i),:)
!# 2367 "xc_vdW_DF.f90"
     DO q2_i = 1, Nqs
        DO q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i) = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        END DO
     END DO
!# 2373 "xc_vdW_DF.f90"
  END Do
!# 2375 "xc_vdW_DF.f90"
  IF ( gamma_only ) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))
!# 2377 "xc_vdW_DF.f90"
  DEALLOCATE( kernel_of_k )
!# 2379 "xc_vdW_DF.f90"
  END SUBROUTINE thetas_to_uk
!# 2388 "xc_vdW_DF.f90"
  ! ####################################################################
  !                           |                 |
  !                           | GENERATE_KERNEL |
  !                           |_________________|
!# 2393 "xc_vdW_DF.f90"
  SUBROUTINE generate_kernel
!# 2395 "xc_vdW_DF.f90"
  !! This routine calculates the vdW-DF kernel.
  !
  ! The original definition of the kernel function is given in DION
  ! equations 14-16. The Soler method makes the kernel function a
  ! function of only 1 variable (r) by first putting it in the form
  ! phi(q1*r, q2*r). Then, the q-dependence is removed by expanding the
  ! function in a special way (see SOLER equation 3). This yields a
  ! separate function for each pair of q points that is a function of r
  ! alone. There are (Nqs^2+Nqs)/2 unique functions, where Nqs is the
  ! number of q points used. In the Soler method, the kernel is first
  ! made in the form phi(d1, d2) but this is not done here. It was found
  ! that, with q's chosen judiciously ahead of time, the kernel and the
  ! second derivatives required for interpolation could be tabulated
  ! ahead of time for faster use of the vdW-DF functional. Through
  ! testing we found no need to soften the kernel and correct for this
  ! later (see SOLER eqations 6-7).
  !
  ! The algorithm employed here is "embarrassingly parallel," meaning
  ! that it parallelizes very well up to (Nqs^2+Nqs)/2 processors,
  ! where, again, Nqs is the number of q points chosen. However,
  ! parallelization on this scale is unnecessary. In testing the code
  ! runs in under a minute on 16 Intel Xeon processors.
  !
  ! IMPORTANT NOTICE: Results are very sensitive to compilation details.
  ! In particular, the usage of FMA (Fused Multiply-and-Add)
  ! instructions used by modern CPUs such as AMD Interlagos (Bulldozer)
  ! and Intel Ivy Bridge may affect quite heavily some components of the
  ! kernel (communication by Ake Sandberg, Umea University). In practice
  ! this should not be a problem, since most affected elements are the
  ! less relevant ones.
  !
  ! For the calculation of the kernel we have benefited from access to
  ! earlier vdW-DF implementation into PWscf and ABINIT, written by Timo
  ! Thonhauser, Valentino Cooper, and David Langreth. These codes, in
  ! turn, benefited from earlier codes written by Maxime Dion and Henrik
  ! Rydberg.
!# 2433 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 2435 "xc_vdW_DF.f90"
  INTEGER  :: a_i, b_i, q1_i, q2_i, r_i
  ! Indexing variables
!# 2438 "xc_vdW_DF.f90"
  REAL(DP) :: weights( Nintegration_points )
  ! Array to hold dx values for the Gaussian-Legendre integration of the kernel.
!# 2441 "xc_vdW_DF.f90"
  REAL(DP) :: sin_a( Nintegration_points ), cos_a( Nintegration_points )
  ! Sine and cosine values of the aforementioned points a.
!# 2444 "xc_vdW_DF.f90"
  REAL(DP) :: d1, d2, d, integral
  ! Intermediate values
!# 2447 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! The following variables control the parallel environment.
!# 2450 "xc_vdW_DF.f90"
  INTEGER :: my_start_q, my_end_q, Ntotal
  ! Starting and ending q value for each  processor, also the total
  ! number of calculations to do, i.e. (Nqs^2 + Nqs)/2.
!# 2454 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE :: phi(:,:), phi_deriv(:,:)
  ! Arrays to store the kernel functions and their second derivatives.
  ! They are stored as phi(radial_point, idx).
!# 2458 "xc_vdW_DF.f90"
  INTEGER, ALLOCATABLE  :: indices(:,:), proc_indices(:,:)
  ! Indices holds the values of q1 and q2 as partitioned out to the
  ! processors. It is an Ntotal x 2 array stored as indices(index of
  ! point number, q1:q2). Proc_indices holds the section of the indices
  ! array that is assigned to each processor. This is a Nproc x 2
  ! array, stored as proc_indices(processor_number,
  ! starting_index:ending_index)
!# 2466 "xc_vdW_DF.f90"
  INTEGER :: Nper, Nextra, start_q, end_q
  ! Baseline number of jobs per processor, number of processors that
  ! get an extra job in case the number of jobs doesn't split evenly
  ! over the number of processors, starting index into the indices
  ! array, ending index into the indices array.
!# 2472 "xc_vdW_DF.f90"
  INTEGER :: nproc, mpime
  ! Number or procs, rank of current processor.
!# 2475 "xc_vdW_DF.f90"
  INTEGER :: proc_i, my_Nqs
!# 2480 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Start the timer.
!# 2483 "xc_vdW_DF.f90"
  CALL start_clock ( 'vdW_kernel' )
!# 2486 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! The total number of phi_alpha_beta functions that have to be
  ! calculated.
!# 2490 "xc_vdW_DF.f90"
  Ntotal = (Nqs**2 + Nqs)/2
  ALLOCATE ( indices(Ntotal, 2) )
!# 2494 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! This part fills in the indices array. It just loops through the q1
  ! and q2 values and stores them. Sections of this array will be
  ! assigned to each of the processors later.
!# 2499 "xc_vdW_DF.f90"
  idx = 1
!# 2501 "xc_vdW_DF.f90"
  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i
        indices(idx, 1) = q1_i
        indices(idx, 2) = q2_i
        idx = idx + 1
     END DO
  END DO
!# 2510 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Figure out the baseline number of functions to be calculated by each
  ! processor and how many processors get one extra job.
!# 2514 "xc_vdW_DF.f90"
  nproc  = mp_size( intra_image_comm )
  mpime  = mp_rank( intra_image_comm )
  Nper   = Ntotal/nproc
  Nextra = MOD(Ntotal, nproc)
!# 2519 "xc_vdW_DF.f90"
  ALLOCATE( proc_indices(nproc, 2) )
!# 2521 "xc_vdW_DF.f90"
  start_q = 0
  end_q   = 0
!# 2525 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Loop over all the processors and figure out which section of the
  ! indices array each processor should do. All processors figure this
  ! out for every processor so there is no need to communicate results.
!# 2530 "xc_vdW_DF.f90"
  DO proc_i = 1, nproc
!# 2532 "xc_vdW_DF.f90"
     start_q = end_q + 1
     end_q   = start_q + (Nper - 1)
     IF (proc_i <= Nextra) end_q = end_q + 1
!# 2536 "xc_vdW_DF.f90"
     ! This is to prevent trouble if number of processors exceeds Ntotal.
     IF ( proc_i > Ntotal ) THEN
        start_q    = Ntotal
        end_q      = Ntotal
     END IF
!# 2542 "xc_vdW_DF.f90"
     IF ( proc_i == (mpime+1) ) THEN
        my_start_q = start_q
        my_end_q   = end_q
     END IF
!# 2547 "xc_vdW_DF.f90"
     proc_indices(proc_i, 1) = start_q
     proc_indices(proc_i, 2) = end_q
!# 2550 "xc_vdW_DF.f90"
  END DO
!# 2553 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Store how many jobs are assigned to me.
!# 2556 "xc_vdW_DF.f90"
  my_Nqs    = my_end_q - my_start_q + 1
  ALLOCATE( phi( 0:Nr_points, my_Nqs ), phi_deriv( 0:Nr_points, my_Nqs ) )
!# 2559 "xc_vdW_DF.f90"
  phi       = 0.0D0
  phi_deriv = 0.0D0
  kernel    = 0.0D0
  d2phi_dk2 = 0.0D0
!# 2565 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Find the integration points we are going to use in the
  ! Gaussian-Legendre integration.
!# 2569 "xc_vdW_DF.f90"
  CALL prep_gaussian_quadrature( weights )
!# 2572 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Get a, a^2, sin(a), cos(a) and the weights for the Gaussian-Legendre
  ! integration.
!# 2576 "xc_vdW_DF.f90"
  DO a_i=1, Nintegration_points
     a_points (a_i) = TAN( a_points(a_i) )
     a_points2(a_i) = a_points(a_i)**2
     weights(a_i)   = weights(a_i)*( 1 + a_points2(a_i) )
     cos_a(a_i)     = COS( a_points(a_i) )
     sin_a(a_i)     = SIN( a_points(a_i) )
  END DO
!# 2585 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Calculate the value of the W function defined in DION equation 16
  ! for each value of a and b.
!# 2589 "xc_vdW_DF.f90"
  DO a_i = 1, Nintegration_points
  DO b_i = 1, Nintegration_points
     W_ab(a_i, b_i) = 2.0D0 * weights(a_i)*weights(b_i) * (           &
        (3.0D0-a_points2(a_i))*a_points(b_i) *sin_a(a_i)*cos_a(b_i) + &
        (3.0D0-a_points2(b_i))*a_points(a_i) *cos_a(a_i)*sin_a(b_i) + &
        (a_points2(a_i)+a_points2(b_i)-3.0D0)*sin_a(a_i)*sin_a(b_i) - &
        3.0D0*a_points(a_i)*a_points(b_i)*cos_a(a_i)*cos_a(b_i) )   / &
        (a_points(a_i)*a_points(b_i) )
  END DO
  END DO
!# 2601 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! vdW-DF analysis tool as described in PRB 97, 085115 (2018).
!# 2604 "xc_vdW_DF.f90"
  IF      ( vdW_DF_analysis == 1 ) THEN
!# 2606 "xc_vdW_DF.f90"
     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = weights(a_i)*weights(b_i) *                  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO
!# 2613 "xc_vdW_DF.f90"
  ELSE IF ( vdW_DF_analysis == 2 ) THEN
!# 2615 "xc_vdW_DF.f90"
     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = W_ab(a_i, b_i) - weights(a_i)*weights(b_i) *  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO
!# 2622 "xc_vdW_DF.f90"
  END IF
!# 2625 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Now, we loop over all the pairs (q1,q2) that are assigned to us and
  ! perform our calculations.
!# 2629 "xc_vdW_DF.f90"
  DO idx = 1, my_Nqs
!# 2631 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! First, get the value of phi(q1*r, q2*r) for each r and the
     ! particular values of q1 and q2 we are using.
!# 2635 "xc_vdW_DF.f90"
     DO r_i = 1, Nr_points
        d1  = q_mesh( indices(idx+my_start_q-1, 1) ) * dr * r_i
        d2  = q_mesh( indices(idx+my_start_q-1, 2) ) * dr * r_i
        phi(r_i, idx) = phi_value(d1, d2)
     END DO
!# 2642 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Now, perform a radial FFT to turn our phi_alpha_beta(r) into
     ! phi_alpha_beta(k) needed for SOLER equation 8.
!# 2646 "xc_vdW_DF.f90"
     CALL radial_fft( phi(:,idx) )
!# 2649 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Determine the spline interpolation coefficients for the Fourier
     ! transformed kernel function.
!# 2653 "xc_vdW_DF.f90"
     CALL set_up_splines( phi(:, idx), phi_deriv(:, idx) )
!# 2655 "xc_vdW_DF.f90"
  END DO
!# 2658 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Finally, we collect the results after letting everybody catch up.
!# 2661 "xc_vdW_DF.f90"
  CALL mp_barrier( intra_image_comm )
!# 2663 "xc_vdW_DF.f90"
  DO proc_i = 0, nproc-1
!# 2665 "xc_vdW_DF.f90"
     IF ( proc_i >= Ntotal ) EXIT
!# 2667 "xc_vdW_DF.f90"
     CALL mp_get ( phi      , phi      , mpime, 0, proc_i, 0, intra_image_comm )
     CALL mp_get ( phi_deriv, phi_deriv, mpime, 0, proc_i, 0, intra_image_comm )
!# 2670 "xc_vdW_DF.f90"
     IF ( mpime == 0 ) THEN
!# 2672 "xc_vdW_DF.f90"
        DO idx = proc_indices(proc_i+1,1), proc_indices(proc_i+1,2)
           q1_i = indices(idx, 1)
           q2_i = indices(idx, 2)
           kernel    (:, q1_i, q2_i) = phi       (:, idx - proc_indices(proc_i+1,1) + 1)
           d2phi_dk2 (:, q1_i, q2_i) = phi_deriv (:, idx - proc_indices(proc_i+1,1) + 1)
           kernel    (:, q2_i, q1_i) = kernel    (:, q1_i, q2_i)
           d2phi_dk2 (:, q2_i, q1_i) = d2phi_dk2 (:, q1_i, q2_i)
        END DO
!# 2681 "xc_vdW_DF.f90"
     END IF
!# 2683 "xc_vdW_DF.f90"
  END DO
!# 2685 "xc_vdW_DF.f90"
  CALL mp_bcast ( kernel   , 0, intra_image_comm )
  CALL mp_bcast ( d2phi_dk2, 0, intra_image_comm )
!# 2689 "xc_vdW_DF.f90"
  DEALLOCATE( indices, proc_indices, phi, phi_deriv )
!# 2692 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Stop the timer.
!# 2695 "xc_vdW_DF.f90"
  CALL stop_clock ( 'vdW_kernel' )
!# 2697 "xc_vdW_DF.f90"
  END SUBROUTINE generate_kernel
!# 2706 "xc_vdW_DF.f90"
  ! ####################################################################
  !                    |                            |
  !                    |  PREP_GAUSSIAN_QUADRATURE  |
  !                    |____________________________|
!# 2711 "xc_vdW_DF.f90"
  SUBROUTINE prep_gaussian_quadrature( weights )
!# 2713 "xc_vdW_DF.f90"
  !! Routine to calculate the points and weights for the
  !! Gaussian-Legendre integration. This routine is modeled after the
  !! routine GAULEG from NUMERICAL_RECIPES.
!# 2718 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: weights(:)
  ! The points and weights for the Gaussian-Legendre integration.
!# 2721 "xc_vdW_DF.f90"
  INTEGER  :: Npoints
  ! The number of points we actually have to calculate. The rest will
  ! be obtained from symmetry.
!# 2725 "xc_vdW_DF.f90"
  REAL(DP) :: poly_1, poly_2, poly_3
  ! Temporary storage for Legendre polynomials.
!# 2728 "xc_vdW_DF.f90"
  INTEGER  :: i_point, i_poly
  ! Indexing variables
!# 2731 "xc_vdW_DF.f90"
  REAL(DP) :: root, dp_dx, last_root
  ! The value of the root of a given Legendre polynomial, the derivative
  ! of the polynomial at that root and the value of the root in the last
  ! iteration (to check for convergence of Newton's method).
!# 2736 "xc_vdW_DF.f90"
  real(dp) :: midpoint, length
  ! The middle of the x-range and the length to that point.
!# 2742 "xc_vdW_DF.f90"
  Npoints  = (Nintegration_points + 1)/2
  midpoint = 0.5D0 * ( ATAN(a_min) + ATAN(a_max) )
  length   = 0.5D0 * ( ATAN(a_max) - ATAN(a_min) )
!# 2746 "xc_vdW_DF.f90"
  DO i_point = 1, Npoints
     ! -----------------------------------------------------------------
     ! Make an initial guess for the root.
!# 2750 "xc_vdW_DF.f90"
     root = COS(DBLE(pi*(i_point - 0.25D0)/(Nintegration_points + 0.5D0)))
!# 2752 "xc_vdW_DF.f90"
     DO
        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial,
        ! evaluated at the approximate root. See NUMERICAL_RECIPES.
!# 2757 "xc_vdW_DF.f90"
        poly_1 = 1.0D0
        poly_2 = 0.0D0
!# 2760 "xc_vdW_DF.f90"
        DO i_poly = 1, Nintegration_points
!# 2762 "xc_vdW_DF.f90"
           poly_3 = poly_2
           poly_2 = poly_1
           poly_1 = ((2.0D0 * i_poly - 1.0D0)*root*poly_2 - (i_poly-1.0D0)*poly_3)/i_poly
!# 2766 "xc_vdW_DF.f90"
        END DO
!# 2769 "xc_vdW_DF.f90"
        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial.
        ! Find the derivative of the polynomial and use it in Newton's
        ! method to refine our guess for the root.
!# 2774 "xc_vdW_DF.f90"
        dp_dx = Nintegration_points * (root*poly_1 - poly_2)/(root**2 - 1.0D0)
!# 2776 "xc_vdW_DF.f90"
        last_root = root
        root      = last_root - poly_1/dp_dx
!# 2780 "xc_vdW_DF.f90"
        ! --------------------------------------------------------------
        ! Check for convergence.
!# 2783 "xc_vdW_DF.f90"
        IF (abs(root - last_root) <= 1.0D-14) EXIT
!# 2785 "xc_vdW_DF.f90"
     END DO
!# 2788 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Fill in the array of evaluation points.
!# 2791 "xc_vdW_DF.f90"
     a_points(i_point) = midpoint - length*root
     a_points(Nintegration_points + 1 - i_point) = midpoint + length*root
!# 2795 "xc_vdW_DF.f90"
     ! -----------------------------------------------------------------
     ! Fill in the array of weights.
!# 2798 "xc_vdW_DF.f90"
     weights(i_point) = 2.0D0 * length/((1.0D0 - root**2)*dp_dx**2)
     weights(Nintegration_points + 1 - i_point) = weights(i_point)
!# 2801 "xc_vdW_DF.f90"
  END DO
!# 2803 "xc_vdW_DF.f90"
  END SUBROUTINE prep_gaussian_quadrature
!# 2812 "xc_vdW_DF.f90"
  ! ####################################################################
  !                            |             |
  !                            |  PHI_VALUE  |
  !                            |_____________|
!# 2817 "xc_vdW_DF.f90"
  REAL(DP) FUNCTION phi_value(d1, d2)
!# 2819 "xc_vdW_DF.f90"
  !! This function returns the value of the kernel calculated via DION
  !! equation 14.
!# 2823 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN) :: d1, d2
  ! The point at which to evaluate the kernel. d1 = q1*r and d2 = q2*r.
!# 2826 "xc_vdW_DF.f90"
  REAL(DP) :: w, x, y, z, T
  ! Intermediate values
!# 2829 "xc_vdW_DF.f90"
  REAL(DP) :: nu(Nintegration_points), nu1(Nintegration_points)
  ! Defined in the discussio below equation 16 of DION.
!# 2832 "xc_vdW_DF.f90"
  INTEGER  :: a_i, b_i
  ! Indexing variables
!# 2838 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Loop over all integration points and calculate the value of the nu
  ! functions defined in the discussion below equation 16 in DION.
!# 2842 "xc_vdW_DF.f90"
  DO a_i = 1, Nintegration_points
     nu(a_i)  = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d1 ))
     nu1(a_i) = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d2 ))
  END DO
!# 2848 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Carry out the integration of DION equation 13.
!# 2851 "xc_vdW_DF.f90"
  phi_value = 0.0D0
!# 2853 "xc_vdW_DF.f90"
  DO a_i = 1, Nintegration_points
     w = nu(a_i)
     y = nu1(a_i)
     DO b_i = 1, Nintegration_points
        x = nu(b_i)
        z = nu1(b_i)
        T = (1.0D0/(w+x) + 1.0D0/(y+z))*(1.0D0/((w+y)*(x+z)) + 1.0D0/((w+z)*(y+x)))
        phi_value = phi_value + T * W_ab(a_i, b_i)
     END DO
  END DO
!# 2864 "xc_vdW_DF.f90"
  phi_value = 1.0D0/pi**2*phi_value
!# 2866 "xc_vdW_DF.f90"
  END FUNCTION phi_value
!# 2875 "xc_vdW_DF.f90"
  ! ####################################################################
  !                            |              |
  !                            |  RADIAL_FFT  |
  !                            |______________|
!# 2880 "xc_vdW_DF.f90"
  SUBROUTINE radial_fft(phi)
!# 2882 "xc_vdW_DF.f90"
  !! This subroutine performs a radial Fourier transform on the
  !! real-space kernel functions.
  !!
  !! Basically, this is just `4*pi*r^2*phi*sin(k*r)/(k*r) dr`
  !! integrated from `0` to `r_max`. That is, it is the kernel function
  !! `phi` integrated with the 0^th spherical Bessel function radially,
  !! with a `4*pi` assumed from angular integration since we have
  !! spherical symmetry. The spherical symmetry comes in because the
  !! kernel function depends only on the magnitude of the vector between
  !! two points. The integration is done using the trapezoid rule.
!# 2894 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  ! On input holds the real-space function phi_q1_q2(r).
  ! On output hold the reciprocal-space function phi_q1_q2(k).
!# 2898 "xc_vdW_DF.f90"
  REAL(DP) :: phi_k(0:Nr_points)
  ! Temporary storage for phi_q1_q2(k).
!# 2901 "xc_vdW_DF.f90"
  INTEGER  :: k_i, r_i
  ! Indexing variables
!# 2904 "xc_vdW_DF.f90"
  REAL(DP) :: r, k
  ! The real and reciprocal space points.
!# 2910 "xc_vdW_DF.f90"
  phi_k = 0.0D0
!# 2912 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Handle the k=0 point separately.
!# 2915 "xc_vdW_DF.f90"
  DO r_i = 1, Nr_points
     r        = r_i * dr
     phi_k(0) = phi_k(0) + phi(r_i)*r**2
  END DO
!# 2921 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Subtract half of the last value off because of the trapezoid rule.
!# 2924 "xc_vdW_DF.f90"
  phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)
!# 2927 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Integration for the rest of the k-points.
!# 2930 "xc_vdW_DF.f90"
  DO k_i = 1, Nr_points
     k = k_i * dk
     DO r_i = 1, Nr_points
        r          = r_i * dr
        phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * SIN(k*r) / k
     END DO
     phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r * SIN(k*r) / k
  END DO
!# 2940 "xc_vdW_DF.f90"
  ! --------------------------------------------------------------------
  ! Add in the 4*pi and the dr factor for the integration.
!# 2943 "xc_vdW_DF.f90"
  phi = 4.0D0 * pi * phi_k * dr
!# 2945 "xc_vdW_DF.f90"
  END SUBROUTINE radial_fft
!# 2954 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |                  |
  !                          |  SET UP SPLINES  |
  !                          |__________________|
!# 2959 "xc_vdW_DF.f90"
  SUBROUTINE set_up_splines(phi, D2)
!# 2961 "xc_vdW_DF.f90"
  !! This subroutine accepts a function (`phi`) and finds at each point
  !! the second derivative (`D2`) for use with spline interpolation. This
  !! function assumes we are using the expansion described in SOLER
  !! equation 3. That is, the derivatives are those needed to interpolate
  !! Kronecker delta functions at each of the `q` values. Other than some
  !! special modification to speed up the algorithm in our particular
  !! case, this algorithm is taken directly from NUMERICAL_RECIPES.
!# 2970 "xc_vdW_DF.f90"
  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  ! The k-space kernel function for a particular q1 and q2.
!# 2973 "xc_vdW_DF.f90"
  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  ! The second derivatives to be used in the interpolation expansion
  ! (SOLER equation 3).
!# 2977 "xc_vdW_DF.f90"
  REAL(DP), ALLOCATABLE   :: temp_array(:)         ! Temporary storage
  REAL(DP)                :: temp_1, temp_2
!# 2980 "xc_vdW_DF.f90"
  INTEGER  :: r_i
  ! Indexing variable
!# 2986 "xc_vdW_DF.f90"
  ALLOCATE( temp_array(0:Nr_points) )
!# 2988 "xc_vdW_DF.f90"
  D2         = 0
  temp_array = 0
!# 2991 "xc_vdW_DF.f90"
  DO r_i = 1, Nr_points - 1
     temp_1  = DBLE(r_i - (r_i - 1))/DBLE( (r_i + 1) - (r_i - 1) )
     temp_2  = temp_1 * D2(r_i-1) + 2.0D0
     D2(r_i) = (temp_1 - 1.0D0)/temp_2
     temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/DBLE( dk*((r_i+1) - r_i) ) - &
          ( phi(r_i) - phi(r_i-1))/DBLE( dk*(r_i - (r_i-1)) )
     temp_array(r_i) = (6.0D0*temp_array(r_i)/DBLE( dk*((r_i+1) - (r_i-1)) )-&
          temp_1*temp_array(r_i-1))/temp_2
  END DO
!# 3001 "xc_vdW_DF.f90"
  D2(Nr_points) = 0.0D0
  DO  r_i = Nr_points-1, 0, -1
     D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)
  END DO
!# 3006 "xc_vdW_DF.f90"
  DEALLOCATE( temp_array )
!# 3008 "xc_vdW_DF.f90"
  END SUBROUTINE set_up_splines
!# 3017 "xc_vdW_DF.f90"
  ! ####################################################################
  !                          |            |
  !                          |  VDW_INFO  |
  !                          |____________|
!# 3022 "xc_vdW_DF.f90"
  SUBROUTINE vdW_info
!# 3024 "xc_vdW_DF.f90"
  IMPLICIT NONE
!# 3029 "xc_vdW_DF.f90"
  WRITE(stdout,'(/)')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% You are using vdW-DF, which was implemented by the Thonhauser group. %")')
  WRITE(stdout,'(5x,"% Please cite the following two papers that made this development      %")')
  WRITE(stdout,'(5x,"% possible and the two reviews that describe the various versions:     %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   T. Thonhauser et al., PRL 115, 136402 (2015).                      %")')
  WRITE(stdout,'(5x,"%   T. Thonhauser et al., PRB 76,  125112 (2007).                      %")')
  WRITE(stdout,'(5x,"%   K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).             %")')
  WRITE(stdout,'(5x,"%   D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009). %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% If you are calculating stress with vdW-DF, please also cite:         %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   R. Sabatini et al., J. Phys.: Condens. Matter 24, 424209 (2012).   %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   And for spin-polarized stress cases also:                          %")')
  WRITE(stdout,'(5x,"%   C.M. Frostenson et al., Electr. Struct. 4, 014001 (2022).          %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'()')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%  vdW-DF NEWS:                                                        %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%  * vdW-DF3-mc is now available. DOI: 10.1103/hp9d-4kpf               %")')
  WRITE(stdout,'(5x,"%    use with input_dft = ''vdW-DF3-mc''                                 %")')
  WRITE(stdout,'(5x,"%  * vdW-DF3 is now available. DOI: 10.1021/acs.jctc.0c00471           %")')
  WRITE(stdout,'(5x,"%    use with input_dft = ''vdW-DF3-opt1'' or ''vdW-DF3-opt2''             %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%  * Unscreened and range-separated hybrid vdW-DF-cx functionals       %")')
  WRITE(stdout,'(5x,"%    DOI: 10.1063/1.4986522 and 10.1088/1361-648X/ac2ad2               %")')
  WRITE(stdout,'(5x,"%    use with input_dft = ''vdW-DF-cx0'' and ''vdW-DF-ahcx''               %")')
  WRITE(stdout,'(5x,"%  * Unscreened and range-separated hybrid vdW-DF2-b86r functionals    %")')
  WRITE(stdout,'(5x,"%    DOI: 10.1063/1.4986522 and DOI: 10.1103/PhysRevX.12.041003        %")')
  WRITE(stdout,'(5x,"%    use with input_dft = ''vdW-DF2-br0'' and ''vdW-DF2-ahbr''             %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(/)')
!# 3070 "xc_vdW_DF.f90"
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(5x,"Carrying out vdW-DF run using the following parameters:")')
     WRITE(stdout,'(5X,A,I3,A,I5,A,F8.3)' ) "Nqs    = ", Nqs, "  Npoints = ", Nr_points, &
                  "  r_max = ", r_max
     WRITE(stdout,'(5X,"q_mesh =",4F12.8)') (q_mesh(idx), idx=1, 4)
     WRITE(stdout,'(13X,4F12.8)') (q_mesh(idx), idx=5, Nqs)
  END IF
!# 3078 "xc_vdW_DF.f90"
  END SUBROUTINE
!# 3083 "xc_vdW_DF.f90"
END MODULE vdW_DF
