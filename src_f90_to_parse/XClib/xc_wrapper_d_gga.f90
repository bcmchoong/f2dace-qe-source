!# 1 "xc_wrapper_d_gga.f90"
!# 2 "xc_wrapper_d_gga.f90"
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss, gpu_args_ )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sp
  !! number of spin components
  REAL(DP), INTENT(IN) :: r_in(length,sp)
  !! charge density
  REAL(DP), INTENT(IN) :: g_in(length,3,sp)
  !! gradient
  REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                           dvxc_ss(length,sp,sp)
  LOGICAL, OPTIONAL, INTENT(IN) :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  !
  gpu_args = .FALSE.
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    CALL dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( r_in, g_in ), copyout( dvxc_rr, dvxc_sr, dvxc_ss )
    CALL dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    !$acc end data
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE constants_l,          ONLY: e2
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda,   &
                                  ishybrid, exx_started, exx_fraction
  USE qe_drivers_d_gga
!# 71 "xc_wrapper_d_gga.f90"
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sp
  !! number of spin components
  REAL(DP), INTENT(IN) :: r_in(length,sp)
  !! charge density
  REAL(DP), INTENT(IN) :: g_in(length,3,sp)
  !! gradient
  REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                           dvxc_ss(length,sp,sp)
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:)
  REAL(DP), ALLOCATABLE :: vrrc(:,:), vsrc(:,:), vssc(:), vrzc(:,:)
  !
  INTEGER :: fkind
!# 101 "xc_wrapper_d_gga.f90"
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: k, length_dlxc
  REAL(DP) :: rht, zeta, xcoef
  REAL(DP), ALLOCATABLE :: sigma(:)
  REAL(DP), PARAMETER :: small=1.E-10_DP, rho_trash=0.5_DP
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-6
  !
  !$acc data present( r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !
  IF ( ANY(.NOT.is_libxc(3:4)) ) THEN
    rho_threshold_gga = small ;  grho_threshold_gga = small
  ENDIF
  !
  !$acc kernels
  dvxc_rr(:,:,:) = 0._DP
  dvxc_sr(:,:,:) = 0._DP
  dvxc_ss(:,:,:) = 0._DP
  !$acc end kernels
  !
  fkind = -1
  fkind_is_XC = .FALSE.
  !
!# 273 "xc_wrapper_d_gga.f90"
  !
  ! ... QE DERIVATIVE FOR EXCHANGE AND CORRELATION
  !
  IF ( ((.NOT.is_libxc(3).AND.igcx/=0) .OR. (.NOT.is_libxc(4).AND.igcc/=0)) &
        .AND. (.NOT.fkind_is_XC) ) THEN
    !
    ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
    ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
    !$acc data create( vrrx, vsrx, vssx, vrrc, vsrc, vssc )
    !
    IF ( sp == 1 ) THEN
       !
       IF (.NOT. ALLOCATED(sigma)) THEN
         ALLOCATE( sigma(length) )
         !$acc enter data create(sigma)
       ENDIF
       !
       !$acc parallel loop
       DO k = 1, length
         sigma(k) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
       ENDDO
       !
       CALL dgcxc_unpol( length, r_in(:,1), sigma, vrrx(:,1), vsrx(:,1), vssx(:,1), &
                         vrrc(:,1), vsrc(:,1), vssc )
       !
       !$acc parallel loop
       DO k = 1, length
         dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * (vrrx(k,1) + vrrc(k,1))
         dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * (vsrx(k,1) + vsrc(k,1))
         dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * (vssx(k,1) + vssc(k)  )
       ENDDO
       !
    ELSEIF ( sp == 2 ) THEN
       !
       ALLOCATE( vrzc(length,sp) )
       !$acc data create( vrzc )
       !
       CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
       !
       !$acc parallel loop
       DO k = 1, length
         rht = r_in(k,1) + r_in(k,2)
         IF (rht > epsr) THEN
           zeta = (r_in(k,1) - r_in(k,2))/rht
           !
           dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2*(vrrx(k,1) + vrrc(k,1) + &
                                                 vrzc(k,1)*(1.d0 - zeta)/rht)
           dvxc_rr(k,1,2) = dvxc_rr(k,1,2) + e2*(vrrc(k,1) - vrzc(k,1)*(1.d0 + zeta)/rht)
           dvxc_rr(k,2,1) = dvxc_rr(k,2,1) + e2*(vrrc(k,2) + vrzc(k,2)*(1.d0 - zeta)/rht)
           dvxc_rr(k,2,2) = dvxc_rr(k,2,2) + e2*(vrrx(k,2) + vrrc(k,2) - &
                                                 vrzc(k,2)*(1.d0 + zeta)/rht)
         ENDIF
         !
         dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * (vsrx(k,1) + vsrc(k,1))
         dvxc_sr(k,1,2) = dvxc_sr(k,1,2) + e2 * vsrc(k,1)
         dvxc_sr(k,2,1) = dvxc_sr(k,2,1) + e2 * vsrc(k,2)
         dvxc_sr(k,2,2) = dvxc_sr(k,2,2) + e2 * (vsrx(k,2) + vsrc(k,2))
         !
         dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * (vssx(k,1) + vssc(k))
         dvxc_ss(k,1,2) = dvxc_ss(k,1,2) + e2 * vssc(k)
         dvxc_ss(k,2,1) = dvxc_ss(k,2,1) + e2 * vssc(k)
         dvxc_ss(k,2,2) = dvxc_ss(k,2,2) + e2 * (vssx(k,2) + vssc(k))
       ENDDO
       !
       !$acc end data
       DEALLOCATE( vrzc )
       !
    ENDIF
    !
    !$acc end data
    DEALLOCATE( vrrx, vsrx, vssx )
    DEALLOCATE( vrrc, vsrc, vssc )
    !
  ENDIF
  !
!# 354 "xc_wrapper_d_gga.f90"
  IF ( ALLOCATED(sigma) ) THEN
    !$acc exit data delete(sigma)
    DEALLOCATE( sigma )
  ENDIF
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dgcxc_
