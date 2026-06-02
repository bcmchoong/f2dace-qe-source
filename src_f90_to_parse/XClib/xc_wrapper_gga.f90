!# 1 "xc_wrapper_gga.f90"
!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud, &
                   gpu_args_ )
  !-------------------------------------------------------------------------
  !! Wrapper to gpu or non gpu version of \(\texttt{xc_gcx}\).
  !
  USE kind_l,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: ns
  !! spin dimension for input
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! Charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! gradient
  REAL(DP), INTENT(OUT) :: ex(length)
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec(length)
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! exchange potential (density part)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! exchange potential (gradient part)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! correlation potential (density part)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  !! correlation potential (gradient part)
  REAL(DP), INTENT(OUT), OPTIONAL :: v2c_ud(length)
  !! correlation potential, cross term
  LOGICAL, INTENT(IN), OPTIONAL :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  REAL(DP), ALLOCATABLE :: v2c_dummy(:)
  !
  gpu_args = .FALSE.
  !
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF (ns==2 .AND. .NOT. PRESENT(v2c_ud)) CALL xclib_infomsg( 'xc_gcx', 'WARNING: cross &
                &term v2c_ud not found xc_gcx (gga) call with polarized case' )
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( rho, grho, ex, ec, v1x, v2x, v1c, v2c )
    IF (PRESENT(v2c_ud)) THEN
      !$acc data present( v2c_ud )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
      !$acc end data
    ELSE
      ALLOCATE( v2c_dummy(length) )
      !$acc data create( v2c_dummy )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      !$acc end data
      DEALLOCATE( v2c_dummy )
    ENDIF
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho, grho ), copyout( ex, ec, v1x, v2x, v1c, v2c )
    IF (PRESENT(v2c_ud)) THEN
      !$acc data copyout( v2c_ud )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
      !$acc end data
    ELSE
      ALLOCATE( v2c_dummy(length) )
      !$acc data create( v2c_dummy )
      CALL xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_dummy )
      !$acc end data
      DEALLOCATE( v2c_dummy )
    ENDIF
    !$acc end data
    !
  ENDIF  
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gcx_( length, ns, rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !-------------------------------------------------------------------------
  !! GGA wrapper routine - gpu double.
  !
!# 102 "xc_wrapper_gga.f90"
  !
  USE kind_l,               ONLY: DP
  USE xclib_utils_and_para, ONLY: error_msg, nowarning
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda,   &
                                  ishybrid, exx_started, exx_fraction
  USE qe_drivers_gga
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: ns
  REAL(DP), INTENT(IN) :: rho(length,ns)
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  REAL(DP), INTENT(OUT) :: ex(length)
  REAL(DP), INTENT(OUT) :: ec(length)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c(length,ns)
  REAL(DP), INTENT(OUT) :: v2c_ud(length)
  !
  ! ... local variables
  !
!# 144 "xc_wrapper_gga.f90"
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: k, is, ierr, fkind_x
  REAL(DP) :: rho_up, rho_dw, grho_up, grho_dw, sgn1
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  !$acc data present( rho, grho, ex, ec, v1x, v2x, v1c, v2c, v2c_ud )
  !
  ierr = 0
  fkind_x = -1
  fkind_is_XC = .FALSE.
  !
!# 212 "xc_wrapper_gga.f90"
  !
  IF (ANY(.NOT.is_libxc(3:4))) THEN
    !
    ALLOCATE( rh(length), grho2(length,ns) )
    !$acc enter data create( rh, grho2 )
    !$acc parallel loop
    DO k = 1, length
       rh(k) = ABS(rho(k,1))
       grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
    ENDDO
    !
    IF ( ns==1 ) THEN
      CALL gcxc( length, rh, grho2(:,1), ex, ec, v1x(:,1), v2x(:,1), v1c(:,1), &
                 v2c(:,1), ierr )
      !
      !$acc parallel loop
      DO k = 1, length
         sgn1 = SIGN(1._DP, rho(k,1))
         ex(k) = ex(k) * sgn1
         ec(k) = ec(k) * sgn1
      ENDDO
    ENDIF
    !
  ENDIF
  !
  ! ---- GGA CORRELATION
  !
  IF ( is_libxc(4) ) THEN  !lda part of LYP not present in libxc (still so? - check)
    !
!# 304 "xc_wrapper_gga.f90"
    !
  ELSEIF ( (.NOT.is_libxc(4)) .AND. (.NOT.fkind_is_XC) ) THEN
    !
    IF ( ns /= 1 ) THEN
       !
       IF (igcc==3 .OR. igcc==7 .OR. igcc==13 ) THEN
          !
          ALLOCATE( grho_ud(length) )
          !$acc data create(grho_ud)
          !
          !$acc parallel loop
          DO k = 1, length
            grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
            grho_ud(k) = grho(1,k,1) * grho(1,k,2) + grho(2,k,1) * grho(2,k,2) + &
                         grho(3,k,1) * grho(3,k,2)
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
          ENDDO
          !
          CALL gcc_spin_more( length, rho, grho2, grho_ud, ec, v1c, v2c, v2c_ud )
          !
          !$acc end data
          DEALLOCATE( grho_ud )
          !
       ELSE
          !
          ALLOCATE( zeta(length) )
          !$acc data create( zeta )
          !
          !$acc parallel loop
          DO k = 1, length
            rh(k) = rho(k,1) + rho(k,2)
            IF ( rh(k) > rho_threshold_gga ) THEN
              zeta(k) = (rho(k,1)-rho(k,2)) / rh(k)
            ELSE
              zeta(k) = 2.0_DP ! trash value, gcc-routines get rid of it when present
            ENDIF
            grho2(k,1) = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                         ( grho(2,k,1) + grho(2,k,2) )**2 + &
                         ( grho(3,k,1) + grho(3,k,2) )**2
            grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
          ENDDO
          !
          CALL gcc_spin( length, rh, zeta, grho2(:,1), ec, v1c, v2c(:,1) )
          !
          !$acc parallel loop
          DO k = 1, length
            v2c(k,2) = v2c(k,1)
            IF ( ns==2 ) v2c_ud(k) = v2c(k,1)
          ENDDO
          !
          !$acc end data
          DEALLOCATE( zeta )
          !
       ENDIF
       !
    ENDIF
    !
  ENDIF
  !
  ! --- GGA EXCHANGE
  !
  IF ( is_libxc(3) ) THEN
    !
!# 427 "xc_wrapper_gga.f90"
    !
  ELSE
    !
    IF ( ns > 1 ) THEN
      !$acc parallel loop collapse(2)
      DO is = 1, ns
        DO k = 1, length
          grho2(k,is) = grho(1,k,is)**2 + grho(2,k,is)**2 + grho(3,k,is)**2
        ENDDO
      ENDDO
      !
      CALL gcx_spin( length, rho, grho2, ex, v1x, v2x, ierr )
    ENDIF
    !
  ENDIF
  !
  IF (ANY(.NOT.is_libxc(3:4))) THEN
    !$acc exit data delete( rh, grho2 )
    DEALLOCATE( rh, grho2 )
  ENDIF
!# 453 "xc_wrapper_gga.f90"
  !
  !$acc end data
  !
  IF (ierr/=0 .AND. .NOT.nowarning) CALL xclib_error( 'xc_gcx_', error_msg(ierr), 1 )
  !
  RETURN
  !
END SUBROUTINE xc_gcx_
