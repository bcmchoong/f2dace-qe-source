!# 1 "xc_wrapper_lda_lsda.f90"
!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------------
SUBROUTINE xc( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out, gpu_args_ )
  !--------------------------------------------------------------------------------------
  !! Wrapper routine to \(\texttt{xc_}\) or \(\texttt{xc_gpu}\).
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: srd
  !! spin dimension of rho
  INTEGER, INTENT(IN) :: svd
  !! spin dimension of v
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! Charge density
  REAL(DP), INTENT(OUT) :: ex_out(length)
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vx_out(length,svd)
  !! \(dE_x(\text{rho})/d\text{rho}  ( NOT d\epsilon_x(\text{rho})/d\text{rho} )
  REAL(DP), INTENT(OUT) :: ec_out(length)
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vc_out(length,svd)
  !! \(dE_c(\text{rho})/d\text{rho}  ( NOT d\epsilon_c(\text{rho})/d\text{rho} )
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
    !$acc data present( rho_in, ex_out, ec_out, vx_out, vc_out )
    CALL xc_( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho_in ), copyout( ex_out, ec_out, vx_out, vc_out )
    CALL xc_( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
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
SUBROUTINE xc_( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
  !-------------------------------------------------------------------------
  !! Wrapper xc LDA - openACC version.  
  !! See comments in routine \(\texttt{xc}\) for variable explanations.
  !
!# 71 "xc_wrapper_lda_lsda.f90"
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: iexch, icorr, is_libxc, rho_threshold_lda, &
                                  finite_size_cell_volume_set
  USE qe_drivers_lda_lsda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: srd, svd
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  REAL(DP), INTENT(OUT) :: ex_out(length), vx_out(length,svd)
  REAL(DP), INTENT(OUT) :: ec_out(length), vc_out(length,svd)
  !
  ! ... local variables
  !
  LOGICAL :: is_libxc1, is_libxc2
!# 100 "xc_wrapper_lda_lsda.f90"
  REAL(DP), ALLOCATABLE :: zeta(:)
  REAL(DP) :: arho_ir
  INTEGER  :: ir
  !
  !$acc data present( rho_in, ex_out, ec_out, vx_out, vc_out )
  !
  is_libxc1 = is_libxc(1)
  is_libxc2 = is_libxc(2)
  !
!# 183 "xc_wrapper_lda_lsda.f90"
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) ) THEN
     !
     SELECT CASE( srd )
     CASE( 1 )
        !
        IF ((iexch==8 .AND. .NOT.is_libxc1) .OR. (icorr==10 .AND. .NOT.is_libxc2)) THEN
          IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
              'finite size corrected exchange used w/o initialization', 1 )
        ENDIF
        CALL xc_lda( length, rho_in(:,1), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
        !
     CASE( 2 )
        !
        ALLOCATE( zeta(length) )
        !$acc data create( zeta )
        !$acc parallel loop
        DO ir = 1, length
          arho_ir = ABS(rho_in(ir,1))
          IF (arho_ir > rho_threshold_lda) zeta(ir) = rho_in(ir,2) / arho_ir
        ENDDO
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end data
        DEALLOCATE( zeta )
        !
     CASE( 4 )
        !
        ALLOCATE( zeta(length) )
        !$acc data create( zeta )
        !$acc parallel loop
        DO ir = 1, length
          arho_ir = ABS( rho_in(ir,1) )
          IF (arho_ir > rho_threshold_lda) zeta(ir) = SQRT( rho_in(ir,2)**2 + rho_in(ir,3)**2 + &
                                                          rho_in(ir,4)**2 ) / arho_ir ! amag/arho
        ENDDO
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end data
        DEALLOCATE( zeta )
        !
     CASE DEFAULT
        !
        CALL xclib_error( 'xc_LDA', 'Wrong ns input', 2 )
        !
     END SELECT
     !
  ENDIF
  !
  !  ... fill output arrays
  !
!# 267 "xc_wrapper_lda_lsda.f90"
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE xc_
