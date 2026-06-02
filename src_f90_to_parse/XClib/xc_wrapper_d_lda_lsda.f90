!# 1 "xc_wrapper_d_lda_lsda.f90"
!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE dmxc( length, srd, rho_in, dmuxc, gpu_args_ )
  !----------------------------------------------------------------------
  !! Wrapper routine. Calls internal dmxc-driver routines or the external
  !! ones from Libxc, depending on the input choice.
  !
  USE kind_l,   ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: srd
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,srd,srd)
  !! the derivative of the xc potential
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
    !$acc data present( rho_in, dmuxc )
    CALL dmxc_( length, srd, rho_in, dmuxc )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho_in ), copyout( dmuxc )
    CALL dmxc_( length, srd, rho_in, dmuxc )
    !$acc end data
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
!------------------------------------------------------------------------
SUBROUTINE dmxc_( length, srd, rho_in, dmuxc )
  !----------------------------------------------------------------------
  !! Wrapper routine. Calls internal dmxc-driver routines or the external
  !! ones from Libxc, depending on the input choice.
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: iexch, icorr, is_libxc, rho_threshold_lda
  USE qe_drivers_d_lda_lsda
  !
!# 67 "xc_wrapper_d_lda_lsda.f90"
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: srd
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,srd,srd)
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
!# 91 "xc_wrapper_d_lda_lsda.f90"
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: ir, length_lxc, length_dlxc, fkind_x
  REAL(DP), PARAMETER :: small=1.E-10_DP, rho_trash=0.5_DP
  !
  !$acc data present( rho_in, dmuxc )
  !
  !$acc kernels
  dmuxc(:,:,:) = 0.0_DP
  !$acc end kernels
  !
  fkind_x = -1
  fkind_is_XC = .FALSE.
  !
!# 178 "xc_wrapper_d_lda_lsda.f90"
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. (.NOT.fkind_is_XC) ) THEN
    rho_threshold_lda = small
    IF ( srd == 1 ) CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
    IF ( srd == 2 ) CALL dmxc_lsda( length, rho_in, dmuxc )
    IF ( srd == 4 ) CALL dmxc_nc( length, rho_in, dmuxc )
  ENDIF
  !
!# 245 "xc_wrapper_d_lda_lsda.f90"
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dmxc_
