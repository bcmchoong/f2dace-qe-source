!# 1 "xc_wrapper_mgga.f90"
!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                       v2c, v3c, gpu_args_ )
  !----------------------------------------------------------------------------------
  !! Wrapper to gpu or non gpu version of \(\texttt{xc_metagcx}\).
  !
  USE kind_l,               ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER, INTENT(IN) :: ns
  !! spin components
  INTEGER, INTENT(IN) :: np
  !! first dimension of v2c
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! the charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! grho = \nabla rho
  REAL(DP), INTENT(IN) :: tau(length,ns)
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: ex(length)
  !! sx = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: ec(length)
  !! sc = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x(length,ns)
  !! v3x = D(E_x)/D(tau)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c(np,length,ns)
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c(length,ns)
  !! v3c = D(E_c)/D(tau)
  LOGICAL, INTENT(IN), OPTIONAL :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  !
  gpu_args = .FALSE.
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
    CALL xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                      v2c, v3c )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho, grho, tau ), copyout( ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
    CALL xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, &
                      v2c, v3c )
    !$acc end data
    !
  ENDIF  
  !
  RETURN
  !
END SUBROUTINE
!
!
!----------------------------------------------------------------------------------------
SUBROUTINE xc_metagcx_( length, ns, np, rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-------------------------------------------------------------------------------------
  !! Wrapper routine. Calls internal metaGGA drivers or the Libxc ones,
  !! depending on the input choice.
  !
!# 87 "xc_wrapper_mgga.f90"
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: imeta, imetac, is_libxc, rho_threshold_mgga,&
                                  grho2_threshold_mgga, tau_threshold_mgga,   &
                                  ishybrid, exx_started, exx_fraction
  USE qe_drivers_mgga
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER, INTENT(IN) :: ns
  !! spin components
  INTEGER, INTENT(IN) :: np
  !! first dimension of v2c
  REAL(DP), INTENT(IN) :: rho(length,ns)
  !! the charge density
  REAL(DP), INTENT(IN) :: grho(3,length,ns)
  !! grho = \nabla rho
  REAL(DP), INTENT(IN) :: tau(length,ns)
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: ex(length)
  !! sx = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: ec(length)
  !! sc = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x(length,ns)
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x(length,ns)
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x(length,ns)
  !! v3x = D(E_x)/D(tau)
  REAL(DP), INTENT(OUT) :: v1c(length,ns)
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c(np,length,ns)
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c(length,ns)
  !! v3c = D(E_c)/D(tau)
  !
  ! ... local variables
  !
  INTEGER :: k, is, ipol
  REAL(DP), ALLOCATABLE :: grho2(:,:)
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
!# 145 "xc_wrapper_mgga.f90"
  !
  !$acc data present( rho, grho, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !
!# 220 "xc_wrapper_mgga.f90"
  !
  IF ( .NOT.is_libxc(5) .AND. imetac==0 ) THEN
    IF (ns == 1) THEN
      !
      ALLOCATE( grho2(length,ns) )
      !$acc data create( grho2 )
      !
      !$acc parallel loop
      DO k = 1, length
        grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
      ENDDO
      !
      CALL tau_xc( length, rho(:,1), grho2, tau(:,1), ex, ec, v1x(:,1), &
                   v2x(:,1), v3x(:,1), v1c(:,1), v2c, v3c(:,1) )
      !
      !$acc end data
      DEALLOCATE( grho2 )
      !
    ELSEIF (ns == 2) THEN
      !
      CALL tau_xc_spin( length, rho, grho, tau, ex, ec, v1x, v2x, v3x, &
                        v1c, v2c, v3c )
      !
    ENDIF
  ENDIF
  !
!# 387 "xc_wrapper_mgga.f90"
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE xc_metagcx_
