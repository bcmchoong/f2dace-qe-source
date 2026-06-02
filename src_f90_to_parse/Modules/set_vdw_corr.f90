!# 1 "set_vdw_corr.f90"
!# 2 "set_vdw_corr.f90"
! Copyright (C) 2019 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw, mbd_vdw, lxdm )
  USE io_global, ONLY: stdout
  !
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(in) :: vdw_corr
  LOGICAL, INTENT(out) :: llondon, ldftd3, ts_vdw, mbd_vdw, lxdm
  !
!# 17 "set_vdw_corr.f90"
  llondon= .FALSE.
  ldftd3 = .FALSE.
  ts_vdw = .FALSE.
  mbd_vdw= .FALSE.
  lxdm   = .FALSE.
!# 23 "set_vdw_corr.f90"
  SELECT CASE( TRIM( vdw_corr ) )
  CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
     llondon= .TRUE.
!# 27 "set_vdw_corr.f90"
  CASE( 'grimme-d3', 'Grimme-D3', 'DFT-D3', 'dft-d3' )
     ldftd3 = .TRUE.
!# 30 "set_vdw_corr.f90"
  CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
     ts_vdw = .TRUE.
!# 33 "set_vdw_corr.f90"
  CASE( 'MBD', 'mbd', 'many-body-dispersion', 'mbd_vdw' )
     ts_vdw = .TRUE.
     mbd_vdw = .TRUE.
!# 37 "set_vdw_corr.f90"
  CASE( 'XDM', 'xdm' )
     lxdm   = .TRUE.
!# 40 "set_vdw_corr.f90"
  CASE('none','')
!# 42 "set_vdw_corr.f90"
  CASE DEFAULT
     WRITE (stdout,*)
     CALL infomsg('set_vdw_corr','WARNING: unknown vdw correction (vdw_corr): '//TRIM(vdw_corr)//'. No vdw correction used.')
     WRITE (stdout,*)
!# 47 "set_vdw_corr.f90"
  END SELECT
    
!# 50 "set_vdw_corr.f90"
END SUBROUTINE set_vdw_corr
