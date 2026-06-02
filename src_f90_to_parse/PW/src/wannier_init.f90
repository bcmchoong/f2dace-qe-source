!# 1 "wannier_init.f90"
! Copyright (C) 2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 10 "wannier_init.f90"
!----------------------------------------------------------------------
SUBROUTINE wannier_init(hwwa)
  !----------------------------------------------------------------------
  !    
  ! ... This routine ALLOCATEs all dynamically ALLOCATEd arrays for wannier calc
  !
  USE wannier_new 
  USE wvfct, only : nbnd, npwx
  USE input_parameters, only: constrain_pot, wan_data
  USE lsda_mod, only: nspin
  USE ions_base, only : nat
  USE basis, only : natomwfc
  USE constants, only: rytoev
  USE klist, only: nks
  USE io_files
  USE buffers
  USE noncollin_module, ONLY : npol
!# 28 "wannier_init.f90"
  IMPLICIT NONE 
  
  LOGICAL,INTENT(IN) :: hwwa ! have we Wannier already?
  LOGICAL :: exst = .FALSE.,opnd
  INTEGER :: i, io_level
!# 34 "wannier_init.f90"
  ALLOCATE(pp(nwan,nbnd))
  ALLOCATE(wan_in(nwan,nspin))
  ALLOCATE(wannier_energy(nwan,nspin))
  ALLOCATE(wannier_occ(nwan,nwan,nspin))
  ALLOCATE(coef(natomwfc,nwan,nspin))
  
  coef = (0.d0,0.d0)
  wannier_energy = (0.d0,0.d0)
  wannier_occ = (0.d0,0.d0)
!# 44 "wannier_init.f90"
  wan_in(1:nwan,1:nspin) = wan_data(1:nwan,1:nspin)
  
  IF(.NOT. hwwa) THEN
!# 48 "wannier_init.f90"
     IF(use_energy_int) THEN
        do i=1,nwan
           wan_in(i,:)%bands_from = (1.d0/rytoev)*wan_in(i,:)%bands_from
           wan_in(i,:)%bands_to = (1.d0/rytoev)*wan_in(i,:)%bands_to
        end do
     END IF
     
     CALL wannier_check()
  end if
!# 58 "wannier_init.f90"
  ALLOCATE(wan_pot(nwan,nspin))
  wan_pot(1:nwan,1:nspin) = constrain_pot(1:nwan,1:nspin)
  
  !now open files to store projectors and wannier functions
  nwordwpp = nwan*nbnd*npol
  nwordwf = nwan*npwx*npol
  io_level = 1
  CALL open_buffer( iunwpp, 'wproj', nwordwpp, io_level, exst )
  CALL open_buffer( iunwf, 'wwf', nwordwf, io_level, exst )
!# 68 "wannier_init.f90"
  ! For atomic wavefunctions
  nwordatwfc = npwx*natomwfc*npol
  INQUIRE( UNIT = iunsat, OPENED = opnd )
  IF(.NOT. opnd) CALL open_buffer( iunsat,'satwfc',nwordatwfc,io_level,exst )
!# 73 "wannier_init.f90"
  RETURN
  !
END SUBROUTINE wannier_init
