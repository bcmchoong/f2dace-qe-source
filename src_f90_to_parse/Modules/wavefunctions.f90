!# 1 "wavefunctions.f90"
!
! Copyright (C) 2002-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 9 "wavefunctions.f90"
!=----------------------------------------------------------------------------=!
   MODULE wavefunctions
!=----------------------------------------------------------------------------=!
     !! Wavefunction arrays.
     !
     USE kinds, ONLY :  DP
     USE control_flags, ONLY : use_gpu
!# 20 "wavefunctions.f90"
     IMPLICIT NONE
     SAVE
!# 23 "wavefunctions.f90"
     !
!# 27 "wavefunctions.f90"
     COMPLEX(DP), ALLOCATABLE, TARGET :: evc(:,:)
!# 29 "wavefunctions.f90"
       !! wavefunctions in the PW basis set.  
       !! noncolinear case: first index is a combined PW + spin index
       !
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic(:)
     !! additional memory for FFT
     COMPLEX(DP) , ALLOCATABLE, TARGET :: psic_nc(:,:)
     !! additional memory for FFT for the noncolinear case
     !
   CONTAINS
!# 39 "wavefunctions.f90"
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc(npwx, npol, nbnd)
  !----------------------------------------------------------------------------
  !! Dynamical allocation of wavefunctions.  
  !! Requires dimensions: \(\text{npwx}\), \(\text{nbnd}\), \(\text{npol}\)
  !
!# 48 "wavefunctions.f90"
  INTEGER, INTENT(IN) :: npwx, npol, nbnd
  !
  INTEGER :: istat
  !
  !
  ALLOCATE( evc(npwx*npol,nbnd) )
!civn: PIN evc memory here
!# 59 "wavefunctions.f90"
  !
END SUBROUTINE allocate_wfc
!
!----------------------------------------------------------------------------
SUBROUTINE deallocate_wfc()
  !----------------------------------------------------------------------------
!# 68 "wavefunctions.f90"
  !
  IMPLICIT NONE
  INTEGER :: istat
  !
!# 76 "wavefunctions.f90"
  IF( ALLOCATED( evc ) ) DEALLOCATE( evc )
  !
END SUBROUTINE deallocate_wfc
     !
!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions
!=----------------------------------------------------------------------------=!
