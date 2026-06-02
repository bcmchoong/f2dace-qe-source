!# 1 "fft_support.f90"
!# 2 "fft_support.f90"
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 10 "fft_support.f90"
!=----------------------------------------------------------------------=!
   MODULE fft_support
!=----------------------------------------------------------------------=!
!# 14 "fft_support.f90"
       USE fft_param
       IMPLICIT NONE
       SAVE
!# 18 "fft_support.f90"
       PRIVATE
       PUBLIC :: good_fft_dimension, allowed, good_fft_order
!# 21 "fft_support.f90"
!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!
!
!         FFT support Functions/Subroutines
!
!=----------------------------------------------------------------------=!
!
!
integer function good_fft_dimension (n)
  !
  ! Determines the optimal maximum dimensions of fft arrays
  ! Useful on some machines to avoid memory conflicts
  !
  IMPLICIT NONE
  INTEGER :: n, nx
  REAL(DP) :: log2n
  !
  ! this is the default: max dimension = fft dimension
  nx = n
  !
!# 53 "fft_support.f90"
  !
  good_fft_dimension = nx
  return
end function good_fft_dimension
!# 59 "fft_support.f90"
!=----------------------------------------------------------------------=!
!# 61 "fft_support.f90"
function allowed (nr)
!# 64 "fft_support.f90"
  ! find if the fft dimension is a good one
  ! a "bad one" is either not implemented (as on IBM with ESSL)
  ! or implemented but with awful performances (most other cases)
!# 68 "fft_support.f90"
  implicit none
  integer :: nr
!# 71 "fft_support.f90"
  logical :: allowed
  integer :: pwr (5)
  integer :: mr, i, fac, p, maxpwr
  integer :: factors( 5 ) = (/ 2, 3, 5, 7, 11 /)
!# 76 "fft_support.f90"
  ! find the factors of the fft dimension
!# 78 "fft_support.f90"
  mr  = nr
  pwr = 0
  factors_loop: do i = 1, 5
     fac = factors (i)
     maxpwr = NINT ( LOG( DBLE (mr) ) / LOG( DBLE (fac) ) ) + 1
     do p = 1, maxpwr
        if ( mr == 1 ) EXIT factors_loop
        if ( MOD (mr, fac) == 0 ) then
           mr = mr / fac
           pwr (i) = pwr (i) + 1
        endif
     enddo
  end do factors_loop
!# 92 "fft_support.f90"
  IF ( nr /= ( mr * 2**pwr (1) * 3**pwr (2) * 5**pwr (3) * 7**pwr (4) * 11**pwr (5) ) ) &
     CALL fftx_error__ (' allowed ', ' what ?!? ', 1 )
!# 95 "fft_support.f90"
  if ( mr /= 1 ) then
!# 97 "fft_support.f90"
     ! fft dimension contains factors > 11 : no good in any case
!# 99 "fft_support.f90"
     allowed = .false.
!# 101 "fft_support.f90"
  else
!# 113 "fft_support.f90"
     ! fftw and all other cases: no factors 7 and 11
!# 115 "fft_support.f90"
     allowed = ( ( pwr(4) == 0 ) .and. ( pwr(5) == 0 ) )
!# 119 "fft_support.f90"
  endif
!# 121 "fft_support.f90"
  return
end function allowed
!# 124 "fft_support.f90"
!=----------------------------------------------------------------------=!
!# 126 "fft_support.f90"
   INTEGER FUNCTION good_fft_order( nr, np )
!# 128 "fft_support.f90"
!
!    This function find a "good" fft order value greater or equal to "nr"
!
!    nr  (input) tentative order n of a fft
!
!    np  (optional input) if present restrict the search of the order
!        in the ensemble of multiples of np
!
!    Output: the same if n is a good number
!         the closest higher number that is good
!         an fft order is not good if not implemented (as on IBM with ESSL)
!         or implemented but with awful performances (most other cases)
!
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: nr
     INTEGER, OPTIONAL, INTENT(IN) :: np
     INTEGER :: new
!# 146 "fft_support.f90"
     new = nr
     IF( PRESENT( np ) ) THEN
       IF (np <= 0 .OR. np > nr) &
           CALL fftx_error__( ' good_fft_order ', ' invalid np ', 1 )
       DO WHILE( ( ( .NOT. allowed( new ) ) .OR. ( MOD( new, np ) /= 0 ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     ELSE
       DO WHILE( ( .NOT. allowed( new ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     END IF
!# 159 "fft_support.f90"
     IF( new > nfftx ) &
       CALL fftx_error__( ' good_fft_order ', ' fft order too large ', new )
!# 162 "fft_support.f90"
     good_fft_order = new
!# 164 "fft_support.f90"
     RETURN
   END FUNCTION good_fft_order
!# 168 "fft_support.f90"
!=----------------------------------------------------------------------=!
   END MODULE fft_support
!=----------------------------------------------------------------------=!
