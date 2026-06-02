!# 1 "random_numbers_gpu.f90"
!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE random_numbers_gpum
  !----------------------------------------------------------------------------
  !! Module for random numbers generation - GPU double.
  !
  USE kinds, ONLY : DP
!# 17 "random_numbers_gpu.f90"
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    FUNCTION randy_gpu ( irand )
      !------------------------------------------------------------------------
      REAL(DP) :: randy_gpu
      INTEGER, optional    :: irand
!# 30 "random_numbers_gpu.f90"
      call errore('randy','use randy_vect_gpu on GPUs',1)
    END FUNCTION randy_gpu
    !------------------------------------------------------------------------
    SUBROUTINE randy_vect_gpu ( r_d, n, irand )
      !------------------------------------------------------------------------
      !
      ! randy_vect_gpu(r, n, irand): reseed with initial seed idum=irand ( 0 <= irand <= ic, see below)
      !                     if randyv is not explicitly initialized, it will be
      !                     initialized with seed idum=0 the first time it is called
      ! randy_vect_gpu(r, n) : generate uniform real(DP) numbers x in [0,1]
      !
      USE random_numbers, ONLY : randy
!# 45 "random_numbers_gpu.f90"
      REAL(DP) :: r_d(n)
!# 49 "random_numbers_gpu.f90"
      INTEGER              :: i, n
      INTEGER, optional    :: irand
      !
      INTEGER              :: ist
      INTEGER, SAVE        :: idum=0
!# 76 "random_numbers_gpu.f90"
      ! randy_vect_gpu is not a GPU array in this case
      !
      ! ist means starting index here
      ist = 1
      IF ( present(irand) ) THEN
         r_d(1) = randy(irand)
         ist = 2
      END IF
      DO i = ist, n
         r_d(i) = randy()
      END DO
!# 88 "random_numbers_gpu.f90"
      RETURN
      !
    END SUBROUTINE randy_vect_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE randy_vect_debug_gpu (r_d, n, irand )
      !------------------------------------------------------------------------
      !
      ! randy_vect_debug_gpu(r, n, irand): reseed with initial seed idum=irand ( 0 <= irand <= ic, see below)
      !                           if randyv is not explicitly initialized, it will be
      !                           initialized with seed idum=0 the first time it is called
      ! randy_vect_debug_gpu(r, n) : generate uniform real(DP) numbers x in [0,1]
      !
      USE random_numbers, ONLY : randy
      !
      REAL(DP) :: r_d(n)
      INTEGER, optional    :: irand
!# 108 "random_numbers_gpu.f90"
      INTEGER :: n, i, ist
      REAL(DP), ALLOCATABLE :: aux_v(:)
      !
      ALLOCATE(aux_v(n))
      !
      ist = 1
      IF ( present(irand) ) THEN
         aux_v(1) = randy(irand)
         ist = 2
      END IF
      !
      DO i = ist, n
         aux_v(i) = randy()
      END DO
      !
      r_d(1:n) = aux_v(1:n)
      !
      DEALLOCATE(aux_v)
    END SUBROUTINE randy_vect_debug_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_random_seed ( )
      !------------------------------------------------------------------------
      !
      ! poor-man random seed for randy
      !
      INTEGER, DIMENSION (8) :: itime
      INTEGER  :: iseed
      REAL(DP) :: drand(1)
!# 140 "random_numbers_gpu.f90"
      !
      CALL date_and_time ( values = itime ) 
      ! itime contains: year, month, day, time difference in minutes, hours,
      !                 minutes, seconds and milliseconds. 
      iseed = ( itime(8) + itime(6) ) * ( itime(7) + itime(4) )
      CALL randy_vect_gpu ( drand, 1, iseed )
      CALL randy_vect_debug_gpu (drand, 1, iseed )
      !
    END SUBROUTINE set_random_seed
    !
END MODULE random_numbers_gpum
