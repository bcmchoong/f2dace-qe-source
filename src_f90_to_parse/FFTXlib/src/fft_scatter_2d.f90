!# 1 "fft_scatter_2d.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni, modified by Paolo Giannozzi
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_scatter_2d
!=----------------------------------------------------------------------=!
!# 18 "fft_scatter_2d.f90"
        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param
!# 21 "fft_scatter_2d.f90"
        IMPLICIT NONE
!# 23 "fft_scatter_2d.f90"
        SAVE
!# 25 "fft_scatter_2d.f90"
        PRIVATE
!# 27 "fft_scatter_2d.f90"
        PUBLIC :: fft_scatter
!# 29 "fft_scatter_2d.f90"
!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!   ALLTOALL based SCATTER, should be better on network
!   with a defined topology, like on bluegene and cray machine
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn )
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (isgn > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nr3x elements for a fft of order nr3
  !    nr3x can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !    f_in  contains input columns, is destroyed on output
  !    f_aux contains output planes
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains input planes, is destroyed on output
  !    f_in  contains output columns
  !
  IMPLICIT NONE
!# 69 "fft_scatter_2d.f90"
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
!# 259 "fft_scatter_2d.f90"
  RETURN
!# 261 "fft_scatter_2d.f90"
END SUBROUTINE fft_scatter
!
!=----------------------------------------------------------------------=!
   END MODULE fft_scatter_2d
!=----------------------------------------------------------------------=!
