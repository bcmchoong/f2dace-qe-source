!# 1 "fft_scatter.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! Basic scatter operations needed by parallel FFT
! Written by Carlo Cavazzoni and Stefano de Gironcoli, modified by PG
!
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_scatter
!=----------------------------------------------------------------------=!
!# 19 "fft_scatter.f90"
        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param
!# 22 "fft_scatter.f90"
        IMPLICIT NONE
!# 24 "fft_scatter.f90"
        SAVE
!# 26 "fft_scatter.f90"
        PRIVATE
!# 28 "fft_scatter.f90"
        PUBLIC :: fft_scatter_xy, fft_scatter_yz, fft_scatter_many_xy, fft_scatter_many_yz
        PUBLIC :: fft_scatter_tg, fft_scatter_tg_opt
!# 31 "fft_scatter.f90"
!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_xy ( desc, f_in, f_aux, nxx_, isgn, comm )
  !-----------------------------------------------------------------------
  !
  ! Transpose of the fft xy planes across the comm communicator.
  ! If the optional comm is not provided as input, the transpose is made
  ! across desc%comm2 communicator.
  !
  ! a) From Y-oriented columns to X-oriented partial slices (isgn > 0)
  !    Active columns along the Y direction corresponding to a subset of the
  !    active X values and a range of Z values (in this order) are stored
  !    consecutively for each processor and are such that the subgroup owns
  !    all data for a range of Z values.
  !
  !    The Y pencil -> X-oriented partial slices transposition is performed
  !    in the subgroup of processors (desc%comm2) owning this range of Z values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Y
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc2" are desc%nr2p(iproc2)*desc%nr1p(me2)*desc%my_nr3p big.
  !    2) all processors communicate to exchange slices (all sectin of columns with
  !       Y in the slice belonging to "me" must be received, all the others
  !       must be sent to "iproc2")
  !
  !    Finally one gets the "partial slice" representation: each processor has
  !    all the X values of desc%my_nr2p Y and desc%my_nr3p Z values.
  !    Data are organized with the X index running fastest, then Y, then Z.
  !
  !    f_in  contains the input Y columns, is destroyed on output
  !    f_aux contains the output X-oriented partial slices.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input X-oriented partial slices, is destroyed on output
  !    f_in  contains the output Y columns.
  !
  IMPLICIT NONE
!# 76 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  INTEGER, OPTIONAL, INTENT(in) :: comm
  INTEGER :: nr1_temp(1), comm_
!# 295 "fft_scatter.f90"
END SUBROUTINE fft_scatter_xy
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_many_xy ( desc, f_in, f_aux, isgn, howmany)
  !-----------------------------------------------------------------------
  !
  ! transpose of the fft xy planes across the desc%comm2 communicator
  !
  ! a) From Y-oriented columns to X-oriented partial slices (isgn > 0)
  !    Active columns along the Y direction corresponding to a subset of the
  !    active X values and a range of Z values (in this order) are stored
  !    consecutively for each processor and are such that the subgroup owns
  !    all data for a range of Z values.
  !
  !    The Y pencil -> X-oriented partial slices transposition is performed
  !    in the subgroup of processors (desc%comm2) owning this range of Z values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Y
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc2" are desc%nr2p(iproc2)*desc%nr1p(me2)*desc%my_nr3p big.
  !    2) all processors communicate to exchange slices (all sectin of columns with
  !       Y in the slice belonging to "me" must be received, all the others
  !       must be sent to "iproc2")
  !
  !    Finally one gets the "partial slice" representation: each processor has
  !    all the X values of desc%my_nr2p Y and desc%my_nr3p Z values.
  !    Data are organized with the X index running fastest, then Y, then Z.
  !
  !    f_in  contains the input Y columns, is destroyed on output
  !    f_aux contains the output X-oriented partial slices.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input X-oriented partial slices, is destroyed on output
  !    f_in  contains the output Y columns.
  !
  IMPLICIT NONE
!# 335 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in), TARGET :: desc
  INTEGER, INTENT(in)                            :: isgn
  INTEGER, INTENT(in)                            :: howmany
!# 339 "fft_scatter.f90"
  COMPLEX (DP), INTENT(inout)   :: f_in (:), f_aux (:)
!# 575 "fft_scatter.f90"
  RETURN
98 format ( 10 ('(',2f12.9,')') )
99 format ( 20 ('(',2f12.9,')') )
!# 579 "fft_scatter.f90"
END SUBROUTINE fft_scatter_many_xy
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_yz ( desc, f_in, f_aux, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! transpose of the fft yz planes across the desc%comm3 communicator
  !
  ! a) From Z-oriented columns to Y-oriented colums (isgn > 0)
  !    Active columns (or sticks or pencils) along the Z direction for each
  !    processor are stored consecutively and are such that they correspond
  !    to a subset of the active X values.
  !
  !    The pencil -> slices transposition is performed in the subgroup
  !    of processors (desc%comm3) owning these X values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Z
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc3" are desc%nr3p(iproc3)*desc%nsw/nsp(me) big.
  !    2) all processors communicate to exchange slices (all columns with
  !       Z in the slice belonging to "me" must be received, all the others
  !       must be sent to "iproc3")
  !
  !    Finally one gets the "slice" representation: each processor has
  !    desc%nr3p(mype3) Z values of all the active pencils along Y for the
  !    X values of the current group. Data are organized with the Y index
  !    running fastest, then the reordered X values, then Z.
  !
  !    f_in  contains the input Z columns, is destroyed on output
  !    f_aux contains the output Y colums.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input Y columns, is destroyed on output
  !    f_in  contains the output Z columns.
  !
  IMPLICIT NONE
!# 619 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
!# 841 "fft_scatter.f90"
END SUBROUTINE fft_scatter_yz
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_many_yz ( desc, f_in, f_aux, isgn, howmany )
  !-----------------------------------------------------------------------
  !
  ! transpose of the fft yz planes across the desc%comm3 communicator
  !
  ! a) From Z-oriented columns to Y-oriented colums (isgn > 0)
  !    Active columns (or sticks or pencils) along the Z direction for each
  !    processor are stored consecutively and are such that they correspond
  !    to a subset of the active X values.
  !
  !    The pencil -> slices transposition is performed in the subgroup
  !    of processors (desc%comm3) owning these X values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Z
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc3" are desc%nr3p(iproc3)*desc%nsw/nsp(me) big.
  !    2) all processors communicate to exchange slices (all columns with
  !       Z in the slice belonging to "me" must be received, all the others
  !       must be sent to "iproc3")
  !
  !    Finally one gets the "slice" representation: each processor has
  !    desc%nr3p(mype3) Z values of all the active pencils along Y for the
  !    X values of the current group. Data are organized with the Y index
  !    running fastest, then the reordered X values, then Z.
  !
  !    f_in  contains the input Z columns, is destroyed on output
  !    f_aux contains the output Y colums.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input Y columns, is destroyed on output
  !    f_in  contains the output Z columns.
  !
  IMPLICIT NONE
!# 881 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in), TARGET :: desc
  COMPLEX (DP), INTENT(inout)                    :: f_in(:), f_aux(:)
  INTEGER, INTENT(in)                            :: isgn, howmany
!# 1110 "fft_scatter.f90"
  RETURN
98 format ( 10 ('(',2f12.9,')') )
99 format ( 20 ('(',2f12.9,')') )
!# 1114 "fft_scatter.f90"
END SUBROUTINE fft_scatter_many_yz
!# 1116 "fft_scatter.f90"
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg ( desc, f_in, f_aux, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! task group wavefunction redistribution
  !
  ! a) (isgn >0 ) From many-wfc partial-plane arrangement to single-wfc whole-plane one
  !
  ! b) (isgn <0 ) From single-wfc whole-plane arrangement to many-wfc partial-plane one
  !
  ! in both cases:
  !    f_in  contains the input data, is overwritten with the desired output
  !    f_aux is used as working array, may contain garbage in output
  !
  IMPLICIT NONE
!# 1132 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
!# 1136 "fft_scatter.f90"
  INTEGER :: ierr
!# 1138 "fft_scatter.f90"
  CALL start_clock ('fft_scatt_tg')
!# 1140 "fft_scatter.f90"
  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )
!# 1159 "fft_scatter.f90"
  CALL stop_clock ('fft_scatt_tg')
!# 1161 "fft_scatter.f90"
  RETURN
!# 1163 "fft_scatter.f90"
END SUBROUTINE fft_scatter_tg
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg_opt ( desc, f_in, f_out, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! task group wavefunction redistribution
  !
  ! a) (isgn >0 ) From many-wfc partial-plane arrangement to single-wfc whole-plane one
  !
  ! b) (isgn <0 ) From single-wfc whole-plane arrangement to many-wfc partial-plane one
  !
  ! in both cases:
  !    f_in  contains the input data
  !    f_out contains the output data
  !
  IMPLICIT NONE
!# 1181 "fft_scatter.f90"
  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_out (nxx_)
!# 1185 "fft_scatter.f90"
  INTEGER :: ierr
!# 1187 "fft_scatter.f90"
  CALL start_clock ('fft_scatt_tg')
!# 1189 "fft_scatter.f90"
  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )
!# 1207 "fft_scatter.f90"
  CALL stop_clock ('fft_scatt_tg')
!# 1209 "fft_scatter.f90"
  RETURN
!# 1211 "fft_scatter.f90"
END SUBROUTINE fft_scatter_tg_opt
!
!=----------------------------------------------------------------------=!
   END MODULE fft_scatter
!=----------------------------------------------------------------------=!
