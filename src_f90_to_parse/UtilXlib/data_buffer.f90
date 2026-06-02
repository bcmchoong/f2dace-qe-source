!# 1 "data_buffer.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE data_buffer
    USE util_param,  ONLY : DP
!# 14 "data_buffer.f90"
    !
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE, dimension(:)  :: mp_buff_r
    INTEGER, ALLOCATABLE, dimension(:)  :: mp_buff_i
    PUBLIC :: mp_buff_r, mp_buff_i
    !
!# 28 "data_buffer.f90"
END MODULE data_buffer
