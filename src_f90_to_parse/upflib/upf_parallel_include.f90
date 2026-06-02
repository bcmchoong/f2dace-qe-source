!# 1 "upf_parallel_include.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE upf_parallel_include
!# 20 "upf_parallel_include.f90"
        ! dummy world and null communicator
        INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
        INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
        INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
!# 25 "upf_parallel_include.f90"
END MODULE upf_parallel_include
