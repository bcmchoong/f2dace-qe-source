!# 1 "fft_param.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE fft_param
  use iso_fortran_env, only : stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT
!# 18 "fft_param.f90"
  INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
  INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
  INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
!# 22 "fft_param.f90"
  
  INTEGER, PARAMETER :: ndims = 20
  !! Number of different FFT tables that the module
  !!could keep into memory without reinitialization
!# 27 "fft_param.f90"
  INTEGER, PARAMETER :: nfftx = 16385
  !!Max allowed fft dimension
!# 30 "fft_param.f90"
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
!# 32 "fft_param.f90"
  REAL(DP), PARAMETER :: eps8  = 1.0E-8_DP
!# 34 "fft_param.f90"
END MODULE fft_param
