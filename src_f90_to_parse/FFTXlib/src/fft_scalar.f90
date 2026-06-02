!# 1 "fft_scalar.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI
! (both 3d for serial execution and 1d+2d FFTs for parallel execution);    !
! legacy NEC ASL libraries (3d only, no parallel execution)                !
! CUDA FFT for NVidiia GPUs
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood                !
! Last update Feb 2021
!--------------------------------------------------------------------------!
!# 21 "fft_scalar.f90"
!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!
!# 25 "fft_scalar.f90"
     USE fft_param
!# 35 "fft_scalar.f90"
     USE fft_scalar_fftw
!# 42 "fft_scalar.f90"
     IMPLICIT NONE
     SAVE
!# 45 "fft_scalar.f90"
     PRIVATE
     PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds
!# 51 "fft_scalar.f90"
   END MODULE fft_scalar
