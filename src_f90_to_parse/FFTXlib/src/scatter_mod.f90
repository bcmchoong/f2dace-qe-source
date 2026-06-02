!# 1 "scatter_mod.f90"
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
! Rewritten by Stefano de Gironcoli
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE scatter_mod
!=----------------------------------------------------------------------=!
!# 19 "scatter_mod.f90"
        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param
!# 22 "scatter_mod.f90"
        IMPLICIT NONE
!# 24 "scatter_mod.f90"
        INTERFACE gather_grid
           MODULE PROCEDURE gather_real_grid, gather_complex_grid
        END INTERFACE
!# 28 "scatter_mod.f90"
        INTERFACE scatter_grid
           MODULE PROCEDURE scatter_real_grid, scatter_complex_grid
        END INTERFACE
!# 32 "scatter_mod.f90"
        SAVE
!# 34 "scatter_mod.f90"
        PRIVATE
!# 36 "scatter_mod.f90"
        PUBLIC :: gather_grid, scatter_grid
        PUBLIC :: cgather_sym, cgather_sym_many, cscatter_sym_many
!# 39 "scatter_mod.f90"
!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!----------------------------------------------------------------------------
SUBROUTINE gather_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - version for real arrays
  !
  ! ... REAL*8  f_in  = distributed variable (dfft%nnr)
  ! ... REAL*8  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  !
!# 119 "scatter_mod.f90"
  CALL fftx_error__(' gather_real_grid', 'do not use in serial execution', 1)
!# 121 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE gather_real_grid
!# 126 "scatter_mod.f90"
!----------------------------------------------------------------------------
SUBROUTINE gather_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - complex arrays
  !
  ! ... COMPLEX*16  f_in  = distributed variable (dfft%nnr)
  ! ... COMPLEX*16  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
!# 204 "scatter_mod.f90"
  CALL fftx_error__('gather_complex_grid', 'do not use in serial execution', 1)
!# 206 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE gather_complex_grid
!# 211 "scatter_mod.f90"
!----------------------------------------------------------------------------
SUBROUTINE scatter_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... REAL*8  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... REAL*8  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  REAL(DP), ALLOCATABLE ::  f_aux(:)
  !
!# 282 "scatter_mod.f90"
  CALL fftx_error__('scatter_real_grid', 'do not use in serial execution', 1)
!# 284 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE scatter_real_grid
!----------------------------------------------------------------------------
SUBROUTINE scatter_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... COMPLEX*16  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... COMPLEX*16  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
!# 367 "scatter_mod.f90"
  CALL fftx_error__('scatter_complex_grid', 'do not use in serial execution', 1)
!# 369 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE scatter_complex_grid
!
! ... "gather"-like subroutines
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( dfftp, f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (used in phonon code)
  ! ... Differs from gather_grid because mpi_allgatherv is used instead
  ! ... of mpi_gatherv - all data is gathered on ALL processors
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  COMPLEX(DP) :: f_in( : ), f_out(:)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
!# 441 "scatter_mod.f90"
  CALL fftx_error__('cgather_sym', 'do not use in serial execution', 1)
!# 443 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym_many( dfftp, f_in, f_out, nbnd, nbnd_proc, start_nbnd_proc )
  !-----------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... This routine generalizes cgather_sym, receiveng nbnd complex
  ! ... distributed functions and collecting nbnd_proc(dfftp%mype+1)
  ! ... functions in each processor.
  ! ... start_nbnd_proc(dfftp%mype+1), says where the data for each processor
  ! ... start in the distributed variable
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx,nbnd)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x,nbnd_proc(dfftp%mype+1))
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nnr,nbnd)
  COMPLEX(DP) :: f_out(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
!# 576 "scatter_mod.f90"
  CALL fftx_error__('cgather_sym_many', 'do not use in serial execution', 1)
!# 578 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE cgather_sym_many
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym_many( dfftp, f_in, f_out, target_ibnd, nbnd, nbnd_proc, &
                              start_nbnd_proc   )
  !----------------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... generalizes cscatter_sym. It assumes that each processor has
  ! ... a certain number of bands (nbnd_proc(dfftp%mype+1)). The processor
  ! ... that has target_ibnd scatters it to all the other processors
  ! ... that receive a distributed part of the target function.
  ! ... start_nbnd_proc(dfftp%mype+1) is used to identify the processor
  ! ... that has the required band
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1x*nr2x*nr3x, nbnd_proc(dfftp%mype+1) )
  ! ... COMPLEX*16  f_out = distributed variable (nrxx)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP) :: f_out(dfftp%nnr)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  INTEGER :: target_ibnd
  !
!# 667 "scatter_mod.f90"
  CALL fftx_error__('cscatter_sym_many', 'do not use in serial execution', 1)
!# 669 "scatter_mod.f90"
  !
  RETURN
  !
END SUBROUTINE cscatter_sym_many
!# 674 "scatter_mod.f90"
!=----------------------------------------------------------------------=!
   END MODULE scatter_mod
!=----------------------------------------------------------------------=!
!
