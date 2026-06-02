!# 1 "fft_parallel_2d.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni
!     Last update Apr. 2009
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel_2d
!
!# 25 "fft_parallel_2d.f90"
   !
   USE fft_param
   IMPLICIT NONE
   SAVE
   !
!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE fft_scatter_2d,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : )  ! array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  !
  !
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  !LOGICAL                    :: use_tg
  !
  !
  IF (dfft%has_task_groups) CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  ALLOCATE( aux( dfft%nnr ) )
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        CALL cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL fw_scatter( isgn ) ! forward scatter from stick to planes
     !
     CALL cft_2xy( f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes )
     !
  ELSE
     !
     IF ( isgn == -1 ) THEN
        !
        planes = dfft%iplp
        !
     ELSE IF ( isgn == -2 ) THEN
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL cft_2xy( f, dfft%my_nr3p, n1, n2, nx1, nx2, isgn, planes )
     !
     CALL bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
        !
     ELSE
        !
        CALL cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE fw_scatter( iopt )
!# 155 "fft_parallel_2d.f90"
     !Transpose data for the 2-D FFT on the x-y plane
     !
     !NOGRP*dfft%nnr: The length of aux and f
     !nr3x: The length of each Z-stick
     !aux: input - output
     !f: working space
     !isgn: type of scatter
     !dfft%nsw(me) holds the number of Z-sticks proc. me has.
     !dfft%nr3p: number of planes per processor
     !
     !
     USE fft_scatter_2d, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE fw_scatter
!# 183 "fft_parallel_2d.f90"
  !
!# 185 "fft_parallel_2d.f90"
  SUBROUTINE bw_scatter( iopt )
     !
     USE fft_scatter_2d, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE bw_scatter
  !
END SUBROUTINE tg_cft3s
!
!
!
!# 592 "fft_parallel_2d.f90"
!
END MODULE fft_parallel_2d
