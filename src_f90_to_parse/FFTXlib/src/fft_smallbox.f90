!# 1 "fft_smallbox.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 10 "fft_smallbox.f90"
!=----------------------------------------------------------------------=!
   MODULE fft_smallbox
!=----------------------------------------------------------------------=!
!# 14 "fft_smallbox.f90"
!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE iso_c_binding
       IMPLICIT NONE
       SAVE
!# 19 "fft_smallbox.f90"
        PRIVATE
        PUBLIC :: cft_b, cft_b_omp_init, cft_b_omp
!# 22 "fft_smallbox.f90"
! ...   Local Parameter
!# 24 "fft_smallbox.f90"
        INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
!# 26 "fft_smallbox.f90"
        !   ndims   Number of different FFT tables that the module
        !           could keep into memory without reinitialization
!# 29 "fft_smallbox.f90"
        INTEGER, PARAMETER :: ndims = 3
!# 31 "fft_smallbox.f90"
        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)
!# 35 "fft_smallbox.f90"
        INTEGER   :: cft_b_dims( 3 )
!$omp threadprivate (cft_b_dims)
        TYPE(C_PTR) :: cft_b_bw_planz = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planz)
        TYPE(C_PTR) :: cft_b_bw_planx = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planx)
        TYPE(C_PTR) :: cft_b_bw_plany = C_NULL_PTR
!$omp threadprivate (cft_b_bw_plany)
!# 44 "fft_smallbox.f90"
!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!
!# 48 "fft_smallbox.f90"
!
!=----------------------------------------------------------------------=!
!
!
!
!         3D parallel FFT on sub-grids
!
!
!
!=----------------------------------------------------------------------=!
!
!# 60 "fft_smallbox.f90"
   SUBROUTINE cft_b ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )
!# 62 "fft_smallbox.f90"
!     driver routine for 3d complex fft's on box grid, parallel case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented for FFTW, only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
      USE fftw_interfaces
      implicit none
      integer nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(dp) :: f(:)
!# 75 "fft_smallbox.f90"
      integer isign, naux, ibid, k
      integer nplanes
      real(DP) :: tscale
!# 79 "fft_smallbox.f90"
      integer :: ip, i, first_index, how_many_y
      integer, save :: icurrent = 1
      integer, save :: dims( 3, ndims ) = -1
!# 83 "fft_smallbox.f90"
      TYPE(C_PTR), save :: bw_planz(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_planx(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_plany(  ndims ) = C_NULL_PTR
!# 87 "fft_smallbox.f90"
      isign = -sgn
      tscale = 1.0_DP
!# 90 "fft_smallbox.f90"
      if ( isign > 0 ) then
         call fftx_error__('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1
!# 99 "fft_smallbox.f90"
      !
      !   Here initialize table only if necessary
      !
!# 103 "fft_smallbox.f90"
      ip = -1
      DO i = 1, ndims
!# 106 "fft_smallbox.f90"
        !   first check if there is already a table initialized
        !   for this combination of parameters
!# 109 "fft_smallbox.f90"
        IF ( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. ( nz == dims(3,i) ) ) THEN
           ip = i
           EXIT
        END IF
!# 114 "fft_smallbox.f90"
      END DO
!# 116 "fft_smallbox.f90"
      IF( ip == -1 ) THEN
!# 118 "fft_smallbox.f90"
        !   no table exist for these parameters
        !   initialize a new one
!# 121 "fft_smallbox.f90"
        if ( C_ASSOCIATED(bw_planz(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planz(icurrent) )
        call CREATE_PLAN_1D( bw_planz(icurrent), nz, 1 )
!# 125 "fft_smallbox.f90"
        if ( C_ASSOCIATED(bw_planx(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planx(icurrent) )
        call CREATE_PLAN_1D( bw_planx(icurrent), nx, 1 )
!# 129 "fft_smallbox.f90"
        if ( C_ASSOCIATED(bw_plany(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_plany(icurrent) )
        call CREATE_PLAN_1D( bw_plany(icurrent), ny, 1 )
!# 133 "fft_smallbox.f90"
!        if ( C_ASSOCIATED(bw_planxy(icurrent)) ) &
!             call DESTROY_PLAN_2D( bw_planxy(icurrent) )
!        call CREATE_PLAN_2D( bw_planxy(icurrent), nx, ny, 1 )
!
        dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1
!# 141 "fft_smallbox.f90"
      END IF
!# 143 "fft_smallbox.f90"
      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1 )
     
      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FFTW_INPLACE_DRV_1D( bw_plany(ip), nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FFTW_INPLACE_DRV_1D( bw_planx(ip), how_many_y, f(first_index), 1, ldx )
!# 160 "fft_smallbox.f90"
      end do   
!# 162 "fft_smallbox.f90"
      RETURN
   END SUBROUTINE cft_b
!# 165 "fft_smallbox.f90"
!
!=----------------------------------------------------------------------=!
!
!
!
!   3D parallel FFT on sub-grids, to be called inside OpenMP region
!
!
!
!=----------------------------------------------------------------------=!
!
!# 177 "fft_smallbox.f90"
   SUBROUTINE cft_b_omp_init ( nx, ny, nz )
!# 179 "fft_smallbox.f90"
!     driver routine for 3d complex fft's on box grid, init subroutine
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz
      !
      !   Here initialize table 
      !
!$omp parallel
!# 189 "fft_smallbox.f90"
      IF( .NOT. C_ASSOCIATED(cft_b_bw_planz) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planz, nz, 1 )
         cft_b_dims(3) = nz
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_planx) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planx, nx, 1 )
         cft_b_dims(1) = nx
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_plany, ny, 1 )
         cft_b_dims(2) = ny
      END IF
!# 202 "fft_smallbox.f90"
!$omp end parallel
!# 204 "fft_smallbox.f90"
     RETURN
   END SUBROUTINE cft_b_omp_init
!# 208 "fft_smallbox.f90"
   SUBROUTINE cft_b_omp ( f, nx, ny, nz, ldx, ldy, ldz, imin2, imax2, imin3, imax3, sgn )
!# 210 "fft_smallbox.f90"
!     driver routine for 3d complex fft's on box grid, parallel (MPI+OpenMP) case
!     fft along z for all xy values 
!     fft along y is done only for the local  z values: i.e. z-planes with imin3 <= nz <= imax3.
!     fft along x is done only for the local yz values: i.e. z-planes with imin3 <= nz <= imax3 
!                                                        and y-planes with imin2 <= ny <= imax2.
!     implemented ONLY for internal fftw, and only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
!     This driver is meant for calls inside parallel OpenMP sections
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz,ldx,ldy,ldz,imin2,imax2,imin3,imax3,sgn
      complex(dp) :: f(:)
!# 225 "fft_smallbox.f90"
      INTEGER, SAVE :: k, first_index, how_many_y
!$omp threadprivate (k,first_index,how_many_y)
!# 228 "fft_smallbox.f90"
      if ( -sgn > 0 ) then
         CALL fftx_error__('cft_b_omp','forward transform not implemented',1)
      end if
!# 232 "fft_smallbox.f90"
      IF ( .NOT. C_ASSOCIATED(cft_b_bw_planz) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_planx) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL fftx_error__('cft_b_omp','plan not initialized',1)
      END IF
!# 238 "fft_smallbox.f90"
      !  consistency check
!# 240 "fft_smallbox.f90"
      IF ( ( nx /= cft_b_dims(1) ) .or. ( ny /= cft_b_dims(2) ) .or. ( nz /= cft_b_dims(3) ) ) THEN
         CALL fftx_error__('cft_b_omp', 'dimensions are inconsistent with the existing plan',1) 
      END IF
!# 244 "fft_smallbox.f90"
      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( cft_b_bw_planz, ldx*ldy, f(1), ldx*ldy, 1 )
!# 249 "fft_smallbox.f90"
      do k = imin3, imax3
      !
      !  fft along Y
      !
        first_index = (k-1)*ldx*ldy + 1
        call FFTW_INPLACE_DRV_1D( cft_b_bw_plany, nx, f(first_index), ldx, 1 )
      !
      !  fft along X
      !
        first_index = first_index + (imin2-1)*ldx ; how_many_y = imax2 + 1 - imin2
        call FFTW_INPLACE_DRV_1D( cft_b_bw_planx, how_many_y, f(first_index), 1, ldx )
      end do   
!# 262 "fft_smallbox.f90"
     RETURN
   END SUBROUTINE cft_b_omp
!# 266 "fft_smallbox.f90"
!=----------------------------------------------------------------------=!
   END MODULE fft_smallbox
!=----------------------------------------------------------------------=!
