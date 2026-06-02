!# 1 "fft_scalar.FFTW.f90"
!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!# 16 "fft_scalar.FFTW.f90"
!=----------------------------------------------------------------------=!
   MODULE fft_scalar_FFTW
!=----------------------------------------------------------------------=!
       USE fft_param
!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE iso_c_binding
       USE fftw_interfaces
!# 24 "fft_scalar.FFTW.f90"
       IMPLICIT NONE
       SAVE
       PRIVATE
       PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds
!# 29 "fft_scalar.FFTW.f90"
!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!
!# 33 "fft_scalar.FFTW.f90"
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z"
!
!
!
!=----------------------------------------------------------------------=!
!
!# 45 "fft_scalar.FFTW.f90"
   SUBROUTINE cft_1z(c, nsl, nz, ldz, isign, cout)
!# 47 "fft_scalar.FFTW.f90"
!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : backward (f(G)=>f(R)), isign < 0 : forward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available
!# 56 "fft_scalar.FFTW.f90"
     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz
!# 59 "fft_scalar.FFTW.f90"
     COMPLEX (DP) :: c(:), cout(:)
!# 61 "fft_scalar.FFTW.f90"
     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: found
!# 71 "fft_scalar.FFTW.f90"
     !   Pointers to the "C" structures containing FFT factors ( PLAN )
!# 73 "fft_scalar.FFTW.f90"
     TYPE(C_PTR), SAVE :: fw_planz( ndims ) = C_NULL_PTR
     TYPE(C_PTR), SAVE :: bw_planz( ndims ) = C_NULL_PTR
!# 76 "fft_scalar.FFTW.f90"
     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF
!# 80 "fft_scalar.FFTW.f90"
     !
     !   Here initialize table only if necessary
     !
     CALL lookup()
!# 85 "fft_scalar.FFTW.f90"
     IF( .NOT. found ) THEN
!# 87 "fft_scalar.FFTW.f90"
       !   no table exist for these parameters
       !   initialize a new one
!# 90 "fft_scalar.FFTW.f90"
       CALL init_plan()
!# 92 "fft_scalar.FFTW.f90"
     END IF
!# 94 "fft_scalar.FFTW.f90"
     !
     !   Now perform the FFTs using machine specific drivers
     !
!# 135 "fft_scalar.FFTW.f90"
     IF (isign < 0) THEN
        CALL FFT_Z_STICK(fw_planz( ip), c(1), ldz, nsl)
        tscale = 1.0_DP / nz
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
        CALL FFT_Z_STICK(bw_planz( ip), c(1), ldz, nsl)
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
     END IF
!# 150 "fft_scalar.FFTW.f90"
     RETURN
!# 152 "fft_scalar.FFTW.f90"
   CONTAINS
!# 154 "fft_scalar.FFTW.f90"
     SUBROUTINE lookup()
     DO ip = 1, ndims
        !   first check if there is already a table initialized
        !   for this combination of parameters
        found = ( nz == zdims(1,ip) )
        IF (found) EXIT
     END DO
     END SUBROUTINE lookup
!# 163 "fft_scalar.FFTW.f90"
     SUBROUTINE init_plan()
       IF( C_ASSOCIATED(fw_planz( icurrent)) ) CALL DESTROY_PLAN_1D( fw_planz( icurrent) )
       IF( C_ASSOCIATED(bw_planz( icurrent)) ) CALL DESTROY_PLAN_1D( bw_planz( icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_planz( icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_planz( icurrent), nz, idir)
       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan
!# 173 "fft_scalar.FFTW.f90"
   END SUBROUTINE cft_1z
!# 175 "fft_scalar.FFTW.f90"
!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!
!# 189 "fft_scalar.FFTW.f90"
   SUBROUTINE cft_2xy(r, nzl, nx, ny, ldx, ldy, isign, pl2ix)
!# 191 "fft_scalar.FFTW.f90"
!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : backward (f(G)=>f(R)), isign < 0 : forward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available
!# 201 "fft_scalar.FFTW.f90"
     IMPLICIT NONE
!# 203 "fft_scalar.FFTW.f90"
     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     COMPLEX (DP) :: r( : )
     INTEGER :: i, k, j, err, idir, ip, kk
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), found
!# 224 "fft_scalar.FFTW.f90"
     TYPE(C_PTR), SAVE :: fw_plan( 2, ndims ) = C_NULL_PTR
     TYPE(C_PTR), SAVE :: bw_plan( 2, ndims ) = C_NULL_PTR
!# 228 "fft_scalar.FFTW.f90"
     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF
!# 237 "fft_scalar.FFTW.f90"
     !
     !   Here initialize table only if necessary
     !
!# 241 "fft_scalar.FFTW.f90"
     CALL lookup()
!# 243 "fft_scalar.FFTW.f90"
     IF( .NOT. found ) THEN
!# 245 "fft_scalar.FFTW.f90"
       !   no table exist for these parameters
       !   initialize a new one
!# 248 "fft_scalar.FFTW.f90"
       CALL init_plan()
!# 250 "fft_scalar.FFTW.f90"
     END IF
!# 252 "fft_scalar.FFTW.f90"
     !
     !   Now perform the FFTs using machine specific drivers
     !
!# 356 "fft_scalar.FFTW.f90"
     IF( isign < 0 ) THEN
!# 358 "fft_scalar.FFTW.f90"
       CALL FFT_X_STICK( fw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )
!# 360 "fft_scalar.FFTW.f90"
       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK(fw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do
!# 369 "fft_scalar.FFTW.f90"
       tscale = 1.0_DP / ( nx * ny )
       r(1:ldx * ldy * nzl) = r(1:ldx * ldy * nzl) * tscale
!# 372 "fft_scalar.FFTW.f90"
     ELSE IF( isign > 0 ) THEN
!# 374 "fft_scalar.FFTW.f90"
       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK( bw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do
!# 383 "fft_scalar.FFTW.f90"
       CALL FFT_X_STICK( bw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )
!# 385 "fft_scalar.FFTW.f90"
    END IF
!# 393 "fft_scalar.FFTW.f90"
     RETURN
!# 395 "fft_scalar.FFTW.f90"
   CONTAINS
!# 397 "fft_scalar.FFTW.f90"
     SUBROUTINE lookup()
     DO ip = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       found = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       found = found .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (found) EXIT
     END DO
     END SUBROUTINE lookup
!# 407 "fft_scalar.FFTW.f90"
     SUBROUTINE init_plan()
!# 414 "fft_scalar.FFTW.f90"
       IF( C_ASSOCIATED(fw_plan( 2,icurrent)) )   CALL DESTROY_PLAN_1D( fw_plan( 2,icurrent) )
       IF( C_ASSOCIATED(bw_plan( 2,icurrent)) )   CALL DESTROY_PLAN_1D( bw_plan( 2,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2,icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2,icurrent), ny, idir)
!# 419 "fft_scalar.FFTW.f90"
       IF( C_ASSOCIATED(fw_plan( 1,icurrent)) ) CALL DESTROY_PLAN_1D( fw_plan( 1,icurrent) )
       IF( C_ASSOCIATED(bw_plan( 1,icurrent)) ) CALL DESTROY_PLAN_1D( bw_plan( 1,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1,icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1,icurrent), nx, idir)
!# 424 "fft_scalar.FFTW.f90"
       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan
!# 430 "fft_scalar.FFTW.f90"
   END SUBROUTINE cft_2xy
!# 432 "fft_scalar.FFTW.f90"
!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs
!
!
!
!=----------------------------------------------------------------------=!
!
!# 444 "fft_scalar.FFTW.f90"
   SUBROUTINE cfft3d( f, nx, ny, nz, ldx, ldy, ldz, howmany, isign )
!# 446 "fft_scalar.FFTW.f90"
  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available
!# 457 "fft_scalar.FFTW.f90"
     IMPLICIT NONE
!# 459 "fft_scalar.FFTW.f90"
     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1
!# 466 "fft_scalar.FFTW.f90"
     TYPE(C_PTR), save :: fw_plan(ndims) = C_NULL_PTR
     TYPE(C_PTR), save :: bw_plan(ndims) = C_NULL_PTR
!# 469 "fft_scalar.FFTW.f90"
     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3',' nz is less than 1 ', 1)
!# 476 "fft_scalar.FFTW.f90"
     !
     !   Here initialize table only if necessary
     !
     CALL lookup()
!# 481 "fft_scalar.FFTW.f90"
     IF( ip == -1 ) THEN
!# 483 "fft_scalar.FFTW.f90"
       !   no table exist for these parameters
       !   initialize a new one
!# 486 "fft_scalar.FFTW.f90"
       CALL init_plan()
!# 488 "fft_scalar.FFTW.f90"
     END IF
!# 490 "fft_scalar.FFTW.f90"
     !
     !   Now perform the 3D FFT using the machine specific driver
     !
!# 494 "fft_scalar.FFTW.f90"
     IF( isign < 0 ) THEN
!# 496 "fft_scalar.FFTW.f90"
       call FFTW_INPLACE_DRV_3D( fw_plan(ip), 1, f(1), 1, 1 )
!# 498 "fft_scalar.FFTW.f90"
       tscale = 1.0_DP / DBLE( nx * ny * nz )
       f(1:nx * ny * nz) = f(1:nx * ny * nz) * tscale
!# 501 "fft_scalar.FFTW.f90"
     ELSE IF( isign > 0 ) THEN
!# 503 "fft_scalar.FFTW.f90"
       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f(1), 1, 1 )
!# 505 "fft_scalar.FFTW.f90"
     END IF
!# 507 "fft_scalar.FFTW.f90"
     RETURN
!# 509 "fft_scalar.FFTW.f90"
   CONTAINS
!# 511 "fft_scalar.FFTW.f90"
     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup
!# 525 "fft_scalar.FFTW.f90"
     SUBROUTINE init_plan()
       IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
         call fftx_error__('cfft3','not implemented',1)
       IF( C_ASSOCIATED (fw_plan(icurrent)) ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       IF( C_ASSOCIATED (bw_plan(icurrent)) ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nx, ny, nz, idir)
       idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nx, ny, nz, idir)
       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan
!# 537 "fft_scalar.FFTW.f90"
   END SUBROUTINE cfft3d
!# 539 "fft_scalar.FFTW.f90"
!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!
!# 551 "fft_scalar.FFTW.f90"
SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, &
     do_fft_z, do_fft_y)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ldy)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none
!# 566 "fft_scalar.FFTW.f90"
  integer :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation
!# 572 "fft_scalar.FFTW.f90"
  complex(DP) :: f ( ldx * ldy * ldz * howmany )
  integer :: do_fft_z(:), do_fft_y(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj, h, ldh
  REAL(DP) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1
!# 581 "fft_scalar.FFTW.f90"
  TYPE(C_PTR), SAVE :: fw_plan ( 3, ndims ) = C_NULL_PTR
  TYPE(C_PTR), SAVE :: bw_plan ( 3, ndims ) = C_NULL_PTR
!# 584 "fft_scalar.FFTW.f90"
  tscale = 1.0_DP
  ldh = ldx * ldy * ldz
!# 587 "fft_scalar.FFTW.f90"
  IF( ny /= ldy ) &
    CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )
  IF( howmany < 1 ) &
    CALL fftx_error__(' cfft3ds ', ' howmany less than one ', 1 )
!# 592 "fft_scalar.FFTW.f90"
     CALL lookup()
!# 594 "fft_scalar.FFTW.f90"
     IF( ip == -1 ) THEN
!# 596 "fft_scalar.FFTW.f90"
       !   no table exist for these parameters
       !   initialize a new one
!# 599 "fft_scalar.FFTW.f90"
       CALL init_plan()
!# 601 "fft_scalar.FFTW.f90"
     END IF
!# 603 "fft_scalar.FFTW.f90"
     IF ( isign > 0 ) THEN
!# 605 "fft_scalar.FFTW.f90"
        DO h = 0, howmany - 1
           !
           !  k-direction ...
           !
!# 610 "fft_scalar.FFTW.f90"
           incx1 = ldx * ldy;  incx2 = 1;  m = 1
!# 612 "fft_scalar.FFTW.f90"
           do i =1, nx
              do j = 1, ny
                 ii = i + ldx * (j-1)
                 if ( do_fft_z(ii) > 0 ) then
                    call FFTW_INPLACE_DRV_1D( bw_plan( 3, ip), m, f( ii + h*ldh ), incx1, incx2 )
                 end if
              end do
           end do
!# 621 "fft_scalar.FFTW.f90"
           !
           !  ... j-direction ...
           !
!# 625 "fft_scalar.FFTW.f90"
           incx1 = ldx;  incx2 = ldx*ldy;  m = nz
!# 627 "fft_scalar.FFTW.f90"
           do i = 1, nx
              if ( do_fft_y( i ) == 1 ) then
                call FFTW_INPLACE_DRV_1D( bw_plan( 2, ip), m, f( i + h*ldh ), incx1, incx2 )
              endif
           enddo
!# 633 "fft_scalar.FFTW.f90"
           !
           !  ... i - direction
           !
!# 637 "fft_scalar.FFTW.f90"
           incx1 = 1;  incx2 = ldx;  m = ldy*nz
!# 639 "fft_scalar.FFTW.f90"
           call FFTW_INPLACE_DRV_1D( bw_plan( 1, ip), m, f( 1 + h*ldh ), incx1, incx2 )
!# 641 "fft_scalar.FFTW.f90"
        END DO
!# 643 "fft_scalar.FFTW.f90"
     ELSE
!# 645 "fft_scalar.FFTW.f90"
        DO h = 0, howmany - 1
           !
           !  i - direction ...
           !
!# 650 "fft_scalar.FFTW.f90"
           incx1 = 1;  incx2 = ldx;  m = ldy*nz
!# 652 "fft_scalar.FFTW.f90"
           call FFTW_INPLACE_DRV_1D( fw_plan( 1, ip), m, f( 1 + h*howmany ), incx1, incx2 )
!# 654 "fft_scalar.FFTW.f90"
           !
           !  ... j-direction ...
           !
!# 658 "fft_scalar.FFTW.f90"
           incx1 = ldx;  incx2 = ldx*ldy;  m = nz
!# 660 "fft_scalar.FFTW.f90"
           do i = 1, nx
              if ( do_fft_y ( i ) == 1 ) then
                call FFTW_INPLACE_DRV_1D( fw_plan( 2, ip), m, f( i + h*howmany ), incx1, incx2 )
              endif
           enddo
!# 666 "fft_scalar.FFTW.f90"
           !
           !  ... k-direction
           !
!# 670 "fft_scalar.FFTW.f90"
           incx1 = ldx * ny;  incx2 = 1;  m = 1
!# 672 "fft_scalar.FFTW.f90"
           do i = 1, nx
              do j = 1, ny
                 ii = i + ldx * (j -1)
                 if ( do_fft_z ( ii ) > 0 ) then
                    call FFTW_INPLACE_DRV_1D( fw_plan( 3, ip), m, f( ii + h*howmany ), incx1, incx2 )
                 end if
              end do
           end do
!# 681 "fft_scalar.FFTW.f90"
           f(h*howmany+1:h*howmany+ldx*ldy*nz) = f(h*howmany+1:h*howmany+ldx*ldy*nz) * (1.0_DP/(nx*ny*nz))
        END DO
!# 684 "fft_scalar.FFTW.f90"
     END IF
     RETURN
!# 687 "fft_scalar.FFTW.f90"
   CONTAINS
!# 689 "fft_scalar.FFTW.f90"
     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
           ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup
!# 702 "fft_scalar.FFTW.f90"
     SUBROUTINE init_plan()
       IF( C_ASSOCIATED (fw_plan( 1, icurrent)) ) CALL DESTROY_PLAN_1D( fw_plan( 1, icurrent) )
       IF( C_ASSOCIATED (bw_plan( 1, icurrent)) ) CALL DESTROY_PLAN_1D( bw_plan( 1, icurrent) )
       IF( C_ASSOCIATED (fw_plan( 2, icurrent)) ) CALL DESTROY_PLAN_1D( fw_plan( 2, icurrent) )
       IF( C_ASSOCIATED (bw_plan( 2, icurrent)) ) CALL DESTROY_PLAN_1D( bw_plan( 2, icurrent) )
       IF( C_ASSOCIATED (fw_plan( 3, icurrent)) ) CALL DESTROY_PLAN_1D( fw_plan( 3, icurrent) )
       IF( C_ASSOCIATED (bw_plan( 3, icurrent)) ) CALL DESTROY_PLAN_1D( bw_plan( 3, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1, icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1, icurrent), nx, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2, icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2, icurrent), ny, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 3, icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 3, icurrent), nz, idir)
!# 716 "fft_scalar.FFTW.f90"
       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan
!# 721 "fft_scalar.FFTW.f90"
   END SUBROUTINE cfft3ds
!=----------------------------------------------------------------------=!
 END MODULE fft_scalar_FFTW
!=----------------------------------------------------------------------=!
