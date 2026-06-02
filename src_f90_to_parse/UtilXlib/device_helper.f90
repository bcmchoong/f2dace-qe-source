!# 1 "device_helper.f90"
!
! Copyright (C) 2003-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file initiated by Carlo Cavazzoni 2020
!
! Purpose: collect miscellaneous subroutines to help dealing with 
!          accelerator devices
!# 13 "device_helper.f90"
SUBROUTINE MYDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!# 18 "device_helper.f90"
!     .. Scalar Arguments ..
    DOUBLE PRECISION ::  ALPHA
    INTEGER          ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    DOUBLE PRECISION :: A( LDA, * ), X( * ), Y( * )
!# 26 "device_helper.f90"
    CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!# 28 "device_helper.f90"
END SUBROUTINE MYDGER
!# 30 "device_helper.f90"
SUBROUTINE MYZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!# 35 "device_helper.f90"
!     .. Scalar Arguments ..
    COMPLEX*16, INTENT(IN) :: ALPHA
    INTEGER,    INTENT(IN) :: INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    COMPLEX*16 :: A( LDA, * ), X( * ), Y( * )
!# 44 "device_helper.f90"
    CALL ZGERC  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!# 47 "device_helper.f90"
END SUBROUTINE MYZGERC
!# 49 "device_helper.f90"
!=----------------------------------------------------------------------------=!
!# 51 "device_helper.f90"
SUBROUTINE MYDGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!# 56 "device_helper.f90"
    CHARACTER*1, INTENT(IN) ::        TRANSA, TRANSB
    INTEGER, INTENT(IN) ::            M, N, K, LDA, LDB, LDC
    DOUBLE PRECISION, INTENT(IN) ::   ALPHA, BETA
    DOUBLE PRECISION  :: A( LDA, * ), B( LDB, * ), C( LDC, * )
!# 64 "device_helper.f90"
    CALL dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!# 67 "device_helper.f90"
END SUBROUTINE MYDGEMM
!# 69 "device_helper.f90"
SUBROUTINE MYZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!# 74 "device_helper.f90"
    CHARACTER*1, INTENT(IN) ::        TRANSA, TRANSB
    INTEGER, INTENT(IN) ::            M, N, K, LDA, LDB, LDC
    COMPLEX*16, INTENT(IN) ::   ALPHA, BETA
    COMPLEX*16  :: A( LDA, * ), B( LDB, * ), C( LDC, * )
!# 82 "device_helper.f90"
    CALL zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
!# 85 "device_helper.f90"
END SUBROUTINE MYZGEMM
!# 87 "device_helper.f90"
SUBROUTINE MYDGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!# 92 "device_helper.f90"
      DOUBLE PRECISION, INTENT(IN) :: ALPHA,BETA
      INTEGER, INTENT(IN) :: INCX,INCY,LDA,M,N
      CHARACTER*1, INTENT(IN) :: TRANS
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!# 100 "device_helper.f90"
    CALL dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!# 102 "device_helper.f90"
END SUBROUTINE MYDGEMV
!# 104 "device_helper.f90"
SUBROUTINE MYZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!# 109 "device_helper.f90"
      COMPLEX*16, INTENT(IN) :: ALPHA,BETA
      INTEGER, INTENT(IN) :: INCX,INCY,LDA,M,N
      CHARACTER*1, INTENT(IN) :: TRANS
      COMPLEX*16 :: A(LDA,*),X(*),Y(*)
!# 117 "device_helper.f90"
    CALL zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!# 119 "device_helper.f90"
END SUBROUTINE MYZGEMV
!# 121 "device_helper.f90"
!=----------------------------------------------------------------------------=
!# 123 "device_helper.f90"
DOUBLE PRECISION FUNCTION MYDDOT(N,DX,INCX,DY,INCY)
!# 128 "device_helper.f90"
    INTEGER, INTENT(IN) :: INCX,INCY,N
    DOUBLE PRECISION, INTENT(IN) :: DX(*),DY(*)
!# 134 "device_helper.f90"
    DOUBLE PRECISION DDOT
    MYDDOT = DDOT(N,DX,INCX,DY,INCY)
!# 138 "device_helper.f90"
END FUNCTION MYDDOT
!# 140 "device_helper.f90"
! this is analogus to MYDDOT, but the result is on device
DOUBLE PRECISION FUNCTION MYDDOT_VECTOR_GPU(N,DX,DY)
!# 145 "device_helper.f90"
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: DX(*),DY(*)
    DOUBLE PRECISION :: RES
!# 172 "device_helper.f90"
    DOUBLE PRECISION DDOT
    MYDDOT_VECTOR_GPU = DDOT(N,DX,1,DY,1)
!# 175 "device_helper.f90"
END FUNCTION MYDDOT_VECTOR_GPU
!# 177 "device_helper.f90"
function MYDDOTv2 (n, dx, incx, dy, incy)
!# 182 "device_helper.f90"
implicit none
  DOUBLE PRECISION :: MYDDOTv2
  integer :: n, incx, incy
  DOUBLE PRECISION, dimension(*)  :: dx, dy
!# 194 "device_helper.f90"
  DOUBLE PRECISION DDOT
  MYDDOTv2=DDOT(n, dx, incx, dy, incy)
!# 198 "device_helper.f90"
  return
end function MYDDOTv2
!# 201 "device_helper.f90"
subroutine MYDDOTv3 (n, dx, incx, dy, incy, result)
!# 206 "device_helper.f90"
  implicit none
  integer :: n, incx, incy
  DOUBLE PRECISION, dimension(*)  :: dx, dy
  DOUBLE PRECISION :: result
!# 218 "device_helper.f90"
  DOUBLE PRECISION DDOT
  result=DDOT(n, dx, incx, dy, incy)
!# 222 "device_helper.f90"
  return
end subroutine MYDDOTv3
!# 225 "device_helper.f90"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MYDTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) 
!# 230 "device_helper.f90"
implicit none
  character*1 :: side, uplo, transa, diag 
  integer :: m, n, lda, ldb 
  DOUBLE PRECISION :: alpha 
  DOUBLE PRECISION, dimension(lda, *) :: a 
  DOUBLE PRECISION, dimension(ldb, *) :: b 
!# 240 "device_helper.f90"
  call DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)  
!# 242 "device_helper.f90"
  return
end subroutine MYDTRSM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE COMPLEX function MYZDOTC(n, zx, incx, zy, incy) 
!# 249 "device_helper.f90"
implicit none
  integer :: n, incx, incy 
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 256 "device_helper.f90"
  DOUBLE COMPLEX, EXTERNAL :: ZDOTC
  MYZDOTC = ZDOTC(n, zx, incx, zy, incy)  
!# 259 "device_helper.f90"
  return
end function MYZDOTC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE COMPLEX function MYZDOTC_VECTOR_GPU(n, zx, zy) 
!# 266 "device_helper.f90"
implicit none
  integer :: n
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 285 "device_helper.f90"
  DOUBLE COMPLEX, EXTERNAL :: ZDOTC
  MYZDOTC_VECTOR_GPU = ZDOTC(n, zx, 1, zy, 1)  
!# 288 "device_helper.f90"
  return
end function MYZDOTC_VECTOR_GPU 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSWAP(n, zx, incx, zy, incy) 
!# 295 "device_helper.f90"
implicit none
  integer :: n, incx, incy 
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 302 "device_helper.f90"
  CALL ZSWAP(n, zx, incx, zy, incy)  
!# 304 "device_helper.f90"
  return
END SUBROUTINE MYZSWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSWAP_VECTOR_GPU(n, zx, zy) 
!# 311 "device_helper.f90"
implicit none
  integer :: n
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 326 "device_helper.f90"
  CALL ZSWAP(n, zx, 1, zy, 1)  
!# 328 "device_helper.f90"
  return
END SUBROUTINE MYZSWAP_VECTOR_GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZCOPY(n, zx, incx, zy, incy)
!# 335 "device_helper.f90"
IMPLICIT NONE 
  INTEGER :: n, incx, incy
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 342 "device_helper.f90"
  CALL ZCOPY(n, zx, incx, zy, incy)  
!# 344 "device_helper.f90"
  RETURN
END SUBROUTINE MYZCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZAXPY(n, za, zx, incx, zy, incy)
!# 351 "device_helper.f90"
IMPLICIT NONE 
  INTEGER :: n, incx, incy
  DOUBLE COMPLEX :: za  
  DOUBLE COMPLEX, dimension(*) :: zx, zy
!# 359 "device_helper.f90"
  CALL ZAXPY(n, za, zx, incx, zy, incy)  
!# 361 "device_helper.f90"
  RETURN
END SUBROUTINE MYZAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZDSCAL(n, da, zx, incx)
!# 368 "device_helper.f90"
IMPLICIT NONE 
  INTEGER :: n, incx
  DOUBLE PRECISION :: da
  DOUBLE COMPLEX, dimension(*) :: zx
!# 376 "device_helper.f90"
  CALL ZDSCAL(n, da, zx, incx)
!# 378 "device_helper.f90"
  RETURN
END SUBROUTINE MYZDSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSCAL(n, za, zx, incx)
!# 385 "device_helper.f90"
IMPLICIT NONE 
  INTEGER :: n, incx
  DOUBLE COMPLEX :: za
  DOUBLE COMPLEX, dimension(*) :: zx
!# 393 "device_helper.f90"
  CALL ZSCAL(n, za, zx, incx)
!# 395 "device_helper.f90"
  RETURN
END SUBROUTINE MYZSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDCOPY(n, x, incx, y, incy)
!# 402 "device_helper.f90"
IMPLICIT NONE
  INTEGER :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)   :: x(*)
  DOUBLE PRECISION, INTENT(OUT)  :: y(*)
!# 410 "device_helper.f90"
  call DCOPY(n, x, incx, y, incy)
!# 412 "device_helper.f90"
  RETURN
END SUBROUTINE MYDCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDAXPY(n, a, x, incx, y, incy)
!# 419 "device_helper.f90"
IMPLICIT NONE
  INTEGER :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)  :: a
  DOUBLE PRECISION, INTENT(IN)  :: x(*) 
  DOUBLE PRECISION, INTENT(OUT) :: y(*) 
!# 428 "device_helper.f90"
  call DAXPY( n, a, x, incx, y, incy)
!# 430 "device_helper.f90"
  RETURN
END SUBROUTINE MYDAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSCAL(n, a, x, incx)
!# 437 "device_helper.f90"
IMPLICIT NONE
  integer :: n, incx
  DOUBLE PRECISION :: a
  DOUBLE PRECISION, dimension(*)  :: x
!# 445 "device_helper.f90"
  call DSCAL(n, a, x, incx)
!# 447 "device_helper.f90"
  RETURN
END SUBROUTINE MYDSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSWAP(n, dx, incx, dy, incy) 
!# 454 "device_helper.f90"
implicit none
  integer :: n, incx, incy 
  REAL(8), dimension(*) :: dx, dy
!# 461 "device_helper.f90"
  CALL DSWAP(n, dx, incx, dy, incy)  
!# 463 "device_helper.f90"
  return
END SUBROUTINE MYDSWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSWAP_VECTOR_GPU(n, dx, dy) 
!# 470 "device_helper.f90"
implicit none
  integer :: n
  DOUBLE PRECISION, dimension(*) :: dx, dy
!# 485 "device_helper.f90"
  CALL DSWAP(n, dx, 1, dy, 1)  
!# 487 "device_helper.f90"
  return
END SUBROUTINE MYDSWAP_VECTOR_GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
