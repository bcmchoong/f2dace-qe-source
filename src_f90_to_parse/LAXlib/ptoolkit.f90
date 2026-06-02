!# 1 "ptoolkit.f90"
!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ---------------------------------------------------------------------------------
!# 10 "ptoolkit.f90"
MODULE laxlib_ptoolkit
  IMPLICIT NONE
  SAVE
CONTAINS
!# 15 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmred_x_x( na, a, lda, desca, nb, b, ldb, descb )
   !
   !! Double precision SQuare Matrix REDistribution
   !! 
   !! Copy a global "na * na" matrix locally stored in "a",
   !! and distributed as described by "desca", into a larger
   !! global "nb * nb" matrix stored in "b" and distributed
   !! as described in "descb".
   !!
   ! 
   ! If you want to read, get prepared for an headache!
   ! Written struggling by Carlo Cavazzoni.
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: na
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a
   REAL(DP)            :: a(lda,lda)  
   !!  matrix to be redistributed into b
   TYPE(la_descriptor), INTENT(IN) :: desca
   !! laxlib descriptor of matrix a
   INTEGER, INTENT(IN) :: nb
   !! global dimension of matrix b   
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of matrix b
   REAL(DP)            :: b(ldb,ldb)
   !! redistributed matrix
   TYPE(la_descriptor), INTENT(IN) :: descb
   !! laxlib descriptor of matrix b
!# 51 "ptoolkit.f90"
   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca%npr )
   INTEGER :: ncsnd( desca%npr )
   INTEGER :: displ( desca%npr )
   INTEGER :: irb_new( desca%npr )
   INTEGER :: ire_new( desca%npr )
   INTEGER :: icb_new( desca%npr )
   INTEGER :: ice_new( desca%npr )
   REAL(DP), ALLOCATABLE :: buf(:)
   REAL(DP), ALLOCATABLE :: ab(:,:)
   REAL(DP), ALLOCATABLE :: tst1(:,:)
   REAL(DP), ALLOCATABLE :: tst2(:,:)
!# 74 "ptoolkit.f90"
   IF( desca%active_node <= 0 ) THEN
      RETURN
   END IF
!# 78 "ptoolkit.f90"
   ! preliminary consistency checks
!# 80 "ptoolkit.f90"
   IF( nb < na ) &
      CALL lax_error__( " dsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb%n ) &
      CALL lax_error__( " dsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca%n ) &
      CALL lax_error__( " dsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb%nrcx ) &
      CALL lax_error__( " dsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca%nrcx ) &
      CALL lax_error__( " dsqmred ", " wrong leading dim lda ", lda )
!# 91 "ptoolkit.f90"
   npr   = desca%npr
   myrow = desca%myr
   npc   = desca%npc
   mycol = desca%myc
   comm  = desca%comm
!# 410 "ptoolkit.f90"
   RETURN
END SUBROUTINE laxlib_dsqmred_x_x
!# 413 "ptoolkit.f90"
SUBROUTINE laxlib_zsqmred_x_x( na, a, lda, desca, nb, b, ldb, descb )
   !
   !! double complex (Z) SQuare Matrix REDistribution
   !! 
   !! Copy a global "na * na" matrix locally stored in "a",
   !! and distributed as described by "desca", into a larger
   !! global "nb * nb" matrix stored in "b" and distributed
   !! as described in "descb".
   !!
   ! 
   ! If you want to read, get prepared for an headache!
   ! Written struggling by Carlo Cavazzoni.
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: na
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a
   COMPLEX(DP)            :: a(lda,lda)
   !!  matrix to be redistributed into b
   TYPE(la_descriptor), INTENT(IN) :: desca
   !! laxlib descriptor of matrix a
   INTEGER, INTENT(IN) :: nb
   !! global dimension of matrix b   
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of matrix b
   COMPLEX(DP)            :: b(ldb,ldb)
   !! redistributed matrix
   TYPE(la_descriptor), INTENT(IN) :: descb
   !! laxlib descriptor matrix b
!# 449 "ptoolkit.f90"
   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca%npr )
   INTEGER :: ncsnd( desca%npr )
   INTEGER :: displ( desca%npr )
   INTEGER :: irb_new( desca%npr )
   INTEGER :: ire_new( desca%npr )
   INTEGER :: icb_new( desca%npr )
   INTEGER :: ice_new( desca%npr )
   COMPLEX(DP), ALLOCATABLE :: buf(:)
   COMPLEX(DP), ALLOCATABLE :: ab(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst1(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst2(:,:)
!# 472 "ptoolkit.f90"
   IF( desca%active_node <= 0 ) THEN
      RETURN
   END IF
!# 476 "ptoolkit.f90"
   ! preliminary consistency checks
!# 478 "ptoolkit.f90"
   IF( nb < na ) &
      CALL lax_error__( " zsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb%n ) &
      CALL lax_error__( " zsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca%n ) &
      CALL lax_error__( " zsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb%nrcx ) &
      CALL lax_error__( " zsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca%nrcx ) &
      CALL lax_error__( " zsqmred ", " wrong leading dim lda ", lda )
!# 489 "ptoolkit.f90"
   npr   = desca%npr
   myrow = desca%myr
   npc   = desca%npc
   mycol = desca%myc
   comm  = desca%comm
!# 788 "ptoolkit.f90"
   RETURN
END SUBROUTINE laxlib_zsqmred_x_x
!# 791 "ptoolkit.f90"
END MODULE laxlib_ptoolkit
!# 794 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 797 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmdst_x( n, ar, ldar, a, lda, desc )
  !
  !!  Double precision SQuare Matrix DiSTribution
  !!  This subroutine take a replicated square matrix "ar" and distribute it
  !!  across processors as described by descriptor "desc"
  !
  USE laxlib_descriptor
  !
  implicit none
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! global dimension
  INTEGER, INTENT(IN) :: ldar
  !! leading dimension of matrix ar
  REAL(DP)            :: ar(ldar,*)  
  !! matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of matrix a
  REAL(DP)            :: a(lda,*)
  !! distributed matrix a
  TYPE(la_descriptor), INTENT(IN) :: desc
  !! laxlib descriptor for matrix a
  !!
  !
  REAL(DP), PARAMETER :: zero = 0_DP
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc%active_node <= 0 ) THEN
     RETURN
  END IF
!# 830 "ptoolkit.f90"
  nx  = desc%nrcx
  ir  = desc%ir
  ic  = desc%ic
  nr  = desc%nr
  nc  = desc%nc
!# 836 "ptoolkit.f90"
  IF( lda < nx ) &
     CALL lax_error__( " dsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc%n ) &
     CALL lax_error__( " dsqmdst ", " inconsistent dimension n ", n )
!# 841 "ptoolkit.f90"
  DO j = 1, nc
     DO i = 1, nr
        a( i, j ) = ar( i + ir - 1, j + ic - 1 )
     END DO
     DO i = nr+1, nx
        a( i, j ) = zero
     END DO
  END DO
  DO j = nc + 1, nx
     DO i = 1, nx
        a( i, j ) = zero
     END DO
  END DO
!# 855 "ptoolkit.f90"
  RETURN
!# 857 "ptoolkit.f90"
END SUBROUTINE laxlib_dsqmdst_x
!# 860 "ptoolkit.f90"
SUBROUTINE laxlib_zsqmdst_x( n, ar, ldar, a, lda, desc )
  !
  !! double complex (Z) SQuare Matrix DiSTribution
  !! This subroutine take a replicated square matrix "ar" and distribute it
  !! across processors as described by descriptor "desc"
  !
  USE laxlib_descriptor
  !
  implicit none
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! global dimension
  INTEGER, INTENT(IN) :: ldar
  !! leading dimension of matrix ar
  COMPLEX(DP)            :: ar(ldar,*)
  !! matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of matrix a
  COMPLEX(DP)            :: a(lda,*)
  !! distributed matrix a
  TYPE(la_descriptor), INTENT(IN) :: desc
  !! laxlib descriptor for matrix a
  !!
  !
  COMPLEX(DP), PARAMETER :: zero = ( 0_DP , 0_DP )
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc%active_node <= 0 ) THEN
     RETURN
  END IF
!# 893 "ptoolkit.f90"
  nx  = desc%nrcx
  ir  = desc%ir
  ic  = desc%ic
  nr  = desc%nr
  nc  = desc%nc
!# 899 "ptoolkit.f90"
  IF( lda < nx ) &
     CALL lax_error__( " zsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc%n ) &
     CALL lax_error__( " zsqmdst ", " inconsistent dimension n ", n )
!# 904 "ptoolkit.f90"
  DO j = 1, nc
     DO i = 1, nr
        a( i, j ) = ar( i + ir - 1, j + ic - 1 )
     END DO
     DO i = nr+1, nx
        a( i, j ) = zero
     END DO
  END DO
  DO j = nc + 1, nx
     DO i = 1, nx
        a( i, j ) = zero
     END DO
  END DO
!# 918 "ptoolkit.f90"
  RETURN
!# 920 "ptoolkit.f90"
END SUBROUTINE laxlib_zsqmdst_x
!# 922 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 924 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmcll_x( n, a, lda, ar, ldar, desc, comm )
  !
  !! Double precision SQuare Matrix CoLLect
  !! This sub. take a distributed square matrix "a" and collect 
  !! the block assigned to processors into a replicated matrix "ar",
  !! matrix is distributed as described by descriptor desc
  !!
  !
  USE laxlib_descriptor
  USE laxlib_parallel_include
  !
  implicit none
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! global dimension
  INTEGER, INTENT(IN) :: ldar
  !! leading dimension of matrix ar
  REAL(DP)            :: ar(ldar,*)
  !! matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of matrix a
  REAL(DP)            :: a(lda,*)
  !! distributed matrix a
  TYPE(la_descriptor), INTENT(IN) :: desc
  !! laxlib descriptor for matrix a
  INTEGER, INTENT(IN) :: comm
  !! mpi communicator
  !
  INTEGER :: i, j
!# 1013 "ptoolkit.f90"
  DO j = 1, n
     DO i = 1, n
        ar( i, j ) = a( i, j )
     END DO
  END DO
!# 1021 "ptoolkit.f90"
  RETURN
END SUBROUTINE laxlib_dsqmcll_x
!# 1025 "ptoolkit.f90"
SUBROUTINE laxlib_zsqmcll_x( n, a, lda, ar, ldar, desc, comm )
  !
  !  double complex (Z) SQuare Matrix CoLLect
  !  This sub. take a distributed square matrix "a" and collect 
  !  the block assigned to processors into a replicated matrix "ar",
  !  matrix is distributed as described by descriptor desc
  !
  USE laxlib_descriptor
  USE laxlib_parallel_include
  !
  implicit none
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! global dimension
  INTEGER, INTENT(IN) :: ldar
  !! leading dimension of matrix ar
  COMPLEX(DP)            :: ar(ldar,*)
  !! matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of matrix a
  COMPLEX(DP)            :: a(lda,*)
  !! distributed matrix a
  TYPE(la_descriptor), INTENT(IN) :: desc
  !! laxlib descriptor for matrix a
  INTEGER, INTENT(IN) :: comm
  !! mpi communicator
  !
  INTEGER :: i, j
!# 1113 "ptoolkit.f90"
  DO j = 1, n
     DO i = 1, n
        ar( i, j ) = a( i, j )
     END DO
  END DO
!# 1121 "ptoolkit.f90"
  RETURN
END SUBROUTINE laxlib_zsqmcll_x
!# 1125 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1127 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmwpb_x( n, a, lda, desc )
   !
   !! Double precision SQuare Matrix WiPe Border subroutine
   !! initialize to zero the distributed matrix border
   !
   USE laxlib_descriptor
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a 
   REAL(DP)            :: a(lda,*)  
   !! distributed matrix a
   TYPE(la_descriptor), INTENT(IN) :: desc
   !! laxlib descriptor
   !
   INTEGER :: i, j
   !
   DO j = 1, desc%nc
      DO i = desc%nr + 1, desc%nrcx
         a( i, j ) = 0_DP
      END DO
   END DO
   DO j = desc%nc + 1, desc%nrcx
      DO i = 1, desc%nrcx
         a( i, j ) = 0_DP
      END DO
   END DO
   !
   RETURN
END SUBROUTINE laxlib_dsqmwpb_x
!# 1162 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1164 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmsym_x( n, a, lda, idesc )
   !
   !! Double precision SQuare Matrix SYMmetrization
   !!
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a 
   REAL(DP)            :: a(lda,*)
   !! distributed matrix a
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
!# 1189 "ptoolkit.f90"
   INTEGER :: i, j
   INTEGER :: comm 
   INTEGER :: nr, nc, dest, sreq, ierr, sour
   REAL(DP) :: atmp
!# 1265 "ptoolkit.f90"
   DO j = 1, n
      !
      DO i = j + 1, n
         !
         a(i,j) = a(j,i)
         !
      END DO
      !
   END DO
!# 1277 "ptoolkit.f90"
   RETURN
END SUBROUTINE laxlib_dsqmsym_x
!# 1280 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1428 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1430 "ptoolkit.f90"
SUBROUTINE laxlib_zsqmher_x( n, a, lda, idesc )
   !
   !! double complex (Z) SQuare Matrix HERmitianize
   !!
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a 
   COMPLEX(DP)            :: a(lda,lda)
   !! distributed matrix a
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
!# 1455 "ptoolkit.f90"
   INTEGER :: i, j
   INTEGER :: comm, myid
   INTEGER :: nr, nc, dest, sreq, ierr, sour
   COMPLEX(DP) :: atmp
   COMPLEX(DP), ALLOCATABLE :: tst1(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst2(:,:)
!# 1570 "ptoolkit.f90"
   DO j = 1, n
      !
      a(j,j) = CMPLX( DBLE( a(j,j) ), 0_DP, KIND=DP )
      !
      DO i = j + 1, n
         !
         a(i,j) = CONJG( a(j,i) )
         !
      END DO
      !
   END DO
!# 1584 "ptoolkit.f90"
   RETURN
END SUBROUTINE laxlib_zsqmher_x
!# 1588 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1590 "ptoolkit.f90"
SUBROUTINE laxlib_dsqmred_x( na, a, lda, idesca, nb, b, ldb, idescb )
   !
   !! Double precision SQuare Matrix REDistribution
   !! 
   !! Copy a global "na * na" matrix locally stored in "a",
   !! and distributed as described by integer "idesca", into a larger
   !! global "nb * nb" matrix stored in "b" and distributed
   !! as described in integer "idescb".
   !
   !
   USE laxlib_descriptor
   USE laxlib_ptoolkit
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: na
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a
   REAL(DP)            :: a(lda,lda)
   !! matrix to be redistributed into b
   INTEGER, INTENT(IN) :: idesca(LAX_DESC_SIZE)
   !! integer laxlib descriptor matrix a
   INTEGER, INTENT(IN) :: nb
   !! global dimension of matrix b   
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of matrix b
   REAL(DP)            :: b(ldb,ldb)
   !! redistributed matrix
   INTEGER, INTENT(IN) :: idescb(LAX_DESC_SIZE)
   !! integer laxlib descriptor matrix b
   !
   TYPE(la_descriptor) :: desca
   TYPE(la_descriptor) :: descb
   !
   CALL laxlib_intarray_to_desc(desca,idesca)
   CALL laxlib_intarray_to_desc(descb,idescb)
   CALL laxlib_dsqmred_x_x( na, a, lda, desca, nb, b, ldb, descb )
END SUBROUTINE
!# 1633 "ptoolkit.f90"
SUBROUTINE laxlib_zsqmred_x( na, a, lda, idesca, nb, b, ldb, idescb )
   !
   !! double complex (Z) SQuare Matrix REDistribution
   !! 
   !! Copy a global "na * na" matrix locally stored in "a",
   !! and distributed as described by integer "idesca", into a larger
   !! global "nb * nb" matrix stored in "b" and distributed
   !! as described in integer "idescb".
   !!
   USE laxlib_descriptor
   USE laxlib_ptoolkit
   !
   IMPLICIT NONE
   !
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: na
   !! global dimension of matrix a
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of matrix a
   COMPLEX(DP)            :: a(lda,lda)
   !! matrix to be redistributed into b
   INTEGER, INTENT(IN) :: idesca(LAX_DESC_SIZE)
   !! integer laxlib descriptor matrix a
   INTEGER, INTENT(IN) :: nb
   !! global dimension of matrix b   
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of matrix b
   COMPLEX(DP)            :: b(ldb,ldb)
   !! redistributed matrix
   INTEGER, INTENT(IN) :: idescb(LAX_DESC_SIZE)
   !! integer laxlib descriptor matrix b
   !
   TYPE(la_descriptor) :: desca
   TYPE(la_descriptor) :: descb
   !
   CALL laxlib_intarray_to_desc(desca,idesca)
   CALL laxlib_intarray_to_desc(descb,idescb)
   CALL laxlib_zsqmred_x_x( na, a, lda, desca, nb, b, ldb, descb )
END SUBROUTINE
!# 1676 "ptoolkit.f90"
! ---------------------------------------------------------------------------------
!# 1679 "ptoolkit.f90"
SUBROUTINE rep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !!  Parallel matrix multiplication with replicated matrix
  !!
  !!  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !!
  !!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !!
  !!  WHERE  OP( X ) IS ONE OF
  !!
  !!     OP( X ) = X   OR   OP( X ) = X',
  !!
  !!  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !!  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !  written by Carlo Cavazzoni
  !
  USE laxlib_parallel_include
  implicit none
  include 'laxlib_kinds.fh'
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa
  !! specifies the form of op( A ) to be used in the matrix multiplication as follows:
  !! 'N' or 'n',  op( A ) = A.
  !! 'T' or 't',  op( A ) = A**T.
  !! 'C' or 'c',  op( A ) = A**T. 
  CHARACTER(LEN=1), INTENT(IN) :: transb
  !! specifies the form of op( B ) to be used in the matrix multiplication as
  !follows:
  !! 'N' or 'n',  op( B ) = B.
  !! 'T' or 't',  op( B ) = B**T.
  !! 'C' or 'c',  op( B ) = B**T.
  INTEGER, INTENT(IN) :: m
  !! number of rows of the matrix A and C
  INTEGER, INTENT(IN) :: n
  !! number of columns of the matrix B and C
  INTEGER, INTENT(IN) :: k
  !! number of columns of A and rows of B
  REAL(DP), INTENT(IN) :: alpha
  !! scalar alpha
  REAL(DP), INTENT(IN) :: beta
  !! scalar beta
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of A
  INTEGER, INTENT(IN) :: ldb
  !! leading dimension of B
  INTEGER, INTENT(IN) :: ldc
  !! leading dimension of C
  REAL(DP) :: a(lda,*)
  !! matrix A
  REAL(DP) :: b(ldb,*)
  !! matrix B
  REAL(DP) :: c(ldc,*)
  !! matrix C
  INTEGER, INTENT(IN) :: comm
  !! mpi communicator
  !
!# 1858 "ptoolkit.f90"
     !  if we are not compiling with __MPI this is equivalent to a blas call
!# 1860 "ptoolkit.f90"
     CALL dgemm(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)
!# 1864 "ptoolkit.f90"
  RETURN
!# 1866 "ptoolkit.f90"
END SUBROUTINE rep_matmul_drv
!# 1869 "ptoolkit.f90"
SUBROUTINE zrep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !!  Parallel matrix multiplication with replicated matrix
  !!
  !!  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !!
  !!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !!
  !!  WHERE  OP( X ) IS ONE OF
  !!
  !!     OP( X ) = X   OR   OP( X ) = X',
  !!
  !!  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !!  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !!
  !  written by Carlo Cavazzoni
  !
  USE laxlib_parallel_include
  implicit none
  include 'laxlib_kinds.fh'
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa
  !! specifies the form of op( A ) to be used in the matrix multiplication as
  !follows:
  !! 'N' or 'n',  op( A ) = A.
  !! 'T' or 't',  op( A ) = A**T.
  !! 'C' or 'c',  op( A ) = A**T. 
  CHARACTER(LEN=1), INTENT(IN) :: transb
  !! specifies the form of op( B ) to be used in the matrix multiplication as
  !follows:
  !! 'N' or 'n',  op( B ) = B.
  !! 'T' or 't',  op( B ) = B**T.
  !! 'C' or 'c',  op( B ) = B**T.
  INTEGER, INTENT(IN) :: m
  !! number of rows of the matrix A and C
  INTEGER, INTENT(IN) :: n
  !! number of columns of the matrix B and C
  INTEGER, INTENT(IN) :: k
  !! number of columns of A and rows of B
  REAL(DP), INTENT(IN) :: alpha
  !! scalar alpha
  REAL(DP), INTENT(IN) :: beta
  !! scalar beta
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of A
  INTEGER, INTENT(IN) :: ldb
  !! leading dimension of B
  INTEGER, INTENT(IN) :: ldc
  !! leading dimension of C
  COMPLEX(DP) :: a(lda,*)
  !! matrix A
  COMPLEX(DP) :: b(ldb,*)
  !! matrix B
  COMPLEX(DP) :: c(ldc,*)
  !! matrix C
  INTEGER, INTENT(IN) :: comm
  !! mpi communicator
  !
!# 2048 "ptoolkit.f90"
     !  if we are not compiling with __MPI this is equivalent to a blas call
!# 2050 "ptoolkit.f90"
     CALL zgemm(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)
!# 2054 "ptoolkit.f90"
  RETURN
!# 2056 "ptoolkit.f90"
END SUBROUTINE zrep_matmul_drv
!# 2058 "ptoolkit.f90"
!
!
!=----------------------------------------------------------------------------=!
!
!
!  Cannon's algorithms for parallel matrix multiplication
!  written by Carlo Cavazzoni
!  
!
!
!# 2069 "ptoolkit.f90"
SUBROUTINE sqr_dmm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   !
   !!  
   !!  Double precision parallel square matrix multiplication with Cannon's algorithm
   !!  performs one of the matrix-matrix operations
   !!
   !!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !!
   !!  where  op( x ) is one of
   !!
   !!     OP( X ) = X   OR   OP( X ) = X',
   !!
   !!  alpha and beta are scalars, and a, b and c are square matrices
   !!
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa
   !! specifies the form of op( A ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( A ) = A.
   !! 'T' or 't',  op( A ) = A**T.
   !! 'C' or 'c',  op( A ) = A**T. 
   CHARACTER(LEN=1), INTENT(IN) :: transb
   !! specifies the form of op( B ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( B ) = B.
   !! 'T' or 't',  op( B ) = B**T.
   !! 'C' or 'c',  op( B ) = B**T.
   INTEGER, INTENT(IN) :: n
   !! global dimension
   REAL(DP), INTENT(IN) :: alpha
   !! scalar alpha
   REAL(DP), INTENT(IN) :: beta
   !! scalar beta
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of A
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of B
   INTEGER, INTENT(IN) :: ldc
   !! leading dimension of C
   REAL(DP) :: a(lda,*)
   !! matrix A
   REAL(DP) :: b(ldb,*)
   !! matrix B
   REAL(DP) :: c(ldc,*)
   !! matrix C
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   !
   integer :: ierr
   integer :: np
   integer :: i, j, nr, nc, nb, iter, rowid, colid
   logical :: ta, tb
   INTEGER :: comm
   !
   !
   real(DP), allocatable :: bblk(:,:), ablk(:,:)
   !
!# 2140 "ptoolkit.f90"
   !
   CALL laxlib_intarray_to_desc(desc,idesc)
   !
   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF
!# 2151 "ptoolkit.f90"
   IF( n < 1 ) THEN
      RETURN
   END IF
!# 2155 "ptoolkit.f90"
   IF( desc%npr == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL dgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF
!# 2165 "ptoolkit.f90"
   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_mm_cannon ', ' works only with square processor mesh ', 1 )
   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc%npr
   comm  = desc%comm
   rowid = desc%myr
   colid = desc%myc
   !
   !  Retrieve the size of the local block
   !
   nr    = desc%nr 
   nc    = desc%nc 
   nb    = desc%nrcx
   !
!# 2186 "ptoolkit.f90"
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'T' .OR. TRANSA == 't' )
   tb = ( TRANSB == 'T' .OR. TRANSB == 't' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, 1.0_DP, c, ldc)
      !
   END DO
!# 2268 "ptoolkit.f90"
   deallocate( ablk, bblk )
   
   RETURN
!# 2272 "ptoolkit.f90"
CONTAINS
!# 2274 "ptoolkit.f90"
   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift length
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN
!# 2302 "ptoolkit.f90"
         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid
!# 2307 "ptoolkit.f90"
      ELSE IF( dir == 'S' ) THEN
!# 2309 "ptoolkit.f90"
         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid
!# 2314 "ptoolkit.f90"
      ELSE
!# 2316 "ptoolkit.f90"
         CALL lax_error__( ' sqr_mm_cannon ', ' unknown shift direction ', 1 )
!# 2318 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 2331 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_block
!# 2334 "ptoolkit.f90"
   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow
!# 2375 "ptoolkit.f90"
      ELSE
!# 2377 "ptoolkit.f90"
         CALL lax_error__( ' sqr_mm_cannon ', ' unknown shift_exch direction ', 1 )
!# 2379 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 2392 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_exch_block
!# 2395 "ptoolkit.f90"
END SUBROUTINE sqr_dmm_cannon_x
!# 2397 "ptoolkit.f90"
!=----------------------------------------------------------------------------=!
!# 2776 "ptoolkit.f90"
!=----------------------------------------------------------------------------=!
!# 2778 "ptoolkit.f90"
SUBROUTINE sqr_smm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   !
   !!  Single precision parallel square matrix multiplication with Cannon's algorithm
   !!  performs one of the matrix-matrix operations 
   !!
   !!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !!
   !!  where  op( x ) is one of
   !!
   !!     OP( X ) = X   OR   OP( X ) = X',
   !!
   !!  alpha and beta are scalars, and a, b and c are square matrices
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa
   !! specifies the form of op( A ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( A ) = A.
   !! 'T' or 't',  op( A ) = A**T.
   !! 'C' or 'c',  op( A ) = A**T. 
   CHARACTER(LEN=1), INTENT(IN) :: transb
   !! specifies the form of op( B ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( B ) = B.
   !! 'T' or 't',  op( B ) = B**T.
   !! 'C' or 'c',  op( B ) = B**T.
   INTEGER, INTENT(IN) :: n
   !! global dimension
   REAL(SP), INTENT(IN) :: alpha
   !! scalar alpha
   REAL(SP), INTENT(IN) :: beta
   !! scalar beta
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of A
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of B
   INTEGER, INTENT(IN) :: ldc
   !! leading dimension of C
   REAL(SP) :: a(lda,*)
   !! matrix A
   REAL(SP) :: b(ldb,*)
   !! matrix B
   REAL(SP) :: c(ldc,*)
   !! matrix C
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   !
   integer :: ierr
   integer :: np
   integer :: i, j, nr, nc, nb, iter, rowid, colid
   logical :: ta, tb
   INTEGER :: comm
   !
   !
   real(SP), allocatable :: bblk(:,:), ablk(:,:)
   !
!# 2847 "ptoolkit.f90"
   !
   CALL laxlib_intarray_to_desc(desc,idesc)
   !
   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF
!# 2858 "ptoolkit.f90"
   IF( n < 1 ) THEN
      RETURN
   END IF
!# 2862 "ptoolkit.f90"
   IF( desc%npr == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL sgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF
!# 2872 "ptoolkit.f90"
   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_smm_cannon ', ' works only with square processor mesh ', 1 )
   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc%npr
   comm  = desc%comm
   rowid = desc%myr
   colid = desc%myc
   !
   !  Retrieve the size of the local block
   !
   nr    = desc%nr 
   nc    = desc%nc 
   nb    = desc%nrcx
   !
!# 2893 "ptoolkit.f90"
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_SP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_SP
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = 0.0_SP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = 0.0_SP
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'T' .OR. TRANSA == 't' )
   tb = ( TRANSB == 'T' .OR. TRANSB == 't' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL sgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL sgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, 1.0_SP, c, ldc)
      !
   END DO
!# 2975 "ptoolkit.f90"
   deallocate( ablk, bblk )
   
   RETURN
!# 2979 "ptoolkit.f90"
CONTAINS
!# 2981 "ptoolkit.f90"
   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      REAL(SP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift length
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN
!# 3009 "ptoolkit.f90"
         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid
!# 3014 "ptoolkit.f90"
      ELSE IF( dir == 'S' ) THEN
!# 3016 "ptoolkit.f90"
         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid
!# 3021 "ptoolkit.f90"
      ELSE
!# 3023 "ptoolkit.f90"
         CALL lax_error__( ' sqr_smm_cannon ', ' unknown shift direction ', 1 )
!# 3025 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3038 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_block
!# 3041 "ptoolkit.f90"
   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      REAL(SP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow
!# 3082 "ptoolkit.f90"
      ELSE
!# 3084 "ptoolkit.f90"
         CALL lax_error__( ' sqr_smm_cannon ', ' unknown shift_exch direction ', 1 )
!# 3086 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3099 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_exch_block
!# 3102 "ptoolkit.f90"
END SUBROUTINE sqr_smm_cannon_x
!# 3105 "ptoolkit.f90"
!=----------------------------------------------------------------------------=!
!# 3107 "ptoolkit.f90"
SUBROUTINE sqr_zmm_cannon_x( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, idesc )
   !
   !!  Double precision complex (Z) parallel square matrix multiplication with Cannon's algorithm
   !!  performs one of the matrix-matrix operations 
   !!
   !!     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !!
   !!  where  op( x ) is one of
   !!
   !!     OP( X ) = X   OR   OP( X ) = X',
   !!
   !!  alpha and beta are scalars, and a, b and c are square matrices
   !
   USE laxlib_descriptor
   !
   USE laxlib_parallel_include
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa
   !! specifies the form of op( A ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( A ) = A.
   !! 'T' or 't',  op( A ) = A**T.
   !! 'C' or 'c',  op( A ) = A**T. 
   CHARACTER(LEN=1), INTENT(IN) :: transb
   !! specifies the form of op( B ) to be used in the matrix multiplication as
   !follows:
   !! 'N' or 'n',  op( B ) = B.
   !! 'T' or 't',  op( B ) = B**T.
   !! 'C' or 'c',  op( B ) = B**T.
   INTEGER, INTENT(IN) :: n
   !! global dimension
   COMPLEX(DP), INTENT(IN) :: alpha
   !! scalar alpha
   COMPLEX(DP), INTENT(IN) :: beta
   !! scalar beta
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of A
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of B
   INTEGER, INTENT(IN) :: ldc
   !! leading dimension of C
   COMPLEX(DP) :: a(lda,*)
   !! matrix A
   COMPLEX(DP) :: b(ldb,*)
   !! matrix B
   COMPLEX(DP) :: c(ldc,*)
   !! matrix C
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   !
   INTEGER :: ierr
   INTEGER :: np
   INTEGER :: i, j, nr, nc, nb, iter, rowid, colid
   LOGICAL :: ta, tb
   INTEGER :: comm
   !
   !
   COMPLEX(DP), ALLOCATABLE :: bblk(:,:), ablk(:,:)
   COMPLEX(DP) :: zone = ( 1.0_DP, 0.0_DP )
   COMPLEX(DP) :: zzero = ( 0.0_DP, 0.0_DP )
   !
!# 3178 "ptoolkit.f90"
   !
   CALL laxlib_intarray_to_desc(desc,idesc)
   !
   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF
!# 3189 "ptoolkit.f90"
   IF( n < 1 ) THEN
      RETURN
   END IF
!# 3193 "ptoolkit.f90"
   IF( desc%npr == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL zgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF
!# 3203 "ptoolkit.f90"
   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_zmm_cannon ', ' works only with square processor mesh ', 1 )
   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc%npr
   comm  = desc%comm
   rowid = desc%myr
   colid = desc%myc
   !
   !  Retrieve the size of the local block
   !
   nr    = desc%nr 
   nc    = desc%nc 
   nb    = desc%nrcx
   !
!# 3224 "ptoolkit.f90"
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'C' .OR. TRANSA == 'c' )
   tb = ( TRANSB == 'C' .OR. TRANSB == 'c' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, zone, c, ldc)
      !
   END DO
!# 3306 "ptoolkit.f90"
   deallocate( ablk, bblk )
   
   RETURN
!# 3310 "ptoolkit.f90"
CONTAINS
!# 3312 "ptoolkit.f90"
   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift length
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN
!# 3340 "ptoolkit.f90"
         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid
!# 3345 "ptoolkit.f90"
      ELSE IF( dir == 'S' ) THEN
!# 3347 "ptoolkit.f90"
         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid
!# 3352 "ptoolkit.f90"
      ELSE
!# 3354 "ptoolkit.f90"
         CALL lax_error__( ' sqr_zmm_cannon ', ' unknown shift direction ', 1 )
!# 3356 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3369 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_block
   !
   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow
!# 3413 "ptoolkit.f90"
      ELSE
!# 3415 "ptoolkit.f90"
         CALL lax_error__( ' sqr_zmm_cannon ', ' unknown shift_exch direction ', 1 )
!# 3417 "ptoolkit.f90"
      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3430 "ptoolkit.f90"
      RETURN
   END SUBROUTINE shift_exch_block
!# 3433 "ptoolkit.f90"
END SUBROUTINE sqr_zmm_cannon_x
!# 3435 "ptoolkit.f90"
!
!
!
!
!# 3440 "ptoolkit.f90"
SUBROUTINE sqr_tr_cannon_x( n, a, lda, b, ldb, idesc )
   !
   !!  Parallel square matrix transposition with Cannon's algorithm
   !!
   !
   USE laxlib_parallel_include
   IMPLICIT NONE
   !
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of A
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of B
   REAL(DP)            :: a(lda,*)
   !! matrix A
   REAL(DP)            :: b(ldb,*)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: i, j, nr, nc, nb
   INTEGER :: comm
   !
   REAL(DP), ALLOCATABLE :: ablk(:,:)
   !
!# 3476 "ptoolkit.f90"
   !
   IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
      RETURN
   END IF
!# 3481 "ptoolkit.f90"
   IF( n < 1 ) THEN
     RETURN
   END IF
!# 3485 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) == 1 ) THEN
      CALL mytranspose( a, lda, b, ldb, n, n )
      RETURN
   END IF
!# 3490 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) /= idesc(LAX_DESC_NPC) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' works only with square processor mesh ', 1 )
   IF( n /= idesc(LAX_DESC_N) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size n  ', 1 )
   IF( lda /= idesc(LAX_DESC_NRCX) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size lda  ', 1 )
   IF( ldb /= idesc(LAX_DESC_NRCX) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size ldb  ', 1 )
!# 3499 "ptoolkit.f90"
   comm = idesc(LAX_DESC_COMM)
!# 3501 "ptoolkit.f90"
   rowid = idesc(LAX_DESC_MYR)
   colid = idesc(LAX_DESC_MYC)
   np    = idesc(LAX_DESC_NPR)
   !
   !  Compute the size of the local block
   !
   nr = idesc(LAX_DESC_NR) 
   nc = idesc(LAX_DESC_NC) 
   nb = idesc(LAX_DESC_NRCX)
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   CALL exchange_block( ablk )
   !
!# 3535 "ptoolkit.f90"
   !
   DO j = 1, nr
      DO i = 1, nc
         b( j, i ) = ablk( i, j )
      END DO
   END DO
   !
   deallocate( ablk )
   
   RETURN
!# 3546 "ptoolkit.f90"
CONTAINS
!# 3548 "ptoolkit.f90"
   SUBROUTINE exchange_block( blk )
      !
      !   Block exchange ( transpose )
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      irdst = colid
      icdst = rowid
      irsrc = colid
      icsrc = rowid
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3574 "ptoolkit.f90"
      RETURN
   END SUBROUTINE
!# 3578 "ptoolkit.f90"
END SUBROUTINE
!# 3580 "ptoolkit.f90"
!
!# 3765 "ptoolkit.f90"
!
!# 3767 "ptoolkit.f90"
SUBROUTINE sqr_tr_cannon_sp_x( n, a, lda, b, ldb, idesc )
   !
   !! Parallel square matrix transposition with Cannon's algorithm
   !! single precision version
   !
   USE laxlib_parallel_include
   IMPLICIT NONE
   !
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of A
   INTEGER, INTENT(IN) :: ldb
   !! leading dimension of B
   REAL(SP)            :: a(lda,*)
   !! matrix A
   REAL(SP)            :: b(ldb,*)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: i, j, nr, nc, nb
   INTEGER :: comm
   !
   REAL(SP), ALLOCATABLE :: ablk(:,:)
   !
!# 3803 "ptoolkit.f90"
   !
   IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
      RETURN
   END IF
!# 3808 "ptoolkit.f90"
   IF( n < 1 ) THEN
     RETURN
   END IF
!# 3812 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) == 1 ) THEN
      CALL mytranspose_sp( a, lda, b, ldb, n, n )
      RETURN
   END IF
!# 3817 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) /= idesc(LAX_DESC_NPC) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' works only with square processor mesh ', 1 )
   IF( n /= idesc(LAX_DESC_N) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size n  ', 1 )
   IF( lda /= idesc(LAX_DESC_NRCX) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size lda  ', 1 )
   IF( ldb /= idesc(LAX_DESC_NRCX) ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size ldb  ', 1 )
!# 3826 "ptoolkit.f90"
   comm = idesc(LAX_DESC_COMM)
!# 3828 "ptoolkit.f90"
   rowid = idesc(LAX_DESC_MYR)
   colid = idesc(LAX_DESC_MYC)
   np    = idesc(LAX_DESC_NPR)
   !
   !  Compute the size of the local block
   !
   nr = idesc(LAX_DESC_NR) 
   nc = idesc(LAX_DESC_NC) 
   nb = idesc(LAX_DESC_NRCX)
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_SP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_SP
      END DO
   END DO
   !
   CALL exchange_block( ablk )
   !
!# 3862 "ptoolkit.f90"
   !
   DO j = 1, nr
      DO i = 1, nc
         b( j, i ) = ablk( i, j )
      END DO
   END DO
   !
   deallocate( ablk )
   
   RETURN
!# 3873 "ptoolkit.f90"
CONTAINS
!# 3875 "ptoolkit.f90"
   SUBROUTINE exchange_block( blk )
      !
      !   Block exchange ( transpose )
      !
      IMPLICIT NONE
      REAL(SP) :: blk( :, : )
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      irdst = colid
      icdst = rowid
      irsrc = colid
      icsrc = rowid
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
!# 3901 "ptoolkit.f90"
      RETURN
   END SUBROUTINE
!# 3905 "ptoolkit.f90"
END SUBROUTINE
!# 3909 "ptoolkit.f90"
SUBROUTINE redist_row2col_x( n, a, b, ldx, nx, idesc )
   !
   !!  redistribute a, array whose second dimension is distributed over processor row,
   !!  to obtain b, with the second dim. distributed over processor column 
   !
   !
   USE laxlib_parallel_include
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: ldx
   !! local rows
   INTEGER, INTENT(IN) :: nx
   !! local columns
   REAL(DP)            :: a(ldx,nx)
   !! matrix A
   REAL(DP)            :: b(ldx,nx)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: comm
   INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
   !
!# 3944 "ptoolkit.f90"
   !
   IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
      RETURN
   END IF
!# 3949 "ptoolkit.f90"
   IF( n < 1 ) THEN
     RETURN
   END IF
!# 3953 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) == 1 ) THEN
      b = a
      RETURN
   END IF
!# 3958 "ptoolkit.f90"
   IF( idesc(LAX_DESC_NPR) /= idesc(LAX_DESC_NPC) ) &
      CALL lax_error__( ' redist_row2col ', ' works only with square processor mesh ', 1 )
   IF( n /= idesc(LAX_DESC_N) ) &
      CALL lax_error__( ' redist_row2col ', ' inconsistent size n  ', 1 )
   IF( nx /= idesc(LAX_DESC_NRCX) ) &
      CALL lax_error__( ' redist_row2col ', ' inconsistent size lda  ', 1 )
!# 3965 "ptoolkit.f90"
   comm = idesc(LAX_DESC_COMM)
!# 3967 "ptoolkit.f90"
   rowid = idesc(LAX_DESC_MYR)
   colid = idesc(LAX_DESC_MYC)
   np    = idesc(LAX_DESC_NPR)
   !
   irdst = colid
   icdst = rowid
   irsrc = colid
   icsrc = rowid
   !
   CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
   CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
   !
!# 3991 "ptoolkit.f90"
   b = a
!# 3993 "ptoolkit.f90"
   !
   RETURN
!# 3996 "ptoolkit.f90"
END SUBROUTINE redist_row2col_x
!# 4114 "ptoolkit.f90"
!
!
!
!# 4118 "ptoolkit.f90"
SUBROUTINE cyc2blk_redist_x( n, a, lda, nca, b, ldb, ncb, idesc )
   !
   !!  Parallel square matrix redistribution. Double precision
   !!  A (input) is cyclically distributed by rows across processors
   !!  B (output) is distributed by block across 2D processors grid
   !!
   !
   USE laxlib_descriptor 
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! local rows of A
   INTEGER, INTENT(IN) :: nca
   !! local columns of A
   INTEGER, INTENT(IN) :: ldb
   !! local rows of B
   INTEGER, INTENT(IN) :: ncb
   !! local columns of B
   REAL(DP)            :: a(lda,nca)
   !! matrix A
   REAL(DP)            :: b(ldb,ncb)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
!# 4265 "ptoolkit.f90"
   b( 1:n, 1:n ) = a( 1:n, 1:n )   
!# 4269 "ptoolkit.f90"
   RETURN
!# 4271 "ptoolkit.f90"
CONTAINS
!# 4273 "ptoolkit.f90"
   SUBROUTINE check_sndbuf_index()
      CHARACTER(LEN=38), SAVE :: msg = ' check_sndbuf_index in cyc2blk_redist '
      IF( j  > SIZE(sndbuf,2) ) CALL lax_error__( msg, ' j > SIZE(sndbuf,2) ', ip+1 )
      IF( il > SIZE(sndbuf,1) ) CALL lax_error__( msg, ' il > SIZE(sndbuf,1) ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 < 1 ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 > lda ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip+1 )
      IF( jj < 1 ) CALL lax_error__( msg, ' jj < 1 ', ip+1 )
      IF( jj > n ) CALL lax_error__( msg, ' jj > n ', ip+1 )
      RETURN
   END SUBROUTINE check_sndbuf_index
!# 4284 "ptoolkit.f90"
   SUBROUTINE check_rcvbuf_index()
      CHARACTER(LEN=38), SAVE :: msg = ' check_rcvbuf_index in cyc2blk_redist '
      IF( i > ldb ) CALL lax_error__( msg, ' i > ldb ', ip+1 )
      IF( j > ldb ) CALL lax_error__( msg, ' j > ldb ', ip+1 )
      IF( j > nb  ) CALL lax_error__( msg, ' j > nb  ', ip+1 )
      IF( il > SIZE( rcvbuf, 1 ) ) CALL lax_error__( msg, ' il too large ', ip+1 )
      RETURN
   END SUBROUTINE check_rcvbuf_index
!# 4293 "ptoolkit.f90"
END SUBROUTINE cyc2blk_redist_x
!# 4296 "ptoolkit.f90"
SUBROUTINE cyc2blk_zredist_x( n, a, lda, nca, b, ldb, ncb, idesc )
   !
   !! Parallel square matrix redistribution. Double precision complex (Z)
   !! A (input) is cyclically distributed by rows across processors
   !! B (output) is distributed by block across 2D processors grid
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! local rows of A
   INTEGER, INTENT(IN) :: nca
   !! local columns of A
   INTEGER, INTENT(IN) :: ldb
   !! local rows of B
   INTEGER, INTENT(IN) :: ncb
   !! local columns of B
   COMPLEX(DP)            :: a(lda,nca)
   !! matrix A
   COMPLEX(DP)            :: b(ldb,ncb)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   !
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
!# 4435 "ptoolkit.f90"
   b( 1:n, 1:n ) = a( 1:n, 1:n )   
!# 4439 "ptoolkit.f90"
   RETURN
!# 4441 "ptoolkit.f90"
CONTAINS
!# 4443 "ptoolkit.f90"
   SUBROUTINE check_sndbuf_index()
      CHARACTER(LEN=40), SAVE :: msg = ' check_sndbuf_index in cyc2blk_zredist  '
      IF( j  > SIZE(sndbuf,2) ) CALL lax_error__( msg, ' j > SIZE(sndbuf,2) ', ip+1 )
      IF( il > SIZE(sndbuf,1) ) CALL lax_error__( msg, ' il > SIZE(sndbuf,1) ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 < 1 ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 > SIZE(a,1) ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip+1 )
      IF( jj < 1 ) CALL lax_error__( msg, ' jj < 1 ', ip+1 )
      IF( jj > n ) CALL lax_error__( msg, ' jj > n ', ip+1 )
      RETURN
   END SUBROUTINE check_sndbuf_index
!# 4454 "ptoolkit.f90"
   SUBROUTINE check_rcvbuf_index()
      CHARACTER(LEN=40), SAVE :: msg = ' check_rcvbuf_index in cyc2blk_zredist  '
      IF( i > ldb ) CALL lax_error__( msg, ' i > ldb ', ip+1 )
      IF( j > ldb ) CALL lax_error__( msg, ' j > ldb ', ip+1 )
      IF( j > nb  ) CALL lax_error__( msg, ' j > nb  ', ip+1 )
      IF( il > SIZE( rcvbuf, 1 ) ) CALL lax_error__( msg, ' il too large ', ip+1 )
      RETURN
   END SUBROUTINE check_rcvbuf_index
!# 4463 "ptoolkit.f90"
END SUBROUTINE cyc2blk_zredist_x
!# 4468 "ptoolkit.f90"
SUBROUTINE blk2cyc_redist_x( n, a, lda, nca, b, ldb, ncb, idesc )
   !
   !!  Parallel square matrix redistribution. Double precision.
   !!  A (output) is cyclically distributed by rows across processors
   !!  B (input) is distributed by block across 2D processors grid
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! local rows of A
   INTEGER, INTENT(IN) :: nca
   !! local columns of A
   INTEGER, INTENT(IN) :: ldb
   !! local rows of B
   INTEGER, INTENT(IN) :: ncb
   !! local columns of B
   REAL(DP)            :: a(lda,nca)
   !! matrix A
   REAL(DP)            :: b(ldb,ncb)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   REAL(DP), allocatable :: rcvbuf(:,:,:)
   REAL(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
!# 4601 "ptoolkit.f90"
   a( 1:n, 1:n ) = b( 1:n, 1:n )   
!# 4605 "ptoolkit.f90"
   RETURN
!# 4607 "ptoolkit.f90"
END SUBROUTINE blk2cyc_redist_x
!# 4610 "ptoolkit.f90"
SUBROUTINE blk2cyc_zredist_x( n, a, lda, nca, b, ldb, ncb, idesc )
   !
   !!  Parallel square matrix redistribution. Double precision complex (Z)
   !!  A (output) is cyclically distributed by rows across processors
   !!  B (input) is distributed by block across 2D processors grid
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   INTEGER, INTENT(IN) :: n
   !! global dimension
   INTEGER, INTENT(IN) :: lda
   !! local rows of A
   INTEGER, INTENT(IN) :: nca
   !! local columns of A
   INTEGER, INTENT(IN) :: ldb
   !! local rows of B
   INTEGER, INTENT(IN) :: ncb
   !! local columns of B
   COMPLEX(DP)            :: a(lda,nca)
   !! matrix A
   COMPLEX(DP)            :: b(ldb,ncb)
   !! matrix B
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   !
   TYPE(la_descriptor) :: desc
   !
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
!# 4744 "ptoolkit.f90"
   a( 1:n, 1:n ) = b( 1:n, 1:n )   
!# 4748 "ptoolkit.f90"
   RETURN
!# 4750 "ptoolkit.f90"
END SUBROUTINE blk2cyc_zredist_x
!
!
!
!  Double Complex and Double Precision Cholesky Factorization of
!  an Hermitan/Symmetric block distributed matrix
!  written by Carlo Cavazzoni
!
!
!# 4760 "ptoolkit.f90"
SUBROUTINE laxlib_pzpotrf_x( sll, ldx, n, idesc )
   !
   !! Double precision Complex (Z) Cholesky Factorization of
   !! an Hermitan/Symmetric block distributed matrix
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   implicit none
   !
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   !
   integer :: n
   !! global dimension
   integer :: ldx
   !! leading dimension of sll   
   integer, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   complex(DP) :: sll( ldx, ldx )
   !! matrix sll
   !
   real(DP)  :: one, zero
   complex(DP) :: cone, czero
   integer :: myrow, mycol, ierr
   integer :: jb, info, ib, kb
   integer :: jnr, jir, jic, jnc
   integer :: inr, iir, iic, inc
   integer :: knr, kir, kic, knc
   integer :: nr, nc
   integer :: rcomm, ccomm, color, key, myid, np
   complex(DP), allocatable :: ssnd( :, : ), srcv( :, : )
   TYPE(la_descriptor) :: desc
!# 4794 "ptoolkit.f90"
   one   = 1.0_DP
   cone  = 1.0_DP
   zero  = 0.0_DP
   czero = 0.0_DP
!# 4993 "ptoolkit.f90"
   CALL ZPOTRF( 'L', n, sll, ldx, info )
!# 4995 "ptoolkit.f90"
   IF( info /= 0 ) &
      CALL lax_error__( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
!# 5000 "ptoolkit.f90"
   return
END SUBROUTINE laxlib_pzpotrf_x
!# 5003 "ptoolkit.f90"
!  now the Double Precision subroutine
!# 5005 "ptoolkit.f90"
SUBROUTINE laxlib_pdpotrf_x( sll, ldx, n, idesc )
   !
   !! Double precision Cholesky Factorization of
   !! an Hermitan/Symmetric block distributed matrix
   !
   USE laxlib_descriptor
   USE laxlib_parallel_include
   !
   implicit none
   !
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   !
   integer  :: n
   !! global dimension
   integer  :: ldx
   !! leading dimension of sll
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   REAL(DP) :: sll( ldx, ldx )
   !! matrix sll
   REAL(DP) :: one, zero
   integer  :: myrow, mycol, ierr
   integer  :: jb, info, ib, kb
   integer  :: jnr, jir, jic, jnc
   integer  :: inr, iir, iic, inc
   integer  :: knr, kir, kic, knc
   integer  :: nr, nc
   integer  :: rcomm, ccomm, color, key, myid, np
   REAL(DP), ALLOCATABLE :: ssnd( :, : ), srcv( :, : )
   TYPE(la_descriptor) :: desc
!# 5037 "ptoolkit.f90"
   one   = 1.0_DP
   zero  = 0.0_DP
!# 5231 "ptoolkit.f90"
   CALL DPOTRF( 'L', n, sll, ldx, info )
!# 5233 "ptoolkit.f90"
   IF( info /= 0 ) &
      CALL lax_error__( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
!# 5238 "ptoolkit.f90"
   return
END SUBROUTINE laxlib_pdpotrf_x
!# 5241 "ptoolkit.f90"
!
!
!
!
!# 5246 "ptoolkit.f90"
SUBROUTINE laxlib_pztrtri_x ( sll, ldx, n, idesc )
    ! 
    !! pztrtri computes the parallel inversion of a lower triangular matrix 
    !! distribuited among the processes using a 2-D block partitioning. 
    !! The algorithm is based on the schema below and executes the model 
    !! recursively to each column C2 under the diagonal.     
    !!
    !!     |-------|-------|      |--------------------|--------------------|
    !!     |   A1  |   0   |      |   C1 = trtri(A1)   |          0         |
    !! A = |-------|-------|  C = |--------------------|--------------------|
    !!     |   A2  |   A3  |      | C2 = -C3 * A2 * C1 |   C3 = trtri(A3)   | 
    !!     |-------|-------|      |--------------------|--------------------|
    !!
    !! The recursive steps of multiplication (C2 = -C3 * A2 * C1) is based on the Cannon's algorithms 
    !! for parallel matrix multiplication and is done with BLACS(dgemm)
    !!
    !!
    !
    ! Arguments
    ! ============
    !
    ! sll   = local block of data
    ! ldx   = leading dimension of one block
    ! n     = size of the global array diributed among the blocks
    ! desc  = descriptor of the matrix distribution
    !
    !
    !  written by Ivan Girotto
    !
!# 5276 "ptoolkit.f90"
    USE laxlib_descriptor
    USE laxlib_parallel_include
!# 5279 "ptoolkit.f90"
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
!# 5283 "ptoolkit.f90"
    INTEGER, INTENT(IN)  :: n
    !! global dimension
    INTEGER, INTENT(IN)  :: ldx
    !! leading dimension of sll
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! integer laxlib descriptor
    COMPLEX(DP), INTENT( INOUT ) :: sll( ldx, ldx )
    !! matrix sll
!# 5292 "ptoolkit.f90"
    COMPLEX(DP), PARAMETER :: ONE = (1.0_DP, 0.0_DP)
    COMPLEX(DP), PARAMETER :: ZERO = (0.0_DP, 0.0_DP)
!# 5298 "ptoolkit.f90"
    INTEGER :: req(2), ierr, col_comm
    INTEGER :: send, recv, group_rank, group_size
    INTEGER :: myrow, mycol, np, myid, comm
!# 5302 "ptoolkit.f90"
    ! counters
    INTEGER :: k, i, j, count, step_count, shiftcount, cicle 
    INTEGER :: C3dim   ! Dimension of submatrix B
    INTEGER :: nc, nr ! Local dimension of block
    INTEGER :: info, sup_recv
    INTEGER :: idrowref, idcolref, idref, idrecv 
!# 5309 "ptoolkit.f90"
    ! B and BUF_RECV are used to overload the computation of matrix multiplication and the shift of the blocks
    COMPLEX(DP), ALLOCATABLE, DIMENSION( :, : ) :: B, C, BUF_RECV 
    COMPLEX(DP) :: first
    TYPE(la_descriptor) :: desc
!# 5314 "ptoolkit.f90"
    CALL laxlib_intarray_to_desc(desc,idesc)
    myrow = desc%myr
    mycol = desc%myc
    myid  = desc%mype
    np    = desc%npr
    comm  = desc%comm
!# 5321 "ptoolkit.f90"
    IF( desc%npr /= desc%npc ) THEN
       CALL lax_error__( ' pztrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc%nrcx ) THEN
       CALL lax_error__( ' pztrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF
!# 5328 "ptoolkit.f90"
    nr = desc%nr
    nc = desc%nc
!# 5331 "ptoolkit.f90"
    !  clear elements outside local meaningful block nr*nc
!# 5333 "ptoolkit.f90"
    DO j = nc+1, ldx
       DO i = 1, ldx
          sll( i, j ) = zero
       END DO
    END DO
    DO j = 1, ldx
       DO i = nr+1, ldx
          sll( i, j ) = zero
       END DO
    END DO
!# 5569 "ptoolkit.f90"
    CALL compute_ztrtri()
!# 5573 "ptoolkit.f90"
  CONTAINS
!# 5575 "ptoolkit.f90"
     SUBROUTINE compute_ztrtri()
       !
       !  clear the upper triangle (excluding diagonal terms) and
       !
       DO j = 1, ldx
          DO i = 1, j-1
             sll ( i, j ) = zero
          END DO
       END DO
       !
       CALL ztrtri( 'L', 'N', nr, sll, ldx, info )
       !
       IF( info /= 0 ) THEN
          CALL lax_error__( ' pztrtri ', ' problem in the local inversion ', info )
       END IF
       !
     END SUBROUTINE compute_ztrtri
!# 5594 "ptoolkit.f90"
     INTEGER FUNCTION shift ( idref, id, pos, size, dir )
!# 5596 "ptoolkit.f90"
       IMPLICIT NONE
   
       INTEGER :: idref, id, pos, size
       CHARACTER ( LEN = 1 ) :: dir
   
       IF( ( dir == 'E' ) .OR. ( dir == 'S' ) ) THEN
          shift = idref + MOD ( ( id - idref ) + pos, size )
       ELSE IF( ( dir == 'W' ) .OR. ( dir == 'N' ) ) THEN
          shift = idref + MOD ( ( id - idref ) - pos + size, size )
       ELSE
          shift = -1
       END IF
   
       RETURN
!# 5611 "ptoolkit.f90"
     END FUNCTION shift
!# 5613 "ptoolkit.f90"
END SUBROUTINE laxlib_pztrtri_x
!# 5615 "ptoolkit.f90"
!  now the Double Precision subroutine
!# 5617 "ptoolkit.f90"
SUBROUTINE laxlib_pdtrtri_x ( sll, ldx, n, idesc )
    
    !
    !!
    !! pdtrtri computes the parallel inversion of a lower triangular matrix 
    !! distribuited among the processes using a 2-D block partitioning. 
    !! The algorithm is based on the schema below and executes the model 
    !! recursively to each column C2 under the diagonal.     
    !!
    !!     |-------|-------|      |--------------------|--------------------|
    !!     |   A1  |   0   |      |   C1 = trtri(A1)   |          0         |
    !! A = |-------|-------|  C = |--------------------|--------------------|
    !!     |   A2  |   A3  |      | C2 = -C3 * A2 * C1 |   C3 = trtri(A3)   | 
    !!     |-------|-------|      |--------------------|--------------------|
    !!
    !! The recursive steps of multiplication (C2 = -C3 * A2 * C1) is based on the Cannon's algorithms 
    !! for parallel matrix multiplication and is done with BLACS(dgemm)
    !!
    !
    ! Arguments
    ! ============
    !
    ! sll   = local block of data
    ! ldx   = leading dimension of one block
    ! n     = size of the global array diributed among the blocks
    ! idesc  = descriptor of the matrix distribution
    !
    !
    !  written by Ivan Girotto
    !
!# 5648 "ptoolkit.f90"
    USE laxlib_descriptor
    USE laxlib_parallel_include
!# 5651 "ptoolkit.f90"
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
!# 5655 "ptoolkit.f90"
    INTEGER, INTENT(IN)  :: n
    !! global dimension
    INTEGER, INTENT(IN)  :: ldx
    !! leading dimension of sll
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! integer laxlib descriptor
    REAL(DP), INTENT( INOUT ) :: sll( ldx, ldx )
    !! matrix sll
!# 5664 "ptoolkit.f90"
    REAL(DP), PARAMETER :: ONE = 1.0_DP
    REAL(DP), PARAMETER :: ZERO = 0.0_DP
!# 5670 "ptoolkit.f90"
    INTEGER :: req(2), ierr, col_comm
    INTEGER :: send, recv, group_rank, group_size
    INTEGER :: myrow, mycol, np, myid, comm
!# 5674 "ptoolkit.f90"
    ! counters
    INTEGER :: k, i, j, count, step_count, shiftcount, cicle 
    INTEGER :: C3dim   ! Dimension of submatrix B
    INTEGER :: nc, nr ! Local dimension of block
    INTEGER :: info, sup_recv
    INTEGER :: idrowref, idcolref, idref, idrecv 
!# 5681 "ptoolkit.f90"
    ! B and BUF_RECV are used to overload the computation of matrix multiplication and the shift of the blocks
    REAL(DP), ALLOCATABLE, DIMENSION( :, : ) :: B, C, BUF_RECV 
    REAL(DP) :: first
    TYPE(la_descriptor) :: desc
!# 5686 "ptoolkit.f90"
    CALL laxlib_intarray_to_desc(desc,idesc)
!# 5688 "ptoolkit.f90"
    myrow = desc%myr
    mycol = desc%myc
    myid  = desc%mype
    np    = desc%npr
    comm  = desc%comm
!# 5694 "ptoolkit.f90"
    IF( desc%npr /= desc%npc ) THEN
       CALL lax_error__( ' pdtrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc%nrcx ) THEN
       CALL lax_error__( ' pdtrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF
!# 5701 "ptoolkit.f90"
    nr = desc%nr
    nc = desc%nc
!# 5704 "ptoolkit.f90"
    !  clear elements outside local meaningful block nr*nc
!# 5706 "ptoolkit.f90"
    DO j = nc+1, ldx
       DO i = 1, ldx
          sll( i, j ) = zero
       END DO
    END DO
    DO j = 1, ldx
       DO i = nr+1, ldx
          sll( i, j ) = zero
       END DO
    END DO
!# 5946 "ptoolkit.f90"
    CALL compute_dtrtri()
!# 5950 "ptoolkit.f90"
  CONTAINS
!# 5952 "ptoolkit.f90"
     SUBROUTINE compute_dtrtri()
       !
       !  clear the upper triangle (excluding diagonal terms) and
       !
       DO j = 1, ldx
          DO i = 1, j-1
             sll ( i, j ) = zero
          END DO
       END DO
       !
       CALL dtrtri( 'L', 'N', nr, sll, ldx, info )
       !
       IF( info /= 0 ) THEN
          CALL lax_error__( ' pdtrtri ', ' problem in the local inversion ', info )
       END IF
       !
     END SUBROUTINE compute_dtrtri
!# 5971 "ptoolkit.f90"
     INTEGER FUNCTION shift ( idref, id, pos, size, dir )
!# 5973 "ptoolkit.f90"
       IMPLICIT NONE
   
       INTEGER :: idref, id, pos, size
       CHARACTER ( LEN = 1 ) :: dir
   
       IF( ( dir == 'E' ) .OR. ( dir == 'S' ) ) THEN
          shift = idref + MOD ( ( id - idref ) + pos, size )
       ELSE IF( ( dir == 'W' ) .OR. ( dir == 'N' ) ) THEN
          shift = idref + MOD ( ( id - idref ) - pos + size, size )
       ELSE
          shift = -1
       END IF
   
       RETURN
!# 5988 "ptoolkit.f90"
     END FUNCTION shift
!# 5990 "ptoolkit.f90"
END SUBROUTINE laxlib_pdtrtri_x
!# 5994 "ptoolkit.f90"
SUBROUTINE laxlib_pdsyevd_x( tv, n, idesc, hh, ldh, e )
   !
   !! Parallel version of the HOUSEHOLDER tridiagonalization Algorithm for simmetric matrix.
   !! double precision version
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   include 'laxlib_low.fh'
   LOGICAL, INTENT(IN) :: tv
   !! if true compute eigenvalues and eigenvectors (not used)
   INTEGER, INTENT(IN) :: n
   !! global dimension 
   INTEGER, INTENT(IN) :: ldh
   !! leading dimension of hh
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   REAL(DP) :: hh( ldh, ldh )
   !! matrix to be diagonalized and output eigenvectors
   REAL(DP) :: e( n )
   !! eigenvalues
!# 6016 "ptoolkit.f90"
   INTEGER :: nrlx, nrl, nproc
   REAL(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
   CHARACTER :: jobv
!# 6020 "ptoolkit.f90"
   nrl   = idesc(LAX_DESC_NRL) 
   nrlx  = idesc(LAX_DESC_NRLX) 
   nproc = idesc(LAX_DESC_NPC) * idesc(LAX_DESC_NPR)
!# 6024 "ptoolkit.f90"
   ALLOCATE( diag( nrlx, n ) )
   ALLOCATE( vv( nrlx, n ) )
!# 6027 "ptoolkit.f90"
   jobv = 'N'
   IF( tv ) jobv = 'V'
   !
   !  Redistribute matrix "hh" into "diag",  
   !  matrix "hh" is block distributed, matrix diag is cyclic distributed
   CALL blk2cyc_redist( n, diag, nrlx, n, hh, ldh, ldh, idesc )
   !
   CALL dspev_drv( jobv, diag, nrlx, e, vv, nrlx, nrl, n, &
        nproc, idesc(LAX_DESC_MYPE), idesc(LAX_DESC_COMM) )
   !
   IF( tv ) CALL cyc2blk_redist( n, vv, nrlx, n, hh, ldh, ldh, idesc )
   !
   DEALLOCATE( vv )
   DEALLOCATE( diag )
!# 6042 "ptoolkit.f90"
   RETURN
END SUBROUTINE
!# 6047 "ptoolkit.f90"
SUBROUTINE laxlib_pzheevd_x( tv, n, idesc, hh, ldh, e )
   !
   !! Parallel version of the HOUSEHOLDER tridiagonalization Algorithm for simmetric matrix.
   !! double precision complex(Z) version
   !
   IMPLICIT NONE
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   include 'laxlib_low.fh'
   LOGICAL, INTENT(IN) :: tv
   !! if true compute eigenvalues and eigenvectors (not used)
   INTEGER, INTENT(IN) :: n
   !! global dimensio of matrix 
   INTEGER, INTENT(IN) :: ldh
   !! leading dimension of hh
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor
   COMPLEX(DP) :: hh( ldh, ldh )
   !! matrix to be diagonalized and output eigenvectors
   REAL(DP) :: e( n )
   !! eigenvalues
!# 6069 "ptoolkit.f90"
   INTEGER :: nrlx, nrl
   COMPLEX(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
   CHARACTER :: jobv
!# 6073 "ptoolkit.f90"
   nrl  = idesc(LAX_DESC_NRL)
   nrlx = idesc(LAX_DESC_NRLX)
   !
   ALLOCATE( diag( nrlx, n ) )
   ALLOCATE( vv( nrlx, n ) )
   !
   jobv = 'N'
   IF( tv ) jobv = 'V'
!# 6082 "ptoolkit.f90"
   CALL blk2cyc_redist( n, diag, nrlx, n, hh, ldh, ldh, idesc )
   !
   CALL zhpev_drv( jobv, diag, nrlx, e, vv, nrlx, nrl, n, &
        idesc(LAX_DESC_NPC) * idesc(LAX_DESC_NPR), idesc(LAX_DESC_MYPE), idesc(LAX_DESC_COMM) )
   !
   if( tv ) CALL cyc2blk_redist( n, vv, nrlx, n, hh, ldh, ldh, idesc )
   !
   DEALLOCATE( vv ) 
   DEALLOCATE( diag )
!# 6092 "ptoolkit.f90"
   RETURN
END SUBROUTINE
!# 6097 "ptoolkit.f90"
SUBROUTINE sqr_dsetmat_x( what, n, alpha, a, lda, idesc )
   !
   !!
   !!  Set the double precision values of a square distributed matrix 
   !!
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   CHARACTER(LEN=1), INTENT(IN) :: what
   !! 'A' set all the values of "a" equal to alpha
   !! 'U' set the values in the upper triangle of "a" equal to alpha
   !! 'L' set the values in the lower triangle of "a" equal to alpha
   !! 'D' set the values in the diagonal of "a" equal to alpha
   INTEGER, INTENT(IN) :: n
   !! global dimension of the matrix
   REAL(DP), INTENT(IN) :: alpha
   !! value to be assigned to elements of "a"
   INTEGER, INTENT(IN) :: lda
   !!! leading dimension of a
   REAL(DP) :: a(lda,*)
   !! matrix whose values have to be set
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor of matrix a
!# 6124 "ptoolkit.f90"
   INTEGER :: i, j
!# 6126 "ptoolkit.f90"
   IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF
!# 6134 "ptoolkit.f90"
   SELECT CASE( what )
     CASE( 'U', 'u' )
        IF( idesc(LAX_DESC_MYC) > idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, j - 1
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'L', 'l' )
        IF( idesc(LAX_DESC_MYC) < idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = j + 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'D', 'd' )
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO i = 1, idesc(LAX_DESC_NR)
              a( i, i ) = alpha
           END DO
        END IF
     CASE DEFAULT
        DO j = 1, idesc(LAX_DESC_NC)
           DO i = 1, idesc(LAX_DESC_NR)
              a( i, j ) = alpha
           END DO
        END DO
   END SELECT
   !
   RETURN
END SUBROUTINE sqr_dsetmat_x
!# 6181 "ptoolkit.f90"
SUBROUTINE sqr_zsetmat_x( what, n, alpha, a, lda, idesc )
   !
   !!  Set the double precision complex(Z) values of a square distributed matrix 
   !
   IMPLICIT NONE
   !
   include 'laxlib_kinds.fh'
   include 'laxlib_param.fh'
   !
   CHARACTER(LEN=1), INTENT(IN) :: what
   !! 'A' set all the values of "a" equal to alpha
   !! 'U' set the values in the upper triangle of "a" equal to alpha
   !! 'L' set the values in the lower triangle of "a" equal to alpha
   !! 'D' set the values in the diagonal of "a" equal to alpha
   !! 'H' clear the imaginary part of the diagonal of "a" 
   INTEGER, INTENT(IN) :: n
   !! global dimension of the matrix
   COMPLEX(DP), INTENT(IN) :: alpha
   !! value to be assigned to elements of "a"
   INTEGER, INTENT(IN) :: lda
   !! leading dimension of a
   COMPLEX(DP) :: a(lda,*)
   !! matrix whose values have to be set
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   !! integer laxlib descriptor of matrix a
!# 6207 "ptoolkit.f90"
   INTEGER :: i, j
!# 6209 "ptoolkit.f90"
   IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF
!# 6217 "ptoolkit.f90"
   SELECT CASE( what )
     CASE( 'U', 'u' )
        IF( idesc(LAX_DESC_MYC) > idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, j - 1
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'L', 'l' )
        IF( idesc(LAX_DESC_MYC) < idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO j = 1, idesc(LAX_DESC_NC)
              DO i = j + 1, idesc(LAX_DESC_NR)
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'D', 'd' )
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO i = 1, idesc(LAX_DESC_NR)
              a( i, i ) = alpha
           END DO
        END IF
     CASE( 'H', 'h' )
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           DO i = 1, idesc(LAX_DESC_NR)
              a( i, i ) = CMPLX( DBLE( a(i,i) ), 0_DP, KIND=DP )
           END DO
        END IF
     CASE DEFAULT
        DO j = 1, idesc(LAX_DESC_NC)
           DO i = 1, idesc(LAX_DESC_NR)
              a( i, j ) = alpha
           END DO
        END DO
   END SELECT
   !
   RETURN
END SUBROUTINE sqr_zsetmat_x
!# 6269 "ptoolkit.f90"
!------------------------------------------------------------------------
    SUBROUTINE distribute_lambda_x( lambda_repl, lambda_dist, idesc )
!------------------------------------------------------------------------
       IMPLICIT NONE
       include 'laxlib_kinds.fh'
       include 'laxlib_param.fh'
       REAL(DP), INTENT(IN)  :: lambda_repl(:,:)
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
       INTEGER :: i, j, ic, ir
       IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
          ir = idesc(LAX_DESC_IR)
          ic = idesc(LAX_DESC_IC)
          DO j = 1, idesc(LAX_DESC_NC)
             DO i = 1, idesc(LAX_DESC_NR)
                lambda_dist( i, j ) = lambda_repl( i + ir - 1, j + ic - 1 )
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_lambda_x
    !
!------------------------------------------------------------------------
    SUBROUTINE collect_lambda_x( lambda_repl, lambda_dist, idesc )
!------------------------------------------------------------------------
       USE laxlib_processors_grid,   ONLY: ortho_parent_comm
       USE laxlib_parallel_include
       IMPLICIT NONE
       include 'laxlib_kinds.fh'
       include 'laxlib_param.fh'
       REAL(DP), INTENT(OUT) :: lambda_repl(:,:)
       REAL(DP), INTENT(IN)  :: lambda_dist(:,:)
       INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
       INTEGER :: i, j, ic, ir, ierr
       lambda_repl = 0.0d0
       IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
          ir = idesc(LAX_DESC_IR)
          ic = idesc(LAX_DESC_IC)
          DO j = 1, idesc(LAX_DESC_NC)
             DO i = 1, idesc(LAX_DESC_NR)
                lambda_repl( i + ir - 1, j + ic - 1 ) = lambda_dist( i, j )
             END DO
          END DO
       END IF
!# 6317 "ptoolkit.f90"
       RETURN
    END SUBROUTINE collect_lambda_x
!# 6320 "ptoolkit.f90"
!------------------------------------------------------------------------
    SUBROUTINE collect_zmat_x( zmat_repl, zmat_dist, idesc )
!------------------------------------------------------------------------
       USE laxlib_processors_grid,   ONLY: ortho_parent_comm
       USE laxlib_parallel_include
       IMPLICIT NONE
       include 'laxlib_kinds.fh'
       include 'laxlib_param.fh'
       REAL(DP), INTENT(OUT) :: zmat_repl(:,:)
       REAL(DP), INTENT(IN)  :: zmat_dist(:,:)
       INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
       INTEGER :: i, ii, j, me, np, nrl, ierr
       zmat_repl = 0.0d0
       me = idesc(LAX_DESC_MYPE)
       np = idesc(LAX_DESC_NPC) * idesc(LAX_DESC_NPR)
       nrl = idesc(LAX_DESC_NRL)
       IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
          DO j = 1, idesc(LAX_DESC_N)
             ii = me + 1
             DO i = 1, nrl
                zmat_repl( ii, j ) = zmat_dist( i, j )
                ii = ii + np
             END DO
          END DO
       END IF
!# 6349 "ptoolkit.f90"
       RETURN
    END SUBROUTINE collect_zmat_x
!# 6352 "ptoolkit.f90"
!------------------------------------------------------------------------
    SUBROUTINE setval_lambda_x( lambda_dist, i, j, val, idesc )
!------------------------------------------------------------------------
       IMPLICIT NONE
       include 'laxlib_kinds.fh'
       include 'laxlib_param.fh'
       REAL(DP), INTENT(OUT) :: lambda_dist(:,:)
       INTEGER,  INTENT(IN)  :: i, j
       REAL(DP), INTENT(IN)  :: val
       INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
       INTEGER :: ir, ic
       IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
          ir = idesc(LAX_DESC_IR)
          ic = idesc(LAX_DESC_IC)
          IF( ( i >= ir ) .AND. ( i - ir + 1 <= idesc(LAX_DESC_NR) ) ) THEN
             IF( ( j >= ic ) .AND. ( j - ic + 1 <= idesc(LAX_DESC_NC) ) ) THEN
                lambda_dist( i - ir + 1, j - ic + 1 ) = val
             END IF
          END IF
       END IF
       RETURN
    END SUBROUTINE setval_lambda_x
!# 6375 "ptoolkit.f90"
!------------------------------------------------------------------------
    SUBROUTINE distribute_zmat_x( zmat_repl, zmat_dist, idesc )
!------------------------------------------------------------------------
       IMPLICIT NONE
       include 'laxlib_kinds.fh'
       include 'laxlib_param.fh'
       REAL(DP), INTENT(IN)  :: zmat_repl(:,:)
       REAL(DP), INTENT(OUT) :: zmat_dist(:,:)
       INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
       INTEGER :: i, ii, j, me, np
       me = idesc(LAX_DESC_MYPE)
       np = idesc(LAX_DESC_NPC) * idesc(LAX_DESC_NPR)
       IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
          DO j = 1, idesc(LAX_DESC_N)
             ii = me + 1
             DO i = 1, idesc(LAX_DESC_NRL)
                zmat_dist( i, j ) = zmat_repl( ii, j )
                ii = ii + np
             END DO
          END DO
       END IF
       RETURN
    END SUBROUTINE distribute_zmat_x
    !
