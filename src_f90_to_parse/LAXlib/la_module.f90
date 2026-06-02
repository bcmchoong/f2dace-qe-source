!# 1 "la_module.f90"
MODULE LAXlib
!# 5 "la_module.f90"
  IMPLICIT NONE
  !
  INTERFACE diaghg
     MODULE PROCEDURE cdiaghg_cpu_, rdiaghg_cpu_
!# 12 "la_module.f90"
  END INTERFACE
  !
  INTERFACE pdiaghg
     MODULE PROCEDURE pcdiaghg_, prdiaghg_
!# 19 "la_module.f90"
  END INTERFACE
  !
  INTERFACE diagh
     MODULE PROCEDURE cdiagh_cpu_, rdiagh_cpu_
  END INTERFACE
  !
  INTERFACE pdiagh
     MODULE PROCEDURE pcdiagh_, prdiagh_
  END INTERFACE
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE cdiaghg_cpu_( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! complex matrices version.
    !! On output both matrix are unchanged.
    !!
    !! LAPACK version - uses both ZHEGV and ZHEGVX
    !
!# 46 "la_module.f90"
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,n)
    !! matrix to be diagonalized
    COMPLEX(DP), INTENT(INOUT) :: s(ldh,n)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    LOGICAL, OPTIONAL ::  offload
    !! optionally solve the eigenvalue problem on the GPU
    LOGICAL :: loffload
      !
!# 79 "la_module.f90"
    !
    loffload = .false.
    !
    ! the following ifdef ensures no offload if not compiling from GPU 
!# 86 "la_module.f90"
    !
    ! ... always false when compiling without CUDA support
    !
    IF ( loffload ) THEN
!# 102 "la_module.f90"
    ELSE
      CALL laxlib_cdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE cdiaghg_cpu_
  !
!# 182 "la_module.f90"
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rdiaghg_cpu_( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by diaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! LAPACK version - uses both DSYGV and DSYGVX
    !!
    !
!# 199 "la_module.f90"
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculate
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,n)
    !! matrix to be diagonalized
    REAL(DP), INTENT(INOUT) :: s(ldh,n)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,m)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    LOGICAL, OPTIONAL ::  offload
    !! optionally solve the eigenvalue problem on the GPU   
    LOGICAL :: loffload
      !
!# 232 "la_module.f90"
    !
    loffload = .false.
    !
    ! the following ifdef ensures no offload if not compiling from GPU 
!# 239 "la_module.f90"
    !
    ! ... always false when compiling without CUDA support
    !
    IF ( loffload ) THEN
!# 255 "la_module.f90"
    ELSE
      CALL laxlib_rdiaghg(n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm)
    END IF
    !
    RETURN
    !
  END SUBROUTINE rdiaghg_cpu_
  !
!# 335 "la_module.f90"
  !
  !  === Parallel diagonalization interface subroutines
  !
  !
  !----------------------------------------------------------------------------
  SUBROUTINE prdiaghg_( n, h, s, ldh, e, v, idesc, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! real matrices version.
    !! On output both matrix are unchanged.
    !!
    !! Parallel version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,ldh)
    !! matrix to be diagonalized
    REAL(DP), INTENT(INOUT) :: s(ldh,ldh)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    LOGICAL, OPTIONAL ::  offload
    !! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload
!# 374 "la_module.f90"
    CALL laxlib_prdiaghg( n, h, s, ldh, e, v, idesc)
      
  END SUBROUTINE
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiaghg_( n, h, s, ldh, e, v, idesc, offload )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiaghg interface.
    !! Calculates eigenvalues and eigenvectors of the generalized problem.
    !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
    !! complex matrices version.
    !! On output both matrix are unchanged.
    !!
    !! Parallel version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh)
    !! matrix to be diagonalized
    COMPLEX(DP), INTENT(INOUT) :: s(ldh,ldh)
    !! overlap matrix
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    LOGICAL, OPTIONAL ::  offload
    !! place-holder, offloading on GPU not implemented yet
    LOGICAL :: loffload
!# 412 "la_module.f90"
    CALL laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc)
      
  END SUBROUTINE
  !
!# 507 "la_module.f90"
  !
  !----------------------------------------------------------------------------
  SUBROUTINE cdiagh_cpu_( n, m, h, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !----------------------------------------------------------------------------
    !
    !! Called by diagh interface.
    !! Calculates eigenvalues and eigenvectors of the standard problem.
    !! Solve Hv = ev, with H Hermitian matrix.
    !! complex matrices version.
    !! On output H matrix is unchanged.
    !!
    !! LAPACK version - uses both ZHEEVD and ZHEEVX
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,n)
    !! matrix to be diagonalized
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    !
    CALL laxlib_cdiagh( n, m, h, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE cdiagh_cpu_
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rdiagh_cpu_( n, m, h, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !----------------------------------------------------------------------------
    !
    !! Called by diagh interface.
    !! Calculates eigenvalues and eigenvectors of the standard problem.
    !! Solve Hv = ev, with H symmetric matrix.
    !! real matrices version.
    !! On output H matrix is unchanged.
    !!
    !! LAPACK version - uses both DSYEVD and DSYEVX
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized
    INTEGER, INTENT(IN) :: m
    !! number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,n)
    !! matrix to be diagonalized
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,m)
    !! eigenvectors (column-wise)
    INTEGER,  INTENT(IN)  :: me_bgrp
    !! index of the processor within a band group
    INTEGER,  INTENT(IN)  :: root_bgrp
    !! index of the root processor within a band group
    INTEGER,  INTENT(IN)  :: intra_bgrp_comm
    !! intra band group communicator
    !
    CALL laxlib_rdiagh( n, m, h, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE rdiagh_cpu_
  !
  !----------------------------------------------------------------------------
  SUBROUTINE pcdiagh_( n, h, ldh, e, v, idesc )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiagh interface.
    !! Calculates eigenvalues and eigenvectors of the standard problem.
    !! Solve Hv = ev, with H Hermitian matrix.
    !! complex matrices version.
    !! On output H matrix is unchanged.
    !!
    !! Parallel version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh)
    !! matrix to be diagonalized
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    !
    CALL laxlib_pcdiagh( n, h, ldh, e, v, idesc )
    !
  END SUBROUTINE pcdiagh_
  !
  !----------------------------------------------------------------------------
  SUBROUTINE prdiagh_( n, h, ldh, e, v, idesc )
    !----------------------------------------------------------------------------
    !
    !! Called by pdiagh interface.
    !! Calculates eigenvalues and eigenvectors of the standard problem.
    !! Solve Hv = ev, with H symmetric matrix.
    !! real matrices version.
    !! On output H matrix is unchanged.
    !!
    !! Parallel version with full data distribution
    !!
    !
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    INTEGER, INTENT(IN) :: n
    !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
    INTEGER, INTENT(IN) :: ldh
    !! leading dimension of h, as declared in the calling pgm unit
    REAL(DP), INTENT(INOUT) :: h(ldh,ldh)
    !! matrix to be diagonalized
    REAL(DP), INTENT(OUT) :: e(n)
    !! eigenvalues
    REAL(DP), INTENT(OUT) :: v(ldh,ldh)
    !! eigenvectors (column-wise)
    INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
    !! laxlib descriptor
    !
    CALL laxlib_prdiagh( n, h, ldh, e, v, idesc )
    !
  END SUBROUTINE prdiagh_
  !
END MODULE LAXlib
