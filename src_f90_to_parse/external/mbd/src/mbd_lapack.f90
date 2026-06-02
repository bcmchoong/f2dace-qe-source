!# 1 "mbd_lapack.f90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_lapack
!# 6 "mbd_lapack.f90"
use mbd_constants
use mbd_utils, only: exception_t, tostr
!# 9 "mbd_lapack.f90"
implicit none
private
public :: mmul, inv, invh, inverse, eig, eigh, eigvals, eigvalsh, det, mode
!# 13 "mbd_lapack.f90"
interface mmul
    module procedure mmul_real
    module procedure mmul_complex
end interface
!# 18 "mbd_lapack.f90"
interface inv
    module procedure inv_real
end interface
!# 22 "mbd_lapack.f90"
interface invh
    module procedure invh_real
end interface
!# 26 "mbd_lapack.f90"
interface eig
    module procedure eig_real
end interface
!# 30 "mbd_lapack.f90"
interface eigh
    module procedure eigh_real
    module procedure eigh_complex
end interface
!# 35 "mbd_lapack.f90"
interface eigvals
    module procedure eigvals_real
end interface
!# 39 "mbd_lapack.f90"
interface eigvalsh
    module procedure eigvalsh_real
    module procedure eigvalsh_complex
end interface
!# 44 "mbd_lapack.f90"
interface
    ! The followinbg interfaces were taken straight from the LAPACK codebase,
    ! replacing COMPLEX*16 for COMPLEX(dp)
    SUBROUTINE ZHEEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
    import :: dp
    CHARACTER JOBZ, UPLO
    INTEGER INFO, LDA, LWORK, N
    DOUBLE PRECISION RWORK(*), W(*)
    COMPLEX(dp) A(LDA, *), WORK(*)
    END
    SUBROUTINE DGEEV(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
    CHARACTER JOBVL, JOBVR
    INTEGER INFO, LDA, LDVL, LDVR, LWORK, N
    DOUBLE PRECISION A(LDA, *), VL(LDVL, *), VR(LDVR, *), WI(*), WORK(*), WR(*)
    END
    SUBROUTINE DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
    CHARACTER JOBZ, UPLO
    INTEGER INFO, LDA, LWORK, N
    DOUBLE PRECISION A(LDA, *), W(*), WORK(*)
    END
    SUBROUTINE DGETRF(M, N, A, LDA, IPIV, INFO)
    INTEGER INFO, LDA, M, N
    INTEGER IPIV(*)
    DOUBLE PRECISION A(LDA, *)
    END
    SUBROUTINE DGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
    INTEGER INFO, LDA, LWORK, N
    INTEGER IPIV(*)
    DOUBLE PRECISION A(LDA, *), WORK(*)
    END
    SUBROUTINE DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
    INTEGER INFO, LDA, LDB, N, NRHS
    INTEGER IPIV(*)
    DOUBLE PRECISION A(LDA, *), B(LDB, *)
    END
    SUBROUTINE ZGETRF(M, N, A, LDA, IPIV, INFO)
    import :: dp
    INTEGER INFO, LDA, M, N
    INTEGER IPIV(*)
    COMPLEX(dp) A(LDA, *)
    END
    SUBROUTINE ZGETRI(N, A, LDA, IPIV, WORK, LWORK, INFO)
    import :: dp
    INTEGER INFO, LDA, LWORK, N
    INTEGER IPIV(*)
    COMPLEX(dp) A(LDA, *), WORK(*)
    END
    SUBROUTINE ZGEEV(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    import :: dp
    CHARACTER JOBVL, JOBVR
    INTEGER INFO, LDA, LDVL, LDVR, LWORK, N
    DOUBLE PRECISION RWORK(*)
    COMPLEX(dp) A(LDA, *), VL(LDVL, *), VR(LDVR, *), W(*), WORK(*)
    END
    SUBROUTINE DSYTRI(UPLO, N, A, LDA, IPIV, WORK, INFO)
    CHARACTER UPLO
    INTEGER INFO, LDA, N
    INTEGER IPIV(*)
    DOUBLE PRECISION A(LDA, *), WORK(*)
    END
    SUBROUTINE DSYTRF(UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO)
    CHARACTER UPLO
    INTEGER INFO, LDA, LWORK, N
    INTEGER IPIV(*)
    DOUBLE PRECISION A(LDA, *), WORK(*)
    END
    SUBROUTINE DGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    DOUBLE PRECISION ALPHA, BETA
    INTEGER K, LDA, LDB, LDC, M, N
    CHARACTER TRANSA, TRANSB
    DOUBLE PRECISION A(LDA, *), B(LDB, *), C(LDC, *)
    END
    SUBROUTINE ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    import :: dp
    COMPLEX(dp) ALPHA, BETA
    INTEGER K, LDA, LDB, LDC, M, N
    CHARACTER TRANSA, TRANSB
    COMPLEX(dp) A(LDA, *), B(LDB, *), C(LDC, *)
    END
end interface
!# 125 "mbd_lapack.f90"
contains
!# 127 "mbd_lapack.f90"
function inverse(A, exc)
    real(dp), intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp) :: inverse(size(A, 1), size(A, 2))
!# 132 "mbd_lapack.f90"
    call inv_real(inverse, exc, src=A)
end function
!# 135 "mbd_lapack.f90"
function eigvalsh_real(A, exc, destroy) result(eigvals)
    real(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvals(size(A, 1))
!# 141 "mbd_lapack.f90"
    real(dp), allocatable, target :: A_work(:, :)
    real(dp), pointer :: A_p(:, :)
!# 144 "mbd_lapack.f90"
    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eigh_real(A_p, eigvals, exc, vals_only=.true.)
end function
!# 157 "mbd_lapack.f90"
function eigvalsh_complex(A, exc, destroy) result(eigvals)
    complex(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigvals(size(A, 1))
!# 163 "mbd_lapack.f90"
    complex(dp), allocatable, target :: A_work(:, :)
    complex(dp), pointer :: A_p(:, :)
!# 166 "mbd_lapack.f90"
    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eigh_complex(A_p, eigvals, exc, vals_only=.true.)
end function
!# 179 "mbd_lapack.f90"
function eigvals_real(A, exc, destroy) result(eigvals)
    real(dp), target, intent(in) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigvals(size(A, 1))
!# 185 "mbd_lapack.f90"
    real(dp), allocatable, target :: A_work(:, :)
    real(dp), pointer :: A_p(:, :)
!# 188 "mbd_lapack.f90"
    nullify (A_p)
    if (present(destroy)) then
        if (destroy) then
            A_p => A
        end if
    end if
    if (.not. associated(A_p)) then
        allocate (A_work(size(A, 1), size(A, 1)), source=A)
        A_p => A_work
    end if
    call eig_real(A_p, eigvals, exc, vals_only=.true.)
end function
!# 201 "mbd_lapack.f90"
function mmul_real(A, B, transA, transB) result(C)
    real(dp), intent(in) :: A(:, :), B(:, :)
    character, intent(in), optional :: transA, transB
    real(dp) :: C(size(A, 1), size(B, 2))
!# 206 "mbd_lapack.f90"
    character :: transA_, transB_
    integer :: n
!# 209 "mbd_lapack.f90"
    transA_ = 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = size(A, 1)
    call DGEMM(transA_, transB_, n, n, n, 1d0, A, n, B, n, 0d0, C, n)
end function
!# 217 "mbd_lapack.f90"
function mmul_complex(A, B, transA, transB) result(C)
    complex(dp), intent(in) :: A(:, :), B(:, :)
    character, intent(in), optional :: transA, transB
    complex(dp) :: C(size(A, 1), size(B, 2))
!# 222 "mbd_lapack.f90"
    character :: transA_, transB_
    integer :: n
!# 225 "mbd_lapack.f90"
    transA_ = 'N'
    transB_ = 'N'
    if (present(transA)) transA_ = transA
    if (present(transB)) transB_ = transB
    n = size(A, 1)
    call ZGEMM(transA_, transB_, n, n, n, (1d0, 0d0), A, n, B, n, (0d0, 0d0), C, n)
end function
!# 233 "mbd_lapack.f90"
subroutine inv_real(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
!# 238 "mbd_lapack.f90"
    real(dp), allocatable :: work_arr(:)
    integer, allocatable :: i_pivot(:)
    integer :: n, n_work_arr, error_flag
    real(dp) :: n_work_arr_optim(1)
!# 243 "mbd_lapack.f90"
    n = size(A, 1)
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call DGETRF(n, n, A, n, i_pivot, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGETRF'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
    call DGETRI(n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim(1))
    allocate (work_arr(n_work_arr))
    call DGETRI(n, A, n, i_pivot, work_arr(1), n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGETRI'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
end subroutine
!# 270 "mbd_lapack.f90"
subroutine invh_real(A, exc, src)
    real(dp), intent(inout) :: A(:, :)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
!# 275 "mbd_lapack.f90"
    integer, allocatable :: i_pivot(:)
    real(dp), allocatable :: work_arr(:)
    integer :: n, n_work_arr, error_flag
    real(dp) :: n_work_arr_optim(1)
!# 280 "mbd_lapack.f90"
    n = size(A, 1)
    if (n == 0) return
    if (present(src)) A = src
    allocate (i_pivot(n))
    call DSYTRF('U', n, A, n, i_pivot, n_work_arr_optim, -1, error_flag)
    n_work_arr = nint(n_work_arr_optim(1))
    allocate (work_arr(n_work_arr))
    call DSYTRF('U', n, A, n, i_pivot, work_arr(1), n_work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYTRF'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
    deallocate (work_arr)
    allocate (work_arr(n))
    call DSYTRI('U', n, A, n, i_pivot, work_arr, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYTRI'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
    call fill_tril(A)
end subroutine
!# 310 "mbd_lapack.f90"
subroutine eigh_real(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only
!# 317 "mbd_lapack.f90"
    real(dp), allocatable :: work_arr(:)
    real(dp) :: n_work_arr(1)
    integer :: error_flag, n
!# 321 "mbd_lapack.f90"
    n = size(A, 1)
    if (present(src)) A = src
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, n_work_arr, -1, error_flag)
    allocate (work_arr(nint(n_work_arr(1))))
    call DSYEV(mode(vals_only), 'U', n, A, n, eigs, work_arr(1), size(work_arr), error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DSYEV'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
end subroutine
!# 336 "mbd_lapack.f90"
subroutine eig_real(A, eigs, exc, src, vals_only)
    real(dp), intent(inout) :: A(:, :)
    complex(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    real(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only
!# 343 "mbd_lapack.f90"
    real(dp) :: n_work_arr(1), dummy(1)
    integer :: error_flag, n
    real(dp), allocatable :: eigs_r(:), eigs_i(:), vectors(:, :), work_arr(:)
!# 347 "mbd_lapack.f90"
    n = size(A, 1)
    if (present(src)) A = src
    allocate (eigs_r(n), eigs_i(n))
    if (mode(vals_only) == 'V') then
        allocate (vectors(n, n))
    else
        allocate (vectors(1, 1))
    end if
    call DGEEV( &
        'N', mode(vals_only), n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, n_work_arr, -1, error_flag &
    )
    allocate (work_arr(nint(n_work_arr(1))))
    call DGEEV( &
        'N', mode(vals_only), n, A, n, eigs_r, eigs_i, dummy, 1, &
        vectors, n, work_arr(1), size(work_arr), error_flag &
    )
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'DGEEV'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
    eigs = cmplx(eigs_r, eigs_i, dp)
    if (mode(vals_only) == 'V') A = vectors
end subroutine
!# 376 "mbd_lapack.f90"
subroutine eigh_complex(A, eigs, exc, src, vals_only)
    complex(dp), intent(inout) :: A(:, :)
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    complex(dp), intent(in), optional :: src(:, :)
    logical, intent(in), optional :: vals_only
!# 383 "mbd_lapack.f90"
    complex(dp), allocatable :: work(:)
    complex(dp) :: lwork_cmplx(1)
    real(dp), allocatable :: rwork(:)
    integer :: n, lwork, error_flag
!# 388 "mbd_lapack.f90"
    n = size(A, 1)
    if (present(src)) A = src
    allocate (rwork(max(1, 3 * n - 2)))
    call ZHEEV(mode(vals_only), 'U', n, A, n, eigs, lwork_cmplx, -1, rwork, error_flag)
    lwork = nint(dble(lwork_cmplx(1)))
    allocate (work(lwork))
    call ZHEEV(mode(vals_only), 'U', n, A, n, eigs, work(1), lwork, rwork, error_flag)
    if (error_flag /= 0) then
        if (present(exc)) then
            exc%code = MBD_EXC_LINALG
            exc%origin = 'ZHEEV'
            exc%msg = 'Failed with code '//trim(tostr(error_flag))
        end if
        return
    end if
end subroutine
!# 405 "mbd_lapack.f90"
real(dp) function det(A) result(D)
    real(dp), intent(in) :: A(:, :)
!# 408 "mbd_lapack.f90"
    integer :: n, i, info
    real(dp), allocatable :: LU(:, :)
    integer, allocatable :: ipiv(:)
!# 412 "mbd_lapack.f90"
    n = size(A, 1)
    allocate (ipiv(n))
    LU = A
    call DGETRF(n, n, LU, n, ipiv, info)
    D = product([(LU(i, i), i=1, n)])
end function
!# 419 "mbd_lapack.f90"
subroutine fill_tril(A)
    real(dp), intent(inout) :: A(:, :)
!# 422 "mbd_lapack.f90"
    integer :: i, j
!# 424 "mbd_lapack.f90"
    do i = 1, size(A, 1)
        do j = i + 1, size(A, 1)
            A(j, i) = A(i, j)
        end do
    end do
end subroutine
!# 431 "mbd_lapack.f90"
character(len=1) function mode(vals_only)
    logical, intent(in), optional :: vals_only
!# 434 "mbd_lapack.f90"
    mode = 'V'
    if (present(vals_only)) then
        if (vals_only) mode = 'N'
    end if
end function
!# 440 "mbd_lapack.f90"
end module
