!# 1 "mbd_matrix.F90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!# 5 "mbd_matrix.F90"
module mbd_matrix
!# 7 "mbd_matrix.F90"
use mbd_constants
use mbd_lapack, only: mmul, invh, invh, eigh, eigvals, eigvalsh
use mbd_utils, only: findval, exception_t, atom_index_t, is_true, clock_t
!# 18 "mbd_matrix.F90"
implicit none
!# 20 "mbd_matrix.F90"
private
public :: contract_cross_33
!# 23 "mbd_matrix.F90"
type, public :: matrix_re_t
    real(dp), allocatable :: val(:, :)
    type(atom_index_t) :: idx
!# 29 "mbd_matrix.F90"
    contains
    procedure :: siz => matrix_re_siz
    procedure :: init => matrix_re_init
    procedure :: add_diag => matrix_re_add_diag
    procedure :: add_diag_scalar => matrix_re_add_diag_scalar
    procedure :: mult_cross => matrix_re_mult_cross
    procedure :: mult_rows => matrix_re_mult_rows
    procedure :: mult_cols_3n => matrix_re_mult_cols_3n
    procedure :: mult_col => matrix_re_mult_col
    procedure :: mmul => matrix_re_mmul
    procedure :: invh => matrix_re_invh
    procedure :: eigh => matrix_re_eigh
    procedure :: eigvals => matrix_re_eigvals
    procedure :: eigvalsh => matrix_re_eigvalsh
    procedure :: sum_all => matrix_re_sum_all
    procedure :: contract_n_transp => matrix_re_contract_n_transp
    procedure :: contract_n33diag_cols => matrix_re_contract_n33diag_cols
    procedure :: contract_n33_rows => matrix_re_contract_n33_rows
    procedure :: copy_from => matrix_re_copy_from
    procedure :: move_from => matrix_re_move_from
    procedure :: init_from => matrix_re_init_from
    procedure :: alloc_from => matrix_re_alloc_from
end type
!# 53 "mbd_matrix.F90"
type, public :: matrix_cplx_t
    complex(dp), allocatable :: val(:, :)
    type(atom_index_t) :: idx
!# 59 "mbd_matrix.F90"
    contains
    procedure :: siz => matrix_cplx_siz
    procedure :: init => matrix_cplx_init
    procedure :: add_diag => matrix_cplx_add_diag
    procedure :: add_diag_scalar => matrix_cplx_add_diag_scalar
    procedure :: mult_cross => matrix_cplx_mult_cross
    procedure :: mult_rows => matrix_cplx_mult_rows
    procedure :: mult_cols_3n => matrix_cplx_mult_cols_3n
    procedure :: mult_col => matrix_cplx_mult_col
    procedure :: mmul => matrix_cplx_mmul
    procedure :: eigh => matrix_cplx_eigh
    procedure :: eigvalsh => matrix_cplx_eigvalsh
    procedure :: sum_all => matrix_cplx_sum_all
    procedure :: contract_n_transp => matrix_cplx_contract_n_transp
    procedure :: contract_n33diag_cols => matrix_cplx_contract_n33diag_cols
    procedure :: contract_n33_rows => matrix_cplx_contract_n33_rows
    procedure :: copy_from => matrix_cplx_copy_from
    procedure :: move_from => matrix_cplx_move_from
    procedure :: init_from => matrix_cplx_init_from
    procedure :: alloc_from => matrix_cplx_alloc_from
end type
!# 81 "mbd_matrix.F90"
interface contract_cross_33
    module procedure contract_cross_33_real
    module procedure contract_cross_33_complex
end interface
!# 86 "mbd_matrix.F90"
contains
!# 91 "mbd_matrix.F90"
integer function matrix_re_siz(this, ndim) result(siz)
    class(matrix_re_t), intent(in) :: this
!# 97 "mbd_matrix.F90"
    integer, intent(in) :: ndim
!# 99 "mbd_matrix.F90"
    siz = size(this%val, ndim)
end function
!# 106 "mbd_matrix.F90"
subroutine matrix_re_init(this, idx)
!# 108 "mbd_matrix.F90"
    class(matrix_re_t), intent(out) :: this
!# 117 "mbd_matrix.F90"
    type(atom_index_t), intent(in) :: idx
!# 122 "mbd_matrix.F90"
    this%idx = idx
!# 126 "mbd_matrix.F90"
end subroutine
!# 129 "mbd_matrix.F90"
subroutine matrix_re_init_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other
!# 138 "mbd_matrix.F90"
    this%idx = other%idx
!# 142 "mbd_matrix.F90"
end subroutine
!# 145 "mbd_matrix.F90"
subroutine matrix_re_copy_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other
!# 154 "mbd_matrix.F90"
    call this%init_from(other)
    this%val = other%val
end subroutine
!# 159 "mbd_matrix.F90"
subroutine matrix_re_move_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(inout) :: other
!# 168 "mbd_matrix.F90"
    call this%init_from(other)
    call move_alloc(other%val, this%val)
end subroutine
!# 173 "mbd_matrix.F90"
subroutine matrix_re_alloc_from(this, other)
    class(matrix_re_t), intent(out) :: this
    type(matrix_re_t), intent(in) :: other
!# 182 "mbd_matrix.F90"
    integer :: n1, n2
!# 184 "mbd_matrix.F90"
    call this%init_from(other)
    n1 = other%siz(1)
    n2 = other%siz(2)
    allocate (this%val(n1, n2))
end subroutine
!# 191 "mbd_matrix.F90"
subroutine matrix_re_add_diag_scalar(this, d)
    class(matrix_re_t), intent(inout) :: this
!# 197 "mbd_matrix.F90"
    real(dp), intent(in) :: d
!# 199 "mbd_matrix.F90"
    integer :: i
!# 201 "mbd_matrix.F90"
    call this%add_diag([(d, i=1, this%idx%n_atoms)])
end subroutine
!# 205 "mbd_matrix.F90"
subroutine matrix_re_add_diag(this, d)
    class(matrix_re_t), intent(inout) :: this
!# 211 "mbd_matrix.F90"
    real(dp), intent(in) :: d(:)
!# 213 "mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom, i
!# 215 "mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_diag => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                if (i_atom /= j_atom) cycle
                do i = 1, 3
                    this_diag(i, i) = this_diag(i, i) + d(i_atom)
                end do
            end associate
        end do
    end do
end subroutine
!# 232 "mbd_matrix.F90"
subroutine matrix_re_mult_cross(this, b, c)
    class(matrix_re_t), intent(inout) :: this
!# 238 "mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)
!# 241 "mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 243 "mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                if (present(c)) then
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (b(i_atom) * c(j_atom) + c(i_atom) * b(j_atom))
                else
                    this_sub(:3, :3) = this_sub(:3, :3) * b(i_atom) * b(j_atom)
                end if
            end associate
        end do
    end do
end subroutine
!# 262 "mbd_matrix.F90"
subroutine matrix_re_mult_rows(this, b)
    class(matrix_re_t), intent(inout) :: this
!# 268 "mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
!# 270 "mbd_matrix.F90"
    integer :: my_i_atom
!# 272 "mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        associate ( &
                i_atom => this%idx%i_atom(my_i_atom), &
                this_sub => this%val(3 * (my_i_atom - 1) + 1:, :) &
        )
            this_sub(:3, :) = this_sub(:3, :) * b(i_atom)
        end associate
    end do
end subroutine
!# 283 "mbd_matrix.F90"
subroutine matrix_re_mult_cols_3n(this, b)
    class(matrix_re_t), intent(inout) :: this
!# 289 "mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
!# 291 "mbd_matrix.F90"
    integer :: my_j_atom, i
!# 293 "mbd_matrix.F90"
    do my_j_atom = 1, size(this%idx%j_atom)
        associate ( &
                b_sub => b(3 * (this%idx%j_atom(my_j_atom) - 1) + 1:), &
                this_sub => this%val(:, 3 * (my_j_atom - 1) + 1:) &
        )
            ! TODO should be do-concurrent, but this crashes IBM XL 16.1.1,
            ! see issue #16
            do i = 1, 3
                this_sub(:, i) = this_sub(:, i) * b_sub(i)
            end do
        end associate
    end do
end subroutine
!# 308 "mbd_matrix.F90"
subroutine matrix_re_mult_col(this, idx, a)
    class(matrix_re_t), intent(inout) :: this
!# 314 "mbd_matrix.F90"
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)
!# 317 "mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 319 "mbd_matrix.F90"
    do my_j_atom = 1, size(this%idx%j_atom)
        if (this%idx%j_atom(my_j_atom) /= idx) cycle
        do my_i_atom = 1, size(this%idx%i_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                this_sub(:3, :3) = this_sub(:3, :3) * a(i_atom)
            end associate
        end do
    end do
end subroutine
!# 333 "mbd_matrix.F90"
subroutine matrix_re_eigh(A, eigs, exc, src, vals_only, clock)
    class(matrix_re_t), intent(inout) :: A
    type(matrix_re_t), intent(in), optional :: src
!# 341 "mbd_matrix.F90"
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only
    type(clock_t), intent(inout), optional :: clock
!# 358 "mbd_matrix.F90"
    call eigh(A%val, eigs, exc, src%val, vals_only)
end subroutine
!# 362 "mbd_matrix.F90"
function matrix_re_eigvalsh(A, exc, destroy, clock) result(eigs)
    class(matrix_re_t), intent(inout) :: A
!# 368 "mbd_matrix.F90"
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3 * A%idx%n_atoms)
    type(clock_t), intent(inout), optional :: clock
!# 383 "mbd_matrix.F90"
    eigs = eigvalsh(A%val, exc, destroy)
end function
!# 387 "mbd_matrix.F90"
real(dp) function matrix_re_sum_all(this) result(res)
    class(matrix_re_t), intent(in) :: this
!# 394 "mbd_matrix.F90"
    res = sum(this%val)
!# 398 "mbd_matrix.F90"
end function
!# 401 "mbd_matrix.F90"
subroutine matrix_re_contract_n_transp(this, dir, res)
    class(matrix_re_t), intent(in) :: this
    real(dp), intent(out), target :: res(:, :)
!# 409 "mbd_matrix.F90"
    character(len=*), intent(in) :: dir
!# 411 "mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 413 "mbd_matrix.F90"
    real(dp), pointer :: res_sub(:, :)
!# 418 "mbd_matrix.F90"
    res(:, :) = 0d0
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3 * (this%idx%i_atom(my_i_atom) - 1) + 1:)
            case ('C')
                res_sub => res(3 * (this%idx%j_atom(my_j_atom) - 1) + 1:, :)
            end select
            associate ( &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
!# 437 "mbd_matrix.F90"
end subroutine
!# 440 "mbd_matrix.F90"
function contract_cross_33_real(k_atom, A, A_prime, B, B_prime) result(res)
    type(matrix_re_t), intent(in) :: A, B
    real(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    real(dp) :: res(A%idx%n_atoms)
!# 450 "mbd_matrix.F90"
    integer, intent(in) :: k_atom
!# 452 "mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom, i_atom, j_atom
!# 454 "mbd_matrix.F90"
    res(:) = 0d0
    my_i_atom = findval(A%idx%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%idx%j_atom)
            j_atom = A%idx%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:), &
                    A_prime_sub => A_prime(:, 3 * (j_atom - 1) + 1:) &
            )
                res(j_atom) = -1d0 / 3 * sum(A_sub(:3, :3) * A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%idx%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%idx%i_atom)
            i_atom = A%idx%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:), &
                    B_prime_sub => B_prime(3 * (i_atom - 1) + 1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0 / 3) * sum(B_prime_sub(:3, :) * B_sub(:3, :3))
            end associate
        end do
    end if
!# 483 "mbd_matrix.F90"
end function
!# 486 "mbd_matrix.F90"
function matrix_re_contract_n33diag_cols(A) result(res)
    class(matrix_re_t), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)
!# 495 "mbd_matrix.F90"
    integer :: i_xyz, my_j_atom, j_atom
!# 497 "mbd_matrix.F90"
    res(:) = 0d0
    do my_j_atom = 1, size(A%idx%j_atom)
        j_atom = A%idx%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%val(i_xyz::3, 3 * (my_j_atom - 1) + i_xyz))
        end do
    end do
    res = res / 3
!# 509 "mbd_matrix.F90"
end function
!# 512 "mbd_matrix.F90"
function matrix_re_contract_n33_rows(A) result(res)
    class(matrix_re_t), intent(in) :: A
    real(dp) :: res(A%idx%n_atoms)
!# 521 "mbd_matrix.F90"
    integer :: my_i_atom, i_atom
!# 523 "mbd_matrix.F90"
    res(:) = 0d0
    do my_i_atom = 1, size(A%idx%i_atom)
        i_atom = A%idx%i_atom(my_i_atom)
        associate (A_sub => A%val(3 * (my_i_atom - 1) + 1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
!# 533 "mbd_matrix.F90"
end function
!# 536 "mbd_matrix.F90"
type(matrix_re_t) function matrix_re_mmul( &
        A, B, transA, transB) result(C)
    class(matrix_re_t), intent(in) :: A
    type(matrix_re_t), intent(in) :: B
!# 546 "mbd_matrix.F90"
    character, intent(in), optional :: transA, transB
!# 548 "mbd_matrix.F90"
    C%idx = A%idx
!# 557 "mbd_matrix.F90"
    C%val = mmul(A%val, B%val, transA, transB)
!# 559 "mbd_matrix.F90"
end function
!# 563 "mbd_matrix.F90"
!# 1 "./mbd_matrix.F90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!# 94 "./mbd_matrix.F90"
integer function matrix_cplx_siz(this, ndim) result(siz)
    class(matrix_cplx_t), intent(in) :: this
!# 97 "./mbd_matrix.F90"
    integer, intent(in) :: ndim
!# 99 "./mbd_matrix.F90"
    siz = size(this%val, ndim)
end function
!# 113 "./mbd_matrix.F90"
subroutine matrix_cplx_init(this, idx)
!# 115 "./mbd_matrix.F90"
    class(matrix_cplx_t), intent(out) :: this
!# 117 "./mbd_matrix.F90"
    type(atom_index_t), intent(in) :: idx
!# 122 "./mbd_matrix.F90"
    this%idx = idx
!# 126 "./mbd_matrix.F90"
end subroutine
!# 133 "./mbd_matrix.F90"
subroutine matrix_cplx_init_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other
!# 138 "./mbd_matrix.F90"
    this%idx = other%idx
!# 142 "./mbd_matrix.F90"
end subroutine
!# 149 "./mbd_matrix.F90"
subroutine matrix_cplx_copy_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other
!# 154 "./mbd_matrix.F90"
    call this%init_from(other)
    this%val = other%val
end subroutine
!# 163 "./mbd_matrix.F90"
subroutine matrix_cplx_move_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(inout) :: other
!# 168 "./mbd_matrix.F90"
    call this%init_from(other)
    call move_alloc(other%val, this%val)
end subroutine
!# 177 "./mbd_matrix.F90"
subroutine matrix_cplx_alloc_from(this, other)
    class(matrix_cplx_t), intent(out) :: this
    type(matrix_cplx_t), intent(in) :: other
!# 182 "./mbd_matrix.F90"
    integer :: n1, n2
!# 184 "./mbd_matrix.F90"
    call this%init_from(other)
    n1 = other%siz(1)
    n2 = other%siz(2)
    allocate (this%val(n1, n2))
end subroutine
!# 194 "./mbd_matrix.F90"
subroutine matrix_cplx_add_diag_scalar(this, d)
    class(matrix_cplx_t), intent(inout) :: this
!# 197 "./mbd_matrix.F90"
    real(dp), intent(in) :: d
!# 199 "./mbd_matrix.F90"
    integer :: i
!# 201 "./mbd_matrix.F90"
    call this%add_diag([(d, i=1, this%idx%n_atoms)])
end subroutine
!# 208 "./mbd_matrix.F90"
subroutine matrix_cplx_add_diag(this, d)
    class(matrix_cplx_t), intent(inout) :: this
!# 211 "./mbd_matrix.F90"
    real(dp), intent(in) :: d(:)
!# 213 "./mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom, i
!# 215 "./mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_diag => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                if (i_atom /= j_atom) cycle
                do i = 1, 3
                    this_diag(i, i) = this_diag(i, i) + d(i_atom)
                end do
            end associate
        end do
    end do
end subroutine
!# 235 "./mbd_matrix.F90"
subroutine matrix_cplx_mult_cross(this, b, c)
    class(matrix_cplx_t), intent(inout) :: this
!# 238 "./mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: c(:)
!# 241 "./mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 243 "./mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    j_atom => this%idx%j_atom(my_j_atom), &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                if (present(c)) then
                    this_sub(:3, :3) = this_sub(:3, :3) * &
                        (b(i_atom) * c(j_atom) + c(i_atom) * b(j_atom))
                else
                    this_sub(:3, :3) = this_sub(:3, :3) * b(i_atom) * b(j_atom)
                end if
            end associate
        end do
    end do
end subroutine
!# 265 "./mbd_matrix.F90"
subroutine matrix_cplx_mult_rows(this, b)
    class(matrix_cplx_t), intent(inout) :: this
!# 268 "./mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
!# 270 "./mbd_matrix.F90"
    integer :: my_i_atom
!# 272 "./mbd_matrix.F90"
    do my_i_atom = 1, size(this%idx%i_atom)
        associate ( &
                i_atom => this%idx%i_atom(my_i_atom), &
                this_sub => this%val(3 * (my_i_atom - 1) + 1:, :) &
        )
            this_sub(:3, :) = this_sub(:3, :) * b(i_atom)
        end associate
    end do
end subroutine
!# 286 "./mbd_matrix.F90"
subroutine matrix_cplx_mult_cols_3n(this, b)
    class(matrix_cplx_t), intent(inout) :: this
!# 289 "./mbd_matrix.F90"
    real(dp), intent(in) :: b(:)
!# 291 "./mbd_matrix.F90"
    integer :: my_j_atom, i
!# 293 "./mbd_matrix.F90"
    do my_j_atom = 1, size(this%idx%j_atom)
        associate ( &
                b_sub => b(3 * (this%idx%j_atom(my_j_atom) - 1) + 1:), &
                this_sub => this%val(:, 3 * (my_j_atom - 1) + 1:) &
        )
            ! TODO should be do-concurrent, but this crashes IBM XL 16.1.1,
            ! see issue #16
            do i = 1, 3
                this_sub(:, i) = this_sub(:, i) * b_sub(i)
            end do
        end associate
    end do
end subroutine
!# 311 "./mbd_matrix.F90"
subroutine matrix_cplx_mult_col(this, idx, a)
    class(matrix_cplx_t), intent(inout) :: this
!# 314 "./mbd_matrix.F90"
    integer, intent(in) :: idx
    real(dp), intent(in) :: a(:)
!# 317 "./mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 319 "./mbd_matrix.F90"
    do my_j_atom = 1, size(this%idx%j_atom)
        if (this%idx%j_atom(my_j_atom) /= idx) cycle
        do my_i_atom = 1, size(this%idx%i_atom)
            associate ( &
                    i_atom => this%idx%i_atom(my_i_atom), &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                this_sub(:3, :3) = this_sub(:3, :3) * a(i_atom)
            end associate
        end do
    end do
end subroutine
!# 337 "./mbd_matrix.F90"
subroutine matrix_cplx_eigh(A, eigs, exc, src, vals_only, clock)
    class(matrix_cplx_t), intent(inout) :: A
    type(matrix_cplx_t), intent(in), optional :: src
!# 341 "./mbd_matrix.F90"
    real(dp), intent(out) :: eigs(:)
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: vals_only
    type(clock_t), intent(inout), optional :: clock
!# 358 "./mbd_matrix.F90"
    call eigh(A%val, eigs, exc, src%val, vals_only)
end subroutine
!# 365 "./mbd_matrix.F90"
function matrix_cplx_eigvalsh(A, exc, destroy, clock) result(eigs)
    class(matrix_cplx_t), intent(inout) :: A
!# 368 "./mbd_matrix.F90"
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    real(dp) :: eigs(3 * A%idx%n_atoms)
    type(clock_t), intent(inout), optional :: clock
!# 383 "./mbd_matrix.F90"
    eigs = eigvalsh(A%val, exc, destroy)
end function
!# 390 "./mbd_matrix.F90"
complex(dp) function matrix_cplx_sum_all(this) result(res)
    class(matrix_cplx_t), intent(in) :: this
!# 394 "./mbd_matrix.F90"
    res = sum(this%val)
!# 398 "./mbd_matrix.F90"
end function
!# 405 "./mbd_matrix.F90"
subroutine matrix_cplx_contract_n_transp(this, dir, res)
    class(matrix_cplx_t), intent(in) :: this
    complex(dp), intent(out), target :: res(:, :)
!# 409 "./mbd_matrix.F90"
    character(len=*), intent(in) :: dir
!# 411 "./mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom
!# 415 "./mbd_matrix.F90"
    complex(dp), pointer :: res_sub(:, :)
!# 418 "./mbd_matrix.F90"
    res(:, :) = 0d0
    do my_i_atom = 1, size(this%idx%i_atom)
        do my_j_atom = 1, size(this%idx%j_atom)
            select case (dir(1:1))
            case ('R')
                res_sub => res(:, 3 * (this%idx%i_atom(my_i_atom) - 1) + 1:)
            case ('C')
                res_sub => res(3 * (this%idx%j_atom(my_j_atom) - 1) + 1:, :)
            end select
            associate ( &
                    this_sub => this%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:) &
            )
                res_sub(:3, :3) = res_sub(:3, :3) + transpose(this_sub(:3, :3))
            end associate
        end do
    end do
!# 437 "./mbd_matrix.F90"
end subroutine
!# 445 "./mbd_matrix.F90"
function contract_cross_33_complex(k_atom, A, A_prime, B, B_prime) result(res)
    type(matrix_cplx_t), intent(in) :: A, B
    complex(dp), intent(in) :: A_prime(:, :), B_prime(:, :)
    complex(dp) :: res(A%idx%n_atoms)
!# 450 "./mbd_matrix.F90"
    integer, intent(in) :: k_atom
!# 452 "./mbd_matrix.F90"
    integer :: my_i_atom, my_j_atom, i_atom, j_atom
!# 454 "./mbd_matrix.F90"
    res(:) = 0d0
    my_i_atom = findval(A%idx%i_atom, k_atom)
    if (my_i_atom > 0) then
        do my_j_atom = 1, size(A%idx%j_atom)
            j_atom = A%idx%j_atom(my_j_atom)
            associate ( &
                    A_sub => A%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:), &
                    A_prime_sub => A_prime(:, 3 * (j_atom - 1) + 1:) &
            )
                res(j_atom) = -1d0 / 3 * sum(A_sub(:3, :3) * A_prime_sub(:, :3))
            end associate
        end do
    end if
    my_j_atom = findval(A%idx%j_atom, k_atom)
    if (my_j_atom > 0) then
        do my_i_atom = 1, size(A%idx%i_atom)
            i_atom = A%idx%i_atom(my_i_atom)
            associate ( &
                    B_sub => B%val(3 * (my_i_atom - 1) + 1:, 3 * (my_j_atom - 1) + 1:), &
                    B_prime_sub => B_prime(3 * (i_atom - 1) + 1:, :) &
            )
                res(i_atom) = res(i_atom) + &
                    (-1d0 / 3) * sum(B_prime_sub(:3, :) * B_sub(:3, :3))
            end associate
        end do
    end if
!# 483 "./mbd_matrix.F90"
end function
!# 490 "./mbd_matrix.F90"
function matrix_cplx_contract_n33diag_cols(A) result(res)
    class(matrix_cplx_t), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)
!# 495 "./mbd_matrix.F90"
    integer :: i_xyz, my_j_atom, j_atom
!# 497 "./mbd_matrix.F90"
    res(:) = 0d0
    do my_j_atom = 1, size(A%idx%j_atom)
        j_atom = A%idx%j_atom(my_j_atom)
        do i_xyz = 1, 3
            res(j_atom) = res(j_atom) + &
                sum(A%val(i_xyz::3, 3 * (my_j_atom - 1) + i_xyz))
        end do
    end do
    res = res / 3
!# 509 "./mbd_matrix.F90"
end function
!# 516 "./mbd_matrix.F90"
function matrix_cplx_contract_n33_rows(A) result(res)
    class(matrix_cplx_t), intent(in) :: A
    complex(dp) :: res(A%idx%n_atoms)
!# 521 "./mbd_matrix.F90"
    integer :: my_i_atom, i_atom
!# 523 "./mbd_matrix.F90"
    res(:) = 0d0
    do my_i_atom = 1, size(A%idx%i_atom)
        i_atom = A%idx%i_atom(my_i_atom)
        associate (A_sub => A%val(3 * (my_i_atom - 1) + 1:, :))
            res(i_atom) = res(i_atom) + sum(A_sub(:3, :))
        end associate
    end do
!# 533 "./mbd_matrix.F90"
end function
!# 541 "./mbd_matrix.F90"
type(matrix_cplx_t) function matrix_cplx_mmul( &
        A, B, transA, transB) result(C)
    class(matrix_cplx_t), intent(in) :: A
    type(matrix_cplx_t), intent(in) :: B
!# 546 "./mbd_matrix.F90"
    character, intent(in), optional :: transA, transB
!# 548 "./mbd_matrix.F90"
    C%idx = A%idx
!# 557 "./mbd_matrix.F90"
    C%val = mmul(A%val, B%val, transA, transB)
!# 559 "./mbd_matrix.F90"
end function
!# 564 "mbd_matrix.F90"
!# 565 "mbd_matrix.F90"
subroutine matrix_re_invh(A, exc, src, clock)
    class(matrix_re_t), intent(inout) :: A
    type(matrix_re_t), intent(in), optional :: src
    type(exception_t), intent(out), optional :: exc
    type(clock_t), intent(inout), optional :: clock
!# 586 "mbd_matrix.F90"
    if (present(src)) then
        call invh(A%val, exc, src%val)
    else
        call invh(A%val, exc)
    end if
!# 592 "mbd_matrix.F90"
end subroutine
!# 594 "mbd_matrix.F90"
function matrix_re_eigvals(A, exc, destroy) result(eigs)
    class(matrix_re_t), target, intent(in) :: A
    type(exception_t), intent(out), optional :: exc
    logical, intent(in), optional :: destroy
    complex(dp) :: eigs(3 * A%idx%n_atoms)
!# 608 "mbd_matrix.F90"
    eigs = eigvals(A%val, exc, destroy)
!# 610 "mbd_matrix.F90"
end function
!# 612 "mbd_matrix.F90"
end module
