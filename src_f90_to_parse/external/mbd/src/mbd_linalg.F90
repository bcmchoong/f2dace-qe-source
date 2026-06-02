!# 1 "mbd_linalg.F90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_linalg
!# 6 "mbd_linalg.F90"
use mbd_constants
!# 8 "mbd_linalg.F90"
implicit none
!# 10 "mbd_linalg.F90"
private
public :: outer, eye, diag
!# 13 "mbd_linalg.F90"
interface diag
    module procedure get_diag_real
    module procedure get_diag_complex
    module procedure make_diag_real
end interface
!# 19 "mbd_linalg.F90"
contains
!# 21 "mbd_linalg.F90"
function outer(a, b) result(c)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: c(size(a), size(b))
!# 25 "mbd_linalg.F90"
    integer :: i, j
!# 27 "mbd_linalg.F90"
    do i = 1, size(a)
        do j = 1, size(b)
            c(i, j) = a(i) * b(j)
        end do
    end do
end function
!# 34 "mbd_linalg.F90"
function eye(n) result(A)
    integer, intent(in) :: n
    real(dp) :: A(n, n)
!# 38 "mbd_linalg.F90"
    integer :: i
!# 40 "mbd_linalg.F90"
    A(:, :) = 0.d0
    do concurrent(i=1:n)
        A(i, i) = 1.d0
    end do
end function
!# 46 "mbd_linalg.F90"
function get_diag_real(A) result(d)
    real(dp), intent(in) :: A(:, :)
    real(dp) :: d(size(A, 1))
!# 50 "mbd_linalg.F90"
    integer :: i
!# 52 "mbd_linalg.F90"
    do concurrent(i=1:size(A, 1))
        d(i) = A(i, i)
    end do
end function
!# 57 "mbd_linalg.F90"
function get_diag_complex(A) result(d)
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: d(size(A, 1))
!# 61 "mbd_linalg.F90"
    integer :: i
!# 63 "mbd_linalg.F90"
    do concurrent(i=1:size(A, 1))
        d(i) = A(i, i)
    end do
end function
!# 68 "mbd_linalg.F90"
function make_diag_real(d) result(A)
    real(dp), intent(in) :: d(:)
    real(dp) :: A(size(d), size(d))
!# 72 "mbd_linalg.F90"
    integer :: i
!# 74 "mbd_linalg.F90"
    A(:, :) = 0.d0
    do concurrent(i=1:size(d))
        A(i, i) = d(i)
    end do
end function
!# 80 "mbd_linalg.F90"
end module
