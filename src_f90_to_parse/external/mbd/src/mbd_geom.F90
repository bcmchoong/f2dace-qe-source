!# 1 "mbd_geom.F90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!# 8 "mbd_geom.F90"
module mbd_geom
!! Representing a molecule or a crystal unit cell.
!# 11 "mbd_geom.F90"
use mbd_constants
use mbd_defaults
use mbd_formulas, only: alpha_dyn_qho, C6_from_alpha, omega_qho
use mbd_gradients, only: grad_t, grad_request_t
use mbd_lapack, only: eigvals, inverse
use mbd_utils, only: &
    shift_idx, atom_index_t, quad_pt_t, exception_t, tostr, clock_t, printer, &
    printer_i, logger_t
use mbd_vdw_param, only: ts_vdw_params
!# 27 "mbd_geom.F90"
implicit none
!# 29 "mbd_geom.F90"
private
public :: supercell_circum, get_freq_grid
!# 32 "mbd_geom.F90"
type, public :: param_t
    !! Calculation-wide paramters.
    real(dp) :: dipole_cutoff = 400d0 * ang  ! used only when Ewald is off
    real(dp) :: ewald_real_cutoff_scaling = 1d0
    real(dp) :: ewald_rec_cutoff_scaling = 1d0
    real(dp) :: k_grid_shift = K_GRID_SHIFT
    logical :: ewald_on = .true.
    logical :: zero_negative_eigvals = .false.
    logical :: rpa_rescale_eigs = .false.
    integer :: rpa_order_max = 10
    integer :: n_freq = N_FREQUENCY_GRID
end type
!# 45 "mbd_geom.F90"
type, public :: geom_t
    !! Represents a molecule or a crystal unit cell.
    !!
    !! The documented variables should be set before calling the initializer.
    real(dp), allocatable :: coords(:, :)
        !! (\(3\times N\), a.u.) Atomic coordinates.
    real(dp), allocatable :: lattice(:, :)
        !! (\(3\times 3\), a.u.) Lattice vectors in columns, unallocated if not
        !! periodic.
    integer, allocatable :: k_grid(:)
        !! Number of \(k\)-points along reciprocal axes.
    real(dp), allocatable :: custom_k_pts(:, :)
        !! Custom \(k\)-point grid.
    character(len=10) :: parallel_mode = 'auto'
        !! Type of parallelization:
        !!
        !! - `atoms`: distribute matrices over all MPI tasks using ScaLAPACK,
        !! solve eigenproblems sequentialy.
        !! - `k_points`: parallelize over k-points (each MPI task solves entire
        !! eigenproblems for its k-points)
    logical :: get_eigs = .false.
        !! Whether to keep MBD eigenvalues
    logical :: get_modes = .false.
        !! Whether to calculate MBD eigenvectors
    logical :: do_rpa = .false.
        !! Whether to calculate MBD energy by frequency integration
    logical :: get_rpa_orders = .false.
        !! Whether to calculate RPA orders
    type(logger_t) :: log
        !! Used for logging
!# 86 "mbd_geom.F90"
    ! The following components are set by the initializer and should be
    ! considered read-only
    type(clock_t) :: timer
    type(exception_t) :: exc
    type(quad_pt_t), allocatable :: freq(:)
    real(dp) :: gamm = 0d0
    real(dp) :: real_space_cutoff
    real(dp) :: rec_space_cutoff
    type(param_t) :: param
    type(atom_index_t) :: idx
!# 105 "mbd_geom.F90"
    contains
    procedure :: init => geom_init
    procedure :: destroy => geom_destroy
    procedure :: siz => geom_siz
    procedure :: has_exc => geom_has_exc
!# 113 "mbd_geom.F90"
    procedure :: clock => geom_clock
end type
!# 116 "mbd_geom.F90"
contains
!# 118 "mbd_geom.F90"
subroutine geom_init(this)
    class(geom_t), intent(inout) :: this
!# 121 "mbd_geom.F90"
    integer :: i_atom
    real(dp) :: volume, freq_grid_err
    logical :: is_parallel
    character(len=10) :: log_level_str
!# 130 "mbd_geom.F90"
    if (.not. associated(this%log%printer)) this%log%printer => printer
    call get_environment_variable('LIBMBD_LOG_LEVEL', log_level_str)
    if (log_level_str /= '') read (log_level_str, *) this%log%level
    associate (n => this%param%n_freq)
        allocate (this%freq(0:n))
        call get_freq_grid(n, this%freq(1:n)%val, this%freq(1:n)%weight)
    end associate
    this%freq(0)%val = 0d0
    this%freq(0)%weight = 0d0
    freq_grid_err = test_frequency_grid(this%freq)
    call this%log%info('Frequency grid relative error: '//tostr(freq_grid_err))
    call this%timer%init(100)
    if (allocated(this%lattice)) then
        volume = abs(dble(product(eigvals(this%lattice))))
        if (this%param%ewald_on) then
            this%gamm = 2.5d0 / volume**(1d0 / 3)
            this%real_space_cutoff = 6d0 / this%gamm * this%param%ewald_real_cutoff_scaling
            this%rec_space_cutoff = 10d0 * this%gamm * this%param%ewald_rec_cutoff_scaling
        else
            this%real_space_cutoff = this%param%dipole_cutoff
        end if
    end if
!# 195 "mbd_geom.F90"
    if (this%parallel_mode == 'auto') this%parallel_mode = 'none'
!# 199 "mbd_geom.F90"
    is_parallel = .false.
!# 201 "mbd_geom.F90"
    if (.not. is_parallel) then
        this%idx%i_atom = [(i_atom, i_atom=1, this%siz())]
        this%idx%j_atom = this%idx%i_atom
    end if
    this%idx%n_atoms = this%siz()
    call this%log%info('Will use parallel mode: '//this%parallel_mode)
!# 216 "mbd_geom.F90"
end subroutine
!# 218 "mbd_geom.F90"
subroutine geom_destroy(this)
    class(geom_t), intent(inout) :: this
!# 224 "mbd_geom.F90"
    deallocate (this%freq)
    deallocate (this%timer%timestamps)
    deallocate (this%timer%counts)
    deallocate (this%timer%levels)
end subroutine
!# 230 "mbd_geom.F90"
integer function geom_siz(this) result(siz)
    class(geom_t), intent(in) :: this
!# 233 "mbd_geom.F90"
    if (allocated(this%coords)) then
        siz = size(this%coords, 2)
    else
        siz = 0
    end if
end function
!# 240 "mbd_geom.F90"
logical function geom_has_exc(this) result(has_exc)
    class(geom_t), intent(in) :: this
!# 243 "mbd_geom.F90"
    has_exc = this%exc%code /= 0
end function
!# 270 "mbd_geom.F90"
function supercell_circum(lattice, radius) result(sc)
    real(dp), intent(in) :: lattice(3, 3)
    real(dp), intent(in) :: radius
    integer :: sc(3)
!# 275 "mbd_geom.F90"
    real(dp) :: ruc(3, 3), layer_sep(3)
    integer :: i
!# 278 "mbd_geom.F90"
    ruc = 2 * pi * inverse(transpose(lattice))
    do concurrent(i=1:3)
        layer_sep(i) = sum(lattice(:, i) * ruc(:, i) / sqrt(sum(ruc(:, i)**2)))
    end do
    sc = ceiling(radius / layer_sep + 0.5d0)
end function
!# 285 "mbd_geom.F90"
subroutine geom_clock(this, id)
    class(geom_t), intent(inout) :: this
    integer, intent(in) :: id
!# 289 "mbd_geom.F90"
    call this%timer%clock(id)
end subroutine
!# 292 "mbd_geom.F90"
subroutine get_freq_grid(n, x, w, L)
    integer, intent(in) :: n
    real(dp), intent(out) :: x(n), w(n)
    real(dp), intent(in), optional :: L
!# 297 "mbd_geom.F90"
    real(dp) :: L_
!# 299 "mbd_geom.F90"
    if (present(L)) then
        L_ = L
    else
        L_ = 0.6d0
    end if
    call gauss_legendre(n, x, w)
    w = 2 * L_ / (1 - x)**2 * w
    x = L_ * (1 + x) / (1 - x)
    w = w(n:1:-1)
    x = x(n:1:-1)
end subroutine
!# 311 "mbd_geom.F90"
subroutine gauss_legendre(n, r, w)
    integer, intent(in) :: n
    real(dp), intent(out) :: r(n), w(n)
!# 315 "mbd_geom.F90"
    integer, parameter :: q = selected_real_kind(15)
    integer, parameter :: n_iter = 1000
    real(q) :: x, f, df, dx
    integer :: k, iter, i
    real(q) :: Pk(0:n), Pk1(0:n - 1), Pk2(0:n - 2)
!# 321 "mbd_geom.F90"
    if (n == 1) then
        r(1) = 0d0
        w(1) = 2d0
        return
    end if
    Pk2(0) = 1._q  ! k = 0
    Pk1(0:1) = [0._q, 1._q]  ! k = 1
    do k = 2, n
        Pk(0:k) = ((2 * k - 1)* &
            [0._q, Pk1(0:k - 1)] - (k - 1)*[Pk2(0:k - 2), 0._q, 0._q]) / k
        if (k < n) then
            Pk2(0:k - 1) = Pk1(0:k - 1)
            Pk1(0:k) = Pk(0:k)
        end if
    end do
    ! now Pk contains k-th Legendre polynomial
    do i = 1, n
        x = cos(pi * (i - 0.25_q) / (n + 0.5_q))
        do iter = 1, n_iter
            df = 0._q
            f = Pk(n)
            do k = n - 1, 0, -1
                df = f + x * df
                f = Pk(k) + x * f
            end do
            dx = f / df
            x = x - dx
            if (abs(dx) < 10 * epsilon(dx)) exit
        end do
        r(i) = dble(x)
        w(i) = dble(2 / ((1 - x**2) * df**2))
    end do
end subroutine
!# 355 "mbd_geom.F90"
real(dp) function test_frequency_grid(freq) result(error)
    !! Calculate relative quadrature error in C6 of a carbon atom
    type(quad_pt_t), intent(in) :: freq(0:)
!# 359 "mbd_geom.F90"
    real(dp) :: alpha(1, 0:ubound(freq, 1)), C6(1), C6_ref(1), w(1), a0(1)
    type(grad_t), allocatable :: dalpha(:)
    type(grad_request_t) :: grad
!# 363 "mbd_geom.F90"
    a0(1) = ts_vdw_params(1, 6)
    C6_ref(1) = ts_vdw_params(2, 6)
    w = omega_qho(C6_ref, a0)
    alpha = alpha_dyn_qho(a0, w, freq, dalpha, grad)
    C6 = C6_from_alpha(alpha, freq)
    error = abs(C6(1) / C6_ref(1) - 1d0)
end function
!# 371 "mbd_geom.F90"
end module
