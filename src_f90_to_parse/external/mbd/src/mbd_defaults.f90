!# 1 "mbd_defaults.f90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!# 5 "mbd_defaults.f90"
module mbd_defaults
!! Defaults used at multiple places.
!# 8 "mbd_defaults.f90"
use mbd_constants
!# 10 "mbd_defaults.f90"
implicit none
!# 12 "mbd_defaults.f90"
real(dp), parameter :: N_FREQUENCY_GRID = 15
real(dp), parameter :: K_GRID_SHIFT = 0.5d0
real(dp), parameter :: TS_DAMPING_D = 20d0
real(dp), parameter :: MBD_DAMPING_A = 6d0
real(dp), parameter :: MAX_ATOMS_PER_BLOCK = 6
!# 18 "mbd_defaults.f90"
end module
