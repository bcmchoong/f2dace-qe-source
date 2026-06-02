!# 1 "mbd_version.f90"
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
module mbd_version
!# 6 "mbd_version.f90"
implicit none
!# 8 "mbd_version.f90"
integer, parameter, public :: MBD_VERSION_MAJOR = 0
integer, parameter, public :: MBD_VERSION_MINOR = 13
integer, parameter, public :: MBD_VERSION_PATCH = 0
character(len=30), parameter, public :: MBD_VERSION_SUFFIX = '0'
!# 13 "mbd_version.f90"
end module
