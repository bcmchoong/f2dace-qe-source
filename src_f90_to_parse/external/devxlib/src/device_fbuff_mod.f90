!# 1 "device_fbuff_mod.f90"
!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contains a CUF KERNEL to perform the selected operation
! cu_memsync are wrappers for cuda_memcpy functions
!
!# 15 "device_fbuff_mod.f90"
!
module device_fbuff_m
  use tb_dev, only : tb_dev_t
  use tb_pin, only : tb_pin_t
!# 20 "device_fbuff_mod.f90"
  implicit none
!# 22 "device_fbuff_mod.f90"
  TYPE(tb_dev_t) :: dev_buf ! A global variable hosting a global device buffer
  TYPE(tb_pin_t) :: pin_buf ! A global variable hosting a global pinned memory buffer
!# 25 "device_fbuff_mod.f90"
!# 1 "./device_fbuff_interf.f90"
!# 26 "device_fbuff_mod.f90"
!# 27 "device_fbuff_mod.f90"
end module device_fbuff_m
