!# 1 "device_memcpy_mod.f90"
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
!# 15 "device_memcpy_mod.f90"
!
module device_memcpy_m
  implicit none
!# 19 "device_memcpy_mod.f90"
!# 1 "./device_memcpy_interf.f90"
!
!# 2 "./device_memcpy_interf.f90"
!# 1 "/workspace/develop/q-e/external/devxlib/include/device_macros.h"
!# 3 "./device_memcpy_interf.f90"
!# 3 "./device_memcpy_interf.f90"
!
interface dev_memcpy
    !
    subroutine sp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
      use iso_fortran_env
      implicit none
      !   
      real(real32), intent(inout) :: array_out(:)
      real(real32), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
!# 18 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_r1d
    !
    subroutine sp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      use iso_fortran_env
      implicit none
      !   
      real(real32), intent(inout) :: array_out(:,:)
      real(real32), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
!# 34 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_r2d
    !
    subroutine sp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      use iso_fortran_env
      implicit none
      !   
      real(real32), intent(inout) :: array_out(:,:,:)
      real(real32), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 51 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_r3d
    !
    subroutine sp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      use iso_fortran_env
      implicit none
      !   
      real(real32), intent(inout) :: array_out(:,:,:,:)
      real(real32), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 69 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_r4d
    !
    subroutine dp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
      use iso_fortran_env
      implicit none
      !   
      real(real64), intent(inout) :: array_out(:)
      real(real64), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
!# 84 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_r1d
    !
    subroutine dp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      use iso_fortran_env
      implicit none
      !   
      real(real64), intent(inout) :: array_out(:,:)
      real(real64), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
!# 100 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_r2d
    !
    subroutine dp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      use iso_fortran_env
      implicit none
      !   
      real(real64), intent(inout) :: array_out(:,:,:)
      real(real64), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 117 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_r3d
    !
    subroutine dp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      use iso_fortran_env
      implicit none
      !   
      real(real64), intent(inout) :: array_out(:,:,:,:)
      real(real64), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 135 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_r4d
    !
    subroutine sp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
      use iso_fortran_env
      implicit none
      !   
      complex(real32), intent(inout) :: array_out(:)
      complex(real32), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
!# 150 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_c1d
    !
    subroutine sp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      use iso_fortran_env
      implicit none
      !   
      complex(real32), intent(inout) :: array_out(:,:)
      complex(real32), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
!# 166 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_c2d
    !
    subroutine sp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      use iso_fortran_env
      implicit none
      !   
      complex(real32), intent(inout) :: array_out(:,:,:)
      complex(real32), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 183 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_c3d
    !
    subroutine sp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      use iso_fortran_env
      implicit none
      !   
      complex(real32), intent(inout) :: array_out(:,:,:,:)
      complex(real32), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 201 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memcpy_c4d
    !
    subroutine dp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
      use iso_fortran_env
      implicit none
      !   
      complex(real64), intent(inout) :: array_out(:)
      complex(real64), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
!# 216 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_c1d
    !
    subroutine dp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      use iso_fortran_env
      implicit none
      !   
      complex(real64), intent(inout) :: array_out(:,:)
      complex(real64), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
!# 232 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_c2d
    !
    subroutine dp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      use iso_fortran_env
      implicit none
      !   
      complex(real64), intent(inout) :: array_out(:,:,:)
      complex(real64), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 249 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_c3d
    !
    subroutine dp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      use iso_fortran_env
      implicit none
      !   
      complex(real64), intent(inout) :: array_out(:,:,:,:)
      complex(real64), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 267 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memcpy_c4d
    !
    subroutine i4_dev_memcpy_i1d(array_out, array_in, &
                                            range1, lbound1 )
      use iso_fortran_env
      implicit none
      !   
      integer(int32), intent(inout) :: array_out(:)
      integer(int32), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
!# 282 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memcpy_i1d
    !
    subroutine i4_dev_memcpy_i2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      use iso_fortran_env
      implicit none
      !   
      integer(int32), intent(inout) :: array_out(:,:)
      integer(int32), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
!# 298 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memcpy_i2d
    !
    subroutine i4_dev_memcpy_i3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      use iso_fortran_env
      implicit none
      !   
      integer(int32), intent(inout) :: array_out(:,:,:)
      integer(int32), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 315 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memcpy_i3d
    !
    subroutine i4_dev_memcpy_i4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      use iso_fortran_env
      implicit none
      !   
      integer(int32), intent(inout) :: array_out(:,:,:,:)
      integer(int32), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 333 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memcpy_i4d
    !
    !
!# 609 "./device_memcpy_interf.f90"
    !
!# 1393 "./device_memcpy_interf.f90"
    !
end interface dev_memcpy
!
interface dev_memcpy_async
    !
    subroutine sp_memcpy_d2h_async_r1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 1403 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:)
       real(real32), intent(in)    :: array_in(:)
!# 1411 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1413 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2)
       integer, optional, intent(in) :: lbound1
!# 1418 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_r1d
    !
    subroutine sp_memcpy_d2h_async_r2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 1427 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:)
       real(real32), intent(in)    :: array_in(:,:)
!# 1435 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1437 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2)
       integer, optional, intent(in) :: lbound1,lbound2
!# 1442 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_r2d
    !
    subroutine sp_memcpy_d2h_async_r3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 1452 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:,:)
       real(real32), intent(in)    :: array_in(:,:,:)
!# 1460 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1462 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3
!# 1467 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_r3d
    !
    subroutine sp_memcpy_d2h_async_r4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 1478 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:,:,:)
       real(real32), intent(in)    :: array_in(:,:,:,:)
!# 1486 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1488 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2),range4(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3,lbound4
!# 1493 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_r4d
    !
    subroutine dp_memcpy_d2h_async_r1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 1501 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:)
       real(real64), intent(in)    :: array_in(:)
!# 1509 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1511 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2)
       integer, optional, intent(in) :: lbound1
!# 1516 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_r1d
    !
    subroutine dp_memcpy_d2h_async_r2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 1525 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:)
       real(real64), intent(in)    :: array_in(:,:)
!# 1533 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1535 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2)
       integer, optional, intent(in) :: lbound1,lbound2
!# 1540 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_r2d
    !
    subroutine dp_memcpy_d2h_async_r3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 1550 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:,:)
       real(real64), intent(in)    :: array_in(:,:,:)
!# 1558 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1560 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3
!# 1565 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_r3d
    !
    subroutine dp_memcpy_d2h_async_r4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 1576 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:,:,:)
       real(real64), intent(in)    :: array_in(:,:,:,:)
!# 1584 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1586 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2),range4(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3,lbound4
!# 1591 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_r4d
    !
    subroutine sp_memcpy_d2h_async_c1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 1599 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:)
       complex(real32), intent(in)    :: array_in(:)
!# 1607 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1609 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2)
       integer, optional, intent(in) :: lbound1
!# 1614 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_c1d
    !
    subroutine sp_memcpy_d2h_async_c2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 1623 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:)
       complex(real32), intent(in)    :: array_in(:,:)
!# 1631 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1633 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2)
       integer, optional, intent(in) :: lbound1,lbound2
!# 1638 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_c2d
    !
    subroutine sp_memcpy_d2h_async_c3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 1648 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:,:)
       complex(real32), intent(in)    :: array_in(:,:,:)
!# 1656 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1658 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3
!# 1663 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_c3d
    !
    subroutine sp_memcpy_d2h_async_c4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 1674 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:,:,:)
       complex(real32), intent(in)    :: array_in(:,:,:,:)
!# 1682 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1684 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2),range4(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3,lbound4
!# 1689 "./device_memcpy_interf.f90"
       !
    end subroutine sp_memcpy_d2h_async_c4d
    !
    subroutine dp_memcpy_d2h_async_c1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 1697 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:)
       complex(real64), intent(in)    :: array_in(:)
!# 1705 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1707 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2)
       integer, optional, intent(in) :: lbound1
!# 1712 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_c1d
    !
    subroutine dp_memcpy_d2h_async_c2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 1721 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:)
       complex(real64), intent(in)    :: array_in(:,:)
!# 1729 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1731 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2)
       integer, optional, intent(in) :: lbound1,lbound2
!# 1736 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_c2d
    !
    subroutine dp_memcpy_d2h_async_c3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 1746 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:,:)
       complex(real64), intent(in)    :: array_in(:,:,:)
!# 1754 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1756 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3
!# 1761 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_c3d
    !
    subroutine dp_memcpy_d2h_async_c4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 1772 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:,:,:)
       complex(real64), intent(in)    :: array_in(:,:,:,:)
!# 1780 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1782 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2),range4(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3,lbound4
!# 1787 "./device_memcpy_interf.f90"
       !
    end subroutine dp_memcpy_d2h_async_c4d
    !
    subroutine i4_memcpy_d2h_async_i1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 1795 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:)
       integer(int32), intent(in)    :: array_in(:)
!# 1803 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1805 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2)
       integer, optional, intent(in) :: lbound1
!# 1810 "./device_memcpy_interf.f90"
       !
    end subroutine i4_memcpy_d2h_async_i1d
    !
    subroutine i4_memcpy_d2h_async_i2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 1819 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:)
       integer(int32), intent(in)    :: array_in(:,:)
!# 1827 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1829 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2)
       integer, optional, intent(in) :: lbound1,lbound2
!# 1834 "./device_memcpy_interf.f90"
       !
    end subroutine i4_memcpy_d2h_async_i2d
    !
    subroutine i4_memcpy_d2h_async_i3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 1844 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:,:)
       integer(int32), intent(in)    :: array_in(:,:,:)
!# 1852 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1854 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3
!# 1859 "./device_memcpy_interf.f90"
       !
    end subroutine i4_memcpy_d2h_async_i3d
    !
    subroutine i4_memcpy_d2h_async_i4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 1870 "./device_memcpy_interf.f90"
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:,:,:)
       integer(int32), intent(in)    :: array_in(:,:,:,:)
!# 1878 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 1880 "./device_memcpy_interf.f90"
       integer, optional, intent(in) :: range1(2),range2(2),range3(2),range4(2)
       integer, optional, intent(in) :: lbound1,lbound2,lbound3,lbound4
!# 1885 "./device_memcpy_interf.f90"
       !
    end subroutine i4_memcpy_d2h_async_i4d
    !
    !
!# 2382 "./device_memcpy_interf.f90"
    !
end interface dev_memcpy_async
!
interface 
    !
    subroutine dev_stream_sync(stream)
!# 2391 "./device_memcpy_interf.f90"
       implicit none
!# 2395 "./device_memcpy_interf.f90"
       integer, intent(in) :: stream
!# 2397 "./device_memcpy_interf.f90"
    !
    end subroutine dev_stream_sync
    !
end interface 
!
interface dev_memset
    !
    subroutine sp_dev_memset_r1d(array_out, val, &
                                            
                                            range1, lbound1 )
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:)
       real(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 2417 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_r1d
    !
    subroutine sp_dev_memset_r2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:)
       real(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 2434 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_r2d
    !
    subroutine sp_dev_memset_r3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:,:)
       real(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2452 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_r3d
    !
    subroutine sp_dev_memset_r4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       use iso_fortran_env
       implicit none
       !
       real(real32), intent(inout) :: array_out(:,:,:,:)
       real(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2471 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_r4d
    !
    subroutine dp_dev_memset_r1d(array_out, val, &
                                            
                                            range1, lbound1 )
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:)
       real(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 2487 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_r1d
    !
    subroutine dp_dev_memset_r2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:)
       real(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 2504 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_r2d
    !
    subroutine dp_dev_memset_r3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:,:)
       real(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2522 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_r3d
    !
    subroutine dp_dev_memset_r4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       use iso_fortran_env
       implicit none
       !
       real(real64), intent(inout) :: array_out(:,:,:,:)
       real(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2541 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_r4d
    !
    subroutine sp_dev_memset_c1d(array_out, val, &
                                            
                                            range1, lbound1 )
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:)
       complex(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 2557 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_c1d
    !
    subroutine sp_dev_memset_c2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:)
       complex(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 2574 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_c2d
    !
    subroutine sp_dev_memset_c3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:,:)
       complex(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2592 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_c3d
    !
    subroutine sp_dev_memset_c4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       use iso_fortran_env
       implicit none
       !
       complex(real32), intent(inout) :: array_out(:,:,:,:)
       complex(real32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2611 "./device_memcpy_interf.f90"
       !
    end subroutine sp_dev_memset_c4d
    !
    subroutine dp_dev_memset_c1d(array_out, val, &
                                            
                                            range1, lbound1 )
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:)
       complex(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 2627 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_c1d
    !
    subroutine dp_dev_memset_c2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:)
       complex(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 2644 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_c2d
    !
    subroutine dp_dev_memset_c3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:,:)
       complex(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2662 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_c3d
    !
    subroutine dp_dev_memset_c4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       use iso_fortran_env
       implicit none
       !
       complex(real64), intent(inout) :: array_out(:,:,:,:)
       complex(real64), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2681 "./device_memcpy_interf.f90"
       !
    end subroutine dp_dev_memset_c4d
    !
    subroutine i4_dev_memset_i1d(array_out, val, &
                                            
                                            range1, lbound1 )
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:)
       integer(int32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 2697 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memset_i1d
    !
    subroutine i4_dev_memset_i2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:)
       integer(int32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 2714 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memset_i2d
    !
    subroutine i4_dev_memset_i3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:,:)
       integer(int32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2732 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memset_i3d
    !
    subroutine i4_dev_memset_i4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       use iso_fortran_env
       implicit none
       !
       integer(int32), intent(inout) :: array_out(:,:,:,:)
       integer(int32), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2751 "./device_memcpy_interf.f90"
       !
    end subroutine i4_dev_memset_i4d
    !
    !
!# 3049 "./device_memcpy_interf.f90"
    !
end interface dev_memset
!
!# 20 "device_memcpy_mod.f90"
!# 21 "device_memcpy_mod.f90"
end module device_memcpy_m
