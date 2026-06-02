!# 1 "device_memcpy.f90"
!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contain a CUF KERNEL to perform the selected operation
! cu_memcpy contain also wrappers for cuda_memcpy (sync and async) functions
!
!# 12 "device_memcpy.f90"
!# 1 "/workspace/develop/q-e/external/devxlib/include/device_macros.h"
!# 13 "device_memcpy.f90"
!# 13 "device_memcpy.f90"
!
!=======================================
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
!# 28 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
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
!# 63 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
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
!# 113 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
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
!# 178 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
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
!# 254 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
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
!# 289 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
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
!# 339 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
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
!# 404 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
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
!# 480 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
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
!# 515 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
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
!# 565 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
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
!# 630 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
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
!# 706 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
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
!# 741 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
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
!# 791 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
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
!# 856 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
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
!# 932 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
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
!# 967 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
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
!# 1017 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
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
!# 1082 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine i4_dev_memcpy_i4d
!
!
!======================
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
!# 1162 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
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
!# 1196 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
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
!# 1243 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
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
!# 1303 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
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
!# 1372 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
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
!# 1406 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
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
!# 1453 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
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
!# 1513 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
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
!# 1582 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
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
!# 1616 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
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
!# 1663 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
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
!# 1723 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
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
!# 1792 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
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
!# 1826 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
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
!# 1873 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
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
!# 1933 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
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
!# 2002 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
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
!# 2036 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
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
!# 2083 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
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
!# 2143 "device_memcpy.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine i4_dev_memset_i4d
!
!
!=======================================
!
!# 4228 "device_memcpy.f90"
!
!======================
!# 5414 "device_memcpy.f90"
!
!======================
!
!# 6591 "device_memcpy.f90"
!======================
!
!# 7897 "device_memcpy.f90"
!
!======================
!
subroutine dev_stream_sync(stream)
!# 7904 "device_memcpy.f90"
    implicit none
!# 7910 "device_memcpy.f90"
    integer, intent(in) :: stream
    return
!# 7913 "device_memcpy.f90"
end subroutine dev_stream_sync
!
!======================
!
!
subroutine sp_memcpy_d2h_async_r1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 7923 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real32), intent(inout) :: array_out(:)
    real(real32), intent(in)    :: array_in(:)
!# 7931 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 7933 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 7939 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
!# 7963 "device_memcpy.f90"
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
!# 7966 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_r1d
!
subroutine sp_memcpy_d2h_async_r2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 7975 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real32), intent(inout) :: array_out(:,:)
    real(real32), intent(in)    :: array_in(:,:)
!# 7983 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 7985 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 7991 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
!# 8026 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
!# 8029 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_r2d
!
subroutine sp_memcpy_d2h_async_r3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 8039 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real32), intent(inout) :: array_out(:,:,:)
    real(real32), intent(in)    :: array_in(:,:,:)
!# 8047 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8049 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 8055 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
!# 8101 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
!# 8104 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_r3d
!
subroutine sp_memcpy_d2h_async_r4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 8115 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real32), intent(inout) :: array_out(:,:,:,:)
    real(real32), intent(in)    :: array_in(:,:,:,:)
!# 8123 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8125 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 8131 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
!# 8188 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
!# 8191 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_r4d
!
subroutine dp_memcpy_d2h_async_r1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 8199 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real64), intent(inout) :: array_out(:)
    real(real64), intent(in)    :: array_in(:)
!# 8207 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8209 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 8215 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
!# 8239 "device_memcpy.f90"
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
!# 8242 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_r1d
!
subroutine dp_memcpy_d2h_async_r2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 8251 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real64), intent(inout) :: array_out(:,:)
    real(real64), intent(in)    :: array_in(:,:)
!# 8259 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8261 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 8267 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
!# 8302 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
!# 8305 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_r2d
!
subroutine dp_memcpy_d2h_async_r3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 8315 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real64), intent(inout) :: array_out(:,:,:)
    real(real64), intent(in)    :: array_in(:,:,:)
!# 8323 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8325 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 8331 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
!# 8377 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
!# 8380 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_r3d
!
subroutine dp_memcpy_d2h_async_r4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 8391 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    real(real64), intent(inout) :: array_out(:,:,:,:)
    real(real64), intent(in)    :: array_in(:,:,:,:)
!# 8399 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8401 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 8407 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
!# 8464 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
!# 8467 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_r4d
!
subroutine sp_memcpy_d2h_async_c1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 8475 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real32), intent(inout) :: array_out(:)
    complex(real32), intent(in)    :: array_in(:)
!# 8483 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8485 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 8491 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
!# 8515 "device_memcpy.f90"
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
!# 8518 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_c1d
!
subroutine sp_memcpy_d2h_async_c2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 8527 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real32), intent(inout) :: array_out(:,:)
    complex(real32), intent(in)    :: array_in(:,:)
!# 8535 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8537 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 8543 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
!# 8578 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
!# 8581 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_c2d
!
subroutine sp_memcpy_d2h_async_c3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 8591 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real32), intent(inout) :: array_out(:,:,:)
    complex(real32), intent(in)    :: array_in(:,:,:)
!# 8599 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8601 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 8607 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
!# 8653 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
!# 8656 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_c3d
!
subroutine sp_memcpy_d2h_async_c4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 8667 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real32), intent(inout) :: array_out(:,:,:,:)
    complex(real32), intent(in)    :: array_in(:,:,:,:)
!# 8675 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8677 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 8683 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
!# 8740 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
!# 8743 "device_memcpy.f90"
    !
end subroutine sp_memcpy_d2h_async_c4d
!
subroutine dp_memcpy_d2h_async_c1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 8751 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real64), intent(inout) :: array_out(:)
    complex(real64), intent(in)    :: array_in(:)
!# 8759 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8761 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 8767 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
!# 8791 "device_memcpy.f90"
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
!# 8794 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_c1d
!
subroutine dp_memcpy_d2h_async_c2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 8803 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real64), intent(inout) :: array_out(:,:)
    complex(real64), intent(in)    :: array_in(:,:)
!# 8811 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8813 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 8819 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
!# 8854 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
!# 8857 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_c2d
!
subroutine dp_memcpy_d2h_async_c3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 8867 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real64), intent(inout) :: array_out(:,:,:)
    complex(real64), intent(in)    :: array_in(:,:,:)
!# 8875 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8877 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 8883 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
!# 8929 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
!# 8932 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_c3d
!
subroutine dp_memcpy_d2h_async_c4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 8943 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    complex(real64), intent(inout) :: array_out(:,:,:,:)
    complex(real64), intent(in)    :: array_in(:,:,:,:)
!# 8951 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 8953 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 8959 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
!# 9016 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
!# 9019 "device_memcpy.f90"
    !
end subroutine dp_memcpy_d2h_async_c4d
!
subroutine i4_memcpy_d2h_async_i1d(array_out, array_in, stream, &
                                             range1, lbound1  )
!# 9027 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    integer(int32), intent(inout) :: array_out(:)
    integer(int32), intent(in)    :: array_in(:)
!# 9035 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 9037 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 9043 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
!# 9067 "device_memcpy.f90"
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
!# 9070 "device_memcpy.f90"
    !
end subroutine i4_memcpy_d2h_async_i1d
!
subroutine i4_memcpy_d2h_async_i2d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
!# 9079 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    integer(int32), intent(inout) :: array_out(:,:)
    integer(int32), intent(in)    :: array_in(:,:)
!# 9087 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 9089 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 9095 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
!# 9130 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
!# 9133 "device_memcpy.f90"
    !
end subroutine i4_memcpy_d2h_async_i2d
!
subroutine i4_memcpy_d2h_async_i3d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
!# 9143 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    integer(int32), intent(inout) :: array_out(:,:,:)
    integer(int32), intent(in)    :: array_in(:,:,:)
!# 9151 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 9153 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 9159 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
!# 9205 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
!# 9208 "device_memcpy.f90"
    !
end subroutine i4_memcpy_d2h_async_i3d
!
subroutine i4_memcpy_d2h_async_i4d(array_out, array_in, stream, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
!# 9219 "device_memcpy.f90"
    use iso_fortran_env
    implicit none
    !
    integer(int32), intent(inout) :: array_out(:,:,:,:)
    integer(int32), intent(in)    :: array_in(:,:,:,:)
!# 9227 "device_memcpy.f90"
    integer, intent(in) :: stream
!# 9229 "device_memcpy.f90"
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 9235 "device_memcpy.f90"
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
!# 9292 "device_memcpy.f90"
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
!# 9295 "device_memcpy.f90"
    !
end subroutine i4_memcpy_d2h_async_i4d
!
