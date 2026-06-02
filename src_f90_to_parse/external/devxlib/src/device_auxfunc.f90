!# 1 "device_auxfunc.f90"
!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Auxiliary functions
!
!# 10 "device_auxfunc.f90"
!# 1 "/workspace/develop/q-e/external/devxlib/include/device_macros.h"
!# 11 "device_auxfunc.f90"
!# 11 "device_auxfunc.f90"
!
!=======================================
!
!# 15 "device_auxfunc.f90"
subroutine dp_dev_conjg_c1d(array_inout, &
                                    
                                   range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 27 "device_auxfunc.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_inout(i1 ) = &
            conjg (  array_inout (i1 ) )
    enddo
    !
end subroutine dp_dev_conjg_c1d
!
subroutine dp_dev_conjg_c2d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 61 "device_auxfunc.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2 ) = &
            conjg (  array_inout (i1,i2 ) )
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c2d
!
subroutine dp_dev_conjg_c3d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 108 "device_auxfunc.f90"
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
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
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
        array_inout(i1,i2,i3 ) = &
            conjg (  array_inout (i1,i2,i3 ) )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c3d
!
subroutine dp_dev_conjg_c4d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3, & 
                                   range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 168 "device_auxfunc.f90"
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
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_inout, 4)/)
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
        array_inout(i1,i2,i3,i4 ) = &
            conjg (  array_inout (i1,i2,i3,i4 ) )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c4d
!
subroutine sp_dev_conjg_c1d(array_inout, &
                                    
                                   range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 237 "device_auxfunc.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_inout(i1 ) = &
            conjg (  array_inout (i1 ) )
    enddo
    !
end subroutine sp_dev_conjg_c1d
!
subroutine sp_dev_conjg_c2d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 271 "device_auxfunc.f90"
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2 ) = &
            conjg (  array_inout (i1,i2 ) )
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c2d
!
subroutine sp_dev_conjg_c3d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 318 "device_auxfunc.f90"
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
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
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
        array_inout(i1,i2,i3 ) = &
            conjg (  array_inout (i1,i2,i3 ) )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c3d
!
subroutine sp_dev_conjg_c4d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3, & 
                                   range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 378 "device_auxfunc.f90"
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
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_inout, 4)/)
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
        array_inout(i1,i2,i3,i4 ) = &
            conjg (  array_inout (i1,i2,i3,i4 ) )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c4d
!
!
!=======================================
!
!# 832 "device_auxfunc.f90"
subroutine dxlib_errore( calling_routine, message, ierr )
    implicit none
    character(len=*), intent(in) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
    integer,          intent(in) :: ierr
    !
    print *, trim(calling_routine), trim(message), ierr
    !
end subroutine dxlib_errore
!
!======================
!
!# 846 "device_auxfunc.f90"
subroutine dp_dev_vec_upd_remap_r1d(ndim, vout, v1, map1, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim
    real(PRCSN), intent(inout) :: vout(:) 
    real(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    real(PRCSN), optional, intent(in)    :: scal
!# 858 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))
      enddo
    endif
end subroutine dp_dev_vec_upd_remap_r1d
!
subroutine dp_dev_vec_upd_remap_c1d(ndim, vout, v1, map1, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), optional, intent(in)    :: scal
!# 885 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))
      enddo
    endif
end subroutine dp_dev_vec_upd_remap_c1d
!
subroutine sp_dev_vec_upd_remap_r1d(ndim, vout, v1, map1, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim
    real(PRCSN), intent(inout) :: vout(:) 
    real(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    real(PRCSN), optional, intent(in)    :: scal
!# 912 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))
      enddo
    endif
end subroutine sp_dev_vec_upd_remap_r1d
!
subroutine sp_dev_vec_upd_remap_c1d(ndim, vout, v1, map1, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), optional, intent(in)    :: scal
!# 939 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))
      enddo
    endif
end subroutine sp_dev_vec_upd_remap_c1d
!
!
!======================
!
!# 958 "device_auxfunc.f90"
subroutine dp_dev_vec_upd_v_remap_v_r1d(ndim, vout, v1, map1, v2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim
    real(PRCSN), intent(inout) :: vout(:) 
    real(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    real(PRCSN), intent(in)    :: v2(:) 
    real(PRCSN), optional, intent(in)    :: scal
!# 971 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)
      enddo
    endif
end subroutine dp_dev_vec_upd_v_remap_v_r1d
!
subroutine dp_dev_vec_upd_v_remap_v_c1d(ndim, vout, v1, map1, v2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), intent(in)    :: v2(:) 
    complex(PRCSN), optional, intent(in)    :: scal
!# 999 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)
      enddo
    endif
end subroutine dp_dev_vec_upd_v_remap_v_c1d
!
subroutine sp_dev_vec_upd_v_remap_v_r1d(ndim, vout, v1, map1, v2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim
    real(PRCSN), intent(inout) :: vout(:) 
    real(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    real(PRCSN), intent(in)    :: v2(:) 
    real(PRCSN), optional, intent(in)    :: scal
!# 1027 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)
      enddo
    endif
end subroutine sp_dev_vec_upd_v_remap_v_r1d
!
subroutine sp_dev_vec_upd_v_remap_v_c1d(ndim, vout, v1, map1, v2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), intent(in)    :: v2(:) 
    complex(PRCSN), optional, intent(in)    :: scal
!# 1055 "device_auxfunc.f90"
    integer :: i
    !   
    if (present(scal)) then
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)*scal
      enddo
    else
      !$cuf kernel do(1)
      do i = 1, ndim
          vout(i) = v1(map1(i))*v2(i)
      enddo
    endif
end subroutine sp_dev_vec_upd_v_remap_v_c1d
!
!
!======================
!
!# 1075 "device_auxfunc.f90"
subroutine dp_dev_vec_upd_v_remap_v_x_c1d(ndim, vout, v1,op1, map1, v2,op2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), intent(in)    :: v2(:) 
    character(1), intent(in)    :: op1, op2
    complex(PRCSN), optional, intent(in)    :: scal
!# 1089 "device_auxfunc.f90"
    integer :: i
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*v2(i)*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*v2(i)
         enddo
       endif
    elseif (op1=="C".and.op2=="N") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*v2(i)*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*v2(i)
         enddo
       endif
    elseif (op1=="N".and.op2=="C") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*conjg(v2(i))*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*conjg(v2(i))
         enddo
       endif
    elseif (op1=="N".and.op2=="C") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*conjg(v2(i))*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*conjg(v2(i))
         enddo
       endif
    else
      call dxlib_errore("dp_dev_vec_upd_v_remap_v_x_c1d","invalid op1/op2",10)
    endif
    !
end subroutine dp_dev_vec_upd_v_remap_v_x_c1d
!
subroutine sp_dev_vec_upd_v_remap_v_x_c1d(ndim, vout, v1,op1, map1, v2,op2, scal)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim
    complex(PRCSN), intent(inout) :: vout(:) 
    complex(PRCSN), intent(in)    :: v1(:) 
    integer,      intent(in)    :: map1(:) 
    complex(PRCSN), intent(in)    :: v2(:) 
    character(1), intent(in)    :: op1, op2
    complex(PRCSN), optional, intent(in)    :: scal
!# 1159 "device_auxfunc.f90"
    integer :: i
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*v2(i)*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*v2(i)
         enddo
       endif
    elseif (op1=="C".and.op2=="N") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*v2(i)*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*v2(i)
         enddo
       endif
    elseif (op1=="N".and.op2=="C") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*conjg(v2(i))*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = v1(map1(i))*conjg(v2(i))
         enddo
       endif
    elseif (op1=="N".and.op2=="C") then
       if (present(scal)) then
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*conjg(v2(i))*scal
         enddo
       else
         !$cuf kernel do(1)
         do i = 1, ndim
             vout(i) = conjg(v1(map1(i)))*conjg(v2(i))
         enddo
       endif
    else
      call dxlib_errore("sp_dev_vec_upd_v_remap_v_x_c1d","invalid op1/op2",10)
    endif
    !
end subroutine sp_dev_vec_upd_v_remap_v_x_c1d
!
!
!======================
!
subroutine dp_dev_mat_upd_dMd_r2d(ndim1, ndim2, mat, v1,op1, v2,op2, scal)
    !
    ! performs: mat(i,j) = scal * op1(v1(i)) * mat(i,j) * op2(v2(j))
    ! op = 'N', 'R', 'C',       'RC'
    !       x   1/x  conjg(x)   conjg(1/x)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim1,ndim2
    real(PRCSN), intent(inout) :: mat(:,:)
    real(PRCSN), intent(in)    :: v1(:)
    real(PRCSN), intent(in)    :: v2(:)
    character(1), intent(in)    :: op1, op2
    real(PRCSN), optional, intent(in)  :: scal
!# 1235 "device_auxfunc.f90"
    integer :: i,j
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
!# 1240 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1244 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1250 "device_auxfunc.f90"
       else
!# 1252 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1256 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1262 "device_auxfunc.f90"
       endif
    elseif (op1=="R".and.op2=="R") then
       if (present(scal)) then
!# 1266 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1270 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal / v1(i) * mat(i,j) / v2(j)
         enddo
         enddo
!# 1276 "device_auxfunc.f90"
       else
!# 1278 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1282 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = mat(i,j) /v1(i) / v2(j)
         enddo
         enddo
!# 1288 "device_auxfunc.f90"
       endif
    else
      call dxlib_errore("dp_mat_upd_dMd_r2d","invalid op1/op2",10)
    endif
    !
end subroutine dp_dev_mat_upd_dMd_r2d
!
subroutine dp_dev_mat_upd_dMd_c2d(ndim1, ndim2, mat, v1,op1, v2,op2, scal)
    !
    ! performs: mat(i,j) = scal * op1(v1(i)) * mat(i,j) * op2(v2(j))
    ! op = 'N', 'R', 'C',       'RC'
    !       x   1/x  conjg(x)   conjg(1/x)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    integer,      intent(in)    :: ndim1,ndim2
    complex(PRCSN), intent(inout) :: mat(:,:)
    complex(PRCSN), intent(in)    :: v1(:)
    complex(PRCSN), intent(in)    :: v2(:)
    character(1), intent(in)    :: op1, op2
    complex(PRCSN), optional, intent(in)  :: scal
!# 1312 "device_auxfunc.f90"
    integer :: i,j
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
!# 1317 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1321 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1327 "device_auxfunc.f90"
       else
!# 1329 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1333 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1339 "device_auxfunc.f90"
       endif
    elseif (op1=="R".and.op2=="R") then
       if (present(scal)) then
!# 1343 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1347 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal / v1(i) * mat(i,j) / v2(j)
         enddo
         enddo
!# 1353 "device_auxfunc.f90"
       else
!# 1355 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1359 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = mat(i,j) /v1(i) / v2(j)
         enddo
         enddo
!# 1365 "device_auxfunc.f90"
       endif
    elseif (op1=="C".and.op2=="C") then
       if (present(scal)) then
!# 1369 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1373 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * conjg(v1(i)) * mat(i,j) * conjg(v2(j))
         enddo
         enddo
!# 1379 "device_auxfunc.f90"
       else
!# 1381 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1385 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = conjg(v1(i)) * mat(i,j) * conjg(v2(j))
         enddo
         enddo
!# 1391 "device_auxfunc.f90"
       endif
    else
      call dxlib_errore("dp_mat_upd_dMd_c2d","invalid op1/op2",10)
    endif
    !
end subroutine dp_dev_mat_upd_dMd_c2d
!
subroutine sp_dev_mat_upd_dMd_r2d(ndim1, ndim2, mat, v1,op1, v2,op2, scal)
    !
    ! performs: mat(i,j) = scal * op1(v1(i)) * mat(i,j) * op2(v2(j))
    ! op = 'N', 'R', 'C',       'RC'
    !       x   1/x  conjg(x)   conjg(1/x)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim1,ndim2
    real(PRCSN), intent(inout) :: mat(:,:)
    real(PRCSN), intent(in)    :: v1(:)
    real(PRCSN), intent(in)    :: v2(:)
    character(1), intent(in)    :: op1, op2
    real(PRCSN), optional, intent(in)  :: scal
!# 1415 "device_auxfunc.f90"
    integer :: i,j
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
!# 1420 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1424 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1430 "device_auxfunc.f90"
       else
!# 1432 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1436 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1442 "device_auxfunc.f90"
       endif
    elseif (op1=="R".and.op2=="R") then
       if (present(scal)) then
!# 1446 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1450 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal / v1(i) * mat(i,j) / v2(j)
         enddo
         enddo
!# 1456 "device_auxfunc.f90"
       else
!# 1458 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1462 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = mat(i,j) /v1(i) / v2(j)
         enddo
         enddo
!# 1468 "device_auxfunc.f90"
       endif
    else
      call dxlib_errore("sp_mat_upd_dMd_r2d","invalid op1/op2",10)
    endif
    !
end subroutine sp_dev_mat_upd_dMd_r2d
!
subroutine sp_dev_mat_upd_dMd_c2d(ndim1, ndim2, mat, v1,op1, v2,op2, scal)
    !
    ! performs: mat(i,j) = scal * op1(v1(i)) * mat(i,j) * op2(v2(j))
    ! op = 'N', 'R', 'C',       'RC'
    !       x   1/x  conjg(x)   conjg(1/x)
    implicit none
    !   
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    integer,      intent(in)    :: ndim1,ndim2
    complex(PRCSN), intent(inout) :: mat(:,:)
    complex(PRCSN), intent(in)    :: v1(:)
    complex(PRCSN), intent(in)    :: v2(:)
    character(1), intent(in)    :: op1, op2
    complex(PRCSN), optional, intent(in)  :: scal
!# 1492 "device_auxfunc.f90"
    integer :: i,j
    !   
    if (op1=="N".and.op2=="N") then
       if (present(scal)) then
!# 1497 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1501 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1507 "device_auxfunc.f90"
       else
!# 1509 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1513 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = v1(i) * mat(i,j) * v2(j)
         enddo
         enddo
!# 1519 "device_auxfunc.f90"
       endif
    elseif (op1=="R".and.op2=="R") then
       if (present(scal)) then
!# 1523 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1527 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal / v1(i) * mat(i,j) / v2(j)
         enddo
         enddo
!# 1533 "device_auxfunc.f90"
       else
!# 1535 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1539 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = mat(i,j) /v1(i) / v2(j)
         enddo
         enddo
!# 1545 "device_auxfunc.f90"
       endif
    elseif (op1=="C".and.op2=="C") then
       if (present(scal)) then
!# 1549 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1553 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = scal * conjg(v1(i)) * mat(i,j) * conjg(v2(j))
         enddo
         enddo
!# 1559 "device_auxfunc.f90"
       else
!# 1561 "device_auxfunc.f90"
         !$omp parallel do default(shared), private(i,j), collapse(2)
!# 1565 "device_auxfunc.f90"
         do j = 1, ndim2
         do i = 1, ndim1
             mat(i,j) = conjg(v1(i)) * mat(i,j) * conjg(v2(j))
         enddo
         enddo
!# 1571 "device_auxfunc.f90"
       endif
    else
      call dxlib_errore("sp_mat_upd_dMd_c2d","invalid op1/op2",10)
    endif
    !
end subroutine sp_dev_mat_upd_dMd_c2d
!
!
!======================
!
subroutine dp_dev_mem_addscal_r1d(array_out, array_in, scal, &
                                            range1, lbound1 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(14,200)), intent(inout) :: array_out(:)
    real(selected_real_kind(14,200)), intent(in)    :: array_in(:)
    real(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 1594 "device_auxfunc.f90"
    !
    real(selected_real_kind(14,200)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = &
            array_out(i1 ) + &
            scal_*array_in(i1 )
    enddo
    !
end subroutine dp_dev_mem_addscal_r1d
!
subroutine dp_dev_mem_addscal_r2d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(14,200)), intent(inout) :: array_out(:,:)
    real(selected_real_kind(14,200)), intent(in)    :: array_in(:,:)
    real(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 1635 "device_auxfunc.f90"
    !
    real(selected_real_kind(14,200)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = &
            array_out(i1,i2 ) + &
            scal_*array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_r2d
!
subroutine dp_dev_mem_addscal_r3d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(14,200)), intent(inout) :: array_out(:,:,:)
    real(selected_real_kind(14,200)), intent(in)    :: array_in(:,:,:)
    real(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 1691 "device_auxfunc.f90"
    !
    real(selected_real_kind(14,200)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = &
            array_out(i1,i2,i3 ) + &
            scal_*array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_r3d
!
subroutine dp_dev_mem_addscal_r4d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(14,200)), intent(inout) :: array_out(:,:,:,:)
    real(selected_real_kind(14,200)), intent(in)    :: array_in(:,:,:,:)
    real(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 1762 "device_auxfunc.f90"
    !
    real(selected_real_kind(14,200)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = &
            array_out(i1,i2,i3,i4 ) + &
            scal_*array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_r4d
!
subroutine sp_dev_mem_addscal_r1d(array_out, array_in, scal, &
                                            range1, lbound1 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(6, 37)), intent(inout) :: array_out(:)
    real(selected_real_kind(6, 37)), intent(in)    :: array_in(:)
    real(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 1844 "device_auxfunc.f90"
    !
    real(selected_real_kind(6, 37)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = &
            array_out(i1 ) + &
            scal_*array_in(i1 )
    enddo
    !
end subroutine sp_dev_mem_addscal_r1d
!
subroutine sp_dev_mem_addscal_r2d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:)
    real(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:)
    real(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 1885 "device_auxfunc.f90"
    !
    real(selected_real_kind(6, 37)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = &
            array_out(i1,i2 ) + &
            scal_*array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_r2d
!
subroutine sp_dev_mem_addscal_r3d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:,:)
    real(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:,:)
    real(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 1941 "device_auxfunc.f90"
    !
    real(selected_real_kind(6, 37)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = &
            array_out(i1,i2,i3 ) + &
            scal_*array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_r3d
!
subroutine sp_dev_mem_addscal_r4d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    use iso_fortran_env
    implicit none
    !
    real(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:,:,:)
    real(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:,:,:)
    real(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2012 "device_auxfunc.f90"
    !
    real(selected_real_kind(6, 37)) :: scal_ = 1.0
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = &
            array_out(i1,i2,i3,i4 ) + &
            scal_*array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_r4d
!
subroutine dp_dev_mem_addscal_c1d(array_out, array_in, scal, &
                                            range1, lbound1 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(14,200)), intent(inout) :: array_out(:)
    complex(selected_real_kind(14,200)), intent(in)    :: array_in(:)
    complex(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 2094 "device_auxfunc.f90"
    !
    complex(selected_real_kind(14,200)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = &
            array_out(i1 ) + &
            scal_*array_in(i1 )
    enddo
    !
end subroutine dp_dev_mem_addscal_c1d
!
subroutine dp_dev_mem_addscal_c2d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(14,200)), intent(inout) :: array_out(:,:)
    complex(selected_real_kind(14,200)), intent(in)    :: array_in(:,:)
    complex(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 2135 "device_auxfunc.f90"
    !
    complex(selected_real_kind(14,200)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = &
            array_out(i1,i2 ) + &
            scal_*array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_c2d
!
subroutine dp_dev_mem_addscal_c3d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(14,200)), intent(inout) :: array_out(:,:,:)
    complex(selected_real_kind(14,200)), intent(in)    :: array_in(:,:,:)
    complex(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2191 "device_auxfunc.f90"
    !
    complex(selected_real_kind(14,200)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = &
            array_out(i1,i2,i3 ) + &
            scal_*array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_c3d
!
subroutine dp_dev_mem_addscal_c4d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(14,200)), intent(inout) :: array_out(:,:,:,:)
    complex(selected_real_kind(14,200)), intent(in)    :: array_in(:,:,:,:)
    complex(selected_real_kind(14,200)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2262 "device_auxfunc.f90"
    !
    complex(selected_real_kind(14,200)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = &
            array_out(i1,i2,i3,i4 ) + &
            scal_*array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_mem_addscal_c4d
!
subroutine sp_dev_mem_addscal_c1d(array_out, array_in, scal, &
                                            range1, lbound1 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(6, 37)), intent(inout) :: array_out(:)
    complex(selected_real_kind(6, 37)), intent(in)    :: array_in(:)
    complex(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
!# 2344 "device_auxfunc.f90"
    !
    complex(selected_real_kind(6, 37)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = &
            array_out(i1 ) + &
            scal_*array_in(i1 )
    enddo
    !
end subroutine sp_dev_mem_addscal_c1d
!
subroutine sp_dev_mem_addscal_c2d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:)
    complex(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:)
    complex(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
!# 2385 "device_auxfunc.f90"
    !
    complex(selected_real_kind(6, 37)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = &
            array_out(i1,i2 ) + &
            scal_*array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_c2d
!
subroutine sp_dev_mem_addscal_c3d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:,:)
    complex(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:,:)
    complex(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 2441 "device_auxfunc.f90"
    !
    complex(selected_real_kind(6, 37)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = &
            array_out(i1,i2,i3 ) + &
            scal_*array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_c3d
!
subroutine sp_dev_mem_addscal_c4d(array_out, array_in, scal, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    use iso_fortran_env
    implicit none
    !
    complex(selected_real_kind(6, 37)), intent(inout) :: array_out(:,:,:,:)
    complex(selected_real_kind(6, 37)), intent(in)    :: array_in(:,:,:,:)
    complex(selected_real_kind(6, 37)), optional, intent(in) :: scal
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 2512 "device_auxfunc.f90"
    !
    complex(selected_real_kind(6, 37)) :: scal_ = (1.0, 0.0)
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
    if (present(scal)) scal_=scal
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = &
            array_out(i1,i2,i3,i4 ) + &
            scal_*array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_mem_addscal_c4d
!
!
!======================
!
