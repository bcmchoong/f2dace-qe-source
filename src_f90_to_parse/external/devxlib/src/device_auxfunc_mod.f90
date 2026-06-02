!# 1 "device_auxfunc_mod.f90"
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
!# 15 "device_auxfunc_mod.f90"
!
module device_auxfunc_m
  implicit none
!# 19 "device_auxfunc_mod.f90"
!# 1 "./device_auxfunc_interf.f90"
!# 2 "./device_auxfunc_interf.f90"
!# 1 "/workspace/develop/q-e/external/devxlib/include/device_macros.h"
!# 3 "./device_auxfunc_interf.f90"
!# 3 "./device_auxfunc_interf.f90"
!
interface dev_conjg
    !
!# 7 "./device_auxfunc_interf.f90"
    subroutine dp_dev_conjg_c1d(array_inout, &
                                   
                                   range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 19 "./device_auxfunc_interf.f90"
       !
    end subroutine dp_dev_conjg_c1d
    !
    subroutine dp_dev_conjg_c2d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 35 "./device_auxfunc_interf.f90"
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
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 52 "./device_auxfunc_interf.f90"
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
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 70 "./device_auxfunc_interf.f90"
       !
    end subroutine dp_dev_conjg_c4d
    !
    subroutine sp_dev_conjg_c1d(array_inout, &
                                   
                                   range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
!# 85 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_conjg_c1d
    !
    subroutine sp_dev_conjg_c2d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
!# 101 "./device_auxfunc_interf.f90"
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
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
!# 118 "./device_auxfunc_interf.f90"
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
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
!# 136 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_conjg_c4d
    !
    !
!# 252 "./device_auxfunc_interf.f90"
    !
end interface dev_conjg
!# 255 "./device_auxfunc_interf.f90"
interface dev_vec_upd_remap
    !
!# 258 "./device_auxfunc_interf.f90"
    subroutine dp_dev_vec_upd_remap_v_r1d(ndim, vout, v1, map1, v2, scal)
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       integer,      intent(in)    :: ndim
       real(PRCSN), intent(inout) :: vout(:) 
       real(PRCSN), intent(in)    :: v1(:) 
       integer,      intent(in)    :: map1(:) 
       real(PRCSN), intent(in)    :: v2(:) 
       real(PRCSN), optional, intent(in)    :: scal
!# 271 "./device_auxfunc_interf.f90"
       !
    end subroutine dp_dev_vec_upd_remap_v_r1d
    !
    subroutine dp_dev_vec_upd_remap_v_c1d(ndim, vout, v1, map1, v2, scal)
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       integer,      intent(in)    :: ndim
       complex(PRCSN), intent(inout) :: vout(:) 
       complex(PRCSN), intent(in)    :: v1(:) 
       integer,      intent(in)    :: map1(:) 
       complex(PRCSN), intent(in)    :: v2(:) 
       complex(PRCSN), optional, intent(in)    :: scal
!# 287 "./device_auxfunc_interf.f90"
       !
    end subroutine dp_dev_vec_upd_remap_v_c1d
    !
    subroutine sp_dev_vec_upd_remap_v_r1d(ndim, vout, v1, map1, v2, scal)
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       integer,      intent(in)    :: ndim
       real(PRCSN), intent(inout) :: vout(:) 
       real(PRCSN), intent(in)    :: v1(:) 
       integer,      intent(in)    :: map1(:) 
       real(PRCSN), intent(in)    :: v2(:) 
       real(PRCSN), optional, intent(in)    :: scal
!# 303 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_vec_upd_remap_v_r1d
    !
    subroutine sp_dev_vec_upd_remap_v_c1d(ndim, vout, v1, map1, v2, scal)
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       integer,      intent(in)    :: ndim
       complex(PRCSN), intent(inout) :: vout(:) 
       complex(PRCSN), intent(in)    :: v1(:) 
       integer,      intent(in)    :: map1(:) 
       complex(PRCSN), intent(in)    :: v2(:) 
       complex(PRCSN), optional, intent(in)    :: scal
!# 319 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_vec_upd_remap_v_c1d
    !
    !
end interface dev_vec_upd_remap
!# 325 "./device_auxfunc_interf.f90"
interface dev_vec_upd_v_remap_v
    !
!# 328 "./device_auxfunc_interf.f90"
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
!# 341 "./device_auxfunc_interf.f90"
    !
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
!# 357 "./device_auxfunc_interf.f90"
    !
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
!# 373 "./device_auxfunc_interf.f90"
    !
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
!# 389 "./device_auxfunc_interf.f90"
    !
    end subroutine sp_dev_vec_upd_v_remap_v_c1d
    !
    !
!# 395 "./device_auxfunc_interf.f90"
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
!# 409 "./device_auxfunc_interf.f90"
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
!# 426 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_vec_upd_v_remap_v_x_c1d
    !
    !
end interface dev_vec_upd_v_remap_v
!# 432 "./device_auxfunc_interf.f90"
interface dev_mat_upd_dMd
    !
!# 435 "./device_auxfunc_interf.f90"
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
!# 452 "./device_auxfunc_interf.f90"
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
!# 472 "./device_auxfunc_interf.f90"
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
!# 492 "./device_auxfunc_interf.f90"
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
!# 512 "./device_auxfunc_interf.f90"
       !
    end subroutine sp_dev_mat_upd_dMd_c2d
    !
    !
end interface dev_mat_upd_dMd
!# 518 "./device_auxfunc_interf.f90"
interface dev_mem_addscal
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
!# 533 "./device_auxfunc_interf.f90"
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
!# 549 "./device_auxfunc_interf.f90"
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
!# 566 "./device_auxfunc_interf.f90"
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
!# 584 "./device_auxfunc_interf.f90"
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
!# 599 "./device_auxfunc_interf.f90"
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
!# 615 "./device_auxfunc_interf.f90"
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
!# 632 "./device_auxfunc_interf.f90"
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
!# 650 "./device_auxfunc_interf.f90"
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
!# 665 "./device_auxfunc_interf.f90"
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
!# 681 "./device_auxfunc_interf.f90"
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
!# 698 "./device_auxfunc_interf.f90"
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
!# 716 "./device_auxfunc_interf.f90"
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
!# 731 "./device_auxfunc_interf.f90"
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
!# 747 "./device_auxfunc_interf.f90"
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
!# 764 "./device_auxfunc_interf.f90"
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
!# 782 "./device_auxfunc_interf.f90"
    end subroutine sp_dev_mem_addscal_c4d
    !
    !
end interface dev_mem_addscal
!# 20 "device_auxfunc_mod.f90"
!# 21 "device_auxfunc_mod.f90"
end module device_auxfunc_m
