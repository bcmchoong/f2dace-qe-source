!# 1 "mp.f90"
!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains interfaces to most low-level MPI operations:
! initialization and stopping, broadcast, parallel sum, etc.
!
!------------------------------------------------------------------------------!
MODULE mp
!------------------------------------------------------------------------------!
  USE util_param,     ONLY : DP, stdout, i8b
!# 20 "mp.f90"
  !
  IMPLICIT NONE
  PRIVATE
  ! 
  PUBLIC :: mp_start, mp_abort, mp_stop, mp_end, &
    mp_bcast, mp_sum, mp_max, mp_min, mp_rank, mp_size, &
    mp_gather, mp_alltoall, mp_get, mp_put, &
    mp_barrier, mp_report, mp_group_free, &
    mp_root_sum, mp_comm_free, mp_comm_create, mp_comm_group, &
    mp_group_create, mp_comm_split, mp_set_displs, &
    mp_circular_shift_left, mp_circular_shift_left_start, &
    mp_get_comm_null, mp_get_comm_self, mp_count_nodes, &
    mp_type_create_column_section, mp_type_create_row_section, mp_type_free, &
    mp_allgather, mp_waitall, mp_testall
  !
  INTERFACE mp_bcast
    MODULE PROCEDURE mp_bcast_i1, mp_bcast_r1, mp_bcast_c1, &
      mp_bcast_z, mp_bcast_zv, &
      mp_bcast_iv, mp_bcast_i8v, mp_bcast_rv, mp_bcast_cv, mp_bcast_l, mp_bcast_rm, &
      mp_bcast_cm, mp_bcast_im, mp_bcast_it, mp_bcast_i4d, mp_bcast_rt, mp_bcast_lv, &
      mp_bcast_lm, mp_bcast_r4d, mp_bcast_r5d, mp_bcast_ct,  mp_bcast_c4d,&
      mp_bcast_c5d, mp_bcast_c6d
!# 50 "mp.f90"
  END INTERFACE
  ! 
  INTERFACE mp_sum
    MODULE PROCEDURE mp_sum_i1, mp_sum_iv, mp_sum_i8v, mp_sum_im, mp_sum_it, mp_sum_i4, mp_sum_i5, &
      mp_sum_r1, mp_sum_rv, mp_sum_rm, mp_sum_rm1_nc, mp_sum_rm2_nc, mp_sum_rt, mp_sum_r4d, &
      mp_sum_c1, mp_sum_cv, mp_sum_cm, mp_sum_cm1_nc, mp_sum_cm2_nc, mp_sum_ct, mp_sum_c4d, &
      mp_sum_c5d, mp_sum_c6d, mp_sum_rmm, mp_sum_cmm, mp_sum_r5d, &
      mp_sum_r6d
!# 65 "mp.f90"
  END INTERFACE
  ! 
  INTERFACE mp_root_sum
    MODULE PROCEDURE mp_root_sum_rm, mp_root_sum_cm
!# 72 "mp.f90"
  END INTERFACE
  ! 
  INTERFACE mp_get
    MODULE PROCEDURE mp_get_r1, mp_get_rv, mp_get_cv, mp_get_i1, mp_get_iv, mp_get_rm, mp_get_cm
!# 80 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_put
     MODULE PROCEDURE mp_put_rv, mp_put_cv, mp_put_i1, mp_put_iv, &
       mp_put_rm
!# 89 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_max
     MODULE PROCEDURE mp_max_i, mp_max_r, mp_max_rv, mp_max_iv
!# 96 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_min
     MODULE PROCEDURE mp_min_i, mp_min_r, mp_min_rv, mp_min_iv
!# 103 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_gather
     MODULE PROCEDURE mp_gather_i1, mp_gather_iv, mp_gatherv_rv, mp_gatherv_iv, &
                      mp_gatherv_rm, mp_gatherv_im, mp_gatherv_cv, &
                      mp_gatherv_inplace_cplx_array
!# 113 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_allgather
     MODULE PROCEDURE mp_allgatherv_inplace_cplx_array
     MODULE PROCEDURE mp_allgatherv_inplace_real_array
!# 122 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_alltoall
     MODULE PROCEDURE mp_alltoall_c3d, mp_alltoall_i3d
!# 129 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_circular_shift_left
     MODULE PROCEDURE mp_circular_shift_left_i0, &
       mp_circular_shift_left_i1, &
       mp_circular_shift_left_i2, &
       mp_circular_shift_left_r2d, &
       mp_circular_shift_left_c2d
!# 144 "mp.f90"
   END INTERFACE
   ! 
   INTERFACE mp_circular_shift_left_start
     MODULE PROCEDURE mp_circular_shift_left_start_i0, &
       mp_circular_shift_left_start_i1, &
       mp_circular_shift_left_start_i2, &
       mp_circular_shift_left_start_r2d, &
       mp_circular_shift_left_start_c2d
   END INTERFACE
   ! 
   INTERFACE mp_type_create_column_section
     MODULE PROCEDURE mp_type_create_cplx_column_section
     MODULE PROCEDURE mp_type_create_real_column_section
!# 161 "mp.f90"
   END INTERFACE
!# 163 "mp.f90"
   INTERFACE mp_type_create_row_section
     MODULE PROCEDURE mp_type_create_cplx_row_section
     MODULE PROCEDURE mp_type_create_real_row_section
!# 170 "mp.f90"
   END INTERFACE
!------------------------------------------------------------------------------!
!
   CONTAINS
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!..mp_gather_i1
      SUBROUTINE mp_gather_i1(mydata, alldata, root, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata, root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:)
        INTEGER :: ierr
!# 194 "mp.f90"
        alldata(1) = mydata
!# 196 "mp.f90"
        RETURN
      END SUBROUTINE mp_gather_i1
!# 199 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gather_iv
!..Carlo Cavazzoni
      SUBROUTINE mp_gather_iv(mydata, alldata, root, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: mydata(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER, INTENT(OUT) :: alldata(:,:)
        INTEGER :: msglen, ierr
!# 219 "mp.f90"
        msglen = SIZE(mydata)
        IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop( 8014 )
        alldata(:,1) = mydata(:)
!# 223 "mp.f90"
        RETURN
      END SUBROUTINE mp_gather_iv
!# 226 "mp.f90"
!
!------------------------------------------------------------------------------!
!..mp_start
      SUBROUTINE mp_start(numtask, taskid, group)
!# 231 "mp.f90"
! ...
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (OUT) :: numtask, taskid
        INTEGER, INTENT (IN)  :: group
        INTEGER :: ierr
! ...
        ierr = 0
        numtask = 1
        taskid = 0
!# 254 "mp.f90"
        RETURN
      END SUBROUTINE mp_start
!
!------------------------------------------------------------------------------!
!..mp_abort
!# 260 "mp.f90"
      SUBROUTINE mp_abort(errorcode,gid)
        IMPLICIT NONE
        INTEGER :: ierr
        INTEGER, INTENT(IN):: errorcode, gid
!# 271 "mp.f90"
      END SUBROUTINE mp_abort
!
!------------------------------------------------------------------------------!
!..mp_end
!# 276 "mp.f90"
      SUBROUTINE mp_end(groupid, cleanup)
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: groupid
        LOGICAL, OPTIONAL, INTENT(IN) :: cleanup
        INTEGER :: ierr, taskid
        LOGICAL :: cleanup_
!# 284 "mp.f90"
        ierr = 0
        taskid = 0
!# 298 "mp.f90"
        RETURN
      END SUBROUTINE mp_end
!# 301 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_group
!# 304 "mp.f90"
      SUBROUTINE mp_comm_group( comm, group )
         USE parallel_include
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: comm
         INTEGER, INTENT (OUT) :: group
         INTEGER :: ierr
         ierr = 0
!# 315 "mp.f90"
         group = 0
!# 317 "mp.f90"
      END SUBROUTINE  mp_comm_group
!# 319 "mp.f90"
      SUBROUTINE mp_comm_split( old_comm, color, key, new_comm )
         USE parallel_include
         IMPLICIT NONE
         INTEGER, INTENT (IN) :: old_comm
         INTEGER, INTENT (IN) :: color, key
         INTEGER, INTENT (OUT) :: new_comm
         INTEGER :: ierr
         ierr = 0
!# 331 "mp.f90"
         new_comm = old_comm
!# 333 "mp.f90"
      END SUBROUTINE  mp_comm_split
!# 336 "mp.f90"
      SUBROUTINE mp_group_create( group_list, group_size, old_grp, new_grp )
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: group_list(:), group_size, old_grp
        INTEGER, INTENT (OUT) :: new_grp
        INTEGER :: ierr
!# 343 "mp.f90"
        ierr = 0
        new_grp = old_grp
!# 349 "mp.f90"
      END SUBROUTINE mp_group_create
!# 351 "mp.f90"
!------------------------------------------------------------------------------!
      SUBROUTINE mp_comm_create( old_comm, new_grp, new_comm )
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: old_comm
        INTEGER, INTENT (IN) :: new_grp
        INTEGER, INTENT (OUT) :: new_comm
        INTEGER :: ierr
!# 360 "mp.f90"
        ierr = 0
        new_comm = old_comm
!# 366 "mp.f90"
      END SUBROUTINE mp_comm_create
!# 368 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_group_free
      SUBROUTINE mp_group_free( group )
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: group
        INTEGER :: ierr
        ierr = 0
!# 380 "mp.f90"
      END SUBROUTINE mp_group_free
!------------------------------------------------------------------------------!
!# 383 "mp.f90"
      SUBROUTINE mp_comm_free( comm )
         USE parallel_include
         IMPLICIT NONE
         INTEGER, INTENT (INOUT) :: comm
         INTEGER :: ierr
         ierr = 0
!# 395 "mp.f90"
         RETURN
      END SUBROUTINE mp_comm_free
!# 398 "mp.f90"
!------------------------------------------------------------------------------!
! non-blocking helpers
! waits till all request are completed
      SUBROUTINE mp_waitall(requests)
! ...
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: requests(:)
        INTEGER :: ierr
!# 410 "mp.f90"
        ierr = 0
!# 415 "mp.f90"
        RETURN
      END SUBROUTINE mp_waitall
!# 418 "mp.f90"
!tests all requests
      SUBROUTINE mp_testall(requests, flag)
      ! ...
         USE parallel_include
         IMPLICIT NONE
          INTEGER, INTENT (INOUT) :: requests(:)
          INTEGER :: ierr
!# 428 "mp.f90"
          LOGICAL, INTENT(OUT):: flag
          !
          ierr = 0
          flag = .FALSE.
!# 436 "mp.f90"
          flag = .TRUE.
!# 438 "mp.f90"
          RETURN
       END SUBROUTINE mp_testall
!# 441 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_bcast
      SUBROUTINE mp_bcast_i1(msg,source,gid)
        IMPLICIT NONE
        INTEGER :: msg
        INTEGER :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
!# 456 "mp.f90"
      END SUBROUTINE mp_bcast_i1
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_iv(msg, source, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! Bcast an integer vector
      !!  
      IMPLICIT NONE
      ! 
      INTEGER :: msg(:)
      INTEGER, INTENT(in) :: source
      INTEGER, INTENT(in) :: gid
!# 474 "mp.f90"
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_bcast_iv
      !------------------------------------------------------------------------------!
      ! 
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_i8v(msg, source, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! Bcast an integer vector of kind i8b. 
      !!  
      IMPLICIT NONE
      ! 
      INTEGER(KIND = i8b) :: msg(:)
      INTEGER, INTENT(in) :: source
      INTEGER, INTENT(in) :: gid
!# 494 "mp.f90"
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_bcast_i8v
      !------------------------------------------------------------------------------!
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_im(msg, source, gid)
      !------------------------------------------------------------------------------!
        IMPLICIT NONE
        INTEGER :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 510 "mp.f90"
      END SUBROUTINE mp_bcast_im
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_it( msg, source, gid )
        IMPLICIT NONE
        INTEGER :: msg(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 526 "mp.f90"
      END SUBROUTINE mp_bcast_it
!
!------------------------------------------------------------------------------!
!
! Samuel Ponce
!
      SUBROUTINE mp_bcast_i4d(msg, source, gid)
        IMPLICIT NONE
        INTEGER :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 542 "mp.f90"
      END SUBROUTINE mp_bcast_i4d
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_r1( msg, source, gid )
        IMPLICIT NONE
        REAL (DP) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 556 "mp.f90"
      END SUBROUTINE mp_bcast_r1
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rv(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 570 "mp.f90"
      END SUBROUTINE mp_bcast_rv
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_rm(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 584 "mp.f90"
      END SUBROUTINE mp_bcast_rm
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_rt(msg,source,gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 600 "mp.f90"
      END SUBROUTINE mp_bcast_rt
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_r4d(msg, source, gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 616 "mp.f90"
      END SUBROUTINE mp_bcast_r4d
!# 618 "mp.f90"
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_r5d(msg, source, gid)
        IMPLICIT NONE
        REAL (DP) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 633 "mp.f90"
      END SUBROUTINE mp_bcast_r5d
!# 635 "mp.f90"
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_c1(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 647 "mp.f90"
      END SUBROUTINE mp_bcast_c1
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cv(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 660 "mp.f90"
      END SUBROUTINE mp_bcast_cv
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_cm(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 673 "mp.f90"
      END SUBROUTINE mp_bcast_cm
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_ct(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 686 "mp.f90"
      END SUBROUTINE mp_bcast_ct
!# 688 "mp.f90"
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_bcast_c4d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 700 "mp.f90"
      END SUBROUTINE mp_bcast_c4d
!# 702 "mp.f90"
      SUBROUTINE mp_bcast_c5d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 712 "mp.f90"
      END SUBROUTINE mp_bcast_c5d
!# 714 "mp.f90"
      SUBROUTINE mp_bcast_c6d(msg,source,gid)
        IMPLICIT NONE
        COMPLEX (DP) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 724 "mp.f90"
      END SUBROUTINE mp_bcast_c6d
!# 726 "mp.f90"
!
!------------------------------------------------------------------------------!
!# 729 "mp.f90"
      SUBROUTINE mp_bcast_l(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 739 "mp.f90"
      END SUBROUTINE mp_bcast_l
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_lv(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 755 "mp.f90"
      END SUBROUTINE mp_bcast_lv
!# 757 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_bcast_lm(msg,source,gid)
        IMPLICIT NONE
        LOGICAL :: msg(:,:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
!# 771 "mp.f90"
      END SUBROUTINE mp_bcast_lm
!# 774 "mp.f90"
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_z(msg,source,gid)
        IMPLICIT NONE
        CHARACTER (len=*) :: msg
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, ierr, i
        INTEGER, ALLOCATABLE :: imsg(:)
!# 801 "mp.f90"
      END SUBROUTINE mp_bcast_z
!
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
!
      SUBROUTINE mp_bcast_zv(msg,source,gid)
        IMPLICIT NONE
        CHARACTER (len=*) :: msg(:)
        INTEGER, INTENT(IN) :: source
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: msglen, m1, m2, ierr, i, j
        INTEGER, ALLOCATABLE :: imsg(:,:)
!# 837 "mp.f90"
      END SUBROUTINE mp_bcast_zv
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_i1(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include 
        IMPLICIT NONE
        INTEGER :: msg_dest
        INTEGER, INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 854 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen = 1
!# 861 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 863 "mp.f90"
        msglen = 0
!# 865 "mp.f90"
        IF(dest .NE. sour) THEN
!# 880 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF
!# 891 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_i1
!# 894 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_iv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER :: msg_dest(:)
        INTEGER, INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 909 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 916 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 918 "mp.f90"
        msglen = 0
!# 920 "mp.f90"
        IF(sour .NE. dest) THEN
!# 934 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 942 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_iv
!# 945 "mp.f90"
!------------------------------------------------------------------------------!
!# 947 "mp.f90"
      SUBROUTINE mp_get_r1(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL (DP)             :: msg_dest
        REAL (DP), INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 958 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 965 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 967 "mp.f90"
        msglen = 0
!# 969 "mp.f90"
        IF(sour .NE. dest) THEN
!# 983 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF
!# 991 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_r1
!# 994 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include 
        IMPLICIT NONE
        REAL (DP)             :: msg_dest(:)
        REAL (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1009 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1016 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1018 "mp.f90"
        msglen = 0
!# 1020 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1034 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 1042 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_rv
!# 1045 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_rm(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include 
        IMPLICIT NONE
        REAL (DP) :: msg_dest(:,:)
        REAL (DP), INTENT(IN) :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1060 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1067 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1069 "mp.f90"
        msglen = 0
!# 1071 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1085 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
!# 1093 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_rm
!# 1097 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_get_cv(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        COMPLEX (DP)             :: msg_dest(:)
        COMPLEX (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1112 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1119 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1121 "mp.f90"
        msglen = 0
!# 1123 "mp.f90"
        IF( dest .NE. sour ) THEN
!# 1137 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 1145 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_cv
!# 1150 "mp.f90"
!------------------------------------------------------------------------------!
!
! Marco Govoni
!
      SUBROUTINE mp_get_cm(msg_dest, msg_sour, mpime, dest, sour, ip, gid)
        USE parallel_include 
        IMPLICIT NONE
        COMPLEX (DP)              :: msg_dest(:,:)
        COMPLEX (DP), INTENT(IN)  :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1165 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1172 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1174 "mp.f90"
        msglen = 0
!# 1176 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1190 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1), 1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE( msg_sour )
        END IF
!# 1198 "mp.f90"
        RETURN
      END SUBROUTINE mp_get_cm
!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!
!# 1206 "mp.f90"
      SUBROUTINE mp_put_i1(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        USE parallel_include 
        IMPLICIT NONE
        INTEGER :: msg_dest
        INTEGER, INTENT(IN) :: msg_sour
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1217 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1224 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1226 "mp.f90"
        msglen = 0
!# 1228 "mp.f90"
        IF(dest .NE. sour) THEN
!# 1242 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest = msg_sour
          msglen = 1
        END IF
!# 1250 "mp.f90"
        RETURN
      END SUBROUTINE mp_put_i1
!# 1253 "mp.f90"
!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_iv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER             :: msg_dest(:)
        INTEGER, INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1267 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1272 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1274 "mp.f90"
        msglen = 0
!# 1276 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1290 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 1298 "mp.f90"
        RETURN
      END SUBROUTINE mp_put_iv
!# 1301 "mp.f90"
!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL (DP)             :: msg_dest(:)
        REAL (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1315 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1320 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1322 "mp.f90"
        msglen = 0
!# 1324 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1338 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 1346 "mp.f90"
        RETURN
      END SUBROUTINE mp_put_rv
!# 1349 "mp.f90"
!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_rm(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL (DP)             :: msg_dest(:,:)
        REAL (DP), INTENT(IN) :: msg_sour(:,:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1363 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1368 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1370 "mp.f90"
        msglen = 0
!# 1372 "mp.f90"
        IF(sour .NE. dest) THEN
!# 1386 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour,1),1:SIZE(msg_sour,2)) = msg_sour(:,:)
          msglen = SIZE(msg_sour)
        END IF
!# 1394 "mp.f90"
        RETURN
      END SUBROUTINE mp_put_rm
!# 1398 "mp.f90"
!------------------------------------------------------------------------------!
!
!
      SUBROUTINE mp_put_cv(msg_dest, msg_sour, mpime, sour, dest, ip, gid)
        USE parallel_include
        IMPLICIT NONE
        COMPLEX (DP)             :: msg_dest(:)
        COMPLEX (DP), INTENT(IN) :: msg_sour(:)
        INTEGER, INTENT(IN) :: dest, sour, ip, mpime
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
!# 1412 "mp.f90"
        INTEGER :: ierr, nrcv
        INTEGER :: msglen
!# 1417 "mp.f90"
        ! processors not taking part in the communication have 0 length message
!# 1419 "mp.f90"
        msglen = 0
!# 1421 "mp.f90"
        IF( dest .NE. sour ) THEN
!# 1435 "mp.f90"
        ELSEIF(mpime .EQ. sour)THEN
          msg_dest(1:SIZE(msg_sour)) = msg_sour(:)
          msglen = SIZE(msg_sour)
        END IF
!# 1443 "mp.f90"
        RETURN
      END SUBROUTINE mp_put_cv
!# 1446 "mp.f90"
!
!------------------------------------------------------------------------------!
!
!..mp_stop
!
      SUBROUTINE mp_stop(code)
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: code
        INTEGER :: ierr
        WRITE( stdout, fmt='( "*** error in Message Passing (mp) module ***")' )
        WRITE( stdout, fmt='( "*** error code: ",I5)' ) code
!# 1462 "mp.f90"
        STOP
      END SUBROUTINE mp_stop
!------------------------------------------------------------------------------!
!
!..mp_sum
      SUBROUTINE mp_sum_i1(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 1476 "mp.f90"
      END SUBROUTINE mp_sum_i1
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_iv(msg, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! MPI sum an integer vector from all cores and bcast the result to all.  
      !! 
      IMPLICIT NONE
      ! 
      INTEGER, INTENT(inout) :: msg(:)
      INTEGER, INTENT(in) :: gid
!# 1493 "mp.f90"
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_sum_iv
      !------------------------------------------------------------------------------!
      ! 
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_i8(msg, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! MPI sum an integer vector from all cores and bcast the result to all.  
      !! 
      IMPLICIT NONE
      ! 
      INTEGER(KIND = i8b), INTENT(inout) :: msg
      INTEGER, INTENT(in) :: gid
!# 1512 "mp.f90"
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_sum_i8
      !------------------------------------------------------------------------------!
      ! 
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_i8v(msg, gid)
      !------------------------------------------------------------------------------!
      !! 
      !! MPI sum an integer vector from all cores and bcast the result to all.  
      !! 
      IMPLICIT NONE
      ! 
      INTEGER(KIND = i8b), INTENT(inout) :: msg(:)
      INTEGER, INTENT(in) :: gid
!# 1531 "mp.f90"
      !------------------------------------------------------------------------------!
      END SUBROUTINE mp_sum_i8v
      !------------------------------------------------------------------------------!
      !
      !------------------------------------------------------------------------------!
      SUBROUTINE mp_sum_im(msg,gid)
      !------------------------------------------------------------------------------!
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT(IN) :: gid
!# 1546 "mp.f90"
      END SUBROUTINE mp_sum_im
!
!------------------------------------------------------------------------------!
!# 1550 "mp.f90"
      SUBROUTINE mp_sum_it(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT (IN) :: gid
!# 1559 "mp.f90"
      END SUBROUTINE mp_sum_it
!# 1561 "mp.f90"
!------------------------------------------------------------------------------!
!# 1563 "mp.f90"
      SUBROUTINE mp_sum_i4(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT (IN) :: gid
!# 1572 "mp.f90"
      END SUBROUTINE mp_sum_i4
!# 1574 "mp.f90"
!------------------------------------------------------------------------------!
!# 1576 "mp.f90"
      SUBROUTINE mp_sum_i5(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT (IN) :: gid
!# 1585 "mp.f90"
      END SUBROUTINE mp_sum_i5
!# 1588 "mp.f90"
!------------------------------------------------------------------------------!
!# 1590 "mp.f90"
      SUBROUTINE mp_sum_r1(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT (IN) :: gid
!# 1599 "mp.f90"
      END SUBROUTINE mp_sum_r1
!# 1601 "mp.f90"
!
!------------------------------------------------------------------------------!
!# 1604 "mp.f90"
      SUBROUTINE mp_sum_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT (IN) :: gid
!# 1613 "mp.f90"
      END SUBROUTINE mp_sum_rv
!
!------------------------------------------------------------------------------!
!# 1618 "mp.f90"
      SUBROUTINE mp_sum_rm(msg, gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
!# 1627 "mp.f90"
      END SUBROUTINE mp_sum_rm
!# 1629 "mp.f90"
      SUBROUTINE mp_sum_rm1_nc(msg, k1, k2, gid)
        ! for non-contiguous 1D arrays  
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT (IN) :: gid
        INTEGER, INTENT (IN) :: k1, k2
!# 1646 "mp.f90"
      END SUBROUTINE mp_sum_rm1_nc
!# 1649 "mp.f90"
      SUBROUTINE mp_sum_rm2_nc(msg, k1, k2, k3, k4, gid)
        ! for non-contiguous 2D arrays  
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
        INTEGER, INTENT (IN) :: k1, k2, k3, k4
!# 1666 "mp.f90"
      END SUBROUTINE mp_sum_rm2_nc
!# 1668 "mp.f90"
      SUBROUTINE mp_root_sum_rm( msg, res, root, gid )
        USE parallel_include
        IMPLICIT NONE
        REAL (DP), INTENT (IN)  :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,   INTENT (IN) :: gid
!# 1691 "mp.f90"
        res = msg
!# 1695 "mp.f90"
      END SUBROUTINE mp_root_sum_rm
!# 1698 "mp.f90"
      SUBROUTINE mp_root_sum_cm( msg, res, root, gid )
        USE parallel_include
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN)  :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER,   INTENT (IN)  :: root
        INTEGER,  INTENT (IN) :: gid
!# 1721 "mp.f90"
        res = msg
!# 1725 "mp.f90"
      END SUBROUTINE mp_root_sum_cm
!# 1727 "mp.f90"
!
!------------------------------------------------------------------------------!
!# 1731 "mp.f90"
!------------------------------------------------------------------------------!
!
!# 1734 "mp.f90"
      SUBROUTINE mp_sum_rmm( msg, res, root, gid )
        USE parallel_include
        IMPLICIT NONE
        REAL (DP), INTENT (IN) :: msg(:,:)
        REAL (DP), INTENT (OUT) :: res(:,:)
        INTEGER, INTENT (IN) :: root
        INTEGER, INTENT (IN) :: gid
        INTEGER :: group
        INTEGER :: msglen
        INTEGER :: taskid, ierr
!# 1745 "mp.f90"
        msglen = size(msg)
!# 1762 "mp.f90"
        res = msg
!# 1765 "mp.f90"
      END SUBROUTINE mp_sum_rmm
!# 1768 "mp.f90"
!
!------------------------------------------------------------------------------!
!# 1772 "mp.f90"
      SUBROUTINE mp_sum_rt( msg, gid )
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1781 "mp.f90"
      END SUBROUTINE mp_sum_rt
!# 1783 "mp.f90"
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r4d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1797 "mp.f90"
      END SUBROUTINE mp_sum_r4d
!# 1801 "mp.f90"
!------------------------------------------------------------------------------!
!# 1803 "mp.f90"
      SUBROUTINE mp_sum_c1(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 1812 "mp.f90"
      END SUBROUTINE mp_sum_c1
!
!------------------------------------------------------------------------------!
!# 1816 "mp.f90"
      SUBROUTINE mp_sum_cv(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
!# 1825 "mp.f90"
      END SUBROUTINE mp_sum_cv
!
!------------------------------------------------------------------------------!
!# 1829 "mp.f90"
      SUBROUTINE mp_sum_cm(msg, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
!# 1838 "mp.f90"
      END SUBROUTINE mp_sum_cm
!
!------------------------------------------------------------------------------!
!# 1843 "mp.f90"
      SUBROUTINE mp_sum_cm1_nc(msg, k1, k2, gid)
        ! for non-contiguous 1D arrays  
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT (IN) :: gid
        INTEGER, INTENT (IN) :: k1, k2
!# 1860 "mp.f90"
      END SUBROUTINE mp_sum_cm1_nc
!# 1862 "mp.f90"
      SUBROUTINE mp_sum_cm2_nc(msg, k1, k2, k3, k4, gid)
        ! for non-contiguous 2D arrays  
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:)
        INTEGER, INTENT (IN) :: gid
        INTEGER, INTENT (IN) :: k1, k2, k3, k4
!# 1879 "mp.f90"
      END SUBROUTINE mp_sum_cm2_nc
!
!------------------------------------------------------------------------------!
!# 1883 "mp.f90"
      SUBROUTINE mp_sum_cmm(msg, res, gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (IN) :: msg(:,:)
        COMPLEX (DP), INTENT (OUT) :: res(:,:)
        INTEGER, INTENT (IN) :: gid
!# 1893 "mp.f90"
        res = msg
!# 1895 "mp.f90"
      END SUBROUTINE mp_sum_cmm
!# 1898 "mp.f90"
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_ct(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1912 "mp.f90"
      END SUBROUTINE mp_sum_ct
!# 1914 "mp.f90"
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c4d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1928 "mp.f90"
      END SUBROUTINE mp_sum_c4d
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c5d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1943 "mp.f90"
      END SUBROUTINE mp_sum_c5d
!# 1945 "mp.f90"
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_r5d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1958 "mp.f90"
      END SUBROUTINE mp_sum_r5d
!# 1961 "mp.f90"
      SUBROUTINE mp_sum_r6d(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1970 "mp.f90"
      END SUBROUTINE mp_sum_r6d
!# 1972 "mp.f90"
!
!------------------------------------------------------------------------------!
!
! Carlo Cavazzoni
!
      SUBROUTINE mp_sum_c6d(msg,gid)
        IMPLICIT NONE
        COMPLEX (DP), INTENT (INOUT) :: msg(:,:,:,:,:,:)
        INTEGER, INTENT(IN) :: gid
!# 1986 "mp.f90"
      END SUBROUTINE mp_sum_c6d
!# 1990 "mp.f90"
!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 2000 "mp.f90"
      END SUBROUTINE mp_max_i
!
!------------------------------------------------------------------------------!
!
!..mp_max_iv
!..Carlo Cavazzoni
!
      SUBROUTINE mp_max_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
!# 2016 "mp.f90"
      END SUBROUTINE mp_max_iv
!
!----------------------------------------------------------------------
!# 2020 "mp.f90"
      SUBROUTINE mp_max_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 2029 "mp.f90"
      END SUBROUTINE mp_max_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_max_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
!# 2041 "mp.f90"
      END SUBROUTINE mp_max_rv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_i(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 2052 "mp.f90"
      END SUBROUTINE mp_min_i
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_iv(msg,gid)
        IMPLICIT NONE
        INTEGER, INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
!# 2063 "mp.f90"
      END SUBROUTINE mp_min_iv
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_r(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg
        INTEGER, INTENT(IN) :: gid
!# 2074 "mp.f90"
      END SUBROUTINE mp_min_r
!
!------------------------------------------------------------------------------!
      SUBROUTINE mp_min_rv(msg,gid)
        IMPLICIT NONE
        REAL (DP), INTENT (INOUT) :: msg(:)
        INTEGER, INTENT(IN) :: gid
!# 2086 "mp.f90"
      END SUBROUTINE mp_min_rv
!# 2088 "mp.f90"
!------------------------------------------------------------------------------!
!# 2090 "mp.f90"
      SUBROUTINE mp_barrier(gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: gid
        INTEGER :: ierr
!# 2099 "mp.f90"
      END SUBROUTINE mp_barrier
!# 2101 "mp.f90"
!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_rank
      FUNCTION mp_rank( comm )
        USE parallel_include
        IMPLICIT NONE
        INTEGER :: mp_rank
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, taskid
!# 2111 "mp.f90"
        ierr = 0
        taskid = 0
!# 2117 "mp.f90"
        mp_rank = taskid
      END FUNCTION mp_rank
!# 2120 "mp.f90"
!------------------------------------------------------------------------------!
!.. Carlo Cavazzoni
!..mp_size
      FUNCTION mp_size( comm )
        USE parallel_include
        IMPLICIT NONE
        INTEGER :: mp_size
        INTEGER, INTENT(IN) :: comm
        INTEGER :: ierr, numtask
!# 2130 "mp.f90"
        ierr = 0
        numtask = 1
!# 2136 "mp.f90"
        mp_size = numtask
      END FUNCTION mp_size
!# 2139 "mp.f90"
      SUBROUTINE mp_report
        INTEGER :: i
        WRITE( stdout, *)
!# 2148 "mp.f90"
        WRITE( stdout, *)
!# 2150 "mp.f90"
        RETURN
      END SUBROUTINE mp_report
!# 2154 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_rv
!..Carlo Cavazzoni
!# 2158 "mp.f90"
      SUBROUTINE mp_gatherv_rv( mydata, alldata, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL(DP) :: mydata(:)
        REAL(DP) :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid
!# 2185 "mp.f90"
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
!# 2190 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_rv
!# 2193 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_cv
!..Carlo Cavazzoni
!# 2197 "mp.f90"
      SUBROUTINE mp_gatherv_cv( mydata, alldata, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        COMPLEX(DP) :: mydata(:)
        COMPLEX(DP) :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid
!# 2224 "mp.f90"
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
!# 2229 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_cv
!# 2232 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_rv
!..Carlo Cavazzoni
!# 2236 "mp.f90"
      SUBROUTINE mp_gatherv_iv( mydata, alldata, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER :: mydata(:)
        INTEGER :: alldata(:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid
!# 2263 "mp.f90"
        IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
!# 2268 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_iv
!# 2272 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_rm
!..Carlo Cavazzoni
!# 2276 "mp.f90"
      SUBROUTINE mp_gatherv_rm( mydata, alldata, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL(DP) :: mydata(:,:)  ! Warning first dimension is supposed constant!
        REAL(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)
!# 2314 "mp.f90"
        IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( alldata, 2 ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata, 2  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
!# 2320 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_rm
!# 2323 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_im
!..Carlo Cavazzoni
!# 2327 "mp.f90"
      SUBROUTINE mp_gatherv_im( mydata, alldata, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        INTEGER :: mydata(:,:)  ! Warning first dimension is supposed constant!
        INTEGER :: alldata(:,:)
        INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
        INTEGER, INTENT(IN) :: gid
        INTEGER :: group
        INTEGER :: ierr, npe, myid, nsiz
        INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)
!# 2365 "mp.f90"
        IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( alldata, 2 ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
        IF ( SIZE( mydata, 2  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
        !
        alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
!# 2371 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_im
!# 2375 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_gatherv_inplace_cplx_array
!..Ye Luo
!# 2379 "mp.f90"
      SUBROUTINE mp_gatherv_inplace_cplx_array(alldata, my_column_type, recvcount, displs, root, gid)
        USE parallel_include
        IMPLICIT NONE
        COMPLEX(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: my_column_type
        INTEGER, INTENT(IN) :: recvcount(:), displs(:)
        INTEGER, INTENT(IN) :: root, gid
        INTEGER :: ierr, npe, myid
!# 2405 "mp.f90"
        RETURN
      END SUBROUTINE mp_gatherv_inplace_cplx_array
!# 2408 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_allgatherv_inplace_cplx_array
!..Ye Luo
!# 2412 "mp.f90"
      SUBROUTINE mp_allgatherv_inplace_cplx_array(alldata, my_element_type, recvcount, displs, gid)
        USE parallel_include
        IMPLICIT NONE
        COMPLEX(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: my_element_type
        INTEGER, INTENT(IN) :: recvcount(:), displs(:)
        INTEGER, INTENT(IN) :: gid
        INTEGER :: ierr, npe, myid
!# 2433 "mp.f90"
        RETURN
      END SUBROUTINE mp_allgatherv_inplace_cplx_array
!# 2436 "mp.f90"
!.. SdG added 16/08/19
      SUBROUTINE mp_allgatherv_inplace_real_array(alldata, my_element_type, recvcount, displs, gid)
        USE parallel_include
        IMPLICIT NONE
        REAL(DP) :: alldata(:,:)
        INTEGER, INTENT(IN) :: my_element_type
        INTEGER, INTENT(IN) :: recvcount(:), displs(:)
        INTEGER, INTENT(IN) :: gid
        INTEGER :: ierr, npe, myid
!# 2458 "mp.f90"
        RETURN
      END SUBROUTINE mp_allgatherv_inplace_real_array
!# 2461 "mp.f90"
!------------------------------------------------------------------------------!
!# 2463 "mp.f90"
      SUBROUTINE mp_set_displs( recvcount, displs, ntot, nproc )
        !  Given the number of elements on each processor (recvcount), this subroutine
        !  sets the correct offsets (displs) to collect them on a single
        !  array with contiguous elemets
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: recvcount(:) ! number of elements on each processor
        INTEGER, INTENT(OUT) :: displs(:)   ! offsets/displacements
        INTEGER, INTENT(OUT) :: ntot
        INTEGER, INTENT(IN) :: nproc
        INTEGER :: i
!# 2474 "mp.f90"
        displs( 1 ) = 0
        !
!# 2483 "mp.f90"
        ntot = recvcount( 1 )
!# 2485 "mp.f90"
        RETURN
      END SUBROUTINE mp_set_displs
!# 2488 "mp.f90"
!------------------------------------------------------------------------------!
!# 2491 "mp.f90"
SUBROUTINE mp_alltoall_c3d( sndbuf, rcvbuf, gid )
   USE parallel_include
   IMPLICIT NONE
   COMPLEX(DP) :: sndbuf( :, :, : )
   COMPLEX(DP) :: rcvbuf( :, :, : )
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe
!# 2518 "mp.f90"
   rcvbuf = sndbuf
!# 2522 "mp.f90"
   RETURN
END SUBROUTINE mp_alltoall_c3d
!# 2526 "mp.f90"
!------------------------------------------------------------------------------!
!# 2528 "mp.f90"
SUBROUTINE mp_alltoall_i3d( sndbuf, rcvbuf, gid )
   USE parallel_include
   IMPLICIT NONE
   INTEGER :: sndbuf( :, :, : )
   INTEGER :: rcvbuf( :, :, : )
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe
!# 2555 "mp.f90"
   rcvbuf = sndbuf
!# 2559 "mp.f90"
   RETURN
END SUBROUTINE mp_alltoall_i3d
!# 2562 "mp.f90"
SUBROUTINE mp_circular_shift_left_i0( buf, itag, gid )
   USE parallel_include
   IMPLICIT NONE
   INTEGER :: buf
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2592 "mp.f90"
   ! do nothing
!# 2594 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_i0
!# 2598 "mp.f90"
SUBROUTINE mp_circular_shift_left_i1( buf, itag, gid )
   USE parallel_include
   IMPLICIT NONE
   INTEGER :: buf(:)
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2628 "mp.f90"
   ! do nothing
!# 2630 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_i1
!# 2634 "mp.f90"
SUBROUTINE mp_circular_shift_left_i2( buf, itag, gid )
   USE parallel_include
   IMPLICIT NONE
   INTEGER :: buf(:,:)
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2664 "mp.f90"
   ! do nothing
!# 2666 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_i2
!# 2670 "mp.f90"
SUBROUTINE mp_circular_shift_left_r2d( buf, itag, gid )
   USE parallel_include
   IMPLICIT NONE
   REAL(DP) :: buf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2700 "mp.f90"
   ! do nothing
!# 2702 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_r2d
!# 2705 "mp.f90"
SUBROUTINE mp_circular_shift_left_c2d( buf, itag, gid )
   USE parallel_include
   IMPLICIT NONE
   COMPLEX(DP) :: buf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2735 "mp.f90"
   ! do nothing
!# 2737 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_c2d
!# 2741 "mp.f90"
!------------------------------------------------------------------------------!
!..mp_circular_shift_left_start
SUBROUTINE mp_circular_shift_left_start_i0( sendbuf, recvbuf, itag, gid, requests)
   USE parallel_include
   IMPLICIT NONE
   INTEGER  :: sendbuf, recvbuf
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER, INTENT(INOUT) :: requests(2)
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2780 "mp.f90"
   recvbuf = sendbuf
!# 2783 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_start_i0
!# 2787 "mp.f90"
SUBROUTINE mp_circular_shift_left_start_i1( sendbuf, recvbuf, itag, gid, requests)
   USE parallel_include
   IMPLICIT NONE
   INTEGER  :: sendbuf( : ), recvbuf( : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER, INTENT(INOUT) :: requests(2)
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2825 "mp.f90"
   recvbuf = sendbuf
!# 2828 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_start_i1
!# 2832 "mp.f90"
SUBROUTINE mp_circular_shift_left_start_i2( sendbuf, recvbuf, itag, gid, requests)
   USE parallel_include
   IMPLICIT NONE
   INTEGER  :: sendbuf( :, : ), recvbuf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER, INTENT(INOUT) :: requests(2)
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2870 "mp.f90"
   recvbuf = sendbuf
!# 2873 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_start_i2
!# 2877 "mp.f90"
SUBROUTINE mp_circular_shift_left_start_r2d( sendbuf, recvbuf, itag, gid, requests)
   USE parallel_include
   IMPLICIT NONE
   REAL(DP) :: sendbuf( :, : ), recvbuf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER, INTENT(INOUT) :: requests(2)
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2915 "mp.f90"
   recvbuf = sendbuf
!# 2918 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_start_r2d
!# 2922 "mp.f90"
SUBROUTINE mp_circular_shift_left_start_c2d( sendbuf, recvbuf, itag, gid, requests)
   USE parallel_include
   IMPLICIT NONE
   COMPLEX(DP) :: sendbuf( :, : ), recvbuf( :, : )
   INTEGER, INTENT(IN) :: itag
   INTEGER, INTENT(IN) :: gid
   INTEGER, INTENT(INOUT) :: requests(2)
   INTEGER :: nsiz, group, ierr, npe, sour, dest, mype
!# 2960 "mp.f90"
   recvbuf = sendbuf
!# 2963 "mp.f90"
   RETURN
END SUBROUTINE mp_circular_shift_left_start_c2d
!
!
!------------------------------------------------------------------------------!
!..mp_count_nodes
SUBROUTINE mp_count_nodes(num_nodes, color, key, group)
  !
  ! ... This routine counts the number of nodes using
  ! ...  MPI_GET_PROCESSOR_NAME in the group specified by `group`.
  ! ...  It returns colors and keys to be used in MPI_COMM_SPLIT.
  ! ...  When running in parallel, the evaluation of color and key
  ! ...  is done by all processors.
  ! ...
  ! ...
  ! ... input:
  ! ...    group      Communicator used to count nodes.
  !
  ! ... output:
  ! ...    num_nodes  Number of unique nodes in the communicator
  ! ...    color      Integer (positive), same for all processes residing on a node.
  ! ...    key        Integer, unique number identifying each process on the same node.
  ! ...
  USE parallel_include
  IMPLICIT NONE
  INTEGER, INTENT (OUT) :: num_nodes
  INTEGER, INTENT (OUT) :: color
  INTEGER, INTENT (OUT) :: key
  INTEGER, INTENT (IN)  :: group
!# 2996 "mp.f90"
  
  LOGICAL, ALLOCATABLE   :: found_list(:)
  INTEGER, ALLOCATABLE   :: color_list(:)
  INTEGER, ALLOCATABLE   :: key_list(:)
  !
  INTEGER :: hostname_len, max_hostname_len, numtask, me, ierr
  !
  ! Loops variables
  INTEGER :: i, j, e, s, c, k
  ! ...
  ierr      = 0
  num_nodes = 1
  color     = 1
  key       = 0
  !
!# 3079 "mp.f90"
  RETURN
END SUBROUTINE mp_count_nodes
!
FUNCTION mp_get_comm_null( )
  USE parallel_include
  IMPLICIT NONE
  INTEGER :: mp_get_comm_null
  mp_get_comm_null = MPI_COMM_NULL
END FUNCTION mp_get_comm_null
!# 3089 "mp.f90"
FUNCTION mp_get_comm_self( )
  USE parallel_include
  IMPLICIT NONE
  INTEGER :: mp_get_comm_self
  mp_get_comm_self = MPI_COMM_SELF
END FUNCTION mp_get_comm_self
!# 3096 "mp.f90"
SUBROUTINE mp_type_create_cplx_column_section(dummy, start, length, stride, mytype)
  USE parallel_include
  IMPLICIT NONE
  !
  COMPLEX (DP), INTENT(IN) :: dummy
  INTEGER, INTENT(IN) :: start, length, stride
  INTEGER, INTENT(OUT) :: mytype
  !
!# 3113 "mp.f90"
  mytype = 0;
!# 3115 "mp.f90"
  !
  RETURN
END SUBROUTINE mp_type_create_cplx_column_section
!# 3119 "mp.f90"
SUBROUTINE mp_type_create_real_column_section(dummy, start, length, stride, mytype)
  USE parallel_include
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: dummy
  INTEGER, INTENT(IN) :: start, length, stride
  INTEGER, INTENT(OUT) :: mytype
  !
!# 3136 "mp.f90"
  mytype = 0;
!# 3138 "mp.f90"
  !
  RETURN
END SUBROUTINE mp_type_create_real_column_section
!# 3142 "mp.f90"
SUBROUTINE mp_type_create_cplx_row_section(dummy, column_start, column_stride, row_length, mytype)
  USE parallel_include
  IMPLICIT NONE
  !
  COMPLEX (DP), INTENT(IN) :: dummy
  INTEGER, INTENT(IN) :: column_start, column_stride, row_length
  INTEGER, INTENT(OUT) :: mytype
  !
!# 3170 "mp.f90"
  mytype = 0;
!# 3172 "mp.f90"
  !
  RETURN
END SUBROUTINE mp_type_create_cplx_row_section
!# 3176 "mp.f90"
SUBROUTINE mp_type_create_real_row_section(dummy, column_start, column_stride, row_length, mytype)
  USE parallel_include
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: dummy
  INTEGER, INTENT(IN) :: column_start, column_stride, row_length
  INTEGER, INTENT(OUT) :: mytype
  !
!# 3204 "mp.f90"
  mytype = 0;
!# 3206 "mp.f90"
  !
  RETURN
END SUBROUTINE mp_type_create_real_row_section
!# 3210 "mp.f90"
SUBROUTINE mp_type_free(mytype)
  USE parallel_include
  IMPLICIT NONE
  INTEGER :: mytype, ierr
  !
!# 3219 "mp.f90"
  !
  RETURN
END SUBROUTINE mp_type_free
!------------------------------------------------------------------------------!
!  GPU specific subroutines (Pietro Bonfa')
!------------------------------------------------------------------------------!
! Before hacking on the CUDA part remember that:
!
! 1. all mp_* interface should be blocking with respect to both MPI and CUDA.
!    MPI will only wait for completion on the default stream therefore device
!    synchronization must be enforced.
! 2. Host -> device memory copies of a memory block of 64 KB or less are
!    asynchronous in the sense that they may return before the data is actually
!    available on the GPU. However, the user is still free to change the buffer
!    as soon as those calls return with no ill effects.
!    (https://devtalk.nvidia.com/default/topic/471866/cuda-programming-and-performance/host-device-memory-copies-up-to-64-kb-are-asynchronous/)
! 3. For transfers from device to either pageable or pinned host memory,
!    the function returns only once the copy has completed.
! 4. GPU synchronization is always enforced even if no communication takes place.
!------------------------------------------------------------------------------!
!# 6596 "mp.f90"
!------------------------------------------------------------------------------!
END MODULE mp
!------------------------------------------------------------------------------!
!
! Script to generate stop messages:
!   # coding: utf-8
!   import re
!   import sys
!   i = 8000
!   def replace(match):
!       global i
!       i += 1
!       return 'mp_stop( {0} )'.format(i)
!   
!   with open(sys.argv[1],'r') as f:
!       data = re.sub(r"mp_stop\(\s?\d+\s?\)", replace, f.read())
!       with open(sys.argv[1]+'.new','w') as fo:
!           fo.write(data)
