!# 1 "mp_base.f90"
!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!  Wrapper for MPI implementations that have problems with large messages
!
!# 13 "mp_base.f90"
!  In some MPI implementation the communication subsystem
!  crashes when message exceeds a given size, so we need
!  to break down MPI communications in smaller pieces
!
!# 20 "mp_base.f90"
!  Some implementation of MPI (OpenMPI) if it is not well tuned for the given
!  network hardware (InfiniBand) tend to lose performance or get stuck inside
!  collective routines if processors are not well synchronized
!  A barrier fixes the problem
!
!# 28 "mp_base.f90"
!=----------------------------------------------------------------------------=!
!
! These routines allocate buffer spaces used in reduce_base_real_gpu.
! These should be in data_buffer.f90 but need to be here because size
! depends on the 100000 definition
SUBROUTINE allocate_buffers
    USE data_buffer
    IMPLICIT NONE
    INTEGER, PARAMETER :: maxb = 100000
    !
    IF (.NOT. ALLOCATED(mp_buff_r)) ALLOCATE(mp_buff_r(maxb))
    IF (.NOT. ALLOCATED(mp_buff_i)) ALLOCATE(mp_buff_i(maxb))
    !
END SUBROUTINE allocate_buffers
!# 43 "mp_base.f90"
SUBROUTINE deallocate_buffers
    USE data_buffer
    IMPLICIT NONE
    !
    DEALLOCATE(mp_buff_r, mp_buff_i)
    !
END SUBROUTINE deallocate_buffers
!# 51 "mp_base.f90"
!=----------------------------------------------------------------------------=!
!
!# 54 "mp_base.f90"
SUBROUTINE mp_synchronize( gid )
   USE parallel_include  
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: gid
!# 63 "mp_base.f90"
   RETURN
END SUBROUTINE mp_synchronize
!# 67 "mp_base.f90"
!=----------------------------------------------------------------------------=!
!
!# 70 "mp_base.f90"
   SUBROUTINE bcast_real( array, n, root, gid )
        USE util_param, ONLY: DP
        USE parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        REAL(DP) :: array( n )
!# 115 "mp_base.f90"
        RETURN
   END SUBROUTINE bcast_real 
!# 118 "mp_base.f90"
   !------------------------------------------------------------------------------!
   SUBROUTINE bcast_integer(array, n, root, gid)
   !------------------------------------------------------------------------------!
   !! 
   !! Broadcast integers
   !!   
   USE parallel_include  
   ! 
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n, root, gid
   INTEGER :: array(n)
!# 164 "mp_base.f90"
   RETURN
   !------------------------------------------------------------------------------!
   END SUBROUTINE bcast_integer
   !------------------------------------------------------------------------------!
   ! 
   !------------------------------------------------------------------------------!
   SUBROUTINE bcast_integer8(array, n, root, gid)
   !------------------------------------------------------------------------------!
   !! 
   !! Broadcast integers
   !!   
   USE util_param,     ONLY : i8b
   USE parallel_include
   ! 
   IMPLICIT NONE
   INTEGER, INTENT(in) :: n, root, gid
   INTEGER(KIND = i8b) :: array(n)
!# 216 "mp_base.f90"
   RETURN
   !------------------------------------------------------------------------------!
   END SUBROUTINE bcast_integer8
   !------------------------------------------------------------------------------!
   ! 
   !------------------------------------------------------------------------------!
   SUBROUTINE bcast_logical( array, n, root, gid )
        USE parallel_include  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n, root, gid
        LOGICAL :: array( n )
!# 265 "mp_base.f90"
        RETURN
   END SUBROUTINE bcast_logical
!# 269 "mp_base.f90"
!
! ... "reduce"-like subroutines
!
!# 336 "mp_base.f90"
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY: buff => mp_buff_r
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim     ! size of the array
  REAL(DP)                :: ps(dim) ! array whose elements have to be reduced
  INTEGER,  INTENT(IN)    :: comm    ! communicator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 427 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE reduce_base_real
!
!# 433 "mp_base.f90"
!
!
!# 553 "mp_base.f90"
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE util_param,  ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: dim
  INTEGER, INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  !
!# 643 "mp_base.f90"
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_integer8(dim, ps, comm, root)
  !----------------------------------------------------------------------------
  !!
  !! Sums a distributed variable ps(dim) over the processors.
  !! This version uses a fixed-length buffer of appropriate (?) dim
  !!
  USE util_param,  ONLY : DP, i8b
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)    :: dim
  INTEGER(KIND = i8b), INTENT(inout) :: ps(dim)
  INTEGER, INTENT(in)    :: comm    ! communicator
  INTEGER, INTENT(in)    :: root    ! if root <  0 perform a reduction to all procs
                                    ! if root >= 0 perform a reduce only to root proc.
  INTEGER(KIND = i8b):: buff(dim)   ! quick and dirty fix: automatic array
  !
!# 738 "mp_base.f90"
  !
  RETURN
  !----------------------------------------------------------------------------
  END SUBROUTINE reduce_base_integer8
  !----------------------------------------------------------------------------
  !
!# 745 "mp_base.f90"
  !
  ! ... "reduce"-like subroutines
  !
  !----------------------------------------------------------------------------
  SUBROUTINE reduce_base_real_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed variable ps(dim) over the processors,
  ! ... and store the results in variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE util_param, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  REAL(DP), INTENT(IN)  :: ps(dim)
  REAL(DP)              :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                   ! if root >= 0 perform a reduce only to root proc.
  !
!# 829 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE reduce_base_real_to
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE reduce_base_integer_to( dim, ps, psout, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... sums a distributed integer variable ps(dim) over the processors, and
  ! ... saves the result on the output variable psout.
  ! ... This version uses a fixed-length buffer of appropriate (?) length
  !
  USE util_param, ONLY : DP
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: dim
  INTEGER,  INTENT(IN)  :: ps(dim)
  INTEGER               :: psout(dim)
  INTEGER,  INTENT(IN)  :: comm    ! communecator
  INTEGER,  INTENT(IN)  :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 917 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE reduce_base_integer_to
!
!
!  Parallel MIN and MAX
!
!# 926 "mp_base.f90"
!----------------------------------------------------------------------------
SUBROUTINE parallel_min_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER                 :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 1015 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE parallel_min_integer
!# 1020 "mp_base.f90"
!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_integer( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param,  ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_i
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  INTEGER                 :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 1108 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE parallel_max_integer
!# 1114 "mp_base.f90"
!----------------------------------------------------------------------------
SUBROUTINE parallel_min_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the minimum value of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_r
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP)                :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 1201 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE parallel_min_real
!# 1206 "mp_base.f90"
!
!----------------------------------------------------------------------------
SUBROUTINE parallel_max_real( dim, ps, comm, root )
  !----------------------------------------------------------------------------
  !
  ! ... compute the maximum value of a distributed variable ps(dim) over the processors.
  ! ... This version uses a fixed-length buffer of appropriate (?) dim
  !
  USE util_param, ONLY : DP
  USE data_buffer, ONLY : buff => mp_buff_r
  USE parallel_include  
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: dim
  REAL(DP)                :: ps(dim)
  INTEGER,  INTENT(IN)    :: comm    ! communecator
  INTEGER,  INTENT(IN)    :: root    ! if root <  0 perform a reduction to all procs
                                     ! if root >= 0 perform a reduce only to root proc.
  !
!# 1296 "mp_base.f90"
  !
  RETURN
  !
END SUBROUTINE parallel_max_real
