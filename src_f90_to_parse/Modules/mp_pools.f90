!# 1 "mp_pools.f90"
!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_pools
  !----------------------------------------------------------------------------
  !! Pool groups (processors within a pool of k-points).  
  !! Subdivision of image group, used for k-point parallelization.
  !
  USE mp, ONLY : mp_barrier, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER :: npool       = 1
  !! number of "k-points"-pools
  INTEGER :: nproc_pool  = 1
  !! number of processors within a pool
  INTEGER :: me_pool     = 0
  !! index of the processor within a pool 
  INTEGER :: root_pool   = 0
  !! index of the root processor within a pool
  INTEGER :: my_pool_id  = 0
  !! index of my pool
  INTEGER :: inter_pool_comm  = 0
  !! inter pool communicator
  INTEGER :: intra_pool_comm  = 0
  !! intra pool communicator
  !
  INTEGER :: kunit = 1
  !! granularity of k-point distribution.
  !! kunit=1 standard case. In phonon k and k+q must
  !! be on the same pool, so kunit=2.
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_pools( npool_, parent_comm )
    !---------------------------------------------------------------------------
    !! Divide processors (of the "parent_comm" group) into "pools".
    !! Requires: \(\text{npool_}\) read from command line, 
    !!           \(\text{parent_comm}\), typically \(\text{world_comm} = 
    !! \text{group}\) of all processors.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npool_, parent_comm
    !
    INTEGER :: parent_nproc = 1, parent_mype  = 0
    !
!# 95 "mp_pools.f90"
    !
    RETURN
  END SUBROUTINE mp_start_pools
  !
END MODULE mp_pools
!# 103 "mp_pools.f90"
!----------------------------------------------------------------------------
MODULE mp_orthopools
  !----------------------------------------------------------------------------
  !! Ortho-pool groups. Each orthopool group collects the (n+1)th CPU of each
  !! pool, i.e.:
  !
  !! * orthopool 0: first CPU of each pool;
  !! * orthopool 1: second CPU of each pool.
  !
  USE mp, ONLY : mp_barrier, mp_size, mp_rank, mp_comm_split
  USE mp_pools
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER :: northopool      = 1
  !! number of "k-points"-orthopools, must be equal to nproc_pool
  INTEGER :: nproc_orthopool = 1
  !! number of processors within a orthopool, must be equal to npool
  INTEGER :: me_orthopool    = 0
  !! index of the processor within a orthopool, 
  !! must be equal to the pool id of that cpu
  INTEGER :: root_orthopool  = 0
  !! index of the root processor within a orthopool
  INTEGER :: my_orthopool_id = 0
  !! index of my orthopool
  INTEGER :: inter_orthopool_comm = 0
  !! inter orthopool communicator
  INTEGER :: intra_orthopool_comm = 0
  !! intra orthopool communicator
  !
  LOGICAL,PRIVATE :: init_orthopools = .false.
  ! 
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_stop_orthopools( )
    !! Free the orthopools communicators (if they had been set up).
    !
    USE mp, ONLY : mp_comm_free
    IMPLICIT NONE
    ! 
    IF(init_orthopools) THEN
      CALL mp_comm_free ( inter_orthopool_comm )
      CALL mp_comm_free ( intra_orthopool_comm )
      init_orthopools = .false.
    ENDIF
    !
    RETURN
  END SUBROUTINE
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_orthopools( parent_comm )
    !---------------------------------------------------------------------------
    !! Divide processors (of the "parent_comm" group) into "orthopools".  
    !! Requires: pools being already initialized,
    !!           \(\text{parent_comm}\), typically world_comm = group of all processors
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: parent_comm
    !
    INTEGER :: parent_nproc = 1, parent_mype  = 0
    !
    ! Only init this once (I put this check because initialisation 
    ! of orthopools is done later, during EXX bootstrap, not at the beginning
    IF(init_orthopools) RETURN
    init_orthopools = .true.
    !
!# 202 "mp_pools.f90"
    !
    RETURN
  END SUBROUTINE mp_start_orthopools
  !
END MODULE mp_orthopools
