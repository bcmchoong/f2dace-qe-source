!# 1 "mp_world.f90"
!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_world
  !----------------------------------------------------------------------------
  !! World group - all QE routines using \(\texttt{mp_world_start}\) to 
  !! start MPI will work in the communicator passed as input to 
  !! \(\texttt{mp_world_start}\)
  !
  USE mp, ONLY : mp_barrier, mp_start, mp_end, mp_stop, mp_count_nodes 
  USE io_global, ONLY : meta_ionode_id, meta_ionode
!# 20 "mp_world.f90"
  !
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... World group - all QE routines using mp_world_start to start MPI
  ! ... will work in the communicator passed as input to mp_world_start
  !
  INTEGER :: nnode = 1
  !! number of nodes
  INTEGER :: nproc = 1
  !! number of processors
  INTEGER :: mpime = 0
  !! processor index (starts from 0 to nproc-1)
  INTEGER :: root  = 0
  !! index of the root processor
  INTEGER :: world_comm = 0
  !! communicator
  !
!# 47 "mp_world.f90"
  !
  PRIVATE
  PUBLIC :: nnode, nproc, mpime, root, world_comm, mp_world_start, mp_world_end
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_start ( my_world_comm )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: my_world_comm
    INTEGER :: color, key, ndev
!# 66 "mp_world.f90"
    !
    world_comm = my_world_comm
    !
    ! ... check if mpi is already initialized (library mode) or not
    ! 
!# 83 "mp_world.f90"
    !
    CALL mp_count_nodes ( nnode, color, key, world_comm )
    !
!# 103 "mp_world.f90"
    !
    CALL mp_start( nproc, mpime, world_comm )
    !
    !
    ! ... meta_ionode is true if this processor is the root processor
    ! ... of the world group - "ionode_world" would be a better name
    ! ... meta_ionode_id is the index of such processor
    !
    meta_ionode = ( mpime == root )
    meta_ionode_id = root
    !
    RETURN
    !
  END SUBROUTINE mp_world_start
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_world_end ( )
    !-----------------------------------------------------------------------
!# 124 "mp_world.f90"
    !
    CALL mp_barrier( world_comm )
    CALL mp_end ( world_comm )
!# 133 "mp_world.f90"
    !
  END SUBROUTINE mp_world_end
  !
END MODULE mp_world
