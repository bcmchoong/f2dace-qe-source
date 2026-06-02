!# 1 "mp_images.f90"
!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_images
  !----------------------------------------------------------------------------
  !! Image groups (processors within an image). Images are used for
  !! coarse-grid parallelization of semi-independent calculations,
  !! e.g. points along the reaction path (NEB) or phonon irreps.
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE io_global, ONLY : ionode, ionode_id
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER :: nimage = 1
  !! number of images
  INTEGER :: nproc_image=1
  !! number of processors within an image
  INTEGER :: me_image  = 0
  !! index of the processor within an image
  INTEGER :: root_image= 0
  !! index of the root processor within an image
  INTEGER :: my_image_id=0
  !! index of my image
  INTEGER :: inter_image_comm = 0
  !! inter image communicator
  INTEGER :: intra_image_comm = 0
  !! intra image communicator
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE mp_start_images ( nimage_, parent_comm )
    !-----------------------------------------------------------------------
    !! Divide processors (of the "parent_comm" group) into "images". 
    !! Requires: \(\text{nimage_}\), read from command line and
    !! \(\text{parent_comm}\), typically world_comm = group of all processors
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nimage_, parent_comm
    !
!# 99 "mp_images.f90"
    RETURN
    !
  END SUBROUTINE mp_start_images
  !
  SUBROUTINE mp_init_image ( parent_comm )
    !
    !! There is just one image: set it to the same as parent_comm (world)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: parent_comm
    !
    intra_image_comm = parent_comm 
    nproc_image = mp_size( parent_comm )
    me_image    = mp_rank( parent_comm )
    !
    ! ... no need to set inter_image_comm,  my_image_id, root_image
    ! ... set processor that performs I/O
    !
    ionode = ( me_image == root_image )
    ionode_id = root_image
    !
  END SUBROUTINE mp_init_image
  !
END MODULE mp_images
