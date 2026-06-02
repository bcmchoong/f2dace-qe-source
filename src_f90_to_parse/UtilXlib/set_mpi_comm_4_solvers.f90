!# 1 "set_mpi_comm_4_solvers.f90"
!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE set_mpi_comm_4_solvers(parent_comm, intra_bgrp_comm_, inter_bgrp_comm_ )
  !----------------------------------------------------------------------------
  !
  USE mp_bands_util
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: parent_comm, intra_bgrp_comm_, inter_bgrp_comm_
  ! local variables
  INTEGER :: parent_nproc, parent_mype, ortho_parent_comm_, ierr
  !
  !write(*,*) ' enter set_mpi_comm_4_davidson'
  intra_bgrp_comm   = intra_bgrp_comm_
  inter_bgrp_comm   = inter_bgrp_comm_
  !
!# 61 "set_mpi_comm_4_solvers.f90"
    parent_nproc = 1
    parent_mype  = 0
    nproc_bgrp   = 1
    nbgrp        = 1 
    use_bgrp_in_hpsi = .false.
    my_bgrp_id   = 0
    me_bgrp      = 0
!# 69 "set_mpi_comm_4_solvers.f90"
    !write(*,*) ' exit set_mpi_comm_4_davidson'
    RETURN
  !
END SUBROUTINE set_mpi_comm_4_solvers
