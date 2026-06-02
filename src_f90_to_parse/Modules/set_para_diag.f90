!# 1 "set_para_diag.f90"
!
! Copyright (C) 2021 Quantum ESPRESSO Fondation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE set_para_diag( nbnd, use_para_diag )
  !-----------------------------------------------------------------------------
  !! Sets up the communicator used for parallel diagonalization in LAXlib.
  !! Merges previous code executed at startup with function "check_para_diag".
  !! To be called after the initialization of variables is completed and
  !! the dimension of matrices to be diagonalized is known
  !
  USE io_global,            ONLY : stdout, ionode
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp_world,             ONLY : world_comm
  USE mp_exx,               ONLY : negrp
  USE command_line_options, ONLY : ndiag_
!# 23 "set_para_diag.f90"
  IMPLICIT NONE
!# 25 "set_para_diag.f90"
  INCLUDE 'laxlib.fh'
!# 27 "set_para_diag.f90"
  INTEGER, INTENT(IN) :: nbnd
  !! dimension of matrices to be diagonalized (number of bands)
  LOGICAL, INTENT(INOUT) :: use_para_diag
  !! true if parallel linear algebra is to be used
  !
  LOGICAL, SAVE :: first = .TRUE.
  LOGICAL :: do_diag_in_band_group = .TRUE.
  INTEGER :: np_ortho(2), ortho_parent_comm
!# 36 "set_para_diag.f90"
  IF( .NOT. first ) RETURN
  first = .FALSE.
  !
  IF( negrp > 1 .OR. do_diag_in_band_group ) THEN
     ! one diag group per bgrp with strict hierarchy: POOL > BAND > DIAG
     ! if using exx groups from mp_exx,  always use this diag method
     CALL laxlib_start ( ndiag_, intra_bgrp_comm, .TRUE. )
  ELSE
     ! one diag group per pool ( individual k-point level )
     ! with band group and diag group both being children of POOL comm
     CALL laxlib_start ( ndiag_, intra_pool_comm, .FALSE. )
  END IF
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
       inter_bgrp_comm )
  !
!# 91 "set_para_diag.f90"
  !
  use_para_diag = .FALSE.
  !
!# 95 "set_para_diag.f90"
  !
  RETURN
END SUBROUTINE set_para_diag
