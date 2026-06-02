!# 1 "add_dmft_occ.f90"
!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! --------------------------------------------------------------
MODULE add_dmft_occ
  ! Written by Sophie D. Beck (2020-2021)
  ! arXiv:2111.10289 (2021)
  !
  ! This module contains a subroutine to read DMFT occupation uddaptes
  ! and to change the fermi weights and the wavefunctions accordingly.
  !
  !
  USE wvfct,                ONLY : wg, nbnd, et
  USE klist,                ONLY : wk, nkstot, nks
  USE kinds,                ONLY : dp
  USE io_global,			ONLY : ionode, ionode_id, stdout
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_images,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE control_flags,        ONLY : restart
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: dmft
  !! if .TRUE.: update occupations/eigenvectors once
  LOGICAL :: dmft_updated = .FALSE.
  !! if .TRUE.: the update was done
  CHARACTER(len=256) :: dmft_prefix
  !! name of hdf5 input archive containing the occupation update and bandwindow
  COMPLEX(DP), ALLOCATABLE :: v_dmft(:,:,:)
  !! transformation matrix from previous to new eigenvectors
  !
  PUBLIC :: dmft, dmft_updated, v_dmft
  PUBLIC :: dmft_update
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE dmft_update()
!# 216 "add_dmft_occ.f90"
  END SUBROUTINE dmft_update
  !
  !------------------------------------------------------------------------
END MODULE add_dmft_occ
