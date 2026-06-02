!# 1 "mp_bands.f90"
!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_bands
  !----------------------------------------------------------------------------
  !! Band groups (processors within a pool of bands).  
  !! Subdivision of pool group, used for parallelization over bands.
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  INTEGER :: nbgrp = 1
  !! number of band groups
  INTEGER :: nproc_bgrp = 1
  !! number of processors within a band group
  INTEGER :: me_bgrp = 0
  !! index of the processor within a band group
  INTEGER :: root_bgrp = 0
  !! index of the root processor within a band group
  INTEGER :: my_bgrp_id = 0
  !! index of my band group
  INTEGER :: root_bgrp_id = 0
  !! index of root band group
  INTEGER :: inter_bgrp_comm = 0
  !! inter band group communicator
  INTEGER :: intra_bgrp_comm = 0
  !! intra band group communicator  
  !
  LOGICAL :: use_bgrp_in_hpsi = .FALSE.
  !! TRUE if band parallelization is performed inside \( H\psi \)
  !! and \( S\psi \), FALSE otherwise (band parallelization can be 
  !! performed outside \( H\psi \) and \( S\psi \) though).
  !
  ! ... "task" groups (for band parallelization of FFT)
  !
  INTEGER :: ntask_groups = 1
  !! number of proc. in an orbital "task group"
  !
  ! ... "nyfft" groups (to push FFT parallelization beyond the nz-planes limit)
  !
  INTEGER :: nyfft = 1
  !! number of y-fft groups. By default =1, i.e. y-ffts are done by a 
  !! single proc.
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_bands( nband_, ntg_, nyfft_, parent_comm )
    !--------------------------------------------------------------------------
    !! Divide processors (of the "parent_comm" group) into nband_ pools.
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nband_
    !! it must have been previously read from command line argument
    !! by a call to routine \(\texttt{get\_command\_line}\).
    INTEGER, INTENT(IN) :: parent_comm
    !! typically processors of a k-point pool (intra\_pool\_comm)
    INTEGER, INTENT(IN), OPTIONAL :: ntg_, nyfft_
    !
    INTEGER :: parent_nproc = 1, parent_mype = 0
    !
!# 121 "mp_bands.f90"
    RETURN
    !
  END SUBROUTINE mp_start_bands
  !
END MODULE mp_bands
!
!     
MODULE mp_bands_TDDFPT
  !! Starting and ending band indexes in bgrp parallelization.
  ! NB: These two variables used to be in mp_bands and are loaded from mp_global in TDDFPT 
  !     I think they would better stay in a TDDFPT specific module but leave them here not to
  !     be too invasive on a code I don't know well. SdG
  !     
  INTEGER :: ibnd_start = 0
  !! starting band index used in bgrp parallelization.
  INTEGER :: ibnd_end = 0
  !! ending band index used in bgrp parallelization.
!     
END MODULE mp_bands_TDDFPT
!     
