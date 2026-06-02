!# 1 "mp_swap_ax_rism.f90"
!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE mp_swap_ax_rism(mp_site, mp_task, lsite, rsite, ltask, rtask, isign)
  !---------------------------------------------------------------------------
  !
  ! ... swap parallel axies of RISM's data (for real numbers)
  ! ... if isign > 0, swap rtask -> rsite.
  ! ... if isign < 0, swap rsite -> rtask.
  !
  USE kinds,   ONLY : DP
  USE mp,      ONLY : mp_max
  USE mp_rism, ONLY : mp_rism_site, mp_rism_task
  USE parallel_include
  !
  IMPLICIT NONE
  !
  TYPE(mp_rism_site), INTENT(IN)    :: mp_site
  TYPE(mp_rism_task), INTENT(IN)    :: mp_task
  INTEGER,            INTENT(IN)    :: lsite
  REAL(DP),           INTENT(INOUT) :: rsite(lsite,1:*)
  INTEGER,            INTENT(IN)    :: ltask
  REAL(DP),           INTENT(INOUT) :: rtask(ltask,1:*)
  INTEGER,            INTENT(IN)    :: isign
  !
  REAL(DP), ALLOCATABLE :: rwork(:)
  !
  ALLOCATE(rwork(mp_task%nvec))
  !
  IF (isign > 0) THEN
    CALL rtask_to_rsite()
  ELSE IF (isign < 0) THEN
    CALL rsite_to_rtask()
  ELSE !IF (isign == 0) THEN
    ! NOP
  END IF
  !
  DEALLOCATE(rwork)
  !
CONTAINS
  !
!# 126 "mp_swap_ax_rism.f90"
  SUBROUTINE rtask_to_rsite()
    IMPLICIT NONE
    !
    rsite(1:mp_task%nvec, 1:mp_site%nsite) = rtask(1:mp_task%nvec, 1:mp_site%nsite)
  END SUBROUTINE rtask_to_rsite
  !
  SUBROUTINE rsite_to_rtask()
    IMPLICIT NONE
    !
    rtask(1:mp_task%nvec, 1:mp_site%nsite) = rsite(1:mp_task%nvec, 1:mp_site%nsite)
  END SUBROUTINE rsite_to_rtask
!# 138 "mp_swap_ax_rism.f90"
  !
END SUBROUTINE mp_swap_ax_rism
