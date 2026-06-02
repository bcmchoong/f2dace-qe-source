!# 1 "stop_run.f90"
!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run( exit_status )
  !----------------------------------------------------------------------------
  !! Close all files and synchronize processes before stopping. 
  !! Remove temporary files needed for restart only if exit_status = 0
  !! (successful execution)
  !
  USE io_global,          ONLY : ionode
  USE mp_global,          ONLY : mp_global_end
  USE environment,        ONLY : environment_end
  USE io_files,           ONLY : iuntmp, seqopn
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: exit_status
  LOGICAL             :: exst, opnd, lflag
  !
  lflag = ( exit_status == 0 ) 
  IF ( lflag ) THEN
     ! 
     ! ... remove files needed only to restart
     !
     CALL seqopn( iuntmp, 'restart', 'UNFORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     !
     IF ( ionode ) THEN
        CALL seqopn( iuntmp, 'update', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
        CALL seqopn( iuntmp, 'para', 'FORMATTED', exst )
        CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     ENDIF
     !
  ENDIF
  !
  CALL close_files( lflag )
  !
  CALL print_clock_pw()
  !
  CALL clean_pw( .TRUE. )
  !
  CALL environment_end( )
  !
  CALL mp_global_end()
  !
END SUBROUTINE stop_run
!
!-----------------------------------------
SUBROUTINE do_stop( exit_status )
  !---------------------------------------
  !! Stop the run. Exit status is returned to the shell only if
  !! preprocessing flag __RETURN_EXIT_STATUS is set (default: no).
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: exit_status
  !
!# 65 "stop_run.f90"
  STOP
!# 94 "stop_run.f90"
  !
END SUBROUTINE do_stop
!
!----------------------------------------------------------------------------
SUBROUTINE closefile()
  !----------------------------------------------------------------------------
  !! Close all files and synchronize processes before stopping.  
  !! Called by "sigcatch" when it receives a signal.
  !
  USE io_global,  ONLY :  stdout
  !
  WRITE( stdout,'(5X,"Signal Received, stopping ... ")')
  !
  CALL stop_run( 255 )
  !
  RETURN
  !
END SUBROUTINE closefile
