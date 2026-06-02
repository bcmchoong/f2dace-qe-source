!# 1 "autopilot.f90"
! autopilot.f90
!********************************************************************************
! autopilot.f90                                 Copyright (c) 2005 Targacept, Inc.
!********************************************************************************
!   The Autopilot Feature suite is a user level enhancement that enables the 
! following features:  
!      automatic restart of a job; 
!      preconfiguration of job parameters; 
!      on-the-fly changes to job parameters; 
!      and pausing of a running job.  
!
! For more information, see README.AUTOPILOT in document directory.
!
! This program is free software; you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY FOR A PARTICULAR 
! PURPOSE.  See the GNU General Public License at www.gnu.or/copyleft/gpl.txt for 
! more details.
! 
! THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  
! EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES 
! PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESS OR IMPLIED, 
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND THE 
! PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, 
! YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
!
! IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING, 
! WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE 
! THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY 
! GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR 
! INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA 
! BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A 
! FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER 
! OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
!
! You should have received a copy of the GNU General Public License along with 
! this program; if not, write to the 
! Free Software Foundation, Inc., 
! 51 Franklin Street, 
! Fifth Floor, 
! Boston, MA  02110-1301, USA.
! 
! Targacept's address is 
! 200 East First Street, Suite 300
! Winston-Salem, North Carolina USA 27101-4165 
! Attn: Molecular Design. 
! Email: atp@targacept.com
!
! This work was supported by the Advanced Technology Program of the 
! National Institute of Standards and Technology (NIST), Award No. 70NANB3H3065 
!
!********************************************************************************
 
!# 58 "autopilot.f90"
MODULE autopilot
  !---------------------------------------------------------------------------
  !! This module handles the Autopilot Feature Suite.
  !
  !! Written by Lee Atkinson, with help from the ATP team at Targacept, Inc.  
  !! Created June 2005.  
  !! Modified by Yonas Abraham, Sept 2006.
  !
  ! The address for Targacept, Inc. is:
  !   200 East First Street, Suite
  !   300, Winston-Salem, North Carolina 27101; 
  ! Attn: Molecular Design.
  !
  !! See README.AUTOPILOT in the Doc directory for more information.
  !---------------------------------------------------------------------------
!# 74 "autopilot.f90"
  USE kinds
  USE parser, ONLY :  read_line
!# 77 "autopilot.f90"
  IMPLICIT NONE
  SAVE
!# 80 "autopilot.f90"
  INTEGER, parameter :: MAX_INT = huge(1)  
  INTEGER, parameter :: max_event_step = 32  !right now there can be upto 32 Autopilot Events
  INTEGER, parameter :: n_auto_vars = 10     !right now there are only 10 Autopilot Variables
!# 84 "autopilot.f90"
  INTEGER   :: n_events 
  INTEGER   :: event_index = 0
  INTEGER   :: max_rules = 320 !(max_event_step * n_auto_vars)
  INTEGER   :: n_rules 
  INTEGER   :: event_step(max_event_step)
  INTEGER   :: current_nfi
  LOGICAL   :: pilot_p   = .FALSE.    ! pilot property
  LOGICAL   :: restart_p = .FALSE.    ! restart property
  LOGICAL   :: pause_p   = .FALSE.    ! pause   property
  INTEGER   :: pilot_unit = 42   ! perhaps move this to io_files
  CHARACTER(LEN=256) :: pilot_type
!# 96 "autopilot.f90"
  ! AUTOPILOT VARIABLES
  ! These are the variable tables which change the actual variable
  ! dynamically during the course of a simulation. There are many 
  ! parameters which govern a simulation, yet only these are allowed 
  ! to be changed using the event rule mechanism.
  ! Each of these tables are typed according to their variable 
  ! and begin with event_
!# 104 "autopilot.f90"
  !     &CONTROL
  INTEGER   :: rule_isave(max_event_step)
  INTEGER   :: rule_iprint(max_event_step)
  LOGICAL   :: rule_tprint(max_event_step)
  REAL(DP) :: rule_dt(max_event_step)
  !     &SYSTEM
!# 111 "autopilot.f90"
  !     &ELECTRONS
  REAL(DP)         :: rule_emass(max_event_step)
  CHARACTER(LEN=80) :: rule_electron_dynamics(max_event_step)
  REAL(DP)         :: rule_electron_damping(max_event_step)
  CHARACTER(LEN=80) :: rule_electron_orthogonalization(max_event_step)
!# 118 "autopilot.f90"
  !     &IONS
  CHARACTER(LEN=80) :: rule_ion_dynamics(max_event_step)
  REAL(DP)         :: rule_ion_damping(max_event_step)
  CHARACTER(LEN=80) :: rule_ion_temperature(max_event_step)
  REAL(DP) :: rule_tempw(max_event_step)
  INTEGER  :: rule_nhpcl(max_event_step)
  REAL(DP) :: rule_fnosep(max_event_step)
  !     &CELL
!# 127 "autopilot.f90"
  !     &PHONON
!# 130 "autopilot.f90"
  ! Each rule also needs to be correlated (flagged) against the event time table
  ! This is used to flag the a given variable is to be changed on an
  ! event. Initially all set to false, a rule against an event makes it true
  ! Each of these flags are LOGICALs and begin with event_
  !     &CONTROL
  LOGICAL :: event_isave(max_event_step)           
  LOGICAL :: event_iprint(max_event_step)          
  LOGICAL :: event_tprint(max_event_step)
  LOGICAL :: event_dt(max_event_step)              
  !     &SYSTEM
!# 141 "autopilot.f90"
  !     &ELECTRONS
  LOGICAL :: event_emass(max_event_step)   
  LOGICAL :: event_electron_dynamics(max_event_step)   
  LOGICAL :: event_electron_damping(max_event_step)
  LOGICAL :: event_electron_orthogonalization(max_event_step)
!# 147 "autopilot.f90"
  !     &IONS
  LOGICAL :: event_ion_dynamics(max_event_step)   
  LOGICAL :: event_ion_damping(max_event_step)
  LOGICAL :: event_ion_temperature(max_event_step)   
  LOGICAL :: event_tempw(max_event_step)           
  LOGICAL :: event_nhpcl(max_event_step)
  LOGICAL :: event_fnosep(max_event_step)
  !     &CELL
!# 156 "autopilot.f90"
  !     &PHONON
!# 159 "autopilot.f90"
  PRIVATE
  PUBLIC :: auto_check, init_autopilot, card_autopilot, add_rule, & 
       & assign_rule, restart_p, max_event_step, event_index, event_step, rule_isave, &
       & rule_iprint, &
       & rule_tprint, &
       & rule_dt, rule_emass, rule_electron_dynamics, rule_electron_damping, &
       & rule_ion_dynamics, rule_ion_damping, rule_ion_temperature, rule_tempw, &
       & rule_electron_orthogonalization, &
       & event_isave, event_iprint, &
       & event_tprint, &
       & event_dt, event_emass, &
       & event_electron_dynamics, event_electron_damping, event_ion_dynamics, &
       & current_nfi, pilot_p, pilot_unit, pause_p,auto_error, parse_mailbox, &
       & event_ion_damping, event_ion_temperature, event_tempw, &
       & event_electron_orthogonalization, &
       & event_nhpcl, event_fnosep, rule_nhpcl, rule_fnosep
!# 177 "autopilot.f90"
CONTAINS
!# 179 "autopilot.f90"
  !----------------------------------------------------------------------------
  SUBROUTINE auto_error( calling_routine, message)
    !----------------------------------------------------------------------------
    !! This routine calls errore based upon the pilot property flag. 
    !! If the flag is TRUE, the simulation will not stop, 
    !! but the pause property flag is set to TRUE, 
    !! causing pilot to force a state of sleep, 
    !! until the user can mail proper commands. 
    !! Otherwise, its assumed that dynamics have not started 
    !! and this call is an indirect result of read_cards. 
    !! Thus the simulation will stop. 
    !! Either way errore will always issues a warning message.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: calling_routine
    !! the name of the calling calling_routine
    CHARACTER(LEN=*), INTENT(IN) :: message
    !! the output message
    !
    INTEGER :: ierr
    ! the error flag
!# 202 "autopilot.f90"
    IF (pilot_p) THEN
       ! if ierr < 0 errore writes the message but does not stop
       ierr = -1
       pause_p = .TRUE.
       !call mp_bcast(pause_p, ionode_id, world_comm)
    ELSE
       ! if ierr > 0 it stops
       ierr = 1
    ENDIF
!# 212 "autopilot.f90"
    CALL errore( calling_routine, message, ierr )
!# 214 "autopilot.f90"
  END SUBROUTINE auto_error
!# 217 "autopilot.f90"
  !-----------------------------------------------------------------------
  ! AUTO (restart) MODE
  !
  ! Checks if restart files are present, just like what readfile_cp does, 
  ! which calls cp_readfile which create a path to restart.xml. 
  ! This could be checked with inquire, as in check_restartfile. 
  ! If restart_mode=auto, and restart.xml is present, 
  ! then restart_mode="restart", else its "from_scratch".
  ! Set other associated vars, appropriately.
  !
  ! Put this in setcontrol_flags on the select statement
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION auto_check (ndr, outdir)
    !---------------------------------------------------------------------
    !! Checks if restart files are present.
    !
    USE io_global, ONLY: ionode, ionode_id
    USE mp,        ONLY : mp_bcast
    USE mp_world,  ONLY : world_comm
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndr
    !! I/O unit number
    CHARACTER(LEN=*), INTENT(IN) :: outdir
    !! output directory
    !
    ! ... local variables
    !
    CHARACTER(LEN=256) :: dirname, filename
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    LOGICAL :: restart_p = .FALSE.
    INTEGER :: strlen
    ! right now cp_readfile is called with outdir = ' '
    ! so, in keeping with the current design,
    ! the responsibility of setting 
    ! ndr and outdir is the calling program
!# 254 "autopilot.f90"
    IF (ionode) THEN
       dirname = 'RESTART' // int_to_char( ndr )
       IF ( LEN( outdir ) > 1 ) THEN
          strlen  = index(outdir,' ') - 1
          dirname = outdir(1:strlen) // '/' // dirname
       END IF
!# 261 "autopilot.f90"
       filename = TRIM( dirname ) // '/' // 'restart.xml'
       INQUIRE( FILE = TRIM( filename ), EXIST = restart_p )
!# 264 "autopilot.f90"
       auto_check = restart_p
    END IF
    CALL mp_bcast(auto_check, ionode_id, world_comm)
!# 268 "autopilot.f90"
    return
!# 270 "autopilot.f90"
  END FUNCTION auto_check
!# 273 "autopilot.f90"
  !-----------------------------------------------------------------------
  SUBROUTINE init_autopilot
    !---------------------------------------------------------------------
    !! INITIALIZE AUTOPILOT: must be done, even if not in use.
    !! Add this as a call in read_cards.f90 sub: card_default_values.
    !
    IMPLICIT NONE
    !
    integer :: event
!# 283 "autopilot.f90"
    pause_p = .FALSE.
!# 285 "autopilot.f90"
    ! Initialize all events to an iteration that should never occur
    DO event=1,max_event_step
       event_step(event) = MAX_INT
    ENDDO
!# 290 "autopilot.f90"
    n_events = 0 
    n_rules  = 0
    event_index = 1
!# 294 "autopilot.f90"
    ! Initialize here
    !     &CONTROL
    event_isave(:)            = .false.
    event_iprint(:)           = .false.
    event_tprint(:)           = .false.
    event_dt(:)               = .false.
    !     &SYSTEM
    !     &ELECTRONS
    event_emass(:)             = .false.
    event_electron_dynamics(:) = .false.
    event_electron_damping(:)  = .false.
    event_electron_orthogonalization(:) = .false.
!# 308 "autopilot.f90"
    !     &IONS
    event_ion_dynamics(:)      = .false.
    event_ion_damping(:)       = .false.
    event_ion_temperature(:)   = .false.
    event_tempw(:)             = .false.
    !     &CELL
    !     &PHONON
!# 316 "autopilot.f90"
    rule_isave(:)             = 0
    rule_iprint(:)            = 0
    rule_tprint(:)            = .FALSE.
    rule_dt(:)                = 0.0_DP
    rule_emass(:)             = 0.0_DP
    rule_electron_dynamics(:) = 'NONE'
    rule_electron_damping(:)  = 0.0_DP
    rule_ion_dynamics(:)      = 'NONE'
    rule_ion_damping(:)       = 0.0_DP
    rule_ion_temperature(:)   = 'NOT_CONTROLLED'
    rule_tempw(:)             = 0.01_DP
!# 328 "autopilot.f90"
  END SUBROUTINE init_autopilot
!# 332 "autopilot.f90"
  !-----------------------------------------------------------------------  
  SUBROUTINE card_autopilot( input_line )
    !--------------------------------------------------------------------
    !! Called in READ_CARDS and in PARSE_MAILBOX.
    !
    USE io_global, ONLY: ionode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: input_line
    !
    ! ... local variables
    !
    INTEGER :: i, j, linelen
    LOGICAL            :: process_this_line = .FALSE.
    LOGICAL            :: endrules = .FALSE.
    LOGICAL            :: tend = .FALSE.
    LOGICAL, SAVE      :: tread = .FALSE.
    LOGICAL, EXTERNAL  :: matches
    CHARACTER(LEN=1), EXTERNAL :: capital
!# 353 "autopilot.f90"
    !ASU: copied this here since it seems not to be executed during each
    !     call of the routine. Needed for multiple rules in same block
    process_this_line = .FALSE.
    endrules = .FALSE.
    tend = .FALSE.
    tread = .FALSE.
    
    ! This is so we do not read a autopilot card twice from the input file
    IF (( .NOT. pilot_p ) .and. tread ) THEN
       CALL errore( ' card_autopilot  ', ' two occurrences', 2 )
    END IF
!# 365 "autopilot.f90"
    ! If this routined has been called from parse_mailbox
    ! the pilot_type should be set
    IF ( pilot_p ) THEN
       ! IF its MANUAL then process this line first! 
       ! other skip this line and get the next
       IF (TRIM(pilot_type) .eq. 'MANUAL') THEN
          process_this_line = .TRUE.
       ELSE IF (TRIM(pilot_type) .eq. 'PILOT') THEN
          process_this_line = .FALSE.
       ELSE IF (TRIM(pilot_type) .eq. 'AUTO') THEN
          process_this_line = .FALSE.
       ELSE
          IF( ionode ) WRITE(*,*) 'AUTOPILOT: UNRECOGNIZED PILOT TYPE!', TRIM(pilot_type), '===='
          GO TO 10
       END IF
    ELSE
       ! this routine is called from read_cards
       pilot_type = 'AUTO'
       process_this_line = .FALSE.
    END IF
!# 386 "autopilot.f90"
    j=0
    ! must use a local (j) since n_rules may not get updated correctly
    DO WHILE ((.NOT. endrules) .and. j<=max_rules)
       j=j+1
!# 391 "autopilot.f90"
       IF (j > max_rules) THEN
          CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rules May Have Been Execced!')
          go to 10
       END IF
!# 396 "autopilot.f90"
       ! Assume that pilot_p is an indicator and that
       ! this call to card_autopilot is from parse_mailbox
       ! and not from read_cards
       IF(pilot_p) THEN
          IF ( .NOT. process_this_line ) THEN
             ! get the next one
             CALL read_line( input_line, end_of_file = tend)
          END IF
       ELSE
          ! from read_cards
          CALL read_line( input_line, end_of_file = tend)
       END IF
       
       linelen = LEN_TRIM( input_line )
!# 411 "autopilot.f90"
       DO i = 1, linelen
          input_line( i : i ) = capital( input_line( i : i ) )
       END DO
!# 415 "autopilot.f90"
       ! If ENDRULES isnt found, add_rule will fail
       ! and we run out of line anyway
       
       IF ( tend .or. matches( 'ENDRULES', input_line ) ) GO TO 10
!# 420 "autopilot.f90"
       ! Assume this is a rule
       CALL ADD_RULE(input_line)
       ! now, do not process this line anymore
       ! make sure we get the next one
       process_this_line = .FALSE.
!# 426 "autopilot.f90"
    END DO
!# 428 "autopilot.f90"
    IF( ionode ) WRITE(*,*) 'AUTOPILOT SET'
!# 430 "autopilot.f90"
10  CONTINUE
!# 432 "autopilot.f90"
  END SUBROUTINE card_autopilot
!# 437 "autopilot.f90"
  !-----------------------------------------------------------------------
  SUBROUTINE add_rule( input_line )
    !---------------------------------------------------------------------
    !! ADD RULE
    !
    USE io_global, ONLY: ionode
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256) :: input_line
    !
    ! ... local variables
    !
    integer :: i, linelen
    integer :: eq1_pos, eq2_pos, plus_pos, colon_pos
    CHARACTER(LEN=32)  :: var_label
    CHARACTER(LEN=32)  :: value_str
    INTEGER            :: on_step, now_step, plus_step
    integer            :: ios
    integer            :: event
!# 458 "autopilot.f90"
    LOGICAL, EXTERNAL  :: matches
    LOGICAL            :: new_event
!# 462 "autopilot.f90"
    ! this is a temporary local variable
    event = n_events
!# 465 "autopilot.f90"
    ! important for parsing
    i=0
    eq1_pos   = 0
    eq2_pos   = 0
    plus_pos  = 0
    colon_pos = 0
!# 472 "autopilot.f90"
    linelen=LEN_TRIM( input_line )
!# 475 "autopilot.f90"
    ! Attempt to get PLUS SYMBOL
    i = 1
    do while( (plus_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. '+') then
          plus_pos = i
       endif
    end do
!# 484 "autopilot.f90"
    ! Attempt to get a colon
    i = 1
    do while( (colon_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. ':') then
          colon_pos = i
       endif
    end do
!# 493 "autopilot.f90"
    ! Attempt to get first equals
    i = 1
    do while( (eq1_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. '=') then
          eq1_pos = i
       endif
    end do
!# 503 "autopilot.f90"
    ! Attempt to get second equals
    if ((eq1_pos .ne. 0) .and. (eq1_pos < colon_pos)) then
       i = colon_pos + 1
       do while( (eq2_pos .eq. 0) .and. (i <= linelen) )
          i = i + 1
          if(input_line( i : i ) .eq. '=') then
             eq2_pos = i
          endif
       end do
    endif
!# 514 "autopilot.f90"
    ! Complain if there is a bad parse
    if (colon_pos .eq. 0) then
       call auto_error( ' AutoPilot ','Missing colon separator')
       go to 20
    endif
!# 520 "autopilot.f90"
    if (eq1_pos .eq. 0) then
       call auto_error( ' AutoPilot ','Missing equals sign')
       go to 20
    endif
!# 525 "autopilot.f90"
    if ((plus_pos > 0) .and. (eq1_pos < colon_pos)) then
       call auto_error( ' AutoPilot ','equals and plus found prior to colon')
       go to 20
    endif
!# 531 "autopilot.f90"
    !================================================================================
    ! Detect events
    IF ( (pilot_type .eq. 'MANUAL') .or. (pilot_type .eq. 'PILOT') ) THEN
       !-------------------------------------------
       !Do the equivalent of the following:
       !READ(input_line, *) now_label, plus_label1, plus_step, colon_label, var_label, eq_label2, value_str
       !Format:
       ! [NOW [+ plus_step] :] var_label = value_str
       !-------------------------------------------
!# 541 "autopilot.f90"
       ! if there is a NOW, get it and try to get plus and plus_step
!# 543 "autopilot.f90"
       IF ( matches( 'NOW', input_line ) ) THEN
          ! Attempt to get PLUS_STEP
          plus_step = 0
          ! if all is non-trivial, read a PLUS_STEP
          if ((plus_pos > 0) .and. (colon_pos > plus_pos)) then
             ! assume a number lies between
             read(input_line(plus_pos+1:colon_pos-1),*,iostat=ios) plus_step
             if(ios .ne. 0) then
                CALL auto_error( ' AutoPilot ','Value Type Mismatch on NOW line!')
                go to 20
             end if
          end if
          ! set NOW event
          now_step = current_nfi + plus_step
       ELSE
          ! set NOW event
          now_step = current_nfi
       END IF
!# 563 "autopilot.f90"
       !================================================================================
       ! set event
       !
       ! Heres where it get interesting
       ! We may have a new event , or not! :)
!# 569 "autopilot.f90"
       IF ((event-1) .gt. 0) THEN
          IF ( now_step .lt. event_step(event-1)) THEN
             IF( ionode ) write(*,*) ' AutoPilot: current input_line', input_line 
             CALL auto_error( ' AutoPilot ','Dynamic Rule Event Out of Order!')
             go to 20
          ENDIF
       ENDIF
!# 577 "autopilot.f90"
       IF (event .eq. 0) THEN
          new_event = .true.
       ELSEIF ( now_step .gt. event_step(event)) THEN
          new_event = .true.
       ELSE
          new_event = .false.
       ENDIF
!# 585 "autopilot.f90"
       IF ( new_event ) THEN
          ! new event
          event = event + 1
!# 589 "autopilot.f90"
          IF (event > max_event_step) THEN
             IF( ionode ) write(*,*) ' AutoPilot: current input_line', input_line 
             CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rule Event Exceeded!')
             go to 20
          ENDIF
!# 595 "autopilot.f90"
          event_step(event) = now_step
          n_events = event       
       ENDIF
!# 600 "autopilot.f90"
    ELSE IF ( matches( 'ON_STEP', input_line ) ) THEN
       ! Assuming pilot_type is AUTO
       ! if it isnt and ON_STEP these rules wont take anyway
!# 604 "autopilot.f90"
       !-------------------------------------------
       !Do the equivalent of the following:
       !READ(input_line, *) on_step_label, eq_label1, on_step, colon_label, var_label, eq_label2, value_str
       !Format:
       ! ON_STEP = on_step : var_label = value_str
       !-------------------------------------------
!# 611 "autopilot.f90"
       IF( ionode ) write(*,*) 'ADD_RULE: POWER STEERING'
!# 613 "autopilot.f90"
       ! Attempt to get ON_STEP
       on_step = MAX_INT
       ! if all is non-trivial, read a PLUS_STEP
       if ((eq1_pos > 0) .and. (colon_pos > eq1_pos)) then
          ! assume a number lies between
          read(input_line(eq1_pos+1:colon_pos-1),*,iostat=ios) on_step
          if(ios .ne. 0) then
             CALL auto_error( ' AutoPilot ','Value Type Mismatch on ON_STEP line!')
             go to 20
          end if
       end if
       
!# 627 "autopilot.f90"
       !================================================================================
       ! set event
       !
       ! Heres where it get interesting
       ! We may have a new event , or not! :)       
!# 634 "autopilot.f90"
       IF ( ((event-1) .gt. 0)) THEN 
          IF ( on_step .lt. event_step(event-1))  THEN
              IF( ionode ) write(*,*) ' AutoPilot: current input_line', input_line 
              CALL auto_error( ' AutoPilot ','Dynamic Rule Event Out of Order!')
              go to 20
          ENDIF
       ENDIF
!# 643 "autopilot.f90"
       IF (event .eq. 0) THEN
           new_event = .true.
       ELSEIF (on_step .gt. event_step(event)) THEN
           new_event = .true.
       ELSE
           new_event = .false.
       ENDIF
!# 651 "autopilot.f90"
       IF (new_event) THEN
          ! new event
          event = event + 1
          IF (event > max_event_step) THEN
             IF( ionode ) write(*,*) ' AutoPilot: current input_line', input_line
             CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rule Event Exceeded!')
             go to 20
          ENDIF
          event_step(event) = on_step
          n_events = event       
       ENDIF
!# 663 "autopilot.f90"
    END IF ! Event Detection Complete
!# 666 "autopilot.f90"
    !-------------------------------------
    ! Now look for a label and a value
    !-------------------------------------
!# 670 "autopilot.f90"
    if (eq2_pos .eq. 0) then
       var_label = input_line(colon_pos+1: eq1_pos-1)    
       read( input_line(eq1_pos+1:linelen), *, iostat=ios) value_str  
       if(ios .ne. 0) then
          CALL auto_error( ' AutoPilot ','Value Type Mismatch on NOW_STEP line!')
          go to 20
       end if
    else
       var_label = input_line(colon_pos+1: eq2_pos-1)    
       read( input_line(eq2_pos+1:linelen), *, iostat=ios) value_str  
       if(ios .ne. 0) then
          CALL auto_error( ' AutoPilot ','Value Type Mismatch on ON_STEP line!')
          go to 20
       end if
    endif
!# 686 "autopilot.f90"
    ! The Assignment must lie outside the new event scope since
    ! there can exists more than one rule per event
!# 689 "autopilot.f90"
    IF ( (n_rules+1) .gt. max_rules) THEN
       IF( ionode ) write(*,*) ' AutoPilot: current n_rules', n_rules
       CALL auto_error( ' AutoPilot ', ' invalid number of rules ')
       go to 20
    END IF
!# 695 "autopilot.f90"
    call assign_rule(event, var_label, value_str)    
!# 697 "autopilot.f90"
    !IF( ionode ) write(*,*) '  Number of rules: ', n_rules
!# 699 "autopilot.f90"
    FLUSH(6)
!# 701 "autopilot.f90"
20  CONTINUE
!# 703 "autopilot.f90"
  END SUBROUTINE add_rule
!# 706 "autopilot.f90"
  !-----------------------------------------------------------------------
  SUBROUTINE assign_rule(event, var, value)
    !---------------------------------------------------------------------
    !! ASSIGN RULE
    !
    USE io_global, ONLY: ionode
    !
    IMPLICIT NONE
    !
    INTEGER :: event
    CHARACTER(LEN=32) :: var
    CHARACTER(LEN=32) :: value
    !
    ! ... local variables
    !
    INTEGER   :: i, varlen
    INTEGER   :: int_value
    LOGICAL   :: logical_value
    REAL      :: real_value
    REAL(DP) :: realDP_value
    LOGICAL   :: assigned
    LOGICAL, EXTERNAL  :: matches
    CHARACTER(LEN=1), EXTERNAL :: capital
!# 731 "autopilot.f90"
    var = TRIM(var)
    varlen = LEN_TRIM(var)
!# 734 "autopilot.f90"
    DO i = 1, varlen
       var( i : i ) = capital( var( i : i ) )
    END DO
!# 739 "autopilot.f90"
    IF( ionode ) write(*,'("   Reading rule: ",A20,A20)' ) var, value
    assigned = .TRUE.
!# 742 "autopilot.f90"
    IF ( matches( "ISAVE", var ) ) THEN
       read(value, *) int_value
       rule_isave(event)  = int_value
       event_isave(event) = .true.
    ELSEIF ( matches( "IPRINT", var ) ) THEN
       read(value, *) int_value
       rule_iprint(event)  = int_value
       event_iprint(event) = .true.
    ELSEIF ( matches( "TPRINT", var ) ) THEN
       read(value, *) logical_value
       rule_tprint(event)  = logical_value
       event_tprint(event) = .true.
    ELSEIF ( matches( "DT", var ) ) THEN
       read(value, *) real_value
       rule_dt(event)  = real_value
       event_dt(event) = .true.
       !IF( ionode ) write(*,*) 'RULE_DT', rule_dt(event), 'EVENT', event
    ELSEIF ( matches( "EMASS", var ) ) THEN
       read(value, *) realDP_value
       rule_emass(event)  = realDP_value
       event_emass(event) = .true.
    ELSEIF ( matches( "ELECTRON_DYNAMICS", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'SD') .and. (value .ne. 'VERLET') .and. (value .ne. 'DAMP') &
         .and. (value .ne. 'NONE') .and. (value .ne. 'CG') ) then
          call auto_error(' autopilot ',' unknown electron_dynamics '//trim(value) ) 
          assigned = .FALSE.
          go to 40
       endif
       rule_electron_dynamics(event)  = value
       event_electron_dynamics(event) = .true.
    ELSEIF ( matches( "ELECTRON_DAMPING", var ) ) THEN
       read(value, *) realDP_value
       rule_electron_damping(event)  = realDP_value
       event_electron_damping(event) = .true.
    ELSEIF ( matches( "ION_DYNAMICS", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'SD') .and. (value .ne. 'VERLET') .and. (value .ne. 'DAMP') .and. (value .ne. 'NONE')) then
          call auto_error(' autopilot ',' unknown ion_dynamics '//trim(value) )
          assigned = .FALSE.
          go to 40
       endif
       rule_ion_dynamics(event)  = value
       event_ion_dynamics(event) = .true.
    ELSEIF ( matches( "ORTHOGONALIZATION", var) ) THEN
       read(value, *) value
       if ( (value .ne. 'ORTHO') .and. (value .ne. 'GRAM-SCHMIDT') ) then
          call auto_error(' autopilot ',' unknown orthogonalization '//trim(value) )
          assigned = .false.
          go to 40
       endif
       rule_electron_orthogonalization(event) = value
       event_electron_orthogonalization(event) = .true.
    ELSEIF ( matches( "ION_DAMPING", var ) ) THEN
       read(value, *) realDP_value
       rule_ion_damping(event)  = realDP_value
       event_ion_damping(event) = .true.
    ELSEIF ( matches( "ION_TEMPERATURE", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'NOSE') .and. (value .ne. 'NOT_CONTROLLED') .and. (value .ne. 'RESCALING')) then
          call auto_error(' autopilot ',' unknown ion_temperature '//trim(value) )
          assigned = .FALSE.
          go to 40
       endif
       rule_ion_temperature(event)  = value
       event_ion_temperature(event) = .true.
    ELSEIF ( matches( "TEMPW", var ) ) THEN
       read(value, *) realDP_value
       rule_tempw(event)  = realDP_value
       event_tempw(event) = .true.
    ELSEIF ( matches( "NHPCL", var ) ) THEN
       read(value, *) int_value
       rule_nhpcl(event)  = int_value
       event_nhpcl(event) = .true.
    ELSEIF ( matches( "FNOSEP", var ) ) THEN
       read(value, *) realDP_value
       rule_fnosep(event)  = realDP_value
       event_fnosep(event) = .true.
    ELSE
       CALL auto_error( 'autopilot', ' ASSIGN_RULE: FAILED  '//trim(var)//' '//trim(value) )
    END IF
!# 824 "autopilot.f90"
40  if (assigned) then
       n_rules   = n_rules + 1
    else
       IF( ionode ) write(*,*) '  Autopilot: Rule Assignment Failure '
       CALL auto_error( 'autopilot', ' ASSIGN_RULE: FAILED  '//trim(var)//' '//trim(value) )
    endif
!# 831 "autopilot.f90"
  END SUBROUTINE assign_rule
!# 835 "autopilot.f90"
  !-----------------------------------------------------------------------
  SUBROUTINE parse_mailbox ()
    !---------------------------------------------------------------------
    !! Read the mailbox with a mailbox parser:
    !
    !! * if it starts with ON_STEP, then apply to event table etc;
    !! * if not the try to establish that its a variable to set right now.
    !
    USE io_global, ONLY: ionode
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    CHARACTER(LEN=256) :: input_line
    LOGICAL            :: tend
!# 851 "autopilot.f90"
    CHARACTER(LEN=1), EXTERNAL :: capital
    LOGICAL, EXTERNAL  :: matches
!# 855 "autopilot.f90"
    ! we can use this parser routine, since parse_unit=pilot_unit
    CALL read_line( input_line, end_of_file=tend )
    IF (tend) GO TO 50
!# 859 "autopilot.f90"
    DO i = 1, LEN_TRIM( input_line )
       input_line( i : i ) = capital( input_line( i : i ) )
    END DO
!# 863 "autopilot.f90"
    ! This conditional implements the PAUSE feature calling init_auto_pilot, 
    ! will reset this modules global PAUSE_P variable to FALSE
    IF ( matches( "PAUSE", input_line ) .or. &
         matches( "SLEEP", input_line ) .or. &
         matches( "HOVER", input_line ) .or. &
         matches( "WAIT",  input_line ) .or. &
         matches( "HOLD",  input_line ) ) THEN
!# 871 "autopilot.f90"
       IF( ionode ) write(*,*) 'SLEEPING'
       IF( ionode ) write(*,*) 'INPUT_LINE=', input_line
       pause_p = .TRUE.
     ! now you can pass continue to resume 
    ELSE IF (matches( "CONTINUE", input_line ) .or. &
             matches( "RESUME", input_line ) ) THEN
       IF( ionode ) write(*,*) 'RUNNING'
       IF( ionode ) write(*,*) 'INPUT_LINE=', input_line
       pause_p = .FALSE.
!# 881 "autopilot.f90"
       ! Now just quit this subroutine
    ELSE
       ! Also, We didnt see a PAUSE cmd!
       pause_p = .FALSE.
!# 886 "autopilot.f90"
       ! now lets detect the mode for card_autopilot
       ! even though this line will be passed to it the first time
!# 889 "autopilot.f90"
       IF ( matches( "AUTOPILOT", TRIM(input_line) ) ) THEN
          IF( ionode ) WRITE(*,*) '  New autopilot course detected' 
          pilot_type ='AUTO'
       ELSE IF (matches( "PILOT", TRIM(input_line) ) ) THEN
          IF( ionode ) WRITE(*,*) '  Relative pilot course correction detected'
          pilot_type ='PILOT'
       ELSE IF (matches( "NOW", TRIM(input_line) ) ) THEN
          IF( ionode ) WRITE(*,*) '  Manual piloting detected'
          pilot_type ='MANUAL'
       ELSE
          ! Well lets just pause since this guys is throwing trash
          IF( ionode ) WRITE(*,*) '  Mailbox contents not understood: pausing'
          pause_p = .TRUE.
       ENDIF
!# 904 "autopilot.f90"
    END IF
!# 906 "autopilot.f90"
    IF (pause_p) GO TO 50
!# 909 "autopilot.f90"
    ! ok if one adds a rule during steering`
    ! event table must be cleared (from steer point) forward
    !
    ! Every nodes gets this (and the call to card_autopilot
    ! which calls add_rule, which calls assign_rule, etc
    ! In this way we sync the event table
    ! Then we shouldn't have to sync employ_rules variable
    ! changes, or their subroutine side effects (like ions_nose_init)  
!# 918 "autopilot.f90"
    CALL init_autopilot()
!# 920 "autopilot.f90"
    CALL card_autopilot( input_line )
!# 922 "autopilot.f90"
50  CONTINUE
!# 924 "autopilot.f90"
  end subroutine parse_mailbox
!# 927 "autopilot.f90"
END MODULE autopilot
