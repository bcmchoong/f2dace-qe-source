!# 1 "error_handler.f90"
!
! Copyright (C) 2002-2025 Quantum ESPRESSO Foundation
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE errore( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output: 
  ! ... if ierr <= 0 it does nothing,
  ! ... if ierr  > 0 it stops.
  !
  ! ...          **** Important note for parallel execution ***
  !
  ! ... in parallel execution unit 6 is written only by the first node;
  ! ... all other nodes have unit 6 redirected to nothing (/dev/null).
  ! ... We write to the "*" unit instead, that appears on all nodes.
  ! ... Effective but annoying!
  !
  USE util_param
!# 28 "error_handler.f90"
  USE mp,        ONLY : mp_abort, mp_rank
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
  INTEGER,          INTENT(IN) :: ierr
    ! the error flag
  INTEGER :: crashunit, mpime
  CHARACTER(LEN=6) :: cerr
  !
  IF( ierr <= 0 ) RETURN
  !
  ! ... the error message is written on the "*" unit
  !
  WRITE( cerr, FMT = '(I6)' ) ierr
  WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
  WRITE( UNIT = *, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
  WRITE( UNIT = *, FMT = '(5X,A)' ) TRIM(message)
  WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
  !
  WRITE( *, '("     stopping ...")' )
  !
  FLUSH( stdout )
  !
!# 64 "error_handler.f90"
  !
  !  .. write the message to a file and close it before exiting
  !  .. this will prevent loss of information on systems that
  !  .. do not flush the open streams
  !  .. added by C.C.
  !
  OPEN( NEWUNIT = crashunit, FILE = crash_file, &
        POSITION = 'APPEND', STATUS = 'UNKNOWN' )
  !
  WRITE( UNIT = crashunit, FMT = '(/,1X,78("%"))' )
!# 78 "error_handler.f90"
  WRITE( UNIT = crashunit, &
         FMT = '(5X,"from ",A," : error #",I10)' ) TRIM(calling_routine), ierr
  WRITE( UNIT = crashunit, FMT = '(5X,A)' ) TRIM(message)
  WRITE( UNIT = crashunit, FMT = '(1X,78("%"),/)' )
  !
  CLOSE( UNIT = crashunit )
  !
!# 91 "error_handler.f90"
  !
  STOP 1
  !
END SUBROUTINE errore
!
!----------------------------------------------------------------------
SUBROUTINE infomsg( routine, message )
  !----------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an info message
  ! ... from a given routine to output.
  !
  USE util_param
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*) :: routine, message
  ! the name of the calling routine
  ! the output message
  !
!  IF ( ionode ) THEN   !if not ionode it is redirected to /dev/null anyway
     !
     WRITE( stdout , '(5X,"Message from routine ",A,":")' ) routine
     WRITE( stdout , '(5X,A)' ) message
     !
!  END IF
  !
  RETURN
  !
END SUBROUTINE infomsg
!
module error_handler
  implicit none
  private
!# 126 "error_handler.f90"
  public :: init_error, add_name, chop_name, error_mem, warning
!# 128 "error_handler.f90"
  type chain
   character (len=35)   :: routine_name
   type(chain), pointer :: previous_link
  end type chain
!# 133 "error_handler.f90"
  type(chain), pointer :: routine_chain
!# 135 "error_handler.f90"
contains
!# 137 "error_handler.f90"
  subroutine init_error(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name
!# 141 "error_handler.f90"
    allocate(routine_chain)
!# 143 "error_handler.f90"
    routine_chain%routine_name  =  routine_name
    nullify(routine_chain%previous_link)
!# 146 "error_handler.f90"
    return
  end subroutine init_error
!# 149 "error_handler.f90"
  subroutine add_name(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name
    type(chain), pointer          :: new_link
!# 154 "error_handler.f90"
    allocate(new_link)
    new_link%routine_name  =  routine_name
    new_link%previous_link => routine_chain
    routine_chain          => new_link
!# 159 "error_handler.f90"
    return
  end subroutine add_name
!# 162 "error_handler.f90"
  subroutine chop_name
    implicit none
    type(chain), pointer :: chopped_chain
!# 166 "error_handler.f90"
    chopped_chain => routine_chain%previous_link
    deallocate(routine_chain)
    routine_chain => chopped_chain
!# 170 "error_handler.f90"
    return
  end subroutine chop_name
!# 173 "error_handler.f90"
  recursive subroutine trace_back(error_code)
!# 175 "error_handler.f90"
    implicit none
    integer :: error_code
!# 178 "error_handler.f90"
    write(unit=*,fmt=*) "   Called by ", routine_chain%routine_name
    if (.not.associated(routine_chain%previous_link)) then
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       write(unit=*,fmt=*) " "
       if( error_code > 0 ) then
          stop
       else
          return
       end if
    end if
!# 190 "error_handler.f90"
    routine_chain => routine_chain%previous_link
    call trace_back(error_code)
!# 193 "error_handler.f90"
  end subroutine trace_back
!# 195 "error_handler.f90"
  subroutine error_mem(message,error_code)
    character (len=*), intent(in) :: message
    integer, intent(in), optional :: error_code
    integer                       :: action_code
    type(chain), pointer          :: save_chain
!# 201 "error_handler.f90"
    if (present(error_code)) then
       action_code = error_code
    else
       action_code = 1
    end if
!# 207 "error_handler.f90"
    if( action_code /= 0 ) then
       write(unit=*,fmt=*) " "
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
!# 212 "error_handler.f90"
       if( action_code > 0 ) then
          write(unit=*,fmt=*) "   Fatal error in routine `", &
               trim(routine_chain%routine_name),"': ",message
       else
          write(unit=*,fmt=*) "   Warning from routine `", &
               trim(routine_chain%routine_name),"': ",message
          save_chain => routine_chain
       end if
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       routine_chain => routine_chain%previous_link
       call trace_back(action_code)
       routine_chain => save_chain
    end if
!# 227 "error_handler.f90"
    return
  end subroutine error_mem
!# 230 "error_handler.f90"
  subroutine warning(message)
    character (len=*), intent(in) :: message
    call error_mem(message,-1)
    return
  end subroutine warning
!# 236 "error_handler.f90"
end module error_handler
