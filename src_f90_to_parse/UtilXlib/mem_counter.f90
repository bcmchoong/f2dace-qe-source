!# 1 "mem_counter.f90"
!
SUBROUTINE mem_counter(valore, sumsign, label, sub)
    ! Poor-man memory counter - to be called after each allocation/deallocation
    ! if sumsign = +/-1, add/remove "valore" from allocated memory counter
    !    "label" should be the name of the allocated/deallocated variable
    ! if sumsign = 0, print the current value of the allocated memory counter
    !    "label" should identify the place of the code where you call from
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: valore
    INTEGER, INTENT(IN) :: sumsign
    CHARACTER (LEN=*), INTENT(IN) :: label
    CHARACTER (LEN=*), INTENT(IN) :: sub
    !
    INTEGER, SAVE :: current_memory = 0
    INTEGER, SAVE :: last_print_val = 0
    INTEGER, PARAMETER :: MB = 1024*1024
    !
    IF ( sumsign /= 0 .AND. ABS(sumsign) /= 1 ) RETURN
!# 29 "mem_counter.f90"
    IF ( sumsign == 0 ) THEN
       WRITE(*, '("Memory allocated in ",a,":", I6, " MB")') &
                trim(label), current_memory/MB
    ELSE
       current_memory = current_memory + valore*sumsign
       ! normal output: print when memory has grown by more than 2 MB
       IF ( ( current_memory - last_print_val ) > 2*MB) THEN
          WRITE(*, '("Max allocated memory: ", I6, " MB, variable: ",a)') &
              current_memory/MB, trim(label)
          last_print_val = current_memory
       END IF
    END IF
    !
END SUBROUTINE
