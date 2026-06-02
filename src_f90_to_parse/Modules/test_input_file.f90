!# 1 "test_input_file.f90"
!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION test_input_xml (myunit)
   !
   !! Check if file opened as unit "myunit" is a xml file or not.
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(in) :: myunit
   !
   CHARACTER(LEN=256) :: dummy
   CHARACTER(LEN=1), EXTERNAL :: capital
   INTEGER :: i, j 
   LOGICAL :: exst
   !
   test_input_xml = .false.
   INQUIRE ( UNIT=myunit, EXIST=exst )
   IF ( .NOT. exst ) GO TO 10
   
   ! read until a non-empty line is found
!# 27 "test_input_file.f90"
   dummy = ' '
   DO WHILE ( LEN_TRIM(dummy) < 1 )
      READ ( myunit,'(A)', ERR=10, END=10) dummy
   END DO
!# 32 "test_input_file.f90"
   ! remove blanks from line, convert to capital, clean trailing characters
!# 34 "test_input_file.f90"
   j=1
   DO i=1, LEN_TRIM(dummy) 
      IF ( dummy(i:i) /= ' ' ) THEN
         dummy(j:j) = capital(dummy(i:i))
         j=j+1
      END IF
   END DO
   DO i=j, LEN_TRIM(dummy) 
      dummy(i:i) = ' '
   END DO
!# 45 "test_input_file.f90"
   ! check for string "<?xml" or "<xml" in the beginning, ">" at the end
!# 47 "test_input_file.f90"
   j = LEN_TRIM (dummy)
   test_input_xml = ( (dummy(1:5) == "<?XML") .OR. (dummy(1:4) == "<XML") ) &
                      .AND. (dummy(j:j) == ">")
   RETURN
!# 52 "test_input_file.f90"
10 WRITE (0,"('from test_input_xml: input file not opened or empty')")
!# 54 "test_input_file.f90"
END FUNCTION test_input_xml
