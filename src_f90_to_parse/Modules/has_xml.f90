!# 1 "has_xml.f90"
!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION has_xml(inp_string)
!
!! This function returns TRUE if the last four characters of inp_string
!! are .xml or .XML. On output the string .xml or .XML is removed from 
!! \(\text{inp_string}\).
!
IMPLICIT NONE
!# 16 "has_xml.f90"
CHARACTER(LEN=*), INTENT(INOUT) :: inp_string
!# 18 "has_xml.f90"
INTEGER :: leng, start 
CHARACTER(LEN=4) :: aux
LOGICAL, EXTERNAL :: matches
!# 22 "has_xml.f90"
has_xml=.FALSE.
leng=LEN_TRIM(inp_string)
!# 25 "has_xml.f90"
!cannot match xml if it is only 1 or 2 chars long
IF(leng<3) RETURN
!# 28 "has_xml.f90"
start=MAX(leng-3,1)
aux=inp_string(start:leng)
!# 31 "has_xml.f90"
IF (matches(aux,'.xml').OR.matches(aux,'.XML')) THEN
   has_xml=.TRUE.
   inp_string(leng-3:leng)=' '
ENDIF
!# 36 "has_xml.f90"
RETURN
END FUNCTION has_xml
