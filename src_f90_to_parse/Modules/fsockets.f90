!# 1 "fsockets.f90"
!F90 ISO_C_BINGING wrapper for socket communication.
!# 3 "fsockets.f90"
!Copyright (C) 2013, Michele Ceriotti
!# 5 "fsockets.f90"
!Permission is hereby granted, free of charge, to any person obtaining
!a copy of this software and associated documentation files (the
!"Software"), to deal in the Software without restriction, including
!without limitation the rights to use, copy, modify, merge, publish,
!distribute, sublicense, and/or sell copies of the Software, and to
!permit persons to whom the Software is furnished to do so, subject to
!the following conditions:
!# 13 "fsockets.f90"
!The above copyright notice and this permission notice shall be included
!in all copies or substantial portions of the Software.
!# 16 "fsockets.f90"
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
!CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
!TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
!SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!# 26 "fsockets.f90"
   MODULE F90SOCKETS
   !! Contains both the functions that transmit data to the socket and read the data
   !! back out again once finished, and the function which opens the socket initially.
   !! Functions:
   !
   !! * \(\texttt{open_socket}\): Opens a socket with the required host server, socket 
   !!   type and port number;
   !! * \(\texttt{write_buffer}\): Writes a string to the socket;
   !! * \(\texttt{read_buffer}\): Reads data from the socket.
!# 36 "fsockets.f90"
   USE ISO_C_BINDING
   
   IMPLICIT NONE
!# 40 "fsockets.f90"
  INTERFACE writebuffer
      MODULE PROCEDURE writebuffer_s, &
                       writebuffer_d, writebuffer_dv, &
                       writebuffer_i
                       
  END INTERFACE 
!# 47 "fsockets.f90"
  INTERFACE readbuffer
      MODULE PROCEDURE readbuffer_s, &
                       readbuffer_dv, readbuffer_d, &
                       readbuffer_i
                       
  END INTERFACE 
!# 54 "fsockets.f90"
  INTERFACE
    SUBROUTINE open_csocket(psockfd, inet, port, host) BIND(C, name="open_socket")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd, inet, port
    CHARACTER(KIND=C_CHAR), DIMENSION(*)     :: host
!# 60 "fsockets.f90"
    END SUBROUTINE open_csocket
!# 62 "fsockets.f90"
    
    SUBROUTINE writebuffer_csocket(psockfd, pdata, plen) BIND(C, name="writebuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen
!# 69 "fsockets.f90"
    END SUBROUTINE writebuffer_csocket       
!# 71 "fsockets.f90"
    SUBROUTINE readbuffer_csocket(psockfd, pdata, plen) BIND(C, name="readbuffer")
      USE ISO_C_BINDING
    INTEGER(KIND=C_INT)                      :: psockfd
    TYPE(C_PTR), VALUE                       :: pdata
    INTEGER(KIND=C_INT)                      :: plen
!# 77 "fsockets.f90"
    END SUBROUTINE readbuffer_csocket   
  END INTERFACE
!# 80 "fsockets.f90"
   CONTAINS
   
   SUBROUTINE open_socket(psockfd, inet, port, host)      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: inet, port
      INTEGER, INTENT(OUT) :: psockfd
      CHARACTER(LEN=1024), INTENT(IN) :: host
      CHARACTER(LEN=1,KIND=C_CHAR) :: chost(1024)
!# 89 "fsockets.f90"
      CALL fstr2cstr(host, chost)
      CALL open_csocket(psockfd, inet, port, host)
   END SUBROUTINE
!# 93 "fsockets.f90"
   SUBROUTINE fstr2cstr(fstr, cstr, plen)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fstr
      CHARACTER(LEN=1,KIND=C_CHAR), INTENT(OUT) :: cstr(:)
      INTEGER, INTENT(IN), OPTIONAL :: plen
      
      INTEGER i,n
      IF (PRESENT(plen)) THEN
         n = plen
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
      ELSE
         n = LEN_TRIM(fstr)
         DO i=1,n
            cstr(i) = fstr(i:i)
         ENDDO
         cstr(n+1) = C_NULL_CHAR
      END IF
   END SUBROUTINE
!# 114 "fsockets.f90"
  SUBROUTINE writebuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(IN)                :: fdata
!# 119 "fsockets.f90"
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
!# 121 "fsockets.f90"
      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 8)
  END SUBROUTINE
!# 125 "fsockets.f90"
  SUBROUTINE writebuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd, fdata
!# 129 "fsockets.f90"
    INTEGER(KIND=C_INT), TARGET              :: cdata
!# 131 "fsockets.f90"
      cdata = fdata
      CALL writebuffer_csocket(psockfd, c_loc(cdata), 4)
  END SUBROUTINE
!# 135 "fsockets.f90"
  SUBROUTINE writebuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(IN)             :: fstring
    INTEGER, INTENT(IN)                      :: plen
!# 141 "fsockets.f90"
    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)
!# 144 "fsockets.f90"
      DO i = 1,plen
         cstring(i) = fstring(i:i)
      ENDDO
      CALL writebuffer_csocket(psockfd, c_loc(cstring(1)), plen)
  END SUBROUTINE
!# 150 "fsockets.f90"
  SUBROUTINE writebuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(IN), TARGET        :: fdata(plen)
!# 155 "fsockets.f90"
      CALL writebuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE
!# 158 "fsockets.f90"
  SUBROUTINE readbuffer_d (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    REAL(KIND=8), INTENT(OUT)               :: fdata
!# 163 "fsockets.f90"
    REAL(KIND=C_DOUBLE), TARGET              :: cdata
!# 165 "fsockets.f90"
      CALL readbuffer_csocket(psockfd, c_loc(cdata), 8)
      fdata=cdata
  END SUBROUTINE
!# 169 "fsockets.f90"
  SUBROUTINE readbuffer_i (psockfd, fdata)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    INTEGER, INTENT(OUT)                     :: fdata
!# 174 "fsockets.f90"
    INTEGER(KIND=C_INT), TARGET              :: cdata
!# 176 "fsockets.f90"
      CALL readbuffer_csocket(psockfd, c_loc(cdata), 4)
      fdata = cdata
  END SUBROUTINE
!# 180 "fsockets.f90"
  SUBROUTINE readbuffer_s (psockfd, fstring, plen)
      USE ISO_C_BINDING
    INTEGER, INTENT(IN)                      :: psockfd
    CHARACTER(LEN=*), INTENT(OUT)            :: fstring
    INTEGER, INTENT(IN)                      :: plen
!# 186 "fsockets.f90"
    INTEGER                                  :: i
    CHARACTER(LEN=1, KIND=C_CHAR), TARGET    :: cstring(plen)
!# 189 "fsockets.f90"
      CALL readbuffer_csocket(psockfd, c_loc(cstring(1)), plen)
      fstring=""   
      DO i = 1,plen
         fstring(i:i) = cstring(i)
      ENDDO
  END SUBROUTINE
!# 196 "fsockets.f90"
  SUBROUTINE readbuffer_dv(psockfd, fdata, plen)
      USE ISO_C_BINDING  
    INTEGER, INTENT(IN)                      :: psockfd, plen
    REAL(KIND=8), INTENT(OUT), TARGET       :: fdata(plen)
!# 201 "fsockets.f90"
      CALL readbuffer_csocket(psockfd, c_loc(fdata(1)), 8*plen)
  END SUBROUTINE
  END MODULE
