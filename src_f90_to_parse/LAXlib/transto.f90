!# 1 "transto.f90"
!
! Copyright (C) 2001 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 10 "transto.f90"
!     OPTIMIZED DRIVER FOR MATRIX TRASPOSITION
!   
!     written by Carlo Cavazzoni 
!
!# 20 "transto.f90"
      SUBROUTINE mytranspose(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!
        IMPLICIT NONE
        include 'laxlib_kinds.fh'
!# 28 "transto.f90"
        INTEGER :: ldx, ldy, n, m, what
        REAL(DP) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = 35
        REAL(DP) :: buf(bsiz, bsiz), bswp
!# 35 "transto.f90"
        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n: ",2I6)') ldx, n
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m: ",2I6)') ldy, m
        end if
!# 42 "transto.f90"
        nb = n / bsiz 
        mb = m / bsiz 
!# 45 "transto.f90"
        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF
!# 51 "transto.f90"
        select case (what)
!# 53 "transto.f90"
          case (1)
!# 55 "transto.f90"
            do i=1,n
              do j=1,m
                y(j,i) = x(i,j)
              enddo
            enddo
!# 61 "transto.f90"
          case (2)
!# 63 "transto.f90"
            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = buf(j,i)
                  enddo
                enddo
              enddo
            enddo
!# 87 "transto.f90"
            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 104 "transto.f90"
            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 121 "transto.f90"
            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = buf(i,j)
                enddo
              enddo
            END IF
!# 136 "transto.f90"
          case default
!# 138 "transto.f90"
            write(6,fmt='("trasponi: undefined method")')
!# 140 "transto.f90"
        end select
!# 142 "transto.f90"
        RETURN
      END SUBROUTINE  mytranspose
!# 147 "transto.f90"
      SUBROUTINE mytransposez(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!
!# 153 "transto.f90"
        IMPLICIT NONE
        include 'laxlib_kinds.fh'
!# 157 "transto.f90"
        INTEGER :: ldx, ldy, n, m, what
        COMPLEX(DP) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = 35 / 2
        COMPLEX(DP) :: buf(bsiz, bsiz), bswp
!# 164 "transto.f90"
        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n")')
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m")')
        end if
!# 171 "transto.f90"
        nb = n / bsiz 
        mb = m / bsiz 
!# 174 "transto.f90"
        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF
!# 180 "transto.f90"
        select case (what)
!# 182 "transto.f90"
          case (1)
!# 184 "transto.f90"
            do i=1,n
              do j=1,m
                y(j,i) = x(i,j)
              enddo
            enddo
!# 190 "transto.f90"
          case (2)
!# 192 "transto.f90"
            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = buf(j,i)
                  enddo
                enddo
              enddo
            enddo
!# 216 "transto.f90"
            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 233 "transto.f90"
            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 250 "transto.f90"
            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = buf(i,j)
                enddo
              enddo
            END IF
!# 265 "transto.f90"
          case default
!# 267 "transto.f90"
            write(6,fmt='("trasponi: undefined method")')
!# 269 "transto.f90"
        end select
!# 271 "transto.f90"
        RETURN
      END SUBROUTINE mytransposez
!# 277 "transto.f90"
      SUBROUTINE mytranspose_sp(x, ldx, y, ldy, n, m)
!
!     x  input matrix (n by m) to be trasposed
!     y  output matrix (m by n), the transpose of x
!
        IMPLICIT NONE
        include 'laxlib_kinds.fh'
!# 285 "transto.f90"
        INTEGER :: ldx, ldy, n, m, what
        REAL(SP) :: x(ldx, m), y(ldy, n)
        INTEGER :: i, j, k, d, nb, mb, ib, jb, ioff, joff
        INTEGER :: iind, jind
        INTEGER,  PARAMETER :: bsiz = 35
        REAL(SP) :: buf(bsiz, bsiz), bswp
!# 292 "transto.f90"
        if( n>ldx ) then
          write(6,fmt='("trasponi: inconsistent ldx and n: ",2I6)') ldx, n
        end if
        if( m>ldy ) then
          write(6,fmt='("trasponi: inconsistent ldy and m: ",2I6)') ldy, m
        end if
!# 299 "transto.f90"
        nb = n / bsiz 
        mb = m / bsiz 
!# 302 "transto.f90"
        IF( nb < 2 .AND. mb < 2 ) THEN
          what = 1
        ELSE
          what = 2
        END IF
!# 308 "transto.f90"
        select case (what)
!# 310 "transto.f90"
          case (1)
!# 312 "transto.f90"
            do i=1,n
              do j=1,m
                y(j,i) = x(i,j)
              enddo
            enddo
!# 318 "transto.f90"
          case (2)
!# 320 "transto.f90"
            do ib = 1, nb
              ioff = (ib-1) * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, bsiz
                    buf(i,j) = x(i+ioff, j+joff)
                  enddo
                enddo
                do j = 1, bsiz
                  do i = 1, j-1
                    bswp = buf(i,j)
                    buf(i,j) = buf(j,i) 
                    buf(j,i) = bswp
                  enddo
                enddo
                do i=1,bsiz
                  do j=1,bsiz
                    y(j+joff, i+ioff) = buf(j,i)
                  enddo
                enddo
              enddo
            enddo
!# 344 "transto.f90"
            IF( MIN(1, MOD(n, bsiz)) > 0 ) THEN
              ioff = nb * bsiz
              do jb = 1, mb
                joff = (jb-1) * bsiz
                do j = 1, bsiz
                  do i = 1, MIN(bsiz, n-ioff)
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, MIN(bsiz, n-ioff)
                  do j = 1, bsiz
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 361 "transto.f90"
            IF( MIN(1, MOD(m, bsiz)) > 0 ) THEN
              joff = mb * bsiz
              do ib = 1, nb
                ioff = (ib-1) * bsiz
                do j = 1, MIN(bsiz, m-joff)
                  do i = 1, bsiz
                    buf(i,j) =  x(i+ioff, j+joff)
                  enddo
                enddo
                do i = 1, bsiz
                  do j = 1, MIN(bsiz, m-joff)
                    y(j+joff,i+ioff) = buf(i,j)
                  enddo
                enddo
              enddo
            END IF
!# 378 "transto.f90"
            IF( MIN(1,MOD(n,bsiz))>0 .AND. MIN(1,MOD(m,bsiz))>0 ) THEN
              joff = mb * bsiz
              ioff = nb * bsiz
              do j = 1, MIN(bsiz, m-joff)
                do i = 1, MIN(bsiz, n-ioff)
                  buf(i,j) =  x(i+ioff, j+joff)
                enddo
              enddo
              do i = 1, MIN(bsiz, n-ioff)
                do j = 1, MIN(bsiz, m-joff)
                  y(j+joff,i+ioff) = buf(i,j)
                enddo
              enddo
            END IF
!# 393 "transto.f90"
          case default
!# 395 "transto.f90"
            write(6,fmt='("trasponi: undefined method")')
!# 397 "transto.f90"
        end select
!# 399 "transto.f90"
        RETURN
      END SUBROUTINE  mytranspose_sp
