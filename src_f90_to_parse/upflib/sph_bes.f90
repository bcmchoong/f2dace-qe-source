!# 1 "sph_bes.f90"
!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine sph_bes (msh, r, q, l, jl)
  !$acc routine vector
  !--------------------------------------------------------------------
  !! Spherical Bessel function.
  !
  USE upf_kinds, only: DP
  USE upf_const, only: eps14
  !
  implicit none
  !
  integer :: msh
  !! number of grid points points
  integer :: l
  !! angular momentum (-1 <= l <= 6)
  real(DP) :: r (msh)
  !! radial grid
  real(DP) :: q
  !! q
  real(DP) :: jl (msh)
  !! Output: Spherical Bessel function \(j_l(q*r(i))\)
  !
  ! xseries = convergence radius of the series for small x of j_l(x)
  real(DP) :: x, xl, xseries = 0.05_dp
  integer :: i, ir, ir0
  integer :: semifact
  !
!# 39 "sph_bes.f90"
 
  !  case q=0
!# 42 "sph_bes.f90"
  if (abs (q) < eps14) then
     if (l == -1) then
        stop !call upf_error ('sph_bes', 'j_{-1}(0) ?!?', 1)
     elseif (l == 0) then
        !$acc loop vector
        do ir = 1, msh
          jl(ir) = 1.d0
        enddo
     else
        !$acc loop vector
        do ir = 1, msh
          jl(ir) = 0.d0
        enddo  
     endif
     return
  end if 
!# 59 "sph_bes.f90"
  !  case l=-1
!# 61 "sph_bes.f90"
  if (l == - 1) then
     if (abs (q * r (1) ) < eps14) stop !call upf_error ('sph_bes', 'j_{-1}(0) ?!?',1)
!# 71 "sph_bes.f90"
     !$acc loop vector
     do ir = 1, msh
       jl (ir) = cos (q * r (ir) ) / (q * r (ir) )
     enddo
!# 77 "sph_bes.f90"
     return
!# 79 "sph_bes.f90"
  end if
!# 81 "sph_bes.f90"
  ! series expansion for small values of the argument
  ! ir0 is the first grid point for which q*r(ir0) > xseries
  ! notice that for small q it may happen that q*r(msh) < xseries !
!# 85 "sph_bes.f90"
  ir0 = msh+1
  !$acc loop vector
  do ir = 1, msh
     if ( abs (q * r (ir) ) > xseries ) then
        ir0 = ir
        exit
     end if
  end do
!# 94 "sph_bes.f90"
  !$acc loop vector
  do ir = 1, ir0 - 1
     x = q * r (ir)
     if ( l == 0 ) then
        xl = 1.0_dp
     else
        xl = x**l
     end if
     !--
     semifact = 1
     !$acc loop seq reduction(*:semifact)
     do i = 2*l+1, 1, -2
       semifact = i*semifact
     enddo
     !---
     jl (ir) = xl/DBLE(semifact) * &
                ( 1.0_dp - x**2/1.0_dp/2.0_dp/DBLE(2*l+3) * &
                ( 1.0_dp - x**2/2.0_dp/2.0_dp/DBLE(2*l+5) * &
                ( 1.0_dp - x**2/3.0_dp/2.0_dp/DBLE(2*l+7) * &
                ( 1.0_dp - x**2/4.0_dp/2.0_dp/DBLE(2*l+9) ) ) ) )
  end do
!# 116 "sph_bes.f90"
  ! the following shouldn't be needed but do you trust compilers
  ! to do the right thing in this special case ? I don't - PG
!# 119 "sph_bes.f90"
  if ( ir0 > msh ) return
!# 121 "sph_bes.f90"
  if (l == 0) then
!# 130 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = sin (q * r (ir) ) / (q * r (ir) )
     enddo
!# 136 "sph_bes.f90"
  elseif (l == 1) then
!# 147 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q * r (ir) ) / (q * r (ir) ) - &
                  cos (q * r (ir) ) ) / (q * r (ir) )
     enddo
!# 154 "sph_bes.f90"
  elseif (l == 2) then
!# 165 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = ( (3.d0 / (q*r(ir)) - (q*r(ir)) ) * sin (q*r(ir)) - &
                    3.d0 * cos (q*r(ir)) ) / (q*r(ir))**2
     enddo
!# 172 "sph_bes.f90"
  elseif (l == 3) then
!# 185 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q*r(ir)) * &
                 (15.d0 / (q*r(ir)) - 6.d0 * (q*r(ir)) ) + &
                 cos (q*r(ir)) * ( (q*r(ir))**2 - 15.d0) ) / &
                 (q*r(ir)) **3
     enddo
!# 194 "sph_bes.f90"
  elseif (l == 4) then
!# 208 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (sin (q*r(ir)) * &
                 (105.d0 - 45.d0 * (q*r(ir))**2 + (q*r(ir))**4) + &
                 cos (q*r(ir)) * &
                 (10.d0 * (q*r(ir))**3 - 105.d0 * (q*r(ir))) ) / &
                    (q*r(ir))**5
     enddo
!# 218 "sph_bes.f90"
  elseif (l == 5) then
!# 231 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = (-cos(q*r(ir)) - &
                  (945.d0*cos(q*r(ir))) / (q*r(ir)) ** 4 + &
                  (105.d0*cos(q*r(ir))) / (q*r(ir)) ** 2 + &
                  (945.d0*sin(q*r(ir))) / (q*r(ir)) ** 5 - &
                  (420.d0*sin(q*r(ir))) / (q*r(ir)) ** 3 + &
                  ( 15.d0*sin(q*r(ir))) / (q*r(ir)) ) / (q*r(ir))
     enddo
!# 242 "sph_bes.f90"
  elseif (l == 6) then
!# 257 "sph_bes.f90"
     !$acc loop vector
     do ir = ir0, msh
       jl (ir) = ((-10395.d0*cos(q*r(ir))) / (q*r(ir))**5 + &
                  (  1260.d0*cos(q*r(ir))) / (q*r(ir))**3 - &
                  (    21.d0*cos(q*r(ir))) / (q*r(ir))    - &
                             sin(q*r(ir))                   + &
                  ( 10395.d0*sin(q*r(ir))) / (q*r(ir))**6 - &
                  (  4725.d0*sin(q*r(ir))) / (q*r(ir))**4 + &
                  (   210.d0*sin(q*r(ir))) / (q*r(ir))**2 ) / (q*r(ir))
     enddo
!# 269 "sph_bes.f90"
  else
!# 271 "sph_bes.f90"
     stop !call upf_error ('sph_bes', 'not implemented', abs(l))
!# 273 "sph_bes.f90"
  endif
  !
  return
end subroutine sph_bes
!# 278 "sph_bes.f90"
integer function semifact(n)
  ! semifact(n) = n!!
  implicit none
  integer :: n, i
!# 283 "sph_bes.f90"
  semifact = 1
  do i = n, 1, -2
     semifact = i*semifact
  end do
  return
end function semifact
!
SUBROUTINE sph_dbes ( nr, r, xg, l, jl, djl )
  !
  !! Calculates \(x*dj_l(x)/dx\) using the recursion formula:
  !! $$ dj_l(x)/dx = l/x j_l(x) - j_{l+1}(x) $$
  !! for \(l=0\), and for \(l>0\):
  !! $$ dj_l(x)/dx = j_{l-1}(x) - (l+1)/x j_l(x) $$
  !! Requires \(j_l(r)\) in input.
  !! Used only in CP. Note that upflib uses numerical differentiation.
  !
  USE upf_kinds, only: DP
  USE upf_const, ONLY : eps8
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l, nr
  REAL (DP), INTENT(IN) :: xg, jl(nr), r(nr)
  REAL (DP), INTENT(OUT):: djl(nr)
  !
  if ( xg < eps8 ) then
     !
     ! special case q=0
     ! note that x*dj_l(x)/dx = 0 for x = 0
     !
     djl(:) = 0.0d0
  else
     !
     if ( l > 0 ) then
        call sph_bes ( nr, r, xg, l-1, djl )
        djl(:) = djl(:) * (xg * r(:) ) - (l+1) * jl(:)
     else if ( l == 0 ) then
        call sph_bes ( nr, r, xg, l+1, djl )
        djl(:) = - djl(:) * (xg * r(:) )
     else
        call upf_error('sph_dbes','l < 0 not implemented', abs(l) )
     end if
  end if
  !
end SUBROUTINE sph_dbes
