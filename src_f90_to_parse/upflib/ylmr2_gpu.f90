!# 1 "ylmr2_gpu.f90"
!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module ylmr2_gpum
!# 113 "ylmr2_gpu.f90"
end module ylmr2_gpum
!# 115 "ylmr2_gpu.f90"
subroutine ylmr2_gpu(lmax2, ng, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax, GPU version
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !     Last modified Jan 2024, by PG
  !
!# 129 "ylmr2_gpu.f90"
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  real(DP), intent(out) :: ylm (ng,lmax2)
!# 161 "ylmr2_gpu.f90"
  call upf_error('ylmr2_gpu','you should not be here, go away!',1)
!# 163 "ylmr2_gpu.f90"
  return
end subroutine ylmr2_gpu
