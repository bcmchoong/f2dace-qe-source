!# 1 "beef.f90"
!
MODULE beef
  !
  ! BEEF-vdW Module
  !
  ! by Gabriel S. Gusmao :: gusmaogabriels@gmail.com
  ! Medford Groups @ ChBE, Georgia Institute of Technology
  !
  ! adapted from 
  ! 
  ! Johannes Voss
  !  - https://github.com/vossjo/q-e 
  !  - https://stanford.edu/~vossj/slac/
  !
  !
  USE kinds
  !
  PRIVATE
  PUBLIC :: beef_energies, beef_print, beefxc, energies
  !
  real(DP), allocatable       :: beefxc(:), energies(:)
  real(DP)                    :: ldaxc
!# 24 "beef.f90"
 CONTAINS
!
!# 27 "beef.f90"
!
! obtain xc ensemble energies from non-selfconsistent calculations
! of the xc energies for perturbed BEEF expansion coefficents
! (provided by libbeef)
!
!-------------------------------------------------------------------------
SUBROUTINE beef_energies( )
!-------------------------------------------------------------------------
!# 36 "beef.f90"
  USE io_global,         ONLY  : stdout, ionode
  USE xc_lib,            ONLY  : xclib_dft_is
  USE ener,                 ONLY : vtxc, etxc
  USE scf,                  ONLY : rho, rho_core, rhog_core, v
  !
  USE beef_interface, ONLY: beefsetmode, beefrandinitdef, beefensemble
  !
  implicit none
  real(DP)                    :: ldaxc
  integer                     :: i
!# 47 "beef.f90"
  if (.not. allocated(beefxc)) allocate(beefxc(32))
  if (.not. allocated(energies)) allocate(energies(2000))
!# 50 "beef.f90"
  if (.not. xclib_dft_is('meta')) then
     do i=1,30
        !calculate exchange contributions in Legendre polynomial
        !basis
        call beefsetmode(i-1)
        CALL v_xc( rho, rho_core, rhog_core, beefxc(i), vtxc, v%of_r)
     enddo
        !calculate lda correlation contribution
        call beefsetmode(-3)
        CALL v_xc( rho, rho_core, rhog_core, beefxc(31), vtxc, v%of_r)
        !calculate pbe correlation contribution
        call beefsetmode(-2)
        CALL v_xc( rho, rho_core, rhog_core, beefxc(32), vtxc, v%of_r)
        !calculate lda xc energy
        call beefsetmode(-4)
        CALL v_xc( rho, rho_core, rhog_core, ldaxc, vtxc, v%of_r )
        !restore original, unperturbed xc potential and energy
        call beefsetmode(-1)
        CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
  else
     do i=1,30
        !calculate exchange contributions in Legendre polynomial
        !basis
        call beefsetmode(i-1)
        CALL v_xc_meta( rho, rho_core, rhog_core, beefxc(i), vtxc,v%of_r,v%kin_r )
     enddo
       !calculate lda correlation contribution
       call beefsetmode(-3)
       CALL v_xc_meta( rho, rho_core, rhog_core, beefxc(31), vtxc,v%of_r,v%kin_r )
       !calculate pbe correlation contribution
       call beefsetmode(-2)
       CALL v_xc_meta( rho, rho_core, rhog_core, beefxc(32), vtxc,v%of_r,v%kin_r )
       !calculate ldaxc energy
       call beefsetmode(-4)
       CALL v_xc_meta( rho, rho_core, rhog_core, ldaxc, vtxc,v%of_r,v%kin_r )
       !restore original, unperturbed xc potential and energy
       call beefsetmode(-1)
       CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc,v%of_r,v%kin_r )
  endif
  call beefrandinitdef
  !subtract LDA xc from exchange contributions
  do i=1,32
     beefxc(i) = beefxc(i)-ldaxc
  enddo
  beefxc(32) = beefxc(32)+beefxc(31)
!# 96 "beef.f90"
  call beefensemble(beefxc, energies)
  if (.NOT. ionode) RETURN
!# 99 "beef.f90"
  if ( ionode ) call beef_print( )
!# 101 "beef.f90"
END SUBROUTINE beef_energies
!# 103 "beef.f90"
!-------------------------------------------------------------------------
SUBROUTINE beef_print( )
!-------------------------------------------------------------------------
!# 107 "beef.f90"
  USE io_global,         ONLY  : stdout, ionode
  
  implicit none
  integer                     :: i
!# 112 "beef.f90"
  if (.NOT. ionode) RETURN
!# 114 "beef.f90"
  WRITE(*,*) "BEEFens 2000 ensemble energies"
  do i=1,2000
     WRITE(*, "(E35.15)") energies(i)
  enddo
  WRITE(*,*)
  WRITE(*,*) "BEEF-vdW xc energy contributions"
  do i=1,32
     WRITE(*,*) i, ": ", beefxc(i)
  enddo
!# 124 "beef.f90"
END SUBROUTINE beef_print
!# 126 "beef.f90"
END MODULE beef
