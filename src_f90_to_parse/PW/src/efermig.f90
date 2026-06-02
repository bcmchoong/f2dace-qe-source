!# 1 "efermig.f90"
!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION efermig( et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk )
  !--------------------------------------------------------------------
  !! Finds the Fermi energy - Gaussian Broadening
  !! (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !! Improved bisection algorithm by Flaviano José dos Santos (EPFL)
  !! Functions not passed as arguments (some compilers don't like it) 
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  USE mp,        ONLY : mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: Ngauss
  !! type of smearing technique
  INTEGER, INTENT(IN) :: is
  !! spin label (0 or 1,2)
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  REAL(DP), INTENT(IN) :: wk(nks)
  !! weight of k points
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: Degauss
  !! smearing parameter
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP) :: efermig
  !! the Fermi energy
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-10, eps_cold_MP = 1.0d-2
  !! tolerance for the number of electrons, important for bisection 
  !! smaller tolerance for the number of electrons, important for M-P and Cold smearings
  INTEGER, PARAMETER :: maxiter = 300
  !
  REAL(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  REAL(DP), EXTERNAL :: sumkg, sumkg1, sumkg2
  !! Function to compute the number of electrons for a given energy
  !! Function to compute the first derivative of the number of electrons
  !! Function to compute the second derivative of the number of electrons
  INTEGER :: i, kpoint, Ngauss_
  INTEGER :: info, maxiter_aux
  !
  !  ... find (very safe) bounds for the Fermi energy:
  !  Elw = lowest, Eup = highest energy among all k-points.
  !  Works with distributed k-points, also if nks=0 on some processor
  !
  Elw = 1.0E+8
  Eup =-1.0E+8
  DO kpoint = 1, nks
    Elw = MIN( Elw, et(1,kpoint) )
    Eup = MAX( Eup, et(nbnd,kpoint) )
  ENDDO
  Eup = Eup + 10 * Degauss
  Elw = Elw - 10 * Degauss
  !
  ! ... find min and max across pools
  !
  CALL mp_max( eup, inter_pool_comm )
  CALL mp_min( elw, inter_pool_comm )
  !
!# 79 "efermig.f90"
  ! For M-P and cold smearings, perform a preliminary determination with the Gaussian broadening
  ! to obtain an initial guess
  if( Ngauss .NE. -99 ) then
    Ngauss_ = 0
  else ! If Fermi-Dirac smearing, the preliminary bisection can run with the F-D smearing itself 
    Ngauss_ = Ngauss
  end if
!# 87 "efermig.f90"
  maxiter_aux = maxiter
  
  call bisection_find_efermi( Elw, Eup, ef, eps, maxiter_aux, info)
!# 91 "efermig.f90"
  ! Error handling
  select case( info )
    case( 1 )
      IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
      WRITE( stdout, '(5x,"Warning: too many iterations in bisection" &
        &      5x,"Ef (eV) = ",f15.6," Num. electrons = ",f10.6)' ) &
        Ef * rytoev, num_electrons(Ef)
    case( 2 )
      call errore( 'efermig', 'internal error, cannot bracket Ef', 1 )
  end select
!# 102 "efermig.f90"
  ! If this initial guess already corresponds to the correct number of electron for the actual occupation function, the Fermi energy is found.
  Ngauss_ = Ngauss
  
  ! In case Ngauss = 0 or -99, the function returns here too.
  if( abs_num_electrons_minus_nelec(ef) < eps .or. Ngauss == 0 .or. Ngauss == -99) then 
    
    efermig = ef
    
    goto 98765
  end if
  
  ! If the initial prospected Ef did not provide the correct number of electrons, use Newton's methods to improve it.
  ! Use the prospected Ef as initial guess.
!# 116 "efermig.f90"
  maxiter_aux = maxiter
  
  if( Ngauss_ > 0  .or.  Ngauss_ == -1 ) then ! If methfessel-paxton method or Cold smearing method
!# 120 "efermig.f90"
!    call newton_minimization(dev1_sq_num_electrons, dev2_sq_num_electrons, ef, eps, maxiter_aux, info)
    call newton_minimization( ef, eps, maxiter_aux, info)
!# 123 "efermig.f90"
  end if
!# 125 "efermig.f90"
  ! Error handling
  select case( info )
    case( 1 )
      IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
      WRITE( stdout, '(5x,"Warning: too many iterations in Newtons minimization"/ &
         &      5x,"Ef (eV) = ",f15.6," Num. electrons = ",f10.6,"  Num. steps = ",i0)' ) &
         Ef * rytoev, num_electrons(Ef), maxiter
  end select
!# 134 "efermig.f90"
  if( (Ngauss_ == -1 .or. Ngauss_ >  0) .and. ( abs_num_electrons_minus_nelec(ef) < eps_cold_MP ) ) then
!# 136 "efermig.f90"
    efermig = ef
!# 138 "efermig.f90"
  else
    ! If Newton's minimization did not help. Just use bisection with the actual smearing, which reproduce the original behavior of this function
    Ngauss_ = Ngauss
    maxiter_aux = maxiter
!# 143 "efermig.f90"
!    call bisection_find_efermi(num_electrons_minus_nelec, Elw, Eup, ef, eps, maxiter_aux, info)
    call bisection_find_efermi( Elw, Eup, ef, eps, maxiter_aux, info)
!# 146 "efermig.f90"
    efermig = ef
!# 148 "efermig.f90"
    IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
    WRITE( stdout, '(5x,"Minimization algorithm failed to find Fermi energy: reverting to bisection",&
     & /,5x,"Possible cause: smearing is larger than the electronic band-gap.")' )
  end if
!# 153 "efermig.f90"
  98765 continue
!# 155 "efermig.f90"
  return
!# 157 "efermig.f90"
  contains
!# 159 "efermig.f90"
  function num_electrons_minus_nelec(x)
    real(DP), intent(in) :: x
    real(DP) :: num_electrons_minus_nelec
!# 163 "efermig.f90"
    num_electrons_minus_nelec = num_electrons(x) - nelec
  end function num_electrons_minus_nelec
!# 166 "efermig.f90"
  function num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: num_electrons
!# 170 "efermig.f90"
    num_electrons = sumkg( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function num_electrons
!# 173 "efermig.f90"
  function abs_num_electrons_minus_nelec(ef)
    real(DP), intent(in) :: ef
    real(DP) :: abs_num_electrons_minus_nelec
!# 177 "efermig.f90"
    abs_num_electrons_minus_nelec = abs(num_electrons_minus_nelec(ef))
  end function abs_num_electrons_minus_nelec
!# 180 "efermig.f90"
  function sq_num_electrons_minus_nelec(ef)
    real(DP), intent(in) :: ef
    real(DP) :: sq_num_electrons_minus_nelec
!# 184 "efermig.f90"
    sq_num_electrons_minus_nelec = (num_electrons_minus_nelec(ef))**2
  end function sq_num_electrons_minus_nelec
!# 187 "efermig.f90"
  function dev1_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev1_num_electrons
!# 191 "efermig.f90"
    dev1_num_electrons = sumkg1( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function dev1_num_electrons
!# 194 "efermig.f90"
  function dev2_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev2_num_electrons
!# 198 "efermig.f90"
    dev2_num_electrons = sumkg2( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function dev2_num_electrons
!# 201 "efermig.f90"
  function dev1_sq_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev1_sq_num_electrons
!# 205 "efermig.f90"
    dev1_sq_num_electrons = 2.d0 * num_electrons_minus_nelec(ef) * dev1_num_electrons(ef)
  end function dev1_sq_num_electrons
!# 208 "efermig.f90"
  function dev2_sq_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev2_sq_num_electrons
!# 212 "efermig.f90"
    dev2_sq_num_electrons = 2.d0 * ( (dev1_num_electrons(ef))**2 + num_electrons_minus_nelec(ef) * dev2_num_electrons(ef) )
  end function dev2_sq_num_electrons
!# 215 "efermig.f90"
  subroutine newton_minimization(x, tol, Nmax, info)
    real(DP),          intent(inout) :: x
    !! Initial guess in the entry. Solution in the exit
    real(DP),          intent(in)    :: tol
    integer,           intent(inout) :: Nmax
    integer,           intent(out)   :: info
    !! 0 = solution found; 1 = max number of step (Nmax) reached; 2 = second derivative is zero
!# 223 "efermig.f90"
    real(DP)                         :: abstol, x0, denominator, numerator, factor
    integer                          :: i
!# 226 "efermig.f90"
    abstol = abs(tol)
!# 228 "efermig.f90"
    x0 = x
!# 230 "efermig.f90"
    factor = 1.0d0
!# 232 "efermig.f90"
    do i = 1, Nmax
       
       numerator   = dev1_sq_num_electrons(x)
       denominator = abs(dev2_sq_num_electrons(x))
!# 237 "efermig.f90"
       ! Checking if the denominator is zero
       if( denominator > abstol ) then
          x = x0 - factor*numerator/denominator
!# 241 "efermig.f90"
          ! Checking if a stationary point was achieved
          if( abs(x0-x) < abstol .or. abs_num_electrons_minus_nelec(x) < abstol ) then
             info = 0
             Nmax = i
             return
          ! If a stationary point was not achieved, continue
          else
             x0 = x
          end if
!# 251 "efermig.f90"
       ! If denominator is zero, return an error
       else 
          info = 2
          return
       end if
    end do
!# 258 "efermig.f90"
    ! Checking if max number of steps was reached
    if( i > Nmax ) then
      info = 1
      return
    end if 
  end subroutine newton_minimization
!# 265 "efermig.f90"
  subroutine bisection_find_efermi( energy_lower_bound, energy_upper_bound, x, tol, Nmax, info)
    real(DP),      intent(in)    :: energy_lower_bound
    real(DP),      intent(in)    :: energy_upper_bound
    real(DP),      intent(out)   :: x
    !! Found Fermi energy at exit
    real(DP),      intent(in)    :: tol
    integer,       intent(inout) :: Nmax
    !! In entry: Max number of steps. In exit: number of step taken.
    integer,       intent(out)   :: info
    !! 0 = solution found; 1 = max number of step (Nmax) reached; 2 = cannot bracket root
!# 276 "efermig.f90"
    real(DP)                     :: abs_tol, fx, Elw_local, Eup_local
    integer                      :: i
!# 279 "efermig.f90"
    abs_tol = abs(tol)
!# 281 "efermig.f90"
    Elw_local = energy_lower_bound
    Eup_local = energy_upper_bound
!# 284 "efermig.f90"
    if( num_electrons_minus_nelec(Elw_local) > abs_tol .or. num_electrons_minus_nelec(Eup_local) < -abs_tol ) then
      info = 2
      return
    end if
!# 289 "efermig.f90"
    do i = 1, Nmax
!# 291 "efermig.f90"
      x = ( Eup_local + Elw_local ) * 0.5d0
      fx = num_electrons_minus_nelec(x)
!# 294 "efermig.f90"
      ! Was the root found?
      if( abs(fx) < abs_tol ) then
        info = 0
        Nmax = i
        return
      else
        ! Choosing new boundaries
        if( fx < -abs_tol ) then
          Elw_local = x
        else
          Eup_local = x
        end if
      end if
    end do
!# 309 "efermig.f90"
    ! Checking if max number of steps was reached
    if( i > Nmax ) then
      info = 1
      return
    end if 
  end subroutine bisection_find_efermi
END FUNCTION efermig
