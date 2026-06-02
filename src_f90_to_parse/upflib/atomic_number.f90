!# 1 "atomic_number.f90"
!
! Copyright (C) 2004-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ------------------------------------------------------------------
function atomic_number(atm)
  ! ------------------------------------------------------------------
  !
  use upf_utils, only: capital, lowercase, isnumeric
  implicit none
  character(len=*) :: atm
  integer :: atomic_number
!# 17 "atomic_number.f90"
  character(len=2) :: elements(120), atom
  data elements/' H',                              'He', &
                'Li','Be',' B',' C',' N',' O',' F','Ne', &
                'Na','Mg','Al','Si',' P',' S','Cl','Ar', &
                ' K','Ca','Sc','Ti',' V','Cr','Mn',      &
                          'Fe','Co','Ni','Cu','Zn',      &
                          'Ga','Ge','As','Se','Br','Kr', &
                'Rb','Sr',' Y','Zr','Nb','Mo','Tc',      &
                          'Ru','Rh','Pd','Ag','Cd',      &
                          'In','Sn','Sb','Te',' I','Xe', &
                'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
                               'Tb','Dy','Ho','Er','Tm','Yb','Lu', &
                               'Hf','Ta',' W','Re','Os', &
                          'Ir','Pt','Au','Hg',           &
                          'Tl','Pb','Bi','Po','At','Rn', &
                'Fr','Ra','Ac','Th','Pa',' U','Np','Pu', &
                'Am','Cm','Bk','Cf','Es','Fm','Md','No', &
                'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
                'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Un','Ub' /
  integer :: n
!# 38 "atomic_number.f90"
  atom='  '
  if ( len(atm) == 1 ) then
!
! Case : atm='X'
!
     atom(2:2)=capital(atm(1:1))
  else if ( ( len_trim(atm) == 1 ) .or. ( isnumeric(atm(2:2)) ) .or. &
          ( atm(2:2) == '-' )    .or. ( atm(2:2) == '_' ) ) then
!
! Case : atm='X ', 'X_*', 'X-*', 'X[0-9]* '
!
     atom(2:2)=capital(atm(1:1))
  else if (atm(1:1) == ' ') then
!
! Case : atm=' X*'
!
     atom(2:2)=capital(atm(2:2))
  else
!
! Case : atm='XY*'
!
     atom(1:1)=capital(atm(1:1))
     atom(2:2)=lowercase(atm(2:2))
  end if
      
  do n=1, 120
     if ( atom == elements(n) ) then
        atomic_number=n
        return
     end if
  end do
!# 70 "atomic_number.f90"
  atomic_number = 0
  print '("WARNING: Atom ",a2," not found")', atom
!# 73 "atomic_number.f90"
end function atomic_number
! ------------------------------------------------------------------
function atom_name(atomic_number)
  ! ------------------------------------------------------------------
  !
  integer :: atomic_number
  character(len=2) :: atom_name
!# 81 "atomic_number.f90"
  character(len=2) :: elements(120)
  data elements/' H',                              'He', &
                'Li','Be',' B',' C',' N',' O',' F','Ne', &
                'Na','Mg','Al','Si',' P',' S','Cl','Ar', &
                ' K','Ca','Sc','Ti',' V','Cr','Mn',      &
                          'Fe','Co','Ni','Cu','Zn',      &
                          'Ga','Ge','As','Se','Br','Kr', &
                'Rb','Sr',' Y','Zr','Nb','Mo','Tc',      &
                          'Ru','Rh','Pd','Ag','Cd',      &
                          'In','Sn','Sb','Te',' I','Xe', &
                'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd', &
                               'Tb','Dy','Ho','Er','Tm','Yb','Lu', &
                               'Hf','Ta',' W','Re','Os', &
                          'Ir','Pt','Au','Hg',           &
                          'Tl','Pb','Bi','Po','At','Rn', &
                'Fr','Ra','Ac','Th','Pa',' U','Np','Pu', &
                'Am','Cm','Bk','Cf','Es','Fm','Md','No', &
                'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
                'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og','Un','Ub' /
!# 101 "atomic_number.f90"
  if (atomic_number < 1 .or. atomic_number > 120) then
     print *, 'Invalid atomic number: ',atomic_number
     atom_name = 'XX'
  else
     atom_name=elements(atomic_number)
  end if
  return
!# 109 "atomic_number.f90"
end function atom_name
