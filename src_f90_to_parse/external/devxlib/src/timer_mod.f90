!# 1 "timer_mod.f90"
!
! interface to timing routines
!
module timer_m
  use iso_c_binding
  implicit none
  private
 
  interface 
    function wallclock() bind(C,name="wallclock_c")
       use iso_c_binding, only: c_double
       real(c_double) :: wallclock
    end function
  end interface
!# 16 "timer_mod.f90"
  public :: wallclock
!# 18 "timer_mod.f90"
end module timer_m
