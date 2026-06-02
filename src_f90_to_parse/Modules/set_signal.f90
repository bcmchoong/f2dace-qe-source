!# 1 "set_signal.f90"
!
!# 168 "set_signal.f90"
MODULE set_signal
!! This module is a Fortran 2003 interface to the customize_signals.c C file
!! Compatible with Intel/PGI/Gcc(>=4.3) compilers
!# 172 "set_signal.f90"
USE io_global, ONLY : stdout
!# 174 "set_signal.f90"
CONTAINS
!# 176 "set_signal.f90"
! Place holders to employ when the signal trapping feature is disabled
SUBROUTINE signal_trap_init
  WRITE(stdout, FMT=*) "signal trapping disabled: compile with "
  WRITE(stdout, FMT=*) "-D__TRAP_SIGUSR1 to enable this feature"
END SUBROUTINE signal_trap_init
!# 182 "set_signal.f90"
FUNCTION signal_detected()
  LOGICAL::signal_detected
  signal_detected = .FALSE.
END FUNCTION signal_detected
!# 187 "set_signal.f90"
END MODULE set_signal
