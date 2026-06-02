!# 1 "xc_beef_interface.f90"
!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE beef_interface
    !----------------------------------------------------------------------
    !! This module contains fortran wrapper to BEEF library functions.
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    !
    PUBLIC :: beefx, beeflocalcorr, beeflocalcorrspin, beefsetmode, &
        beefrandinit, beefrandinitdef, beefensemble, beef_set_type
    !
!# 22 "xc_beef_interface.f90"
    INTERFACE
    !
    SUBROUTINE beefx( r, g, e, dr, dg, addlda ) BIND(C, NAME="beefx_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE)            :: r, g, e, dr, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beefx
    !
    SUBROUTINE beeflocalcorr( r, g, e, dr, dg, addlda) BIND(C, NAME="beeflocalcorr_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE), INTENT(INOUT) :: r, g, e, dr, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beeflocalcorr
    !
    SUBROUTINE beeflocalcorrspin(r, z, g, e, drup, drdown, dg, addlda) BIND(C, NAME="beeflocalcorrspin_")
        USE iso_c_binding
        !$acc routine seq
        REAL (C_DOUBLE), INTENT(INOUT) :: r, z, g, e, drup, drdown, dg
        INTEGER(C_INT), INTENT(IN) :: addlda
    END SUBROUTINE beeflocalcorrspin
    !
    SUBROUTINE beefsetmode(mode) BIND(C, NAME="beefsetmode_")
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: mode
    END SUBROUTINE beefsetmode
    !
    SUBROUTINE beefrandinit(seed) BIND(C, NAME="beefrandinit_")
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: seed
    END SUBROUTINE beefrandinit
    !
    SUBROUTINE beefrandinitdef() BIND(C, NAME="beefrandinitdef_")
    END SUBROUTINE beefrandinitdef
    !
    SUBROUTINE beefensemble(beefxc, ensemble) BIND(C, NAME="beefensemble_")
        USE iso_c_binding
        REAL (C_DOUBLE), INTENT(INOUT) :: beefxc(*), ensemble(*)
    END SUBROUTINE beefensemble
    !
    FUNCTION beef_set_type_interface(tbeef, ionode) &
            BIND(C,name="beef_set_type_") RESULT(r)
        USE iso_c_binding
        INTEGER(C_INT), INTENT(IN) :: tbeef, ionode
        INTEGER(C_INT)             :: r
    END FUNCTION beef_set_type_interface
    !
    END INTERFACE
    !
!# 73 "xc_beef_interface.f90"
    !
    CONTAINS
    ! ====================================================================
    !
!# 78 "xc_beef_interface.f90"
    FUNCTION beef_set_type(tbeef, ionode) RESULT(r)
        INTEGER, INTENT(IN) :: tbeef
        LOGICAL, INTENT(IN) :: ionode
        LOGICAL             :: r
        ! ... local variables ...
        INTEGER             :: ionode_ = 0
        INTEGER             :: r_
        !
        IF ( ionode ) ionode_ = 1
        !
        r_ = beef_set_type_interface(tbeef, ionode_)
        !
        IF ( r_ /= 0 ) THEN
            r = .TRUE.
        ELSE
            r = .FALSE.
        END IF
        !
    END FUNCTION beef_set_type
    !
!# 163 "xc_beef_interface.f90"
    !
END MODULE beef_interface
!
