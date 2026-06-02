!# 1 "oscdft_extended.f90"
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE oscdft_nsg3 (nsg, nsnew)
   !
   !! This routine copies the diagonal components of the complex 
   !! generalized occupation matrix nsg to a real array nsnew
   !! that is used to build the constraint - replaces oscdft_nsg(iflag=3)
   !
   USE kinds,           ONLY : DP
   USE parameters,      ONLY : ntypx
   USE io_global,       ONLY : stdout
   USE ions_base,       ONLY : nat, ityp
   USE lsda_mod,        ONLY : nspin
   USE upf_params,      ONLY : lqmax
   USE ldaU,            ONLY : max_num_neighbors, ldmx_tot, neighood, &
                               Hubbard_l, Hubbard_lmax
!# 25 "oscdft_extended.f90"
   !
   IMPLICIT NONE
   COMPLEX(DP), INTENT(IN)  :: nsg (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   REAL(dp), INTENT(OUT) :: nsnew(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin,nat)
   !
   INTEGER :: na, na1, nt, viz, ldim, is, m1, m2
   !
!# 57 "oscdft_extended.f90"
   RETURN
   !
END SUBROUTINE oscdft_nsg3
!
SUBROUTINE oscdft_nsg (lflag,nsg)
   !
   !! This routine adjusts (modifies) the nsg based on constraints
   !! If lflag=1, then copy ctx%inp%occupation to nsg
   !! If lflag=2, then copy nsg to ctx%inp%occupation
   !! If lflag=4, like lflag=1 but it does not nullify the occupations for 
   !!             Hubbard atoms atoms to which we do not apply the constraints
   !! Case lflag=3 implemented in oscdft_nsg3, no longer here
   !
   USE kinds,           ONLY : DP
   USE parameters,      ONLY : ntypx
   USE io_global,       ONLY : stdout
   USE ions_base,       ONLY : nat, ityp
   USE lsda_mod,        ONLY : nspin
   USE upf_params,      ONLY : lqmax
   USE ldaU,            ONLY : max_num_neighbors, ldmx_tot, neighood, &
                               Hubbard_l, Hubbard_lmax
!# 81 "oscdft_extended.f90"
   !
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: lflag
   COMPLEX(DP), INTENT(INOUT) :: nsg (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   !
   INTEGER :: na, na1, nt, viz, ldim, is, m1, m2
   LOGICAL :: found
   !
!# 136 "oscdft_extended.f90"
   RETURN
   !
END SUBROUTINE oscdft_nsg
!# 140 "oscdft_extended.f90"
SUBROUTINE oscdft_v_constraint_extended (nsg, vtot, etot)
   !
   !! Computes the contribution to the potential from the Lagrange multipliers used
   !! to constrain the occupation matrix to the target (DFT+U+V case).
   !! Here we are applying the constraint only using the diagonal (onsite) component 
   !! of the generalized occupaion matrix nsg.
   !
   USE kinds,           ONLY : DP
   USE parameters,      ONLY : ntypx
   USE ions_base,       ONLY : nat, ityp
   USE lsda_mod,        ONLY : nspin
   USE io_global,       ONLY : stdout
   USE control_flags,   ONLY : iverbosity
   USE ldaU,            ONLY : ldim_u, ldmx_tot, max_num_neighbors, neighood 
!# 158 "oscdft_extended.f90"
   IMPLICIT NONE
!# 160 "oscdft_extended.f90"
   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(INOUT) :: vtot(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   REAL(DP), INTENT(INOUT) :: etot
   REAL(DP) :: ec
   INTEGER :: is, na, na1, na2, viz, nt1, m1, m2
   INTEGER, EXTERNAL :: find_viz
   !
!# 205 "oscdft_extended.f90"
   RETURN
   !
END SUBROUTINE oscdft_v_constraint_extended
