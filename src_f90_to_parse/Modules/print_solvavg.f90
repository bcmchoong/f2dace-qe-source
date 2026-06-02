!# 1 "print_solvavg.f90"
!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE print_solvavg(rismt, ext, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... print 3D-RISM's or Laue-RISM's correlations as planar average.
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_rank, mp_sum, mp_max
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvavg,   ONLY : solvavg_init, solvavg_clear, solvavg_print
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(IN)  :: rismt
  CHARACTER(LEN=*), INTENT(IN)  :: ext
  INTEGER,          INTENT(OUT) :: ierr
  !
  INTEGER            :: io_group_id
  INTEGER            :: my_group_id
  INTEGER            :: ista
  CHARACTER(LEN=256) :: filave
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... get process info.
  my_group_id = mp_rank(rismt%mp_site%inter_sitg_comm)
  !
  ! ... find the index of the group which includes ionode
  io_group_id = 0
  IF (ionode) THEN
    io_group_id = my_group_id
  END IF
  CALL mp_sum(io_group_id, rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(io_group_id, rismt%mp_site%inter_sitg_comm)
  !
  ! ... init solvavg
  IF (my_group_id == io_group_id) THEN
    IF (rismt%itype == ITYPE_3DRISM) THEN
      CALL solvavg_init(rismt%dfft, rismt%mp_site%intra_sitg_comm, .FALSE.)
    ELSE !IF (rismt%itype == ITYPE_LAUERISM) THEN
      CALL solvavg_init(rismt%lfft, rismt%mp_site%intra_sitg_comm, .FALSE.)
    END IF
  END IF
  !
  ! ... put data to solvavg
  IF (rismt%itype == ITYPE_3DRISM) THEN
    CALL print_solvavg_3drism(rismt, io_group_id, my_group_id)
  ELSE !IF (rismt%itype == ITYPE_LAUERISM) THEN
    CALL print_solvavg_lauerism(rismt, io_group_id, my_group_id)
  END IF
  !
  ! ... print solvavg
  IF (my_group_id == io_group_id) THEN
    filave = TRIM(tmp_dir) // TRIM(prefix) // '.' // TRIM(ext)
    CALL solvavg_print(filave, &
    & 'solvent densities and electrostatic potentials which act on electron', ista)
    ista = ABS(ista)
  ELSE
    ista = 0
  END IF
  !
  CALL mp_max(ista, rismt%mp_site%inter_sitg_comm)
  !
  IF (ista /= 0) THEN
    CALL errore('print_solvavg', 'cannot write file' // TRIM(filave), ista)
  END IF
  !
  ! ... finalize solvavg
  IF (my_group_id == io_group_id) THEN
    CALL solvavg_clear()
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE print_solvavg
!
!---------------------------------------------------------------------------
SUBROUTINE print_solvavg_3drism(rismt, io_group_id, my_group_id)
  !---------------------------------------------------------------------------
  !
  USE constants,      ONLY : RYTOEV, K_BOLTZMANN_RY, BOHR_RADIUS_ANGS
  USE control_flags,  ONLY : gamma_only
  USE fft_interfaces, ONLY : invfft
  USE io_global,      ONLY : ionode
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum, mp_get, mp_barrier
  USE rism,           ONLY : rism_type
  USE solvavg,        ONLY : solvavg_size, solvavg_put, solvavg_add
  USE solvmol,        ONLY : solVs, get_nuniq_in_solVs, &
                           & iuniq_to_isite, iuniq_to_nsite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN) :: rismt
  INTEGER,         INTENT(IN) :: io_group_id
  INTEGER,         INTENT(IN) :: my_group_id
  !
  REAL(DP) :: beta
  !
  INTEGER, PARAMETER :: LEN_SATOM = 6
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... put data
  CALL put_solvent()
  CALL put_solute()
  CALL put_vtotal()
  CALL put_solvent_g()
  CALL put_solvent_h()
  CALL put_solvent_c()
  CALL put_solvent_u()
  CALL put_solvent_t()
  !
CONTAINS
  !
  SUBROUTINE put_vtotal()
    IMPLICIT NONE
    !
    INTEGER                  :: ig
    INTEGER                  :: ir
    REAL(DP),    ALLOCATABLE :: vpot(:)
    COMPLEX(DP), ALLOCATABLE :: auxs(:)
    !
    IF (rismt%nr < 1) THEN
      RETURN
    END IF
    !
    IF (rismt%ng < 1 .OR. rismt%ng < rismt%gvec%ngm) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ! ... allocate memory
      ALLOCATE(vpot(rismt%nr))
      IF (rismt%dfft%nnr > 0) THEN
        ALLOCATE(auxs(rismt%dfft%nnr))
      END IF
      !
      ! ... Vsolute
      vpot = -rismt%vsr * RYTOEV  ! acting on electron
      CALL solvavg_put('Avg v_total (eV)', .FALSE., vpot)
      !
      vpot = -rismt%vlr * RYTOEV  ! acting on electron
      CALL solvavg_add(solvavg_size(),     .FALSE., vpot)
      !
      ! ... Vsolv
      IF (rismt%dfft%nnr > 0) THEN
        auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rismt%gvec%ngm
        auxs(rismt%dfft%nl(ig)) = rismt%vpot(ig)
      END DO
      IF (gamma_only) THEN
        DO ig = rismt%gvec%gstart, rismt%gvec%ngm
          auxs(rismt%dfft%nlm(ig)) = CONJG(auxs(rismt%dfft%nl(ig)))
        END DO
      END IF
      !
      IF (rismt%dfft%nnr > 0) THEN
        CALL invfft('Rho', auxs, rismt%dfft)
      END IF
      !
      vpot = 0.0_DP
      DO ir = 1, rismt%dfft%nnr
        vpot(ir) = -DBLE(auxs(ir)) * RYTOEV  ! acting on electron
      END DO
      !
      CALL solvavg_add(solvavg_size(),     .FALSE., vpot)
      !
      ! ... deallocate memory
      DEALLOCATE(vpot)
      IF (rismt%dfft%nnr > 0) THEN
        DEALLOCATE(auxs)
      END IF
      !
    END IF
    !
  END SUBROUTINE put_vtotal
  !
  SUBROUTINE put_solute()
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE :: vpot(:)
    !
    IF (rismt%nr < 1) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ALLOCATE(vpot(rismt%nr))
!# 221 "print_solvavg.f90"
      !
      ! ... Vsolute
      vpot = -rismt%vsr * RYTOEV  ! acting on electron
      CALL solvavg_put('Avg v_solute (eV)', .FALSE., vpot)
      !
      vpot = -rismt%vlr * RYTOEV  ! acting on electron
      CALL solvavg_add(solvavg_size(),       .FALSE., vpot)
      !
      DEALLOCATE(vpot)
      !
    END IF
    !
  END SUBROUTINE put_solute
  !
  SUBROUTINE put_solvent()
    IMPLICIT NONE
    !
    INTEGER                  :: ig
    INTEGER                  :: ir
    REAL(DP),    ALLOCATABLE :: rhor(:)
    REAL(DP),    ALLOCATABLE :: vpor(:)
    COMPLEX(DP), ALLOCATABLE :: auxs(:)
    !
    IF (rismt%ng < 1 .OR. rismt%ng < rismt%gvec%ngm) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ! ... allocate memory
      IF (rismt%dfft%nnr > 0) THEN
        ALLOCATE(rhor(rismt%dfft%nnr))
        ALLOCATE(vpor(rismt%dfft%nnr))
        ALLOCATE(auxs(rismt%dfft%nnr))
      END IF
      !
      ! ... Rho
      IF (rismt%dfft%nnr > 0) THEN
        auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rismt%gvec%ngm
        auxs(rismt%dfft%nl(ig)) = rismt%rhog(ig)
      END DO
      IF (gamma_only) THEN
        DO ig = rismt%gvec%gstart, rismt%gvec%ngm
          auxs(rismt%dfft%nlm(ig)) = CONJG(auxs(rismt%dfft%nl(ig)))
        END DO
      END IF
      !
      IF (rismt%dfft%nnr > 0) THEN
        CALL invfft('Rho', auxs, rismt%dfft)
      END IF
      !
      DO ir = 1, rismt%dfft%nnr
        rhor(ir) = DBLE(auxs(ir)) / BOHR_RADIUS_ANGS
      END DO
      !
      IF (rismt%dfft%nnr > 0) THEN
        CALL solvavg_put('Tot chg (e/A)', .TRUE., rhor)
      END IF
      !
      ! ... Vsolv
      IF (rismt%dfft%nnr > 0) THEN
        auxs = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rismt%gvec%ngm
        auxs(rismt%dfft%nl(ig)) = rismt%vpot(ig)
      END DO
      IF (gamma_only) THEN
        DO ig = rismt%gvec%gstart, rismt%gvec%ngm
          auxs(rismt%dfft%nlm(ig)) = CONJG(auxs(rismt%dfft%nl(ig)))
        END DO
      END IF
      !
      IF (rismt%dfft%nnr > 0) THEN
        CALL invfft('Rho', auxs, rismt%dfft)
      END IF
      !
      DO ir = 1, rismt%dfft%nnr
        vpor(ir) = -DBLE(auxs(ir)) * RYTOEV  ! acting on electron
      END DO
      !
      IF (rismt%dfft%nnr > 0) THEN
        CALL solvavg_put('Avg v_solvent (eV)', .FALSE., vpor)
      END IF
      !
      ! ... deallocate memory
      IF (rismt%dfft%nnr > 0) THEN
        DEALLOCATE(rhor)
        DEALLOCATE(vpor)
        DEALLOCATE(auxs)
      END IF
      !
    END IF
    !
  END SUBROUTINE put_solvent
  !
  SUBROUTINE put_solvent_g()
    IMPLICIT NONE
    !
    INTEGER                  :: nq
    INTEGER                  :: iq
    INTEGER                  :: nv
    INTEGER                  :: iv
    INTEGER                  :: isolV
    INTEGER                  :: iatom
    CHARACTER(LEN=LEN_SATOM) :: satom
    INTEGER                  :: owner_group_id
    REAL(DP)                 :: rhov
    REAL(DP), ALLOCATABLE    :: rhor(:)
    !
    IF (rismt%nr < 1) THEN
      RETURN
    END IF
    !
    nq = get_nuniq_in_solVs()
    !
    ALLOCATE(rhor(rismt%nr))
    !
    ! ... Guv
    DO iq = 1, nq
      iv    = iuniq_to_isite(1, iq)
      nv    = iuniq_to_nsite(iq)
      isolV = isite_to_isolV(iv)
      iatom = isite_to_iatom(iv)
      rhov  = solVs(isolV)%density
      satom = ADJUSTL(solVs(isolV)%aname(iatom))
      !
      IF (rismt%mp_site%isite_start <= iq .AND. iq <= rismt%mp_site%isite_end) THEN
        owner_group_id = my_group_id
        rhor = rismt%gr(:, iq - rismt%mp_site%isite_start + 1) * (DBLE(nv) * rhov / BOHR_RADIUS_ANGS)
      ELSE
        owner_group_id = 0
        rhor = 0.0_DP
      END IF
      !
      CALL mp_sum(owner_group_id, rismt%mp_site%inter_sitg_comm)
      CALL mp_get(rhor, rhor, my_group_id, io_group_id, &
                & owner_group_id, iq, rismt%mp_site%inter_sitg_comm)
      !
      IF (my_group_id == io_group_id) THEN
        CALL solvavg_put('Tot rho_'// TRIM(satom) //' (1/A)', .TRUE., rhor)
      END IF
      !
      CALL mp_barrier(rismt%mp_site%inter_sitg_comm)
    END DO
    !
    DEALLOCATE(rhor)
    !
  END SUBROUTINE put_solvent_g
  !
  SUBROUTINE put_solvent_h()
    IMPLICIT NONE
!# 422 "print_solvavg.f90"
  END SUBROUTINE put_solvent_h
  !
  SUBROUTINE put_solvent_c()
    IMPLICIT NONE
!# 501 "print_solvavg.f90"
  END SUBROUTINE put_solvent_c
  !
  SUBROUTINE put_solvent_u()
    IMPLICIT NONE
!# 632 "print_solvavg.f90"
  END SUBROUTINE put_solvent_u
  !
  SUBROUTINE put_solvent_t()
    IMPLICIT NONE
!# 686 "print_solvavg.f90"
  END SUBROUTINE put_solvent_t
  !
END SUBROUTINE print_solvavg_3drism
!
!---------------------------------------------------------------------------
SUBROUTINE print_solvavg_lauerism(rismt, io_group_id, my_group_id)
  !---------------------------------------------------------------------------
  !
  USE constants,      ONLY : RYTOEV, K_BOLTZMANN_RY, BOHR_RADIUS_ANGS
  USE control_flags,  ONLY : gamma_only
  USE fft_interfaces, ONLY : invfft
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : fw_lauefft_2xy
  USE mp,             ONLY : mp_sum, mp_get, mp_barrier
  USE rism,           ONLY : rism_type
  USE rism1d_facade,  ONLY : rism1t
  USE solvavg,        ONLY : solvavg_size, solvavg_put, solvavg_add
  USE solvmol,        ONLY : solVs, get_nuniq_in_solVs, &
                           & iuniq_to_isite, iuniq_to_nsite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN) :: rismt
  INTEGER,         INTENT(IN) :: io_group_id
  INTEGER,         INTENT(IN) :: my_group_id
  !
  REAL(DP) :: beta
  !
  INTEGER, PARAMETER :: LEN_SATOM = 6
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... put data
  CALL put_solvent()
  CALL put_solvent_pbc()
  CALL put_solute()
  CALL put_vtotal()
  CALL put_solvent_g()
  CALL put_solvent_h()
  CALL put_solvent_c()
  CALL put_solvent_u()
  CALL put_solvent_t()
  !
CONTAINS
  !
  SUBROUTINE put_vtotal()
    IMPLICIT NONE
    !
    REAL(DP),    ALLOCATABLE :: vpot(:)
    COMPLEX(DP), ALLOCATABLE :: vpol(:)
    !
    IF (rismt%nr < 1 .OR. (rismt%nrzl * rismt%ngxy) < 1) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ALLOCATE(vpot(rismt%nr))
      ALLOCATE(vpol(rismt%ngxy * rismt%nrzl))
      !
      ! ... Vsolute, R-space + Laue-rep.
      vpot = -rismt%vsr * RYTOEV  ! acting on electron
      CALL solvavg_put('Avg v_total (eV)', .FALSE., vpot)
      !
      vpol = -rismt%vlgz * RYTOEV  ! acting on electron
      CALL solvavg_add(solvavg_size(),     .FALSE., vpol, rismt%nrzl, .TRUE.)
      !
      ! ... Vsolv, Laue-rep.
      vpol = -rismt%vpot * RYTOEV  ! acting on electron
      CALL solvavg_add(solvavg_size(),     .FALSE., vpol, rismt%nrzl, .TRUE.)
      !
      DEALLOCATE(vpot)
      DEALLOCATE(vpol)
      !
    END IF
    !
  END SUBROUTINE put_vtotal
  !
  SUBROUTINE put_solute()
    IMPLICIT NONE
    !
    REAL(DP),    ALLOCATABLE :: vpot(:)
    COMPLEX(DP), ALLOCATABLE :: vpol(:)
    !
    IF (rismt%nr < 1 .OR. (rismt%nrzl * rismt%ngxy) < 1) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ALLOCATE(vpot(rismt%nr))
      ALLOCATE(vpol(rismt%ngxy * rismt%nrzl))
!# 795 "print_solvavg.f90"
      !
      ! ... Vsolute, R-space + Laue-rep.
      vpot = -rismt%vsr * RYTOEV  ! acting on electron
      CALL solvavg_put('Avg v_solute (eV)', .FALSE., vpot)
      !
      vpol = -rismt%vlgz * RYTOEV  ! acting on electron
      CALL solvavg_add(solvavg_size(),      .FALSE., vpol, rismt%nrzl, .TRUE.)
      !
      DEALLOCATE(vpot)
      DEALLOCATE(vpol)
      !
    END IF
    !
  END SUBROUTINE put_solute
  !
  SUBROUTINE put_solvent()
    IMPLICIT NONE
    !
    COMPLEX(DP), ALLOCATABLE :: tmpl(:)
    !
    IF ((rismt%nrzl * rismt%ngxy) < 1) THEN
      RETURN
    END IF
    !
    IF (my_group_id == io_group_id) THEN
      !
      ALLOCATE(tmpl(rismt%ngxy * rismt%nrzl))
      !
      ! ... Rho, Laue-rep.
      tmpl = rismt%rhog / BOHR_RADIUS_ANGS
      CALL solvavg_put('Tot chg (e/A)', .TRUE.,  tmpl, rismt%nrzl, .TRUE.)
!# 829 "print_solvavg.f90"
      !
      ! ... Vsolv, Laue-rep.
      tmpl = -rismt%vpot * RYTOEV  ! acting on electron
      CALL solvavg_put('Avg v_solvent (eV)',    .FALSE., tmpl, rismt%nrzl, .TRUE.)
!# 836 "print_solvavg.f90"
      !
      DEALLOCATE(tmpl)
      !
    END IF
    !
  END SUBROUTINE put_solvent
  !
  SUBROUTINE put_solvent_pbc()
    IMPLICIT NONE
!# 924 "print_solvavg.f90"
  END SUBROUTINE put_solvent_pbc
  !
  SUBROUTINE put_solvent_g()
    IMPLICIT NONE
    !
    INTEGER                  :: nq
    INTEGER                  :: iq
    INTEGER                  :: iiq
    INTEGER                  :: nv
    INTEGER                  :: iv
    INTEGER                  :: isolV
    INTEGER                  :: iatom
    CHARACTER(LEN=LEN_SATOM) :: satom
    INTEGER                  :: owner_group_id
    INTEGER                  :: iz
    INTEGER                  :: iiz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy
    INTEGER                  :: kgxy
    REAL(DP)                 :: rhov_right
    REAL(DP)                 :: rhov_left
    REAL(DP),    ALLOCATABLE :: rhor(:)
    COMPLEX(DP), ALLOCATABLE :: rhol(:)
    COMPLEX(DP), ALLOCATABLE :: ggz(:,:)
    !
    IF ((rismt%nrzs * rismt%ngxy) < 1) THEN
      RETURN
    END IF
    !
    IF ((rismt%nrzl * rismt%ngxy) < 1) THEN
      RETURN
    END IF
!# 962 "print_solvavg.f90"
    !
    nq = get_nuniq_in_solVs()
    !
    IF (rismt%nsite > 0) THEN
      ALLOCATE(ggz(rismt%nrzs * rismt%ngxy, rismt%nsite))
    END IF
    ALLOCATE(rhol(rismt%nrzl * rismt%ngxy))
!# 972 "print_solvavg.f90"
    !
    ! ... gr -> ggz
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq = iq - rismt%mp_site%isite_start + 1
      ggz(:, iiq) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      IF (rismt%dfft%nnr > 0) THEN
        CALL fw_lauefft_2xy(rismt%lfft, rismt%gr(:, iiq), ggz(:, iiq), rismt%nrzs, 1)
      END IF
    END DO
    !
    ! ... Guv, Laue-rep.
    DO iq = 1, nq
      iv         = iuniq_to_isite(1, iq)
      nv         = iuniq_to_nsite(iq)
      isolV      = isite_to_isolV(iv)
      iatom      = isite_to_iatom(iv)
      rhov_right = solVs(isolV)%density
      rhov_left  = solVs(isolV)%subdensity
      satom      = ADJUSTL(solVs(isolV)%aname(iatom))
      !
      IF (rismt%mp_site%isite_start <= iq .AND. iq <= rismt%mp_site%isite_end) THEN
        owner_group_id = my_group_id
        iiq = iq - rismt%mp_site%isite_start + 1
        !
        rhol = rismt%hsgz(:, iq - rismt%mp_site%isite_start + 1) &
           & + rismt%hlgz(:, iq - rismt%mp_site%isite_start + 1)
        IF (rismt%lfft%gxystart > 1) THEN
          rhol(1:rismt%nrzl) = rhol(1:rismt%nrzl) + CMPLX(1.0_DP, 0.0_DP, kind=DP)
        END IF
        !
        DO igxy = 1, rismt%ngxy
          jgxy = (igxy - 1) * rismt%nrzl
          kgxy = (igxy - 1) * rismt%nrzs
          DO iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
            iiz = iz - rismt%lfft%izcell_start + 1
            rhol(iz + jgxy) = ggz(iiz + kgxy, iiq)
          END DO
        END DO
        !
        DO igxy = 1, rismt%ngxy
          jgxy = (igxy - 1) * rismt%nrzl
          DO iz = 1, rismt%lfft%izleft_gedge
            rhol(iz + jgxy) = rhol(iz + jgxy) * rhov_left
          END DO
          DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
            rhol(iz + jgxy) = rhol(iz + jgxy) * rhov_right
          END DO
        END DO
        !
        rhol = rhol * (DBLE(nv) / BOHR_RADIUS_ANGS)
        !
      ELSE
        owner_group_id = 0
        rhol = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      !
      CALL mp_sum(owner_group_id, rismt%mp_site%inter_sitg_comm)
      CALL mp_get(rhol, rhol, my_group_id, io_group_id, &
                & owner_group_id, iq, rismt%mp_site%inter_sitg_comm)
      !
      IF (my_group_id == io_group_id) THEN
        CALL solvavg_put('Tot rho_'// TRIM(satom) //' (1/A)', .TRUE., rhol, rismt%nrzl, .TRUE.)
!# 1037 "print_solvavg.f90"
      END IF
      !
      CALL mp_barrier(rismt%mp_site%inter_sitg_comm)
    END DO
!# 1069 "print_solvavg.f90"
    !
    IF (rismt%nsite > 0) THEN
      DEALLOCATE(ggz)
    END IF
    DEALLOCATE(rhol)
!# 1077 "print_solvavg.f90"
    !
  END SUBROUTINE put_solvent_g
  !
  SUBROUTINE put_solvent_h()
    IMPLICIT NONE
!# 1279 "print_solvavg.f90"
  END SUBROUTINE put_solvent_h
  !
  SUBROUTINE put_solvent_c()
    IMPLICIT NONE
!# 1686 "print_solvavg.f90"
  END SUBROUTINE put_solvent_c
  !
  SUBROUTINE put_solvent_u()
    IMPLICIT NONE
!# 1804 "print_solvavg.f90"
  END SUBROUTINE put_solvent_u
  !
  SUBROUTINE put_solvent_t()
    IMPLICIT NONE
!# 1858 "print_solvavg.f90"
  END SUBROUTINE put_solvent_t
  !
END SUBROUTINE print_solvavg_lauerism
