!# 1 "esm_stres_mod.f90"
MODULE esm_stres_mod
!# 3 "esm_stres_mod.f90"
  USE kinds,    ONLY : DP
  USE esm_common_mod, ONLY : esm_w, esm_bc, &
                             mill_2d, imill_2d, ngm_2d, exp_erfc
  IMPLICIT NONE
!# 8 "esm_stres_mod.f90"
  ! Workaround for Cray bug - note that exp, cosh, sinh with complex argument
  ! are in the F2008 standard so qe_exp, qe_cosh, qe_sinh are no longer needed
!# 16 "esm_stres_mod.f90"
CONTAINS
!# 18 "esm_stres_mod.f90"
  !-----------------------------------------------------------------------
  !--------------ESM STRESS SUBROUTINE------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE esm_stres_har(sigmahar, rhog)
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    IMPLICIT NONE
    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)
!# 28 "esm_stres_mod.f90"
    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_har must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_har_bc1(sigmahar, rhog)
    CASE ('bc2')
      CALL esm_stres_har_bc2(sigmahar, rhog)
    CASE ('bc3')
      CALL esm_stres_har_bc3(sigmahar, rhog)
    CASE ('bc4')
      STOP 'esm_stres_har has not yet implemented for esm_bc = bc4'
    END SELECT
!# 41 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_har
!# 44 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewa(sigmaewa)
    !-----------------------------------------------------------------------
    !
    ! Calculates Ewald stresswith both G- and R-space terms.
    ! Determines optimal alpha. Should hopefully work for any structure.
    !
    USE kinds,     ONLY : DP
    USE constants, ONLY : tpi
    USE cell_base, ONLY : tpiba2
    USE ions_base, ONLY : zv, nat, ityp
    USE gvect,     ONLY : gcutm
    IMPLICIT NONE
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 58 "esm_stres_mod.f90"
    ! output: the ewald stress
    !
    !    here the local variables
    !
    INTEGER :: ia
    ! counter on atoms
!# 65 "esm_stres_mod.f90"
    REAL(DP) :: charge, alpha, upperbound
    ! total ionic charge in the cell
    ! alpha term in ewald sum
    ! the maximum radius to consider real space sum
    REAL(DP) :: sigmaewg(3, 3), sigmaewr(3, 3)
    ! ewald stress computed in reciprocal space
    ! ewald stress computed in real space
!# 73 "esm_stres_mod.f90"
    charge = sum(zv(ityp(:)))
!# 75 "esm_stres_mod.f90"
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    alpha = 2.9d0
    DO
      alpha = alpha - 0.1d0
      IF (alpha .le. 0.d0) CALL errore('esm_stres_ewa', 'optimal alpha not found', 1)
      upperbound = 2.d0*charge**2*sqrt(2.d0*alpha/tpi)* &
                   erfc(sqrt(tpiba2*gcutm/4.d0/alpha))
      IF (upperbound < 1.0d-7) EXIT
    END DO
!# 86 "esm_stres_mod.f90"
    ! G-space sum here.
    ! Determine if this processor contains G=0 and set the constant term
    CALL esm_stres_ewg(alpha, sigmaewg)
!# 90 "esm_stres_mod.f90"
    ! R-space sum here (only for the processor that contains G=0)
    CALL esm_stres_ewr(alpha, sigmaewr)
!# 93 "esm_stres_mod.f90"
    sigmaewa(:, :) = sigmaewg(:, :) + sigmaewr(:, :)
!# 95 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewa
!# 98 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_loclong(sigmaloclong, rhog)
    USE kinds,    ONLY : DP
    USE gvect,    ONLY : ngm
    IMPLICIT NONE
    REAL(DP), INTENT(out)   :: sigmaloclong(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)
!# 105 "esm_stres_mod.f90"
    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_loclong must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_loclong_bc1(sigmaloclong, rhog)
    CASE ('bc2')
      CALL esm_stres_loclong_bc2(sigmaloclong, rhog)
    CASE ('bc3')
      CALL esm_stres_loclong_bc3(sigmaloclong, rhog)
    CASE ('bc4')
      STOP 'esm_stres_loclong has not yet implemented for esm_bc = bc4'
    END SELECT
!# 118 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_loclong
!# 121 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_har_bc1(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 133 "esm_stres_mod.f90"
    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)
!# 136 "esm_stres_mod.f90"
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 151 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 157 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
!# 162 "esm_stres_mod.f90"
    ! initialize
    sigmahar(:, :) = 0.0d0
!# 165 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 171 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 174 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 183 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 186 "esm_stres_mod.f90"
      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig
!# 196 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 203 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 205 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 214 "esm_stres_mod.f90"
      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 224 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2
      END DO ! igz
!# 231 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 239 "esm_stres_mod.f90"
        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
          - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
          + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m
      END DO ! iz
!# 250 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 257 "esm_stres_mod.f90"
      ! add bare coulomn terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 263 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
!# 265 "esm_stres_mod.f90"
        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz
!# 270 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 275 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp
!# 283 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 294 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz
!# 299 "esm_stres_mod.f90"
      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 307 "esm_stres_mod.f90"
        rg3 = rhog3(1, imill_2d(0, 0))
        Vr(iz) = &
          -tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c
      END DO ! iz
!# 315 "esm_stres_mod.f90"
      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 338 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 348 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 351 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 358 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 366 "esm_stres_mod.f90"
      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 373 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz
!# 377 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 382 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
!# 389 "esm_stres_mod.f90"
    END IF ! imill_2d(0,0) > 0
!# 391 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)
!# 395 "esm_stres_mod.f90"
    CALL mp_sum(sigmahar, intra_bgrp_comm)
!# 397 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 403 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_har_bc1
!# 406 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_har_bc2(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 418 "esm_stres_mod.f90"
    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)
!# 421 "esm_stres_mod.f90"
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum1c, sum2c, sum2p, sum2m
    COMPLEX(DP) :: sum1sp, sum1sm, sum1cp, sum1cm, sum2sp, sum2sm
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 435 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 441 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
!# 447 "esm_stres_mod.f90"
    ! initialize
    sigmahar(:, :) = 0.0d0
!# 450 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 456 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 459 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 468 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 471 "esm_stres_mod.f90"
      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig
!# 481 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 488 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 490 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 499 "esm_stres_mod.f90"
      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      sum1sp = (0.d0, 0.d0)
      sum1sm = (0.d0, 0.d0)
      sum1cp = (0.d0, 0.d0)
      sum1cm = (0.d0, 0.d0)
      sum2sp = (0.d0, 0.d0)
      sum2sm = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 515 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
!# 517 "esm_stres_mod.f90"
        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2
!# 522 "esm_stres_mod.f90"
        sum1sp = sum1sp + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)
        sum1sm = sum1sm + rg3*QE_SINH(gp*z0 - ci*gz*z0)/(gp - ci*gz)
!# 525 "esm_stres_mod.f90"
        sum1cp = sum1cp + rg3*QE_COSH(gp*z0 + ci*gz*z0)/(gp + ci*gz)*z0
        sum1cm = sum1cm + rg3*QE_COSH(gp*z0 - ci*gz*z0)/(gp - ci*gz)*z0
!# 528 "esm_stres_mod.f90"
        sum2sp = sum2sp + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)**2
        sum2sm = sum2sm + rg3*QE_SINH(gp*z0 - ci*gz*z0)/(gp - ci*gz)**2
      END DO ! igz
!# 532 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 540 "esm_stres_mod.f90"
        !! BC1 terms
        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
          - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
          + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m
!# 551 "esm_stres_mod.f90"
        !! BC2 terms
        dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                             + dgp_deps(:, :)*( &
                             - tpi/gp**2*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1) &
                             + tpi/gp*((-z - 2*z1)*EXP(-gp*(z + 2*z1)) - z*EXP(+gp*z))/sinh(2*gp*z1) &
                             - tpi/gp*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)) &
                             * sum1sp &
                             + tpi/gp*(EXP(-gp*(z + 2*z1)) - EXP(+gp*z))/sinh(2*gp*z1) &
                             * (-delta(:, :)*sum1sp + dgp_deps(:, :)*(sum1cp - sum2sp)) &
                             + dgp_deps(:, :)*( &
                             - tpi/gp**2*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1) &
                             + tpi/gp*((+z - 2*z1)*EXP(+gp*(z - 2*z1)) + z*EXP(-gp*z))/sinh(2*gp*z1) &
                             - tpi/gp*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)) &
                             * sum1sm &
                             + tpi/gp*(EXP(+gp*(z - 2*z1)) - EXP(-gp*z))/sinh(2*gp*z1) &
                             * (-delta(:, :)*sum1sm + dgp_deps(:, :)*(sum1cm - sum2sm))
!# 568 "esm_stres_mod.f90"
      END DO ! iz
!# 570 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 577 "esm_stres_mod.f90"
      ! add bare couloum terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 583 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
!# 585 "esm_stres_mod.f90"
        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz
!# 590 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 595 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp
!# 603 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 614 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz
!# 619 "esm_stres_mod.f90"
      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 627 "esm_stres_mod.f90"
        rg3 = rhog3(1, imill_2d(0, 0))
!# 629 "esm_stres_mod.f90"
        !! BC1 terms
        Vr(iz) = &
          - tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c
!# 636 "esm_stres_mod.f90"
        !! BC2 terms
        Vr(iz) = Vr(iz) &
                 + tpi*z1*2*z0*rg3 - tpi*(-z/z1)*2*z0*sum1c
      END DO ! iz
!# 641 "esm_stres_mod.f90"
      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f1 = f1 &
           + tpi*z1*2*z0*rg3 - tpi*(-z_r/z1)*2*z0*sum1c
!# 651 "esm_stres_mod.f90"
      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f2 = f2 &
           + tpi*z1*2*z0*rg3 - tpi*(-z_l/z1)*2*z0*sum1c
!# 658 "esm_stres_mod.f90"
      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f3 = f3 &
           - tpi*(-1.0d0/z1)*2*z0*sum1c
!# 663 "esm_stres_mod.f90"
      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      f4 = f4 &
           - tpi*(-1.0d0/z1)*2*z0*sum1c
!# 668 "esm_stres_mod.f90"
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 676 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 686 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 689 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 696 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 704 "esm_stres_mod.f90"
      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 711 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz
!# 715 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 720 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END IF ! imill_2d(0,0) > 0
!# 728 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)
!# 732 "esm_stres_mod.f90"
    CALL mp_sum(sigmahar, intra_bgrp_comm)
!# 734 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 740 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_har_bc2
!# 743 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_har_bc3(sigmahar, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 755 "esm_stres_mod.f90"
    REAL(DP), INTENT(out)   :: sigmahar(3, 3)
    COMPLEX(DP), INTENT(in) :: rhog(ngm)   !  n(G)
!# 758 "esm_stres_mod.f90"
    INTEGER :: ig, iga, igb, igz, igp, la, mu, iz, jz
    REAL(DP) :: L, S, z0, z1, z
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: sum1p, sum1m, sum2p, sum2m, sum1c, sum2c
    COMPLEX(DP) :: sum1sh, sum1ch, sum2sh
    REAL(DP)    :: z_l, z_r
    COMPLEX(DP) :: f1, f2, f3, f4, a0, a1, a2, a3
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dgp2_deps(2, 2)  !! dgp^2/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 774 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 780 "esm_stres_mod.f90"
    REAL(DP) :: sigmahar_bc1(3, 3)
!# 782 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
!# 788 "esm_stres_mod.f90"
    ! initialize
    sigmahar(:, :) = 0.0d0
!# 791 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 797 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 800 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 809 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 812 "esm_stres_mod.f90"
      ! expand function symmetrically to gz<0
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      END IF
    END DO ! ig
!# 822 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 829 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 831 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dgp2_deps(la, mu) = -g(la)*g(mu)*2.0d0
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 840 "esm_stres_mod.f90"
      ! summations over gz
      sum1p = (0.d0, 0.d0)
      sum1m = (0.d0, 0.d0)
      sum2p = (0.d0, 0.d0)
      sum2m = (0.d0, 0.d0)
      sum1sh = (0.d0, 0.d0)
      sum1ch = (0.d0, 0.d0)
      sum2sh = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 853 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
        sum1p = sum1p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)
        sum1m = sum1m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)
        sum2p = sum2p + rg3*QE_EXP(+ci*gz*z0)/(gp - ci*gz)**2
        sum2m = sum2m + rg3*QE_EXP(-ci*gz*z0)/(gp + ci*gz)**2
        sum1sh = sum1sh + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)
        sum1ch = sum1ch + rg3*QE_COSH(gp*z0 + ci*gz*z0)/(gp + ci*gz)*z0
        sum2sh = sum2sh + rg3*QE_SINH(gp*z0 + ci*gz*z0)/(gp + ci*gz)**2
      END DO ! igz
!# 863 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 871 "esm_stres_mod.f90"
        !! BC1 terms
        dVr_deps(iz, :, :) = &
          -(dgp_deps(:, :)*tpi/gp**2*(gp*(z - z0) - 1.0d0) &
            - delta(:, :)*tpi/gp) &
          *EXP(+gp*(z - z0))*sum1p &
          + dgp_deps(:, :)*tpi/gp*EXP(+gp*(z - z0))*sum2p &
          + (dgp_deps(:, :)*tpi/gp**2*(gp*(z + z0) + 1.0d0) &
             + delta(:, :)*tpi/gp) &
          *EXP(-gp*(z + z0))*sum1m &
          + dgp_deps(:, :)*tpi/gp*EXP(-gp*(z + z0))*sum2m
!# 882 "esm_stres_mod.f90"
        !! BC3 termn
        dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                             - dgp_deps(:, :)*( &
                             -fpi/gp**2*EXP(-gp*(-z + 2*z1)) &
                             - fpi/gp*(-z + 2*z1)*EXP(-gp*(-z + 2*z1)) &
                             )*sum1sh &
                             - fpi/gp*EXP(-gp*(-z + 2*z1))*( &
                             -delta(:, :)*sum1sh &
                             + dgp_deps(:, :)*(sum1ch - sum2sh))
      END DO ! iz
!# 893 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 900 "esm_stres_mod.f90"
      ! add bare coulomn terms to dV(gz)/deps
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 906 "esm_stres_mod.f90"
        rg3 = rhog3(iz, igp)
!# 908 "esm_stres_mod.f90"
        dVg_deps(iz, :, :) = dVg_deps(iz, :, :) &
                              - delta(:, :)*fpi*rg3/(gp**2 + gz**2) &
                              - dgp2_deps(:, :)*fpi*rg3/(gp**2 + gz**2)**2
      END DO ! igz
!# 913 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 918 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, :, :))
      END DO ! igz
    END DO ! igp
!# 926 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! summations over gz
      sum1c = (0.d0, 0.d0)
      sum2c = (0.d0, 0.d0)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 937 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        sum1c = sum1c + rg3*ci*cos(gz*z0)/gz
        sum2c = sum2c + rg3*cos(gz*z0)/gz**2
      END DO ! igz
!# 942 "esm_stres_mod.f90"
      ! calculate V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 950 "esm_stres_mod.f90"
        rg3 = rhog3(1, imill_2d(0, 0))
        !! BC1 terms
        Vr(iz) = &
          - tpi*z**2*rg3 &
          - tpi*z0**2*rg3 &
          - fpi*z*sum1c &
          - fpi*sum2c
!# 958 "esm_stres_mod.f90"
        !! BC3 terms
        Vr(iz) = Vr(iz) - tpi*(z - 2*z1)*2*z0*rg3 + fpi*z0*sum1c
      END DO ! iz
!# 962 "esm_stres_mod.f90"
      ! separation by polynomial
      z_l = -z0
      z_r = +z0
      f1 = -tpi*z_r**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_r*sum1c &
           - fpi*sum2c
      f1 = f1 &
           - tpi*(z_r - 2*z1)*2*z0*rg3 + fpi*z0*sum1c
!# 972 "esm_stres_mod.f90"
      f2 = -tpi*z_l**2*rg3 &
           - tpi*z0**2*rg3 &
           - fpi*z_l*sum1c &
           - fpi*sum2c
      f2 = f2 &
           - tpi*(z_l - 2*z1)*2*z0*rg3 + fpi*z0*sum1c
!# 979 "esm_stres_mod.f90"
      f3 = -fpi*z_r*rg3 &
           - fpi*sum1c
      f3 = f3 &
           - tpi*(1.0d0)*2*z0*rg3
!# 984 "esm_stres_mod.f90"
      f4 = -fpi*z_l*rg3 &
           - fpi*sum1c
      f4 = f4 &
           - tpi*(1.0d0)*2*z0*rg3
!# 989 "esm_stres_mod.f90"
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 997 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 1007 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 1010 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 1017 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 1025 "esm_stres_mod.f90"
      ! add bare coulomn terms to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 1032 "esm_stres_mod.f90"
        rg3 = rhog3(iz, imill_2d(0, 0))
        Vg(iz) = Vg(iz) + fpi*rg3/gz**2
      END DO ! igz
!# 1036 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 1041 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
!# 1045 "esm_stres_mod.f90"
        sigmahar(1:2, 1:2) = sigmahar(1:2, 1:2) &
                             + REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END IF ! imill_2d(0,0) > 0
!# 1050 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmahar(:, :) = sigmahar(:, :)*(-0.5d0*e2)
!# 1054 "esm_stres_mod.f90"
    CALL mp_sum(sigmahar, intra_bgrp_comm)
!# 1056 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 1062 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_har_bc3
!# 1065 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewr(alpha, sigmaewa)
    USE kinds, ONLY : DP
    IMPLICIT NONE
!# 1069 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1072 "esm_stres_mod.f90"
    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc2')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc3')
      CALL esm_stres_ewr_pbc(alpha, sigmaewa)
    CASE ('bc4')
      STOP 'esm_stres_ewa has not yet implemented for esm_bc = bc4'
    END SELECT
!# 1085 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewr
!# 1088 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewr_pbc(alpha, sigmaewa)
    USE kinds,     ONLY : DP
    USE constants, ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base, ONLY : omega, alat, at, tpiba, bg
    USE ions_base, ONLY : zv, nat, tau, ityp
    USE gvect,     ONLY : gstart
    USE mp_bands,  ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    IMPLICIT NONE
!# 1098 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1101 "esm_stres_mod.f90"
    INTEGER, PARAMETER :: mxr = 50
    ! the maximum number of R vectors included in r sum
    INTEGER  :: ia, ib, nr, nrm, la, mu
    REAL(DP) :: Qa, Qb, dtau(3), rmax
    REAL(DP) :: salp, r(3, mxr), r2(mxr), rr, fac
!# 1107 "esm_stres_mod.f90"
    salp = sqrt(alpha)
!# 1109 "esm_stres_mod.f90"
    ! initialize
    sigmaewa(:, :) = 0.d0
!# 1112 "esm_stres_mod.f90"
    !
    ! R-space sum here (only for the processor that contains G=0)
    !
    IF (gstart == 2) THEN
      rmax = 4.0d0/salp/alat
      !
      ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-1
      !
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          !
          !     generates nearest-neighbors shells r(i)=R(i)-dtau(i)
          !
          dtau(:) = tau(:, ib) - tau(:, ia)
          CALL rgen(dtau, rmax, mxr, at, bg, r, r2, nrm)
!# 1130 "esm_stres_mod.f90"
          DO nr = 1, nrm
            rr = sqrt(r2(nr))*alat
            r(:, nr) = r(:, nr)*alat
!# 1134 "esm_stres_mod.f90"
            fac = Qb*Qa/rr**3 &
                  *(erfc(salp*rr) &
                  + rr*2.0d0*salp*sqrtpm1*EXP(-alpha*rr**2))
            DO la = 1, 3
              DO mu = 1, 3
                sigmaewa(la, mu) = sigmaewa(la, mu) + fac*r(la, nr)*r(mu, nr)
              END DO ! mu
            END DO ! la
          END DO ! nr
        END DO ! ia
      END DO ! ib
    END IF
!# 1147 "esm_stres_mod.f90"
    sigmaewa(:, :) = sigmaewa(:, :)*(e2/2.0d0/omega)
!# 1149 "esm_stres_mod.f90"
    CALL mp_sum(sigmaewa, intra_bgrp_comm)
!# 1151 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewr_pbc
!# 1154 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewg(alpha, sigmaewa)
    USE kinds, ONLY : DP
    IMPLICIT NONE
!# 1158 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1161 "esm_stres_mod.f90"
    SELECT CASE (esm_bc)
    CASE ('pbc')
      STOP 'esm_stres_ewa must not be called for esm_bc = pbc'
    CASE ('bc1')
      CALL esm_stres_ewg_bc1(alpha, sigmaewa)
    CASE ('bc2')
      CALL esm_stres_ewg_bc2(alpha, sigmaewa)
    CASE ('bc3')
      CALL esm_stres_ewg_bc3(alpha, sigmaewa)
    CASE ('bc4')
      STOP 'esm_stres_ewa must not be called for esm_bc = bc4'
    END SELECT
!# 1174 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewg
!# 1177 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewg_bc1(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 1188 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1191 "esm_stres_mod.f90"
    INTEGER  :: ia, ib, igp, iga, igb, la, mu, iz
    REAL(DP) :: L, S, salp
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact
!# 1198 "esm_stres_mod.f90"
    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 1203 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    salp  = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp
!# 1210 "esm_stres_mod.f90"
    ! initialize
    sigmaewa(:, :) = 0.0d0
!# 1213 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF
!# 1222 "esm_stres_mod.f90"
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 1230 "esm_stres_mod.f90"
        ! distance between atoms is defined
        zbza = zb - za
!# 1233 "esm_stres_mod.f90"
        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 1241 "esm_stres_mod.f90"
          IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 1243 "esm_stres_mod.f90"
          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza
!# 1249 "esm_stres_mod.f90"
          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO
!# 1257 "esm_stres_mod.f90"
          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          !
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp
        END DO ! igp
!# 1286 "esm_stres_mod.f90"
        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF
!# 1291 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
!# 1294 "esm_stres_mod.f90"
      END DO ! ia
    END DO ! ib
!# 1297 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF
!# 1307 "esm_stres_mod.f90"
        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
!# 1316 "esm_stres_mod.f90"
          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
                 + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)
        END DO ! ia
!# 1321 "esm_stres_mod.f90"
        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb
!# 1323 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF
!# 1328 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)
!# 1332 "esm_stres_mod.f90"
    CALL mp_sum(sigmaewa, intra_bgrp_comm)
!# 1334 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewg_bc1
!# 1337 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewg_bc2(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 1348 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1351 "esm_stres_mod.f90"
    INTEGER  :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP) :: exph1, exph2, exph3
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact
!# 1359 "esm_stres_mod.f90"
    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 1364 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    salp = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp
!# 1373 "esm_stres_mod.f90"
    ! initialize
    sigmaewa(:, :) = 0.0d0
!# 1376 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF
!# 1385 "esm_stres_mod.f90"
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 1393 "esm_stres_mod.f90"
        ! distance between atoms is defined
        zbza = zb - za
!# 1396 "esm_stres_mod.f90"
        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 1404 "esm_stres_mod.f90"
          IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 1406 "esm_stres_mod.f90"
          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza
!# 1412 "esm_stres_mod.f90"
          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO
!# 1420 "esm_stres_mod.f90"
          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          !
          exph1 = (cosh(gp*(zb - za))*EXP(-2*gp*z1) - cosh(gp*(zb + za)))/sinh(2*gp*z1)
          exph2 = ((zb - za)*sinh(gp*(zb - za))*EXP(-2*gp*z1) &
                  - 2*z1*cosh(gp*(zb - za))*EXP(-2*gp*z1) &
                  - (zb + za)*sinh(gp*(zb + za)))/sinh(2*gp*z1)
          exph3 = -(cosh(gp*(zb - za))*EXP(-2*gp*z1) - cosh(gp*(zb + za)))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)
!# 1446 "esm_stres_mod.f90"
          !! BC1 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp
!# 1455 "esm_stres_mod.f90"
          !! BC2 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*tpi/gp*Qb*Qa/S*cosgpr*exph1 &
                          - tpi/gp*delta(:, :)*Qb*Qa/S*cosgpr*exph1 &
                          + tpi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*(exph2 + exph3)
        END DO ! igp
!# 1462 "esm_stres_mod.f90"
        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF
!# 1467 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
!# 1470 "esm_stres_mod.f90"
      END DO ! ia
    END DO ! ib
!# 1473 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF
!# 1483 "esm_stres_mod.f90"
        ! [note] this Vr does not contain a term due to efield z*efield
        ! because it vanishes in the differentiation with respect to strain.
        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
!# 1494 "esm_stres_mod.f90"
          !! BC1 terms
          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
                 + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)
!# 1499 "esm_stres_mod.f90"
          !! BC2 terms
          Vr = Vr + tpi*Qa/S*(-zb*za + z1*z1)/z1
        END DO ! ia
!# 1503 "esm_stres_mod.f90"
        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb
!# 1505 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF
!# 1510 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)
!# 1514 "esm_stres_mod.f90"
    CALL mp_sum(sigmaewa, intra_bgrp_comm)
!# 1516 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewg_bc2
!# 1519 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_ewg_bc3(alpha, sigmaewa)
    USE kinds,         ONLY : DP
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE gvect,         ONLY : gstart
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 1530 "esm_stres_mod.f90"
    REAL(DP), INTENT(in)  :: alpha
    REAL(DP), INTENT(out) :: sigmaewa(3, 3)
!# 1533 "esm_stres_mod.f90"
    INTEGER  :: ia, ib, igp, iga, igb, la, mu
    REAL(DP) :: L, S, salp, z0, z1
    REAL(DP) :: Qa, Qb, ra(2), rb(2), za, zb
    REAL(DP) :: g(2), gp, Vr
    REAL(DP) :: cosgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp, expm
    REAL(DP) :: zbza, isalp, gpzbza, gp2a, mgazz, pgazz, fact
!# 1540 "esm_stres_mod.f90"
    REAL(DP) :: dE_deps(2, 2)
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
!# 1545 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    salp = sqrt(alpha)
    isalp = 1.0_DP/salp
    fact  = 0.5_DP * isalp
!# 1554 "esm_stres_mod.f90"
    ! initialize
    sigmaewa(:, :) = 0.0d0
!# 1557 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO ib = 1, nat
      Qb = (-1.0d0)*zv(ityp(ib))
      rb(1:2) = tau(1:2, ib)*alat
      zb = tau(3, ib)*alat
      IF (zb > L*0.5d0) THEN
        zb = zb - L
      END IF
!# 1566 "esm_stres_mod.f90"
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 1574 "esm_stres_mod.f90"
        ! distance between atoms is defined
        zbza = zb - za
!# 1577 "esm_stres_mod.f90"
        ! summations over gp
        dE_deps(:, :) = 0.0d0
        DO igp = 1, ngm_2d
          iga = mill_2d(1, igp)
          igb = mill_2d(2, igp)
          g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
          gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 1585 "esm_stres_mod.f90"
          IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 1587 "esm_stres_mod.f90"
          ! define exp phases
          gpzbza = gp * zbza
          gp2a   = gp * 0.5_DP * isalp
          mgazz  = gp2a - salp * zbza
          pgazz  = gp2a + salp * zbza
!# 1593 "esm_stres_mod.f90"
          ! derivatives by strain tensor
          DO la = 1, 2
            DO mu = 1, 2
              dgp_deps(la, mu) = -g(la)*g(mu)/gp
              dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
            END DO
          END DO
!# 1601 "esm_stres_mod.f90"
          ! coefficients
          cosgpr = cos(g(1)*(rb(1) - ra(1)) + g(2)*(rb(2) - ra(2)))
          experfcm = exp_erfc(-gpzbza, mgazz)
          experfcp = exp_erfc(+gpzbza, pgazz)
          dexperfcm_dgp = -zbza*experfcm &
                          -exp_gauss( -gpzbza, mgazz ) * fact
          dexperfcp_dgp = +zbza*experfcp &
                          -exp_gauss( +gpzbza, pgazz ) * fact
          !
          ! Old code is not safe, because diverged terms are included.
          ! However, this code is a faithful for original formula.
          ! For this reason, we leave following old codes as comment.
          !
          ! experfcm = exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za))
          ! experfcp = exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za))
          ! dexperfcm_dgp = -(zb - za)*exp_erfc(-gp*(zb - za), gp/2.d0/salp - salp*(zb - za)) &
          !                 - EXP(-gp*(zb - za))*qe_gauss(gp/2.d0/salp - salp*(zb - za))/2.d0/salp
          ! dexperfcp_dgp = +(zb - za)*exp_erfc(+gp*(zb - za), gp/2.d0/salp + salp*(zb - za)) &
          !                 - EXP(+gp*(zb - za))*qe_gauss(gp/2.d0/salp + salp*(zb - za))/2.d0/salp
          expm = EXP(-gp*(-zb + 2*z1 - za))
          !
          !! BC1 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcm &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcm &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcm_dgp &
                          + gp*dinvgp_deps(:, :)*pi/gp*Qb*Qa/S*cosgpr*experfcp &
                          - pi/gp*delta(:, :)*Qb*Qa/S*cosgpr*experfcp &
                          + pi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*dexperfcp_dgp
!# 1631 "esm_stres_mod.f90"
          !! BC3 terms
          dE_deps(:, :) = dE_deps(:, :) &
                          - gp*dinvgp_deps(:, :)*tpi/gp*Qb*Qa/S*cosgpr*expm &
                          + tpi/gp*delta(:, :)*Qb*Qa/S*cosgpr*expm &
                          + tpi/gp*Qb*Qa/S*cosgpr*dgp_deps(:, :)*(-zb + 2*z1 - za)*expm
        END DO ! igp
!# 1638 "esm_stres_mod.f90"
        ! modifications
        IF (gamma_only) THEN
          dE_deps(:, :) = dE_deps(:, :)*2.0d0
        END IF
!# 1643 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
!# 1646 "esm_stres_mod.f90"
      END DO ! ia
    END DO ! ib
!# 1649 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (gstart == 2) THEN
      DO ib = 1, nat
        Qb = (-1.0d0)*zv(ityp(ib))
        rb(1:2) = tau(1:2, ib)*alat
        zb = tau(3, ib)*alat
        IF (zb > L*0.5d0) THEN
          zb = zb - L
        END IF
!# 1659 "esm_stres_mod.f90"
        Vr = 0.0d0
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
!# 1668 "esm_stres_mod.f90"
          !! BC1 terms
          Vr = Vr - tpi*Qa/S &
               *((zb - za)*erf(salp*(zb - za)) &
               + EXP(-alpha*(zb - za)**2)*sqrtpm1/salp)
!# 1673 "esm_stres_mod.f90"
          !! BC3 terms
          Vr = Vr + tpi*Qa/S*(-zb + 2*z1 - za)
        END DO ! ia
!# 1677 "esm_stres_mod.f90"
        dE_deps(1:2, 1:2) = -delta(1:2, 1:2)*Vr*Qb
!# 1679 "esm_stres_mod.f90"
        ! calculate stress tensor
        sigmaewa(1:2, 1:2) = sigmaewa(1:2, 1:2) - dE_deps(1:2, 1:2)/omega
      END DO ! ib
    END IF
!# 1684 "esm_stres_mod.f90"
    ! half means removing duplications.
    ! e2 means hartree -> Ry.
    sigmaewa(:, :) = sigmaewa(:, :)*(0.5d0*e2)
!# 1688 "esm_stres_mod.f90"
    CALL mp_sum(sigmaewa, intra_bgrp_comm)
!# 1690 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_ewg_bc3
!# 1693 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_loclong_bc1(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 1706 "esm_stres_mod.f90"
    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)
!# 1709 "esm_stres_mod.f90"
    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    COMPLEX(DP) :: poly_fr, poly_fl, poly_dfr, poly_dfl
    COMPLEX(DP) :: poly_a, poly_b, poly_c, poly_d
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza
!# 1726 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 1732 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 1738 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 1741 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 1750 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 1753 "esm_stres_mod.f90"
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig
!# 1762 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S     = omega/L
    z0    = L/2.d0
    alpha = 1.0d0
    salp  = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    ! initialize
    sigmaloclong(:, :) = 0.0d0
!# 1774 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 1781 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 1783 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 1791 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 1799 "esm_stres_mod.f90"
        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp
!# 1840 "esm_stres_mod.f90"
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp
        END DO ! ia
      END DO ! iz
!# 1847 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 1854 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 1859 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp
!# 1867 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 1883 "esm_stres_mod.f90"
        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 1890 "esm_stres_mod.f90"
          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                     + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)
        END DO ! iz
!# 1895 "esm_stres_mod.f90"
        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
               + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
               + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
      END DO ! ia
!# 1907 "esm_stres_mod.f90"
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
                                              + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
            - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
                                                                   + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 1915 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 1925 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 1928 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 1935 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 1943 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 1948 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0
!# 1956 "esm_stres_mod.f90"
    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)
!# 1959 "esm_stres_mod.f90"
    CALL mp_sum(sigmaloclong, intra_bgrp_comm)
!# 1961 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 1967 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_loclong_bc1
!# 1970 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_loclong_bc2(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 1983 "esm_stres_mod.f90"
    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)
!# 1986 "esm_stres_mod.f90"
    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: exph1, exph2, exph3
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza
!# 2002 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 2008 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 2014 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 2017 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 2026 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 2029 "esm_stres_mod.f90"
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig
!# 2038 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    ! initialize
    sigmaloclong(:, :) = 0.0d0
!# 2051 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 2058 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 2060 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 2068 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 2076 "esm_stres_mod.f90"
        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          exph1 = (cosh(gp*(z - za))*EXP(-2*gp*z1) - cosh(gp*(z + za)))/sinh(2*gp*z1)
          exph2 = ((z - za)*sinh(gp*(z - za))*EXP(-2*gp*z1) &
                   - 2*z1*cosh(gp*(z - za))*EXP(-2*gp*z1) &
                   - (z + za)*sinh(gp*(z + za)))/sinh(2*gp*z1)
          exph3 = -(cosh(gp*(z - za))*EXP(-2*gp*z1) - cosh(gp*(z + za)))/sinh(2*gp*z1)**2*2*z1*cosh(2*gp*z1)
!# 2118 "esm_stres_mod.f90"
          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp
!# 2124 "esm_stres_mod.f90"
          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp
!# 2130 "esm_stres_mod.f90"
          !! BC2 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*tpi/gp*Qa/S*expimgpr*exph1 &
                               - tpi/gp*delta(:, :)*Qa/S*expimgpr*exph1 &
                               + tpi/gp*Qa/S*expimgpr*dgp_deps(:, :)*(exph2 + exph3)
        END DO ! ia
      END DO ! iz
!# 2138 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 2145 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 2150 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp
!# 2158 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      ! [note] this Vr does not contain a term due to efield z*efield
      ! because it vanishes in the differentiation with respect to strain.
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 2176 "esm_stres_mod.f90"
        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 2183 "esm_stres_mod.f90"
          !! BC1 terms
          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                     + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)
!# 2188 "esm_stres_mod.f90"
          !! BC2 terms
          Vr(iz) = Vr(iz) + tpi*Qa/S*(-z*za + z1*z1)/z1
        END DO ! iz
!# 2192 "esm_stres_mod.f90"
        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
               + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f1 = f1 + tpi*Qa/S*(-z_r*za + z1*z1)/z1
!# 2197 "esm_stres_mod.f90"
        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
               + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f2 = f2 + tpi*Qa/S*(-z_l*za + z1*z1)/z1
!# 2202 "esm_stres_mod.f90"
        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f3 = f3 + tpi*Qa/S*(-za)/z1
!# 2206 "esm_stres_mod.f90"
        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
        f4 = f4 + tpi*Qa/S*(-za)/z1
!# 2210 "esm_stres_mod.f90"
      END DO ! ia
!# 2212 "esm_stres_mod.f90"
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 2220 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 2230 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 2233 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 2240 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 2248 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 2253 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0
!# 2261 "esm_stres_mod.f90"
    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)
!# 2264 "esm_stres_mod.f90"
    CALL mp_sum(sigmaloclong, intra_bgrp_comm)
!# 2266 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 2272 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_loclong_bc2
!# 2275 "esm_stres_mod.f90"
  SUBROUTINE esm_stres_loclong_bc3(sigmaloclong, rhog)
    USE kinds,         ONLY : DP
    USE gvect,         ONLY : ngm, mill
    USE constants,     ONLY : pi, sqrtpm1, tpi, fpi, e2
    USE cell_base,     ONLY : omega, alat, at, tpiba, bg
    USE ions_base,     ONLY : zv, nat, tau, ityp
    USE control_flags, ONLY : gamma_only
    USE fft_base,      ONLY : dfftp
    USE fft_scalar,    ONLY : cft_1z
    USE mp_bands,      ONLY : intra_bgrp_comm
    USE mp,            ONLY : mp_sum
    IMPLICIT NONE
!# 2288 "esm_stres_mod.f90"
    REAL(DP), INTENT(out) :: sigmaloclong(3, 3)
    COMPLEX(DP) :: rhog(ngm)   !  n(G)
!# 2291 "esm_stres_mod.f90"
    INTEGER  :: ig, iga, igb, igz, igp, la, mu, iz, jz, ia
    REAL(DP) :: L, S, z0, z1, alpha, salp, z
    REAL(DP) :: Qa, ra(2), za
    REAL(DP) :: g(2), gp, gz
    COMPLEX(DP), PARAMETER :: ci = dcmplx(0.0d0, 1.0d0)
    COMPLEX(DP) :: rg3
    COMPLEX(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp
    COMPLEX(DP) :: expm
    REAL(DP)    :: z_r, z_l
    COMPLEX(DP) :: a0, a1, a2, a3, f1, f2, f3, f4
    REAL(DP), PARAMETER :: delta(2, 2) = reshape((/1.0d0, 0.0d0, 0.0d0, 1.0d0/), (/2, 2/))
    REAL(DP) :: dgp_deps(2, 2)  !! dgp/deps
    REAL(DP) :: dinvgp_deps(2, 2)  !! dgp^-1/deps
    REAL(DP) :: isalp, fact
    REAL(DP) :: zza, mgza, pgza, g2a_maza, g2a_paza, tz1
!# 2307 "esm_stres_mod.f90"
    COMPLEX(DP), ALLOCATABLE :: rhog3(:, :)
    COMPLEX(DP), ALLOCATABLE :: dVr_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: dVg_deps(:, :, :)
    COMPLEX(DP), ALLOCATABLE :: Vr(:)
    COMPLEX(DP), ALLOCATABLE :: Vg(:)
!# 2313 "esm_stres_mod.f90"
    REAL(DP) :: sigmaloclong_bc1(3, 3)
!# 2315 "esm_stres_mod.f90"
    ALLOCATE (rhog3(dfftp%nr3, ngm_2d))
    ALLOCATE (dVr_deps(dfftp%nr3, 2, 2))
    ALLOCATE (dVg_deps(dfftp%nr3, 2, 2))
    ALLOCATE (Vr(dfftp%nr3))
    ALLOCATE (Vg(dfftp%nr3))
!# 2321 "esm_stres_mod.f90"
    ! reconstruct rho(gz,gp)
    rhog3(:, :) = (0.d0, 0.d0)
!# 2324 "esm_stres_mod.f90"
    DO ig = 1, ngm
      iga = mill(1, ig)
      igb = mill(2, ig)
      igz = mill(3, ig) + 1
      igp = imill_2d(iga, igb)
      IF (igz < 1) THEN
        igz = igz + dfftp%nr3
      END IF
!# 2333 "esm_stres_mod.f90"
      rg3 = rhog(ig)
      rhog3(igz, igp) = rg3
!# 2336 "esm_stres_mod.f90"
      IF (gamma_only .and. iga == 0 .and. igb == 0) THEN
        igz = 1 - mill(3, ig)
        IF (igz < 1) THEN
          igz = igz + dfftp%nr3
        END IF
        rhog3(igz, igp) = CONJG(rg3)
      ENDIF
    END DO ! ig
!# 2345 "esm_stres_mod.f90"
    ! cell settings
    L = at(3, 3)*alat
    S = omega/L
    z0 = L/2.d0
    z1 = z0 + esm_w
    alpha = 1.0d0
    salp = sqrt(alpha)
    ! useful values are setted
    isalp = 1.d0/sqrt(alpha)
    fact  = 1.d0/2.d0/salp
    tz1 = 2.d0 * z1
    !
    ! initialize
    sigmaloclong(:, :) = 0.0d0
!# 2360 "esm_stres_mod.f90"
    !****For gp!=0 case ********************
    DO igp = 1, ngm_2d
      iga = mill_2d(1, igp)
      igb = mill_2d(2, igp)
      g(1:2) = (iga*bg(1:2, 1) + igb*bg(1:2, 2))*tpiba
      gp = sqrt(g(1)*g(1) + g(2)*g(2))
!# 2367 "esm_stres_mod.f90"
      IF (gp == 0.0d0) CYCLE ! skip gp=0
!# 2369 "esm_stres_mod.f90"
      ! derivatives by strain tensor
      DO la = 1, 2
        DO mu = 1, 2
          dgp_deps(la, mu) = -g(la)*g(mu)/gp
          dinvgp_deps(la, mu) = +g(la)*g(mu)/gp**3
        END DO
      END DO
!# 2377 "esm_stres_mod.f90"
      ! calculate dV(z)/deps
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 2385 "esm_stres_mod.f90"
        ! summations over all atoms
        dVr_deps(iz, :, :) = (0.0d0, 0.0d0)
        DO ia = 1, nat
          Qa = (-1.0d0)*zv(ityp(ia))
          ra(1:2) = tau(1:2, ia)*alat
          za = tau(3, ia)*alat
          IF (za > L*0.5d0) THEN
            za = za - L
          END IF
! --------------------------------------------------------------------------------------------------
!         Following code is old version. However, this code explicitly shows the formulation of 
!         stress tensor within ESM scheme. For this reason, we left the code as comment.
!          expimgpr = qe_exp(-ci*(g(1)*ra(1) + g(2)*ra(2)))
!          experfcm = exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za))
!          experfcp = exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za))
!          dexperfcm_dgp = -(z - za)*exp_erfc(-gp*(z - za), gp/2.d0/salp - salp*(z - za)) &
!                          - EXP(-gp*(z - za))*qe_gauss(gp/2.d0/salp - salp*(z - za))/2.d0/salp
!          dexperfcp_dgp = +(z - za)*exp_erfc(+gp*(z - za), gp/2.d0/salp + salp*(z - za)) &
!                          - EXP(+gp*(z - za))*qe_gauss(gp/2.d0/salp + salp*(z - za))/2.d0/salp
!---------------------------------------------------------------------------------------------------
          !
          ! ... Set useful values
          zza = z - za
          mgza = -gp*zza
          pgza = +gp*zza
          g2a_maza = gp*0.5d0*isalp - salp * zza
          g2a_paza = gp*0.5d0*isalp + salp * zza
          !
          expimgpr = QE_EXP(-ci*(g(1)*ra(1) + g(2)*ra(2)))
          experfcm = exp_erfc(mgza, g2a_maza)
          experfcp = exp_erfc(pgza, g2a_paza)
          dexperfcm_dgp = -zza * experfcm &
                          -exp_gauss( mgza, g2a_maza ) * fact
          dexperfcp_dgp = +zza * experfcp &
                          -exp_gauss( pgza, g2a_paza ) * fact
          !
          expm = EXP(-gp*(-z + tz1 - za))
!# 2423 "esm_stres_mod.f90"
          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcm &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcm &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcm_dgp
!# 2429 "esm_stres_mod.f90"
          !! BC1 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               + gp*dinvgp_deps(:, :)*pi/gp*Qa/S*expimgpr*experfcp &
                               - pi/gp*delta(:, :)*Qa/S*expimgpr*experfcp &
                               + pi/gp*Qa/S*expimgpr*dgp_deps(:, :)*dexperfcp_dgp
!# 2435 "esm_stres_mod.f90"
          !! BC3 terms
          dVr_deps(iz, :, :) = dVr_deps(iz, :, :) &
                               - gp*dinvgp_deps(:, :)*tpi/gp*Qa/S*expimgpr*expm &
                               + tpi/gp*delta(:, :)*Qa/S*expimgpr*expm &
                               + tpi/gp*Qa/S*expimgpr*dgp_deps(:, :)*(-z + 2*z1 - za)*expm
        END DO ! ia
      END DO ! iz
!# 2443 "esm_stres_mod.f90"
      ! convert dV(z)/deps to dV(gz)/deps
      DO la = 1, 2
        DO mu = 1, 2
          CALL cft_1z(dVr_deps(:, la, mu), 1, dfftp%nr3, dfftp%nr3, -1, dVg_deps(:, la, mu))
        END DO
      END DO
!# 2450 "esm_stres_mod.f90"
      ! modifications
      IF (gamma_only) THEN
        dVg_deps(:, :, :) = dVg_deps(:, :, :)*2.0d0
      END IF
!# 2455 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, igp)
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    END DO ! igp
!# 2463 "esm_stres_mod.f90"
    !****For gp=0 case ********************
    IF (imill_2d(0, 0) > 0) THEN
      ! calculate V(z)
      Vr(:) = 0.0d0
      ! separation by polynomial
      f1 = (0.d0, 0.d0); f2 = (0.d0, 0.d0); f3 = (0.d0, 0.d0); f4 = (0.d0, 0.d0)
      z_l = -z0
      z_r = +z0
      DO ia = 1, nat
        Qa = (-1.0d0)*zv(ityp(ia))
        ra(1:2) = tau(1:2, ia)*alat
        za = tau(3, ia)*alat
        IF (za > L*0.5d0) THEN
          za = za - L
        END IF
!# 2479 "esm_stres_mod.f90"
        DO iz = 1, dfftp%nr3
          jz = iz - 1
          IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
            jz = jz - dfftp%nr3
          END IF
          z = DBLE(jz) / DBLE(dfftp%nr3) * L
!# 2486 "esm_stres_mod.f90"
          !! BC1 terms
          Vr(iz) = Vr(iz) - tpi*Qa/S &
                   *((z - za)*erf(salp*(z - za)) &
                   + EXP(-alpha*(z - za)**2)*sqrtpm1/salp)
!# 2491 "esm_stres_mod.f90"
          !! BC3 terms
          Vr(iz) = Vr(iz) + tpi*Qa/S*(-z + 2*z1 - za)
        END DO ! iz
!# 2495 "esm_stres_mod.f90"
        f1 = f1 - tpi*Qa/S &
             *((z_r - za)*erf(salp*(z_r - za)) &
             + EXP(-alpha*(z_r - za)**2)*sqrtpm1/salp)
        f1 = f1 + tpi*Qa/S*(-z_r + 2*z1 - za)
!# 2500 "esm_stres_mod.f90"
        f2 = f2 - tpi*Qa/S &
             *((z_l - za)*erf(salp*(z_l - za)) &
             + EXP(-alpha*(z_l - za)**2)*sqrtpm1/salp)
        f2 = f2 + tpi*Qa/S*(-z_l + 2*z1 - za)
!# 2505 "esm_stres_mod.f90"
        f3 = f3 - tpi*Qa/S &
             *erf(salp*(z_r - za))
        f3 = f3 + tpi*Qa/S*(-1.0d0)
!# 2509 "esm_stres_mod.f90"
        f4 = f4 - tpi*Qa/S &
             *erf(salp*(z_l - za))
        f4 = f4 + tpi*Qa/S*(-1.0d0)
      END DO ! ia
!# 2514 "esm_stres_mod.f90"
      a0 = (f1*z_l**2*(z_l - 3.d0*z_r) + z_r*(f3*z_l**2*(-z_l + z_r) &
           + z_r*(f2*(3.d0*z_l - z_r) + f4*z_l*(-z_l + z_r))))/(z_l - z_r)**3
      a1 = (f3*z_l**3 + z_l*(6.d0*f1 - 6.d0*f2 + (f3 + 2.d0*f4)*z_l)*z_r &
           - (2*f3 + f4)*z_l*z_r**2 - f4*z_r**3)/(z_l - z_r)**3
      a2 = (-3*f1*(z_l + z_r) + 3.d0*f2*(z_l + z_r) - (z_l - z_r)*(2*f3*z_l &
           + f4*z_l + f3*z_r + 2*f4*z_r))/(z_l - z_r)**3
      a3 = (2.d0*f1 - 2.d0*f2 + (f3 + f4)*(z_l - z_r))/(z_l - z_r)**3
!# 2522 "esm_stres_mod.f90"
      ! remove polynomial from V(z)
      DO iz = 1, dfftp%nr3
        jz = iz - 1
        IF (jz >= (dfftp%nr3 - dfftp%nr3/2)) THEN
          jz = jz - dfftp%nr3
        END IF
        z = DBLE(jz) / DBLE(dfftp%nr3) * L
        Vr(iz) = Vr(iz) - (a0 + a1*z + a2*z**2 + a3*z**3)
      ENDDO
!# 2532 "esm_stres_mod.f90"
      ! convert V(z) to V(gz) without polynomial
      CALL cft_1z(Vr, 1, dfftp%nr3, dfftp%nr3, -1, Vg)
!# 2535 "esm_stres_mod.f90"
      ! add polynomial to V(gz)
      DO igz = -(dfftp%nr3 - 1)/2, (dfftp%nr3 - 1)/2
        IF (igz == 0) CYCLE
        iz = igz + 1
        IF (iz < 1) iz = iz + dfftp%nr3
        gz = dble(igz)*tpi/L
!# 2542 "esm_stres_mod.f90"
        Vg(iz) = Vg(iz) &
                  + a1*ci*cos(gz*z0)/gz &
                  + a2*2.0d0*cos(gz*z0)/gz**2 &
                  + a3*ci*z0**2*cos(gz*z0)/gz &
                  - a3*ci*6.0d0*cos(gz*z0)/gz**3
      END DO
      Vg(1) = Vg(1) + a0*1.0d0 + a2*z0**2/3.0d0
!# 2550 "esm_stres_mod.f90"
      ! calculate dV/deps(gz)
      DO igz = 1, dfftp%nr3
        dVg_deps(igz, :, :) = -delta(:, :)*Vg(igz)
      END DO ! igz
!# 2555 "esm_stres_mod.f90"
      ! calculate stress tensor
      DO igz = 1, dfftp%nr3
        rg3 = rhog3(igz, imill_2d(0, 0))
        sigmaloclong(1:2, 1:2) = sigmaloclong(1:2, 1:2) &
                                 - REAL(CONJG(rg3)*dVg_deps(igz, 1:2, 1:2))
      END DO ! igz
    ENDIF ! imill_2d(0,0) > 0
!# 2563 "esm_stres_mod.f90"
    ! e2 means hartree -> Ry.
    sigmaloclong(:, :) = sigmaloclong(:, :)*(e2)
!# 2566 "esm_stres_mod.f90"
    CALL mp_sum(sigmaloclong, intra_bgrp_comm)
!# 2568 "esm_stres_mod.f90"
    DEALLOCATE (rhog3)
    DEALLOCATE (dVr_deps)
    DEALLOCATE (dVg_deps)
    DEALLOCATE (Vr)
    DEALLOCATE (Vg)
!# 2574 "esm_stres_mod.f90"
    RETURN
  END SUBROUTINE esm_stres_loclong_bc3
!# 2577 "esm_stres_mod.f90"
  COMPLEX(DP) FUNCTION qe_exp(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s
!# 2581 "esm_stres_mod.f90"
    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)
!# 2586 "esm_stres_mod.f90"
    qe_exp = EXP(r)*cmplx(c, s, kind=DP)
!# 2588 "esm_stres_mod.f90"
  END FUNCTION qe_exp
!# 2590 "esm_stres_mod.f90"
  COMPLEX(DP) FUNCTION qe_sinh(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s
!# 2594 "esm_stres_mod.f90"
    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)
!# 2599 "esm_stres_mod.f90"
    qe_sinh = 0.5d0*(EXP(r)*cmplx(c, s, kind=DP) - EXP(-r)*cmplx(c, -s, kind=DP))
!# 2601 "esm_stres_mod.f90"
  END FUNCTION qe_sinh
!# 2603 "esm_stres_mod.f90"
  COMPLEX(DP) FUNCTION qe_cosh(x)
    COMPLEX(DP), INTENT(in) :: x
    REAL(DP) :: r, i, c, s
!# 2607 "esm_stres_mod.f90"
    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)
!# 2612 "esm_stres_mod.f90"
    qe_cosh = 0.5d0*(EXP(r)*cmplx(c, s, kind=DP) + EXP(-r)*cmplx(c, -s, kind=DP))
!# 2614 "esm_stres_mod.f90"
  END FUNCTION qe_cosh
!# 2616 "esm_stres_mod.f90"
  FUNCTION qe_gauss(x) result(gauss)
    USE kinds,     ONLY : DP
    USE constants, ONLY : sqrtpm1  ! 1/sqrt(pi)
    IMPLICIT NONE
!# 2621 "esm_stres_mod.f90"
    REAL(DP), INTENT(in) :: x
    REAL(DP) :: gauss
!# 2624 "esm_stres_mod.f90"
    gauss = 2.0d0*sqrtpm1*EXP(-x*x)
!# 2626 "esm_stres_mod.f90"
  END FUNCTION qe_gauss
!# 2628 "esm_stres_mod.f90"
  FUNCTION exp_gauss( x, y )
    USE kinds,     ONLY : DP
    USE constants, ONLY : sqrtpm1 !1/sqrt(pi)
!# 2632 "esm_stres_mod.f90"
    REAL(DP), INTENT(IN) :: x, y
    REAL(DP) :: exp_gauss
!# 2635 "esm_stres_mod.f90"
    exp_gauss = 2._DP*sqrtpm1*EXP( x - y*y )
!# 2637 "esm_stres_mod.f90"
  END FUNCTION exp_gauss
END MODULE esm_stres_mod
