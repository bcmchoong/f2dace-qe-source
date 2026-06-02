!# 1 "paro_gamma_new.f90"
!
! Copyright (C) 2015-2016 Aihui Zhou's group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
!
! We propose some parallel orbital updating based plane wave basis methods
! for electronic structure calculations, which aims to the solution of the corresponding eigenvalue
! problems. Compared to the traditional plane wave methods, our methods have the feature of two level
! parallelization, which make them have great advantage in large-scale parallelization.
!
! The approach following Algorithm is the parallel orbital updating algorithm:
! 1. Choose initial $E_{\mathrm{cut}}^{(0)}$ and then obtain $V_{N_G^{0}}$, use the SCF method to solve
!    the Kohn-Sham equation in $V_{G_0}$ and get the initial $(\lambda_i^{0},u_i^{0}), i=1, \cdots, N$ 
!    and let $n=0$.
! 2. For $i=1,2,\ldots,N$, find $e_i^{n+1/2}\in V_{G_n}$ satisfying
!    $$a(\rho_{in}^{n}; e_i^{n+1/2}, v) = -[(a(\rho_{in}^{n}; u_i^{n}, v) - \lambda_i^{n} (u_i^{n}, v))]  $$
!    in parallel , where $\rho_{in}^{n}$ is the input charge density obtained by the orbits obtained in the 
!    $n$-th iteration or the former iterations.
! 3. Find $\{\lambda_i^{n+1},u_i^{n+1}\} \in \mathbf{R}\times \tilde{V}_N$   satisfying
!      $$a(\tilde{\rho}; u_i^{n+1}, v) = ( \lambda_i^{n+1}u_i^{n+1}, v) \quad  \forall v \in \tilde{V}_N$$
!      where $\tilde{V}_N = \mathrm{span}\{e_1^{n+1/2},\ldots,e_N^{n+1/2},u_1^{n},\ldots,u_N^{n}\}$, 
!      $\tilde{\rho}(x)$ is the input charge density obtained from the previous orbits.
! 4. Convergence check: if not converged, set $n=n+1$, go to step 2; else,  stop.
!
! You can see the detailed information through
!  X. Dai, X. Gong, A. Zhou, J. Zhu,
!   A parallel orbital-updating approach for electronic structure calculations, arXiv:1405.0260 (2014).
! X. Dai, Z. Liu, X. Zhang, A. Zhou,
!  A Parallel Orbital-updating Based Optimization Method for Electronic Structure Calculations, 
!   arXiv:1510.07230 (2015).
! Yan Pan, Xiaoying Dai, Xingao Gong, Stefano de Gironcoli, Gian-Marco Rignanese, and Aihui Zhou,
!  A Parallel Orbital-updating Based Plane Wave Basis Method. J. Comp. Phys. 348, 482-492 (2017).
!
! The file is written mainly by Stefano de Gironcoli and Yan Pan.
! GPU porting by Ivan Carnimeo
!
!NOTE (Ivan Carnimeo, May, 30th, 2022): 
!   paro_k_new and paro_gamma_new have been ported to GPU with OpenACC, 
!   the previous CUF versions (paro_k_new_gpu and paro_gamma_new_gpu) have been removed, 
!   and now paro_k_new and paro_gamma_new are used for both CPU and GPU execution.
!   If you want to see the previous code checkout to commit: 55c4e48ba650745f74bad43175f65f5449fd1273 (on Fri May 13 10:57:23 2022 +0000)
!
!-------------------------------------------------------------------------------
SUBROUTINE paro_gamma_new( h_psi_ptr, s_psi_ptr, hs_psi_ptr, g_psi_ptr, overlap, &
                 npwx, npw, nbnd, evc, eig, btype, ethr, notconv, nhpsi )
  !-------------------------------------------------------------------------------
  !paro_flag = 1: modified parallel orbital-updating method
!# 54 "paro_gamma_new.f90"
  ! global variables
  USE util_param,          ONLY : DP, stdout
  USE mp_bands_util,       ONLY : inter_bgrp_comm, nbgrp, my_bgrp_id
  USE mp,                  ONLY : mp_sum, mp_allgather, mp_barrier, &
                                  mp_type_create_column_section, mp_type_free
!# 60 "paro_gamma_new.f90"
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
!# 64 "paro_gamma_new.f90"
  ! I/O variables
  LOGICAL, INTENT(IN)        :: overlap
  INTEGER, INTENT(IN)        :: npw, npwx, nbnd
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nbnd)
  REAL(DP), INTENT(IN)       :: ethr
  REAL(DP), INTENT(INOUT)    :: eig(nbnd)
  INTEGER, INTENT(IN)        :: btype(nbnd)
  INTEGER, INTENT(OUT)       :: notconv, nhpsi
!  INTEGER, INTENT(IN)        :: paro_flag
  
  ! local variables (used in the call to cegterg )
  !------------------------------------------------------------------------
  EXTERNAL h_psi, s_psi, hs_psi, g_psi
  ! subroutine h_psi_ptr (npwx,npw,nvec,evc,hpsi)  computes H*evc  using band parallelization
  ! subroutine s_psi_ptr (npwx,npw,nvec,evc,spsi)  computes S*evc  using band parallelization
  ! subroutine hs_psi_ptr(npwx,npw,evc,hpsi,spsi)  computes H*evc and S*evc for a single band
  ! subroutine g_psi_ptr (npwx,npw,npol,m,psi,eig) computes g*psi  for m bands
  !
  ! ... local variables
  !
  INTEGER :: itry, paro_ntr, nconv, nextra, nactive, nbase, ntrust, ndiag, nvecx, nproc_ortho
  REAL(DP), ALLOCATABLE    :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
  LOGICAL, ALLOCATABLE     :: conv(:)
!# 89 "paro_gamma_new.f90"
  REAL(DP), PARAMETER      :: extra_factor = 0.5 ! workspace is at most this factor larger than nbnd
  INTEGER, PARAMETER       :: min_extra = 4      ! but at least this lager
!# 92 "paro_gamma_new.f90"
  INTEGER :: ibnd, ibnd_start, ibnd_end, how_many, lbnd, kbnd, last_unconverged, &
             recv_counts(nbgrp), displs(nbgrp), column_type
  !
  !$acc data deviceptr(eig)
  !
  ! ... init local variables
  !
  CALL laxlib_getval( nproc_ortho = nproc_ortho )
  paro_ntr = 20
!# 102 "paro_gamma_new.f90"
  nvecx = nbnd + max ( nint ( extra_factor * nbnd ), min_extra )
  !
  CALL start_clock( 'paro_gamma' ); !write (6,*) ' enter paro diag'
!# 106 "paro_gamma_new.f90"
  !$acc host_data use_device(evc)
  CALL mp_type_create_column_section(evc(1,1), 0, npwx, npwx, column_type)
  !$acc end host_data
!# 110 "paro_gamma_new.f90"
  ALLOCATE ( psi(npwx,nvecx), hpsi(npwx,nvecx), spsi(npwx,nvecx), ew(nvecx) )
  ALLOCATE ( conv(nbnd) )
  !$acc enter data create(psi, hpsi, spsi, ew)
!# 114 "paro_gamma_new.f90"
  CALL start_clock( 'paro:init' ); 
  conv(:) =  .FALSE. ; nconv = COUNT ( conv(:) )
  !$acc kernels
  psi(:,1:nbnd) = evc(:,1:nbnd) ! copy input evc into work vector
  !$acc end kernels
  !$acc data present(psi, spsi, hpsi)
  call h_psi (npwx,npw,nbnd,psi,hpsi) ! computes H*psi
  call s_psi (npwx,npw,nbnd,psi,spsi) ! computes S*psi
  !$acc end data
!# 124 "paro_gamma_new.f90"
  nhpsi = 0 ; IF (my_bgrp_id==0) nhpsi = nbnd
  CALL stop_clock( 'paro:init' ); 
!# 130 "paro_gamma_new.f90"
     CALL rotate_HSpsi_gamma (  npwx, npw, nbnd, nbnd, psi, hpsi, overlap, spsi, eig )
!# 140 "paro_gamma_new.f90"
  !write (6,'(10f10.4)') psi(1:5,1:3)
!# 142 "paro_gamma_new.f90"
  !write (6,*) eig(1:nbnd)
!# 144 "paro_gamma_new.f90"
  ParO_loop : &
  DO itry = 1,paro_ntr
!# 147 "paro_gamma_new.f90"
     !write (6,*) ' paro_itry =', itry, ethr
!# 149 "paro_gamma_new.f90"
     !----------------------------
!# 151 "paro_gamma_new.f90"
     nactive = nbnd - (nconv+1)/2 ! number of correction vectors to be computed (<nbnd)
     notconv = nbnd - nconv       ! number of needed roots
     nextra  = nactive - notconv  ! number of extra vectors
     nbase   = nconv + nactive    ! number of orbitals the correction should be orthogonal to (<2*nbnd)
     ndiag   = nbase + nactive    ! dimension of the matrix to be diagonalized at this iteration (<2*nbnd)
!# 157 "paro_gamma_new.f90"
     !----------------------------
!# 159 "paro_gamma_new.f90"
     nactive = min ( (nvecx-nconv)/2, nvecx-nbnd) ! number of corrections there is space for
     notconv = nbnd - nconv                       ! number of needed roots
     nextra  = max ( nactive - notconv, 0 )       ! number of extra vectors, if any
     nbase   = max ( nconv + nactive , nbnd )     ! number of orbitals to be orthogonal to  (<nvecx)
     ntrust  = min ( nconv + nactive , nbnd )     ! number of orbitals that will be actually corrected
     ndiag   = nbase + nactive       ! dimension of the matrix to be diagonalized at this iteration (<nvecx)
!# 166 "paro_gamma_new.f90"
     !write (6,*) itry, notconv, conv
     !write (6,*) ' nvecx, nbnd, nconv, notconv, nextra, nactive, nbase, ntrust, ndiag  =', nvecx, nbnd, nconv, notconv, nextra, nactive, nbase, ntrust, ndiag
     
     CALL divide_all(inter_bgrp_comm,nactive,ibnd_start,ibnd_end,recv_counts,displs)
     how_many = ibnd_end - ibnd_start + 1
     !write (6,*) nactive, ibnd_start, ibnd_end, recv_counts, displs
!# 173 "paro_gamma_new.f90"
     CALL start_clock( 'paro:pack' ); 
     lbnd = 1; kbnd = 1
     DO ibnd = 1, ntrust ! pack unconverged roots in the available space
        IF (.NOT.conv(ibnd) ) THEN
           !$acc kernels 
           psi (:,nbase+kbnd)  = psi(:,ibnd)
           hpsi(:,nbase+kbnd) = hpsi(:,ibnd)
           spsi(:,nbase+kbnd) = spsi(:,ibnd)
           ew(kbnd) = eig(ibnd) 
           !$acc end kernels
           last_unconverged = ibnd
           lbnd=lbnd+1 ; kbnd=kbnd+recv_counts(mod(lbnd-2,nbgrp)+1); if (kbnd>nactive) kbnd=kbnd+1-nactive
        END IF
     END DO
     DO ibnd = nbnd+1, nbase   ! add extra vectors if it is the case
        !$acc kernels 
        psi (:,nbase+kbnd)  = psi(:,ibnd)
        hpsi(:,nbase+kbnd) = hpsi(:,ibnd)
        spsi(:,nbase+kbnd) = spsi(:,ibnd)
        ew(kbnd) = eig(last_unconverged)
        !$acc end kernels
        lbnd=lbnd+1 ; kbnd=kbnd+recv_counts(mod(lbnd-2,nbgrp)+1); if (kbnd>nactive) kbnd=kbnd+1-nactive
     END DO
     CALL stop_clock( 'paro:pack' ); 
   
!     write (6,*) ' check nactive = ', lbnd, nactive
     if (lbnd .ne. nactive+1 ) stop ' nactive check FAILED '
!# 201 "paro_gamma_new.f90"
     CALL bpcg_gamma(hs_psi, g_psi, psi, spsi, npw, npwx, nbnd, how_many, &
                psi(:,nbase+1), hpsi(:,nbase+1), spsi(:,nbase+1), ethr, ew(1), nhpsi)
!# 204 "paro_gamma_new.f90"
     CALL start_clock( 'paro:mp_bar' ); 
     CALL mp_barrier(inter_bgrp_comm)
     CALL stop_clock( 'paro:mp_bar' ); 
     CALL start_clock( 'paro:mp_sum' ); 
!# 209 "paro_gamma_new.f90"
     !$acc host_data use_device(psi, hpsi, spsi)
     CALL mp_allgather(psi (:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     CALL mp_allgather(hpsi(:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     CALL mp_allgather(spsi(:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     !$acc end host_data
     CALL stop_clock( 'paro:mp_sum' ); 
!# 219 "paro_gamma_new.f90"
        !$acc host_data use_device(ew)
        CALL rotate_HSpsi_gamma (  npwx, npw, ndiag, ndiag, psi, hpsi, overlap, spsi, ew )
        !$acc end host_data
!# 232 "paro_gamma_new.f90"
     !write (6,*) ' ew : ', ew(1:nbnd)
     ! only the first nbnd eigenvalues are relevant for convergence
     ! but only those that have actually been corrected should be trusted
     conv(1:nbnd) = .FALSE.
!# 237 "paro_gamma_new.f90"
     !$acc kernels copy(conv) 
     conv(1:ntrust) = ABS(ew(1:ntrust)-eig(1:ntrust)).LT.ethr 
     !$acc end kernels
!# 241 "paro_gamma_new.f90"
     nconv = COUNT(conv(1:ntrust)) ; notconv = nbnd - nconv
     !$acc kernels 
     eig(1:nbnd)  = ew(1:nbnd)
     !$acc end kernels
     IF ( nconv == nbnd ) EXIT ParO_loop
!# 247 "paro_gamma_new.f90"
  END DO ParO_loop
!# 249 "paro_gamma_new.f90"
  !$acc kernels
  evc(:,1:nbnd) = psi(:,1:nbnd)
  !$acc end kernels
  !
  !$acc end data
  !
  CALL mp_sum(nhpsi,inter_bgrp_comm)
!# 257 "paro_gamma_new.f90"
  !$acc exit data delete(psi, hpsi, spsi, ew)
  DEALLOCATE ( ew, conv, psi, hpsi, spsi )
  CALL mp_type_free( column_type )
!# 261 "paro_gamma_new.f90"
  CALL stop_clock( 'paro_gamma' ); !write (6,*) ' exit paro diag'
!# 263 "paro_gamma_new.f90"
END SUBROUTINE paro_gamma_new
