!# 1 "casino_pp.f90"
!# 2 "casino_pp.f90"
MODULE casino_pp
!# 4 "casino_pp.f90"
  !
  ! All variables read from CASINO file format
  !
  ! trailing underscore means that a variable with the same name
  ! is used in module 'upf' containing variables to be written
  !
  USE upf_kinds, ONLY : dp
!# 12 "casino_pp.f90"
  CHARACTER(len=20) :: dft_
  CHARACTER(len=2)  :: psd_
  REAL(dp) :: zp_
  INTEGER nlc, nnl, lmax_, lloc, nchi, rel_
  LOGICAL :: numeric, bhstype, nlcc_
  CHARACTER(len=2), ALLOCATABLE :: els_(:)
  REAL(dp) :: zmesh
  REAL(dp) :: xmin      = -7.0_dp
  REAL(dp) :: dx        = 20.0_dp/1500.0_dp
  REAL(dp) :: tn_prefac = 0.75E-6_dp
  LOGICAL  :: tn_grid   = .true.
!# 25 "casino_pp.f90"
  REAL(dp), ALLOCATABLE::  r_(:)
  INTEGER :: mesh_
!# 28 "casino_pp.f90"
  REAL(dp), ALLOCATABLE::  vnl(:,:)
  INTEGER, ALLOCATABLE:: lchi_(:), nns_(:)
  REAL(dp), ALLOCATABLE:: chi_(:,:),  oc_(:)
!# 32 "casino_pp.f90"
CONTAINS
  !
  !     ----------------------------------------------------------
  SUBROUTINE read_casino(iunps,nofiles, waveunit)
    !     ----------------------------------------------------------
    !
    !     Reads in a CASINO tabulated pp file and it's associated
    !     awfn files. Some basic processing such as removing the
    !     r factors from the potentials is also performed.
!# 43 "casino_pp.f90"
    USE upf_kinds, ONLY : dp
    IMPLICIT NONE
    TYPE :: wavfun_list
       INTEGER :: occ,eup,edwn, nquant, lquant
       CHARACTER(len=2) :: label
       REAL(dp), ALLOCATABLE :: wavefunc(:)
       TYPE (wavfun_list), POINTER :: p
!# 51 "casino_pp.f90"
    END TYPE wavfun_list
!# 53 "casino_pp.f90"
    TYPE :: channel_list
       INTEGER :: lquant
       REAL(dp), ALLOCATABLE :: channel(:)
       TYPE (channel_list), POINTER :: p
!# 58 "casino_pp.f90"
    END TYPE channel_list
!# 61 "casino_pp.f90"
    TYPE (channel_list), POINTER :: phead
    TYPE (channel_list), POINTER :: pptr
    TYPE (channel_list), POINTER :: ptail
!# 65 "casino_pp.f90"
    TYPE (wavfun_list), POINTER :: mhead
    TYPE (wavfun_list), POINTER :: mptr
    TYPE (wavfun_list), POINTER :: mtail
!# 69 "casino_pp.f90"
    INTEGER :: iunps, nofiles, ios
    !
    LOGICAL :: groundstate, found
    CHARACTER(len=2) :: label, rellab
!# 74 "casino_pp.f90"
    INTEGER :: l, i, ir, nb, gsorbs, j,k,m,tmp, lquant, orbs, nquant
    INTEGER, ALLOCATABLE :: gs(:,:)
    INTEGER, INTENT(in) :: waveunit(nofiles)
!# 78 "casino_pp.f90"
    NULLIFY (  mhead, mptr, mtail )
    dft_ = 'HF'   !Hardcoded at the moment should eventually be HF anyway
!# 81 "casino_pp.f90"
    nlc = 0              !These two values are always 0 for numeric pps
    nnl = 0              !so lets just hard code them
!# 84 "casino_pp.f90"
    nlcc_ = .false.       !Again these two are alwas false for CASINO pps
    bhstype = .false.
!# 89 "casino_pp.f90"
    READ(iunps,'(a2,35x,a2)') rellab, psd_
    READ(iunps,*)
    IF ( rellab == 'DF' ) THEN
       rel_=1
    ELSE
       rel_=0
    ENDIF
!# 97 "casino_pp.f90"
    READ(iunps,*) zmesh,zp_  !Here we are reading zmesh (atomic #) and
    DO i=1,3                 !zp_ (pseudo charge)
       READ(iunps,*)
    ENDDO
    READ(iunps,*) lloc               !reading in lloc
    IF ( zp_<=0d0 ) &
         CALL upf_error( 'read_casino','Wrong zp ',1 )
    IF ( lloc>3.or.lloc<0 ) &
         CALL upf_error( 'read_casino','Wrong lloc ',1 )
!# 108 "casino_pp.f90"
    !
    !    compute the radial mesh
    !
!# 112 "casino_pp.f90"
    DO i=1,3
       READ(iunps,*)
    ENDDO
    READ(iunps,*) mesh_   !Reading in total no. of mesh points
!# 118 "casino_pp.f90"
    ALLOCATE(  r_(mesh_))
!# 120 "casino_pp.f90"
    READ(iunps,*)
    DO i=1,mesh_
       READ(iunps,*) r_(i)
    ENDDO
!# 126 "casino_pp.f90"
    ! Read in the different channels of V_nl
    ALLOCATE(phead)
    ptail => phead
    pptr  => phead
!# 131 "casino_pp.f90"
    ALLOCATE( pptr%channel(mesh_) )
    READ(iunps, '(15x,I1,7x)') l
    pptr%lquant=l
    READ(iunps, *)  (pptr%channel(ir),ir=1,mesh_)
!# 137 "casino_pp.f90"
    DO
       READ(iunps, '(15x,I1,7x)', IOSTAT=ios) l
!# 140 "casino_pp.f90"
       IF (ios /= 0 ) THEN
          exit
       ENDIF
!# 144 "casino_pp.f90"
       ALLOCATE(pptr%p)
       pptr=> pptr%p
       ptail=> pptr
       ALLOCATE( pptr%channel(mesh_) )
       pptr%lquant=l
       READ(iunps, *)  (pptr%channel(ir),ir=1,mesh_)
!# 151 "casino_pp.f90"
    ENDDO
!# 153 "casino_pp.f90"
    !Compute the number of channels read in.
    lmax_ =-1
    pptr => phead
    DO
       IF ( .not. associated(pptr) )exit
       lmax_=lmax_+1
!# 160 "casino_pp.f90"
       pptr =>pptr%p
    ENDDO
!# 163 "casino_pp.f90"
    ALLOCATE(vnl(mesh_,0:lmax_))
    i=0
    pptr => phead
    DO
       IF ( .not. associated(pptr) )exit
       !         lchi_(i) = pptr%lquant
!# 170 "casino_pp.f90"
       DO ir=1,mesh_
          vnl(ir,i) = pptr%channel(ir)
       ENDDO
       DEALLOCATE( pptr%channel )
       pptr =>pptr%p
       i=i+1
    ENDDO
!# 178 "casino_pp.f90"
    !Clean up the linked list (deallocate it)
    DO
       IF ( .not. associated(phead) )exit
       pptr => phead
       phead => phead%p
       DEALLOCATE( pptr )
    ENDDO
!# 186 "casino_pp.f90"
    DO l = 0, lmax_
       DO ir = 1, mesh_
          vnl(ir,l) = vnl(ir,l)/r_(ir) !Removing the factor of r CASINO has
       ENDDO
       ! correcting for possible divide by zero
       IF ( r_(1) == 0 ) THEN
          vnl(1,l) = 0
       ENDIF
    ENDDO
!# 196 "casino_pp.f90"
    ALLOCATE(mhead)
    mtail => mhead
!# 199 "casino_pp.f90"
    mptr => mhead
!# 201 "casino_pp.f90"
    NULLIFY(mtail%p)
    groundstate=.true.
    DO j=1,nofiles
!# 205 "casino_pp.f90"
       DO i=1,4
!# 207 "casino_pp.f90"
          READ(waveunit(j),*)
       ENDDO
!# 210 "casino_pp.f90"
       READ(waveunit(j),*) orbs
!# 212 "casino_pp.f90"
       IF ( groundstate ) THEN
!# 214 "casino_pp.f90"
          ALLOCATE( gs(orbs,3) )
!# 216 "casino_pp.f90"
          gs = 0
          gsorbs = orbs
       ENDIF
!# 220 "casino_pp.f90"
       DO i=1,2
          READ(waveunit(j),*)
       ENDDO
!# 224 "casino_pp.f90"
       READ(waveunit(j),*) mtail%eup, mtail%edwn
       READ(waveunit(j),*)
!# 227 "casino_pp.f90"
       DO i=1,mtail%eup+mtail%edwn
          READ(waveunit(j),*) tmp, nquant, lquant
!# 230 "casino_pp.f90"
          IF ( groundstate ) THEN
             found = .true.
!# 233 "casino_pp.f90"
             DO m=1,orbs
!# 235 "casino_pp.f90"
                IF ( (nquant==gs(m,1) .and. lquant==gs(m,2)) ) THEN
                   gs(m,3) = gs(m,3) + 1
                   exit
                ENDIF
!# 240 "casino_pp.f90"
                found = .false.
!# 242 "casino_pp.f90"
             ENDDO
!# 244 "casino_pp.f90"
             IF (.not. found ) THEN
!# 246 "casino_pp.f90"
                DO m=1,orbs
!# 248 "casino_pp.f90"
                   IF ( gs(m,1) == 0 ) THEN
                      gs(m,1) = nquant
                      gs(m,2) = lquant
                      gs(m,3) = 1
!# 253 "casino_pp.f90"
                      exit
                   ENDIF
!# 256 "casino_pp.f90"
                ENDDO
!# 258 "casino_pp.f90"
             ENDIF
!# 260 "casino_pp.f90"
          ENDIF
!# 262 "casino_pp.f90"
       ENDDO
!# 264 "casino_pp.f90"
       READ(waveunit(j),*)
       READ(waveunit(j),*)
!# 267 "casino_pp.f90"
       DO i=1,mesh_
          READ(waveunit(j),*)
       ENDDO
!# 271 "casino_pp.f90"
       DO k=1,orbs
          READ(waveunit(j),'(13x,a2)', err=300) label
          READ(waveunit(j),*) tmp, nquant, lquant
!# 275 "casino_pp.f90"
          IF ( .not. groundstate ) THEN
             found = .false.
!# 278 "casino_pp.f90"
             DO m = 1,gsorbs
!# 280 "casino_pp.f90"
                IF ( nquant == gs(m,1) .and. lquant == gs(m,2) ) THEN
                   found = .true.
                   exit
                ENDIF
             ENDDO
             mptr => mhead
             DO
                IF ( .not. associated(mptr) )exit
                IF ( nquant == mptr%nquant .and. lquant == mptr%lquant ) found = .true.
                mptr =>mptr%p
             ENDDO
             IF ( found ) THEN
                DO i=1,mesh_
                   READ(waveunit(j),*)
                ENDDO
!# 296 "casino_pp.f90"
                CYCLE
             ENDIF
          ENDIF
            IF ( allocated(mtail%wavefunc) ) THEN
                ALLOCATE(mtail%p)
                mtail=>mtail%p
                NULLIFY(mtail%p)
                ALLOCATE( mtail%wavefunc(mesh_) )
             ELSE
                ALLOCATE( mtail%wavefunc(mesh_) )
             ENDIF
             mtail%label = label
             mtail%nquant = nquant
             mtail%lquant = lquant
!# 312 "casino_pp.f90"
             READ(waveunit(j), *, err=300) (mtail%wavefunc(ir),ir=1,mesh_)
          ENDDO
          groundstate = .false.
       ENDDO
!# 317 "casino_pp.f90"
       nchi =0
       mptr => mhead
       DO
          IF ( .not. associated(mptr) )exit
          nchi=nchi+1
!# 323 "casino_pp.f90"
          mptr =>mptr%p
       ENDDO
!# 326 "casino_pp.f90"
       ALLOCATE(lchi_(nchi), els_(nchi), nns_(nchi))
       ALLOCATE(oc_(nchi))
       ALLOCATE(chi_(mesh_,nchi))
       oc_ = 0
!# 331 "casino_pp.f90"
       !Sort out the occupation numbers
       DO i=1,gsorbs
          oc_(i)=gs(i,3)
       ENDDO
       DEALLOCATE( gs )
!# 337 "casino_pp.f90"
       i=1
       mptr => mhead
       DO
          IF ( .not. associated(mptr) )exit
          nns_(i) = mptr%nquant
          lchi_(i) = mptr%lquant
          els_(i) = mptr%label
!# 345 "casino_pp.f90"
          DO ir=1,mesh_
!# 347 "casino_pp.f90"
             chi_(ir:,i) = mptr%wavefunc(ir)
          ENDDO
          DEALLOCATE( mptr%wavefunc )
          mptr =>mptr%p
          i=i+1
       ENDDO
!# 354 "casino_pp.f90"
       !Clean up the linked list (deallocate it)
       DO
          IF ( .not. associated(mhead) )exit
          mptr => mhead
          mhead => mhead%p
          DEALLOCATE( mptr )
       ENDDO
!# 363 "casino_pp.f90"
       !     ----------------------------------------------------------
       WRITE (0,'(a)') 'Pseudopotential successfully read'
       !     ----------------------------------------------------------
       RETURN
!# 368 "casino_pp.f90"
300    CALL upf_error('read_casino','pseudo file is empty or wrong',1)
!# 370 "casino_pp.f90"
     END SUBROUTINE read_casino
!# 372 "casino_pp.f90"
     !     ----------------------------------------------------------
     SUBROUTINE convert_casino(upf_out)
       !     ----------------------------------------------------------
       USE upf_kinds,    ONLY : dp
       USE pseudo_types, ONLY : pseudo_upf
!# 378 "casino_pp.f90"
       IMPLICIT NONE
!# 380 "casino_pp.f90"
       TYPE(pseudo_upf), INTENT(inout)       :: upf_out
!# 382 "casino_pp.f90"
       REAL(dp), ALLOCATABLE :: aux(:)
       REAL(dp) :: vll
       INTEGER :: kkbeta, l, iv, ir, i, nb
!# 386 "casino_pp.f90"
       !
       upf_out%nv       = "2.0.1"
       upf_out%tvanp    = .false.
       upf_out%tpawp    = .false.
       upf_out%tcoulombp= .false.
       upf_out%has_so   = .false.
       upf_out%has_wfc  = .false.
       upf_out%has_gipaw= .false.
       upf_out%paw_as_gipaw = .false.
       upf_out%with_metagga_info = .false.
       !
       WRITE(upf_out%generated, '("From a Trail & Needs tabulated &
            &PP for CASINO")')
       WRITE(upf_out%author,'("unknown")')
       WRITE(upf_out%date,'("unknown")')
       upf_out%comment = 'Info: automatically converted from CASINO &
            &Tabulated format'
!# 404 "casino_pp.f90"
       IF (rel_== 0) THEN
          upf_out%rel = 'no'
       ELSEIF (rel_==1 ) THEN
          upf_out%rel = 'scalar'
       ELSE
          upf_out%rel = 'full'
       ENDIF
!# 412 "casino_pp.f90"
       IF (xmin == 0 ) THEN
          xmin= log(zmesh * r_(2) )
       ENDIF
!# 416 "casino_pp.f90"
       ! Allocate and assign the radial grid
!# 418 "casino_pp.f90"
       upf_out%mesh  = mesh_
       upf_out%zmesh = zmesh
       upf_out%dx    = dx
       upf_out%xmin  = xmin
!# 423 "casino_pp.f90"
       ALLOCATE(upf_out%rab(upf_out%mesh))
       ALLOCATE(  upf_out%r(upf_out%mesh))
!# 426 "casino_pp.f90"
       upf_out%r = r_
       DEALLOCATE( r_ )
!# 429 "casino_pp.f90"
       upf_out%rmax = maxval(upf_out%r)
!# 432 "casino_pp.f90"
       !
       ! subtract out the local part from the different
       ! potential channels
       !
!# 437 "casino_pp.f90"
       DO l = 0, lmax_
          IF ( l/=lloc ) vnl(:,l) = vnl(:,l) - vnl(:,lloc)
       ENDDO
!# 441 "casino_pp.f90"
       ALLOCATE (upf_out%vloc(upf_out%mesh))
       upf_out%vloc(:) = vnl(:,lloc)
!# 445 "casino_pp.f90"
       ! Compute the derivatives of the grid. The Trail and Needs
       ! grids use r(i) = (tn_prefac / zmesh)*( exp(i*dx) - 1 ) so
       ! must be treated differently to standard QE grids.
!# 449 "casino_pp.f90"
       IF ( tn_grid ) THEN
          DO ir = 1, upf_out%mesh
             upf_out%rab(ir) = dx * ( upf_out%r(ir) + tn_prefac / zmesh )
          ENDDO
       ELSE
          DO ir = 1, upf_out%mesh
             upf_out%rab(ir) = dx  * upf_out%r(ir)
          ENDDO
       ENDIF
!# 460 "casino_pp.f90"
       !
       !    compute the atomic charges
       !
       ALLOCATE (upf_out%rho_at(upf_out%mesh))
       upf_out%rho_at(:) = 0.d0
!# 466 "casino_pp.f90"
       DO nb = 1, nchi
          IF( oc_(nb)/=0.d0) THEN
             upf_out%rho_at(:) = upf_out%rho_at(:) +&
               &  oc_(nb)*chi_(:,nb)**2
          ENDIF
       ENDDO
!# 473 "casino_pp.f90"
       ! This section deals with the pseudo wavefunctions.
       ! These values are just given directly to the pseudo_upf structure
       upf_out%nwfc  = nchi
!# 477 "casino_pp.f90"
       ALLOCATE( upf_out%oc(upf_out%nwfc), upf_out%epseu(upf_out%nwfc) )
       ALLOCATE( upf_out%lchi(upf_out%nwfc), upf_out%nchi(upf_out%nwfc) )
       ALLOCATE( upf_out%els(upf_out%nwfc) )
       ALLOCATE( upf_out%rcut_chi(upf_out%nwfc) )
       ALLOCATE( upf_out%rcutus_chi (upf_out%nwfc) )
!# 483 "casino_pp.f90"
       DO i=1, upf_out%nwfc
          upf_out%nchi(i)  = nns_(i)
          upf_out%lchi(i)  = lchi_(i)
          upf_out%rcut_chi(i)  = 0.0d0
          upf_out%rcutus_chi(i)= 0.0d0
          upf_out%oc (i)   = oc_(i)
          upf_out%els(i) = els_(i)
          upf_out%epseu(i) = 0.0d0
       ENDDO
       DEALLOCATE (lchi_, oc_, nns_)
!# 494 "casino_pp.f90"
       upf_out%psd = psd_
       upf_out%typ = 'NC'
       upf_out%nlcc = nlcc_
       upf_out%zp = zp_
       upf_out%etotps = 0.0d0
       upf_out%ecutrho=0.0d0
       upf_out%ecutwfc=0.0d0
       upf_out%lloc=lloc
!# 503 "casino_pp.f90"
       IF ( lmax_ == lloc) THEN
          upf_out%lmax = lmax_-1
       ELSE
          upf_out%lmax = lmax_
       ENDIF
       upf_out%nbeta = lmax_
!# 510 "casino_pp.f90"
       ALLOCATE ( upf_out%els_beta(upf_out%nbeta) )
       ALLOCATE ( upf_out%rcut(upf_out%nbeta) )
       ALLOCATE ( upf_out%rcutus(upf_out%nbeta) )
!# 514 "casino_pp.f90"
       upf_out%rcut=0.0d0
       upf_out%rcutus=0.0d0
       upf_out%dft =dft_
!# 519 "casino_pp.f90"
       IF (upf_out%nbeta > 0) THEN
!# 521 "casino_pp.f90"
          ALLOCATE(upf_out%kbeta(upf_out%nbeta), upf_out%lll(upf_out%nbeta))
          upf_out%kkbeta=upf_out%mesh
          DO ir = 1,upf_out%mesh
             IF ( upf_out%r(ir) > upf_out%rmax ) THEN
                upf_out%kkbeta=ir
                exit
             ENDIF
          ENDDO
!# 530 "casino_pp.f90"
          ! make sure kkbeta is odd as required for simpson
          IF(mod(upf_out%kkbeta,2) == 0) upf_out%kkbeta=upf_out%kkbeta-1
          upf_out%kbeta(:) = upf_out%kkbeta
          ALLOCATE(aux(upf_out%kkbeta))
          ALLOCATE(upf_out%beta(upf_out%mesh,upf_out%nbeta))
          ALLOCATE(upf_out%dion(upf_out%nbeta,upf_out%nbeta))
!# 537 "casino_pp.f90"
          upf_out%dion(:,:) =0.d0
!# 539 "casino_pp.f90"
          iv=0
          DO i=1,upf_out%nwfc
             l=upf_out%lchi(i)
             IF (l/=upf_out%lloc) THEN
                iv=iv+1
                upf_out%els_beta(iv)=upf_out%els(i)
                upf_out%lll(iv)=l
                DO ir=1,upf_out%kkbeta
!# 548 "casino_pp.f90"
                   upf_out%beta(ir,iv)=chi_(ir,i)*vnl(ir,l)
                   aux(ir) = chi_(ir,i)**2*vnl(ir,l)
!# 551 "casino_pp.f90"
                ENDDO
                CALL simpson(upf_out%kkbeta,aux,upf_out%rab,vll)
                upf_out%dion(iv,iv) = 1.0d0/vll
             ENDIF
!# 556 "casino_pp.f90"
             IF(iv >= upf_out%nbeta) exit  ! skip additional pseudo wfns
          ENDDO
!# 560 "casino_pp.f90"
          DEALLOCATE (vnl, aux)
!# 562 "casino_pp.f90"
          !
          !   redetermine ikk2
          !
          DO iv=1,upf_out%nbeta
             upf_out%kbeta(iv)=upf_out%kkbeta
             DO ir = upf_out%kkbeta,1,-1
                IF ( abs(upf_out%beta(ir,iv)) > 1.d-12 ) THEN
                   upf_out%kbeta(iv)=ir
                   exit
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!# 576 "casino_pp.f90"
       ALLOCATE (upf_out%chi(upf_out%mesh,upf_out%nwfc))
       upf_out%chi = chi_
       DEALLOCATE (chi_)
!# 580 "casino_pp.f90"
       RETURN
     END SUBROUTINE convert_casino
!# 584 "casino_pp.f90"
     SUBROUTINE write_casino_tab(upf_in, fileout)
!# 586 "casino_pp.f90"
       USE pseudo_types, ONLY : pseudo_upf
!# 588 "casino_pp.f90"
       IMPLICIT NONE
!# 590 "casino_pp.f90"
       CHARACTER(LEN=*), INTENT(in)       :: fileout
       TYPE(pseudo_upf), INTENT(in)       :: upf_in
       INTEGER :: i, lp1, unout_
!# 594 "casino_pp.f90"
       INTEGER, EXTERNAL :: atomic_number
!# 596 "casino_pp.f90"
       OPEN ( NEWUNIT=unout_, FILE=TRIM(fileout), ACTION = 'WRITE') 
!# 598 "casino_pp.f90"
       WRITE(unout_,*) "Converted Pseudopotential in REAL space for ", upf_in%psd
       WRITE(unout_,*) "Atomic number and pseudo-charge"
       WRITE(unout_,"(I3,F8.2)") atomic_number( upf_in%psd ),upf_in%zp
       WRITE(unout_,*) "Energy units (rydberg/hartree/ev):"
       WRITE(unout_,*) "rydberg"
       WRITE(unout_,*) "Angular momentum of local component (0=s,1=p,2=d..)"
       WRITE(unout_,"(I2)") upf_in%lloc
       WRITE(unout_,*) "NLRULE override (1) VMC/DMC (2) config gen (0 ==> &
            &input/default VALUE)"
       WRITE(unout_,*) "0 0"
       WRITE(unout_,*) "Number of grid points"
       WRITE(unout_,*) upf_in%mesh
       WRITE(unout_,*) "R(i) in atomic units"
       WRITE(unout_, "(T4,E22.15)") upf_in%r(:)
!# 613 "casino_pp.f90"
       lp1 = size ( vnl, 2 )
       DO i=1,lp1
          WRITE(unout_, "(A,I1,A)") 'r*potential (L=',i-1,') in Ry'
          WRITE(unout_, "(T4,E22.15)") vnl(:,i)
       ENDDO
       CLOSE (unout_) 
       DEALLOCATE(vnl)
!# 621 "casino_pp.f90"
     END SUBROUTINE write_casino_tab
!# 623 "casino_pp.f90"
     SUBROUTINE conv_upf2casino(upf_in)
!# 625 "casino_pp.f90"
       USE pseudo_types, ONLY : pseudo_upf
!# 628 "casino_pp.f90"
       IMPLICIT NONE
!# 630 "casino_pp.f90"
       TYPE(pseudo_upf), INTENT(in)       :: upf_in
       INTEGER :: i, l, channels
!# 633 "casino_pp.f90"
       REAL(dp), PARAMETER :: offset=1E-20_dp
       !This is an offset added to the wavefunctions to
       !eliminate any divide by zeros that may be caused by
       !zeroed wavefunction terms.
!# 638 "casino_pp.f90"
       IF (upf_in%typ /= 'NC') THEN
          WRITE(0,*) ''
          WRITE(0,*) 'WRONG PSEUDOPOTENTIAL!'
          WRITE(0,*) 'Only norm-conserving pps can be used in CASINO!'
          STOP
       ENDIF
       
       WRITE(0,*) "Number of grid points: ", upf_in%mesh
       WRITE(0,*) "Number of KB projectors: ", upf_in%nbeta
       WRITE(0,*) "Channel(s) of KB projectors: ", upf_in%lll
       WRITE(0,*) "Number of channels to be re-constructed: ", upf_in%nbeta+1
!# 650 "casino_pp.f90"
       channels=upf_in%nbeta+1
       ALLOCATE ( vnl(upf_in%mesh,channels) )
!# 653 "casino_pp.f90"
       !Set up the local component of each channel
       DO i=1,channels
          vnl(:,i)=upf_in%r(:)*upf_in%vloc(:)
       ENDDO
!# 659 "casino_pp.f90"
       DO i=1,upf_in%nbeta
          l=upf_in%lll(i)+1
!# 662 "casino_pp.f90"
          !Check if any wfc components have been zeroed
          !and apply the offset IF they have
!# 665 "casino_pp.f90"
          IF ( minval(abs(upf_in%chi(:,l))) /= 0 ) THEN
             vnl(:,l)= (upf_in%beta(:,i)/(upf_in%chi(:,l)) &
                  *upf_in%r(:)) + vnl(:,l)
          ELSE
             WRITE(0,"(A,ES10.3,A)") 'Applying ',offset , ' offset to &
                  &wavefunction to avoid divide by zero'
             vnl(:,l)= (upf_in%beta(:,i)/(upf_in%chi(:,l)+offset) &
                  *upf_in%r(:)) + vnl(:,l)
          ENDIF
!# 675 "casino_pp.f90"
       ENDDO
!# 677 "casino_pp.f90"
     END SUBROUTINE conv_upf2casino
!# 679 "casino_pp.f90"
END MODULE casino_pp
