!# 1 "io_base.f90"
!
! Copyright (C) 2016-2017 Quantum ESPRESSO Foundation 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_base
  !----------------------------------------------------------------------------
  !
  !! Subroutines used to read and write binary data produced by QE.  
  !! Author: Paolo Giannozzi, based on previous work by Carlo Cavazzoni
  !
  USE kinds,     ONLY : dp
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: write_wfc, read_wfc, write_rhog, read_rhog
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, filename, root_in_group, intra_group_comm, &
         ik, xk, ispin, nspin, wfc, ngw, gamma_only, nbnd, igl, ngwl, &
         b1,b2,b3, mill_k, scalef )
      !------------------------------------------------------------------------
      !
      !! Collects wfc, distributed on "intra_group_comm", writes them
      !! together with related information to file "filename.*"
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Only processor "root_in_group" collects data and writes to file
      !!
      USE mp_wave,    ONLY : mergewf, mergekg
      USE mp,         ONLY : mp_size, mp_rank, mp_max
      !
!# 44 "io_base.f90"
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      CHARACTER(LEN=*),   INTENT(IN) :: filename
      INTEGER,            INTENT(IN) :: ik, ispin, nspin
      REAL(DP),           INTENT(IN) :: xk(:)
      COMPLEX(DP),        INTENT(IN) :: wfc(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      INTEGER,            INTENT(IN) :: mill_k(:,:)
      REAL(DP),           INTENT(IN) :: b1(3), b2(3), b3(3)    
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm
      !
      LOGICAL                  :: ionode_in_group
      INTEGER                  :: igwx, npwx, npol, j
      INTEGER                  :: me_in_group, nproc_in_group, my_group
      INTEGER, ALLOCATABLE     :: itmp(:,:)
      COMPLEX(DP), ALLOCATABLE, TARGET :: wtmp(:)
      COMPLEX(DP), POINTER             :: wtmp2(:)
      !
!# 74 "io_base.f90"
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      npol = 1
      IF ( nspin == 4 ) npol = 2
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN
!# 101 "io_base.f90"
         OPEN ( UNIT = iuni, FILE = TRIM(filename)//'.dat', &
              FORM='unformatted', STATUS = 'unknown' )
         WRITE(iuni) ik, xk, ispin, gamma_only, scalef
         WRITE(iuni) ngw, igwx, npol, nbnd
!# 106 "io_base.f90"
         !
      END IF
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( itmp( 3, MAX (igwx,1) ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      itmp (:,:) = 0
      CALL mergekg( mill_k, itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      IF ( ionode_in_group ) THEN
!# 132 "io_base.f90"
         WRITE(iuni) b1, b2, b3
         WRITE(iuni) itmp(1:3,1:igwx)
!# 135 "io_base.f90"
      END IF
      DEALLOCATE( itmp )
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( wtmp( MAX( npol*igwx, 1 ) ) )
         IF ( npol == 2 ) wtmp2 => wtmp( igwx+1:2*igwx )
      ELSE
         ALLOCATE( wtmp( 1 ) )
         IF ( npol == 2 ) wtmp2 => wtmp( 1:1 )
      ENDIF
      wtmp = 0.0_DP
      !
!# 157 "io_base.f90"
      DO j = 1, nbnd
         !
         IF ( npol == 2 ) THEN
            !
            ! Quick-and-dirty noncolinear case - mergewf should be modified
            ! Collect into wtmp(1:igwx) the first set of plane waves components
            !
            CALL mergewf( wfc(1:npwx,       j), wtmp , ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
            ! Collect into wtmp(igwx+1:2*igwx) the second set of plane waves
            ! components - pointer wtmp2 is used instead of wtmp(igwx+1:2*igwx)
            ! in order to avoid a bogus out-of-bound error
            !
            CALL mergewf( wfc(npwx+1:2*npwx,j), wtmp2, ngwl, igl,&
                 me_in_group, nproc_in_group, root_in_group, intra_group_comm )
            !
         ELSE
            !
            CALL mergewf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                 nproc_in_group, root_in_group, intra_group_comm )
            !
         END IF
         !
         IF ( ionode_in_group ) THEN
!# 186 "io_base.f90"
            WRITE(iuni) wtmp(1:npol*igwx)
!# 188 "io_base.f90"
         END IF
         !
      END DO
      IF ( ionode_in_group ) THEN
!# 196 "io_base.f90"
         CLOSE (UNIT = iuni, STATUS = 'keep' )
!# 198 "io_base.f90"
      END IF
      !
      IF ( npol == 2 ) NULLIFY ( wtmp2 )
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, filename, root_in_group, intra_group_comm,  &
         ik, xk, ispin, npol, wfc, ngw, gamma_only, nbnd, igl, ngwl, &
         b1, b2, b3, mill_k, scalef, ierr )
      !
      !! Processor "root_in_group" reads wfc and related information from file 
      !! "filename.*" (* = dat if fortran binary, * = hdf5 if HDF5),
      !! distributes wfc on "intra_group_comm"
      !! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp_wave,     ONLY : splitwf, splitkg
      USE mp,          ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
!# 225 "io_base.f90"
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER,            INTENT(IN)    :: ik
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, npol
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      REAL(DP),           INTENT(OUT)   :: xk(3)
      REAL(DP),           INTENT(OUT)   :: b1(3), b2(3), b3(3)
      INTEGER,            INTENT(OUT)   :: mill_k(:,:)
      LOGICAL,            INTENT(OUT)   :: gamma_only
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      INTEGER, ALLOCATABLE              :: itmp(:,:)
      COMPLEX(DP), ALLOCATABLE, TARGET  :: wtmp(:)
      COMPLEX(DP), POINTER              :: wtmp2(:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, ik_, nbnd_
      INTEGER                           :: me_in_group, nproc_in_group
      LOGICAL                           :: ionode_in_group
!# 255 "io_base.f90"
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      !
      IF ( ionode_in_group ) THEN
!# 265 "io_base.f90"
         OPEN ( UNIT = iuni, FILE=TRIM(filename)//'.dat', &
                FORM='unformatted', STATUS = 'old', IOSTAT = ierr_)
!# 270 "io_base.f90"
      END IF
      CALL mp_bcast( ierr_, root_in_group, intra_group_comm )
      IF ( PRESENT(ierr) ) THEN
         ierr = ierr_
         IF ( ierr /= 0 ) RETURN
      ELSE
         CALL errore( 'read_wfc ', &
              'cannot open restart file ' // TRIM(filename) //' for reading', ierr_ )
      END IF
      !
      IF ( ionode_in_group ) THEN
!# 297 "io_base.f90"
         READ (iuni) ik_, xk, ispin, gamma_only, scalef
         READ (iuni) ngw, igwx_, npol, nbnd_
!# 300 "io_base.f90"
      END IF
      !
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( xk,     root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only, root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( npol,   root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd_,   root_in_group, intra_group_comm )
      !
      npwx = SIZE( wfc, 1 ) / npol
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( itmp( 3,MAX( igwx_, igwx ) ) )
!# 324 "io_base.f90"
         READ (iuni) b1, b2, b3
         READ (iuni) itmp(1:3,1:igwx_)
!# 327 "io_base.f90"
         IF ( igwx > igwx_ ) itmp(1:3,igwx_+1:igwx) = 0
      ELSE
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      CALL splitkg( mill_k(:,:), itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      DEALLOCATE (itmp)
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( wtmp( npol*MAX( igwx_, igwx ) ) )
         IF ( npol == 2 ) wtmp2 => wtmp(igwx_+1:2*igwx_)
!# 342 "io_base.f90"
      ELSE
         ALLOCATE( wtmp(1) )
         IF ( npol == 2 ) wtmp2 => wtmp( 1:1 )
      ENDIF
      nbnd = nbnd_ 
      DO j = 1, nbnd_ 
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
!# 357 "io_base.f90"
               READ (iuni) wtmp(1:npol*igwx_) 
!# 359 "io_base.f90"
               IF ( igwx > igwx_ ) wtmp((npol*igwx_+1):npol*igwx) = 0.0_DP
               !
            END IF
            !
            IF ( npol == 2 ) THEN
               !
               ! Quick-and-dirty noncolinear case - mergewf should be modified
               ! Collect into wtmp(1:igwx_) first set of plane wave components
               !
               CALL splitwf( wfc(1:npwx,       j), wtmp ,   &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
               !
               ! Collect into wtmp(igwx_+1:2*igwx_) the second set of plane wave
               ! components - instead of wtmp(igwx_+1:2*igwx_), pointer wtmp2
               ! is used, in order to prevent a bogus out-of-bound error
               !
               CALL splitwf( wfc(npwx+1:2*npwx,j), wtmp2,  &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
            ELSE
               CALL splitwf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                    nproc_in_group, root_in_group, intra_group_comm )
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode_in_group ) THEN
!# 393 "io_base.f90"
         CLOSE ( UNIT = iuni, STATUS = 'keep' )
!# 395 "io_base.f90"
      END IF
      !
      IF ( npol == 2 ) NULLIFY ( wtmp2 )
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rhog ( filename, root_in_group, intra_group_comm, &
         b1, b2, b3, gamma_only, mill, ig_l2g, rho )
      !------------------------------------------------------------------------
      !! Collects rho(G), distributed on "intra_group_comm", writes it
      !! together with related information to file "filename".*
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Processor "root_in_group" collects data and writes to file
      !
      USE mp,                   ONLY : mp_sum, mp_bcast, mp_size, mp_rank
      USE mp_wave,              ONLY : mergewf, mergekg
!# 418 "io_base.f90"
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !! name of file written (to which a suffix is added)
      INTEGER,            INTENT(IN) :: root_in_group
      !! root processor that collects and writes
      INTEGER,            INTENT(IN) :: intra_group_comm
      !! rho(G) is distributed over this group of processors
      REAL(dp),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !!  b1, b2, b3 are the three primitive vectors in a.u.
      INTEGER,          INTENT(IN) :: mill(:,:)
      !! Miller indices for local G-vectors
      !! G = mill(1)*b1 + mill(2)*b2 + mill(3)*b3
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      LOGICAL,          INTENT(IN) :: gamma_only
      !! if true, only the upper half of G-vectors (z >=0) is present
      COMPLEX(dp),      INTENT(IN) :: rho(:,:)
      !! rho(G) on this processor
      !
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      !! Local rho(G), with LSDA workaround
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      !! Global rho(G) collected on root proc
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !! Global Miller indices collected on root proc
      INTEGER                  :: me_in_group, nproc_in_group
      LOGICAL                  :: ionode_in_group
      INTEGER                  :: ngm, nspin, ngm_g, igwx
      INTEGER                  :: iun, ns, ig, ierr
      !
!# 457 "io_base.f90"
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (mill, 2) .OR. ngm /= SIZE (ig_l2g, 1) ) &
         CALL errore('write_rhog', 'inconsistent input dimensions', 1)
      nspin= SIZE (rho, 2)
!# 474 "io_base.f90"
      iun  = 4
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_group_comm )
      !
      ierr = 0
!# 486 "io_base.f90"
      IF ( ionode_in_group ) OPEN ( UNIT = iun, FILE = TRIM(filename)//'.dat', &
                FORM = 'unformatted', STATUS = 'unknown', iostat = ierr )
!# 489 "io_base.f90"
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error opening file ' &
           & // TRIM( filename ), 1 )
      IF ( ionode_in_group ) THEN
!# 499 "io_base.f90"
          WRITE (iun, iostat=ierr) gamma_only, ngm_g, nspin
          WRITE (iun, iostat=ierr) b1, b2, b3
!# 502 "io_base.f90"
      END IF
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 1 )
      !
      ! ... collect all G-vectors across processors within the band group
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( mill_g( 3, ngm_g ) )
      ELSE
         ! not used: some compiler do not like passing unallocated arrays
         ALLOCATE( mill_g( 3, 1 ) )
      END IF
      !
      ! ... mergekg collects distributed array mill(1:3,ig) where ig is the
      ! ... local index, into array mill_g(1:3,ig_g), where ig_g=ig_l2g(ig)
      ! ... is the global index. mill_g is collected on root_bgrp only
      !
      CALL mergekg( mill, mill_g, ngm, ig_l2g, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      !
      ! ... write G-vectors
      !
      IF ( ionode_in_group ) THEN
!# 538 "io_base.f90"
         WRITE (iun, iostat=ierr) mill_g(1:3,1:ngm_g)
!# 540 "io_base.f90"
      END IF
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... deallocate to save memory
      !
      DEALLOCATE( mill_g )
      !
      ! ... now collect all G-vector components of the charge density
      ! ... (one spin at the time to save memory) using the same logic
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( rho_g( ngm_g ) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
      !
      DO ns = 1, nspin
         !
         DO ig = 1, ngm
               rhoaux(ig) = rho(ig,ns)
         END DO
         !
         rho_g = 0
         CALL mergewf( rhoaux, rho_g, ngm, ig_l2g, me_in_group, &
              nproc_in_group, root_in_group, intra_group_comm )
         !
         IF ( ionode_in_group ) THEN
!# 577 "io_base.f90"
            WRITE (iun, iostat=ierr) rho_g(1:ngm_g)
!# 579 "io_base.f90"
         END IF
         CALL mp_bcast( ierr, root_in_group, intra_group_comm )
         IF ( ierr > 0 ) CALL errore ( 'write_rhog','error writing file ' &
              & // TRIM( filename ), 2+ns )
         !
      END DO
      !
!# 589 "io_base.f90"
      IF (ionode_in_group) CLOSE (UNIT = iun, status ='keep' )
!# 591 "io_base.f90"
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      !
      RETURN
      !
    END SUBROUTINE write_rhog
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rhog ( filename, root_in_group, intra_group_comm, &
         ig_l2g, nspin, rho, gamma_only, ier_ )
      !------------------------------------------------------------------------
      !! Read and distribute rho(G) from file  "filename".* 
      !! (* = dat if fortran binary, * = hdf5 if HDF5)
      !! Processor "root_in_group" reads from file, distributes to
      !! all processors in the intra_group_comm communicator 
      !
      USE mp,         ONLY : mp_size, mp_rank, mp_bcast
      USE mp_wave,    ONLY : splitwf
      USE gvect,      ONLY : ngm_g
      !
!# 615 "io_base.f90"
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      !! name of file read (to which a suffix is added)
      INTEGER,          INTENT(IN) :: root_in_group
      !! root processor that reads and distributes
      INTEGER,          INTENT(IN) :: intra_group_comm
      !! rho(G) is distributed over this group of processors
      INTEGER,          INTENT(IN) :: ig_l2g(:)
      !! local-to-global indices, for machine- and mpi-independent ordering
      !! on this processor, G(ig) maps to G(ig_l2g(ig)) in global ordering
      INTEGER,          INTENT(IN) :: nspin
      !! read up to nspin components
      COMPLEX(dp),  INTENT(INOUT) :: rho(:,:)
      !! temporary check while waiting for more definitive solutions
      LOGICAL, OPTIONAL, INTENT(IN) :: gamma_only
      !! if present, don't stop in case of open error, return a nonzero value
      INTEGER, OPTIONAL, INTENT(OUT):: ier_
      !
      COMPLEX(dp), ALLOCATABLE :: rho_g(:)
      COMPLEX(dp), ALLOCATABLE :: rhoaux(:)
      REAL(dp)                 :: b1(3), b2(3), b3(3)
      INTEGER                  :: ngm, nspin_, isup, isdw
      INTEGER                  :: iun, ns, ig, ierr
      INTEGER                  :: me_in_group, nproc_in_group
      LOGICAL                  :: ionode_in_group, gamma_only_, readmill
      INTEGER                  :: ngm_g_
      INTEGER, ALLOCATABLE     :: mill_g(:,:)
      !
!# 649 "io_base.f90"
      !
      ngm  = SIZE (rho, 1)
      IF (ngm /= SIZE (ig_l2g, 1) ) &
           CALL errore('read_rhog', 'inconsistent input dimensions', 1)
      !
      iun  = 4
      ierr = 0
      IF ( PRESENT(ier_) ) ier_ = 0
      !
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
      !
      IF ( ionode_in_group ) THEN
!# 680 "io_base.f90"
         OPEN ( UNIT = iun, FILE = TRIM( filename ) // '.dat', &
              FORM = 'unformatted', STATUS = 'old', iostat = ierr )
         IF ( ierr /= 0 ) THEN
            IF ( PRESENT(ier_) ) THEN
               ier_ = ierr
               RETURN
            END IF
            ierr = 1
            GO TO 10
         END IF
         READ (iun, iostat=ierr) gamma_only_, ngm_g_, nspin_
         IF ( ierr /= 0 ) THEN
            ierr = 2
            GO TO 10
         END IF
         READ (iun, iostat=ierr) b1, b2, b3
         IF ( ierr /= 0 ) ierr = 3
!# 698 "io_base.f90"
10       CONTINUE 
      END IF
      !
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), ierr )
      CALL mp_bcast( ngm_g_, root_in_group, intra_group_comm )
      CALL mp_bcast( nspin_, root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only_, root_in_group, intra_group_comm )
      !
      IF ( nspin > nspin_ ) &
         CALL infomsg('read_rhog', 'some spin components not found')
      IF ( ngm_g < MAXVAL (ig_l2g(:)) ) &
           CALL infomsg('read_rhog', 'some G-vectors are missing, zero-padding' )
      !
      ! ... if required and if there is a mismatch between input gamma tricks
      ! ... and gamma tricks read from file: allocate and read Miller indices
      !
      readmill = PRESENT(gamma_only) 
      IF ( readmill ) readmill = ( gamma_only .NEQV. gamma_only_ ) 
      !
      IF (readmill .AND. ionode_in_group) THEN
         ALLOCATE (mill_g(3,ngm_g_))
!# 727 "io_base.f90"
         READ (iun, iostat=ierr) mill_g(1:3,1:ngm_g_)
!# 729 "io_base.f90"
      ELSE
         ALLOCATE (mill_g(1,1))
!# 732 "io_base.f90"
         ! .. skip record containing G-vector indices
         IF ( ionode_in_group) READ (iun, iostat=ierr) mill_g(1,1)
!# 735 "io_base.f90"
      END IF
      !
      CALL mp_bcast( ierr, root_in_group, intra_group_comm )
      IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
           & // TRIM( filename ), 2 )
      !
      ! ... now read, broadcast and re-order G-vector components
      ! ... of the charge density (one spin at the time to save memory)
      !
      IF ( ionode_in_group ) THEN
         ALLOCATE( rho_g(MAX(ngm_g_,ngm_g)) )
      ELSE
         ALLOCATE( rho_g( 1 ) )
      END IF
      ALLOCATE (rhoaux(ngm))
!# 757 "io_base.f90"
      !
      DO ns = 1, nspin_
         !
         IF ( ionode_in_group ) THEN
!# 766 "io_base.f90"
           READ (iun, iostat=ierr) rho_g(1:ngm_g_)
!# 768 "io_base.f90"
           IF ( ngm_g > ngm_g_) rho_g(ngm_g_+1:ngm_g) = cmplx(0.d0,0.d0, KIND = DP) 
         END IF
         CALL mp_bcast( ierr, root_in_group, intra_group_comm )
         IF ( ierr > 0 ) CALL errore ( 'read_rhog','error reading file ' &
              & // TRIM( filename ), 2+ns )
         !
         ! ... Convert charge from full G-vector to half G-vector format
         !
         IF ( readmill ) CALL charge_k_to_g (ngm_g_, rho_g, mill_g, &
              root_in_group,intra_group_comm, gamma_only)
         !
         CALL splitwf( rhoaux, rho_g, ngm, ig_l2g, me_in_group, &
              nproc_in_group, root_in_group, intra_group_comm )
         DO ig = 1, ngm
            rho(ig,ns) = rhoaux(ig)
         END DO
         ! 
      END DO
      !
!# 790 "io_base.f90"
      IF ( ionode_in_group ) CLOSE (UNIT = iun, status ='keep' )
!# 792 "io_base.f90"
      !
      DEALLOCATE( rhoaux )
      DEALLOCATE( rho_g )
      IF (ALLOCATED(mill_g))  DEALLOCATE( mill_g )
      !
      RETURN
      !
    END SUBROUTINE read_rhog
    !
    SUBROUTINE charge_k_to_g ( ngm_g_file, rho_g, mill_g_file, root_in_group, &
         intra_group_comm , this_run_is_gamma_only)
      !
      !! This routine reorders G-vectors for the charge density on global mesh
      !! from the k case to the gamma-only one.
      !
      USE io_global,     ONLY : stdout
      USE gvect,         ONLY : ngm, ngm_g, ig_l2g, mill
      USE mp,            ONLY : mp_size,mp_rank
      USE mp_wave,       ONLY : mergekg
     
      IMPLICIT NONE
      INTEGER, INTENT(in) :: intra_group_comm,root_in_group
      INTEGER, INTENT(in) :: ngm_g_file  
      !! number of g vectors found in file 
      INTEGER, INTENT(in) :: mill_g_file(:,:)
      COMPLEX(kind=DP), INTENT(inout) :: rho_g(:)!relative to k case  in input, gamma case in output
      LOGICAL, OPTIONAL, INTENT(in) :: this_run_is_gamma_only 
      INTEGER                  :: me_in_group, npr
      COMPLEX(kind=DP), ALLOCATABLE :: rho_aux(:)
      LOGICAL                  :: ionode_in_group
      INTEGER :: nproc_in_group
      INTEGER, ALLOCATABLE :: mill_g(:,:), grid(:,:,:) 
      CHARACTER(len=256)   :: mesg
      INTEGER :: ig, jg, nr1b2,nr2b2,nr3b2  
      IF ( .NOT. PRESENT (this_run_is_gamma_only) ) RETURN 
      IF ( this_run_is_gamma_only) THEN 
         call infomsg('read_rhog','Conversion: K charge Gamma charge') 
      ELSE 
         call infomsg ('read_rhog', 'Conversion: Gamma charge to K charge') 
      ENDIF 
!# 833 "io_base.f90"
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      ionode_in_group = ( me_in_group == root_in_group )
!# 837 "io_base.f90"
      IF(ionode_in_group) THEN
         allocate(rho_aux(MAX(ngm_g_file, ngm_g) ))
         allocate(mill_g(3,ngm_g))
         rho_aux(1:ngm_g_file)=rho_g(1:ngm_g_file)
      ELSE
         allocate(rho_aux(1))
         allocate(mill_g(1,1))
      ENDIF
!# 846 "io_base.f90"
      CALL mergekg( mill, mill_g, ngm, ig_l2g, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
!# 849 "io_base.f90"
      IF(ionode_in_group) THEN
         rho_g(:)= cmplx(0.d0, 0.d0,KIND = DP) 
         IF ( this_run_is_gamma_only ) THEN 
            ig = 1 
            DO jg=1,ngm_g_file
               if(  mill_g(1,ig)==mill_g_file(1,jg) .and. &
                  mill_g(2,ig)==mill_g_file(2,jg) .and. &
                  mill_g(3,ig)==mill_g_file(3,jg) ) then
                  rho_g(ig)=rho_aux(jg)
                  ig = ig + 1 
               endif
               IF ( ig .GE. ngm_g ) EXIT  
               END DO
         ELSE ! this run uses full fft mesh 
            nr1b2 = MAX(MAXVAL(ABS(mill_g(1,:))),MAXVAL(ABS(mill_g_file(1,:)))) 
            nr2b2 = MAX(MAXVAL(ABS(mill_g(2,:))),MAXVAL(ABS(mill_g_file(2,:)))) 
            nr3b2 = MAX(MAXVAL(ABS(mill_g(3,:))),MAXVAL(ABS(mill_g_file(3,:))))   
            ALLOCATE(grid(-nr1b2-1:nr1b2+1,-nr2b2-1:nr2b2+1,-nr3b2-1:nr3b2+1)) 
            grid = 10 * ngm_g_file
            !$omp do   
            DO ig = 1, ngm_g_file
               grid ( mill_g_file(1,ig), mill_g_file(2,ig),mill_g_file(3,ig))     = ig
               grid(-mill_g_file(1,ig),-mill_g_file(2,ig),-mill_g_file(3,ig))     = -ig  
            END DO 
            !$omp do private(ig) 
            DO jg =1, ngm_g
               ig = grid(mill_g(1,jg),mill_g(2,jg),mill_g(3,jg))
               IF (ig .LE. ngm_g_file) THEN  
                 IF (ig .GE. 0 ) THEN 
                   rho_g(jg) = rho_aux(ig) 
                 ELSE 
                   rho_g(jg) = CONJG(rho_aux(-ig)) 
                 END IF 
               END IF
            END DO 
            deallocate(grid) 
         END IF 
      END IF
      deallocate(rho_aux,mill_g)  
      return
 
    END SUBROUTINE charge_k_to_g
    !
  END MODULE io_base
