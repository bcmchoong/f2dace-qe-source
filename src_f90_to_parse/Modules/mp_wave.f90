!# 1 "mp_wave.f90"
!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
    MODULE mp_wave
      !
      !! MPI management of wave function related arrays.
      !
      IMPLICIT NONE
      SAVE
!# 16 "mp_wave.f90"
    CONTAINS
!# 18 "mp_wave.f90"
      SUBROUTINE mergewf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )
!# 20 "mp_wave.f90"
      !! This subroutine merges the pieces of a wave functions (pw) splitted across 
      !! processors into a total wave function (pwt) containing al the components
      !! in a pre-defined order (the same as if only one processor is used).
!# 24 "mp_wave.f90"
      USE kinds
      USE parallel_include
!# 27 "mp_wave.f90"
      IMPLICIT NONE
!# 29 "mp_wave.f90"
      COMPLEX(DP), intent(in) :: PW(:)
      !! piece of wave function
      COMPLEX(DP), intent(out) :: PWT(:)
      !! total wave function
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor ( the one that should receive the data )
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl
!# 44 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)
!# 47 "mp_wave.f90"
      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid
!# 53 "mp_wave.f90"
!
! ... Subroutine Body
!
!# 57 "mp_wave.f90"
      igwx = MAXVAL( ig_l2g(1:ngwl) )
!# 70 "mp_wave.f90"
      IF ( mpime == root .AND. igwx > SIZE( pwt ) ) &
        CALL errore(' mergewf ',' wrong size for pwt ',SIZE(pwt) )
!# 114 "mp_wave.f90"
      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
        PWT( ig_l2g(i) ) = pw(i)
      END DO
!# 125 "mp_wave.f90"
      RETURN
      END SUBROUTINE mergewf
!# 128 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
      
      SUBROUTINE mergekg ( mill, millt, ngwl, ig_l2g, mpime, nproc, root, comm )
!# 132 "mp_wave.f90"
      !! Same logic as for \(\texttt{mergewf}\), for Miller indices.
!# 134 "mp_wave.f90"
      USE kinds
      USE parallel_include
!# 137 "mp_wave.f90"
      IMPLICIT NONE
!# 139 "mp_wave.f90"
      INTEGER, intent(in) :: mill(:,:)
      !! Miller indices: distributed input
      INTEGER, intent(out):: millt(:,:)
      !! Miller indices: collected output
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl
!# 154 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
      INTEGER, ALLOCATABLE :: mill_ip(:,:)
!# 157 "mp_wave.f90"
      INTEGER :: ierr, i, ip, ngw_ip, ngw_lmax, itmp, igwx, gid
!# 163 "mp_wave.f90"
!
! ... Subroutine Body
!
!# 167 "mp_wave.f90"
      igwx = MAXVAL( ig_l2g(1:ngwl) )
!# 179 "mp_wave.f90"
      IF ( mpime == root .AND. igwx > SIZE( millt, 2 ) ) &
        CALL errore(' mergekg',' wrong size for millt ',SIZE(millt,2) )
!# 223 "mp_wave.f90"
      DO I = 1, ngwl
        ! WRITE( stdout,*) 'MW ', ig_l2g(i), i
         millt(:,ig_l2g(i) ) = mill(:,i)
      END DO
!# 234 "mp_wave.f90"
      RETURN
    END SUBROUTINE mergekg
!# 237 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 239 "mp_wave.f90"
      SUBROUTINE splitwf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )
!# 241 "mp_wave.f90"
      !! This subroutine splits a total wave function (PWT) containing al the components
      !! in a pre-defined order (the same as if only one processor is used), across 
      !! processors (PW).
!# 245 "mp_wave.f90"
      USE kinds
      USE parallel_include
      IMPLICIT NONE
!# 249 "mp_wave.f90"
      COMPLEX(DP), INTENT(OUT) :: PW(:)
      !! piece of wave function
      COMPLEX(DP), INTENT(IN) :: PWT(:)
      !! total wave function
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl
!# 264 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
      COMPLEX(DP), ALLOCATABLE :: pw_ip(:)
!# 267 "mp_wave.f90"
      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp, size_pwt
!# 273 "mp_wave.f90"
!
! ... Subroutine Body
!
!# 277 "mp_wave.f90"
      igwx = MAXVAL( ig_l2g(1:ngwl) )
!# 290 "mp_wave.f90"
      IF ( mpime == root .AND. igwx > SIZE(pwt )) &
        CALL errore (' splitwf ',' wrong size for pwt', SIZE(pwt) )
!# 327 "mp_wave.f90"
      DO I = 1, ngwl
           pw(i) = pwt( ig_l2g(i) ) 
      END DO
!# 337 "mp_wave.f90"
      RETURN
      END SUBROUTINE splitwf
!# 340 "mp_wave.f90"
      !=----------------------------------------------------------------------------=!
!# 342 "mp_wave.f90"
      SUBROUTINE splitkg ( mill, millt, ngwl, ig_l2g, mpime, nproc, root, comm )
!# 344 "mp_wave.f90"
      !! Same logic as for \(\texttt{splitwf}\), for Miller indices.
!# 346 "mp_wave.f90"
      USE kinds
      USE parallel_include
      IMPLICIT NONE
!# 350 "mp_wave.f90"
      INTEGER, INTENT(OUT):: mill(:,:)
      !! Miller indices: distributed output
      INTEGER, INTENT(IN) :: millt(:,:)
      !! Miller indices: collected input
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ig_l2g(:)
      INTEGER, INTENT(IN) :: ngwl
!# 365 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
      INTEGER, ALLOCATABLE :: mill_ip(:,:)
!# 368 "mp_wave.f90"
      INTEGER ierr, i, ngw_ip, ip, ngw_lmax, gid, igwx, itmp
!# 374 "mp_wave.f90"
!
! ... Subroutine Body
!
!# 378 "mp_wave.f90"
      igwx = MAXVAL( ig_l2g(1:ngwl) )
!# 391 "mp_wave.f90"
      IF ( mpime == root .AND. igwx > SIZE( millt, 2 ) ) &
        CALL errore(' splitkg ',' wrong size for millt ',SIZE(millt,2) )
!# 428 "mp_wave.f90"
      DO I = 1, ngwl
         mill(:,i) = millt(:,ig_l2g(i)) 
      END DO
!# 438 "mp_wave.f90"
      RETURN
    END SUBROUTINE splitkg
!# 441 "mp_wave.f90"
      SUBROUTINE mergeig(igl, igtot, ngl, mpime, nproc, root, comm)
!# 443 "mp_wave.f90"
      !! This subroutine merges the pieces of a vector splitted across 
      !! processors into a total vector (igtot) containing al the components
      !! in a pre-defined order (the same as if only one processor is used).
!# 447 "mp_wave.f90"
      USE kinds
      USE parallel_include
!# 450 "mp_wave.f90"
      IMPLICIT NONE
!# 452 "mp_wave.f90"
      INTEGER, intent(in)  :: igl(:)
      !! piece of splitted vector
      INTEGER, intent(out) :: igtot(:)
      !! total vector
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ngl
!# 466 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
!# 468 "mp_wave.f90"
      INTEGER :: ierr, i, ip, ng_ip, ng_lmax, ng_g, gid, igs
!# 527 "mp_wave.f90"
      igtot( 1:ngl ) = igl( 1:ngl )
!# 535 "mp_wave.f90"
      RETURN
      END SUBROUTINE mergeig
!# 538 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 540 "mp_wave.f90"
      SUBROUTINE splitig(igl, igtot, ngl, mpime, nproc, root, comm)
!# 542 "mp_wave.f90"
      !! This subroutine splits a replicated vector (\(\text{igtot}\)) stored on
      !! the \(\text{root}\) proc across processors (\(\text{igl}\)).
!# 545 "mp_wave.f90"
      USE kinds
      USE parallel_include
      IMPLICIT NONE
!# 549 "mp_wave.f90"
      INTEGER, INTENT(OUT) :: igl(:)
      !! vector splitted across procs
      INTEGER, INTENT(IN)  :: igtot(:)
      !! replicated vector on root proc
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: root
      !! root processor
      INTEGER, INTENT(IN) :: comm
      !! communicator
      INTEGER, INTENT(IN) :: ngl
!# 563 "mp_wave.f90"
      INTEGER ierr, i, ng_ip, ip, ng_lmax, ng_g, gid, igs
!# 569 "mp_wave.f90"
      INTEGER, ALLOCATABLE :: ig_ip(:)
!# 625 "mp_wave.f90"
      igl( 1:ngl ) = igtot( 1:ngl )
!# 633 "mp_wave.f90"
      RETURN
      END SUBROUTINE splitig
!# 636 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 638 "mp_wave.f90"
   SUBROUTINE pwscatter( c, ctmp, ngw, indi_l, sour_indi, dest_indi, &
      n_indi_rcv, n_indi_snd, icntix, mpime, nproc, group )
!# 641 "mp_wave.f90"
      USE kinds
      USE parallel_include
!# 644 "mp_wave.f90"
      implicit none
!# 646 "mp_wave.f90"
      integer :: indi_l(:)
      !! list of G-vec index to be exchanged
      integer :: sour_indi(:)
      !! the list of source processors
      integer :: dest_indi(:)
      !! the list of destination processors
      integer :: n_indi_rcv
      !! number of G-vectors to be received
      integer :: n_indi_snd
      !! number of G-vectors to be sent
      integer :: icntix
      !! total number of G-vec to be exchanged
      INTEGER, INTENT(IN) :: mpime
      !! index of the calling processor ( starting from 0 )
      INTEGER, INTENT(IN) :: nproc
      !! number of processors
      INTEGER, INTENT(IN) :: group
      
      COMPLEX(DP) :: c(:)
      COMPLEX(DP) :: ctmp(:)
      integer  ::  ngw
!# 668 "mp_wave.f90"
      integer :: ig, icsize
      INTEGER :: me, idest, isour, ierr
!# 671 "mp_wave.f90"
      COMPLEX(DP), ALLOCATABLE :: my_buffer( : )
      COMPLEX(DP), ALLOCATABLE :: mp_snd_buffer( : )
      COMPLEX(DP), ALLOCATABLE :: mp_rcv_buffer( : )
      INTEGER, ALLOCATABLE :: ibuf(:)
!# 676 "mp_wave.f90"
      !
      ! ... SUBROUTINE BODY
      !
!# 680 "mp_wave.f90"
      me = mpime + 1
!# 682 "mp_wave.f90"
      if( icntix .lt. 1 ) then
        icsize = 1
      else
        icsize = icntix
      endif
!# 688 "mp_wave.f90"
      ALLOCATE( mp_snd_buffer( icsize * nproc ) )
      ALLOCATE( mp_rcv_buffer( icsize * nproc ) )
      ALLOCATE( my_buffer( ngw ) )
      ALLOCATE( ibuf( nproc ) )
      ctmp = ( 0.0_DP, 0.0_DP )
!# 694 "mp_wave.f90"
      ! WRITE( stdout,*) 'D: ', nproc, mpime, group
!# 696 "mp_wave.f90"
      ibuf = 0
      DO IG = 1, n_indi_snd
        idest = dest_indi(ig)
        ibuf(idest) = ibuf(idest) + 1;
        if(idest .ne. me) then
          mp_snd_buffer( ibuf(idest) + (idest-1)*icsize ) = C( indi_l( ig ) )
        else
          my_buffer(ibuf(idest)) = C(indi_l(ig))
        end if
      end do
!# 713 "mp_wave.f90"
      CALL errore(' pwscatter ',' no communication protocol ',0)
!# 717 "mp_wave.f90"
      ibuf = 0
      DO IG = 1, n_indi_rcv
        isour = sour_indi(ig)
        if(isour.gt.0 .and. isour.ne.me) then
          ibuf(isour) = ibuf(isour) + 1
          CTMP(ig) = mp_rcv_buffer(ibuf(isour) + (isour-1)*icsize)
        else if(isour.gt.0) then
          ibuf(isour) = ibuf(isour) + 1
          CTMP(ig) = my_buffer(ibuf(isour))
        else
          CTMP(ig) = (0.0_DP,0.0_DP)
        end if
      end do
!# 731 "mp_wave.f90"
      DEALLOCATE( mp_snd_buffer )
      DEALLOCATE( mp_rcv_buffer )
      DEALLOCATE( my_buffer )
      DEALLOCATE( ibuf )
!# 736 "mp_wave.f90"
      RETURN
    END SUBROUTINE pwscatter
!# 742 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 744 "mp_wave.f90"
SUBROUTINE redistwf( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !! Redistribute wave function.
   !
   USE kinds
   USE parallel_include
!# 751 "mp_wave.f90"
   implicit none
!# 753 "mp_wave.f90"
   COMPLEX(DP) :: c_dist_pw(:,:)
   !! the wave functions with plane waves distributed over processors 
   COMPLEX(DP) :: c_dist_st(:,:)
   !! the wave functions with electronic states distributed over processors 
   INTEGER, INTENT(IN) :: npw_p(:)
   !! the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)
   !! the number of states on each processor
   INTEGER, INTENT(IN) :: comm
   !! group communicator
   INTEGER, INTENT(IN) :: idir
   !! direction of the redistribution:  
   !! \(\text{idir}>0\):  \(\text{c_dist_pw}\rightarrow\text{c_dist_st}\)  
   !! \(\text{idir}<0\):  \(\text{c_dist_pw}\leftarrow\text{c_dist_st}\)
!# 768 "mp_wave.f90"
   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww, ii
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   COMPLEX(DP), ALLOCATABLE :: ctmp( : )
!# 857 "mp_wave.f90"
   RETURN
END SUBROUTINE redistwf
!# 860 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 862 "mp_wave.f90"
SUBROUTINE redistwfr( c_dist_pw, c_dist_st, npw_p, nst_p, comm, idir )
   !
   !!  Redistribute wave function.
   !
   USE kinds
   USE parallel_include
!# 869 "mp_wave.f90"
   implicit none
!# 871 "mp_wave.f90"
   REAL(DP) :: c_dist_pw(:,:)
   !! the wave functions with plane waves distributed over processors 
   REAL(DP) :: c_dist_st(:,:)
   !! the wave functions with electronic states distributed over processors 
   INTEGER, INTENT(IN) :: npw_p(:)
   !! the number of plane wave on each processor
   INTEGER, INTENT(IN) :: nst_p(:)
   !! the number of states on each processor
   INTEGER, INTENT(IN) :: comm
   !! group communicator
   INTEGER, INTENT(IN) :: idir
   !! direction of the redistribution:  
   !! \(\text{idir}>0\):  \(\text{c_dist_pw}\rightarrow\text{c_dist_st}\)  
   !! \(\text{idir}<0\):  \(\text{c_dist_pw}\leftarrow\text{c_dist_st}\)
!# 886 "mp_wave.f90"
   INTEGER :: mpime, nproc, ierr, npw_t, nst_t, proc, i, j, ngpww
   INTEGER, ALLOCATABLE :: rdispls(:),  recvcount(:)
   INTEGER, ALLOCATABLE :: sendcount(:),  sdispls(:)
   REAL(DP), ALLOCATABLE :: ctmp( : )
!# 973 "mp_wave.f90"
   RETURN
END SUBROUTINE redistwfr
!# 976 "mp_wave.f90"
!=----------------------------------------------------------------------------=!
!# 978 "mp_wave.f90"
    END MODULE mp_wave
