!# 1 "distools.f90"
!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
!
SUBROUTINE block_distribute( nat, me, nproc, ia_s, ia_e, mykey )
    !
    ! Distribute "nat" objects (.e.g atoms) among "nproc" processors
    ! Atoms "ia_s" to "ia_e" are assigned to this ("me") processor
    ! If nproc > nat, atoms are assigned more than once to processors,
    ! "mykey" labels how many times the same atom appears
    !
    INTEGER, INTENT(IN)  :: nat, me, nproc  
    INTEGER, INTENT(OUT) :: ia_s, ia_e, mykey
    INTEGER :: na_loc, r, nproc_ia
!# 20 "distools.f90"
    INTEGER, EXTERNAL :: ldim_block, gind_block
!# 22 "distools.f90"
    ! compute how many processors we have for a given atom
    !
    nproc_ia = nproc / nat
    !
    IF( nproc_ia == 0 ) THEN
       !       
       ! here we have less than one processor per atom
       !       
       mykey  = 0
       na_loc = ldim_block( nat, nproc, me)
       ia_s   = gind_block( 1, nat, nproc, me )
       ia_e   = ia_s + na_loc - 1
       !       
    ELSE 
       !       
       ! here we have more than one proc per atom
       !
       r = MOD( nproc, nat )
       !
       IF( me < (nproc_ia + 1)*r ) THEN
          ! processors that do the work, more procs work on a single atom
          ia_s  = me/(nproc_ia + 1) + 1
          mykey = MOD( me, nproc_ia + 1 )
       ELSE
          ia_s  = ( me - (nproc_ia + 1)*r ) / nproc_ia + 1 + r
          mykey = MOD( me - (nproc_ia + 1)*r , nproc_ia )
       END IF
       !
       ia_e = ia_s
       !
    END IF
!# 54 "distools.f90"
    RETURN
!# 56 "distools.f90"
END SUBROUTINE
!
!
FUNCTION block_size(ia, nat, nproc) RESULT (res) 
   !! counts how many procs have the same ia   
   IMPLICIT NONE 
   !
   INTEGER   :: res
   INTEGER,INTENT(IN) :: ia, nat, nproc 
   !            
   res  = nproc/nat
   IF (ia <= MOD(nproc,nat)) res = res + 1   
END FUNCTION block_size  
!
!
SUBROUTINE GRID2D_DIMS( grid_shape, nproc, nprow, npcol )
   !
   ! This subroutine factorizes the number of processors (NPROC)
   ! into NPROW and NPCOL according to the shape
   !
   !    Written by Carlo Cavazzoni
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: grid_shape
   INTEGER, INTENT(IN)  :: nproc
   INTEGER, INTENT(OUT) :: nprow, npcol
   INTEGER :: sqrtnp, i
   !
   sqrtnp = INT( SQRT( REAL( nproc ) + 0.1 ) )
   !
   IF( grid_shape == 'S' ) THEN
      ! Square grid
      nprow = sqrtnp
      npcol = sqrtnp
   ELSE
      ! Rectangular grid
      DO i = 1, sqrtnp + 1
         IF( MOD( nproc, i ) == 0 ) nprow = i
      end do
      npcol = nproc / nprow
   END IF
   RETURN
END SUBROUTINE
!# 100 "distools.f90"
SUBROUTINE GRID2D_COORDS( order, rank, nprow, npcol, row, col )
   !
   !  this subroutine computes the cartesian coordinates "row" and "col"
   !  of the processor whose MPI task id is "rank". 
   !  Note that if the rank is larger that the grid size
   !  all processors whose MPI task id is greather or equal 
   !  than nprow * npcol are placed on the diagonal extension of the grid itself
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: order
   INTEGER, INTENT(IN)  ::  rank          ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol  ! dimensions of the processor grid
   INTEGER, INTENT(OUT) ::  row, col
   IF( rank >= 0 .AND. rank < nprow * npcol ) THEN
      IF( order == 'C' .OR. order == 'c' ) THEN
         !  grid in COLUMN MAJOR ORDER
         row = MOD( rank, nprow )
         col = rank / nprow
      ELSE
         !  grid in ROW MAJOR ORDER
         row = rank / npcol
         col = MOD( rank, npcol )
      END IF
   ELSE
      row = rank
      col = rank
   END IF
   RETURN
END SUBROUTINE
!# 130 "distools.f90"
SUBROUTINE GRID2D_RANK( order, nprow, npcol, row, col, rank )
   !
   !  this subroutine computes the processor MPI task id "rank" of the processor  
   !  whose cartesian coordinate are "row" and "col".
   !  Note that the subroutine assumes cyclic indexing ( row = nprow = 0 )
   !
   IMPLICIT NONE
   CHARACTER, INTENT(IN) :: order
   INTEGER, INTENT(OUT) ::  rank         ! process index starting from 0
   INTEGER, INTENT(IN)  ::  nprow, npcol ! dimensions of the processor grid
   INTEGER, INTENT(IN)  ::  row, col
   
   IF( order == 'C' .OR. order == 'c' ) THEN
     !  grid in COLUMN MAJOR ORDER
     rank = MOD( row + nprow, nprow ) + MOD( col + npcol, npcol ) * nprow
   ELSE
     !  grid in ROW MAJOR ORDER
     rank = MOD( col + npcol, npcol ) + MOD( row + nprow, nprow ) * npcol
   END IF
   !
   RETURN
END SUBROUTINE
!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!# 161 "distools.f90"
       INTEGER FUNCTION ldim_cyclic(gdim, np, me)
!# 163 "distools.f90"
!   gdim = global dimension of distributed array
!   np   = number of processors
!   me   = index of the calling processor (starting from 0)
!  
!   this function returns the number of elements of the distributed array
!   stored in the local memory of the processor "me" for a cyclic 
!   data distribution.
!   Example of the cyclic distribution of a 10 elements array on 4 processors
!   array elements  |  PEs
!    a(1)           |   0
!    a(2)           |   1
!    a(3)           |   2
!    a(4)           |   3
!    a(5)           |   0
!    a(6)           |   1
!    a(7)           |   2
!    a(8)           |   3
!    a(9)           |   0
!    a(10)          |   1
!# 183 "distools.f90"
       IMPLICIT NONE
       INTEGER :: gdim, np, me, r, q
!# 186 "distools.f90"
       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_cyclic: arg no. 3 out of range '
         STOP
       END IF
!# 191 "distools.f90"
       q = INT(gdim / np)
       r = MOD(gdim, np)
!# 194 "distools.f90"
       IF( me .LT. r ) THEN
! ...    if my index is less than the reminder I got an extra element
         ldim_cyclic = q+1
       ELSE
         ldim_cyclic = q
       END IF
 
       RETURN
       END FUNCTION ldim_cyclic
!# 204 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 206 "distools.f90"
       INTEGER FUNCTION ldim_block(gdim, np, me)
!# 208 "distools.f90"
!   gdim = global dimension of distributed array
!   np   = number of processors
!   me   = index of the calling processor (starting from 0)
!  
!   this function returns the number of elements of the distributed array
!   stored in the local memory of the processor "me" for a balanced block 
!   data distribution, with the larger block on the lower index processors.
!   Example of the block distribution of 10 elements array a on 4 processors
!   array elements  |  PEs
!    a(1)           |   0
!    a(2)           |   0
!    a(3)           |   0
!    a(4)           |   1
!    a(5)           |   1
!    a(6)           |   1
!    a(7)           |   2
!    a(8)           |   2
!    a(9)           |   3
!    a(10)          |   3
!# 228 "distools.f90"
       IMPLICIT NONE
       INTEGER :: gdim, np, me, r, q
!# 231 "distools.f90"
       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_block: arg no. 3 out of range '
         STOP
       END IF
!# 236 "distools.f90"
       q = INT(gdim / np)
       r = MOD(gdim, np)
!# 239 "distools.f90"
       IF( me .LT. r ) THEN
! ...    if my index is less than the reminder I got an extra element
         ldim_block = q+1
       ELSE
         ldim_block = q
       END IF
 
       RETURN
       END FUNCTION ldim_block
!# 249 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 251 "distools.f90"
       INTEGER FUNCTION ldim_block_sca( gdim, np, me )
!# 253 "distools.f90"
!   gdim = global dimension of distributed array
!   np   = number of processors
!   me   = index of the calling processor (starting from 0)
!  
!   this function returns the number of elements of the distributed array
!   stored in the local memory of the processor "me" for equal block
!   data distribution, all block have the same size but the last one.
!   Example of the block distribution of 10 elements array a on 4 processors
!   array elements  |  PEs
!    a(1)           |   0
!    a(2)           |   0
!    a(3)           |   0
!    a(4)           |   1
!    a(5)           |   1
!    a(6)           |   1
!    a(7)           |   2
!    a(8)           |   2
!    a(9)           |   2
!    a(10)          |   3
!# 273 "distools.f90"
       IMPLICIT NONE
       INTEGER :: gdim, np, me, nb
!# 276 "distools.f90"
       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_block: arg no. 3 out of range '
         STOP
       END IF
!# 281 "distools.f90"
       nb = INT( gdim / np )
       IF( MOD( gdim,  np ) /= 0 ) THEN
         nb = nb+1
         ! ... last processor take the rest
         IF( me == ( np - 1 ) ) nb = gdim - (np-1)*nb
       END IF
!# 288 "distools.f90"
       ldim_block_sca = nb
!# 291 "distools.f90"
       RETURN
       END FUNCTION ldim_block_sca
!# 294 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 296 "distools.f90"
       INTEGER FUNCTION lind_cyclic(ig, nx, np, me)
!
!   INPUT :
!      ig  global index of the x dimension of array element
!      nx  dimension of the global array
!      np  number of processor in the x dimension of the processors grid
!      me  index of the local processor in the processor grid
!                (starting from zero)
!
!   OUTPUT :
!
!      lind_cyclic return the local index corresponding to the
!      global index "ig" for a cyclic distribution
!   
!# 311 "distools.f90"
       IMPLICIT NONE
!# 313 "distools.f90"
       INTEGER :: ig, nx, np, me
!# 315 "distools.f90"
       lind_cyclic = (ig-1)/np + 1
!# 317 "distools.f90"
       RETURN 
       END FUNCTION lind_cyclic
!# 321 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 324 "distools.f90"
      INTEGER FUNCTION gind_cyclic( lind, n, np, me )
!# 326 "distools.f90"
!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.
!
!# 336 "distools.f90"
            INTEGER, INTENT(IN) :: lind, n, me, np
            INTEGER r, q
!# 339 "distools.f90"
            gind_cyclic = (lind-1) * np + me + 1
!# 341 "distools.f90"
            RETURN
      END FUNCTION gind_cyclic
!# 345 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 348 "distools.f90"
      INTEGER FUNCTION gind_block( lind, n, np, me )
!# 350 "distools.f90"
!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.
!# 360 "distools.f90"
            INTEGER, INTENT(IN) :: lind, n, me, np
            INTEGER r, q
!# 363 "distools.f90"
              q = INT(n/np)
              r = MOD(n,np)
              IF( me < r ) THEN
                gind_block = (Q+1)*me + lind
              ELSE
                gind_block = Q*me + R + lind
              END IF
!# 371 "distools.f90"
         RETURN
      END FUNCTION gind_block
!# 374 "distools.f90"
!=----------------------------------------------------------------------------=!
!# 376 "distools.f90"
      INTEGER FUNCTION gind_block_sca( lind, n, np, me )
!# 378 "distools.f90"
!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  np        The total number processes over which the distributed
!            matrix is distributed.
!# 388 "distools.f90"
       INTEGER, INTENT(IN) :: lind, n, me, np
       INTEGER nb
!# 391 "distools.f90"
       IF( me >= np .OR. me < 0 ) THEN
         WRITE(6,*) ' ** ldim_block: arg no. 3 out of range '
         STOP
       END IF
!# 396 "distools.f90"
       nb = INT( n / np )
       IF( MOD( n,  np ) /= 0 ) nb = nb+1
!# 399 "distools.f90"
       gind_block_sca = lind + me * nb
!# 401 "distools.f90"
       RETURN
!# 403 "distools.f90"
    END FUNCTION gind_block_sca
