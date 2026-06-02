!# 1 "zhpev_drv.f90"
!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE zhpev_module
!# 10 "zhpev_drv.f90"
   IMPLICIT NONE
   SAVE
!# 13 "zhpev_drv.f90"
CONTAINS
   !
   !-------------------------------------------------------------------------
   SUBROUTINE pzhptrd( n, nrl, ap, lda, d, e, tau, nproc, me, comm )
     !-------------------------------------------------------------------------
      !
      !  Parallel MPI version of the LAPACK routine ZHPTRD
      !
      !     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
      !     Dicember 12, 1999
      !
      !  REFERENCES :
      !
      !     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
      !     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
      !     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
      !
      !     PARALLEL NUMERICAL ALGORITHMS,
      !     T.L. FREEMAN AND C.PHILLIPS,
      !     PRENTICE HALL INTERNATIONAL (1992).
      !
      !     LAPACK routine (version 2.0) --
      !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
      !     Courant Institute, Argonne National Lab, and Rice University
      !
!# 39 "zhpev_drv.f90"
      USE laxlib_parallel_include
!# 41 "zhpev_drv.f90"
      IMPLICIT NONE
!# 43 "zhpev_drv.f90"
      include 'laxlib_kinds.fh'
!# 45 "zhpev_drv.f90"
!     .. __SCALAR Arguments ..
      INTEGER            LDA, N, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      REAL(DP)             D( * ), E( * )
      COMPLEX(DP)         AP(LDA, * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  PZHPTRD reduces a complex Hermitian distributed matrix AP  to
!  real symmetric tridiagonal form T by a unitary similarity
!  transformation: Q**H * A * Q = T.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  AP      (input/output) COMPLEX(DP) array, dimension (LDA,N)
!          On entry, the Hermitian matrix AP.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1.
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..
!# 83 "zhpev_drv.f90"
!          On exit, the diagonal and first subdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements below the first
!          subdiagonal, with the array TAU, represent the unitary
!          matrix Q as a product of elementary reflectors;
!
!  LDA     (input) INTEGER
!          Leading dimension of the local matrix AP, LDA > NRL
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = AP(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i+1,i)
!
!  TAU     (output) COMPLEX(DP) array, dimension (N-1)
!          The __SCALAR factors of the elementary reflectors (see Further
!          Details).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )
!# 110 "zhpev_drv.f90"
!
!  Further Details
!  ===============
!
!  the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a complex __SCALAR, and v is a complex vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in AP,
!  overwriting A(i+2:n,i), and tau is stored in TAU(i).
!
!  =====================================================================
!
!     .. Parameters ..
!# 131 "zhpev_drv.f90"
      COMPLEX(DP)  ONE, ZERO, HALF
      PARAMETER   ( ONE = ( 1.0_DP, 0.0_DP ),ZERO = ( 0.0_DP, 0.0_DP ),  &
     &             HALF = ( 0.5_DP, 0.0_DP ) )
      REAL(DP)      RONE, RZERO
      PARAMETER   ( RONE = 1.0_DP, RZERO = 0.0_DP )
!# 137 "zhpev_drv.f90"
      INTEGER QI
      INTEGER IL(N+1)
      INTEGER OW(N+1)
      COMPLEX(DP) CTMP
      COMPLEX(DP) CTMPV(N+1)
      COMPLEX(DP) TAUL(N+1)
      COMPLEX(DP) APKI(N+1)
      REAL(DP)     TMP
      REAL(DP)     TMPV(N+1)
!# 147 "zhpev_drv.f90"
!     ..
!     .. Local __SCALARs ..
      INTEGER            J, I, I1, K, I2, NI1, JL
      INTEGER            KL, J1
      COMPLEX(DP)         ALPHA, TAUI
      INTEGER            KNT, IERR
      REAL(DP)             ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Subroutines ..
      EXTERNAL           zaxpy
      EXTERNAL           zdscal, zscal
!     ..
!     .. External Functions ..
      ! some compiler don't like complex functions
      !COMPLEX(DP)         zdotc
      !EXTERNAL           zdotc
      !COMPLEX(DP)         ZLADIV
      !EXTERNAL            ZLADIV
      REAL(DP)             DLAMCH, DLAPY3, DZNRM2
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DABS, DBLE, AIMAG, SIGN
!
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF(N.LE.0) THEN
        RETURN
      END IF
!# 179 "zhpev_drv.f90"
      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC)
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO
!
!        Reduce the lower triangle of A.
!
         IF (OW(1).EQ.ME) THEN
           AP( IL(1), 1 ) = DBLE( AP( IL(1), 1 ) )
         END IF
!# 195 "zhpev_drv.f90"
         DO I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            IF (OW(I+1).EQ.ME) THEN
              ALPHA = AP( IL(I+1), I )
            END IF
!# 209 "zhpev_drv.f90"
            IF( (N-I).LE.0 ) THEN
              TAUI = RZERO
            ELSE
              IF(OW(I+2).EQ.ME) THEN
                I2 = IL(I+2)
              ELSE
                I2 = IL(I+2) + 1          ! I+2
              ENDIF
              NI1 = NRL - I2 + 1          ! N-I-1
!# 219 "zhpev_drv.f90"
              IF((N-I-1).GT.0) THEN
                IF( NI1 .GT. 0 ) THEN
                   XNORM = DZNRM2( NI1, AP( I2, I ), 1 )
                ELSE
                   XNORM = 0.0_DP
                END IF
!# 231 "zhpev_drv.f90"
              ELSE
                XNORM = 0.0_DP
              ENDIF
!# 235 "zhpev_drv.f90"
              ALPHR = DBLE( ALPHA )
              ALPHI = AIMAG( ALPHA )
              IF( XNORM.EQ.RZERO .AND. ALPHI.EQ.RZERO ) THEN
                TAUI = RZERO
              ELSE
                BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
                SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
                RSAFMN = RONE / SAFMIN
                IF( DABS( BETA ).LT.SAFMIN ) THEN
                  KNT = 0
   10             CONTINUE
                  KNT = KNT + 1
!# 248 "zhpev_drv.f90"
                  IF(NI1.GT.0) THEN
                    CALL zdscal( NI1, RSAFMN, AP( I2, I ), 1 )
                  ENDIF
!# 252 "zhpev_drv.f90"
                  BETA = BETA*RSAFMN
                  ALPHI = ALPHI*RSAFMN
                  ALPHR = ALPHR*RSAFMN
                  IF( DABS( BETA ).LT.SAFMIN ) GO TO 10
!# 257 "zhpev_drv.f90"
                  IF((N-I-1).GT.0) THEN
                    XNORM = DZNRM2( NI1, AP( I2, I ), 1 )
!# 265 "zhpev_drv.f90"
                  ELSE
                    XNORM = 0.0_DP
                  ENDIF
!# 269 "zhpev_drv.f90"
                  ALPHA = CMPLX( ALPHR, ALPHI, KIND=DP )
                  BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
                  TAUI = CMPLX( (BETA-ALPHR)/BETA, -ALPHI/BETA, KIND=DP )
                  ! next line yields problems on some compilers
                  ! ALPHA = ZLADIV( ONE, ALPHA-BETA )
                  ALPHA = ONE / (ALPHA-BETA)
!# 276 "zhpev_drv.f90"
                  IF(NI1.GT.0) THEN
                    CALL zscal( NI1, ALPHA, AP( I2, I ), 1 )
                  ENDIF
!# 280 "zhpev_drv.f90"
                  ALPHA = BETA
                  DO J = 1, KNT
                    ALPHA = ALPHA*SAFMIN
                  END DO
!# 285 "zhpev_drv.f90"
                ELSE
!# 287 "zhpev_drv.f90"
                  TAUI = CMPLX( (BETA-ALPHR)/BETA, -ALPHI/BETA, KIND=DP )
                  ! next line yields problems on some compilers
                  ! ALPHA = ZLADIV( ONE, ALPHA-BETA )
                  ALPHA = ONE / (ALPHA-BETA)
!# 292 "zhpev_drv.f90"
                  IF(NI1.GT.0) THEN
                    CALL zscal( NI1, ALPHA, AP( I2, I ), 1 )
                  ENDIF
!# 296 "zhpev_drv.f90"
                  ALPHA = BETA
                END IF
              END IF
            ENDIF
!
            E( I ) = ALPHA
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               ! ... AP( I+1, I ) = ONE
               IF (OW(I+1).EQ.ME) THEN
                 AP( IL(I+1), I ) = ONE
               END IF
!
!              Compute  y := tau * A * v  storing y in TAU(i:n-1)
!
!# 315 "zhpev_drv.f90"
               ! ... broadcast A(K,I)
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+2
               ENDIF
!# 333 "zhpev_drv.f90"
               DO J = I+1,N
                 APKI(J) = AP(J,I)
               ENDDO
!# 337 "zhpev_drv.f90"
               DO J = I+1, N+1
                 TAU(J-1) = ZERO
               END DO
               DO JL = I1, NRL
                 J = ME + (JL-1)*NPROC + 1
                 TAU(J-1) = ZERO
                 DO K = I+1, J
                   TAU(J-1) = TAU(J-1) + TAUI * AP(JL,K) * APKI(K)
                 END DO
               END DO
               DO J = I+1, N
                 IF(OW(J+1).EQ.ME) THEN
                   J1 = IL(J+1)
                 ELSE
                   J1 = IL(J+1) + 1          ! I+2
                 ENDIF
                 DO KL = J1, NRL
                   K = ME + (KL-1)*NPROC + 1
                   TAU(J-1) = TAU(J-1) + TAUI * CONJG(AP(KL,J)) * APKI(K)
                 END DO
               END DO
!# 365 "zhpev_drv.f90"
!
!              Compute  w := y - 1/2 * tau * (y'*v) * v
!
               ! ... ALPHA = -HALF*TAUI*zdotc(N-I,TAU(I),1,AP(I+1,I),1)
!# 370 "zhpev_drv.f90"
               JL = 1
               DO J = I, N
                 IF(OW(J+1).EQ.ME) THEN
                   TAUL(JL) = TAU(J)
                   JL = JL + 1
                 END IF
               END DO
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+1
               ENDIF
               NI1 = NRL - I1 + 1          ! N-I
               IF ( NI1 > 0 ) THEN
                  ! next line yields problems on some compilers
                  !ALPHA = -HALF*TAUI*zdotc(NI1,TAUL(1),1,AP(I1,I),1)
                  ALPHA = -HALF*TAUI*dot_product(TAUL(1:NI1),AP(I1:I1+NI1-1,I))
               ELSE
                  ALPHA = 0.0_DP
               END IF
!# 411 "zhpev_drv.f90"
               CALL zaxpy(N-I,ALPHA,AP(I+1,I),1,TAU(I),1)
!# 414 "zhpev_drv.f90"
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               ! ... broadcast A(K,I)
               IF(OW(I+1).EQ.ME) THEN
                 I1 = IL(I+1)
               ELSE
                 I1 = IL(I+1) + 1          ! I+2
               ENDIF
!# 436 "zhpev_drv.f90"
               DO J = I+1, N
                 APKI(J) = AP(J,I)
               END DO
!# 441 "zhpev_drv.f90"
               DO K = I+1,N
                 DO JL = I1,NRL
                   J = ME + (JL-1)*NPROC + 1
                   AP(JL,K) = AP(JL,K) - ONE * AP(JL,I) * CONJG(TAU(K-1)) - &
     &                CONJG(ONE) * TAU(J-1) * CONJG(APKI(K))
                 END DO
               END DO
!
            END IF
            IF(OW(I+1).EQ.ME) THEN
              AP(IL(I+1),I) = E( I )
            END IF
            IF(OW(I).EQ.ME) THEN
              D( I ) = DBLE(AP( IL(I),I ))
            END IF
!# 460 "zhpev_drv.f90"
            TAU( I ) = TAUI
         END DO
         IF(OW(I).EQ.ME) THEN
            D( N ) = DBLE(AP( IL(I),I ))
         END IF
!# 469 "zhpev_drv.f90"
!
      RETURN
!# 472 "zhpev_drv.f90"
!
!     End of ZHPTRD
!
      END SUBROUTINE pzhptrd
!# 477 "zhpev_drv.f90"
!==----------------------------------------------==!
!# 479 "zhpev_drv.f90"
   SUBROUTINE pzupgtr( n, nrl, ap, lda, tau, q, ldq, nproc, me, comm)
!# 481 "zhpev_drv.f90"
!
!  Parallel MPI version of the LAPACK routine ZUPGTR
!
!     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
!     Dicember 12, 1999
!
!  REFERENCES :
!
!     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
!     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
!     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
!
!     PARALLEL NUMERICAL ALGORITHMS,
!     T.L. FREEMAN AND C.PHILLIPS,
!     PRENTICE HALL INTERNATIONAL (1992).
!
!     LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!# 501 "zhpev_drv.f90"
      USE laxlib_parallel_include
!# 503 "zhpev_drv.f90"
      IMPLICIT NONE
!# 505 "zhpev_drv.f90"
      include 'laxlib_kinds.fh'
!
!     .. __SCALAR Arguments ..
!# 509 "zhpev_drv.f90"
      INTEGER            INFO, LDQ, N, LDA, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      COMPLEX(DP)         AP(LDA, * ), Q( LDQ, * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  PZUPGTR generates a complex unitary matrix Q which is defined as the
!  product of n-1 elementary reflectors H(i) of order n, as returned by
!  PZHPTRD :
!
!  Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  AP      (input) COMPLEX(DP) array, dimension (LDA,N)
!          The vectors which define the elementary reflectors, as
!          returned by PZHPTRD.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1.
!              Example for NPROC = 4 :
!              ROW | PE
!              1   | 0
!              2   | 1
!              3   | 2
!              4   | 3
!              5   | 0
!              6   | 1
!              ..  | ..
!
!  LDA     (input) INTEGER
!          Leading dimension of the local matrix AP, LDA > NRL
!
!  TAU     (input) COMPLEX(DP) array, dimension (N-1)
!          TAU(i) must contain the __SCALAR factor of the elementary
!          reflector H(i), as returned by PZHPTRD.
!
!  Q       (output) COMPLEX(DP) array, dimension (LDQ,N)
!          The N-by-N unitary matrix Q.
!          The rows of the matrix are distributed among processors
!          in the same way of the matrix AP
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q. LDQ >= max(1,NRL).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )
!
!  =====================================================================
!
!     .. Parameters ..
!# 573 "zhpev_drv.f90"
      COMPLEX(DP)         ONE, ZERO
      PARAMETER          ( ONE = (1.0_DP,0.0_DP), ZERO = (0.0_DP,0.0_DP) )
!# 576 "zhpev_drv.f90"
      !  change the following parameters to tune the performances
      !
      INTEGER, PARAMETER :: opt_zgemv = 40
      INTEGER, PARAMETER :: opt_zgerc = 40
!# 581 "zhpev_drv.f90"
      INTEGER QI
      INTEGER IL(N+1)
      INTEGER OW(N+1)
      COMPLEX(DP) CTMP
      COMPLEX(DP) WORK(N+1)
!# 587 "zhpev_drv.f90"
!     ..
!     .. Local __SCALARs ..
      INTEGER :: I, IINFO, J, K, JL, KL, J1, I1, I2, NI1, L, IERR
      INTEGER :: ibeg, iend, nr
      INTEGER, EXTERNAL :: ldim_cyclic, lind_cyclic
!     ..
!# 594 "zhpev_drv.f90"
!     .. Executable Statements ..
!
!     Test the input arguments
!
!     Quick return if possible
!
      IF( N == 0 ) THEN
        RETURN
      END IF
!# 604 "zhpev_drv.f90"
      nr = ldim_cyclic( n, nproc, me )
      !
      IF( nr /= nrl ) &
         CALL lax_error__( " pzupgtr ", " inconsistent dimensions ", nrl )
      !
      ibeg = lind_cyclic( 1,  n, nproc, me )
      iend = lind_cyclic( nr, n, nproc, me )
!
      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC)
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO
!
!        Unpack the vectors which define the elementary reflectors and
!        set the first row and column of Q equal to those of the unit
!        matrix
!
      IF(OW(1).EQ.ME) THEN
        Q( IL(1), 1 ) = ONE
        DO KL = 2, NRL
          Q( KL, 1 ) = ZERO
        END DO
        DO J = 2, N
          Q( IL(1), J ) = ZERO
        END DO
      ELSE
        DO KL = 1, NRL
          Q( KL, 1 ) = ZERO
        END DO
      ENDIF
!# 640 "zhpev_drv.f90"
      DO J = 2, N
        IF(OW(J+1).EQ.ME) THEN
          J1 = IL(J+1)
        ELSE
          J1 = IL(J+1) + 1
        ENDIF
        DO KL = J1, NRL
          Q( KL, J ) = AP( KL, J-1 )
        END DO
      END DO
!# 651 "zhpev_drv.f90"
      IF( N.GT.1 ) THEN
!
!           Generate Q(2:n,2:n)
!
        DO I = N-1, 1, -1
!
!         Apply H(i) to A(i:m,i:n) from the left
!
          IF( I.LT.(N-1) ) THEN
!# 661 "zhpev_drv.f90"
            IF(OW(I+1).EQ.ME) THEN
              Q( IL(I+1), I+1 ) = ONE
            END IF
!
!           Form  H * C
!
            IF( TAU(I).NE.ZERO ) THEN
!
!             w := C' * v
!
              IF(OW(I+1).EQ.ME) THEN
                I1 = IL(I+1)
              ELSE
                I1 = IL(I+1) + 1
              ENDIF
              !
              IF( N-1-I > OPT_ZGEMV ) THEN
                 IF( NRL-I1+1 > 0 ) THEN
                    CALL zgemv( 'C', NRL-I1+1, N-1-I, one, Q( I1, I+1+1 ), ldq, Q( I1, I+1 ), 1, zero, work, 1 )
                 ELSE
                    work( 1 : N-1-I ) = 0.0_DP
                 END IF
              ELSE
                 DO J = 1, N-1-I
                    CTMP = ZERO
                    DO KL = I1, NRL
                       CTMP = CTMP + CONJG( Q( KL, J+I+1 ) ) * Q( KL,I+1 )
                     END DO
                     WORK(J) = CTMP
                  END DO
              END IF
!# 697 "zhpev_drv.f90"
              !
              !  C := C - v * w'
              !
              IF( N-1-I > opt_zgerc ) THEN
                 IF( NRL-I1+1 > 0 ) THEN
                    CALL zgerc( NRL-I1+1, N-1-I, -TAU(I), Q(I1, I+1), 1, work, 1, Q( I1, 1+I+1 ), ldq )
                 END IF
              ELSE
                 DO J = 1, N-1-I
                   CTMP = -TAU(I) * CONJG( WORK( J ) )
                   DO KL = I1, NRL
                     Q( KL, J+I+1 ) = Q( KL, J+I+1 ) + CTMP * Q(KL, I+1)
                   END DO
                 END DO
              END IF
            END IF
          END IF
!# 715 "zhpev_drv.f90"
          IF( I.LT.(N-1) ) THEN
            IF(OW(I+2).EQ.ME) THEN
              I2 = IL(I+2)              ! I+2
            ELSE
              I2 = IL(I+2) + 1          ! local ind. of the first element > I+2
            ENDIF
            NI1 = NRL - I2 + 1          ! N-I-1
            IF ( NI1 > 0 ) CALL zscal( NI1, -TAU( I ), Q( I2, I+1 ), 1 )
          END IF
!# 725 "zhpev_drv.f90"
          IF(OW(I+1).EQ.ME) THEN
            Q( IL(I+1), I+1 ) = ONE - TAU( I )
          END IF
!
!             Set A(1:i-1,i) to zero
!
          DO L = 1, I - 1
            IF(OW(L+1).EQ.ME) THEN
              Q( IL(L+1), I+1 ) = ZERO
            END IF
          END DO
        END DO
      END IF
!# 740 "zhpev_drv.f90"
      RETURN
!# 742 "zhpev_drv.f90"
!
!     End of ZUPGTR
!
      END SUBROUTINE pzupgtr
!# 747 "zhpev_drv.f90"
!==----------------------------------------------==!
!# 749 "zhpev_drv.f90"
      SUBROUTINE pzsteqr( compz, n, nrl, d, e, z, ldz, nproc, me, comm )
!
!  Parallel MPI version of the LAPACK routine ZHPTRD
!
!     Carlo Cavazzoni (carlo.cavazzoni@cineca.it) -- CINECA
!     Dicember 12, 1999
!
!  REFERENCES :
!
!     NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING.
!     W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, AND W.T. VETTERLING,
!     CAMBRIDGE UNIVERSITY PRESS, CAMBRIDGE.
!
!     PARALLEL NUMERICAL ALGORITHMS,
!     T.L. FREEMAN AND C.PHILLIPS,
!     PRENTICE HALL INTERNATIONAL (1992).
!
!     LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!
      USE laxlib_parallel_include
!# 772 "zhpev_drv.f90"
      IMPLICIT NONE
!# 774 "zhpev_drv.f90"
      include 'laxlib_kinds.fh'
!# 776 "zhpev_drv.f90"
!     .. __SCALAR Arguments ..
      CHARACTER          COMPZ
      INTEGER            LDZ, N, NRL, NPROC, ME, comm
!     ..
!     .. Array Arguments ..
      REAL(DP)   D( * ), E( * )
      COMPLEX(DP)         Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  PZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band complex Hermitian matrix can also
!  be found if PZHPTRD has been used to reduce this
!  matrix to tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  Hermitian matrix.  On entry, Z must contain the
!                  unitary matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the mglobal atrix AP.  N >= 0.
!
!  NRL     (input) INTEGER
!          The number of local rows of the matrix AP. NRL >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) COMPLEX(DP) array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the unitary
!          matrix used in the reduction to tridiagonal form.
!          On exit if COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original Hermitian matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!          The rows of the matrix are distributed among processors
!          with blocking factor 1, i.e. for NPROC = 4 :
!              ROW Index | Processor index owning the row
!                    1   |    0
!                    2   |    1
!                    3   |    2
!                    4   |    3
!                    5   |    0
!                    6   |    1
!                    ..  |    ..
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,NRL).
!
!  NPROC   (input) INTEGER
!          Number of processors
!
!  ME      (input) INTEGER
!          Index of the local processor  ( 0, 1, 2, ..., NPROC-1 )
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(DP)  RZERO, RONE, TWO, THREE, CTEMP, STEMP
      PARAMETER          ( RZERO = 0.0_DP, RONE = 1.0_DP, TWO = 2.0_DP, &
     &                   THREE = 3.0_DP )
      COMPLEX(DP)         ZERO, ONE,ZTEMP
      PARAMETER          ( ZERO = ( 0.0_DP, 0.0_DP ), ONE = ( 1.0_DP, 0.0_DP ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!# 863 "zhpev_drv.f90"
      INTEGER  :: QI, KL, INFO
      INTEGER  :: IL(N+1)
      INTEGER  :: OW(N+1)
      REAL(DP) :: WORK(2*N)
      REAL(DP) :: dvar(6)
!# 869 "zhpev_drv.f90"
!     .. Local __SCALARs ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
     &                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,   &
     &                   NM1, NMAXIT, IERR
      REAL(DP)   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,   &
     &                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL(DP)   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA
      EXTERNAL           ZLASET, ZLASR, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!# 894 "zhpev_drv.f90"
! DEBUG START
!      if( n > 400 ) then
!      write( 4000 + me, * ) LDZ, N, NRL, NPROC, ME, comm
!      do i = 1, n
!         write( 4000 + me, * ) d( i )
!      end do
!      do i = 1, n
!         write( 4000 + me, * ) e( i )
!      end do
!      do j = 1, n
!      do i = 1, nrl
!         write( 4000 + me, * ) z( i, j )
!      end do
!      end do
!      close( 4000 + me )
!      call mpi_barrier( comm, i )
!      stop 'qui'
!      end if
! DEBUG END
!# 914 "zhpev_drv.f90"
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( (LDZ.LT.1) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX(1,NRL) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF(N.LE.0) THEN
        RETURN
      END IF
!
      DO I = 1,N+1
        QI     = (I-1)/NPROC
        OW(I)  = MOD((I-1),NPROC)
        IF(ME .le. OW(I) ) then
          IL(I) = QI + 1
        ELSE
          IL(I) = QI
        END IF
      END DO
!# 952 "zhpev_drv.f90"
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 .AND. OW(1).EQ.ME ) Z( IL(1), 1 ) = ONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!     We ensure that all procs have the same data!
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = RONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
      !
      dvar(1) = EPS
      dvar(2) = EPS2
      dvar(3) = SAFMIN
      dvar(4) = SAFMAX
      dvar(5) = SSFMAX
      dvar(6) = SSFMIN
      !
!# 978 "zhpev_drv.f90"
      !
      EPS     = dvar(1)
      EPS2    = dvar(2)
      SAFMIN  = dvar(3)
      SAFMAX  = dvar(4)
      SSFMAX  = dvar(5)
      SSFMIN  = dvar(6)
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ.EQ.2 ) THEN
        CALL ZLASET( 'Full', NRL, N, ZERO, ZERO, Z, LDZ )
        DO J = 1, N
          IF(OW(J).EQ.ME) THEN
            Z( IL(J), J ) = ONE
          END IF
        END DO
      END IF
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
!# 1010 "zhpev_drv.f90"
      IF( L1 .GT. N )   GO TO 160
!# 1012 "zhpev_drv.f90"
      IF( L1 .GT. 1 )   E( L1-1 ) = RZERO
!# 1014 "zhpev_drv.f90"
      IF( me == 0 ) THEN
!# 1016 "zhpev_drv.f90"
         IF( L1.LE.NM1 ) THEN
            DO M = L1, NM1
               TST = DABS( E( M ) )
               IF( TST .EQ. RZERO )        GO TO 30
               IF( TST .LE. ( SQRT(DABS(D(M)))*SQRT(DABS(D(M+1))) ) * EPS ) THEN
                  E( M ) = RZERO
                  GO TO 30
               END IF
            END DO
         END IF
         M = N
!
   30    CONTINUE
!# 1030 "zhpev_drv.f90"
      END IF
!# 1039 "zhpev_drv.f90"
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )  GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.RZERO )   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( DABS( D( LEND ) ).LT.DABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GT.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
!# 1076 "zhpev_drv.f90"
         IF( me == 0 ) THEN
!# 1078 "zhpev_drv.f90"
            IF( L.NE.LEND ) THEN
               LENDM1 = LEND - 1
               DO M = L, LENDM1
                  TST = DABS( E( M ) )**2
                  IF( TST.LE.( EPS2*DABS(D(M)) )*DABS(D(M+1))+ SAFMIN )GO TO 60
               END DO
            END IF
!
            M = LEND
!
   60       CONTINUE
!# 1090 "zhpev_drv.f90"
         END IF
!# 1097 "zhpev_drv.f90"
         IF( M.LT.LEND )  E( M ) = RZERO
         P = D( L )
         IF( M.EQ.L ) THEN
!
!        Eigenvalue found.
!
            D( L ) = P
            L = L + 1
            IF( L.LE.LEND )   GO TO 40
            GO TO 140
         END IF
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CTEMP = WORK( L )
               STEMP = WORK( N-1+L )
               IF( ( CTEMP.NE.RONE ) .OR. ( STEMP.NE.RZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, 1+L )
                     Z( KL, 1+L ) = CTEMP*ZTEMP - STEMP*Z( KL, L )
                     Z( KL, L )   = STEMP*ZTEMP + CTEMP*Z( KL, L )
                  END DO
               END IF
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = RZERO
            L = L + 2
            IF( L.LE.LEND )     GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT )   GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         !
         ! iteration is performed on one processor and results broadcast
         ! to all others to prevent potential problems if all processors
         ! do not behave in exactly the same way (even with the same data!)
         !
         if ( me == 0 ) then
!# 1149 "zhpev_drv.f90"
            G = ( D( L+1 )-P ) / ( TWO*E( L ) )
            R = DLAPY2( G, RONE )
            G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
            S = RONE
            C = RONE
            P = RZERO
!
!        Inner loop
!
            MM1 = M - 1
            DO I = MM1, L, -1
               F = S*E( I )
               B = C*E( I )
               CALL DLARTG( G, F, C, S, R )
               IF( I.NE.M-1 )  E( I+1 ) = R
               G = D( I+1 ) - P
               R = ( D( I )-G )*S + TWO*C*B
               P = S*R
               D( I+1 ) = G + P
               G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
               IF( ICOMPZ.GT.0 ) THEN
                  WORK( I ) = C
                  WORK( N-1+I ) = -S
               END IF
            END DO
            D( L ) = D( L ) - P
            E( L ) = G
         END IF
!# 1187 "zhpev_drv.f90"
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
!# 1195 "zhpev_drv.f90"
           DO J = M - L + 1 - 1, 1, -1
             CTEMP =  WORK( L + J -1)
             STEMP =  WORK( N-1+L + J-1)
             IF( ( CTEMP.NE.RONE ) .OR. ( STEMP.NE.RZERO ) ) THEN
               DO KL = 1, NRL
                 ZTEMP = Z( KL, J+1+L-1 )
                 Z( KL, J+1+L-1 ) = CTEMP*ZTEMP - STEMP*Z( KL, J+L-1 )
                 Z( KL, J+L-1 ) = STEMP*ZTEMP + CTEMP*Z( KL, J+L-1 )
               END DO
             END IF
           END DO
         END IF
!
         GO TO 40
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
!# 1218 "zhpev_drv.f90"
         IF( me == 0 ) THEN
!# 1220 "zhpev_drv.f90"
            IF( L.NE.LEND ) THEN
               LENDP1 = LEND + 1
               DO 100 M = L, LENDP1, -1
                  TST = DABS( E( M-1 ) )**2
                  IF( TST.LE.(EPS2*DABS(D(M)))*DABS(D(M-1))+ SAFMIN )GO TO 110
  100          CONTINUE
            END IF
!
            M = LEND
!
  110       CONTINUE
!# 1232 "zhpev_drv.f90"
         END IF
!# 1239 "zhpev_drv.f90"
         IF( M.GT.LEND )   E( M-1 ) = RZERO
         P = D( L )
         IF( M.EQ.L ) THEN
!
!        Eigenvalue found.
!
            D( L ) = P
            L = L - 1
            IF( L.GE.LEND )   GO TO 90
            GO TO 140
         END IF
!
!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CTEMP = WORK( M )
               STEMP = WORK( N-1+M )
               IF( ( CTEMP.NE.RONE ) .OR. ( STEMP.NE.RZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, L)
                     Z( KL, L )   = CTEMP*ZTEMP - STEMP*Z( KL, L-1 )
                     Z( KL, L-1 ) = STEMP*ZTEMP + CTEMP*Z( KL, L-1 )
                  END DO
               END IF
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = RZERO
            L = L - 2
            IF( L.GE.LEND )    GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT.EQ.NMAXIT )  GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         !
         ! iteration is performed on one processor and results broadcast
         ! to all others to prevent potential problems if all processors
         ! do not behave in exactly the same way (even with the same data!)
         !
         if ( me == 0 ) then
!# 1291 "zhpev_drv.f90"
            G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
            R = DLAPY2( G, RONE )
            G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
            S = RONE
            C = RONE
            P = RZERO
!
!        Inner loop
!
            LM1 = L - 1
            DO I = M, LM1
               F = S*E( I )
               B = C*E( I )
               CALL DLARTG( G, F, C, S, R )
               IF( I.NE.M )     E( I-1 ) = R
               G = D( I ) - P
               R = ( D( I+1 )-G )*S + TWO*C*B
               P = S*R
               D( I ) = G + P
               G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
               IF( ICOMPZ.GT.0 ) THEN
                  WORK( I ) = C
                  WORK( N-1+I ) = S
               END IF
            END DO
            D( L ) = D( L ) - P
            E( LM1 ) = G
         END IF
!# 1329 "zhpev_drv.f90"
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ.GT.0 ) THEN
!# 1337 "zhpev_drv.f90"
            DO J = 1, L - M
               CTEMP = WORK( M+J-1 )
               STEMP = WORK( N-1+M+J-1 )
               IF( ( CTEMP.NE.RONE ) .OR. ( STEMP.NE.RZERO ) ) THEN
                  DO KL = 1, NRL
                     ZTEMP = Z( KL, J+M )
                     Z( KL, J+M )   = CTEMP*ZTEMP - STEMP*Z(KL, J+M-1)
                     Z( KL, J+M-1 ) = STEMP*ZTEMP + CTEMP*Z(KL, J+M-1)
                  END DO
               END IF
            END DO
         END IF
!
         GO TO 90
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
!# 1358 "zhpev_drv.f90"
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
     &                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
     &                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
     &                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
     &                N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT .EQ. NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ) .NE. RZERO )  INFO = INFO + 1
  150    CONTINUE
         WRITE(6,*) 'WARNING pzsteqr, convergence not achieved INFO = ', INFO
         RETURN
      END IF
      GO TO 10
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
!# 1386 "zhpev_drv.f90"
      IF( ICOMPZ.EQ.0 ) THEN
!
!        Use Quick Sort
!
         CALL DLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( NRL, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
!# 1414 "zhpev_drv.f90"
      RETURN
!
!     End of ZSTEQR
!
      END SUBROUTINE pzsteqr
!# 1420 "zhpev_drv.f90"
!==----------------------------------------------==!
!# 1577 "zhpev_drv.f90"
END MODULE zhpev_module
!# 1580 "zhpev_drv.f90"
!==----------------------------------------------==!
!# 1583 "zhpev_drv.f90"
   SUBROUTINE zhpev_drv_x( JOBZ, UPLO, N, AP, W, Z, LDZ )
!# 1585 "zhpev_drv.f90"
        use zhpev_module
        IMPLICIT NONE
        include 'laxlib_kinds.fh'
!# 1589 "zhpev_drv.f90"
        CHARACTER ::       JOBZ, UPLO
        INTEGER   ::       IOPT, INFO, LDZ, N
        COMPLEX(DP) ::  AP( * ), Z( LDZ, * )
        REAL(DP) ::  W( * )
        REAL(DP), ALLOCATABLE :: RWORK(:)
        COMPLEX(DP), ALLOCATABLE :: ZWORK(:)
!# 1596 "zhpev_drv.f90"
        ALLOCATE( rwork( MAX(1, 3*n-2) ), zwork( MAX(1, 2*n-1)) )
        CALL ZHPEV(jobz, uplo, n, ap, w, z, ldz, zwork, rwork, INFO)
        DEALLOCATE( rwork, zwork )
        IF( INFO .NE. 0 ) THEN
          CALL lax_error__( ' zhpev_drv ', ' diagonalization failed ',INFO )
        END IF
!# 1603 "zhpev_drv.f90"
        RETURN
   END SUBROUTINE
!# 1606 "zhpev_drv.f90"
!==----------------------------------------------==!
!# 1608 "zhpev_drv.f90"
   SUBROUTINE pzhpev_drv_x( jobz, ap, lda, w, z, ldz, nrl, n, nproc, mpime, comm )
!# 1610 "zhpev_drv.f90"
     use zhpev_module
!# 1612 "zhpev_drv.f90"
     IMPLICIT NONE
     include 'laxlib_kinds.fh'
     CHARACTER :: JOBZ
     INTEGER, INTENT(IN) :: lda, ldz, nrl, n, nproc, mpime
     INTEGER, INTENT(IN) :: comm
     COMPLEX(DP) :: ap( lda, * ), z( ldz, * )
     REAL(DP) :: w( * )
     REAL(DP), ALLOCATABLE :: rwork( : )
     COMPLEX(DP), ALLOCATABLE :: cwork( : )
     !
     ALLOCATE( rwork( n ) )
     ALLOCATE( cwork( n ) )
     !
     CALL pzhptrd( n, nrl, ap, lda, w, rwork, cwork, nproc, mpime, comm)
!# 1627 "zhpev_drv.f90"
     IF( jobz == 'V' .OR. jobz == 'v' ) THEN
        CALL pzupgtr( n, nrl, ap, lda, cwork, z, ldz, nproc, mpime, comm)
     END IF
!# 1631 "zhpev_drv.f90"
     CALL pzsteqr( jobz, n, nrl, w, rwork, z, ldz, nproc, mpime, comm)
!# 1633 "zhpev_drv.f90"
     DEALLOCATE( cwork )
     DEALLOCATE( rwork )
!# 1636 "zhpev_drv.f90"
     RETURN
   END SUBROUTINE
