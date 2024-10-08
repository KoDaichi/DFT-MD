      SUBROUTINE DLAED4X( N, I, D, Z, DELTA, RHO, DLAM, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
*     Courant Institute, NAG Ltd., and Rice University
*     December 23, 1999
*
*     .. Scalar Arguments ..
      INTEGER            I, INFO, N
      DOUBLE PRECISION   DLAM, RHO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), DELTA( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  This subroutine computes the I-th updated eigenvalue of a symmetric
*  rank-one modification to a diagonal matrix whose elements are
*  given in the array d, and that
*
*             D(i) < D(j)  for  i < j
*
*  and that RHO > 0.  This is arranged by the calling routine, and is
*  no loss in generality.  The rank-one modified system is thus
*
*             diag( D )  +  RHO *  Z * Z_transpose.
*
*  where we assume the Euclidean norm of Z is 1.
*
*  The method consists of approximating the rational functions in the
*  secular equation by simpler interpolating rational functions.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The length of all arrays.
*
*  I      (input) INTEGER
*         The index of the eigenvalue to be computed.  1 <= I <= N.
*
*  D      (input) DOUBLE PRECISION array, dimension (N)
*         The original eigenvalues.  It is assumed that they are in
*         order, D(I) < D(J)  for I < J.
*
*  Z      (input) DOUBLE PRECISION array, dimension (N)
*         The components of the updating vector.
*
*  DELTA  (output) DOUBLE PRECISION array, dimension (N)
*         If N .ne. 1, DELTA contains (D(j) - lambda_I) in its  j-th
*         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
*         contains the information necessary to construct the
*         eigenvectors.
*
*  RHO    (input) DOUBLE PRECISION
*         The scalar in the symmetric updating formula.
*
*  DLAM   (output) DOUBLE PRECISION
*         The computed lambda_I, the I-th updated eigenvalue.
*
*  INFO   (output) INTEGER
*         = 0:  successful exit
*         > 0:  if INFO = 1, the updating process failed.
*
*  Internal Parameters
*  ===================
*
*  Logical variable ORGATI (origin-at-i?) is used for distinguishing
*  whether D(i) or D(i+1) is treated as the origin.
*
*            ORGATI = .true.    origin at i
*            ORGATI = .false.   origin at i+1
*
*   Logical variable SWTCH3 (switch-for-3-poles?) is for noting
*   if we are working with THREE poles!
*
*   MAXIT is the maximum number of iterations allowed for each
*   eigenvalue.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ren-Cang Li, Computer Science Division, University of California
*     at Berkeley, USA
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE, FOUR, EIGHT, TEN
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0,
     $                   TEN = 10.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ORGATI, SWTCH, SWTCH3
      INTEGER            II, IIM1, IIP1, IP1, ITER, J, NITER
      DOUBLE PRECISION   A, B, C, DEL, DLTLB, DLTUB, DPHI, DPSI, DW,
     $                   EPS, ERRETM, ETA, MIDPT, PHI, PREW, PSI,
     $                   RHOINV, TAU, TEMP, TEMP1, W, WK
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   ZZ( 3 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAED5, DLAED6X
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
!$      REAL(8) ZZZ1(16)
!$      REAL(8) ZZZ2(16)
!$      REAL(8) ZZZ3(16)
!$      INTEGER, EXTERNAL :: OMP_GET_NUM_THREADS
!$      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
*     .. Executable Statements ..
*
*     Since this routine is called in an inner loop, we do no argument
*     checking.
*
*     Quick return for N=1 and 2.
*
      INFO = 0
      IF( N.EQ.1 ) THEN
*
*         Presumably, I=1 upon entry
*
         DLAM = D( 1 ) + RHO*Z( 1 )*Z( 1 )
         DELTA( 1 ) = ONE
         RETURN
      END IF
      IF( N.EQ.2 ) THEN
         CALL DLAED5( I, D, Z, DELTA, RHO, DLAM )
         RETURN
      END IF
*
*     Compute machine epsilon
*
      EPS = DLAMCH( 'Epsilon' )
      RHOINV = ONE / RHO
*
*     The case I = N
*
      IF( I.EQ.N ) THEN
*
*        Initialize some basic variables
*
         II = N - 1
         NITER = 1
*
*        Calculate initial guess
*
         MIDPT = RHO / TWO
*
*        If ||Z||_2 is not one, then TEMP should be set to
*        RHO * ||Z||_2^2 / TWO
*
!$OMP PARALLEL DO
         DO 10 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) )
   10    CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         DO 11 J = 1, N
            DELTA( J ) = DELTA( J ) - MIDPT
   11    CONTINUE
!$OMP END PARALLEL DO
*
         PSI = ZERO
!$    ZZZ1(1)=PSI
!$OMP PARALLEL PRIVATE(PSI)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$OMP DO
         DO 20 J = 1, N - 2
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
   20    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
*
         C = RHOINV + PSI
         W = C + Z( II )*Z( II ) / DELTA( II ) +
     $       Z( N )*Z( N ) / DELTA( N )
*
         IF( W.LE.ZERO ) THEN
            TEMP = Z( N-1 )*Z( N-1 ) / ( D( N )-D( N-1 )+RHO ) +
     $             Z( N )*Z( N ) / RHO
            IF( C.LE.TEMP ) THEN
               TAU = RHO
            ELSE
               DEL = D( N ) - D( N-1 )
               A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
               B = Z( N )*Z( N )*DEL
               IF( A.LT.ZERO ) THEN
                  TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
               ELSE
                  TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
               END IF
            END IF
*
*           It can be proved that
*               D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
*
            DLTLB = MIDPT
            DLTUB = RHO
         ELSE
            DEL = D( N ) - D( N-1 )
            A = -C*DEL + Z( N-1 )*Z( N-1 ) + Z( N )*Z( N )
            B = Z( N )*Z( N )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( SQRT( A*A+FOUR*B*C )-A )
            ELSE
               TAU = ( A+SQRT( A*A+FOUR*B*C ) ) / ( TWO*C )
            END IF
*
*           It can be proved that
*               D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2
*
            DLTLB = ZERO
            DLTUB = MIDPT
         END IF
*
!$OMP PARALLEL DO
         DO 30 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) )
   30    CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         DO 31 J = 1, N
            DELTA( J ) = DELTA( J ) - TAU
   31    CONTINUE
!$OMP END PARALLEL DO
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 40 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PSI = PSI + WK
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + WK*(II+1-J)
   40    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            DLAM = D( I ) + TAU
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
         A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $       DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
         B = DELTA( N-1 )*DELTA( N )*W
         IF( C.LT.ZERO )
     $      C = ABS( C )
         IF( C.EQ.ZERO ) THEN
*          ETA = B/A
*           ETA = RHO - TAU
            ETA = DLTUB - TAU
         ELSE IF( A.GE.ZERO ) THEN
            ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
         ELSE
            ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GT.ZERO )
     $      ETA = -W / ( DPSI+DPHI )
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
!$OMP PARALLEL DO
         DO 50 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
   50    CONTINUE
!$OMP END PARALLEL DO
*
         TAU = TAU + ETA
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 60 J = 1, II
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PSI = PSI + WK
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + WK*(II+1-J)
   60    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         TEMP = Z( N ) / DELTA( N )
         PHI = Z( N )*TEMP
         DPHI = TEMP*TEMP
         ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $            ABS( TAU )*( DPSI+DPHI )
*
         W = RHOINV + PHI + PSI
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 90 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               DLAM = D( I ) + TAU
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            C = W - DELTA( N-1 )*DPSI - DELTA( N )*DPHI
            A = ( DELTA( N-1 )+DELTA( N ) )*W -
     $          DELTA( N-1 )*DELTA( N )*( DPSI+DPHI )
            B = DELTA( N-1 )*DELTA( N )*W
            IF( A.GE.ZERO ) THEN
               ETA = ( A+SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A-SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GT.ZERO )
     $         ETA = -W / ( DPSI+DPHI )
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
!$OMP PARALLEL DO
            DO 70 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
   70       CONTINUE
!$OMP END PARALLEL DO
*
            TAU = TAU + ETA
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
            DO 80 J = 1, II
               TEMP = Z( J ) / DELTA( J )
               WK = Z( J )*TEMP
               PSI = PSI + WK
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + WK*(II+1-J)
   80       CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            TEMP = Z( N ) / DELTA( N )
            PHI = Z( N )*TEMP
            DPHI = TEMP*TEMP
            ERRETM = EIGHT*( -PHI-PSI ) + ERRETM - PHI + RHOINV +
     $               ABS( TAU )*( DPSI+DPHI )
*
            W = RHOINV + PHI + PSI
   90    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         DLAM = D( I ) + TAU
         GO TO 250
*
*        End for the case I = N
*
      ELSE
*
*        The case for I < N
*
         NITER = 1
         IP1 = I + 1
*
*        Calculate initial guess
*
         DEL = D( IP1 ) - D( I )
         MIDPT = DEL / TWO
!$OMP PARALLEL DO
         DO 100 J = 1, N
            DELTA( J ) = ( D( J )-D( I ) )
  100    CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
         DO 101 J = 1, N
            DELTA( J ) = DELTA( J ) - MIDPT
  101    CONTINUE
!$OMP END PARALLEL DO
*
         PSI = ZERO
!$    ZZZ1(1)=PSI
!$OMP PARALLEL PRIVATE(PSI)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$OMP DO
         DO 110 J = 1, I - 1
            PSI = PSI + Z( J )*Z( J ) / DELTA( J )
  110    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
*
         PHI = ZERO
!$    ZZZ1(1)=PHI
!$OMP PARALLEL PRIVATE(PHI)
!$    PHI=0.0D+00
!$OMP MASTER
!$    PHI=ZZZ1(1)
!$OMP END MASTER
!$OMP DO
         DO 120 J = N, I + 2, -1
            PHI = PHI + Z( J )*Z( J ) / DELTA( J )
  120    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PHI
!$OMP BARRIER
!$OMP MASTER
!$    PHI=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PHI=PHI+ZZZ1(J__)
!$    ENDDO
!$     ZZZ1(1)=PHI
!$OMP END MASTER
!$OMP END PARALLEL
!$     PHI=ZZZ1(1)
         C = RHOINV + PSI + PHI
         W = C + Z( I )*Z( I ) / DELTA( I ) +
     $       Z( IP1 )*Z( IP1 ) / DELTA( IP1 )
*
         IF( W.GT.ZERO ) THEN
*
*           d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
*
*           We choose d(i) as origin.
*
            ORGATI = .TRUE.
            A = C*DEL + Z( I )*Z( I ) + Z( IP1 )*Z( IP1 )
            B = Z( I )*Z( I )*DEL
            IF( A.GT.ZERO ) THEN
               TAU = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            ELSE
               TAU = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = ZERO
            DLTUB = MIDPT
         ELSE
*
*           (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
*
*           We choose d(i+1) as origin.
*
            ORGATI = .FALSE.
            A = C*DEL - Z( I )*Z( I ) - Z( IP1 )*Z( IP1 )
            B = Z( IP1 )*Z( IP1 )*DEL
            IF( A.LT.ZERO ) THEN
               TAU = TWO*B / ( A-SQRT( ABS( A*A+FOUR*B*C ) ) )
            ELSE
               TAU = -( A+SQRT( ABS( A*A+FOUR*B*C ) ) ) / ( TWO*C )
            END IF
            DLTLB = -MIDPT
            DLTUB = ZERO
         END IF
*
         IF( ORGATI ) THEN
!$OMP PARALLEL DO
            DO 130 J = 1, N
               DELTA( J ) = ( D( J )-D( I ) )
  130       CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
            DO 131 J = 1, N
               DELTA( J ) = DELTA( J ) - TAU
  131       CONTINUE
!$OMP END PARALLEL DO
         ELSE
!$OMP PARALLEL DO
            DO 140 J = 1, N
               DELTA( J ) = ( D( J )-D( IP1 ) )
  140       CONTINUE
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
            DO 141 J = 1, N
               DELTA( J ) = DELTA( J ) - TAU
  141       CONTINUE
!$OMP END PARALLEL DO
         END IF
         IF( ORGATI ) THEN
            II = I
         ELSE
            II = I + 1
         END IF
         IIM1 = II - 1
         IIP1 = II + 1
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 150 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PSI = PSI + WK
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + WK*(IIM1+1-J)
  150    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
!$    ZZZ1(1)=PHI
!$    ZZZ2(1)=DPHI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PHI,DPHI,ERRETM)
!$    PHI=0.0D+00
!$OMP MASTER
!$    PHI=ZZZ1(1)
!$OMP END MASTER
!$    DPHI=0.0D+00
!$OMP MASTER
!$    DPHI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 160 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PHI = PHI + WK
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + WK*(J+1-IIP1)
  160    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PHI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPHI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PHI=0.0D+00
!$    DPHI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PHI=PHI+ZZZ1(J__)
!$     DPHI=DPHI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PHI
!$     ZZZ2(1)=DPHI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PHI=ZZZ1(1)
!$     DPHI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
*
         W = RHOINV + PHI + PSI
*
*        W is the value of the secular function with
*        its ii-th element removed.
*
         SWTCH3 = .FALSE.
         IF( ORGATI ) THEN
            IF( W.LT.ZERO )
     $         SWTCH3 = .TRUE.
         ELSE
            IF( W.GT.ZERO )
     $         SWTCH3 = .TRUE.
         END IF
         IF( II.EQ.1 .OR. II.EQ.N )
     $      SWTCH3 = .FALSE.
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = W + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU )*DW
*
*        Test for convergence
*
         IF( ABS( W ).LE.EPS*ERRETM ) THEN
            IF( ORGATI ) THEN
               DLAM = D( I ) + TAU
            ELSE
               DLAM = D( IP1 ) + TAU
            END IF
            GO TO 250
         END IF
*
         IF( W.LE.ZERO ) THEN
            DLTLB = MAX( DLTLB, TAU )
         ELSE
            DLTUB = MIN( DLTUB, TAU )
         END IF
*
*        Calculate the new step
*
         NITER = NITER + 1
         IF( .NOT.SWTCH3 ) THEN
            IF( ORGATI ) THEN
               C = W - DELTA( IP1 )*DW - ( D( I )-D( IP1 ) )*
     $             ( Z( I ) / DELTA( I ) )**2
            ELSE
               C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $             ( Z( IP1 ) / DELTA( IP1 ) )**2
            END IF
            A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $          DELTA( I )*DELTA( IP1 )*DW
            B = DELTA( I )*DELTA( IP1 )*W
            IF( C.EQ.ZERO ) THEN
               IF( A.EQ.ZERO ) THEN
                  IF( ORGATI ) THEN
                     A = Z( I )*Z( I ) + DELTA( IP1 )*DELTA( IP1 )*
     $                   ( DPSI+DPHI )
                  ELSE
                     A = Z( IP1 )*Z( IP1 ) + DELTA( I )*DELTA( I )*
     $                   ( DPSI+DPHI )
                  END IF
               END IF
               ETA = B / A
            ELSE IF( A.LE.ZERO ) THEN
               ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
            ELSE
               ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
            END IF
         ELSE
*
*           Interpolation using THREE most relevant poles
*
            TEMP = RHOINV + PSI + PHI
            IF( ORGATI ) THEN
               TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $             ( D( IIM1 )-D( IIP1 ) )*TEMP1
               ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
               ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                   ( ( DPSI-TEMP1 )+DPHI )
            ELSE
               TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
               TEMP1 = TEMP1*TEMP1
               C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $             ( D( IIP1 )-D( IIM1 ) )*TEMP1
               ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                   ( DPSI+( DPHI-TEMP1 ) )
               ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
            END IF
            ZZ( 2 ) = Z( II )*Z( II )
            CALL DLAED6X( NITER, ORGATI, C, DELTA( IIM1 ), ZZ, W, ETA,
     $                   INFO )
            IF( INFO.NE.0 )
     $         GO TO 250
         END IF
*
*        Note, eta should be positive if w is negative, and
*        eta should be negative otherwise. However,
*        if for some reason caused by roundoff, eta*w > 0,
*        we simply use one Newton step instead. This way
*        will guarantee eta*w < 0.
*
         IF( W*ETA.GE.ZERO )
     $      ETA = -W / DW
         TEMP = TAU + ETA
         IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
            IF( W.LT.ZERO ) THEN
               ETA = ( DLTUB-TAU ) / TWO
            ELSE
               ETA = ( DLTLB-TAU ) / TWO
            END IF
         END IF
*
         PREW = W
*
  170    CONTINUE
!$OMP PARALLEL DO
         DO 180 J = 1, N
            DELTA( J ) = DELTA( J ) - ETA
  180    CONTINUE
!$OMP END PARALLEL DO
*
*        Evaluate PSI and the derivative DPSI
*
         DPSI = ZERO
         PSI = ZERO
         ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 190 J = 1, IIM1
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PSI = PSI + WK
            DPSI = DPSI + TEMP*TEMP
            ERRETM = ERRETM + WK*(IIM1+1-J)
  190    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
         ERRETM = ABS( ERRETM )
*
*        Evaluate PHI and the derivative DPHI
*
         DPHI = ZERO
         PHI = ZERO
!$    ZZZ1(1)=PHI
!$    ZZZ2(1)=DPHI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PHI,DPHI,ERRETM)
!$    PHI=0.0D+00
!$OMP MASTER
!$    PHI=ZZZ1(1)
!$OMP END MASTER
!$    DPHI=0.0D+00
!$OMP MASTER
!$    DPHI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
         DO 200 J = N, IIP1, -1
            TEMP = Z( J ) / DELTA( J )
            WK = Z( J )*TEMP
            PHI = PHI + WK
            DPHI = DPHI + TEMP*TEMP
            ERRETM = ERRETM + WK*(J+1-IIP1)
  200    CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PHI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPHI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PHI=0.0D+00
!$    DPHI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PHI=PHI+ZZZ1(J__)
!$     DPHI=DPHI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PHI
!$     ZZZ2(1)=DPHI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PHI=ZZZ1(1)
!$     DPHI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
*
         TEMP = Z( II ) / DELTA( II )
         DW = DPSI + DPHI + TEMP*TEMP
         TEMP = Z( II )*TEMP
         W = RHOINV + PHI + PSI + TEMP
         ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $            THREE*ABS( TEMP ) + ABS( TAU+ETA )*DW
*
         SWTCH = .FALSE.
         IF( ORGATI ) THEN
            IF( -W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         ELSE
            IF( W.GT.ABS( PREW ) / TEN )
     $         SWTCH = .TRUE.
         END IF
*
         TAU = TAU + ETA
*
*        Main loop to update the values of the array   DELTA
*
         ITER = NITER + 1
*
         DO 240 NITER = ITER, MAXIT
*
*           Test for convergence
*
            IF( ABS( W ).LE.EPS*ERRETM ) THEN
               IF( ORGATI ) THEN
                  DLAM = D( I ) + TAU
               ELSE
                  DLAM = D( IP1 ) + TAU
               END IF
               GO TO 250
            END IF
*
            IF( W.LE.ZERO ) THEN
               DLTLB = MAX( DLTLB, TAU )
            ELSE
               DLTUB = MIN( DLTUB, TAU )
            END IF
*
*           Calculate the new step
*
            IF( .NOT.SWTCH3 ) THEN
               IF( .NOT.SWTCH ) THEN
                  IF( ORGATI ) THEN
                     C = W - DELTA( IP1 )*DW -
     $                   ( D( I )-D( IP1 ) )*( Z( I ) / DELTA( I ) )**2
                  ELSE
                     C = W - DELTA( I )*DW - ( D( IP1 )-D( I ) )*
     $                   ( Z( IP1 ) / DELTA( IP1 ) )**2
                  END IF
               ELSE
                  TEMP = Z( II ) / DELTA( II )
                  IF( ORGATI ) THEN
                     DPSI = DPSI + TEMP*TEMP
                  ELSE
                     DPHI = DPHI + TEMP*TEMP
                  END IF
                  C = W - DELTA( I )*DPSI - DELTA( IP1 )*DPHI
               END IF
               A = ( DELTA( I )+DELTA( IP1 ) )*W -
     $             DELTA( I )*DELTA( IP1 )*DW
               B = DELTA( I )*DELTA( IP1 )*W
               IF( C.EQ.ZERO ) THEN
                  IF( A.EQ.ZERO ) THEN
                     IF( .NOT.SWTCH ) THEN
                        IF( ORGATI ) THEN
                           A = Z( I )*Z( I ) + DELTA( IP1 )*
     $                         DELTA( IP1 )*( DPSI+DPHI )
                        ELSE
                           A = Z( IP1 )*Z( IP1 ) +
     $                         DELTA( I )*DELTA( I )*( DPSI+DPHI )
                        END IF
                     ELSE
                        A = DELTA( I )*DELTA( I )*DPSI +
     $                      DELTA( IP1 )*DELTA( IP1 )*DPHI
                     END IF
                  END IF
                  ETA = B / A
               ELSE IF( A.LE.ZERO ) THEN
                  ETA = ( A-SQRT( ABS( A*A-FOUR*B*C ) ) ) / ( TWO*C )
               ELSE
                  ETA = TWO*B / ( A+SQRT( ABS( A*A-FOUR*B*C ) ) )
               END IF
            ELSE
*
*              Interpolation using THREE most relevant poles
*
               TEMP = RHOINV + PSI + PHI
               IF( SWTCH ) THEN
                  C = TEMP - DELTA( IIM1 )*DPSI - DELTA( IIP1 )*DPHI
                  ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*DPSI
                  ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*DPHI
               ELSE
                  IF( ORGATI ) THEN
                     TEMP1 = Z( IIM1 ) / DELTA( IIM1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIP1 )*( DPSI+DPHI ) -
     $                   ( D( IIM1 )-D( IIP1 ) )*TEMP1
                     ZZ( 1 ) = Z( IIM1 )*Z( IIM1 )
                     ZZ( 3 ) = DELTA( IIP1 )*DELTA( IIP1 )*
     $                         ( ( DPSI-TEMP1 )+DPHI )
                  ELSE
                     TEMP1 = Z( IIP1 ) / DELTA( IIP1 )
                     TEMP1 = TEMP1*TEMP1
                     C = TEMP - DELTA( IIM1 )*( DPSI+DPHI ) -
     $                   ( D( IIP1 )-D( IIM1 ) )*TEMP1
                     ZZ( 1 ) = DELTA( IIM1 )*DELTA( IIM1 )*
     $                         ( DPSI+( DPHI-TEMP1 ) )
                     ZZ( 3 ) = Z( IIP1 )*Z( IIP1 )
                  END IF
               END IF
               CALL DLAED6X( NITER, ORGATI, C, DELTA( IIM1 ), ZZ,W,ETA,
     $                      INFO )
               IF( INFO.NE.0 )
     $            GO TO 250
            END IF
*
*           Note, eta should be positive if w is negative, and
*           eta should be negative otherwise. However,
*           if for some reason caused by roundoff, eta*w > 0,
*           we simply use one Newton step instead. This way
*           will guarantee eta*w < 0.
*
            IF( W*ETA.GE.ZERO )
     $         ETA = -W / DW
            TEMP = TAU + ETA
            IF( TEMP.GT.DLTUB .OR. TEMP.LT.DLTLB ) THEN
               IF( W.LT.ZERO ) THEN
                  ETA = ( DLTUB-TAU ) / TWO
               ELSE
                  ETA = ( DLTLB-TAU ) / TWO
               END IF
            END IF
*
!$OMP PARALLEL DO
            DO 210 J = 1, N
               DELTA( J ) = DELTA( J ) - ETA
  210       CONTINUE
!$OMP END PARALLEL DO
*
            TAU = TAU + ETA
            PREW = W
*
*           Evaluate PSI and the derivative DPSI
*
            DPSI = ZERO
            PSI = ZERO
            ERRETM = ZERO
!$    ZZZ1(1)=PSI
!$    ZZZ2(1)=DPSI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PSI,DPSI,ERRETM)
!$    PSI=0.0D+00
!$OMP MASTER
!$    PSI=ZZZ1(1)
!$OMP END MASTER
!$    DPSI=0.0D+00
!$OMP MASTER
!$    DPSI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
            DO 220 J = 1, IIM1
               TEMP = Z( J ) / DELTA( J )
               WK = Z( J )*TEMP
               PSI = PSI + WK
               DPSI = DPSI + TEMP*TEMP
               ERRETM = ERRETM + WK*(J-IIM1+1)
  220       CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PSI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPSI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PSI=0.0D+00
!$    DPSI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PSI=PSI+ZZZ1(J__)
!$     DPSI=DPSI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PSI
!$     ZZZ2(1)=DPSI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PSI=ZZZ1(1)
!$     DPSI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
            ERRETM = ABS( ERRETM )
*
*           Evaluate PHI and the derivative DPHI
*
            DPHI = ZERO
            PHI = ZERO
!$    ZZZ1(1)=PHI
!$    ZZZ2(1)=DPHI
!$    ZZZ3(1)=ERRETM
!$OMP PARALLEL PRIVATE(TEMP,WK,PHI,DPHI,ERRETM)
!$    PHI=0.0D+00
!$OMP MASTER
!$    PHI=ZZZ1(1)
!$OMP END MASTER
!$    DPHI=0.0D+00
!$OMP MASTER
!$    DPHI=ZZZ2(1)
!$OMP END MASTER
!$    ERRETM=0.0D+00
!$OMP MASTER
!$    ERRETM=ZZZ3(1)
!$OMP END MASTER
!$OMP DO
            DO 230 J = N, IIP1, -1
               TEMP = Z( J ) / DELTA( J )
               WK = Z( J )*TEMP
               PHI = PHI + WK
               DPHI = DPHI + TEMP*TEMP
               ERRETM = ERRETM + WK*(J+1-IIP1)
  230       CONTINUE
!$OMP ENDDO
!$    ZZZ1(OMP_GET_THREAD_NUM()+1)=PHI
!$    ZZZ2(OMP_GET_THREAD_NUM()+1)=DPHI
!$    ZZZ3(OMP_GET_THREAD_NUM()+1)=ERRETM
!$OMP BARRIER
!$OMP MASTER
!$    PHI=0.0D+00
!$    DPHI=0.0D+00
!$    ERRETM=0.0D+00
!$    DO J__=1,OMP_GET_NUM_THREADS()
!$     PHI=PHI+ZZZ1(J__)
!$     DPHI=DPHI+ZZZ2(J__)
!$     ERRETM=ERRETM+ZZZ3(J__)
!$    ENDDO
!$     ZZZ1(1)=PHI
!$     ZZZ2(1)=DPHI
!$     ZZZ3(1)=ERRETM
!$OMP END MASTER
!$OMP END PARALLEL
!$     PHI=ZZZ1(1)
!$     DPHI=ZZZ2(1)
!$     ERRETM=ZZZ3(1)
*
            TEMP = Z( II ) / DELTA( II )
            DW = DPSI + DPHI + TEMP*TEMP
            TEMP = Z( II )*TEMP
            W = RHOINV + PHI + PSI + TEMP
            ERRETM = EIGHT*( PHI-PSI ) + ERRETM + TWO*RHOINV +
     $               THREE*ABS( TEMP ) + ABS( TAU )*DW
            IF( W*PREW.GT.ZERO .AND. ABS( W ).GT.ABS( PREW ) / TEN )
     $         SWTCH = .NOT.SWTCH
*
  240    CONTINUE
*
*        Return with INFO = 1, NITER = MAXIT and not converged
*
         INFO = 1
         IF( ORGATI ) THEN
            DLAM = D( I ) + TAU
         ELSE
            DLAM = D( IP1 ) + TAU
         END IF
*
      END IF
*
  250 CONTINUE
*
      RETURN
*
*     End of DLAED4X
*
      END
