      SUBROUTINE PDLAEDZ( N, N1, ID, Q, IQ, JQ, LDQ, DESCQ, Z, WORK )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     December 31, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ID, IQ, JQ, LDQ, N, N1
*     ..
*     .. Array Arguments ..
      INTEGER            DESCQ( * )
      DOUBLE PRECISION   Q( LDQ, * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLAEDZ Form the z-vector which consists of the last row of Q_1
*  and the first row of Q_2.
*  =====================================================================
*
*     .. Parameters ..
*
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DTYPE_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
*
      INTEGER            COL, I, IBUF, ICTXT, IIQ, IIZ1, IIZ2, IQCOL,
     $                   IQROW, IZ, IZ1, IZ1COL, IZ1ROW, IZ2, IZ2COL,
     $                   IZ2ROW, J, JJQ, JJZ1, JJZ2, MYCOL, MYROW, N2,
     $                   NB, NBLOC, NPCOL, NPROW, NQ1, NQ2, ZSIZ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DCOPY, DGEBR2D, DGEBS2D,
     $                   DGERV2D, DGESD2D, INFOG2L
*     ..
*     .. External Functions ..
      INTEGER            NUMROC
      EXTERNAL           NUMROC
*     ..

      REAL(8), POINTER :: TEMP(:)
      include 'mpif.h'

*     .. Executable Statements ..
*
*       This is just to keep ftnchek and toolpack/1 happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DTYPE_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      ICTXT = DESCQ( CTXT_ )
      NB = DESCQ( NB_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      CALL INFOG2L( ID, ID, DESCQ, NPROW, NPCOL, MYROW, MYCOL, IIQ, JJQ,
     $              IQROW, IQCOL )
      N2 = N - N1

      allocate( TEMP(1:N) )
      TEMP(1:N) = 0D0
*
*     Form z1 which consist of the last row of Q1
*
      CALL INFOG2L( IQ-1+( ID+N1-1 ), JQ-1+ID, DESCQ, NPROW, NPCOL,
     $              MYROW, MYCOL, IIZ1, JJZ1, IZ1ROW, IZ1COL )
      NQ1 = NUMROC( N1, NB, MYCOL, IZ1COL, NPCOL )
      IF( ( MYROW.EQ.IZ1ROW ) .AND. ( NQ1.NE.0 ) ) THEN
         CALL DCOPY( NQ1, Q( IIZ1, JJZ1 ), LDQ, WORK, 1 )
         I = MOD( MYCOL + NPCOL - MOD( IZ1COL, NPCOL ), NPCOL )
         NBLOC = ( NQ1-1 ) / NB + 1
!$OMP PARALLEL DO PRIVATE(J,ZSIZ,IZ1,IZ)
         DO J = 1, NBLOC
            IZ1 = (J-1)*NB
            IZ = I*NB + 1 + (J-1)*NB*NPCOL
            ZSIZ = MIN( NB, NQ1-IZ1 )
            IF ( ZSIZ > 0 ) THEN
               CALL DCOPY( ZSIZ, WORK(1+IZ1), 1, TEMP(IZ), 1)
            END IF
         END DO
!$OMP END PARALLEL DO
      END IF
*
*     Form z2 which consist of the first row of Q2
*
      CALL INFOG2L( IQ-1+( ID+N1 ), JQ-1+( ID+N1 ), DESCQ, NPROW, NPCOL,
     $              MYROW, MYCOL, IIZ2, JJZ2, IZ2ROW, IZ2COL )
      NQ2 = NUMROC( N2, NB, MYCOL, IZ2COL, NPCOL )
      IF( ( MYROW.EQ.IZ2ROW ) .AND. ( NQ2.NE.0 ) ) THEN
         CALL DCOPY( NQ2, Q( IIZ2, JJZ2 ), LDQ, WORK, 1 )
         I = MOD( MYCOL + NPCOL - MOD( IZ2COL, NPCOL ), NPCOL )
         NBLOC = ( NQ2-1 ) / NB + 1
!$OMP PARALLEL DO PRIVATE(J,ZSIZ,IZ2,IZ)
         DO J = 1, NBLOC
            IZ2 = (J-1)*NB
            IZ = NB*I + N1 + 1 + (J-1)*NB*NPCOL
            ZSIZ = MIN( NB, NQ2-IZ2 )
            IF ( ZSIZ > 0 ) THEN
              CALL DCOPY( ZSIZ, WORK(1+IZ2), 1, TEMP(IZ), 1)
            END IF
         END DO
!$OMP END PARALLEL DO
      END IF
*
*     proc(IQROW,IQCOL) broadcast Z=(Z1,Z2)
*
      call MPI_Allreduce( TEMP( 1 ), Z( 1 ), N, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MPI_COMM_EIGEN, J)

      deallocate( TEMP )
*
      RETURN
*
*     End of PDLAEDZ
*
      END
