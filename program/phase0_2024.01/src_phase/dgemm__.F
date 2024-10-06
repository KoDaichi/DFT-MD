C *******************************
C *                             *
C *  Coded by NEC CORPORATION.  *
C *                             *
C *******************************
C
      SUBROUTINE dgemm__
     &(TRANSA,TRANSB,NM,NL,NN,ALPHA,A,LMA,B,LNB,BETA,C,LMC)
      IMPLICIT NONE
C
C ::: ARGUMENTS :::
      CHARACTER TRANSA,TRANSB
      INTEGER  LMA,NM,NN,LNB,NL,LMC
      DOUBLE PRECISION A(LMA,*),B(LNB,*),C(LMC,*),ALPHA,BETA
C***********************************************************************
C
      CHARACTER TRANSAW,TRANSBW
      INTEGER MODE
      INTEGER MPANT,MMDNT,IPAR,IST,MM
      INTEGER NT,isw,ierr
      CHARACTER*80 F_RSVTASK
      ierr=0
      isw=1
C
CCCCC SETUP
C
      TRANSAW=TRANSA
      TRANSBW=TRANSB
      IF (TRANSA.EQ.'n') THEN
         TRANSAW='N'
      ELSE IF (TRANSA.EQ.'t') THEN
         TRANSAW='T'
      ELSE IF (TRANSA.EQ.'c') THEN
         TRANSAW='T'
      ELSE IF (TRANSA.EQ.'C') THEN
         TRANSAW='T'
      ENDIF
      IF (TRANSB.EQ.'n') THEN
         TRANSBW='N'
      ELSE IF (TRANSB.EQ.'t') THEN
         TRANSBW='T'
      ELSE IF (TRANSB.EQ.'c') THEN
         TRANSBW='T'
      ELSE IF (TRANSB.EQ.'C') THEN
         TRANSBW='T'
      ENDIF
C
C     GET NUMBER OF TASK
C
      CALL GETENV( "F_RSVTASK", F_RSVTASK)
      READ (F_RSVTASK,'(I8)') NT
      IF (NT.LE.0)THEN
         NT=8
c$$$         NT=1
      ENDIF

CCCCCCCCCCCCCCCCC
C
C   Quickly return if possible
C

      IF( ( NM.EQ.0 ).OR.( NL.EQ.0 ).OR.
     $    ((( ALPHA.EQ.0.0d0 ).OR.( NN.EQ.0 )).AND.( BETA.EQ.1.0d0)))
     $   RETURN



C
CCCCC MULTIPLICATION
C
C.... A * B
      IF(TRANSAW.EQ.'N'.AND.TRANSBW.EQ.'N')THEN
         IF(NM.GE.NL) THEN
            MODE=1
            IF(NL.LT.NT) THEN
               MODE=2
            ENDIF
         ELSE
            MODE=2
            IF(MOD(LNB,2).EQ.0 .AND. MOD(LMC,2).EQ.0) THEN
               MODE=1
            ENDIF
            IF(NM.LT.NT) THEN
               MODE=1
            ENDIF
         ENDIF
C     
         IF(MODE.EQ.1) THEN
            MPANT = NL/NT
            MMDNT = MOD(NL,NT)-1
!CDIR PARDO BY=1
            DO 110 IPAR = 0 , MIN(NT-1,NL-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!          CALL DAM1MM(A,LMA,NM,NN,B(1,IST),LNB,MM,C(1,IST),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,NM,MM,NN,
     &              ALPHA,A,LMA,B(1,IST),LNB,BETA,C(1,IST),LMC)
 110        CONTINUE
         ELSE
            MPANT = NM/NT
            MMDNT = MOD(NM,NT)-1
!CDIR PARDO BY=1
            DO 120 IPAR = 0 , MIN(NT-1,NM-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!          CALL DAM1MM(A(IST,1),LMA,MM,NN,B,LNB,NL,C(IST,1),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,MM,NL,NN,
     &              ALPHA,A(IST,1),LMA,B,LNB,BETA,C(IST,1),LMC)
 120        CONTINUE
         ENDIF
C.... A * B'
      ELSE IF(TRANSAW.EQ.'N'.AND.TRANSBW.EQ.'T')THEN
         IF(NM.GE.NL) THEN
            MODE=1
            IF(NL.LT.NT) THEN
               MODE=2
            ENDIF
         ELSE
            MODE=2
            IF(MOD(LMC,2).EQ.0) THEN
               MODE=1
            ENDIF
            IF(NM.LT.NT) THEN
               MODE=1
            ENDIF
         ENDIF
C     
         IF(MODE.EQ.1) THEN
            MPANT = NL/NT
            MMDNT = MOD(NL,NT)-1
!CDIR PARDO BY=1
            DO 210 IPAR = 0 , MIN(NT-1,NL-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1MT(A,LMA,NM,NN,B(IST,1),LNB,MM,C(1,IST),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,NM,MM,NN,
     &              ALPHA,A,LMA,B(IST,1),LNB,BETA,C(1,IST),LMC)
 210        CONTINUE
         ELSE
            MPANT = NM/NT
            MMDNT = MOD(NM,NT)-1
               IF (MIN(NT-1,NM-1) .GE. 0) THEN
!CDIR PARDO BY=1
            DO 220 IPAR = 0 , MIN(NT-1,NM-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1MT(A(IST,1),LMA,MM,NN,B,LNB,NL,C(IST,1),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,MM,NL,NN,
     &              ALPHA,A(IST,1),LMA,B,LNB,BETA,C(IST,1),LMC)
 220        CONTINUE
               ENDIF
         ENDIF
C.... A' * B
      ELSE IF(TRANSAW.EQ.'T'.AND.TRANSBW.EQ.'N')THEN
         IF(NM.GE.NL) THEN
            MODE=1
            IF(MOD(LMA,2).EQ.0 .AND. MOD(LNB,2).EQ.1) THEN
               MODE=2
            ENDIF
            IF(NL.LT.NT) THEN
               MODE=2
            ENDIF
         ELSE
            MODE=2
            IF(MOD(LNB,2).EQ.0 .AND. MOD(LMA,2).EQ.1) THEN
               MODE=1
            ENDIF
            IF(NM.LT.NT) THEN
               MODE=1
            ENDIF
         ENDIF
C     
         IF(MODE.EQ.1) THEN
            MPANT = NL/NT
            MMDNT = MOD(NL,NT)-1
!CDIR PARDO BY=1
            DO 310 IPAR = 0 , MIN(NT-1,NL-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1TM(A,LMA,NM,NN,B(1,IST),LNB,MM,C(1,IST),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,NM,MM,NN,
     &              ALPHA,A,LMA,B(1,IST),LNB,BETA,C(1,IST),LMC)
 310        CONTINUE
         ELSE
            MPANT = NM/NT
            MMDNT = MOD(NM,NT)-1
!CDIR PARDO BY=1
            DO 320 IPAR = 0 , MIN(NT-1,NM-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1TM(A(1,IST),LMA,MM,NN,B,LNB,NL,C(IST,1),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,MM,NL,NN,
     &              ALPHA,A(1,IST),LMA,B,LNB,BETA,C(IST,1),LMC)
 320        CONTINUE
         ENDIF
C.... A' * B'
      ELSE IF(TRANSAW.EQ.'T'.AND.TRANSBW.EQ.'T')THEN
         IF(NM.GE.NL) THEN
            MODE=1
            IF(MOD(LMA,2).EQ.0 .AND. MOD(LMC,2).EQ.1) THEN
               MODE=2
            ENDIF
            IF(NL.LT.NT) THEN
               MODE=2
            ENDIF
         ELSE
            MODE=2
            IF(MOD(LMA,2).EQ.1 .AND. MOD(LMC,2).EQ.0) THEN
               MODE=1
            ENDIF
            IF(NM.LT.NT) THEN
               MODE=1
            ENDIF
         ENDIF
C     
         IF(MODE.EQ.1) THEN
            MPANT = NL/NT
            MMDNT = MOD(NL,NT)-1
!CDIR PARDO BY=1
            DO 410 IPAR = 0 , MIN(NT-1,NL-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1TT(A,LMA,NM,NN,B(IST,1),LNB,MM,C(1,IST),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,NM,MM,NN,
     &              ALPHA,A,LMA,B(IST,1),LNB,BETA,C(1,IST),LMC)
 410        CONTINUE
         ELSE
            MPANT = NM/NT
            MMDNT = MOD(NM,NT)-1
!CDIR PARDO BY=1
            DO 420 IPAR = 0 , MIN(NT-1,NM-1)
               IF (IPAR.LE.MMDNT) THEN
                  MM = MPANT+1
               ELSE
                  MM = MPANT
               ENDIF
               IST = IPAR*MPANT + MIN(MMDNT+1,IPAR) + 1
!               CALL DAM1TT(A(1,IST),LMA,MM,NN,B,LNB,NL,C(IST,1),LMC,ISW,IERR)
               CALL DGEMM
     &              (TRANSA,TRANSB,MM,NL,NN,
     &              ALPHA,A(1,IST),LMA,B,LNB,BETA,C(IST,1),LMC)
 420        CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
