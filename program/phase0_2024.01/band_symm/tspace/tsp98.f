C SUBROUTINE TSPACE ====*====3====*====4====*====5====*====6====*====7
C
C    PREPARE THE ROTATION OPERATIOS IT(3,48)
C    AND MULTIPLICATION TABLE IM(48,48)
C    INVERSE ELEMENTS ARE PREPARED IN IV(48)
C           IF(IL.LE.0) THEN D6H GROUP
C           IF(IL.GE.1) THEN OH GROUP
C
C                MODIFIED AT 1988/10/10
C                BY AKIRA YANASE
C                CALCULATION FOR ROTATION MTRIX FOR S=1/2(SPINMC)
C                AND MULTIPRICATION TABLE IN SPIN SPACE(SPNMLT)
C                ARE CALCULATED HERE.
C
C                MODIFIED AT 1994/11/29
C                MODIFIRD AT !995/09/06
C                BY AKIRA YANASE
C                FOR V4.1                 
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSPACE(IILL)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG1/,/SPG2/,IMAT,IC,IH,IK
      INTEGER II(2),IJ(3),IMAT(2,9)
      INTEGER IC(3,24),IH(3,12),IK(3,12)
      DATA IMAT/
     &  -1,1, 9,9, 0,-1, -1,0, 9,9,
     &   1,0, 0,1, 9,9,   1,-1/
      DATA IC/
     &   1, 2, 3,  1,-2,-3, -1, 2,-3, -1,-2, 3,
     &   3, 1, 2, -3, 1,-2, -3,-1, 2,  3,-1,-2,
     &   2, 3, 1,  2,-3,-1, -2, 3,-1, -2,-3, 1,
     &   2, 1,-3, -2,-1,-3,  3,-2, 1, -1, 3, 2,
     &  -3,-2,-1, -1,-3,-2,  1,-3, 2,  3, 2,-1,
     &  -2, 1, 3,  1, 3,-2, -3, 2, 1,  2,-1, 3/
      DATA IH/
     &   1, 2, 3,  4, 1, 3, -2, 4, 3, -1,-2, 3,
     &  -4,-1, 3,  2,-4, 3, -4, 2,-3,  1, 4,-3,
     &  -2,-1,-3,  4,-2,-3, -1,-4,-3,  2, 1,-3/
      DATA IK/
     &   1, 2, 3, -2, 4, 3, -4, 1, 3, -1,-2, 3,
     &   2,-4, 3,  4,-1, 3, -1, 4,-3,  4,-2,-3,
     &  -2,-1,-3,  1,-4,-3, -4, 2,-3,  2, 1,-3/
      IL=IILL
      IF(-1.GT.IL.OR.IL.GT.4) THEN
         WRITE(6,900) IL
  900 FORMAT(' IL FOR LATTICE TYPE SHOULD BE -1.GE.IL.LE.4'
     &      /' HOWEVER YOUR IL=',I6,' THEN STOP HERE')
         STOP
      END IF
C-----------------------------------
      WRITE(6,*) ' ----- WELCOME TO TSPACE V4.1 1995/09/06 -----'
      CALL SPINMC
      CALL SPNMLT
C     CALL TSPHDS
C
C------------------------------------
      IF(IL.GT.0) THEN
         DO 1 I=1,24
         DO 1 K=1,3
         IT(K,I)=IC(K,I)
         IT(K,I+24)=IC(K,I)*(-1)
    1    CONTINUE
      ELSE
         DO 2 I=1,12
         DO 2 K=1,3
         IT(K,I)=IH(K,I)
         IT(K,I+24)=IK(K,I)
         IT(K,I+12)=IH(K,I)*(-1)
    2    IT(K,I+36)=IK(K,I)*(-1)
      ENDIF
      DO 3 I=1,24
      DO 3 J=1,24
      DO 4 K=1,3
         IF(IABS(IT(K,I)).EQ.4) THEN
            N1=IT(1,J)+5
            N2=IT(2,J)+5
            DO 5 L=1,2
            II(L)=IMAT(L,N1)-IMAT(L,N2)
            IF(IT(K,I).EQ.-4) II(L)=-II(L)
    5       CONTINUE
            DO 6 L=1,9
            IF(II(1).EQ.IMAT(1,L).AND.
     &      II(2).EQ.IMAT(2,L)) GO TO 7
    6       CONTINUE
    7       N=L-5
         ELSE
            M=IABS(IT(K,I))
            N=IT(M,J)
            IF(IT(K,I).LT.0) N=-N
         ENDIF
         IJ(K)=N
    4 CONTINUE
      DO 8 K=1,24
      IF(IJ(1).EQ.IT(1,K).AND.
     &   IJ(2).EQ.IT(2,K).AND.
     &   IJ(3).EQ.IT(3,K)) GO TO 9
    8 CONTINUE
    9 IM(I,J)=K
    3 CONTINUE
      IF(IL.GT.0)THEN
      NN=48
      DO 19 I=1,24
      DO 19 J=1,24
      IM(I,J+24)=IM(I,J)+24
      IM(I+24,J)=IM(I,J)+24
   19 IM(I+24,J+24)=IM(I,J)
      ELSE
      NN=24
      ENDIF
      DO 10 I=1,NN
      DO 11 J=1,NN
      IF(IM(I,J).EQ.1) GO TO 12
   11 CONTINUE
   12 IV(I)=J
   10 CONTINUE
      NG=1
      IG(1)=1
      DO 13 I=1,3
      JV(1,I,1)=0
   13 JV(2,I,1)=1
      RETURN
      END
C SUBROUTINE SPINMC ====*====3====*====4====*====5====*====6====*====7
C
C      ROTATION MATRIIS FOR SPIN SPACE ARE PREPARED
C      FOR CUBIC OR HEXAGONAL OPERATIONS
C
C      1988.10.10 :  A. YANASE
C
C     CALL TSPACE(1)
C     CALL SPINMC
C     CALL SPNMLT
C     STOP
C     END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE SPINMC
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SN
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG2/,/SPG6/
      DIMENSION  EULER(3,24)
      INTEGER ICREAL(2,3,24)
      INTEGER IC(16),ICM(16,2,2),K1M(2,2),K2M(2,2),K3M(2,2),ISM(2,2)
C
      IF(IL.GE.1) THEN
         CALL EULERC(ICREAL,EULER)
         NN=24
      ELSE
         CALL EULERH(ICREAL,EULER)
         NN=12
      END IF
      IPR=0
      DO 31 J=1,2
      DO 31 K=1,2
         MA=3-2*J
         MB=3-2*K
         CALL TSRMHI(1,MA,MB,K1,K2,K3,IS,IC,IPR)
         DO 32 II=1,16
           ICM(II,J,K)=IC(II)
   32    CONTINUE
         K1M(J,K)=K1
         K2M(J,K)=K2
         K3M(J,K)=K3
         ISM(J,K)=IS
   31 CONTINUE
      DO 10 I=1,NN
C       WRITE(6,660) (EULER(K,I),K=1,3)
C 660   FORMAT(1H ,3F10.5)
        CA=COS(EULER(1,I)*0.5D0)
        CC=COS(EULER(3,I)*0.5D0)
        SA=SIN(EULER(1,I)*0.5D0)
        SC=SIN(EULER(3,I)*0.5D0)
        CB=COS(EULER(2,I))
        W1=COS(EULER(2,I)*0.5D0)
        W3=SIN(EULER(2,I)*0.5D0)
        DO 20 J=1,2
        DO 21 K=1,2
          MA=3-2*J
          MB=3-2*K
          WW=0.0
          DO 22 II=1,12
            IF(ICM(II,J,K).NE.0) THEN
            WW=WW+ICM(II,J,K)*(CB**(II-1))
            END IF
   22     CONTINUE
          WW=((WW*K1M(J,K))/K2M(J,K))*SQRT(DBLE(K3M(J,K)))
          IF(ISM(J,K).EQ.1) THEN
             WW=WW*W1
          ELSE
             WW=WW*W3
          END IF
          SN(J,K,I)=DCMPLX(CA,-MA*SA)*DCMPLX(CC,-MB*SC)*WW
          IF(IL.LE.0) THEN
             SN(J,K,I+12)=SN(J,K,I)
          END IF
   21   CONTINUE
C       WRITE(6,600) I,MA,(SN(J,K,I),K=1,2)
C 600   FORMAT(2I3,'/2',2(' (',2F10.7,')'))
   20   CONTINUE
C     WRITE(6,600)
   10 CONTINUE
      RETURN
      END
C SUBROUTINE SPNMLT ====*====3====*====4====*====5====*====6====*====7
C
C      MULITIPLICATIO TABLE OF SPIN HALF MATRIX
C      FOR CUBIC OR HEXAGONAL OPERATIONS
C
C      1988.10.10 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SPNMLT
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 SN
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG2/,/SPG6/
      COMPLEX*16 WW,WM(2,2)
      IF(IL.GE.1) NN=24
      IF(IL.LE.0) NN=12
      DO 10 I=1,NN
      DO 20 J=1,NN
C---
         DO 31 K1=1,2
         DO 33 K2=1,2
            WW=0.0
            DO 32 K=1,2
            WW=WW+SN(K1,K,I)*SN(K,K2,J)
   32       CONTINUE
         WM(K1,K2)=WW
   33    CONTINUE
C        WRITE(6,600) I,J,3-2*K1,(WM(K1,K2),K2=1,2)
C 600   FORMAT(2I3,I2,'/2',2(' (',2F10.7,')'))
   31    CONTINUE
C        WRITE(6,600)
C---
      DO 41 K=1,NN
         INIT=0
      DO 42 K1=1,2
      DO 43 K2=1,2
         IF(ABS(WM(K1,K2)).GT.0.000001D0) THEN
           IF(ABS(WM(K1,K2)-SN(K1,K2,K)).LT.0.000001D0) THEN
              INDS=1
           ELSE IF(ABS(WM(K1,K2)+SN(K1,K2,K)).LT.0.000001D0) THEN
              INDS=-1
           ELSE
              INDS=0
           END IF
           IF(INDS.EQ.0) GO TO 41
           IF(INIT.EQ.1.AND.INDS.NE.INDO) GO TO 41
           INDO=INDS
           INIT=1
         END IF
   43 CONTINUE
   42 CONTINUE
         IDF(I,J)=INDS*K
         IF(IL.LE.0) THEN
             IDF(I+12,J)=INDS*(K+12)
             IDF(I+12,J+12)=IDF(I,J)
             IDF(I,J+12)=INDS*(K+12)
         END IF
         GO TO 20
   41 CONTINUE
         WRITE(6,900) I,J
  900    FORMAT(' STOP AT 41 IN SPNMLT FOR',2I5)
         STOP
   20 CONTINUE
C     WRITE(6,610) (IDF(I,J),J=1,NN)
C 610 FORMAT(1H ,24I3)
   10 CONTINUE
      RETURN
      END
C SUBROUTINE EULERC ====*====3====*====4====*====5====*====6====*====7
C
C      EULERIAN ANGLES OF CUBIC OPERATION.
C
C      1983.7.11. :  N. HAMADA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE EULERC(ICREAL,EULER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE IC1,IC2
C
C      INTEGER ICREAL(2,3,24),ICR(2,3,24)
      INTEGER ICREAL(2,3,24)
      DIMENSION EULER(3,24)
      INTEGER IC1(2,3,12),IC2(2,3,12)
      DATA IC1
     &          /  0,  1,  0,  1,  0,  1,
     &            -1,  2,  1,  1,  1,  2,
     &             0,  1,  1,  1,  0,  1,
     &             0,  1,  0,  1,  1,  1,
     &             0,  1,  1,  2,  1,  2,
     &             0,  1, -1,  2,  1,  2,
C
     &             0,  1, -1,  2, -1,  2,
     &             0,  1,  1,  2, -1,  2,
     &            -1,  2, -1,  2,  0,  1,
     &            -1,  2,  1,  2,  0,  1,
     &             1,  2,  1,  2,  0,  1,
     &             1,  2, -1,  2,  0,  1/
C
      DATA IC2
     &          /  0,  1,  1,  1,  1,  2,
     &             1,  2,  1,  1,  0,  1,
     &             1,  1, -1,  2,  0,  1,
     &             1,  2,  1,  2,  1,  2,
     &             1,  1,  1,  2,  0,  1,
     &             1,  2, -1,  2,  1,  2,
C
     &            -1,  2,  1,  2,  1,  2,
     &             0,  1,  1,  2,  0,  1,
     &             1,  2,  0,  1,  0,  1,
     &            -1,  2, -1,  2,  1,  2,
     &             0,  1, -1,  2,  0,  1,
     &            -1,  2,  0,  1,  0,  1  /
      DATA PAI / 3.14159265358979312D0 /
      DO 10 I=1,24
      DO 10 J=1,3
        DO 20 K=1,2
        IF(I.LE.12) THEN
           ICREAL(K,J,I)=IC1(K,J,I)
        ELSE
           ICREAL(K,J,I)=IC2(K,J,I-12)
        END IF
   20   CONTINUE
        EULER(J,I)=PAI*ICREAL(1,J,I)/ICREAL(2,J,I)
   10 CONTINUE
      RETURN
      END
C SUBROUTINE EULERH ====*====3====*====4====*====5====*====6====*====7
C
C     EULERIAN ANGLES OF HEXAGONAL OPERATION
C
C      1983.7.11. :  N. HAMADA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE EULERH(ICREAL,EULER)
C
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE ICR
C
      DIMENSION EULER(3,12)
      INTEGER ICREAL(2,3,12)
      INTEGER ICR(2,3,12)
      DATA ICR
     &          /  0,  1,  0,  1,  0,  1,
     &             1,  3,  0,  1,  0,  1,
     &             2,  3,  0,  1,  0,  1,
     &             1,  1,  0,  1,  0,  1,
     &            -2,  3,  0,  1,  0,  1,
     &            -1,  3,  0,  1,  0,  1,
     &             0,  1,  1,  1,  0,  1,
     &             2,  3,  1,  1, -2,  3,
     &            -2,  3,  1,  1,  2,  3,
     &            -1,  2,  1,  1,  1,  2,
     &             1,  6,  1,  1, -1,  6,
     &             5,  6,  1,  1, -5,  6 /
      DATA PAI / 3.14159265358979312D0 /
      DO 10 I=1,12
      DO 10 J=1,3
        DO 20 K=1,2
          ICREAL(K,J,I)=ICR(K,J,I)
   20   CONTINUE
        H1=DBLE(ICREAL(1,J,I))
        H2=DBLE(ICREAL(2,J,I))
        EULER(J,I)=PAI*H1/H2
   10 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY26 ====*====3====*====4====*====5====*====6====*====7
C
C       CALCULATION FOR THE FACTOR SYSTEM
C       CALLED FROM TSIREP
C           PHAI(I,J))=EXP(I*2PAI*(IFA(I,J)/IFC)
C           EXTRA FACTOR DOUE TO THE DOUBLE GROUP WAS STORED
C           AS THE SIGN OF IDF(I,J) IN SPNMLT
C
C                             MODIFIED AT 1988/10/10
C                             BY AKIRA YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE ZZZY26
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 SN
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48),
     &   IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG6/
      KW=1
      DO 1 J=1,MG
      DO 2 I=1,MG
      IW=0
      JW=1
      IF(IZ.NE.0) THEN
         DO 3 K=1,3
         IW=JK(K,I)*JV(1,K,JG(J))*JW
     &             +JV(2,K,JG(J))*IW
         JW=JV(2,K,JG(J))*JW
         CALL ZZZY24(IW,JW)
    3    CONTINUE
         IF(IW.LT.0) IW=IW+JW
         IF(IW.LT.0) IW=IW+JW
         IF(IW.LT.0) IW=IW+JW
         IF(IW.LT.0) IW=IW+JW
      ENDIF
      IF(IDOUB.EQ.1) THEN
         IW1=MOD(IG(JG(I))-1,24)+1
         IW2=MOD(IG(JG(J))-1,24)+1
         IF(IDF(IW1,IW2).LT.0) THEN
            IW=IW*2+JW
            JW=JW*2
            CALL ZZZY24(IW,JW)
         ENDIF
      ENDIF
      IW=MOD(IW,JW)
      IFA(I,J)=IW*256+JW
      KW1=KW
      CALL ZZZY24(KW,JW)
      KW=KW1*JW
    2 CONTINUE
    1 CONTINUE
      IFC=KW
      DO 4 J=1,MG
      DO 5 I=1,MG
         IW=IFA(I,J)/256
         JW=MOD(IFA(I,J),256)
         IFA(I,J)=IW*(KW/JW)
    5 CONTINUE
    4 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY35 ====*====3====*====4====*====5====*====6====*====7
C
C       FIND THE TIME REVERSAL ELEMENTS AND  CALCULATION OF HERING'S
C       SUM WHICH SUM IS STORED IN  COMMON/SPG4/-------,IATR,---
C       RESULTS CAN BE SEEN BY CALLING TSTRDS(ENTRY IN TSPRDS)
C
C           CALLED FROM TSIREP.
C           EXTRA FACTOR DOUE TO THE DOUBLE GROUP WAS STORED
C           AS THE SIGN OF IDF(I,J) IN SPNMLT
C
C                             MODIFIED AT 1988/10/10
C                             BY AKIRA YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY35
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW(12),C,SN
      LOGICAL LOGIC
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48),
     &   IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/SPG6/
      INTEGER KA(3),KC(3)
      DO 31 IO=1,NR
   31 CW(IO)=0.D0
      MTRG=0
      DO 1 I=1,NG
      ITA=IV(IG(I))
      ITB=IG(I)
      IF(IL.LE.0) ITA=ITA+24
      IND=0
      DO 2 J=1,3
      MMA=IT(J,ITA)
      IF(IABS(MMA).NE.4) GO TO 11
      KW=KB(1)+KB(2)
      IF(MMA.EQ.-4) KA(J)=-KB(J)+KW
      IF(MMA.EQ.4) KA(J)=-KB(J)-KW
      GO TO 12
   11 CONTINUE
      MB=IABS(MMA)
      IF(MMA.GT.0) KA(J)=-KB(J)-KB(MB)
      IF(MMA.LT.0) KA(J)=-KB(J)+KB(MB)
   12 CONTINUE
      IF(MOD(IABS(KA(J)),ICB).NE.0) GO TO 1
      KA(J)=KA(J)/ICB
      KC(J)=IABS(KA(J))
      IF(KA(J).NE.0) IND=1
    2 CONTINUE
      IF(IND.EQ.0) GO TO 3
      IF(IL.EQ.0) GO TO 3
      IF(IL.EQ.1) GO TO 3
      IF(IL.EQ.2) GO TO 7
      IF(IL.EQ.3) GO TO 4
      IF(IL.EQ.-1) GO TO 9
      IF(MOD(KC(1)+KC(2),2).NE.0) GO TO 1
      GO TO 3
    9 KSUM=-KA(1)+KA(2)+KA(3)
      IF(MOD(IABS(KSUM),3).NE.0) GO TO 1
      GO TO 3
    7 CONTINUE
      IF(MOD(KC(1),2).NE.MOD(KC(2),2)) GO TO 1
      IF(MOD(KC(1),2).NE.MOD(KC(3),2)) GO TO 1
      GO TO 3
    4 CONTINUE
      IF(MOD(KC(1)+KC(2)+KC(3),2).NE.0) GO TO 1
    3 ITC=IM(IG(I),IG(I))
      DO 5 J=1,MG
      IF(ITC.EQ.IG(JG(J))) GO TO 6
    5 CONTINUE
    6 ITD=MOD(ITB-1,24)+1
      LOGIC=IDF(ITD,ITD).LT.0.AND.IDOUB.EQ.1
      MTRG=MTRG+1
      JTRG(4,MTRG)=I
      DO 21 K=1,3
      JTRG(K,MTRG)=KA(K)
   21 CONTINUE
      IW=0
      JW=1
      DO 8 K=1,3
      KD=JV(1,K,I)
      KE=JV(2,K,I)
      IW=IW*KE+JW*KA(K)*KD
      JW=JW*KE
      CALL ZZZY24(IW,JW)
    8 CONTINUE
   40 IF(IW.GE.0) GO TO 41
      IW=IW+JW
      GO TO 40
   41 IF(IW.LT.JW) GO TO 42
      IW=IW-JW
      GO TO 41
   42 IP=IW*(24/JW)
      ITRC(1,MTRG)=ITC
      ITRC(2,MTRG)=IW
      ITRC(3,MTRG)=JW
      DO 32 IO=1,NR
      C=CR(J,IO)
      IF(LOGIC) C=-C
      CALL ZZZY27(C,IP)
      CW(IO)=CW(IO)+C
   32 CONTINUE
    1 CONTINUE
      DO 33 IO=1,NR
      AW=DBLE(CW(IO))
      IF(AW.GE.0.D0) IATR(IO)=AW+0.5
      IF(AW.LT.0.D0) IATR(IO)=AW-0.5
C   IF NO TIME REVERSAL ELEMENT THEN
      IF(MTRG.EQ.0)  IATR(IO)=99
   33 CONTINUE
      RETURN
      END
C SUBROUTINE TSGENR ====*====3====*====4====*====5====*====6====*====7
C
C      GENERATOR OF THE SPACE GROUP IS STORED IN COMMON/SPG2/
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSGENR(JA,JB)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
      INTEGER JB(2,3)
      DO 1 I=1,NG
      IF(IG(I).EQ.JA) RETURN
    1 CONTINUE
      NG=NG+1
      IG(NG)=JA
      DO 2 I=1,3
      JV(1,I,NG)=JB(1,I)
      CALL CHKDNM(JB(2,I))
    2 JV(2,I,NG)=JB(2,I)
      RETURN
      END
C SUBROUTINE ZZZY22 ====*====3====*====4====*====5====*====6====*====7
C
C       VECTOR GIVEN BY JB IS TRANSFORMED THROUGH THE OPERATION
C       GIVEN BY IA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY22(IA,JB)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG1/,/SPG2/
      INTEGER JA(2,3),JB(2,3),JC(2,3),JD(2,3)
      DO 1 I=1,3
      MA=IT(I,IA)
      IF(IABS(MA).EQ.4) THEN
      IW1=JB(1,1)*JB(2,2)-JB(1,2)*JB(2,1)
      IW2=JB(2,1)*JB(2,2)
      IF(MA.EQ.-4) IW1=-IW1
      IF(IW1.LT.0) IW1=IW1+IW2
      CALL ZZZY24(IW1,IW2)
      IW1=MOD(IW1,IW2)
      JA(1,I)=IW1
      JA(2,I)=IW2
      ELSE
      MB=IABS(MA)
      JA(1,I)=JB(1,MB)
      IF(MA.LT.0.AND.JB(1,MB).NE.0) JA(1,I)=JB(2,MB)-JB(1,MB)
      JA(2,I)=JB(2,MB)
      ENDIF
    1 CONTINUE
      DO 2 I=1,3
      JB(1,I)=JA(1,I)
    2 JB(2,I)=JA(2,I)
      RETURN
C ENTRY ZZZY23 ====2====*====3====*====4====*====5====*====6====*====7
C
C         (VECTOR JC) + (VECTOR JD)*IS
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY ZZZY23(JC,JD,IS)
      DO 3 I=1,3
      IW1=JC(1,I)*JD(2,I)+JC(2,I)*JD(1,I)*IS
      IW2=JC(2,I)*JD(2,I)
      IF(IW1.LT.0) IW1=IW1+IW2
      IW1=MOD(IW1,IW2)
      CALL ZZZY24(IW1,IW2)
      JC(1,I)=IW1
    3 JC(2,I)=IW2
      RETURN
      END
C SUBROUTINE ZZZY24 ====*====3====*====4====*====5====*====6====*====7
C
C        REDUCTION OF IW/JW
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY24(IW,JW)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      IF(IW.EQ.0) GO TO 2
      I=IABS(IW)
      J=IABS(JW)
    1 K=MOD(J,I)
      IF(K.EQ.0) GO TO 3
      J=I
      I=K
      GO TO 1
    3 IW=IW/I
      JW=JW/I
      RETURN
    2 IW=0
      JW=1
      RETURN
      END
C SUBROUTINE TSPHDS ====*====3====*====4====*====5====*====6====*====7
C
C   PRINT OUT THE ROTATION MATRIC SPIN HALF
C   FOR CUBIC OR HEXAGONAL OPERATION
C             1988/10/10 AKIRA YANASE
C
C---*----1----*----0----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSPHDS
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 SN
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG1/,/SPG2/,/SPG6/,XYZ,CMN,HMN,AC
      CHARACTER*4 XYZ(4),CMN(24)
      CHARACTER*4 HMN(12)
      CHARACTER A*1,B*4,GN(2,3)*4
      DIMENSION  EULER(3,24)
      INTEGER ICREAL(2,3,24)
      CHARACTER*1 AC(3)
	  DATA XYZ/'X','Y','Z','W'/
      DATA CMN/'E','C2X','C2Y','C2Z',
     &   'C31+','C32+','C33+','C34+',
     &   'C31-','C32-','C33-','C34-',
     &   'C2A','C2B','C2C','C2D','C2E','C2F',
     &   'C4X+','C4Y+','C4Z+','C4X-','C4Y-','C4Z-'/
	  DATA HMN/
     &  'E','C6+','C3+','C2','C3-','C6-',
     &  'C211','C221','C231','C212','C222','C232'/
	  DATA AC/' ','I','-'/
      IPS=0
      WRITE(6,615)
  615 FORMAT(1H //' TABLE OF OPERATION ')
      WRITE(6,609)
      WRITE(6,605)
  609   FORMAT(4X,'NAME ','JONES REP.   ','EULER ANGLE '
     &      ,14X,'ROT. MAT, OF S=1/2')
  605   FORMAT(18X,'ALPHA',2X,'BETA',3X,'GAMMA',2X
     &          ,' MA',5X,'MB=1/2',12X,'-1/2')
      IF(IL.GE.1) THEN
         CALL EULERC(ICREAL,EULER)
         MNG=24
      ELSE
         CALL EULERH(ICREAL,EULER)
         MNG=12
      END IF
      DO 1 I=1,MNG
         ITA=I
         A=AC(1)
      IF(IL.LE.0) B=HMN(MOD(ITA-1,12)+1)
      IF(IL.GT.0) B=CMN(MOD(ITA-1,24)+1)
      DO 2 J=1,3
      MA=IT(J,ITA)
      MB=IABS(MA)
      IF(MA.GT.0) GN(1,J)=AC(1)
      IF(MA.LT.0) GN(1,J)=AC(3)
      GN(2,J)=XYZ(MB)
    2 CONTINUE
        WRITE(6,606) ITA,A,B,GN,(EULER(K,I),K=1,3),(SN(1,K,I),K=1,2)
        WRITE(6,607) (SN(2,K,I),K=1,2)
C       WRITE(6,608)
  606   FORMAT(I3,A1,A4,3(1X,2A1),3F7.4,' 1/2',2(' (',2F7.4,')'))
  607   FORMAT(                    38X, '-1/2',2(' (',2F7.4,')'))
C 608   FORMAT(1X)
    1 CONTINUE
      IF(IL.LE.0) WRITE(6,663)
  663 FORMAT(1H //' GROUP TABLE FOR D6 IN SPIN SPACE')
      IF(IL.GT.0) WRITE(6,662)
  662 FORMAT(1H //' GROUP TABLE FOR O-GROUP IN SPIN SPACE')
      DO 7 I=1,MNG
    7 WRITE(6,660) (IDF(I,J),J=1,MNG)
  660 FORMAT(24I3)
      RETURN
      END
C SUBROUTINE ZZZY34 ====*====3====*====4====*====5====*====6====*====7
C
C            WAVE VECTOR KBB IS TRANSFORMED THROUGH
C            THE SPECE GROUP OPERATION GIVEN IN IG(48).
C            WHEN THE TRASFORMED VECTOR IS EQUIVALENT WITH KBB,
C            TRANSFORMATION IS REGISTERED IN JG(48) AS A ELEMENT
C            OF P-K GROUP AND THE RECIPROCAL LATTICE VECTOR IS
C            REGISTERD IN JK(3.48).
C            IZ BECOMES NONE ZERO, WHEN KBB ON THE ZONE BOUNDARY
C            AND GROUP IS NON SYMMORPHIC FOR THIS KBB EFFCTIVELY.
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY34(KBB,ICBB,IDOUBB)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      SAVE /SPG1/,/SPG2/,/SPG3/
      INTEGER KBB(3),KA(3),KC(3)
      IDOUB=IDOUBB
      ICB=ICBB
      DO 1 I=1,3
    1 KB(I)=KBB(I)
      MG=0
      INDA=0
      INDB=0
      DO 2 I=1,NG
      IA=IV(IG(I))
      IF(IL.LE.0) IA=IA+24
      INDC=0
      DO 3 J=1,3
      MA=IT(J,IA)
      IF(IABS(MA).EQ.4) KW=KB(1)+KB(2)
      IF(IABS(MA).NE.4) KW=KB(IABS(MA))
      IF(MA.GT.0) KA(J)=KB(J)-KW
      IF(MA.LT.0) KA(J)=KB(J)+KW
      IF(MOD(IABS(KA(J)),ICB).NE.0) GO TO 2
      KA(J)=KA(J)/ICB
      KC(J)=IABS(KA(J))
      IF(KA(J).NE.0) INDC=1
    3 CONTINUE
      IF(INDC.EQ.1) THEN
      IF(IL.EQ.-1) GO TO 4
      IF(IL.EQ.0) GO TO 5
      IF(IL.EQ.1) GO TO 5
      IF(IL.EQ.2) GO TO 6
      IF(IL.EQ.3) GO TO 7
      IF(IL.EQ.4) GO TO 8
    4 IIWW=-KA(1)+KA(2)+KA(3)
      IF(MOD(IABS(IIWW),3).NE.0) GO TO 2
      GO TO 5
    6 IF(MOD(KC(1),2).NE.MOD(KC(2),2)) GO TO 2
      IF(MOD(KC(2),2).NE.MOD(KC(3),2)) GO TO 2
      GO TO 5
    7 IF(MOD(KC(1)+KC(2)+KC(3),2).NE.0) GO TO 2
      GO TO 5
    8 IF(MOD(KC(1)+KC(2),2).NE.0) GO TO 2
    5 INDA=1
      ENDIF
      MG=MG+1
      DO 9 J=1,3
      JK(J,MG)=KA(J)
      IF(JV(1,J,I).NE.0) INDB=1
    9 CONTINUE
      JG(MG)=I
    2 CONTINUE
      IZ=INDA*INDB
      RETURN
      END
C SUBROUTINE TSRMHI ====*====3====*====4====*====5====*====6====*====7
C
C   EXPRESSION OF THE ROTATION MATRIX IN HALF INTEGER J/2 SPACE
C   FOR EULARIAN ANGLE BETA
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSRMHI(J,MA,MB,K1,K2,K3,IS,IC,IPR)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*13 CR(4)
      CHARACTER*2  CD
      COMMON/IPRC/IPA(10,2)
      SAVE /IPRC/,INDA,IPRIM
      INTEGER IPRIM(10)
      DIMENSION IC(16),ID(16)
      DATA IPRIM/2,3,5,7,11,13,17,19,23,29/
      DATA INDA/0/
C     DATA CR/' ','SQRT((1+C)/2)','SQRT(1-C**2)','SQRT((1-C)/2)'/
      DATA CR/' ','COS(BETA/2)','SIN(BETA)','SIN(BETA/2)'/
      IF(J.LT.1.OR.MOD(J,2).EQ.0) THEN
          WRITE(6,900) J
  900 FORMAT(' J SHOLD BE ODD AND POSITIVE FOR TSRMHI'
     &   /' HOWEVER J=',I3,' THEN STOP HERE')
          STOP
      END IF
      IF(MA.LT.-J.OR.J.LT.MA) THEN
          WRITE(6,901) J,J,MA
  901 FORMAT(' MA SHOULD BE -',I2,'.GE.MA.LE.',I2
     &  /' HOWEVER MA=',I3,' THEN STOP HERE')
          STOP
      END IF
      IF(MB.LT.-J.OR.J.LT.MB) THEN
          WRITE(6,902) J,J,MB
  902 FORMAT(' MB SHOULD BE -',I2,'.GE.MB.LE.',I2
     &  /' HOWEVER MB=',I3,' THEN STOP HERE')
          STOP
      END IF
      JAP=(J+MA)/2
      JAM=(J-MA)/2
      JBP=(J+MB)/2
      JBM=(J-MB)/2
      MAB=(MA-MB)/2
      JHI=1
      GO TO 20
C ENTRY TSRMI =====2====*====3====*====4====*====5====*====6====*====7
C     FOR INTEGER J
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY TSRMI(J,MA,MB,K1,K2,K3,IS,IC,IPR)
      IF(J.LT.0) THEN
          WRITE(6,903) J
  903 FORMAT(' J SHOLD NOT BE NEGATIVE FOR TSRMI'
     &   /' HOWEVER J=',I3,' THEN STOP HERE')
          STOP
      END IF
      IF(MA.LT.-J.OR.J.LT.MA) THEN
          WRITE(6,901) J,J,MA
          STOP
      END IF
      IF(MB.LT.-J.OR.J.LT.MB) THEN
          WRITE(6,902) J,J,MB
          STOP
      END IF
      JAP=J+MA
      JAM=J-MA
      JBP=J+MB
      JBM=J-MB
      MAB=MA-MB
      JHI=0
   20 CONTINUE
      ITI=MAX0(0,MAB)+1
      ITA=MIN0(JAP,JBM)+1
      NC=(JAP+JBM-2*(ITI-1))/2
     &  +(2*(ITI-1)-MAB)/2+1
      IF(INDA.EQ.1) GO TO 19
      CALL ZZZY54
      INDA=1
   19 CONTINUE
      DO 2 I=1,10
    2 IPA(I,1)=0
      DO 1 I=1,16
    1 IC(I)=0
      DO 4 ITT=ITI,ITA
      IT=ITT-1
      DO 3 I=1,10
    3 IPA(I,2)=0
      CALL ZZZY51(2,JAP-IT)
      CALL ZZZY51(2,JBM-IT)
      CALL ZZZY51(2,IT)
      CALL ZZZY51(2,IT-MAB)
      IW=1
      JW=1
      DO 6 I=1,10
      IF(ITT.EQ.ITI) GO TO 16
      IF(IPA(I,1).EQ.IPA(I,2)) GO TO 6
      IF(IPA(I,1).LT.IPA(I,2)) GO TO 7
      IW=IW*(IPRIM(I)**(IPA(I,1)-IPA(I,2)))
      GO TO 6
    7 JW=JW*(IPRIM(I)**(IPA(I,2)-IPA(I,1)))
   16 IPA(I,1)=IPA(I,2)
    6 CONTINUE
      IF(JW.EQ.1) GO TO 5
      DO 8 I=1,NC
    8 IC(I)=IC(I)*JW
    5 IW=IW*(1-2*MOD(IT,2))
      I1=JAP+JBM-2*IT
      I2=2*IT-MAB
      J1=I1/2
      J2=I2/2
      CALL ZZZY52(ID,J1,J2)
      DO 9 I=1,NC
    9 IC(I)=IC(I)+IW*ID(I)
    4 CONTINUE
      I1=MOD(I1,2)
      I2=MOD(I2,2)
      IS=I1+I2
      IF(I1.EQ.0.AND.I2.EQ.1) IS=3
      IW=0
      DO 11 I=1,NC
      CALL ZZZY53(IW,IC(I),JW)
      IW=JW
      IF(JW.EQ.1) GO TO 12
   11 CONTINUE
      DO 13 I=1,NC
   13 IC(I)=IC(I)/JW
      DO 17 I=1,10
   18 IF(MOD(JW,IPRIM(I)).NE.0) GO TO 17
      IPA(I,1)=IPA(I,1)-1
      JW=JW/IPRIM(I)
      IF(JW.EQ.1) GO TO 12
      GO TO 18
   17 CONTINUE
   12 CONTINUE
      DO 10 I=1,10
   10 IPA(I,1)=IPA(I,1)*(-2)
      CALL ZZZY51(1,JAP)
      CALL ZZZY51(1,JAM)
      CALL ZZZY51(1,JBP)
      CALL ZZZY51(1,JBM)
      IF(JHI.EQ.1) IPA(1,1)=IPA(1,1)-J+1
      IF(JHI.EQ.0) IPA(1,1)=IPA(1,1)-2*J
      K1=1
      K2=1
      K3=1
      K4=1
      DO 14 I=1,10
      IF(IPA(I,1).EQ.0) GO TO 14
      IF(IPA(I,1).LT.0) GO TO 15
      K1=K1*(IPRIM(I)**(IPA(I,1)/2))
      IF(MOD(IPA(I,1),2).EQ.1) K3=K3*IPRIM(I)
      GO TO 14
   15 K2=K2*(IPRIM(I)**(IABS(IPA(I,1))/2))
      IF(MOD(-IPA(I,1),2).EQ.1) K4=K4*IPRIM(I)
   14 CONTINUE
      IF(K4.EQ.1) GO TO 21
      K3=K3*K4
      K2=K2*K4
   21 CONTINUE
      IF(IPR.NE.1) RETURN
      CD='  '
      IF(JHI.EQ.1) CD='/2'
      IF(JHI.EQ.0.AND.J.GT.8) GO TO 22
      IF(JHI.EQ.1.AND.J.GT.15) GO TO 22
      WRITE(6,600) J,CD,MA,CD,MB,CD,K1,K2,K3,CR(IS+1)
  600 FORMAT(' J=',I3,A2,2X,'MA=',I3,A2,'  MB=',I3,A2,'  (',
     &  I3,'/',I5,')SQRT(',I5,')',A13)
      GO TO 23
   22 CONTINUE
      WRITE(6,601) J,CD,MA,CD,MB,CD,K1,K2,K3,CR(IS+1)
  601 FORMAT(' J=',I3,A2,2X,'MA=',I3,A2,'  MB=',I3,A2,'  (',
     &  I3,'/',I5,')SQRT(',I10,')',A13)
   23 CONTINUE
      WRITE(6,660) (IC(I),I-1,I=1,NC)
  660 FORMAT(2X,I9,'C**',I2,I9,'C**',I2,I9,'C**',I2,
     &           I9,'C**',I2,I9,'C**',I2)
      RETURN
      END
C SUBROUTINE ZZZY51 ====*====3====*====4====*====5====*====6====*====7
C
C        INTEGER J IS RESOLVE TO A PRODUCT OF PRIMITIVE INTEGER
C        RESULT IS STORED IN IPA(10,K).
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY51(K,J)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/IPRC/IPA(10,2)
      SAVE /IPRC/,IPRIM
      INTEGER IPRIM(10)
      DATA IPRIM/2,3,5,7,11,13,17,19,23,29/
      IF(J.LE.1) RETURN
      DO 1 I=2,J
      N=I
      DO 2 IP=1,10
    3 IF(MOD(N,IPRIM(IP)).NE.0) GO TO 2
      IPA(IP,K)=IPA(IP,K)+1
      N=N/IPRIM(IP)
      IF(N.EQ.1) GO TO 1
      GO TO 3
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY52 ====*====3====*====4====*====5====*====6====*====7
C
C     SUM OF COMBINATIONS IN WIGNER'S EXPRESSION IS CALCULATED
C     ENTRY ZZZY54 PREPARE THE COMBINATION IC(16,16).
C     BE CAREFUL!! OCCATIONARY, CALCULATION FOR IC HAS A TRUBLE
C     BY THE VECTORIZATION.
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY52(IB,M,N)
C
        IMPLICIT REAL*8(A-H,O-Z)
      SAVE IC
C
      DIMENSION IC(16,16),IB(16)
      GO TO 5
C ENTRY ZZZY54 ====2====*====3====*====4====*====5====*====6====*====7
C    PREPERATION OF THE COMBINATIONS
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY ZZZY54
      IC(1,1)=1
      IC(1,2)=1
      IC(2,2)=1
      DO 1 I=3,16
      IC(1,I)=1
      IC(I,I)=1
      DO 2 J=1,I-2
    2 IC(J+1,I)=IC(J,I-1)+IC(J+1,I-1)
    1 CONTINUE
      RETURN
    5 NX=M+N+1
      DO 3 I=1,NX
      KI=MAX0(1,I-N)
      KA=MIN0(M+1,I)
      IW=0
      DO 4 K=KI,KA
    4 IW=IW+IC(K,M+1)*IC(I-K+1,N+1)
     &     *(1-2*MOD(I-K,2))
      IB(I)=IW
    3 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY53 ====*====3====*====4====*====5====*====6====*====7
C
C    K=G.C.M. OF ABS(I) AND ABS(J)
C      IF I.EQ.0  THEN K=ABS(J),   IF J.EQ.0 THEN K=ABS(I)
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY53(I,J,K)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      II=IABS(I)
      JJ=IABS(J)
      IF(II.EQ.0) GO TO 1
      IF(JJ.EQ.0) GO TO 4
    3 KK=MOD(II,JJ)
      IF(KK.EQ.0) GO TO 2
      II=JJ
      JJ=KK
      GO TO 3
    2 K=JJ
      RETURN
    1 K=JJ
      RETURN
    4 K=II
      RETURN
      END
C SUBROUTINE TSIREP ====*====3====*====4====*====5====*====6====*====7
C
C     MAIN ROUTINE TO CALCULTE THE IRREDUCIBLEREPRESENTATION.
C     THE METHOD IS GIVEN BY V.C.SAHNI AND G.VENKATARAMAN
C     PHYS. KONDENS. MATERIE 11,199-211(1970)
C     DIMENSION OF THE REPRESENTATION IS IN ND(12)
C     MATRIX OF THE IRREDUCIBLE REPRESENTATION IS STORED IN
C     IR(7,48,12)=IR(ROW,GROUP ELEMENT,NUMBER OF REPRESENTATION).
C     AT MOST 6 ELEMENTS IN ONE ROW IS STORED IN ONE WORD, BY USING
C     5 BITS FOR EACH.  ALL ELEMNETS CAN BE EXPRESSED BY 0 OR
C     EXP(2PAI*I*(N/24)) --N=1,2,3,,,,24--.
C     THE VALUES STORED IN THESE 5 BITS ARE THESE N OR 0.
C     WHEN YOU FIND A EXPRESSION WITH NUMBER 32, DECOMPOSITION TO
C     EACH ELEMENT IS DONE THERE
C     ROW=7 IS NOT FOR THE ELEMENTS BUT KEEP THE NORMARIZATION OF THE
C     MATRIX.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSIREP(KKBB,IICC,IDDOUB)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,C
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/SPG5/IH(48),KH(24,2),IO,NN,IPA(48,2)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/SPG5/
      INTEGER IB(6),LH(2),IO1(12),ION(12,6)
     &         ,IND(12),JH(2),M(7),M2(7),M3(7)
     &         ,JR(7,48,12),JX(7),JY(7),JZ(7),MD(12)
      DIMENSION KKBB(3)
      CALL CHKDNM(IICC)
C---------------------------------------------
C     FACTORSYSTEM
      CALL ZZZY34(KKBB,IICC,IDDOUB)
      CALL ZZZY26
C---------------------------------------------
C     CALL TSPKDS
C---------------------------------------------
C   GROUP OF K HAS ONLY THE IDENTITY ELEMENT
      IF(MG.GT.1) GO TO 85
      NR=1
      NH=1
      ND(1)=1
      IR(1,1,1)=24*(32**5)
      IR(7,1,1)=1
      CR(1,1)=1.0
      CALL ZZZY35
      RETURN
C---------------------------------------------
C     STARTING CYCLIC GROUP
   85 CONTINUE
      ITA=IG(JG(2))
      NH=2
      IF(IL.LE.0) GO TO 83
      IF(ITA.GE.5.AND.ITA.LE.12) NH=3
      ITB=0
      IF(MG.GE.4) ITB=IG(JG(3))
      IF(ITB.GE.19.AND.ITB.LE.24) NH=4
      IF(ITB.GE.43.AND.ITB.LE.48) NH=4
      GO TO 84
   83 CONTINUE
      IF(ITA.EQ.2) NH=6
      IF(ITA.EQ.3) NH=3
   84 CONTINUE
      NR=NH
      IF(NH.LT.6) GO TO 81
      JX(1)=24
      IP=24
      DO 82 J=2,6
      JX(J)=IP
      CALL ZZZY29(IP,IFA(2,J)*(24/IFC),1)
   82 CONTINUE
      GO TO 6
   81 CONTINUE
      JX(1)=24
      JX(2)=24
      IF(NH.LE.3) GO TO 1
      JX(3)=IFA(3,3)*(24/IFC)
      JX(4)=IFA(3,2)*(24/IFC)
      CALL ZZZY29(JX(4),JX(3),1)
      IP=JX(4)
      CALL ZZZY29(IP,IFA(3,4)*(24/IFC),1)
      GO TO 6
    1 CONTINUE
      IP=IFA(2,2)*(24/IFC)
      IQ=IFA(2,3)*(24/IFC)
      JX(3)=IP
      IF(NH.EQ.3) CALL ZZZY29(IP,IQ,1)
    6 CALL ZZZY28(IP,NH,IB)
      DO 2 I=1,NR
      ND(I)=1
      IA=24
      DO 3 J=1,NH
      JJ=J
      IF(NH.EQ.4.AND.J.EQ.2) JJ=3
      IF(NH.EQ.4.AND.J.EQ.3) JJ=2
      JA=IA
      CALL ZZZY29(JA,JX(J),-1)
      IR(1,JJ,I)=JA*32**5
      IR(7,JJ,I)=1
      DO 4 K=2,6
    4 IR(K,JJ,I)=0
      CALL ZZZY29(IA,IB(I),1)
    3 CONTINUE
      IH(I)=IG(JG(I))
    2 CONTINUE
C----------------------------------------------------------
C  CALCULATION OF THE CHARACTER
      CALL ZZZY32
C----------------------------------------------------------
    5 CONTINUE
C     CALL TSIRDS
      IF(NH.GE.MG) GO TO 75
C----------------------------------------------------------
C     ORBIT CLASSIFICATION
C
      JH(1)=IG(JG(NH+1))
      JH(2)=IM(JH(1),JH(1))
      N=2
      IF(JH(2).NE.1) N=3
      IF(N.EQ.2) GO TO 91
      IF(N.EQ.3.AND.IM(JH(2),JH(1)).EQ.1) GO TO 91
      JH(1)=IG(JG(NH+2))
      LH(1)=NH+2
      N=2
      GO TO 8
   91 CONTINUE
      LH(1)=NH+1
      IF(N.EQ.2) GO TO 8
      DO 10 I=NH+2,MG
      IF(IG(JG(I)).EQ.JH(2)) GO TO 9
   10 CONTINUE
    9 LH(2)=I
    8 CONTINUE
C     WRITE(6,651)
C 651 FORMAT(//' CONJUGATE ELEMENT')
      INDI=0
      DO 11 I=1,N-1
      JHV=IV(JH(I))
      DO 12 J=1,NH
      IHA=IM(JHV,IM(IH(J),JH(I)))
      DO 13 K=1,NH
      IF(IHA.EQ.IH(K)) GO TO 14
   13 CONTINUE
   14 KH(J,I)=K
      IF(K.NE.J) INDI=1
   12 CONTINUE
C     WRITE(6,602) (KH(J,I),J=1,NH)
C 602 FORMAT(24I3)
   11 CONTINUE
      DO 19 I=1,N-1
      IW=LH(I)
      DO 15 J=1,NH
      IPA(J,I)=IFA(J,IW)*(24/IFC)
      IP=IFA(IW,KH(J,I))*(24/IFC)
      CALL ZZZY29(IPA(J,I),IP,-1)
   15 CONTINUE
   19 CONTINUE
      DO 16 I=1,NR
   16 IND(I)=0
      I1=0
      IN=0
      I=1
   27 DO 17 J=1,NH
      C=CR(KH(J,1),I)
      CALL ZZZY27(C,IPA(J,1))
      IF(CDABS(C-CR(J,I)).GT.1.0D-4) GO TO 18
   17 CONTINUE
      I1=I1+1
      IO1(I1)=I
      IND(I)=1
      GO TO 23
   18 IN=IN+1
      ION(1,IN)=I
      IND(I)=1
      DO 21 JN=1,N-1
      DO 20 IW=I+1,NR
      IF(IND(IW).EQ.1) GO TO 20
      DO 22 J=1,NH
      C=CR(KH(J,JN),I)
      CALL ZZZY27(C,IPA(J,JN))
      IF(CDABS(C-CR(J,IW)).GT.1.0D-4) GO TO 20
   22 CONTINUE
      ION(JN+1,IN)=IW
      IND(IW)=1
      GO TO 21
   20 CONTINUE
   21 CONTINUE
   23 IF(I.GE.NR) GO TO 26
      DO 24 IW=I+1,NR
      IF(IND(IW).EQ.0) GO TO 25
   24 CONTINUE
      GO TO 26
   25 I=IW
      GO TO 27
   26 CONTINUE
C     WRITE(6,611) I1
C 611 FORMAT(' NUMBER OF ORBIT OF ORDER ONE=',I3)
C     IF(I1.GT.0) WRITE(6,603) (IO1(I),I=1,I1)
C     WRITE(6,613) N,IN
C 613 FORMAT(' NUMBER OF ORBIT OF ORDER',I3,2H =,I3)
C     IF(IN.GT.0) WRITE(6,603)((ION(I,J),I=1,N),J=1,IN)
C 603 FORMAT(20I5)
      JJ=0
C
C     ORBIT OF ORDER ONE
C
      IF(I1.EQ.0) GO TO 61
      DO 31 IW=1,I1
      IO=IO1(IW)
      NN=ND(IO)
      IF(NN.EQ.1) GO TO 42
      IF(INDI.EQ.0.AND.NN.EQ.3) GO TO 42
      CALL ZZZY33(M)
      CALL ZZZY31(M,M,M2,NN)
      IF(N.EQ.3) CALL ZZZY31(M2,M,M3,NN)
C     WRITE(6,604) M2
C     IF(N.EQ.3) WRITE(6,604) M3
C 604 FORMAT(6(3X,O12))
      LAM=M2(1)/(32**5)
      IF(N.EQ.3) LAM=M3(1)/(32**5)
      GO TO 41
   42 LAM=24
      DO 37 I=1,6
   37 M(I)=0
      DO 38 I=1,NN
      M(I)=LAM*(32**(6-I))
   38 CONTINUE
      M(7)=1
      CALL ZZZY31(M,M,M2,NN)
   41 IP=IFA(LH(1),LH(1))*(24/IFC)
      IF(N.EQ.3) CALL ZZZY29(IP,IFA(NH+1,LH(2))*(24/IFC),1)
      IX=IP
      CALL ZZZY29(IX,LAM,-1)
      CALL ZZZY28(IX,N,IB)
      DO 43 J=1,N
      JJ=JJ+1
      MD(JJ)=NN
      DO 44 I=1,NH
      DO 45 K=1,7
      JR(K,I,JJ)=IR(K,I,IO)
   45 CONTINUE
   44 CONTINUE
      DO 46 MA=1,N-1
      IF(MA.EQ.2) JW=IFA(NH+1,NH+1)*(24/IFC)
      DO 47 K=1,7
      JX(K)=M(K)
      IF(MA.EQ.2) JX(K)=M2(K)
   47 CONTINUE
      CALL ZZZY30(JX,IB(J),NN)
      IF(MA.EQ.2) CALL ZZZY30(JX,IB(J),NN)
      IF(MA.EQ.2) CALL ZZZY30(JX,JW,NN)
      DO 48 MB=1,NH
      JA=IM(JH(MA),IG(JG(MB)))
      DO 49 MC=NH+1,MG
      IF(JA.EQ.IG(JG(MC))) GO TO 50
   49 CONTINUE
   50 DO 51 K=1,7
      JY(K)=IR(K,MB,IO)
   51 CONTINUE
      CALL ZZZY30(JY,(IFC-IFA(LH(MA),MB))*(24/IFC),NN)
      CALL ZZZY31(JX,JY,JZ,NN)
      DO 52 K=1,7
      JR(K,MC,JJ)=JZ(K)
   52 CONTINUE
   48 CONTINUE
   46 CONTINUE
   43 CONTINUE
   31 CONTINUE
C
C     ORBIT OF ORDER N
C
   61 CONTINUE
      IF(IN.EQ.0) GO TO 71
      DO 62 IW=1,IN
      IO=ION(1,IW)
      NN=ND(IO)
      JJ=JJ+1
      MD(JJ)=NN*N
      DO 63 K=1,NH*N
      DO 69 K1=1,6
   69 JR(K1,K,JJ)=0
      JA=IG(JG(K))
      DO 64 I=1,N
      JB=JA
      IF(I.NE.1) JB=IM(IV(JH(I-1)),JA)
      DO 65 J=1,N
      JC=JB
      IF(J.NE.1) JC=IM(JB,JH(J-1))
      DO 66 KA=1,NH
      IF(JC.EQ.IG(JG(KA))) GO TO 67
   66 CONTINUE
   65 CONTINUE
   67 JR(7,K,JJ)=IR(7,KA,IO)
      DO 68 J1=1,NN
      DO 60 J2=1,NN
      IW1=32**(6-J2)
      MA=IR(J1,KA,IO)/IW1
      MA=MA-(MA/32)*32
      IF(MA.EQ.0) GO TO 60
      IF(J.NE.1) CALL ZZZY29(MA,IFA(K,LH(J-1))*(24/IFC),1)
      IF(I.NE.1) CALL ZZZY29(MA,IFA(LH(I-1),KA)*(24/IFC),-1)
      MB=(I-1)*NN+J1
      MC=(J-1)*NN+J2
      JR(MB,K,JJ)=JR(MB,K,JJ)+MA*(32**(6-MC))
   60 CONTINUE
   68 CONTINUE
   64 CONTINUE
   63 CONTINUE
   62 CONTINUE
   71 NR=JJ
      NH=NH*N
      DO 72 I=1,NR
      ND(I)=MD(I)
      DO 73 J=1,NH
      IH(J)=IG(JG(J))
      DO 74 K=1,7
      IR(K,J,I)=JR(K,J,I)
   74 CONTINUE
   73 CONTINUE
   72 CONTINUE
C-----------------------------------------------------
C    CALCULATION OF CHARACTER
      CALL ZZZY32
C-----------------------------------------------------
      GO TO 5
   75 CONTINUE
C-----------------------------------------------------
C   TIME REVERSAL ELEMENTS AND HERRING'S SUM
      CALL ZZZY35
C-----------------------------------------------------
      RETURN
      END
C SUBROUTINE TSIRDS ====*====3====*====4====*====5====*====6====*====7
C
C    PRINT OUT THE MATRIX ELEMENTS OF THE IRREDUCIBLE REPRESENTATION
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSIRDS
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &    ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG3/,/SPG4/,CAR
      CHARACTER*1 CAR(25)
      CHARACTER*1 CA(48)
	  DATA CAR/
     &   'A','K','P','X','B','I',
     &   'C','V','Q','L','D','-',
     &   'E','M','R','W','F','J',
     &   'G','Y','S','N','H','+',
     &   '0'/
      IF(IFC.LE.1) GO TO 5
      WRITE(6,601) IFC
  601 FORMAT(//' FACTOR SYSTEM'/' COMMON DENOMINATOR=',I3)
      WRITE(6,602) (J,J=1,MG)
  602 FORMAT(3X,48I2)
      DO 7 I=1,MG
      WRITE(6,600) I,(IFA(I,J),J=1,MG)
    7 CONTINUE
  600 FORMAT(I3,48I2)
    5 WRITE(6,660) NR,NH
  660 FORMAT(//' NUMBER OF IMR=',I7/
     &         ' NUMBER OF ELEMENT=',I3)
      DO 4 I=1,NR
      WRITE(6,661) I,ND(I)
  661 FORMAT(/' IMR NO',I3,'  DIMENSION=',I3)
C    &       ' CHARACTER TABLE')
C     WRITE(6,662) (CR(J,I),J=1,NH)
C 662 FORMAT(6(2X,2F7.3))
      WRITE(6,663)
  663 FORMAT(' MATRIX REPRESENTATION')
      WRITE(6,668) (J,J=1,NH)
  668 FORMAT(4X,48I2)
      NN=ND(I)
      DO 1 I1=1,NN
      DO 2 I2=1,NN
      DO 3 K=1,NH
      IX=IR(I1,K,I)/(32**(6-I2))
      IX=IX-(IX/32)*32
      IF(IX.EQ.0) IX=25
      CA(K)=CAR(IX)
    3 CONTINUE
      WRITE(6,669) I1,I2,(CA(K),K=1,NH)
  669 FORMAT(2I2,48A2)
    2 CONTINUE
    1 CONTINUE
    4 CONTINUE
      RETURN
      END
C SUBROUTINE TSIRNR ====*====3====*====4====*====5====*====6====*====7
C
C     NNRR           :NUMBER OF REPRESENTATION
C     NNHH           :ORDER OF QROUP OFPOINT K
C     NNDD(12)       :DIMENSION OF REPRESENTATION
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSIRNR(NNRR,NNHH,NNDD)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG4/
      DIMENSION NNDD(12)
      DO 1 I=1,NR
    1 NNDD(I)=ND(I)
      NNHH=NH
      NNRR=NR
      RETURN
      END
C SUBROUTINE ZZZY27 ====*====3====*====4====*====5====*====6====*====7
C     SLAVE ROUTINE FOR TSIREP
C       C*EXP(2PAI*I*(IP/24)) --------> C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY27(C,IP)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 C
      T=2.0D0*3.14159265358979312D0
     &   *(DBLE(MOD(IP,32))/DBLE(24))
      C=C*DCMPLX(COS(T),SIN(T))
      RETURN
      END
C SUBROUTINE ZZZY28 ====*====3====*====4====*====5====*====6====*====7
C    SLAVE ROUTINE OF TSIREP
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
C
      SUBROUTINE ZZZY28(IP,N,IB)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IB(6)
      JA=MOD(IP,32)
      IF(JA.EQ.0) JA=24
      IA=JA/N
      IF(JA.NE.IA*N) WRITE(6,600) IA,JA,N
  600 FORMAT(' ERROR IN CP',3I5)
      NI=24/N
      IB(1)=IA+(N-1)*NI
      DO 1 I=2,N
      IB(I)=IA+(I-2)*NI
    1 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY29 ====*====3====*====4====*====5====*====6====*====7
C    SLAVE ROUTINE OF TSIREP
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY29(IA,IB,IS)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      JW1=MOD(IA,32)
      JW2=MOD(IB,32)
      IW=JW1+JW2*IS
      IF(IW.LT.0) IW=IW+24
      IW=MOD(IW,24)
      IF(IW.EQ.0) IW=24
      IA=IA-JW1+IW
      RETURN
      END
C SUBROUTINE ZZZY30 ====*====3====*====4====*====5====*====6====*====7
C   SLAVE ROUTINE OF TSIREP
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY30(IX,JA,N)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IX(7),IY(7)
      J2=MOD(JA,32)
      DO 1 I=1,6
    1 IY(I)=0
      IY(7)=IX(7)
      DO 2 I=1,N
      DO 3 J=1,N
      J3=IX(I)/(32**(6-J))
      J3=J3-(J3/32)*32
      IF(J3.EQ.0) GO TO 3
      IW=J3+J2
      IW=MOD(IW,24)
      IF(IW.EQ.0) IW=24
      IY(I)=IY(I)+IW*(32**(6-J))
    3 CONTINUE
    2 CONTINUE
      DO 4 I=1,7
    4 IX(I)=IY(I)
      RETURN
      END
C SUBROUTINE ZZZY31 ====*====3====*====4====*====5====*====6====*====7
C   SLAVE ROUTINE OF TSIREP
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY31(IX,IY,IZ,N)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 C
      DIMENSION IX(7),IY(7),IZ(7)
      DO 1 I=1,6
      IZ(I)=0
    1 CONTINUE
      IND=0
      DO 2 I=1,N
      DO 3 J=1,N
      C=(0.D0,0.D0)
      DO 4 K=1,N
      J1=IX(I)/(32**(6-K))
      J1=J1-(J1/32)*32
      IF(J1.EQ.0) GO TO 4
      T1=2.0*3.14159265358979312D0*(DBLE(J1)/DBLE(24))
      J3=IY(K)/(32**(6-J))
      J3=J3-(J3/32)*32
      IF(J3.EQ.0) GO TO 4
      T2=2.0D0*3.14159265358979312D0*(DBLE(J3)/DBLE(24))
      C=C+DCMPLX(COS(T1),SIN(T1))*DCMPLX(COS(T2),SIN(T2))
    4 CONTINUE
      IF(CDABS(C).LT.1.0D-4) GO TO 3
      W=ATAN2(DIMAG(C),DBLE(C))
      W=W/3.14159265358979312D0
      IF(W.LT.0.D0) W=2.0D0+W
      IW=W*12.0+0.5
      IF(IW.EQ.0) IW=24
      IZ(I)=IZ(I)+IW*(32**(6-J))
      M=CDABS(C)**2+0.5
      IF(IND.EQ.1) GO TO 5
      IZ(7)=(IX(7)*IY(7))/M
      IND=1
      GO TO 3
    5 IF(IZ(7)*M.NE.IX(7)*IY(7)) GO TO 9
    3 CONTINUE
    2 CONTINUE
      RETURN
    9 WRITE(6,600)
  600 FORMAT(' ERROR IN CMM')
      RETURN
      END
C SUBROUTINE ZZZY32 ====*====3====*====4====*====5====*====6====*====7
C  CHARACTERS OF I.R. ARE CALCULATED
C--*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY32
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 C,CR
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG4/
      DO 2 I=1,NR
      N=ND(I)
      DO 3 J=1,NH
      C=(0.D0,0.D0)
      DO 1 K=1,N
      IW=32**(6-K)
      JA=IR(K,J,I)/IW
      JA=JA-(JA/32)*32
      IF(JA.EQ.0) GO TO 1
      IF(JA.EQ.24) JA=0
      T=2.0D0*3.14159265358979312D0*(DBLE(JA)/DBLE(24))
      C=C+DCMPLX(COS(T),SIN(T))
    1 CONTINUE
      CR(J,I)=C
     & /DCMPLX(SQRT(DBLE(IR(7,J,I))),0.D0)
    3 CONTINUE
    2 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY33 ====*====3====*====4====*====5====*====6====*====7
C  SLAVE ROUTINE OF TSIEP.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY33(M)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/SPG5/IH(48),KH(24,2),IO,NN,IPA(48,2)
      SAVE /SPG4/,/SPG5/
      DIMENSION IN(4,4),IE(6),M(7)
      NE=0
      DO 1 J1=1,NN
      DO 1 J2=1,NN
      IN(J1,J2)=63
      IF(NN.EQ.4.AND.J1.NE.J2) IN(J1,J2)=0
    1 CONTINUE
      DO 2 I=2,NH
      K=KH(I,1)
      IF(IR(7,I,IO).NE.1) GO TO 2
      IF(IR(7,K,IO).NE.1) GO TO 2
      DO 3 J1=1,NN
      DO 4 J4=1,NN
      JA=IR(J1,I,IO)/(32**(6-J4))
      JA=JA-(JA/32)**32
      IF(JA.NE.0) GO TO 5
    4 CONTINUE
    5 DO 6 J2=1,NN
      JC=JA
      IJ1=J4*16+J2
      DO 7 J3=1,NN
      JB=IR(J3,K,IO)/(32**(6-J2))
      JB=JB-(JB/32)*32
      IF(JB.NE.0) GO TO 8
    7 CONTINUE
    8 IF(IN(J1,J3).EQ.0) GO TO 16
      IF(IN(J4,J2).EQ.0) GO TO 15
      IJ2=J1*16+J3
      CALL ZZZY29(JB,IPA(I,1),1)
      IF(IJ1-IJ2) 9,10,11
    9 IF(J4.NE.1) GO TO 6
      CALL ZZZY29(JC,JB,-1)
C
C                    65536 = 2**16
C
      NEW=IJ1*65536+IJ2*256+JC
      GO TO 12
   11 IF(J1.NE.1) GO TO 6
      CALL ZZZY29(JB,JC,-1)
      NEW=IJ2*65536+IJ1*256+JB
   12 CONTINUE
C     WRITE(6,612) I,J1,J2,NEW
C 612 FORMAT(3I5,2X,Z8)
      IF(NE.EQ.0) GO TO 13
      DO 14 J=1,NE
      IF(IE(J).EQ.NEW) GO TO 6
   14 CONTINUE
   13 NE=NE+1
C     WRITE(6,612) NE,I,J1,J2,NEW
C 612 FORMAT(4I5,2X,Z8)
      IE(NE)=NEW
      GO TO 6
   10 IF(JB.EQ.JA) GO TO 6
   15 IN(J1,J3)=0
      GO TO 6
   16 IN(J4,J2)=0
    6 CONTINUE
    3 CONTINUE
    2 CONTINUE
C     WRITE(6,604) (IE(I),I=1,NE)
C 604 FORMAT(4(3X,Z8))
C     WRITE(6,605) ((IN(I,J),J=1,NN),I=1,NN)
C 605 FORMAT(16I4)
      DO 21 J2=1,NN
      IF(IN(1,J2).NE.0) GO TO 22
   21 CONTINUE
   22 JC=24
      IN(1,J2)=JC
      IJ=16+J2
      DO 23 I=1,NE
      IF(IJ.NE.IE(I)/65536) GO TO 23
      J3=IE(I)/4096
      J3=J3-(J3/16)*16
      J4=IE(I)/256
      J4=J4-(J4/16)*16
      IF(IN(J3,J4).NE.63) GO TO 23
      JA=IE(I)-(IE(I)/256)*256
      CALL ZZZY29(JA,JC,1)
      IN(J3,J4)=JA
   23 CONTINUE
C     WRITE(6,606) ((IN(I,J),J=1,NN),I=1,NN)
C 606 FORMAT(9(2X,I2))
      DO 27 I=1,6
   27 M(I)=0
      DO 28 J1=1,NN
      DO 26 J2=1,NN
      M(J1)=M(J1)+IN(J1,J2)*(32**(6-J2))
   26 CONTINUE
   28 CONTINUE
      M(7)=1
      IF(NN.EQ.2.AND.IN(1,1).NE.0.AND.IN(1,2).NE.0) M(7)=2
      RETURN
      END
C SUBROUTINE TSCRST ====*====3====*====4====*====5====*====6====*====7
C
C   CRYSTAL STRUCTURE STORED IN /ATT/
C     ATOMIC POSITIONS TRASFORMAED BY THE SPACE GROUP OPERATIONS
C     ARE TABLED IN ISITR(50,48)
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSCRST(VA,KKAA,NKAT,NAT,KJ)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/ATT   /ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &   ,NKATOM,NATOM,KS(11,LMNKAT),IAA(LMNATM),JRCH
      SAVE /SPG1/,/SPG2/,/ATT   /
      DIMENSION VA(3,LMNATM),KKAA(2,LMNKAT),KJ(11,LMNKAT),VD(3)
      IF(IL.LE.0) JRCH=2
      IF(IL.GE.1) JRCH=1
      NKATOM=NKAT
      NATOM=NAT
      DO 31 I=1,NKATOM
      DO 33 K=KKAA(1,I),KKAA(2,I)
   33 IAA(K)=I
   31 KION(I)=KKAA(2,I)
      DO 32 I=1,NATOM
      DO 32 K=1,3
   32 VATOM(K,I)=VA(K,I)
      DO 34 K=1,NKATOM
      DO 34 J=1,11
   34 KS(J,K)=KJ(J,K)
      DO 6 J=1,NATOM
      DO 7 I=1,3
      IF(VATOM(I,J).LT.0.0D0) VATOM(I,J)=VATOM(I,J)+1.0D0
    7 CONTINUE
    6 CONTINUE
      INDIC=0
      INDIC2=0
      WRITE(6,602)
  602 FORMAT(1H ,/' TABLE OF THE TRANSLATION OF ATOMIC POSITION')
C
C
 9000 CONTINUE
C
C---*
      DO 11 J=1,NG
      JA=IG(J)
                          IF(INDIC2.EQ.2) THEN
                          WRITE(6,*) ' J=',J,'    IG(J)=',JA
                          END IF
      DO 12 IA=1,NATOM
      DO 1 K=1,3
      MA=IT(K,JA)
      IF(IABS(MA).EQ.4) THEN
C
C     X-Y
C
      VV=VATOM(1,IA)-VATOM(2,IA)
      IF(MA.EQ.-4) VV=-VV
      IF(VV.LT.0.0D0) VV=1.0D0+VV
      VD(K)=VV
      ELSE
      MB=IABS(MA)
      VD(K)=VATOM(MB,IA)
      IF(MA.LT.0) VD(K)=1.0D0-VD(K)
      END IF
      KE=JV(1,K,J)
      KC=JV(2,K,J)
      VD(K)=VD(K)+DBLE(KE)/DBLE(KC)
      IF(VD(K).GE.0.9999) VD(K)=VD(K)-1.0D0
    1 CONTINUE
                          IF(INDIC2.EQ.2) THEN
                          WRITE(6,9100) IA, (VATOM(K,IA),K=1,3),
     &                                      (VD(K),K=1,3)
 9100                     FORMAT(1H ,5X,'IA=',I2,'(',3F10.5,' ) ---->',
     &                               '(',3F10.5,' )'     )
                          END IF
C
      DO 2 IB=1,NATOM
        IBB=IB
      IF(IL.NE.-1) GO TO 9
        WV1=VATOM(1,IB)-VD(1)
        WV2=VATOM(2,IB)-VD(2)
        WV3=VATOM(3,IB)-VD(3)
        IF(DABS(WV1).LE.1.0D-4.AND.
     &     DABS(WV2).LE.1.0D-4.AND.
     &     DABS(WV3).LE.1.0D-4) GO TO 4
        IF(WV1.LT.0.0D0) WV1=WV1+1.0D0
        IF(WV2.LT.0.0D0) WV2=WV2+1.0D0
        IF(WV3.LT.0.0D0) WV3=WV3+1.0D0
        IF(DABS(WV1-1.0D0/3.0D0).LE.1.0D-4.AND.
     &     DABS(WV2-2.0D0/3.0D0).LE.1.0D-4.AND.
     &     DABS(WV3-2.0D0/3.0D0).LE.1.0D-4) GO TO 4
        IF(DABS(WV1-2.0D0/3.0D0).LE.1.0D-4.AND.
     &     DABS(WV2-1.0D0/3.0D0).LE.1.0D-4.AND.
     &     DABS(WV3-1.0D0/3.0D0).LE.1.0D-4) GO TO 4
      GO TO 2
C
    9 CONTINUE
        NHALF=0
        DO 3 K=1,3
          WV=ABS(VATOM(K,IB)-VD(K))
          IF(WV.LE.1.0D-4) GO TO 3
          IF(ABS(WV-0.5).GT.1.0D-4) GO TO 2
          NHALF=NHALF+1
    3   CONTINUE
        IF(NHALF.EQ.0) GO TO 4
        IF(IL.EQ.2.AND.NHALF.EQ.2) GO TO 4
        IF(IL.EQ.3.AND.NHALF.EQ.3) GO TO 4
        IF(IL.EQ.4.AND.NHALF.EQ.2.AND.
     &    ABS(VATOM(3,IB)-VD(3)).LE.1.0D-4) GO TO 4
    2 CONTINUE
C
      ISITR(IA,JA)=0
      INDIC=1
      GO TO 12
C
    4 ISITR(IA,JA)=IBB
      IF(IAA(IBB).NE.IAA(IA)) INDIC=2
C
   12 CONTINUE
      IF(NATOM.LE.25) WRITE(6,600) J,JA,(ISITR(IA,JA),IA=1,NATOM)
      IF(NATOM.GT.25) WRITE(6,600) J,JA,(ISITR(IA,JA),IA=1,25)
      IF(NATOM.GT.25) WRITE(6,603) (ISITR(IA,JA),IA=26,NATOM)
   11 CONTINUE
C---*
  600 FORMAT(27I3)
  603 FORMAT(6X,27I3)
      IF(INDIC.EQ.0) THEN
        RETURN
      ELSE IF(INDIC2.EQ.0) THEN
        WRITE(6,*) ' === CHECK THE POSITIONS OF ATOMS . ==='
        WRITE(6,*) ' === CHECK THE SYMMETRY OPERATIONS .==='
        INDIC2=2
        GO TO 9000
      ELSE
        STOP ' === STOP IN SUB.TSCRST. (POSITIONS OF ATOMS) ==='
      END IF
C
      END
C SUBROUTINE TSTARK ====*====3====*====4====*====5====*====6====*====7
C   PRINT OUT THE TABLE OF STAR OF K
C   CALCULATIONS ARE DONE IN ZZZY38 AND TABLED IN /STK/
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSTARK
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CW
      CHARACTER*4 A,B,GN(2,3)
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/STK/,CMN,HMN,AC,XYZ
      CHARACTER*4 CMN(24),HMN(12),AC(3),XYZ(4)
      INTEGER KA(2,3)
      DATA CMN/'E','C2X','C2Y','C2Z','C31+','C32+'
     &                   ,'C33+','C34+','C31-','C32-','C33-','C34-'
     &                   ,'C2A','C2B','C2C','C2D','C2E','C2F'
     &                   ,'C4X+','C4Y+','C4Z+','C4X-','C4Y-','C4Z-'/
     &           ,HMN/'E','C6+','C3+','C2','C3-','C6-'
     &                   ,'C211','C221','C231','C212','C222','C232'/
     &           ,AC/' ','I','-'/
     &           ,XYZ/'X','Y','Z','W'/
C---------------------------------------------------------
C  CALCULATIONS ARE DONE IN ZZZY38 AND TABLED IN /STK/
      CALL ZZZY38
C---------------------------------------------------------
      WRITE(6,601)
  601 FORMAT(//' TABLE FOR STAR OF K')
      DO 1 I=1,NS
      ITA=IG(JS(I))
      ITC=MOD(ITA-1,24)+1
      IF(IL.LE.0) ITC=MOD(ITA-1,12)+1
      A=AC(1)
      IF(IL.GT.0.AND.ITA.GT.24) A=AC(2)
      IF(IL.LE.0.AND.ITA.GT.12) A=AC(2)
      IF(IL.GT.0) B=CMN(ITC)
      IF(IL.LE.0) B=HMN(ITC)
      DO 2 J=1,3
      MA=IT(J,ITA)
      IF(IL.LE.0) MA=IT(J,ITA+24)
      MB=IABS(MA)
      IF(MA.GT.0) GN(1,J)=AC(1)
      IF(MA.LT.0) GN(1,J)=AC(3)
      GN(2,J)=XYZ(MB)
      KA(1,J)=JV(1,J,JS(I))
      KA(2,J)=JV(2,J,JS(I))
    2 CONTINUE
      WRITE(6,602) I,(KS(K,I),K=1,3),ICBB,ITA,A,B,GN,KA
  602 FORMAT(' K',I2,'=',3I3,1H/,I3,I5,2X,A1,A4,3(1X,2A1),3(I3,1H/,I1))
    1 CONTINUE
      RETURN
      END
C SUBROUTINE TSPGCR ====*====3====*====4====*====5====*====6====*====7
C
C   PRINT OUT THE TABLE OF CHARACTER OF FULL REPRESENTATION
C   CALCULATIONS ARE DONE IN ZZZY39 AND TABLED IN /STK/----,CW(48,12)
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSPGCR
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CW,CR
      CHARACTER*4 A,B
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG4/,/STK/,CMN,HMN,AC
      CHARACTER*4 CMN(24),HMN(12),AC(3)
      DIMENSION REALP(12),AIMGP(12)
	  DATA CMN/'E','C2X','C2Y','C2Z','C31+','C32+'
     &                   ,'C33+','C34+','C31-','C32-','C33-','C34-'
     &                   ,'C2A','C2B','C2C','C2D','C2E','C2F'
     &                   ,'C4X+','C4Y+','C4Z+','C4X-','C4Y-','C4Z-'/
      DATA HMN/'E','C6+','C3+','C2','C3-','C6-'
     &     ,'C211','C221','C231','C212','C222','C232'/
      DATA AC/' ','I','-'/
C-----------------------------------------------------------
C   CALCULATIONS ARE DONE IN ZZZY39
      CALL ZZZY39
C-----------------------------------------------------------
      WRITE(6,601)
  601 FORMAT(//' CHARACTER TABLE OF THE FULL SPACE GROUP IR.')
      WRITE(6,602) (I,I=1,NR)
  602 FORMAT(5X,12I10)
      DO 1 I=1,NG
      ITA=IG(I)
      ITC=MOD(ITA-1,24)+1
      IF(IL.LE.0) ITC=MOD(ITA-1,12)+1
      A=AC(1)
      IF(IL.GT.0.AND.ITA.GT.24) A=AC(2)
      IF(IL.LE.0.AND.ITA.GT.12) A=AC(2)
      IF(IL.GT.0) B=CMN(ITC)
      IF(IL.LE.0) B=HMN(ITC)
      DO 2 J=1,NR
      REALP(J)=DBLE(CW(I,J))
    2 AIMGP(J)=DIMAG(CW(I,J))
      WRITE(6,603) I,A,B,(REALP(J),J=1,NR)
      WRITE(6,604) (AIMGP(J),J=1,NR)
  603 FORMAT(I3,2X,A1,A4,12F10.5)
  604 FORMAT(10X,12F10.5)
    1 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY38 ====*====3====*====4====*====5====*====6====*====7
C
C  CALCULATION FOR STAR OF K
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY38
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CW
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/STK/
      INTEGER KA(3),KC(3),KT(3)
      ICBB=ICB
      NS=1
      DO 1 I=1,3
    1 KS(I,NS)=KB(I)
      JS(NS)=1
      DO 2 I=1,NG
      IA=IG(I)
      IF(IL.LE.0) IA=IA+24
      DO 3 J=1,3
      MA=IT(J,IA)
      IF(IABS(MA).EQ.4) KW=KB(1)+KB(2)
      IF(IABS(MA).NE.4) KW=KB(IABS(MA))
      IF(MA.GT.0) KT(J)=KW
      IF(MA.LT.0) KT(J)=-KW
    3 CONTINUE
      DO 4 J=1,NS
      INDC=0
      DO 5 K=1,3
      KA(K)=KS(K,J)-KT(K)
      IF(MOD(IABS(KA(K)),ICB).NE.0) GO TO 4
      KA(K)=KA(K)/ICB
      KC(K)=IABS(KA(K))
      IF(KA(K).NE.0) INDC=1
    5 CONTINUE
      IF(INDC.EQ.0) GO TO 13
      IF(IL.EQ.-1) GO TO 7
      IF(IL.EQ.0) GO TO 13
      IF(IL.EQ.1) GO TO 13
      IF(IL.EQ.2) GO TO 8
      IF(IL.EQ.3) GO TO 9
      IF(IL.EQ.4) GO TO 10
    7 IIWW=-KA(1)+KA(2)+KA(3)
      IF(MOD(IABS(IIWW),3).NE.0) GO TO 4
      GO TO 13
    8 IF(MOD(KC(1),2).NE.MOD(KC(2),2)) GO TO 4
      IF(MOD(KC(2),2).NE.MOD(KC(3),2)) GO TO 4
      GO TO 13
    9 IF(MOD(KC(1)+KC(2)+KC(3),2).NE.0) GO TO 4
      GO TO 13
   10 IF(MOD(KC(1)+KC(2),2).NE.0) GO TO 4
      GO TO 13
    4 CONTINUE
      NS=NS+1
      DO 11 K=1,3
   11 KS(K,NS)=KT(K)
      JS(NS)=I
      GO TO 2
   13 DO 14 K=1,3
      IF(JV(1,K,JS(J)).LT.JV(1,K,I)) GO TO 2
      IF(JV(1,K,JS(J)).GT.JV(1,K,I)) GO TO 19
   14 CONTINUE
      GO TO 2
   19 IF(J.EQ.NS) GO TO 15
      DO 16 JJ=J,NS-1
      JS(JJ)=JS(JJ+1)
      DO 17 K=1,3
   17 KS(K,JJ)=KS(K,JJ+1)
   16 CONTINUE
   15 JS(NS)=I
      DO 18 K=1,3
   18 KS(K,NS)=KT(K)
    2 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY39 ====*====3====*====4====*====5====*====6====*====7
C
C CALCULATION FOR THE CHARACTER OF FULL REPRESENTATION OF SPACE GROUP
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY39
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW,WW
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/STK/
      CALL ZZZY38
      DO 3 I=1,NR
      DO 3 J=1,NG
    3 CW(J,I)=(0.D0,0.D0)
      DO 4 I=1,NG
      DO 5 J=1,NS
      KW=IM(IV(IG(JS(J))),IM(IG(I),IG(JS(J))))
      DO 6 K=1,MG
      IF(KW.EQ.IG(JG(K))) GO TO 7
    6 CONTINUE
      GO TO 5
    7 IA=IG(JS(J))
      IF(IL.LE.0) IA=IA+24
      IW=0
      JW=1
      DO 8 L=1,3
      IF(JV(1,L,JS(J)).EQ.0) GO TO 9
      MA=IT(L,IA)
      IF(IABS(MA).EQ.4) KW=JK(1,K)+JK(2,K)
      IF(IABS(MA).NE.4) KW=JK(IABS(MA),K)
      IF(MA.LT.0) KW=-KW
      IW=IW*JV(2,L,JS(J))+JW*KW*JV(1,L,JS(J))
      JW=JW*JV(2,L,JS(J))
    9 IW=IW*JV(2,L,I)*ICB-JW*JV(1,L,I)*KS(L,J)
      JW=JW*JV(2,L,I)*ICB
    8 CONTINUE
      W=2.0D0*3.14159265358979312D0*(DBLE(IW)/DBLE(JW))
      WW=DCMPLX(COS(W),SIN(W))
      DO 10 L=1,NR
      CW(I,L)=CW(I,L)+CR(K,L)*WW
   10 CONTINUE
    5 CONTINUE
C
      DO 1 L=1,NR
      IF(ABS(DBLE(CW(I,L))).LT.1.0D-4)
     &              CW(I,L)=DCMPLX(0.D0,DIMAG(CW(I,L)))
    1 IF(ABS(DIMAG(CW(I,L))).LT.1.0D-4)
     &              CW(I,L)=DCMPLX(DBLE(CW(I,L)),0.D0)
C
    4 CONTINUE
C
      RETURN
      END
C SUBROUTINE ZZZY40 ====*====3====*====4====*====5====*====6====*====7
C
C   WAWE VECTOR GIVEN KM(4,JU) IS TRANSFORMED BY OPERATION WITH CODE JA
C   TO KK(3).
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY40(JA,JU,KK)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 U
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPW/,/SPG2/
      DIMENSION IT(3,48),IH(3,24),KK(3)
      DATA IT/
     &  4,5,6, 4,1,0, 2,5,0, 2,1,6, 6,4,5, 0,4,1,
     &  0,2,5, 6,2,1, 5,6,4, 5,0,2, 1,6,2, 1,0,4,
     &  5,4,0, 1,2,0, 6,1,4, 2,6,5, 0,1,2, 2,0,1,
     &  4,0,5, 6,5,2, 1,4,6, 4,6,1, 0,5,4, 5,2,6,
     &  2,1,0, 2,5,6, 4,1,6, 4,5,0, 0,2,1, 6,2,5,
     &  6,4,1, 0,4,5, 1,0,2, 1,6,4, 5,0,4, 5,6,2,
     &  1,2,6, 5,4,6, 0,5,2, 4,0,1, 6,5,4, 4,6,5,
     &  2,6,1, 0,1,4, 5,2,0, 2,0,5, 6,1,2, 1,4,0/
      DATA IH/
     &  4,5,6, 1,3,6, 7,4,6, 2,1,6, 5,7,6, 3,2,6,
     &  2,3,0, 3,1,0, 1,2,0, 4,7,0, 7,5,0, 5,4,0,
     &  2,1,0, 5,7,0, 3,2,0, 4,5,0, 1,3,0, 7,4,0,
     &  4,7,6, 7,5,6, 5,4,6, 2,3,6, 3,1,6, 1,2,6/
      DO 1 K=1,3
      MA=IT(K,JA)-3
      IF(IL.LE.0) MA=IH(K,JA)-3
      IF(MA.EQ.0) GO TO 2
      IF(MA.EQ.4) GO TO 3
      MB=IABS(MA)
      KK(K)=KM(MB,JU)
      IF(MA.LT.0) KK(K)=-KK(K)
      GO TO 1
    2 KK(K)=KM(1,JU)+KM(2,JU)
      GO TO 1
    3 KK(K)=-(KM(1,JU)+KM(2,JU))
    1 CONTINUE
      RETURN
      END
C SUBROUTINE ZZZY41 ====*====3====*====4====*====5====*====6====*====7
C   WAVE VECTORS OF KB(3)/IC+G ARE GENERATED IN KM(,)
C   A() HAS A LENGTH OF WAVE VECTOR
C   KT() HAS A NUMBER WHICH INDICATES CHANGE POINT OF THE LENGTH
C   IN THE TABLE OF KM(,)
C   IK :NUMBER OF WAVE VECTOR IN KM(,)
C   IE :NUMBER OF SHELL
C
C   BUG FIX by A.YANASE     1998/01/26
C           
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY41(KKB,ICC,AM)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 U
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      COMMON/LAT/AA,B,C,CA,CB,CC,A1,B1,C1
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPW/,/LAT/,/SPG2/
      DIMENSION KB(3),KKB(3),KBB(3,10),IRECP(3,10)
      CALL CHKDNM(ICC)
      CALL NEAREC(KKB,ICC,KBB,IRECP,NNG)
      KB(1)=KBB(1,1)
      KB(2)=KBB(2,1)
      KB(3)=KBB(3,1)
      IC=ICC
      KXM=AM+4.5
      KYM=AM*(A1/B1)*(B/AA)+4.5
      KZM=AM*(A1/C1)*(C/AA)+4.5
      IK=1
      KM(1,1)=KB(1)
      KM(2,1)=KB(2)
      KM(3,1)=KB(3)
      CALL ZZZY37(KB(1),KB(2),KB(3),IC,S)
      A(1)=S
      DO 41 IX=1,KXM
      KX=(IX-1)*IC
      DO 42 I1=1,2
      DO 43 IY=1,KYM
      KY=(IY-1)*IC
      DO 44 I2=1,2
      DO 45 IZ=1,KZM
      IF(IX*IY*IZ.EQ.1) GO TO 45
      KZ=(IZ-1)*IC
      DO 46 I3=1,2
      IF(IL.EQ.0) GO TO 460
      IF(IL.EQ.1) GO TO 460
      IF(IL.EQ.2) GO TO 462
      IF(IL.EQ.3) GO TO 463
      IF(IL.EQ.-1) GO TO 461
      IF(MOD(IABS(KX+KY),2*IC).NE.0) GO TO 47
      GO TO 460
  461 KSUM=-KX+KY+KZ
      IF(MOD(IABS(KSUM),3*IC).NE.0) GO TO 47
      GO TO 460
  463 IF(MOD(IABS(KX+KY+KZ),2*IC).NE.0) GO TO 47
      GO TO 460
  462 IF(MOD(IABS(KX),2*IC).NE.MOD(IABS(KY),2*IC)) GO TO 47
      IF(MOD(IABS(KY),2*IC).NE.MOD(IABS(KZ),2*IC)) GO TO 47
  460 CONTINUE
      K1=KX+KB(1)
      K2=KY+KB(2)
      K3=KZ+KB(3)
      CALL ZZZY37(K1,K2,K3,IC,W)
      IF(W.GT.AM) GO TO 47
      DO 48 JK=1,IK
      JJK=IK-JK+1
      IF(A(JJK).LT.W) GO TO 49
      DO 50 I=1,3
   50 KM(I,JJK+1)=KM(I,JJK)
      A(JJK+1)=A(JJK)
   48 CONTINUE
      JJK=0
   49 KM(1,JJK+1)=K1
      KM(2,JJK+1)=K2
      KM(3,JJK+1)=K3
      A(JJK+1)=W
      IK=IK+1
      IF(IK.GT.MAXNPW) GO TO 53
   47 IF(KZ.EQ.0) GO TO 45
      KZ=-KZ
   46 CONTINUE
   45 CONTINUE
      IF(KY.EQ.0) GO TO 43
      KY=-KY
   44 CONTINUE
   43 CONTINUE
      IF(KX.EQ.0) GO TO 41
      KX=-KX
   42 CONTINUE
   41 CONTINUE
C     DO 51 II=1,IK
C     WRITE(6,601) II, (KM(K,II),K=1,3),IC,A(II)
C 601 FORMAT(I5,2H (,3I4,2H)/,I3,F10.5)
C  51 CONTINUE
      IE=1
      W=A(1)
      KM(4,1)=0
      DO 52 II=2,IK
      KM(4,II)=0
      IF(ABS(A(II)-W).LT.1.0D-4) GO TO 52
      KT(IE)=II-1
      IE=IE+1
      W=A(II)
   52 CONTINUE
      KT(IE)=IK
      RETURN
   53 WRITE(6,660)
  660 FORMAT(' STOP AT 53 IN ZZZY41 AM IS TOO LARGE')
      STOP
      END

C SUBROUTINE TSPWA  ====*====3====*====4====*====5====*====6====*====7
C
C      SYMMETRIZE THE PLANE WAVES.
C      ALL COMPONENTS ARE GENRATED
C      1985.7.20. :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSPWA(KB,IC,IR,NDI,KO,V,B,D,IND,ND1,ND2,NNWW)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 V,U,WD,CR,CW,WW
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      SAVE /SPW/
      DIMENSION KO(1),V(1),B(1),IND(4,1),KB(3),D(1)
      DIMENSION JGA(2,3,48),WD(48),JG(48),CR(48),KK(3),NIP(100)
CC
CC
CC
      PAI2=8.0D0*DATAN(1.0D0)
CC
      CALL CHKDNM(IC)
CC
      DO 13 I=1,IK
   13 KM(4,I)=0
CC
CC
C     WRITE(6,690) IE,IK
C 690 FORMAT(10I5)
      NK=IK
      CALL TSIRMR(IR,1,1,MG,ND,JG,JGA,CR,WD)
C     WRITE(6,690) ND,MG
      NDI=ND
      JV=0
      ND1=0
      ND2=0
      IW=0
      DO 21 ID=1,IE
      NN=ND1+1
      JW=JV+1
      JV=KT(ID)
C
C     CHARACTER TEST
C
      CW=0.0D0
      DO 22 II=1,MG
      IF(CDABS(CR(II)).LT.1.0D-4) GO TO 22
      JA=JG(II)
      WW=0.0D0
      DO 23 JU=JW,JV
      JJ=JU
      CALL ZZZY40(JA,JJ,KK)
      WA=0.0D0
      DO 24 K=1,3
      IF(KK(K).NE.KM(K,JU)) GO TO 23
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
      WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
   24 CONTINUE
      X=PAI2*WA
      WW=WW+DCMPLX(COS(X),SIN(X))
   23 CONTINUE
      CW=CW+CONJG(CR(II))*WW
   22 CONTINUE
      CW=CW/MG
      NI=CW+0.5
C     WRITE(6,690) ID,JW,JV,NI
      IF(NI.EQ.0) GO TO 21
      IP=JW
      NP=0
      NB=1
      NDR=1
C      CALL TSIRMR(IR,1,NB,MG,ND,JG,JGA,CR,WD)
      CALL TSIRME(1,NB,WD)
      INDIC=0
C
C     PROJECTION OPERATOR
C
   47 DO 27 JU=JW,JV
   27 U(JU)=0.0
      DO 28 II=1,MG
      IF(CDABS(WD(II)).LT.1.0D-4) GO TO 28
      JA=JG(II)
      CALL ZZZY40(JA,IP,KK)
      WA=0.0D0
      DO 29 K=1,3
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
   29 WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
      X=PAI2*WA
      WW=DCMPLX(COS(X),SIN(X))
      DO 30 JU=JW,JV
      DO 20 K=1,3
      IF(KM(K,JU).NE.KK(K)) GO TO 30
   20 CONTINUE
      KM(4,JU)=IP
      U(JU)=U(JU)+WW*CONJG(WD(II))
      GO TO 28
   30 CONTINUE
      WRITE(6,602)
  602 FORMAT(' STOP AT 30 IN TSPWA')
      STOP
   28 CONTINUE
C
C     ORTHOGONALITY TEST
C
C     WRITE(6,699) NN,ND1,ND2,NDR,NB,IP
C 699 FORMAT(' ORTHGONAL TEST',6I5)
      IF(NN.GT.ND1) GO TO 31
      DO 32 II=NN,ND1
C     WRITE(6,699) II
      ID1=IND(1,II)
      ID2=IND(2,II)
      WW=0.0D0
      DO 33 I=ID1,ID2
      JU=JW
   36 IF(KO(I).EQ.JU) GO TO 34
      GO TO 35
   34 CONTINUE
      WW=WW+CONJG(V(I))*U(JU)
C     WRITE(6,697) WW,I,JU
C 697 FORMAT(' ORTH WW',2F10.6,2I5)
      GO TO 33
   35 JU=JU+1
      IF(JU.LE.JV) GO TO 36
      WRITE(6,603)
  603 FORMAT(' STOP AT 36 IN TSPWA')
      STOP
   33 CONTINUE
      IF(CDABS(WW).LE.1.0D-4) GO TO 32
      IF(NDR.EQ.1) GO TO 37
      WRITE(6,606)
  606 FORMAT(' STOP AT 32 IN TSPWA')
      STOP
   32 CONTINUE
C
C     REGISTRATION
C
   31 ND1=ND1+1
      IND(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW+1
      DO 38 JU=JW,JV
      IF(CDABS(U(JU)).GT.1.0D-4) GO TO 381
      IF(JU.EQ.IP) GO TO 381
      IF(KM(4,JU).EQ.IP) KM(4,JU)=0
      GO TO 38
  381 IW=IW+1
      KO(IW)=JU
      V(IW)=U(JU)*(DBLE(ND)/DBLE(MG))
      IF(ABS(DIMAG(V(IW))).LT.1.D-4) V(IW)=DBLE(V(IW))
      IF(ABS(DBLE(V(IW))).LT.1.D-4) V(IW)=DCMPLX(0.0D0,DIMAG(V(IW)))
      WA=WA+DBLE(V(IW)*CONJG(V(IW)))
      IF(JU.NE.IP) GO TO 38
      IF(IW .EQ.IIWW) GO TO 38
      KKW=KO(IW)
      KO(IW)=KO(IIWW)
      KO(IIWW)=KKW
      WW=V(IW)
      V(IW)=V(IIWW)
      V(IIWW)=WW
   38 CONTINUE
      IND(2,ND1)=IW
      B(ND1)=SQRT(WA)
      IF(WA.GT.1.0D-4) GO TO 39
      IW=IIWW-1
      ND1=ND1-1
      IF(NDR.EQ.1) GO TO 37
      WRITE(6,605)
  605 FORMAT(' STOP AT 39 IN TSPWA')
      STOP
   39 IND(3,ND1)=NDR
      IND(4,ND1)=NB
      D(ND1)=A(IP)
      IF(NDR.NE.1) GO TO 40
      IF(INDIC.EQ.0) THEN
      NP=NP+1
      IF(NP.GT.100)STOP ' REASON NP.GT.100 IN TSPWA'
      NIP(NP)=IP
      ENDIF
      ND2=ND2+1
   40 IF(NDR.GE.ND) GO TO 44
      NDR=NDR+1
      CALL TSIRME(NDR,NB,WD)
      GO TO 47
   44 NI=NI-1
      IF(NI.LE.0) GO TO 21
C     WRITE(6,698) ID,NI
C 698 FORMAT(' POINT 44',2I5)
C
C SET THE NEXT PIVOT
C
   37 IF(INDIC.EQ.1) GO TO 48
      JJWW=IP+1
      IP=IP+1
      IF(IP.GT.JV) GO TO 41
C     DO 41 JU=JJWW,JV
C     IF(KM(4,JU).NE.0) GO TO 41
C     IP=JU
      NDR=1
      CALL TSIRME(1,NB,WD)
      GO TO 47
   41 CONTINUE
      NNP=NP
      IF(ND.NE.2.AND.ND.NE.4) GO TO 46
      IF(ND.EQ.2) NB=2
      IF(ND.EQ.4) NB=4
C     WRITE(6,690) NNP,NP,NDR,NB
      NP=0
   48 NP=NP+1
      IF(NP.GT.NNP) GO TO 46
      IP=NIP(NP)
      INDIC=1
      NDR=1
      CALL TSIRME(1,NB,WD)
      GO TO 47
   46 WRITE(6,604)
  604 FORMAT(' STOP AT 46 IN TSPWA')
      STOP
   21 CONTINUE
      NW=IW
      NNWW=IW
C
      DO 14 I=1,IK
   14 KM(4,I)=IC
C
      RETURN
      ENTRY TSPWDS
      WRITE(6,663) IR,NDI,NK,ND1,ND2
  663 FORMAT(//5I5/)
      IF(NW.EQ.0) GO TO 492
      ID1=0
      DO 49 IW=1,NW
      IF(IW.NE.IND(1,ID1+1).OR.ID1.GE.ND1) GO TO 491
      ID1=ID1+1
      WRITE(6,661) ID1,(IND(K,ID1),K=1,4),D(ID1),B(ID1)
  661 FORMAT(I4,4I5,2F10.5)
      ICOUNT=0
  491 ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KM(K,KO(IW)),K=1,3),IC,V(IW)
  660 FORMAT(4X,2I5,2H (,3I5,2H)/,I3,3H  (,2F10.5,1H))
   49 CONTINUE
      RETURN
  492 WRITE(6,662)
  662 FORMAT(' NO STATE FOR THIS REPRESENTATION')
      RETURN
      END
C SUBROUTINE TSPWAM ====*====3====*====4====*====5====*====6====*====7
C
C      SYMMETRIZE THE PLANE WAVES.
C      ALL COMPONENTS ARE GENRATED.
C      AT 1989.5.17. :  A. YANASE
C      MODIFIED FROM TSPWA
C      MATRIX ELEMENTS OF THE I.R. ARE DIRECTLY OBTAINED FROM
C      IR(7,48,12) IN COMMON/SPG4/ FOR AN ADJUSTMENT TO TSRDSD
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSPWAM(KB,IC,IR,NDI,KO,V,B,D,IND,ND1,ND2,NNWW)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 V,U,WD,CR,CW,WW
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW   / KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      SAVE /SPW/
      DIMENSION KO(1),V(1),B(1),IND(4,1),KB(3),D(1)
      DIMENSION JGA(2,3,48),WD(48),JG(48),CR(48),KK(3),NIP(100)
CC
CC
CC
      PAI2=8.0D0*DATAN(1.0D0)
CC
      CALL CHKDNM(IC)
CC
CC
      DO 13 I=1,IK
   13 KM(4,I)=0
CC
CC
C     WRITE(6,690) IE,IK
C 690 FORMAT(10I5)
      NK=IK
      CALL IRMATA(IR,1,1,MG,ND,JG,JGA,CR,WD)
C     WRITE(6,690) ND,MG
      NDI=ND
      JV=0
      ND1=0
      ND2=0
      IW=0
      DO 21 ID=1,IE
      NN=ND1+1
      JW=JV+1
      JV=KT(ID)
C
C     CHARACTER TEST
C
      CW=0.0D0
      DO 22 II=1,MG
      IF(CDABS(CR(II)).LT.1.0D-4) GO TO 22
      JA=JG(II)
      WW=0.0D0
      DO 23 JU=JW,JV
      JJ=JU
      CALL ZZZY40(JA,JJ,KK)
      WA=0.0D0
      DO 24 K=1,3
      IF(KK(K).NE.KM(K,JU)) GO TO 23
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
      WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
   24 CONTINUE
      X=PAI2*WA
      WW=WW+DCMPLX(COS(X),SIN(X))
   23 CONTINUE
      CW=CW+DCONJG(CR(II))*WW
   22 CONTINUE
      CW=CW/MG
      NI=CW+0.5
C     WRITE(6,690) ID,JW,JV,NI
      IF(NI.EQ.0) GO TO 21
      IP=JW
      NP=0
      NB=1
      NDR=1
      CALL IRMATE(1,NB,WD)
      INDIC=0
C
C     PROJECTION OPERATOR
C
   47 DO 27 JU=JW,JV
   27 U(JU)=0.0
      DO 28 II=1,MG
      IF(CDABS(WD(II)).LT.1.0D-4) GO TO 28
      JA=JG(II)
      CALL ZZZY40(JA,IP,KK)
      WA=0.0D0
      DO 29 K=1,3
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
   29 WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
      X=PAI2*WA
      WW=DCMPLX(COS(X),SIN(X))
      DO 30 JU=JW,JV
      DO 20 K=1,3
      IF(KM(K,JU).NE.KK(K)) GO TO 30
   20 CONTINUE
      KM(4,JU)=IP
      U(JU)=U(JU)+WW*DCONJG(WD(II))
      GO TO 28
   30 CONTINUE
      WRITE(6,602)
  602 FORMAT(' STOP AT 30 IN TSPWAM')
      STOP
   28 CONTINUE
C
C     ORTHOGONALITY TEST
C
C     WRITE(6,699) NN,ND1,ND2,NDR,NB,IP
C 699 FORMAT(' ORTHGONAL TEST',6I5)
      IF(NN.GT.ND1) GO TO 31
      DO 32 II=NN,ND1
C     WRITE(6,699) II
      ID1=IND(1,II)
      ID2=IND(2,II)
      WW=0.0D0
      DO 33 I=ID1,ID2
      JU=JW
   36 IF(KO(I).EQ.JU) GO TO 34
      GO TO 35
   34 CONTINUE
      WW=WW+DCONJG(V(I))*U(JU)
C     WRITE(6,697) WW,I,JU
C 697 FORMAT(' ORTH WW',2F10.6,2I5)
      GO TO 33
   35 JU=JU+1
      IF(JU.LE.JV) GO TO 36
      WRITE(6,603)
  603 FORMAT(' STOP AT 36 IN TSPWAM')
      STOP
   33 CONTINUE
      IF(CDABS(WW).LE.1.0D-4) GO TO 32
      IF(NDR.EQ.1) GO TO 37
      WRITE(6,606)
  606 FORMAT(' STOP AT 32 IN TSPWAM')
      STOP
   32 CONTINUE
C
C     REGISTRATION
C
   31 ND1=ND1+1
      IND(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW+1
      DO 38 JU=JW,JV
      IF(CDABS(U(JU)).GT.1.0D-4) GO TO 381
      IF(JU.EQ.IP) GO TO 381
      IF(KM(4,JU).EQ.IP) KM(4,JU)=0
      GO TO 38
  381 IW=IW+1
      KO(IW)=JU
      V(IW)=U(JU)*(DBLE(ND)/DBLE(MG))
      IF(ABS(DIMAG(V(IW))).LT.1.D-4) V(IW)=DBLE(V(IW))
      IF(ABS(DBLE(V(IW))).LT.1.D-4) V(IW)=DCMPLX(0.0D0,DIMAG(V(IW)))
      WA=WA+DBLE(V(IW)*DCONJG(V(IW)))
      IF(JU.NE.IP) GO TO 38
      IF(IW .EQ.IIWW) GO TO 38
      KKW=KO(IW)
      KO(IW)=KO(IIWW)
      KO(IIWW)=KKW
      WW=V(IW)
      V(IW)=V(IIWW)
      V(IIWW)=WW
   38 CONTINUE
      IND(2,ND1)=IW
      B(ND1)=SQRT(WA)
      IF(WA.GT.1.0D-4) GO TO 39
      IW=IIWW-1
      ND1=ND1-1
      IF(NDR.EQ.1) GO TO 37
      WRITE(6,605)
  605 FORMAT(' STOP AT 39 IN TSPWAM')
      STOP
   39 IND(3,ND1)=NDR
      IND(4,ND1)=NB
      D(ND1)=A(IP)
      IF(NDR.NE.1) GO TO 40
      IF(INDIC.EQ.0) THEN
      NP=NP+1
      IF(NP.GT.100)STOP ' REASON NP.GT.100 IN TSPWAM'
      NIP(NP)=IP
      ENDIF
      ND2=ND2+1
   40 IF(NDR.GE.ND) GO TO 44
      NDR=NDR+1
      CALL IRMATE(NDR,NB,WD)
      GO TO 47
   44 NI=NI-1
      IF(NI.LE.0) GO TO 21
C     WRITE(6,698) ID,NI
C 698 FORMAT(' POINT 44',2I5)
C
C SET THE NEXT PIVOT
C
   37 IF(INDIC.EQ.1) GO TO 48
      JJWW=IP+1
      IP=IP+1
      IF(IP.GT.JV) GO TO 41
C     DO 41 JU=JJWW,JV
C     IF(KM(4,JU).NE.0) GO TO 41
C     IP=JU
      NDR=1
      CALL IRMATE(1,NB,WD)
      GO TO 47
   41 CONTINUE
      NNP=NP
      IF(ND.NE.2.AND.ND.NE.4) GO TO 46
      IF(ND.EQ.2) NB=2
      IF(ND.EQ.4) NB=4
C     WRITE(6,690) NNP,NP,NDR,NB
      NP=0
   48 NP=NP+1
      IF(NP.GT.NNP) GO TO 46
      IP=NIP(NP)
      INDIC=1
      NDR=1
      CALL IRMATE(1,NB,WD)
      GO TO 47
   46 WRITE(6,604)
  604 FORMAT(' STOP AT 46 IN TSPWAM')
      STOP
   21 CONTINUE
      NW=IW
      NNWW=IW
C
      DO 14 I=1,IK
   14 KM(4,I)=IC
C
      RETURN
      ENTRY TSPMDS
      WRITE(6,663) IR,NDI,NK,ND1,ND2
  663 FORMAT(//5I5/)
      IF(NW.EQ.0) GO TO 492
      ID1=0
      DO 49 IW=1,NW
      IF(IW.NE.IND(1,ID1+1).OR.ID1.GE.ND1) GO TO 491
      ID1=ID1+1
      WRITE(6,661) ID1,(IND(K,ID1),K=1,4),D(ID1),B(ID1)
  661 FORMAT(I4,4I5,2F10.5)
      ICOUNT=0
  491 ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KM(K,KO(IW)),K=1,3),IC,V(IW)
  660 FORMAT(4X,2I5,2H (,3I5,2H)/,I3,3H  (,2F10.5,1H))
   49 CONTINUE
      RETURN
  492 WRITE(6,662)
  662 FORMAT(' NO STATE FOR THIS REPRESENTATION')
      RETURN
      END
C SUBROUTINE IRMATA ====*====3====*====4====*====5====*====6====*====7
C
C  NA,NB ELEMENT OF IRREDUCIBLE REPRESENTATION IIR IS GIVEN IN WD
C  NND        :DIMENSION OF REPRESENTATION IIR
C  JJG(48)    :OPERATION CODE
C  JGA(2,3,48):TRANSLATION FOR THE CORRESPONDING JJG
C  CCR(48)    :CHARACTER OF REPRESENTATION IIR
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE IRMATA(IIRR,NA,NB,MMG,NND,JJG,JGA,CCR,WD)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,WD,CCR
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &    ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG2/,/SPG3/,/SPG4/,IIR
      DIMENSION JGA(2,3,48),WD(48),JJG(48),CCR(48)
      MMG=MG
      IIR=IIRR
      CALL CHKNIR(IIR,NR)
      NND=ND(IIR)
      DO 1 I=1,MG
      JJG(I)=IG(JG(I))
      CCR(I)=CR(I,IIR)
      DO 2 K=1,3
      JGA(1,K,I)=JV(1,K,JG(I))
      JGA(2,K,I)=JV(2,K,JG(I))
    2 CONTINUE
    1 CONTINUE
      GO TO 3
C ENTRY IRMATE ====2====*====3====*====4====*====5====*====6====*====7
C   SHORT FORM TO GET NA,NB ELEMENTS, AFTER YOU CALL IRMATA ONCE.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY IRMATE(NA,NB,WD)
    3 CONTINUE
      CALL CHKNST(NA,ND(IIR))
      CALL CHKNST(NB,ND(IIR))
      IWA=32**(6-NB)
      DO 4 I=1,MG
      IW=IR(NA,I,IIR)/IWA
      IW=IW-(IW/32)*32
      IF(IW.NE.0) GO TO 7
      WD(I)=0.D0
      GO TO 4
    7 X=(DBLE(IW)/24.0D0)*2.0D0*3.14159265358979312D0
      WD(I)=DCMPLX(COS(X),SIN(X))/SQRT(DBLE(IR(7,I,IIR)))
    4 CONTINUE
      RETURN
      END
C SUBROUTINE TSFIPW ====*====3====*====4====*====5====*====6====*====7
C
C      FIND THE IRREDUCIBLE REPRESENTATION,JJR, FOR THE EIGEN VECTOR
C      GIVEN BY VR AND VI IN THE CASE OF DOUBLE REPRESENTATION.
C      1985.10.23 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSFIPW(KB,IC,VR,VI,JJR)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 U,WD,CR,WW,SN
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      SAVE /SPW/,/SPG6/
      COMPLEX*16 UU(2,96)
      DIMENSION VR(2,1),VI(2,1)
      DIMENSION KB(3)
      DIMENSION JGA(2,3,48),WD(48),JG(48),CR(48),KK(3)
      DIMENSION ND(12)
CC
      PAI2=8.0D0*DATAN(1.0D0)
CC
C     WRITE(6,690) IE,IK
C 690 FORMAT(10I5)
      CALL CHKDNM(IC)
      CALL TSIRNR(NR,NH,ND)
C     WRITE(6,*) NR,NH,'NR,NH'
CC
      JV=0
      DO 21 ID=1,IE
      JW=JV+1
      JV=KT(ID)
      IF(JV-JW+1.GT.96) THEN
         WRITE(6,*) ' STOP AT DO 21 IN TSFIPW'
         STOP
      END IF
C
      DO 25 IIR=1,NR
      IR=IIR
C     WRITE(6,*) IR
C----------------------------------------------------
C    PROJECTION
      IU=0
      DO 27 JU=JW,JV
      IU=IU+1
      UU(1,IU)=0.0
      UU(2,IU)=0.0
   27 CONTINUE
      CALL TSIRMR(IR,1,1,MG,NND,JG,JGA,CR,WD)
      DO 22 II=1,MG
      IF(CDABS(CR(II)).LT.1.0D-4) GO TO 22
      JA=JG(II)
      JB=MOD(JA-1,24)+1
      DO 23 JU=JW,JV
      JJ=JU
      IF((ABS(VR(1,JU)).LT.0.000001D0).AND.
     &   (ABS(VI(1,JU)).LT.0.000001D0).AND.
     &   (ABS(VR(2,JU)).LT.0.000001D0).AND.
     &   (ABS(VI(2,JU)).LT.0.000001D0)) GO TO 23
      CALL ZZZY40(JA,JJ,KK)
      DO 28 JJU=JW,JV
      IU=JJU-JW+1
      WA=0.0
      DO 24 K=1,3
      IF(KM(K,JJU).NE.KK(K)) GO TO 28
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
      WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
   24 CONTINUE
      X=PAI2*WA
      WW=DCONJG(CR(II))*DCMPLX(COS(X),SIN(X))
      UU(1,IU)=UU(1,IU)+WW*SN(1,1,JB)*DCMPLX(VR(1,JU),VI(1,JU))
     &                 +WW*SN(1,2,JB)*DCMPLX(VR(2,JU),VI(2,JU))
      UU(2,IU)=UU(2,IU)+WW*SN(2,1,JB)*DCMPLX(VR(1,JU),VI(1,JU))
     &                 +WW*SN(2,2,JB)*DCMPLX(VR(2,JU),VI(2,JU))
      GO TO 23
   28 CONTINUE
       WRITE(6,*) ' STOP AT 28 IN TSFIPW'
       STOP
   23 CONTINUE
   22 CONTINUE
C
C    END OF PROJECTION
C-------------------------------------------
C  IF THE PROJECTION GIVES NON ZERO STATE
C  THEN I.R. IS IR
C
      IU=0
      DO 29 JU=JW,JV
      IU=IU+1
      IF(CDABS(UU(1,IU)).GT.0.0001D0) GO TO 30
      IF(CDABS(UU(2,IU)).GT.0.0001D0) GO TO 30
   29 CONTINUE
   25 CONTINUE
   21 CONTINUE
C----------------------------------------------------
C  LOOP END. PROGRAM DOES NOT EXPECT REACH HERE
C
      WRITE(6,*) ' STOP AT 21 IN TSFIPW'
      STOP
C----------------------------------------------------
   30 CONTINUE
      JJR=IR
      RETURN
      END
C SUBROUTINE DGTRST ====*====3====*====4====*====5====*====6====*====7
C
C   JDUB=0 WITHOUT SPIN-ORBIT
C   JDUB=1 WITH SPIN-ORBIT
C   NUMBER OF I.R.     NR
C   NUMBER OF ELEMNTS OF KP-GROUP MG
C   NUMBER OF STAR     NSTAR
C   DEGENERACY OF I.R. NDEG(12)
C   HERING'S SUM       JTR(12)
C   PARTNER            IPART(12)
C   ARE OBTAINED
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE DGTRST(JDUB,NNR,MMG,NSTR,NDEG,JTR,IPART)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/STK/
      DIMENSION IPART(12),NDEG(12),JTR(12)
      CALL FIDPRT(IPART)
      JDUB=IDOUB
      NNR=NR
      NSTR=NS
      MMG=MG
      DO 1 I=1,NR
      NDEG(I)=ND(I)
      JTR(I)=IATR(I)
    1 CONTINUE
C     WRITE(6,601) JDUB,NNR,NSTR
C 601 FORMAT(' JDUB=',I2,' NR=',I3,' NSTR=',I3)
C     WRITE(6,602) (I,I=1,NR)
C 602 FORMAT(' NO  ',12I5)
C     WRITE(6,603) (NDEG(I),I=1,NR)
C 603 FORMAT(' NDEG',12I5)
C     WRITE(6,604) (JTR(I),I=1,NR)
C 604 FORMAT(' JTR ',12I5)
C     WRITE(6,605) (IPART(I),I=1,NR)
C 605 FORMAT(' PRTN',12I5)
      RETURN
      END
C SUBROUTINE GENLTP ====*====3====*====4====*====5====*====6====*====7
C
C    SLAVE ROUTINE FOR FIDPRT. LATTICE VECTORS ARE GENERATED
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE GENLTP(LAT,ICD)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/,JP,JPER,IPER
      DIMENSION JP(3)
      INTEGER IPER(3,6),LAT(3)
      DATA IPER/ 1,2,3, 2,3,1, 3,1,2, 1,3,2, 3,2,1, 2,1,3/
      IF(ICD.EQ.0) THEN
         JP(1)=0
         JP(2)=0
         JP(3)=0
         JPER=1
      END IF
      IX=JP(IPER(1,JPER))
      IY=JP(IPER(2,JPER))
      IZ=JP(IPER(3,JPER))
      IF(IL.EQ.-1) THEN
          LAT(1)=IX*2-IY-IZ
          LAT(2)=IX+IY-IZ*2
          LAT(3)=IX+IY+IZ
          ICD=3
      ELSE IF(IL.EQ.0.OR.IL.EQ.1) THEN
          LAT(1)=IX
          LAT(2)=IY
          LAT(3)=IZ
          ICD=1
      ELSE IF(IL.EQ.2) THEN
          LAT(1)=IY+IZ
          LAT(2)=IZ+IX
          LAT(3)=IX+IY
          ICD=2
      ELSE IF(IL.EQ.3) THEN
          LAT(1)=-IX+IY+IZ
          LAT(2)= IX-IY+IZ
          LAT(3)= IX+IY-IZ
          ICD=2
      ELSE IF(IL.EQ.4) THEN
          LAT(1)=IX-IY
          LAT(2)=IX+IY
          LAT(3)=IZ*2
          ICD=2
      END IF
      IF(JP(2).EQ.JP(3)) THEN
          IF(JP(1).EQ.JP(2)) THEN
               JP(1)=JP(1)+1
               JP(2)=0
               JP(3)=0
               JPER=1
          ELSE
               IF(JPER.LT.3) THEN
                  JPER=JPER+1
               ELSE
                  JP(2)=JP(2)+1
                  JP(3)=0
                  JPER=1
               END IF
          END IF
      ELSE
          IF(JP(1).EQ.JP(2)) THEN
               IF(JPER.LT.3) THEN
                  JPER=JPER+1
               ELSE
                  JPER=1
                  JP(3)=JP(3)+1
               END IF
          ELSE
               IF(JPER.LT.6) THEN
                  JPER=JPER+1
               ELSE
                  JP(3)=JP(3)+1
                  JPER=1
               END IF
          END IF
      END IF
C     WRITE(6,*) ' JP,JPER',JP,JPER
      RETURN
      END
C SUBROUTINE TSBZPL ====*====3====*====4====*====5====*====6====*====7
C
C   PLOT PROGRAM OF THE  FIRST B.Z.
C
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSBZPL
        IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NRECPT(3,30),RECTAX(4,30)
      DIMENSION CORNER(3,2,100)
      CALL FIBZPL(NRECPT,NRP)
      CALL RECTAG(NRECPT,RECTAX,NRP)
      WRITE(6,901)
  901 FORMAT(/' RECIPROCAL LATTICE VECTOR FOR THE B.Z. FORMATION'
     &       /' UNIT IS A*, X-AXIS//A* Y-AXIS IS IN A* B* PLANE')
      WRITE(6,902)
  902 FORMAT('  NO ','   A* B* C*',
     &     '   RECTANGULAR COORDINATE',7X,'LENGTH')
      DO 41 J=1,NRP
      WRITE(6,900) J,(NRECPT(I,J),I=1,3),(RECTAX(I,J),I=1,4)
  900 FORMAT(I4,2X,3I3,4F10.5)
   41 CONTINUE
      NLIN=0
      DO 31 N1=1,NRP-1
      NN1=N1
      DO 32 N2=N1+1,NRP
      NN2=N2
      NNCO=0
      DO 33 N3=1,NRP
      NN3=N3
C     WRITE(6,*) N1,N2,N3
      IF(N3.EQ.N1.OR.N3.EQ.N2) GO TO 33
         CALL FICORN(NN1,NN2,NN3,RECTAX,NRP,X,Y,Z,IND)
         IF(IND.EQ.0) THEN
            IF(NNCO.EQ.0) THEN
               NNCO=1
               NLIN=NLIN+1
               CORNER(1,1,NLIN)=X
               CORNER(2,1,NLIN)=Y
               CORNER(3,1,NLIN)=Z
            ELSE
               WA=(CORNER(1,1,NLIN)-X)**2
     &           +(CORNER(2,1,NLIN)-Y)**2
     &           +(CORNER(3,1,NLIN)-Z)**2
               IF(WA.GT.0.000001D0) THEN
                   CORNER(1,2,NLIN)=X
                   CORNER(2,2,NLIN)=Y
                   CORNER(3,2,NLIN)=Z
                 WRITE(6,601) NLIN,((CORNER(K,I,NLIN),K=1,3),I=1,2)
  601              FORMAT(1H ,I4,' FROM(',3F9.5,') TO(',3F9.5,')')
                   GO TO 32
               END IF
            END IF
         END IF
   33 CONTINUE
      IF(NNCO.EQ.1) THEN
         WRITE(6,*) ' STOP AT 33 IN TSBZPL',N1,N2
      END IF
   32 CONTINUE
   31 CONTINUE
      RETURN
      END
C SUBROUTINE RECTAG ====*====3====*====4====*====5====*====6====*====7
C
C   RECTANGULAR COORDINATES FOR THE RECIPROCAL LATTICE POINTS
C   GIVEN BY NRECPT(3,) ARE CALCULATED IN RECTAX(3,).
C   X-AXIS IS IN A*-AXIS.
C   Y-AXIS IS IN A*,B* PLANE
C   UNIT OF LENGTH IS LATTICE PARAMETER A*
C
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE RECTAG(NRECPT,RECTAX,NRP)
        IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NRECPT(3,30),RECTAX(4,30)
      CALL GETRCP(AK,BK,CK,CAK,CBK,CCK)
      F1=1.0
      F2=(BK/AK)*CCK
      F3=(CK/AK)*CBK
      WSIN=SQRT(1.0D0-CCK*CCK)
      F4=(BK/AK)*WSIN
      F5=(CK/AK)*(CAK-CBK*CCK)/WSIN
      F6=(CK/AK)*SQRT(1.0D0-CAK*CAK-CBK*CBK-CCK*CCK
     &   +2.0D0*CAK*CBK*CCK)/WSIN
      DO 1 I=1,NRP
         RECTAX(1,I)=F1*DBLE(NRECPT(1,I))+F2*DBLE(NRECPT(2,I))
     &              +F3*DBLE(NRECPT(3,I))
         RECTAX(2,I)=F4*DBLE(NRECPT(2,I))+F5*DBLE(NRECPT(3,I))
         RECTAX(3,I)=F6*DBLE(NRECPT(3,I))
         RECTAX(4,I)=SQRT(RECTAX(1,I)**2
     &              +RECTAX(2,I)**2+RECTAX(3,I)**2)
    1 CONTINUE
      RETURN
      END
C SUBROUTINE FICORN ====*====3====*====4====*====5====*====6====*====7
C
C   COMMON POINT OF THE THREE PLANE GIVEN BY N1,N2,N3 IS CALCULATED.
C   WHEN EQUATION IS SINGULAR       THEN IND=2
C   WHEN THE POINT IS OUT OF B.Z,   THEN IND=1
C   ELSE                                 IND=0
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FICORN(N1,N2,N3,RECTAX,NRP,X,Y,Z,IND)
        IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RECTAX(4,30)
      DIMENSION A(3,3),B(3),IPVT(3)
C     WRITE(6,*) N1,N2,N3
         A(1,1)=RECTAX(1,N1)
         A(1,2)=RECTAX(2,N1)
         A(1,3)=RECTAX(3,N1)
         B(1)=0.5*(A(1,1)**2+A(1,2)**2+A(1,3)**2)
         A(2,1)=RECTAX(1,N2)
         A(2,2)=RECTAX(2,N2)
         A(2,3)=RECTAX(3,N2)
         B(2)=0.5*(A(2,1)**2+A(2,2)**2+A(2,3)**2)
         A(3,1)=RECTAX(1,N3)
         A(3,2)=RECTAX(2,N3)
         A(3,3)=RECTAX(3,N3)
         B(3)=0.5*(A(3,1)**2+A(3,2)**2+A(3,3)**2)
      CALL LINEA3(A,B,IPVT,IERR)
c      WRITE(6,*) IERR,n1,n2,n3
      IF(IERR.NE.0) THEN
         IND=2
         RETURN
      ELSE
c      write(6,*) B
         WA=B(1)**2+B(2)**2+B(3)**2
c         INDFAR=1
         DO 1 N=1,NRP
c         IF(N.EQ.N1.OR.N.EQ.N2.OR.N.EQ.N3) GO TO 11
         IF(N.EQ.N1.OR.N.EQ.N2.OR.N.EQ.N3) GO TO 1
         WB=(B(1)-RECTAX(1,N))**2
     &     +(B(2)-RECTAX(2,N))**2
     &     +(B(3)-RECTAX(3,N))**2
c       WRITE(6,*) WA,WB
        IF(WA.GT.WB+0.00001D0) THEN
           IND=1
           RETURN
        END IF
C   11   IF(RECTAX(4,N).GT.WA) INDFAR=0
    1   CONTINUE
c        IF(INDFAR.EQ.1) THEN
c           IND=1
c           RETURN
c        END IF
      END IF
      IND=0
      X=B(1)
      Y=B(2)
      Z=B(3)
      RETURN
      END
C SUBROUTINE FIBZPL ====*====3====*====4====*====5====*====6====*====7
C
C   LIST UP IN NRECPT(3,) THE RECIPROCAL POINTS WHICH CORRESPOND
C   TO A PLAN OF THE FIRST B.Z.
C
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FIBZPL(NRECPT,NRP)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
      DIMENSION NRECPT(3,1)
      IF(IL.EQ.-1) THEN
         KXM=3
         KYM=3
         KZM=3
      ELSE IF(IL.EQ.1.OR.IL.EQ.0) THEN
         KXM=1
         KYM=1
         KZM=1
      ELSE IF(IL.EQ.2.OR.IL.EQ.3) THEN
         KXM=2
         KYM=2
         KZM=2
      ELSE IF(IL.EQ.4) THEN
         KXM=2
         KYM=2
         KZM=1
      END IF
      IF(NG.LE.2) THEN
         KXM=2
         KYM=2
         KZM=2
      END IF
      NRP=0
      DO 45 IZ=0,KZM
          IF(IZ.EQ.0) I3M=1
          IF(IZ.NE.0) I3M=2
      DO 43 IY=0,KYM
          IF(IY.EQ.0) I2M=1
          IF(IY.NE.0) I2M=2
      DO 41 IX=0,KXM
          IF(IX.EQ.0) I1M=1
          IF(IX.NE.0) I1M=2
      IF(IX.EQ.0.AND.IY.EQ.0.AND.IZ.EQ.0) GO TO 41
      DO 42 I3=1,I3M
      DO 44 I2=1,I2M
      DO 46 I1=1,I1M
         KX=IX*(3-2*I1)
         KY=IY*(3-2*I2)
         KZ=IZ*(3-2*I3)
         IF(IL.EQ.0) GO TO 460
         IF(IL.EQ.1) GO TO 460
         IF(IL.EQ.2) GO TO 462
         IF(IL.EQ.3) GO TO 463
         IF(IL.EQ.-1) GO TO 461
         IF(MOD(IABS(KX+KY),2).NE.0) GO TO 46
         GO TO 460
  461    KSUM=-KX+KY+KZ
         IF(MOD(IABS(KSUM),3).NE.0) GO TO 46
         GO TO 460
  463    IF(MOD(IABS(KX+KY+KZ),2).NE.0) GO TO 46
         GO TO 460
  462    IF(MOD(IABS(KX),2).NE.MOD(IABS(KY),2)) GO TO 46
         IF(MOD(IABS(KY),2).NE.MOD(IABS(KZ),2)) GO TO 46
  460 CONTINUE
      CALL TSKFBZ(KX,KY,KZ,2,IND)
C------------------------------------------------------
C   MID POINT TO KX,KY,KZ IS ON A FACE OF THE B.Z.
      IF(IND.EQ.2) THEN
C------------------------------------------------------
C        WRITE(6,*) KX,KY,KZ
         NRP=NRP+1
         NRECPT(1,NRP)=KX
         NRECPT(2,NRP)=KY
         NRECPT(3,NRP)=KZ
      END IF
   46 CONTINUE
   44 CONTINUE
   42 CONTINUE
   41 CONTINUE
   43 CONTINUE
   45 CONTINUE
      RETURN
      END
C SUBROUTINE TSKFBZ ====*====3====*====4====*====5====*====6====*====7
C
C    IF POINT (KX/ICC, KY/ICC, KZ/ICC) IS
C                    OUT OF FIRST B.Z.    IND=0
C                    IN THE FIRST B.Z.    IND=1
C                    ON THE FACE OF B.Z.  IND=2
C                    ON THE EDGE OR CONER IND>3
C          IND IS THE NUMBER OF THE RECIPROCAL LATTICE POINTS
C          WHICH HAVE THE SAME DISTANCE TO THE POINT
C               1988.10.18  AKIRA YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSKFBZ(KKX,KKY,KKZ,IC,IND)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
      CALL CHKDNM(IC)
      IF(IL.EQ.-1) THEN
         KXM=3
         KYM=3
         KZM=3
      ELSE IF(IL.EQ.1.OR.IL.EQ.0) THEN
         KXM=1
         KYM=1
         KZM=1
      ELSE IF(IL.EQ.2.OR.IL.EQ.3) THEN
         KXM=2
         KYM=2
         KZM=2
      ELSE IF(IL.EQ.4) THEN
         KXM=2
         KYM=2
         KZM=1
      END IF
      IF(NG.LE.2) THEN
         KXM=2
         KYM=2
         KZM=2
      END IF
      CALL ZZZY37(KKX,KKY,KKZ,IC,S0)
      NPOINT=0
      DO 41 IX=0,KXM
          IF(IX.EQ.0) I1M=1
          IF(IX.NE.0) I1M=2
      DO 42 I1=1,I1M
         KX=IX*IC*(3-2*I1)
      DO 43 IY=0,KYM
          IF(IY.EQ.0) I2M=1
          IF(IY.NE.0) I2M=2
      DO 44 I2=1,I2M
         KY=IY*IC*(3-2*I2)
      DO 45 IZ=0,KZM
          IF(IZ.EQ.0) I3M=1
          IF(IZ.NE.0) I3M=2
      IF(IX.EQ.0.AND.IY.EQ.0.AND.IZ.EQ.0) GO TO 45
      DO 46 I3=1,I3M
         KZ=IZ*IC*(3-2*I3)
         IF(IL.EQ.0) GO TO 460
         IF(IL.EQ.1) GO TO 460
         IF(IL.EQ.2) GO TO 462
         IF(IL.EQ.3) GO TO 463
         IF(IL.EQ.-1) GO TO 461
         IF(MOD(IABS(KX+KY),2*IC).NE.0) GO TO 46
         GO TO 460
  461    KSUM=-KX+KY+KZ
         IF(MOD(IABS(KSUM),3*IC).NE.0) GO TO 46
         GO TO 460
  463    IF(MOD(IABS(KX+KY+KZ),2*IC).NE.0) GO TO 46
         GO TO 460
  462    IF(MOD(IABS(KX),2*IC).NE.MOD(IABS(KY),2*IC)) GO TO 46
         IF(MOD(IABS(KY),2*IC).NE.MOD(IABS(KZ),2*IC)) GO TO 46
  460 CONTINUE
         KDX=KX-KKX
         KDY=KY-KKY
         KDZ=KZ-KKZ
      CALL ZZZY37(KDX,KDY,KDZ,IC,W)
      IF(W+0.000001D0.LT.S0) GO TO 47
      IF(ABS(W-S0).LT.0.000001D0) NPOINT=NPOINT+1
   46 CONTINUE
   45 CONTINUE
   44 CONTINUE
   43 CONTINUE
   42 CONTINUE
   41 CONTINUE
      IND=NPOINT+1
      RETURN
   47 IND=0
      RETURN
      END
C SUBROUTINE TSKPGH ====*====3====*====4====*====5====*====6====*====7
C
C   KPOINTS IN THE FIRST B.Z. ARE GENERATED IN KKM.
C   FOR HEXAGONAL LATTICE WITH FULL SYMMTERY
C   ICC IS COMMON DENOMINATOR AND NK IS NUMBER OF KPOINT.
C   AXIS OF GAMMA TO M IN THE X-DIRECTION
C        IS DIVIDED BY NX.
C   AXIS OF M TO K IN THE DIRECTION PERPENDICULAR TO X
C        IS DIVIDED BY NY=NX.
C   AXIS OF GAMMA TO Z IN THE Z-DIRECTION
C        IS DIVIDED BY NZ.
C
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSKPGH(NX,NY,NZ,KKM,ICC,NK)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
C
      DIMENSION KKM(3,1)
      CALL ADDINV(NTLATC)
      IF(IL.NE.0.OR.NTLATC.NE.0) THEN
         WRITE(6,*) ' TSKPGH IS SPECIALLY DESIGNED TO GENERATE KPOINT'
         WRITE(6,*) ' FOR HEXAGONAL LATTICE WITH FULL SYMMETRY.'
         WRITE(6,*) ' PLEASE CALL TSKPGN FOR OTHER LATTICE.'
         WRITE(6,*) ' STOP IN TSKPGH'
         STOP
      END IF
      IF(NX.NE.NY) THEN
         WRITE(6,*) ' NX SHOULD BR EQUAL TO NY FOR TSKPGH'
         NY=NX
      END IF
      ICXY=NX*2*3
      CALL ZZZY53(ICXY,NZ,KK)
      LCM=ICXY*(NZ/KK)
      IF(MOD((LCM/NZ),2).EQ.0) THEN
           KDX=LCM/(NX*2)
           KDY=LCM/(NY*3)
           KDXX=KDY/2
           KDZ=(LCM/NZ)/2
           ICC=LCM
      ELSE
           ICC=2*LCM
           KDX=LCM/NX
           KDXX=LCM/(NY*3)
           KDY=KDXX*2
           KDZ=LCM/NZ
      END IF
C     WRITE(6,*) NX,NY,NZ,LCM
C     WRITE(6,901)
C 901 FORMAT(/' GENERATED KPOINT'
C    &       /'    NO  KX, KY, KZ')
      NK=0
      IND32=0
      DO 1 IZ=0,NZ
         KKZ=IZ*KDZ
      DO 2 IX=0,NX
      DO 3 IY=0,IX
         KKX=IX*KDX-IY*KDXX
         KKY=IY*KDY
         CALL TSKFBZ(KKX,KKY,KKZ,ICC,IND)
C        WRITE(6,601) KKX,KKY,KKZ,ICC,IND
C 601    FORMAT(' B.Z. TEST',4I6)
         IF(IND.EQ.0) GO TO 3
         CALL KALRST(KKX,KKY,KKZ,ICC,KKM,NK,IND)
C        WRITE(6,602) KKX,KKY,KKZ,ICC,IND
C 602    FORMAT(' TEST IF ALREADY REGISTERED',4I5)
         IF(IND.EQ.1) THEN
            NK=NK+1
            KKM(1,NK)=KKX
            KKM(2,NK)=KKY
            KKM(3,NK)=KKZ
C           WRITE(6,900) NK,KKX,KKY,KKZ,ICC,IND
C 900 FORMAT(1H ,I5,3I4,'/',I4,' IND=',I2)
            IF(IND32.EQ.0) THEN
                IF(3*KKX.EQ.ICC.AND.6*KKY.EQ.ICC) IND32=1
            END IF
         END IF
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
      IF(IND32.EQ.0) THEN
         KKX=ICC/3
         KKY=ICC/6
         DO 4 IZ=0,NZ
         KKZ=IZ*KDZ
         NK=NK+1
         KKM(1,NK)=KKX
         KKM(2,NK)=KKY
         KKM(3,NK)=KKZ
C           WRITE(6,900) NK,KKX,KKY,KKZ,ICC,IND
    4    CONTINUE
      END IF
      CALL REMINV
      RETURN
      END
C SUBROUTINE TSNMKP ====*====3====*====4====*====5====*====6====*====7
C
C  THE NAME OF KPOINT IN THE B.Z. IS GIVEN
C  
C                  BY A. YANASE
C  MODIFIED BY A.YANASE    1944/03/18
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSNMKP(KX,KY,KZ,ICC,CIR)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      CHARACTER*2 CIR
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/LAT/A,B,C,CA,CB,CC,A1,B1,C1
      SAVE /SPG2/,/LAT/
      CALL CHKDNM(ICC)
      CALL ADDINV(NLATIC)
      IF(IL.EQ.-1) THEN
         CALL TRIGON(KX,KY,KZ,ICC,CIR)
      ELSE IF(IL.EQ.0) THEN
         CALL HEXAGN(KX,KY,KZ,ICC,CIR)
      ELSE IF(IL.EQ.1) THEN
          IF(NLATIC.EQ.8.OR.NLATIC.EQ.7) THEN
                CALL SIMCUB(KX,KY,KZ,ICC,CIR)
C          ELSE IF(NLATIC.EQ.7) THEN
C                CALL SIMPTH(KX,KY,KZ,ICC,CIR)
          ELSE IF(NLATIC.EQ.6.OR.NLATIC.EQ.5) THEN
                CALL SIMTET(KX,KY,KZ,ICC,CIR)
          ELSE IF(NLATIC.EQ.4) THEN
                CALL SIMORT(KX,KY,KZ,ICC,CIR)
          ELSE IF(NLATIC.EQ.3.OR.NLATIC.EQ.2.OR.NLATIC.EQ.9) THEN
                CALL SIMMON(KX,KY,KZ,ICC,NLATIC,CIR)
          ELSE IF(NLATIC.EQ.1) THEN
                CALL SIMTRI(KX,KY,KZ,ICC,CIR)
          ELSE
                WRITE(6,*) ' WE ARE VERY SORRY '
                WRITE(6,*) ' WE HAVE NOT YET PREPARED THE SUBROUTINE '
                CIR='GN'
          END IF
      ELSE IF(IL.EQ.2) THEN
          IF(NLATIC.EQ.8.OR.NLATIC.EQ.7) THEN
                CALL FACCUB(KX,KY,KZ,ICC,CIR)
C          ELSE IF(NLATIC.EQ.7) THEN
C                CALL FACDTH(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.6.OR.NLATIC.EQ.5) THEN
                 CALL FACTET(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.4) THEN
                 CALL FACORT(KX,KY,KZ,ICC,CIR)
           ELSE
                 WRITE(6,*) ' WE ARE VERY SORRY '
                 WRITE(6,*) ' WE HAVE NOT YET PREPARED THE SUBROUTINE '
                 CIR='GN'
           END IF
       ELSE IF(IL.EQ.3) THEN
           IF(NLATIC.EQ.8.OR.NLATIC.EQ.7) THEN
                 CALL BDCCUB(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.6.OR.NLATIC.EQ.5) THEN
                 CALL BDCTET(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.4) THEN
                 CALL BDCORT(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.9.OR.NLATIC.EQ.3.OR.NLATIC.EQ.2) THEN
                 CALL BDCMON(KX,KY,KZ,ICC,NLATIC,CIR)
           ELSE
                 WRITE(6,*) ' WE ARE VERY SORRY '
                 WRITE(6,*) ' WE HAVE NOT YET PREPARED THE SUBROUTINE '
                 CIR='GN'
           END IF
       ELSE IF(IL.EQ.4) THEN
           IF(NLATIC.EQ.4) THEN
                 CALL BACORT(KX,KY,KZ,ICC,CIR)
           ELSE IF(NLATIC.EQ.9.OR.NLATIC.EQ.3) THEN
                 CALL BACMON(KX,KY,KZ,ICC,NLATIC,CIR) 
           ELSE 
                 WRITE(6,*) ' WE ARE VERY SORRY '
                 WRITE(6,*) ' WE HAVE NOT YET PREPARED THE SUBROUTINE '
                 CIR='GN'
          END IF
       END IF
       CALL REMINV
       RETURN
       END
C SUBROUTINE FACTET ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR FACE CENTERED TETRAGONAL LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 2002.06.24 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FACTET(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
        IF(KY.EQ.0) THEN 
          IF(KX.EQ.0) THEN
                  CN='GM'
          ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
          ELSE
                  CN='DT'
          END IF
        ELSE IF(KY.EQ.KX) THEN
          IF(KX.EQ.ICC) THEN
                  CN='Z '
          ELSE
                  CN='SM'
          END IF
        ELSE IF(KX.EQ.ICC) THEN
                  CN='Y '
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KX.EQ.0) THEN
                  CN='Z '
          ELSE IF(KY.EQ.0) THEN
                  CN='U '
          ELSE IF(KX.EQ.KY) THEN
                  CN='F '
          ELSE
                  CN='ZB'
          END IF
      ELSE IF(KZ*2.EQ.ICC.AND.KX+KY.EQ.ICC) THEN
          IF(KY.EQ.KX) THEN
                  CN='N '
          ELSE IF(KX.EQ.ICC) THEN
                  CN='P '
          ELSE
                  CN='Q '
          END IF
      ELSE
          IF(KX.EQ.0) THEN
                  CN='LD'
          ELSE IF(KY.EQ.0.AND.KX.EQ.ICC) THEN
                  CN='W '
          ELSE IF(KX.EQ.ICC.AND.KY.EQ.KX) THEN
                  CN='V '
          ELSE IF(KX.EQ.KY) THEN
                  CN='XY'
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='XB'
          ELSE IF(KY.EQ.0) THEN
                  CN='YP'
          ELSE
                  CN='GN'
          END IF
      END IF
C         write(6,*) KX,KY,KZ,ICC,CN
      RETURN
      END
CSUBROUTINE SIMTRI ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR SIMPLE TRICLINIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMTRI(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      CN='GN'
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='X '
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='S '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Y '
              END IF
          ELSE IF(KY.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='X '
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='U '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Z '
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='T '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='R '
              END IF
          ELSE IF(KY.EQ.ICC.AND.KX.EQ.0) THEN
                  CN='Z '
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KY.EQ.0.AND.KX*2.EQ.ICC) THEN
                 CN='X '
          ELSE IF(KY*2.EQ.ICC.AND.KX.EQ.0) THEN
                 CN='Y '
          END IF
      END IF
      RETURN
      END
C SUBROUTINE SIMMON ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR SIMPLE MONOCLINIC LATTICE
C   ICC :COMMON DENOMINATOR
C   NLATIC=9   TWO-FOLD AXIS // A
C   NLATIC=3   TWO-FOLD AXIS // B
C   NLATIC=2   TWO-FOLD AXIS // C
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMMON(KKX,KKY,KKZ,ICC,NLATIC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      CN='GN'
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='X '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='LD'
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='A '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='V '
              END IF
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='LD'
              ELSE IF(KX*2.EQ.ICC.AND.NLATIC.EQ.3) THEN
                  CN='V '
              ELSE IF(NLATIC.EQ.2) THEN
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='D '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='W '
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='C '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='E '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='U '
              END IF
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='W '
              ELSE IF(KX*2.EQ.ICC.AND.NLATIC.EQ.3) THEN
                  CN='U '
              ELSE IF(NLATIC.EQ.2) THEN
                  CN='ZB'
              END IF
          END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0.AND.NLATIC.EQ.2) THEN
                  CN='LD'
              ELSE IF(KX*2.EQ.ICC.AND.NLATIC.EQ.2) THEN
                  CN='V '
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='YP'
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0.AND.NLATIC.EQ.2) THEN
                  CN='W '
              ELSE IF(KX*2.EQ.ICC.AND.NLATIC.EQ.2) THEN
                  CN='U '
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='YB'
              END IF
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.9) THEN
                  CN='XP'
              ELSE IF(KX*2.EQ.ICC.AND.NLATIC.EQ.9) THEN
                  CN='XB'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE BACMON ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BASE CENTERED MONOCLINIC LATTICE
C   ICC :COMMON DENOMINATOR
C   NLATIC=9   TWO-FOLD AXIS // A
C   NLATIC=3   TWO-FOLD AXIS // B
C   
C   CN  :NAME
C                 1994.03.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BACMON(KKX,KKY,KKZ,ICC,NLATIC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      CN='GN'
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Y '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='LD'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='C '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='V '
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC.AND.NLATIC.EQ.3) THEN
                  CN='C '
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='M '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='U '
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='M '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='E '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='L '
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='U '
              ELSE IF(KX.EQ.ICC.AND.NLATIC.EQ.3) THEN
                  CN='E '
              END IF
          END IF
      ELSE IF(KY.EQ.0.AND.NLATIC.EQ.3) THEN
         CN='YP'
      ELSE IF(KY.EQ.ICC.AND.NLATIC.EQ.3) THEN
         CN='YB'
      ELSE IF(KX.EQ.0.AND.NLATIC.EQ.9) THEN
         CN='XP'
      ELSE IF(KX.EQ.ICC.AND.NLATIC.EQ.9) THEN
         CN='XB' 
      END IF
      RETURN
      END
C SUBROUTINE BDCMON ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BODY CENTERED ORTHORHOMBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   NLATIC=9   TWO-FOLD AXIS // A
C   NLATIC=3   TWO-FOLD AXIS // B
C   NLATIC=2   TWO-FOLD AXIS // C
C   CN  :NAME
C                 1994.03.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BDCMON(KKX,KKY,KKZ,ICC,NLATIC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      CN='GN'
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='LD'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='F '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='T '
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC.AND.NLATIC.EQ.3) THEN
                  CN='F '
              ELSE IF(NLATIC.EQ.2) THEN
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='F '
              END IF
          ELSE
              IF(KX.EQ.0.AND.NLATIC.EQ.3) THEN
                  CN='F '
              ELSE IF(NLATIC.EQ.2) THEN
                  CN='ZB'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY*2.EQ.ICC) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='W '
              ELSE IF(KX.EQ.0) THEN
                  CN='S '
              ELSE IF(NLATIC.EQ.9) THEN
                  CN='D '
              END IF
           ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='R '
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='D '
              END IF
           ELSE IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='YP'
              END IF
           ELSE IF(KX.EQ.0.AND.NLATIC.EQ.9) THEN
                  CN='XP'
           END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0.AND.NLATIC.EQ.2) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC.AND.NLATIC.EQ.2) THEN
                  CN='F '
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='YP'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0.AND.NLATIC.EQ.2) THEN
                  CN='F '
              ELSE IF(NLATIC.EQ.3) THEN
                  CN='YB'
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC.AND.NLATIC.EQ.2) THEN
                  CN='D '
          ELSE IF(NLATIC.EQ.9) THEN
              IF(KX.EQ.0) THEN
                  CN='XP'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='XB'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE SIMTET ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BODY SIMPLE TETRAGONAL LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMTET(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
          IF(KX.EQ.0) THEN
                  CN='GM'
          ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='X '
              ELSE IF(KY*2.EQ.ICC) THEN
                  CN='M '
              ELSE
                  CN='Y '
              END IF
          ELSE IF(KY.EQ.0) THEN
                  CN='DT'
          ELSE IF(KX.EQ.KY) THEN
                  CN='SM'
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KX.EQ.0) THEN
                  CN='Z '
          ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='R '
              ELSE IF(KY*2.EQ.ICC) THEN
                  CN='A '
              ELSE
                  CN='T '
              END IF
          ELSE IF(KX.EQ.KY) THEN
                  CN='S '
          ELSE
              IF(KY.EQ.0) THEN
                  CN='U '
              ELSE
                  CN='ZB'
              END IF
          END IF
      ELSE
          IF(KX.EQ.0) THEN
                  CN='LD'
          ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='W '
              ELSE IF(KY*2.EQ.ICC) THEN
                  CN='V '
              ELSE
                  CN='XB'
              END IF
          ELSE IF(KX.EQ.KY) THEN
                  CN='XY'
          ELSE
              IF(KY.EQ.0) THEN
                  CN='XP'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE SIMPTH ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR SIMPLE CUBIC LATTICE WITH TH
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMPTH(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      DO 1 N=1,2
          IF(KZ.GT.KY.OR.KZ.GT.KX) THEN
              IW=KZ
              KZ=KY
              KY=KX
              KX=IW
          ELSE
              GO TO 2
          END IF
    1 CONTINUE
    2 CONTINUE
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='DT'
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='M '
              ELSE
                  CN='Z1'
              END IF
          ELSE IF(KX.EQ.KY) THEN
                  CN='SM'
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='Z2'
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.KY) THEN
          IF(KY.EQ.KX) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='R '
               ELSE
                  CN='LD'
              END IF
          ELSE IF(KX*2.EQ.ICC) THEN
                  CN='S '
          ELSE
                  CN='XY'
          END IF
      ELSE IF(KZ.EQ.KX) THEN
          IF(KY*2.EQ.ICC) THEN
                  CN='S '
          ELSE
                  CN='XY'
          END IF
      ELSE
          IF(KY*2.EQ.ICC) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='T '
              ELSE
                  CN='B1'
              END IF
          ELSE IF(KY.EQ.KX) THEN
                  CN='XY'
          ELSE
              IF(KX*2.EQ.ICC) THEN
                  CN='B2'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE SIMORT ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR SIMPLE ORTHORHOMBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMORT(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='S '
              ELSE
                  CN='C '
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='D '
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='U '
              ELSE
                  CN='A '
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='T '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='R '
              ELSE
                  CN='E '
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='B '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='P '
              ELSE
                  CN='ZB'
              END IF
          END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='G '
              ELSE
                  CN='YP'
              END IF
          ELSE IF(KY*2.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='H '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='Q '
              ELSE
                  CN='YB'
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='XP'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='XB'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE SIMCUB ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR SIMPLE CUBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE SIMCUB(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      DO 1 N=1,2
          IF(KZ.GT.KY.OR.KZ.GT.KX) THEN
              IW=KZ
              KZ=KY
              KY=KX
              KX=IW
          ELSE
              GO TO 2
          END IF
    1 CONTINUE
    2 CONTINUE
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='DT'
              END IF
          ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='X '
              ELSE IF(KY*2.EQ.ICC) THEN
                  CN='M '
              ELSE
                  CN='Z '
              END IF
          ELSE IF(KX.EQ.KY) THEN
                  CN='SM'
          ELSE
              IF(KY.EQ.0) THEN
                  CN='DT'
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.KY) THEN
          IF(KY.EQ.KX) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='R '
              ELSE
                  CN='LD'
              END IF
          ELSE IF(KX*2.EQ.ICC) THEN
                  CN='S '
          ELSE
                  CN='XY'
          END IF
      ELSE
          IF(KX*2.EQ.ICC) THEN
              IF(KY*2.EQ.ICC) THEN
                  CN='T '
              ELSE
                  CN='ZB'
              END IF
          ELSE IF(KY.EQ.KX) THEN
                  CN='XY'
          ELSE
                  CN='GN'
          END IF
      END IF
      RETURN
      END
C SUBROUTINE FACORT ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR FACE CENTERED ORTHORHOMBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FACORT(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Z '
              ELSE
                  CN='C '
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='D '
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Y '
              ELSE
                  CN='A '
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE
                  CN='U '
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='B '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='R '
              ELSE
                  CN='ZB'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC.AND.KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='L '
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='G '
              ELSE
                  CN='YP'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='H '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Q '
              ELSE
                  CN='YB'
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='XP'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='XB'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE FACCUB ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR FACE CENTERED CUBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FACCUB(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      DO 1 N=1,2
          IF(KZ.GT.KY.OR.KZ.GT.KX) THEN
              IW=KZ
              KZ=KY
              KY=KX
              KX=IW
          ELSE
              GO TO 2
          END IF
    1 CONTINUE
    2 CONTINUE
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='DT'
              END IF
          ELSE IF(KX.EQ.ICC) THEN
              IF(KY*2.EQ.ICC) THEN
                  CN='W '
              ELSE
                  CN='Z '
              END IF
          ELSE IF(KX.EQ.KY) THEN
              IF(KX*4.EQ.ICC*3) THEN
                  CN='K '
              ELSE
                  CN='SM'
              END IF
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ.EQ.KY) THEN
          IF(KY.EQ.KX) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='L '
              ELSE
                  CN='LD'
              END IF
          ELSE IF(KX.EQ.ICC) THEN
              IF(KY*4.EQ.ICC) THEN
                  CN='U '
              ELSE
                  CN='S '
              END IF
          ELSE
                  CN='XY'
          END IF
      ELSE
          IF(KX.EQ.ICC) THEN
                  CN='ZB'
          ELSE IF(KY.EQ.KX) THEN
                  CN='XY'
          ELSE IF(KY*2.EQ.ICC.AND.KX+KZ.EQ.ICC) THEN
                  CN='Q '
          ELSE
                  CN='GN'
          END IF
      END IF
      RETURN
      END
C SUBROUTINE BACORT ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BASE CENTERED ORTHORHOMBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BACORT(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Y '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='Y '
              ELSE
                  CN='C '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='S '
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='F '
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='Z '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='T '
              ELSE
                  CN='A '
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='T '
              ELSE
                  CN='E '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='R '
          ELSE
              IF(KX.EQ.0) THEN
                  CN='B '
              ELSE IF(KX.EQ.ICC) THEN
                  CN='G '
              ELSE
                  CN='ZB'
              END IF
          END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='H '
              ELSE
                  CN='YP'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='H '
              ELSE
                  CN='YB'
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='D '
          ELSE
              IF(KX.EQ.0) THEN
                  CN='XP'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='XB'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE BDCORT ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BODY CENTERED ORTHORHOMBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BDCORT(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE
                  CN='F '
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='T '
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='U '
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE
                  CN='F '
              END IF
          ELSE
              IF(KX.EQ.0) THEN
                  CN='U '
              ELSE
                  CN='ZB'
              END IF
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY*2.EQ.ICC) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='W '
              ELSE IF(KX.EQ.0) THEN
                  CN='S '
              ELSE
                  CN='D '
              END IF
           ELSE IF(KX*2.EQ.ICC) THEN
              IF(KY.EQ.0) THEN
                  CN='R '
              ELSE
                  CN='Q '
              END IF
           ELSE IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE
                  CN='YP'
              END IF
           ELSE IF(KX.EQ.0) THEN
                  CN='XP'
           ELSE
                  CN='GN'
           END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='G '
              ELSE
                  CN='YP'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='G '
              ELSE
                  CN='YB'
              END IF
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='P '
          ELSE
              IF(KX.EQ.0) THEN
                  CN='XP'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='XB'
              ELSE
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE BDCCUB ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BODY CENTERED CUBIC LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BDCCUB(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      DO 1 N=1,2
          IF(KZ.GT.KY.OR.KZ.GT.KX) THEN
              IW=KZ
              KZ=KY
              KY=KX
              KX=IW
          ELSE
              GO TO 2
          END IF
    1 CONTINUE
    2 CONTINUE
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='H '
              ELSE
                  CN='DT'
              END IF
          ELSE IF(KX.EQ.KY) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='N '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='G '
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ.EQ.KY) THEN
          IF(KY.EQ.KX) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='P '
              ELSE
                  CN='LD'
              END IF
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='F '
          ELSE
                  CN='XY'
          END IF
      ELSE
          IF(KX.EQ.KY) THEN
              IF(KY*2.EQ.ICC) THEN
                  CN='D '
              ELSE
                  CN='XY'
              END IF
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='ZB'
          ELSE
                  CN='GN'
          END IF
      END IF
      RETURN
      END
C SUBROUTINE BDCTET ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR BODY CENTERED TETRAGONAL LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE BDCTET(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
         IF(KY.GT.KX) THEN
              IW=KY
              KY=KX
              KX=IW
         END IF
      IF(KZ.EQ.0) THEN
          IF(KX.EQ.0) THEN
                  CN='GM'
          ELSE IF(KX.EQ.ICC) THEN
                  CN='Z '
          ELSE IF(KY.EQ.0) THEN
                  CN='SM'
          ELSE IF(KX.EQ.KY) THEN
             IF(KX*2.EQ.ICC) THEN
                  CN='X '
             ELSE
                  CN='DT'
             END IF
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='Y '
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ.EQ.ICC) THEN
          IF(KX.EQ.0) THEN
                  CN='Z '
          ELSE IF(KY.EQ.0) THEN
                  CN='F '
          ELSE IF(KX.EQ.KY) THEN
                  CN='U '
          ELSE
                  CN='ZB'
          END IF
      ELSE IF(KZ*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
                  CN='N '
          ELSE IF(KY*2.EQ.ICC) THEN
                  CN='P '
          ELSE
                  CN='Q '
          END IF
      ELSE
          IF(KX.EQ.0) THEN
                  CN='LD'
          ELSE IF(KY*2.EQ.ICC.AND.KX*2.EQ.ICC) THEN
                  CN='W '
          ELSE IF(KX.EQ.ICC.AND.KY.EQ.0) THEN
                  CN='V '
          ELSE IF(KX.EQ.KY) THEN
                  CN='XY'
          ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='XB'
          ELSE IF(KY.EQ.0) THEN
                  CN='YP'
          ELSE
                  CN='GN'
          END IF
      END IF
      RETURN
      END
C SUBROUTINE HEXAGN ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR HEXAGONAL LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE HEXAGN(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      IF(KKX.GE.0.AND.KKY.GE.0) THEN
             KX=KKX
             KY=KKY
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.LT.0) THEN
             KX=-KKX
             KY=-KKY
             KZ=-KKZ
      ELSE IF(KKX.GE.0.AND.KKY.LT.0.AND.KKX+KKY.GE.0) THEN
             KX=KKX+KKY
             KY=-KKY
             KZ=KKZ
      ELSE IF(KKX.GE.0.AND.KKY.LT.0.AND.KKX+KKY.LT.0) THEN
             KX=-(KKX+KKY)
             KY=KKX
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.GE.0.AND.KKX+KKY.GE.0) THEN
             KX=-KKX
             KY=KKX+KKY
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.GE.0.AND.KKX+KKY.LT.0) THEN
             KX=KKY
             KY=-(KKX+KKY)
             KZ=KKZ
      END IF
      IF(KZ.LT.0) THEN
             IW=KX
             KX=KY
             KY=IW
             KZ=-KZ
      END IF
      IF(KY.GT.KX) THEN
             IW=KY
             KY=KX
             KX=IW
      END IF
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='M '
              ELSE
                  CN='SM'
              END IF
          ELSE IF(KX.EQ.KY) THEN
              IF(KX*3.EQ.ICC) THEN
                  CN='K '
              ELSE
                  CN='T '
              END IF
          ELSE IF(KX*2+KY.EQ.ICC) THEN
                  CN='TP'
          ELSE
                  CN='ZP'
          END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='A '
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='L '
              ELSE
                  CN='R '
              END IF
          ELSE IF(KX.EQ.KY) THEN
              IF(KX*3.EQ.ICC) THEN
                  CN='H '
              ELSE
                  CN='S '
              END IF
          ELSE IF(KX*2+KY.EQ.ICC) THEN
                  CN='SP'
          ELSE
                  CN='ZB'
          END IF
      ELSE
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='LD'
              ELSE IF(KX*2.EQ.ICC) THEN
                  CN='U '
              ELSE
                  CN='YP'
              END IF
          ELSE IF(KX.EQ.KY) THEN
              IF(KX*3.EQ.ICC) THEN
                  CN='P '
              ELSE
                  CN='XY'
              END IF
          ELSE IF(KX*2+KY.EQ.ICC) THEN
                  CN='XB'
          ELSE
                  CN='GN'
          END IF
      END IF
      RETURN
      END
C SUBROUTINE TRIGON ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR TRIGONAL LATTICE
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TRIGON(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      IF(KKX.GE.0.AND.KKY.GE.0) THEN
             KX=KKX
             KY=KKY
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.LT.0) THEN
             KX=-KKX
             KY=-KKY
             KZ=-KKZ
      ELSE IF(KKX.GE.0.AND.KKY.LT.0.AND.KKX+KKY.GE.0) THEN
             KX=KKX+KKY
             KY=-KKY
             KZ=KKZ
      ELSE IF(KKX.GE.0.AND.KKY.LT.0.AND.KKX+KKY.LT.0) THEN
             KX=-(KKX+KKY)
             KY=KKX
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.GE.0.AND.KKX+KKY.GE.0) THEN
             KX=-KKX
             KY=KKX+KKY
             KZ=KKZ
      ELSE IF(KKX.LT.0.AND.KKY.GE.0.AND.KKX+KKY.LT.0) THEN
             KX=KKY
             KY=-(KKX+KKY)
             KZ=KKZ
      END IF
      IF(KZ.LT.0) THEN
             IW=KX
             KX=KY
             KY=IW
             KZ=-KZ
      END IF
      IF(KZ.EQ.0) THEN
         IF(KX.EQ.KY) THEN
             IF(KX.EQ.0) THEN
                  CN='GM'
             ELSE IF(KX*2.EQ.ICC) THEN
                  CN='F '
             ELSE
                  CN='SM'
             END IF
         ELSE IF(KY.EQ.0) THEN
             IF(KX.EQ.ICC) THEN
                  CN='P '
             ELSE
                  CN='YP'
             END IF
         ELSE IF(KX.EQ.0) THEN
             IF(KY.EQ.ICC) THEN
                  CN='P '
             ELSE
                  CN='XP'
             END IF
         ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='ZB'
         ELSE
                  CN='GN'
         END IF
      ELSE IF(KZ*2.EQ.ICC) THEN
         IF(KX*2+KY.EQ.ICC) THEN
             IF(KX.EQ.0) THEN
                 CN='Z '
             ELSE IF(KY.EQ.0) THEN
                 CN='L '
             ELSE
                 CN='Y '
             END IF
         ELSE IF(KY.EQ.0) THEN
             IF(KX.EQ.0) THEN
                  CN='LD'
             ELSE
                  CN='YP'
             END IF
         ELSE IF(KX.EQ.0) THEN
                  CN='XP'
         ELSE
                  CN='GN'
         END IF
      ELSE IF(KZ.EQ.ICC) THEN
         IF(KX+KY*2.EQ.ICC) THEN
             IF(KX.EQ.0) THEN
                  CN='F '
             ELSE
                  CN='Q '
             END IF
         ELSE IF(KY.EQ.0) THEN
             IF(KX.EQ.0) THEN
                  CN='LD'
             ELSE
                  CN='YP'
             END IF
         ELSE IF(KX.EQ.0) THEN
                  CN='XP'
         ELSE
                  CN='GN'
         END IF
      ELSE IF(KZ*2.EQ.3*ICC) THEN
         IF(KX.EQ.KY) THEN
             IF(KX.EQ.0) THEN
                  CN='Z '
             ELSE
                  CN='B '
             END IF
         ELSE IF(KY.EQ.0) THEN
                  CN='YP'
         ELSE IF(KX.EQ.0) THEN
                  CN='XP'
         ELSE
                  CN='GN'
         END IF
      ELSE
         IF(KY.EQ.0) THEN
             IF(KX.EQ.0) THEN
                  CN='LD'
             ELSE IF(KX.EQ.ICC) THEN
                  CN='P '
             ELSE
                  CN='YP'
             END IF
         ELSE IF(KX.EQ.0) THEN
             IF(KY.EQ.ICC) THEN
                  CN='P '
             ELSE
                  CN='XP'
             END IF
         ELSE IF(KX+KY.EQ.ICC) THEN
                  CN='ZB'
         ELSE
                  CN='GN'
         END IF
      END IF
      RETURN
      END
C SUBROUTINE TSRDSD ====*====3====*====4====*====5====*====6====*====7
C
C  REDUCED BASIS FOR THE PRODUCT REPRESENTATION OF SINGLE REPRESENTATION
C  AND SPIN MATRIX ARE OBTAINED IN DB(6,2,6,2)
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSRDSD(JRD,JRS,DB,NI)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 DB,CR,CRS,CW,WWD,WD,WA,CE
      COMPLEX*16 SN
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/SPG5/IH(48),KH(24,2),IO,NN,IPA(48,2)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      SAVE /SPG2/,/SPG3/,/SPG4/,/SPG5/,/SPG6/
      DIMENSION KKB(3),INDA(6,2),INDB(6,2)
      DIMENSION IRS(7,48),CRS(48)
      DIMENSION WWD(48,6),WD(48),DB(6,2,6,2),WA(6,2)
      DIMENSION JJG(48),JGA(2,3,48),CE(48)
      DO 11 K=1,3
      KKB(K)=KB(K)
   11 CONTINUE
      ICC=ICB
      CALL TSIREP(KKB,ICC,0)
      CALL CHKNIR(JRS,NR)
      NDS=ND(JRS)
      DO 12 J=1,NH
          CRS(J)=CR(J,JRS)
          DO 13 K=1,7
            IRS(K,J)=IR(K,J,JRS)
   13     CONTINUE
   12 CONTINUE
      CALL TSIREP(KKB,ICC,1)
      CALL CHKNIR(JRD,NR)
C
C   CHARACTER TEST
C
        CW=0
        DO 23 I=1,NH
          IF(IL.GE.1) THEN
             JJ=MOD(IG(JG(I))-1,24)+1
          ELSE
             JJ=MOD(IG(JG(I))-1,12)+1
          END IF
            CW=CW+DCONJG(CR(I,JRD))*CRS(I)
     &        *(SN(1,1,JJ)+SN(2,2,JJ))
C     WRITE(6,*) CW,I,JRD,JRS
   23   CONTINUE
        CW=CW/NH
        NI=CW+0.5
C     WRITE(6,*) NI
        IF(NI.GT.2) THEN
           WRITE(6,*) ' DIMENSION OF DB IS NOT ENUGH NI=',NI
           STOP ' STOP AT 23 IN TSRDSD '
        END IF
        IF(NI.EQ.0) GO TO 80
      DO 81 NA=1,NDS
      DO 81 NBS=1,2
        INDA(NA,NBS)=0
   81 CONTINUE
      DO 82 NNI=1,NI
C
C   PROJECTION OPERATOR
C
C   PIVOT IS SET AT THE COMPONENT OF SINGLE REPRESENTATIO
C   WHICH IS NOT APPERED IN THE PREVIOUS PROJECTION FOR THE
C   FIRST COMPONENT OF DUBLE REPRESENTATION
C
      DO 83 IB=1,NDS
      DO 831 IBS=1,2
         IF(INDA(IB,IBS).EQ.1) GO TO 831
            NB=IB
            NBS=IBS
C---------------------------------------------------
C  MATRIX ELEMENTS OF SINGLE REPRESENTATION
      CALL GETMXE(IRS,NDS,NH,NB,WWD)
C----------------------------------------------------
      DO 84 N=1,ND(JRD)
C----------------------------------------------------
C  MATRIX ELEMENT OF DOUBLE REPRESENTATION
C  SECOND SUFFIX OF PROJECTION OPERATOR IS FIXED AT NBD
      NNN=N
      CALL TSIRMR(JRD,NNN,1,MMG,NDND,JJG,JGA,CE,WD)
C-----------------------------------------------------
C     WRITE(6,*) MMG,(WD(NA),NA=1,MMG)
      DO 31 NA=1,NDS
      DO 32 NAS=1,2
         WA(NA,NAS)=0.0
         INDB(NA,NAS)=0
   32 CONTINUE
   31 CONTINUE
      DO 33 J=1,NH
      IF(ABS(WD(J)).GT.1.0D-7) THEN
         IF(IL.GE.1) THEN
            JJ=MOD(IG(JG(J))-1,24)+1
         ELSE
            JJ=MOD(IG(JG(J))-1,12)+1
         END IF
         DO 34 NA=1,NDS
         IF(ABS(WWD(J,NA)).LT.1.0D-7) GO TO 34
         DO 35 NAS=1,2
         IF(ABS(SN(NAS,NBS,JJ)).LT.1.0D-7) GO TO 35
         IF(N.EQ.1) INDB(NA,NAS)=1
             WA(NA,NAS)=WA(NA,NAS)
     &           +DCONJG(WD(J))*WWD(J,NA)*SN(NAS,NBS,JJ)
   35    CONTINUE
   34    CONTINUE
      END IF
   33 CONTINUE
      WC=0.0
      DO 36 NA=1,NDS
      DO 36 NAS=1,2
         WC=WC+DCONJG(WA(NA,NAS))*WA(NA,NAS)
   36 CONTINUE
C     WRITE(6,*) 'NB=',NB,' NBS=',NBS,' N=',N,' NNI=',NNI,' WC=',WC
      IF(WC.LT.1.0D-5) THEN
          IF(N.EQ.1) GO TO 831
          IF(N.NE.1) WRITE(6,*) ' THIS CAN NOT BE HAPPEN '
          STOP ' STOP AT 36 IN TSRDSD '
      END IF
      DO 92 NA=1,NDS
      DO 92 NAS=1,2
         IF(INDB(NA,NAS).EQ.1) INDA(NA,NAS)=1
   92 CONTINUE
      WC=DSQRT(WC)
      DO 37 NA=1,NDS
      DO 37 NAS=1,2
         DB(NA,NAS,N,NNI)=WA(NA,NAS)/WC
   37 CONTINUE
   84 CONTINUE
      GO TO 82
  831 CONTINUE
   83 CONTINUE
      WRITE(6,*) ' PIVOT SET ALGORISM IS NOT ENUGH '
      STOP ' STOP AT 83 IN TSRDSD '
   82 CONTINUE
   80 CONTINUE
      RETURN
      END
C SUBROUTINE GETMXE ====*====3====*====4====*====5====*====6====*====7
C
C  SLAVE ROUTINE FOR TSRDSD
C  NA=1->NDS,NB ELEMENT OF IRREDUCIBLE REPRESENTATION
C  IS GIVEN IN WWD(48,6)
C  IRS        :INTEGER FORM MATRIXRLEMENT
C  NDS        :DIMENSION OF REPRESENTATION IRS
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE GETMXE(IRS,NDS,NH,NB,WWD)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 WWD
      DIMENSION WWD(48,6),IRS(7,48)
      DO 41 NA=1,NDS
      IWA=32**(6-NB)
      DO 42 I=1,NH
      IW=IRS(NA,I)/IWA
      IW=IW-(IW/32)*32
      IF(IW.EQ.0) THEN
        WWD(I,NA)=0.0D0
      ELSE
        X=(DBLE(IW)/24.0D0)*2.0D0*3.14159265358979312D0
        WWD(I,NA)=DCMPLX(COS(X),SIN(X))/SQRT(DBLE(IRS(7,I)))
      END IF
   42 CONTINUE
   41 CONTINUE
      RETURN
      END
C SUBROUTINE TSLCLS ====*====3====*====4====*====5====*====6====*====7
C
C        SYMMETRIZED COEFFICIENTS FOR LOCAL SYMMETRY STATES(JR)
C        AT ATOMS WITH KINDS NUMBER IA.
C        NA COMPONENT OF I.R IIR IS OBTAINED IN THIS SUBROUTINE
C                 BY A.YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSLCLS(IIR,NA,JR,IA,KP,U,INS,A,ND1,ND2,NND)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 W,WB,V,U,WW,WD,CR
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IIGG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &     ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/ATT/ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &   ,NKATOM,NATOM,MKS(11,LMNKAT),IIAA(LMNATM),JRCH
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/ATT/
     &    ,NDIMC,NDIMH,ISIGMC,ISIGMH,IEMATC,IEMATH
      DIMENSION IND(150),MSIG(48),MIEM(48),MIA(150),MIT(150),
     &      V(150),MLL(3,LMNATM),NIQ(150),NJQ(150)
      DIMENSION KP(2,1),U(1),INS(3,1),A(1),NIP(150),WD(48)
      INTEGER NDIMC(10),NDIMH(12),ISIGMC(4,10),ISIGMH(12,6),IEMATC(24)
     &       ,IEMATH(12,2)
      REAL*8 EMAT(2,2,6)
      DATA NDIMC/
     &   1,1,1,1,2,2,3,3,3,3/
      DATA NDIMH/
     &   1,1,1,1,1,1,1,1,2,2,2,2/
      DATA ISIGMC/
     &   1,1,1,1, 1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1,
     &   1,1,1,1, 1,1,-1,-1,
     &   1,1,-1,-1, 1,1,1,1, 1,-1,-1,1, 1,-1,1,-1/
      DATA ISIGMH/
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &   1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1,-1,-1,-1, 1, 1, 1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
      DATA IEMATC/
     &   1,1,1,1,2,2,2,2,3,3,3,3,
     &   4,4,5,6,5,6,6,5,4,6,5,4/
      DATA IEMATH/
     &   1,3,2,1,3,2,4,6,5,4,6,5,
     &   1,2,3,1,2,3,4,5,6,4,5,6/
C    &  /1.0,0.0,0.0,1.0,
C    &   -0.5,0.8660254,-0.8660254,-0.5,
C    &   -0.5,-0.8660254,0.8660254,-0.5,
C    &   1.0,0.0,0.0,-1.0,
C    &   -0.5,0.8660254,0.8660254,0.5,
C    &   -0.5,-0.8660254,-0.8660254,0.5/
C
      IF(JRCH.EQ.1) CALL CHKNIR(JR,10)
      IF(JRCH.EQ.2) CALL CHKNIR(JR,12)
      CALL CHKNKA(IA,NKATOM)
      CALL CHKNIR(IIR,NR)
      CALL EMAT0(EMAT)
      PAI2=8.D0*DATAN(1.D0)
      NND=ND(IIR)
      IF(JRCH.EQ.2) NDIM=NDIMH(JR)
      IF(JRCH.EQ.1) NDIM=NDIMC(JR)
      IAS=1
      IF(IA.NE.1) IAS=KION(IA-1)+1
      IAL=KION(IA)
      NSI=IAL-IAS+1
      NNN=NDIM*NSI
      L=0
      DO 10 IAA=IAS,IAL
      DO 10 I=1,NDIM
      L=L+1
      MLL(I,IAA)=L
      MIA(L)=IAA
      MIT(L)=I
   10 CONTINUE
C
C     CHARACTER TEST
C
      W=0.0D0
      DO 2 IG=1,MG
      J=IIGG(JG(IG))
      IF(JRCH.EQ.2) ISG=ISIGMH(MOD(J-1,12)+1,(JR+1)/2)
      IF(JRCH.EQ.2.AND.J.GE.13.AND.MOD(JR,2).EQ.0) ISG=-ISG
      IF(JRCH.EQ.1) ISG=ISIGMC(((J-1)/12)+1,JR)
      MSIG(IG)=ISG
      IF(NDIM.NE.2) GO TO 27
      IF(JRCH.EQ.2) IEM=IEMATH(MOD(J-1,12)+1,(JR-9)/2+1)
      IF(JRCH.EQ.1) IEM=IEMATC(MOD(J-1,24)+1)
      MIEM(IG)=IEM
   27 WB=0.0D0
      DO 23 IAA=IAS,IAL
      IF(ISITR(IAA,J).NE.IAA) GO TO 23
      DO 24 I=1,NDIM
      SIG=ISG
      IF(NDIM.EQ.1) GO TO 25
      IF(NDIM.EQ.2) GO TO 26
      MA=IT(I,IV(J))
      IF(IABS(MA).NE.I) GO TO 24
      IF(MA.LT.0) SIG=-SIG
      GO TO 25
   26 SIG=SIG*EMAT(I,I,IEM)
   25 WA=0.0D0
      DO 19 K=1,3
      KK=JK(K,IG)
   19 WA=WA+VATOM(K,IAA)*KK
      WA=WA*PAI2
      WB=WB+DCMPLX(COS(WA),SIN(WA))*SIG
   24 CONTINUE
   23 CONTINUE
      W=W+CONJG(CR(IG,IIR))*WB
    2 CONTINUE
      W=W/MG
      NI=W+0.5
      NNII=NI
      ND1=0
      ND2=0
      IW=0
      IF(NI.EQ.0) GO TO 1
      NB=NA
      NP=0
      NQ=0
      INDIC=0
      INDIQ=0
      CALL ZZZY45(IIR)
      CALL TSIRME(NA,NA,WD)
      DO 3 I=1,NNN
    3 IND(I)=0
      IN=1
C
C     PROJECTION OPERATOR
C
   42 INT=MIT(IN)
      INA=MIA(IN)
      DO  8 J=1,NNN
    8 V(J)=0.0D0
      DO 9 IG=1,MG
      IF(ABS(WD(IG)).LT.1.D-4) GO TO 9
      WA=0.0D0
      DO 20 K=1,3
      KK=JK(K,IG)
   20 WA=WA+VATOM(K,INA)*DBLE(KK)
      X=PAI2*WA
      W=DCMPLX(COS(X),SIN(X))*CONJG(WD(IG))
      J=IIGG(JG(IG))
      IAA=ISITR(INA,J)
      SIG=MSIG(IG)
      IF(NDIM.EQ.1) GO TO 50
      IF(NDIM.EQ.2) GO TO 51
      MA=IT(INT,IV(J))
      MB=IABS(MA)
      IF(MA.LT.0) SIG=-SIG
      JJL=MLL(MB,IAA)
      V(JJL)=V(JJL)+W*SIG
      IND(JJL)=IN
      GO TO 9
   51 JJL=MLL(1,IAA)
      V(JJL)=V(JJL)+W*SIG*EMAT(1,INT,MIEM(IG))
      V(JJL+1)=V(JJL+1)+W*SIG*EMAT(2,INT,MIEM(IG))
      IND(JJL)=IN
      IND(JJL+1)=IN
      GO TO 9
   50 JJL=MLL(1,IAA)
      V(JJL)=V(JJL)+W*SIG
      IND(JJL)=IN
    9 CONTINUE
C
C     ORTHOGONALITY TEST
C
      IF(ND1.EQ.0) GO TO 31
      DO 32 J=1,ND1
      ID1=INS(1,J)
      ID2=INS(2,J)
      WW=0.0D0
      DO 33 I=ID1,ID2
      JU=MLL(KP(2,I),KP(1,I))
      WW=WW+CONJG(U(I))*V(JU)
   33 CONTINUE
      IF(ABS(WW).GT.1.D-4) GO TO 40
   32 CONTINUE
C
C     REGISTRATION
C
   31 ND1=ND1+1
      INS(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW+1
      DO 36 J=1,NNN
      IF(ABS(V(J)).GE.1.D-4) GO TO 362
      IF(J.EQ.IN) GO TO 362
      IF(IND(J).EQ.IN) IND(J)=0
      GO TO 36
  362 IW=IW+1
      KP(1,IW)=MIA(J)
      KP(2,IW)=MIT(J)
      U(IW)=V(J)*(DBLE(ND(IIR))/DBLE(MG))
      IF(ABS(DIMAG(U(IW))).LT.1.D-4) U(IW)=DBLE(U(IW))
      IF(ABS(DBLE(U(IW))).LT.1.D-4) U(IW)=DCMPLX(0.D0,DIMAG(U(IW)))
      WA=WA+DBLE(U(IW)*CONJG(U(IW)))
      IF(J.NE.IN) GO TO 36
      IF(IW.EQ.IIWW) GO TO 36
      DO 361 K=1,2
      KKW=KP(K,IW)
      KP(K,IW)=KP(K,IIWW)
  361 KP(K,IIWW)=KKW
      W=U(IW)
      U(IW)=U(IIWW)
      U(IIWW)=W
   36 CONTINUE
      INS(2,ND1)=IW
      A(ND1)=SQRT(WA)
      IF(WA.GT.1.D-4) GO TO 37
      ND1=ND1-1
      IW=IIWW-1
      IF(INDIC.NE.0) GO TO 40
      NQ=NQ+1
      NIQ(NQ)=IN
      IF(MOD(IN,2).EQ.0) NJQ(NQ)=IN-1
      IF(MOD(IN,2).EQ.1) NJQ(NQ)=IN+1
      GO TO 40
   37 IF(INDIC.EQ.1) GO TO 38
      IF(INDIC.EQ.2) GO TO 39
      IF(INDIC.EQ.3) GO TO 66
      NP=NP+1
      NIP(NP)=IN
   66 NI=NI-1
      INS(3,ND1)=0
      IF(NI.LE.0) GO TO 1
      IF(INDIC.EQ.3) GO TO 65
      GO TO 41
   38 INS(3,ND1)=ND1+1
      INDIC=2
      CALL TSIRME(NB,NB,WD)
      ND2=ND2+1
      GO TO 42
   39 INS(3,ND1)=-(ND1-1)
      NI=NI-1
      IF(NI.LE.0) GO TO 1
      GO TO 43
C
C
C
   40 IF(INDIC.EQ.2) GO TO 52
      IF(INDIC.EQ.3) GO TO 65
      IF(NA.NE.NB) GO TO 43
   41 JJWW=IN+1
      IF(JJWW.GT.NNN) GO TO 450
      DO 45 J=JJWW,NNN
      IF(IND(J).NE.0) GO TO 45
      IN=J
      GO TO 42
   45 CONTINUE
  450 CONTINUE
      NNP=NP
   64 NM=ND(IIR)
      IF(NM.NE.2.AND.NM.NE.4) GO TO 46
      IF(NM.EQ.2) NB=3-NA
      IF(NM.EQ.4) NB=5-NA
      NP=0
   43 NP=NP+1
      IF(NP.GT.NNP) GO TO 61
      IN=NIP(NP)
      INDIC=1
      CALL TSIRME(NA,NB,WD)
      GO TO 42
   61 IF(INDIQ.EQ.1) GO TO 46
      INDIC=3
      NNP=NQ
      DO 62 NQ=1,NNP
   62 NIP(NQ)=NIQ(NQ)
      NB=NA
      CALL TSIRME(NA,NA,WD)
      ND1=0
      ND2=0
      NI=NNII
      IW=0
      NQ=0
   65 NQ=NQ+1
      IF(NQ.GT.NNP) GO TO 63
      IN=NJQ(NQ)
      GO TO 42
   63 INDIC=1
      INDIQ=1
      GO TO 64
   46 WRITE(6,604)
  604 FORMAT(' STOP AT 46 IN TSLCLS')
      STOP
   52 WRITE(6,605)
  605 FORMAT(' STOP AT 52 IN TSLCLS')
      STOP
    1 NW=IW
      RETURN
      ENTRY TSLSDS
      WRITE(6,663) IIR,NA,JR,IA,NNN,ND1,ND2,ND1-ND2
  663 FORMAT(/8I5)
      IF(NW.EQ.0) GO TO 47
      ID1=0
      DO 48 IW=1,NW
      IF(IW.NE.INS(1,ID1+1).OR.ID1.GE.ND1) GO TO 49
      ID1=ID1+1
      WRITE(6,661) ID1,(INS(K,ID1),K=1,3),A(ID1)
  661 FORMAT(I4,3I5,F10.5)
      ICOUNT=0
   49 ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KP(K,IW),K=1,2),U(IW)
  660 FORMAT(4X,2I5,2H (,2I5,4H)  (,2F10.5,1H))
   48 CONTINUE
   47 RETURN
      END
C SUBROUTINE TSLCLA ====*====3====*====4====*====5====*====6====*====7
C
C        SYMMETRIZED COEFFICIENTS FOR LOCAL SYMMETRY STATES(JR)
C        AT ATOMS WITH KINDS NUMBER IA.
C        ALL COMPONENT OF I.R IIR IS OBTAINED IN THIS SUBROUTINE
C                 BY A.YANASE
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSLCLA(IIR,JR,IA,KP,U,INS,A,ND1,NND)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 W,WB,V,U,WW,WD,CR
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IIGG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &     ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/ATT/ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &    ,NKATOM,NATOM,MKS(11,LMNKAT),IIAA(LMNATM),JRCH
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/ATT/
     &  ,NDIMC,NDIMH,ISIGMC,ISIGMH,IEMATC,IEMATH
      DIMENSION IND(150),MSIG(48),MIEM(48),MIA(150),MIT(150),
     &      V(150),MLL(3,LMNATM),NIQ(150),NJQ(150)
      DIMENSION KP(2,1),U(1),INS(4,1),A(1),NIP(150),WD(48)
      INTEGER NDIMC(10),NDIMH(12),ISIGMC(4,10),ISIGMH(12,6)
     &       ,IEMATC(24),IEMATH(12,2)
      REAL*8 EMAT(2,2,6)      
      DATA NDIMC/
     &   1,1,1,1,2,2,3,3,3,3/
      DATA NDIMH/
     &   1,1,1,1,1,1,1,1,2,2,2,2/
      DATA ISIGMC/
     &   1,1,1,1, 1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1,
     &   1,1,1,1, 1,1,-1,-1,
     &   1,1,-1,-1, 1,1,1,1, 1,-1,-1,1, 1,-1,1,-1/
      DATA ISIGMH/
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &   1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1,-1,-1,-1, 1, 1, 1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
      DATA IEMATC/
     &   1,1,1,1,2,2,2,2,3,3,3,3,
     &   4,4,5,6,5,6,6,5,4,6,5,4/
      DATA IEMATH/
     &   1,3,2,1,3,2,4,6,5,4,6,5,
     &   1,2,3,1,2,3,4,5,6,4,5,6/
C    &  /1.0,0.0,0.0,1.0,
C    &   -0.5,0.8660254,-0.8660254,-0.5,
C    &   -0.5,-0.8660254,0.8660254,-0.5,
C    &   1.0,0.0,0.0,-1.0,
C    &   -0.5,0.8660254,0.8660254,0.5,
C    &   -0.5,-0.8660254,-0.8660254,0.5/
C
      IF(JRCH.EQ.1) CALL CHKNIR(JR,10)
      IF(JRCH.EQ.2) CALL CHKNIR(JR,12)
      CALL CHKNKA(IA,NKATOM)
      CALL CHKNIR(IIR,NR)
      CALL EMAT0(EMAT)
      PAI2=8.D0*DATAN(1.D0)
      NND=ND(IIR)
      IF(JRCH.EQ.2) NDIM=NDIMH(JR)
      IF(JRCH.EQ.1) NDIM=NDIMC(JR)
      IAS=1
      IF(IA.NE.1) IAS=KION(IA-1)+1
      IAL=KION(IA)
      NSI=IAL-IAS+1
      NNN=NDIM*NSI
      L=0
      DO 10 IAA=IAS,IAL
      DO 10 I=1,NDIM
      L=L+1
      MLL(I,IAA)=L
      MIA(L)=IAA
      MIT(L)=I
   10 CONTINUE
C
C     CHARACTER TEST
C
      W=0.0D0
      DO 2 IG=1,MG
      J=IIGG(JG(IG))
      IF(JRCH.EQ.2) ISG=ISIGMH(MOD(J-1,12)+1,(JR+1)/2)
      IF(JRCH.EQ.2.AND.J.GE.13.AND.MOD(JR,2).EQ.0) ISG=-ISG
      IF(JRCH.EQ.1) ISG=ISIGMC(((J-1)/12)+1,JR)
      MSIG(IG)=ISG
      IF(NDIM.NE.2) GO TO 27
      IF(JRCH.EQ.2) IEM=IEMATH(MOD(J-1,12)+1,(JR-9)/2+1)
      IF(JRCH.EQ.1) IEM=IEMATC(MOD(J-1,24)+1)
      MIEM(IG)=IEM
   27 WB=0.0D0
      DO 23 IAA=IAS,IAL
      IF(ISITR(IAA,J).NE.IAA) GO TO 23
      DO 24 I=1,NDIM
      SIG=ISG
      IF(NDIM.EQ.1) GO TO 25
      IF(NDIM.EQ.2) GO TO 26
      MA=IT(I,IV(J))
      IF(IABS(MA).NE.I) GO TO 24
      IF(MA.LT.0) SIG=-SIG
      GO TO 25
   26 SIG=SIG*EMAT(I,I,IEM)
   25 WA=0.0D0
      DO 19 K=1,3
      KK=JK(K,IG)
   19 WA=WA+VATOM(K,IAA)*KK
      WA=WA*PAI2
      WB=WB+DCMPLX(COS(WA),SIN(WA))*SIG
   24 CONTINUE
   23 CONTINUE
      W=W+CONJG(CR(IG,IIR))*WB
    2 CONTINUE
      W=W/MG
      NI=W+0.5
      NNII=NI
      NA=1
      NC=1
      ND1=0
      IW=0
      IF(NI.EQ.0) GO TO 1
      NB=NA
      NP=0
      NQ=0
      INDIC=0
      INDIQ=0
      CALL ZZZY45(IIR)
      CALL TSIRME(1,NA,WD)
      DO 3 I=1,NNN
    3 IND(I)=0
      IN=1
C
C     PROJECTION OPERATOR
C
   42 INT=MIT(IN)
      INA=MIA(IN)
      DO  8 J=1,NNN
    8 V(J)=0.0D0
      DO 9 IG=1,MG
      IF(ABS(WD(IG)).LT.1.D-4) GO TO 9
      WA=0.0D0
      DO 20 K=1,3
      KK=JK(K,IG)
   20 WA=WA+VATOM(K,INA)*DBLE(KK)
      X=PAI2*WA
      W=DCMPLX(COS(X),SIN(X))*CONJG(WD(IG))
      J=IIGG(JG(IG))
      IAA=ISITR(INA,J)
      SIG=MSIG(IG)
      IF(NDIM.EQ.1) GO TO 50
      IF(NDIM.EQ.2) GO TO 51
      MA=IT(INT,IV(J))
      MB=IABS(MA)
      IF(MA.LT.0) SIG=-SIG
      JJL=MLL(MB,IAA)
      V(JJL)=V(JJL)+W*SIG
      IF(NC.EQ.1) IND(JJL)=IN
      GO TO 9
   51 JJL=MLL(1,IAA)
      V(JJL)=V(JJL)+W*SIG*EMAT(1,INT,MIEM(IG))
      V(JJL+1)=V(JJL+1)+W*SIG*EMAT(2,INT,MIEM(IG))
      IF(NC.EQ.1) IND(JJL)=IN
      IF(NC.EQ.1) IND(JJL+1)=IN
      GO TO 9
   50 JJL=MLL(1,IAA)
      V(JJL)=V(JJL)+W*SIG
      IF(NC.EQ.1) IND(JJL)=IN
    9 CONTINUE
C
C     ORTHOGONALITY TEST
C
      IF(ND1.EQ.0) GO TO 31
      DO 32 J=1,ND1
      ID1=INS(1,J)
      ID2=INS(2,J)
      WW=0.0D0
      DO 33 I=ID1,ID2
      JU=MLL(KP(2,I),KP(1,I))
      WW=WW+CONJG(U(I))*V(JU)
   33 CONTINUE
      IF(ABS(WW).GT.1.D-4) GO TO 40
   32 CONTINUE
C
C     REGISTRATION
C
   31 ND1=ND1+1
      INS(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW+1
      DO 36 J=1,NNN
      IF(ABS(V(J)).GE.1.D-4) GO TO 362
      IF(J.EQ.IN) GO TO 362
      IF(IND(J).EQ.IN.AND.NC.EQ.1) IND(J)=0
      GO TO 36
  362 IW=IW+1
      KP(1,IW)=MIA(J)
      KP(2,IW)=MIT(J)
      U(IW)=V(J)*(DBLE(ND(IIR))/DBLE(MG))
      IF(ABS(DIMAG(U(IW))).LT.1.D-4) U(IW)=DBLE(U(IW))
      IF(ABS(DBLE(U(IW))).LT.1.D-4) U(IW)=DCMPLX(0.D0,DIMAG(U(IW)))
      WA=WA+DBLE(U(IW)*CONJG(U(IW)))
      IF(J.NE.IN) GO TO 36
      IF(IW.EQ.IIWW) GO TO 36
      DO 361 K=1,2
      KKW=KP(K,IW)
      KP(K,IW)=KP(K,IIWW)
  361 KP(K,IIWW)=KKW
      W=U(IW)
      U(IW)=U(IIWW)
      U(IIWW)=W
   36 CONTINUE
      INS(2,ND1)=IW
      A(ND1)=SQRT(WA)
      IF(WA.GT.1.D-4) GO TO 37
      ND1=ND1-1
      IW=IIWW-1
      IF(INDIC.NE.0) GO TO 40
      NQ=NQ+1
      NIQ(NQ)=IN
      IF(MOD(IN,2).EQ.0) NJQ(NQ)=IN-1
      IF(MOD(IN,2).EQ.1) NJQ(NQ)=IN+1
      GO TO 40
   37 INS(3,ND1)=NC
      INS(4,ND1)=NB
      IF(NC.GE.ND(IIR)) GO TO 38
      NC=NC+1
      CALL TSIRME(NC,NB,WD)
      GO TO 42
   38 NI=NI-1
      IF(NI.LE.0) GO TO 1
      NC=1
      IF(INDIC.EQ.1) GO TO 43
      IF(INDIC.EQ.3) GO TO 65
      NP=NP+1
      NIP(NP)=IN
      GO TO 41
C
C
C
   40 IF(NC.NE.1) GO TO 52
      NC=1
      IF(INDIC.EQ.3) GO TO 65
      IF(INDIC.EQ.1) GO TO 43
   41 JJWW=IN+1
      IF(JJWW.GT.NNN) GO TO 450
      DO 45 J=JJWW,NNN
      IF(IND(J).NE.0) GO TO 45
      IN=J
      CALL TSIRME(1,NA,WD)
      GO TO 42
   45 CONTINUE
  450 CONTINUE
      NNP=NP
   64 NM=ND(IIR)
      IF(NM.NE.2.AND.NM.NE.4) GO TO 46
      IF(NM.EQ.2) NB=3-NA
      IF(NM.EQ.4) NB=5-NA
      NP=0
   43 NP=NP+1
      IF(NP.GT.NNP) GO TO 61
      IN=NIP(NP)
      INDIC=1
      CALL TSIRME(1,NB,WD)
      GO TO 42
   61 IF(INDIQ.EQ.1) GO TO 46
      INDIC=3
      NNP=NQ
      DO 62 NQ=1,NNP
   62 NIP(NQ)=NIQ(NQ)
      NB=NA
      ND1=0
      NI=NNII
      IW=0
      NQ=0
   65 NQ=NQ+1
      IF(NQ.GT.NNP) GO TO 63
      IN=NJQ(NQ)
      CALL TSIRME(1,NA,WD)
      GO TO 42
   63 INDIC=1
      INDIQ=1
      GO TO 64
   46 WRITE(6,604)
  604 FORMAT(' STOP AT 46 IN TSLCLA')
      STOP
   52 WRITE(6,605)
  605 FORMAT(' STOP AT 52 IN TSLCLA')
      STOP
    1 NW=IW
      RETURN
      ENTRY TSLADS
      WRITE(6,663) IIR,NA,JR,IA,NNN,ND1,ND1/ND(IIR)
  663 FORMAT(/8I5)
      IF(NW.EQ.0) GO TO 47
      ID1=0
      DO 48 IW=1,NW
      IF(IW.NE.INS(1,ID1+1).OR.ID1.GE.ND1) GO TO 49
      ID1=ID1+1
      WRITE(6,661) ID1,(INS(K,ID1),K=1,4),A(ID1)
  661 FORMAT(I4,4I5,F10.5)
      ICOUNT=0
   49 ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KP(K,IW),K=1,2),U(IW)
  660 FORMAT(4X,2I5,2H (,2I5,4H)  (,2F10.5,1H))
   48 CONTINUE
   47 RETURN
      END
C SUBROUTINE ZZZY46 ====*====3====*====4====*====5====*====6====*====7
C
C     SPHERICAL HARMONICS FOR(L=0,1,2 AND 3) ARE TRANSFORMED TO
C     CUBIC(IJ.EQ.1) OR HEXAGONAL(IJ.EQ.2) HARMONICS BY TSTRLM
C     16 BASIS ARE LISTED IN THE MANUAL OR IN THE DATA OF LM,JC,
C     JH,NC AND NH
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE ZZZY46(IJ)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 TT,U
      COMMON/SCC/TT(16,7)
      SAVE /SCC/,LM,JC,JH,NC,NH
      DIMENSION U(150),KP(150),INS(4,21)
      INTEGER LM(16),JC(16),JH(16),NC(16),NH(16)
      DATA LM/0,1, 1, 1,2,2,2, 2, 2,3, 3, 3, 3, 3, 3, 3/
      DATA JC/1,8, 8, 8,5,5,9, 9, 9,4, 8, 8, 8,10,10,10/
      DATA JH/1,4,10,10,1,9,9,11,11,4,10,10,12,12, 6, 8/
      DATA NC/1,1, 2, 3,1,2,1, 2, 3,1, 1, 2, 3, 1, 2, 3/
      DATA NH/1,1, 1, 2,1,1,2, 1, 2,1, 1, 2, 1, 2, 1, 1/
      DO 1 LL=1,4
      L=LL-1
      JRCH=IJ
      CALL TSTRLM(L,U,KP,INS)
C     CALL TSLMDS
      DO 2 J=1,16
      IF(LM(J).NE.L) GO TO 2
      DO 3 IS=1,2*L+1
      IF(JRCH.EQ.1.AND.JC(J).NE.INS(4,IS)) GO TO 3
      IF(JRCH.EQ.1.AND.NC(J).NE.INS(3,IS)) GO TO 3
      IF(JRCH.EQ.2.AND.JH(J).NE.INS(4,IS)) GO TO 3
      IF(JRCH.EQ.2.AND.NH(J).NE.INS(3,IS)) GO TO 3
      DO 4 I=1,7
    4 TT(J,I)=0.0D0
      DO 5 I=INS(1,IS),INS(2,IS)
      TT(J,KP(I))=U(I)
    5 CONTINUE
      GO TO 2
    3 CONTINUE
      WRITE(6,600) J
  600 FORMAT(' STOP AT 3 IN ZZZY46 FOR',I2)
      STOP
    2 CONTINUE
    1 CONTINUE
C     DO 60 J=1,16
C     WRITE(6,660) (TT(J,I),I=1,7)
C 660 FORMAT(8F9.5)
C  60 CONTINUE
      RETURN
      END
C SUBROUTINE TSTRLM ====*====3====*====4====*====5====*====6====*====7
C
C    SPHERICAL HARMONICS(L) ARE TRANSFORMED TO CUBIC(IL.GE.1)
C    OR HEXAGONAL(IL.LE.0) HARMINICS
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSTRLM(L,U,KP,INS)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NGNG,IIGG(48),JV(2,3,48)
      SAVE /SPG2/,RC,RH,TC,TH,NDIMC,NDIMH
      COMPLEX*16 U(150),WD(21),WE(21),CW,WW
      DIMENSION IC(16),WM(3),IM(3),W(3),KP(150),INS(4,21)
      INTEGER RC(3,24),RH(3,12),TC(24),TH(12),NDIMC(10),NDIMH(12)
      DATA RC/
     &  0,0,0, 0,2,2, 0,2,0, 0,0,2,
     &  0,1,1, 2,1,3, 2,1,1, 0,1,3,
     &  1,1,2, 3,1,0, 1,1,0, 3,1,2,
     &  0,2,1, 0,2,3, 0,1,2, 1,1,1, 2,1,0, 3,1,3,
     &  3,1,1, 0,1,0, 1,0,0, 1,1,3, 2,1,2, 3,0,0/
      DATA RH/
     &  0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0,
     &  0,2,0, 4,2,0, 2,2,0, 3,2,0, 1,2,0, 5,2,0/
      DATA TC/
     &  1,2,2,2,3,3,3,3,3,3,3,3,
     &  2,2,2,2,2,2,4,4,4,4,4,4/
      DATA TH/
     &  1,6,3,2,3,6,2,2,2,2,2,2/
      DATA NDIMC/
     &    1,1,1,1,2,2,3,3,3,3/
      DATA NDIMH/
     &    1,1,1,1,1,1,1,1,2,2,2,2/
      PAI=4.D0*DATAN(1.D0)
      ND1=0
      IW=0
      JRCH=1
      IF(IL.LE.0) JRCH=2
      NG=24
      IF(JRCH.EQ.2) NG=12
      IPR=0
      NJR=9
      IF(JRCH.EQ.2) NJR=11
      DO 801 JR=1,NJR,2
      NDIM=NDIMC(JR)
      IF(JRCH.EQ.2) NDIM=NDIMH(JR)
C
C   CHARCTER TEST
C
      WB=0.0D0
      DO 204 IG=1,NG
      CL=2*L+1
      IF(IG.EQ.1) GO TO 201
      IF(JRCH.EQ.1) XX=PAI/TC(IG)
      IF(JRCH.EQ.2) XX=PAI/TH(IG)
      CL=SIN((2*L+1)*XX)/SIN(XX)
  201 CONTINUE
C     WRITE(6,693) L,IG,CL
C 693 FORMAT(2I5,F8.4)
      IF(ABS(CL).LT.0.001D0) GO TO 204
      J=IG
      WA=0.0D0
      DO 205 IN=1,NDIM
      I=IN
      JRJR=JR
      CALL TSTRLS(JRJR,I,J,WM,JRCH)
      WA=WA+WM(IN)
  205 CONTINUE
C     WRITE(6,691) WA,CL
C 691 FORMAT(8F10.4)
      IF(ABS(WA).LT.0.001D0) GO TO 204
      WB=WB+CL*WA
  204 CONTINUE
      WB=WB/NG
      NI=WB+0.5D0
C     WRITE(6,690) L,JR,NI
C 690 FORMAT(16I5)
      IF(NI.EQ.0) GO TO 801
C
C  PROJECTION OPERATOR
C
      NC=1
      NB=1
      MB=L
  150 CONTINUE
C     WRITE(6,695) L,JR,NI,MB,NC,NB
C 695 FORMAT(6I4)
      DO 101 IA=1,2*L+1
      MA=L-IA+1
      CALL TSRMI(L,MA,MB,K1,K2,K3,IS,IC,IPR)
      IM(1)=0
      IM(2)=IC(1)
      IM(3)=0
      IF(IS.EQ.2) GO TO 102
      DO 103 I=1,16
      IM(1)=IM(1)+IC(I)
      IM(3)=IM(3)+IC(I)*(MOD(I,2)*2-1)
  103 CONTINUE
  102 WW=(DBLE(K1)/DBLE(K2))*SQRT(DBLE(K3))
      W(1)=WW*IM(1)
      W(2)=WW*IM(2)
      W(3)=WW*IM(3)
      WD(IA)=0.0D0
      DO 104 IG=1,NG
      IF(JRCH.EQ.1) I2=RC(2,IG)+1
      IF(JRCH.EQ.2) I2=RH(2,IG)+1
      IF(IM(I2).EQ.0) GO TO 104
      J=IG
      JJR=JR
      CALL TSTRLS(JJR,NB,J,WM,JRCH)
      IF(ABS(WM(NC)).LT.0.001) GO TO 104
      IF(JRCH.EQ.2) GO TO 105
      JW=RC(1,IG)*MA+RC(3,IG)*MB
      XX=DBLE(MOD(JW,4))*(PAI/2.0D0)
      GO TO 106
  105 JW=RH(1,IG)*MA
      XX=DBLE(MOD(JW,6))*(PAI/3.0D0)
  106 CW=DCMPLX(COS(XX),-SIN(XX))
      WD(IA)=WD(IA)+CW*W(I2)*WM(NC)
  104 CONTINUE
  101 CONTINUE
C
C    ORTHOGONALITY TEST
C
      IF(ND1.EQ.0) GO TO 111
      DO 112 ID=1,ND1
      ID1=INS(1,ID)
      ID2=INS(2,ID)
      WW=0.0D0
      DO 113 I=ID1,ID2
      JU=KP(I)
      WW=WW+CONJG(U(I))*WD(JU)
  113 CONTINUE
      IF(ABS(WW).LE.1.0D-4) GO TO 112
C
C    ORTHOGONALIZATION
C
      DO 121 I=INS(1,ID),INS(2,ID)
      JU=KP(I)
      WD(JU)=WD(JU)-U(I)*WW
  121 CONTINUE
C
C
  112 CONTINUE
C
C   REGISTRATION
C
  111 ND1=ND1+1
      INS(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW
      IMC=99
      DO 114 IA=1,2*L+1
      IF(ABS(WD(IA)).LT.1.0D-4) GO TO 114
      IW=IW+1
      IF(IMC.EQ.-(L-IA+1)) IMD=IW
      IF(IMC.EQ.99) IMC=L-IA+1
      KP(IW)=IA
      U(IW)=WD(IA)
      WA=WA+DBLE(U(IW)*CONJG(U(IW)))
  114 CONTINUE
      WA=SQRT(WA)
      IF(WA.GT.1.0D-4) GO TO 116
      IW=IIWW
      ND1=ND1-1
      GO TO 303
  116 INS(2,ND1)=IW
      INS(3,ND1)=NC
      INS(4,ND1)=JR
      IF(MOD(L,2).EQ.1) INS(4,ND1)=JR+1
      CW=1.0/WA
      IF(IMC.EQ.0) GO TO 115
      IF(MOD(IMC,2).EQ.1.AND.DBLE(U(IMD)).GT.0.D0) CW=(0.D0,-1.D0)/WA
      IF(MOD(IMC,2).EQ.0.AND.DBLE(U(IMD)).LT.0.D0) CW=(0.D0,-1.D0)/WA
  115 DO 117 I=INS(1,ND1),INS(2,ND1)
  117 U(I)=U(I)*CW
C
C    PARTNER
C
      IF(NDIM.EQ.1) GO TO 130
      DO 131 IA=1,2*L+1
      MA=L-IA+1
      WD(IA)=0
      WE(IA)=0
      DO 137 I=INS(1,ND1),INS(2,ND1)
      MBB=L-KP(I)+1
      CALL TSRMI(L,MA,MBB,K1,K2,K3,IS,IC,IPR)
      IM(1)=0
      IM(2)=IC(1)
      IM(3)=0
      IF(IS.EQ.2) GO TO 132
      DO 133 J=1,16
      IM(1)=IM(1)+IC(J)
      IM(3)=IM(3)+IC(J)*(MOD(J,2)*2-1)
  133 CONTINUE
  132 WW=(DBLE(K1)/DBLE(K2))*SQRT(DBLE(K3))
      W(1)=WW*IM(1)
      W(2)=WW*IM(2)
      W(3)=WW*IM(3)
      DO 134 IG=8,9
      IF(JRCH.EQ.1) I2=RC(2,IG)+1
      IF(JRCH.EQ.2) I2=RH(2,IG)+1
      IF(IM(I2).EQ.0) GO TO 134
      IF(JRCH.EQ.2) GO TO 135
      JW=RC(1,IG)*MA+RC(3,IG)*MBB
      XX=DBLE(MOD(JW,4))*(PAI/2.0)
      GO TO 136
  135 JW=RH(1,IG)*MA
      XX=DBLE(MOD(JW,6))*(PAI/3.0D0)
  136 CW=DCMPLX(COS(XX),-SIN(XX))
      IF(NDIM.EQ.2.AND.IG.EQ.8) WD(IA)=WD(IA)+CW*W(I2)*U(I)
      IF(NDIM.EQ.2.AND.IG.EQ.9) WD(IA)=WD(IA)-CW*W(I2)*U(I)
      IF(NDIM.EQ.3.AND.IG.EQ.9) WE(IA)=WE(IA)+CW*W(I2)*U(I)
      IF(NDIM.EQ.3.AND.IG.EQ.8) WD(IA)=WD(IA)-CW*W(I2)*U(I)
  134 CONTINUE
  137 CONTINUE
  131 CONTINUE
      DO 138 N=2,NDIM
      ND1=ND1+1
      INS(1,ND1)=IW+1
      DO 139 IA=1,2*L+1
      IF(N.EQ.2.AND.ABS(WD(IA)).LT.1.D-4) GO TO 139
      IF(N.EQ.3.AND.ABS(WE(IA)).LT.1.D-4) GO TO 139
      IW=IW+1
      KP(IW)=IA
      IF(NDIM.EQ.2) U(IW)=WD(IA)/SQRT(3.0D0)
      IF(NDIM.EQ.2.AND.JR.EQ.9) U(IW)=-U(IW)
      IF(NDIM.EQ.3.AND.N.EQ.2) U(IW)=WD(IA)
      IF(NDIM.EQ.3.AND.N.EQ.3) U(IW)=WE(IA)
  139 CONTINUE
      INS(2,ND1)=IW
      INS(3,ND1)=N
      INS(4,ND1)=JR
      IF(MOD(L,2).EQ.1) INS(4,ND1)=JR+1
  138 CONTINUE
  130 CONTINUE
      NI=NI-1
      IF(NI.EQ.0) GO TO 801
      GO TO 302
C
C    FIND A NEW PIVOT
C
  303 CONTINUE
C     WRITE(6,693)
C 693 FORMAT(' ZERO VECTOR')
  302 CONTINUE
      MB=MB-1
      IF(MB.GE.-L) GO TO 150
      WRITE(6,601)
  601 FORMAT(' STOP NEAR 301 IN TSTRLM')
      STOP
  801 CONTINUE
      NW=IW
      RETURN
      ENTRY TSLMDS
      WRITE(6,661) L,JRCH,NW
  661 FORMAT(' L=',I3,' JRCH=',I3,' NW=',I5)
      ID1=0
      DO 811 IW=1,NW
      IF(IW.NE.INS(1,ID1+1).OR.ID1.GE.ND1) GO TO 812
      ID1=ID1+1
      WRITE(6,662) ID1,(INS(K,ID1),K=1,4)
  662 FORMAT(I4,4I5)
      ICOT=0
  812 ICOT=ICOT+1
      WRITE(6,663) IW,ICOT,L-KP(IW)+1,U(IW)
  663 FORMAT(4X,2I5,I4,3H  (,2F10.5,1H))
  811 CONTINUE
      RETURN
      END
C SUBROUTINE TSTRLS ====*====3====*====4====*====5====*====6====*====7
C
C    CUBIC(JRCH=1) OR HEXAGONAL(JRCH=2) HARMONICS DIFFINEDD IN THE MANUAL
C    OR ZZZY46 ARE TRANSFORMED BY ROTATION J AND COEFFICENTS ARE GIVEN IN
C    WD
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSTRLS(JR,NR,J,WD,JRCH)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      SAVE /SPG1/
     &  ,NDIMC,NDIMH,ISIGMC,ISIGMH,IEMATC,IEMATH      
      INTEGER NDIMC(10),NDIMH(12),ISIGMC(4,10),ISIGMH(12,6)
     &       ,IEMATC(24),IEMATH(12,2)
      DIMENSION EMAT(2,2,6)
      DIMENSION WD(3)
      DATA NDIMC/
     &   1,1,1,1,2,2,3,3,3,3/
      DATA NDIMH/
     &   1,1,1,1,1,1,1,1,2,2,2,2/
      DATA ISIGMC/
     &   1,1,1,1, 1,1,-1,-1, 1,-1,1,-1, 1,-1,-1,1,
     &   1,1,1,1, 1,1,-1,-1,
     &   1,1,-1,-1, 1,1,1,1, 1,-1,-1,1, 1,-1,1,-1/
      DATA ISIGMH/
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     &   1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1,-1, 1,-1, 1,-1,-1,-1,-1, 1, 1, 1,
     &   1,-1, 1,-1, 1,-1, 1, 1, 1,-1,-1,-1,
     &   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
      DATA IEMATC/
     &   1,1,1,1,2,2,2,2,3,3,3,3,
     &   4,4,5,6,5,6,6,5,4,6,5,4/
      DATA IEMATH/
     &   1,3,2,1,3,2,4,6,5,4,6,5,
     &   1,2,3,1,2,3,4,5,6,4,5,6/
C    &  /1.0,0.0,0.0,1.0,
C    &   -0.5,0.8660254,-0.8660254,-0.5,
C    &   -0.5,-0.8660254,0.8660254,-0.5,
C    &   1.0,0.0,0.0,-1.0,
C    &   -0.5,0.8660254,0.8660254,0.5,
C    &   -0.5,-0.8660254,-0.8660254,0.5/
C
      CALL EMAT0(EMAT)
      IF(JRCH.EQ.2) NDIM=NDIMH(JR)
      IF(JRCH.EQ.1) NDIM=NDIMC(JR)
      IF(JRCH.EQ.2) SIG=ISIGMH(MOD(J-1,12)+1,(JR+1)/2)
      IF(JRCH.EQ.2.AND.J.GE.13.AND.MOD(JR,2).EQ.0) SIG=-SIG
      IF(JRCH.EQ.1) SIG=ISIGMC(((J-1)/12)+1,JR)
      IF(NDIM.EQ.1) GO TO 50
      IF(NDIM.EQ.3) GO TO 51
      IF(JRCH.EQ.2) IEM=IEMATH(MOD(J-1,12)+1,(JR-9)/2+1)
      IF(JRCH.EQ.1) IEM=IEMATC(MOD(J-1,24)+1)
      WD(1)=SIG*EMAT(1,NR,IEM)
      WD(2)=SIG*EMAT(2,NR,IEM)
      GO TO 59
   51 MA=IT(NR,IV(J))
      MB=IABS(MA)
      IF(MA.LT.0) SIG=-SIG
      DO 55 I=1,3
   55 WD(I)=0.0D0
      WD(MB)=SIG
      GO TO 59
   50 WD(1)=SIG
   59 CONTINUE
      RETURN
      END
C SUBROUTINE TSTRLB ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSTRLB(L,U,KP,INS)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NNGG,IIGG(48),JV(2,3,48)
      SAVE /SPG2/,RC,RH
      COMPLEX*16 WD(21,24),CW,U(150),WC
      DIMENSION KP(150),INS(4,21),IM(3),W(3),IC(16)
      INTEGER RC(3,24),RH(3,12)
      DATA RC/
     &  0,0,0, 0,2,2, 0,2,0, 0,0,2,
     &  0,1,1, 2,1,3, 2,1,1, 0,1,3,
     &  1,1,2, 3,1,0, 1,1,0, 3,1,2,
     &  0,2,1, 0,2,3, 0,1,2, 1,1,1, 2,1,0, 3,1,3,
     &  3,1,1, 0,1,0, 1,0,0, 1,1,3, 2,1,2, 3,0,0/
      DATA RH/
     &  0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0,
     &  0,2,0, 4,2,0, 2,2,0, 3,2,0, 1,2,0, 5,2,0/
      PAI=4.D0*DATAN(1.D0)
      JRCH=1
      IF(IL.LE.0) JRCH=2
      NG=24
      IF(JRCH.EQ.2) NG=12
      DO 130 ND1=1,2*L+1
      DO 131 IA=1,2*L+1
      MA=L-IA+1
      DO 151 IG=1,NG
  151 WD(IA,IG)=0.0D0
      DO 137 I=INS(1,ND1),INS(2,ND1)
      MBB=L-KP(I)+1
      CALL TSRMI(L,MA,MBB,K1,K2,K3,IS,IC,IPR)
      IM(1)=0
      IM(2)=IC(1)
      IM(3)=0
      IF(IS.EQ.2) GO TO 132
      DO 133 J=1,16
      IM(1)=IM(1)+IC(J)
      IM(3)=IM(3)+IC(J)*(MOD(J,2)*2-1)
  133 CONTINUE
  132 WW=(DBLE(K1)/DBLE(K2))*SQRT(DBLE(K3))
      W(1)=WW*IM(1)
      W(2)=WW*IM(2)
      W(3)=WW*IM(3)
      DO 134 IG=1,NG
      IF(JRCH.EQ.1) I2=RC(2,IG)+1
      IF(JRCH.EQ.2) I2=RH(2,IG)+1
      IF(IM(I2).EQ.0) GO TO 134
      IF(JRCH.EQ.2) GO TO 135
      JW=RC(1,IG)*MA+RC(3,IG)*MBB
      XX=DBLE(MOD(JW,4))*(PAI/2.0D0)
      GO TO 136
  135 JW=RH(1,IG)*MA
      XX=DBLE(MOD(JW,6))*(PAI/3.0D0)
  136 CW=DCMPLX(COS(XX),-SIN(XX))
      WD(IA,IG)=WD(IA,IG)+CW*W(I2)*U(I)
  134 CONTINUE
  137 CONTINUE
  131 CONTINUE
      DO 138 ND2=1,2*L+1
      DO 140 IG=1,NG
      WC=0.0D0
      DO 139 I=INS(1,ND2),INS(2,ND2)
      WC=WC+CONJG(U(I))*WD(KP(I),IG)
  139 CONTINUE
      IF(ABS(WC).GT.1.D-4) WRITE(6,660) ND2,ND1,IG,WC
  660 FORMAT(3I5,2F10.5)
  140 CONTINUE
  138 CONTINUE
  130 CONTINUE
      RETURN
      END
C SUBROUTINE EMAT0 2====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE EMAT0(EMAT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION EMAT(2,2,6)
C
      P2=DSQRT(3.0D0)/2.0D0
      EMAT(1,1,1) = 1.0D0
      EMAT(2,1,1) = 0.0D0
      EMAT(1,2,1) = 0.0D0
      EMAT(2,2,1) = 1.0D0
C
      EMAT(1,1,2) =-0.5D0
      EMAT(2,1,2) = P2
      EMAT(1,2,2) =-P2
      EMAT(2,2,2) =-0.5D0
C
      EMAT(1,1,3) =-0.5D0
      EMAT(2,1,3) =-P2
      EMAT(1,2,3) = P2
      EMAT(2,2,3) =-0.5D0
C
      EMAT(1,1,4) = 1.0D0
      EMAT(2,1,4) = 0.0D0
      EMAT(1,2,4) = 0.0D0
      EMAT(2,2,4) =-1.0D0
C
      EMAT(1,1,5) =-0.5D0
      EMAT(2,1,5) = P2
      EMAT(1,2,5) = P2
      EMAT(2,2,5) = 0.5D0
C
      EMAT(1,1,6) =-0.5D0
      EMAT(2,1,6) =-P2
      EMAT(1,2,6) =-P2
      EMAT(2,2,6) = 0.5D0
C
      RETURN
      END
C SUBROUTINE TSBASE ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSBASE(IR,NA)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 V,U,WWST
      COMMON/BON/KBO(2,30),BOL(30),NBO,KBOND(4,200),NBOND,WWST(4000)
     &       ,IVEC(10000),VST(3,4000),NVEC,NDIM
     &       ,KA(2,1500),V(1500),INS(3,500),A(500),NT,NS
      COMMON/ATT   /ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &      ,NKATOM,NATOM,KS(11,LMNKAT),IA(LMNATM),JRCH
      SAVE /BON/,/ATT/,ICH,ICHN
      DIMENSION KP(2,500),U(500),IN(3,500),B(500)
      INTEGER ICH(11,2),ICHN(11,2)
      DATA ICH/1,8,5,9,4,8,10,0,0,0,0,
     &                    1,4,10,1,9,11,4,10,12,6,8/
      DATA ICHN/0,1,4,6,9,10,13,0,0,0,0,
     &                     0,1,2,4,5,7,9,10,12,14,15/
      NDIM=0
      NT=0
      NS=0
      DO 1 I=1,NKATOM
      IF(JRCH.EQ.1) JN=7
      IF(JRCH.EQ.2) JN=11
      DO 2 J=1,JN
      IF(KS(J,I).EQ.0) GO TO 2
      JR=ICH(J,JRCH)
      JJ=ICHN(J,JRCH)
      II=I
      CALL TSLCLS(IR,NA,JR,II,KP,U,IN,B,ND1,ND2,NND)
      IF(ND1.EQ.0) GO TO 2
      NDIM=NDIM+ND1-ND2
      NSS=NS
      DO 3 K=1,ND1
      NS=NS+1
      A(NS)=B(K)
      INS(1,NS)=NT+1
      INS(2,NS)=NT+IN(2,K)-IN(1,K)+1
      INS(3,NS)=0
      IF(IN(3,K).GT.0) INS(3,NS)=NSS+IN(3,K)
      IF(IN(3,K).LT.0) INS(3,NS)=-(NSS-IN(3,K))
      DO 4 I1=IN(1,K),IN(2,K)
      NT=NT+1
      IF(NT.LE.1500) GO TO 91
      WRITE(6,901)
  901 FORMAT(' STOP AT 91 IN TSBASE')
      STOP
   91 CONTINUE
      V(NT)=U(I1)
      KA(1,NT)=KP(1,I1)
    4 KA(2,NT)=KP(2,I1)+JJ
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
      RETURN
      ENTRY TSBSDS
      IF(NT.EQ.0) GO TO 5
      WRITE(6,663) IR,NA,NDIM,NS,NS-NDIM,NT
  663 FORMAT(//6I5)
      IS=0
      DO 6 I=1,NT
      IF(I.NE.INS(1,IS+1).OR.IS.GE.NS) GO TO 7
      IS=IS+1
      WRITE(6,661) IS,(INS(K,IS),K=1,3),A(IS)
  661 FORMAT(I4,3I5,F10.5)
      IC=0
    7 IC=IC+1
      WRITE(6,660) I,IC,(KA(K,I),K=1,2),V(I)
  660 FORMAT(4X,2I5,2H (,2I5,4H)  (,2F10.5,1H))
    6 CONTINUE
    5 CONTINUE
      RETURN
      END
C SUBROUTINE TSBOGN ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSBOGN(JBO,BO,NB,JFBOND)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 WWST,V
      COMMON/BON/KBO(2,30),BOL(30),NBO,KBOND(4,200),NBOND,WWST(4000)
     &      ,IVEC(10000),VST(3,4000),NVEC,NDIM
     &      ,KA(2,1500),V(1500),INS(3,500),A(500),NT,NS
      COMMON/ATT   /ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &   ,NKATOM,NATOM,KS(11,LMNKAT),IA(LMNATM),JRCH
      SAVE /BON/,/ATT/,NSCH
      DIMENSION JBO(2,10),BO(10)
      INTEGER NSCH(11,2)
      DATA NSCH/1,2,3,3,4,4,4,0,0,0,0,
     &                     1,2,2,3,3,3,4,4,4,4,4/
      WRITE(6,602)
  602 FORMAT(1H ,' TABLE OF GENERATED BONDS')
      NBO=NB
      DO 1 I=1,NBO
      BOL(I)=BO(I)
      IF(JBO(2,I).LT.JBO(1,I)) GO TO 2
      KBO(1,I)=JBO(1,I)
      KBO(2,I)=JBO(2,I)
      GO TO 1
    2 KBO(1,I)=JBO(2,I)
      KBO(2,I)=JBO(1,I)
    1 CONTINUE
      NBOND=0
      DO 3 I=1,NBO
      N1=7
      IF(JRCH.EQ.2) N1=11
      DO 4 I1=1,N1
      IF(KS(I1,KBO(1,I)).EQ.0) GO TO 4
      N2=1
      IF(KBO(1,I).EQ.KBO(2,I)) N2=I1
      DO 5 I2=N2,N1
      IF(KS(I2,KBO(2,I)).EQ.0) GO TO 5
      N3=MIN0(NSCH(I1,JRCH),NSCH(I2,JRCH))
      IF(JFBOND.EQ.0.AND.N3.EQ.4) GO TO 5
      IF(JFBOND.EQ.9.AND.N3.EQ.4
     &   .AND.I.EQ.1) GO TO 5
      DO 6 J=1,N3
      NBOND=NBOND+1
      IF(NBOND.LE.200) GO TO 91
      WRITE(6,901)
  901 FORMAT(' STOP AT 91 IN TSBOGN')
      STOP
   91 CONTINUE
      KBOND(1,NBOND)=I
      KBOND(2,NBOND)=I1
      KBOND(3,NBOND)=I2
      KBOND(4,NBOND)=J-1
      WRITE(6,600) NBOND,(KBO(K,I),K=1,2),(KBOND(K,NBOND),K=1,4)
  600 FORMAT(7I5)
    6 CONTINUE
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
C     NNAA=NKATOM*7
C     IF(JRCH.EQ.2) NNAA=NKATOM*11
C     WRITE(6,601) NNAA,NBOND
C     WRITE(1) NNAA,NBOND
C 601 FORMAT(2I5)
      RETURN
      END
C SUBROUTINE TSBOCF ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSBOCF(KB,ICB,NKB,IIR,JF,KF)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 WWST,V,WA,WC
      COMMON/BON/KBO(2,30),BOL(30),NBO,KBOND(4,200),NBOND,WWST(4000)
     &       ,IVEC(10000),VST(3,4000),NVEC,NDIM
     &       ,KA(2,1500),V(1500),INS(3,500),A(500),NT,NS
      COMMON/ATT   /ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &     ,NKATOM,NATOM,KS(11,LMNKAT),IA(LMNATM),JRCH
      SAVE /BON/,/ATT/,NNKB,ISM
      INTEGER ISM(16,2)
      DIMENSION VD(3),KB(3),VDC(3)
      DIMENSION NSNS(40),JAM(40),JSM(40)
      DIMENSION NDEG(12),JTR(12),IPART(12)
      DATA ISM/1,2,2,2,3,3,4,4,4,5,6,6,6,7,7,7,
     &         1,2,3,3,4,5,5,6,6,7,8,8,9,9,10,11/
      DATA NNKB/0/
      NAB=NATOM*NATOM
      IF(NNKB.EQ.NKB) GO TO 3
      CALL DGTRST(JDUB,NNR,MMG,NSTR,NDEG,JTR,IPART)
      IF(KF.NE.0) WRITE(KF,611) KB,ICB,MMG,NNR,NSTR
  611 FORMAT(/' TSBOCF FOR K=',3I3,'/',I3
     &   /' ORDER OF K-POINT GROUP=',I3
     &   /' NUMBER OF REPRESENTATION=',I2
     &   /' NUMBER OF STAR=',I2)
      IF(KF.NE.0) WRITE(KF,612) (I,I=1,NNR)
  612 FORMAT(' NO  ',12I5)
      IF(KF.NE.0) WRITE(KF,613) (NDEG(I),I=1,NNR)
  613 FORMAT(' NDEG',12I5)
      IF(KF.NE.0) WRITE(KF,614) (JTR(I),I=1,NNR)
  614 FORMAT(' JTR ',12I5)
      IF(KF.NE.0) WRITE(KF,615) (IPART(I),I=1,NNR)
  615 FORMAT(' PRTN',12I5)
      DO 1 I=1,NVEC
      WW=0.0D0
      DO 2 K=1,3
    2 WW=WW+(DBLE(KB(K))/DBLE(ICB))*VST(K,I)
      WW=2.0D0*3.14159265358979312D0*WW
    1 WWST(I)=DCMPLX(COS(WW),SIN(WW))
      NNKB=NKB
    3 CONTINUE
      IF(NS.EQ.0) RETURN
      N711=7
      IF(JRCH.EQ.2) N711=11
      J=0
      DO 21 I=1,NS
      IF(INS(3,I).LT.0) GO TO 21
      J=J+1
      JA=IA(KA(1,INS(1,I)))
      JS=ISM(KA(2,INS(1,I)),JRCH)
      NSNS(J)=N711*(JA-1)+JS
      JAM(J)=JA
      JSM(J)=JS
   21 CONTINUE
      IF(KF.NE.0) WRITE(KF,603) NKB,IIR,J
  603 FORMAT(' NKB=',I3,' REPRESENTAION NO=',I2
     &     ,' DIMENSION OF MATRIX=',I3)
      IF(KF.NE.0) WRITE(KF,602) (' (',JAM(I),JSM(I),')',I=1,J)
  602 FORMAT(8(A2,2I3,A1))
      WRITE(JF) NKB,IIR,J,MMG,NSTR,NDEG(IIR),JTR(IIR),IPART(IIR)
      WRITE(JF) (NSNS(I),I=1,J)
      IJ=JRCH
      IT=0
      IZ=0
      DO 4 I=1,NS
      IF(INS(3,I).LT.0) GO TO 4
      IZ=IZ+1
      JZ=IZ-1
      DO 5 J=I,NS
      IF(INS(3,J).LT.0) GO TO 5
      JZ=JZ+1
      IX=I
      JX=J
      IF(INS(3,IX).EQ.0.OR.INS(3,JX).GT.0) GO TO 13
      IX=J
      JX=I
      IW=IZ
      IZ=JZ
      JZ=IW
   13 J1=INS(1,JX)
      J2=INS(2,JX)
      I1=INS(1,IX)
      IF(INS(3,IX).EQ.0) GO TO 6
      J1=INS(1,INS(3,JX))
      J2=INS(2,INS(3,JX))
    6 JA=IA(KA(1,INS(1,IX)))
      JB=IA(KA(1,INS(1,JX)))
      IS=ISM(KA(2,INS(1,IX)),IJ)
      JS=ISM(KA(2,INS(1,JX)),IJ)
      IF(JA.LE.JB) GO TO 11
      IW=JB
      JB=JA
      JA=IW
      IW=JS
      JS=IS
      IS=IW
   11 IF(JA.NE.JB.OR.IS.LE.JS) GO TO 12
       IW=IS
       JS=IS
      IS=IW
   12 DO 7 K=1,NBOND
      IF(KBO(1,KBOND(1,K)).NE.JA) GO TO 7
      IF(KBO(2,KBOND(1,K)).NE.JB) GO TO 7
      IF(KBOND(2,K).NE.IS) GO TO 7
      IF(KBOND(3,K).NE.JS) GO TO 7
      WA=0.0D0
      DO 8 JJ=J1,J2
      IN=(KBOND(1,K)-1)*NAB+(KA(1,I1)-1)*NATOM+KA(1,JJ)
      IN1=MOD(IVEC(IN),10000)
      IN2=IVEC(IN)/10000
      WC=0.0D0
      IF(IN1.EQ.0) GO TO 8
      DO 9 IIN=IN1,IN2
      DO 10 KK=1,3
   10 VD(KK)=VST(KK,IIN)
      CALL ZZZY43(VD,VDC,SSS)
      WC=WC+WWST(IIN)*TSKCOF(KBOND(4,K),KA(2,I1),KA(2,JJ),VDC,IJ,0)
    9 CONTINUE
      WA=WA+WC*V(JJ)
    8 CONTINUE
      WA=WA/(A(IX)*A(JX))
      IF(CDABS(WA).LT.1.0D-4) GO TO 7
      IT=IT+1
      IF(KF.NE.0) WRITE(KF,601) IT,IZ,JZ
     &     ,K,(KBO(II,KBOND(1,K)),II=1,2),(KBOND(II,K),II=1,4),WA
  601 FORMAT(I5,3I3,2H (,6I3,4H)  (,2F10.5,1H))
      WRITE(JF) IZ,JZ,K,WA
    7 CONTINUE
      IF(IZ.LE.JZ) GO TO 5
      IW=IZ
      IZ=JZ
      JZ=IW
    5 CONTINUE
    4 CONTINUE
      WRITE(JF) 0,0,0,(0.0D0,0.0D0)
      RETURN
      END
C FUNCTION TSKCOF =2====*====3====*====4====*====5====*====6====*====7
C
C   SLATTER AND KOSTER COEFICENT
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      REAL*8 FUNCTION TSKCOF(IS,IT,JT,V,IJ,IPR)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 TT,W,WC,WB,WD
      COMMON/SCC/TT(16,7)
      SAVE /SCC/,ISTA,IJST,LM
      DIMENSION V(3),IC(16),ICM(4,7,4,3),AK(7,4,3),ISM(7,4,3)
      INTEGER LM(16)
      DATA LM/0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3/
      DATA ISTA,IJST/0,0/
      IF(IJST.EQ.IJ) GO TO 10
      IJST=IJ
      CALL ZZZY46(IJST)
   10 CONTINUE
      IF(ISTA.EQ.1) GO TO 1
      ISTA=1
      DO 2 I=1,3
      L=I
      DO 2 J=1,2*L+1
      MA=I-J+1
      DO 2 K=1,L+1
      MB=I-K+1
      CALL TSRMI(L,MA,MB,K1,K2,K3,ISA,IC,0)
      ISM(J,K,I)=ISA
      AK(J,K,I)=(DBLE(K1)/DBLE(K2))*SQRT(DBLE(K3))
      DO 3 M=1,4
    3 ICM(M,J,K,I)=IC(M)
    2 CONTINUE
    1 CONTINUE
      SB=SQRT(1.0D0-V(3)*V(3))
      SX=1.0D0
      SY=0.0D0
      IF(SB.GT.1.0D-4) SX=V(1)/SB
      IF(SB.GT.1.0D-4) SY=V(2)/SB
      WD=DCMPLX(SX,SY)
      II=IT
      IGO=1
      L1=LM(IT)
    9 K1=L1-IS+1
      W=1.0D0
      IF(L1.EQ.0) GO TO 4
      W=0.0D0
      DO 5 K=1,2*L1+1
      IF(CDABS(TT(II,K)).LT.1.0D-4) GO TO 5
      KK=L1-K+1
      WA=ICM(1,K,K1,L1)
      DO 6 I=2,L1+1
    6 WA=WA+DBLE(ICM(I,K,K1,L1))*(V(3)**(I-1))
      WA=WA*AK(K,K1,L1)
      IF(ISM(K,K1,L1).EQ.2) WA=WA*SB
      IF(KK.EQ.0) WB=WA
      IF(KK.GT.0) WB=WA*(CONJG(WD)**KK)
      IF(KK.LT.0) WB=WA*(WD**(-KK))
      W=W+WB*CONJG(TT(II,K))
    5 CONTINUE
    4 GO TO (7,8), IGO
    7 WC=W
      IF(IT.EQ.JT) GO TO 8
      L1=LM(JT)
      II=JT
      IGO=2
      GO TO 9
    8 WA=WC*CONJG(W)
      IF(IS.NE.0) WA=WA*2.0D0
      IF(LM(IT).GT.LM(JT).AND.
     &   MOD(LM(IT)+LM(JT),2).EQ.1) WA=-WA
      TSKCOF=WA
      IF(IPR.EQ.1) WRITE(6,600) IS,IT,JT,V,WA
  600 FORMAT(' TSKCOF',3I5,3F9.5,F9.5)
      RETURN
      END
C SUBROUTINE LATT01 ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE LATT01(JF)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/LAT   / A,B,C,CA,CB,CC,A1,B1,C1
      COMMON/RMTRC0/ RM(3,3)
      COMMON/GLATT0/ AK,BK,CK,BCK,CAK,ABK
      COMMON/GMTRC0/ GM(3,3)
      PAI2=8.D0*DATAN(1.D0)
      CAA=CA**2
      CBB=CB**2
      CCC=CC**2
      CW=1.0-CAA-CBB-CCC+2.0*CA*CB*CC
      SA=DSQRT(1.D0-CAA)
      SB=DSQRT(1.D0-CBB)
      SC=DSQRT(1.D0-CCC)
      AK=PAI2*A1/A
      BK=PAI2*B1/B
      CK=PAI2*C1/C
      BCK=(CB*CC-CA)/(SB*SC)
      CAK=(CC*CA-CB)/(SC*SA)
      ABK=(CA*CB-CC)/(SA*SB)
      IF(JF.GE.1 .OR. JF.LE.99) THEN
        WRITE(JF,100) A,B,C,CA,CB,CC
        WRITE(JF,200) AK,BK,CK,BCK,CAK,ABK
        WRITE(JF,300)
      END IF
  100 FORMAT(1H ,'--- LATTICE PARAMETERS IN REAL SPACE ---------'/
     &       1H ,' A, B, C=',3F16.10/
     &       1H ,'CA,CB,CC=',3F17.10                         )
  200 FORMAT(1H ,'--- LATTICE PARAMETERS IN RECIPROCAL SPACE ---'/
     &       1H ,' A, B, C=',3F16.10/
     &       1H ,'CA,CB,CC=',3F17.10                         )
  300 FORMAT(1H ,'----------------------------------------------')
C
      RM(1,1)=A*A
      RM(2,2)=B*B
      RM(3,3)=C*C
      RM(1,2)=A*B*CC
      RM(2,3)=B*C*CA
      RM(3,1)=C*A*CB
      RM(2,1)=RM(1,2)
      RM(3,2)=RM(2,3)
      RM(1,3)=RM(3,1)
      GM(1,1)=AK*AK
      GM(2,2)=BK*BK
      GM(3,3)=CK*CK
      GM(1,2)=AK*BK*ABK
      GM(2,3)=BK*CK*BCK
      GM(3,1)=CK*AK*CAK
      GM(2,1)=GM(1,2)
      GM(3,2)=GM(2,3)
      GM(1,3)=GM(3,1)
      RETURN
      END
C SUBROUTINE LATT01 ====*====3====*====4====*====5====*====6====*====7
C
C
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE GLATTC(NATM)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 UU
      INCLUDE 'prmtsp.f'
      SAVE KP1
C     PARAMETER(MAXNPW=4854)
      COMMON/SPW   / KM(4,MAXNPW),AA(MAXNPW),UU(MAXNPW),KT(MAXNPW),IE,IK
      INTEGER KP1(4) 
      DATA KP1/0, 0, 0, 1/
C
      IKMIN=NATM*150
      IF(IKMIN.GT.MAXNPW) IKMIN=MAXNPW
      AMX=1.D0
    1 AMX=AMX+0.5D0
      CALL ZZZY41(KP1(1),KP1(4),AMX)
      WRITE(6,*) 'IK=',IK,'   IE=',IE
      IF(IK.LT.IKMIN) GO TO 1
      DO 10 I=1,IE
   10 WRITE(6,100) (KM(K,KT(I)),K=1,3),KT(I),AA(KT(I))
  100 FORMAT(1H ,'(',3I5,' )','  (',I5,' )   ',D14.5  )
      RETURN
      END
C SUBROUTINE TSBZEG ====*====3====*====4====*====5====*====6====*====7
C
C   EDGES OF THE  FIRST B.Z.
C
C                 1988.10.18 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSBZEG(NRECPT,RECTAX,NRP,CO,NLIN)
        IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NRECPT(3,30),RECTAX(4,30)
      DIMENSION CO(3,2,100)
      CALL FIBZPL(NRECPT,NRP)
      CALL RECTAG(NRECPT,RECTAX,NRP)
      NLIN=0
      DO 31 N1=1,NRP-1
      NN1=N1
      DO 32 N2=N1+1,NRP
      NN2=N2
      NNCO=0
      DO 33 N3=1,NRP
      NN3=N3
c      WRITE(6,*) N1,N2,N3
      IF(N3.EQ.N1.OR.N3.EQ.N2) GO TO 33
         CALL FICORN(NN1,NN2,NN3,RECTAX,NRP,X,Y,Z,IND)
         IF(IND.EQ.0) THEN
c         write(6,*) x,y,z
            IF(NNCO.EQ.0) THEN
               NNCO=1
               NLIN=NLIN+1
               CO(1,1,NLIN)=X
               CO(2,1,NLIN)=Y
               CO(3,1,NLIN)=Z
            ELSE
               WA=(CO(1,1,NLIN)-X)**2
     &           +(CO(2,1,NLIN)-Y)**2
     &           +(CO(3,1,NLIN)-Z)**2
               IF(WA.GT.0.000001D0) THEN
                   CO(1,2,NLIN)=X
                   CO(2,2,NLIN)=Y
                   CO(3,2,NLIN)=Z
                   GO TO 32
               END IF
            END IF
         END IF
   33 CONTINUE
      IF(NNCO.EQ.1) THEN
          NLIN=NLIN-1 
      END IF
   32 CONTINUE
   31 CONTINUE
      RETURN
      END
C SUBROUTINE TSCSDT ====*====3====*====4====*====5====*====6====*====7
C
C          DATA FOR CRYSTAL STRUCTURE PLOT
C IN-PUT
C   ATOMIC POSITIONS FOR NA*NB*NC UNIT CELLS
C   IF(NA.EQ.0) THEN NC HEXAGONAL UNIT CELLS
C   KABC=1,2,3,4,5,6 DETERMINE THE DIRECTION. SEE DATA OF I123
C   XOUT(3,30)    POSITION OF EXTRA ATOMS,
C   NKO(30)       NUMBERS OF EXTRA ATOMS IN TSCRST
C   NOUT          NUMBER OF EXTRA ATOMS
C OUT-PUT
C   XFU(3,2)     FIGURE SIZE IN UNIT OF LATTICE CONSTANT I123(1,KABC)
C                IN RECTANGULAR COORDINATE.
C   XC(3,300)    CENTER OF SPHERES
C   IK(300)      NUMBERS OF ATOMS
C   NNS          NUMBER OF SPHERES
C   JB(2,300)    NUMBERS OF END SPHERES FOR BARS
C   IBB(300)     KINDS OF BONDS
C   NBB          NUMBER OF BARS
C   X1(3,300),X2(3,300)  END POINTS OF LINES
C   NLINE        NUMBER OF LINES
C
C               1989/12/19      AKIRA YANASE
C                  ORIGINAL PROGRAM IS TSLSPL
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSCSDT(NA,NB,NC,KABC,XOUT,NKO,NOUT
     &    ,XFU,XC,IK,NNS,JB,IBB,NBB,X1,X2,NLINE)
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 WWST,V
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/LAT/A,B,C,CA,CB,CC,A1,B1,C1
      COMMON/ATT/ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &      ,NKATOM,NATOM,KS(11,LMNKAT),IA(LMNATM),JRCH
      COMMON/BON/KBO(2,30),BOL(30),NBO,KBOND(4,200),NBOND,WWST(4000)
     &       ,IVEC(10000),VST(3,4000),NVEC,NDIM
     &       ,KA(2,1500),V(1500),INS(3,500),AA(500),NT,NS
      SAVE /SPG2/,/LAT/,/BON/,/ATT/,I123
      INTEGER I123(3,6)
      DIMENSION XOUT(3,30),NKO(30)
      REAL*8 XFU(3,2),XC(3,500),X1(3,500),X2(3,500)
      INTEGER IK(500),JB(2,500),IBB(500)
      DIMENSION N(3),II(3),X(3),XS(3,500),INNS(3,50)
      DIMENSION WCOS(3),ABC(3),XA(3),XX(2,6)
      DATA I123/
     &  1,2,3, 1,3,2, 2,3,1, 2,1,3, 3,1,2, 3,2,1/
C     ATOM POSITIONS
C
      WRITE(6,*) ' POSITION OF SPHERES'
      WRITE(6,*) ' NO  COORDINATE IN A,B,C UNIT      KIND NATOM'
      NNS=0
      N(1)=NA
      N(2)=NB
      N(3)=NC
      IHEX=0
      IF(NA.GT.0) GO TO 1001
      IF(IL.GT.0) RETURN
      IHEX=1
      N(1)=2
      N(2)=2
 1001 CONTINUE
      WRITE(6,6006) IL,IHEX,N
 6006 FORMAT(10I5)
      DO 100 IS=1,NATOM
      INNS(1,IS)=NNS+1
      DO 101 IAA=1,N(1)+2
      II(1)=IAA-2
      DO 102 IB=1,N(2)+2
      II(2)=IB-2
      DO 103 IC=1,N(3)+2
      II(3)=IC-2
      K=1
      IF(IL.EQ.-1) K=3
      IF(IL.EQ.2) K=4
      IF(IL.EQ.3) K=2
      IF(IL.EQ.4) K=2
      DO 110 J=1, K
      DO 111 I=1,3
      W=VATOM(I,IS)
      IF(J.EQ.1) GO TO 112
      IF(IL.NE.-1) GO TO 113
      W=W+1.0/3.0
      IF(J.EQ.2.AND.I.EQ.1) W=W+1.0/3.0
      IF(J.EQ.3.AND.I.NE.1) W=W+1.0/3.0
      GO TO 112
  113 IF(IL.EQ.3) W=W+0.5
      IF(IL.EQ.2.AND.I.NE.J-1) W=W+0.5
      IF(IL.EQ.4.AND.I.NE.3) W=W+0.5
  112 W=W+II(I)
      IF(W-N(I).GT.1.0E-4) GO TO 110
      IF(W.LT.-1.0E-4) GO TO 110
      X(I)=W
  111 CONTINUE
      IF(IHEX.EQ.0) GO TO 1002
      IF(X(1)-X(2)-1.0.GT.1.0E-4) GO TO 110
      IF(X(1)-X(2)+1.0.LT.-1.0E-4) GO TO 110
 1002 CONTINUE
      IF(NNS.GE.300) STOP
      NNS=NNS+1
      DO 114 I=1,3
  114 XS(I,NNS)=X(I)
      IK(NNS)=IS
      WRITE(6,6001) NNS,X,IA(IS),IS
 6001 FORMAT(I5,3F10.5,2I5)
  110 CONTINUE
  103 CONTINUE
  102 CONTINUE
  101 CONTINUE
      INNS(3,IS)=NNS
      IF(NOUT.EQ.0) GO TO 116
      DO 115 I=1,NOUT
      IF(NKO(I).NE.IS) GO TO 115
      IF(NNS.GE.300) STOP
      NNS=NNS+1
      DO 117 K=1,3
  117 XS(K,NNS)=XOUT(K,I)
      IK(NNS)=IS
      WRITE(6,6001) NNS,(XS(K,NNS),K=1,3),IA(IS),IS
  115 CONTINUE
  116 INNS(2,IS)=NNS
  100 CONTINUE
      NSS=NNS
C
C     BARS BETWEEN ATOMS
C
      NAB=NATOM*NATOM
      NBB=0
      WRITE(6,*) ' BARS BETWEEN SPHERE'
      WRITE(6,*) ' NO   END SPHERE  BOND  END POINT WITHOUT SPHERE'
      DO 202 IS=1,NATOM
      INS1=INNS(1,IS)
      INS2=INNS(2,IS)
      DO 203 IT=1,NATOM
      ISIG=1
      ISA=IS
      ITA=IT
      IF(IS.LE.IT) GO TO 214
      ISIG=-1
      ISA=IT
      ITA=IS
  214 CONTINUE
      DO 201 IB=1,NBO
      IF(KBO(1,IB).NE.IA(ISA)) GO TO 201
      IF(KBO(2,IB).NE.IA(ITA)) GO TO 201
      JNS1=INNS(1,IT)
      JNS2=INNS(2,IT)
      IN=(IB-1)*NAB+(ISA-1)*NATOM+ITA
      IN1=MOD(IVEC(IN),10000)
      IN2=IVEC(IN)/10000
      IF(IN1.EQ.0) GO TO 201
      DO 204 INSS=INS1,INS2
      DO 205 I1N=IN1,IN2
      DO 206 KK=1,3
      XA(KK)=XS(KK,INSS)+VST(KK,I1N)*ISIG
  206 CONTINUE
C     WRITE(6,6002) I1N,INSS,IT,IS,IB,XA
C6002 FORMAT(5I5,3F10.5)
      DO 207 JNS=JNS1,JNS2
      DO 208 KK=1,3
      IF(ABS(XA(KK)-XS(KK,JNS)).GT.1.0E-5) GO TO 207
  208 CONTINUE
      IF(IT.LT.IS) GO TO 205
      IF(NBB.GE.300) STOP
      NBB=NBB+1
      JB(1,NBB)=INSS
      JB(2,NBB)=JNS
      IBB(NBB)=IB
      WRITE(6,6003) NBB,(JB(I,NBB),I=1,2),IBB(NBB)
 6003 FORMAT(4I5)
      GO TO 205
  207 CONTINUE
      DO 209 KK=1,3
      IF(XA(KK).GT.N(KK)+1.0E-4) GO TO 210
      IF(XA(KK).GT.-1.0E-4) GO TO 209
      IF(XS(KK,INSS).LT.1.0E-4) GO TO 205
      WA=XS(KK,INSS)/(XS(KK,INSS)-XA(KK))
      GO TO 211
  210 IF(XS(KK,INSS).GT.N(KK)-1.0E-4) GO TO 205
      WA=(N(KK)-XS(KK,INSS))/(XA(KK)-XS(KK,INSS))
  211 DO 212 K=1,3
  212 XA(K)=XS(K,INSS)+(XA(K)-XS(K,INSS))*WA
  209 CONTINUE
      IF(IHEX.EQ.0) GO TO 1003
      IF(XA(1)-XA(2)-1.0.GT.1.0E-4) GO TO 2004
      IF(XA(1)-XA(2)+1.0.GE.-1.0E-4) GO TO 1003
      S=1.0
      GO TO 2005
 2004 S=-1.0
 2005 WA=XA(2)-XS(2,INSS)-(XA(1)-XS(1,INSS))
      W1=(XS(1,INSS)*XA(2)-XS(2,INSS)*XA(1)
     &  +S*(XA(1)-XS(1,INSS)))/WA
      W2=W1+S
      XA(3)=((XS(1,INSS)-XS(2,INSS))*XA(3)
     &     -(XA(1)-XA(2))*XS(3,INSS)+S*(XA(3)-XS(3,INSS)))/WA
      XA(1)=W1
      XA(2)=W2
 1003 CONTINUE
      IF(NNS.GE.300) STOP
      NNS=NNS+1
      DO 213 K=1,3
  213 XS(K,NNS)=XA(K)
      IK(NNS)=0
      IF(NBB.GE.300) STOP
      NBB=NBB+1
      JB(1,NBB)=INSS
      JB(2,NBB)=NNS
      IBB(NBB)=IB
      WRITE(6,6004) NBB,(JB(I,NBB),I=1,2),IBB(NBB),XA
 6004 FORMAT(4I5,3F10.5)
  205 CONTINUE
  204 CONTINUE
  201 CONTINUE
  203 CONTINUE
  202 CONTINUE
C
C     COEFFICIENTS OF TRANSFORMATION
C
      WCOS(1)=CA
      WCOS(2)=CB
      WCOS(3)=CC
      ABC(1)=A
      ABC(2)=B
      ABC(3)=C
      I1=I123(1,KABC)
      I2=I123(2,KABC)
      I3=I123(3,KABC)
      WA=ABC(I1)
      ABC(1)=ABC(1)/WA
      ABC(2)=ABC(2)/WA
      ABC(3)=ABC(3)/WA
      WSIN=SQRT(1.0-WCOS(I3)**2)
      F1=ABC(I1)
      F2=ABC(I2)*WCOS(I3)
      F3=ABC(I3)*WCOS(I2)
      F4=ABC(I2)*WSIN
      F5=ABC(I3)*(WCOS(I1)-WCOS(I2)*WCOS(I3))/WSIN
      F6=ABC(I3)*SQRT(1.0-CA*CA-CB*CB-CC*CC
     &                  +2.0*CA*CB*CC)/WSIN
C
C     FIGURE SIZE
C
      XFU(1,2)=DMAX1(F1*N(I1),F1*N(I1)+F2*N(I2),
     &  F1*N(I1)+F3*N(I3),F1*N(I1)+F2*N(I2)+F3*N(I3))
      XFU(1,1)=DMIN1(0.0D0,F2*N(I2),F3*N(I3),F2*N(I2)+F3*N(I3))
      XFU(2,2)=DMAX1(F4*N(I2),F4*N(I2)+F5*N(I3))
      XFU(2,1)=DMIN1(0.0D0,F5*N(I3))
      XFU(3,1)=0.0
      XFU(3,2)=F6*N(I3)
      IF(NOUT.EQ.0) GO TO 12
      DO 11 I=1,NOUT
      XO=F1*XOUT(I1,I)+F2*XOUT(I2,I)+F3*XOUT(I3,I)
      IF(XO.GT.XFU(1,2)) XFU(1,2)=XO
      IF(XO.LT.XFU(1,1)) XFU(1,1)=XO
      YO=F4*XOUT(I2,I)+F5*XOUT(I3,I)
      IF(YO.GT.XFU(2,2)) XFU(2,2)=YO
      IF(YO.LT.XFU(2,1)) XFU(2,1)=YO
      ZO=F6*XOUT(I3,I)
      IF(ZO.GT.XFU(3,2)) XFU(3,2)=ZO
      IF(ZO.LT.XFU(3,1)) XFU(3,1)=ZO
   11 CONTINUE
   12 CONTINUE
C
C     TPERSP
C
      DO 51 I=1,NNS
      XC(1,I)=F1*XS(I1,I)+F2*XS(I2,I)+F3*XS(I3,I)
      XC(2,I)=F4*XS(I2,I)+F5*XS(I3,I)
      XC(3,I)=F6*XS(I3,I)
   51 CONTINUE
      NLINE=0
      IF(IHEX.EQ.1) GO TO 1014
      DO 1004 J1=1,N(I1)+1
      DO 1005 J2=1,N(I2)+1
      NLINE=NLINE+1
      X1(1,NLINE)=F1*(J1-1)+F2*(J2-1)
      X1(2,NLINE)=F4*(J2-1)
      X1(3,NLINE)=0.0
      X2(1,NLINE)=X1(1,NLINE)+F3*N(I3)
      X2(2,NLINE)=X1(2,NLINE)+F5*N(I3)
      X2(3,NLINE)=F6*N(I3)
 1005 CONTINUE
 1004 CONTINUE
      DO 1006 J2=1,N(I2)+1
      DO 1007 J3=1,N(I3)+1
      NLINE=NLINE+1
      X1(1,NLINE)=F2*(J2-1)+F3*(J3-1)
      X1(2,NLINE)=F4*(J2-1)+F5*(J3-1)
      X1(3,NLINE)=F6*(J3-1)
      X2(1,NLINE)=X1(1,NLINE)+F1*N(I1)
      X2(2,NLINE)=X1(2,NLINE)
      X2(3,NLINE)=X1(3,NLINE)
 1007 CONTINUE
 1006 CONTINUE
      DO 1008 J3=1,N(I3)+1
      DO 1009 J1=1,N(I1)+1
      NLINE=NLINE+1
      X1(1,NLINE)=F1*(J1-1)+F3*(J3-1)
      X1(2,NLINE)=F5*(J3-1)
      X1(3,NLINE)=F6*(J3-1)
      X2(1,NLINE)=X1(1,NLINE)+F2*N(I2)
      X2(2,NLINE)=X1(2,NLINE)+F4*N(I2)
      X2(3,NLINE)=X1(3,NLINE)
 1009 CONTINUE
 1008 CONTINUE
      IF(IHEX.EQ.0) GO TO 1010
 1014 CONTINUE
      XX(1,1)=0.0
      XX(2,1)=0.0
      XX(1,2)=F1
      XX(2,2)=0.0
      XX(1,3)=F1*2.0+F2
      XX(2,3)=F4
      XX(1,4)=F1*2.0+F2*2.0
      XX(2,4)=F4*2.0
      XX(1,5)=F1+F2*2.0
      XX(2,5)=F4*2.0
      XX(1,6)=F2
      XX(2,6)=F4
      DO 1011 J1=1,N(3)+1
      DO 1012 J2=1,6
      JW=J2+1
      IF(J2.EQ.6) JW=1
      ZZ=(J1-1)*F6
      NLINE=NLINE+1
      X1(1,NLINE)=XX(1,J2)
      X1(2,NLINE)=XX(2,J2)
      X1(3,NLINE)=ZZ
      X2(1,NLINE)=XX(1,JW)
      X2(2,NLINE)=XX(2,JW)
      X2(3,NLINE)=ZZ
 1012 CONTINUE
 1011 CONTINUE
      DO 1013 J2=1,6
      NLINE=NLINE+1
      X1(1,NLINE)=XX(1,J2)
      X1(2,NLINE)=XX(2,J2)
      X1(3,NLINE)=0.0
      X2(1,NLINE)=XX(1,J2)
      X2(2,NLINE)=XX(2,J2)
      X2(3,NLINE)=F6*N(I3)
 1013 CONTINUE
 1010 CONTINUE
      RETURN
      END
C SUBROUTINE CHKNIR ====*====3====*====4====*====5====*====6====*====7
C
C     ROUTINE FOR IR NUMBER CHECK
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CHKNIR(IR,NR)
      IF(IR.LT.1) THEN
         WRITE(6,601) IR
  601 FORMAT(' NUMBER IR OF I.R. SHOULD BE POSITVE'
     &      /' HOWEVER NOW IR=',I3,' THEN STOP HERE')
      END IF
      IF(IR.GT.NR) THEN
         WRITE(6,602) NR,IR
  602 FORMAT(' NUMBER IR OF I.R. SHOULD BE .LE.',I2
     &      /' HOWEVER NOW IR=',I3,' THEN STOP HERE')
      STOP
      END IF
      RETURN
      END
C SUBROUTINE CHKNST ====*====3====*====4====*====5====*====6====*====7
C
C     ROUTINE FOR STATE NUMBER IN I.R. CHECK
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CHKNST(NA,ND)
      IF(NA.LT.1) THEN
         WRITE(6,601) NA
  601 FORMAT(' STATE NUMBER OF NA IN I.R. SHOULD BE POSITVE'
     &      /' HOWEVER NOW NA=',I5,' THEN STOP HERE')
      STOP
      END IF
      IF(NA.GT.ND) THEN
         WRITE(6,602) ND,NA
  602 FORMAT(' STATE NUMBER OF NA IN I.R SHOULD BE .LE.',I2
     &      /' HOWEVER NOW NA=',I3,' THEN STOP HERE')
      STOP
      END IF
      RETURN
      END
C SUBROUTINE CHKNKA ====*====3====*====4====*====5====*====6====*====7
C
C     ROUTINE FOR ATOM KIND NUMBER CHECK
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CHKNKA(IA,NK)
      IF(IA.LT.1) THEN
         WRITE(6,601) IA
  601 FORMAT(' KIND NUMBER OF IA SHOULD BE POSITVE'
     &      /' HOWEVER NOW IA=',I3,' THEN STOP HERE')
      STOP
      END IF
      IF(IA.GT.NK) THEN
         WRITE(6,602) IA,NK
  602 FORMAT(' KIND NUMBER IA SHOULD BE .LE.',I2
     &      /' HOWEVER NOW NA=',I3,' THEN STOP HERE')
      STOP
      END IF
      RETURN
      END
C SUBROUTINE CHKDNM ====*====3====*====4====*====5====*====6====*====7
C
C     ROUTINE FOR ATOM KIND NUMBER CHECK
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CHKDNM(IC)
      IF(IC.LT.1) THEN
         WRITE(6,601) IC
  601 FORMAT(' DENOMINATOR IC SHOULD BE POSITVE'
     &      /' HOWEVER NOW IC=',I3,' THEN STOP HERE')
      STOP
      END IF
      RETURN
      END
C SUBROUTINE TSPGRP ====*====3====*====4====*====5====*====6====*====7
C
C    SPACE GOUP IS CONSTRUCTRED FROM GENERATOR GIVEN BY CALLING 
C    TSGENR. IF INV00.EQ.1 THEN ORIGIN OF THE COORDINATE IS MOVED
C    SO THAT THE TRANSLATION FOR THE INVERSION OPERATION BECOMES 0.
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSPGRP(INV00)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG1/,/SPG2/,IB
      INTEGER JA(2,3,48),JB(2,3),JC(2,3),JD(2,3) 
      INTEGER IB(2,3,6)
      DATA IB/
     &  1,2, 1,2, 1,2,  0,1, 1,2, 1,2,  1,2, 0,1, 1,2,
     &  1,2, 1,2, 0,1,  2,3, 1,3, 1,3,  1,3, 2,3, 2,3/
      N=0
    1 NN=N
      N=NG
      DO 2 I=1,N
      DO 3 J=1,N
      IW=IM(IG(I),IG(J))
      DO 4 K=1,NG
      IF(IW.EQ.IG(K)) THEN
         KSAME=K
         ISAME=1
         GO TO 42
      END IF
    4 CONTINUE
      ISAME=0 
      NG=NG+1
      IG(NG)=IW
   42 CONTINUE  
      DO 41 K=1,3
      JB(1,K)=JV(1,K,J)
      JB(2,K)=JV(2,K,J)
      JC(1,K)=JV(1,K,I)
   41 JC(2,K)=JV(2,K,I)
      CALL ZZZY22(IG(I),JB)
      CALL ZZZY23(JB,JC,1)
      IF(ISAME.EQ.1) THEN
         DO 43 K=1,3
C         IF(JB(1,K).NE.JV(1,K,KSAME).OR.
C     &      JB(2,K).NE.JV(2,K,KSAME)) GO TO 44
         IF(JB(1,K)*JV(2,K,KSAME).NE.JB(2,K)*JV(1,K,KSAME)) GO TO 44
   43    CONTINUE
         GO TO 3
   44    IF(IL.EQ.0.OR.IL.EQ.1) GO TO 48 
         IF(IL.EQ.-1) KS=5
         IF(IL.EQ.2) KS=2
         IF(IL.EQ.3) KS=1
         IF(IL.EQ.4) KS=4
         IF(IL.EQ.-1) KL=6
         IF(IL.EQ.2) KL=4
         IF(IL.EQ.3) KL=1
         IF(IL.EQ.4) KL=4
         DO 45 KK=KS,KL
         DO 46 L=1,3
            JC(1,L)=IB(1,L,KK)
            JC(2,L)=IB(2,L,KK)
   46    CONTINUE
         CALL ZZZY23(JC,JB,1)
         DO 47 K=1,3
C         IF(JC(1,K).NE.JV(1,K,KSAME).OR.
C     &      JC(2,K).NE.JV(2,K,KSAME)) GO TO 45
           IF(JC(1,K)*JV(2,k,KSAME).NE.
     &       JC(2,K)*JV(1,K,KSAME)) GO TO 45
   47    CONTINUE
         GO TO 3
   45    CONTINUE
   48    WRITE(6,*) 'YOUR GENERATOR(S) ARE NOT CORRECT'
         WRITE(6,*) 'STOP AT 45 IN TSPGRP'
         STOP          
      ELSE
         DO 5 L=1,3
         JV(1,L,NG)=JB(1,L)
    5    JV(2,L,NG)=JB(2,L)
      END IF
    3 CONTINUE
    2 CONTINUE
      IF(N.LT.NG) GO TO 1
      DO 6 I=2,NG
      DO 9 J=1,3
      JB(1,J)=JV(1,J,I)
    9 JB(2,J)=JV(2,J,I)
      IW=IG(I)
      DO 7 J=1,I-1
      JJ=I-J
      IF(IW.GT.IG(JJ)) GO TO 11
      IG(JJ+1)=IG(JJ)
      DO 8 K=1,3
      JV(1,K,JJ+1)=JV(1,K,JJ)
    8 JV(2,K,JJ+1)=JV(2,K,JJ)
    7 CONTINUE
      JJ=0
   11 IG(JJ+1)=IW
      DO 10 K=1,3
      JV(1,K,JJ+1)=JB(1,K)
   10 JV(2,K,JJ+1)=JB(2,K)
    6 CONTINUE
      IF(INV00.NE.1) GO TO 12
      DO 13 I=1,NG
      IF(IL.LE.0.AND.IG(I).EQ.13) GO TO 14
      IF(IL.GT.0.AND.IG(I).EQ.25) GO TO 14
   13 CONTINUE
      GO TO 12
   14 DO 15 J=1,3
      IW=JV(1,J,I)
      JW=JV(2,J,I)*2
      CALL ZZZY24(IW,JW)
      JB(1,J)=IW
      JB(2,J)=JW
   15 CONTINUE
      DO 16 I=1,NG
      DO 17 J=1,3
      JD(1,J)=JV(1,J,I)
      JD(2,J)=JV(2,J,I)
      JC(1,J)=JB(1,J)
   17 JC(2,J)=JB(2,J)
      CALL ZZZY22(IG(I),JC)
      CALL ZZZY23(JC,JB,-1)
      CALL ZZZY23(JD,JC,1)
      DO 18 J=1,3
      JV(1,J,I)=JD(1,J)
   18 JV(2,J,I)=JD(2,J)
   16 CONTINUE
   12 IF(IL.EQ.0.OR.IL.EQ.1) RETURN
      DO 19 I=1,3
      JA(1,I,1)=0
   19 JA(2,I,1)=1
      NV=1
      DO 20 I=1,NG
      DO 21 J=1,NV
      DO 23 K=1,3
      IF(JV(1,K,I).NE.JA(1,K,J).OR.
     &   JV(2,K,I).NE.JA(2,K,J)) GO TO 21
   23 CONTINUE
      GO TO 20
   21 CONTINUE
      IF(IL.EQ.-1) KS=5
      IF(IL.EQ.2) KS=2
      IF(IL.EQ.3) KS=1
      IF(IL.EQ.4) KS=4
      IF(IL.EQ.-1) KL=6
      IF(IL.EQ.2) KL=4
      IF(IL.EQ.3) KL=1
      IF(IL.EQ.4) KL=4
      DO 24 K=KS,KL
      DO 25 L=1,3
      JC(1,L)=IB(1,L,K)
      JC(2,L)=IB(2,L,K)
      JD(1,L)=JV(1,L,I)
   25 JD(2,L)=JV(2,L,I)
      CALL ZZZY23(JC,JD,1)
      DO 26 L=1,NV
      DO 27 N=1,3
      IF(JC(1,N).NE.JA(1,N,L).OR.
     &   JC(2,N).NE.JA(2,N,L)) GO TO 26
   27 CONTINUE
      DO 28 N=1,3
      JV(1,N,I)=JA(1,N,L)
   28 JV(2,N,I)=JA(2,N,L)
      GO TO 20
   26 CONTINUE
   24 CONTINUE
      NV=NV+1
      DO 29 K=1,3
      JB(1,K)=JD(1,K)
   29 JB(2,K)=JD(2,K)
      DO 30 K=KS,KL
      DO 31 L=1,3
      JC(1,L)=IB(1,L,K)
   31 JC(2,L)=IB(2,L,K)
      CALL ZZZY23(JC,JD,1)
      DO 32 L=1,3
      IF(JB(1,L).LT.JC(1,L)) GO TO 30
      IF(JB(1,L).GT.JC(1,L)) GO TO 33
   32 CONTINUE
   33 DO 34 L=1,3
      JB(1,L)=JC(1,L)
   34 JB(2,L)=JC(2,L)
   30 CONTINUE
      DO 35 L=1,3
      JV(1,L,I)=JB(1,L)
      JV(2,L,I)=JB(2,L)
      JA(1,L,NV)=JB(1,L)
   35 JA(2,L,NV)=JB(2,L)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE TSNTNM(NN,NF,NC,SCH,HMN)
      CHARACTER*5 SCHNAM,SCH
      CHARACTER*10 HMNAME,HMN
      IF(NN.LE.0.OR.NN.GT.230) GO TO 4
      REWIND NF
      NC=0
    2 READ(NF,100,END=3) NO,IL,NG,SCHNAM,HMNAME
  100 FORMAT(I3,1X,I2,I2,1X,A5,1X,A10)
      DO 1 I=1,NG
         READ(NF,100)
    1 CONTINUE
      IF(NO.EQ.NN) THEN
         SCH=SCHNAM
         HMN=HMNAME
         NC=NC+1
      ELSE IF(NO.GT.NN) THEN
         RETURN
      END IF
      GO TO 2
    3 IF(NC.NE.0) RETURN
    4 WRITE(6,*) ' YOU MAY GIVE INCORRECT SPACE GROUP NUMBER'
      RETURN
      END    
      SUBROUTINE TSCHTN(SCH,NF,NC,NN,HMN)
      CHARACTER*5 SCHNAM,SCH
      CHARACTER*10 HMNAME,HMN
      REWIND NF
      NC=0
    2 READ(NF,100,END=3) NO,IL,NG,SCHNAM,HMNAME
  100 FORMAT(I3,1X,I2,I2,1X,A5,1X,A10)
      DO 5 I=1,NG
         READ(NF,100)
    5 CONTINUE
      IF(SCH.EQ.SCHNAM) THEN
         NN=NO
         HMN=HMNAME
         NC=NC+1
      ELSE
         IF(NC.NE.0) RETURN
      END IF
      GO TO 2
    3 IF(NC.NE.0.AND.NO.EQ.230) RETURN
      WRITE(6,*) ' YOU MAY GIVE INCORRECT SCHNAME'
      RETURN
      END    
      SUBROUTINE TSHMTN(HMN,NF,NC,NN,SCH)
      CHARACTER*5 SCHNAM,SCH
      CHARACTER*10 HMNAME,HMN
      REWIND NF
      NC=0
    2 READ(NF,100,END=3) NO,IL,NG,SCHNAM,HMNAME
  100 FORMAT(I3,1X,I2,I2,1X,A5,1X,A10)
      DO 5 I=1,NG
         READ(NF,100)
    5 CONTINUE
      IF(HMN.EQ.HMNAME) THEN
         SCH=SCHNAM
         NN=NO
         NC=NC+1
      ELSE
         IF(NC.NE.0) RETURN
      END IF
      GO TO 2
    3 IF(NC.NE.0.AND.NO.EQ.230) RETURN
      WRITE(6,*) ' YOU MAY GIVE INCORRECT INTERNATIONAL NAME'
      RETURN
      END    
      SUBROUTINE TSPNGE(NN,NC,NF)
      CHARACTER*5 SCHNAM
      CHARACTER*10 HMNAME
      INTEGER JB(2,3)
      REWIND NF
      NCHOI=0
    2 READ(NF,100,END=3) NO,IL,NG,SCHNAM,HMNAME
  100 FORMAT(I3,1X,I2,I2,1X,A5,1X,A10)
      IF(NO.EQ.NN) THEN
         NCHOI=NCHOI+1
         IF(NC.EQ.NCHOI) THEN
            WRITE(6,601)
  601       FORMAT(//1X)
            CALL TSPACE(IL)
            WRITE(6,600) NN,SCHNAM,HMNAME,NCHOI
  600       FORMAT(I4,3X,A5,2X,A10,' CHOICE',I2)
            DO 1 K=1,NG
               READ(NF,200) JA,((JB(I,J),I=1,2),J=1,3)
  200          FORMAT(I2,6I2)
               CALL TSGENR(JA,JB)
    1       CONTINUE
            CALL TSPGRP(0)
            CALL TSPGDS
         ELSE
            DO 5 I=1,NG
               READ(NF,200)
    5       CONTINUE 
            GO TO 2   
         END IF
      ELSE
         DO 4 I=1,NG
            READ(NF,200)
    4    CONTINUE
      END IF
      IF(NO.GT.NN) RETURN
      GO TO 2
    3 RETURN
      END
C     SUBROUTINE TSWYRD *====3====*====4====*====5====*====6====*====7   
C     
C     WYCOFF POSITIONS ARE READ FROM FILE NF
C     NSPG: GROUP NUMBER
C     NC  : CHOICE NUMBER
C     NF  : FILE NUMBER
C     NUC : NUMBER OF POSITION
C     NPOS: 1,2, ........,192
C     MPOS: a,b, ........
C     XYZ : REPRESENTATIVE POSITION
C     JB  : 
C     
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C     
      SUBROUTINE TSWYRD(NSPG,NC,NF,NUC,NPOS,MPOS,XYZ,JB)
C     
      INTEGER JB(2,3,30),XYZ(3,30),NPOS(30)
      CHARACTER*1 MPOS(30),SYMB(3,30)
C     
      NCOUT=0
      REWIND NF
   12 READ(NF,100) NN,NUC
  100 FORMAT(2I3)
      IF(NN.NE.NSPG) THEN
         DO 13 K=1,NUC
            READ(NF,200)
   13    CONTINUE
         IF(NCOUT.NE.0) THEN
            WRITE(6,*) ' CHOICE NUMBER IS INCORRECT'
            WRITE(6,*) ' STOP IN TSWYRD'
            STOP
         END IF
         GO TO 12
      ELSE
         NCOUT=NCOUT+1
         DO 14 K=1,NUC
            READ(NF,200) NPOS(K),MPOS(K),(SYMB(I,K),I=1,3)
  200       FORMAT(I2,A1,3(1X,A1))
            IF(NPOS(K).EQ.92) NPOS(K)=192
   14    CONTINUE
      END IF
      IF(NCOUT.LT.NC) GO TO 12
      DO 16 K=1,NUC
         DO 15 I=1,3
            IF(SYMB(I,K).EQ.'0') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'h') THEN
               JB(1,I,K)=1
               JB(2,I,K)=2
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'q') THEN
               JB(1,I,K)=1
               JB(2,I,K)=4
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'t') THEN
               JB(1,I,K)=3
               JB(2,I,K)=4
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'1') THEN
               JB(1,I,K)=1
               JB(2,I,K)=3
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'2') THEN
               JB(1,I,K)=2
               JB(2,I,K)=3
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'3') THEN
               JB(1,I,K)=3
               JB(2,I,K)=8
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'5') THEN
               JB(1,I,K)=5
               JB(2,I,K)=6
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'6') THEN
               JB(1,I,K)=1
               JB(2,I,K)=6
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'7') THEN
               JB(1,I,K)=7
               JB(2,I,K)=8
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'8') THEN
               JB(1,I,K)=1
               JB(2,I,K)=8
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'f') THEN
               JB(1,I,K)=5
               JB(2,I,K)=8
               XYZ(I,K)=0
            ELSE IF(SYMB(I,K).EQ.'x') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=1
            ELSE IF(SYMB(I,K).EQ.'y') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=2
            ELSE IF(SYMB(I,K).EQ.'z') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=3
            ELSE IF(SYMB(I,K).EQ.'-') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=-1
            ELSE IF(SYMB(I,K).EQ.'w') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=-2
            ELSE IF(SYMB(I,K).EQ.'d') THEN
               JB(1,I,K)=0
               JB(2,I,K)=1
               XYZ(I,K)=5
            ELSE IF(SYMB(I,K).EQ.'m') THEN
               JB(1,I,K)=1
               JB(2,I,K)=2
               XYZ(I,K)=-2
            ELSE IF(SYMB(I,K).EQ.'p') THEN
               JB(1,I,K)=1
               JB(2,I,K)=2
               XYZ(I,K)=2
            ELSE IF(SYMB(I,K).EQ.'n') THEN
               JB(1,I,K)=1
               JB(2,I,K)=4
               XYZ(I,K)=-2
            ELSE IF(SYMB(I,K).EQ.'r') THEN
               JB(1,I,K)=1
               JB(2,I,K)=4
               XYZ(I,K)=2
            ELSE IF(SYMB(I,K).EQ.'s') THEN
               JB(1,I,K)=1
               JB(2,I,K)=2
               XYZ(I,K)=1
            ELSE IF(SYMB(I,K).EQ.'u') THEN
               JB(1,I,K)=1
               JB(2,I,K)=2
               XYZ(I,K)=-1
            ELSE IF(SYMB(I,K).EQ.'v') THEN
               JB(1,I,K)=1
               JB(2,I,K)=4
               XYZ(I,K)=1
            END IF
   15    CONTINUE
   16 CONTINUE
      RETURN
      END
C     SUBROUTINE TSWYCF *====3====*====4====*====5====*====6====*====7   
C     
C     
C     
C     
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C     
      SUBROUTINE TSWYCF(A,NUC,NPOS,MPOS,XYZ,JB,NSITE,XYZWM,JAM)
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)                           
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)                             
      SAVE /SPG1/,/SPG2/,IB                                              
      CHARACTER*1 A,MPOS(30)
      INTEGER IB(2,3,6)
      INTEGER XYZ(3,30),JB(2,3,30),NPOS(30)
      INTEGER JA(2,3),XYZW(3),XYZWM(3,48),JAM(2,3,48)
      INTEGER JC(2,3)     
      DATA IB/
     &     1,2, 1,2, 1,2,  0,1, 1,2, 1,2,  1,2, 0,1, 1,2,
     &     1,2, 1,2, 0,1,  2,3, 1,3, 1,3,  1,3, 2,3, 2,3/
      DO 21 I=1,NUC
         IF(MPOS(I).EQ.A) GO TO 22
   21 CONTINUE
      WRITE(6,*) ' POSITION NAME IS INCORRECT'
      WRITE(6,*) ' STOP IN TSWYCF'
      STOP
   22 K=I
      IF(IL.EQ.-1)           NSIT=NPOS(K)/3
      IF(IL.EQ.0.OR.IL.EQ.1) NSIT=NPOS(K)
      IF(IL.EQ.2)            NSIT=NPOS(K)/4
      IF(IL.GE.3)            NSIT=NPOS(K)/2
      NSITE=0
      DO 11 JJ=1,NG
         DO 23 I=1,3
            DO 24 J=1,2
               JA(J,I)=JB(J,I,K)
               JC(J,I)=JV(J,I,JJ)
   24       CONTINUE
   23    CONTINUE
         CALL ZZZY22(IG(JJ),JA)
         CALL ZZZY23(JA,JC,1)
         DO 25 I=1,3
            MA=IT(I,IG(JJ))                            
            IF(IABS(MA).EQ.4) THEN                    
               IF(XYZ(1,K).EQ.1) THEN
                  IF(XYZ(2,K).EQ.0) XYZW(I)=1
                  IF(XYZ(2,K).EQ.2) XYZW(I)=4
                  IF(XYZ(2,K).EQ.-1) XYZW(I)=5
                  IF(XYZ(2,K).EQ.5) XYZW(I)=-1
               ELSE 
                  XYZW(I)=0
               END IF     
               IF(MA.EQ.-4) XYZW(I)=-XYZW(I)
            ELSE                                      
               MB=IABS(MA)                            
               IF(MA.GT.0) XYZW(I)=XYZ(MB,K)
               IF(MA.LT.0) XYZW(I)=-XYZ(MB,K)
            END IF               
   25    CONTINUE               
         IF(NSITE.EQ.0) THEN
            NSITE=NSITE+1
            GO TO 12
         END IF
         DO 13 L=1,NSITE
            DO 14 I=1,3
               IF(XYZWM(I,L).NE.XYZW(I)) GO TO 13
   14       CONTINUE
            DO 43 I=1,3
               IF(JAM(1,I,L).NE.JA(1,I).OR.
     &              JAM(2,I,L).NE.JA(2,I)) GO TO 44
   43       CONTINUE
            GO TO 11
   44       IF(IL.EQ.0.OR.IL.EQ.1) GO TO 13 
            IF(IL.EQ.-1) KS=5
            IF(IL.EQ.2) KS=2
            IF(IL.EQ.3) KS=1
            IF(IL.EQ.4) KS=4
            IF(IL.EQ.-1) KL=6
            IF(IL.EQ.2) KL=4
            IF(IL.EQ.3) KL=1
            IF(IL.EQ.4) KL=4
            DO 45 KK=KS,KL
               DO 46 I=1,3
                  JC(1,I)=IB(1,I,KK)
                  JC(2,I)=IB(2,I,KK)
   46          CONTINUE
               CALL ZZZY23(JC,JA,1)
               DO 47 I=1,3
                  IF(JC(1,I).NE.JAM(1,I,L).OR.
     &                 JC(2,I).NE.JAM(2,I,L)) GO TO 45
   47          CONTINUE
               GO TO 11
   45       CONTINUE
   13    CONTINUE
         NSITE=NSITE+1
   12    DO 2 I=1,3
            XYZWM(I,NSITE)=XYZW(I)
            DO 3 J=1,2
               JAM(J,I,NSITE)=JA(J,I)
    3       CONTINUE
    2    CONTINUE
   11 CONTINUE
      IF(NSIT.NE.NSITE) THEN
         WRITE(6,601) NSITE,NSIT
  601    FORMAT(' NSITE=',I3,' BUT EXPECTED VALUE=',I3)
C         do 61 n=1,nsite
C            write(6,602) (xyzwm(i,n),(jam(j,i,n),j=1,2),i=1,3)
C  602       format(3(i3,i2,'/',i1))
C   61    continue
         STOP
      END IF
      RETURN              
      END
C SUBROUTINE TSPGDS ====*====3====*====4====*====5====*====6====*====7     
C                                                                          
C   PRINT OUT THE SPACE GROUP ELEMENT                                      
C                                                                          
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7     
C                                                                          
      SUBROUTINE TSPGDS                                                    
C                                                                          
        IMPLICIT REAL*8(A-H,O-Z)                                           
C                                                                          
      COMPLEX*16 CR                                                        
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)                                
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)                                  
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)                             
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)                
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)              
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,XYZ,CMN,HMN,AC                           
          
      CHARACTER*4 XYZ(4),HMN(12),CMN(24)
      CHARACTER*1 AC(3)
      CHARACTER A*1,B*4,GN(2,3)*4                                          
      INTEGER KK(3),KA(2,3) 
      DATA XYZ/'X','Y','Z','W'/                                 
      DATA CMN/'E','C2X','C2Y','C2Z',                                    
     &   'C31+','C32+','C33+','C34+',                                      
     &   'C31-','C32-','C33-','C34-',                                      
     &   'C2A','C2B','C2C','C2D','C2E','C2F',                              
     &   'C4X+','C4Y+','C4Z+','C4X-','C4Y-','C4Z-'/                        
      DATA HMN/                                                 
     &  'E','C6+','C3+','C2','C3-','C6-',                                  
     &  'C211','C221','C231','C212','C222','C232'/                         
      DATA AC/' ','I','-'/                                                      
                                
      MNG=NG                                                               
      IPS=1                                                                
      IF(IL.EQ.-1) WRITE(6,*) ' TRIGONAL LATTICE'
      IF(IL.EQ.0)  WRITE(6,*) ' HEXAGONAL LATTICE'
      IF(IL.EQ.1)  WRITE(6,*) ' SIMPLE LATTICE'
      IF(IL.EQ.2)  WRITE(6,*) ' FACE CENTERED LATTICE'
      IF(IL.EQ.3)  WRITE(6,*) ' BODY CENTERED LATTICE'
      IF(IL.EQ.4)  WRITE(6,*) ' C- CENTERED LATTICE'
      WRITE(6,602)                                                         
  602 FORMAT('   GROUP ELEMENTS')                                         
      GO TO 11                                                             
C ENTRY TSTRDS ====2====*====3====*====4====*====5====*====6====*====7     
C   TIME REVERSAL ELEMENTS                                                 
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7     
      ENTRY TSTRDS                                                         
      IPS=3                                                                
      MNG=MTRG                                                             
      WRITE(6,611)                                                         
  611 FORMAT(//' TIME REVERSAL ELEMENT')                                   
      IF(MTRG.GT.0) GO TO 11                                               
      WRITE(6,612)                                                         
  612 FORMAT(' THERE IS NO TIME REVERSAL ELEMENT')                         
      RETURN                                                               
C ENTRY TSOPDS ====2====*====3====*====4====*====5====*====6====*====7     
C  TABLE OF OPERATION CODE AND MULTIPLICATION TABLE                        
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7     
      ENTRY TSOPDS                                                         
      IPS=0                                                                
      WRITE(6,605)                                                         
  605 FORMAT('   TABLE OF OPERATION CODE')                                
      IF(IL.LE.0) MNG=24                                                   
      IF(IL.GT.0) MNG=48                                                   
      GO TO 11                                                             
C ENTRY TSPKDS ====2====*====3====*====4====*====5====*====6====*====7     
C     GROUP OF POINT K                                                     
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7     
      ENTRY TSPKDS                                                         
      MNG=MG                                                               
      IPS=2                                                                
      WRITE(6,603) KB,ICB                                                  
  603 FORMAT('   ELEMENTS OF PK-GROUP K=',3I3,1H/,I3)                     
      IF(IDOUB.EQ.1) WRITE(6,604)                                          
  604 FORMAT(' DOUBLE REPRESENTATION')                                     
   11 CONTINUE                                                             
      DO 1 I=1,MNG                                                         
      IF(IPS.EQ.0) ITA=I                                                   
      IF(IPS.EQ.1) ITA=IG(I)                                               
      IF(IPS.EQ.2) ITA=IG(JG(I))                                           
      IF(IPS.EQ.3) ITA=IG(JTRG(4,I))                                       
      A=AC(1)                                                              
      IF((IL.LE.0.AND.ITA.GT.12).OR.                                       
     &   (IL.GT.0.AND.ITA.GT.24))  A=AC(2)                                 
      IF(IL.LE.0) B=HMN(MOD(ITA-1,12)+1)                                   
      IF(IL.GT.0) B=CMN(MOD(ITA-1,24)+1)                                   
      DO 2 J=1,3                                                           
      MA=IT(J,ITA)                                                         
      IF((IPS.EQ.2.OR.IPS.EQ.3).AND.IL.LE.0)                               
     &   MA=IT(J,ITA+24)                                                   
      MB=IABS(MA)                                                          
      IF(MA.GT.0) GN(1,J)=AC(1)                                            
      IF(MA.LT.0) GN(1,J)=AC(3)                                            
      GN(2,J)=XYZ(MB)                                                      
    2 CONTINUE                                                             
      IF(IPS.EQ.0) GO TO 13                                                
      IF(IPS.EQ.1) JGA=I                                                   
      IF(IPS.EQ.2) JGA=JG(I)                                               
      IF(IPS.EQ.3) JGA=JTRG(4,I)                                           
      DO 3 J=1,3                                                           
      KA(1,J)=JV(1,J,JGA)                                                  
    3 KA(2,J)=JV(2,J,JGA)                                                  
      IF(IPS.EQ.1) GO TO 12                                                
      IF(IPS.EQ.2) GO TO 5                                                 
      DO 6 J=1,3                                                           
    6 KK(J)=JTRG(J,I)                                                      
      IITR=ITRC(1,I)                                                       
      IIWW=ITRC(2,I)                                                       
      JJWW=ITRC(3,I)                                                       
      WRITE(6,616) I,ITA,A,B,GN,KA,KK,IITR,IIWW,JJWW                       
  616 FORMAT(I3,I4,2X,A1,A4,3(1X,2A1),3(I3,1H/,I1),3I3,                    
     &       I4,2H (,I2,1H/,I2,1H))                                        
      GO TO 1                                                              
    5 DO 4 J=1,3                                                           
    4 KK(J)=JK(J,I)                                                        
      WRITE(6,601) I,ITA,A,B,GN,KA,KK                                      
  601 FORMAT(I3,I4,2X,A1,A4,3(1X,2A1),                                     
     &   3(I3,1H/,I1),3I3)                                                 
      GO TO 1                                                              
   12 WRITE(6,601) I,ITA,A,B,GN,KA                                         
      GO TO 1                                                              
   13 WRITE(6,606) ITA,A,B,GN                                              
  606 FORMAT(I4,2X,A1,A4,3(1X,2A1))                                        
    1 CONTINUE                                                             
      IF(IPS.EQ.0) GO TO 9                                                 
      IF(IPS.NE.3) RETURN                                                  
      WRITE(6,613) (I,I=1,NR)                                              
      WRITE(6,614) (IATR(I),I=1,NR)                                        
  613 FORMAT(//' TIME REVERSAL SUM'/12I4)                                  
  614 FORMAT(12I4)                                                         
      RETURN                                                               
    9 IF(IL.LE.0) WRITE(6,663)                                             
  663 FORMAT(1H //' GROUP TABLE FOR D6H')                                  
      IF(IL.GT.0) WRITE(6,662)                                             
  662 FORMAT(1H //' GROUP TABLE FOR O-GROUP')                              
      DO 7 I=1,24                                                          
    7 WRITE(6,660) (IM(I,J),J=1,24)                                        
  660 FORMAT(24I3)                                                         
      WRITE(6,661)                                                         
  661 FORMAT(//' INVERS ELEMENTS')                                         
      IF(IL.LE.0) NN=24                                                    
      IF(IL.GT.0) NN=48                                                    
      WRITE(6,660) (IV(I),I=1,NN)                                          
      RETURN                                                               
      END                                                                  
C SUBROUTINE TSSLPW ====*====3====*====4====*====5====*====6====*====7
C
C      SYMMETRIZED PLANE WAVE.
C
C      1983.4.11. :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSSLPW(AM,KB,IC,IR,NA,NDI,KO,V,B,D,IND,NK,ND1,ND2
     &                  ,NNWW)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER(MAXNPW=4854)
      COMPLEX*16 V,U,WD,CR,CW,WW
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      SAVE /SPW/,KKBB,ICIC
      DIMENSION KO(3,1),V(1),B(1),IND(3,1),KB(3),D(1)
      DIMENSION JGA(2,3,48),WD(48),JG(48),CR(48),KK(3),NIP(100)
      DIMENSION KKBB(3)
      DATA KKBB/0, 0, 0/
      DATA ICIC/0/
      CALL CHKDNM(IC)
      PAI2=8.D0*DATAN(1.D0)
      IF(ICIC.NE.IC) GO TO 11
      DO 12 K=1,3
      IF(KKBB(K).NE.KB(K)) GO TO 11
   12 CONTINUE
      DO 13 I=1,IK
   13 KM(4,I)=0
      GO TO 15
   11 CONTINUE
      CALL ZZZY41(KB,IC,AM)
      ICIC=IC
      DO 14 K=1,3
   14 KKBB(K)=KB(K)
   15 CONTINUE
C      WRITE(6,690) IE,IK
C  690 FORMAT(10I5)
      NK=IK
      CALL TSIRMR(IR,NA,NA,MG,ND,JG,JGA,CR,WD)
C      WRITE(6,690) ND,MG
      NDI=ND
      JV=0
      ND1=0
      ND2=0
      IW=0
      DO 21 ID=1,IE
      NN=ND1+1
      JW=JV+1
      JV=KT(ID)
C
C     CHARACTER TEST
C
      CW=0.0D0
      DO 22 II=1,MG
      IF(CDABS(CR(II)).LT.1.0D-4) GO TO 22
      JA=JG(II)
      WW=0.0D0
      DO 23 JU=JW,JV
      JJ=JU
      CALL ZZZY40(JA,JJ,KK)
      WA=0.0D0
      DO 24 K=1,3
      IF(KK(K).NE.KM(K,JU)) GO TO 23
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
      WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
   24 CONTINUE
      X=PAI2*WA
      WW=WW+DCMPLX(COS(X),SIN(X))
   23 CONTINUE
      CW=CW+CONJG(CR(II))*WW
   22 CONTINUE
      CW=CW/MG
      NI=CW+0.5
C      WRITE(6,690) ID,JW,JV,NI,nd1,nd2,iw
      IF(NI.EQ.0) GO TO 21
      IP=JW
      NP=0
      NB=NA
      CALL TSIRME(NA,NA,WD)
      INDIC=0
C
C     PROJECTION OPERATOR
C
   47 DO 27 JU=JW,JV
   27 U(JU)=0.0D0
      DO 28 II=1,MG
      IF(CDABS(WD(II)).LT.1.0D-4) GO TO 28
      JA=JG(II)
      CALL ZZZY40(JA,IP,KK)
      WA=0.0D0
      DO 29 K=1,3
      KD=JGA(1,K,II)
      KC=JGA(2,K,II)
   29 WA=WA+DBLE((KB(K)-KK(K))*KD)/DBLE(IC*KC)
      X=PAI2*WA
      WW=DCMPLX(COS(X),SIN(X))
      DO 30 JU=JW,JV
      DO 20 K=1,3
      IF(KM(K,JU).NE.KK(K)) GO TO 30
   20 CONTINUE
      KM(4,JU)=IP
      U(JU)=U(JU)+WW*CONJG(WD(II))
      GO TO 28
   30 CONTINUE
      WRITE(6,602)
  602 FORMAT(' STOP AT 30 IN TSSLPW')
      STOP
   28 CONTINUE
C
C     ORTHOGONALITY TEST
C
C     WRITE(6,690) NN,ND1,ND2,NA,NB,IP
      IF(INDIC.EQ.2) GO TO 31
      IF(NN.GT.ND1) GO TO 31
      DO 32 II=NN,ND1
      ID1=IND(1,II)
      ID2=IND(2,II)
      WW=0.0D0
      DO 33 I=ID1,ID2
      JU=JW
   36 DO 34 K=1,3
      IF(KO(K,I).NE.KM(K,JU)) GO TO 35
   34 CONTINUE
      WW=WW+CONJG(V(I))*U(JU)
      GO TO 33
   35 JU=JU+1
      IF(JU.LE.JV) GO TO 36
      WRITE(6,603)
  603 FORMAT(' STOP AT 36 IN TSSLPW')
      WRITE(6,*) NN,ND1,ND2,NA,NB,IP
      write(6,*) iw,ie,ik,ic,ii,id1,id2,nn,nd1
      write(6,*) i,(ko(k,i),k=1,3)
     &     ,(km(k,jw),k=1,3),(km(k,jv),k=1,3)
      STOP
   33 CONTINUE
C      write(6,*) ii,ip,ww,nd1,nd2,na,nb,indic
      IF(CDABS(WW).GT.1.0D-4) GO TO 37
   32 CONTINUE
C
C     REGISTRATION
C
   31 ND1=ND1+1
      IND(1,ND1)=IW+1
      WA=0.0D0
      IIWW=IW+1
      DO 38 JU=JW,JV
      IF(CDABS(U(JU)).GT.1.0D-4) GO TO 381
      IF(JU.EQ.IP) GO TO 381
      IF(KM(4,JU).EQ.IP) KM(4,JU)=0
      GO TO 38
  381 IW=IW+1
      KO(1,IW)=KM(1,JU)
      KO(2,IW)=KM(2,JU)
      KO(3,IW)=KM(3,JU)
      V(IW)=U(JU)*(DBLE(ND)/DBLE(MG))
      IF(ABS(DIMAG(V(IW))).LT.1.0D-4) V(IW)=DBLE(V(IW))
      IF(ABS(DBLE(V(IW))).LT.1.0D-4) V(IW)=DCMPLX(0.0D0,DIMAG(V(IW)))
      WA=WA+DBLE(V(IW)*CONJG(V(IW)))
      IF(JU.NE.IP) GO TO 38
      IF(IW .EQ.IIWW) GO TO 38
      DO 382 K=1,3
      KKW=KO(K,IW)
      KO(K,IW)=KO(K,IIWW)
  382 KO(K,IIWW)=KKW
      WW=V(IW)
      V(IW)=V(IIWW)
      V(IIWW)=WW
   38 CONTINUE
      IND(2,ND1)=IW
      B(ND1)=SQRT(WA)
      IF(WA.GT.1.0D-4) GO TO 39
      IW=IIWW-1
      ND1=ND1-1
C      write(6,*) ' 0-vector',iiww,nd1+1
      GO TO 37
   39 IF(INDIC.EQ.1) GO TO 40
      IF(INDIC.EQ.2) GO TO 44
      NP=NP+1
      NIP(NP)=IP
      IND(3,ND1)=0
      D(ND1)=A(IP)
      NI=NI-1
      IF(NI.LE.0) GO TO 21
      GO TO 25
   40 IND(3,ND1)=ND1+1
      D(ND1)=A(IP)
      INDIC=2
      ND2=ND2+1
C      write(6,*) 'nd2+1',nd2
      CALL TSIRME(NB,NB,WD)
      GO TO 47
   44 IND(3,ND1)=-(ND1-1)
      D(ND1)=A(IP)
      NI=NI-1
      IF(NI.LE.0) GO TO 21
      GO TO 48
C
C
C
   37 IF(NA.NE.NB) GO TO 48
   25 JJWW=IP+1
C      IP=IP+1
C      IF(IP.GT.JV) GO TO 41
      DO 41 JU=JJWW,JV
      IF(KM(4,JU).NE.0) GO TO 41
      IP=JU
      GO TO 47
   41 CONTINUE
      NNP=NP
      IF(ND.NE.2.AND.ND.NE.4) GO TO 46
      IF(ND.EQ.2) NB=3-NA
      IF(ND.EQ.4) NB=5-NA
C     WRITE(6,690) NNP,NP,NA,NB
      NP=0
   48 NP=NP+1
      IF(NP.GT.NNP) GO TO 46
      IP=NIP(NP)
      INDIC=1
      CALL TSIRME(NA,NB,WD)
      GO TO 47
   46 WRITE(6,604)
  604 FORMAT(' STOP AT 46 IN TSSLPW')
      STOP
   21 CONTINUE
      NW=IW
      NNWW=IW
      RETURN
      ENTRY TSSWDS
      WRITE(6,663) IR,NA,NDI,NK,ND1,ND2,ND1-ND2
  663 FORMAT(//7I5/)
      IF(NW.EQ.0) GO TO 492
      ID1=0
      DO 49 IW=1,NW
      IF(IW.NE.IND(1,ID1+1).OR.ID1.GE.ND1) GO TO 491
      ID1=ID1+1
      WRITE(6,661) ID1,(IND(K,ID1),K=1,3),D(ID1),B(ID1)
  661 FORMAT(I4,3I5,2F10.5)
      ICOUNT=0
  491 ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KO(K,IW),K=1,3),IC,V(IW)
  660 FORMAT(4X,2I5,2H (,3I5,2H)/,I3,3H  (,2F10.5,1H))
   49 CONTINUE
      RETURN
  492 WRITE(6,662)
  662 FORMAT(' NO STATE FOR THIS REPRESENTATION')
      RETURN
      END
      SUBROUTINE LINEA3(A,B,IPV,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(3,3),B(3),IPV(3),D(3),X(3)
      N3=3
      DO 10 K=1,N3
         IPV(K)=K
   10 CONTINUE    
      DO 1 K=1,N3-1
         DO 2 L=K,N3
            IF(DABS(A(IPV(L),K)).GT.1.0D-4) GO TO 3
    2    CONTINUE
         IERR=1
         RETURN       
    3    IF(K.NE.L) THEN
            IW=IPV(K)
            IPV(K)=IPV(L)
            IPV(L)=IW
         END IF
         D(IPV(K))=1.0D0/A(IPV(K),K)
         DO 4 I=K+1,N3
            WW=A(IPV(I),K)*D(IPV(K))
            DO 5 J=k+1,N3
               A(IPV(I),J)=A(IPV(I),J)-WW*A(IPV(K),J)
    5       CONTINUE   
            B(IPV(I))=B(IPV(I))-WW*B(IPV(K))
    4   CONTINUE
    1 CONTINUE
      IF(DABS(A(IPV(N3),N3)).GT.1.0D-4) GO TO 6
      IERR=1
      RETURN
    6 D(IPV(N3))=1.0D0/A(IPV(N3),N3)
      X(N3)=B(IPV(N3))*D(IPV(N3))
      DO 7 KK=1,N3-1
         K=N3-KK
         WW=0.0D0
         DO 8 I=K+1,N3
            WW=WW+A(IPV(K),I)*X(I)
    8    CONTINUE   
         X(K)=(B(IPV(K))-WW)*D(IPV(K))
    7    CONTINUE
      DO 9 K=1,N3
         B(K)=X(K)
    9 CONTINUE
      IERR=0
      RETURN
      END   
C SUBROUTINE TSLATC ====*====3====*====4====*====5====*====6====*====7
C
C  LATTICE PARAMETERS ARE STORED IN COMMON/LAT/
C  IF THE SYMMETRY IS CUBIC,      THIS PROGRAM MAKES A=B=C AND CA=CB=CC=0.0
C  IF THE SYMMETRY IS TETRAGONAL, THIS PROGRAM MAKES A=B AND CA=CB=CC=0.0
C  IF THE SYMMETRY IS ORTHRHOMBIC, THIS PROGRAM MAKES CA=CB=CC=0.0
C  IF THE SYMMETRY IS MONOCLINIC, CA=CB=0.0, CB=CC=0.0 OR CC=CA=0.0
C  IF THE SYMMETRY IS HEXAOGONAL OR TRIGONAL,
C                    THIS PROGRAM MAKES A=B, CA=CB=0.0 AND CC=-0.5.
C    IN SPITE OF THE PARAMETERS GIVEN BY AX,BX, E.T.C HAVE NOT THESE
C    SYMMETRY.
C
C               MODIFIED ON 1988/10/21
C                 ENTRY GETRCP TO GET THE PARAMETERS FOR RECIPROCAL
C                 LATTICE IS MADE IN THIS TIME.
C                  BY A. YANASE
C               MODIFIED ON 1993/12/24
C                 ENTRY GETKVC IS MADE
C                   BY A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSLATC(AX,BX,CX,CAX,CBX,CCX)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      LOGICAL LA,LB,LC,LCUB,LTET
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/LAT/A,B,C,CA,CB,CC,A1,B1,C1
      SAVE /LAT/,/SPG2/,AAK,BBK,CCK,BCK,CAK,ABK,W1,W2,W3,W4,W5
     &   ,AS,BS,CS,COAK,COBK,COCK,X1,X2,X3,X4,X5
      DIMENSION V(3),DC(3)
      IF(AX.LT.1.0D-5.OR.BX.LT.1.0D-5.OR.CX.LT.1.0D-5) THEN
        WRITE(6,*) ' LATTICE CONSTANTS SHOULD BE FINITE POSITIVIE'
        WRITE(6,*) ' STOP IN TSLATC'
        STOP
      END IF
      IF(ABS(CAX).GT.1.0D0.OR.ABS(CBX).GT.1.0D0
     &   .OR.ABS(CCX).GT.1.0D0) THEN
        WRITE(6,*) ' THE ABSOLUTE VALUE OF COS SHOULD BE LESS THAN 1'
        WRITE(6,*) ' STOP IN TSLATC'
        STOP
      END IF
      IF(AX*ABS(CCX).GT.BX.OR.BX*ABS(CCX).GT.AX.OR.
     &   BX*ABS(CAX).GT.BX.OR.CX*ABS(CAX).GT.BX.OR.
     &   CX*ABS(CBX).GT.AX.OR.AX*ABS(CBX).GT.CX) THEN
        WRITE(6,*) ' YOUR LATTICE HAS BETTER CHOICE OF LAT. CONST.'
        WRITE(6,*) ' STOP IN TSLATC'
        STOP
      END IF  
      A=AX
      B=BX
      C=CX
      CA=CAX
      CB=CBX
      CC=CCX
      IF(IL.GE.1) GO TO 1
      B=A
      CC=-0.5D0
      CA=0.0D0
      CB=0.0D0
      GO TO 2
    1 CONTINUE
      LA=.FALSE.
      LB=.FALSE.
      LC=.FALSE.
      LCUB=.FALSE.
      LTET=.FALSE.
      DO 3 I=1,NG
      JJ=IG(I)
      IF(JJ.EQ.4) LC=.TRUE.
      IF(JJ.EQ.28) LC=.TRUE.
      IF(JJ.EQ.2) LA=.TRUE.
      IF(JJ.EQ.26) LA=.TRUE.
      IF(JJ.EQ.27) LB=.TRUE.
      IF(JJ.EQ.3) LB=.TRUE.
      IF(JJ.EQ.5) LCUB=.TRUE.
      IF(JJ.EQ.21) LTET=.TRUE.
      IF(JJ.EQ.48) LTET=.TRUE.
    3 CONTINUE
      IF(LC) CA=0.0D0
      IF(LC) CB=0.0D0
      IF(LA) CB=0.0D0
      IF(LA) CC=0.0D0
      IF(LB) CC=0.0D0
	  IF(LTET) CC=0.0D0
C    MODIFIED AT 2006/12/11 by A.YANASE	  
      IF(LB) CA=0.0D0
      IF(LTET) B=A
      IF(LCUB) B=A
      IF(LCUB) C=A
    2 CONTINUE
      WRITE(6,600) A,B,C,CA,CB,CC
  600 FORMAT(' LATTICE CONSTANTS ARE SET AS'
     &      /'  A=',F9.5,'  B=',F9.5,'  C=',F9.5
     &      /' CA=',F9.5,' CB=',F9.5,' CC=',F9.5)
      CAA=CA*CA
      CBB=CB*CB
      CCC=CC*CC
      CW=1.0D0-CAA-CBB-CCC+2.0D0*CA*CB*CC
      IF(CW.LT.1.0D-5) THEN
         WRITE(6,*) ' THE ANGLES BETWEEN AXIS CAN NOT GIVE ANY LATTICE'
         WRITE(6,*) ' STOP IN TSLATC'
         STOP
      END IF 
      A2=(1.0D0-CAA)/CW
      B2=(1.0D0-CBB)/CW
      C2=(1.0D0-CCC)/CW
      IF(A2.LT.1.0D-5.OR.B2.LT.1.0D-5.OR.C2.LT.1.0D-5) THEN
         WRITE(6,*) ' THE ANGLES BETWEEN AXIS CAN NOT GIVE ANY LATTICE'
         WRITE(6,*) ' STOP IN TSLATC'
         STOP
      END IF 
      A1=SQRT(A2)
      B1=SQRT(B2)
      C1=SQRT(C2)
      AAK=1.0D0
      BBK=((1.0D0-CBB)/(1.0D0-CAA))*((A/B)**2)
      CCK=((1.0D0-CCC)/(1.0D0-CAA))*((A/C)**2)
      BCK=2.0D0*((CB*CC-CA)/(1.0D0-CAA))*((A*A)/(C*B))
      CAK=2.0D0*((CC*CA-CB)/(1.0D0-CAA))*(A/C)
      ABK=2.0D0*((CA*CB-CC)/(1.0D0-CAA))*(A/B)
      SC=SQRT(1.0D0-CCC)
      W1=(B/A)*CC
      W2=(C/A)*CB
      W3=(B/A)*SC
      W4=(C/A)*((CA-CB*CC)/SC)
      W5=(C/A)*(SQRT(CW)/SC)
      AS=SQRT(AAK)
      BS=SQRT(BBK)
      CS=SQRT(CCK)
      COAK=(0.5*BCK)/(BS*CS)
      COBK=(0.5*CAK)/(CS*AS)
      COCK=(0.5*ABK)/(AS*BS)
      CWK=1.0D0-COAK*COAK-COBK*COBK-COCK*COCK
     &      +2.0D0*COAK*COBK*COCK
      SCK=SQRT(1.0D0-COCK*COCK)
      X1=(BS/AS)*COCK
      X2=(CS/AS)*COBK
      X3=(BS/AS)*SCK
      X4=(CS/AS)*((COAK-COBK*COCK)/SCK)
      X5=(CS/AS)*(SQRT(CWK)/SCK)
      RETURN
C ENTRY GETKVC ====2====*====3====*====4====*====5====*====6====*====7
C     XK,YK,ZK ARE IN RECTANGLE COORDINATE
C     XXK,YYK,ZZK ARE IN A*,B*,C* COORDINATE 
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY GETKVC(XK,YK,ZK,XXK,YYK,ZZK)
      ZZK=ZK/X5
      YYK=(YK-ZZK*X4)/X3
      XXK=XK-YYK*X1-ZZK*X2
      RETURN
C ENTRY ZZZY37 ====2====*====3====*====4====*====5====*====6====*====7
C    S=LENGTH OF RECIPROCAL VECTOR(KX/IA,KY/IA,KZ/IA) IN THE UNIT OF A*
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY ZZZY37(KX,KY,KZ,IA,S)
      X=DBLE(KX)/DBLE(IA)
      Y=DBLE(KY)/DBLE(IA)
      Z=DBLE(KZ)/DBLE(IA)
      S=SQRT(X*X*AAK+Y*Y*BBK+Z*Z*CCK
     &      +X*Y*ABK+Y*Z*BCK+Z*X*CAK)
      RETURN
C ENTRY ZZZY43 ====2====*====3====*====4====*====5====*====6====*====7
C   S=LEMGTH OF VECTOR V IN UNIT OF A
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY ZZZY43(V,DC,S)
      DC(1)=V(1)+V(2)*W1+V(3)*W2
      DC(2)=V(2)*W3+V(3)*W4
      DC(3)=V(3)*W5
      S=SQRT(DC(1)*DC(1)+DC(2)*DC(2)+DC(3)*DC(3))
      DC(1)=DC(1)/S
      DC(2)=DC(2)/S
      DC(3)=DC(3)/S
      RETURN
C ENTRY GETRCP ====2====*====3====*====4====*====5====*====6====*====7
C    LATTICE PARAMETERS OF THE RECIPROCAL LATTICE ARE GIVEN HERE
C    IN THE UNIT OF A*. THE VALUE OF ASTAR IS ALWAYS 1.0D0 NOW.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY GETRCP(ASTAR,BSTAR,CSTAR,COSAK,COSBK,COSCK)
      ASTAR=AS
      BSTAR=BS
      CSTAR=CS
      COSAK=COAK
      COSBK=COBK
      COSCK=COCK
      END
C SUBROUTINE TSCRFR ====*====3====*====4====*====5====*====6====*====7
C
C    THE CALCULATION FOR THE CHARACTER CCW OF FULL REPRESENTATION
C    OF SPACE GROUP ELEMENT WITH LATTICE TRANSRATION OF LAT(2,3)
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSCRFR(LAT,CCW)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW,WW,SN,CC
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/SPG6/SN(2,2,24),IDF(24,24)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/STK/,/SPG6/
      COMPLEX*16 CCW(48,12)
      DIMENSION LAT(2,3)
      REAL*8 IW,JW
      CALL ZZZY38
      DO 3 I=1,NR
      DO 3 J=1,NG
    3 CCW(J,I)=(0.D0,0.D0)
      DO 4 I=1,NG
      DO 5 J=1,NS
      KW=IM(IV(IG(JS(J))),IM(IG(I),IG(JS(J))))
      DO 6 K=1,MG
      IF(KW.EQ.IG(JG(K))) GO TO 7
    6 CONTINUE
      GO TO 5
    7 IA=IG(JS(J))
C--------------------------------------
C      DOUBLE GROUP FACTOR   1993/12/28
C--------------------------------------
      IEXTF=1
      IF(IDOUB.EQ.0) GO TO 11
      IB=IG(I)
      IC=IV(IA)
      ID=IA
      IE=IM(IB,ID) 
      IF(ID.GT.24) ID=ID-24
      IF(IB.GT.24) IB=IB-24
      IF(IC.GT.24) IC=IC-24
      IF(IE.GT.24) IE=IE-24
      IEXTF=IDF(IC,IE)*IDF(IB,ID)*IDF(IC,ID)
   11 CONTINUE
      IF(IL.LE.0) IA=IA+24
      IW=0.D0
      JW=1.D0
      DO 8 L=1,3
      CALL CHKDNM(LAT(2,L))
      IF(JV(1,L,JS(J)).EQ.0) GO TO 9
      MA=IT(L,IA)
      IF(IABS(MA).EQ.4) KW=JK(1,K)+JK(2,K)
      IF(IABS(MA).NE.4) KW=JK(IABS(MA),K)
      IF(MA.LT.0) KW=-KW
      IW=IW*JV(2,L,JS(J))+JW*KW*JV(1,L,JS(J))
      JW=JW*JV(2,L,JS(J))
    9 IW=IW*JV(2,L,I)*ICB-JW*JV(1,L,I)*KS(L,J)
      JW=JW*JV(2,L,I)*ICB
      IW=IW*LAT(2,L)*ICB-JW*LAT(1,L)*KS(L,J)
      JW=JW*LAT(2,L)*ICB
    8 CONTINUE
      W=2.0D0*3.14159265358979312D0*IW/JW
      WW=DCMPLX(COS(W),SIN(W))
      DO 10 L=1,NR
      CC=CR(K,L)
      IF(IEXTF.LT.0) CC=-CC
      CCW(I,L)=CCW(I,L)+CC*WW
   10 CONTINUE
C      write(6,*) i,j,k,(ccw(i,l),l=1,nr)
    5 CONTINUE
C
      DO 1 L=1,NR
      IF(ABS(DBLE(CCW(I,L))).LT.1.0D-4)
     &              CCW(I,L)=DCMPLX(0.D0,DIMAG(CCW(I,L)))
    1 IF(ABS(DIMAG(CCW(I,L))).LT.1.0D-4)
     &              CCW(I,L)=DCMPLX(DBLE(CCW(I,L)),0.D0)
C
    4 CONTINUE
      RETURN
      END
C SUBROUTINE FIDPRT ====*====3====*====4====*====5====*====6====*====7
C
C    WHEN THE TIME REVERSAL SYMMETRY IS PAIRING, THIS ROUTINE FIND THE
C    PARTNER BY THE CALCULATION FOR THE CHARACTER OF FULL REPRESENTATION
C    OF SPACE GROUP
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE FIDPRT(IPART)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/STK/
      COMPLEX*16 CCW(48,12),CWW(48,8,12)
      DIMENSION IPART(12),LAT(2,3),LATG(3)
      DIMENSION LATT(2,3,8)
C-----------------------------------------
C   STAR OF KB
C
      CALL ZZZY35
C------------------------------------------
C   TIME REVERSAL
C
      CALL ZZZY38
C-------------------------------------------
      ICD=0
      DO 100 NREP=1,8
      CALL GENLTP(LATG,ICD)
      LAT(1,1)=LATG(1)
      LAT(1,2)=LATG(2)
      LAT(1,3)=LATG(3)
      LAT(2,1)=ICD
      LAT(2,2)=ICD
      LAT(2,3)=ICD
      CALL TSCRFR(LAT,CCW)
      DO 37 I=1,3
      DO 37 J=1,2
      LATT(J,I,NREP)=LAT(J,I)
   37 CONTINUE
      DO 31 I=1,NG
      DO 31 J=1,NR
      CWW(I,NREP,J)=CCW(I,J)
   31 CONTINUE
Cc      DO 42 LT=1,NREP
Cc      DO 41 I=1,NG
Cc          DO 43 J=1,NR
Cc             IF(ABS(DBLE(CWW(I,LT,J))).GT.1.0D-5) GO TO 44
Cc             IF(ABS(DIMAG(CWW(I,LT,J))).GT.1.0D-5) GO TO 44
Cc   43     CONTINUE
Cc          GO TO 41
Cc   44 CONTINUE
Cc      WRITE(6,603) LT,IG(I),(LATT(1,K,LT),K=1,3),LATT(2,1,LT)
Cc     &             ,(DBLE(CWW(I,LT,J)),J=1,NR)
Cc      WRITE(6,604) (DIMAG(CWW(I,LT,J)),J=1,NR)
Cc  603 FORMAT(2I3,3I2,'/',I1,12F7.4)
Cc  604 FORMAT(14X,12F7.4)
Cc   41 CONTINUE
Cc   42 CONTINUE
C
      DO 34 L=1,NR
         IPART(L)=0
   34 CONTINUE
      INDR=0
      DO 32 L=1,NR
      IF(IATR(L).NE.0.OR.IPART(L).NE.0) GO TO 32
        INDC=0
        IF(L.EQ.NR) GO TO 36
        DO 33 LL=L+1,NR
        IF(IATR(LL).NE.0.OR.IPART(LL).NE.0) GO TO 33
        DO 35 LT=1,NREP
        DO 35 J=1,NG
        IF(ABS(DBLE(CWW(J,LT,L))-DBLE(CWW(J,LT,LL))).GT.1.0D-5
     &  .OR.ABS(DIMAG(CWW(J,LT,L))+DIMAG(CWW(J,LT,LL))).GT.1.0D-5)
     &      GO TO 33
   35   CONTINUE
        IF(INDC.EQ.0) THEN
           IPART(L)=LL
           IPART(LL)=L
           INDC=1
        ELSE
C          WRITE(6,*) ' PARTNER IS NOT UNIQUE '
           INDR=1
           GO TO 32
        END IF
   33   CONTINUE
           IF(INDC.EQ.1) GO TO 32
   36      WRITE(6,*) ' PARTNER IS NOT FOUND '
Cc           return
           STOP
   32 CONTINUE
      IF(INDR.EQ.0) GO TO 101
  100 CONTINUE
           WRITE(6,*) ' PARTNER IS NOT UNIQUELY DETERMINED '
           STOP
  101 CONTINUE
C     WRITE(6,*) (I,I=1,NR)
C     WRITE(6,*) (IPART(I),I=1,NR)
      RETURN
      END
C SUBROUTINE TSNBTB ====*====3====*====4====*====5====*====6====*====7
C
C      NEIGHBORS OF A SPACE LATTICE.
C
C      1983.4.11. :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSNBTB(RMRM,JF)
C     NEIGHBOURS OF A SPACE LATTICE
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
C     PARAMETER (LMNATM=50,LMNKAT=10)
      COMPLEX*16 WWST,V
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/LAT/A,B,C,CA,CB,CC,A1,B1,C1
      COMMON/ATT/ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &       ,NKATOM,NATOM,KS(11,LMNKAT),IA(LMNATM),JRCH
      COMMON/BON/KBO(2,30),BOL(30),NBO,KBOND(4,200),NBOND,WWST(4000)
     &      ,IVEC(10000),VST(3,4000),NVEC,NDIM
     &      ,KA(2,1500),V(1500),INS(3,500),AA(500),NT,NS
      SAVE /SPG2/,/LAT/,/ATT/,/BON/
      DIMENSION XS(3),XA(8,201)
      RM=RMRM
      ITSBO=0
      GO TO 41
      ENTRY TSBOVR(JF)
      NAB=NATOM*NATOM
      IF(NBO*NAB.LE.10000) GO TO 91
      WRITE(6,901)
  901 FORMAT(' STOP AT 91 IN TSBOVR')
      STOP
   91 CONTINUE
      NBONAB=NBO*NAB
      DO 94 I=1,NBONAB
   94 IVEC(I)=0
      ITSBO=1
      RM=0.0D0
      DO 42 I=1,NBO
      IF(RM.LT.BOL(I)) RM=BOL(I)
   42 CONTINUE
      NVEC=0
      RM=RM+1.0D-4
   41 CONTINUE
      P=B/A
      Q=C/A
      XM=RM*A1
      IW=XM
      NX=2*IW+3
      MX=IW+2
      YM=RM*B1*(A/B)
      IW=YM
      NY=2*IW+3
      MY=IW+2
      ZM=RM*C1*(A/C)
      IW=ZM
      NZ=2*IW+3
      MZ=IW+2
      DO 11 I1=1,NATOM
      NN=0
      DO 12 IS=1,NATOM
      K=1
      IF(IL.EQ.-1) K=3
      IF(IL.EQ.2) K=4
      IF(IL.EQ.3) K=2
      IF(IL.EQ.4) K=2
      DO 13 J=1,K
      DO 14 I=1,3
      W=VATOM(I,IS)
      IF(J.EQ.1) GO TO 15
      IF(IL.NE.-1) GO TO 16
      W=W+1.0D0/3.0D0
      IF(J.EQ.2.AND.I.NE.1) W=W+1.0D0/3.0D0
      IF(J.EQ.3.AND.I.EQ.1) W=W+1.0D0/3.0D0
      GO TO 15
   16 IF(IL.EQ.3) W=W+0.5D0
      IF(IL.EQ.2.AND.I.NE.J-1) W=W+0.5D0
      IF(IL.EQ.4.AND.I.NE.3) W=W+0.5D0
   15 XS(I)=W
   14 CONTINUE
      DO 21 IX=1,NX
      X=XS(1)+DBLE(IX-MX)
      XD=X-VATOM(1,I1)
      IF(ABS(XD).GT.XM) GO TO 21
      DO 22 IY=1,NY
      Y=XS(2)+DBLE(IY-MY)
      YD=Y-VATOM(2,I1)
      IF(ABS(YD).GT.YM) GO TO 22
      DO 23 IZ=1,NZ
      Z=XS(3)+DBLE(IZ-MZ)
      ZD=Z-VATOM(3,I1)
      IF(ABS(ZD).GT.ZM) GO TO 23
       X8=IS
      R=SQRT(XD*XD+(YD*P)**2+(ZD*Q)**2
     & +2.0D0*P*Q*YD*ZD*CA+2.0D0*Q*ZD*XD*CB+2.0D0*P*XD*YD*CC)
      IF(R.LT.1.0D-5) GO TO 23
      IF(R.GT.RM) GO TO 23
      IF(NN.EQ.0) GO TO 34
      DO 31 IP=1,NN
      M=NN-IP+1
      IF((XA(7,M)-R).LT.-1.0D-5) GO TO 32
      IF((XA(7,M)-R).LT.1.0D-5
     &  .AND.XA(8,M).LE.X8) GO TO 32
      MP=M+1
      IF(MP.GT.201) GO TO 31
      DO 33 I=1,8
   33 XA(I,MP)=XA(I,M)
   31 CONTINUE
   34 M=0
   32 MP=M+1
      IF(MP.GT.201) GO TO 23
      XA(1,MP)=X
      XA(2,MP)=Y
      XA(3,MP)=Z
      XA(4,MP)=XD
      XA(5,MP)=YD
      XA(6,MP)=ZD
      XA(7,MP)=R
      XA(8,MP)=X8
      IF(NN.LT.201) NN=NN+1
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
   13 CONTINUE
   12 CONTINUE
      IF(ITSBO.EQ.0) GO TO 43
      RX=0.0D0
      N8=0
      IF(NN.EQ.0) GO TO 61
      DO 44 J=1,NN
      N=XA(8,J)+0.5D0
      IF(ABS(RX-XA(7,J)).GT.1.0D-5.OR.N.NE.N8) GO TO 50
      IF(INDA.EQ.0) GO TO 44
      GO TO 45
   50 RX=XA(7,J)
      N8=N
      IA1=IA(I1)
      IA2=IA(N8)
      IF(IA1.LE.IA2) GO TO 47
      IW=IA1
      IA1=IA2
      IA2=IW
   47 DO 46 K=1,NBO
      IF(KBO(1,K).NE.IA1) GO TO 46
      IF(KBO(2,K).NE.IA2) GO TO 46
      IF(ABS(BOL(K)-RX).LT.1.0D-5) GO TO 48
   46 CONTINUE
      INDA=0
      GO TO 44
   48 IN=(K-1)*NAB+(I1-1)*NATOM+N8
      IVEC(IN)=NVEC+1
C     WRITE(6,999) IN,NVEC,IVEC(IN),K,I1,N8
      INDA=1
   45 NVEC=NVEC+1
      IF(NVEC.LE.4000) GO TO 92
      WRITE(6,902)
  902 FORMAT(' STOP AT 92 IN TSBOVR')
   92 CONTINUE
      IVEC(IN)=IVEC(IN)-(IVEC(IN)/10000)*10000+NVEC*10000
C     WRITE(6,999) IN,NVEC,IVEC(IN),K,I1,N8
      DO 49 I=1,3
   49 VST(I,NVEC)=XA(I+3,J)
   44 CONTINUE
C     WRITE(6,999) (IVEC(IN),IN=1,16)
C 999 FORMAT(8I10)
   43 IF(JF.EQ.0) GO TO 11
      WRITE(JF,610) RM
  610 FORMAT(/' TABLE OF NEIGHBOURS IN THE RANGE OF',F9.5
     &   ,' FROM THE ATOM AT')
      WRITE(JF,611) (VATOM(I,I1),I=1,3),IA(I1),I1
  611 FORMAT(3F9.5,' WITH KIND NUMBER OF',I2,' AND ATOM NUMBER OF',I3)
      WRITE(JF,612) A,B,C,CA,CB,CC
  612 FORMAT(' LATTICE CONSTANTS  A=',F9.5,'  B=',F9.5,'  C=',F9.5/
     &       '                   CA=',F9.5,' CB=',F9.5,' CC=',F9.5)
      IF(NN.EQ.0) GO TO 61
      WRITE(JF,602)
  602 FORMAT(3X,1HX,7X,1HY,7X,1HZ,6X,2HXD,6X,2HYD,6X,2HZD
     &   ,4X,' DISTANCE/A')
      RX=0.0D0
      M8=0
      K=0
      DO 62 J=1,NN
      N=XA(8,J)+0.5D0
      M=IA(N)
      IF(ABS(RX-XA(7,J)).GT.1.0D-5) I=0
      IF(M.GT.M8) I=0
      IF(I.EQ.0) RX=XA(7,J)
      IF(I.EQ.0) M8=M
      IF(I.EQ.0) K=K+1
      I=I+1
      WRITE(JF,601) (XA(L,J),L=1,7),M,N,J,K,I
  601 FORMAT(6F8.4,F9.5,2I3,I4,2I3)
   62 CONTINUE
      IF(NN.EQ.201) WRITE(JF,*) ' THERE MAY EXIST MORE NEIGHBOURS'
      GO TO 11
   61 IF(JF.NE.0) WRITE(JF,603)
  603 FORMAT(' WE CAN NOT FIND ANY NEIGHBOUR')
   11 CONTINUE
      IF(ITSBO.EQ.0) RETURN
C     DO 51 J=1,NATOM
C     DO 51 I=1,NATOM
C     DO 51 K=1,NBO
C     IN=(K-1)*NAB+(J-1)*NATOM+I
C     IF(IVEC(1,IN).EQ.0) GO TO 51
C      WRITE(6,604) I,J,K,IN,(IVEC(N,IN),N=1,2)
C 604 FORMAT(3I3,3I5)
C      WRITE(6,605) ((VST(N,M),N=1,3),M=IVEC(1,IN),IVEC(2,IN))
C 605 FORMAT(9F8.4)
C  51 CONTINUE
      RETURN
      END
C SUBROUTINE FACDTH ====*====3====*====4====*====5====*====6====*====7
C
C   SLAVE ROUTINE OF TSNMKP
C   FOR FACE CENTERED CUBIC LATTICE WITH TH
C   ICC :COMMON DENOMINATOR
C   CN  :NAME
C                 1988.11.17 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FACDTH(KKX,KKY,KKZ,ICC,CN)
      CHARACTER*2 CN
      KZ=IABS(KKZ)
      KY=IABS(KKY)
      KX=IABS(KKX)
      DO 1 N=1,2
          IF(KZ.GT.KY.OR.KZ.GT.KX) THEN
              IW=KZ
              KZ=KY
              KY=KX
              KX=IW
          ELSE
              GO TO 2
          END IF
    1 CONTINUE
    2 CONTINUE
      IF(KZ.EQ.0) THEN
          IF(KY.EQ.0) THEN
              IF(KX.EQ.0) THEN
                  CN='GM'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='X '
              ELSE
                  CN='DT'
              END IF
          ELSE IF(KY.EQ.ICC) THEN
              IF(KX.EQ.0) THEN
                  CN='X '
              ELSE
                  CN='Z1'
              END IF
          ELSE IF(KX.EQ.ICC) THEN
              IF(KY*2.EQ.ICC) THEN
                  CN='W '
              ELSE
                  CN='Z2'
              END IF
          ELSE IF(KX.EQ.KY) THEN
                  CN='SM'
          ELSE
              IF(KX.EQ.0) THEN
                  CN='DT'
              ELSE IF(KX.EQ.ICC) THEN
                  CN='Z2'
              ELSE
                  CN='ZP'
              END IF
          END IF
      ELSE IF(KZ.EQ.KY) THEN
          IF(KY.EQ.KX) THEN
              IF(KX*2.EQ.ICC) THEN
                  CN='L '
              ELSE
                  CN='LD'
              END IF
          ELSE IF(KX.EQ.ICC) THEN
                  CN='S '
          ELSE
                  CN='XY'
          END IF
      ELSE IF(KZ.EQ.KX) THEN
          IF(KY.EQ.ICC) THEN
                  CN='S '
          ELSE
                  CN='XY'
          END IF
      ELSE
          IF(KY.EQ.ICC) THEN
                  CN='B1'
          ELSE IF(KY.EQ.KX) THEN
                  CN='XY'
          ELSE IF(KY*2.EQ.ICC.AND.KX+KZ.EQ.ICC) THEN
                  CN='Q '
          ELSE IF(KX*2.EQ.ICC.AND.KY+KZ.EQ.ICC) THEN
                  CN='Q '
          ELSE
              IF(KX.EQ.ICC) THEN
                  CN='B2'
              ELSE 
                  CN='GN'
              END IF
          END IF
      END IF
      RETURN
      END
C SUBROUTINE EQUIKK ====*====3====*====4====*====5====*====6====*====7
C
C    IF KB(3)/IC IS EQUIVALENT TO KBB(3)/ICC THEN
C       KG(K)=KBB(K)/ICC-KS(K,IS)/ICBB
C       IND=JS(IS) :OPERATION CODE TRANSLATING KS(K,IS) FROM KB/IC 
C       WHERE KS(K,IS)/ICBB IS ONE OF THE STAR OF KB(K)/IC
C    ELSE IF KB(3)/IC IS EQUIVALENT TO -KBB(3)/ICC THEN
C       KG(K)=-KBB(K)/ICC-KS(K,IS)/ICBB
C       IND=-JS(IS) :OPERATION CODE TRANSLATING KS(K,IS) FROM KB/IC 
C       WHERE KS(K,IS)/ICBB IS ONE OF THE STAR OF KB(K)/IC
C    ELSE 
C       KG(K)=999, IND=0
C                1994/2/15
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE EQUIKK(KB,IC,KBB,ICC,KG,IND)
      COMPLEX*16 CW
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      DIMENSION KBB(3),KB(3),KG(3)
C
      CALL CHKDNM(IC)
      CALL CHKDNM(ICC)
      CALL TSIREP(KB,IC,0)
      CALL ZZZY38
      DO 1 IS=1,NS
      DO 2 K=1,3
      KG(K)=KBB(K)*ICBB-KS(K,IS)*ICC
      IF(MOD(KG(K),ICBB*ICC).NE.0) GO TO 1
      KG(K)=KG(K)/(ICBB*ICC)
    2 CONTINUE
      IF(IL.EQ.-1) THEN
          IF(MOD(IABS(-KG(1)+KG(2)+KG(3)),3).NE.0) GO TO 1
      ELSE IF(IL.EQ.2) THEN
          IF(MOD(IABS(KG(1)),2).NE.MOD(IABS(KG(2)),2)) GO TO 1
          IF(MOD(IABS(KG(2)),2).NE.MOD(IABS(KG(3)),2)) GO TO 1
      ELSE IF(IL.EQ.3) THEN
          IF(MOD(IABS(KG(1)+KG(2)+KG(3)),2).NE.0) GO TO 1
      ELSE IF(IL.EQ.4) THEN 
          IF(MOD(IABS(KG(1)+KG(2)),2).NE.0) GO TO 1
      END IF
      IND=JS(IS)
      RETURN
    1 CONTINUE
      DO 3 IS=1,NS
      DO 4 K=1,3
      KG(K)=-KBB(K)*ICBB-KS(K,IS)*ICC
      IF(MOD(IABS(KG(K)),ICBB*ICC).NE.0) GO TO 3
      KG(K)=KG(K)/(ICBB*ICC)
    4 CONTINUE
      IF(IL.EQ.-1) THEN
          IF(MOD(IABS(-KG(1)+KG(2)+KG(3)),3).NE.0) GO TO 3
      ELSE IF(IL.EQ.2) THEN
          IF(MOD(IABS(KG(1)),2).NE.MOD(IABS(KG(2)),2)) GO TO 3
          IF(MOD(IABS(KG(2)),2).NE.MOD(IABS(KG(3)),2)) GO TO 3
      ELSE IF(IL.EQ.3) THEN
          IF(MOD(IABS(KG(1)+KG(2)+KG(3)),2).NE.0) GO TO 3
      ELSE IF(IL.EQ.4) THEN 
          IF(MOD(IABS(KG(1)+KG(2)),2).NE.0) GO TO 3
      END IF
      IND=-JS(IS)
      RETURN
    3 CONTINUE
      KG(1)=999
      KG(2)=999
      KG(3)=999
      IND=0
      RETURN
      END   
C SUBROUTINE KPNAME ====*====3====*====4====*====5====*====6====*====7
C
C   NAME OF KPOINT 
C   THIS SUBROUTINE GIVES THE NAME OF K-POINT : NP
C   EVEN WHEN THE POINT IS OUTSIDE OF FIRST B.Z.
C   THE RECIPROCAL LATTICE FOR THE POINT IS GIVEN IN IRECP(3)
C
C                 1994.02.20 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE KPNAME(KX,KY,KZ,IC,NP)
      DIMENSION KB(3),KBB(3,10),IRECP(3,10)
      CHARACTER*2 NP,NAME
      CALL CHKDNM(IC)
      KB(1)=KX
      KB(2)=KY
      KB(3)=KZ
      ICC=IC
      CALL NEAREC(KB,ICC,KBB,IRECP,NG)
      KKX=KBB(1,1)
      KKY=KBB(2,1)
      KKZ=KBB(3,1)
      CALL TSNMKP(KKX,KKY,KKZ,ICC,NAME)
C      write(6,*) kkx,kky,kkz,icc,name
      NP=NAME
      RETURN
      END
C SUBROUTINE TSKPGN ====*====3====*====4====*====5====*====6====*====7
C
C   KPOINTS IN THE FIRST B.Z. ARE GENERATED IN KKM.
C   ICC IS COMMON DENOMINATOR AND NK IS NUMBER OF KPOINT.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE X-DIRECTION
C        IS DIVIDED BY NX.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE Y-DIRECTION
C        IS DIVIDED BY NY.
C   AXIS OF GAMMA TO ZONE BOUNDARY CENTER IN THE Z-DIRECTION
C        IS DIVIDED BY NZ.
C   ZONE BOUNDARY CENTER MEANS, HERE, MID-POINT BETWEEN GAMMA AND
C   THE NEAREST REXIPROCAL LATTICE POINT IN THE RESPECTIVE DIRECTION.
C
C                 1988.10.18 :  A. YANASE
C   MODIFIED TO INCLUDE NX=0 OR/AND NY=0 OR/AND NZ=0
C                 1994.03.17 :  A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE TSKPGN(NX,NY,NZ,KKM,ICC,NK)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
C
      DIMENSION KKM(3,1)
      IF(NX.LT.0.OR.NY.LT.0.OR.NZ.LT.0) THEN
         WRITE(6,900) NX,NY,NZ
  900    FORMAT(' ALL OF NX,NY,NZ SHOULD NOT BE NEGATIVE',
     &       ' HOWEVER NX=',I3,' NY=',I3,' NZ=',I3,' THEN STOP HERE')
         STOP
      END IF
      CALL ADDINV(NTLATC)
      CALL ZZZY53(NX,NY,KK)
      IF(NX.NE.0.AND.NY.NE.0) THEN
         LCM=NX*(NY/KK)
      ELSE IF(NX.NE.0) THEN
         LCM=NX
      ELSE
         LCM=NY
      END IF
      CALL ZZZY53(LCM,NZ,KK)
      IF(LCM.NE.0.AND.NZ.NE.0) THEN
         LCM=LCM*(NZ/KK)
      ELSE IF(LCM.EQ.0) THEN
         LCM=NZ
      END IF
C      WRITE(6,*) NX,NY,NZ,LCM
      KDX=0
      KDY=0
      KDZ=0
      IF(IL.EQ.-1) THEN
         IF(MOD(LCM,3).NE.0) THEN
           IF(NX.NE.0) KDX=3*(LCM/NX)
           IF(NY.NE.0) KDY=3*(LCM/NY)
           IF(NZ.NE.0) KDZ=3*(LCM/NZ)
           ICC=2*LCM
           IF(ICC.EQ.0) ICC=1
         ELSE
           IF(NX.NE.0) KDX=LCM/NX
           IF(NY.NE.0) KDY=LCM/NY
           IF(NZ.NE.0) KDZ=LCM/NZ
           ICC=2*(LCM/3)
           IF(ICC.EQ.0) ICC=1
         END IF
      ELSE IF(IL.EQ.0.OR.IL.EQ.1) THEN
           IF(NX.NE.0) KDX=LCM/NX
           IF(NY.NE.0) KDY=LCM/NY
           IF(NZ.NE.0) KDZ=LCM/NZ
           ICC=LCM*2
           IF(ICC.EQ.0) ICC=1
      ELSE IF(IL.EQ.2.OR.IL.EQ.3) THEN
           IF(NX.NE.0) KDX=LCM/NX
           IF(NY.NE.0) KDY=LCM/NY
           IF(NZ.NE.0) KDZ=LCM/NZ
           ICC=LCM
           IF(ICC.EQ.0) ICC=1
      ELSE IF(IL.EQ.4) THEN
           IF(NX.NE.0) KDX=2*LCM/NX
           IF(NY.NE.0) KDY=2*LCM/NY
           IF(NZ.NE.0) KDZ=LCM/NZ
           ICC=LCM*2
           IF(MOD(KDZ,2).EQ.0) THEN
               KDZ=KDZ/2
               KDY=KDY/2
               KDX=KDX/2
               ICC=ICC/2
           END IF
           IF(ICC.EQ.0) ICC=1
      END IF
C     WRITE(6,901)
C 901 FORMAT(/' GENERATED KPOINT'
C    &       /'    NO  KX, KY, KZ')
      NK=0
      NNY=NY
      NNX=NX
      NNZ=NZ
      IF(IL.LE.0) THEN
         NNX=NX+NY/3
         NNY=NY+NX/3
      ELSE IF((NTLATC.LE.3.AND.NTLATC.GE.1).OR.NTLATC.EQ.9) THEN
         NNX=2*NX
         NNY=2*NY
         NNZ=2*NZ
      END IF
      DO 1 IZ=0,NNZ
         KKZ=IZ*KDZ
         DO 4 JX=1,2
            DO 5 JY=1,2
               DO 2 IY=0,NNY
                  KY=IY*KDY
                  KKY=KY*(3-2*JY)
                  DO 3 IX=0,NNX
                     KX=IX*KDX
                     KKX=KX*(3-2*JX)
                     CALL TSKFBZ(KKX,KKY,KKZ,ICC,IND)
                     IF(IND.EQ.0) GO TO 3
                     CALL KALRST(KKX,KKY,KKZ,ICC,KKM,NK,IND)
                     IF(IND.EQ.1) THEN
                        NK=NK+1
                        KKM(1,NK)=KKX
                        KKM(2,NK)=KKY
                        KKM(3,NK)=KKZ
C                        WRITE(6,902) NK,KKX,KKY,KKZ,ICC,IND
C  902                   FORMAT(1H ,I5,3I4,'/',I4,' IND=',I2)
                     END IF
    3             CONTINUE
    2          CONTINUE
    5       CONTINUE
    4    CONTINUE
    1 CONTINUE
      CALL REMINV
      RETURN
      END
C SUBROUTINE ADDINV ====*====3====*====4====*====5====*====6====*====7
C
C   THE INVERSION SYMMETRY IS INCLUDED TENTATIVILY TO IMITATE
C   THE TIME REVERSAL SYMMETRY
C   NTLATC=8     OH  SIMMETRY
C   NTLATC=7     TH  SIMMETRY
C   NTLATC=6     D4H SIMMETRY
C   NTLATC=5     C4H SIMMETRY
C   NTLATC=4     D2H SIMMETRY
C   NTLATC=9     C2H SIMMETRY 2-FOLD AXIS//A-AXIS
C   NTLATC=3     C2H SIMMETRY 2-FOLD AXIS//B-AXIS
C   NTLATC=2     C2H SIMMETRY 2-FOLD AXIS//C-AXIS
C   NTLATC=1     CI  SIMMETRY
C   NTLATC=0     D6H SIMMETRY
C   NTLATC=-1    D3D SIMMETRY
C   NTLATC=-2    C6H SIMMETRY
C   LTLATC=-3    C3I SIMMETRY
C                 1988.11.4 :  A. YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE ADDINV(NTLATC)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/,IIG,JJV,ILL,NNG
C
      LOGICAL LA,LB,LC,LCUB,LTET,LSIX
C
      DIMENSION JB(2,3),JJV(2,3,48),IIG(48)
      DO 10 I=1,NG
      IIG(I)=IG(I)
      DO 11 K=1,3
      DO 12 J=1,2
      JJV(J,K,I)=JV(J,K,I)
      IF(J.EQ.1) JV(J,K,I)=0
      IF(J.EQ.2) JV(J,K,I)=1
   12 CONTINUE
   11 CONTINUE
   10 CONTINUE
      NNG=NG
      ILL=IL
      DO 20 K=1,3
      JB(1,K)=0
      JB(2,K)=1
   20 CONTINUE
      IF(IL.GE.1) NA=25
      IF(IL.LE.0) NA=13
      CALL TSGENR(NA,JB)
      CALL TSPGRP(0)
      IF(IL.GE.1) THEN
         LA=.FALSE.
         LB=.FALSE.
         LC=.FALSE.
         LCUB=.FALSE.
         LTET=.FALSE.
         DO 30 I=1,NG
         IF(IG(I).EQ.2) LA=.TRUE.
         IF(IG(I).EQ.3) LB=.TRUE.
         IF(IG(I).EQ.4) LC=.TRUE.
         IF(IG(I).EQ.5) LCUB=.TRUE.
         IF(IG(I).EQ.21) LTET=.TRUE.
   30    CONTINUE
         IF(LCUB.AND.LTET)        NTLATC=8
         IF(LCUB.AND.(.NOT.LTET)) NTLATC=7
         IF(LA.AND.LB.AND.LC.AND.(.NOT.LCUB).AND.LTET) NTLATC=6
         IF((.NOT.LA).AND.(.NOT.LB).AND.LTET)          NTLATC=5
         IF(LA.AND.LB.AND.LC.AND.(.NOT.LCUB).AND.(.NOT.LTET)) NTLATC=4
         IF(LA.AND.(.NOT.LB).AND.(.NOT.LC)) NTLATC=9
         IF((.NOT.LA).AND.LB.AND.(.NOT.LC)) NTLATC=3
         IF((.NOT.LA).AND.(.NOT.LB).AND.LC) NTLATC=2
         IF((.NOT.LA).AND.(.NOT.LB).AND.(.NOT.LC)) NTLATC=1
      ELSE IF(IL.LE.0) THEN
         LA=.FALSE.
         LB=.FALSE.
         LSIX=.FALSE.
         DO 31 I=1,NG
         IF(IG(I).EQ.2) LSIX=.TRUE.
         IF(IG(I).EQ.7) LA=.TRUE.
         IF(IG(I).EQ.10) LB=.TRUE.
   31    CONTINUE
         IF(LSIX.AND.LA)        NTLATC=0
         IF((.NOT.LSIX).AND.(LA.OR.LB)) NTLATC=-1
         IF(LSIX.AND.(.NOT.LA).AND.(.NOT.LB)) NTLATC=-2
         IF((.NOT.LSIX).AND.(.NOT.LA).AND.(.NOT.LB)) NTLATC=-3
      END IF
      RETURN
C ENTRY REMINV ====2====*====3====*====4====*====5====*====6====*====7
C
C         RECOVER TO ORIGINAL SPACE GROUP
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY REMINV
      DO 15 I=1,NNG
      IG(I)=IIG(I)
      DO 16 K=1,3
      DO 17 J=1,2
      JV(J,K,I)=JJV(J,K,I)
   17 CONTINUE
   16 CONTINUE
   15 CONTINUE
      NG=NNG
      IL=ILL
      RETURN
      END
C SUBROUTINE KALRST ===*====3====*====4====*====5====*====6====*====7
C
C    IF STAR OF (KX,KY,KZ)/IC IS ALREADY IN KKM, IND=0
C    (KX,KY,KZ)/IC IS NEW                      , IND=1
C               1988.10.18  AKIRA YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE KALRST(KX,KY,KZ,IC,KKM,NK,IND)
        IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 CW
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &  ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG2/,/SPG3/,/STK/
      DIMENSION KKM(3,1),KKB(3),KA(3),KC(3)
      KKB(1)=KB(1)
      KKB(2)=KB(2)
      KKB(3)=KB(3)
      ICBB=ICB
      KB(1)=KX
      KB(2)=KY
      KB(3)=KZ
      ICB=IC
C---------------------------------------------
C   CALCULATION OF STAR
      CALL ZZZY38
C---------------------------------------------
      IF(NK.EQ.0) GO TO 11
      DO 1 IS=1,NS
      DO 2 J=1,NK
         INDC=0
      DO 3 K=1,3
         KA(K)=KS(K,IS)-KKM(K,J)
         IF(MOD(IABS(KA(K)),ICB).NE.0) GO TO 2
         KA(K)=KA(K)/ICB
         KC(K)=IABS(KA(K))
         IF(KA(K).NE.0) INDC=1
    3 CONTINUE
      IF(INDC.EQ.1) THEN
         IF(IL.EQ.-1) GO TO 4
         IF(IL.EQ.0) GO TO 5
         IF(IL.EQ.1) GO TO 5
         IF(IL.EQ.2) GO TO 6
         IF(IL.EQ.3) GO TO 7
         IF(IL.EQ.4) GO TO 8
    4    IIWW=-KA(1)+KA(2)+KA(3)
         IF(MOD(IABS(IIWW),3).NE.0) GO TO 2
         GO TO 5
    6    IF(MOD(KC(1),2).NE.MOD(KC(2),2)) GO TO 2
         IF(MOD(KC(2),2).NE.MOD(KC(3),2)) GO TO 2
         GO TO 5
    7    IF(MOD(KC(1)+KC(2)+KC(3),2).NE.0) GO TO 2
         GO TO 5
    8    IF(MOD(KC(1)+KC(2),2).NE.0) GO TO 2
         GO TO 5
       ELSE
         GO TO 5
       END IF
    2 CONTINUE
    1 CONTINUE
   11    IND=1
         GO TO 13
    5    IND=0
   13 KB(1)=KKB(1)
      KB(2)=KKB(2)
      KB(3)=KKB(3)
      ICB=ICBB
      RETURN
      END
C SUBROUTINE CORRES ====*====3====*====4====*====5====*====6====*====7
C
C    WHEN KB(3)/IC IS EQUIVALENT TO KBB(3)/ICC, CORRESPONDENCE TABLE
C    ITCR(12) IS GIVEN BY THE CALCULATION FOR THE CHARACTER OF 
C    FULL REPRESENTATION OF SPACE GROUP
C    IW=0 SINGLE REPRESENTATION
C    IW=1 DOUBLE REPRESENTATION
C    NRR: NUMBER OF REPRESENTATION
C    IND=1,  ITCR(IR) FOR KB IS THE SAME AS IR FOR KBB
C    IND=-1, ITCR(IR) FOR KB IS THE TIME REVERSAL OF IR FOR KBB
C    IND=0 NOT
C                 1994/2/28
C                    BY S. TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CORRES(KB,IC,KBB,ICC,IW,ITCR,NRR,IND)
       IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 CCW1(48,12),CCW2(48,12)
      COMPLEX*16 CCW(48,12)     
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/
      DIMENSION LAT(2,3),ND(12),LATG(3),JTR(12),IPA(12)
      DIMENSION KG(3),JTCR(12)
      DIMENSION KB(3),KBB(3),itcr(12)
C---------------------------------
C   input data check
C---------------------------------
      CALL CHKDNM(IC)
      CALL CHKDNM(ICC)
      IF(IW.NE.0.AND.IW.NE.1) THEN
        WRITE(6,*) ' IW SHOLD BE 0 OR 1, HOWEVER IW=',IW
        WRITE(6,*) ' STOP IN CORRES'
        STOP
      END IF
      CALL EQUIKK(KB,IC,KBB,ICC,KG,INDE)
      IF(INDE.EQ.0) THEN
        IND=0
        RETURN
      END IF
      CALL TSIREP(KB,IC,IW)
      CALL DGTRST(JDOB,NRR,NH,NSTR,ND,JTR,IPA)
C      write(6,*) kb,ic,nrr,(jtr(it),it=1,nrr),(ipa(it),it=1,nrr)
C---------------------------
C  CLEAR ITCR AND JTCR
C---------------------------
      DO 10 IT=1,12
      ITCR(IT)=0
      JTCR(IT)=0
   10 CONTINUE
C-------------------------
C  COMPAIR KB WITH KBB
C-------------------------
      DO 11 IP=1,3
      IF((KB(IP)*ICC).NE.(KBB(IP)*IC)) GO TO 12
   11 CONTINUE
C--------------------------------
C   KB/IC is the same as KBB/ICC
C--------------------------------
      DO 15 IT=1,NRR
      ITCR(IT)=IT
   15 CONTINUE
      IND=1
      RETURN
C----------------------------
   12 CONTINUE
      ICD=0
      DO 100 NREPT=1,32
      CALL GENLTP(LATG,ICD)
      LAT(1,1)=LATG(1)
      LAT(1,2)=LATG(2)
      LAT(1,3)=LATG(3)
      LAT(2,1)=ICD
      LAT(2,2)=ICD
      LAT(2,3)=ICD
C-------------------------
C  CAL. FOR CCW1
C-------------------------
      CALL TSIREP(KB,IC,IW)
      CALL TSCRFR(LAT,CCW)
      DO 2 J=1,NRR
      DO 3 K=1,NG
      CCW1(K,J)=CCW(K,J)
    3 CONTINUE
    2 CONTINUE
C-------------------------
C  CAL. FOR CCW2
C-------------------------
      CALL TSIREP(KBB,ICC,IW)
      CALL TSIRNR(NR,NH,ND)
      IF(NR.NE.NRR) GO TO 911
      CALL TSCRFR(LAT,CCW)
      DO 5 J=1,NR
      DO 6 K=1,NG
      CCW2(K,J)=CCW(K,J)
    6 CONTINUE
    5 CONTINUE
C---------------------------
C  COMPAIR CCW1 WITH CCW2
C---------------------------       
      DO 7 IR2=1,NR
      IF(ITCR(IR2).NE.0) GO TO 7
      DO 8 IR1=1,NRR
      IF(JTCR(IR1).NE.0) GO TO 8
      JHATA=0
      DO 9 IGG=1,NG
      if(abs(ccw2(igg,ir2)).gt.1.0d-4) then
        jhata=1
        IF(INDE.GT.0) THEN
          IF(ABS(CCW1(IGG,IR1)-CCW2(IGG,IR2)).GT.1.0D-5) GO TO 8
        else
          IF(ABS(CCW1(IGG,IR1)-DCONJG(CCW2(IGG,IR2))).GT.1.0D-5) 
     &           GO TO 8
        end if
      else
        if(abs(ccw1(igg,ir1)).gt.1.0d-4) go to 8 
      end if
    9 CONTINUE
      if(jhata.eq.0) go to 7
      IF(IR2.EQ.1) GO TO 16
      DO 13 IRR=1,IR2-1
         IF(ITCR(IRR).EQ.IR1) then
             itcr(ir2)=0
             itcr(irr)=0
             go to 7
         end if
   13 CONTINUE  
   16 continue
        ITCR(IR2)=IR1
        GO TO 7
    8 CONTINUE
      if(inde.gt.0) then
         WRITE(6,*) ' WE CANNOT FIND FOR',IR2
      else
         WRITE(6,*) ' WE CANNOT FIND TIME REVERSAL FOR',IR2
      end if
      write(6,*) ' STOP AT 8 IN CORRES'
      STOP
    7 CONTINUE
C      write(6,*) ' kb=',kb,ic,' kbb=',kbb,ic
C      write(6,*) (itcr(iiii),iiii=1,nr)
      ihata=1
      do 71 ir2=1,nr
        if(itcr(ir2).ne.0) then
           jtcr(itcr(ir2))=ir2
        else
           ihata=0
        end if 
   71 continue
      if(ihata.eq.0) go to 100
      if(inde.gt.0) IND=1
      if(inde.lt.0) IND=-1
      RETURN
  100 CONTINUE
      write(6,*) ' kb=',kb,ic,' kbb=',kbb,ic
      WRITE(6,*) ' WE CANNOT DETERMINE THE CORRESPONDENCE'
      WRITE(6,*) ' STOP AT 100 IN CORRES'
      STOP
  911 WRITE(6,*) ' NUMBERS OF REPRESENTATION ARE DIFFERENT'
      WRITE(6,*) ' FROM EACH OTHER NRR=',NRR,' NR=',NR 
      WRITE(6,*) ' STOP AT 911 IN CORRES'
      STOP
      END    
C SUBROUTINE COMPAT ====*====3====*====4====*====5====*====6====*====7
C
C    COMPATIBILITY OF THE REPRESENTATION OF KB(3)/IC
C    TO THAT OF KBB(3)/ICC, WHEN KP-GROUP OF KB IS A SUBGROUP
C    OF KP-GROUP OF KBB
C    NR:  NUMBER OF REPRESENTATION OF KB/IC
C    NRT: NUMBER OF REPRESENTATION OF KBB/ICC
C    IWW: 0 SINGLE REP., 1 DOUBLE REP.
C    ICP(1-NRT,1-NR): TABLE OF COMPATIBILITY
C
C                1994/02/28
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE COMPAT(KB,IC,NR,KBB,ICC,NRT,IWW,ICP,INDC)
      DIMENSION KB(3),KBB(3),ICP(12,12)
      DIMENSION IFA(12),ND(12),NDT(12),JG(48),JGT(48),JGA(2,3,48)
      DIMENSION NTAB(12)
      COMPLEX*16 WD(48),CW(48,12),CWT(48),CR(48)
C
      CALL CHKDNM(IC)
      CALL CHKDNM(ICC)
      IF(IWW.NE.0.AND.IWW.NE.1) THEN
        WRITE(6,*) ' IW SHOLD BE 0 OR 1, HOWEVER IWW=',IWW
        WRITE(6,*) ' STOP IN COMPAT'
        STOP
      END IF
      CALL TSIREP(KB,IC,IWW)
      CALL TSIRNR(NR,NK,ND)
      DO 41 IR=1,NR
      CALL IRMATA(IR,1,1,MG,NND,JG,JGA,CR,WD)
      DO 42 I=1,MG
      CW(I,IR)=CR(I)
   42 CONTINUE
   41 CONTINUE
      CALL TSIREP(KBB,ICC,IWW)
      CALL TSIRNR(NRT,NKT,NDT)
      CALL IRMATA(1,1,1,MGT,NND,JGT,JGA,CR,WD)
      CALL SUBGRP(MG,JG,MGT,JGT,NTAB,IND)
      IF(IND.NE.0) GO TO 99
C-----------------------------------------------
C   MAIN PROCESS STARTS
C-----------------------------------------------
      DO 100 IRT=1,NRT
      IRRT=IRT
      CALL IRMATA(IRRT,1,1,MGT,NND,JGT,JGA,CWT,WD)
C-----------------------------------------------
C   DIMESNSION TEST
C-----------------------------------------------
      DO 1 I=1,NR
    1 IFA(I)=0
      NCOUNT=0
   12 IRP=1
   14 IFA(IRP)=IFA(IRP)+1
      NCOUNT=NCOUNT+1
      IF(NCOUNT.GT.500) GO TO 98
   11 CALL UPTEST(IFA,ND,NND,NR,IND)
c      WRITE(6,*) (IFA(I),I=1,NR),IND,IRP,IRT
      IF(IND.EQ.0) GO TO 20
      IF(IND.LT.0) GO TO 12
   13 IF(IRP.GE.NR) GO TO 98
      DO 2 I=1,IRP
      IFA(I)=0
    2 CONTINUE
      IRP=IRP+1
      GO TO 14
C--------------------------------------------
C    COMPATIBILITY CHECK
C-------------------------------------------- 
   20 CALL TESTCR(IFA,NR,CW,CWT,NTAB,MG,IND)
      IF(IND.EQ.0) GO TO 21
      GO TO 13
   21 CONTINUE
      DO 40 IR=1,NR
         ICP(IRRT,IR)=IFA(IR)
   40 CONTINUE
  100 CONTINUE
      INDC=0
      RETURN
   98 WRITE(6,*) ' COMPATIBILITY TEST IS FAILED THEN STOP HERE'
      STOP
   99 WRITE(6,*) ' K-POINT GROUP OF KB IS NOT SUBGROUP 
     &  OF K-POINT GROUP OF KBB'
      INDC=1
      RETURN
      END
C SUBROUTINE UPTEST ====*====3====*====4====*====5====*====6====*====7
C
C    IF IF NDSUM(IFA(I)*ND(I),I=1,NR)=NDT THEN
C          IND=0
C    ELSE IF NDSUM < NDT THEN
C          IND=-1
C    ELSE IF NDSUM > NDT THEN
C          IND=1
C                 1993/12/25
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE UPTEST(IFA,ND,NDT,NR,IND)
      DIMENSION ND(12),IFA(12)
      NDSUM=0
      DO 10 I=1,NR
      NDSUM=NDSUM+IFA(I)*ND(I)
   10 CONTINUE
      IF(NDSUM.LT.NDT) GO TO 1 
      IF(NDSUM.GT.NDT) GO TO 2 
      IND=0
      RETURN
    1 IND=-1
      RETURN
    2 IND=1
      RETURN
      END
C SUBROUTINE SUBGRP ====*====3====*====4====*====5====*====6====*====7
C
C    IF (JG(I),I=1,MG) IS A SUBGROUP OF (JGT(J),J=1,MGT) THEN 
C          TABLE (NTAB(I),I=1,MG) IS MADE HERE AND IND=0
C    ELSE 
C          IND=-1
C
C                 1993/12/25
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE SUBGRP(MG,JG,MGT,JGT,NTAB,IND)
      DIMENSION NTAB(48),JG(48),JGT(48)
C
      DO 1 I=1,MG
      DO 2 J=1,MGT
      IF(JG(I).EQ.JGT(J)) GO TO 3 
    2 CONTINUE
      IND=-1
      RETURN
    3 NTAB(I)=J
    1 CONTINUE
      IND=0
      RETURN
      END
C SUBROUTINE TESTCR ====*====3====*====4====*====5====*====6====*====7
C
C    IF THE CHARACTERS SATISFY THE COMPATIVILITY RELATION
C          IND=0
C    ELSE 
C          IND=1
C
C                 1993/12/25
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TESTCR(IFA,NR,CW,CWT,NTAB,MG,IND)
      COMPLEX*16 CW(48,12),CWT(48),CWSUM
      DIMENSION IFA(12),NTAB(12)      
C
c      WRITE(6,*) ' We are now in TEST'
      DO 20 J=1,MG
      CWSUM=0.0
      DO 10 I=1,NR
      CWSUM=CWSUM+IFA(I)*CW(J,I)
   10 CONTINUE
      IF(ABS(CWSUM-CWT(NTAB(J))).GT.0.0001D0) GO TO 1
   20 CONTINUE
      IND=0
c      WRITE(6,*) ' TEST IS END',ind
      RETURN
    1 IND=1
c      WRITE(6,*) ' TEST IS END',ind
      RETURN    
      END   
C SUBROUTINE CMPTRV ====*====3====*====4====*====5====*====6====*====7
C
C    COMPATIBILITY TABLE ICPP OBTAINED BY COMPAT IS REFORMED TO
C    ICCPP0 BY TAKING ACOOUNT THE TIME REVERSAL SYMMETRY.
C
C                 1994/02/28
C                   BY  S.TANAKA AND A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE CMPTRV(KB,IC,KBB,ICC,IWW,ICPP,ICPP0)   
C
      DIMENSION KB(3),KBB(3),ICPP(12,12),ICPP0(12,12)
      DIMENSION ICPTR(12,12)
      DIMENSION NDEGT(12),JTT(12),IPARTT(12)
      DIMENSION NDEGP(12),JTP(12),IPARTP(12)
C
C------------------------------------------------------------
      CALL CHKDNM(IC)
      CALL CHKDNM(ICC)
      IF(IWW.NE.0.AND.IWW.NE.1) THEN
        WRITE(6,*) ' IW SHOLD BE 0 OR 1, HOWEVER IWW=',IWW
        WRITE(6,*) ' STOP IN CORRES'
        STOP
      END IF
      CALL TSIREP(KBB,ICC,IWW)
      CALL DGTRST(JDUB,NRT,MMGT,NSTRT,NDEGT,JTT,IPARTT)
      NRRT=NRT
C      WRITE(6,1100) MMGT,NSTRT,NRT
C 1100 FORMAT(' ORDER OF K-POINT GROUP=',I3,' NUMBER OF STAR=',I3
C     &      /' NUMBER OF REPRESENTATION=',I3)
C      WRITE(6,1110) (I,I=1,NRT)
C      WRITE(6,1120) (NDEGT(I),I=1,NRT)
C      WRITE(6,1130) (JTT(I),I=1,NRT)
C      WRITE(6,1140) (IPARTT(I),I=1,NRT)
C      WRITE(6,*) '                        '
C 1110 FORMAT('    NO      ',12I4)
C 1120 FORMAT(' DEGENERACY ',12I4)
C 1130 FORMAT(' HERRING SUM',12I4)
C 1140 FORMAT(' PARTNER NO ',12I4)
C-------------------------------------------------------------
      CALL TSIREP(KB,IC,IWW)
      CALL DGTRST(JDUB,NRP,MMGP,NSTRP,NDEGP,JTP,IPARTP)
      NRRP=NRP
C      WRITE(6,2100) MMGP,NSTRP,NRRP
C 2100 FORMAT(' ORDER OF K-POINT GROUP=',I3,' NUMBER OF STAR=',I3
C     &      /' NUMBER OF REPRESENTATION=',I3)
C      WRITE(6,2110) (I,I=1,NRRP)
C      WRITE(6,2120) (NDEGP(I),I=1,NRRP)
C      WRITE(6,2130) (JTP(I),I=1,NRRP)
C      WRITE(6,2140) (IPARTP(I),I=1,NRRP)
C      WRITE(6,*) '                         '
C 2110 FORMAT('    NO      ',12I4)
C 2120 FORMAT(' DEGENERACY ',12I4)
C 2130 FORMAT(' HERRING SUM',12I4)
C 2140 FORMAT(' PARTNER NO ',12I4)
      DO 1 I=1,12
      DO 2 J=1,12
      ICPTR(I,J)=0
    2 CONTINUE
    1 CONTINUE
C--------------------------------------------------------------
      IF(IWW.EQ.0) THEN
         DO 10 I=1,NRRT
            DO 20 J=1,NRRP
               IF((JTT(I).GT.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=0
               IF((JTT(I).GT.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=1
               IF((JTT(I).GT.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=2
               IF((JTT(I).EQ.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=3
               IF((JTT(I).EQ.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=4
               IF((JTT(I).EQ.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=5
               IF((JTT(I).LT.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=6
               IF((JTT(I).LT.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=7
               IF((JTT(I).LT.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=8
   20       CONTINUE
   10    CONTINUE
      ELSE
         DO 101 I=1,NRRT
            IF(JTT(I).EQ.99) JTT(I)=-99
            DO 201 J=1,NRRP
            IF(JTP(J).EQ.99) JTP(J)=-99
               IF((JTT(I).LT.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=0
               IF((JTT(I).LT.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=1
               IF((JTT(I).LT.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=2
               IF((JTT(I).EQ.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=3
               IF((JTT(I).EQ.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=4
               IF((JTT(I).EQ.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=5
               IF((JTT(I).GT.0).AND.(JTP(J).LT.0)) ICPTR(I,J)=6
               IF((JTT(I).GT.0).AND.(JTP(J).EQ.0)) ICPTR(I,J)=7
               IF((JTT(I).GT.0).AND.(JTP(J).GT.0)) ICPTR(I,J)=8
  201       CONTINUE
  101    CONTINUE
      END IF
C--------------------------------------------------------------
      DO 21 I=1,NRRT
      DO 22 J=1,NRRP
      IF(ICPTR(I,J).EQ.0) ICPP0(I,J)=ICPP(I,J)
      IF(ICPTR(I,J).EQ.1) THEN
          IF(J.LT.IPARTP(J)) THEN
            ICPP0(I,J)=(ICPP(I,J)+ICPP(I,IPARTP(J)))/2
            ICPP0(I,IPARTP(J))=ICPP0(I,J)
          ENDIF
      ELSE IF(ICPTR(I,J).EQ.2) THEN
            ICPP0(I,J)=ICPP(I,J)/2
      ELSE IF(ICPTR(I,J).EQ.3) THEN
        IF(I.LT.IPARTT(I)) THEN
            ICPP0(I,J)=ICPP(I,J)+ICPP(IPARTT(I),J)
            ICPP0(IPARTT(I),J)=ICPP0(I,J)
        END IF
      ELSE IF(ICPTR(I,J).EQ.4) THEN
        IF((I.LT.IPARTT(I)).AND.(J.LT.IPARTP(J))) THEN  
            ICPP0(I,J)=(ICPP(I,J)+ICPP(I,IPARTP(J))
     &        +ICPP(IPARTT(I),J)+ICPP(IPARTT(I),IPARTP(J)))/2
            ICPP0(IPARTT(I),J)=ICPP0(I,J)
            ICPP0(I,IPARTP(J))=ICPP0(I,J)
            ICPP0(IPARTT(I),IPARTP(J))=ICPP0(I,J)
        ENDIF
      ELSE IF(ICPTR(I,J).EQ.5) THEN
        IF(I.LT.IPARTT(I)) THEN
            ICPP0(I,J)=(ICPP(I,J)+ICPP(IPARTT(I),J))/2
            ICPP0(IPARTT(I),J)=ICPP0(I,J)
        END IF
      ELSE IF(ICPTR(I,J).EQ.6) THEN
            ICPP0(I,J)=ICPP(I,J)*2
      ELSE IF(ICPTR(I,J).EQ.7) THEN
        IF(J.LT.IPARTP(J)) THEN
            ICPP0(I,J)=ICPP(I,J)+ICPP(I,IPARTP(J))
            ICPP0(I,IPARTP(J))=ICPP0(I,J)
        END IF
      ELSE IF(ICPTR(I,J).EQ.8) THEN
            ICPP0(I,J)=ICPP(I,J)
      END IF
   22 CONTINUE
   21 CONTINUE
C--------------------------------------------------------------           
C      DO 30 J=1,NRRP
C      WRITE(6,1300) (ICPTR(I,J),I=1,NRRT)
C   30 CONTINUE
C 1300 FORMAT(' ICPTR=',12I4)
C
C      DO 40 J=1,NRRP
C      WRITE(6,1401) (ICPP(I,J),I=1,NRRT)
C   40 CONTINUE
C      DO 41 J=1,NRRP
C      WRITE(6,1400) (ICPP0(I,J),I=1,NRRT)
C   41 CONTINUE
C 1400 FORMAT(' ICPP0=',12I4)
C 1401 FORMAT(' ICPP =',12I4)
      RETURN
      END
C SUBROUTINE NEAREC ====*====3====*====4====*====5====*====6====*====7
C
C    KB(K)/IC=KBB(K)/IC+KG(K) K=1,3
C       WHERE KBB(K)/IC IS IN THE FIRST B.Z.
C       AND KG(K) IS THE RECIPROCAL LATTICE VECTOR
C
C                1994/2/20
C                   BY  A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE NEAREC(KB,IC,KBB,KG,NG)
      DIMENSION KB(3),KBB(3,10),KG(3,10),KW(3),KKG(3)
      CALL CHKDNM(IC)
      DO 1 K=1,3
         KW(K)=MOD(KB(K),IC)
         KKG(K)=(KB(K)-KW(K))/IC
    1 CONTINUE
      JJ=0
      INDF=0
    2 CALL FINDRP(KKG,INDF)
      KX=KB(1)-KKG(1)*IC
      KY=KB(2)-KKG(2)*IC
      KZ=KB(3)-KKG(3)*IC
      ICC=IC
      CALL TSKFBZ(KX,KY,KZ,ICC,IND)
c      write(6,*) kx,ky,kz,icc,ind
      IF(IND.EQ.0) GO TO 2
      IF(JJ.GT.0) THEN
         DO 4 J=1,JJ
            DO 5 K=1,3
               IF(KKG(K).NE.KG(K,J)) GO TO 4
    5       CONTINUE
            GO TO 2
    4    CONTINUE
      END IF
      JJ=JJ+1
      DO 3 K=1,3
         KG(K,JJ)=KKG(K)      
         KBB(K,JJ)=KB(K)-KKG(K)*ICC
    3 CONTINUE
      IF(IND.GT.JJ) GO TO 2
      NG=JJ
      RETURN
      END
C SUBROUTINE FINDRP ====*====3====*====4====*====5====*====6====*====7
C
C        KG(K) IS THE RECIPROCAL LATTICE VECTOR
C
C                1994/2/20
C                   BY  A. YANASE
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      SUBROUTINE FINDRP(KG,INDF)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      SAVE /SPG2/,KO,JPER,NPER,IPER,KAD
      INTEGER IPER(3,6),KAD(3,44),KG(3),KO(3)
      DATA IPER/ 1,2,3, 2,3,1, 3,1,2, 1,3,2, 3,2,1, 2,1,3/
      DATA KAD/ -1,0,0, 1,0,0, 0,-1,-1, -1,1,0, 0,1,1
     &   ,-1,-1,-1, 1,-1,-1, -1,1,1, 1,1,1, -2,0,0,  2,0,0
     &   ,3,0,0,    -3,0,0,   2,1,0,-2,1,0,  2,-1,0,-2,-1,0
     &   ,4,0,0,    -4,0,0,   2,2,0,-2,2,0, -2,-2,0  ,3,1,0
     &  ,-3,1,0,    3,-1,0, -3,-1,0 ,2,1,1,  2,-1,1,  2,-1,-1 
     &   ,-2,1,1,  -2,-1,1, -2,-1,-1 
     &   ,3,2,0, -3,2,0, 3,-2,0,-3,-2,0
     &   ,2,2,1, -2,2,1,-2,-2,1, 2,2,-1, 2,-2,-1, -2,-2,-1
     &   ,5,0,0, -5,0,0/
      IF(INDF.EQ.0) THEN
         DO 3 K=1,3
            KO(K)=KG(K)
    3    CONTINUE
         JPER=0
         NPER=1
         GO TO 1
      END IF
    2 JPER=JPER+1
      IF(JPER.EQ.7) THEN
         NPER=NPER+1
         IF(NPER.GT.44) GO TO 5
         JPER=1
      END IF
      DO 4 K=1,3
         KG(K)=KO(K)+KAD(IPER(K,JPER),NPER)
    4 CONTINUE
    1 CONTINUE
      IF(IL.EQ.-1) THEN
         IF(MOD(-KG(1)+KG(2)+KG(3),3).NE.0) GO TO 2
      ELSE IF(IL.EQ.2) THEN
         IF(MOD(IABS(KG(1)),2).NE.MOD(IABS(KG(2)),2)) GO TO 2
         IF(MOD(IABS(KG(2)),2).NE.MOD(IABS(KG(3)),2)) GO TO 2
      ELSE IF(IL.EQ.3) THEN
         IF(MOD(KG(1)+KG(2)+KG(3),2).NE.0) GO TO 2
      ELSE IF(IL.EQ.4) THEN
         IF(MOD(KG(1)+KG(2),2).NE.0) GO TO 2
      END IF
      INDF=INDF+1
      RETURN
    5 WRITE(6,*) ' STOP AT 5 IN FINDRP'
      STOP
      END
C SUBROUTINE DGTRMD ====*====3====*====4====*====5====*====6====*====7
C
C   JDUB=0    WITHOUT SPIN-ORBIT
C   JDUB=1    WITH SPIN-ORBIT
C   NUMBER OF I.R.                   NNR 
C   NUMBER OF ELEMNTS OF KP-GROUP    MMG
C   NUMBER OF STAR                   NSTAR
C   DEGENERACY OF I.R.               NDEG(12)
C   NUMBER OF TIME REVERSAL ELEMENTS MTR
C   HERING'S SUM                     JTR(12)
C   PARTNER                          IPART(12)
C   ARE OBTAINED
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE DGTRMD(JDUB,NNR,MMG,NSTR,NDEG,MTR,JTR,IPART)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,CW
      COMMON/SPG1/IT(3,48),IM(48,48),IV(48)
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &   ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      COMMON/STK/KS(3,48),JS(48),NS,ICBB,CW(48,12)
      SAVE /SPG1/,/SPG2/,/SPG3/,/SPG4/,/STK/
      DIMENSION IPART(12),NDEG(12),JTR(12)
      CALL FIDPRT(IPART)
      MTR=MTRG
      JDUB=IDOUB
      NNR=NR
      NSTR=NS
      MMG=MG
      DO 1 I=1,NR
      NDEG(I)=ND(I)
      JTR(I)=IATR(I)
    1 CONTINUE
C     WRITE(6,601) JDUB,NNR,NSTR
C 601 FORMAT(' JDUB=',I2,' NR=',I3,' NSTR=',I3)
C     WRITE(6,602) (I,I=1,NR)
C 602 FORMAT(' NO  ',12I5)
C     WRITE(6,603) (NDEG(I),I=1,NR)
C 603 FORMAT(' NDEG',12I5)
C     WRITE(6,604) (JTR(I),I=1,NR)
C 604 FORMAT(' JTR ',12I5)
C     WRITE(6,605) (IPART(I),I=1,NR)
C 605 FORMAT(' PRTN',12I5)
      RETURN
      END
C SUBROUTINE TSIRMR ====*====3====*====4====*====5====*====6====*====7
C
C  NA,NB ELEMENT OF IRREDUCIBLE REPRESENTATION IIR IS GIVEN IN WD
C  NND        :DIMENSION OF REPRESENTATION IIR
C  JJG(48)    :OPERATION CODE
C  JGA(2,3,48):TRANSLATION FOR THE CORRESPONDING JJG
C  CCR(48)    :CHARACTER OF REPRESENTATION IIR
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE TSIRMR(IIR,NA,NB,MMG,NND,JJG,JGA,CCR,WD)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,WD,CCR
      COMMON/SPG2/IL,NG,IG(48),JV(2,3,48)
      COMMON/SPG3/KB(3),ICB,MG,JG(48),JK(3,48)
     &    ,IZ,IFA(48,48),IFC,IDOUB,MTRG,JTRG(4,48),ITRC(3,48)
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG2/,/SPG3/,/SPG4/,JR,INDA,RD
      DIMENSION RD(6,6,48),JGA(2,3,48),WD(48),JJG(48),CCR(48)
      CALL CHKNIR(IIR,NR)
      CALL CHKNST(NA,ND(IIR))
      CALL CHKNST(NB,ND(IIR))
      INDT=0
      GO TO 8
C ENTRY ZZZY45 ====2====*====3====*====4====*====5====*====6====*====7
C  IF YOU ENTER HERE, ONLY THE TRANSFORMATION TO REAL MATRIX IS DONE.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY ZZZY45(IIR)
      CALL CHKNIR(IIR,NR)
      INDT=1
    8 JR=IIR
C-----------------------------------------------------------------
C   TRANSFORMATION TO REAL MATRIX, IF POSSIBLE (INDA.EQ.0)
      CALL ZZZY36(JR,RD,INDA)
C-----------------------------------------------------------------
      IF(INDT.EQ.1) RETURN
      MMG=MG
      NND=ND(IIR)
      DO 1 I=1,MG
      JJG(I)=IG(JG(I))
      CCR(I)=CR(I,IIR)
      DO 2 K=1,3
      JGA(1,K,I)=JV(1,K,JG(I))
      JGA(2,K,I)=JV(2,K,JG(I))
    2 CONTINUE
    1 CONTINUE
      GO TO 3
C ENTRY TSIRME ====2====*====3====*====4====*====5====*====6====*====7
C   SHORT FORM TO GET NA,NB ELEMENTS, AFTER YOU CALL TSIRMR ONCE.
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      ENTRY TSIRME(NA,NB,WD)
      CALL CHKNST(NA,ND(JR))
      CALL CHKNST(NB,ND(JR))
    3 IF(INDA.EQ.0) GO TO 5
      IWA=32**(6-NB)
      DO 4 I=1,MG
      IW=IR(NA,I,JR)/IWA
      IW=IW-(IW/32)*32
      IF(IW.NE.0) GO TO 7
      WD(I)=0.D0
      GO TO 4
    7 X=(DBLE(IW)/24.0D0)*2.0D0*3.14159265358979312D0
      WD(I)=DCMPLX(COS(X),SIN(X))/SQRT(DBLE(IR(7,I,JR)))
    4 CONTINUE
      RETURN
    5 DO 6 I=1,MG
    6 WD(I)=DCMPLX(RD(NA,NB,I),0.D0)
      RETURN
      END
C SUBROUTINE ZZZY36 ====*====3====*====4====*====5====*====6====*====7
C
C       MATRIX ELEMENTS OF THE I.R. ARE TANSFORMED TO REAL FORM,
C       IF POSSIBLE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      SUBROUTINE ZZZY36(JR,RD,INDA)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
      COMPLEX*16 CR,TR,CW,CD,CE,CWK
      COMMON/SPG4/NR,NH,ND(12),IR(7,48,12),CR(48,12),IATR(12)
      SAVE /SPG4/,IT,ITR
      DIMENSION TR(6,6),RD(6,6,48),CD(6,6),CE(6,6)
      INTEGER IT(2,3,3),ITR(2,2,2)
      DATA IT/1,2, 0,0, 0,0,  2,3, 1,4, 0,0,  1,4, 2,5, 3,6/
      DATA ITR/0,0,2,6, 1,7,7,1/
      IND=0
      INDC=0
      DO 1 I=1,NH
      IF(ABS(DIMAG(CR(I,JR))).GE.1.0D-4) GO TO 21
      III=IR(1,I,JR)/(32**5)
      IF(III.NE.0.AND.III.NE.24.AND.III.NE.12) INDC=1
      IF(IND.NE.0) GO TO 1
      IF(III.EQ.8) IND=1
      IF(III.EQ.6) IND=2
      IF(III.EQ.18) IND=2
    1 CONTINUE
      IF(INDC.EQ.0) GO TO 24
      IF(IND.EQ.0) IND=1
      INDB=0
      NND=ND(JR)
      N=NND/2
   23 CONTINUE
      DO 6 I=1,6
      DO 6 J=1,6
    6 TR(I,J)=0.D0
      DO 3 I=1,2
      DO 4 J=1,2
      W=3.14159265358979312D0*(DBLE(ITR(I,J,IND))/4.0D0)
      CWK=DCMPLX(COS(W),SIN(W))/DSQRT(2.0D0)
      DO 5 K=1,N
      TR(IT(I,K,N),IT(J,K,N))=CWK
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
      INDA=0
      DO 10 K=1,NH
      XW=IR(7,K,JR)
      XW=SQRT(XW)
      DO 11 I=1,NND
      DO 12 J=1,NND
      CW=0.D0
      IW=32**(6-J)
      IW=IR(I,K,JR)/IW
      IW=IW-(IW/32)*32
      IF(IW.EQ.0) GO TO 13
      W=3.14159265358979312D0*(DBLE(IW)/12.0D0)
      CW=DCMPLX(COS(W),SIN(W))/XW
   13 CD(I,J)=CW
   12 CONTINUE
   11 CONTINUE
      DO 14 I=1,NND
      DO 14 J=1,NND
      CW=0.D0
      DO 15 II=1,NND
   15 CW=CW+CD(I,II)*TR(II,J)
      CE(I,J)=CW
   14 CONTINUE
      DO 16 I=1,NND
      DO 16 J=1,NND
      CW=0.D0
      DO 17 II=1,NND
   17 CW=CW+CONJG(TR(II,I))*CE(II,J)
      CD(I,J)=CW
      RD(I,J,K)=DBLE(CW)
      IF(ABS(DIMAG(CW)).GE.1.0D-4) INDA=1
   16 CONTINUE
C     DO 18 I=1,NND
C     WRITE(6,600) (CD(I,J),J=1,NND)
C 600 FORMAT(8F8.4)
C  18 CONTINUE
   10 CONTINUE
      IF(INDA.EQ.0) RETURN
      IF(INDB.EQ.1) GO TO 22
      IF(IND.EQ.1) IND=2
      IF(IND.EQ.2) IND=1
      INDB=1
      GO TO 23
   24 INDA=2
C      WRITE(6,603) JR
C  603 FORMAT(' 1,1 ELEMENTS ARE REAL' ,I5)
      RETURN
   21 CONTINUE
C      WRITE(6,601) I,JR
C  601 FORMAT(' COMPLEX*16 CHARACTER',2I5)
      INDA=3
      RETURN
   22 CONTINUE
C      WRITE(6,602) JR
C  602 FORMAT(' WE DO NOT KNOW THE METHOD',I5)
      RETURN
      END
      SUBROUTINE DSLCLA(IR,JR,IA,KP,U,INS,AA,ND1,NND)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
      COMPLEX*16 U,CR
      COMMON/SPG4/NR,NH,ND(12),IIR(7,48,12),CR(48,12),IATR(12)
      COMMON/ATT/ISITR(LMNATM,48),KION(LMNKAT),VATOM(3,LMNATM)
     &    ,NKATOM,NATOM,MKS(11,LMNKAT),IIAA(LMNATM),JRCH
      SAVE /SPG4/,/ATT/,NDIMC,NDIMH
      DIMENSION KP(2,1),U(1),INS(4,1),AA(1)
      INTEGER NDIMC(10),NDIMH(12)
      DATA NDIMC/
     &   1,1,1,1,2,2,3,3,3,3/
      DATA NDIMH/
     &   1,1,1,1,1,1,1,1,2,2,2,2/
      IAS=1
      IF(IA.NE.1) IAS=KION(IA-1)+1
      IAL=KION(IA)
      NSI=IAL-IAS+1
      IF(JRCH.EQ.2) NDIM=NDIMH(JR)
      IF(JRCH.EQ.1) NDIM=NDIMC(JR)
      WRITE(6,663) IR,ND(IR),IA,NSI,JR,NDIM
  663 FORMAT(' SIMMETRIZED STATES FOR IR=',I2,' ND(IR)=',I2
     &       /' GENERATED FROM IA=',I2,' (NSI=',I2,') JR=',I2,
     &       ' (NDIM=',I2,')') 
      IF(ND1.EQ.0) GO TO 47
      WRITE(6,664) ND1,ND1/ND(IR)
  664 FORMAT(' NUMBER OF STATES=',I2,',NUMBER OF SET=',I2)
      DO 48 ID1=1,ND1
      WRITE(6,661) ID1,(INS(K,ID1),K=1,4),AA(ID1)
  661 FORMAT(I4,I3,' -->',I3,' MU=',I2,' NU=',I2,F10.5)
      ICOUNT=0
      DO 49 ID2=INS(1,ID1),INS(2,ID1)
      ICOUNT=ICOUNT+1
      WRITE(6,660) ID2,ICOUNT,(KP(K,ID2),K=1,2),U(ID2)
  660 FORMAT(4X,2I5,2H (,2I5,4H)  (,2F10.5,1H))
   49 CONTINUE
   48 CONTINUE
      RETURN
   47 WRITE(6,662)
  662 FORMAT(' THERE IS NO STATE')
      RETURN
      END
      SUBROUTINE DSPWA(KB,ICC,IR,NDI,KO,V,BB,D,IND,ND1,ND2,NNWW)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'prmtsp.f'
      COMPLEX*16 U,V
      COMMON/SPW/KM(4,MAXNPW),A(MAXNPW),U(MAXNPW),KT(MAXNPW),IE,IK
      SAVE /SPW/
      DIMENSION KB(3),KO(1),IND(4,1),V(1),BB(1),D(1)
      WRITE(6,663) IR,NDI,KB,ICC,A(IK),ND1,ND2,IK
  663 FORMAT(' SIMMETRIZED PLANE WAVES FOR IR=',I2,' ND(IR)=',I2
     &       /' GENERATED FOR KB=(',3I5,')/',I6,'  .LE.',F10.5, 
     &       /' NUMBER OF STATES=',I3,',NUMBER OF SET=',I3,
     &       ' NUMBER OF USED PLANE WAVES=',I4) 
      IF(ND1.EQ.0) GO TO 492
      DO 49 ID1=1,ND1
      WRITE(6,661) ID1,(IND(K,ID1),K=1,4),D(ID1),BB(ID1)
  661 FORMAT(I4,I3,' -->',I3,' MU=',I2,' NU=',I2,2F10.5)
      ICOUNT=0
      DO 48 IW=IND(1,ID1),IND(2,ID1)
      ICOUNT=ICOUNT+1
      WRITE(6,660) IW,ICOUNT,(KM(K,KO(IW)),K=1,3),ICC,V(IW)
  660 FORMAT(4X,2I5,2H (,3I5,2H)/,I3,3H  (,2F10.5,1H))
   48 CONTINUE
   49 CONTINUE
      RETURN
  492 WRITE(6,662)
  662 FORMAT(' NO STATE FOR THIS REPRESENTATION')
      RETURN
      END
