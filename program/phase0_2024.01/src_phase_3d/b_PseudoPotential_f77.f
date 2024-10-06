!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: cnvrtp, cnvrtc, xcpot, sphset2, gegant, yset2, crrsph, 
!             s101, gaunt0
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation 
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan. 
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
c $Id: b_PseudoPotential_f77.f 376 2014-06-17 07:48:31Z jkoga $
*---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine cnvrtp
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: cnvrtp
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     > (K,N1,N2,L,M,CO,R,
     <  A)
C       a(i) = sum_mm=0^m co(mm)*r(i)^(2mm+l) , for n1=<i=<n2
C***********************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER   K,N1,N2,L,M,NN,MM
      REAL*8 A(K),R(K),CO(0:*)
      DO 200 NN = N1,N2
          A(NN) =                  CO(M)
  200 CONTINUE
      DO 300 MM = M-1,0,-1
        DO 400 NN = N1,N2
          A(NN) = A(NN)*R(NN)**2 + CO(MM)
  400   CONTINUE
  300 CONTINUE
      DO 500 NN = N1,N2
          A(NN) = A(NN)*R(NN)**L
  500 CONTINUE
      RETURN
      END

*---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine cnvrtc
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: cnvrtc
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     > (K,N1,N2,L,M,CO,R,E,PHI,VLOCR,
     <  A)
C  chi(i) = sum_mm=1^m co(mm)*(2mm+2l+1)*mm*r(i)^(2mm+l-1)
C           +(E-Vloc(i)*phi(i)                  ,for n1=<i=<n2
C***********************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER   K,N1,N2,L,M,NN,MM
      REAL*8 A(K),R(K),CO(0:*),PHI(K),E,FAC,VLOCR(K)
      FAC = DFLOAT((2*M+2*L+1)*M)
      DO 200 NN = N1,N2
          A(NN) =                  CO(M) *FAC
  200 CONTINUE
      DO 300 MM = M-1,1,-1
        FAC = DFLOAT((2*MM+2*L+1)*MM)
        DO 400 NN = N1,N2
          A(NN) = A(NN)*R(NN)**2 + CO(MM)*FAC
  400   CONTINUE
  300 CONTINUE
      DO 500 NN = N1,N2
          A(NN) = A(NN)*R(NN)**(L+1)+(E-VLOCR(NN))*PHI(NN)
  500 CONTINUE
      RETURN
      END

C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine xcpot(name,xc1,xc2,rs,zeta,exc)
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: xcpot
!
!  AUTHOR(S): T. Yamasaki, Stefan Bl"ugel,   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

C********************************************************************
C   This subroutine calculates the xc-energy and potential in Rydberg
C   units. Wigner type xc-energy and potential has been added by
C                                 Stefan Bl"ugel, ISSP, Aug. 1990
c  #1) VWN potential is added          94/12/17 by T.Yamasaki
C********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      character*7 name
C     C13     = 1/3
C     C43     = 4/3
C     C213    = 2**(1/3)-1
C     C243    = 2**(4/3)-2 = 2*C213
      DATA C13 /0.333333333333333D0/
      DATA C43 /1.333333333333333D0/
      DATA C213/0.259921050D0/
      DATA C243/0.519842100D0/
      DATA WX0 /0.916330586D0/
C
      F1(X)   = ( (1.0D0+X)**C43 + (1.0D0-X)**C43 - 2.0D0 ) / C243
      F2(X)   = ( (1.0D0+X)**C13 -1.0D0 ) / C213
      G(X)    = (1.0D0+X**3)*LOG(1.0D0+1.0D0/X)
     /         - X**2 + X/2.0D0 - C13
C
      one     = 1
      two     = 2
      three   = 3
      seven   = 7
c
      WXP     = - WX0 / RS
c
      if(      name.eq.'vwn    '.or.name.eq.'VWN    '
     &     .or.name.eq.'pz     '.or.name.eq.'PZ     '
     &     .or.name.eq.'PERZUN '.or.name.eq.'perzun ' ) then
        pai=4*datan(one)
        c23  = two/three
        c76  = seven/6
        f20  = 4/(9*c213)
c
        f00 = ( (1+zeta)**c43 + (1-zeta)**c43 - 2 ) / c243             
        wx  = wxp*(1+c213*f00)                                         
        vx1 = c43*wxp*(1+zeta)**c13                                    
        vx2 = c43*wxp*(1-zeta)**c13                                    
        if(name.eq.'vwn    '.or.name.eq.'VWN    ') then
          x = sqrt(rs)                                                    
          a = -1/(3*pai*pai)
          b =  1.13107d0                                                  
          c = 13.0045d0                                                   
          d = -0.0047584d0                                                
          q = sqrt(4*c-b**2)                                           
          p = x**2 + b*x + c                                              
          p0= d**2 + b*d + c                                              
          r = datan(q/(2*x+b))                                          
          h00= (a/f20) * ( dlog((x**2)/p)+2*b*r/q
     &         -b*d*(dlog(((x-d)**2)/p)+2*(b+2*d)*r/q)/p0 )
          h10= (a/f20) * ( 2/x - (1-b*d/p0)*(2*x+b)/p -                  
     &         2*b*d/((x-d)*p0) -                                
     &         4*b*(1-(b+2*d)*d/p0)/(q**2+(2*x+b)**2) )
c
          a =  0.0621814d0                                                
          b =  3.72744d0                                                  
          c = 12.9352d0                                                   
          d = -0.10498d0                                                  
          q = sqrt(4*c-b**2)                                           
          p = x**2 + b*x + c                                              
          p0= d**2 + b*d + c                                              
          r = datan(q/(2*x+b))                                          
          h0p= a * ( dlog((x**2)/p)+2*b*r/q
     &         -b*d*(dlog(((x-d)**2)/p)+2*(b+2*d)*r/q)/p0 )
          h1p= a * ( 2/x - (1-b*d/p0)*(2*x+b)/p -                  
     &         2*b*d/((x-d)*p0) -                                
     &         4*b*(1-(b+2*d)*d/p0)/(q**2+(2*x+b)**2) )
c
          a =  0.0310907d0                                                
          b =  7.06042d0                                                  
          c = 18.0578d0                                                   
          d = -0.32500d0                                                  
          q = sqrt(4*c-b**2)                                           
          p = x**2 + b*x + c                                              
          p0= d**2 + b*d + c                                              
          r = datan(q/(2*x+b))                                          
          h0f= a * ( dlog((x**2)/p)+2*b*r/q
     &         -b*d*(dlog(((x-d)**2)/p)+2*(b+2*d)*r/q)/p0 )       
          h1f= a * ( 2/x - (1-b*d/p0)*(2*x+b)/p -                  
     &         2*b*d/((x-d)*p0) -                                
     &         4*b*(1-(b+2*d)*d/p0)/(q**2+(2*x+b)**2) )
c
          z4=zeta**4                                                      
          h000=h0f-h0p-h00                                                
          wc=h0p+h00*f00+h000*z4*f00                                      
          vc=wc-x*((1-z4*f00)*h1p+z4*f00*h1f+(1-z4)*f00*h10)/6
          g00 = 4*((1+zeta)**c13-(1-zeta)**c13) / (6*c213)
          vc0 = 4* h000*(zeta**3)*f00 + (h000*z4+h00)*g00  
          vc1=vc+(1-zeta)*vc0                                          
          vc2=vc-(1+zeta)*vc0                                          
        else
          x=rs
          z=zeta
          if(x.gt.0.9999d0) then
            gp=  -0.1423d0
            b1p=  1.0529d0
            b2p=  0.3334d0
            gf=  -0.0843d0
            b1f=  1.3981d0
            b2f=  0.2611d0
            b1px=b1p*sqrt(x)
            b2px=b2p*x
            b1fx=b1f*sqrt(x)
            b2fx=b2f*x
            wcp=gp/(1+b1px+b2px)
            wcf=gf/(1+b1fx+b2fx)
            vcp=wcp*(1+c76*b1px+c43*b2px)/(1+b1px+b2px)
            vcf=wcf*(1+c76*b1fx+c43*b2fx)/(1+b1fx+b2fx)
          else
            ap=   0.0311d0
            bp=  -0.0480d0
            cp=   0.0020d0
            dp=  -0.0116d0
            af=   0.01555d0
            bf=  -0.0269d0
            cf=   0.0007d0
            df=  -0.0048d0
            wcp=ap*dlog(x)+bp+cp*x*dlog(x)+dp*x
            wcf=af*dlog(x)+bf+cf*x*dlog(x)+df*x
            vcp=ap*dlog(x)+bp-c13*ap+c23*cp*x*dlog(x)+c13*(2*dp-cp)*x
            vcf=af*dlog(x)+bf-c13*af+c23*cf*x*dlog(x)+c13*(2*df-cf)*x
          endif
          wc=wcp+f00*(wcf-wcp)
          t1=vcp+f00*(vcf-vcp)
          f11=c23/c213*((1+z)**c13-(1-z)**c13)
          t2=f11*(wcf-wcp)
          vc1=t1+t2*(1-z)
          vc2=t1-t2*(1+z)
c-- change from hartree unit to rydberg unit
          wc=wc*two
          vc1=vc1*two
          vc2=vc2*two
c
        endif
        exc = wx + wc
        xc1 = vx1+vc1                                                     
        xc2 = vx2+vc2                                                     
        return
      endif
      WXF     =   WXP * (2.0D0**(1.0D0/3.0D0))
      VXP     =   WXP * C43
C
C------> Wigner Type xc-energy and potential
C
      IF (     NAME .EQ. 'WIGN   '.or.name.eq.'WIGNER '
     &     .or.name .eq. 'wign   '.or.name.eq.'wigner ') then
         IF ( ABS(ZETA) .GT. 1.D-3 ) THEN
             WRITE(6,*)' WIGNER XC-ENERGY AND POTENTIAL NOT SUITABLE',
     >                 ' FOR SPIN-POLARIZED CALCULATION'
             STOP 'WIGNER'
         END IF
         WCP =  - 0.88D0 / ( RS + 7.8D0 )
         VCP =    WCP * C43 + 2.288D0 / ( RS + 7.8D0 ) ** 2
         XC1 =  VXP + VCP
         XC2 =  XC1
         EXC =  WXP + WCP
         RETURN
      END IF
C
C
C------> Perdew Zunger type xc--potential
C
      IF (  NAME.EQ.'PZ     '.or.name.eq.'pz     '
     &     .or.name.eq.'PERZUN '.or.name.eq.'perzun ' ) then
        IF (RS.GT.1.0D0) THEN
          GAMP    = - 0.1432D0
          BET1P   =   1.0529D0
          BET2P   =   0.3334D0
          GAMF    = - 0.0843D0
          BET1F   =   1.3981D0
          BET2F   =   0.2611D0
          WCP     =   2.0D0*GAMP/(1.0D0+BET1P*SQRT(RS)+BET2P*RS)
          WCF     =   2.0D0*GAMF/(1.0D0+BET1F*SQRT(RS)+BET2F*RS)
          VCP     =   WCP*(1.0D0+7.0D0*BET1P*SQRT(RS)/6.0D0
     /         +C43*BET2P*RS)/(1.0D0+BET1P*SQRT(RS)+BET2P*RS)
          VCF     =   WCF*(1.0D0+7.0D0*BET1F*SQRT(RS)/6.0D0
     /         +C43*BET2F*RS)/(1.0D0+BET1F*SQRT(RS)+BET2F*RS)
        ELSE
          AP      =   0.0311D0
          BP      = - 0.048D0
          CP      =   0.0020D0
          DP      = - 0.0116D0
          AF      =   0.01555D0
          BF      = - 0.0269D0
          CF      =   0.0007D0
          DF      = - 0.0048D0
          WCP     =   2.0D0*(AP*LOG(RS)+BP+CP*RS*LOG(RS)+DP*RS)
          WCF     =   2.0D0*(AF*LOG(RS)+BF+CF*RS*LOG(RS)+DF*RS)
          VCP     =   2.0D0*(AP*LOG(RS)+(BP-C13*AP)
     /         +2.0D0*C13*CP*RS*LOG(RS)
     /         +C13*(2.0D0*DP-CP)*RS)
          VCF     =   2.0D0*(AF*LOG(RS)+(BF-C13*AF)
     /         +2.0D0*C13*CF*RS*LOG(RS)
     /         +C13*(2.0D0*DF-CF)*RS)
        END IF
      ELSE
        IF (NAME. EQ. 'BH     '.or.NAME. EQ. 'bh     ') THEN
C---*   VON BARTH-HEDIN
          CP    = 0.0504D0
          CF    = 0.0254D0
          RP    = 30.0D0
          RF    = 75.0D0
C
        ELSE IF (NAME. EQ. 'MJW    '.or.NAME. EQ. 'mjw    ') THEN
C---*   MORUZZI-JANAK-WILLIAMS
          CP    = 0.0450D0
          CF    = 0.0225D0
          RP    = 21.0D0
          RF    = 52.91668408D0
        ELSE IF (NAME .EQ. 'GL     '.or.NAME .EQ. 'gl     ') THEN
C---*   GUNNARSSON-LUNDQVIST
          CP    = 0.0504D0
          CF    = 0.0225D0
          RP    = 11.4D0
          RF    = 15.9D0
        ELSE
          WRITE (6,*) ' ==== STOP IN XCPOT. (NO DATA FOR ',NAME,') =='
          STOP ' ===== STOP IN XCOT. ==='
        END IF
C
        WCP     = - CP * G(RS/RP)
        WCF     = - CF * G(RS/RF)
        VCP     = - CP * LOG(1.0D0+RP/RS)
        VCF     = - CF * LOG(1.0D0+RF/RS)
      END IF
C
      T1      = VCF - VCP - C43*(WCF-WCP)
      T2      = C43 * ((WXF-WXP)+(WCF-WCP))
C
      FF1     = F1(ZETA)
C
C------> XC--POTENTIAL
C
      FF21    = F2(ZETA)
      FF22    = F2(-ZETA)
C     
      XC2     = VXP + VCP + T1*FF1
      XC1     = XC2 + T2*FF21
      XC2     = XC2 + T2*FF22
C
C------> XC--ENERGY
C
      EXC     = WXP + WCP + ((WXF-WXP)+(WCF-WCP))*FF1
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine sphset2
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: sphset2
!
!  AUTHOR(S): Y. Morikawa   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     >     (NFOUT,IPRI,
c f-non-locality included. 22nd Jan. 1996 Y.M
     >      lcmax,
c f-non-locality included. 22nd Jan. 1996 Y.M
     <      CR,ISPH,MMT)
c                           @(#)sphset2.f 9.1 97/05/08 14:49:23
c     index for spherical harmonics setting up.
c     extended to l=3          23rd Jan. 1996. Y.M
C************************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER  LM1(16),LM2(16),LM3(16,16,20),LMLMAX,LMAX
     &        ,I,J,K,L,II,JJ,MM
      INTEGER  L1(2),L2(2),M1(2),M2(2),M3,L3S,L3L
     &        ,NLM1,NLM2,MMT(16,16),lcmax
      INTEGER  LMT,L1T,L2T,M1T,M2T,ISPH(16,16,6),NFOUT,IPRI
      PARAMETER (LMLMAX=7)
      COMPLEX*16 AC(16,2),BC(16,2),CC(16,16,20),CCT
      REAL*8 GAUNT((LMLMAX+1)**2,(LMLMAX+1)**2,0:6),CR(16,16,6),TMP
      LMAX  = 3
c     lcmax = 4
c--%--%--%D --> Added by T.Yamasaki being suggested by Y.Morikawa
c                              at 2nd May 1994
      call rsreal((lmlmax+1)**4*7,gaunt(1,1,0))
c--%--%--%D <--
C---*---*---*---*---*---*---*---*---*---*---*---*
      CALL  GEGANT
     >                 (LMAX,NFOUT,
     <                  GAUNT)
C---*---*---*---*---*---*---*---*---*---*---*---*
      CALL  YSET2
     <     (AC,BC,LM1,LM2)
C---*---*---*---*---*---*---*---*---*---*---*---*
      DO 100 I=1,16
        L1(1) = LM1(I)/100
        M1(1) = LM1(I)-L1(1)*100-L1(1)
        L1(2) = LM2(I)/100
        M1(2) = LM2(I)-L1(2)*100-L1(2)
        DO 110 J=1,16
          L2(1) = LM1(J)/100
          M2(1) = LM1(J)-L2(1)*100-L2(1)
          L2(2) = LM2(J)/100
          M2(2) = LM2(J)-L2(2)*100-L2(2)
          MM=0
          DO 120 II=1,2
            DO 130 JJ=1,2
              L3S = ABS(L1(II)-L2(JJ))
              L3L =     L1(II)+L2(JJ)
c
c  l expantion of the charge density is truncated here.
c
  133         continue
              if(l3l.gt.lcmax) then
                 l3l = l3l-2
                 goto 133
              end if

              M3  = M1(II)-M2(JJ)
              NLM1=(L1(II)+1)**2-L1(II)+M1(II)
              NLM2=(L2(JJ)+1)**2-L2(JJ)+M2(JJ)
              DO 140 K=L3S,L3L,2
                TMP =  GAUNT(NLM1,NLM2,K)
                IF(ABS(TMP).GT.1.D-15) THEN
                    LMT = K*100+M3+K
                    CCT = BC(I,II)*AC(J,JJ)*TMP
                    IF(MM.EQ.0) GOTO 170
                    DO 150 L = 1,MM
                      IF(LMT.EQ.LM3(I,J,L)) THEN
                          CC(I,J,L) = CC(I,J,L)+CCT
                          GOTO 160
                      END IF
  150               CONTINUE
  170               MM=MM+1
                    LM3(I,J,MM) = LMT
                    CC(I,J,MM)  = CCT
  160               CONTINUE
                END IF
  140         CONTINUE
  130       CONTINUE
  120     CONTINUE
          MMT(I,J)=MM
          MM=0
          DO 200 II=1,MMT(I,J)
            IF(ABS(CC(I,J,II)).LT.1.D-15) LM3(I,J,II)=-123
            L1T = LM3(I,J,II)/100
            M1T = LM3(I,J,II)-100*L1T-L1T
            IF(L1T.LT.0) GOTO 200
            IF(M1T.EQ.0) THEN
                MM=MM+1
C---*---*---*---*---*---*---*---*---*---*---*---*
                CALL CRRSPH
     >              (L1T,M1T,CC(I,J,II),CC(I,J,II),NFOUT,
     <               ISPH(I,J,MM),CR(I,J,MM))
C---*---*---*---*---*---*---*---*---*---*---*---*
              ELSE
                DO 210 JJ=II,MMT(I,J)
                  L2T = LM3(I,J,JJ)/100
                  M2T = LM3(I,J,JJ)-100*L2T-L2T
                  IF((L1T.EQ.L2T).AND.(M1T.EQ.-M2T)) THEN
                    MM=MM+1
C---*---*---*---*---*---*---*---*---*---*---*---*
                    CALL CRRSPH
     >                (L1T,M1T,CC(I,J,II),CC(I,J,JJ),NFOUT,
     <                 ISPH(I,J,MM),CR(I,J,MM))
C---*---*---*---*---*---*---*---*---*---*---*---*
                    LM3(I,J,JJ) = -123
                    GOTO 200
                  END IF
  210           CONTINUE
            WRITE(NFOUT,*) ' THERE IS NO PAIR LM3(I,J,II)=',LM3(I,J,II)
            STOP
            END IF
  200     CONTINUE
          MMT(I,J)=MM
        IF(IPRI.GE.2) THEN
          WRITE(NFOUT,300) I,J,(CR(I,J,II),ISPH(I,J,II),II=1,MMT(I,J))
  300     FORMAT(' ',2I2,F10.6,'*ISPH(',I2,')',
     &              4(3X,F10.6,'*ISPH(',I2,')'))
        END IF
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
c                           @(#)gegant.f 9.1 97/05/08 14:48:23 
C---*-SUBROUTINE GEGANT-*----3----*----4----*----5----*----6----*----7
C
C        GENERATE THE GAUNT COEFFCIENT
C                  1986.8.8   :  A.YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine gegant
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: gegant
!
!  AUTHOR(S): A.YANASE   August/08/1986
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     >                 (LMAX,NFOUT,
     <                  GAUNT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (LMLMAX=7)
      REAL*8 GAUNT((LMLMAX+1)**2,(LMLMAX+1)**2,0:6)
      DIMENSION FAC(50)
      CALL S101(FAC,PAI4)
      DO 1 L0=0,6
      DO 2 L1=0,LMAX
      DO 3 L2=0,LMAX
      DO 4 M1=-L1,L1
      DO 5 M2=-L2,L2
      M0=M1-M2
      IF(IABS(M0).LE.L0) THEN
         CALL GAUNT0(FAC,PAI4,L1,M1,L2,M2,L0,GA,NFOUT)
         NLM1=(L1+1)**2-L1+M1
         NLM2=(L2+1)**2-L2+M2
         GAUNT(NLM1,NLM2,L0)=GA
      END IF
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine yset2
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 376 $)
!
!  SUBROUINE: yset2
!
!  AUTHOR(S): Y. Morikawa   January/22/1996
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     <     (AC,BC,LM1,LM2)
c                           @(#)yset2.f 9.1 97/05/08 14:49:44
c
c     Y_lm(i) = a(i,1)*Y_lm1(i)  + a(i,2)*Y_lm2(i)
c             = b(i,1)*Y_lm1(i)^*+ b(i,2)*Y_lm2(i)^*
c     modified from yset2.
c     extended to l=3                22nd Jan. 1996 Y. M
C***************************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER  LM1(16),LM2(16)
      COMPLEX*16 AC(16,2),BC(16,2),RR,RI,RE
      RR = DCMPLX( 1.D0/SQRT(2.D0) , 0.D0            )
      RI = DCMPLX( 0.D0            , 1.D0/SQRT(2.D0) )
      RE = DCMPLX( 1.D0            , 0.D0            )
      AC( 1,1)  =  RE*0.5D0
      AC( 2,1)  = -RR
      AC( 3,1)  =  RI
      AC( 4,1)  =  RE*0.5D0
      AC( 5,1)  =  RE*0.5D0
      AC( 6,1)  =  RR
      AC( 7,1)  = -RI
      AC( 8,1)  =  RI
      AC( 9,1)  = -RR
      AC(10,1)  =  RE*0.5d0
      AC(11,1)  = -RR
      AC(12,1)  =  RI
      AC(13,1)  =  RR
      AC(14,1)  = -RI
      AC(15,1)  = -RR
      AC(16,1)  =  RI
      AC( 1,2)  =  RE*0.5D0
      AC( 2,2)  =  RR
      AC( 3,2)  =  RI
      AC( 4,2)  =  RE*0.5D0
      AC( 5,2)  =  RE*0.5D0
      AC( 6,2)  =  RR
      AC( 7,2)  =  RI
      AC( 8,2)  =  RI
      AC( 9,2)  =  RR
      AC(10,2)  =  RE*0.5d0
      AC(11,2)  =  RR
      AC(12,2)  =  RI
      AC(13,2)  =  RR
      AC(14,2)  =  RI
      AC(15,2)  =  RR
      AC(16,2)  =  RI
      BC( 1,1)  =  RE*0.5D0
      BC( 2,1)  = -RR
      BC( 3,1)  = -RI
      BC( 4,1)  =  RE*0.5D0
      BC( 5,1)  =  RE*0.5D0
      BC( 6,1)  =  RR
      BC( 7,1)  =  RI
      BC( 8,1)  = -RI
      BC( 9,1)  = -RR
      BC(10,1)  =  RE*0.5d0
      BC(11,1)  = -RR
      BC(12,1)  = -RI
      BC(13,1)  =  RR
      BC(14,1)  =  RI
      BC(15,1)  = -RR
      BC(16,1)  = -RI
      BC( 1,2)  =  RE*0.5D0
      BC( 2,2)  =  RR
      BC( 3,2)  = -RI
      BC( 4,2)  =  RE*0.5D0
      BC( 5,2)  =  RE*0.5D0
      BC( 6,2)  =  RR
      BC( 7,2)  = -RI
      BC( 8,2)  = -RI
      BC( 9,2)  =  RR
      BC(10,2)  =  RE*0.5d0
      BC(11,2)  =  RR
      BC(12,2)  = -RI
      BC(13,2)  =  RR
      BC(14,2)  = -RI
      BC(15,2)  =  RR
      BC(16,2)  = -RI
      LM1( 1) = 000
      LM1( 2) = 102
      LM1( 3) = 102
      LM1( 4) = 101
      LM1( 5) = 202
      LM1( 6) = 204
      LM1( 7) = 204
      LM1( 8) = 203
      LM1( 9) = 203
      LM1(10) = 303
      LM1(11) = 304
      LM1(12) = 304
      LM1(13) = 305
      LM1(14) = 305
      LM1(15) = 306
      LM1(16) = 306
      LM2( 1) = 000
      LM2( 2) = 100
      LM2( 3) = 100
      LM2( 4) = 101
      LM2( 5) = 202
      LM2( 6) = 200
      LM2( 7) = 200
      LM2( 8) = 201
      LM2( 9) = 201
      LM2(10) = 303
      LM2(11) = 302
      LM2(12) = 302
      LM2(13) = 301
      LM2(14) = 301
      LM2(15) = 300
      LM2(16) = 300
      RETURN
      END

C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine crrsph
     >     (L1T,M1T,CC1,CC2,NFOUT,
     <      ISPH,CR)
c                           @(#)crrsph.f 9.1 97/05/08 14:47:58 
C*************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER  L1T,M1T,ISPH,IPM,IPMC,NFOUT
      COMPLEX*16  CC1,CC2,CCT
      REAL*8   CR,DDT
      IF(M1T.EQ.0) THEN
         IF(ABS(DIMAG(CC1)).GT.1.D-15) THEN
            WRITE(NFOUT,*) ' M1T.EQ.0 AND DIMAG(CC1).GT.0 '
            STOP
         END IF
         CR = DREAL(CC1)
         IF(L1T.EQ.0) THEN
            ISPH=1
         ELSE IF(L1T.EQ.1) THEN
            ISPH=4
         ELSE IF(L1T.EQ.2) THEN
            ISPH=5
         ELSE IF(L1T.EQ.3) THEN
            ISPH=10
         ELSE IF(L1T.EQ.4) THEN
            ISPH=17
         ELSE
            WRITE(NFOUT,*) ' L1T.GT.4 '
            STOP
         END IF
      ELSE
         IPMC = M1T/ABS(M1T)
         CCT = CC1+CC2
         DDT = ABS(CCT)
         IF(DDT.GT.1.D-15)   THEN
            IPM = 1
         ELSE IF(DDT.LT.1.D-15)  THEN
            IPM =-1
         END IF
         IF(L1T.EQ.1.AND.IPM.EQ.1) THEN
            ISPH = 3
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.1.AND.IPM.EQ.-1) THEN
            ISPH = 2
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.2.AND.ABS(M1T).EQ.1.AND.IPM.EQ.+1) THEN
            ISPH = 8
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.2.AND.ABS(M1T).EQ.1.AND.IPM.EQ.-1) THEN
            ISPH = 9
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.2.AND.ABS(M1T).EQ.2.AND.IPM.EQ.+1) THEN
            ISPH = 6
            CR   = DREAL(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.2.AND.ABS(M1T).EQ.2.AND.IPM.EQ.-1) THEN
            ISPH = 7
            CR   = DIMAG(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.1.AND.IPM.EQ.+1) THEN
            ISPH =12
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.1.AND.IPM.EQ.-1) THEN
            ISPH =11
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.2.AND.IPM.EQ.+1) THEN
            ISPH =13
            CR   = DREAL(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.2.AND.IPM.EQ.-1) THEN
            ISPH =14
            CR   = DIMAG(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.3.AND.IPM.EQ.+1) THEN
            ISPH =16
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.3.AND.ABS(M1T).EQ.3.AND.IPM.EQ.-1) THEN
            ISPH =15
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.1.AND.IPM.EQ.+1) THEN
            ISPH =19
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.1.AND.IPM.EQ.-1) THEN
            ISPH =18
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.2.AND.IPM.EQ.+1) THEN
            ISPH =20
            CR   = DREAL(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.2.AND.IPM.EQ.-1) THEN
            ISPH =21
            CR   = DIMAG(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.3.AND.IPM.EQ.+1) THEN
            ISPH =23
            CR   =-DIMAG(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.3.AND.IPM.EQ.-1) THEN
            ISPH =22
            CR   =-DREAL(CC1)*SQRT(2.D0)*IPMC
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.4.AND.IPM.EQ.+1) THEN
            ISPH =24
            CR   = DREAL(CC1)*SQRT(2.D0)
         ELSE IF(L1T.EQ.4.AND.ABS(M1T).EQ.4.AND.IPM.EQ.-1) THEN
            ISPH =25
            CR   = DIMAG(CC1)*SQRT(2.D0)*IPMC
         ELSE
            WRITE(NFOUT,*) ' L1T,M1T,IPM = ',L1T,M1T,IPM
            STOP
         END IF
      END IF
      IF(ABS(DREAL(CC1)).GT.1.d-15.AND.ABS(DIMAG(CC1)).GT.1.d-15) THEN
        WRITE(NFOUT,111) DREAL(CC1),DIMAG(CC1)
        STOP
  111   FORMAT(' ','DREAL,DIMAG IN SUB C R R S P H=',2F20.10)
      END IF
      IF(ABS(CR).LT.1.d-15) THEN
        WRITE(NFOUT,112) CR
  112   FORMAT(' ',' ABS(CR) IS 0 IN C R R S P H',F20.10)
        STOP
      END IF
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
C
      subroutine s101(fac,pai4)
c                           @(#)s101.f 9.1 97/05/08 14:49:20 
C
      IMPLICIT REAL*8(A-H,O-Z)
C  S101(N)=(N-1) FACTORIAL
C
      REAL*8 FAC(*)
      FAC(1)=1.D0
      DO 1 I=1,49
  1   FAC(I+1)=DBLE(I)*FAC(I)
      PAI4=16.0D0*DATAN(1.0D0)
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      subroutine gaunt0(fac,pai4,l1,m1,l2,m2,l,gaunt,nfout)
c                           @(#)gaunt0.f 9.1 97/05/08 14:48:23 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FAC(50)
      LL1=IABS(L1-L2)
      LL2=L1+L2
      DO 10 I=LL1,LL2,2
        IF(L.EQ.I) GO TO 20
   10 CONTINUE
      GAUNT=0.0D0
      RETURN
C
   20 CONTINUE
      X1=DBLE(L1)
      X2=DBLE(L2)
      X=DBLE(L)
      I=L1+L2+L+M2
      M3=M1-M2
      IF(IABS(M1).GT.L1 .OR. IABS(M2).GT.L2 .OR. IABS(M3).GT.L) THEN
        WRITE(NFOUT,*) 'L1,M1,L2,M2,L,M=',L1,M1,L2,M2,L,M3
        STOP ' === STOP IN SUB.GAUNT0 (INVALID M1,M2) ==='
      END IF
C     GAUNT=SQRT((2.D0*X1+1.D0)*(2.D0*X2+1.D0))*DBLE(1+4*(I/2)-2*I)
C    & *F100(FAC,L1,L2,L,M1,-M2,M3)*F102(FAC,L1,L2,L)/(2.D0*X+1.D0)
      GAUNT=SQRT((2.D0*X1+1.D0)*(2.D0*X2+1.D0)/((2.D0*X+1.D0)*PAI4))
     &     * DBLE(1+4*(I/2)-2*I)
     &     * F100(FAC,L1,L2,L,M1,-M2,M3)*F102(FAC,L1,L2,L)
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      function f100(fac,j1,j2,j3,m1,m2,m3)
c                           @(#)f100.f 9.1 97/05/08 14:48:15 
C
C      F100 : WIGNER COEFFICIENT <J1M1J2M2]J3M3>
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 FAC(*)
      IF(M3.NE.M1+M2) GO TO 2
      K1=J1+J2-J3+1
      K2=J3+J1-J2+1
      K3=J3+J2-J1+1
      K4=J1+J2+J3+2
      T=DBLE(2*J3+1)*FAC(K1)*FAC(K2)*FAC(K3)/FAC(K4)
      K1=J1+M1+1
      K2=J1-M1+1
      K3=J2+M2+1
      K4=J2-M2+1
      K5=J3+M3+1
      K6=J3-M3+1
      T=SQRT(T*FAC(K1)*FAC(K2)*FAC(K3)*FAC(K4)*FAC(K5)*FAC(K6))
      N1=MAX0(J2-J3-M1,J1-J3+M2,0)+1
      N2=MIN0(J1+J2-J3,J1-M1,J2+M2)+1
      IF(N1.GT.N2) GO TO 2
      T1=0.D0
*VDIR NOVECTOR
      DO 1 M=N1,N2
      N=M-1
      K1=J1+J2-J3-N+1
       K2=J1-M1-N+1
      K3=J2+M2-N+1
      K4=J3-J2+M1+N+1
      K5=J3-J1-M2+N+1
 1    T1=T1+DBLE(1+4*(N/2)-2*N)/(FAC(M)*FAC(K1)*FAC(K2)*FAC(K3)
     1  *FAC(K4)*FAC(K5))
      F100=T*T1
      RETURN
 2    F100=0.D0
      RETURN
      END
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
C
      function f102(fac,l1,l2,l3)
c                           @(#)f102.f 9.1 97/05/08 14:48:16 
C
C      F102 = <L1,0,L2,0]L3,0>
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER X,P
      REAL*8 FAC(*)
      LT=L1+L2+L3
      P=LT/2
      IF(2*P.NE.LT) GO TO 1
      F102=SQRT(DBLE(2*L3+1)/DBLE(LT+1))
      F102=F102*FAC(P+1)/SQRT(FAC(2*P+1))
      X=P-L1
      F102=F102*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L2
      F102=F102*SQRT(FAC(2*X+1))/FAC(X+1)
      X=P-L3
      F102=F102*SQRT(FAC(2*X+1))/FAC(X+1)
      IF(X.GT.2*(X/2)) F102=-F102
      RETURN
 1    F102=0.D0
      RETURN
      END

C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine sphset3
!=======================================================================
!
!  SOFTWARE NAME : PHASE version 2.0
!
!  SUBROUINE: sphset3
!
!  AUTHOR(S): Y. Morikawa   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     >     (NFOUT,IPRI,CR,ISPH,MMT)
c                           @(#)sphset2.f 9.1 97/05/08 14:49:23
c     index for spherical harmonics setting up.
c     extended to l=3          23rd Jan. 1996. Y.M
C************************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER  LM1(25),LM2(25),LM3(25,25,20),LMLMAX,LMAX
     &        ,I,J,K,L,II,JJ,MM
      INTEGER  L1(2),L2(2),M1(2),M2(2),M3,L3S,L3L
     &        ,NLM1,NLM2,MMT(25,25),lcmax
      INTEGER  LMT,L1T,L2T,M1T,M2T,ISPH(25,25,6),NFOUT,IPRI
      PARAMETER (LMLMAX=7)
      COMPLEX*16 AC(25,2),BC(25,2),CC(25,25,20),CCT
      REAL*8 GAUNT(25,25,0:6),CR(25,25,6),TMP
      LMAX  = 4
      lcmax = 4
c--%--%--%D --> Added by T.Yamasaki being suggested by Y.Morikawa
c                              at 2nd May 1994
      call rsreal(5**4*7,gaunt(1,1,0))
c--%--%--%D <--
C---*---*---*---*---*---*---*---*---*---*---*---*
      CALL  GEGANT2
     >                 (LMAX,NFOUT,
     <                  GAUNT)
C---*---*---*---*---*---*---*---*---*---*---*---*
      CALL  YSET3
     <     (AC,BC,LM1,LM2)
C---*---*---*---*---*---*---*---*---*---*---*---*
      DO 100 I=1,25
        L1(1) = LM1(I)/100
        M1(1) = LM1(I)-L1(1)*100-L1(1)
        L1(2) = LM2(I)/100
        M1(2) = LM2(I)-L1(2)*100-L1(2)
        DO 110 J=1,25
          L2(1) = LM1(J)/100
          M2(1) = LM1(J)-L2(1)*100-L2(1)
          L2(2) = LM2(J)/100
          M2(2) = LM2(J)-L2(2)*100-L2(2)
          MM=0
          DO 120 II=1,2
            DO 130 JJ=1,2
              L3S = ABS(L1(II)-L2(JJ))
              L3L =     L1(II)+L2(JJ)
c
c  l expantion of the charge density is truncated here.
c
  133         continue
              if(l3l.gt.lcmax) then
                 l3l = l3l-2
                 goto 133
              end if

              M3  = M1(II)-M2(JJ)
              NLM1=(L1(II)+1)**2-L1(II)+M1(II)
              NLM2=(L2(JJ)+1)**2-L2(JJ)+M2(JJ)
              DO 140 K=L3S,L3L,2
                TMP =  GAUNT(NLM1,NLM2,K)
                IF(ABS(TMP).GT.1.D-15) THEN
                    LMT = K*100+M3+K
                    CCT = BC(I,II)*AC(J,JJ)*TMP
                    IF(MM.EQ.0) GOTO 170
                    DO 150 L = 1,MM
                      IF(LMT.EQ.LM3(I,J,L)) THEN
                          CC(I,J,L) = CC(I,J,L)+CCT
                          GOTO 160
                      END IF
  150               CONTINUE
  170               MM=MM+1
                    LM3(I,J,MM) = LMT
                    CC(I,J,MM)  = CCT
  160               CONTINUE
                END IF
  140         CONTINUE
  130       CONTINUE
  120     CONTINUE
          MMT(I,J)=MM
          MM=0
          DO 200 II=1,MMT(I,J)
            IF(ABS(CC(I,J,II)).LT.1.D-15) LM3(I,J,II)=-123
            L1T = LM3(I,J,II)/100
            M1T = LM3(I,J,II)-100*L1T-L1T
            IF(L1T.LT.0) GOTO 200
            IF(M1T.EQ.0) THEN
                MM=MM+1
C---*---*---*---*---*---*---*---*---*---*---*---*
                CALL CRRSPH
     >              (L1T,M1T,CC(I,J,II),CC(I,J,II),NFOUT,
     <               ISPH(I,J,MM),CR(I,J,MM))
C---*---*---*---*---*---*---*---*---*---*---*---*
              ELSE
                DO 210 JJ=II,MMT(I,J)
                  L2T = LM3(I,J,JJ)/100
                  M2T = LM3(I,J,JJ)-100*L2T-L2T
                  IF((L1T.EQ.L2T).AND.(M1T.EQ.-M2T)) THEN
                    MM=MM+1
C---*---*---*---*---*---*---*---*---*---*---*---*
                    CALL CRRSPH
     >                (L1T,M1T,CC(I,J,II),CC(I,J,JJ),NFOUT,
     <                 ISPH(I,J,MM),CR(I,J,MM))
C---*---*---*---*---*---*---*---*---*---*---*---*
                    LM3(I,J,JJ) = -123
                    GOTO 200
                  END IF
  210           CONTINUE
            WRITE(NFOUT,*) ' THERE IS NO PAIR LM3(I,J,II)=',LM3(I,J,II)
            STOP
            END IF
  200     CONTINUE
          MMT(I,J)=MM
        IF(IPRI.GE.2) THEN
          WRITE(NFOUT,300) I,J,(CR(I,J,II),ISPH(I,J,II),II=1,MMT(I,J))
  300     FORMAT(' ',2I2,F10.6,'*ISPH(',I2,')',
     &              4(3X,F10.6,'*ISPH(',I2,')'))
        END IF
  110   CONTINUE
  100 CONTINUE
      RETURN
      END
      
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine yset3
!=======================================================================
!
!  SOFTWARE NAME : PHASE version 2.0
!
!  SUBROUINE: yset3
!
!  AUTHOR(S): Y. Morikawa   January/22/1996
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     <     (AC,BC,LM1,LM2)
c                           @(#)yset2.f 9.1 97/05/08 14:49:44
c
c     Y_lm(i) = a(i,1)*Y_lm1(i)  + a(i,2)*Y_lm2(i)
c             = b(i,1)*Y_lm1(i)^*+ b(i,2)*Y_lm2(i)^*
c     modified from yset2.
c     extended to l=3                22nd Jan. 1996 Y. M
C***************************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER  LM1(25),LM2(25)
      COMPLEX*16 AC(25,2),BC(25,2),RR,RI,RE
      RR = DCMPLX( 1.D0/SQRT(2.D0) , 0.D0            )
      RI = DCMPLX( 0.D0            , 1.D0/SQRT(2.D0) )
      RE = DCMPLX( 1.D0            , 0.D0            )
      AC( 1,1)  =  RE*0.5D0
      AC( 2,1)  = -RR
      AC( 3,1)  =  RI
      AC( 4,1)  =  RE*0.5D0
      AC( 5,1)  =  RE*0.5D0
      AC( 6,1)  =  RR
      AC( 7,1)  = -RI
      AC( 8,1)  =  RI
      AC( 9,1)  = -RR
      AC(10,1)  =  RE*0.5d0
      AC(11,1)  = -RR
      AC(12,1)  =  RI
      AC(13,1)  =  RR
      AC(14,1)  = -RI
      AC(15,1)  = -RR
      AC(16,1)  =  RI
      AC(17,1)  =  RE*0.5D0
      AC(18,1)  = -RR
      AC(19,1)  =  RI
      AC(20,1)  =  RR
      AC(21,1)  = -RI
      AC(22,1)  = -RR
      AC(23,1)  =  RI
      AC(24,1)  =  RR
      AC(25,1)  = -RI
      
      AC( 1,2)  =  RE*0.5D0
      AC( 2,2)  =  RR
      AC( 3,2)  =  RI
      AC( 4,2)  =  RE*0.5D0
      AC( 5,2)  =  RE*0.5D0
      AC( 6,2)  =  RR
      AC( 7,2)  =  RI
      AC( 8,2)  =  RI
      AC( 9,2)  =  RR
      AC(10,2)  =  RE*0.5d0
      AC(11,2)  =  RR
      AC(12,2)  =  RI
      AC(13,2)  =  RR
      AC(14,2)  =  RI
      AC(15,2)  =  RR
      AC(16,2)  =  RI
      AC(17,2)  =  RE*0.5D0
      AC(18,2)  =  RR
      AC(19,2)  =  RI
      AC(20,2)  =  RR
      AC(21,2)  =  RI
      AC(22,2)  =  RR
      AC(23,2)  =  RI
      AC(24,2)  =  RR
      AC(25,2)  =  RI
      
      BC( 1,1)  =  RE*0.5D0
      BC( 2,1)  = -RR
      BC( 3,1)  = -RI
      BC( 4,1)  =  RE*0.5D0
      BC( 5,1)  =  RE*0.5D0
      BC( 6,1)  =  RR
      BC( 7,1)  =  RI
      BC( 8,1)  = -RI
      BC( 9,1)  = -RR
      BC(10,1)  =  RE*0.5d0
      BC(11,1)  = -RR
      BC(12,1)  = -RI
      BC(13,1)  =  RR
      BC(14,1)  =  RI
      BC(15,1)  = -RR
      BC(16,1)  = -RI
      BC(17,1)  =  RE*0.5D0
      BC(18,1)  = -RR
      BC(19,1)  = -RI
      BC(20,1)  =  RR
      BC(21,1)  =  RI
      BC(22,1)  = -RR
      BC(23,1)  = -RI
      BC(24,1)  =  RR
      BC(25,1)  =  RI
      
      BC( 1,2)  =  RE*0.5D0
      BC( 2,2)  =  RR
      BC( 3,2)  = -RI
      BC( 4,2)  =  RE*0.5D0
      BC( 5,2)  =  RE*0.5D0
      BC( 6,2)  =  RR
      BC( 7,2)  = -RI
      BC( 8,2)  = -RI
      BC( 9,2)  =  RR
      BC(10,2)  =  RE*0.5d0
      BC(11,2)  =  RR
      BC(12,2)  = -RI
      BC(13,2)  =  RR
      BC(14,2)  = -RI
      BC(15,2)  =  RR
      BC(16,2)  = -RI
      BC(17,2)  =  RE*0.5D0
      BC(18,2)  =  RR
      BC(19,2)  = -RI
      BC(20,2)  =  RR
      BC(21,2)  = -RI
      BC(22,2)  =  RR
      BC(23,2)  = -RI
      BC(24,2)  =  RR
      BC(25,2)  = -RI
      
      LM1( 1) = 000
      LM1( 2) = 102
      LM1( 3) = 102
      LM1( 4) = 101
      LM1( 5) = 202
      LM1( 6) = 204
      LM1( 7) = 204
      LM1( 8) = 203
      LM1( 9) = 203
      LM1(10) = 303
      LM1(11) = 304
      LM1(12) = 304
      LM1(13) = 305
      LM1(14) = 305
      LM1(15) = 306
      LM1(16) = 306
      LM1(17) = 404
      LM1(18) = 405
      LM1(19) = 405
      LM1(20) = 406
      LM1(21) = 406
      LM1(22) = 407
      LM1(23) = 407 
      LM1(24) = 408
      LM1(25) = 408
      
      LM2( 1) = 000
      LM2( 2) = 100
      LM2( 3) = 100
      LM2( 4) = 101
      LM2( 5) = 202
      LM2( 6) = 200
      LM2( 7) = 200
      LM2( 8) = 201
      LM2( 9) = 201
      LM2(10) = 303
      LM2(11) = 302
      LM2(12) = 302
      LM2(13) = 301
      LM2(14) = 301
      LM2(15) = 300
      LM2(16) = 300
      LM2(17) = 404
      LM2(18) = 403
      LM2(19) = 403
      LM2(20) = 402
      LM2(21) = 402
      LM2(22) = 401
      LM2(23) = 401
      LM2(24) = 400
      LM2(25) = 400
      
      RETURN
      END
c                           @(#)gegant.f 9.1 97/05/08 14:48:23 
C---*-SUBROUTINE GEGANT-*----3----*----4----*----5----*----6----*----7
C
C        GENERATE THE GAUNT COEFFCIENT
C                  1986.8.8   :  A.YANASE
C
C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine gegant2
!=======================================================================
!
!  SOFTWARE NAME : PHASE version 2.0
!
!  SUBROUINE: gegant2
!
!  AUTHOR(S): A.YANASE   August/08/1986
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
!
!=======================================================================

     >                 (LMAX,NFOUT,
     <                  GAUNT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (LMLMAX=7)
      REAL*8 GAUNT(5**2,5**2,0:6)
      DIMENSION FAC(50)
      CALL S101(FAC,PAI4)
      DO 1 L0=0,6
      DO 2 L1=0,LMAX
      DO 3 L2=0,LMAX
      DO 4 M1=-L1,L1
      DO 5 M2=-L2,L2
      M0=M1-M2
      IF(IABS(M0).LE.L0) THEN
         CALL GAUNT0(FAC,PAI4,L1,M1,L2,M2,L0,GA,NFOUT)
         NLM1=(L1+1)**2-L1+M1
         NLM2=(L2+1)**2-L2+M2
         GAUNT(NLM1,NLM2,L0)=GA
      END IF
    5 CONTINUE
    4 CONTINUE
    3 CONTINUE
    2 CONTINUE
    1 CONTINUE
      RETURN
      END
