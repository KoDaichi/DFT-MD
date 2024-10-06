!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: xc_gga_rad, EXCHPBE_wb_rad, EXCHrevPBE_wb_rad, 
!             EXCHRPBE_wb_rad, CORPBE_wb_rad, NEWTO2, 
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
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
!
c
c%%%%%%%%%%%%%%%%%%%%% Header of hex.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c $Id: gncpp_xc_gga_rad.F 570 2017-04-21 20:34:50Z yamasaki $
c #1) gga introduction
c      1994/11/22 t.yamasaki
c
c     Modified by Y.Morikawa(JRCAT-NAIR), 7th Aug. 1998
c
c     A relativistic self-consistent field computer program
c     for atoms and ions using local potentials.
c     ( spin-orbit interaction is neglected. )
c     jfp = 1  spin-polarized state.
c           0  paramagnetic state.
c     made by akira hasegawa (1978).
c     ref:  d. d. koelling and b. n. harmon,
c           j. phys. c: solid state physics 10, 3107-14.
c**********************************************************************
c   input                                                             *
c        jrh  : if jrh=0 , starting charge density will be prepared.  *
c        jfp  : 0----para , 1----ferro                                *
c     i_rela  : 1: scalar rela.
c             : 2: scalar rela. small comp. neglected for valence states
c   output                                                            *
c        xe   : energy eigen value (hartree atomic unit)              *
c        rhc  : charge density of core states * r**2 * 4pi            *
c        rhv  : charge density of valence states * r**2 * 4pi         *
c        rh   : total charge density * r**2 * 4pi                     *
c        inder: 0---converged , 1---not converged                     *
c        ga   : radial wave function * r .   pl(r) in k.h.            *
c        gb   :                              ql(r) in k.h.            *
c        vr   : atomic self-consistent potential * r                  *
c**********************************************************************


C---*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine xc_gga_rad
     > (kmesh,mesh,kspin,nspin,nfout,rad,h,omo,rh,xctype
     < ,vxc,exc
     w ,grad,wwk
     & )
c
c                           @(#)gncpp_xc_gga_rad.F 1.1 03/02/19 01:30:00
c Exchange and Correlation Potential based on a Generalized Gradient
c Approximation formed by using FFTs.
c
c GGA functional:
c  John P. Perdew, Kieron Burke, and Matthias Ernzerhof,
c  Phys. Rev. Lett. vol. 77, 3865, (1996)
c
c Gradients of the soft charge density are calculated by using FFT.
c Gradients of the augmentation charge are calculated on fine grid
c meshs.
c Real space calculation is implemented by T. Sanada @ Thinking Machines
c
c              26th Feb. 1998, Y. Morikawa @ JRCAT-NAIR
C***********************************************************
      use m_Const_Parameters, only : PAI
      implicit none
      integer kmesh,mesh,kspin,nspin,ispin,nfout
     & ,i1,lgga,nr,ig,ncut
c$$$     & ,lpot
      real*8 rh(kmesh,kspin),vxc(kmesh,kspin),exc,eps_chg,eps_grad
     & ,rad(kmesh),omo(kmesh),h
     w      ,grad(kmesh,0:2),wwk(kmesh,7),pi,pi4
      character*7 xctype
      parameter(pi=PAI)
      parameter(pi4 = 4.d0*pi)
c
c set up charge cutoff and gradient cutoff
c
      eps_chg  = 1.d-25
      eps_grad = 1.d-40
c     xctype='ldapw91'
c      xctype='ggapw91'
c       xctype='ggapbe '
c
c Initialization
c
      do ispin = 1,nspin
       do nr=1,mesh
        vxc(nr,ispin) = 0.d0
       end do
      end do

      exc  = 0.d0
c$$$      lpot = 1
      if(xctype.eq.'ldapw91'.or.xctype.eq.'LDAPW91') then
        lgga=0
       else if(xctype.eq.'ggapw91'.or.xctype.eq.'ggapbe '
     &     .or.xctype.eq.'GGAPW91'.or.xctype.eq.'GGAPBE '
     &     .or.xctype.eq.'revpbe '.or.xctype.eq.'REVPBE '
     &     .or.xctype.eq.'rpbe   '.or.xctype.eq.'RPBE   ') then
        lgga=1
      end if
c
c Exchange energy and potential
c
c The spin-scaling relationship
c Ex[up,down] = 0.5*(Ex[2*up]+Ex[2*down])
c is used.
c
c     Spin loop start
c
      do ispin = 1,nspin
c
c Initialization
c
       do i1 = 0,2
        do nr=1,mesh
         grad(nr,i1) = 0.d0
        end do
       end do
c
c Form the spin density.
c
       do nr=1,mesh
        grad(nr,0) = rh(nr,ispin)/(pi4*rad(nr)**2)*dble(nspin)
       enddo
       if(xctype.eq.'ldapw91'.or.xctype.eq.'LDAPW91') then
c
c LSDA: set the gradients to zero
c
         do i1=1,2
          do nr=1,mesh
           grad(nr,i1) = 0.d0
          end do
         end do
        else if(xctype.eq.'ggapw91'.or.xctype.eq.'ggapbe '
     &      .or.xctype.eq.'GGAPW91'.or.xctype.eq.'GGAPBE '
     &      .or.xctype.eq.'revpbe '.or.xctype.eq.'REVPBE '
     &      .or.xctype.eq.'rpbe   '.or.xctype.eq.'RPBE   ') then
c
c Form the (spin) density gradients in the real space.
c
c
c The first derivatives of the (spin) density.
c
         call newto2
     >   (kmesh,mesh,rad,grad(1,0)
     <   ,grad(1,1)
     w ,wwk(1,1),wwk(1,2),wwk(1,3),wwk(1,4),wwk(1,5),wwk(1,6),wwk(1,7))
c
c Spin scaling for spin-polarized system.
c
         if(nspin.eq.2) then
           do nr=1,mesh
             grad(nr,1) = grad(nr,1)*2.d0
           end do
         end if
       end if

c$$$cdebug
c$$$      write(6,1919)(nr,grad(nr,0),grad(nr,1),nr=1400,mesh)
c$$$ 1919 format('GRDEX:',i4,2d20.10)
c$$$cdebug

c
c Form the Exchange potential and energy
c
c grad(*,0) = rho * ex
c grad(*,1) = Grad(rho)/|Grad rho| * d(fx)/d(|Grad rho|)
c
       if(xctype.eq.'rpbe   '.or.xctype.eq.'RPBE   ') then
         call EXCHRPBE_wb_rad
     >   (kmesh,mesh,lgga,eps_chg
     =   ,grad,vxc(1,ispin))
        else if(xctype.eq.'revpbe '.or.xctype.eq.'REVPBE ') then
         call EXCHrevPBE_wb_rad
     >   (kmesh,mesh,lgga,eps_chg
     =   ,grad,vxc(1,ispin))
        else
         call EXCHPBE_wb_rad
     >   (kmesh,mesh,lgga,eps_chg
     =   ,grad,vxc(1,ispin))
       endif
c
c Construct GGA correction for the exchange potential.
c
       call newto2
     > (kmesh,mesh,rad,grad(1,1)
     < ,grad(1,2)
     w ,wwk(1,1),wwk(1,2),wwk(1,3),wwk(1,4),wwk(1,5),wwk(1,6),wwk(1,7))
c$$$cdebug
c$$$      write(6,1918)(nr,vxc(nr,1),grad(nr,1),grad(nr,2),nr=1400,mesh)
c$$$ 1918 format('WBEX:',i4,3d20.10)
c$$$cdebug
c
c Form the potential
c
       do ig=1,mesh
        vxc(ig,ispin)=vxc(ig,ispin)+grad(ig,2)+2.d0/rad(ig)*grad(ig,1)
       enddo

c$$$
#ifdef GGA_CHECK
       write(6,'(" after gradient correction for ex ( xc_gga_rad )")')
       write(6,'(" rad, rho, vxc")')
       do ig = 1, mesh
         write(6,'(" ( ",i4,") rad, rho, vxc = ",3d20.8," (exch)")')
     &        ig, rad(ig), rh(ig,1)/(pi4*rad(ig)**2), vxc(ig,1)
       end do
#endif
c
c Sum up the exchange energy
c
       do ig = 1,mesh
          exc = exc + pi4*grad(ig,0)*rad(ig)**2*omo(ig)/dble(nspin)
       enddo
c
c     Spin loop end
c
      end do
c
c Correlation part.
c
c
c Initialization
c
      call rsreal(kmesh*3,grad)
c
c Form the total charge density
c
      do ispin=1,nspin
       do nr=1,mesh
        grad(nr,0) = grad(nr,0)+rh(nr,ispin)/(pi4*rad(nr)**2)
       enddo
      enddo

      if( xctype.eq.'ggapw91'.or.xctype.eq.'ggapbe '
     &.or.xctype.eq.'GGAPW91'.or.xctype.eq.'GGAPBE '
     &.or.xctype.eq.'revpbe '.or.xctype.eq.'REVPBE '
     &.or.xctype.eq.'rpbe   '.or.xctype.eq.'RPBE   ') then
c
c Form the density gradients in the real space.
c
c
c The first derivatives of the charge density.
c
         call newto2
     >   (kmesh,mesh,rad,grad(1,0)
     <   ,grad(1,1)
     w ,wwk(1,1),wwk(1,2),wwk(1,3),wwk(1,4),wwk(1,5),wwk(1,6),wwk(1,7))
      end if
#ifdef GGA_CHECK
      write(6,'(" after newto2 for cr ( xc_gga_rad )")')
      write(6,'(" rad, rho, grad(*,1)")')
      do ig = 1, mesh
         write(6,'(" ( ",i4,") rad, rho, grad = ",3d20.8," (corl)")')
     &        ig, rad(ig), grad(ig,0), grad(ig,1)
      end do
#endif
c
c Form the spin polarization.
c
      if(nspin.eq.2) then
        do ig=1,mesh
         grad(ig,2) = (rh(ig,1)-rh(ig,2))/(rh(ig,1)+rh(ig,2))
        enddo
      endif
c
c Form the Correlation potential and energy
c
      call CORPBE_wb_rad
     > (kmesh,mesh,kspin,nspin,lgga,eps_chg,eps_grad
     = ,grad,vxc)
c
c Construct GGA correction for the correlation potential.
c
       call newto2
     > (kmesh,mesh,rad,grad(1,1)
     < ,grad(1,2)
     w ,wwk(1,1),wwk(1,2),wwk(1,3),wwk(1,4),wwk(1,5),wwk(1,6),wwk(1,7))
c
c Sumup
c

c$$$cdebug
c$$$      write(6,1991)(nr,vxc(nr,1),grad(nr,2)+2.d0/rad(nr)*grad(nr,1)
c$$$     &             ,nr=1400,mesh)
c$$$ 1991 format('EXCR:',i4,2d20.10)
c$$$cdebug

c
c Form the potential
c
c$$$      write(6,'(" before gradient correction for cr ( xc_gga_rad )")')
c$$$      write(6,'(" vxc, rho, rad")')
c$$$
c$$$      do ig = 1, mesh
c$$$         write(6,'(" ( ",i4,") rad, rho, vxc = ",3d20.8," (corl)")')
c$$$     &        ig, rad(ig), rh(ig,1)/(pi4*rad(ig)**2), vxc(ig,1)
c$$$      end do
c$$$      write(6,'(" ---")')
c$$$      do ig = mesh-10+1, mesh
c$$$         write(6,'(" ( ",i4,") rad, rho, vxc = ",3d20.8)')
c$$$     &        ig, rad(ig), rh(ig,1)/(pi4*rad(ig)**2), vxc(ig,1)
c$$$      end do

      do ispin=1,nspin
       do ig=1,mesh
         vxc(ig,ispin)=vxc(ig,ispin)+grad(ig,2)+2.d0/rad(ig)*grad(ig,1)
       enddo
      enddo
#ifdef GGA_CHECK
      write(6,'(" after gradient correction for cr ( xc_gga_rad )")')
      write(6,'(" rad, rho, vxc")')
      do ig = 1, mesh
         write(6,'(" ( ",i4,") rad,rho,vxc,v4,v5 = "
     &        ,5d18.6," (corl)")')
     &        ig, rad(ig), rh(ig,1)/(pi4*rad(ig)**2), vxc(ig,1)
     &        , 2/rad(ig)*grad(ig,1), grad(ig,2)
      end do
#endif
c
c Sum up the correlation energy
c
      do ig = 1,mesh
         exc = exc + pi4*grad(ig,0)*rad(ig)**2*omo(ig)
      enddo
c
c Cut error in low density region
c
      do ispin = 1,nspin
       do nr=mesh,1,-1
        if(rh(nr,ispin)/(pi4*rad(nr)**2).gt.eps_chg) then
          ncut = nr
          goto 1001
        endif
       enddo
 1001  continue
       do nr=ncut+1,mesh
        vxc(nr,ispin) = 0.d0
       enddo
      enddo
      return
      end
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHPBE_wb_rad
     > (kmesh,mesh,lgga,eps_chg
     = ,grad,vxc)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c
c----------------------------------------------------------------------
c  Modified to adopt White and Bird scheme(PRB50, 4954 (1994))
c  for GGA potential construction. by Y. Morikawa(JRCAT-NAIR),August 3, 1998.
c
c  Modified to built into gncpp program code. Y.Morikawa August 7th, 1998.
c
c  Inputs:
c   grad(kmesh,0) : rho
c   grad(kmesh,1) : d(rho)/dr
c  Outputs:
c   grad(kmesh,0) : rho * e_x
c   grad(kmesh,1) : d(rho)/dr * partial(f_x)/partial(|Grad rho|) / |Grad rho|
c   vxc( kmesh)   : ++ partial(f_x)/partial(rho)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
cYM
      use m_Const_Parameters, only : PAI
      implicit none
      real*8 thrd,thrd4,pi,ax,um,uk,ul,exunif,rho,grad_rho,s,s2,k_f2
     & ,thrpi2,P0,FxPBE,Fs,vxg,eps_chg
cYM
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=PAI)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
c      parameter(um=0.2195149727645171d0,uk=1.2450d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
cYM
      parameter(thrpi2 = 3.d0*pi*pi)
      integer kmesh,mesh,lgga,nr
      real*8 grad(kmesh,0:*), vxc(kmesh)
cYM
#ifdef GGA_CHECK
      write(6,'(" -- EXCHPBE_wb_rad --")')
#endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      if(lgga.eq.0)then
        do nr = 1,mesh
         rho        = max(grad(nr,0),eps_chg)
         exunif     = AX*rho**THRD
         grad(nr,0) = exunif*rho
         vxc( nr)   = vxc( nr) + exunif*thrd4
#ifdef GGA_CHECK
         write(6,'(" (",i4,") rho, vxc = ",2d20.8)') nr, rho, vxc(nr)
#endif
        end do
	return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do nr=1,mesh
       rho    = max(grad(nr,0),eps_chg)
c construct LDA exchange energy density
       exunif = AX*rho**THRD
c
c Form S = ABS(Grad rho)/(2*k_F*rho)
c
       grad_rho = abs(  grad(nr,1))
       k_f2     = (thrpi2*rho)**thrd*2.d0
       s        = grad_rho/(k_f2*rho)
c
c construct PBE enhancement factor
c
       S2 = S*S
       P0=1.d0+ul*S2
       FxPBE = 1d0+uk-uk/P0
c       EX = exunif*FxPBE
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) FxPBE = 1.d0
ccheck
       grad(nr,0) =  exunif*FxPBE*rho
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c
       Fs=2.d0*uk*ul/(P0*P0)
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) Fs = 0.d0
ccheck
       vxc(nr) = vxc(nr) + THRD4*exunif*(FxPBE-FS*S2)
c
c
c
       vxg           = 0.375d0*Fs/(pi*k_f2*rho)
       grad(nr,1)    = grad(nr,1) * vxg

#ifdef GGA_CHECK
       write(6,'(" (",i4,") rho, vxc = ",2d20.8," (EXCHPBE_wb_rad)")')
     &      nr, rho, vxc(nr)
#endif

      end do


      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHrevPBE_wb_rad
     > (kmesh,mesh,lgga,eps_chg
     = ,grad,vxc)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c
c----------------------------------------------------------------------
c  Modified to adopt White and Bird scheme(PRB50, 4954 (1994))
c  for GGA potential construction. by Y. Morikawa(JRCAT-NAIR),August 3, 1998.
c
c  Modified to built into gncpp program code. Y.Morikawa August 7th, 1998.
c
c  Inputs:
c   grad(kmesh,0) : rho
c   grad(kmesh,1) : d(rho)/dr
c  Outputs:
c   grad(kmesh,0) : rho * e_x
c   grad(kmesh,1) : d(rho)/dr * partial(f_x)/partial(|Grad rho|) / |Grad rho|
c   vxc( kmesh)   : ++ partial(f_x)/partial(rho)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
cYM
      use m_Const_Parameters, only : PAI
      implicit none
      real*8 thrd,thrd4,pi,ax,um,uk,ul,exunif,rho,grad_rho,s,s2,k_f2
     & ,thrpi2,P0,FxPBE,Fs,vxg,eps_chg
cYM
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=PAI)
      parameter(ax=-0.738558766382022405884230032680836d0)
c$$$  parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
      parameter(um=0.2195149727645171d0,uk=1.2450d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
cYM
      parameter(thrpi2 = 3.d0*pi*pi)
      integer kmesh,mesh,lgga,nr
      real*8 grad(kmesh,0:*), vxc(kmesh)
cYM
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      if(lgga.eq.0)then
        do nr = 1,mesh
         rho        = max(grad(nr,0),eps_chg)
         exunif     = AX*rho**THRD
         grad(nr,0) = exunif*rho
         vxc( nr)   = vxc( nr) + exunif*thrd4
        end do
	return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do nr=1,mesh
       rho    = max(grad(nr,0),eps_chg)
c construct LDA exchange energy density
       exunif = AX*rho**THRD
c
c Form S = ABS(Grad rho)/(2*k_F*rho)
c
       grad_rho = abs(  grad(nr,1))
       k_f2     = (thrpi2*rho)**thrd*2.d0
       s        = grad_rho/(k_f2*rho)
c
c construct PBE enhancement factor
c
       S2 = S*S
       P0=1.d0+ul*S2
       FxPBE = 1d0+uk-uk/P0
c       EX = exunif*FxPBE
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) FxPBE = 1.d0
ccheck
       grad(nr,0) =  exunif*FxPBE*rho
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c
       Fs=2.d0*uk*ul/(P0*P0)
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) Fs = 0.d0
ccheck
       vxc(nr) = vxc(nr) + THRD4*exunif*(FxPBE-FS*S2)
c
c
c
       vxg           = 0.375d0*Fs/(pi*k_f2*rho)
       grad(nr,1)    = grad(nr,1) * vxg
      end do
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE EXCHRPBE_wb_rad
     > (kmesh,mesh,lgga,eps_chg
     = ,grad,vxc)
c----------------------------------------------------------------------
C  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
c  K Burke's modification of PW91 codes, May 14, 1996
c  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
c
c  Revised PBE functional after B.Hammer, L.B. Hansen and J.K. Noerskov
c
c----------------------------------------------------------------------
c  Modified to adopt White and Bird scheme(PRB50, 4954 (1994))
c  for GGA potential construction. by Y. Morikawa(JRCAT-NAIR),August 3, 1998.
c
c  Modified to built into gncpp program code. Y.Morikawa August 7th, 1998.
c
c  Inputs:
c   grad(kmesh,0) : rho
c   grad(kmesh,1) : d(rho)/dr
c  Outputs:
c   grad(kmesh,0) : rho * e_x
c   grad(kmesh,1) : d(rho)/dr * partial(f_x)/partial(|Grad rho|) / |Grad rho|
c   vxc( kmesh)   : ++ partial(f_x)/partial(rho)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  INPUT rho : DENSITY
C  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
C  INPUT U:  (GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN rho)/(rho*(2*KF)**2)
c   (for U,V, see PW86(24))
c  input lgga:  (=0=>don't put in gradient corrections, just LDA)
c  input lpot:  (=0=>don't get potential and don't need U and V)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
c [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
c     {\bf 40},  3399  (1989) (E).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c Formulas:
c   	e_x[unif]=ax*rho^(4/3)  [LDA]
c ax = -0.75*(3/pi)^(1/3)
c	e_x[PBE]=e_x[unif]*FxPBE(s)
c$$$	FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
c	FxPBE(s)=1+uk-uk*exp(-ul*s*s)
c uk, ul defined after [a](13) 
c----------------------------------------------------------------------
c----------------------------------------------------------------------
cYM
      use m_Const_Parameters, only : PAI
      implicit none
      real*8 thrd,thrd4,pi,ax,um,uk,ul,exunif,rho,grad_rho,s,s2,k_f2
     & ,thrpi2,FxPBE,Fs,vxg,eps_chg
c$$$     & ,P0
cYM
      parameter(thrd=1.d0/3.d0,thrd4=4.d0/3.d0)
      parameter(pi=PAI)
      parameter(ax=-0.738558766382022405884230032680836d0)
      parameter(um=0.2195149727645171d0,uk=0.8040d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
c      parameter(um=0.2195149727645171d0,uk=1.2450d0,ul=um/uk)
c revPBE parameter Y.M 1998/08/14
cYM
      parameter(thrpi2 = 3.d0*pi*pi)
      integer kmesh,mesh,lgga,nr
      real*8 grad(kmesh,0:*), vxc(kmesh)
cYM
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c construct LDA exchange energy density
      if(lgga.eq.0)then
        do nr = 1,mesh
         rho        = max(grad(nr,0),eps_chg)
         exunif     = AX*rho**THRD
         grad(nr,0) = exunif*rho
         vxc( nr)   = vxc( nr) + exunif*thrd4
        end do
	return
      endif
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      do nr=1,mesh
       rho    = max(grad(nr,0),eps_chg)
c construct LDA exchange energy density
       exunif = AX*rho**THRD
c
c Form S = ABS(Grad rho)/(2*k_F*rho)
c
       grad_rho = abs(  grad(nr,1))
       k_f2     = (thrpi2*rho)**thrd*2.d0
       s        = grad_rho/(k_f2*rho)
c
c construct PBE enhancement factor
c
       S2 = S*S
c$$$       P0=1.d0+ul*S2
c$$$   FxPBE = 1d0+uk-uk/P0
       FxPBE = 1d0+uk-uk*exp(-ul*s2)
c       EX = exunif*FxPBE
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) FxPBE = 1.d0
ccheck
       grad(nr,0) =  exunif*FxPBE*rho
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C  ENERGY DONE. NOW THE POTENTIAL:
c  find first and second derivatives of Fx w.r.t s.
c  Fs=(1/s)*d FxPBE/ ds
c
c$$$   Fs=2.d0*uk*ul/(P0*P0)
       Fs=2.d0*uk*ul*exp(-ul*s2)
ccheck
       if(abs(grad(nr,0)).lt.eps_chg) Fs = 0.d0
ccheck
       vxc(nr) = vxc(nr) + THRD4*exunif*(FxPBE-FS*S2)
c
c
c
       vxg           = 0.375d0*Fs/(pi*k_f2*rho)
       grad(nr,1)    = grad(nr,1) * vxg
      end do
      RETURN
      END
c----------------------------------------------------------------------
c######################################################################
c----------------------------------------------------------------------
      SUBROUTINE CORPBE_wb_rad
     > (kmesh,mesh,kspin,nspin,lgga,eps_chg,eps_grad
     = ,grad,vxc)
c----------------------------------------------------------------------
c  Official PBE correlation code. K. Burke, May 14, 1996.
C  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
C       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
C       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
C       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3)
C       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2)
C       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2
c       :  UU,VV,WW, only needed for PBE potential
c       : lgga=flag to do gga (0=>LSD only)
c       : lpot=flag to do potential (0=>energy only)
c  output: ec=lsd correlation energy from [a]
c        : vcup=lsd up correlation potential
c        : vcdn=lsd dn correlation potential
c        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
c        : dvcup=nonlocal correction to vcup
c        : dvcdn=nonlocal correction to vcdn
c----------------------------------------------------------------------
c  Modified to adopt White and Bird scheme(PRB50, 4954 (1994))
c  for GGA potential construction. by Y. Morikawa(JRCAT-NAIR),August 4, 1998.
c  Inputs:
c   grad(kmesh,0) : rho=rho_up+rho_down
c   grad(kmesh,1) : d(rho)/dr
c   grad(kmesh,2) : (rho_up-rho_down)/(rho_up+rho_dowm)
c  Outputs:
c   grad(kmesh,0) : rho * e_c
c   grad(kmesh,1) : d(rho)/dr * d(f_c)/d(|Grad rho|) / |Grad rho|
c   vxc( kmesh,ispin) : ++ partial(f_c)/partial(rho)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c References:
c [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
c     {\sl Generalized gradient approximation made simple}, sub.
c     to Phys. Rev.Lett. May 1996.
c [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
c     construction of a generalized gradient approximation:  The PW91
c     density functional}, submitted to Phys. Rev. B, Feb. 1996.
c [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
c----------------------------------------------------------------------
c----------------------------------------------------------------------
cYM
      use m_Const_Parameters, only : PAI
      implicit none
      integer kmesh,mesh,kspin,nspin,ig,lgga
      real*8 grad(kmesh,0:2), vxc( kmesh,kspin  )
cYM
c thrd*=various multiples of 1/3
c numbers for use in LSD energy spin-interpolation formula, [c](9).
c      GAM= 2^(4/3)-2
c      FZZ=f''(0)= 8/(9*GAM)
c numbers for construction of PBE
c      gamma=(1-log(2))/pi^2
c      bet=coefficient in gradient expansion for correlation, [a](4).
c      eta=small number to stop d phi/ dzeta from blowing up at 
c          |zeta|=1.
      real*8 thrd,thrdm,thrd2,sixthm,thrd4,GAM,fzz,gamma,bet,delt,eta
      parameter(thrd=1.d0/3.d0,thrdm=-thrd,thrd2=2.d0*thrd)
      parameter(sixthm=thrdm/2.d0)
      parameter(thrd4=4.d0*thrd)
      parameter(GAM=0.5198420997897463295344212145565d0)
      parameter(fzz=8.d0/(9.d0*GAM))
      parameter(gamma=0.03109069086965489503494086371273d0)
      parameter(bet=0.06672455060314922d0,delt=bet/gamma)
      parameter(eta=1.d-12)
cYM
      real*8 rs, rtrs, zet, pi,thopi,tho4pi,third,sixth,Q0,Q1,Q2,Q3
      parameter(pi=PAI)
      parameter(thopi=3.d0/pi,tho4pi=3.d0/(pi*4.d0)
     & ,third=1.d0/3.d0,sixth=1.d0/6.d0)
      real*8 Au,a1u,b1u,b2u,b3u,b4u,eu,eurs
     &      ,Ap,a1p,b1p,b2p,b3p,b4p,ep,eprs
     &      ,Aa,a1a,b1a,b2a,b3a,b4a,alfm,alfrsm
      parameter(Au  = 0.0310907D0 , a1u = 0.21370D0, b1u = 7.5957D0
     &         ,b2u = 3.5876D0    , b3u = 1.6382D0 , b4u = 0.49294D0)
      parameter(Ap  = 0.01554535D0, a1p = 0.20548D0, b1p =14.1189D0
     &         ,b2p = 6.1977D0    , b3p = 3.3662D0 , b4p = 0.62517D0)
      parameter(Aa  = 0.0168869D0 , a1a = 0.11125D0, b1a =10.357D0
     &         ,b2a = 3.6231D0    , b3a = 0.88026D0, b4a = 0.49671D0)
      real*8 rho,grad_rho,eps_grad,eps_chg,z4,f,ec,ecrs,fz,eczet
     & ,comm,g,g3,k_s,k_s2phi,T,PON,B,B2,T2,T4,Q4,Q5,H,G4,T6
     & ,RSTHRD,GZ,FAC,BG,BEC,Q8,Q9,hB,hRS,hZ,hT,pref,vcg
c$$$     & ,alfc
cYM
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c find LSD energy contributions, using [c](10) and Table I[c].
c EU=unpolarized LSD correlation energy
c EURS=dEU/drs
c EP=fully polarized LSD correlation energy
c EPRS=dEP/drs
c ALFM=-spin stiffness, [c](3).
c ALFRSM=-dalpha/drs
c F=spin-scaling factor from [c](9).
c construct ec, using [c](8)

#ifdef GGA_CHECK
      write(6,'(" -- CORPBE_wb_rad --")')
#endif

      if(lgga.eq.0) then
c
c Only LSDA
c
      do ig = 1,mesh

       rho    = max(grad(ig,0),eps_chg)
       zet    = grad(ig,2)
       rs     = (1.d0/rho*tho4pi)**third
       rtrs=dsqrt(rs)

       Q0 = -2.D0*Au*(1.D0+A1u*rtrs*rtrs)
       Q1 =  2.D0*Au*rtrs*(B1u+rtrs*(B2u+rtrs*(B3u+B4u*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       eu =  Q0*Q2
       Q3 =  Au*(B1u/rtrs+2.D0*B2u+rtrs*(3.D0*B3u+4.D0*B4u*rtrs))
       eurs = -2.D0*Au*A1u*Q2-Q0*Q3/(Q1*(1.d0+Q1))

       Q0 = -2.D0*Ap*(1.D0+A1p*rtrs*rtrs)
       Q1 =  2.D0*Ap*rtrs*(B1p+rtrs*(B2p+rtrs*(B3p+B4p*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       ep =  Q0*Q2
       Q3 =  Ap*(B1p/rtrs+2.D0*B2p+rtrs*(3.D0*B3p+4.D0*B4p*rtrs))
       eprs = -2.D0*Ap*A1p*Q2-Q0*Q3/(Q1*(1.d0+Q1))

       Q0 = -2.D0*Aa*(1.D0+A1a*rtrs*rtrs)
       Q1 =  2.D0*Aa*rtrs*(B1a+rtrs*(B2a+rtrs*(B3a+B4a*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       alfm =  Q0*Q2
       Q3 =  Aa*(B1a/rtrs+2.D0*B2a+rtrs*(3.D0*B3a+4.D0*B4a*rtrs))
       alfrsm = -2.D0*Ap*A1p*Q2-Q0*Q3/(Q1*(1.d0+Q1))

c$$$       ALFC = -ALFM
       Z4 = ZET**4
       F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
       EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
       ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
       FZ   = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
       ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
       COMM = EC -RS*ECRS/3.D0-ZET*ECZET

       grad(ig,0)     = rho*ec
       vxc( ig,1    ) = vxc(ig,1    ) + COMM + ECZET
       if(nspin.eq.2) 
     & vxc( ig,nspin) = vxc(ig,nspin) + COMM - ECZET

#ifdef GGA_CHECK
       write(6,'(" (",i4,") rho, vxc = ",2d18.8)') ig, rho, vxc(ig,1)
#endif
      enddo
      return
      endif
c----------------------------------------------------------------------

c
c LSDA + PBE Gradient Correction
c
      do ig = 1,mesh

       rho    = max(grad(ig,0),eps_chg)
       zet    = grad(ig,2)
       rs     = (1.d0/rho*tho4pi)**third
       rtrs=dsqrt(rs)

       Q0 = -2.D0*Au*(1.D0+A1u*rtrs*rtrs)
       Q1 =  2.D0*Au*rtrs*(B1u+rtrs*(B2u+rtrs*(B3u+B4u*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       eu =  Q0*Q2
       Q3 =  Au*(B1u/rtrs+2.D0*B2u+rtrs*(3.D0*B3u+4.D0*B4u*rtrs))
       eurs = -2.D0*Au*A1u*Q2-Q0*Q3/(Q1*(1.d0+Q1))

       Q0 = -2.D0*Ap*(1.D0+A1p*rtrs*rtrs)
       Q1 =  2.D0*Ap*rtrs*(B1p+rtrs*(B2p+rtrs*(B3p+B4p*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       ep =  Q0*Q2
       Q3 =  Ap*(B1p/rtrs+2.D0*B2p+rtrs*(3.D0*B3p+4.D0*B4p*rtrs))
       eprs = -2.D0*Ap*A1p*Q2-Q0*Q3/(Q1*(1.d0+Q1))

       Q0 = -2.D0*Aa*(1.D0+A1a*rtrs*rtrs)
       Q1 =  2.D0*Aa*rtrs*(B1a+rtrs*(B2a+rtrs*(B3a+B4a*rtrs)))
       Q2 =  DLOG(1.D0+1.D0/Q1)
       alfm =  Q0*Q2
       Q3 =  Aa*(B1a/rtrs+2.D0*B2a+rtrs*(3.D0*B3a+4.D0*B4a*rtrs))
       alfrsm = -2.D0*Ap*A1p*Q2-Q0*Q3/(Q1*(1.d0+Q1))

c$$$       ALFC = -ALFM
       Z4 = ZET**4
       F=((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
       EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c LSD potential from [c](A1)
c ECRS = dEc/drs [c](A2)
c ECZET=dEc/dzeta [c](A3)
c FZ = dF/dzeta [c](A4)
       ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
       FZ   = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
       ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
       COMM = EC -RS*ECRS/3.D0-ZET*ECZET

       grad(ig,0)     = rho*ec
       vxc( ig,1    ) = vxc(ig,1    ) + COMM + ECZET
       if(nspin.eq.2) 
     & vxc( ig,nspin) = vxc(ig,nspin) + COMM - ECZET
c----------------------------------------------------------------------
c PBE correlation energy
c G=phi(zeta), given after [a](3)
c DELT=bet/gamma
c B=A of [a](8)

       G  =((1.d0+ZET)**thrd2+(1.d0-ZET)**thrd2)/2.d0
       G3 = G**3

       grad_rho = max(abs(grad(ig,1)),eps_grad)
c$$$       grad_rho = abs(grad(ig,1))

       k_s      = 2.d0*(thopi*rho)**sixth
       k_s2phi  = 2.d0*k_s*g
       t  = grad_rho/(rho*k_s2phi)

       PON=-EC/(G3*gamma)
       B  = DELT/(DEXP(PON)-1.D0)
       B2 = B*B
       T2 = T*T
       T4 = T2*T2
c$$$       RS2 = RS*RS
       Q4 = 1.D0+B*T2
       Q5 = 1.D0+B*T2+B2*T4
       H = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
ccheck
       if(abs(grad(ig,0)).lt.eps_chg) H = 0.d0
ccheck
c----------------------------------------------------------------------
c----------------------------------------------------------------------
C ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
       G4 = G3*G
       T6 = T4*T2
       RSTHRD = RS/3.D0
       GZ=(((1.d0+zet)**2+eta)**sixthm-
     1     ((1.d0-zet)**2+eta)**sixthm)/3.d0
       FAC = DELT/B+1.D0
       BG = -3.D0*B2*EC*FAC/(BET*G4)
       BEC = B2*FAC/(BET*G3)
       Q8 = Q5*Q5+DELT*Q4*Q5*T2
       Q9 = 1.D0+2.D0*B*T2
       hB = -BET*G3*B*T6*(2.D0+B*T2)/Q8
       hRS = -RSTHRD*hB*BEC*ECRS
       hZ = 3.D0*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
       hT = 2.d0*BET*G3*Q9/Q8

       COMM = H+HRS-7.d0*T2*HT/6.D0
       PREF = HZ-GZ*T2*HT/G
       COMM = COMM-PREF*ZET

       vcg  = -t*hT/(grad_rho*k_s2phi)
c
c GGA energy correction
c
       grad(ig,0)     = grad(ig,0)+rho*H
c
c GGA potential correction
c
ccheck
       if(abs(grad(ig,0)).lt.eps_chg) COMM = 0.d0
       if(abs(grad(ig,0)).lt.eps_chg) PREF = 0.d0
       if(abs(grad(ig,0)).lt.eps_chg) vcg = 0.d0
ccheck
       vxc(ig,1    ) = vxc(ig,1    ) + COMM + PREF
       if(nspin.eq.2)
     & vxc(ig,nspin) = vxc(ig,nspin) + COMM - PREF
c
c GGA potential correction
c
       grad(ig,1)    = grad(ig,1) * vcg
#ifdef GGA_CHECK
       write(6,'(" (",i4,") rho, vxc = ",2d18.8)') ig, rho, vxc(ig,1)
#endif
      enddo
      RETURN
      END


      SUBROUTINE NEWTO2
     >                 (KMESH,MESH,X,A,
     <                  B,
     W                  A1,A2,A3,A4,A5,A6,A7)
C    FUNCTION A(X)
C    B(X) = DA/DX
C    7TH ORDER NEWTON'S FUNDAMENTAL FORMULA
C***********************************************************
      IMPLICIT LOGICAL(A-Z)
      INTEGER*4 KMESH,MESH,I,M
      REAL*8 X(KMESH),A(KMESH),B(KMESH),A1(KMESH),A2(KMESH),
     & A3(KMESH),A4(KMESH),A5(KMESH),A6(KMESH),A7(KMESH)
      DO 100 I=MESH-6,MESH
        A1(I) = 0.D0
        A2(I) = 0.D0
        A3(I) = 0.D0
        A4(I) = 0.D0
        A5(I) = 0.D0
        A6(I) = 0.D0
* a6 is changed to a7. 24th Apr. 1996 Y.M
c       A6(I) = 0.D0
        A7(I) = 0.D0
* a6 is changed to a7. 24th Apr. 1996 Y.M
  100 CONTINUE
      DO 200 I=1,MESH-1
        A1(I) = (A(I+1)-A(I))/(X(I+1)-X(I))
  200 CONTINUE
      DO 210 I=1,MESH-2
        A2(I) = (A1(I+1)-A1(I))/(X(I+2)-X(I))
  210 CONTINUE
      DO 220 I=1,MESH-3
        A3(I) = (A2(I+1)-A2(I))/(X(I+3)-X(I))
  220 CONTINUE
      DO 230 I=1,MESH-4
        A4(I) = (A3(I+1)-A3(I))/(X(I+4)-X(I))
  230 CONTINUE
      DO 240 I=1,MESH-5
        A5(I) = (A4(I+1)-A4(I))/(X(I+5)-X(I))
  240 CONTINUE
      DO 250 I=1,MESH-6
        A6(I) = (A5(I+1)-A5(I))/(X(I+6)-X(I))
  250 CONTINUE
      DO 260 I=1,MESH-7
        A7(I) = (A6(I+1)-A6(I))/(X(I+7)-X(I))
  260 CONTINUE
C--*--*--* EDGE POINTS
      M=MESH
      B(1) =    A1(1)
     &       + (X(1)-X(2))*(A2(1)
     &       + (X(1)-X(3))*(A3(1)
     &       + (X(1)-X(4))*(A4(1)
     &       + (X(1)-X(5))*(A5(1)
     &       + (X(1)-X(6))*(A6(1)
     &       + (X(1)-X(7))* A7(1))))))
      B(2) =    A1(1)
     &       + (X(2)-X(1))*(A2(1)
     &       + (X(2)-X(3))*(A3(1)
     &       + (X(2)-X(4))*(A4(1)
     &       + (X(2)-X(5))*(A5(1)
     &       + (X(2)-X(6))*(A6(1)
     &       + (X(2)-X(7))* A7(1))))))
      B(3) =    A1(1) + (X(3)-X(2))*A2(1)
     &       + (X(3)-X(1))*(A2(1)
     &       + (X(3)-X(2))*(A3(1)
     &       + (X(3)-X(4))*(A4(1)
     &       + (X(3)-X(5))*(A5(1)
     &       + (X(3)-X(6))*(A6(1)
     &       + (X(3)-X(7))* A7(1))))))
      B(4) =    A1(1) + (X(4)-X(2))*(A2(1)+(X(4)-X(3))*A3(1))
     &       + (X(4)-X(1))*(A2(1)+(X(4)-X(3))*A3(1)
     &       + (X(4)-X(2))*(A3(1)
     &       + (X(4)-X(3))*(A4(1)
     &       + (X(4)-X(5))*(A5(1)
     &       + (X(4)-X(6))*(A6(1)
     &       + (X(4)-X(7))* A7(1))))))
      B(M-2) = A1(M-7)
     &       + (X(M-2)-X(M-6))*(A2(M-7)
     &       + (X(M-2)-X(M-5))*(A3(M-7)
     &       + (X(M-2)-X(M-4))*(A4(M-7)
     &       + (X(M-2)-X(M-3))* A5(M-7))))
     &       + (X(M-2)-X(M-7))*(A2(M-7)
     &       + (X(M-2)-X(M-5))*(A3(M-7)
     &       + (X(M-2)-X(M-4))*(A4(M-7)
     &       + (X(M-2)-X(M-3))* A5(M-7)))
     &       + (X(M-2)-X(M-6))*(A3(M-7)
     &       + (X(M-2)-X(M-4))*(A4(M-7)
     &       + (X(M-2)-X(M-3))* A5(M-7))
     &       + (X(M-2)-X(M-5))*(A4(M-7)
     &       + (X(M-2)-X(M-3))* A5(M-7)
     &       + (X(M-2)-X(M-4))*(A5(M-7)
     &       + (X(M-2)-X(M-3))*(A6(M-7)
     &       + (X(M-2)-X(M-1))* A7(M-7))))))
      B(M-1) = A1(M-7)
     &       + (X(M-1)-X(M-6))*(A2(M-7) + (X(M-1)-X(M-5))*(A3(M-7)
     &       + (X(M-1)-X(M-4))*(A4(M-7) + (X(M-1)-X(M-3))*(A5(M-7)
     &       + (X(M-1)-X(M-2))* A6(M-7)))))
     &       + (X(M-1)-X(M-7))*(A2(M-7) + (X(M-1)-X(M-5))*(A3(M-7)
     &       + (X(M-1)-X(M-4))*(A4(M-7) + (X(M-1)-X(M-3))*(A5(M-7)
     &       + (X(M-1)-X(M-2))* A6(M-7))))
     &       + (X(M-1)-X(M-6))*(A3(M-7) + (X(M-1)-X(M-4))*(A4(M-7)
     &       + (X(M-1)-X(M-3))*(A5(M-7) + (X(M-1)-X(M-2))* A6(M-7)))
     &       + (X(M-1)-X(M-5))*(A4(M-7) + (X(M-1)-X(M-3))*(A5(M-7)
     &       + (X(M-1)-X(M-2))* A6(M-7))+ (X(M-1)-X(M-4))*(A5(M-7)
     &       + (X(M-1)-X(M-2))* A6(M-7) + (X(M-1)-X(M-3))*(A6(M-7)
     &       + (X(M-1)-X(M-2))* A7(M-7))))))
      B(M  ) = A1(M-7)
     &       + (X(M  )-X(M-6))*(A2(M-7) + (X(M  )-X(M-5))*(A3(M-7)
     &       + (X(M  )-X(M-4))*(A4(M-7) + (X(M  )-X(M-3))*(A5(M-7)
     &       + (X(M  )-X(M-2))*(A6(M-7) + (X(M  )-X(M-1))* A7(M-7))))))
     &       + (X(M  )-X(M-7))*(A2(M-7) + (X(M  )-X(M-5))*(A3(M-7)
     &       + (X(M  )-X(M-4))*(A4(M-7) + (X(M  )-X(M-3))*(A5(M-7)
     &       + (X(M  )-X(M-2))*(A6(M-7) + (X(M  )-X(M-1))* A7(M-7)))))
     &       + (X(M  )-X(M-6))*(A3(M-7) + (X(M  )-X(M-4))*(A4(M-7)
     &       + (X(M  )-X(M-3))*(A5(M-7) + (X(M  )-X(M-2))*(A6(M-7)
     &       + (X(M  )-X(M-1))* A7(M-7))))
     &       + (X(M  )-X(M-5))*(A4(M-7) + (X(M  )-X(M-3))*(A5(M-7)
     &       + (X(M  )-X(M-2))*(A6(M-7) + (X(M  )-X(M-1))* A7(M-7)))
     &       + (X(M  )-X(M-4))*(A5(M-7) + (X(M  )-X(M-2))*(A6(M-7)
     &       + (X(M  )-X(M-1))* A7(M-7))+ (X(M  )-X(M-3))*(A6(M-7)
     &       + (X(M  )-X(M-1))* A7(M-7)
     &       + (X(M  )-X(M-2))* A7(M-7))))))
C--*--*--* INNER POINTS
      DO 300 I=5,MESH-3
        B(I) = A1(I-4)
     &        +(X(I)-X(I-3))*(A2(I-4)
     &        +(X(I)-X(I-2))*(A3(I-4)
     &        +(X(I)-X(I-1))* A4(I-4)))
     &        +(X(I)-X(I-4))*(A2(I-4)
     &        +(X(I)-X(I-2))*(A3(I-4)
     &        +(X(I)-X(I-1))* A4(I-4))
     &        +(X(I)-X(I-3))*(A3(I-4)
     &        +(X(I)-X(I-1))* A4(I-4)
     &        +(X(I)-X(I-2))*(A4(I-4)
     &        +(X(I)-X(I-1))*(A5(I-4)
     &        +(X(I)-X(I+1))*(A6(I-4)
     &        +(X(I)-X(I+2))* A7(I-4))))))
*bug fixed    29th Apr. 1996 Y.M
*    &        +(X(I)-X(I-1))* A5(I-4)
*    &        +(X(I)-X(I+1))*(A6(I-4)
*    &        +(X(I)-X(I+2))* A7(I-4)))))
*bug fixed
  300 CONTINUE
      RETURN
      END

