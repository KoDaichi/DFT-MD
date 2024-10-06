!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: forcaf2
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
!  $Id: b_Force_f77.F 570 2017-04-21 20:34:50Z yamasaki $
      subroutine forcaf2
     >     (katm,katm2,kopr,natm,nopr,napt,iwei,op,ipri,af
     =     ,fnlxyz_l
     w     ,forc_up,forc_dn )
c $Id: b_Force_f77.F 570 2017-04-21 20:34:50Z yamasaki $
c The forces of the down-spin are calculated by operating
c the anti-ferromagnetic operation to those of the up-spin.
c This subroutine is called from vnlsum.
c
c antiferromagnetic calculation is added on 9th Jul. 1996
c                                          by H.Sawada
C**********************************************************************
      implicit none

#ifdef VPP
#ifdef PARA
      include 'nproc.inc'
!xocl processor pe(npe)
!xocl subprocessor pes(npe)=pe(1:npe)
      include 'inddef_p.inc'
#endif
#endif
      integer katm,katm2,natm,kopr,nopr,iopr
     &     ,iwei(katm),napt(katm,kopr)
      integer ipri,af
      integer ia,mm,nap, j
      real*8 op(3,3,kopr)
     w     ,forc_up(katm2,3),forc_dn(katm2,3),px,py,pz
c
#ifdef PARA
      real*8 fnlxyz_l(katm_p,3)
!xocl local fnlxyz_l(/ind_katm,:)
#else
      real*8 fnlxyz_l(katm,3)
#endif

      if(ipri.ge.2) then
         print *, ' forces only for up-spin (vnlsum)'
         do 610 ia = 1, natm
!xocl spread do/ind_katm
!xocl index ia
         print '(" ia=",i3,3d15.7)',ia,(fnlxyz_l(ia,j),j=1,3)
!xocl end spread
 610     continue
      endif

      do 620 ia = 1, katm2
      forc_up(ia,1)=0.d0
      forc_up(ia,2)=0.d0
      forc_up(ia,3)=0.d0
 620  continue

!xocl spread do/ind_katm
      do 500 ia = 1, natm
      do j = 1, 3
      forc_up(ia,j) = fnlxyz_l(ia,j)
      enddo
 500  continue
!xocl end spread sum(forc_up)

      mm=0
      do 300 ia = 1, natm
      if (iwei(ia).eq.2) then
         mm=mm+1
         forc_up(natm+mm,1) = - forc_up(ia,1)
         forc_up(natm+mm,2) = - forc_up(ia,2)
         forc_up(natm+mm,3) = - forc_up(ia,3)
      end if
 300  continue

      iopr = nopr +af
      do 100 ia = 1, natm
      nap = napt(ia,iopr)
      px = forc_up(nap,1)
      py = forc_up(nap,2)
      pz = forc_up(nap,3)
      forc_dn(ia,1) = op(1,1,iopr)*px+op(2,1,iopr)*py+op(3,1,iopr)*pz
      forc_dn(ia,2) = op(1,2,iopr)*px+op(2,2,iopr)*py+op(3,2,iopr)*pz
      forc_dn(ia,3) = op(1,3,iopr)*px+op(2,3,iopr)*py+op(3,3,iopr)*pz
 100  continue

!xocl spread do/ind_katm
      do 400 ia = 1, natm
      do j = 1, 3
      fnlxyz_l(ia,j) = fnlxyz_l(ia,j) + forc_dn(ia,j)
      enddo
 400  continue
!xocl end spread

      if(ipri.ge.2) then
         print *, ' forces for up- and down-spin (vnlsum)'
         do 600 ia = 1, natm
!xocl spread do/ind_katm
!xocl index ia
         print '(" ia=",i3,3d15.7)',ia,(fnlxyz_l(ia,j),j=1,3)
!xocl end spread
 600     continue
      endif

      return
      end
