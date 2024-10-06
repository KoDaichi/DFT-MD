!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: anlmes packuc initid setowh initai mvnbr sortsa sortsa_2
!             bsort dsort swapio getmin getmin_2 caldst calttl distbs
!             patmbr gtrcmt gtbrmt posatm gtnabc gtwscr wdmesp dwinly2
!             exstak
!
!  AUTHOR(S): T. Yamasaki   February/08/2004
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
! $Id: b_Ldos_f77.F 633 2020-12-01 05:11:03Z jkoga $
C   *********************************
C     by T.Yamasaki
C           10th Jun 1992
C     revised 14th Nov 1992
C   *********************************
C --------------------------------------------
C
C       modified by T. Yamasaki, 17th Nov 1992 
C       modified by T. Yamasaki, 18th Jan 2004
C       modified by T. Yamasaki, 08th Feb 2004
C
C         this subroutine is specified for the case of the system
C        which has the inversion symmetry.
C
C
      subroutine anlmes_ac(ista_atm,iend_atm,ista_atm2,iend_atm2,nfout
     &              ,ipri,nmesh,ac_mesh 
     &              ,nac_mesh,ifxrh,ifyrh,ifxr,ifyr,ifzr
     &              ,altv,rltv,catxyz,katm2,crtdst_aldos,delta_dist
     &              ,catoms,ioddst,dstnc,cnmesh,nmmax)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit none
      integer ista_atm,iend_atm,ista_atm2,iend_atm2
      integer nfout,ipri,nmesh,ifxrh, ifyrh, ifxr, ifyr, ifzr, katm2
      integer ac_mesh(ista_atm2:iend_atm2,nmesh,3)
      integer nac_mesh(ista_atm2:iend_atm2)
      real*8 altv(3,3), rltv(3,3)
      real*8 catxyz(katm2,3),catoms(katm2,3)
      real*8 crtdst_aldos, delta_dist
      integer ioddst(katm2+1, *)
      real*8 dstnc(katm2)
      logical cnmesh
      integer nmmax
c
      integer INITIA, IAAX, IBAX, ICAX
      parameter (INITIA = 1)
      parameter (IAAX = 1, IBAX = 2, ICAX = 3 )
      real*8 ttl(6)

      integer nabc(3), nmabc(3), onwhca(3,2)
      real*8 abcmes(3)
      real*8 dface(3),dline(2,3)
      real*8 offset(3)

      real*8 dx,dy,dz
      real*8 cella,cellb,cellc
      integer ix,iy,iz

      real*8 wscr
      real*8 dm
      integer i,  ia, iatom, icount, ii
      integer :: ixx,iyy,izz,iix,iiy,iiz
      real*8, allocatable, dimension(:,:) :: catomst

      allocate(catomst(katm2,3))
      wscr = crtdst_aldos
      dx = 1.d0/dble(ifxr)
      dy = 1.d0/dble(ifyr)
      dz = 1.d0/dble(ifzr)
      cella = dx*sqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
      cellb = dy*sqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
      cellc = dz*sqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
      ixx=nint(wscr/cella)
      iyy=nint(wscr/cellb)
      izz=nint(wscr/cellc)
      if(mod(ixx,2)/=0) ixx = ixx+1
      if(mod(iyy,2)/=0) iyy = iyy+1
      if(mod(izz,2)/=0) izz = izz+1
      if(ixx.gt.ifxr/2) ixx=ifxr/2
      if(iyy.gt.ifyr/2) iyy=ifyr/2
      if(izz.gt.ifzr/2) izz=ifzr/2
      if(ipri>=2) write(nfout,'(a,7f10.5)') 'wscr cellabc w/c '
     & ,wscr,cella,cellb,cellc,wscr/cella,wscr/cellb,wscr/cellc
      if(ipri>=2) write(nfout,'(a,3i8)') 
     &            ' !** ix   iy   iz   ',ixx,iyy,izz
      if(ipri>=2) write(nfout,'(a,3i8)') 
     &            ' !** ifxr ifyr ifzr ',ifxr,ifyr,ifzr
      call calttl(ipri,altv,ttl)
      call distbs(nfout,ipri,altv,rltv,dface,dline)
C
      call posatm(catxyz,katm2,katm2,rltv,catoms)
c
      call packuc(catoms,katm2,katm2)

      nmmax = 0
      do ia=ista_atm,iend_atm
         icount = 0
         do ix=-ixx,ixx
            iix = ix
            if(ix<=0) iix = iix+ifxr
            abcmes(1) = ix*dx
            do iy=-iyy,iyy
               iiy = iy
               if(iy<=0) iiy = iiy+ifyr
               abcmes(2) = iy*dy
               do iz=-izz,izz
                  iiz = iz
                  if(iz<=0) iiz = iiz+ifzr
                  abcmes(3) = iz*dz
                  offset(:) = catoms(ia,:)
                  do i=1,3
                     if(abcmes(i)>1) abcmes(i) = abcmes(i)-1
                     if(abcmes(i)<0) abcmes(i) = abcmes(i)+1
                  enddo
                  do i=1,3
                  catomst(1:katm2,i) = catoms(1:katm2,i)-offset(i)
                  enddo
                  call packuc(catomst,katm2,katm2)
                  call getmin_2(nfout,catomst,katm2,abcmes,ttl,wscr
     &                ,delta_dist,iatom)
                  if(iatom.eq.ia ) then
                     icount = icount+1
                     if(.not.cnmesh)then
                        do i=1,3
                           ac_mesh(ia,icount,1) = iix
                           ac_mesh(ia,icount,2) = iiy
                           ac_mesh(ia,icount,3) = iiz
                        enddo
                     endif
                  endif
               enddo
            enddo
         enddo
         if(icount>nmmax) then
            nmmax = icount
         endif
         if(.not.cnmesh)then
            nac_mesh(ia) = icount
         endif
      enddo
      deallocate(catomst)
c
      return
      end
C --------------------------------------------------------

      subroutine anlmes(nfout,ipri,mesh,ifxrh,ifyrh,ifxr,ifyr,ifzr
     &              ,altv,rltv,catxyz,katm2,crtdst_aldos,delta_dist
     &              ,catoms,ioddst,dstnc)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit none
      integer nfout,ipri,ifxrh, ifyrh, ifxr, ifyr, ifzr, katm2
      integer mesh(ifxrh,ifyrh,ifzr)
      real*8 altv(3,3), rltv(3,3)
      real*8 catxyz(katm2,3),catoms(katm2,3)
      real*8 crtdst_aldos, delta_dist
      integer ioddst(katm2+1, *)
      real*8 dstnc(katm2)
c
      integer INITIA, IAAX, IBAX, ICAX
      parameter (INITIA = 1)
      parameter (IAAX = 1, IBAX = 2, ICAX = 3 )
      real*8 ttl(6)
c$$$      dimension nabc(3), abcmes(3), onwhca(3,2), nmabc(3)
      integer nabc(3), nmabc(3), onwhca(3,2)
      real*8 abcmes(3)
      real*8 dface(3),dline(2,3)
c$$$      common/bdist/dface(3),dline(2,3)
c$$$      common/bcunt/icount_eq, icount_h, icount_1, icount_2

      real*8 wscr, dx1,dx2,dx3
      integer i, n1ax, n2ax, n3ax, ic, ib, ia, iatom, i_all, isp, icount

      wscr = crtdst_aldos

      if(ipri >= 1) 
     &     write(nfout,'(" !!ldos, dimension of mesh = ",3i8)')
     &          ifxrh, ifyrh, ifzr

      call calttl(ipri,altv,ttl)
      call distbs(nfout,ipri,altv,rltv,dface,dline)
C
      call posatm(catxyz,katm2,katm2,rltv,catoms)
c
      call packuc(catoms,katm2,katm2)
      if(ipri >= 1) then
         write(nfout,'(" !!ldos *** packuc ***")')
         do i = 1, katm2
            write(nfout,9003) i
     &           ,catxyz(i,1), catxyz(i,2), catxyz(i,3)
     &           ,catoms(i,1),catoms(i,2),catoms(i,3)
 9003       format(i4,2('(',3f16.8,')'))
         end do
      end if
C
      call gtnabc(nabc,ifxr,ifyr,ifzr)
      if(ipri >= 2) then
         write(nfout,*) ' !!ldos nabc(1 - 3) = ',nabc(1),nabc(2),nabc(3)
         write(nfout,*) ' !!ldos *** gtnabc ***'
      end if

c$$$      call gtwscr(nfout,altv,wscr)
      if(ipri >= 1) 
     &     write(nfout,*) ' !!ldos *** gtwscr *** , wscr = ', wscr
C
      n1ax = ICAX
      dx1 = 1.0/dble(nabc(n1ax))
      n2ax = IBAX
      dx2 = 1.0/dble(nabc(n2ax))
      n3ax = IAAX
      dx3 = 1.0/dble(ifxr)
C
      if(ipri >= 2) 
     &     write(nfout,'(" !!ldos: dx1,n1ax,dx2,n2ax,dx3,n3ax = "
     &                  ,3(f8.4,i4))') dx1,n1ax,dx2,n2ax,dx3,n3ax
C
#ifndef _FAST_WAY_
      do 1001 ic = 1, nabc(n1ax)
         nmabc(n1ax) = ic
         abcmes(n1ax) = (ic-1)*dx1
         do 1002 ib = 1, nabc(n2ax)
            nmabc(n2ax) = ib
            abcmes(n2ax) = (ib-1)*dx2
            do 1003 ia = 1, nabc(n3ax)
               nmabc(n3ax) = ia
               abcmes(n3ax) = (ia-1)*dx3
#ifdef DEBUG
               write(nfout,'(" !!ldos ia,ib,ic = ", 3i5)') ia,ib,ic
#endif
               call getmin_2(nfout,catoms, katm2, abcmes, ttl
     &              ,wscr, delta_dist
     &              ,iatom)
               mesh(nmabc(1),nmabc(2),nmabc(3)) = iatom
#ifdef DEBUG
               write(nfout,9005) ia, ib, ic, iatom
     &              ,abcmes(n1ax),abcmes(n2ax),abcmes(n3ax)
     &              ,catoms(iatom,1),catoms(iatom,2),catoms(iatom,3)
#endif
 1003       continue
 1002    continue
 1001 continue

#else
      call rsreal(ICAX, abcmes)

      call initid(nfout,ioddst(1,INITIA), katm2+1, katm2)

      call setowh(onwhca, n1ax, n2ax)

      do 1 ic = 1, nabc(n1ax)
         nmabc(n1ax) = ic
         abcmes(n1ax) = (ic-1)*dx1
c         write(nfout,*) ' *** ( ic = ',ic,' ) ***'
         call mvnbr(catoms(1,n1ax), katm2,abcmes(n1ax), 
     &        ioddst(1,INITIA), katm2+1)
c$$$         call sortsa(nfout,catoms(1,1),katm2,abcmes,ttl,onwhca(1,1),
c$$$     &               wscr,ioddst(1,INITIA),katm2+1,
c$$$     &               ioddst(1,INITIA+n1ax),dstnc)
C$$ASASASASAS
C$$$         call sortsa_2(nfout,catoms(1,1),katm2,abcmes,ttl,onwhca(1,1)
C$$$    &        ,wscr,dface,dline,ioddst(1,INITIA),katm2+1
C$$$    &        ,ioddst(1,INITIA+n1ax),dstnc)
         call sortsa_2(nfout,ipri,catoms(1,1),katm2,abcmes,ttl   !USAMI
     &        ,onwhca(1,1),wscr,dface,dline,ioddst(1,INITIA)
     &        ,katm2+1,ioddst(1,INITIA+n1ax),dstnc)
C$$ASASASASAS
C
         do 2 ib = 1, nabc(n2ax)
C           write(nfout,*) '   ***** ( ib = ', ib,' ) *****'
            nmabc(n2ax) = ib
            abcmes(n2ax) = (ib-1)*dx2
            call mvnbr(catoms(1,n2ax), katm2, abcmes(n2ax),
     &           ioddst(1,INITIA+n1ax),katm2+1)
c$$$            call sortsa(nfout,catoms(1,1),katm2,abcmes,ttl,onwhca(1,2),
c$$$     &                  wscr,ioddst(1,INITIA+n1ax),katm2+1,
c$$$     &                  ioddst(1,INITIA+n2ax),dstnc)
C$$ASASASASAS            call sortsa_2(nfout,catoms(1,1),katm2,abcmes,ttl,onwhca(1,2)
C$$ASASASASAS     &           ,wscr,dface,dline,ioddst(1,INITIA+n1ax),katm2+1
C$$ASASASASAS     &           ,ioddst(1,INITIA+n2ax),dstnc)
            call sortsa_2(nfout,ipri,catoms(1,1),katm2,abcmes,ttl !USAMI
     &           ,onwhca(1,2),wscr,dface,dline,ioddst(1,INITIA+n1ax)
     &           ,katm2+1,ioddst(1,INITIA+n2ax),dstnc)
C$$ASASASASAS
            do 3 ia = 1, nabc(n3ax)
               nmabc(n3ax) = ia
               abcmes(n3ax) = (ia-1)*dx3
               call getmin(nfout,n3ax,catoms(1,1), katm2, abcmes, ttl,
     &              wscr,ioddst(1,INITIA+n2ax),
     &              katm2+1,dstnc, iatom)
               mesh(nmabc(1),nmabc(2),nmabc(3)) = iatom
#ifdef DEBUG
               write(nfout,9005) ia, ib, ic, iatom
     &              ,abcmes(n1ax),abcmes(n2ax),abcmes(n3ax)
     &              ,catoms(iatom,1),catoms(iatom,2),catoms(iatom,3)
 9005          format('(',3i3,') (',i4,',',3f6.3,' )',3f6.3)
#endif
 3          continue
 2       continue
 1    continue
#endif
c
      if(ipri >= 3) then
         i_all = 0
         do 2001 isp = 1, katm2 + 1
            icount = 0
            do 2002 ic = 1, nabc(3)
               do 2003 ib = 1, nabc(2)
                  do 2004 ia = 1, nabc(1)
                     if(mesh(ia,ib,ic).eq.isp) icount = icount + 1
 2004             continue
 2003          continue
 2002       continue
            write(nfout,*) ' !D (isp ) ',isp,' : ', icount
            i_all = i_all + icount
 2001    continue
         write(nfout,*) ' !D i_all = ', i_all
#ifdef DEBUG_MESH
         write(nfout,*) ' !D mesh '
         do 11 ic = 1, nabc(3)
            do 12 ib = 1, nabc(2)
               write(nfout,9002) (mesh(ia,ib,ic),ia = 1, nabc(1))
 12         continue
 11      continue
 9002    format(16i4)
#endif
      end if
      return
      end
C --------------------------------------------------------
C
C        anlmesl: coded by T.Yamasaki 6th Apr. 2017
C 
C --------------------------------------------------------
      subroutine anlmesl(nfout,ipri,meshl,ifxrh,ifyrh,ifxr,ifyr,ifzr
     &              ,winlay,mlayer,normal_axis,nmeshlay)
      implicit none
      integer nfout,ipri,ifxrh, ifyrh, ifxr, ifyr, ifzr, mlayer
      integer normal_axis
      integer meshl(ifxrh,ifyrh,ifzr),nmeshlay(mlayer)
      real*8  winlay(mlayer,2)
c
      integer IAAX, IBAX, ICAX
      parameter (IAAX = 1, IBAX = 2, ICAX = 3 )
      integer nabc(3), nmabc(3)
      real*8 abcmes(3)

      real*8  dx1,dx2,dx3
      integer i, n1ax, n2ax, n3ax, ic, ib, ia, i_all, isp, icount,ilayer

      if(ipri >= 1) 
     &     write(nfout,'(" !!ldos, dimension of meshl = ",3i8)')
     &          ifxrh, ifyrh, ifzr

      call gtnabc(nabc,ifxr,ifyr,ifzr)

      n1ax = ICAX
      dx1 = 1.0/dble(nabc(n1ax))
      n2ax = IBAX
      dx2 = 1.0/dble(nabc(n2ax))
      n3ax = IAAX
      dx3 = 1.0/dble(ifxr)
C
      if(ipri >= 2) 
     &     write(nfout,'(" !!ldos: dx1,n1ax,dx2,n2ax,dx3,n3ax = "
     &                  ,3(f8.4,i4))') dx1,n1ax,dx2,n2ax,dx3,n3ax
C
      do 1001 ic = 1, nabc(n1ax)
         nmabc(n1ax) = ic
         abcmes(n1ax) = (ic-1)*dx1
         do 1002 ib = 1, nabc(n2ax)
            nmabc(n2ax) = ib
            abcmes(n2ax) = (ib-1)*dx2
            do 1003 ia = 1, nabc(n3ax)
               nmabc(n3ax) = ia
               abcmes(n3ax) = (ia-1)*dx3
#ifdef DEBUG
               write(nfout,'(" !!ldos ia,ib,ic = ", 3i5)') ia,ib,ic
#endif
               call getmin_layer(nfout,winlay, mlayer
     &              ,normal_axis,abcmes,ilayer)
               meshl(nmabc(1),nmabc(2),nmabc(3)) = ilayer
#ifdef DEBUG
               write(nfout,9005) ia, ib, ic, ilayer
     &              ,abcmes(n1ax),abcmes(n2ax),abcmes(n3ax)
 9005          format('(',3i3,') (',i4,',',3f6.3,' )',3f6.3)
#endif
 1003       continue
 1002    continue
 1001 continue

      i_all = 0
      do 2001 isp = 1, mlayer
         icount = 0
         do 2002 ic = 1, nabc(3)
            do 2003 ib = 1, nabc(2)
               do 2004 ia = 1, nabc(1)
                  if(meshl(ia,ib,ic).eq.isp) icount = icount + 1
 2004          continue
 2003       continue
 2002    continue
         if(ipri >= 1) write(nfout,*) ' !D (isp ) ',isp,' : ', icount
         i_all = i_all + icount
         nmeshlay(isp) = icount
 2001 continue
         
      if(ipri>=1) then
         write(nfout,*) ' !D i_all = ', i_all
#ifdef DEBUG_MESH
         write(nfout,*) ' !D meshl '
         do 11 ic = 1, nabc(3)
            do 12 ib = 1, nabc(2)
               write(nfout,9003) ib,ic
               write(nfout,9002) (meshl(ia,ib,ic),ia = 1, nabc(1))
 12         continue
 11      continue
 9003    format("ib,ic = ",2i4)
 9002    format(16i4)
#endif
      end if
      return
      end
C --------------------------------------------------------
C
C        packuc: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine packuc(catom,matom,natom)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      dimension catom(matom,3)
      do 1 k = 1, 3
         do 2 i = 1, natom
            n = catom(i,k)
            catom(i,k) = catom(i,k) - n
            if(catom(i,k).lt.0.0) catom(i,k) = catom(i,k) + 1.0
 2       continue
 1    continue
      return
      end
C --------------------------------------------------------
C
C        initid: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine initid(nfout,ioddst,m,n)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit none
      integer nfout,m,n, ioddst(m)
      integer i
      if(m.lt.n+1) then
         write(nfout,*) ' **** invalid array size (at sub. initid) ****'
         write(nfout,*) '         ( m = ',m,', n = ',n,')'
         call phase_error_with_msg(nfout,'invalid array size'//
     $    '(at sub. initidi)',__LINE__,__FILE__)
      endif
      do 1 i = 1, n
         ioddst(i) = i
 1    continue
      ioddst(m) = n
      return
      end
C --------------------------------------------------------
C
C        setowh: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine setowh(onwhca, n1, n2)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      integer ON, OFF
      parameter (ON = 1, OFF = 0)
      integer onwhca
      dimension onwhca(3,2)
      dimension nord(3)
      data nord/1,2,3/
      call initai(onwhca,6, OFF)
      onwhca(nord(n1),1) = ON
      onwhca(nord(n1),2) = ON
      onwhca(nord(n2),2) = ON
      return
      end
C --------------------------------------------------------
C
C        initai: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine initai(narray,m,nvalue)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      dimension narray(m)
      do 1 i = 1, m
         narray(i) = nvalue
 1    continue
      return
      end
C --------------------------------------------------------
C
C        mvnbr: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine mvnbr(catom,m,pmesh,ioddst,modst)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      dimension catom(m)
      dimension ioddst(modst)
C
      do 1 i = 1, ioddst(modst)
         ip = ioddst(i)
         d = catom(ip) - pmesh
         if(d.gt.0.5) then
            catom(ip) = catom(ip) - 1.0
         else if(d.lt.-0.5) then
            catom(ip) = catom(ip) + 1.0
         endif
 1    continue
      return
      end
C --------------------------------------------------------
C
C        sortsa: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine sortsa(nfout,catom,mcatm,abcmes,ttl,onwhca,wscr,
     &                  iodfrm, moddst,iodto,dstnc)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      integer nfout
      integer ON, OFF, onwhca
      parameter (ON = 1, OFF = 0)
      parameter (IAAX = 1, IBAX = 2, ICAX = 3 )
      dimension catom(mcatm,3), abcmes(3), ttl(6)
      dimension onwhca(3), iodfrm(moddst), iodto(moddst)
      dimension dstnc(mcatm)
      data icount/-20/
C
      icount = icount + 1
C
      katom = 0
      do 1 i = 1, iodfrm(moddst)
         jatom = iodfrm(i)
         if(onwhca(IAAX).eq.OFF) then
            x = 0.0
         else
            x = catom(jatom,IAAX) - abcmes(IAAX)
         endif
C
         if(onwhca(IBAX).eq.OFF) then
            y = 0.0
         else
            y = catom(jatom,IBAX) - abcmes(IBAX)
         endif
C
         if(onwhca(ICAX).eq.OFF) then
            z = 0.0
         else
            z = catom(jatom,ICAX) - abcmes(ICAX)
         endif
C
         xdst = ttl(1)*x*x + ttl(2)*y*y + ttl(3)*z*z
     &        + ttl(4)*x*y + ttl(5)*y*z + ttl(6)*z*x
         dstnc(jatom) = dsqrt(xdst)
         if(dstnc(jatom) .lt. wscr) then
            katom = katom + 1
            iodto(katom) = jatom
         endif
 1    continue
      iodto(moddst) = katom
C
c      if(icount.lt.4) then
      if(icount.lt.0) then
         write(nfout,*) ' *** iodto(before dsort) = ',
     &        iodto(moddst), ' ***'
         do 11 j = 1, iodto(moddst)
            jatom = iodto(j)
            write(nfout,9001) jatom, (catom(jatom,i), i = 1, 3),
     &           (abcmes(i),i=1,3),dstnc(jatom)
 11      continue
 9001    format(' !!ldos ','(',3f8.4,') (',3f8.4,')',f8.4)
      endif
C
      call bsort(iodto, moddst, dstnc, mcatm)
C
C     call dsort(nfout,iodto, moddst, dstnc, mcatm)
C
c      if(icount.lt.4) then
      if(icount.lt.0) then
         write(nfout,'(" !!ldos --- iodto(after dsort) = ",i6," ---")')
     &        iodto(moddst)
         do 21 j = 1, iodto(moddst)
            jatom = iodto(j)
            write(nfout,9001) jatom,(catom(jatom,i),i = 1, 3),
     &           (abcmes(i),i=1,3),dstnc(jatom)
 21      continue
      endif
C
      return
      end
C --------------------------------------------------------
C
C        sortsa: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
C$$ASASASASAS
C$$      subroutine sortsa_2(nfout,catom,mcatm,abcmes,ttl,onwhca,wscr
C$$     &     ,dface,dline,iodfrm, moddst,iodto,dstnc)
      subroutine sortsa_2(nfout,ipri,catom,mcatm,abcmes,ttl,onwhca,wscr
     &     ,dface,dline,iodfrm, moddst,iodto,dstnc)
C$$ASASASASAS
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8(a-h,o-z)
      integer ON, OFF, onwhca
      parameter (ON = 1, OFF = 0)
      parameter (IAAX = 1, IBAX = 2, ICAX = 3 )
      integer Case_D, Case_L, case
      parameter (Case_D = 1, Case_L = 2)
      dimension catom(mcatm,3), abcmes(3), ttl(6)
      dimension onwhca(3), iodfrm(moddst), iodto(moddst)
      dimension dstnc(mcatm), dface(3), dline(2,3)
      data icount/-20/
C
      icount = icount + 1
C
      if(onwhca(1)+onwhca(2)+onwhca(3).eq.ON) then
         case = Case_D
         if(onwhca(IAAX).eq.ON) then
            ipnt = IAAX
         else if(onwhca(IBAX).eq.ON) then
            ipnt = IBAX
         else
            ipnt = ICAX
         endif
      else if(onwhca(1)+onwhca(2)+onwhca(3).eq.2*ON) then
         case = Case_L
      else
         call phase_error_with_msg(nfout,' !E onwhca error'
     $                            ,__LINE__,__FILE__)
      endif
c
      katom = 0
      if(case.eq.Case_D) then
         do 1 i = 1, iodfrm(moddst)
         jatom = iodfrm(i)
         alpha = catom(jatom,ipnt) - abcmes(ipnt)
         xdst = alpha*dface(ipnt)
         dstnc(jatom) = dabs(xdst)
         if(dstnc(jatom) .lt. wscr) then
            katom = katom + 1
            iodto(katom) = jatom
         endif
 1       continue
      else if(case.eq.Case_L) then
         do 11 i = 1, iodfrm(moddst)
         jatom = iodfrm(i)
         alpha = catom(jatom,1) - abcmes(1)
         beta  = catom(jatom,2) - abcmes(2)
         gamma = catom(jatom,3) - abcmes(3)
         if(onwhca(1).eq.OFF) then
            dh = beta*dline(1,1) + gamma*dline(2,1)
            x = -dh
            y = beta
            z = gamma
         else if(onwhca(2).eq.OFF) then
            dh = alpha*dline(1,2) + gamma*dline(2,2)
            x = alpha
            y = -dh
            z = gamma
         else if(onwhca(3).eq.OFF) then
            dh = alpha*dline(1,3) + beta*dline(2,3)
            x = alpha
            y = beta
            z = -dh
         endif
         xdst = ttl(1)*x*x + ttl(2)*y*y + ttl(3)*z*z
     &        + ttl(4)*x*y + ttl(5)*y*z + ttl(6)*z*x
         dstnc(jatom) = dsqrt(xdst)
         if(dstnc(jatom) .lt. wscr) then
            katom = katom + 1
            iodto(katom) = jatom
         endif
 11      continue
      endif
      iodto(moddst) = katom
C
C$$ASASASASAS
C$$      if(icount.lt.0) then
      if(ipri >= 2 .and. icount.lt.0) then
C$$ASASASASAS
         write(nfout,*) ' *** iodto(before dsort) = ',
     &        iodto(moddst), ' ***'
         do 21 j = 1, iodto(moddst)
            jatom = iodto(j)
            write(nfout,9001) jatom, (catom(jatom,i), i = 1, 3),
     &           (abcmes(i),i=1,3),dstnc(jatom)
 21      continue
 9001    format(' !!ldos ',i4,'(',3f8.4,') (',3f8.4,')',f8.4)
      endif
C
      call bsort(iodto, moddst, dstnc, mcatm)
C
C     call dsort(nfout,iodto, moddst, dstnc, mcatm)
C
c      if(icount.lt.4) then
C$$ASASASASAS
C$$      if(icount.lt.0) then
      if(ipri >= 2 .and. icount.lt.0) then
C$$ASASASASAS
         write(nfout,*) ' --- iodto(after dsort) = ',
     &        iodto(moddst),' ---'
         do 31 j = 1, iodto(moddst)
            jatom = iodto(j)
            write(nfout,9001) jatom,(catom(jatom,i),i = 1, 3),
     &           (abcmes(i),i=1,3),dstnc(jatom)
 31      continue
      endif
C
      return
      end
C --------------------------------------------------------
C
C        bsort: coded by T.Yamasaki 13th Oct. 1992
C 
C --------------------------------------------------------
      subroutine bsort(iodto, moddst, dstnc, mcatm)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension iodto(moddst), dstnc(mcatm)
C
      nelemt = iodto(moddst)
C
      do 1 i = 1, nelemt-1
         pbas = dstnc(iodto(i))
         pmin = pbas
         do 2 j = i+1, nelemt
            p = dstnc(iodto(j))
            if(p.lt.pmin) then
               pmin = p
               jmin = j
            endif
 2       continue
         if(pmin.lt.pbas) then
            k = iodto(i)
            iodto(i) = iodto(jmin)
            iodto(jmin) = k
         endif
 1    continue
      return
      end
C --------------------------------------------------------
C
C        dsort: coded by T.Yamasaki 8th Oct. 1992
C 
C --------------------------------------------------------
      subroutine dsort(nfout,iodto, moddst, dstnc, mcatm)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer h, t
      parameter (mstack = 60)
      dimension iodto(moddst), dstnc(mcatm)
C
      dimension ISH(0:MSTACK), IST(0:MSTACK)
c
      mmax = mcatm
C
      nstack = 0
      is = 0
      h = 1
C
      nelemn = iodto(moddst)
      t = nelemn
 1001 continue
      if(h.ge.t) then
         if(is.ge.0) then
            if(is.gt.mstack) call exstak
            h = ish(is)
            t = ish(is)
         endif
         is = is - 1
      else
         np = (h+t)/2 + 0.5
         if(np.gt.mmax) write(nfout,*) ' !!ldos warning *** np > nmax'
         p = dstnc(iodto(np))
         i = h
         j = t
C
         if(np.ne.t) then
            call swapio(iodto,np,t)
         endif
 1000    continue
         if(j.lt.i) goto 2000
 1002    continue
         if(dstnc(iodto(j)).ge.p) then
            j = j - 1
            if(j.lt.h) goto 1000
            goto 1002
         endif
 1003    continue
         if(dstnc(iodto(i)).lt.p) then
            i = i + 1
            goto 1003
         endif
         if(i.lt.j) then
            call swapio(iodto,i,j)
            i = i + 1
            j = j - 1
         endif
         goto 1000
 2000    continue
C
         if(i.ne.t) then
            call swapio(iodto,i,t)
         endif
C
         if(i.lt.t) then
            is = is + 1
            if(is.gt.nstack) nstack = nstack + 1
            if(is.gt.mstack) call exstak
            ish(is) = i + 1
            ist(is) = t
         endif
         t = j
      endif
      if(is.ge.0) goto 1001
C
      return
      end
C --------------------------------------------------------
C
C        swapio: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine swapio(iodto,i,j)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension iodto(*)
      k = iodto(i)
      iodto(i) = iodto(j)
      iodto(j) = k
      return
      end
C --------------------------------------------------------
C
C        getmin: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine getmin(nfout,n3ax,catom,matom,abcmes,ttl,wscr,
     &                  ioddst,moddst,dstnc,iatom)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer nfout
      dimension catom(matom,3),abcmes(3),ttl(6)
      dimension ioddst(moddst),dstnc(matom)
      data delta/1.d-5/
      data icycle/0/
c
      if(ioddst(moddst).lt.1) then
         iatom = matom+1
         return
      endif
      call mvnbr(catom(1,n3ax),matom,abcmes(n3ax)
     &     ,ioddst,moddst)
#ifdef DEBUG
      write(nfout,*) ' !!ldos === ioddst(after dsort) = ',
     &     ioddst(moddst),' ==='
 9001 format(' !!ldos ',i4,'(',3f8.4,') (',3f8.4,')',f8.4)
      do 21 j = 1, ioddst(moddst)
      jatom = ioddst(j)
      write(nfout,9001) jatom,(catom(jatom,i),i = 1, 3),
     &     (abcmes(i),i=1,3),dstnc(jatom)
 21   continue
#endif
      call caldst(nfout,catom,matom,ioddst(1),abcmes,ttl,dst)
      dmin = dst
      iatom = ioddst(1)
      if(ioddst(moddst).ge.2) then
         do 1 i = 2, ioddst(moddst)
         ip = ioddst(i)
         if(dmin.lt.dstnc(ip)) then
            goto 1001
         endif
         call caldst(nfout,catom,matom,ip,abcmes,ttl,dst)
         if(dst.lt.dmin) then
            dmin = dst
            iatom = ip
         endif
C           if(dst.lt.dstnc(ip)) then
C             return
C           endif
 1       continue
 1001    continue
      endif
      if(dmin.gt.wscr) iatom = matom+1
      return
      end
C --------------------------------------------------------
C
C        getmin: coded by T.Yamasaki 10th Jun 1992
C 
C --------------------------------------------------------
      subroutine getmin_2(nfout,catom,matom,abcmes,ttl
     &     ,wscr,delta_dist,iatom)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer nfout
      dimension catom(matom,3),abcmes(3),ttl(6)
      real*8 dst

c$$$      real*8 rnd
      data icount/0/
c$$$      integer save icount = 0
c
      call caldst(nfout,catom,matom,1,abcmes,ttl,dst)
      dmin = dst

#ifdef DEBUG
      write(nfout,'(" !!ldos dist(",i4,") = ",f12.8)') 1, dst
#endif
      iatom = 1
      do 1 i = 2, matom
      ip = i
      call caldst(nfout,catom,matom,ip,abcmes,ttl,dst)

#ifdef DEBUG
      write(nfout,'(" !!ldos dist(",i4,") = ",f12.8)') i, dst
#endif
      if(dabs(dst-dmin) .lt. delta_dist) then
         icount = icount + 1
         if(mod(icount,2) == 0) then
c$$$         rnd = rand(0)
c$$$         write(nfout,'(" !!ldos icount = ",i8," rnd = ",f12.8)')
c$$$     &        icount, rnd
c$$$         if(rnd < 0.5d0) then
            dmin = dst
            iatom = ip
         end if
      else if(dst.lt.dmin) then
         dmin = dst
         iatom = ip
      endif
 1    continue
      if(dmin.gt.wscr) iatom = matom + 1
      return
      end

      subroutine getmin_3(nfout,catom,matom,abcmes,ttl
     &     ,wscr,delta_dist,iatom,dm)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer nfout
      dimension catom(matom,3),abcmes(3),ttl(6)
      real*8 dst

c$$$      real*8 rnd
      data icount/0/
c$$$      integer save icount = 0
c
      call caldst(nfout,catom,matom,1,abcmes,ttl,dst)
      dmin = dst

#ifdef DEBUG
      write(nfout,'(" !!ldos dist(",i4,") = ",f12.8)') 1, dst
#endif
      iatom = 1
      do 1 i = 2, matom
      ip = i
      call caldst(nfout,catom,matom,ip,abcmes,ttl,dst)

#ifdef DEBUG
      write(nfout,'(" !!ldos dist(",i4,") = ",f12.8)') i, dst
#endif
      if(dabs(dst-dmin) .lt. delta_dist) then
         icount = icount + 1
         if(mod(icount,2) == 0) then
c$$$         rnd = rand(0)
c$$$         write(nfout,'(" !!ldos icount = ",i8," rnd = ",f12.8)')
c$$$     &        icount, rnd
c$$$         if(rnd < 0.5d0) then
            dmin = dst
            iatom = ip
         end if
      else if(dst.lt.dmin) then
         dmin = dst
         iatom = ip
      endif
 1    continue
      if(dmin.gt.wscr) iatom = matom + 1
      dm = dmin
      return
      end

C --------------------------------------------------------
C
C        getmin_layer: coded by T.Yamasaki 7th Apr 2017
C 
C --------------------------------------------------------
      subroutine getmin_layer(nfout,winlay,mlayer,normal_axis
     &     ,abcmes,ilayer)
c
      implicit real*8 (a-h,o-z)
      integer nfout,mlayer
      dimension winlay(mlayer,2),abcmes(3)

      ilayer = 0
      dd = abcmes(normal_axis)
      do 1 i = 1, mlayer
         zmin = winlay(i,1)
         zmax = winlay(i,2)
         if(dd-zmin>0.5) zmin=zmin+1.0
         if(dd-zmin<-0.5) zmin=zmin-1.0
         if(dd-zmax>0.5) zmax=zmax+1.0
         if(dd-zmax<-0.5) zmax=zmax-1.0
         if(zmin<=dd .and. dd< zmax) then
            ilayer = i
            goto 2
         end if
 1    continue
 2    continue
      if(ilayer > mlayer) ilayer = 0

      return
      end
C --------------------------------------------------------
C
C        caldst: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine caldst(nfout,catom,matom,ip,abcmes,ttl,dst)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer nfout
      dimension catom(matom,3),ttl(6),abcmes(3)
      dimension xyz(3)
C      write(nfout,*) ' *** abcmes(1-3) = '
C     &      ,abcmes(1),abcmes(2),abcmes(3),' ***'
C      write(nfout,*) ' *** ip = ', ip
c$$$      do 1 i = 1, 3
c$$$         xyz(i) = catom(ip,i) - abcmes(i)
c$$$ 1    continue
c$$$      xdst =  ttl(1)*xyz(1)*xyz(1) 
c$$$     &     + ttl(2)*xyz(2)*xyz(2)
c$$$     &     + ttl(3)*xyz(3)*xyz(3)
c$$$     &     + ttl(4)*xyz(1)*xyz(2)
c$$$     &     + ttl(5)*xyz(2)*xyz(3)
c$$$     &     + ttl(6)*xyz(3)*xyz(1)
c$$$C
c$$$      dst = dsqrt(xdst)
c
      dst = 1.d9
      do 11 k = -1, 1
      xyz(3) = catom(ip,3) - abcmes(3) + k
      do 12 j = -1, 1
      xyz(2) = catom(ip,2) - abcmes(2) + j
      do 13 i = -1, 1
      xyz(1) = catom(ip,1) - abcmes(1) + i
c$$$      if(i*j*k.ne.0) then
#ifdef _ORTHOGONAL_
      xdst = ttl(1)*xyz(1)**2
     &     + ttl(2)*xyz(2)**2
     &     + ttl(3)*xyz(3)**2
#else
      xdst = ttl(1)*xyz(1)**2
     &     + ttl(2)*xyz(2)**2
     &     + ttl(3)*xyz(3)**2
     &     + ttl(4)*xyz(1)*xyz(2)
     &     + ttl(5)*xyz(2)*xyz(3)
     &     + ttl(6)*xyz(3)*xyz(1)
#endif
      xdst = dsqrt(xdst)
      if(xdst.lt.dst) dst = xdst
 13   continue
 12   continue
 11   continue
c      
      return
      end
C --------------------------------------------------------
C
C        calttl: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine calttl(ipri,brmt,ttl)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension brmt(3,3),ttl(6)
      do 11 j = 1, 6
         ttl(j) = 0.d0
 11   continue
      do 1 j = 1, 3
         jj = mod(j+1,3)
         if(jj.eq.0) jj = 3
         do 2 i = 1, 3
            ttl(j) = ttl(j) + brmt(i,j)*brmt(i,j)
            ttl(j+3) = ttl(j+3) + 2.d0*brmt(i,j)*brmt(i,jj)
 2       continue
 1    continue
      if(ipri>=1) write(6,9001) (ttl(j),j = 1, 6)
 9001 format(' !!ldos (ttl(1-6)) = (',6f12.4,')')
C
      return
      end
C --------------------------------------------------------
C
C        distbs: coded by T.Yamasaki 21st Apr 1996
C 
C --------------------------------------------------------
      subroutine distbs(nfout,ipri,brmt, rcmt, dface, dline)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      integer nfout, ipri
      dimension brmt(3,3),rcmt(3,3), dface(3), dline(2,3)
c
c   checking of orthogonality
      if(ipri >= 2) then
         write(nfout,'(" !!ldos: -- checking of orthogonality --")')
         do i = 1, 3
            do j = 1, 3
               prod = 0.d0
               do 3 k = 1, 3
                  prod = prod + brmt(k,i)*rcmt(k,j)
 3             continue
               write(nfout,'(" !!ldos - product of brmt and rcmt -")')
               write(nfout,'(" !!ldos i = ",i5," j = ",i5
     &                                 ," prod = ",f20.12)') i, j, prod
            end do
         end do
      end if
c
      do 11 i = 1, 3
         p = 0.d0
         d = 0.d0
         do 12 j = 1, 3
            p = p + brmt(j,i)*rcmt(j,i)
            d = d + rcmt(j,i)*rcmt(j,i)
 12      continue
         dface(i) = p/dsqrt(d)
 11   continue
c
      if(ipri >= 2) then
         write(nfout,9001) (dface(j),j = 1, 3)
 9001    format(' !!ldos (dface(1-3)) = (',3f12.4,')')
      end if
c
      do 21 i = 1, 3
         d = 0.d0
         do 32 k = 1, 3
            d = d + brmt(k,i)*brmt(k,i)
 32      continue
         do 22 j = 1, 2
            jj = i+j
            jj = mod(jj,3)
            if(jj .eq. 0) jj = 3
            p = 0.d0
            do 23 k = 1, 3
               p = p + brmt(k,jj)*brmt(k,i)
 23         continue
            if(i.eq.2) then
               dline(3-j,i) = p/d
            else
               dline(j,i) = p/d
            endif
 22      continue
 21   continue

      if(ipri >= 1) then
         write(nfout,9002) ((dline(j,i),j=1,2),i=1,3)
 9002    format(' !!ldos (dline(1-2,1-3)) = ,',3('(',2f12.4,')'))
      end if
C
      return
      end
C --------------------------------------------------------
C
C        patmbr: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine patmbr(catxyz,matom,natom,rcmt,catom)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      parameter (LBRAVS = 1, LCARTS = 2, MCORDS = 2)
      dimension catxyz(matom,3),catom(matom,3)
      dimension rcmt(3,3),work(3)
C
      pi2 = 8.d0*datan(1.d0)
C
      do 1 i = 1, natom
         do 2 j = 1, 3
            work(j) = 0.d0
            do 3 k = 1, 3
               work(j) = work(j) + rcmt(k,j)*catxyz(i,k)/pi2
 3          continue
 2       continue
         catom(i,1) = work(1)
         catom(i,2) = work(2)
         catom(i,3) = work(3)
C
 1    continue
C
      return
      end
C --------------------------------------------------------
C
C        gtrcmt: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine gtrcmt(brmt,rcmt)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension brmt(3,3),rcmt(3,3)
C
      do 1 i = 1, 3
         j = mod(i,3) + 1
         k = mod(i+1,3) + 1
         do 2 ii = 1, 3
            jj = mod(ii,3) + 1
            kk = mod(ii+1,3) + 1
            rcmt(ii,i) = brmt(jj,j)*brmt(kk,k)
     &                 - brmt(kk,j)*brmt(jj,k)
 2       continue
 1    continue
C
      volunt = 0.d0
      do 11 i = 1, 3
         volunt = volunt + brmt(i,1)*rcmt(i,1)
 11   continue
C
      do 21 i = 1, 3
         do 22 j = 1, 3
            rcmt(j,i) = rcmt(j,i)/volunt
 22      continue
 21   continue
C
      write(6,9212) ((rcmt(j,i),i=1,3),j=1,3)
 9212 format(/' !!ldos rcmt-matrix',3(/3f10.5))
      return
      end
C --------------------------------------------------------
C
C        gtbrmt: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine gtbrmt(brmt)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension brmt(3,3)
      do 1 j = 1, 3
         do 2 i = 1, 3
            brmt(i,j) = 0.d0
 2       continue
 1    continue
      brmt(1,1) = 14.49568901
      brmt(2,2) = 14.49568901
      brmt(3,3) = 51.25
C
      return
      end
C --------------------------------------------------------
C
C        posatm: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine posatm(catxyz,matom,natom,rcmt,catom)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension catxyz(matom,3), catom(matom,3)
      dimension rcmt(3,3)
C
      call patmbr(catxyz,matom,natom,rcmt,catom)
C
c      write(6,9001) ((catxyz(ia,j),j = 1, 3),ia = 1, natom)
c 9001 format(' *** catxyz ***',99(/3f12.4))
C
c      write(6,9002) ((catom(ia,j),j = 1, 3),ia = 1, natom)
c 9002 format(' *** catoms(in BRAVaiS) ***',99(/3f12.4))

      return
      end
C --------------------------------------------------------
C
C        gtnabc: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
      subroutine gtnabc(nabc,mesha,meshb,meshc)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension nabc(3)
C
      nabc(1) = mesha
      nabc(2) = meshb
      nabc(3) = meshc
C
      return
      end
C --------------------------------------------------------
C
C        gtwscr: coded by T.Yamasaki 8th Oct 1992
C 
C --------------------------------------------------------
c$$$      subroutine gtwscr(nfout,brmt,wscr)
c$$$c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
c$$$      implicit real*8 (a-h,o-z)
c$$$      integer nfout
c$$$      dimension brmt(3,3)
c$$$C
c$$$c$$$      z = 0.2
c$$$c$$$      wscr = z * brmt(1,1)
c$$$c$$$      if(wscr.lt.5.0) then
c$$$c$$$         wscr = 5.d0
c$$$c$$$      endif
c$$$      wscr = 6.0d0
c$$$      write(nfout,*) ' !!ldos **** polyhedral critical radius = ', wscr
c$$$c        unit is a.u.
c$$$c
c$$$      return
c$$$      end
C -------------
      subroutine wdmesp(mesh,ifxr,ifyr,ifzr,nfile)
c                           @(#)anlmes.f 9.2 01/05/14 20:51:42
      implicit real*8 (a-h,o-z)
      dimension mesh(ifxr,ifyr,ifzr)
C
      mt = ifxr*ifyr*ifzr
C
      write(nfile,*) ' CHARGE DENSITY NE = ', mt,' ( ',
     &     ifyr,ifzr,ifxr, ' )'
C
      write(nfile,9001)
     &    (((mesh(i,j,k), j = 1, ifyr), k = 1, ifzr),i = 1, ifxr)
 9001 format(20i4)
C
      return
      end
C ----------------------------------------------
C       Nov. 17th 1992 by T.Yamasaki
C
C       modified by T. Yamasaki, 09th Feb 2004
C ----------------------------------------------
C                           %Z%%M% %I% %E% %U%
      subroutine dwinly2(kimg,catxyz,matom,natom,winlay,mlayer,nlayer
     &                 ,numlay,nhatom,maxhv,slicing_selection,deltaz
     &                 ,crtdst_winlay,nfout)
      implicit none
      integer kimg,matom,natom,mlayer,nlayer,nhatom,slicing_selection

c$$$      parameter ( crtdst = 3.5 )
      real*8 catxyz(matom,3),winlay(mlayer,2)
      integer numlay(nhatom)
      real*8 maxhv
      real*8 deltaz, crtdst_winlay
      integer nfout

      real*8 tmp, maxv, zmin
      integer i, nl

      integer REGULAR_INTERVALS, BY_ATOMIC_POSITION
      parameter ( REGULAR_INTERVALS = 1, BY_ATOMIC_POSITION = 2)

      if(kimg .eq.1) maxv = maxhv*2
      if(kimg .eq.2) maxv = maxhv

      call rsreal(mlayer*2, winlay)
      if(slicing_selection .eq. REGULAR_INTERVALS) then
         write(nfout,'(" !!ldos slicing_selection = EAZY")')
         if(kimg .eq. 1) then
            winlay(1,1) = 0.d0
            winlay(mlayer,2) = maxhv
            do i = 1, nlayer
              winlay(i,2) = deltaz*i
              if(winlay(i,2) .gt. maxhv) winlay(i,2) = maxhv
              if(i+1 .le. mlayer) winlay(i+1,1) = winlay(i,2)
            enddo
         else if(kimg .eq. 2) then
            zmin = maxv*30
            do i = 1, nhatom
              if(zmin .gt. catxyz(i,1)) zmin = catxyz(i,1)
            enddo
            winlay(1,1) = zmin - crtdst_winlay
            winlay(mlayer,2) = winlay(1,1) + maxv
            do i = 1, nlayer
              winlay(i,2) = winlay(i,1) + deltaz
              if(winlay(i,2) .gt. winlay(mlayer,2))
     &             winlay(i,2) = winlay(mlayer,2)
              if(i+1. le. mlayer) winlay(i+1,1) = winlay(i,2)
            enddo
         endif
      else if(slicing_selection .eq. BY_ATOMIC_POSITION) then
         write(nfout,'(" !!ldos slicing_selection = BY_ATOMS")')
            
         do i = 1, mlayer
           winlay(i,1) = +maxv*30
           winlay(i,2) = -maxv*30
         enddo

         do i = 1, nhatom
           nl = numlay(i)
           if(winlay(nl,1) .gt. catxyz(i,1)) winlay(nl,1) = catxyz(i,1)
           if(winlay(nl,2) .lt. catxyz(i,1)) winlay(nl,2) = catxyz(i,1)
         enddo
c           winlay(nl,1) : minimum z among nl-th layer atomic coordinates
c           winlay(nl,2) : maxmum  z among nl-th layer atomic coordinates
c
         write(nfout,*)
     &        ' !!ldos a range of atomic positions of each layer ****'
         do i = 1, nlayer
           write(nfout,9001) i, winlay(i,1), winlay(i,2)
         enddo
C
         tmp = winlay(1,1) - crtdst_winlay
         if(kimg .eq. 2) then
            if(tmp .lt. winlay(nlayer,2) - maxhv ) then
c     ( a unit cell has no vacuum region )
               winlay(1,1) = (winlay(1,1) + winlay(nlayer,2)-maxv)*0.5
               winlay(nlayer,2) = winlay(1,1)+maxv
               winlay(mlayer,1) = winlay(1,1)
               winlay(mlayer,2) = winlay(1,1)
            else
c     ( a unit cell has a vacuum region)
               winlay(1,1) = winlay(1,1) - crtdst_winlay
               winlay(nlayer,2) = winlay(nlayer,2) + crtdst_winlay
               winlay(mlayer,1) = winlay(nlayer,2)
               winlay(mlayer,2) = winlay(1,1) + maxv
            endif
         else if(kimg .eq. 1) then
            winlay(1,1)      = 0
            winlay(mlayer,2) = maxhv
            winlay(nlayer,2) = winlay(nlayer,2) + crtdst_winlay
            winlay(mlayer,1) = winlay(nlayer,2)
         endif

         do i = 1, nlayer-1
           winlay(i,2) = (winlay(i,2) + winlay(i+1,1))*0.5
           winlay(i+1,1) = winlay(i,2)
         enddo
      endif

      write(nfout,*) ' !!ldos     no,        min,           max '
      do i = 1, nlayer
      write(nfout,9001) i, winlay(i,1), winlay(i,2)
      enddo
      if(mlayer.gt.nlayer) then
         write(nfout,9001) mlayer, winlay(mlayer,1), winlay(mlayer,2)
      endif
 9001 format(' !!ldos ',i4,2f20.8)
C
      return
      end
c -- * -- 1 -- * -- 2 -- * -- 3 -- * -- 4 -- * -- 5 -- * -- 6 -- * -- 7
c                           @(#)exstak.f 9.1 97/09/10 20:34:10
c
c ----- subroutine exstak ---------------------------------------------
c
c                        1994 6/26    by yamasaki takahiro
c
      subroutine exstak
c
      write(6,*) ' !!ldos == error   at exstak == '
      write(6,*) ' !!ldos   shortage of stack size '
      write(6,*) ' !!ldos   expande stack size '
      call phase_error_with_msg(nfout,'error   at exstak'
     $                         ,__LINE__,__FILE__)
      end
