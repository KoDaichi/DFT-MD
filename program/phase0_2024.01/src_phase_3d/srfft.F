!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: hitachi_z3fft, set_length
!
!  AUTHOR(S): O. Sugino  December/22/2003
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
#ifdef SRFFT
      subroutine hitachi_z3fft(afft,ftw, id,nl,nm,nn, ip1,ip2,ip3,
     & cw1,cw2,cw3, inverse_or_direct, itable_flag)
! $Id: srfft.F 570 2017-04-21 20:34:50Z yamasaki $
      implicit none
      complex*16 afft(*), ftw(*)
      integer id, nl, nm, nn
      integer ip1(3), ip2(3), ip3(3)
      complex*16 cw1(*), cw2(*), cw3(*)
      integer inverse_or_direct, itable_flag
      integer ier, nwkw
      integer ndim1, ndim2
#ifndef SRFTOUTSMP
      integer iopt(2)
      integer i
#else
      integer iopt(3)
      integer i,itop,nloops
#endif

      ndim1 = id
      ndim2 = nm
      iopt(1) = inverse_or_direct
      iopt(2) = itable_flag
      nwkw = 1
#ifndef SRFTOUTSMP
      call hzft7m(afft,ndim1,ndim2,nl,nm,nn,ip1,ip2,ip3,
     & iopt,cw1,cw2,cw3,ftw,nwkw,ier)
#else
      if(itable_flag.eq.0)then
        nloops=1
      else
        nloops=nn
      endif

*poption parallel
      do i=1,nloops
        itop=(ndim1*ndim2)*(i-1)+1
        call hzft9m(afft(itop),ndim1,ndim2,nl,nm,ip1,ip2,
     &         iopt,cw1,cw2,ftw(itop),nwkw,ier)
      enddo

      call trans_axis(ftw,afft,nl,nm,nn,nwkw)

      iopt(3)=2
      if(itable_flag.eq.0)then
        nloops=1
      else
        nloops=nm
      endif

*poption parallel
      do i=1,nloops
        itop=((nl+nwkw)*nn)*(i-1)+1
        call hzft5m(ftw(itop),nl+nwkw,nn,ip3,nl,
     &         iopt,cw3,afft(itop),ier)
      enddo

      call trans_axis(afft,ftw,nl,nn,nm,nwkw)
#endif

!ccc This loop is redundant.
      do i = 1, nm*nn
          afft(i*id) = (0.Q0,0.Q0)
      end do

      return
      end

      subroutine set_length(n, ip)
      implicit none
      integer n, ip(3)
      integer nn, lp2, lp3, lp5

      nn =n
      lp5=0
      lp3=0
      lp2=0
      do while(mod(nn,5).eq.0)
          lp5=lp5+1
          nn=nn/5
      end do
      do while(mod(nn,3).eq.0)
          lp3=lp3+1
          nn=nn/3
      end do
      do while(mod(nn,2).eq.0)
          lp2=lp2+1
          nn=nn/2
      end do
      ip(1) = lp2
      ip(2) = lp3
      ip(3) = lp5

      if(nn.ne.1) stop 'n is an improper number.'

      return
      end
#ifdef SRFTOUTSMP
      subroutine trans_axis(c1,c2,n1,n2,n3,nwkw)
      complex*16 c1(n1+nwkw,n3,n2),c2(n1+nwkw,n2,n3)

      do k=1,n3
      do j=1,n2
      do i=1,n1
        c1(i,k,j)=c2(i,j,k) 
      enddo
      enddo
      enddo

      return
      end
#endif
#endif

