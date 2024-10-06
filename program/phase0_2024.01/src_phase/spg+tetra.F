      module m_spg_tetra
        implicit none

!       for nskpb0
        integer :: np1s,np2s
        real*8,  allocatable, dimension(:,:) :: pb0s
        integer, allocatable, dimension(:)   :: ip10s,ip20s
        integer, allocatable, dimension(:)   :: ip01s,ip21s
        integer, allocatable, dimension(:)   :: ip02s,ip12s
        integer, allocatable, dimension(:)   :: iu21s,iv21s
        integer :: lmnp0s=0, lmnp1s=0, lmnp2s=0

!       for wtetra
        integer, allocatable, dimension(:)   :: ip2cubs
        integer :: ncubs=0

        contains

        logical function wtetra_vars_ready(ncub)
          integer, intent(in) :: ncub
          wtetra_vars_ready = ncub .eq. ncubs .and.
     &                        allocated(ip2cubs)
          return
        end function wtetra_vars_ready

        logical function nspb0_vars_ready(lmnp0, lmnp1, lmnp2)
          integer, intent(in) :: lmnp0, lmnp1, lmnp2
          nspb0_vars_ready = lmnp0s .eq. lmnp0 .and.
     &                       lmnp1s .eq. lmnp1 .and.
     &                       lmnp2s .eq. lmnp2 .and.
     &                       allocated(pb0s)
          return
        end function nspb0_vars_ready

        subroutine get_wtetra_vars(ip2cub)
          integer, dimension(:), intent(out) :: ip2cub
          ip2cub = ip2cubs
          return
        end subroutine get_wtetra_vars

        subroutine set_wtetra_vars(ncub,ip2cub)
          integer, intent(in) :: ncub
          integer, dimension(ncub), intent(in) :: ip2cub
          if(allocated(ip2cubs)) deallocate(ip2cubs)
          ncubs = ncub
          allocate(ip2cubs(ncubs))
          ip2cubs = ip2cub
          return
        end subroutine set_wtetra_vars

        subroutine get_nskpb0_vars(np1, np2, pb0,
     &             ip10, ip20, ip01,ip21,ip02,ip12,iu21,iv21)
          integer, intent(out) :: np1,np2
          real*8,  dimension(:,:), intent(out) :: pb0
          integer, dimension(:),   intent(out) :: ip10, ip20
          integer, dimension(:),   intent(out) :: ip01,ip21
          integer, dimension(:),   intent(out) :: ip02,ip12
          integer, dimension(:),   intent(out) :: iu21,iv21
          np1 = np1s
          np2 = np2s
          pb0 = pb0s
          ip10=ip10s
          ip20=ip20s
          ip01=ip01s
          ip21=ip21s
          ip02=ip02s
          ip12=ip12s
          iu21=iu21s
          iv21=iv21s
          return
        end subroutine get_nskpb0_vars

        subroutine set_nskpb0_vars(lmnp0,lmnp1,lmnp2,np1, np2, pb0,
     &             ip10, ip20, ip01,ip21,ip02,ip12,iu21,iv21)
          integer, intent(in) :: lmnp0, lmnp1, lmnp2
          integer, intent(in) :: np1,np2
          real*8,  dimension(3,lmnp0), intent(in) :: pb0
          integer, dimension(lmnp0), intent(in) :: ip10, ip20
          integer, dimension(lmnp1), intent(in) :: ip01,ip21
          integer, dimension(lmnp2), intent(in) :: ip02,ip12
          integer, dimension(lmnp1), intent(in) :: iu21,iv21
          np1s = np1
          np2s = np2
          if(allocated(pb0s))  deallocate(pb0s)
          if(allocated(ip10s)) deallocate(ip10s)
          if(allocated(ip20s)) deallocate(ip20s)
          if(allocated(ip01s)) deallocate(ip01s)
          if(allocated(ip21s)) deallocate(ip21s)
          if(allocated(ip02s)) deallocate(ip02s)
          if(allocated(ip12s)) deallocate(ip12s)
          if(allocated(iu21s)) deallocate(iu21s)
          if(allocated(iv21s)) deallocate(iv21s)

          allocate(pb0s(3,lmnp0));pb0s=pb0
          allocate(ip10s(lmnp0));ip10s=ip10
          allocate(ip20s(lmnp0));ip20s=ip20
          allocate(ip01s(lmnp1));ip01s=ip01
          allocate(ip21s(lmnp1));ip21s=ip21
          allocate(ip02s(lmnp2));ip02s=ip02
          allocate(ip12s(lmnp2));ip12s=ip12
          allocate(iu21s(lmnp1));iu21s=iu21
          allocate(iv21s(lmnp1));iv21s=iv21
          lmnp0s = lmnp0
          lmnp1s = lmnp1
          lmnp2s = lmnp2
          return
        end subroutine set_nskpb0_vars

      end module m_spg_tetra

      subroutine wd_spg_tetra_vars()
        use m_spg_tetra
        use m_Control_Parameters, only : sw_output_kpoint_info
        use m_Const_Parameters, only : OFF
        use m_Files, only : m_Files_open_kpoint_bin,
     &      m_Files_close_kpoint_bin, nfkpoint_bin, F_KPOINT_BIN
        use m_Parallelization, only : mype, MPI_CommGroup

        implicit none
        if(sw_output_kpoint_info==OFF) return
        if(allocated(pb0s)) then
          call m_Files_open_kpoint_bin()
          if(mype .eq. 0) then
            write(nfkpoint_bin) lmnp0s,lmnp1s, lmnp2s
            write(nfkpoint_bin) np1s,np2s,pb0s
            write(nfkpoint_bin) ip10s,ip20s,ip01s,ip21s,ip02s
     &           ,ip12s,iu21s,iv21s
          endif
        endif
        if(allocated(ip2cubs)) then
          if(mype .eq. 0) then
            write(nfkpoint_bin) ncubs
            write(nfkpoint_bin) ip2cubs
          endif
        endif

        call m_Files_close_kpoint_bin()
      end subroutine wd_spg_tetra_vars

      subroutine rd_spg_tetra_vars()
        use m_spg_tetra
        use m_Control_Parameters, only : sw_output_kpoint_info
        use m_Const_Parameters, only : OFF
        use m_Parallelization, only : mype, MPI_CommGroup
        use m_Files, only : m_Files_open_kpoint_bin,
     &      m_Files_close_kpoint_bin, nfkpoint_bin, F_KPOINT_BIN

        implicit none

        include 'mpif.h'
        logical :: exi, ierr

        if(sw_output_kpoint_info==OFF) return

        if(mype .eq. 0) inquire(file=F_KPOINT_BIN, exist = exi)
        call mpi_bcast(exi,1,mpi_logical,0, MPI_CommGroup,ierr)
        if(.not. exi) return

         call m_Files_open_kpoint_bin()
         if(mype .eq. 0) read(nfkpoint_bin,err=5) lmnp0s,lmnp1s,lmnp2s
         call mpi_bcast(lmnp0s,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(lmnp1s,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(lmnp2s,1,mpi_integer,0,MPI_CommGroup,ierr)
         if(allocated(pb0s))  deallocate(pb0s)
         if(allocated(ip10s)) deallocate(ip10s)
         if(allocated(ip20s)) deallocate(ip20s)
         if(allocated(ip01s)) deallocate(ip01s)
         if(allocated(ip21s)) deallocate(ip21s)
         if(allocated(ip02s)) deallocate(ip02s)
         if(allocated(ip12s)) deallocate(ip12s)
         if(allocated(iu21s)) deallocate(iu21s)
         if(allocated(iv21s)) deallocate(iv21s)

         allocate(pb0s(3,lmnp0s))
         allocate(ip10s(lmnp0s))
         allocate(ip20s(lmnp0s))
         allocate(ip01s(lmnp1s))
         allocate(ip21s(lmnp1s))
         allocate(ip02s(lmnp2s))
         allocate(ip12s(lmnp2s))
         allocate(iu21s(lmnp1s))
         allocate(iv21s(lmnp1s))

         if(mype .eq. 0) then
           read(nfkpoint_bin,err=5) np1s,np2s,pb0s
           read(nfkpoint_bin,err=5) ip10s,ip20s,ip01s,ip21s,ip02s
     &     ,ip12s,iu21s,iv21s
         endif

         call mpi_bcast(np1s,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(np2s,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(pb0s,3*lmnp0s,mpi_double_precision,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip10s,lmnp0s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip20s,lmnp0s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip01s,lmnp1s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip21s,lmnp1s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip02s,lmnp2s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(ip12s,lmnp2s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(iu21s,lmnp1s,mpi_integer,0,
     &        MPI_CommGroup,ierr)
         call mpi_bcast(iv21s,lmnp1s,mpi_integer,0,
     &        MPI_CommGroup,ierr)

5        continue

         if(mype .eq. 0) read(nfkpoint_bin,err=6) ncubs
         call mpi_bcast(ncubs,1,mpi_integer,0,MPI_CommGroup,ierr)
         if(allocated(ip2cubs)) deallocate(ip2cubs)
         allocate(ip2cubs(ncubs))
         if(mype .eq. 0) read(nfkpoint_bin,err=6) ip2cubs
         call mpi_bcast(ip2cubs,ncubs,mpi_integer,0,MPI_CommGroup,ierr)
6        continue
         call m_Files_close_kpoint_bin()

         return
      end subroutine rd_spg_tetra_vars

!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINES: fermi1, fermi2, fermi4, fermi5, lsatpz, lsdos0, lsnbat
!     lsnbst, lstt1i, msatpz, msdos0, msdosi, msnbat, msnbst, mstt1i,
!     msttie, msttii, msttio, nrrpc2, nrrpj0, nrrpl1, nrrpl2, nrskc0,
!     nrskc1, nrtrc0, nrtrc2, nrylm1, nrylm2, nsatpz, nsdos0, nsdos2,
!     nsdos3, nsdosi, nseulc, nseulh, nsgcm2, nsgcm3, nsgrp1, nsgrpb
!     nsgrpw, nsig10, nsjonh, nskg00, nskg01, nskgs0, nskgs1, nskma0,
!     nskp00, nskpb0, nskpbm, nskpr0, nskpw0, nslat3, nslata, nslatb,
!     nslatc, nslatr, nslatz, nsm3iw, nsmetr, nsmlt1, nsmltw, nsmult,
!     nsnbat, nsnbrd, nsnbst, nspace, nspbg0, nspbge, nspgr1, nspgrp,
!     nspgw1, nsprmv, nsrduc, nsrmxc, nsrmxh, nsrot1, nsrota, nsrotc,
!     nsroth, nsrotk, nsrotl, nsrotr, nss0, nss00r, nss01r, nss02r,
!     nssdjg, nssum1, nstrsh, nstt0c, nstt0i, nstt0r, nstt1i, nstt1r,
!     nstt2i, nstt3i, nsttod, nstts1, nsxlog, osatpz, osnbat, osnbst,
!     rdprp, setkp0, setkp0_n, ka00, setkp0_default, setkp0_default_n,
!     setspg, setspg_n, setspg_default, setspg_default_n, getspgtab,
!     tbspg, tspaca, tspacb, tsrmhi, wtetra, zzzy51, zzzy52, zzzy53
!
!  AUTHORS: N. Hamada, H. Mizouchi, K. Mae    August/20/2003
!           A. Yanase, and K. Terakura   June/07/1986
!
!  FURTHER MODIFICATION: T. Yamasaki     September/16/2007
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
c sub.fermi1       2====*====3====*====4====*====5====*====6====*====7
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c

! ========================= modified by K. Tagami ================== 11.0
!      subroutine fermi1(nxx,nyy,nzz,np2,lmnp2,neig,lmneig,nspin,
!     &                  eig2,ip20,np0,valenc,eferm,eband,valud,
!     &                  iwt,ip2cub,ipri )
!
      subroutine fermi1(nxx,nyy,nzz,np2,lmnp2,neig,lmneig,nspin,
     &                  eig2,ip20,np0,valenc,eferm,eband,valud,
     &                  iwt,ip2cub,ipri, flag_noncol )
! =================================================================== 11.0

c     $Id: spg+tetra.F 633 2020-12-01 05:11:03Z jkoga $
c
c      eferm :  fermi energy
c      dos   :  density of states at fermi energy for each spin
c      sidos :  integrated dos at fermi energy for each spin
c      ddd = sum of dos
c      sss = sum of sidos
c
      implicit real*8(a-h,o-z)
      integer ip20(np0)
      integer iwt(np2),ip2cub(nxx*nyy*nzz)
      integer ipri
      real*8  eig2(lmnp2,lmneig,nspin)
      !!$real*8  dos(2),sidos(2),eawk(20000),valud(nspin)
      real*8  dos(2),sidos(2),eawk(np0),valud(nspin)
      real*8  efm(200),qfm(200)

c ================================== added by K. Tagami ================- 11.0
      real*8 valenc
! === Optional argument needs interface!!! by T.Kato 2013/07/02 ================
!     logical, optional:: flag_noncol
      logical :: flag_noncol
! ==============================================================================
c ====================================================================== 11.0
c
#ifdef __TIMER_SUB__
      call timer_sta(710)
#endif

      idim=3
*     eps=dfloat(10)**(-6)
      eps=dfloat(10)**(-10)
      neig1=1
      neig2=neig
      !!$if(np0.gt.20000) then
      !!$  write(6,*) ' np0=',np0,' > 20000 : error at sub.fermi1.'
      !!$  stop 'error at sub.fermi1 (np0).'
      !!$end if
      emin=eig2(1,1,1)
      emax=eig2(1,1,1)
#ifdef __TIMER_DO__
      call timer_sta(811)
#endif
      do 10 is=1,nspin
      do 10 i=1,neig
      do 10 ip2=1,np2
        if(emin.gt.eig2(ip2,i,is)) emin=eig2(ip2,i,is)
        if(emax.lt.eig2(ip2,i,is)) emax=eig2(ip2,i,is)
   10 continue
#ifdef __TIMER_DO__
      call timer_end(811)
#endif
      if(ipri >= 2) write(6,*) ' emax=',emax,'   emin=',emin
      wei=2.d0
      if(nspin.eq.2) wei=1.d0

! ===========================- added by K. Tagami ==================== 11.0
      if ( flag_noncol ) wei = 1.0d0
! ==================================================================== 11.0

c modified by H.Sawada on May 1, 1997
      efermi = eferm
c modified by H.Sawada on May 1, 1997
      e1     = emin
      e2     = emax
c
      iii=0
      instts1 = 0
 55   continue
 33   continue
        ddd=0.d0
        sss=0.d0
        ebnd=0.d0
#ifdef __TIMER_DO__
      call timer_sta(812)
#endif
        do 32 ispin=1,nspin
        call nsdos2(idim,efermi,nxx,nyy,nzz,lmnp2,neig1,neig2,
     &              eig2(1,1,ispin),
     &              ip20,np0,eawk,dos(ispin),sidos(ispin),dm1,
     &              instts1,np2,iwt,ip2cub)
        ddd=ddd+dos(ispin)*wei
        sss=sss+sidos(ispin)*wei
        ebnd=ebnd+dm1*wei
        valud(ispin)=sidos(ispin)*wei
 32     continue
#ifdef __TIMER_DO__
      call timer_end(812)
#endif
        iii=iii+1
        if(iii.eq.100 .and. instts1.eq.0) then
           instts1 = 1
           efermi = eferm
           emin=eig2(1,1,1)
           emax=eig2(1,1,1)
#ifdef __TIMER_DO__
      call timer_sta(813)
#endif
           do 11 is=1,nspin
           do 11 i=1,neig
           do 11 ip2=1,np2
           if(emin.gt.eig2(ip2,i,is)) emin=eig2(ip2,i,is)
           if(emax.lt.eig2(ip2,i,is)) emax=eig2(ip2,i,is)
 11        continue
#ifdef __TIMER_DO__
      call timer_end(813)
#endif
           e1     = emin
           e2     = emax
           iii = 0
           goto 55
        elseif(iii.gt.100) then
           write(6,*) 'iii=',iii,'   efermi=',efermi
           write(6,*) 'sss=',sss,'  val=',valenc
           do 40 i=1,100
 40        write(6,*) ' efermi=',efm(i),'   q=',qfm(i)
           stop ' === stop sub.fermi1. (iii>200) ==='
        end if
        efm(iii)=efermi
        qfm(iii)=sss
        if(abs( sss-valenc        ).lt.eps) go to 34
c       if(abs((sss-valenc)/valenc).lt.eps) go to 34
c
        if(sss.lt.valenc) then
          e1     = efermi
          efermi = efermi + (e2-efermi)/2
        else
          e2     = efermi
          efermi = efermi + (e1-efermi)/2
        end if
      go to 33
c
   34 continue
c
      if(abs(ddd).lt.1.d-20) then
        call fermi2(np2,lmnp2,lmneig,neig,nspin,eig2,efermi,ipri)
      end if
c      write( 6,111) efermi,ddd,sss
c      write(16,111) efermi,ddd,sss
c  111 format(5x,'fermi-en=',f12.6,5x,'dos=',f14.6,5x,'int.dos=',f14.6)
c --*
      if(nspin.eq.2) then
        smom=sidos(1)-sidos(2)
c        write( 6,120) smom
c        write(16,120) smom
      end if
c  120 format(' magnetic moment=',f12.6,' per cell')
      eferm=efermi
      eband=ebnd
#ifdef __TIMER_SUB__
      call timer_end(710)
#endif
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine fermi2(np2,lmnp2,lmneig,neig,nspin,eig2,eferm,ipri)
c
      implicit real*8(a-h,o-z)
      dimension eig2(lmnp2,lmneig,nspin)
      integer ipri
c
#ifdef __TIMER_SUB__
      call timer_sta(712)
#endif
      eps=dfloat(10)**(-20)
      zero=0
      e1=-10000
      e2= 10000
#ifdef __TIMER_DO__
      call timer_sta(815)
#endif
      do 10 is=1,nspin
      do 10 i=1,neig
      do 10 ip2=1,np2
        e0=eig2(ip2,i,is)-eferm
        if(e0.lt.zero .and. e0.gt.(e1-eferm)) then
          e1=eig2(ip2,i,is)
        else if(e0.gt.zero .and. e0.lt.(e2-eferm)) then
          e2=eig2(ip2,i,is)
        end if
   10 continue
#ifdef __TIMER_DO__
      call timer_end(815)
#endif
c     eferm= e1+eps
      eferm=(e1+e2)/2
c     eferm= e2-eps
      if(ipri .ge.2) write( 6
     &     ,'('' e(val.top)='',f12.6,5x,''e(cond.bottom)='',f12.6)')
     &      e1,e2
#ifdef __TIMER_SUB__
      call timer_end(712)
#endif
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine fermi4(np2,lmnp2,lmneig,neig,eig2,eferm)
c
      implicit real*8(a-h,o-z)
      dimension eig2(lmnp2,lmneig)
c
      eps=dfloat(10)**(-20)
      zero=0
      e1=-10000
      e2= 10000
c$$$      do 10 is=1,nspin
      do 10 i=1,neig
      do 10 ip2=1,np2
        e0=eig2(ip2,i)-eferm
        if(e0.lt.zero .and. e0.gt.(e1-eferm)) then
          e1=eig2(ip2,i)
        else if(e0.gt.zero .and. e0.lt.(e2-eferm)) then
          e2=eig2(ip2,i)
        end if
   10 continue
c     eferm= e1+eps
      eferm=(e1+e2)/2
c     eferm= e2-eps
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine lsatpz(ng1,tb1,lrb1,nat,y,kat,lmnat,iatsym)
c
      implicit real*8(a-h,o-z)
      real*8    tb1(3,48),y(3,nat),z(3),u(3)
      integer lrb1(3,3,48),kat(nat),iatsym(lmnat,48)
c
      eps=dfloat(10)**(-5)
      one=1
c
      write(6,*) ' '
      write(6,*) ' symmetry operations on atoms'
      do 10 k=1,ng1
        do 20 i=1,nat
          z(1)=lrb1(1,1,k)*y(1,i)+lrb1(1,2,k)*y(2,i)+
     &         lrb1(1,3,k)*y(3,i)+tb1(1,k)
          z(2)=lrb1(2,1,k)*y(1,i)+lrb1(2,2,k)*y(2,i)+
     &         lrb1(2,3,k)*y(3,i)+tb1(2,k)
          z(3)=lrb1(3,1,k)*y(1,i)+lrb1(3,2,k)*y(2,i)+
     &         lrb1(3,3,k)*y(3,i)+tb1(3,k)
          do 30 j=1,nat
            jj=j
            u(1)=dabs(z(1)-y(1,j))
            u(2)=dabs(z(2)-y(2,j))
            u(3)=dmod(dabs(z(3)-y(3,j))+eps/2,one)
            if(u(1).lt.eps .and. u(2).lt.eps .and. u(3).lt.eps) then
              if(kat(i).eq.kat(j)) then
                go to 34
              else
                go to 32
              end if
            end if
   30     continue
          write(6,*) ' k=',k
          write(6,300) i,(y(l,i),l=1,3),(z(l),l=1,3)
          write(6,*) ' === sub.lsatpz. (no atom) ==='
          stop '=lsatpz(atomic position)='
c
   32     continue
          write(6,*) ' k=',k
          write(6,320) i,(y(l,i),l=1,3),jj,(z(l),l=1,3)
          write(6,*) ' === sub.lsatpz. (different kind) ==='
          stop '=lsatpz(atomic position)='
c
   34     continue
          iatsym(i,k)=jj
   20   continue
        write(6,100) k,(iatsym(i,k),i=1,nat)
   10 continue
      return
  300 format(' iatom=',i5,'  (',3f7.3,')  --->  no atom',
     &                 5x,'  (',3f7.3,')')
  320 format(' iatom=',i5,'  (',3f7.3,')  --->   iatom=',
     &                 i5,'  (',3f7.3,')')
  100 format(' k=',i2,2x,20i3/(7x,20i3))
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine lsdos0(engy,nzz,eig,ddd,sss,eee)
c
      implicit real*8(a-h,o-z)
      real*8 eig(nzz+1)
c
      eps=dfloat(10)**(-3)
      ddd=0
      sss=0
      eee=0
      do 10 i=1,nzz
        if(eig(i+1).ge.eig(i)) then
          if(engy.le.eig(i))then
          else if(engy.le.eig(i+1)) then
            ww=eig(i+1)-eig(i)
            w1=engy-eig(i)
            if(ww.ge.eps) then
              ddd=ddd+1/ww
              sss=sss+w1/ww
            else
              ddd=ddd+1/eps
              sss=sss+1
            end if
          else
            sss=sss+1
          end if
c
        else
          if(engy.le.eig(i+1))then
          else if(engy.le.eig(i)) then
            ww=eig(i)-eig(i+1)
            w1=engy-eig(i+1)
            if(ww.ge.eps) then
              ddd=ddd+1/ww
              sss=sss+w1/ww
            else
              ddd=ddd+1/eps
              sss=sss+1
            end if
          else
            sss=sss+1
          end if
        end if
c
   10 continue
      ddd=ddd/nzz
      sss=sss/nzz
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.lsnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia)
c
c#12  input:     jpr     : print control
c#12             rm      : maximum of neighbor distance (standard unit)
c#12             grb(3,3): metric tensor
c#12             nat     : # of atoms in a unit cell
c#12             u(3,nat): atomic positions in a unit cell
c#12
c#12  output:    nb(nat) : # of neighbors
c#12             lmnb    : limit of # of neighbors
c#12             r (lmnb,nat): distance to the neighbors (standard unit)
c#12             ia(lmnb,nat): atom index of the neighbor
c#12  noexternal:
c
c#21  to get a distance from each atom to the neighboring atoms
c
c#31  1991.10.16.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine lsnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia,w)
c
      implicit real*8(a-h,o-z)
      real*8 grb(3,3)
      real*8 u(3,nat), r(lmnb,nat), w(3,lmnb,nat)
      integer ia(lmnb,nat), nb(nat)
c
      eps=dfloat(10)**(-5)
      one=1
      pai2=8*datan(one)
      pai22=pai2**2

      nz=int(rm/sqrt(grb(3,3)))+3
c
      do 10 j0=1,nat
        v1=u(1,j0)
        v2=u(2,j0)
        v3=u(3,j0)
        k=0
        do 20 j1=1,nat
          ww1=u(1,j1)
          ww2=u(2,j1)
          ww3=u(3,j1)
          w1=ww1-v1
          w2=ww2-v2
          do 22 iz=-nz,nz
            w3=ww3+iz-v3
            rr=sqrt(w1*grb(1,1)*w1+w1*grb(1,2)*w2+w1*grb(1,3)*w3
     &             +w2*grb(2,1)*w1+w2*grb(2,2)*w2+w2*grb(2,3)*w3
     &             +w3*grb(3,1)*w1+w3*grb(3,2)*w2+w3*grb(3,3)*w3)
c                                                 registration
            if(rr.gt.eps .and. rr.le.rm) then
              k=k+1
              if(k.gt.lmnb) then
                write(6,*) ' k=',k,'  > lmnb=',lmnb
                stop 'sub.lsnbat'
              end if
              r (k,j0)=rr
              ia(k,j0)=j1
              w(1,k,j0)=w1
              w(2,k,j0)=w2
              w(3,k,j0)=w3
c             write(6,900) j0,k,ia(k,j0),r(k,j0)
            end if
   24     continue
   22     continue
   20   continue
        nb(j0)=k
c
   10 continue
c
      do 90 j0=1,nat
        call lsnbst(nb(j0),r(1,j0),ia(1,j0),w(1,1,j0))
        if(jpr.ge.2) then
          write(6,*) ' atom0=',j0
          do 92  k=1,nb(j0)
          write(6,900) k,ia(k,j0),(w(i,k,j0),i=1,3),r(k,j0)
   92     continue
        end if
   90 continue
  900 format(1h ,'atom1=',i5,' (',i5,')','  w=',3f8.3,'   r=',f10.3)
      return
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.lsnbst(n,ra,ia,wa)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  ia(n): an array to be rearranged correspondingly
c#12                  wa(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#21  to sort ra(n) in ascending order with an index ia(n)
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine lsnbst(n,ra,ia,wa)
c
      implicit real*8(a-h,o-z)
      real*8  ra(n),wa(3,n)
      integer ia(n)
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          iia=ia(l)
          w1a=wa(1,l)
          w2a=wa(2,l)
          w3a=wa(3,l)
        else
          rra=ra(ir)
          iia=ia(ir)
          w1a=wa(1,ir)
          w2a=wa(2,ir)
          w3a=wa(3,ir)
          ra(ir)=ra(1)
          ia(ir)=ia(1)
          wa(1,ir)=wa(1,1)
          wa(2,ir)=wa(2,1)
          wa(3,ir)=wa(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            ia(1)=iia
            wa(1,1)=w1a
            wa(2,1)=w2a
            wa(3,1)=w3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            ia(i)=ia(j)
            wa(1,i)=wa(1,j)
            wa(2,i)=wa(2,j)
            wa(3,i)=wa(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        ia(i)=iia
        wa(1,i)=w1a
        wa(2,i)=w2a
        wa(3,i)=w3a
c
      go to 10
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine lstt1i(engy,nzz,eig,cdos,cind)
c
      implicit real*8(a-h,o-z)
      real*8 eig(nzz+1), cdos(nzz+1), cind(nzz+1)
c
      eps=dfloat(10)**(-3)
      one=1
      half=one/2
      do 10 i=1,nzz+1
        cdos(i)=0
        cind(i)=0
   10 continue
c
      do 20 i=1,nzz
        if(eig(i+1).ge.eig(i)) then
          if(engy.le.eig(i))then
          else if(engy.le.eig(i+1)) then
            ww=eig(i+1)-eig(i)
            if(ww.ge.eps) then
              w1=engy-eig(i)
              w2=eig(i+1)-engy
              www=ww**2
              ww2=w2/www
              ww1=w1/www
              cdos(i  )=cdos(i  )+ww2
              cdos(i+1)=cdos(i+1)+ww1
              cind(i  )=cind(i  )+(ww+w2)*ww1/2
              cind(i+1)=cind(i+1)+w1*ww1/2
            else
              cdos(i  )=cdos(i  )+half/eps
              cdos(i+1)=cdos(i+1)+half/eps
              cind(i  )=cind(i  )+half
              cind(i+1)=cind(i+1)+half
            end if
          else
            cind(i  )=cind(i  )+half
            cind(i+1)=cind(i+1)+half
          end if
c
        else
          if(engy.le.eig(i+1))then
          else if(engy.le.eig(i)) then
            ww=eig(i)-eig(i+1)
            if(ww.ge.eps) then
              w1=engy-eig(i+1)
              w2=eig(i)-engy
              www=ww**2
              ww2=w2/www
              ww1=w1/www
              cdos(i  )=cdos(i  )+ww1
              cdos(i+1)=cdos(i+1)+ww2
              cind(i  )=cind(i  )+w1*ww1/2
              cind(i+1)=cind(i+1)+(ww+w2)*ww1/2
            else
              cdos(i  )=cdos(i  )+half/eps
              cdos(i+1)=cdos(i+1)+half/eps
              cind(i  )=cind(i  )+half
              cind(i+1)=cind(i+1)+half
            end if
          else
            cind(i  )=cind(i  )+half
            cind(i+1)=cind(i+1)+half
          end if
        end if
   20 continue
c
      do 30 i=1,nzz+1
        cdos(i)=cdos(i)/nzz
        cind(i)=cind(i)/nzz
   30 continue
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msatpz(ng1,tb1,lrb1,nat,y,kat,lmnat,iatsym)
c
      implicit real*8(a-h,o-z)
      real*8  tb1(3,48),y(3,nat),z(3),u(3)
      integer lrb1(3,3,48),kat(nat),iatsym(lmnat,48)
c
      eps=dfloat(10)**(-5)
      one=1
c
      write(6,*) ' '
      write(6,*) ' symmetry operations on atoms'
      do 10 k=1,ng1
        do 20 i=1,nat
          z(1)=lrb1(1,1,k)*y(1,i)+lrb1(1,2,k)*y(2,i)+
     &         lrb1(1,3,k)*y(3,i)+tb1(1,k)
          z(2)=lrb1(2,1,k)*y(1,i)+lrb1(2,2,k)*y(2,i)+
     &         lrb1(2,3,k)*y(3,i)+tb1(2,k)
          z(3)=lrb1(3,1,k)*y(1,i)+lrb1(3,2,k)*y(2,i)+
     &         lrb1(3,3,k)*y(3,i)+tb1(3,k)
          do 30 j=1,nat
            jj=j
            u(1)=dmod(dabs(z(1)-y(1,j))+eps/2,one)
            u(2)=dmod(dabs(z(2)-y(2,j))+eps/2,one)
            u(3)=dabs(z(3)-y(3,j))
            if(u(1).lt.eps .and. u(2).lt.eps .and. u(3).lt.eps) then
              if(kat(i).eq.kat(j)) then
                go to 34
              else
                go to 32
              end if
            end if
   30     continue
          write(6,*) ' k=',k
          write(6,300) i,(y(l,i),l=1,3),(z(l),l=1,3)
          write(6,*) ' === sub.msatpz. (no atom) ==='
          stop '=msatpz(atomic position)='
c
   32     continue
          write(6,*) ' k=',k
          write(6,320) i,(y(l,i),l=1,3),jj,(z(l),l=1,3)
          write(6,*) ' === sub.msatpz. (different kind) ==='
          stop '=msatpz(atomic position)='
c
   34     continue
          iatsym(i,k)=jj
   20   continue
        write(6,100) k,(iatsym(i,k),i=1,nat)
   10 continue
      return
  300 format(' iatom=',i5,'  (',3f7.3,')  --->  no atom',
     &                 5x,'  (',3f7.3,')')
  320 format(' iatom=',i5,'  (',3f7.3,')  --->   iatom=',
     &                 i5,'  (',3f7.3,')')
  100 format(' k=',i2,2x,20i3/(7x,20i3))
      end

c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msdos0(e,nx,ny,ea,dos,dm0,dm1)
c
c     ***  by linear interpolation in triangles             ***
c     ***
c     ***
c     ***  nx  number of mesh points in b.z. in x-direction
c     ***  ny  number of mesh points in b.z. in y-direction
c     ***  dos density of states
c     ***  dm0 0th moment of dos (integrated up to e)
c     ***  dm1 1st moment of dos (integrated up to e)
c
      implicit  real*8(a-h,o-z)
c
      dimension ea(*)
      dimension ecu(2,2), ec(4), et(4), eb(4)
      dimension icu(2,2), ic(4), it(4), ib(4)
      equivalence(ec(1),ecu(1,1))
      equivalence(ic(1),icu(1,1))
      integer iqt(3,4)
      data    iqt/1,2,4, 1,3,4,  2,1,3, 2,4,3/
c
      one=1
      c3=one/3
      npx=nx+1
      npy=ny+1
      np=npx*npy
      nc=nx*ny
      nt=4*nc
c
c     np: the number of mesh points
c     nc: the number of squares
c     nt: the number of triangles (consider four triangles in a square)
c
      dos=0
      dm0=0
      dm1=0
c
      emax0=-1.e30
      emin0= 1.e30
      do 12 ip=1,np
        if(ea(ip).gt.emax0) emax0=ea(ip)
        if(ea(ip).lt.emin0) emin0=ea(ip)
   12 continue
c                                   ============ if 1 ==
      if(e.gt.emin0) then
c
c     ***       sampling over squares      ***
c
      jc=0
      do 20 iy=0,ny-1
      do 20 ix=0,nx-1
        jc=jc+1
c       ***  energies at square corners      ***
        ni=npx*iy+ix
        emax=-1.e30
        emin= 1.e30
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(ky-1)+kx
        ecu(kx,ky)=ea(ip0)
        icu(kx,ky)=ip0
        if(ea(ip0).gt.emax) emax=ea(ip0)
        if(ea(ip0).lt.emin) emin=ea(ip0)
   30   continue
c                                   ============ if 2 ==
        if(e.gt.emin) then
c         *** two triangles  for each type   ***
c         *** sampling over triangles        ***
          do 40 jt=1,4
            et(1)=ec(iqt(1,jt))
            et(2)=ec(iqt(2,jt))
            et(3)=ec(iqt(3,jt))
            it(1)=ic(iqt(1,jt))
            it(2)=ic(iqt(2,jt))
            it(3)=ic(iqt(3,jt))
            do 42 m=1,3
            eb(m)=et(m)
            ib(m)=it(m)
   42       continue
c           ***  eb(1).le.eb(2).le.eb(3)  ***
            call msttio(eb,ib)
            e1=eb(1)
            e2=eb(2)
            e3=eb(3)
            call msttie(e1,e2,e3)
            call msdosi(e,e1,e2,e3,d,d0,d1)
            dos=dos+d
            dm0=dm0+d0
            dm1=dm1+d1
   40     continue
        end if
c                                   ============ if 2 ==
   20 continue
c
      dos=dos/nt
      dm0=dm0/nt
      dm1=dm1/nt
c
      end if
c                                   ============ if 1 ==
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msdosi(e,e1,e2,e3,d,d0,d1)
c
c  ** triangle corners at one energy           **  c
c  ** according to lambin and vigneron,        **  c
c  ** phys. rev. b29, 3430 (1984)              **  c
c
      implicit real*8(a-h,o-z)
      d21=e2-e1
      d31=e3-e1
      d32=e3-e2
      d1=e-e1
      d3=e3-e
      if(d21.lt.0 .or. d31.lt.0 .or. d32.lt.0) then
        write(6,*) ' e1,e2,e3=',e1,e2,e3
        write(6,*) ' energy order error  in sub.msttii.'
        stop ' = sub.msttii ='
      end if
c
      if(e.le.e1) then
        d =0
        d0=0
        d1=0
c
      else if(e.le.e2) then
        b2=d1/d21
        b3=d1/d31
        d =2*b2/d31
        d0=  b2*b3
        d1=  b2*b3*(2*e+e1)/3
c
      else if(e.lt.e3) then
        d3=e3-e
        b1=d3/d31
        b2=d3/d32
        d =2*b1/d32
        d0=1-b1*b2
        d1=(e1+e2+e3-b1*b2*(2*e+e3))/3
c
      else
        d =0
        d0=1
        d1=(e1+e2+e3)/3
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.msnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia)
c
c#12  input:     jpr     : print control
c#12             rm      : maximum of neighbor distance (standard unit)
c#12             grb(3,3): metric tensor
c#12             nat     : # of atoms in a unit cell
c#12             u(3,nat): atomic positions in a unit cell
c#12
c#12  output:    nb(nat) : # of neighbors
c#12             lmnb    : limit of # of neighbors
c#12             r (lmnb,nat): distance to the neighbors (standard unit)
c#12             ia(lmnb,nat): atom index of the neighbor
c#12  noexternal:
c
c#21  to get a distance from each atom to the neighboring atoms
c
c#31  1990.06.05.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia,w)
c
      implicit real*8(a-h,o-z)
      real*8  grb(3,3)
      real*8  u(3,nat), r(lmnb,nat), w(3,lmnb,nat)
      integer ia(lmnb,nat), nb(nat)
c
      eps=dfloat(10)**(-5)
      one=1
      pai2=8*datan(one)
      pai22=pai2**2
c
      nx=0
      ny=0
      do 32 iy=-1,1
      do 34 ix=-1,1
      if(ix.ne.0 .and. iy.ne.0) then
        p1=ix
        p2=iy
        rr=sqrt(p1*grb(1,1)*p1+p1*grb(1,2)*p2
     &         +p2*grb(2,1)*p1+p2*grb(2,2)*p2)
        nn=(rm/rr)+5
        nnx=iabs(ix)*nn
        nny=iabs(iy)*nn
        if(nnx.gt.nx) nx=nnx
        if(nny.gt.ny) ny=nny
      endif
   34 continue
   32 continue
c
c     write(6,'('' nx,ny='',3i5)') nx,ny
c
c
      do 10 j0=1,nat
        v1=u(1,j0)
        v2=u(2,j0)
        v3=u(3,j0)
        k=0
        do 20 j1=1,nat
          ww1=u(1,j1)
          ww2=u(2,j1)
          ww3=u(3,j1)
          w3=ww3-v3
          do 22 iy=-ny,ny
            w2=ww2+iy-v2
          do 24 ix=-nx,nx
            w1=ww1+ix-v1
            rr=sqrt(w1*grb(1,1)*w1+w1*grb(1,2)*w2+w1*grb(1,3)*w3
     &             +w2*grb(2,1)*w1+w2*grb(2,2)*w2+w2*grb(2,3)*w3
     &             +w3*grb(3,1)*w1+w3*grb(3,2)*w2+w3*grb(3,3)*w3)
c                                                 registration
            if(rr.gt.eps .and. rr.le.rm) then
              k=k+1
              if(k.gt.lmnb) then
                write(6,*) ' k=',k,'  > lmnb=',lmnb
                stop 'sub.msnbat'
              end if
              r (k,j0)=rr
              ia(k,j0)=j1
              w(1,k,j0)=w1
              w(2,k,j0)=w2
              w(3,k,j0)=w3
c             write(6,900) j0,k,ia(k,j0),r(k,j0)
            end if
   24     continue
   22     continue
   20   continue
        nb(j0)=k
c
   10 continue
c
      do 90 j0=1,nat
        call msnbst(nb(j0),r(1,j0),ia(1,j0),w(1,1,j0))
        if(jpr.ge.2) then
          write(6,*) ' atom0=',j0
          do 92  k=1,nb(j0)
          write(6,900) k,ia(k,j0),(w(i,k,j0),i=1,3),r(k,j0)
   92     continue
        end if
   90 continue
  900 format(1h ,'atom1=',i5,' (',i5,')','  w=',3f8.3,'   r=',f10.6)
      return
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.msnbst(n,ra,ia,wa)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  ia(n): an array to be rearranged correspondingly
c#12                  wa(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#21  to sort ra(n) in ascending order with an index ia(n)
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msnbst(n,ra,ia,wa)
c
      implicit real*8(a-h,o-z)
      real*8  ra(n),wa(3,n)
      integer ia(n)
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          iia=ia(l)
          w1a=wa(1,l)
          w2a=wa(2,l)
          w3a=wa(3,l)
        else
          rra=ra(ir)
          iia=ia(ir)
          w1a=wa(1,ir)
          w2a=wa(2,ir)
          w3a=wa(3,ir)
          ra(ir)=ra(1)
          ia(ir)=ia(1)
          wa(1,ir)=wa(1,1)
          wa(2,ir)=wa(2,1)
          wa(3,ir)=wa(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            ia(1)=iia
            wa(1,1)=w1a
            wa(2,1)=w2a
            wa(3,1)=w3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            ia(i)=ia(j)
            wa(1,i)=wa(1,j)
            wa(2,i)=wa(2,j)
            wa(3,i)=wa(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        ia(i)=iia
        wa(1,i)=w1a
        wa(2,i)=w2a
        wa(3,i)=w3a
c
      go to 10
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine mstt1i(e,nx,ny,ea,cdos,css0)
c
c     ***  In order to obtain 'the density of states and    ***
c     ***  the integrated density of states',               ***
c     ***  weights for each k point are calculated          ***
c     ***  for a single energy                              ***
c     ***  by linear interpolation in triangles             ***
c     ***
c     ***
c     ***  npx  number of mesh points in b.z.
c     ***         in x-direction
c     ***  npy  number of mesh points in b.z.
c     ***         in y-direction
c     ***  cdos density of states 'coefficient'
c     ***  cind number  of states 'coefficient'
c     ***  ref.)  ph. lambin and j. p. vigneron,
c                 phys. rev. b29 (1984) 3430.
c
      implicit  real*8(a-h,o-z)
c
      dimension ea(*),cdos(*),css0(*)
      dimension ecu(2,2), ec(4), et(4), eb(4)
      dimension icu(2,2), ic(4), it(4), ib(4)
      dimension dos(3),ss0(3)
      equivalence(ec(1),ecu(1,1))
      equivalence(ic(1),icu(1,1))
      integer iqt(3,4)
      data    iqt/1,2,4, 1,3,4,  2,1,3, 2,4,3/
c
      one=1
      c3=one/3
      npx=nx+1
      npy=ny+1
      np=npx*npy
      nc=nx*ny
      nt=4*nc
c
c     np: the number of mesh points
c     nc: the number of squares
c     nt: the number of triangles (consider four triangles in a square)
c
      do 10 ip=1,np
        cdos(ip)=0
        css0(ip)=0
   10 continue
c
      emax0=-1.d30
      emin0= 1.d30
      do 12 ip=1,np
        if(ea(ip).gt.emax0) emax0=ea(ip)
        if(ea(ip).lt.emin0) emin0=ea(ip)
   12 continue
c                                   ============ if 1 ==
      if(e.gt.emin0) then
c
c     ***       sampling over squares      ***
c
      jc=0
      do 20 iy=0,ny-1
      do 20 ix=0,nx-1
        jc=jc+1
c       ***  energies at square corners      ***
        ni=npx*iy+ix
        emax=-1.d30
        emin= 1.d30
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(ky-1)+kx
        ecu(kx,ky)=ea(ip0)
        icu(kx,ky)=ip0
        if(ea(ip0).gt.emax) emax=ea(ip0)
        if(ea(ip0).lt.emin) emin=ea(ip0)
   30   continue
c                                   ============ if 2 ==
        if(e.ge.emax) then
          css0(ic(1))=css0(ic(1))+1
          css0(ic(2))=css0(ic(2))+1
          css0(ic(3))=css0(ic(3))+1
          css0(ic(4))=css0(ic(4))+1
c
        else if(e.gt.emin) then
c         *** two triangles  for each type   ***
c         *** sampling over triangles        ***
          do 40 jt=1,4
            et(1)=ec(iqt(1,jt))
            et(2)=ec(iqt(2,jt))
            et(3)=ec(iqt(3,jt))
            it(1)=ic(iqt(1,jt))
            it(2)=ic(iqt(2,jt))
            it(3)=ic(iqt(3,jt))
            do 42 m=1,3
            eb(m)=et(m)
            ib(m)=it(m)
   42       continue
c           ***  eb(1).le.eb(2).le.eb(3)  ***
            call msttio(eb,ib)
            e1=eb(1)
            e2=eb(2)
            e3=eb(3)
            call msttie(e1,e2,e3)
            call msttii(e,e1,e2,e3,dos,ss0)
            cdos(ib(1))=cdos(ib(1))+dos(1)
            cdos(ib(2))=cdos(ib(2))+dos(2)
            cdos(ib(3))=cdos(ib(3))+dos(3)
            css0(ib(1))=css0(ib(1))+ss0(1)
            css0(ib(2))=css0(ib(2))+ss0(2)
            css0(ib(3))=css0(ib(3))+ss0(3)
   40     continue
        end if
c                                   ============ if 2 ==
   20 continue
c
      do 50 i=1,np
      cdos(i)=cdos(i)/nt
      css0(i)=css0(i)/nt
   50 continue
c
      end if
c                                   ============ if 1 ==
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msttie(e1,e2,e3)
c
      implicit  real*8(a-h,o-z)
c
      eps=dfloat(10)**(-3)
c
      if(abs(e2-e1).lt.eps) e1=e2-eps
      if(abs(e3-e2).lt.eps) e3=e2+eps
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msttii(e,e1,e2,e3,dos,ss0)
c
c  ** dos and integrated dos 'coefficients' at **  c
c  ** triangle corners at one energy           **  c
c  ** according to lambin and vigneron,        **  c
c  ** phys. rev. b29, 3430 (1984)              **  c
c
      implicit real*8(a-h,o-z)
      dimension dos(3),ss0(3)
c
      one=1
      d21=e2-e1
      d31=e3-e1
      d32=e3-e2
      d1=e-e1
      d3=e3-e
      if(d21.lt.0.0 .or. d31.lt.0.0 .or. d32.lt.0.0) then
        write(6,*) ' e1,e2,e3=',e1,e2,e3
        write(6,*) ' energy order error  in sub.msttii.'
        stop ' = sub.msttii ='
      end if
c
      if(e.le.e1) then
        do 10 i=1,3
          dos(i)=0
          ss0(i)=0
   10   continue
c
      else if(e.le.e2) then
        d2=e2-e
        d3=e3-e
        a2=d2/d21
        a3=d3/d31
        aa=d1/(d21*d31)
        b2=d1/d21
        b3=d1/d31
        dos(1)=(a2+a3)*aa
        dos(2)= b2    *aa
        dos(3)=    b3 *aa
        ss0(1)=(1+a2+a3)*d1*aa/3
        ss0(2)=   b2    *d1*aa/3
        ss0(3)=      b3 *d1*aa/3
c
      else if(e.lt.e3) then
        d2=e-e2
        d3=e3-e
        a1=d1/d31
        a2=d2/d32
        aa=d3/(d31*d32)
        b1=d3/d31
        b2=d3/d32
        dos(1)= b1    *aa
        dos(2)=    b2 *aa
        dos(3)=(a1+a2)*aa
        ss0(1)=(1-b1*d3*aa)/3
        ss0(2)=(1-b2*d3*aa)/3
        ss0(3)=(1-(1+a1+a2)*d3*aa)/3
c
      else
        c3=one/3
        do 20 i=1,3
          dos(i)=0
          ss0(i)=c3
   20   continue
      end if
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine msttio(eb,ib)
c
c     ***  eb(1).le.eb(2).le.eb(3)  ***
c
      implicit real*8(a-h,o-z)
      real*8 eb(3)
      integer ib(3)
          do 44 m=2,3
          if(eb(m).lt.eb(1)) then
            ebm=eb(m)
            ibm=ib(m)
            eb(m)=eb(1)
            ib(m)=ib(1)
            eb(1)=ebm
            ib(1)=ibm
          end if
   44     continue
          if(eb(3).le.eb(2)) then
            ebm=eb(3)
            ibm=ib(3)
            eb(3)=eb(2)
            ib(3)=ib(2)
            eb(2)=ebm
            ib(2)=ibm
          end if
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrrpc2(euler,lmax,lmlmax,dc)
c
      implicit real*8(a-h,o-z)
      parameter (lmaxm=2)
c
      real*8     euler(3)
      complex*16 dl ((lmaxm+1)**2,(lmaxm+1)**2)
      complex*16 tlc((lmaxm+1)**2,(lmaxm+1)**2),wk
      real*8     dc ((lmlmax+1)**2,(lmlmax+1)**2)
c
      if(lmax.gt.lmaxm) then
        write(6,*) ' lmax=',lmax,'  : too large'
      end if
      call nrrpl2(euler,lmax,lmaxm,dl)
      call nrtrc2(lmax,lmaxm,tlc)
c
      zero=0
      imax=(lmax+1)**2
      do 10 j=1,imax
      do 10 i=1,imax
        wk=dcmplx(zero,zero)
        do 20 k1=1,imax
        do 20 k2=1,imax
   20   wk=wk+conjg(tlc(k1,i))*dl(k1,k2)*tlc(k2,j)
        if(abs(dimag(wk)).gt.1.e-5) then
          write(6,*) ' d(',i,',',j,') =',wk,'  should be real.'
          stop '=nrrpc2(dc)='
        end if
        dc(i,j)=dreal(wk)
   10 continue
c
      return
      end
c member name  nrrpj0
c     implicit real*8(a-h,o-z)
c     real*8     euler(3)
c     complex*16 wd
c     one=1
c     pai = datan(one)*4
c     j=1
c     jd=1
c     euler(1)=pai
c     euler(2)=0.0
c     euler(3)=0.0
c     do 10 m1=j,-j,-jd
c     do 20 m2=j,-j,-jd
c       call nrrpj0(euler,j,jd,m1,m2,wd)
c       write(6,*) wd
c  20 continue
c  10 continue
c     stop
c     end
c subroutine nrrpj0 ===*====3====*====4====*====5====*====6====*====7
c
c     representation of rotation group.
c          j/jd :  angular momentum.
c
c      1983.7.15. :  n. hamada.
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrrpj0(euler,j,jd,m1,m2,wd)
c
      implicit real*8(a-h,o-z)
c
      real*8     euler(3)
      complex*16 wd, ea, ec
      integer    ic(0:15)
c
      eps=dfloat(10)**(-5)
      zero=0
c
      ipr=0
      if(jd.eq.1) then
        call tsrmi(j,m1,m2,k1,k2,k3,is,ic,ipr)
      else if(jd.eq.2) then
        call tsrmhi(j,m1,m2,k1,k2,k3,is,ic,ipr)
      else
        stop ' === stop in sub.pgchjr. (jd.ne. 1 nor 2) ==='
      end if
c
      cb=cos(euler(2))
      sb=sin(euler(2))
      if(k1.eq.0 .or. k3.eq.0) then
        wd=dcmplx(zero,zero)
      else if(is.eq.2 .and. dabs(sb).lt.eps) then
        wd=dcmplx(zero,zero)
      else if(is.eq.1 .and. dabs(1+cb).lt.eps) then
        wd=dcmplx(zero,zero)
      else if(is.eq.3 .and. dabs(1-cb).lt.eps) then
        wd=dcmplx(zero,zero)
      else
        r1=dfloat(k3)
        r1=(dfloat(k1)/k2)*sqrt(r1)
        r2=dfloat(ic(0))
        do 20 i=1,15
          if(ic(i).ne.0) then
            r2=r2+ic(i)*(cb**i)
          end if
   20   continue
        if(is.eq.0) then
          r3=1
        else if(is.eq.2) then
          r3=sb
        else if(is.eq.1) then
          r3=sqrt((1+cb)/2)
        else if(is.eq.3) then
          r3=sqrt((1-cb)/2)
        end if
        ca=cos((euler(1)*m1)/jd)
        cc=cos((euler(3)*m2)/jd)
        sa=-sin((euler(1)*m1)/jd)
        sc=-sin((euler(3)*m2)/jd)
        ea=dcmplx(ca,sa)
        ec=dcmplx(cc,sc)
        wd=ea*(r1*r2*r3)*ec
      end if
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrrpl1(euler,inver,lmax,lmlmax,d)
c
      implicit real*8(a-h,o-z)
c
      real*8     euler(3)
      complex*16 wd,d((lmlmax+1)**2,(lmlmax+1)**2)
c
      zero=0
      do 10 i=1,(lmax+1)**2
      do 10 j=1,(lmax+1)**2
   10 d(i,j)=dcmplx(zero,zero)
c
c     write(6,*) 'euler=',(euler(i),i=1,3)
      jd=1
      do 20 l=0,lmax
        do 22 m2=-l,l
          lm2=l**2+m2+l+1
        do 22 m1=-l,l
          lm1=l**2+m1+l+1
          call nrrpj0(euler,l,jd,m1,m2,wd)
          d(lm1,lm2)=wd*inver
c         write(6,*) lm1,lm2,d(lm1,lm2)
   22   continue
   20 continue
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrrpl2(euler,lmax,lmlmax,d)
c
      implicit real*8(a-h,o-z)
c
      real*8     euler(3)
      complex*16 wd,d((lmlmax+1)**2,(lmlmax+1)**2)
c
      zero=0
      do 10 i=1,(lmax+1)**2
      do 10 j=1,(lmax+1)**2
   10 d(i,j)=dcmplx(zero,zero)
c
c     write(6,*) 'euler=',(euler(i),i=1,3)
      jd=1
      do 20 l=0,lmax
        do 22 m2=-l,l
          lm2=l**2+m2+l+1
        do 22 m1=-l,l
          lm1=l**2+m1+l+1
          call nrrpj0(euler,l,jd,m1,m2,wd)
          d(lm1,lm2)=wd
c         write(6,*) lm1,lm2,wd
   22   continue
   20 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nrskc0(jpr,mm0,mm1,xx,v,hmx)
c
c#12  input :         jpr   : output control
c#12                  mm0   : orbital index
c#12                  mm1   : orbital index
c#12                  xx(3) : position of (mm1) relative to (mm0)
c#12                          (Cartesian coordinate)
c#12  output:         hmx   : transfer matrix element
c#13  noexternal
c
c#21  to get Koster-Slater coefficients.
c#22
c#22  orbital index  1: s
c#22                 2: x,  3: y,  4: z
c#22                 5: xy, 6: yz, 7: zx, 8: x**2-y**2, 9: 3*z**2-r**2
c#22
c#22  bond index     1: ss, 2: sp, 3:pp1, 4:pp2,
c#22                 5: sd, 6: pd1, 7: pd2, 8: dd1, 9: dd2, 10: dd3
c
c#31  1990.06.07.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrskc0(jpr,mm0,mm1,xx,v,hmx)
c
      implicit real*8(a-h,o-z)
      real*8 xx(3),v(10)
      integer ip(9)
c
      data ip/1, -1,-1,-1, 1,1,1,1,1/
c
      three=3
      if(mm0.lt.1 .or. mm0.gt.9 .or. mm1.lt.1 .or. mm1.gt.9) then
        write(6,*) ' mm0,mm1=',mm0,mm1
        stop '=sub.nrskc0='
      end if
      if(mm0.le.mm1) then
        m0=mm0
        m1=mm1
      else
        m0=mm1
        m1=mm0
      end if
c
      a3=sqrt(three)
      xxx=sqrt(xx(1)**2+xx(2)**2+xx(3)**2)
      x1=xx(1)/xxx
      x2=xx(2)/xxx
      x3=xx(3)/xxx
c
      if(m0.eq.1) then
        if(m1.eq.1) then
          hmx=v(1)
        else if(m1.eq.2) then
          hmx=v(2)*x1
        else if(m1.eq.3) then
          hmx=v(2)*x2
        else if(m1.eq.4) then
          hmx=v(2)*x3
        else if(m1.eq.5) then
          hmx=v(5)*x1*x2*a3
        else if(m1.eq.6) then
          hmx=v(5)*x2*x3*a3
        else if(m1.eq.7) then
          hmx=v(5)*x3*x1*a3
        else if(m1.eq.8) then
          hmx=v(5)*(x1*x1-x2*x2)*a3/2
        else if(m1.eq.9) then
          hmx=v(5)*   (x3*x3-(x1*x1+x2*x2)/2)
        end if
      else if(m0.eq.2) then
        if(m1.eq.2) then
          hmx=v(3)*x1*x1+v(4)*(1-x1*x1)
        else if(m1.eq.3) then
          hmx=(v(3)-v(4))*x1*x2
        else if(m1.eq.4) then
          hmx=(v(3)-v(4))*x1*x3
        else if(m1.eq.5) then
          hmx=v(6)*x1*x1*x2*a3+v(7)*(1-2*x1*x1)*x2
        else if(m1.eq.6) then
          hmx=v(6)*x1*x2*x3*a3-v(7)*x1*x2*x3*2
        else if(m1.eq.7) then
          hmx=v(6)*x1*x3*x1*a3+v(7)*(1-2*x1*x1)*x3
        else if(m1.eq.8) then
          hmx=v(6)*x1*(x1*x1-x2*x2)*a3/2
     &       +v(7)*x1*(1-x1*x1+x2*x2)
        else if(m1.eq.9) then
          hmx=v(6)*x1*(x3*x3-(x1*x1+x2*x2)/2)
     &       -v(7)*x1*x3*x3*a3
        end if
      else if(m0.eq.3) then
        if(m1.eq.3) then
          hmx=v(3)*x2*x2+v(4)*(1-x2*x2)
        else if(m1.eq.4) then
          hmx=(v(3)-v(4))*x2*x3
        else if(m1.eq.5) then
          hmx=v(6)*x2*x1*x2*a3+v(7)*(1-2*x2*x2)*x1
        else if(m1.eq.6) then
          hmx=v(6)*x2*x2*x3*a3+v(7)*(1-2*x2*x2)*x3
        else if(m1.eq.7) then
          hmx=v(6)*x2*x3*x1*a3-v(7)*x2*x3*x1*2
        else if(m1.eq.8) then
          hmx=v(6)*x2*(x1*x1-x2*x2)*a3/2
     &       -v(7)*x2*(1+x1*x1-x2*x2)
        else if(m1.eq.9) then
          hmx=v(6)*x2*(x3*x3-(x1*x1+x2*x2)/2)
     &       -v(7)*x2*x3*x3*a3
        end if
      else if(m0.eq.4) then
        if(m1.eq.4) then
          hmx=v(3)*x3*x3+v(4)*(1-x3*x3)
        else if(m1.eq.5) then
          hmx=v(6)*x3*x1*x2*a3-v(7)*x3*x1*x2*2
        else if(m1.eq.6) then
          hmx=v(6)*x3*x2*x3*a3+v(7)*(1-2*x3*x3)*x2
        else if(m1.eq.7) then
          hmx=v(6)*x3*x3*x1*a3+v(7)*(1-2*x3*x3)*x1
        else if(m1.eq.8) then
          hmx=v(6)*x3*(x1*x1-x2*x2)*a3/2
     &       -v(7)*x3*(x1*x1-x2*x2)
        else if(m1.eq.9) then
          hmx=v(6)*x3*(x3*x3-(x1*x1+x2*x2)/2)
     &       +v(7)*x3*(x1*x1+x2*x2)*a3
        end if
      else if(m0.eq.5) then
        if(m1.eq.5) then
          hmx=v(8)*x1*x2*x1*x2*3
     &       +v(9)*(x1*x1+x2*x2-x1*x1*x2*x2*4)
     &       +v(10)*(x3*x3+x1*x1*x2*x2)
        else if(m1.eq.6) then
          hmx=v(8)*x1*x2*x2*x3*3
     &       +v(9)*x1*(1-x2*x2*4)*x3
     &       +v(10)*x1*(x2*x2-1)*x3
        else if(m1.eq.7) then
          hmx=v(8)*x1*x2*x3*x1*3
     &       +v(9)*(1-x1*x1*4)*x2*x3
     &       +v(10)*(x1*x1-1)*x2*x3
        else if(m1.eq.8) then
          hmx=(v(8)*3/2-v(9)*2+v(10)/2)*x1*x2*(x1*x1-x2*x2)
        else if(m1.eq.9) then
          hmx=v(8)*x1*x2*(x3*x3-(x1*x1+x2*x2)/2)*a3
     &       -v(9)*x1*x2*x3*x3*a3
     &       +v(10)*x1*x2*(1+x3*x3)*a3/2
        end if
      else if(m0.eq.6) then
        if(m1.eq.6) then
          hmx=v(8)*x2*x3*x2*x3*3
     &       +v(9)*(x2*x2+x3*x3-x2*x2*x3*x3*4)
     &       +v(10)*(x1*x1+x2*x2*x3*x3)
        else if(m1.eq.7) then
          hmx=v(8)*x2*x3*x3*x1*3
     &       +v(9)*x2*(1-x3*x3*4)*x1
     &       +v(10)*x2*(x3*x3-1)*x1
        else if(m1.eq.8) then
          hmx=v(8)*x2*x3*(x1*x1-x2*x2)*3/2
     &       -v(9)*x2*x3*(1+(x1*x1-x2*x2)*2)
     &       +v(10)*x2*x3*(1+(x1*x1-x2*x2)/2)
        else if(m1.eq.9) then
          hmx=v(8)*x2*x3*(x3*x3-(x1*x1+x2*x2)/2)*a3
     &       +v(9)*x2*x3*(x1*x1+x2*x2-x3*x3)*a3
     &       -v(10)*x2*x3*(x1*x1+x2*x2)*a3/2
        end if
      else if(m0.eq.7) then
        if(m1.eq.7) then
          hmx=v(8)*x3*x1*x3*x1*3
     &       +v(9)*(x3*x3+x1*x1-x3*x3*x1*x1*4)
     &       +v(10)*(x2*x2+x3*x3*x1*x1)
        else if(m1.eq.8) then
          hmx=v(8)*x3*x1*(x1*x1-x2*x2)*3/2
     &       +v(9)*x3*x1*(1-(x1*x1-x2*x2)*2)
     &       -v(10)*x3*x1*(1-(x1*x1-x2*x2)/2)
        else if(m1.eq.9) then
          hmx=v(8)*x3*x1*(x3*x3-(x1*x1+x2*x2)/2)*a3
     &       +v(9)*x3*x1*(x1*x1+x2*x2-x3*x3)*a3
     &       -v(10)*x3*x1*(x1*x1+x2*x2)*a3/2
        end if
      else if(m0.eq.8) then
        if(m1.eq.8) then
          hmx=v(8)*(x1*x1-x2*x2)**2*3/4
     &       +v(9)*(x1*x1+x2*x2-(x1*x1-x2*x2)**2)
     &       +v(10)*(x3*x3+(x1*x1-x2*x2)**2/4)
        else if(m1.eq.9) then
          hmx=v(8)*(x1*x1-x2*x2)*(x3*x3-(x1*x1+x2*x2)/2)*a3/2
     &       -v(9)*(x1*x1-x2*x2)*x3*x3*a3
     &       +v(10)*(x1*x1-x2*x2)*(1+x3*x3)*a3/4
        end if
      else if(m0.eq.9) then
        if(m1.eq.9) then
          hmx=v(8)*(x3*x3-(x1*x1+x2*x2)/2)**2
     &       +v(9)*x3*x3*(x1*x1+x2*x2)*3
     &       +v(10)*(x1*x1+x2*x2)**2*3/4
        end if
      end if
c
      if(mm0.gt.mm1) then
        hmx=hmx*ip(mm0)*ip(mm1)
      end if
      if(jpr.ge.3) write(6,100) mm0,mm1,(xx(i),i=1,3),hmx
  100 format(1h ,2i5,3f6.3,f9.3)
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrskc1(jpr,mm0,mm1,xx,v1,v2,hmx)
c
      implicit real*8(a-h,o-z)
      real*8 xx(3),v1(10),v2(10)
      integer ip(9)
c
      data ip/1, -1,-1,-1, 1,1,1,1,1/
c
      if(mm0.le.mm1) then
        call nrskc0(jpr,mm0,mm1,xx,v1,hmx)
      else
        call nrskc0(jpr,mm1,mm0,xx,v2,hmx)
        hmx=ip(mm0)*ip(mm1)*hmx
      end if
      if(jpr.ge.3) write(6,100) mm0,mm1,(xx(i),i=1,3),hmx
  100 format(1h ,2i5,3f6.3,f9.3)
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nrtrc0(c0,c1,c2)
c
c#12  nput   :   lmax
c#12  output :   c0, c1, c2     : transformation matrix
c#13  noexternal
c
c#21  to get cubic harmonics from sherical harmonics
c
c#22
c#22  orbital  1: s                                              for l=0
c#22           1: x,  2: y,  3: z                                for l=1
c#22           1: xy, 2: yz, 3: zx, 4: x**2-y**2, 5: 3*z**2-r**2 for l=2
c#22
c
c#31  1991.04.16.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrtrc0(c0,c1,c2)
c
      implicit real*8(a-h,o-z)
      complex*16 c0, c1(-1:1,3), c2(-2:2,5), s2r, s2i
c
      zero=0
      one=1
      two=2
c
      s2=one/sqrt(two)
      s2r=dcmplx(s2,zero)
      s2i=dcmplx(zero,s2)
c
c     l=0
      c0=dcmplx(one,zero)
c
c     l=1
      c1(-1,1)= s2r
      c1( 0,1)= 0
      c1( 1,1)=-s2r
c
      c1(-1,2)= s2i
      c1( 0,2)= 0
      c1( 1,2)= s2i
c
      c1(-1,3)= 0
      c1( 0,3)= 1
      c1( 1,3)= 0
c
c     l=2
      c2(-2,1)= s2i
      c2(-1,1)= 0
      c2( 0,1)= 0
      c2( 1,1)= 0
      c2( 2,1)=-s2i
c                      2<-->3  2001.5.1. N.Hamada
      c2(-2,2)= 0
      c2(-1,2)= s2i
      c2( 0,2)= 0
      c2( 1,2)= s2i
      c2( 2,2)= 0
c
      c2(-2,3)= 0
      c2(-1,3)= s2r
      c2( 0,3)= 0
      c2( 1,3)=-s2r
      c2( 2,3)= 0
c
      c2(-2,4)= s2r
      c2(-1,4)= 0
      c2( 0,4)= 0
      c2( 1,4)= 0
      c2( 2,4)= s2r
c
      c2(-2,5)= 0
      c2(-1,5)= 0
      c2( 0,5)= 1
      c2( 1,5)= 0
      c2( 2,5)= 0
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nrtrc2(lmax,lmlmax,c)
c
c#12  input  :   lmax, lmlmax
c#12  output :   c    : transformation matrix
c#13  noexternal
c
c#21  to get cubic harmonics from sherical harmonics
c
c#22  orbital  1:(0, 0), 2:(1,-1), 3:(1, 0), 4:(1, 1)
c#22           5:(2,-2), 6:(2,-1), 7:(2, 0), 8:(2, 1), 9(2, 2)
c#22
c#22  orbital  1: s                                              for l=0
c#22           2: x,  3: y,  4: z                                for l=1
c#22           5: xy, 6: yz, 7: zx, 8: x**2-y**2, 9: 3*z**2-r**2 for l=2
c#22
c
c#31  1991.04.16.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrtrc2(lmax,lmlmax,c)
c
      implicit real*8(a-h,o-z)
      complex*16 c((lmlmax+1)**2,(lmlmax+1)**2)
      complex*16 c0, c1(-1:1,3), c2(-2:2,5)
c
      zero=0
c
      if(lmax.lt.0 .or. lmax.gt.2 .or. lmax.gt.lmlmax) then
        write(6,*) ' lmax=',lmax,'   lmlmax=',lmlmax
        write(6,*) ' lmax should be less than 3.'
        stop '=nrtrc2(lmax)='
      end if
c
      call nrtrc0(c0,c1,c2)
c
      do 10 j=1,(lmax+1)**2
      do 10 i=1,(lmax+1)**2
   10 c(i,j)=dcmplx(zero,zero)
c
      l=0
      c(1,1)=c0
      if(lmax.eq.0) return
c
      l=1
      do 20 m=-l,l
        lm1=l**2+m+l+1
      do 20 i=1,2*l+1
        lm2=l**2+i
        c(lm1,lm2)=c1(m,i)
   20 continue
      if(lmax.eq.1) return
c
      l=2
      do 30 m=-l,l
        lm1=l**2+m+l+1
      do 30 i=1,2*l+1
        lm2=l**2+i
        c(lm1,lm2)=c2(m,i)
   30 continue
      if(lmax.eq.2) return
c
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrylm1(lmax,lmlmax,theta,phi,ylm)
c
      implicit real*8(a-h,o-z)
      real*8     euler(3)
      complex*16 wd,ylm(-lmlmax:lmlmax,0:lmlmax)
c
      one=1
      pai4=16*datan(one)
      jd=1
      euler(1)=phi
      euler(2)=theta
      euler(3)=0
c
      do 10 l=0,lmax
      do 10 m=-l,l,jd
        call nrrpj0(euler,l,jd,m,0,wd)
        ylm(m,l)=sqrt((2*l+1)/(pai4))*dconjg(wd)
        write(6,200) l,m,ylm(m,l)
   10 continue
  200 format(' l=',i2,'   m=',i2,'   ylm=',2f12.6)
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nrylm2(lmax,theta,phi,ylm)
c
      implicit real*8(a-h,o-z)
      real*8     euler(3)
      complex*16 wd,ylm((lmax+1)**2)
c
      one=1
      pai4=16*datan(one)
      jd=1
      euler(1)=phi
      euler(2)=theta
      euler(3)=0.0
c
      do 10 l=0,lmax
      do 10 m=-l,l
        lm=l**2+m+l+1
        call nrrpj0(euler,l,jd,m,0,wd)
        ylm(lm)=sqrt((2*l+1)/(pai4))*dconjg(wd)
        write(6,200) l,m,ylm(lm)
   10 continue
  200 format(' l=',i2,'   m=',i2,'   ylm=',2f12.6)
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsatpz(ng1,tb1,lrb1,nat,y,kat,lmnat,iatsym)
c
      implicit real*8(a-h,o-z)
      real*8  tb1(3,48),y(3,nat),z(3),u(3)
      integer lrb1(3,3,48),kat(nat),iatsym(lmnat,48)
c
      eps=dfloat(10)**(-5)
      one=dfloat(1)
c
      write(6,*) ' '
      write(6,*) ' symmetry operations on atoms'
      do 10 k=1,ng1
        do 20 i=1,nat
          z(1)=lrb1(1,1,k)*y(1,i)+lrb1(1,2,k)*y(2,i)+
     &         lrb1(1,3,k)*y(3,i)+tb1(1,k)
          z(2)=lrb1(2,1,k)*y(1,i)+lrb1(2,2,k)*y(2,i)+
     &         lrb1(2,3,k)*y(3,i)+tb1(2,k)
          z(3)=lrb1(3,1,k)*y(1,i)+lrb1(3,2,k)*y(2,i)+
     &         lrb1(3,3,k)*y(3,i)+tb1(3,k)
          do 30 j=1,nat
            jj=j
            u(1)=dmod(dabs(z(1)-y(1,j))+eps/2,one)
            u(2)=dmod(dabs(z(2)-y(2,j))+eps/2,one)
            u(3)=dmod(dabs(z(3)-y(3,j))+eps/2,one)
            if(u(1).lt.eps .and. u(2).lt.eps .and. u(3).lt.eps) then
              if(kat(i).eq.kat(j)) then
                go to 34
              else
                go to 32
              end if
            end if
   30     continue
          write(6,*) ' k=',k
          write(6,300) i,(y(l,i),l=1,3),(z(l),l=1,3)
          write(6,*) ' === sub.nsatpz. (no atom) ==='
          stop '=nsatpz(atomic position)='
c
   32     continue
          write(6,*) ' k=',k
          write(6,320) i,(y(l,i),l=1,3),jj,(z(l),l=1,3)
          write(6,*) ' === sub.nsatpz. (different kind) ==='
          stop '=nsatpz(atomic position)='
c
   34     continue
          iatsym(i,k)=jj
   20   continue
        write(6,100) k,(iatsym(i,k),i=1,nat)
   10 continue
      return
  300 format(' iatom=',i5,'  (',3f7.3,')  --->  no atom',
     &                 5x,'  (',3f7.3,')')
  320 format(' iatom=',i5,'  (',3f7.3,')  --->   iatom=',
     &                 i5,'  (',3f7.3,')')
  100 format(' k=',i2,2x,20i3/(7x,20i3))
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsdos0(e,nxx,nyy,nzz,ea,dos,dm0,dm1,
     &     lmnp2,eig2,instts1,np2,iwt,ip2cub)
c
c     ***  by linear interporation in tetrahedrons
c     ***
c     ***
c     ***  nxx  number of mesh points in x-direction
c     ***  nyy  number of mesh points in y-direction
c     ***  nzz  number of mesh points in z-direction
c     ***
c     ***  dos  density of states
c     ***  dm0  0th moment of dos (interated up to e)
c     ***  dm1  1st moment of dos (interated up to e)
c
      implicit none
c
      real*8  ea(*)
      real*8  ecub(2,2,2), ec(8), et(4), eb(4)
      integer iecub(2,2,2),iec(8),iet(4),ieb(4)
      equivalence(ec(1),ecub(1,1,1))
      equivalence(iec(1),iecub(1,1,1))
      integer iqmat(6,2)
      data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
      integer np2,lmnp2,iwt(np2)
      integer instts1
      integer nxx,nyy,nzz,ip2cub(nxx*nyy*nzz)
ccc  commented by K.Mae 2003.8.5
c     real*8 doscub(5000),dm0cub(5000),dm1cub(5000)
ccc  added by K.Mae 2003.8.5
      real*8 doscub(nxx*nyy*nzz)
      real*8 dm0cub(nxx*nyy*nzz)
      real*8 dm1cub(nxx*nyy*nzz)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real*8 eig2(lmnp2)

      real*8 d,d0,d1,dm0,dm1,dos,e,e1,e2,e3,e4
      real*8 emax,emax0,emin,emin0,eps
      integer icub,ip,ip0,iq,it,ix,iy,iz,kx,ky,kz,m
      integer ncub,ni,np,npx,npy,npz,ntet
c  definition of eps  <- must be consistent with <nstts1>
**    eps=dfloat(10)**(-4)
      eps=dfloat(10)**(-4)
c
      npx=nxx+1
      npy=nyy+1
      npz=nzz+1
      np=npx*npy*npz
      ncub=nxx*nyy*nzz
      ntet=6*ncub
c!!! commented by K. Mae 2003.8.5 !!!!
c!      if(ncub.gt.5000) then
c!         write(6,*) ' ncub= ',ncub,' > 5000'
c!         write(6,*) ' : error at sub.nsdos0.'
c!         stop 'error at sub.nsdos0 (ncub).'
c!      end if
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c
      dos=0.d0
      dm0=0.d0
      dm1=0.d0
c
      emax0=-dfloat(10)**30
      emin0= dfloat(10)**30
      do 12 ip=1,np
        if(ea(ip).gt.emax0) emax0=ea(ip)
        if(ea(ip).lt.emin0) emin0=ea(ip)
   12 continue
c                                   ============ if 1 ==
**    if(e.gt.emin0-eps*2) then
      if(e.gt.emax0 .and. instts1.eq.0) then
         do 13 ip = 1, np2
         dm1 = dm1 + eig2(ip) * iwt(ip)
 13      continue
         dm1 = dm1 / ntet / 4.d0
         dm0 = 1.d0
      elseif(e.gt.emin0-eps*2) then
c
c     ***  integration over b.z. starts    ***
c
c     ***       sampling over cubes        ***
c
      icub=0
      do 20 iz=0,nzz-1
      do 20 iy=0,nyy-1
      do 20 ix=0,nxx-1
      icub=icub+1
      if(icub.ne.ip2cub(icub)) then
         dos=dos+doscub(ip2cub(icub))
         dm0=dm0+dm0cub(ip2cub(icub))
         dm1=dm1+dm1cub(ip2cub(icub))
      else
         doscub(icub)=0.d0
         dm0cub(icub)=0.d0
         dm1cub(icub)=0.d0
c     ***  energies at cube corners  ***
        ni=npx*(npy*iz+iy)+ix
        emax=-dfloat(10)**30
        emin= dfloat(10)**30
        do 30 kz=1,2
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
        ecub(kx,ky,kz)=ea(ip0)
        iecub(kx,ky,kz)=ip0
        if(ea(ip0).gt.emax) emax=ea(ip0)
        if(ea(ip0).lt.emin) emin=ea(ip0)
   30   continue
c                                   ============ if 2 ==
        if(e.gt.emin-eps*2) then
**      if(e.gt.emin) then
c         ***      six tetrahedrons      ***
c         *** sampling over tetrahedrons ***
          et(1)=ec(1)
          et(4)=ec(8)
          iet(1)=iec(1)
          iet(4)=iec(8)
          do 40 it=1,6
            do 42 ip=1,2
              iq=iqmat(it,ip)
              et(ip+1)=ec(iq)
              iet(ip+1)=iec(iq)
   42       continue
            do 44 m=1,4
              eb(m)=et(m)
              ieb(m)=iet(m)
   44       continue
c
c           ***  eb(1).le.eb(2).le.eb(3).le.eb(4)  ***
            call nsttod(eb,ieb)
c
            e1=eb(1)
            e2=eb(2)
            e3=eb(3)
            e4=eb(4)
c                                   ============ if 3 ==
c$$$           if(e.ge.e4) then
c$$$               dm0=dm0+1
c$$$               dm1=dm1+(e1+e2+e3+e4)/4.d0
c
c$$$           else if(e.gt.e1) then
c
              call nstts1(e1,e2,e3,e4)
              call nsdosi(e,e1,e2,e3,e4,d,d0,d1)
              dos=dos+d
              dm0=dm0+d0
              dm1=dm1+d1
              doscub(icub)=doscub(icub)+d
              dm0cub(icub)=dm0cub(icub)+d0
              dm1cub(icub)=dm1cub(icub)+d1
c$$$           end if
c                                   ============ if 3 ==
   40   continue
      end if
      endif
c                                   ============ if 2 ==
   20 continue
c
      dos=dos/ntet
      dm0=dm0/ntet
      dm1=dm1/ntet
c
      end if
c                                   ============ if 1 ==

      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsdos2(idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2,
     &     ip20,np0,ea,dos,dm0,dm1,instts1,np2,iwt,ip2cub)
c
      implicit  real*8(a-h,o-z)
c
      real*8 eig2(lmnp2,neig2), ea(np0)
      integer ip20(np0)
      integer iwt(np2),ip2cub(nx*ny*nz)
      integer instts1
c
#ifdef __TIMER_SUB__
      call timer_sta(711)
#endif
      dos=0.d0
      dm0=0.d0
      dm1=0.d0
#ifdef __TIMER_DO__
      call timer_sta(814)
#endif
      do 10 ieig=neig1,neig2
        do 20 k0=1,np0
        ea(k0)=eig2(ip20(k0),ieig)
c       write(6,*) ' k0,k1=',k0,ip20(k0),'   ea=',ea(k0)
   20   continue
        if(idim.eq.1) then
          call lsdos0(e,nz,ea,d,d0,d1)
        else if(idim.eq.2) then
          call msdos0(e,nx,ny,ea,d,d0,d1)
        else if(idim.eq.3) then
          call nsdos0(e,nx,ny,nz,ea,d,d0,d1,
     &          lmnp2,eig2(1,ieig),instts1,np2,iwt,ip2cub)
        else
          write(6,*) ' idim=',idim,' : error'
          stop '=error(idim) sub.nsdos2='
        end if
        dos  =dos  +d
        dm0  =dm0  +d0
        dm1  =dm1  +d1
   10 continue
#ifdef __TIMER_DO__
      call timer_end(814)
#endif
c
#ifdef __TIMER_SUB__
      call timer_end(711)
#endif
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
! ====================== modified by K. Tagami ==================== 11.0
!      subroutine nsdos3(jf,idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2,
!     &                  ip20,np0,ea,instts1,np2,iwt,ip2cub)
      subroutine nsdos3(jf,idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2,
     &                  ip20,np0,ea,instts1,np2,iwt,ip2cub, ipri )
! ================================================================== 11.0
c
      implicit  real*8(a-h,o-z)
c
      real*8 eig2(lmnp2,neig2), ea(np0)
      real*8 dos(100),dm0(100)
      integer ip20(np0)
      integer iwt(np2),ip2cub(nx*ny*nz)
      integer instts1

! ==================== added by K. Tagami ========= 11.0
      integer ipri
! ================================================= 11.0
c
      one=1
      eps=dfloat(10)**(-6)
      iflag=0
      do 10 ieig=neig1,neig2
        do 20 k0=1,np0
        ea(k0)=eig2(ip20(k0),ieig)
c$$$        write(6,*) ' k0,k1=',k0,ip20(k0),'   ea=',ea(k0)
   20   continue
        if(idim.eq.1) then
          call lsdos0(e,nz,ea,d,d0,d1)
        else if(idim.eq.2) then
          call msdos0(e,nx,ny,ea,d,d0,d1)
        else if(idim.eq.3) then
          call nsdos0(e,nx,ny,nz,ea,d,d0,d1,
     &          lmnp2,eig2(1,ieig),instts1,np2,iwt,ip2cub)
        else
          write(6,*) ' idim=',idim,' : error'
          stop '=error(idim) sub.nsdos3='
        end if
        if(iflag.eq.0) then
          if((one-d0).lt.eps) then
            dos(1)=d
            dm0(1)=d0
          else
            dos(2)=d
            dm0(2)=d0
            i0=ieig-2
            j =2
            iflag=1
          end if
        else
          j=j+1
          dos(j)=d
          dm0(j)=d0
          if(d0.lt.eps .or. j.ge.100) go to 30
        end if
   10 continue
c
   30 continue
      if(ipri >= 1) then
         j1=j
         write(jf,*) ' ---        dos for each band, each spin ---'
         write(jf,100) (i0+j,dos(j),j=1,j1)
         write(jf,*) ' --- occupation for each band, each spin ---'
         write(jf,100) (i0+j,dm0(j),j=1,j1)
         write(jf,*) ' -------------------------------------------'
         write(jf,'(" i0 = ",i8)') i0
      end if
  100 format((5(i6,f9.6)))
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsdosi(e,e1,e2,e3,e4,dos,dm0,dm1)
c
c  ** tetrahedron corners at one energy        **  c
c  ** according to lambin and vigneron,        **  c
c  ** phys. rev. b29, 3430 (1984)              **  c
c
      implicit real*8(a-h,o-z)
c     zero=0.d0   change by Tsuyoshi Miyazaki 94.8.28
      zero=1.d-10
      d21=e2-e1
      d31=e3-e1
      d41=e4-e1
      d32=e3-e2
      d42=e4-e2
      d43=e4-e3
      if(d21.le.zero .or. d32.le.zero .or. d43.le.zero) then
        write(6,*) ' Warning!! in = sub.nsdosi ='
        write(6,*) ' e1,e2,e3,e4=',e1,e2,e3,e4
        write(6,*) ' energy order error  in sub.nsdosi.'
c       stop ' = sub.nsdosi ='
      end if
c
      if(e.le.e1) then
        dos=0.d0
        dm0=0.d0
        dm1=0.d0
c
      else if(e.le.e2) then
        d1=e-e1
        dd=d21*d31*d41
        dos=3.d0*(d1**2)/dd
        dm0=  (d1**3)/dd
        dm1=dm0*(3.d0*e+e1)/4.d0
c
      else if(e.lt.e3) then
        d1=e-e1
        d2=e-e2
        d3=e3-e
        d4=e4-e
        dd1=d42*d32*d31
        dd2=d42*d41*d31
        e32=(3.d0*e3-e2)/2.d0-e
        e41=(3.d0*e4-e1)/2.d0-e
        dos=3.d0*(d3*d2/dd1+d4*d1/dd2)
        dm0=(d2**2)*e32/dd1+(d1**2)*e41/dd2-(d21**2)/(2.d0*d31*d42)
        o=(d21**2)*(3.d0*e2+e1)/(4.d0*d31*d41)
        p1=d2*(3.d0*e-4.d0*e3+e2)/4.d0
        p2=e2*(2.d0*e-3.d0*e3+e2)/2.d0
        p=-(d2**2)*(p1+p2)/dd1
        q1=(d2**2)*(-6.d0*d4*d1-2.d0*d2*(2.d0*e-e4-e1)+d2**2)/4.d0
        q2=e2*((d1**2)*(2.d0*e-3.d0*e4+e1)
     &              +(d21**2)*(3.d0*e4-2.d0*e2-e1))/2.d0
        q=-(q1+q2)/dd2
        dm1=o+p+q
      else if(e.lt.e4) then
        d4=e4-e
        dd=d43*d42*d41
        dos=3.d0*(d4**2)/dd
        dm0=1.d0-(d4**3)/dd
        dm1=(e1+e2+e3+e4)/4.d0 - (d4**3)*(3.d0*e+e4)/(4.d0*dd)
      else
        dos =0
        dm0=1.d0
        dm1=(e1+e2+e3+e4)/4.d0
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nseulc(ieuler,euler)
c
c#12  noinput
c#12  output: ieuler(3,48)*2*pai/4: euler angles for cubic
c#12           euler(3,48)        : euler angles (alpha,beta, gamma)
c#13  noexternal
c
c#21  to get euler angles for point-group operations of cubic lattice
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nseulc(ieuler,euler)
c
      implicit real*8(a-h,o-z)
c
      integer  ieuler(3,24)
      dimension euler(3,24)
      integer ie4(3,24)
      data    ie4/
     &  0,0,0, 0,2,2, 0,2,0, 0,0,2,
     &  0,1,1, 2,1,3, 2,1,1, 0,1,3,
     &  1,1,2, 3,1,0, 1,1,0, 3,1,2,
     &  0,2,1, 0,2,3, 0,1,2, 1,1,1, 2,1,0, 3,1,3,
     &  3,1,1, 0,1,0, 1,0,0, 1,1,3, 2,1,2, 3,0,0/
c
      one=1
      pai2=8*datan(one)
      do 10 i=1,24
      do 10 j=1,3
        ieuler(j,i)=ie4(j,i)
        euler(j,i)=(pai2*ie4(j,i))/4
   10 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nseulh(ieuler,euler)
c
c#12  noinput
c#12  output: ieuler(3,12)*2*pai/4: euler angles for hexagonal
c#12           euler(3,12)        : euler angles (alpha,beta, gamma)
c#13  noexternal
c
c#21  to get euler angles for hexagonal operations
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nseulh(ieuler,euler)
c
      implicit real*8(a-h,o-z)
c
      dimension euler(3,12)
      integer ieuler(3,12)
      integer ie6(3,12)
      data    ie6/
     &  0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0,
     &  0,3,0, 4,3,0, 2,3,0, 3,3,0, 1,3,0, 5,3,0/
c
      one=1
      pai2=8*datan(one)
      do 10 i=1,12
      do 10 j=1,3
        ieuler(j,i)=ie6(j,i)
        euler(j,i)=(pai2*dble(ieuler(j,i)))/6
   10 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsgcm2(i1,i2,m)
c
c#12  input   :  i1,i2
c#12  output  :  m     : the greatest common measure
c
c#21  to get the greatest common measure of i1 and i2
c
c#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsgcm2(i1,i2,m)
c
      implicit real*8(a-h,o-z)
c
      if(i1.eq.0 .or. i2.eq.0) then
        write(6,*) '=== error in sub.nsgcm2 ==='
        write(6,*) 'i1,i2=',i1,i2,' : They must not be zero.'
        stop '=== error in sub.nsgcm2 ==='
      else
        i=iabs(i1)
        j=iabs(i2)
    1   continue
        k=mod(j,i)
        if(k.eq.0) go to 2
        j=i
        i=k
        go to 1
    2   m=i
      end if
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsgcm3(i1,i2,i3,m)
c
c#12  input   :  i1,i2,i3
c#12  output  :  m     : the greatest common measure
c#13  external:  sub.nsgcm2
c
c#21  to get the greatest common measure of i1, i2 and i3
c
c#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsgcm3(i1,i2,i3,m)
c
      implicit real*8(a-h,o-z)
      external nsgcm2
c
      call nsgcm2(i1,i2,m1)
      call nsgcm2(i2,i3,m2)
      call nsgcm2(m1,m2,m)
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
c#11 &           schoen,jones,euler,rot,ieuler,irot,ir1234,
c#11 &           im0,iv0,ng0,ng1,ig10,ig01,jg1,
c#11 &           irotr1,irotk1,im1,iv1,movo)
c
c#12  input:      jf : output file
c#12             jpr : print control
c#12              il : lattice type
c#12             inv : parameter (0,1) for moving the origin
c#12         igen(3) : generator (rotation part)
c#12     jgen(2,3,3) : generator (nonprimitive translation part)
c#12  output:  schoen(48),jones(3,48),euler(3,24),rot(3,3,48)
c#12           ieuler(3,24),irot(3,3,48),ir1234(3,48)
c#12           im0(48,48),iv0(48),ng0
c#12           ng1,ig10(48),ig01(48),jg1(2,3,48),
c#12           irotr1(3,3,48),irotk1(3,3,48),im1(48,48),iv1(48)
c#12           mo(2,3): translation of the origin from the initial
c#13  external:  nsrotl, nsmult, nspgrp
c
c#21  to get space group
c
c#31  1990.01.09.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
c ======================= modified by K. Tagami ===================== 12.0A
c      subroutine nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
c     &                  schoen,jones,euler,rot,ieuler,irot,ir1234,
c     &                  im0,iv0,ng0,ng1,ig10,ig01,jg1,
c     &                  irotr1,irotk1,im1,iv1,movo)
      subroutine nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
     &                  schoen,jones,euler,rot,ieuler,irot,ir1234,
     &                  im0,iv0,ng0,ng1,ig10,ig01,jg1,
     &                  irotr1,irotk1,im1,iv1,movo,
     &                  use_trs, tac, tca, tab, tba,
     &                  gen_name_in_carts )
c ==================================================================== 12.0A
c
      implicit real*8(a-h,o-z)
      integer igen(ngen),jgen(2,3,ngen),movo(2,3)

c ========================== added by K. Tagami ====================== 12.0A
      logical use_trs, gen_name_in_carts
      real*8 tac(3,3), tca(3,3), tab(3,3), tba(3,3)
c ==================================================================== 12.0A
c
      character*2 jones(3,48)
      character*5 schoen(48)
      real*8  euler(3,24), rot(3,3,48)
      integer ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &        im0(48,48),iv0(48),ng0
      integer il,ng1,ig10(48),ig01(48),jg1(2,3,48),
     &        irotr1(3,3,48),irotk1(3,3,48),im1(48,48),iv1(48)

c
      call nsrotl(il,jpr,jf,ng0,ieuler,euler,irot,rot,ir1234,jones,
     &            schoen)
c
      call nsmult(jf,jpr,ng0,irot,im0,iv0)
c
c ======================= modified by K. Tagami ===================== 12.0A
c      call nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im0,irot,schoen,
c     &            movo,ng1,ig01,jg1)
c
      call nspgrp_kt(jf,jpr,il,inv,ngen,igen,jgen,im0,irot,schoen,
     &            movo,ng1,ig01,jg1,
     &            use_trs, tac, tca, tab, tba,
     &            gen_name_in_carts )
c ==================================================================== 12.0A
c
      do 10 k=1,ng1
      do 10 j=1,3
      do 10 i=1,3
   10 irotr1(i,j,k)=irot(i,j,ig01(k))
c
      if(il.le.0) then
        do 12 k=1,ng1
        do 12 j=1,3
        do 12 i=1,3
   12   irotk1(i,j,k)=irot(i,j,ig01(k)+24)
      else
        do 14 k=1,ng1
        do 14 j=1,3
        do 14 i=1,3
   14   irotk1(i,j,k)=irot(i,j,iv0(ig01(k)))
      end if
c
      do 20 i=1,48
   20 ig10(i)=0
      do 22 i=1,ng1
   22 ig10(ig01(i))=i
c
      do 24 j=1,ng1
      do 24 i=1,ng1
   24 im1(i,j)=ig10(im0(ig01(i),ig01(j)))
c
      do 26 i=1,ng1
   26 iv1(i)=ig10(iv0(ig01(i)))
c
      if(jpr.ge.0) then
        write(jf,*) ' '
        write(jf,*) ' ----- group elements -------'
        do 50 i=1,ng1
   50  write(jf,100) i,ig01(i),schoen(ig01(i)),(jones(j,ig01(i)),j=1,3),
     &              ((jg1(j,k,i),j=1,2),k=1,3)
  100  format(i5,'   (',i2,')',3x,a5,5x,'(',2(a2,','),a2,')',5x,
     &         '(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )
c
        write(jf,*) ' '
        write(jf,*) ' index ig10(i) :'
        write(jf,120) (i,i= 1,24)
        write(jf,140) (ig10(i),i= 1,24)
        write(jf,*) ' '
        write(jf,120) (i,i=25,48)
        write(jf,140) (ig10(i),i=25,48)
  120 format('    i=',24i3)
  140 format(' ig10=',24i3)
      end if
c
      if(jpr.ge.1) then
        no=ng1/6
        mo=mod(ng1,6)
        write(jf,*) '   '
        write(jf,*) ' matrix representation of operation '
        write(jf,*) ' (for real-space coordinate)'
        do 60 ii=1,no
        write(jf,*) '   '
        i1=(ii-1)*6+1
        i2=ii*6
        write(jf,200) (i,(irotr1(1,k,i),k=1,3),i=i1,i2)
        do 62 j=2,3
   62   write(jf,220) (  (irotr1(j,k,i),k=1,3),i=i1,i2)
   60   continue
c
        if(mo.ne.0) then
        i1=no*6+1
        i2=no*6+mo
        write(jf,*) '   '
        write(jf,200) (i,(irotr1(1,k,i),k=1,3),i=i1,i2)
        do 64 j=2,3
   64   write(jf,220) (  (irotr1(j,k,i),k=1,3),i=i1,i2)
        end if
c
        write(jf,*) '   '
        write(jf,*) ' matrix representation of operation '
        write(jf,*) ' (for reciprocal-space coordinate)'
        do 66 ii=1,no
        write(jf,*) '   '
        i1=(ii-1)*6+1
        i2=ii*6
        write(jf,200) (i,(irotk1(1,k,i),k=1,3),i=i1,i2)
        do 68 j=2,3
   68   write(jf,220) (  (irotk1(j,k,i),k=1,3),i=i1,i2)
   66   continue
c
        if(mo.ne.0) then
        i1=no*6+1
        i2=no*6+mo
        write(jf,*) '   '
        write(jf,200) (i,(irotk1(1,k,i),k=1,3),i=i1,i2)
        do 70 j=2,3
   70   write(jf,220) (  (irotk1(j,k,i),k=1,3),i=i1,i2)
        end if
      end if
  200 format(1h ,6('(',i2,')',3i2,2x))
  220 format(1h ,6(4x        ,3i2,2x))
c
      if(jpr.ge.1) then
        write(jf,*) '   '
        write(jf,*) '--- group multiplication table ---'
        if(ng1.le.24) then
          write(jf,320) (j,j=1,ng1)
          write(jf,340) ('---',j=1,ng1)
          do 80 i=1,ng1
   80     write(jf,300) i,(im1(i,j),j=1,ng1)
        else
          write(jf,*) ' '
          write(jf,320) (j,j=1,24)
          write(jf,340) ('---',j=1,24)
          do 82 i=1,ng1
   82     write(jf,300) i,(im1(i,j),j=1,24)
          write(jf,*) ' '
          write(jf,320) (j,j=25,ng1)
          write(jf,340) ('---',j=25,ng1)
          do 84 i=1,ng1
   84     write(jf,300) i,(im1(i,j),j=25,ng1)
  300     format((i3,2x,24i3))
  320     format((5x,24i3))
  340     format(5x,24a3)
        end if
        write(jf,*) ' '
        write(jf,*) '--- invers elements ---'
        write(jf,400) (iv1(j),j=1,ng1)
  400   format((5x,24i3))
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsgrpb(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
c
c#12  input:     ng1 : order of group
c#12        tab(3,3) : transformation matix
c#12        tba(3,3) : transformation matix
c#12             ta1 : nonprimitive translation vector (A system)
c#12             ra1 : rotation matrix in real space (A system)
c#12             sa1 : rotation matrix in reciprocal space (A system)
c#12                   ra1*sa1=unit matrix
c#12  output:    tb1 : nonprimitive translation vector (B system)
c#12             rb1 : rotation matrix in real space (B system)
c#12             sb1 : rotation matrix in reciprocal space (B system)
c#12                   rb1*sb1=unit matrix
c
c#21  to get space-group matrix notation in B system
c
c#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsgrpb(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
c
      implicit real*8(a-h,o-z)
      real*8 tab(3,3),tba(3,3)
      real*8 ta1(3,48),ra1(3,3,48),sa1(3,3,48)
      real*8 tb1(3,48),rb1(3,3,48),sb1(3,3,48)
c
      do 20 k=1,ng1
        tb1(1,k)=tba(1,1)*ta1(1,k)+tba(1,2)*ta1(2,k)+tba(1,3)*ta1(3,k)
        tb1(2,k)=tba(2,1)*ta1(1,k)+tba(2,2)*ta1(2,k)+tba(2,3)*ta1(3,k)
        tb1(3,k)=tba(3,1)*ta1(1,k)+tba(3,2)*ta1(2,k)+tba(3,3)*ta1(3,k)
   20 continue
c
      do 22 k=1,ng1
        do 24 i=1,3
        do 24 j=1,3
          r=0
          s=0
          do 26 l=1,3
          do 26 m=1,3
            r=r+tba(i,l)*ra1(l,m,k)*tab(m,j)
            s=s+tba(i,l)*sa1(l,m,k)*tab(m,j)
   26     continue
          rb1(i,j,k)=r
          sb1(i,j,k)=s
   24  continue
   22 continue
c
      return
      end

c ================================= added by K. Tagami =============== 12.0A
      subroutine nsgrpb_kt(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
      implicit none
c inout
      integer ng1
      real*8 tab(3,3),tba(3,3)
      real*8 ta1(3,48),ra1(3,3,48),sa1(3,3,48)
      real*8 tb1(3,48),rb1(3,3,48),sb1(3,3,48)
c local
      integer i, j, k, l, m
      real*8 r, s
c begin
c      do k=1,ng1
c         tb1(1,k)=tba(1,1)*ta1(1,k)+tba(1,2)*ta1(2,k)+tba(1,3)*ta1(3,k)
c         tb1(2,k)=tba(2,1)*ta1(1,k)+tba(2,2)*ta1(2,k)+tba(2,3)*ta1(3,k)
c         tb1(3,k)=tba(3,1)*ta1(1,k)+tba(3,2)*ta1(2,k)+tba(3,3)*ta1(3,k)
c      end do
      do k=1,ng1
         tb1(1,k)=tab(1,1)*ta1(1,k)+tab(1,2)*ta1(2,k)+tab(1,3)*ta1(3,k)
         tb1(2,k)=tab(2,1)*ta1(1,k)+tab(2,2)*ta1(2,k)+tab(2,3)*ta1(3,k)
         tb1(3,k)=tab(3,1)*ta1(1,k)+tab(3,2)*ta1(2,k)+tab(3,3)*ta1(3,k)
      end do
c
      do k=1,ng1
         do i=1,3
            do j=1,3
               r=0
               s=0
               do l=1,3
                  do m=1,3
                     r = r +tab(i,l) *ra1(l,m,k) *tba(m,j)
                     s = s +tab(i,l) *sa1(l,m,k) *tba(m,j)
                  end do
               end do
               rb1(i,j,k)=r
               sb1(i,j,k)=s
            end do
         end do
      end do
      return
      end
c ==================================================================== 12.0A

c ===================== added by K. Tagami =========================== 12.0A
      subroutine nsgrpa( ng1,tac,tca, rc1, sc1, ra1, sa1 )
      implicit none
c inout
      integer ng1
      real*8 tac(3,3),tca(3,3)
      real*8 ra1(3,3,48),sa1(3,3,48)
      real*8 rc1(3,3,48),sc1(3,3,48)
c local
      integer i, j, k, m1, m2
      real*8 c1, c2
c begin
c ---------------------------------------------
c Note :
c          ra1 is similar to m_CS_op_in_PUCD
c                             in  m_Crystal_Structure.F90
c
c ----------------------------------------------
      do k=1,ng1
         do i=1,3
            do j=1,3
               c1 = 0.0d0
               c2 = 0.0d0

               Do m1=1, 3
                  Do m2=1, 3
                     c1 = c1 + tac(i,m1) *rc1(m1,m2,k) *tca(m2,j)
                     c2 = c2 + tac(i,m1) *sc1(m1,m2,k) *tca(m2,j)
                  End do
               End do
               ra1(i,j,k) = c1
               sa1(i,j,k) = c2
            End do
         End do
      End do
      return
      end
c ==================================================================== 12.0A



c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsgrpw(jf,schoen,ng1,ig01,t)
c
c#12  input:           jf : output file
c#12           schoen(48) : Schoenflies notation
c#12                 ng1  : order of the group
c#12                ig01  : list vector to get the element code
c#12                   t  : nonprimitive translation
c#12  nooutput
c#13  noexternal
c
c#21  to write group elements
c
c#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsgrpw(jf,schoen,ng1,ig01,t)
c
      implicit real*8(a-h,o-z)
      character*5 schoen(48)
      real*8  t(3,48)
      integer ng1,ig01(48)
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsig10(ng1,ig01,ig10)
c
      implicit real*8(a-h,o-z)
      dimension ig01(48),ig10(48)
c
      do 10 ig0=1,48
   10 ig10(ig0)=99
      do 20 ig1=1,ng1
   20 ig10(ig01(ig1))=ig1
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsjonh(irot,ir1234,jones)
c
c#12  input : irot(3,3,24)    (i) : rotation matrix in integer
c#12  output: irot(3,3,24)    (i) : rotation matrix in integer
c#12          ir1234(3,24)    (i) : jones faithful representation
c#12           jones(3,24)    (a2): jones faithful representation
c#13  noexternal
c
c#21  to get rotation matrix, and jones faithfull representation
c#21         for point-group operations of hexagonal lattice
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsjonh(irot,ir1234,jones)
c
      implicit real*8(a-h,o-z)
      character*2 jones(3,24)
      integer irot(3,3,24),ir1234(3,24)
c
c
      do 20 i=1,24
        do 24 j=1,2
          k1=irot(j,1,i)
          k2=irot(j,2,i)
          k3=irot(j,3,i)
          if(k3.ne.0) then
            write(6,9000) i,j,(irot(j,k,i),k=1,3)
            stop ' === stop in sub.nsjonh. (invalid r, 1)'
          end if
          if(k1.eq.1) then
            if(k2.eq.1) then
              ir1234(j,i)=4
              jones(j,i)=' w'
c             --- for reciprocal space ---
            else if(k2.eq.0) then
              ir1234(j,i)=1
              jones(j,i)=' x'
            else if(k2.eq.-1) then
              ir1234(j,i)=4
              jones(j,i)=' w'
            else
              write(6,9000) i,j,(irot(j,k,i),k=1,3)
              stop ' === stop in sub.nsjonh. (invalid r, 1)'
            end if
          else if(k1.eq.0) then
            if(k2.eq.1) then
              ir1234(j,i)=2
              jones(j,i)=' y'
            else if(k2.eq.-1) then
              ir1234(j,i)=-2
              jones(j,i)='-y'
            else
              write(6,9000) i,j,(irot(j,k,i),k=1,3)
              stop ' === stop in sub.nsjonh. (invalid r, 3) ==='
            end if
          else if(k1.eq.-1) then
            if(k2.eq.1) then
              ir1234(j,i)=-4
              jones(j,i)='-w'
            else if(k2.eq.0) then
              ir1234(j,i)=-1
              jones(j,i)='-x'
            else if(k2.eq.-1) then
              ir1234(j,i)=-4
              jones(j,i)='-w'
c             --- for reciprocal space ---
            else
              write(6,9000) i,j,(irot(j,k,i),k=1,3)
              stop ' === stop in sub.nsjonh. (invalid r, 4) ==='
            end if
          end if
   24   continue
c
        j=3
        if(irot(3,1,i).ne.0 .or. irot(3,2,i).ne.0) then
          write(6,9000) i,j,(irot(j,k,i),k=1,3)
          stop ' === stop in sub.nsjonh. (invalid r, 5) ==='
        else if(irot(3,3,i).eq.1)  then
          ir1234(j,i)=3
          jones(j,i)=' z'
        else if(irot(3,3,i).eq.-1) then
          ir1234(j,i)=-3
          jones(j,i)='-z'
        else
          write(6,9000) i,j,(irot(j,k,i),k=1,3)
          stop ' === stop in sub.nsjonh. (invalid r, 6) ==='
        end if
   20 continue
c
 9000 format(1h ,'operation=',i2,'  j=',i1,'   (r(j,k),k=1,3)',3i3)
      return
      end
c     implicit real*8(a-h,o-z)
c     parameter(lmnkg=10000,lmne=10000)
c     real*8  gkb(3,3)
c     real*8  ekg(lmnkg)
c     integer kg(3,lmnkg),jkge(lmne)
c     data gkb /1.d0, 0.d0, 0.d0,
c    &          0.d0, 1.d0, 0.d0,
c    &          0.d0, 0.d0, 1.d0/
c     open(unit= 6, file='a.out6')
c     jpr=3
c     emax=5.1d0/2
c     call nskg00(jpr,emax,gkb, nkg,lmnkg,kg,ekg, ne,lmne,jkge)
c     stop
c     end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nskg00(jpr,emax,gkb, nkg,lmnkg,kg,ekg, ne,lmne,jkge)
c
c#12  input:     jpr     : print control
c#12             emax    : maximum kinetic energy (Hartree unit)
c#12             gkb(3,3): metric tensor (atomic unit in k-space)
c#12             lmnkg   : limit of nkg
c#12             lmne    : limit of ne
c#12
c#12  output:    nkg     : # of g vectors
c#12             kg(3,lmnkg): g vectors
c#12             ekg( lmnkg): kinetic energy of g vector (Hartree unit)
c#12             ne         : # of energy shells
c#12             jkge(lmne) : total # of vectors inside the shell
c#12  noexternal:
c
c#21  to get g vectors (reciprocal lattice)
c
c#31  1993.04.21.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskg00(jpr,emax,gkb, nkg,lmnkg,kg,ekg,ne,lmne,jkge)
c
      implicit real*8(a-h,o-z)
      real*8  gkb(3,3)
      real*8  ekg(lmnkg)
      integer kg(3,lmnkg), jkge(lmne)
c
      eps=dfloat(10)**(-8)
c
      rm=sqrt(emax*2)
c
      nx=0
      ny=0
      nz=0
      do 30 iz=-1,1
      do 32 iy=-1,1
      do 34 ix=-1,1
      if(ix.ne.0 .or. iy.ne.0 .or. iz.ne.0) then
        p1=ix
        p2=iy
        p3=iz
        rr=sqrt(p1*gkb(1,1)*p1+p1*gkb(1,2)*p2+p1*gkb(1,3)*p3
     &         +p2*gkb(2,1)*p1+p2*gkb(2,2)*p2+p2*gkb(2,3)*p3
     &         +p3*gkb(3,1)*p1+p3*gkb(3,2)*p2+p3*gkb(3,3)*p3)
        nn=(rm/rr)+5
        nnx=iabs(ix)*nn
        nny=iabs(iy)*nn
        nnz=iabs(iz)*nn
        if(nnx.gt.nx) nx=nnx
        if(nny.gt.ny) ny=nny
        if(nnz.gt.nz) nz=nnz
      endif
   34 continue
   32 continue
   30 continue
c
c     write(6,'('' nx,ny,nz='',3i5)') nx,ny,nz
c
c
      k=0
      do 10 iz=-nz,nz
      do 12 iy=-ny,ny
      do 14 ix=-nx,nx
        ee=(ix*gkb(1,1)*ix+ix*gkb(1,2)*iy+ix*gkb(1,3)*iz
     &     +iy*gkb(2,1)*ix+iy*gkb(2,2)*iy+iy*gkb(2,3)*iz
     &     +iz*gkb(3,1)*ix+iz*gkb(3,2)*iy+iz*gkb(3,3)*iz)/2
c
c                                                registration
        if(ee.le.emax) then
          k=k+1
          if(k.gt.lmnkg) then
            write(6,*) ' k=',k,'  > lmnkg=',lmnkg,' in sub.nskg00'
            stop 'error in sub.nskg00 (lmnkg)'
          end if
          ekg(k)=ee
          kg(1,k)=ix
          kg(2,k)=iy
          kg(3,k)=iz
        end if
   14 continue
   12 continue
   10 continue
      nkg=k
c
c                                                sorting
      call nskgs0(nkg,ekg,kg)
c

      m=0
      do 20 k=2,nkg
        if(abs(ekg(k-1)-ekg(k)).gt.eps) then
          m=m+1
          if(m.gt.lmne) then
            write(6,*) ' m=',m,' > lmne=',lmne,'  in sub.nskg00'
            stop 'error in sub.nskg00 (lmne)'
          end if
          jkge(m)=k-1
        end if
   20 continue
      ne=m+1
      jkge(ne)=nkg
c
      if(jpr.ge.3) then
        write(6,900) nkg
        do 90  k=1,nkg
          write(6,920) k,(kg(i,k),i=1,3),ekg(k),ekg(k)*2
   90   continue
      end if
c
      if(jpr.ge.1) then
        write(6,940) ne
        write(6,900) nkg
        do 92  m=1,ne
          k=jkge(m)
          write(6,920) k,(kg(i,k),i=1,3),ekg(k),ekg(k)*2
   92   continue
      end if
c
  900 format(i10,' g vectors', 15x,'energy(Hartree)',5x,' (Rydberg)')
  920 format(i10,'  (',3i5,' )',5x,2f15.6)
  940 format(i10,' shells of g vectors were produced.')
      return
      end
c     implicit real*8(a-h,o-z)
c     parameter(lmnpg=10000,lmne=10000)
c     real*8  gkb(3,3), p(3)
c     real*8  epg(lmnpg)
c     real*8  pg(3,lmnpg)
c     integer jpge(lmne)
c     data gkb /1.d0, 0.d0, 0.d0,
c    &          0.d0, 1.d0, 0.d0,
c    &          0.d0, 0.d0, 1.d0/
c     data p   /0.5d0, 0.5d0, 0.5d0/
c     open(unit= 6, file='a.out6')
c     jpr=3
c     emax=5.1d0/2
c     call nskg01(jpr,emax,gkb,p, npg,lmnpg,pg,epg, ne,lmne,jpge)
c     stop
c     end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nskg01(jpr,emax,gkb,p, npg,lmnpg,pg,epg, ne,lmne,jpge)
c
c#12  input:     jpr     : print control
c#12             emax    : maximum kinetic energy (Hartree unit)
c#12             gkb(3,3): metric tensor (atomic unit in k-space)
c#12             p(3)    : k vector
c#12             lmnpg   : limit of npg
c#12             lmne    : limit of ne
c#12
c#12  output:    npg     : # of k+g vectors
c#12             pg(3,lmnpg): k+g vectors
c#12             epg( lmnpg): kinetic energy of p+g vector (Hartree unit)
c#12             ne         : # of energy shells
c#12             jpge(lmne) : total # of vectors inside the shell
c#12  noexternal:
c
c#21  to get k+g vectors (reciprocal lattice)
c
c#31  1993.04.26.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskg01(jpr,emax,gkb,p0,npg,lmnpg,pg,epg,ne,lmne,jpge)
c
      implicit real*8(a-h,o-z)
      real*8  gkb(3,3), p0(3), p(3)
      real*8  epg(lmnpg), pg(3,lmnpg)
      integer jpge(lmne)
c
      eps=dfloat(10)**(-8)
      one=1
c
      rm=sqrt(emax*2)
c
      p(1)=dmod(p0(1),one)
      p(2)=dmod(p0(2),one)
      p(3)=dmod(p0(3),one)
c
      nx=0
      ny=0
      nz=0
      do 30 iz=-1,1
      do 32 iy=-1,1
      do 34 ix=-1,1
      if(ix.ne.0 .or. iy.ne.0 .or. iz.ne.0) then
        p1=ix
        p2=iy
        p3=iz
        rr=sqrt(p1*gkb(1,1)*p1+p1*gkb(1,2)*p2+p1*gkb(1,3)*p3
     &         +p2*gkb(2,1)*p1+p2*gkb(2,2)*p2+p2*gkb(2,3)*p3
     &         +p3*gkb(3,1)*p1+p3*gkb(3,2)*p2+p3*gkb(3,3)*p3)
        nn=(rm/rr)+5
        nnx=iabs(ix)*nn
        nny=iabs(iy)*nn
        nnz=iabs(iz)*nn
        if(nnx.gt.nx) nx=nnx
        if(nny.gt.ny) ny=nny
        if(nnz.gt.nz) nz=nnz
      endif
   34 continue
   32 continue
   30 continue
c
c     write(6,'('' nx,ny,nz='',3i5)') nx,ny,nz
c
      k=0
      nnx=0
      nny=0
      nnz=0
      do 10 iz=-nz,nz
      do 12 iy=-ny,ny
      do 14 ix=-nx,nx
        p1=p(1)+ix
        p2=p(2)+iy
        p3=p(3)+iz
        ee=(p1*gkb(1,1)*p1+p1*gkb(1,2)*p2+p1*gkb(1,3)*p3
     &     +p2*gkb(2,1)*p1+p2*gkb(2,2)*p2+p2*gkb(2,3)*p3
     &     +p3*gkb(3,1)*p1+p3*gkb(3,2)*p2+p3*gkb(3,3)*p3)/2
c
c                                                registration
        if(ee.le.emax) then
          k=k+1
          if(k.gt.lmnpg) then
            write(6,*) ' k=',k,'  > lmnpg=',lmnpg,' in sub.nskg01'
            stop 'error in sub.nskg01 (lmnpg)'
          end if
          epg(k)=ee
          pg(1,k)=p1
          pg(2,k)=p2
          pg(3,k)=p3
          if(ix.gt.nnx) nnx=ix
          if(iy.gt.nny) nny=iy
          if(iz.gt.nnz) nnz=iz
        end if
   14 continue
   12 continue
   10 continue
      npg=k
c
c     write(6,'('' nnx,nny,nnz='',3i5)') nnx,nny,nnz
c
c                                                sorting
      call nskgs1(npg,epg,pg)
c

      m=0
      do 20 k=2,npg
        if(abs(epg(k-1)-epg(k)).gt.eps) then
          m=m+1
          if(m.gt.lmne) then
            write(6,*) ' m=',m,' > lmne=',lmne,'  in sub.nskg01'
            stop 'error in sub.nskg01 (lmne)'
          end if
          jpge(m)=k-1
        end if
   20 continue
      ne=m+1
      jpge(ne)=npg
c
      if(jpr.ge.1) then
        write(6,940) ne
        write(6,900) npg
      end if
c
      if(jpr.ge.3) then
        write(6,900) npg
        do 90  k=1,npg
          write(6,920) k,(pg(i,k),i=1,3),epg(k),epg(k)*2
   90   continue
      end if
c
      if(jpr.ge.2) then
        write(6,940) ne
        write(6,900) npg
        do 92  m=1,ne
          k=jpge(m)
          write(6,920) k,(pg(i,k),i=1,3),epg(k),epg(k)*2
   92   continue
      end if
c
  900 format(i10,' g vectors', 24x,'energy(Hartree)',5x,' (Rydberg)')
  920 format(i10,'  (',3f8.3,' )',5x,2f15.6)
  940 format(i10,' shells of g vectors were produced.')
      return
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsnbs0(n,ra,ka)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  ka(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskgs0(n,ra,ka)
c
      implicit real*8(a-h,o-z)
      real*8  ra(n)
      integer ka(3,n)
c
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          k1a=ka(1,l)
          k2a=ka(2,l)
          k3a=ka(3,l)
        else
          rra=ra(ir)
          k1a=ka(1,ir)
          k2a=ka(2,ir)
          k3a=ka(3,ir)
          ra(ir)=ra(1)
          ka(1,ir)=ka(1,1)
          ka(2,ir)=ka(2,1)
          ka(3,ir)=ka(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            ka(1,1)=k1a
            ka(2,1)=k2a
            ka(3,1)=k3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            ka(1,i)=ka(1,j)
            ka(2,i)=ka(2,j)
            ka(3,i)=ka(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        ka(1,i)=k1a
        ka(2,i)=k2a
        ka(3,i)=k3a
c
      go to 10
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nskgs1(n,ra,pa)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  pa(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskgs1(n,ra,pa)
c
      implicit real*8(a-h,o-z)
      real*8  ra(n), pa(3,n)
c
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          p1a=pa(1,l)
          p2a=pa(2,l)
          p3a=pa(3,l)
        else
          rra=ra(ir)
          p1a=pa(1,ir)
          p2a=pa(2,ir)
          p3a=pa(3,ir)
          ra(ir)=ra(1)
          pa(1,ir)=pa(1,1)
          pa(2,ir)=pa(2,1)
          pa(3,ir)=pa(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            pa(1,1)=p1a
            pa(2,1)=p2a
            pa(3,1)=p3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            pa(1,i)=pa(1,j)
            pa(2,i)=pa(2,j)
            pa(3,i)=pa(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        pa(1,i)=p1a
        pa(2,i)=p2a
        pa(3,i)=p3a
c
      go to 10
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
c
c     (nxx+1)*(nyy+1)*(nzz+1) :  no. of mesh points
c     nxx*nyy*nzz :  no. of cubes
c     nx1,ny1,nz1 :  for constructing k point
c     nd :  denominator for constructing k point
c
c
      implicit real*8(a-h,o-z)
      data id/1/
c
      if(il.lt.-1 .or. il.gt.4) then
c        write(6,*) ' il=',il
        stop ' === stop in sub.nskma0. (il) ==='
      end if
c
      if(id.eq.1) then
        if(il.eq.-1) then
          nxx=nx
          nyy=ny
          nzz=nz*3
        else if(il.le.1) then
          nxx=nx
          nyy=ny
          nzz=nz
        else if(il.eq.2) then
          nxx=nx
          nyy=ny*2
          nzz=nz*2
        else if(il.eq.3) then
          nxx=nx
          nyy=ny
          nzz=nz*2
        else if(il.eq.4) then
          nxx=nx
          nyy=ny*2
          nzz=nz
        else if(il.eq.5) then
          nxx=nx
          nyy=ny
          nzz=nz*2
        else if(il.eq.6) then
          nxx=nx
          nyy=ny
          nzz=nz*2
        end if
      else if(id.eq.2) then
        if(il.eq.-1) then
          nxx=nx*3
          nyy=ny*3
          nzz=nz*3
        else if(il.le.1) then
          nxx=nx
          nyy=ny
          nzz=nz
        else if(il.le.3) then
          nxx=nx*2
          nyy=ny*2
          nzz=nz*2
        else if(il.eq.4) then
          nxx=nx*2
          nyy=ny*2
          nzz=nz
        else if(il.eq.5) then
          nxx=nx
          nyy=ny*2
          nzz=nz*2
        else if(il.eq.6) then
          nxx=nx*2
          nyy=ny
          nzz=nz*2
        end if
      else
c        write(6,*) ' id=',id
        stop ' === stop in sub.nskma0. (id) ==='
      end if
c
      if(nx.eq.0) then
         nx1=1
      else
         nx1=nx
      end if
      if(ny.eq.0) then
         ny1=1
      else
         ny1=ny
      end if
      if(nz.eq.0) then
         nz1=1
      else
         nz1=nz
      end if
      nd=nx1*ny1*nz1
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp, np,p)
c
      implicit real*8(a-h,o-z)
c
      real*8 p(3,lmnp)
c
      np=0
      do 10 iz=0,nzz
      do 10 iy=0,nyy
      do 10 ix=0,nxx
        np=np+1
        if(np.gt.lmnp) then
          write(6,*) 'ix,iy,iz=',ix,iy,iz,'   np,lmnp=',np,lmnp
          write(6,*) 'nxx,nyy,nzz=',nxx,nyy,nzz
          stop ' === stop in sub.nskp00. (np>lmnp) ==='
        end if
        p(1,np)=ix *ny1*nz1/dfloat(nd)
        p(2,np)=nx1*iy *nz1/dfloat(nd)
        p(3,np)=nx1*ny1*iz /dfloat(nd)
   10 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
c#11 &           pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
c#11 &           iu21,iv21)
c
c#12  input:
c#12
c#12
c#12
c#12  output:
c#12
c#12
c#13  external:  nspbge
c
c#21  to get coordination systems (C,A,B)
c
c#31  1990.11.20.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &                  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &                  iu21,iv21, itrs )
c
      implicit real*8(a-h,o-z)
c
      real*8  tab(3,3)
      integer lsb1(3,3,ng1),iv1(ng1)
c
      real*8  pa0(3,lmnp0),pb0(3,lmnp0),p(3)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer, allocatable, dimension(:,:) :: pb05i
      integer :: pi(3)
c
      do 10 i=1,np0
        pb0(1,i)=pa0(1,i)*tab(1,1)+pa0(2,i)*tab(2,1)+pa0(3,i)*tab(3,1)
        pb0(2,i)=pa0(1,i)*tab(1,2)+pa0(2,i)*tab(2,2)+pa0(3,i)*tab(3,2)
        pb0(3,i)=pa0(1,i)*tab(1,3)+pa0(2,i)*tab(2,3)+pa0(3,i)*tab(3,3)
   10 continue
      allocate(pb05i(3,lmnp0))
      pb05i(:,:) = nint(100000*pb0(:,:))
c
      np1=1
      ip10(1)=1
      ip01(1)=1
      do 20 i=2,np0
        do 22 j=1,np1
          jj=j
!          call nspbge(pb0(1,i),pb0(1,ip01(j)),ind)
          call nspbgei(pb05i(1,i),pb05i(1,ip01(j)),ind)
          if(ind.eq.0) go to 24
   22   continue
        np1=np1+1
        if(np1.gt.lmnp1) then
          write(6,*) ' np1=',np1,' > lmnp1=',lmnp1
          stop '=nskpb0(np1)='
        end if
        ip10(i)=np1
        ip01(np1)=i
        go to 20
   24   continue
        ip10(i)=jj
   20 continue
c
      np2=1
      ip21(1)=1
      ip12(1)=1
      ip02(1)=1
      iu21(1)=1
      iv21(1)=1
      do 30 i=2,np1
        do k=1,ng1
           do m=0, itrs
              kk=k
              p(1)=pb0(1,ip01(i))*lsb1(1,1,k)
     &             +pb0(2,ip01(i))*lsb1(2,1,k)
     &             +pb0(3,ip01(i))*lsb1(3,1,k)
              p(2)=pb0(1,ip01(i))*lsb1(1,2,k)
     &             +pb0(2,ip01(i))*lsb1(2,2,k)
     &             +pb0(3,ip01(i))*lsb1(3,2,k)
              p(3)=pb0(1,ip01(i))*lsb1(1,3,k)
     &             +pb0(2,ip01(i))*lsb1(2,3,k)
     &             +pb0(3,ip01(i))*lsb1(3,3,k)
              p(1) = p(1) *( 1 -2*m )
              p(2) = p(2) *( 1 -2*m )
              p(3) = p(3) *( 1 -2*m )
              pi = nint(100000*p)

              do j=1,np2
                 jj=j
!                 call nspbge(p,pb0(1,ip02(j)),ind)
                 call nspbgei(pi,pb05i(1,ip02(j)),ind)
                 if(ind.eq.0) go to 34
              end do
           end do
        end do
 32     continue

        np2=np2+1
        if(np2.gt.lmnp2) then
          write(6,*) ' np2=',np2,' > lmnp2=',lmnp2
          stop '=nskpb0(np2)='
        end if
        ip21(i)=np2
        ip12(np2)=i
        ip02(np2)=ip01(i)
        iu21(i)=1
        iv21(i)=1
        go to 30
   34   continue
        ip21(i)=jj
        iu21(i)=kk
        iv21(i)=iv1(kk)
   30 continue
c
      do 40 i=1,np0
   40 ip20(i)=ip21(ip10(i))
c
      if(jpr.ge.3) then
         write(6,*) ' === k points === '
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
 120     format((1h ,2('     (',3f9.5,' ) ')))
      end if
      if(jpr.ge.2) then
         write(6,*) '=== k-point code === '
         write(6,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
 140     format((1h ,2(i4,' (',3f9.5,' ) ')))
c
      write(6,*) '----------------------------------------------------'
         write(6,180)
 180     format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
         write(6,200) (ip0,ip10(ip0),iu21(ip10(ip0)),ip21(ip10(ip0)),
     &        ip0=1,np0)
 200     format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
      write(6,*) '----------------------------------------------------'
      end if
      return
      end

      subroutine nskpb0_s(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,
     &                  lmnp2,pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &                  iu21,iv21, itrs )
      use m_Timing, only : tstatc0_begin, tstatc0_end
      use m_spg_tetra

c
      implicit real*8(a-h,o-z)
      include 'mpif.h'
c
      real*8  tab(3,3)
      integer lsb1(3,3,ng1),iv1(ng1)
c
      real*8  pa0(3,lmnp0),pb0(3,lmnp0),p(3)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer, allocatable, dimension(:,:) :: pb05i
      integer :: pi(3)
      logical :: exi
      integer :: id_sname = -1
c
      if (nspb0_vars_ready(lmnp0,lmnp1,lmnp2)) then
        call get_nskpb0_vars(np1, np2, pb0,
     &       ip10, ip20, ip01,ip21,ip02,ip12,iu21,iv21)
        return
      endif

      call tstatc0_begin('nskpb0_s ',id_sname)
      do 10 i=1,np0
        pb0(1,i)=pa0(1,i)*tab(1,1)+pa0(2,i)*tab(2,1)+pa0(3,i)*tab(3,1)
        pb0(2,i)=pa0(1,i)*tab(1,2)+pa0(2,i)*tab(2,2)+pa0(3,i)*tab(3,2)
        pb0(3,i)=pa0(1,i)*tab(1,3)+pa0(2,i)*tab(2,3)+pa0(3,i)*tab(3,3)
   10 continue
      allocate(pb05i(3,lmnp0))
      pb05i(:,:) = nint(100000*pb0(:,:))
c
      np1=1
      ip10(1)=1
      ip01(1)=1
      do 20 i=2,np0
        do 22 j=1,np1
          jj=j
          call nspbgei(pb05i(1,i),pb05i(1,ip01(j)),ind)
          if(ind.eq.0) go to 24
   22   continue
        np1=np1+1
        if(np1.gt.lmnp1) then
          write(6,*) ' np1=',np1,' > lmnp1=',lmnp1
          stop '=nskpb0(np1)='
        end if
        ip10(i)=np1
        ip01(np1)=i
        go to 20
   24   continue
        ip10(i)=jj
   20 continue
c
      np2=1
      ip21(1)=1
      ip12(1)=1
      ip02(1)=1
      iu21(1)=1
      iv21(1)=1
      do 30 i=2,np1
        do k=1,ng1
           do m=0, itrs
              kk=k
              p(1)=pb0(1,ip01(i))*lsb1(1,1,k)
     &             +pb0(2,ip01(i))*lsb1(2,1,k)
     &             +pb0(3,ip01(i))*lsb1(3,1,k)
              p(2)=pb0(1,ip01(i))*lsb1(1,2,k)
     &             +pb0(2,ip01(i))*lsb1(2,2,k)
     &             +pb0(3,ip01(i))*lsb1(3,2,k)
              p(3)=pb0(1,ip01(i))*lsb1(1,3,k)
     &             +pb0(2,ip01(i))*lsb1(2,3,k)
     &             +pb0(3,ip01(i))*lsb1(3,3,k)
              p(1) = p(1) *( 1 -2*m )
              p(2) = p(2) *( 1 -2*m )
              p(3) = p(3) *( 1 -2*m )

              pi = nint(100000*p)
              do j=1,np2
                 jj=j
                 call nspbgei(pi,pb05i(1,ip02(j)),ind)
                 if(ind.eq.0) go to 34
              end do
           end do
        end do
 32     continue

        np2=np2+1
        if(np2.gt.lmnp2) then
          write(6,*) ' np2=',np2,' > lmnp2=',lmnp2
          stop '=nskpb0(np2)='
        end if
        ip21(i)=np2
        ip12(np2)=i
        ip02(np2)=ip01(i)
        iu21(i)=1
        iv21(i)=1
        go to 30
   34   continue
        ip21(i)=jj
        iu21(i)=kk
        iv21(i)=iv1(kk)
   30 continue
c
      do 40 i=1,np0
   40 ip20(i)=ip21(ip10(i))
c
      if(jpr.ge.3) then
         write(6,*) ' === k points === '
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
 120     format((1h ,2('     (',3f9.5,' ) ')))
      end if
      if(jpr.ge.2) then
         write(6,*) '=== k-point code === '
         write(6,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
 140     format((1h ,2(i4,' (',3f9.5,' ) ')))
c
      write(6,*) '----------------------------------------------------'
         write(6,180)
 180     format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
         write(6,200) (ip0,ip10(ip0),iu21(ip10(ip0)),ip21(ip10(ip0)),
     &        ip0=1,np0)
 200     format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
      write(6,*) '----------------------------------------------------'
      end if
      deallocate(pb05i)

      call set_nskpb0_vars(lmnp0,lmnp1,lmnp2,np1, np2, pb0,
     &     ip10, ip20, ip01,ip21,ip02,ip12,iu21,iv21)

      call tstatc0_end(id_sname)

      return
      end

c ================================ added by K. Tagami ================= 12.0A
      subroutine nskpa0_kt(jpr,tab,ng1,lsa1,iv1,np0,pa0,lmnp0,
     &                  lmnp1,lmnp2,
     &                  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &                  iu21,iv21, itrs )
c
      implicit real*8(a-h,o-z)
c
      real*8  tab(3,3)
      integer lsa1(3,3,ng1),iv1(ng1)
c
      real*8  pa0(3,lmnp0),pb0(3,lmnp0),p(3)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
c
      pb0 = pa0
c
      do 10 i=1,np0
        pa0(1,i)=pb0(1,i)*tab(1,1)+pb0(2,i)*tab(2,1)+pb0(3,i)*tab(3,1)
        pa0(2,i)=pb0(1,i)*tab(1,2)+pb0(2,i)*tab(2,2)+pb0(3,i)*tab(3,2)
        pa0(3,i)=pb0(1,i)*tab(1,3)+pb0(2,i)*tab(2,3)+pb0(3,i)*tab(3,3)
 10   continue
c
      np1=1
      ip10(1)=1
      ip01(1)=1
      do 20 i=2,np0
         do 22 j=1,np1
            jj=j
            call nspbge(pa0(1,i),pa0(1,ip01(j)),ind)
            if(ind.eq.0) go to 24
 22      continue
         np1=np1+1
         if(np1.gt.lmnp1) then
            write(6,*) ' np1=',np1,' > lmnp1=',lmnp1
            stop '=nskpb0(np1)='
         end if
         ip10(i)=np1
         ip01(np1)=i
         go to 20
 24      continue
        ip10(i)=jj
 20   continue
c
      np2=1
      ip21(1)=1
      ip12(1)=1
      ip02(1)=1
      iu21(1)=1
      iv21(1)=1
      do 30 i=2,np1
         do k=1,ng1
            do m=0, itrs
               kk=k
               p(1)=pa0(1,ip01(i))*lsa1(1,1,k)
     &              +pa0(2,ip01(i))*lsa1(2,1,k)
     &              +pa0(3,ip01(i))*lsa1(3,1,k)
               p(2)=pa0(1,ip01(i))*lsa1(1,2,k)
     &              +pa0(2,ip01(i))*lsa1(2,2,k)
     &              +pa0(3,ip01(i))*lsa1(3,2,k)
               p(3)=pa0(1,ip01(i))*lsa1(1,3,k)
     &              +pa0(2,ip01(i))*lsa1(2,3,k)
     &              +pa0(3,ip01(i))*lsa1(3,3,k)
               p(1) = p(1)* (1-2*m)
               p(2) = p(2)* (1-2*m)
               p(3) = p(3)* (1-2*m)

               do j=1,np2
                  jj=j
                  call nspbge(p,pa0(1,ip02(j)),ind)
                  if(ind.eq.0) go to 34
               end do
            end do
         end do
 32      continue

         np2=np2+1
         if(np2.gt.lmnp2) then
            write(6,*) ' np2=',np2,' > lmnp2=',lmnp2
            stop '=nskpb0(np2)='
         end if
         ip21(i)=np2
         ip12(np2)=i
         ip02(np2)=ip01(i)
         iu21(i)=1
         iv21(i)=1
         go to 30
 34      continue
         ip21(i)=jj
         iu21(i)=kk
         iv21(i)=iv1(kk)
 30   continue
c
      do 40 i=1,np0
 40       ip20(i)=ip21(ip10(i))
c
      if(jpr.ge.3) then
         write(6,*) ' === k points === '
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
         write(6,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
         write(6,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
 120          format((1h ,2('     (',3f9.5,' ) ')))
      end if
      if(jpr.ge.2) then
         write(6,*) '=== k-point code === '
         write(6,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
 140     format((1h ,2(i4,' (',3f9.5,' ) ')))
c
         write(6,*) '------------------------------------------------'
         write(6,180)
 180     format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
         write(6,200) (ip0,ip10(ip0),iu21(ip10(ip0)),ip21(ip10(ip0)),
     &        ip0=1,np0)
 200     format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
         write(6,*) '------------------------------------------------'
      end if
      return
      end
c ==================================================================== 12.0A

c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11        subroutine nskpbm(np0,pa0,np1,ip10,ip01)
c#12  input:
c#12
c#12
c#12
c#12  output:
c#12
c#12
c#13  external:  nspbge
c
c
c#41  2002.12.10 This subroutine is modified from nskpb0 for tetrahedron method in the program 'phase'.
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskpbm(np0,lmnp0,lmnp1,pa0,np1,ip10,ip01)
c
      implicit real*8(a-h,o-z)
c
      real*8  pa0(3,lmnp0)
      integer ip10(lmnp0),ip01(lmnp1)
c

      np1=1
      ip10(1)=1
      ip01(1)=1
      do 20 i=2,np0
        do 22 j=1,np1
          jj=j
          call nspbge(pa0(1,i),pa0(1,ip01(j)),ind)
          if(ind.eq.0) go to 24
   22   continue
        np1=np1+1
        ip10(i)=np1
        ip01(np1)=i
        go to 20
   24   continue
        ip10(i)=jj
   20 continue
c
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskpr0(jf,nxx,nyy,nzz,nx1,ny1,nz1,np0,np1,np2,pa0,pb0,
     &                  ip10,ip20,ip01,ip21,ip02,ip12,iu21,iv21)
c
      implicit real*8(a-h,o-z)
c     real*8  pa0(3,np0),pb0(3,np0)
c     integer ip10(np0),ip20(np0)
c     integer ip01(np1),ip21(np1)
c     integer ip02(np2),ip12(np2)
c     integer iu21(np1),iv21(np1)
      real*8  pa0(3,*),pb0(3,*)
      integer ip10(*),ip20(*)
      integer ip01(*),ip21(*)
      integer ip02(*),ip12(*)
      integer iu21(*),iv21(*)
c
      write(6,*) ' nskpr0 start!'
      read(jf,*) nxx,nyy,nzz,nx1,ny1,nz1
      read(jf,*) np0,np1,np2
      read(jf,200) ((pa0(j,i),j=1,3),i=1,np0)
      read(jf,200) ((pb0(j,i),j=1,3),i=1,np0)
      read(jf,*) (ip10(i),i=1,np0)
      read(jf,*) (ip20(i),i=1,np0)
      read(jf,*) (ip01(i),i=1,np1)
      read(jf,*) (ip21(i),i=1,np1)
      read(jf,*) (ip02(i),i=1,np2)
      read(jf,*) (ip12(i),i=1,np2)
      read(jf,*) (iu21(i),i=1,np1)
      read(jf,*) (iv21(i),i=1,np1)
  100 format(16i5)
  200 format(3d24.16)
      write(6,*) ' nskpr0 end!'
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nskpw0(jf,nxx,nyy,nzz,nx1,ny1,nz1,np0,np1,np2,pa0,pb0,
     &                  ip10,ip20,ip01,ip21,ip02,ip12,iu21,iv21)
c
      implicit real*8(a-h,o-z)
      real*8  pa0(3,np0),pb0(3,np0)
      integer ip10(np0),ip20(np0)
      integer ip01(np1),ip21(np1)
      integer ip02(np2),ip12(np2)
      integer iu21(np1),iv21(np1)
c
      write(jf,100) nxx,nyy,nzz,nx1,ny1,nz1
      write(jf,100) np0,np1,np2
      write(jf,200) ((pa0(j,i),j=1,3),i=1,np0)
      write(jf,200) ((pb0(j,i),j=1,3),i=1,np0)
      write(jf,100) (ip10(i),i=1,np0)
      write(jf,100) (ip20(i),i=1,np0)
      write(jf,100) (ip01(i),i=1,np1)
      write(jf,100) (ip21(i),i=1,np1)
      write(jf,100) (ip02(i),i=1,np2)
      write(jf,100) (ip12(i),i=1,np2)
      write(jf,100) (iu21(i),i=1,np1)
      write(jf,100) (iv21(i),i=1,np1)
  100 format(8i10)
  200 format(3d24.16)
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
c#11 &                                 grc,gkc,gra,gka,grb,gkb)
c
c#12  input:     jpr : print control
c#12              il : lattice type
c#12           a,b,c : lattice parameter (length)
c#12        ca,cb,cc : lattice parameter (cos(angle))
c#12  output:  tca(3,3),tac(3,3),tab(3,3),tba(3,3),tcb(3,3),tbc(3,3)
c#12           grc(3,3),gkc(3,3),gra(3,3),gka(3,3),grb(3,3),gkb(3,3)
c#13  external:  nslatc,nslata,nslatb
c
c#21  to get coordination systems (C,A,B)
c
c#31  1990.01.20.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
c ================================= modified by K. Tagami =============== 12.0A
c      subroutine nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
c     &                  grc,gkc,gra,gka,grb,gkb)
      subroutine nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
     &                  grc,gkc,gra,gka,grb,gkb,
     &                  use_altv_rltv, altv, rltv )
c ======================================================================== 12.0A
c
      implicit real*8(a-h,o-z)
      real*8 tca(3,3),tac(3,3),tab(3,3),tba(3,3),tcb(3,3),tbc(3,3)
      real*8 grc(3,3),gkc(3,3),gra(3,3),gka(3,3),grb(3,3),gkb(3,3)

c ===================== added by K. Tagami ============================= 12.0A
      logical use_altv_rltv
      real*8 altv(3,3), rltv(3,3)
c ====================================================================== 12.0A
c
      call nslatc(grc,gkc)

c =============================- modified by K. Tagami =============== 12.0A
c      call nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
c
      if ( use_altv_rltv ) then
         call nslata_kt( a,b,c,ca,cb,cc,tca,tac,gra,gka, altv, rltv )
      else
         call nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
      endif
c ======================================================================= 12.0A

      call nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)
c
      if(jpr.ge.3) then
      write(6,*) ' '
      write(6,*) ' === C coordinate system ==='
      write(6,*) ' (Cartesian Coordinate System)'
      write(6,*) ' metric tensors :  (real space) (l**2)      ',
     &           ' (reciprocal space) (1/l**2)'
      do 10 i=1,3
   10 write(6,100) (grc(i,j),j=1,3),(gkc(i,j),j=1,3)
      write(6,*) ' =================================================='
      end if
c
      if(jpr.ge.1) then
      write(6,*) ' '
      write(6,*) ' === A Coordinate system ==='
      write(6,*) ' (Conventional Crystalline Coordinate System)'
      write(6,*) ' transformation matrices Tca and Tac'
      do 20 i=1,3
   20 write(6,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
      write(6,*) ' metric tensors :  (real space) (l**2)      ',
     &           ' (reciprocal space) (1/l**2)'
      do 30 i=1,3
   30 write(6,100) (gra(i,j),j=1,3),(gka(i,j),j=1,3)
      write(6,*) ' ================================================'
        if(il.eq.-1) then
         write(6,*) ' < Trigonal lattice (R) >'
        else if(il.eq.0) then
          write(6,*) ' < Hexagonal lattice > '
        else if(il.eq.0 .or. il.eq.1) then
          write(6,*) ' < Primitive lattice >'
        else if(il.eq.2) then
          write(6,*) ' < face centered lattice (F) >'
        else if(il.eq.3) then
          write(6,*) ' < body centered lattice (I) >'
        else if(il.eq.4) then
          write(6,*) ' < one(ab)-face centered lattice (C) >'
        else if(il.eq.5) then
          write(6,*) ' < one(bc)-face centered lattice (C) >'
        else if(il.eq.6) then
          write(6,*) ' < one(ca)-face centered lattice (C) >'
        end if
      end if
c
      if(jpr.ge.2) then
      write(6,*) ' === B Coordinate system ==='
      write(6,*) ' (Primitive Crystalline Coordinate System)'
      write(6,*) ' transformation matrices Tab and Tba'
      do 40 i=1,3
   40 write(6,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
      write(6,*) ' transformation matrices Tcb and Tbc'
      do 42 i=1,3
   42 write(6,100) (tcb(i,j),j=1,3), (tbc(i,j),j=1,3)
      write(6,*) ' metric tensors :  (real space) (l**2)      ',
     &           ' (reciprocal space) (1/l**2)'
      do 50 i=1,3
   50 write(6,100) (grb(i,j),j=1,3),(gkb(i,j),j=1,3)
      write(6,*) ' ================================================'
      end if
      if(jpr.ge.0) then
          write(6,*) ' transformation matrices Tab and Tba'
          do 51 i=1,3
 51          write(6,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
      end if
  100 format(1h ,3f12.6,3x,3f12.6)
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
c
c#12  input:     a,b,c,ca,cb,cc: lattice parameters
c#12  output:    tca(3,3), tac(3,3): transformation matix
c#12                                 between C and A systems
c#12             gra(3,3),gka(3,3): metric tensors in A system
c#12  external:  nsmetr
c
c#21  to define the A coordinate system
c
c#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
c
      implicit real*8(a-h,o-z)
c
      real*8 tca(3,3),tac(3,3),gra(3,3),gka(3,3)
c
      eps=dfloat(10)**(-8)
c
      vv=sqrt(1-ca*ca-cb*cb-cc*cc+2*ca*cb*cc)
      sc=sqrt(1-cc*cc)
      p=(ca-cb*cc)/sc
      q=vv/sc
c
      tca(1,1)=a
      tca(2,1)=0
      tca(3,1)=0
      tca(1,2)=b*cc
      tca(2,2)=b*sc
      tca(3,2)=0
      tca(1,3)=c*cb
      tca(2,3)=c*p
      tca(3,3)=c*q
c
      tac(1,1)= 1/a
      tac(2,1)= 0
      tac(3,1)= 0
      tac(1,2)=-cc/(a*sc)
      tac(2,2)=  1/(b*sc)
      tac(3,2)= 0
      tac(1,3)= (ca*cc-cb)/(a*vv*sc)
      tac(2,3)=-p/(b*vv)
      tac(3,3)= sc/(c*vv)
c
c       error check.
      ier=0
      do 20 j=1,3
      do 20 i=1,3
        er1=0
        eg1=0
      do 22 k=1,3
        er1=er1+ tca(i,k)*tac(k,j)
        eg1=eg1+ tac(i,k)*tca(k,j)
   22 continue

        if(i.eq.j) then
          if(dabs(er1-1).gt.eps .or.
     &       dabs(eg1-1).gt.eps     ) then
             write(6,*) 'i,j=',i,j,'  e=',er1,eg1
             ier=1
          end if
        else
          if(dabs(er1).gt.eps .or.
     &       dabs(eg1).gt.eps     ) then
             write(6,*) 'i,j=',i,j,'  e=',er1,eg1
             ier=1
          end if
        end if
   20 continue
c
      if(ier.ne.0) then
      write(6,*) ' '
      write(6,*) ' === A Coordinate system ==='
      write(6,*) ' (Conventional Crystalline Coordinate System)'
      write(6,*) ' transformation matrices Tca and Tac'
      do 60 i=1,3
   60 write(6,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
  100 format(1h ,3f12.5,3x,3f12.5)
      stop ' === stop in sub.nslata ==='
      end if
c
      call nsmetr(tca,tac,gra,gka)
c
      return
      end

c =========================== added by K. Tagami ======================= 12.0A
      subroutine nslata_kt( a,b,c,ca,cb,cc,tca,tac,gra,gka,
     &                      altv, rltv )
      implicit none
c inout
      real*8 a, b, c, ca, cb, cc
      real*8 tca(3,3),tac(3,3),gra(3,3),gka(3,3)
      real*8 rltv(3,3), altv(3,3)
c local
      real*8 eps, PAI2, eg1, er1
      integer i, j, k, ier
c begin
      eps=dfloat(10)**(-8)
c
      Do i=1, 3
         Do k=1, 3
            tca(k,i) = altv(k,i)
         End do
      End do
c
      PAI2 = atan2(1.0d0,1.0d0) *8.0d0
      Do i=1, 3
         Do k=1, 3
            tac(k,i) = rltv(i,k) /PAI2     ! confirmed in the case of graphene
         End do
      End do

c ---   error check.
      ier=0
      do j=1,3
         do i=1,3
            er1=0
            eg1=0
            do k=1,3
               er1=er1+ tca(i,k)*tac(k,j)
               eg1=eg1+ tac(i,k)*tca(k,j)
            End do

            if(i.eq.j) then
               if(dabs(er1-1).gt.eps .or.
     &              dabs(eg1-1).gt.eps     ) then
                  write(6,*) 'i,j=',i,j,'  e=',er1,eg1
                  ier=1
               end if
            else
               if(dabs(er1).gt.eps .or.
     &              dabs(eg1).gt.eps     ) then
                  write(6,*) 'i,j=',i,j,'  e=',er1,eg1
                  ier=1
               end if
            end if
         End do
      End do
c
      if(ier.ne.0) then
         write(6,*) ' '
         write(6,*) ' === A Coordinate system ==='
         write(6,*) ' (Conventional Crystalline Coordinate System)'
         write(6,*) ' transformation matrices Tca and Tac'
         do i=1,3
            write(6,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
         End do
 100     format(1h ,3f12.5,3x,3f12.5)
         stop ' === stop in sub.nslata ==='
      end if
c
      call nsmetr(tca,tac,gra,gka)
c
      return
      end
c ==================================================================== 12.0A

c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)
c
c#12  input:     il                 : lattice type
c#12             tca(3,3), tac(3,3) : transformation matrix (C:A)
c#12  output:    tab(3,3), tba(3,3) : transformation matrix (A:B)
c#12             tcb(3,3), tbc(3,3) : transformation matrix (C:B)
c#12             grb(3,3), gkb(3,3) : metric tensors
c#12  external:  nsmetr
c
c#21  to define the B coordinate system
c
c#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)
c
      implicit real*8(a-h,o-z)
c
      real*8 tca(3,3),tac(3,3)
      real*8 tcb(3,3),tbc(3,3),grb(3,3),gkb(3,3)
      real*8 tab(3,3),tba(3,3)
c
      one=1
      two=2
c
      p1=one/2
      q1=one/3
      q2=two/3
      if(il.eq.-1) then
c       write(6,*) ' < Trigonal lattice (R) >'
        tab(1,1)= q2
        tab(2,1)= q1
        tab(3,1)= q1
        tab(1,2)=-q1
        tab(2,2)= q1
        tab(3,2)= q1
        tab(1,3)=-q1
        tab(2,3)=-q2
        tab(3,3)= q1
        tba(1,1)= 1
        tba(2,1)=-1
        tba(3,1)= 0
        tba(1,2)= 0
        tba(2,2)= 1
        tba(3,2)=-1
        tba(1,3)= 1
        tba(2,3)= 1
        tba(3,3)= 1
      else if(il.eq.0) then
c       write(6,*) ' < Hexagonal lattice > '
        do 10 j=1,3
        do 10 i=1,3
          tab(i,j)= 0
          tba(i,j)= 0
   10   continue
        do 12 i=1,3
          tab(i,i)= 1
          tba(i,i)= 1
   12   continue
      else if(il.eq.0 .or. il.eq.1) then
c       write(6,*) ' < Primitive lattice >'
        do 20 j=1,3
        do 20 i=1,3
          tab(i,j)= 0
          tba(i,j)= 0
   20   continue
        do 22 i=1,3
          tab(i,i)= 1
          tba(i,i)= 1
   22   continue
      else if(il.eq.2) then
c       write(6,*) ' < face centered lattice (F) >'
        tab(1,1)=  0
        tab(2,1)= p1
        tab(3,1)= p1
        tab(1,2)= p1
        tab(2,2)=  0
        tab(3,2)= p1
        tab(1,3)= p1
        tab(2,3)= p1
        tab(3,3)=  0
        tba(1,1)=-1
        tba(2,1)= 1
        tba(3,1)= 1
        tba(1,2)= 1
        tba(2,2)=-1
        tba(3,2)= 1
        tba(1,3)= 1
        tba(2,3)= 1
        tba(3,3)=-1
      else if(il.eq.3) then
c       write(6,*) ' < body centered lattice (I) >'
        tab(1,1)=-p1
        tab(2,1)= p1
        tab(3,1)= p1
        tab(1,2)= p1
        tab(2,2)=-p1
        tab(3,2)= p1
        tab(1,3)= p1
        tab(2,3)= p1
        tab(3,3)=-p1
        tba(1,1)= 0
        tba(2,1)= 1
        tba(3,1)= 1
        tba(1,2)= 1
        tba(2,2)= 0
        tba(3,2)= 1
        tba(1,3)= 1
        tba(2,3)= 1
        tba(3,3)= 0
      else if(il.eq.4) then
c       write(6,*) ' < one(ab)-face centered lattice (C) >'
        tab(1,1)= p1
        tab(2,1)=-p1
        tab(3,1)=  0
        tab(1,2)= p1
        tab(2,2)= p1
        tab(3,2)=  0
        tab(1,3)=  0
        tab(2,3)=  0
        tab(3,3)=  1
        tba(1,1)= 1
        tba(2,1)= 1
        tba(3,1)= 0
        tba(1,2)=-1
        tba(2,2)= 1
        tba(3,2)= 0
        tba(1,3)= 0
        tba(2,3)= 0
        tba(3,3)= 1
      else if(il.eq.5) then
c       write(6,*) ' < one(bc)-face centered lattice (C) >'
        tab(1,1)=  1
        tab(2,1)=  0
        tab(3,1)=  0
        tab(1,2)=  0
        tab(2,2)= p1
        tab(3,2)=-p1
        tab(1,3)=  0
        tab(2,3)= p1
        tab(3,3)= p1
        tba(1,1)= 1
        tba(2,1)= 0
        tba(3,1)= 0
        tba(1,2)= 0
        tba(2,2)= 1
        tba(3,2)= 1
        tba(1,3)= 0
        tba(2,3)=-1
        tba(3,3)= 1
      else if(il.eq.6) then
c       write(6,*) ' < one(ca)-face centered lattice (C) >'
        tab(1,1)= p1
        tab(2,1)=  0
        tab(3,1)= p1
        tab(1,2)=  0
        tab(2,2)=  1
        tab(3,2)=  0
        tab(1,3)=-p1
        tab(2,3)=  0
        tab(3,3)= p1
        tba(1,1)= 1
        tba(2,1)= 0
        tba(3,1)=-1
        tba(1,2)= 0
        tba(2,2)= 1
        tba(3,2)= 0
        tba(1,3)= 1
        tba(2,3)= 0
        tba(3,3)= 1
      else
        write(6,*) ' il=',il
        write(6,*) ' === stop in sub.nslatb. (il) ==='
        stop 'error in sub.nslatb'
      end if
c
      do 30 i=1,3
      do 30 j=1,3
        x=0.0
        y=0.0
        do 32 k=1,3
          x=x+tca(i,k)*tab(k,j)
          y=y+tba(i,k)*tac(k,j)
   32   continue
        tcb(i,j)=x
        tbc(i,j)=y
   30 continue
c
c     write(6,*) ' === B Coordinate system ==='
c     write(6,*) ' (Primitive Crystalline Coordinate System)'
c     write(6,*) ' transformation matrices Tab and Tba'
c     do 40 i=1,3
c  40 write(6,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
c     write(6,*) ' transformation matrices Tcb and Tbc'
c     do 42 i=1,3
c  42 write(6,100) (tcb(i,j),j=1,3), (tbc(i,j),j=1,3)
c 100 format(1h ,3f12.5,3x,3f12.5)
c
      call nsmetr(tcb,tbc,grb,gkb)
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nslatc(grc,gkc)
c
c#12  noinput:
c#12  output:    grc(3,3): metric tensor in real space
c#12             gkc(3,3): metric tensor in reciprocal space
c#12  noexternal:
c
c#21  to get metric tensors in C (Cartesian) coordinate system
c
c#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nslatc(grc,gkc)
c
      implicit real*8(a-h,o-z)
      real*8 grc(3,3),gkc(3,3)
c
      one=1
      pai2=8*datan(one)
      pai22=pai2**2
      do 10 j=1,3
      do 10 i=1,3
        grc(i,j)=0
        gkc(i,j)=0
   10 continue
c
      do 20 i=1,3
        grc(i,i)=1
        gkc(i,i)=pai22
   20 continue
c
c     write(6,*) ' '
c     write(6,*) ' === C coordinate system ==='
c     write(6,*) ' (Cartesian Coordinate System)'
c     write(6,*) ' metric tensors :  (real space) (l**2)      ',
c    &           ' (reciprocal space) (1/l**2)'
c     do 30 i=1,3
c  30 write(6,100) (grc(i,j),j=1,3),(gkc(i,j),j=1,3)
c 100 format(1h ,3f12.6,3x,3f12.6)
c     write(6,*) ' =================================================='
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nslatr(il,r,ind)
c
c     il :  lattice type
c     ind = 0 :  r is lattice vector.
c           1 :  r in not lattice vector.
c
      implicit real*8(a-h,o-z)
      dimension r(3), j(3)
      data one/1.0/, half/0.5/, eps/1.e-4/
c
      thrd=1.0/3.0
      if(il.le.-2 .or. il.ge.5) then
        write(6,*) ' il=',il
        stop '=== stop in sub.nslatr. (il) ==='
      end if
      ind=1
      do 10 i=1,3
        aaa=mod(abs(r(i)),one)
        if(aaa.gt.eps .and. abs(one-aaa).gt.eps) go to 2
   10 continue
      ind=0
      return
c
    2 continue
      if(il.eq.-1) then
        do 20 i=1,3
          aaa=mod(abs(r(i)),thrd)
          if(aaa.gt.eps .and. abs(thrd-aaa).gt.eps) return
   20   continue
        aaa=mod(abs(-r(1)+r(2)+r(3)),one)
        if(aaa.lt.eps .or. abs(one-aaa).lt.eps) ind=0
      else if(il.le.1) then
        return
      else if(il.eq.2 .or. il.eq.3) then
        do 30 i=1,3
          aaa=mod(abs(r(i)),half)
          if(aaa.gt.eps .and. abs(half-aaa).gt.eps) return
   30   continue
        if(il.eq.2) then
          aaa=mod(abs(r(1)+r(2)+r(3)),one)
          if(aaa.lt.eps .or. abs(one-aaa).lt.eps) ind=0
        else if(il.eq.3) then
          do 40 i=1,3
            s=sign(0.1d0,r(i))
            j(i)=iabs(int(2.0*r(i)+s))
   40     continue
          if((mod(j(1),2).eq.0 .and. mod(j(2),2).eq.0 .and.
     &        mod(j(3),2).eq.0) .or. (mod(j(1),2).eq.1 .and.
     &        mod(j(2),2).eq.1 .and. mod(j(3),2).eq.1)) ind=0
        end if
      else if(il.eq.4) then
        do 50 i=1,2
          aaa=mod(abs(r(i)),half)
          if(aaa.gt.eps .and. abs(half-aaa).gt.eps) return
   50   continue
          aaa=mod(abs(r(3)),one)
          if(aaa.gt.eps .and. abs(one-aaa).gt.eps) return
          aaa=mod(abs(r(1)+r(2)),one)
          if(aaa.lt.eps .or. abs(one-aaa).lt.eps) ind=0
      end if
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
c
c#12  input:       il : lattice type
c#12              ng1 : # of group elements (order of the space group)
c#12             ig01 : list vector for getting the element code
c#12  in-output:   a,b,c,ca,cb,cc: lattice parameters
c#12  noexternal:
c
c#21  to check the consistency between the space group and
c#21                                  the lattice parameters
c
c#31  1990.04.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
c
      implicit real*8(a-h,o-z)
c
      integer il,ng1,ig01(48)
      logical la,lb,lc,lcub,ltet
c
      one=1
      if(il.le.0) then
        b=a
        cc=-one/2
        ca= 0
        cb= 0
      else
        la=.false.
        lb=.false.
        lc=.false.
        lcub=.false.
        ltet=.false.
        do 3 i=1,ng1
          jj=ig01(i)
          if(jj.eq.4)  lc=.true.
          if(jj.eq.28) lc=.true.
          if(jj.eq.2)  la=.true.
          if(jj.eq.26) la=.true.
          if(jj.eq.27) lb=.true.
          if(jj.eq.3)  lb=.true.
          if(jj.eq.5)  lcub=.true.
          if(jj.eq.21) ltet=.true.
          if(jj.eq.48) ltet=.true.
    3   continue
        if(lc) ca=0
        if(lc) cb=0
        if(la) cb=0
        if(la) cc=0
        if(lb) cc=0
        if(lb) ca=0
        if(ltet) b=a
        if(lcub) b=a
        if(lcub) c=a
      end if
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsm3iw(jf,ng1,ir)
c
c#12  input:      jf : output file
c#12             ng1 : # of matrices (3x3)
c#12              ir : matrix
c#12  nooutput
c#13  noexternal
c
c#21  to write 3x3 matrices
c
c#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsm3iw(jf,ng1,ir)
c
      implicit real*8(a-h,o-z)
c
      integer ir(3,3,ng1)
c
      no=ng1/6
      mo=mod(ng1,6)
      do 60 ii=1,no
      write(jf,*) '   '
      i1=(ii-1)*6+1
      i2=ii*6
      write(jf,200) (i,(ir(1,k,i),k=1,3),i=i1,i2)
      do 62 j=2,3
   62 write(jf,220) (  (ir(j,k,i),k=1,3),i=i1,i2)
   60 continue
c
      if(mo.ne.0) then
        i1=no*6+1
        i2=no*6+mo
        write(jf,*) '   '
        write(jf,200) (i,(ir(1,k,i),k=1,3),i=i1,i2)
        do 64 j=2,3
   64   write(jf,220) (  (ir(j,k,i),k=1,3),i=i1,i2)
      end if
c
  200 format(1h ,6('(',i2,')',3i2,2x))
  220 format(1h ,6(4x        ,3i2,2x))
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsmetr(tca,tac,gra,gka)
c
c#12  input:     tca(3,3), tac(3,3) : transformation matrix
c#12  output:    gra(3,3), gka(3,3)  : metric tensors
c#12  noexternal:
c
c#21  to get metric tensors
c
c#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsmetr(tca,tac,gra,gka)
c
      implicit real*8(a-h,o-z)
      real*8 tca(3,3),tac(3,3),gra(3,3),gka(3,3)
c
      one=1
      pai2=8*datan(one)
      pai22=pai2**2
c
      do 10 i=1,3
      do 10 j=1,3
        x=0
        y=0
        do 20 k=1,3
          x=x+tca(k,i)*tca(k,j)
          y=y+tac(i,k)*tac(j,k)
   20   continue
        gra(i,j)=x
        gka(i,j)=y*pai22
   10 continue
c     write(6,*) ' metric tensors :  (real space) (l**2)      ',
c    &           ' (reciprocal space) (1/l**2)'
c     do 30 i=1,3
c  30 write(6,100) (gra(i,j),j=1,3),(gka(i,j),j=1,3)
c 100 format(1h ,3f12.6,3x,3f12.6)
c     write(6,*) ' ================================================'
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsmlt1(nn,n,a,b,c)
c
c#12  input :  nn
c#12            n        : n=<nn
c#12            a(nn,n)  : n*n matrix
c#12            b(nn,n)  : n*n matrix
c#12  output    c(nn,n)  :  c=a*b
c#13  noexternal
c
c#21  to get matrix (n*n)  multiplication
c
c#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsmlt1(nn,n,a,b,c)
c
      implicit real*8(a-h,o-z)
c
      dimension a(nn,n),b(nn,n),c(nn,n)
c
      do 10 i=1,n
      do 10 j=1,n
         c(i,j)=0
      do 10 k=1,n
        c(i,j)=c(i,j)+a(i,k)*b(k,j)
   10 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsmltw(jf,ng1,im1,iv1)
c
c#12  input:      jf : output file
c#12             ng1 : order of the group
c#12      im1(48,48) : multiplication table
c#12         iv1(48) : inverse element
c#12  nooutput
c#13  noexternal
c
c#21  to write multiplication table
c
c#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsmltw(jf,ng1,im1,iv1)
c
      implicit real*8(a-h,o-z)
c
      integer ng1,im1(48,48),iv1(48)
c
      write(jf,*) '   '
      write(jf,*) '--- group multiplication table ---'
      if(ng1.le.24) then
        write(jf,320) (j,j=1,ng1)
        write(jf,340) ('---',j=1,ng1)
        do 80 i=1,ng1
   80   write(jf,300) i,(im1(i,j),j=1,ng1)
      else
        write(jf,*) ' '
        write(jf,320) (j,j=1,24)
        write(jf,340) ('---',j=1,24)
        do 82 i=1,ng1
   82   write(jf,300) i,(im1(i,j),j=1,24)
        write(jf,*) ' '
        write(jf,320) (j,j=25,ng1)
        write(jf,340) ('---',j=25,ng1)
        do 84 i=1,ng1
   84   write(jf,300) i,(im1(i,j),j=25,ng1)
      end if
  300 format((i3,2x,24i3))
  320 format((5x,24i3))
  340 format(5x,24a3)
c
      write(jf,*) ' '
      write(jf,*) '--- invers elements ---'
      write(jf,400) (iv1(j),j=1,ng1)
  400 format((5x,24i3))
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsmult(jf,jpr,ng0,irot,im,iv)
c
c#12  input :          jf   : output file
c#12                  jpr   : output control
c#12                  ng0   : # of group elements
c#12         irot(3,3,48)   : rotation matrix
c#12  output:     im(48,48) : multiplication table
c#12              iv(48)    : inverse elements
c#13  noexternal
c
c#21  to get multiplication table from rotation matrices
c
c#31  1989.12.26.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsmult(jf,jpr,ng0,irot,im,iv)
c
      implicit real*8(a-h,o-z)
c
      integer irot(3,3,48),ia(3,3),im(48,48),iv(48)
c
      do 20 i=1,ng0
      do 20 j=1,ng0
        do 22 k=1,3
        do 22 l=1,3
          ia(k,l)=0
          do 24 m=1,3
   24     ia(k,l)=ia(k,l)+irot(k,m,i)*irot(m,l,j)
   22   continue
        do 26 k=1,ng0
          if(ia(1,1).eq.irot(1,1,k) .and.
     &       ia(2,1).eq.irot(2,1,k) .and.
     &       ia(3,1).eq.irot(3,1,k) .and.
     &       ia(1,2).eq.irot(1,2,k) .and.
     &       ia(2,2).eq.irot(2,2,k) .and.
     &       ia(3,2).eq.irot(3,2,k) .and.
     &       ia(1,3).eq.irot(1,3,k) .and.
     &       ia(2,3).eq.irot(2,3,k) .and.
     &       ia(3,3).eq.irot(3,3,k)      ) go to 28
   26   continue
        stop '=== stop in sub.nsmult. (im) ==='
   28   im(i,j)=k
   20 continue
c
      do 30 i=1,ng0
        do 32 j=1,ng0
          jj=j
          if(im(i,j).eq.1) go to 34
   32   continue
        stop '=== stop in sub.nsmult. (iv) ==='
   34   iv(i)=jj
   30 continue
c
      if(jpr.ge.3) then
        write(jf,*) '   '
        write(jf,*) '--- group multiplication table ---'
        write(jf,120) (j,j=1,24)
        write(jf,140)
        do 40 i=1,24
   40   write(jf,100) i,(im(i,j),j=1,24)
c
        if (ng0.gt.24) then
        write(jf,*) ' '
        write(jf,120) (j,j=1,24)
        write(jf,140)
        do 42 i=25,ng0
   42   write(jf,100) i,(im(i,j),j=1,24)
        write(jf,*) ' '
        write(jf,120) (j,j=25,ng0)
        write(jf,140)
        do 44 i=1,24
   44   write(jf,100) i,(im(i,j),j=25,ng0)
        write(jf,*) ' '
        write(jf,120) (j,j=25,ng0)
        write(jf,140)
        do 46 i=25,ng0
   46   write(jf,100) i,(im(i,j),j=25,ng0)
        end if
c
        write(jf,*) ' '
        write(jf,*) '--- invers elements ---'
        write(jf,120) (iv(j),j=1,ng0)
  100   format((i3,2x,24i3))
  120   format((5x,24i3))
  140   format(5x,72('-'))
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia)
c
c#12  input:     jpr     : print control
c#12             rm      : maximum of neighbor distance (standard unit)
c#12             grb(3,3): metric tensor
c#12             nat     : # of atoms in a unit cell
c#12             u(3,nat): atomic positions in a unit cell
c#12
c#12  output:    nb(nat) : # of neighbors
c#12             lmnb    : limit of # of neighbors
c#12             r (lmnb,nat): distance to the neighbors (standard unit)
c#12             ia(lmnb,nat): atom index of the neighbor
c#12  noexternal:
c
c#21  to get a distance from each atom to the neighboring atoms
c
c#31  1990.11.29.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia,w)
c
      implicit real*8(a-h,o-z)
      real*8  grb(3,3)
      real*8  u(3,nat), r(lmnb,nat), w(3,lmnb,nat)
      integer ia(lmnb,nat), nb(nat)
c
      one=1
      pai2=8*datan(one)
      pai22=pai2**2
c
      nx=0
      ny=0
      nz=0
      do 30 iz=-1,1
      do 32 iy=-1,1
      do 34 ix=-1,1
      if(ix.ne.0 .and. iy.ne.0 .and. iz.ne.0) then
        p1=ix
        p2=iy
        p3=iz
        rr=sqrt(p1*grb(1,1)*p1+p1*grb(1,2)*p2+p1*grb(1,3)*p3
     &         +p2*grb(2,1)*p1+p2*grb(2,2)*p2+p2*grb(2,3)*p3
     &         +p3*grb(3,1)*p1+p3*grb(3,2)*p2+p3*grb(3,3)*p3)
        nn=(rm/rr)+5
        nnx=iabs(ix)*nn
        nny=iabs(iy)*nn
        nnz=iabs(iz)*nn
        if(nnx.gt.nx) nx=nnx
        if(nny.gt.ny) ny=nny
        if(nnz.gt.nz) nz=nnz
      endif
   34 continue
   32 continue
   30 continue
c
c     write(6,'('' nx,ny,nz='',3i5)') nx,ny,nz
c
c
      do 10 j0=1,nat
        v1=u(1,j0)
        v2=u(2,j0)
        v3=u(3,j0)
        k=0
        do 20 j1=1,nat
          ww1=u(1,j1)
          ww2=u(2,j1)
          ww3=u(3,j1)
          do 21 iz=-nz,nz
            w3=ww3+iz-v3
          do 22 iy=-ny,ny
            w2=ww2+iy-v2
          do 24 ix=-nx,nx
            w1=ww1+ix-v1
            rr=sqrt(w1*grb(1,1)*w1+w1*grb(1,2)*w2+w1*grb(1,3)*w3
     &             +w2*grb(2,1)*w1+w2*grb(2,2)*w2+w2*grb(2,3)*w3
     &             +w3*grb(3,1)*w1+w3*grb(3,2)*w2+w3*grb(3,3)*w3)
c                                                 registration
            if(rr.gt.1.e-5 .and. rr.le.rm) then
              k=k+1
              if(k.gt.lmnb) then
                write(6,*) ' k=',k,'  > lmnb=',lmnb
                stop 'sub.nsnbat'
              end if
              r (k,j0)=rr
              ia(k,j0)=j1
              w(1,k,j0)=w1
              w(2,k,j0)=w2
              w(3,k,j0)=w3
c             write(6,900) j0,k,ia(k,j0),r(k,j0)
            end if
   24     continue
   22     continue
   21     continue
   20   continue
        nb(j0)=k
c
   10 continue
c
      do 90 j0=1,nat
        call nsnbst(nb(j0),r(1,j0),ia(1,j0),w(1,1,j0))
        if(jpr.ge.2) then
          write(6,*) ' atom0=',j0
          do 92  k=1,nb(j0)
          write(6,900) k,ia(k,j0),(w(i,k,j0),i=1,3),r(k,j0)
   92     continue
        end if
   90 continue
  900 format(1h ,'atom1=',i5,' (',i5,')','  w=',3f8.3,'   r=',f10.6)
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsnbrd(jf,jpr,nat,nb,lmnb,wb)
c
      implicit real*8(a-h,o-z)
c
      integer nb(nat)
      real*8  wb(3,lmnb,nat)
c
      read(jf) nat
      read(jf) (nb(ia0),ia0=1,nat)
      do 10 ia0=1,nat
        read(jf) ((wb(j,ib1,ia0),j=1,3),ib1=1,nb(ia0))
   10 continue
c
      if(jpr.gt.3) then
        write(6,*) ' nb:'
        write(6,700) (nb(i),i=1,nat)
        write(6,*) ' ib1  ia0 :   wb'
        do 60 ia0=1,nat
        do 60 ib1=1,nb(ia0)
          write(6,780) ib1,ia0,(wb(j,ib1,ia0),j=1,3)
   60   continue
      end if
c
  700 format(20i3)
  780 format(2i4,2x,3f12.6)
c
      return
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsnbst(n,ra,ia,wa)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  ia(n): an array to be rearranged correspondingly
c#12                  wa(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#21  to sort ra(n) in ascending order with an index ia(n)
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsnbst(n,ra,ia,wa)
c
      implicit real*8(a-h,o-z)
      real*8  ra(n),wa(3,n)
      integer ia(n)
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          iia=ia(l)
          w1a=wa(1,l)
          w2a=wa(2,l)
          w3a=wa(3,l)
        else
          rra=ra(ir)
          iia=ia(ir)
          w1a=wa(1,ir)
          w2a=wa(2,ir)
          w3a=wa(3,ir)
          ra(ir)=ra(1)
          ia(ir)=ia(1)
          wa(1,ir)=wa(1,1)
          wa(2,ir)=wa(2,1)
          wa(3,ir)=wa(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            ia(1)=iia
            wa(1,1)=w1a
            wa(2,1)=w2a
            wa(3,1)=w3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            ia(i)=ia(j)
            wa(1,i)=wa(1,j)
            wa(2,i)=wa(2,j)
            wa(3,i)=wa(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        ia(i)=iia
        wa(1,i)=w1a
        wa(2,i)=w2a
        wa(3,i)=w3a
c
      go to 10
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nspace(jf,jpr,il,ngen,inv,igen,jgen,
c#11 &           ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1,
c#11 &           omove,euler1,inver1)
c
c#12  input:      jf : output file
c#12             jpr : print control
c#12              il : lattice type
c#12             inv : parameter (0,1) for moving the origin
c#12         igen(3) : generator (rotation part)
c#12     jgen(2,3,3) : generator (nonprimitive translation part)
c#12  output:       ng00 : # of elements in Schoenflies notation
c#12          schoe0(48) : Schoenflies notation (character*5)
c#12                 ng1 : # of group elements
c#12                ig01 : element code
c#12                 ta1 : nonprimitive translation vector (A system)
c#12                 ra1 : rotation matrix in real space (A system)
c#12                 sa1 : rotation matrix in reciprocal space (A system)
c#12                       ra1*sa1=unit matrix
c#12          im1(48,48) : multiplication table
c#12             iv1(48) : inverse element
c#12            omove(3) : translation of the origin from the initial
c#12        euler1(3,48) : Euler angle
c#12        inver1(  48) : index for inversion operation (=1 or -1)
c#13  external:  nsgrp1, tspaca
c
c#21  to get space group
c
c#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
c ================================== modified by K. Tagami ============ 12.0A
c      subroutine nspace(jf,jpr,il,ngen,inv,igen,jgen,
c     &                  ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1,
c     &                  omove,euler1,inver1)
      subroutine nspace(jf,jpr,il,ngen,inv,igen,jgen,
     &                  ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1,
     &                  omove,euler1,inver1,
     &                  use_trs, tac, tca, tab, tba,
     &                  gen_name_in_carts )
c ====================================================================== 12.0A
c
      implicit real*8(a-h,o-z)
      integer igen(ngen),jgen(2,3,ngen)
      real*8  ta1(3,48),ra1(3,3,48),sa1(3,3,48)
      real*8  omove(3)
c ====================================== added by K. Tagami ============ 12.0A
      logical use_trs, gen_name_in_carts
      real*8 tac(3,3), tca(3,3), tab(3,3), tba(3,3)
c ====================================================================== 12.0A
c
c
      character*2 jones
      character*5 schoen
      common/nspg00/ schoen(48),jones(3,48)
      common/nspg0 / euler(3,24), rot(3,3,48),
     &              ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &               im0(48,48),iv0(48),ng0
c     common/nspg1 / il,ng1,ig10(48),ig01(48),jg1(2,3,48),
c    &               irotr1(3,3,48),irotk1(3,3,48),im1(48,48),iv1(48)
c
c     real*8  euler(3,24), rot(3,3,48)
c     integer ieuler(3,24),irot(3,3,48),ir1234(3,48),
c    &        im0(48,48),iv0(48),ng0
c
      character*5 schoe0(48)
      integer il,ng1,ig10(48),ig01(48),jg1(2,3,48),
     &        irotr1(3,3,48),irotk1(3,3,48),im1(48,48),iv1(48)
      integer movo(2,3)
      real*8  euler1(3,48)
      integer inver1(48)
c
c ====================================== modified by K. Tagami ========= 12.0A
c      call nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
c     &            schoen,jones,euler,rot,ieuler,irot,ir1234,
c     &            im0,iv0,ng0,ng1,ig10,ig01,jg1,
c     &            irotr1,irotk1,im1,iv1,movo)
      call nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
     &            schoen,jones,euler,rot,ieuler,irot,ir1234,
     &            im0,iv0,ng0,ng1,ig10,ig01,jg1,
     &            irotr1,irotk1,im1,iv1,movo,
     &            use_trs, tac, tca, tab, tba,
     &            gen_name_in_carts )
c ======================================================================= 12.0A
c
      do 10 k=1,ng1
        ta1(1,k)=dfloat(jg1(1,1,k))/jg1(2,1,k)
        ta1(2,k)=dfloat(jg1(1,2,k))/jg1(2,2,k)
        ta1(3,k)=dfloat(jg1(1,3,k))/jg1(2,3,k)
   10 continue
c
      do 12 k=1,ng1
        do 14 j=1,3
        do 14 i=1,3
          ra1(i,j,k)=irotr1(i,j,k)
          sa1(i,j,k)=irotk1(i,j,k)
   14   continue
   12 continue
c
      do 20 k=1,ng1
        ig0=ig01(k)
        if(il.le.0) then
          if(ig0.le.12) then
            ii=ig0
            inver1(k)= 1
          else
            ii=ig0-12
            inver1(k)=-1
          end if
        else
          if(ig0.le.24) then
            ii=ig0
            inver1(k)= 1
          else
            ii=ig0-24
            inver1(k)=-1
          end if
        end if
        euler1(1,k)=euler(1,ii)
        euler1(2,k)=euler(2,ii)
        euler1(3,k)=euler(3,ii)
   20 continue
c
      ng00=ng0
      do 30 k=1,ng0
        schoe0(k)=schoen(k)
   30 continue

c
      omove(1)=dfloat(movo(1,1))/movo(2,1)
      omove(2)=dfloat(movo(1,2))/movo(2,2)
      omove(3)=dfloat(movo(1,3))/movo(2,3)
c
c     for being compatible with the tspace package
      call tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nspbg0(pb,np0,pb0,ip0)
c
      implicit real*8(a-h,o-z)
      real*8 pb(3),pb0(3,np0)
      do 10 ip=1,np0
        ip0=ip
        call nspbge(pb,pb0(1,ip0),ind)
        if(ind.eq.0) go to 12
   10 continue
      write(6,*) ' pb=',pb,' : not suitable (sub.nspbg0)'
      stop 'error(sub.nspbg0)'
   12 return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nspbge(p1,p2,ind)
c
      implicit real*8(a-h,o-z)
c
      real*8 p1(3),p2(3)
c
      eps=dfloat(10)**(-5)
      one=1
      x=abs(p1(1)-p2(1))+eps/2
      y=abs(p1(2)-p2(2))+eps/2
      z=abs(p1(3)-p2(3))+eps/2
      if(dmod(x,one).le.eps .and.
     &   dmod(y,one).le.eps .and.
     &   dmod(z,one).le.eps      ) then
        ind=0
      else
        ind=1
      end if
      return
      end

      subroutine nspbgei(p1,p2,ind)
c
      implicit real*8(a-h,o-z)
c
      integer :: p1(3),p2(3)
      integer :: x,y,z
      integer :: one = 100000
c
      x=abs(p1(1)-p2(1))
      y=abs(p1(2)-p2(2))
      z=abs(p1(3)-p2(3))
      if(mod(x,one).eq.0 .and.
     &   mod(y,one).eq.0 .and.
     &   mod(z,one).eq.0      ) then
        ind=0
      else
        ind=1
      end if
      return
      end

c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nspgr1(jf,tca,tac,tab,tba,tcb,tbc,
     &                     grc,gkc,gra,gka,grb,gkb,schoen,
     &           il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1,
     &           euler1,inver1)
c
      implicit real*8(a-h,o-z)
      real*8 tca(3,3),tac(3,3),tab(3,3),tba(3,3),tcb(3,3),tbc(3,3)
      real*8 grc(3,3),gkc(3,3),gra(3,3),gka(3,3),grb(3,3),gkb(3,3)
      character*5 schoen(48)
      integer il,ng0,ng1,ig01(48),im1(48,48),iv1(48)
      integer lra1(3,3,48),lsa1(3,3,48),lrb1(3,3,48),lsb1(3,3,48)
      real*8  ta1(3,48),tb1(3,48)
      real*8  euler1(3,48)
      integer inver1(  48)
c
      read(jf,100) ((tca(i,j),i=1,3),j=1,3)
      read(jf,100) ((tac(i,j),i=1,3),j=1,3)
      read(jf,100) ((tab(i,j),i=1,3),j=1,3)
      read(jf,100) ((tba(i,j),i=1,3),j=1,3)
      read(jf,100) ((tcb(i,j),i=1,3),j=1,3)
      read(jf,100) ((tbc(i,j),i=1,3),j=1,3)
c
      read(jf,100) ((grc(i,j),i=1,3),j=1,3)
      read(jf,100) ((gkc(i,j),i=1,3),j=1,3)
      read(jf,100) ((gra(i,j),i=1,3),j=1,3)
      read(jf,100) ((gka(i,j),i=1,3),j=1,3)
      read(jf,100) ((grb(i,j),i=1,3),j=1,3)
      read(jf,100) ((gkb(i,j),i=1,3),j=1,3)
c
      read(jf,200) il,ng0
      read(jf,300) (schoen(i),i=1,ng0)
      read(jf,200) ng1
      read(jf,200) (ig01(i),i=1,ng1)
      read(jf,200) ((im1(i,j),i=1,ng1),j=1,ng1)
      read(jf,200) ( iv1(i),i=1,ng1)
      read(jf,100) ((ta1(i,j),i=1,3),j=1,ng1)
      read(jf,220) (((lra1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      read(jf,220) (((lsa1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      read(jf,100) ((tb1(i,j),i=1,3),j=1,ng1)
      read(jf,220) (((lrb1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      read(jf,220) (((lsb1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      read(jf,100) ((euler1(i,j),i=1,3),j=1,ng1)
      read(jf,200) ( inver1(i),i=1,ng1)
c
  100 format((3d24.16))
  200 format((24i3))
  220 format((18i3))
  300 format((10(2x,a5)))
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im,ir,schoen,mo,ng,ig,jg)
c
c#12  input:    jf  : output file
c#12           jpr  : output control
c#12            il  : lattice type
c#12           inv  : flug for moving the origin (0:nomove, 1:move)
c#12          ngen  : # of generator
c#12          igen  : rotation of generater
c#12          jgen  : nonprimitive translation vector of generator
c#12            im  : multiplication table
c#12            ir  : rotation matrix
c#12        schoen  : schoenfries index
c#12  output:   mo(2,3)    : translation of the origin
c#12            ng         : # of operations (order of group)
c#12            ig(48)     : rotation of group element
c#12            jg(2,3,48) : nonprimitive translation of the element
c#13  external: sub.nsrot1, nssum1
c
c#21  to get space group
c
c#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im,ir,schoen,
     &                  mo,ng,ig,jg)
c
c
      implicit real*8(a-h,o-z)
      integer igen(*),jgen(2,3,*),ig(*),jg(2,3,*)
      integer im(48,48),ir(3,3,48)
      integer iflg(48),jw(2,3,48),jww(2,3),mo(2,3),mo1(2,3),jw1(2,3)
      integer lprm(2,3,4)
      character*5 schoen(48)
c
c --- register the generator ---
c
      iflg(1)=1
      do 10 j=1,3
        jw(1,j,1)=0
        jw(2,j,1)=1
   10 continue
      do 12 i=2,48
   12 iflg(i)=0
      do 14 i=1,ngen
        if(il.ge.1 .and. igen(i).gt.48) then
          write(6,*) ' igen=',igen(i),' > 48 : error'
          stop 'error:sub.nspgrp (igen)'
        end if
        if(il.le.0 .and. igen(i).gt.24) then
          write(6,*) ' igen=',igen(i),' > 24 : error'
          stop 'error:sub.nspgrp (igen)'
        end if
        iflg(igen(i))=1
        do 16 j=1,3
          jw(1,j,igen(i))=jgen(1,j,i)
          jw(2,j,igen(i))=jgen(2,j,i)
   16   continue
   14 continue
      n=0
      do 18 i=1,48
   18 n=n+iflg(i)
      if(jpr.ge.0)      write(jf,*) ' '
      if(jpr.ge.0)      write(jf,*) ' ----- generator -------'
      if(jpr.ge.0) then
      nn=0
      do 80 i=1,48
        if(iflg(i).eq.1) then
          nn=nn+1
          write(jf,100) nn,i,schoen(i),((jw(j,k,i),j=1,2),k=1,3)
        end if
   80 continue
      end if
c
c --- produce group elements by multiplications
c
    1 ng=n
      do 20 i=1,48
      if(iflg(i).eq.1) then
        do 22 j=1,48
        if(iflg(j).eq.1) then
c         iw=im(i,j)
          iw=im(j,i)
          if(iflg(iw).ne.1) then
            iflg(iw)=1
            call nsrot1(jw(1,1,i),ir(1,1,j),jww)
            call nssum1(1,jww,jw(1,1,j),jw(1,1,iw))
          end if
        end if
   22   continue
      end if
   20 continue
      n=0
      do 24 i=1,48
   24 n=n+iflg(i)
      if(n.gt.ng) go to 1
c
c --- register group elements
c
      n=0
      do 30 i=1,48
      if(iflg(i).eq.1) then
        n=n+1
        ig(n)=i
        do 32 k=1,3
        do 32 j=1,2
   32   jg(j,k,n)=jw(j,k,i)
      end if
   30 continue
c
c --- move the origin
c
      do 38 i=1,3
        mo(1,i)=0
        mo(2,i)=1
   38 continue
      if(inv.ne.0) then
        i1=0
        do 40 i=1,ng
          if(il.le.0 .and. ig(i).eq.13) i1=i
          if(il.ge.1 .and. ig(i).eq.25) i1=i
   40   continue
        if(i1.ne.0) then
          do 42 i=1,3
            mo(1,i)=jg(1,i,i1)
            mo(2,i)=jg(2,i,i1)*2
            call nsrduc(mo(1,i),mo(2,i))
   42     continue
          if(jpr.ge.1) write(6,200) ((mo(i,j),i=1,2),j=1,3)
          do 44 i=1,ng
            call nsrot1(mo,ir(1,1,ig(i)),mo1)
            call nssum1( 1,jg(1,1,i),mo1,jw1)
            call nssum1(-1,jw1,mo,jg(1,1,i))
   44     continue
        end if
      end if
c
      call nsprmv(il,nprm,lprm)
      call nssdjg(jf,jpr,nprm,lprm,ng,jg)
c
c --- print
c

      if(jpr.ge.2) then
        write(jf,*) ' '
        write(jf,*) ' ----- group elements -------'
        do 90 i=1,ng
   90   write(jf,100) i,ig(i),schoen(ig(i)),
     &                ((jg(j,k,i),j=1,2),k=1,3)
  100   format(i5,'   (',i2,')',3x,a5,5x,
     &         '(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )
      end if
  200   format(/'    move the origin by ','  (',2(i3,' /',i3,'  ,'),
     &                                         i3,' /',i3,' )'      )
c
      return
      end

c ==================================== added by K. Tagami =============== 12.0A
      subroutine nspgrp_kt( jf, jpr, il, inv, ngen, igen, jgen,
     &                      im, ir, schoen,
     &                      mo,ng,ig,jg,
     &                      use_trs, tac, tca, tab, tba,
     &                      gen_name_in_carts )
c
      implicit real*8(a-h,o-z)
      integer igen(*),jgen(2,3,*),ig(*),jg(2,3,*)
      integer im(48,48),ir(3,3,48)
      integer iflg(48),jw(2,3,48),jww(2,3),mo(2,3),mo1(2,3),jw1(2,3)
      integer lprm(2,3,4)
      character*5 schoen(48)

      logical use_trs, gen_name_in_carts
      logical with_inverse
      integer iflg_tmp(48)
c --
      integer i, j, k
      integer m1, m2
      integer ir_sysA(3,3,48), ir_sysB(3,3,48)
      real*8 c1
      real*8 tac(3,3), tca(3,3), tab(3,3), tba(3,3)
c
      integer i_shift
c -
c -- make op matrix in the A system
      if ( gen_name_in_carts ) then
         Do k=1, 48
            Do i=1, 3
               Do j=1, 3
                  c1 = 0.0d0
                  Do m1=1, 3
                     Do m2=1, 3
                        c1 = c1 + tac(i,m1)*ir(m1,m2,k)*tca(m2,j)
                     End do
                  End do
                  ir_sysA(i,j,k) = nint(c1)
               End do
            End do
         End do
c -- make op matrix in the B system
         Do k=1, 48
            Do i=1, 3
               Do j=1, 3
                  c1 = 0.0d0
                  Do m1=1, 3
                     Do m2=1, 3
                        c1 = c1 + tab(i,m1)*ir_sysA(m1,m2,k)*tba(m2,j)
                     End do
                  End do
                  ir_sysB(i,j,k) = nint(c1)
               End do
            End do
         End do

      else
         ir_sysB = ir
      endif
c
c --- register the generator ---
c
      iflg(1)=1
      do j=1,3
        jw(1,j,1)=0
        jw(2,j,1)=1
      End do

      do i=2,48
         iflg(i)=0
      End do

      do i=1,ngen
         if(il.ge.1 .and. igen(i).gt.48) then
            write(6,*) ' igen=',igen(i),' > 48 : error'
            stop 'error:sub.nspgrp (igen)'
         end if
         if(il.le.0 .and. igen(i).gt.24) then
            write(6,*) ' igen=',igen(i),' > 24 : error'
            stop 'error:sub.nspgrp (igen)'
         end if
         iflg(igen(i))=1
         do j=1,3
            jw(1,j,igen(i))=jgen(1,j,i)
            jw(2,j,igen(i))=jgen(2,j,i)
         End do
      End do

      n=0
      do 18 i=1,48
   18 n=n+iflg(i)

      if(jpr.ge.0)      write(jf,*) ' '
      if(jpr.ge.0)      write(jf,*) ' ----- generator -------'
      if(jpr.ge.0) then
      nn=0
      do 80 i=1,48
        if(iflg(i).eq.1) then
          nn=nn+1
          write(jf,100) nn,i,schoen(i),((jw(j,k,i),j=1,2),k=1,3)
        end if
   80 continue
      end if
c

c --- produce group elements by multiplications
c
    1 ng=n
      do 20 i=1,48
      if(iflg(i).eq.1) then
        do 22 j=1,48
        if(iflg(j).eq.1) then
c         iw=im(i,j)
          iw=im(j,i)
          if(iflg(iw).ne.1) then

            iflg(iw)=1
            call nsrot1(jw(1,1,i),ir_sysB(1,1,j),jww)
            call nssum1(1,jww,jw(1,1,j),jw(1,1,iw))
          end if
        end if
   22   continue
      end if
   20 continue
      n=0
      do 24 i=1,48
   24 n=n+iflg(i)
      if(n.gt.ng) go to 1

c
c --- register group elements
c
      n=0
      do 30 i=1,48
      if(iflg(i).eq.1) then
        n=n+1
        ig(n)=i
        do 32 k=1,3
        do 32 j=1,2
   32   jg(j,k,n)=jw(j,k,i)
      end if
   30 continue
c
c --
c      write(*,*) 'use_trs now = ',use_trs
c
c
c ----------------------------- check if opr has inv-symm --
      with_inverse = .false.
      iflg_tmp = 0
c
      if ( il .ge. 1 ) then
         if ( iflg(25) .eq. 1 ) then
            with_inverse = .true.
         endif
      else if ( il .le. 0 ) then
         if ( iflg(13) .eq. 1 ) then
            with_inverse = .true.
         endif
      endif
c ---
      if ( (.not. with_inverse) .and. use_trs ) then
         iflg_tmp = 0

         if ( il. ge. 1 ) then
            Do i=1, 48
               if ( i .le. 24 ) then
                  i_shift = 24
               else
                  i_shift = -24
               endif
               if ( iflg(i) .eq. 1 ) then
                  iflg_tmp(i) = 1
                  iflg_tmp(i+i_shift) = 1
               endif
            End do
         else if ( il .le. 0 ) then
            Do i=1, 24
               if ( i .le. 12 ) then
                  i_shift = 12
               else
                  i_shift = -12
               endif
               if ( iflg(i) .eq. 1 ) then
                  iflg_tmp(i) = 1
                  iflg_tmp(i+i_shift) = 1
               endif
            End do
         endif

         n=0
         do i=1,48
            if ( iflg_tmp(i).eq.1 ) then
               n=n+1
               ig(n)=i
               if ( iflg(i).eq.1 ) then
                  Do k=1, 3
                     Do j=1, 2
                        jg(j,k,n)=jw(j,k,i)
                     End do
                  End do
               else
                  if ( il.ge.1 ) then
                     if ( i .le. 24 ) then
                        i_shift = 24
                     else
                        i_shift = -24
                     endif
                  else
                     if ( i .le. 12 ) then
                        i_shift = 12
                     else
                        i_shift = -12
                     endif
                  endif
                  Do k=1, 3
                     Do j=1, 2
                        jg(j,k,n)=jw(j,k,i+i_shift)
                     End do
                  End do
               endif
            endif
         end do

        iflg = iflg_tmp

      end if
c

c --- move the origin
c
      do 38 i=1,3
        mo(1,i)=0
        mo(2,i)=1
   38 continue
      if(inv.ne.0) then
        i1=0
        do 40 i=1,ng
          if(il.le.0 .and. ig(i).eq.13) i1=i
          if(il.ge.1 .and. ig(i).eq.25) i1=i
   40   continue
        if(i1.ne.0) then
          do 42 i=1,3
            mo(1,i)=jg(1,i,i1)
            mo(2,i)=jg(2,i,i1)*2
            call nsrduc(mo(1,i),mo(2,i))
   42     continue
          if(jpr.ge.1) write(6,200) ((mo(i,j),i=1,2),j=1,3)
          do 44 i=1,ng
            call nsrot1(mo,ir_sysB(1,1,ig(i)),mo1)
            call nssum1( 1,jg(1,1,i),mo1,jw1)
            call nssum1(-1,jw1,mo,jg(1,1,i))
   44     continue
        end if
      end if
c
      call nsprmv(il,nprm,lprm)
      call nssdjg(jf,jpr,nprm,lprm,ng,jg)
c
c --- print
c

      if(jpr.ge.2) then
        write(jf,*) ' '
        write(jf,*) ' ----- group elements -------'
        do 90 i=1,ng
   90   write(jf,100) i,ig(i),schoen(ig(i)),
     &                ((jg(j,k,i),j=1,2),k=1,3)
  100   format(i5,'   (',i2,')',3x,a5,5x,
     &         '(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )
      end if
  200   format(/'    move the origin by ','  (',2(i3,' /',i3,'  ,'),
     &                                         i3,' /',i3,' )'      )
c
      return
      end
c ===================================================================== 12.0A

c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nspgw1(jf,tca,tac,tab,tba,tcb,tbc,
     &                     grc,gkc,gra,gka,grb,gkb,schoen,
     &           il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1,
     &           euler1,inver1)
c
      implicit real*8(a-h,o-z)
      real*8 tca(3,3),tac(3,3),tab(3,3),tba(3,3),tcb(3,3),tbc(3,3)
      real*8 grc(3,3),gkc(3,3),gra(3,3),gka(3,3),grb(3,3),gkb(3,3)
      character*5 schoen(48)
      integer il,ng0,ng1,ig01(48),im1(48,48),iv1(48)
      integer lra1(3,3,48),lsa1(3,3,48),lrb1(3,3,48),lsb1(3,3,48)
      real*8  ta1(3,48),tb1(3,48)
      real*8  euler1(3,48)
      integer inver1(  48)
c
      write(jf,100) ((tca(i,j),i=1,3),j=1,3)
      write(jf,100) ((tac(i,j),i=1,3),j=1,3)
      write(jf,100) ((tab(i,j),i=1,3),j=1,3)
      write(jf,100) ((tba(i,j),i=1,3),j=1,3)
      write(jf,100) ((tcb(i,j),i=1,3),j=1,3)
      write(jf,100) ((tbc(i,j),i=1,3),j=1,3)
c
      write(jf,100) ((grc(i,j),i=1,3),j=1,3)
      write(jf,100) ((gkc(i,j),i=1,3),j=1,3)
      write(jf,100) ((gra(i,j),i=1,3),j=1,3)
      write(jf,100) ((gka(i,j),i=1,3),j=1,3)
      write(jf,100) ((grb(i,j),i=1,3),j=1,3)
      write(jf,100) ((gkb(i,j),i=1,3),j=1,3)
c
      write(jf,200) il,ng0
      write(jf,300) (schoen(i),i=1,ng0)
      write(jf,200) ng1
      write(jf,200) (ig01(i),i=1,ng1)
      write(jf,200) ((im1(i,j),i=1,ng1),j=1,ng1)
      write(jf,200) (iv1(i),i=1,ng1)
      write(jf,100) ((ta1(i,j),i=1,3),j=1,ng1)
      write(jf,220) (((lra1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      write(jf,220) (((lsa1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      write(jf,100) ((tb1(i,j),i=1,3),j=1,ng1)
      write(jf,220) (((lrb1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      write(jf,220) (((lsb1(i,j,k),i=1,3),j=1,3),k=1,ng1)
      write(jf,100) ((euler1(i,j),i=1,3),j=1,ng1)
      write(jf,200) ( inver1(i),i=1,ng1)
c
  100 format((3e24.16))
  200 format((24i3))
  220 format((18i3))
  300 format((10(2x,a5)))
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsprmv(il,nv,lv)
c
c#12  input   :       il : lattice type
c#12  output  :       nv : # of vectors
c#12           lv(2,3,4) : primitive translation vector
c#12  noexternal
c
c#21  to get primitive translation vectors
c     e.g., (0/1,0/1,0/1), (1/2,1/2,1/2) for body centered lattice
c
c#31  1990.01.06.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsprmv(il,nv,lv)
c
      implicit real*8(a-h,o-z)
      dimension lv(2,3,4)
c
      do 10 i=1,3
        lv(1,i,1)=0
        lv(2,i,1)=1
   10 continue
c
      if(il.eq.-1) then
        nv=3
        lv(1,1,2)=2
        lv(2,1,2)=3
        lv(1,2,2)=1
        lv(2,2,2)=3
        lv(1,3,2)=1
        lv(2,3,2)=3
c
        lv(1,1,3)=1
        lv(2,1,3)=3
        lv(1,2,3)=2
        lv(2,2,3)=3
        lv(1,3,3)=2
        lv(2,3,3)=3
      else if(il.eq.0 .or. il.eq.1) then
        nv=1
      else if(il.eq.2) then
        nv=4
        lv(1,1,2)=1
        lv(2,1,2)=2
        lv(1,2,2)=1
        lv(2,2,2)=2
        lv(1,3,2)=0
        lv(2,3,2)=1
c
        lv(1,1,3)=1
        lv(2,1,3)=2
        lv(1,2,3)=0
        lv(2,2,3)=1
        lv(1,3,3)=1
        lv(2,3,3)=2
c
        lv(1,1,4)=0
        lv(2,1,4)=1
        lv(1,2,4)=1
        lv(2,2,4)=2
        lv(1,3,4)=1
        lv(2,3,4)=2
      else if(il.eq.3) then
        nv=2
        lv(1,1,2)=1
        lv(2,1,2)=2
        lv(1,2,2)=1
        lv(2,2,2)=2
        lv(1,3,2)=1
        lv(2,3,2)=2
      else if(il.eq.4) then
        nv=2
        lv(1,1,2)=1
        lv(2,1,2)=2
        lv(1,2,2)=1
        lv(2,2,2)=2
        lv(1,3,2)=0
        lv(2,3,2)=1
      else if(il.eq.5) then
        nv=2
        lv(1,1,2)=0
        lv(2,1,2)=1
        lv(1,2,2)=1
        lv(2,2,2)=2
        lv(1,3,2)=1
        lv(2,3,2)=2
      else if(il.eq.6) then
        nv=2
        lv(1,1,2)=1
        lv(2,1,2)=2
        lv(1,2,2)=0
        lv(2,2,2)=1
        lv(1,3,2)=1
        lv(2,3,2)=2
      else
        write(6,*) 'il=',il
        write(6,*) '=== stop at sub.nsprmv. (il) ==='
        stop 'error:sub.nsprmv'
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrduc(in,id)
c
c#12  in-out put :    in : numerator
c#12                  id : denominator
c
c#21  to reduce a fraction in/id
c
c#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrduc(in,id)
c
      implicit real*8(a-h,o-z)
c
      if(in.eq.0) then
        id=1
      else
        call nsgcm2(in,id,m)
        in=in/m
        id=id/m
      end if
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrmxc(euler,r)
c
c#12  input :  euler(3) : euler angles (alpha,beta, gamma)
c#12  output     r(3,3) : rotation matrix in real
c#13  noexternal
c
c#21  to get rotation matrix for cubic operations
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrmxc(euler,r)
c
      implicit real*8(a-h,o-z)
      dimension euler(3),r(3,3)
c
c                        euler(1) :  alpha
c                        euler(2) :  beta
c                        euler(3) :  gamma
c
c                  point : x         u     x
c                          y  ---->  v = r y
c                          z    r    w     z
c
      ca=cos(euler(1))
      cb=cos(euler(2))
      cc=cos(euler(3))
      sa=sin(euler(1))
      sb=sin(euler(2))
      sc=sin(euler(3))
      r(1,1) = ca*cb*cc - sa*   sc
      r(2,1) = sa*cb*cc + ca*   sc
      r(3,1) =          -    sb*cc
      r(1,2) =-ca*cb*sc - sa*   cc
      r(2,2) =-sa*cb*sc + ca*   cc
      r(3,2) =               sb*sc
      r(1,3) =            ca*sb
      r(2,3) =            sa*sb
      r(3,3) =               cb
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrmxh(euler,r)
c
c#12  input :  euler(3) : euler angles (alpha,beta, gamma)
c#12  output     r(3,3) : rotation matrix for hexagonal
c#13  external : sub.nsmlt1
c
c#21  to get a rotation matrix for a hexagonal operation
c
c#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrmxh(euler,r)
c
      implicit real*8(a-h,o-z)
      dimension euler(3),r(3,3),rr(3,3,3),rrr(3,3)
c
c                        euler(1) :  alpha
c                        euler(2) :  beta
c                        euler(3) :  gamma
c
c                  point : x         u     x
c                          y  ---->  v = r y
c                          z    r    w     z
c
      three=3
      ca=cos(euler(1))
      cb=cos(euler(2))
      cc=cos(euler(3))
      sa=sin(euler(1))
      sb=sin(euler(2))
      sc=sin(euler(3))
      p2=2
      p3=sqrt(three)
      rr(1,1,1) =  ca + sa/p3
      rr(2,1,1) =  p2*sa/p3
      rr(3,1,1) =  0
      rr(1,2,1) = -rr(2,1,1)
      rr(2,2,1) =  ca -sa/p3
      rr(3,2,1) =  0
      rr(1,3,1) =  0
      rr(2,3,1) =  0
      rr(3,3,1) =  1
c
      rr(1,1,2) =  cb
      rr(2,1,2) =  0
      rr(3,1,2) = -sb
      rr(1,2,2) = (1-cb)/p2
      rr(2,2,2) =  1
      rr(3,2,2) =  sb/p2
      rr(1,3,2) =  sb
      rr(2,3,2) =  0
      rr(3,3,2) =  cb
c
      rr(1,1,3) =  cc+sc/p3
      rr(2,1,3) =  p2*sc/p3
      rr(3,1,3) =  0
      rr(1,2,3) = -rr(2,1,3)
      rr(2,2,3) =  cc-sc/p3
      rr(3,2,3) =  0
      rr(1,3,3) =  0
      rr(2,3,3) =  0
      rr(3,3,3) =  1
c
      call nsmlt1(3,3,rr(1,1,1),rr(1,1,2),rrr)
      call nsmlt1(3,3,rrr,rr(1,1,3),r)
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrot1(ja,ir,jb)
c
c#12  input   : ja(2,3)  : a vector
c#12            ir(3,3)  : rotation matrix
c#12  output  : jb(2,3)  : a vector rotated
c#13  external: nsgcm3
c
c#21  to get a rotated vector
c#21  (jb(1,1)/jb(2,1), jb(1,2)/jb(2,2), jb(1,3)/jb(2,3))
c
c#31  1989.12.28.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrot1(ja,ir,jb)
c
      implicit real*8(a-h,o-z)
      integer ja(2,3),ir(3,3),jb(2,3)
c
      call nsgcm3(ja(2,1),ja(2,2),ja(2,3),m)
      mm=m*m
      id=ja(2,1)*ja(2,2)*ja(2,3)/mm
      i1=ja(1,1)*(id/ja(2,1))
      i2=ja(1,2)*(id/ja(2,2))
      i3=ja(1,3)*(id/ja(2,3))
c
c     ja=(i1,i2,i3)/id
c
      jb(1,1)=ir(1,1)*i1+ir(1,2)*i2+ir(1,3)*i3
      jb(1,2)=ir(2,1)*i1+ir(2,2)*i2+ir(2,3)*i3
      jb(1,3)=ir(3,1)*i1+ir(3,2)*i2+ir(3,3)*i3
      jb(2,1)=id
      jb(2,2)=id
      jb(2,3)=id
c
      call nsrduc(jb(1,1),jb(2,1))
      call nsrduc(jb(1,2),jb(2,2))
      call nsrduc(jb(1,3),jb(2,3))
      jb(1,1)=mod(jb(1,1),jb(2,1))
      jb(1,2)=mod(jb(1,2),jb(2,2))
      jb(1,3)=mod(jb(1,3),jb(2,3))
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrota(natm,lmnatm,vatm,iatm21,ratm21)
c
      implicit real*8(a-h,o-z)
      common/spg2  / il,ng,ig0(48),jv(2,3,48)
      dimension vatm(3,lmnatm),iatm21(lmnatm,48),ratm21(3,lmnatm,48)
      dimension v(3),r(3)
c
      do 10 ig1=1,ng
      do 10 j1=1,natm
        call nsrotr(ig1,vatm(1,j1),v)
        do 20 j2=1,natm
          jj2=j2
          do 22 i=1,3
   22     r(i)=vatm(i,j2)-v(i)
          call nslatr(il,r,ind)
        if(ind.eq.0) go to 24
   20   continue
        stop ' === stop in sub.nsrota. (no atom) ==='
   24   continue
        iatm21(j1,ig1)=jj2
        do 26 i=1,3
   26   ratm21(i,j1,ig1)=r(i)
   10 continue
c
      return
      end
c
c
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
c
c#12  noinput
c#12  output: ng0                 : # of group elements
c#12          ieuler(3,24)*2*pai/4: euler angles for cubic
c#12           euler(3,24)        : euler angles (alpha,beta, gamma)
c#12          irot(3,3,48)        : rotation matrix in integer
c#12           rot(3,3,48)        : rotation matrix in real
c#12          ir1234(3,48)        : jones faithful representation
c#12           jones(3,48)   (a2) : jones faithful representation
c#12          schoen(3,48)   (a5) : schoenflies notation
c#13  external : sub.nseulc, nsrmxc
c
c#21  to get euler angles, rotation matrix,
c#21         jones faithfull representation, and
c#21         schoenflies notation
c#21         for point-group operations
c#21         of cubic lattice
c
c#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
c
      implicit real*8(a-h,o-z)
      character*2 jones(3,48)
      integer ieuler(3,24),irot(3,3,48),ir1234(3,48)
      real*8  euler(3,24), rot(3,3,48), r(3,3)
      character*1 cha0(3)
      character*5 schoen(48)
      character*4 cmn(24)
      data cmn            /'e   ','c2x ','c2y ','c2z ',
     &                     'c31+','c32+','c33+','c34+',
     &                     'c31-','c32-','c33-','c34-',
     &       'c2a ','c2b ','c2c ','c2d ','c2e ','c2f ',
     &       'c4x+','c4y+','c4z+','c4x-','c4y-','c4z-'/
      data cha0 /'x','y','z'/
      data one /1.00001/
c
      ng0=48
      call nseulc(ieuler,euler)
      do 10 i=1,24
        call nsrmxc(euler(1,i),r)
        do 12 j=1,3
        do 12 k=1,3
          irot(j,k,i)=int(one*r(j,k))
          rot(j,k,i)=r(j,k)
   12   continue
        do 20 j=1,3
        do 20 k=1,3
          irr=irot(j,k,i)
          if(irr.eq.-1) then
              ir1234(j,i)=-k
          else if(irr.eq.1) then
              ir1234(j,i)= k
          end if
   20   continue
   10 continue
c
      do 30 i=1,24
        schoen(i   )=' '//cmn(i)
        schoen(i+24)='i'//cmn(i)
      do 30 j=1,3
        ir1234(j,i+24)=-ir1234(j,i)
      do 30 k=1,3
        irot(j,k,i+24)=-irot(j,k,i)
         rot(j,k,i+24)=- rot(j,k,i)
   30 continue
c
      do 40 i=1,48
      do 40 j=1,3
        k=ir1234(j,i)
        if(k.lt.0) then
        jones(j,i)='-'//cha0(iabs(k))
      else
        jones(j,i)=' '//cha0(iabs(k))
      end if
   40 continue
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
c
c#12  noinput
c#12  output: ng0                 : # of group elements
c#12          ieuler(3,24)*2*pai/4: euler angles for hexagonal
c#12           euler(3,24)        : euler angles (alpha,beta, gamma)
c#12          irot(3,3,48)        : rotation matrix in integer
c#12           rot(3,3,48)        : rotation matrix in real
c#12          ir1234(3,48)        : jones faithful representation
c#12           jones(3,48)   (a2) : jones faithful representation
c#12          schoen(3,48)   (a5) : schoenflies notation
c#13  external : sub.nseulh, nsrmxh, nsjonh, nstrsh, nsmlt1
c
c#21  to get euler angles, rotation matrix,
c#21         jones faithfull representation, and
c#21         schoenflies notation
c#21         for point-group operations
c#21         of hexagonal lattice
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
c
      implicit real*8(a-h,o-z)
      character*2 jones(3,48)
      character*5 schoen(48)
      integer ieuler(3,24),irot(3,3,48),ir1234(3,48),iw(3,3)
      real*8  euler(3,24), rot(3,3,48), r(3,3), u(3,3), v(3,3)
      character*4 hmn(12)
      data        hmn    /'e   ','c6+ ','c3+ ','c2  ','c3- ','c6- ',
     &                    'c211','c221','c231','c212','c222','c232'/
      data onep /1.00001/
c
      ng0=24
      call nseulh(ieuler,euler)
      do 10 i=1,12
        call nsrmxh(euler(1,i),r)
        do 12 k=1,3
        do 12 j=1,3
          rot(j,k,i)=r(j,k)
          irot(j,k,i)=int(onep*r(j,k))
          rot(j,k,i+12)=-r(j,k)
          irot(j,k,i+12)=-int(onep*r(j,k))
   12   continue
   10 continue
c
      call nsjonh(irot,ir1234,jones)
c
      do 30 i=1,12
        schoen(i   )=' '//hmn(i)
        schoen(i+12)='i'//hmn(i)
   30 continue
c
      call nstrsh(u,v)
      do 40 i=1,24
        call nsmlt1(3,3,v,rot(1,1,i),r)
        call nsmlt1(3,3,r,u,rot(1,1,i+24))
   40 continue
      do 50 i=25,48
      do 50 k=1,3
      do 50 j=1,3
        irot(j,k,i)=int(onep*rot(j,k,i))
   50 continue
      call nsjonh(irot(1,1,25),ir1234(1,25),jones(1,25))
c
      do 60 i=25,48
        do 62 k=1,3
        do 62 j=1,3
   62   iw(j,k)=irot(k,j,i)
        do 64 k=1,3
        do 64 j=1,3
          irot(j,k,i)=iw(j,k)
           rot(j,k,i)=iw(j,k)
   64   continue
   60 continue
c
c     check
c
      do 70 i=1,24
        do 72 j=1,3
        do 72 k=1,3
          i1=0
          i2=0
          do 74 l=1,3
            i1=i1+irot(j,l,i)*irot(l,k,i+24)
            i2=i2+irot(j,l,i+24)*irot(l,k,i)
   74     continue
          if((j.eq.k .and. (i1.ne.1 .or. i2.ne.1)) .or.
     &       (j.ne.k .and. (i1.ne.0 .or. i2.ne.0))     ) then
            write(6,*) ' j,k=',j,k,'   i1,i2=',i1,i2,'  : i=',i
            write(6,*) ' === stop in sub.nsroth. ',
     &                 ' (inverse matrix condition) ==='
            stop 'error in sub.nsroth'
          end if
   72   continue
   70 continue
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrotk(ig,k1,k2)
c
      implicit real*8(a-h,o-z)
      common/nspg0 / euler(3,24), rot(3,3,48),
     &              ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &               im0(48,48),iv0(48),ng0
      common/spg2  / il,ng,ig0(48),jv(2,3,48)
      integer k1(4),k2(4)
c
      if(il.le.0) then
        igg=ig0(ig)+24
      else
        igg=ig0(ig)
      end if
      k2(4)=k1(4)
      do 10 i=1,3
        k=0
        do 12 j=1,3
        k=k+irot(i,j,igg)*k1(j)
   12   continue
        k2(i)=k
   10 continue
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nsrotl(jf,jpr,il,ng0,ieuler,euler,irot,rot,
c#11             ir1234,jones,schoen                 )
c
c#12  input: jf    : output fileil=<0 for hexagonal, il>0 for cubic.
c#12         jpr   : output control
c#12         il    : lattice typei: il=<0 for hexagonal, il>0 for cubic.
c
c#12  output: ng0                 : # of group elements
c#12          ieuler(3,48)*2*pai/6: euler angles for hexagonal
c#12          ieuler(3,48)*2*pai/4: euler angles for cubic
c#12           euler(3,48)        : euler angles (alpha,beta, gamma)
c#12          irot(3,3,48)        : rotation matrix in integer
c#12           rot(3,3,48)        : rotation matrix in real
c#12           jones(3,48)   (a2) : jones faithful representaion
c#12          schoen(3,48)   (a5) : schoenflies notation
c#13  external : sub.nsrotc, nsroth
c
c#21  to get euler angles, rotation matrix,
c#21         jones faithfull representation, and
c#21         schoenflies notation
c#21         for point-group operations
c#21         of hexagonal and cubic lattices
c
c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrotl(il,jpr,jf,ng0,ieuler,euler,irot,rot,
     &                  ir1234,jones,schoen)
c
      implicit real*8(a-h,o-z)
      character*2 jones(3,48)
      character*5 schoen(48)
      integer ieuler(3,24),irot(3,3,48),ir1234(3,48)
      real*8  euler(3,24), rot(3,3,48)
      if(il.le.0) then
        call nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
      else
        call nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
      end if
      no=ng0
      if(jpr.ge.2) then
      write(jf,*) '     '
      if(il.le.0) then
        write(jf,*) ' table of operation code for the point group',
     &              ' 6/mmm(d6h) '
      else
        write(jf,*) ' table of operation code for the point group',
     &              ' m3m(oh) '
      end if
      if(il.le.0) write(jf,*) ' (for real-space coordinate, w=x-y)'
      write(jf,100)
        do 10 i=1,no/2
          write(jf,120) i,schoen(i),(jones(j,i),j=1,3),
     &                  (ir1234(j,i),j=1,3),(euler(j,i),j=1,3)
   10   continue
        do 12 i=no/2+1,no
          write(jf,140) i,schoen(i),(jones(j,i),j=1,3),
     &                  (ir1234(j,i),j=1,3)
   12   continue
      if(il.le.0) then
        write(jf,*) '    '
        write(jf,*) ' (for reciprocal-space coordinate, w=x+y) '
        do 14 i=25,48
          write(jf,160) i,          (jones(j,i),j=1,3),
     &                  (ir1234(j,i),j=1,3)
   14   continue
      end if
      end if
  100 format(1h ,' no.',' schoenflies   jones       ir1234    ',
     &           '   euler-angle(alpha,beta,gamma)')
  120 format(1h ,'(',i2,')',5x,a5,' (',3(1x,a2),' )  (',3(1x,i2),' )',
     &           5x,'(',3f8.3,' )')
  140 format(1h ,'(',i2,')',5x,a5,' (',3(1x,a2),' )  (',3(1x,i2),' )')
  160 format(1h ,'(',i2,')',5x,5x,' (',3(1x,a2),' )  (',3(1x,i2),' )')
c
      if(jpr.ge.3) then
      write(jf,*) '   '
      write(jf,*) ' matrix representation of operation '
      if(il.le.0) write(jf,*) ' (for real-space coordinate)'
      do 20 ii=1,no,6
      write(jf,*) '   '
      write(jf,200) (i,(irot(1,k,i),k=1,3),i=ii,ii+5)
      do 22 j=2,3
   22 write(jf,220) (  (irot(j,k,i),k=1,3),i=ii,ii+5)
   20 continue
      if(il.le.0) then
      write(jf,*) '   '
      write(jf,*) ' matrix representation of operation '
      write(jf,*) ' (for reciprocal-space coordinate)'
      do 30 ii=25,48,6
      write(jf,*) '   '
      write(jf,200) (i,(irot(1,k,i),k=1,3),i=ii,ii+5)
      do 32 j=2,3
   32 write(jf,220) (  (irot(j,k,i),k=1,3),i=ii,ii+5)
   30 continue
      end if
      end if
  200 format(1h ,6('(',i2,')',3i2,2x))
  220 format(1h ,6(4x        ,3i2,2x))
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsrotr(ig,r1,r2)
c
      implicit real*8(a-h,o-z)
      common/nspg0 / euler(3,24), rot(3,3,48),
     &              ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &               im0(48,48),iv0(48),ng0
      common/spg2  / il,ng,ig0(48),jv(2,3,48)
      dimension r1(3),r2(3)
c
      do 10 i=1,3
        a=0.0
        do 12 j=1,3
        a=a+rot(i,j,ig0(ig))*r1(j)
   12   continue
        r2(i)=a+dfloat(jv(1,i,ig))/jv(2,i,ig)
   10 continue
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nss0(es,ef,nxx,nyy,nzz,ev,ew,s0)
c
c     ***  by linear interporation in tetrahedrons
c     ***
c     ***
c     ***  nxx  number of mesh points in x-direction
c     ***  nyy  number of mesh points in y-direction
c     ***  nzz  number of mesh points in z-direction
c     ***
c     ***  s0   susceptibility with a constant matrix element(=one)
c     ***  es   energy argument for s0
c     ***  ef   fermi energy
c     ***  ev   energy band (expected to be   occupied)
c     ***  ew   energy band (expected to be unoccupied)
c     ***
c
      implicit  real*8(a-h,o-z)
c
      real*8  ev(*),ew(*)
      real*8  evcub(2,2,2), evc(8), evt(4), evb(4)
      real*8  ewcub(2,2,2), ewc(8), ewt(4), ewb(4)
      integer iecub(2,2,2), iec(8), iet(4), ieb(4)
      equivalence(evc(1),evcub(1,1,1))
      equivalence(ewc(1),ewcub(1,1,1))
      equivalence(iec(1),iecub(1,1,1))
      integer iqmat(6,2)
      data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
c
      npx=nxx+1
      npy=nyy+1
      npz=nzz+1
      np=npx*npy*npz
      ncub=nxx*nyy*nzz
      ntet=6*ncub
c
      s0=0
c
      emin0= 1.e30
      emax0=-1.e30
      do 12 ip=1,np
        if(ev(ip).lt.emin0) emin0=ev(ip)
        if(ew(ip).gt.emax0) emax0=ew(ip)
   12 continue
c                                   ============ if 1 ==
      if(ef.gt.emin0 .and. ef.lt.emax0) then
c
c     ***  integration over b.z. starts    ***
c
c     ***       sampling over cubes        ***
c
      icub=0
      do 20 iz=0,nzz-1
      do 20 iy=0,nyy-1
      do 20 ix=0,nxx-1
        icub=icub+1
c     ***  energies at cube corners  ***
        ni=npx*(npy*iz+iy)+ix
        emax=-1.e30
        emin= 1.e30
        do 30 kz=1,2
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
        evcub(kx,ky,kz)=ev(ip0)
        ewcub(kx,ky,kz)=ew(ip0)
        iecub(kx,ky,kz)=ip0
        if(ew(ip0).gt.emax) emax=ew(ip0)
        if(ev(ip0).lt.emin) emin=ev(ip0)
   30   continue
c                                   ============ if 2 ==
        if(ef.gt.emin .and. ef.lt.emax) then
c         ***      six tetrahedrons      ***
c         *** sampling over tetrahedrons ***
          evt(1)=evc(1)
          evt(4)=evc(8)
          ewt(1)=ewc(1)
          ewt(4)=ewc(8)
          iet(1)=iec(1)
          iet(4)=iec(8)
          do 40 it=1,6
            do 42 ip=1,2
              iq=iqmat(it,ip)
              evt(ip+1)=evc(iq)
              ewt(ip+1)=ewc(iq)
              iet(ip+1)=iec(iq)
   42       continue
            do 44 m=1,4
              evb(m)=evt(m)
              ewb(m)=ewt(m)
              ieb(m)=iet(m)
   44       continue
c
            call nss02r(es,ef,evb,ewb,s)
            s0=s0+s
c
   40     continue
        end if
c                                   ============ if 2 ==
   20 continue
c
      s0=s0/ntet
c
      end if
c                                   ============ if 1 ==
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nss00r(ed,s0)
      implicit real*8 (a-h,o-z)
      real*8 ed(4), r(4)
      e=0
      e1=ed(1)
      e2=ed(2)
      e3=ed(3)
      e4=ed(4)
      call nstts1(e1,e2,e3,e4)
      call nstt0r(e,e1,e2,e3,e4,r)
      s0=-r(1)-r(2)-r(3)-r(4)
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nss01r(es0,ef,eb,eq,s0)
c
c     susceptibility (constant matrix element)
c     by linear interporation in a tetrahedron
c     es : energy argument s0
c     ef : fermi energy for eb and eq
c     eb : occupied
c     eq : expected to be unoccupied
c     ref) J.Rath & A.J.Freeman, Phys.Rev.B11(1975)2109.
c
      implicit  real*8(a-h,o-z)
c
      real*8  eb(4), ec(4), ed(4), eq(4), er(4), es(4)
      integer iq(4), id(4)
c
      e=ef
      if(eb(1).gt.eb(2) .or. eb(2).gt.eb(3) .or. eb(3).gt.eb(4) .or.
     &   eb(4).gt.e) then
        write(6,*) ' eb=',eb,'    e=',e
        stop 'error(eb,e): sub.nss01r'
      end if
      eps=1.e-6
      do 10 i=1,4
        er(i)=eq(i)
        iq(i)=i
   10 continue
c     ***  er(1).le.er(2).le.er(3).le.er(4)  ***
      call nsttod(er,iq)
c
      e1=er(1)
      e2=er(2)
      e3=er(3)
      e4=er(4)
      do 12 i=1,4
   12 ec(i)=er(i)-eb(iq(i))-es0
c
      if(e.ge.e4) then
        s0=0
c       write(6,*) ' 0 '
      else if(e.ge.e3) then
        es(1)=e
        es(2)=e
        es(3)=e
        es(4)=e4
        dd=e4-e
        d41=e4-e1
        d42=e4-e2
        d43=e4-e3
        vv1=(dd**3)/(d41*d42*d43)
        ed(1)=(ec(1)-ec(4))*dd/d41+ec(4)
        ed(2)=(ec(2)-ec(4))*dd/d42+ec(4)
        ed(3)=(ec(3)-ec(4))*dd/d42+ec(4)
        ed(4)= ec(4)
        call nsttod(ed,id)
        call nss00r(ed,s01)
        s0=vv1*s01
c       write(6,*) ' 1 : vv1=',vv1
c
      else if(e.ge.e2) then
c
        es(1)=e
        es(2)=e
        es(3)=e3
        es(4)=e4
        dd=e4-e
        d41=e4-e1
        d42=e4-e2
        vv1=(dd**2)/(d41*d42)
        ed(1)=(ec(1)-ec(4))*dd/d41+ec(4)
        ed(2)=(ec(2)-ec(4))*dd/d42+ec(4)
        ed(3)= ec(3)
        ed(4)= ec(4)
        call nsttod(ed,id)
        call nss00r(ed,s01)
c
        es(1)=e
        es(2)=e
        es(3)=e
        es(4)=e3
        dd =e4-e
        d41=e4-e1
        d42=e4-e2
        d32=e3-e2
        vv2=(e4-es(1))*(es(2)-e2)*(e3-es(3))/(d41*d42*d32)
        ed(1)=(ec(1)-ec(4))*dd/d41+ec(4)
        ed(2)=(ec(2)-ec(4))*dd/d42+ec(4)
        ed(3)=(ec(2)-ec(3))*(e3-es(3))/d32+ec(3)
        ed(4)= ec(3)
        call nsttod(ed,id)
        call nss00r(ed,s02)
c
        es(1)=e
        es(2)=e
        es(3)=e
        es(4)=e3
        d41=e4-e1
        d32=e3-e2
        d31=e3-e1
        vv3=(es(1)-e1)*(e3-es(3))*(e3-es(2))/(d41*d32*d31)
        ed(1)=(ec(1)-ec(4))*(e4-es(1))/d41+ec(4)
        ed(2)=(ec(1)-ec(3))*(e3-es(2))/d31+ec(3)
        ed(3)=(ec(2)-ec(3))*(e3-es(3))/d32+ec(3)
        ed(4)= ec(3)
        call nsttod(ed,id)
        call nss00r(ed,s03)
c
        s0=s01*vv1+s02*vv2+s03*vv3
c       write(6,*) ' 2 : vv1=',vv1,'   vv2=',vv2,'   vv3=',vv3
c
      else if(e.ge.e1) then
c
        es(1)=e
        es(2)=e2
        es(3)=e3
        es(4)=e4
        dd=e4-e
        d41=e4-e1
        d42=e4-e2
        vv1=dd/d41
        ed(1)=(ec(1)-ec(4))*dd/d41+ec(4)
        ed(2)= ec(2)
        ed(3)= ec(3)
        ed(4)= ec(4)
        call nsttod(ed,id)
        call nss00r(ed,s01)
c
        es(1)=e
        es(2)=e
        es(3)=e2
        es(4)=e3
        d31=e3-e1
        d41=e4-e1
        d42=e4-e2
        vv2=(e3-es(1))*(es(2)-e1)/(d31*d41)
        ed(1)=(ec(1)-ec(3))*(e3-es(1))/d31+ec(3)
        ed(2)=(ec(1)-ec(4))*(e4-es(2))/d41+ec(4)
        ed(3)= ec(2)
        ed(4)= ec(3)
        call nsttod(ed,id)
        call nss00r(ed,s02)
c
        es(1)=e
        es(2)=e
        es(3)=e
        es(4)=e2
        d21=e2-e1
        d31=e3-e1
        d41=e4-e1
        vv3=(es(1)-e1)*(es(3)-e1)*(e2-es(2))/(d41*d31*d21)
        ed(1)=(ec(1)-ec(4))*(e4-es(1))/d41+ec(4)
        ed(2)=(ec(1)-ec(2))*(e2-es(2))/d21+ec(2)
        ed(3)=(ec(1)-ec(3))*(e3-es(3))/d31+ec(3)
        ed(4)= ec(2)
        call nsttod(ed,id)
        call nss00r(ed,s03)
c
        s0=s01*vv1+s02*vv2+s03*vv3
c       write(6,*) ' 3 : vv1=',vv1,'   vv2=',vv2,'   vv3=',vv3
c
      else
        es(1)=er(1)
        es(2)=er(2)
        es(3)=er(3)
        es(4)=er(4)
        ed(1)=ec(1)
        ed(2)=ec(2)
        ed(3)=ec(3)
        ed(4)=ec(4)
        call nsttod(ed,id)
        call nss00r(ed,s01)
        s0=s01
c       write(6,*) ' 4 '
c
      end if
  200 format(4f9.5)
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nss02r(es,ef,ev,ew,s0)
c
c     susceptibility (constant matrix element)
c     by linear interporation in a tetrahedron
c     es : energy argument of s0(susceptibility)
c     ef : fermi level for ev and ew
c     ev : expected to be   occupied
c     ew : expected to be unoccupied
c     ref) J.Rath & A.J.Freeman, Phys.Rev.B11(1975)2109.
c
      implicit  real*8(a-h,o-z)
c
      real*8  ev(4), ew(4), ea(4), eb(4), ep(4), eq(4)
      integer iv(4)
c
      e=ef
      do 10 i=1,4
        ea(i)=ev(i)
        iv(i)=i
   10 continue
c     ***  ea(1).le.ea(2).le.ea(3).le.ea(4)  ***
      call nsttod(ea,iv)
c
      e1=ea(1)
      e2=ea(2)
      e3=ea(3)
      e4=ea(4)
      do 12 i=1,4
   12 ep(i)=ew(iv(i))
c
      if(e.le.e1) then
        s0=0
c       write(6,*) ' 0 '
      else if(e.le.e2) then
        eb(1)=ea(1)
        eb(2)=e
        eb(3)=e
        eb(4)=e
        dd=e-ea(1)
        d21=ea(2)-ea(1)
        d31=ea(3)-ea(1)
        d41=ea(4)-ea(1)
        vv1=(dd**3)/(d21*d31*d41)
        eq(1)=ep(1)
        eq(2)=(ep(2)-ep(1))*dd/d21+ep(1)
        eq(3)=(ep(3)-ep(1))*dd/d31+ep(1)
        eq(4)=(ep(4)-ep(1))*dd/d41+ep(1)
        call nss01r(es,e,eb,eq,s01)
c       write(6,200) eb
c       write(6,200) eq
        s0=vv1*s01
c       write(6,*) ' 1 : vv1=',vv1
c
      else if(e.le.e3) then
        eb(1)=ea(1)
        eb(2)=ea(2)
        eb(3)=e
        eb(4)=e
        dd=e-ea(1)
        d31=ea(3)-ea(1)
        d41=ea(4)-ea(1)
        vv1=dd*dd/(d31*d41)
        eq(1)= ep(1)
        eq(2)= ep(2)
        eq(3)=(ep(3)-ep(1))*dd/d31+ep(1)
        eq(4)=(ep(4)-ep(1))*dd/d41+ep(1)
        call nss01r(es,e,eb,eq,s01)
c       write(6,200) eb
c       write(6,200) eq
c
        eb(1)=ea(2)
        eb(2)=e
        eb(3)=e
        eb(4)=e
        d31=ea(3)-ea(1)
        d32=ea(3)-ea(2)
        d41=ea(4)-ea(1)
        vv2=(ea(3)-eb(2))*(eb(3)-ea(2))*(eb(4)-ea(1))/(d31*d32*d41)
        eq(1)= ep(2)
        eq(2)=(ep(3)-ep(1))*(eb(2)-ea(1))/d31+ep(1)
        eq(3)=(ep(3)-ep(2))*(eb(3)-ea(2))/d32+ep(2)
        eq(4)=(ep(4)-ep(1))*(eb(4)-ea(1))/d41+ep(1)
        call nss01r(es,e,eb,eq,s02)
c       write(6,200) eb
c       write(6,200) eq
c
        eb(1)=ea(2)
        eb(2)=e
        eb(3)=e
        eb(4)=e
        d41=ea(4)-ea(1)
        d32=ea(3)-ea(2)
        d42=ea(4)-ea(2)
        vv3=(ea(4)-eb(2))*(eb(3)-ea(2))*(eb(4)-ea(2))/(d41*d32*d42)
        eq(1)= ep(2)
        eq(2)=(ep(4)-ep(1))*(eb(2)-ea(1))/d41+ep(1)
        eq(3)=(ep(3)-ep(2))*(eb(3)-ea(2))/d32+ep(2)
        eq(4)=(ep(4)-ep(2))*(eb(4)-ea(2))/d42+ep(2)
        call nss01r(es,e,eb,eq,s03)
c       write(6,200) eb
c       write(6,200) eq
c
        s0=s01*vv1+s02*vv2+s03*vv3
c       write(6,*) ' 2 : vv1=',vv1,'   vv2=',vv2,'   vv3=',vv3
c
      else if(e.lt.e4) then
        eb(1)=ea(1)
        eb(2)=ea(2)
        eb(3)=ea(3)
        eb(4)=e
        d41=ea(4)-ea(1)
        vv1=(eb(4)-ea(1))/d41
        eq(1)= ep(1)
        eq(2)= ep(2)
        eq(3)= ep(3)
        eq(4)=(ep(4)-ep(1))*(eb(4)-ea(1))/d41+ep(1)
        call nss01r(es,e,eb,eq,s01)
c       write(6,200) eb
c       write(6,200) eq
c
        eb(1)=ea(2)
        eb(2)=ea(3)
        eb(3)=e
        eb(4)=e
        d41=ea(4)-ea(1)
        d42=ea(4)-ea(2)
        vv2=(ea(4)-eb(3))*(eb(4)-ea(2))/(d41*d42)
        eq(1)= ep(2)
        eq(2)= ep(3)
        eq(3)=(ep(4)-ep(1))*(eb(3)-ea(1))/d41+ep(1)
        eq(4)=(ep(4)-ep(2))*(eb(4)-ea(2))/d42+ep(2)
        call nss01r(es,e,eb,eq,s02)
c       write(6,200) eb
c       write(6,200) eq
c
        eb(1)=ea(3)
        eb(2)=e
        eb(3)=e
        eb(4)=e
        d41=ea(4)-ea(1)
        d42=ea(4)-ea(2)
        d43=ea(4)-ea(3)
        vv3=(ea(4)-eb(3))*(ea(4)-eb(4))*(eb(2)-ea(3))/(d41*d42*d43)
        eq(1)= ep(3)
        eq(2)=(ep(4)-ep(3))*(eb(2)-ea(3))/d43+ep(3)
        eq(3)=(ep(4)-ep(1))*(eb(3)-ea(1))/d41+ep(1)
        eq(4)=(ep(4)-ep(2))*(eb(4)-ea(2))/d42+ep(2)
        call nss01r(es,e,eb,eq,s03)
        s0=s01*vv1+s02*vv2+s03*vv3
c       write(6,200) eb
c       write(6,200) eq
c       write(6,*) ' 3 : vv1=',vv1,'   vv2=',vv2,'   vv3=',vv3
c
      else
        eb(1)=ea(1)
        eb(2)=ea(2)
        eb(3)=ea(3)
        eb(4)=ea(4)
        eq(1)=ep(1)
        eq(2)=ep(2)
        eq(3)=ep(3)
        eq(4)=ep(4)
        call nss01r(es,e,eb,eq,s01)
        s0=s01
c       write(6,200) eb
c       write(6,200) eq
c       write(6,*) ' 4 '
c
      end if
  200 format(4f9.5)
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nssdjg(jf,jpr,nv,lv,ng,jg)
c
c#12  input     :        jf : output file
c#12                    jpr : print control
c#12                     nv : # of lv vectors
c#12              lv(2,3,4) : primitive translation vector
c#12                     ng : oder of group
c#12  in-output : jg(2,3,48): translation vector
c#12  noexternal
c
c#21  to imodify nonprimitive translation vectors
c#21  (positive, smaller elements)
c
c#31  1990.01.06.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nssdjg(jf,jpr,nv,lv,ng,jg)
c
      implicit real*8(a-h,o-z)
      integer jg(2,3,48),lv(2,3,4)
      integer ja(2,3,48),jb(2,3),jc(2,3)
c
      do 10 j=1,3
      do 10 k=1,ng
        jg(1,j,k)=mod(jg(1,j,k),jg(2,j,k))
        if(jg(1,j,k).lt.0) jg(1,j,k)=jg(1,j,k)+jg(2,j,k)
   10 continue
c
      do 20 j=1,3
        ja(1,j,1)=0
        ja(2,j,1)=1
   20 continue
      na=1
c
      do 30 m=1,ng
c
        do 50 ia=1,na
        do 50  k=1,nv
          call nssum1(1,jg(1,1,m),lv(1,1,k),jb)
          if(jb(1,1).eq.ja(1,1,ia) .and.
     &       jb(2,1).eq.ja(2,1,ia) .and.
     &       jb(1,2).eq.ja(1,2,ia) .and.
     &       jb(2,2).eq.ja(2,2,ia) .and.
     &       jb(1,3).eq.ja(1,3,ia) .and.
     &       jb(2,3).eq.ja(2,3,ia)      ) then
            do 52 i=1,2
            do 52 j=1,3
   52       jg(i,j,m)=ja(i,j,ia)
            go to 30
          end if
   50   continue
c
c       registration
c
        na=na+1
        do 60 i=1,2
        do 60 j=1,3
   60   jb(i,j)=jg(i,j,m)
        do 70 k=2,nv
          call nssum1(1,jg(1,1,m),lv(1,1,k),jc)
          if(jb(1,1).gt.jc(1,1) .or.
     &       jb(1,2).gt.jc(1,2) .or.
     &       jb(1,3).gt.jc(1,3)     ) then
            do 72 j=1,3
              jb(1,j)=jc(1,j)
              jb(2,j)=jc(2,j)
   72       continue
          end if
   70   continue
        do 80 i=1,2
        do 80 j=1,3
          jg(i,j,m) =jb(i,j)
          ja(i,j,na)=jb(i,j)
   80   continue
c
   30 continue
c
      if(jpr.gt.1) then
        write(jf,*) ' '
        write(jf,*) '= nonprimitive translation vectors =   na=',na
        do 90 k=1,na
   90   write(jf,100) ((ja(i,j,k),i=1,2),j=1,3)
  100   format(1h ,'( ',3(i2,'/',i2,2x),')')
      end if
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nssum1(is,ja,jb,jc)
c
c#12  input   : is       : add(1) or subtract(-1)
c#12            ja(2,3)  : a vector
c#12            jb(2,3)  : a vector
c#12  output  : jc(2,3)  : sum of ja and jb, and take modulus
c#13  external: nsgcm3
c
c#21  to get a summation of vectors
c#21  (ja(1,1)/ja(2,1), ja(1,2)/ja(2,2), ja(1,3)/ja(2,3))  etc.
c
c#31  1989.12.28.:  n. hamada, a. yanase, and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nssum1(is,ja,jb,jc)
c
      implicit real*8(a-h,o-z)
      integer ja(2,3),jb(2,3),jc(2,3)
c
      call nsgcm3(ja(2,1),ja(2,2),ja(2,3),m)
      call nsgcm3(jb(2,1),jb(2,2),jb(2,3),n)
      id1=ja(2,1)*ja(2,2)*ja(2,3)/m*m
      id2=jb(2,1)*jb(2,2)*jb(2,3)/n*n
      call nsgcm2(id1,id2,k)
      id=id1*id2/k
      i1=ja(1,1)*(id/ja(2,1))
      i2=ja(1,2)*(id/ja(2,2))
      i3=ja(1,3)*(id/ja(2,3))
      j1=jb(1,1)*(id/jb(2,1))
      j2=jb(1,2)*(id/jb(2,2))
      j3=jb(1,3)*(id/jb(2,3))
c
      jc(1,1)=i1+j1*is
      jc(1,2)=i2+j2*is
      jc(1,3)=i3+j3*is
      jc(2,1)=id
      jc(2,2)=id
      jc(2,3)=id
c
      call nsrduc(jc(1,1),jc(2,1))
      call nsrduc(jc(1,2),jc(2,2))
      call nsrduc(jc(1,3),jc(2,3))
      jc(1,1)=mod(jc(1,1),jc(2,1))
      jc(1,2)=mod(jc(1,2),jc(2,2))
      jc(1,3)=mod(jc(1,3),jc(2,3))
c
      return
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.nstrsh(u,v)
c
c#12  noinput
c#12  output: u(3,3)
c#12          v(3,3)
c#13  noexternal
c
c#21  to get matrix u and v
c
c#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstrsh(u,v)
c
      implicit real*8(a-h,o-z)
      real*8 u(3,3),v(3,3)
c
      three=3
      do 10 j=1,3
      do 10 i=1,3
        u(i,j)=0
        v(i,j)=0
   10 continue
      u(3,3)=1.0
      v(3,3)=1.0
      a=sqrt(three)
      u(1,1)= 2/a
      u(1,2)= 1/a
      u(2,1)= 1/a
      u(2,2)= 2/a
      v(1,1)= 2/a
      v(1,2)=-1/a
      v(2,1)=-1/a
      v(2,2)= 2/a
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt0c(e,e1,e2,e3,e4,r)
c
c     (keep implicit 8-byte real)
      implicit real * 8 (a-h,o-z)
c
      real*8 e1,e2,e3,e4
      complex*16 r(4),e
c
c     (keep 16-byte complex)
      complex * 16 p1,p2,p3,p4,c1,c2,c3,c4
      complex * 16 d1,d2,d3,d4,a1,a2,a3,a4,b1,b2,b3,b4
c$$$      complex * 16 dc log
      complex * 16 cd log
c
      one=1
      d12=e1-e2
      d13=e1-e3
      d14=e1-e4
      d23=e2-e3
      d24=e2-e4
      d34=e3-e4
      d21=-d12
      d31=-d13
      d41=-d14
      d32=-d23
      d42=-d24
      d43=-d34
      d1=e-e1
      d2=e-e2
      d3=e-e3
      d4=e-e4
      a1=(d1**2)/(d21*d31*d41)
      a2=(d2**2)/(d12*d32*d42)
      a3=(d3**2)/(d13*d23*d43)
      a4=(d4**2)/(d14*d24*d34)
c     (keep dc log)
c$$$      p1= dc log(d1)
c$$$      p2= dc log(d2)
c$$$      p3= dc log(d3)
c$$$      p4= dc log(d4)
      p1= cd log(d1)
      p2= cd log(d2)
      p3= cd log(d3)
      p4= cd log(d4)
      b1=a1*d1*p1
      b2=a2*d2*p2
      b3=a3*d3*p3
      b4=a4*d4*p4
      c1=d2/d21+d3/d31+d4/d41
      c2=d1/d12+d3/d32+d4/d42
      c3=d1/d13+d2/d23+d4/d43
      c4=d1/d14+d2/d24+d3/d34
      r(1)=a1*(one-c1*p1)-b2/d21-b3/d31-b4/d41
      r(2)=a2*(one-c2*p2)-b1/d12-b3/d32-b4/d42
      r(3)=a3*(one-c3*p3)-b1/d13-b2/d23-b4/d43
      r(4)=a4*(one-c4*p4)-b1/d14-b2/d24-b3/d34
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt0i(e,e1,e2,e3,e4,dos,dosin)
c
c  ** dos and integrated dos 'coefficients' at **  c
c  ** tetrahedron corners at one energy        **  c
c  ** according to lambin and vigneron,        **  c
c  ** phys. rev. b29, 3430 (1984)              **  c
c
      implicit real*8(a-h,o-z)
      dimension dos(4),dosin(4)
      common/ntcntr/ncounter,ncounter1,ncounter2,ncounter3
      d21=e2-e1
      d31=e3-e1
      d41=e4-e1
      d32=e3-e2
      d42=e4-e2
      d43=e4-e3
      d1=e-e1
      d4=e4-e
      if(e.le.e2) then
         ncounter1 = ncounter1+1
         d2=e2-e
         d3=e3-e
         yy=d41*d31*d21
         x=d2/d21+d3/d31+d4/d41
         y=(d1*d1)/yy
         dos(1)=x*y
         dosin(1)=0.25*d1*y*(x+1.0)
         xx=d1*d1*d1
         x=xx/(d21*yy)
         dos(2)=x
         dosin(2)=0.25*d1*x
         x=xx/(d31*yy)
         dos(3)=x
         dosin(3)=0.25*d1*x
         x=xx/(d41*yy)
         dos(4)=x
         dosin(4)=0.25*d1*x
      else if(e.ge.e3) then
         ncounter3 = ncounter3+1
        d2=e-e2
        d3=e-e3
        xx=d4*d4*d4
        yy=d41*d42*d43
        x=xx/(d41*yy)
        dos(1)=x
        dosin(1)=0.25*(1.0-d4*x)
        x=xx/(d42*yy)
        dos(2)=x
        dosin(2)=0.25*(1.0-d4*x)
        x=xx/(d43*yy)
        dos(3)=x
        dosin(3)=0.25*(1.0-d4*x)
        x=d3/d43+d2/d42+d1/d41
        y=(d4*d4)/yy
        dos(4)=x*y
        dosin(4)=0.25*(1.0-d4*y*(x+1.0))
      else
         ncounter2 = ncounter2+1
         d2=e-e2
         d3=e3-e
         y=1.0/(d31*d42)+1.0/(d41*d32)
         x1=(d3*d3)/(d31*d31*d32)*(d2/d42+d1/d41)
         x2=(d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
         x3=(d3*d4*d1)/(d31*d41)*y
         dos(1)=0.5*(x1+x2+x3)
         x1=d2*d2*(d32*d2+3.0*d3*(d32+d3))/12.0
         y1=x1
         x1=x1/(d31*d31*d32*d42)
         x2=d2*(d2*d2*(d31+3.0*d21)+3.0*d3*(d2*d3+d32*(3.0*d21+d1)))
         x2=x2/12.0
         y2=x2
         x2=x2/(d31*d31*d32*d41)
         x3=d2*d2*(d42*d2+3.0*d4*(d42+d4))/12.0
         y3=x3
         x3=x3/(d41*d41*d42*d32)
         x4=d2*(d2*d2*(d41+3.0*d21)+3.0*d4*(d2*d4+d42*(3.0*d21+d1)))
         x4=x4/12.0
         y4=x4
         x4=x4/(d41*d41*d42*d31)
         x5=0.5*d2*d3*d4*(d1+d21)
         x5=x5+d2*d2*(2.0*d21*(d3+d42)+(d1+d21)*(2.0*d3+d4+d42))/12.0
         x5=x5*y/(d31*d41)
         x6=0.25*d21*d21*(d42/d41+d32/d31+1.0)/(d41*d31)
         dosin(1)=0.5*(x1+x2+x3+x4+x5)+x6
         x1=(d3*d3)/(d32*d32*d31)*(d2/d42+d1/d41)
         x2=(d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
         x3=(d3*d4*d2)/(d42*d32)*y
         dos(2)=0.5*(x1+x2+x3)
         x1=y1/(d32*d32*d31*d42)
         x2=y2/(d32*d32*d31*d41)
         x3=y3/(d42*d42*d41*d32)
         x4=y4/(d42*d42*d41*d31)
         x5=d2*d2*(d3*(d42+3.0*d4)+d32*(d42+d4))/12
         x5=x5*y/(d42*d32)
         x6=0.25*d21/d31*d21/d41
         dosin(2)=0.5*(x1+x2+x3+x4+x5)+x6
         x1=(d2*d2)/(d32*d32*d42)*(d3/d31+d4/d41)
         x2=(d1*d1)/(d31*d31*d41)*(d3/d32+d4/d42)
         x3=(d1*d2*d3)/(d32*d31)*y
         dos(3)=0.5*(x1+x2+x3)
         x1=d2*d2*d2*(3*d3+d32)/12
         y1=x1
         x1=x1/(d32*d32*d42*d31)
         x2=d2*d2*d2*(3*d4+d42)/12
         y2=x2
         x2=x2/(d32*d32*d42*d41)
         x3=d2*(d2*d31*(d2+3*d21)+3*d3*(d2*d2+3*d21*d1)+
     &        3*d21*d21*d32)/12
         y3=x3
         x3=x3/(d31*d31*d41*d32)
         x4=d2*(d2*d41*(d2+3*d21)+3*d4*(d2*d2+3*d21*d1)+3*
     &        d21*d21*d42)/12
         y4=x4
         x4=x4/(d31*d31*d41*d42)
         x5=(d2*d2)*(d3*(d21+3*d1)+d32*(d21+d1))/12
         x5=x5/(d32*d31)*y
         x6=0.25*(d21/d31*d21/d31*d21/d41)
         dosin(3)=0.5*(x1+x2+x3+x4+x5)+x6
         x1=(d2*d2)/(d42*d42*d32)*(d3/d31+d4/d41)
         x2=(d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
         x3=(d1*d2*d4)/(d41*d42)*y
         dos(4)=0.5*(x1+x2+x3)
         x1=y1/(d42*d42*d32*d31)
         x2=y2/(d42*d42*d32*d41)
         x3=y3/(d41*d41*d31*d32)
         x4=y4/(d41*d41*d31*d42)
         x5=(d2*d2)*(d4*(d21+3*d1)+d42*(d21+d1))/12
         x5=x5/(d41*d42)*y
         x6=0.25*(d21/d31*d21/d41*d21/d41)
         dosin(4)=0.5*(x1+x2+x3+x4+x5)+x6
      end if
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsttim(e1,e2,e3,e4,dos,dosin)
      implicit real * 8 (a-h,o-z)
      real*8  e1,e2,e3,e4,dos(4),dosin(4)
      real*8  tdos,e(4),d(4), esum

      e(1) = e1
      e(2) = e2
      e(3) = e3
      e(4) = e4
      esum = e1+e2+e3+e4
c$$$      tdos = 0.025d0*sum(dos(1:4))
      tdos = 0.025d0*(dos(1)+dos(2)+dos(3)+dos(4))
      do i=1,4
c$$$         dosin(i) = dosin(i) + tdos*sum(e(1:4)-e(i))
         dosin(i) = dosin(i) + tdos*(esum-4.d0*e(i))
      end do
      end subroutine nsttim
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt0r(e,e1,e2,e3,e4,r)
c
c     (keep implicit 8-byte real)
      implicit real * 8 (a-h,o-z)
c
      real*8  e,e1,e2,e3,e4,r(4)
c
      one=1
      d12=e1-e2
      d13=e1-e3
      d14=e1-e4
      d23=e2-e3
      d24=e2-e4
      d34=e3-e4
      d21=-d12
      d31=-d13
      d41=-d14
      d32=-d23
      d42=-d24
      d43=-d34
      d1=e-e1
      d2=e-e2
      d3=e-e3
      d4=e-e4
      a1=(d1**2)/(d21*d31*d41)
      a2=(d2**2)/(d12*d32*d42)
      a3=(d3**2)/(d13*d23*d43)
      a4=(d4**2)/(d14*d24*d34)
      call nsxlog(d1,p1)
      call nsxlog(d2,p2)
      call nsxlog(d3,p3)
      call nsxlog(d4,p4)
      b1=a1*d1*p1
      b2=a2*d2*p2
      b3=a3*d3*p3
      b4=a4*d4*p4
      c1=d2/d21+d3/d31+d4/d41
      c2=d1/d12+d3/d32+d4/d42
      c3=d1/d13+d2/d23+d4/d43
      c4=d1/d14+d2/d24+d3/d34
      r(1)=a1*(one-c1*p1)-b2/d21-b3/d31-b4/d41
      r(2)=a2*(one-c2*p2)-b1/d12-b3/d32-b4/d42
      r(3)=a3*(one-c3*p3)-b1/d13-b2/d23-b4/d43
      r(4)=a4*(one-c4*p4)-b1/d14-b2/d24-b3/d34
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt1i(e,nxx,nyy,nzz,ea,cdos,cind,flag)
c
c     ***  'densities of states and integrated densities of ***
c     ***  coefficients' are obtained at k-mesh points      ***
c     ***  for a single energy                              ***
c     ***  by linear interpolation in tetrahedorons         ***
c     ***
c     ***
c     ***  npx  number of mesh points in b.z.             ***
c     ***         in x-direction
c     ***  npy  number of mesh points in b.z.
c     ***         in y-direction
c     ***  npz  number of mesh points in b.z.             ***
c     ***         in z-direction
c     ***  cdos density of states 'coefficient'
c     ***  cind number  of states 'coefficient'
c     ***  flag If flag is .true., quadric correction
c     ***       will be performed
c     ***  ref.)  ph. lambin and j. p. vigneron,
c                 phys. rev. b29 (1984) 3430.
c
      implicit real*8(a-h,o-z)
c
      dimension ea(*),cdos(*),cind(*)
      dimension  ecub(2,2,2), ec(8), et(4), eb(4)
      dimension iecub(2,2,2),iec(8),iet(4),ieb(4)
      dimension dos(4),dosin(4)
      equivalence(ec(1),ecub(1,1,1))
      equivalence(iec(1),iecub(1,1,1))
      dimension iqmat(6,2)
      common/ntcntr/ncounter,ncounter1,ncounter2,ncounter3
      logical   flag
      data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
c  definition of eps  <- must be consistent with <nstts1>
      eps=dfloat(10)**(-4)
c
      npx=nxx+1
      npy=nyy+1
      npz=nzz+1
      np=npx*npy*npz
      ncub=nxx*nyy*nzz
      ntet=6*ncub
c
      do 10 ip=1,np
        cdos(ip)=0
        cind(ip)=0
   10 continue
c
      emax0=-dfloat(10)**30
      emin0= dfloat(10)**30
      do 12 ip=1,np
        if(ea(ip).gt.emax0) emax0=ea(ip)
        if(ea(ip).lt.emin0) emin0=ea(ip)
   12 continue
c                                   ============ if 1 ==
      if(e.gt.emin0) then
c
c     ***  integration over b.z. starts    ***
c
c     ***       sampling over cubes        ***
c
      icub=0
      do 20 iz=0,nzz-1
      do 20 iy=0,nyy-1
      do 20 ix=0,nxx-1
        icub=icub+1
c     ***  energies at cube corners  ***
        ni=npx*(npy*iz+iy)+ix
        emax=-dfloat(10)**30
        emin= dfloat(10)**30
        do 30 kz=1,2
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
        ecub(kx,ky,kz)=ea(ip0)
        iecub(kx,ky,kz)=ip0
        if(ea(ip0).gt.emax) emax=ea(ip0)
        if(ea(ip0).lt.emin) emin=ea(ip0)
   30   continue
c                                   ============ if 2 ==
**    if(e.ge.emax+2*eps) then
      if(e.ge.emax) then
        cind(iec(1))=cind(iec(1))+dfloat(3)/2
        cind(iec(8))=cind(iec(8))+dfloat(3)/2
        do 36 i=2,7
        cind(iec(i))=cind(iec(i))+dfloat(1)/2
   36   continue
c
      else if(e.gt.emin-2*eps) then
c     ***      six tetrahedrons      ***
c*
c     *** sampling over tetrahedrons ***
c*
        et(1)=ec(1)
        et(4)=ec(8)
        iet(1)=iec(1)
        iet(4)=iec(8)
        do 40 it=1,6
          do 42 ip=1,2
             iq=iqmat(it,ip)
             et(ip+1)=ec(iq)
             iet(ip+1)=iec(iq)
 42       continue
          do 44 m=1,4
             eb(m)=et(m)
             ieb(m)=iet(m)
 44       continue
c     ***  eb(1).le.eb(2).le.eb(3).le.eb(4)  ***
c
          call nsttod(eb,ieb)
c
          e1=eb(1)
          e2=eb(2)
          e3=eb(3)
          e4=eb(4)
c                                   ============ if 3 ==
          call nstts1(e1,e2,e3,e4)

          if(e.ge.e4) then
            do 46 i=1,4
            cind(ieb(i))=cind(ieb(i))+dfloat(1)/4
   46       continue
c
          else if(e.gt.e1) then
c
*           call nstts1(e1,e2,e3,e4)
            call nstt0i(e,e1,e2,e3,e4,dos,dosin)
            if(flag) call nsttim(e1,e2,e3,e4,dos,dosin)
c
            do 48 i=1,4
            cdos(ieb(i))=cdos(ieb(i))+dos(i)
            cind(ieb(i))=cind(ieb(i))+dosin(i)
   48       continue
          end if
c                                   ============ if 3 ==
   40   continue
      end if
c                                   ============ if 2 ==
   20 continue
c
      do 50 i=1,np
        cdos(i)=cdos(i)/ntet
        cind(i)=cind(i)/ntet
   50 continue
c
      end if
c                                   ============ if 1 ==
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt1r(e,nxx,nyy,nzz,ea,r)
c
c     ***  for a single energy                            ***
c     ***  by linear interpolation in tetrahedorons       ***
c     ***                                                 ***
c     ***  npx  number of mesh points in b.z.
c     ***         in x-direction
c     ***  npy  number of mesh points in b.z.
c     ***         in y-direction
c     ***  npz  number of mesh points in b.z.
c     ***         in z-direction
c     ***  r
c     ***
c     ***  ref.)  ph. lambin and j. p. vigneron,
c                 phys. rev. b29 (1984) 3430.
c
      implicit real*8(a-h,o-z)
c
      dimension ea(*),r(*)
      dimension  ecub(2,2,2), ec(8), et(4), eb(4)
      dimension iecub(2,2,2),iec(8),iet(4),ieb(4)
      dimension rr(4)
      equivalence(ec(1),ecub(1,1,1))
      equivalence(iec(1),iecub(1,1,1))
      dimension iqmat(6,2)
      data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
c
      npx=nxx+1
      npy=nyy+1
      npz=nzz+1
      np=npx*npy*npz
      ncub=nxx*nyy*nzz
      ntet=6*ncub
c
      do 10 ip=1,np
        r(ip)=0
   10 continue
c
c     ***  integration over b.z. starts    ***
c
c     ***       sampling over cubes        ***
c
      icub=0
      do 20 iz=0,nzz-1
      do 20 iy=0,nyy-1
      do 20 ix=0,nxx-1
        icub=icub+1
c     ***  energies at cube corners  ***
        ni=npx*(npy*iz+iy)+ix
        do 30 kz=1,2
        do 30 ky=1,2
        do 30 kx=1,2
        ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
        ecub(kx,ky,kz)=ea(ip0)
        iecub(kx,ky,kz)=ip0
   30   continue
c
c     ***      six tetrahedrons      ***
c*
c     *** sampling over tetrahedrons ***
c*
        et(1)=ec(1)
        et(4)=ec(8)
        iet(1)=iec(1)
        iet(4)=iec(8)
        do 40 it=1,6
          do 42 ip=1,2
          iq=iqmat(it,ip)
          et(ip+1)=ec(iq)
          iet(ip+1)=iec(iq)
   42     continue
          do 44 m=1,4
          eb(m)=et(m)
          ieb(m)=iet(m)
   44     continue
c     ***  eb(1).le.eb(2).le.eb(3).le.eb(4)  ***
c
          call nsttod(eb,ieb)
c
          e1=eb(1)
          e2=eb(2)
          e3=eb(3)
          e4=eb(4)
c
            call nstts1(e1,e2,e3,e4)
c           if(e1.le.0) then
c             e2=e2-e1+1.e-3
c             e3=e3-e1+1.e-3
c             e4=e4-e1+1.e-3
c             e1=1.e-3
c           end if
            call nstt0r(e,e1,e2,e3,e4,rr)
c
            do 48 i=1,4
            r(ieb(i))=r(ieb(i))+rr(i)
   48       continue
   40   continue
   20 continue
c
      do 50 i=1,np
        r(i)=r(i)/ntet
   50 continue
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt2i(idim,e,nx,ny,nz,np2,eig2,ip20,np0,ea,cd,cs,
     &                  cdos,cind)
c
      implicit  real*8(a-h,o-z)
c
      real*8  eig2(np2),cdos(np2),cind(np2)
      integer ip20(np0)
      real*8  ea(np0),cd(np0),cs(np0)
c
#ifdef __TIMER_SUB__
      call timer_sta(714)
#endif
#ifdef __TIMER_DO__
      call timer_sta(819)
#endif
      do 10 k0=1,np0
        ea(k0)=eig2(ip20(k0))
c$$$        write(6,*) ' k0,k1=',k0,ip20(k0),'   ea=',ea(k0),' <<nstt2i>>'
   10 continue
#ifdef __TIMER_DO__
      call timer_end(819)
#endif
      if(idim.eq.1) then
        call lstt1i(e,nz,ea,cd,cs)
      else if(idim.eq.2) then
        call mstt1i(e,nx,ny,ea,cd,cs)
      else if(idim.eq.3) then
        call nstt1i(e,nx,ny,nz,ea,cd,cs,.false.)
      else if(idim.eq.-3) then
        call nstt1i(e,nx,ny,nz,ea,cd,cs,.true.)
      else
        write(6,*) ' idim=',idim,' : error'
        stop '=error(idim) sub.nstt2i='
      end if
      do 20 k2=1,np2
        cdos(k2)=0
        cind(k2)=0
   20 continue
#ifdef __TIMER_DO__
      call timer_sta(820)
#endif
      do 30 k0=1,np0
        cdos(ip20(k0))=cdos(ip20(k0))+cd(k0)
        cind(ip20(k0))=cind(ip20(k0))+cs(k0)
   30 continue
#ifdef __TIMER_DO__
      call timer_end(820)
#endif
#ifdef __TIMER_SUB__
      call timer_end(714)
#endif
c
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt3i(idim,ne,e,nxx,nyy,nzz,
     &                  np2,lmnp2e,neig,eeig,ip20,np0,eawk,cdwk,cswk,
     &                  lmnp2c,lmneig,cdos,cind,eps )
c
      implicit real*8(a-h,o-z)
c
      real*8  e(0:ne)
      integer ip20(np0)
      real*8  eeig(lmnp2e,neig)
      real*8  eawk(np0), cdwk(np0), cswk(np0)
      real*8  cdos(lmnp2c,lmneig,0:ne), cind(lmnp2c,lmneig,0:ne)
      real*8  eps
      common/ntcntr/ncounter,ncounter1,ncounter2,ncounter3
c$$$      real*8 etime,wct_now,wct_start
c
c$$$      etime = 0.0d0
c$$$      call gettod(wct_start)
#ifdef __TIMER_SUB__
      call timer_sta(713)
#endif
C     eps=dfloat(10)**(-5)
c
      iloop = 0
      ncounter = 0
#ifdef __TIMER_DO__
      call timer_sta(816)
#endif
      do 10 ie=0,ne
         do 20 ieig=1,neig
           call nstt2i(idim,e(ie),nxx,nyy,nzz,np2,eeig(1,ieig),ip20,np0,
     &           eawk,cdwk,cswk,cdos(1,ieig,ie),cind(1,ieig,ie))
            iloop = iloop + (nxx*nyy*nzz*6)
 20      continue
 10   continue
#ifdef __TIMER_DO__
      call timer_end(816)
#endif

c$$$      write(6,'(" !dos iloop = ",i12)') iloop
c$$$      write(6,'(" !dos ncounter1 = ",i12)') ncounter1
c$$$      write(6,'(" !dos ncounter2 = ",i12)') ncounter2
c$$$      write(6,'(" !dos ncounter3 = ",i12)') ncounter3
c$$$      goto 31
c
c     take care of a weight on a degenerate state
c
      do 30 k2=1,np2
        ieig=1
   40   continue
        n=1
        !!$do 42 i=1,20
#ifdef __TIMER_DO__
      call timer_sta(817)
#endif
        do 42 i=1,neig
          if(ieig+i.gt.neig) go to 44
          if(dabs(eeig(k2,ieig+i)-eeig(k2,ieig)).lt.eps) then
            n=n+1
            go to 42
          end if
          go to 44
   42   continue
#ifdef __TIMER_DO__
      call timer_end(817)
#endif
        !!$write(6,100) (eeig(k2,ieig+i),i=0,20)
        !!$write(6,*) ' 21 states are degenerate! (error)'
        !!$stop 'error === in sub.nstt3i ==='

c
   44   continue
c$$$        if(n >= 2) then
c$$$           write(6,'(" !nstt3i (k2,ieig,n) = ",3i5)') k2, ieig, n
c$$$           write(6,'(" !nstt3i ",10x," eeig = ",6f8.4)')
c$$$     &          (eeig(k2,ieig+i),i=0,n-1)
c$$$        end if
#ifdef __TIMER_DO__
      call timer_sta(818)
#endif
        do 50 ie=0,ne
          c1=0
          c2=0
          do 52 i=0,n-1
            c1=c1+cdos(k2,ieig+i,ie)
            c2=c2+cind(k2,ieig+i,ie)
   52     continue
          c1=c1/n
          c2=c2/n
          do 54 i=0,n-1
            cdos(k2,ieig+i,ie)=c1
            cind(k2,ieig+i,ie)=c2
   54     continue
   50   continue
#ifdef __TIMER_DO__
      call timer_end(818)
#endif
        ieig=ieig+n
        if(ieig.lt.neig) go to 40
c
   30 continue
c$$$ 31   call gettod(wct_now)
c$$$      etime = etime + (wct_now-wct_start)*1.d-6
c$$$      write(6,'(" !nstt3i etime = ",f16.8)') etime
#ifdef __TIMER_SUB__
      call timer_end(713)
#endif
      return
c$$$  100 format(' e=',4e16.8/(3x,4e16.8))
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsttod(eb,ieb)
c
      implicit real*8(a-h,o-z)
      dimension eb(*),ieb(*)
      do 1 k=1,3
      a=eb(k)
      ih=ieb(k)
      ind=k
      do 2 i=k+1,4
      if(eb(i).lt.a) then
      a=eb(i)
      ih=ieb(i)
      ind=i
      end if
    2 continue
      eb(ind)=eb(k)
      ieb(ind)=ieb(k)
      eb(k)=a
      ieb(k)=ih
    1 continue
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstts1(e1,e2,e3,e4)
c
      implicit real*8(a-h,o-z)
c
*     eps=dfloat(10)**(-3)
      eps=dfloat(10)**(-4)
      eps2=eps/2
      eps3=(eps*3)/2
      a21=abs(e2-e1)
      a32=abs(e3-e2)
      a43=abs(e4-e3)
      if(a21.lt.eps) then
        if(a32.lt.eps) then
          if(a43.lt.eps) then
            o=(e1+e2+e3+e4)/4
            e1=o-eps3
            e2=o-eps2
            e3=o+eps2
            e4=o+eps3
          else
            e1=e3-eps*2
            e2=e3-eps
          end if
        else
          if(a43.lt.eps) then
            e1=e2-eps
            e4=e3+eps
          else
            e1=e2-eps
          end if
        end if
      else
        if(a32.lt.eps) then
          if(a43.lt.eps) then
            e3=e2+eps
            e4=e2+eps*2
          else
            o=(e2+e3)/2
            e2=o-eps2
            e3=o+eps2
          end if
        else
          if(a43.lt.eps) then
            e4=e3+eps
          end if
        end if
      end if
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nsxlog(x,p)
c
c     (keep implicit 8-byte real)
      implicit real * 8 (a-h,o-z)
c
      eps=dfloat(10)**(-4)
      y=abs(x)
c     (keep d log)
      if(y.le.eps) then
        p= d log(eps)
      else
        p= d log(y)
      end if
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine osatpz(ng1,tb1,lrb1,nat,y,kat,lmnat,iatsym)
c
      implicit real*8(a-h,o-z)
      real*8    tb1(3,48),y(3,nat),z(3),u(3)
      integer lrb1(3,3,48),kat(nat),iatsym(lmnat,48)
c
      eps=dfloat(10)**(-5)
c
      write(6,*) ' '
      write(6,*) ' symmetry operations on atoms'
      do 10 k=1,ng1
        do 20 i=1,nat
          z(1)=lrb1(1,1,k)*y(1,i)+lrb1(1,2,k)*y(2,i)+
     &         lrb1(1,3,k)*y(3,i)+tb1(1,k)
          z(2)=lrb1(2,1,k)*y(1,i)+lrb1(2,2,k)*y(2,i)+
     &         lrb1(2,3,k)*y(3,i)+tb1(2,k)
          z(3)=lrb1(3,1,k)*y(1,i)+lrb1(3,2,k)*y(2,i)+
     &         lrb1(3,3,k)*y(3,i)+tb1(3,k)
          do 30 j=1,nat
            jj=j
            u(1)=dabs(z(1)-y(1,j))
            u(2)=dabs(z(2)-y(2,j))
            u(3)=dabs(z(3)-y(3,j))
            if(u(1).lt.eps .and. u(2).lt.eps .and. u(3).lt.eps) then
              if(kat(i).eq.kat(j)) then
                go to 34
              else
                go to 32
              end if
            end if
   30     continue
          write(6,*) ' k=',k
          write(6,300) i,(y(l,i),l=1,3),(z(l),l=1,3)
          write(6,*) ' === sub.osatpz. (no atom) ==='
          stop '=osatpz(atomic position)='
c
   32     continue
          write(6,*) ' k=',k
          write(6,320) i,(y(l,i),l=1,3),jj,(z(l),l=1,3)
          write(6,*) ' === sub.osatpz. (different kind) ==='
          stop '=osatpz(atomic position)='
c
   34     continue
          iatsym(i,k)=jj
   20   continue
        write(6,100) k,(iatsym(i,k),i=1,nat)
   10 continue
      return
  300 format(' iatom=',i5,'  (',3f7.3,')  --->  no atom',
     &                 5x,'  (',3f7.3,')')
  320 format(' iatom=',i5,'  (',3f7.3,')  --->   iatom=',
     &                 i5,'  (',3f7.3,')')
  100 format(' k=',i2,2x,20i3/(7x,20i3))
      end
c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.osnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia)
c
c#12  input:     jpr     : print control
c#12             rm      : maximum of neighbor distance (standard unit)
c#12             grb(3,3): metric tensor
c#12             nat     : # of atoms in a unit cell
c#12             u(3,nat): atomic positions in a unit cell
c#12
c#12  output:    nb(nat) : # of neighbors
c#12             lmnb    : limit of # of neighbors
c#12             r (lmnb,nat): distance to the neighbors (standard unit)
c#12             ia(lmnb,nat): atom index of the neighbor
c#12  noexternal:
c
c#21  to get a distance from each atom to the neighboring atoms
c
c#31  1991.10.16.:  n. hamada
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine osnbat(jpr,rm,grb,nat,u,nb,lmnb,r,ia,w)
c
      implicit real*8(a-h,o-z)
      real*8 grb(3,3)
      real*8 u(3,nat), r(lmnb,nat), w(3,lmnb,nat)
      integer ia(lmnb,nat), nb(nat)
c
      eps=dfloat(10)**(-5)
      one=1
      pai2=8*datan(one)
      pai22=pai2**2
c
      do 10 j0=1,nat
        v1=u(1,j0)
        v2=u(2,j0)
        v3=u(3,j0)
        k=0
        do 20 j1=1,nat
          ww1=u(1,j1)
          ww2=u(2,j1)
          ww3=u(3,j1)
          w1=ww1-v1
          w2=ww2-v2
          w3=ww3-v3
            rr=sqrt(w1*grb(1,1)*w1+w1*grb(1,2)*w2+w1*grb(1,3)*w3
     &             +w2*grb(2,1)*w1+w2*grb(2,2)*w2+w2*grb(2,3)*w3
     &             +w3*grb(3,1)*w1+w3*grb(3,2)*w2+w3*grb(3,3)*w3)
c                                                 registration
            if(rr.gt.eps .and. rr.le.rm) then
              k=k+1
              if(k.gt.lmnb) then
                write(6,*) ' k=',k,'  > lmnb=',lmnb
                stop 'sub.osnbat'
              end if
              r (k,j0)=rr
              ia(k,j0)=j1
              w(1,k,j0)=w1
              w(2,k,j0)=w2
              w(3,k,j0)=w3
c             write(6,900) j0,k,ia(k,j0),r(k,j0)
            end if
   20   continue
        nb(j0)=k
c
   10 continue
c
      do 90 j0=1,nat
        call osnbst(nb(j0),r(1,j0),ia(1,j0),w(1,1,j0))
        if(jpr.ge.2) then
          write(6,*) ' atom0=',j0
          do 92  k=1,nb(j0)
          write(6,900) k,ia(k,j0),(w(i,k,j0),i=1,3),r(k,j0)
   92     continue
        end if
   90 continue
  900 format(1h ,'atom1=',i5,' (',i5,')','  w=',3f8.3,'   r=',f10.6)
      return
      end
c ==*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.osnbst(n,ra,ia,wa)
c
c#12  input :          n   : array length
c#12  in-output       ra(n): an array to be sorted
c#12                  ia(n): an array to be rearranged correspondingly
c#12                  wa(3,n): an array to be rearranged correspondingly
c#13  noexternal
c
c#21  to sort ra(n) in ascending order with an index ia(n)
c
c#31  1990.5.31.:  n. hamada (ref. 'Numerical recipes' Press et al )
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine osnbst(n,ra,ia,wa)
c
      implicit real*8(a-h,o-z)
      real*8    ra(n),wa(3,n)
      integer ia(n)
      if(n.le.1) then
        return
      end if
      l=n/2+1
      ir=n
c
   10 continue
        if(l.gt.1) then
          l=l-1
          rra=ra(l)
          iia=ia(l)
          w1a=wa(1,l)
          w2a=wa(2,l)
          w3a=wa(3,l)
        else
          rra=ra(ir)
          iia=ia(ir)
          w1a=wa(1,ir)
          w2a=wa(2,ir)
          w3a=wa(3,ir)
          ra(ir)=ra(1)
          ia(ir)=ia(1)
          wa(1,ir)=wa(1,1)
          wa(2,ir)=wa(2,1)
          wa(3,ir)=wa(3,1)
          ir=ir-1
          if(ir.eq.1) then
            ra(1)=rra
            ia(1)=iia
            wa(1,1)=w1a
            wa(2,1)=w2a
            wa(3,1)=w3a
            return
          end if
        end if
c
        i=l
        j=l+l
   20   if(j.le.ir) then
          if(j.lt.ir) then
            if(ra(j).lt.ra(j+1)) j=j+1
          end if
          if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            ia(i)=ia(j)
            wa(1,i)=wa(1,j)
            wa(2,i)=wa(2,j)
            wa(3,i)=wa(3,j)
            i=j
            j=j+j
          else
            j=ir+1
          end if
        go to 20
        end if
c
        ra(i)=rra
        ia(i)=iia
        wa(1,i)=w1a
        wa(2,i)=w2a
        wa(3,i)=w3a
c
      go to 10
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine rdprp(jpr,cname,idim,il,ngen,inv,igen,jgen,
     &                  imag,ianti,janti,
     &                  a,b,c,ca,cb,cc,nfspg,jpri_spg)
      implicit real*8(a-h,o-z)
      integer igen(ngen),jgen(2,3,ngen)
      integer janti(2,3)
      character*(60) cname
c
      rewind nfspg
      read(nfspg,*) jpr
      if(jpr.ge.jpri_spg) jpr = jpri_spg
      read(nfspg,800) cname
      read(nfspg,*) idim, il, ngen, inv

c      write(6,820) cname,idim,il,ngen,inv
c      write(6,*) 'ngen=',ngen

      if(jpri_spg.ge.1) then
          write(6,820) cname,idim,il,ngen,inv
          write(6,*) 'ngen=',ngen
      end if

      do 80 i=1,ngen
c        read (5,*)   igen(i),((jgen(j,k,i),j=1,2),k=1,3)
        read (nfspg,*)   igen(i),((jgen(j,k,i),j=1,2),k=1,3)
   80 continue
  800 format(a60)
  820 format(' == ',a60,' ==' /
     &  ' dimension=',i2,'      il=',i2,'   ngen=',i2,'   inv=',i2)
  840 format(i5,3(3x,i3,' /',i3))
      read(nfspg,*) imag

      if(imag.eq.1) then
        write(6,*) ' antiferromagnetic calculation'
      endif
      if(imag.eq.1) then
        read (nfspg,*)   ianti,((janti(j,k),j=1,2),k=1,3)
      endif
c
c --* lattice paremeters
c
c      read(5,*) a,b,c,ca,cb,cc
      read(nfspg,*) a,b,c,ca,cb,cc

c      write(6,*) '   '
c      write(6,860) a,b,c,ca,cb,cc
c  860 format(' a, b, c =',3f12.6/
c     &       'ca,cb,cc =',3f12.6)

      if(jpri_spg.ge.1) then
          write(6,*) '   '
          write(6,860) a,b,c,ca,cb,cc
 860      format(' a, b, c =',3f12.6/
     &       'ca,cb,cc =',3f12.6)
       end if
      return
      end

c
      subroutine setkp0(np2,np1,np0,lmnp0,lmnp1,lmnp2
     & ,nx,ny,nz,nxx,nyy,nzz
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     & ,nstar2,pa0,pb0,pb,ka0,ka2
     & ,nfspg,ipri_kp,ipri_spg, itrs )
      implicit real*8(a-h,o-z)
      parameter(lmnl=30)
c      character*5 schoen(48)
      character*60 cname
      integer igen(3),jgen(2,3,3)
      integer janti(2,3)
c     integer kanti(2,3)
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c      real*8  tabn(3,3),ta1n(3,48)
c      integer  lra1n(3,3,48)
c
c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      integer ka(4)
      integer ka2(4,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)
      integer ncounter, jpri_spg
      data    ncounter/0/
      save    ncounter

      ncounter = ncounter + 1
      if(jpri_spg .ge. 1 .and. ncounter .le. 1)then
         jpri_spg = 1
      else
         jpri_spg = 0
      end if
c
      eps=dfloat(10)**(-5)
c
      call rdprp(jpr,cname,idim,ill,ngen,inv,igen,jgen,
     &            imag,ianti,janti,
     &            a,b,c,ca,cb,cc,nfspg,jpri_spg)

      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,jpri_spg)

c
c
      if(idim.eq.1) then
        if(nx.ge.0) nx=0
        ny=0
      end if
      if(idim.eq.2) nz=0

      if(ipri_kp .ge. 1 .and. ncounter .le. 1) then
         jpr = 3
      else
         jpr = -1
      end if
c
        call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
        call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)
c
!        call nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &              iu21,iv21, itrs )
        call nskpb0_s(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &              iu21,iv21, itrs )
c
        if(idim.eq.1) then
        do 90 i=1,np0
        if(abs(pa0(1,i)).gt.eps .or. abs(pb0(1,i)).gt.eps .or.
     &     abs(pa0(2,i)).gt.eps .or. abs(pb0(2,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   90   continue
        end if
c
        if(idim.eq.2) then
        do 92 i=1,np0
        if(abs(pa0(3,i)).gt.eps .or. abs(pb0(3,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   92   continue
  900   format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
        end if
c
c
c
        do 10 j=1,np2
        pb(1,j)=pb0(1,ip02(j))
        pb(2,j)=pb0(2,ip02(j))
        pb(3,j)=pb0(3,ip02(j))
        nstar2(j)=0
   10   continue
        do 12 j=1,np1
   12   nstar2(ip21(j))=nstar2(ip21(j))+1
c
      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
c========================

c  100 format(8i10)
  120 format(i5,5x,3f9.3)
  200 format((3d24.16))
  220 format((3d24.16,i6))
c  300 format(5i12)
c
      return
      end

c ---------------------
c     subroutine setkp0_n is rewritten from <setkp0> by T. Yamasaki,
c                                     31th May 2003
      subroutine setkp0_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc
     & ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     & ,nx,ny,nz,nxx,nyy,nzz
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     & ,nstar2,pa0,pb0,pb,ka0,ka2
     & ,ipri_kp, itrs )

      implicit real*8(a-h,o-z)
c      integer il1
      integer ngen,inv
      integer igen(ngen),jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc

      parameter(lmnl=30)
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c      real*8  tabn(3,3),ta1n(3,48)
c      integer  lra1n(3,3,48)
c
c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      integer ka(4)
      integer ka2(4,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)


      eps=dfloat(10)**(-5)

      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,ipri_spg)


        call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
        call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)
c
!        call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &              iu21,iv21, itrs )
        call nskpb0_s(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &              iu21,iv21, itrs )
c
        if(idim.eq.1) then
        do 90 i=1,np0
        if(abs(pa0(1,i)).gt.eps .or. abs(pb0(1,i)).gt.eps .or.
     &     abs(pa0(2,i)).gt.eps .or. abs(pb0(2,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   90   continue
        end if
c
        if(idim.eq.2) then
        do 92 i=1,np0
        if(abs(pa0(3,i)).gt.eps .or. abs(pb0(3,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   92   continue
  900   format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
        end if
c
        do 10 j=1,np2
        pb(1,j)=pb0(1,ip02(j))
        pb(2,j)=pb0(2,ip02(j))
        pb(3,j)=pb0(3,ip02(j))
        nstar2(j)=0
   10   continue
        do 12 j=1,np1
   12   nstar2(ip21(j))=nstar2(ip21(j))+1

      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

  120 format(i5,5x,3f9.3)
  200 format((3d24.16))
  220 format((3d24.16,i6))
      return
      end

      subroutine setkp0_n_kt(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     &     ,nx,ny,nz,nxx,nyy,nzz
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     &     ,nstar2,pa0,pb0,pb,ka0,ka2
     &     ,ipri_kp,
     &     use_altv_rltv, altv, rltv, itrs,
     &     gen_name_in_carts )

      implicit real*8(a-h,o-z)
c      integer il1
      integer ngen,inv
      integer igen(ngen),jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc

      parameter(lmnl=30)
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c      real*8  tabn(3,3),ta1n(3,48)
c      integer  lra1n(3,3,48)
c
c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      integer ka(4)
      integer ka2(4,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)

c ----
      integer itrs
      logical use_altv_rltv, use_trs, gen_name_in_carts
      real*8 altv(3,3), rltv(3,3)
c -------
c init
      if(itrs == 1) then
         use_trs = .true.
      else
         use_trs = .false.
      end if

      eps=dfloat(10)**(-5)

      call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &     omove,jpri_spg, use_altv_rltv, altv, rltv,
     &     use_trs, gen_name_in_carts )

      call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
      call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)
c
!      call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &     iu21,iv21, itrs )
      call nskpb0_s(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &     iu21,iv21, itrs )
c
      if(idim.eq.1) then
         do 90 i=1,np0
            if(abs(pa0(1,i)).gt.eps .or. abs(pb0(1,i)).gt.eps .or.
     &           abs(pa0(2,i)).gt.eps .or. abs(pb0(2,i)).gt.eps) then
               write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
               stop '=bnkpgn(pa0,pb0)='
            end if
 90      continue
      end if
c
      if(idim.eq.2) then
         do 92 i=1,np0
            if(abs(pa0(3,i)).gt.eps .or. abs(pb0(3,i)).gt.eps) then
               write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
               stop '=bnkpgn(pa0,pb0)='
            end if
 92      continue
 900     format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
      end if
c
      do 10 j=1,np2
         pb(1,j)=pb0(1,ip02(j))
         pb(2,j)=pb0(2,ip02(j))
         pb(3,j)=pb0(3,ip02(j))
         nstar2(j)=0
 10   continue
      do 12 j=1,np1
 12      nstar2(ip21(j))=nstar2(ip21(j))+1

      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

 120  format(i5,5x,3f9.3)
 200  format((3d24.16))
 220  format((3d24.16,i6))
      return
      end

c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp, np,ka0)
c
      implicit real*8(a-h,o-z)
c
      integer ka0(4,lmnp)
c
      np=0
      do 10 iz=0,nzz
      do 10 iy=0,nyy
      do 10 ix=0,nxx
        np=np+1
        if(np.gt.lmnp) then
          write(6,*) 'ix,iy,iz=',ix,iy,iz,'   np,lmnp=',np,lmnp
          stop ' === stop in sub.ka00. (np>lmnp) ==='
        end if
        ka0(1,np)=ix *ny1*nz1
        ka0(2,np)=nx1*iy *nz1
        ka0(3,np)=nx1*ny1*iz
        ka0(4,np)=nd
c$$$        write(6,'(" !spg  (",3i5,") : ka0 (1:4",i8,") = ",4i8)')
c$$$     &       ix,iy,iz,np,ka0(1:4,np)
   10 continue

      return
      end
c
      subroutine setkp0_default(nbztyp1,altv,nx,ny,nz
     & ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     & ,nxx,nyy,nzz
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     & ,nstar2,pa0,pb0,pb,ka0,ka2
     & ,ipri_kp,ipri_spg, itrs )
c
      implicit real*8(a-h,o-z)
      parameter(lmnl=30)
c      character*5 schoen(48)
c      character*60 cname
      integer igen(3),jgen(2,3,3)
c      integer janti(2,3),kanti(2,3)
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c      real*8  tabn(3,3),ta1n(3,48)
c      integer  lra1n(3,3,48)
      real*8  altv(3,3)
c
c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      integer ka(4)
      integer ka2(4,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)
      integer ncounter
      data    ncounter/0/
      save    ncounter

      ncounter = ncounter + 1
c
      one=1
      eps=dfloat(10)**(-5)
c

      write(6,'(" << setkp0_default >>")')
      if(ipri_spg .ge. 1. and. ncounter .le. 1)
     &     write(6,*) 'nbztyp(nbztyp_spg) ', nbztyp1

      a=dsqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
      b=dsqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
      c=dsqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
      ca=(altv(1,2)*altv(1,3)+altv(2,2)*altv(2,3)+altv(3,2)*altv(3,3))
     & /(b*c)
      cb=(altv(1,3)*altv(1,1)+altv(2,3)*altv(2,1)+altv(3,3)*altv(3,1))
     & /(c*a)
      cc=(altv(1,1)*altv(1,2)+altv(2,1)*altv(2,2)+altv(3,1)*altv(3,2))
     & /(a*b)

      if(jpri_spg.ge.1) then
          write(6,*) '   '
          write(6,860) a,b,c,ca,cb,cc
  860     format(' a, b, c =',3f12.6/
     &       'ca,cb,cc =',3f12.6)
       end if

      if(nbztyp1.eq.2) then
         idim=3
         ill=1
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.3) then
         idim=3
         ill=3
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.4) then
         idim=3
         ill=2
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.5) then
         idim=3
         ill=2
         ngen=3
         inv=0
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1

cccccccc 2nd choice ccccccccccccccccccccc
         igen(2)=19
         jgen(1,1,2)=1
         jgen(2,1,2)=4
         jgen(1,2,2)=1
         jgen(2,2,2)=2
         jgen(1,3,2)=3
         jgen(2,3,2)=4
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1
cccccccc 2nd choice ccccccccccccccccccccc

      else if(nbztyp1.eq.6) then
         idim=3
         ill=0
         ngen=2
         inv=0
         igen(1)=3
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=2
         jgen(2,3,1)=3
         igen(2)=10
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1

      end if

      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,jpri_spg)
c
      if(idim.eq.1) then
        if(nx.ge.0) nx=0
        ny=0
      end if
      if(idim.eq.2) nz=0
      if(ipri_kp .ge. 1 .and. ncounter .le. 1) then
          jpr = 3
      else
          jpr = -1
      end if
c
        call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
        call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)
c
!        call nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &              iu21,iv21, itrs )
        call nskpb0_s(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &              pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &              iu21,iv21, itrs )
c
        if(idim.eq.1) then
        do 90 i=1,np0
        if(abs(pa0(1,i)).gt.eps .or. abs(pb0(1,i)).gt.eps .or.
     &     abs(pa0(2,i)).gt.eps .or. abs(pb0(2,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   90   continue
        end if
c
        if(idim.eq.2) then
        do 92 i=1,np0
        if(abs(pa0(3,i)).gt.eps .or. abs(pb0(3,i)).gt.eps) then
          write(6,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
          stop '=bnkpgn(pa0,pb0)='
        end if
   92   continue
  900   format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
        end if
c
c
c
        do 10 j=1,np2
        pb(1,j)=pb0(1,ip02(j))
        pb(2,j)=pb0(2,ip02(j))
        pb(3,j)=pb0(3,ip02(j))
        nstar2(j)=0
   10   continue
        do 12 j=1,np1
   12   nstar2(ip21(j))=nstar2(ip21(j))+1
c
      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

c========================

c  100 format(8i10)
  120 format(i5,5x,3f9.3)
  200 format((3d24.16))
  220 format((3d24.16,i6))
c  300 format(5i12)
c
      return
      end

c --- subroutine setkp0_default_n ---
c     This subroutine is revised from <setkp0_default> by T. Yamasaki,
c                                                      31th May 2003
      subroutine setkp0_default_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc
     &     ,nx,ny,nz
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     &     ,nxx,nyy,nzz
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     &     ,nstar2,pa0,pb0,pb,ka0
     &     ,ipri_kp, itrs )

      implicit real*8(a-h,o-z)
      parameter(lmnl=30)

      integer ill,ngen,inv
      integer igen(ngen), jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc
      integer nx,ny,nz,np2,np1,np0,lmnp0,lmnp1,lmnp2
      integer nxx,nyy,nzz
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)

c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)

      if(ipri_kp >= 1) then
         write(6,'(" << setkp0_default_n >>")')

         write(6,'(" !! ill = ",2i6)') ill
         write(6,'(" !! ngen = ",i6)') ngen
      end if

      if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
      if(ipri_kp >= 1) then
         do j = 1, ngen
            write(6,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j)
     &           ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j)
     &           ,jgen(2,3,j)
         enddo
      end if
      idim = 3
      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,jpri_spg)
      if(ipri_kp >= 1) write(6,'(" !! ill, il = ",2i6)') ill,il

      call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
      call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)

!      call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &     iu21,iv21, itrs )
      call nskpb0_s(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &     iu21,iv21, itrs )

      do 10 j=1,np2
        pb(1,j)=pb0(1,ip02(j))
        pb(2,j)=pb0(2,ip02(j))
        pb(3,j)=pb0(3,ip02(j))
        nstar2(j)=0
 10   continue
      do 12 j=1,np1
 12      nstar2(ip21(j))=nstar2(ip21(j))+1

      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

      if(ipri_kp >= 1)
     &     write(6,'(" !! ngen = ",i6," <<- setkp0_default_n")') ngen
      return
      end

c ======================================= added by K. Tagami ================ 12.0A
      subroutine setkp0_default_n_kt2(ill,ngen,inv,igen,jgen,
     &                               a,b,c,ca,cb,cc
     &     ,nx,ny,nz
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     &     ,nxx,nyy,nzz
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     &     ,nstar2,pa0,pb0, pa, ka0
     &     ,ipri_kp,
     &      use_altv_rltv, altv, rltv, itrs,
     &      gen_name_in_carts )

      implicit real*8(a-h,o-z)
      parameter(lmnl=30)

      integer ill,ngen,inv
      integer igen(ngen), jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc
      integer nx,ny,nz,np2,np1,np0,lmnp0,lmnp1,lmnp2
      integer nxx,nyy,nzz
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)

c      real*8  pa(3)
      real*8  pa (3,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)
c ----
      integer itrs
      logical use_altv_rltv, use_trs, gen_name_in_carts
      real*8 altv(3,3), rltv(3,3)
c -------
c init
      if(itrs == 1) then
         use_trs = .true.
      else
         use_trs = .false.
      end if

      if(ipri_kp >= 1) then
         write(6,'(" << setkp0_default_n >>")')

         write(6,'(" !! ill = ",2i6)') ill
         write(6,'(" !! ngen = ",i6)') ngen
      end if

      if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
      if(ipri_kp >= 1) then
         do j = 1, ngen
            write(6,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j)
     &           ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j)
     &           ,jgen(2,3,j)
         enddo
      end if
      idim = 3
      call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &              omove,jpri_spg, use_altv_rltv, altv, rltv,
     &              use_trs, gen_name_in_carts )

      if(ipri_kp >= 1) write(6,'(" !! ill, il = ",2i6)') ill,il

      call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
      call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)

c ----------------------------------------------
c      if ( kmode .eq. 1 ) then
c         call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,
c     &               lmnp1,lmnp2,
c     &        pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
c     &        iu21,iv21 )
c      else if ( kmode .eq. 2 ) then
c         call nskpa0_kt(ipri_kp,tab,ng1,lsa1,iv1,np0,pa0,lmnp0,
c     &                  lmnp1,lmnp2,
c     &        pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
c     &        iu21,iv21 )
c      endif
c  --
         call nskpa0_kt(ipri_kp,tab,ng1,lsa1,iv1,np0,pa0,lmnp0,
     &     lmnp1,lmnp2,
     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &     iu21,iv21, itrs )
c ------------------------------------------------
      do j=1,np2
        pa(1,j) = pa0(1,ip02(j))
        pa(2,j) = pa0(2,ip02(j))
        pa(3,j) = pa0(3,ip02(j))
        nstar2(j)=0
      end do

      do j=1,np1
         nstar2(ip21(j))=nstar2(ip21(j))+1
      end do
      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

      if(ipri_kp >= 1)
     &     write(6,'(" !! ngen = ",i6," <<- setkp0_default_n")') ngen
      return
      end

      subroutine setkp0_default_n_kt(ill,ngen,inv,igen,jgen,
     &                               a,b,c,ca,cb,cc
     &     ,nx,ny,nz
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2
     &     ,nxx,nyy,nzz
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21
     &     ,nstar2,pa0,pb0, pb, ka0
     &     ,ipri_kp,
     &      use_altv_rltv, altv, rltv, itrs,
     &      gen_name_in_carts )

      implicit real*8(a-h,o-z)
      parameter(lmnl=30)

      integer ill,ngen,inv
      integer igen(ngen), jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc
      integer nx,ny,nz,np2,np1,np0,lmnp0,lmnp1,lmnp2
      integer nxx,nyy,nzz
      real*8  omove(3)
      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)

c      real*8  pa(3)
      real*8  pb (3,lmnp2)
c      real*8  pn(3,0:lmnl)
c      integer kn(4,0:lmnl)
      real*8  pa0(3,lmnp0),pb0(3,lmnp0)
      integer ip10(lmnp0),ip20(lmnp0)
      integer ip01(lmnp1),ip21(lmnp1)
      integer ip02(lmnp2),ip12(lmnp2)
      integer iu21(lmnp1),iv21(lmnp1)
      integer nstar2(lmnp2)
      integer ka0(4,lmnp0)
c ----
      integer itrs
      logical use_altv_rltv, use_trs, gen_name_in_carts
      real*8 altv(3,3), rltv(3,3)
c -------
c init
      if(itrs == 1) then
         use_trs = .true.
      else
         use_trs = .false.
      end if

      if(ipri_kp >= 1) then
         write(6,'(" << setkp0_default_n >>")')

         write(6,'(" !! ill = ",2i6)') ill
         write(6,'(" !! ngen = ",i6)') ngen
      end if

      if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
      if(ipri_kp >= 1) then
         do j = 1, ngen
            write(6,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j)
     &           ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j)
     &           ,jgen(2,3,j)
         enddo
      end if
      idim = 3
      call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &              omove,jpri_spg, use_altv_rltv, altv, rltv,
     &              use_trs, gen_name_in_carts )

      if(ipri_kp >= 1) write(6,'(" !! ill, il = ",2i6)') ill,il

      call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
      call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,pa0)

!      call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &     iu21,iv21, itrs )
      call nskpb0_s(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
     &     iu21,iv21, itrs )

      do 10 j=1,np2
        pb(1,j)=pb0(1,ip02(j))
        pb(2,j)=pb0(2,ip02(j))
        pb(3,j)=pb0(3,ip02(j))
        nstar2(j)=0
 10   continue

      do j=1,np1
         nstar2(ip21(j))=nstar2(ip21(j))+1
      end do
      call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)

      if(ipri_kp >= 1)
     &     write(6,'(" !! ngen = ",i6," <<- setkp0_default_n")') ngen
      return
      end
c =========================================================================== 12.0A

      subroutine setspg(tabn,ng1n,ta1n,lra1n,nfspg,ipri_spg)
c
      implicit real*8(a-h,o-z)
      character*60 cname
      integer igen(3),jgen(2,3,3)
      integer janti(2,3)
c     integer kanti(2,3)
      real*8  omove(3)

      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02

      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
      real*8  tabn(3,3),ta1n(3,48)
      integer  lra1n(3,3,48)
      integer  ncounter, jpri_spg
      data     ncounter/0/
      save     ncounter

      ncounter = ncounter + 1
      if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
         jpri_spg = 1
      else
         jpri_spg = 0
      end if
c
      one=1
      eps=dble(10)**(-5)

      call rdprp(jpr,cname,idim,ill,ngen,inv,igen,jgen,
     &            imag,ianti,janti,
     &            a,b,c,ca,cb,cc,nfspg,jpri_spg)

      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,jpri_spg)
c

      do ii = 1,3
         tabn(ii,1) = tab(ii,1)
         tabn(ii,2) = tab(ii,2)
         tabn(ii,3) = tab(ii,3)
      end do
      ng1n = ng1
      do iop = 1,ng1n
         ta1n(1,iop) = ta1(1,iop)
         ta1n(2,iop) = ta1(2,iop)
         ta1n(3,iop) = ta1(3,iop)
      end do
      do iop = 1,ng1n
         do jj =1,3
            do ii = 1,3
                lra1n(ii,jj,iop) = lra1(ii,jj,iop)
            end do
         end do
      end do
c
      return
      end

      subroutine setspg_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc
     &     ,tabn,ng1n,ta1n,lra1n,ig01n)
c        Modified by T. Yamasaki(FUJITSU Lab.), 31 May 2003
c
      implicit real*8(a-h,o-z)
      integer ill,ngen,inv
      integer igen(ngen),jgen(2,3,ngen)
      real*8 a,b,c,ca,cb,cc
      real*8  omove(3)

      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02

      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c$$$      real*8  tabn(3,3),ta1n(3,48)
c$$$      integer  lra1n(3,3,48),ig01n(48)
      real*8  tabn(3,3),ta1n(3,*)
      integer  lra1n(3,3,*),ig01n(48)

      idim = 3
      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,jpri_spg)

      do ii = 1,3
         tabn(ii,1) = tab(ii,1)
         tabn(ii,2) = tab(ii,2)
         tabn(ii,3) = tab(ii,3)
      end do
      ng1n = ng1
      do iop = 1,ng1n
         ta1n(1,iop) = ta1(1,iop)
         ta1n(2,iop) = ta1(2,iop)
         ta1n(3,iop) = ta1(3,iop)
      end do
      do iop = 1,ng1n
         do jj =1,3
            do ii = 1,3
                lra1n(ii,jj,iop) = lra1(ii,jj,iop)
            end do
         end do
      end do
      do iop = 1,ng1n
         ig01n(iop) = ig01(iop)
      end do
c
      return
      end

      subroutine setspg_default(nbztyp1,altv,tabn,ng1n,ta1n,lra1n
     & ,ipri_spg)
c
      implicit real*8(a-h,o-z)
      integer igen(3),jgen(2,3,3)
c      integer janti(2,3),kanti(2,3)
      real*8  omove(3)

      integer irotr2, irotk2
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02

      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
      real*8  tabn(3,3),ta1n(3,48)
      integer  lra1n(3,3,48)
      real*8  altv(3,3)
      integer  ncounter
      data     ncounter /0/
      save     ncounter

      ncounter = ncounter + 1
c
      one=1
      eps=dble(10)**(-5)

      if(ipri_spg .ge. 1 .and. ncounter .le. 1)
     &     write(6,*) 'nbztyp ', nbztyp1
      a=dsqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
      b=dsqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
      c=dsqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
      ca=(altv(1,2)*altv(1,3)+altv(2,2)*altv(2,3)+altv(3,2)*altv(3,3))
     & /(b*c)
      cb=(altv(1,3)*altv(1,1)+altv(2,3)*altv(2,1)+altv(3,3)*altv(3,1))
     & /(c*a)
      cc=(altv(1,1)*altv(1,2)+altv(2,1)*altv(2,2)+altv(3,1)*altv(3,2))
     & /(a*b)

      if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
          write(6,*) '   '
          write(6,860) a,b,c,ca,cb,cc
  860     format(' a, b, c =',3f12.6/
     &       'ca,cb,cc =',3f12.6)
       end if

      if(nbztyp1.eq.2) then
         idim=3
         ill=1
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.3) then
         idim=3
         ill=3
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.4) then
         idim=3
         ill=2
         ngen=3
         inv=1
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1
         igen(2)=19
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1

      else if(nbztyp1.eq.5) then
         idim=3
         ill=2
         ngen=3
         inv=0
         igen(1)=5
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=0
         jgen(2,3,1)=1

cccccccc 2nd choice ccccccccccccccccccccc
         igen(2)=19
         jgen(1,1,2)=1
         jgen(2,1,2)=4
         jgen(1,2,2)=1
         jgen(2,2,2)=2
         jgen(1,3,2)=3
         jgen(2,3,2)=4
         igen(3)=25
         jgen(1,1,3)=0
         jgen(2,1,3)=1
         jgen(1,2,3)=0
         jgen(2,2,3)=1
         jgen(1,3,3)=0
         jgen(2,3,3)=1
cccccccc 2nd choice ccccccccccccccccccccc

      else if(nbztyp1.eq.6) then
         idim=3
         ill=0
         ngen=2
         inv=0
         igen(1)=3
         jgen(1,1,1)=0
         jgen(2,1,1)=1
         jgen(1,2,1)=0
         jgen(2,2,1)=1
         jgen(1,3,1)=2
         jgen(2,3,1)=3
         igen(2)=10
         jgen(1,1,2)=0
         jgen(2,1,2)=1
         jgen(1,2,2)=0
         jgen(2,2,2)=1
         jgen(1,3,2)=0
         jgen(2,3,2)=1

      end if

c
      call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &           omove,ipri_spg)
c

      do ii = 1,3
         tabn(ii,1) = tab(ii,1)
         tabn(ii,2) = tab(ii,2)
         tabn(ii,3) = tab(ii,3)
      end do
      ng1n = ng1
      do iop = 1,ng1n
         ta1n(1,iop) = ta1(1,iop)
         ta1n(2,iop) = ta1(2,iop)
         ta1n(3,iop) = ta1(3,iop)
      end do
      do iop = 1,ng1n
         do jj =1,3
            do ii = 1,3
                lra1n(ii,jj,iop) = lra1(ii,jj,iop)
            end do
         end do
      end do
c
      return
      end

      subroutine setspg_default_n(ill,ngen,inv,igen,jgen,imag,iaf,jaf
     &     ,a,b,c,ca,cb,cc
     &     ,tabn,ng1n,ta1n,lra1n,ipri_spg,ig01n)
c !                    Modified by T. Yamasaki(FUJITSU Lab.) 31 May 2003.
c
      implicit real*8(a-h,o-z)
      integer ill,ngen,inv
      integer igen(ngen),jgen(2,3,ngen)
      integer imag,iaf,jaf(2,3)
      real*8 a,b,c,ca,cb,cc
c      integer janti(2,3),kanti(2,3)
      real*8  omove(3)

      integer irotr2, irotk2
      common/nspg0 / euler(3,24), rot(3,3,48),
     &              ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &               im0(48,48),iv0(48),ng00
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02

      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c$$$      real*8  tabn(3,3),ta1n(3,48)
c$$$      integer  lra1n(3,3,48),ig01n(48)
      real*8  tabn(3,3),ta1n(3,49)
      integer  lra1n(3,3,49),ig01n(48)

      if(ipri_spg .ge. 1 ) then
          write(6,*) '   '
          write(6,860) a,b,c,ca,cb,cc
  860     format(' a, b, c =',3f12.6/
     &       'ca,cb,cc =',3f12.6)
       end if

       idim = 3
       call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc,
     &      omove,ipri_spg)

      do ii = 1,3
         tabn(ii,1) = tab(ii,1)
         tabn(ii,2) = tab(ii,2)
         tabn(ii,3) = tab(ii,3)
      end do
      ng1n = ng1
      do iop = 1,ng1n
         ta1n(1,iop) = ta1(1,iop)
         ta1n(2,iop) = ta1(2,iop)
         ta1n(3,iop) = ta1(3,iop)
      end do
      do iop = 1,ng1n
         do jj =1,3
            do ii = 1,3
                lra1n(ii,jj,iop) = lra1(ii,jj,iop)
            end do
         end do
      end do
      do iop = 1,ng1n
         ig01n(iop) = ig01(iop)
      end do

      ! Antiferro case
      if(imag==1) then
         lra1n(:,:,ng1n+1) = irot(:,:,iaf)
         do jj=1,3
            ta1n(jj,ng1n+1) = dble(jaf(1,jj))/dble(jaf(2,jj))
         end do
      end if
c
      return
      end

! =====================================- added by K. Tagami ================= 12.0A
      subroutine setspg_default_n_kt( ill, ngen, inv, igen, jgen,
     &                                imag, iaf, jaf,
     &                                a, b, c, ca, cb, cc,
     &                                tabn, ng1n, ta1n, lra1n,
     &                                ipri_spg, ig01n,
     &                                use_altv_rltv, altv, rltv,
     &                                gen_name_in_carts )
c
      implicit real*8(a-h,o-z)
      integer ill,ngen,inv
      integer igen(ngen),jgen(2,3,ngen)
      integer imag,iaf,jaf(2,3)
      real*8 a,b,c,ca,cb,cc
c      integer janti(2,3),kanti(2,3)
      real*8  omove(3)

c ---
      logical use_altv_rltv, gen_name_in_carts
      real*8 altv(3,3), rltv(3,3)
c ---

      integer irotr2, irotk2
      common/nspg0 / euler(3,24), rot(3,3,48),
     &              ieuler(3,24),irot(3,3,48),ir1234(3,48),
     &               im0(48,48),iv0(48),ng00
      common/nspg2 / eeule2(3),eeulv2(3),ieule2(3),ieulv2(3),
     &               irotr2(3,3),irotk2(3,3),ig02,iv02

      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
c      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
c      integer lrb2(3,3),lrbv(3,3)
c      real*8  rrb2(3,3),rrbv(3,3)
c      real*8  rtc2(3,3),rtcv(3,3)
c$$$      real*8  tabn(3,3),ta1n(3,48)
c$$$      integer  lra1n(3,3,48),ig01n(48)
      real*8  tabn(3,3),ta1n(3,49)
      integer  lra1n(3,3,49),ig01n(48)

      if(ipri_spg .ge. 1 ) then
          write(6,*) '   '
          write(6,860) a,b,c,ca,cb,cc
  860     format(' a, b, c =',3f12.6/
     &       'ca,cb,cc =',3f12.6)
       end if

       idim = 3
c ----------
       call tbspg_kt( idim, ill, ngen, inv, igen, jgen,
     &                a, b, c, ca, cb, cc, omove, ipri_spg,
     &                use_altv_rltv, altv, rltv, .false.,
     &                gen_name_in_carts )

c ----------

      do ii = 1,3
         tabn(ii,1) = tab(ii,1)
         tabn(ii,2) = tab(ii,2)
         tabn(ii,3) = tab(ii,3)
      end do
      ng1n = ng1
      do iop = 1,ng1n
         ta1n(1,iop) = ta1(1,iop)
         ta1n(2,iop) = ta1(2,iop)
         ta1n(3,iop) = ta1(3,iop)
      end do
      do iop = 1,ng1n
         do jj =1,3
            do ii = 1,3
                lra1n(ii,jj,iop) = lra1(ii,jj,iop)
            end do
         end do
      end do
      do iop = 1,ng1n
         ig01n(iop) = ig01(iop)
      end do

      ! Antiferro case
      if(imag==1) then
         lra1n(:,:,ng1n+1) = irot(:,:,iaf)
         do jj=1,3
            ta1n(jj,ng1n+1) = dble(jaf(1,jj))/dble(jaf(2,jj))
         end do
      end if
c
      return
      end
! ====================================================================== 12.0A

      subroutine getspgtab(tabn)
      implicit real*8(a-h,o-z)
      real*8 tabn(3,3)
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      integer i,j

      do j=1,3
         do i=1,3
            tabn(i,j) = tab(i,j)
         end do
      end do

      return
      end

c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine tbspg(idim,ill,ngen,inv,igen,jgen,
     &                 a,b,c,ca,cb,cc, omove,ipri_spg)
c
      implicit real*8(a-h,o-z)
      integer igen(ngen),jgen(2,3,ngen)
      character*5  schoen(48)
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
      real*8 ra1(3,3,48),sa1(3,3,48),rb1(3,3,48),sb1(3,3,48)
      real*8 omove(3)
      integer ncounter
      data    ncounter/0/
      save    ncounter

c ==================================== added by K. Tagami ========== 12.0A
      logical gen_name_in_carts
c
      real*8 dummy1(3,3), dummy2(3,3)
c
      dummy1 = 0.0d0
      dummy2 = 0.0d0
c
      gen_name_in_carts = .false.
c ================================================================== 12.0A
c
      ncounter = ncounter + 1

      one=1
      il=ill
      if(idim.eq.2) then
      if(il.ne.0 .and. il.ne.1 .and. il.ne.4) then
        write(6,*) ' il=',il,' : stop in sub.tbspg.'
        stop '=tbspg (il, 2D)='
      end if
      if(il.ge.1) then
        do 10 i=1,ngen
        if((igen(i).ge. 5 .and. igen(i).le.12) .or.
     &     (igen(i).ge.15 .and. igen(i).le.20) .or.
     &     (igen(i).ge.22 .and. igen(i).le.23) .or.
     &     (igen(i).ge.29 .and. igen(i).le.36) .or.
     &     (igen(i).ge.39 .and. igen(i).le.44) .or.
     &     (igen(i).ge.46 .and. igen(i).le.47)     ) then
          write(6,*) ' igen(',i,')=',igen(i),'  in sub.tbspg'
          stop '=tbspg (igen, 2D)='
        end if
   10   continue
      end if
c
      do 12 i=1,ngen
        if(jgen(1,3,i).ne.0) then
          write(6,*) ' jgen(1,3,',i,')=',jgen(1,3,i),'  in sub.tbspg'
          stop '=tbspg (jgen, 2D)='
        end if
   12 continue
      end if
c
      jf=6
       if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
          jpr = 0
       else
          jpr = -1
       end if

c ======================================== modified by K. Tagami ======== 12.0A
c      call nspace(jf,jpr,il,ngen,inv,igen,jgen,
c     &            ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1,
c     &            omove,euler1,inver1)
      call nspace(jf,jpr,il,ngen,inv,igen,jgen,
     &            ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1,
     &            omove,euler1,inver1,
     &            .false., tac, tca, tab, tba,
     &             gen_name_in_carts )
c ======================================================================= 12.0A

      if(ipri_spg .ge. 1 .and. ncounter .le. 1)
     &     write(6,120) (omove(i),i=1,3)
  120 format(/'omove=(',3f9.6,' )')
c
      call nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
      if(ipri_spg .ge. 2 .and. ncounter  .le. 1)
     & write(6,140) a,b,c,ca,cb,cc
  140 format(/'  a, b, c=',3f12.6/' ca,cb,cc=',3f12.6)

      if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
         jpr3 = 0
      else
         jpr3 = -1
      end if

c ======================================= modified by K. Tagami ========== 12.0A
c      call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
c     &                                   grc,gkc,gra,gka,grb,gkb)
      call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
     &                                   grc,gkc,gra,gka,grb,gkb,
     &                                   .false., dummy1, dummy2 )
c ======================================================================== 12.0A
      call nsgrpb(ng1,tab,tba,ta1,ra1,sa1,tb1,rb1,sb1)
c
      eps=dble(10)**(-5)
      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(ra1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' ra1(',i,j,k,') =',ra1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ='
        end if
        lra1(i,j,k)=nint(ra1(i,j,k))
      enddo
      enddo
      enddo
c
      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(sa1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' sa1(',i,j,k,') =',sa1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ='
        end if
        lsa1(i,j,k)=nint(sa1(i,j,k))
      enddo
      enddo
      enddo
c
      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(rb1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' rb1(',i,j,k,') =',rb1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ='
        end if
        lrb1(i,j,k)=nint(rb1(i,j,k))
      enddo
      enddo
      enddo

      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(sb1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' sb1(',i,j,k,') =',sb1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ='
        end if
        lsb1(i,j,k)=nint(sb1(i,j,k))
      enddo
      enddo
      enddo

      return
      end

c ============================================ added by K. Tagami ========= 12.0
      subroutine tbspg_kt(idim,ill,ngen,inv,igen,jgen,
     &                 a,b,c,ca,cb,cc, omove,ipri_spg,
     &                 use_altv_rltv, altv, rltv,
     &                 use_trs, gen_name_in_carts )
c
      implicit real*8(a-h,o-z)
      integer igen(ngen),jgen(2,3,ngen)
      character*5  schoen(48)
      common/nspg03/ tca(3,3),tac(3,3),tab(3,3),tba(3,3),
     &               tcb(3,3),tbc(3,3)
      common/nspg04/ grc(3,3),gkc(3,3),gra(3,3),gka(3,3),
     &               grb(3,3),gkb(3,3)
      common/nspg06/ id0,il,ng0,ng1,ig01(48),im1(48,48),iv1(48),
     &               ta1(3,48),lra1(3,3,48),lsa1(3,3,48),
     &               tb1(3,48),lrb1(3,3,48),lsb1(3,3,48)
      common/nspg07/ euler1(3,48),inver1(48)
      real*8 ra1(3,3,48),sa1(3,3,48),rb1(3,3,48),sb1(3,3,48)
      real*8 omove(3)
      integer ncounter
      data    ncounter/0/
      save    ncounter
c -------------------
      logical use_altv_rltv, use_trs, gen_name_in_carts
      real*8 rltv(3,3), altv(3,3)
c local
      real*8 rc1(3,3,48), sc1(3,3,48)
c -------------------
c
      ncounter = ncounter + 1

      one=1
      il=ill
      if(idim.eq.2) then
      if(il.ne.0 .and. il.ne.1 .and. il.ne.4) then
        write(6,*) ' il=',il,' : stop in sub.tbspg.'
        stop '=tbspg (il, 2D)='
      end if
      if(il.ge.1) then
        do 10 i=1,ngen
        if((igen(i).ge. 5 .and. igen(i).le.12) .or.
     &     (igen(i).ge.15 .and. igen(i).le.20) .or.
     &     (igen(i).ge.22 .and. igen(i).le.23) .or.
     &     (igen(i).ge.29 .and. igen(i).le.36) .or.
     &     (igen(i).ge.39 .and. igen(i).le.44) .or.
     &     (igen(i).ge.46 .and. igen(i).le.47)     ) then
          write(6,*) ' igen(',i,')=',igen(i),'  in sub.tbspg'
          stop '=tbspg (igen, 2D)='
        end if
 10        continue
      end if
c
      do 12 i=1,ngen
        if(jgen(1,3,i).ne.0) then
          write(6,*) ' jgen(1,3,',i,')=',jgen(1,3,i),'  in sub.tbspg'
          stop '=tbspg (jgen, 2D)='
        end if
 12      continue
      end if
c
      jf=6
       if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
          jpr = 0
       else
          jpr = -1
       end if

      call nspace(jf,jpr,il,ngen,inv,igen,jgen,
     &             ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1,
     &             omove,euler1,inver1,
     &             use_trs, tac, tca, tab, tba,
     &             gen_name_in_carts )

      if(ipri_spg .ge. 1 .and. ncounter .le. 1)
     &     write(6,120) (omove(i),i=1,3)
 120   format(/'omove=(',3f9.6,' )')
c
! --
      if ( .not. use_altv_rltv ) then
         call nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
      endif

      if(ipri_spg .ge. 2 .and. ncounter  .le. 1)
     & write(6,140) a,b,c,ca,cb,cc
 140   format(/'  a, b, c=',3f12.6/' ca,cb,cc=',3f12.6)

      if(ipri_spg .ge. 1 .and. ncounter .le. 1) then
         jpr3 = 0
      else
         jpr3 = -1
      end if

c ------------------------------------------
      call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
     &                                   grc,gkc,gra,gka,grb,gkb,
     &                                   use_altv_rltv, altv, rltv )

c
      if ( gen_name_in_carts ) then
         rc1 = ra1
         sc1 = sa1
         call nsgrpa( ng1, tac, tca, rc1, sc1, ra1,sa1 )
      endif

      call nsgrpb_kt(ng1,tab,tba,ta1,ra1,sa1,tb1,rb1,sb1)

      eps=dble(10)**(-5)

      if ( gen_name_in_carts ) then
         call check_if_unit_matirx( ng1, rc1, sc1, 'rc1 or sc1' )
      endif
      call check_if_unit_matirx( ng1, ra1, sa1, 'ra1 or sa1' )
      call check_if_unit_matirx( ng1, rb1, sb1, 'rb1 or sb1' )
c -------------------------------

      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(ra1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' ra1(',i,j,k,') =',ra1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ph1='
        end if
        lra1(i,j,k)=nint(ra1(i,j,k))
      enddo
      enddo
      enddo
c
      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(sa1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' sa1(',i,j,k,') =',sa1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ph2='
        end if
        lsa1(i,j,k)=nint(sa1(i,j,k))
      enddo
      enddo
      enddo
c
      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(rb1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' rb1(',i,j,k,') =',rb1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ph3='
        end if
        lrb1(i,j,k)=nint(rb1(i,j,k))
      enddo
      enddo
      enddo

      do k=1,ng1
      do i=1,3
      do j=1,3
        x=dabs(sb1(i,j,k))+eps/2
        if(dmod(x,one).gt.eps) then
          ind=1
          write(6,*) ' sb1(',i,j,k,') =',sb1(i,j,k)
          write(6,*) ' stop in sub.tbspg (noninteger matrix el.)'
          stop '=tbspg ph4='
        end if
        lsb1(i,j,k)=nint(sb1(i,j,k))
      enddo
      enddo
      enddo

      return
      end

      subroutine check_if_unit_matirx( ng1, rmat_in, smat_in, comment1 )
      implicit none
c inout
      integer ng1
      real*8 rmat_in(3,3,ng1), smat_in(3,3,ng1)
      character*10 comment1
c local
      integer i, j, k, l
      real*8 c1, c2, eps
c init
      eps = 1.0D-5
c begin
      Do k=1, ng1
         Do i=1, 3
            Do j=1, 3
               c1 = 0.0d0
               c2 = 0.0d0
               Do l=1, 3
                  c1 = c1 + rmat_in(i,l,k) *smat_in(l,j,k)
                  c2 = c2 + rmat_in(l,i,k) *smat_in(j,l,k)
               End do
               if (  ( i.eq.j .and. abs(c1-1.0)>eps ) .or.
     &               ( i.ne.j .and. abs(c1) > eps ) ) then
                  write(6,*) trim(comment1)
     &                              // ' is not properly obtained.'
                  write(6,*) 'i j c1 = ', i, j, c1
               endif
               if (  ( i.eq.j .and. abs(c2-1.0)>eps ) .or.
     &               ( i.ne.j .and. abs(c2) > eps ) ) then
                  write(6,*) trim(comment1)
     &                              // ' is not properly obtained.'
                  write(6,*) 'i j c2 = ', i, j, c2
               endif

            End do
         End do
      End do

      return
      end
c ======================================================================== 12.0A


c===*====1====*====2====*====3====*====4====*====5====*====6====*====7
c
c#11  sub.tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
c
c#12  input:      il, ng1, ir1234, ig01, iv0, im0, jg1
c#12  output: commons
c#13  noexternal:
c
c#21  to be compatible with the tspace package
c
c#31  1990.01.09.:  n. hamada, a. yanase and k. terakura
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
c
      implicit real*8(a-h,o-z)
c
      integer ir1234(3,48),im0(48,48),iv0(48),ig01(48),jg1(2,3,48)
c
c     for being compatible with tspace program package
c
      common/spg1  / it(3,48),im(48,48),iv(48)
      common/spg2  / il9,ng9,ig0(48),jv(2,3,48)
c
      do 90 j=1,48
      do 90 i=1,3
   90 it(i,j)=ir1234(i,j)
c$$$      do 92 j=1,48
c$$$      ig0(j)=ig01(j)
      do 92 j=1,48
         if(j<=ng1) then
            ig0(j)=ig01(j)
         else
            ig0(j) = j
         end if
   92 iv(j)=iv0(j)
      do 94 j=1,48
      do 94 i=1,48
   94 im(i,j)=im0(i,j)
      il9=il
      ng9=ng1
      do 96 i=1,2
      do 96 j=1,3
      do 96 k=1,ng1
   96 jv(i,j,k)=jg1(i,j,k)
c
c      open(unit= 2, file='a.spg', status='unknown')
c      write( 2,200) il,ng1
c      write( 2,200) ((it(i,j),i=1, 3),j=1,48)
c      write( 2,200) ((im(i,j),i=1,48),j=1,48)
c      write( 2,200) (iv  (j),j=1,48)
c      write( 2,200) (ig01(j),j=1,ng1)
c      write( 2,200) (((jg1(i,j,k),i=1,2),j=1,3),k=1,ng1)
c      close(unit= 2)
c
      return
c  200 format(24i3)
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine tspacb
c
      implicit real*8(a-h,o-z)
c
c     for being compatible with tspace program package
c
      common/spg1  / it(3,48),im(48,48),iv(48)
      common/spg2  / il,ng1,ig01(48),jg1(2,3,48)
c
      open(unit= 2, file='a.spg', status='old')
      read( 2,200) il,ng1
      read( 2,200) ((it(i,j),i=1, 3),j=1,48)
      read( 2,200) ((im(i,j),i=1,48),j=1,48)
      read( 2,200) (iv  (j),j=1,48)
      read( 2,200) (ig01(j),j=1,ng1)
      read( 2,200) (((jg1(i,j,k),i=1,2),j=1,3),k=1,ng1)
      close(unit= 2)
c
      return
  200 format(24i3)
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine tsrmhi(j,ma,mb,k1,k2,k3,is,ic,ipr)
c
      implicit real*8(a-h,o-z)
c
      character*13 cr(4)
      character*2  cd
      common/ipr/ipa(10,2)
      dimension ic(16),id(16)
      integer iprim(10)
      data    iprim /2,3,5,7,11,13,17,19,23,29/
c     data inda/0/
      data cr/' ','sqrt((1+c)/2)','sqrt(1-c**2)','sqrt((1-c)/2)'/
      jap=(j+ma)/2
      jam=(j-ma)/2
      jbp=(j+mb)/2
      jbm=(j-mb)/2
      mab=(ma-mb)/2
      jhi=1
      go to 20
      entry tsrmi(j,ma,mb,k1,k2,k3,is,ic,ipr)
      jap=j+ma
      jam=j-ma
      jbp=j+mb
      jbm=j-mb
      mab=ma-mb
      jhi=0
   20 continue
      iti=max0(0,mab)+1
      ita=min0(jap,jbm)+1
      nc=(jap+jbm-2*(iti-1))/2
     &  +(2*(iti-1)-mab)/2+1
c     if(inda.eq.1) go to 19
c     call zzzy54
c     inda=1
c  19 continue
      do 2 i=1,10
    2 ipa(i,1)=0
      do 1 i=1,16
    1 ic(i)=0
      do 4 itt=iti,ita
      it=itt-1
      do 3 i=1,10
    3 ipa(i,2)=0
      call zzzy51(2,jap-it)
      call zzzy51(2,jbm-it)
      call zzzy51(2,it)
      call zzzy51(2,it-mab)
      iw=1
      jw=1
      do 6 i=1,10
      if(itt.eq.iti) go to 16
      if(ipa(i,1).eq.ipa(i,2)) go to 6
      if(ipa(i,1).lt.ipa(i,2)) go to 7
      iw=iw*(iprim(i)**(ipa(i,1)-ipa(i,2)))
      go to 6
    7 jw=jw*(iprim(i)**(ipa(i,2)-ipa(i,1)))
   16 ipa(i,1)=ipa(i,2)
    6 continue
      if(jw.eq.1) go to 5
      do 8 i=1,nc
    8 ic(i)=ic(i)*jw
    5 iw=iw*(1-2*mod(it,2))
      i1=jap+jbm-2*it
      i2=2*it-mab
      j1=i1/2
      j2=i2/2
      call zzzy52(id,j1,j2)
      do 9 i=1,nc
    9 ic(i)=ic(i)+iw*id(i)
    4 continue
      i1=mod(i1,2)
      i2=mod(i2,2)
      is=i1+i2
      if(i1.eq.0.and.i2.eq.1) is=3
      iw=0
      do 11 i=1,nc
      call zzzy53(iw,ic(i),jw)
      iw=jw
      if(jw.eq.1) go to 12
   11 continue
      do 13 i=1,nc
   13 ic(i)=ic(i)/jw
      do 17 i=1,10
   18 if(mod(jw,iprim(i)).ne.0) go to 17
      ipa(i,1)=ipa(i,1)-1
      jw=jw/iprim(i)
      if(jw.eq.1) go to 12
      go to 18
   17 continue
   12 continue
      do 10 i=1,10
   10 ipa(i,1)=ipa(i,1)*(-2)
      call zzzy51(1,jap)
      call zzzy51(1,jam)
      call zzzy51(1,jbp)
      call zzzy51(1,jbm)
      if(jhi.eq.1) ipa(1,1)=ipa(1,1)-j+1
      if(jhi.eq.0) ipa(1,1)=ipa(1,1)-2*j
      k1=1
      k2=1
      k3=1
      k4=1
      do 14 i=1,10
      if(ipa(i,1).eq.0) go to 14
      if(ipa(i,1).lt.0) go to 15
      k1=k1*(iprim(i)**(ipa(i,1)/2))
      if(mod(ipa(i,1),2).eq.1) k3=k3*iprim(i)
      go to 14
   15 k2=k2*(iprim(i)**(iabs(ipa(i,1))/2))
      if(mod(-ipa(i,1),2).eq.1) k4=k4*iprim(i)
   14 continue
      if(k4.eq.1) go to 21
      k3=k3*k4
      k2=k2*k4
   21 continue
      if(ipr.ne.1) return
      cd='  '
      if(jhi.eq.1) cd='/2'
      if(jhi.eq.0.and.j.gt.8) go to 22
      if(jhi.eq.1.and.j.gt.15) go to 22
      write(6,600) j,cd,ma,cd,mb,cd,k1,k2,k3,cr(is+1)
  600 format(' j=',i3,a2,2x,'ma=',i3,a2,'  mb=',i3,a2,'  (',
     &  i3,'/',i5,')sqrt(',i5,')',a13)
      go to 23
   22 continue
      write(6,601) j,cd,ma,cd,mb,cd,k1,k2,k3,cr(is+1)
  601 format(' j=',i3,a2,2x,'ma=',i3,a2,'  mb=',i3,a2,'  (',
     &  i3,'/',i5,')sqrt(',i10,')',a13)
   23 continue
      write(6,660) (ic(i),i-1,i=1,nc)
  660 format(2x,i9,'c**',i2,i9,'c**',i2,i9,'c**',i2,
     &           i9,'c**',i2,i9,'c**',i2)
      return
      end
      subroutine wtetra(nxx,nyy,nzz,np0,np2,ip20,iwt,ip2cub,ip2cub_wk)
      use m_Timing, only : tstatc0_begin, tstatc0_end
      use m_spg_tetra
c
      implicit none
c
      real*8  ecub(2,2,2), ec(8)
      integer iecub(2,2,2),iec(8),iet(4),ieb(4)
      equivalence(ec(1),ecub(1,1,1))
      equivalence(iec(1),iecub(1,1,1))
      integer iqmat(6,2)
      data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
      integer np0,np2
      integer ip20(np0),iwt(np2)

      integer npx,npy,npz,nxx,nyy,nzz,np,ntet,ncub
      integer ip,ip0,ip1,ip2,ix,iy,iz,ni,kx,ky,kz,it,iq,m,icub,jcub
      integer ip2cub_wk(9,nxx*nyy*nzz),ip2cub(nxx*nyy*nzz)
      integer :: id_sname = -1
c
      if(wtetra_vars_ready(nxx*nyy*nzz)) then
        call get_wtetra_vars(ip2cub)
        return
      endif
      npx=nxx+1
      npy=nyy+1
      npz=nzz+1
      np=npx*npy*npz
      ncub=nxx*nyy*nzz
      ntet=6*ncub

      do 10 ip=1,np2
      iwt(ip)=0
 10   continue
      do 20 iz=0,nzz-1
      do 20 iy=0,nyy-1
      do 20 ix=0,nxx-1
      ni=npx*(npy*iz+iy)+ix
      do 30 kz=1,2
      do 30 ky=1,2
      do 30 kx=1,2
      ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
      if(ip0.gt.np0) then
         write(6,*) ' ip0, np0  ',ip0,np0
         stop ' wtetra -- ip0.ne.np0'
      endif
      iecub(kx,ky,kz)=ip0
 30   continue
      iet(1)=iec(1)
      iet(4)=iec(8)
      do 40 it=1,6
      do 42 ip=1,2
      iq=iqmat(it,ip)
      iet(ip+1)=iec(iq)
 42   continue
      do 44 m=1,4
      ieb(m)=iet(m)
      if(ip20(ieb(m)).gt.np2) then
         write(6,*) ' ip20, np2  ',ip20(ieb(m)),np2
         stop ' wtetra -- ip20.ne.np2'
      endif
      iwt(ip20(ieb(m))) = iwt(ip20(ieb(m))) + 1
 44   continue
 40   continue
 20   continue

      icub=0
      do 50 iz=0,nzz-1
      do 50 iy=0,nyy-1
      do 50 ix=0,nxx-1
      icub=icub+1
      ni=npx*(npy*iz+iy)+ix
      ip=0
      do 60 kz=1,2
      do 60 ky=1,2
      do 60 kx=1,2
      ip=ip+1
      ip0=ni+npx*(npy*(kz-1)+ky-1)+kx
      iec(ip)=ip20(ip0)
 60   continue
      do 71 ip1=2,6
      do 71 ip2=ip1+1,7
      if(iec(ip1).gt.iec(ip2)) then
         ip=iec(ip1)
         iec(ip1)=iec(ip2)
         iec(ip2)=ip
      endif
 71   continue
      if(iec(1).lt.iec(8)) then
         ip2cub_wk(1,icub)=iec(1)
         ip2cub_wk(8,icub)=iec(8)
      else
         ip2cub_wk(8,icub)=iec(1)
         ip2cub_wk(1,icub)=iec(8)
      endif
      ip0=iec(1)+iec(8)
      do 72 ip=2,7
      ip2cub_wk(ip,icub)=iec(ip)
      ip0=ip0+iec(ip)
 72   continue
      ip2cub(icub)=icub
      ip2cub_wk(9,icub)=ip0
 50   continue

      call tstatc0_begin('wtetra0 ',id_sname)
      ncub=nxx*nyy*nzz
      do 80 icub=2,ncub
      do 81 jcub=1,icub-1
      if(ip2cub_wk(9,icub).eq.ip2cub_wk(9,jcub)) then
         do 82 ip=1,8
         if(ip2cub_wk(ip,icub).ne.ip2cub_wk(ip,jcub)) goto 81
 82      continue
         ip2cub(icub)=jcub
         goto 80
      endif
 81   continue
 80   continue
      call tstatc0_end(id_sname)

      jcub=0
      do 90 icub=1,ncub
c$$$      write(6,'(2x,i5,''  -->'',i5)') icub,ip2cub(icub)
      if(icub.eq.ip2cub(icub)) jcub=jcub+1
 90   continue
c$$$      write(6,'(2x,''number of cube which should be calculated'',i5
c$$$     &     ,"<<wtetra>>")')  jcub

      call set_wtetra_vars(nxx*nyy*nzz, ip2cub)

      return
      end
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine zzzy51(k,j)
c
      implicit real*8(a-h,o-z)
c
      common/ipr/ipa(10,2)
      integer iprim(10)
      data    iprim /2,3,5,7,11,13,17,19,23,29/
      if(j.le.1) return
      do 1 i=2,j
      n=i
      do 2 ip=1,10
    3 if(mod(n,iprim(ip)).ne.0) go to 2
      ipa(ip,k)=ipa(ip,k)+1
      n=n/iprim(ip)
      if(n.eq.1) go to 1
      go to 3
    2 continue
    1 continue
      return
      end
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine zzzy52(ib,m,n)
c
      implicit real*8(a-h,o-z)
c
      dimension ic(16,16),ib(16)
c     go to 5
c     entry zzzy54
      ic(1,1)=1
      ic(1,2)=1
      ic(2,2)=1
      do 1 i=3,16
      ic(1,i)=1
      ic(i,i)=1
      do 2 j=1,i-2
    2 ic(j+1,i)=ic(j,i-1)+ic(j+1,i-1)
    1 continue
c     return
c   5 continue
      nx=m+n+1
      do 3 i=1,nx
      ki=max0(1,i-n)
      ka=min0(m+1,i)
      iw=0
      do 4 k=ki,ka
    4 iw=iw+ic(k,m+1)*ic(i-k+1,n+1)
     &     *(1-2*mod(i-k,2))
      ib(i)=iw
    3 continue
      return
      end
c
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine zzzy53(i,j,k)
c
      implicit real*8(a-h,o-z)
c
      ii=iabs(i)
      jj=iabs(j)
      if(ii.eq.0) go to 1
      if(jj.eq.0) go to 4
    3 kk=mod(ii,jj)
      if(kk.eq.0) go to 2
      ii=jj
      jj=kk
      go to 3
    2 k=jj
      return
    1 k=jj
      return
    4 k=ii
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt4i(ipri,idim,newindows,nxyz_tetra,np2,lmnp2e,neig,eeig,
     &     ip20,np0,lmnp2c,lmneig,mtetra,nttra,deltae,
     &     nep,dos,dosin,cdos,cind)
c
c     nstt0i, nstt1i, nstt2i are nstt3i are merged into
c  this subroutine nstt4i
c               by T. Yamasaki, Aug 2007
c
      implicit real*8(a-h,o-z)
      real*8 eeig(lmnp2e,lmneig)
      integer nxyz_tetra(3)
      integer ip20(np0), ieb(4),nttra(mtetra,4)
      real*8 dos(0:nep,4),dosin(0:nep,4)
      real*8 cdos(0:newindows,lmnp2c,lmneig)
     &     , cind(0:newindows,lmnp2c,lmneig)
      real*8 eb(4)
      integer ncounter,ncounter1,ncounter2,ncounter3
      real*8 etime1,etime2,etime3, wct_now, wct_start,etime4,etime5
     &     ,etime0

      eps = dfloat(10)**(-5)

      nxx = nxyz_tetra(1); nyy = nxyz_tetra(2); nzz = nxyz_tetra(3)
      npx = nxx+1
      npy = nyy+1
      npz = nzz+1
      np = npx*npy*npz
      ncub = nxx*nyy*nzz
      ntet = 6*ncub

      if(ipri>=2) then
         write(6,'(" lmnp2e, lmneig, mtetra, newindows, lmnp2c,"
     &,"nep = ",6i8)') lmnp2e,lmneig,mtetra,newindows,lmnp2e,nep
         write(6,'(" np = ",i5," np2 = ",i5)')np,np2
      end if
      es =  dfloat(10)**30
      ee = -dfloat(10)**30
      do ib = 1, neig
      do ip = 1, lmnp2e
         if(eeig(ip,ib) < es) es = eeig(ip,ib)
         if(eeig(ip,ib) > ee) ee = eeig(ip,ib)
      end do
      end do
      if(ipri>=2) write(6,'(" es,ee = ",2d12.4)') es,ee

      ncounter = 0
      ncounter1 = 0
      ncounter2 = 0
      ncounter3 = 0
      etime0 = 0.d0
      etime1 = 0.d0
      etime2 = 0.d0
      etime3 = 0.d0
      etime4 = 0.d0
      etime5 = 0.d0

      iloop = 0
      do iv = 1, ntet
c$$$         write(6,'(" iv = ",i8)') iv
c$$$         write(6,'(" -- nttra  = ",4i8)') (nttra(iv,ie),ie=1,4)
         do ib = 1, neig
            call gettod(wct_start)
            do ie = 1,4
c$$$               if(nttra(iv,ie).gt.np0 .or. nttra(iv,ie).le.0) then
c$$$                  write(6,'(" nttra(",i5,",",i5,") = ",i5)')
c$$$     &                 iv,ie,nttra(iv,ie)
c$$$                  stop 'nttra illegal'
c$$$               end if
               ieb(ie) = ip20(nttra(iv,ie))
               eb(ie) = eeig(ieb(ie),ib)
            end do

            call nsttod(eb,ieb)

            e1 = eb(1)
            e2 = eb(2)
            e3 = eb(3)
            e4 = eb(4)
            call nstts1(e1,e2,e3,e4)
            if(e1 > ee) cycle
            ns = (e1 - es)/deltae+1
            ne = (e4 - es)/deltae+1
            if(ns < 0) ns = 0
            if(ne > nep) ne = nep
c$$$            write(6,'(" iv,ib = ",2i5," e = ",4f8.4," ns,ne=",2i8)')
c$$$     &           iv,ib, e1,e2,e3,e4,ns,ne

            d21=e2-e1
            d31=e3-e1
            d41=e4-e1
            d32=e3-e2
            d42=e4-e2
            d43=e4-e3

            iloop = iloop + (ne-ns+1)
            ne1 = (e2-es)/deltae
            if(es+deltae*(ne1-0.5)>e2) ne1 = ne1-1
            ns2 = ne1+1
            ne2 = (e3-es)/deltae
            if(es+deltae*(ne2-0.5)>e3) ne2 = ne2-1
            ns3 = ne2+1

            if(ipri>=2) then
               write(6,'(" ns, ne1, ns2, ne2, ns3, ne = ",6i9)')
     &              ns ,ne1, ns2, ne2, ns3, ne
            end if

c$$$            do idos = ns, ne
            yy1 = d41*d31*d21
            d21yy1 = d21*yy1
            d31yy1 = d31*yy1
            d41yy1 = d41*yy1
            call gettod(wct_now)
            etime0 = etime0 + (wct_now-wct_start)*1.d-6

            call gettod(wct_start)
            do idos = ns, ne1
               e = es+deltae*(idos-0.5)
               d1=e-e1
               d4=e4-e
               ncounter1 = ncounter1+1
               d2=e2-e
               d3=e3-e
c$$$                  yy=d41*d31*d21
               x=d2/d21+d3/d31+d4/d41
               y=(d1*d1)/yy1
               dos(idos,1)=x*y
               dosin(idos,1)=0.25*d1*y*(x+1.0)
               xx=d1*d1*d1
               x=xx/d21yy1
               dos(idos,2)=x
               dosin(idos,2)=0.25*d1*x
               x=xx/d31yy1
               dos(idos,3)=x
               dosin(idos,3)=0.25*d1*x
               x=xx/d41yy1
               dos(idos,4)=x
               dosin(idos,4)=0.25*d1*x
            end do
            call gettod(wct_now)
            etime1 = etime1 + (wct_now-wct_start)*1.d-6

            yy2 = d41*d42*d43
            d41yy2 = d41*yy2
            d42yy2 = d42*yy2
            d43yy2 = d43*yy2

            call gettod(wct_start)
            do idos = ns3,ne
               e = es+deltae*(idos-0.5)
               d1=e-e1
               d4=e4-e
               ncounter3 = ncounter3+1
               d2=e-e2
               d3=e-e3
               xx=d4*d4*d4
c$$$                  yy=d41*d42*d43
               x=xx/d41yy2
               dos(idos,1)=x
               dosin(idos,1)=0.25*(1.0-d4*x)
               x=xx/d42yy2
               dos(idos,2)=x
               dosin(idos,2)=0.25*(1.0-d4*x)
               x=xx/d43yy2
               dos(idos,3)=x
               dosin(idos,3)=0.25*(1.0-d4*x)
               x=d3/d43+d2/d42+d1/d41
               y=(d4*d4)/yy2
               dos(idos,4)=x*y
               dosin(idos,4)=0.25*(1.0-d4*y*(x+1.0))
            end do
            call gettod(wct_now)
            etime3 = etime3 + (wct_now-wct_start)*1.d-6
            yy3 = 1.0/(d31*d42)+1.0/(d41*d32)
            d31d31d32d41 = d31*d31*d32*d41
            d31d41       = d31*d41
            d42d32       = d42*d32
            d32d32d31    = d32*d32*d31

            call gettod(wct_start)
            do idos=ns2,ne2
               e = es+deltae*(idos-0.5)
               d1=e-e1
               d4=e4-e
               ncounter2 = ncounter2+1
               d2=e-e2
               d3=e3-e
c$$$                  y=1.0/(d31*d42)+1.0/(d41*d32)
               x1=(d3*d3)/(d31*d31*d32)*(d2/d42+d1/d41)
               x2=(d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
               x3=(d3*d4*d1)/(d31d41)*yy3
               dos(idos,1)=0.5*(x1+x2+x3)
               x1=d2*d2*(d32*d2+3.0*d3*(d32+d3))/12.0
               y1=x1
               x1=x1/(d31*d31*d32*d42)
               x2=d2*(d2*d2*(d31+3.0*d21)+
     &              3.0*d3*(d2*d3+d32*(3.0*d21+d1)))
               x2=x2/12.0
               y2=x2
               x2 = x2/d31d31d32d41
               x3=d2*d2*(d42*d2+3.0*d4*(d42+d4))/12.0
               y3=x3
               x3=x3/(d41*d41*d42d32)
               x4=d2*(d2*d2*(d41+3.0*d21)+
     &              3.0*d4*(d2*d4+d42*(3.0*d21+d1)))
               x4=x4/12.0
               y4=x4
               x4=x4/(d41*d41*d42*d31)
               x5=0.5*d2*d3*d4*(d1+d21)
               x5=x5+d2*d2*(2.0*d21*(d3+d42)+
     &              (d1+d21)*(2.0*d3+d4+d42))/12.0
               x5=x5*yy3/(d31d41)
               x6=0.25*d21*d21*(d42/d41+d32/d31+1.0)/(d41*d31)
               dosin(idos,1)=0.5*(x1+x2+x3+x4+x5)+x6
               x1=(d3*d3)/(d32d32d31)*(d2/d42+d1/d41)
               x2=(d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
               x3=(d3*d4*d2)/(d42d32)*yy3
               dos(idos,2)=0.5*(x1+x2+x3)
               x1=y1/(d32d32d31*d42)
               x2=y2/(d32d32d31*d41)
               x3=y3/(d42*d42*d41*d32)
               x4=y4/(d42*d42*d41*d31)
               x5=d2*d2*(d3*(d42+3.0*d4)+d32*(d42+d4))/12
               x5=x5*yy3/(d42d32)
               x6=0.25*d21/d31*d21/d41
               dosin(idos,2)=0.5*(x1+x2+x3+x4+x5)+x6
               x1=(d2*d2)/(d32*d32*d42)*(d3/d31+d4/d41)
               x2=(d1*d1)/(d31*d31d41)*(d3/d32+d4/d42)
               x3=(d1*d2*d3)/(d32*d31)*yy3
               dos(idos,3)=0.5*(x1+x2+x3)
               x1=d2*d2*d2*(3*d3+d32)/12
               y1=x1
               x1=x1/(d32*d32*d42*d31)
               x2=d2*d2*d2*(3*d4+d42)/12
               y2=x2
               x2=x2/(d32*d32*d42*d41)
               x3=d2*(d2*d31*(d2+3*d21)+3*d3*(d2*d2+3*d21*d1)+
     &              3*d21*d21*d32)/12
               y3=x3
               x3=x3/(d31*d31d41*d32)
               x4=d2*(d2*d41*(d2+3*d21)+3*d4*(d2*d2+3*d21*d1)+3*
     &              d21*d21*d42)/12
               y4=x4
               x4=x4/(d31*d31d41*d42)
               x5=(d2*d2)*(d3*(d21+3*d1)+d32*(d21+d1))/12
               x5=x5/(d32*d31)*yy3
               x6=0.25*(d21/d31*d21/d31*d21/d41)
               dosin(idos,3)=0.5*(x1+x2+x3+x4+x5)+x6

               x1=(d2*d2)/(d42*d42d32)*(d3/d31+d4/d41)
               x2=(d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
               x3=(d1*d2*d4)/(d41*d42)*yy3
               dos(idos,4)=0.5*(x1+x2+x3)
               x1=y1/(d42*d42d32*d31)
               x2=y2/(d42*d42d32*d41)
               x3=y3/(d41*d41*d31*d32)
               x4=y4/(d41*d41*d31*d42)
               x5=(d2*d2)*(d4*(d21+3*d1)+d42*(d21+d1))/12
               x5=x5/(d41*d42)*yy3
               x6=0.25*(d21/d31*d21/d41*d21/d41)
               dosin(idos,4)=0.5*(x1+x2+x3+x4+x5)+x6
            end do
            call gettod(wct_now)
            etime2 = etime2 + (wct_now-wct_start)*1.d-6
            call gettod(wct_start)
            if(idim .eq. -3) then
               esum = e1+e2+e4+e4
               do idos = ns, ne
                  tdos = 0.025d0*(dos(idos,1)+dos(idos,2)
     &                 +dos(idos,3)+dos(idos,4))
                  dosin(idos,1)=dosin(idos,1)+tdos*(esum-4.d0*e1)
                  dosin(idos,2)=dosin(idos,2)+tdos*(esum-4.d0*e2)
                  dosin(idos,3)=dosin(idos,3)+tdos*(esum-4.d0*e1)
                  dosin(idos,4)=dosin(idos,4)+tdos*(esum-4.d0*e2)
               end do
            end if
            do ie = ns, ne
               cdos(ie,ieb(1),ib)=cdos(ie,ieb(1),ib)+dos(ie,1)
               cdos(ie,ieb(2),ib)=cdos(ie,ieb(2),ib)+dos(ie,2)
               cdos(ie,ieb(3),ib)=cdos(ie,ieb(3),ib)+dos(ie,3)
               cdos(ie,ieb(4),ib)=cdos(ie,ieb(4),ib)+dos(ie,4)
               cind(ie,ieb(1),ib)=cind(ie,ieb(1),ib)+dosin(ie,1)
               cind(ie,ieb(2),ib)=cind(ie,ieb(2),ib)+dosin(ie,2)
               cind(ie,ieb(3),ib)=cind(ie,ieb(3),ib)+dosin(ie,3)
               cind(ie,ieb(4),ib)=cind(ie,ieb(4),ib)+dosin(ie,4)
            end do
            do idos = ne+1, nEwindows
               cind(idos,ieb(1),ib) = cind(idos,ieb(1),ib)+0.25d0
               cind(idos,ieb(2),ib) = cind(idos,ieb(2),ib)+0.25d0
               cind(idos,ieb(3),ib) = cind(idos,ieb(3),ib)+0.25d0
               cind(idos,ieb(4),ib) = cind(idos,ieb(4),ib)+0.25d0
            end do
            call gettod(wct_now)
            etime4 = etime4 + (wct_now-wct_start)*1.d-6
         end do
      end do
      call gettod(wct_start)
      if(ipri>=2) write(6,'(" before cdos, cind substitution")')

      do ib = 1, neig
         do ip = 1, lmnp2c
            do idos = 0, newindows
               cdos(idos,ip,ib) = cdos(idos,ip,ib)/ntet
               cind(idos,ip,ib) = cind(idos,ip,ib)/ntet
            end do
         end do
      end do
c
c ---- following lines are the later part of subroutine nstt3i --
c     take care of a weight on a degenerate state
c
      do 30 k2=1,np2
        ieig=1
   40   continue
        n=1
        !!$do 42 i=1,20
        do 42 i=1,neig
          if(ieig+i.gt.neig) go to 44
          if(dabs(eeig(k2,ieig+i)-eeig(k2,ieig)).lt.eps) then
            n=n+1
            go to 42
          end if
          go to 44
   42   continue

   44   continue
        do 50 ie=0,ne
          c1=0
          c2=0
          do 52 i=0,n-1
            c1=c1+cdos(ie,k2,ieig+i)
            c2=c2+cind(ie,k2,ieig+i)
   52     continue
          c1=c1/n
          c2=c2/n
          do 54 i=0,n-1
            cdos(ie,k2,ieig+i)=c1
            cind(ie,k2,ieig+i)=c2
   54     continue
   50   continue
        ieig=ieig+n
        if(ieig.lt.neig) go to 40
c
   30 continue
c -------------------------------

      call gettod(wct_now)
      if(ipri>=2) then
         etime5 = etime5 + (wct_now-wct_start)*1.d-6
         write(6,'(" !dos iloop    = ",i12)') iloop
         write(6,'(" !dos ncounter1 = ",i12)') ncounter1
         write(6,'(" !dos ncounter2 = ",i12)') ncounter2
         write(6,'(" !dos ncounter3 = ",i12)') ncounter3
         write(6,'(" !dos etime0    = ",f16.8)') etime0
         write(6,'(" !dos etime1    = ",f16.8)') etime1
         write(6,'(" !dos etime2    = ",f16.8)') etime2
         write(6,'(" !dos etime3    = ",f16.8)') etime3
         write(6,'(" !dos etime4    = ",f16.8)') etime4
         write(6,'(" !dos etime5    = ",f16.8)') etime5
      end if
      return
      end
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
c
      subroutine nstt5i(ipri,idim,es,ee,newindows,nxyz_tetra,
     &     np2,lmnp2e,neig,
     &     eeig,ip20,np0,lmnp2c,lmneig,mtetra,nttra,deltae,dos_weight,
     &     nep,dos,dosin,cdos,csumdos)
c
c     nstt0i, nstt1i, nstt2i and nstt3i are merged into
c  this subroutine nstt5i
c               by T. Yamasaki, Aug 2007
c
      use m_Parallelization, only : npes,mype,MPI_CommGroup
      implicit real*8(a-h,o-z)
      include 'mpif.h'
      real*8 eeig(lmnp2e,lmneig)
      integer nxyz_tetra(3)
      integer ip20(np0), ieb(4),nttra(mtetra,4)
      real*8 dos(0:nep,4),dosin(0:nep,4),dos_weight(lmneig,lmnp2c)
      real*8 cdos(0:newindows)
     &     , csumdos(0:newindows)
!MPI
      real*8,allocatable, dimension(:,:) :: dos_mpi,dosin_mpi
      real*8,allocatable, dimension(:) :: cdos_mpi,csumdos_mpi
!MPI

      real*8 eb(4)
      integer ncounter,ncounter1,ncounter2,ncounter3
      real*8 etime1,etime2,etime3, wct_now, wct_start,etime4,etime5
     &     ,etime0
      integer ierr,itmp

      eps = dfloat(10)**(-4)

      nxx = nxyz_tetra(1); nyy = nxyz_tetra(2); nzz = nxyz_tetra(3)
      npx = nxx+1
      npy = nyy+1
      npz = nzz+1
      np = npx*npy*npz
      ncub = nxx*nyy*nzz
      ntet = 6*ncub

      dos = 0.0d0
      dosin = 0.0d0
      cdos = 0.0d0
      csumdos = 0.0d0
      if(ipri.ge.2) write(6,'(" np = ",i5," np2 = ",i5," idim = ",i5)')
     &     np,np2,idim
c$$$      es =  dfloat(10)**30
c$$$      ee = -dfloat(10)**30
c$$$      do ib = 1, neig
c$$$      do ip = 1, lmnp2e
c$$$         if(eeig(ip,ib) < es) es = eeig(ip,ib)
c$$$         if(eeig(ip,ib) > ee) ee = eeig(ip,ib)
c$$$      end do
c$$$      end do
c$$$      es = es - 0.005
c$$$      ee = ee + 0.005
c$$$      write(6,'(" es,ee = ",2f10.6)') es,ee

      ncounter = 0
      ncounter1 = 0
      ncounter2 = 0
      ncounter3 = 0
      etime0 = 0.d0
      etime1 = 0.d0
      etime2 = 0.d0
      etime3 = 0.d0
      etime4 = 0.d0
      etime5 = 0.d0

      iloop = 0
      do iv = 1, ntet
         if(mod(iv,npes)/=mype) cycle
         if(ipri.ge.3) then
            write(6,'(" iv = ",i8)') iv
            write(6,'(" -- nttra  = ",4i8)') (nttra(iv,ie),ie=1,4)
         endif
         do ib = 1, neig
            call gettod(wct_start)
            do ie = 1,4
               ieb(ie) = ip20(nttra(iv,ie))
               eb(ie) = eeig(ieb(ie),ib)
            end do

            call nsttod(eb,ieb)

            e1 = eb(1)
            e2 = eb(2)
            e3 = eb(3)
            e4 = eb(4)
            call nstts1(e1,e2,e3,e4)
            if(ipri.ge.3) then
               write(6,'(" (iv,ib) = (",2i5,") e(1:4) = ",4f9.5
     &              , " ieb(1:4) = ",4i3)') iv,ib,e1,e2,e3,e4
     &              ,ieb(1),ieb(2),ieb(3),ieb(4)
            end if
            if(e1 > ee) cycle
            ns = (e1 - es)/deltae
            if(es+deltae*ns .lt. e1) ns = ns + 1
            ne = (e4 - es)/deltae+1
            if(es+deltae*ne .ge. e4) ne = ne - 1
            if(es+deltae*ne .ge. e4) ne = ne - 1

            iloop = iloop + (ne-ns+1)
            ne1 = (e2-es)/deltae
            if(es+deltae*ne1>e2) ne1 = ne1-1
            ns2 = ne1+1
            ns3 = (e3-es)/deltae+1
            if(es+deltae*ns3 < e3) ns3 = ns3+1
            ne2 = ns3-1

            if(ipri .ge. 3) then
               write(6,'(" ns, ne1 = ",2i8," e = ",f9.5," - ",f9.5)')
     &              ns,ne1,es+deltae*ns,es+deltae*ne1
               write(6,'(" ns2,ne2 = ",2i8," e = ",f9.5," - ",f9.5)')
     &              ns2,ne2,es+deltae*ns2,es+deltae*ne2
               write(6,'(" ns3,ne  = ",2i8," e = ",f9.5," - ",f9.5)')
     &              ns3,ne,es+deltae*ns3,es+deltae*ne
               if(es+deltae*ns <e1)write(6,'(" !! es+deltae+ns  < e1")')
               if(es+deltae*ne1>e2)write(6,'(" !! es+deltae*ne1 > e2")')
               if(es+deltae*ns2<e2)write(6,'(" !! es+deltae*ns2 < e2")')
               if(es+deltae*ne2>e3)write(6,'(" !! es+deltae*ne2 > e3")')
               if(es+deltae*ns3<e3)write(6,'(" !! es+deltae*ns3 < e3")')
               if(es+deltae*ne >e4)write(6,'(" !! es+deltae*ne  > e4")')
            end if

            call gettod(wct_now)
            etime0 = etime0 + (wct_now-wct_start)*1.d-6

            call gettod(wct_start)

            d21=e2-e1
            d31=e3-e1
            d41=e4-e1
            d32=e3-e2
            d42=e4-e2
            d43=e4-e3
            ncounter1 = ncounter1+(ne1-ns+1)
            do idos = ns, ne1
               e = es + deltae*idos
               d1=e-e1
               d4=e4-e
               d2=e2-e
               d3=e3-e
               yy=d41*d31*d21
               x=d2/d21+d3/d31+d4/d41
               y=(d1*d1)/yy
               dos(idos,1)=x*y
               dosin(idos,1)=0.25*d1*y*(x+1.0)
               xx=d1*d1*d1
               x=xx/(d21*yy)
               dos(idos,2)=x
               dosin(idos,2)=0.25*d1*x
               x=xx/(d31*yy)
               dos(idos,3)=x
               dosin(idos,3)=0.25*d1*x
               x=xx/(d41*yy)
               dos(idos,4)=x
               dosin(idos,4)=0.25*d1*x
            end do
            call gettod(wct_now)
            etime1 = etime1 + (wct_now-wct_start)*1.d-6

            call gettod(wct_start)
            ncounter3 = ncounter3+(ne-ns3+1)
            do idos = ns3,ne
               e = es + deltae*idos
               d1=e-e1
               d4=e4-e
               d2=e-e2
               d3=e-e3
               xx=d4*d4*d4
               yy=d41*d42*d43
               x=xx/(d41*yy)
               dos(idos,1)=x
               dosin(idos,1)=0.25*(1.0-d4*x)
               x=xx/(d42*yy)
               dos(idos,2)=x
               dosin(idos,2)=0.25*(1.0-d4*x)
               x=xx/(d43*yy)
               dos(idos,3)=x
               dosin(idos,3)=0.25*(1.0-d4*x)
               x=d3/d43+d2/d42+d1/d41
               y=(d4*d4)/yy
               dos(idos,4)=x*y
               dosin(idos,4)=0.25*(1.0-d4*y*(x+1.0))
            end do
            call gettod(wct_now)
            etime3 = etime3 + (wct_now-wct_start)*1.d-6

            call gettod(wct_start)
            ncounter2 = ncounter2+(ne2-ns2+1)
            td21 = 3.0*d21
            yy3=1.0/(d31*d42)+1.0/(d41*d32)
            do idos=ns2,ne2
               e = es + deltae*idos
c$$$               e = es+deltae*(idos-0.5)
               d1=e-e1
               d4=e4-e
               d2=e-e2
               d3=e3-e
               d2d2 = d2*d2
               d2d2d2 = d2d2*d2
               d3d3 = d3*d3
               x1=(d3d3)/(d31*d31*d32)*(d2/d42+d1/d41)
               x2=(d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
               x3=(d3*d4*d1)/(d31*d41)*yy3
               dos(idos,1)=0.5*(x1+x2+x3)
               x1=d2d2*(d32*d2+3.0*d3*(d32+d3))/12.0
               y1=x1
               x1=x1/(d31*d31*d32*d42)
               x2=d2*(d2d2*(d31+td21)+3.0*d3*(d2*d3+d32*(td21+d1)))
               x2=x2/12.0
               y2=x2
               x2=x2/(d31*d31*d32*d41)
               x3=d2d2*(d42*d2+3.0*d4*(d42+d4))/12.0
               y3=x3
               x3=x3/(d41*d41*d42*d32)
               x4=d2*(d2d2*(d41+td21)+3.0*d4*(d2*d4+d42*(td21+d1)))
               x4=x4/12.0
               y4=x4
               x4=x4/(d41*d41*d42*d31)
               x5=0.5*d2*d3*d4*(d1+d21)
               x5=x5+d2d2*(2.0*d21*(d3+d42)
     &              +(d1+d21)*(2.0*d3+d4+d42))/12.0
               x5=x5*yy3/(d31*d41)
               x6=0.25*d21*d21*(d42/d41+d32/d31+1.0)/(d41*d31)
               dosin(idos,1)=0.5*(x1+x2+x3+x4+x5)+x6

               x1=(d3d3)/(d32*d32*d31)*(d2/d42+d1/d41)
               x2=(d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
               x3=(d3*d4*d2)/(d42*d32)*yy3
               dos(idos,2)=0.5*(x1+x2+x3)
               x1=y1/(d32*d32*d31*d42)
               x2=y2/(d32*d32*d31*d41)
               x3=y3/(d42*d42*d41*d32)
               x4=y4/(d42*d42*d41*d31)
               x5=d2d2*(d3*(d42+3.0*d4)+d32*(d42+d4))/12
               x5=x5*yy3/(d42*d32)
               x6=0.25*d21/d31*d21/d41
               dosin(idos,2)=0.5*(x1+x2+x3+x4+x5)+x6

               x1=(d2d2)/(d32*d32*d42)*(d3/d31+d4/d41)
               x2=(d1*d1)/(d31*d31*d41)*(d3/d32+d4/d42)
               x3=(d1*d2*d3)/(d32*d31)*yy3
               dos(idos,3)=0.5*(x1+x2+x3)
               x1=d2d2d2*(3*d3+d32)/12
               y1=x1
               x1=x1/(d32*d32*d42*d31)
               x2=d2d2d2*(3*d4+d42)/12
               y2=x2
               x2=x2/(d32*d32*d42*d41)
               x3=d2*(d2*d31*(d2+3*d21)+3*d3*(d2d2+3*d21*d1)+
     &              3*d21*d21*d32)/12
               y3=x3
               x3=x3/(d31*d31*d41*d32)
               x4=d2*(d2*d41*(d2+3*d21)+3*d4*(d2d2+3*d21*d1)+3*
     &              d21*d21*d42)/12
               y4=x4
               x4=x4/(d31*d31*d41*d42)
               x5=(d2d2)*(d3*(d21+3*d1)+d32*(d21+d1))/12
               x5=x5/(d32*d31)*yy3
               x6=0.25*(d21/d31*d21/d31*d21/d41)
               dosin(idos,3)=0.5*(x1+x2+x3+x4+x5)+x6

               x1=(d2d2)/(d42*d42*d32)*(d3/d31+d4/d41)
               x2=(d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
               x3=(d1*d2*d4)/(d41*d42)*yy3
               dos(idos,4)=0.5*(x1+x2+x3)
               x1=y1/(d42*d42*d32*d31)
               x2=y2/(d42*d42*d32*d41)
               x3=y3/(d41*d41*d31*d32)
               x4=y4/(d41*d41*d31*d42)
               x5=(d2d2)*(d4*(d21+3*d1)+d42*(d21+d1))/12
               x5=x5/(d41*d42)*yy3
               x6=0.25*(d21/d31*d21/d41*d21/d41)
               dosin(idos,4)=0.5*(x1+x2+x3+x4+x5)+x6
            end do
            call gettod(wct_now)
            etime2 = etime2 + (wct_now-wct_start)*1.d-6
            call gettod(wct_start)
            if(idim .eq. -3) then
               esum = e1+e2+e3+e4
               do idos = ns, ne
                  tdos = 0.025d0*(dos(idos,1)+dos(idos,2)
     &                 +dos(idos,3)+dos(idos,4))
                  dosin(idos,1)=dosin(idos,1)+tdos*(esum-4.d0*e1)
                  dosin(idos,2)=dosin(idos,2)+tdos*(esum-4.d0*e2)
                  dosin(idos,3)=dosin(idos,3)+tdos*(esum-4.d0*e3)
                  dosin(idos,4)=dosin(idos,4)+tdos*(esum-4.d0*e4)
               end do
            end if
            do ip = 1,4
               do ie = ns, ne
                  cdos(ie) = cdos(ie)
     &                 +dos(ie,ip)*dos_weight(ib,ieb(ip))
                  csumdos(ie) = csumdos(ie)
     &                 + dosin(ie,ip)*dos_weight(ib,ieb(ip))
               end do
            end do
            w = (dos_weight(ib,ieb(1))+dos_weight(ib,ieb(2))
     &           +dos_weight(ib,ieb(3))+dos_weight(ib,ieb(4)))*0.25d0
            do ie = ne+1,nEwindows
               csumdos(ie) = csumdos(ie) + 1.d0*w
            end do
            call gettod(wct_now)
            etime4 = etime4 + (wct_now-wct_start)*1.d-6
         end do
      end do
      call gettod(wct_start)
      do ie = 0, nEwindows
         cdos(ie) = cdos(ie)/ntet
         csumdos(ie) = csumdos(ie)/ntet
      end do

      if (npes .gt. 1) then
         itmp = (nep+1)*4
         allocate(dos_mpi(0:nep,4));dos_mpi=0.d0
         allocate(dosin_mpi(0:nep,4));dosin_mpi=0.d0
         allocate(cdos_mpi(0:newindows));cdos_mpi=0.d0
         allocate(csumdos_mpi(0:newindows));csumdos_mpi=0.d0

         call mpi_allreduce(dos,dos_mpi,itmp,mpi_double_precision,
     &   mpi_sum,MPI_CommGroup,ierr)
         dos = dos_mpi

         call mpi_allreduce(dosin,dosin_mpi,itmp,mpi_double_precision,
     &   mpi_sum,MPI_CommGroup,ierr)
         dosin = dosin_mpi

         call mpi_allreduce(cdos,cdos_mpi,(newindows+1),
     &   mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         cdos = cdos_mpi

         call mpi_allreduce(csumdos,csumdos_mpi,(newindows+1),
     &   mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         csumdos = csumdos_mpi

         deallocate(dos_mpi);deallocate(dosin_mpi)
         deallocate(cdos_mpi)
         deallocate(csumdos_mpi)
      endif

      call gettod(wct_now)
      etime5 = etime5 + (wct_now-wct_start)*1.d-6

      if(ipri .ge. 2) then
         write(6,'(" !dos iloop    = ",i12)') iloop
         write(6,'(" !dos ncounter1 = ",i12)') ncounter1
         write(6,'(" !dos ncounter2 = ",i12)') ncounter2
         write(6,'(" !dos ncounter3 = ",i12)') ncounter3
         write(6,'(" !dos etime0    = ",f16.8)') etime0
         write(6,'(" !dos etime1    = ",f16.8)') etime1
         write(6,'(" !dos etime2    = ",f16.8)') etime2
         write(6,'(" !dos etime3    = ",f16.8)') etime3
         write(6,'(" !dos etime4    = ",f16.8)') etime4
         write(6,'(" !dos etime5    = ",f16.8)') etime5
      end if
      return
      end
