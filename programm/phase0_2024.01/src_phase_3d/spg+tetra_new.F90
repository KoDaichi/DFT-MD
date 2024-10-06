!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  The spg+tetra subroutine file(Fortran90 version 5.3)
!
! Last Updetaed 2022/8/15
!

!  SUBROUINES:
!     fermi1, fermi2, nsdos0, nsdos2, nsdos3, nsdosi, nseulc
!     nseulh, nsgcm2, nsgcm3, nsgrp1, nsgrpb, nsgrpb_kt, nsgrpa
!     nsjonh, nskma0, nskp00, nskpb0, nskpb0_s, nskpa0_kt, nskpbm
!     nslat3, nslata, nslata_kt, nslatb, nslatc, nslatr, nslatz,
!     nsmetr, nsmlt1, nsmult, nspace, nspgrp, nspgrp_kt, nsprmv, nsrduc
!     nsrmxc, nsrmxh, nsrot1, nsrota, nsrotc, nsroth, nsrotk
!     nsrotl, nsrotr, nssdjg, nssum1, nstrsh, nstt0i, nstt1i
!     nstt2i, nstt3i, nsttod, nstts1, rdprp, setkp0, setkp0_n
!     setkp0_n_kt, ka00, setkp0_default, setkp0_default_n,
!     setkp0_default_n_kt2, setkp0_default_n_kt, setspg,
!     setspg_n, setspg_default, setspg_default_n,
!    setspg_default_n_kt, getspgtab, tbspg, bspg_kt, check_if_unit_matrix,
!    tspaca, wtetra, nstt4i, nstt5i
!
!  HISTORY of DEVELOPMENTS and AUTHORS:
!    1. The original spg+tetra subroutine codes were written in Fortran77 by
!       A. Yanase, and K. Terakura  June/07/1986
!    2. The codes were modifined for PHASE/0 by
!        N. Hamada, H. Mizouchi, K. Mae  August/20/2003
!    3. The codes were further modified by
!        T. Yamasaki  September/16/2007
!    4. The codes were re-written in Fortran90 and were modified  by
!       T. Hamada  January/23/2022
!
!  REMARKS and NOTES
!    1. This Fortran90 file contains 76 spg+tetra subroutines used
!       by PHASE/0.
!    2. This file  must be complied along with m_Spg_plus_Tetra_Common_Blocks.F90.
!    3. The following subroutines are largely modified for the performance
!       improvement of PHASE/0: fermi2, nsdos2, nsdos0, nsdosi, nstt2i,
!       nstt1i nstt0i, nskpb0, nskpb0_s, nskpb0_kt, nskpbm.
!    4. Some subroutines have OpenMP directives for paralllel compiler.
!    5. The spg+tetra subroutines not used by PHASE/0 are eliminated.
!    6. The original Fortran77 codes of the spg+tetra subroutines are
!       in spg+tetra_F77.F.
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
! ========================= modified by K. Tagami ================== 11.0
!     subroutine fermi1(nxx,nyy,nzz,np2,lmnp2,neig,lmneig,nspin,
!     &                  eig2,ip20,np0,valenc,eferm,eband,valud,
!     &                  iwt,ip2cub,ipri )
!
      subroutine fermi1(nxx,nyy,nzz,np2,lmnp2,neig,lmneig,nspin, &
     &                  eig2,ip20,np0,valenc,eferm,eband,valud, &
     &                  iwt,ip2cub,ipri, flag_noncol )
! =================================================================== 11.0

!     $Id: spg+tetra.F 633 2020-12-01 05:11:03Z jkoga $
!
!      eferm :  fermi energy
!      dos   :  density of states at fermi energy for each spin
!      sidos :  integrated dos at fermi energy for each spin
!      ddd = sum of dos
!      sss = sum of sidos
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in)      :: nxx, nyy, nzz, np2, lmnp2, neig, lmneig &
                                   & ,nspin, np0, ipri
         integer, intent(in), dimension(np0)  ::  ip20
         integer, intent(in), dimension(np2)  ::  iwt
         integer, intent(in), dimension(nxx*nyy*nzz) :: ip2cub
         real(DP), intent(inout) :: eferm, eband
         real(DP), intent(in), dimension(lmnp2,lmneig,nspin) :: eig2
         integer :: idim, ispin, neig1, neig2, is, i, ip2, iii, instts1
         real(DP)  :: emin, emax, wei, efermi, e1, e2, ddd, sss, dm1, smom , ebnd
         real(DP), dimension(nspin)  :: dos, sidos
         real(DP), dimension(np0) :: eawk
         real(DP), dimension(nspin) :: valud
         real(DP), dimension(200)  :: efm, qfm
         real(DP) :: eps = 1.0d-10
! ================================== added by K. Tagami ================- 11.0
         real(DP)  :: valenc
! === Optional argument needs interface!!! by T.Kato 2013/07/02 ================
!     logical, optional:: flag_noncol
         logical  :: flag_noncol
! ==============================================================================
! ====================================================================== 11.0

#ifdef __TIMER_SUB__
         call timer_sta(710)
#endif

         idim = 3
!     eps=dfloat(10)**(-6)
!         eps = 1.0d-10
         neig1 = 1
         neig2 = neig
      !!$if(np0.gt.20000) then
      !!$  write(6,*) ' np0=',np0,' > 20000 : error at sub.fermi1.'
      !!$  stop 'error at sub.fermi1 (np0).'
      !!$end if
         emin = eig2(1,1,1)
         emax = eig2(1,1,1)
#ifdef __TIMER_DO__
         call timer_sta(811)
#endif
         do  is=1, nspin
            do  i=1, neig
               do  ip2=1,np2
                 if(emin > eig2(ip2,i,is)) emin = eig2(ip2,i,is)
                 if(emax < eig2(ip2,i,is)) emax = eig2(ip2,i,is)
               end do
            end do
         end do
#ifdef __TIMER_DO__
         call timer_end(811)
#endif
         if(ipri >= 2) then
            if(printable) write(nfout,*) ' emax=',emax,'   emin=',emin
         end if
         wei=2.d0
         if(nspin == 2) wei=1.d0

! ===========================- added by K. Tagami ==================== 11.0
         if ( flag_noncol ) wei = 1.0d0
! ==================================================================== 11.0

! modified by H.Sawada on May 1, 1997
         efermi = eferm
! modified by H.Sawada on May 1, 1997
         e1     = emin
         e2     = emax
!
         iii = 0
         instts1 = 0
 55      continue
 33      continue
         ddd = 0.d0
         sss = 0.d0
         ebnd = 0.d0
#ifdef __TIMER_DO__
         call timer_sta(812)
#endif
         do ispin=1, nspin
           call nsdos2(idim,efermi,nxx,nyy,nzz,lmnp2,neig1,neig2, &
        &  eig2(1,1,ispin), &
        &  ip20,np0,eawk,dos(ispin),sidos(ispin),dm1, &
        &  instts1,np2,iwt,ip2cub)
           ddd=ddd+dos(ispin)*wei
           sss=sss+sidos(ispin)*wei
           ebnd=ebnd+dm1*wei
           valud(ispin)=sidos(ispin)*wei
         end do
#ifdef __TIMER_DO__
         call timer_end(812)
#endif
         iii=iii+1
         if(iii == 100 .and. instts1==0) then
            instts1 = 1
            efermi = eferm
            emin = eig2(1,1,1)
            emax = eig2(1,1,1)
#ifdef __TIMER_DO__
            call timer_sta(813)
#endif
            do  is=1, nspin
               do  i=1, neig
                  do  ip2=1, np2
                    if(emin > eig2(ip2,i,is)) emin = eig2(ip2,i,is)
                    if(emax < eig2(ip2,i,is)) emax = eig2(ip2,i,is)
                  end do
               end do
            end do
#ifdef __TIMER_DO__
            call timer_end(813)
#endif
            e1     = emin
            e2     = emax
            iii = 0
            goto 55
         elseif(iii > 100) then
            if(printable) write(nfout,*) 'iii=',iii,'   efermi=',efermi
            if(printable) write(nfout,*) 'sss=',sss,'  val=',valenc
            do  i=1,100
               if(printable) write(nfout,*) ' efermi=',efm(i),'   q=',qfm(i)
            end do
            stop ' === stop sub.fermi1. (iii>200) ==='
         end if

         efm(iii) = efermi
         qfm(iii) = sss
         if(dabs( sss-valenc)< eps) goto 34
!       if(abs((sss-valenc)/valenc).lt.eps) go to 34
!
         if(sss < valenc) then
            e1     = efermi
            efermi = efermi + (e2-efermi)*0.50d0
         else
            e2     = efermi
            efermi = efermi + (e1-efermi)*0.50d0
         end if
         go to 33

   34    continue

         if(dabs(ddd) < 1.d-20) then
            call fermi2(np2,lmnp2,lmneig,neig,nspin,eig2,efermi,ipri)
         end if
!      write( 6,111) efermi,ddd,sss
!      write(16,111) efermi,ddd,sss
!  111 format(5x,'fermi-en=',f12.6,5x,'dos=',f14.6,5x,'int.dos=',f14.6)
! --*
         if(nspin == 2) then
            smom=sidos(1)-sidos(2)
!        write( 6,120) smom
!        write(16,120) smom
         end if
!  120 format(' magnetic moment=',f12.6,' per cell')
         eferm=efermi
         eband=ebnd

#ifdef __TIMER_SUB__
         call timer_end(710)
#endif
      end subroutine fermi1
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine fermi2(np2,lmnp2,lmneig,neig,nspin,eig2,eferm,ipri)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: np2, lmnp2, lmneig, neig, nspin, ipri
         integer :: i, ip2, is
         real(DP), intent(in), dimension(lmnp2,lmneig,nspin) :: eig2
         real(DP), intent(inout) :: eferm
         real(DP) :: eps, zero, e1, e2, e0

#ifdef __TIMER_SUB__
         call timer_sta(712)
#endif
         eps = 1.0d-20
         e1 = -1.0d4
         e2 =  1.0d4
#ifdef __TIMER_DO__
         call timer_sta(815)
#endif
         do  is=1, nspin
            do  i=1, neig
               do ip2=1, np2
                  e0 = eig2(ip2,i,is)-eferm
                  if(e0 < 0.0d0 .and. e0 > (e1-eferm)) then
                     e1 = eig2(ip2,i,is)
                  else if(e0 > 0.0d0 .and. e0 <(e2-eferm)) then
                     e2=eig2(ip2,i,is)
                  end if
               end do
            end do
         end do
#ifdef __TIMER_DO__
        call timer_end(815)
#endif
!     eferm= e1+eps
        eferm=(e1+e2)*0.5d0
!     eferm= e2-eps
        if(ipri >= 2) then
           if(printable) then
              write(nfout &
     &        ,'('' e(val.top)='',f12.6,5x,''e(cond.bottom)='',f12.6)') &
     &       e1,e2
           end if
        end if
#ifdef __TIMER_SUB__
        call timer_end(712)
#endif
      end subroutine fermi2
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nsdos0(e,nxx,nyy,nzz,ea,dos,dm0,dm1, &
          & lmnp2,eig2,instts1,np0,np2,iwt,ip2cub)

!     ***  by linear interporation in tetrahedrons
!     ***
!     ***
!     ***  nxx  number of mesh points in x-direction
!     ***  nyy  number of mesh points in y-direction
!     ***  nzz  number of mesh points in z-direction
!     ***
!     ***  dos  density of states
!     ***  dm0  0th moment of dos (interated up to e)
!     ***  dm1  1st moment of dos (interated up to e)
!
!     This subroutine is modified for PHASE/0 by T. Hamada
!     December/20/2021
!
         use m_Const_Parameters, only :  DP
         use m_Control_Parameters, only : width_tetra
         use m_Spg_plus_Tetra_Common_Blocks, only : npx, npy, npz, np, &
         & ntet, ncub, rntet_1, tetra_eps, ip0_index1, ip0_index2, &
         &  tet_sub_called
         use m_Parallelization, only : npes,mype, MPI_CommGroup, &
       & mpi_double_precision, mpi_integer, istatus, ierr
         use m_Timing, only : tstatc0_begin, tstatc0_end
         implicit none
!         include 'mpif.h'
          integer, intent(in) :: nxx, nyy, nzz, np0,np2, lmnp2,  instts1
          integer, intent(in), dimension(np2) :: iwt
          integer, intent(in), dimension(nxx*nyy*nzz) :: ip2cub
          integer :: ni, icub, ip, ip0, iq, it, ix, iy, iz, kx, ky, kz, m
          integer :: ip2cub_1
          integer, dimension(2,2,2) :: iecub
          integer, dimension(8) :: iec
          integer, dimension(4) :: iet, ieb
          integer, dimension(6,2) :: iqmat
          real(DP) , intent(in) :: e
          real(DP), intent(in), dimension(np0) :: ea
          real(DP), intent(in), dimension(lmnp2) :: eig2
          real(DP), intent(out) :: dos, dm0, dm1
          real(DP) :: d, d0, d1, e1, e2, e3, e4, emax, emin, eps2
          real(DP) :: emin0 =  1.0d30
          real(DP) :: emax0 = -1.0d30
          real(DP) :: eps = 1.0d-4
          real(DP), dimension(2,2,2) :: ecub
          real(DP), dimension(8) :: ec
          real(DP), dimension(4) :: et, eb

          equivalence(ec(1),ecub(1,1,1))
          equivalence(iec(1),iecub(1,1,1))

          data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/

!  commented by K.Mae 2003.8.5
!     real*8 doscub(5000),dm0cub(5000),dm1cub(5000)
!   added by K.Mae 2003.8.5
         real(DP), dimension(nxx*nyy*nzz) :: doscub, dm0cub, dm1cub
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         integer :: id_sname = -1

         call tstatc0_begin('nsdos0 ', id_sname)

!  definition of eps  <- must be consistent with <nstts1>
!*    eps=dfloat(10)**(-4)
!        eps = 1.0d-4
         eps = tetra_eps !_ defined in m_Spg_plus_Tetra_Common_Blocks
!         eps = width_tetra*0.5d0 ! can be configured from input
         eps2 = eps*2.0d0

         if(.not. tet_sub_called) then
            npx = nxx+1
            npy = nyy+1
            npz = nzz+1
            np = npx*npy*npz
            ncub = nxx*nyy*nzz
            ntet = 6*ncub
            rntet_1 = 1.0d0/real(ntet,kind=DP)

            allocate(ip0_index1(ncub))

            icub = 0
            do iz=0, nzz-1
               do iy=0, nyy-1
                  do ix=0, nxx-1
                     icub = icub+1
                     ip0_index1(icub)  = npx*(npy*iz+iy)+ix
                  end do
               end do
            end do

            do  kz=1, 2
               do ky=1, 2
                  do  kx=1, 2
                     ip0_index2(kx,ky,kz) = npx*(npy*(kz-1)+ky-1)+kx
                  end do
               end do
            end do
            if(npes > 1) then
               call mpi_bcast(npx,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(npy,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(npz,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(np,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ncub,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ntet,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(rntet_1,1,mpi_double_precision,0,MPI_CommGroup,ierr)
               call mpi_bcast(ip0_index1,ncub,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ip0_index2,8,mpi_integer, 0,MPI_CommGroup,ierr)
            end if
            tet_sub_called = .true.
         end if

!!!! commented by K. Mae 2003.8.5 !!!!
!!      if(ncub.gt.5000) then
!!         write(6,*) ' ncub= ',ncub,' > 5000'
!!         write(6,*) ' : error at sub.nsdos0.'
!!         stop 'error at sub.nsdos0 (ncub).'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         dos = 0.0d0
         dm0 = 0.0d0
         dm1 = 0.0d0

         do ip=1, np
            if(ea(ip) > emax0) emax0 = ea(ip)
            if(ea(ip) < emin0) emin0 = ea(ip)
         end do
!                                   ============ if 1 ==
!**    if(e.gt.emin0-eps*2) then
         if(e > emax0 .and. instts1 == 0) then
            do  ip = 1, np2
               dm1 = dm1 + eig2(ip) * real(iwt(ip),kind=DP)
            end do
            dm1 = dm1 / real(ntet,kind=DP) / 4.d0
            dm0 = 1.0d0
         elseif(e > (emin0-eps2)) then

!     ***  integration over b.z. starts    ***
!
!     ***       sampling over cubes        ***

            icub = 0
            do iz=0, nzz-1
               do iy=0, nyy-1
                  do ix=0, nxx-1
                     icub = icub+1
                     if(icub /= ip2cub(icub)) then
                        dos = dos+doscub(ip2cub(icub))
                        dm0 = dm0+dm0cub(ip2cub(icub))
                        dm1 = dm1+dm1cub(ip2cub(icub))
                     else
!     ***  energies at cube corners  ***
                        doscub(icub) =  0.0d0
                        dm0cub(icub) =  0.0d0
                        dm1cub(icub)  = 0.0d0
!!                      emax =-1.0d30
                        emin = 1.0d30
                        do  kz=1, 2
                           do ky=1, 2
                              do  kx=1, 2
                                 ip0 = ip0_index1(icub) + &
                                  & ip0_index2(kx,ky,kz)
                                  ecub(kx,ky,kz) = ea(ip0)
                                  iecub(kx,ky,kz) = ip0
!!                               if(ea(ip0) > emax) emax = ea(ip0)
                                 if(ea(ip0) < emin) emin = ea(ip0)
                              end do
                           end do
                        end do
!                                   ============ if 2 ==
!                        if(e > (emin-eps2)) then
                        if(e > (emin-eps2)) then
!**      if(e.gt.emin) then
!         ***      six tetrahedrons      ***
!         *** sampling over tetrahedrons ***
                           et(1) = ec(1)
                           et(4) = ec(8)
                           iet(1) = iec(1)
                           iet(4) = iec(8)
                           do it=1, 6
                              iq = iqmat(it,1)
                              et(2) = ec(iq)
                              iet(2) = iec(iq)
                              iq = iqmat(it,2)
                              et(3) = ec(iq)
                              iet(3) = iec(iq)
                              eb(1:4) = et(1:4)
                              ieb(1:4) = iet(1:4)
!           ***  eb(1).le.eb(2).le.eb(3).le.eb(4)  ***
                              call nsttod(eb,ieb)
                              e1 = eb(1)
                              e2 = eb(2)
                              e3 = eb(3)
                              e4 = eb(4)
!                                   ============ if 3 ==
!$$$           if(e.ge.e4) then
!$$$               dm0=dm0+1
!$$$               dm1=dm1+(e1+e2+e3+e4)/4.d0
!$$$           else if(e.gt.e1) then
                              call nstts1(e1,e2,e3,e4)
                              d=0.0d0
                              d0 = 0.0d0
                              d1 = 0.0d0
                              if(e <= e1) cycle
                              if(e >= e4) then
                                 d  = 0.0d0
                                 d0 = 1.d0
                                 d1=(e1+e2+e3+e4)*0.25d0
                              else
                                 call nsdosi(e,e1,e2,e3,e4,d,d0,d1)
                              end if
                              dos = dos+d
                              dm0 = dm0+d0
                              dm1 = dm1+d1
                              doscub(icub) = doscub(icub)+d
                              dm0cub(icub) = dm0cub(icub)+d0
                              dm1cub(icub) = dm1cub(icub)+d1
!$$$           end if
!                                   ============ if 3 ==
                           end do
                       end if
                     end if
!                                   ============ if 2 ==
                  end do
               end do
            end do
            dos=dos*rntet_1
            dm0=dm0*rntet_1
            dm1=dm1*rntet_1
         end if

         call tstatc0_end(id_sname)
      end subroutine nsdos0
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsdos2(idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2, &
         &ip20,np0,ea,dos,dm0,dm1,instts1,np2,iwt,ip2cub)
         use m_Const_Parameters, only :  DP, CMPLDP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
!
! This subroutine was modified for PHASE/ by T. Gamada
! December/20/2021
!
         implicit none
         integer, intent(in) :: idim, nx, ny, nz, np0, np2, lmnp2, neig1, neig2, instts1
         integer, intent(in), dimension(np0) :: ip20
         integer, intent(in), dimension(np2) :: iwt
         integer :: ieig, k0
         integer, intent(in), dimension(nx*ny*nz) :: ip2cub
         real(DP), intent(in) :: e
         real(DP), intent(in), dimension(lmnp2, neig2) :: eig2
         real(DP), intent(out) :: dos, dm0, dm1
         real(DP), intent(out),dimension(np0) :: ea
         real(DP) :: d, d0, d1


#ifdef __TIMER_SUB__
         call timer_sta(711)
#endif
         dos = 0.0d0
         dm0 = 0.0d0
         dm1 = 0.0d0
#ifdef __TIMER_DO__
         call timer_sta(814)
#endif
         do  ieig=neig1, neig2
            ea(1:np0) = eig2(ip20(1:np0),ieig)
!       write(6,*) ' k0,k1=',k0,ip20(k0),'   ea=',ea(k0)
            call nsdos0(e,nx,ny,nz,ea,d,d0,d1, &
              & lmnp2,eig2(1,ieig),instts1,np0,np2,iwt,ip2cub)
            dos = dos + d
            dm0 = dm0 + d0
            dm1 = dm1 + d1
         end do
#ifdef __TIMER_DO__
           call timer_end(814)
#endif

#ifdef __TIMER_SUB__
           call timer_end(711)
#endif
      end  subroutine nsdos2
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
! ====================== modified by K. Tagami ==================== 11.0
!      subroutine nsdos3(jf,idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2,
!     &                  ip20,np0,ea,instts1,np2,iwt,ip2cub)
      subroutine nsdos3(jf,idim,e,nx,ny,nz,lmnp2,neig1,neig2,eig2, &
             & ip20,np0,ea,instts1,np2,iwt,ip2cub, ipri )
! ================================================================== 11.0
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
!$ use omp_lib
         implicit none
         integer, intent(in) :: jf, idim, nx, ny, nz, np0, np2, lmnp2, neig1, neig2, instts1
         integer, intent(in), dimension(np0) :: ip20
         integer, intent(in), dimension(np2) :: iwt
         integer :: ieig, k0, iflag, i0, j , j1
         integer, intent(in), dimension(nx*ny*nz) :: ip2cub
         real(DP), intent(in) :: e
         real(DP), intent(in), dimension(lmnp2, neig2) :: eig2
         real(DP), intent(out),dimension(np0) :: ea
         real(DP) :: one, eps, d, d0, d1
         real(DP), dimension(100) :: dos, dm0

! ==================== added by K. Tagami ========= 11.0
         integer, intent(in) ::  ipri
! ================================================= 11.0

         one = 1.0d0
         eps = 1.0d-6
         iflag = 0

         do ieig=neig1, neig2
            ea(1:np0)=eig2(ip20(1:np0),ieig)
!$$$        write(6,*) ' k0,k1=',k0,ip20(k0),'   ea=',ea(k0)
               call nsdos0(e,nx,ny,nz,ea,d,d0,d1, &
              & lmnp2,eig2(1,ieig),instts1,np0,np2,iwt,ip2cub)
            if(iflag == 0) then
               if((1.0d0-d0) < eps) then
                  dos(1) = d
                  dm0(1) = d0
               else
                  dos(2) = d
                  dm0(2) = d0
                  i0=ieig - 2
                  j =2
                  iflag=1
               end if
            else
               j = j + 1
               dos(j) = d
               dm0(j) = d0
               if(d0 < eps .or. j >= 100) then
                  exit
               end if
            end if
         end do

         if(ipri >= 1) then
            j1=j
            write(jf,*) ' ---        dos for each band, each spin ---'
            write(jf,100) (i0+j,dos(j),j=1,j1)
            write(jf,*) ' --- occupation for each band, each spin ---'
            write(jf,100) (i0+j,dm0(j),j=1,j1)
            write(jf,*) ' -------------------------------------------'
            write(jf,'(" i0 = ",i8)') i0
         end if
  100    format((5(i6,f9.6)))
      end subroutine nsdos3
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nsdosi(e,e1,e2,e3,e4,dos,dm0,dm1)
!
!  ** tetrahedron corners at one energy        **  c
!  ** according to lambin and vigneron,        **  c
!  ** phys. rev. b29, 3430 (1984)              **  c
!
!      implicit real*8(a-h,o-z)
!     zero=0.d0   change by Tsuyoshi Miyazaki 94.8.28
!     This subroutine was modified for PHASE/0 by T. Hamada
!     December/20/2021
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         real(DP), intent(in) :: e, e1,e2, e3, e4
         real(DP), intent(out) :: dos, dm0, dm1
         real(DP) :: de, d21, d21d21, d31, d32, d41, d42, d43, d1, d1d1, &
        & d2, d2d2,d3, d4, dd, dd1_1, dd2_1, te4
         real(DP) :: e32, e41, o, p, p1, p2, q, q1, q2
         real(DP) :: zero = 1.0d-10

         d21 = e2-e1
         d32 = e3-e2
         d43 = e4-e3

         if(d21.le.zero .or. d32.le.zero .or. d43.le.zero) then
            if(printable) then
               write(nfout,*) ' Warning!! in = sub.nsdosi ='
               write(nfout,*) ' e1,e2,e3,e4=',e1,e2,e3,e4
               write(6,*) ' energy order error  in sub.nsdosi.'
            end if
         end if

         if(e.le.e2) then
            d31 = e3-e1
            d41 = e4-e1
            d1=e-e1
            dd=d21*d31*d41
            dos = (d1*d1)/dd
            dm0 = d1*dos
            dos = 3.0d0*dos
            dm1=dm0*(3.d0*e+e1)*0.25d0
         else if(e.lt.e3) then
            d31 = e3-e1
            d41 = e4-e1
            d42 = e4-e2
            d1=e-e1
            d2=e-e2
            d3=e3-e
            d4=e4-e
            d1d1 = d1*d1
            d2d2 = d2*d2
            d21d21 = d21*d21
            de = 2.0d0*e
            te4 = 3.0d0*e4
            e32=(3.0d0*e3-e2)*0.5d0-e
            e41=(te4-e1)*0.5d0-e
            dd1_1=1.0d0/(d42*d32*d31)
            dd2_1=1.0d0/(d42*d41*d31)
            dos=3.d0*(d3*d2*dd1_1+d4*d1*dd2_1)
            dm0=(d2d2)*e32*dd1_1+(d1d1)*e41*dd2_1-0.50d0*d32*(d21d21)*dd1_1
            o=0.25d0*(d21d21)*(3.d0*e2+e1)*d42*dd2_1
            p1=d2*(de+e-4.0d0*e3+e2)*0.25d0
            p2=e2*(de-3.0d0*e3+e2)*0.50d0
            p=-(d2d2)*(p1+p2)*dd1_1
            q1=(d2d2)*(-6.d0*d4*d1-2.d0*d2*(de-e4-e1)+d2d2)*0.25d0
            q2=e2*((d1d1)*(de-te4+e1) &
         &   +(d21d21)*(te4-2.d0*e2-e1))*0.50d0
        q=-(q1+q2)*dd2_1
        dm1=o+p+q
         else if(e.lt.e4) then
            d41 = e4-e1
            d42 = e4-e2
            d4 = e4-e
            dd=d43*d42*d41
            dos= (d4*d4)/dd
            dm0=1.d0-d4*dos
            dm1=((e1+e2+e3+e4) - (d4*dos)*(3.d0*e+e4))*0.25d0
            dos = 3.0d0*dos
         endif
      end subroutine nsdosi
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nseulc(ieuler,euler)

!#12  noinput
!#12  output: ieuler(3,48)*2*pai/4: euler angles for cubic
!#12           euler(3,48)        : euler angles (alpha,beta, gamma)
!#13  noexternal
!
!#21  to get euler angles for point-group operations of cubic lattice
!
!#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nseulc(ieuler,euler)
         use m_Const_Parameters, only :  DP, PAI2
         implicit none
         integer, intent(out), dimension(3,24) :: ieuler
         integer :: i, j
         integer, dimension(3,24) :: ie4
         real(DP), intent(out), dimension(3,24) :: euler

         data    ie4/ &
      &  0,0,0, 0,2,2, 0,2,0, 0,0,2, &
      &  0,1,1, 2,1,3, 2,1,1, 0,1,3, &
      &  1,1,2, 3,1,0, 1,1,0, 3,1,2, &
      &  0,2,1, 0,2,3, 0,1,2, 1,1,1, 2,1,0, 3,1,3, &
      &  3,1,1, 0,1,0, 1,0,0, 1,1,3, 2,1,2, 3,0,0/

        ieuler(1:3,1:24)=ie4(1:3,1:24)
        euler(1:3,1:24)=(PAI2*real(ie4(1:3,1:24),kind=DP))*0.25d0
      end subroutine nseulc
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nseulh(ieuler,euler)
!
!#12  noinput
!#12  output: ieuler(3,12)*2*pai/4: euler angles for hexagonal
!#12           euler(3,12)        : euler angles (alpha,beta, gamma)
!#13  noexternal
!
!#21  to get euler angles for hexagonal operations
!
!#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nseulh(ieuler,euler)
         use m_Const_Parameters, only :  DP, PAI
         use m_Spg_plus_Tetra_Common_Blocks, only : factor_1_3
         implicit none
         integer, intent(out), dimension(3,12) :: ieuler
         integer :: i, j
         integer, dimension(3,12) :: ie6
         real(DP), intent(out), dimension(3,12) :: euler
         data    ie6/  &
        &  0,0,0, 1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0, &
        &  0,3,0, 4,3,0, 2,3,0, 3,3,0, 1,3,0, 5,3,0/


         ieuler(1:3,1:12) = ie6(1:3,1:12)
         euler(1:3,1:12)=PAI*real(ieuler(1:3,1:12),kind=DP)*factor_1_3
      end subroutine nseulh
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsgcm2(i1,i2,m)
!
!#12  input   :  i1,i2
!#12  output  :  m     : the greatest common measure
!
!#21  to get the greatest common measure of i1 and i2
!
!#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsgcm2(i1,i2,m)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: i1, i2
         integer, intent(out) :: m
         integer :: i, j, k

         if(i1 == 0 .or. i2 == 0) then
            if(printable) then
               write(nfout,*) '=== error in sub.nsgcm2 ==='
               write(nfout,*) 'i1,i2=',i1,i2,' : They must not be zero.'
            end if
            stop '=== error in sub.nsgcm2 ==='
         else
            i = iabs(i1)
            j = iabs(i2)
    1       continue
            k=mod(j,i)
            if(k == 0) goto 2
            j = i
            i = k
            goto 1
    2       m = i
         end if
      end subroutine nsgcm2
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!c
!c#11  sub.nsgcm3(i1,i2,i3,m)
!c
!c#12  input   :  i1,i2,i3
!c#12  output  :  m     : the greatest common measure
!c#13  external:  sub.nsgcm2
!c
!c#21  to get the greatest common measure of i1, i2 and i3
!c
!c#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
!c
!c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsgcm3(i1,i2,i3,m)
         implicit none
         integer, intent(in) :: i1, i2, i3
         integer, intent(out) :: m
         integer :: m1, m2
         external nsgcm2

         call nsgcm2(i1,i2,m1)
         call nsgcm2(i2,i3,m2)
         call nsgcm2(m1,m2,m)
      end subroutine nsgcm3
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
!#11 &           schoen,jones,euler,rot,ieuler,irot,ir1234,
!#11 &           im0,iv0,ng0,ng1,ig10,ig01,jg1,
!#11 &           irotr1,irotk1,im1,iv1,movo)
!
!#12  input:      jf : output file
!#12             jpr : print control
!#12              il : lattice type
!#12             inv : parameter (0,1) for moving the origin
!#12         igen(3) : generator (rotation part)
!#12     jgen(2,3,3) : generator (nonprimitive translation part)
!#12  output:  schoen(48),jones(3,48),euler(3,24),rot(3,3,48)
!#12           ieuler(3,24),irot(3,3,48),ir1234(3,48)
!#12           im0(48,48),iv0(48),ng0
!#12           ng1,ig10(48),ig01(48),jg1(2,3,48),
!c#12           irotr1(3,3,48),irotk1(3,3,48),im1(48,48),iv1(48)
!#12           mo(2,3): translation of the origin from the initial
!#13  external:  nsrotl, nsmult, nspgrp
!
!#21  to get space group
!
!#31  1990.01.09.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!c
!c ======================= modified by K. Tagami ===================== 12.0A
!c      subroutine nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
!c     &                  schoen,jones,euler,rot,ieuler,irot,ir1234,
!c     &                  im0,iv0,ng0,ng1,ig10,ig01,jg1,
!c     &                  irotr1,irotk1,im1,iv1,movo)
      subroutine nsgrp1(jf,jpr,il,ngen,inv,igen,jgen, &
                      & schoen,jones,euler,rot,ieuler,irot,ir1234, &
                      & im0,iv0,ng0,ng1,ig10,ig01,jg1, &
                      & irotr1,irotk1,im1,iv1,movo, &
                      & use_trs, tac, tca, tab, tba, &
                      & gen_name_in_carts )
!==================================================================== 12.0A
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: jf, jpr, il, ngen, inv
         integer, intent(in), dimension(3) :: igen
         integer, intent(in), dimension(2,3,3) :: jgen
         integer, intent(out) :: ng0, ng1
         integer, intent(out), dimension(48) :: ig10, ig01, iv1
         integer, dimension(48) :: iv0
         integer, intent(out), dimension(2,3,48) :: jg1
         integer, intent(out), dimension(3,3,48) ::irotr1, irotk1
         integer, intent(out), dimension(48,48) :: im0, im1
         integer, intent(out), dimension(3,48) :: ir1234
         integer, intent(out), dimension(3,24) :: ieuler
         integer, intent(out), dimension(3,3,48) :: irot
         integer, intent(out), dimension(2,3) :: movo
         integer :: i, i1, i2, ii, j, k, no, mo
         real(DP), intent(out), dimension(3,24) :: euler
         real(DP), intent(out), dimension(3,3,48) :: rot
         character(len =2), intent(out), dimension(3,48) :: jones
         character(len =5), intent(out), dimension(48) :: schoen
! ========================== added by K. Tagami ====================== 12.0A
         logical use_trs, gen_name_in_carts
         real(DP), dimension(3,3) ::  tac, tca, tab, tba
! ==================================================================== 12.0A

         call nsrotl(il,jpr,jf,ng0,ieuler,euler,irot,rot,ir1234,jones, &
                  & schoen)

         call nsmult(jf,jpr,ng0,irot,im0,iv0)

! ======================= modified by K. Tagami ===================== 12.0A
!      call nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im0,irot,schoen,
!     &            movo,ng1,ig01,jg1)
!
         call nspgrp_kt(jf,jpr,il,inv,ngen,igen,jgen,im0,irot,schoen, &
                      & movo,ng1,ig01,jg1, &
                      & use_trs, tac, tca, tab, tba, &
                      & gen_name_in_carts )
! ==================================================================== 12.0A

         irotr1(1:3,1:3,1:ng1) = irot(1:3,1:3,ig01(1:ng1))

         if(il <= 0) then
            do k=1, ng1
               irotk1(1:3,1:3,k) = irot(1:3,1:3,ig01(k)+24)
            end do
         else
            irotk1(1:3,1:3,1:ng1) = irot(1:3,1:3,iv0(ig01(1:ng1)))
         end if

         ig10(1:48) = 0

         do i=1, ng1
            ig10(ig01(i)) = i
         end do

         do i= 1, ng1
            do j = 1, ng1
               im1(i,j) = ig10(im0(ig01(i),ig01(j)))
            end do
         end do

         iv1(1:ng1) = ig10(iv0(ig01(1:ng1)))

         if(jpr >= 0) then
            write(jf,*) ' '
            write(jf,*) ' ----- group elements -------'
            do i=1, ng1
               write(jf,100) i,ig01(i),schoen(ig01(i)),(jones(j,ig01(i)),j=1,3), &
                   & ((jg1(j,k,i),j=1,2),k=1,3)
            end do
  100       format(i5,'   (',i2,')',3x,a5,5x,'(',2(a2,','),a2,')',5x, &
                 &'(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )

            write(jf,*) ' '
            write(jf,*) ' index ig10(i) :'
            write(jf,120) (i,i= 1,24)
            write(jf,140) (ig10(i),i= 1,24)
            write(jf,*) ' '
            write(jf,120) (i,i=25,48)
            write(jf,140) (ig10(i),i=25,48)
  120       format('    i=',24i3)
  140       format(' ig10=',24i3)
         end if

         if(jpr >= 1) then
            no = ng1/6
            mo = mod(ng1,6)
            write(jf,*) '   '
            write(jf,*) ' matrix representation of operation '
            write(jf,*) ' (for real-space coordinate)'
            do ii=1, no
               write(jf,*) '   '
               i1 = (ii-1)*6+1
               i2 = ii*6
               write(jf,200) (i,(irotr1(1,k,i),k=1,3),i=i1,i2)
               do j=2, 3
                  write(jf,220) (  (irotr1(j,k,i),k=1,3),i=i1,i2)
               end do
            end do

            if(mo /= 0) then
               i1 = no*6+1
               i2 = no*6+mo
               write(jf,*) '   '
               write(jf,200) (i,(irotr1(1,k,i),k=1,3),i=i1,i2)
               do j=2, 3
                  write(jf,220) (  (irotr1(j,k,i),k=1,3),i=i1,i2)
               end do
            end if
            write(jf,*) '   '
            write(jf,*) ' matrix representation of operation '
            write(jf,*) ' (for reciprocal-space coordinate)'
            do  ii=1, no
               write(jf,*) '   '
               i1 = (ii-1)*6+1
               i2 = ii*6
               write(jf,200) (i,(irotk1(1,k,i),k=1,3),i=i1,i2)
               do j=2, 3
                  write(jf,220) (  (irotk1(j,k,i),k=1,3),i=i1,i2)
               end do
            end do

            if(mo /= 0) then
               i1 = no*6+1
               i2 = no*6 + mo
               write(jf,*) '   '
               write(jf,200) (i,(irotk1(1,k,i),k=1,3),i=i1,i2)
               do j=2, 3
                  write(jf,220) (  (irotk1(j,k,i),k=1,3),i=i1,i2)
               end do
            end if
         end if
  200    format(1h ,6('(',i2,')',3i2,2x))
  220    format(1h ,6(4x        ,3i2,2x))

         if(jpr >= 1) then
            write(jf,*) '   '
            write(jf,*) '--- group multiplication table ---'
            if(ng1 <= 24) then
               write(jf,320) (j,j=1,ng1)
               write(jf,340) ('---',j=1,ng1)
               do i=1, ng1
                  write(jf,300) i,(im1(i,j),j=1,ng1)
               end do
            else
               write(jf,*) ' '
               write(jf,320) (j,j=1,24)
               write(jf,340) ('---',j=1,24)
               do  i=1, ng1
                  write(jf,300) i,(im1(i,j),j=1,24)
               end do
               write(jf,*) ' '
               write(jf,320) (j,j=25,ng1)
               write(jf,340) ('---',j=25,ng1)
               do i=1,ng1
                  write(jf,300) i,(im1(i,j),j=25,ng1)
               end do
  300          format((i3,2x,24i3))
  320          format((5x,24i3))
  340          format(5x,24a3)
            end if
            write(jf,*) ' '
            write(jf,*) '--- invers elements ---'
            write(jf,400) (iv1(j),j=1,ng1)
  400       format((5x,24i3))
         end if
      end  subroutine nsgrp1
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsgrpb(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
!
!#12  input:     ng1 : order of group
!#12        tab(3,3) : transformation matix
!#12        tba(3,3) : transformation matix
!#12             ta1 : nonprimitive translation vector (A system)
!#12             ra1 : rotation matrix in real space (A system)
!#12             sa1 : rotation matrix in reciprocal space (A system)
!#12                   ra1*sa1=unit matrix
!#12  output:    tb1 : nonprimitive translation vector (B system)
!#12             rb1 : rotation matrix in real space (B system)
!#12             sb1 : rotation matrix in reciprocal space (B system)
!#12                   rb1*sb1=unit matrix
!
!#21  to get space-group matrix notation in B system
!
!#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
!     2020.12.20 ; OMP directives for parallel compiler by T. Hamada
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsgrpb(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
         use m_Const_Parameters, only :  DP
!$ use omp_lib
         implicit none
         integer , intent(in) :: ng1
         integer :: i, j, k, l, m
         real(DP), intent(in), dimension(3,3) :: tba, tab
         real(DP), intent(in), dimension(3,48) :: ta1
         real(DP), intent(in), dimension(3,3,48) :: ra1, sa1
         real(DP), intent(out), dimension(3,48) :: tb1
         real(DP), intent(out), dimension(3,3,48) :: rb1, sb1
         real(DP) :: r, s

         do i = 1, ng1
            tb1(1,i) = tba(1,1)*ta1(1,i)+tba(1,2)*ta1(2,i)+tba(1,3)*ta1(3,i)
            tb1(2,i) = tba(2,1)*ta1(1,i)+tba(2,2)*ta1(2,i)+tba(2,3)*ta1(3,i)
            tb1(3,i) = tba(3,1)*ta1(1,i)+tba(3,2)*ta1(2,i)+tba(3,3)*ta1(3,i)
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  r = 0.0d0
                  s = 0.0d0
                  do l=1, 3
                     do m=1, 3
                        r = r+tba(i,l)*ra1(l,m,k)*tab(m,j)
                        s = s+tba(i,l)*sa1(l,m,k)*tab(m,j)
                     end do
                  end do
                  rb1(i,j,k) = r
                  sb1(i,j,k) = s
               end do
            end do
         end do
      end subroutine nsgrpb

! ================================= added by K. Tagami =============== 12.0A
      subroutine nsgrpb_kt(ng1,tab,tba,ta1,ra1,sa1, tb1,rb1,sb1)
         use m_Const_Parameters, only :  DP
!   2021.12.20 OMP Directives for parallel compilers by T. Hamada
!$ use omp_lib
         implicit none
!  inout
         integer, intent(in)  :: ng1
         real(DP), intent(in), dimension(3,3) :: tab, tba
         real(DP), intent(in), dimension(3,48) :: ta1
         real(DP), intent(in), dimension(3,3,48) :: ra1, sa1
         real(DP), intent(out), dimension(3,48) :: tb1
         real(DP), intent(out), dimension(3,3,48) :: rb1, sb1
! local
         integer :: i, j, k, l, m
         real(DP) ::  r, s
! begin

         do k=1, ng1
            tb1(1,k)=tab(1,1)*ta1(1,k)+tab(1,2)*ta1(2,k)+tab(1,3)*ta1(3,k)
            tb1(2,k)=tab(2,1)*ta1(1,k)+tab(2,2)*ta1(2,k)+tab(2,3)*ta1(3,k)
            tb1(3,k)=tab(3,1)*ta1(1,k)+tab(3,2)*ta1(2,k)+tab(3,3)*ta1(3,k)
         end do

         do k=1,ng1
            do i=1,3
               do j=1,3
                  r=0.0d0
                  s=0.0d0
                  do l=1,3
                     do m=1,3
                        r = r +tab(i,l) *ra1(l,m,k) *tba(m,j)
                        s = s +tab(i,l) *sa1(l,m,k) *tba(m,j)
                     end do
                  end do
                  rb1(i,j,k) = r
                  sb1(i,j,k) = s
               end do
            end do
         end do
      end subroutine nsgrpb_kt
! ==================================================================== 12.0A

! ===================== added by K. Tagami =========================== 12.0A
      subroutine nsgrpa( ng1,tac,tca, rc1, sc1, ra1, sa1 )
         use m_Const_Parameters, only :  DP
         implicit none
! inout
         integer, intent(in) ::  ng1
         real(DP), intent(in), dimension(3,3) :: tac, tca
         real(DP), intent(in), dimension(3,3,48) :: rc1, sc1
         real(DP), intent(out), dimension(3,3,48) :: ra1, sa1
! local
         integer :: i, j, k, m1, m2
         real(DP) ::  c1, c2
! begin
! ---------------------------------------------
! Note :
!          ra1 is similar to m_CS_op_in_PUCD
!                             in  m_Crystal_Structure.F90
!
!    ----------------------------------------------
         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  c1 = 0.0d0
                  c2 = 0.0d0
                  do m1=1,  3
                     do m2=1, 3
                        c1 = c1 + tac(i,m1) *rc1(m1,m2,k) *tca(m2,j)
                        c2 = c2 + tac(i,m1) *sc1(m1,m2,k) *tca(m2,j)
                     end do
                  end do
                  ra1(i,j,k) = c1
                  sa1(i,j,k) = c2
               end do
            end do
         end do

      end subroutine nsgrpa
! ==================================================================== 12.0A
!
!#11  sub.nsjonh(irot,ir1234,jones)
!
!#12  input : irot(3,3,24)    (i) : rotation matrix in integer
!#12  output: irot(3,3,24)    (i) : rotation matrix in integer
!#12          ir1234(3,24)    (i) : jones faithful representation
!#12           jones(3,24)    (a2): jones faithful representation
!#13  noexternal
!
!#21  to get rotation matrix, and jones faithfull representation
!#21         for point-group operations of hexagonal lattice
!
!#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsjonh(irot,ir1234,jones)
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in), dimension(3,3,24) :: irot
         integer, intent(out), dimension(3,24) :: ir1234
         integer :: i, j, k, k1, k2, k3
         character(len=2), intent(out), dimension(3,24) :: jones

         do i=1, 24
            do  j=1, 2
               k1 = irot(j,1,i)
               k2 = irot(j,2,i)
               k3 = irot(j,3,i)
               if(k3 /= 0) then
                  if(printable) write(nfout,9000) i,j,(irot(j,k,i),k=1,3)
                  stop ' === stop in sub.nsjonh. (invalid r, 1)'
               end if
               if(k1 == 1) then
                  if(k2 == 1) then
                     ir1234(j,i)=4
                     jones(j,i)=' w'
!             --- for reciprocal space ---
                  else if(k2 == 0) then
                     ir1234(j,i)=1
                     jones(j,i)=' x'
                  else if(k2 == -1) then
                     ir1234(j,i)=4
                     jones(j,i)=' w'
                  else
                    if(printable) write(nfout,9000) i,j,(irot(j,k,i),k=1,3)
                    stop ' === stop in sub.nsjonh. (invalid r, 1)'
                  end if
               else if(k1 == 0) then
                  if(k2 == 1) then
                     ir1234(j,i)=2
                     jones(j,i)=' y'
                  else if(k2 == -1) then
                     ir1234(j,i)=-2
                     jones(j,i)='-y'
                  else
                     if(printable) write(nfout,9000) i,j,(irot(j,k,i),k=1,3)
                     stop ' === stop in sub.nsjonh. (invalid r, 3) ==='
                  end if
               else if(k1 == -1) then
                  if(k2 == 1) then
                     ir1234(j,i)=-4
                     jones(j,i)='-w'
                  else if(k2 == 0) then
                     ir1234(j,i)=-1
                     jones(j,i)='-x'
                  else if(k2 == -1) then
                     ir1234(j,i)=-4
                     jones(j,i)='-w'
!             --- for reciprocal space ---
                  else
                     if(printable) write(nfout,9000) i,j,(irot(j,k,i),k=1,3)
                     stop ' === stop in sub.nsjonh. (invalid r, 4) ==='
                  end if
               end if
            end do

            j=3
            if(irot(3,1,i) /= 0 .or. irot(3,2,i) /= 0) then
               if(printable) write(nfout,9000) i,j,(irot(j,k,i),k=1,3)
               stop ' === stop in sub.nsjonh. (invalid r, 5) ==='
            else if(irot(3,3,i) == 1)  then
               ir1234(j,i)=3
               jones(j,i)=' z'
            else if(irot(3,3,i) == -1) then
               ir1234(j,i)=-3
               jones(j,i)='-z'
            else
               if(printable) write(6,9000) i,j,(irot(j,k,i),k=1,3)
               stop ' === stop in sub.nsjonh. (invalid r, 6) ==='
            end if
         end do
 9000  format(1h ,'operation=',i2,'  j=',i1,'   (r(j,k),k=1,3)',3i3)
      end subroutine nsjonh
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
!
!     (nxx+1)*(nyy+1)*(nzz+1) :  no. of mesh points
!     nxx*nyy*nzz :  no. of cubes
!     nx1,ny1,nz1 :  for constructing k point
!     nd :  denominator for constructing k point
!
         use m_Const_Parameters, only :  DP!
         implicit none
         integer, intent(in) :: il, nx, ny,nz
         integer, intent(out) :: nxx, nyy, nzz, nx1, ny1, nz1, nd
         integer, save :: nxxs, nyys, nzzs, nx1s, ny1s, nz1s, nds
         integer :: id
         data id/1/

         if(il < -1 .or. il > 4) then
!        write(6,*) ' il=',il
            stop ' === stop in sub.nskma0. (il) ==='
         end if

         if(id == 1) then
            if(il == -1) then
               nxx = nx
               nyy = ny
               nzz = nz*3
            else if(il <= 1) then
               nxx = nx
               nyy = ny
               nzz = nz
            else if(il == 2) then
               nxx = nx
               nyy = ny*2
               nzz = nz*2
            else if(il == 3) then
               nxx = nx
               nyy = ny
               nzz = nz*2
            else if(il == 4) then
               nxx = nx
               nyy = ny*2
               nzz = nz
            else if(il == 5) then
               nxx = nx
               nyy = ny
               nzz = nz*2
            else if(il == 6) then
               nxx = nx
               nyy = ny
               nzz = nz*2
            end if
         else if(id == 2) then
            if(il == -1) then
               nxx = nx*3
               nyy = ny*3
               nzz = nz*3
            else if(il <=1) then
               nxx = nx
               nyy = ny
               nzz = nz
            else if(il <= 3) then
               nxx = nx*2
               nyy = ny*2
               nzz = nz*2
            else if(il == 4) then
               nxx = nx*2
               nyy = ny*2
               nzz = nz
            else if(il == 5) then
               nxx = nx
               nyy = ny*2
               nzz = nz*2
            else if(il == 6) then
               nxx = nx*2
               nyy = ny
               nzz = nz*2
            end if
         else
!        write(6,*) ' id=',id
            stop ' === stop in sub.nskma0. (id) ==='
         end if

         if(nx == 0) then
            nx1 = 1
         else
            nx1 = nx
         end if
         if(ny == 0) then
            ny1 = 1
         else
            ny1 = ny
         end if
         if(nz == 0) then
            nz1 = 1
         else
            nz1 = nz
         end if
         nd = nx1*ny1*nz1
      end subroutine nskma0
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp, np,k,p)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : nskp00_called
         use m_Timing, only : tstatc0_begin, tstatc0_end
! 2021.12.20: modifued by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: nxx, nyy,nzz, nx1, ny1, nz1, nd, lmnp
         integer, intent(in), dimension(4,lmnp) :: k
         integer, intent(in) :: np
         integer :: i
         real(DP) :: rnd_1
         real(DP), dimension(3,lmnp), intent(out) :: p
         real(DP), allocatable, dimension(:,:),save :: ps
         integer :: id_sname = -1

         if(nskp00_called) then
!$OMP simd
            do i = 1, np
                 p(1:3,i) = ps(1:3,i)
            end do
!$OMP end simd
            return
         end if
        call tstatc0_begin('nskp00 ',id_sname)
        rnd_1 = 1.0/real(nd,kind=DP)

!$OMP simd
        do i = 1, np
           p(1:3,i) = real(k(1:3,i),kind=DP)*rnd_1
           enddo
!$OMP end simd
        allocate(ps(3,lmnp))
!$OMP simd
        do i = 1, np
           ps(1:3,i) = p(1:3,i)
        end do
!$OMP end simd
        nskp00_called = .true.
        call tstatc0_end(id_sname)
      end subroutine nskp00
!==*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!#11 &           pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!#11 &           iu21,iv21)
!!c
!#12  input:
!!#12
!#12
!#12
!#12  output:
!#12
!#12
!#13  external:  nspbge
!
!#21  to get coordination systems (C,A,B)
!
!#31  1990.11.20.:  n. hamada, a. yanase and k. terakura
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,&
     &                  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
     &                  iu21,iv21, itrs )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
! 2021.12.20 : OMP directives and modifucations by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jpr, ng1, np0, lmnp0, lmnp1, lmnp2, itrs
         integer, intent(in), dimension(3,3,ng1) :: lsb1
         integer, intent(in), dimension(ng1)  :: iv1
         integer, intent(out) :: np1, np2
         integer, intent(out), dimension(lmnp0) :: ip10, ip20
         integer, intent(out), dimension(lmnp1) :: ip01, ip21
         integer, intent(inout), dimension(lmnp2) :: ip02, ip12
         integer, intent(out), dimension(lmnp1) :: iu21, iv21
         integer, allocatable, dimension(:,:) :: pb05i
         integer, dimension(3,3,ng1) :: lsb1_1
         integer, dimension(3) :: pi
         integer :: i, j, jj, k, kk, m, ind, ip0, ip2, id
         real(DP), intent(in), dimension(3,3) :: tab
         real(DP), intent(in), dimension(3,lmnp0) :: pa0
         real(DP), intent(out), dimension(3,lmnp0) :: pb0
         real(DP),dimension(3) :: pb0_0
         real(DP) ,dimension(3,3) :: rlsb1
         real(DP), dimension(3) :: p

        do i = 1, np0
           pb0(1,i)=pa0(1,i)*tab(1,1)+pa0(2,i)*tab(2,1)+pa0(3,i)*tab(3,1)
           pb0(2,i)=pa0(1,i)*tab(1,2)+pa0(2,i)*tab(2,2)+pa0(3,i)*tab(3,2)
           pb0(3,i)=pa0(1,i)*tab(1,3)+pa0(2,i)*tab(2,3)+pa0(3,i)*tab(3,3)
        end do

        allocate(pb05i(3,lmnp0))

        pb05i(1,1:lmnp0) = nint(100000*pb0(1,1:lmnp0))
        pb05i(2,1:lmnp0) = nint(100000*pb0(2,1:lmnp0))
        pb05i(3,1:lmnp0) = nint(100000*pb0(3,1:lmnp0))

        np1 = 1
        ip10(1)= 1
        ip01(1) = 1
        do 1  i=2, np0
           do  j=1, np1
              jj = j
! ---------- in-line expansion of sub. nspbgei start --------
              id = abs(pb05i(1,i)-pb05i(1,ip01(j)))
              if(mod(id,100000)/=0) cycle
              id = abs(pb05i(2,i)-pb05i(2,ip01(j)))
              if(mod(id,100000)/=0) cycle
              id = abs(pb05i(3,i)-pb05i(3,ip01(j)))
              if(mod(id,100000)/=0) cycle
! ---------- in-line expansion of sub. nspbgei end ----------
              ip10(i) = jj
              goto 1
           end do
           np1 = np1+1
           if(np1 > lmnp1) then
              if(printable) write(nfout,*) ' np1=',np1,' > lmnp1=',lmnp1
              stop '=nskpb0(np1)='
           end if
           ip10(i) = np1
           ip01(np1) = i
     1  end do
        np2 = 1
        ip21(1) = 1
        ip12(1) = 1
        ip02(1) = 1
        iu21(1) = 1
        iv21(1) = 1
        do 2 i=2,np1
           pb0_0(1:3) =  pb0(1:3,ip01(i))
           do k=1,ng1
              rlsb1(1:3,1:3) = real(lsb1(1:3,1:3,k),kind=DP)
              p(1:3) = pb0_0(1)*rlsb1(1,1:3) &
                &     +pb0_0(2)*rlsb1(2,1:3) &
                &     +pb0_0(3)*rlsb1(3,1:3)
              do m=0, itrs
                 kk = k
                 pi(1:3) = nint(100000*(1-2*m)*p(1:3))
                 do j=1, np2
                    jj = j
!---------- in-line expansoin of sub. nsobgei start ----------
                    id = abs(pi(1)-pb05i(1,ip02(j)))
                    if(mod(id,100000)/=0) cycle
                    id = abs(pi(2)-pb05i(2,ip02(j)))
                    if(mod(id,100000)/=0) cycle
                    id = abs(pi(3)-pb05i(3,ip02(j)))
                    if(mod(id,100000)/=0) cycle
!---------- in -line expansion of sub. nspbgei end ----------
                    ip21(i) = jj
                    iu21(i) = kk
                    iv21(i) = iv1(kk)
                    goto 2
                end do
              end do
           end do
           np2 = np2+1
           if(np2 > lmnp2) then
              if(printable) write(nfout,*) ' np2=',np2,' > lmnp2=',lmnp2
              stop '=nskpb0(np2)='
           end if
           ip21(i) = np2
           ip12(np2) = i
           ip02(np2) = ip01(i)
           iu21(i) = 1
           iv21(i) = 1
     2  end do

        ip20(1:np0)=ip21(ip10(1:np0))

        deallocate(pb05i)

        if(jpr >= .3) then
           if(printable) then
              write(nfout,*) ' === k points === '
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
           end if
        end if
 120    format((1h ,2('     (',3f9.5,' ) ')))

        if(jpr >= 2) then
           if(printable) then
              write(nfout,*) '=== k-point code === '
              write(nfout,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
              write(nfout,*) &
             & '----------------------------------------------------'
              write(nfout,180)
              write(nfout,200) (ip0,ip10(ip0),iu21(ip10(ip0)),&
              & ip21(ip10(ip0)),   ip0=1,np0)
              write(nfout,*) &
             & '----------------------------------------------------'
           end if
        end if
 140    format((1h ,2(i4,' (',3f9.5,' ) ')))
 180    format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
 200    format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
      end subroutine nskpb0

      subroutine nskpb0_s(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1, &
     &                  lmnp2,pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,&
     &                  iu21,iv21, itrs)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : nskpb0_s_called
! 2021.12.20: OpenMP directives and modifications by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jpr, ng1, np0, lmnp0, lmnp1, lmnp2,itrs
         integer, intent(in), dimension(3,3,ng1) :: lsb1
         integer, intent(in), dimension(ng1) :: iv1
         integer, intent(out) :: np1, np2
         integer, intent(out), dimension(lmnp0) :: ip10, ip20
         integer, intent(out), dimension(lmnp1) :: ip01, ip21
         integer, intent(out), dimension(lmnp2) :: ip02, ip12
         integer, intent(out), dimension(lmnp1) :: iu21, iv21
         integer, allocatable, dimension(:,:) :: pb05i
         real(DP), allocatable, dimension(:,:,:) :: pl
         real(DP),dimension(3) :: p
         integer, dimension(3) :: pi
         integer, dimension(lmnp1,3,ng1) :: pil
         integer, save :: np1s,np2s
         integer i, j, k, m, id, ip0, ip2
         integer, allocatable, dimension(:), save :: ip10s,ip20s
         integer, allocatable, dimension(:), save :: ip01s,ip21s
         integer, allocatable, dimension(:), save :: ip02s,ip12s
         integer, allocatable, dimension(:), save :: iu21s,iv21s
         real(DP), intent(in), dimension(3,3) :: tab
         real(DP), intent(in), dimension(3,lmnp0) :: pa0
         real(DP), intent(out), dimension(3,lmnp0) :: pb0
         real(DP), allocatable, dimension(:,:) :: pb0_1
         real(DP) , allocatable, dimension(:,:), save :: pb0s
         integer :: ifact

         if(nskpb0_s_called) then
           np1 = np1s
           np2 = np2s
           pb0 = pb0s
           ip10 = ip10s
           ip20 = ip20s
           ip01 = ip01s
           ip21 = ip21s
           ip02 = ip02s
           ip12 = ip12s
           iu21 = iu21s
           iv21 = iv21s
           return
        endif
        do i= 1, np0
           pb0(1,i)=pa0(1,i)*tab(1,1)+pa0(2,i)*tab(2,1)+pa0(3,i)*tab(3,1)
           pb0(2,i)=pa0(1,i)*tab(1,2)+pa0(2,i)*tab(2,2)+pa0(3,i)*tab(3,2)
           pb0(3,i)=pa0(1,i)*tab(1,3)+pa0(2,i)*tab(2,3)+pa0(3,i)*tab(3,3)
        end do

        allocate(pb0_1(3,lmnp0))
        allocate(pb05i(3,lmnp0))

        pb0_1 = 100000*pb0

        pb05i = nint(pb0_1)

        np1 = 1
        ip10(1) = 1
        ip01(1) = 1

        do 1  i=2, np0
           do  j=1, np1
! ---------- in line expansion of modified sub. nspbgei start ----------
              id = abs(pb05i(1,i)-pb05i(1,ip01(j)))
              if(mod(id,100000)/=0) cycle
              id = abs(pb05i(2,i)-pb05i(2,ip01(j)))
              if(mod(id,100000)/=0) cycle
              id = abs(pb05i(3,i)-pb05i(3,ip01(j)))
              if(mod(id,100000)/=0) cycle
! ---------- in-line expasion of modified sub. nsbpgei end ----------
              ip10(i) = j
              goto 1
           end do
           np1 = np1+1
           ip10(i)=np1
           ip01(np1)=i
     1  end do

        np2 = 1
        ip21(1) = 1
        ip12(1) = 1
        ip02(1) = 1
        iu21(1) = 1
        iv21(1) = 1

        allocate(pl(np1,3,ng1))

        do i = 2, np1
           do j = 1, 3
              do k = 1, ng1
                  pl(i,j,k) = &
               &  pb0_1(1,ip01(i))*lsb1(1,j,k) &
               & +pb0_1(2,ip01(i))*lsb1(2,j,k) &
               & +pb0_1(3,ip01(i))*lsb1(3,j,k)
              end do
           end do
        end do

        deallocate(pb0_1)

        do 2 i = 2, np1
           do k = 1, ng1
              p(1:3) = pl(i,1:3,k)
              do m  = 0, itrs
                 pi=  nint((1-2*m)*p)
                 do  j=1, np2
! ---------- in-line expansion of modified sub. nspbgei start ----------
                    id = iabs(pi(1)-pb05i(1,ip02(j)))
                    if(mod(id,100000)/=0) cycle
                    id = iabs(pi(2)-pb05i(2,ip02(j)))
                    if(mod(id,100000)/=0) cycle
                    id = iabs(pi(3)-pb05i(3,ip02(j)))
                    if(mod(id,100000)/=0) cycle
! ---------- in-line expansion of modified sub npsbgei end ----------
                    ip21(i) = j
                    iu21(i) = k
                    iv21(i) = iv1(k)
                    goto 2
                 end do
              end do
           end do
           np2 = np2 + 1
           ip21(i) = np2
           ip12(np2) = i
           ip02(np2) = ip01(i)
           iu21(i) = 1
           iv21(i) = 1
     2  end do

        deallocate(pb05i)
        deallocate(pl)

        do i = 1,np0
           ip20(i)=ip21(ip10(i))
        end do

        if(jpr >= 3) then
           if(printable) then
              write(nfout,*) ' === k points === '
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
           end if
        end if
 120    format((1h ,2('     (',3f9.5,' ) ')))

        if(jpr >= 2) then
           if(printable ) then
              write(nfout,*) '=== k-point code === '
              write(nfout,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
              write(nfout,180)
              write(nfout,200) (ip0,ip10(ip0),iu21(ip10(ip0)),&
              & ip21(ip10(ip0)),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
           end if
        end if
 140    format((1h ,2(i4,' (',3f9.5,' ) ')))
 180    format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
 200    format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
        np1s = np1
        np2s = np2
        allocate(pb0s(3,lmnp0))
        allocate(ip10s(lmnp0))
        allocate(ip20s(lmnp0))
        allocate(ip01s(lmnp1))
        allocate(ip21s(lmnp1))
        allocate(ip02s(lmnp2))
        allocate(ip12s(lmnp2))
        allocate(iu21s(lmnp1))
        allocate(iv21s(lmnp1))
        pb0s= pb0
        ip10s=ip10
        ip20s=ip20
        ip01s=ip01
        ip21s=ip21
        ip02s=ip02
        ip12s=ip12
        iu21s=iu21
        iv21s=iv21
        nskpb0_s_called = .true.
      end subroutine nskpb0_s

      subroutine nskpb0_s_1(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1, &
        &    lmnp2,pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,&
        &    iu21,iv21, itrs)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Crystal_Structure, only     : inversion_symmetry
         use m_Spg_plus_Tetra_Common_Blocks, only : nskpb0_s_called
         use m_Timing, only : tstatc0_begin, tstatc0_end
! 2021.12.20: OpenMP directives and modifications by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jpr, ng1, np0, lmnp0, lmnp1, lmnp2,itrs
         integer, intent(in), dimension(3,3,ng1) :: lsb1
         integer, intent(in), dimension(ng1) :: iv1
         integer, intent(out) :: np1, np2
         integer, intent(out), dimension(lmnp0) :: ip10, ip20
         integer, intent(out), dimension(lmnp1) :: ip01, ip21
         integer, intent(out), dimension(lmnp2) :: ip02, ip12
         integer, intent(out), dimension(lmnp1) :: iu21, iv21
         integer, allocatable, dimension(:,:) :: pb05i
         integer, allocatable, dimension(:,:,:) :: pil
         integer, save :: np1s,np2s
         integer i, j, k, id, m, ip0, ip2
         integer, allocatable, dimension(:), save :: ip10s,ip20s
         integer, allocatable, dimension(:), save :: ip01s,ip21s
         integer, allocatable, dimension(:), save :: ip02s,ip12s
         integer, allocatable, dimension(:), save :: iu21s,iv21s
         real(DP), intent(in), dimension(3,3) :: tab
         real(DP), intent(in), dimension(3,lmnp0) :: pa0
         real(DP), intent(out), dimension(3,lmnp0) :: pb0
         real(DP), allocatable, dimension(:,:) :: pb0_1
         real(DP) , allocatable, dimension(:,:), save :: pb0s
         real(DP) :: x, y, z
         integer :: id_sname = -1

         if(nskpb0_s_called) then
           np1 = np1s
           np2 = np2s
!$OMP simd
           do i = 1, np0
              pb0(1:3,i) = pb0s(1:3,i)
           end do
!$OMP end simd
           ip10 = ip10s
           ip20 = ip20s
           ip01 = ip01s
           ip21 = ip21s
           ip02 = ip02s
           ip12 = ip12s
           iu21 = iu21s
           iv21 = iv21s
           return
        endif
        call tstatc0_begin('nskpb0_s_1 ',id_sname)
        allocate(pb0_1(3,np0))
        allocate(pb05i(3,np0))

        call set_pb0

        call set_pb0_1

        call set_pb05i

        np1 = 1
        ip10(1) = 1
        ip01(1) = 1

!$OMP ordered
           do 1  i=2, np0
              do  j=1, np1
! ---------- in line expansion of modified sub. nspbgei start ----------
                 id = pb05i(1,i)-pb05i(1,ip01(j))
                 id = mod(id,100000)
                 if(id /=0) cycle
                 id = pb05i(2,i)-pb05i(2,ip01(j))
                 id = mod(id,100000)
                 if(id/=0) cycle
                 id = pb05i(3,i)-pb05i(3,ip01(j))
                 id = mod(id,100000)
                 if(id/=0) cycle
! ---------- in-line expasion of modified sub. nsbpgei end ----------
                 ip10(i) = j
                 goto 1
              end do
              np1 = np1+1
              ip10(i)=np1
              ip01(np1)=i
     1     end do

!$OMP end ordered

        np2 = 1
        ip21(1) = 1
        ip12(1) = 1
        ip02(1) = 1
        iu21(1) = 1
        iv21(1) = 1

!$OMP ordered
       if(ng1>1) then
          allocate(pil(2:np1,3,ng1))
           call set_pil
           do 2 i = 2, np1
              do k = 1, ng1
                 do  j=1, np2
! ---------- in-line expansion of modified sub. nspbgei start ----------
                    id = pil(i,1,k)-pb05i(1,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
                    id = pil(i,2,k)-pb05i(2,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
                    id = pil(i,3,k)-pb05i(3,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
! ---------- in-line expansion of modified sub npsbgei end ----------
                    ip21(i) = j
                    iu21(i) = k
                    iv21(i) = iv1(k)
                    goto 2
                 end do
              end do
              np2 = np2 + 1
              ip21(i) = np2
              ip12(np2) = i
              ip02(np2) = ip01(i)
              iu21(i) = 1
              iv21(i) = 1
     2     end do
           deallocate(pil)
      else
           np2 = np1
           do i = 1, np2
              ip02(i) =  ip01(i)
              ip21(i) =  i
              ip12(i) =  i
           end do
           iu21 = 1
           iv21 = iv1(1)
        end if
!$OMP end ordered

        deallocate(pb0_1)
        deallocate(pb05i)

        call set_ip20

        if(jpr >= 3) then
           if(printable) then
              write(nfout,*) ' === k points === '
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
           end if
        end if
 120    format((1h ,2('     (',3f9.5,' ) ')))

        if(jpr >= 2) then
           if(printable ) then
              write(nfout,*) '=== k-point code === '
              write(nfout,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
              write(nfout,180)
              write(nfout,200) (ip0,ip10(ip0),iu21(ip10(ip0)),&
              & ip21(ip10(ip0)),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
           end if
        end if
 140    format((1h ,2(i4,' (',3f9.5,' ) ')))
 180    format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
 200    format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
        np1s = np1
        np2s = np2
        allocate(pb0s(3,lmnp0))
        allocate(ip10s(lmnp0))
        allocate(ip20s(lmnp0))
        allocate(ip01s(lmnp1))
        allocate(ip21s(lmnp1))
        allocate(ip02s(lmnp2))
        allocate(ip12s(lmnp2))
        allocate(iu21s(lmnp1))
        allocate(iv21s(lmnp1))
!$OMP simd
        do i = 1, np0
           pb0s(1:3,i)= pb0(1:3,i)
        end do
!$OMP end simd
        ip10s=ip10
        ip20s=ip20
        ip01s=ip01
        ip21s=ip21
        ip02s=ip02
        ip12s=ip12
        iu21s=iu21
        iv21s=iv21
        nskpb0_s_called = .true.
        call tstatc0_end(id_sname)
        contains
         subroutine set_pb0
!$OMP simd
           do i  = 1, np0
              pb0(1:3,i) =  pa0(1,i)*tab(1,1:3)
              pb0(1:3,i) =  pb0(1:3,i)+pa0(2,i)*tab(2,1:3)
              pb0(1:3,i) =  pb0(1:3,i)+pa0(3,i)*tab(3,1:3)
           end do
!$OMP end simd
         end subroutine set_pb0

         subroutine set_pb0_1
!$OMP simd
           do i = 1, np0
              pb0_1(1:3,i) = 1.0d5*pb0(1:3,i)
           end do
!$OMP end simd
         end subroutine set_pb0_1

         subroutine set_pb05i
!$OMP simd
           do i = 1, np0
              pb05i(1:3,i)   = nint(pb0_1(1:3,i))
           end do
!$OMP end simd
         end subroutine set_pb05i

         subroutine set_pil
!$OMP simd
          do i = 2, np1
             do k = 1, ng1
                 pil(i,1:3,k) = &
             &   nint(pb0_1(1,ip01(i))*lsb1(1,1:3,k) &
             &       +pb0_1(2,ip01(i))*lsb1(2,1:3,k) &
             &       +pb0_1(3,ip01(i))*lsb1(3,1:3,k))
             end do
          enddo
!$OMP end simd
         end subroutine set_pil

         subroutine set_ip20
          if(ng1>1) then
!$OMP simd
             do i = 1, np0
                ip20(i) = ip21(ip10(i))
             end do
!$OMP end simd
         else
!$OMP simd
            do i=1, np0
               ip20(i) = ip10(i)
            end do
!$OMP end simd
        end if
      end subroutine set_ip20
      end subroutine nskpb0_s_1

      subroutine nskpb0_s_2(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1, &
       &    lmnp2,pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,&
       &    iu21,iv21, itrs)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : nskpb0_s_called
         use m_Timing, only : tstatc0_begin, tstatc0_end
! 2021.12.20: OpenMP directives and modifications by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jpr, ng1, np0, lmnp0, lmnp1, lmnp2,itrs
         integer, intent(in), dimension(3,3,ng1) :: lsb1
         integer, intent(in), dimension(ng1) :: iv1
         integer, intent(out) :: np1, np2
         integer, intent(out), dimension(lmnp0) :: ip10, ip20
         integer, intent(out), dimension(lmnp1) :: ip01, ip21
         integer, intent(out), dimension(lmnp2) :: ip02, ip12
         integer, intent(out), dimension(lmnp1) :: iu21, iv21
         integer, allocatable, dimension(:,:) :: pb05i
         integer, allocatable, dimension(:,:,:,:) :: pil
         integer, save :: np1s,np2s
         integer i, j, k, m, id, ip0, ip2
         integer, allocatable, dimension(:), save :: ip10s,ip20s
         integer, allocatable, dimension(:), save :: ip01s,ip21s
         integer, allocatable, dimension(:), save :: ip02s,ip12s
         integer, allocatable, dimension(:), save :: iu21s,iv21s
         real(DP), intent(in), dimension(3,3) :: tab
         real(DP), intent(in), dimension(3,lmnp0) :: pa0
         real(DP), intent(out), dimension(3,lmnp0) :: pb0
         real(DP), allocatable, dimension(:,:) :: pb0_1
         real(DP) , allocatable, dimension(:,:), save :: pb0s
         integer :: id_sname = -1

         if(nskpb0_s_called) then
           np1 = np1s
           np2 = np2s
           pb0 = pb0s
           ip10 = ip10s
           ip20 = ip20s
           ip01 = ip01s
           ip21 = ip21s
           ip02 = ip02s
           ip12 = ip12s
           iu21 = iu21s
           iv21 = iv21s
           return
        endif
        call tstatc0_begin('nskpb0_s_2 ',id_sname)

        allocate(pb05i(3,np0))
        allocate(pb0_1(3,np0))

        call set_pb0

        call set_pb0_1

        call set_pb05i

        np1 = 1
        ip10(1) = 1
        ip01(1) = 1

        do 1  i=2, np0
           do  j=1, np1
! ---------- in line expansion of modified sub. nspbgei start ----------
              id = pb05i(1,i)-pb05i(1,ip01(j))
              id = mod(id,100000)
              if(id/=0) cycle
              id = pb05i(2,i)-pb05i(2,ip01(j))
              id = mod(id,100000)
              if(id/=0) cycle
              id = pb05i(3,i)-pb05i(3,ip01(j))
              id = mod(id,100000)
              if(id/=0) cycle
              ip10(i) = j
              goto 1
           end do
           np1 = np1+1
           ip10(i)=np1
           ip01(np1)=i
     1  end do

        np2 = 1
        ip21(1) = 1
        ip12(1) = 1
        ip02(1) = 1
        iu21(1) = 1
        iv21(1) = 1

        allocate(pil(2:np1,3,ng1,0:itrs))

        call set_pil

        deallocate(pb0_1)

        do 2 i = 2, np1
           do k = 1, ng1
              do m  = 0, itrs
                  do  j=1, np2
! ---------- in-line expansion of modified sub. nspbgei start ----------
                    id = pil(i,1,k,m)-pb05i(1,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
                    id = pil(i,2,k,m)-pb05i(2,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
                    id = pil(i,3,k,m)-pb05i(3,ip02(j))
                    id = mod(id,100000)
                    if(id/=0) cycle
! ---------- in-line expansion of modified sub npsbgei end ----------
                    ip21(i) = j
                    iu21(i) = k
                    iv21(i) = iv1(k)
                    goto 2
                 end do
              end do
           end do
           np2 = np2 + 1
           ip21(i) = np2
           ip12(np2) = i
           ip02(np2) = ip01(i)
           iu21(i) = 1
           iv21(i) = 1
     2  end do

        deallocate(pil)
        deallocate(pb05i)

        call set_ip20

        if(jpr >= 3) then
           if(printable) then
              write(nfout,*) ' === k points === '
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
              write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
              write(nfout,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
           end if
        end if
 120    format((1h ,2('     (',3f9.5,' ) ')))

        if(jpr >= 2) then
           if(printable ) then
              write(nfout,*) '=== k-point code === '
              write(nfout,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
              write(nfout,180)
              write(nfout,200) (ip0,ip10(ip0),iu21(ip10(ip0)),&
              & ip21(ip10(ip0)),ip0=1,np0)
              write(nfout,*) &
              & '----------------------------------------------------'
           end if
        end if
 140    format((1h ,2(i4,' (',3f9.5,' ) ')))
 180    format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
 200    format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
        np1s = np1
        np2s = np2
        allocate(pb0s(3,lmnp0))
        allocate(ip10s(lmnp0))
        allocate(ip20s(lmnp0))
        allocate(ip01s(lmnp1))
        allocate(ip21s(lmnp1))
        allocate(ip02s(lmnp2))
        allocate(ip12s(lmnp2))
        allocate(iu21s(lmnp1))
        allocate(iv21s(lmnp1))
        pb0s= pb0
        ip10s=ip10
        ip20s=ip20
        ip01s=ip01
        ip21s=ip21
        ip02s=ip02
        ip12s=ip12
       iu21s=iu21
        iv21s=iv21
        nskpb0_s_called = .true.
        call tstatc0_end(id_sname)
        contains
         subroutine set_pb0
!$OMP simd
           do i  = 1, np0
              pb0(1:3,i) = pa0(1,i)*tab(1,1:3)
              pb0(1:3,i) = pb0(1:3,i)+pa0(2,i)*tab(2,1:3)
              pb0(1:3,i) = pb0(1:3,i)+pa0(3,i)*tab(3,1:3)
           end do
!$OMP end simd
         end subroutine set_pb0

         subroutine set_pb0_1
!$OMP simd
           do i = 1, np0
              pb0_1(1:3,i) = 1.0d5*pb0(1:3,i)
           end do
!$OMP end simd
          end subroutine set_pb0_1

         subroutine set_pb05i
!$OMP simd
           do i = 1, np0
              pb05i(1:3,i)   = nint(pb0_1(1:3,i))
           end do
!$OMP end simd
         end subroutine set_pb05i

         subroutine set_pil
!$OMP simd
           do i = 2, np1
              do k = 1, ng1
                    pil(i,1:3,k,0) = nint(   &
                 & (pb0_1(1,ip01(i))*lsb1(1,1:3,k) &
                 & +pb0_1(2,ip01(i))*lsb1(2,1:3,k) &
                 & +pb0_1(3,ip01(i))*lsb1(3,1:3,k)))
                    pil(i,1:3,k,1) = nint(-1*   &
                 & (pb0_1(1,ip01(i))*lsb1(1,1:3,k) &
                 & +pb0_1(2,ip01(i))*lsb1(2,1:3,k) &
                 & +pb0_1(3,ip01(i))*lsb1(3,1:3,k)))
              end do
           end do
!$OMP end simd
         end subroutine set_pil

         subroutine set_ip20
!$OMP simd
           do i = 1, np0
              ip20(i) = ip21(ip10(i))
           end do
!$OMP end simd
         end subroutine set_ip20
      end subroutine nskpb0_s_2

! ================================ added by K. Tagami ================= 12.0A
      subroutine nskpa0_kt(jpr,tab,ng1,lsa1,iv1,np0,pa0,lmnp0, &
     &                  lmnp1,lmnp2, &
     &                  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
     &                  iu21,iv21, itrs )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
!2021.12.20 : OpenMP directives and modifictions by T. Hmmada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jpr, ng1, np0, lmnp0, lmnp1, lmnp2,itrs
         integer, intent(in), dimension(3,3,lmnp0) :: lsa1
         integer, intent(in), dimension(ng1) :: iv1
         integer, intent(out) :: np1, np2
         integer, intent(out), dimension(lmnp0) :: ip10, ip20
         integer, intent(out), dimension(lmnp1) :: ip01, ip21
         integer, intent(out), dimension(lmnp2) :: ip02, ip12
         integer, intent(out), dimension(lmnp1) :: iu21, iv21
         integer :: i, j, jj, k, kk, m, ind, ip0, ip2
         real(DP), intent(in), dimension(3,3) :: tab
         real(DP), intent(inout), dimension(3,lmnp0) :: pa0
         real(DP), intent(out), dimension(3,lmnp0) :: pb0
         real(DP) :: eps, d, factor
         real(DP), dimension(3,3) :: rlsa1
         real(DP), dimension(3) :: pa0_0, p

         eps = 1.0d-5
         pb0(1:3,1:lmnp0) = pa0(1:3,1:lmnp0)

         do i = 1, np0
            pa0(1,i)=pb0(1,i)*tab(1,1)+pb0(2,i)*tab(2,1)+pb0(3,i)*tab(3,1)
            pa0(2,i)=pb0(1,i)*tab(1,2)+pb0(2,i)*tab(2,2)+pb0(3,i)*tab(3,2)
            pa0(3,i)=pb0(1,i)*tab(1,3)+pb0(2,i)*tab(2,3)+pb0(3,i)*tab(3,3)
         end do

         np1 = 1
         ip10(1) = 1
         ip01(1) = 1

         do 1 i=2, np0
            do j=1, np1
               jj = j
! ---------- in-line expansion of sub. nspbge start ----------
               d = dabs(pa0(1,i)-pa0(1,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
               d = dabs(pa0(2,i)-pa0(2,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
               d = dabs(pa0(3,i)-pa0(3,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
               ind = 0
! ---------- in-line expansion of sub. nspbge end ----------
               ip10(i) = jj
               goto 1
            end do
            np1 = np1+1
            if(np1 > lmnp1) then
               if(printable) write(6,*) ' np1=',np1,' > lmnp1=',lmnp1
               stop '=nskpb0(np1)='
            end if
            ip10(i) = np1
            ip01(np1) = i
     1   end do
         np2 = 1
         ip21(1) = 1
         ip12(1) = 1
         ip02(1) = 1
         iu21(1) = 1
         iv21(1) = 1
         do 2 i=2,np1
            pa0_0(1:3) = pa0(1:3,ip01(i))
            do k=1,ng1
               p(1:3) = pa0_0(1)*real(lsa1(1,1:3,k),kind=DP) &
                     & +pa0_0(2)*real(lsa1(2,1:3,k),kind=DP) &
                     & +pa0_0(3)*real(lsa1(3,1:3,k),kind=DP)
               do m=0, itrs
                  kk =k
                  factor = real((1-2*m),kind=DP)
                  p(1:3) = p(1:3)*factor
                  do j=1, np2
                     jj = j
! ---------- in-line expansion of sub. nspbge start ----------
                     d = dabs(p(1)-pa0(1,ip02(j))) + 0.5d-5
                     if(dmod(d,1.0d0) <= 1.0d-5) cycle
                     d = dabs(p(2)-pa0(2,ip02(j))) + 0.5d-5
                     if(dmod(d,1.0d0) <= 1.0d-5) cycle
                     d = dabs(p(3)-pa0(3,ip02(j))) + 0.5d-5
                     if(dmod(d,1.0d0) <= 1.0d-5) cycle
! ---------- in0line expansion of sub. nspbge end ---------
                     ip21(i) = jj
                     iu21(i) = kk
                     iv21(i) = iv1(kk)
                     goto 2
                  end do
               end do
            end do
            np2 = np2+1
            if(np2 > lmnp2) then
               if(printable)  write(nfout,*) ' np2=',np2,' > lmnp2=',lmnp2
               stop '=nskpb0(np2)='
            end if
            ip21(i) = np2
            ip12(np2) = i
            ip02(np2) = ip01(i)
            iu21(i) = 1
            iv21(i) = 1
      2  end do

         ip20(1:np0)=ip21(ip10(1:np0))

         if(jpr >= 3) then
            if(printable) then
               write(nfout,*) ' === k points === '
               write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
               write(nfout,120) ((pa0(i,ip02(ip2)),i=1,3),ip2=1,np2)
               write(nfout,*) ' np0=',np0,'   np1=',np1,'   np2=',np2
               write(nfout,120) ((pb0(i,ip02(ip2)),i=1,3),ip2=1,np2)
            end if
         end if
 120     format((1h ,2('     (',3f9.5,' ) ')))

         if(jpr >= 2) then
            if(printable) then
               write(nfout,*) '=== k-point code === '
               write(nfout,140) (ip0,(pa0(i,ip0),i=1,3),ip0=1,np0)
               write(nfout,*) &
              & '------------------------------------------------'
               write(nfout,180)
               write(nfout,200) (ip0,ip10(ip0),iu21(ip10(ip0)),&
              &  ip21(ip10(ip0)), ip0=1,np0)
               write(6,*) &
              & '------------------------------------------------'
            end if
         end if
 140     format((1h ,2(i4,' (',3f9.5,' ) ')))
 180     format(1h ,2(' ip0',4x,' ip1',' (iu21)  ip2',5x))
 200     format((1h ,2(i4,' -->',i4,' -(',i2,')->',i4,5x)))
      end subroutine nskpa0_kt
! ==================================================================== 12.0A
!
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!c
!#11        subroutine nskpbm(np0,pa0,np1,ip10,ip01)
!#12  input:
!#12
!#12
!#12
!#12  output:
!#12
!#12
!#13  external:  nspbge
!
!
!#41  2002.12.10 This subroutine is modified from nskpb0 for tetrahedron method in the program 'phase'.
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nskpbm(np0,lmnp0,lmnp1,pa0,np1,ip10,ip01)
         use m_Const_Parameters, only :  DP
! 2021.12.20: Open MP directives and modifications by T. Hamada
         implicit none
         integer, intent(in) :: np0, lmnp0, lmnp1
         integer, intent(out) :: np1
         integer, intent(out), dimension(lmnp0) :: ip10
         integer, intent(out), dimension(lmnp1) :: ip01
         integer :: i, j, k, jj, ip01_j
         real(DP), intent(in), dimension(3,lmnp0) :: pa0
         real(DP) :: d

         np1 = 1
         ip10(1) = 1
         ip01(1) = 1
         do 1 i=2, np0
            do  j=1, np1
               jj = j
! ---------- in-line expansion of sub. nspbge start ----------
               d = dabs(pa0(1,i)-pa0(1,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
               d = dabs(pa0(2,i)-pa0(2,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
               d = dabs(pa0(3,i)-pa0(3,ip01(j))) + 0.5d-5
               if(dmod(d,1.0d0) <= 1.0d-5) cycle
! ---------- in-line expansion of sub. nspbge end ----------
               ip10(i) = jj
               goto 1
            end do
            np1 = np1+1
            ip10(i) = np1
            ip01(np1) = i
      1  end do
      end  subroutine nskpbm
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
!#11 &                                 grc,gkc,gra,gka,grb,gkb)
!
!#12  input:     jpr : print control
!#12              il : lattice type
!#12           a,b,c : lattice parameter (length)
!#12        ca,cb,cc : lattice parameter (cos(angle))
!#12  output:  tca(3,3),tac(3,3),tab(3,3),tba(3,3),tcb(3,3),tbc(3,3)
!#12           grc(3,3),gkc(3,3),gra(3,3),gka(3,3),grb(3,3),gkb(3,3)
!#13  external:  nslatc,nslata,nslatb
!
!#21  to get coordination systems (C,A,B)
!
!#31  1990.01.20.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
! ================================= modified by K. Tagami =============== 12.0A
!      subroutine nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
!     &                  grc,gkc,gra,gka,grb,gkb)
      subroutine nslat3(jpr,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc, &
     &                  grc,gkc,gra,gka,grb,gkb, &
     &                  use_altv_rltv, altv, rltv )
! ======================================================================== 12.0A
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: jpr, il
         integer :: i, j
         real(DP), intent(in) :: a, b, c, ca, cb, cc
         real(DP), intent(out), dimension(3,3) :: tca, tac, tab, tba, &
        & tcb, tbc, grc, gkc, gra, gka, grb, gkb
! ===================== added by K. Tagami ============================= 12.0A
         logical  :: use_altv_rltv
         real(DP), dimension(3,3) :: altv, rltv
! ====================================================================== 12.0A

         call nslatc(grc,gkc)

! =============================- modified by K. Tagami =============== 12.0A
!        call nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)

         if ( use_altv_rltv ) then
           call nslata_kt( a,b,c,ca,cb,cc,tca,tac,gra,gka, altv, rltv )
         else
            call nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
         endif
! ======================================================================= 12.0A

         call nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)

         if(jpr >= 3) then
            if(printable) then
               write(nfout,*) ' '
               write(nfout,*) ' === C coordinate system ==='
               write(nfout,*) ' (Cartesian Coordinate System)'
               write(nfout,*) ' metric tensors :  (real space) (l**2) ',&
             & ' (reciprocal space) (1/l**2)'
              do  i=1,3
                 write(nfout,100) (grc(i,j),j=1,3),(gkc(i,j),j=1,3)
              end do
              write(nfout,*) &
             & ' =================================================='
            end if
         end if

         if(jpr >= 1) then
            if(printable) then
               write(nfout,*) ' '
               write(nfout,*) ' === A Coordinate system ==='
               write(nfout,*) &
              & ' (Conventional Crystalline Coordinate System)'
               write(nfout,*) ' transformation matrices Tca and Tac'
               do i=1, 3
                  write(nfout,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
               end do
               write(nfout,*) &
              & ' metric tensors :  (real space) (l**2)', &
              & ' (reciprocal space) (1/l**2)'
               do i=1, 3
                  write(nfout,100) (gra(i,j),j=1,3),(gka(i,j),j=1,3)
               end do
               write(nfout,*) &
              &' ================================================'
            end if
            if(il == -1) then
               if(printable) write(nfout,*) ' < Trigonal lattice (R) >'
            else if(il == 0) then
               if(printable) write(nfout,*) ' < Hexagonal lattice > '
            else if(il == 0 .or. il == 1) then
               if(printable) write(nfout,*) ' < Primitive lattice >'
            else if(il == 2) then
               if(printable) write(nfout,*) &
              & ' < face centered lattice (F) >'
            else if(il == 3) then
               if(printable) write(nfout,*) &
              & ' < body centered lattice (I) >'
            else if(il == 4) then
               if(printable) write(nfout,*) &
              &' < one(ab)-face centered lattice (C) >'
            else if(il == 5) then
               if(printable) write(nfout,*) &
              & ' < one(bc)-face centered lattice (C) >'
            else if(il == 6) then
               if(printable) write(nfout,*) &
              & ' < one(ca)-face centered lattice (C) >'
            end if
         end if

         if(jpr >= 2) then
            if(printable) then
               write(nfout,*) ' === B Coordinate system ==='
               write(nfout,*) &
              & ' (Primitive Crystalline Coordinate System)'
               write(nfout,*) ' transformation matrices Tab and Tba'
               do  i=1, 3
                  write(nfout,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
               end do
               write(nfout,*) ' transformation matrices Tcb and Tbc'
               do i=1, 3
                  write(nfout,100) (tcb(i,j),j=1,3), (tbc(i,j),j=1,3)
               end do
               write(nfout,*) &
             & ' metric tensors :  (real space) (l**2)      ', &
             & ' (reciprocal space) (1/l**2)'
               do i=1, 3
                  write(nfout,100) (grb(i,j),j=1,3),(gkb(i,j),j=1,3)
               end do
               write(nfout,*) &
             & ' ================================================'
            end if
         end if

        if(jpr >= 0) then
           if(printable) then
              write(nfout,*) ' transformation matrices Tab and Tba'
              do  i=1, 3
                 write(nfout,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
              end do
           end if
        end if
  100   format(1h ,3f12.6,3x,3f12.6)
      end subroutine nslat3
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
!
!#12  input:     a,b,c,ca,cb,cc: lattice parameters
!#12  output:    tca(3,3), tac(3,3): transformation matix
!#12                                 between C and A systems
!#12             gra(3,3),gka(3,3): metric tensors in A system
!#12  external:  nsmetr
!
!#21  to define the A coordinate system
!
!#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nslata(a,b,c,ca,cb,cc,tca,tac,gra,gka)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer :: i, j, k, ier
         real(DP), intent(in) :: a, b, c, ca, cb, cc
         real(DP), intent(out), dimension(3,3) :: tca, tac, gra, gka
         real(DP) :: eps, vv, sc, p, q, er1, eg1

         eps = 1.0d-8

         vv = dsqrt(1.0d0-ca*ca-cb*cb-cc*cc+2.0d0*ca*cb*cc)
         sc = dsqrt(1.0d0-cc*cc)
         p = (ca-cb*cc)/sc
         q = vv/sc

         tca(1,1) = a
         tca(2,1) = 0.0d0
         tca(3,1) = 0.0d0
         tca(1,2) = b*cc
         tca(2,2) = b*sc
         tca(3,2) = 0.0d0
         tca(1,3) = c*cb
         tca(2,3) = c*p
         tca(3,3) = c*q

         tac(1,1) = 1.0d0/a
         tac(2,1) = 0.0d0
         tac(3,1) = 0.0d0
         tac(1,2) = -1.0d0*cc/(a*sc)
         tac(2,2) = 1.0d0/(b*sc)
         tac(3,2) = 0.0d0
         tac(1,3) = (ca*cc-cb)/(a*vv*sc)
         tac(2,3) =-p/(b*vv)
         tac(3,3) = sc/(c*vv)

!       error check.
         ier=0
         do j=1, 3
            do i=1, 3
               er1 = 0.0d0
               eg1 = 0.0d0
               do k=1, 3
                  er1 = er1 + tca(i,k)*tac(k,j)
                  eg1 = eg1 + tac(i,k)*tca(k,j)
               end do
               if(i == j) then
                  if(dabs(er1-1.0d0) > eps .or. &
                 &   dabs(eg1-1.0d0) > eps) then
                     if(printable) write(nfout,*) 'i,j=',i,j,'  e=',er1,eg1
                     ier=1
                  end if
               else
                  if(dabs(er1) > eps .or.  &
                 &   dabs(eg1) > eps) then
                     if(printable)  write(nfout,*) 'i,j=',i,j,'  e=',er1,eg1
                     ier=1
                  end if
               end if
            end do
         end do

         if(ier /= 0) then
            if(printable) then
               write(nfout,*) ' '
               write(nfout,*) ' === A Coordinate system ==='
               write(nfout,*) &
              &' (Conventional Crystalline Coordinate System)'
               write(nfout,*) ' transformation matrices Tca and Tac'
               do i=1, 3
                 write(nfout,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
               end do
            end if
            stop ' === stop in sub.nslata ==='
         end if

 100 format(1h ,3f12.5,3x,3f12.5)

         call nsmetr(tca,tac,gra,gka)

      end subroutine nslata
! =========================== added by K. Tagami ======================= 12.0A
      subroutine nslata_kt( a,b,c,ca,cb,cc,tca,tac,gra,gka, &
     &                      altv, rltv )
         use m_Const_Parameters, only :  DP, PAI2
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
!$ use omp_lib
         implicit none
! inout
         real(DP), intent(in) :: a, b, c, ca, cb, cc
         real(DP), intent(out), dimension(3,3) :: tca, tac, gra, gka
         real(DP), intent(in), dimension(3,3) :: rltv, altv
!local
         real(DP) :: eps, eg1, er1
         integer  :: i, j, k, ier
! begin
         eps = 1.0d-8

         tca(1,1:3) = altv(1,1:3)
         tca(2,1:3) = altv(2,1:3)
         tca(3,1:3) = altv(3,1:3)

         do i=1, 3
            do k=1, 3
               tac(k,i) = rltv(i,k) /PAI2     ! confirmed in the case of graphene
            end do
         end do

! ---   error check.
        ier=0

        do j=1, 3
           do i=1, 3
              er1 = 0.0d0
              eg1 = 0.0d0
              do k=1, 3
                 er1 = er1+tca(i,k)*tac(k,j)
                 eg1 = eg1+tac(i,k)*tca(k,j)
              end do

              if(i == j) then
                 if(dabs(er1-1.0d0) > eps .or. &
                 &  dabs(eg1-1.0d0) > eps) then
                    if(printable) write(nfout,*) 'i,j=',i,j,'  e=',er1,eg1
                    ier = 1
                 end if
              else
                 if(dabs(er1) > eps .or. &
                 &  dabs(eg1) > eps) then
                    if(printable) write(nfout,*) 'i,j=',i,j,'  e=',er1,eg1
                    ier = 1
                 end if
              end if
           end do
        end do

        if(ier /= 0) then
           if(printable) then
              write(nfout,*) ' '
              write(nfout,*) &
             & ' === A Coordinate system ==='
              write(nfout,*) &
             & ' (Conventional Crystalline Coordinate System)'
              write(nfout,*) &
             & ' transformation matrices Tca and Tac'
              do i=1, 3
                write(nfout,100) (tca(i,j),j=1,3), (tac(i,j),j=1,3)
              end do
           end if
           stop ' === stop in sub.nslata ==='
        end if
 100    format(1h ,3f12.5,3x,3f12.5)

        call nsmetr(tca,tac,gra,gka)

      end subroutine nslata_kt
!c ==================================================================== 12.0A

!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)
!
!#12  input:     il                 : lattice type
!#12             tca(3,3), tac(3,3) : transformation matrix (C:A)
!#12  output:    tab(3,3), tba(3,3) : transformation matrix (A:B)
!#12             tcb(3,3), tbc(3,3) : transformation matrix (C:B)
!#12             grb(3,3), gkb(3,3) : metric tensors
!#12  external:  nsmetr
!
!#21  to define the B coordinate system
!
!c#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nslatb(il,tca,tac,tab,tba,tcb,tbc,grb,gkb)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: il
         integer :: i, j, k
         real(DP), intent(in), dimension(3,3) :: tca, tac
         real(DP), intent(out), dimension(3,3) :: tcb, tbc, grb, gkb, tab, tba
         real(DP) :: one, two, p1, q1, q2, x, y

         one = 1.0d0
         two = 2.0d0

         p1 = one/2.0d0
         q1 = one/3.0d0
         q2 = two/3.0d0

         if(il == -1) then
            if(printable)   write(nfout,*) ' < Trigonal lattice (R) >'
            tab(1,1) = q2
            tab(2,1) = q1
            tab(3,1) = q1
            tab(1,2) =-1.0d0*q1
            tab(2,2) = q1
            tab(3,2) = q1
            tab(1,3) =-1.0d0*q1
            tab(2,3) =-1.0d0*q2
            tab(3,3) = q1
            tba(1,1) = 1.0d0
            tba(2,1) =-1.0d0
            tba(3,1) = 0.0d0
            tba(1,2) = 0.0d0
            tba(2,2) = 1.0d0
            tba(3,2) =-1.0d0
            tba(1,3) = 1.0d0
            tba(2,3) = 1.0d0
            tba(3,3) = 1.0d0
         else if(il == 0) then
!       write(6,*) ' < Hexagonal lattice > '
            tab(1:3,1:3) = 0.0d0
            tba(1:3,1:3) = 0.0d0
            do i=1, 3
               tab(i,i)= 1.0d0
               tba(i,i)= 1.0d0
            end do
         else if(il == 0 .or. il == 1) then
!       write(6,*) ' < Primitive lattice >'
            tab(1:3,1:3) = 0.0d0
            tba(1:3,1:3) = 0.0d0
            do i=1, 3
               tab(i,i)= 1.0d0
               tba(i,i)= 1.0d0
            end do
         else if(il == 2) then
!       write(6,*) ' < face centered lattice (F) >'
            tab(1,1) = 0.0d0
            tab(2,1) = p1
            tab(3,1) = p1
            tab(1,2) = p1
            tab(2,2) = 0.0d0
            tab(3,2) = p1
            tab(1,3) = p1
            tab(2,3) = p1
            tab(3,3) = 0.0d0
            tba(1,1) =-1.0d0
            tba(2,1) = 1.0d0
            tba(3,1) = 1.0d0
            tba(1,2) = 1.0d0
            tba(2,2) =-1.0d0
            tba(3,2) = 1.0d0
            tba(1,3) = 1.0d0
            tba(2,3) = 1.0d0
            tba(3,3) =-1.0d0
         else if(il == 3) then
!       write(6,*) ' < body centered lattice (I) >'
            tab(1,1) =-1.0d0*p1
            tab(2,1) = p1
            tab(3,1) = p1
            tab(1,2) = p1
            tab(2,2) =-p1
            tab(3,2) = p1
            tab(1,3) = p1
            tab(2,3) = p1
            tab(3,3) =-p1
            tba(1,1) = 0.0d0
            tba(2,1) = 1.0d0
            tba(3,1) = 1.0d0
            tba(1,2) = 1.0d0
            tba(2,2) = 0.0d0
            tba(3,2) = 1.0d0
            tba(1,3) = 1.0d0
            tba(2,3) = 1.0d0
            tba(3,3) = 0.0d0
         else if(il == 4) then
!       write(6,*) ' < one(ab)-face centered lattice (C) >'
            tab(1,1) = p1
            tab(2,1) =-1.0d0*p1
            tab(3,1) = 0.0d0
            tab(1,2) = p1
            tab(2,2) = p1
            tab(3,2) = 0.0d0
            tab(1,3) = 0.0d0
            tab(2,3) = 0.0d0
            tab(3,3) = 1.0d0
            tba(1,1) = 1.0d0
            tba(2,1) = 1.0d0
            tba(3,1) = 0.0d0
            tba(1,2) =-1.0d0
            tba(2,2) = 1.0d0
            tba(3,2) = 0.0d0
            tba(1,3) = 0.0d0
            tba(2,3) = 0.0d0
            tba(3,3) = 1.0d0
         else if(il == 5) then
!       write(6,*) ' < one(bc)-face centered lattice (C) >'
            tab(1,1) = 1.0d0
            tab(2,1) = 0.0d0
            tab(3,1) = 0.0d0
            tab(1,2) = 0.0d0
            tab(2,2) = p1
            tab(3,2) =-1.0d0*p1
            tab(1,3) = 0.0d0
            tab(2,3) = p1
            tab(3,3) = p1
            tba(1,1) = 1.0d0
            tba(2,1) = 0.0d0
            tba(3,1) = 0.0d0
            tba(1,2) = 0.0d0
            tba(2,2) = 1.0d0
            tba(3,2) = 1.0d0
            tba(1,3) = 0.0d0
            tba(2,3) =-1.0d0
            tba(3,3) = 1.0d0
         else if(il == 6) then
!       write(6,*) ' < one(ca)-face centered lattice (C) >'
            tab(1,1) = p1
            tab(2,1) = 0.0d0
            tab(3,1) = p1
            tab(1,2) = 0.0d0
            tab(2,2) = 1.0d0
            tab(3,2) = 0.0d0
            tab(1,3) =-1.0d0*p1
            tab(2,3) = 0.0d0
            tab(3,3) = p1
            tba(1,1) = 1.0d0
            tba(2,1) = 0.0d0
            tba(3,1) =-1.0d0
            tba(1,2) = 0.0d0
            tba(2,2) = 1.0d0
            tba(3,2) = 0.0d0
            tba(1,3) = 1.0d0
            tba(2,3) = 0.0d0
            tba(3,3) = 1.0d0
         else
            if(printable) then
               write(nfout,*) ' il=',il
               write(nfout,*) ' === stop in sub.nslatb. (il) ==='
            end if
            stop 'error in sub.nslatb'
         end if

         do i=1, 3
            do j=1, 3
               x = 0.0d0
               y = 0.0d0
               do  k =1, 3
                  x = x+tca(i,k)*tab(k,j)
                  y = y+tba(i,k)*tac(k,j)
               end do
               tcb(i,j) = x
               tbc(i,j) = y
            end do
         end do
!     write(6,*) ' === B Coordinate system ==='
!     write(6,*) ' (Primitive Crystalline Coordinate System)'
!     write(6,*) ' transformation matrices Tab and Tba'
!     do 40 i=1,3
!  40 write(6,100) (tab(i,j),j=1,3), (tba(i,j),j=1,3)
!     write(6,*) ' transformation matrices Tcb and Tbc'
!     do 42 i=1,3
!  42 write(6,100) (tcb(i,j),j=1,3), (tbc(i,j),j=1,3)
! 100 format(1h ,3f12.5,3x,3f12.5)
!
         call nsmetr(tcb,tbc,grb,gkb)
      end subroutine nslatb
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nslatc(grc,gkc)
!
!#12  noinput:
!#12  output:    grc(3,3): metric tensor in real space
!#12             gkc(3,3): metric tensor in reciprocal space
!#12  noexternal:
!
!#21  to get metric tensors in C (Cartesian) coordinate system
!
!#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nslatc(grc,gkc)
        use m_Const_Parameters, only :  DP, PAI2
        implicit none
        integer :: i, j
        real(DP), intent(out), dimension(3,3) :: grc, gkc
        real(DP) ::  pai22

        pai22 = PAI2**2

        grc(1:3,1:3) = 0.0d0
        gkc(1:3,1:3) = 0.0d0

        do i=1,3
           grc(i,i) = 1.0d0
           gkc(i,i) = pai22
        end do
!
!     write(6,*) ' '
!     write(6,*) ' === C coordinate system ==='
!     write(6,*) ' (Cartesian Coordinate System)'
!     write(6,*) ' metric tensors :  (real space) (l**2)      ',
!    &           ' (reciprocal space) (1/l**2)'
!     do 30 i=1,3
!  30 write(6,100) (grc(i,j),j=1,3),(gkc(i,j),j=1,3)
! 100 forma0t(1h ,3f12.6,3x,3f12.6)
!     write(6,*) ' =================================================='
!
      end subroutine nslatc
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nslatr(il,r,ind)
!
!     il :  lattice type
!     ind = 0 :  r is lattice vector.
!           1 :  r in not lattice vector.
!
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : factor_1_3
         implicit none
         integer, intent(in) :: il
         integer, intent(out) :: ind
         integer :: i, flag
         integer, dimension(3) :: j
         real(DP), intent(in), dimension(3) :: r
         real(DP) :: one, half, thrd, eps, aaa, s
         data one/1.0d0/, half/0.5d0/, eps/1.0d-4/


         if(il <= -2 .or. il >= 5) then
            if(printable) write(nfout,*) ' il=',il
            stop '=== stop in sub.nslatr. (il) ==='
         end if

         ind = 1
         do i=1, 3
            flag = 0
            aaa = dmod(dabs(r(i)),1.0d0)
            if(aaa > eps .and. abs(1.0d0-aaa) > eps) then
               flag = 1
               exit
            end if
         end do
         if(flag == 0) then
            ind = 0
            return
         end if

         if(il == -1) then
            do i=1, 3
               aaa = dmod(dabs(r(i)),factor_1_3)
               if(aaa > eps .and. dabs(factor_1_3-aaa) > eps) then
                  return
               end if
            end do
            aaa =dmod(dabs(-r(1)+r(2)+r(3)),1.0d0)
            if(aaa < eps .or. dabs(1.0d0-aaa) < eps) then
               ind = 0
            end if
         else if(il <= 1) then
            return
         else if(il == 2 .or. il == 3) then
            do i=1, 3
               aaa = dmod(dabs(r(i)),0.5d0)
               if(aaa.gt.eps .and. dabs(0.5d0-aaa) > eps) then
                  return
               endif
            end do
            if(il == 2) then
               aaa = dmod(dabs(r(1)+r(2)+r(3)),1.0d0)
               if(aaa < eps .or. dabs(1.0d0-aaa) < eps) then
                  ind = 0
               end if
            else if(il == 3) then
               do i=1, 3
                  s = sign(0.1d0,r(i))
                  j(i) = iabs(int(2.0d0*r(i)+s))
               end do
               if((mod(j(1),2) == 0 .and. mod(j(2),2) == 0 .and. &
                & mod(j(3),2) == 0) .or. (mod(j(1),2) == 1 .and. &
                & mod(j(2),2) == 1 .and. mod(j(3),2) ==1)) then
                  ind = 0
               end if
            end if
         else if(il == 4) then
            do i=1, 2
               aaa = dmod(dabs(r(i)),0.5d0)
               if(aaa > eps .and. dabs(0.5d0-aaa) > eps) then
                  return
               end if
            end do
            aaa = dmod(dabs(r(3)),1.0d0)
            if(aaa > eps .and. dabs(1.0d0-aaa) > eps) then
               return
            end if
            aaa = dmod(dabs(r(1)+r(2)),1.0d0)
            if(aaa < eps .or. dabs(1.0d0-aaa) < eps) then
               ind = 0
            end if
         end if
      end subroutine nslatr
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
!
!#12  input:       il : lattice type
!#12              ng1 : # of group elements (order of the space group)
!#12             ig01 : list vector for getting the element code
!#12  in-output:   a,b,c,ca,cb,cc: lattice parameters
!#12  noexternal:
!
!!#21  to check the consistency between the space group and
!#21                                  the lattice parameters
!!
!!#31  1990.04.12.:  n. hamada, a. yanase and k. terakura

! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(in) :: il, ng1
         integer, intent(in), dimension(48) :: ig01
         integer :: i, jj
         real(DP), intent(in) :: a
         real(DP), intent(out) :: b, c, ca, cb, cc
         logical :: la,lb,lc,lcub,ltet

         if(il <= 0) then
            b = a
            cc =-0.5d0
            ca = 0.0d0
            cb = 0.0d0
         else
            la = .false.
            lb = .false.
            lc = .false.
            lcub = .false.
            ltet = .false.
            do  i=1, ng1
               jj=ig01(i)
               if(jj == 4)  lc = .true.
               if(jj == 28) lc = .true.
               if(jj == 2)  la = .true.
               if(jj == 26) la = .true.
               if(jj == 27) lb = .true.
               if(jj == 3)  lb = .true.
               if(jj == 5)  lcub = .true.
               if(jj == 21) ltet = .true.
               if(jj == 48) ltet = .true.
            end do
            if(lc) ca = 0.0d0
            if(lc) cb = 0.0d0
            if(la) cb = 0.0d0
            if(la) cc = 0.0d0
            if(lb) cc = 0.0d0
            if(lb) ca = 0.0d0
            if(ltet) b = a
            if(lcub) b = a
            if(lcub) c = a
         end if
      end  subroutine nslatz
! ===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsmetr(tca,tac,gra,gka)
!
!#12  input:     tca(3,3), tac(3,3) : transformation matrix
!#12  output:    gra(3,3), gka(3,3)  : metric tensors
!#12  noexternal:
!
!#21  to get metric tensors
!
!#31  1990.01.18.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsmetr(tca,tac,gra,gka)
         use m_Const_Parameters, only :  DP, PAI2
         implicit none
         integer :: i, j, k
         real(DP), intent(in), dimension(3,3) :: tca, tac
         real(DP), intent(out), dimension(3,3):: gra, gka
         real(DP) ::  pai22, x, y

         pai22 = PAI2**2

         do i=1, 3
            do j=1, 3
               x = 0.0d0
               y = 0.0d0
               do k=1, 3
                  x = x+tca(k,i)*tca(k,j)
                  y = y+tac(i,k)*tac(j,k)
               end do
               gra(i,j) = x
               gka(i,j) = y*pai22
            end do
         end do
      end subroutine nsmetr
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsmlt1(nn,n,a,b,c)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(in) :: n, nn
         integer :: i, j, k
         real(DP), intent(in), dimension(nn,n) :: a, b
         real(DP), intent(out), dimension(nn,n) :: c
         real(DP) :: sum

         do i=1, n
            do j= 1, n
               sum = 0.0d0
              do k=1, n
                  sum = sum + a(i,k)*b(k,j)
               end do
               c(i,j) =sum
            end do
         end do
      end subroutine nsmlt1

! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!==*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsmult(jf,jpr,ng0,irot,im,iv)
!
!#12  input :          jf   : output file
!#12                  jpr   : output control
!#12                  ng0   : # of group elements
!#12         irot(3,3,48)   : rotation matrix
!#12  output:     im(48,48) : multiplication table
!#12              iv(48)    : inverse elements
!#13  noexternal
!
!#21  to get multiplication table from rotation matrices
!
!#31  1989.12.26.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsmult(jf,jpr,ng0,irot,im,iv)
         use m_Const_Parameters, only :  DP
!$ use omp_lib
         implicit none
         integer, intent(in) :: jf, jpr, ng0
         integer, intent(in), dimension(3,3,48) :: irot
         integer, intent(inout), dimension(48,48) :: im
         integer, intent(out), dimension(48) :: iv
         integer :: i, j, k, l, m, jj
         integer :: flag
         integer, dimension(3,3) :: ia

         do i=1, ng0
            do j=1, ng0
               do  k=1, 3
                  do l=1, 3
                     ia(k,l) = 0
                     do m=1, 3
                        ia(k,l) = ia(k,l)+irot(k,m,i)*irot(m,l,j)
                     end do
                  end do
               end do
               do  k=1, ng0
                  flag = 0
                  if(ia(1,1) == irot(1,1,k) .and. &
                  &  ia(2,1) == irot(2,1,k) .and. &
                  &  ia(3,1) == irot(3,1,k) .and. &
                  &  ia(1,2) == irot(1,2,k) .and. &
                  &  ia(2,2) == irot(2,2,k) .and. &
                  &  ia(3,2) == irot(3,2,k) .and. &
                  &  ia(1,3) == irot(1,3,k) .and. &
                  &  ia(2,3) == irot(2,3,k) .and. &
                  &  ia(3,3) == irot(3,3,k)) then
                     flag = 1
                     im(i,j) = k
                     exit
                  end if
               end do
               if(flag == 0) then
                  stop '=== stop in sub.nsmult. (im) ==='
               else
                  cycle
              end if
            end do
         end do

         do i=1, ng0
            do j=1, ng0
               jj = j
                flag = 0
               if(im(i,j) == 1) then
                  flag = 1
                  iv(i)=jj
                  exit
               end if
            end do
            if(flag == 0) then
               stop '=== stop in sub.nsmult. (iv) ==='
            else
               cycle
            end if
         end do

         if(jpr >= 3) then
            write(jf,*) '   '
            write(jf,*) '--- group multiplication table ---'
            write(jf,120) (j,j=1,24)
            write(jf,140)
            do i=1, 24
               write(jf,100) i,(im(i,j),j=1,24)
            end do
            if(ng0 > 24) then
               write(jf,*) ' '
               write(jf,120) (j,j=1,24)
               write(jf,140)
               do i=25, ng0
                  write(jf,100) i,(im(i,j),j=1,24)
               end do
               write(jf,*) ' '
               write(jf,120) (j,j=25,ng0)
               write(jf,140)
               do i=1, 24
                  write(jf,100) i,(im(i,j),j=25,ng0)
               end do
               write(jf,*) ' '
               write(jf,120) (j,j=25,ng0)
               write(jf,140)
               do i=25, ng0
                  write(jf,100) i,(im(i,j),j=25,ng0)
               end do
            end if

            write(jf,*) ' '
            write(jf,*) '--- invers elements ---'
            write(jf,120) (iv(j),j=1,ng0)
  100       format((i3,2x,24i3))
  120       format((5x,24i3))
  140       format(5x,72('-'))
         end if
      end subroutine nsmult
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nspace(jf,jpr,il,ngen,inv,igen,jgen,
!#11 &           ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1,
!#11 &           omove,euler1,inver1)
!
!#12  input:      jf : output file
!#12             jpr : print control
!#12              il : lattice type
!!#12             inv : parameter (0,1) for moving the origin
!#12         igen(3) : generator (rotation part)
!#12     jgen(2,3,3) : generator (nonprimitive translation part)
!#12  output:       ng00 : # of elements in Schoenflies notation
!#12          schoe0(48) : Schoenflies notation (character*5)
!#12                 ng1 : # of group elements
!#12                ig01 : element code
!#12                 ta1 : nonprimitive translation vector (A system)
!#12                 ra1 : rotation matrix in real space (A system)
!#12                 sa1 : rotation matrix in reciprocal space (A system)
!#12                       ra1*sa1=unit matrix
!#12          im1(48,48) : multiplication table
!#12             iv1(48) : inverse element
!#12            omove(3) : translation of the origin from the initial
!#12        euler1(3,48) : Euler angle
!#12        inver1(  48) : index for inversion operation (=1 or -1)
!#13  external:  nsgrp1, tspaca
!
!#21  to get space group
!
!#31  1990.11.12.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
! ================================== modified by K. Tagami ============ 12.0A
!     subroutine nspace(jf,jpr,il,ngen,inv,igen,jgen,
!    &                  ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1,
!    &                  omove,euler1,inver1)
      subroutine nspace(jf,jpr,il,ngen,inv,igen,jgen, &
     &                  ng00,schoe0,ng1,ig01,ta1,ra1,sa1,im1,iv1, &
     &                  omove,euler1,inver1, &
     &                  use_trs, tac, tca, tab, tba, &
     &                  gen_name_in_carts )
! ====================================================================== 12.0A
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : schoen, jones, euler, rot, &
        & ieuler, irot, ir1234, im0, iv0, ng0_0
! 2021.12.20 : OpenMP directives by T. Hamada
!$ use omp_lib
         implicit none
         integer, intent(in) :: jf, jpr, il, ngen, inv
         integer, intent(in), dimension(3) :: igen
         integer, intent(in), dimension(2,3,3) :: jgen
         integer, intent(out) :: ng00, ng1
         integer, intent(out), dimension(48)  :: ig01, iv1, inver1
         real(DP), intent(out), dimension(3,48) :: ta1
         integer, intent(out), dimension(48,48) :: im1
         real(DP), intent(out), dimension(3,3,48) :: ra1, sa1
         integer :: i, j, k, ig0, ii
         integer, dimension(48) :: ig10
         integer, dimension(2,3,48) :: jg1
         integer,  dimension(3,3,48) :: irotr1, irotk1
         integer, dimension(2,3) :: movo
         real(DP), intent(out), dimension(3) :: omove
         real(DP), intent(out), dimension(3,48) :: euler1
         character(len=5), dimension(48) :: schoe0

! ====================================== added by K. Tagami ============ 12.0A
         logical use_trs, gen_name_in_carts
         real(DP), dimension(3,3) ::  tac, tca, tab, tba
! ====================================================================== 12.0A
!====================================== modified by K. Tagami ========= 12.0A
!      call nsgrp1(jf,jpr,il,ngen,inv,igen,jgen,
!     &            schoen,jones,euler,rot,ieuler,irot,ir1234,
!     &            im0,iv0,ng0,ng1,ig10,ig01,jg1,
!     &            irotr1,irotk1,im1,iv1,movo)
         call nsgrp1(jf,jpr,il,ngen,inv,igen,jgen, &
     &            schoen,jones,euler,rot,ieuler,irot,ir1234, &
     &            im0,iv0,ng0_0,ng1,ig10,ig01,jg1, &
     &            irotr1,irotk1,im1,iv1,movo, &
     &            use_trs, tac, tca, tab, tba, &
     &            gen_name_in_carts )
! ======================================================================= 12.0A
!$omp parallel
         ta1(1:3,1:ng1) = real(jg1(1,1:3,1:ng1),kind=DP)/real(jg1(2,1:3,1:ng1),kind=DP)
         ra1(1:3,1:3,1:ng1) = real(irotr1(1:3,1:3,1:ng1),kind=DP)
         sa1(1:3,1:3,1:ng1) = real(irotk1(1:3,1:3,1:ng1),kind=DP)
!$omp end parallel

         do k=1, ng1
            ig0 = ig01(k)
            if(il <= 0) then
               if(ig0 <= 12) then
                  ii = ig0
                  inver1(k) = 1
               else
                  ii = ig0-12
                  inver1(k) =-1
               end if
            else
               if(ig0 <= 24) then
                  ii = ig0
                  inver1(k) = 1
               else
                  ii =ig0-24
                  inver1(k) = -1
               end if
            end if
            euler1(1:3,k) = euler(1:3,ii)
         end do

         ng00 = ng0_0
         schoe0(1:ng0_0) = schoen(1:ng0_0)

         omove(1:3) = real(movo(1,1:3),kind=DP)/real(movo(2,1:3),kind=DP)
!         omove(2) = real(movo(1,2),kind=DP)/real(movo(2,2),kind=DP)
!         omove(3) = real(movo(1,3),kind=DP)/real(movo(2,3),kind=DP)

!     for being compatible with the tspace package
         call tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
      end subroutine nspace
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nspbge(p1,p2,ind)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(out) :: ind
         real(DP), dimension(3) ::  p1, p2
         real(DP), dimension(3) :: d
         real(DP) :: one, eps, x, y, z

         eps = 1.0d-5
         ind = 1
         d(1) = dabs(p1(1)-p2(1)) + 0.5d-5
         if(dmod(d(1),1.0d0) <= 1.0d-5) return
         d(2) = dabs(p1(2)-p2(2)) + 0.5d-5
         if(dmod(d(2),1.0d0) <= 1.0d-5) return
         d(3) = dabs(p1(3)-p2(3)) + 0.5d-5
         if(dmod(d(3),1.0d0) <= 1.0d-5) return
         ind = 0
!         x = dabs(p1(1)-p2(1))+0.5d-5
!         y = dabs(p1(2)-p2(2))+0.5d-5
!         z = dabs(p1(3)-p2(3))+0.5d-5
!      d = dabs(p1 - p2) + 0.5d-5
!      if(dmod(d(1),1.0d0).le. 1.0d-5 .and. &
!      &  dmod(d(2),1.0d0).le. 1.0d-5 .and. &
!      &  dmod(d(3),1.0d0).le. 1.0d-5      ) then
!            ind = 0
!         else
!            ind = 1
!         end if
      end subroutine nspbge

      subroutine nspbgei(p1,p2,ind)
         implicit none
         integer, intent(in), dimension(3) :: p1, p2
         integer, intent(out) :: ind
         integer, dimension(3) :: d
         integer :: x, y, z
         integer :: one = 100000

!         x = iabs(p1(1)-p2(1))
!         y = iabs(p1(2)-p2(2))
!         z = iabs(p1(3)-p2(3))
         ind = 1
         d(1) = iabs(p1(1)-p2(1))
         if(mod(d(1),100000) /=0) return
         d(2) = iabs(p1(2)-p2(2))
         if(mod(d(2),100000) /=0) return
         d(3) = iabs(p1(3)-p2(3))
         if(mod(d(3),100000) /=0) return
         ind = 0
!        if(mod(d(1),100000).eq.0 .and. &
!       &   mod(d(2),100000).eq.0 .and. &
!       &   mod(d(3),100000).eq.0      ) then
!            ind = 0
!         else
!            ind = 1
!         end if
      end subroutine nspbgei

! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im,ir,schoen,mo,ng,ig,jg)
!!c
!#12  input:    jf  : output file
!#12           jpr  : output control
!#12            il  : lattice type
!#12           inv  : flug for moving the origin (0:nomove, 1:move)
!#12          ngen  : # of generator
!#12          igen  : rotation of generater
!#12          jgen  : nonprimitive translation vector of generator
!#12            im  : multiplication table
!#12            ir  : rotation matrix
!#12        schoen  : schoenfries index
!#12  output:   mo(2,3)    : translation of the origin
!#12            ng         : # of operations (order of group)
!#12            ig(48)     : rotation of group element
!#12            jg(2,3,48) : nonprimitive translation of the element
!#13  external: sub.nsrot1, nssum1
!
!#21  to get space group
!
!#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nspgrp(jf,jpr,il,inv,ngen,igen,jgen,im,ir,schoen, &
     &                  mo,ng,ig,jg)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: jf, jpr, il, inv, ngen
         integer, intent(in), dimension(ngen) :: igen
         integer, intent(in), dimension(2,3,ngen) :: jgen
         integer, intent(in), dimension(48,48) :: im
         integer, intent(in), dimension(2,3,48) :: ir
         integer, intent(out) :: ng
         integer, intent(out), dimension(2,3) :: mo
         integer, intent(out), dimension(:) :: ig
         integer, intent(out), dimension(2,3,48) :: jg
         integer :: i, j, k, n, nn, i1, iw, nprm
         integer, dimension(48) :: iflg
         integer, dimension(2,3,48) :: jw
         integer, dimension(2,3) :: jww, mo1, jw1
         integer, dimension(2,3,4) :: lprm
         character(len=5), dimension(48) :: schoen

! --- register the generator ---

         iflg(1) = 1
         jw(1,1:3,1) = 0
         jw(2,1:3,1) = 1

         iflg(2:48) = 0

         do i=1, ngen
            if(il >= 1 .and. igen(i) > 48) then
               if(printable) then
                  write(nfout,*) ' igen=',igen(i),' > 48 : error'
               end if
               stop 'error:sub.nspgrp (igen)'
            end if
            if(il <= 0 .and. igen(i) > 24) then
               if(printable) then
                  write(nfout,*) ' igen=',igen(i),' > 24 : error'
               end if
               stop 'error:sub.nspgrp (igen)'
            end if
            iflg(igen(i))=1
            jw(1:2,1:3,igen(i)) = jgen(1:2,1:3,i)
         end do
         n = 0
         do i=1, 48
            n=n+iflg(i)
         end do
         if(jpr >= 0)      write(jf,*) ' '
         if(jpr >= 0)      write(jf,*) ' ----- generator -------'
         if(jpr >= 0) then
             nn=0
             do  i=1, 48
                if(iflg(i) == 1) then
                   nn=nn+1
                   write(jf,100) nn,i,schoen(i),((jw(j,k,i),j=1,2),k=1,3)
                end if
             end do
         end if

!--- produce group elements by multiplications

    1    ng = n
         do i=1, 48
            if(iflg(i) == 1) then
               do j=1, 48
                  if(iflg(j) == 1) then
!                   iw=im(i,j)
                     iw = im(j,i)
                     if(iflg(iw) /= 1) then
                        iflg(iw) = 1
                        call nsrot1(jw(1,1,i),ir(1,1,j),jww)
                        call nssum1(1,jww,jw(1,1,j),jw(1,1,iw))
                     end if
                  end if
               end do
            end if
         end do
         n = 0
         do i=1, 48
            n = n+iflg(i)
         end do
         if(n > ng) go to 1

! --- register group elements

         n = 0
         do  i=1, 48
            if(iflg(i) == 1) then
               n = n+1
               ig(n) = i
               jg(1:2,1:3,n) = jw(1:2,1:3,i)
            end if
         end do
! --- move the origin

         mo(1,1:3) = 0
         mo(2,1:3) = 1
         if(inv /= 0) then
            i1 = 0
            do i=1, ng
               if(il <= 0 .and. ig(i) == 13) i1 = i
               if(il >= 1 .and. ig(i) == 25) i1 = i
            end do
            if(i1 /= 0) then
               do i=1, 3
                  mo(1,i) = jg(1,i,i1)
                  mo(2,i) = jg(2,i,i1)**2
                  call nsrduc(mo(1,i),mo(2,i))
               end do
               if(jpr >= 1) then
                  if(printable) write(nfout,200) ((mo(i,j),i=1,2),j=1,3)
               end if
               do i=1, ng
                  call nsrot1(mo,ir(1,1,ig(i)),mo1)
                  call nssum1( 1,jg(1,1,i),mo1,jw1)
                  call nssum1(-1,jw1,mo,jg(1,1,i))
               end do
            end if
         end if

         call nsprmv(il,nprm,lprm)
         call nssdjg(jf,jpr,nprm,lprm,ng,jg)

! --- print


         if(jpr >= 2) then
            write(jf,*) ' '
           write(jf,*) ' ----- group elements -------'
           do i=1, ng
              write(jf,100) i,ig(i),schoen(ig(i)), &
             & ((jg(j,k,i),j=1,2),k=1,3)
           end do
        end if
  100   format(i5,'   (',i2,')',3x,a5,5x, &
       &         '(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )
  200   format(/'    move the origin by ','  (',2(i3,' /',i3,'  ,'), &
       &                                         i3,' /',i3,' )'      )
      end subroutine nspgrp

! ==================================== added by K. Tagami =============== 12.0A
      subroutine nspgrp_kt( jf, jpr, il, inv, ngen, igen, jgen, &
     &                      im, ir, schoen, &
     &                      mo,ng,ig,jg, &
     &                      use_trs, tac, tca, tab, tba, &
     &                      gen_name_in_carts)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: jf, jpr, il, inv, ngen
         integer, intent(in), dimension(3) :: igen
         integer, intent(in), dimension(2,3,3) :: jgen
         integer, intent(in), dimension(48,48) :: im
         integer, intent(in), dimension(3,3,48) :: ir
         integer, intent(out) :: ng
         integer, intent(out), dimension(2,3) :: mo
         integer, intent(out), dimension(48) :: ig
         integer, intent(out), dimension(2,3,48) :: jg
         integer, dimension(48) :: iflg
         integer, dimension(2,3,48) :: jw
         integer, dimension(2,3) :: jww, mo1,jw1
         integer, dimension(2,3,4) ::  lprm
         character(len=5), dimension(48) :: schoen
         logical :: use_trs, gen_name_in_carts
         logical :: with_inverse
         integer, dimension(48) :: iflg_tmp
         integer :: i, j, k, n, nn, iw, i1, nprm
         integer :: m1, m2
         integer, dimension(3,3,48) :: ir_sysA, ir_sysB
         real(DP) ::  c1
         real(DP), dimension(3,3) ::  tac, tca, tab, tba

         integer i_shift

! -- make op matrix in the A system
         if ( gen_name_in_carts ) then
            do k=1, 48
               do i=1, 3
                  do j=1, 3
                     c1 = 0.0d0
                     do m1=1, 3
                        do m2=1, 3
                           c1 = c1 + tac(i,m1)*ir(m1,m2,k)*tca(m2,j)
                        end do
                     end do
                     ir_sysA(i,j,k) = nint(c1)
                  end do
               end do
            end do
! -- make op matrix in the B system
            do k=1, 48
               do i=1, 3
                  do j=1, 3
                     c1 = 0.0d0
                     do m1=1, 3
                        do m2=1, 3
                           c1 = c1 + tab(i,m1)*ir_sysA(m1,m2,k)*tba(m2,j)
                        end do
                     end do
                     ir_sysB(i,j,k) = nint(c1)
                  end do
               end do
            end do
         else
            ir_sysB = ir
         endif

! --- register the generator ---

         iflg(1)=1
         jw(1,1:3,1) = 0
         jw(2,1:3,1) = 1

         iflg(2:48) = 0

         do i=1, ngen
            if(il >=1 .and. igen(i) > 48) then
               if(printable) then
                  write(nfout,*) ' igen=',igen(i),' > 48 : error'
               end if
               stop 'error:sub.nspgrp (igen)'
            end if
            if(il <= 0 .and. igen(i) > 24) then
               if(printable) then
                  write(nfout,*) ' igen=',igen(i),' > 24 : error'
               end if
               stop 'error:sub.nspgrp (igen)'
            end if
            iflg(igen(i)) = 1
            jw(1:2,1:3,igen(i))=jgen(1:2,1:3,i)
         end do

         n = 0
         do i=1, 48
            n = n+iflg(i)
         end do
         if(jpr >= 0) write(jf,*) ' '
         if(jpr >= 0) write(jf,*) ' ----- generator -------'
         if(jpr >= 0) then
            nn = 0
            do i=1, 48
               if(iflg(i) == 1) then
                  nn = nn+1
                  write(jf,100) nn,i,schoen(i),((jw(j,k,i),j=1,2),k=1,3)
                end if
            end do
         end if
! --- produce group elements by multiplications
    1    ng=n
         do i=1,48
            if(iflg(i) == 1) then
               do j=1, 48
                  if(iflg(j) == 1) then
!                    iw=im(i,j)
                     iw = im(j,i)
                     if(iflg(iw) /= 1) then
                        iflg(iw) = 1
                        call nsrot1(jw(1,1,i),ir_sysB(1,1,j),jww)
                        call nssum1(1,jww,jw(1,1,j),jw(1,1,iw))
                     end if
                  end if
               end do
            end if
         end do
         n = 0
         do i=1, 48
            n = n+iflg(i)
         end do
         if(n > ng) go to 1
! --- register group elements
         n = 0
         do i=1, 48
            if(iflg(i) == 1) then
               n = n+1
               ig(n) = i
               jg(1:2,1:3,n)=jw(1:2,1:3,i)
            end if
         end do

!      write(*,*) 'use_trs now = ',use_trs

! ----------------------------- check if opr has inv-symm --
         with_inverse = .false.
         iflg_tmp = 0

         if ( il >= 1 ) then
            if ( iflg(25)  == 1 ) then
               with_inverse = .true.
            endif
         else if ( il <= 0 ) then
            if ( iflg(13) == 1 ) then
               with_inverse = .true.
            endif
         endif

         if ( (.not. with_inverse) .and. use_trs ) then
            iflg_tmp = 0
            if ( il >= 1 ) then
               do i=1, 48
                  if ( i <= 24 ) then
                     i_shift = 24
                  else
                     i_shift = -24
                  endif
                  if ( iflg(i)  == 1 ) then
                     iflg_tmp(i) = 1
                     iflg_tmp(i+i_shift) = 1
                  endif
               end do
            else if ( il <= 0 ) then
               do i=1, 24
                  if ( i <= 12 ) then
                     i_shift = 12
                  else
                     i_shift = -12
                  endif
                  if ( iflg(i) == 1 ) then
                     iflg_tmp(i) = 1
                     iflg_tmp(i+i_shift) = 1
                  endif
               end do
            end if

            n = 0
            do i=1, 48
               if ( iflg_tmp(i) == 1 ) then
                  n = n+1
                  ig(n) = i
                  if ( iflg(i) == 1 ) then
                     jg(1:2,1:3,n) = jw(1:2,1:3,i)
                  else
                     if ( il >= 1 ) then
                        if ( i <= 24 ) then
                           i_shift = 24
                        else
                           i_shift = -24
                        endif
                     else
                        if ( i <= 12 ) then
                           i_shift = 12
                        else
                           i_shift = -12
                        endif
                     endif
                     jg(1:2,1:3,n)=jw(1:2,1:3,i+i_shift)
                  end if
               end if
            end do
            iflg = iflg_tmp
         end if

! --- move the origin

         mo(1,1:3) =0
         mo(2,1:3) =1
         if(inv /= 0) then
            i1=0
            do  i=1, ng
               if(il <= 0 .and. ig(i) == 13) i1=i
               if(il >= 1 .and. ig(i) == 25) i1=i
            end do
            if(i1 /= 0) then
               do i=1, 3
                  mo(1,i) = jg(1,i,i1)
                  mo(2,i) = jg(2,i,i1)*2
                  call nsrduc(mo(1,i),mo(2,i))
               end do
               if(jpr >= 1) then
                  if(printable) write(nfout,200) ((mo(i,j),i=1,2),j=1,3)
               end if
               do i=1, ng
                  call nsrot1(mo,ir_sysB(1,1,ig(i)),mo1)
                  call nssum1( 1,jg(1,1,i),mo1,jw1)
                  call nssum1(-1,jw1,mo,jg(1,1,i))
               end do
            end if
         end if

         call nsprmv(il,nprm,lprm)
         call nssdjg(jf,jpr,nprm,lprm,ng,jg)

! --- print

         if(jpr.ge.2) then
            write(jf,*) ' '
            write(jf,*) ' ----- group elements -------'
            do i=1, ng
               write(jf,100) i,ig(i),schoen(ig(i)), &
              & ((jg(j,k,i),j=1,2),k=1,3)
            end do
         end if
  100    format(i5,'   (',i2,')',3x,a5,5x, &
       &'(',2(i3,' /',i3,'  ,'),i3,' /',i3,' )'      )
  200    format(/'    move the origin by ','  (',2(i3,' /',i3,'  ,'), &
       &                                         i3,' /',i3,' )'      )
      end subroutine nspgrp_kt
! ===================================================================== 12.0A
!
!#11  sub.nsprmv(il,nv,lv)
!
!#12  input   :       il : lattice type
!#12  output  :       nv : # of vectors
!#12           lv(2,3,4) : primitive translation vector
!#12  noexternal
!
!#21  to get primitive translation vectors
!     e.g., (0/1,0/1,0/1), (1/2,1/2,1/2) for body centered lattice
!
!#31  1990.01.06.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsprmv(il,nv,lv)
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: il
         integer, intent(out) :: nv
         integer, intent(out), dimension(2,3,4) :: lv

         lv(1,1:3,1) = 0
         lv(2,1:3,1) = 1

         if(il == -1) then
            nv = 3
            lv(1,1,2) = 2
            lv(2,1,2) = 3
            lv(1,2,2) = 1
            lv(2,2,2) = 3
            lv(1,3,2) = 1
            lv(2,3,2) = 3
            lv(1,1,3) = 1
            lv(2,1,3) = 3
            lv(1,2,3) = 2
            lv(2,2,3) = 3
            lv(1,3,3) = 2
            lv(2,3,3) = 3
         else if(il== 0 .or. il == 1) then
            nv=1
         else if(il == 2) then
            nv = 4
            lv(1,1,2) = 1
            lv(2,1,2) = 2
            lv(1,2,2) = 1
            lv(2,2,2) = 2
            lv(1,3,2) = 0
            lv(2,3,2) = 1
            lv(1,1,3) = 1
            lv(2,1,3) = 2
            lv(1,2,3) = 0
            lv(2,2,3) = 1
            lv(1,3,3) = 1
            lv(2,3,3) = 2
            lv(1,1,4) = 0
            lv(2,1,4) = 1
            lv(1,2,4) = 1
            lv(2,2,4) = 2
            lv(1,3,4) = 1
            lv(2,3,4) = 2
         else if(il == 3) then
            nv = 2
            lv(1,1,2) = 1
            lv(2,1,2) = 2
            lv(1,2,2) = 1
            lv(2,2,2) = 2
            lv(1,3,2) = 1
            lv(2,3,2) = 2
         else if(il.eq.4) then
            nv = 2
            lv(1,1,2) = 1
            lv(2,1,2) = 2
            lv(1,2,2) = 1
            lv(2,2,2) = 2
            lv(1,3,2) = 0
            lv(2,3,2) = 1
         else if(il == 5) then
            nv = 2
            lv(1,1,2) = 0
            lv(2,1,2) = 1
            lv(1,2,2) = 1
            lv(2,2,2) = 2
            lv(1,3,2) = 1
            lv(2,3,2) = 2
         else if(il == 6) then
            nv = 2
            lv(1,1,2) = 1
            lv(2,1,2) = 2
            lv(1,2,2) = 0
            lv(2,2,2) = 1
            lv(1,3,2) = 1
            lv(2,3,2) = 2
         else
            if(printable) then
               write(nfout,*) 'il=',il
               write(nfout,*) '=== stop at sub.nsprmv. (il) ==='
            end if
            stop 'error:sub.nsprmv'
         end if
      end subroutine nsprmv
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrduc(in,id)
!
!#12  in-out put :    in : numerator
!#12                  id : denominator
!
!#21  to reduce a fraction in/id
!
!#31  1989.12.28.:  n. hamada, a. yanase and k. terakura
!
!c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrduc(in,id)
         implicit none
         integer, intent(inout) :: in, id
         integer :: m

         if(in == 0) then
            id = 1
         else
            call nsgcm2(in,id,m)
            in=in/m
            id=id/m
         end if
      end subroutine nsrduc
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrmxc(euler,r)
!
!#12  input :  euler(3) : euler angles (alpha,beta, gamma)
!#12  output     r(3,3) : rotation matrix in real
!#13  noexternal
!
!#21  to get rotation matrix for cubic operations
!
!#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrmxc(euler,r)
use m_Const_Parameters, only :  DP
         implicit none
         real(DP), intent(in), dimension(3) :: euler
         real(DP), intent(out), dimension(3,3) :: r
         real(DP) :: ca, cb, cc,sa, sb, sc
!
!                        euler(1) :  alpha
!                        euler(2) :  beta
!                        euler(3) :  gamma
!
!                  point : x         u     x
!                          y  ---->  v = r y
!                          z    r    w     z
!
         ca = dcos(euler(1))
         cb = dcos(euler(2))
         cc = dcos(euler(3))
         sa = dsin(euler(1))
         sb = dsin(euler(2))
         sc = dsin(euler(3))
         r(1,1) = ca*cb*cc - sa*   sc
         r(2,1) = sa*cb*cc + ca*   sc
         r(3,1) =          -    sb*cc
         r(1,2) =-ca*cb*sc - sa*   cc
         r(2,2) =-sa*cb*sc + ca*   cc
         r(3,2) =               sb*sc
         r(1,3) =            ca*sb
         r(2,3) =            sa*sb
         r(3,3) =               cb
      end subroutine nsrmxc
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrmxh(euler,r)
!
!#12  input :  euler(3) : euler angles (alpha,beta, gamma)
!#12  output     r(3,3) : rotation matrix for hexagonal
!#13  external : sub.nsmlt1
!
!#21  to get a rotation matrix for a hexagonal operation
!
!#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrmxh(euler,r)
         use m_Const_Parameters, only :  DP
         implicit none
         real(DP), intent(in), dimension(3)  :: euler
         real(DP), intent(out), dimension(3,3) :: r
         real(DP) :: three, ca, cb, cc, sa, sb, sc, p2, p3
         real(DP), dimension(3,3,3) :: rr
         real(DP), dimension(3,3) :: rrr
!
!                        euler(1) :  alpha
!                        euler(2) :  beta
!                        euler(3) :  gamma
!
!                  point : x         u     x
!                          y  ---->  v = r y
!                          z    r    w     z
!
         three = 3.0d0
         ca = dcos(euler(1))
         cb = dcos(euler(2))
         cc = dcos(euler(3))
         sa = dsin(euler(1))
         sb = dsin(euler(2))
         sc = dsin(euler(3))
         p2 = 2.0d0
         p3 = dsqrt(three)
         rr(1,1,1) =  ca + sa/p3
         rr(2,1,1) =  p2*sa/p3
         rr(3,1,1) =  0.0d0
         rr(1,2,1) = -rr(2,1,1)
         rr(2,2,1) =  ca -sa/p3
         rr(3,2,1) =  0.0d0
         rr(1,3,1) =  0.0d0
         rr(2,3,1) =  0.0d0
         rr(3,3,1) =  1.0d0

         rr(1,1,2) =  cb
         rr(2,1,2) =  0.0d0
         rr(3,1,2) = -sb
         rr(1,2,2) = (1.0d0-cb)/p2
         rr(2,2,2) =  1.0d0
         rr(3,2,2) =  sb/p2
         rr(1,3,2) =  sb
         rr(2,3,2) =  0.0d0
         rr(3,3,2) =  cb

         rr(1,1,3) =  cc+sc/p3
         rr(2,1,3) =  p2*sc/p3
         rr(3,1,3) =  0.0d0
         rr(1,2,3) = -p2*sc/p3
         rr(2,2,3) =  cc-sc/p3
         rr(3,2,3) =  0.0d0
         rr(1,3,3) =  0.0d0
         rr(2,3,3) =  0.0d0
         rr(3,3,3) =  1.0d0

         call nsmlt1(3,3,rr(1,1,1),rr(1,1,2),rrr)
         call nsmlt1(3,3,rrr,rr(1,1,3),r)
      end subroutine nsrmxh
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrot1(ja,ir,jb)
!
!12  input   : ja(2,3)  : a vector
!#12            ir(3,3)  : rotation matrix
!#12  output  : jb(2,3)  : a vector rotated
!!c#13  external: nsgcm3
!
!#21  to get a rotated vector
!#21  (jb(1,1)/jb(2,1), jb(1,2)/jb(2,2), jb(1,3)/jb(2,3))
!
!#31  1989.12.28.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrot1(ja,ir,jb)
         implicit none
         integer, intent(in), dimension(2,3) :: ja
         integer, intent(in), dimension(3,3) :: ir
         integer, intent(out), dimension(2,3) :: jb
         integer ::  m, mm, id, i1, i2, i3

         call nsgcm3(ja(2,1),ja(2,2),ja(2,3),m)
         mm = m*m
         id = ja(2,1)*ja(2,2)*ja(2,3)/mm
         i1 = ja(1,1)*(id/ja(2,1))
         i2 = ja(1,2)*(id/ja(2,2))
         i3 = ja(1,3)*(id/ja(2,3))

!    ja=(i1,i2,i3)/id

         jb(1,1) = ir(1,1)*i1+ir(1,2)*i2+ir(1,3)*i3
         jb(1,2) = ir(2,1)*i1+ir(2,2)*i2+ir(2,3)*i3
         jb(1,3) = ir(3,1)*i1+ir(3,2)*i2+ir(3,3)*i3
         jb(2,1) = id
         jb(2,2) = id
         jb(2,3) = id

         call nsrduc(jb(1,1),jb(2,1))
         call nsrduc(jb(1,2),jb(2,2))
         call nsrduc(jb(1,3),jb(2,3))
         jb(1,1)=mod(jb(1,1),jb(2,1))
         jb(1,2)=mod(jb(1,2),jb(2,2))
         jb(1,3)=mod(jb(1,3),jb(2,3))
      end subroutine nsrot1
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nsrota(natm,lmnatm,vatm,iatm21,ratm21)
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : il9, ng9, ig0, jv
!$ use omp_lib
         implicit none
         integer, intent(in) :: natm, lmnatm
         integer, intent(out), dimension(lmnatm,3) :: iatm21
         integer :: i, ig1, j1, j2, jj2
         integer :: ind
         integer :: nthreads
         real(DP), intent(in), dimension(3,lmnatm) :: vatm
         real(DP), intent(out), dimension(3,lmnatm,48) :: ratm21
         real(DP), dimension(3) :: v, r

         nthreads = min(ng9,3)
         do ig1=1, ng9
            do j1=1,natm
               call nsrotr(ig1,vatm(1,j1),v)
               do  j2=1 ,natm
                  jj2 = j2
                  r(1:3) = vatm(1:3,j2)-v(1:3)
                  call nslatr(il9,r,ind)
                  if(ind == 0) then
                     exit
                  end if
               end do
               if(ind /= 0) then
                  stop ' === stop in sub.nsrota. (no atom) ==='
               end if
               iatm21(j1,ig1) = jj2
               ratm21(1:3,j1,ig1) = r(1:3)
            end do
         end do
      end subroutine nsrota
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
!
!#12  noinput
!#12  output: ng0                 : # of group elements
!#12          ieuler(3,24)*2*pai/4: euler angles for cubic
!#12           euler(3,24)        : euler angles (alpha,beta, gamma)
!#12          irot(3,3,48)        : rotation matrix in integer
!#12           rot(3,3,48)        : rotation matrix in real
!#12          ir1234(3,48)        : jones faithful representation
!#12           jones(3,48)   (a2) : jones faithful representation
!#12          schoen(3,48)   (a5) : schoenflies notation
!#13  external : sub.nseulc, nsrmxc
!
!#21  to get euler angles, rotation matrix,
!#21         jones faithfull representation, and
!#21         schoenflies notation
!#21         for point-group operations
!#21         of cubic lattice
!
!#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(out) :: ng0
         integer, intent(out), dimension(3,24) :: ieuler
         integer, intent(out), dimension(3,3,48) :: irot
         integer, intent(out), dimension(3,48) :: ir1234
         integer :: i, j, k, irr
         real(DP), intent(out), dimension(3,24) :: euler
         real(DP), intent(out), dimension(3,3,48) :: rot
         real(DP), dimension(3,3) :: r
         real(DP) :: one
         character(len=2), dimension(3,48) :: jones
         character(len=5), dimension(48) :: schoen
         character(len=1), dimension(3) :: cha0
         character(len=4), dimension(24) :: cmn

         data cmn            /'e   ','c2x ','c2y ','c2z ', &
         &                    'c31+','c32+','c33+','c34+', &
         &                    'c31-','c32-','c33-','c34-', &
         &      'c2a ','c2b ','c2c ','c2d ','c2e ','c2f ', &
         &      'c4x+','c4y+','c4z+','c4x-','c4y-','c4z-'/
         data cha0 /'x','y','z'/
         data one /1.00001d0/

         ng0 = 48
         call nseulc(ieuler,euler)
         do i=1, 24
            call nsrmxc(euler(1,i),r)
            irot(1:3,1:3,i)=int(one*r(1:3,1:3))
            rot(1:3,1:3,i)=r(1:3,1:3)
            do j=1, 3
               do k=1, 3
                  irr = irot(j,k,i)
                  if(irr == -1) then
                     ir1234(j,i) = -k
                  else if(irr == 1) then
                     ir1234(j,i) = k
                  end if
               end do
            end do
         end do
         schoen(1:24) = ' '//cmn(1:24)
         schoen(25:48) = 'i'//cmn(1:24)
         ir1234(1:3,25:48) = -ir1234(1:3,1:24)
         irot(1:3,1:3,25:48) = -irot(1:3,1:3,1:24)
         rot(1:3,1:3,25:48) = -rot(1:3,1:3,1:24)
         do i=1, 48
            do j=1, 3
               k = ir1234(j,i)
               if(k < 0) then
                  jones(j,i) ='-'//cha0(iabs(k))
               else
                  jones(j,i)=' '//cha0(iabs(k))
               end if
            end do
         end do
      end subroutine nsrotc
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
!
!#12  noinput
!#12  output: ng0                 : # of group elements
!#12          ieuler(3,24)*2*pai/4: euler angles for hexagonal
!#12           euler(3,24)        : euler angles (alpha,beta, gamma)
!#12          irot(3,3,48)        : rotation matrix in integer
!#12           rot(3,3,48)        : rotation matrix in real
!#12          ir1234(3,48)        : jones faithful representation
!#12           jones(3,48)   (a2) : jones faithful representation
!#12          schoen(3,48)   (a5) : schoenflies notation
!#13  external : sub.nseulh, nsrmxh, nsjonh, nstrsh, nsmlt1
!
!#21  to get euler angles, rotation matrix,
!#21         jones faithfull representation, and
!#21         schoenflies notation
!#21         for point-group operations
!#21         of hexagonal lattice
!
!c#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(out) :: ng0
         integer, intent(out), dimension(3,24) :: ieuler
         integer, intent(out), dimension(3,3,48) :: irot
         integer, intent(out), dimension(3,48) :: ir1234
         integer :: i, j, k, l, i1, i2
         integer, dimension(3,3) :: iw
         real(DP), intent(out), dimension(3,24) :: euler
         real(DP), intent(out), dimension(3,3,48) :: rot
         real(DP), dimension(3,3) :: r, u, v, wk1, wk2
         real(DP) :: onep
         character(len=2), intent(out), dimension(3,48) :: jones
         character(len=5), intent(out), dimension(48) :: schoen
         character(len=4), dimension(12) :: hmn
         data        hmn /'e   ','c6+ ','c3+ ','c2  ','c3- ','c6- ', &
         &                'c211','c221','c231','c212','c222','c232'/
        data onep /1.00001/

        ng0 = 24
        call nseulh(ieuler,euler)
        do  i=1, 12
           call nsrmxh(euler(1,i),r)
           do k = 1,3
              do j = 1,3
                 rot(j,k,i) = r(j,k)
                 irot(j,k,i) = int(onep*r(j,k))
                 rot(j,k,i+12) = -r(j,k)
                 irot(j,k,i+12) = -int(onep*r(j,k))
              end do
           end do
        end do
        call nsjonh(irot,ir1234,jones)

        schoen(1:12) = ' '//hmn(1:12)
        schoen(13:24) = 'i'//hmn(1:12)

        call nstrsh(u,v)
        do i=1, 24
           wk1(1:3,1:3) = rot(1:3,1:3,i)
           call nsmlt1(3,3,v,wk1,r)
           call nsmlt1(3,3,r,u,wk2)
           rot(1:3,1:3,i+24) = wk2(1:3,1:3)
        end do

        irot(1:3,1:3,25:48)=int(onep*rot(1:3,1:3,25:48))

        call nsjonh(irot(1,1,25),ir1234(1,25),jones(1,25))

        do i=25, 48
           do j = 1, 3
              do k = 1,3
                 iw(j,k)=irot(k,j,i)
              end do
           end do
           irot(1:3,1:3,i) = iw(1:3,1:3)
           rot(1:3,1:3,i) = real(iw(1:3,1:3),DP)
        end do

!     check

        do i=1, 24
            do j=1, 3
               do k=1, 3
                  i1 = 0
                  i2 = 0
                  do l=1, 3
                     i1 = i1+irot(j,l,i)*irot(l,k,i+24)
                     i2 = i2+irot(j,l,i+24)*irot(l,k,i)
                  end do
                  if((j == k .and. (i1 /= 1 .or. i2 /= 1)) .or. &
                   & (j /= k .and. (i1 /= 0 .or. i2 /= 0))) then
                     if(printable) then
                        write(nfout,*) &
                       & ' j,k=',j,k,'   i1,i2=',i1,i2,'  : i=',i
                        write(nfout,*) ' === stop in sub.nsroth. ', &
                       & ' (inverse matrix condition) ==='
                     end if
                     stop 'error in sub.nsroth'
                  end if
               end do
            end do
         end do
      end  subroutine nsroth
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrotk(ig,k1,k2)
         use m_Spg_plus_Tetra_Common_Blocks, only : il9, ig0, irot
         implicit none
         integer, intent(in) :: ig
         integer, intent(in), dimension(4) :: k1
         integer, intent(out), dimension(4) :: k2
         integer :: igg, i, j, k

         if(il9 <= 0) then
            igg = ig0(ig)+24
         else
            igg = ig0(ig)
         end if
         k2(4)=k1(4)
         do i=1, 3
            k = 0
            do j=1, 3
               k = k+irot(i,j,igg)*k1(j)
            end do
            k2(i) = k
         end do
      end subroutine nsrotk
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nsrotl(jf,jpr,il,ng0,ieuler,euler,irot,rot,
!#11             ir1234,jones,schoen                 )
!
!#12  input: jf    : output fileil=<0 for hexagonal, il>0 for cubic.
!#12         jpr   : output control
!#12         il    : lattice typei: il=<0 for hexagonal, il>0 for cubic.
!
!#12  output: ng0                 : # of group elements
!#12          ieuler(3,48)*2*pai/6: euler angles for hexagonal
!#12          ieuler(3,48)*2*pai/4: euler angles for cubic
!#12           euler(3,48)        : euler angles (alpha,beta, gamma)
!#12          irot(3,3,48)        : rotation matrix in integer
!#12           rot(3,3,48)        : rotation matrix in real
!#12           jones(3,48)   (a2) : jones faithful representaion
!#12          schoen(3,48)   (a5) : schoenflies notation
!#13  external : sub.nsrotc, nsroth
!
!#21  to get euler angles, rotation matrix,
!#21         jones faithfull representation, and
!#21         schoenflies notation
!#21         for point-group operations
!#21         of hexagonal and cubic lattices
!
!#31  1986.06.07.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrotl(il,jpr,jf,ng0,ieuler,euler,irot,rot, &
     &                  ir1234,jones,schoen)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(in) :: il, jpr, jf
         integer, intent(out) :: ng0
         integer, intent(out), dimension(3,24) :: ieuler
         integer, intent(out), dimension(3,3,48) :: irot
         integer, intent(out), dimension(3,48) :: ir1234
         integer :: no, n0, i, j, k, ii
         real(DP), intent(out), dimension(3,24) :: euler
         real(DP), intent(out), dimension(3,3,48) :: rot
         character(len=2), intent(out), dimension(3,48) :: jones
         character(len=5), intent(out), dimension(48) :: schoen

         if(il <= 0) then
            call nsroth(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
         else
            call nsrotc(ng0,ieuler,euler,irot,rot,ir1234,jones,schoen)
         end if

         no = ng0
         if(jpr >= 2) then
            write(jf,*) '     '
            if(il <= 0) then
               write(jf,*) &
              & ' table of operation code for the point group', &
              &   ' 6/mmm(d6h) '
            else
               write(jf,*) &
              & ' table of operation code for the point group',&
              & ' m3m(oh) '
            end if
            if(il <= 0) write(jf,*) &
           & ' (for real-space coordinate, w=x-y)'
            write(jf,100)
            do i=1, no/2
               write(jf,120) i,schoen(i),(jones(j,i),j=1,3), &
              & (ir1234(j,i),j=1,3),(euler(j,i),j=1,3)
            end do
            do i=no/2+1, no
               write(jf,140) i,schoen(i),(jones(j,i),j=1,3), &
             & (ir1234(j,i),j=1,3)
            end do
            if(il <= 0) then
               write(jf,*) '    '
               write(jf,*) ' (for reciprocal-space coordinate, w=x+y) '
               do i=25, 48
                  write(jf,160) i,          (jones(j,i),j=1,3), &
                & (ir1234(j,i),j=1,3)
               end do
            end if
         end if
  100    format(1h ,' no.',' schoenflies   jones       ir1234    ',&
        & '   euler-angle(alpha,beta,gamma)')
  120    format(1h ,'(',i2,')',5x,a5,' (',3(1x,a2),' )  (',3(1x,i2),' )', &
        &       5x,'(',3f8.3,' )')
  140    format(1h ,'(',i2,')',5x,a5,' (',3(1x,a2),' )  (',3(1x,i2),' )')
  160    format(1h ,'(',i2,')',5x,5x,' (',3(1x,a2),' )  (',3(1x,i2),' )')

         if(jpr >= 3) then
            write(jf,*) '   '
            write(jf,*) ' matrix representation of operation '
            if(il <= 0) write(jf,*) ' (for real-space coordinate)'
            do  ii=1, no, 6
               write(jf,*) '   '
               write(jf,200) (i,(irot(1,k,i),k=1,3),i=ii,ii+5)
               do  j=2,3
                  write(jf,220) (  (irot(j,k,i),k=1,3),i=ii,ii+5)
               end do
            end do
            if(il <= 0) then
               write(jf,*) '   '
               write(jf,*) ' matrix representation of operation '
               write(jf,*) ' (for reciprocal-space coordinate)'
               do ii=25, 48, 6
                  write(jf,*) '   '
                  write(jf,200) (i,(irot(1,k,i),k=1,3),i=ii,ii+5)
                  do j=2, 3
                     write(jf,220) (  (irot(j,k,i),k=1,3),i=ii,ii+5)
                  end do
               end do
            end if
         end if
  200    format(1h ,6('(',i2,')',3i2,2x))
  220    format(1h ,6(4x        ,3i2,2x))
      end  subroutine nsrotl
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsrotr(ig,r1,r2)
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : rot,ig0, jv
         implicit none
         integer, intent(in) :: ig
         real(DP), intent(in), dimension(3) :: r1
         real(DP), intent(out), dimension(3) :: r2

         integer :: i, j
         real(DP) :: a

         do i=1, 3
            a = 0.0d0
            do j=1, 3
               a = a+rot(i,j,ig0(ig))*r1(j)
            end do
            r2(i) = a+real(jv(1,i,ig),kind=DP)/jv(2,i,ig)
         end do
      end  subroutine nsrotr
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nssdjg(jf,jpr,nv,lv,ng,jg)
!
!#12  input     :        jf : output file
!#12                    jpr : print control
!#12                     nv : # of lv vectors
!#12              lv(2,3,4) : primitive translation vector
!#12                     ng : oder of group
!#12  in-output : jg(2,3,48): translation vector
!#12  noexternal
!
!#21  to imodify nonprimitive translation vectors
!#21  (positive, smaller elements)
!
!#31  1990.01.06.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nssdjg(jf,jpr,nv,lv,ng,jg)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(in) :: jf, jpr, nv, ng
         integer, intent(in), dimension(2,3,4) :: lv
         integer, intent(out), dimension(2,3,48) :: jg
         integer :: iflag
         integer :: na, m, ia, i, j, k
         integer, dimension(2,3,48) :: ja
         integer, dimension(2,3) :: jb, jc

         do  j=1, 3
            do  k=1, ng
              jg(1,j,k) = mod(jg(1,j,k),jg(2,j,k))
             if(jg(1,j,k)< 0) jg(1,j,k) = jg(1,j,k)+jg(2,j,k)
            end do
         end do

         ja(1,1:3,1) = 0
         ja(2,1:3,1) = 1

         na = 1

         do m=1,ng
            do ia=1, na
               do  k=1, nv
                  iflag = 0
                  call nssum1(1,jg(1,1,m),lv(1,1,k),jb)
                  if(jb(1,1) == ja(1,1,ia) .and. &
                 &   jb(2,1) == ja(2,1,ia) .and. &
                 &   jb(1,2) == ja(1,2,ia) .and. &
                 &   jb(2,2) == ja(2,2,ia) .and. &
                 &   jb(1,3) == ja(1,3,ia) .and. &
                 &   jb(2,3) == ja(2,3,ia)) then
                     jg(1:2,1:3,m) = ja(1:2,1:3,ia)
                     iflag = 1
                     exit
                 end if
               end do
               if(iflag == 1) exit
            end do
            if(iflag /= 1 ) then
!       registration
               na = na+1
               jb(1:2,1:3)=jg(1:2,1:3,m)
               do k=2, nv
                  call nssum1(1,jg(1,1,m),lv(1,1,k),jc)
                  if(jb(1,1) > jc(1,1) .or. &
                 &   jb(1,2) > jc(1,2) .or. &
                 &   jb(1,3) > jc(1,3)) then
                     jb(1:2,1:3) = jc(1:2,1:3)
                  end if
               end do
               jg(1:2,1:3,m) =jb(1:2,1:3)
               ja(1:2,1:3,na)=jb(1:2,1:3)
            else
               cycle
            end if
         end do

         if(jpr > 1) then
            write(jf,*) ' '
            write(jf,*) '= nonprimitive translation vectors =   na=',na
            do k=1, na
               write(jf,100) ((ja(i,j,k),i=1,2),j=1,3)
            end do
         end if
  100   format(1h ,'( ',3(i2,'/',i2,2x),')')
      end subroutine nssdjg
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nssum1(is,ja,jb,jc)
!
!#12  input   : is       : add(1) or subtract(-1)
!#12            ja(2,3)  : a vector
!#12            jb(2,3)  : a vector
!#12  output  : jc(2,3)  : sum of ja and jb, and take modulus
!#13  external: nsgcm3
!
!#21  to get a summation of vectors
!#21  (ja(1,1)/ja(2,1), ja(1,2)/ja(2,2), ja(1,3)/ja(2,3))  etc.
!
!#31  1989.12.28.:  n. hamada, a. yanase, and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nssum1(is,ja,jb,jc)
         implicit none
         integer, intent(in) :: is
         integer, intent(in), dimension(2,3) :: ja, jb
         integer, intent(out), dimension(2,3) :: jc
         integer :: m, n, id, k, id1, id2, i1, i2, i3, j1, j2, j3

         call nsgcm3(ja(2,1),ja(2,2),ja(2,3),m)
         call nsgcm3(jb(2,1),jb(2,2),jb(2,3),n)
         id1 = ja(2,1)*ja(2,2)*ja(2,3)/m*m
         id2 = jb(2,1)*jb(2,2)*jb(2,3)/n*n
         call nsgcm2(id1,id2,k)
         id = id1*id2/k
         i1 = ja(1,1)*(id/ja(2,1))
         i2 = ja(1,2)*(id/ja(2,2))
         i3 = ja(1,3)*(id/ja(2,3))
         j1 = jb(1,1)*(id/jb(2,1))
         j2 = jb(1,2)*(id/jb(2,2))
         j3 = jb(1,3)*(id/jb(2,3))
         jc(1,1) = i1+j1*is
         jc(1,2) = i2+j2*is
         jc(1,3) = i3+j3*is
         jc(2,1) = id
         jc(2,2) = id
         jc(2,3) = id
         call nsrduc(jc(1,1),jc(2,1))
         call nsrduc(jc(1,2),jc(2,2))
         call nsrduc(jc(1,3),jc(2,3))
         jc(1,1)=mod(jc(1,1),jc(2,1))
         jc(1,2)=mod(jc(1,2),jc(2,2))
         jc(1,3)=mod(jc(1,3),jc(2,3))
      end subroutine nssum1
!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!#11  sub.nstrsh(u,v)
!
!#12  noinput
!#12  output: u(3,3)
!#12          v(3,3)
!#13  noexternal
!
!#21  to get matrix u and v
!
!#31  1986.06.07.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstrsh(u,v)
         use m_Const_Parameters, only :  DP
         implicit none
         real(DP), dimension(3,3) :: u, v
         real(DP) :: three, a, a1, a2


         u(1:3,1:3) = 0.0d0
         v(1:3,1:3) = 0.0d0
         u(3,3) = 1.0d0
         v(3,3) = 1.0d0

         a = dsqrt(3.0d0)
         a1 = 1.0d0/a
         a2 = 2.0d0/a

         u(1,1) = a2
         u(1,2) = a1
         u(2,1) = a1
         u(2,2) = a2
         v(1,1) = a2
         v(1,2) = -a1
         v(2,1) = -a1
         v(2,2) = a2
      end subroutine nstrsh
!c --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nstt0i(e,e1,e2,e3,e4,dos,dosin)
!
!  ** dos and integrated dos 'coefficients' at **  c
!  ** tetrahedron corners at one energy        **  c
!  ** according to lambin and vigneron,        **  c
!  ** phys. rev. b29, 3430 (1984)              **  c
! 2021.12.20: modified by T> Hamada for perfromnce improvement
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : ncounter, ncounter1,ncounter2, ncounter3, factor_1_12
         implicit none
         real(DP), intent(in) :: e, e1, e2, e3, e4
         real(DP), intent(out), dimension(4) :: dos, dosin
         real(DP) :: d1, d2, d3, d4, d21, d31, d41, d32, d42, d43
         real(DP) :: d1d1, d2d2, d3d3, d4d4, d21d21,d31d31,d32d32,d41d41,d42d42, &
                & d31d32, d31d41, d32d41, d31d42, d32d42, d41d42, &
                & d1_3, d3_3, d4_3, d21_3
         real(DP) :: x, y, xx, zz, yy, x1, x2, x3, x4, x5, x6, &
        & y1, y2, y3, y4


         d21 = e2-e1
         d31 = e3-e1
         d41 = e4-e1
         d32 = e3-e2
         d42 = e4-e2
         d43 = e4-e3
         d1 = e-e1
         d4 = e4-e

         if(e <= e2) then
            ncounter1 = ncounter1+1
            d2=e2-e
            d3=e3-e
            yy=d41*d31*d21
            x=d2/d21+d3/d31+d4/d41
            d1d1 = d1**2
            y=(d1d1)/yy
            dos(1)=x*y
            dosin(1)=0.25*d1*y*(x+1.0d0)
            xx=d1d1*d1
            x=xx/(d21*yy)
            dos(2)=x
            dosin(2)=0.25d0*d1*x
            x=xx/(d31*yy)
            dos(3)=x
            dosin(3)=0.25d0*d1*x
            x=xx/(d41*yy)
            dos(4)=x
            dosin(4)=0.25d0*d1*x
         else if(e >= e3) then
            ncounter3 = ncounter3+1
            d2=e-e2
            d3=e-e3
            xx=d4**3
            yy=d41*d42*d43
            x=xx/(d41*yy)
            dos(1)=x
            dosin(1)=0.25d0*(1.0d0-d4*x)
            x=xx/(d42*yy)
            dos(2)=x
            dosin(2)=0.25d0*(1.0d0-d4*x)
            x=xx/(d43*yy)
            dos(3)=x
            dosin(3)=0.25d0*(1.0d0-d4*x)
            x=d3/d43+d2/d42+d1/d41
            y=(d4*d4)/yy
            dos(4)=x*y
            dosin(4)=0.25*(1.0d0-d4*y*(x+1.0d0))
         else
            ncounter2 = ncounter2+1
            d2 = e-e2
            d3 = e3-e
            d1d1 = d1**2
            d2d2 = d2**2
            d3d3 = d3**2
            d4d4 = d4**2
            d21d21 = d21**2
            d31d31 = d31**2
            d32d32 = d32**2
            d41d41 = d41**2
            d42d42 = d42**2
            d31d32 = d31*d32
            d31d41 = d31*d41
            d31d42 = d31*d42
            d32d41 = d32*d41
            d32d42 = d32*d42
            d41d42 = d41*d42
            d1_3 = d1*3.0d0
            d3_3 = d3*3.0d0
            d4_3 = d4*3.0d0
            d21_3 = d21*3.0d0

            y = 1.0d0/(d31d42)+1.0d0/(d32d41)
            x1 = (d3d3)/(d31d31*d32)*(d2/d42+d1/d41)
            x2 = (d4d4)/(d41d41*d42)*(d2/d32+d1/d31)
            x3 = (d3*d4*d1)/(d31d41)*y
            dos(1) = 0.5d0*(x1+x2+x3)
            x1 = d2d2*(d32*d2+d3_3*(d32+d3))*factor_1_12
            y1 = x1
            x1 = x1/(d31d31*d32d42)
            x2 = d2*(d2d2*(d31+d21_3)+d3_3*(d2*d3+d32*(d21_3+d1)))
            x2 = x2*factor_1_12
            y2 = x2
            x2 = x2/(d31d31*d32d41)
            x3 = d2d2*(d42*d2+d4_3*(d42+d4))*factor_1_12
            y3 = x3
            x3 = x3/(d41d41*d32d42)
            x4 = d2*(d2d2*(d41+d21_3)+d4_3*(d2*d4+d42*(d21_3+d1)))
            x4 = x4*factor_1_12
            y4 = x4
            x4 = x4/(d41d41*d31d42)
            x5 = 0.5d0*d2*d3*d4*(d1+d21)
            x5 = x5+d2d2*((d21_3-d21)*(d3+d42)+(d1+d21)*((d3_3-d3)+d4+d42))*factor_1_12
            x5 = x5*y/(d31d41)
            x6 = 0.25d0*d21d21*(d42/d41+d32/d31+1.0d0)/(d31d41)
            dosin(1) = 0.5d0*(x1+x2+x3+x4+x5)+x6

            x1 = (d3d3)/(d32d32*d31)*(d2/d42+d1/d41)
            x2 = (d4d4)/(d42d42*d41)*(d2/d32+d1/d31)
            x3 = (d3*d4*d2)/(d32d42)*y
            dos(2) = 0.5d0*(x1+x2+x3)
            x1 = y1/(d32d32*d31d42)
            x2 = y2/(d32d32*d31d41)
            x3 = y3/(d42d42*d32d41)
            x4 = y4/(d42d42*d31d41)
            x5 = d2d2*(d3*(d42+d4_3)+d32*(d42+d4))*factor_1_12
            x5 = x5*y/(d32d42)
            x6 = 0.25d0*d21/d31*d21/d41
            dosin(2) = 0.5d0*(x1+x2+x3+x4+x5)+x6
            x1 = (d2d2)/(d32d32*d42)*(d3/d31+d4/d41)
            x2 = (d1d1)/(d31d31*d41)*(d3/d32+d4/d42)
            x3 = (d1*d2*d3)/(d31d32)*y

            dos(3) = 0.5d0*(x1+x2+x3)
            x1 = d2d2*d2*(d3_3+d32)*factor_1_12
            y1 = x1
            x1 = x1/(d32d32*d31d42)
            x2 = d2d2*d2*(d4_3+d42)*factor_1_12
            y2 = x2
            x2 = x2/(d32d32*d41d42)
            x3 = d2*(d2*d31*(d2+d21_3)+d3_3*(d2d2+d21_3*d1)+ &
         &  3.0d0*d21d21*d32)*factor_1_12
            y3 = x3
            x3 = x3/(d31d31*d32d41)
            x4= d2*(d2*d41*(d2+d21_3)+d4_3*(d2d2+d21_3*d1)+3.0d0* &
         &  d21d21*d42)*factor_1_12
            y4 = x4
            x4 = x4/(d31d31*d41d42)
            x5 = (d2d2)*(d3*(d21+d1_3)+d32*(d21+d1))*factor_1_12
            x5 = x5/(d31d32)*y
            x6 = 0.25d0*(d21/d31*d21/d31*d21/d41)
            dosin(3) = 0.5d0*(x1+x2+x3+x4+x5)+x6

            x1 = (d2d2)/(d42d42*d32)*(d3/d31+d4/d41)
            x2 = (d1d1)/(d41d41*d31)*(d3/d32+d4/d42)
            x3 = (d1*d2*d4)/(d41d42)*y
            dos(4) = 0.5d0*(x1+x2+x3)
            x1 = y1/(d42d42*d31d32)
            x2 = y2/(d42d42*d32d41)
            x3 = y3/(d41d41*d31d32)
            x4 = y4/(d41d41*d31d42)
            x5 = (d2d2)*(d4*(d21+d1_3)+d42*(d21+d1))*factor_1_12
            x5 = x5/(d41d42)*y
            x6 = 0.25d0*(d21/d31*d21/d41*d21/d41)
            dosin(4) = 0.5d0*(x1+x2+x3+x4+x5)+x6
         end if
      end subroutine nstt0i
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!      subroutine nsttim(e1,e2,e3,e4,dos,dosin)
!         use m_Const_Parameters, only :  DP
!         implicit none
!         real(DP), intent(in) ::  e1,e2,e3,e4
!         real(DP), intent(in), dimension(4) :: dos
!         real(DP), intent(out), dimension(4) :: dosin
!         real(DP) :: tdos, esum
!         real(DP), dimension(4) :: e, d
!
!         e(1) = e1
!         e(2) = e2
!         e(3) = e3
!         e(4) = e4
!         esum = e1+e2+e3+e4
!!$$$      tdos = 0.025d0*sum(dos(1:4))
!         tdos = 0.025d0*(dos(1)+dos(2)+dos(3)+dos(4))
!$$$         dosin(i) = dosin(i) + tdos*sum(e(1:4)-e(i))
!         dosin = 0.0d0
!         dosin(1) = dosin(1) + tdos*(esum-4.d0*e(1))
!         dosin(2) = dosin(2) + tdos*(esum-4.d0*e(2))
!         dosin(3) = dosin(3) + tdos*(esum-4.d0*e(3))
!jj         dosin(4) = dosin(4) + tdos*(esum-4.d0*e(4))
!      end subroutine nsttim
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
!
      subroutine nstt1i(e,nxx,nyy,nzz,ea,cdos,cind)
!
!     ***  'densities of states and integrated densities of ***
!     ***  coefficients' are obtained at k-mesh points      ***
!     ***  for a single energy                              ***
!     ***  by linear interpolation in tetrahedorons         ***
!     ***
!     ***
!     ***  npx  number of mesh points in b.z.             ***
!     ***         in x-direction
!     ***  npy  number of mesh points in b.z.
!     ***         in y-direction
!     ***  npz  number of mesh points in b.z.             ***
!     ***         in z-direction
!     ***  cdos density of states 'coefficient'
!     ***  cind number  of states 'coefficient'
!     ***  flag If flag is .true., quadric correction
!     ***       will be performed
!     ***  ref.)  ph. lambin and j. p. vigneron,
!                 phys. rev. b29 (1984) 3430.
! 2o20.12.20:modified by T. Hamada
         use m_Const_Parameters, only :  DP
         use m_Control_Parameters, only : width_tetra
         use m_Spg_plus_Tetra_Common_Blocks, only : ncounter,ncounter1,ncounter2,ncounter3,  &
        & npx, npy, npz, np, ncub, ntet, rntet_1, tetra_eps, ip0_index1,ip0_index2, &
        & tet_sub_called, nstt1i_flag
         use m_Parallelization, only : npes,mype, MPI_CommGroup, &
       & mpi_double_precision, mpi_integer, istatus, ierr
         use m_Timing, only : tstatc0_begin, tstatc0_end
!$ use omp_lib
         implicit none
!         include 'mpif.h'
         integer, intent(in) :: nxx, nyy, nzz
         integer :: ni, ip, ip0, iq, it, icub, &
                 & i, ix, iy, iz, kx, ky, kz
         integer, dimension(4) :: iet, ieb
         integer, dimension(8) :: iec
         integer, dimension(2,2,2) :: iecub
         integer, dimension(6,2) :: iqmat
         real(DP), intent(in) :: e
         real(DP), intent(in), dimension((nxx+1)*(nyy+1)*(nzz+1)) :: ea
         real(DP), intent(out), dimension((nxx+1)*(nyy+1)*(nzz+1)) :: cdos, cind
         real(DP) :: eps, emax0, emin0, emax, emin, e1, e2, e3, e4, tdos, esum
         real(DP), dimension(4) :: et, eb, dos, dosin
         real(DP), dimension(8) :: ec
         real(DP), dimension(2,2,2) :: ecub
         logical :: flag

         equivalence(ec(1),ecub(1,1,1))
         equivalence(iec(1),iecub(1,1,1))

         data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/
         integer :: id_sname = -1
!  definition of eps  <- must be consistent with <nstts1>
!        eps = 1.0d-4
         eps = tetra_eps
!         eps = width_tetra*0.5d0

         call tstatc0_begin('nstt1i ',id_sname)
         if( .not. tet_sub_called) then
            npx = nxx+1
            npy = nyy+1
            npz = nzz+1
            np = npx*npy*npz
            ncub = nxx*nyy*nzz
            ntet = 6*ncub
            rntet_1 = 1.0d0/real(ntet,kind=DP)

            allocate(ip0_index1(ncub))
            icub = 0
            do iz=0, nzz-1
               do iy=0, nyy-1
                  do ix=0, nxx-1
                     icub = icub+1
                     ip0_index1(icub)  = npx*(npy*iz+iy)+ix
                  end do
               end do
            end do

            do  kz=1, 2
               do ky=1, 2
                  do  kx=1, 2
                     ip0_index2(kx,ky,kz) = npx*(npy*(kz-1)+ky-1)+kx
                  end do
               end do
            end do

            if(npes > 1) then
               call mpi_bcast(npx,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(npy,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(npz,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(np,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ncub,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ntet,1,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(rntet_1,1,mpi_double_precision,0,MPI_CommGroup,ierr)
               call mpi_bcast(ip0_index1,ncub,mpi_integer,0,MPI_CommGroup,ierr)
               call mpi_bcast(ip0_index2,8,mpi_integer, 0,MPI_CommGroup,ierr)
            end if
             tet_sub_called = .true.
         end if

         cdos(1:np) = 0.0d0
         cind(1:np) = 0.0d0

!        emax0 =-1.0d30
         emin0 = 1.0d30
         do ip=1,np
!           if(ea(ip) > emax0) emax0 = ea(ip)
            if(ea(ip) < emin0) emin0 = ea(ip)
         end do
!                                   ============ if 1 ==
         if(e > emin0) then
!     ***  integration over b.z. starts    ***
!
!     ***       sampling over cubes        ***
            icub = 0
            do iz=0, nzz-1
               do iy=0, nyy-1
                  do ix=0, nxx-1
                     icub = icub+1
!     ***  energies at cube corners  ***
                     emax =-1.0d30
                     emin = 1.0d30
                     do kz=1, 2
                        do ky=1, 2
                           do kx=1, 2
                              ip0 = ip0_index1(icub) + &
                               &    ip0_index2(kx,ky,kz)
                              ecub(kx,ky,kz) = ea(ip0)
                              iecub(kx,ky,kz) = ip0
                              if(ea(ip0) > emax) emax=ea(ip0)
                              if(ea(ip0) < emin) emin=ea(ip0)
                           end do
                        end do
                     end do
!**    if(e.ge.emax+2*eps) then
                     if(e >= emax) then
                        cind(iec(1)) = cind(iec(1))+1.5d0
                        cind(iec(8)) = cind(iec(8))+1.5d0
                        cind(iec(2:7)) = cind(iec(2:7))+0.5d0
                     else if(e> (emin-2.0d0*eps)) then
!     ***      six tetrahedrons      ***
!*
!     *** sampling over tetrahedrons ***
!*
                        et(1) = ec(1)
                        et(4) = ec(8)
                        iet(1) = iec(1)
                        iet(4) = iec(8)
                        do  it=1, 6
                           iq = iqmat(it,1)
                           et(2) = ec(iq)
                           iet(2) = iec(iq)
                           iq = iqmat(it,2)
                           et(3) = ec(iq)
                           iet(3) = iec(iq)
                           eb(1:4) = et(1:4)
                           ieb(1:4) = iet(1:4)
!     ***  eb(1).le.eb(2).le.eb(3).le.eb(4)  ***

                           call nsttod(eb,ieb)

                           e1 = eb(1)
                           e2 = eb(2)
                           e3 = eb(3)
                           e4 = eb(4)
!                                   ============ if 3 ==
                           call nstts1(e1,e2,e3,e4)

                           if(e<=e1) cycle
                           if(e >= e4) then
                               cind(ieb(1:4))=cind(ieb(1:4))+0.25d0
!!*           call nstts1(e1,e2,e3,e4)
                            else
                               call nstt0i(e,e1,e2,e3,e4,dos,dosin)
                               if(nstt1i_flag) call nsttim ! contained here
                               cdos(ieb(1:4))=cdos(ieb(1:4))+dos(1:4)
                               cind(ieb(1:4))=cind(ieb(1:4))+dosin(1:4)
                           end if
                        end do
                     end if
!                                   ============ if 2 ==
                  end do
               end do
            end do
!$OMP simd
            do i = 1, np
               cdos(i) = cdos(i)*rntet_1
               cind(i) = cind(i)*rntet_1
            end do
!$OMP end simd
         end if

         call tstatc0_end(id_sname)
         contains
          subroutine nsttim
! 2021.12.20 : nsttim as a sub origram   T. Hamada
             esum = e1+e2+e3+e4
             tdos = 0.025d0*(dos(1)+dos(2)+dos(3)+dos(4))
             dosin(1) = dosin(1) + tdos*(esum-4.d0*e1)
             dosin(2) = dosin(2) + tdos*(esum-4.d0*e2)
             dosin(3) = dosin(3) + tdos*(esum-4.d0*e3)
             dosin(4) = dosin(4) + tdos*(esum-4.d0*e4)
          end subroutine nsttim
      end  subroutine nstt1i
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstt2i(idim,e,nx,ny,nz,np2,eig2,ip20,np0,ea,cd,cs, &
     &                  cdos,cind)
! 2021.12.20 : modified by T. Hamada
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: idim, nx, ny, nz, np0, np2
         integer, intent(in), dimension(np0) :: ip20
         integer :: i
         real(DP), intent(in) :: e
         real(DP), intent(in), dimension(np2) :: eig2
         real(DP), intent(out), dimension(np0) :: ea
         real(DP), intent(out), dimension(np2) :: cdos, cind
         real(DP), intent(out), dimension(np0) :: cd, cs

#ifdef __TIMER_SUB__
         call timer_sta(714)
#endif
#ifdef __TIMER_DO__
         call timer_sta(819)
#endif
         ea(1:np0)=eig2(ip20(1:np0))
#ifdef __TIMER_DO__
         call timer_end(819)
#endif
         call nstt1i(e,nx,ny,nz,ea,cd,cs)

         cdos(1:np2) = 0.0d0
         cind(1:np2) = 0.0d0
#ifdef __TIMER_DO__
         call timer_sta(820)
#endif
         do i= 1, np0
            cdos(ip20(i)) = cdos(ip20(i))+cd(i)
            cind(ip20(i)) = cind(ip20(i))+cs(i)
         end do
!         do i = 1, np0
!            cdos(ip20(i)) = cdos(ip20(i))+cd(i)
!            cind(ip20(i)) = cind(ip20(i))+cs(i)
!         end do
#ifdef __TIMER_DO__
         call timer_end(820)
#endif
#ifdef __TIMER_SUB__
         call timer_end(714)
#endif
      end  subroutine nstt2i
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstt3i(idim,ne,e,nxx,nyy,nzz, &
     &                  np2,lmnp2e,neig,eeig,ip20,np0,eawk,cdwk,cswk, &
     &                  lmnp2c,lmneig,cdos,cind,eps )
! modified by T. Hamada
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : ncounter,ncounter2,ncounter3, &
         & nstt3i_called, nstt1i_flag
         use m_Parallelization, only : npes,mype, MPI_CommGroup, &
       & mpi_logical, istatus, ierr
!$ use omp_lib

         use m_Timing, only : tstatc0_begin, tstatc0_end
         implicit none
!         include 'mpif.h'
         integer, intent(in) :: idim, ne, nxx, nyy, nzz, np2, lmnp2e, &
         & neig, np0, lmnp2c, lmneig
         integer, intent(in), dimension(np0) :: ip20
         integer :: iloop, ie, ieig, k2, i, n
         real(DP), intent(in) :: eps
         real(DP), intent(in), dimension(0:ne) :: e
         real(DP), intent(in), dimension(lmnp2e,neig) :: eeig
         real(DP), intent(out), dimension(lmnp2c,lmneig,0:ne) :: cdos, cind
         real(DP) :: c1, c2, rn_1
         real(DP), dimension(np0) :: eawk, cdwk, cswk
         integer :: id_sname = -1
!$$$      real*8 etime,wct_now,wct_start
!
!$$$      etime = 0.0d0
!$$$      call gettod(wct_start)
#ifdef __TIMER_SUB__
         call timer_sta(713)
#endif
         if(.not.nstt3i_called) then
            if(idim == 3) then
               nstt1i_flag = .false.
            else
               nstt1i_flag = .true.
            end if
            call mpi_bcast(nstt1i_flag,1,mpi_logical,0,MPI_CommGroup,ierr)
            nstt3i_called = .true.
         end if
!        eps=dfloat(10)**(-5)

         call tstatc0_begin('nstt3i ',id_sname)
         iloop = 0
         ncounter = 0
#ifdef __TIMER_DO__
         call timer_sta(816)
#endif
!$OMP parallel private(ieig,n, rn_1,c1, c2)
!$OMP do
         do  ie=0, ne
            do ieig=1, neig
               call nstt2i(idim,e(ie),nxx,nyy,nzz,np2,eeig(1,ieig),ip20,np0, &
              & eawk,cdwk,cswk,cdos(1,ieig,ie),cind(1,ieig,ie))
               iloop = iloop + (nxx*nyy*nzz*6)
            end do
         end do
!$OMP end do
#ifdef __TIMER_DO__
         call timer_end(816)
#endif
!$$$      write(6,'(" !dos iloop = ",i12)') iloop
!$$$      write(6,'(" !dos ncounter1 = ",i12)') ncounter1
!$$$      write(6,'(" !dos ncounter2 = ",i12)') ncounter2
!$$$      write(6,'(" !dos ncounter3 = ",i12)') ncounter3
!$$$      goto 31

!     take care of a weight on a degenerate state
!$OMP do
         do k2=1, np2
            ieig = 1
   40       continue
            n = 1
        !!$do 42 i=1,20
#ifdef __TIMER_DO__
            call timer_sta(817)
#endif
            do i=1, neig
               if(ieig+i > neig) exit
               if(dabs(eeig(k2,ieig+i)-eeig(k2,ieig)) < eps) then
                   n = n+1
               else
                   exit
               end if
            end do
#ifdef __TIMER_DO__
            call timer_end(817)
#endif
        !!$write(6,100) (eeig(k2,ieig+i),i=0,20)
        !!$write(6,*) ' 21 states are degenerate! (error)'
        !!$stop 'error === in sub.nstt3i ==='

!c$$$        if(n >= 2) then
!c$$$           write(6,'(" !nstt3i (k2,ieig,n) = ",3i5)') k2, ieig, n
!c$$$           write(6,'(" !nstt3i ",10x," eeig = ",6f8.4)')
!c$$$     &          (eeig(k2,ieig+i),i=0,n-1)
!c$$$        end if
#ifdef __TIMER_DO__
            call timer_sta(818)
#endif
            rn_1 = 1.0d0/real(n,kind=DP)
            do ie=0, ne
               c1 = 0.0d0
               c2 = 0.0d0
               do  i=0, n-1
                  c1 = c1+cdos(k2,ieig+i,ie)
                  c2 = c2+cind(k2,ieig+i,ie)
               end do
               c1 = c1*rn_1
               c2 = c2*rn_1
               do i=0, n-1
                  cdos(k2,ieig+i,ie) = c1
                  cind(k2,ieig+i,ie) = c2
               end  do
            end do
#ifdef __TIMER_DO__
            call timer_end(818)
#endif
            ieig = ieig+n
            if(ieig < neig) goto 40
         end do
!$OMP end do
!$OMP end parallel
!$$$ 31   call gettod(wct_now)
!$$$      etime = etime + (wct_now-wct_start)*1.d-6
!$$$      write(6,'(" !nstt3i etime = ",f16.8)') etime
#ifdef __TIMER_SUB__
         call timer_end(713)
#endif
!$$$  100 format(' e=',4e16.8/(3x,4e16.8))
         call tstatc0_end(id_sname)
      end subroutine nstt3i
!--*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nsttod(eb,ieb)
         use m_Const_Parameters, only :  DP
         implicit none
         integer, intent(inout), dimension(4) :: ieb
         real(DP), intent(inout), dimension(4) :: eb
         integer :: i, k, ih, ind
         real(DP) :: a

         do  k=1, 3
            a = eb(k)
            ih = ieb(k)
            ind = k
            do  i=k+1, 4
               if(eb(i) < a) then
                  a = eb(i)
                  ih = ieb(i)
                  ind = i
               end if
            end do
            eb(ind) = eb(k)
            ieb(ind) = ieb(k)
            eb(k) = a
            ieb(k) = ih
         end do
      end  subroutine nsttod
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstts1(e1,e2,e3,e4)
         use m_Const_Parameters, only :  DP
         use m_Control_Parameters, only : width_tetra
         use m_Spg_plus_Tetra_Common_Blocks, only : tetra_eps
          implicit none
          real(DP), intent(inout) :: e1, e2, e3, e4
          real(DP) :: eps, eps2, eps3
          real(DP) :: a21, a32, a43, o

!         eps = 1.0d-4
          eps = tetra_eps

          a21 = dabs(e2-e1)
          a32 = dabs(e3-e2)
          a43 = dabs(e4-e3)


          eps2 = eps*0.5d0
          eps3 = eps*1.5d0

          if(a21 < eps) then
             if(a32 < eps) then
                if(a43 < eps) then
                   o = (e1+e2+e3+e4)*0.25d0
                   e1 = o-eps3
                   e2 = o-eps2
                   e3 = o+eps2
                   e4 = o+eps3
                else
                   e1 = e3-eps*2.0d0
                   e2 = e3-eps
                end if
             else
                if(a43 < eps) then
                   e1 = e2-eps
                   e4 = e3+eps
                else
                   e1 = e2-eps
                end if
             end if
          else
             if(a32 < eps) then
                if(a43 < eps) then
                   e3 = e2+eps
                   e4 = e2+eps*2.0d0
                else
                   o = (e2+e3)*0.5d0
                   e2 = o-eps2
                   e3 = o+eps2
                end if
             else
                if(a43 < eps) then
                   e4 = e3+eps
                end if
             end if
          end if
      end  subroutine nstts1

      subroutine rdprp(jpr,cname,idim,il,ngen,inv,igen,jgen, &
     &                  imag,ianti,janti, &
     &                  a,b,c,ca,cb,cc,nfspg,jpri_spg)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
         integer, intent(in) :: nfspg, jpri_spg
         integer, intent(out) :: jpr, idim, il, ngen, inv, imag
         integer, intent(out), dimension(3) :: igen
         integer, intent(out), dimension(2,3,3) :: jgen
         integer, intent(out) :: ianti
         integer, intent(out), dimension(2,3) :: janti
         integer :: i, j, k
         real(DP) :: a, b, c, ca, cb, cc
         character(len=60), intent(out)  :: cname

         rewind nfspg

         read(nfspg,*) jpr
         if(jpr >= jpri_spg) jpr = jpri_spg
         read(nfspg,800) cname
         read(nfspg,*) idim, il, ngen, inv

!      write(6,820) cname,idim,il,ngen,inv
!      write(6,*) 'ngen=',ngen

         if(jpri_spg >= 1) then
            if(printable) then
               write(nfout,820) cname,idim,il,ngen,inv
               write(nfout,*) 'ngen=',ngen
            end if
         end if

         do i=1, 3
            read (nfspg,*)  igen(i),((jgen(j,k,i),j=1,2),k=1,3)
         end do

  800    format(a60)
  820    format(' == ',a60,' ==' / &
        &  ' dimension=',i2,'      il=',i2,'   ngen=',i2,'   inv=',i2)
  840    format(i5,3(3x,i3,' /',i3))

         read(nfspg,*) imag

         if(imag == 1) then
            if(printable) then
               write(nfout,*) ' antiferromagnetic calculation'
            end if
         end if
         if(imag == 1) then
           read (nfspg,*)  ianti,((janti(j,k),j=1,2),k=1,3)
         endif

! --* lattice paremeters
!      read(5,*) a,b,c,ca,cb,cc
        read(nfspg,*) a,b,c,ca,cb,cc

!      write(6,*) '   '
!      write(6,860) a,b,c,ca,cb,cc
!      860 format(' a, b, c =',3f12.6/
!      &       'ca,cb,cc =',3f12.6)

         if(jpri_spg >= 1) then
            if(printable) then
               write(nfout,*) '   '
               write(nfout,860) a,b,c,ca,cb,cc
            end if
         end if
 860     format(' a, b, c =',3f12.6/ &
        &       'ca,cb,cc =',3f12.6)
      end  subroutine rdprp
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine setkp0(np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     & ,nx,ny,nz,nxx,nyy,nzz &
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     & ,nstar2,pa0,pb0,pb,ka0,ka2 &
     & ,nfspg,ipri_kp,ipri_spg, itrs )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : eeule2,eeulv2,ieule2,ieulv2,ig02,iv02 ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : tca,tac,tab,tba,tcb,tbc ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : grc,gkc,gra,gka,grb,gkb ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1
         implicit none
         integer, intent(in) :: lmnp0, lmnp1, lmnp2, ipri_kp, ipri_spg,&
                              &  nfspg, itrs
         integer :: np2, np1, np0
         integer :: nx, ny, nz, nxx, nyy, nzz
         integer :: jpr, idim, ill, ngen, ianti
         integer :: i, j, inv, nx1, ny1, nz1, nd, imag
         integer, dimension(3) :: igen
         integer, dimension(2,3,3) :: jgen
         integer, dimension(2,3) :: janti
         real(DP) :: eps
         real(DP) :: a, b, c, ca, cb, cc
         real(DP), dimension(3) :: omove
         character(len=60) :: cname
!      character*5 schoen(48)
!     integer kanti(2,3)
!
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!      real*8  tabn(3,3),ta1n(3,48)
!      integer  lra1n(3,3,48)
!
!      real*8  pa(3)
!      real*8  pb (3,lmnp2)
!      integer ka(4)
         integer :: jpri_spg
         integer :: ncounter
         integer, dimension(4,lmnp0) :: ka0
         integer, dimension(4,lmnp2) :: ka2
!     integer kn(4,0:lmnl)
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) :: ip01, ip21, iu21, iv21
         integer, dimension (lmnp2) :: ip02, ip12, nstar2
         real(DP), dimension(3,lmnp0) :: pa0, pb0
         real(DP), dimension(3,lmnp2) :: pb

        data    ncounter/0/
        save    ncounter

         ncounter = ncounter + 1

         if(jpri_spg  >= 1 .and. ncounter  <= 1)then
            jpri_spg = 1
          else
            jpri_spg = 0
         end if

         eps = 1.0d-5

         call rdprp(jpr,cname,idim,ill,ngen,inv,igen,jgen, &
     &            imag,ianti,janti, &
     &            a,b,c,ca,cb,cc,nfspg,jpri_spg)

         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
     &           omove,jpri_spg)

         if(idim == 1) then
            if(nx >= 0) nx=0
            ny=0
         end if

         if(idim == 2) nz=0

         if(ipri_kp >= 1 .and. ncounter <= 1) then
            jpr = 3
         else
           jpr = -1
         end if

         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)

         call nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
     &             pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
     &             iu21,iv21, itrs )
         if(idim == 1) then
            do i=1, np0
               if(dabs(pa0(1,i)) > eps .or. dabs(pb0(1,i)) > eps .or. &
               &  dabs(pa0(2,i)) > eps .or. dabs(pb0(2,i)) > eps) then
                  if(printable) write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if

         if(idim == 2) then
            do i=1, np0
               if(dabs(pa0(3,i)) > eps .or. dabs(pb0(3,i)) > eps) then
                  if(printable) write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if
  900    format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
         pb(1:3,1:np2) = pb0(1:3,ip02(1:np2))
         nstar2(1:np2) = 0
         do i = 1, np1
            nstar2(ip21(i)) = nstar2(ip21(i))+1
         end do
!========================

!  100 format(8i10)
  120    format(i5,5x,3f9.3)
  200    format((3d24.16))
  220    format((3d24.16,i6))
! 300 format(5i12)
      end  subroutine setkp0
! ---------------------
!     subroutine setkp0_n is rewritten from <setkp0> by T. Yamasaki,
      subroutine setkp0_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
     & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     & ,nx,ny,nz,nxx,nyy,nzz &
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     & ,nstar2,pa0,pb0,pb,ka0,ka2 &
     & ,ipri_kp, itrs )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         implicit none
!      integer il1
         integer, intent(in) :: ill, ngen, inv
         integer, intent(in) :: lmnp0, lmnp1, lmnp2, itrs
         integer, intent(in) :: ipri_kp
         integer, intent(in), dimension(ngen) :: igen
         integer :: np0, np1, np2
         integer, intent(in), dimension(2,3,ngen) :: jgen
         integer :: nx, ny, nz, nxx, nyy, nzz
         integer ::  nx1, ny1, nz1, nd
         integer :: idim
         integer :: ipri_spg
         integer :: i, j
         real(DP), intent(in) :: a, b, c, ca, cb, cc
         real(DP) :: eps
         real(DP), dimension(3) :: omove

!        real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!        integer lrb2(3,3),lrbv(3,3)
!        real*8  rrb2(3,3),rrbv(3,3)
!        real*8  rtc2(3,3),rtcv(3,3)
!        real*8  tabn(3,3),ta1n(3,48)
!        integer  lra1n(3,3,48)

!      real*8  pa(3)
         real(DP), dimension(3,lmnp2) :: pb
!      integer ka(4)
         integer, dimension(4,lmnp0) :: ka0
         integer, dimension(4,lmnp2) :: ka2
!      real*8  pn(3,0:lmnl)
!      integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) :: pa0,  pb0
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) :: ip01, ip21
         integer, dimension(lmnp2) :: ip02, ip12
         integer, dimension(lmnp1) :: iu21, iv21
         integer, dimension(lmnp2) ::  nstar2

         eps = 1.0d-5

         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
       & omove,ipri_spg)


         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
          call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)

         call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
       & pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
       & iu21,iv21, itrs )
         if(idim == 1) then
            do i=1, np0
               if(dabs(pa0(1,i)) > eps .or. abs(pb0(1,i)) > eps .or. &
               &  dabs(pa0(2,i)) > eps .or. abs(pb0(2,i)) > eps) then
                  if(printable) then
                     write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  end if
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if

         if(idim == 2) then
            do i=1, np0
               if(dabs(pa0(3,i)) > eps .or. abs(pb0(3,i)) > eps) then
                  if(printable) then
                    write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  end if
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if
  900    format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
         pb(1:3,1:np2) = pb0(1:3,ip02(1:np2))
         nstar2(1:np2) = 0
         do i = 1, np1
            nstar2(ip21(i)) = nstar2(ip21(i))+1
         end do

  120    format(i5,5x,3f9.3)
  200    format((3d24.16))
  220    format((3d24.16,i6))
      end subroutine setkp0_n

      subroutine setkp0_n_kt(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     &     ,nx,ny,nz,nxx,nyy,nzz &
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     &     ,nstar2,pa0,pb0,pb,ka0,ka2 &
     &     ,ipri_kp, &
     &     use_altv_rltv, altv, rltv, itrs, &
     &     gen_name_in_carts )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         implicit none
!     integer il1
         integer, intent(in) :: ill, ngen, inv
         integer, intent(in) :: lmnp0, lmnp1, lmnp2
         integer :: nx, ny, nz, nxx, nyy, nzz
         integer :: nx1, ny1, nz1, nd
         integer :: idim
         integer :: np0, np1, np2
         integer :: ipri_kp
         integer :: jpri_spg
         integer :: i, j
         integer, dimension(ngen) :: igen
         integer, dimension(2,3,ngen) :: jgen
         real(DP) :: a, b, c, ca, cb, cc
         real(DP) :: eps
         real(DP), dimension(3) :: omove
         integer :: irotr2, irotk2

!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!      real*8  tabn(3,3),ta1n(3,48)
!      integer  lra1n(3,3,48)

!      real*8  pa(3)
         real(DP), dimension(3,lmnp2)  :: pb
!      integer ka(4)
         integer, dimension(4,lmnp0)  :: ka0
         integer, dimension(4,lmnp2)  :: ka2
!     real*8  pn(3,0:lmnl)
!     integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) :: pa0, pb0
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) ::  ip01, ip21
         integer, dimension(lmnp2) ::  ip02, ip12
         integer, dimension(lmnp1) ::  iu21, iv21
         integer, dimension(lmnp2) ::  nstar2
! ---
         integer :: itrs
         logical :: use_altv_rltv, use_trs, gen_name_in_carts
         real(DP), dimension(3,3) :: altv, rltv
! -------
!c init
         if(itrs == 1) then
            use_trs = .true.
         else
            use_trs = .false.
         end if

         eps = 1.0d-5

         call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
       & omove,jpri_spg, use_altv_rltv, altv, rltv, &
       & use_trs, gen_name_in_carts )

         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)

         call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
       & pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
       & iu21,iv21, itrs )

         if(idim == 1) then
            do i=1, np0
               if(dabs(pa0(1,i)) > eps .or. dabs(pb0(1,i)) > eps .or. &
                & dabs(pa0(2,i)) > eps .or. dabs(pb0(2,i)) > eps) then
                  if(printable) then
                     write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  end if
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if

         if(idim == 2) then
            do  i=1, np0
               if(dabs(pa0(3,i)) > eps .or. dabs(pb0(3,i)) > eps) then
                   if(printable) then
                      write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                   end if
                  stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if
 900     format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)
         pb(1:3,1:np2) = pb0(1:3,ip02(1:np2))
         nstar2(1:np2)=0
         do i = 1, np1
            nstar2(ip21(i)) = nstar2(ip21(i))+1
         end do

 120     format(i5,5x,3f9.3)
 200     format((3d24.16))
 220     format((3d24.16,i6))
      end  subroutine setkp0_n_kt

! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp, np,ka0)
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : ka00_called
!$ use omp_lib
         implicit none
         integer, intent(in) :: nxx, nyy, nzz, nx1, ny1, nz1
         integer, intent(in) :: nd, lmnp
         integer, intent(inout) :: np
         integer, intent(out), dimension(4,lmnp) :: ka0
         integer :: i, ix, iy, iz, ip
         integer :: ny1nz1,  nx1nz1, nx1ny1
         integer,save :: nps
         integer, pointer, dimension(:,:), save :: ka0s


         if(ka00_called) then
            np = nps
!$OMP simd
            do i = 1, np
               ka0(1:4,i) = ka0s(1:4,i)
            end do
!$OMP end simd
            return
         end if

!        ka0(4,:) = nd
         ka0=nd

         np = 0
         do iz=0, nzz
            do iy=0, nyy
                do ix=0, nxx
                 np = np+1
!                  if(np > lmnp) then
!                    if(printable) then
!                       write(nfout,*) 'ix,iy,iz=',ix,iy,iz,'   np,lmnp=',np,lmnp
!                    end if
!                    stop ' === stop in sub.ka00. (np>lmnp) ==='
!                  end if
                  ka0(1,np) = ix
                  ka0(2,np) = iy
                  ka0(3,np) = iz
!                 ka0(4,np) = nd
!                  ka0(4,np) = nd
!$$$        write(6,'(" !spg  (",3i5,") : ka0 (1:4",i8,") = ",4i8)')
!$$$     &       ix,iy,iz,np,ka0(1:4,np)
              end do
           end do
        end do

         ny1nz1 = ny1*nz1
         nx1nz1 = nx1*nz1
         nx1ny1 = nx1*ny1

!$OMP simd
        do i = 1, np
           ka0(1,i) = ka0(1,i)*ny1nz1
           ka0(2,i) = ka0(2,i)*nx1nz1
           ka0(3,i) = ka0(3,i)*nx1ny1
        end do
!$OMP end simd

        allocate(ka0s(4,lmnp))

        nps = np
!$OMP simd
        do i = 1, np
           ka0s(1:4,i) = ka0(1:4,i)
        end do
!$OMP end simd
        ka00_called = .true.
      end subroutine ka00

      subroutine setkp0_default(nbztyp1,altv,nx,ny,nz &
     & ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     & ,nxx,nyy,nzz &
     & ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     & ,nstar2,pa0,pb0,pb,ka0,ka2 &
     & ,ipri_kp,ipri_spg, itrs )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         implicit none
!      character*5 schoen(48)
!      character*60 cname
         integer, intent(in) :: nbztyp1
         integer, intent(in) :: lmnp0, lmnp1, lmnp2
         integer, intent(in) :: ipri_kp, ipri_spg, itrs
         integer :: ill, idim, inv, jpr
         integer :: nx, ny, nz, nxx, nyy, nzz, nx1, ny1, nz1, nd
         integer :: ngen
         integer :: np0, np1, np2
         integer :: i, j
         integer, dimension(3) :: igen
         integer, dimension(2,3,3) :: jgen
         real(DP) one, eps
         real(DP) :: a, b, c, ca, cb, cc
         real(DP), dimension(3,3) :: altv
         real(DP), dimension(3) :: omove
!      integer janti(2,3),kanti(2,3)
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!      real*8  tabn(3,3),ta1n(3,48)
!      integer  lra1n(3,3,48)

!      real*8  pa(3)
         real(DP), dimension(3,lmnp2) ::  pb
!       integer ka(4)
         integer, dimension(4,lmnp2) :: ka2
!     real*8  pn(3,0:lmnl)
!     integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) :: pa0, pb0
         integer, dimension(lmnp0) ::  ip10, ip20
         integer, dimension(lmnp1) :: ip01, ip21
         integer, dimension(lmnp2) :: ip02, ip12
         integer, dimension(lmnp1) :: iu21, iv21
         integer, dimension(lmnp2) ::  nstar2
         integer, dimension(4,lmnp0) ::  ka0
         integer :: ncounter
         data    ncounter/0/
         save    ncounter

         ncounter = ncounter + 1

         one = 1.0d0
         eps = 1.0d-5

         if(printable) write(nfout,'(" << setkp0_default >>")')
         if(ipri_spg  >= 1 .and. ncounter <= 1 ) &
        &  write(nfout,*) 'nbztyp(nbztyp_spg) ', nbztyp1

         a = dsqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
         b = dsqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
         c = dsqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
         ca = (altv(1,2)*altv(1,3)+altv(2,2)*altv(2,3)+altv(3,2)*altv(3,3)) &
          & /(b*c)
         cb = (altv(1,3)*altv(1,1)+altv(2,3)*altv(2,1)+altv(3,3)*altv(3,1)) &
         & /(c*a)
         cc =(altv(1,1)*altv(1,2)+altv(2,1)*altv(2,2)+altv(3,1)*altv(3,2)) &
         & /(a*b)

         if(ipri_spg >= 1) then
            if(printable) Then
               write(nfout,*) '   '
               write(nfout,860) a,b,c,ca,cb,cc
            end if
         end if
  860    format(' a, b, c =',3f12.6/ &
        &  'ca,cb,cc =',3f12.6)

         if(nbztyp1 == 2) then
            idim = 3
            ill = 1
            ngen = 3
            inv = 1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
         else if(nbztyp1 == 3) then
            idim =3
            ill =3
            ngen = 3
            inv = 1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
           jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
      else if(nbztyp1 == 4) then
            idim =3
            ill =2
            ngen =3
            inv =1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
         else if(nbztyp1 == 5) then
            idim = 3
            ill = 2
            ngen = 3
            inv = 0
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
!cccccccc 2nd choice ccccccccccccccccccccc
            igen(2) = 19
            jgen(1,1,2) = 1
            jgen(2,1,2) = 4
            jgen(1,2,2) = 1
            jgen(2,2,2) = 2
            jgen(1,3,2) = 3
            jgen(2,3,2) = 4
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
!ccccccc 2nd choice ccccccccccccccccccccc
         else if(nbztyp1 == 6) then
            idim = 3
            ill = 0
            ngen = 2
            inv = 0
            igen(1) = 3
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 2
            jgen(2,3,1) = 3
            igen(2) = 10
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
         end if

         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
       & omove,ipri_spg)

         if(idim == 1) then
            if(nx >= 0) nx = 0
            ny = 0
         end if
         if(idim == 2) nz = 0
         if(ipri_kp  >= 1 .and. ncounter <= 1) then
            jpr = 3
         else
            jpr = -1
         end if

         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
          call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)

         call nskpb0(jpr,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
       &             pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
       &            iu21,iv21, itrs )

         if(idim == 1) then
            do i=1,  np0
               if(dabs(pa0(1,i)) > eps .or. dabs(pb0(1,i)) > eps .or. &
               &  dabs(pa0(2,i)) > eps .or. abs(pb0(2,i)) > eps) then
                  if(printable) then
                     write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  end if
                   stop '=bnkpgn(pa0,pb0)='
               end if
            end do
         end if

         if(idim == 2) then
            do i=1, np0
               if(dabs(pa0(3,i)) > eps .or. abs(pb0(3,i)) > eps) then
                  if(printable) then
                     write(nfout,900) i,(pa0(j,i),j=1,3),(pb0(j,i),j=1,3)
                  end if
                  stop '=bnkpgn(pa0,pb0)='
              end if
            end do
         end if
  900   format(i5,'   pa0=',3f12.6/5x,'   pb0=',3f12.6)

         pb(1:3,1:np2) = pb0(1:3,ip02(1:np2))
         nstar2(1:np2) = 0
        do i = 1, np1
           nstar2(ip21(i)) = nstar2(ip21(i))+1
        end do

!========================

! 100 format(8i10)
  120    format(i5,5x,3f9.3)
  200    format((3d24.16))
  220    format((3d24.16,i6))
!  300 format(5i12)
      end subroutine setkp0_default

! --- subroutine setkp0_default_n ---
!     This subroutine is revised from <setkp0_default> by T. Yamasaki,
!                                                      31th May 2003
      subroutine setkp0_default_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
        &     ,nx,ny,nz &
        &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
        &     ,nxx,nyy,nzz &
        &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
        &     ,nstar2,pa0,pb0,pb,ka0 &
        &     ,ipri_kp, itrs )

         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         use m_Spg_plus_Tetra_Common_Blocks, only : &
         setkp0_default_n_called
         use m_Timing, only : tstatc0_begin, tstatc0_end
!$ use omp_lib
         implicit  none
         integer :: ill,ngen,inv
         integer :: idim, ipri_spg
         integer, dimension(ngen) :: igen
         integer, dimension(2,3,ngen)  :: jgen
         real(DP) ::  a,b,c,ca,cb,cc
         integer  :: nx,ny,nz,nxx, nyy, nzz, np2, np1, np0, lmnp0, lmnp1,lmnp2
         integer  :: nx1, ny1, nz1, nd
         integer  :: i, j
         real(DP) ::  omove(3)
         integer  :: ipri_kp
         integer  :: itrs

!      real*8  pa(3)
         real(DP), dimension(3,lmnp2) ::  pb
!      real*8  pn(3,0:lmnl)
!      integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) ::   pa0, pb0
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) ::  ip01, ip21
         integer, dimension(lmnp2) :: ip02, ip12
         integer, dimension(lmnp1) ::  iu21, iv21
         integer, dimension(lmnp2) ::  nstar2
         integer, dimension(4,lmnp0) ::  ka0

         real(DP), allocatable, dimension(:,:),save :: pbs
         integer, allocatable, dimension(:),save :: nstar2s
         integer :: id_sname = -1
         call tstatc0_begin('setkp0_default_n ',id_sname)
         if(ipri_kp >= 1) then
            if(printable) then
               write(6,'(" << setkp0_default_n >>")')
               write(6,'(" !! ill = ",2i6)') ill
               write(6,'(" !! ngen = ",i6)') ngen
            end if
         end if

         if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
         if(ipri_kp >= 1) then
            if(printable) then
               do j = 1, ngen
                  write(6,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j) &
                 & ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j ) &
                 & ,jgen(2,3,j)
               enddo
            end if
         end if

         idim = 3
         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        &           omove,ipri_spg)
         if(ipri_kp >= 1) then
            if(printable) write(6,'(" !! ill, il = ",2i6)') ill,il
         end if

         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
        call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)

!$o,p para;;e;
!      call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2,
!     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &     iu21,iv21, itrs )
         if(itrs == 0) then
            call nskpb0_s_1(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
           &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
           &     iu21,iv21, itrs)
         else
            call nskpb0_s_2(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
           &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
           &     iu21,iv21, itrs)
         end if

!----- OMP directives and modifications T. Hamada 2022/1/2 -----
         if(setkp0_default_n_called) then
!$OMP simd
            do i = 1, np2
               pb(1:3,i) = pbs(1:3,i)
            end do
!$OMP end simd
             nstar2 = nstar2s
            return
         end if
         call set_pb
         call set_nstar2
         allocate(pbs(3,lmnp2))
         allocate(nstar2s(lmnp2))
!$OMP simd
         do i = 1, np2
            pbs(1:3,i) = pb(1:3,i)
         end do
!$OMP end simd
         nstar2s = nstar2
         setkp0_default_n_called = .true.
!-----------------------------------------------------------------
         if(ipri_kp >= 1) then
            if(printable) then
               write(nfout,'(" !! ngen = ",i6," <<- s,etkp0_default_n")') ngen

            end if
         end if

         call tstatc0_end(id_sname)

         contains
          subroutine set_pb
!$OMP simd
            do i = 1, np2
               pb(1:3,i)  = pb0(1:3,ip02(i))
            end do
!$OMP end simd
          end subroutine set_pb

          subroutine set_nstar2
             nstar2 = 0
!$OMP simd
             do i = 1, np1
               nstar2(ip21(i)) = nstar2(ip21(i))+1
            end do
!$OMP end simd
         end subroutine set_nstar2
        end subroutine setkp0_default_n

! ======================================= added by K. Tagami ================ 12.0A
      subroutine setkp0_default_n_kt2(ill,ngen,inv,igen,jgen &
     &     , a,b,c,ca,cb,cc &
     &     ,nx,ny,nz &
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     &     ,nxx,nyy,nzz &
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     &     ,nstar2,pa0,pb0, pa, ka0 &
     &     ,ipri_kp &
     &     ,use_altv_rltv, altv, rltv, itrs &
     &     ,gen_name_in_carts )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         implicit none
         integer :: ill, ngen, inv
         integer, dimension(ngen) :: igen
         integer, dimension(2,3,ngen) :: jgen
         real(DP)  :: a,b,c,ca,cb,cc
         integer   :: idim, ipri_spg
         integer  :: nx, ny, nz, nxx, nyy, nzz, np2, np1, np0, lmnp0, lmnp1, lmnp2
         integer  :: nx1, ny1, nz1, nd
         integer  :: i, j
         integer  :: ipri_kp
         real(DP), dimension(3) ::  omove
!      real*8  pa(3)
         real(DP), dimension(3,lmnp2) :: pa
!     real*8  pn(3,0:lmnl)
!     integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) :: pa0, pb0
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) :: ip01, ip21
         integer, dimension(lmnp2) :: ip02, ip12
         integer, dimension(lmnp1) :: iu21, iv21
         integer, dimension(lmnp2) :: nstar2
         integer, dimension(4,lmnp0) ::  ka0
! ----
         integer :: itrs
         logical use_altv_rltv, use_trs, gen_name_in_carts
         real(DP), dimension(3,3) :: altv, rltv
! init
         if(itrs == 1) then
            use_trs = .true.
         else
            use_trs = .false.
         end if

         if(ipri_kp >= 1) then
            if(printable) then
               write(nfout,'(" << setkp0_default_n >>")')
               write(nfout,'(" !! ill = ",2i6)') ill
               write(nfout,'(" !! ngen = ",i6)') ngen
            end if
         end if

         if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
         if(ipri_kp >= 1) then
            if(printable) then
               do j = 1, ngen
                  write(nfout,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j) &
                & ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j) &
                & ,jgen(2,3,j)
               enddo
            end if
         end if
         idim = 3
         call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
            &          omove,ipri_spg, use_altv_rltv, altv, rltv, &
           &           use_trs, gen_name_in_carts )

         if(ipri_kp >= 1) then
            if(printable) write(nfout,'(" !! ill, il = ",2i6)') ill,il
         end if
         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)
! ----------------------------------------------
!      if ( kmode .eq. 1 ) then
!         call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,
!     &               lmnp1,lmnp2,
!     &        pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &        iu21,iv21 )
!      else if ( kmode .eq. 2 ) then
!         call nskpa0_kt(ipri_kp,tab,ng1,lsa1,iv1,np0,pa0,lmnp0,
!     &                  lmnp1,lmnp2,
!     &        pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12,
!     &        iu21,iv21 )
!      endif
! ----------------- T. Hamada 2021.12.23 --------------------------
         call nskpa0_kt(ipri_kp,tab,ng1,lsa1,iv1,np0,pa0,lmnp0, &
     &     lmnp1,lmnp2, &
     &     pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
     &     iu21,iv21, itrs )
!         if(itrs == 0) then
!            call nskpb0_s_1(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,&
!           &  lmnp2, &
!           &  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
!           &  iu21,iv21, itrs )
!         else
!            call nskpb0_s_2(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,&
!           &  lmnp2, &
!           &  pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
!           &  iu21,iv21, itrs )
!         end if
! -----------------------------------------------------------------
         pa(1:3,1:np2) = pa0(1:3,ip02(1:np2))
         nstar2(1:np2) = 0
! ----------------------------------------------------------------
         do i = 1, np1
            nstar2(ip21(i)) = nstar2(ip21(i))+1
         end do

         if(ipri_kp >= 1 ) then
            if(printable) then
               write(nfout,'(" !! ngen = ",i6," <<- setkp0_default_n")') ngen
            end if
         end  if
      end subroutine setkp0_default_n_kt2

      subroutine setkp0_default_n_kt(ill,ngen,inv,igen,jgen &
     &     ,a,b,c,ca,cb,cc &
     &     ,nx,ny,nz &
     &     ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
     &     ,nxx,nyy,nzz &
     &     ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
     &     ,nstar2,pa0,pb0, pb, ka0 &
     &     ,ipri_kp &
     &     ,use_altv_rltv, altv, rltv, itrs &
     &     ,gen_name_in_carts )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & setkp0_default_n_kt_called
         use m_Timing, only : tstatc0_begin, tstatc0_end
!$ use omp_lib
         implicit none
         integer, intent(in) ::  ill, ngen, inv
         integer, dimension(ngen) ::  igen
         integer :: ipri_kp
         integer :: idim, jpri_spg
         integer :: i, j
         integer, dimension(2,3,ngen) :: jgen
         real(DP) :: a, b, c, ca, cb, cc
         integer :: nx, ny, nz, nxx, nyy, nzz, np2, np1, np0, lmnp0, lmnp1, lmnp2
         integer :: nx1, ny1, nz1, nd
         real(DP), dimension(3) :: omove
         integer :: irotr2, irotk2
!      real*8  pa(3)
         real(DP), dimension(3,lmnp2)  :: pb
!      real*8  pn(3,0:lmnl)
!      integer kn(4,0:lmnl)
         real(DP), dimension(3,lmnp0) :: pa0, pb0
         integer, dimension(lmnp0) :: ip10, ip20
         integer, dimension(lmnp1) :: ip01, ip21
         integer, dimension(lmnp2) :: ip02, ip12
         integer, dimension(lmnp1) :: iu21, iv21
         integer, dimension(lmnp2) :: nstar2
         integer, dimension(4,lmnp0) :: ka0

         integer itrs
         logical :: use_altv_rltv, use_trs, gen_name_in_carts
         real(DP), dimension(3,3) :: altv, rltv

         real(DP), allocatable, dimension(:,:), save :: pbs
         integer, allocatable, dimension(:), save :: nstar2s
         integer :: id_sname = -1
! -------
! init

         call tstatc0_begin('setkp0_default_n_kt ',id_sname)

         if(itrs == 1) then
            use_trs = .true.
         else
            use_trs = .false.
         end if

         if(ipri_kp >= 1) then
            if(printable) then
               write(nfout,'(" << setkp0_default_n >>")')
               write(nfout,'(" !! ill = ",2i6)') ill
               write(nfout,'(" !! ngen = ",i6)') ngen
            end if
         end if

         if(ngen > 3) stop ' ngen >3 <<setkp0_default_n>>'
         if(ipri_kp >= 1) then
            if(printable) then
               do j = 1, ngen
                  write(nfout,'(" !!  igen, jgen = ",7i6)') igen(j), jgen(1,1,j) &
                  & ,jgen(2,1,j),jgen(1,2,j),jgen(2,2,j),jgen(1,3,j) &
                  & ,jgen(2,3,j)
               enddo
            end if
         end if
         idim = 3
         call tbspg_kt(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        &             omove,jpri_spg, use_altv_rltv, altv, rltv, &
        &             use_trs, gen_name_in_carts )

         if(ipri_kp >= 1) then
            if(printable) write(nfout,'(" !! ill, il = ",2i6)') ill,il
         end if
         call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)
         call ka00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0)
         call nskp00(nxx,nyy,nzz,nx1,ny1,nz1,nd,lmnp0, np0,ka0,pa0)
! ------------------------ T. Hamada 2021.12.23 -----------------------
         if(itrs == 0) then
            call nskpb0_s_1(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
           & pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
           & iu21,iv21, itrs)
         else
            call nskpb0_s_2(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
           & pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
           & iu21,iv21, itrs)
         end if

!        call nskpb0(ipri_kp,tab,ng1,lsb1,iv1,np0,pa0,lmnp0,lmnp1,lmnp2, &
!       & pb0,np1,np2,ip10,ip20,ip01,ip21,ip02,ip12, &
!       & iu21,iv21, itrs )
!       pb(1:3,1:np2)=pb0(1:3,ip02(1:np2))
! ---------------------------------------------------------------------
! ----- OMP directives and modifications by T. Hamada 2022/1/2 ------
         if(setkp0_default_n_kt_called) then
            pb = pbs
            nstar2 = nstar2s
            return
         end if
!$OMP parallel
         call set_pb
         call set_nstar2
!$OMP end parallel
         allocate(pbs(3,lmnp2))
         allocate(nstar2s(lmnp2))
         pbs = pb
         nstar2s = nstar2
         setkp0_default_n_kt_called = .true.
! ---------------------------------------------------------------
         if(ipri_kp >= 1) then
            if(printable) then
               write(nfout,'(" !! ngen = ",i6," <<- setkp0_default_n")') ngen
            end if
         end if

         call tstatc0_end(id_sname)

         contains
          subroutine set_pb
!$OMP simd
             do i = 1, np2
               pb(1,i)   = pb0(1,ip02(i))
               pb(2,i)   = pb0(2,ip02(i))
               pb(3,i)   = pb0(3,ip02(i))
            end do
!$OMP end simd
       end subroutine set_pb

          subroutine set_nstar2
            nstar2(1:np2) = 0
!$OMP simd
            do i = 1, np1
                nstar2(ip21(i)) = nstar2(ip21(i))+1
            end do
!$OMP end simd
      end subroutine set_nstar2

      end subroutine setkp0_default_n_kt
! =========================================================================== 12.0A
      subroutine setspg(tabn,ng1n,ta1n,lra1n,nfspg,ipri_spg)
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 !  common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 !  common/nspg07
         implicit none
         character(len=60) ::  cname
         integer, intent(in) :: nfspg, ipri_spg
         integer :: jpr, ill, idim, ngen, inv, imag, ianti
         integer, dimension(3) :: igen
         integer, dimension(2,3,3) :: jgen
         integer, dimension(2,3) ::  janti
         integer :: ng1n
!     integer kanti(2,3)
         real(DP) :: one, eps
         real(DP) :: a, b, c, ca, cb, cc
         real(DP), dimension(3,3) :: tabn
         real(DP), dimension(3,48) :: ta1n
         real(DP), dimension(3) :: omove
         integer  :: irotr2, irotk2
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!      real*8  tabn(3,3),ta1n(3,48)
         integer, dimension(3,3,48) ::  lra1n
         integer  ::  ncounter
         integer :: jpri_spg
         data  ncounter/0/
         save ncounter

         ncounter = ncounter + 1
         if(ipri_spg  >= 1 .and. ncounter  <= 1) then
            jpri_spg = 1
         else
            jpri_spg = 0
         end if

         one = 1.0d0
         eps = 1.0d-5

        call rdprp(jpr,cname,idim,ill,ngen,inv,igen,jgen, &
        &            imag,ianti,janti, &
        &            a,b,c,ca,cb,cc,nfspg,jpri_spg)

         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        & omove,jpri_spg)


         tabn(1:3,1:3) = tab(1:3,1:3)
         ng1n = ng1
         ta1n(1:3,1:ng1n) = ta1(1:3,1:ng1n)
         lra1n(1:3,1:3,1:ng1n) = lra1(1:3,1:3,1:ng1n)
      end  subroutine setspg

      subroutine setspg_n(ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
     &     ,tabn,ng1n,ta1n,lra1n,ig01n)
!        Modified by T. Yamasaki(FUJITSU Lab.), 31 May 2003
         use m_Const_Parameters, only :  DP
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07

         implicit none
         integer :: ill, ngen,inv
         integer :: ng1n
         integer, intent(in), dimension(ngen) :: igen
         integer, intent(in), dimension(2,3,ngen) :: jgen
         integer :: idim
         integer :: jpri_spg
         integer, intent(out), dimension(48) :: ig01n
         integer, intent(out), dimension(3,3,ng1) :: lra1n
         real(DP), intent(in) :: a, b, c, ca, cb, cc
         real(DP), intent(out), dimension(3,3) :: tabn
         real(DP), intent(out), dimension(3,ng1) :: ta1n
         real(DP), dimension(8) :: omove

!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!$$$      real*8  tabn(3,3),ta1n(3,48)
!$$$      integer  lra1n(3,3,48),ig01n(48)

         idim = 3
         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        & omove,jpri_spg)

         tabn(1:3,1:3) = tab(1:3,1:3)
         ng1n = ng1
         ta1n(1:3,1:ng1n) = ta1(1:3,1:ng1n)
         lra1n(1:3,1:3,1:ng1n) = lra1(1:3,1:3,1:ng1n)
         ig01n(1:ng1n) = ig01(1:ng1n)
      end subroutine setspg_n

      subroutine setspg_default(nbztyp1,altv,tabn,ng1n,ta1n,lra1n &
     & ,ipri_spg)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         implicit none
         integer, intent(in) :: nbztyp1, ipri_spg
         integer :: ng1n
         integer, intent(out), dimension(3,3,3) :: lra1n
         integer :: irotr2, irotk2
         integer :: idim, ill, ngen, inv
         integer, dimension(3) :: igen
         integer, dimension(2,3,3) :: jgen
         real(DP) :: one, eps, a, b, c, ca, cb, cc
         real(DP), intent(in), dimension(3,3) :: altv
         real(DP), intent(out), dimension(3,3) :: tabn
         real(DP), intent(out), dimension(3,48) :: ta1n
         real(DP), dimension(3) :: omove
!      integer janti(2,3),kanti(2,3)
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
         integer  ::   ncounter
         data     ncounter /0/
         save     ncounter

         ncounter = ncounter + 1

         one = 1.0d0
         eps = 1.0d-5

         if(ipri_spg  >= 1 .and. ncounter <= 1) then
            if(printable) then
               write(nfout,*) 'nbztyp ', nbztyp1
            end if
         end if
         a=dsqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
         b=dsqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
         c=dsqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
         ca=(altv(1,2)*altv(1,3)+altv(2,2)*altv(2,3)+altv(3,2)*altv(3,3))   &
        & /(b*c)
         cb=(altv(1,3)*altv(1,1)+altv(2,3)*altv(2,1)+altv(3,3)*altv(3,1)) &
        & /(c*a)
         cc=(altv(1,1)*altv(1,2)+altv(2,1)*altv(2,2)+altv(3,1)*altv(3,2)) &
        & /(a*b)

         if(ipri_spg  >= 1 .and. ncounter  <= 1) then
            if(printable) then
               write(nfout,*) '   '
               write(nfout,860) a,b,c,ca,cb,cc
            end if
         end if
  860    format(' a, b, c =',3f12.6/ &
        &     'ca,cb,cc =',3f12.6)

         if(nbztyp1 == 2) then
            idim = 3
            ill = 1
            ngen = 3
            inv = 1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) =25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
         else if(nbztyp1 == 3) then
            idim =3
            ill = 3
            ngen = 3
            inv = 1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
         else if(nbztyp1 == 4) then
            idim = 3
            ill = 2
            ngen = 3
            inv = 1
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
            igen(2) = 19
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
         else if(nbztyp1 == 5) then
            idim = 3
            ill = 2
            ngen = 3
            inv = 0
            igen(1) = 5
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 0
            jgen(2,3,1) = 1
!ccccccc 2nd choice ccccccccccccccccccccc
            igen(2) = 19
            jgen(1,1,2) = 1
            jgen(2,1,2) = 4
            jgen(1,2,2) = 1
            jgen(2,2,2) = 2
            jgen(1,3,2) = 3
            jgen(2,3,2) = 4
            igen(3) = 25
            jgen(1,1,3) = 0
            jgen(2,1,3) = 1
            jgen(1,2,3) = 0
            jgen(2,2,3) = 1
            jgen(1,3,3) = 0
            jgen(2,3,3) = 1
!ccccccc 2nd choice ccccccccccccccccccccc
         else if(nbztyp1 == 6) then
            idim = 3
            ill = 0
            ngen = 2
            inv = 0
            igen(1) = 3
            jgen(1,1,1) = 0
            jgen(2,1,1) = 1
            jgen(1,2,1) = 0
            jgen(2,2,1) = 1
            jgen(1,3,1) = 2
            jgen(2,3,1) = 3
            igen(2) = 10
            jgen(1,1,2) = 0
            jgen(2,1,2) = 1
            jgen(1,2,2) = 0
            jgen(2,2,2) = 1
            jgen(1,3,2) = 0
            jgen(2,3,2) = 1
         end if

         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        & omove,ipri_spg)

         ng1n = ng1
         tabn(1:3,1:3) = tab(1:3,1:3)
         ta1n(1:3,1:ng1n) = ta1(1:3,1:ng1n)
         lra1n(1:3,1:3,1:ng1n) = lra1(1:3,1:3,1:ng1n)
      end subroutine setspg_default

      subroutine setspg_default_n(ill,ngen,inv,igen,jgen,imag,iaf,jaf &
     &     ,a,b,c,ca,cb,cc &
     &     ,tabn,ng1n,ta1n,lra1n,ipri_spg,ig01n)
!                     Modified by T. Yamasaki(FUJITSU Lab.) 31 May 2003.
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : irot ! common/nspg0/
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
!$ use omp_lib
         implicit none
         integer ill, ngen,inv, ipri_spg
         integer :: ng1n
         integer :: idim
         integer :: jj
         integer, dimension(ngen) :: igen
         integer, dimension(2,3,ngen) :: jgen
         integer :: imag,iaf
         integer, dimension(2,3) :: jaf
         real(DP) ::  a, b, c, ca, cb, cc
!      integer janti(2,3),kanti(2,3)
         real(DP), dimension(3) :: omove
         integer :: rotr2, irotk2
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!$$$      real*8  tabn(3,3),ta1n(3,48)
!$$$      integer  lra1n(3,3,48),ig01n(48)
         real(DP), dimension(3,3) :: tabn
         real(DP), dimension(3,49) :: ta1n
         integer, dimension(3,3,49) :: lra1n
         integer, dimension(48) :: ig01n

         if(ipri_spg >=  1 ) then
            if(printable) then
               write(nfout,*) '   '
               write(nfout,860) a,b,c,ca,cb,cc
             end if
         end if
  860    format(' a, b, c =',3f12.6/ &
         &       'ca,cb,cc =',3f12.6)

         idim = 3
         call tbspg(idim,ill,ngen,inv,igen,jgen,a,b,c,ca,cb,cc, &
        & omove,ipri_spg)
         ng1n = ng1
!$OMP parallel
         tabn(1:3,1:3) = tab(1:3,1:3)
         ta1n(1:3,1:ng1n) = ta1(1:3,1:ng1n)
         lra1n(1:3,1:3,1:ng1n) = lra1(1:3,1:3,1:ng1n)
         ig01n(1:ng1n) = ig01(1:ng1n)
!$OMP end parallel

!     ! Antiferro case
         if(imag==1) then
         lra1n(:,:,ng1n+1) = irot(:,:,iaf)
            do jj=1, 3
               ta1n(jj,ng1n+1) = real(jaf(1,jj),kind=DP)/real(jaf(2,jj),kind=DP)
            end do
         end if
      end subroutine setspg_default_n
! =====================================- added by K. Tagami ================= 12.0A
      subroutine setspg_default_n_kt( ill, ngen, inv, igen, jgen, &
     &                                imag, iaf, jaf, &
     &                                a, b, c, ca, cb, cc, &
     &                                tabn, ng1n, ta1n, lra1n, &
     &                                ipri_spg, ig01n, &
     &                                use_altv_rltv, altv, rltv, &
     &                                gen_name_in_carts )

         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : irot ! common/nspg0/
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & eeule2,eeulv2,ieule2,ieulv2,ig02,iv02         ! common/nspg2 /
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
!$ use omp_lib
         implicit none
         integer :: ill, ngen, inv, ipri_spg
         integer, dimension(ngen) :: igen
         integer, dimension(2,3,ngen) :: jgen
         integer :: idim, imag, iaf, ng1n
         integer :: jj
         integer, dimension(2,3) ::  jaf
         real(DP) ::  a, b, c, ca, cb, cc
!c      integer janti(2,3),kanti(2,3)
         real(DP), dimension(3) :: omove
! ---
         logical :: use_altv_rltv, gen_name_in_carts
         real(DP), dimension(3,3) :: altv, rltv
! ---
         integer :: irotr2, irotk2
!      real*8  ta2(3),tb2(3),tbv(3),tc2(3),tcv(3)
!      integer lrb2(3,3),lrbv(3,3)
!      real*8  rrb2(3,3),rrbv(3,3)
!      real*8  rtc2(3,3),rtcv(3,3)
!$$$      real*8  tabn(3,3),ta1n(3,48)
!$$$      integer  lra1n(3,3,48),ig01n(48)
         real(DP), dimension(3,3) :: tabn
         real(DP), dimension(3,49) :: ta1n
         integer, dimension(3,3,49) :: lra1n
         integer, dimension(48) :: ig01n

         if(ipri_spg  >= 1 ) then
            if(printable) then
               write(6,*) '   '
               write(6,860) a,b,c,ca,cb,cc
            end if
         end if
  860    format(' a, b, c =',3f12.6/ &
         &       'ca,cb,cc =',3f12.6)

         idim = 3
! ----------
          call tbspg_kt( idim, ill, ngen, inv, igen, jgen, &
        &               a, b, c, ca, cb, cc, omove, ipri_spg, &
        &               use_altv_rltv, altv, rltv, .false., &
        &               gen_name_in_carts )
! ----------

         ng1n = ng1
!$OMP parallel
         tabn(1:3,1:3) = tab(1:3,1:3)
         ta1n(1:3,1:ng1n) = ta1(1:3,1:ng1n)
         lra1n(1:3,1:3,1:ng1n) = lra1(1:3,1:3,1:ng1n)
         ig01n(1:ng1n) = ig01(1:ng1n)
!$OMP end parallel

!       Antiferro case
         if(imag==1) then
            lra1n(:,:,ng1n+1) = irot(:,:,iaf)
            do jj=1,3
               ta1n(jj,ng1n+1) = real(jaf(1,jj),kind=DP)/real(jaf(2,jj),kind=DP)
            end do
         end if
      end subroutine setspg_default_n_kt
! ====================================================================== 12.0A
      subroutine getspgtab(tabn)
          use m_Const_Parameters, only :  DP
          use m_Spg_Plus_Tetra_Common_Blocks, only : tab
         implicit none
         real(DP), intent(out), dimension(3,3) :: tabn

         tabn(1:3,1:3) = tab(1:3,1:3)

      end subroutine getspgtab
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine tbspg(idim,ill,ngen,inv,igen,jgen, &
     &   a,b,c,ca,cb,cc, omove,ipri_spg)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         use m_spg_plus_Tetra_Common_Blocks, only : tbspg_called
!$ use omp_lib
         implicit none
         integer, intent(in) :: idim,ill, ngen, inv
         integer, intent(in), dimension(ngen) :: igen
         integer, intent(in), dimension(2,3,ngen) :: jgen
         integer, intent(in) :: ipri_spg
         integer :: i, j, k, jf, jpr, jpr3, ind
         real(DP) :: a, b, c, ca, cb, cc
         real(DP) :: one, eps, x
         real(DP), dimension(3,3,48) ::  ra1, sa1, rb1, sb1
         real(DP), dimension(3) :: omove
         character(len=5), dimension(48) :: schoen
         integer ::  ncounter
         data    ncounter/0/
         save    ncounter

! ==================================== added by K. Tagami ========== 12.0A
         logical gen_name_in_carts

         real(DP), dimension(3,3) :: dummy1, dummy2

         dummy1 = 0.0d0
         dummy2 = 0.0d0

         gen_name_in_carts = .false.
! ================================================================== 12.0A

         if(tbspg_called) return

         ncounter = ncounter + 1

         one = 1.0d0
         il = ill

         if(idim == 2) then
            if(il /= 0 .and. il /= 1 .and. il /= 4) then
               if(printable) then
                  write(nfout,*) ' il=',il,' : stop in sub.tbspg.'
               end if
               stop '=tbspg (il, 2D)='
            end if

            if(il >= 1) then
               do i=1, ngen
                  if((igen(i) >= 5 .and. igen(i) <= 12) .or. &
                 &  (igen(i) >= 15 .and. igen(i) <= 20) .or. &
                 &  (igen(i) >= 22 .and. igen(i) <= 23) .or. &
                 &  (igen(i) >= 29 .and. igen(i) <= 36) .or. &
                 &  (igen(i) >= 39 .and. igen(i) <= 44) .or. &
                 &  (igen(i) >=.46 .and. igen(i) <= 47)) then
                     if(printable ) then
                        write(nfout,*) &
                       & ' igen(',i,')=',igen(i),'  in sub.tbspg'
                     end if
                     stop '=tbspg (igen, 2D)='
                  end if
               end do
            end if

            do i=1, ngen
               if(jgen(1,3,i) /= 0) then
                  if(printable) then
                     write(nfout,*) &
                    & ' jgen(1,3,',i,')=',jgen(1,3,i),'  in sub.tbspg'
                  end if
                  stop '=tbspg (jgen, 2D)='
               end if
            end do
         end if

         jf=6
         if(ipri_spg >= 1 .and. ncounter  <= 1) then
            jpr = 0
         else
            jpr = -1
         end if

! ======================================== modified by K. Tagami ======== 12.0A
!      call nspace(jf,jpr,il,ngen,inv,igen,jgen,
!     &            ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1,
!     &            omove,euler1,inver1)
         call nspace(jf,jpr,il,ngen,inv,igen,jgen, &
             &            ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1, &
             &            omove,euler1,inver1, &
             &            .false., tac, tca, tab, tba, &
             &             gen_name_in_carts )
! ======================================================================= 12.0A

         if(ipri_spg >= 1 .and. ncounter <= 1) then
            if(printable) then
                write(6,120) (omove(i),i=1,3)
            end if
         end if
  120    format(/'omove=(',3f9.6,' )')

         call nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
         if(ipri_spg >= 2 .and. ncounter  <= 1) then
            if(printable) write(nfout,140) a,b,c,ca,cb,cc
         end if
  140    format(/'  a, b, c=',3f12.6/' ca,cb,cc=',3f12.6)

         if(ipri_spg  >= 1 .and. ncounter <= 1) then
            jpr3 = 0
         else
            jpr3 = -1
         end if

! ======================================= modified by K. Tagami ========== 12.0A
!      call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc,
!     &                                   grc,gkc,gra,gka,grb,gkb)
         call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc, &
             &       grc,gkc,gra,gka,grb,gkb, &
             &       .false., dummy1, dummy2 )
! ======================================================================== 12.0A
         call nsgrpb(ng1,tab,tba,ta1,ra1,sa1,tb1,rb1,sb1)

         eps = 1.0d-5

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(ra1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind=1
                     if(printable) then
                        write(nfout,*) ' ra1(',i,j,k,') =',ra1(i,j,k)
                        write(nfout,*) &
                      & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ='
                  end if
                  lra1(i,j,k) = nint(ra1(i,j,k))
               end do
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(sa1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind=1
                     if(printable) then
                        write(nfout,*) ' sa1(',i,j,k,') =',sa1(i,j,k)
                        write(nfout,*) &
                       & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ='
                  end if
                  lsa1(i,j,k) = nint(sa1(i,j,k))
               enddo
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(rb1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind=1
                     if(printable) then
                        write(nfout,*) ' rb1(',i,j,k,') =',rb1(i,j,k)
                        write(nfout,*) &
                       & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ='
                  end if
                  lrb1(i,j,k) = nint(rb1(i,j,k))
               end do
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(sb1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind=1
                     if(printable) then
                        write(nfout,*) ' sb1(',i,j,k,') =',sb1(i,j,k)
                        write(6,*) &
                       & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ='
                  end if
                  lsb1(i,j,k)=nint(sb1(i,j,k))
               end do
            end do
         end do

         tbspg_called = .true.
      end  subroutine tbspg
! ============================================ added by K. Tagami ========= 12.0
      subroutine tbspg_kt(idim,ill,ngen,inv,igen,jgen, &
      &                 a,b,c,ca,cb,cc, omove,ipri_spg, &
      &                 use_altv_rltv, altv, rltv, &
      &                 use_trs, gen_name_in_carts )
         use m_Const_Parameters, only :  DP
        use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & tca,tac,tab,tba,tcb,tbc                       ! common/nspg03
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & grc,gkc,gra,gka,grb,gkb                       ! common/nspg04
         use m_Spg_plus_Tetra_Common_Blocks, only : &
        & id0,il,ng0,ng1,ig01,im1,iv1,ta1,lra1,lsa1,tb1,lrb1,lsb1 ! common/nspg06
         use m_Spg_plus_Tetra_Common_Blocks, only : euler1,inver1 ! common/nspg07
         use m_Spg_plus_Tetra_Common_Blocks, only : tbspg_kt_called
         implicit none
         integer, intent(in)  :: idim , ill, ngen, inv, ipri_spg
         integer, intent(in), dimension(3) :: igen
         integer, intent(in), dimension(2,3,3) :: jgen
         integer :: i, j, k, jf, jpr, jpr3, ind
         real(DP) :: a, b, c, ca, cb, cc
         real(DP) :: one, eps, x
         real(DP), dimension(3,3,48) :: ra1, sa1, rb1, sb1, rc1, sc1
         real(DP), dimension(3) :: omove
         character(len=5), dimension(48) :: schoen
         integer ::  ncounter
         data    ncounter/0/
         save ncounter
! -------------------
         logical :: use_altv_rltv, use_trs, gen_name_in_carts
         real(DP), dimension(3,3) :: rltv, altv
! local
! -------------------

         if(tbspg_kt_called) return

         ncounter = ncounter + 1

         il = ill
         if(idim == 2) then
            if(il /= 0 .and. il /= 1 .and. il /= 4) then
               if(printable) then
                  write(nfout,*) ' il=',il,' : stop in sub.tbspg.'
               end if
               stop '=tbspg (il, 2D)='
            end if
            if(il >= 1) then
               do i=1, ngen
                  if((igen(i) >=  5 .and. igen(i) <= 12) .or. &
                 &   (igen(i) >= 15 .and. igen(i) <= 20) .or. &
                 &   (igen(i) >= 22 .and. igen(i) <= 23) .or. &
                 &   (igen(i) >= 29 .and. igen(i) <= 36) .or. &
                 &   (igen(i) >- 39 .and. igen(i) <= 44) .or. &
                 &   (igen(i) >= 46 .and. igen(i) <= 47)) then
                     if(printable) then
                        write(nfout,*) &
                       &' igen(',i,')=',igen(i),'  in sub.tbspg'
                     end if
                     stop '=tbspg (igen, 2D)='
                  end if
               end do
            end if
            do i=1, ngen
               if(jgen(1,3,i) /= 0) then
                  if(printable) then
                     write(nfout,*) &
                   & ' jgen(1,3,',i,')=',jgen(1,3,i),'  in sub.tbspg'
                  end if
                  stop '=tbspg (jgen, 2D)='
               end if
            end do
         end if

         jf = 6
         if(ipri_spg  >= 1 .and. ncounter  <= 1) then
            jpr = 0
         else
            jpr = -1
         end if

         call nspace(jf,jpr,il,ngen,inv,igen,jgen, &
        &             ng0,schoen,ng1,ig01,ta1,ra1,sa1,im1,iv1, &
        &             omove,euler1,inver1, &
        &             use_trs, tac, tca, tab, tba, &
        &             gen_name_in_carts )

         if(ipri_spg  >= 1 .and. ncounter <= 1) then
            if(printable) write(nfout,120) (omove(i),i=1,3)
         end if
 120     format(/'omove=(',3f9.6,' )')

! --
         if( .not. use_altv_rltv ) then
            call nslatz(il,ng1,ig01,a,b,c,ca,cb,cc)
         endif

         if(ipri_spg >= 2 .and. ncounter  <= 1) then
            if(printable)  write(nfout,140) a,b,c,ca,cb,cc
         end if
 140     format(/'  a, b, c=',3f12.6/' ca,cb,cc=',3f12.6)

         if(ipri_spg  >= 1 .and. ncounter  <= 1) then
            jpr3 = 0
         else
            jpr3 = -1
         end if

! ------------------------------------------
         call nslat3(jpr3,il,a,b,c,ca,cb,cc,tca,tac,tab,tba,tcb,tbc, &
        &           grc,gkc,gra,gka,grb,gkb, &
        &           use_altv_rltv, altv, rltv )
         if ( gen_name_in_carts ) then
            rc1 = ra1
            sc1 = sa1
            call nsgrpa( ng1, tac, tca, rc1, sc1, ra1,sa1 )
         endif

         call nsgrpb_kt(ng1,tab,tba,ta1,ra1,sa1,tb1,rb1,sb1)

         eps = 1.0d-5

         if ( gen_name_in_carts ) then
            call check_if_unit_matirx( ng1, rc1, sc1, 'rc1 or sc1' )
         endif
         call check_if_unit_matirx( ng1, ra1, sa1, 'ra1 or sa1' )
         call check_if_unit_matirx( ng1, rb1, sb1, 'rb1 or sb1' )
!  -------------------------------

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(ra1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind = 1
                     if(printable) then
                        write(nfout,*) ' ra1(',i,j,k,') =',ra1(i,j,k)
                        write(nfout,*) &
                       & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ph1='
                  end if
                  lra1(i,j,k)=nint(ra1(i,j,k))
               end do
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(sa1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind = 1
                     if(printable) then
                        write(nfout,*) ' sa1(',i,j,k,') =',sa1(i,j,k)
                        write(nfout,*) &
                      & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ph2='
                  end if
                  lsa1(i,j,k)=nint(sa1(i,j,k))
               end do
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(rb1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind = 1
                     if(printable) then
                        write(nfout,*) ' rb1(',i,j,k,') =',rb1(i,j,k)
                        write(nfout,*) &
                      & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ph3='
                  end if
                  lrb1(i,j,k)=nint(rb1(i,j,k))
               end do
            end do
         end do

         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  x = dabs(sb1(i,j,k))+0.5d0*eps
                  if(dmod(x,1.0d0) > eps) then
                     ind = 1
                     if(printable) then
                        write(nfout,*) ' sb1(',i,j,k,') =',sb1(i,j,k)
                        write(nfout,*) &
                       & ' stop in sub.tbspg (noninteger matrix el.)'
                     end if
                     stop '=tbspg ph4='
                 end if
                 lsb1(i,j,k) = nint(sb1(i,j,k))
               end do
            end do
         end do
         tbspg_kt_called = .true.
      end subroutine tbspg_kt

      subroutine check_if_unit_matirx( ng1, rmat_in, smat_in, comment1 )
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         implicit none
! inout
         integer :: ng1
         real(DP), dimension(3,3,ng1)  :: rmat_in, smat_in
         character(len=10) :: comment1
! local
         integer :: i, j, k, l
         real(DP)  ::  c1, c2, eps
! init
         eps = 1.0D-5
! begin
!$OMP parallel private(c1,c2)
!$OMP do
         do k=1, ng1
            do i=1, 3
               do j=1, 3
                  c1 = 0.0d0
                  c2 = 0.0d0
                  do l=1, 3
                     c1 = c1 + rmat_in(i,l,k) *smat_in(l,j,k)
                     c2 = c2 + rmat_in(l,i,k) *smat_in(j,l,k)
                  end do
                  if (  ( i == j .and. dabs(c1-1.0d0) > eps ) .or. &
                 &      ( i /= j .and. dabs(c1) > eps ) ) then
                     if(printable) then
                        write(nfout,*) trim(comment1) &
                       &  // ' is not properly obtained.'
                        write(nfout,*) 'i j c1 = ', i, j, c1
                     endif
                  end if
                  if (  ( i == j .and. dabs(c2-1.0d0) > eps ) .or. &
                      & ( i /= j .and. abs(c2) > eps ) ) then
                     if(printable) then
                        write(nfout,*) trim(comment1) &
                      &  // ' is not properly obtained.'
                       write(6,*) 'i j c2 = ', i, j, c2
                     endif
                  end if
               end do
            end do
         end do
!$OMP end do
!$OMP end parallel
      end subroutine check_if_unit_matirx
! ======================================================================== 12.0A


!===*====1====*====2====*====3====*====4====*====5====*====6====*====7
!
!11  sub.tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
!
!#12  input:      il, ng1, ir1234, ig01, iv0, im0, jg1
!#12  output: commons
!#13  noexternal:
!
!21  to be compatible with the tspace package

!#31  1990.01.09.:  n. hamada, a. yanase and k. terakura
!
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine tspaca(il,ng1,ir1234,ig01,iv0,im0,jg1)
         use m_Spg_plus_Tetra_Common_Blocks, only : it, im, iv, il9, ng9, ig0, jv
!$ use omp_lib
         implicit none
         integer, intent(in) :: il, ng1
         integer :: j
         integer, intent(in), dimension(3,48) :: ir1234
         integer, intent(in), dimension(48,48) :: im0
         integer, intent(in), dimension(48) :: iv0, ig01
         integer, intent(in), dimension(2,3,48) :: jg1

!     for being compatible with tspace program package

         it(1:3,1:48)=ir1234(1:3,1:48)
!$$$      do 92 j=1,48
!$$$      ig0(j)=ig01(j)
         do  j=1, 48
            if(j <= ng1) then
               ig0(j) = ig01(j)
            else
               ig0(j) = j
            end if
            iv(j) = iv0(j)
         end do
!$OMP parallel
         im(1:48,1:48)=im0(1:48,1:48)
         il9 = il
         ng9 = ng1
         jv(1:2,1:3,1:ng1)=jg1(1:2,1:3,1:ng1)
!$OMP end parallel

!      open(unit= 2, file='a.spg', status='unknown')
!      write( 2,200) il,ng1
!      write( 2,200) ((it(i,j),i=1, 3),j=1,48)
!      write( 2,200) ((im(i,j),i=1,48),j=1,48)
!      write( 2,200) (iv  (j),j=1,48)
!      write( 2,200) (ig01(j),j=1,ng1)
!      write( 2,200) (((jg1(i,j,k),i=1,2),j=1,3),k=1,ng1)
!      close(unit= 2)
!  200 format(24i3)
      end subroutine tspaca

      subroutine wtetra(nxx,nyy,nzz,np0,np2,ip20,iwt,ip2cub,ip2cub_wk)
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
!$ use omp_lib
         implicit none
         integer, intent(in) :: nxx, nyy, nzz, np0, np2
         integer, intent(in), dimension(np0) :: ip20
         integer, intent(out), dimension(np2) :: iwt
         integer, intent(out), dimension(nxx*nyy*nzz) :: ip2cub
         integer, dimension(9,nxx*nyy*nzz) :: ip2cub_wk
         integer ::  npx,npy,npz, np,ntet,ncub
         integer :: ip,ip0,ip1,ip2,ix,iy,iz,kx,ky,kz,it,iq,m,icub,jcub
         integer :: flag
         integer, dimension(2,2,2)  :: iecub
         integer, dimension(8) :: iec
         integer, dimension(4) :: iet, ieb
         integer, dimension(6,2) :: iqmat
         integer, allocatable, dimension(:,:,:) :: ni, ni0
         real(DP), dimension(2,2,2) ::  ecub
         real(DP), dimension(8) :: ec

         equivalence(ec(1),ecub(1,1,1))
         equivalence(iec(1),iecub(1,1,1))
         data iqmat/2,2,5,3,3,5, 4,6,6,4,7,7/

         npx = nxx+1
         npy = nyy+1
         npz = nzz+1
         np = npx*npy*npz
         ncub = nxx*nyy*nzz
         ntet = 6*ncub

         iwt(1:np2) = 0
         allocate(ni(0:nxx-1,0:nyy-1, 0:nzz-1))
         allocate(ni0(2,2,2))


!$OMP simd
         do iz=0, nzz-1
            do iy=0, nyy-1
               do ix=0, nxx-1
                  ni(ix,iy,iz) = npx*(npy*iz+iy)+ix
               end do
            end do
         end do
!$OMP end simd

!$OMP simd
         do kz = 1, 2
            do ky = 1, 2
               do kx = 1,2
                  ni0(kx,ky,kz) = npx*(npy*(kz-1)+ky-1)+kx
               end do
            end do
         end do
!$OMP end simd

         do iz=0, nzz-1
            do iy=0, nyy-1
               do ix=0, nxx-1
!                 ni = npx*(npy*iz+iy)+ix
                  do kz=1,2
                     do ky=1, 2
                        do kx=1, 2
!                          ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
!                           if(ip0 > np0) then
!                              if(printable) write(nfout,*) ' ip0, np0  ',ip0,np0
!                              stop ' wtetra -- ip0.ne.np0'
!                           endif
!                          iecub(kx,ky,kz) = ip0
                           iecub(kx,ky,kz) =  ni(ix,iy,iz) + ni0(kx,ky,kz)
                        end do
                     end do
                  end do
                  iet(1) = iec(1)
                  iet(4) = iec(8)
                  do it=1, 6
                     do ip=1, 2
                        iq = iqmat(it,ip)
                        iet(ip+1) = iec(iq)
                     end do
                     do m=1, 4
                        ieb(m) = iet(m)
                        if(ip20(ieb(m)) > np2) then
                           if(printable) then
                              write(nfout,*) ' ip20, np2  ',ip20(ieb(m)),np2
                           end if
                           stop ' wtetra -- ip20.ne.np2'
                        endif
                        iwt(ip20(ieb(m))) = iwt(ip20(ieb(m))) + 1
                     end do
                  end do
               end do
            end do
         end do

         icub = 0
         do iz=0, nzz-1
            do iy=0, nyy-1
               do ix=0, nxx-1
                  icub = icub+1
!                 ni = npx*(npy*iz+iy)+ix
                  ip = 0
                  do kz=1,2
                     do ky=1, 2
                        do kx=1, 2
                           ip = ip+1
!                          ip0 = ni+npx*(npy*(kz-1)+ky-1)+kx
                           ip0 = ni(ix,iy,iz)+ ni0(kx,ky,kz)
                           iec(ip) = ip20(ip0)
                        end do
                     end do
                  end do
                  do ip1=2, 6
                     do ip2=ip1+1, 7
                        if(iec(ip1) > iec(ip2)) then
                           ip = iec(ip1)
                           iec(ip1) = iec(ip2)
                           iec(ip2) = ip
                        endif
                     end do
                  end do
                  if(iec(1) < iec(8)) then
                     ip2cub_wk(1,icub) = iec(1)
                     ip2cub_wk(8,icub) = iec(8)
                  else
                     ip2cub_wk(8,icub) = iec(1)
                      ip2cub_wk(1,icub)=iec(8)
                  end if
                  ip0 = iec(1)+iec(8)
                  do ip=2, 7
                     ip2cub_wk(ip,icub) = iec(ip)
                     ip0 = ip0+iec(ip)
                  end do
                  ip2cub(icub) = icub
                  ip2cub_wk(9,icub) = ip0
               end do
            end do
         end do

         deallocate(ni)
         deallocate(ni0)

         ncub=nxx*nyy*nzz

!$OMP parallel
!$OMP do
         do 80 icub=2, ncub
            do 81 jcub=1, icub-1
               if(ip2cub_wk(9,icub) == ip2cub_wk(9,jcub)) then
                  do 82  ip=1,8
                     if(ip2cub_wk(ip,icub) /= ip2cub_wk(ip,jcub)) then
                        goto 81
                     end if
     82           end do
                  ip2cub(icub)=jcub
                  goto 80
               end if
     81       end do
     80   end do
!$OMP end do
!$OMP end parallel

         jcub = 0
!$OMP parallel reduction(+:jcub)
!$OMP do
         do icub=1, ncub
!$$$      write(6,'(2x,i5,''  -->'',i5)') icub,ip2cub(icub)
            if(icub == ip2cub(icub)) jcub=jcub+1
         end do
!$OMP end do
!$OMP end parallel
!$$$      write(6,'(2x,''number of cube which should be calculated'',i5
!$$$     &     ,"<<wtetra>>")')  jcub
      end subroutine wtetra
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstt4i(ipri,idim,newindows,nxyz_tetra,np2,lmnp2e,neig,eeig, &
        &     ip20,np0,lmnp2c,lmneig,mtetra,nttra,deltae, &
        &     nep,dos,dosin,cdos,cind)
!
!     nstt0i, nstt1i, nstt2i are nstt3i are merged into
!  this subroutine nstt4i
!               by T. Yamasaki, Aug 2007
!
         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Spg_plus_Tetra_Common_Blocks, only : factor_1_12
         implicit none
         integer, intent(in) :: ipri
         integer, intent(in) :: idim, newindows, np2, lmnp2e, neig,&
       & np0, lmnp2c, lmneig, mtetra, nep
         integer, intent(in), dimension(3) :: nxyz_tetra
         integer, intent(in), dimension(np0) :: ip20
         integer, intent(in), dimension(mtetra,4) :: nttra
         integer, dimension(4) :: ieb
         integer :: nxx, nyy, nzz, npx, npy, npz, np, ncub, ntet
         integer :: ib, ip, iloop, idos, iv, ie, ns, ns1, ns2, ns3, ne, ne1, ne2, ne3
         integer :: k2, ieig, i, n
         real(DP), intent(in) :: deltae
         real(DP), intent(in), dimension(lmnp2e,lmneig) :: eeig
         real(DP), intent(out), dimension(0:nep,4) :: dos, dosin
         real(DP), intent(out), dimension(0:newindows,lmnp2c,lmneig) :: cdos, cind
         real(DP), dimension(4) :: eb
         integer :: ncounter,ncounter1,ncounter2,ncounter3
         real(DP) :: eps, es, ee, e, e1, e2, e3, e4, esum, tdos, d21, d31, d41, d32,&
       &  d42, d43, d1, d2, d3, d4, x, xx, y, yy1, yy2, yy3, &
       &  d43yy2, d21yy1, d31yy1, d41yy1, d21yy2, d31yy2, d41yy2, d42yy2, &
       & x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, &
       & d31d31d32d41, d31d41, d42d32, d32d32d31
         real(DP) :: c1, c2
         real(DP) ::  etime1,etime2,etime3, wct_now, wct_start,etime4, &
       & etime5, etime0

         eps = 1.0d-5

         nxx = nxyz_tetra(1); nyy = nxyz_tetra(2); nzz = nxyz_tetra(3)
         npx = nxx+1
         npy = nyy+1
         npz = nzz+1
         np = npx*npy*npz
         ncub = nxx*nyy*nzz
         ntet = 6*ncub

         if(ipri>=2) then
            if(printable) then
               write(nfout,'(" lmnp2e, lmneig, mtetra, newindows, lmnp2c," &
              &,"nep = ",6i8)') lmnp2e,lmneig,mtetra,newindows,lmnp2e,nep
               write(nfout,'(" np = ",i5," np2 = ",i5)')np,np2
            end if
         end if
         es =  1.0d30
         ee = -1.0d30
         do ib = 1, neig
            do ip = 1, lmnp2e
               if(eeig(ip,ib) < es) es = eeig(ip,ib)
               if(eeig(ip,ib) > ee) ee = eeig(ip,ib)
            end do
         end do
         if(ipri>=2) then
            if(printable) write(nfout,'(" es,ee = ",2d12.4)') es,ee
         end if

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
!$$$         write(6,'(" iv = ",i8)') iv
!$$$         write(6,'(" -- nttra  = ",4i8)') (nttra(iv,ie),ie=1,4)
            do ib = 1, neig
               call gettod(wct_start)
!$$$               if(nttra(iv,ie).gt.np0 .or. nttra(iv,ie).le.0) then
!$$$                  write(6,'(" nttra(",i5,",",i5,") = ",i5)')
!$$$     &                 iv,ie,nttra(iv,ie)
!$$$                  stop 'nttra illegal'
!$$$               end if
               do ie = 1, 4
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
               ns = int((e1 - es)/deltae)+1
               ne = int((e4 - es)/deltae)+1
               if(ns < 0) ns = 0
               if(ne > nep) ne = nep
!$$$            write(6,'(" iv,ib = ",2i5," e = ",4f8.4," ns,ne=",2i8)')
!$$$     &           iv,ib, e1,e2,e3,e4,ns,ne

               d21 = e2-e1
               d31 = e3-e1
               d41 = e4-e1
               d32 = e3-e2
               d42 = e4-e2
               d43 = e4-e3

               iloop = iloop + (ne-ns+1)
               ne1 = int((e2-es)/deltae)
               if((es+deltae*(real(ne1,kind=DP)-0.5d0))>e2) ne1 = ne1-1
               ns2 = ne1+1
               ne2 = int((e3-es)/deltae)
               if((es+deltae*(real(ne2,kind=DP)-0.5d0))>e3) ne2 = ne2-1
               ns3 = ne2+1

               if(ipri >=2 ) then
                  if(printable) then
                     write(nfout,&
                    & '(" ns, ne1, ns2, ne2, ns3, ne = ",6i9)') &
                    &     ns ,ne1, ns2, ne2, ns3, ne
                  end if
               end if

!$$$            do idos = ns, ne
               yy1 = d41*d31*d21
               d21yy1 = d21*yy1
               d31yy1 = d31*yy1
               d41yy1 = d41*yy1
               call gettod(wct_now)
               etime0 = etime0 + (wct_now-wct_start)*1.d-6

               call gettod(wct_start)
               do idos = ns, ne1
                  e = es+deltae*(real(idos,kind=DP)-0.5d0)
                  d1 = e-e1
                  d4 = e4-e
                  ncounter1 = ncounter1+1
                  d2 = e2-e
                  d3 = e3-e
!$$$                  yy=d41*d31*d21
                  x = d2/d21+d3/d31+d4/d41
                  y = (d1*d1)/yy1
                  dos(idos,1) = x*y
                  dosin(idos,1) = 0.25d0*d1*y*(x+1.0d0)
                  xx = d1*d1*d1
                  x = xx/d21yy1
                  dos(idos,2) = x
                  dosin(idos,2) = 0.25d0*d1*x
                  x = xx/d31yy1
                  dos(idos,3) = x
                  dosin(idos,3) = 0.25d0*d1*x
                  x = xx/d41yy1
                  dos(idos,4) = x
                  dosin(idos,4) = 0.25d0*d1*x
               end do
               call gettod(wct_now)
               etime1 = etime1 + (wct_now-wct_start)*1.d-6

               yy2 = d41*d42*d43
               d41yy2 = d41*yy2
               d42yy2 = d42*yy2
               d43yy2 = d43*yy2

               call gettod(wct_start)
               do idos = ns3, ne
                  e = es+deltae*(real(idos,kind=DP)-0.5d0)
                  d1 = e-e1
                  d4 = e4-e
                  ncounter3 = ncounter3+1
                  d2 = e-e2
                  d3 = e-e3
                  xx = d4*d4*d4
!$$$                  yy=d41*d42*d43
                  x = xx/d41yy2
                  dos(idos,1) = x
                  dosin(idos,1) = 0.25d0*(1.0d0-d4*x)
                  x = xx/d42yy2
                  dos(idos,2) = x
                  dosin(idos,2) = 0.25d0*(1.0d0-d4*x)
                  x = xx/d43yy2
                  dos(idos,3) = x
                  dosin(idos,3) =0.25d0*(1.0d0-d4*x)
                  x = d3/d43+d2/d42+d1/d41
                  y = (d4*d4)/yy2
                  dos(idos,4) = x*y
                  dosin(idos,4) =0.25d0*(1.0d0-d4*y*(x+1.0d0))
               end do
               call gettod(wct_now)
               etime3 = etime3 + (wct_now-wct_start)*1.d-6
               yy3 = 1.0d0/(d31*d42)+1.0d0/(d41*d32)
               d31d31d32d41 = d31*d31*d32*d41
               d31d41       = d31*d41
               d42d32       = d42*d32
               d32d32d31    = d32*d32*d31

               call gettod(wct_start)
               do idos=ns2, ne2
                  e = es+deltae*(real(idos,kind=DP)-0.5d0)
                  d1 = e-e1
                  d4 = e4-e
                  ncounter2 = ncounter2+1
                  d2 = e-e2
                  d3 = e3-e
!$$$                  y=1.0/(d31*d42)+1.0/(d41*d32)
                  x1 = (d3*d3)/(d31*d31*d32)*(d2/d42+d1/d41)
                  x2 = (d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
                  x3 = (d3*d4*d1)/(d31d41)*yy3
                  dos(idos,1)=0.5d0*(x1+x2+x3)
                  x1 = d2*d2*(d32*d2+3.0d0*d3*(d32+d3))*factor_1_12
                  y1 = x1
                  x1 = x1/(d31*d31*d32*d42)
                  x2 = d2*(d2*d2*(d31+3.0d0*d21)+ &
                &   3.0d0*d3*(d2*d3+d32*(3.0d0*d21+d1)))
                  x2 = x2*factor_1_12
                  y2 = x2
                  x2 = x2/d31d31d32d41
                  x3 = d2*d2*(d42*d2+3.0d0*d4*(d42+d4))*factor_1_12
                  y3 = x3
                  x3 = x3/(d41*d41*d42d32)
                  x4 = d2*(d2*d2*(d41+3.0d0*d21)+ &
                &   3.0d0*d4*(d2*d4+d42*(3.0d0*d21+d1)))
                  x4 = x4*factor_1_12
                  y4 = x4
                  x4 = x4/(d41*d41*d42*d31)
                  x5 = 0.5d0*d2*d3*d4*(d1+d21)
                  x5 = x5+d2*d2*(2.0d0*d21*(d3+d42)+ &
                &   (d1+d21)*(2.0d0*d3+d4+d42))*factor_1_12
                  x5 = x5*yy3/(d31d41)
                  x6 = 0.25d0*d21*d21*(d42/d41+d32/d31+1.0d0)/(d41*d31)
                  dosin(idos,1) = 0.5d0*(x1+x2+x3+x4+x5)+x6
                  x1 = (d3*d3)/(d32d32d31)*(d2/d42+d1/d41)
                  x2 = (d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
                  x3 = (d3*d4*d2)/(d42d32)*yy3
                  dos(idos,2) = 0.5d0*(x1+x2+x3)
                  x1 = y1/(d32d32d31*d42)
                  x2 = y2/(d32d32d31*d41)
                  x3 = y3/(d42*d42*d41*d32)
                  x4 = y4/(d42*d42*d41*d31)
                  x5 = d2*d2*(d3*(d42+3.0d0*d4)+d32*(d42+d4))*factor_1_12
                  x5 = x5*yy3/(d42d32)
                  x6 =0.25d0*d21/d31*d21/d41
                  dosin(idos,2) =0.5d0*(x1+x2+x3+x4+x5)+x6
                  x1 = (d2*d2)/(d32*d32*d42)*(d3/d31+d4/d41)
                  x2 = (d1*d1)/(d31*d31d41)*(d3/d32+d4/d42)
                  x3 = (d1*d2*d3)/(d32*d31)*yy3
                  dos(idos,3) = 0.5d0*(x1+x2+x3)
                  x1 = d2*d2*d2*(3.0d0*d3+d32)*factor_1_12
                  y1 = x1
                  x1 = x1/(d32*d32*d42*d31)
                  x2 = d2*d2*d2*(3.0d0*d4+d42)*factor_1_12
                  y2 = x2
                  x2 = x2/(d32*d32*d42*d41)
                  x3 = d2*(d2*d31*(d2+3.0d0*d21)+3.0d0*d3*(d2*d2+3.0d0*d21*d1)+ &
                &   3.0d0*d21*d21*d32)*factor_1_12
                  y3 = x3
                  x3 = x3/(d31*d31d41*d32)
                  x4 = d2*(d2*d41*(d2+3.0d0*d21)+3.0d0*d4*(d2*d2+3.0d0*d21*d1)+3.0d0* &
                &   d21*d21*d42)*factor_1_12
                  y4 = x4
                  x4 = x4/(d31*d31d41*d42)
                  x5 = (d2*d2)*(d3*(d21+3.0d0*d1)+d32*(d21+d1))*factor_1_12
                  x5 = x5/(d32*d31)*yy3
                  x6 =0.25d0*(d21/d31*d21/d31*d21/d41)
                  dosin(idos,3)=0.5d0*(x1+x2+x3+x4+x5)+x6
                  x1 = (d2*d2)/(d42*d42d32)*(d3/d31+d4/d41)
                  x2 = (d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
                  x3 = (d1*d2*d4)/(d41*d42)*yy3
                  dos(idos,4)=0.5d0*(x1+x2+x3)
                  x1 = y1/(d42*d42d32*d31)
                  x2 = y2/(d42*d42d32*d41)
                  x3 = y3/(d41*d41*d31*d32)
                  x4 = y4/(d41*d41*d31*d42)
                  x5 = (d2*d2)*(d4*(d21+3.0d0*d1)+d42*(d21+d1))*factor_1_12
                  x5 = x5/(d41*d42)*yy3
                  x6 =0.25d0*(d21/d31*d21/d41*d21/d41)
                  dosin(idos,4) = 0.5d0*(x1+x2+x3+x4+x5)+x6
               end do
               call gettod(wct_now)
               etime2 = etime2 + (wct_now-wct_start)*1.d-6
               call gettod(wct_start)
               if(idim == -3) then
                  esum = e1+e2+e4+e4
                  do idos = ns, ne
                     tdos = 0.025d0*(dos(idos,1)+dos(idos,2) &
                    &   +dos(idos,3)+dos(idos,4))
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
         if(ipri>=2) then
            if(printable) then
               write(nfout,'(" before cdos, cind substitution")')
            end if
         end if

         cdos(0:newindows,1:lmnp2c,1:neig) = &
       & cdos(0:newindows,1:lmnp2c,1:neig)/real(ntet,kind=DP)
         cind(0:newindows,1:lmnp2c,1:neig) = &
       & cind(0:newindows,1:lmnp2c,1:neig)/real(ntet,kind=DP)

! ---- following lines are the later part of subroutine nstt3i --
!     take care of a weight on a degenerate state
         do k2=1, np2
            ieig=1
   40       continue
            n = 1
!        !!$do 42 i=1,20
            do i=1, neig
               if(ieig+i > neig) then
                  exit
               end if
               if(dabs(eeig(k2,ieig+i)-eeig(k2,ieig)) < eps) then
                  n = n+1
               else
                  exit
               end if
            end do

            do ie=0, ne
               c1 = 0.0d0
               c2 = 0.0d0
               do i=0, n-1
                  c1 = c1+cdos(ie,k2,ieig+i)
                  c2 = c2+cind(ie,k2,ieig+i)
               end do
               c1 = c1/real(n,kind=DP)
               c2 = c2/real(n,kind=DP)
               do i=0, n-1
                  cdos(ie,k2,ieig+i)=c1
                  cind(ie,k2,ieig+i)=c2
               end do
            end do
            ieig = ieig+n
            if(ieig < neig) go to 40
         end do
! -------------------------------
         call gettod(wct_now)
         if(ipri>=2) then
            etime5 = etime5 + (wct_now-wct_start)*1.d-6
            if(printable) then
               write(nfout,'(" !dos iloop    = ",i12)') iloop
               write(nfout,'(" !dos ncounter1 = ",i12)') ncounter1
               write(nfout,'(" !dos ncounter2 = ",i12)') ncounter2
               write(nfout,'(" !dos ncounter3 = ",i12)') ncounter3
               write(nfout,'(" !dos etime0    = ",f16.8)') etime0
               write(nfout,'(" !dos etim      = ",f16.8)') etime1
               write(nfout,'(" !dos etime2    = ",f16.8)') etime2
               write(nfout,'(" !dos etime3    = ",f16.8)') etime3
               write(nfout,'(" !dos etime4    = ",f16.8)') etime4
               write(nfout,'(" !dos etime5    = ",f16.8)') etime5
            end if
         end if
      end subroutine nstt4i
! --*----1----*----2----*----3----*----4----*----5----*----6----*----7
      subroutine nstt5i(ipri,idim,es,ee,newindows,nxyz_tetra, &
     &     np2,lmnp2e,neig, &
     &     eeig,ip20,np0,lmnp2c,lmneig,mtetra,nttra,deltae,dos_weight, &
     &     nep,dos,dosin,cdos,csumdos)
!
!     nstt0i, nstt1i, nstt2i and nstt3i are merged into
!  this subroutine nstt5i
!               by T. Yamasaki, Aug 2007

         use m_Const_Parameters, only :  DP
         use m_Files,            only :  nfout
         use m_Control_Parameters, only    : printable
         use m_Parallelization, only : npes,mype, MPI_CommGroup, &
       & mpi_double_precision, mpi_sum
         use m_Spg_plus_Tetra_Common_Blocks, only: factor_1_12
         implicit none
!         include 'mpif.h'
         integer, intent(in) :: ipri, idim, newindows
         integer, intent(in) :: neig, np0, np2, lmneig, lmnp2e, lmnp2c, mtetra, nep
         integer, intent(in), dimension(3) :: nxyz_tetra
         integer, intent(in), dimension(np0) :: ip20
         integer, intent(in), dimension(mtetra,4) :: nttra
         integer :: ip
         integer, dimension(4) :: ieb
         real(DP), intent(in) ::  es, ee, deltae
         real(DP), intent(in), dimension(lmnp2e,lmneig) :: eeig
         real(DP), intent(in), dimension(lmneig,lmnp2c) :: dos_weight
         real(DP), intent(out), dimension(0:nep,4) :: dos, dosin
         real(DP), intent(out), dimension(0:newindows) :: cdos, csumdos
!MPI
         real(DP), allocatable, dimension(:,:) :: dos_mpi,dosin_mpi
         real(DP), allocatable, dimension(:) :: cdos_mpi,csumdos_mpi
!MPI
         real(DP), dimension(4) ::  eb
         integer :: ncounter,ncounter1,ncounter2,ncounter3
         real(DP) ::  etime1,etime2,etime3, wct_now, wct_start,etime4,etime5 &
       & ,etime0
         integer :: ierr,itmp

         integer :: nxx, nyy, nzz, npx, npy, npz, np, ncub, ntet
         integer :: iloop, iv, ib, ie, idos
         integer :: ns, ns1, ns2, ns3, ne, ne1, ne2, ne3
         real(DP) :: eps, e, e1, e2, e3, e4, d1, d2, d3, d4, d21, d31, d32, d41, d42, dr2, d43, &
        & x, y, xx, yy, td21, yy2, yy3, w
         real(DP) :: d2d2, d2d2d2, d3d3, x1, x2, x3, x4, x5, x6, y1, y2, y3, y4
         real(DP) :: esum, tdos


         eps = 1.0d-4

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
         if(ipri /= 2) then
            if(printable) then
               write(nfout,'(" np = ",i5," np2 = ",i5," idim = ",i5)') &
             & np,np2,idim
            end if
         end if
!$$$      es =  dfloat(10)**30
!$$$      ee = -dfloat(10)**30
!$$$      do ib = 1, neig
!$$$      do ip = 1, lmnp2e
!$$$         if(eeig(ip,ib) < es) es = eeig(ip,ib)
!$$$         if(eeig(ip,ib) > ee) ee = eeig(ip,ib)
!$$$      end do
!$$$      end do
!$$$      es = es - 0.005
!$$$      ee = ee + 0.005
!$$$      write(6,'(" es,ee = ",2f10.6)') es,ee

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
            if(ipri >= 3) then
               if(printable) then
                  write(nfout,'(" iv = ",i8)') iv
                  write(nfout,'(" -- nttra  = ",4i8)') (nttra(iv,ie),ie=1,4)
               end if
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
               if(ipri >= 3) then
                  if(printable) then
                     write(nfout,'(" (iv,ib) = (",2i5,") e(1:4) = ",4f9.5, &
                    & "ieb(1:4) = ",4i3)') iv,ib,e1,e2,e3,e4 &
                    &,ieb(1),ieb(2),ieb(3),ieb(4)
                 end if
               end if
               if(e1 > ee) cycle
               ns = (e1 - es)/deltae
               if(es+deltae*ns <  e1) ns = ns + 1
               ne = (e4 - es)/deltae+1
               if(es+deltae*ne >= e4) ne = ne - 1
               if(es+deltae*ne >= e4) ne = ne - 1

               iloop = iloop + (ne-ns+1)
               ne1 = (e2-es)/deltae
               if(es+deltae*ne1 > e2) ne1 = ne1-1
               ns2 = ne1+1
               ns3 = (e3-es)/deltae+1
               if(es+deltae*ns3 < e3) ns3 = ns3+1
               ne2 = ns3-1

               if(ipri >= 3 .and. printable) then
                  write(nfout,&
                &  '(" ns, ne1 = ",2i8," e = ",f9.5," - ",f9.5)') &
                &   ns,ne1,es+deltae*ns,es+deltae*ne1
                  write(nfout,&
                &  '(" ns2,ne2 = ",2i8," e = ",f9.5," - ",f9.5)') &
                &   ns2,ne2,es+deltae*ns2,es+deltae*ne2
                 write(nfout,&
                &  '(" ns3,ne  = ",2i8," e = ",f9.5," - ",f9.5)') &
                &   ns3,ne,es+deltae*ns3,es+deltae*ne
                 if(es+deltae*ns  < e1) &
                & write(nfout,'(" !! es+deltae+ns  < e1")')
                 if(es+deltae*ne1 > e2) &
                & write(nfout,'(" !! es+deltae*ne1 > e2")')
                 if(es+deltae*ns2 < e2) &
                & write(nfout,'(" !! es+deltae*ns2 < e2")')
                 if(es+deltae*ne2 > e3) &
                & write(nfout,'(" !! es+deltae*ne2 > e3")')
                 if(es+deltae*ns3 < e3) &
                & write(nfout,'(" !! es+deltae*ns3 < e3")')
                 if(es+deltae*ne > e4) &
                & write(nfout,'(" !! es+deltae*ne  > e4")')
               end if

               call gettod(wct_now)
               etime0 = etime0 + (wct_now-wct_start)*1.d-6

               call gettod(wct_start)

               d21 = e2-e1
               d31 = e3-e1
               d41 = e4-e1
               d32 = e3-e2
               d42 = e4-e2
               d43 = e4-e3
               ncounter1 = ncounter1+(ne1-ns+1)
               do idos = ns, ne1
                  e = es + deltae*real(idos,kind=DP)
                  d1 = e-e1
                  d4 = e4-e
                  d2 = e2-e
                  d3 = e3-e
                  yy = d41*d31*d21
                  x = d2/d21+d3/d31+d4/d41
                  y = (d1*d1)/yy
                  dos(idos,1) = x*y
                  dosin(idos,1) = 0.25d0*d1*y*(x+1.0d0)
                  xx = d1*d1*d1
                  x = xx/(d21*yy)
                  dos(idos,2) = x
                  dosin(idos,2) = 0.25d0*d1*x
                  x = xx/(d31*yy)
                  dos(idos,3) = x
                  dosin(idos,3) = 0.25d0*d1*x
                  x =xx/(d41*yy)
                  dos(idos,4) = x
                  dosin(idos,4) =0.25d0*d1*x
               end do
               call gettod(wct_now)
               etime1 = etime1 + (wct_now-wct_start)*1.d-6

               call gettod(wct_start)
               ncounter3 = ncounter3+(ne-ns3+1)
               do idos = ns3,ne
                  e = es + deltae*idos
                  d1 = e-e1
                  d4 = e4-e
                  d2 = e-e2
                  d3 = e-e3
                  xx = d4*d4*d4
                  yy = d41*d42*d43
                  x = xx/(d41*yy)
                  dos(idos,1) = x
                  dosin(idos,1) = 0.25d0*(1.0d0-d4*x)
                  x = xx/(d42*yy)
                  dos(idos,2) = x
                  dosin(idos,2) = 0.25d0*(1.0d0-d4*x)
                  x = xx/(d43*yy)
                  dos(idos,3) = x
                  dosin(idos,3) = 0.25d0*(1.0d0-d4*x)
                  x = d3/d43+d2/d42+d1/d41
                  y = (d4*d4)/yy
                  dos(idos,4) = x*y
                  dosin(idos,4) = 0.25d0*(1.0d0-d4*y*(x+1.0d0))
               end do
               call gettod(wct_now)
               etime3 = etime3 + (wct_now-wct_start)*1.d-6

               call gettod(wct_start)
               ncounter2 = ncounter2+(ne2-ns2+1)
               td21 = 3.0d0*d21
               yy3=1.0d0/(d31*d42)+1.0d0/(d41*d32)
               do idos = ns2, ne2
                  e = es + deltae*real(idos,kind=DP)
!$$$               e = es+deltae*(idos-0.5)
                  d1 =e-e1
                  d4 =e4-e
                  d2 =e-e2
                  d3 =e3-e
                  d2d2 = d2*d2
                  d2d2d2 = d2d2*d2
                  d3d3 = d3*d3
                  x1 = (d3d3)/(d31*d31*d32)*(d2/d42+d1/d41)
                  x2 = (d4*d4)/(d41*d41*d42)*(d2/d32+d1/d31)
                  x3 = (d3*d4*d1)/(d31*d41)*yy3
                  dos(idos,1) = 0.5d0*(x1+x2+x3)
                  x1 = d2d2*(d32*d2+3.0d0*d3*(d32+d3))*factor_1_12
                  y1 = x1
                  x1 = x1/(d31*d31*d32*d42)
                  x2 = d2*(d2d2*(d31+td21)+3.0d0*d3*(d2*d3+d32*(td21+d1)))
                  x2 = x2*factor_1_12
                  y2 = x2
                  x2 = x2/(d31*d31*d32*d41)
                  x3 = d2d2*(d42*d2+3.0d0*d4*(d42+d4))*factor_1_12
                  y3 = x3
                  x3 = x3/(d41*d41*d42*d32)
                  x4 = d2*(d2d2*(d41+td21)+3.0d0*d4*(d2*d4+d42*(td21+d1)))
                  x4 = x4*factor_1_12
                  y4 = x4
                  x4 = x4/(d41*d41*d42*d31)
                  x5 = 0.5d0*d2*d3*d4*(d1+d21)
                  x5 = x5+d2d2*(2.0d0*d21*(d3+d42) &
                &   + (d1+d21)*(2.0d0*d3+d4+d42))*factor_1_12
                  x5 = x5*yy3/(d31*d41)
                  x6 = 0.25d0*d21*d21*(d42/d41+d32/d31+1.0d0)/(d41*d31)
                  dosin(idos,1) = 0.5d0*(x1+x2+x3+x4+x5)+x6

                  x1 = (d3d3)/(d32*d32*d31)*(d2/d42+d1/d41)
                  x2 = (d4*d4)/(d42*d42*d41)*(d2/d32+d1/d31)
                  x3 = (d3*d4*d2)/(d42*d32)*yy3
                  dos(idos,2) = 0.5d0*(x1+x2+x3)
                  x1 = y1/(d32*d32*d31*d42)
                  x2 = y2/(d32*d32*d31*d41)
                  x3 = y3/(d42*d42*d41*d32)
                  x4 = y4/(d42*d42*d41*d31)
                  x5 = d2d2*(d3*(d42+3.0d0*d4)+d32*(d42+d4))*factor_1_12
                  x5 =x5*yy3/(d42*d32)
                  x6 =0.25d0*d21/d31*d21/d41
                  dosin(idos,2) = 0.5d0*(x1+x2+x3+x4+x5)+x6

                  x1 = (d2d2)/(d32*d32*d42)*(d3/d31+d4/d41)
                  x2 = (d1*d1)/(d31*d31*d41)*(d3/d32+d4/d42)
                  x3 = (d1*d2*d3)/(d32*d31)*yy3
                  dos(idos,3) = 0.5d0*(x1+x2+x3)
                  x1 = d2d2d2*(3.0d0*d3+d32)*factor_1_12
                  y1 = x1
                  x1 = x1/(d32*d32*d42*d31)
                  x2 = d2d2d2*(3.0d0*d4+d42)*factor_1_12
                  y2 = x2
                  x2 = x2/(d32*d32*d42*d41)
                  x3 = d2*(d2*d31*(d2+td21)+3.0d0*d3*(d2d2+td21*d1)+ &
                & 3.0d0*d21*d21*d32)*factor_1_12
                  y3 = x3
                  x3 = x3/(d31*d31*d41*d32)
                  x4 = d2*(d2*d41*(d2+td21)+3.0d0*d4*(d2d2+td21*d1)+3.0d0* &
                  d21*d21*d42)*factor_1_12
                  y4 = x4
                  x4 = x4/(d31*d31*d41*d42)
                  x5 = (d2d2)*(d3*(d21+3.0d0*d1)+d32*(d21+d1))*factor_1_12
                  x5 = x5/(d32*d31)*yy3
                  x6 = 0.25d0*(d21/d31*d21/d31*d21/d41)
                  dosin(idos,3) = 0.5d0*(x1+x2+x3+x4+x5)+x6

                  x1 = (d2d2)/(d42*d42*d32)*(d3/d31+d4/d41)
                  x2 = (d1*d1)/(d41*d41*d31)*(d3/d32+d4/d42)
                  x3 = (d1*d2*d4)/(d41*d42)*yy3
                  dos(idos,4) = 0.5d0*(x1+x2+x3)
                  x1 = y1/(d42*d42*d32*d31)
                  x2 = y2/(d42*d42*d32*d41)
                  x3 = y3/(d41*d41*d31*d32)
                  x4 = y4/(d41*d41*d31*d42)
                  x5 =(d2d2)*(d4*(d21+3.0d0*d1)+d42*(d21+d1))*factor_1_12
                  x5 = x5/(d41*d42)*yy3
                  x6= 0.25d0*(d21/d31*d21/d41*d21/d41)
                  dosin(idos,4) = 0.5d0*(x1+x2+x3+x4+x5)+x6
               end do
               call gettod(wct_now)
               etime2 = etime2 + (wct_now-wct_start)*1.d-6
               call gettod(wct_start)
               if(idim == -3) then
                  esum = e1+e2+e3+e4
                  do idos = ns, ne
                     tdos = 0.025d0*(dos(idos,1)+dos(idos,2) &
                    & +dos(idos,3)+dos(idos,4))
                     dosin(idos,1) = dosin(idos,1)+tdos*(esum-4.d0*e1)
                     dosin(idos,2) = dosin(idos,2)+tdos*(esum-4.d0*e2)
                     dosin(idos,3) = dosin(idos,3)+tdos*(esum-4.d0*e3)
                     dosin(idos,4) = dosin(idos,4)+tdos*(esum-4.d0*e4)
                  end do
               end if
               do ip = 1,4
                  do ie = ns, ne
                     cdos(ie) = cdos(ie) &
                  & +dos(ie,ip)*dos_weight(ib,ieb(ip))
                     csumdos(ie) = csumdos(ie) &
                  & +dosin(ie,ip)*dos_weight(ib,ieb(ip))
                 end do
               end do
               w = (dos_weight(ib,ieb(1))+dos_weight(ib,ieb(2)) &
               &   +dos_weight(ib,ieb(3))+dos_weight(ib,ieb(4)))*0.25d0
               do ie = ne+1,nEwindows
                  csumdos(ie) = csumdos(ie) + 1.d0*w
               end do
               call gettod(wct_now)
               etime4 = etime4 + (wct_now-wct_start)*1.d-6
            end do
         end do
         call gettod(wct_start)
         cdos(0:newindows) = cdos(0:newindows)/real(ntet,kind=DP)
         csumdos(0:newindows) = csumdos(0:newindows)/real(ntet,kind=DP)

         if (npes  > 1) then
            itmp = (nep+1)*4
            allocate(dos_mpi(0:nep,4));dos_mpi=0.d0
            allocate(dosin_mpi(0:nep,4));dosin_mpi=0.d0
            allocate(cdos_mpi(0:newindows));cdos_mpi=0.d0
            allocate(csumdos_mpi(0:newindows));csumdos_mpi=0.d0

            call mpi_allreduce(dos,dos_mpi,itmp,mpi_double_precision, &
           & mpi_sum,MPI_CommGroup,ierr)
            dos = dos_mpi

            call mpi_allreduce(dosin,dosin_mpi,itmp,mpi_double_precision, &
           & mpi_sum,MPI_CommGroup,ierr)
            dosin = dosin_mpi

            call mpi_allreduce(cdos,cdos_mpi,(newindows+1), &
           & mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            cdos = cdos_mpi

            call mpi_allreduce(csumdos,csumdos_mpi,(newindows+1), &
           & mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
            csumdos = csumdos_mpi

            deallocate(dos_mpi) ; deallocate(dosin_mpi)
            deallocate(cdos_mpi)
            deallocate(csumdos_mpi)
         endif

         call gettod(wct_now)
         etime5 = etime5 + (wct_now-wct_start)*1.d-6

         if(ipri > 2) then
            write(nfout,'(" !dos iloop    = ",i12)') iloop
            write(nfout,'(" !dos ncounter1 = ",i12)') ncounter1
            write(nfout,'(" !dos ncounter2 = ",i12)') ncounter2
            write(nfout,'(" !dos ncounter3 = ",i12)') ncounter3
            write(nfout,'(" !dos etime0    = ",f16.8)') etime0
            write(nfout,'(" !dos etime1    = ",f16.8)') etime1
            write(nfout,'(" !dos etime2    = ",f16.8)') etime2
            write(nfout,'(" !dos etime3    = ",f16.8)') etime3
            write(nfout,'(" !dos etime4    = ",f16.8)') etime4
            write(nfout,'(" !dos etime5    = ",f16.8)') etime5
         end if
      end subroutine nstt5i
