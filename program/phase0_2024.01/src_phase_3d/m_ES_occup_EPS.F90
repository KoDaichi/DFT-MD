!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Module Name : m_ES_occup_EPS                                     #
!#                                                                       #
!#                                Written by T. Hamada 2007/5/8          #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
module m_ES_occup_EPS
!     (m_ESoc)
! $Id: m_ES_occup_EPS.F90 410 2014-10-30 05:20:45Z jkoga $
  use m_Electronic_Structure, only : efermi, occup_l, eko_l, totch, neordr, eko_ek &
       &                            , vbm, metalic_system, check_if_metalic_flag
  use m_Kpoints,              only : kv3, kv3_ek, qwgt, qwgt_ek &
       &                            ,np0,np2,ip20,iwt,ip2cub,nxyz_tetra
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,   only : ipri, iprioccup, width, nspin, neg, af, ekmode, width_tetra
!!$       &                            ,deltaE_dos_Gaussdistrib, variance_dos_Gaussdistrib
  use m_Const_Parameters,     only : DP, DELTA_FermiSearchRange, ON, OFF
  use m_Parallelization,      only : npes,map_ek,mype,map_e,map_k,myrank_e,myrank_k &
       &                            ,ierr,np_e,map_z,ista_e,ista_k,iend_k
!! UVSOR
  use m_ES_occup,             only : check_totch, check_if_metalic
!! UVSOR

! ==== KT_add ==== 2014/09/22
  use m_Control_Parameters,    only : noncol, ndim_spinor
! ================ 2014/09/22
  use mpi

  implicit none
!   1. m_ESoc_fermi_parabolic       <-(ChargeDensity_Construction) ->(3)
!      - get_occup_l_and_tot, - find_eko_minimum_and_maximum
!   2. m_ESoc_fermi_tetrahedron     <-(ChargeDensity_Construction) ->(3)
!   3. check_totch           <-(1),(2)
!
!   The subroutine "m_ESoc_fermi_parabolic" was originally coded by
!   Noriaki Hamada, 1984.05.17
!
!! UVSOR
  real(kind=DP), allocatable, dimension(:,:) :: occup_l_ek
!! UVSOR
!  include 'mpif.h'
contains

!! UVSOR
  subroutine m_ESoc_EPS_alloc_occup_l_ek
   allocate(occup_l_ek(neg,kv3_ek)) ; occup_l_ek=0.0d0
  end subroutine m_ESoc_EPS_alloc_occup_l_ek

  subroutine m_ESoc_EPS_dealloc_occup_l_ek
   deallocate(occup_l_ek)
  end subroutine m_ESoc_EPS_dealloc_occup_l_ek
!! UVSOR

  subroutine m_ESoc_EPS_fermi_parabolic_ek(nfout)
!
!   Coded from <m_ESoc_fermi_parabolic> 
!                      by T. Yamasaki (FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    integer :: id_sname = -1
!! UVSOR
    integer             :: ik, ie
    real(kind=DP),allocatable,dimension(:,:) :: occup_mpi, eko_tmp

    allocate(occup_mpi(neg,kv3_ek)); occup_mpi =0.0d0
!! UVSOR

    call tstatc0_begin('m_ESoc_EPS_fermi_parabolic_ek ', id_sname)

    call check_totch(nfout)

! ==== KT_mod ====== 2014/09/22
!    emin = minval(eko_ek)
!    emax = maxval(eko_ek)
!
    if ( noncol ) then
       allocate( eko_tmp( neg, kv3_ek/ndim_spinor ) )
       Do ik=1, kv3_ek /ndim_spinor
          eko_tmp(:,ik) = eko_ek(:,(ik-1)*ndim_spinor+1)
       End do
       emin = minval( eko_tmp );  emax = maxval( eko_tmp )
       deallocate( eko_tmp )
    else
       emin = minval(eko_ek);     emax = maxval(eko_ek)
    endif
! =================== 2014/09/22

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

! === KT_mod ===== 2014/09/22
!       call get_tot_eps()             ! -(contained here) ->tot
!          ~~~~~~~~~~~~~~~~~~~~
       if ( noncol ) then
          call get_tot_eps_noncl()
       else
          call get_tot_eps()
       endif
! =============== 2014/09/22

       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)
       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch
    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    call check_if_metalic(nfout,metalic_system)
    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band enrgy = ", f10.4)') vbm

   call tstatc0_end(id_sname)

! ==== KT_mod ==== 2014/09/22
!    do ik = 1, kv3_ek, af+1                             ! MPI
!       do ie = 1, neg                                   ! MPI
!          if(map_e(ie) == myrank_e) then                ! MPI
!             occup_l_ek(map_z(ie),ik) = occup_mpi(ie,ik)! MPI
!          end if                                        ! MPI
!       end do                                           ! MPI
!    end do                                              ! MPI
!   deallocate(occup_mpi)
!
   if ( noncol ) then
      do ik = 1, kv3_ek, ndim_spinor
         do ie = 1, neg
            if(map_e(ie) == myrank_e) then
               occup_l_ek(map_z(ie),ik)   = occup_mpi(ie,ik)
               occup_l_ek(map_z(ie),ik+1) = occup_mpi(ie,ik)
            end if
         end do
      end do
   else
      do ik = 1, kv3_ek, af+1                             ! MPI
         do ie = 1, neg                                   ! MPI
            if(map_e(ie) == myrank_e) then                ! MPI
               occup_l_ek(map_z(ie),ik) = occup_mpi(ie,ik)! MPI
            end if                                        ! MPI
         end do                                           ! MPI
      end do                                              ! MPI
   endif

   deallocate(occup_mpi)
! ================ 2014/09/22

  contains

    subroutine get_tot_eps
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw
      tot = 0.d0
      do k = 1, kv3_ek, af+1
         do i = 1, neg
            e = eko_ek(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight*wspin*kv3_ek*qwgt_ek(k)
!! UVSOR
            occup_mpi(i,k) = totw
!! UVSOR
            tot = tot + 2*totw
         end do
      end do
      tot = tot/kv3_ek * (af+1)
      if(iprioccup >= 2) write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') efermi,tot, jcount
    end subroutine get_tot_eps

! === KT_add === 2014/09/22
    subroutine get_tot_eps_noncl
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw, wktmp

      tot = 0.d0
      wktmp = kv3_ek / ndim_spinor

      do k = 1, kv3_ek, ndim_spinor
         do i = 1, neg
            e = eko_ek(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight *wspin *wktmp *qwgt_ek(k)
            occup_mpi(i,k) = totw
            tot = tot + totw
         end do
      end do

      tot = tot /wktmp
      if(iprioccup >= 2) then
         write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') &
              &          efermi,tot, jcount
      endif

    end subroutine get_tot_eps_noncl
! ========= 2014/09/22

  end subroutine m_ESoc_EPS_fermi_parabolic_ek

  subroutine m_ESoc_EPS_fermi_tetra_ek(nfout)
!  Coded from <m_ESoc_fermi_tetrahedron>
!          by Takahiro Yamasaki(FUJITSU LABORATORIES Ltd.), 28th Jun. 2003
!
    integer, intent(in) :: nfout
#ifndef NO_TETRAHEDRON
    real(kind=DP), parameter :: delta = 1.d-12
    integer, parameter       :: idim = 3
    integer        :: neig,nengy,ispin,ip2,ik,instts1,ikee,nxx,nyy,nzz, ib, jb
    real(kind=DP)  :: efermi2,eval,totind, et
    real(kind=DP), pointer, dimension(:,:,:) :: eig2
    real(kind=DP), pointer, dimension(:) :: eawk,cdwk,cswk,cdos,cind,valud
    integer             :: id_sname = -1
!! UVSOR
    integer        :: ip_mpi, ieig
    integer, allocatable, dimension(:,:) :: neordr_ek
    real(kind=DP), pointer, dimension(:,:,:) :: occup2
!! UVSOR

    call tstatc0_begin('m_ESoc_EPS_fermi_tetra_ek ', id_sname)

    if ( noncol ) stop "PPPPPPPPPP"

    allocate(eig2(np2,neg,nspin)); eig2 = 0.d0
    allocate(eawk(np0)); eawk = 0.d0
    allocate(cdwk(np0)); cdwk = 0.d0
    allocate(cswk(np0)); cswk = 0.d0
    allocate(cdos(np2*neg)); cdos = 0.d0
    allocate(cind(np2*neg)); cind = 0.d0
    allocate(valud(nspin)); valud = 0.d0
!! UVSOR
    allocate(neordr_ek(neg,kv3_ek)) ; neordr_ek = 0
    allocate(occup2(neg,np2,nspin)) ; occup2 = 0.0d0
!! UVSOR

    call check_totch(nfout)           ! totch is checked

    nxx = nxyz_tetra(1)
    nyy = nxyz_tetra(2)
    nzz = nxyz_tetra(3)
    neig=neg
    nengy=0
    do ispin=1,nspin
       do ip2=1,np2
          ik=nspin*(ip2-1)+ispin
!! UVSOR
          neordr_ek(1:neg,ik) = (/(ib,ib = 1, neg)/)
!! UVSOR
          do ib=1,neg
             eig2(ip2,ib,ispin)=eko_ek(ib,ik)
          enddo
          do ib = 1,neg-1
             do jb = ib+1, neg
                if(eig2(ip2,jb,ispin) < eig2(ip2,ib,ispin)-delta) then
                   et = eig2(ip2,ib,ispin)
                   eig2(ip2,ib,ispin) = eig2(ip2,jb,ispin)
                   eig2(ip2,jb,ispin) = et
!! UVSOR
                   neordr_ek(ib,ik) =jb
                   neordr_ek(jb,ik) =ib
!! UVSOR
                end if
             end do
          end do
       enddo
    enddo

    if(iprioccup>=2) write(nfout,*) ' === tetrahedron method', &
                &  ' for k-space integration ==='
    if(iprioccup>=2) then
       write(nfout,'(" --- eig2 ---")')
       do ip2 = 1, np2
! === DEBUG by tkato 2013/11/19 ================================================
!         write(nfout,'(" (ip2 = ",i8,")",99(/,8d12.4))') ip2,(eig2(ip2,ib,ispin),ib=1,neg)
          write(nfout,'(" (ip2 = ",i8,")",99(/,8d12.4))') ip2,(eig2(ip2,ib,nspin),ib=1,neg)
! ==============================================================================
       end do
    end if
    efermi2 = efermi
! === Optional argument needs interface!!! by T.Kato 2013/07/02 ================
!   call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
!            &  eig2,ip20,np0,totch, efermi2,eval,valud, &
!            &  iwt,ip2cub,iprioccup)
    call fermi1(nxx,nyy,nzz,np2,np2,neig,neg,nspin, &
             &  eig2,ip20,np0,totch, efermi2,eval,valud, &
             &  iwt,ip2cub,iprioccup,.false.)
! ==============================================================================
    if(iprioccup>=2) write(nfout,*) 'eval=',eval
    do ispin=1,nspin
!!$      write(nfout,*) ' ispin=',ispin
      instts1 = 0
      call nsdos3(nfout,idim,efermi2,nxx,nyy,nzz,np2,1,neig, &
               &  eig2(1,1,ispin),ip20,np0,eawk,instts1,np2,iwt,ip2cub,iprioccup)
      call nstt3i(idim,nengy,efermi2,nxx,nyy,nzz, &
               &  np2,np2,neig,eig2(1,1,ispin), &
               &  ip20,np0,eawk,cdwk,cswk,np2,neig,cdos,cind,width_tetra )
      efermi=efermi2

      totind=0.d0
      do ip2=1,np2
         do ib=1,neig
            ikee=ip2+np2*(ib-1)
!! UVSOR
            ik = nspin * (ip2-1) + ispin
            ip_mpi = neordr_ek(ib,ik)
            occup2(ip_mpi,ip2,ispin) = cind(ikee)*dble(np2)
!! UVSOR
            totind=totind+cind(ikee)
         enddo
      enddo
      if(iprioccup>=2) write(nfout,*) 'for spin=',ispin, &
                  &  ' ** TOTAL CHARGE after fermi1 = ',totind
    enddo

!! UVSOR
    do ieig = 1, neg
       if(map_e(ieig) /= myrank_e) cycle                          ! MPI
       do ispin=1,nspin
          do ip2=1,np2
             ik=nspin*(ip2-1)+ispin
!             if(map_k(ik) == myrank_k) then                      ! MPI
                occup_l_ek(map_z(ieig),ik)=occup2(ieig,ip2,ispin) ! MPI
!             end if
          enddo
       enddo
    enddo
    deallocate(neordr_ek)
    deallocate(occup2)
!! UVSOR

    deallocate(valud); deallocate(cind); deallocate(cdos); deallocate(cswk)
    deallocate(cdwk);  deallocate(eawk)
    deallocate(eig2)
    
    call check_if_metalic(nfout,metalic_system)
    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band enrgy = ", f10.4)') vbm

    call tstatc0_end(id_sname)
#endif
  end subroutine m_ESoc_EPS_fermi_tetra_ek

! =========================== KT_add ======================= 13.0E
  subroutine m_ESoc_EPS_fermi_dirac_ek(nfout)
    integer, intent(in) :: nfout

    integer             :: jcount
    real(kind=DP)       :: emin, emax, e1, e2, tot
    real(kind=DP), parameter :: DELTA_TOTCH = 1.d-10
    integer,       parameter :: MAXITR = 100000
    integer :: id_sname = -1

    integer             :: ik, ie
    real(kind=DP),allocatable,dimension(:,:) :: occup_mpi, eko_tmp

    allocate(occup_mpi(neg,kv3_ek)); occup_mpi =0.0d0

    call tstatc0_begin('m_ESoc_EPS_fermi_dirac_ek ', id_sname)

    call check_totch(nfout)

! ==== KT_mod ====== 2014/09/22
!    emin = minval(eko_ek)
!    emax = maxval(eko_ek)
!
    if ( noncol ) then
       allocate( eko_tmp( neg, kv3_ek/ndim_spinor ) )
       Do ik=1, kv3_ek /ndim_spinor
          eko_tmp(:,ik) = eko_ek(:,(ik-1)*ndim_spinor+1)
       End do
       emin = minval( eko_tmp );  emax = maxval( eko_tmp )
       deallocate( eko_tmp )
    else
       emin = minval(eko_ek);     emax = maxval(eko_ek)
    endif
! =================== 2014/09/22

    efermi = emax
    e1 = emin - dabs(DELTA_FermiSearchRange)
    e2 = emax + dabs(DELTA_FermiSearchRange)

    jcount = 1
    occupied_ch_equals_totch : do

! === KT_mod ===== 2014/09/22
!       call get_tot_eps()             ! -(contained here) ->tot
!          ~~~~~~~~~~~~~~~~~~~~
       if ( noncol ) then
          call get_tot_eps_noncl()
       else
          call get_tot_eps()
       endif
! =============== 2014/09/22

       if(jcount == 1 .and. tot < totch) &
            call wd_fermi_error1(nfout,emin,emax,tot,totch) ! -(b_Fermi)

       if(dabs(tot - totch) < DELTA_TOTCH) exit occupied_ch_equals_totch
       if( tot < totch) then
          e1 = efermi; efermi = efermi + (e2-efermi)/2
       else
          e2 = efermi; efermi = efermi + (e1-efermi)/2
       end if
       jcount = jcount + 1
       if(jcount > MAXITR) then
          call wd_fermi_error2 &   !-(b_Fermi)
               & (nfout,e1,e2,efermi,emin,emax,tot,totch,neg,MAXITR)
       end if
    end do occupied_ch_equals_totch

    if(iprioccup >= 2) &
         &write(nfout,'( " efermi = ", f10.4, " tot = ", f10.4)') efermi, tot

    call check_if_metalic(nfout,metalic_system)
    if(iprioccup >= 2 .and. .not.metalic_system) &
         &write(nfout,'( " The highest occupied band enrgy = ", f10.4)') vbm

   call tstatc0_end(id_sname)

! ==== KT_mod ==== 2014/09/22
!    do ik = 1, kv3_ek, af+1                             ! MPI
!       do ie = 1, neg                                   ! MPI
!          if(map_e(ie) == myrank_e) then                ! MPI
!             occup_l_ek(map_z(ie),ik) = occup_mpi(ie,ik)! MPI
!          end if                                        ! MPI
!       end do                                           ! MPI
!    end do                                              ! MPI
!   deallocate(occup_mpi)
!
   if ( noncol ) then
      do ik = 1, kv3_ek, ndim_spinor
         do ie = 1, neg
            if(map_e(ie) == myrank_e) then
               occup_l_ek(map_z(ie),ik)   = occup_mpi(ie,ik)
               occup_l_ek(map_z(ie),ik+1) = occup_mpi(ie,ik)
            end if
         end do
      end do
   else
      do ik = 1, kv3_ek, af+1                             ! MPI
         do ie = 1, neg                                   ! MPI
            if(map_e(ie) == myrank_e) then                ! MPI
               occup_l_ek(map_z(ie),ik) = occup_mpi(ie,ik)! MPI
            end if                                        ! MPI
         end do                                           ! MPI
      end do                                              ! MPI
   endif

   deallocate(occup_mpi)
! ================ 2014/09/22

  contains

    subroutine get_tot_eps
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw

      tot = 0.d0
      do k = 1, kv3_ek, af+1
         do i = 1, neg
            e = eko_ek(i,k)

            call width_fermi_dirac(e,efermi,width,dos,weight)  ! -(b_Fermi)

            totw = weight*wspin*kv3_ek*qwgt_ek(k)
            occup_mpi(i,k) = totw
            tot = tot + 2*totw
         end do
      end do

      tot = tot/kv3_ek * (af+1)
      if (iprioccup >= 2) then
         write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') &
              &                         efermi,tot, jcount
      endif

    end subroutine get_tot_eps

! === KT_add === 2014/09/22
    subroutine get_tot_eps_noncl
      integer       :: k, i
      real(kind=DP) :: wspin = 1.d0, e, dos, weight, totw, wktmp

      tot = 0.d0
      wktmp = kv3_ek / ndim_spinor

      do k = 1, kv3_ek, ndim_spinor
         do i = 1, neg
            e = eko_ek(i,k)
            call width2(e,efermi,width,dos,weight)  ! -(b_Fermi)
            totw = weight *wspin *wktmp *qwgt_ek(k)
            occup_mpi(i,k) = totw
            tot = tot + totw
         end do
      end do

      tot = tot /wktmp
      if(iprioccup >= 2) then
         write(nfout,'(" efermi, tot = ",d20.12,d20.12," jcount = ",i7)') &
              &          efermi,tot, jcount
      endif

    end subroutine get_tot_eps_noncl
! ========= 2014/09/22

  end subroutine m_ESoc_EPS_fermi_dirac_ek
! ============================================================ 13.0E
end module m_ES_occup_EPS
