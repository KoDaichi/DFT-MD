!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 594 $)
!
!  MODULE: m_ES_Intgr_VlhxcQlm
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki   May/17/2005
!      Further modification by T. Yamasaki   Aug/31/2007
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
!

#ifdef VPP
#define _VECTOR_TUNING_
#endif
#ifdef SX
#define _VECTOR_TUNING_
#endif
#ifdef HIUX
#define _VECTOR_TUNING_
#endif

module m_ES_Intgr_VlhxcQlm
! $Id: m_ES_Intgr_VlhxcQlm.F90 594 2019-07-01 07:15:23Z jkoga $
  use m_Electronic_Structure,only : vlhxcQ, vlhxc_l, dhub
  use m_PlaneWaveBasisSet,   only : kgp,ngabc,m_pwBS_sphrp2,ylm_l
  use m_PseudoPotential,     only : ilmt,ltp,taup,isph,il2p,iqitg,qitg_l,dl2p &
       &                           ,modnrm,m_PP_find_maximum_l, nqitg, nlmt &
       &                           ,m_PP_include_vanderbilt_pot, ivanl &
       &                           ,m_PP_set_index_arrays1 &
       &                           ,m_PP_set_index_arrays2
  use m_Ionic_System,        only : ntyp, natm, ityp, pos, ihubbard
  use m_Crystal_Structure,   only : af, rltv, univol
  use m_Timing,              only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters,  only : nspin,iprivlhxcq,kimg,nel_Ylm, m_CtrlP_cachesize, sw_communicator_for_chg
  use m_Const_Parameters,    only : DP, PAI2, zi, SKIP, OFF
  use m_Parallelization,     only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr, np_e, mpi_chg_world
!$$#ifdef PARA3D
  use m_PlaneWaveBasisSet,   only : m_pwBS_sphrp2_3D
  use m_PlaneWaveBasisSet,   only : ngabc_kngp_l, ngabc_kngp_B_l
! ==== DEBUG by tkato 2012/04/03 ===============================================
  use m_Parallelization,     only : mpi_k_world, myrank_k
! ==============================================================================
!$$#endif

! ========================== added by K. Tagami ================== 11.0
  use m_Control_Parameters,  only : ndim_magmom
! ================================================================ 11.0
  use mpi

  implicit none

!  include 'mpif.h'
!  52. m_ESiVQ_integrate_VlhxcQlm
#ifdef _VECTOR_TUNING_
!      - symmetrize_vlhxcq, - integrate_qitg_x_ylm, - integration_of_VlhxcQlm_core3,
#else
!      - symmetrize_vlhxcq, - integration_of_VlhxcQlm_core2,
#endif
!      - integration_of_VlhxcQlm_core
contains


  subroutine m_ESiVQ_add_dhub_to_vlhxcQ(nfout)
    integer, intent(in) :: nfout
    integer :: is,ia
! ========================================== modified by K. Tagami ======== 11.0
!    do is=1,nspin
    do is=1, ndim_magmom
! ========================================================================= 11.0
       do ia=1,natm
          if(ihubbard(ia) > 0) then
             vlhxcQ(:,:,ia,is) = vlhxcQ(:,:,ia,is)+dhub(:,:,ia,is)
          end if
       end do
    end do
  end subroutine m_ESiVQ_add_dhub_to_vlhxcQ

!$$#ifdef PARA3D
!=======================================================================
!=======================================================================
!=======================================================================
!! #define _NOPARA_KNGP_B_
  subroutine m_ESiVQ_integrate_VlhxcQlm_3D(nfout)
   use m_Parallelization,     only : ista_kngp_B , iend_kngp_B    &
  &                                , ista_atm_B, iend_atm_B, mem_atm_B &
  &                                , mpi_ke_world

    ! ***
    ! (Rev) T. Yamaskai, 31, Aug, 2007
    !     1. 'call set_index_arrays1' that included a bug is replaced
    !       by 'call m_PP_set_index_arrays1', whose bug is fixed.
    !     2. 'call set_index_arrays2' is also replaced by 'call
    !       m_PP_set_index_arrays2' that can be referred from other modules.
    !     3. contained subroutines, set_index_arrays1 and set_index_arrays2 were
    !       deleted.
    ! ***
    integer, intent(in)   :: nfout

    real(kind=DP), allocatable, dimension(:)           :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext
    real(kind=DP), allocatable, target, dimension(:)   :: zfcos, zfsin

    integer :: ncache, ibsize, ibl1, ibl2
    real(kind=DP), allocatable, dimension(:)           :: veq ! d(ista_kngp:iend_kngp)
    real(kind=DP)                                      :: vq_ia
    real(kind=DP), allocatable, dimension(:,:)         :: qitg_red, ylm_red, vlhxc_red
    integer,       allocatable, dimension(:)           :: il3

    integer :: is,it,lmt1,lmt2,n,ilm3,l3,ia

    integer :: m, maxm, ip, np, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)

    real(kind=DP), allocatable, dimension(:,:,:,:) :: vlhxcQ_mpi ! d(nlmt,nlmt,natm,nspin)
    integer :: iab
    integer :: ista,iend
    integer :: id_sname = -1

    integer :: i, ilm
    if(modnrm == SKIP) return
    if(sw_communicator_for_chg == OFF)then
      ista = ista_kngp_B
      iend = iend_kngp_B
    else
      ista = ista_kngp
      iend = iend_kngp
    endif
!--- dummy timer --->
#ifdef __TIMER_SUB__
  call timer_sta(1061)
  call timer_end(1061)
  call timer_sta(1062)
  call timer_end(1062)
#endif

#ifdef __TIMER_DO__
  call timer_sta(1085)
  call timer_end(1085)
#endif
!<--- dummy timer ---

    call tstatc0_begin('m_ESiVQ_integrate_VlhxcQlm_3D ',id_sname)

#ifdef __TIMER_SUB__
  call timer_sta(1057)
#endif

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)

    allocate(nqitg_sp(ntyp)); allocate(nqitg_sp0(ntyp))
    allocate(iq2l3(nqitg))
    allocate(nc(mcritical,nqitg));nc=0
    call m_PP_set_index_arrays1(nfout,ntyp,nqitg,mcritical,n**2,il3 &
         & ,maxm,mc,nqitg_sp,nqitg_sp0,iq2l3,nc)
    allocate(nc2lmt1(mc,maxm,nqitg))
    allocate(nc2lmt2(mc,maxm,nqitg))
    allocate(nc2n(mc,maxm,nqitg))
    call m_PP_set_index_arrays2(nfout,mc,maxm,nqitg,mcritical,n**2,il3,iq2l3 &
         & ,nc2lmt1,nc2lmt2,nc2n,nc) ! -> nc2lmt1, nc2lmt2, nc2n, nc

    if(n**2 > nel_Ylm) then
#ifndef _NOPARA_KNGP_B_
       allocate(ylm_ext(ista:iend,nel_Ylm+1:n**2));  ylm_ext = 0.d0
       allocate(ylm_t(ista:iend)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2_3D(ilm3,rltv,ista,iend,ylm_t)  !(ilm3,rltv,ngabc,gr_l)->(ylm)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
#else
       allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2));  ylm_ext = 0.d0
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2_3D(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)  !(ilm3,rltv,ngabc,gr_l)->(ylm)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
#endif
       deallocate(ylm_t)
    end if

    vlhxcQ = 0.d0

#ifndef _NOPARA_KNGP_B_
    do is = 1, nspin, af+1

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          ibsize = iend-ista+1
       else
          ibsize=ncache/(8*(n**2))
       end if
       if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ibsize = ", i8)') ibsize
!f       allocate(zfcos(ibsize)); zfcos = 0.d0
!f       allocate(zfsin(ibsize)); zfsin = 0.d0
!f       allocate(veq(ibsize))
       allocate(qitg_red(ibsize,nqitg))
       allocate(vlhxc_red(ibsize,kimg))
       allocate(ylm_red(ibsize,n**2))
!$omp parallel default(none)                                                    &
!$omp   shared(ityp,natm,ivanl,nqitg_sp0,nqitg_sp,iq2l3,nc,nc2lmt1,             &
!$omp          nc2lmt2,nc2n,vlhxcQ,is,univol,dl2p,pos,ista_kngp_B,iend_kngp_B,  &
!$omp          ngabc_kngp_B_l,qitg_red,vlhxc_red,ylm_red,kimg,ibsize,mype,      &
!$omp          qitg_l,vlhxc_l,ylm_l,ylm_ext,nqitg,nel_Ylm,n,ista,iend,          &
!$omp          sw_communicator_for_chg)                                         &
!$omp   private(ibl1,ibl2,ia,it,iq,l3,m,ilm3,vq_ia,ip,lmt1,lmt2,np,zfcos,zfsin,veq,i,ilm)
#ifdef __TIMER_OMP1__
  call timer_sta(1088)
#endif
       allocate(zfcos(ibsize)); zfcos = 0.d0
       allocate(zfsin(ibsize)); zfsin = 0.d0
       allocate(veq(ibsize))
#ifdef __TIMER_OMP1__
  call timer_end(1088)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1076)
#endif
!#ifdef USE_PROF
!  call timer_sta(1076)
!#endif
!-!OCL PREFETCH_INDIRECT
       do ibl1 = ista, iend, ibsize
          ibl2 = ibl1+ibsize-1
          if(ibl2 .gt. iend) ibl2 = iend
!oo       if(iprivlhxcq >= 2) then
!oo          write(nfout,'(" !mESiVQ:    ibl1,ibl2,ibsize = ", 3i8)') ibl1, ibl2, ibsize
!oo       end if
!f          call substitute_qitgred(ibl1,ibl2)  ! qitg_l -> qitg_red
!---
!$omp do
          do iq = 1, nqitg
             do i = 1, ibl2-ibl1+1
                qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
             end do
          end do
!$omp end do
!---
!f          call substitute_vlhxcred(ibl1,ibl2) ! vlhxc_l -> vlhxc_red
!---
          if(kimg==1) then
!$omp do
             do i = 1, ibl2-ibl1+1
                vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
             end do
!$omp end do
          else if(kimg==2) then
!$omp do
             do i = 1, ibl2-ibl1+1
                vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
                vlhxc_red(i,2) = vlhxc_l(i+ibl1-1,2,is)
             end do
!$omp end do
          end if
!---
!f          call substitute_ylmred(ibl1,ibl2)   ! ylm_l, ylm_ext -> ylm_red
!$omp do
      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do
!$omp end do
      if(n**2 > nel_Ylm) then
         do ilm = nel_ylm+1, n**2
!$omp do
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_ext(i+ibl1-1,ilm)
            end do
!$omp end do
         end do
      end if

#ifdef __TIMER_OMP1__
  call timer_sta(1087)
#endif
#ifdef __TIMER_OMP1__
  call timer_sta(1093)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1084)
#endif
!$omp do schedule(runtime)
          do ia = 1, natm
             it = ityp(ia)
             if(sum(ivanl(:,it)) == 0) cycle
!             call calc_phase_div(ia)  !   -> zfcos, zfsin
!f             call calc_phase_div_B_3D(ia)  !   -> zfcos, zfsin
#ifdef __TIMER_OMP3__
  call timer_sta(1089)
#endif
             if(sw_communicator_for_chg == OFF)then
             call calc_phase_div_B_3D(ia,zfcos,zfsin,ibsize,ibl1,ibl2)  !   -> zfcos, zfsin
             else
             call calc_phase_div_3D(ia)  !   -> zfcos, zfsin
             endif
#ifdef __TIMER_OMP3__
  call timer_end(1089)
#endif
!oo          if(iprivlhxcq >= 2) then
!oo             write(nfout,'(" !mESiVQ:    ia,  iq,  l3,   m,ilm3,  ip,lmt1,lmt2,  np")')
!oo          end if
#ifdef __TIMER_OMP2__
  call timer_sta(1095)
#endif
             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
!f                call dp_Vlhxcq_exp_Q_div(iq) ! vlhxc, exp, qitg_l  -> veq
#ifdef __TIMER_OMP3__
  call timer_sta(1090)
#endif
                call dp_Vlhxcq_exp_Q_div(iq,zfcos,zfsin,veq,ibsize,l3,ibl1,ibl2) ! vlhxc, exp, qitg_l  -> veq
#ifdef __TIMER_OMP3__
  call timer_end(1090)
#endif
#ifdef __TIMER_OMP2__
  call timer_sta(1096)
#endif
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
!f                   call veQ_dot_ylm_div(vq_ia) ! veq, ylm -> vq_ia
#ifdef __TIMER_OMP3__
  call timer_sta(1091)
#endif
                   call veQ_dot_ylm_div(vq_ia,veq,ibsize,ilm3,l3,ibl1,ibl2) ! veq, ylm -> vq_ia
#ifdef __TIMER_OMP3__
  call timer_end(1091)
#endif
!f--#ifdef __TIMER_DO__
!f--  call timer_sta(1076)
!f--#endif
#ifdef __TIMER_OMP2__
  call timer_sta(1092)
#endif
                   do ip = 1, nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq)
                      lmt2 = nc2lmt2(ip,m,iq)
                      np = nc2n(ip,m,iq)
!oo                   if(iprivlhxcq >= 2) then
!oo                      write(nfout,'(" !mESiVQ: ",9i5)') ia, iq, l3, m, ilm3, ip, lmt1,lmt2,np
!oo                   end if
                      vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is) &
                     &                          + univol*vq_ia*dl2p(lmt1,lmt2,np,it)
                   end do
#ifdef __TIMER_OMP2__
  call timer_end(1092)
#endif
!f--#ifdef __TIMER_DO__
!f--  call timer_end(1076)
!f--#endif
                end do
#ifdef __TIMER_OMP2__
  call timer_end(1096)
#endif
             end do
#ifdef __TIMER_OMP2__
  call timer_end(1095)
#endif
          end do
!$omp end do
#ifdef __TIMER_DO__
  call timer_end(1084)
#endif
#ifdef __TIMER_OMP1__
  call timer_end(1093)
#endif
#ifdef __TIMER_OMP1__
  call timer_end(1087)
#endif
       end do
!#ifdef USE_PROF
!  call timer_end(1076)
!#endif
#ifdef __TIMER_DO__
  call timer_end(1076)
#endif
#ifdef __TIMER_OMP1__
  call timer_sta(1094)
#endif
       deallocate(veq)
       deallocate(zfsin);deallocate(zfcos)
#ifdef __TIMER_OMP1__
  call timer_end(1094)
#endif
!$omp end parallel
       deallocate(ylm_red)
       deallocate(vlhxc_red)
       deallocate(qitg_red)
!f       deallocate(veq)

!f       deallocate(zfsin);deallocate(zfcos)
    end do

    if(npes > 1) then
       allocate(vlhxcQ_mpi(nlmt,nlmt,natm,nspin)); vlhxcQ_mpi = 0.d0
! ==== DEBUG by tkato 2012/04/03 ===============================================
#if 0
#ifdef __TIMER_COMM__
  call timer_barrier(MPI_CommGroup)
  call timer_sta(1077)
#endif
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, MPI_CommGroup, ierr)
#else
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_k_world(myrank_k))
  call timer_sta(1077)
#endif
       if(sw_communicator_for_chg == OFF)then
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, mpi_k_world(myrank_k), ierr)
       else
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, mpi_chg_world, ierr)
       endif
#endif
! ==============================================================================
#ifdef __TIMER_COMM__
  call timer_end(1077)
#endif
       vlhxcQ = vlhxcQ_mpi
       deallocate(vlhxcQ_mpi)
    end if

#else

    do is = 1, nspin, af+1

       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          ibsize = iend_kngp-ista_kngp+1
       else
          ibsize=ncache/(8*(n**2))
       end if
       if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ibsize = ", i8)') ibsize
       allocate(zfcos(ibsize)); zfcos = 0.d0
       allocate(zfsin(ibsize)); zfsin = 0.d0
       allocate(veq(ibsize))
       allocate(qitg_red(ibsize,nqitg))
       allocate(vlhxc_red(ibsize,kimg))
       allocate(ylm_red(ibsize,n**2))
       do ibl1 = ista_kngp, iend_kngp, ibsize
          ibl2 = ibl1+ibsize-1
          if(ibl2 .gt. iend_kngp) ibl2 = iend_kngp
!oo       if(iprivlhxcq >= 2) then
!oo          write(nfout,'(" !mESiVQ:    ibl1,ibl2,ibsize = ", 3i8)') ibl1, ibl2, ibsize
!oo       end if
          call substitute_qitgred(ibl1,ibl2)  ! qitg_l -> qitg_red
          call substitute_vlhxcred(ibl1,ibl2) ! vlhxc_l -> vlhxc_red
          call substitute_ylmred(ibl1,ibl2)   ! ylm_l, ylm_ext -> ylm_red

!!        do ia = ista_atm_B, iend_atm_B
          do iab = ista_atm_B, iend_atm_B
             ia = mem_atm_B(iab)
             it = ityp(ia)
             if(sum(ivanl(:,it)) == 0) cycle
!             call calc_phase_div(ia)  !   -> zfcos, zfsin
             call calc_phase_div_3D(ia)  !   -> zfcos, zfsin
!oo          if(iprivlhxcq >= 2) then
!oo             write(nfout,'(" !mESiVQ:    ia,  iq,  l3,   m,ilm3,  ip,lmt1,lmt2,  np")')
!oo          end if
             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
!f                call dp_Vlhxcq_exp_Q_div(iq) ! vlhxc, exp, qitg_l  -> veq
                call dp_Vlhxcq_exp_Q_div(iq,zfcos,zfsin,veq,ibsize,l3,ibl1,ibl2) ! vlhxc, exp, qitg_l  -> veq
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
!f                   call veQ_dot_ylm_div(vq_ia) ! veq, ylm -> vq_ia
                   call veQ_dot_ylm_div(vq_ia,veq,ibsize,ilm3,l3,ibl1,ibl2) ! veq, ylm -> vq_ia
#ifdef __TIMER_DO__
  call timer_sta(1078)
#endif
                   do ip = 1, nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq)
                      lmt2 = nc2lmt2(ip,m,iq)
                      np = nc2n(ip,m,iq)
!oo                   if(iprivlhxcq >= 2) then
!oo                      write(nfout,'(" !mESiVQ: ",9i5)') ia, iq, l3, m, ilm3, ip, lmt1,lmt2,np
!oo                   end if
                      vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is) &
                     &                          + univol*vq_ia*dl2p(lmt1,lmt2,np,it)
                   end do
#ifdef __TIMER_DO__
  call timer_end(1078)
#endif
                end do
             end do
          end do
       end do
       deallocate(ylm_red)
       deallocate(vlhxc_red)
       deallocate(qitg_red)
       deallocate(veq)

       deallocate(zfsin);deallocate(zfcos)
    end do

    if(npes > 1) then
       allocate(vlhxcQ_mpi(nlmt,nlmt,natm,nspin)); vlhxcQ_mpi = 0.d0
! ==== DEBUG by tkato 2012/04/03 ===============================================
#if 0
#ifdef __TIMER_COMM__
  call timer_barrier(MPI_CommGroup)
  call timer_sta(1079)
#endif
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, MPI_CommGroup, ierr)
#else
#ifdef __TIMER_COMM__
  call timer_barrier(mpi_k_world(myrank_k))
  call timer_sta(1079)
#endif
       if(sw_communicator_for_chg == OFF)then
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, mpi_k_world(myrank_k), ierr)
       else
       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
      &                   mpi_double_precision, mpi_sum, mpi_chg_world, ierr)
       endif
#endif
! ==============================================================================
#ifdef __TIMER_COMM__
  call timer_end(1079)
#endif
       vlhxcQ = vlhxcQ_mpi
       deallocate(vlhxcQ_mpi)
    end if
#endif

!   if(nrank_g > 1) then
!      allocate(vlhxcQ_mpi(nlmt,nlmt,natm,nspin)); vlhxcQ_mpi = 0.d0
!      call mpi_allreduce(vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*nspin,  &
!     &                   mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
!      call mpi_allgather(vlhxcQ_l_3D, nlmt*nlmt*nspin*np_atm_B, MPI_DOUBLE_PRECISION &
!     &                   vlhxcQ_mpi,  nlmt*nlmt*nspin*mp_atm_B, MPI_DOUBLE_PRECISION &
!     &                   mpi_ke_world, ierr )
!      vlhxcQ = vlhxcQ_mpi
!      deallocate(vlhxcQ_mpi)
!   end if

    call symmetrize_vlhxcq
    if(allocated(ylm_ext)) deallocate(ylm_ext)
    deallocate(il3)

    deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)

    if(iprivlhxcq >= 2) then
       write(nfout,*) ' <<< m_ESiVQ_integrate_VlhxcQlm >>> vlhxcQ start'
       do is = 1, nspin, af+1
          do it = 1, natm
             write(nfout,*) ' (vlhxcQ) -- #atom = ', it
             do lmt1 = 1, ilmt(ityp(it))
                write(nfout,*) ' --- lmt1 = ', lmt1
                write(nfout,'(99(5d16.8,/))')&
                     & (vlhxcq(lmt1,lmt2,it,is),lmt2=1,ilmt(ityp(it)))
             end do
          end do
       end do
       write(nfout,*) ' <<< m_ESiVQ_integrate_VlhxcQlm >>> vlhxcQ end'
    endif
#ifdef __TIMER_SUB__
  call timer_end(1057)
#endif
    call tstatc0_end(id_sname)
  contains

    subroutine symmetrize_vlhxcq
      integer :: is, ia, it, lmt1, lmt2
#ifdef __TIMER_SUB__
  call timer_sta(1063)
#endif
#ifdef __TIMER_DO__
  call timer_sta(1086)
#endif
      do is = 1, nspin, af+1
         do ia = 1, natm
            it = ityp(ia)
            do lmt1 = 1, ilmt(it)
               do lmt2 = lmt1+1, ilmt(it)
                  vlhxcQ(lmt2,lmt1,ia,is) = vlhxcQ(lmt1,lmt2,ia,is)
               end do
            end do
         end do
      end do
#ifdef __TIMER_DO__
  call timer_end(1086)
#endif
#ifdef __TIMER_SUB__
  call timer_end(1063)
#endif
    end subroutine symmetrize_vlhxcq

!f    subroutine substitute_qitgred()
    subroutine substitute_qitgred(ib1,ib2)
      integer, intent(in) :: ib1, ib2
      integer :: iq, i
      do iq = 1, nqitg
!f         do i = 1, ibl2-ibl1+1
!f            qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
         do i = 1, ib2-ib1+1
            qitg_red(i, iq) = qitg_l(i+ib1-1,iq)
         end do
      end do
    end subroutine substitute_qitgred

!f    subroutine substitute_vlhxcred
    subroutine substitute_vlhxcred(ib1,ib2)
      integer, intent(in) :: ib1, ib2
      integer :: i
#ifdef __TIMER_SUB__
  call timer_sta(1060)
#endif

#ifdef __TIMER_DO__
  call timer_sta(1083)
#endif
      if(kimg==1) then
!f         do i = 1, ibl2-ibl1+1
!f            vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
         do i = 1, ib2-ib1+1
            vlhxc_red(i,1) = vlhxc_l(i+ib1-1,1,is)
         end do
      else if(kimg==2) then
!f         do i = 1, ibl2-ibl1+1
!f            vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
!f            vlhxc_red(i,2) = vlhxc_l(i+ibl1-1,2,is)
         do i = 1, ib2-ib1+1
            vlhxc_red(i,1) = vlhxc_l(i+ib1-1,1,is)
            vlhxc_red(i,2) = vlhxc_l(i+ib1-1,2,is)
         end do
      end if
#ifdef __TIMER_DO__
  call timer_end(1083)
#endif
#ifdef __TIMER_SUB__
  call timer_end(1060)
#endif
    end subroutine substitute_vlhxcred

!f    subroutine substitute_ylmred
    subroutine substitute_ylmred(ib1,ib2)
      integer, intent(in) :: ib1, ib2
      integer :: ilm, i
      do ilm = 1, nel_Ylm
!f         do i = 1, ibl2-ibl1+1
!f            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         do i = 1, ib2-ib1+1
            ylm_red(i,ilm) = ylm_l(i+ib1-1,ilm)
         end do
      end do
      if(n**2 > nel_Ylm) then
         do ilm = nel_ylm+1, n**2
!f            do i = 1, ibl2-ibl1+1
!f               ylm_red(i,ilm) = ylm_ext(i+ibl1-1,ilm)
            do i = 1, ib2-ib1+1
               ylm_red(i,ilm) = ylm_ext(i+ib1-1,ilm)
            end do
         end do
      end if
    end subroutine substitute_ylmred

    subroutine calc_phase_div(ia)
      integer, intent(in) :: ia
      integer :: i, ig
      real(kind=DP) :: fx, fy, fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do i = 1, ibl2-ibl1+1
         ig = i + ibl1 - 1
         ph = ngabc(ig,1)*fx + ngabc(ig,2)*fy + ngabc(ig,3)*fz
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
    end subroutine calc_phase_div

    subroutine calc_phase_div_3D(ia)
      integer, intent(in) :: ia
      integer :: i, ig
      real(kind=DP) :: fx, fy, fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
      do i = 1, ibl2-ibl1+1
         ig = i + ibl1 - 1
         ph = ngabc_kngp_l(ig,1)*fx + ngabc_kngp_l(ig,2)*fy + ngabc_kngp_l(ig,3)*fz
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
    end subroutine calc_phase_div_3D

!f    subroutine calc_phase_div_B_3D(ia)
    subroutine calc_phase_div_B_3D(ia,zfcos0,zfsin0,ib,ib1,ib2)
      integer, intent(in) :: ia, ib, ib1, ib2
      real(kind=DP), intent(out) :: zfcos0(ib), zfsin0(ib)
      integer :: i, ig
      real(kind=DP) :: fx, fy, fz, ph
      fx = pos(ia,1)*PAI2
      fy = pos(ia,2)*PAI2
      fz = pos(ia,3)*PAI2
!f      do i = 1, ibl2-ibl1+1
!f         ig = i + ibl1 - 1
      do i = 1, ib2-ib1+1
         ig = i + ib1 - 1
         ph = ngabc_kngp_B_l(ig,1)*fx + ngabc_kngp_B_l(ig,2)*fy + ngabc_kngp_B_l(ig,3)*fz
         zfcos0(i) = dcos(ph)
         zfsin0(i) = dsin(ph)
      end do
    end subroutine calc_phase_div_B_3D

!f    subroutine dp_Vlhxcq_exp_Q_div(iq)
    subroutine dp_Vlhxcq_exp_Q_div(iq,zfcos0,zfsin0,veq0,ib,l,ib1,ib2)
      integer, intent(in) :: iq, ib, l, ib1, ib2
      real(kind=DP), intent(out) :: veq0(ib)
      real(kind=DP), target, intent(in) :: zfcos0(ib), zfsin0(ib)
      integer       :: i
      real(kind=DP) :: f
      real(kind=DP), pointer, dimension(:) :: zfsc1, zfsc2

!f--#ifdef __TIMER_SUB__
!f--  call timer_sta(1061)
!f--#endif

!      if(mod(l3,2) == 0) then
      if(mod(l,2) == 0) then
         zfsc1 => zfcos0
         if(kimg == 2) zfsc2 => zfsin0
         f = -1.d0
!      else if(mod(l3,2) == 1) then
      else if(mod(l,2) == 1) then
         zfsc1 => zfsin0
         if(kimg == 2) zfsc2 => zfcos0
         f = 1.d0
      end if

!f--#ifdef __TIMER_DO__
!f--  call timer_sta(1084)
!f--#endif
      if(kimg == 1) then
!f         do i = 1, ibl2-ibl1+1
         do i = 1, ib2-ib1+1
            veq0(i) = qitg_red(i,iq)*vlhxc_red(i,1)*zfsc1(i)
         end do
      else if(kimg == 2) then
!f         do i = 1, ibl2-ibl1+1
         do i = 1, ib2-ib1+1
            veq0(i)  = qitg_red(i,iq)*(vlhxc_red(i,1)*zfsc1(i)+f*vlhxc_red(i,2)*zfsc2(i))
         end do
      end if
!f--#ifdef __TIMER_DO__
!f--  call timer_end(1084)
!f--#endif
!f--#ifdef __TIMER_SUB__
!f--  call timer_end(1061)
!f--#endif
    end subroutine dp_Vlhxcq_exp_Q_div

    subroutine veQ_dot_ylm_div(vq,veq0,ib,il,l,ib1,ib2)
      real(kind=DP), intent(out) :: vq
      integer, intent(in) :: ib, il, l, ib1, ib2
      real(kind=DP), intent(in) :: veq0(ib)

      integer       :: i
!f--#ifdef __TIMER_SUB__
!f--  call timer_sta(1062)
!f--#endif
      vq = 0.d0
!f--#ifdef __TIMER_DO__
!f--  call timer_sta(1085)
!f--#endif
!      do i = 1, ibl2-ibl1+1
!         vq = vq + veq(i)*ylm_red(i,ilm3)
      do i = 1, ib2-ib1+1
         vq = vq + veq0(i)*ylm_red(i,il)
      end do
!f--#ifdef __TIMER_DO__
!f--  call timer_end(1085)
!f--#endif

!      if(mod(l3,4) == 0) then
      if(mod(l,4) == 0) then
         vq = vq
!      else if(mod(l3,4) == 1) then
      else if(mod(l,4) == 1) then
         vq = -vq
!      else if(mod(l3,4) == 2) then
      else if(mod(l,4) == 2) then
         vq = -vq
!      else if(mod(l3,4) == 3) then
      else if(mod(l,4) == 3) then
         vq = vq
      end if
!f--#ifdef __TIMER_SUB__
!f--  call timer_end(1062)
!f--#endif
    end subroutine veQ_dot_ylm_div

  end subroutine m_ESiVQ_integrate_VlhxcQlm_3D

!$$#endif

end module m_ES_Intgr_VlhxcQlm
