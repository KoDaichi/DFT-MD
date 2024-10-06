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
  use m_Parallelization,     only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype,ierr, np_e

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

!$$#ifndef PARA3D
  subroutine m_ESiVQ_integrate_VlhxcQlm(nfout)
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

    real(kind=DP), pointer, dimension(:)               :: ylm
    real(kind=DP), allocatable, dimension(:)           :: ylm_t
    real(kind=DP), allocatable, target, dimension(:,:) :: ylm_ext
#ifdef _VECTOR_TUNING_
    real(kind=DP), allocatable, target, dimension(:,:) :: zfcos_x, zfsin_x  ! d(ista_kngp:iend_kngp_adj,n_ialist0)
#ifndef _m_ESIV_no_loop_exchange_
    integer                                            :: ia_s
    real(kind=DP), allocatable, dimension(:,:)         :: veq_x ! d(ista_kngp:iend_kngp_adj,n_ialist0)
    real(kind=DP), allocatable, dimension(:)           :: vq_x ! d(n_ialist0)
#endif
    integer :: iend_kngp_adj, n_ialist, n_ialist0, ia_start, n_iagroup, n_ia, ia_g, l3_odd, l3_even
    integer, allocatable, dimension(:) :: ia_list
    real(kind=DP), allocatable, dimension(:,:) ::         qitgylm
#else
    real(kind=DP), allocatable, target, dimension(:)   :: zfcos, zfsin
#ifndef _m_ESIV_no_loop_exchange_
    integer :: ncache, ibsize, ibl1, ibl2
    real(kind=DP), allocatable, dimension(:)           :: veq ! d(ista_kngp:iend_kngp)
    real(kind=DP)                                      :: vq_ia
    real(kind=DP), allocatable, dimension(:,:)         :: qitg_red, ylm_red, vlhxc_red
#endif
#endif
    integer,       allocatable, dimension(:)           :: il3

    integer :: is,it,lmt1,lmt2,il1,tau1,il2,tau2,n,ilm3,l3,iiqitg,ia

#ifndef _m_ESIV_no_loop_exchange_
    integer :: m, maxm, ip, np, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
#endif

    real(kind=DP), allocatable, dimension(:,:,:,:) :: vlhxcQ_mpi ! d(nlmt,nlmt,natm,nspin)

    integer :: id_sname = -1
    if(modnrm == SKIP) return

    call tstatc0_begin('m_ESiVQ_integrate_VlhxcQlm ',id_sname)

    call m_PP_find_maximum_l(n)    ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)

#ifndef _m_ESIV_no_loop_exchange_
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
#endif

#ifdef _VECTOR_TUNING_
    n_ialist = 1
#ifdef HIUX
    n_ialist = 16
#endif
#ifdef VPP
    n_ialist = 8
#endif
#ifdef SX
    n_ialist = 8
#endif
    iend_kngp_adj = iend_kngp
#ifdef HIUX
#else
    if(mod(iend_kngp-ista_kngp+1,2) == 0) iend_kngp_adj = iend_kngp_adj+1
#endif
    allocate(ia_list(n_ialist))
#endif

    if(n**2 > nel_Ylm) then
       allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2));  ylm_ext = 0.d0
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)  ! (ilm3,rltv,ngabc,gr_l)->(ylm)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
       deallocate(ylm_t)
    end if

    vlhxcQ = 0.d0

! ====================================== modified by K. Tagami ============== 11.0
!    do is = 1, nspin, af+1
    do is = 1, ndim_magmom, af+1
! =========================================================================== 11.0

#ifdef _VECTOR_TUNING_
       do it = 1, ntyp
          if(sum(ivanl(:,it)) == 0) cycle
          allocate(qitgylm(ista_kngp:iend_kngp,2))

          n_ia = 0
          do ia = 1, natm
             if(ityp(ia) == it) n_ia = n_ia + 1
          end do
          if(n_ialist <=0) call phase_error_with_msg(nfout, 'n_ialist is illegal <<m_ESiVQ_integrate_VlhxcQlm>>'&
                                                    ,__LINE__,__FILE__)
          n_iagroup = n_ia/n_ialist + 1
          ia_start = 1
          if(iprivlhxcq >= 2) write(nfout,'(" !m_ESiVQ_integrate_VlhxcQlm: n_iagroup = ",i8, " ityp = ",i8)') n_iagroup,it
          do ia_g = 1, n_iagroup
             n_ialist0 = 0
             ia_list = 0
             AtomcountLoop: do ia = ia_start, natm
                if(ityp(ia) == it) then
                   n_ialist0 = n_ialist0 + 1
                   ia_list(n_ialist0) = ia
                end if
                if(n_ialist0 >= n_ialist) exit AtomcountLoop
             end do AtomcountLoop
             ia_start = ia+1
             if(n_ialist0 >= 1) then
                if(iprivlhxcq >= 2) write(nfout,'(" !m_ESiVQ_integrate_VlhxcQlm: ia_list = ",8i8)') (ia_list(ia),ia=1,n_ialist0)
                allocate(zfcos_x(ista_kngp:iend_kngp_adj,n_ialist0))
                allocate(zfsin_x(ista_kngp:iend_kngp_adj,n_ialist0))
                call calc_phase_b(natm,pos,ia_list,n_ialist0,kgp,ngabc,ista_kngp,iend_kngp &
                     &                                                ,ista_kngp,iend_kngp_adj,zfcos_x,zfsin_x)
#ifndef _m_ESIV_no_loop_exchange_
                allocate(veq_x(ista_kngp:iend_kngp_adj,n_ialist0))
                allocate(vq_x(n_ialist0))
                if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ia,  iq,  l3,   m,ilm3,  ip,lmt1,lmt2,  np")')
                do iq = nqitg_sp0(it), nqitg_sp(it)
                   l3 = iq2l3(iq)
                   if(l3 == 0) then
                      ilm3 = 1
                      call integration_of_VlhxcQlm_core4(iq,vq_x) ! iq, vlhxc,exp,qitg_l,ylm -> vq_x
                   else
                      call dp_Vlhxcq_exp_Q_x(iq) ! vlhxc, exp, qitg_l  -> veq_x
                   end if
                   do m = 1, 2*l3+1
                      ilm3 = l3*l3+m
                      if(l3 /= 0) call veQ_dot_ylm_x(vq_x) ! veq_x, ylm -> vq_x
                      do ia_s = 1, n_ialist0
                         ia = ia_list(ia_s)
                         do ip = 1, nc(m,iq)
                            lmt1 = nc2lmt1(ip,m,iq)
                            lmt2 = nc2lmt2(ip,m,iq)
                            np = nc2n(ip,m,iq)
                            if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ: ",9i5)') ia, iq, l3, m, ilm3, ip, lmt1,lmt2,np
                            vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is) + univol*vq_x(ia_s)*dl2p(lmt1,lmt2,np,it)
                         end do
                      end do
                   end do
                end do
                deallocate(vq_x,veq_x)
#else
                if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ia,lmt1, il1,tau1,lmt2, il2,tau2,il2p,ilm3,  l3,iiqitg ")')

                do lmt1 = 1, ilmt(it)
                   il1  = ltp( lmt1, it);   tau1 = taup(lmt1, it)
                   do lmt2 = lmt1, ilmt(it)
                      il2  = ltp( lmt2, it);   tau2 = taup(lmt2, it)
                      l3_odd = 0; l3_even = 0

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
                      if ( il1 /= il2 ) cycle
                      if ( im1 /= im2 ) cycle
#endif
! ---------------------------------------------------- 11.0S


                      qitgylm = 0.d0
                      do n = 1, il2p(lmt1,lmt2,it)
                         ilm3 = isph(lmt1,lmt2,n,it);   l3 = il3(ilm3)
                         iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                         if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ: ",11i5)') ia,lmt1,il1,tau1,lmt2,il2,tau2,n,ilm3,l3,iiqitg
                         if(iiqitg == 0) cycle
                         call integrate_qitg_x_ylm() ! -> qitgylm
                      end do
                      if(l3_odd + l3_even >=1) call integration_of_VlhxcQlm_core3()
                      !                             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   end do
                end do
#endif
                deallocate(zfsin_x,zfcos_x)
             end if
          end do
          deallocate(qitgylm)
       end do
#else

#ifndef _m_ESIV_no_loop_exchange_
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
          if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ibl1,ibl2,ibsize = ", 3i8)') ibl1, ibl2, ibsize

          call substitute_qitgred()  ! qitg_l -> qitg_red
          call substitute_vlhxcred() ! vlhxc_l -> vlhxc_red
          call substitute_ylmred()   ! ylm_l, ylm_ext -> ylm_red
          do ia = 1, natm
             it = ityp(ia)
             if(sum(ivanl(:,it)) == 0) cycle
             call calc_phase_div(ia)  !   -> zfcos, zfsin
             if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ia,  iq,  l3,   m,ilm3,  ip,lmt1,lmt2,  np")')
             do iq = nqitg_sp0(it), nqitg_sp(it)
                l3 = iq2l3(iq)
                call dp_Vlhxcq_exp_Q_div(iq) ! vlhxc, exp, qitg_l  -> veq
                do m = 1, 2*l3+1
                   ilm3 = l3*l3+m
!!$                   if(ilm3 <= nel_Ylm) then
!!$                      ylm => ylm_l(:,ilm3)
!!$                   else
!!$                      ylm => ylm_ext(:,ilm3)
!!$                   end if
                   call veQ_dot_ylm_div(vq_ia) ! veq, ylm -> vq_ia
                   do ip = 1, nc(m,iq)
                      lmt1 = nc2lmt1(ip,m,iq)
                      lmt2 = nc2lmt2(ip,m,iq)
                      np = nc2n(ip,m,iq)
                      if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ: ",9i5)') ia, iq, l3, m, ilm3, ip, lmt1,lmt2,np
                      vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is) + univol*vq_ia*dl2p(lmt1,lmt2,np,it)
                   end do
                end do
             end do
          end do
       end do
       deallocate(ylm_red)
       deallocate(vlhxc_red)
       deallocate(qitg_red)
       deallocate(veq)
#else
       allocate(zfcos(ista_kngp:iend_kngp)); zfcos = 0.d0
       allocate(zfsin(ista_kngp:iend_kngp)); zfsin = 0.d0
       do ia = 1, natm
          it = ityp(ia)
          if(sum(ivanl(:,it)) == 0) cycle
          call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos,zfsin)
          !                -> zfcos, zfsin
          if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ:    ia,lmt1, il1,tau1,lmt2, il2,tau2,il2p,ilm3,  l3,iiqitg ")')
          do lmt1 = 1, ilmt(it)
             il1  = ltp( lmt1, it);   tau1 = taup(lmt1, it)
             do lmt2 = lmt1, ilmt(it)
                il2  = ltp( lmt2, it);   tau2 = taup(lmt2, it)

#ifdef SKIP_TEST
! ---------------------------------------------------- 11.0S
                if ( il1 /= il2 ) cycle
                if ( im1 /= im2 ) cycle
! ---------------------------------------------------- 11.0S
#endif

                do n = 1, il2p(lmt1,lmt2,it)
                   ilm3 = isph(lmt1,lmt2,n,it);   l3 = il3(ilm3)
                   iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                   if(iprivlhxcq >= 2) write(nfout,'(" !mESiVQ: ",11i5)') ia,lmt1,il1,tau1,lmt2,il2,tau2,n,ilm3,l3,iiqitg
                   if(iiqitg == 0) cycle
                   if(ilm3 <= nel_Ylm) then
                      ylm => ylm_l(:,ilm3)
                   else
                      ylm => ylm_ext(:,ilm3)
                   end if
                   call integration_of_VlhxcQlm_core2 ! -> vlhxcQ(lmt1,lmt2,ia,is)
                   !    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                end do
             end do
          end do
       end do
#endif
#endif

#ifndef _VECTOR_TUNING_
       deallocate(zfsin);deallocate(zfcos)
#endif
    end do
    if(npes > 1) then
! ================================= modified by K. Tagami ================ 11.0
!       allocate(vlhxcQ_mpi(nlmt,nlmt,natm,nspin)); vlhxcQ_mpi = 0.d0
!       call mpi_allreduce(vlhxcQ, vlhxcQ_mpi,nlmt*nlmt*natm*nspin, mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

       allocate( vlhxcQ_mpi(nlmt,nlmt,natm,ndim_magmom) ); vlhxcQ_mpi = 0.d0
       call mpi_allreduce( vlhxcQ, vlhxcQ_mpi, nlmt*nlmt*natm*ndim_magmom, &
        &                  mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
! ======================================================================= 11.0
       vlhxcQ = vlhxcQ_mpi
       deallocate(vlhxcQ_mpi)
    end if
    call symmetrize_vlhxcq
    if(allocated(ylm_ext)) deallocate(ylm_ext)
#ifdef _VECTOR_TUNING_
    deallocate(ia_list); deallocate(il3)
#else
!!$    deallocate(zfsin); deallocate(zfcos); deallocate(il3)
    deallocate(il3)
#endif

#ifndef _m_ESIV_no_loop_exchange_
    deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif

    if(iprivlhxcq >= 2) then
       write(nfout,*) ' <<< m_ESiVQ_integrate_VlhxcQlm >>> vlhxcQ start'
! =================================== modified by K. Tagami =============== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ========================================================================= 11.0

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
    call tstatc0_end(id_sname)
  contains

    subroutine symmetrize_vlhxcq
      integer :: is, ia, it, lmt1, lmt2
! =================================== modified by K. Tagami =============== 11.0
!       do is = 1, nspin, af+1
       do is = 1, ndim_magmom, af+1
! ========================================================================= 11.0

         do ia = 1, natm
            it = ityp(ia)
            do lmt1 = 1, ilmt(it)
               do lmt2 = lmt1+1, ilmt(it)
                  vlhxcQ(lmt2,lmt1,ia,is) = vlhxcQ(lmt1,lmt2,ia,is)
               end do
            end do
         end do
      end do
    end subroutine symmetrize_vlhxcq

#ifdef _VECTOR_TUNING_
#ifndef _m_ESIV_no_loop_exchange_
    subroutine integration_of_VlhxcQlm_core4(iq,vq)
      integer, intent(in) :: iq
      real(kind=DP),dimension(n_ialist0), intent(out) :: vq
      integer       :: i,ia_s,ia
      real(kind=DP) ::  f


      vq = 0.d0
      if(ilm3 <= nel_Ylm) then
         if(kimg == 1) then
            if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_l(i,ilm3)*vlhxc_l(i,1,is)*zfcos_x(i,ia_s)
                  end do
               end do
            else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_l(i,ilm3)*vlhxc_l(i,1,is)*zfsin_x(i,ia_s)
                  end do
               end do
            end if
         else if(kimg == 2) then
            if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_l(i,ilm3) &
                          & *(vlhxc_l(i,1,is)*zfcos_x(i,ia_s)-vlhxc_l(i,2,is)*zfsin_x(i,ia_s))
                  end do
               end do
            else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  !for mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_l(i,ilm3) &
                          & *(vlhxc_l(i,1,is)*zfsin_x(i,ia_s)+vlhxc_l(i,2,is)*zfcos_x(i,ia_s))
                  end do
               end do
            end if
         end if
      else if(ilm3 > nel_Ylm) then
         if(kimg == 1) then
            if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_ext(i,ilm3)*vlhxc_l(i,1,is)*zfcos_x(i,ia_s)
                  end do
               end do
            else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_ext(i,ilm3)*vlhxc_l(i,1,is)*zfsin_x(i,ia_s)
                  end do
               end do
            end if
         else if(kimg == 2) then
            if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  ! mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_ext(i,ilm3) &
                          & *(vlhxc_l(i,1,is)*zfcos_x(i,ia_s)-vlhxc_l(i,2,is)*zfsin_x(i,ia_s))
                  end do
               end do
            else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
               do ia_s = 1, n_ialist0
                  do i = ista_kngp, iend_kngp  !for mpi
                     vq(ia_s) = vq(ia_s) + qitg_l(i,iq)*ylm_ext(i,ilm3) &
                          & *(vlhxc_l(i,1,is)*zfsin_x(i,ia_s)+vlhxc_l(i,2,is)*zfcos_x(i,ia_s))
                  end do
               end do
            end if
         end if
      end if

      if(mod(l3,2) == 0) then
         vq(:) = real(zi**l3)*vq(:)
      else if(mod(l3,2) == 1) then
         vq(:) = -aimag(zi**l3)*vq(:)
      end if
    end subroutine integration_of_VlhxcQlm_core4

    subroutine dp_Vlhxcq_exp_Q_x(iq)
      integer, intent(in) :: iq
      integer       :: i, ia_s

      if(kimg == 1) then
         if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  veq_x(i,ia_s) = qitg_l(i,iq)*vlhxc_l(i,1,is)*zfcos_x(i,ia_s)
               end do
            end do
         else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  veq_x(i,ia_s) = qitg_l(i,iq)*vlhxc_l(i,1,is)*zfsin_x(i,ia_s)
               end do
            end do
         end if
      else if(kimg == 2) then
         if(mod(l3,2) == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  veq_x(i,ia_s) = qitg_l(i,iq)*(vlhxc_l(i,1,is)*zfcos_x(i,ia_s)-vlhxc_l(i,2,is)*zfsin_x(i,ia_s))
               end do
            end do
         else if(mod(l3,2) == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  !for mpi
                  veq_x(i,ia_s) = qitg_l(i,iq)*(vlhxc_l(i,1,is)*zfsin_x(i,ia_s)+vlhxc_l(i,2,is)*zfcos_x(i,ia_s))
               end do
            end do
         end if
      end if
    end subroutine dp_Vlhxcq_exp_Q_x

    subroutine veQ_dot_ylm_x(vq_x)
      real(kind=DP), intent(out) :: vq_x(n_ialist0)

      integer       :: i,iy, ia

      if(ilm3 <= nel_Ylm) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
#ifdef SX
!CDIR PARALLEL DO PRIVATE(ia,i)
#endif
         do ia = 1, n_ialist0
            vq_x(ia) = 0.d0
            do i = ista_kngp, iend_kngp  ! mpi
               vq_x(ia) = vq_x(ia) + veq_x(i,ia)*ylm_l(i,ilm3)
            end do
         end do
      else
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
#ifdef SX
!CDIR PARALLEL DO PRIVATE(ia,i)
#endif
         do ia = 1, n_ialist0
            vq_x(ia) = 0.d0
            do i = ista_kngp, iend_kngp  ! mpi
               vq_x(ia) = vq_x(ia) + veq_x(i,ia)*ylm_ext(i,ilm3)
            end do
         end do
      end if

      if(mod(l3,4) == 0) then
         vq_x = vq_x
      else if(mod(l3,4) == 1) then
         vq_x = -vq_x
      else if(mod(l3,4) == 2) then
         vq_x = -vq_x
      else if(mod(l3,4) == 3) then
         vq_x = vq_x
      end if
    end subroutine veQ_dot_ylm_x
#else
    subroutine integrate_qitg_x_ylm()
      integer :: i
      real(kind=DP) :: f

      if(ilm3 <= nel_Ylm) then
         if(mod(l3,2) == 0) then
            l3_even = l3_even+1
            f = real(zi**l3)*dl2p(lmt1,lmt2,n,it)
            do i = ista_kngp, iend_kngp  ! mpi
               qitgylm(i,1) = qitgylm(i,1) + f*qitg_l(i,iiqitg)*ylm_l(i,ilm3)
            end do
         else if(mod(l3,2) == 1) then
            l3_odd = l3_odd+1
            f = -aimag(zi**l3)*dl2p(lmt1,lmt2,n,it)
            do i = ista_kngp, iend_kngp  !mpi
               qitgylm(i,2) = qitgylm(i,2) + f*qitg_l(i,iiqitg)*ylm_l(i,ilm3)
            end do
         end if
      else
         if(mod(l3,2) == 0) then
            l3_even = l3_even+1
            f = real(zi**l3)*dl2p(lmt1,lmt2,n,it)
            do i = ista_kngp, iend_kngp  ! mpi
               qitgylm(i,1) = qitgylm(i,1) + f*qitg_l(i,iiqitg)*ylm_ext(i,ilm3)
            end do
         else if(mod(l3,2) == 1) then
            l3_odd = l3_odd+1
            f = -aimag(zi**l3)*dl2p(lmt1,lmt2,n,it)
            do i = ista_kngp, iend_kngp  !mpi
               qitgylm(i,2) = qitgylm(i,2) + f*qitg_l(i,iiqitg)*ylm_ext(i,ilm3)
            end do
         end if
      end if
    end subroutine integrate_qitg_x_ylm

    subroutine integration_of_VlhxcQlm_core3
      integer       :: i,ia_s,ia
      real(kind=DP),allocatable,dimension(:) :: vq

      allocate(vq(n_ialist0)); vq = 0.d0

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
         do ia_s = 1, n_ialist0
            do i = ista_kngp, iend_kngp  ! mpi
               vq(ia_s) = vq(ia_s) + vlhxc_l(i,1,is)*(qitgylm(i,1)*zfcos_x(i,ia_s) &
                    &                               + qitgylm(i,2)*zfsin_x(i,ia_s))
            end do
         end do
      else if(kimg == 2) then
         if(l3_even > 0 .and. l3_odd >0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  vq(ia_s) = vq(ia_s) &
                       & + qitgylm(i,1)*(vlhxc_l(i,1,is)*zfcos_x(i,ia_s)-vlhxc_l(i,2,is)*zfsin_x(i,ia_s)) &
                       & + qitgylm(i,2)*(vlhxc_l(i,1,is)*zfsin_x(i,ia_s)+vlhxc_l(i,2,is)*zfcos_x(i,ia_s))
               end do
            end do
        else if(l3_even > 0 .and. l3_odd == 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  vq(ia_s) = vq(ia_s) &
                       & + qitgylm(i,1)*(vlhxc_l(i,1,is)*zfcos_x(i,ia_s)-vlhxc_l(i,2,is)*zfsin_x(i,ia_s))
               end do
            end do
         else if(l3_even == 0 .and. l3_odd > 0) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef HIUX
*poption indep
*poption parallel
#endif
            do ia_s = 1, n_ialist0
               do i = ista_kngp, iend_kngp  ! mpi
                  vq(ia_s) = vq(ia_s) &
                       & + qitgylm(i,2)*(vlhxc_l(i,1,is)*zfsin_x(i,ia_s)+vlhxc_l(i,2,is)*zfcos_x(i,ia_s))
               end do
            end do
         end if
      end if

      do ia_s = 1, n_ialist0
         ia = ia_list(ia_s)
         vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is)+univol*vq(ia_s)
      end do
      deallocate(vq)
    end subroutine integration_of_VlhxcQlm_core3
#endif


!!!-->  #else if .not. _VECTOR_TUNING_
#else

    subroutine substitute_qitgred()
      integer :: iq, i
      do iq = 1, nqitg
         do i = 1, ibl2-ibl1+1
            qitg_red(i, iq) = qitg_l(i+ibl1-1,iq)
         end do
      end do
    end subroutine substitute_qitgred

    subroutine substitute_vlhxcred
      integer :: i
      if(kimg==1) then
         do i = 1, ibl2-ibl1+1
            vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
         end do
      else if(kimg==2) then
         do i = 1, ibl2-ibl1+1
            vlhxc_red(i,1) = vlhxc_l(i+ibl1-1,1,is)
            vlhxc_red(i,2) = vlhxc_l(i+ibl1-1,2,is)
         end do
      end if
    end subroutine substitute_vlhxcred

    subroutine substitute_ylmred
      integer :: ilm, i
      do ilm = 1, nel_Ylm
         do i = 1, ibl2-ibl1+1
            ylm_red(i,ilm) = ylm_l(i+ibl1-1,ilm)
         end do
      end do
      if(n**2 > nel_Ylm) then
         do ilm = nel_ylm+1, n**2
            do i = 1, ibl2-ibl1+1
               ylm_red(i,ilm) = ylm_ext(i+ibl1-1,ilm)
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

#ifndef _m_ESIV_no_loop_exchange_
    subroutine dp_Vlhxcq_exp_Q(iq)
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f
      real(kind=DP), pointer, dimension(:) :: zfsc1, zfsc2

      if(mod(l3,2) == 0) then
         zfsc1 => zfcos
         if(kimg == 2) zfsc2 => zfsin
         f = -1.d0
      else if(mod(l3,2) == 1) then
         zfsc1 => zfsin
         if(kimg == 2) zfsc2 => zfcos
         f = 1.d0
      end if

      if(kimg == 1) then
         do i = ista_kngp, iend_kngp  ! mpi
            veq(i) = qitg_l(i,iq)*vlhxc_l(i,1,is)*zfsc1(i)
         end do
      else if(kimg == 2) then
         do i = ista_kngp, iend_kngp  ! mpi
            veq(i)  = qitg_l(i,iq)*(vlhxc_l(i,1,is)*zfsc1(i)+f*vlhxc_l(i,2,is)*zfsc2(i))
         end do
      end if
    end subroutine dp_Vlhxcq_exp_Q

    subroutine veQ_dot_ylm(vq)
      real(kind=DP), intent(out) :: vq

      integer       :: i,iy
      vq = 0.d0
      do i = ista_kngp, iend_kngp  ! mpi
         iy = i-ista_kngp+1
         vq = vq + veq(i)*ylm(iy)
      end do

      if(mod(l3,4) == 0) then
         vq = vq
      else if(mod(l3,4) == 1) then
         vq = -vq
      else if(mod(l3,4) == 2) then
         vq = -vq
      else if(mod(l3,4) == 3) then
         vq = vq
      end if
    end subroutine veQ_dot_ylm

    subroutine dp_Vlhxcq_exp_Q_div(iq)
      integer, intent(in) :: iq
      integer       :: i
      real(kind=DP) :: f
      real(kind=DP), pointer, dimension(:) :: zfsc1, zfsc2

      if(mod(l3,2) == 0) then
         zfsc1 => zfcos
         if(kimg == 2) zfsc2 => zfsin
         f = -1.d0
      else if(mod(l3,2) == 1) then
         zfsc1 => zfsin
         if(kimg == 2) zfsc2 => zfcos
         f = 1.d0
      end if

      if(kimg == 1) then
         do i = 1, ibl2-ibl1+1
            veq(i) = qitg_red(i,iq)*vlhxc_red(i,1)*zfsc1(i)
         end do
      else if(kimg == 2) then
         do i = 1, ibl2-ibl1+1
            veq(i)  = qitg_red(i,iq)*(vlhxc_red(i,1)*zfsc1(i)+f*vlhxc_red(i,2)*zfsc2(i))
         end do
      end if
    end subroutine dp_Vlhxcq_exp_Q_div

    subroutine veQ_dot_ylm_div(vq)
      real(kind=DP), intent(out) :: vq

      integer       :: i,iy
      vq = 0.d0
      do i = 1, ibl2-ibl1+1
         vq = vq + veq(i)*ylm_red(i,ilm3)
      end do

      if(mod(l3,4) == 0) then
         vq = vq
      else if(mod(l3,4) == 1) then
         vq = -vq
      else if(mod(l3,4) == 2) then
         vq = -vq
      else if(mod(l3,4) == 3) then
         vq = vq
      end if
    end subroutine veQ_dot_ylm_div

#else
    subroutine integration_of_VlhxcQlm_core2
      integer       :: i,iy
      real(kind=DP) :: vq
      real(kind=DP) :: f
      real(kind=DP), pointer, dimension(:) :: zfsc1, zfsc2

      vq = 0.d0
      if(mod(l3,2) == 0) then
         zfsc1 => zfcos
         if(kimg == 2) zfsc2 => zfsin
         f = -1.d0
      else if(mod(l3,2) == 1) then
         zfsc1 => zfsin
         if(kimg == 2) zfsc2 => zfcos
         f = 1.d0
      end if

      if(kimg == 1) then
         do i = ista_kngp, iend_kngp  ! mpi
            iy = i-ista_kngp+1
            vq = vq + qitg_l(i,iiqitg)*ylm(iy)*vlhxc_l(i,1,is)*zfsc1(i)
         end do
      else if(kimg == 2) then
         do i = ista_kngp, iend_kngp  ! mpi
            iy = i-ista_kngp+1
            vq  = vq + qitg_l(i,iiqitg)*ylm(iy)*(vlhxc_l(i,1,is)*zfsc1(i)+f*vlhxc_l(i,2,is)*zfsc2(i))
         end do
      end if
      if(mod(l3,2) == 0) then
         vq = real(zi**l3)*dl2p(lmt1,lmt2,n,it)*vq
      else if(mod(l3,2) == 1) then
         vq = -aimag(zi**l3)*dl2p(lmt1,lmt2,n,it)*vq
      end if

      vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is)+univol*vq
    end subroutine integration_of_VlhxcQlm_core2
#endif
#endif

!!$    subroutine integration_of_VlhxcQlm_core
!!$      integer       :: ia, i, iy
!!$      real(kind=DP) :: dga, vq, grt, qy, zfc, zfs
!!$
!!$      dga = dl2p(lmt1,lmt2,n,it)
!!$      do ia = 1, natm
!!$         if(ityp(ia) /= it) cycle
!!$         vq = 0.d0
!!$         if(mod(il1+il2,2) == 0) then
!!$            do i = ista_kngp, iend_kngp  !for mpi
!!$               iy = i-ista_kngp+1
!!$               grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2)&
!!$                    &+pos(ia,3)*ngabc(i,3))*PAI2
!!$               zfc = dcos(grt); if(kimg == 2) zfs = dsin(grt)
!!$               qy  = qitg_l(i,iiqitg)*ylm(iy)
!!$               vq  = vq+ qy*vlhxc_l(i,1,is)*zfc
!!$               if(kimg==2) vq = vq - qy*vlhxc_l(i,kimg,is)*zfs
!!$            end do
!!$            vq = real(zi**l3)*dga*vq
!!$         else if(mod(il1+il2,2) == 1) then
!!$            do i = ista_kngp, iend_kngp  !for mpi
!!$               iy = i-ista_kngp+1
!!$               grt = (pos(ia,1)*ngabc(i,1) + pos(ia,2)*ngabc(i,2)&
!!$                    &+pos(ia,3)*ngabc(i,3))*PAI2
!!$               zfs = dsin(grt); if(kimg==2) zfc = dcos(grt)
!!$               qy  = qitg_l(i,iiqitg)*ylm(iy)
!!$               vq  = vq + qy*vlhxc_l(i,1,is)*zfs
!!$               if(kimg==2) vq = vq + qy*vlhxc_l(i,kimg,is)*zfc
!!$            end do
!!$            vq = -aimag(zi**l3)*dga*vq
!!$         endif
!!$         vlhxcQ(lmt1,lmt2,ia,is) = vlhxcQ(lmt1,lmt2,ia,is)+univol*vq
!!$      end do
!!$    end subroutine integration_of_VlhxcQlm_core
  end subroutine m_ESiVQ_integrate_VlhxcQlm
!$$#endif

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


end module m_ES_Intgr_VlhxcQlm
