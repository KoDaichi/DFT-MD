!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_Stress
!
!  AUTHOR(S): H. Sawada   August/20/2003
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
! *************************************************************
!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.0B]
!                 stress terms are avaiable in the case of DFT-D2 calculations
!
! =============================================================

#ifndef NO_NONLOCAL_DGEMM
#define NONLOCAL_DGEMM
#endif

module m_Stress

!$$#ifndef PARA3D

!   This module is for the calculation of stress tensor.
!   Original program was coded with fortan77 by H. Sawada in 1997,
!  and translated into this module with fortran90 and mpi by H. Sawada
!  and T. Yamasaki in 1999.
!
! $Id: m_Stress.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters,     only : DP, PAI, PAI2, PAI4, SKIP, zi &
                                  &, Valence_plus_PC_Charge &
!!$ASASASASAS
!!$                                  &, EXC_ONLY,STRESS_, CRDTYP, GAMMA
                                  &, EXC_ONLY,STRESS_, CRDTYP, GAMMA &
                                  &, MAPPED,NOTMAPPED, ON
!!$ASASASASAS
  use m_Crystal_Structure,    only : altv,rltv,univol,nopr,op, p2bmat, b2pmat
  use m_Ionic_System,         only : ntyp, natm,ityp,iwei,pos,zfm3_l,s_ew
  use m_PlaneWaveBasisSet,    only : kg1,kg,kgp,ngabc,iba,nbase,nbmx,gr_l,ylm_l &
                                  &, m_pwBS_sphrp2,m_pwBS_sphrp2_diff
  use m_Electronic_Structure, only : zaj_l,occup_l,eko_l,totch &
                                  &, fsr_l,fsi_l
!!$ASASASASAS
!!$                                  &, m_ES_alloc_arai, m_ES_dealloc_arai
!!$ASASASASAS
  use m_ES_nonlocal,          only :                                            &
 &                                   m_ES_betar_dot_WFs_4_lmta_k                &
 &                                 , G_dot_R_map,G_dot_R_mpi                    &
 &                                 , m_ES_alloc_zfsincos, m_ES_dealloc_zfsincos &
 &                                 , m_ES_alloc_arai, m_ES_dealloc_arai         &
 &                                 , alloc_zfsincos_mpi, dealloc_zfsincos_mpi
#ifdef NONLOCAL_DGEMM
  use m_ES_nonlocal,          only : m_ES_betar_dot_Psi_4_each_k_snl
#endif
  use m_Kpoints,              only : k_symmetry
  use m_PseudoPotential,      only : nlmt,nlmta,lmta,ilmt,lmtt,ltp,mtp,il2p,isph &
                                  &, iqitg,taup,dl2p &
                                  &, q,dion,psc_l,psc_diff_l &
                                  &, qitg_l,qitg_diff_l, modnrm &
                                  &, itpcc,rhpcg_l,rhpcg_diff_l,epc,etot1,ival &
                                  &, m_PP_find_maximum_l, m_PP_include_vanderbilt_pot &
                                  &, ipaw,dion_paw, flg_paw &
                                  &, nqitg &
                                  &, m_PP_set_index_arrays1 &
                                  &, m_PP_set_index_arrays2
  use m_NonLocal_Potential,   only : snld
  use m_XC_Potential,         only : vxc_l, vxcpc_l, exc, s_gga1, s_gga2

#ifndef DISABLE_VDWDF
  use m_vdWDF,  only : s_cnl1, s_cnl2, s_cnl1_pc, s_cnl2_pc, ecnl_vdwdf
#endif

  use m_Charge_Density,       only : chgq_l, hsr, chgsoft
  use m_Control_Parameters,   only : nspin,ipri,kimg,af,nel_Ylm,sw_screening_correction,nhistory_stress, printable&
                                  &, sw_optimize_lattice,sw_uniform &
 &                                 , ke_e0,ke_sigma,ke_a, xctype &
 &                                 , m_CtrlP_cachesize
  use m_Timing,               only : tstatc0_begin, tstatc0_end
  use m_Parallelization,      only : MPI_CommGroup,map_k,myrank_k,np_e,ierr,ista_k,iend_k &
                                  &, ista_kngp,iend_kngp,ista_snl,iend_snl,npes,mype &
                                  &, ista_spin, iend_spin
  use m_Screening,            only : screening

! =================================== Added by K. Tagami ============== p1100
  use m_Electronic_Structure,  only : dhub
  use m_Ionic_System,          only : ihubbard
  use m_Control_Parameters,    only : proj_attribute
! ===================================================================== p1100

  use m_IterationNumbers,      only : iteration_unit_cell

! ================================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor, ndim_magmom, ndim_chgpot
  use m_PseudoPotential,     only : dion0_noncl, q_noncl
  use m_ES_NonCollinear,       only : m_ES_MagMom_To_DensMat_hsr
  use m_Electronic_Structure,  only : dhub_aimag
! =================================================================== 11.0

! =================================== KT_add ============ 13.0B &13.1AS
  use m_Control_Parameters,    only : external_stress, &
       &                              sw_fix_lattice_angles, sw_vdw_correction, &
       &                              sw_fix_lattice_lengths, &
       &                              fix_angle_alpha, fix_angle_beta, fix_angle_gamma, &
       &                              fix_length_a, fix_length_b, fix_length_c, &
       &                              sw_neglect_stress_offdiagonal, &
       &                              sw_fix_lattice_shapes, fix_shape_ab_plane, &
       &                              fix_shape_bc_plane, fix_shape_ac_plane, &
       &                              estress_has_been_set, sw_add_vlocal_gzero
  use m_Ionic_System,          only : s_vdw
! ======================================================= 13.0B &13.1AS

  use m_Charge_Density,      only : m_CD_symmtrz_of_ff_noncl_C
  use mpi

  implicit none
!  include 'mpif.h'

  integer,      private,pointer, dimension(:)   :: il3
  real(kind=DP),private,pointer, dimension(:)   :: vlength,wkc
  real(kind=DP),private,pointer, dimension(:)   :: grinv,zfsin_tmp,zfcos_tmp,mzfsin_tmp
  real(kind=DP),private,pointer, dimension(:)   :: flchgq
  real(kind=DP),private,pointer, dimension(:,:) :: g,alinvt,vloc
  real(kind=DP),private,pointer, dimension(:,:) :: s_ke,s_or,s_nl
  real(kind=DP),private,pointer, dimension(:,:) :: s_h1,s_h2,s_loc1,s_loc2
  real(kind=DP),private,pointer, dimension(:,:) :: s_xc1,s_xc2
  real(kind=DP),private,pointer, dimension(:,:) :: s_pc
  real(kind=DP),private,allocatable, dimension(:,:,:) :: fsr_diff_l,fsi_diff_l,ylmd
  real(kind=DP),private,pointer, dimension(:)   :: ylm
  real(kind=DP),private,pointer, dimension(:,:,:) :: flchgqd
  real(kind=DP),public, pointer, dimension(:,:,:,:,:,:) :: hsrd,esrd

! ========================= added by K. Tagami =========================== 11.0
  real(kind=DP),public, allocatable, dimension(:,:,:,:,:,:) :: hsid, esid
! ======================================================================== 11.0

  real(kind=DP),private,pointer, dimension(:,:) :: s_0,s_1,s_2
  real(kind=DP),private,pointer, dimension(:,:) :: s_tmp,s_up,s_dn,s
  character(len=2) ,private                     :: stress_name

! ====================== Added by K. Tagami ========================= p1100
  real(kind=DP),private,allocatable, dimension(:,:) :: s_hub
! =================================================================== p1100
  real(kind=DP),private,allocatable, dimension(:,:) :: s_vdwdf

  real(kind=DP),private :: stressmx=-1
  real(kind=DP),dimension(3,3),private :: currstress=0.d0

  real(kind=DP), dimension(3), private :: etot3=0.d0
  real(kind=DP), private :: univol0
  real(kind=DP) :: stress_correction = 0.d0
  character(len('stress_correction')), parameter, private :: tag_stress_correction='stress_correction'
  real(kind=DP), allocatable, target, dimension(:,:) :: ylm_red,ylm_ext
  integer :: nblk = 1000

  interface m_Stress_in_correction
    module procedure in_correction0
    module procedure in_correction1
  end interface m_Stress_in_correction
contains
  subroutine m_Stress_alloc
    allocate(s_ke(3,3)); s_ke = 0.d0
    allocate(s_or(3,3)); s_or = 0.d0
    allocate(s_nl(3,3)); s_nl = 0.d0
    allocate(s_xc1(3,3)); s_xc1 = 0.d0
    allocate(s_xc2(3,3)); s_xc2 = 0.d0
    allocate(s_h1 (3,3)); s_h1  = 0.d0
    allocate(s_h2 (3,3)); s_h2  = 0.d0
    allocate(s_loc1(3,3)); s_loc1 = 0.d0
    allocate(s_loc2(3,3)); s_loc2 = 0.d0
    allocate(s_pc(3,3)); s_pc = 0.d0

! ======================= modified by K. Tagami ================ 11.0
!    allocate(hsrd(natm,nlmt,nlmt,nspin,3,3)); hsrd = 0.d0
!    allocate(esrd(natm,nlmt,nlmt,nspin,3,3)); esrd = 0.d0
!
    if ( noncol ) then
       allocate(hsrd(natm,nlmt,nlmt,ndim_magmom,3,3)); hsrd = 0.d0
       allocate(esrd(natm,nlmt,nlmt,ndim_magmom,3,3)); esrd = 0.d0
       allocate(hsid(natm,nlmt,nlmt,ndim_magmom,3,3)); hsid = 0.d0
       allocate(esid(natm,nlmt,nlmt,ndim_magmom,3,3)); esid = 0.d0
    else
       allocate(hsrd(natm,nlmt,nlmt,nspin,3,3)); hsrd = 0.d0
       allocate(esrd(natm,nlmt,nlmt,nspin,3,3)); esrd = 0.d0
    endif
! ============================================================== 11.0
    allocate(s_0(3,3)); s_0 = 0.d0
    allocate(s_1(3,3)); s_1 = 0.d0
    allocate(s_2(3,3)); s_2 = 0.d0
    allocate(s_tmp(3,3)); s_tmp = 0.d0
    allocate(s_up(3,3)); s_up = 0.d0
    allocate(s_dn(3,3)); s_dn = 0.d0
    allocate(s(3,3)); s = 0.d0
    allocate(alinvt(3,3)); alinvt = rltv / PAI2

! =========================================== Added by K.Tagami ==== p1100
    allocate(s_hub(3,3)); s_hub = 0.d0
! ================================================================== p1100
    allocate(s_vdwdf(3,3)); s_vdwdf = 0.d0

  end subroutine m_Stress_alloc

  subroutine m_Stress_dealloc
    deallocate(s_ke)
    deallocate(s_or)
    deallocate(s_nl)
    deallocate(s_xc1)
    deallocate(s_xc2)
    deallocate(s_h1)
    deallocate(s_h2)
    deallocate(s_loc1)
    deallocate(s_loc2)
    deallocate(s_pc)
    deallocate(hsrd)
    deallocate(esrd)

! ============================== added by K. Tagami ============= 11.0
    if ( allocated(hsid) ) deallocate(hsid)
    if ( allocated(esid) ) deallocate(esid)
! =============================================================== 11.0

    deallocate(s_0)
    deallocate(s_1)
    deallocate(s_2)
    deallocate(s_tmp)
    deallocate(s_up)
    deallocate(s_dn)
    deallocate(s)
    deallocate(alinvt)

! =========================================== Added by K.Tagami ==== p1100
    deallocate(s_hub)
! ================================================================== p1100
    deallocate(s_vdwdf)

  end subroutine m_Stress_dealloc

  subroutine m_Stress_alloc_grinv_g
    integer :: i
    real(kind=DP) :: ga,gb,gc
    allocate(grinv(ista_kngp:iend_kngp)); grinv = 0.d0
    allocate(g(ista_kngp:iend_kngp,3)); g = 0.d0
    do i = ista_kngp, iend_kngp  !for mpi
      if(gr_l(i).lt.1.0d-20) then
        grinv(i)=0.d0
      else
        grinv(i)=1.d0/gr_l(i)
      endif
    enddo
    do i = ista_kngp, iend_kngp  !for mpi
      ga = ngabc(i,1)
      gb = ngabc(i,2)
      gc = ngabc(i,3)
      g(i,1) = rltv(1,1)*ga+rltv(1,2)*gb+rltv(1,3)*gc
      g(i,2) = rltv(2,1)*ga+rltv(2,2)*gb+rltv(2,3)*gc
      g(i,3) = rltv(3,1)*ga+rltv(3,2)*gb+rltv(3,3)*gc
    enddo
  end subroutine m_Stress_alloc_grinv_g

  subroutine m_Stress_dealloc_grinv_g
    deallocate(grinv)
    deallocate(g)
  end subroutine m_Stress_dealloc_grinv_g

  subroutine m_Stress_ke(nfout,kv3,kp_crdtyp,vkxyz)
    integer, intent(in) :: nfout, kv3, kp_crdtyp
    real(kind=DP), intent(in), dimension(kv3,3,kp_crdtyp) :: vkxyz
    real(kind=DP)       :: fac
    integer             :: ik,i,j,ib,ig, ispin
    real(kind=DP),pointer,dimension(:,:) :: s_ke_mpi   ! MPI
    integer             :: id_sname = -1

    call tstatc0_begin('m_Stress_ke ',id_sname)

    allocate(s_ke_mpi(3,3))
    allocate(g(kg1,3)); g = 0.d0
    allocate(vlength(kg1)); vlength = 0.d0
    allocate(wkc(kg1));  wkc  = 0.d0

    do ispin = ista_spin, iend_spin, af+1
    do ik = ispin, kv3-nspin+ispin, nspin
    !do ik = 1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle               ! MPI
       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            & ,g(1,1),g(1,2),g(1,3),vlength)         ! ->(bottom_Subr.)
       do ig = 1, iba(ik)
          wkc(ig) = 1.d0 &
               &  + 2.d0*ke_a/sqrt(PAI)/ke_sigma  &
               &  * exp(-((vlength(ig)**2/2.d0-ke_e0)/ke_sigma)**2)
       end do
       if(k_symmetry(ik) == GAMMA) then
          do ig = 2, iba(ik)
             wkc(ig) = wkc(ig)*2.d0
          end do
       end if

       do i = 1, 3
       do j = 1, 3
       if(kimg.eq.1) then
          do ib = 1, np_e                           ! MPI
             do ig =1,iba(ik)
                fac = wkc(ig)*g(ig,i)*(g(ig,1)*alinvt(1,j) &
                     &                +g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                s_ke(i,j) = s_ke(i,j) &
                     &        - fac*occup_l(ib,ik)*zaj_l(ig,ib,ik,1)**2
             end do
          end do
       elseif(kimg.eq.2) then
          do ib = 1, np_e                           ! MPI
             do ig =1,iba(ik)
                fac = wkc(ig)*g(ig,i)*(g(ig,1)*alinvt(1,j) &
                     &       +g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                s_ke(i,j) = s_ke(i,j) &
                     &        - fac*occup_l(ib,ik) &
                     &        *(zaj_l(ig,ib,ik,1)**2+zaj_l(ig,ib,ik,2)**2)
             end do
          end do
       endif
       enddo
       enddo
    enddo
    enddo

    call mpi_allreduce(s_ke,s_ke_mpi,3*3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) !MPI
    s_ke = s_ke_mpi  ! MPI

    s_ke=s_ke/dble(kv3)*2.d0

    stress_name = 'KI'
    if(af == 0) then
       s_0 = s_ke
       call check_stress(nfout)
    else
       s_up = s_ke
       call check_af(nfout)
    endif

    deallocate(s_ke_mpi)
    deallocate(g)
    deallocate(vlength)
    deallocate(wkc)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_ke

! ========================= added by K. Tagami ======================= 11.0
  subroutine m_Stress_ke_noncl( nfout, kv3, kp_crdtyp, vkxyz )
    integer, intent(in) :: nfout, kv3, kp_crdtyp
    real(kind=DP), intent(in), dimension(kv3,3,kp_crdtyp) :: vkxyz
    real(kind=DP)       :: fac
    integer             :: ik,i,j,ib,ig
    real(kind=DP),pointer,dimension(:,:) :: s_ke_mpi   ! MPI
    integer             :: id_sname = -1

    integer :: is, k1

    call tstatc0_begin('m_Stress_ke_noncl ',id_sname)

    allocate(s_ke_mpi(3,3))
    allocate(g(kg1,3)); g = 0.d0
    allocate(vlength(kg1)); vlength = 0.d0
    allocate(wkc(kg1));  wkc  = 0.d0

    do ik = 1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle               ! MPI
       call k_plus_G_vectors(ik,kgp,kg1,kv3,iba,nbase,vkxyz,ngabc,rltv&
            & ,g(1,1),g(1,2),g(1,3),vlength)         ! ->(bottom_Subr.)

       do ig = 1, iba(ik)
          wkc(ig) = 1.d0 &
               &  + 2.d0*ke_a/sqrt(PAI)/ke_sigma  &
               &  * exp(-((vlength(ig)**2/2.d0-ke_e0)/ke_sigma)**2)
       end do
       if(k_symmetry(ik) == GAMMA) then
          do ig = 2, iba(ik)
             wkc(ig) = wkc(ig)*2.d0
          end do
       end if

       Do is=1, ndim_spinor
          k1 = ik + is -1

          do i = 1, 3
             do j = 1, 3

                if(kimg.eq.1) then
                   do ib = 1, np_e                           ! MPI
                      do ig =1,iba(ik)
                         fac = wkc(ig)*g(ig,i)*(g(ig,1)*alinvt(1,j) &
                              &                +g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                         s_ke(i,j) = s_ke(i,j) &
                              &        - fac*occup_l(ib,ik)*zaj_l(ig,ib,k1,1)**2
                      end do
                   end do
                elseif(kimg.eq.2) then
                   do ib = 1, np_e                           ! MPI
                      do ig =1,iba(ik)
                         fac = wkc(ig)*g(ig,i)*(g(ig,1)*alinvt(1,j) &
                              &       +g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                         s_ke(i,j) = s_ke(i,j) &
                              &        - fac*occup_l(ib,ik) &
                              &        *(zaj_l(ig,ib,k1,1)**2+zaj_l(ig,ib,k1,2)**2)
                      end do
                   end do
                endif
             enddo
          enddo

       End do

    enddo


    if ( npes > 1 ) then
       call mpi_allreduce( s_ke, s_ke_mpi, 3*3, mpi_double_precision, &
            &              mpi_sum, MPI_CommGroup, ierr ) !MPI
       s_ke = s_ke_mpi  ! MPI
    endif

    s_ke = s_ke/ dble(kv3 /ndim_spinor)

    stress_name = 'KI'
    s_0 = s_ke
    call check_stress(nfout)

    deallocate(s_ke_mpi); deallocate(g);  deallocate(vlength); deallocate(wkc)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_ke_noncl
! ================================================================== 11.0

  subroutine m_Stress_hsr_diff(nfout,kv3)
    integer, intent(in) :: nfout,kv3
!!$ASASASASAS
    !!$integer    :: i1,i2,ia,it,ispin,ik,i,lmt1,lmt2,lmta1,lmta2,iksnl
    integer    :: i1,i2,ia,it,ispin,ik,i,lmt1,lmt2,lmta1,lmta2,iksnl,mapmode
!!$ASASASASAS
    integer    :: id_sname = -1
    real(kind=DP) :: fac, w_n, ek_n, c
    real(kind=DP),pointer,dimension(:,:,:,:,:,:) :: hsrd_mpi, esrd_mpi  ! MPI
    real(kind=DP),pointer,dimension(:,:,:)       :: snld_mpi            ! MPI
!!$ASASASASAS
    real(kind=DP) :: dtmp1=0.d0
    real(kind=DP) :: dtmp2=0.d0
!!$ASASASASAS

    call tstatc0_begin('m_Stress_hsr_diff ',id_sname)

    allocate(hsrd_mpi(natm,nlmt,nlmt,nspin,3,3))   ! MPI
    allocate(esrd_mpi(natm,nlmt,nlmt,nspin,3,3))   ! MPI

    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai(0)
!!$ASASASASAS
    call alloc_zfsincos_mpi()
!!$ASASASASAS
    allocate(fsr_diff_l(np_e,nlmta,ista_k:iend_k)); fsr_diff_l = 0.d0
    allocate(fsi_diff_l(np_e,nlmta,ista_k:iend_k)); fsi_diff_l = 0.d0

    fac = 2.d0/dble(kv3)
    hsrd = 0.d0
    esrd = 0.d0

#ifdef NONLOCAL_DGEMM
    do i1 = 1, 3
    do i2 = 1, 3
!    do ispin = 1, nspin, af+1
    do ispin = ista_spin, iend_spin, af+1
       do ik = ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle         ! MPI
            snld_mpi => snld(:,:,:,i1,i2)   ! MPI
            call m_ES_betar_dot_Psi_4_each_k_snl &
            &   (nfout,snld_mpi,zaj_l,ista_k,iend_k,ik,fsr_diff_l,fsi_diff_l,precalculate_phase=.false.)
            do ia = 1, natm
              it = ityp(ia)
              do i = 1, np_e                          ! MPI
                 w_n = occup_l(i,ik) * fac
                 ek_n = eko_l(i,ik)
                 do lmt1 = 1,ilmt(it)
                   lmta1    = lmta(lmt1,ia)
                   do lmt2 = lmt1,ilmt(it)
                     lmta2    = lmta(lmt2,ia)
                     if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
                       dtmp1 = fsi_l(i,lmta2,ik)
                       dtmp2 = fsi_l(i,lmta1,ik)
                     endif
                     c= w_n * &
                      &     (  fsr_diff_l(i,lmta1,ik)*fsr_l(i,lmta2,ik) &
                      &     +  fsi_diff_l(i,lmta1,ik)*dtmp1 &
                      &     +  fsr_diff_l(i,lmta2,ik)*fsr_l(i,lmta1,ik) &
                      &     +  fsi_diff_l(i,lmta2,ik)*dtmp2 )
                     hsrd(ia,lmt1,lmt2,ispin,i1,i2)=hsrd(ia,lmt1,lmt2,ispin,i1,i2) &
                                   &     + c
                     esrd(ia,lmt1,lmt2,ispin,i1,i2)=esrd(ia,lmt1,lmt2,ispin,i1,i2) &
                                   &     + ek_n * c
                   enddo
                 enddo
              enddo
            enddo
       enddo
    enddo
    enddo
    enddo
#else
    do i1 = 1, 3
    do i2 = i1, 3
      do ia = 1, natm
        it = ityp(ia)
!!$ASASASASAS
!!$        call G_dot_R(natm,ia,pos,kgp,nbmx,ngabc,zfcos,zfsin)
        if ( kv3/nspin == 1 ) then
          mapmode = MAPPED
          call G_dot_R_map(ia,1)
        else
          mapmode = NOTMAPPED
          call G_dot_R_mpi(ia)
        endif
!!$ASASASASAS
!        do ispin = 1, nspin, af+1
        do ispin = ista_spin, iend_spin, af+1
        do ik = ispin, kv3+ispin-nspin, nspin
          if(map_k(ik) /= myrank_k) cycle         ! MPI
          iksnl = (ik-1)/nspin + 1
          snld_mpi => snld(:,:,ista_snl:iend_snl,i1,i2)   ! MPI
          call m_ES_betar_dot_WFs_4_lmta_k &
!!$ASASASASAS
!!$           & (ista_k,iend_k,ik,zaj_l,ia,iksnl,snld_mpi,fsr_diff_l,fsi_diff_l,0)  ! MPI
           & (ista_k,iend_k,ik,zaj_l,ia,iksnl,snld_mpi,fsr_diff_l,fsi_diff_l,mapmode)  ! MPI
!!$ASASASASAS
!!$           & (ista_k,iend_k,ik,zaj_l,ia,snld(1,1,ista_snl,i1,i2),fsr_diff_l,fsi_diff_l)  ! MPI
          do i = 1, np_e                          ! MPI
             w_n = occup_l(i,ik) * fac
             ek_n = eko_l(i,ik)
             do lmt1 = 1,ilmt(it)
               lmta1    = lmta(lmt1,ia)
               do lmt2 = lmt1,ilmt(it)
                 lmta2    = lmta(lmt2,ia)
!!$ASASASASAS
                 if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
                   dtmp1 = fsi_l(i,lmta2,ik)
                   dtmp2 = fsi_l(i,lmta1,ik)
                 endif
!!$                 c= w_n * &
!!$                  &     (  fsr_diff_l(i,lmta1,ik)*fsr_l(i,lmta2,ik) &
!!$                  &     +  fsi_diff_l(i,lmta1,ik)*fsi_l(i,lmta2,ik) &
!!$                  &     +  fsr_diff_l(i,lmta2,ik)*fsr_l(i,lmta1,ik) &
!!$                  &     +  fsi_diff_l(i,lmta2,ik)*fsi_l(i,lmta1,ik) )
                 c= w_n * &
                  &     (  fsr_diff_l(i,lmta1,ik)*fsr_l(i,lmta2,ik) &
                  &     +  fsi_diff_l(i,lmta1,ik)*dtmp1 &
                  &     +  fsr_diff_l(i,lmta2,ik)*fsr_l(i,lmta1,ik) &
                  &     +  fsi_diff_l(i,lmta2,ik)*dtmp2 )
!!$ASASASASAS
                 hsrd(ia,lmt1,lmt2,ispin,i1,i2)=hsrd(ia,lmt1,lmt2,ispin,i1,i2) &
                                  &     + c
                 esrd(ia,lmt1,lmt2,ispin,i1,i2)=esrd(ia,lmt1,lmt2,ispin,i1,i2) &
                                  &     + ek_n * c
               enddo    ! lmt2
             enddo    ! lmt1
          enddo    ! i
       enddo    ! ik
       enddo    ! ispin
    enddo    !ia
    enddo    !i2
    enddo    !i1
#endif

    call mpi_allreduce(hsrd,hsrd_mpi,natm*nlmt*nlmt*nspin*3*3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
    call mpi_allreduce(esrd,esrd_mpi,natm*nlmt*nlmt*nspin*3*3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr) ! MPI
    hsrd = hsrd_mpi  ! MPI
    esrd = esrd_mpi  ! MPI

!!$ASASASASAS
    call dealloc_zfsincos_mpi()
!!$ASASASASAS
    call m_ES_dealloc_arai()
    call m_ES_dealloc_zfsincos()
    deallocate(fsr_diff_l); deallocate(fsi_diff_l)
    deallocate(hsrd_mpi); deallocate(esrd_mpi)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_hsr_diff

! ========================== added by K. Tagami =========== 11.0
  subroutine m_Stress_hsr_diff_noncl(nfout,kv3)
    integer, intent(in) :: nfout,kv3
    integer    :: i1,i2,ia,it,ispin,ik,i,lmt1,lmt2,lmta1,lmta2,iksnl,mapmode

    integer    :: id_sname = -1
    real(kind=DP) :: fac, w_n, ek_n, c1, c2
    real(kind=DP),pointer,dimension(:,:,:,:,:,:) :: hsrd_mpi, esrd_mpi  ! MPI
    real(kind=DP),pointer,dimension(:,:,:)       :: snld_mpi            ! MPI

    real(kind=DP) :: rho_r_magmom(ndim_magmom)
    real(kind=DP) :: rho_i_magmom(ndim_magmom)
    real(kind=DP) :: c_nr( ndim_spinor,ndim_spinor )
    real(kind=DP) :: c_ni( ndim_spinor,ndim_spinor )
    real(kind=DP), allocatable, dimension(:,:,:) :: fsr_diff_l,fsi_diff_l

    integer :: is1, is2, istmp, k1, k2

    call tstatc0_begin('m_Stress_hsr_diff_noncl ',id_sname)

    allocate( hsrd_mpi(natm,nlmt,nlmt,ndim_magmom,3,3) )   ! MPI
    allocate( esrd_mpi(natm,nlmt,nlmt,ndim_magmom,3,3) )   ! MPI

    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai(0)

    call alloc_zfsincos_mpi()

    allocate(fsr_diff_l(np_e,nlmta,ista_k:iend_k)); fsr_diff_l = 0.d0
    allocate(fsi_diff_l(np_e,nlmta,ista_k:iend_k)); fsi_diff_l = 0.d0

    fac = 1.d0 /dble(kv3/ndim_spinor)

    hsrd = 0.d0;  esrd = 0.d0
    hsid = 0.d0;  esid = 0.d0

    do i1 = 1, 3
       do i2 = 1, 3
          do ia = 1, natm
             it = ityp(ia)

             if ( kv3 / ndim_spinor == 1 ) then
                mapmode = MAPPED
                call G_dot_R_map(ia,1)
             else
                mapmode = NOTMAPPED
                call G_dot_R_mpi(ia)
             endif

             Do ik=1, kv3, ndim_spinor
                if (map_k(ik) /= myrank_k) cycle         ! MPI
                iksnl = (ik-1) /ndim_spinor + 1

                snld_mpi => snld(:,:,ista_snl:iend_snl,i1,i2)   ! MPI

                Do is1=1, ndim_spinor
                   k1 = ik + is1 -1
                   call m_ES_betar_dot_WFs_4_lmta_k( ista_k, iend_k, k1, &
                        &                            zaj_l, ia, iksnl, snld_mpi,&
                        &                            fsr_diff_l, fsi_diff_l, mapmode )
                                                                   ! MPI
                End do

                Do i = 1, np_e                          ! MPI
                   w_n = occup_l(i,ik) * fac
                   ek_n = eko_l(i,ik)

                   do lmt1 = 1,ilmt(it)
                      lmta1    = lmta(lmt1,ia)

                      do lmt2 = lmt1,ilmt(it)
!!!!!!!                      do lmt2 = 1,ilmt(it)
                         lmta2    = lmta(lmt2,ia)

                         Do is1=1, ndim_spinor
                            Do is2=1, ndim_spinor
                               k1 = ik + is1 -1
                               k2 = ik + is2 -1
                               istmp = ( is1 -1 )*ndim_spinor + is2

                               c1 = w_n *( fsr_diff_l(i,lmta1,k1)*fsr_l(i,lmta2,k2) &
                                    &     +fsi_diff_l(i,lmta1,k1)*fsi_l(i,lmta2,k2) &
                                    &     +fsr_diff_l(i,lmta2,k2)*fsr_l(i,lmta1,k1) &
                                    &     +fsi_diff_l(i,lmta2,k2)*fsi_l(i,lmta1,k1) )
                               c2 = w_n *(-fsr_diff_l(i,lmta1,k1)*fsi_l(i,lmta2,k2) &
                                    &     +fsi_diff_l(i,lmta1,k1)*fsr_l(i,lmta2,k2) &
                                    &     +fsr_diff_l(i,lmta2,k2)*fsi_l(i,lmta1,k1) &
                                    &     -fsi_diff_l(i,lmta2,k2)*fsr_l(i,lmta1,k1) )

                               c_nr(is1,is2) = c1;      c_ni(is1,is2) = c2
                            End do
                         End do

! -- convert from ss to magmom representaion --
                         rho_r_magmom(1) = c_nr(1,1) + c_nr(2,2);         ! ctot
                         rho_r_magmom(2) = c_nr(1,2) + c_nr(2,1);         ! mx
                         rho_r_magmom(3) = c_ni(2,1) - c_ni(1,2);         ! my
                         rho_r_magmom(4) = c_nr(1,1) - c_nr(2,2);          ! mz
!
                         rho_i_magmom(1) = c_ni(1,1) + c_ni(2,2);         ! ctot
                         rho_i_magmom(2) = c_ni(1,2) + c_ni(2,1);         ! mx
                         rho_i_magmom(3) =-c_nr(2,1) + c_nr(1,2);         ! my
                         rho_i_magmom(4) = c_ni(1,1) - c_ni(2,2);          ! mz
! -----------------------------------------------------------------------

                         Do istmp=1, ndim_magmom
                            hsrd(ia,lmt1,lmt2,istmp,i1,i2) &
                                 & = hsrd(ia,lmt1,lmt2,istmp,i1,i2) &
                                 &  +rho_r_magmom(istmp)
                            esrd(ia,lmt1,lmt2,istmp,i1,i2) &
                                 & = esrd(ia,lmt1,lmt2,istmp,i1,i2) &
                                 &  +ek_n *rho_r_magmom(istmp)

                            hsid(ia,lmt1,lmt2,istmp,i1,i2) &
                                 & = hsid(ia,lmt1,lmt2,istmp,i1,i2) &
                                 &  +rho_i_magmom(istmp)
                            esid(ia,lmt1,lmt2,istmp,i1,i2) &
                                 & = esid(ia,lmt1,lmt2,istmp,i1,i2) &
                                 &  +ek_n *rho_i_magmom(istmp)
                         End do

                      enddo    ! lmt2
                   enddo    ! lmt1
                End do    ! i

             enddo    ! ik
          enddo    !ia
       enddo    !i2
    enddo    !i1

    if ( npes > 1 ) then
       call mpi_allreduce( hsrd, hsrd_mpi, natm*nlmt*nlmt*ndim_magmom*3*3, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr ) ! MPI
       call mpi_allreduce( esrd, esrd_mpi, natm*nlmt*nlmt*ndim_magmom*3*3,&
            &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr ) ! MPI
       hsrd = hsrd_mpi  ! MPI
       esrd = esrd_mpi  ! MPI
!
       call mpi_allreduce( hsid, hsrd_mpi, natm*nlmt*nlmt*ndim_magmom*3*3, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr ) ! MPI
       call mpi_allreduce( esid, esrd_mpi, natm*nlmt*nlmt*ndim_magmom*3*3,&
            &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr ) ! MPI
       hsid = hsrd_mpi  ! MPI
       esid = esrd_mpi  ! MPI

    endif

    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai()
    call m_ES_dealloc_zfsincos()

    deallocate(fsr_diff_l); deallocate(fsi_diff_l)
    deallocate(hsrd_mpi); deallocate(esrd_mpi)

    if ( flg_paw ) then
       Do i1=1, 3
          Do i2=1, 3
             call m_CD_symmtrz_of_ff_noncl_C( hsrd(:,:,:,:,i1,i2), hsid(:,:,:,:,i1,i2) )
             call m_CD_symmtrz_of_ff_noncl_C( esrd(:,:,:,:,i1,i2), esid(:,:,:,:,i1,i2) )
          End Do
       End Do
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_hsr_diff_noncl
! ================================================================ 11.0

!!$  subroutine m_Stress_xcfft_diff(nfout)
!!$    integer, intent(in) :: nfout
!!$    integer    :: id_sname = -1
!!$
!!$    call tstatc0_begin('m_Stress_xcfft_diff ',id_sname)
!!$
!!$    call xc_potential(nfout,Valence_plus_PC_Charge,chgq_l,STRESS &
!!$                 &   ,chgsoft,hsr,hsrd)
!!$
!!$    if(af == 0) then
!!$       s_0 = s_gga1
!!$       stress_name = 'G1'
!!$       call check_stress(nfout)
!!$       s_0 = s_gga2
!!$       stress_name = 'G2'
!!$       call check_stress(nfout)
!!$    else
!!$       s_up = s_gga1
!!$       stress_name = 'G1'
!!$       call check_af(nfout)
!!$       s_up = s_gga2
!!$       stress_name = 'G2'
!!$       call check_af(nfout)
!!$    endif
!!$
!!$    call tstatc0_end(id_sname)
!!$  end subroutine m_Stress_xcfft_diff

  subroutine m_Stress_check_g12(nfout)
    integer, intent(in) :: nfout
    if(af == 0) then
       s_0 = s_gga1
       stress_name = 'G1'
       call check_stress(nfout)
       s_0 = s_gga2
       stress_name = 'G2'
       call check_stress(nfout)
    else
       s_up = s_gga1
       stress_name = 'G1'
       call check_af(nfout)
       s_up = s_gga2
       stress_name = 'G2'
       call check_af(nfout)
    endif
  end subroutine m_Stress_check_g12

  subroutine m_Stress_lclchg_diff(nfout)
    integer, intent(in) :: nfout
    integer    :: n,ilm3,i,it,ik,mdvdb,ia,is,mspin,i1,i2
    integer    :: lmt1,lmt2,il1,il2,tau1,tau2,l3,iiqitg
    real(kind=DP) :: fac,dga,c1,c2
    real(kind=DP),allocatable,target,dimension(:) :: ylm_t
    real(kind=DP), pointer, dimension(:,:)        :: ylmd_mpi
    integer :: ibl1,ibl2
#ifndef _m_Stress_no_loop_exchange_
    integer :: m, maxm, ip, np, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
#endif
    integer    :: id_sname = -1

    if(modnrm == SKIP) return

    call tstatc0_begin('m_Stress_lclchg_diff ',id_sname,level=1)

    do i1=1,3
    do i2=1,3

!    allocate(zfcos_tmp(ista_kngp:iend_kngp))
!    allocate(zfsin_tmp(ista_kngp:iend_kngp))
!    zfcos_tmp = 0.d0; zfsin_tmp = 0.d0
    allocate(vloc(ista_kngp:iend_kngp,kimg)); vloc = 0.d0

    allocate(flchgq(nspin)); flchgq = 0.d0
    allocate(flchgqd(3,3,nspin)); flchgqd = 0.d0

    do ik = 1, kimg
      do it = 1, ntyp
         do i = ista_kngp, iend_kngp  !for mpi
            vloc(i,ik)=vloc(i,ik)+psc_l(i,it)*zfm3_l(i,it,ik)
            if(sw_screening_correction==ON) then
               vloc(i,ik)=vloc(i,ik)-ival(it)*screening%phik(i)/univol*zfm3_l(i,it,ik)
            end if
        enddo
      enddo
    enddo

    call m_PP_find_maximum_l(n)   ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
#ifndef _m_Stress_no_loop_exchange_
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
!    allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
    allocate(ylmd(ista_kngp:iend_kngp,3,n**2)); ylmd = 0.d0
    allocate(ylmd_mpi(ista_kngp:iend_kngp,3)); ylmd_mpi = 0.d0
    do ilm3 = 1, n**2
!mpi       call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm2(1,ilm3))
!!$       call m_pwBS_sphrp2(ilm3,rltv,ylm2_mpi)
!!$       do i = ista_kngp, iend_kngp  !for mpi
!!$          ylm2(i,ilm3) = ylm2_mpi(i)
!!$       enddo
!mpi       call m_pwBS_sphrp2_diff(ilm3,rltv,ylmd(1,1,ilm3))
       call m_pwBS_sphrp2_diff(ilm3,rltv,ylmd_mpi)
       do i = ista_kngp, iend_kngp  !for mpi
          ylmd(i,1,ilm3) = ylmd_mpi(i,1)
          ylmd(i,2,ilm3) = ylmd_mpi(i,2)
          ylmd(i,3,ilm3) = ylmd_mpi(i,3)
       enddo
    end do
    deallocate(ylmd_mpi)
    if(n**2 > nel_Ylm) then
       allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2));  ylm_ext = 0.d0
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)  ! (ilm3,rltv,ngabc,gr_l)->(ylm)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
       deallocate(ylm_t)
    end if

    call blksize(n)
    do ibl1 = ista_kngp, iend_kngp, nblk
    ibl2 = ibl1+nblk-1
    if(ibl2 .gt. iend_kngp) ibl2 = iend_kngp
#ifndef _m_Stress_no_loop_exchange_
    allocate(zfsin_tmp(ibl1:ibl2));zfsin_tmp=0.d0
    allocate(zfcos_tmp(ibl1:ibl2));zfcos_tmp=0.d0
    allocate(mzfsin_tmp(ibl1:ibl2));mzfsin_tmp = 0.d0
    do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle
       call calc_phase2(natm,pos,ia,kgp,ngabc,ibl1,ibl2,zfcos_tmp,zfsin_tmp)
       mzfsin_tmp = -zfsin_tmp
       do iiqitg = nqitg_sp0(it), nqitg_sp(it)
          l3 = iq2l3(iiqitg)
          do m = 1, 2*l3+1
             ilm3 = l3*l3+m
             if(ilm3 <= nel_Ylm) then
                ylm => ylm_l(ibl1:ibl2,ilm3)
             else
                ylm => ylm_ext(ibl1:ibl2,ilm3)
             end if
             do ip = 1, nc(m,iiqitg)
                lmt1 = nc2lmt1(ip,m,iiqitg)
                lmt2 = nc2lmt2(ip,m,iiqitg)
                np = nc2n(ip,m,iiqitg)
                dga = dl2p(lmt1,lmt2,np,it)
                fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
                if(mod(l3,2) == 0) then
                   call even_case
                else
                   call odd_case
                endif
             enddo
          enddo
       enddo
    enddo
    deallocate(mzfsin_tmp)
    deallocate(zfsin_tmp)
    deallocate(zfcos_tmp)
#else
    do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle
       call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos_tmp,zfsin_tmp)
                         ! -(b_Elec.)  -> zfcos_tmp, zfsin_tmp
       do lmt1 = 1,ilmt(it)
          il1 = ltp(lmt1,it); tau1 = taup(lmt1,it)
          do lmt2 = lmt1, ilmt(it)
            il2 = ltp(lmt2,it); tau2 = taup(lmt2,it)
            fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
            do n = 1, il2p(lmt1,lmt2,it)
              ilm3 = isph(lmt1,lmt2,n,it)
              l3   =  il3(ilm3)
              iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
              if(iiqitg == 0) cycle
              if(ilm3 <= nel_Ylm) then
                 ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
              else
                 call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
                 ylm => ylm_t(ista_kngp:iend_kngp)
              end if
              dga = dl2p(lmt1,lmt2,n,it)

              if (mod(il1+il2,2) == 0) then
                 call even_case
              else
                 call odd_case
              endif

            enddo
          enddo
       enddo
    enddo
#endif
    enddo

!    deallocate(zfsin_tmp); deallocate(zfcos_tmp)
    deallocate(il3)
!    deallocate(ylm_t)
    deallocate(ylmd)
    if(allocated(ylm_ext)) deallocate(ylm_ext)
#ifndef _m_Stress_no_loop_exchange_
    deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif


    deallocate(vloc)
    deallocate(flchgq); deallocate(flchgqd)

    enddo
    enddo

    call mpi_allreduce(mpi_in_place,s_xc1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(mpi_in_place,s_h1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(mpi_in_place,s_loc1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

    s_xc1  = s_xc1  * univol
    s_h1   = s_h1   * univol * PAI4
    s_loc1 = s_loc1 * univol

    if(af == 0) then
       s_0 = s_xc1
       stress_name = 'X1'
       call check_stress(nfout)
       s_0 = s_h1
       stress_name = 'H1'
       call check_stress(nfout)
       s_0 = s_loc1
       stress_name = 'L1'
       call check_stress(nfout)
    else
       s_up = s_xc1
       stress_name = 'X1'
       call check_af(nfout)
       s_up = s_h1
       stress_name = 'H1'
       call check_af(nfout)
       s_up = s_loc1
       stress_name = 'L1'
       call check_af(nfout)
    endif

    call tstatc0_end(id_sname)

  contains

    subroutine blksize(n)
      integer, intent(in) :: n
      integer :: ncache
      ncache = (m_CtrlP_cachesize()*1024)*3/4
      if(ncache == 0) then
         nblk = iend_kngp-ista_kngp+1
      else
         nblk=ncache/(8*(n**2))
      end if
      if (nblk<32) nblk = 1000
      !write(nfout,'(a,i8)') ' !** block size : ',nblk
    end subroutine blksize

    subroutine even_case
      do is=1,nspin,af+1
         flchgq(is) = fac*real(zi**(-l3))*dga &
              & *hsr(ia,lmt1,lmt2,is)
         flchgqd(i1,i2,is) = fac*real(zi**(-l3))*dga &
         & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if(kimg == 1) then
         call real_case(ilm3,iiqitg,zfcos_tmp)
      else
         call complex_case(ilm3,iiqitg,zfcos_tmp,mzfsin_tmp)
      endif
    end subroutine even_case

    subroutine odd_case
      do is=1,nspin,af+1
         flchgq(is) = fac*aimag(zi**(-l3))*dga &
              & *hsr(ia,lmt1,lmt2,is)
         flchgqd(i1,i2,is) = fac*aimag(zi**(-l3))*dga &
         & *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if(kimg == 1) then
         call real_case(ilm3,iiqitg,zfsin_tmp)
      else
         call complex_case(ilm3,iiqitg,zfsin_tmp,zfcos_tmp)
      endif
    end subroutine odd_case

    subroutine real_case(ilm3,iiqitg,zf)
      integer,       intent(in) :: ilm3,iiqitg
      real(kind=DP), intent(in), dimension(ibl1:ibl2) :: zf
      real(kind=DP) :: sxc1,sh1,sloc1,c1,c2,c3up,c3down
      integer :: iy
      integer :: id_sname = -1

      sxc1 = 0.d0;sh1=0.d0;sloc1=0.d0
      do is=1,nspin,af+1
        flchgq(is) = flchgq(is) * dble(iwei(ia))
        flchgqd(i1,i2,is) = flchgqd(i1,i2,is) * dble(iwei(ia))
      enddo

      if(nspin==1 .or. af==1) then
         mspin=(nspin+1)/nspin
         do i = ibl1, ibl2  !for mpi
            iy = i - ibl1+1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)

            sxc1=sxc1 &
                 &  +c3up*vxc_l(i,1,1)*zf(i)
            sh1=sh1 &
                 &  +grinv(i)**2*c3up &
                 &  *(chgq_l(i,1,1)+chgq_l(i,1,nspin))/mspin &
                 &  *zf(i)
            sloc1=sloc1 &
                 &  +c3up*zf(i)*vloc(i,1)
         enddo
!xocl end spread sum(s_xc1),sum(s_loc1),sum(s_h1)
      else  ! nspin ==2 and af /= 1
         do i = ibl1, ibl2  !for mpi
            iy = i - ibl1 + 1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            c3down = &
             &  - flchgq(nspin) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,nspin)*ylm(iy)*qitg_l(i,iiqitg)

            sxc1=sxc1 &
                 &  +(c3up*vxc_l(i,1,1) &
                 &   +c3down*vxc_l(i,1,nspin)) &
                 &   *zf(i)
            sh1=sh1 &
                 &  +grinv(i)**2*(c3up+c3down) &
                 &  *(chgq_l(i,1,1)+chgq_l(i,1,nspin))*zf(i)
            sloc1=sloc1 &
                 &  +(c3up+c3down) &
                 &  *zf(i)*vloc(i,1)
         enddo
!xocl end spread sum(s_xc1),sum(s_loc1),sum(s_h1)
      endif
      s_xc1(i1,i2)  =  s_xc1(i1,i2)  + sxc1
      s_h1(i1,i2)   =  s_h1(i1,i2)   + sh1
      s_loc1(i1,i2) =  s_loc1(i1,i2) + sloc1
    end subroutine real_case

    subroutine complex_case(ilm3,iiqitg,zf1,zf2)
      integer,       intent(in) :: ilm3,iiqitg
      real(kind=DP), intent(in), dimension(ibl1:ibl2) :: zf1,zf2
      real(kind=DP) :: sxc1,sh1,sloc1,c1,c2,c3up,c3down
      integer :: iy
      real(kind=DP) :: factor
      integer :: id_sname = -1
      sxc1 = 0.d0;sh1=0.d0;sloc1=0.d0
      if(nspin.eq.1 .or. af.eq.1) then
!         mspin=(nspin+1)/nspin
!         do i = ista_kngp, iend_kngp  !for mpi
         factor = 0.5d0
         if(af .eq. 1) factor = 1.d0
         do i = ibl1, ibl2  !for mpi
!            iy = i - ista_kngp + 1
            iy = i - ibl1 + 1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up = &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)

            sxc1=sxc1 &
                 &  +c3up &
                 &  *(vxc_l(i,1,1)*zf1(i)+vxc_l(i,2,1)*zf2(i))
            sh1=sh1 &
                 &  +grinv(i)**2*c3up &
!                 &  *((chgq_l(i,1,1)+chgq_l(i,1,nspin))/mspin*zf1(i) &
!                 &   +(chgq_l(i,2,1)+chgq_l(i,2,nspin))/mspin*zf2(i))
                 &  *factor*((chgq_l(i,1,1)+chgq_l(i,1,nspin))*zf1(i) &
                 &   +(chgq_l(i,2,1)+chgq_l(i,2,nspin))*zf2(i))
            sloc1=sloc1 &
                 &  +c3up &
                 &  *(zf1(i)*vloc(i,1)+zf2(i)*vloc(i,2))
         enddo
!xocl end spread sum(s_xc1),sum(s_loc1),sum(s_h1)
      else  ! nspin ==2 and af /= 1
!         do i = ista_kngp, iend_kngp  !for mpi
!            iy = i - ista_kngp + 1
         do i = ibl1, ibl2  !for mpi
            iy = i - ibl1 + 1
            c1=(  g(i,1)*alinvt(1,i2) &
             &  + g(i,2)*alinvt(2,i2) &
             &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
            c2=(  ylmd(i,1,ilm3)*alinvt(1,i2) &
             &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
             &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)
            c3up =  &
             &  - flchgq(1) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,1)*ylm(iy)*qitg_l(i,iiqitg)
            c3down = &
             &  - flchgq(nspin) &
             &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
             &  +qitg_diff_l(i,iiqitg)*c1) &
             &  +c2*qitg_l(i,iiqitg)) &
             &  + flchgqd(i1,i2,nspin)*ylm(iy)*qitg_l(i,iiqitg)

            sxc1=sxc1 &
                 &  +c3up &
                 &  *(vxc_l(i,1,1)*zf1(i)+vxc_l(i,2,1)*zf2(i)) &
                 &  +c3down &
                 &  *(vxc_l(i,1,nspin)*zf1(i)+vxc_l(i,2,nspin)*zf2(i))
            sh1=sh1 &
                 &  +grinv(i)**2*(c3up+c3down) &
                 &  *((chgq_l(i,1,1)+chgq_l(i,1,nspin))*zf1(i) &
                 &   +(chgq_l(i,2,1)+chgq_l(i,2,nspin))*zf2(i))
            sloc1=sloc1 &
                 &  +(c3up+c3down) &
                 &  *(zf1(i)*vloc(i,1)+zf2(i)*vloc(i,2))
         enddo
!xocl end spread sum(s_xc1),sum(s_loc1),sum(s_h1)
      endif
      s_xc1(i1,i2)  =  s_xc1(i1,i2)  + sxc1
      s_h1(i1,i2)   =  s_h1(i1,i2)   + sh1
      s_loc1(i1,i2) =  s_loc1(i1,i2) + sloc1
    end subroutine complex_case

  end subroutine m_Stress_lclchg_diff

! ========================= added by K. Tagami ======================= 11.0
  subroutine m_Stress_lclchg_diff_noncl(nfout)

    integer, intent(in) :: nfout
    integer    :: n,ilm3,i,it,ik,mdvdb,ia,is,mspin,i1,i2
    integer    :: lmt1,lmt2,il1,il2,tau1,tau2,l3,iiqitg
    real(kind=DP) :: fac,dga,c1,c2
    real(kind=DP),allocatable,target,dimension(:) :: ylm_t
    real(kind=DP), pointer, dimension(:,:)        :: ylmd_mpi
    integer :: ibl1,ibl2
#ifndef _m_Stress_no_loop_exchange_
    integer :: m, maxm, ip, np, iq
    integer, parameter :: mcritical = 4*2+1
    integer, allocatable, dimension(:) :: nqitg_sp, nqitg_sp0 !d(ntyp)
    integer, allocatable, dimension(:) :: iq2l3 ! d(nqitg)
    integer, allocatable, dimension(:,:) :: nc  ! d(maxm,nqitg)
    integer :: mc ! maxval(nc)
    integer, allocatable, dimension(:,:,:) :: nc2lmt1, nc2lmt2, nc2n ! d(mc,maxm,nqitg)
#endif
    integer    :: id_sname = -1

    if(modnrm == SKIP) return

    call tstatc0_begin('m_Stress_lclchg_diff_noncol ',id_sname,level=1)

    do i1=1,3
    do i2=1,3

    allocate(vloc(ista_kngp:iend_kngp,kimg)); vloc = 0.d0

    allocate(flchgq(ndim_magmom)); flchgq = 0.d0
    allocate(flchgqd(3,3,ndim_magmom)); flchgqd = 0.d0

    do ik = 1, kimg
      do it = 1, ntyp
         do i = ista_kngp, iend_kngp  !for mpi
            vloc(i,ik)=vloc(i,ik)+psc_l(i,it)*zfm3_l(i,it,ik)
            if(sw_screening_correction==ON) then
               vloc(i,ik)=vloc(i,ik)-ival(it)*screening%phik(i)/univol*zfm3_l(i,it,ik)
            end if
        enddo
      enddo
    enddo


    call m_PP_find_maximum_l(n)   ! n-1: maximum l
    n = (n-1) + (n-1) + 1
    allocate(il3(n**2)); call substitute_il3(n**2,il3) ! -(b_Elec..)
#ifndef _m_Stress_no_loop_exchange_
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
!    allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
    allocate(ylmd(ista_kngp:iend_kngp,3,n**2)); ylmd = 0.d0
    allocate(ylmd_mpi(ista_kngp:iend_kngp,3)); ylmd_mpi = 0.d0
    do ilm3 = 1, n**2
!mpi       call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm2(1,ilm3))
!!$       call m_pwBS_sphrp2(ilm3,rltv,ylm2_mpi)
!!$       do i = ista_kngp, iend_kngp  !for mpi
!!$          ylm2(i,ilm3) = ylm2_mpi(i)
!!$       enddo
!mpi       call m_pwBS_sphrp2_diff(ilm3,rltv,ylmd(1,1,ilm3))
       call m_pwBS_sphrp2_diff(ilm3,rltv,ylmd_mpi)
       do i = ista_kngp, iend_kngp  !for mpi
          ylmd(i,1,ilm3) = ylmd_mpi(i,1)
          ylmd(i,2,ilm3) = ylmd_mpi(i,2)
          ylmd(i,3,ilm3) = ylmd_mpi(i,3)
       enddo
    end do
    deallocate(ylmd_mpi)
    if(n**2 > nel_Ylm) then
       allocate(ylm_ext(ista_kngp:iend_kngp,nel_Ylm+1:n**2));  ylm_ext = 0.d0
       allocate(ylm_t(ista_kngp:iend_kngp)); ylm_t = 0.d0
       do ilm3 = nel_Ylm+1, n**2
          call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)  ! (ilm3,rltv,ngabc,gr_l)->(ylm)
          ylm_ext(:,ilm3) = ylm_t(:)
       end do
       deallocate(ylm_t)
    end if

    call blksize(n)
    do ibl1 = ista_kngp, iend_kngp, nblk
    ibl2 = ibl1+nblk-1
    if(ibl2 .gt. iend_kngp) ibl2 = iend_kngp
#ifndef _m_Stress_no_loop_exchange_
    allocate(zfsin_tmp(ibl1:ibl2));zfsin_tmp=0.d0
    allocate(zfcos_tmp(ibl1:ibl2));zfcos_tmp=0.d0
    allocate(mzfsin_tmp(ibl1:ibl2));mzfsin_tmp = 0.d0
    do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle
       call calc_phase2(natm,pos,ia,kgp,ngabc,ibl1,ibl2,zfcos_tmp,zfsin_tmp)
       mzfsin_tmp = -zfsin_tmp
       do iiqitg = nqitg_sp0(it), nqitg_sp(it)
          l3 = iq2l3(iiqitg)
          do m = 1, 2*l3+1
             ilm3 = l3*l3+m
             if(ilm3 <= nel_Ylm) then
                ylm => ylm_l(ibl1:ibl2,ilm3)
             else
                ylm => ylm_ext(ibl1:ibl2,ilm3)
             end if
             do ip = 1, nc(m,iiqitg)
                lmt1 = nc2lmt1(ip,m,iiqitg)
                lmt2 = nc2lmt2(ip,m,iiqitg)
                np = nc2n(ip,m,iiqitg)
                dga = dl2p(lmt1,lmt2,np,it)
                fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
                if(mod(l3,2) == 0) then
                   call even_case_noncl
                else
                   call odd_case_noncl
                endif
             enddo
          enddo
       enddo
    enddo
    deallocate(mzfsin_tmp)
    deallocate(zfsin_tmp)
    deallocate(zfcos_tmp)
#else
    do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(mdvdb == SKIP) cycle
       call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos_tmp,zfsin_tmp)
                         ! -(b_Elec.)  -> zfcos_tmp, zfsin_tmp
       do lmt1 = 1,ilmt(it)
          il1 = ltp(lmt1,it); tau1 = taup(lmt1,it)
          do lmt2 = lmt1, ilmt(it)
            il2 = ltp(lmt2,it); tau2 = taup(lmt2,it)
            fac = 2.d0 ; if(lmt1 == lmt2) fac = 1.d0
            do n = 1, il2p(lmt1,lmt2,it)
              ilm3 = isph(lmt1,lmt2,n,it)
              l3   =  il3(ilm3)
              iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
              if(iiqitg == 0) cycle
              if(ilm3 <= nel_Ylm) then
                 ylm => ylm_l(ista_kngp:iend_kngp,ilm3)
              else
                 call m_pwBS_sphrp2(ilm3,rltv,ista_kngp,iend_kngp,ylm_t)
                 ylm => ylm_t(ista_kngp:iend_kngp)
              end if
              dga = dl2p(lmt1,lmt2,n,it)

              if (mod(il1+il2,2) == 0) then
                 call even_case_noncl
              else
                 call odd_case_noncl
              endif

            enddo
          enddo
       enddo
    enddo
#endif
    enddo

!    deallocate(zfsin_tmp); deallocate(zfcos_tmp)
    deallocate(il3)
!    deallocate(ylm_t)
    deallocate(ylmd)
    if(allocated(ylm_ext)) deallocate(ylm_ext)
#ifndef _m_Stress_no_loop_exchange_
    deallocate(nc2n,nc2lmt2,nc2lmt1,nc,iq2l3,nqitg_sp,nqitg_sp0)
#endif


    deallocate(vloc)
    deallocate(flchgq); deallocate(flchgqd)

    enddo
    enddo

    call mpi_allreduce(mpi_in_place,s_xc1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(mpi_in_place,s_h1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(mpi_in_place,s_loc1,9 &
                 &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

    s_xc1  = s_xc1  * univol
    s_h1   = s_h1   * univol * PAI4
    s_loc1 = s_loc1 * univol

    if(af == 0) then
       s_0 = s_xc1
       stress_name = 'X1'
       call check_stress(nfout)
       s_0 = s_h1
       stress_name = 'H1'
       call check_stress(nfout)
       s_0 = s_loc1
       stress_name = 'L1'
       call check_stress(nfout)
    else
       s_up = s_xc1
       stress_name = 'X1'
       call check_af(nfout)
       s_up = s_h1
       stress_name = 'H1'
       call check_af(nfout)
       s_up = s_loc1
       stress_name = 'L1'
       call check_af(nfout)
    endif

    call tstatc0_end(id_sname)

  contains

    subroutine blksize(n)
      integer, intent(in) :: n
      integer :: ncache
      ncache = (m_CtrlP_cachesize()*1024)*3/4
      if(ncache == 0) then
         nblk = iend_kngp-ista_kngp+1
      else
         nblk=ncache/(8*(n**2))
      end if
      if (nblk<32) nblk = 1000
    end subroutine blksize

    subroutine even_case_noncl
      integer :: is

      do is = 1,ndim_magmom
         flchgq(is) = fac *real(zi**(-l3)) *dga &
              &           *hsr(ia,lmt1,lmt2,is)
         flchgqd(i1,i2,is) = fac *real(zi**(-l3)) *dga &
                    &                  *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if (kimg == 1) then
         call real_case_noncl(ilm3,iiqitg,zfcos_tmp)
      else
         call complex_case_noncl(ilm3,iiqitg,zfcos_tmp,mzfsin_tmp)
      endif
    end subroutine even_case_noncl

    subroutine odd_case_noncl
      integer :: is

      do is = 1,ndim_magmom
         flchgq(is) = fac *aimag(zi**(-l3)) *dga &
              &           *hsr(ia,lmt1,lmt2,is)
         flchgqd(i1,i2,is) = fac * aimag(zi**(-l3)) *dga &
              &                  *hsrd(ia,lmt1,lmt2,is,i1,i2)
      enddo
      if (kimg == 1) then
         call real_case_noncl(ilm3,iiqitg,zfsin_tmp)
      else
         call complex_case_noncl(ilm3,iiqitg,zfsin_tmp,zfcos_tmp)
      endif
    end subroutine odd_case_noncl

    subroutine real_case_noncl( ilm3,iiqitg, zf )
      integer,       intent(in) :: ilm3,iiqitg
      real(kind=DP), intent(in), dimension(ibl1:ibl2) :: zf
      real(kind=DP) :: sxc1,sh1,sloc1,c1,c2
      real(kind=DP), allocatable, dimension(:) :: drhodh
      integer :: iy

      integer :: is
      allocate(drhodh(ndim_magmom))
      do is=1,ndim_magmom
         flchgq(is) = flchgq(is) * dble(iwei(ia))
         do i1=1,3
            do i2=1,3
               flchgqd(i1,i2,is) = flchgqd(i1,i2,is) * dble(iwei(ia))
            enddo
         enddo
      enddo

      sxc1=0.d0;sh1=0.d0;sloc1=0.d0
      do i = ibl1, ibl2  !for mpi
         iy = i - ibl1+1
         c1 =  ( g(i,1)*alinvt(1,i2) &
              &  + g(i,2)*alinvt(2,i2) &
              &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
         c2 =  ( ylmd(i,1,ilm3)*alinvt(1,i2) &
              &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
              &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)

         Do is=1, ndim_magmom
            drhodh(is) =  &
                 &  - flchgq(is) &
                 &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
                 &  +qitg_diff_l(i,iiqitg)*c1) &
                 &  +c2*qitg_l(i,iiqitg)) &
                 &  + flchgqd(i1,i2,is)*ylm(iy)*qitg_l(i,iiqitg)
         End do
         Do is=1, ndim_magmom
            sxc1 = sxc1 &
                 &      + drhodh(is)*vxc_l(i,1,is)*zf(i)
         End do

!               s_xc1_mpi1(i1,i2) = s_xc1_mpi1(i1,i2) &
!                    &            + drhodh(1)*vxc_l(i,1,1)*zf(i)
!
         sh1 = sh1 &
               & + grinv(i)**2 *drhodh(1) *chgq_l(i,1,1) *zf(i)
         sloc1 = sloc1 &
               & + drhodh(1)*zf(i)*vloc(i,1)
      end do

      s_xc1(i1,i2)  =  s_xc1(i1,i2)  + sxc1
      s_h1(i1,i2)   =  s_h1(i1,i2)   + sh1
      s_loc1(i1,i2) =  s_loc1(i1,i2) + sloc1
      deallocate(drhodh)

    end subroutine real_case_noncl

    subroutine complex_case_noncl(ilm3,iiqitg,zf1,zf2)
      integer,       intent(in) :: ilm3,iiqitg
      real(kind=DP), intent(in), dimension(ibl1:ibl2) &
           &  :: zf1,zf2
      integer :: iy
      integer :: is
      real(kind=DP) :: sxc1,sh1,sloc1,c1,c2
      real(kind=DP), allocatable, dimension(:) :: drhodh

      allocate(drhodh(ndim_magmom))
      sxc1=0.d0;sh1=0.d0;sloc1=0.d0
      do i = ibl1, ibl2  !for mpi
         iy = i - ibl1 + 1
         c1= ( g(i,1)*alinvt(1,i2) &
              &  + g(i,2)*alinvt(2,i2) &
              &  + g(i,3)*alinvt(3,i2) )  *g(i,i1) * grinv(i)
         c2= ( ylmd(i,1,ilm3)*alinvt(1,i2) &
              &  + ylmd(i,2,ilm3)*alinvt(2,i2) &
              &  + ylmd(i,3,ilm3)*alinvt(3,i2) )  *g(i,i1)

         Do is=1, ndim_magmom
            drhodh(is) =  &
              &  - flchgq(is) &
              &  *(ylm(iy)*(qitg_l(i,iiqitg)*alinvt(i1,i2) &
              &  +qitg_diff_l(i,iiqitg)*c1) &
              &  +c2*qitg_l(i,iiqitg)) &
              &  + flchgqd(i1,i2,is)*ylm(iy)*qitg_l(i,iiqitg)
         End do
         Do is=1, ndim_magmom
            sxc1 = sxc1 &
                 &    + drhodh(is) &
                 &    *( vxc_l(i,1,is) *zf1(i) &
                 &    +vxc_l(i,2,is) *zf2(i) )
         End do

         sh1 = sh1 &
              &    + grinv(i)**2 *drhodh(1) &
              &    *( chgq_l(i,1,1)*zf1(i) &
              &    +chgq_l(i,2,1)*zf2(i) )
         sloc1 = sloc1 &
              &    + drhodh(1) &
              &    *(zf1(i)*vloc(i,1)+zf2(i)*vloc(i,2))
      enddo
      s_xc1(i1,i2)  =  s_xc1(i1,i2)  + sxc1
      s_h1(i1,i2)   =  s_h1(i1,i2)   + sh1
      s_loc1(i1,i2) =  s_loc1(i1,i2) + sloc1
      deallocate(drhodh)
    end subroutine complex_case_noncl

  end subroutine m_Stress_lclchg_diff_noncl
! =========================================================================== 11.0

  subroutine m_Stress_xc(nfout)
    integer, intent(in) :: nfout
    integer    :: is,ig,i,j,ia,it,ipc
    real(kind=DP) :: vxcp,pc
    integer    :: id_sname = -1
    integer       :: iend  !mpi
    real(kind=DP) :: vxcp_mpi
    real(kind=DP), pointer, dimension(:,:) :: s_xc2_mpi1, s_xc2_mpi2
    real(kind=DP), pointer, dimension(:,:) :: s_pc_mpi1, s_pc_mpi2

    call tstatc0_begin('m_Stress_xc ',id_sname)

    allocate(zfcos_tmp(ista_kngp:iend_kngp))
    allocate(zfsin_tmp(ista_kngp:iend_kngp))
    zfcos_tmp = 0.d0; zfsin_tmp = 0.d0

    vxcp_mpi = 0.d0
    vxcp=0.d0
    if(kimg==1) then
       do is=1,nspin
          iend = iend_kngp
          if( iend_kngp > kg ) iend = kg
          if( ista_kngp <= iend ) then
             do i = ista_kngp, iend  !for mpi
                vxcp_mpi=vxcp_mpi+vxc_l(i,1,is)*chgsoft(i,1,is)
             enddo
          endif
!xocl end spread sum(vxcp)
       enddo
    else
       do is=1,nspin
          iend = iend_kngp
          if( iend_kngp > kg ) iend = kg
          if( ista_kngp <= iend ) then
             do i = ista_kngp, iend  !for mpi
                vxcp_mpi=vxcp_mpi+vxc_l(i,1,is)*chgsoft(i,1,is) &
                     &    +vxc_l(i,2,is)*chgsoft(i,2,is)
             enddo
          endif
!xocl end spread sum(vxcp)
       enddo
    endif
    call mpi_allreduce(vxcp_mpi,vxcp,1 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

#ifndef DISABLE_VDWDF
    s_xc2 = ( exc +ecnl_vdwdf -vxcp *univol ) *alinvt
#else
    s_xc2 = ( exc -vxcp *univol ) *alinvt
#endif

    allocate(s_xc2_mpi1(3,3)); s_xc2_mpi1 = 0.d0
    allocate(s_xc2_mpi2(3,3)); s_xc2_mpi2 = 0.d0
    allocate(s_pc_mpi1(3,3)); s_pc_mpi1 = 0.d0
    allocate(s_pc_mpi2(3,3)); s_pc_mpi2 = 0.d0

    do ia = 1, natm
    it=ityp(ia)
    if(itpcc(it).ne.0) then
       call calc_phase2(natm,pos,ia,kgp,ngabc,ista_kngp,iend_kngp,zfcos_tmp,zfsin_tmp)
                         ! -(b_Elec.)  -> zfcos_tmp, zfsin_tmp
       ipc=itpcc(it)
       if(kimg==1) then
          do is=1,nspin
          do i=1,3
          do j=1,3
             do ig = ista_kngp, iend_kngp  !for mpi
                pc = (rhpcg_l(ig,ipc)*alinvt(i,j) &
                     &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i) &
                     &   *(g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j) &
                     &    +g(ig,3)*alinvt(3,j))) &
                     &   *univol &
                     &   * dble(iwei(ia))
                s_xc2_mpi1(i,j) = s_xc2_mpi1(i,j) &
                     &   - vxc_l(ig,1,is) * zfcos_tmp(ig) * pc / nspin
                s_pc_mpi1(i,j) = s_pc_mpi1(i,j) &
                     &   - vxcpc_l(ig,1) * zfcos_tmp(ig) * pc / nspin
             enddo
          enddo
          enddo
          enddo
!xocl end spread sum(s_xc2),sum(s_pc)
       else
          do is=1,nspin
          do i=1,3
          do j=1,3
             do ig = ista_kngp, iend_kngp  !for mpi
                pc = (rhpcg_l(ig,ipc)*alinvt(i,j) &
                     &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i) &
                     &   *(g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j) &
                     &    +g(ig,3)*alinvt(3,j))) &
                     &   *univol &
                     &   * dble(iwei(ia))
                s_xc2_mpi1(i,j) = s_xc2_mpi1(i,j) &
                     &   - (vxc_l(ig,1,is) * zfcos_tmp(ig) &
                     &    - vxc_l(ig,2,is) * zfsin_tmp(ig) ) * pc / nspin
                s_pc_mpi1(i,j) = s_pc_mpi1(i,j) &
                     &   - (vxcpc_l(ig,1) * zfcos_tmp(ig) &
                     &    - vxcpc_l(ig,2) * zfsin_tmp(ig) ) * pc / nspin
             enddo
          enddo
          enddo
          enddo
!xocl end spread sum(s_xc2),sum(s_pc)
       endif
    endif
    enddo

    call mpi_allreduce(s_xc2_mpi1,s_xc2_mpi2,9 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(s_pc_mpi1,s_pc_mpi2,9 &
         &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)

    s_xc2 = s_xc2 + s_xc2_mpi2
    s_pc = s_pc + s_pc_mpi2 + epc * alinvt
    deallocate(s_xc2_mpi1); deallocate(s_xc2_mpi2)
    deallocate(s_pc_mpi1); deallocate(s_pc_mpi2)

    s_0 = s_pc
    stress_name = 'PC'
    call check_stress(nfout)
    s_0 = s_xc1
    stress_name = 'X1'
    call check_stress(nfout)
    s_0 = s_xc2
    stress_name = 'X2'
    call check_stress(nfout)
    if(af == 0) then
       s_0 = s_xc1 + s_xc2 + s_gga1 + s_gga2
       stress_name = 'XC'
       call check_stress(nfout)
    else
       s_0 = s_xc2
       stress_name = 'X2'
       call check_stress(nfout)
    endif

    deallocate(zfsin_tmp); deallocate(zfcos_tmp)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_xc

! ==================================== added by K. Tagami ============= 11.0
  subroutine m_Stress_xc_noncl(nfout)
    integer, intent(in) :: nfout
    integer    :: is,ig,i,j,ia,it,ipc
    real(kind=DP) :: vxcp,pc
    integer    :: id_sname = -1
    integer       :: iend  !mpi
    real(kind=DP) :: vxcp_mpi
    real(kind=DP), pointer, dimension(:,:) :: s_xc2_mpi1, s_xc2_mpi2
    real(kind=DP), pointer, dimension(:,:) :: s_pc_mpi1, s_pc_mpi2

    call tstatc0_begin('m_Stress_xc_noncl ',id_sname)

    allocate(zfcos_tmp(ista_kngp:iend_kngp))
    allocate(zfsin_tmp(ista_kngp:iend_kngp))
    zfcos_tmp = 0.d0; zfsin_tmp = 0.d0

    vxcp_mpi = 0.d0
    vxcp=0.d0

    if (kimg==1) return

    do is=1,ndim_magmom
       iend = iend_kngp
       if( iend_kngp > kg ) iend = kg
       if( ista_kngp <= iend ) then
          do i = ista_kngp, iend  !for mpi
             vxcp_mpi = vxcp_mpi +vxc_l(i,1,is)*chgsoft(i,1,is) &
                  &              +vxc_l(i,2,is)*chgsoft(i,2,is)
          enddo
       endif
    enddo

    call mpi_allreduce( vxcp_mpi, vxcp, 1, &
         &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )

#ifndef DISABLE_VDWDF
    s_xc2 = ( exc +ecnl_vdwdf -vxcp *univol ) *alinvt
#else
    s_xc2 = ( exc -vxcp *univol ) *alinvt
#endif

    allocate(s_xc2_mpi1(3,3)); s_xc2_mpi1 = 0.d0
    allocate(s_xc2_mpi2(3,3)); s_xc2_mpi2 = 0.d0
    allocate(s_pc_mpi1(3,3)); s_pc_mpi1 = 0.d0
    allocate(s_pc_mpi2(3,3)); s_pc_mpi2 = 0.d0

    do ia = 1, natm
       it = ityp(ia)

       if (itpcc(it).eq.0 ) cycle
       call calc_phase2( natm, pos, ia, kgp, ngabc, ista_kngp, iend_kngp, &
            &            zfcos_tmp,zfsin_tmp )
                         ! -(b_Elec.)  -> zfcos_tmp, zfsin_tmp
       ipc = itpcc(it)
       do i=1,3
          do j=1,3
             do ig = ista_kngp, iend_kngp  !for mpi
                pc = (rhpcg_l(ig,ipc)*alinvt(i,j) &
                     &   +rhpcg_diff_l(ig,ipc)*grinv(ig)*g(ig,i) &
                     &   *(g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j) &
                     &    +g(ig,3)*alinvt(3,j))) &
                     &   *univol &
                     &   * dble(iwei(ia))

                s_xc2_mpi1(i,j) = s_xc2_mpi1(i,j) &
                     &    -(  vxc_l(ig,1,1) * zfcos_tmp(ig) &
                     &      - vxc_l(ig,2,1) * zfsin_tmp(ig) ) * pc

                s_pc_mpi1(i,j) = s_pc_mpi1(i,j) &
                     &     -(  vxcpc_l(ig,1) * zfcos_tmp(ig) &
                     &       - vxcpc_l(ig,2) * zfsin_tmp(ig) ) * pc
             enddo
          enddo
       enddo
    enddo

    call mpi_allreduce( s_xc2_mpi1, s_xc2_mpi2, 9, &
         &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
    call mpi_allreduce( s_pc_mpi1, s_pc_mpi2, 9, &
         &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )

    s_xc2 = s_xc2 + s_xc2_mpi2
    s_pc = s_pc + s_pc_mpi2 + epc * alinvt
    deallocate(s_xc2_mpi1); deallocate(s_xc2_mpi2)
    deallocate(s_pc_mpi1); deallocate(s_pc_mpi2)


! ---
    s_0 = s_pc
    stress_name = 'PC'
    call check_stress(nfout)

    s_0 = s_xc1 + s_xc2 + s_gga1 + s_gga2;
    stress_name = 'XC'
    call check_stress(nfout)

    deallocate(zfsin_tmp); deallocate(zfcos_tmp)

    call tstatc0_end(id_sname)

  end subroutine m_Stress_xc_noncl
! ============================================================= 11.0

  subroutine m_Stress_loc(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,it,is,ik,ig
    real(kind=DP) :: rsv
    integer    :: id_sname = -1
    integer    :: ist  !mpi
    real(kind=DP), pointer, dimension(:,:) :: s_loc2_mpi
    real(kind=DP) :: rsv_mpi

! ========================= added by K. Tagami =========== 11.0
    integer :: ismax
! ======================================================== 11.0

    call tstatc0_begin('m_Stress_loc ',id_sname)

    allocate(s_loc2_mpi(3,3)); s_loc2_mpi = 0.d0

! ========================= added by K. Tagami =========== 11.0
    if ( noncol ) then
       ismax = 1
    else
       ismax = nspin
    endif
! ======================================================== 11.0

    ist = ista_kngp
    if(ista_kngp==1) ist = 2
    if ( sw_add_vlocal_gzero == ON ) ist = ista_kngp

! ======================== modified by K. Tagami ======== 11.0
!    do is=1,nspin
    do is=1, ismax
! ======================================================== 11.0
       do i=1,3
          do j=1,3
             do it=1,ntyp
                do ik=1,kimg
!                   ist = ista_kngp
!                   if(ista_kngp==1) ist = 2
                   do ig = ist, iend_kngp  !for mpi
                      s_loc2_mpi(i,j) = s_loc2_mpi(i,j) &
                           &  - univol  &
                           &    * chgq_l(ig,ik,is)*zfm3_l(ig,it,ik)*psc_diff_l(ig,it) &
                           &    * 2.d0 * g(ig,i) &
                           &    * (g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j) &
                           &      +g(ig,3)*alinvt(3,j))
                   enddo
!xocl end spread sum(s_loc2)
                enddo
             enddo
          enddo
       enddo
    enddo

    call mpi_allreduce( s_loc2_mpi, s_loc2, 9, &
         &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
    deallocate(s_loc2_mpi)

    rsv_mpi = 0.d0
    rsv = 0.d0

    ist = ista_kngp
    if(ista_kngp==1) ist = 2
    if ( sw_add_vlocal_gzero == ON ) ist = ista_kngp

! ======================== modified by K. Tagami ======== 11.0
!    do is=1,nspin
    do is=1, ismax
! ======================================================== 11.0

      do it=1,ntyp
         do ik=1,kimg
            do ig = ist, iend_kngp  !for mpi
               rsv_mpi = rsv_mpi + chgsoft(ig,ik,is)*zfm3_l(ig,it,ik)*psc_l(ig,it)
              if(sw_screening_correction==ON) then
                  rsv_mpi = rsv_mpi - chgsoft(ig,ik,is)*zfm3_l(ig,it,ik)*ival(it)*Screening%phik(i)/univol
              end if
            enddo
         enddo
      enddo
    enddo

    call mpi_allreduce( rsv_mpi, rsv, 1, &
         &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )

    s_loc2 = s_loc2 - (etot1*totch + univol*rsv) * alinvt

    if(af == 0) then
       s_0 = s_loc1 + s_loc2
       stress_name = 'LO'
       call check_stress(nfout)
    else
       s_0 = s_loc2
       stress_name = 'L2'
       call check_stress(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_loc

  subroutine m_Stress_h(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ik,ig
    integer    :: id_sname = -1
    real(kind=DP) :: eh
    integer       :: ist  !mpi
    real(kind=DP) :: eh_mpi
    real(kind=DP), pointer, dimension(:,:) :: s_h2_mpi1, s_h2_mpi2

    call tstatc0_begin('m_Stress_h ',id_sname)

    ist = ista_kngp
    if(ista_kngp==1) ist = 2

    eh = 0.d0
    do ik = 1,kimg
! ========================== modified by K. Tagami ========== 11.0
!       if(nspin.eq.1) then
       if ( nspin.eq.1 .or. noncol ) then
! ============================================================ 11.0
          do i = ist, iend_kngp  !for mpi
             eh = eh + (chgq_l(i,ik,1)*grinv(i))**2
          enddo
       else if(nspin.eq.2) then
          do i = ist, iend_kngp  !for mpi
             eh = eh + ((chgq_l(i,ik,1)+chgq_l(i,ik,nspin))*grinv(i))**2
          enddo
       endif
    enddo

    do ik = 1, kimg
! ========================== modified by K. Tagami ========== 11.0
!       if(nspin.eq.1) then
       if ( nspin.eq.1 .or. noncol ) then
! ============================================================ 11.0
          do i = ist, iend_kngp  !for mpi
             eh = eh - 2 * chgq_l(i,ik,1)*chgsoft(i,ik,1)*grinv(i)**2
          enddo
       else if(nspin.eq.2) then
          do i = ist, iend_kngp  !for mpi
             eh = eh - 2 * (chgq_l(i,ik,1)+chgq_l(i,ik,nspin)) &
                  &  * (chgsoft(i,ik,1)+chgsoft(i,ik,nspin)) * grinv(i)**2
          enddo
      endif
    enddo

    if(npes > 1) then
       call mpi_allreduce(eh,eh_mpi,1 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
       eh = eh_mpi
    end if

    eh = univol*pai4*eh*0.5d0
    s_h2 = s_h2 + eh * alinvt

    allocate(s_h2_mpi1(3,3)); s_h2_mpi1 = 0.d0
    allocate(s_h2_mpi2(3,3)); s_h2_mpi2 = 0.d0

    do i=1,3
       do j=1,3
          do ik=1,kimg
! ========================== modified by K. Tagami ========== 11.0
!             if(nspin.eq.1) then
             if ( nspin.eq.1 .or. noncol ) then
! ============================================================ 11.0
                do ig = ist, iend_kngp  !for mpi
                   s_h2_mpi1(i,j) = s_h2_mpi1(i,j) &
                        &   + univol*pai4 &
                        &   *chgq_l(ig,ik,1)**2*grinv(ig)**4 &
                        &   *g(ig,i) &
                        &   *(g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                enddo
             else if(nspin.eq.2) then
                do ig = ist, iend_kngp  !for mpi
                   s_h2_mpi1(i,j) = s_h2_mpi1(i,j) &
                        &   + univol*pai4 &
                        &   *(chgq_l(ig,ik,1)+chgq_l(ig,ik,nspin))**2*grinv(ig)**4 &
                        &   *g(ig,i) &
                        &   *(g(ig,1)*alinvt(1,j)+g(ig,2)*alinvt(2,j)+g(ig,3)*alinvt(3,j))
                enddo
             endif
          enddo
       enddo
    enddo

    if(npes > 1) then
       call mpi_allreduce(s_h2_mpi1,s_h2_mpi2,9 &
            &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    else
       s_h2_mpi2 = s_h2_mpi1
    end if

    s_h2 = s_h2 + s_h2_mpi2
    deallocate(s_h2_mpi1); deallocate(s_h2_mpi2)

    if(af == 0) then
       s_0 = s_h1 +s_h2
       stress_name = 'HA'
       call check_stress(nfout)
    else
       s_0 = s_h2
       stress_name = 'H2'
       call check_stress(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_h

  subroutine m_Stress_ew(nfout)
    integer, intent(in) :: nfout

    s_0 = s_ew
    stress_name = 'EW'
    call check_stress(nfout)

  end subroutine m_Stress_ew

! =================================== KT_add ====================== 13.0B
  subroutine m_Stress_vdw(nfout)
    integer, intent(in) :: nfout
    integer :: i, j

    s_0 = s_vdw
    stress_name = 'VDW'

!    call check_stress(nfout)
    if(ipri >= 1) write(nfout,'(''  STRESS TENSOR'',5x,A2,/,(3f20.10))') &
         &     stress_name,((s_0(i,j),j=1,3),i=1,3)

  end subroutine m_Stress_vdw
! ================================================================= 13.0B

#ifndef DISABLE_VDWDF
  subroutine m_Stress_vdwdf(nfout)
    integer, intent(in) :: nfout
    integer :: i, j
!    s_vdwdf = s_cnl1 +s_cnl2 -s_cnl1_pc -s_cnl2_pc
    s_vdwdf = s_cnl1 +s_cnl2
    s_0 = s_vdwdf
!
    stress_name = 'VDWDF'

!    call check_stress(nfout)
    if(ipri >= 1) write(nfout,'(''  STRESS TENSOR'',5x,A2,/,(3f20.10))') &
         &     stress_name,((s_0(i,j),j=1,3),i=1,3)
  end subroutine m_Stress_vdwdf
#endif

  subroutine m_Stress_nl(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_nl ',id_sname)

    do i = 1, 3
    do j = 1, 3
    do ispin = 1, nspin, af+1
    do it = 1, ntyp
    do lmt1 = 1, ilmt(it)
    il1 = ltp(lmt1,it)
    im1 = mtp(lmt1,it)
    do lmt2 = lmt1, ilmt(it)
    il2 = ltp(lmt2,it)
    im2 = mtp(lmt2,it)
    if(ipaw(it)==0) then
        if(il1==il2 .and. im1==im2) then
           do ia = 1, natm
           if(ityp(ia)==it) then
              fac = 2.d0*dble(iwei(ia))
              if(lmt1==lmt2) fac = fac*0.5d0
              s_nl(i,j) = s_nl(i,j) &
                & + fac * dion(lmt1,lmt2,it) * hsrd(ia,lmt1,lmt2,ispin,i,j)

!              write(980,*) 'i j it = ', i, j, it
!              write(980,*) 'lmt1 kmt2 ia = ', lmt1, lmt2, ia
!              write(980,*) 'dion= ', dion(lmt1,lmt2,it)
!              write(980,*) 'hsrd= ', hsrd(ia,lmt1,lmt2,1,i,j), &
!                   &                 hsrd(ia,lmt1,lmt2,2,i,j)
           endif
           enddo
!
        endif
    else
! --- Check ---
       do ia = 1, natm
       if(ityp(ia)==it) then
          fac = 2.d0*dble(iwei(ia))
          if(lmt1==lmt2) fac = fac*0.5d0
          s_nl(i,j) = s_nl(i,j) &
            & + fac * dion_paw(lmt1,lmt2,ispin,ia) * hsrd(ia,lmt1,lmt2,ispin,i,j)
       endif
       enddo
    endif
    enddo  ! lmt2
    enddo  ! lmt1
    enddo  ! it
    enddo  ! ispin
    enddo  ! j
    enddo  ! i

    if(af == 0) then
       s_0 = s_nl
       stress_name = 'NL'
       call check_stress(nfout)
    else
       s_up = s_nl
       stress_name = 'NL'
       call check_af(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_nl

! =========================== added by K. Tagami ============ 11.0
  subroutine m_Stress_nl_nonclA(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    integer:: is

    call tstatc0_begin('m_Stress_nl_noncl ',id_sname)

    do i = 1, 3
       do j = 1, 3
          do it = 1, ntyp
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it)
                im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it)
                   im2 = mtp(lmt2,it)

                   if (ipaw(it)==0) then
                      if(il1==il2 .and. im1==im2) then
                         do ia = 1, natm
                            if(ityp(ia)==it) then
                               fac = 2.d0*dble(iwei(ia))
                               if(lmt1==lmt2) fac = fac*0.5d0
                               s_nl(i,j) = s_nl(i,j) +fac *dion(lmt1,lmt2,it) &
                                    &                     *hsrd(ia,lmt1,lmt2,1,i,j)

!                               write(980,*) 'i j it = ', i, j, it
!                               write(980,*) 'lmt1 kmt2 ia = ', lmt1, lmt2, ia
!                               write(980,*) 'dion= ', dion(lmt1,lmt2,it)
!                               write(980,*) 'hsrd= ', hsrd(ia,lmt1,lmt2,1,i,j), &
!                                    &                 hsrd(ia,lmt1,lmt2,2,i,j)

                            endif
                         enddo
                      endif
                   else

                      do ia = 1, natm
                         if(ityp(ia)==it) then
                            fac = 2.d0*dble(iwei(ia))
                            if(lmt1==lmt2) fac = fac*0.5d0

                            Do is=1, ndim_magmom
                               s_nl(i,j) = s_nl(i,j) +fac *dion_paw(lmt1,lmt2,is,ia)&
                                    &                     *hsrd(ia,lmt1,lmt2,is,i,j)
                            End do
                         endif
                      enddo
                   endif
                enddo  ! lmt2
             enddo  ! lmt1
          enddo  ! it
       enddo  ! j
    enddo  ! i

    s_0 = s_nl;        stress_name = 'NL'
    call check_stress(nfout)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_nl_nonclA

  subroutine m_Stress_nl_nonclB(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1
    integer:: is

    real(kind=DP) :: fac, c1, c2
    real(kind=DP), allocatable ::  hsr_ssrep( :,:,:,: )
    real(kind=DP), allocatable ::  hsi_ssrep( :,:,:,: )

    call tstatc0_begin('m_Stress_nl_noncl ',id_sname)

    allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) )
    allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) )
!
    hsr_ssrep = 0.0d0;  hsi_ssrep = 0.0d0

    do i = 1, 3
       do j = 1, 3
!
          call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, &
               &                           hsrd(:,:,:,:,i,j), hsid(:,:,:,:,i,j ), &
               &                           hsr_ssrep, hsi_ssrep )
!
          do it = 1, ntyp
!
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it)
                im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it)
                   im2 = mtp(lmt2,it)

                   do ia = 1, natm
                      if ( ityp(ia) /= it ) cycle
                      fac = 2.d0*dble(iwei(ia))
                      if(lmt1==lmt2) fac = fac*0.5d0

                      Do is=1, ndim_chgpot
                         c1 = real(  dion0_noncl( lmt1,lmt2,is,ia ) ) &
                              &  *hsr_ssrep( ia,lmt1,lmt2,is )
                         c2 = aimag(  dion0_noncl( lmt1,lmt2,is,ia ) ) &
                              &  *hsi_ssrep( ia,lmt1,lmt2,is )
                         s_nl(i,j) = s_nl(i,j) +fac *( c1+c2 )
                      End do
                   enddo
                enddo  ! lmt2
             enddo  ! lmt1
          enddo  ! it
       enddo  ! j
    enddo  ! i
!
    deallocate( hsr_ssrep );   deallocate( hsi_ssrep )
!
    s_0 = s_nl;        stress_name = 'NL'
    call check_stress(nfout)

    call tstatc0_end(id_sname)

  end subroutine m_Stress_nl_nonclB
! ================================================================= 11.0

  subroutine m_Stress_or(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_or ',id_sname)

    do i = 1, 3
    do j = 1, 3
    do ispin = 1, nspin, af+1
    do it = 1, ntyp
    do lmt1 = 1, ilmt(it)
    il1 = ltp(lmt1,it)
    im1 = mtp(lmt1,it)
    do lmt2 = lmt1, ilmt(it)
    il2 = ltp(lmt2,it)
    im2 = mtp(lmt2,it)
    if(il1==il2 .and. im1==im2) then
       do ia = 1, natm
       if(ityp(ia)==it) then
          fac = 2.d0*dble(iwei(ia))
          if(lmt1==lmt2) fac = fac*0.5d0
          s_or(i,j) = s_or(i,j) &
            & + fac * q(lmt1,lmt2,it) * esrd(ia,lmt1,lmt2,ispin,i,j)
       endif
       enddo
    endif
    enddo  ! lmt2
    enddo  ! lmt1
    enddo  ! it
    enddo  ! ispin
    enddo  ! j
    enddo  ! i

    if(af == 0) then
       s_0 = s_or
       stress_name = 'OR'
       call check_stress(nfout)
    else
       s_up = s_or
       stress_name = 'OR'
       call check_af(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_or

! ================================= added by K. Tagami ============= 11.0
  subroutine m_Stress_or_nonclA(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_or_noncl ',id_sname)

    do i = 1, 3
       do j = 1, 3

          do it = 1, ntyp
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it)
                im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it)
                   im2 = mtp(lmt2,it)

                   if(il1==il2 .and. im1==im2) then
                      do ia = 1, natm
                         if(ityp(ia)==it) then
                            fac = 2.d0*dble(iwei(ia))
                            if(lmt1==lmt2) fac = fac*0.5d0
                            s_or(i,j) = s_or(i,j) &
                                 & + fac * q(lmt1,lmt2,it) * esrd(ia,lmt1,lmt2,1,i,j)
                         endif
                      enddo
                   endif
                enddo  ! lmt2
             enddo  ! lmt1
          enddo  ! it
       enddo  ! j
    enddo  ! i

    s_0 = s_or;        stress_name = 'OR'
    call check_stress(nfout)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_or_nonclA

  subroutine m_Stress_or_nonclB(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
    integer    :: id_sname = -1

    integer :: is
    real(kind=DP) :: fac, c1, c2
    real(kind=DP), allocatable ::  hsr_ssrep( :,:,:,: )
    real(kind=DP), allocatable ::  hsi_ssrep( :,:,:,: )

    call tstatc0_begin('m_Stress_or_noncl ',id_sname)

    allocate( hsr_ssrep( natm, nlmt, nlmt, ndim_chgpot ) )
    allocate( hsi_ssrep( natm, nlmt, nlmt, ndim_chgpot ) )
!
    hsr_ssrep = 0.0d0;  hsi_ssrep = 0.0d0

    do i = 1, 3
       do j = 1, 3
!
          call m_ES_MagMom_To_DensMat_hsr( natm, nlmt, &
               &                           esrd(:,:,:,:,i,j), esid(:,:,:,:,i,j), &
               &                           hsr_ssrep, hsi_ssrep )
!
          do it = 1, ntyp
             do lmt1 = 1, ilmt(it)
                il1 = ltp(lmt1,it)
                im1 = mtp(lmt1,it)
                do lmt2 = lmt1, ilmt(it)
                   il2 = ltp(lmt2,it)
                   im2 = mtp(lmt2,it)

                   if (il1/=il2) cycle

                   do ia = 1, natm
                      if (ityp(ia)/=it) cycle

                      if ( ityp(ia) /= it ) cycle
                      fac = 2.d0*dble(iwei(ia))
                      if(lmt1==lmt2) fac = fac*0.5d0

                      Do is=1, ndim_chgpot
                         c1 = real(  q_noncl( lmt1,lmt2,is,it ) ) &
                              &  *hsr_ssrep( ia,lmt1,lmt2,is )
                         c2 = aimag(  q_noncl( lmt1,lmt2,is,it ) ) &
                              &  *hsi_ssrep( ia,lmt1,lmt2,is )
                         s_or(i,j) = s_or(i,j) +fac *( c1+c2 )
                      End do

                   end do
                enddo  ! lmt2
             enddo  ! lmt1
          enddo  ! it
       enddo  ! j
    enddo  ! i
!
    deallocate( hsr_ssrep );   deallocate( hsi_ssrep )
!
    s_0 = s_or;        stress_name = 'OR'
    call check_stress(nfout)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_or_nonclB
! =============================================================== 11.0

! =========================== Added by K. Tagami =========== p1100
  subroutine m_Stress_Hub(nfout)
    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,ispin,lmt1,lmt2,il1,il2,im1,im2
!
    integer ih, l1p
!
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_Hub ',id_sname)

    do i = 1, 3
       do j = 1, 3
          do ispin = 1, nspin, af+1
             do ia = 1, natm
                ih = ihubbard(ia)
                if(ih == 0) cycle
                it = ityp(ia)
                l1p = proj_attribute(ih)%l+1
                do lmt1 = 1, ilmt(it)
                   il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                   if(il1 /= l1p) cycle
                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                      if(il2 /= l1p) cycle
      !!$if(im1 /= im2) cycle
                      fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                      s_hub(i,j) = s_hub(i,j) &
                           & + fac * dhub(lmt1,lmt2,ia,ispin)  &
                           &       * hsrd(ia,lmt1,lmt2,ispin,i,j)
                   end do
                end do
             end do
          end do
       enddo
    enddo
! --
!    if ( nspin==1 ) then
!       s_hub =s_hub * 2.0d0
!    endif
!    s_hub = 0.0d0
    ! --
    if(af == 0) then
       s_0 = s_hub
       stress_name = 'HU'
       call check_stress(nfout)
    else
       s_up = s_hub
       stress_name = 'HU'
       call check_af(nfout)
    endif

    call tstatc0_end(id_sname)
  end subroutine m_Stress_Hub

! ====================================================== p1100

! ============================ added by K. Tagami ============ 11.0
  subroutine m_Stress_Hub_nonclA(nfout)

    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,is,lmt1,lmt2,il1,il2,im1,im2
!
    integer ih, l1p
!
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_Hub_noncl ',id_sname)

    do i = 1, 3
       do j = 1, 3

          Do is=1, ndim_magmom
             do ia = 1, natm
                ih = ihubbard(ia)
                if(ih == 0) cycle
                it = ityp(ia)
                l1p = proj_attribute(ih)%l+1

                do lmt1 = 1, ilmt(it)
                   il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                   if(il1 /= l1p) cycle

                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                      if(il2 /= l1p) cycle
      !!$if(im1 /= im2) cycle
                      fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                      s_hub(i,j) = s_hub(i,j) &
                           & + fac * dhub(lmt1,lmt2,ia,is)  &
                           &       * hsrd(ia,lmt1,lmt2,is,i,j)
                   end do
                end do
             end do
          End do

       enddo
    enddo
!
    s_0 = s_hub;     stress_name = 'HU'
    call check_stress(nfout)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_Hub_nonclA

  subroutine m_Stress_Hub_nonclB(nfout)

    integer, intent(in) :: nfout
    integer    :: i,j,ia,it,is,lmt1,lmt2,il1,il2,im1,im2
!
    integer ih, l1p
!
    integer    :: id_sname = -1
    real(kind=DP) :: fac

    call tstatc0_begin('m_Stress_Hub_noncl ',id_sname)

    do i = 1, 3
       do j = 1, 3

          Do is=1, ndim_magmom
             do ia = 1, natm
                ih = ihubbard(ia)
                if(ih == 0) cycle
                it = ityp(ia)
                l1p = proj_attribute(ih)%l+1

                do lmt1 = 1, ilmt(it)
                   il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)
                   if(il1 /= l1p) cycle

                   do lmt2 = lmt1, ilmt(it)
                      il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)
                      if(il2 /= l1p) cycle
      !!$if(im1 /= im2) cycle
                      fac = 2.d0*iwei(ia); if(lmt1 == lmt2) fac = fac*0.5d0
                      s_hub(i,j) = s_hub(i,j) &
                           & + fac * dhub(lmt1,lmt2,ia,is)  &
                           &       * hsrd(ia,lmt1,lmt2,is,i,j) &
                           & + fac * dhub_aimag(lmt1,lmt2,ia,is)  &
                           &       * hsid(ia,lmt1,lmt2,is,i,j)
                   end do
                end do
             end do
          End do

       enddo
    enddo
!
    s_0 = s_hub;     stress_name = 'HU'
    call check_stress(nfout)

    call tstatc0_end(id_sname)
  end subroutine m_Stress_Hub_nonclB
! ================================================================ 11.0

  subroutine check_stress(nfout)
    integer, intent(in) :: nfout
    integer             :: i,j,k,l,iopr

    s_1 = 0.d0
    s_2 = 0.d0
    do i=1,3
    do j=1,3
      do k=1,3
      s_1(i,j)=s_1(i,j)-s_0(i,k)*altv(j,k)/univol
      enddo
    enddo
    enddo

    do i=1,3
    do j=1,3
      do k=1,3
      do l=1,3
        do iopr=1,nopr
        s_2(i,j)=s_2(i,j) + op(i,k,iopr)*op(j,l,iopr)*s_1(k,l)
        enddo
      enddo
      enddo
    s_2(i,j)=s_2(i,j)/nopr
    enddo
    enddo

    if(ipri >= 1) write(nfout,'(''  STRESS TENSOR'',5x,A2,/,(3f20.10))') &
         &     stress_name,((s_2(i,j),j=1,3),i=1,3)
  end subroutine check_stress

  subroutine check_af(nfout)
    integer, intent(in) :: nfout
    integer             :: i,j,k,l,iopr

    s_1 = 0.d0
    s_2 = 0.d0
    s_dn = 0.d0

    do i=1,3
    do j=1,3
      do k=1,3
      s_1(i,j)=s_1(i,j)-s_up(i,k)*altv(j,k)/univol
      enddo
    enddo
    enddo

    iopr = nopr +af
    do i=1,3
    do j=1,3
      do k=1,3
      do l=1,3
        s_dn(i,j)=s_dn(i,j) + op(i,k,iopr)*op(j,l,iopr)*s_1(k,l)
      enddo
      enddo
    enddo
    enddo

    s_1 = s_1 + s_dn

    do i=1,3
    do j=1,3
      do k=1,3
      do l=1,3
        do iopr=1,nopr
        s_2(i,j)=s_2(i,j) + op(i,k,iopr)*op(j,l,iopr)*s_1(k,l)
        enddo
      enddo
      enddo
    s_2(i,j)=s_2(i,j)/nopr
    enddo
    enddo

    if(ipri >= 1) write(nfout,'(''  STRESS TENSOR'',5x,A2,/,(3f20.10))') &
         &     stress_name,((s_2(i,j),j=1,3),i=1,3)
  end subroutine check_af

  subroutine m_Stress_sum(nfout)
    integer, intent(in) :: nfout
    integer             :: i,j,k,l,iopr
    real(kind=DP) :: c1
    real(kind=DP) :: altv_brav(3,3)

    s_tmp=0.d0
    s=0.d0
    if(af==0) then
       s = s_ke + s_xc1 + s_xc2 + s_gga1 + s_gga2 + s_loc1 + s_loc2 &
         &  + s_h1 + s_h2 + s_nl - s_or + s_ew

! ================================== Added by K. Tagami ======= p1100
       s = s + s_hub
! ============================================================== p1100

       do i=1,3
         do j=1,3
           do k=1,3
             s_tmp(i,j) = s_tmp(i,j) - s(i,k) * altv(j,k) / univol
           enddo
         enddo
       enddo
    else
       s_up = s_ke + s_xc1 + s_gga1 + s_gga2 + s_loc1 &
         &  + s_h1 + s_nl - s_or

! ================================== Added by K. Tagami ======= p1100
       s_up = s_up + s_hub
! ============================================================= p1100

       do i=1,3
         do j=1,3
           do k=1,3
             s_tmp(i,j) = s_tmp(i,j) - s_up(i,k) * altv(j,k) / univol
           enddo
         enddo
       enddo
       iopr = nopr + af
       s_dn = 0.d0
       do i=1,3
         do j=1,3
           do k=1,3
             do l=1,3
               s_dn(i,j) = s_dn(i,j) + op(i,k,iopr)*op(j,l,iopr)*s_tmp(k,l)
             enddo
           enddo
         enddo
       enddo
       s_tmp = s_tmp + s_dn
       s = s_xc2 + s_loc2 + s_h2 + s_ew

       do i=1,3
         do j=1,3
           do k=1,3
             s_tmp(i,j) = s_tmp(i,j) - s(i,k) * altv(j,k) / univol
           enddo
         enddo
       enddo
    endif

! ================================== KT_add  =============== 13.0B
    if ( sw_vdw_correction == ON ) s_tmp = s_tmp +s_vdw
! ========================================================== 13.0B
!    if ( xctype == "vdwdf" ) s = s +s_vdwdf

    s = 0.d0
    do i=1,3
      do j=1,3
        do k=1,3
          do l=1,3
            do iopr=1,nopr
              s(i,j) = s(i,j) + op(i,k,iopr)*op(j,l,iopr)*s_tmp(k,l)
            enddo
          enddo
        enddo
      enddo
    enddo
    s = s / nopr

    if ( xctype == "vdwdf" ) s = s +s_vdwdf

    if(ipri >= 1) write(nfout,'(''  STRESS TENSOR'',/,(3f20.10))') &
         &     ((s(i,j),j=1,3),i=1,3)

! ============================= KT_add ======================== 13.0B
    if (ipri >= 1) then
       write(nfout,*) '-------------------------- '
       write(nfout,'(''  External STRESS TENSOR'',/,(3f20.10))') &
            &            ( ( -external_stress(i,j),j=1,3),i=1,3)
       write(nfout,'(''  Total STRESS TENSOR'',/,(3f20.10))') &
            &            ( ( s(i,j)-external_stress(i,j),j=1,3),i=1,3)
       write(nfout,*) '-------------------------- '
    endif
! ============================================================= 13.0B

! ================================ KT_mod ================== 13.0B
!    stressmx = abs(s(1,1))
!    do i=1,3
!       do j=1,3
!          if(abs(s(j,i))>stressmx) stressmx = abs(s(j,i))
!       enddo
!    enddo
!    currstress(:,:) = s(:,:)
!    if(sw_uniform==ON)then
!       stressmx = abs(currstress(1,1)+currstress(2,2)+currstress(3,3))/3.0d0
!    endif

    if(.not.m_Stress_in_correction()) then
      currstress = s -external_stress
      call cal_stressmx()
    endif
!

! ================================================================= 13.0B

  end subroutine m_Stress_sum

  subroutine cal_stressmx()
    real(kind=DP) :: s_tmp(3,3)
    real(kind=DP) :: altv_brav(3,3)
    real(kind=DP) :: c1
    integer :: i,j
    stressmx = abs(currstress(1,1))
    do i=1,3
       do j=1,3
          if(abs(currstress(j,i))>stressmx) stressmx = abs(currstress(j,i))
       enddo
    enddo
    if ( sw_uniform == ON ) then
       stressmx = abs(currstress(1,1)+currstress(2,2)+currstress(3,3))/3.0d0
    endif

! === KT_add === 13.1AS
    if ( sw_neglect_stress_offdiagonal == ON ) then
       s_tmp = 0.0d0
       s_tmp(1,1) = currstress(1,1)
       s_tmp(2,2) = currstress(2,2)
       s_tmp(3,3) = currstress(3,3)

       stressmx = 0.0d0
       do i=1,3
          c1 = sqrt( s_tmp(i,1)**2 +s_tmp(i,2)**2 +s_tmp(i,3)**2 )
          if ( c1 > stressmx ) stressmx = c1
       enddo
    endif

    if ( sw_fix_lattice_angles == ON .or. sw_fix_lattice_lengths == ON ) then
       s_tmp = currstress

       if ( sw_fix_lattice_shapes == ON ) then
          if ( fix_shape_ab_plane ) then
             s_tmp(1,1) = ( currstress(1,1)+currstress(2,2) )/2.0d0
             s_tmp(2,2) = ( currstress(1,1)+currstress(2,2) )/2.0d0
             s_tmp(1,2) = 0.0d0
             s_tmp(2,1) = 0.0d0
          endif
          if ( fix_shape_ac_plane ) then
             s_tmp(1,1) = ( currstress(1,1)+currstress(3,3) )/2.0d0
             s_tmp(3,3) = ( currstress(1,1)+currstress(3,3) )/2.0d0
             s_tmp(1,3) = 0.0d0
             s_tmp(3,1) = 0.0d0
          endif
          if ( fix_shape_bc_plane ) then
             s_tmp(2,2) = ( currstress(2,2)+currstress(3,3) )/2.0d0
             s_tmp(3,3) = ( currstress(2,2)+currstress(3,3) )/2.0d0
             s_tmp(2,3) = 0.0d0
             s_tmp(3,2) = 0.0d0
          endif
       endif

       call get_latvec_from_prim_to_brav( p2bmat, altv, altv_brav )

       s_tmp = projected_stress_along_latvecs( altv_brav, s_tmp )
       stressmx = 0.0d0

       if ( sw_fix_lattice_angles == ON ) then
          s_tmp = constrain_on_stress_type1( altv_brav, s_tmp )
       endif
       if ( sw_fix_lattice_lengths == ON ) then
          s_tmp = constrain_on_stress_type2( altv_brav, s_tmp )
       endif

       do i=1,3
          c1 = sqrt( s_tmp(i,1)**2 +s_tmp(i,2)**2 +s_tmp(i,3)**2 )
          if ( c1 > stressmx ) stressmx = c1
       enddo
    endif
! ============= 13.1AS
  end subroutine cal_stressmx


  function m_Stress_get_stressmx() result(res)
     real(kind=DP) :: res
     res = stressmx
  end function m_Stress_get_stressmx

  subroutine m_Stress_set_stressmx(strmx)
    real(kind=DP), intent(in) :: strmx
    stressmx = strmx
  end subroutine m_Stress_set_stressmx

  function m_Stress_get_curr_stress() result(res)
     real(kind=DP), dimension(3,3) :: res
     res = currstress
  end function m_Stress_get_curr_stress
!$$#endif

! === KT_add ==== 13.1AS
function projected_stress_along_latvecs(altv_in,stress) result(ret)
  real(kind=DP), dimension(3,3), intent(in) :: stress, altv_in
  real(kind=DP), dimension(3,3) :: ret

  integer :: i
  real(kind=DP) :: normalized_lat_vec(3,3), c1
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

  ret = 0.0d0

  do i=1,3
     ret(1,i) = dot_product(stress(i,:),normalized_lat_vec(:,1)) ! a-vec
     ret(2,i) = dot_product(stress(i,:),normalized_lat_vec(:,2)) ! b-vec
     ret(3,i) = dot_product(stress(i,:),normalized_lat_vec(:,3)) ! c-vec
  enddo
end function projected_stress_along_latvecs

function constrain_on_stress_type1(altv_in,s_in) result(s_out)
  real(kind=DP), intent(in) :: s_in(3,3), altv_in(3,3)
  real(kind=DP) :: s_out(3,3)

  integer :: i
  real(kind=DP) :: c1, c2, c3
  real(kind=DP) :: normalized_lat_vec(3,3), vectmp(3), vectmp2(3)
!
  s_out = s_in
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

! -----------------------------------
  if ( fix_angle_alpha .and. fix_angle_beta .and. fix_angle_gamma ) then
     Do i=1, 3
        c1 = dot_product( s_out(i,:), normalized_lat_vec(:,i) )
        s_out(i,:) = c1 *normalized_lat_vec(:,i)
     End do
     return
  endif

! -----------------------------------
  if ( fix_angle_alpha ) then
     vectmp(:) = normalized_lat_vec(:,2) + normalized_lat_vec(:,3)

     Do i=2, 3
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( s_out(i,:), vectmp2(:) )
        s_out(i,:) = s_out(i,:) - c3 *vectmp2(:)
     End do
  endif

! -----------------------------------
  if ( fix_angle_beta ) then
     vectmp(:) = normalized_lat_vec(:,1) + normalized_lat_vec(:,3)

     Do i=1, 3, 2
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( s_out(i,:), vectmp2(:) )
        s_out(i,:) = s_out(i,:) - c3 *vectmp2(:)
     End do
  endif

! -----------------------------------
  if ( fix_angle_gamma ) then
     vectmp(:) = normalized_lat_vec(:,1) + normalized_lat_vec(:,2)

     Do i=1, 2
        c1 = dot_product( vectmp(:), normalized_lat_vec(:,i) )
        vectmp2(:) = vectmp(:) - c1*normalized_lat_vec(:,i)

        c2 = sqrt( vectmp2(1)**2 + vectmp2(2)**2 + vectmp2(3)**2 )
        vectmp2 = vectmp2 /c2

        c3 = dot_product( s_out(i,:), vectmp2(:) )
        s_out(i,:) = s_out(i,:) - c3 *vectmp2(:)
     End do

  endif

end function constrain_on_stress_type1

function constrain_on_stress_type2(altv_in,s_in) result(s_out)
  real(kind=DP), intent(in) :: s_in(3,3), altv_in(3,3)
  real(kind=DP) :: s_out(3,3)

  integer :: i
  real(kind=DP) :: c1, c2, c3
  real(kind=DP) :: normalized_lat_vec(3,3), vectmp(3), vectmp2(3)
!
  s_out = s_in
!
  Do i=1, 3
     c1 = sqrt( altv_in(1,i)**2 +altv_in(2,i)**2 + altv_in(3,i)**2 )
     normalized_lat_vec(:,i) = altv_in(:,i) / c1
  End do

! -----------------------------------
  if ( fix_length_a ) then
     c1 = dot_product( s_out(1,:), normalized_lat_vec(:,1) )
     s_out(1,:) = s_out(1,:) -c1 *normalized_lat_vec(:,1)
  endif
  if ( fix_length_b ) then
     c1 = dot_product( s_out(2,:), normalized_lat_vec(:,2) )
     s_out(2,:) = s_out(2,:) -c1 *normalized_lat_vec(:,2)
  endif
  if ( fix_length_c ) then
     c1 = dot_product( s_out(3,:), normalized_lat_vec(:,3) )
     s_out(3,:) = s_out(3,:) -c1 *normalized_lat_vec(:,3)
  endif
end function constrain_on_stress_type2
! ====== 13.1AS


subroutine m_Stress_correction(nfout,cntn)
  use m_IterationNumbers, only : iteration_stress_correction
  use m_Total_Energy, only: etotal0
  use m_Crystal_Structure, only : univol
  use m_Control_Parameters, only : gmax_org,decut_stress_correction,external_stress
  implicit none
  integer, intent(in) :: nfout
  logical, intent(in), optional :: cntn
  real(kind=DP) :: ecut,idelta_v,e1,e2,s1,s2
  logical :: from_cntn
  logical, save :: evaluated_correction = .false.
  if (evaluated_correction) return
  from_cntn = .false.
  if(present(cntn)) from_cntn = cntn
  if(iteration_stress_correction>=1 .and. iteration_stress_correction <= 2)then
    etot3(iteration_stress_correction) = etotal0
    univol0 = univol
  endif
  if(iteration_stress_correction == 2 .and. .not. from_cntn .or. &
  & (from_cntn.and.iteration_stress_correction>=3)) then !eval stress correction
    ecut = gmax_org*gmax_org
    idelta_v = -ecut/(1.5d0*univol0*decut_stress_correction)
    e1 = -etot3(1)
    e2 =  etot3(2)
    s1 = -e1*idelta_v
    s2 = -e2*idelta_v
    stress_correction = (s1+s2)*0.5d0
    if(printable) then
       write(nfout,'(a,f20.15)')  ' !** Pulay stress : ',stress_correction
    endif
    external_stress(1,1) = stress_correction
    external_stress(2,2) = stress_correction
    external_stress(3,3) = stress_correction
    estress_has_been_set = .true.
    evaluated_correction = .true.
  endif
end subroutine m_Stress_correction

function m_Stress_get_correction() result(res)
  implicit none
  real(kind=DP) :: res
  res = stress_correction
  return
end function m_Stress_get_correction

logical function in_correction0()
  use m_Const_Parameters, only : OFF
  use m_Control_Parameters, only : sw_stress_correction
  use m_IterationNumbers, only : iteration_stress_correction
  implicit none
  in_correction0 = in_correction1(2)
end function in_correction0

logical function in_correction1(ncount)
  use m_Const_Parameters, only : OFF
  use m_Control_Parameters, only : sw_stress_correction
  use m_IterationNumbers, only : iteration_stress_correction
  implicit none
  integer, intent(in) :: ncount
  if(sw_stress_correction == OFF) then
    in_correction1 = .false.
    return
  endif
  in_correction1 = iteration_stress_correction <= ncount
  return
end function in_correction1

function m_Stress_decut_for_correction() result(ret)
  use m_Const_Parameters, only : OFF,DP
  use m_Control_Parameters, only : sw_stress_correction, decut_stress_correction
  use m_IterationNumbers, only : iteration_stress_correction
  implicit none
  real(kind=DP) :: ret
  if(sw_stress_correction == OFF)then
     ret = 0.d0
     return
  endif
  if(iteration_stress_correction == 1) then
     ret = -decut_stress_correction
  else if (iteration_stress_correction == 2) then
     ret = +decut_stress_correction
  else
     ret = 0.d0
  endif
  return
end function m_Stress_decut_for_correction

subroutine m_Stress_wd_cdata_4_correction(nfcntn)
  use m_IterationNumbers, only : iteration_stress_correction
  implicit none
  integer, intent(in) :: nfcntn
  if(mype==0)then
     write(nfcntn,*) tag_stress_correction
     write(nfcntn,'(i8)') iteration_stress_correction
     write(nfcntn,'(4f20.10)') etot3(1),etot3(2),etot3(3),univol0
     write(nfcntn,'(l2)') estress_has_been_set
     write(nfcntn,'(3f20.10)') currstress(1,1),currstress(1,2),currstress(1,3)
     write(nfcntn,'(3f20.10)') currstress(2,1),currstress(2,2),currstress(2,3)
     write(nfcntn,'(3f20.10)') currstress(3,1),currstress(3,2),currstress(3,3)
  endif
end subroutine m_Stress_wd_cdata_4_correction

subroutine m_Stress_rd_cdata_4_correction(nfout,nfcntn)
  use m_IterationNumbers, only : iteration_stress_correction
  implicit none
  integer, intent(in) :: nfout,nfcntn
  integer,      parameter   :: len_str = 132
  character(len=len_str)    :: str
  logical             :: tag_is_found, EOF_reach
  integer :: ierr
  if(mype==0)then
     call rewind_to_tag0(nfcntn,len(tag_stress_correction),tag_stress_correction &
          &, EOF_reach, tag_is_found, str,len_str)
     if(.not.tag_is_found)then
       write(nfout,'(a)') '!** '//trim(tag_stress_correction)//' does not exist'
     else
       read(nfcntn,*) iteration_stress_correction
       read(nfcntn,*) etot3(1),etot3(2),etot3(3),univol0
       read(nfcntn,*) estress_has_been_set
       read(nfcntn,*) currstress(1,1),currstress(1,2),currstress(1,3)
       read(nfcntn,*) currstress(2,1),currstress(2,2),currstress(2,3)
       read(nfcntn,*) currstress(3,1),currstress(3,2),currstress(3,3)
     endif
  endif
  if(npes>1)then
     call mpi_bcast(tag_is_found,1,mpi_logical,0,MPI_CommGroup,ierr)
     if(tag_is_found)then
        call mpi_bcast(iteration_stress_correction,1,mpi_integer,0,MPI_CommGroup,ierr)
        call mpi_bcast(etot3,3,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(univol0,1,mpi_double_precision,0,MPI_CommGroup,ierr)
        call mpi_bcast(estress_has_been_set,1,mpi_logical,0,MPI_CommGroup,ierr)
        call mpi_bcast(currstress,9,mpi_double_precision,0,MPI_CommGroup,ierr)
     endif
  endif
end subroutine m_Stress_rd_cdata_4_correction

end module m_Stress
