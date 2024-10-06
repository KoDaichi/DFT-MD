!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver: 801)
!
!  MODULE: m_ES_nonlocal
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, April/15/2006, September/02/2008
!  FURTHER MODIFICATION: T. Yamasaki, T. Uda and T. Ohno, September 2009 (MGS_DGEMM)
!  FURTHER MODIFICATION: T. Yamasaki and T. Yamamoto,   October 2009  (NONLOCAL_DGEMM)
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

!   This module has been revised for the GAMMA point (k=(0,0,0)) by T. Yamasaki
!  in April 2006. Number of operations for the Gamma point have been tremendously
!  reduced in subroutines of m_ES_betar_dot_Wfs, m_ES_Vnonlocal_W, and
!  m_ES_modified_gram_schmidt.
!
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
#endif
#ifdef __TIMER_DGEMM__
#   define __TIMER_DGEMM_START(a)  call timer_sta(a)
#   define __TIMER_DGEMM_STOP(a)   call timer_end(a)
#else
#   define __TIMER_DGEMM_START(a)
#   define __TIMER_DGEMM_STOP(a)
#endif
#ifdef __TIMER_DO__
#   define __TIMER_DO_START(a)   call timer_sta(a)
#   define __TIMER_DO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_DO_START(a)
#   define __TIMER_DO_STOP(a)
#endif
#ifdef __TIMER_COMM__
#   define __TIMER_COMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_COMM_START(a)       call timer_sta(a)
#   define __TIMER_COMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_COMM_START_w_BARRIER(str,a)
#   define __TIMER_COMM_START(a)
#   define __TIMER_COMM_STOP(a)
#endif

#ifndef SX
#define DGEMM__       DGEMM
#endif

#ifndef NO_NONLOCAL_DGEMM
#define NONLOCAL_DGEMM
#endif

module m_ES_nonlocal
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process &
       &                         , iteration, nkgroup
  use m_NonLocal_Potential, only : snl,snl_add,phig,paog, norm_phig &
!       &                         , nmesh_rs,nmesh_rs_max,meshx_rs,meshy_rs,meshz_rs,meshxyz_rs,snl_rs &
!       &                         , meshxyz_rs_conjg,nmesh_rs_h,meshxyz_rs_h,map_h,nmesh_rs_max_h,snl_rs_h
       &                         , snl_rs, snl_rs_h &
       &                         , i_exp_snl, AtaulmaG, BtaulmaG
  use m_Realspace,          only : nmesh_rs,nmesh_rs_max,meshx_rs,meshy_rs,meshz_rs,meshxyz_rs,meshxyz_rs_conjg &
       &                         , nmesh_rs_h,meshxyz_rs_h,map_h,nmesh_rs_max_h
  use m_PlaneWaveBasisSet,  only : ngabc,kg1,nbase,nbmx,iba,igf
  use m_PseudoPotential,    only : ilmt,nlmtt,nlmta,lmta,lmtt,ltp,mtp,q,dion,taup &
       &                         , ilmt_phi,nlmt_phi,nlmtt_phi,nlmta_phi &
       &                         , lmta_phi,lmtt_phi,ltp_phi, mtp_phi, taup_phi &
       &                         , ilmt_add,nlmt_add,nlmtt_add,nlmta_add &
       &                         , lmta_add,lmtt_add,ltp_add,mtp_add &
       &                         , ilmt_pao,nlmt_pao,nlmtt_pao,nlmta_pao &
       &                         , lmta_pao,lmtt_pao,ltp_pao,mtp_pao,ibpao &
       &                         , modnrm,nac,fqwei,nlmta1,nlmta2 &
       &                         , m_PP_include_vanderbilt_pot &
       &                         , ipaw, dion_paw, q_phirt_pw, ia2ia_symmtry_op_inv &
       &                         , irorb, nrorb, crorb
  use m_Crystal_Structure,  only : op, nopr, univol
  use m_Kpoints,            only : kv3, kv3_ek,k_symmetry,vkxyz, sw_force_kpt_inside_bz
  use m_Ionic_System,       only : ntyp, natm, iwei, ityp, pos
  use m_FFT,                only : nfft &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work &
       &                         , fft_box_size_WF, m_FFT_WF
  use m_Files,              only : nfout
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only : nspin,ipri, ipribetar &
       &                         , kimg, neg, meg, af &
       &                         , m_CtrlP_cachesize &
       &                         , sw_use_add_proj &
       &                         , flag_mpi_g_dot_r, flag_mpi_g_dot_r_k &
       &                         , nblocksize_betar_dot_wfs_nlmta  &
       &                         , nblocksize_betar_nlmta_is_given &
       &                         , nblocksize_betar_dot_wfs, nblocksize_vnonlocal_w &
       &                         , nblocksize_betar_is_given &
       &                         , nblocksize_vnonlocal_is_given &
       &                         , nblocksize_betar_npe_is_given &
       &                         , nblocksize_betar_dot_wfs_npe  &
       &                         , sw_rspace,sw_save_memory, m_CtrlP_realspace_integ_OK &
       &                         , nblocksize_rspace_betar, sw_rspace_v &
       &                         , nblocksize_rspace_v &
       &                         , sw_betar_dot_wfs_exp, sw_precalculate_phase_vnonlocal
  use m_Const_Parameters,   only : DP, CMPLDP, SKIP, EXECUT, ON, OFF, PAI2 &
       &                         , GAMMA, MAPPED, NOTMAPPED, BUCS, DIRECT, ELECTRON, INVERSE
  use m_Parallelization,    only : mype &
       &                         , myrank_e,map_e,ista_e,iend_e,istep_e &
       &                         , map_z,np_e,myrank_k,map_k,ista_k,iend_k &
       &                         , ista_snl, iend_snl &
       &                         , ng_nbmx, ista_nbmx, mpi_nbmx_world &
       &                         , np_nbmx, mp_nbmx &
       &                         , ng_nbmx_k, ista_nbmx_k, mpi_nbmx_world_k &
!!$       &                         , np_nbmx_k, mp_nbmx_k &
       &                         , ista_kg1_k, np_kg1_k, mp_kg1_k &
       &                         , nbmx_ext, ierr, ista_atm,iend_atm
  use m_Electronic_Structure, only : zaj_l,  fsr_l, fsi_l , fsr_add_l, fsi_add_l &
       &                           , eko_l, vnlph_l, vlhxcQ, compr_l, compi_l    &
       &                           , fsr_l_2D_k, fsi_l_2D_k, neordr, nrvf_ordr

  use m_Control_Parameters,   only : nblocksize_vnonlocal_w_nlmta_is_given            &
 &                                 , nblocksize_vnonlocal_w_nlmta                     &
 &                                 , nblocksize_vnonlocal_w_f_is_given                &
 &                                 , nblocksize_vnonlocal_w_f
  use m_FFT,                  only : m_FFT_W_Vlocal_W_3D
  use m_Parallelization,      only : lrank, nbsn, nbsn_sta, nbsn_end, neg_g           &
 &                                 , nbs_num, nbsn_num, nbs_sta, nbs_end              &
 &                                 , mpi_kg_world, mpi_ke_world           &
 &                                 , myrank_g, nrank_g &
 &                                 , ista_g1k, iend_g1k, np_g1k, mp_g1k   &
 &                                 , nis_g1k, nie_g1k, nel_g1k               &
 &                                 , nis_e, nel_fs, nis_fs, np_fs        &
 &                                 , ista_fs, iend_fs     
  use m_PlaneWaveBasisSet,    only : m_pwBS_kinetic_energies_3D
  use m_Electronic_Structure, only : fsr_l_2D, fsi_l_2D         &
 &                                 , fsr_ball, fsi_ball         &
 &                                 , m_ES_gather_f_3d_to_2d_blk &
 &                                 , m_ES_gather_f_3d_to_2d     &
 &                                 , m_ES_alloc_fsr_l_2d        &
 &                                 , m_ES_alloc_fsi_l_2d        &
 &                                 , m_ES_dealloc_fsr_l_2d      &
 &                                 , m_ES_dealloc_fsi_l_2d
  use m_PlaneWaveBasisSet,   only : ngabc_kngp_l, ngabc_kngp_B_l
  use m_Electronic_Structure,only :  afft,bfft 

  use m_Electronic_Structure, only : m_ES_WF_in_Rspace_3D, m_ES_WF_2D_psi

!<<----

! ============================ added by K. Tagami =================== 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor, SpinOrbit_mode, sw_hubbard, &
       &                             orb_popu_method, sw_band_unfolding, &
       &                             band_unfolding_active
  use m_PseudoPotential,      only : dion_scr_noncl, q_noncl, taup_pao
  use m_ES_NonCollinear,      only : m_ES_alloc_spinor_eigenwfn_0, &
       &                             m_ES_dealloc_spinor_eigenwfn_0, &
       &                             m_ES_set_spinor_eigenwfn_0, &
       &                             Spinor_EigenWfn0_atomtyp
  use m_SpinOrbit_Potential,  only : EigenWfns_MatLS_L0, &
       &                            EigenWfns_MatLS_L1, &
       &                            EigenWfns_MatLS_L2, &
       &                            EigenWfns_MatLS_L3
  use m_SpinOrbit_Potential,   only : EigenVals_MatLS_L0, &
       &                            EigenVals_MatLS_L1, &
       &                            EigenVals_MatLS_L2, &
       &                            EigenVals_MatLS_L3
  use m_Const_Parameters,    only : Neglected, BuiltIn, ByProjector, ByPawPot, &
       &                            ZeffApprox, ReadFromPP
! ================================================================== 11.0

  use m_PlaneWaveBasisSet,   only :  GVec_on_refcell
  use mpi

  implicit none

! === DEBUG by tkato for restart!!! 2012/02/12 =================================
! integer                                             :: npesize
  integer                                             :: npesize = 0
! ==============================================================================

  real(kind=DP),        allocatable,dimension(:)      :: zfcos, zfsin  !d(nbmx_ext)
  real(kind=DP),        allocatable,dimension(:)      :: zfcos_mpi,zfsin_mpi !d(mp_nbmx)
!!$  real(kind=DP),        allocatable,dimension(:)      :: zfcos_mpi_k, zfsin_mpi_k !d(mp_nbmx_k)
  real(kind=DP),private,allocatable,dimension(:)      :: ar, ai

  real(kind=DP), allocatable, target, dimension(:)    :: sc, ss, qc, qs
  real(kind=DP), allocatable, target, dimension(:)    :: sc_l

  integer, private, parameter                         :: sw_timing_2ndlevel = ON

#ifdef NONLOCAL_DGEMM
  real(kind=DP),        allocatable,dimension(:,:)    :: wk_zfcos, wk_zfsin  !d(nbmx_ext)
!!$  integer      ,allocatable,dimension(:)              :: lmtt_tmp
  real(kind=DP),private,allocatable,dimension(:,:)    :: wk_ar,  wk_ai
  real(kind=DP),allocatable,dimension(:,:)            :: bp_tmp1,  bp_tmp2
!!$  real(kind=DP),allocatable,dimension(:)              :: ia_tmp
  logical :: DGEMM_DEBUG = .false.
#endif

  real(kind=DP),          allocatable, dimension(:,:) :: pre_sc_without, pre_ss_without

#ifdef SX
  integer,private                                     :: nb_vnonlocal_default = 5000  ! TY 26Aug2009
  integer,private                                     :: nb_betar_default     = 10000 ! TY
#else
  integer,private                                     :: nb_vnonlocal_default = 1000  ! TY 26Aug2009
  integer,private                                     :: nb_betar_default     = 32    ! TY
#endif

  integer, allocatable, dimension(:)                  :: full_to_dist,dist_to_full

  logical, allocatable, dimension(:)                  :: mapl
  integer                                             :: nonzero_fft_elements
  integer, allocatable, dimension(:)                  :: revmap

  integer,private                                     :: nlmta_us
!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  contains


  subroutine m_ES_alloc_scss_etc_3D()
    allocate(sc_l(1:maxval(np_g1k(:))))
  end subroutine m_ES_alloc_scss_etc_3D

  subroutine m_ES_dealloc_scss_etc()
    if(allocated(qs)) deallocate(qs)
    if(allocated(qc)) deallocate(qc)
    if(allocated(ss)) deallocate(ss)
    if(allocated(sc)) deallocate(sc)
    if(allocated(sc_l)) deallocate(sc_l)
    call m_FFT_dealloc_WF_work()
  end subroutine m_ES_dealloc_scss_etc


  subroutine m_ES_alloc_afft_scss_etc_3D()
! === DEBUG by tkato 2011/07/12 ================================================
!   allocate(afft(nfft))
!   allocate(bfft(nfft))
! Why are these lines comments???
    allocate(afft(nfft))
    allocate(bfft(nfft))
! ==============================================================================
    call m_ES_alloc_scss_etc_3D()
  end subroutine m_ES_alloc_afft_scss_etc_3D

  subroutine m_ES_dealloc_afft_scss_etc()
    call m_ES_dealloc_scss_etc()
    if(allocated(bfft)) deallocate(bfft)
    if(allocated(afft)) deallocate(afft)
  end subroutine m_ES_dealloc_afft_scss_etc

  subroutine alloc_zfsincos_mpi()
    if(flag_mpi_g_dot_r) then
       allocate(zfcos_mpi(mp_nbmx))
       allocate(zfsin_mpi(mp_nbmx))
    end if
  end subroutine alloc_zfsincos_mpi

  subroutine dealloc_zfsincos_mpi()
    if(flag_mpi_g_dot_r) deallocate(zfsin_mpi,zfcos_mpi)
  end subroutine dealloc_zfsincos_mpi

!!$  subroutine alloc_zfsincos_mpi_k()
!!$    if(flag_mpi_g_dot_r_k) then
!!$       allocate(zfcos_mpi_k(mp_nbmx_k))
!!$       allocate(zfsin_mpi_k(mp_nbmx_k))
!!$    end if
!!$  end subroutine alloc_zfsincos_mpi_k

!!$  subroutine dealloc_zfsincos_mpi_k
!!$    if(flag_mpi_g_dot_r_k) then
!!$       deallocate(zfcos_mpi_k)
!!$       deallocate(zfsin_mpi_k)
!!$    end if
!!$  end subroutine dealloc_zfsincos_mpi_k

  subroutine alloc_zfsincos(ibsize)
    integer, intent(in) :: ibsize
    allocate(zfsin(ibsize))
    allocate(zfcos(ibsize))
  end subroutine alloc_zfsincos

  subroutine dealloc_zfsincos
    deallocate(zfcos,zfsin)
  end subroutine dealloc_zfsincos


  subroutine dealloc_arai
    deallocate(ai,ar)
  end subroutine dealloc_arai

  subroutine m_ES_alloc_zfsincos(ik)
    integer, intent(in) :: ik
    if(ik == 0) then
! ============================== modified by K. Tagami =============== 11.0
!       if(kv3/nspin == 1) then
!          allocate(zfcos(iba(1)))
!          allocate(zfsin(iba(1)))
!       else
!          allocate(zfcos(nbmx_ext))
!          allocate(zfsin(nbmx_ext))
!       end if

       if ( noncol ) then
         if ( kv3/ndim_spinor == 1 ) then
            allocate(zfcos(iba(1)));  allocate(zfsin(iba(1)))
         else
            allocate(zfcos(nbmx_ext));  allocate(zfsin(nbmx_ext))
         endif
       else
         if ( kv3/nspin == 1 ) then
            allocate(zfcos(iba(1)));  allocate(zfsin(iba(1)))
         else
            allocate(zfcos(nbmx_ext));  allocate(zfsin(nbmx_ext))
         endif
       end if
! ==================================================================== 11.0
    else
       allocate(zfcos(iba(ik)))
       allocate(zfsin(iba(ik)))
    end if
! ================================= added by K. Tagami ============= 11.0
    zfcos = 0.0d0;  zfsin = 0.0d0
! ================================================================== 11.0
  end subroutine m_ES_alloc_zfsincos

  subroutine m_ES_dealloc_zfsincos()
    if(allocated(zfsin)) deallocate(zfsin)
    if(allocated(zfcos)) deallocate(zfcos)
  end subroutine m_ES_dealloc_zfsincos

  subroutine m_ES_alloc_arai_3D(ik)
    integer, intent(in) :: ik
    allocate(ar(maxval(np_g1k)))
    allocate(ai(maxval(np_g1k)))
  end subroutine m_ES_alloc_arai_3D

  subroutine m_ES_dealloc_arai_3D()
    deallocate(ar,ai)
  end subroutine m_ES_dealloc_arai_3D

#ifdef NONLOCAL_DGEMM
  subroutine alloc_wkzfsincos(ibsize)
    integer, intent(in) :: ibsize
    integer             :: ichkalloc
    allocate(wk_zfsin(ibsize,natm), stat=ichkalloc)
    allocate(wk_zfcos(ibsize,natm), stat=ichkalloc)
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos', ibsize, natm
      call phase_error_with_msg(nfout,'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos',__LINE__,__FILE__)
    endif
  end subroutine alloc_wkzfsincos

  subroutine alloc_wkzfsincos_red(ibsize,natmsize)
    integer, intent(in) :: ibsize, natmsize
    integer             :: ichkalloc
    allocate(wk_zfsin(ibsize,natmsize), stat=ichkalloc)
    allocate(wk_zfcos(ibsize,natmsize), stat=ichkalloc)
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos_red', ibsize, natm
      call phase_error_with_msg(nfout,'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos_red',__LINE__,__FILE__)
    endif
    wk_zfsin = 0.d0; wk_zfcos = 0.d0
  end subroutine alloc_wkzfsincos_red

  subroutine dealloc_wkzfsincos
    deallocate(wk_zfcos,wk_zfsin)
  end subroutine dealloc_wkzfsincos

  subroutine alloc_wkarai(ibsize,ibsize2)
    integer, intent(in) :: ibsize, ibsize2
    integer             :: ichkalloc
    allocate( wk_ar(ibsize,ibsize2),stat=ichkalloc )
    allocate( wk_ai(ibsize,ibsize2),stat=ichkalloc )
    if( ichkalloc /= 0 ) then
      write(nfout,*) ' could not allocate wk_ar or wk_ai in alloc_wkarai ', ibsize, ibsize2
      call phase_error_with_msg(nfout,' could not allocate wk_ar or wk_ai in alloc_wkarai ',__LINE__,__FILE__)
    endif
    wk_ar=0.d0
    wk_ai=0.d0
  end subroutine alloc_wkarai

  subroutine dealloc_wkarai
    deallocate( wk_ai, wk_ar)
  end subroutine dealloc_wkarai

  subroutine alloc_wkbp(ik,LD11,LD12)
    integer , intent(in) :: ik,LD11,LD12
    integer             :: ichkalloc
    if(k_symmetry(ik) == GAMMA) then
       allocate( bp_tmp1(LD11,LD12)  ,stat=ichkalloc )
    else
       allocate( bp_tmp1(LD11,LD12)  ,stat=ichkalloc )
       allocate( bp_tmp2(LD11,LD12)  ,stat=ichkalloc )
    end if
    if( ichkalloc /= 0 ) then
      write(nfout,*) ' could not allocate bp_tmp in alloc_wkbp ', LD11, LD12
      call phase_error_with_msg(nfout,' could not allocate bp_tmp in alloc_wkbp ',__LINE__,__FILE__)
    endif
  end subroutine alloc_wkbp

  subroutine dealloc_wkbp(ik)
    integer, intent(in) :: ik
    if(k_symmetry(ik) == GAMMA) then
       deallocate( bp_tmp1)
    else
       deallocate( bp_tmp1,bp_tmp2)
    end if
  end subroutine dealloc_wkbp

!!$  subroutine alloc_wkother(icnt)
!!$    integer,intent(in) :: icnt
!!$    integer             :: ichkalloc
!!$    allocate(lmtt_tmp(icnt), stat=ichkalloc )
!!$    allocate(   ia_tmp(icnt), stat=ichkalloc )
!!$    if( ichkalloc /= 0 ) then
!!$      write(nfout,*) 'could not allocate work-array in alloc_wkother', icnt
!!$      stop
!!$    endif
!!$  end subroutine alloc_wkother

!!$  subroutine dealloc_wkother
!!$    deallocate( ia_tmp )
!!$    deallocate( lmtt_tmp)
!!$  end subroutine dealloc_wkother
#endif


#ifdef NONLOCAL_DGEMM

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  subroutine m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part,map)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part
    logical, intent(in), dimension(np_e), optional :: map

    integer :: mdvdb, it, ia, lmt2
    integer :: ibl1,ibl2,ibsize
    integer :: ichkalloc
    integer :: id_sname = -1
!!  real(kind=DP),allocatable,dimension(:,:) :: wk_sc_without,  wk_ss_without
    real(kind=DP),allocatable,dimension(:,:) :: wk_sc_with,     wk_ss_with,     wk_qc_with,     wk_qs_with
    real(kind=DP),allocatable,dimension(:,:) :: fsr_tmp_without, fsi_tmp_without
    real(kind=DP),allocatable,dimension(:,:) :: fsr_tmp_with   , fsi_tmp_with
    integer,allocatable,dimension(:) :: lmt2_tmp_without, ia_tmp_without, lmta2_tmp_without
    integer,allocatable,dimension(:) :: lmt2_tmp_with,    ia_tmp_with,    lmta2_tmp_with
    integer,allocatable,dimension(:) :: il2_tmp_without,  im2_tmp_without
    integer,allocatable,dimension(:) :: il2_tmp_with,     im2_tmp_with
    integer :: icnt_without  , icnt_with , icblk , icbl1 , icbl2
    logical :: mdvdb_EXECUT = .false.
    integer, save :: ibsize_print = OFF
    integer,save :: pre_proc_switch = 999
    integer,save :: pre_proc_ik     = 0
    logical :: chg_cond
!!  integer :: ist
    integer :: lnblck, lnsize, lnb1, lnb2
! === rspace on 3D_Parallel ==============================================================
    logical :: flag_in_realspace

#ifdef __FAPP__
    call fapp_start('Vnonlocal_W',1,1)
#endif
    flag_in_realspace = .false.
    if(sw_rspace == ON) then
       if(m_CtrlP_realspace_integ_OK()) flag_in_realspace = .true.
    end if

    if(ipribetar >= 2 .and. sw_rspace == ON .and. .not. flag_in_realspace) write(nfout,&
       & '(" sw_rspace = ON, but the realspace integration is not installed in the solver now applied. ")')

    if(flag_in_realspace) then
       if(allocated(vnlph_l)) then
          if(size(vnlph_l,1) .ne. np_g1k(ik)) then
             deallocate(vnlph_l)
             allocate(vnlph_l(np_g1k(ik),np_e,kimg))
          end if
       else
          allocate(vnlph_l(np_g1k(ik),np_e,kimg))
       end if
       call Vnonlocal_W_in_realspace(ik,ispin,switch_of_eko_part)
    else
    if(sw_precalculate_phase_vnonlocal==ON) then
      if(present(map)) then
        call m_ES_Vnonlocal_W_precphase(ik,iksnl,ispin,switch_of_eko_part,map)
      else
        call m_ES_Vnonlocal_W_precphase(ik,iksnl,ispin,switch_of_eko_part)
      endif
      return
    endif
! ========================================================================================
    call tstatc0_begin('m_ES_Vnonlocal_W_3D ',id_sname,level=1)
! --> T. Yamasaki, 26th Aug. 2009
!!$    ibsize=5000
                                                                __TIMER_SUB_START(201)
    if(nblocksize_vnonlocal_is_given) then
       ibsize = nblocksize_vnonlocal_w
    else
       ibsize= nb_vnonlocal_default
       if(ipribetar >= 2 .and. ibsize_print == OFF) then
          write(nfout,'(" !ibsize(=nblocksize_vnonlocal_w) (m_ES_Vnonlocal_W) = ",i8)') ibsize
!!$             ibsize_print = ON
       end if
    end if
! <--
    if(nblocksize_vnonlocal_w_nlmta_is_given) then
       icblk = nblocksize_vnonlocal_w_nlmta
    else
       icblk = nlmta
       if(ipribetar >= 2 .and. ibsize_print == OFF) then
          write(nfout,'(" !icblk(=nblocksize_vnonlocal_w_nlmta) (m_ES_Vnonlocal_W) = ",i8)') icblk
       end if
    end if

#ifdef _NONLOCAL_MEMORY_KEEP_
    if ((pre_proc_switch .eq. switch_of_eko_part) .and. (pre_proc_ik .eq. ik)) then
       chg_cond = .false.
    else
       chg_cond = .true.
    endif
#else
    chg_cond = .true.
#endif
    if(nblocksize_vnonlocal_w_f_is_given) then
       lnblck = nblocksize_vnonlocal_w_f
       if(np_e < lnblck) then
          lnblck = np_e
       end if
       if(1 > lnblck) then
          lnblck = np_e
       end if
    else
       lnblck = np_e
       if(ipribetar >= 2 .and. ibsize_print == OFF) then
          write(nfout,'(" !lnblck(=nblocksize_vnonlocal_w_f) (m_ES_Vnonlocal_W) = ",i8)') lnblck
       end if
    end if

    icnt_without=0
    icnt_with=0

    allocate( lmt2_tmp_without(nlmta), stat=ichkalloc ) ;  lmt2_tmp_without=0
    allocate(   ia_tmp_without(nlmta), stat=ichkalloc ) ;    ia_tmp_without=0
    allocate(lmta2_tmp_without(nlmta), stat=ichkalloc ) ; lmta2_tmp_without=0
    allocate( lmt2_tmp_with(nlmta),    stat=ichkalloc ) ;     lmt2_tmp_with=0
    allocate(   ia_tmp_with(nlmta),    stat=ichkalloc ) ;       ia_tmp_with=0
    allocate(lmta2_tmp_with(nlmta),    stat=ichkalloc ) ;    lmta2_tmp_with=0
    allocate( il2_tmp_without(nlmta),  stat=ichkalloc ) ;   il2_tmp_without=0
    allocate( im2_tmp_without(nlmta),  stat=ichkalloc ) ;   im2_tmp_without=0
    allocate( il2_tmp_with(nlmta),     stat=ichkalloc ) ;      il2_tmp_with=0
    allocate( im2_tmp_with(nlmta),     stat=ichkalloc ) ;      im2_tmp_with=0
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate work-array in m_ES_Vnonlocal_W', nlmta
      call phase_error_with_msg(nfout,'could not allocate work-array in m_ES_Vnonlocal_W',__LINE__,__FILE__)
    endif

    call pre_m_ES_Vnonlocal_W( icnt_without, icnt_with )

    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       if(mdvdb == EXECUT) then
         mdvdb_EXECUT = .true.
         exit
       endif
    end do

#ifdef _NONLOCAL_MEMORY_KEEP_
    if (chg_cond) then
       if (allocated(pre_sc_without)) deallocate(pre_sc_without)
       if (allocated(pre_ss_without)) deallocate(pre_ss_without)
       allocate( pre_sc_without(np_g1k(ik),icnt_without), stat=ichkalloc ) ;   pre_sc_without=0.d0
       allocate( pre_ss_without(np_g1k(ik),icnt_without), stat=ichkalloc ) ;   pre_ss_without=0.d0
       if( ichkalloc /= 0 ) then
         write(nfout,*) 'could not allocate work-array1 in m_ES_Vnonlocal_W', ibsize, icnt_without, icnt_with, np_e
         call phase_error_with_msg(nfout,'could not allocate work-array1 in m_ES_Vnonlocal_W',__LINE__,__FILE__)
       endif
    endif
#else
    allocate( pre_sc_without(np_g1k(ik),icnt_without), stat=ichkalloc ) ;   pre_sc_without=0.d0
    allocate( pre_ss_without(np_g1k(ik),icnt_without), stat=ichkalloc ) ;   pre_ss_without=0.d0
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate work-array1 in m_ES_Vnonlocal_W', ibsize, icnt_without, icnt_with, np_e
      call phase_error_with_msg(nfout,'could not allocate work-array1 in m_ES_Vnonlocal_W',__LINE__,__FILE__)
    endif
#endif

    if (ibsize .gt. np_g1k(ik)) then
       ibsize = np_g1k(ik)
    endif
    call alloc_wkzfsincos(ibsize)

!!  allocate( wk_sc_without(ibsize,icnt_without), stat=ichkalloc ) ;   wk_sc_without=0.d0
!!  allocate( wk_ss_without(ibsize,icnt_without), stat=ichkalloc ) ;   wk_ss_without=0.d0
!x  allocate( fsr_tmp_without(np_e,min(icblk,icnt_without)) ,stat=ichkalloc ) ; fsr_tmp_without=0.d0
    allocate( fsr_tmp_without(lnblck,min(icblk,icnt_without)) ,stat=ichkalloc ) ; fsr_tmp_without=0.d0
    allocate( wk_sc_with(ibsize,min(icblk,icnt_with)),       stat=ichkalloc ) ;      wk_sc_with=0.d0
    allocate( wk_ss_with(ibsize,min(icblk,icnt_with)),       stat=ichkalloc ) ;      wk_ss_with=0.d0
!x  allocate( fsr_tmp_with(np_e,min(icblk,icnt_with)),       stat=ichkalloc ) ;    fsr_tmp_with=0.d0
    allocate( fsr_tmp_with(lnblck,min(icblk,icnt_with)),       stat=ichkalloc ) ;    fsr_tmp_with=0.d0
!x  call m_ES_alloc_fsr_l_2d()
    call m_ES_alloc_fsr_l_2d(lnblck,nlmta)

    if( ichkalloc /= 0 ) then
       write(nfout,'(" ibsize, icblk, icnt_without, icnt_with = ", 4i8)') ibsize, icblk, icnt_without, icnt_with
       write(nfout,'(" np_e, nlmta = ", 2i8)') np_e, nlmta
      write(nfout,*) 'could not allocate work-array2 in m_ES_Vnonlocal_W', ibsize, icnt_without, icnt_with, np_e
      call phase_error_with_msg(nfout,' could not allocate work-array2 in m_ES_Vnonlocal_W (2)',__LINE__,__FILE__)
    endif

    if( k_symmetry(ik) /= GAMMA ) then
!x    allocate( fsi_tmp_without(np_e,icnt_without), stat=ichkalloc )
!x    allocate( fsi_tmp_with(np_e,min(icblk,icnt_with))      , stat=ichkalloc )
!x    call m_ES_alloc_fsi_l_2d()
      allocate( fsi_tmp_without(lnblck,icnt_without), stat=ichkalloc )
      allocate( fsi_tmp_with(lnblck,min(icblk,icnt_with))      , stat=ichkalloc )
      call m_ES_alloc_fsi_l_2d(lnblck,nlmta)
      fsi_tmp_without=0.d0
      fsi_tmp_with=0.d0
      if( ichkalloc /= 0 ) then
        write(nfout,*) 'could not allocate work-array3 in m_ES_Vnonlocal_W', np_e,icnt_without,icnt_with
        call phase_error_with_msg(nfout,'could not allocate work-array3 in m_ES_Vnonlocal_W',__LINE__,__FILE__)
      endif
    endif

    if(mdvdb_EXECUT) then
      allocate( wk_qc_with(ibsize,min(icblk,icnt_with)), stat=ichkalloc ) ; wk_qc_with=0.d0
      allocate( wk_qs_with(ibsize,min(icblk,icnt_with)), stat=ichkalloc ) ; wk_qs_with=0.d0
      if( ichkalloc /= 0 ) then
        write(nfout,*) 'could not allocate work-array4 in m_ES_Vnonlocal_W', icnt_without, icnt_with, natm
        call phase_error_with_msg(nfout,'could not allocate work-array4 in m_ES_Vnonlocal_W',__LINE__,__FILE__)
      endif
    endif
!
!
! Revised by T. Yamasaki, September 2008
!
   if(DGEMM_DEBUG) write(nfout,*)' DGEMM_debug Vnonlocal :', &
                   icnt_without,icnt_with, nlmta, ibsize, kimg, k_symmetry(ik), GAMMA
   if(ibsize_print == OFF .and. ipribetar >= 2) then
     write(nfout,*)'<<m_ES_Vnonlocal_W>>'
     write(nfout,'(" mype = ",i3, " ibsize, np_e, icnt_with, icnt_without, iba(ik)/ibsize = " &
    &   ,5i8," <<m_ES_Vnonlocal_W>>")') mype, ibsize, np_e, icnt_with, icnt_without, iba(ik)/ibsize
     write(nfout,*)'myrank_e=',myrank_e,'myrank_e=',myrank_e,'myrank_g=',myrank_g, &
    &             'ik=',ik,'ispin=',ispin,'switch_of_eko_part=',switch_of_eko_part
     write(nfout,*)'np_e=',np_e,'iba(ik)=',iba(ik),'ibsize=',ibsize,'nlmta=',nlmta,&
    &             'iksnl=',iksnl,'kimg=',kimg,' chg_cond=',chg_cond
     write(nfout,*)'ista_g1k(ik)=',ista_g1k(ik),'iend_g1k(ik)=',iend_g1k(ik), &
    &             'np_g1k(ik)=',np_g1k(ik),    'mp_g1k(ik)=',mp_g1k(ik)
     write(nfout,*)'ista_e=',ista_e,'iend_e=',iend_e, 'neg=',neg,' lnblck=',lnblck
     write(nfout,*)'icnt_with=',icnt_with,'icnt_without=',icnt_without,'icblk=',icblk, &
    &             'k_symmetry(ik)=',k_symmetry(ik),'GAMMA=',GAMMA,'iteration=',iteration
     ibsize_print = ON
   end if

!!   call gather_f_3d_to_2d(fsr_l_3D, fsr_l_2D)
!x   call m_ES_gather_f_3d_to_2d(fsr_l, fsr_l_2D, ik)
!x   if( k_symmetry(ik) /= GAMMA ) then
!!      call gather_f_3d_to_2d(fsi_l_3D, fsi_l_2D)
!x      call m_ES_gather_f_3d_to_2d(fsi_l, fsi_l_2D, ik)
!x   endif
     if (allocated(vnlph_l)) then
       if (size(vnlph_l,1) .ne. np_g1k(ik)) then
          deallocate(vnlph_l)   
          allocate(vnlph_l(np_g1k(ik),np_e,kimg))
       endif
     else
        allocate(vnlph_l(np_g1k(ik),np_e,kimg))
     endif

     vnlph_l = 0.0d0

  do lnb1 = 1, np_e, lnblck
     lnb2=min( lnb1+lnblck-1,np_e )
     lnsize = lnb2-lnb1 + 1

     call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, ik, lnblck, lnb1, lnsize)
     if( k_symmetry(ik) /= GAMMA ) then
        call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi_l_2D, ik, lnblck, lnb1, lnsize)
     endif

     do ibl1=ista_g1k(ik),iend_g1k(ik),ibsize
        ibl2=min( ibl1+ibsize-1,iend_g1k(ik) )
        call calc_phase_blk(ik,ibl1,ibl2)
        if( icnt_with .gt. 0 ) then
           do icbl1 = 1 , icnt_with , icblk
              icbl2 = min( icbl1+icblk-1,icnt_with )
              call Vnonlocal_W_part_with_blk_3D( ibl1,ibl2,icnt_with )
              call add_vnlph_l_with_eko_blk_3D(  ibsize,ibl1,ibl2,icnt_with, vnlph_l )
           enddo
        endif
        if( icnt_without .gt. 0 ) then
           do icbl1 = 1 , icnt_without , icblk
              icbl2 = min( icbl1+icblk-1,icnt_without )
              call Vnonlocal_W_part_without_blk_3D( ibl1,ibl2,icnt_without)
              call add_vnlph_l_without_eko_blk_3D(  ibsize,ibl1,ibl2,icnt_without,vnlph_l)
           enddo
        endif
     enddo

   enddo

    deallocate( lmt2_tmp_without,ia_tmp_without,lmta2_tmp_without,il2_tmp_without,im2_tmp_without )
    deallocate( lmt2_tmp_with,   ia_tmp_with,   lmta2_tmp_with,   il2_tmp_with,   im2_tmp_with    )
!!  deallocate( wk_sc_without,   wk_ss_without                   )
    deallocate( wk_sc_with,      wk_ss_with                      )
    deallocate( fsr_tmp_without, fsr_tmp_with )
    call m_ES_dealloc_fsr_l_2d()
    if(mdvdb_EXECUT) then
       deallocate( wk_qc_with,  wk_qs_with )
    endif
    if( k_symmetry(ik) /= GAMMA ) then
       deallocate( fsi_tmp_without, fsi_tmp_with )
       call m_ES_dealloc_fsi_l_2d()
    endif
    call dealloc_wkzfsincos()
#ifndef _NONLOCAL_MEMORY_KEEP_
    deallocate(pre_sc_without)
    deallocate(pre_ss_without)
#endif
    pre_proc_switch = switch_of_eko_part
    pre_proc_ik     = ik

    call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(201)
! === rspace on 3D_Parallel ==============================================================
    end if
! ========================================================================================
#ifdef __FAPP__
   call fapp_stop('Vnonlocal_W',1,1)
#endif
  contains
    subroutine pre_m_ES_Vnonlocal_W( icnt_without, icnt_with )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer,intent(inout) :: icnt_without, icnt_with
                                                                __TIMER_SUB_START(202)
                                                                __TIMER_DO_START(209)
      do ia = 1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(switch_of_eko_part == OFF) mdvdb= SKIP
         if( mdvdb==SKIP ) then
            do lmt2 = 1, ilmt(it)
               icnt_without = icnt_without + 1
               lmt2_tmp_without(icnt_without)  = lmt2
               ia_tmp_without(icnt_without)    =   ia
               lmta2_tmp_without(icnt_without) = lmta(lmt2,ia)
               il2_tmp_without(icnt_without)   =  ltp(lmt2,it)
               im2_tmp_without(icnt_without)   =  mtp(lmt2,it)
            enddo
         elseif( mdvdb==EXECUT ) then
            do lmt2 = 1, ilmt(it)
               icnt_with = icnt_with + 1
               lmt2_tmp_with(icnt_with)  = lmt2
               ia_tmp_with(icnt_with)    =   ia
               lmta2_tmp_with(icnt_with) = lmta(lmt2,ia)
               il2_tmp_with(icnt_with)   =  ltp(lmt2,it)
               im2_tmp_with(icnt_with)   =  mtp(lmt2,it)
            enddo
         endif
      enddo
                                                                __TIMER_DO_STOP(209)
                                                                __TIMER_SUB_STOP(202)
    end subroutine pre_m_ES_Vnonlocal_W

    subroutine calc_phase_blk(ik,ibl1,ibl2)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in) :: ik, ibl1, ibl2
      integer :: i, nb, ia
      real(kind=DP) :: ph, f1, f2, f3
      integer :: id_sname0 = -1
                                                                __TIMER_SUB_START(204)
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('calc_phase_blk ',id_sname0)

!cdir outerunroll=16
                                                                __TIMER_DO_START(213)
!fj!OCL SERIAL
!OCL NOFLTLD
      do i = ibl1, ibl2
         nb = nbase(i,ik)
!fj!OCL PARALLEL
         do ia=1,natm
            f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
            ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
            wk_zfcos(i-ibl1+1,ia) = dcos(ph)
            wk_zfsin(i-ibl1+1,ia) = dsin(ph)
         enddo
      enddo
                                                                __TIMER_DO_STOP(213)
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname0)
                                                                __TIMER_SUB_STOP(204)
    end subroutine calc_phase_blk

    subroutine Vnonlocal_W_part_without_blk_3D(iblk1,iblk2,icnt_without)
!
! (FUJITSU), December 2009
!

      integer, intent(in) :: iblk1,iblk2, icnt_without
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic, iadd, icadd
      real(kind=DP) :: tmp
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(207)
      call tstatc0_begin('Vnonlocal_W_part_without_blk_3D ',id_sname)

!!      wk_sc_without = 0.d0
!!      wk_ss_without = 0.d0

!!      do ic = 1, icnt_without
                                                                __TIMER_DO_START(226)
!$OMP PARALLEL DO SCHEDULE(RUNTIME) DEFAULT(none)    &
!$OMP             SHARED(icbl1,icbl2,ia_tmp_without,ityp,lmta2_tmp_without,il2_tmp_without,  &
!$OMP                    im2_tmp_without,lmt2_tmp_without,fsr_tmp_without,fsr_l_2D,k_symmetry, &
!$OMP                    ik,fsi_tmp_without,fsi_l_2D,lnsize,ilmt,lmtt,ltp,mtp,dion,dion_paw,vlhxcQ, &
!$OMP                    ispin,iwei,iblk1,iblk2,pre_sc_without,wk_zfcos,snl,iksnl,pre_ss_without, &
!$OMP                    wk_zfsin,ista_g1k,q,chg_cond,lnblck,np_e,lnb1,ipaw)   &
!$OMP             PRIVATE(ic,ia,it,lmta2,il2,im2,lmt2,icadd,lmt1,lmtt1,il1,im1,il11,  &
!$OMP                     mdl,tmp,i,iadd)
!OCL NOFLTLD
        do ic = icbl1, icbl2
           ia   = ia_tmp_without(ic)
           it   = ityp(ia)
           lmta2= lmta2_tmp_without(ic)
           il2  = il2_tmp_without(ic)
           im2  = im2_tmp_without(ic)
           lmt2 = lmt2_tmp_without(ic)
           icadd = ic-icbl1+1
!x         fsr_tmp_without(1:np_e,icadd) = fsr_l_2D(1:np_e,lmta2)
           fsr_tmp_without(1:lnsize,icadd) = fsr_l_2D(1:lnsize,lmta2)
           if(k_symmetry(ik) /= GAMMA) then
!x            fsi_tmp_without(1:np_e,icadd) = fsi_l_2D(1:np_e,lmta2)
              fsi_tmp_without(1:lnsize,icadd) = fsi_l_2D(1:lnsize,lmta2)
           endif

           if (.not. chg_cond) then
              cycle
           endif

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
           if (lnblck /= np_e) then
              if (lnb1 /= 1) then
                 cycle
              end if
           end if
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)
              if(il1 == il2 .and. im1 == im2) then
! === DEBUG by tkato 2011/10/05 ================================================
!                tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 endif
! ==============================================================================
              else
! === DEBUG by tkato 2011/10/05 ================================================
!                tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 end if
! ==============================================================================
              endif
              tmp = tmp * iwei(ia)
              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
              if(mdl == 0 .or. mdl == 2) then
!fj                                                                __TIMER_DO_START(227)
                 do i =  iblk1, iblk2
!!                  wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) + tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
!!                  wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    iadd = i-ista_g1k(ik)+1
                    pre_sc_without(iadd,ic) = pre_sc_without(iadd,ic) + tmp*wk_zfcos(i-iblk1+1,ia)*snl(iadd,lmtt1,iksnl)
                    pre_ss_without(iadd,ic) = pre_ss_without(iadd,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(iadd,lmtt1,iksnl)
                 enddo
!fj                                                                __TIMER_DO_STOP(227)
              else if(mdl == 1 .or. mdl == 3) then
!fj                                                                __TIMER_DO_START(228)
                 do i =  iblk1, iblk2
!!                  wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
!!                  wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) - tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    iadd = i-ista_g1k(ik)+1
                    pre_sc_without(iadd,ic) = pre_sc_without(iadd,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(iadd,lmtt1,iksnl)
                    pre_ss_without(iadd,ic) = pre_ss_without(iadd,ic) - tmp*wk_zfcos(i-iblk1+1,ia)*snl(iadd,lmtt1,iksnl)
                 enddo
!fj                                                                __TIMER_DO_STOP(228)
              endif
           end do
        end do
                                                                __TIMER_DO_STOP(226)

!!$                                                                __TIMER_DO_START(227)
!!$                                                                __TIMER_DO_START(228)
!!$                                                                __TIMER_DO_STOP(228)
!!$                                                                __TIMER_DO_STOP(227)

        call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(207)
    end subroutine Vnonlocal_W_part_without_blk_3D

    subroutine Vnonlocal_W_part_with_blk_3D(iblk1,iblk2,icnt_with)
!
! (FUJITSU), December 2009
!

      integer, intent(in) :: iblk1,iblk2, icnt_with
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic, icadd
      real(kind=DP) :: tmp
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(205)
      call tstatc0_begin('Vnonlocal_W_part_with_blk_3D ',id_sname)

        wk_sc_with = 0.d0
        wk_ss_with = 0.d0
        wk_qc_with = 0.d0
        wk_qs_with = 0.d0

!!      do ic = 1, icnt_with
                                                                __TIMER_DO_START(214)
!#ifdef USE_PROF
!  call timer_sta(214)
!#endif
!$OMP PARALLEL DO SCHEDULE(RUNTIME) DEFAULT(none)    &
!$OMP             SHARED(icbl1,icbl2,ia_tmp_with,ityp,lmta2_tmp_with,il2_tmp_with,  &
!$OMP                    im2_tmp_with,lmt2_tmp_with,fsr_tmp_with,fsr_l_2D,k_symmetry, &
!$OMP                    ik,fsi_tmp_with,fsi_l_2D,lnsize,ilmt,lmtt,ltp,mtp,dion,dion_paw,ipaw,vlhxcQ, &
!$OMP                    ispin,iwei,iblk1,iblk2,wk_sc_with,wk_zfcos,snl,iksnl,wk_ss_with, &
!$OMP                    wk_zfsin,ista_g1k,q,wk_qc_with,wk_qs_with)   &
!$OMP             PRIVATE(ic,ia,it,lmta2,il2,im2,lmt2,icadd,lmt1,lmtt1,il1,im1,il11,  &
!$OMP                     mdl,tmp,i)
!OCL NOFLTLD
        do ic = icbl1, icbl2
           ia   = ia_tmp_with(ic)
           it   = ityp(ia)
           lmta2= lmta2_tmp_with(ic)
           il2  = il2_tmp_with(ic)
           im2  = im2_tmp_with(ic)
           lmt2 = lmt2_tmp_with(ic)
           icadd = ic-icbl1+1
!x         fsr_tmp_with(1:np_e,icadd) = fsr_l_2D(1:np_e,lmta2)
           fsr_tmp_with(1:lnsize,icadd) = fsr_l_2D(1:lnsize,lmta2)
           if(k_symmetry(ik) /= GAMMA) then
!x            fsi_tmp_with(1:np_e,icadd) = fsi_l_2D(1:np_e,lmta2)
              fsi_tmp_with(1:lnsize,icadd) = fsi_l_2D(1:lnsize,lmta2)
           endif

           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)
              if(il1 == il2 .and. im1 == im2) then
! === DEBUG by tkato 2011/10/05 ================================================
!                tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 endif
! ==============================================================================
              else
! === DEBUG by tkato 2011/10/05 ================================================
!                tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 end if
! ==============================================================================
              endif
              tmp = tmp * iwei(ia)
              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
!fj                                                              __TIMER_DO_START(215)
              if(mdl == 0 .or. mdl == 2) then
                 do i =  iblk1, iblk2
                    wk_sc_with(i-iblk1+1,icadd) = wk_sc_with(i-iblk1+1,icadd) &
                     &+ tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,icadd) = wk_ss_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                 enddo
              else if(mdl == 1 .or. mdl == 3) then
                 do i =  iblk1, iblk2
                    wk_sc_with(i-iblk1+1,icadd) = wk_sc_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,icadd) = wk_ss_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                 enddo
              endif
!fj                                                              __TIMER_DO_STOP(215)
!fj                                                              __TIMER_DO_START(216)
              if( il1 == il2 .and. im1 == im2) then
                 tmp = q(lmt1,lmt2,it)*iwei(ia)
                 if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                 if(mdl == 0 .or. mdl == 2) then
                    do i =  iblk1, iblk2
                       wk_qc_with(i-iblk1+1,icadd) = wk_qc_with(i-iblk1+1,icadd) &
                     &+ tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                       wk_qs_with(i-iblk1+1,icadd) = wk_qs_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    enddo
                 else if(mdl == 1 .or. mdl == 3) then
                    do i =  iblk1, iblk2
                       wk_qc_with(i-iblk1+1,icadd) = wk_qc_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfsin(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                       wk_qs_with(i-iblk1+1,icadd) = wk_qs_with(i-iblk1+1,icadd) &
                     &- tmp*wk_zfcos(i-iblk1+1,ia)*snl(i-ista_g1k(ik)+1,lmtt1,iksnl)
                    enddo
                 end if
              end if
!fj                                                              __TIMER_DO_STOP(216)
           end do
        end do
                                                                __TIMER_DO_STOP(214)

!!$                                                                __TIMER_DO_START(215)
!!$                                                                __TIMER_DO_START(216)
!!$                                                                __TIMER_DO_STOP(216)
!!$                                                                __TIMER_DO_STOP(215)
        call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(205)
    end subroutine Vnonlocal_W_part_with_blk_3D

    subroutine add_vnlph_l_with_eko_blk_3D(ibsize,ibl1,ibl2,icnt_with,vnlph)
!
! (FUJITSU), December 2009
!
      integer, intent(in) :: ibsize, ibl1, ibl2, icnt_with
      real(kind=DP), intent(inout), dimension(np_g1k(ik),np_e,kimg) :: vnlph
      integer          :: ic, ib
      integer       :: icsize , ista , icnt_size
      real(kind=DP) :: alpha, beta
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(206)
      call tstatc0_begin('add_vnlph_l_with_eko_blk_3D ',id_sname)

      ista = ibl1-ista_g1k(ik)+1
      icnt_size = icbl2-icbl1+1

      if(kimg == 1) then
         icsize=ibl2-ibl1+1
         alpha= 1.d0;  beta= 1.d0

                                                                __TIMER_DGEMM_START(217)
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik) )
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
        &              wk_sc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik) )

         alpha=-1.d0;  beta= 1.d0
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik) )
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
        &              wk_ss_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik) )
                                                                __TIMER_DGEMM_STOP(217)

!!       do ic = 1, icnt_with
                                                                __TIMER_DO_START(223)
         do ic = icbl1, icbl2
!x          do ib = 1, np_e                                           ! MPI
            do ib = 1, lnsize                                         ! MPI
!x             fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib,ik)
!x             fsi_tmp_with(ib,ic-icbl1+1) = fsi_tmp_with(ib,ic-icbl1+1)*eko_l(ib,ik)
               fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib+lnb1-1,ik)
               fsi_tmp_with(ib,ic-icbl1+1) = fsi_tmp_with(ib,ic-icbl1+1)*eko_l(ib+lnb1-1,ik)
            enddo
         enddo
                                                                __TIMER_DO_STOP(223)

         alpha=-1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(218)
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik) )
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
        &              wk_qc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik) )
         alpha= 1.d0;  beta= 1.d0
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,wk_qs_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik) )
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
        &              wk_qs_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik) )
                                                                __TIMER_DGEMM_STOP(218)
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            icsize=ibl2-ibl1+1
            alpha= 1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(219)
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &             wk_sc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &             wk_ss_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(219)

!!          do ic = 1, icnt_with
                                                                __TIMER_DO_START(224)
!OCL NOFLTLD
            do ic = icbl1, icbl2
!x             do ib = 1, np_e                                           ! MPI
               do ib = 1, lnsize                                         ! MPI
!x                fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib,ik)
                  fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib+lnb1-1,ik)
               enddo
            enddo
                                                                __TIMER_DO_STOP(224)
            alpha=-1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(220)
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qs_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qs_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(220)
         else
            icsize=ibl2-ibl1+1
            alpha= 1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(221)
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_sc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha=-1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_ss_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha= 1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_sc_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_sc_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
            alpha= 1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_ss_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(221)

!!          do ic = 1, icnt_with
                                                                __TIMER_DO_START(225)
            do ic = icbl1, icbl2
!x             do ib = 1, np_e                                           ! MPI
               do ib = 1, lnsize                                         ! MPI
!x                fsi_tmp_with(ib,ic-icbl1+1) = fsi_tmp_with(ib,ic-icbl1+1)*eko_l(ib,ik)
!x                fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib,ik)
                  fsi_tmp_with(ib,ic-icbl1+1) = fsi_tmp_with(ib,ic-icbl1+1)*eko_l(ib+lnb1-1,ik)
                  fsr_tmp_with(ib,ic-icbl1+1) = fsr_tmp_with(ib,ic-icbl1+1)*eko_l(ib+lnb1-1,ik)
               enddo
            enddo
                                                                __TIMER_DO_STOP(225)

            alpha=-1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(222)
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qc_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha= 1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qs_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qs_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha=-1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qc_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qc_with,ibsize,fsi_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
            alpha=-1.d0;  beta= 1.d0
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x         &              alpha,wk_qs_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha, &
           &              wk_qs_with,ibsize,fsr_tmp_with,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(222)
         end if
      end if
      call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(206)
    end subroutine add_vnlph_l_with_eko_blk_3D

    subroutine add_vnlph_l_without_eko_blk_3D(ibsize,ibl1,ibl2,icnt_without,vnlph)
!
! (FUJITSU), December 2009
!
      integer, intent(in) :: ibsize, ibl1, ibl2, icnt_without
      real(kind=DP), intent(inout), dimension(np_g1k(ik),np_e,kimg) :: vnlph
      integer       :: icsize, ista, icnt_size
      real(kind=DP) :: alpha, beta
      integer :: id_sname = -1
      call tstatc0_begin('add_vnlph_l_without_eko_blk_3D ',id_sname)
                                                                __TIMER_SUB_START(208)

      ista = ibl1+1-ista_g1k(ik)
      icnt_size = icbl2-icbl1+1

!!    wk_sc_without(1:ibsize,1:icnt_without) = pre_sc_without(ista:ista+ibsize-1,1:icnt_without)
!!    wk_ss_without(1:ibsize,1:icnt_without) = pre_ss_without(ista:ista+ibsize-1,1:icnt_without)

      if(kimg == 1) then
         icsize=ibl2-ibl1+1
         alpha= 1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(229)
!!       call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,pre_sc_without(ista,icbl1),np_g1k(ik), fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_sc_without(ista,icbl1), &
        &              np_g1k(ik),fsr_tmp_without,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))

         alpha=-1.d0;  beta= 1.d0
!!       call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &              alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
!x       call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &              alpha,pre_ss_without(ista,icbl1),np_g1k(ik), fsi_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
         call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_ss_without(ista,icbl1), &
        &              np_g1k(ik),fsi_tmp_without,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(229)
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            icsize=ibl2-ibl1+1
            alpha= 1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(230)
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_sc_without(ista,icbl1),np_g1k(ik), fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_sc_without(ista,icbl1), &
        &                 np_g1k(ik),fsr_tmp_without,lnblck,beta,vnlph(ista,lnb1,1),np_g1k(ik))
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_ss_without(ista,icbl1),np_g1k(ik), fsr_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_ss_without(ista,icbl1), &
        &                 np_g1k(ik),fsr_tmp_without,lnblck,beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(230)
         else
            icsize=ibl2-ibl1+1
            alpha= 1.d0;  beta= 1.d0
                                                                __TIMER_DGEMM_START(231)
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_sc_without(ista,icbl1),np_g1k(ik), fsr_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_sc_without(ista,icbl1), &
        &                 np_g1k(ik), fsr_tmp_without,lnblck, beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha=-1.d0;  beta= 1.d0
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_ss_without(ista,icbl1),np_g1k(ik), fsi_tmp_without,np_e, beta,vnlph(ista,1,1),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_ss_without(ista,icbl1), &
        &                 np_g1k(ik), fsi_tmp_without,lnblck, beta,vnlph(ista,lnb1,1),np_g1k(ik))
            alpha= 1.d0;  beta= 1.d0
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_sc_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_sc_without(ista,icbl1),np_g1k(ik), fsi_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_sc_without(ista,icbl1), &
        &                 np_g1k(ik), fsi_tmp_without,lnblck, beta,vnlph(ista,lnb1,2),np_g1k(ik))
            alpha= 1.d0;  beta= 1.d0
!!          call DGEMM__('N','T', icsize,np_e,icnt_without, &
!!      &                 alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
!x          call DGEMM__('N','T', icsize,np_e,icnt_size, &
!x      &                 alpha,pre_ss_without(ista,icbl1),np_g1k(ik), fsr_tmp_without,np_e, beta,vnlph(ista,1,2),np_g1k(ik))
            call DGEMM__('N','T', icsize,lnsize,icnt_size,alpha,pre_ss_without(ista,icbl1), &
        &                 np_g1k(ik), fsr_tmp_without,lnblck, beta,vnlph(ista,lnb1,2),np_g1k(ik))
                                                                __TIMER_DGEMM_STOP(231)
         end if
      end if
      call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(208)
    end subroutine add_vnlph_l_without_eko_blk_3D

  end subroutine m_ES_Vnonlocal_W_3D

  subroutine m_ES_AtaulmnaG(hardpart)
    logical, intent(in), optional :: hardpart
    complex(kind=DP), allocatable, dimension(:,:) :: wkexp
    integer :: ia,it,ig,igg,iksnl,ik,lmt1,lmt2,im1,im2,lmtt1,lmtt2,lmta1,lmta2,il1,il2,mil
    integer :: ispin
    complex(kind=DP) :: ctmp,cfac,i1=(0.d0,1.d0)
    real(kind=DP) :: tmp
    integer :: mdvdb
    integer :: icount
    logical :: hp
    integer :: ibsize,ibl,ibl1,ibl2
    integer :: id_sname = -1
    call tstatc0_begin('m_ES_AtaulmnaG ',id_sname,1)
    if(nblocksize_vnonlocal_is_given) then
       ibsize = nblocksize_vnonlocal_w
    else
       ibsize= nb_vnonlocal_default
    end if
    allocate(wkexp(maxval(np_g1k),natm))
    hp = .false.
    if(present(hardpart)) hp = hardpart
    if(modnrm==ON .and. hp) BtaulmaG = 0.d0
    do ispin = 1, nspin, af+1
    do ik=ispin, kv3-nspin+ispin,nspin
      if(map_k(ik) /= myrank_k) cycle
      if ( noncol ) then
         iksnl = (ik-1)/ndim_spinor + 1   
      else
         iksnl = (ik-1)/nspin + 1         
      endif
      call cal_phase()
!      do ibl1=1,np_g1k(ik),ibsize
      icount = 0
!      ibl2=min( ibl1+ibsize-1,np_g1k(ik) )
      AtaulmaG(1:np_g1k(ik),:,ik,:) = 0.d0
      do ia=1,natm
        it = ityp(ia)
        mdvdb = m_PP_include_vanderbilt_pot(it)
        do lmt2 = 1, ilmt(it)
          lmtt2 = lmtt(lmt2,it)
          lmta2 = lmta(lmt2,ia)
          il2   = ltp(lmt2,it)
          im2   = mtp(lmt2,it)
          if(mdvdb==ON) icount = icount + 1
          do lmt1 = 1, ilmt(it)
            lmtt1 = lmtt(lmt1,it)
            lmta1 = lmta(lmt1,ia)
            il1   = ltp(lmt1,it)
            im1   = mtp(lmt1,it)
            mil   = mod(il1-1,4)
            cfac  = i1**(-mil)
            if(il1 == il2 .and. im1 == im2) then
!!$               tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
               if(ipaw(it)==0) then
                  tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
               else
                  tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
               endif
            else
!!$               tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
               if(ipaw(it)==0) then
                  tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
               else
                  tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
               end if
            endif
            tmp = tmp * iwei(ia)
!            do ig=ibl1,ibl2
            do ig=1,np_g1k(ik)
              ctmp  = cfac*wkexp(ig,ia)*snl(ig,lmtt1,iksnl)
              AtaulmaG(ig,lmta2,ik,1) = AtaulmaG(ig,lmta2,ik,1) + tmp *  real(ctmp)
              AtaulmaG(ig,lmta2,ik,2) = AtaulmaG(ig,lmta2,ik,2) + tmp * dimag(ctmp)
            enddo
            if(hp .and. mdvdb==ON .and. il1 == il2 .and. im1 == im2) then
              tmp = q(lmt1,lmt2,it)*iwei(ia)
              !do ig=ibl1,ibl2
              do ig=1,np_g1k(ik)
                ctmp  = cfac*wkexp(ig,ia)*snl(ig,lmtt1,iksnl)
                BtaulmaG(ig,icount,ik,1) = BtaulmaG(ig,icount,ik,1) + tmp *  real(ctmp)
                BtaulmaG(ig,icount,ik,2) = BtaulmaG(ig,icount,ik,2) + tmp * dimag(ctmp)
              enddo
            endif
          enddo
        enddo
      !enddo
      enddo
      nlmta_us = icount
    enddo
    enddo
    deallocate(wkexp)
    call tstatc0_end(id_sname)

    contains

    subroutine cal_phase()
      real(kind=DP) :: f1,f2,f3,ph
      integer :: nb
      do ia=1, natm
        f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
        do ig=ista_g1k(ik),iend_g1k(ik)
          igg = ig-ista_g1k(ik)+1
          nb = nbase(ig,ik)
          ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
          wkexp(igg,ia) = dcmplx(dcos(ph),-dsin(ph))
        enddo
      enddo
    end subroutine cal_phase

  end subroutine m_ES_AtaulmnaG

  subroutine m_ES_Vnonlocal_W_precphase(ik,iksnl,ispin,switch_of_eko_part,map)
    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part
    logical, intent(in), dimension(np_e), optional :: map
    integer :: i, ik_for_pointing_eko, ib, icount, ia, it, mdvdb, lmt2, lmta2
    integer :: ibsize,ibl,ibl1,ibl2,mg,ndata,icountb
    logical :: umap
    real(kind=DP) :: e
    real(kind=DP), allocatable, dimension(:,:) :: efsr_l, efsi_l, fsr, fsi
    real(kind=DP), allocatable, dimension(:,:) :: fsrtmp,fsitmp,vnltmpr,vnltmpi
    integer :: id_sname = -1
    if(nblocksize_vnonlocal_is_given) then
       ibsize = nblocksize_vnonlocal_w
    else
       ibsize= nb_vnonlocal_default
    end if
    if (allocated(vnlph_l)) then
      if (size(vnlph_l,1) .ne. np_g1k(ik)) then
         deallocate(vnlph_l)   
         allocate(vnlph_l(np_g1k(ik),np_e,kimg))
      endif
    else
       allocate(vnlph_l(np_g1k(ik),np_e,kimg))
    endif
    call tstatc0_begin('m_ES_Vnonlocal_W_precphase ',id_sname,1)
    allocate(fsr(np_e,nlmta))
    if(k_symmetry(ik)/=GAMMA) allocate(fsi(np_e,nlmta))
    call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr, ik, np_e, 1, np_e)
    if( k_symmetry(ik)/=GAMMA ) then
       call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi, ik, np_e, 1, np_e)
    endif

    mg = np_g1k(ik)
    umap = .false.
    if(present(map)) then
      umap = .true.
    endif
    if(umap) then
      ndata = 0
      do i=1,np_e 
        if(.not. map(i)) ndata = ndata+1
      enddo
    endif
    if(ndata == np_e) umap = .false.
    if(umap) then
      allocate(fsrtmp(ndata,nlmta))
      allocate(vnltmpr(mg,ndata))
      if(kimg==2) allocate(vnltmpi(mg,ndata))
      if(k_symmetry(ik) /= GAMMA) allocate(fsitmp(ndata,nlmta))
      icount = 0
      do i=1,np_e
        if(.not.map(i)) then
          icount = icount+1
          fsrtmp(icount,:) = fsr(i,:)
          if(k_symmetry(ik) /= GAMMA) fsitmp(icount,:) = fsi(i,:)
        endif
      enddo
    endif
    if(.not.umap) then
      if(kimg==1) then
        call DGEMM__('N','T',mg,np_e,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsr,np_e,0.d0,vnlph_l(1,1,1),mg)
        call DGEMM__('N','T',mg,np_e,nlmta,-1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsi,np_e,1.d0,vnlph_l(1,1,1),mg)
      else
        call DGEMM__('N','T',mg,np_e,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsr,np_e,0.d0,vnlph_l(1,1,1),mg)
        if(k_symmetry(ik)/=GAMMA) &
        call DGEMM__('N','T',mg,np_e,nlmta,-1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsi,np_e,1.d0,vnlph_l(1,1,1),mg)
        call DGEMM__('N','T',mg,np_e,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsr,np_e,0.d0,vnlph_l(1,1,2),mg)
        if(k_symmetry(ik)/=GAMMA) &
        call DGEMM__('N','T',mg,np_e,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsi,np_e,1.d0,vnlph_l(1,1,2),mg)
      endif
    else
      if(kimg==1) then
        call DGEMM__('N','T',mg,ndata,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsrtmp,ndata,0.d0,vnltmpr,mg)
        call DGEMM__('N','T',mg,ndata,nlmta,-1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsitmp,ndata,1.d0,vnltmpr,mg)
      else
        call DGEMM__('N','T',mg,ndata,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsrtmp,ndata,0.d0,vnltmpr,mg)
        if(k_symmetry(ik)/=GAMMA) &
        call DGEMM__('N','T',mg,ndata,nlmta,-1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsitmp,ndata,1.d0,vnltmpr,mg)
        call DGEMM__('N','T',mg,ndata,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,2),mg,fsrtmp,ndata,0.d0,vnltmpi,mg)
        if(k_symmetry(ik)/=GAMMA) &
        call DGEMM__('N','T',mg,ndata,nlmta, 1.d0,AtaulmaG(1:mg,:,ik,1),mg,fsitmp,ndata,1.d0,vnltmpi,mg)
      endif
    endif
    if(modnrm==ON .and. switch_of_eko_part==ON) then
      if ( noncol ) then
         ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1
      else
         ik_for_pointing_eko = ik
      endif
      icount = 0
      if(.not.umap) then
         allocate(efsr_l(np_e,nlmta_us))
         if(k_symmetry(ik) /= GAMMA) then
           allocate(efsi_l(np_e,nlmta_us))
         endif
         do ia=1,natm
           it = ityp(ia)
           mdvdb = m_PP_include_vanderbilt_pot(it)
           do lmt2 = 1, ilmt(it)
             if(mdvdb==ON) then
               icount = icount+1
              lmta2 = lmta(lmt2,ia)
              do ib=1,np_e
                efsr_l(ib,icount) = -eko_l(ib,ik_for_pointing_eko) * fsr(ib,lmta2)
                if(k_symmetry(ik) /=  GAMMA) &
                efsi_l(ib,icount) = -eko_l(ib,ik_for_pointing_eko) * fsi(ib,lmta2)
              enddo
            endif
          enddo
        enddo

        if(kimg==1) then
          call DGEMM__('N','T',mg,np_e,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsr_l,np_e,1.d0,vnlph_l(1,1,1),mg)
          call DGEMM__('N','T',mg,np_e,nlmta_us,-1.d0,BtaulmaG(1:mg,:,ik,2),mg,efsi_l,np_e,1.d0,vnlph_l(1,1,1),mg)
        else
          call DGEMM__('N','T',mg,np_e,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsr_l,np_e,1.d0,vnlph_l(1,1,1),mg)
          if(k_symmetry(ik) /= GAMMA) &
          call DGEMM__('N','T',mg,np_e,nlmta_us,-1.d0,BtaulmaG(1:mg,:,ik,2),mg,efsi_l,np_e,1.d0,vnlph_l(1,1,1),mg)
          call DGEMM__('N','T',mg,np_e,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,2),mg,efsr_l,np_e,1.d0,vnlph_l(1,1,2),mg)
          if(k_symmetry(ik) /= GAMMA) &
          call DGEMM__('N','T',mg,np_e,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsi_l,np_e,1.d0,vnlph_l(1,1,2),mg)
        endif
        deallocate(efsr_l)
        if(k_symmetry(ik) /= GAMMA) then
          deallocate(efsi_l)
        endif
      else
         allocate(efsr_l(ndata,nlmta_us))
         if(k_symmetry(ik) /= GAMMA) then
           allocate(efsi_l(ndata,nlmta_us))
         endif
         do ia=1,natm
           it = ityp(ia)
           mdvdb = m_PP_include_vanderbilt_pot(it)
           do lmt2 = 1, ilmt(it)
             if(mdvdb==ON) then
               icount = icount+1
              lmta2 = lmta(lmt2,ia)
              icountb = 0
              do ib=1,np_e
                if(.not.map(ib)) then
                  icountb = icountb+1
                  efsr_l(icountb,icount) = -eko_l(ib,ik_for_pointing_eko) * fsr(ib,lmta2)
                  if(k_symmetry(ik) /=  GAMMA) &
                  efsi_l(icountb,icount) = -eko_l(ib,ik_for_pointing_eko) * fsi(ib,lmta2)
                endif
              enddo
            endif
          enddo
        enddo

        if(kimg==1) then
          call DGEMM__('N','T',mg,ndata,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsr_l,ndata,1.d0,vnltmpr,mg)
          call DGEMM__('N','T',mg,ndata,nlmta_us,-1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsi_l,ndata,1.d0,vnltmpr,mg)
        else
          call DGEMM__('N','T',mg,ndata,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsr_l,ndata,1.d0,vnltmpr,mg)
          if(k_symmetry(ik) /= GAMMA) &
          call DGEMM__('N','T',mg,ndata,nlmta_us,-1.d0,BtaulmaG(1:mg,:,ik,2),mg,efsi_l,ndata,1.d0,vnltmpr,mg)
          call DGEMM__('N','T',mg,ndata,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,2),mg,efsr_l,ndata,1.d0,vnltmpi,mg)
          if(k_symmetry(ik) /= GAMMA) &
          call DGEMM__('N','T',mg,ndata,nlmta_us, 1.d0,BtaulmaG(1:mg,:,ik,1),mg,efsi_l,ndata,1.d0,vnltmpi,mg)
        endif
        deallocate(efsr_l)
        if(k_symmetry(ik) /= GAMMA) then
          deallocate(efsi_l)
        endif
      endif
    endif
    if(umap) then
      icountb = 0
      do ib=1,np_e
        if(.not.map(ib)) then
          icountb = icountb+1
          vnlph_l(1:mg,ib,1) = vnltmpr(1:mg,icountb)
          if(kimg==2) vnlph_l(1:mg,ib,2) = vnltmpi(1:mg,icountb)
        endif
      enddo
    endif
    deallocate(fsr)
    if(k_symmetry(ik) /= GAMMA) then
      deallocate(fsi)
    endif
    if(umap) then
      deallocate(fsrtmp)
      if(k_symmetry(ik)/=GAMMA) deallocate(fsitmp)
      deallocate(vnltmpr)
      if(kimg==2) deallocate(vnltmpi)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_ES_Vnonlocal_W_precphase
!!! else NONLOCAL_DGEMM
#else

  subroutine m_ES_Vnonlocal_W_3D(ik,iksnl,ispin,switch_of_eko_part)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!

    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part

    print *,"ERROR : Not Support Multi Dimension Para on DGEM"
    MPI_ABORT(mpi_comm_world, 140, ik)

  end subroutine m_ES_Vnonlocal_W_3D
#endif
!! #ifdef NONLOCAL_DGEMM end



  subroutine m_ES_wd_fsr_fsi()
    integer :: ik
    if(ipribetar >= 2) then
       write(nfout,'(" --- fsr_l, fsi_l ---")')
       do ik = ista_k, iend_k                              ! MPI
          call wd_fsr_fsi_3D(ik)    ! MPI
       end do
    end if
  end subroutine m_ES_wd_fsr_fsi


  subroutine G_dot_R(ia,mapmode,ik)
    integer, intent(in) :: ia,mapmode
    integer, intent(in), optional :: ik
    integer :: i, i1, ip
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1

    call tstatc0_begin('G_dot_R ',id_sname)
    f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
    if(.not.allocated(nbase)) call phase_error_with_msg(nfout,' nbase is not allocated',__LINE__,__FILE__)
    if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
       do i = 1, np_g1k(ik)
          i1 = nbase(i+ista_g1k(ik)-1,ik)
          grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    else
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
       do i = 1, np_nbmx
          ip = i + ista_nbmx -1
          grt = ngabc(ip,1)*f1+ngabc(ip,2)*f2+ngabc(ip,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine G_dot_R

  subroutine G_dot_R_map_div(ia,ik,ibl1,ibl2)
    integer, intent(in) :: ia,ik,ibl1,ibl2
    integer :: i, i1
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1
    call tstatc0_begin('G_dot_R_map_div ',id_sname)

    f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
    if(.not.allocated(nbase)) call phase_error_with_msg(nfout,' nbase is not allocated',__LINE__,__FILE__)

    if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       zfcos = 0.0D0;   zfsin = 0.0D0
       do i = 1, ibl2-ibl1+1
          i1 = nbase(i+ibl1-1,ik)
          if ( sw_force_kpt_inside_bz == ON ) then
             if ( GVec_on_refcell(i1,ik) == 0 ) cycle
          else
             if ( GVec_on_refcell(i1,1) == 0 ) cycle
          endif
          grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    else
!!$    do i = 1, iba(ik)
       do i = 1, ibl2-ibl1+1
          i1 = nbase(i+ibl1-1,ik)
          grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    endif
    call tstatc0_end(id_sname)
  end subroutine G_dot_R_map_div

#ifdef NONLOCAL_DGEMM
  subroutine G_dot_R_map_blk_3D(ik,ibl1,ibl2,ia1,ia2)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!

    integer, intent(in) :: ik,ibl1,ibl2, ia1,ia2
    integer :: ia, i, i1
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1
                                                                __TIMER_SUB_START(410)

    if(.not.allocated(nbase)) call phase_error_with_msg(nfout,' nbase is not allocated',__LINE__,__FILE__)
    call tstatc0_begin('G_dot_R_map_blk ',id_sname)

   if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       wk_zfcos = 0.0D0;     wk_zfsin = 0.0D0
       do i = ibl1, ibl2
          i1 = nbase(i+ista_g1k(ik)-1,ik)
          if ( sw_force_kpt_inside_bz == ON ) then
             if ( GVec_on_refcell(i1,ik) == 0 ) cycle
          else
             if ( GVec_on_refcell(i1,1) == 0 ) cycle
          endif
          do ia=ia1,ia2
             f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
             grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
             wk_zfcos(i-ibl1+1,ia-ia1+1) = dcos(grt)
             wk_zfsin(i-ibl1+1,ia-ia1+1) = dsin(grt)
          end do
       end do
    else
!cdir outerunroll=16
       do i = ibl1, ibl2
          !f       i1 = nbase(i,ik)
          i1 = nbase(i+ista_g1k(ik)-1,ik)
!!$       do ia=1,natm
          do ia=ia1,ia2
             f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
             grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
             wk_zfcos(i-ibl1+1,ia-ia1+1) = dcos(grt)
             wk_zfsin(i-ibl1+1,ia-ia1+1) = dsin(grt)
          end do
       end do
    endif
    call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(410)
  end subroutine G_dot_R_map_blk_3D
#endif
!! #ifdef NONLOCAL_DGEMM end



!kukan 4
  subroutine m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik,dealloc_fs2d)
    integer, intent(in) :: nfout,ik
    logical, intent(in), optional :: dealloc_fs2d
    integer :: mod_ball = ON
    call m_ES_betar_dot_Psi_4_each_k_3D(nfout,zaj_l,ista_k,iend_k,ik,fsr_l,fsi_l,mod_ball,dealloc_fs2d)
  end subroutine m_ES_betar_dot_WFs_4_each_k_3D

#ifdef NONLOCAL_DGEMM
  subroutine m_ES_betar_dot_Psi_4_each_k_3D(nfout,psi_l,k1,k2,ik,bpr_l,bpi_l,mod_ball,dealloc_fs2d,map)
    integer, intent(in)    :: nfout, k1,k2,ik
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg):: psi_l ! MPI
    real(kind=DP),intent(out),dimension(np_e,np_fs,k1:k2) :: bpr_l, bpi_l !MPI
    integer, intent(in), optional :: mod_ball
    logical, intent(in), optional :: dealloc_fs2d
    logical, intent(in), dimension(np_e), optional :: map
    if(present(map)) then
      call m_ES_betar_dot_Psi_k_3D_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l,mod_ball,dealloc_fs2d,map=map)
    else
      call m_ES_betar_dot_Psi_k_3D_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l,mod_ball,dealloc_fs2d)
    endif
  end subroutine m_ES_betar_dot_Psi_4_each_k_3D

  subroutine m_ES_betar_dot_Psi_k_3D_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l,mod_ball,dealloc_fs2d,precalculate_phase,map)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, October 2009 : Mblock
!
    use m_Parallelization,only      : ball_buff, ball_addr
! === DEBUG by tkato 2014/08/08 ==========================================================
    use m_Electronic_Structure, only : m_ES_gather_f_3d_to_2d
! ========================================================================================

    integer, intent(in)    :: nfout, k1,k2,ik
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt,ista_snl:iend_snl) :: snl
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg):: psi_l ! MPI
    real(kind=DP),intent(out),dimension(np_e,np_fs,k1:k2) :: bpr_l, bpi_l !MPI
    logical, intent(in), optional :: precalculate_phase
    integer, intent(in), optional :: mod_ball
    logical, intent(in), optional :: dealloc_fs2d
    logical, intent(in), dimension(np_e), optional :: map
    logical :: deallocfs2d

    integer :: iksnl, ibsize, iend
    integer :: datasize

    integer :: ia, it, lmt1, msize, msize_target, msizemax, msizesum, natm_redmax
    integer :: ibl1,ibl2, ia1, ia2
    integer :: LD11,  LD12
    logical :: tran1
!!$    integer, allocatable, dimension(:) :: mil
    real(kind=DP), allocatable, dimension(:,:,:) :: psi_ri
!f
!   real(kind=DP),allocatable, dimension(:,:) :: fsr_l_2D,fsi_l_2D
    real(kind=DP),allocatable, dimension(:,:) :: wk_fsr_l,wk_fsi_l
    real(kind=DP),allocatable, dimension(:,:) :: wk_fsr_ball,wk_fsi_ball
    real(kind=DP),allocatable, dimension(:,:) :: wk_fsr_2D,wk_fsi_2D
    real(kind=DP),allocatable, dimension(:,:) :: fsrt,fsit
    integer :: i, j, k, n
    logical :: prec

    integer :: id_sname = -1
! === rspace on 3D_Parallel ==============================================================
    logical :: flag_in_realspace


    prec = sw_betar_dot_wfs_exp==ON
    if(present(precalculate_phase)) then
      prec = precalculate_phase
    endif

    flag_in_realspace = .false.
    if(sw_rspace == ON) then
       if(m_CtrlP_realspace_integ_OK()) flag_in_realspace = .true.
    end if

    if(flag_in_realspace) then
       allocate(fsrt(np_e,nlmta));fsrt=0.d0
!       if(k_symmetry(ik) /= GAMMA)then
          allocate(fsit(np_e,nlmta));fsit=0.d0
!       endif
       call betar_dot_Psi_4_each_k_in_rs(nfout,k1,k2,ik,psi_l,fsrt,fsit)
       do i=ista_fs,iend_fs
          bpr_l(1:np_e,i-ista_fs+1,ik) = fsrt(1:np_e,i)
          if(k_symmetry(ik) /= Gamma) then
             bpi_l(1:np_e,i-ista_fs+1,ik) = fsit(1:np_e,i)
          endif
       enddo
       deallocate(fsrt)
       !if(k_symmetry(ik) /= GAMMA) deallocate(fsit)
       deallocate(fsit)
       if(present(mod_ball)) then
       if(mod_ball == ON) then
          allocate(wk_fsr_l(np_fs,np_e))
          allocate(wk_fsr_ball(np_fs,neg))
          if(k_symmetry(ik) /= GAMMA) then
             allocate(wk_fsi_l(np_fs,np_e))
             allocate(wk_fsi_ball(np_fs,neg))
          endif
          do i=1,np_e
            do j=1,np_fs
              wk_fsr_l(j,i)=bpr_l(i,j,ik)
              if(k_symmetry(ik) /= GAMMA) wk_fsi_l(j,i)=bpi_l(i,j,ik)
            enddo
          enddo
          call mpi_allgatherv(wk_fsr_l,np_fs*np_e,mpi_double_precision,wk_fsr_ball(1,1),ball_buff,&
          & ball_addr,mpi_double_precision,mpi_kg_world,ierr)
          if(k_symmetry(ik) /= GAMMA)call mpi_allgatherv(wk_fsi_l,np_fs*np_e,mpi_double_precision,&
                                 &wk_fsi_ball,ball_buff,ball_addr,mpi_double_precision,mpi_kg_world,ierr)
!OCL NOFLTLD
          do i=1,neg
            do j=1,np_fs
              fsr_ball(i,j,ik)=wk_fsr_ball(j,i)
              if(k_symmetry(ik) /= GAMMA) fsi_ball(i,j,ik)=wk_fsi_ball(j,i)
            enddo
          enddo
          deallocate(wk_fsr_l,wk_fsr_ball)
          !if(k_symmetry(ik) /= GAMMA) deallocate(wk_fsi_l,wk_fsi_ball)
          if(k_symmetry(ik) /= GAMMA) deallocate(wk_fsi_l,wk_fsi_ball)
       end if
       end if

       call m_ES_alloc_fsr_l_2d(np_e,nlmta)
       call m_ES_gather_f_3d_to_2d(fsr_l, fsr_l_2D, ik)

       if(k_symmetry(ik) /= GAMMA) then
         call m_ES_alloc_fsi_l_2d(np_e, nlmta)
         call m_ES_gather_f_3d_to_2d(fsi_l, fsi_l_2D, ik)
       endif
    else

    if (sw_betar_dot_wfs_exp==ON .and. prec) then
      if(.not.allocated(fsr_l_2D)) allocate(fsr_l_2D(np_e,nlmta))
      if(k_symmetry(ik) /= GAMMA .and. .not.allocated(fsi_l_2D)) allocate(fsi_l_2D(np_e,nlmta))
      if(present(map)) then
        call m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,fsr_l_2D,fsi_l_2D,map=map)
      else
        call m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,fsr_l_2D,fsi_l_2D)
      endif
      if(nrank_g > 1) then
        call mpi_reduce_scatter(fsr_l_2D,bpr_l(:,:,ik),nel_fs*np_e,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
        if(k_symmetry(ik) /= GAMMA) then
          call mpi_reduce_scatter(fsi_l_2D,bpi_l(:,:,ik),nel_fs*np_e, mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
        endif
      else
        bpr_l(:,:,ik) = fsr_l_2D(:,:)
        if(k_symmetry(ik) /= GAMMA) bpi_l(:,:,ik) = fsi_l_2D(:,:)
      endif
      return
    endif

! ========================================================================================
                                                                __TIMER_SUB_START(401)
#ifdef __FAPP__
    call fapp_start('betar_dot_WFs_gspace',1,1)
#endif
    call tstatc0_begin('betar_dot_WFs ',id_sname,level=1)
    iksnl = (ik-1)/nspin + 1

!f    fsr_l(:,:,ik) = 0.d0
!f    if(k_symmetry(ik) /= GAMMA) fsi_l(:,:,ik) = 0.d0

!---
    if( nblocksize_betar_npe_is_given) then
       npesize = nblocksize_betar_dot_wfs_npe
       if( npesize > np_e) npesize = np_e
    else
       npesize = np_e
    end if

!FF    call m_ES_alloc_fsr_l_2d()
    if(allocated(fsr_l_2D)) deallocate(fsr_l_2D)
    call m_ES_alloc_fsr_l_2d(npesize,nlmta)
    allocate(wk_fsr_2D(npesize,np_fs))
    fsr_l_2D = 0.0d0
!----------

    allocate(wk_fsr_l(np_fs,np_e))
! === for np_fs == 0 ===========================================================
!   allocate(wk_fsr_ball(np_fs,neg))
    allocate(wk_fsr_ball(max(np_fs,1),neg))
! ==============================================================================
    if(k_symmetry(ik) /= GAMMA) then
!FF        call m_ES_alloc_fsi_l_2d()
      if(allocated(fsi_l_2D)) deallocate(fsi_l_2D)
      call m_ES_alloc_fsi_l_2d(npesize,nlmta)
      fsi_l_2D = 0.0d0
      allocate(wk_fsi_l(np_fs,np_e))
      allocate(wk_fsi_ball(np_fs,neg))
      allocate(wk_fsi_2D(npesize,np_fs))
    endif

    call set_msize() ! -> ibsize, msize_target, msizemax, natm_redmax

   do n = 1,np_e,npesize
    call alloc_wkzfsincos_red(ibsize,natm_redmax)
    call alloc_wkarai(ibsize,msizemax)
!!$    call alloc_wkother(msizemax)
!!$    allocate(mil(msizemax))
!FF    allocate(psi_ri(ibsize,np_e,kimg))
    allocate(psi_ri(ibsize,npesize,kimg))

!----
!   do ibl1 = 1, np_g1k(ik), ibsize
!          ibl2 = min( ibl1+ibsize-1, np_g1k(ik) )
!          call cpzaj_l_to_psi_ri(ibl1,ibl2)
!-----
    ia1 = 1
    msizesum = 0
    Mblock: do
       msize = ilmt(ityp(ia1))
       do ia = ia1+1, natm
          it = ityp(ia)
          msize = msize + ilmt(it)
          if(msize > msize_target) exit
       end do
       ia2 = ia-1

       call pre_lmta_k_blk(msize,LD11,LD12) !  -> msize, LD11, LD12
       call alloc_wkbp(ik,LD11,LD12)
       bp_tmp1=0.d0
       if(k_symmetry(ik) /= GAMMA) bp_tmp2=0.d0

       do ibl1 = 1, np_g1k(ik), ibsize
          ibl2 = min( ibl1+ibsize-1, np_g1k(ik) )
          call cppsi_l_to_psi_ri(ibl1,ibl2,n)
          call G_dot_R_map_blk_3D(ik,ibl1,ibl2, ia1,ia2) ! pos(ia,1:3), ngabc -> wk_zfcos,wk_zfsin
          call m_ES_betar_dot_WFs_4_lmta_k_blk_3D(ik,ibsize,ibl1,ibl2,psi_ri &
               & ,iksnl,snl, msize,LD11,tran1,ia1,ia2) !  --> bp_tmp1, bp_tmp2
       end do

! === DEBUG by tkato 2011/07/12 ================================================
!      call post_lmta_k_blk(ista_k,iend_k, fsr_l_2D,fsi_l_2D)
       if(k_symmetry(ik) /= GAMMA) then
          call post_lmta_k_blk(ista_k,iend_k,n,fsr_l_2D,fsi_l_2D)
       else
          call post_lmta_k_blk(ista_k,iend_k,n,fsr_l_2D)
       endif
! ==============================================================================
       msizesum = msizesum+msize

       call dealloc_wkbp(ik)
       ia1 = ia2+1
       if(ia1 > natm) exit Mblock
    end do Mblock
!----
!   end do
!----

    if(ipribetar >= 2) call wd_fsr_fsi_3D(ik)

    deallocate(psi_ri)
!!$    call dealloc_wkother()
    call dealloc_wkarai()
    call dealloc_wkzfsincos()


!f reduce_scatter f
!modify2010
   if(nrank_g > 1) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,419)
!    call mpi_reduce_scatter(fsr_l_2D,fsr_l(n,1,ik),nel_fs*npesize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
!    if(k_symmetry(ik) /= GAMMA) call mpi_reduce_scatter(fsi_l_2D,fsi_l(n,1,ik),nel_fs*npesize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    call mpi_reduce_scatter(fsr_l_2D,wk_fsr_2D,nel_fs*npesize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    if(k_symmetry(ik) /= GAMMA) call mpi_reduce_scatter(fsi_l_2D,wk_fsi_2D,nel_fs*npesize,&
                                    & mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(419)
     iend = n+npesize-1
     if(iend >= np_e) iend=np_e
!OCL NOFLTLD
     do i = n,iend
     do j = 1,np_fs
        bpr_l(i,j,ik) = wk_fsr_2D(i-n+1,j)
        if(k_symmetry(ik) /= GAMMA)bpi_l(i,j,ik) = wk_fsi_2D(i-n+1,j)
     enddo;enddo
   else
     iend = n+npesize-1
     if(iend >= np_e) iend=np_e
     do i = n,iend
     do j = 1,np_fs
        bpr_l(i,j,ik) = fsr_l_2D(i-n+1,j)
        if(k_symmetry(ik) /= GAMMA)bpi_l(i,j,ik) = fsi_l_2D(i-n+1,j)
     enddo;enddo
!     fsr_l(:,:,ik)=fsr_l_2D(:,:)
!     if(k_symmetry(ik) /= GAMMA) fsi_l(:,:,ik)=fsi_l_2D(:,:)
   endif

  enddo ! npesize

  call tstatc0_end(id_sname)

! for eigen value sort 
  if(present(mod_ball)) then
  if(mod_ball == ON) then
                                                 __TIMER_COMM_START(420)
  do i=1,np_e
    do j=1,np_fs
      wk_fsr_l(j,i)=bpr_l(i,j,ik)
      if(k_symmetry(ik) /= GAMMA) wk_fsi_l(j,i)=bpi_l(i,j,ik)
    enddo
  enddo
                                                 __TIMER_COMM_STOP(420)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,421)
  call mpi_allgatherv(&
       & wk_fsr_l,np_fs*np_e,mpi_double_precision,wk_fsr_ball(1,1),ball_buff,ball_addr,mpi_double_precision,mpi_kg_world,ierr)
  if(k_symmetry(ik) /= GAMMA)call mpi_allgatherv(wk_fsi_l,np_fs*np_e,mpi_double_precision,&
! === DEBUG for -check all =====================================================
!                                 &wk_fsi_ball(1,1),ball_buff,ball_addr,mpi_double_precision,mpi_kg_world,ierr)
                                  &wk_fsi_ball,ball_buff,ball_addr,mpi_double_precision,mpi_kg_world,ierr)
! ==============================================================================
                                                 __TIMER_COMM_STOP(421)
                                                 __TIMER_COMM_START(422)
!OCL NOFLTLD
  do i=1,neg
    do j=1,np_fs
      fsr_ball(i,j,ik)=wk_fsr_ball(j,i)
      if(k_symmetry(ik) /= GAMMA) fsi_ball(i,j,ik)=wk_fsi_ball(j,i)
    enddo
  enddo
                                                 __TIMER_COMM_STOP(422)

  end if
  end if

   deallocate(wk_fsr_l,wk_fsr_ball)
   if(k_symmetry(ik) /= GAMMA) deallocate(wk_fsi_l,wk_fsi_ball)
   deallocfs2d = .true.
   if(present(dealloc_fs2d)) deallocfs2d = dealloc_fs2d

   if(deallocfs2d)then
    deallocate(fsr_l_2D)
    if(k_symmetry(ik) /= GAMMA) deallocate(fsi_l_2D)
   end if
#ifdef __FAPP__
   call fapp_stop('betar_dot_WFs_gspace',1,1)
#endif
                                                                __TIMER_SUB_STOP(401)
! === rspace on 3D_Parallel ==============================================================
    end if
! ========================================================================================

   contains
     subroutine cppsi_l_to_psi_ri(ibl1,ibl2,np)
       integer, intent(in) :: ibl1,ibl2,np
       integer :: ip, i, ri, iend
                                                                __TIMER_SUB_START(402)
                                                                __TIMER_DO_START(423)
       iend = np+npesize-1
       if(iend >= np_e ) iend = np_e
       do ri = 1, kimg
          do i = np, iend
             do ip = ibl1, ibl2
                psi_ri(ip-ibl1+1,i-np+1,ri) = psi_l(ip,i,ik,ri)
             end do
          end do
       end do
                                                                __TIMER_DO_STOP(423)
                                                                __TIMER_SUB_STOP(402)
     end subroutine cppsi_l_to_psi_ri

     subroutine set_msize()
       integer :: ncache, N, K, M, ia_i, ia_f, msize
       integer, save :: msize_print = OFF
                                                                __TIMER_SUB_START(403)

       if( nblocksize_betar_is_given) then
          ibsize = nblocksize_betar_dot_wfs
       else
          ibsize = nb_betar_default
          if(ipribetar >= 2 .and. msize_print == OFF) then
             write(nfout,'(" !ibsize(=nblocksize_betar_w) (m_ES_betar_dot_Psi_4_each_k_3D) = ",i8)') ibsize
          end if
       end if
!f
       if( ibsize == 0 ) ibsize = np_g1k(ik)
!f
!!$    ibsize=10000
    if( nblocksize_betar_nlmta_is_given) then
       M = nblocksize_betar_dot_wfs_nlmta
       if(M > nlmta) M = nlmta
    else
       ncache = (m_CtrlP_cachesize()*1024)*3/4
!F       if(ncache == 0) then
!F          M = nlmta
!F       else
          K = ibsize
          N = np_e
          if(kimg == 1) then
             if(k_symmetry(ik) == GAMMA) then
                M = (ncache - K*N*8)/((K+N)*8)
             else 
                M = (ncache - K*N*8)/((K+N)*16)
             endif
          else
             if(k_symmetry(ik) == GAMMA) then
                M = (ncache - K*N*16)/((2*K+N)*8)
             else
                M = (ncache - K*N*16)/((K+N)*16)
             end if
          endif
          if(M > nlmta) M = nlmta
!F       endif
!!$$!!$$!!
!F       M = nlmta
!!$$!!$$!!
    endif
       msize_target = M
       if(msize_target<np_e) msize_target = np_e
       if(msize_target < ilmt(1)) msize_target = ilmt(1)
       if(ipribetar>=2) then
          if(msize_print == OFF) then
             write(nfout,'(" mype = ",i3, " np_e, msize_target, nlmta, ibsize, np_g1k(ik)/ibsize = ",5i8 &
                  & ," <<m_ES_betar_dot_Psi_4_each_k_3D>>")') &
                  & mype, np_e, msize_target, nlmta, ibsize, np_g1k(ik)/ibsize
          end if
       end if

!F       if(ncache == 0) then
!F          msizemax = nlmta
!F          natm_redmax = natm
!F       else
          natm_redmax = 0
          msizemax = 0
          ia_i = 1
!!$          if(msize_target < ilmt(1)) then
!!$             msize = nlmta
!!$             msizemax = msize
!!$          else
             Mblock: do
                msize = ilmt(ityp(ia_i))
                do ia = ia_i+1, natm
                   it = ityp(ia)
                   msize = msize + ilmt(it)
                   if(msize > msize_target) exit
                end do
                if(ia <= natm) then
                   msize = msize - ilmt(it)
                end if
                ia_f = ia - 1
                if(ia_f - ia_i + 1 > natm_redmax) natm_redmax = ia_f -ia_i + 1
                if(msize > msizemax) msizemax = msize
                if(ipribetar >= 2) then
                   if(msize_print == OFF) &
                        & write(nfout,'(" ia_i, ia_f, msize = ",3i8)') ia_i, ia_f, msize
                end if
                ia_i = ia_f+1
                if(ia_i > natm) exit Mblock
             end do Mblock
!!$          end if
          if(ipribetar>=2) then
             if(msize_print == OFF) &
                  & write(nfout,'(" mype = ",i3, " msizemax, natm_redmax = ",2i8)') mype, msizemax, natm_redmax
          end if
!F       end if
       msize_print = ON
!!$$!!$$!!
!F       msizemax = nlmta
!F       natm_redmax = natm
!!$$!!$$!!
                                                                __TIMER_SUB_STOP(403)
     end subroutine set_msize

     subroutine pre_lmta_k_blk(msize,LD11,LD12)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
       integer,intent(out) :: msize,LD11, LD12
       integer :: icnt, ia, it, lmt1
                                                                __TIMER_SUB_START(404)

       icnt = 0
!!$       do ia = 1, natm
       do ia = ia1, ia2
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             icnt = icnt+1
!!$             mil(icnt) = mod( ltp(lmt1,it),4)
!!$             lmtt_tmp(icnt) = lmtt(lmt1,it)
!!$             ia_tmp(icnt) = ia-ia1+1
          end do
       end do
       msize = icnt
       if(npesize .gt. msize) then
          LD11 = npesize ; LD12 = msize
          tran1 = .true.
       else
          LD11 = msize ; LD12 = npesize
          tran1 = .false.
       end if
#ifdef SX
       if( mod(LD11,2) .eq. 0) LD11 = LD11+1
#endif
                                                                __TIMER_SUB_STOP(404)
     end subroutine pre_lmta_k_blk

    subroutine post_lmta_k_blk( k1,k2,np,bpr_l,bpi_l)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, Octover 2009 : multipication of i^L is removed
      integer,intent(in)  :: k1,k2,np
!f      real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2)  ::  bpr_l,bpi_l  ! MPI
! === DEBUG by tkato 2011/07/12 ================================================
!     real(kind=DP),intent(out),dimension(np_e,nlmta)  ::  bpr_l,bpi_l  ! MPI
      real(kind=DP),intent(out),dimension(npesize,nlmta)          :: bpr_l
      real(kind=DP),intent(out),dimension(npesize,nlmta),optional :: bpi_l
! ==============================================================================
      integer :: j,ib, iend
                                                                __TIMER_SUB_START(405)
                                                                __TIMER_DO_START(424)
       iend = np+npesize-1
       if(iend >= np_e ) iend = np_e

      if(k_symmetry(ik) == GAMMA) then  !! betar_dot_WFs_core2_blk
         if(tran1) then
!cdir nodep
            do j = 1, msize
               do ib = 1, npesize
!f                  bpr_l(ib,j+msizesum,ik) = bp_tmp1(ib,j)
                  bpr_l(ib,j+msizesum) = bp_tmp1(ib,j)
!!                  bpr_l(ib,j+msizesum) = bpr_l(ib,j+msizesum) + bp_tmp1(ib,j)
               enddo
            enddo
         else
!cdir nodep
            do j = 1, msize
               do ib = 1, npesize
!f                  bpr_l(ib,j+msizesum,ik) = bp_tmp1(j,ib)
                  bpr_l(ib,j+msizesum) = bp_tmp1(j,ib)
!!                  bpr_l(ib,j+msizesum) = bpr_l(ib,j+msizesum) + bp_tmp1(j,ib)
               enddo
            enddo
         endif
      else
         if(tran1) then
!cdir nodep
            do j = 1, msize
               do ib = 1, npesize
!f                  bpr_l(ib,j+msizesum,ik) =  bp_tmp1(ib,j)
!f                  bpi_l(ib,j+msizesum,ik) =  bp_tmp2(ib,j)
                  bpr_l(ib,j+msizesum) =  bp_tmp1(ib,j)
                  bpi_l(ib,j+msizesum) =  bp_tmp2(ib,j)
!!                  bpr_l(ib,j+msizesum) =  bpr_l(ib,j+msizesum) + bp_tmp1(ib,j)
!!                  bpi_l(ib,j+msizesum) =  bpi_l(ib,j+msizesum) + bp_tmp2(ib,j)
               end do
            end do
          else
!cdir nodep
            do j = 1, msize
               do ib = 1, npesize
!f                  bpr_l(ib,j+msizesum,ik) =  bp_tmp1(j,ib)
!f                  bpi_l(ib,j+msizesum,ik) =  bp_tmp2(j,ib)
                  bpr_l(ib,j+msizesum) =  bp_tmp1(j,ib)
                  bpi_l(ib,j+msizesum) =  bp_tmp2(j,ib)
!!                  bpr_l(ib,j+msizesum) =  bpr_l(ib,j+msizesum) + bp_tmp1(j,ib)
!!                  bpi_l(ib,j+msizesum) =  bpi_l(ib,j+msizesum) + bp_tmp2(j,ib)
               end do
            end do
          endif
       end if
                                                                __TIMER_DO_STOP(424)
                                                                __TIMER_SUB_STOP(405)
    end subroutine post_lmta_k_blk

  end subroutine m_ES_betar_dot_Psi_k_3D_snl

  subroutine m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,bpr_l,bpi_l,map)
    integer, intent(in) :: ik,k1,k2
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(in),dimension(maxval(np_g1k),nlmta,ista_snl:iend_snl,2) :: i_exp_snl
    real(kind=DP),intent(inout),dimension(np_e,nlmta) ::  bpr_l,bpi_l !MPI
    logical, intent(in), dimension(np_e), optional :: map
    real(kind=DP), allocatable, dimension(:,:,:) :: zajtmp
    real(kind=DP), allocatable, dimension(:,:) :: resr,resi
    integer :: ig,ib,it,iksnl,lmt1,ia,lmtt1,lmta1,icount,ng,mg
    integer :: ndata
    logical :: umap
    integer :: id_sname = -1
    call tstatc0_begin('m_ES_betar_dot_WFs_exp ',id_sname,1)
    iksnl = (ik-1)/nspin + 1         

    umap = .false.
    if(present(map)) umap = .true.
    if(umap)then
      ndata = 0
      do ib=1,np_e
        if(.not.map(ib)) ndata = ndata+1
      enddo
    endif
    if(ndata==np_e) umap = .false.

    if(umap) then
      allocate(zajtmp(np_g1k(ik),ndata,kimg))
      allocate(resr(nlmta,ndata))
      if(k_symmetry(ik)/=GAMMA) allocate(resi(nlmta,ndata))
      icount = 0
      do ib=1,np_e
        if(.not.map(ib))then 
          icount = icount+1
          zajtmp(1:np_g1k(ik),icount,1:kimg) = psi_l(1:np_g1k(ik),ib,ik,1:kimg)
        endif
      enddo 
    else
      allocate(zajtmp(np_g1k(ik),np_e,kimg))
      zajtmp(1:np_g1k(ik),1:np_e,1:kimg) = psi_l(1:np_g1k(ik),1:np_e,ik,1:kimg)
      allocate(resr(nlmta,np_e))
      if(k_symmetry(ik)/=GAMMA) allocate(resi(nlmta,np_e))
    endif

    ng = np_g1k(ik);mg=maxval(np_g1k)
    if(.not.umap) then
      if(kimg==1) then
        call DGEMM__('T','N',nlmta,np_e,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,1),ng, 0.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) &
        call DGEMM__('T','N',nlmta,np_e,ng, 1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,1),ng, 0.d0,resi,nlmta)
      else
        call DGEMM__('T','N',nlmta,np_e,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,1),ng, 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,np_e,ng,-1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,2),ng, 1.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) then
          call DGEMM__('T','N',nlmta,np_e,ng, 1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,1),ng, 0.d0,resi,nlmta)
          call DGEMM__('T','N',nlmta,np_e,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,2),ng, 1.d0,resi,nlmta)
        endif
      endif
      do lmt1=1,nlmta
        do ib=1,np_e
          bpr_l(ib,lmt1) = resr(lmt1,ib)
          if(k_symmetry(ik) /= GAMMA) bpi_l(ib,lmt1) = resi(lmt1,ib)
        enddo
      enddo
    else
      if(kimg==1) then
        call DGEMM__('T','N',nlmta,ndata,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,1),ng, 0.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) &
        call DGEMM__('T','N',nlmta,ndata,ng, 1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,1),ng, 0.d0,resi,nlmta)
      else
        call DGEMM__('T','N',nlmta,ndata,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,1),ng, 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,ndata,ng,-1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,2),ng, 1.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) then
          call DGEMM__('T','N',nlmta,ndata,ng, 1.d0, i_exp_snl(:,:,iksnl,2),mg,zajtmp(:,:,1),ng, 0.d0,resi,nlmta)
          call DGEMM__('T','N',nlmta,ndata,ng, 1.d0, i_exp_snl(:,:,iksnl,1),mg,zajtmp(:,:,2),ng, 1.d0,resi,nlmta)
        endif
      endif
      do lmt1=1,nlmta
        icount = 0
        do ib=1,np_e
          if(.not.map(ib)) then
            icount = icount+1
            bpr_l(ib,lmt1) = resr(lmt1,icount)
            if(k_symmetry(ik) /= GAMMA) bpi_l(ib,lmt1) = resi(lmt1,icount)
          endif
        enddo
      enddo
    endif
    deallocate(zajtmp)
    deallocate(resr)
    if(k_symmetry(ik) /= GAMMA) deallocate(resi)
    call tstatc0_end(id_sname)
  end subroutine m_ES_betar_dot_WFs_exp

#else
  subroutine m_ES_betar_dot_Psi_4_each_k_3D(nfout,ik)
    integer, intent(in)    :: nfout, ik
    integer :: ia, iksnl, mapmode

    print *,"ERROR : Not Support Multi Dimension Para on DGEM"
    MPI_ABORT(mpi_comm_world, 141, ik)

  end subroutine m_ES_betar_dot_Psi_4_each_k_3D
#endif
!! #ifdef NONLOCAL_DGEMM end


  subroutine m_ES_betar_dot_WFs_4_lmta_k_3D(k1,k2,ik,psi_l,ia,iksnl,snl_or_snld,bpr_l,bpi_l,mapmode)
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l ! MPI
!!$    real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l ! MPI
!!$    real(kind=DP),intent(inout),dimension(np_e,np_fs,k1:k2) :: bpr_l,bpi_l ! MPI
    real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l ! MPI
! === DEBUG by tkato 2011/08/31 ================================================
!   real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt,ista_snl:iend_snl) :: snl_or_snld
! ==============================================================================

    integer                    ::it, lmt1, lmtt1, lmta1, il1
    integer :: id_sname = -1
#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) then
       ibsize = iend_g1k(ik) - ista_g1k(ik) + 1
    else
    if(k_symmetry(ik) == GAMMA) then ! core2
      if(kimg == 1) then ! ar(i):1,psi_l(i,ib,ik,1):np_e or
                         ! ar(i):1,psi_l(i,ib,ik,2):np_e
        ibsize=ncache/(8*(np_e+1))
      else !  ar(i),ai(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,2):np_e*2  or
           !  ai(i),ar(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,2):np_e*2
        ibsize=ncache/(8*(2+np_e*2))
      endif
    else ! ar(i),ai(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,1):np_e*2
      ibsize=ncache/(8*(2+np_e*2))
    endif
    end if
! nec debug
!  write(6,887) 'ik=',ik,' np_e=',np_e,' iba(ik)=',iba(ik),' ibsize=',ibsize,&
!      ' ncache=',ncache
! 887 format(a,i3,a,i4,a,i6,a,i6,a,i8)
! NEC tune ------------------------------->
#endif

    call tstatc0_begin('m_ES_betar_dot_WFs_4_lmta_k_3D ',id_sname)

    it    = ityp(ia)
#ifndef SX
! NEC tune ------------------------------->
    do ibl1=ista_g1k(ik),iend_g1k(ik),ibsize
    ibl2=ibl1+ibsize-1
    if(ibl2.gt.iend_g1k(ik)) ibl2=iend_g1k(ik)
! nec debug
!  write(6,888) 'ibl1,ibl2=',ibl1,ibl2,' iba(ik)=',iba(ik)
! 888 format(a,2i6,a,i6)
! NEC tune <-------------------------------
#endif
    do lmt1 = 1, ilmt(it)
       lmtt1 = lmtt(lmt1,it)
       lmta1 = lmta(lmt1,ia)
       il1   = ltp(lmt1,it)
!!$       do so = 1, mso+1
!!$          lmtt1 = lmtt0*mso - (mso-so)
!!$          lmta1 = lmta1*mso - (mso-so)
!!$          lmtt1 = nlmtt*(so-1) + lmtt0
!!$          lmta1 = nlmta*(so-1) + lmta0
       call G_dot_R_mult_snl(snl_or_snld)   ! lmtt1, exp(iGR), zfcos,zfsin,snl_or_snld ->ar,ai
       if(k_symmetry(ik) == GAMMA) then
          call betar_dot_WFs_core2(il1,psi_l,bpr_l)  ! lmta1, ar,ai,psi_l -> bpr_l
       else
          call betar_dot_WFs_core(psi_l,bpr_l,bpi_l) ! lmta1, sum(c(k+G)exp(iGR)*snl() ,ar,ai,psi_l -> bpr_l,bpi_l
#ifdef SX
! NEC tune ------------------------------->
          call multiple_i_l()                 !   i**l*( )
#else
          if(ibl2.eq.iend_g1k(ik))then
            call multiple_i_l()                 !   i**l*( )
          endif
! NEC tune <-------------------------------
#endif
       end if
    end do
#ifndef SX
! NEC tune ------------------------------->
    end do
! NEC tune <-------------------------------
#endif
    call tstatc0_end(id_sname)
  contains
    subroutine G_dot_R_mult_snl(snl_or_snld)
! === DEBUG by tkato 2011/08/31 ================================================
!     real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
      real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt,ista_snl:iend_snl) :: snl_or_snld
! ==============================================================================
      integer :: i, i1, iadd
#ifdef BETAR_GAMMA_TUNE
      integer :: mil
#endif
      integer :: id_sname = -1
      call tstatc0_begin('G_dot_R_mult_snl ',id_sname)

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
#ifdef BETAR_GAMMA_TUNE
         mil = mod(il1,4)
         if(mapmode == MAPPED) then
            if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
               if(mil == 1 .or. mil == 3) then
                  if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
                     ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
                  endif
                  do i = ibl1, ibl2
                     iadd = i - ista_g1k(ik) + 1
! === DEBUG by tkato 2011/08/31 ================================================
!                    ar(iadd) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
!                    ai(iadd) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
                     ar(iadd) = zfcos(i)*snl_or_snld(iadd,lmtt1,iksnl)
                     ai(iadd) = zfsin(i)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
                  end do
               else if(mil == 2 .or. mil == 0) then
                  if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
                     ar(1) = 0.d0
                  endif
                  do i = ibl1, ibl2
                     iadd = i - ista_g1k(ik) + 1
! === DEBUG by tkato 2011/08/31 ================================================
!                    ai(iadd) = -zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
!                    ar(iadd) =  zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
                     ai(iadd) = -zfcos(i)*snl_or_snld(iadd,lmtt1,iksnl)
                     ar(iadd) =  zfsin(i)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
                  end do
               end if
            else
               if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
                  ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
               endif
               do i = ibl1, ibl2
                  iadd = i - ista_g1k(ik) + 1
! === DEBUG by tkato 2011/08/31 ================================================
!                 ar(iadd) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
!                 ai(iadd) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
                  ar(iadd) = zfcos(i)*snl_or_snld(iadd,lmtt1,iksnl)
                  ai(iadd) = zfsin(i)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
               end do
            end if
         else
            if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
               if(mil == 1 .or. mil ==3) then
                  if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
                     i1    = nbase(1,ik)
                     ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
                  endif
                  do i = ibl1, ibl2
                     iadd = i - ista_g1k(ik) + 1
                     i1    = nbase(i,ik)
! === DEBUG by tkato 2011/08/31 ================================================
!                    ar(iadd) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
!                    ai(iadd) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
                     ar(iadd) = zfcos(i1)*snl_or_snld(iadd,lmtt1,iksnl)
                     ai(iadd) = zfsin(i1)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
                  end do
               else if(mil == 2 .or. mil == 0) then
                  if(ista_g1k(ik).eq.1)then
                     ar(1) = 0.d0
                  endif
                  do i = ibl1, ibl2
                     iadd = i - ista_g1k(ik) + 1
                     i1    = nbase(i,ik)
! === DEBUG by tkato 2011/08/31 ================================================
!                    ai(iadd) = -zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
!                    ar(iadd) =  zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
                     ai(iadd) = -zfcos(i1)*snl_or_snld(iadd,lmtt1,iksnl)
                     ar(iadd) =  zfsin(i1)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
                  end do
               else
               end if
            else
               if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
                  i1    = nbase(1,ik)
                  ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
               endif
               do i = ibl1, ibl2
                  iadd = i - ista_g1k(ik) + 1
                  i1    = nbase(i,ik)
! === DEBUG by tkato 2011/08/31 ================================================
!                 ar(iadd) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
!                 ai(iadd) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
                  ar(iadd) = zfcos(i1)*snl_or_snld(iadd,lmtt1,iksnl)
                  ai(iadd) = zfsin(i1)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
               end do
            end if
         end if
#else
         if(mapmode == MAPPED) then
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1) then
               ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
! === DEBUG by tkato 2011/08/31 ================================================
!              ar(iadd) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
!              ai(iadd) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
               ar(iadd) = zfcos(i)*snl_or_snld(iadd,lmtt1,iksnl)
               ai(iadd) = zfsin(i)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
            end do
         else
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
               i1    = nbase(1,ik)
               ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
               i1    = nbase(i,ik)
! === DEBUG by tkato 2011/08/31 ================================================
!              ar(iadd) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
!              ai(iadd) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
               ar(iadd) = zfcos(i1)*snl_or_snld(iadd,lmtt1,iksnl)
               ai(iadd) = zfsin(i1)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
            end do
         end if
#endif
         if(ipribetar>=3 ) then
            write(nfout,'(" --- ar --- <<G_dot_R_mult_snl>>")')
            write(nfout,'(5f20.10)') (ar(i),i=1,iend_g1k(ik)-ista_g1k(ik)+1)
            write(nfout,'(" --- ai --- <<G_dot_R_mult_snl>>")')
            write(nfout,'(5f20.10)') (ai(i),i=1,iend_g1k(ik)-ista_g1k(ik)+1)
         end if
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
               ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
! === DEBUG by tkato 2011/08/31 ================================================
!              ar(iadd) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
               ar(iadd) = zfcos(i)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
            end do
         else
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iend_g1k(ik) .and. ista_g1k(ik).eq.1)then
               i1    = nbase(1,ik)
               ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
               i1    = nbase(i,ik)
! === DEBUG by tkato 2011/08/31 ================================================
!              ar(iadd) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
               ar(iadd) = zfcos(i1)*snl_or_snld(iadd,lmtt1,iksnl)
! ==============================================================================
            end do
         end if
      end if
      call tstatc0_end(id_sname)
    end subroutine G_dot_R_mult_snl

    subroutine betar_dot_WFs_core(psi_l,bpr_l,bpi_l)
      real(kind=DP),intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l ! MPI
!      real(kind=DP),intent(out),dimension(np_e,np_fs,k1:k2) :: bpr_l, bpi_l ! MPI
      real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) :: bpr_l, bpi_l ! MPI
      real(kind=DP) :: bpr, bpi
      integer       :: ib, i, iadd

      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core ',id_sname)

#ifndef SX
! NEC tune --------------------------->
!      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
!      bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      if(ibl1.eq.ista_g1k(ik))then
         bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
         bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      endif
! NEC tune <---------------------------
#endif

! nec debug
!      write(6,999) 'np_e=',np_e,' ik=',ik,' iba(ik)=',iba(ik),' k1,k2=',k1,&
!         ',',k2,' kimg=',kimg,' ibl1,ibl2=',ibl1,',',ibl2,' lmta1=',lmta1
! 999 format(a,i3,a,i1,a,i5,a,i1,a,i1,a,i1,a,i5,a,i5,a,i5)

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
         do ib = 1, np_e ! MPI
#ifdef SX
! NEC tune ------------------------------------------>
            bpr = 0.d0; bpi = 0.d0
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = ibl1, ibl2
! NEC tune <------------------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
               bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,1)
               bpi = bpi + ai(iadd)*psi_l(iadd,ib,ik,1)
!!$               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)
!!$               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
#ifndef SX
         if(ibl2 .eq. iend_g1k(ik)) then
#endif
            call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
            call mpi_allreduce(MPI_IN_PLACE,bpi_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
         endif
#endif
      else if(kimg == 2) then
!!$#ifdef _VECTOR_TUNING_
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
         do ib = 1, np_e ! MPI
#ifdef SX
! NEC tune ------------------------------------------>
            bpr = 0.d0; bpi = 0.d0
            do i = ista_g1k(ik), iend_g1k(ik)
#else
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = ibl1, ibl2
! NEC tune <------------------------------------------
#endif
               iadd = i - ista_g1k(ik) + 1
               bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,1)-ai(iadd)*psi_l(iadd,ib,ik,2)
               bpi = bpi + ai(iadd)*psi_l(iadd,ib,ik,1)+ar(iadd)*psi_l(iadd,ib,ik,2)
!!$               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
!!$               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
#ifndef SX
         if(ibl2 .eq. iend_g1k(ik)) then
#endif
            call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
            call mpi_allreduce(MPI_IN_PLACE,bpi_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
         endif
#endif
!$#else
!!$         do ib = 1, np_e                   ! MPI
!!$            bpr = 0.d0; bpi = 0.d0
!!$            do i = 1, iba(ik)
!!$               bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
!!$               bpi = bpi + ai(i)*psi_l(i,ib,ik,1)
!!$            end do
!!$            do i = 1, iba(ik)
!!$               bpr = bpr - ai(i)*psi_l(i,ib,ik,2)
!!$               bpi = bpi + ar(i)*psi_l(i,ib,ik,2)
!!$            end do
!!$            bpr_l(ib,lmta1,ik) = bpr
!!$            bpi_l(ib,lmta1,ik) = bpi
!!$         end do
!!$#endif
      end if
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core

    subroutine multiple_i_l
      integer mil, ib
      real(kind=DP) :: tempo

      mil = mod(il1,4)
      if(mil == 2) then
         do ib = 1, np_e ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = -tempo
         end do
      else if(mil == 3) then
         do ib = 1, np_e ! MPI
            bpr_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpi_l(ib,lmta1,ik)
         end do
      else if(mil == 0) then
         do ib = 1, np_e ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = tempo
         end do
      end if
    end subroutine multiple_i_l

    subroutine betar_dot_WFs_core2(il1,psi_l,bpr_l)
      integer, intent(in) :: il1
      real(kind=DP),intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l ! MPI
!!$      real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) :: bpr_l ! MPI
!      real(kind=DP),intent(out),dimension(np_e,np_fs,k1:k2) :: bpr_l ! MPI
      real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) :: bpr_l ! MPI
      integer       :: ib, i, mil, iadd
      real(kind=DP) :: bpr, fp

#ifndef SX
! NEC tune ------------------------------->
      integer :: ibl1_2
#endif
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core2 ',id_sname)

#ifdef SX
      bpr_l(1:np_e,lmta1,ik) = 0.d0 ! MPI
#else
      if(ibl1.eq.1) then
        ibl1_2=2
      else
        ibl1_2=ibl1
      endif
      if(ibl1.eq.ista_g1k(ik))then
        bpr_l(1:np_e,lmta1,ik) = 0.d0 ! MPI
      endif
! NEC tune <-------------------------------
#endif
      mil = mod(il1,4)
      if(mil == 0 .or. mil == 1) then
         fp = 1.d0
      else if(mil == 2 .or. mil == 3) then
         fp = -1.d0
      end if
      if(kimg == 1) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = max(2,ista_g1k(ik)), iend_g1k(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  iadd = i - ista_g1k(ik) + 1
                  bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,1)
               end do
#ifdef SX
! NEC tune ------------------------------->
               if(ista_g1k(ik) == 1) then
                  bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
               else
                  bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr)
               endif
#else
               if(ibl2.eq.iend_g1k(ik))then
                  if(ista_g1k(ik) == 1) then
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
                  else
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr)
                  endif
               else
                  bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
#ifndef SX
            if(ibl2 .eq. iend_g1k(ik)) then
#endif
               call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
            endif
#endif
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = max(2,ista_g1k(ik)), iend_g1k(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  iadd = i - ista_g1k(ik) + 1
                  bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
#else
               if(ibl2.eq.iend_g1k(ik))then
                  bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
               else
                  bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
#ifndef SX
            if(ibl2 .eq. iend_g1k(ik)) then
#endif
               call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
            endif
#endif
         end if
      else if(kimg == 2) then
#ifdef BETAR_GAMMA_TUNE

            do ib = 1, np_e
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
                  iadd = i - ista_g1k(ik) + 1
                  bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,1) - ai(iadd)*psi_l(iadd,ib,ik,2)
               end do
               if(ibl2.eq.iend_g1k(ik))then
                  if(ista_g1k(ik) == 1) then
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
                  else
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr)
                  endif
               else
                  bpr_l(ib,lmta1,ik) = bpr
               endif

            end do
            if(ibl2 .eq. iend_g1k(ik)) then
               call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
            endif
#else
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = max(2,ista_g1k(ik)), iend_g1k(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  iadd = i - ista_g1k(ik) + 1
                  bpr = bpr + ar(iadd)*psi_l(iadd,ib,ik,1) - ai(iadd)*psi_l(iadd,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               if(ista_g1k(ik) == 1) then
                  bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
               else
                  bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr)
               endif
#else
               if(ibl2.eq.iend_g1k(ik))then
                  if(ista_g1k(ik) == 1) then
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
                  else
                     bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr)
                  endif
               else
                  bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
#ifndef SX
            if(ibl2 .eq. iend_g1k(ik)) then
#endif
               call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
            endif
#endif
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = max(2,ista_g1k(ik)), iend_g1k(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  iadd = i - ista_g1k(ik) + 1
                  bpr = bpr + ai(iadd)*psi_l(iadd,ib,ik,1) + ar(iadd)*psi_l(iadd,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
#else
               if(ibl2.eq.iend_g1k(ik))then
                  bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
               else
                  bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
#ifndef SX
            if(ibl2 .eq. iend_g1k(ik)) then
#endif
               call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
#ifndef SX
            endif
#endif
         end if
#endif
      end if
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core2
  end subroutine m_ES_betar_dot_WFs_4_lmta_k_3D

  subroutine m_ES_betar_dot_WFs_4_lmta_k_div(k1,k2,ik,ibsize,ibl1,ibl2,psi_l,ia,iksnl,snl_or_snld,bpr_l,bpi_l,mapmode)
    integer, intent(in)    :: k1,k2, ik, ibsize,ibl1,ibl2, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) ::  psi_l        ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l  ! MPI
    real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld

    integer                    ::it, lmt1, lmtt1, lmta1, il1
    integer :: id_sname = -1
    call tstatc0_begin('m_ES_betar_dot_WFs_4_lmta_k_div ',id_sname)

    it    = ityp(ia)

    do lmt1 = 1, ilmt(it)
       lmtt1 = lmtt(lmt1,it)
       lmta1 = lmta(lmt1,ia)
       il1   = ltp(lmt1,it)
       call G_dot_R_mult_snl(snl_or_snld)   ! lmtt1, exp(iGR), zfcos,zfsin,snl_or_snld ->ar,ai
       if(k_symmetry(ik) == GAMMA) then
          call betar_dot_WFs_core2(il1,psi_l,bpr_l)  ! lmta1, ar,ai,psi_l -> bpr_l
       else
          call betar_dot_WFs_core(psi_l,bpr_l,bpi_l) ! lmta1, sum(c(k+G)exp(iGR)*snl() ,ar,ai,psi_l -> bpr_l,bpi_l
#ifndef BETAR_GAMMA_TUNE
          if(ibl2.eq.iba(ik))then
            call multiple_i_l()                 !   i**l*( )
          endif
#endif
       end if
    end do
    call tstatc0_end(id_sname)

  contains
    subroutine G_dot_R_mult_snl(snl_or_snld)
      real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
      integer :: i
#ifdef BETAR_GAMMA_TUNE
      integer :: mil
      real(kind=DP) :: fp
#endif
      integer :: id_sname = -1
      call tstatc0_begin('G_dot_R_mult_snl ',id_sname)

      if(mapmode /= MAPPED) then
         call phase_error_with_msg(nfout,' mapmode /= MAPPED <<G_dot_R_mult_snl in m_ES_betar_dot_WFs_4_lmta_k_div>>'&
                                  ,__LINE__,__FILE__)
      end if

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
#ifdef BETAR_GAMMA_TUNE
         mil = mod(il1,4)
         if(k_symmetry(ik) == GAMMA) then
            if(mil == 0 .or. mil == 1) then
               fp = 2.d0
            else if(mil == 2 .or. mil == 3) then
               fp = -2.d0
            end if
            if(mil == 0 .or. mil == 2) then
               do i = 1, ibl2-ibl1+1
                  ar(i) = fp*zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
                  ai(i) = fp*zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
               end do
               if(ibl1 == 1) then
                  ar(1) = 0.d0; ai(1) = 0.d0
               end if
            else if(mil == 1 .or. mil == 3) then
               do i = 1, ibl2-ibl1+1
                  ar(i) = fp*zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
                  ai(i) = -fp*zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
               end do
               if(ibl1 == 1) then
                  ar(1) = ar(1)*0.5d0; ai(1) = 0.d0
               end if
            end if
         else
            if(mil == 0 .or. mil == 1) then
               fp = 1.d0
            else if(mil == 2 .or. mil == 3) then
               fp = -1.d0
            end if
            if(mil == 0 .or. mil == 2) then
               do i = 1, ibl2-ibl1+1
                  ar(i) =  fp*zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
                  ai(i) = -fp*zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
               end do
            else if(mil == 1 .or. mil == 3) then
               do i = 1, ibl2-ibl1+1
                  ar(i) = fp*zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
                  ai(i) = fp*zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
               end do
            end if
!!$            do i = 1, ibl2-ibl1+1
!!$               ar(i) = zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
!!$               ai(i) = zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
!!$            end do
         end if
#else
         do i = 1, ibl2-ibl1+1
            ar(i) = zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
            ai(i) = zfsin(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
         end do
#endif
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         do i = 1, ibl2-ibl1+1
            ar(i) = zfcos(i)*snl_or_snld(i+ibl1-1,lmtt1,iksnl)
         end do
      end if
      call tstatc0_end(id_sname)
    end subroutine G_dot_R_mult_snl

    subroutine betar_dot_WFs_core(psi_l,bpr_l,bpi_l)
      real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) :: psi_l        ! MPI
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) :: bpr_l,bpi_l  ! MPI
      real(kind=DP) :: bpr, bpi
      integer       :: ib, i
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core ',id_sname)

      if(ibl1.eq.1)then
         bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
         bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      endif

      if(kimg == 1) then
         do ib = 1, np_e                   ! MPI
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = 1, ibl2-ibl1+1
               bpr = bpr + ar(i)*psi_l(i,ib,1)
               bpi = bpi + ai(i)*psi_l(i,ib,1)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
      else if(kimg == 2) then
         do ib = 1, np_e                   ! MPI
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = 1, ibl2-ibl1+1
               bpr = bpr + ar(i)*psi_l(i,ib,1)-ai(i)*psi_l(i,ib,2)
               bpi = bpi + ai(i)*psi_l(i,ib,1)+ar(i)*psi_l(i,ib,2)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
      end if
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core

    subroutine multiple_i_l
      integer mil, ib
      real(kind=DP) :: tempo
      integer :: id_sname = -1
      call tstatc0_begin('multiple_i_l ',id_sname)

      mil = mod(il1,4)
      if(mil == 2) then
         do ib = 1, np_e                 ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = -tempo
         end do
      else if(mil == 3) then
         do ib = 1, np_e                 ! MPI
            bpr_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpi_l(ib,lmta1,ik)
         end do
      else if(mil == 0) then
         do ib = 1, np_e                 ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = tempo
         end do
      end if
      call tstatc0_end(id_sname)
    end subroutine multiple_i_l

    subroutine betar_dot_WFs_core2(il1,psi_l,bpr_l)
      integer, intent(in) :: il1
      real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) :: psi_l        ! MPI
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) :: bpr_l  ! MPI
      integer       :: ib, i, mil
      real(kind=DP) :: bpr, fp

      integer :: ibl1_2
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core2 ',id_sname)

      if(ibl1.eq.1) then
         bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      endif

#ifdef BETAR_GAMMA_TUNE
      if(kimg == 1) then
         do ib = 1, np_e
            bpr = 0.d0
            do i = 1, ibl2-ibl1+1
               bpr = bpr + ar(i)*psi_l(i,ib,1)
            end do
            bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + bpr
         end do
      else if(kimg == 2) then
         do ib = 1, np_e
            bpr = 0.d0
            do i = 1, ibl2-ibl1+1
               bpr = bpr + ar(i)*psi_l(i,ib,1) + ai(i)*psi_l(i,ib,2)
            end do
            bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + bpr
         end do
      end if
#else
      if(ibl1.eq.1) then
        ibl1_2=2
      else
        ibl1_2=1
      endif

      mil = mod(il1,4)
      if(mil == 0 .or. mil == 1) then
         fp = 1.d0
      else if(mil == 2 .or. mil == 3) then
         fp = -1.d0
      end if

      if(kimg == 1) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
            do ib = 1, np_e
               bpr = 0.d0
               do i = ibl1_2, ibl2-ibl1+1
                  bpr = bpr + ar(i)*psi_l(i,ib,1)
               end do
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*2.d0*bpr
               if(ibl1 == 1) bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*psi_l(1,ib,1)*ar(1)
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
            do ib = 1, np_e
               bpr = 0.d0
               do i = ibl1_2, ibl2-ibl1+1
                  bpr = bpr + ar(i)*psi_l(i,ib,2)
               end do
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*2.d0*bpr
            end do
         end if
      else if(kimg == 2) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
            do ib = 1, np_e
               bpr = 0.d0
               do i = ibl1_2, ibl2-ibl1+1
                  bpr = bpr + ar(i)*psi_l(i,ib,1) - ai(i)*psi_l(i,ib,2)
               end do
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*2.d0*bpr
               if(ibl1 == 1) bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*psi_l(1,ib,1)*ar(1)
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
            do ib = 1, np_e
               bpr = bpr_l(ib,lmta1,ik)
               bpr = 0.d0
               do i = ibl1_2, ibl2-ibl1+1
                  bpr = bpr + ai(i)*psi_l(i,ib,1) + ar(i)*psi_l(i,ib,2)
               end do
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + fp*2.d0*bpr
            end do
         end if
      end if
#endif
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core2
  end subroutine m_ES_betar_dot_WFs_4_lmta_k_div

#ifdef NONLOCAL_DGEMM

  subroutine m_ES_betar_dot_WFs_4_lmta_k_blk_3D(ik,ibsize,ibl1,ibl2,psi_l,iksnl,snl_or_snld &
       &                               ,msize,LD11, tran1, ia1, ia2 )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
    integer, intent(in) :: ik, ibsize,ibl1,ibl2, iksnl
    integer, intent(in) :: msize,LD11,ia1,ia2
    logical, intent(in) :: tran1
    real(kind=DP),intent(in), dimension(ibsize,npesize,kimg) ::  psi_l
!f    real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) ::  snl_or_snld
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt,ista_snl:iend_snl) ::  snl_or_snld
!    real(kind=DP),intent(in), dimension(np_g1k(ik),nlmtt,iksnl:iksnl) ::  snl_or_snld
!    real(kind=DP),intent(in), dimension(kg1,nlmtt,iksnl:iksnl) ::  snl_or_snld
    integer :: id_sname = -1
                                                                __TIMER_SUB_START(406)
!!$    integer :: id_sname2 = -1
    call tstatc0_begin('m_ES_betar_dot_WFs_4_lmta_k_blk ',id_sname)

!!$    if(npes > 1) then
!!$       call tstatc0_begin('mpi_barrier(betar_dot_WFs) ',id_sname2)
!!$       call mpi_barrier(mpi_k_world(myrank_k),ierr)
!!$       call tstatc0_end(id_sname2)
!!$    end if

       call G_dot_R_mult_snl_blk(snl_or_snld)      ! lmtt1, exp(iGR), zfcos,zfsin,snl_or_snld ->ar,ai
       if(k_symmetry(ik) == GAMMA) then
          call betar_dot_WFs_core2_blk(ibsize,tran1,psi_l )  ! lmta1, ar,ai,psi_l -> bp_tmp1
       else
          call betar_dot_WFs_core_blk (ibsize,tran1, psi_l )  ! lmta1, sum(c(k+G)exp(iGR)*snl() ,ar,ai,psi_l -> bp_tmp1, bp_tmp2
       end if

    call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(406)
  contains
    subroutine G_dot_R_mult_snl_blk(snl_or_snld)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, Octover 2009 : mapping algorithm with the subroutine "post_lmta_k_blk"
!
!f      real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
      real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt,ista_snl:iend_snl) :: snl_or_snld
!    real(kind=DP),intent(in), dimension(np_g1k(ik),nlmtt,iksnl:iksnl) ::  snl_or_snld
!    real(kind=DP),intent(in), dimension(kg1,nlmtt,iksnl:iksnl) ::  snl_or_snld
      integer :: j, i
      integer :: mil0, ia, it, lmt1, lt
      real(kind=DP) :: fp, fp0
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(407)
      call tstatc0_begin('G_dot_R_mult_snl_blk ',id_sname)

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
         if(k_symmetry(ik) == GAMMA) then
            fp0 = 2.d0
         else
            fp0 = 1.d0
         end if

         j = 0
                                                                __TIMER_DO_START(425)
         do ia = 1, ia2-ia1+1
            it = ityp(ia+ia1-1)
            do lmt1 = 1, ilmt(it)
               j = j + 1
               mil0 = mod(ltp(lmt1,it),4)
               lt = lmtt(lmt1,it)
               if(mil0 == 0 .or. mil0 == 1) then
                  fp = fp0
               else if(mil0 == 2 .or. mil0 == 3) then
                  fp = -fp0
               end if
               if(k_symmetry(ik) == GAMMA) then
                  if(mil0 == 0 .or. mil0 == 2) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) = fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                     if(ibl1 == 1 .and. myrank_g == 0) then
                        wk_ar(1, j) = 0.d0; wk_ai(1, j) = 0.d0
                     end if
                  else if(mil0 == 1 .or. mil0 == 3) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) = fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = -fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                     if(ibl1 == 1 .and. myrank_g == 0) then
                        wk_ar(1, j) = wk_ar(1,j)*0.5d0; wk_ai(1, j) = 0.d0
                     end if
                  end if
               else
                  if(mil0 == 0 .or. mil0 == 2) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) =  fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = -fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                  else if(mil0 == 1 .or. mil0 == 3) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) = fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                  end if
               end if
            end do
         end do
                                                                __TIMER_DO_STOP(425)
      else ! if(kimg== 1 .and. k_symmetry(ik) == GAMMA) then
         j = 0
                                                                __TIMER_DO_START(426)
         do ia = 1, ia2-ia1+1
            it = ityp(ia+ia1-1)
            do lmt1 = 1, ilmt(it)
               j = j+1
               lt = lmtt(lmt1,it)
               do i = 1, ibl2-ibl1+1
                  wk_ar(i,j) = wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
               end do
            end do
         end do
                                                                __TIMER_DO_STOP(426)
      end if
      call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(427)
    end subroutine G_dot_R_mult_snl_blk

    subroutine betar_dot_WFs_core_blk( ibsize, tran1,  psi_l )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer,intent(in)  :: ibsize
      logical, intent(in)  :: tran1
      real(kind=DP),intent(in), dimension(ibsize,npesize,kimg) :: psi_l
      integer       :: icsize, M, LDA, LDB
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(408)
      call tstatc0_begin('betar_dot_WFs_core_blk ',id_sname)

      M = msize
      icsize = ibl2-ibl1+1
      LDA = ibsize;   LDB = ibsize
      if(kimg == 1) then
         if(tran1) then
                                                              __TIMER_DGEMM_START(411)
            call DGEMM__('T','N',npesize,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',npesize,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
                                                              __TIMER_DGEMM_STOP(411)
         else
                                                              __TIMER_DGEMM_START(412)
            call DGEMM__('T','N',M,npesize,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,npesize,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
                                                              __TIMER_DGEMM_STOP(412)
         end if
      else if(kimg == 2) then
         if(tran1) then
                                                              __TIMER_DGEMM_START(413)
            call DGEMM__('T','N',npesize,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',npesize,M,icsize,-1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',npesize,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
            call DGEMM__('T','N',npesize,M,icsize, 1.d0,psi_l(1,1,2),LDB,wk_ar,LDA,1.d0,bp_tmp2,LD11)
                                                              __TIMER_DGEMM_STOP(413)
         else
                                                              __TIMER_DGEMM_START(414)
            call DGEMM__('T','N',M,npesize,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,npesize,icsize,-1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,npesize,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
            call DGEMM__('T','N',M,npesize,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp2,LD11)
                                                              __TIMER_DGEMM_STOP(414)
         endif
      end if
      call tstatc0_end(id_sname)
                                                                __TIMER_SUB_STOP(408)
    end subroutine betar_dot_WFs_core_blk

    subroutine betar_dot_WFs_core2_blk( ibsize, tran1,  psi_l )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in)  :: ibsize
      logical, intent(in)  :: tran1
      real(kind=DP),intent(in), dimension(ibsize,npesize,kimg) :: psi_l  ! MPI

      integer       :: icsize, M, LDA, LDB
      integer :: id_sname = -1
                                                                __TIMER_SUB_START(409)
      call tstatc0_begin('betar_dot_WFs_core2_blk ',id_sname)

      icsize = ibl2-ibl1+1
      M = msize
      LDA = ibsize; LDB = ibsize
      if(kimg == 1) then
         if(tran1) then
                                                              __TIMER_DGEMM_START(415)
            call DGEMM__('T','N',npesize,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
                                                              __TIMER_DGEMM_STOP(415)
         else
                                                              __TIMER_DGEMM_START(416)
            call DGEMM__('T','N',M,npesize,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
                                                              __TIMER_DGEMM_STOP(416)
         endif
      else if(kimg==2) then
         if(tran1) then
                                                              __TIMER_DGEMM_START(417)
            call DGEMM__('T','N',npesize,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',npesize,M,icsize,1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
                                                              __TIMER_DGEMM_STOP(417)
         else
                                                              __TIMER_DGEMM_START(418)
            call DGEMM__('T','N',M,npesize,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,npesize,icsize,1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
                                                              __TIMER_DGEMM_STOP(418)
         end if
      end if
      call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(409)
    end subroutine betar_dot_WFs_core2_blk
  end subroutine m_ES_betar_dot_WFs_4_lmta_k_blk_3D


#endif
!! #ifdef NONLOCAL_DGEMM end

  subroutine m_ES_add_betar_dot_WFs(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, iksnl, mapmode
    integer     :: id_sname = -1

    if(sw_use_add_proj == OFF) return

    call tstatc0_begin('m_ES_add_betar_dot_WFs ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_add_betar_dot_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai_3D(0)
    call alloc_zfsincos_mpi()
    do ia = 1, natm
! ================================== modified by K. Tagami ============= 11.0
!       if(kv3/nspin == 1) then
!          call G_dot_R_map(ia,1)
!          mapmode = MAPPED
!       else
!          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
!          mapmode = NOTMAPPED
!       end if
          if ( kv3/nspin == 1) then
            mapmode = MAPPED
            call G_dot_R(ia,mapmode,1)
         else
            mapmode = NOTMAPPED
            call G_dot_R(ia,mapmode)
         end if
! ====================================================================== 11.0

       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI

! ================================= modified by K. Tagami ============== 11.0
!          iksnl = (ik-1)/nspin + 1
!
          if ( noncol ) then
            iksnl = ( ik-1 ) / ndim_spinor + 1
          else
            iksnl = (ik-1)/nspin + 1
          endif
! ====================================================================== 11.0

          call m_ES_add_betar_dot_WFs_4_lmta_k(ista_k,iend_k,ik,zaj_l,ia,iksnl,snl_add,fsr_add_l,fsi_add_l,mapmode)
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai_3D()
    call m_ES_dealloc_zfsincos()
    if(ipribetar >= 2) then
       write(nfout,'(" --- fsr_add_l, fsi_add_l ---")')
       do ik = ista_k, iend_k                              ! MPI
          call wd_fsr_fsi_add(ik)                          ! MPI
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine m_ES_add_betar_dot_WFs

  subroutine wd_fsr_fsi_add(ik)
    integer, intent(in) :: ik
    integer :: it, ia
    if(k_symmetry(ik) == GAMMA) then
       write(nfout,'(" -- fsr_add --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (fsr_add_l(map_z(it),ia,ik),ia=1,nlmta_add)
       end do
    else
       write(nfout,'(" -- fsr_add, fsi_add --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (fsr_add_l(map_z(it),ia,ik),fsi_add_l(map_z(it),ia,ik),ia=1,nlmta_add)!MPI
       end do
    end if
  end subroutine wd_fsr_fsi_add

  subroutine m_ES_add_betar_dot_WFs_4_lmta_k(k1,k2,ik,psi_l,ia,iksnl,snl_add,bpr_l,bpi_l,mapmode)
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l        ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta_add,k1:k2) ::    bpr_l,bpi_l  ! MPI
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt_add,ista_snl:iend_snl) :: snl_add

    integer                    ::it, lmt1, lmtt1, lmta1, il1

    it    = ityp(ia)
    do lmt1 = 1, ilmt_add(it)
       lmtt1 = lmtt_add(lmt1,it)
       lmta1 = lmta_add(lmt1,ia)
       il1   = ltp_add(lmt1,it)
       call G_dot_R_mult_snl()             ! exp(iGR)
       if(k_symmetry(ik) == GAMMA) then
          call betar_dot_WFS_core2(il1)
       else
          call betar_dot_WFs_core() ! sum(c(k+G)exp(iGR)*snl_add()
          call multiple_i_l()                 !   i**l*( )
       end if
    end do
  contains
    subroutine G_dot_R_mult_snl()
      integer :: i, i1, iadd

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA)) then
         if(mapmode == MAPPED) then
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               ar(iadd) = zfcos(i)*snl_add(iadd,lmtt1,iksnl)
               ai(iadd) = zfsin(i)*snl_add(iadd,lmtt1,iksnl)
            end do
         else
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               i1 = nbase(i,ik)
               ar(iadd) = zfcos(i1)*snl_add(iadd,lmtt1,iksnl)
               ai(iadd) = zfsin(i1)*snl_add(iadd,lmtt1,iksnl)
            end do
         end if
      else ! if(kimg == 1 .and. k_symmetry(ik) == GAMMA) then
         if(mapmode == MAPPED) then
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               ar(iadd) = zfcos(i)*snl_add(iadd,lmtt1,iksnl)
            end do
         else
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               i1 = nbase(i,ik)
               ar(iadd) = zfcos(i1)*snl_add(iadd,lmtt1,iksnl)
            end do
         end if
      end if
    end subroutine G_dot_R_mult_snl

    subroutine betar_dot_WFs_core()
      integer       :: ib, i, iadd

      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI

      if(kimg == 1) then
         do ib = 1, np_e                   ! MPI
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(iadd)*psi_l(iadd,ib,ik,1)
               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(iadd)*psi_l(iadd,ib,ik,1)
            end do
         end do
      else if(kimg == 2) then
         do ib = 1, np_e                   ! MPI
            do i = ista_g1k(ik), iend_g1k(ik)
               iadd = i - ista_g1k(ik) + 1
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(iadd)*psi_l(iadd,ib,ik,1)-ai(iadd)*psi_l(iadd,ib,ik,2)
               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(iadd)*psi_l(iadd,ib,ik,1)+ar(iadd)*psi_l(iadd,ib,ik,2)
            end do
         end do
      end if
      call mpi_allreduce(MPI_IN_PLACE,bpr_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
      call mpi_allreduce(MPI_IN_PLACE,bpi_l(1,lmta1,ik),np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
    end subroutine betar_dot_WFs_core

    subroutine multiple_i_l
      integer mil, ib
      real(kind=DP) :: tempo

      mil = mod(il1,4)
      if(mil == 2) then
         do ib = 1, np_e                 ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = -tempo
         end do
      else if(mil == 3) then
         do ib = 1, np_e                 ! MPI
            bpr_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpi_l(ib,lmta1,ik)
         end do
      else if(mil == 0) then
         do ib = 1, np_e                 ! MPI
            tempo = bpi_l(ib,lmta1,ik)
            bpi_l(ib,lmta1,ik) = -bpr_l(ib,lmta1,ik)
            bpr_l(ib,lmta1,ik) = tempo
         end do
      end if
    end subroutine multiple_i_l

    subroutine betar_dot_WFs_core2(il1)
      integer, intent(in) :: il1

      integer       :: ib, i, mil
      real(kind=DP) :: bpr, fp

      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI

      mil = mod(il1,4)
      if(mil == 0 .or. mil == 1) then
         fp = 1.d0
      else if(mil == 2 .or. mil == 3) then
         fp = -1.d0
      end if
      if(kimg == 1) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
               end do
               bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,2)
               end do
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
            end do
         end if
      else if(kimg == 2) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1) - ai(i)*psi_l(i,ib,ik,2)
               end do
               bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ai(i)*psi_l(i,ib,ik,1) + ar(i)*psi_l(i,ib,ik,2)
               end do
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
            end do
         end if
      end if
    end subroutine betar_dot_WFs_core2
  end subroutine m_ES_add_betar_dot_WFs_4_lmta_k


  subroutine add_overlap_phirt_psirpw
    integer :: ik, ia, it, lmt1, lmt2, p, q, s
    integer :: il1, il2, im1, im2, t1, t2
    integer :: ikphig
    real(kind=DP) :: ctmp, c1

    Do ik=1, kv3
       if ( map_k(ik) /= myrank_k ) cycle

       ikphig = ( ik -1 )/nspin +1

       Do ia=1, natm
          it = ityp(ia)

          Do lmt1=1, ilmt_phi(it)
             il1=ltp_phi(lmt1,it);  im1 = mtp_phi(lmt1,it);  t1 = taup_phi(lmt1,it)
             p = lmta_phi(lmt1,ia);  s = lmtt_phi(lmt1,it)
             c1 = sqrt( norm_phig( s,ikphig ) )

             Do lmt2=1, ilmt(it)
                il2=ltp(lmt2,it);  im2 = mtp(lmt2,it);  t2 = taup(lmt2,it)
                if ( il1 /= il2 ) cycle
                if ( im1 /= im2 ) cycle

                q = lmta(lmt2,ia)
                ctmp = q_phirt_pw(il1,t1,t2,it) /c1

                if ( k_symmetry(ik) == GAMMA ) then
                   compr_l(:,p,1,ik) = compr_l(:,p,1,ik) &
                        &              +ctmp *fsr_l(:,q,ik) 
                else
                   compr_l(:,p,1,ik) = compr_l(:,p,1,ik) +ctmp *fsr_l(:,q,ik)
                   compi_l(:,p,1,ik) = compi_l(:,p,1,ik) +ctmp *fsi_l(:,q,ik)
                endif
             End Do
          End Do
       End Do
    End Do
  end subroutine add_overlap_phirt_psirpw

  subroutine set_compri_if_anfiferro
    integer :: ik, iorb, jorb, mm, iopr

    iopr = nopr +af

    Do ik=1, kv3, af+1
       if ( map_k(ik) /= myrank_k ) cycle

       compr_l( :,:,1,ik+1 ) = 0.d0
       compi_l( :,:,1,ik+1 ) = 0.d0

       if (k_symmetry(ik) == GAMMA) then
          do iorb=1,nlmta_phi
             do mm=1,nrorb(iorb,iopr)
                jorb = irorb(mm,iorb,iopr)
                compr_l( :,iorb,1,ik+1 ) = compr_l( :,iorb,1,ik+1 ) &
                     &                   + compr_l( :,jorb,1,ik )*crorb(mm,iorb,iopr)
             end do
          end do
       else
          do iorb=1,nlmta_phi
             do mm=1,nrorb(iorb,iopr)
                jorb=irorb(mm,iorb,iopr)
                compr_l( :,iorb,1,ik+1 ) = compr_l( :,iorb,1,ik+1 ) &
                     &                   + compr_l( :,jorb,1,ik )*crorb(mm,iorb,iopr)
                compi_l( :,iorb,1,ik+1 ) = compi_l( :,iorb,1,ik+1 ) &
                     &                   + compi_l( :,jorb,1,ik )*crorb(mm,iorb,iopr)
             end do
          end do
      endif
    end Do

  end subroutine set_compri_if_anfiferro

  subroutine wd_compr_compi(ik)
    integer, intent(in) :: ik
    integer :: it, ia
    if(k_symmetry(ik) == GAMMA) then
       write(nfout,'(" -- compr --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (compr_l(map_z(it),ia,1,ik),ia=1,nlmta_phi)!MPI
       end do
    else
       write(nfout,'(" -- compr, compi --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (compr_l(map_z(it),ia,1,ik),compi_l(map_z(it),ia,1,ik),ia=1,nlmta_phi)!MPI
       end do
    end if
  end subroutine wd_compr_compi

  subroutine m_ES_phir_dot_WFs_3D(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, ikphig, mapmode
    integer     :: id_sname = -1
    call tstatc0_begin('m_ES_phir_dot_WFs_3D ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_phir_dot_WFs_3D --")')
    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai_3D(0)
    call alloc_zfsincos_mpi()
    do ia = 1, natm
! ============================ modified by K. Tagami ================ 11.0
!       if(kv3/nspin == 1) then
!          call G_dot_R_map(ia,1)
!          mapmode = MAPPED
!       else
!          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
!          mapmode = NOTMAPPED
!       end if
!
       if(kv3/nspin == 1) then
!!$          call G_dot_R_map(ia,1)
          mapmode = MAPPED
          call G_dot_R(ia,mapmode,1)
       else
!!$          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = NOTMAPPED
          call G_dot_R(ia,mapmode)
       end if
! ====================================================================== 11.0

       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI
! ==================================== modified by K. Tagami ============ 11.0
!          ikphig = (ik-1)/nspin + 1
          if ( noncol ) then
            ikphig = (ik-1)/ndim_spinor + 1
          else
            ikphig = (ik-1)/nspin + 1
          endif
! ====================================================================== 11.0          

          call m_ES_phir_dot_WFs_4_lmta_k_3D( ista_k, iend_k, ik, zaj_l, ia, ikphig, &
               &                           phig, nopr, compr_l, compi_l, mapmode )
       end do
    end do

    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai_3D()
    call m_ES_dealloc_zfsincos()

    if ( orb_popu_method == 2 )  call add_overlap_phirt_psirpw
    if ( af == 1 ) call set_compri_if_anfiferro

    if(ipribetar >= 2) then
       write(nfout,'(" --- compr_l, compi_l ---")')
       do ik = ista_k, iend_k                            ! MPI
          call wd_compr_compi(ik)                        ! MPI
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine m_ES_phir_dot_WFs_3D

  subroutine m_ES_phir_dot_WFs_4_lmta_k_3D( k1, k2, ik, psi_l, ia, ikphig, phig, &
       &                                 nopr, cr_l, ci_l, mapmode )
    integer, intent(in)    :: k1,k2, ik, ia, ikphig, mapmode, nopr
    real(kind=DP),intent(in), dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l        ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta_phi,nopr,k1:k2) ::    cr_l,ci_l  ! MPI
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt_phi,ista_snl:iend_snl) :: phig

    integer                    ::it, lmt1, lmtt1, lmta1, il1, ierr

    it    = ityp(ia)
    do lmt1 = 1, ilmt_phi(it)
       lmtt1 = lmtt_phi(lmt1,it)
       lmta1 = lmta_phi(lmt1,ia)
       il1   = ltp_phi(lmt1,it)

       if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
          call G_dot_R_mult_phig_unfolding()             ! exp(iGR)
       else
          call G_dot_R_mult_phig()             ! exp(iGR)
       endif

       if(k_symmetry(ik) == GAMMA) then
          call phir_dot_WFs_core2(il1) ! sum(c(k+G)exp(iGR)*phig()
       else
          call phir_dot_WFs_core() ! sum(c(k+G)exp(iGR)*phig()
          call multiple_i_l()                 !   i**l*( )
       end if
    end do

  contains

    subroutine G_dot_R_mult_phig()
      integer :: i, i1, iadd

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, np_g1k(ik)
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, np_g1k(ik)
               i1    = nbase(i+ista_g1k(ik)-1,ik)
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, np_g1k(ik)
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, np_g1k(ik)
               i1    = nbase(i+ista_g1k(ik)-1,ik)
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      end if
    end subroutine G_dot_R_mult_phig

    subroutine G_dot_R_mult_phig_unfolding()
      integer :: i, i1, iadd

      ar = 0.0d0;  ai = 0.0d0

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, np_g1k(ik)
               i1 = i +ista_g1k(ik) -1
               if ( sw_force_kpt_inside_bz == ON ) then
                  if ( GVec_on_refcell(i1,ik) == 0 ) cycle
               else
                  if ( GVec_on_refcell(i1,1) == 0 ) cycle
               endif
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, np_g1k(ik)
               i1    = nbase(i+ista_g1k(ik)-1,ik)
               if ( sw_force_kpt_inside_bz == ON ) then
                  if ( GVec_on_refcell(i1,ik) == 0 ) cycle
               else
                  if ( GVec_on_refcell(i1,1) == 0 ) cycle
               endif
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, np_g1k(ik)
               i1 = i +ista_g1k(ik) -1
               if ( sw_force_kpt_inside_bz == ON ) then
                  if ( GVec_on_refcell(i1,ik) == 0 ) cycle
               else
                  if ( GVec_on_refcell(i1,1) == 0 ) cycle
               endif
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, np_g1k(ik)
               i1    = nbase(i+ista_g1k(ik)-1,ik)
               if ( sw_force_kpt_inside_bz == ON ) then
                  if ( GVec_on_refcell(i1,ik) == 0 ) cycle
               else
                  if ( GVec_on_refcell(i1,1) == 0 ) cycle
               endif
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      end if
    end subroutine G_dot_R_mult_phig_unfolding

    subroutine phir_dot_WFs_core()
      integer       :: ib, i
      real(kind=DP), allocatable, dimension(:) :: cr,ci
      allocate(cr(np_e));cr=0.0d0
      allocate(ci(np_e));ci=0.d0
      cr_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI
      ci_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
         do ib = 1, np_e                   ! MPI
            do i = 1, np_g1k(ik)
               cr(ib) = cr(ib) + ar(i)*psi_l(i,ib,ik,1)
               ci(ib) = ci(ib) + ai(i)*psi_l(i,ib,ik,1)
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
            do ib = 1, np_e                   ! MPI
               do i = 2, np_g1k(ik)
                  cr(ib) = cr(ib) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
               end do
               cr(ib) = cr(ib)*2.d0 + ar(1)*psi_l(1,ib,ik,1)-ai(1)*psi_l(1,ib,ik,2)
            end do
         else
#ifdef VPP
*vocl loop, unroll(4)
#endif
            do ib = 1, np_e                   ! MPI
               do i = 1, np_g1k(ik)
                  cr(ib) = cr(ib) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
                  ci(ib) = ci(ib) + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
               end do
            end do
         end if
      end if
      call mpi_allreduce(MPI_IN_PLACE,cr,np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
      call mpi_allreduce(MPI_IN_PLACE,ci,np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
      cr_l(1:np_e,lmta1,1,ik) = cr(1:np_e)
      ci_l(1:np_e,lmta1,1,ik) = ci(1:np_e)
      deallocate(cr,ci)
    end subroutine phir_dot_WFs_core

    subroutine multiple_i_l
      integer mil, ib
      real(kind=DP) :: tempo

      mil = mod(il1,4)
      if(mil == 2) then
         do ib = 1, np_e                 ! MPI
            tempo = ci_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik)
            cr_l(ib,lmta1,1,ik) = -tempo
         end do
      else if(mil == 3) then
         do ib = 1, np_e                 ! MPI
            cr_l(ib,lmta1,1,ik) = -cr_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = -ci_l(ib,lmta1,1,ik)
         end do
      else if(mil == 0) then
         do ib = 1, np_e                 ! MPI
            tempo = ci_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = -cr_l(ib,lmta1,1,ik)
            cr_l(ib,lmta1,1,ik) = tempo
         end do
      end if
    end subroutine multiple_i_l

    subroutine phir_dot_WFs_core2(il1)
      integer, intent(in) :: il1
      integer       :: ib, i, mil
      real(kind=DP) :: bpr, fp
      real(kind=DP), allocatable, dimension(:) :: cr
      allocate(cr(np_e));cr=0.d0
      cr_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI

      mil = mod(il1,4)
      if(mil == 0 .or. mil == 1) then
         fp = 1.d0
      else if(mil == 2 .or. mil == 3) then
         fp = -1.d0
      end if
      if(kimg == 1) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e                   ! MPI
               bpr = 0.d0
               do i = 2, np_g1k(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
               end do
               cr(ib) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, np_g1k(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,2)
               end do
               cr(ib) = fp*2.d0*bpr
            end do
         end if
      else if(kimg == 2) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, np_g1k(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1) - ai(i)*psi_l(i,ib,ik,2)
               end do
               cr(ib) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, np_g1k(ik)
                  bpr = bpr + ai(i)*psi_l(i,ib,ik,1) + ar(i)*psi_l(i,ib,ik,2)
               end do
               cr(ib) = fp*2.d0*bpr
            end do
         end if
      end if
      call mpi_allreduce(MPI_IN_PLACE,cr,np_e,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_ke_world,ierr)
      cr_l(:,lmta1,1,ik) = cr(:)
      deallocate(cr)
    end subroutine phir_dot_WFs_core2
  end subroutine m_ES_phir_dot_WFs_4_lmta_k_3D

  subroutine m_ES_phir_dot_WFs_4_lmta_k( k1, k2, ik, psi_l, ia, ikphig, phig, &
       &                                 nopr, cr_l, ci_l, mapmode )
    integer, intent(in)    :: k1,k2, ik, ia, ikphig, mapmode, nopr
    real(kind=DP),intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l        ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta_phi,nopr,k1:k2) ::    cr_l,ci_l  ! MPI
    real(kind=DP),intent(in), dimension(kg1,nlmtt_phi,ista_snl:iend_snl) :: phig

    integer                    ::it, lmt1, lmtt1, lmta1, il1

    it    = ityp(ia)
    do lmt1 = 1, ilmt_phi(it)
       lmtt1 = lmtt_phi(lmt1,it)
       lmta1 = lmta_phi(lmt1,ia)
       il1   = ltp_phi(lmt1,it)
       call G_dot_R_mult_phig()             ! exp(iGR)
       if(k_symmetry(ik) == GAMMA) then
          call phir_dot_WFs_core2(il1) ! sum(c(k+G)exp(iGR)*phig()
       else
          call phir_dot_WFs_core() ! sum(c(k+G)exp(iGR)*phig()
          call multiple_i_l()                 !   i**l*( )
       end if
    end do
  contains
    subroutine G_dot_R_mult_phig()
      integer :: i, i1

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
               ai(i) = zfsin(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
               ar(i) = zfcos(i)*phig(i,lmtt1,ikphig)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
               ar(i) = zfcos(i1)*phig(i,lmtt1,ikphig)
            end do
         end if
      end if
    end subroutine G_dot_R_mult_phig

    subroutine phir_dot_WFs_core()
      integer       :: ib, i

      cr_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI
      ci_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
         do ib = 1, np_e                   ! MPI
            do i = 1, iba(ik)
               cr_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik) + ar(i)*psi_l(i,ib,ik,1)
               ci_l(ib,lmta1,1,ik) = ci_l(ib,lmta1,1,ik) + ai(i)*psi_l(i,ib,ik,1)
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
            do ib = 1, np_e                   ! MPI
               do i = 2, iba(ik)
                  cr_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
               end do
               cr_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik)*2.d0 + ar(1)*psi_l(1,ib,ik,1)-ai(1)*psi_l(1,ib,ik,2)
            end do
         else
#ifdef VPP
*vocl loop, unroll(4)
#endif
            do ib = 1, np_e                   ! MPI
               do i = 1, iba(ik)
                  cr_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
                  ci_l(ib,lmta1,1,ik) = ci_l(ib,lmta1,1,ik) + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
               end do
            end do
         end if
      end if
    end subroutine phir_dot_WFs_core

    subroutine multiple_i_l
      integer mil, ib
      real(kind=DP) :: tempo

      mil = mod(il1,4)
      if(mil == 2) then
         do ib = 1, np_e                 ! MPI
            tempo = ci_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = cr_l(ib,lmta1,1,ik)
            cr_l(ib,lmta1,1,ik) = -tempo
         end do
      else if(mil == 3) then
         do ib = 1, np_e                 ! MPI
            cr_l(ib,lmta1,1,ik) = -cr_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = -ci_l(ib,lmta1,1,ik)
         end do
      else if(mil == 0) then
         do ib = 1, np_e                 ! MPI
            tempo = ci_l(ib,lmta1,1,ik)
            ci_l(ib,lmta1,1,ik) = -cr_l(ib,lmta1,1,ik)
            cr_l(ib,lmta1,1,ik) = tempo
         end do
      end if
    end subroutine multiple_i_l

    subroutine phir_dot_WFs_core2(il1)
      integer, intent(in) :: il1
      integer       :: ib, i, mil
      real(kind=DP) :: bpr, fp

      cr_l(1:np_e,lmta1,1,ik) = 0.d0        ! MPI

      mil = mod(il1,4)
      if(mil == 0 .or. mil == 1) then
         fp = 1.d0
      else if(mil == 2 .or. mil == 3) then
         fp = -1.d0
      end if
      if(kimg == 1) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e                   ! MPI
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
               end do
               cr_l(ib,lmta1,1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,2)
               end do
               cr_l(ib,lmta1,1,ik) = fp*2.d0*bpr
            end do
         end if
      else if(kimg == 2) then
         if(mil == 1 .or. mil == 3) then ! l = 0, 4, 8, ... or 2, 6, 10, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1) - ai(i)*psi_l(i,ib,ik,2)
               end do
               cr_l(ib,lmta1,1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
            end do
         else if(mil == 2 .or. mil == 0) then ! l = 1, 5, 9, ... or 3, 7, 11, ...
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
#ifdef NEC_TUNE_SMP
!CDIR SELECT(CONCUR)
#endif
            do ib = 1, np_e
               bpr = 0.d0
               do i = 2, iba(ik)
                  bpr = bpr + ai(i)*psi_l(i,ib,ik,1) + ar(i)*psi_l(i,ib,ik,2)
               end do
               cr_l(ib,lmta1,1,ik) = fp*2.d0*bpr
            end do
         end if
      end if
    end subroutine phir_dot_WFs_core2
  end subroutine m_ES_phir_dot_WFs_4_lmta_k


  subroutine m_ES_PAO_WFs_3D(nfout)

    integer, intent(in)    :: nfout

    integer ia, ik, iksnl, mapmode
    integer     :: id_sname = -1

    call tstatc0_begin('m_ES_PAO_WFs ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_PAO_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call alloc_zfsincos_mpi()
    do ia = 1, natm
       if(kv3/nspin == 1) then
!!$          call G_dot_R_map(ia,1)
          mapmode = MAPPED
          call G_dot_R(ia,mapmode,1)
       else
!!$          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = NOTMAPPED
          call G_dot_R(ia,mapmode)
       end if
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI
          iksnl = (ik-1)/nspin + 1
          call m_ES_PAO_WFs_4_lmta_k_3D(ista_k,iend_k,ik,ia,iksnl,paog,zaj_l,mapmode)
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_zfsincos()
    call tstatc0_end(id_sname)
  end subroutine m_ES_PAO_WFs_3D

  subroutine m_ES_PAO_WFs_4_lmta_k_3D(k1,k2,ik,ia,iksnl,paog,psi_l,mapmode)

    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(maxval(np_g1k),nlmtt_pao,ista_snl:iend_snl) :: paog ! MPI
    real(kind=DP),intent(out), dimension(maxval(np_g1k),ista_e:iend_e,k1:k2,kimg) :: psi_l        ! MPI

    integer                    ::it, lmt1, lmtt1, lmta1, il1, ib, ib2

    it    = ityp(ia)
    do lmt1 = 1, ilmt_pao(it)
       lmtt1 = lmtt_pao(lmt1,it)
       lmta1 = lmta_pao(lmt1,ia)
       il1   = ltp_pao(lmt1,it)
       ib    = ibpao(lmta1)
       ib2   = 0
       if(ib < 0) then
          ib = abs(ib)
          ib2 = ib + 1
       end if
       if(ib >= ista_e .and. ib <= iend_e) then
          call G_dot_R_mult_paog(ib,.true.)             ! exp(iGR)
       end if
       if(ib2 >= ista_e .and. ib2 <= iend_e) then
          call G_dot_R_mult_paog(ib2,.false.)             ! exp(iGR)
       end if
    end do
  contains
    subroutine G_dot_R_mult_paog(ib,fcos)
      integer, intent(in) :: ib
      logical, intent(in) :: fcos
      integer :: i, i1, mil
      if(kimg==1) then
         if(fcos) then
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, np_g1k(ik)
                  psi_l(i,ib,ik,1) = zfcos(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, np_g1k(ik)
                  i1    = nbase(i+ista_g1k(ik)-1,ik)
                  psi_l(i,ib,ik,1) = zfcos(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, np_g1k(ik)
                  psi_l(i,ib,ik,1) = zfsin(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, np_g1k(ik)
                  i1    = nbase(i+ista_g1k(ik)-1,ik)
                  psi_l(i,ib,ik,1) = zfsin(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         end if
      else
         mil = mod(il1,4)
         if(mil == 1 .or. mil == 3) then
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, np_g1k(ik)
                  psi_l(i,ib,ik,1) = zfcos(i)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) = zfsin(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, np_g1k(ik)
                  i1    = nbase(i+ista_g1k(ik)-1,ik)
                  psi_l(i,ib,ik,1) = zfcos(i1)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) = zfsin(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, np_g1k(ik)
                  psi_l(i,ib,ik,1) = -zfsin(i)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) =  zfcos(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, np_g1k(ik)
                  i1    = nbase(i+ista_g1k(ik)-1,ik)
                  psi_l(i,ib,ik,1) = -zfsin(i1)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) =  zfcos(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         end if
      end if
    end subroutine G_dot_R_mult_paog

  end subroutine m_ES_PAO_WFs_4_lmta_k_3D

  subroutine wd_fsr_fsi_3D(ik)

    integer, intent(in) :: ik
    integer :: it, ia
    if(k_symmetry(ik) == GAMMA) then
       write(nfout,'(" -- fsr --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
! === DEBUG by tkato 2013/08/28 ================================================
!         write(nfout,'(6d12.4)') (fsr_l(map_z(it),ia,ik),ia=1,nlmta)!MPI
          write(nfout,'(6d12.4)') (fsr_l(map_z(it),ia,ik),ia=1,np_fs)!MPI
! ==============================================================================
       end do
    else
       write(nfout,'(" -- fsr, fsi --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
! === DEBUG by tkato 2013/08/28 ================================================
!         write(nfout,'(6d12.4)') (fsr_l(map_z(it),ia,ik),fsi_l(map_z(it),ia,ik),ia=1,nlmta)!MPI
          write(nfout,'(6d12.4)') (fsr_l(map_z(it),ia,ik),fsi_l(map_z(it),ia,ik),ia=1,np_fs)!MPI
! ==============================================================================
       end do
    end if
    call flush(nfout)
  end subroutine wd_fsr_fsi_3D

  subroutine m_ES_betar_dot_WFs_3D(nfout,ik)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, October 2009 : Mblock
!
    use m_Parallelization,only      : ball_buff, ball_addr

    integer, intent(in)    :: nfout, ik
    integer :: iksnl, ibsize
    integer :: datasize

    integer :: ia, it, lmt1, msize, msize_target, msizemax, msizesum, natm_redmax
    integer :: ibl1,ibl2, ia1, ia2
    integer :: LD11,  LD12
    logical :: tran1
!!$    integer, allocatable, dimension(:) :: mil
    real(kind=DP), allocatable, dimension(:,:,:) :: psi_ri
!f
!   real(kind=DP),allocatable, dimension(:,:) :: fsr_l_2D,fsi_l_2D
    real(kind=DP),allocatable, dimension(:,:) :: wk_fsr_l,wk_fsi_l
    real(kind=DP),allocatable, dimension(:,:) :: wk_fsr_ball,wk_fsi_ball
    integer :: i, j, k 

    integer :: id_sname = -1
! === rspace on 3D_Parallel ==============================================================
    logical :: flag_in_realspace
    call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)
  end subroutine m_ES_betar_dot_WFs_3D

  subroutine build_fft_map(nfft,kimg)
    integer, intent(in) :: nfft,kimg
    integer :: ia,nma,im
    if(.not.allocated(mapl)) allocate(mapl(nfft/2))
    mapl = .false.
    do ia=ista_atm,iend_atm
       nma = nmesh_rs(ia) 
       do im=1,nma
          mapl(meshxyz_rs(im,ia)) = .true.
       enddo
    enddo
    nonzero_fft_elements = 0
    do ia=ista_atm,iend_atm
       nma = nmesh_rs(ia) 
       do im=1,nma
          if(mapl(meshxyz_rs(im,ia))) nonzero_fft_elements = nonzero_fft_elements + 1
       enddo
    enddo
    write(0,*) 'nfft, nonzero_fft_elements ',nfft,nonzero_fft_elements
  end subroutine build_fft_map

  subroutine allreduce_wf(lsize,kimg,nfft,ibsize,psi_l,bff)
    integer, intent(in) :: lsize,kimg,nfft,ibsize
    real(kind=DP), dimension(lsize*kimg,ibsize), intent(in) :: psi_l
    real(kind=DP), dimension(nfft,ibsize), intent(out) :: bff
    integer :: imesh,i,j,ierr,ib
    integer :: id_sname = -1
    call tstatc0_begin('allreduce_wf ',id_sname)
    bff = 0.d0
    do ib=1,ibsize
    do imesh=1,lsize
       i = 2*dist_to_full(imesh)
       j = 2*imesh
       bff(i-1,ib) = psi_l(j-1,ib)
       bff(i,ib)   = psi_l(j,ib)
    enddo
    enddo
    call mpi_allreduce(mpi_in_place,bff,nfft*ibsize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    call tstatc0_end(id_sname)
  end subroutine allreduce_wf

  subroutine betar_dot_Psi_4_each_k_in_rs(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l)
    integer,intent(in) :: nfout 
    integer,intent(in) :: k1,k2,ik
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l
    real(kind=DP),intent(out),dimension(np_e,nlmta) :: bpr_l, bpi_l
    call betar_dot_Psi_4_each_k_in_rs0(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l,nmesh_rs_max,snl_rs)
  end subroutine betar_dot_Psi_4_each_k_in_rs

  subroutine betar_dot_Psi_4_each_k_in_rs0(nfout,k1,k2,ik,psi_l,bprt,bpit,nsnl,snlsnl)
    use m_Parallelization,    only : nel_fft_z, nel_fft_y, nel_fft_x,ista_lmta,iend_lmta
    integer,intent(in) :: nfout 
    integer,intent(in) :: k1,k2,ik
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l
    real(kind=DP),intent(out),dimension(np_e,nlmta) :: bprt, bpit
    integer, intent(in) :: nsnl
    real(kind=DP), intent(in), dimension(nsnl,ista_lmta:iend_lmta) :: snlsnl
    integer :: ia,ib
    real(kind=DP),allocatable,dimension(:,:) :: cos_kr,sin_kr
#ifdef MULT_PHASE_RSPACE
    real(kind=DP),allocatable,dimension(:) :: cos_a,sin_a
#endif
    real(kind=DP), allocatable, dimension(:,:) :: bff
    real(kind=DP), allocatable, dimension(:,:) :: wk_bfft_l
    real(kind=DP), allocatable, dimension(:,:) :: rrr, iii
    integer :: id_sname = -1
    integer :: ib1, ib2, ibsize, isrsize, lsize, iesize, ierr,i,j
    logical, save :: initialized = .false.
    call tstatc0_begin('betar_dot_Psi_4_each_k_in_rs ',id_sname,level=1)
    if(.not.initialized)then
    call build_fft_map(nfft,kimg)
    initialized = .true.
    endif
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
!    if(.not.initialized) then
!       call build_fft_map(lsize,nfft,kimg)
!       initialized = .true.
!    endif
    ibsize = nblocksize_rspace_betar
    if(ibsize<1) ibsize = 1
    if(ibsize>np_e) ibsize = np_e
!    if (nblocksize_fftw_is_given) then
!       ibsize = nblocksize_fftw
!       if(ibsize < 1) ibsize = 1
!    endif

    call alloc_work_arrays()
    call k_dot_r(ik,cos_kr,sin_kr)
#ifdef MULT_PHASE_RSPACE
    call k_dot_pos(ik,cos_a,sin_a)
#endif
    do ib1=1,np_e, ibsize
       !ib2 = min(ib1+ibsize-1,np_e)
       ib2 = min(ibsize,np_e+1-ib1)
       call m_ES_WF_2D_psi(ik,wk_bfft_l,psi_l,min(ibsize+ib1-1,np_e),ib1,k1,k2,ibsize,lsize,INVERSE,bff)
!       call m_ES_WF_in_Rspace_3D(k1,k2,ik,ib1,ib2,ibsize,lsize,psi_l,wk_bfft_l)
!       call allreduce_wf(lsize,kimg,nfft,ibsize,wk_bfft_l,bff)
       do ia=ista_atm,iend_atm
          call betar_dot_Psi_atm_band_pre(bff,ia)
          call betar_dot_Psi_atm_band(ia)
       enddo
    enddo
    
    if(nrank_g>1)then
       call mpi_allreduce(mpi_in_place,bprt,np_e*nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       if(k_symmetry(ik) /= GAMMA) then
          call mpi_allreduce(mpi_in_place,bpit,np_e*nlmta,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       endif
    endif

    call dealloc_work_arrays()
    call tstatc0_end(id_sname)

    !write(nfout,'(a)') 'fsr fsr '
    !do i=1,np_e
    !   write(nfout,'(a,i8)') 'band ',i
    !   write(nfout,'(5f20.10)') (bprt(i,j),j=1,nlmta)
    !enddo
    contains

    subroutine alloc_work_arrays()
       allocate(wk_bfft_l(lsize*kimg,ibsize))
       allocate(bff(nfft,ibsize));bff=0.d0
       allocate(cos_kr(nmesh_rs_max,natm));cos_kr=0.d0
       allocate(sin_kr(nmesh_rs_max,natm));sin_kr=0.d0
#ifdef MULT_PHASE_RSPACE
       allocate(cos_a(natm));cos_a=0.d0
       allocate(sin_a(natm));sin_a=0.d0
#endif
       allocate(rrr(nmesh_rs_max,ibsize)); rrr = 0.0d0
       allocate(iii(nmesh_rs_max,ibsize)); iii = 0.0d0
    end subroutine alloc_work_arrays

    subroutine dealloc_work_arrays()
       deallocate(wk_bfft_l)
       deallocate(bff)
       deallocate(cos_kr)
       deallocate(sin_kr)
#ifdef MULT_PHASE_RSPACE
       deallocate(cos_a)
       deallocate(sin_a)
#endif
       deallocate(rrr)
       deallocate(iii)
    end subroutine dealloc_work_arrays

    subroutine betar_dot_Psi_atm_band_pre(bff,iatm)
       real(kind=DP),dimension(nfft,ibsize),intent(in) :: bff
       integer,intent(in) :: iatm
       integer :: lmt1,i,lmta1,imesh,nma,it,ib
       real(kind=DP) :: s,bpr,bpi
       real(kind=DP) :: fac,cosa,sina,rr,ii,psr,psi
       real(kind=DP), allocatable, dimension(:) :: cos_krt,sin_krt
#ifdef USE_DDOT
       real(kind=DP) :: ddot
#endif
       integer :: id_sname=-1
       call tstatc0_begin('betar_dot_Psi_atm_band_pre ',id_sname,level=1)
       fac = dsqrt(univol)/dble(fft_box_size_WF(1,1)*fft_box_size_WF(2,1)*fft_box_size_WF(3,1))
       nma = nmesh_rs(iatm)
       it = ityp(iatm)
#ifdef MULT_PHASE_RSPACE
       cosa = cos_a(iatm)
       sina = sin_a(iatm)
#endif
       allocate(cos_krt(nma));cos_krt(:) = cos_kr(:,iatm)
       allocate(sin_krt(nma));sin_krt(:) = sin_kr(:,iatm)
!       do ib=ib1,ib2
       do ib=1,ib2
       if(kimg==1)then
          do i=1,nma
             imesh   = 2*meshxyz_rs(i,iatm)
#ifdef MULT_PHASE_RSPACE
             psr = bff(imesh-1,ib)
             psi = bff(imesh,ib)*meshxyz_rs_conjg(i,iatm)
             rr = psr*cosa-psi*sina
             ii = psi*cosa+psr*sina
#else
             rr = bff(imesh-1,ib)
             ii = bff(imesh,ib)*meshxyz_rs_conjg(i,iatm)
#endif
             rrr(i,ib) = rr*cos_krt(i)-ii*sin_krt(i) 
             iii(i,ib) = ii*cos_krt(i)+rr*sin_krt(i)
          enddo
       else
          do i=1,nma
             imesh   = 2*meshxyz_rs(i,iatm)
#ifdef MULT_PHASE_RSPACE
             psr = bff(imesh-1,ib)
             psi = bff(imesh,ib)
             rr = psr*cosa-psi*sina
             ii = psi*cosa+psr*sina
#else
             rr = bff(imesh-1,ib)
             ii = bff(imesh,ib)
#endif
             rrr(i,ib) = rr*cos_krt(i)-ii*sin_krt(i) 
             iii(i,ib) = ii*cos_krt(i)+rr*sin_krt(i)
          enddo
       endif
       enddo
       deallocate(cos_krt,sin_krt)
       call tstatc0_end(id_sname)
    end subroutine betar_dot_Psi_atm_band_pre

    subroutine betar_dot_Psi_atm_band(iatm)
       integer,intent(in) :: iatm
       integer :: lmt1,i,lmta1,imesh,nma,it,ib,iband
       real(kind=DP) :: s,bpr,bpi
       real(kind=DP) :: fac,cosa,sina,rr,ii,psr,psi
#ifdef USE_DDOT
       real(kind=DP) :: ddot
#endif
       integer :: id_sname=-1
       call tstatc0_begin('betar_dot_Psi_atm_band ',id_sname,level=1)
       fac = dsqrt(univol)/dble(fft_box_size_WF(1,1)*fft_box_size_WF(2,1)*fft_box_size_WF(3,1))
       nma = nmesh_rs(iatm)
       it = ityp(iatm)
#ifdef NO_DGEMM_RSPACE
       do ib=1,ib2
          iband = ib+ib1-1
          do lmt1=1,ilmt(it)
             lmta1 = lmta(lmt1,iatm)
#ifdef USE_DDOT
             bprt(iband,lmta1) = fac*ddot(nma,snlsnl(1:nma,lmta1),1,rrr(1:nma,ib),1)
#else
             bprt(iband,lmta1) = fac*dot_product(snlsnl(1:nma,lmta1),rrr(1:nma,ib))
#endif
          enddo
       enddo
#else
       call dgemm('T','N',ib2,ilmt(it),nma,fac,rrr,nmesh_rs_max,snlsnl(1:nma,lmta(1,iatm)) &
                     ,nmesh_rs_max,0.d0,bprt(ib1,lmta(1,iatm)),np_e)
#endif
       if(k_symmetry(ik) /= GAMMA)then
#ifdef NO_DGEMM_RSPACE
       do ib=1,ib2
          iband = ib+ib1-1
          do lmt1=1,ilmt(it)
             lmta1 = lmta(lmt1,iatm)
#ifdef USE_DDOT
             bpit(iband,lmta1) = fac*ddot(nma,snlsnl(1:nma,lmta1),1,iii(1:nma,ib),1)
#else
             bpit(iband,lmta1) = fac*dot_product(snlsnl(1:nma,lmta1),iii(1:nma,ib))
#endif
          enddo
       enddo
#else
       call dgemm('T','N',ib2,ilmt(it),nma,fac,iii,nmesh_rs_max,snlsnl(1:nma,lmta(1,iatm)) &
                     ,nmesh_rs_max,0.d0,bpit(ib1,lmta(1,iatm)),np_e)
#endif
       endif
       call tstatc0_end(id_sname)
    end subroutine betar_dot_Psi_atm_band

  end subroutine betar_dot_Psi_4_each_k_in_rs0

  subroutine Vnonlocal_W_in_realspace(ik,ispin,switch_of_eko_part)
    use m_Electronic_Structure, only : m_ES_gather_f_3d_to_2d
    integer, intent(in) :: ik,ispin
    integer, intent(in) :: switch_of_eko_part
    integer :: mdvdb, it, ia, lmt2, lmta2, il2, im2,nma
    real(kind=DP),allocatable,dimension(:,:) :: vnlr1d,vnli1d
    real(kind=DP),allocatable,dimension(:) :: sc,qc
    real(kind=DP),allocatable,dimension(:,:) :: sc_lmta,qc_lmta
    real(kind=DP),allocatable,dimension(:,:) :: cos_kr,sin_kr
#ifdef MULT_PHASE_RSPACE
    real(kind=DP),allocatable,dimension(:) :: cos_a,sin_a
#endif
    real(kind=DP),allocatable,dimension(:,:) :: bff
    real(kind=DP),allocatable,dimension(:,:) :: bffb,fsrt,fsit
    integer :: ib,ierr,ii,ib2
    integer :: ibl1,ibl2,nbsize,ibsize
    integer :: id_sname = -1,id_sname2=-1,id_sname3=-1
    call tstatc0_begin('Vnonlocal_W_in_realspace ',id_sname,level=1)
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       if(mdvdb == EXECUT) exit
    end do
    
    ibsize = nblocksize_rspace_v
    if(ibsize<1) ibsize = 1
    if(ibsize>np_e) ibsize = np_e
    call alloc_arrays()
    call m_ES_gather_f_3d_to_2d(fsr_l, fsrt, ik)
    if(k_symmetry(ik) /= GAMMA) then
      call m_ES_gather_f_3d_to_2d(fsi_l, fsit, ik)
    endif
    call k_dot_r(ik,cos_kr,sin_kr)
#ifdef MULT_PHASE_RSPACE
    call k_dot_pos(ik,cos_a,sin_a)
#endif
    Loop_natm0: do ia=ista_atm,iend_atm
       if(kimg==1)then
          nma = nmesh_rs_h(ia)
       else
          nma = nmesh_rs(ia)
       endif
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       do lmt2 = 1, ilmt(it)
          lmta2 = lmta(lmt2,ia)
          il2   = ltp(lmt2,it)
          im2   = mtp(lmt2,it)
          call part_sum_over_lmt1_rs()
       enddo
    end do Loop_natm0

    do ib=1,np_e,ibsize
       bff=0.d0
       ib2 = min(ibsize,np_e+1-ib)
       Loop_natm1: do ia=ista_atm,iend_atm
          it = ityp(ia)
          mdvdb = m_PP_include_vanderbilt_pot(it)
          if(switch_of_eko_part == OFF) mdvdb= SKIP
          vnlr1d = 0.d0
          vnli1d = 0.d0
          if(kimg==1)then
             nma = nmesh_rs_h(ia)
          else
             nma = nmesh_rs(ia)
          endif
          if(mdvdb==SKIP)then
             call add_vnlph_l_without_eko_part_rs()
          else
             call add_vnlph_l_with_eko_part_rs()
          endif
          call multiply_phase()
          call map_vnl_to_bff()
       end do Loop_natm1
       if(nrank_g>1)then
         call mpi_allreduce(mpi_in_place,bff,nfft*ibsize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       endif
       call Vnonlocal_W_to_Gspace()
    enddo
    call dealloc_arrays()

    call tstatc0_end(id_sname)

    contains

    subroutine alloc_arrays()
       integer :: nmm
       nmm = nmesh_rs_max
       if(kimg==1) nmm = nmesh_rs_max_h
       allocate(vnlr1d(nmm,ibsize));vnlr1d=0.d0
       allocate(vnli1d(nmm,ibsize));vnli1d=0.d0
       allocate(sc_lmta(nmm,nlmta));sc_lmta=0.d0
       if(mdvdb==EXECUT)then
          allocate(qc_lmta(nmm,nlmta));qc_lmta=0.d0
       endif
       allocate(cos_kr(nmesh_rs_max,natm));cos_kr=0.d0
       allocate(sin_kr(nmesh_rs_max,natm));sin_kr=0.d0
#ifdef MULT_PHASE_RSPACE
       allocate(cos_a(natm));cos_a=0.d0
       allocate(sin_a(natm));sin_a=0.d0
#endif
       if(sw_save_memory==ON)then
          allocate(bff(nfft,ibsize));bff=0.d0
       else
          allocate(bffb(nfft,np_e));bffb=0.d0
       endif
       allocate(fsrt(np_e,nlmta));fsrt=0.d0
       allocate(fsit(np_e,nlmta));fsit=0.d0
    end subroutine alloc_arrays

    subroutine dealloc_arrays()
       deallocate(vnlr1d)
       deallocate(vnli1d)
       deallocate(sc_lmta)
       if(allocated(qc_lmta)) deallocate(qc_lmta)
       deallocate(bff)
       deallocate(cos_kr)
       deallocate(sin_kr)
#ifdef MULT_PHASE_RSPACE
       deallocate(cos_a)
       deallocate(sin_a)
#endif
       deallocate(fsrt)
       deallocate(fsit)
    end subroutine dealloc_arrays

    subroutine part_sum_over_lmt1_rs()
       integer       :: lmta1,lmt1, il1, im1, i
       real(kind=DP) :: tmp,rr,ii
       real(kind=DP),pointer,dimension(:,:) :: snlt
       integer :: id_sname = -1
       call tstatc0_begin('part_sum_over_lmt1_rs ',id_sname,level=1)
       if(sw_save_memory/=ON)then
          sc = 0.d0
          if(mdvdb == EXECUT) then
             qc = 0.d0
          endif
       else
          sc_lmta(:,lmta2) = 0.d0
          if(mdvdb == EXECUT) then
             qc_lmta(:,lmta2) = 0.d0
          endif
       endif
       if(kimg==1)then
          snlt => snl_rs_h
       else
          snlt => snl_rs
       endif
       do lmt1 = 1,ilmt(it)
          lmta1 = lmta(lmt1,ia)
          il1   = ltp(lmt1,it)
          im1   = mtp(lmt1,it)
          if(il1 == il2 .and. im1 == im2) then
             if(ipaw(it)==0) then
                 tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
             else
                 tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
             end if
          else
             if(ipaw(it)==0) then
                 tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
             else
                 tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
             end if
          endif
          do i=1,nma
             sc_lmta(i,lmta2) = sc_lmta(i,lmta2)+tmp*snlt(i,lmta1)
          enddo
          if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
             tmp = q(lmt1,lmt2,it)
             do i=1,nma
                qc_lmta(i,lmta2) = qc_lmta(i,lmta2) + tmp*snlt(i,lmta1)
             enddo
          end if
       end do

       call tstatc0_end(id_sname)
    end subroutine part_sum_over_lmt1_rs

    subroutine add_vnlph_l_with_eko_part_rs()
       integer :: i,iband,ib1,nmm
       real(kind=DP) :: fr,fi,e
       real(kind=DP), allocatable, dimension(:,:) :: fre,fie
       integer :: id_sname = -1
       call tstatc0_begin('add_vnlph_l_with_eko_part_rs ',id_sname,level=1)
#ifdef NO_DGEMM_RSPACE
       do ib1=1,ib2
       iband = ib+ib1-1
       do lmt2 = 1, ilmt(it)
          lmta2 = lmta(lmt2,ia)
          fr = fsrt(iband,lmta2)
          e = eko_l(iband,ik)
          do i=1,nma
             vnlr1d(i,ib1) = vnlr1d(i,ib1)+fr*(sc_lmta(i,lmta2)-e*qc_lmta(i,lmta2))
          enddo
          if(k_symmetry(ik)/=GAMMA) then
             fi = fsit(iband,lmta2)
             do i=1,nma
                vnli1d(i,ib1) = vnli1d(i,ib1)+fi*(sc_lmta(i,lmta2)-e*qc_lmta(i,lmta2))
             enddo
          endif
       enddo
       enddo
#else
       nmm = nmesh_rs_max
       if(kimg==1) nmm = nmesh_rs_max_h
       allocate(fre(ib2,ilmt(it)));fre=0.d0
       if(k_symmetry(ik) /= GAMMA) allocate(fie(ib2,ilmt(it)));fie=0.d0
       do ib1=1,ib2
          iband = ib+ib1-1
          do lmt2 = 1, ilmt(it)
             lmta2 = lmta(lmt2,ia)
             fre(ib1,lmt2) = -fsrt(iband,lmta2)*eko_l(iband,ik)
             if(k_symmetry(ik) /= GAMMA) fie(ib1,lmt2) = -fsit(iband,lmta2)*eko_l(iband,ik)
          enddo
       enddo
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsrt(ib,lmta(1,ia)),np_e,0.d0,vnlr1d,nmm)
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,qc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fre,ib2,1.d0,vnlr1d,nmm)
       deallocate(fre)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsit(ib,lmta(1,ia)),np_e,0.d0,vnli1d,nmm)
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,qc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fie,ib2,1.d0,vnli1d,nmm)
          deallocate(fie)
       endif
#endif
       call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_with_eko_part_rs

    subroutine add_vnlph_l_without_eko_part_rs()
       real(kind=DP) :: fr,fi
       integer :: i,iband,ib1,nmm
       integer :: id_sname = -1
       call tstatc0_begin('add_vnlph_l_without_eko_part_rs ',id_sname,level=1)
#ifdef NO_DGEMM_RSPACE
       do ib1=1,ib2
       iband = ib1+ib-1
       do lmt2 = 1, ilmt(it)
          lmta2 = lmta(lmt2,ia)
          fr = fsrt(iband,lmta2)
          do i=1,nma
             vnlr1d(i,ib1) = vnlr1d(i,ib1)+fr*sc_lmta(i,lmta2)
          enddo
          if(k_symmetry(ik)/=GAMMA) then
             fi = fsit(iband,lmta2)
             do i=1,nma
                vnli1d(i,ib1) = vnli1d(i,ib1)+fi*sc_lmta(i,lmta2)
             enddo
          endif
       enddo
       enddo
#else
       nmm = nmesh_rs_max
       if(kimg==1) nmm = nmesh_rs_max_h
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsrt(ib,lmta(1,ia)),np_e,0.d0,vnlr1d,nmm)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsit(ib,lmta(1,ia)),np_e,0.d0,vnli1d,nmm)
       endif
#endif
       call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_without_eko_part_rs

    subroutine multiply_phase()
       integer :: i,iband,ib1
       real(kind=DP) :: rr,ii
       real(kind=DP) :: cosa,sina,fac
       integer :: id_sname = -1
       call tstatc0_begin('multiply_phase ',id_sname,level=1)
       fac = dsqrt(univol)/dble(fft_box_size_WF(1,1)*fft_box_size_WF(2,1)*fft_box_size_WF(3,1))
#ifdef MULT_PHASE_RSPACE
       cosa =  cos_a(ia)
       sina = -sin_a(ia)
#endif

       if(kimg==1)then
          do ib1=1,ib2
          do i=1,nma
#ifdef MULT_PHASE_RSPACE
             rr = cosa*vnlr1d(i,ib1)-sina*vnli1d(i,ib1)
             ii = sina*vnlr1d(i,ib1)+cosa*vnli1d(i,ib1)
#else
             rr = vnlr1d(i,ib1)
             ii = vnli1d(i,ib1)
#endif
             vnlr1d(i,ib1) = (rr*cos_kr(map_h(i,ia),ia)+ii*sin_kr(map_h(i,ia),ia))*fac
             vnli1d(i,ib1) = (ii*cos_kr(map_h(i,ia),ia)-rr*sin_kr(map_h(i,ia),ia))*fac
          enddo
          enddo
       else
          do ib1=1,ib2
          do i=1,nma
#ifdef MULT_PHASE_RSPACE
             rr = cosa*vnlr1d(i,ib1)-sina*vnli1d(i,ib1)
             ii = sina*vnlr1d(i,ib1)+cosa*vnli1d(i,ib1)
#else
             rr = vnlr1d(i,ib1)
             ii = vnli1d(i,ib1)
#endif
             vnlr1d(i,ib1) = (rr*cos_kr(i,ia)+ii*sin_kr(i,ia))*fac
             vnli1d(i,ib1) = (ii*cos_kr(i,ia)-rr*sin_kr(i,ia))*fac
          enddo
          enddo
       endif
       call tstatc0_end(id_sname)
    end subroutine multiply_phase

    subroutine map_vnl_to_bff()
       integer :: i,i1,iband,ib1
       integer :: id_sname=-1
       call tstatc0_begin('map_vnl_to_bff ',id_sname,level=1)
       if(kimg==1)then 
          do ib1=1,ib2
          do i=1,nma
             i1 = 2*meshxyz_rs_h(i,ia)
             bff(i1-1,ib1) = bff(i1-1,ib1) + vnlr1d(i,ib1)
             bff(i1,ib1)   = bff(i1,ib1)   - vnli1d(i,ib1)
          enddo
          enddo
       else
          do ib1=1,ib2
          do i=1,nma
             i1 = 2*meshxyz_rs(i,ia)
             bff(i1-1,ib1) = bff(i1-1,ib1) + vnlr1d(i,ib1)
             bff(i1,ib1)   = bff(i1,ib1)   + vnli1d(i,ib1)
          enddo
          enddo
       endif
       call tstatc0_end(id_sname)
    end subroutine map_vnl_to_bff

    subroutine Vnonlocal_W_to_Gspace()
       integer :: i,i1,ri,iband,ib1
       integer :: id_sname = -1
       call tstatc0_begin('Vnonlocal_W_to_Gspace ',id_sname)
       do ib1=1,ib2
       iband = ib1+ib-1
       vnlph_l(:,iband,:) = 0.d0
       call m_FFT_WF(ELECTRON,nfout,bff(1:nfft,ib1),DIRECT,ON)
       do ri=1,kimg
          do i=ista_g1k(ik),iend_g1k(ik)
             i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
             vnlph_l(i-ista_g1k(ik)+1,iband,ri) = bff(i1,ib1)
          enddo
       enddo
       enddo
       call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_to_Gspace

  end subroutine Vnonlocal_W_in_realspace

  subroutine k_dot_r(ik,zc_ar,zs_ar)
    integer, intent(in) :: ik
    real(kind=DP),dimension(nmesh_rs_max,natm),intent(out) :: zc_ar,zs_ar
    integer :: ia,i
    real(kind=DP) :: inl,inm,inn
    real(kind=DP) :: rx,ry,rz,kdr
    integer :: id_sname = -1
    call tstatc0_begin('k_dot_r ',id_sname,level=1)
    inl = 1.d0/dble(fft_box_size_WF(1,1))
    inm = 1.d0/dble(fft_box_size_WF(2,1))
    inn = 1.d0/dble(fft_box_size_WF(3,1))
    do ia=ista_atm,iend_atm
       do i=1,nmesh_rs(ia)
          rx = dble(meshx_rs(i,ia))*inl
          ry = dble(meshy_rs(i,ia))*inm
          rz = dble(meshz_rs(i,ia))*inn
          kdr = (rx*vkxyz(ik,1,BUCS)+ry*vkxyz(ik,2,BUCS)+rz*vkxyz(ik,3,BUCS))*PAI2
          zc_ar(i,ia) = dcos(kdr)
          zs_ar(i,ia) = dsin(kdr)
       enddo
    enddo
  end subroutine k_dot_r

#ifdef MULT_PHASE_RSPACE
  subroutine k_dot_pos(ik,zc,zs)
    integer, intent(in) :: ik
    real(kind=DP),dimension(natm),intent(out) :: zc,zs
    real(kind=DP) :: kdp
    integer :: ia
    do ia=1,natm
       kdp = (pos(ia,1)*vkxyz(ik,1,BUCS)+pos(ia,2)*vkxyz(ik,2,BUCS)+pos(ia,3)*vkxyz(ik,3,BUCS))*PAI2
       zc(ia) =  dcos(kdp)
       zs(ia) = -dsin(kdp)
    enddo
  end subroutine k_dot_pos

#endif

end module m_ES_nonlocal
