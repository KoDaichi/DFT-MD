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
       &                         , irorb, nrorb, crorb, nlmt
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
       &                           , neordr, nrvf_ordr

  use m_Electronic_Structure,only :  afft,bfft 

  use m_Electronic_Structure, only : m_ES_WF_in_Rspace

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

#ifdef _CUDA_
  use cublas, only : cublasDgemm
#endif
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

  integer, private, parameter                         :: sw_timing_2ndlevel = ON

#ifdef NONLOCAL_DGEMM
  real(kind=DP),        allocatable,dimension(:,:)    :: wk_zfcos, wk_zfsin  !d(nbmx_ext)
!!$  integer      ,allocatable,dimension(:)              :: lmtt_tmp
  real(kind=DP),private,allocatable,dimension(:,:)    :: wk_ar,  wk_ai
  real(kind=DP),allocatable,dimension(:,:)            :: bp_tmp1,  bp_tmp2
!!$  real(kind=DP),allocatable,dimension(:)              :: ia_tmp
  logical :: DGEMM_DEBUG = .false.
#endif


#ifdef SX
  integer,private                                     :: nb_vnonlocal_default = 5000  ! TY 26Aug2009
  integer,private                                     :: nb_betar_default     = 10000 ! TY
#else
  integer,private                                     :: nb_vnonlocal_default = 1000  ! TY 26Aug2009
  integer,private                                     :: nb_betar_default     = 32    ! TY
#endif

  integer,private                                     :: nlmta_us
!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  contains

  subroutine m_ES_alloc_scss_etc()
    call m_FFT_alloc_WF_work()  ! allocate(ftw)
    allocate(sc(kg1))
    allocate(ss(kg1))
    allocate(qc(kg1))
    allocate(qs(kg1))
  end subroutine m_ES_alloc_scss_etc


  subroutine m_ES_dealloc_scss_etc()
    if(allocated(qs)) deallocate(qs)
    if(allocated(qc)) deallocate(qc)
    if(allocated(ss)) deallocate(ss)
    if(allocated(sc)) deallocate(sc)
    call m_FFT_dealloc_WF_work()
  end subroutine m_ES_dealloc_scss_etc

  subroutine m_ES_alloc_afft_scss_etc()
    allocate(afft(nfft))
    allocate(bfft(nfft))
    call m_ES_alloc_scss_etc()
  end subroutine m_ES_alloc_afft_scss_etc


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

  subroutine alloc_arai(ibsize)
    integer, intent(in) :: ibsize
    allocate(ar(ibsize),ai(ibsize))
  end subroutine alloc_arai

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

  subroutine m_ES_alloc_arai(ik)
    integer, intent(in) :: ik
    if(ik == 0) then
! ====================================== modified by K. Tagami ============ 11.0
!       if(kv3/nspin == 1) then
!          allocate(ar(iba(1)))
!          allocate(ai(iba(1)))
!       else
!          allocate(ar(nbmx))
!          allocate(ai(nbmx))
!       end if

      if ( noncol ) then
        if ( kv3/ndim_spinor == 1 ) then
           allocate(ar(iba(1)));       allocate(ai(iba(1)))
        else
           allocate(ar(nbmx));         allocate(ai(nbmx))
        end if
      else
        if ( kv3/nspin == 1 ) then
           allocate(ar(iba(1)));       allocate(ai(iba(1)))
        else
           allocate(ar(nbmx));         allocate(ai(nbmx))
        end if
      endif
! ======================================================================== 11.0
    else
       allocate(ar(iba(ik)))
       allocate(ai(iba(ik)))
    end if
! ==================================== added by K. Tagami =============== 11.0
    ar = 0.0d0;  ai = 0.0d0
! ======================================================================== 11.0
  end subroutine m_ES_alloc_arai

  subroutine m_ES_dealloc_arai()
    deallocate(ar,ai)
  end subroutine m_ES_dealloc_arai

#ifdef NONLOCAL_DGEMM
  subroutine alloc_wkzfsincos(ibsize)
    integer, intent(in) :: ibsize
    integer             :: ichkalloc
    allocate(wk_zfsin(ibsize,natm), stat=ichkalloc)
    allocate(wk_zfcos(ibsize,natm), stat=ichkalloc)
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos', ibsize, natm
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
    endif
  end subroutine alloc_wkzfsincos

  subroutine alloc_wkzfsincos_red(ibsize,natmsize)
    integer, intent(in) :: ibsize, natmsize
    integer             :: ichkalloc
    allocate(wk_zfsin(ibsize,natmsize), stat=ichkalloc)
    allocate(wk_zfcos(ibsize,natmsize), stat=ichkalloc)
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate wk_zfsin or wk_zfcos in alloc_wkzfsincos_red', ibsize, natm
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
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
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
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
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
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
  subroutine m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part,map)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part
    logical, intent(in), optional, dimension(np_e) :: map

    integer :: mdvdb, it, ia, lmt2
    integer :: ibl1,ibl2,ibsize
    integer :: ichkalloc
    integer :: id_sname = -1
    real(kind=DP),allocatable,dimension(:,:) :: wk_sc_without,  wk_ss_without
    real(kind=DP),allocatable,dimension(:,:) :: wk_sc_with,     wk_ss_with,     wk_qc_with,     wk_qs_with
    real(kind=DP),allocatable,dimension(:,:) :: fsr_tmp_without, fsi_tmp_without
    real(kind=DP),allocatable,dimension(:,:) :: fsr_tmp_with   , fsi_tmp_with
#ifdef _CUDA_
    real(kind=DP),allocatable,dimension(:,:) :: fsr_tmp_with_e   , fsi_tmp_with_e
#endif
    integer,allocatable,dimension(:) :: lmt2_tmp_without, ia_tmp_without, lmta2_tmp_without
    integer,allocatable,dimension(:) :: lmt2_tmp_with,    ia_tmp_with,    lmta2_tmp_with
    integer,allocatable,dimension(:) :: il2_tmp_without,  im2_tmp_without
    integer,allocatable,dimension(:) :: il2_tmp_with,     im2_tmp_with
    integer :: icnt_without  , icnt_with
    integer,allocatable,dimension(:) :: mdvdba
    logical :: mdvdb_EXECUT = .false.
    integer, save :: ibsize_print = OFF

! ================================== added by K. Tagami =============== 11.0
    integer :: ik_for_pointing_eko
! ===================================================================== 11.0
    integer :: ic, ib
! ----------- Revised by T.Yamasaki, 1 Aug. 2014 ---------->>
    logical :: flag_in_realspace

    flag_in_realspace = .false.
    if(sw_rspace == ON .and. sw_rspace_v == ON) then
       if(m_CtrlP_realspace_integ_OK()) flag_in_realspace = .true.
    end if

    if(ipribetar>=2 .and. sw_rspace == ON .and. .not.flag_in_realspace) write(nfout,&
         & '(" sw_rspace = ON, but the realspace integration is not installed in the solver now applied. ")')

!!$    if(sw_rspace==ON)then
    if(flag_in_realspace)then
! <<---------------------------------------------------------
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
    call tstatc0_begin('m_ES_Vnonlocal_W ',id_sname,level=1)
! --> T. Yamasaki, 26th Aug. 2009
!!$    ibsize=5000
    if(nblocksize_vnonlocal_is_given) then
       ibsize = nblocksize_vnonlocal_w
    else
       ibsize= nb_vnonlocal_default
       if(ipribetar >= 2 .and. ibsize_print == OFF) then
          write(nfout,'(" !ibsize(=nblocksize_vnonlocal_w) (m_ES_Vnonlocal_W) = ",i8)') ibsize
          ibsize_print = ON
       end if
    end if
! <--
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
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
    endif

    call pre_m_ES_Vnonlocal_W_count( icnt_without, icnt_with )
    allocate(mdvdba(ntyp))
    do it = 1, ntyp
       mdvdba(it) = m_PP_include_vanderbilt_pot(it)
    enddo
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       if(mdvdb == EXECUT) then
         mdvdb_EXECUT = .true.
         exit
       endif
    end do
    call alloc_wkzfsincos(ibsize)
    allocate( wk_sc_without(ibsize,icnt_without), stat=ichkalloc ) ;   wk_sc_without=0.d0
    allocate( wk_ss_without(ibsize,icnt_without), stat=ichkalloc ) ;   wk_ss_without=0.d0
    allocate( fsr_tmp_without(np_e,icnt_without), stat=ichkalloc ) ; fsr_tmp_without=0.d0
    allocate( wk_sc_with(ibsize,icnt_with),       stat=ichkalloc ) ;      wk_sc_with=0.d0
    allocate( wk_ss_with(ibsize,icnt_with),       stat=ichkalloc ) ;      wk_ss_with=0.d0
    allocate( fsr_tmp_with(np_e,icnt_with),       stat=ichkalloc ) ;    fsr_tmp_with=0.d0
#ifdef _CUDA_
    allocate( fsr_tmp_with_e(np_e,icnt_with), stat=ichkalloc )
    allocate( fsi_tmp_with_e(np_e,icnt_with), stat=ichkalloc )
#endif
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate work-array2 in m_ES_Vnonlocal_W', ibsize, icnt_without, icnt_with, np_e
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
    endif

    if( k_symmetry(ik) /= GAMMA ) then
      allocate( fsi_tmp_without(np_e,icnt_without), stat=ichkalloc ) ; fsi_tmp_without=0.d0
      allocate( fsi_tmp_with(np_e,icnt_with)      , stat=ichkalloc ) ;    fsi_tmp_with=0.d0
      if( ichkalloc /= 0 ) then
        write(nfout,*) 'could not allocate work-array3 in m_ES_Vnonlocal_W', np_e,icnt_without,icnt_with
        call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
      endif
    endif

    if(mdvdb_EXECUT) then
      allocate( wk_qc_with(ibsize,icnt_with), stat=ichkalloc ) ; wk_qc_with=0.d0
      allocate( wk_qs_with(ibsize,icnt_with), stat=ichkalloc ) ; wk_qs_with=0.d0
      if( ichkalloc /= 0 ) then
        write(nfout,*) 'could not allocate work-array4 in m_ES_Vnonlocal_W', icnt_without, icnt_with, natm
        call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
      endif
    endif
!
!
! Revised by T. Yamasaki, September 2008
!
     if(DGEMM_DEBUG) write(nfout,*)' DGEMM_debug Vnonlocal :', icnt_without,icnt_with, nlmta, ibsize, kimg, k_symmetry(ik), GAMMA
   if(ibsize_print == OFF .and. ipribetar >= 2) then
        write(nfout,'(" mype = ",i3, " ibsize, np_e, icnt_with, icnt_without, iba(ik)/ibsize = ",5i8," <<m_ES_Vnonlocal_W>>")') &
             & mype, ibsize, np_e, icnt_with, icnt_without, iba(ik)/ibsize
        ibsize_print = ON
     end if


! ======================================== added by K. Tagami ========= 11.0
    if ( noncol ) then
       ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1
    else
       ik_for_pointing_eko = ik
    endif
! ====================================================================== 11.0

! ======================================= modified by K. Tagami ======== 11.0
!     do ibl1=1,iba(ik),ibsize
!        ibl2=min( ibl1+ibsize-1,iba(ik) )
!        call calc_phase_blk(ik,ibl1,ibl2)
!        if( icnt_with .gt. 0 ) then
!          call Vnonlocal_W_part_with_blk   ( ibl1,ibl2,icnt_with   )
!          call add_vnlph_l_with_eko_blk   (  ibsize,ibl1,ibl2,icnt_with,   vnlph_l )
!        endif
!        if( icnt_without .gt. 0 ) then
!          call Vnonlocal_W_part_without_blk( ibl1,ibl2,icnt_without)
!          call add_vnlph_l_without_eko_blk(  ibsize,ibl1,ibl2,icnt_without,vnlph_l )
!        endif
!     enddo
!
    if( ichkalloc /= 0 ) then
      write(nfout,*) 'could not allocate work-array in m_ES_Vnonlocal_W', nlmta
      call phase_error_with_msg(nfout,'allocation error at m_ES_nonlocal',__LINE__,__FILE__)
    endif
    !$acc data copyout(vnlph_l) &
    !$acc copyin(fsr_l, fsi_l, vlhxcQ, eko_l, k_symmetry, iwei, mdvdba, pos) &
    !$acc create(fsr_tmp_without, fsi_tmp_without, wk_sc_without, wk_ss_without, &
    !$acc        fsr_tmp_with, fsi_tmp_with, wk_sc_with, wk_ss_with, wk_qc_with, &
    !$acc        wk_qs_with, fsr_tmp_with_e, fsi_tmp_with_e, wk_zfcos, wk_zfsin, &
    !$acc        ia_tmp_with, lmta2_tmp_with, il2_tmp_with, im2_tmp_with, lmt2_tmp_with, &
    !$acc        ia_tmp_without, lmta2_tmp_without, il2_tmp_without, im2_tmp_without, lmt2_tmp_without) &
    !$acc present(snl, ngabc, nbase, q, dion, dion_paw, ityp, ilmt, lmtt, lmta, ltp, mtp, ipaw)
    call pre_m_ES_Vnonlocal_W()
    !$acc kernels
    vnlph_l = 0.d0                      ! vnlph_l d(1:kg1,1:np_e,1:kimg)
    !$acc end kernels
    if ( noncol ) then 
       do ibl1=1,iba(ik),ibsize
          ibl2=min( ibl1+ibsize-1,iba(ik) )
          call calc_phase_blk(ik,ibl1,ibl2)
          if( icnt_with .gt. 0 ) then
             call Vnonlocal_W_part_with_blk_noncl( ibl1,ibl2,icnt_with   )
             call add_vnlph_l_with_eko_blk   (  ibsize,ibl1,ibl2,icnt_with, vnlph_l )

          endif
          if( icnt_without .gt. 0 ) then
             call Vnonlocal_W_part_no_blk_noncl( ibl1,ibl2,icnt_without )
             call add_vnlph_l_without_eko_blk(  ibsize,ibl1,ibl2,icnt_without,vnlph_l )

          endif
       enddo
    else 
       do ibl1=1,iba(ik),ibsize
          ibl2=min( ibl1+ibsize-1,iba(ik) )
          call calc_phase_blk(ik,ibl1,ibl2)
          if( icnt_with .gt. 0 ) then
             call Vnonlocal_W_part_with_blk   ( ibl1,ibl2,icnt_with   )
             call add_vnlph_l_with_eko_blk   (  ibsize,ibl1,ibl2,icnt_with, vnlph_l )

          endif
          if( icnt_without .gt. 0 ) then
             call Vnonlocal_W_part_without_blk( ibl1,ibl2,icnt_without)
             call add_vnlph_l_without_eko_blk(  ibsize,ibl1,ibl2,icnt_without,vnlph_l )

          endif
       enddo
    endif
! ========================================================================= 11.0

    !$acc end data

    deallocate( lmt2_tmp_without,ia_tmp_without,lmta2_tmp_without,il2_tmp_without,im2_tmp_without )
    deallocate( lmt2_tmp_with,   ia_tmp_with,   lmta2_tmp_with,   il2_tmp_with,   im2_tmp_with    )
    deallocate( wk_sc_without,   wk_ss_without, wk_sc_with,       wk_ss_with                      )
    deallocate( fsr_tmp_without, fsr_tmp_with )
#ifdef _CUDA_
    deallocate( fsr_tmp_with_e, fsi_tmp_with_e )
#endif
    deallocate(mdvdba)
    if(mdvdb_EXECUT) then
       deallocate( wk_qc_with,  wk_qs_with )
    endif
    if( k_symmetry(ik) /= GAMMA ) then
       deallocate( fsi_tmp_without, fsi_tmp_with )
    endif
    call dealloc_wkzfsincos()
    call tstatc0_end(id_sname)
    endif
  contains
    subroutine pre_m_ES_Vnonlocal_W_count( icnt_without, icnt_with )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer,intent(inout) :: icnt_without, icnt_with
      do ia = 1, natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(switch_of_eko_part == OFF) mdvdb= SKIP
         if( mdvdb==SKIP ) then
            do lmt2 = 1, ilmt(it)
               icnt_without = icnt_without + 1
            enddo
         elseif( mdvdb==EXECUT ) then
            do lmt2 = 1, ilmt(it)
               icnt_with = icnt_with + 1
            enddo
         endif
      enddo
    end subroutine pre_m_ES_Vnonlocal_W_count

    subroutine pre_m_ES_Vnonlocal_W()
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer :: icwithout, icwith
      !$acc data present(ityp, mdvdba, ilmt, lmta, ltp, mtp, &
      !$acc              ia_tmp_with, lmta2_tmp_with, il2_tmp_with, im2_tmp_with, lmt2_tmp_with, &
      !$acc              ia_tmp_without, lmta2_tmp_without, il2_tmp_without, im2_tmp_without, lmt2_tmp_without)
      !$acc serial
      icwithout = 0
      icwith = 0
      do ia = 1, natm
         it = ityp(ia)
         mdvdb = mdvdba(it)
         if(switch_of_eko_part == OFF) mdvdb= SKIP
         if( mdvdb==SKIP ) then
            do lmt2 = 1, ilmt(it)
               icwithout = icwithout + 1
               lmt2_tmp_without(icwithout)  = lmt2
               ia_tmp_without(icwithout)    =   ia
               lmta2_tmp_without(icwithout) = lmta(lmt2,ia)
               il2_tmp_without(icwithout)   =  ltp(lmt2,it)
               im2_tmp_without(icwithout)   =  mtp(lmt2,it)
            enddo
         elseif( mdvdb==EXECUT ) then
            do lmt2 = 1, ilmt(it)
               icwith = icwith + 1
               lmt2_tmp_with(icwith)  = lmt2
               ia_tmp_with(icwith)    =   ia
               lmta2_tmp_with(icwith) = lmta(lmt2,ia)
               il2_tmp_with(icwith)   =  ltp(lmt2,it)
               im2_tmp_with(icwith)   =  mtp(lmt2,it)
            enddo
         endif
      enddo
      !$acc end serial
      !$acc end data
    end subroutine pre_m_ES_Vnonlocal_W

    subroutine calc_phase_blk(ik,ibl1,ibl2)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in) :: ik, ibl1, ibl2
      integer :: i, nb, ia
      real(kind=DP) :: ph, f1, f2, f3
      integer :: id_sname0 = -1
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('calc_phase_blk ',id_sname0)
      !$acc data present(nbase, pos, ngabc, wk_zfcos, wk_zfsin)
      !$acc parallel
!cdir outerunroll=16
      !$acc loop gang
      do i = ibl1, ibl2
         nb = nbase(i,ik)
         !$acc loop vector
         do ia=1,natm
            f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
            ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
            wk_zfcos(i-ibl1+1,ia) = dcos(ph)
            wk_zfsin(i-ibl1+1,ia) = dsin(ph)
         enddo
      enddo
      !$acc end parallel
      !$acc end data
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname0)
    end subroutine calc_phase_blk

    subroutine Vnonlocal_W_part_without_blk(iblk1,iblk2,icnt_without)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by J. Koga, March 2010
!
      integer, intent(in) :: iblk1,iblk2, icnt_without
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic
      integer       :: ii,jj
      real(kind=DP) :: tmp
      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_without_blk ',id_sname)
      !$acc data present(snl, fsr_l, fsi_l, vlhxcQ, k_symmetry, iwei, &
      !$acc              ityp, dion, dion_paw, ilmt, lmtt, ltp, mtp, ipaw, &
      !$acc              wk_sc_without, wk_ss_without, fsr_tmp_without, fsi_tmp_without, wk_zfcos, wk_zfsin, &
      !$acc              ia_tmp_without, lmta2_tmp_without, il2_tmp_without, im2_tmp_without, lmt2_tmp_without)
      !$acc parallel
      !$acc loop gang
        do ii=1, icnt_without
      !$acc loop vector
        do jj=1, ibsize
        wk_sc_without(jj,ii) = 0.d0
        wk_ss_without(jj,ii) = 0.d0
        enddo
        enddo
       
      !$acc loop gang
        do ic = 1, icnt_without
           ia   = ia_tmp_without(ic)
           it   = ityp(ia)
           lmta2= lmta2_tmp_without(ic)
           il2  = il2_tmp_without(ic)
           im2  = im2_tmp_without(ic)
           lmt2 = lmt2_tmp_without(ic)
           fsr_tmp_without(1:np_e,ic) = fsr_l(1:np_e,lmta2,ik)
           if(k_symmetry(ik) /= GAMMA) then
              fsi_tmp_without(1:np_e,ic) = fsi_l(1:np_e,lmta2,ik)
           endif
           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)

              if(il1 == il2 .and. im1 == im2) then
!!$                 tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 endif
              else
!!$                 tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 end if
              endif
              tmp = tmp * iwei(ia)
              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
              if(mdl == 0 .or. mdl == 2) then
      !$acc loop vector
                 do i =  iblk1, iblk2
                    wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) + tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                 enddo
              else if(mdl == 1 .or. mdl == 3) then
      !$acc loop vector
                 do i =  iblk1, iblk2
                    wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) - tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                 enddo
              endif
           end do
        end do
        !$acc end parallel
        !$acc end data
        call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_without_blk

! ================================= added by K. Tagami ================ 11.0
    subroutine Vnonlocal_W_part_no_blk_noncl(iblk1,iblk2,icnt_without)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by J. Koga, March 2010
!
      integer, intent(in) :: iblk1,iblk2, icnt_without
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic
      integer :: mdvdb

      complex(kind=CMPLDP) :: tmp
      real(kind=DP) :: c1, c2

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_no_blk_noncl ',id_sname)

        wk_sc_without = 0.d0
        wk_ss_without = 0.d0

        do ic = 1, icnt_without
           ia   = ia_tmp_without(ic)
           it   = ityp(ia)

! ---------------------------------------------------- 11.0S
!#ifdef SKIP_TEST
           mdvdb = m_PP_include_vanderbilt_pot(it)
!#endif
! ---------------------------------------------------- 11.0S

           lmta2= lmta2_tmp_without(ic)
           il2  = il2_tmp_without(ic)
           im2  = im2_tmp_without(ic)
           lmt2 = lmt2_tmp_without(ic)
           fsr_tmp_without(1:np_e,ic) = fsr_l(1:np_e,lmta2,ik)
           if(k_symmetry(ik) /= GAMMA) then
              fsi_tmp_without(1:np_e,ic) = fsi_l(1:np_e,lmta2,ik)
           endif
           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
              if ( mdvdb == SKIP ) then
                 if ( il1 /= il2 ) cycle
                 if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
                    if ( im1 /= im2 ) cycle
                 endif
              endif
#endif
! ---------------------------------------------------- 11.0S

              tmp = dion_scr_noncl( lmt1, lmt2, ispin, ia )
              tmp = tmp * iwei(ia)

              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
              if(mdl == 0 .or. mdl == 2) then
                 do i =  iblk1, iblk2
                    c1 =  real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
	&                +aimag(tmp) *wk_zfsin(i-iblk1+1,ia) 
                    c2 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
	&                +aimag(tmp) *wk_zfcos(i-iblk1+1,ia) 

                    wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) &
	&                                        + c1 *snl(i,lmtt1,iksnl)
                    wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) &
	&                                        + c2 *snl(i,lmtt1,iksnl)
                 enddo
              else if(mdl == 1 .or. mdl == 3) then
                 do i =  iblk1, iblk2
                    c1 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
	&                +aimag(tmp) *wk_zfcos(i-iblk1+1,ia) 
                    c2 = -real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
	&                -aimag(tmp) *wk_zfsin(i-iblk1+1,ia) 

                    wk_sc_without(i-iblk1+1,ic) = wk_sc_without(i-iblk1+1,ic) &
	&                                       + c1 *snl(i,lmtt1,iksnl)
                    wk_ss_without(i-iblk1+1,ic) = wk_ss_without(i-iblk1+1,ic) &	
        &                                       + c2 *snl(i,lmtt1,iksnl)
                 enddo
              endif
           end do
        end do
        call tstatc0_end(id_sname)

    end subroutine Vnonlocal_W_part_no_blk_noncl
! ========================================================================== 11.0

    subroutine Vnonlocal_W_part_with_blk(iblk1,iblk2,icnt_with)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by J. Koga, March 2010
!
      integer, intent(in) :: iblk1,iblk2, icnt_with
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic
      real(kind=DP) :: tmp
      integer :: ii,jj
      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_with_blk ',id_sname)

      !$acc data present(snl, fsr_l, fsi_l, vlhxcQ, eko_l, q, k_symmetry, iwei, &
      !$acc              ityp, dion, dion_paw, ilmt, lmtt, ltp, mtp, ipaw, &
      !$acc              wk_sc_with, wk_ss_with, wk_qs_with, wk_qc_with, &
      !$acc              wk_zfcos, wk_zfsin, &
      !$acc              fsr_tmp_with, fsi_tmp_with, &
      !$acc              fsr_tmp_with_e, fsi_tmp_with_e, &
      !$acc              ia_tmp_with, lmta2_tmp_with, il2_tmp_with, im2_tmp_with, lmt2_tmp_with)
      !$acc parallel
      !$acc loop gang
      do ii=1, icnt_with
      !$acc loop vector
      do jj=1, ibsize
        wk_sc_with(jj,ii) = 0.d0
        wk_ss_with(jj,ii) = 0.d0
        wk_qc_with(jj,ii) = 0.d0
        wk_qs_with(jj,ii) = 0.d0
      enddo
      enddo

      !$acc loop gang
        do ic = 1, icnt_with
           ia   = ia_tmp_with(ic)
           it   = ityp(ia)
           lmta2= lmta2_tmp_with(ic)
           il2  = il2_tmp_with(ic)
           im2  = im2_tmp_with(ic)
           lmt2 = lmt2_tmp_with(ic)
           fsr_tmp_with(1:np_e,ic) = fsr_l(1:np_e,lmta2,ik)
           if(k_symmetry(ik) /= GAMMA) then
              fsi_tmp_with(1:np_e,ic) = fsi_l(1:np_e,lmta2,ik)
           endif
#ifdef _CUDA_
           fsr_tmp_with_e(1:np_e,ic) = fsr_tmp_with(1:np_e,ic)*eko_l(1:np_e,ik_for_pointing_eko)
           if(k_symmetry(ik) /= GAMMA) then
              fsi_tmp_with_e(1:np_e,ic) = fsi_tmp_with(1:np_e,ic)*eko_l(1:np_e,ik_for_pointing_eko)
           endif
#endif
           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)
              if(il1 == il2 .and. im1 == im2) then
!!$                 tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 endif
              else
!!$                 tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 if(ipaw(it)==0) then
                    tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                 else
                    tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                 end if
              endif
              tmp = tmp * iwei(ia)
              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
              if(mdl == 0 .or. mdl == 2) then
                 !$acc loop vector
                 do i =  iblk1, iblk2
                    wk_sc_with(i-iblk1+1,ic) = wk_sc_with(i-iblk1+1,ic) + tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,ic) = wk_ss_with(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                 enddo
              else if(mdl == 1 .or. mdl == 3) then
                 !$acc loop vector
                 do i =  iblk1, iblk2
                    wk_sc_with(i-iblk1+1,ic) = wk_sc_with(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,ic) = wk_ss_with(i-iblk1+1,ic) - tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                 enddo
              endif
              if( il1 == il2 .and. im1 == im2) then
                 tmp = q(lmt1,lmt2,it)*iwei(ia)
                 if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                 if(mdl == 0 .or. mdl == 2) then
                 !$acc loop vector
                    do i =  iblk1, iblk2
                       wk_qc_with(i-iblk1+1,ic) = wk_qc_with(i-iblk1+1,ic) + tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                       wk_qs_with(i-iblk1+1,ic) = wk_qs_with(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    enddo
                 else if(mdl == 1 .or. mdl == 3) then
                 !$acc loop vector
                    do i =  iblk1, iblk2
                       wk_qc_with(i-iblk1+1,ic) = wk_qc_with(i-iblk1+1,ic) - tmp*wk_zfsin(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                       wk_qs_with(i-iblk1+1,ic) = wk_qs_with(i-iblk1+1,ic) - tmp*wk_zfcos(i-iblk1+1,ia)*snl(i,lmtt1,iksnl)
                    enddo
                 end if
              end if
           end do
        end do
        !$acc end parallel
        !$acc end data
        call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_with_blk

! ================================== added by K. Tagami =================== 11.0
    subroutine Vnonlocal_W_part_with_blk_noncl(iblk1,iblk2,icnt_with)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by J. Koga, March 2010
!
      integer, intent(in) :: iblk1,iblk2, icnt_with
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      integer       :: ia, i, lmta2, il2, im2, lmt2, ic
      integer :: mdvdb

      complex(kind=CMPLDP) :: tmp
      real(kind=DP) :: c1, c2

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_with_blk_noncl ',id_sname)


        wk_sc_with = 0.d0
        wk_ss_with = 0.d0
        wk_qc_with = 0.d0
        wk_qs_with = 0.d0

        do ic = 1, icnt_with
           ia   = ia_tmp_with(ic)
           it   = ityp(ia)

! ---------------------------------------------------- 11.0S
!#ifdef SKIP_TEST
           mdvdb = m_PP_include_vanderbilt_pot(it)
!#endif
! ---------------------------------------------------- 11.0S

           lmta2= lmta2_tmp_with(ic)
           il2  = il2_tmp_with(ic)
           im2  = im2_tmp_with(ic)
           lmt2 = lmt2_tmp_with(ic)
           fsr_tmp_with(1:np_e,ic) = fsr_l(1:np_e,lmta2,ik)
           if(k_symmetry(ik) /= GAMMA) then
              fsi_tmp_with(1:np_e,ic) = fsi_l(1:np_e,lmta2,ik)
           endif
           do lmt1 = 1,ilmt(it)
              lmtt1 = lmtt(lmt1,it)
              il1   = ltp(lmt1,it)
              im1   = mtp(lmt1,it)
              il11  = il1 - 1
              mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
              if ( mdvdb == SKIP ) then
                 if ( il1 /= il2 ) cycle
                 if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
                    if ( im1 /= im2 ) cycle
                 endif
              endif
#endif
! ---------------------------------------------------- 11.0S

              tmp = dion_scr_noncl( lmt1, lmt2, ispin, ia )
              tmp = tmp * iwei(ia)

              if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
              if(mdl == 0 .or. mdl == 2) then
                 do i =  iblk1, iblk2
                    c1 =  real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
        &                +aimag(tmp) *wk_zfsin(i-iblk1+1,ia)
                    c2 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
        &                +aimag(tmp) *wk_zfcos(i-iblk1+1,ia)

                    wk_sc_with(i-iblk1+1,ic) = wk_sc_with(i-iblk1+1,ic) &
	&                                    + c1 *snl(i,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,ic) = wk_ss_with(i-iblk1+1,ic) &
	&                                    + c2 *snl(i,lmtt1,iksnl)
                 enddo
              else if(mdl == 1 .or. mdl == 3) then
                 do i =  iblk1, iblk2

                    c1 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
        &                +aimag(tmp) *wk_zfcos(i-iblk1+1,ia)
                    c2 = -real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
        &                -aimag(tmp) *wk_zfsin(i-iblk1+1,ia)

                    wk_sc_with(i-iblk1+1,ic) = wk_sc_with(i-iblk1+1,ic) &
	&                                     + c1 *snl(i,lmtt1,iksnl)
                    wk_ss_with(i-iblk1+1,ic) = wk_ss_with(i-iblk1+1,ic) &
	&                                     + c2 *snl(i,lmtt1,iksnl)
                 enddo
              endif

! ----------------------------------------------------- 11.0S
!              if( il1 == il2 .and. im1 == im2) then

              if ( mdvdb == EXECUT .and. il1 == il2  ) then
                 if ( SpinOrbit_mode /= BuiltIn ) then
                    if ( im1 /= im2 ) cycle
                 endif
! ------------------------------------------------------ 11.0S

                 tmp = q_noncl(lmt1,lmt2,ispin,it)*iwei(ia)

                 if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
                 if(mdl == 0 .or. mdl == 2) then
                    do i =  iblk1, iblk2

                      c1 =  real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
        &                  +aimag(tmp) *wk_zfsin(i-iblk1+1,ia)
                      c2 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
        &                  +aimag(tmp) *wk_zfcos(i-iblk1+1,ia)

                      wk_qc_with(i-iblk1+1,ic) = wk_qc_with(i-iblk1+1,ic) &
	&                                       + c1 *snl(i,lmtt1,iksnl)
                      wk_qs_with(i-iblk1+1,ic) = wk_qs_with(i-iblk1+1,ic) &
	&                                       + c2 *snl(i,lmtt1,iksnl)
                    enddo
                 else if(mdl == 1 .or. mdl == 3) then
                    do i =  iblk1, iblk2

                      c1 = -real(tmp)  *wk_zfsin(i-iblk1+1,ia) &
        &                  +aimag(tmp) *wk_zfcos(i-iblk1+1,ia)
                      c2 = -real(tmp)  *wk_zfcos(i-iblk1+1,ia) &
        &                  -aimag(tmp) *wk_zfsin(i-iblk1+1,ia)

                      wk_qc_with(i-iblk1+1,ic) = wk_qc_with(i-iblk1+1,ic) &
	&                                      + c1 *snl(i,lmtt1,iksnl)
                      wk_qs_with(i-iblk1+1,ic) = wk_qs_with(i-iblk1+1,ic) &
	&                                      + c2 *snl(i,lmtt1,iksnl)
                    enddo
                 end if
              end if
           end do
        end do
        call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_with_blk_noncl
! ============================================================================== 11.0

    subroutine add_vnlph_l_with_eko_blk(ibsize,ibl1,ibl2,icnt_with,vnlph)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in) :: ibsize, ibl1, ibl2, icnt_with
      real(kind=DP), intent(inout), dimension(kg1,np_e,kimg) :: vnlph
      integer          :: ic, ib
      integer       :: icsize
      real(kind=DP) :: alpha, beta
      integer :: id_sname = -1
      !$acc host_data use_device(wk_sc_with, fsr_tmp_with, wk_ss_with, fsi_tmp_with, &
      !$acc wk_qc_with, wk_qs_with, fsr_tmp_with_e, fsi_tmp_with_e, &
      !$acc vnlph)
      call tstatc0_begin('add_vnlph_l_with_eko_blk ',id_sname)
      if(kimg == 1) then
         icsize=ibl2-ibl1+1
#ifdef _CUDA_
         alpha= 1.d0;  beta= 1.d0
         call  cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
        &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
         alpha=-1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
        &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
#else
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_with, &
        &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
         alpha=-1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_with, &
        &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
#endif

! ========================== modified by K. Tagami =============== 11.0
!         do ic = 1, icnt_with
!            do ib = 1, np_e                                              ! MPI
!               fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik)
!               fsi_tmp_with(ib,ic) = fsi_tmp_with(ib,ic)*eko_l(ib,ik)
!            enddo
!         enddo
#ifndef _CUDA_
         do ic = 1, icnt_with
            do ib = 1, np_e                                              ! MPI
               fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik_for_pointing_eko)
               fsi_tmp_with(ib,ic) = fsi_tmp_with(ib,ic)*eko_l(ib,ik_for_pointing_eko)
            enddo
         enddo
#endif
! ================================================================ 11.0

#ifdef _CUDA_
         alpha=-1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
        &              alpha,wk_qc_with,ibsize, fsr_tmp_with_e,np_e, beta,vnlph(ibl1,1,1),kg1 )
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
        &              alpha,wk_qs_with,ibsize, fsi_tmp_with_e,np_e, beta,vnlph(ibl1,1,1),kg1 )
#else
         alpha=-1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_with, &
        &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_with, &
        &              alpha,wk_qs_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1 )
#endif
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
            icsize=ibl2-ibl1+1
#ifdef _CUDA_
            alpha= 1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#else
            alpha= 1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif

! ======================================== modified by K. Tagami ============== 11.0
!            do ic = 1, icnt_with
!               do ib = 1, np_e                                              ! MPI
!                  fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik)
!               enddo
!            enddo
#ifndef _CUDA_
            do ic = 1, icnt_with
               do ib = 1, np_e                                              ! MPI
                  fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik_for_pointing_eko)
               enddo
            enddo
#endif
! =============================================================================== 11.0

#ifdef _CUDA_
            alpha=-1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsr_tmp_with_e,np_e, beta,vnlph(ibl1,1,1),kg1)
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsr_tmp_with_e,np_e, beta,vnlph(ibl1,1,2),kg1)
#else
            alpha=-1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif
         else
            icsize=ibl2-ibl1+1
#ifdef _CUDA_
            alpha= 1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha=-1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha= 1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
            alpha= 1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#else
            alpha= 1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha=-1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha= 1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_sc_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
            alpha= 1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_ss_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif
! ==================================== modified by K. Tagami ================ 11.0
!            do ic = 1, icnt_with
!               do ib = 1, np_e                                              ! MPI
!                  fsi_tmp_with(ib,ic) = fsi_tmp_with(ib,ic)*eko_l(ib,ik)
!                  fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik)
!               enddo
!            enddo
#ifndef _CUDA_
            do ic = 1, icnt_with
               do ib = 1, np_e                                              ! MPI
                  fsi_tmp_with(ib,ic) = fsi_tmp_with(ib,ic)*eko_l(ib,ik_for_pointing_eko)
                  fsr_tmp_with(ib,ic) = fsr_tmp_with(ib,ic)*eko_l(ib,ik_for_pointing_eko)
               enddo
            enddo
#endif
! ============================================================================ 11.0
#ifdef _CUDA_
            alpha=-1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsr_tmp_with_e,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha= 1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsi_tmp_with_e,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha=-1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsi_tmp_with_e,np_e, beta,vnlph(ibl1,1,2),kg1)
            alpha=-1.d0;  beta= 1.d0
            call cublasDgemm('N', 'T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsr_tmp_with_e,np_e, beta,vnlph(ibl1,1,2),kg1)
#else
            alpha=-1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha= 1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,1),kg1)
            alpha=-1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qc_with,ibsize, fsi_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
            alpha=-1.d0;  beta= 1.d0
            call DGEMM__('N','T', icsize,np_e,icnt_with, &
           &              alpha,wk_qs_with,ibsize, fsr_tmp_with,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif
         end if
      end if
      !$acc end host_data
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_with_eko_blk

    subroutine add_vnlph_l_without_eko_blk(ibsize,ibl1,ibl2,icnt_without,vnlph)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in) :: ibsize, ibl1, ibl2, icnt_without
      real(kind=DP), intent(inout), dimension(kg1,np_e,kimg) :: vnlph
      integer       :: icsize
      real(kind=DP) :: alpha, beta
      integer :: id_sname = -1
      !$acc host_data use_device(wk_sc_without, fsr_tmp_without, wk_ss_without, fsi_tmp_without, vnlph_l)
      call tstatc0_begin('add_vnlph_l_without_eko_blk ',id_sname)

      if(kimg == 1) then
         icsize=ibl2-ibl1+1
#ifdef _CUDA_
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph_l(ibl1,1,1),kg1 )
         alpha=-1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph_l(ibl1,1,1),kg1 )
#else
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ibl1,1,1),kg1 )
         alpha=-1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ibl1,1,1),kg1 )
#endif
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
         icsize=ibl2-ibl1+1
#ifdef _CUDA_
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph_l(ibl1,1,1),kg1)
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph_l(ibl1,1,2),kg1)
#else
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ibl1,1,1),kg1)
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif
         else
         icsize=ibl2-ibl1+1
#ifdef _CUDA_
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph_l(ibl1,1,1),kg1)
         alpha=-1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph_l(ibl1,1,1),kg1)
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsi_tmp_without,np_e, beta,vnlph_l(ibl1,1,2),kg1)
         alpha= 1.d0;  beta= 1.d0
         call cublasDgemm('N', 'T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph_l(ibl1,1,2),kg1)
#else
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ibl1,1,1),kg1)
         alpha=-1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ibl1,1,1),kg1)
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_sc_without,ibsize, fsi_tmp_without,np_e, beta,vnlph(ibl1,1,2),kg1)
         alpha= 1.d0;  beta= 1.d0
         call DGEMM__('N','T', icsize,np_e,icnt_without, &
        &              alpha,wk_ss_without,ibsize, fsr_tmp_without,np_e, beta,vnlph(ibl1,1,2),kg1)
#endif
         end if
      end if
      call tstatc0_end(id_sname)
      !$acc end host_data
    end subroutine add_vnlph_l_without_eko_blk
  end subroutine m_ES_Vnonlocal_W

  subroutine m_ES_AtaulmnaG(hardpart)
    logical, intent(in), optional :: hardpart
    complex(kind=DP), allocatable, dimension(:,:) :: wkexp
    integer :: ia,it,ig,iksnl,ik,lmt1,lmt2,im1,im2,lmtt1,lmtt2,lmta1,lmta2,il1,il2,mil
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
    allocate(wkexp(maxval(iba),natm))
    hp = .false.
    if(present(hardpart)) hp = hardpart
    AtaulmaG = 0.d0
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
      do ibl1=1,iba(ik),ibsize
      icount = 0
      ibl2=min( ibl1+ibsize-1,iba(ik) )
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
!            do ig=1,iba(ik)
            do ig=ibl1,ibl2
              ctmp  = cfac*wkexp(ig,ia)*snl(ig,lmtt1,iksnl)
              AtaulmaG(ig,lmta2,ik,1) = AtaulmaG(ig,lmta2,ik,1) + tmp *  real(ctmp)
              AtaulmaG(ig,lmta2,ik,2) = AtaulmaG(ig,lmta2,ik,2) + tmp * dimag(ctmp)
            enddo
            if(hp .and. mdvdb==ON .and. il1 == il2 .and. im1 == im2) then
              tmp = q(lmt1,lmt2,it)*iwei(ia)
!              do ig=1,iba(ik)
              do ig=ibl1,ibl2
                ctmp  = cfac*wkexp(ig,ia)*snl(ig,lmtt1,iksnl)
                BtaulmaG(ig,icount,ik,1) = BtaulmaG(ig,icount,ik,1) + tmp *  real(ctmp)
                BtaulmaG(ig,icount,ik,2) = BtaulmaG(ig,icount,ik,2) + tmp * dimag(ctmp)
              enddo
            endif
          enddo
        enddo
      enddo
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
        do ig=1,iba(ik)
          nb = nbase(ig,ik)
          ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
          wkexp(ig,ia) = dcmplx(dcos(ph),-dsin(ph))
        enddo
      enddo
    end subroutine cal_phase

  end subroutine m_ES_AtaulmnaG

  subroutine m_ES_Vnonlocal_W_precphase(ik,iksnl,ispin,switch_of_eko_part,map)
    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part
    logical, intent(in), dimension(np_e), optional :: map
    integer :: i, ik_for_pointing_eko, ib, icount, ia, it, mdvdb, lmt2, lmta2
    integer :: ibsize,ibl,ibl1,ibl2,icountb,ndata
    logical :: umap
    real(kind=DP) :: e
    real(kind=DP), allocatable, dimension(:,:) :: efsr_l, efsi_l
    real(kind=DP), allocatable, dimension(:,:) :: atmpr,atmpi
    real(kind=DP), allocatable, dimension(:,:) :: fsrtmp,fsitmp,vnltmpr,vnltmpi
    integer :: id_sname = -1
    if(nblocksize_vnonlocal_is_given) then
       ibsize = nblocksize_vnonlocal_w
    else
       ibsize= nb_vnonlocal_default
    end if
    call tstatc0_begin('m_ES_Vnonlocal_W_precphase ',id_sname,1)
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
    if(umap)then
      allocate(fsrtmp(ndata,nlmta))
      allocate(vnltmpr(ibsize,ndata))
      if(kimg==2) allocate(vnltmpi(ibsize,ndata))
      if(k_symmetry(ik) /= GAMMA) allocate(fsitmp(ndata,nlmta))
      icount = 0
      do i=1,np_e
        if(.not.map(i)) then
          icount = icount+1
          fsrtmp(icount,:) = fsr_l(i,:,ik)
          if(k_symmetry(ik) /= GAMMA) fsitmp(icount,:) = fsi_l(i,:,ik)
        endif
      enddo
    endif
    allocate(atmpr(ibsize,nlmta))
    allocate(atmpi(ibsize,nlmta))
    if(.not.umap) then
      do ibl1=1,iba(ik),ibsize
        ibl2=min( ibl1+ibsize-1,iba(ik) )
        atmpr(1:ibl2-ibl1+1,1:nlmta) = AtaulmaG(ibl1:ibl2,1:nlmta,ik,1)
        atmpi(1:ibl2-ibl1+1,1:nlmta) = AtaulmaG(ibl1:ibl2,1:nlmta,ik,2)
        if(kimg==1) then
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta, 1.d0,atmpr,ibsize,fsr_l(:,:,ik),np_e,0.d0,vnlph_l(ibl1,1,1),kg1)
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta,-1.d0,atmpi,ibsize,fsi_l(:,:,ik),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
        else
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta, 1.d0,atmpr,ibsize,fsr_l(:,:,ik),np_e,0.d0,vnlph_l(ibl1,1,1),kg1)
          if(k_symmetry(ik)/=GAMMA) &
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta,-1.d0,atmpi,ibsize,fsi_l(:,:,ik),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta, 1.d0,atmpi,ibsize,fsr_l(:,:,ik),np_e,0.d0,vnlph_l(ibl1,1,2),kg1)
          if(k_symmetry(ik)/=GAMMA) &
          call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta, 1.d0,atmpr,ibsize,fsi_l(:,:,ik),np_e,1.d0,vnlph_l(ibl1,1,2),kg1)
        endif
      enddo
    else
      do ibl1=1,iba(ik),ibsize
        ibl2=min( ibl1+ibsize-1,iba(ik) )
        atmpr(1:ibl2-ibl1+1,1:nlmta) = AtaulmaG(ibl1:ibl2,1:nlmta,ik,1)
        atmpi(1:ibl2-ibl1+1,1:nlmta) = AtaulmaG(ibl1:ibl2,1:nlmta,ik,2)
        if(kimg==1) then
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta, 1.d0,atmpr,ibsize,fsrtmp(:,:),ndata,0.d0,vnltmpr,ibsize)
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta,-1.d0,atmpi,ibsize,fsitmp(:,:),ndata,1.d0,vnltmpr,ibsize)
        else
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta, 1.d0,atmpr,ibsize,fsrtmp(:,:),ndata,0.d0,vnltmpr,ibsize)
          if(k_symmetry(ik)/=GAMMA) &
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta,-1.d0,atmpi,ibsize,fsitmp(:,:),ndata,1.d0,vnltmpr,ibsize)
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta, 1.d0,atmpi,ibsize,fsrtmp(:,:),ndata,0.d0,vnltmpi,ibsize)
          if(k_symmetry(ik)/=GAMMA) &
          call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta, 1.d0,atmpr,ibsize,fsitmp(:,:),ndata,1.d0,vnltmpi,ibsize)
        endif
      enddo
      icount = 0
      do i=1,np_e
        if(.not.map(i)) then
          icount = icount+1
          vnlph_l(ibl1:ibl2,i,1) = vnltmpr(1:ibl2-ibl1+1,icount)
          if(kimg==2) &
          vnlph_l(ibl1:ibl2,i,2) = vnltmpi(1:ibl2-ibl1+1,icount)
        endif
      enddo
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
        if(k_symmetry(ik) /= GAMMA) allocate(efsi_l(np_e,nlmta_us))
        do ia=1,natm
          it = ityp(ia)
          mdvdb = m_PP_include_vanderbilt_pot(it)
          do lmt2 = 1, ilmt(it)
            if(mdvdb==ON) then
              icount = icount+1
              lmta2 = lmta(lmt2,ia)
              do ib=1,np_e
                efsr_l(ib,icount) = -eko_l(ib,ik_for_pointing_eko) * fsr_l(ib,lmta2,ik)
                if(k_symmetry(ik) /=  GAMMA) &
                efsi_l(ib,icount) = -eko_l(ib,ik_for_pointing_eko) * fsi_l(ib,lmta2,ik)
              enddo
            endif
          enddo
        enddo
   
        do ibl1=1,iba(ik),ibsize
          ibl2=min( ibl1+ibsize-1,iba(ik) )
          atmpr(1:ibl2-ibl1+1,1:nlmta_us) = BtaulmaG(ibl1:ibl2,1:nlmta_us,ik,1)
          atmpi(1:ibl2-ibl1+1,1:nlmta_us) = BtaulmaG(ibl1:ibl2,1:nlmta_us,ik,2)
          if(kimg==1) then
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us, 1.d0,atmpr,ibsize,efsr_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us,-1.d0,atmpi,ibsize,efsi_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
          else
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us, 1.d0,atmpr,ibsize,efsr_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
            if(k_symmetry(ik) /= GAMMA) &
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us,-1.d0,atmpi,ibsize,efsi_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,1),kg1)
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us, 1.d0,atmpi,ibsize,efsr_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,2),kg1)
            if(k_symmetry(ik) /= GAMMA) &
            call DGEMM__('N','T',ibl2-ibl1+1,np_e,nlmta_us, 1.d0,atmpr,ibsize,efsi_l(:,:),np_e,1.d0,vnlph_l(ibl1,1,2),kg1)
          endif
        enddo
        deallocate(efsr_l)
        if(k_symmetry(ik) /= GAMMA) deallocate(efsi_l)
      else
        allocate(efsr_l(ndata,nlmta_us))
        if(k_symmetry(ik) /= GAMMA) allocate(efsi_l(ndata,nlmta_us))
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
                  efsr_l(icountb,icount) = -eko_l(ib,ik_for_pointing_eko) * fsr_l(ib,lmta2,ik)
                  if(k_symmetry(ik) /=  GAMMA) &
                  efsi_l(icountb,icount) = -eko_l(ib,ik_for_pointing_eko) * fsi_l(ib,lmta2,ik)
                endif
              enddo
            endif
          enddo
        enddo
   
        do ibl1=1,iba(ik),ibsize
          ibl2=min( ibl1+ibsize-1,iba(ik) )
          atmpr(1:ibl2-ibl1+1,1:nlmta_us) = BtaulmaG(ibl1:ibl2,1:nlmta_us,ik,1)
          atmpi(1:ibl2-ibl1+1,1:nlmta_us) = BtaulmaG(ibl1:ibl2,1:nlmta_us,ik,2)
          if(kimg==1) then
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us, 1.d0,atmpr,ibsize,efsr_l(:,:),ndata,1.d0,vnltmpr,ibsize)
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us,-1.d0,atmpi,ibsize,efsi_l(:,:),ndata,1.d0,vnltmpr,ibsize)
          else
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us, 1.d0,atmpr,ibsize,efsr_l(:,:),ndata,1.d0,vnltmpr,ibsize)
            if(k_symmetry(ik) /= GAMMA) &
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us,-1.d0,atmpi,ibsize,efsi_l(:,:),ndata,1.d0,vnltmpr,ibsize)
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us, 1.d0,atmpi,ibsize,efsr_l(:,:),ndata,1.d0,vnltmpi,ibsize)
            if(k_symmetry(ik) /= GAMMA) &
            call DGEMM__('N','T',ibl2-ibl1+1,ndata,nlmta_us, 1.d0,atmpr,ibsize,efsi_l(:,:),ndata,1.d0,vnltmpi,ibsize)
          endif
          icount = 0
          do i=1,np_e
            if(.not.map(i)) then
              icount = icount+1
              vnlph_l(ibl1:ibl2,i,1) = vnlph_l(ibl1:ibl2,i,1) + vnltmpr(1:ibl2-ibl1+1,icount)
              if(kimg==2) &
              vnlph_l(ibl1:ibl2,i,2) = vnlph_l(ibl1:ibl2,i,2) + vnltmpi(1:ibl2-ibl1+1,icount)
            endif
          enddo
        enddo
        deallocate(efsr_l)
        if(k_symmetry(ik) /= GAMMA) deallocate(efsi_l)
      endif
    endif
    deallocate(atmpr)
    deallocate(atmpi)
    if(umap) then
      deallocate(fsrtmp)
      deallocate(vnltmpr)
      if(kimg==2) deallocate(vnltmpi)
      if(k_symmetry(ik) /= GAMMA) deallocate(fsitmp)
    endif
    call tstatc0_end(id_sname)
  end subroutine m_ES_Vnonlocal_W_precphase

#else
  subroutine m_ES_Vnonlocal_W(ik,iksnl,ispin,switch_of_eko_part)
    integer, intent(in) :: ik,iksnl,ispin
    integer, intent(in) :: switch_of_eko_part

    real(kind=DP), allocatable,dimension(:):: zfcos_mpi0,zfsin_mpi0 !d(np_kg1_k)
    integer :: mdvdb, it, ia, lmt2, lmta2, il2, im2, isize
    integer :: id_sname = -1
    real(kind=DP),allocatable,dimension(:) :: sc, ss, qc, qs

! ================================== added by K. Tagami =============== 11.0
    integer :: ik_for_pointing_eko
! ===================================================================== 11.0

#ifndef SX
    complex(kind=CMPLDP),allocatable,dimension(:,:) :: vnlph_cmplx_l
    real(kind=DP),allocatable,dimension(:,:,:) :: vnlph_red_l
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
!!$    ncache = (cachesize(3)*1024)*3/4
    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) then
       ibsize = iba(ik)
    else
    if(switch_of_eko_part == OFF) then ! without_eko_part
      if(kimg == 1) then ! vnlph_l(i,ib,1):np_e,sc(i):1,ss(i):1
        ibsize=ncache/(8*(np_e+2))
      else if(kimg == 2) then
        if(k_symmetry(ik) == GAMMA) then ! vnlph_l(i,ib,1):np_e
          ibsize=ncache/(8*(np_e*2+2))   ! vnlph_l(i,ib,2):np_e,sc(i):1,ss(i):1
        else ! vnlph_l(i,ib,1):np_e,vnlph_l(i,ib,2):np_e,sc(i):1,ss(i):1
          ibsize=ncache/(8*(np_e*2+2))
        endif
      endif
    else ! with_eko_part
      if(kimg == 1) then
        if(k_symmetry(ik) == GAMMA) then ! vnlph_l(i,ib,1):np_e
          ibsize=ncache/(8*(np_e+4))     ! sc(i):1,qc(i):1,ss(i):1,qs(i):1
        else ! vnlph_l(i,ib,1):np_e,vnlph_l(i,ib,2):np_e
          ibsize=ncache/(8*(np_e*2+4))   ! sc(i):1,qc(i):1,ss(i):1,qs(i):1
        endif
      else ! vnlph_l(i,ib,1):np_e,vnlph_l(i,ib,2):np_e
          ibsize=ncache/(8*(np_e*2+4))   ! sc(i):1,qc(i):1,ss(i):1,qs(i):1
      endif
    endif
    end if
! nec debug
    if(ipribetar >= 2) then
       write(nfout,885) 'ik=',ik,' np_e=',np_e,' iba(ik)=',iba(ik),' ibsize=',ibsize,&
            & ' ncache=',ncache
    end if
885 format(a,i3,a,i4,a,i6,a,i6,a,i8)
! NEC tune <-------------------------------
#endif
    call tstatc0_begin('m_ES_Vnonlocal_W ',id_sname,level=1)

#ifdef SX
    call m_ES_alloc_zfsincos(ik)
    if(flag_mpi_g_dot_r_k) then
       allocate(zfcos_mpi0(mp_kg1_k))
       allocate(zfsin_mpi0(mp_kg1_k))
    end if
    isize = kg1
#else
    call alloc_zfsincos(ibsize)
    isize = ibsize
    allocate(vnlph_red_l(ibsize,np_e,kimg))
#endif

    allocate(sc(isize),ss(isize))
    do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       if(mdvdb == EXECUT) goto 1001
    end do
    goto 1002
1001 allocate(qc(isize),qs(isize))
1002 continue

    vnlph_l = 0.d0                      ! vnlph_l d(1:kg1,1:np_e,1:kimg)

! ======================================== added by K. Tagami ========= 11.0
    if ( noncol ) then
       ik_for_pointing_eko = ( iksnl -1 )*ndim_spinor +1
    else
       ik_for_pointing_eko = ik
    endif
! ====================================================================== 11.0

    Loop_ntyp: do it = 1, ntyp
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
#ifdef SX
       Loop_natm : do ia = 1, natm
          if(ityp(ia) /= it) cycle
          call calc_phase_mpi(ik,ia) ! pos(ia,1:3),ngabc -> zfcos,zfsin
          do lmt2 = 1, ilmt(it)
             lmta2 = lmta(lmt2,ia)
             il2   = ltp(lmt2,it)
             im2   = mtp(lmt2,it)

! ==================================== modified by K. Tagami ============ 11.0
!             call Vnonlocal_W_part_sum_over_lmt1(iba(ik))
!
             if ( noncol ) then
               call Vnonlocal_W_part_sum_lmt1_noncl(iba(ik))
             else
               call Vnonlocal_W_part_sum_over_lmt1(iba(ik))
             endif
! ====================================================================== 11.0

             if(mdvdb == SKIP) then
                call add_vnlph_l_without_eko_part(kg1,iba(ik),vnlph_l)
             else if(mdvdb == EXECUT) then
                call add_vnlph_l_with_eko_part(kg1,iba(ik),vnlph_l)
             endif
          end do
       end do Loop_natm
#else
! Revised by T. Yamasaki, September 2008 
       Gdiv_Loop: do ibl1=1,iba(ik),ibsize
          ibl2=ibl1+ibsize-1
          if(ibl2.gt.iba(ik)) ibl2=iba(ik)
          if(ipribetar >= 2) write(nfout,'(" Gdiv_loop = ",i8," ibl1, ibl2 = ",2i8)') ceiling(dble(ibl1)/ibsize),ibl1,ibl2
!!$          if(kimg==2) then
!!$             allocate(vnlph_cmplx_l(ibl2-ibl1+1,np_e))
!!$             vnlph_cmplx_l = cmplx(0.d0,0.d0)
!!$          end if
          vnlph_red_l = 0.d0
! NEC tune <-------------------------------
          Loop_natm : do ia = 1, natm
             if(ityp(ia) /= it) cycle
             call calc_phase_div(ik,ia) ! pos(ia,1:3),ngabc -> zfcos,zfsin
             do lmt2 = 1, ilmt(it)
                lmta2 = lmta(lmt2,ia)
                il2   = ltp(lmt2,it)
                im2   = mtp(lmt2,it)

! ======================================== modified by K. Tagami ============ 11.0
!                call Vnonlocal_W_part_sum_over_lmt1(ibl1,ibl2)
!
                if ( noncol ) then
                  call Vnonlocal_W_part_sum_lmt1_noncl(ibl1,ibl2)
                else
                  call Vnonlocal_W_part_sum_over_lmt1(ibl1,ibl2)
                endif
! ============================================================================ 11.0

                if(mdvdb == SKIP) then
                   call add_vnlph_l_without_eko_part(ibsize,ibl2-ibl1+1,vnlph_red_l)
                else if(mdvdb == EXECUT) then
                   call add_vnlph_l_with_eko_part(ibsize,ibl2-ibl1+1,vnlph_red_l)
                endif
             end do
          end do Loop_natm
!!$          if(kimg==2) call cp_vnlphcmplx2vnlph()
!!$          if(kimg==2) deallocate(vnlph_cmplx_l)
          call cp_vnlphred2vnlph()
! NEC tune ------------------------------->
       end do Gdiv_Loop
! NEC tune <-------------------------------
#endif
    end do Loop_ntyp
    if(allocated(qc)) deallocate(qs,qc)
    deallocate(ss,sc)
#ifdef SX
    if(flag_mpi_g_dot_r_k) deallocate(zfsin_mpi0,zfcos_mpi0)
    call m_ES_dealloc_zfsincos()
#else
    deallocate(vnlph_red_l)
    call dealloc_zfsincos()
#endif
    call tstatc0_end(id_sname)
  contains
!!$    subroutine cp_vnlphcmplx2vnlph()
!!$      integer :: ib, i
!!$      integer :: id_sname0 = -1
!!$      call tstatc0_begin('cp_vnlphcmplx2vnlph ',id_sname0)
!!$
!!$      do ib = 1, np_e
!!$         do i = ibl1, ibl2
!!$            vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + dreal(vnlph_cmplx_l(i-ibl1+1,ib))
!!$            vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + dimag(vnlph_cmplx_l(i-ibl1+1,ib))
!!$         end do
!!$      end do
!!$      call tstatc0_end(id_sname0)
!!$    end subroutine cp_vnlphcmplx2vnlph

#ifndef SX
    subroutine cp_vnlphred2vnlph()
      integer :: ib, i
      integer :: id_sname0 = -1
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('cp_vnlphred2vnlph ',id_sname0)

      if(kimg==1) then
         do ib = 1, np_e
            do i = ibl1, ibl2
               vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + vnlph_red_l(i-ibl1+1,ib,1)
            end do
         end do
      else if(kimg==2) then
         do ib = 1, np_e
            do i = ibl1, ibl2
               vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + vnlph_red_l(i-ibl1+1,ib,1)
               vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + vnlph_red_l(i-ibl1+1,ib,2)
            end do
         end do
      end if
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname0)
    end subroutine cp_vnlphred2vnlph
#endif

    subroutine calc_phase_mpi(ik,ia)
      integer, intent(in) :: ik, ia
      integer :: i, nb, ip
      real(kind=DP) :: ph, f1, f2, f3
      integer :: id_sname0 = -1
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('calc_phase_mpi ',id_sname0)

      f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
      if(flag_mpi_g_dot_r_k) then
         zfcos_mpi0(np_kg1_k+1:mp_kg1_k) = 0.d0; zfsin_mpi0(np_kg1_k+1:mp_kg1_k) = 0.d0
!!$         imax = min(iba(ik),iend_kg1_k)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
         do i = 1, np_kg1_k
            ip = i + ista_kg1_k - 1
!!$         do i = ista_kg1_k, imax
            nb = nbase(ip,ik)
            ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
            zfcos_mpi0(i) = dcos(ph)
            zfsin_mpi0(i) = dsin(ph)
         end do
         call mpi_allgather(zfcos_mpi0, mp_kg1_k, mpi_double_precision &
              &  , zfcos,mp_kg1_k,mpi_double_precision, mpi_nbmx_world_k(ng_nbmx_k),ierr)
         call mpi_allgather(zfsin_mpi0, mp_kg1_k, mpi_double_precision &
              &  , zfsin,mp_kg1_k,mpi_double_precision, mpi_nbmx_world_k(ng_nbmx_k),ierr)
!!$         call mpi_allgatherv(zfcos_mpi0,nel_kg1_k(mype),mpi_double_precision &
!!$              & ,zfcos, nel_kg1_k, idisp_kg1_k, mpi_double_precision, mpi_nbmx_world_k(ng_nbmx_k),ierr)
!!$         call mpi_allgatherv(zfsin_mpi0,nel_kg1_k(mype),mpi_double_precision &
!!$              & ,zfsin, nel_kg1_k, idisp_kg1_k, mpi_double_precision, mpi_nbmx_world_k(ng_nbmx_k),ierr)
!!$         call mpi_allreduce(zfcos_mpi,zfcos,iba(ik),mpi_double_precision,mpi_sum,mpi_nbmx_world_k(ng_nbmx_k),ierr)
!!$         call mpi_allreduce(zfsin_mpi,zfsin,iba(ik),mpi_double_precision,mpi_sum,mpi_nbmx_world_k(ng_nbmx_k),ierr)
      else
         if(ipribetar >= 2) write(nfout,'(" -- ik = ",i5, " <<calc_phase_mpi>>, iba(ik) = ",i8)') ik, iba(ik)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
         do i = 1, iba(ik)
            nb = nbase(i,ik)
            ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
            zfcos(i) = dcos(ph)
            zfsin(i) = dsin(ph)
         end do
      end if
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname0)
    end subroutine calc_phase_mpi

#ifndef SX
    subroutine calc_phase_div(ik,ia)
      integer, intent(in) :: ik, ia
      integer :: i, nb, ip
      real(kind=DP) :: ph, f1, f2, f3
      integer :: id_sname0 = -1
      if(sw_timing_2ndlevel == ON) call tstatc0_begin('calc_phase_div ',id_sname0)

      f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
      do i = 1, ibl2-ibl1+1
         nb = nbase(i+ibl1-1,ik)
         ph = ngabc(nb,1)*f1+ngabc(nb,2)*f2+ngabc(nb,3)*f3
         zfcos(i) = dcos(ph)
         zfsin(i) = dsin(ph)
      end do
      if(sw_timing_2ndlevel == ON) call tstatc0_end(id_sname0)
    end subroutine calc_phase_div
#endif

#ifdef SX
    subroutine Vnonlocal_W_part_sum_over_lmt1(ibaik)
      integer, intent(in) :: ibaik
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      real(kind=DP) :: tmp
      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_over_lmt1 ',id_sname)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif

      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)
         if(il1 == il2 .and. im1 == im2) then
!!$            tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
            if(ipaw(it)==0) then
                tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
            else
                tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
            end if
         else
!!$            tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
            if(ipaw(it)==0) then
                tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
            else
                tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
            end if
         endif
         tmp = tmp * iwei(ia)
!!$         if(tmp <DELTA) cycle
         if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
         if(mdl == 0 .or. mdl == 2) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            sc(1:ibaik) = sc(1:ibaik) + tmp*zfcos(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            ss(1:ibaik) = ss(1:ibaik) - tmp*zfsin(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
         else if(mdl == 1 .or. mdl == 3) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            sc(1:ibaik) = sc(1:ibaik) - tmp*zfsin(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            ss(1:ibaik) = ss(1:ibaik) - tmp*zfcos(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
         endif
         if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
            tmp = q(lmt1,lmt2,it)*iwei(ia)
            if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if(mdl == 0 .or. mdl == 2) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qc(1:ibaik) = qc(1:ibaik) + tmp*zfcos(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qs(1:ibaik) = qs(1:ibaik) - tmp*zfsin(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
            else if(mdl == 1 .or. mdl == 3) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qc(1:ibaik) = qc(1:ibaik) - tmp*zfsin(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qs(1:ibaik) = qs(1:ibaik) - tmp*zfcos(1:ibaik)*snl(1:ibaik,lmtt1,iksnl)
            end if
         end if
      end do
      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_over_lmt1
#else
    subroutine Vnonlocal_W_part_sum_over_lmt1(ibl1,ibl2)
      integer, intent(in) :: ibl1,ibl2
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl, i
      real(kind=DP) :: tmp
      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_over_lmt1 ',id_sname)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif


      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)
         if(il1 == il2 .and. im1 == im2) then
!!$            tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
            if(ipaw(it)==0) then
                tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
            else
                tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
            end if
         else
!!$            tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
            if(ipaw(it)==0) then
                tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
            else
                tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
            end if
         endif
         tmp = tmp * iwei(ia)
!!$         if(tmp <DELTA) cycle
         if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
         if(mdl == 0 .or. mdl == 2) then
            do i = 1, ibl2-ibl1+1
               sc(i) = sc(i) + tmp*zfcos(i)*snl(i+ibl1-1,lmtt1,iksnl)
               ss(i) = ss(i) - tmp*zfsin(i)*snl(i+ibl1-1,lmtt1,iksnl)
            end do
         else if(mdl == 1 .or. mdl == 3) then
            do i = 1, ibl2-ibl1+1
               sc(i) = sc(i) - tmp*zfsin(i)*snl(i+ibl1-1,lmtt1,iksnl)
               ss(i) = ss(i) - tmp*zfcos(i)*snl(i+ibl1-1,lmtt1,iksnl)
            end do
         endif
         if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
            tmp = q(lmt1,lmt2,it)*iwei(ia)
            if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if(mdl == 0 .or. mdl == 2) then
               do i = 1, ibl2-ibl1+1
                  qc(i) = qc(i) + tmp*zfcos(i)*snl(i+ibl1-1,lmtt1,iksnl)
                  qs(i) = qs(i) - tmp*zfsin(i)*snl(i+ibl1-1,lmtt1,iksnl)
               end do
            else if(mdl == 1 .or. mdl == 3) then
               do i = 1, ibl2-ibl1+1
                  qc(i) = qc(i) - tmp*zfsin(i)*snl(i+ibl1-1,lmtt1,iksnl)
                  qs(i) = qs(i) - tmp*zfcos(i)*snl(i+ibl1-1,lmtt1,iksnl)
               end do
            end if
         end if
      end do
      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_over_lmt1
#endif
    
! =================================== added by K. Tagami ===================== 11.0
#ifdef SX
    subroutine Vnonlocal_W_part_sum_lmt1_noncl(ibaik)
      integer, intent(in) :: ibaik
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl
      complex(kind=CMPLDP) :: tmp
      real(kind=DP):: c1, c2

      integer :: id_sname = -1
      call tstatc0_begin('Vnonlocal_W_part_sum_lmt1_noncl ',id_sname)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif

      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
         if ( mdvdb == SKIP ) then
            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif
         endif
#endif
! ---------------------------------------------------- 11.0S

         tmp = dion_scr_noncl( lmt1, lmt2, ispin, ia )
         tmp = tmp * iwei(ia)

!!$         if(tmp <DELTA) cycle
         if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
         if(mdl == 0 .or. mdl == 2) then

            c1 = real(tmp);  c2 = aimag(tmp)
 
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            sc(1:ibaik) = sc(1:ibaik) + ( c1*zfcos(1:ibaik) +c2*zfsin(1:ibaik) )&
	&                               *snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            ss(1:ibaik) = ss(1:ibaik) + ( -c1*zfsin(1:ibaik) +c2*zfcos(1:ibaik) ) &
	&                               *snl(1:ibaik,lmtt1,iksnl)
         else if(mdl == 1 .or. mdl == 3) then
            c1 = real(tmp);  c2 = aimag(tmp)

#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            sc(1:ibaik) = sc(1:ibaik) +( -c1*zfsin(1:ibaik) +c2*zfcos(1:ibaik) ) &
	&                              *snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
            ss(1:ibaik) = ss(1:ibaik) + ( -c1*zfcos(1:ibaik) -c2*zfsin(1:ibaik) ) &
	&                              *snl(1:ibaik,lmtt1,iksnl)
         endif

! ----------------------------------------------------- 11.0S
!          if (mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
         if ( mdvdb == EXECUT .and. il1 == il2 ) then
            if ( SpinOrbit_mode /= BuiltIn ) then
               if ( im1 /= im2 ) cycle
            endif
! ------------------------------------------------------ 11.0S

            tmp = q_noncl(lmt1,lmt2,ispin,it)*iwei(ia)

            if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if(mdl == 0 .or. mdl == 2) then

               c1 = real(tmp);  c2 = aimag(tmp)

#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qc(1:ibaik) = qc(1:ibaik) + ( c1*zfcos(1:ibaik) +c2*zfsin(1:ibaik) )&
	&                                *snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qs(1:ibaik) = qs(1:ibaik) + ( -c1*zfsin(1:ibaik) +c2*zfcos(1:ibaik) ) &
	&                                *snl(1:ibaik,lmtt1,iksnl)

            else if(mdl == 1 .or. mdl == 3) then
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qc(1:ibaik) = qc(1:ibaik) +( -c1*zfsin(1:ibaik) +c2*zfcos(1:ibaik) ) &
	&                                *snl(1:ibaik,lmtt1,iksnl)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
               qs(1:ibaik) = qs(1:ibaik) + ( -c1*zfcos(1:ibaik) -c2*zfsin(1:ibaik) ) &
	&                                *snl(1:ibaik,lmtt1,iksnl)

            end if
         end if
      end do
      call tstatc0_end(id_sname)

    end subroutine Vnonlocal_W_part_sum_lmt1_noncl
#else
    subroutine Vnonlocal_W_part_sum_lmt1_noncl(ibl1,ibl2)
      integer, intent(in) :: ibl1,ibl2
      integer       :: lmt1, lmtt1, il1, im1, il11, mdl, i
      complex(kind=CMPLDP) :: tmp
      real(kind=DP) :: c1, c2
      integer :: id_sname = -1

      call tstatc0_begin('Vnonlocal_W_part_sum_lmt1_noncl ',id_sname)

      sc = 0.d0; ss = 0.d0
      if(mdvdb == EXECUT) then
         qc = 0.d0; qs = 0.d0
      endif


      do lmt1 = 1,ilmt(it)
         lmtt1 = lmtt(lmt1,it)
         il1   = ltp(lmt1,it)
         im1   = mtp(lmt1,it)
         il11  = il1 - 1
         mdl   = mod(il11,4)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
         if ( mdvdb == SKIP ) then
            if ( il1 /= il2 ) cycle
            if ( SpinOrbit_mode == Neglected .and. sw_hubbard == OFF ) then
               if ( im1 /= im2 ) cycle
            endif
         endif
#endif
! ---------------------------------------------------- 11.0S

         tmp = dion_scr_noncl( lmt1,lmt2,ispin,ia )
         tmp = tmp * iwei(ia)

!!$         if(tmp <DELTA) cycle
         if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
         if(mdl == 0 .or. mdl == 2) then
            c1 = real(tmp);  c2 = aimag(tmp)

            do i = 1, ibl2-ibl1+1
               sc(i) = sc(i) + ( c1*zfcos(i)+c2*zfsin(i) )*snl(i+ibl1-1,lmtt1,iksnl)
               ss(i) = ss(i) + (-c1*zfsin(i)+c2*zfcos(i) )*snl(i+ibl1-1,lmtt1,iksnl)
            end do
         else if(mdl == 1 .or. mdl == 3) then
            c1 = real(tmp);  c2 = aimag(tmp)

            do i = 1, ibl2-ibl1+1
               sc(i) = sc(i) + (-c1*zfsin(i)+c2*zfcos(i) )*snl(i+ibl1-1,lmtt1,iksnl)
               ss(i) = ss(i) + (-c1*zfcos(i)-c2*zfsin(i) )*snl(i+ibl1-1,lmtt1,iksnl)
            end do
         endif


! ----------------------------------------------------- 11.0S
!          if (mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
         if ( mdvdb == EXECUT .and. il1 == il2 ) then
            if ( SpinOrbit_mode /= BuiltIn ) then
               if ( im1 /= im2 ) cycle
            endif
! ------------------------------------------------------ 11.0S

            tmp = q_noncl(lmt1,lmt2,ispin,it)*iwei(ia)

            if(mdl == 2 .or. mdl == 3) tmp = -1*tmp
            if(mdl == 0 .or. mdl == 2) then
               c1 = real(tmp);  c2 = aimag(tmp)

               do i = 1, ibl2-ibl1+1

                  qc(i) = qc(i) + ( c1*zfcos(i)+c2*zfsin(i) )*snl(i+ibl1-1,lmtt1,iksnl)
                  qs(i) = qs(i) + (-c1*zfsin(i)+c2*zfcos(i) )*snl(i+ibl1-1,lmtt1,iksnl)
               end do
            else if(mdl == 1 .or. mdl == 3) then
               c1 = real(tmp);  c2 = aimag(tmp)

               do i = 1, ibl2-ibl1+1
                  qc(i) = qc(i) + (-c1*zfsin(i)+c2*zfcos(i) )*snl(i+ibl1-1,lmtt1,iksnl)
                  qs(i) = qs(i) + (-c1*zfcos(i)-c2*zfsin(i) )*snl(i+ibl1-1,lmtt1,iksnl)
               end do
            end if
         end if
      end do
      call tstatc0_end(id_sname)
    end subroutine Vnonlocal_W_part_sum_lmt1_noncl
#endif

! ================================================================== 11.0

    subroutine add_vnlph_l_with_eko_part(nsize,isize,vnlph)
      integer, intent(in) :: nsize, isize
      real(kind=DP), intent(inout), dimension(nsize,np_e,kimg) :: vnlph
      integer          :: ib, i
      real(kind=DP)    :: fr,fi,e
      integer :: id_sname = -1
      call tstatc0_begin('add_vnlph_l_with_eko_part ',id_sname)

      if(isize > nsize) call phase_error_with_msg(nfout, ' isize > nsize, <<add_vnlph_l_with_eko_part>>'&
                                                 , __LINE__, __FILE__)

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
         do ib = 1, np_e                                              ! MPI
            fr = fsr_l(ib,lmta2,ik);    fi = fsi_l(ib,lmta2,ik)

! ============================== modified by K. Tagami ============== 11.0
!            e = eko_l(ib,ik)
            e = eko_l(ib,ik_for_pointing_eko)
! ==================================================================== 11.0

            do i = 1, isize
               vnlph(i,ib,1) = vnlph(i,ib,1) + fr*(sc(i)-e*qc(i))-fi*(ss(i)-e*qs(i))
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do ib = 1, np_e                                              ! MPI
               fr = fsr_l(ib,lmta2,ik)
! ========================================== modified by K. Tagami ========= 11.0
!               e = eko_l(ib,ik)
               e = eko_l(ib,ik_for_pointing_eko)
! ========================================================================= 11.0

               do i = 1, isize
                  vnlph(i,ib,1) = vnlph(i,ib,1)+fr*(sc(i)-e*qc(i))
                  vnlph(i,ib,2) = vnlph(i,ib,2)+fr*(ss(i)-e*qs(i))
               end do
            end do
         else
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do ib = 1, np_e                                              ! MPI
! ====================================== modified by K. Tagami ========== 11.0
!               e = eko_l(ib,ik)
               e = eko_l(ib,ik_for_pointing_eko)
! ====================================================================== 11.0

               fr = fsr_l(ib,lmta2,ik);    fi = fsi_l(ib,lmta2,ik)
               do i = 1, isize
                  vnlph(i,ib,1) = vnlph(i,ib,1)+fr*(sc(i)-e*qc(i))-fi*(ss(i)-e*qs(i))
                  vnlph(i,ib,2) = vnlph(i,ib,2)+fi*(sc(i)-e*qc(i))+fr*(ss(i)-e*qs(i))
               end do
            end do
         end if
      end if
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_with_eko_part

    subroutine add_vnlph_l_without_eko_part(nsize,isize,vnlph)
      integer, intent(in) :: nsize, isize
      real(kind=DP), intent(inout), dimension(nsize,np_e,kimg) :: vnlph
      integer       :: ib, i
      real(kind=DP) :: fr,fi
      integer :: id_sname = -1
      call tstatc0_begin('add_vnlph_l_without_eko_part ',id_sname)
#ifndef SX
!!$      complex(kind=CMPLDP),allocatable, dimension(:) :: scss_cmplx
!!$      complex(kind=CMPLDP) :: frfi
#endif

      if(isize > nsize) call phase_error_with_msg(nfout, ' isize > nsize, <<add_vnlph_l_without_eko_part>>'&
                                                 , __LINE__, __FILE__)

      if(kimg == 1) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
         do ib = 1, np_e                                              ! MPI
            fr = fsr_l(ib,lmta2,ik); fi = fsi_l(ib,lmta2,ik)
            do i = 1, isize
               vnlph(i,ib,1) = vnlph(i,ib,1) + fr*sc(i) - fi*ss(i)
            end do
         end do
      else if(kimg == 2) then
         if(k_symmetry(ik) == GAMMA) then
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do ib = 1, np_e
               fr = fsr_l(ib,lmta2,ik)
               do i = 1, isize
                  vnlph(i,ib,1) = vnlph(i,ib,1) + fr*sc(i)
                  vnlph(i,ib,2) = vnlph(i,ib,2) + fr*ss(i)
               end do
!!$               do i = ibl1, ibl2
!!$                  vnlph_l(i,ib,1) = vnlph_l(i,ib,1) + fr*sc(i)
!!$                  vnlph_l(i,ib,2) = vnlph_l(i,ib,2) + fr*ss(i)
!!$               end do
            enddo
         else
#ifdef VPP
*vocl loop, unroll(4)
#endif
#ifdef NEC_TUNE2
!CDIR OUTERUNROLL=4
#endif
            do ib = 1, np_e
               fr = fsr_l(ib,lmta2,ik); fi = fsi_l(ib,lmta2,ik)
               do i = 1, isize
                  vnlph(i,ib,1) = vnlph(i,ib,1) + fr*sc(i) - fi*ss(i)
                  vnlph(i,ib,2) = vnlph(i,ib,2) + fi*sc(i) + fr*ss(i)
               end do
            enddo
         end if
      end if
      call tstatc0_end(id_sname)
    end subroutine add_vnlph_l_without_eko_part
  end subroutine m_ES_Vnonlocal_W
#endif


  subroutine G_dot_R_mpi(ia)
    integer, intent(in) :: ia
    integer :: i, ip
    real(kind=DP) :: grt, f1, f2, f3
!!$    integer :: id_sname = -1
!!$    call tstatc0_begin('G_dot_R_mpi ',id_sname)

    f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
    if(flag_mpi_g_dot_r) then
       zfcos_mpi(np_nbmx+1:mp_nbmx) = 0.d0; zfsin_mpi(np_nbmx+1:mp_nbmx) = 0.d0
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
       do i = 1, np_nbmx
          ip = i + ista_nbmx - 1
          grt = ngabc(ip,1)*f1+ngabc(ip,2)*f2+ngabc(ip,3)*f3
          zfcos_mpi(i) = dcos(grt)
          zfsin_mpi(i) = dsin(grt)
       end do
       call mpi_allgather(zfcos_mpi,mp_nbmx,mpi_double_precision &
            & , zfcos,mp_nbmx,mpi_double_precision, mpi_nbmx_world(ng_nbmx),ierr)
       call mpi_allgather(zfsin_mpi,mp_nbmx,mpi_double_precision &
            & , zfsin,mp_nbmx,mpi_double_precision, mpi_nbmx_world(ng_nbmx),ierr)
!!$       call mpi_allreduce(zfcos_mpi,zfcos,nbmx,mpi_double_precision,mpi_sum,mpi_nbmx_world(ng_nbmx),ierr)
!!$       call mpi_allreduce(zfsin_mpi,zfsin,nbmx,mpi_double_precision,mpi_sum,mpi_nbmx_world(ng_nbmx),ierr)
    else
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
       do i = 1, nbmx
          grt = ngabc(i,1)*f1+ngabc(i,2)*f2+ngabc(i,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    end if
!!$    call tstatc0_end(id_sname)
  end subroutine G_dot_R_mpi

  subroutine m_ES_betar_dot_WFs(nfout)
    integer, intent(in)    :: nfout
! ----------- Written by T.Yamasaki, 1 Aug. 2014 ---------->>
    logical :: flag_in_realspace

    flag_in_realspace = .false.
    if(sw_rspace==ON) then
       if(m_CtrlP_realspace_integ_OK()) flag_in_realspace = .true.
    end if

!!$    if(sw_rspace==ON)then
    if(flag_in_realspace)then
! <<--------------------------------------------------------
       call betar_dot_WFs_in_rspace(nfout)
    else
       call betar_dot_WFs_in_gspace(nfout)
    endif
  end subroutine m_ES_betar_dot_WFs

  subroutine betar_dot_WFs_in_gspace(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, iksnl, mapmode
    integer     :: id_sname = -1
    call tstatc0_begin('betar_dot_WFs (gspace) ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_betar_dot_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai(0)
    call alloc_zfsincos_mpi()
    do ia = 1, natm
! ====================================== modified by K. Tagami ============== 11.0
!       if(kv3/nspin == 1) then
!          call G_dot_R_map(ia,1)
!          mapmode = MAPPED
!       else
!          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
!          mapmode = NOTMAPPED
!       end if

       if ( noncol ) then 
          if ( kv3/ndim_spinor == 1 ) then
            call G_dot_R_map(ia,1);      mapmode = MAPPED
          else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
          endif
       else
          if ( kv3/nspin == 1 ) then
            call G_dot_R_map(ia,1);      mapmode = MAPPED
          else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
          endif
       end if
! ======================================================================= 11.0

       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI

! ========================================= modified by K. Tagami ============ 11.0
!          iksnl = (ik-1)/nspin + 1

          if ( noncol ) then
            iksnl = (ik-1)/ndim_spinor + 1
          else
            iksnl = (ik-1)/nspin + 1
          endif
! ============================================================================ 11.0
          call m_ES_betar_dot_WFs_4_lmta_k(ista_k,iend_k,ik,zaj_l,ia,iksnl,snl,fsr_l,fsi_l,mapmode)
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai()
    call m_ES_dealloc_zfsincos()
    if(ipribetar >= 2) then
       write(nfout,'(" --- fsr_l, fsi_l ---")')
       do ik = ista_k, iend_k                              ! MPI
          call wd_fsr_fsi(ista_k,iend_k,ik,fsr_l,fsi_l)    ! MPI
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine betar_dot_WFs_in_gspace

  subroutine m_ES_wd_fsr_fsi()
    integer :: ik
    if(ipribetar >= 2) then
       write(nfout,'(" --- fsr_l, fsi_l ---")')
       do ik = ista_k, iend_k                              ! MPI
          call wd_fsr_fsi(ista_k,iend_k,ik,fsr_l,fsi_l)    ! MPI
       end do
    end if
  end subroutine m_ES_wd_fsr_fsi

  subroutine wd_fsr_fsi(k1,k2,ik,bwr_l,bwi_l)
    integer, intent(in) :: k1,k2,ik
    real(kind=DP), intent(in) :: bwr_l(np_e,nlmta,k1:k2),bwi_l(np_e,nlmta,k1:k2)
    integer :: it, ia
    if(k_symmetry(ik) == GAMMA) then
       write(nfout,'(" -- fsr --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (bwr_l(map_z(it),ia,ik),ia=1,nlmta)!MPI
       end do
    else
       write(nfout,'(" -- fsr, fsi --")')
       do it = 1, neg
          if(map_e(it) /= myrank_e) cycle                  ! MPI
          write(nfout,'(" ik = ",i8," ib = ",i8)') ik, it
          write(nfout,'(6d12.4)') (bwr_l(map_z(it),ia,ik),bwi_l(map_z(it),ia,ik),ia=1,nlmta)!MPI
       end do
    end if
    call flush(nfout)
  end subroutine wd_fsr_fsi

  subroutine G_dot_R_map(ia,ik)
    integer, intent(in) :: ia,ik
    integer :: i, i1
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1

    call tstatc0_begin('G_dot_R_map ',id_sname)
    f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
    if(.not.allocated(nbase)) call phase_error_with_msg(nfout, ' nbase is not allocated'&
                                                       , __LINE__, __FILE__)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
    if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       zfcos = 0.0D0;   zfsin = 0.0D0
       do i = 1, iba(ik)
          i1 = nbase(i,ik)
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
       do i = 1, iba(ik)
          i1 = nbase(i,ik)
          grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
          zfcos(i) = dcos(grt)
          zfsin(i) = dsin(grt)
       end do
    endif
    call tstatc0_end(id_sname)
  end subroutine G_dot_R_map


  subroutine G_dot_R_map_div(ia,ik,ibl1,ibl2)
    integer, intent(in) :: ia,ik,ibl1,ibl2
    integer :: i, i1
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1
    call tstatc0_begin('G_dot_R_map_div ',id_sname)

    f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
    if(.not.allocated(nbase)) call phase_error_with_msg(nfout, ' nbase is not allocated',__LINE__,__FILE__)

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
  subroutine G_dot_R_map_blk(ik,ibl1,ibl2,ia1,ia2)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
    integer, intent(in) :: ik,ibl1,ibl2, ia1,ia2
    integer :: ia, i, i1
    real(kind=DP) :: grt, f1, f2, f3
    integer :: id_sname = -1

    if(.not.allocated(nbase)) call phase_error_with_msg(nfout, ' nbase is not allocated'&
                                                       , __LINE__, __FILE__)
    call tstatc0_begin('G_dot_R_map_blk ',id_sname)
    !$acc data present(GVec_on_refcell, nbase, wk_zfcos, wk_zfsin, pos, ngabc)
    !$acc kernels
    if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       wk_zfcos = 0.0D0;    wk_zfsin = 0.0D0
       do i = ibl1, ibl2
          i1 = nbase(i,ik)
          if ( sw_force_kpt_inside_bz == ON ) then
             if ( GVec_on_refcell(i1,ik) == 0 ) cycle
          else
             if ( GVec_on_refcell(i1,1) == 0 ) cycle
          endif
!!$       do ia=1,natm
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
          i1 = nbase(i,ik)
!!$       do ia=1,natm
          do ia=ia1,ia2
             f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
             grt = ngabc(i1,1)*f1+ngabc(i1,2)*f2+ngabc(i1,3)*f3
             wk_zfcos(i-ibl1+1,ia-ia1+1) = dcos(grt)
             wk_zfsin(i-ibl1+1,ia-ia1+1) = dsin(grt)
          end do
       end do
    endif
    !$acc end kernels
    !$acc end data
    call tstatc0_end(id_sname)
  end subroutine G_dot_R_map_blk
#endif
!! #ifdef NONLOCAL_DGEMM end


  subroutine m_ES_betar_dot_WFs_4_each_k(nfout,ik)
    integer, intent(in) :: nfout,ik
    call m_ES_betar_dot_Psi_4_each_k(nfout,zaj_l,ista_k,iend_k,ik,fsr_l,fsi_l)
  end subroutine m_ES_betar_dot_WFs_4_each_k
#ifdef NONLOCAL_DGEMM
  subroutine m_ES_betar_dot_Psi_4_each_k(nfout,psi_l,k1,k2,ik,bpr_l,bpi_l,map)
    integer, intent(in) :: nfout,k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l !MPI
    logical, intent(in), optional, dimension(np_e) :: map
    if(present(map)) then
      call m_ES_betar_dot_Psi_4_each_k_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l,map=map)
    else
      call m_ES_betar_dot_Psi_4_each_k_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l)
    endif
  end subroutine m_ES_betar_dot_Psi_4_each_k

  subroutine m_ES_betar_dot_Psi_4_each_k_snl(nfout,snl,psi_l,k1,k2,ik,bpr_l,bpi_l,precalculate_phase,map)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, October 2009 : Mblock
!
    integer, intent(in) :: nfout,k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l !MPI
    logical, intent(in), optional :: precalculate_phase
    logical, intent(in), optional, dimension(np_e) :: map
    integer :: iksnl, ibsize

    integer :: ia, it, lmt1, msize, msize_target, msizemax, msizesum, natm_redmax
    integer :: ibl1,ibl2, ia1, ia2
    integer :: LD11,  LD12
    logical :: tran1
!!$    integer, allocatable, dimension(:) :: mil
    real(kind=DP), allocatable, dimension(:,:,:) :: psi_ri
    logical :: prec
#ifdef _CUDA_
    real(kind=DP),allocatable,dimension(:,:) :: bprt, bpit
#endif

    integer :: id_sname = -1
! ----------- Written by T.Yamasaki, 1 Aug. 2014 ---------->>
    logical :: flag_in_realspace

    prec = sw_betar_dot_wfs_exp==ON
    if(present(precalculate_phase)) then
      prec = precalculate_phase
    endif

    flag_in_realspace = .false.
    if(sw_rspace==ON) then
       if(m_CtrlP_realspace_integ_OK()) flag_in_realspace = .true.
    end if

!!$    if(sw_rspace==ON)then
    if(flag_in_realspace)then
! <<--------------------------------------------------------
       bpr_l(:,:,ik)=0.d0;if(k_symmetry(ik)/=GAMMA) bpi_l(:,:,ik)=0.d0
       call betar_dot_Psi_4_each_k_in_rs(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l)
       if(ipribetar >= 2) call wd_fsr_fsi(k1,k2,ik,bpr_l,bpi_l)
    else

    if (sw_betar_dot_wfs_exp==ON .and. prec) then
      if(present(map)) then
        call m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,bpr_l,bpi_l,map)
      else
        call m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,bpr_l,bpi_l)
      endif
      return
    endif
    call tstatc0_begin('betar_dot_Psi ',id_sname,level=1)

! ======================================= modified by K. Tagami ============= 11.0
!    iksnl = (ik-1)/nspin + 1
!
     if ( noncol ) then
       iksnl = (ik-1)/ndim_spinor + 1
     else
       iksnl = (ik-1)/nspin + 1
     endif

! =========================================================================== 11.0

!!$    if(ik < ista_k .or. iend_k < ik) then
    if(ik < k1 .or. k2 < ik) then
       call phase_error_with_msg(nfout, ' ik < ista_k .or. iend_k <ik (m_ES_betar_dot_Psi_4_each_k)',__LINE__,__FILE__)
    end if

    bpr_l(:,:,ik) = 0.d0
    if(k_symmetry(ik) /= GAMMA) bpi_l(:,:,ik) = 0.d0

    call set_msize() ! -> ibsize, msize_target, msizemax, natm_redmax

    call alloc_wkzfsincos_red(ibsize,natm_redmax)
    call alloc_wkarai(ibsize,msizemax)
!!$    call alloc_wkother(msizemax)
!!$    allocate(mil(msizemax))
#ifdef _CUDA_
    allocate(bprt(np_e,nlmta))
    allocate(bpit(np_e,nlmta))
#endif
    allocate(psi_ri(ibsize,np_e,kimg))
    ia1 = 1
    msizesum = 0
    !$acc data copyin(psi_l, k_symmetry, wk_zfsin, wk_zfcos, &
    !$acc             GVec_on_refcell, pos) &
    !$acc      create( wk_ar, wk_ai, psi_ri) &
    !$acc        copyout(bprt, bpit) &
    !$acc     present(snl, ngabc, nbase, ityp, ilmt, ltp, lmtt)
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

       !$acc data create(bp_tmp1, bp_tmp2)
       !$acc kernels
       bp_tmp1=0.d0
       if(k_symmetry(ik) /= GAMMA) bp_tmp2=0.d0
       !$acc end kernels

       do ibl1 = 1, iba(ik), ibsize
          ibl2 = min( ibl1+ibsize-1, iba(ik) )
          call cppsi_l_to_psi_ri(ibl1,ibl2)
          call G_dot_R_map_blk(ik,ibl1,ibl2, ia1,ia2) ! pos(ia,1:3), ngabc -> wk_zfcos,wk_zfsin
          call m_ES_betar_dot_WFs_4_lmta_k_blk(ik,ibsize,ibl1,ibl2,psi_ri &
               & ,iksnl,snl, msize,LD11,tran1,ia1,ia2) !  --> bp_tmp1, bp_tmp2
       end do

       if(k_symmetry(ik) == GAMMA) then
#ifdef _CUDA_
          call post_lmta_k_blk_GAMMA(k1,k2,bprt)
#else
          call post_lmta_k_blk_GAMMA(k1,k2,bpr_l)
#endif
       else
#ifdef _CUDA_
          call post_lmta_k_blk(k1,k2, bprt,bpit)
#else
          call post_lmta_k_blk(k1,k2, bpr_l,bpi_l)
#endif
       end if
       !$acc end data
       msizesum = msizesum+msize

       call dealloc_wkbp(ik)
       ia1 = ia2+1
       if(ia1 > natm) exit Mblock
    end do Mblock
    !$acc end data

    if(ipribetar >= 2) call wd_fsr_fsi(k1,k2,ik,bpr_l,bpi_l)
#ifdef _CUDA_
    bpr_l(:,:,ik) = bprt(:,:)
    if(k_symmetry(ik)/=GAMMA) bpi_l(:,:,ik) = bpit(:,:)
    deallocate(bprt, bpit)
#endif
    deallocate(psi_ri)
!!$    call dealloc_wkother()
    call dealloc_wkarai()
    call dealloc_wkzfsincos()

    call tstatc0_end(id_sname)

    endif

   contains
     subroutine cppsi_l_to_psi_ri(ibl1,ibl2)
       integer, intent(in) :: ibl1,ibl2
       integer :: ip, n, ri
       !$acc data present(psi_ri, psi_l)
       !$acc kernels
       do ri = 1, kimg
          do n = 1, np_e
             do ip = ibl1, ibl2
                psi_ri(ip-ibl1+1,n,ri) = psi_l(ip,n,ik,ri)
             end do
          end do
       end do
       !$acc end kernels
       !$acc end data
     end subroutine cppsi_l_to_psi_ri

     subroutine set_msize()

       integer :: ncache, N, K, M, ia_i, ia_f, msize
       integer, save :: msize_print = OFF

       if( nblocksize_betar_is_given) then
          ibsize = nblocksize_betar_dot_wfs
       else
          ibsize = nb_betar_default
          if(ipribetar >= 2 .and. msize_print == OFF) then
             write(nfout,'(" !ibsize(=nblocksize_betar_w) (m_ES_betar_dot_WFs_4_each_k) = ",i8)') ibsize
          end if
       end if
! === DEBUG by tkato 2011/06/29 ================================================
       if( ibsize == 0 ) ibsize = iba(ik)
! ==============================================================================
!!$    ibsize=10000
       ncache = (m_CtrlP_cachesize()*1024)*3/4
       if(ncache == 0) then
          M = nlmta
       else
          K = ibsize
          N = np_e
#ifdef _BLOCK_NLMTA_
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
#else
          M = nlmta
#endif
       endif
       msize_target = M
       if(msize_target < ilmt(1)) msize_target = ilmt(1)
       if(ipribetar>=2) then
          if(msize_print == OFF) then
!            write(nfout,'(" mype = ",i3, " np_e, msize_target, nlmta, ibsize, iba(ik)/ibsize = ",5i8 &
!                 & ," <<m_ES_betar_dot_WFs_4_each_k>>")') &
!                 & mype, np_e, msize_target, nlmta, ibsize, iba(ik)/ibsize
          end if
       end if

       if(ncache == 0) then
          msizemax = nlmta
          natm_redmax = natm
       else
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
       end if
       msize_print = ON

     end subroutine set_msize

     subroutine pre_lmta_k_blk(msize,LD11,LD12)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
       integer,intent(out) :: msize,LD11, LD12
       integer :: icnt, ia, it, lmt1

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
       if(np_e .gt. msize) then
          LD11 = np_e ; LD12 = msize
          tran1 = .true.
       else
          LD11 = msize ; LD12 = np_e
          tran1 = .false.
       end if
#ifdef SX
       if( mod(LD11,2) .eq. 0) LD11 = LD11+1
#endif
     end subroutine pre_lmta_k_blk

#ifdef _CUDA_
    subroutine post_lmta_k_blk_Gamma( k1,k2, bprt)
#else
    subroutine post_lmta_k_blk_Gamma( k1,k2, bpr_l)
#endif
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, Octover 2009 : multipication of i^L is removed
      integer,intent(in)  :: k1,k2
#ifdef _CUDA_
      real(kind=DP),intent(inout),dimension(np_e,nlmta)  ::  bprt ! MPI
#else
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2)  ::  bpr_l ! MPI
#endif
      integer :: j,ib

      if(msize+msizesum > nlmta) then
         write(nfout,'(" msize+msizesum > nlmta")')
         if(tran1) write(nfout,'(" tran1 = .true.")')
         if(.not.tran1) write(nfout,'(" tran1 = .false.")')
!!$         call mpi_stop(__FILE__,__LINE__,"post_lmta_k_blk_Gamma")
         call phase_error_with_msg(nfout, ' msize+msizesum > nlmta',__LINE__,__FILE__)
      end if
      !$acc data present(bprt, bp_tmp1)
      !$acc kernels
      if(tran1) then
!cdir nodep
         do j = 1, msize
            do ib = 1, np_e
#ifdef _CUDA_
               bprt(ib,j+msizesum) = bp_tmp1(ib,j)
#else
               bpr_l(ib,j+msizesum,ik) = bp_tmp1(ib,j)
#endif
            enddo
         enddo
      else
!cdir nodep
         do j = 1, msize
            do ib = 1, np_e
#ifdef _CUDA_
               bprt(ib,j+msizesum) = bp_tmp1(j,ib)
#else
               bpr_l(ib,j+msizesum,ik) = bp_tmp1(j,ib)
#endif
            enddo
         enddo
      endif
      !$acc end kernels
      !$acc end data
    end subroutine post_lmta_k_blk_Gamma

#ifdef _CUDA_
    subroutine post_lmta_k_blk( k1,k2, bprt,bpit)
#else
    subroutine post_lmta_k_blk( k1,k2, bpr_l,bpi_l)
#endif
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, Octover 2009 : multipication of i^L is removed
      integer,intent(in)  :: k1,k2
#ifdef _CUDA_
      real(kind=DP),intent(inout),dimension(np_e,nlmta)  ::  bprt,bpit  ! MPI
#else
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2)  ::  bpr_l,bpi_l  ! MPI
#endif
      integer :: j,ib

      !$acc data present(bprt, bpit, bp_tmp1, bp_tmp2)
      !$acc kernels
      if(k_symmetry(ik) == GAMMA) then  !! betar_dot_WFs_core2_blk
         if(tran1) then
!cdir nodep
            do j = 1, msize
               do ib = 1, np_e
#ifdef _CUDA_
                  bprt(ib,j+msizesum) = bp_tmp1(ib,j)
#else
                  bpr_l(ib,j+msizesum,ik) = bp_tmp1(ib,j)
#endif
               enddo
            enddo
         else
!cdir nodep
            do j = 1, msize
               do ib = 1, np_e
#ifdef _CUDA_
                  bprt(ib,j+msizesum) = bp_tmp1(j,ib)
#else
                  bpr_l(ib,j+msizesum,ik) = bp_tmp1(j,ib)
#endif
               enddo
            enddo
         endif
      else
         if(tran1) then
!cdir nodep
            do j = 1, msize
               do ib = 1, np_e
#ifdef _CUDA_
                  bprt(ib,j+msizesum) =  bp_tmp1(ib,j)
                  bpit(ib,j+msizesum) =  bp_tmp2(ib,j)
#else
                  bpr_l(ib,j+msizesum,ik) =  bp_tmp1(ib,j)
                  bpi_l(ib,j+msizesum,ik) =  bp_tmp2(ib,j)
#endif
               end do
            end do
          else
!cdir nodep
            do j = 1, msize
               do ib = 1, np_e
#ifdef _CUDA_
                  bprt(ib,j+msizesum) =  bp_tmp1(j,ib)
                  bpit(ib,j+msizesum) =  bp_tmp2(j,ib)
#else
                  bpr_l(ib,j+msizesum,ik) =  bp_tmp1(j,ib)
                  bpi_l(ib,j+msizesum,ik) =  bp_tmp2(j,ib)
#endif
               end do
            end do
          endif
       end if
      !$acc end kernels
      !$acc end data
    end subroutine post_lmta_k_blk

  end subroutine m_ES_betar_dot_Psi_4_each_k_snl
#else
!! #ifdef NONLOCAL_DGEMM else
  subroutine m_ES_betar_dot_Psi_4_each_k(nfout,psi_l,k1,k2,ik,bpr_l,bpi_l)
    integer, intent(in)    :: nfout, k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l !MPI
    
    integer :: ia, iksnl, mapmode
#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache, ie, ig, ri
    real(kind=DP),allocatable,dimension(:,:,:) :: zaj_red ! d(ibsize,np_e,kimg)
#endif
    integer :: id_sname = -1
    call tstatc0_begin('betar_dot_Psi ',id_sname,level=1)

! ====================================== modified by K. Tagami ============ 11.0
!    iksnl = (ik-1)/nspin + 1
    if ( noncol ) then
      iksnl = (ik-1)/ndim_spinor + 1
    else
      iksnl = (ik-1)/nspin + 1
    endif
! ========================================================================== 11.0

!!$    call alloc_zfsincos_mpi_k()
#ifdef SX
    call m_ES_alloc_zfsincos(ik)
    call m_ES_alloc_arai(ik)

    do ia = 1, natm
       if(ipribetar >= 2) write(nfout,'("(m_ES_betar_dot_Psi_4_each_k) ia = ", i5)') ia
       call G_dot_R_map(ia,ik) ! pos(ia,1:3), ngabc -> zfcos,zfsin
       mapmode = MAPPED
!!$       call G_dot_R_mpi_k(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
!!$       mapmode = NOTMAPPED
!!$       call m_ES_betar_dot_WFs_4_lmta_k(ista_k,iend_k,ik,zaj_l,ia,iksnl,snl,fsr_l,fsi_l,mapmode)
       call m_ES_betar_dot_WFs_4_lmta_k(k1,k2,ik,psi_l,ia,iksnl,snl,bpr_l,bpi_l,mapmode)
       !                                             ->fsr_l,fsi_l
    end do
    if(ipribetar >= 2) call wd_fsr_fsi(k1,k2,ik,bpr_l,bpi_l)
    call m_ES_dealloc_arai
    call m_ES_dealloc_zfsincos
!!$    call dealloc_zfsincos_mpi_k()
#else
! NEC tune ------------------------------->
! Revised by T. Yamasaki, September 2008 
    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) then
       ibsize = iba(ik)
    else
       if(k_symmetry(ik) == GAMMA) then ! core2
          if(kimg == 1) then ! ar(i):1,psi_l(i,ib,ik,1):np_e or
             !                 ar(i):1,psi_l(i,ib,ik,2):np_e
             ibsize=ncache/(8*(np_e+1))
          else !  ar(i),ai(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,2):np_e*2  or
             !    ai(i),ar(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,2):np_e*2
             ibsize=ncache/(8*(2+np_e*2))
          endif
       else ! ar(i),ai(i):2,psi_l(i,ib,ik,1),psi_l(i,ib,ik,1):np_e*2
          ibsize=ncache/(8*(2+np_e*2))
       end if
    end if

    allocate(zaj_red(ibsize,np_e,kimg))
    call alloc_zfsincos(ibsize)
    call alloc_arai(ibsize)

    bpr_l(:,:,ik) = 0.d0
    if(k_symmetry(ik) /= GAMMA) bpi_l(:,:,ik) = 0.d0

    do ibl1 = 1, iba(ik), ibsize
       ibl2 = min( ibl1+ibsize-1, iba(ik))
       do ri = 1, kimg
          do ie = 1, np_e
             do ig = 1, ibl2-ibl1+1
                zaj_red(ig,ie,ri) = psi_l(ig+ibl1-1,ie,ik,ri)
             end do
          end do
       end do

       do ia = 1, natm
          if(ipribetar >= 2) write(nfout,'("(m_ES_betar_dot_Psi_4_each_k) ia = ", i5)') ia
          call G_dot_R_map_div(ia,ik,ibl1,ibl2) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = MAPPED
          call m_ES_betar_dot_WFs_4_lmta_k_div(k1,k2,ik,ibsize,ibl1,ibl2,zaj_red &
               & ,ia,iksnl,snl,bpr_l,bpi_l,mapmode)   !      ->fsr_l,fsi_l
       end do
    end do
    if(ipribetar >= 2) call wd_fsr_fsi(k1,k2,ik,bpr_l,bpi_l)
    call dealloc_arai()
    call dealloc_zfsincos()
    deallocate(zaj_red)
#endif
    call tstatc0_end(id_sname)
  end subroutine m_ES_betar_dot_Psi_4_each_k
#endif
!! #ifdef NONLOCAL_DGEMM end

!kukan 4

  subroutine m_ES_betar_dot_Px_4_each_k(psi_l,k1,k2,ik,bpr_l,bpi_l)
    integer, intent(in)                                     :: k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l !MPI

    integer :: ia, iksnl, mapmode
    integer :: id_sname = -1
    call tstatc0_begin('betar_dot_Psi ',id_sname,level=1)

    call m_ES_alloc_zfsincos(ik)
    call m_ES_alloc_arai(ik)

! ============================== modified by K. Tagami ===================== 11.0
!    iksnl = (ik-1)/nspin + 1
    if ( noncol ) then
      iksnl = (ik-1)/ndim_spinor + 1
    else
      iksnl = (ik-1)/nspin + 1
    endif
! ========================================================================= 11.0

    do ia = 1, natm
!!$       call G_dot_R(natm,ia,pos,kgp,nbmx,ngabc,zfcos,zfsin) ! -(bottom_Subroutines)
!!$       mapmode = NOTMAPPED
       call G_dot_R_map(ia,ik)
       mapmode = MAPPED
       call m_ES_betar_dot_WFs_4_lmta_k(k1,k2,ik,psi_l,ia,iksnl,snl,bpr_l,bpi_l,mapmode)
    end do
    call m_ES_dealloc_arai
    call m_ES_dealloc_zfsincos

    call tstatc0_end(id_sname)
  end subroutine m_ES_betar_dot_Px_4_each_k

  subroutine m_ES_betar_dot_WFs_4_lmta_k(k1,k2,ik,psi_l,ia,iksnl,snl_or_snld,bpr_l,bpi_l,mapmode)
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l        ! MPI
    real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) ::    bpr_l,bpi_l  ! MPI
    real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld

    integer                    ::it, lmt1, lmtt1, lmta1, il1
    integer :: id_sname = -1
#ifndef SX
! NEC tune ------------------------------->
    integer :: ibl1,ibl2,ibsize,ncache
    ncache = (m_CtrlP_cachesize()*1024)*3/4
    if(ncache == 0) then
       ibsize = iba(ik)
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

    call tstatc0_begin('m_ES_betar_dot_WFs_4_lmta_k ',id_sname)

    it    = ityp(ia)
#ifndef SX
! NEC tune ------------------------------->
    do ibl1=1,iba(ik),ibsize
    ibl2=ibl1+ibsize-1
    if(ibl2.gt.iba(ik)) ibl2=iba(ik)
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
          if(ibl2.eq.iba(ik))then
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
      real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
      integer :: i, i1
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
                  if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
                     ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
                  endif
                  do i = ibl1, ibl2
                     ar(i) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
                     ai(i) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
                  end do
               else if(mil == 2 .or. mil == 0) then
                  if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
                     ar(1) = 0.d0
                  endif
                  do i = ibl1, ibl2
                     ai(i) = -zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
                     ar(i) =  zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
                  end do
               end if
            else
               if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
                  ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
               endif
               do i = ibl1, ibl2
                  ar(i) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
                  ai(i) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
               end do
            end if
         else
            if(kimg==2 .and. k_symmetry(ik) == GAMMA) then
               if(mil == 1 .or. mil ==3) then
                  if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
                     i1    = nbase(1,ik)
                     ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
                  endif
                  do i = ibl1, ibl2
                     i1    = nbase(i,ik)
                     ar(i) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
                     ai(i) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
                  end do
               else if(mil == 2 .or. mil == 0) then
                  ar(1) = 0.d0
                  do i = ibl1, ibl2
                     i1    = nbase(i,ik)
                     ai(i) = -zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
                     ar(i) =  zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
                  end do
               else
               end if
            else
               if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
                  i1    = nbase(1,ik)
                  ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
               endif
               do i = ibl1, ibl2
                  i1    = nbase(i,ik)
                  ar(i) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
                  ai(i) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
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
            do i = 1, iba(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
               ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               ar(i) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
               ai(i) = zfsin(i)*snl_or_snld(i,lmtt1,iksnl)
            end do
         else
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
               i1    = nbase(1,ik)
               ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               i1    = nbase(i,ik)
               ar(i) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
               ai(i) = zfsin(i1)*snl_or_snld(i,lmtt1,iksnl)
            end do
         end if
#endif
         if(ipribetar>=3 ) then
            write(nfout,'(" --- ar --- <<G_dot_R_mult_snl>>")')
            write(nfout,'(5f20.10)') (ar(i),i=1,iba(ik))
            write(nfout,'(" --- ai --- <<G_dot_R_mult_snl>>")')
            write(nfout,'(5f20.10)') (ai(i),i=1,iba(ik))
         end if
      else  ! if(kimg == 1 .and. k_symmetry == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
               ar(1) = zfcos(1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               ar(i) = zfcos(i)*snl_or_snld(i,lmtt1,iksnl)
            end do
         else
#ifdef SX
! NEC tune ------------------------------->
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
#else
            if(ibl1.ne.1 .and. ibl2.eq.iba(ik))then
               i1    = nbase(1,ik)
               ar(1) = zfcos(i1)*snl_or_snld(1,lmtt1,iksnl)
            endif
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = ibl1, ibl2
! NEC tune <-------------------------------
#endif
               i1    = nbase(i,ik)
               ar(i) = zfcos(i1)*snl_or_snld(i,lmtt1,iksnl)
            end do
         end if
      end if
      call tstatc0_end(id_sname)
    end subroutine G_dot_R_mult_snl

    subroutine betar_dot_WFs_core(psi_l,bpr_l,bpi_l)
      real(kind=DP),intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l        ! MPI
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) ::    bpr_l,bpi_l  ! MPI
      real(kind=DP) :: bpr, bpi
      integer       :: ib, i

      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core ',id_sname)

#ifndef SX
! NEC tune --------------------------->
!      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
!      bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      if(ibl1.eq.1)then
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
         do ib = 1, np_e                   ! MPI
#ifdef SX
! NEC tune ------------------------------------------>
            bpr = 0.d0; bpi = 0.d0
            do i = 1, iba(ik)
#else
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = ibl1, ibl2
! NEC tune <------------------------------------------
#endif
               bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
               bpi = bpi + ai(i)*psi_l(i,ib,ik,1)
!!$               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)
!!$               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
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
         do ib = 1, np_e                   ! MPI
#ifdef SX
! NEC tune ------------------------------------------>
            bpr = 0.d0; bpi = 0.d0
            do i = 1, iba(ik)
#else
            bpr = bpr_l(ib,lmta1,ik)
            bpi = bpi_l(ib,lmta1,ik)
            do i = ibl1, ibl2
! NEC tune <------------------------------------------
#endif
               bpr = bpr + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
               bpi = bpi + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
!!$               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
!!$               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
            end do
            bpr_l(ib,lmta1,ik) = bpr
            bpi_l(ib,lmta1,ik) = bpi
         end do
!!$#else
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

    subroutine betar_dot_WFs_core2(il1,psi_l,bpr_l)
      integer, intent(in) :: il1
      real(kind=DP),intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l        ! MPI
      real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) ::    bpr_l  ! MPI
      integer       :: ib, i, mil
      real(kind=DP) :: bpr, fp

#ifndef SX
! NEC tune ------------------------------->
      integer :: ibl1_2
#endif
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core2 ',id_sname)

#ifdef SX
      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
#else
      if(ibl1.eq.1) then
        ibl1_2=2
        bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      else
        ibl1_2=ibl1
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
               do i = 2, iba(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
#else
               if(ibl2.eq.iba(ik))then
                 bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
               else
                 bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
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
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = 2, iba(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
#else
               if(ibl2.eq.iba(ik))then
                 bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
               else
                 bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
         end if
      else if(kimg == 2) then
#ifdef BETAR_GAMMA_TUNE

            do ib = 1, np_e
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1) - ai(i)*psi_l(i,ib,ik,2)
               end do
               if(ibl2.eq.iba(ik))then
                 bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
               else
                 bpr_l(ib,lmta1,ik) = bpr
               endif

            end do
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
               do i = 2, iba(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  bpr = bpr + ar(i)*psi_l(i,ib,ik,1) - ai(i)*psi_l(i,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
#else
               if(ibl2.eq.iba(ik))then
                 bpr_l(ib,lmta1,ik) = fp*(2.d0*bpr + psi_l(1,ib,ik,1)*ar(1))
               else
                 bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
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
#ifdef SX
! NEC tune ------------------------------->
               bpr = 0.d0
               do i = 2, iba(ik)
#else
               bpr = bpr_l(ib,lmta1,ik)
               do i = ibl1_2, ibl2
! NEC tune <-------------------------------
#endif
                  bpr = bpr + ai(i)*psi_l(i,ib,ik,1) + ar(i)*psi_l(i,ib,ik,2)
               end do
#ifdef SX
! NEC tune ------------------------------->
               bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
#else
               if(ibl2.eq.iba(ik))then
                 bpr_l(ib,lmta1,ik) = fp*2.d0*bpr
               else
                 bpr_l(ib,lmta1,ik) = bpr
               endif
! NEC tune <-------------------------------
#endif
            end do
         end if
#endif
      end if
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core2
  end subroutine m_ES_betar_dot_WFs_4_lmta_k

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
         call phase_error_with_msg(nfout, ' mapmode /= MAPPED <<G_dot_R_mult_snl in m_ES_betar_dot_WFs_4_lmta_k_div>>'&
                                  , __LINE__, __FILE__)
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
  subroutine m_ES_betar_dot_WFs_4_lmta_k_blk(ik,ibsize,ibl1,ibl2,psi_l,iksnl,snl_or_snld &
       &                               ,msize,LD11, tran1, ia1, ia2 )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
    integer, intent(in) :: ik, ibsize,ibl1,ibl2, iksnl
    integer, intent(in) :: msize,LD11,ia1,ia2
    logical, intent(in) :: tran1
    real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) ::  psi_l
    real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) ::  snl_or_snld
    integer, allocatable, dimension(:,:) :: lmtmap
    real(kind=DP), allocatable, dimension(:,:) :: fpar
    integer :: ia, it, lmt1, j, mil0, lt
    real(kind=DP) :: fp, fp0
    integer :: id_sname = -1
     
!!$    integer :: id_sname2 = -1
    call tstatc0_begin('m_ES_betar_dot_WFs_4_lmta_k_blk ',id_sname)

!!$    if(npes > 1) then
!!$       call tstatc0_begin('mpi_barrier(betar_dot_WFs) ',id_sname2)
!!$       call mpi_barrier(mpi_k_world(myrank_k),ierr)
!!$       call tstatc0_end(id_sname2)
!!$    end if
      allocate(lmtmap(nlmt,natm))
      allocate(fpar(nlmt,natm))
      if(k_symmetry(ik) == GAMMA) then
         fp0 = 2.d0
      else
         fp0 = 1.d0
      end if
      j = 0
      do ia = 1, ia2-ia1+1
         it = ityp(ia+ia1-1)
         do lmt1 = 1, ilmt(it)
            j = j + 1
            lmtmap(lmt1, ia) = j
            mil0 = mod(ltp(lmt1,it),4)
            if(mil0 == 0 .or. mil0 == 1) then
               fpar(lmt1, ia) = fp0
            else if(mil0 == 2 .or. mil0 == 3) then
               fpar(lmt1, ia) = -fp0
            end if
         enddo
       enddo
      !$acc data present(psi_l, k_symmetry, ityp, ilmt, ltp, lmtt, wk_zfsin, wk_zfcos, snl_or_snld, &
      !$acc      wk_ar, wk_ai, bp_tmp1, bp_tmp2) &
      !$acc      copyin(lmtmap, fpar)
       call G_dot_R_mult_snl_blk(snl_or_snld)      ! lmtt1, exp(iGR), zfcos,zfsin,snl_or_snld ->ar,ai
       if(k_symmetry(ik) == GAMMA) then
          call betar_dot_WFs_core2_blk(ibsize,tran1,psi_l )  ! lmta1, ar,ai,psi_l -> bp_tmp1
       else
          call betar_dot_WFs_core_blk (ibsize,tran1, psi_l )  ! lmta1, sum(c(k+G)exp(iGR)*snl() ,ar,ai,psi_l -> bp_tmp1, bp_tmp2
       end if

       !$acc end data

       deallocate(lmtmap)
       deallocate(fpar)

    call tstatc0_end(id_sname)
  contains
    subroutine G_dot_R_mult_snl_blk(snl_or_snld)
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
! Revised by T. Yamasaki, Octover 2009 : mapping algorithm with the subroutine "post_lmta_k_blk"
!
      real(kind=DP),intent(in), dimension(kg1,nlmtt,ista_snl:iend_snl) :: snl_or_snld
      integer :: j, i
      integer :: mil0, ia, it, lmt1, lt
      real(kind=DP) :: fp, fp0
      integer :: id_sname = -1
      call tstatc0_begin('G_dot_R_mult_snl_blk ',id_sname)

      !$acc data present(k_symmetry, ityp, ilmt, ltp, lmtt, wk_ar, wk_ai, wk_zfsin, wk_zfcos, snl_or_snld, lmtmap, fpar)
      !$acc kernels
      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
!!$         if(.not.allocated(wk_ar)) stop 'wk_ar is not allocated'
!!$         if(.not.allocated(wk_ai)) stop 'wk_ai is not allocated'
!!$         if(.not.allocated(wk_zfsin)) stop 'wk_zfsin is not allocated'
!!$         if(.not.allocated(wk_zfcos)) stop 'wk_zfcos is not allocated'

         if(k_symmetry(ik) == GAMMA) then
            fp0 = 2.d0
         else
            fp0 = 1.d0
         end if

         j = 0
         !$acc loop independent
         do ia = 1, ia2-ia1+1
            it = ityp(ia+ia1-1)
            do lmt1 = 1, ilmt(it)
               !j = j + 1
               j = lmtmap(lmt1, ia)
               mil0 = mod(ltp(lmt1,it),4)
               lt = lmtt(lmt1,it)
!               if(mil0 == 0 .or. mil0 == 1) then
!                  fp = fp0
!               else if(mil0 == 2 .or. mil0 == 3) then
!                  fp = -fp0
!               end if
               fp = fpar(lmt1, ia)

               if(k_symmetry(ik) == GAMMA) then
                  if(mil0 == 0 .or. mil0 == 2) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) = fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                     if(ibl1 == 1) then
                        wk_ar(1, j) = 0.d0; wk_ai(1, j) = 0.d0
                     end if
                  else if(mil0 == 1 .or. mil0 == 3) then
                     do i = 1, ibl2-ibl1+1
                        wk_ar(i,j) = fp*wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                        wk_ai(i,j) = -fp*wk_zfsin(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
                     end do
                     if(ibl1 == 1) then
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
      else ! if(kimg== 1 .and. k_symmetry(ik) == GAMMA) then
!!$         if(.not.allocated(wk_ar)) stop 'wk_ar is not allocated'
!!$         if(.not.allocated(wk_zfcos)) stop 'wk_zfcos is not allocated'
         j = 0
         !$acc loop independent
         do ia = 1, ia2-ia1+1
            it = ityp(ia+ia1-1)
            do lmt1 = 1, ilmt(it)
               !j = j+1
               j = lmtmap(lmt1, ia)
               lt = lmtt(lmt1,it)
               do i = 1, ibl2-ibl1+1
                  wk_ar(i,j) = wk_zfcos(i,ia)*snl_or_snld(i+ibl1-1,lt,iksnl)
               end do
            end do
         end do
      end if
      !$acc end kernels
      !$acc end data
      call tstatc0_end(id_sname)
    end subroutine G_dot_R_mult_snl_blk

    subroutine betar_dot_WFs_core_blk( ibsize, tran1,  psi_l )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer,intent(in)  :: ibsize
      logical, intent(in)  :: tran1
      real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) :: psi_l
      integer       :: icsize, M, LDA, LDB
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core_blk ',id_sname)

      M = msize
      icsize = ibl2-ibl1+1
      LDA = ibsize;   LDB = ibsize
      !$acc data present(psi_l, wk_ar, wk_ai, bp_tmp1, bp_tmp2)
      !$acc host_data use_device(psi_l, wk_ar, wk_ai, bp_tmp1, bp_tmp2)
#ifdef _CUDA_
      if(kimg == 1) then
         if(tran1) then
            call cublasDgemm('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
         else
            call cublasDgemm('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',M,np_e,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
         end if
      else if(kimg == 2) then
         if(tran1) then
            call cublasDgemm('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',np_e,M,icsize,-1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
            call cublasDgemm('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,2),LDB,wk_ar,LDA,1.d0,bp_tmp2,LD11)
         else
            call cublasDgemm('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',M,np_e,icsize,-1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',M,np_e,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
            call cublasDgemm('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp2,LD11)
         endif
      end if
#else
      if(kimg == 1) then
         if(tran1) then
            call DGEMM__('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
         else
            call DGEMM__('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,np_e,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
         end if
      else if(kimg == 2) then
         if(tran1) then
            call DGEMM__('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',np_e,M,icsize,-1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,1),LDB,wk_ai,LDA,1.d0,bp_tmp2,LD11)
            call DGEMM__('T','N',np_e,M,icsize, 1.d0,psi_l(1,1,2),LDB,wk_ar,LDA,1.d0,bp_tmp2,LD11)
         else
            call DGEMM__('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,np_e,icsize,-1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,np_e,icsize, 1.d0,wk_ai,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp2,LD11)
            call DGEMM__('T','N',M,np_e,icsize, 1.d0,wk_ar,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp2,LD11)
         endif
      end if
#endif
      !$acc end host_data
      !$acc end data
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core_blk

    subroutine betar_dot_WFs_core2_blk( ibsize, tran1,  psi_l )
!
! Revised by T. Kokubo & D. Fukata (NEC), September 2009
!
      integer, intent(in)  :: ibsize
      logical, intent(in)  :: tran1
      real(kind=DP),intent(in), dimension(ibsize,np_e,kimg) :: psi_l  ! MPI

      integer       :: icsize, M, LDA, LDB
      integer :: id_sname = -1
      call tstatc0_begin('betar_dot_WFs_core2_blk ',id_sname)

      icsize = ibl2-ibl1+1
      M = msize
      LDA = ibsize; LDB = ibsize
      !$acc data present(psi_l, wk_ar, wk_ai, bp_tmp1)
      !$acc host_data use_device(psi_l, wk_ar, wk_ai, bp_tmp1)
#ifdef _CUDA_
      if(kimg == 1) then
         if(tran1) then
            call cublasDgemm('T','N',np_e,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
         else
            call cublasDgemm('T','N',M,np_e,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
         endif
      else if(kimg==2) then
         if(tran1) then
            call cublasDgemm('T','N',np_e,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',np_e,M,icsize,1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
         else
            call cublasDgemm('T','N',M,np_e,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call cublasDgemm('T','N',M,np_e,icsize,1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
         end if
      end if
#else
      if(kimg == 1) then
         if(tran1) then
            call DGEMM__('T','N',np_e,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
         else
            call DGEMM__('T','N',M,np_e,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
         endif
      else if(kimg==2) then
         if(tran1) then
            call DGEMM__('T','N',np_e,M,icsize,1.d0,psi_l(1,1,1),LDB,wk_ar,LDA,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',np_e,M,icsize,1.d0,psi_l(1,1,2),LDB,wk_ai,LDA,1.d0,bp_tmp1,LD11)
         else
            call DGEMM__('T','N',M,np_e,icsize,1.d0,wk_ar,LDA,psi_l(1,1,1),LDB,1.d0,bp_tmp1,LD11)
            call DGEMM__('T','N',M,np_e,icsize,1.d0,wk_ai,LDA,psi_l(1,1,2),LDB,1.d0,bp_tmp1,LD11)
         end if
      end if
#endif
      !$acc end host_data
      !$acc end data
      call tstatc0_end(id_sname)
    end subroutine betar_dot_WFs_core2_blk
 end subroutine m_ES_betar_dot_WFs_4_lmta_k_blk

  subroutine m_ES_betar_dot_WFs_exp(ik,k1,k2,psi_l,i_exp_snl,bpr_l,bpi_l,map)
    integer, intent(in) :: ik,k1,k2
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg):: psi_l       ! MPI
    real(kind=DP),intent(in),dimension(kg1,nlmta,ista_snl:iend_snl,2) :: i_exp_snl
    real(kind=DP),intent(inout),dimension(np_e,nlmta,k1:k2) ::  bpr_l,bpi_l !MPI
    logical, intent(in), dimension(np_e), optional :: map
    real(kind=DP), allocatable, dimension(:,:,:) :: zajtmp
    real(kind=DP), allocatable, dimension(:,:) :: resr,resi
    integer :: ig,ib,it,iksnl,lmt1,ia,lmtt1,lmta1,icount,icountb,ik_for_pointing_eko
    integer :: ndata
    logical :: umap
    real(kind=DP) :: e
    integer :: id_sname = -1
    call tstatc0_begin('m_ES_betar_dot_WFs_exp ',id_sname,1)
    if ( noncol ) then
       iksnl = (ik-1)/ndim_spinor + 1   
    else
       iksnl = (ik-1)/nspin + 1         
    endif
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
      allocate(zajtmp(iba(ik),ndata,kimg))
      icount = 0
      do ib=1,np_e
        if(.not.map(ib)) then
          icount = icount+1
          zajtmp(1:iba(ik),icount,1:kimg) = psi_l(1:iba(ik),ib,ik,1:kimg)
        endif
      enddo
      allocate(resr(nlmta,ndata))
      if(k_symmetry(ik)/=GAMMA) allocate(resi(nlmta,ndata))
    else
      allocate(zajtmp(iba(ik),np_e,kimg))
      zajtmp(1:iba(ik),1:np_e,1:kimg) = psi_l(1:iba(ik),1:np_e,ik,1:kimg)
      allocate(resr(nlmta,np_e))
      if(k_symmetry(ik)/=GAMMA) allocate(resi(nlmta,np_e))
    endif
    if(.not.umap) then
      if(kimg==1) then
        call DGEMM__('T','N',nlmta,np_e,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,np_e,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resi,nlmta)
      else
        call DGEMM__('T','N',nlmta,np_e,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,np_e,iba(ik),-1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,2),iba(ik), 1.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) then
          call DGEMM__('T','N',nlmta,np_e,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resi,nlmta)
          call DGEMM__('T','N',nlmta,np_e,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,2),iba(ik), 1.d0,resi,nlmta)
        endif
      endif
      icount=0
      do ia=1,natm
        it = ityp(ia)
        do lmt1 = 1, ilmt(it)
          lmta1 = lmta(lmt1,ia)
          icount = icount+1
          do ib=1,np_e
            bpr_l(ib,lmta1,ik) = resr(icount,ib)
            if(k_symmetry(ik) /= GAMMA) bpi_l(ib,lmta1,ik) = resi(icount,ib)
          enddo
        enddo
      enddo
    else
      if(kimg==1) then
        call DGEMM__('T','N',nlmta,ndata,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,ndata,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resi,nlmta)
      else
        call DGEMM__('T','N',nlmta,ndata,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resr,nlmta)
        call DGEMM__('T','N',nlmta,ndata,iba(ik),-1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,2),iba(ik), 1.d0,resr,nlmta)
        if(k_symmetry(ik) /= GAMMA) then
          call DGEMM__('T','N',nlmta,ndata,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,2),kg1,zajtmp(:,:,1),iba(ik), 0.d0,resi,nlmta)
          call DGEMM__('T','N',nlmta,ndata,iba(ik), 1.d0, i_exp_snl(:,:,iksnl,1),kg1,zajtmp(:,:,2),iba(ik), 1.d0,resi,nlmta)
        endif
      endif
      icount=0
      do ia=1,natm
        it = ityp(ia)
        do lmt1 = 1, ilmt(it)
          lmta1 = lmta(lmt1,ia)
          icount = icount+1
          icountb = 0
          do ib=1,np_e
            if(.not.map(ib)) then
              icountb = icountb+1
              bpr_l(ib,lmta1,ik) = resr(icount,icountb)
              if(k_symmetry(ik) /= GAMMA) bpi_l(ib,lmta1,ik) = resi(icount,icountb)
            endif
          enddo
        enddo
      enddo
    endif
    deallocate(zajtmp)
    deallocate(resr)
    if(k_symmetry(ik)/=GAMMA) deallocate(resi)
    call tstatc0_end(id_sname)
  end subroutine m_ES_betar_dot_WFs_exp

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
    call m_ES_alloc_arai(0)
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
       if ( noncol ) then
          if ( kv3/ndim_spinor == 1) then
            call G_dot_R_map(ia,1)
            mapmode = MAPPED
         else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
         end if
       else
          if ( kv3/nspin == 1) then
            call G_dot_R_map(ia,1)
            mapmode = MAPPED
         else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
         end if
       endif
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
    call m_ES_dealloc_arai()
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
    real(kind=DP),intent(in), dimension(kg1,np_e,k1:k2,kimg) :: psi_l        ! MPI
    real(kind=DP),intent(out),dimension(np_e,nlmta_add,k1:k2) ::    bpr_l,bpi_l  ! MPI
    real(kind=DP),intent(in), dimension(kg1,nlmtt_add,ista_snl:iend_snl) :: snl_add

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
      integer :: i, i1

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA)) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
               ar(i) = zfcos(i)*snl_add(i,lmtt1,iksnl)
               ai(i) = zfsin(i)*snl_add(i,lmtt1,iksnl)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
               i1 = nbase(i,ik)
               ar(i) = zfcos(i1)*snl_add(i,lmtt1,iksnl)
               ai(i) = zfsin(i1)*snl_add(i,lmtt1,iksnl)
            end do
         end if
      else ! if(kimg == 1 .and. k_symmetry(ik) == GAMMA) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
               ar(i) = zfcos(i)*snl_add(i,lmtt1,iksnl)
            end do
         else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
               ar(i) = zfcos(i1)*snl_add(i,lmtt1,iksnl)
            end do
         end if
      end if
    end subroutine G_dot_R_mult_snl

    subroutine betar_dot_WFs_core()
      integer       :: ib, i

      bpr_l(1:np_e,lmta1,ik) = 0.d0        ! MPI
      bpi_l(1:np_e,lmta1,ik) = 0.d0        ! MPI

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
         do ib = 1, np_e                   ! MPI
            do i = 1, iba(ik)
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)
               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)
            end do
         end do
      else if(kimg == 2) then
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
            do i = 1, iba(ik)
               bpr_l(ib,lmta1,ik) = bpr_l(ib,lmta1,ik) + ar(i)*psi_l(i,ib,ik,1)-ai(i)*psi_l(i,ib,ik,2)
               bpi_l(ib,lmta1,ik) = bpi_l(ib,lmta1,ik) + ai(i)*psi_l(i,ib,ik,1)+ar(i)*psi_l(i,ib,ik,2)
            end do
         end do
      end if
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

  subroutine m_ES_phir_dot_WFs(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, ikphig, mapmode
    integer     :: id_sname = -1
    call tstatc0_begin('m_ES_phir_dot_WFs ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_phir_dot_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call m_ES_alloc_arai(0)
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
       if ( noncol ) then
         if( kv3/ndim_spinor == 1) then
            call G_dot_R_map(ia,1)
            mapmode = MAPPED
         else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
         end if
       else
         if(kv3/nspin == 1) then
            call G_dot_R_map(ia,1)
            mapmode = MAPPED
         else
            call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
            mapmode = NOTMAPPED
         end if
       endif
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

          call m_ES_phir_dot_WFs_4_lmta_k( ista_k, iend_k, ik, zaj_l, ia, ikphig, &
               &                           phig, nopr, compr_l, compi_l, mapmode )
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_arai()
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
  end subroutine m_ES_phir_dot_WFs

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

    subroutine G_dot_R_mult_phig_unfolding()
      integer :: i, i1

      ar = 0.0d0;  ai = 0.0d0

      if(kimg == 2 .or. (kimg ==1 .and. k_symmetry(ik) /= GAMMA) ) then
         if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
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
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
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
            do i = 1, iba(ik)
               i1 = nbase(i,ik)
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
            do i = 1, iba(ik)
               i1    = nbase(i,ik)
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

  subroutine m_ES_PAO_WFs(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, iksnl, mapmode
    integer     :: id_sname = -1

    call tstatc0_begin('m_ES_PAO_WFs ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_PAO_WFs --")')
    call m_ES_alloc_zfsincos(0)
    call alloc_zfsincos_mpi()
    do ia = 1, natm
       if(kv3/nspin == 1) then
          call G_dot_R_map(ia,1)
          mapmode = MAPPED
       else
          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = NOTMAPPED
       end if
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle         ! MPI
          iksnl = (ik-1)/nspin + 1
          call m_ES_PAO_WFs_4_lmta_k(ista_k,iend_k,ik,ia,iksnl,paog,zaj_l,mapmode)
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_zfsincos()
    call tstatc0_end(id_sname)
  end subroutine m_ES_PAO_WFs

! ===================================== added by K. Tagami ============== 11.0
  subroutine m_ES_PAO_WFs_noncl(nfout)
    integer, intent(in)    :: nfout

    integer ia, ik, iksnl, mapmode
    integer :: ikskip

    integer     :: id_sname = -1

    call tstatc0_begin('m_ES_PAO_WFs_noncl ',id_sname,level=1)

    if(ipribetar >= 2) write(nfout,'(" -- m_ES_PAO_WFs_noncl --")')
    call m_ES_alloc_zfsincos(0)
    call alloc_zfsincos_mpi()

    call m_ES_alloc_spinor_eigenwfn_0
    call m_ES_set_spinor_eigenwfn_0

    do ia = 1, natm
       if(kv3/ndim_spinor == 1) then
          call G_dot_R_map(ia,1)
          mapmode = MAPPED
       else
          call G_dot_R_mpi(ia) ! pos(ia,1:3), ngabc -> zfcos,zfsin
          mapmode = NOTMAPPED
       end if

       ikskip = ndim_spinor
       do ik = 1, kv3, ikskip
          if(map_k(ik) /= myrank_k) cycle         ! MPI
          iksnl = (ik-1)/ndim_spinor + 1

          if ( SpinOrbit_mode == Neglected ) then
             call m_ES_PAO_WFs_4_lmta_k_nonclA( ista_k, iend_k, ik, ia, iksnl,&
                  &                              paog, zaj_l, mapmode )
          else if ( SpinOrbit_mode == ByProjector .or. &
                    SpinOrbit_mode == ZeffApprox ) then
             call m_ES_PAO_WFs_4_lmta_k_nonclB( ista_k, iend_k, ik, ia, iksnl,&
                  &                              paog, zaj_l, mapmode )
          else if ( SpinOrbit_mode == ByPawPot .or. &
               &    SpinOrbit_mode == ReadFromPP ) then
             call m_ES_PAO_WFs_4_lmta_k_nonclB( ista_k, iend_k, ik, ia, iksnl,&
                  &                              paog, zaj_l, mapmode )
          else if ( SpinOrbit_mode == BuiltIn ) then
!             call m_ES_PAO_WFs_4_lmta_k_nonclC( ista_k, iend_k, ik, ia, iksnl,&
!                  &                             paog_soc, zaj_l, mapmode )
          endif
       end do
    end do
    call dealloc_zfsincos_mpi()
    call m_ES_dealloc_zfsincos()
    call m_ES_dealloc_spinor_eigenwfn_0

    call tstatc0_end(id_sname)

  end subroutine m_ES_PAO_WFs_noncl
! =========================================================================== 11.0

  subroutine m_ES_PAO_WFs_4_lmta_k(k1,k2,ik,ia,iksnl,paog,psi_l,mapmode)
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(kg1,nlmtt_pao,ista_snl:iend_snl) :: paog ! MPI
    real(kind=DP),intent(out), dimension(kg1,ista_e:iend_e,k1:k2,kimg) :: psi_l        ! MPI

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
               do i = 1, iba(ik)
                  psi_l(i,ib,ik,1) = zfcos(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, iba(ik)
                  i1    = nbase(i,ik)
                  psi_l(i,ib,ik,1) = zfcos(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, iba(ik)
                  psi_l(i,ib,ik,1) = zfsin(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, iba(ik)
                  i1    = nbase(i,ik)
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
               do i = 1, iba(ik)
                  psi_l(i,ib,ik,1) = zfcos(i)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) = zfsin(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, iba(ik)
                  i1    = nbase(i,ik)
                  psi_l(i,ib,ik,1) = zfcos(i1)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) = zfsin(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               do i = 1, iba(ik)
                  psi_l(i,ib,ik,1) = -zfsin(i)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) =  zfcos(i)*paog(i,lmtt1,iksnl)
               end do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               do i = 1, iba(ik)
                  i1    = nbase(i,ik)
                  psi_l(i,ib,ik,1) = -zfsin(i1)*paog(i,lmtt1,iksnl)
                  psi_l(i,ib,ik,2) =  zfcos(i1)*paog(i,lmtt1,iksnl)
               end do
            end if
         end if
      end if
    end subroutine G_dot_R_mult_paog

  end subroutine m_ES_PAO_WFs_4_lmta_k

! ============================== added by K. Tagami ==================== 11.0
  subroutine m_ES_PAO_WFs_4_lmta_k_nonclA(k1,k2,ik,ia,iksnl,paog,psi_l,mapmode)
! 
!                                       in the case of no spin-orbit 
! 
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(kg1,nlmtt_pao,ista_snl:iend_snl) :: paog
                                                                      ! MPI
    real(kind=DP),intent(out), dimension(kg1,ista_e:iend_e,k1:k2,kimg) :: psi_l
                                                                      ! MPI

    integer                    ::it, lmt1, lmtt1, lmta1, il1, ib, ib2
    integer :: is, ib_kt, ib2_kt

! ---------------------------------
    it    = ityp(ia)
! 
    do lmt1 = 1, ilmt_pao(it)
       lmtt1 = lmtt_pao(lmt1,it)
       lmta1 = lmta_pao(lmt1,ia)
       il1   = ltp_pao(lmt1,it)
       ib    = ibpao(lmta1)
       ib2   = 0
       ib2_kt   = 0
! ------------
       Do is=1, ndim_spinor
          ib_kt = ndim_spinor *( ib -1 )+is
          if (ib < 0) then
             ib = abs(ib)
             ib_kt  = ndim_spinor *(ib-1) +is
             ib2_kt = ndim_spinor * ib + is
          end if

!          write(*,*) 'is, ib_kt = ', is, ib_kt, ib2_kt

          if ( ib_kt >= ista_e .and. ib_kt <= iend_e) then
             call G_dot_R_mult_paog_noncl( ib_kt,.true. )             ! exp(iGR)
          endif
          if ( ib2_kt >= ista_e .and. ib2_kt <= iend_e) then
             call G_dot_R_mult_paog_noncl( ib2_kt,.false.)             ! exp(iGR)
          end if
       End do
    end do

  contains

    subroutine G_dot_R_mult_paog_noncl(ib,fcos)
      integer, intent(in) :: ib
      logical, intent(in) :: fcos
      integer :: i, i1, mil
!
      integer :: is1
!
      real(kind=DP) :: c1, s1, ctmp1, stmp1
      complex(kind=CMPLDP) :: zfac( ndim_spinor )
!
      if ( mod(ib, ndim_spinor ) == 1 ) then
         zfac(:) = Spinor_EigenWfn0_atomtyp( it,:,1 )
      else
         zfac(:) = Spinor_EigenWfn0_atomtyp( it,:,2 )
      endif
!
      if (kimg==1) then
         call phase_error_with_msg(nfout, 'Not supported : kimg = 1 for non-collinear system',__LINE__,__FILE__)
      else
         mil = mod(il1,4)
         if(mil == 1 .or. mil == 3) then
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  c1 = real( zfac(is1) );  s1 = aimag( zfac(is1) )
                  do i = 1, iba(ik)
                     ctmp1 = c1 *zfcos(i) -s1 *zfsin(i)
                     stmp1 = c1 *zfsin(i) +s1 *zfcos(i)
                     psi_l(i,ib,ik+is1-1,1) = ctmp1 *paog(i,lmtt1,iksnl)
                     psi_l(i,ib,ik+is1-1,2) = stmp1 *paog(i,lmtt1,iksnl)
                  end do
               End do

!               write(*,*) 'A lmtt1 ', lmtt1, paog(1,lmtt1,iksnl)
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  c1 = real( zfac(is1) );  s1 = aimag( zfac(is1) )
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ctmp1 = c1 *zfcos(i1) -s1 *zfsin(i1)
                     stmp1 = c1 *zfsin(i1) +s1 *zfcos(i1)

                     psi_l(i,ib,ik+is1-1,1) = ctmp1 *paog(i,lmtt1,iksnl)
                     psi_l(i,ib,ik+is1-1,2) = stmp1 *paog(i,lmtt1,iksnl)
                  end do
               End do

!               write(*,*) 'B lmtt1 ', lmtt1, paog(1,lmtt1,iksnl)
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  c1 = real( zfac(is1) );  s1 = aimag( zfac(is1) )
                  do i = 1, iba(ik)
                     ctmp1 = -c1 *zfsin(i) -s1 *zfcos(i)
                     stmp1 =  c1 *zfcos(i) -s1 *zfsin(i)

                     psi_l(i,ib,ik+is1-1,1) = ctmp1 *paog(i,lmtt1,iksnl)
                     psi_l(i,ib,ik+is1-1,2) = stmp1 *paog(i,lmtt1,iksnl)

                  end do
               End do
!               write(*,*) 'C lmtt1 ', lmtt1, paog(1,lmtt1,iksnl)
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  c1 = real( zfac(is1) );  s1 = aimag( zfac(is1) )
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ctmp1 = -c1 *zfsin(i1) -s1 *zfcos(i1)
                     stmp1 =  c1 *zfcos(i1) -s1 *zfsin(i1)

                     psi_l(i,ib,ik+is1-1,1) = ctmp1 *paog(i,lmtt1,iksnl)
                     psi_l(i,ib,ik+is1-1,2) = stmp1 *paog(i,lmtt1,iksnl)
                  end do
               End do
!               write(*,*) 'D lmtt1 ', lmtt1, paog(1,lmtt1,iksnl)
            end if
         end if
      end if
    end subroutine G_dot_R_mult_paog_noncl

  end subroutine m_ES_PAO_WFs_4_lmta_k_nonclA

  subroutine m_ES_PAO_WFs_4_lmta_k_nonclB(k1,k2,ik,ia,iksnl,paog,psi_l,mapmode)
! 
!                                       in the case of spin-orbit 
! 
    integer, intent(in)    :: k1,k2, ik, ia, iksnl, mapmode
    real(kind=DP),intent(in), dimension(kg1,nlmtt_pao,ista_snl:iend_snl) :: paog 
    real(kind=DP),intent(out), dimension(kg1,ista_e:iend_e,k1:k2,kimg) :: psi_l  

    integer :: it, lmt1, lmtt1, lmta1, il1, ib, ib2
    integer :: is1, is2, ib_kt, ib2_kt
    integer :: lmt2, lmtt2, il2, im1, im2
    integer :: itau1, itau2
    integer :: nsize
!
    complex(kind=CMPLDP), allocatable :: paog_tmp( :,: )
    complex(kind=CMPLDP) :: z1

! ---------------------------------
    allocate( paog_tmp( kg1, ndim_spinor ) ); paog_tmp = 0.d0
! ----------
!    Do ib=1, 1*2
!       Do ib2=1, 1*2
!          write(*,*) ib, ib2, EigenWfns_MatLS_L0( ib, ib2 )
!       End do
!    End do
!    Do ib2=1, 3*2
!       write(*,*) 'Eigen Val ib2 ', ib2, EigenVals_MatLS_L1( ib2 )
!       Do ib=1, 3*2
!          write(*,*) ib2, ib, EigenWfns_MatLS_L1( ib, ib2 )
!       End do
!    End do
!    stop

! --------
    it    = ityp(ia)
! 
    do lmt1 = 1, ilmt_pao(it)
       lmtt1 = lmtt_pao(lmt1,it)
       lmta1 = lmta_pao(lmt1,ia)

       il1   = ltp_pao(lmt1,it)
       im1   = mtp_pao(lmt1,it)
       itau1  = taup_pao(lmt1,it)

       ib    = ibpao(lmta1)
       ib2   = 0
       ib2_kt   = 0

       Do is1=1, ndim_spinor
          paog_tmp = 0.0d0

          Do lmt2=1, ilmt_pao(it)

             lmtt2 = lmtt_pao( lmt2,it )
             il2 = ltp_pao( lmt2,it )
             im2 = mtp_pao( lmt2,it )
             itau2 = taup_pao( lmt2,it )
          
             if ( itau1 /= itau2 ) cycle
             if ( il1 /= il2 ) cycle

             nsize = 2* ( il1 -1 ) +1

             Do is2=1, ndim_spinor
                if ( il1 == 1 ) then                 ! s-orbital
                   z1 = EigenWfns_MatLS_L0( im2 +nsize*( is2 -1 ), &
                        &                   im1 +nsize*( is1 -1 ) )
                   paog_tmp(:,is2) = paog_tmp(:,is2) + z1 *paog( :,lmtt2,iksnl )

                else if ( il1 == 2 ) then               ! p-orbital
                   z1 = EigenWfns_MatLS_L1( im2 +nsize*( is2 -1 ), &
                        &                   im1 +nsize*( is1 -1 ) )
                   paog_tmp(:,is2) = paog_tmp(:,is2) + z1 *paog( :,lmtt2,iksnl )

                else if ( il1 == 3 ) then               ! d-orbital
                   z1 = EigenWfns_MatLS_L2( im2 +nsize*( is2 -1 ), &
                        &                   im1 +nsize*( is1 -1 ) )
                   paog_tmp(:,is2) = paog_tmp(:,is2) + z1 *paog( :,lmtt2,iksnl )

                else if ( il1 == 4 ) then               ! f-orbital
                   z1 = EigenWfns_MatLS_L3( im2 +nsize*( is2 -1 ), &
                        &                   im1 +nsize*( is1 -1 ) )
                   paog_tmp(:,is2) = paog_tmp(:,is2) + z1 *paog( :,lmtt2,iksnl )
                endif
             End Do
          End do

! --
          ib_kt = ndim_spinor *( ib -1 )+is1

          if (ib < 0) then
             ib = abs(ib)
             ib_kt  = ndim_spinor *(ib-1) +is1
             ib2_kt = ndim_spinor * ib + is1
          end if

!!!!!!!          write(*,*) 'is, ib_kt = ', is, ib_kt, ib2_kt

          if ( ib_kt >= ista_e .and. ib_kt <= iend_e) then
             call G_dot_R_mult_paog_noncl( ib_kt,.true. )             ! exp(iGR)
!!             call G_dot_R_mult_paog_noncl0( ib_kt,.true. )             ! exp(iGR)
          endif
          if ( ib2_kt >= ista_e .and. ib2_kt <= iend_e) then
             call G_dot_R_mult_paog_noncl( ib2_kt,.false.)             ! exp(iGR)
!!             call G_dot_R_mult_paog_noncl0( ib2_kt,.false.)             ! exp(iGR)
          end if
       End do
    end do

    deallocate( paog_tmp )

  contains

    subroutine G_dot_R_mult_paog_noncl(ib,fcos)
      integer, intent(in) :: ib
      logical, intent(in) :: fcos
      integer :: i, i1, mil
!
      integer :: is1
!
      complex(kind=CMPLDP) :: zfac( ndim_spinor, ndim_spinor ), ztmp1
!
      if ( mod(ib, ndim_spinor ) == 1 ) then
         zfac(:,:) = Spinor_EigenWfn0_atomtyp( it,:,: )
      else
         zfac(:,:) = Spinor_EigenWfn0_atomtyp( it,:,: )
      endif
!
!      write(*,*) 'zfac = ', zfac(1,1)
!      write(*,*) 'zfac = ', zfac(1,2)
!      write(*,*) 'zfac = ', zfac(2,1)
!      write(*,*) 'zfac = ', zfac(2,2)
      !      stop
!
      if (kimg==1) then
         call phase_error_with_msg(nfout, 'Not supported : kimg = 1 for non-collinear system',__LINE__,__FILE__)
      else
         mil = mod(il1,4)
         if(mil == 1 .or. mil == 3) then
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     ztmp1 = dcmplx( zfcos(i),zfsin(i) ) &
                          &  *( zfac(is1,1) *paog_tmp(i,1) &
                          &   + zfac(is1,2) *paog_tmp(i,2) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ztmp1 = dcmplx( zfcos(i1),zfsin(i1) ) &
                          &  *( zfac(is1,1) *paog_tmp(i,1) &
                          &   + zfac(is1,2) *paog_tmp(i,2) )
                     psi_l(i,ib,ik+is1-1,1) = real(  ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     ztmp1 = dcmplx( -zfsin(i),zfcos(i) ) &
                          &  *( zfac(is1,1) *paog_tmp(i,1) &
                          &   + zfac(is1,2) *paog_tmp(i,2) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ztmp1 = dcmplx( -zfsin(i1),zfcos(i1) ) &
                          &  *( zfac(is1,1) *paog_tmp(i,1) &
                          &   + zfac(is1,2) *paog_tmp(i,2) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            end if
         end if
      end if
    end subroutine G_dot_R_mult_paog_noncl

    subroutine G_dot_R_mult_paog_noncl0(ib,fcos)
      integer, intent(in) :: ib
      logical, intent(in) :: fcos
      integer :: i, i1, mil
!
      integer :: is1
!
      complex(kind=CMPLDP) :: zfac( ndim_spinor ), ztmp1
!
      if ( mod(ib, ndim_spinor ) == 1 ) then
         zfac(:) = Spinor_EigenWfn0_atomtyp( it,:,1 )
      else
         zfac(:) = Spinor_EigenWfn0_atomtyp( it,:,2 )
      endif
!
      write(*,*) 'zfac = ', zfac(1)
      write(*,*) 'zfac = ', zfac(2)
!      stop
!
      if (kimg==1) then
         call phase_error_with_msg(nfout, 'Not supported : kimg = 1 for non-collinear system',__LINE__,__FILE__)
      else
         mil = mod(il1,4)
         if(mil == 1 .or. mil == 3) then
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     ztmp1 = dcmplx( zfcos(i),zfsin(i) ) &
                          &  *( zfac(is1) *paog_tmp(i,is1) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ztmp1 = dcmplx( zfcos(i1),zfsin(i1) ) &
                          &  *( zfac(is1) *paog_tmp(i,is1) )
                     psi_l(i,ib,ik+is1-1,1) = real(  ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            end if
         else
            if(mapmode == MAPPED) then
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     ztmp1 = dcmplx( -zfsin(i),zfcos(i) ) &
                          &  *( zfac(is1) *paog_tmp(i,is1) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            else
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(i1)
#endif
               Do is1=1, ndim_spinor
                  do i = 1, iba(ik)
                     i1    = nbase(i,ik)
                     ztmp1 = dcmplx( -zfsin(i1),zfcos(i1) ) &
                          &  *( zfac(is1) *paog_tmp(i,is1) )
                     psi_l(i,ib,ik+is1-1,1) = real( ztmp1 )
                     psi_l(i,ib,ik+is1-1,2) = aimag( ztmp1 )
                  end do
               End do
            end if
         end if
      end if
    end subroutine G_dot_R_mult_paog_noncl0

  end subroutine m_ES_PAO_WFs_4_lmta_k_nonclB
! ============================================================================ 11.0

  subroutine betar_dot_WFs_in_rspace(nfout)
    integer, intent(in) :: nfout
    integer :: ik
    integer :: id_sname = -1
    call tstatc0_begin('betar_dot_WFs_in_rspace ',id_sname,level=1)
    do ik=1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle
       call betar_dot_Psi_4_each_k_in_rs(nfout,ista_k,iend_k,ik,zaj_l,fsr_l,fsi_l)
    enddo

    if(ipribetar >= 2) then
       write(nfout,'(" --- fsr_l, fsi_l ---")')
       do ik = ista_k, iend_k                              ! MPI
          call wd_fsr_fsi(ista_k,iend_k,ik,fsr_l,fsi_l)    ! MPI
       end do
    end if
    call tstatc0_end(id_sname)
  end subroutine betar_dot_WFs_in_rspace

  subroutine betar_dot_Psi_4_each_k_in_rs(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l)
    integer,intent(in) :: nfout 
    integer,intent(in) :: k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(out),dimension(np_e,nlmta,k1:k2) ::    bpr_l,bpi_l  ! MPI
    call betar_dot_Psi_4_each_k_in_rs0(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l,nmesh_rs_max,snl_rs)
  end subroutine betar_dot_Psi_4_each_k_in_rs

  subroutine betar_dot_Psi_4_each_k_in_rs0(nfout,k1,k2,ik,psi_l,bpr_l,bpi_l,nsnl,snlsnl)
    integer,intent(in) :: nfout 
    integer,intent(in) :: k1,k2,ik
    real(kind=DP),intent(in),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(out),dimension(np_e,nlmta,k1:k2) ::    bpr_l,bpi_l  ! MPI
    integer, intent(in) :: nsnl
    real(kind=DP), intent(in), dimension(nsnl,nlmta) :: snlsnl
    integer :: ia,ib,ibsize,ib1,ib2
    real(kind=DP),allocatable,dimension(:,:) :: cos_kr,sin_kr
#ifdef MULT_PHASE_RSPACE
    real(kind=DP),allocatable,dimension(:) :: cos_a,sin_a
#endif
    real(kind=DP),allocatable,dimension(:,:) :: bffb
    real(kind=DP),allocatable,dimension(:,:) :: rrrb,iiib
    integer :: id_sname = -1
    call tstatc0_begin('betar_dot_Psi_4_each_k_in_rs ',id_sname,level=1)

    ibsize = nblocksize_rspace_betar
    if(ibsize<1) ibsize = 1
    if(ibsize>np_e) ibsize = np_e

    call alloc_work_arrays()
    call k_dot_r(ik,cos_kr,sin_kr)
#ifdef MULT_PHASE_RSPACE
    call k_dot_pos(ik,cos_a,sin_a)
#endif

    do ib=1,np_e,ibsize
       ib2 = min(ibsize,np_e+1-ib)
       do ib1=1,ib2
          call m_ES_WF_in_Rspace(k1,k2,ik,ista_e+ib+ib1-2,psi_l,bffb(1:nfft,ib1))
       enddo
       do ia=1,natm
          call betar_dot_Psi_atm_band(bffb,ia)
       enddo
    enddo

    call dealloc_work_arrays()
    call tstatc0_end(id_sname)

    contains

    subroutine alloc_work_arrays()
       allocate(bffb(nfft,ibsize));bffb=0.d0
       allocate(cos_kr(nmesh_rs_max,natm));cos_kr=0.d0
       allocate(sin_kr(nmesh_rs_max,natm));sin_kr=0.d0
#ifdef MULT_PHASE_RSPACE
       allocate(cos_a(natm));cos_a=0.d0
       allocate(sin_a(natm));sin_a=0.d0
#endif
       allocate(rrrb(nmesh_rs_max,ibsize));rrrb=0.d0
       allocate(iiib(nmesh_rs_max,ibsize));iiib=0.d0
    end subroutine alloc_work_arrays

    subroutine dealloc_work_arrays()
       deallocate(bffb)
       deallocate(cos_kr)
       deallocate(sin_kr)
#ifdef MULT_PHASE_RSPACE
       deallocate(cos_a)
       deallocate(sin_a)
#endif
       deallocate(rrrb)
       deallocate(iiib)
    end subroutine dealloc_work_arrays

    subroutine betar_dot_Psi_atm_band(bff,iatm)
       real(kind=DP),dimension(nfft,ibsize),intent(in) :: bff
       integer,intent(in) :: iatm
       integer :: lmt1,i,lmta1,imesh,nma,it,iband,ib1
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
#ifdef MULT_PHASE_RSPACE
       cosa = cos_a(iatm)
       sina = sin_a(iatm)
#endif
       if(kimg==1)then
          do ib1=1,ib2
          do i=1,nma
             imesh   = 2*meshxyz_rs(i,iatm)
#ifdef MULT_PHASE_RSPACE
             psr = bff(imesh-1,ib1)
             psi = bff(imesh,ib1)*meshxyz_rs_conjg(i,iatm)
             rr = psr*cosa-psi*sina
             ii = psi*cosa+psr*sina
#else
             rr = bff(imesh-1,ib1)
             ii = bff(imesh,ib1)*meshxyz_rs_conjg(i,iatm)
#endif
             rrrb(i,ib1) = rr*cos_kr(i,iatm)-ii*sin_kr(i,iatm) 
             iiib(i,ib1) = ii*cos_kr(i,iatm)+rr*sin_kr(i,iatm)
          enddo
          enddo
       else
          do ib1=1,ib2
          do i=1,nma
             imesh   = 2*meshxyz_rs(i,iatm)
#ifdef MULT_PHASE_RSPACE
             psr = bff(imesh-1,ib1)
             psi = bff(imesh,ib1)
             rr = psr*cosa-psi*sina
             ii = psi*cosa+psr*sina
#else
             rr = bff(imesh-1,ib1)
             ii = bff(imesh,ib1)
#endif
             rrrb(i,ib1) = rr*cos_kr(i,iatm)-ii*sin_kr(i,iatm) 
             iiib(i,ib1) = ii*cos_kr(i,iatm)+rr*sin_kr(i,iatm)
          enddo
          enddo
       endif
       
#ifdef NO_DGEMM_RSPACE
       do ib1=1,ib2
          iband = ib+ib1-1
          do lmt1=1,ilmt(it)
             lmta1 = lmta(lmt1,iatm)
#ifdef USE_DDOT
             bpr_l(iband,lmta1,ik) = fac*ddot(nma,snlsnl(1:nma,lmta1),1,rrrb(1:nma,ib1),1)
#else
             bpr_l(iband,lmta1,ik) = fac*dot_product(snlsnl(1:nma,lmta1),rrrb(1:nma,ib1))
#endif
          enddo
       enddo
#else
       call dgemm('T','N',ib2,ilmt(it),nma,fac,rrrb,nmesh_rs_max,snlsnl(1:nma,lmta(1,iatm)) &
                     ,nmesh_rs_max,0.d0,bpr_l(ib,lmta(1,iatm),ik),np_e)
#endif

       if(k_symmetry(ik)/=GAMMA)then
#ifdef NO_DGEMM_RSPACE
       do ib1=1,ib2
          iband = ib+ib1-1
          do lmt1=1,ilmt(it)
            lmta1 = lmta(lmt1,iatm)
#ifdef USE_DDOT
            bpi_l(iband,lmta1,ik) = fac*ddot(nma,snlsnl(1:nma,lmta1),1,iiib(1:nma,ib1),1)
#else
            bpi_l(iband,lmta1,ik) = fac*dot_product(snlsnl(1:nma,lmta1),iiib(1:nma,ib1))
#endif
          enddo
       enddo
#else
       call dgemm('T','N',ib2,ilmt(it),nma,fac,iiib,nmesh_rs_max,snlsnl(1:nma,lmta(1,iatm)) &
                     ,nmesh_rs_max,0.d0,bpi_l(ib,lmta(1,iatm),ik),np_e)
#endif
       endif
       call tstatc0_end(id_sname)
    end subroutine betar_dot_Psi_atm_band

  end subroutine betar_dot_Psi_4_each_k_in_rs0

  subroutine Vnonlocal_W_in_realspace(ik,ispin,switch_of_eko_part)
    integer, intent(in) :: ik,ispin
    integer, intent(in) :: switch_of_eko_part
    integer :: mdvdb, it, ia, lmt2, lmta2, il2, im2,nma
#ifndef RSPACE_DGEMM
    real(kind=DP),allocatable,dimension(:,:) :: vnlr,vnli
#else
    real(kind=DP),allocatable,dimension(:,:,:) :: vnlr,vnli
#endif
    real(kind=DP),allocatable,dimension(:) :: vnlr1d,vnli1d
#ifndef RSPACE_DGEMM
    real(kind=DP),allocatable,dimension(:) :: sc,qc
#else
    real(kind=DP),allocatable,dimension(:,:,:) :: sc,qc
#endif
    real(kind=DP),allocatable,dimension(:,:) :: sc_lmta,qc_lmta
    real(kind=DP),allocatable,dimension(:,:) :: cos_kr,sin_kr
#ifdef MULT_PHASE_RSPACE
    real(kind=DP),allocatable,dimension(:) :: cos_a,sin_a
#endif
    real(kind=DP),allocatable,dimension(:,:) :: bffb
    integer :: ib,ii,ib2
    integer :: ibl1,ibl2,nbsize,ibsize
    integer :: id_sname = -1
#ifdef RSPACE_DGEMM
    integer :: i, i1, lmt1, lmta1, il1, im1, nfft_, ialmt2, maxialmt2, iband
    real(kind=DP) :: fac, tmp
    real(kind=DP), pointer, dimension(:,:) :: snlt
    real(kind=DP), allocatable, dimension(:,:,:) :: atmp, btmp
#endif
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
    call k_dot_r(ik,cos_kr,sin_kr)
#ifdef MULT_PHASE_RSPACE
    call k_dot_pos(ik,cos_a,sin_a)
#endif

#ifndef RSPACE_DGEMM
    do ia=1,natm
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
    end do

    do ib=1,np_e,ibsize
       bffb=0.d0
       ib2 = min(ibsize,np_e+1-ib)
       do ia=1,natm
         it = ityp(ia)
         mdvdb = m_PP_include_vanderbilt_pot(it)
         if(switch_of_eko_part == OFF) mdvdb= SKIP
         vnlr = 0.d0
         vnli = 0.d0
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
       end do
!      vnl(RS) --> FT --> vnl(GS)
       call Vnonlocal_W_to_Gspace()
    enddo
#else
    fac = dsqrt(univol)/dble(fft_box_size_WF(1,1)*fft_box_size_WF(2,1)*fft_box_size_WF(3,1))
    nfft_ = nfft/2
    if(kimg == 1) then
       snlt => snl_rs_h
    else
       snlt => snl_rs
    end if

    ialmt2 = 0
    Loop_count: do ia = 1,natm
       it = ityp(ia)
       do lmt2 = 1, ilmt(it)
          ialmt2 = ialmt2 + 1
       end do
    end do Loop_count
    maxialmt2 = ialmt2

    allocate(sc(nfft_,maxialmt2,2)); sc = 0.0d0
    if(mdvdb == EXECUT) then
       allocate(qc(nfft_,maxialmt2,2)); qc = 0.0d0
    end if
    allocate(atmp(np_e,maxialmt2,4)); atmp = 0.0d0
    allocate(btmp(np_e,maxialmt2,4)); btmp = 0.0d0

    ialmt2 = 0
    do ia = 1, natm
       it = ityp(ia)
       mdvdb = m_PP_include_vanderbilt_pot(it)
       if(switch_of_eko_part == OFF) mdvdb= SKIP
       if(kimg == 1) then
          nma = nmesh_rs_h(ia)
       else
          nma = nmesh_rs(ia)
       end if
       do lmt2 = 1, ilmt(it)
          lmta2 = lmta(lmt2,ia)
          il2   = ltp(lmt2,it)
          im2   = mtp(lmt2,it)
          ialmt2 = ialmt2 + 1
          do lmt1 = 1,ilmt(it)
             lmta1 = lmta(lmt1,ia)
             il1   = ltp(lmt1,it)
             im1   = mtp(lmt1,it)
             if(il1 == il2 .and. im1 == im2) then
                if(ipaw(it) == 0) then
                   tmp = dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,ispin)
                else
                   tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                end if
             else
                if(ipaw(it) == 0) then
                   tmp = vlhxcQ(lmt1,lmt2,ia,ispin)
                else
                   tmp = dion_paw(lmt1,lmt2,ispin,ia) + vlhxcQ(lmt1,lmt2,ia,ispin)
                end if
             end if
             if(kimg == 1) then
                do i = 1, nma
                   i1 = meshxyz_rs_h(i,ia)
                   sc(i1,ialmt2,1) = sc(i1,ialmt2,1) + tmp*snlt(i,lmta1)*cos_kr(map_h(i,ia),ia)
                   sc(i1,ialmt2,2) = sc(i1,ialmt2,2) + tmp*snlt(i,lmta1)*sin_kr(map_h(i,ia),ia)
                end do
                if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
                   tmp = q(lmt1,lmt2,it)
                   do i = 1, nma
                      i1 = meshxyz_rs_h(i,ia)
                      qc(i1,ialmt2,1) = qc(i1,ialmt2,1) + tmp*snlt(i,lmta1)*cos_kr(map_h(i,ia),ia)
                      qc(i1,ialmt2,2) = qc(i1,ialmt2,2) + tmp*snlt(i,lmta1)*sin_kr(map_h(i,ia),ia)
                   end do
                end if
             else
                do i = 1, nma
                   i1 = meshxyz_rs(i,ia)
                   sc(i1,ialmt2,1) = sc(i1,ialmt2,1) + tmp*snlt(i,lmta1)*cos_kr(i,ia)
                   sc(i1,ialmt2,2) = sc(i1,ialmt2,2) + tmp*snlt(i,lmta1)*sin_kr(i,ia)
                end do
                if(mdvdb == EXECUT .and. il1 == il2 .and. im1 == im2) then
                   tmp = q(lmt1,lmt2,it)
                   do i = 1, nma
                      i1 = meshxyz_rs(i,ia)
                      qc(i1,ialmt2,1) = qc(i1,ialmt2,1) + tmp*snlt(i,lmta1)*cos_kr(i,ia)
                      qc(i1,ialmt2,2) = qc(i1,ialmt2,2) + tmp*snlt(i,lmta1)*sin_kr(i,ia)
                   end do
                end if
             end if
          end do
          if(mdvdb == SKIP)then
             do iband = 1, np_e
                atmp(iband,ialmt2,1) = fac*cos_a(ia)*fsr_l(iband,lmta2,ik)
                atmp(iband,ialmt2,2) = fac*sin_a(ia)*fsr_l(iband,lmta2,ik)
             end do
             if(k_symmetry(ik) /= GAMMA) then
                do iband = 1, np_e
                   atmp(iband,ialmt2,3) = fac*cos_a(ia)*fsi_l(iband,lmta2,ik)
                   atmp(iband,ialmt2,4) = fac*sin_a(ia)*fsi_l(iband,lmta2,ik)
                end do
             end if
          else
             do iband = 1, np_e
                atmp(iband,ialmt2,1) = fac*cos_a(ia)*fsr_l(iband,lmta2,ik)
                atmp(iband,ialmt2,2) = fac*sin_a(ia)*fsr_l(iband,lmta2,ik)
                btmp(iband,ialmt2,1) = eko_l(iband,ik)*atmp(iband,ialmt2,1)
                btmp(iband,ialmt2,2) = eko_l(iband,ik)*atmp(iband,ialmt2,2)
             end do
             if(k_symmetry(ik) /= GAMMA) then
                do iband = 1, np_e
                   atmp(iband,ialmt2,3) = fac*cos_a(ia)*fsi_l(iband,lmta2,ik)
                   atmp(iband,ialmt2,4) = fac*sin_a(ia)*fsi_l(iband,lmta2,ik)
                   btmp(iband,ialmt2,3) = eko_l(iband,ik)*atmp(iband,ialmt2,3)
                   btmp(iband,ialmt2,4) = eko_l(iband,ik)*atmp(iband,ialmt2,4)
                end do
             end if
          end if
       end do
    end do

    if(mdvdb == SKIP)then
       call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,1),nfft_,atmp(1,1,1),np_e,0.0d0,vnlr(1,1,1),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,2),nfft_,atmp(1,1,1),np_e,0.0d0,vnlr(1,1,2),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,1),nfft_,atmp(1,1,2),np_e,0.0d0,vnlr(1,1,3),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,2),nfft_,atmp(1,1,2),np_e,0.0d0,vnlr(1,1,4),nfft_)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,1),nfft_,atmp(1,1,4),np_e,0.0d0,vnli(1,1,1),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,2),nfft_,atmp(1,1,3),np_e,0.0d0,vnli(1,1,2),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,1),nfft_,atmp(1,1,3),np_e,0.0d0,vnli(1,1,3),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,1.0d0,sc(1,1,2),nfft_,atmp(1,1,4),np_e,0.0d0,vnli(1,1,4),nfft_)
       end if
    else
       call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,1),nfft_,atmp(1,1,1),np_e,0.0d0,vnlr(1,1,1),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,1),nfft_,btmp(1,1,1),np_e,1.0d0,vnlr(1,1,1),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,2),nfft_,atmp(1,1,1),np_e,0.0d0,vnlr(1,1,2),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,2),nfft_,btmp(1,1,1),np_e,1.0d0,vnlr(1,1,2),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,1),nfft_,atmp(1,1,2),np_e,0.0d0,vnlr(1,1,3),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,1),nfft_,btmp(1,1,2),np_e,1.0d0,vnlr(1,1,3),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,2),nfft_,atmp(1,1,2),np_e,0.0d0,vnlr(1,1,4),nfft_)
       call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,2),nfft_,btmp(1,1,2),np_e,1.0d0,vnlr(1,1,4),nfft_)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,1),nfft_,atmp(1,1,4),np_e,0.0d0,vnli(1,1,1),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,1),nfft_,btmp(1,1,4),np_e,1.0d0,vnli(1,1,1),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,2),nfft_,atmp(1,1,3),np_e,0.0d0,vnli(1,1,2),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,2),nfft_,btmp(1,1,3),np_e,1.0d0,vnli(1,1,2),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,1),nfft_,atmp(1,1,3),np_e,0.0d0,vnli(1,1,3),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,1),nfft_,btmp(1,1,3),np_e,1.0d0,vnli(1,1,3),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2, 1.0d0,sc(1,1,2),nfft_,atmp(1,1,4),np_e,0.0d0,vnli(1,1,4),nfft_)
          call dgemm('N','T',nfft_,np_e,maxialmt2,-1.0d0,qc(1,1,2),nfft_,btmp(1,1,4),np_e,1.0d0,vnli(1,1,4),nfft_)
       end if
    end if

    if(kimg == 1) then
       do iband = 1, np_e
          do i = 1, nfft_
             bffb(2*i-1,iband) = bffb(2*i-1,iband) + vnlr(i,iband,1) + vnli(i,iband,1)
             bffb(2*i-1,iband) = bffb(2*i-1,iband) - vnlr(i,iband,4) + vnli(i,iband,2)
             bffb(2*i,  iband) = bffb(2*i,  iband) + vnlr(i,iband,3) - vnli(i,iband,3)
             bffb(2*i,  iband) = bffb(2*i,  iband) + vnlr(i,iband,2) + vnli(i,iband,4)
          end do
       end do
    else
       do iband = 1, np_e
          do i = 1, nfft_
             bffb(2*i-1,iband) = bffb(2*i-1,iband) + vnlr(i,iband,1) + vnli(i,iband,1)
             bffb(2*i-1,iband) = bffb(2*i-1,iband) - vnlr(i,iband,4) + vnli(i,iband,2)
             bffb(2*i,  iband) = bffb(2*i,  iband) - vnlr(i,iband,3) + vnli(i,iband,3)
             bffb(2*i,  iband) = bffb(2*i,  iband) - vnlr(i,iband,2) - vnli(i,iband,4)
          end do
       end do
    end if

    deallocate(sc)
    if(mdvdb == EXECUT) deallocate(qc)
    deallocate(atmp)
    deallocate(btmp)

    call Vnonlocal_W_to_Gspace()
#endif

    call dealloc_arrays()

    call tstatc0_end(id_sname)

    contains

    subroutine alloc_arrays()
       integer :: nmm
       nmm = nmesh_rs_max
       if(kimg==1) nmm = nmesh_rs_max_h
#ifndef RSPACE_DGEMM
       allocate(vnlr(nmm,ibsize));vnlr=0.d0
       allocate(vnli(nmm,ibsize));vnli=0.d0
       !allocate(sc(nmm));sc=0.d0
       !if(mdvdb==EXECUT)then
       !   allocate(qc(nmm));qc=0.d0
       !endif
       allocate(sc_lmta(nmm,nlmta));sc_lmta=0.d0
       if(mdvdb==EXECUT)then
          allocate(qc_lmta(nmm,nlmta));qc_lmta=0.d0
       endif
#else
       allocate(vnlr(nfft/2,np_e,4));vnlr=0.d0
       allocate(vnli(nfft/2,np_e,4));vnli=0.d0
#endif
       allocate(cos_kr(nmesh_rs_max,natm));cos_kr=0.d0
       allocate(sin_kr(nmesh_rs_max,natm));sin_kr=0.d0
#ifdef MULT_PHASE_RSPACE
       allocate(cos_a(natm));cos_a=0.d0
       allocate(sin_a(natm));sin_a=0.d0
#endif
       allocate(bffb(nfft,ibsize));bffb=0.d0
    end subroutine alloc_arrays

    subroutine dealloc_arrays()
#ifndef RSPACE_DGEMM
       deallocate(vnlr)
       deallocate(vnli)
       !deallocate(sc)
       !if(allocated(qc)) deallocate(qc)
       deallocate(sc_lmta)
       if(allocated(qc_lmta)) deallocate(qc_lmta)
       deallocate(bffb)
#else
       deallocate(vnlr)
       deallocate(vnli)
       deallocate(bffb)
#endif
       deallocate(cos_kr)
       deallocate(sin_kr)
#ifdef MULT_PHASE_RSPACE
       deallocate(cos_a)
       deallocate(sin_a)
#endif
    end subroutine dealloc_arrays

#ifndef RSPACE_DGEMM
    subroutine part_sum_over_lmt1_rs()
       integer       :: lmta1,lmt1, il1, im1, i
       real(kind=DP) :: tmp,rr,ii
       real(kind=DP),pointer,dimension(:,:) :: snlt
       integer :: id_sname = -1
       call tstatc0_begin('part_sum_over_lmt1_rs ',id_sname,level=1)
       sc_lmta(:,lmta2) = 0.d0
       if(mdvdb == EXECUT) then
          qc_lmta(:,lmta2) = 0.d0
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
             fr = fsr_l(iband,lmta2,ik)
             e = eko_l(iband,ik)
             do i=1,nma
                vnlr(i,ib1) = vnlr(i,ib1)+fr*(sc_lmta(i,lmta2)-e*qc_lmta(i,lmta2))
             enddo
             if(k_symmetry(ik)/=GAMMA) then
                fi = fsi_l(iband,lmta2,ik)
                do i=1,nma
                   vnli(i,ib1) = vnli(i,ib1)+fi*(sc_lmta(i,lmta2)-e*qc_lmta(i,lmta2))
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
             fre(ib1,lmt2) = -fsr_l(iband,lmta2,ik)*eko_l(iband,ik)
             if(k_symmetry(ik) /= GAMMA) fie(ib1,lmt2) = -fsi_l(iband,lmta2,ik)*eko_l(iband,ik)
          enddo
       enddo
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsr_l(ib,lmta(1,ia),ik),np_e,0.d0,vnlr,nmm)
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,qc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fre,ib2,1.d0,vnlr,nmm)
       deallocate(fre)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsi_l(ib,lmta(1,ia),ik),np_e,0.d0,vnli,nmm)
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,qc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fie,ib2,1.d0,vnli,nmm)
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
             fr = fsr_l(iband,lmta2,ik)
             do i=1,nma
                vnlr(i,ib1) = vnlr(i,ib1)+fr*sc_lmta(i,lmta2)
             enddo
             if(k_symmetry(ik)/=GAMMA) then
                fi = fsi_l(iband,lmta2,ik)
                do i=1,nma
                   vnli(i,ib1) = vnli(i,ib1)+fi*sc_lmta(i,lmta2)
                enddo
             endif
          enddo
       enddo
#else
       nmm = nmesh_rs_max
       if(kimg==1) nmm = nmesh_rs_max_h
       call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsr_l(ib,lmta(1,ia),ik),np_e,0.d0,vnlr,nmm)
       if(k_symmetry(ik) /= GAMMA) then
          call dgemm('N','T',nma,ib2,ilmt(it),1.d0,sc_lmta(1,lmta(1,ia)),nmesh_rs_max &
                     , fsi_l(ib,lmta(1,ia),ik),np_e,0.d0,vnli,nmm)
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
                rr = cosa*vnlr(i,ib1)-sina*vnli(i,ib1)
                ii = sina*vnlr(i,ib1)+cosa*vnli(i,ib1)
#else
                rr = vnlr(i,ib1)
                ii = vnli(i,ib1)
#endif
                vnlr(i,ib1) = (rr*cos_kr(map_h(i,ia),ia)+ii*sin_kr(map_h(i,ia),ia))*fac
                vnli(i,ib1) = (ii*cos_kr(map_h(i,ia),ia)-rr*sin_kr(map_h(i,ia),ia))*fac
             enddo
          enddo
       else
          do ib1=1,ib2
             do i=1,nma
#ifdef MULT_PHASE_RSPACE
                rr = cosa*vnlr(i,ib1)-sina*vnli(i,ib1)
                ii = sina*vnlr(i,ib1)+cosa*vnli(i,ib1)
#else
                rr = vnlr(i,ib1)
                ii = vnli(i,ib1)
#endif
                vnlr(i,ib1) = (rr*cos_kr(i,ia)+ii*sin_kr(i,ia))*fac
                vnli(i,ib1) = (ii*cos_kr(i,ia)-rr*sin_kr(i,ia))*fac
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
                bffb(i1-1,ib1) = bffb(i1-1,ib1) + vnlr(i,ib1)
                bffb(i1,ib1)   = bffb(i1,ib1)   - vnli(i,ib1)
             enddo
          enddo
       else
          do ib1=1,ib2
             do i=1,nma
                i1 = 2*meshxyz_rs(i,ia)
                bffb(i1-1,ib1) = bffb(i1-1,ib1) + vnlr(i,ib1)
                bffb(i1,ib1)   = bffb(i1,ib1)   + vnli(i,ib1)
             enddo
          enddo
       endif
       call tstatc0_end(id_sname)
    end subroutine map_vnl_to_bff
#endif

    subroutine Vnonlocal_W_to_Gspace()
       integer :: i,i1,ri,iband,ib1
       integer :: ii
       integer :: id_sname = -1
       call tstatc0_begin('Vnonlocal_W_to_Gspace ',id_sname,level=1)
       do ib1=1,ib2
          iband = ib1+ib-1
          vnlph_l(:,iband,:) = 0.d0
          call m_FFT_WF(ELECTRON,nfout,bffb(1:nfft,ib1),DIRECT,ON)
          do ri=1,kimg
             do i=1,iba(ik)
                i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
                vnlph_l(i,iband,ri) = bffb(i1,ib1)
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
    do ia=1,natm
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
