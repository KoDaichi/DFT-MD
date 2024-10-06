!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_Electronic_Structure
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
#ifdef SAVE_FFT_TIMES
#define DEBUG_SAVE_FFT_TIMES
#endif
#ifdef __TIMER_SUB__
#   define __TIMER_SUB_START(a)  call timer_sta(a)
#   define __TIMER_SUB_STOP(a)   call timer_end(a)
#else
#   define __TIMER_SUB_START(a)
#   define __TIMER_SUB_STOP(a)
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

module m_Electronic_Structure
! $Id: m_Electronic_Structure.F90 633 2020-12-01 05:11:03Z jkoga $
#ifdef MPI_FFTW
  use, intrinsic :: iso_c_binding
#endif
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, nkgroup, iteration
!!!!!BRANCH_P ORG_Parallel
! === DEBUG for icond=2,3 & one_by_one by tkato 2014/01/24 =====================
  use m_IterationNumbers,   only : first_kpoint_in_this_job
! ==============================================================================
!!!!!BRANCH_P_END ORG_Parallel
  use m_NonLocal_Potential, only : snl,norm_phig
  use m_PlaneWaveBasisSet,  only : igf,kg1,kg,kgp,nbase,nbmx,iba &
       &                         , nbase_gamma,kg1_prev, iba_prev &
       &                         , nbase_prev, nbase_gamma_prev, igf_prev &
       &                         , GVec_on_refcell
  use m_PlaneWaveBasisSet,   only : m_pwBS_kinetic_energies
#ifdef __EDA__
  use m_PlaneWaveBasisSet,   only : ngabc, gr_l, m_pwBS_kinetic_energies
#endif
  use m_PseudoPotential,    only : ival,ilmt,nlmt,nlmtt,nlmta,lmta,lmtt,ltp,mtp,q,dion &
       &                         , ilmt_phi,nlmt_phi,nlmtt_phi,nlmta_phi &
       &                         , lmta_phi,lmtt_phi,ltp_phi,mtp_phi,taup_phi &
       &                         , iproj_phi, nlmta_add, nac_p &
       &                         , modnrm,nac,fqwei,nlmta1,nlmta2 &
       &                         , porb, qorb, m_PP_tell_iorb_lmtt &
       &                         , nrorb, irorb, crorb &
       &                         , nlmta1_p, nlmta2_p &
       &                         , ipaw, dion_paw
  use m_Crystal_Structure,  only : op, nopr, nlpnt, additional_charge, altv_refcell
  use m_Kpoints,            only : kv3,vkxyz, kv3_ek,vkxyz_ek,k_symmetry, vkxyz_refcell
  use m_Ionic_System,       only : ntyp,iatom, natm, iwei, ityp, pos, cps, qex &
       &                         , if_pdos, speciesname, iproj_group
  use m_FFT,                only : nfft &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work &
#ifdef MPI_FFTW
       &                         , fft_box_size_WF, fft_box_size_WF_prev &
       &                         , afft_mpifftw, afft_mpifftw_kimg1
#else
       &                         , fft_box_size_WF, fft_box_size_WF_prev
#endif
  use m_Files,              only : nfout
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only : nspin,ipri, iprievdff, iprieigenvalue, ipribetar &
       &                         , kimg, neg, meg, neg_previous, num_extra_bands, printable &
       &                         , af, ekmode, delta_eigenvalue, evaluation_eko_diff &
       &                         , delta_eigenvalue_conduction, delta_eigenvalue_cond_is_given &
       &                         , m_CtrlP_ntcnvg_incre &
       &                         , m_CtrlP_ntcnvg_clear &
       &                         , m_CtrlP_ntcnvg_reset &
       &                         , m_CtrlP_cachesize &
       &                         , m_CntrlP_set_neg, m_CntrlP_set_meg &
       &                         , sw_orb_popu, sw_use_add_proj, sw_rsb &
       &                         , sw_serial_fft &
       &                         , sw_modified_kpoint_increment &
       &                         , sw_communicator_for_chg &
       &                         , sw_rspace &
#ifdef MPI_FFTW
       &                         , sw_mpi_fftw &
#endif
#ifdef SAVE_FFT_TIMES
       &                         , sw_save_fft &
#endif
       &                         , sw_hybrid_functional, sw_fft_xzy &
       &                         , sw_interpolate_wfs &
       &                         , sw_keep_hloc_phi &
#ifndef ENABLE_ESM_PACK
       &                         , nblocksize_mgs, nblocksize_mgs_is_given
#else
       &                         , nblocksize_mgs, nblocksize_mgs_is_given, esm_qbac
#endif
  use m_Const_Parameters,   only : DP, CMPLDP, SKIP, EXECUT, ON, OFF, INVERSE, DELTA, PAI2 &
       &                         , DELTAevdff, BUCS, SCF, EK, CARTS, DIRECT &
       &                         , ELECTRON, GAMMA, EK_CONVERGED, STORED_AND_NEW, OLD
  use m_Parallelization,    only : MPI_CommGroup,ista_kngp,iend_kngp,npes,mype &
       &                         , np_kngp, mp_kngp, nel_kngp, myrank_g &
       &                         , nrank_e,myrank_e,map_e,ista_e,iend_e,istep_e,idisp_e &
       &                         , map_z,np_e,mpi_k_world,myrank_k,map_k,ista_k,iend_k &
       &                         , ista_atm, iend_atm &
       &                         , ierr,mp_e,nel_e &
       &                         , ista_g1k,iend_g1k,np_g1k,mp_g1k &
       &                         , np_fs, mp_fs,  nel_fs &
       &                         , ista_kngp_gw, iend_kngp_gw, mpi_g_world &
       &                         , np_kngp_gw, mp_kngp_gw, nel_kngp_gw &
       &                         , ista_g1k_prev, iend_g1k_prev, nrank_k, nis_kv3_ek &
       &                         , np_fft_y, np_fft_z, ista_spin, iend_spin
  use m_Control_Parameters,  only : nblocksize_fftw, nblocksize_fftw_is_given        &
       &                          , nblocksize_gather_f_is_given                     &
       &                          , nblocksize_gather_f
  use m_Parallelization,     only : mpi_kg_world, mpi_ke_world, mpi_chg_world, mpi_ske_world &
       &                          , nrank_g &
       &                          , neg_g, neg_g_all                            &
       &                          , nel_fft_z, nel_fft_y, nel_fft_x &
       &                          , fft_X_x_dim, fft_X_y_dim, fft_X_z_dim &
       &                          , fft_X_x_nel, fft_X_y_nel, fft_X_z_nel &
       &                          , xyz_fft_x, xyz_fft_y, xyz_fft_z,lrank
! === For epsmain by tkato 2013/11/14 ==========================================
  use m_Parallelization,     only : mpi_ge_world
! ==============================================================================
  use m_Parallelization,     only : fft_Y_x_dim,nis_fft_Y_x,nie_fft_Y_x
  use m_PlaneWaveBasisSet,   only : m_pwBS_kinetic_energies_3D
#ifdef FFT_3D_DIVISION
  use m_FFT,                 only : m_FFT_WF_3DIV_3D, m_FFT_Inverse_3DIV_3D &
       &                          , m_FFT_W_Vlocal_W_3DIV_3D
#else
  use m_FFT,                 only : m_FFT_WF_3D,     m_FFT_Inverse_3D       &
       &                          , m_FFT_WF_XYZ_3D, m_FFT_Inverse_XYZ_3D   &
       &                          , m_FFT_Inverse_XYZ_3D_oo_place           &
       &                          , m_FFT_W_Vlocal_W_3D                     &
       &                          , m_FFT_W_Vlocal_W_XYZ_3D
#endif
  use m_FFT,                 only : m_FFT_WF
! ============================== added by K. Tagami ========== 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                              SpinOrbit_mode, sw_hubbard
  use m_Const_Parameters,     only : BuiltIn, Neglected

  use m_PseudoPotential,      only : dion_scr_noncl, fqwei_noncl
!!!!!  use m_FFT,  only : m_FFT_W_Vlocal_W_noncl
  use m_ES_NonCollinear,      only : m_ES_MagMom_to_DensMat_Gspace, &
       &                            m_ES_MagMom_to_DensMat_vlhxcl, &
       &                             m_ES_DensMat_To_MagMom_dm, &
       &                             m_ES_DensMat_To_MagMom_porb
! ============================================================ 11.0

! ============================== added by K. Tagami ========== 12.0Exp
  use m_Control_Parameters,    only :   fixed_charge_k_parallel, sw_band_unfolding, &
       &                               band_unfolding_active
  use m_Const_Parameters,      only :   ONE_BY_ONE
! ============================== added by K. Tagami ========== 12.0Exp

! ========== KT_add ======= 13.0U2
  use m_Control_Parameters, only : sw_potential_mixing, sw_mix_charge_hardpart, &
       &                           sw_modified_TFW_functional, use_metagga, vtau_exists
  use m_PseudoPotential,    only : flg_paw, nloc
! ========================= 13.0U2
  use m_Kpoints,       only : sw_force_kpt_inside_bz

! === KT_add ==== 2015/02/21
  use m_Ionic_System,  only : ionic_charge_atomtyp, ionic_charge_atoms, &
       &                      mag_moment0_atoms_is_defined, sw_set_initial_magmom_by_atom
! =============== 2015/02/21

#ifdef MPI_FFTW
  use m_FFT, only : m_FFT_W_Vlocal_W_mpifftw,m_FFT_W_Vlocal_W_mpifftw3d
#endif
  use mpi

  implicit none

#ifdef SAVE_FFT_TIMES
#ifdef DEBUG_SAVE_FFT_TIMES
  integer, parameter :: ipri_save_fft = 1
#else
  integer, parameter :: ipri_save_fft = 2
#endif
#endif

  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: hlocphi_l   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l_buf   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l_prev
#ifdef SAVE_FFT_TIMES
  real(kind=DP),allocatable, dimension(:,:,:)  :: Phifftr_l
  !                                      d(nfft,np_e,ista_k:iend_k) or d(lsize,np_e,ista_k:iend_k)
  integer, allocatable, dimension(:,:) :: status_saved_phifftr ! d(np_e,ista_k:iend_k) 0: nothing or old, 1: stored and new
#endif
  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_ball
  real(kind=DP),allocatable, dimension(:,:,:)  :: zah_ball
!!$  integer,      allocatable, dimension(:)     :: symmetric_or_antisymmetric
!!$                               !d(np_e), this is allocated when kimg == 1 and the system is symmetric
!!$  real(kind=DP),allocatable, dimension(:,:)   :: zaj_gamma_neg
!!$                               !d(kg1,np_e), zaj -G part for Gamma points when kimg = 1 and the system is symmetric
  integer,      allocatable, dimension(:,:)   :: neordr, nrvf_ordr !d(neg,ista_k:iend_k)
  real(kind=DP),allocatable, dimension(:,:)   :: eko_l             !d(np_e,ista_k:iend_k)
  real(kind=DP),allocatable, dimension(:,:)   :: eko_ek            !d(neg,kv3_ek)
  integer,  allocatable, dimension(:)         :: iconv_ek          !d(kv3_ek)
  real(kind=DP),allocatable, dimension(:,:)   :: occup_l           !d(np_e,ista_k:iend_k)
  real(kind=DP)                               :: efermi, efermi_spin(2), vbm = -9.99d10
  logical ::                                     metalic_system = .false.
  logical ::                                     check_if_metalic_flag = .false.
  real(kind=DP)                               :: totch
  real(kind=DP)                               :: band_entropy = 0.d0

  real(kind=DP),allocatable, dimension(:,:,:) :: fsr_l,fsi_l!d(np_e,nlmta,ista_k:iend_k)
  integer,      allocatable, dimension(:,:)   :: neordr_old  !d(neg,ista_k:iend_k)
  real(kind=DP),allocatable, dimension(:,:,:) :: fsr_ball,fsi_ball
  real(kind=DP),allocatable, dimension(:,:,:) :: fsr_add_l,fsi_add_l!d(np_e,nlmta_add,ista_k:iend_k)
  real(kind=DP),allocatable, dimension(:,:,:,:) :: compr_l,compi_l!d(np_e,nlmta_phi,nopr,ista_k:iend_k)
  real(kind=DP),allocatable, dimension(:,:,:,:) :: compr_l_ek,compi_l_ek!d(np_e,nlmta_phi,nopr,ista_k:iend_k)

!  real(kind=DP), allocatable, dimension(:,:,:) :: vlhxc_l !d(ista_kngp:iend_kngp,kimg,nspin)
  real(kind=DP), target, allocatable, dimension(:,:,:) :: vlhxc_l
                                               !d(ista_kngp:iend_kngp,kimg,nspin)
  real(kind=DP), target, allocatable, dimension(:,:,:) :: vlhxc_l_old

! ------------
! meta-gga
  real(kind=DP), target, allocatable, dimension(:,:,:) :: vtau_l
  real(kind=DP), target, allocatable, dimension(:,:) :: vtau_phl
! ------------

! =========================================== added by K. Tagami ============ 11.0
  real(kind=DP), allocatable, dimension(:,:,:) :: vlhxc_ssrep
                                     !d(ista_kngp:iend_kngp,kimg,ndim_chgpot)
  real(kind=DP), allocatable, dimension(:,:,:,:):: dhub_aimag
! =========================================================================== 11.0

  real(kind=DP), allocatable, dimension(:,:,:,:):: vlhxcQ  !d(nlmt,nlmt,natm,nspin)
! ===== KT_add ========= 13.0U2
  real(kind=DP), allocatable, dimension(:,:,:,:) :: vlhxcQ_old
! ====================== 13.0U2

  real(kind=DP), allocatable, dimension(:,:,:,:):: dhub  !d(nlmt,nlmt,natm,nspin)
  real(kind=DP), allocatable, dimension(:,:,:) :: vnlph_l  !d(kg1,np_e,kimg) work array

  real(kind=DP),private,allocatable,dimension(:,:) :: eko1_l  !d(np_e,ista_k:iend_k)
  real(kind=DP),private,allocatable,dimension(:)   :: evdff   !d(3)
  real(kind=DP),private,allocatable,dimension(:)   :: evdffr  !d(3)

  real(kind=DP), allocatable, dimension(:,:)       :: chg_softpart
  logical                                          :: chg_has_been_calculated = .false.
  ! evdff(1) : dsqrt((sum((e_old - e_new)**2))/(sum))
  ! evdff(2) : (1/sum)sum(dabs(e_old - e_new))
  ! evdff(3) : (1/sum)sum(dsqrt(dabs(e_old**2 - e_new**2)))

  integer,private                                     :: NB
#ifdef SX
  integer,        parameter                           :: nblocksize_mgs_default = 200
#else
!f  integer,        parameter                           :: nblocksize_mgs_default = 8
  integer,parameter                                   :: nblocksize_mgs_default = 8
#endif

  real(kind=DP),private,allocatable,dimension(:)      :: ar, ai

  real(kind=DP),        allocatable, dimension(:)     :: afft, bfft
  integer, private, parameter                         :: sw_timing_2ndlevel = ON

#ifndef NO_NONLOCAL_DGEMM
  logical :: DGEMM_DEBUG = .false.
  real(kind=DP), allocatable, dimension(:,:) :: fsr_l_2D,fsi_l_2D
  real(kind=DP), allocatable, dimension(:,:) :: fsr_diff_l_2D,fsi_diff_l_2D
  real(kind=DP), allocatable, dimension(:,:,:) :: fsr_gall,fsi_gall
  real(kind=DP), allocatable, dimension(:,:) :: pre_sc_without, pre_ss_without
  real(kind=DP), allocatable, dimension(:,:,:) :: fsr_l_2D_k,fsi_l_2D_k
#endif

  complex(kind=CMPLDP), allocatable, dimension(:) :: vloc_esm

  logical                   :: firstcall = .true.
#ifdef MPI_FFTW
  logical                   :: firstcall_mpifftw = .true.
#endif
  integer, allocatable,  dimension(:)   :: fftscnt0, fftrcnt0, fftindex0, fftdist0
  integer, allocatable,  dimension(:,:) :: fftsend0, fftrecv0
  integer :: fftmaxrecv0, fftmaxsend0
#ifdef MPI_FFTW
  integer, allocatable,  dimension(:)   :: fftscnt0_mpifftw, fftrcnt0_mpifftw, fftindex0_mpifftw, fftdist0_mpifftw
  integer, allocatable,  dimension(:,:) :: fftsend0_mpifftw, fftrecv0_mpifftw
  integer :: fftmaxrecv0_mpifftw, fftmaxsend0_mpifftw
  integer,  allocatable, dimension(:) :: mmx,mmy,mmz
  integer :: maxrecv
#endif

#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),        allocatable, dimension(:)     :: bfft_kin
  real(kind=DP),        allocatable, dimension(:)     :: bfft_hartr
  real(kind=DP), dimension(6)                         :: ttr
! -----  ascat ceases modifying  -----
#endif

!  1-1.  m_ES_alloc_zaj_etc
!  1-2.  m_ES_alloc_eko1
!  1-3.  m_ES_alloc_eko_ek
!  1-4.  m_ES_alloc_vlhxc
!  1-5.  m_ES_dealloc_vlhxc
!  1-6.  m_ES_alloc_vlhxcQ
!  1-7.  m_ES_dealloc_vlhxcQ
!  1-8.  m_ES_alloc_Dhub
!  1-9.  m_ES_dealloc_Dhub
!  1-10. m_ES_alloc_scss_etc
!  1-11. m_ES_dealloc_scss_etc
!  1-12. m_ES_alloc_fft_related
!  1-13. m_ES_dealloc_fft_related
!  1-14. m_ES_alloc_afft_scss_etc
!  1-15. m_ES_dealloc_afft_scss_etc
!  1-16. alloc_zfsincos_mpi
!  1-17. dealloc_zfsincos_mpi
!  1-18. m_ES_mgs_alloc
!        - mgs_vdb_alloc, - mgs_nrc_alloc
!  1-19. m_ES_mgs_dealloc
!        - mgs_vdb_dealloc, - mgs_nrc_dealloc
!  1-20. m_ES_alloc_zfsincos
!  1-21. m_ES_dealloc_zfsincos
!  1-22. m_ES_alloc_arai
!  1-23. m_ES_dealloc_arai
!  2-1.  m_ES_gtotch
!  2-2.  m_ES_wd_zaj_small_portion
!  2-3.  m_ES_decide_precon_factor
!  2-4.  kinetic_energy               <- (2-3)
!  2-5.  m_ES_Vnonlocal_W
!        - calc_phase_mpi, - Vnonlocal_W_part_sum_over_lmt1
!        - add_vnlph_l_with_eko_part, - add_vnlph_l_without_eko_part
!  2-6.  m_ES_sort_eigen_values
!        - cp_eigen_values_for_af, - expand_neordr_and_nrvf_ordr
!        - heap_sorting
!  2-7.  m_ES_energy_eigen_values
!        - get_ipri0
!  2-8.  m_ES_energy_eigen_values_ext
!        - get_ipri0
!  2-9.  m_ES_eigen_values_for_each_k
!        - W_T_W, - W_Vnonlocal_W
!  2-10. m_ES_eigen_values_for_each_kex
!        - W_T_W, - W_Vnonlocal_W
!  2-11. m_ES_Vlocal_in_Rspace
!        - map_vlhxc_l_onto_afft
!  2-12. m_ES_WF_in_Rspace
!  2-13. G_dot_R_mpi                  <- (2-14), (2-20), (2-23), (2-26)
!  2-14. m_ES_betar_dot_WFs           -> (2-15), (2-16)
!  2-15. wd_fsr_fsi                   <- (2-14), (2-17)
!  2-16. G_dot_R_map                  <- (2-14), (2-17), (2-18), (2-20), (2-23), (2-26)
!  2-17. m_ES_betar_dot_WFs_4_each_k  -> (2-15), (2-16)
!  2-18. m_ES_betar_dot_Psi_4_each_k  -> (2-16)
!  2-19. m_ES_betar_dot_WFs_4_lmta_k
!       - G_dot_R_mult_snl, - betar_dot_WFs_core, - multiple_i_l
!       - betar_dot_WFs_core2
!  2-20. m_ES_phir_dot_WFs            -> (2-13), (2-16)
!  2-21. wd_compr_compi
!  2-22. m_ES_phir_dot_WFs_4_lmta_k
!       - G_dot_R_mult_phig, - phir_dot_WFs_core, - multiple_i_l, - phir_dot_WFs_core2
!  2-23. m_ES_add_betar_dot_WFs       -> (2-13), (2-16), (2-24)
!  2-24. wd_fsr_fsi_add               <- (2-23)
!  2-25. m_ES_add_betar_dot_WFs_4_lmta_k
!        - G_dot_R_mult_snl, - betar_dot_WFs_core, - multiple_i_l, - betar_dot_WFs_core2
!  2-26. m_ES_PAO_WFs                 -> (2-13), (2-16)
!  2-27. m_ES_PAO_WFs_4_lmta_k
!        - G_dot_R_mult_paog
!  2-28. m_ES_modified_gram_schmidt
!  2-29. m_ES_orthogonalize_SD_to_WFs -> (2-30)
!  2-30. mgs_sd2wf_each_k_G           <- (2-29)  -> (2-34), (2-44)
!        - broadcast_fs, - Psi1SPhi2_t, - modify_bsd_and_phi_t, - alloc_phi_w_and_brd_w
!        - dealloc_phi_w_and_brd_w, - WSW, - normalize_bsd_and_phi, - Psi1SPhi2
!        - modify_bsd_and_phi
!  2-31. orthogonalize_SD             -> (2-33)
!  2-32. m_ES_MGS_4_each_k            -> (2-17), (2-33)
!        - wd_title_of_the_operation
!  2-33. mgs_4_each_k_G               <- (2-31), (2-32)  -> (2-34), (2-44)
!        - WSW_t, - normalize_bp_and_psi_t, - W1SW2_t_r
!        - modify_bp_and_psi_t_r, - substitute_jto_ib2back, - W1SW2_t
!        - modify_bp_and_psi_t, - alloc_and_brd_w, - dealloc_and_brd_w
!        - WSW, - normalize_bp_and_psi, - W1SW2, - modify_bp_and_psi
!  2-34. set_npzri                    <- (2-30), (2-33)
!  2-35. m_ES_W_transpose
!  2-36. m_ES_W_transpose_back
!  2-37. m_ES_W_transpose_back2
!  2-38. m_ES_W_transpose_r
!  2-39. m_ES_W_transpose_back_r
!  2-40. m_ES_F_transpose
!  2-41. m_ES_F_transpose_back
!  2-42. m_ES_F_transpose_r
!  2-43. m_ES_F_transpose_back_r
!  2-44. cp_bp_and_psi_2_brd_w2      <- (2-30), (2-33)
!  2-45. m_ES_wd_eko
!  2-46. m_ES_wd_eko_cond
!  2-47. m_ES_wd_eko2
!  2-48. m_ES_sum_of_LocalPart
!  2-48. m_ES_sum_of_LocalPart2
!  2-49. m_ES_cpeko
!  2-50. logical fnct. m_ES_eekdif
!        - get_ipri0
!  2-51. logical fnct. m_ES_eekdif_cond
!  2-52. m_ES_cp_eko_l_to_eko_ek
!  2-53. real(DP) fnct. m_ES_what_is_evdff_now
!  2-54. real(kind=DP) fnct. m_ES_get_energy
!  2-55. m_ES_orbital_population
!  2-56. m_ES_sym_comp
!  2-57. m_ES_orbital_den_mat
!  2-58. m_ES_set_num_bands_super
!

#ifdef MPI_FFTW
  include 'fftw3-mpi.f03'
#endif
!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  interface m_ES_WF_in_Rspace_3D
    module procedure m_ES_WF_in_Rspace_3D0
    module procedure m_ES_WF_in_Rspace_3D1
    module procedure m_ES_WF_in_Rspace_3D2
  end interface m_ES_WF_in_Rspace_3D

  interface m_ES_WF_2D
    module procedure m_ES_WF_2D0
    module procedure m_ES_WF_2D1
  end interface m_ES_WF_2D

contains

  subroutine m_ES_alloc_fsri_l_2D(ik)
    integer, intent(in) :: ik
    if(.not.allocated(fsr_l_2D)) allocate(fsr_l_2D(np_e,nlmta)); fsr_l_2D = 0.0d0
    if( k_symmetry(ik) /= GAMMA ) then
       if(.not.allocated(fsi_l_2D)) allocate(fsi_l_2D(np_e,nlmta)); fsi_l_2D = 0.0d0
    end if
    if(sw_rspace == ON)then
       if(allocated(fsr_l_2D_k)) deallocate(fsr_l_2D_k)
       if( k_symmetry(ik) /= GAMMA ) then
         if(allocated(fsi_l_2D_k)) deallocate(fsi_l_2D_k)
       endif
    endif
  end subroutine m_ES_alloc_fsri_l_2D

  subroutine m_ES_dealloc_fsri_l_2D(ik)
    integer, intent(in) :: ik
    if(allocated(fsr_l_2D)) deallocate(fsr_l_2D)
    if( k_symmetry(ik) /= GAMMA) then
       if(allocated(fsi_l_2D)) deallocate(fsi_l_2D)
    end if
    if(sw_rspace == ON)then
       if(allocated(fsr_l_2D_k)) deallocate(fsr_l_2D_k)
       if( k_symmetry(ik) /= GAMMA ) then
         if(allocated(fsi_l_2D_k)) deallocate(fsi_l_2D_k)
       endif
    endif
  end subroutine m_ES_dealloc_fsri_l_2D

  subroutine m_ES_alloc_fsri_diff_l_2D(ik)
    integer, intent(in) :: ik
    if(.not.allocated(fsr_diff_l_2D)) allocate(fsr_diff_l_2D(np_e,nlmta)); fsr_diff_l_2D = 0.0d0
    if( k_symmetry(ik) /= GAMMA ) then
       if(.not.allocated(fsi_diff_l_2D)) allocate(fsi_diff_l_2D(np_e,nlmta)); fsi_diff_l_2D = 0.0d0
    end if
  end subroutine m_ES_alloc_fsri_diff_l_2D

  subroutine m_ES_dealloc_fsri_diff_l_2D(ik)
    integer, intent(in) :: ik
    if(allocated(fsr_diff_l_2D)) deallocate(fsr_diff_l_2D)
    if( k_symmetry(ik) /= GAMMA) then
       if(allocated(fsi_diff_l_2D)) deallocate(fsi_diff_l_2D)
    end if
  end subroutine m_ES_dealloc_fsri_diff_l_2D


  subroutine m_ES_realloc_zaj()
    if(allocated(zaj_l)) deallocate(zaj_l)
    allocate( zaj_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg) );zaj_l = 0.0d0
    if(sw_keep_hloc_phi==ON) then
      if(allocated(hlocphi_l)) deallocate(hlocphi_l)
      allocate(hlocphi_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg));hlocphi_l = 0.d0
    endif
  end subroutine m_ES_realloc_zaj

  subroutine m_ES_alloc_zaj_etc()
    integer :: ik,ib
#ifdef SAVE_FFT_TIMES
    integer :: lsize
#endif
    if(iend_k - ista_k < 0) call phase_error_with_msg(nfout,' iend_k - ista_k < 0 (in m_ES_alloc_zaj_etc)', &
    __LINE__,__FILE__)
    if(nlmta <= 0) call phase_error_with_msg(nfout,' nlmta <= 0 (in m_ES_alloc_zaj_etc)',__LINE__,__FILE__)
    if(nblocksize_mgs_is_given) then
       NB = nblocksize_mgs
    else
       NB = nblocksize_mgs_default
    end if
    allocate( zaj_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg) );zaj_l = 0.0d0
    if(sw_keep_hloc_phi==ON) then
      allocate(hlocphi_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg)); hlocphi_l = 0.d0
    endif
    allocate( zaj_ball(maxval(np_g1k),neg,ista_k:iend_k,kimg) )
    allocate( zah_ball(maxval(np_g1k),neg,kimg) )
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
!!$    ibsize = 1
!!$    if (nblocksize_fftw_is_given) then
!!$       ibsize = nblocksize_fftw
!!$       if (ibsize < 1) ibsize = 1
!!$    endif
#ifdef FFT_3D_DIVISION
       lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
       lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
       allocate(Phifftr_l(lsize*kimg,np_e,ista_k:iend_k)); Phifftr_l = 0.d0
       allocate(status_saved_phifftr(np_e,ista_k:iend_k)); status_saved_phifftr = 0
    end if
#endif

    allocate(neordr(neg,ista_k:iend_k))
    allocate(nrvf_ordr(neg,ista_k:iend_k))
    do ik = ista_k, iend_k
       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)
       nrvf_ordr(1:neg,ik) = (/(ib,ib=1,neg)/)
    end do
    allocate(neordr_old(neg,ista_k:iend_k))
    allocate(fsr_l(np_e,np_fs,ista_k:iend_k));fsr_l =0.0d0
    allocate(fsr_ball(neg,np_fs,ista_k:iend_k));fsr_ball = 0.0d0
    allocate(occup_l(np_e,ista_k:iend_k)); occup_l = 0.d0
    if(sw_rspace == ON)then
       allocate(fsr_l_2D_k(np_e,nlmta,ista_k:iend_k)); fsr_l_2D_k = 0.0d0
       allocate(fsi_l_2D_k(np_e,nlmta,ista_k:iend_k)); fsi_l_2D_k = 0.0d0
    endif
!!$ASASASAS
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi_l(np_e,np_fs,ista_k:iend_k)); fsi_l = 0.d0
       allocate(fsi_ball(neg,np_fs,ista_k:iend_k));fsi_ball = 0.0d0
    else
       allocate(fsi_l(1,1,1))
       allocate(fsi_ball(1,1,1))
    end if
    fsi_l = 0.d0
    allocate(vnlph_l(maxval(np_g1k),np_e,kimg)) ; vnlph_l = 0.0d0
    allocate(eko_l(np_e,ista_k:iend_k))      ; eko_l = 0.0d0
    ! -- for PDOS
    if(sw_orb_popu == ON) then
       allocate(compr_l(np_e,nlmta_phi,nopr,ista_k:iend_k))
       compr_l=0.d0         !ASMS 2016/09/08

       if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          allocate(compi_l(np_e,nlmta_phi,nopr,ista_k:iend_k))
!!$ASASASAS
!!%       end if
       else
          allocate(compi_l(1,1,1,1))
       end if
       compi_l = 0.d0
!!$ASASASAS
       if(ekmode == ON) then
          allocate(compr_l_ek(np_e,nlmta_phi,nopr,kv3_ek))
          if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
             allocate(compi_l_ek(np_e,nlmta_phi,nopr,kv3_ek))
          else
             allocate(compi_l_ek(1,1,1,1))
          end if
       endif
    end if
    if(sw_use_add_proj == ON) then
       allocate(fsr_add_l(np_e,nlmta_add,ista_k:iend_k))
       if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
          allocate(fsi_add_l(np_e,nlmta_add,ista_k:iend_k))
       end if
    end if

  end subroutine m_ES_alloc_zaj_etc


  subroutine m_ES_alloc_eko1()
    allocate(eko1_l(np_e,ista_k:iend_k)); eko1_l = 0.d0
    allocate(evdff(3));evdff=1.d99
    allocate(evdffr(3))
  end subroutine m_ES_alloc_eko1

  subroutine m_ES_alloc_eko_ek()
    if ( allocated( eko_ek ) ) deallocate( eko_ek )
    allocate(eko_ek(neg,kv3_ek)); eko_ek = 0.0d0

    if(.not.allocated(iconv_ek)) then
       allocate(iconv_ek(kv3_ek)); iconv_ek = 0
    end if
  end subroutine m_ES_alloc_eko_ek

  subroutine m_ES_alloc_fsr_l_2d(n1, n2)
    integer, intent(in) :: n1, n2
!FNS--debug
    if(allocated(fsr_l_2D)) deallocate(fsr_l_2D)
!----------
    allocate(fsr_l_2D(n1,n2))
    fsr_l_2D = 0.0d0
  end subroutine m_ES_alloc_fsr_l_2d

  subroutine m_ES_alloc_fsi_l_2d(n1, n2)
    integer, intent(in) :: n1, n2
!FNS--debug
    if(allocated(fsi_l_2D)) deallocate(fsi_l_2D)
!----------
    allocate(fsi_l_2D(n1, n2))
    fsi_l_2D = 0.0d0
  end subroutine m_ES_alloc_fsi_l_2d

  subroutine m_ES_alloc_fsr_diff_l_2d(n1, n2)
    integer, intent(in) :: n1, n2
!FNS--debug
    if(allocated(fsr_diff_l_2D)) deallocate(fsr_diff_l_2D)
!----------
    allocate(fsr_diff_l_2D(n1,n2))
    fsr_l_2D = 0.0d0
  end subroutine m_ES_alloc_fsr_diff_l_2d

  subroutine m_ES_alloc_fsi_diff_l_2d(n1, n2)
    integer, intent(in) :: n1, n2
!FNS--debug
    if(allocated(fsi_diff_l_2D)) deallocate(fsi_diff_l_2D)
!----------
    allocate(fsi_diff_l_2D(n1, n2))
    fsi_diff_l_2D = 0.0d0
  end subroutine m_ES_alloc_fsi_diff_l_2d

  subroutine m_ES_dealloc_fsr_l_2d()
    if(allocated(fsr_l_2D)) deallocate(fsr_l_2D)
  end subroutine m_ES_dealloc_fsr_l_2d

  subroutine m_ES_dealloc_fsi_l_2d()
    if(allocated(fsi_l_2D)) deallocate(fsi_l_2D)
  end subroutine m_ES_dealloc_fsi_l_2d

  subroutine m_ES_dealloc_fsr_diff_l_2d()
    if(allocated(fsr_diff_l_2D)) deallocate(fsr_diff_l_2D)
  end subroutine m_ES_dealloc_fsr_diff_l_2d

  subroutine m_ES_dealloc_fsi_diff_l_2d()
    if(allocated(fsi_diff_l_2D)) deallocate(fsi_diff_l_2D)
  end subroutine m_ES_dealloc_fsi_diff_l_2d

  subroutine m_ES_cp_iconv(numk,iconv_ek_tmp)
    integer, intent(in) :: numk
    integer, intent(in), dimension(numk) :: iconv_ek_tmp
    integer :: i
    if(numk>=1) then
       if(.not.allocated(iconv_ek)) then
          allocate(iconv_ek(kv3_ek))
          iconv_ek = 0
       end if
       do i = 1, numk
          iconv_ek(i) = iconv_ek_tmp(i)
       end do
       if(ipri >= 1) then
          write(nfout,'(" ! -- iconv_ek -- <<m_ES_cp_iconv>>")')
          write(nfout,'(" ! ",8i6)') (iconv_ek(i),i=1,numk)
          write(nfout,'(" ! ",8i6)') (iconv_ek(i),i=numk+1,kv3_ek)
       end if
    end if
  end subroutine m_ES_cp_iconv

  subroutine m_ES_alloc_vlhxc()
    allocate(vlhxc_l(ista_kngp:iend_kngp,kimg,nspin)); vlhxc_l = 0.d0
   if ( use_metagga .and. vtau_exists ) then
      allocate(vtau_l(ista_kngp:iend_kngp,kimg,nspin)); vtau_l = 0.d0
   endif
  end subroutine m_ES_alloc_vlhxc

  subroutine m_ES_dealloc_vlhxc()
    if (allocated(vlhxc_l)) deallocate(vlhxc_l)
    if (allocated(vlhxc_l_old)) deallocate(vlhxc_l_old)
    if (allocated(vtau_l))    deallocate(vtau_l)
  end subroutine m_ES_dealloc_vlhxc

  subroutine m_ES_alloc_vlhxcQ()
    if(nlmt <= 0) call phase_error_with_msg(nfout,' nlmt <=  0 (in m_ES_alloc_vlhxcQ)',__LINE__,__FILE__)
! =================================== modified by K. Tagami ================ 11.0
!    allocate(vlhxcQ(nlmt,nlmt,natm,nspin)); vlhxcQ = 0.d0
!
    if ( noncol ) then
       allocate(vlhxcQ(nlmt,nlmt,natm,ndim_magmom))
    else
       if(.not.allocated(vlhxcQ)) allocate(vlhxcQ(nlmt,nlmt,natm,nspin))
    endif
    vlhxcQ = 0.d0
! ========================================================================== 11.0

#if 0
! ====== KT_add ====== 13.0U2
    if ( sw_potential_mixing == ON ) then
       if ( flg_paw .or. sw_mix_charge_hardpart == ON ) then
          allocate(vlhxcQ_old(nlmt,nlmt,natm,nspin))
          vlhxcQ_old = 0.d0
       endif
   endif
! ==================== 13.0U2
#endif

  end subroutine m_ES_alloc_vlhxcQ

  subroutine m_ES_dealloc_vlhxcQ()
    if(allocated(vlhxcQ))  deallocate(vlhxcQ)
#if 0
! ======== KT_add ==== 13.0U2
    if(allocated(vlhxcQ_old))  deallocate(vlhxcQ_old)
! ==================== 13.0U2
#endif
  end subroutine m_ES_dealloc_vlhxcQ

! ====== KT_add =================== 13.0U2
  subroutine m_ES_cp_vlhxc_to_old
    vlhxc_l_old = vlhxc_l
    if ( flg_paw .or. sw_mix_charge_hardpart == ON) then
       vlhxcQ_old = vlhxcQ
#if 0
       if ( noncol ) vlhxcQ_i_old = vlhxcQ_i
#endif
    endif
  end subroutine m_ES_cp_vlhxc_to_old

  subroutine m_ES_cp_vlhxcQ_to_old
    if ( flg_paw .or. sw_mix_charge_hardpart == ON) then
       vlhxcQ_old = vlhxcQ
#if 0
       if ( noncol ) vlhxcQ_i_old = vlhxcQ_i
#endif
    endif
  end subroutine m_ES_cp_vlhxcQ_to_old
! ==================================== 13.0U2

  subroutine m_ES_wd_vlhxcQ()
    integer :: i, j, na
    write(nfout,'(" --- m_ES_wd_vlhxcQ ---")')
    do na = 1, natm
       write(nfout,'("  na = ", i8)') na
       do j = 1, nlmt
          write(nfout,'(" j = ", i8)') j
          write(nfout,'(8f10.4)') (vlhxcQ(i,j,na,1),i=1,nlmt)
       end do
    end do
  end subroutine m_ES_wd_vlhxcQ

  subroutine m_ES_alloc_Dhub()
    if(nlmt <= 0) call phase_error_with_msg(nfout,' nlmt <=  0 (in m_ES_alloc_Dhub)',__LINE__,__FILE__)
! ======================================== modified by K. Tagami ========= 11.0
!!    allocate(dhub(nlmt,nlmt,natm,nspin)); dhub = 0.d0
!
    if ( noncol ) then
       if(allocated(dhub)) deallocate(dhub)
       if(allocated(dhub_aimag)) deallocate(dhub_aimag)
       allocate(dhub(nlmt,nlmt,natm,ndim_magmom))
       allocate(dhub_aimag(nlmt,nlmt,natm,ndim_magmom))
       dhub = 0.d0; dhub_aimag = 0.0d0
    else
       if(allocated(dhub)) deallocate(dhub)
       allocate(dhub(nlmt,nlmt,natm,nspin))
       dhub = 0.d0
    endif
! =========================================================================11.0

  end subroutine m_ES_alloc_Dhub

  subroutine m_ES_dealloc_Dhub()
    if(allocated(dhub)) deallocate(dhub)
! ======================================== added by K. Tagami ========= 11.0
    if(allocated(dhub_aimag)) deallocate(dhub_aimag)
! ===================================================================== 11.0
  end subroutine m_ES_dealloc_Dhub

  subroutine m_ES_alloc_fft_related()
    allocate(afft(nfft))
    allocate(bfft(nfft))
    call m_FFT_alloc_WF_work() ! allocate(ftw)
  end subroutine m_ES_alloc_fft_related

  subroutine m_ES_dealloc_fft_related()
    call m_FFT_dealloc_WF_work()
    deallocate(bfft)
    deallocate(afft)
  end subroutine m_ES_dealloc_fft_related

  subroutine m_ES_gtotch(nfout)
    integer, intent(in) :: nfout
    integer :: it, ia

    totch = 0.d0
    do it = 1, ntyp
       totch = totch + ival(it)*iatom(it) + qex(it)
    end do
#ifdef ENABLE_ESM_PACK
    totch = totch - esm_qbac
#endif
! ===== KT_add === 2014/06/08
    totch = totch - additional_charge      ! totch is num. of electrons
! ================ 2014/06/08

! === KT_add === 2015/02/21
    if ( sw_set_initial_magmom_by_atom == ON ) then
       if ( mag_moment0_atoms_is_defined ) then
          Do ia=1, natm
             totch = totch -ionic_charge_atoms(ia) *iwei(ia)
          End do
       else
          Do ia=1, natm
             it = ityp(ia)
             totch = totch -ionic_charge_atomtyp(it) *iwei(ia)
          end do
       endif
    else
       Do ia=1, natm
          it = ityp(ia)
          totch = totch -ionic_charge_atomtyp(it) *iwei(ia)
       end do
    endif
! ============ 2015/02/21

    if(printable) write(nfout,'(" TOTCH (total charge) = ",d25.12)') totch
    if(totch <= 1.d-20) call phase_error_with_msg(nfout,' ! illegal TOTCH value (m_ES_gtotch)',__LINE__,__FILE__)
    if(totch > natm*500.0) then
       if(printable) then
          do it=1,ntyp
             write(nfout,'(" !! it, ival, iatom, qex = ",i4,f8.4,i4,f8.4)') it,ival(it),iatom(it),qex(it)
          end do
          write(nfout,'(" ! illegal TOTCH value (m_ES_gtotch)")')
       end if
       call phase_error_with_msg(nfout,' ! illegal TOTCH value (m_ES_gtotch)',__LINE__,__FILE__)
    end if
  end subroutine m_ES_gtotch

  subroutine m_ES_wd_zaj_small_portion_3D(nfout,ik,comment,nc)

    integer,        intent(in) :: nfout, ik, nc
    character(len=nc), intent(in) :: comment

    character(len=5) :: a
    integer :: i, ib, ri, j
    real(kind=DP), allocatable, dimension(:,:) :: zaj_tmp
    integer, parameter :: NZAJSIZE = 20
    integer :: nelm

    write(nfout,*) comment
    a = "     "
    do ib = ista_e, iend_e, istep_e               ! MPI
       if(k_symmetry(ik) == GAMMA) then
          allocate(zaj_tmp(5,2)); zaj_tmp = 0.d0
          do i = 1, 5
             find_j: do j = 1, 10
                if(nbase_gamma(j,1) == i) then
                   zaj_tmp(i,1) = zaj_l(j,map_z(ib),ik,1)
                   if(kimg == 2) zaj_tmp(i,2) = zaj_l(j,map_z(ib),ik,2)
                   goto 1001
                end if
             end do find_j
             find_j2: do j = 2, 10
                if(nbase_gamma(j,2) == i) then
                   zaj_tmp(i,1) = zaj_l(j,map_z(ib),ik,1)
                   if(kimg == 2) zaj_tmp(i,2) = -zaj_l(j,map_z(ib),ik,2)
                   exit find_j2
                end if
             end do find_j2
1001         continue
          end do

          do ri = 1, kimg
             if(ri == 1 .and. kimg == 2) a = "(Re) "
             if(ri == 2) a = "(Im) "

             if(ri == 1) write(nfout,'(" eko(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
                  & ib,ik,eko_l(map_z(ib),ik),a,(zaj_tmp(i,ri),i=1,5)
             if(ri == 2) write(nfout,'(31x,a4,5e14.6)') a,(zaj_tmp(i,ri),i=1,5)
          end do
          deallocate(zaj_tmp)
       else
          do ri = 1, kimg
             if(ri == 1 .and. kimg == 2) a = "(Re) "
             if(ri == 2) a = "(Im) "

             if(ri == 1) write(nfout,'(" eko(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
                  & ib,ik,eko_l(map_z(ib),ik),a,(zaj_l(i,map_z(ib),ik,ri),i=1,5)
             if(ri == 2) write(nfout,'(31x,a4,5e14.6)') a,(zaj_l(i,map_z(ib),ik,ri),i=1,5)
          end do
       end if
       if(ipri >= 3) then
          nelm = min(iba(ik),NZAJSIZE)
          do ri = 1, kimg
             if(ri == 1 .and. kimg == 2) a = "(Re) "
             if(ri == 2) a = "(Im) "

             if(ri == 1) write(nfout,'(" eko(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
                  & ib,ik,eko_l(map_z(ib),ik),a,(zaj_l(i,map_z(ib),ik,ri),i=1,5)
             if(ri == 2) write(nfout,'(31x,a4,5e14.6)') a,(zaj_l(i,map_z(ib),ik,ri),i=1,5)
             write(nfout,'(35x,5e14.6)') (zaj_l(i,map_z(ib),ik,ri),i=6,nelm)
!!$             write(nfout,'(35x,5e14.6)') (zaj_l(i,map_z(ib),ik,ri),i=6,kg1)
          end do
       end if
    end do
!!$       if(kimg==2) a = "(Im) "
!!$       write(nfout,'(" eko(",i4,",",i3,")= ",e14.6," ",a4,5e14.6)') &
!!$            &  ib,ik,eko_l(map_z(ib),ik),a,(zaj_l(i,map_z(ib),ik,1),i=1,5)
!!$       write(nfout,455) ik,ib,eko_l(map_z(ib),ik) ! MPI
!!$       do ri = 1, kimg
!!$          if(ri == 1 .and. kimg == 2) write(nfout,*) ' (zaj real part)'
!!$          if(ri == 2)                 write(nfout,*) ' (zaj imag part)'
!!$          if(ri == 1 .and. kimg == 2) a = "(Re)"
!!$          if(ri == 2)                 a = "(Im)"
!!$          write(nfout,'(a4,5e14.6)') a,(zaj_l(i,map_z(ib),ik,ri),i=1,5) ! MPI
!!$       end do
!!$    end do
!!$455 format(' ',' ik, ib = ',2i5,' e = ',d20.12)
  end subroutine m_ES_wd_zaj_small_portion_3D

  subroutine m_ES_decide_precon_factor_3D(precon,ik,ib1,ib2,ibesize,ekin,p)
    integer, intent(in)                         :: precon,ik,ib1,ib2,ibesize
    real(kind=DP), intent(in),  dimension(np_g1k(ik)) :: ekin
    real(kind=DP), intent(out), dimension(mp_g1k(ik),ibesize) :: p

    integer       :: i, ib
    real(kind=DP) :: x, x1, x2
    real(kind=DP), dimension(ibesize) :: ektot, d_ektot
                                                  __TIMER_SUB_START(302)
    if(precon == ON) then
       call kinetic_energy_3D(ik,ib1,ib2,ibesize,ekin,ektot)   ! -(m_E.S.)
       d_ektot = 1.d0/ektot
                                                  __TIMER_DO_START(310)
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             x = ekin(i-ista_g1k(ik)+1)*d_ektot(ib-ib1+1)
             x1 = 27 + ( 18 + (12 + 8*x) *x) *x
             x2 = 16*(x*x)*(x*x)
             p(i-ista_g1k(ik)+1,ib-ib1+1)  = x1/(x1 + x2 )
          end do
       end do
                                                  __TIMER_DO_STOP(310)
    else
       p = 1.d0
    end if
                                                  __TIMER_SUB_STOP(302)
  end subroutine m_ES_decide_precon_factor_3D

  subroutine kinetic_energy_3D(ik,ib1,ib2,ibesize,dekin,ektot)
    integer, intent(in) :: ik, ib1, ib2, ibesize
    real(kind=DP), intent(in), dimension(np_g1k(ik)) :: dekin
    real(kind=DP), intent(out), dimension(ibesize) :: ektot
    real(kind=DP)             , dimension(ibesize) :: ektot_mpi
    integer  :: i, ib
                                                  __TIMER_SUB_START(303)
    ektot(:) = 0.d0
                                                  __TIMER_DO_START(311)
    if(kimg == 1) then
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             ektot(ib-ib1+1) = ektot(ib-ib1+1) + dekin(i-ista_g1k(ik)+1)* &
            &                                    zaj_l(i-ista_g1k(ik)+1,ib,ik,1)**2
          end do
       end do
    else
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             ektot(ib-ib1+1) = ektot(ib-ib1+1) + dekin(i-ista_g1k(ik)+1)* &
            &                                    ( zaj_l(i-ista_g1k(ik)+1,ib,ik,1)**2 &
            &                                    + zaj_l(i-ista_g1k(ik)+1,ib,ik,2)**2 )
          end do
       end do
    end if
                                                 __TIMER_DO_STOP(311)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,312)
    call mpi_allreduce(ektot,ektot_mpi,ibesize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
                                                 __TIMER_COMM_STOP(312)
    ektot = ektot_mpi

    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
                                                 __TIMER_SUB_STOP(303)
  end subroutine kinetic_energy_3D








!!$  integer function cachesize(level)
!!$    integer, intent(in) :: level
!!$    if(level == 1) then
!!$       cachesize = 64
!!$    else if(level == 2) then
!!$       cachesize = 512
!!$    else if(level == 3) then
!!$       cachesize = 1024
!!$    else
!!$       cachesize = 100000
!!$    end if
!!$  end function cachesize

#ifdef __EDA__
! -----  ascat starts modifying  -----
  subroutine m_ES_WF_kin_in_Rspace(ik,ib,bfft_kin,zajtmp)
    use m_Crystal_Structure, only : rltv

    integer, intent(in)                             :: ik, ib
    integer                                         :: i, i1, ri, nb
    integer :: ii, j, i2
    real(kind=DP), intent(inout), dimension(nfft)   :: bfft_kin
    real(kind=DP), intent(in), dimension(kg1,kimg) :: zajtmp
    real(kind=DP)                                   :: ga, gb, gc, ekin

    bfft_kin = 0.d0
    call getttr(rltv,ttr)

!!!    call map_zaj_to_fft_box(ik,kg1,nfft,zajtmp,bfft_kin)  ! ASMS mod 2024/03/19

    if (k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf(1)
          bfft_kin(i1) = 0.0d0

          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             ga = vkxyz(ik,1,BUCS) + ngabc(i,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(i,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(i,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i1 = igf(i)
             bfft_kin(i1) = ekin *zajtmp(ii,1)

             j = nbase_gamma(ii,2)
             ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i2 = igf(j)
             bfft_kin(i2) = ekin *zajtmp(ii,1)
          end do
       else if ( kimg == 2 ) then
          i1 = 2*igf(1) - 1
          bfft_kin(i1)   = 0.0d0
          bfft_kin(i1+1) = 0.0d0

          do ii = 2, iba(ik)
             i = nbase(ii,ik)

             ga = vkxyz(ik,1,BUCS) + ngabc(i,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(i,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(i,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i1 = 2 *igf(i) -1
             bfft_kin(i1  ) = ekin *zajtmp(ii,1)
             bfft_kin(i1+1) = ekin *zajtmp(ii,2)

             j = nbase_gamma(ii,2)
             ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i2 = 2*igf(j)-1
             bfft_kin(i2  ) = ekin *zajtmp(ii,1)
             bfft_kin(i2+1) = -ekin *zajtmp(ii,2)
          end do
       endif
    else
       do ri = 1, kimg
          do i = 1, iba(ik)
             nb = nbase(i,ik)
             ga = vkxyz(ik,1,BUCS) + ngabc(nb,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(nb,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(nb,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                  & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0
             
!       ekin = (ga*ga + gb*gb + gc*gc)*0.5d0
             i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)

             bfft_kin(i1) = ekin*zajtmp(i,ri)   ! MPI
!!             bfft_kin(i1) = ekin*bfft_kin(i1)
          end do
       end do
    endif

    call m_FFT_WF(ELECTRON,nfout,bfft_kin,INVERSE,ON)

  end subroutine m_ES_WF_kin_in_Rspace

#endif


  subroutine m_ES_wd_zaj_small_portion0(str,nc)
    integer, intent(in) :: nc
    character(len=nc), intent(in) :: str

    integer :: ik
    do ik = ista_k, iend_k, af+1                            ! MPI
       if(ipri>=1) call m_ES_wd_zaj_small_portion_3D(nfout,ik,str,nc)
    end do
  end subroutine m_ES_wd_zaj_small_portion0


  subroutine m_ES_wd_eko_cond(nfout,ibcm,mode)
    integer, intent(in) ::                  nfout
    integer, intent(in) ::                  mode, ibcm
    real(kind=DP), pointer, dimension(:) :: eko, eko_t ! d(neg) MPI

    integer :: ik, ib, ikp

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    if(ibcm > neg) return

    allocate(eko(neg)); allocate(eko_t(neg))           ! MPI
    write(nfout,*) '=== energy_eigen_values ==='

    ikp = 0
! =========================== modified by K. Tagami =================== 12.0Exp
!    if(mode == EK) ikp = nk_in_the_process - 1
    if (mode == EK) then
       if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
          ikp = nk_in_the_process - 1
       endif
    endif
! ======================================================================= 12.0Exp

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
    do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                 ! MPI
       eko_t = 0                                       ! MPI
       do ib = 1, neg                                  ! MPI
          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik) ! MPI
       end do                                          ! MPI
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
       if(mode == EK .or. ik == 1) then
          if(iprieigenvalue >= 1) then
             write(nfout,'(" ik = ",i7," ( ",3f10.6," ) [conduction bands]" &
                  & ,/99(4f18.10,/))') ik+ikp,vkxyz(ik,1:3,BUCS) &
                  & , (eko(neordr(ib,ik)),ib=ibcm,neg-num_extra_bands)
          end if
          if(iprieigenvalue >= 2 .and. num_extra_bands >= 1) &
               & write(nfout,'(" -- extra_bands --",/99(4f18.10,/))') &
               & (eko(neordr(ib,ik)),ib=neg-num_extra_bands+1, neg)
       else
          if(iprieigenvalue >= 2) then
             write(nfout,'(" ik = ",i7," (",3f10.6," )",/99(10f8.4,/))')&
                  &ik,vkxyz(ik,1:3,BUCS), (eko(neordr(ib,ik)),ib=ibcm,neg)
          end if
       endif
    end do
    deallocate(eko); deallocate(eko_t)                 ! MPI
  end subroutine m_ES_wd_eko_cond

  subroutine m_ES_wd_eko2(nfout,mode)
    integer, intent(in) ::                  nfout
    integer, intent(in) ::                  mode
    real(kind=DP), pointer, dimension(:) :: eko, eko_t ! d(neg) MPI

    integer :: ik, ib, ikp

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    allocate(eko(neg)); allocate(eko_t(neg))           ! MPI
    write(nfout,*) '=== energy_eigen_values ==='

    ikp = 0

! =========================== modified by K. Tagami =================== 12.0Exp
!    if(mode == EK) ikp = nk_in_the_process - 1
    if (mode == EK) then
       if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
          ikp = nk_in_the_process - 1
       endif
    endif
! ======================================================================= 12.0Exp

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
    do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle                 ! MPI
       eko_t = 0                                       ! MPI
       do ib = 1, neg                                  ! MPI
          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik) ! MPI
       end do                                          ! MPI
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
       if(mode == EK .or. ik == 1) then
          write(nfout,'(" ik = ",i7," ( ",3f10.6," )",/99(4f18.10,/))')&
               &ik+ikp,vkxyz(ik,1:3,BUCS), (eko(neordr(ib,ik)),ib=1,neg-num_extra_bands)
!!$          write(nfout,'(" ik = ",i4,99(4f18.10,/))')&
!!$               &ik,(eko(neordr(ib,ik)),ib=1,neg)
          if(iprieigenvalue >= 2 .and. num_extra_bands >= 1) &
               & write(nfout,'(" -- extra_bands --",/99(4d20.8,/))') &
               & (eko(neordr(ib,ik)),ib=neg-num_extra_bands+1, neg)
       else
          write(nfout,'(" ik = ",i7," (",3f10.6," )",/99(5f20.8,/))')&
               &ik,vkxyz(ik,1:3,BUCS), (eko(neordr(ib,ik)),ib=1,neg)
       endif
    end do
    deallocate(eko); deallocate(eko_t)                 ! MPI
  end subroutine m_ES_wd_eko2


  subroutine m_ES_sum_of_LocalPart(ik,ibo,bpr_l,bpi_l,dz)
    integer, intent(in)                                 :: ik,ibo
    real(kind=DP),intent(in),dimension(np_e,np_fs,ik:ik):: bpr_l,bpi_l !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q

    ib = ibo
    dz = 0.d0
    do ia = 1, nac_p
       p = nlmta1_p(ia); q = nlmta2_p(ia)
       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr_l(ib,q,ik) + bpi_l(ib,p,ik)*bpi_l(ib,q,ik))
    end do
  end subroutine m_ES_sum_of_LocalPart


!!$  subroutine m_ES_sum_of_LocalPart_3D(ik,ib,bpr_l,bpi_l,dz)
!!$
!!$    integer, intent(in)                                 :: ik,ib
!!$    real(kind=DP),intent(in),dimension(np_e,nlmta,ik:ik):: bpr_l,bpi_l !MPI
!!$    real(kind=DP),intent(out)                           :: dz
!!$
!!$    integer  :: ia, p, q
!!$
!!$    dz = 0.d0
!!$    do ia = 1, nac
!!$       p = nlmta1(ia); q = nlmta2(ia)
!!$       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr_l(ib,q,ik) + bpi_l(ib,p,ik)*bpi_l(ib,q,ik))
!!$    end do
!!$  end subroutine m_ES_sum_of_LocalPart_3D

  subroutine m_ES_sum_of_LocalPart2(ik,ibo,bpr_l,bpi_l,bpr_np,bpi_np,dz)
    integer, intent(in)                                 :: ik,ibo
    real(kind=DP),intent(in),dimension(np_e,nlmta,ik:ik):: bpr_l,bpi_l,bpr_np,bpi_np !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q

    ib = map_z(ibo)
    dz = 0.d0
    do ia = 1, nac
       p = nlmta1(ia); q = nlmta2(ia)
       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr_np(ib,q,ik) + bpi_l(ib,p,ik)*bpi_np(ib,q,ik))
    end do
  end subroutine m_ES_sum_of_LocalPart2

! ============================== added by K. Tagami ======================== 11.0
  subroutine m_ES_sum_of_LocalPart2_noncl( ik, ibo, kst, ken, &
       &                                   bpr_l, bpi_l, &
       &                                   bpr_np, bpi_np, dz )
    integer, intent(in)       :: ik,ibo
    integer, intent(in)       :: kst, ken
    real(kind=DP), intent(in) :: bpr_l(np_e,nlmta,kst:ken)
    real(kind=DP), intent(in) :: bpi_l(np_e,nlmta,kst:ken)
    real(kind=DP), intent(in) :: bpr_np(np_e,nlmta,kst:ken)
    real(kind=DP), intent(in) :: bpi_np(np_e,nlmta,kst:ken)

    real(kind=DP),intent(out)      :: dz

    integer  :: ib, ia, p, q
    integer :: k1, k2
    integer :: is1, is2, is_tmp
    real(kind=DP) :: c1, c2

    ib = map_z(ibo)
    dz = 0.d0
    do ia = 1, nac
       p = nlmta1(ia); q = nlmta2(ia)
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             is_tmp = (is1 -1 )*ndim_spinor + is2
             k1 = ik + is1 -1
             k2 = ik + is2 -1

             c1 = real( fqwei_noncl(ia,is_tmp) ) &
                  &   *( bpr_l(ib,p,k1) *bpr_np(ib,q,k2) &
                  &    + bpi_l(ib,p,k1) *bpi_np(ib,q,k2))
             c2 =-aimag(fqwei_noncl(ia,is_tmp) ) &
                  &   *( bpr_l(ib,p,k1) *bpi_np(ib,q,k2) &
                  &     -bpi_l(ib,p,k1) *bpr_np(ib,q,k2))

             dz = dz + c1 + c2
          Enddo
       End do
    end do
  end subroutine m_ES_sum_of_LocalPart2_noncl
! ==================================================================== 11.0


  subroutine m_ES_sum_of_LocalPart3(ik,ibo,bpr_l,bpi_l,bpr1_l,bpi1_l,dz)
    integer, intent(in)                                 :: ik,ibo
    real(kind=DP),intent(in),dimension(np_e,np_fs,ik:ik):: bpr_l,bpi_l !MPI
    real(kind=DP),intent(in),dimension(np_e,np_fs,ista_k:iend_k):: bpr1_l,bpi1_l !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q

    ib = ibo
    dz = 0.d0
    do ia = 1, nac_p
       p = nlmta1_p(ia); q = nlmta2_p(ia)
       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr1_l(ib,q,ik) + bpi_l(ib,p,ik)*bpi1_l(ib,q,ik))
    end do
  end subroutine m_ES_sum_of_LocalPart3

! ============================== added by K. Tagami ======================== 11.0
  subroutine m_ES_sum_of_LocalPart3_noncl( ik, ibo, kst, ken, &
       &                                   bpr_l, bpi_l, &
       &                                   bpr1_l, bpi1_l, dz )
    integer, intent(in) :: ik,ibo
    integer, intent(in) :: kst, ken
    real(kind=DP),intent(in),dimension(np_e,nlmta,kst:ken):: bpr_l,bpi_l !MPI
    real(kind=DP),intent(in),dimension(np_e,nlmta,ista_k:iend_k):: bpr1_l,bpi1_l !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q
    integer :: k1, k2
    integer :: is1, is2, is_tmp
    real(kind=DP) :: c1, c2

    ib = map_z(ibo)
    dz = 0.d0
    do ia = 1, nac
       p = nlmta1(ia); q = nlmta2(ia)
       Do is1=1, ndim_spinor
          Do is2=1, ndim_spinor
             is_tmp = (is1 -1 )*ndim_spinor + is2
             k1 = ik + is1 -1
             k2 = ik + is2 -1

             c1 = real( fqwei_noncl(ia,is_tmp) ) &
                  &   *( bpr_l(ib,p,k1) *bpr1_l(ib,q,k2) &
                  &    + bpi_l(ib,p,k1) *bpi1_l(ib,q,k2))
             c2 =-aimag(fqwei_noncl(ia,is_tmp) ) &
                  &   *( bpr_l(ib,p,k1) *bpi1_l(ib,q,k2) &
                  &     -bpi_l(ib,p,k1) *bpr1_l(ib,q,k2))

             dz = dz + c1 + c2
          Enddo
       End do
    end do
  end subroutine m_ES_sum_of_LocalPart3_noncl
! =========================================================================== 11.0

  subroutine m_ES_cpeko()
    eko1_l = eko_l
  end subroutine m_ES_cpeko

  logical function m_ES_eekdif()
    integer       :: ik, ib
    real(kind=DP) :: fac, tmp, ekosum

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    if(ekmode == ON .and. evaluation_eko_diff == OFF) then
       m_ES_eekdif = .false.
       return
    end if

    evdff = 0.d0; evdffr = 0.d0

! =================================== modified by K. Tagami ============ 11.0
!!    fac = 1.d0/dble(kv3*neg)
!
!    if ( noncol ) then
!      fac = 1.0d0 / dble( kv3/ndim_spinor *neg )
!    else
!      fac = 1.d0/dble(kv3*neg)
!    endif
    fac = 1.d0/dble(neg)
!    endif
! ====================================================================== 11.0

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

    ekosum = 0.d0
! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
    do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle
! === For epsmain by tkato 2013/11/14 ==========================================
       do ib = 1, np_e
          if(nrvf_ordr(neg_g(ib),ik) > neg - num_extra_bands) cycle
          tmp = eko_l(ib,ik) - eko1_l(ib,ik)
! ==============================================================================
          evdff(1) = evdff(1) + tmp*tmp
          evdff(2) = evdff(2) + dabs(tmp)
          ekosum = ekosum + eko_l(ib,ik)
          if(dabs(tmp).gt.DELTAevdff) then
! === For epsmain by tkato 2013/11/14 ==========================================
             evdff(3) = evdff(3) + dsqrt(dabs(eko1_l(ib,ik)**2 - eko_l(ib,ik)**2))
! ==============================================================================
          end if
       end do
    end do

    if(npes > 1) then
! === For epsmain by tkato 2013/11/14 ==========================================
       call mpi_allreduce(MPI_IN_PLACE,evdff,3 &
         & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)  ! MPI
!       call mpi_allreduce(evdff,evdffr,3 &
!         & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)  ! MPI
       call mpi_allreduce(evdff(1),evdffr(1),1 &
         & ,mpi_double_precision,mpi_max,mpi_ge_world,ierr)  ! MPI
       call mpi_allreduce(evdff(2),evdffr(2),1 &
         & ,mpi_double_precision,mpi_max,mpi_ge_world,ierr)  ! MPI
       call mpi_allreduce(evdff(3),evdffr(3),1 &
         & ,mpi_double_precision,mpi_max,mpi_ge_world,ierr)  ! MPI
       call mpi_allreduce(MPI_IN_PLACE,ekosum,1,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)
!       call mpi_allreduce(MPI_IN_PLACE,ekosum,1,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)
! ==============================================================================
       evdff = evdffr
    end if

!!$    if(iprievdff >= 2) write(6,'(" ! evdff = ",3d16.8)') evdff(1:3)
    evdff(1) = dsqrt(fac*evdff(1))*(af+1)
    evdff(2) = fac*evdff(2)*(af+1)
    evdff(3) = fac*evdff(3)*(af+1)

    eko1_l = eko_l

    if(iprievdff >= 1) &
!!$         & write(nfout,'(" >> (",i6,") <eko_old-eko_new>:(",3d13.5,")")') &
!!$         &    iteration_electronic, evdff(1), evdff(2), evdff(3)
         & write(nfout,'(" ENERGY EIGEN VALUE SUM. FOR",i6," -TH ITER (K=",i6,":",i6 &
         & ,") = ",F16.8, " <eko_old-eko_new>:(",3d13.5,")")') &
         &   iteration_electronic,nk_in_the_process, nk_in_the_process+kv3-1 &
         & , ekosum, evdff(1), evdff(2), evdff(3)

    if(max(evdff(1),evdff(2)) < delta_eigenvalue) then
       ib = m_CtrlP_ntcnvg_incre()
       if(printable) then
          write(nfout,'(" !iter = ",i7," ntcnvg = ",i7)') iteration_electronic, ib
          if(iprievdff < 1) &
               & write(nfout,'(" >> (",i6,") <eko_old-eko_new>:(",3d13.5,")")') &
               &    iteration_electronic, evdff(1), evdff(2), evdff(3)
       end if
    else
       call m_CtrlP_ntcnvg_reset()    ! K.Mae 030808
    end if
    m_ES_eekdif = m_CtrlP_ntcnvg_clear()

  end function m_ES_eekdif

  logical function m_ES_eekdif2()
    integer       :: ik, ib, it, kv3_e
    real(kind=DP), allocatable, dimension(:,:) :: evdff2, evdffr2 ! d(3,kv3)
    integer, allocatable, dimension(:)         :: iconv           ! d(kv3)
    real(kind=DP) :: fac, tmp

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    if(evaluation_eko_diff == OFF) then
       m_ES_eekdif2 = .false.
       return
    end if

    allocate(evdff2(3,kv3));  evdff2  = 0.d0
    allocate(evdffr2(3,kv3)); evdffr2 = 0.d0
    allocate(iconv(kv3)); iconv = 0

    fac = 1.d0/dble(neg)

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
    do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle
       do ib = ista_e, iend_e, istep_e
          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          tmp = eko_l(map_z(ib),ik) - eko1_l(map_z(ib),ik)
          evdff2(1,ik) = evdff2(1,ik) + tmp*tmp
          evdff2(2,ik) = evdff2(2,ik) + dabs(tmp)
          if(dabs(tmp).gt.DELTAevdff) then
             evdff2(3,ik) = evdff2(3,ik)  &
                  & + dsqrt(dabs(eko1_l(map_z(ib),ik)**2 - eko_l(map_z(ib),ik)**2))
          end if
       end do
    end do
    if(npes > 1) then
       call mpi_allreduce(mpi_in_place,evdff2,3*kv3 &
            & ,mpi_double_precision,mpi_sum,mpi_kg_world,ierr)  ! MPI
       call mpi_allreduce(mpi_in_place,evdff2,3*kv3 &
            & ,mpi_double_precision,mpi_sum,mpi_ge_world,ierr)  ! MPI
    end if

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
    do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       evdff2(1,ik) = dsqrt(fac*evdff2(1,ik))*(af+1)
       evdff2(2,ik) = fac*evdff2(2,ik)*(af+1)
       evdff2(3,ik) = fac*evdff2(3,ik)*(af+1)
    end do

    eko1_l = eko_l

    if(iprievdff >= 1) then
! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
      do ik = 1, kv3, ikskip
! ====================================================================== 11.0
          write(nfout,'(" >> (",i6,") <eko_old-eko_new>(ik=",i6,") :(",3d13.5 &
               & ,")")') iteration_electronic,ik &
               & ,evdff2(1,ik),evdff2(2,ik),evdff2(3,ik)
       end do
    end if

    it = 0
! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
      do ik = 1, kv3, ikskip
! ====================================================================== 11.0
       if(max(evdff2(1,ik),evdff2(2,ik)) < delta_eigenvalue) then
          iconv(ik) = EK_CONVERGED
       else
          it = it + 1
       end if
    end do
!!!!!BRANCH_P 3D_Parallel
!!!!    kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1))
!!!!!BRANCH_P_END 3D_Parallel
!!!!!BRANCH_P ORG_Parallel
! === DEBUG for icond=2,3 & one_by_one by tkato 2014/01/24 =====================
    if(first_kpoint_in_this_job == 0) then
       kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1))
    else
       kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1)-first_kpoint_in_this_job+1)
    end if
! ==============================================================================
!!!!!BRANCH_P_END ORG_Parallel
    if(iprievdff >= 2) then
       write(nfout,'(" ! -- iconv_ek(before) -- <<m_ES_eekdif2>>")')
       write(nfout,'(" ! ",10i8)') iconv_ek(1:kv3_ek)
    end if
    do ik = 1, kv3_e
       if(iconv_ek(nk_in_the_process+ik-1) /= EK_CONVERGED) then
          iconv_ek(nk_in_the_process+ik-1) = iconv(ik)
       end if
    end do
!!$    iconv_ek(nk_in_the_process:nk_in_the_process+kv3_e-1) = iconv(1:kv3_e)
    if(iprievdff >= 2) then
       write(nfout,'(" ! -- iconv_ek(after) -- <<m_ES_eekdif2>>")')
       write(nfout,'(" ! ",10i8)') iconv_ek(1:kv3_ek)
    end if

    if(it == 0) then
       ib = m_CtrlP_ntcnvg_incre()
       if(printable) then
          write(nfout,'(" !iter = ",i7," ntcnvg = ",i7)') iteration_electronic, ib
          if(iprievdff < 1) then
! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
             do ik = 1, kv3, ikskip
! ====================================================================== 11.0

                write(nfout,'(" >> (",i5,") <eko_(old-new)>(ik=",i4,")=("&
                     & ,3d12.4 ,", iconv = ",i2,")")') iteration_electronic,ik&
                     & ,evdff2(1,ik),evdff2(2,ik),evdff2(3,ik),iconv(ik)
             end do
          end if
       end if
    else
       call m_CtrlP_ntcnvg_reset()    ! K.Mae 030808
    end if
    m_ES_eekdif2 = m_CtrlP_ntcnvg_clear()
    deallocate(evdff2,evdffr2)
    deallocate(iconv)
  end function m_ES_eekdif2

  logical function m_ES_eekdif_cond()
    integer       :: ik, ib, ibcm, ipri0
    real(kind=DP) :: fac, tmp

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

!!$    if(evaluation_eko_diff == OFF) then
    if(.not.delta_eigenvalue_cond_is_given) then
       m_ES_eekdif_cond = .true.
       return
    end if

    ibcm = totch/2 + 1.1
    if(printable) write(nfout,'("!m_ES_eekdif_cond: ibcm = ",i8)') ibcm

    evdff = 0.d0; evdffr = 0.d0
    if(neg-ibcm > 0) then
       fac = 1.d0/dble(kv3*(neg-ibcm))
    else
       fac = 1.d0/dble(kv3*neg)
    end if

! ======================= added by K. Tagami ========================= 11.0
    if ( noncol ) fac = fac * dble(ndim_spinor)
! =================================================================== 11.0

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
      do ik = 1, kv3, ikskip
! ====================================================================== 11.0

       if(map_k(ik) /= myrank_k) cycle
       do ib = ista_e, iend_e, istep_e
          if(nrvf_ordr(ib,ik) < ibcm .or. nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          tmp = eko_l(map_z(ib),ik) - eko1_l(map_z(ib),ik)
          evdff(1) = evdff(1) + tmp*tmp
          evdff(2) = evdff(2) + dabs(tmp)
          if(dabs(tmp).gt.DELTAevdff) then
             evdff(3) = evdff(3) + dsqrt(dabs(eko1_l(map_z(ib),ik)**2 - eko_l(map_z(ib),ik)**2))
          end if
       end do
    end do

    if(npes > 1) then
       call mpi_allreduce(evdff,evdffr,3 &
         & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)  ! MPI
       evdff = evdffr
    end if

    evdff(1) = dsqrt(fac*evdff(1))*(af+1)
    evdff(2) = fac*evdff(2)*(af+1)
    evdff(3) = fac*evdff(3)*(af+1)

    eko1_l = eko_l


    if(iprievdff >= 1) &
         & write(nfout,'(" >> (",i6,") <eko_old-eko_new>:(",3d13.5,")")') &
         &    iteration_electronic, evdff(1), evdff(2), evdff(3)

    if(max(evdff(1),evdff(2)) < delta_eigenvalue_conduction) then
       m_ES_eekdif_cond = .true.
       if(printable) then
          if(iprievdff < 1) &
               & write(nfout,'(" >> (",i6,") <eko_old-eko_new>:(",3d13.5,") : m_ES_eekdif_cond")') &
               &    iteration_electronic, evdff(1), evdff(2), evdff(3)
       end if
    else
       m_ES_eekdif_cond = .false.
    end if

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 1 ) call m_ES_wd_eko_cond(nfout,ibcm,mode=SCF)

  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0

  end function m_ES_eekdif_cond

  subroutine m_ES_cp_eko_l_to_eko_ek()
    real(kind=DP),pointer,dimension(:,:) :: eko_t, eko_t2 ! d(neg,kv3)
    integer :: ik, ib, kv3_e

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    allocate(eko_t(neg,kv3)); allocate(eko_t2(neg,kv3))
    eko_t = 0.d0 ; eko_t2 = 0.d0

! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = af +1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3, af+1
      do ik = 1, kv3, ikskip
! ====================================================================== 11.0
       if(map_k(ik) /= myrank_k) cycle
       do ib = 1, neg
          if(map_e(ib) == myrank_e) eko_t(ib,ik) = eko_l(map_z(ib),ik)
       end do
    end do
    call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
         & ,mpi_sum,mpi_kg_world,ierr)
    kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1))
    if(nk_in_the_process+kv3_e-1>kv3_ek) kv3_e = kv3_ek+1-nk_in_the_process
    eko_ek(:,nk_in_the_process:nk_in_the_process+kv3_e-1) = eko_t2(:,1:kv3_e)
    deallocate(eko_t); deallocate(eko_t2)
  end subroutine m_ES_cp_eko_l_to_eko_ek

  subroutine m_ES_cp_eko_l_to_eko_ek2()
    real(kind=DP),allocatable,dimension(:,:) :: eko_t, eko_t2 ! d(neg,kv3)
    integer :: ik, ib, kv3_e

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0
    integer :: ikt,is, kv3t

    allocate(eko_t(neg,kv3)); allocate(eko_t2(neg,kv3))

    eko_t = 0.d0
! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = 1
    endif
! ====================================================================== 11.0

    do is=ista_spin, iend_spin, ikskip
! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3
!      do ik = 1, kv3, ikskip
! ====================================================================== 11.0
      do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle
       do ib = 1, neg
          if(map_e(ib) /= myrank_e) cycle
          eko_t(ib,ik) = eko_l(map_z(ib),ik)
       end do
      end do
    end do
    if(npes >= 2) then
       call mpi_allreduce(MPI_IN_PLACE,eko_t,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, mpi_kg_world,ierr)
       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, mpi_ge_world,ierr)
    else
       eko_t2 = eko_t
    end if

    kv3_e = min(kv3,kv3_ek-kv3*(nkgroup-1))
    if(nk_in_the_process+kv3_e-1 > kv3_ek) kv3_e = kv3_ek - nk_in_the_process + 1
    if(kv3_e <=0) kv3_e = 1

    if(sw_modified_kpoint_increment == ON)then
      do ik=0,nrank_k-1
         do is=1,nspin
            ikt = nspin*(nis_kv3_ek(ik)-1)+(nkgroup-1)*nspin+is
            if(ikt>kv3_ek) ikt = kv3_ek
            !eko_ek(:,ikt) = eko_t2(:,(ik-1)*nspin+is)
            eko_ek(:,ikt) = eko_t2(:,ik*nspin+is)
         enddo
      enddo
    else
      eko_ek(:,nk_in_the_process:nk_in_the_process+kv3_e-1) = eko_t2(:,1:kv3_e)
    endif

    deallocate(eko_t); deallocate(eko_t2)
  end subroutine m_ES_cp_eko_l_to_eko_ek2

  subroutine m_ES_cp_eko_l_to_eko_ek0()
    real(kind=DP),allocatable,dimension(:,:) :: eko_t, eko_t2 ! d(neg,kv3)
    integer :: ik, ib, kv3_e

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    allocate(eko_t(neg,kv3)); allocate(eko_t2(neg,kv3))

    eko_t = 0.d0
! ==================================== added by K. Tagami =============== 11.0
    if ( noncol ) then
      ikskip = ndim_spinor
    else
      ikskip = 1
    endif
! ====================================================================== 11.0

! ================================== modified by K. Tagami ============= 11.0
!!!    do ik = 1, kv3
      do ik = 1, kv3, ikskip
! ====================================================================== 11.0
       if(map_k(ik) /= myrank_k) cycle
       do ib = 1, neg
          if(map_e(ib) /= myrank_e) cycle
          eko_t(ib,ik) = eko_l(map_z(ib),ik)
       end do
    end do
    if(npes >= 2) then
       call mpi_allreduce(MPI_IN_PLACE,eko_t,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, mpi_kg_world,ierr)
       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, mpi_ge_world,ierr)
    else
       eko_t2 = eko_t
    end if

    eko_ek(:,:) = eko_t2(:,:)

    deallocate(eko_t); deallocate(eko_t2)
  end subroutine m_ES_cp_eko_l_to_eko_ek0

  real(DP) function m_ES_what_is_evdff_now()
    m_ES_what_is_evdff_now = max(evdff(1),evdff(2))
  end function m_ES_what_is_evdff_now

  real(kind=DP) function m_ES_get_energy(ilevel)
    integer, intent(in) :: ilevel
    real(kind=DP), pointer, dimension(:) :: eko, eko_t

    integer :: ib
    allocate(eko(neg)); allocate(eko_t(neg))

    do ib = 1, neg
       if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),1)
    end do
    call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI

    m_ES_get_energy = eko(neordr(ilevel,1))
  end function m_ES_get_energy

  subroutine m_ES_orbital_population(nfout)
    integer, intent(in) :: nfout

    integer :: is,ik,ib,iopr,i,it,ia,il,im,tau,lmt,ilmta,iksnl, lmtt
    real(kind=DP) :: temp, pup, pdn, fac
    real(kind=DP), dimension(nlmta_phi,2) :: porb_mpi

    if(.not.allocated(porb)) allocate(porb(nlmta_phi,2))

    porb_mpi = 0.d0
    porb   = 0.d0

    is = 1
    do ik = 1, kv3
       if(map_k(ik) /= myrank_k) cycle
       if(nspin==2) is = mod(ik-1,2)+1
       iksnl = (ik-1)/nspin + 1
       do ib = 1, neg
          if(map_e(ib) == myrank_e) then
             do i=1,nlmta_phi
                call m_PP_tell_iorb_lmtt(i,lmtt)
                temp = 0.d0
                if(k_symmetry(ik) == GAMMA) then
                   do iopr = 1,nopr
! =============================== Modified by K. Tagami =====  0.2b
!                      temp = temp + (  compr_l(map_z(ib),i,iopr,ik)**2) &
!                           &      *(1.d0+qorb(i)/norm_phig(lmt,iksnl))
                      temp = temp + compr_l(map_z(ib),i,iopr,ik)**2 /2.0 &
                           &      *(1.d0+qorb(i)/(norm_phig(lmtt,iksnl)*2.))
! ==========================================================
                   end do
                else
                   do iopr = 1,nopr
                      temp = temp + (  compr_l(map_z(ib),i,iopr,ik)**2 &
                           &         + compi_l(map_z(ib),i,iopr,ik)**2) &
                           &      *(1.d0+qorb(i)/norm_phig(lmtt,iksnl))
                   end do
                end if
                porb_mpi(i,is) = porb_mpi(i,is) + occup_l(map_z(ib),ik)*temp/dble(nopr)
             end do
          end if
       end do
    end do
    call mpi_allreduce(porb_mpi,porb,nlmta_phi*2,mpi_double_precision &
         &            ,mpi_sum,mpi_kg_world,ierr)
    porb_mpi = porb
    call mpi_allreduce(porb_mpi,porb,nlmta_phi*2,mpi_double_precision &
         &            ,mpi_sum,mpi_ge_world,ierr)

    fac = 2.d0/kv3
    do is=1,nspin
       do i=1,nlmta_phi
          porb(i,is) = fac*porb(i,is)
       end do
    end do
    if(nspin == 1) then
       do i=1,nlmta_phi
          porb(i,1)=0.5d0*porb(i,1)
          porb(i,2)=porb(i,1)
       end do
    end if

    pup=0.d0; pdn=0.d0
    do i=1,nlmta_phi
       pup = pup +porb(i,1);    pdn = pdn +porb(i,2)
    end do

    ! Write down Porb
    if(ipri >= 1) then  ! -->

       write(nfout,'(" --------- Orbital population ---------")')
       write(nfout,'("  ia   l   m   t    Porb(UP)   Porb(DN) element")')
       do ia=1,natm
          if(if_pdos(ia) == OFF) cycle
! ======================================== Added by K. Tagami == 0.2a =
          if(iproj_group(ia) == 0) cycle
! ================================================================
          it = ityp(ia)
          do lmt=1,ilmt_phi(it)
             ilmta = lmta_phi(lmt,ia)
             il  = ltp_phi(lmt,it)
             im  = mtp_phi(lmt,it)
             tau = taup_phi(lmt,it)
             write(nfout,'(4(1x,i3),2(1x,f10.5),4x,4a)') &
               &    ia, il-1, im, tau, porb(ilmta,1:2), speciesname(it)
          end do
       end do
       write(nfout,'(8x,"Total : ",3(1x,f10.5))') pup,pdn,pup+pdn

    end if      ! <--

  end subroutine m_ES_orbital_population

! ============================== added by K. Tagami =================== 11.0
  subroutine m_ES_orbital_population_noncl(nfout)
    integer, intent(in) :: nfout

    integer :: ik,ib,iopr,i,it,ia,il,im,tau,lmt,ilmta,iksnl
    real(kind=DP) :: fac
    integer :: is1, is2, istmp, lmtt

    real(kind=DP), allocatable :: porb_mpi( :,: )
    real(kind=DP) :: porb_sum( ndim_magmom )
    complex(kind=CMPLDP) :: porb_ssrep( nlmta_phi, ndim_magmom )
    complex(kind=CMPLDP) :: z1, z2, ztemp

    if(.not.allocated(porb)) allocate(porb(nlmta_phi,ndim_magmom))

    porb_ssrep = 0.d0;      porb   = 0.d0

    do ik = 1, kv3, ndim_spinor

       if (map_k(ik) /= myrank_k) cycle
       iksnl = (ik-1)/ndim_spinor + 1

       do ib = 1, neg
          if(map_e(ib) == myrank_e) then
             do i=1,nlmta_phi
                call m_PP_tell_iorb_lmtt(i,lmtt)

                if(k_symmetry(ik) == GAMMA) then
                   call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in noncollinear system.', &
                   __LINE__,__FILE__)
                else
                   Do is1=1, ndim_spinor
                      Do is2=1, ndim_spinor
                         istmp = ( is1 -1 )*ndim_spinor + is2

                         ztemp = 0.0d0
                         do iopr = 1,nopr
                            z1 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is1-1 ), &
                                 &       compi_l(map_z(ib),i,iopr,ik+is1-1 ) )
                            z2 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is2-1 ), &
                                 &       compi_l(map_z(ib),i,iopr,ik+is2-1 ) )
                            ztemp = ztemp + z1 *conjg(z2) &
                                 &      *( 1.d0+qorb(i)/norm_phig(lmtt,iksnl) )
                         End do

                         porb_ssrep(i,istmp) = porb_ssrep(i,istmp) &
                              &  + occup_l(map_z(ib),ik) *ztemp /dble(nopr)

                      End do
                   End do

                end if

             end do
          end if
       end do
    end do
! --

    fac = 1.d0 / ( kv3 /ndim_spinor )
    porb_ssrep = porb_ssrep *fac

!------------
    if ( SpinOrbit_mode == BuiltIn ) then
       call phase_error_with_msg(nfout,'kt : Under construction porb ',__LINE__,__FILE__)
    endif
! -----

    call m_ES_DensMat_To_MagMom_porb( nlmta_phi, porb_ssrep, porb )

    if ( npes >=2 ) then
       allocate( porb_mpi(nlmta_phi,ndim_magmom ) )
       porb_mpi = 0.0d0
       call mpi_allreduce( porb, porb_mpi, nlmta_phi*ndim_magmom, &
            &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
       porb = porb_mpi
       deallocate( porb_mpi )
    endif

    porb_sum = 0.0d0
    do i=1,nlmta_phi
       porb_sum(:) = porb_sum(:) + porb(i,:)
    end do

! Write down Porb
    if (ipri >= 1) then

       write(nfout,'(" --------- Orbital population ---------")')
       write(nfout,'("  ia   l   m   t    Porb(tot)  Porb(mx)   Porb(my)   Porb(mz) element")')
       do ia=1,natm
          if(if_pdos(ia) == OFF) cycle

          if(iproj_group(ia) == 0) cycle

          it = ityp(ia)
          do lmt=1,ilmt_phi(it)
             ilmta = lmta_phi(lmt,ia)
             il  = ltp_phi(lmt,it)
             im  = mtp_phi(lmt,it)
             tau = taup_phi(lmt,it)
             write(nfout,'(4(1x,i3),4(1x,f10.5),4x,4a)') &
                  & ia,il,im,tau,porb(ilmta,1:ndim_magmom),speciesname(it)
          end do
       end do

       write(nfout,'(" -------- summation ---------")')
       write(nfout,'("             total       mx         my         mz")')
       write(nfout,'(8x,4(1x,f10.5))') porb_sum(1:ndim_magmom)

    end if

  end subroutine m_ES_orbital_population_noncl
! ===================================================================== 11.0

  subroutine m_ES_sym_comp(nfout)
    integer, intent(in) :: nfout

    integer :: ik,iopr,mm,jorb,iorb,ib,is

    do is=ista_spin, iend_spin
!    do ik=1,kv3
    do ik = is, kv3-nspin+is, nspin
       if(map_k(ik) /= myrank_k) cycle
       if(k_symmetry(ik) == GAMMA) then
          compr_l(:,:,2:nopr,ik) = 0.d0
          do iopr=2,nopr
             do iorb=1,nlmta_phi
                do mm=1,nrorb(iorb,iopr)
                   jorb=irorb(mm,iorb,iopr)
                   do ib=1,np_e
                      compr_l(ib,iorb,iopr,ik) = compr_l(ib,iorb,iopr,ik) &
                           & + compr_l(ib,jorb,1,ik)*crorb(mm,iorb,iopr)
                   end do
                end do
             end do
          end do
       else
          compr_l(:,:,2:nopr,ik) = 0.d0
          compi_l(:,:,2:nopr,ik) = 0.d0
          do iopr=2,nopr
             do iorb=1,nlmta_phi
                do mm=1,nrorb(iorb,iopr)
                   jorb=irorb(mm,iorb,iopr)
                   do ib=1,np_e
                      compr_l(ib,iorb,iopr,ik) = compr_l(ib,iorb,iopr,ik) &
                           & + compr_l(ib,jorb,1,ik)*crorb(mm,iorb,iopr)
                      compi_l(ib,iorb,iopr,ik) = compi_l(ib,iorb,iopr,ik) &
                           & + compi_l(ib,jorb,1,ik)*crorb(mm,iorb,iopr)
                   end do
                end do
             end do
          end do
       end if
    end do
    end do

    if(ipri>=2) write(nfout,*) '** compr and compi were symmetrized. **'

  end subroutine m_ES_sym_comp

  subroutine m_ES_orbital_den_mat(nfout,dm,max2lp,max_projs)
    integer, intent(in) :: nfout,max2lp,max_projs
    real(kind=DP), intent(out) :: dm(max2lp,max2lp,max_projs,natm,nspin)

    integer :: is,ik,ib,iopr,i,it,ia,im,lmt,iksnl
    integer :: im1,im2,i1,i2,lmt1,lmt2,ip
    integer :: lmtt, lmtt1, lmtt2

    real(kind=DP) :: temp, fac
    real(kind=DP), dimension(max2lp,max2lp,max_projs,natm,nspin) :: dm_mpi

    dm_mpi = 0.d0
    dm     = 0.d0

    is = 1
    do ik = 1, kv3
       if(map_k(ik) /= myrank_k) cycle
       if(nspin==2) is = mod(ik-1,2)+1
       iksnl = (ik-1)/nspin + 1
       do ib = 1, neg
          if(map_e(ib) == myrank_e) then
             do ia=1,natm
                it = ityp(ia)
                !!if(if_pdos(ia) == OFF) cycle
                if(iproj_group(ia) == 0) cycle

                ! diagonal part
                do lmt=1,ilmt_phi(it)
                   im  = mtp_phi(lmt,it)
                   ip  = iproj_phi(lmt,it)
                   i = lmta_phi(lmt,ia)
                   lmtt = lmtt_phi(lmt,it)

                   temp = 0.d0
                   if(k_symmetry(ik) == GAMMA) then
                      do iopr = 1,nopr
! ============================ Modified by K. Tagami ==== 0.2b
!                         temp = temp + (  compr_l(map_z(ib),i,iopr,ik)**2) &
!                             &      *(1.d0+qorb(i)/norm_phig(lmt,iksnl))
                         temp = temp + compr_l(map_z(ib),i,iopr,ik)**2 /2.0 &
                             &      *(1.d0+qorb(i)/(norm_phig(lmtt,iksnl)*2.0))
! ======================================================
                      end do
                   else
                      do iopr = 1,nopr
                         temp = temp + (  compr_l(map_z(ib),i,iopr,ik)**2 &
                              &         + compi_l(map_z(ib),i,iopr,ik)**2) &
                              &      *(1.d0+qorb(i)/norm_phig(lmtt,iksnl))
                      end do
                   end if
                   dm_mpi(im,im,ip,ia,is) = dm_mpi(im,im,ip,ia,is) &
                        &                  + occup_l(map_z(ib),ik)*temp/dble(nopr)
                end do

                ! non-diagonal part
                do lmt2=1,ilmt_phi(it)
                   do lmt1=1,lmt2-1
                      if(ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle
                      ip  = iproj_phi(lmt1,it)
                      im1  = mtp_phi(lmt1,it)
                      i1 = lmta_phi(lmt1,ia)
                      im2  = mtp_phi(lmt2,it)
                      i2 = lmta_phi(lmt2,ia)

                      lmtt1 = lmtt_phi(lmt1,it);     lmtt2 = lmtt_phi(lmt2,it)

                      temp = 0.d0
                      if(k_symmetry(ik) == GAMMA) then
                         do iopr = 1,nopr
! ===================== Modified by K. Tagami ==== ???? === 0.2b
!                            temp = temp + ( &
!                                 &   compr_l(map_z(ib),i1,iopr,ik)*compr_l(map_z(ib),i2,iopr,ik)) &
!                                 & /sqrt(norm_phig(lmt1,iksnl)*norm_phig(lmt2,iksnl))
                            temp = temp + ( &
                                 &   compr_l(map_z(ib),i1,iopr,ik)*compr_l(map_z(ib),i2,iopr,ik)) &
                                 & /sqrt(norm_phig(lmtt1,iksnl)*norm_phig(lmtt2,iksnl)) / 4.0
! ========================================================
                         end do
                      else
                         do iopr = 1,nopr
                            temp = temp + ( &
                                 &   compr_l(map_z(ib),i1,iopr,ik)*compr_l(map_z(ib),i2,iopr,ik) &
                                 & + compi_l(map_z(ib),i1,iopr,ik)*compi_l(map_z(ib),i2,iopr,ik)) &
                                 & /sqrt(norm_phig(lmtt1,iksnl)*norm_phig(lmtt2,iksnl))
                         end do
                      end if
                      dm_mpi(im1,im2,ip,ia,is) = dm_mpi(im1,im2,ip,ia,is) + occup_l(map_z(ib),ik)*temp/dble(nopr)
                      dm_mpi(im2,im1,ip,ia,is) = dm_mpi(im1,im2,ip,ia,is)
                   end do
                end do
             end do
          end if
       end do
    end do
    call mpi_allreduce(dm_mpi,dm,max2lp**2*max_projs*natm*nspin,mpi_double_precision &
                    & ,mpi_sum,MPI_CommGroup,ierr)

    fac = 2.d0/kv3
    dm = fac*dm
    if(nspin == 1) dm = 0.5d0*dm

  end subroutine m_ES_orbital_den_mat

! ================================= added by K. Tagami ================= 11.0
  subroutine m_ES_orbital_den_mat_noncl( nfout, dm, max2lp, max_projs )
!
    integer, intent(in) :: nfout,max2lp,max_projs
    real(kind=DP), intent(out) :: dm(max2lp,max2lp,max_projs,natm,ndim_magmom)

    integer :: ik,ib,iopr,i,it,ia,im,lmt,iksnl
    integer :: im1,im2,i1,i2,lmt1,lmt2,ip
    integer :: lmtt, lmtt1, lmtt2
    integer :: is1, is2, istmp

    real(kind=DP) :: fac
    complex(kind=CMPLDP) :: z1, z2, ztemp
    complex(kind=CMPLDP), allocatable :: dm_ssrep( :, :, :, :, : )

    real(kind=DP), allocatable :: dm_mpi( :, :, :, :, : )

    allocate( dm_ssrep( max2lp, max2lp, max_projs, natm, ndim_chgpot ) )
    dm_ssrep = 0.d0; dm  = 0.d0

    do ik = 1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle

       iksnl = (ik-1)/ndim_spinor + 1

       do ib = 1, neg
          if(map_e(ib) == myrank_e) then
             do ia=1,natm
                it = ityp(ia)
                !!if(if_pdos(ia) == OFF) cycle
                if(iproj_group(ia) == 0) cycle

! ------------------------------- diagonal part ----------
                do lmt=1,ilmt_phi(it)
                   im  = mtp_phi(lmt,it)
                   ip  = iproj_phi(lmt,it)
                   i = lmta_phi(lmt,ia)
                   lmtt = lmtt_phi(lmt,it)

                   if (k_symmetry(ik) == GAMMA) then
                      call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in noncollinear system.', &
                      __LINE__,__FILE__)
                   else
                      Do is1=1, ndim_spinor
                         Do is2=1, ndim_spinor
                            istmp = ( is1 -1 )*ndim_spinor + is2

                            ztemp = 0.0d0
                            do iopr = 1,nopr
                               z1 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is1-1 ), &
                                    &       compi_l(map_z(ib),i,iopr,ik+is1-1 ) )
                               z2 = dcmplx( compr_l(map_z(ib),i,iopr,ik+is2-1 ), &
                                    &       compi_l(map_z(ib),i,iopr,ik+is2-1 ) )
                               ztemp = ztemp + z1 *conjg(z2) &
                                    &      *( 1.d0+qorb(i)/norm_phig(lmtt,iksnl) )
                            End do

                            dm_ssrep( im, im, ip, ia, istmp ) &
                                 & = dm_ssrep( im, im, ip, ia, istmp ) &
                                 &  + occup_l(map_z(ib),ik) *ztemp /dble(nopr)
                         End do
                      End do
                   end if

                end do

! --------------------  non-diagonal part --------------
                do lmt2 = 1,ilmt_phi(it)
!!!!                   do lmt1=1,lmt2-1

                   do lmt1 = 1, ilmt_phi(it)

                      if ( lmt1 == lmt2 ) cycle
                      if (ltp_phi(lmt1,it) /= ltp_phi(lmt2,it)) cycle

                      ip  = iproj_phi(lmt1,it)
                      im1  = mtp_phi(lmt1,it)
                      i1 = lmta_phi(lmt1,ia)
                      im2  = mtp_phi(lmt2,it)
                      i2 = lmta_phi(lmt2,ia)

                      lmtt1 = lmtt_phi(lmt1,it);   lmtt2 = lmtt_phi(lmt2,it)

                      if (k_symmetry(ik) == GAMMA) then
                         call phase_error_with_msg(nfout,'Not supported : Gamma symmetry in noncollinear system.', &
                         __LINE__,__FILE__)
                      else
                         Do is1=1, ndim_spinor
                            Do is2=1, ndim_spinor
                               istmp = ( is1 -1 )*ndim_spinor + is2

                               ztemp = 0.0d0
                               do iopr = 1,nopr
                                  z1 = dcmplx( compr_l(map_z(ib),i1,iopr,ik+is1-1 ), &
                                       &       compi_l(map_z(ib),i1,iopr,ik+is1-1 ) )
                                  z2 = dcmplx( compr_l(map_z(ib),i2,iopr,ik+is2-1 ), &
                                       &       compi_l(map_z(ib),i2,iopr,ik+is2-1 ) )
                                  ztemp = ztemp + z1 *conjg(z2) &
                                       & / sqrt( norm_phig(lmtt1,iksnl) &
                                       &         *norm_phig(lmtt2,iksnl) )
                               End do

                               dm_ssrep( im1, im2, ip, ia, istmp ) &
                                    & = dm_ssrep( im1, im2, ip, ia, istmp ) &
                                    &  + occup_l(map_z(ib),ik) *ztemp /dble(nopr)
                            End do
                         End do
                      end if
                   end do
                end do
! ------------------------------------------------------
             end do
          end if
       end do
    end do
! ----------
    fac = 1.d0 / ( kv3 /ndim_spinor )
    dm_ssrep = dm_ssrep *fac

!------------
    if ( SpinOrbit_mode == BuiltIn ) then
       call phase_error_with_msg(nfout,'kt : Under construction porb ',__LINE__,__FILE__)
    endif
! -----

    call m_ES_DensMat_To_MagMom_dm( dm_ssrep, dm, max2lp, max_projs )

    if ( npes > 1 ) then
       allocate( dm_mpi( max2lp, max2lp, max_projs, natm, ndim_magmom ) )
       dm_mpi = 0.0d0

       call mpi_allreduce( dm, dm_mpi, max2lp**2*max_projs*natm*ndim_magmom, &
         &                 mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       dm = dm_mpi
       deallocate( dm_mpi )
    endif
! --

  end subroutine m_ES_orbital_den_mat_noncl
! ================================================================ 11.0

  subroutine m_ES_set_num_bands_super()
    integer :: i
    i = neg*nlpnt
    call m_CntrlP_set_neg(i,nfout)
    call m_CntrlP_set_meg(i)
    if(printable) then
       write(nfout,*) 'num_bands will be changed.'
       write(nfout,'("neg,meg=",2i7)') neg,meg
    end if
  end subroutine m_ES_set_num_bands_super

  subroutine m_ES_add_neg(ne)
     integer, intent(in) :: ne
     integer :: newneg
     if(ne==0) return
     newneg = neg+ne
     if(printable) write(nfout,'(a,i6,a,i6)') ' !** the number of bands will change from ',neg,' to ',newneg
     call m_CntrlP_set_neg(newneg,neg)
     call m_CntrlP_set_meg(newneg)
     neg_previous = newneg
     !call m_Parallel_dealloc_mpi_elec()
     !call m_Parallel_init_mpi_elec(nfout,0,printable,neg,kv3,nspin,kg1)
  end subroutine m_ES_add_neg

  subroutine m_ES_cp_zaj_prev_to_zaj()
    if(mype == 0) then
      write(nfout,'(a)') ' !** copied previous wavefunctions to current wavefunctions'
    endif
    zaj_l = zaj_l_prev
  end subroutine m_ES_cp_zaj_prev_to_zaj

  subroutine m_ES_cp_zaj_to_zaj_prev()
    if(mype == 0) then
      write(nfout,'(a)') ' !** copied current wavefunctions to previous wavefunctions'
    endif
    zaj_l_prev = zaj_l
  end subroutine m_ES_cp_zaj_to_zaj_prev

  subroutine map_zaj_to_fft_box(ik,kg1,nfft,psi_l,bfft)
    integer, intent(in):: ik,kg1,nfft
    real(kind=DP), dimension(kg1,kimg), intent(in) :: psi_l
    real(kind=DP), dimension(nfft), intent(out) :: bfft
    integer :: i,i1,ri, j, i2, ii
    integer :: id_sname2=-1
    call tstatc0_begin('map_zaj_to_fft_box ',id_sname2)

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf(1)
          bfft(i1) = psi_l(1,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = igf(i)
             bfft(i1) = psi_l(ii,1)
             j = nbase_gamma(ii,2)
             i2 = igf(j)
             bfft(i2) =   psi_l(ii,1)
          end do
       else if(kimg == 2) then
          i1 = 2*igf(1) - 1
          bfft(i1)   = psi_l(1,1)
          bfft(i1+1) = psi_l(1,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = 2*igf(i)-1
             bfft(i1  ) = psi_l(ii,1)
             bfft(i1+1) = psi_l(ii,2)
             j = nbase_gamma(ii,2)
             i2 = 2*igf(j)-1
             bfft(i2  ) = psi_l(ii,1)
             bfft(i2+1) = -psi_l(ii,2)
          end do
       end if
    else
#ifdef NEC_TUNE_SMP
!CDIR NOLOOPCHG
#endif
       do ri = 1, kimg
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do i = 1, iba(ik)
             i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
             bfft(i1) = psi_l(i,ri)   ! MPI
          end do
       end do
    end if
  end subroutine map_zaj_to_fft_box

#ifdef FFTW3
  subroutine m_ES_gen_zaj_from_prev_zaj(nfout)
    integer, intent(in) :: nfout
    integer :: ib,ik,ri
    integer :: i,j,k, id, nl, nm, nn, nlhf,inew,jnew,knew,ip,mm,is,nl2,nm2,nn2,ig
    integer :: ii,jj,kk,ierr
    integer :: nfft_prev
    real(kind=DP), allocatable, dimension(:) :: xmat,ymat,zmat
    real(kind=DP), allocatable, dimension(:,:,:) :: cmatr, cmati, cmatr2, cmati2
    real(kind=DP) :: trilinear_interpolation
    integer(kind=DP) :: plan
    real(kind=DP), allocatable, dimension(:,:) :: zajtmp
    integer, parameter :: FFTW_MEASURE=0
    integer,save                   :: id_sname = -1
    call tstatc0_begin('m_ES_gen_zaj_from_prev_zaj ', id_sname, 1)
    allocate(bfft(nfft))
    id = fft_box_size_WF_prev(1,0)
    mm = fft_box_size_WF_prev(2,0)
    nl = fft_box_size_WF_prev(1,1)
    nm = fft_box_size_WF_prev(2,1)
    nn = fft_box_size_WF_prev(3,1)
    nl2=nl+2;nm2=nm+2;nn2=nn+2
    nfft_prev = product(fft_box_size_WF_prev(1:3,0)) * kimg
    allocate(afft(nfft_prev))
    allocate(xmat(nl2));allocate(ymat(nm2));allocate(zmat(nn2))
    allocate(cmatr(nl2,nm2,nn2))
    allocate(cmati(nl2,nm2,nn2))
    if(kimg == 1) then
       nlhf = id/2
    else
       nlhf = id
    end if

    if(kimg == 1)then
      call dfftw_plan_dft_r2c_3d(plan,nl,nm,nn,afft(1),afft(1),FFTW_MEASURE)
    else
      call dfftw_plan_dft_3d(plan,nl,nm,nn,afft(1),afft(1),+1,FFTW_MEASURE)
    endif

    zaj_l = 0.d0
    allocate(zajtmp(kg1_prev,kimg))
    do ik=1,kv3
       if(map_k(ik) /= myrank_k) cycle
!       do ib=ista_e,iend_e,istep_e
       do ib=1,np_e
          zajtmp = 0.d0
          do ig=ista_g1k_prev(ik),iend_g1k_prev(ik)
!             zajtmp(ig,:) = zaj_l_prev(ig-ista_g1k_prev(ik)+1,map_z(ib),ik,:)
             zajtmp(ig,:) = zaj_l_prev(ig-ista_g1k_prev(ik)+1,ib,ik,:)
          enddo
          call mpi_allreduce(MPI_IN_PLACE,zajtmp,kg1_prev*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
          call map_zaj_to_fft_box_prev(ik,kg1_prev,nfft_prev,zajtmp,afft)
          if(kimg==1)then
            call dfftw_execute_dft_r2c(plan,afft(1),afft(1))
          else
            call dfftw_execute_dft(plan,afft(1),afft(1))
          endif
          id = fft_box_size_WF_prev(1,0)
          mm = fft_box_size_WF_prev(2,0)
          nl = fft_box_size_WF_prev(1,1)
          nm = fft_box_size_WF_prev(2,1)
          nn = fft_box_size_WF_prev(3,1)
          if(kimg == 1) then
             nlhf = id/2
          else
             nlhf = id
          end if
          cmatr=0.d0
          cmati=0.d0
          do i = 0, nm+1
            do j = 0, nn+1
              do k = 0, nl+1
                ii=i;jj=j;kk=k
                if(i.eq.0) ii = nm
                if(j.eq.0) jj = nn
                if(k.eq.0) kk = nl
                if(i.eq.nm+1) ii = 1
                if(j.eq.nn+1) jj = 1
                if(k.eq.nl+1) kk = 1
                if(kimg == 1 .and. kk > nlhf) then
                   knew = id - kk
                   jnew = nn+2 - jj
                   inew = nm+2 - ii
                   if(jnew > nn) then
                      jnew = jnew - nn
                   end if
                   if(inew > nm) then
                      inew = inew - nm
                   end if
                else
                   knew = kk; jnew = jj; inew = ii
                end if
                ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
                xmat(k+1) = dble(k)/dble(nl)
                ymat(i+1) = dble(i)/dble(nm)
                zmat(j+1) = dble(j)/dble(nn)
                cmatr(k+1,i+1,j+1) = afft(ip*2-1)
                cmati(k+1,i+1,j+1) = afft(ip*2)
              enddo
            enddo
          enddo

          id = fft_box_size_WF(1,0)
          mm = fft_box_size_WF(2,0)
          nl = fft_box_size_WF(1,1)
          nm = fft_box_size_WF(2,1)
          nn = fft_box_size_WF(3,1)
          if(kimg == 1) then
             nlhf = id/2
          else
             nlhf = id
          end if
          bfft=0.d0
          do i = 1, nm
             do j = 1, nn
                do k = 1, nl
                  if(kimg == 1 .and. k > nlhf) then
                     knew = id - k
                     jnew = nn+2 - j
                     inew = nm+2 - i
                     if(jnew > nn) then
                        jnew = jnew - nn
                     end if
                     if(inew > nm) then
                        inew = inew - nm
                     end if
                  else
                     knew = k; jnew = j; inew = i
                  end if
                  ip = nlhf*mm*(jnew-1) + nlhf*(inew-1) + knew
                  bfft(ip*2-1) = trilinear_interpolation( &
               &  nl2,nm2,nn2,cmatr,xmat,ymat,zmat, &
               &  dble(k)/dble(nl),dble(i)/dble(nm),dble(j)/dble(nn))

                  bfft(ip*2) = trilinear_interpolation( &
               &  nl2,nm2,nn2,cmati,xmat,ymat,zmat, &
               &  dble(k)/dble(nl),dble(i)/dble(nm),dble(j)/dble(nn))
                enddo
             enddo
          enddo
          call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
          do ri=1,kimg
            do ig=ista_g1k(ik),iend_g1k(ik)
!              zaj_l(ig-ista_g1k(ik)+1,map_z(ib),ik,ri) = bfft(kimg*igf(nbase(ig,ik))+(ri-kimg))
              zaj_l(ig-ista_g1k(ik)+1,ib,ik,ri) = bfft(kimg*igf(nbase(ig,ik))+(ri-kimg))
            enddo
          enddo
       enddo
    enddo
    deallocate(zaj_l_prev)
    deallocate(zajtmp)
    deallocate(afft)
    deallocate(bfft)
    deallocate(xmat)
    deallocate(ymat)
    deallocate(zmat)
    deallocate(cmatr)
    deallocate(cmati)
    call dfftw_destroy_plan(plan)
    call tstatc0_end(id_sname)

    contains

    subroutine map_zaj_to_fft_box_prev(ik,kg1,nfft,psi_l,bfft)
    integer, intent(in):: ik,kg1,nfft
    real(kind=DP), dimension(kg1,kimg), intent(in) :: psi_l
    real(kind=DP), dimension(nfft), intent(out) :: bfft
    integer :: i,i1,ri, j, i2, ii
    integer :: id_sname2=-1
    call tstatc0_begin('map_zaj_to_fft_box ',id_sname2)

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf_prev(1)
          bfft(i1) = psi_l(1,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_prev(ik)
!!$             i = nbase(ii,1)
             i = nbase_prev(ii,ik)
             i1 = igf_prev(i)
             bfft(i1) = psi_l(ii,1)
             j = nbase_gamma_prev(ii,2)
             i2 = igf_prev(j)
             bfft(i2) =   psi_l(ii,1)
          end do
       else if(kimg == 2) then
          i1 = 2*igf_prev(1) - 1
          bfft(i1)   = psi_l(1,1)
          bfft(i1+1) = psi_l(1,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba_prev(ik)
!!$             i = nbase(ii,1)
             i = nbase_prev(ii,ik)
             i1 = 2*igf_prev(i)-1
             bfft(i1  ) = psi_l(ii,1)
             bfft(i1+1) = psi_l(ii,2)
             j = nbase_gamma_prev(ii,2)
             i2 = 2*igf_prev(j)-1
             bfft(i2  ) = psi_l(ii,1)
             bfft(i2+1) = -psi_l(ii,2)
          end do
       end if
    else
#ifdef NEC_TUNE_SMP
!CDIR NOLOOPCHG
#endif
       do ri = 1, kimg
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do i = 1, iba_prev(ik)
             i1 = kimg*igf_prev(nbase_prev(i,ik)) + (ri - kimg)
             bfft(i1) = psi_l(i,ri)   ! MPI
          end do
       end do
    end if
    end subroutine map_zaj_to_fft_box_prev

  end subroutine m_ES_gen_zaj_from_prev_zaj
#endif


  subroutine m_ES_alloc_zaj_l_prev(np_g1k)
    integer, dimension(kv3), intent(in) :: np_g1k
      if(.not.allocated(zaj_l_prev)) allocate(zaj_l_prev(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
  end subroutine m_ES_alloc_zaj_l_prev

  subroutine m_ES_dealloc_zaj_l_prev()
    if(allocated(zaj_l_prev)) deallocate(zaj_l_prev)
  end subroutine m_ES_dealloc_zaj_l_prev

  subroutine m_ES_dealloc(store_prev)
    logical, intent(in), optional :: store_prev
    logical :: unitcell_can_change
    logical :: store
    store = unitcell_can_change() .and. sw_interpolate_wfs == ON
    if(present(store_prev)) then
      if(store_prev) then
        store = .true.
      endif
    endif
    if(allocated(zaj_l_prev)) deallocate(zaj_l_prev)
    if(store) then
      allocate(zaj_l_prev(maxval(np_g1k),np_e,ista_k:iend_k,kimg))
      zaj_l_prev = zaj_l
    endif

    if(allocated(vlhxcQ)) deallocate(vlhxcQ)
    if(allocated(vlhxc_l)) deallocate(vlhxc_l)
    if(allocated(vnlph_l)) deallocate(vnlph_l)

    if(allocated(zaj_l)) deallocate(zaj_l)
    if(allocated(hlocphi_l)) deallocate(hlocphi_l)
    if(allocated(zaj_l_buf)) deallocate(zaj_l_buf)

    if(allocated(hlocphi_l)) deallocate(hlocphi_l)
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       if(allocated(Phifftr_l)) deallocate(Phifftr_l)
       if(allocated(status_saved_phifftr)) deallocate(status_saved_phifftr)
    end if
#endif
    if(allocated(fsr_l)) deallocate(fsr_l)
    if(allocated(fsi_l)) deallocate(fsi_l)
    if(allocated(eko_l)) deallocate(eko_l)
    if(allocated(occup_l)) deallocate(occup_l)
    if(allocated(fsr_l_2D_k)) deallocate(fsr_l_2D_k)
    if(allocated(fsi_l_2D_k)) deallocate(fsi_l_2D_k)
    if(allocated(zaj_ball)) deallocate(zaj_ball)
    if(allocated(zah_ball)) deallocate(zah_ball)
    if(allocated(neordr)) deallocate(neordr)
    if(allocated(nrvf_ordr)) deallocate(nrvf_ordr)
    if(allocated(fsr_ball)) deallocate(fsr_ball)
    if(allocated(fsi_ball)) deallocate(fsi_ball)
    if(allocated(eko_ek)) deallocate(eko_ek)
    if(allocated(eko1_l)) deallocate(eko1_l)
    if(allocated(evdff)) deallocate(evdff)
    if(allocated(evdffr)) deallocate(evdffr)
    if(allocated(compr_l)) deallocate(compr_l)
    if(allocated(compi_l)) deallocate(compi_l)
    if(allocated(fsr_add_l)) deallocate(fsr_add_l)
    if(allocated(fsi_add_l)) deallocate(fsi_add_l)
    if(allocated(fsr_l_2d)) deallocate(fsr_l_2d)
    if(allocated(fsi_l_2d)) deallocate(fsi_l_2d)
    if(allocated(pre_sc_without)) deallocate(pre_sc_without)
    if(allocated(pre_ss_without)) deallocate(pre_ss_without)
    if(allocated(neordr_old)) deallocate(neordr_old)
    if(allocated(fftindex0)) deallocate(fftindex0)
    if(allocated(fftdist0)) deallocate(fftdist0)
    if(allocated(fftsend0)) deallocate(fftsend0)
    if(allocated(fftrecv0)) deallocate(fftrecv0)
    if(allocated(fftscnt0)) deallocate(fftscnt0)
    if(allocated(fftrcnt0)) deallocate(fftrcnt0)
#ifdef MPI_FFTW
    if(allocated(fftindex0_mpifftw)) deallocate(fftindex0_mpifftw)
    if(allocated(fftdist0_mpifftw)) deallocate(fftdist0_mpifftw)
    if(allocated(fftsend0_mpifftw)) deallocate(fftsend0_mpifftw)
    if(allocated(fftrecv0_mpifftw)) deallocate(fftrecv0_mpifftw)
    if(allocated(fftscnt0_mpifftw)) deallocate(fftscnt0_mpifftw)
    if(allocated(fftrcnt0_mpifftw)) deallocate(fftrcnt0_mpifftw)
#endif
    firstcall = .true.
#ifdef MPI_FFTW
    firstcall_mpifftw = .true.
#endif

  end subroutine m_ES_dealloc


  subroutine m_ES_reset_first_call()
    firstcall = .true.
#ifdef MPI_FFTW
    firstcall_mpifftw = .true.
#endif
  end subroutine m_ES_reset_first_call

!===============================================================================
  subroutine m_ES_Vlocal_in_Rspace_3D(is,afft_l,lsize,iesize,flag,pot_l)
    integer, intent(in) :: is, lsize, iesize, flag
#ifdef FFT_3D_DIVISION
    real(kind=DP), dimension(lsize*2,iesize) :: afft_l
#else
    real(kind=DP), dimension(lsize*kimg,iesize) :: afft_l
#endif
    real(kind=DP), intent(in), optional :: pot_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

    integer,save                   :: id_sname = -1

    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
! === DEBUG by tkato 2012/06/04 ================================================
#ifndef USE_NONBLK_COMM
     integer, allocatable, dimension(:,:) :: rbuf
#endif
! ==============================================================================

#ifdef __FAPP__
    call fapp_start('Vlocal_in_Rspace',1,1)
#endif
                                                  __TIMER_SUB_START(101)
    call tstatc0_begin('m_ES_Vlocal_in_Rspace_3D ', id_sname)

    if (firstcall) then
       call cnt_vlhxc_to_fft_box_3D()
       firstcall = .false.
    endif
    if ( present(pot_l) ) then    !-(contained here) vlhxc_l  --> afft ;using (igf)
       call map_vlhxc_l_onto_afft_3D(is,pot_l)
    else
       call map_vlhxc_l_onto_afft_3D(is,vlhxc_l)
    endif
                                                  __TIMER_COMM_START(190)
#ifdef FFT_3D_DIVISION
    call m_FFT_WF_3DIV_3D(ELECTRON,nfout,afft_l,lsize,iesize,INVERSE) ! afft -> afft
#else
    if (sw_fft_xzy > 0) then
       call m_FFT_WF_3D(ELECTRON,nfout,afft_l,lsize,iesize,INVERSE) ! afft -> afft
    else
       call m_FFT_WF_XYZ_3D(ELECTRON,nfout,afft_l,lsize,iesize,INVERSE) ! afft -> afft
    end if
#endif
                                                  __TIMER_COMM_STOP(190)
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(101)
#ifdef __FAPP__
    call fapp_stop('Vlocal_in_Rspace',1,1)
#endif
  contains

    subroutine cnt_vlhxc_to_fft_box_3D()
      integer, allocatable, dimension(:,:,:) :: xyz
      integer :: max_fft_x, lx, ly, lz, mx, my, mz, mm, i, j, i1
      integer :: iadd, ladd, lrank, iend
      integer, parameter :: itag = 18
#ifdef FFT_3D_DIVISION
      integer :: kx1p, kx2p, kx3p
#endif
                                                  __TIMER_SUB_START(102)
       max_fft_x = maxval(nel_fft_x(:))
       allocate(xyz(2,3,0:nrank_g-1))
#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,111)
       do i = 0, nrank_g - 1
          call mpi_irecv(xyz(1,1,i), 6, mpi_integer, &
         &               i, itag, mpi_ke_world, req_r(i), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,154,ierr)
           endif
       enddo
       do i = 0, nrank_g - 1
          call mpi_isend(xyz_fft_x, 6, mpi_integer, &
         &               i, itag, mpi_ke_world, req_s(i), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,155,ierr)
           endif
       enddo
       call mpi_waitall(nrank_g, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,156,ierr)
        endif
       call mpi_waitall(nrank_g, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(111)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,157,ierr)
        endif
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,266)
! === DEBUG by tkato 2012/06/04 ================================================
!    integer, allocatable, dimension(:,:) :: rbuf
! ==============================================================================
     allocate(rbuf(6,0:nrank_g-1))
     call MPI_ALLGATHER(xyz_fft_x, 6, mpi_integer, &
    &                   rbuf,6, mpi_integer, &
    &                                mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_allgather error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 174, ierr)
     endif
     do i = 1,2
       do j = 1,3
         xyz(i,j,:)=rbuf(i+(j-1)*2,:)
       enddo
     enddo
     deallocate(rbuf)
                                                 __TIMER_COMM_STOP(266)
#endif

       lx = fft_box_size_WF(1,0)
       ly = fft_box_size_WF(2,0)
       lz = fft_box_size_WF(3,0)
#ifdef FFT_3D_DIVISION
       kx1p = fft_X_x_nel
       kx2p = fft_X_y_nel
       kx3p = fft_X_z_nel
#endif
       allocate(fftscnt0(0:nrank_g-1))
       allocate(fftrcnt0(0:nrank_g-1))
       fftscnt0(:) = 0
       fftrcnt0(:) = 0

       allocate(fftindex0(np_kngp_gw))
       allocate(fftdist0(np_kngp_gw))
       allocate(fftsend0(mp_kngp_gw,0:nrank_g-1))
       allocate(fftrecv0(mp_kngp_gw,0:nrank_g-1))
       fftdist0(:) = -1
       fftsend0(:,:) = 0
       fftrecv0(:,:) = 0
       iend = iend_kngp_gw
       if (iend > kg) iend = kg
                                                 __TIMER_DO_START(112)
       do j = ista_kngp_gw, iend
          iadd = j-ista_kngp_gw+1
          i1 = igf(j)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx
          do i = 0, nrank_g-1
#ifdef FFT_3D_DIVISION
             if ((xyz(1,1,i)<=mx).and.(mx<=xyz(2,1,i)).and.    &
            &    (xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
            &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                ladd = mx-xyz(1,1,i)+1+kx1p*(my-xyz(1,2,i))+kx1p*kx2p*(mz-xyz(1,3,i))
#else
             if ((xyz(1,2,i)<=my).and.(my<=xyz(2,2,i)).and.    &
            &    (xyz(1,3,i)<=mz).and.(mz<=xyz(2,3,i))) then
                ladd = mx+lx*(my-xyz(1,2,i))+lx*(xyz(2,2,i)-xyz(1,2,i)+1)*(mz-xyz(1,3,i))
#endif
                fftscnt0(i) = fftscnt0(i) + 1
                fftindex0(iadd) = fftscnt0(i)
                fftdist0 (iadd) = i
                fftsend0(fftscnt0(i),i) = ladd
                exit
             endif
          enddo
       enddo
                                                 __TIMER_DO_STOP(112)
       deallocate(xyz)

#ifdef USE_NONBLK_COMM
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,113)
       do lrank = 0, nrank_g-1
          call mpi_irecv(fftrecv0(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,158,ierr)
           endif
       enddo
       do lrank = 0, nrank_g - 1
          call mpi_isend(fftsend0(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,159,ierr)
           endif
       enddo
       call mpi_waitall(nrank_g, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,160,ierr)
        endif
       call mpi_waitall(nrank_g, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(113)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,161,ierr)
        endif
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,267)
     call MPI_ALLTOALL(fftsend0, mp_kngp_gw, mpi_integer, &
    &                  fftrecv0, mp_kngp_gw, mpi_integer, &
    &                                     mpi_ke_world, ierr )
     if (ierr /= 0) then
        write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 162, ierr)
     endif
                                                 __TIMER_COMM_STOP(267)
#endif

                                                 __TIMER_DO_START(114)
!OCL NOFLTLD
       do i = 0, nrank_g - 1
!OCL NOFLTLD
          do j = 1, nel_kngp_gw(i)
             if (fftrecv0(j,i) == 0) then
                exit
             end if
             fftrcnt0(i) = fftrcnt0(i) + 1
          enddo
       end do
                                                 __TIMER_DO_STOP(114)
       fftmaxsend0 = maxval(fftscnt0)
       fftmaxrecv0 = maxval(fftrcnt0)
                                                 __TIMER_SUB_STOP(102)
    end subroutine cnt_vlhxc_to_fft_box_3D

    subroutine map_vlhxc_l_onto_afft_3D( is, pot_l )
      integer, intent(in) :: is
!      real(kind=DP), intent(in), target :: pot_l(ista_kngp:iend_kngp,kimg,ndim_magmom)
      real(kind=DP), intent(in), target :: pot_l(ista_kngp:iend_kngp,kimg,nspin)

      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, lrank
      integer :: iend, i, j, iadd, mpi_comm_all
      integer, parameter :: itag = 19
!!!ifdef FFT_USE_SSL2
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
      integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================
      real(kind=DP), pointer, dimension(:,:,:) :: vlhxc_p
      real(kind=DP), allocatable, target, dimension(:,:,:) :: vlhxct
      integer :: ist, ien
      vlhxc_p => null()
      if (sw_communicator_for_chg == ON)then
        ist = ista_kngp_gw
        ien = iend_kngp_gw
        allocate(vlhxct(ist:ien,kimg,nspin));vlhxct=0.d0
        do i=ista_kngp,iend_kngp
           vlhxct(i,:,:) = pot_l(i,:,:)
        enddo
        call mpi_allreduce(mpi_in_place,vlhxct,np_kngp_gw*kimg*nspin, &
        &    mpi_double_precision,mpi_sum,mpi_g_world,ierr)
        vlhxc_p => vlhxct
      else
        vlhxc_p => pot_l
        ist = ista_kngp
        ien = iend_kngp
      endif
                                                 __TIMER_SUB_START(103)
!      mpi_comm_all = mpi_k_world(myrank_k)
       mpi_comm_all = mpi_ke_world
                                                 __TIMER_DO_START(115)
       if (fftmaxsend0 /= 0) then
          allocate(sendbuf(fftmaxsend0*kimg,0:nrank_g-1))
          sendbuf = 0.0d0
          if (kimg == 1) then
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0(iadd),fftdist0(iadd)) = vlhxc_p(i,1,is)
             enddo
          else
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0(iadd)*2-1,fftdist0(iadd)) = vlhxc_p(i,1,is)
                sendbuf(fftindex0(iadd)*2,  fftdist0(iadd)) = vlhxc_p(i,2,is)
             enddo
          endif
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(sendbuf(1,1))
! ==============================================================================
       endif
                                                 __TIMER_DO_STOP(115)
       if (fftmaxrecv0 /= 0) then
          allocate(recvbuf(fftmaxrecv0*kimg,0:nrank_g-1))
          recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(recvbuf(1,1))
! ==============================================================================
       endif

#ifndef USE_ALLTOALLV
       icnt_recv = 0
                                                 __TIMER_COMM_START_w_BARRIER(mpi_comm_all,116)
       do lrank = 0, nrank_g - 1
          if (fftrcnt0(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fftrcnt0(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,162,ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
       icnt_send = 0
       do lrank = 0, nrank_g - 1
          if (fftscnt0(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fftscnt0(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,163,ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,164,ierr)
        endif
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(116)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,165,ierr)
        endif
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_comm_all,268)
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=fftmaxsend0*kimg*i
          rdsp(i)=fftmaxrecv0*kimg*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, fftscnt0*kimg, sdsp, &
      &   mpi_double_precision, recvbuf, fftrcnt0*kimg, rdsp, &
      &   mpi_double_precision, mpi_comm_all, ierr )
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 166, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                 __TIMER_COMM_STOP(268)
#endif

       afft_l = 0.0d0
                                                 __TIMER_DO_START(117)
#ifdef FFT_3D_DIVISION
       if (kimg == 1) then
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
!               afft_l(fftrecv0(j,i),1) = recvbuf(j,i)
                afft_l(fftrecv0(j,i)*2-1,1) = recvbuf(j,i)
                afft_l(fftrecv0(j,i)*2  ,1) = 0.0d0
             enddo
          enddo
       else
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
                afft_l(fftrecv0(j,i)*2-1,1) = recvbuf(j*2-1,i)
                afft_l(fftrecv0(j,i)*2  ,1) = recvbuf(j*2  ,i)
             enddo
          enddo
       endif
#else
#ifdef FFT_USE_SSL2
       nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
       ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
       nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
       nxp = nx

       if (kimg == 1) then
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
                iadd = fftrecv0(j,i)
                iz = (iadd-1)/(nx*ny)+1
                nn = mod(iadd,nx*ny)
                if (nn==0) then
                   iy = ny
                else
                   iy = (nn-1)/nx+1
                end if
                ix = mod(nn,nx)
                if(ix==0) ix = nx
                afft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,1) = recvbuf(j,i)
             enddo
          enddo
       else
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
                iadd = fftrecv0(j,i)
                iz = (iadd-1)/(nx*ny)+1
                nn = mod(iadd,nx*ny)
                if (nn==0) then
                   iy = ny
                else
                   iy = (nn-1)/nx+1
                end if
                ix = mod(nn,nx)
                if(ix==0) ix = nx
                afft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,1) = recvbuf(j*2-1,i)
                afft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,1) = recvbuf(j*2  ,i)
             enddo
          enddo
       endif
#else
       if (kimg == 1) then
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
                afft_l(fftrecv0(j,i),1) = recvbuf(j,i)
             enddo
          enddo
       else
          do i = 0, nrank_g - 1
             do j = 1, fftrcnt0(i)
                afft_l(fftrecv0(j,i)*2-1,1) = recvbuf(j*2-1,i)
                afft_l(fftrecv0(j,i)*2  ,1) = recvbuf(j*2  ,i)
             enddo
          enddo
       endif
#endif
#endif
                                                 __TIMER_DO_STOP(117)
       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)
                                                 __TIMER_SUB_STOP(103)
       if(sw_communicator_for_chg == ON) deallocate(vlhxct)
    end subroutine map_vlhxc_l_onto_afft_3D

  end subroutine m_ES_Vlocal_in_Rspace_3D

!------------------------------------------------------------------------------
#ifdef MPI_FFTW
  subroutine m_ES_Vlocal_in_Rspace_mpifftw3d(is,lx,local_n,lz,afft_l_3d,pot_l)
    use m_FFT, only : m_FFT_Inverse_MPI_FFTW, bfft_mpifftw, afft_mpifftw, bfft_mpifftw_kimg1, afft_mpifftw_kimg1

    integer, intent(in) :: is
    integer(C_INTPTR_T), intent(in) :: lx, local_n, lz
    real(kind=DP), dimension(lx,lz,local_n), intent(out) :: afft_l_3d
    real(kind=DP), intent(in), optional :: pot_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    integer :: mm
    integer,save                   :: id_sname = -1

    call tstatc0_begin('m_ES_Vlocal_in_Rspace_mpifftw ', id_sname)

    if (firstcall_mpifftw) then
       call cnt_vlhxc_to_fft_box_3D()
       firstcall_mpifftw = .false.
    endif
    if ( present(pot_l) ) then    !-(contained here) vlhxc_l  --> afft ;using (igf)
       call map_vlhxc_l_onto_afft_3D(is,pot_l)
    else
       call map_vlhxc_l_onto_afft_3D(is,vlhxc_l)
    endif

    call m_FFT_Inverse_MPI_FFTW(nfout)
    call copy_result()
    call tstatc0_end(id_sname)

  contains

    subroutine copy_result()
      integer :: ix,iy,iz
      integer :: lxh

      if(kimg==2) then
        afft_l_3d = real(afft_mpifftw)
      else
        lxh = lx/2
        do iy=1,local_n
        do iz=1,lz
        do ix=1,lxh
          afft_l_3d(ix,iz,iy) = real (afft_mpifftw_kimg1(ix,iz,iy))
        enddo
        enddo
        enddo
      endif
    end subroutine copy_result

    subroutine cnt_vlhxc_to_fft_box_3D()
      integer :: max_fft_x, i, j, i1
      integer :: iadd, ladd, lrank, iend
      integer, parameter :: itag = 18
      integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz
      integer, allocatable, dimension(:) :: local_ns, local_n_offsets
      integer :: mx,my,mz,mm

       lx = fft_box_size_WF(1,0)
       ly = fft_box_size_WF(2,0)
       lz = fft_box_size_WF(3,0)

       alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
       allocate(local_ns(0:nrank_g-1));local_ns=0
       allocate(local_n_offsets(0:nrank_g-1));local_n_offsets=0
       local_ns(myrank_g) = local_n
       local_n_offsets(myrank_g) = local_n_offset
       call mpi_allreduce(mpi_in_place, local_ns, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)
       call mpi_allreduce(mpi_in_place, local_n_offsets, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)

       allocate(fftscnt0_mpifftw(0:nrank_g-1))
       allocate(fftrcnt0_mpifftw(0:nrank_g-1))
       fftscnt0_mpifftw(:) = 0
       fftrcnt0_mpifftw(:) = 0

       allocate(fftindex0_mpifftw(np_kngp_gw))
       allocate(fftdist0_mpifftw(np_kngp_gw))
       allocate(fftsend0_mpifftw(mp_kngp_gw,0:nrank_g-1))
       allocate(fftrecv0_mpifftw(mp_kngp_gw,0:nrank_g-1))
       fftdist0_mpifftw(:) = -1
       fftsend0_mpifftw(:,:) = 0
       fftrecv0_mpifftw(:,:) = 0
       iend = iend_kngp_gw
       if (iend > kg) iend = kg
                                                 __TIMER_DO_START(112)
       do j = ista_kngp_gw, iend
          iadd = j-ista_kngp_gw+1
          i1 = igf(j)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx
          do i = 0, nrank_g-1
             if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                fftscnt0_mpifftw(i) = fftscnt0_mpifftw(i) + 1
                fftindex0_mpifftw(iadd) = fftscnt0_mpifftw(i)
                fftdist0_mpifftw (iadd) = i
                fftsend0_mpifftw(fftscnt0_mpifftw(i),i) = ladd
                exit
             endif
          enddo
       enddo

       deallocate(local_ns)
       deallocate(local_n_offsets)

       do lrank = 0, nrank_g-1
          call mpi_irecv(fftrecv0_mpifftw(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,158,ierr)
           endif
       enddo
       do lrank = 0, nrank_g - 1
          call mpi_isend(fftsend0_mpifftw(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,159,ierr)
           endif
       enddo
       call mpi_waitall(nrank_g, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,160,ierr)
        endif
       call mpi_waitall(nrank_g, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(113)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,161,ierr)
        endif

!OCL NOFLTLD
       do i = 0, nrank_g - 1
!OCL NOFLTLD
          do j = 1, nel_kngp_gw(i)
             if (fftrecv0_mpifftw(j,i) == 0) then
                exit
             end if
             fftrcnt0_mpifftw(i) = fftrcnt0_mpifftw(i) + 1
          enddo
       end do
                                                 __TIMER_DO_STOP(114)
       fftmaxsend0_mpifftw = maxval(fftscnt0_mpifftw)
       fftmaxrecv0_mpifftw = maxval(fftrcnt0_mpifftw)
                                                 __TIMER_SUB_STOP(102)
    end subroutine cnt_vlhxc_to_fft_box_3D

    subroutine map_vlhxc_l_onto_afft_3D(is,pot_l)
      integer, intent(in) :: is
      real(kind=DP), target, intent(in) :: pot_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, lrank
      integer :: iend, i, j, iadd, mpi_comm_all
      integer, parameter :: itag = 19
!!!ifdef FFT_USE_SSL2
      integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
      real(kind=DP), pointer, dimension(:,:,:) :: vlhxc_p
      real(kind=DP), allocatable, target, dimension(:,:,:) :: vlhxct
      integer :: ist, ien
      integer :: i1
      integer :: mx,my,mz,mm
      integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz
      vlhxc_p => null()
      lx = fft_box_size_WF(1,0)
      ly = fft_box_size_WF(2,0)
      lz = fft_box_size_WF(3,0)
      if (sw_communicator_for_chg == ON)then
        ist = ista_kngp_gw
        ien = iend_kngp_gw
        allocate(vlhxct(ist:ien,kimg,nspin));vlhxct=0.d0
        do i=ista_kngp,iend_kngp
           vlhxct(i,:,:) = pot_l(i,:,:)
        enddo
        call mpi_allreduce(mpi_in_place,vlhxct,np_kngp_gw*kimg*nspin, &
        &    mpi_double_precision,mpi_sum,mpi_g_world,ierr)
        vlhxc_p => vlhxct
      else
        vlhxc_p => pot_l
        ist = ista_kngp
        ien = iend_kngp
      endif

       mpi_comm_all = mpi_ke_world

       if (fftmaxsend0_mpifftw /= 0) then
          allocate(sendbuf(fftmaxsend0_mpifftw*kimg,0:nrank_g-1))
          sendbuf = 0.0d0
          if (kimg == 1) then
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0_mpifftw(iadd),fftdist0_mpifftw(iadd)) = vlhxc_p(i,1,is)
             enddo
          else
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0_mpifftw(iadd)*2-1,fftdist0_mpifftw(iadd)) = vlhxc_p(i,1,is)
                sendbuf(fftindex0_mpifftw(iadd)*2,  fftdist0_mpifftw(iadd)) = vlhxc_p(i,2,is)
             enddo
          endif
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(sendbuf(1,1))
! ==============================================================================
       endif
                                                 __TIMER_DO_STOP(115)
       if (fftmaxrecv0_mpifftw /= 0) then
          allocate(recvbuf(fftmaxrecv0_mpifftw*kimg,0:nrank_g-1))
          recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(recvbuf(1,1))
! ==============================================================================
       endif

       icnt_recv = 0
       do lrank = 0, nrank_g - 1
          if (fftrcnt0_mpifftw(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fftrcnt0_mpifftw(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,162,ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
       icnt_send = 0
       do lrank = 0, nrank_g - 1
          if (fftscnt0_mpifftw(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fftscnt0_mpifftw(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,163,ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,164,ierr)
        endif
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(116)
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,165,ierr)
       endif

       if(kimg==2) then
         bfft_mpifftw = (0.0d0,0.0d0)
         do i = 0, nrank_g - 1
            do j = 1, fftrcnt0_mpifftw(i)
              i1 = fftrecv0_mpifftw(j,i)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              bfft_mpifftw(mx,my,mz) = dcmplx(recvbuf(j*2-1,i),recvbuf(j*2,i))
            enddo
         enddo
       else
         bfft_mpifftw_kimg1=0.d0
         do i = 0, nrank_g - 1
            do j = 1, fftrcnt0_mpifftw(i)
              i1 = fftrecv0_mpifftw(j,i)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,lx*ly)
              if (mm==0) then
                 my = ly
              else
                 my = (mm-1)/lx+1
              end if
              mx = mod(mm,lx)
              if(mx==0) mx = lx
              bfft_mpifftw_kimg1(mx,my,mz) = recvbuf(j,i)
            enddo
         enddo
       endif
       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)
       if(sw_communicator_for_chg == ON) deallocate(vlhxct)
    end subroutine map_vlhxc_l_onto_afft_3D

  end subroutine m_ES_Vlocal_in_Rspace_mpifftw3d

!===============================================================================
  subroutine m_ES_Vlocal_in_Rspace_mpifftw(is,afft_l,lsize,iesize,flag)
    use m_FFT, only : m_FFT_Inverse_MPI_FFTW, bfft_mpifftw, afft_mpifftw, bfft_mpifftw_kimg1, afft_mpifftw_kimg1

    integer, intent(in) :: is, lsize, iesize, flag
    real(kind=DP), dimension(lsize*kimg,iesize) :: afft_l

    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    integer :: mm
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz,lxh
    integer,save                   :: id_sname = -1

    call tstatc0_begin('m_ES_Vlocal_in_Rspace_mpifftw ', id_sname)

    if (firstcall_mpifftw) then
       call cnt_vlhxc_to_fft_box_3D()
       firstcall_mpifftw = .false.
    endif
    call map_vlhxc_l_onto_afft_3D(is)   !-(contained here) vlhxc_l  --> afft ;using (igf)
    call m_FFT_Inverse_MPI_FFTW(nfout)
    call to_1d()
    call tstatc0_end(id_sname)

  contains

    subroutine to_1d()
      integer :: ix,iy,iz
      integer :: i1

      lx = fft_box_size_WF(1,0)
      ly = fft_box_size_WF(2,0)
      lz = fft_box_size_WF(3,0)

      if(kimg==2) then
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
        afft_l = 0.d0
        do iy=1,local_n
        do iz=1,lz
        do ix=1,lx
          i1 = (iy-1)*lx*lz+(iz-1)*lx+ix
          afft_l(2*i1-1,1) = real (afft_mpifftw(ix,iz,iy))
          afft_l(2*i1,1)   = aimag(afft_mpifftw(ix,iz,iy))
        enddo
        enddo
        enddo
      else
        lxh = lx/2
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lxh,mpi_ke_world,local_n,local_n_offset)
        afft_l = 0.d0
        do iy=1,local_n
        do iz=1,lz
        do ix=1,lxh
          i1 = (iy-1)*lxh*lz+(iz-1)*lxh+ix
          afft_l(2*i1-1,1) = real (afft_mpifftw_kimg1(ix,iz,iy))
          afft_l(2*i1,1)   = aimag(afft_mpifftw_kimg1(ix,iz,iy))
        enddo
        enddo
        enddo
      endif
    end subroutine to_1d

    subroutine cnt_vlhxc_to_fft_box_3D()
      integer :: max_fft_x, i, j, i1
      integer :: iadd, ladd, lrank, iend
      integer, parameter :: itag = 18
      integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz
      integer, allocatable, dimension(:) :: local_ns, local_n_offsets

       lx = fft_box_size_WF(1,0)
       ly = fft_box_size_WF(2,0)
       lz = fft_box_size_WF(3,0)

       alloc_local = fftw_mpi_local_size_3d(lz,ly,lx,mpi_ke_world,local_n,local_n_offset)
       allocate(local_ns(0:nrank_g-1));local_ns=0
       allocate(local_n_offsets(0:nrank_g-1));local_n_offsets=0
       local_ns(myrank_g) = local_n
       local_n_offsets(myrank_g) = local_n_offset
       call mpi_allreduce(mpi_in_place, local_ns, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)
       call mpi_allreduce(mpi_in_place, local_n_offsets, nrank_g, mpi_integer, mpi_sum, mpi_ke_world, ierr)

       allocate(fftscnt0_mpifftw(0:nrank_g-1))
       allocate(fftrcnt0_mpifftw(0:nrank_g-1))
       fftscnt0_mpifftw(:) = 0
       fftrcnt0_mpifftw(:) = 0

       allocate(fftindex0_mpifftw(np_kngp_gw))
       allocate(fftdist0_mpifftw(np_kngp_gw))
       allocate(fftsend0_mpifftw(mp_kngp_gw,0:nrank_g-1))
       allocate(fftrecv0_mpifftw(mp_kngp_gw,0:nrank_g-1))
       fftdist0_mpifftw(:) = -1
       fftsend0_mpifftw(:,:) = 0
       fftrecv0_mpifftw(:,:) = 0
       iend = iend_kngp_gw
       if (iend > kg) iend = kg
                                                 __TIMER_DO_START(112)
       do j = ista_kngp_gw, iend
          iadd = j-ista_kngp_gw+1
          i1 = igf(j)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx
          do i = 0, nrank_g-1
             if (mz>local_n_offsets(i) .and. mz<=local_n_offsets(i)+local_ns(i)) then
                ladd = mx+lx*(my-1)+lx*ly*(mz-local_n_offsets(i)-1)
                fftscnt0_mpifftw(i) = fftscnt0_mpifftw(i) + 1
                fftindex0_mpifftw(iadd) = fftscnt0_mpifftw(i)
                fftdist0_mpifftw (iadd) = i
                fftsend0_mpifftw(fftscnt0_mpifftw(i),i) = ladd
                exit
             endif
          enddo
       enddo

       deallocate(local_ns)
       deallocate(local_n_offsets)

       do lrank = 0, nrank_g-1
          call mpi_irecv(fftrecv0_mpifftw(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_r(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,158,ierr)
           endif
       enddo
       do lrank = 0, nrank_g - 1
          call mpi_isend(fftsend0_mpifftw(1,lrank), mp_kngp_gw, mpi_integer, &
         &               lrank, itag, mpi_ke_world, req_s(lrank), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,159,ierr)
           endif
       enddo
       call mpi_waitall(nrank_g, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,160,ierr)
        endif
       call mpi_waitall(nrank_g, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(113)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,161,ierr)
        endif

!OCL NOFLTLD
       do i = 0, nrank_g - 1
!OCL NOFLTLD
          do j = 1, nel_kngp_gw(i)
             if (fftrecv0_mpifftw(j,i) == 0) then
                exit
             end if
             fftrcnt0_mpifftw(i) = fftrcnt0_mpifftw(i) + 1
          enddo
       end do
                                                 __TIMER_DO_STOP(114)
       fftmaxsend0_mpifftw = maxval(fftscnt0_mpifftw)
       fftmaxrecv0_mpifftw = maxval(fftrcnt0_mpifftw)
                                                 __TIMER_SUB_STOP(102)
    end subroutine cnt_vlhxc_to_fft_box_3D

    subroutine map_vlhxc_l_onto_afft_3D(is)
      integer, intent(in) :: is
      real(kind=DP), allocatable, dimension(:,:) :: sendbuf, recvbuf
      integer :: icnt_send, icnt_recv, lrank
      integer :: iend, i, j, iadd, mpi_comm_all
      integer, parameter :: itag = 19
!!!ifdef FFT_USE_SSL2
      integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
      real(kind=DP), pointer, dimension(:,:,:) :: vlhxc_p
      real(kind=DP), allocatable, target, dimension(:,:,:) :: vlhxct
      integer :: ist, ien
      integer :: i1
      vlhxc_p => null()
      if (sw_communicator_for_chg == ON)then
        ist = ista_kngp_gw
        ien = iend_kngp_gw
        allocate(vlhxct(ist:ien,kimg,nspin));vlhxct=0.d0
        do i=ista_kngp,iend_kngp
           vlhxct(i,:,:) = vlhxc_l(i,:,:)
        enddo
        call mpi_allreduce(mpi_in_place,vlhxct,np_kngp_gw*kimg*nspin, &
        &    mpi_double_precision,mpi_sum,mpi_g_world,ierr)
        vlhxc_p => vlhxct
      else
        vlhxc_p => vlhxc_l
        ist = ista_kngp
        ien = iend_kngp
      endif

       mpi_comm_all = mpi_ke_world

       if (fftmaxsend0_mpifftw /= 0) then
          allocate(sendbuf(fftmaxsend0_mpifftw*kimg,0:nrank_g-1))
          sendbuf = 0.0d0
          if (kimg == 1) then
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0_mpifftw(iadd),fftdist0_mpifftw(iadd)) = vlhxc_p(i,1,is)
             enddo
          else
             iend = ien
             if (iend > kg) iend = kg
             do i = ist, iend
                iadd = i-ist+1
                sendbuf(fftindex0_mpifftw(iadd)*2-1,fftdist0_mpifftw(iadd)) = vlhxc_p(i,1,is)
                sendbuf(fftindex0_mpifftw(iadd)*2,  fftdist0_mpifftw(iadd)) = vlhxc_p(i,2,is)
             enddo
          endif
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(sendbuf(1,1))
! ==============================================================================
       endif
                                                 __TIMER_DO_STOP(115)
       if (fftmaxrecv0_mpifftw /= 0) then
          allocate(recvbuf(fftmaxrecv0_mpifftw*kimg,0:nrank_g-1))
          recvbuf = 0.0d0
! === DEBUG by tkato 2012/06/05 ================================================
       else
          allocate(recvbuf(1,1))
! ==============================================================================
       endif

       icnt_recv = 0
       do lrank = 0, nrank_g - 1
          if (fftrcnt0_mpifftw(lrank) /= 0) then
             call mpi_irecv(recvbuf(1,lrank), fftrcnt0_mpifftw(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_r(icnt_recv), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_irecv error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,162,ierr)
              endif
             icnt_recv = icnt_recv + 1
          endif
       enddo
       icnt_send = 0
       do lrank = 0, nrank_g - 1
          if (fftscnt0_mpifftw(lrank) /= 0) then
             call mpi_isend(sendbuf(1,lrank), fftscnt0_mpifftw(lrank)*kimg, mpi_double_precision, &
         &                  lrank, itag, mpi_comm_all, req_s(icnt_send), ierr)
              if (ierr /= 0) then
                 write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_isend error'
                 call flush(nfout)
                 call mpi_abort(mpi_comm_world,163,ierr)
              endif
             icnt_send = icnt_send + 1
          endif
       enddo
       call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
        if (ierr /= 0) then
           write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
           call flush(nfout)
           call mpi_abort(mpi_comm_world,164,ierr)
        endif
       call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(116)
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D :  mpi_waitall error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world,165,ierr)
       endif

       lx = fft_box_size_WF(1,0)
       ly = fft_box_size_WF(2,0)
       lz = fft_box_size_WF(3,0)

       if(kimg==2) then
         bfft_mpifftw = (0.0d0,0.0d0)
         do i = 0, nrank_g - 1
            do j = 1, fftrcnt0_mpifftw(i)
              i1 = fftrecv0_mpifftw(j,i)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,(lx*ly))
              if (mm==0) mm=lx*ly
              my = (mm-1)/lx+1
              mx = mod(mm,lx)
              if (mx==0) mx = lx
              bfft_mpifftw(mx,my,mz) = dcmplx(recvbuf(j*2-1,i),recvbuf(j*2,i))
            enddo
         enddo
       else
         bfft_mpifftw_kimg1=0.d0
         do i = 0, nrank_g - 1
            do j = 1, fftrcnt0_mpifftw(i)
              i1 = fftrecv0_mpifftw(j,i)
              mz = (i1-1)/(lx*ly)+1
              mm = mod(i1,lx*ly)
              if (mm==0) then
                 my = ly
              else
                 my = (mm-1)/lx+1
              end if
              mx = mod(mm,lx)
              if(mx==0) mx = lx
              bfft_mpifftw_kimg1(mx,my,mz) = recvbuf(j,i)
            enddo
         enddo
       endif
       if (allocated(sendbuf)) deallocate(sendbuf)
       if (allocated(recvbuf)) deallocate(recvbuf)
       if(sw_communicator_for_chg == ON) deallocate(vlhxct)
    end subroutine map_vlhxc_l_onto_afft_3D

  end subroutine m_ES_Vlocal_in_Rspace_mpifftw

#endif

!------------------------------------------------------------------------------
  subroutine m_ES_init_chgsoft(nfft)
    integer, intent(in) :: nfft
    if(.not.allocated(chg_softpart)) allocate(chg_softpart(nfft,nspin));chg_softpart=0.d0
  end subroutine m_ES_init_chgsoft

  subroutine m_ES_dealloc_chgsoft()
    if(allocated(chg_softpart)) deallocate(chg_softpart)
    chg_has_been_calculated = .false.
  end subroutine m_ES_dealloc_chgsoft
!------------------------------------------------------------------------------

  subroutine m_ES_eigen_values_for_each_k_3D(is,ik,ekin_l,afft_l,lsize,eval_charge)
    integer,       intent(in)                  :: is,ik,lsize    ! is: spin
    real(kind=DP), intent(in), dimension(np_g1k(ik)) :: ekin_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in), dimension(lsize*2) :: afft_l
#else
    real(kind=DP), intent(in), dimension(lsize*kimg) :: afft_l
#endif
    logical,       intent(in), optional        :: eval_charge

    real(kind=DP), allocatable,dimension(:,:) :: bfft_l, bfft_l_in

    integer       :: ib, ib1, ib2, ibsize, ibesize
    real(kind=DP),dimension(np_e) :: eg, eko
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    real(kind=DP),dimension(np_e) :: eko_tmp
! ==============================================================================
    real(kind=DP), allocatable, dimension(:) :: wk_mpi
    real(kind=DP)                    :: prd
    logical       :: chg
    integer       :: id_sname = -1
                                                 __TIMER_SUB_START(601)
    call tstatc0_begin('m_ES_eigen_values_for_each_k_3D ', id_sname,1)

    chg = .false.
    if(present(eval_charge)) chg = eval_charge

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif

#ifdef FFT_3D_DIVISION
    allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
#else
    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(bfft_l_in(lsize*kimg,ibsize) ,stat=ierr);bfft_l_in=0.d0
#endif
     if (ierr /= 0) then
        write(nfout,*)' m_ES_eigen_values_for_each_k_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 178, ierr)
     endif

    eko(:) = 0.0d0
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    if(sw_hybrid_functional == ON) eko_tmp(:) = eko_l(:,ik)
! ==============================================================================
    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part

    do ib1 = 1, np_e, ibsize   ! MPI
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l, 6)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l)
!       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l_in,bfft_l)
#endif
       if(chg) then
         call add_occupied_densities(bfft_l,ib1,ib2,is,ik)
       endif
#ifdef FFT_3D_DIVISION
       call m_FFT_W_Vlocal_W_3DIV_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_W_Vlocal_W_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
       else
          call m_FFT_W_Vlocal_W_XYZ_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
       end if
#endif
    end do

    deallocate(bfft_l)
    deallocate(bfft_l_in)

    prd = product(fft_box_size_WF(1:3,1))

    eko(:) = eko(:) + eg(:) / prd

    call W_Vnonlocal_W(ik)

    if(chg) chg_has_been_calculated = .true.

    if (nrank_g > 1) then
       allocate(wk_mpi(1:np_e))
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,608)
!       call mpi_allreduce(eko, wk_mpi, np_e, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
       call mpi_allreduce(eko, wk_mpi, np_e, mpi_double_precision, mpi_sum, mpi_ske_world, ierr)
                                                 __TIMER_COMM_STOP(608)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + wk_mpi(1:np_e)
       else
          eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       end if
! ==============================================================================
       deallocate(wk_mpi)
    else
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = eko(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + eko(1:np_e)
       else
          eko_l(1:np_e,ik) = eko(1:np_e)
       end if
! ==============================================================================
    end if

    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(601)
  contains
    subroutine W_T_W()
      integer :: ib, i
                                                 __TIMER_SUB_START(602)
                                                 __TIMER_DO_START(609)
      if(kimg==1) then
         do ib = 1, np_e                     ! MPI
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i)*zaj_l(i,ib,ik,1)**2
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      else if(kimg==2) then
!OCL NOFLTLD
         do ib = 1, np_e                     ! MPI
!OCL NOFLTLD
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i) * (zaj_l(i,ib,ik,1)**2+zaj_l(i,ib,ik,2)**2)
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      end if
                                                 __TIMER_DO_STOP(609)
      if(ipri >= 2) then
         write(6,'(" -- eko (W_T_W) --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko(ib),ib=1,np_e) ! MPI
      endif
                                                 __TIMER_SUB_STOP(602)
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)

      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib
      real(kind=DP) :: fac
      integer         :: lnblck, lnb1, lnb2, lnsize
      integer :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)
                                                 __TIMER_SUB_START(603)
#if 1

    if(nblocksize_gather_f_is_given) then
       lnblck = nblocksize_gather_f
       if(np_e < lnblck) then
          lnblck = np_e
       end if
       if(1 > lnblck) then
          lnblck = np_e
       end if
    else
       lnblck = np_e
    end if

    call m_ES_alloc_fsr_l_2d(lnblck, nlmta)
    if(k_symmetry(ik) /= GAMMA) then
       call m_ES_alloc_fsi_l_2d(lnblck, nlmta)
    end if
                                                 __TIMER_DO_START(610)
    do lnb1 = 1, np_e , lnblck
       lnb2=min( lnb1+lnblck-1,np_e )
       lnsize = lnb2-lnb1 + 1
       call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, ik, lnblck, lnb1, lnsize)
       if( k_symmetry(ik) /= GAMMA ) then
          call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi_l_2D, ik, lnblck, lnb1, lnsize)
       endif

      do ia = ista_atm, iend_atm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q))
                  end do
               end do
            else
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q)  &
                    &                       + fsi_l_2D(ib-lnb1+1,p)*fsi_l_2D(ib-lnb1+1,q))
                  end do
               end do
            end if
         end do
      end do
    end do
                                                 __TIMER_DO_STOP(610)
#else

      call m_ES_alloc_fsr_l_2d(np_e, nlmta)
      call m_ES_gather_f_3d_to_2d(fsr_l, fsr_l_2D, ik)
      if(k_symmetry(ik) /= GAMMA) then
         call m_ES_alloc_fsi_l_2d(np_e, nlmta)
         call m_ES_gather_f_3d_to_2d(fsi_l, fsi_l_2D, ik)
      endif
                                                 __TIMER_DO_START(610)
      do ia = ista_atm, iend_atm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = 1, np_e                               ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib,p)*fsr_l_2D(ib,q))
                  end do
               end do
            else
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = 1, np_e                               ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib,p)*fsr_l_2D(ib,q)  &
                    &                       + fsi_l_2D(ib,p)*fsi_l_2D(ib,q))
                  end do
               end do
            end if
         end do
      end do
                                                 __TIMER_DO_STOP(610)
#endif
      call m_ES_dealloc_fsr_l_2d()
      call m_ES_dealloc_fsi_l_2d()
                                                 __TIMER_SUB_STOP(603)
      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W

    subroutine add_occupied_densities(bfft_l,ib1,ib2,is,ik)
      real(kind=DP), intent(in), dimension(lsize*kimg,ibsize) :: bfft_l
      integer, intent(in) :: ib1,ib2,is,ik
      integer  :: i, ib
      real(kind=DP) :: occupation
      do ib = ib1, ib2
         occupation = occup_l(ib,ik)
         if(abs(occupation) < DELTA) cycle
#ifdef FFT_3D_DIVISION
         do i = 1, lsize*2-1, 2
            chg_softpart(i,is) = chg_softpart(i,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
         end do
#else
         if (sw_fft_xzy > 0) then
            do i = 1, np_fft_y*kimg-1, 2
               chg_softpart(i,is) = chg_softpart(i,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         else
            do i = 1, np_fft_z*kimg-1, 2
               chg_softpart(i,is) = chg_softpart(i,is) + occupation*(bfft_l(i,ib-ib1+1)**2+bfft_l(i+1,ib-ib1+1)**2) ! MPI
            end do
         end if
#endif
         if(ipri >= 2) then
            write(nfout,'(" !cdsoft    ik, ib1 = ",2i8)') ik, ib
            write(nfout,'(" !cdsoft     occupation = ",d16.8)') occupation
            write(nfout,'(" !cdsoft     afft : ",5d16.8)') (chg_softpart(i,1),i=1,15)
         end if
      end do
    end subroutine add_occupied_densities

  end subroutine m_ES_eigen_values_for_each_k_3D

!------------------------------------------------------------------------------
#ifdef MPI_FFTW
  subroutine m_ES_eigen_values_for_each_k_mpifftw(is,ik,ekin_l,afft_l,lsize)
    integer,       intent(in)                  :: is,ik,lsize    ! is: spin
    real(kind=DP), intent(in), dimension(np_g1k(ik)) :: ekin_l
    real(kind=DP), intent(in), dimension(lsize*kimg) :: afft_l

    integer       :: ib, ib1, ib2, ibsize, ibesize
    real(kind=DP),dimension(np_e) :: eg, eko
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    real(kind=DP),dimension(np_e) :: eko_tmp
! ==============================================================================
    real(kind=DP), allocatable, dimension(:) :: wk_mpi
    real(kind=DP)                    :: prd
    integer       :: id_sname = -1

    call tstatc0_begin('m_ES_eigen_values_for_each_k_mpifftw ', id_sname,1)

    ibsize = 1
!    if (nblocksize_fftw_is_given) then
!       ibsize = nblocksize_fftw
!       if (ibsize < 1) ibsize = 1
!    endif

    eko(:) = 0.0d0
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    if(sw_hybrid_functional == ON) eko_tmp(:) = eko_l(:,ik)
! ==============================================================================
    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part
    call m_ES_fftbox_map(ik)
    do ib1 = 1, np_e, ibsize   ! MPI
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

       call m_ES_WF_in_Rspace_mpifftw(ista_k,iend_k,ik,ib1,zaj_l)
       call m_FFT_W_Vlocal_W_mpifftw(afft_l,lsize,ibesize,eg(ib1))
    end do

    prd = product(fft_box_size_WF(1:3,1))

    eko(:) = eko(:) + eg(:) / prd

    call W_Vnonlocal_W(ik)

    if (nrank_g > 1) then
       allocate(wk_mpi(1:np_e))
       call mpi_allreduce(eko, wk_mpi, np_e, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + wk_mpi(1:np_e)
       else
          eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       end if
! ==============================================================================
       deallocate(wk_mpi)
    else
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = eko(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + eko(1:np_e)
       else
          eko_l(1:np_e,ik) = eko(1:np_e)
       end if
! ==============================================================================
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine W_T_W()
      integer :: ib, i
      if(kimg==1) then
         do ib = 1, np_e                     ! MPI
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i)*zaj_l(i,ib,ik,1)**2
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      else if(kimg==2) then
!OCL NOFLTLD
         do ib = 1, np_e                     ! MPI
!OCL NOFLTLD
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i) * (zaj_l(i,ib,ik,1)**2+zaj_l(i,ib,ik,2)**2)
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      end if
      if(ipri >= 2) then
         write(6,'(" -- eko (W_T_W) --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko(ib),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)

      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib
      real(kind=DP) :: fac
      integer         :: lnblck, lnb1, lnb2, lnsize
      integer :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)
    if(nblocksize_gather_f_is_given) then
       lnblck = nblocksize_gather_f
       if(np_e < lnblck) then
          lnblck = np_e
       end if
       if(1 > lnblck) then
          lnblck = np_e
       end if
    else
       lnblck = np_e
    end if

    call m_ES_alloc_fsr_l_2d(lnblck, nlmta)
    if(k_symmetry(ik) /= GAMMA) then
       call m_ES_alloc_fsi_l_2d(lnblck, nlmta)
    end if
                                                 __TIMER_DO_START(610)
    do lnb1 = 1, np_e , lnblck
       lnb2=min( lnb1+lnblck-1,np_e )
       lnsize = lnb2-lnb1 + 1
       call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, ik, lnblck, lnb1, lnsize)
       if( k_symmetry(ik) /= GAMMA ) then
          call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi_l_2D, ik, lnblck, lnb1, lnsize)
       endif

      do ia = ista_atm, iend_atm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q))
                  end do
               end do
            else
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q)  &
                    &                       + fsi_l_2D(ib-lnb1+1,p)*fsi_l_2D(ib-lnb1+1,q))
                  end do
               end do
            end if
         end do
      end do
    end do
      call m_ES_dealloc_fsr_l_2d()
      call m_ES_dealloc_fsi_l_2d()
      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W

  end subroutine m_ES_eigen_values_for_each_k_mpifftw

  subroutine m_ES_eigen_values_for_each_k_mpifftw3d(is,ik,lx,local_n,lz,ekin_l,afft_l,eval_charge)
    integer,             intent(in)                     :: is,ik
    integer(C_INTPTR_T), intent(in)                     :: lx,local_n,lz
    real(kind=DP), intent(in), dimension(np_g1k(ik))    :: ekin_l
    real(kind=DP), intent(in), dimension(lx,lz,local_n) :: afft_l
    logical,             intent(in), optional           :: eval_charge

    integer       :: ib, ib1, ib2, ibsize, ibesize
    real(kind=DP),dimension(np_e) :: eg, eko
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    real(kind=DP),dimension(np_e) :: eko_tmp
! ==============================================================================
    real(kind=DP), allocatable, dimension(:) :: wk_mpi
    real(kind=DP)                    :: prd
    logical :: chg
    integer       :: id_sname = -1

    call tstatc0_begin('m_ES_eigen_values_for_each_k_mpifftw ', id_sname,1)
    chg = .false.
    if(present(eval_charge)) chg = eval_charge
    ibsize = 1
!    if (nblocksize_fftw_is_given) then
!       ibsize = nblocksize_fftw
!       if (ibsize < 1) ibsize = 1
!    endif

    eko(:) = 0.0d0
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
    if(sw_hybrid_functional == ON) eko_tmp(:) = eko_l(:,ik)
! ==============================================================================
    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part

    do ib1 = 1, np_e, ibsize   ! MPI
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

       call m_ES_WF_in_Rspace_mpifftw(ista_k,iend_k,ik,ib1,zaj_l)
       if(chg) then
         call add_occupied_densities(ib1,ib2,is,ik)
       endif
       call m_FFT_W_Vlocal_W_mpifftw3d(lx,local_n,lz,afft_l,ibesize,eg(ib1))
    end do

    prd = product(fft_box_size_WF(1:3,1))

    eko(:) = eko(:) + eg(:) / prd

    call W_Vnonlocal_W(ik)

    if(chg) chg_has_been_calculated = .true.

    if (nrank_g > 1) then
       allocate(wk_mpi(1:np_e))
       call mpi_allreduce(eko, wk_mpi, np_e, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + wk_mpi(1:np_e)
       else
          eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       end if
! ==============================================================================
       deallocate(wk_mpi)
    else
! === Support Hybrid on 3D_Parallel by tkato 2013/02/14 ========================
!      eko_l(1:np_e,ik) = eko(1:np_e)
       if(sw_hybrid_functional == ON) then
          eko_l(1:np_e,ik) = eko_tmp(1:np_e) + eko(1:np_e)
       else
          eko_l(1:np_e,ik) = eko(1:np_e)
       end if
! ==============================================================================
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine W_T_W()
      integer :: ib, i
      if(kimg==1) then
         do ib = 1, np_e                     ! MPI
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i)*zaj_l(i,ib,ik,1)**2
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      else if(kimg==2) then
!OCL NOFLTLD
         do ib = 1, np_e                     ! MPI
!OCL NOFLTLD
            do i = 1, np_g1k(ik)
               eko(ib) = eko(ib) + ekin_l(i) * (zaj_l(i,ib,ik,1)**2+zaj_l(i,ib,ik,2)**2)
            end do
            if(k_symmetry(ik) == GAMMA) eko(ib) = 2.d0*eko(ib)
         end do
      end if
      if(ipri >= 2) then
         write(6,'(" -- eko (W_T_W) --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko(ib),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)

      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib
      real(kind=DP) :: fac
      integer         :: lnblck, lnb1, lnb2, lnsize
      integer :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)
    if(nblocksize_gather_f_is_given) then
       lnblck = nblocksize_gather_f
       if(np_e < lnblck) then
          lnblck = np_e
       end if
       if(1 > lnblck) then
          lnblck = np_e
       end if
    else
       lnblck = np_e
    end if

    call m_ES_alloc_fsr_l_2d(lnblck, nlmta)
    if(k_symmetry(ik) /= GAMMA) then
       call m_ES_alloc_fsi_l_2d(lnblck, nlmta)
    end if
                                                 __TIMER_DO_START(610)
    do lnb1 = 1, np_e , lnblck
       lnb2=min( lnb1+lnblck-1,np_e )
       lnsize = lnb2-lnb1 + 1
       call m_ES_gather_f_3d_to_2d_blk(fsr_l, fsr_l_2D, ik, lnblck, lnb1, lnsize)
       if( k_symmetry(ik) /= GAMMA ) then
          call m_ES_gather_f_3d_to_2d_blk(fsi_l, fsi_l_2D, ik, lnblck, lnb1, lnsize)
       endif

      do ia = ista_atm, iend_atm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q))
                  end do
               end do
            else
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = lnb1, lnb2                            ! MPI
                     eko(ib) = eko(ib) + fac*(fsr_l_2D(ib-lnb1+1,p)*fsr_l_2D(ib-lnb1+1,q)  &
                    &                       + fsi_l_2D(ib-lnb1+1,p)*fsi_l_2D(ib-lnb1+1,q))
                  end do
               end do
            end if
         end do
      end do
    end do
      call m_ES_dealloc_fsr_l_2d()
      call m_ES_dealloc_fsi_l_2d()
      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W

    subroutine add_occupied_densities(ib1,ib2,is,ik)
      integer, intent(in) :: ib1,ib2,is,ik
      integer  :: i1, ib, mm
      integer  :: ix, iy, iz
      real(kind=DP) :: occupation
      integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
      integer :: lxh
      lx = fft_box_size_WF(1,0)
      ly = fft_box_size_WF(2,0)
      lz = fft_box_size_WF(3,0)
      if(kimg==2) then
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
      else
        alloc_local = fftw_mpi_local_size_3d(ly,lz,lx/2,mpi_ke_world,local_n,local_n_offset)
      endif
      do ib = ib1, ib2
         occupation = occup_l(ib,ik)
         if(abs(occupation) < DELTA) cycle
         if(kimg==2) then
           do iy=1,local_n
           do iz=1,lz
           do ix=1,lx
             i1=(iy-1)*lx*lz+(iz-1)*lx+ix
             chg_softpart(i1*2-1,is) = chg_softpart(i1*2-1,is)+occupation*(real(afft_mpifftw(ix,iz,iy))**2 &
                                     +                                    aimag(afft_mpifftw(ix,iz,iy))**2)
           enddo
           enddo
           enddo
         else
           lxh = lx/2
           do iy=1,local_n
           do iz=1,lz
           do ix=1,lxh
             i1=(iy-1)*lxh*lz+(iz-1)*lxh+ix
             chg_softpart(i1*2-1,is) = chg_softpart(i1*2-1,is)+occupation*(real(afft_mpifftw_kimg1(ix,iz,iy))**2 &
                                     +                                    aimag(afft_mpifftw_kimg1(ix,iz,iy))**2)
           enddo
           enddo
           enddo
         endif
      end do
    end subroutine add_occupied_densities

  end subroutine m_ES_eigen_values_for_each_k_mpifftw3d
#endif
!------------------------------------------------------------------------------

  subroutine m_ES_gather_f_3d_to_2d_k(inn,out)
    real(kind=DP), dimension(:,:), intent(in)  :: inn    !fsr_3D
    real(kind=DP), dimension(:,:)  , intent(out) :: out    !fsr_2D
!   real(kind=DP), dimension(np_e*mp_fs,0:nrank_g-1) :: recvbuf
    integer(kind=4) :: ierr , i, j, k

    integer(kind=4),dimension(0:nrank_g-1) :: recvcnt, recvdsp
                                                 __TIMER_SUB_START(203)

!   call mpi_allgather(inn(1,1,ik), np_e*np_fs, MPI_DOUBLE_PRECISION, &
!  &                   recvbuf, np_e*mp_fs, MPI_DOUBLE_PRECISION, mpi_ke_world, ierr)

!   do k = 0, nrank_g - 1
!      do j = 1, nel_fs(k)
!         do i = 1, np_e
!            out(i,nis_fs(k)+j-1) = recvbuf(np_e*(j-1)+i,k)
!         enddo
!      enddo
!   enddo
                                                 __TIMER_DO_START(210)
    do i = 0, nrank_g - 1
       recvcnt(i) = nel_fs(i) * np_e
    end do
                                                 __TIMER_DO_STOP(210)
    recvdsp(0) = 0
                                                 __TIMER_DO_START(211)
    do i = 1, nrank_g - 1
       recvdsp(i) = recvdsp(i-1) + recvcnt(i-1)
    end do
                                                 __TIMER_DO_STOP(211)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,212)
    call mpi_allgatherv(inn(1,1), np_e*np_fs, MPI_DOUBLE_PRECISION, out, &
   &                    recvcnt, recvdsp, MPI_DOUBLE_PRECISION, mpi_ke_world, ierr)
                                                 __TIMER_COMM_STOP(212)

                                                 __TIMER_SUB_STOP(203)
  end subroutine m_ES_gather_f_3d_to_2d_k

  subroutine m_ES_gather_f_3d_to_2d(inn,out,ik)
    real(kind=DP), dimension(np_e,np_fs,ista_k:iend_k), intent(in)  :: inn    !fsr_3D
    real(kind=DP), dimension(:,:)  , intent(out) :: out    !fsr_2D
    integer(kind=4), intent(in) :: ik
!   real(kind=DP), dimension(np_e*mp_fs,0:nrank_g-1) :: recvbuf
    integer(kind=4) :: ierr , i, j, k

    integer(kind=4),dimension(0:nrank_g-1) :: recvcnt, recvdsp
                                                 __TIMER_SUB_START(203)

!   call mpi_allgather(inn(1,1,ik), np_e*np_fs, MPI_DOUBLE_PRECISION, &
!  &                   recvbuf, np_e*mp_fs, MPI_DOUBLE_PRECISION, mpi_ke_world, ierr)

!   do k = 0, nrank_g - 1
!      do j = 1, nel_fs(k)
!         do i = 1, np_e
!            out(i,nis_fs(k)+j-1) = recvbuf(np_e*(j-1)+i,k)
!         enddo
!      enddo
!   enddo
                                                 __TIMER_DO_START(210)
    do i = 0, nrank_g - 1
       recvcnt(i) = nel_fs(i) * np_e
    end do
                                                 __TIMER_DO_STOP(210)
    recvdsp(0) = 0
                                                 __TIMER_DO_START(211)
    do i = 1, nrank_g - 1
       recvdsp(i) = recvdsp(i-1) + recvcnt(i-1)
    end do
                                                 __TIMER_DO_STOP(211)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,212)
    call mpi_allgatherv(inn(1,1,ik), np_e*np_fs, MPI_DOUBLE_PRECISION, out, &
   &                    recvcnt, recvdsp, MPI_DOUBLE_PRECISION, mpi_ke_world, ierr)
                                                 __TIMER_COMM_STOP(212)

                                                 __TIMER_SUB_STOP(203)
  end subroutine m_ES_gather_f_3d_to_2d

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  subroutine m_ES_gather_f_3d_to_2d_blk(inn,out,ik,iasize,ista,ibsize)
! ==== DEBUG by tkato 2012/04/03 ===============================================
!   real(kind=DP), dimension(:,:,:), intent(in)  :: inn    !fsr_3D
!   real(kind=DP), dimension(:,:)  , intent(inout) :: out    !fsr_2D
    real(kind=DP), dimension(np_e, np_fs, ista_k:iend_k), intent(in)  :: inn ! fsr_l
    real(kind=DP), dimension(iasize, nlmta), intent(inout) :: out            ! fsr_l_2D
! ==============================================================================
    integer(kind=4), intent(in) :: ik, iasize, ista, ibsize
! === for np_fs == 0 ===========================================================
!   real(kind=DP), dimension(iasize,np_fs) :: buf
    real(kind=DP), dimension(iasize,max(np_fs,1)) :: buf
! ==============================================================================
    integer(kind=4) :: ierr , i, j, k

    integer(kind=4),dimension(0:nrank_g-1) :: recvcnt, recvdsp
                                                 __TIMER_SUB_START(203)
                                                 __TIMER_DO_START(210)
    do j = 1, np_fs
       do i = 1, ibsize
          buf(i,j) = inn(ista+i-1,j,ik)
       end do
    end do
                                                 __TIMER_DO_STOP(210)
                                                 __TIMER_DO_START(211)
    do i = 0, nrank_g - 1
       recvcnt(i) = nel_fs(i) * iasize
    end do
    recvdsp(0) = 0
    do i = 1, nrank_g - 1
       recvdsp(i) = recvdsp(i-1) + recvcnt(i-1)
    end do
                                                 __TIMER_DO_STOP(211)
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,212)
    call mpi_allgatherv(buf(1,1), iasize*np_fs, MPI_DOUBLE_PRECISION, out, &
   &                    recvcnt, recvdsp, MPI_DOUBLE_PRECISION, mpi_ke_world, ierr)
                                                 __TIMER_COMM_STOP(212)
                                                 __TIMER_SUB_STOP(203)
  end subroutine m_ES_gather_f_3d_to_2d_blk

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  subroutine m_ES_sort_eigen_values_3D()
    integer             :: ik, ib
    real(kind=DP), allocatable, dimension(:)     :: eko_t
    real(kind=DP),              dimension(neg)   :: eko_s
    integer,                    dimension(neg)   :: k_ordr, tr_neg
    integer                            :: neg_noskip, neg_noskip_local
    integer, allocatable, dimension(:) :: ib_noskip, ib_noskip_sort ! d(neg_noskip)
!!$    integer, allocatable, dimension(:) :: ib_noskip_local           ! d(neg_noskip_local)

    integer                  :: i, id_sname = -1
                                                 __TIMER_SUB_START(604)
    call tstatc0_begin('m_ES_sort_eigen_values_3D ',id_sname,1)
                                                 __TIMER_DO_START(612)
    do i = 1, neg
       tr_neg(neg_g_all(i)) = i
    end do
                                                 __TIMER_DO_STOP(612)
    allocate(eko_t(mp_e*nrank_e), stat=ierr)
    if (ierr/=0) then
      call phase_error_with_msg('failed to allocate temporary array at m_ES_sort_eigen_values_3D', &
      &                        __LINE__, __FILE__)
    endif

    do ik = 1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle    ! MPI

       call expand_eko_l_to_eko_s(ik) ! -> eko_s

       k_ordr(1:neg) = neordr(1:neg,ik)
       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)

       if(ipri >= 2) call wd_eko_s(ik,1)

#ifdef _NO_HEAP_SORT_EIGENVALUES_
       call bubble_sorting(neg,eko_s,neordr(1,ik))
                                                 __TIMER_DO_START(618)
       do ib = 1, neg
          eko_t(ib) = eko_s(neordr(ib,ik))
       enddo
                                                 __TIMER_DO_STOP(618)
#else
       call heap_sorting(neg,eko_s,neordr(1,ik))
       eko_t(1:neg) = eko_s(1:neg)
#endif

       call set_nrvf_ordr(ik)           ! neordr(:,ik) -> nrvf_ordr(:,ik)
       call set_neg_noskip(ik)          ! ib_noskip, ib_noskip_sort
       call sort_eko_l_using_eko_t(ik)  ! eko_t -> eko_s -> eko_l
       call sort_zaj(ik)
#ifdef SAVE_FFT_TIMES
       if(sw_save_fft == ON) then
          call sort_Phifftr_l(ik)
          call sort_status_saved_phifftr(ik)
       endif
#endif
       call sort_fsrfsi(ik)
       call dealloc_ib_noskips()          ! ib_noskip, ib_noskip_sort

       neordr_old(:,ik) = neordr(:,ik)
       call set_neordr_and_nrvf_ordr(ik)

       if(ipri >= 2) call wd_eko_s(ik,2)
    end do ! do-loop of ik

    if(af /= 0) then
       call cp_eigen_values_for_af()
       call expand_neordr_and_nrvf_ordr()
    end if

    deallocate(eko_t)
    call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(604)
  contains
!!$#ifdef SAVE_FFT_TIMES
    subroutine set_nrvf_ordr(ik)
      integer, intent(in) :: ik
      integer :: ib,jb

      do ib = 1, neg
         do jb = 1, neg
            if(ib == neordr(jb,ik)) then
               nrvf_ordr(ib,ik) = jb
               exit
            endif
         enddo
      enddo
    end subroutine set_nrvf_ordr
!!$#endif

    subroutine set_neordr_and_nrvf_ordr(ik)
      integer, intent(in) :: ik
      integer, allocatable, dimension(:) :: t_ordr
      integer :: ib, jb

      allocate(t_ordr(neg))
                                                 __TIMER_DO_START(630)
      do ib = 1, neg
         t_ordr(ib) = k_ordr(neordr(ib,ik))
      end do
                                                 __TIMER_DO_STOP(630)
      neordr(1:neg,ik) = t_ordr(1:neg)
                                                 __TIMER_DO_START(631)
      do ib = 1, neg
         do jb = 1, neg
            if(ib == neordr(jb,ik)) then
               nrvf_ordr(ib,ik) = jb
               exit
            endif
         enddo
      enddo
                                                 __TIMER_DO_STOP(631)
      deallocate(t_ordr)
    end subroutine set_neordr_and_nrvf_ordr

    subroutine wd_eko_s(ik,mode)
      integer, intent(in) :: ik,mode
      integer :: i
      if(mode == 1) then
         write(nfout,'(" !Esort  eko_s before sorting ")')
      else if(mode == 2) then
         write(nfout,'(" !Esort  eko_s in order ")')
      else
         return
      end if
      write(nfout,'(" !Esort ",10f8.4)') (eko_s(neordr(i,ik)),i=1,neg)
    end subroutine wd_eko_s

    subroutine sort_eko_l_using_eko_t(ik)
      integer, intent(in) :: ik
      integer :: ib
                                                 __TIMER_DO_START(619)
      do ib = 1, neg
         eko_s(tr_neg(ib)) = eko_t(ib)
      enddo
                                                 __TIMER_DO_STOP(619)
                                                 __TIMER_DO_START(620)
      do ib = ista_e, iend_e
         eko_l(ib-ista_e+1,ik) = eko_s(ib)
      enddo
                                                 __TIMER_DO_STOP(620)
    end subroutine sort_eko_l_using_eko_t

    subroutine expand_eko_l_to_eko_s(ik)
      integer, intent(in) :: ik
      real(kind=DP), allocatable, dimension(:)     :: eko_t, eko_t2
      integer :: ib, kb, jb, iadd

      allocate(eko_t(mp_e*nrank_e), stat=ierr)
      allocate(eko_t2(mp_e), stat=ierr); eko_t2 = 0.d0
      if (ierr /= 0) then
         write(nfout,'("Not allocate error")')
         call mpi_abort(mpi_comm_world, 179, ierr)
      end if

                                                 __TIMER_DO_START(613)
      do ib = 1, np_e
         eko_t2(ib) = eko_l(ib,ik)
      end do
                                                 __TIMER_DO_STOP(613)
      if(nrank_e > 1) then
                                                 __TIMER_COMM_START_w_BARRIER(mpi_kg_world,614)
         call mpi_allgather(eko_t2, mp_e, mpi_double_precision &
              &                 , eko_t,  mp_e, mpi_double_precision, mpi_kg_world, ierr)
                                                 __TIMER_COMM_STOP(614)
         iadd = 0
                                                 __TIMER_DO_START(615)
         do kb = 0, nrank_e-1
            if (nel_e(kb) == mp_e) then
               iadd = iadd + mp_e
            else
               do jb = 1, nel_e(kb)
                  iadd = iadd + 1
                  eko_t(iadd) = eko_t(mp_e*kb+jb)
               enddo
            end if
         enddo
                                                 __TIMER_DO_STOP(615)
      else
         eko_t = eko_t2
      end if
                                                 __TIMER_DO_START(616)
      do ib = 1, neg
         eko_s(neg_g_all(ib)) = eko_t(ib)
      enddo
                                                 __TIMER_DO_STOP(616)
      deallocate(eko_t)
      deallocate(eko_t2)
    end subroutine expand_eko_l_to_eko_s

    subroutine set_neg_noskip(ik)
      integer, intent(in) :: ik
      integer :: icount, ib, jb
      ! --- global ---
      icount = 0
      do ib = 1, neg
         jb = neg_g_all(ib)
         if(neordr(jb,ik) == jb) cycle
         icount = icount+1
      end do
      neg_noskip = icount

      if(neg_noskip >= 1) then
         allocate(ib_noskip(neg_noskip))
         allocate(ib_noskip_sort(neg_noskip))
      end if

      icount = 0
      do ib = 1, neg
         jb = neg_g_all(ib)
         if(neordr(jb,ik) == jb) cycle
         icount = icount+1
         ib_noskip(icount) = ib
         ib_noskip_sort(icount) = tr_neg(nrvf_ordr(jb,ik))
      end do

!!$      ! --- local ---
!!$      icount = 0
!!$      do ib = ista_e, iend_e
!!$         jb = neg_g_all(ib)
!!$         if(neordr(jb,ik) == jb) cycle
!!$         icount = icount+1
!!$      end do
!!$      neg_noskip_local = icount
!!$
!!$      if(neg_noskip_local >= 1) then
!!$         allocate(ib_noskip_local(neg_noskip_local))
!!$      end if
!!$      icount = 0
!!$      do ib = ista_e, iend_e
!!$         jb = neg_g_all(ib)
!!$         if(neordr(jb,ik) == jb) cycle
!!$         icount = icount + 1
!!$         ib_noskip_local(icount) = ib
!!$      end do

    end subroutine set_neg_noskip

    subroutine dealloc_ib_noskips()
      if(neg_noskip >= 1) then
         deallocate(ib_noskip)
         deallocate(ib_noskip_sort)
      end if
!!$      if(neg_noskip_local >= 1) then
!!$         deallocate(ib_noskip_local)
!!$      end if
    end subroutine dealloc_ib_noskips

    subroutine sort_zaj(ik)
      integer, intent(in) :: ik
      real(kind=DP), allocatable, dimension(:,:,:) :: wk_zaj
      integer :: kb, jb, ib, ibr, jb2, i

!!$      max_g1k = maxval(np_g1k(:))
      allocate(wk_zaj(maxval(np_g1k(:)),neg,kimg), stat=ierr)
      do kb = 1, kimg
                                                 __TIMER_DO_START(621)
!OCL NOFLTLD
         do jb = 1, neg_noskip
            ib  = ib_noskip(jb)
            jb2 = ib_noskip_sort(jb)
!OCL NOFLTLD
            do i = 1, np_g1k(ik)
               wk_zaj(i,jb2,kb) = zaj_ball(i,ib,ik,kb)
            end do
         end do
                                                 __TIMER_DO_STOP(621)
                                                 __TIMER_DO_START(622)
         do jb = 1, neg_noskip
            ib = ib_noskip(jb)
            do i = 1, np_g1k(ik)
               zaj_ball(i,ib,ik,kb) = wk_zaj(i,ib,kb)
            end do
         end do
                                                 __TIMER_DO_STOP(622)
                                                 __TIMER_DO_START(623)
!!$         do ib = 1, neg_noskip_local
!!$            jb = ib_noskip_local(ib)
!!$            do i = 1, np_g1k(ik)
!!$               zaj_l(i,jb-ista_e+1,ik,kb) = zaj_ball(i,jb,ik,kb)
!!$            end do
!!$         end do
         do jb = ista_e, iend_e
!!!$            ib = jb + ista_e - 1
            ibr = neg_g_all(jb)
            if(neordr(ibr,ik) == ibr) cycle
            do ib = 1, np_g1k(ik)
!!!$               zaj_l(ib,jb-ista_e+1,ik,kb) = zaj_ball(ib,jb,ik,kb)
               zaj_l(ib,jb-ista_e+1,ik,kb) = zaj_ball(ib,jb,ik,kb)
            end do
         end do
                                                 __TIMER_DO_STOP(623)
      end do
      deallocate(wk_zaj)
    end subroutine sort_zaj

    subroutine sort_fsrfsi(ik)
      integer, intent(in) :: ik
#ifdef SORT_FSRFSI_WITH_SKIP_ARRAYS
      integer :: jb, ib, ibr, jb2, i
      real(kind=DP), allocatable, dimension(:,:)   :: wk_fsr, wk_fsi
#else
      integer :: ib,jb
      real(kind=DP), allocatable, dimension(:,:)   :: wk_fsr
#endif
      allocate(wk_fsr(neg,np_fs), stat=ierr)
#ifndef SORT_FSRFSI_WITH_SKIP_ARRAYS
                                                 __TIMER_DO_START(624)
!OCL NOFLTLD
      do jb = 1, np_fs
         do ib = 1, neg
            wk_fsr(neg_g_all(ib),jb) = fsr_ball(ib,jb,ik)
         end do
      end do
                                                 __TIMER_DO_STOP(624)
                                                 __TIMER_DO_START(625)
!OCL NOFLTLD
      do jb = 1, np_fs
         do ib = 1, neg
            fsr_ball(tr_neg(ib),jb,ik) = wk_fsr(neordr(ib,ik),jb)
         end do
      end do
                                                 __TIMER_DO_STOP(625)
                                                 __TIMER_DO_START(626)
      do jb = 1, np_fs
         do ib = ista_e, iend_e
            fsr_l(ib-ista_e+1,jb,ik) = fsr_ball(ib,jb,ik)
         end do
      end do
                                                 __TIMER_DO_STOP(626)
      if(k_symmetry(ik) /= GAMMA) then
!!$         allocate(wk_fsr(neg,np_fs), stat=ierr)
                                                 __TIMER_DO_START(627)
         do jb = 1, np_fs
            do ib = 1, neg
               wk_fsr(neg_g_all(ib),jb) = fsi_ball(ib,jb,ik)
            end do
         end do
                                                 __TIMER_DO_STOP(627)
                                                 __TIMER_DO_START(628)
         do jb = 1, np_fs
            do ib = 1, neg
               fsi_ball(tr_neg(ib),jb,ik) = wk_fsr(neordr(ib,ik),jb)
            end do
         end do
                                                 __TIMER_DO_STOP(628)
                                                 __TIMER_DO_START(629)
         do jb = 1, np_fs
            do ib = ista_e, iend_e
               fsi_l(ib-ista_e+1,jb,ik) = fsi_ball(ib,jb,ik)
            end do
         end do
                                                 __TIMER_DO_STOP(629)
      endif
#endif
#ifdef SORT_FSRFSI_WITH_SKIP_ARRAYS
      if(k_symmetry(ik) /= GAMMA) then
         allocate(wk_fsi(neg,np_fs), stat=ierr)
!OCL NOFLTLD
         do i = 1, np_fs
!!$            do ib = 1, neg
            do ib = 1, neg_noskip
!!$               jb = neg_g_all(ib)
!!$               if(neordr(jb,ik) == jb) cycle
!!$               jb2 = tr_neg(nrvf_ordr(jb,ik))
               jb  = ib_noskip(ib)
               jb2 = ib_noskip_sort(ib)
               wk_fsr(jb2,i) = fsr_ball(jb,i,ik)
               wk_fsi(jb2,i) = fsi_ball(jb,i,ik)
            end do
         end do
!OCL NOFLTLD
         do i = 1, np_fs
!!$            do ib = 1, neg
!!$               jb = neg_g_all(ib)
!!$               if(neordr(jb,ik) == jb) cycle
            do ib = 1, neg_noskip
               jb = ib_noskip(ib)
               fsr_ball(jb,i,ik) = wk_fsr(jb,i)
               fsi_ball(jb,i,ik) = wk_fsi(jb,i)
            end do
         end do
         do i = 1, np_fs
!!$            do jb = 1, neg_noskip_local
!!$               ib = ib_noskip_local(jb)
            do ib = ista_e, iend_e
               jb = neg_g_all(ib)
               if(neordr(jb,ik) == jb) cycle
               fsr_l(ib-ista_e+1,i,ik) = fsr_ball(ib,i,ik)
               fsi_l(ib-ista_e+1,i,ik) = fsi_ball(ib,i,ik)
            end do
         end do
         deallocate(wk_fsi)
      else
!OCL NOFLTLD
         do i = 1, np_fs
!!$            do ib = 1, neg
!!$               jb = neg_g_all(ib)
!!$               if(neordr(jb,ik) == jb) cycle
!!$               ibr = nrvf_ordr(jb,ik)
!!$               jb2 = tr_neg(ibr)
            do ib = 1, neg_noskip
               jb  = ib_noskip(ib)
               jb2 = ib_noskip_sort(ib)
               wk_fsr(jb2,i) = fsr_ball(jb,i,ik)
            end do
         end do
!OCL NOFLTLD
         do i = 1, np_fs
!!$            do ib = 1, neg
!!$               jb = neg_g_all(ib)
!!$               if(neordr(jb,ik) == jb) cycle
            do ib = 1, neg_noskip
               jb = ib_noskip(ib)
               fsr_ball(jb,i,ik) = wk_fsr(jb,i)
            end do
         end do
         do i = 1, np_fs
            do jb = ista_e, iend_e
               ibr = neg_g_all(jb)
               if(neordr(ibr,ik) == ibr) cycle
               fsr_l(ib-ista_e+1,i,ik) = fsr_ball(ib,i,ik)
            end do
!!$            do jb = 1, neg_noskip_local
!!$               ib = ib_noskip_local(jb)
!!$               fsr_l(ib-ista_e+1,i,ik) = fsr_ball(ib,i,ik)
!!$            end do
         end do
      endif
#endif
      deallocate(wk_fsr)
    end subroutine sort_fsrfsi

#ifdef SAVE_FFT_TIMES
    subroutine sort_Phifftr_l(ik)
#ifdef SAVE_FFT_TIMES02
      use m_Parallelization, only :  nis_e, nie_e
#endif
      integer, intent(in) :: ik
      integer                                      :: lsize, ib, jb
#ifdef SAVE_FFT_TIMES00
      real(kind=DP), allocatable, dimension(:,:)   :: Phifftr_ball  !d(lsize*kimg,neg)
      integer, parameter :: SIZE_OF_MPIALLREDUCE_per_neg = 500
      integer         :: nloop,msize,istart,iend,il
#else
#ifdef SAVE_FFT_TIMES02
      real(kind=DP), allocatable, dimension(:,:)   :: Phifftr_ball  !d(lsize*kimg,neg)
      real(kind=DP), allocatable, dimension(:,:)   :: Phifftr_wk  !d(lsize*kimg,mp_e)
      real(kind=DP), allocatable, dimension(:,:,:) :: Phifftr_mpi  !d(lsize*kimg,mp_e,nrank_e)
      integer, parameter :: SIZE_OF_MPIALLREDUCE_per_neg = 1000
      integer         :: nloop,msize,istart,iend,il
      integer         :: nb
#else
      real(kind=DP), allocatable, dimension(:)     :: Phifftr_wk  !d(lsize*kimg)
      integer :: nbs
#ifdef SAVE_FFT_TIMES01
      real(kind=DP), allocatable, dimension(:,:)   :: Phifftr_ball  !d(lsize*kimg,neg)
#else
      real(kind=DP), allocatable, dimension(:,:)   :: Phifftr_trans  !d(lsize*kimg,np_e)
      integer :: ibr, nbsr
!!$      integer  , allocatable, dimension(:)   :: nrvf_ordr_t ! d(neg)
#endif
#endif
#endif

#ifdef FFT_3D_DIVISION
       lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
       lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif

#ifdef SAVE_FFT_TIMES00
       msize = min(SIZE_OF_MPIALLREDUCE_per_neg,lsize)
       nloop = lsize*kimg/msize
       if(nloop*msize < lsize*kimg) nloop = nloop + 1
       istart = 1; iend = min(msize,lsize*kimg)
       allocate(Phifftr_ball(msize,neg)); Phifftr_ball = 0.d0
       do il = 1, nloop
          Phifftr_ball = 0.d0
          do ib = 1, np_e
             do i = istart, iend
                Phifftr_ball(i-istart+1,neg_g(ib)) = Phifftr_l(i,ib,ik)
             end do
          end do
          call mpi_allreduce(MPI_IN_PLACE,Phifftr_ball,msize*neg,MPI_DOUBLE_PRECISION &
               &                      ,MPI_SUM,mpi_kg_world,ierr)

          do ib = 1, neg
             if(neordr(ib,ik) == ib) cycle
             if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
                jb = tr_neg(ib)
                do i = istart, iend
                   Phifftr_l(i,jb-ista_e+1,ik) = Phifftr_ball(i-istart+1,neordr(ib,ik))
                end do
             end if
          end do
          istart = istart+msize; iend = min(iend+msize,lsize*kimg)
       end do
       deallocate(Phifftr_ball)
#else
#ifdef SAVE_FFT_TIMES02

       msize = min(SIZE_OF_MPIALLREDUCE_per_neg,lsize)
       nloop = lsize*kimg/msize
       if(nloop*msize < lsize*kimg) nloop = nloop + 1

       istart = 1; iend = min(msize,lsize*kimg)
       allocate(Phifftr_ball(msize,neg), stat=ierr)
       do il = 1, nloop
          allocate(Phifftr_wk(msize,mp_e),stat=ierr);      Phifftr_wk = 0.d0
          allocate(Phifftr_mpi(msize,mp_e,0:nrank_e-1),stat=ierr)
          do ib = 1, np_e
             do i = istart, iend
                Phifftr_wk(i-istart+1,ib) = Phifftr_l(i,ib,ik)
             end do
          end do
          call mpi_allgather(Phifftr_wk, msize*mp_e, MPI_DOUBLE_PRECISION &
               &            ,Phifftr_mpi,msize*mp_e, MPI_DOUBLE_PRECISION, mpi_kg_world,ierr)
          deallocate(Phifftr_wk)
          do nb = 0, nrank_e-1
             do ib = nis_e(nb), nie_e(nb)
                do i = 1, msize
                   Phifftr_ball(i,neg_g_all(ib)) = Phifftr_mpi(i,ib-nis_e(nb)+1,nb)
                end do
             end do
          end do
          deallocate(Phifftr_mpi)
          do ib = 1, neg
             if(neordr(ib,ik) == ib) cycle
             if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
                jb = tr_neg(ib)
                do i = istart, iend
                   Phifftr_l(i,jb-ista_e+1,ik) = Phifftr_ball(i-istart+1,neordr(ib,ik))
                end do
             end if
          end do
          istart = istart+msize; iend = min(iend+msize,lsize*kimg)
       end do
       deallocate(Phifftr_ball)
#else
#ifdef SAVE_FFT_TIMES01
       allocate(Phifftr_ball(lsize*kimg,neg)); Phifftr_ball = 0.d0
       allocate(Phifftr_wk(lsize*kimg))
       do ib = 1, neg
          if(neordr(ib,ik) == ib) cycle
          nbs = (ib-1)/NB + 1
          if(lrank(nbs) == myrank_e) Phifftr_wk = Phifftr_l(:,map_z(ib),ik)
          call mpi_bcast(Phifftr_wk,lsize*kimg,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
!!$          Phifftr_ball(:,neg_g_all(ib)) = Phifftr_wk  ! <-----
          Phifftr_ball(:,ib) = Phifftr_wk  ! <-----
       end do
       deallocate(Phifftr_wk)
       do ib = 1, neg
!!$          jb = tr_neg(ib)
          if(neordr(ib,ik) == ib) cycle
          if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
             jb = tr_neg(ib)
             Phifftr_l(:,jb-ista_e+1,ik) = Phifftr_ball(:,neordr(ib,ik))
          end if
       end do
       deallocate(Phifftr_ball)
#else
       allocate(Phifftr_trans(lsize*kimg,np_e)); Phifftr_trans = 0.d0
       allocate(Phifftr_wk(lsize*kimg))

       do ib = 1, neg
          if(neordr(ib,ik)==ib) cycle
          ibr = nrvf_ordr(ib,ik)
          nbs = (ib-1)/NB + 1
          nbsr = (ibr-1)/NB + 1
          if(lrank(nbs) == lrank(nbsr)) cycle
          if(lrank(nbs) == myrank_e) Phifftr_wk = Phifftr_l(:,map_z(ib),ik)
          call mpi_bcast(Phifftr_wk,lsize*kimg,mpi_double_precision,lrank(nbs),mpi_kg_world,ierr)
          if(ista_e <= tr_neg(ibr) .and. tr_neg(ibr) <= iend_e) then
             jb = tr_neg(ibr)
             Phifftr_trans(:,jb-ista_e+1) = Phifftr_wk
          end if
       end do
       do ib = 1, neg
          if(neordr(ib,ik)==ib) cycle
          ibr = nrvf_ordr(ib,ik)
          nbs = (ib-1)/NB + 1
          nbsr = (ibr-1)/NB + 1
          if(lrank(nbs) == lrank(nbsr) .and. lrank(nbs) == myrank_e) then
             if(ista_e <= tr_neg(ibr) .and. tr_neg(ibr) <= iend_e) then
                jb = tr_neg(ibr)
                Phifftr_trans(:,jb-ista_e+1) = Phifftr_l(:,map_z(ib),ik)
             end if
          end if
       end do

!!$       do ib = 1, neg
!!$          ibr = neordr(ib,ik)
!!$          if(ibr==ib) cycle
!!$          if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
!!$             jb = tr_neg(ib)
!!$             Phifftr_l(:,jb-ista_e+1,ik) = Phifftr_trans(:,jb-ista_e+1)
!!$          end if
!!$       end do
       do jb = 1, np_e
          ib = jb + ista_e - 1
          ibr = neg_g_all(ib)
          if(neordr(ibr,ik) == ibr) cycle
          Phifftr_l(:,jb,ik) = Phifftr_trans(:,jb)
       end do

       deallocate(Phifftr_wk)
       deallocate(Phifftr_trans)
#endif
#endif
#endif
    end subroutine sort_Phifftr_l

    subroutine sort_status_saved_phifftr(ik)
      integer, intent(in) :: ik
      integer                                      :: status_wk
      integer    :: ib
#ifdef SAVE_FFT_TIMES00
      integer,       allocatable, dimension(:)     :: status_ball   !d(neg)
#else
      integer    :: nbs,jb
#ifdef SAVE_FFT_TIMES01
      integer,       allocatable, dimension(:)     :: status_ball   !d(neg)
#else
      integer,       allocatable, dimension(:)     :: status_trans   !d(np_e)
      integer :: ibr, nbsr
#endif
#endif

#ifdef SAVE_FFT_TIMES00
      allocate(status_ball(neg)); status_ball = 0
      do ib = 1, np_e
         status_ball(neg_g(ib)) = status_saved_phifftr(ib,ik)
      end do
      call mpi_allreduce(MPI_IN_PLACE,status_ball,neg,MPI_INTEGER,MPI_SUM,mpi_kg_world,ierr)
      do ib = 1, neg
         if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
            jb = tr_neg(ib)
            if(neordr(jb,ik) /= jb) then
               status_saved_phifftr(jb-ista_e+1,ik) = status_ball(neordr(jb,ik))
            end if
         end if
      end do
      deallocate(status_ball)
#else
#ifdef SAVE_FFT_TIMES01
      allocate(status_ball(neg)); status_ball = 0
      do ib = 1, neg
         if(neordr(ib,ik) == ib) cycle
         nbs = (ib-1)/NB + 1
         if(lrank(nbs) == myrank_e) status_wk = status_saved_phifftr(map_z(ib),ik)
         call mpi_bcast(status_wk,1,mpi_integer,lrank(nbs),mpi_kg_world,ierr)
         status_ball(ib) = status_wk
      end do
      do ib = 1, neg
         if(neordr(ib,ik) == ib) cycle
         if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
            jb = tr_neg(ib)
            status_saved_phifftr(jb-ista_e+1,ik) = status_ball(neordr(ib,ik))
         end if
      end do
      deallocate(status_ball)
#else
      allocate(status_trans(np_e))
      do ib = 1, neg
          if(neordr(ib,ik)==ib) cycle
          ibr = nrvf_ordr(ib,ik)
          nbs = (ib-1)/NB + 1
          nbsr = (ibr-1)/NB + 1
          if(lrank(nbs) == lrank(nbsr)) cycle
          if(lrank(nbs) == myrank_e) status_wk = status_saved_phifftr(map_z(ib),ik)
          call mpi_bcast(status_wk,1,mpi_integer,lrank(nbs),mpi_kg_world,ierr)
          if(ista_e <= tr_neg(ibr) .and. tr_neg(ibr) <= iend_e) then
             jb = tr_neg(ibr)
             status_trans(jb-ista_e+1) = status_wk
!!$          status_ball(ib) = status_wk
          end if
       end do
       do ib = 1, neg
          if(neordr(ib,ik)==ib) cycle
          ibr = nrvf_ordr(ib,ik)
          nbs = (ib-1)/NB + 1
          nbsr = (ibr-1)/NB + 1
          if(lrank(nbs) == lrank(nbsr) .and. lrank(nbs) == myrank_e) then
             if(ista_e <= tr_neg(ibr) .and. tr_neg(ibr) <= iend_e) then
                jb = tr_neg(ibr)
                status_trans(jb-ista_e+1) = status_saved_phifftr(map_z(ib),ik)
             end if
          end if
       end do
!!$       do ib = 1, neg
!!$          ibr = neordr(ib,ik)
!!$          if(ibr==ib) cycle
!!$          if(ista_e <= tr_neg(ib) .and. tr_neg(ib) <= iend_e) then
!!$             jb = tr_neg(ib)
!!$             status_saved_phifftr(jb-ista_e+1,ik) = status_ball(ibr)
!!$          end if
!!$       end do
       do jb = 1, np_e
          ib = jb + ista_e - 1
          ibr = neg_g_all(ib)
          if(neordr(ibr,ik) == ibr) cycle
          status_saved_phifftr(jb,ik) = status_trans(jb)
       end do
       deallocate(status_trans)
#endif
#endif

    end subroutine sort_status_saved_phifftr
#endif

    subroutine cp_eigen_values_for_af()
                                                 __TIMER_SUB_START(606)
                                                 __TIMER_DO_START(634)
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle                        ! MPI
         do ib = 1, np_e                                           ! MPI
            eko_l(ib,ik+af) = eko_l(ib,ik)
         enddo
      enddo
                                                 __TIMER_DO_STOP(634)
                                                 __TIMER_SUB_STOP(606)
    end subroutine cp_eigen_values_for_af

    subroutine expand_neordr_and_nrvf_ordr()
      integer :: ik
                                                 __TIMER_SUB_START(607)
                                                 __TIMER_DO_START(635)
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle                        ! MPI
         neordr(1:neg,ik+af) = neordr(1:neg,ik)
         nrvf_ordr(1:neg,ik+af) = nrvf_ordr(1:neg,ik)
      end do
                                                 __TIMER_DO_STOP(635)
                                                 __TIMER_SUB_STOP(607)
    end subroutine expand_neordr_and_nrvf_ordr

#ifdef _NO_HEAP_SORT_EIGENVALUES_
    subroutine bubble_sorting(n,eig,iord)
      integer, intent(in) :: n
      real(kind=DP), intent(in) :: eig(n)
      integer, intent(inout) :: iord(n)
      real(kind=DP), parameter :: delta = 1.d-12
      integer :: ib, jb, ibo, jbo

                                                 __TIMER_DO_START(617)
      do ib = 1, neg-1
         do jb = ib+1, neg
            ibo = iord(ib)
            jbo = iord(jb)
            if(eig(jbo)  < eig(ibo)-delta) then        ! MPI
               iord(jb) = ibo
               iord(ib) = jbo
            end if
         end do
      end do
                                                 __TIMER_DO_STOP(617)
     end subroutine bubble_sorting
#else
    subroutine heap_sorting(n,eig,iord)
      integer, intent(in) :: n
      real(kind=DP), intent(inout) :: eig(n)
      integer, intent(inout) :: iord(n)

      integer :: i,j,k,itmp,it
      real(kind=DP) :: rt
                                                 __TIMER_SUB_START(605)
      !*** initialization: heapfy eig and iord ***
                                                 __TIMER_DO_START(632)
      do i=2,n
         rt = eig(i)
         it = iord(i)
         j=i
10       itmp=j/2
         if(eig(itmp).ge.rt) go to 20
         eig(j) = eig(itmp)
         iord(j) = iord(itmp)
         j=itmp
         if(j.gt.1) go to 10
20       eig(j) = rt
         iord(j) = it
      end do
                                                 __TIMER_DO_STOP(632)
      !*** to be fully sort ***
                                                 __TIMER_DO_START(633)
      do k=n-1,1,-1
         rt = eig(1)
         eig(1) = eig(k+1)
         eig(k+1) = rt
         rt = eig(1)
         it = iord(1)
         iord(1) = iord(k+1)
         iord(k+1) = it
         it = iord(1)
         j=1
         itmp=2
30       if(itmp.gt.k) go to 40
         if(itmp.lt.k) then
            if(eig(itmp+1).gt.eig(itmp)) itmp=itmp+1
         end if
         if(rt.ge.eig(itmp)) go to 40
         eig(j) = eig(itmp)
         iord(j) = iord(itmp)
         j=itmp
         itmp=j*2
         go to 30
40       eig(j) = rt
         iord(j) = it
      end do
                                                 __TIMER_DO_STOP(633)
                                                 __TIMER_SUB_STOP(605)
    end subroutine heap_sorting
#endif
  end subroutine m_ES_sort_eigen_values_3D

! === DEBUG by tkato 2013/02/24 ================================================
! In RTP-TDDFT, such a sorting brings on difference to ORG_Parallel in current.
! The order of band elements on zaj_l is not match to that on occup_l.
! But, this sort is necessary for other case.
!subroutine m_ES_energy_eigen_values_3D(nfout)
 subroutine m_ES_energy_eigen_values_3D(nfout,no_sort)
! ==============================================================================
    integer, intent(in) :: nfout
! === DEBUG by tkato 2013/02/24 ================================================
! In RTP-TDDFT, such a sorting brings on difference to ORG_Parallel in current.
! The order of band elements on zaj_l is not match to that on occup_l.
! But, this sort is necessary for other case.
    logical, optional, intent(in) :: no_sort
! ==============================================================================

    integer             :: is, ik, ipri0
    real(kind=DP), allocatable, dimension(:) :: ekin_l
    real(kind=DP), allocatable, dimension(:) :: afft_l
    integer             :: lsize, ierr, ii
#ifdef MPI_FFTW
    integer(C_INTPTR_T)  :: local_n, local_n_offset, alloc_local, lx, ly, lz, mx,my,mz
    lx = fft_box_size_WF(1,0)
    ly = fft_box_size_WF(2,0)
    lz = fft_box_size_WF(3,0)
    if(kimg==2) then
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lx,mpi_ke_world,local_n,local_n_offset)
    else
      alloc_local = fftw_mpi_local_size_3d(ly,lz,lx/2,mpi_ke_world,local_n,local_n_offset)
    endif
#endif

    call m_FFT_alloc_WF_work()
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
    allocate(afft_l(lsize*2), stat=ierr)
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#ifdef MPI_FFTW
    if(sw_mpi_fftw==ON) then
      lsize = local_n*lx*lz
    endif
#endif
    allocate(afft_l(lsize*kimg), stat=ierr)
#endif
     if(ierr /= 0) then
        write(nfout,*)' m_ES_energy_eigen_values_3D : Not allocated afft_l array'
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 201, ierr)
     endif
    allocate(ekin_l(1:maxval(np_g1k(:))))

    do is = 1, nspin, af+1
#ifdef MPI_FFTW
       if(sw_mpi_fftw==ON) then
         call m_ES_Vlocal_in_Rspace_mpifftw(is,afft_l,lsize,1,OFF)      ! (ptfft1) vlhxc_l->afft
       else
         call m_ES_Vlocal_in_Rspace_3D(is,afft_l,lsize,1,OFF)
       endif
#else
       call m_ES_Vlocal_in_Rspace_3D(is,afft_l,lsize,1,OFF)
#endif
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle                   ! MPI
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l)
#ifdef MPI_FFTW
          if(sw_mpi_fftw==ON) then
            call m_ES_eigen_values_for_each_k_mpifftw(is,ik,ekin_l,afft_l,lsize)
          else
            call m_ES_eigen_values_for_each_k_3D(is,ik,ekin_l,afft_l,lsize)
          endif
#else
          call m_ES_eigen_values_for_each_k_3D(is,ik,ekin_l,afft_l,lsize)
#endif
       end do
    end do
! === DEBUG by tkato 2013/02/24 ================================================
! In RTP-TDDFT, such a sorting brings on difference to ORG_Parallel in current.
! The order of band elements on zaj_l is not match to that on occup_l.
! But, this sort is necessary for other case.
!   call m_ES_sort_eigen_values_3D()
    if(.not. present(no_sort)) then
       call m_ES_sort_eigen_values_3D()
    end if
! ==============================================================================

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko_3D(nfout,mode=SCF)
    deallocate(ekin_l)
    deallocate(afft_l)
    call m_FFT_dealloc_WF_work()
  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
  end subroutine m_ES_energy_eigen_values_3D

  subroutine m_ES_wd_eko_3D(nfout,mode)

    integer, intent(in) ::                  nfout
    integer, intent(in) ::                  mode
    real(kind=DP), pointer, dimension(:) :: eko, eko_t ! d(neg) MPI

    integer :: ik, ib, ikp
    integer :: j, k
    real(kind=DP) :: c1, vkxyz_wk(3)

    allocate(eko(neg)); allocate(eko_t(neg))           ! MPI
    if(printable) write(nfout,*) '=== energy_eigen_values ==='
    ikp = 0
    if(mode == EK) ikp = nk_in_the_process - 1

    do ik = 1, kv3, af+1
       if ( sw_band_unfolding == ON .and. mode == EK ) then
          Do j=1, 3
             c1 = 0.0d0
             Do k=1, 3
                c1 = c1 +altv_refcell(k,j) *vkxyz(ik,k,CARTS)
             End Do
             vkxyz_wk(j) = c1 /PAI2
          End Do
       else
          vkxyz_wk = vkxyz(ik,1:3,BUCS)
       endif

       if(map_k(ik) /= myrank_k) cycle                 ! MPI
       eko_t = 0                                       ! MPI
       do ib = 1, neg                                  ! MPI
          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik) ! MPI
       end do                                          ! MPI
! === DEBUG by tkato 2013/11/19 ================================================
!      call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_kg_world,ierr) ! MPI
! ==============================================================================
       if(mode == EK .or. ik == 1) then
          if(printable) write(nfout,'(" ik = ",i7," ( ",3f10.6," )",/99(4f18.10,/))')&
               &ik+ikp,vkxyz(ik,1:3,BUCS), (eko(neordr(ib,ik)),ib=1,neg-num_extra_bands)
!!$          write(nfout,'(" ik = ",i4,99(4f18.10,/))')&
!!$               &ik,(eko(neordr(ib,ik)),ib=1,neg)
          if(iprieigenvalue >= 2 .and. num_extra_bands >= 1) &
               & write(nfout,'(" -- extra_bands --",/99(4f18.10,/))') &
               & (eko(neordr(ib,ik)),ib=neg-num_extra_bands+1, neg)
       else
          if(printable) write(nfout,'(" ik = ",i7," (",3f10.6," )",/99(10f8.4,/))')&
               & ik,vkxyz_wk(1:3), (eko(neordr(ib,ik)),ib=1,neg)
       endif
    end do
    deallocate(eko); deallocate(eko_t)                 ! MPI
  end subroutine m_ES_wd_eko_3D
!----
  subroutine m_ES_energy_eigen_values_ext_3D(nfout)
    integer, intent(in) :: nfout

    integer             :: is, ik, ipri0
    real(kind=DP), allocatable, dimension(:) :: ekin_l
    real(kind=DP), allocatable, dimension(:) :: afft_l
    integer             :: lsize, ierr, ii

    call m_FFT_alloc_WF_work()
#ifdef FFT_3D_DIVISION
    lsize = fft_X_x_nel*fft_X_y_nel*fft_X_z_nel
#else
    lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
#endif
    allocate(afft_l(lsize*kimg), stat=ierr)
    allocate(ekin_l(1:maxval(np_g1k(:))))

    do is = 1, nspin, af+1
       call m_ES_Vlocal_in_Rspace_3D(is,afft_l,lsize,1,OFF)
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle                   ! MPI
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin_l)

          if(ipri >= 3) then
             write(nfout,'(" -- ekin <<m_ES_eigen_values_ext>> -- ik = ",i8)') ik
             write(nfout,'(8f8.4)') (ekin_l(ipri0),ipri0=1,iba(ik))
          end if

          call m_ES_eigen_values_for_each_kex_3D(nfout,is,ik,ekin_l,afft_l,lsize)
       end do
    end do
    call m_ES_sort_eigen_values_3D()

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko2(nfout,mode=SCF)

    deallocate(afft_l)
    deallocate(ekin_l)
    call m_FFT_dealloc_WF_work()
  contains
    subroutine get_ipri0(ipri_in, ipri_out)
      integer, intent(in)  :: ipri_in
      integer, intent(out) :: ipri_out
      if(npes > 1) then
         if(mype == 0) ipri_out = ipri_in
         call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
      else
         ipri_out = ipri_in
      end if
    end subroutine get_ipri0
  end subroutine m_ES_energy_eigen_values_ext_3D

  subroutine m_ES_eigen_values_for_each_kex_3D(nfout,is,ik,ekin_l,afft_l,lsize)
    integer,       intent(in)                  :: nfout,is,ik,lsize    ! is: spin
    real(kind=DP), intent(in), dimension(np_g1k(ik)) :: ekin_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(in), dimension(2*kimg) :: afft_l
#else
    real(kind=DP), intent(in), dimension(lsize*kimg) :: afft_l
#endif

    real(kind=DP), allocatable,dimension(:,:) :: bfft_l, bfft_l_in

    integer       :: ib, ib1, ib2, ibsize, ibesize
    real(kind=DP),dimension(np_e) :: eg, eko
    real(kind=DP), allocatable, dimension(:) :: wk_mpi
    real(kind=DP)                    :: prd

    integer,save  :: id_sname = -1
    call tstatc0_begin('energy_eigen_values ', id_sname,1)

    ibsize = 1
    if (nblocksize_fftw_is_given) then
       ibsize = nblocksize_fftw
       if (ibsize < 1) ibsize = 1
    endif

#ifdef FFT_3D_DIVISION
    allocate(bfft_l(lsize*2,ibsize) ,stat=ierr)
#else
    allocate(bfft_l(lsize*kimg,ibsize) ,stat=ierr)
    allocate(bfft_l_in(lsize*kimg,ibsize) ,stat=ierr)
#endif
     if (ierr /= 0) then
        write(nfout,*)' m_ES_eigen_values_for_each_k_3D :  Not allocate '
        call flush(nfout)
        call mpi_abort(mpi_comm_world, 178, ierr)
     endif

    eko(:) = 0.0d0

    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part

    do ib1 = 1, np_e, ibsize   ! MPI
       ib2 = min(ib1+ibsize-1,np_e)
       ibesize = ib2 - ib1 + 1

#ifdef __TIMER_COMM__
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l, 0)
#else
       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l)
!       call m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l_in,bfft_l)
#endif
#ifdef FFT_3D_DIVISION
       call m_FFT_W_Vlocal_W_3DIV_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
#else
       if (sw_fft_xzy > 0) then
          call m_FFT_W_Vlocal_W_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
       else
          call m_FFT_W_Vlocal_W_XYZ_3D(ELECTRON,afft_l,bfft_l,lsize,ibesize,eg(ib1)) ! (eigens) --> eg
       end if
#endif
    end do

    deallocate(bfft_l)
    deallocate(bfft_l_in)

    prd = product(fft_box_size_WF(1:3,1))

    eko(:) = eko(:) + eg(:) / prd

    call W_Vnonlocal_W(ik)

    if (nrank_g > 1) then
       allocate(wk_mpi(1:np_e))
       call mpi_allreduce(eko, wk_mpi, np_e, mpi_double_precision, mpi_sum, mpi_ke_world, ierr)
       eko_l(1:np_e,ik) = wk_mpi(1:np_e)
       deallocate(wk_mpi)
    else
       eko_l(1:np_e,ik) = eko(1:np_e)
    end if

    call tstatc0_end(id_sname)
  contains
    subroutine W_T_W()
      integer :: ib, i, ip, ibt

!!$      eko_l(1:np_e,ik) = 0.d0                   ! MPI

      if(kimg==1) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
!         do ib = neg_previous+1, neg                   ! MPI
      do ib = 1, np_e
         if(ista_e+ib-1 < neg_previous+1) cycle
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, np_g1k(ik)
                  eko(ibt) = eko(ibt) + ekin_l(i)*zaj_l(i,ibt,ik,1)**2
               end do
               if(k_symmetry(ik) == GAMMA) eko(ibt) = eko(ibt)*2.d0
            end if
         end do
      else if(kimg==2) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
!!$         do ib = 1, np_e                        ! MPI
!         do ib = neg_previous+1, neg
      do ib = 1, np_e
         if(ista_e+ib-1 < neg_previous+1) cycle
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, np_g1k(ik)
                  eko(ibt) = eko(ibt) + ekin_l(i)*(zaj_l(i,ibt,ik,1)**2+zaj_l(i,ibt,ik,2)**2)
               end do
               if(k_symmetry(ik) == GAMMA) eko(ibt) = eko(ibt)*2.d0
            end if
         end do
      end if

      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_T_W) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko(ib),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)

      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib, ibo
      real(kind=DP) :: fac

      integer,save :: id_sname = -1
                                                 __TIMER_SUB_START(603)
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)

      call m_ES_alloc_fsr_l_2d(np_e, nlmta)
      call m_ES_gather_f_3d_to_2d(fsr_l, fsr_l_2D, ik)
      if(k_symmetry(ik) /= GAMMA) then
         call m_ES_alloc_fsi_l_2d(np_e, nlmta)
         call m_ES_gather_f_3d_to_2d(fsi_l, fsi_l_2D, ik)
      endif

                                                 __TIMER_DO_START(610)
      do ia = 1, natm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            do lmt2 = lmt1, ilmt(it)
               q = lmta(lmt2,ia)
               fac   = 2.d0 * iwei(ia)
               if(lmt1 == lmt2) fac = iwei(ia)
!!$               fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
               if(ipaw(it).eq.0) then
                  fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
               else
                  fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
               end if
                                                 __TIMER_DO_START(611)
               if(k_symmetry(ik) == GAMMA) then
                  do ib = 1, np_e                                  ! MPI
                     ibo = ista_e + ista_e*(ib-1)
                     if(nrvf_ordr(ibo,ik) <= neg_previous) cycle
                     eko(ib) = eko(ib) &
                          & + fac*(fsr_l_2D(ib,p)*fsr_l_2D(ib,q))
                  end do
               else
                  do ib = 1, np_e                                  ! MPI
                     ibo = ista_e + ista_e*(ib-1)
                     if(nrvf_ordr(ibo,ik) <= neg_previous) cycle
                     eko(ib) = eko(ib) &
                          & + fac*(fsr_l_2D(ib,p)*fsr_l_2D(ib,q)&
                          &      + fsi_l_2D(ib,p)*fsi_l_2D(ib,q))
                  end do
               end if
                                                 __TIMER_DO_STOP(611)
            end do
         end do
      end do
                                                 __TIMER_DO_STOP(610)
      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_Vnonlocal_W) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif

      call m_ES_dealloc_fsr_l_2d()
      call m_ES_dealloc_fsi_l_2d()

      call tstatc0_end(id_sname)
                                                 __TIMER_SUB_STOP(603)
    end subroutine W_Vnonlocal_W
  end subroutine m_ES_eigen_values_for_each_kex_3D


!---
!!$#ifdef __TIMER_COMM__
!!$  subroutine m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l, kukanNo)
!!$#else
!!$  subroutine m_ES_WF_in_Rspace_3D(ik,ib1,ib2,ibsize,lsize,bfft_l)
!!$#endif
!!$    use m_Parallelization,     only : wf_fft_scnt, wf_fft_rcnt &
!!$   &                                , wf_fft_recv &
!!$   &                                , wf_fft_index, wf_fft_dist &
!!$   &                                , wf_fft_maxrecv, wf_fft_maxsend
!!$    integer, intent(in)                           :: ik, ib1,ib2, ibsize,lsize
!!$#ifdef __TIMER_COMM__
!!$    integer, intent(in)                           :: kukanNo
!!$#endif
!!$#ifdef FFT_3D_DIVISION
!!$    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
!!$#else
!!$    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
!!$#endif
!!$!   real(kind=DP), dimension(nfft,ibsize) :: bfft
!!$! === DEBUG by tkato 2013/08/28 ================================================
!!$!   real(kind=DP), dimension(1,1) :: bfft
!!$! ==============================================================================
!!$    integer :: i, j, k, ii, jj, iesize, iadd, ierr
!!$    integer, dimension(0:nrank_g-1) ::req_r,req_s
!!$    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
!!$    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
!!$    integer :: icnt_send, icnt_recv, lrank
!!$    integer, parameter :: itag = 21
!!$    integer, save :: savesize = 0, savesend=0, saverecv=0
!!$!!!ifdef FFT_USE_SSL2
!!$    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!$!!!endif
!!$    integer,save  :: id_sname = -1
!!$#ifdef SAVE_FFT_TIMES
!!$    logical :: fft_stored
!!$    integer,save :: id_sname2 = -1
!!$#endif
!!$
!!$! === DEBUG by tkato 2012/06/05 ================================================
!!$#ifdef USE_ALLTOALLV
!!$    integer, allocatable, dimension(:) :: sdsp, rdsp
!!$#endif
!!$! ==============================================================================
!!$
!!$#ifdef SAVE_FFT_TIMES
!!$    if(sw_save_fft == ON) then
!!$       fft_stored = .true.
!!$       do jj = ib1, ib2
!!$          if(status_saved_phifftr(jj,ik) /= STORED_AND_NEW) fft_stored = .false.
!!$       end do
!!$    else
!!$       fft_stored = .false.
!!$    end if
!!$!!$!!x!!$    if(status_saved_phifftr(map_z(ib),ik) == STORED_AND_NEW) then
!!$    if(fft_stored) then
!!$       call tstatc0_begin('m_ES_WF_in_Rspace_3D(2) ',id_sname2)
!!$       if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_fftr(stored) --> bfft")')
!!$       iesize = ib2 - ib1 + 1
!!$       bfft_l(:,1:iesize) = Phifftr_l(:,ib1:ib2,ik)
!!$       call tstatc0_end(id_sname2)
!!$    else
!!$#endif
!!$
!!$#ifdef __TIMER_SUB__
!!$    call timer_sta(304)
!!$#endif
!!$
!!$    call tstatc0_begin('m_ES_WF_in_Rspace_3D(1) ',id_sname)
!!$
!!$    if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_l --(FFT)--> bfft")')
!!$
!!$    iesize = ib2 - ib1 + 1
!!$
!!$!   allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!!$!   allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!!$!   if ((iesize /= savesize) .or. (wf_fft_maxsend(ik) /= savesend)) then
!!$       if (allocated(sendbuf)) deallocate(sendbuf)
!!$       allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!!$!      savesend = wf_fft_maxsend(ik)
!!$!   end if
!!$!   if ((iesize /= savesize) .or. (wf_fft_maxrecv(ik) /= savesend)) then
!!$       if (allocated(recvbuf)) deallocate(recvbuf)
!!$       allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!!$!      saverecv = wf_fft_maxrecv(ik)
!!$!   end if
!!$    savesize = iesize
!!$    sendbuf = 0.0d0
!!$    recvbuf = 0.0d0
!!$
!!$#ifdef __TIMER_DO__
!!$  call timer_sta(313)
!!$#endif
!!$    if(k_symmetry(ik) == GAMMA) then
!!$       if (kimg == 1) then
!!$          do jj = 1, iesize
!!$             do ii = 1, np_g1k(ik)
!!$                sendbuf(iesize*(wf_fft_index(ii*2-1,ik)-1)+1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
!!$                sendbuf(iesize*(wf_fft_index(ii*2  ,ik)-1)+1,wf_fft_dist(ii*2  ,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
!!$             enddo
!!$          enddo
!!$       else
!!$!OCL NORECURRENCE
!!$          do jj = 1, iesize
!!$             do ii = 1, np_g1k(ik)
!!$                iadd = iesize*2*(wf_fft_index(ii*2-1,ik)-1)+jj*2
!!$                sendbuf(iadd-1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
!!$                sendbuf(iadd,  wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
!!$                iadd = iesize*2*(wf_fft_index(ii*2,ik)-1)+jj*2
!!$                sendbuf(iadd-1,wf_fft_dist(ii*2  ,ik)) =  zaj_l(ii,ib1+jj-1,ik,1)
!!$                sendbuf(iadd,  wf_fft_dist(ii*2  ,ik)) = -zaj_l(ii,ib1+jj-1,ik,2)
!!$             enddo
!!$          enddo
!!$       endif
!!$    else
!!$       if (kimg == 1) then
!!$!OCL NORECURRENCE
!!$          do jj = 1, iesize
!!$             do ii = 1, np_g1k(ik)
!!$                iadd = iesize*(wf_fft_index(ii,ik)-1)+jj
!!$                sendbuf(iadd,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
!!$             enddo
!!$          enddo
!!$       else
!!$!OCL NORECURRENCE
!!$          do jj = 1, iesize
!!$             do ii = 1, np_g1k(ik)
!!$                iadd = iesize*2*(wf_fft_index(ii,ik)-1)+jj*2
!!$                sendbuf(iadd-1,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
!!$                sendbuf(iadd  ,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
!!$             enddo
!!$          enddo
!!$       endif
!!$    endif
!!$#ifdef __TIMER_DO__
!!$  call timer_end(313)
!!$#endif
!!$
!!$#ifndef USE_ALLTOALLV
!!$    icnt_recv = 0
!!$#ifdef __TIMER_COMM__
!!$  call timer_barrier(mpi_ke_world)
!!$  call timer_sta(314)
!!$#endif
!!$    do i = 0, nrank_g - 1
!!$       if (wf_fft_rcnt(i,ik) /= 0) then
!!$          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt(i,ik)*kimg*iesize, mpi_double_precision, &
!!$      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
!!$           if (ierr /= 0) then
!!$              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_irecv error'
!!$              call flush(nfout)
!!$              call mpi_abort(mpi_comm_world,166,ierr)
!!$           endif
!!$          icnt_recv = icnt_recv + 1
!!$       endif
!!$    enddo
!!$    icnt_send = 0
!!$    do i = 0, nrank_g - 1
!!$       if (wf_fft_scnt(i,ik) /= 0) then
!!$          call mpi_isend(sendbuf(1,i), wf_fft_scnt(i,ik)*kimg*iesize, mpi_double_precision, &
!!$      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
!!$           if (ierr /= 0) then
!!$              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_isend error'
!!$              call flush(nfout)
!!$              call mpi_abort(mpi_comm_world,167,ierr)
!!$           endif
!!$          icnt_send = icnt_send + 1
!!$       endif
!!$    enddo
!!$    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
!!$     if (ierr /= 0) then
!!$        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
!!$        call flush(nfout)
!!$        call mpi_abort(mpi_comm_world,168,ierr)
!!$     endif
!!$    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
!!$#ifdef __TIMER_COMM__
!!$  call timer_end(314)
!!$#endif
!!$     if (ierr /= 0) then
!!$        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
!!$        call flush(nfout)
!!$        call mpi_abort(mpi_comm_world,169,ierr)
!!$     endif
!!$#else
!!$#ifdef __TIMER_COMM__
!!$  call timer_barrier(mpi_ke_world)
!!$  call timer_sta(329)
!!$#endif
!!$! === DEBUG by tkato 2012/06/05 ================================================
!!$!      integer, allocatable, dimension(:) :: sdsp, rdsp
!!$! ==============================================================================
!!$       allocate(sdsp(0:nrank_g-1), stat=ierr)
!!$       allocate(rdsp(0:nrank_g-1), stat=ierr)
!!$       do i = 0, nrank_g - 1
!!$          sdsp(i)=wf_fft_maxsend(ik)*kimg*iesize*i
!!$          rdsp(i)=wf_fft_maxrecv(ik)*kimg*iesize*i
!!$       enddo
!!$       call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt(:,ik)*kimg*iesize, sdsp,&
!!$      &   mpi_double_precision, recvbuf, wf_fft_rcnt(:,ik)*kimg*iesize, rdsp,&
!!$      &   mpi_double_precision, mpi_ke_world, ierr )
!!$       if (ierr /= 0) then
!!$          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
!!$          call flush(nfout)
!!$          call mpi_abort(mpi_comm_world, 170, ierr)
!!$       endif
!!$       deallocate(sdsp)
!!$       deallocate(rdsp)
!!$#ifdef __TIMER_COMM__
!!$  call timer_end(329)
!!$#endif
!!$#endif
!!$
!!$    bfft_l = 0.0d0
!!$#ifdef __TIMER_DO__
!!$  call timer_sta(315)
!!$#endif
!!$#ifdef FFT_3D_DIVISION
!!$    if (kimg == 1) then
!!$!OCL NORECURRENCE
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             do k = 1, iesize
!!$!               bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
!!$                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*(j-1)+k,i)
!!$                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = 0.0d0
!!$             enddo
!!$          enddo
!!$       enddo
!!$    else
!!$!OCL NORECURRENCE
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             do k = 1, iesize
!!$                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
!!$                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    endif
!!$#else
!!$#ifdef FFT_USE_SSL2
!!$    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
!!$    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
!!$    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
!!$    nxp = nx
!!$    if (kimg == 1) then
!!$!OCL NORECURRENCE
!!$!OCL NOFLTLD
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             iadd = wf_fft_recv(j,ik,i)
!!$             iz = (iadd-1)/(nx*ny)+1
!!$             nn = mod(iadd,nx*ny)
!!$             if (nn==0) then
!!$                iy = ny
!!$             else
!!$                iy = (nn-1)/nx+1
!!$             end if
!!$             ix = mod(nn,nx)
!!$             if(ix==0) ix = nx
!!$             do k = 1, iesize
!!$                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    else
!!$!OCL NORECURRENCE
!!$!OCL NOFLTLD
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             iadd = wf_fft_recv(j,ik,i)
!!$             iz = (iadd-1)/(nx*ny)+1
!!$             nn = mod(iadd,nx*ny)
!!$             if (nn==0) then
!!$                iy = ny
!!$             else
!!$                iy = (nn-1)/nx+1
!!$             end if
!!$             ix = mod(nn,nx)
!!$             if(ix==0) ix = nx
!!$             do k = 1, iesize
!!$                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
!!$                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    endif
!!$#else
!!$    if (kimg == 1) then
!!$!OCL NORECURRENCE
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             do k = 1, iesize
!!$                bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    else
!!$!OCL NORECURRENCE
!!$       do i = 0, nrank_g - 1
!!$          do j = 1, wf_fft_rcnt(i,ik)
!!$             do k = 1, iesize
!!$                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
!!$                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
!!$             enddo
!!$          enddo
!!$       enddo
!!$    endif
!!$#endif
!!$#endif
!!$#ifdef __TIMER_DO__
!!$  call timer_end(315)
!!$#endif
!!$!!  deallocate(sendbuf)
!!$!!  deallocate(recvbuf)
!!$
!!$#ifdef __TIMER_COMM__
!!$    if(kukanNo==3) call timer_sta(327)
!!$    if(kukanNo==6) call timer_sta(636)
!!$    if(kukanNo==7) call timer_sta(691)
!!$    if(kukanNo==8) call timer_sta(961)
!!$#endif
!!$#ifdef FFT_3D_DIVISION
!!$    call m_FFT_Inverse_3DIV_3D(nfout, bfft_l, lsize, iesize)
!!$#else
!!$    if (sw_fft_xzy > 0) then
!!$       call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
!!$    else
!!$       call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
!!$    end if
!!$#endif
!!$#ifdef __TIMER_COMM__
!!$    if(kukanNo==3) call timer_end(327)
!!$    if(kukanNo==6) call timer_end(636)
!!$    if(kukanNo==7) call timer_end(691)
!!$    if(kukanNo==8) call timer_end(961)
!!$#endif
!!$
!!$    if(ipri >= 2) then
!!$       if(ik <=  2 .and. ib1 <= 1) then
!!$          write(6,'(" ! bfft R-space ik = ",i3," ib1 = ",i3," <<m_ES_WF_in_Rspace>>")') ik, ib1
!!$! === DEBUG by tkato 2013/08/28 ================================================
!!$!         write(6,'(8f8.4)') (bfft(i,1),i=1,120)
!!$          write(6,'(8f8.4)') (bfft_l(i,1),i=1,120)
!!$! ==============================================================================
!!$       end if
!!$    end if
!!$
!!$#ifdef __TIMER_SUB__
!!$    call timer_end(304)
!!$#endif
!!$
!!$#ifdef SAVE_FFT_TIMES
!!$    if(sw_save_fft == ON) then
!!$       iesize = ib2 - ib1 + 1
!!$       Phifftr_l(:,ib1:ib2,ik) = bfft_l(:,1:iesize)
!!$       status_saved_phifftr(ib1:ib2,ik) = STORED_AND_NEW
!!$    end if
!!$#endif
!!$    call tstatc0_end(id_sname)
!!$#ifdef SAVE_FFT_TIMES
!!$    end if
!!$#endif
!!$
!!$  end subroutine m_ES_WF_in_Rspace_3D

  subroutine m_ES_WF_2D0(ik,bfft_l,ib2,ib1,ibsize,lsize,inverse_or_direct)
    integer, intent(in) :: ik,ib2,ib1,ibsize,lsize,inverse_or_direct
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
    call m_ES_WF_2D_psi(ik,bfft_l,zaj_l,ib2,ib1,ista_k,iend_k,ibsize,lsize,inverse_or_direct)
  end subroutine m_ES_WF_2D0

  subroutine m_ES_WF_2D1(ik,bfft_l,ib2,ib1,ibsize,lsize,inverse_or_direct,bff)
    integer, intent(in) :: ik,ib2,ib1,ibsize,lsize,inverse_or_direct
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
    real(kind=DP), intent(out), dimension(nfft,ibsize) :: bff
    call m_ES_WF_2D_psi(ik,bfft_l,zaj_l,ib2,ib1,ista_k,iend_k,ibsize,lsize,inverse_or_direct,bff)
  end subroutine m_ES_WF_2D1

  subroutine m_ES_WF_2D_psi(ik,bfft_l,psi_l,ib2,ib1,k1,k2,ibsize,lsize,inverse_or_direct,bff)
    integer, intent(in) :: ik,ib2,ib1,k1,k2,ibsize,lsize,inverse_or_direct
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
    real(kind=DP),intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(out), dimension(nfft,ibsize), optional :: bff
    integer :: ib,ig, ir, gk,gj,gl,gi,iadd,iadd1,gk_,gj_,gl_,nl,nm,nn
    real(kind=DP), allocatable, dimension(:,:), save :: zajbuf
    real(kind=DP), allocatable, dimension(:), save :: bfft
    integer :: nx, ny, nz, iesize,i
    logical :: distrib = .true.
    integer(kind=4),dimension(0:nrank_g-1) :: recvcnt, recvdsp
    integer :: id_sname = -1
    call tstatc0_begin('m_ES_WF_2D_psi ',id_sname)
    distrib = .true.
    if(present(bff)) distrib = .false.

    if(.not.allocated(bfft)) allocate(bfft(nfft))
    if(size(bfft) .ne. nfft) then
       deallocate(bfft)
       allocate(bfft(nfft))
    endif
    iesize = ib2-ib1+1
    nx = fft_box_size_WF(1,1)
    ny = fft_box_size_WF(2,1)
    nz = fft_box_size_WF(3,1)

    if(inverse_or_direct == INVERSE) then
       if(.not.allocated(zajbuf).and.nrank_g>1) allocate(zajbuf(kg1,kimg))
       if(nrank_g>1 .and. size(zajbuf) .ne. kg1*kimg)then
         deallocate(zajbuf)
         allocate(zajbuf(kg1,kimg))
       endif

       nl = xyz_fft_y(2,1) - xyz_fft_y(1,1) + 1
       nm = xyz_fft_y(2,2) - xyz_fft_y(1,2) + 1
       nn = xyz_fft_y(2,3) - xyz_fft_y(1,3) + 1
       bfft_l = 0.d0
       do ib=1,iesize
          if(nrank_g>1) then
            zajbuf = 0.d0
            do ir=1,kimg
               do ig=ista_g1k(ik),iend_g1k(ik)
                  zajbuf(ig,ir) = psi_l(ig-ista_g1k(ik)+1,ib+ib1-1,ik,ir)
               enddo
            enddo
            call mpi_allreduce(MPI_IN_PLACE,zajbuf,kg1*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
            call wf_in_rspace_2D(ik,zajbuf,bfft)
            if(distrib) then
               iadd = 0
               if(kimg==1)then
                  do gj = xyz_fft_y(1,2),xyz_fft_y(2,2)
                     do gl = xyz_fft_y(1,3),xyz_fft_y(2,3)
                        do gk = xyz_fft_y(1,1),xyz_fft_y(2,1)
                           iadd1 = (gk+(gl-1)*nx+(gj-1)*nx*nz)
                           iadd  = iadd + 1
                           bfft_l(iadd,ib) = bfft(iadd1)
                        end do
                     end do
                  enddo
               else
                  do gj = xyz_fft_z(1,3),xyz_fft_z(2,3)
                     do gl = xyz_fft_z(1,2),xyz_fft_z(2,2)
                        do gk = xyz_fft_z(1,1),xyz_fft_z(2,1)
                           iadd1 = (gk+(gl-1)*nx+(gj-1)*nx*ny)
                           iadd  = iadd + 1
                           bfft_l(2*iadd,ib)   = bfft(2*iadd1)
                           bfft_l(2*iadd-1,ib) = bfft(2*iadd1-1)
                        end do
                     end do
                  enddo
               endif
            else
               bff(1:nfft,ib) = bfft(1:nfft)
            endif
          else
            if(distrib) then
               call wf_in_rspace_2D(ik,psi_l(:,ib+ib1-1,ik,:),bfft_l(:,ib))
            else
               call wf_in_rspace_2D(ik,psi_l(:,ib+ib1-1,ik,:),bff(:,ib))
            endif
          endif
       enddo
    else
       do ib=1,iesize
          if(nrank_g>1)then
            iadd = 0
            bfft = 0.d0
            if(kimg==1)then
               do gj = xyz_fft_z(1,3),xyz_fft_z(2,3)
                  do gl = xyz_fft_z(1,2),xyz_fft_z(2,2)
                     do gk = xyz_fft_z(1,1),xyz_fft_z(2,1)
                        iadd1 = (gk+(gl-1)*nx+(gj-1)*nx*ny)
                        iadd  = iadd + 1
                        bfft(iadd1) = bfft_l(iadd,ib)
                     end do
                  end do
               enddo
            else
               do gj = xyz_fft_z(1,3),xyz_fft_z(2,3)
                  do gl = xyz_fft_z(1,2),xyz_fft_z(2,2)
                     do gk = xyz_fft_z(1,1),xyz_fft_z(2,1)
                        iadd1 = (gk+(gl-1)*nx+(gj-1)*nx*ny)
                        iadd  = iadd + 1
                        bfft(2*iadd1) = bfft_l(2*iadd,ib)
                        bfft(2*iadd1-1) = bfft_l(2*iadd-1,ib)
                     end do
                  end do
               enddo
            endif
            call mpi_allreduce(MPI_IN_PLACE,bfft,nfft,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
            call m_FFT_WF(ELECTRON,nfout,bfft,DIRECT,ON)
            if(distrib) then
               iadd = 0
               bfft_l = 0.d0
               if(kimg==1)then
                  do gl = xyz_fft_x(1,3),xyz_fft_x(2,3)
                     do gj = xyz_fft_x(1,2),xyz_fft_x(2,2)
                        do gk = xyz_fft_x(1,1),xyz_fft_x(2,1)
                           iadd1 = (gk+(gj-1)*nx+(gl-1)*nx*ny)
                           iadd  = iadd + 1
                           bfft_l(iadd,ib) = bfft(iadd1)
                        end do
                     end do
                  enddo
               else
                  do gl = xyz_fft_x(1,3),xyz_fft_x(2,3)
                     do gj = xyz_fft_x(1,2),xyz_fft_x(2,2)
                        do gk = xyz_fft_x(1,1),xyz_fft_x(2,1)
                           iadd1 = (gk+(gj-1)*nx+(gl-1)*nx*ny)
                           iadd  = iadd + 1
                           bfft_l(2*iadd,ib) = bfft(2*iadd1)
                           bfft_l(2*iadd-1,ib) = bfft(2*iadd1-1)
                        end do
                     end do
                  enddo
               endif
            else
               bff(1:nfft,ib) = bfft(1:nfft)
            endif
          else
            if(distrib)then
               call m_FFT_WF(ELECTRON,nfout,bfft_l(:,ib),DIRECT,ON)
            else
               call m_FFT_WF(ELECTRON,nfout,bff(:,ib),DIRECT,ON)
            endif
          endif
       enddo
    endif

    call tstatc0_end(id_sname)
    contains

    subroutine wf_in_rspace_2D(ik,psi_l,bfft)
      integer, intent(in):: ik
      real(kind=DP), dimension(kg1,kimg), intent(in) :: psi_l
      real(kind=DP), dimension(nfft), intent(out) :: bfft
      integer :: i,i1,ri, j, i2, ii
      integer :: id_sname2=-1
      call tstatc0_begin('wf_in_rspace_2D ',id_sname2)

      bfft = 0.d0
      if(k_symmetry(ik) == GAMMA) then
         if(kimg == 1) then
            i1 = igf(1)
            bfft(i1) = psi_l(1,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
               i = nbase(ii,ik)
               i1 = igf(i)
               bfft(i1) = psi_l(ii,1)
               j = nbase_gamma(ii,2)
               i2 = igf(j)
               bfft(i2) =   psi_l(ii,1)
            end do
         else if(kimg == 2) then
            i1 = 2*igf(1) - 1
            bfft(i1)   = psi_l(1,1)
            bfft(i1+1) = psi_l(1,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
               i = nbase(ii,ik)
               i1 = 2*igf(i)-1
               bfft(i1  ) = psi_l(ii,1)
               bfft(i1+1) = psi_l(ii,2)
               j = nbase_gamma(ii,2)
               i2 = 2*igf(j)-1
               bfft(i2  ) = psi_l(ii,1)
               bfft(i2+1) = -psi_l(ii,2)
            end do
         end if
      else
#ifdef NEC_TUNE_SMP
!CDIR NOLOOPCHG
#endif
         do ri = 1, kimg
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
            do i = 1, iba(ik)
               i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
               bfft(i1) = psi_l(i,ri)   ! MPI
            end do
         end do
      end if

      call m_FFT_WF(ELECTRON,nfout,bfft,inverse_or_direct,ON)

      call tstatc0_end(id_sname2)
  end subroutine wf_in_rspace_2D

  end subroutine m_ES_WF_2D_psi

#ifdef __TIMER_COMM__
  subroutine m_ES_WF_in_Rspace_3D0(ik,ib1,ib2,ibsize,lsize,bfft_l, kukanNo)
#else
  subroutine m_ES_WF_in_Rspace_3D0(ik,ib1,ib2,ibsize,lsize,bfft_l)
#endif
    integer, intent(in)                           :: ik, ib1,ib2, ibsize,lsize
#ifdef __TIMER_COMM__
    integer, intent(in)                           :: kukanNo
#endif
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
#ifdef __TIMER_COMM__
    call m_ES_WF_in_Rspace_3D1(ista_k,iend_k,ik,ib1,ib2,ibsize,lsize,zaj_l,bfft_l,kukanNo)
#else
    call m_ES_WF_in_Rspace_3D1(ista_k,iend_k,ik,ib1,ib2,ibsize,lsize,zaj_l,bfft_l)
#endif
  end subroutine m_ES_WF_in_Rspace_3D0

#ifdef __TIMER_COMM__
  subroutine m_ES_WF_in_Rspace_3D1(k1,k2,ik,ib1,ib2,ibsize,lsize,zaj_l,bfft_l, kukanNo)
#else
  subroutine m_ES_WF_in_Rspace_3D1(k1,k2,ik,ib1,ib2,ibsize,lsize,zaj_l,bfft_l)
#endif
    use m_Parallelization,     only : wf_fft_scnt, wf_fft_rcnt &
   &                                , wf_fft_recv &
   &                                , wf_fft_index, wf_fft_dist &
   &                                , wf_fft_maxrecv, wf_fft_maxsend
    integer, intent(in)                           :: k1,k2,ik, ib1,ib2, ibsize,lsize
#ifdef __TIMER_COMM__
    integer, intent(in)                           :: kukanNo
#endif
    real(kind=DP), intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: zaj_l
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
!   real(kind=DP), dimension(nfft,ibsize) :: bfft
! === DEBUG by tkato 2013/08/28 ================================================
!   real(kind=DP), dimension(1,1) :: bfft
! ==============================================================================
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
!!!ifdef FFT_USE_SSL2
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
    integer,save  :: id_sname = -1
#ifdef SAVE_FFT_TIMES
    logical :: fft_stored
    integer,save :: id_sname2 = -1
#endif
    integer :: ib,ig, ir, gk,gj,gl,gi,nnx,iadd0,iadd1,nis,nie
    real(kind=DP), allocatable, dimension(:,:) :: zajbuf
    real(kind=DP), allocatable, dimension(:) :: bfft
    real(kind=DP) :: t

! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       fft_stored = .true.
       do jj = ib1, ib2
          if(status_saved_phifftr(jj,ik) /= STORED_AND_NEW) fft_stored = .false.
       end do
       if(ipri >= ipri_save_fft) then
          if(ib1 == 1 .and. ik == 1) &
               & write(nfout,'(" fft_stored (",i4,":",i4,",",i3,") = ",L)') ib1,ib2,ik,fft_stored
       end if
    else
       fft_stored = .false.
    end if
!!$!!x!!$    if(status_saved_phifftr(map_z(ib),ik) == STORED_AND_NEW) then
    if(fft_stored) then
       call tstatc0_begin('m_ES_WF_in_Rspace_3D(2) ',id_sname2)
       if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_fftr(stored) --> bfft")')
       iesize = ib2 - ib1 + 1
       bfft_l(:,1:iesize) = Phifftr_l(:,ib1:ib2,ik)
       call tstatc0_end(id_sname2)
    else
#endif
                                                 __TIMER_SUB_START(304)
    call tstatc0_begin('m_ES_WF_in_Rspace_3D(1) ',id_sname)
! serial or paralle FFT
    if(sw_serial_fft == ON)then
!serial FFT
       call m_ES_WF_2D(ik,bfft_l,ib2,ib1,ibsize,lsize,INVERSE)
    else
#ifdef __FAPP__
    call fapp_start('wf_in_rspace_pre',1,1)
#endif
!parallel FFT

    if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_l --(FFT)--> bfft")')

    bfft_l = 0.0d0

    iesize = ib2 - ib1 + 1

!   allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!   allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!   if ((iesize /= savesize) .or. (wf_fft_maxsend(ik) /= savesend)) then
       if (allocated(sendbuf)) deallocate(sendbuf)
       allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!      savesend = wf_fft_maxsend(ik)
!   end if
!   if ((iesize /= savesize) .or. (wf_fft_maxrecv(ik) /= savesend)) then
       if (allocated(recvbuf)) deallocate(recvbuf)
       allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!      saverecv = wf_fft_maxrecv(ik)
!   end if
    savesize = iesize
!    sendbuf = 0.0d0
!    recvbuf = 0.0d0
                                                 __TIMER_DO_START(313)
    if(k_symmetry(ik) == GAMMA) then
       if (kimg == 1) then
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                sendbuf(iesize*(wf_fft_index(ii*2-1,ik)-1)+1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iesize*(wf_fft_index(ii*2  ,ik)-1)+1,wf_fft_dist(ii*2  ,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii*2-1,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
                iadd = iesize*2*(wf_fft_index(ii*2,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2  ,ik)) =  zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2  ,ik)) = -zaj_l(ii,ib1+jj-1,ik,2)
             enddo
          enddo
       endif
    else
       if (kimg == 1) then
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*(wf_fft_index(ii,ik)-1)+jj
                sendbuf(iadd,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd  ,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
             enddo
          enddo
       endif
    endif
                                                 __TIMER_DO_STOP(313)

#ifndef USE_ALLTOALLV
    icnt_recv = 0
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,314)
    do i = 0, nrank_g - 1
       if (wf_fft_rcnt(i,ik) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,166,ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if (wf_fft_scnt(i,ik) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,167,ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,168,ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(314)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,169,ierr)
     endif
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,329)
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=wf_fft_maxsend(ik)*kimg*iesize*i
          rdsp(i)=wf_fft_maxrecv(ik)*kimg*iesize*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt(:,ik)*kimg*iesize, sdsp,&
      &   mpi_double_precision, recvbuf, wf_fft_rcnt(:,ik)*kimg*iesize, rdsp,&
      &   mpi_double_precision, mpi_ke_world, ierr )
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 170, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                 __TIMER_COMM_STOP(329)
#endif
                                                 __TIMER_DO_START(315)
#ifdef FFT_3D_DIVISION
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
!               bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = 0.0d0
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif
#endif
                                                 __TIMER_DO_STOP(315)
!!  deallocate(sendbuf)
!!  deallocate(recvbuf)

#ifdef __FAPP__
    call fapp_stop('wf_in_rspace_pre',1,1)
#endif

#ifdef __TIMER_COMM__
    if(kukanNo==3) call timer_sta(327)
    if(kukanNo==6) call timer_sta(636)
    if(kukanNo==7) call timer_sta(691)
    if(kukanNo==8) call timer_sta(961)
#endif
#ifdef __FAPP__
    call fapp_start('wf_in_rspace_fft',1,1)
#endif
#ifdef FFT_3D_DIVISION
    call m_FFT_Inverse_3DIV_3D(nfout, bfft_l, lsize, iesize)
#else
    if (sw_fft_xzy > 0) then
       call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
    else
       call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
    end if
#endif
#ifdef __TIMER_COMM__
    if(kukanNo==3) call timer_end(327)
    if(kukanNo==6) call timer_end(636)
    if(kukanNo==7) call timer_end(691)
    if(kukanNo==8) call timer_end(961)
#endif
#ifdef __FAPP__
    call fapp_stop('wf_in_rspace_fft',1,1)
#endif

    if(ipri >= 2) then
       if(ik <=  2 .and. ib1 <= 1) then
          write(6,'(" ! bfft R-space ik = ",i3," ib1 = ",i3," <<m_ES_WF_in_Rspace>>")') ik, ib1
! === DEBUG by tkato 2013/08/28 ================================================
!         write(6,'(8f8.4)') (bfft(i,1),i=1,120)
          write(6,'(8f8.4)') (bfft_l(i,1),i=1,120)
! ==============================================================================
       end if
    end if

                                                 __TIMER_SUB_STOP(304)

! serial or parallel FFT
    endif

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       iesize = ib2 - ib1 + 1
       Phifftr_l(:,ib1:ib2,ik) = bfft_l(:,1:iesize)
       status_saved_phifftr(ib1:ib2,ik) = STORED_AND_NEW
    end if
#endif
    call tstatc0_end(id_sname)
#ifdef SAVE_FFT_TIMES
    end if
#endif

  end subroutine m_ES_WF_in_Rspace_3D1

  ! use out-of-place FFT
  subroutine m_ES_WF_in_Rspace_3D2(ik,ib1,ib2,ibsize,lsize,bfft_l_in,bfft_l)
    use m_Parallelization,     only : wf_fft_scnt, wf_fft_rcnt &
   &                                , wf_fft_recv &
   &                                , wf_fft_index, wf_fft_dist &
   &                                , wf_fft_maxrecv, wf_fft_maxsend
    integer, intent(in)                           :: ik, ib1,ib2, ibsize,lsize
#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l_in
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l_in
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif

!   real(kind=DP), dimension(nfft,ibsize) :: bfft
! === DEBUG by tkato 2013/08/28 ================================================
!   real(kind=DP), dimension(1,1) :: bfft
! ==============================================================================
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
!!!ifdef FFT_USE_SSL2
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
    integer,save  :: id_sname = -1
#ifdef SAVE_FFT_TIMES
    logical :: fft_stored
    integer,save :: id_sname2 = -1
#endif
    integer :: ib,ig, ir, gk,gj,gl,gi,nnx,iadd0,iadd1,nis,nie
    real(kind=DP), allocatable, dimension(:,:) :: zajbuf
    real(kind=DP), allocatable, dimension(:) :: bfft
    real(kind=DP) :: t

! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       fft_stored = .true.
       do jj = ib1, ib2
          if(status_saved_phifftr(jj,ik) /= STORED_AND_NEW) fft_stored = .false.
       end do
       if(ipri >= ipri_save_fft) then
          if(ib1 == 1 .and. ik == 1) &
               & write(nfout,'(" fft_stored (",i4,":",i4,",",i3,") = ",L)') ib1,ib2,ik,fft_stored
       end if
    else
       fft_stored = .false.
    end if
!!$!!x!!$    if(status_saved_phifftr(map_z(ib),ik) == STORED_AND_NEW) then
    if(fft_stored) then
       call tstatc0_begin('m_ES_WF_in_Rspace_3D(2) ',id_sname2)
       if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_fftr(stored) --> bfft")')
       iesize = ib2 - ib1 + 1
       bfft_l(:,1:iesize) = Phifftr_l(:,ib1:ib2,ik)
       call tstatc0_end(id_sname2)
    else
#endif
    call tstatc0_begin('m_ES_WF_in_Rspace_3D(1) ',id_sname)
! serial or paralle FFT
    if(sw_serial_fft == ON)then
!serial FFT
       call m_ES_WF_2D(ik,bfft_l,ib2,ib1,ibsize,lsize,INVERSE)
    else
!parallel FFT

    if(ipri>=2 .and. ik==1 .and. ib1==1) write(nfout,'(" !### zaj_l --(FFT)--> bfft")')

    iesize = ib2 - ib1 + 1

!   allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!   allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!   if ((iesize /= savesize) .or. (wf_fft_maxsend(ik) /= savesend)) then
       if (allocated(sendbuf)) deallocate(sendbuf)
       allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!      savesend = wf_fft_maxsend(ik)
!   end if
!   if ((iesize /= savesize) .or. (wf_fft_maxrecv(ik) /= savesend)) then
       if (allocated(recvbuf)) deallocate(recvbuf)
       allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!      saverecv = wf_fft_maxrecv(ik)
!   end if
    savesize = iesize
!    sendbuf = 0.0d0
!    recvbuf = 0.0d0
    if(k_symmetry(ik) == GAMMA) then
       if (kimg == 1) then
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                sendbuf(iesize*(wf_fft_index(ii*2-1,ik)-1)+1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iesize*(wf_fft_index(ii*2  ,ik)-1)+1,wf_fft_dist(ii*2  ,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii*2-1,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
                iadd = iesize*2*(wf_fft_index(ii*2,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2  ,ik)) =  zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2  ,ik)) = -zaj_l(ii,ib1+jj-1,ik,2)
             enddo
          enddo
       endif
    else
       if (kimg == 1) then
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*(wf_fft_index(ii,ik)-1)+jj
                sendbuf(iadd,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                iadd = iesize*2*(wf_fft_index(ii,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,1)
                sendbuf(iadd  ,wf_fft_dist(ii,ik)) = zaj_l(ii,ib1+jj-1,ik,2)
             enddo
          enddo
       endif
    endif

#ifndef USE_ALLTOALLV
    icnt_recv = 0
    do i = 0, nrank_g - 1
       if (wf_fft_rcnt(i,ik) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,166,ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if (wf_fft_scnt(i,ik) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,167,ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,168,ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,169,ierr)
     endif
#else
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=wf_fft_maxsend(ik)*kimg*iesize*i
          rdsp(i)=wf_fft_maxrecv(ik)*kimg*iesize*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt(:,ik)*kimg*iesize, sdsp,&
      &   mpi_double_precision, recvbuf, wf_fft_rcnt(:,ik)*kimg*iesize, rdsp,&
      &   mpi_double_precision, mpi_ke_world, ierr )
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 170, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
#endif

#ifdef FFT_3D_DIVISION
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
!               bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l_in(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l_in(wf_fft_recv(j,ik,i)*2  ,k) = 0.0d0
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l_in(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l_in(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l_in(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l_in((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l_in((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l_in(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l_in(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l_in(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif
#endif

!!!#ifdef FFT_3D_DIVISION
!!!    call m_FFT_Inverse_3DIV_3D(nfout, bfft_l, lsize, iesize)
!!!#else
!!!    if (sw_fft_xzy > 0) then
!!!       call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
!!!    else
       call m_FFT_Inverse_XYZ_3D_oo_place(nfout, bfft_l_in, bfft_l, lsize, iesize)
!!!    end if

    if(ipri >= 2) then
       if(ik <=  2 .and. ib1 <= 1) then
          write(6,'(" ! bfft R-space ik = ",i3," ib1 = ",i3," <<m_ES_WF_in_Rspace>>")') ik, ib1
! === DEBUG by tkato 2013/08/28 ================================================
!         write(6,'(8f8.4)') (bfft(i,1),i=1,120)
          write(6,'(8f8.4)') (bfft_l(i,1),i=1,120)
! ==============================================================================
       end if
    end if

! serial or paralle FFT
    endif

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       iesize = ib2 - ib1 + 1
       Phifftr_l(:,ib1:ib2,ik) = bfft_l(:,1:iesize)
       status_saved_phifftr(ib1:ib2,ik) = STORED_AND_NEW
    end if
#endif
    call tstatc0_end(id_sname)
#ifdef SAVE_FFT_TIMES
    end if
#endif
  end subroutine m_ES_WF_in_Rspace_3D2

  subroutine m_ES_WF_in_Rspace_kt_3D(k1,k2,ik,ibsize,lsize,zaj_l,bfft_l)
    use m_Parallelization,     only : wf_fft_scnt, wf_fft_rcnt &
   &                                , wf_fft_recv &
   &                                , wf_fft_index, wf_fft_dist &
   &                                , wf_fft_maxrecv, wf_fft_maxsend
    integer, intent(in)                           :: k1,k2,ik, ibsize,lsize
    real(kind=DP), intent(in),dimension(maxval(np_g1k),1,k1:k2,kimg) :: zaj_l

#ifdef FFT_3D_DIVISION
    real(kind=DP), intent(inout), dimension(lsize*2,ibsize) :: bfft_l
#else
    real(kind=DP), intent(inout), dimension(lsize*kimg,ibsize) :: bfft_l
#endif
!   real(kind=DP), dimension(nfft,ibsize) :: bfft
! === DEBUG by tkato 2013/08/28 ================================================
!   real(kind=DP), dimension(1,1) :: bfft
! ==============================================================================
    integer :: i, j, k, ii, jj, iesize, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
    integer, save :: savesize = 0, savesend=0, saverecv=0
!!!ifdef FFT_USE_SSL2
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
    integer,save  :: id_sname = -1

    integer :: ib,ig, ir, gk,gj,gl,gi,nnx,iadd0,iadd1,nis,nie
    real(kind=DP), allocatable, dimension(:,:) :: zajbuf
    real(kind=DP), allocatable, dimension(:) :: bfft
    real(kind=DP) :: t

! === DEBUG by tkato 2012/06/05 ================================================
#ifdef USE_ALLTOALLV
    integer, allocatable, dimension(:) :: sdsp, rdsp
#endif
! ==============================================================================

    integer :: skip_gvec_unfolding, i1

    call tstatc0_begin('m_ES_WF_in_Rspace_kt_3D ',id_sname)

!parallel FFT

    iesize = 1

!   allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!   allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!   if ((iesize /= savesize) .or. (wf_fft_maxsend(ik) /= savesend)) then
       if (allocated(sendbuf)) deallocate(sendbuf)
       allocate(sendbuf(wf_fft_maxsend(ik)*kimg*iesize,0:nrank_g-1))
!      savesend = wf_fft_maxsend(ik)
!   end if
!   if ((iesize /= savesize) .or. (wf_fft_maxrecv(ik) /= savesend)) then
       if (allocated(recvbuf)) deallocate(recvbuf)
       allocate(recvbuf(wf_fft_maxrecv(ik)*kimg*iesize,0:nrank_g-1))
!      saverecv = wf_fft_maxrecv(ik)
!   end if
    savesize = iesize
    sendbuf = 0.0d0
    recvbuf = 0.0d0
                                                 __TIMER_DO_START(313)

    if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       skip_gvec_unfolding = 1
    else
       skip_gvec_unfolding = 0
    endif

    if(k_symmetry(ik) == GAMMA) then
       if (kimg == 1) then
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                if ( skip_gvec_unfolding == 1 ) then
                   i1 = ii +ista_g1k(ik) -1
                   if ( sw_force_kpt_inside_bz == ON ) then
                      if ( GVec_on_refcell(i1,ik) == 0 ) cycle
                   else
                      if ( GVec_on_refcell(i1,1) == 0 ) cycle
                   endif
                endif
                sendbuf(iesize*(wf_fft_index(ii*2-1,ik)-1)+1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,jj,ik,1)
                sendbuf(iesize*(wf_fft_index(ii*2  ,ik)-1)+1,wf_fft_dist(ii*2  ,ik)) = zaj_l(ii,jj,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                if ( skip_gvec_unfolding == 1 ) then
                   i1 = ii +ista_g1k(ik) -1
                   if ( sw_force_kpt_inside_bz == ON ) then
                      if ( GVec_on_refcell(i1,ik) == 0 ) cycle
                   else
                      if ( GVec_on_refcell(i1,1) == 0 ) cycle
                   endif
                endif
                iadd = iesize*2*(wf_fft_index(ii*2-1,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,jj,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2-1,ik)) = zaj_l(ii,jj,ik,2)
                iadd = iesize*2*(wf_fft_index(ii*2,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii*2  ,ik)) =  zaj_l(ii,jj,ik,1)
                sendbuf(iadd,  wf_fft_dist(ii*2  ,ik)) = -zaj_l(ii,jj,ik,2)
             enddo
          enddo
       endif
    else
       if (kimg == 1) then
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                if ( skip_gvec_unfolding == 1 ) then
                   i1 = ii +ista_g1k(ik) -1
                   if ( sw_force_kpt_inside_bz == ON ) then
                      if ( GVec_on_refcell(i1,ik) == 0 ) cycle
                   else
                      if ( GVec_on_refcell(i1,1) == 0 ) cycle
                   endif
                endif
                iadd = iesize*(wf_fft_index(ii,ik)-1)+jj
                sendbuf(iadd,wf_fft_dist(ii,ik)) = zaj_l(ii,jj,ik,1)
             enddo
          enddo
       else
!OCL NORECURRENCE
          do jj = 1, iesize
             do ii = 1, np_g1k(ik)
                if ( skip_gvec_unfolding == 1 ) then
                   i1 = ii +ista_g1k(ik) -1
                   if ( sw_force_kpt_inside_bz == ON ) then
                      if ( GVec_on_refcell(i1,ik) == 0 ) cycle
                   else
                      if ( GVec_on_refcell(i1,1) == 0 ) cycle
                   endif
                endif
                iadd = iesize*2*(wf_fft_index(ii,ik)-1)+jj*2
                sendbuf(iadd-1,wf_fft_dist(ii,ik)) = zaj_l(ii,jj,ik,1)
                sendbuf(iadd  ,wf_fft_dist(ii,ik)) = zaj_l(ii,jj,ik,2)
             enddo
          enddo
       endif
    endif
                                                 __TIMER_DO_STOP(313)

#ifndef USE_ALLTOALLV
    icnt_recv = 0
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,314)
    do i = 0, nrank_g - 1
       if (wf_fft_rcnt(i,ik) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,166,ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if (wf_fft_scnt(i,ik) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt(i,ik)*kimg*iesize, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,167,ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,168,ierr)
     endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(314)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,169,ierr)
     endif
#else
                                                 __TIMER_COMM_START_w_BARRIER(mpi_ke_world,329)
! === DEBUG by tkato 2012/06/05 ================================================
!      integer, allocatable, dimension(:) :: sdsp, rdsp
! ==============================================================================
       allocate(sdsp(0:nrank_g-1), stat=ierr)
       allocate(rdsp(0:nrank_g-1), stat=ierr)
       do i = 0, nrank_g - 1
          sdsp(i)=wf_fft_maxsend(ik)*kimg*iesize*i
          rdsp(i)=wf_fft_maxrecv(ik)*kimg*iesize*i
       enddo
       call MPI_ALLTOALLV(      sendbuf, wf_fft_scnt(:,ik)*kimg*iesize, sdsp,&
      &   mpi_double_precision, recvbuf, wf_fft_rcnt(:,ik)*kimg*iesize, rdsp,&
      &   mpi_double_precision, mpi_ke_world, ierr )
       if (ierr /= 0) then
          write(nfout,*)' m_ES_Vlocal_in_Rspace_3D : mpi_alltoallv error'
          call flush(nfout)
          call mpi_abort(mpi_comm_world, 170, ierr)
       endif
       deallocate(sdsp)
       deallocate(rdsp)
                                                 __TIMER_COMM_STOP(329)
#endif

    bfft_l = 0.0d0
                                                 __TIMER_DO_START(315)
#ifdef FFT_3D_DIVISION
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
!               bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*(j-1)+k,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = 0.0d0
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
#ifdef FFT_USE_SSL2
    nx = xyz_fft_x(2,1) - xyz_fft_x(1,1) + 1
    ny = xyz_fft_x(2,2) - xyz_fft_x(1,2) + 1
    nz = xyz_fft_x(2,3) - xyz_fft_x(1,3) + 1
    nxp = nx
    if (kimg == 1) then
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l(ix+(iy-1)*nx+(iz-1)*nx*ny,k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
!OCL NOFLTLD
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             iadd = wf_fft_recv(j,ik,i)
             iz = (iadd-1)/(nx*ny)+1
             nn = mod(iadd,nx*ny)
             if (nn==0) then
                iy = ny
             else
                iy = (nn-1)/nx+1
             end if
             ix = mod(nn,nx)
             if(ix==0) ix = nx
             do k = 1, iesize
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l((ix+(iy-1)*nxp+(iz-1)*nxp*ny)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#else
    if (kimg == 1) then
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i),k) = recvbuf(iesize*(j-1)+k,i)
             enddo
          enddo
       enddo
    else
!OCL NORECURRENCE
       do i = 0, nrank_g - 1
          do j = 1, wf_fft_rcnt(i,ik)
             do k = 1, iesize
                bfft_l(wf_fft_recv(j,ik,i)*2-1,k) = recvbuf(iesize*2*(j-1)+k*2-1,i)
                bfft_l(wf_fft_recv(j,ik,i)*2  ,k) = recvbuf(iesize*2*(j-1)+k*2  ,i)
             enddo
          enddo
       enddo
    endif
#endif
#endif
                                                 __TIMER_DO_STOP(315)
!!  deallocate(sendbuf)
!!  deallocate(recvbuf)

#ifdef FFT_3D_DIVISION
    call m_FFT_Inverse_3DIV_3D(nfout, bfft_l, lsize, iesize)
#else
    if (sw_fft_xzy > 0) then
       call m_FFT_Inverse_3D(nfout, bfft_l, lsize, iesize)
    else
       call m_FFT_Inverse_XYZ_3D(nfout, bfft_l, lsize, iesize)
    end if
#endif
                                                 __TIMER_SUB_STOP(304)

    call tstatc0_end(id_sname)

  end subroutine m_ES_WF_in_Rspace_kt_3D

#ifdef MPI_FFTW
  subroutine m_ES_fftbox_map(ik)
    use m_Parallelization, only : wf_fft_rcnt_mfftw, wf_fft_recv_mfftw
    integer, intent(in) :: ik
    integer :: i,j,icount,i1,mm,my,mx,mz,lx,ly,lz
   lx = fft_box_size_WF(1,0)
   ly = fft_box_size_WF(2,0)
   lz = fft_box_size_WF(3,0)
   icount = 0
   do i = 0, nrank_g - 1
     do j = 1, wf_fft_rcnt_mfftw(i,ik)
        icount = icount+1
     enddo
   enddo
   maxrecv = icount
   if(allocated(mmx)) deallocate(mmx)
   if(allocated(mmy)) deallocate(mmy)
   if(allocated(mmz)) deallocate(mmz)
   allocate(mmx(maxrecv))
   allocate(mmy(maxrecv))
   allocate(mmz(maxrecv))
   icount = 0
   if(kimg==2) then
     do i = 0, nrank_g - 1
       do j = 1, wf_fft_rcnt_mfftw(i,ik)
          i1 = wf_fft_recv_mfftw(j,ik,i)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,(lx*ly))
          if (mm==0) mm=lx*ly
          my = (mm-1)/lx+1
          mx = mod(mm,lx)
          if (mx==0) mx = lx
          icount = icount+1
          mmx(icount) = mx
          mmy(icount) = my
          mmz(icount) = mz
       enddo
     enddo
   else
     do i = 0, nrank_g - 1
       do j = 1, wf_fft_rcnt_mfftw(i,ik)
          i1 = wf_fft_recv_mfftw(j,ik,i)
          mz = (i1-1)/(lx*ly)+1
          mm = mod(i1,lx*ly)
          if (mm==0) then
             my = ly
          else
             my = (mm-1)/lx+1
          end if
          mx = mod(mm,lx)
          if(mx==0) mx = lx
          icount = icount+1
          mmx(icount) = mx
          mmy(icount) = my
          mmz(icount) = mz
       enddo
     enddo
   endif
  end subroutine m_ES_fftbox_map

  subroutine m_ES_WF_in_Rspace_mpifftw(k1,k2,ik,ib1,zaj_l)
    use m_FFT,                 only : m_FFT_Inverse_MPI_FFTW, bfft_mpifftw, afft_mpifftw, bfft_mpifftw_kimg1, afft_mpifftw_kimg1
    use m_Parallelization,     only : wf_fft_scnt_mfftw, wf_fft_rcnt_mfftw &
   &                                , wf_fft_recv_mfftw &
   &                                , wf_fft_index_mfftw, wf_fft_dist_mfftw &
   &                                , wf_fft_maxrecv_mfftw, wf_fft_maxsend_mfftw
    integer, intent(in)                           :: k1,k2,ik,ib1
    real(kind=DP), intent(in),dimension(maxval(np_g1k),np_e,k1:k2,kimg) :: zaj_l

!   real(kind=DP), dimension(nfft,ibsize) :: bfft
! === DEBUG by tkato 2013/08/28 ================================================
!   real(kind=DP), dimension(1,1) :: bfft
! ==============================================================================
    integer :: i, j, k, ii, iadd, ierr
    integer, dimension(0:nrank_g-1) ::req_r,req_s
    integer, dimension(MPI_STATUS_SIZE,0:nrank_g-1)::sta_r, sta_s
    real(kind=DP), allocatable, dimension(:,:),save :: sendbuf, recvbuf
    integer :: icnt_send, icnt_recv, lrank
    integer, parameter :: itag = 21
!!!ifdef FFT_USE_SSL2
    integer :: nx, ny, nz, nxp, nn, ix, iy, iz
!!!endif
    integer,save  :: id_sname = -1, id_sname1=-1, id_sname2=-1

    integer :: ib,ig, ir, gk,gj,gl,gi,nnx,iadd0,iadd1,nis,nie
    real(kind=DP), allocatable, dimension(:,:) :: zajbuf
    real(kind=DP), allocatable, dimension(:) :: bfft
    real(kind=DP) :: t
    complex(kind=DP), parameter :: ZERO = (0.d0,0.d0)
    integer :: i1,lx,ly,lz,mm
    integer :: icount
    integer(C_INTPTR_T) :: mz,mx,my

    call tstatc0_begin('m_ES_WF_in_Rspace_mpifft ',id_sname)

!parallel FFT

    call tstatc0_begin('m_ES_WF_in_Rspace_mpifft_1 ',id_sname1)
    if (allocated(sendbuf)) deallocate(sendbuf)
    allocate(sendbuf(wf_fft_maxsend_mfftw(ik)*kimg,0:nrank_g-1))
    if (allocated(recvbuf)) deallocate(recvbuf)
    allocate(recvbuf(wf_fft_maxrecv_mfftw(ik)*kimg,0:nrank_g-1))

    if(k_symmetry(ik) == GAMMA) then
       if (kimg == 1) then
             do ii = 1, np_g1k(ik)
                sendbuf(wf_fft_index_mfftw(ii*2-1,ik),wf_fft_dist_mfftw(ii*2-1,ik)) = zaj_l(ii,ib1,ik,1)
                sendbuf(wf_fft_index_mfftw(ii*2  ,ik),wf_fft_dist_mfftw(ii*2  ,ik)) = zaj_l(ii,ib1,ik,1)
             enddo
       else
!OCL NORECURRENCE
             do ii = 1, np_g1k(ik)
                iadd = 2*(wf_fft_index_mfftw(ii*2-1,ik))
                sendbuf(iadd-1,wf_fft_dist_mfftw(ii*2-1,ik)) = zaj_l(ii,ib1,ik,1)
                sendbuf(iadd,  wf_fft_dist_mfftw(ii*2-1,ik)) = zaj_l(ii,ib1,ik,2)
                iadd = 2*(wf_fft_index_mfftw(ii*2,ik))
                sendbuf(iadd-1,wf_fft_dist_mfftw(ii*2  ,ik)) =  zaj_l(ii,ib1,ik,1)
                sendbuf(iadd,  wf_fft_dist_mfftw(ii*2  ,ik)) = -zaj_l(ii,ib1,ik,2)
             enddo
       endif
    else
       if (kimg == 1) then
!OCL NORECURRENCE
             do ii = 1, np_g1k(ik)
                iadd = wf_fft_index_mfftw(ii,ik)
                sendbuf(iadd,wf_fft_dist_mfftw(ii,ik)) = zaj_l(ii,ib1,ik,1)
             enddo
       else
!OCL NORECURRENCE
             do ii = 1, np_g1k(ik)
                iadd = 2*(wf_fft_index_mfftw(ii,ik))
                sendbuf(iadd-1,wf_fft_dist_mfftw(ii,ik)) = zaj_l(ii,ib1,ik,1)
                sendbuf(iadd  ,wf_fft_dist_mfftw(ii,ik)) = zaj_l(ii,ib1,ik,2)
             enddo
       endif
    endif

    icnt_recv = 0
    do i = 0, nrank_g - 1
       if (wf_fft_rcnt_mfftw(i,ik) /= 0) then
          call mpi_irecv(recvbuf(1,i), wf_fft_rcnt_mfftw(i,ik)*kimg, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_r(icnt_recv), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_irecv error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,166,ierr)
           endif
          icnt_recv = icnt_recv + 1
       endif
    enddo
    icnt_send = 0
    do i = 0, nrank_g - 1
       if (wf_fft_scnt_mfftw(i,ik) /= 0) then
          call mpi_isend(sendbuf(1,i), wf_fft_scnt_mfftw(i,ik)*kimg, mpi_double_precision, &
      &                  i, itag, mpi_ke_world, req_s(icnt_send), ierr)
           if (ierr /= 0) then
              write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_isend error'
              call flush(nfout)
              call mpi_abort(mpi_comm_world,167,ierr)
           endif
          icnt_send = icnt_send + 1
       endif
    enddo
    call mpi_waitall(icnt_recv, req_r, sta_r, ierr)
    if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,168,ierr)
    endif
    call mpi_waitall(icnt_send, req_s, sta_s, ierr)
                                                 __TIMER_COMM_STOP(314)
     if (ierr /= 0) then
        write(nfout,*)' m_ES_WF_in_Rspace_3D :  mpi_waitall error'
        call flush(nfout)
        call mpi_abort(mpi_comm_world,169,ierr)
     endif
     call tstatc0_end(id_sname1)

     call tstatc0_begin('m_ES_WF_in_Rspace_mpifftw_2',id_sname2)

     lx = fft_box_size_WF(1,0)
     ly = fft_box_size_WF(2,0)
     lz = fft_box_size_WF(3,0)
     icount = 0
     if(kimg==2) then

       bfft_mpifftw = 0.d0

       do i = 0, nrank_g - 1
         do j = 1, wf_fft_rcnt_mfftw(i,ik)
!           i1 = wf_fft_recv_mfftw(j,ik,i)
!           mz = (i1-1)/(lx*ly)+1
!           mm = mod(i1,(lx*ly))
!           if (mm==0) mm=lx*ly
!           my = (mm-1)/lx+1
!           mx = mod(mm,lx)
!           if (mx==0) mx = lx
            icount = icount+1
            mx = mmx(icount)
            my = mmy(icount)
            mz = mmz(icount)
            bfft_mpifftw(mx,my,mz) = dcmplx(recvbuf(2*j-1,i),recvbuf(2*j,i))
         enddo
       enddo
     else
       bfft_mpifftw_kimg1 = 0.d0
       do i = 0, nrank_g - 1
         do j = 1, wf_fft_rcnt_mfftw(i,ik)
!           i1 = wf_fft_recv_mfftw(j,ik,i)
!           mz = (i1-1)/(lx*ly)+1
!           mm = mod(i1,lx*ly)
!           if (mm==0) then
!              my = ly
!           else
!              my = (mm-1)/lx+1
!           end if
!           mx = mod(mm,lx)
!           if(mx==0) mx = lx
            icount = icount+1
            mx = mmx(icount)
            my = mmy(icount)
            mz = mmz(icount)
            bfft_mpifftw_kimg1(mx,my,mz) = recvbuf(j,i)
         enddo
       enddo
     endif
     call tstatc0_end(id_sname2)

     call m_FFT_Inverse_MPI_FFTW(nfout)

     call tstatc0_end(id_sname)

  end subroutine m_ES_WF_in_Rspace_mpifftw
#endif

  subroutine m_ES_decide_precon_factor_3D_new(precon,ik,ib1,ib2,ibesize,ekin,p)
    integer, intent(in)                         :: precon,ik,ib1,ib2,ibesize
    real(kind=DP), intent(in),  dimension(np_g1k(ik)) :: ekin
    real(kind=DP), intent(out), dimension(mp_g1k(ik),ibesize) :: p

    integer       :: i, ib
    real(kind=DP) :: x, x1, x2
    real(kind=DP), dimension(ibesize) :: ektot, d_ektot
                                                 __TIMER_SUB_START(302)
    if(precon == ON) then
       call kinetic_energy_3D_new(ik,ib1,ib2,ibesize,ekin,ektot)   ! -(m_E.S.)
       d_ektot = 1.d0/ektot
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             x = ekin(i-ista_g1k(ik)+1)*d_ektot(ib-ib1+1)
             x1 = 27 + ( 18 + (12 + 8*x) *x) *x
             x2 = 16*(x*x)*(x*x)
             p(i-ista_g1k(ik)+1,ib-ib1+1)  = x1/(x1 + x2 )
          end do
       end do
    else
       p = 1.d0
    end if
                                                 __TIMER_SUB_STOP(302)
  end subroutine m_ES_decide_precon_factor_3D_new

  subroutine kinetic_energy_3D_new(ik,ib1,ib2,ibesize,dekin,ektot)
    integer, intent(in) :: ik, ib1, ib2, ibesize
    real(kind=DP), intent(in), dimension(np_g1k(ik)) :: dekin
    real(kind=DP), intent(out), dimension(ibesize) :: ektot
    real(kind=DP)             , dimension(ibesize) :: ektot_mpi
    integer  :: i, ib
                                                 __TIMER_SUB_START(303)

    ektot(:) = 0.d0
    if(kimg == 1) then
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             ektot(ib-ib1+1) = ektot(ib-ib1+1) + dekin(i-ista_g1k(ik)+1)* &
            &                                    zaj_l(i-ista_g1k(ik)+1,ib,ik,1)**2
          end do
       end do
    else
       do ib = ib1, ib2
          do i = ista_g1k(ik), iend_g1k(ik)
             ektot(ib-ib1+1) = ektot(ib-ib1+1) + dekin(i-ista_g1k(ik)+1)* &
            &                                    ( zaj_l(i-ista_g1k(ik)+1,ib,ik,1)**2 &
            &                                    + zaj_l(i-ista_g1k(ik)+1,ib,ik,2)**2 )
          end do
       end do
    end if

    call mpi_allreduce(ektot,ektot_mpi,ibesize,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    ektot = ektot_mpi

    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
                                                 __TIMER_SUB_STOP(303)
  end subroutine kinetic_energy_3D_new


  subroutine m_ES_add_it_to_vnlph(ik,ib,v)
    integer, intent(in) :: ik,ib
    real(kind=DP), intent(in) :: v(maxval(np_g1k),kimg)

    integer :: ig

    if(kimg==1) then
       do ig=1,np_g1k(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + v(ig,1)
       end do
    else
       do ig=1,np_g1k(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + v(ig,1)
          vnlph_l(ig,ib,2) = vnlph_l(ig,ib,2) + v(ig,2)
       end do
    end if
  end subroutine m_ES_add_it_to_vnlph

  subroutine m_ES_keep_retrieve_fsr( keep )
    logical, intent(in) :: keep
    real(kind=DP),allocatable,dimension(:,:,:),save :: fsr_tmp
    real(kind=DP),allocatable,dimension(:,:,:),save :: fsi_tmp

    if (keep) then
       if ( .not.allocated(fsr_tmp) ) allocate( fsr_tmp(np_e,np_fs,ista_k:iend_k))
       fsr_tmp = fsr_l
       if ( .not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2) ) then
          if ( .not.allocated(fsi_tmp) ) allocate(fsi_tmp(np_e,np_fs,ista_k:iend_k));
          fsi_tmp = fsi_l
       endif
    else
       fsr_l = fsr_tmp
       deallocate( fsr_tmp )
       if ( allocated(fsi_tmp) ) then
          fsi_l = fsi_tmp;   deallocate( fsi_tmp )
       endif
    endif
  end subroutine m_ES_keep_retrieve_fsr

! === For epsmain by tkato 2013/11/14 ==========================================
  subroutine m_ES_epsmain_reallocate()
     real(kind=DP), allocatable, dimension(:,:,:,:) :: zaj_t
     integer :: i,prevl,maxs
     allocate(zaj_t(maxval(np_g1k),np_e,ista_k:iend_k,kimg));zaj_t=0.d0
     prevl = size(zaj_l,1)
     maxs = maxval(np_g1k)
     if(prevl<maxs) maxs = prevl
     zaj_t(1:maxs,1:np_e,ista_k:iend_k,1:kimg) = zaj_l(1:maxs,1:np_e,ista_k:iend_k,1:kimg)
     if(allocated(zaj_l)) deallocate(zaj_l)
     allocate(zaj_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg)); zaj_l = zaj_t
     if(allocated(zaj_ball)) deallocate(zaj_ball)
     allocate(zaj_ball(maxval(np_g1k),neg,ista_k:iend_k,kimg));
     if(allocated(zah_ball)) deallocate(zah_ball)
     allocate(zah_ball(maxval(np_g1k),neg,kimg));
     deallocate(zaj_t)
     if(sw_keep_hloc_phi==ON) then
       if(allocated(hlocphi_l)) deallocate(hlocphi_l)
       allocate(hlocphi_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg));hlocphi_l = 0.d0
     endif
  end subroutine m_ES_epsmain_reallocate
! ==============================================================================

end module m_Electronic_Structure
