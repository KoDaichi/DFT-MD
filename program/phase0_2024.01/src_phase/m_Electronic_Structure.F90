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
  use m_IterationNumbers,   only : iteration_electronic, nk_in_the_process, nkgroup, iteration
!!!!!BRANCH_P ORG_Parallel
! === DEBUG for icond=2,3 & one_by_one by tkato 2014/01/24 =====================
  use m_IterationNumbers,   only : first_kpoint_in_this_job
! ==============================================================================
!!!!!BRANCH_P_END ORG_Parallel
  use m_NonLocal_Potential, only : snl,norm_phig,AtaulmaG,BtaulmaG
  use m_PlaneWaveBasisSet,  only : igf,kg1,kg,kgp,nbase,nbmx,iba &
       &                         , nbase_gamma,kg1_prev, iba_prev &
       &                         , nbase_prev, nbase_gamma_prev, igf_prev &
       &                         , GVec_on_refcell
  use m_PlaneWaveBasisSet,   only : ngabc, m_pwBS_kinetic_energies
#ifdef __EDA__
  use m_PlaneWaveBasisSet,   only : gr_l, m_pwBS_kinetic_energies
#endif
  use m_PseudoPotential,    only : ival,ilmt,nlmt,nlmtt,nlmta,lmta,lmtt,ltp,mtp,q,dion &
       &                         , ilmt_phi,nlmt_phi,nlmtt_phi,nlmta_phi &
       &                         , lmta_phi,lmtt_phi,ltp_phi,mtp_phi,taup_phi &
       &                         , iproj_phi, nlmta_add, nac_p &
       &                         , modnrm,nac,fqwei,nlmta1,nlmta2 &
       &                         , porb, qorb, m_PP_tell_iorb_lmtt &
       &                         , nrorb, irorb, crorb &
       &                         , ipaw, dion_paw
  use m_Crystal_Structure,  only : op, nopr, nlpnt, additional_charge, altv_refcell
  use m_Kpoints,            only : kv3,vkxyz, kv3_ek,vkxyz_ek,k_symmetry, vkxyz_refcell
  use m_Ionic_System,       only : ntyp,iatom, natm, iwei, ityp, pos, cps, qex &
       &                         , if_pdos, speciesname, iproj_group
  use m_FFT,                only : nfft &
       &                         , m_FFT_alloc_WF_work &
       &                         , m_FFT_dealloc_WF_work &
       &                         , fft_box_size_WF, fft_box_size_WF_prev
  use m_FFT,                only : m_FFT_WF, m_FFT_W_Vlocal_W
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
       &                         , sw_betar_dot_wfs_exp &
       &                         , sw_keep_hloc_phi &
#ifdef SAVE_FFT_TIMES
       &                         , sw_save_fft &
#endif
       &                         , sw_hybrid_functional &
       &                         , sw_interpolate_wfs &
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
       &                         , mpi_e_world &
       &                         , ista_atm, iend_atm &
       &                         , ierr,mp_e,nel_e &
       &                         , ista_g1k,iend_g1k,np_g1k,mp_g1k &
       &                         , np_fs, mp_fs,  nel_fs, nrank_k, nis_kv3_ek &
       &                         , mpi_spin_group, ista_spin, iend_spin
! ============================== added by K. Tagami ========== 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor, ndim_chgpot, ndim_magmom, &
       &                              SpinOrbit_mode, sw_hubbard, sw_band_unfolding, &
       &                              band_unfolding_active
  use m_Const_Parameters,     only : BuiltIn, Neglected

  use m_PseudoPotential,      only : dion_scr_noncl, fqwei_noncl
  use m_FFT,  only : m_FFT_W_Vlocal_W_noncl
  use m_ES_NonCollinear,      only : m_ES_MagMom_to_DensMat_Gspace, &
       &                            m_ES_MagMom_to_DensMat_vlhxcl, &
       &                             m_ES_DensMat_To_MagMom_dm, &
       &                             m_ES_DensMat_To_MagMom_porb
! ============================================================ 11.0

! ============================== added by K. Tagami ========== 12.0Exp
  use m_Control_Parameters,    only :   fixed_charge_k_parallel, sw_band_unfolding
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
  use mpi

  implicit none

  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: hlocphi_l   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l_buf   ! d(kg1,np_e,ista_k:iend_k,kimg) wave functions
  real(kind=DP),allocatable, dimension(:,:,:,:):: zaj_l_prev
#ifdef SAVE_FFT_TIMES
  real(kind=DP),allocatable, dimension(:,:,:)  :: Phifftr_l
  !                                      d(nfft,np_e,ista_k:iend_k) or d(lsize,np_e,ista_k:iend_k)
  integer, allocatable, dimension(:,:) :: status_saved_phifftr ! d(np_e,ista_k:iend_k) 0: nothing or old, 1: stored and new
#endif
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
  integer,        parameter                           :: nblocksize_mgs_default = 8
#endif

  real(kind=DP),private,allocatable,dimension(:)      :: ar, ai

  real(kind=DP),        allocatable, dimension(:)     :: afft, bfft
  integer, private, parameter                         :: sw_timing_2ndlevel = ON

#ifndef NO_NONLOCAL_DGEMM
  logical :: DGEMM_DEBUG = .false.
#endif

  complex(kind=CMPLDP), allocatable, dimension(:) :: vloc_esm

  logical                   :: firstcall = .true.
  integer, allocatable,  dimension(:)   :: fftscnt0, fftrcnt0, fftindex0, fftdist0
  integer, allocatable,  dimension(:,:) :: fftsend0, fftrecv0
  integer :: fftmaxrecv0, fftmaxsend0
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

!  include 'mpif.h'                                      ! MPI
  integer istatus(mpi_status_size)                      ! MPI

  interface m_ES_WF_in_Rspace
    module procedure m_ES_WF_in_Rspace0
    module procedure m_ES_WF_in_Rspace1
  end interface m_ES_WF_in_Rspace

contains


  subroutine m_ES_realloc_zaj()
    if(allocated(zaj_l)) deallocate(zaj_l)
    allocate(zaj_l(kg1,np_e,ista_k:iend_k,kimg));zaj_l = 0.d0
    if(sw_keep_hloc_phi==ON) then
      if(allocated(hlocphi_l)) deallocate(hlocphi_l)
      allocate(hlocphi_l(kg1,np_e,ista_k:iend_k,kimg));zaj_l = 0.d0
    endif
  end subroutine m_ES_realloc_zaj

  subroutine m_ES_alloc_zaj_etc()
    integer :: ik,ib
#ifdef SAVE_FFT_TIMES
    integer :: lsize
#endif
    if(iend_k - ista_k < 0) call phase_error_with_msg(nfout, ' iend_k - ista_k < 0 (in m_ES_alloc_zaj_etc)', &
    __LINE__,__FILE__)
    if(nlmta <= 0) call phase_error_with_msg(nfout, ' nlmta <= 0 (in m_ES_alloc_zaj_etc)',__LINE__,__FILE__)
     allocate(zaj_l(kg1,np_e,ista_k:iend_k,kimg))
     if(sw_keep_hloc_phi==ON) allocate(hlocphi_l(kg1,np_e,ista_k:iend_k,kimg))
     if(sw_rsb==ON) allocate(zaj_l_buf(kg1,np_e,ista_k:iend_k,kimg))
#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       allocate(Phifftr_l(nfft,np_e,ista_k:iend_k)); Phifftr_l = 0.d0
       allocate(status_saved_phifftr(np_e,ista_k:iend_k)); status_saved_phifftr = 0
    end if
#endif
    allocate(neordr(neg,ista_k:iend_k))
    allocate(nrvf_ordr(neg,ista_k:iend_k))
    do ik = ista_k, iend_k
       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)
       nrvf_ordr(1:neg,ik) = (/(ib,ib=1,neg)/)
    end do
    allocate(eko_l(np_e,ista_k:iend_k)); eko_l = 0.d0
    allocate(occup_l(np_e,ista_k:iend_k)); occup_l = 0.d0
    allocate(fsr_l(np_e,nlmta,ista_k:iend_k))
    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
       allocate(fsi_l(np_e,nlmta,ista_k:iend_k)); fsi_l = 0.d0
!!$ASASASAS
!!$    end if
    else
       allocate(fsi_l(1,1,1))
    end if
    fsi_l = 0.d0
!!$ASASASAS
    allocate(vnlph_l(kg1,np_e,kimg))
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

! ======================================= added by K. Tagami =============== 11.0
  subroutine m_ES_alloc_zaj_etc_noncl()
    integer :: ik,ib
    if(iend_k - ista_k < 0) call phase_error_with_msg(nfout, ' iend_k - ista_k < 0 (in m_ES_alloc_zaj_etc_noncl)', &
    __LINE__,__FILE__)
    if(nlmta <= 0) call phase_error_with_msg(nfout, ' nlmta <= 0 (in m_ES_alloc_zaj_etc_noncl)',__LINE__,__FILE__)
    allocate(zaj_l(kg1,np_e,ista_k:iend_k,kimg))
    if(sw_rsb==ON) allocate(zaj_l_buf(kg1,np_e,ista_k:iend_k,kimg))
     if(sw_keep_hloc_phi==ON) then
       if(allocated(hlocphi_l)) deallocate(hlocphi_l)
       allocate(hlocphi_l(maxval(np_g1k),np_e,ista_k:iend_k,kimg));hlocphi_l = 0.d0
     endif
    allocate(neordr(neg,ista_k:iend_k))
    allocate(nrvf_ordr(neg,ista_k:iend_k))
    do ik = ista_k, iend_k
       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)
       nrvf_ordr(1:neg,ik) = (/(ib,ib=1,neg)/)
    end do
    allocate(eko_l(np_e,ista_k:iend_k)); eko_l = 0.d0
    allocate(occup_l(np_e,ista_k:iend_k)); occup_l = 0.d0
    allocate(fsr_l(np_e,nlmta,ista_k:iend_k)); fsr_l = 0.0d0
    allocate(fsi_l(np_e,nlmta,ista_k:iend_k));  fsi_l = 0.d0

    allocate(vnlph_l(kg1,np_e,kimg)); vnlph_l = 0.0d0
    ! -- for PDOS
    if(sw_orb_popu == ON) then
       allocate(compr_l(np_e,nlmta_phi,nopr,ista_k:iend_k))
       allocate(compi_l(np_e,nlmta_phi,nopr,ista_k:iend_k))
       compr_l = 0.d0; compi_l = 0.0d0
    end if
!
    if(sw_use_add_proj == ON) then
       allocate(fsr_add_l(np_e,nlmta_add,ista_k:iend_k))
       allocate(fsi_add_l(np_e,nlmta_add,ista_k:iend_k))
       fsr_add_l = 0.0d0; fsi_add_l = 0.0d0
    end if

  end subroutine m_ES_alloc_zaj_etc_noncl
! =======================================================================-11.0

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
! =================================== modified by K. Tagami ================ 11.0
!!    allocate(vlhxc_l(ista_kngp:iend_kngp,kimg,nspin)); vlhxc_l = 0.d0
!
    if ( noncol ) then
      allocate(vlhxc_l(ista_kngp:iend_kngp,kimg,ndim_magmom))
    else
      if(.not.allocated(vlhxc_l)) allocate(vlhxc_l(ista_kngp:iend_kngp,kimg,nspin))
    endif
    vlhxc_l = 0.d0
! ========================================================================== 11.0

! ====== KT_add ====== 13.0U2
    if ( sw_potential_mixing == ON .or. sw_modified_tfw_functional == ON ) then
      allocate(vlhxc_l_old(ista_kngp:iend_kngp,kimg,ndim_magmom))
      vlhxc_l_old = 0.0d0
   endif
! ==================== 13.0U2

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
    if(nlmt <= 0) call phase_error_with_msg(nfout, ' nlmt <=  0 (in m_ES_alloc_vlhxcQ)',__LINE__,__FILE__)
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
    if(nlmt <= 0) call phase_error_with_msg(nfout, ' nlmt <=  0 (in m_ES_alloc_Dhub)',__LINE__,__FILE__)
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
    if(totch <= 1.d-20) call phase_error_with_msg(nfout, ' ! illegal TOTCH value (m_ES_gtotch)',__LINE__,__FILE__)
    if(totch > natm*500.0) then
       if(printable) then
          do it=1,ntyp
             write(nfout,'(" !! it, ival, iatom, qex = ",i4,f8.4,i4,f8.4)') it,ival(it),iatom(it),qex(it)
          end do
          write(nfout,'(" ! illegal TOTCH value (m_ES_gtotch)")')
       end if
       call phase_error_with_msg(nfout, ' ! illegal TOTCH value (m_ES_gtotch)',__LINE__,__FILE__)
    end if
  end subroutine m_ES_gtotch

  subroutine m_ES_wd_zaj_small_portion(nfout,ik,comment,nc)
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
  end subroutine m_ES_wd_zaj_small_portion

  subroutine m_ES_decide_precon_factor(precon,ik,ib,ekin,p)
    integer, intent(in)                         :: precon,ik,ib
    real(kind=DP), intent(in),  dimension(kg1) :: ekin
    real(kind=DP), intent(out), dimension(kg1) :: p

    integer       :: i
    real(kind=DP) :: ektot, x, x1, x2, d_ektot

! ======================== added by K. Tagami ============== 11.0
    real(kind=DP) :: ctmp1
    integer :: is
! =========================================================== 11.0

    if(precon == ON) then
! ======================== modiifed by K. Tagami ============ 11.0
!       call kinetic_energy(ik,ib,ekin,ektot)   ! -(m_E.S.)

       if ( noncol ) then
          ektot = 0.0d0
          Do is=1, ndim_spinor
             call kinetic_energy( ik+is-1, ib, ekin, ctmp1 )
             ektot = ektot + ctmp1
          End do
       else
          call kinetic_energy(ik,ib,ekin,ektot)   ! -(m_E.S.)
       endif

! =========================================================== 11.0

       d_ektot = 1.d0/ektot
       do i = 1, iba(ik)
          x = ekin(i)*d_ektot
          x1 = 27 + ( 18 + (12 + 8*x) *x) *x
          x2 = 16*(x*x)*(x*x)
          p(i)  = x1/(x1 + x2 )
       end do
    else
       p = 1.d0
    end if
  end subroutine m_ES_decide_precon_factor

  subroutine kinetic_energy(ik,ib,dekin,ektot)
    integer, intent(in) :: ik, ib
    real(kind=DP), intent(in), dimension(kg1) :: dekin
    real(kind=DP), intent(out)                 :: ektot
    integer  :: i
    ektot = 0.d0
    if(kimg == 1) then
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*zaj_l(i,map_z(ib),ik,1)**2
       end do
    else
       do i = 1, iba(ik)
          ektot = ektot + dekin(i)*( zaj_l(i,map_z(ib),ik,1)**2 &
               &                   + zaj_l(i,map_z(ib),ik,2)**2)
       end do
    end if
    if(k_symmetry(ik) == GAMMA)  ektot = ektot*2.d0
  end subroutine kinetic_energy

  subroutine m_ES_sort_eigen_values()
    integer             :: ik, ib, jb
#ifdef _NO_HEAP_SORT_EIGENVALUES_
    integer             :: ibo, jbo
#endif
    real(kind=DP), parameter :: delta = 1.d-12
    real(kind=DP), allocatable,dimension(:) :: eko_t, eko_t2
#ifdef _ESORT_ALLGATHER_
    integer :: n_little_elements, n_big_elements, ns, in, ip, ip0
    real(kind=DP), allocatable,dimension(:) :: work_subst
#endif

    integer                  :: i, id_sname = -1, id_sname2 = -1, id_sname3 = -1
    call tstatc0_begin('m_ES_sort_eigen_values ',id_sname,1)


#ifdef _ESORT_ALLGATHER_
    allocate(eko_t(mp_e*nrank_e))
    allocate(eko_t2(mp_e)); eko_t2 = 0.d0
#else
    allocate(eko_t(neg))
    allocate(eko_t2(np_e))   ! MPI
#endif

    do ik = 1, kv3, af+1
       if(map_k(ik) /= myrank_k) cycle    ! MPI

       do ib = 1, np_e
          eko_t2(ib) = eko_l(ib,ik)
       end do
       if(nrank_e > 1) then
          call tstatc0_begin('mpi_barrier ',id_sname2)
          call mpi_barrier(mpi_k_world(myrank_k),ierr)
          call tstatc0_end(id_sname2)
#ifdef _ESORT_ALLGATHER_
          call tstatc0_begin('mpi_allgather ',id_sname3)
          call mpi_allgather(eko_t2,mp_e, mpi_double_precision &
               &           , eko_t, mp_e, mpi_double_precision,mpi_k_world(myrank_k),ierr)

          n_little_elements = mp_e*nrank_e - neg
          if(ipri >= 2) write(nfout,'("!esort n_little_elements = ",i6)') n_little_elements
          if(ipri >= 2) write(nfout,'("!esort eko_t  = ",10(/,10f8.4))') (eko_t(i),i=1,mp_e*nrank_e)
          if(n_little_elements > 1) then
             n_big_elements = nrank_e - n_little_elements
             ns = mp_e*(n_big_elements + 1) - 1
             allocate(work_subst(neg-ns))
             do in = n_big_elements+2, nrank_e
                ip0 = mp_e*(in-1)
                ip  = (mp_e-1)*(in - n_big_elements - 2)
                do i = 1, nel_e(in-1)
                   work_subst(ip+i) = eko_t(ip0+i)
                end do
             end do
             do i = 1, neg-ns
                eko_t(ns+i-1) = work_subst(i)
             end do
             deallocate(work_subst)
          end if
#else
          call tstatc0_begin('mpi_allgatherv ',id_sname3)
          call mpi_allgatherv(eko_t2,nel_e(myrank_e),mpi_double_precision, eko_t, nel_e &
               & , idisp_e,mpi_double_precision,mpi_k_world(myrank_k),ierr)
#endif
          call tstatc0_end(id_sname3)
       else
          eko_t = eko_t2
       end if
!!$#else
!!$       eko_t = 0.d0                       ! MPI
!!$       do ib = 1, neg                     ! MPI
!!$          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik)   ! MPI
!!$       end do                             ! MPI
!!$
!!$       if(nrank_e > 1) then
!!$          call tstatc0_begin('mpi_allreduce ',id_sname3)
!!$          call mpi_allreduce(eko_t,eko_t2,neg,mpi_double_precision &
!!$               & ,mpi_sum,mpi_k_world(myrank_k),ierr)       ! MPI
!!$          eko_t = eko_t2
!!$          call tstatc0_end(id_sname3)
!!$       end if
!!$#endif

       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)

       if(ipri >=2 ) then
          write(nfout,'(" !Esort  eko_t before sorting ")')
          write(nfout,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
       end if

#ifdef _NO_HEAP_SORT_EIGENVALUES_
       do ib = 1, neg-1
          do jb = ib+1, neg
             ibo = neordr(ib,ik)
             jbo = neordr(jb,ik)
             if(eko_t(jbo)  < eko_t(ibo)-delta) then        ! MPI
                neordr(jb,ik) = ibo
                neordr(ib,ik) = jbo
             endif
          enddo
       enddo
#else
       call heap_sorting(neg,eko_t,neordr(1,ik))
#endif

       do ib = 1, neg
          do jb = 1, neg
             if(ib == neordr(jb,ik)) then
                nrvf_ordr(ib,ik) = jb
                exit
             endif
          enddo
       enddo

       if(ipri >=2 ) then
          write(nfout,'(" !Esort  eko_t in order ")')
          write(nfout,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
       end if
    end do ! do-loop of ik

    if(af /= 0) then
       call cp_eigen_values_for_af()
       call expand_neordr_and_nrvf_ordr()
    end if
    deallocate(eko_t)
!!$#ifdef TRANSPOSE
    deallocate(eko_t2)
!!$#else
!!$    if(nrank_e > 1) deallocate(eko_t2)  ! MPI
!!$#endif
!!$    call mpi_barrier(MPI_CommGroup,ierr)
    call tstatc0_end(id_sname)
  contains
    subroutine cp_eigen_values_for_af()
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle                              ! MPI
         do ib = 1, np_e                                              ! MPI
            eko_l(ib,ik+af) = eko_l(ib,ik)
         enddo
      enddo
    end subroutine cp_eigen_values_for_af

    subroutine expand_neordr_and_nrvf_ordr()
      integer :: ik
      do ik = 1, kv3, af+1
         if(map_k(ik) /= myrank_k) cycle                             ! MPI
         neordr(1:neg,ik+af) = neordr(1:neg,ik)
         nrvf_ordr(1:neg,ik+af) = nrvf_ordr(1:neg,ik)
      end do
    end subroutine expand_neordr_and_nrvf_ordr

#ifndef _NO_HEAP_SORT_EIGENVALUES_
    subroutine heap_sorting(n,eig,iord)
      integer, intent(in) :: n
      real(kind=DP), intent(inout) :: eig(n)
      integer, intent(inout) :: iord(n)

      integer :: i,j,k,itmp,it
      real(kind=DP) :: rt

      !*** initialization: heapfy eig and iord ***
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
      !*** to be fully sort ***
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
    end subroutine heap_sorting
#endif
  end subroutine m_ES_sort_eigen_values

! ===================================== added by K. Tagami ============== 11.0
  subroutine m_ES_sort_eigen_vals_noncl()
    integer             :: ik, ib, jb
#ifdef _NO_HEAP_SORT_EIGENVALUES_
    integer             :: ibo, jbo
#endif
    real(kind=DP), parameter :: delta = 1.d-12
    real(kind=DP), allocatable,dimension(:) :: eko_t, eko_t2
#ifdef _ESORT_ALLGATHER_
    integer :: n_little_elements, n_big_elements, ns, in, ip, ip0
    real(kind=DP), allocatable,dimension(:) :: work_subst
#endif

! ===================================== added by K. Tagami ================ 11.0
    integer :: ikskip, is
! ========================================================================= 11.0

    integer                  :: i, id_sname = -1, id_sname2 = -1, id_sname3 = -1
    call tstatc0_begin('m_ES_sort_eigen_vals_noncl ',id_sname,1)


!!$#ifdef TRANSPOSE
#ifdef _ESORT_ALLGATHER_
    allocate(eko_t(mp_e*nrank_e))
    allocate(eko_t2(mp_e)); eko_t2 = 0.d0
#else
    allocate(eko_t(neg))
    allocate(eko_t2(np_e))   ! MPI
#endif
!!$#else
!!$    allocate(eko_t(neg))
!!$    if(nrank_e > 1) allocate(eko_t2(neg))    ! MPI
!!$#endif

    ikskip = ndim_spinor

    do ik = 1, kv3, ikskip
       if(map_k(ik) /= myrank_k) cycle    ! MPI
!!$#ifdef TRANSPOSE
       do ib = 1, np_e
          eko_t2(ib) = eko_l(ib,ik)
       end do
       if(nrank_e > 1) then
          call tstatc0_begin('mpi_barrier ',id_sname2)
          call mpi_barrier(mpi_k_world(myrank_k),ierr)
          call tstatc0_end(id_sname2)
#ifdef _ESORT_ALLGATHER_
          call tstatc0_begin('mpi_allgather ',id_sname3)
          call mpi_allgather(eko_t2,mp_e, mpi_double_precision &
               &           , eko_t, mp_e, mpi_double_precision,mpi_k_world(myrank_k),ierr)

          n_little_elements = mp_e*nrank_e - neg
          if(ipri >= 2) write(nfout,'("!esort n_little_elements = ",i6)') n_little_elements
          if(ipri >= 2) write(nfout,'("!esort eko_t  = ",10(/,10f8.4))') (eko_t(i),i=1,mp_e*nrank_e)
          if(n_little_elements > 1) then
             n_big_elements = nrank_e - n_little_elements
             ns = mp_e*(n_big_elements + 1) - 1
             allocate(work_subst(neg-ns))
             do in = n_big_elements+2, nrank_e
                ip0 = mp_e*(in-1)
                ip  = (mp_e-1)*(in - n_big_elements - 2)
                do i = 1, nel_e(in-1)
                   work_subst(ip+i) = eko_t(ip0+i)
                end do
             end do
             do i = 1, neg-ns
                eko_t(ns+i-1) = work_subst(i)
             end do
             deallocate(work_subst)
          end if
#else
          call tstatc0_begin('mpi_allgatherv ',id_sname3)
          call mpi_allgatherv(eko_t2,nel_e(myrank_e),mpi_double_precision, eko_t, nel_e &
               & , idisp_e,mpi_double_precision,mpi_k_world(myrank_k),ierr)
#endif
          call tstatc0_end(id_sname3)
       else
          eko_t = eko_t2
       end if
!!$#else
!!$       eko_t = 0.d0                       ! MPI
!!$       do ib = 1, neg                     ! MPI
!!$          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik)   ! MPI
!!$       end do                             ! MPI
!!$
!!$       if(nrank_e > 1) then
!!$          call tstatc0_begin('mpi_allreduce ',id_sname3)
!!$          call mpi_allreduce(eko_t,eko_t2,neg,mpi_double_precision &
!!$               & ,mpi_sum,mpi_k_world(myrank_k),ierr)       ! MPI
!!$          eko_t = eko_t2
!!$          call tstatc0_end(id_sname3)
!!$       end if
!!$#endif

       neordr(1:neg,ik) = (/(ib,ib=1,neg)/)

! ----------------- debug --
!       write(*,'(" !Esort  eko_t before sorting ")')
!       write(*,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
! ------------------------

       if(ipri >=2 ) then
          write(nfout,'(" !Esort  eko_t before sorting ")')
          write(nfout,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
       end if

#ifdef _NO_HEAP_SORT_EIGENVALUES_
       do ib = 1, neg-1
          do jb = ib+1, neg
             ibo = neordr(ib,ik)
             jbo = neordr(jb,ik)
             if(eko_t(jbo)  < eko_t(ibo)-delta) then        ! MPI
                neordr(jb,ik) = ibo
                neordr(ib,ik) = jbo
             endif
          enddo
       enddo
#else
       call heap_sorting(neg,eko_t,neordr(1,ik))
#endif

       do ib = 1, neg
          do jb = 1, neg
             if(ib == neordr(jb,ik)) then
                nrvf_ordr(ib,ik) = jb
                exit
             endif
          enddo
       enddo

       if(ipri >=2 ) then
          write(nfout,'(" !Esort  eko_t in order ")')
          write(nfout,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
       end if

! ------- debug --
!       write(*,'(" !Esort  eko_t in order ")')
!       write(*,'(" !Esort ",10f8.4)') (eko_t(neordr(i,ik)),i=1,neg)
! ---------------

    end do ! do-loop of ik

! ------------------- important- ---------------
!
    Do ik=1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle

       Do is=1, ndim_spinor
          Do ib=1, neg
             neordr( ib,ik+is-1 ) = neordr( ib,ik )
!!!             write(*,*) 'ik, ib ordr : ', ik, ib, neordr( ib,ik )
          End do
       End do
    End do
! --------------------------- unknown
    Do ik=1, kv3, ndim_spinor
       if(map_k(ik) /= myrank_k) cycle

       Do is=1, ndim_spinor
          Do ib=1, neg
              nrvf_ordr( ib,ik+is-1 ) =  nrvf_ordr( ib,ik )
          End do
       End do
    End do
! ---------------------

    deallocate(eko_t)
!!$#ifdef TRANSPOSE
    deallocate(eko_t2)
!!$#else
!!$    if(nrank_e > 1) deallocate(eko_t2)  ! MPI
!!$#endif
!!$    call mpi_barrier(MPI_CommGroup,ierr)
    call tstatc0_end(id_sname)

  contains

    subroutine dummy
    end subroutine dummy

#ifndef _NO_HEAP_SORT_EIGENVALUES_
    subroutine heap_sorting(n,eig,iord)
      integer, intent(in) :: n
      real(kind=DP), intent(inout) :: eig(n)
      integer, intent(inout) :: iord(n)

      integer :: i,j,k,itmp,it
      real(kind=DP) :: rt

      !*** initialization: heapfy eig and iord ***
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
      !*** to be fully sort ***
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
    end subroutine heap_sorting
#endif
  end subroutine m_ES_sort_eigen_vals_noncl
! ========================================================================= 11.0

  subroutine m_ES_energy_eigen_values(nfout,sort)
    integer, intent(in) :: nfout
    logical, intent(in), optional :: sort

    integer             :: is, ik, ipri0
!!$    real(kind=DP), pointer, dimension(:) :: ekin
    real(kind=DP), allocatable, dimension(:) :: ekin
    logical :: sor

    sor = .true.
    if(present(sort)) sor=sort
    allocate(afft(nfft)); allocate(bfft(nfft))
!!$    allocate(sc(kg1))
!!$    ekin => sc(1:kg1)
    allocate(ekin(kg1))
    call m_FFT_alloc_WF_work()

    do is = 1, nspin, af+1
       call m_ES_Vlocal_in_Rspace(is,afft)                  ! (ptfft1)
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle                   ! MPI
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
          call m_ES_eigen_values_for_each_k(is,ik,ekin,afft)
       end do
    end do
    if(sor) call m_ES_sort_eigen_values()

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

    deallocate(afft); deallocate(bfft)
!!$    deallocate(sc)
    deallocate(ekin)
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
  end subroutine m_ES_energy_eigen_values

! =============================== added by K. Tagami ================ 11.0
  subroutine m_ES_energy_eigen_vals_noncl(nfout)
    integer, intent(in) :: nfout
    integer             :: ipri0
    integer :: ik

    real(kind=DP), allocatable, dimension(:) :: ekin
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)

    allocate(ekin(kg1))
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0      ! pot
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0      ! wfn
    call m_FFT_alloc_WF_work()

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )   ! vlhxc_ss -> afft_kt

    Do ik=1, kv3, ndim_spinor
       if (map_k(ik) /= myrank_k) cycle                   ! MPI
       call m_pwBS_kinetic_energies(ik,vkxyz,ekin)
       call m_ES_eigen_vals_each_k_noncl( ik,ekin,afft_kt, bfft_kt )
    End do

    call m_ES_sort_eigen_vals_noncl()

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko(nfout,mode=SCF)

    deallocate( afft_kt, bfft_kt )
    deallocate( ekin)

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
  end subroutine m_ES_energy_eigen_vals_noncl
! ========================================================================= 11.0

  subroutine m_ES_energy_eigen_values_ext(nfout)
    integer, intent(in) :: nfout

    integer             :: is, ik, ipri0
!!$    real(kind=DP), pointer, dimension(:) :: ekin
    real(kind=DP), allocatable, dimension(:) :: ekin

    allocate(afft(nfft)); allocate(bfft(nfft))
!!$    allocate(sc(kg1))
!!$    ekin => sc(1:kg1)
    allocate(ekin(kg1))
    call m_FFT_alloc_WF_work()

    do is = 1, nspin, af+1
       call m_ES_Vlocal_in_Rspace(is,afft)                  ! (ptfft1)
       do ik = is, kv3-nspin+is, nspin
          if(map_k(ik) /= myrank_k) cycle                   ! MPI
          call m_pwBS_kinetic_energies(ik,vkxyz,ekin)

          if(ipri >= 3) then
             write(nfout,'(" -- ekin <<m_ES_eigen_values_ext>> -- ik = ",i8)') ik
             write(nfout,'(8f8.4)') (ekin(ipri0),ipri0=1,iba(ik))
          end if

          call m_ES_eigen_values_for_each_kex(nfout,is,ik,ekin,afft)
       end do
    end do
    call m_ES_sort_eigen_values()

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko2(nfout,mode=SCF)

    deallocate(afft); deallocate(bfft)
!!$    deallocate(sc)
    deallocate(ekin)
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
  end subroutine m_ES_energy_eigen_values_ext

! ==================================== added by K. Tagami ============== 11.0
  subroutine m_ES_energy_eigenvals_ext_noncl(nfout)
    integer, intent(in) :: nfout
    integer             :: ipri0
    integer :: ik

    real(kind=DP), allocatable, dimension(:) :: ekin
    real(kind=DP), allocatable :: afft_kt(:,:)
    real(kind=DP), allocatable :: bfft_kt(:,:)

    allocate(ekin(kg1))
    allocate(afft_kt(nfft,ndim_chgpot)); afft_kt = 0.0d0      ! pot
    allocate(bfft_kt(nfft,ndim_spinor)); bfft_kt = 0.0d0      ! wfn
    call m_FFT_alloc_WF_work()

    call m_ES_Vlocal_in_Rspace_noncl( afft_kt )

    Do ik=1, kv3, ndim_spinor
       if (map_k(ik) /= myrank_k) cycle                   ! MPI
       call m_pwBS_kinetic_energies(ik,vkxyz,ekin)

       if(ipri >= 3) then
          write(nfout,'(" -- ekin <<m_ES_eigen_values_ext>> -- ik = ",i8)') ik
          write(nfout,'(8f8.4)') (ekin(ipri0),ipri0=1,iba(ik))
       end if

       call m_ES_eigen_vals_each_kex_noncl( nfout, ik, ekin, afft_kt, bfft_kt )
    End do

    call m_ES_sort_eigen_vals_noncl()

    call get_ipri0(iprieigenvalue,ipri0)
    if(ipri0 >= 2) call m_ES_wd_eko2(nfout,mode=SCF)

    deallocate( afft_kt, bfft_kt )
    deallocate(ekin)
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
  end subroutine m_ES_energy_eigenvals_ext_noncl
! ========================================================================== 11.0

  subroutine m_ES_init_chgsoft()
    if(.not.allocated(chg_softpart)) allocate(chg_softpart(nfft,nspin));chg_softpart=0.d0
  end subroutine m_ES_init_chgsoft

  subroutine m_ES_dealloc_chgsoft()
    if(allocated(chg_softpart)) deallocate(chg_softpart)
    chg_has_been_calculated = .false.
  end subroutine m_ES_dealloc_chgsoft

  subroutine m_ES_eigen_values_for_each_k(is,ik,ekin,afft,eval_charge)
    integer,       intent(in)                  :: is,ik    ! is: spin
    real(kind=DP), intent(in), dimension(kg1)  :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: afft
    logical,       intent(in), optional        :: eval_charge
#ifdef NEC_TUNE_SMP
    real(kind=DP), dimension(nfft) :: bfft
#endif
    logical       :: chg
    integer       :: ib
    real(kind=DP) :: eg
    integer       :: id_sname = -1
!!$    write(nfout,'(" --- energy_eigen_vlues ---")')
    call tstatc0_begin('energy_eigen_values ', id_sname,1)

    if(ipri >= 3) then
       write(nfout,'(" -- eko_l before <<m_ES_eigen_values_for_each_k>> --, ik = ",i3)') ik
       write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
    endif
    chg = .false.
    if(present(eval_charge)) chg = eval_charge
    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(bfft,eg)
#endif
    do ib = ista_e, iend_e, istep_e  ! MPI
       call m_ES_WF_in_Rspace(ik,ib,bfft) ! -(m_E.S.); (swffft)
       if(chg) then
         call add_occupied_densities(bfft,ib,is,ik)
       endif
       call m_FFT_W_Vlocal_W(ELECTRON,nfft,afft,bfft,eg) ! (eigens) --> eg
       eko_l(map_z(ib),ik) = eko_l(map_z(ib),ik) + eg ! MPI
    end do

    if(ipri >= 3) then
       write(nfout,'(" -- eko_l before W_Vnonlocal_W <<m_ES_eigen_values_for_each_k>> --, ik = ",i3)') ik
       write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
    endif

    call W_Vnonlocal_W(ik)

    if(chg) chg_has_been_calculated = .true.

    call tstatc0_end(id_sname)
  contains
    subroutine W_T_W()
      integer :: ib, i
      real(kind=DP) :: eko

      if(sw_hybrid_functional == OFF) eko_l(1:np_e,ik) = 0.d0                   ! MPI
      if(kimg==1) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = 1, np_e                        ! MPI
            eko = 0.d0
            do i = 1, iba(ik)
               !eko_l(ib,ik) = eko_l(ib,ik) + ekin(i)*zaj_l(i,ib,ik,1)**2
               eko = eko + ekin(i)*zaj_l(i,ib,ik,1)**2
            end do
            if(k_symmetry(ik) == GAMMA) eko = 2.d0*eko
            eko_l(ib,ik) = eko_l(ib,ik) + eko
         end do
      else if(kimg==2) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = 1, np_e                        ! MPI
            eko = 0.d0
            do i = 1, iba(ik)
               !eko_l(ib,ik) = eko_l(ib,ik) + ekin(i)*(zaj_l(i,ib,ik,1)**2+zaj_l(i,ib,ik,2)**2)
               eko = eko + ekin(i)*(zaj_l(i,ib,ik,1)**2+zaj_l(i,ib,ik,2)**2)
            end do
            if(k_symmetry(ik) == GAMMA) eko = 2.d0*eko
            eko_l(ib,ik) = eko_l(ib,ik) + eko
         end do
      end if

      if(ipri >= 2) then
         write(6,'(" -- eko_l (W_T_W) --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)
      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib
      real(kind=DP) :: fac

      integer :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)

      do ia = 1, natm
         it = ityp(ia)
         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)
            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
!!$                  fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = 1, np_e                                  ! MPI
                     eko_l(ib,ik) = eko_l(ib,ik) &
                          & + fac*(fsr_l(ib,p,ik)*fsr_l(ib,q,ik))
                  end do
               end do
            else
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac = iwei(ia)
!!$                  fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  if(ipaw(it).eq.0) then
                     fac   = fac*(dion(lmt1,lmt2,it) + vlhxcQ(lmt1,lmt2,ia,is))
                  else
                     fac   = fac*(dion_paw(lmt1,lmt2,is,ia) + vlhxcQ(lmt1,lmt2,ia,is))
                  end if
                  do ib = 1, np_e                                  ! MPI
                     eko_l(ib,ik) = eko_l(ib,ik) &
                          & + fac*(fsr_l(ib,p,ik)*fsr_l(ib,q,ik)&
                          &      + fsi_l(ib,p,ik)*fsi_l(ib,q,ik))
                  end do
               end do
            end if
         end do
      end do
      if(ipri >= 2) then
         write(6,'(" -- eko_l <<W_Vnonlocal_W>> --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif

      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W

    subroutine add_occupied_densities(bfft,ib1,is,ik)
      real(kind=DP), intent(in), dimension(nfft) :: bfft
      integer, intent(in) :: ib1,is,ik
      integer :: i
      real(kind=DP) :: occupation
      occupation = occup_l(map_z(ib1),ik)
      do i = 1, nfft-1, 2
         chg_softpart(i,is) = chg_softpart(i,is) + occupation*(bfft(i)**2+bfft(i+1)**2) ! MPI
      end do
    end subroutine add_occupied_densities

  end subroutine m_ES_eigen_values_for_each_k

! ============================= added by K. Tagami ======================= 11.0
  subroutine m_ES_eigen_vals_each_k_noncl( ik, ekin, afft_kt, bfft_kt )
    integer,       intent(in)                  :: ik
    real(kind=DP), intent(in), dimension(kg1) :: ekin
    real(kind=DP), intent(in)  :: afft_kt( nfft,ndim_chgpot )
!
    real(kind=DP), intent(inout)  :: bfft_kt( nfft,ndim_spinor )

    integer       :: ib, is
    real(kind=DP) :: eg
    integer       :: id_sname = -1
!!$    write(nfout,'(" --- energy_eigen_vlues ---")')

    call tstatc0_begin('energy_eigen_values_k_noncl ', id_sname,1)

    if(ipri >= 3) then
       write(nfout,'(" -- eko_l before <<m_ES_eigen_vals_each_k_noncl>> --, ik = ",i3)') ik
       write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
    endif

    call W_T_W_noncl()                     ! (eigen1) --> eko_l , kinetic part

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(bfft,eg)
#endif
    do ib = ista_e, iend_e, istep_e  ! MPI
       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik +is -1  , ib, bfft_kt(1:nfft,is) )
       End do

       call m_FFT_W_Vlocal_W_noncl( ELECTRON, nfft, afft_kt, bfft_kt, eg, &
	&                           ndim_spinor, ndim_chgpot ) ! (eigens) --> eg
       eko_l(map_z(ib),ik) = eko_l(map_z(ib),ik) + eg ! MPI

    end do

    if(ipri >= 3) then
       write(nfout,'(" -- eko_l before W_Vnonlocal_W <<m_ES_eigen_values_for_each_k>> --, ik = ",i3)') ik
       write(nfout,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
    endif


!    write(938,*) iwei
!    write(940,*) fsr_l
!    write(942,*) dion_Scr_noncl
!    stop

    call W_Vnonlocal_W_noncl_A(ik)

!    stop

    call tstatc0_end(id_sname)

  contains

    subroutine W_T_W_noncl()
      integer :: ib, i, is
      real(kind=DP) :: ctmp, meko

      if(sw_hybrid_functional == OFF) eko_l(1:np_e,ik) = 0.d0           ! MPI
      if(kimg==1) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = 1, np_e                        ! MPI
            meko = 0.0d0
            do i = 1, iba(ik)
	       ctmp = 0.0d0
               Do is=1, ndim_spinor
                  ctmp = ctmp + zaj_l( i, ib, ik +is-1, 1 )**2
               End do
               meko = meko + ekin(i) *ctmp
!!!!!               eko_l(ib,ik) = eko_l(ib,ik) + ekin(i) *ctmp
            end do
!!!!            if(k_symmetry(ik) == GAMMA) eko_l(ib,ik) = 2.d0*eko_l(ib,ik)

            if ( k_symmetry(ik) == GAMMA ) meko = 2.0d0 *meko
            eko_l(ib,ik) = eko_l(ib,ik) + meko

         end do

      else if(kimg==2) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = 1, np_e                        ! MPI
            meko = 0.0d0
            do i = 1, iba(ik)
	       ctmp = 0.0d0
               Do is=1, ndim_spinor
                  ctmp = ctmp + zaj_l( i, ib, ik +is-1, 1 )**2 &
	&                     + zaj_l( i, ib, ik +is-1, 2 )**2
	       End do
!!               eko_l(ib,ik) = eko_l(ib,ik) + ekin(i) *ctmp
               meko = meko + ekin(i) *ctmp
            end do
!!            if(k_symmetry(ik) == GAMMA) eko_l(ib,ik) = 2.d0*eko_l(ib,ik)

            if ( k_symmetry(ik) == GAMMA ) meko = 2.0d0 *meko
            eko_l(ib,ik) = eko_l(ib,ik) + meko

         end do
      end if

      if(ipri >= 2) then
         write(6,'(" -- eko_l (W_T_W) --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W_noncl

    subroutine W_Vnonlocal_W_noncl_A(ik)
      integer, intent(in) :: ik

      integer :: ia, lmt1, lmt2, it, p, q, ib
      integer :: is1, is2, is_tmp
      integer :: il1, il2, im1, im2
      integer :: k1, k2
      integer :: mdvdb

      real(kind=DP) :: c1, c2, fac0
      complex(kind=CMPLDP) :: fac

      real*8 csum_1, csum_2

      integer :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W_noncl ',id_sname)

!      write(960,*) 'FF'
!      write(960,*) dion_scr_noncl

      csum_1 = 0.0d0;       csum_2 = 0.0d0

      do ia = 1, natm
         it = ityp(ia)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
         mdvdb = m_PP_include_vanderbilt_pot(it)
#endif
! ---------------------------------------------------- 11.0S

         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)

            il1 = ltp(lmt1,it); im1 = mtp(lmt1,it)

            if(k_symmetry(ik) == GAMMA) then
               do lmt2 = lmt1, ilmt(it)
                  q = lmta(lmt2,ia)
                  fac0   = 2.d0 * iwei(ia)
                  if(lmt1 == lmt2) fac0 = iwei(ia)

                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ndim_spinor*( is1 -1 ) +is2
                        fac = fac0 *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )

                        k1 = ik + is1 - 1
                        k2 = ik + is2 - 1

                        Do ib = 1, np_e                                  ! MPI
                           eko_l(ib,ik) = eko_l(ib,ik) &
                                & + fac*(fsr_l(ib,p,k1)*fsr_l(ib,q,k2) )
                        End do
                     End do
                  End do
               End do
            else
! ----------------------------------- Which is better ? --------
!!!!!               Do lmt2 = lmt1, ilmt(it)
               Do lmt2 = 1, ilmt(it)
                  q = lmta(lmt2,ia)

                  il2 = ltp(lmt2,it); im2 = mtp(lmt2,it)


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

!                  fac   = 2.d0 * iwei(ia)
!                  if(lmt1 == lmt2) fac = iwei(ia)

                  fac0 = iwei(ia)

                  Do is1=1, ndim_spinor
                     Do is2=1, ndim_spinor
                        is_tmp = ndim_spinor *( is1 -1 ) +is2
                        fac = fac0 *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )

!                        write(950,*) lmt1, lmt2, is_tmp, ia, fac, dion( lmt1, lmt2, it )

                        k1 = ik + is1 - 1
                        k2 = ik + is2 - 1

                        Do ib = 1, np_e                                  ! MPI
                           c1 =  real(fac) *( fsr_l(ib,p,k1)*fsr_l(ib,q,k2) &
                                &            +fsi_l(ib,p,k1)*fsi_l(ib,q,k2) )
#if 1
! -- orig -
                           c2 = -aimag(fac) *( fsr_l(ib,p,k1)*fsi_l(ib,q,k2) &
	&                                     -fsi_l(ib,p,k1)*fsr_l(ib,q,k2) )
#else
! -- kt : debug --
!                           c2 = aimag(fac) *( fsr_l(ib,p,k1)*fsi_l(ib,q,k2) &
!                                &            -fsi_l(ib,p,k1)*fsr_l(ib,q,k2) )
!--
#endif
                           eko_l(ib,ik) = eko_l(ib,ik) &
                                & + c1 + c2
                        End do
                     End Do
                  End Do
               end do
            end if

         end do
      end do

      if(ipri >= 2) then
         write(6,'(" -- eko_l <<W_Vnonlocal_W_noncl A *>> --, ik = ",i3)') ik
         write(6,'(8f8.4)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif

      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W_noncl_A

  end subroutine m_ES_eigen_vals_each_k_noncl
! ==================================================================== 11.0

  subroutine m_ES_eigen_values_for_each_kex(nfout,is,ik,ekin,afft)
    integer,       intent(in)                  :: nfout,is,ik    ! is: spin
    real(kind=DP), intent(in), dimension(kg1) :: ekin
    real(kind=DP), intent(in), dimension(nfft) :: afft

    integer       :: ib
    real(kind=DP) :: eg
    integer,save  :: id_sname = -1
    call tstatc0_begin('energy_eigen_values ', id_sname,1)

    call W_T_W()                     ! (eigen1) --> eko_l , kinetic part
    do ib = ista_e, iend_e, istep_e  ! MPI
       if(nrvf_ordr(ib,ik) <= neg_previous) cycle
       call m_ES_WF_in_Rspace(ik,ib,bfft) ! -(m_E.S.); (swffft)
       call m_FFT_W_Vlocal_W(ELECTRON,nfft,afft,bfft,eg) ! (eigens) --> eg
       eko_l(map_z(ib),ik) = eko_l(map_z(ib),ik) + eg ! MPI
    end do

    call W_Vnonlocal_W(ik)

    call tstatc0_end(id_sname)
  contains
    subroutine W_T_W()
      integer :: ib, i, ip, ibt

!!$      eko_l(1:np_e,ik) = 0.d0                   ! MPI
      do ib = neg_previous+1, neg
         ip = neordr(ib,ik)
         if(map_e(ip) == myrank_e) eko_l(map_z(ip),ik) = 0.d0
      end do

      if(kimg==1) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = neg_previous+1, neg                   ! MPI
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, iba(ik)
                  eko_l(ibt,ik) = eko_l(ibt,ik) + ekin(i)*zaj_l(i,ibt,ik,1)**2
               end do
               if(k_symmetry(ik) == GAMMA) eko_l(ibt,ik) = eko_l(ibt,ik)*2.d0
            end if
         end do
      else if(kimg==2) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
!!$         do ib = 1, np_e                        ! MPI
         do ib = neg_previous+1, neg
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, iba(ik)
                  eko_l(ibt,ik) = eko_l(ibt,ik) + ekin(i)*(zaj_l(i,ibt,ik,1)**2+zaj_l(i,ibt,ik,2)**2)
               end do
               if(k_symmetry(ik) == GAMMA) eko_l(ibt,ik) = eko_l(ibt,ik)*2.d0
            end if
         end do
      end if

      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_T_W) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W

    subroutine W_Vnonlocal_W(ik)
      integer, intent(in) :: ik

      integer       :: ia, lmt1, lmt2, it, p, q, ib, ibo
      real(kind=DP) :: fac

      integer,save :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W ',id_sname)

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
               if(k_symmetry(ik) == GAMMA) then
                  do ib = 1, np_e                                  ! MPI
                     ibo = ista_e + ista_e*(ib-1)
                     if(nrvf_ordr(ibo,ik) <= neg_previous) cycle
                     eko_l(ib,ik) = eko_l(ib,ik) &
                          & + fac*(fsr_l(ib,p,ik)*fsr_l(ib,q,ik))
                  end do
               else
                  do ib = 1, np_e                                  ! MPI
                     ibo = ista_e + ista_e*(ib-1)
                     if(nrvf_ordr(ibo,ik) <= neg_previous) cycle
                     eko_l(ib,ik) = eko_l(ib,ik) &
                          & + fac*(fsr_l(ib,p,ik)*fsr_l(ib,q,ik)&
                          &      + fsi_l(ib,p,ik)*fsi_l(ib,q,ik))
                  end do
               end if
            end do
         end do
      end do
      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_Vnonlocal_W) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif

      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W
  end subroutine m_ES_eigen_values_for_each_kex

!================================ added by K. Tagami ==================== 11.0
  subroutine m_ES_eigen_vals_each_kex_noncl( nfout,ik,ekin,afft_kt, bfft_kt )
    integer,       intent(in)                  :: nfout, ik
    real(kind=DP), intent(in), dimension(kg1) :: ekin

    real(kind=DP), intent(in)  :: afft_kt( nfft,ndim_chgpot )
    real(kind=DP), intent(inout)  :: bfft_kt( nfft,ndim_spinor )

    integer       :: ib, is
    real(kind=DP) :: eg
    integer,save  :: id_sname = -1

    call tstatc0_begin('energy_eigen_vals_kex_noncl ', id_sname,1)

    call W_T_W_noncl()                     ! (eigen1) --> eko_l , kinetic part

    do ib = ista_e, iend_e, istep_e  ! MPI
       if(nrvf_ordr(ib,ik) <= neg_previous) cycle

       Do is=1, ndim_spinor
         call m_ES_WF_in_Rspace( ik +is-1, ib, bfft_kt(:,is) )
       End do

       call m_FFT_W_Vlocal_W_noncl( ELECTRON,nfft,afft_kt,bfft_kt,eg, &
            &                       ndim_spinor, ndim_chgpot )  ! (eigens) --> eg

       eko_l(map_z(ib),ik) = eko_l(map_z(ib),ik) + eg ! MPI
    end do

    call W_Vnonlocal_W_noncl_A(ik)

    call tstatc0_end(id_sname)
  contains

    subroutine W_T_W_noncl()
      integer :: ib, i, ip, ibt, is
      real(kind=DP) :: ctmp

!!$      eko_l(1:np_e,ik) = 0.d0                   ! MPI
      do ib = neg_previous+1, neg
         ip = neordr(ib,ik)
         if(map_e(ip) == myrank_e) eko_l(map_z(ip),ik) = 0.d0
      end do

      if(kimg==1) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
         do ib = neg_previous+1, neg                   ! MPI
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, iba(ik)
	          ctmp = 0.0d0
         	  Do is=1, ndim_spinor
	            ctmp = ctmp + zaj_l( i, ibt, ik +is-1, 1 )**2
	          End do
                  eko_l(ibt,ik) = eko_l(ibt,ik) + ekin(i) *ctmp
               end do
               if ( k_symmetry(ik)==GAMMA ) eko_l(ibt,ik) = eko_l(ibt,ik)*2.d0
            end if
         end do
      else if(kimg==2) then
#ifdef VPP
*vocl loop,unroll(4)
#endif
!!$         do ib = 1, np_e                        ! MPI
         do ib = neg_previous+1, neg
            ip = neordr(ib,ik)
            if(map_e(ip) == myrank_e) then
               ibt = map_z(ip)
               do i = 1, iba(ik)
	          ctmp = 0.0d0
	          Do is=1, ndim_spinor
                    ctmp = ctmp + zaj_l( i, ib, ik +is-1, 1 )**2 &
	&                       + zaj_l( i, ib, ik +is-1, 2 )**2
	          End do
                  eko_l(ibt,ik) = eko_l(ibt,ik) + ekin(i) *ctmp
               end do
               if ( k_symmetry(ik)==GAMMA ) eko_l(ibt,ik) = eko_l(ibt,ik)*2.d0
            end if
         end do
      end if

      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_T_W_noncl) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif
    end subroutine W_T_W_noncl

    subroutine W_Vnonlocal_W_noncl_A(ik)
      integer, intent(in) :: ik

      integer :: ia, lmt1, lmt2, it, p, q, ib, ibo
      integer :: is1, is2, is_tmp
      integer :: il1, il2, im1, im2
      integer :: k1, k2
      integer :: mdvdb

      real(kind=DP) :: c1, c2, fac0
      complex(kind=CMPLDP) :: fac

      integer,save :: id_sname = -1
      call tstatc0_begin('W_Vnonlocal_W_noncl ',id_sname)

!      write(960,*) 'FF'
!      write(960,*) dion_scr_noncl

      do ia = 1, natm
         it = ityp(ia)

! ---------------------------------------------------- 11.0S
#ifdef SKIP_TEST
         mdvdb = m_PP_include_vanderbilt_pot(it)
#endif
! ---------------------------------------------------- 11.0S

         do lmt1 = 1, ilmt(it)
            p = lmta(lmt1,ia)

            il1 = ltp(lmt1,it);  im1 = mtp(lmt1,it)

!!!!!!            do lmt2 = lmt1, ilmt(it)
            do lmt2 = 1, ilmt(it)
               q = lmta(lmt2,ia)

               il2 = ltp(lmt2,it);  im2 = mtp(lmt2,it)

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

!               fac0   = 2.d0 * iwei(ia)
!               if(lmt1 == lmt2) fac0 = iwei(ia)

               fac0   = iwei(ia)

               Do is1=1, ndim_spinor
                  Do is2=1, ndim_spinor
                     is_tmp = ndim_spinor*( is1 -1 ) +is2
                     fac = fac0 *dion_scr_noncl( lmt1, lmt2, is_tmp, ia )

                     k1 = ik + is1 - 1
                     k2 = ik + is2 - 1

                     if ( k_symmetry(ik) == GAMMA ) then
                        do ib = 1, np_e                                  ! MPI
                           ibo = ista_e + ista_e*(ib-1)
                           if(nrvf_ordr(ibo,ik) <= neg_previous) cycle
                           eko_l(ib,ik) = eko_l(ib,ik) &
                                & + fac*( fsr_l(ib,p,k1)*fsr_l(ib,q,k2) )
                        end do
                     else
                        do ib = 1, np_e                                  ! MPI
                           ibo = ista_e + ista_e*(ib-1)
                           if(nrvf_ordr(ibo,ik) <= neg_previous) cycle

                           c1 =  real(fac) *( fsr_l(ib,p,k1)*fsr_l(ib,q,k2) &
                                &            +fsi_l(ib,p,k1)*fsi_l(ib,q,k2) )
#if 1
! ----------- orig --
                           c2 = -aimag(fac) *( fsr_l(ib,p,k1)*fsi_l(ib,q,k2) &
                                &	     -fsi_l(ib,p,k1)*fsr_l(ib,q,k2) )
#else
! --- kt :debug --
!                           c2 = aimag(fac) *( fsr_l(ib,p,k1)*fsi_l(ib,q,k2) &
!                                &	     -fsi_l(ib,p,k1)*fsr_l(ib,q,k2) )
! ---
#endif
                           eko_l(ib,ik) = eko_l(ib,ik) + c1 +c2

                        end do
                     end if
                  End Do      ! is2
               End do       !is1

            end do
         end do
      end do
      if(ipri >= 3 .and. printable) then
         write(nfout,'(" -- eko_l (W_Vnonlocal_W_noncl :) --, ik = ",i3)') ik
         write(nfout,'(5d20.8)') (eko_l(ib,ik),ib=1,np_e) ! MPI
      endif

      call tstatc0_end(id_sname)
    end subroutine W_Vnonlocal_W_noncl_A

  end subroutine m_ES_eigen_vals_each_kex_noncl
! ========================================================= 11.0

! =============
  subroutine m_ES_Vlocal_in_Rspace(is,afft,pot_l)
    integer, intent(in) :: is
    real(kind=DP), dimension(nfft) :: afft
    real(kind=DP), intent(in), optional :: pot_l(ista_kngp:iend_kngp,kimg,ndim_magmom)

    integer,save                   :: id_sname = -1
!!$    call tstatc0_begin('m_ES_Vlocal_in_Rspace ', id_sname,1)
    call tstatc0_begin('m_ES_Vlocal_in_Rspace ', id_sname)

    call map_vlhxc_l_onto_afft(is)      !-(contained here) vlhxc_l  --> afft ;using (igf)
    call m_FFT_WF(ELECTRON,nfout,afft,INVERSE,OFF) ! afft -> afft

    call tstatc0_end(id_sname)
  contains
    subroutine map_vlhxc_l_onto_afft(is)
      integer, intent(in) :: is
      integer :: i,i1,ri
      integer :: iend
      real(kind=DP), pointer, dimension(:) :: afft_mpi
      if(npes > 1) allocate(afft_mpi(nfft))

      afft = 0.d0
      iend = iend_kngp
      if( iend > kg ) iend = kg
      if ( present(pot_l) ) then
         do ri = 1, kimg
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
               afft(i1) = pot_l(i,ri,is)
            end do
         end do
      else
         do ri = 1, kimg
            do i = ista_kngp, iend  !for mpi
               i1 = kimg*igf(i) + (ri - kimg)
               afft(i1) = vlhxc_l(i,ri,is)
            end do
         end do
      endif

      if(npes > 1) then
         call mpi_allreduce(afft,afft_mpi,nfft &
              &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         afft = afft_mpi
         deallocate(afft_mpi)
      end if

    end subroutine map_vlhxc_l_onto_afft
  end subroutine m_ES_Vlocal_in_Rspace

! ==================================== added by K. Tagami ================== 11.0
  subroutine m_ES_Vlocal_in_Rspace_noncl( afft_kt )
    real(kind=DP), intent(out) :: afft_kt( nfft, ndim_chgpot )

    integer :: is
    integer,save                   :: id_sname = -1
!!$    call tstatc0_begin('m_ES_Vlocal_in_Rspace_noncl ', id_sname,1)
    call tstatc0_begin('m_ES_Vlocal_in_Rspace_noncl ', id_sname)

    allocate( vlhxc_ssrep( ista_kngp:iend_kngp,kimg,ndim_chgpot) )
    vlhxc_ssrep = 0.0d0

    call m_ES_MagMom_to_DensMat_vlhxcl( vlhxc_l, vlhxc_ssrep )

    afft_kt = 0.d0
    Do is=1, ndim_chgpot
      call map_vlhxc_ss_onto_afft(is)
                    !-(contained here) vlhxc_l  --> afft ;using (igf)
      call m_FFT_WF( ELECTRON, nfout, afft_kt(:,is), INVERSE, OFF ) ! afft -> afft
    End do

    deallocate( vlhxc_ssrep )
    call tstatc0_end(id_sname)

  contains

    subroutine map_vlhxc_ss_onto_afft(is)
      integer, intent(in) :: is
      integer :: i,i1,ri
      integer :: iend
      real(kind=DP), allocatable, dimension(:) :: afft_mpi

      if(npes > 1) allocate(afft_mpi(nfft))

      iend = iend_kngp
      if( iend > kg ) iend = kg
      do ri = 1, kimg
         do i = ista_kngp, iend  !for mpi
            i1 = kimg*igf(i) + (ri - kimg)
            afft_kt(i1,is) = vlhxc_ssrep(i,ri,is)
         end do
      end do

      if(npes > 1) then
         call mpi_allreduce( afft_kt(:,is), afft_mpi,nfft &
              &  ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
         afft_kt(:,is) = afft_mpi(:)
         deallocate(afft_mpi)
      end if

    end subroutine map_vlhxc_ss_onto_afft

  end subroutine m_ES_Vlocal_in_Rspace_noncl
! ========================================================================= 11.0

  subroutine m_ES_WF_in_Rspace_kt(k1,k2,ik,psi_l,bfft)
    integer, intent(in) :: k1,k2,ik
    real(kind=DP), intent(in),dimension(kg1,1,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    if ( sw_band_unfolding == ON .and. band_unfolding_active ) then
       call m_ES_map_WF_on_fftmesh_unfold(k1,k2,ik,psi_l,bfft)
    else
       call m_ES_map_WF_on_fftmesh(k1,k2,ik,psi_l,bfft)
    endif
    call m_FFT_WF(ELECTRON,nfout,bfft,INVERSE,ON)

  end subroutine m_ES_WF_in_Rspace_kt

  subroutine m_ES_map_WF_on_fftmesh(k1,k2,ik,psi_l,bfft)
    integer, intent(in) :: k1,k2,ik
    real(kind=DP), intent(in),dimension(kg1,1,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf(1)
          bfft(i1) = psi_l(1,1,ik,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = igf(i)
             bfft(i1) = psi_l(ii,1,ik,1)
             j = nbase_gamma(ii,2)
             i2 = igf(j)
             bfft(i2) =   psi_l(ii,1,ik,1)
          end do
       else if(kimg == 2) then
          i1 = 2*igf(1) - 1
          bfft(i1)   = psi_l(1,1,ik,1)
          bfft(i1+1) = psi_l(1,1,ik,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = 2*igf(i)-1
             bfft(i1  ) = psi_l(ii,1,ik,1)
             bfft(i1+1) = psi_l(ii,1,ik,2)
             j = nbase_gamma(ii,2)
             i2 = 2*igf(j)-1
             bfft(i2  ) = psi_l(ii,1,ik,1)
             bfft(i2+1) = -psi_l(ii,1,ik,2)
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
             bfft(i1) = psi_l(i,1,ik,ri)   ! MPI
          end do
       end do
    end if
  end subroutine m_ES_map_WF_on_fftmesh

  subroutine m_ES_map_WF_on_fftmesh_unfold(k1,k2,ik,psi_l,bfft)
    integer, intent(in) :: k1,k2,ik
    real(kind=DP), intent(in),dimension(kg1,1,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf(1)
          bfft(i1) = psi_l(1,1,ik,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             if ( sw_force_kpt_inside_bz == ON ) then
                if ( GVec_on_refcell(i,ik) == 0 ) cycle
             else
                if ( GVec_on_refcell(i,1) == 0 ) cycle
             endif
             i1 = igf(i)
             bfft(i1) = psi_l(ii,1,ik,1)
          end do

          do ii = 2, iba(ik)
             j = nbase_gamma(ii,2)
             if ( sw_force_kpt_inside_bz == ON ) then
                if ( GVec_on_refcell(j,ik) == 0 ) cycle
             else
                if ( GVec_on_refcell(j,1) == 0 ) cycle
             endif
             i2 = igf(j)
             bfft(i2) =   psi_l(ii,1,ik,1)
          end do

       else if(kimg == 2) then
          i1 = 2*igf(1) - 1
          bfft(i1)   = psi_l(1,1,ik,1)
          bfft(i1+1) = psi_l(1,1,ik,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             if ( sw_force_kpt_inside_bz == ON ) then
                if ( GVec_on_refcell(i,ik) == 0 ) cycle
             else
                if ( GVec_on_refcell(i,1) == 0 ) cycle
             endif
             i1 = 2*igf(i)-1
             bfft(i1  ) = psi_l(ii,1,ik,1)
             bfft(i1+1) = psi_l(ii,1,ik,2)
          end do
          do ii = 2, iba(ik)
             j = nbase_gamma(ii,2)
             if ( sw_force_kpt_inside_bz == ON ) then
                if ( GVec_on_refcell(j,ik) == 0 ) cycle
             else
                if ( GVec_on_refcell(j,1) == 0 ) cycle
             endif
             i2 = 2*igf(j)-1
             bfft(i2  ) = psi_l(ii,1,ik,1)
             bfft(i2+1) = -psi_l(ii,1,ik,2)
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
             ii = nbase(i,ik)
             if ( sw_force_kpt_inside_bz == ON ) then
                if ( GVec_on_refcell(ii,ik) == 0 ) cycle
             else
                if ( GVec_on_refcell(ii,1) == 0 ) cycle
             endif
             i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
             bfft(i1) = psi_l(i,1,ik,ri)   ! MPI
          end do
       end do
    end if
  end subroutine m_ES_map_WF_on_fftmesh_unfold

  subroutine m_ES_WF_in_Rspace1(k1,k2,ik,ib,psi_l,bfft)
    integer, intent(in) :: k1,k2,ik,ib
    real(kind=DP), intent(in),dimension(kg1,np_e,k1:k2,kimg) :: psi_l
    real(kind=DP), intent(inout), dimension(nfft) :: bfft

    integer :: i,i1,ri, j, i2, ii
#ifndef NEC_TUNE_SMP
    integer,save  :: id_sname = -1
#ifdef SAVE_FFT_TIMES
    integer,save :: id_sname2 = -1
#endif
#endif

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       if (status_saved_phifftr(map_z(ib),ik) == STORED_AND_NEW) then
#ifndef NEC_TUNE_SMP
          call tstatc0_begin('m_ES_WF_in_Rspace(2) ',id_sname2,1)
#endif
          if(ipri>=2 .and. ik==1 .and. ib==1) write(nfout,'(" !### zaj_fftr(stored) --> bfft")')
          bfft(:) = Phifftr_l(:,map_z(ib),ik)
#ifndef NEC_TUNE_SMP
          call tstatc0_end(id_sname2)
#endif
       endif
    else
#endif

#ifndef NEC_TUNE_SMP
       call tstatc0_begin('m_ES_WF_in_Rspace(1) ',id_sname,1)
#endif

    if(ipri>=2 .and. ik==1 .and. ib==1) write(nfout,'(" !### zaj_l --(FFT)--> bfft")')

    bfft = 0.d0
    if(k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          i1 = igf(1)
          bfft(i1) = psi_l(1,map_z(ib),ik,1)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = igf(i)
             bfft(i1) = psi_l(ii,map_z(ib),ik,1)
             j = nbase_gamma(ii,2)
             i2 = igf(j)
             bfft(i2) =   psi_l(ii,map_z(ib),ik,1)
          end do
       else if(kimg == 2) then
          i1 = 2*igf(1) - 1
          bfft(i1)   = psi_l(1,map_z(ib),ik,1)
          bfft(i1+1) = psi_l(1,map_z(ib),ik,2)
#ifdef NEC_TUNE_SMP
!CDIR NODEP
#endif
          do ii = 2, iba(ik)
!!$             i = nbase(ii,1)
             i = nbase(ii,ik)
             i1 = 2*igf(i)-1
             bfft(i1  ) = psi_l(ii,map_z(ib),ik,1)
             bfft(i1+1) = psi_l(ii,map_z(ib),ik,2)
             j = nbase_gamma(ii,2)
             i2 = 2*igf(j)-1
             bfft(i2  ) = psi_l(ii,map_z(ib),ik,1)
             bfft(i2+1) = -psi_l(ii,map_z(ib),ik,2)
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
             bfft(i1) = psi_l(i,map_z(ib),ik,ri)   ! MPI
          end do
       end do
    end if
    if(ipri >= 2) then
       if(ik <= 2 .and. ib <= 3) then
          write(6,'(" ! bfft G-space ik = ",i3," ib = ",i3," <<m_ES_WF_in_Rspace>>")') ik, ib
          write(6,'(8f8.4)') (bfft(i),i=1,120)
       end if
    end if
    call m_FFT_WF(ELECTRON,nfout,bfft,INVERSE,ON)
    if(ipri >= 2) then
       if(ik <=  2 .and. ib <= 3) then
          write(6,'(" ! bfft R-space ik = ",i3," ib = ",i3," <<m_ES_WF_in_Rspace>>")') ik, ib
          write(6,'(8f8.4)') (bfft(i),i=1,120)
       end if
    end if

#ifdef SAVE_FFT_TIMES
    if(sw_save_fft == ON) then
       Phifftr_l(:,map_z(ib),ik) = bfft(:)
       status_saved_phifftr(map_z(ib),ik) = STORED_AND_NEW
    end if
#endif
#ifndef NEC_TUNE_SMP
    call tstatc0_end(id_sname)
#endif
#ifdef SAVE_FFT_TIMES
    end if
#endif

  end subroutine m_ES_WF_in_Rspace1

  subroutine m_ES_WF_in_Rspace0(ik,ib,bfft)
    integer, intent(in)                           :: ik, ib
    real(kind=DP), intent(inout), dimension(nfft) :: bfft
    call m_ES_WF_in_Rspace1(ista_k,iend_k,ik,ib,zaj_l,bfft)
  end subroutine m_ES_WF_in_Rspace0

! -----  ascat starts modifying  -----
#ifdef __EDA__
  subroutine m_ES_WF_kin_in_Rspace(ik,ib,bfft_kin)
    use m_Crystal_Structure, only : rltv

    integer, intent(in)                           :: ik, ib
    integer                                       :: i, i1, ri, nb
    integer :: ii, j, i2
    real(kind=DP), intent(inout), dimension(nfft) :: bfft_kin
    real(kind=DP)                                 :: ga, gb, gc, ekin

    bfft_kin = 0.d0
    call getttr(rltv,ttr)

    if (k_symmetry(ik) == GAMMA) then
       if(kimg == 1) then
          do ii = 2, iba(ik)
             i = nbase(ii,ik)
             ga = vkxyz(ik,1,BUCS) + ngabc(i,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(i,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(i,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i1 = igf(i)
             bfft_kin(i1) = ekin *zaj_l(ii,map_z(ib),ik,1)

             j = nbase_gamma(ii,2)
             ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i2 = igf(j)
             bfft_kin(i2) = ekin *zaj_l(ii,map_z(ib),ik,1)
          end do

       else if ( kimg == 2 ) then
          do ii = 2, iba(ik)
             i = nbase(ii,ik)

             ga = vkxyz(ik,1,BUCS) + ngabc(i,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(i,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(i,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i1 = 2 *igf(i) -1
             bfft_kin(i1  ) = ekin *zaj_l(ii,map_z(ib),ik,1)
             bfft_kin(i1+1) = ekin *zaj_l(ii,map_z(ib),ik,2)

             j = nbase_gamma(ii,2)
             ga = vkxyz(ik,1,BUCS) + ngabc(j,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(j,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(j,3)
             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
               & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0

             i2 = 2*igf(j)-1
             bfft_kin(i2  ) = ekin *zaj_l(ii,map_z(ib),ik,1)
             bfft_kin(i2+1) = -ekin *zaj_l(ii,map_z(ib),ik,2)
          end do
       endif

    else
       do ri = 1, kimg
          do i = 1, iba(ik)
             nb = nbase(i,ik)
             ga = vkxyz(ik,1,BUCS) + ngabc(nb,1)
             gb = vkxyz(ik,2,BUCS) + ngabc(nb,2)
             gc = vkxyz(ik,3,BUCS) + ngabc(nb,3)
!       ekin = (ga*ga + gb*gb + gc*gc)*0.5d0

             ekin = ( ttr(1)*ga*ga + ttr(2)*gb*gb + ttr(3)*gc*gc &
                  & +   ttr(4)*ga*gb + ttr(5)*gb*gc + ttr(6)*gc*ga)*0.5d0
             
             i1 = kimg*igf(nbase(i,ik)) + (ri - kimg)
             bfft_kin(i1) = ekin *zaj_l(i,map_z(ib),ik,ri)   ! MPI
          end do
       end do
    endif

    call m_FFT_WF(ELECTRON,nfout,bfft_kin,INVERSE,ON)

  end subroutine m_ES_WF_kin_in_Rspace
#endif
! -----  ascat ceases modifying  -----

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

  subroutine m_ES_wd_zaj_small_portion0(str,nc)
    integer, intent(in) :: nc
    character(len=nc), intent(in) :: str

    integer :: ik
    do ik = ista_k, iend_k, af+1                            ! MPI
       if(ipri>=1) call m_ES_wd_zaj_small_portion(nfout,ik,str,nc)
    end do
  end subroutine m_ES_wd_zaj_small_portion0

  subroutine m_ES_wd_eko(nfout,mode)
    integer, intent(in) ::                  nfout
    integer, intent(in) ::                  mode
    real(kind=DP), pointer, dimension(:) :: eko, eko_t ! d(neg) MPI

    integer :: ik, ib, ikp
    integer :: j, k
    real(kind=DP) :: c1, vkxyz_wk(3)

! ====================================== added by K. Tagami =============== 11.0
    integer :: ikskip
! ========================================================================= 11.0

    allocate(eko(neg)); allocate(eko_t(neg))           ! MPI
    if(printable) write(nfout,*) '=== energy_eigen_values ==='

    ikp = 0
! =========================== modified by K. Tagami =================== 12.0Exp
!    if(mode == EK) ikp = nk_in_the_process - 1
    if (mode == EK) then
       if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
          ikp = nk_in_the_process - 1
       endif
    endif
! ===================================================================== 12.0Exp

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

      if ( sw_band_unfolding == ON .and. mode == EK ) then
         vkxyz_wk(1:3) = vkxyz_refcell(ikp+ik,1:3,BUCS)
      else
         vkxyz_wk(1:3) = vkxyz(ik,1:3,BUCS)
      endif

       if(map_k(ik) /= myrank_k) cycle                 ! MPI
       eko_t = 0                                       ! MPI
       do ib = 1, neg                                  ! MPI
          if(map_e(ib) == myrank_e) eko_t(ib) = eko_l(map_z(ib),ik) ! MPI
       end do                                          ! MPI
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
       if(mode == EK .or. ik == 1) then
          if(printable) write(nfout,'(" ik = ",i7," ( ",3f10.6," )",/99(4f18.10,/))')&
               &ik+ikp,vkxyz_wk(1:3), (eko(neordr(ib,ik)),ib=1,neg-num_extra_bands)
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
  end subroutine m_ES_wd_eko

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
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
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
       call mpi_allreduce(eko_t,eko,neg,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr) ! MPI
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
    real(kind=DP),intent(in),dimension(np_e,nlmta,ik:ik):: bpr_l,bpi_l !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q

    ib = map_z(ibo)
    dz = 0.d0
    do ia = 1, nac
       p = nlmta1(ia); q = nlmta2(ia)
       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr_l(ib,q,ik) + bpi_l(ib,p,ik)*bpi_l(ib,q,ik))
    end do
  end subroutine m_ES_sum_of_LocalPart


!!$  subroutine m_ES_sum_of_LocalPart(ik,ibo,bpr_l,bpi_l,dz)
!!$    integer, intent(in)                                 :: ik,ibo
!!$    real(kind=DP),intent(in),dimension(np_e,nlmta,ik:ik):: bpr_l,bpi_l !MPI
!!$    real(kind=DP),intent(out)                           :: dz
!!$
!!$    integer  :: ib, ia, p, q
!!$
!!$    ib = map_z(ibo)
!!$    dz = 0.d0
!!$    do ia = 1, nac
!!$       p = nlmta1(ia); q = nlmta2(ia)
!!$       dz = dz + fqwei(ia)*(bpr_l(ib,p,ik)*bpr_l(ib,q,ik) + bpi_l(ib,p,ik)*bpi_l(ib,q,ik))
!!$    end do
!!$  end subroutine m_ES_sum_of_LocalPart

! =============================== added by K. Tagami ================= 11.0
  subroutine m_ES_sum_of_LocalPart_noncl( ik,ibo,kst, ken, bpr_l,bpi_l,dz)
    integer, intent(in)  :: ik,ibo
    integer, intent(in)  :: kst, ken
    real(kind=DP), intent(in) :: bpr_l(np_e,nlmta,kst:ken)
    real(kind=DP), intent(in) :: bpi_l(np_e,nlmta,kst:ken)
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q
    integer :: is1, is2, is_tmp, k1, k2
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
                  &   *( bpr_l(ib,p,k1) *bpr_l(ib,q,k2) &
                  &    + bpi_l(ib,p,k1) *bpi_l(ib,q,k2))
             c2 =-aimag(fqwei_noncl(ia,is_tmp) ) &
                  &   *( bpr_l(ib,p,k1) *bpi_l(ib,q,k2) &
                  &     -bpi_l(ib,p,k1) *bpr_l(ib,q,k2))

             dz = dz + c1 + c2

          End do
       End Do
    end do
  end subroutine m_ES_sum_of_LocalPart_noncl
! ================================================================ 11.0


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
    real(kind=DP),intent(in),dimension(np_e,nlmta,ik:ik):: bpr_l,bpi_l !MPI
    real(kind=DP),intent(in),dimension(np_e,nlmta,ista_k:iend_k):: bpr1_l,bpi1_l !MPI
    real(kind=DP),intent(out)                           :: dz

    integer  :: ib, ia, p, q

    ib = map_z(ibo)
    dz = 0.d0
    do ia = 1, nac
       p = nlmta1(ia); q = nlmta2(ia)
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
       do ib = ista_e, iend_e, istep_e
          if(nrvf_ordr(ib,ik) > neg - num_extra_bands) cycle
          tmp = eko_l(map_z(ib),ik) - eko1_l(map_z(ib),ik)
          evdff(1) = evdff(1) + tmp*tmp
          evdff(2) = evdff(2) + dabs(tmp)
          ekosum = ekosum + eko_l(map_z(ib),ik)
          if(dabs(tmp).gt.DELTAevdff) then
             evdff(3) = evdff(3) + dsqrt(dabs(eko1_l(map_z(ib),ik)**2 - eko_l(map_z(ib),ik)**2))
          end if
       end do
    end do

    if(npes > 1) then
       call mpi_allreduce(mpi_in_place,evdff,3 &
         & ,mpi_double_precision,mpi_sum,mpi_e_world(myrank_e),ierr)  ! MPI
       call mpi_allreduce(evdff(1),evdffr(1),1 &
         & ,mpi_double_precision,mpi_max,mpi_k_world(myrank_k),ierr)  ! MPI
       call mpi_allreduce(evdff(2),evdffr(2),1 &
         & ,mpi_double_precision,mpi_max,mpi_k_world(myrank_k),ierr)  ! MPI
       call mpi_allreduce(evdff(3),evdffr(3),1 &
         & ,mpi_double_precision,mpi_max,mpi_k_world(myrank_k),ierr)  ! MPI
       call mpi_allreduce(MPI_IN_PLACE,ekosum,1,mpi_double_precision,mpi_sum,mpi_k_world(myrank_k),ierr)
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
!       call mpi_allreduce(evdff2,evdffr2,3*kv3 &
!            & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)  ! MPI
       call mpi_allreduce(evdff2,evdffr2,3*kv3 &
            & ,mpi_double_precision,mpi_sum,mpi_spin_group,ierr)  ! MPI
       evdff2 = evdffr2
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
!       call mpi_allreduce(evdff,evdffr,3 &
!         & ,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)  ! MPI
       call mpi_allreduce(evdff,evdffr,3 &
         & ,mpi_double_precision,mpi_sum,mpi_spin_group,ierr)  ! MPI
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
         & ,mpi_sum,mpi_k_world(myrank_k),ierr)
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
       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, MPI_CommGroup,ierr)
!       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
!            &              ,mpi_sum, mpi_spin_group,ierr)
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
!       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
!            &              ,mpi_sum, MPI_CommGroup,ierr)
       call mpi_allreduce(eko_t,eko_t2,neg*kv3,mpi_double_precision &
            &              ,mpi_sum, mpi_spin_group,ierr)
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
         &            ,mpi_sum,MPI_CommGroup,ierr)
!    call mpi_allreduce(porb_mpi,porb,nlmta_phi*2,mpi_double_precision &
!         &            ,mpi_sum,mpi_spin_group,ierr)

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
                   call phase_error_with_msg(nfout, 'Not supported : Gamma symmetry in noncollinear system.', &
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
       call phase_error_with_msg(nfout, 'kt : Under construction porb ',__LINE__,__FILE__)
    endif
! -----

    call m_ES_DensMat_To_MagMom_porb( nlmta_phi, porb_ssrep, porb )

    if ( npes >=2 ) then
       allocate( porb_mpi(nlmta_phi,ndim_magmom ) )
       porb_mpi = 0.0d0
!       call mpi_allreduce( porb, porb_mpi, nlmta_phi*ndim_magmom, &
!            &              mpi_double_precision, mpi_sum,MPI_CommGroup, ierr )
       call mpi_allreduce( porb, porb_mpi, nlmta_phi*ndim_magmom, &
            &              mpi_double_precision, mpi_sum,mpi_spin_group, ierr )
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
!    call mpi_allreduce(dm_mpi,dm,max2lp**2*max_projs*natm*nspin,mpi_double_precision &
!                    & ,mpi_sum,MPI_CommGroup,ierr)
    call mpi_allreduce(dm_mpi,dm,max2lp**2*max_projs*natm*nspin,mpi_double_precision &
                    & ,mpi_sum,mpi_spin_group,ierr)

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

!       call mpi_allreduce( dm, dm_mpi, max2lp**2*max_projs*natm*ndim_magmom, &
!         &                 mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       call mpi_allreduce( dm, dm_mpi, max2lp**2*max_projs*natm*ndim_magmom, &
         &                 mpi_double_precision, mpi_sum, mpi_spin_group, ierr )
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
     call m_CntrlP_set_neg(newneg,nfout)
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
          do ig=1,iba_prev(ik)
!             zajtmp(ig,:) = zaj_l_prev(ig,map_z(ib),ik,:)
             zajtmp(ig,:) = zaj_l_prev(ig,ib,ik,:)
          enddo
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
            do ig=1,iba(ik)
!              zaj_l(ig,map_z(ib),ik,ri) = bfft(kimg*igf(nbase(ig,ik))+(ri-kimg))
              zaj_l(ig,ib,ik,ri) = bfft(kimg*igf(nbase(ig,ik))+(ri-kimg))
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

  subroutine m_ES_alloc_zaj_l_prev(kg1)
    integer, intent(in) :: kg1
    if(.not.allocated(zaj_l_prev)) allocate(zaj_l_prev(kg1,np_e,ista_k:iend_k,kimg))
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
      allocate(zaj_l_prev(kg1,np_e,ista_k:iend_k,kimg))
      zaj_l_prev = zaj_l
    endif

    if(allocated(vlhxcQ)) deallocate(vlhxcQ)
    if(allocated(vlhxc_l)) deallocate(vlhxc_l)
    if(allocated(vnlph_l)) deallocate(vnlph_l)

    if(allocated(zaj_l)) deallocate(zaj_l)
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
    if(allocated(neordr)) deallocate(neordr)
    if(allocated(nrvf_ordr)) deallocate(nrvf_ordr)
    if(allocated(eko_ek)) deallocate(eko_ek)
    if(allocated(eko1_l)) deallocate(eko1_l)
    if(allocated(evdff)) deallocate(evdff)
    if(allocated(evdffr)) deallocate(evdffr)
    if(allocated(compr_l)) deallocate(compr_l)
    if(allocated(compi_l)) deallocate(compi_l)
    if(allocated(fsr_add_l)) deallocate(fsr_add_l)
    if(allocated(fsi_add_l)) deallocate(fsi_add_l)

  end subroutine m_ES_dealloc

  subroutine m_ES_add_it_to_vnlph(ik,ib,v)
    integer, intent(in) :: ik,ib
    real(kind=DP), intent(in) :: v(kg1,kimg)

    integer :: ig

    if(kimg==1) then
       do ig=1,iba(ik)
          vnlph_l(ig,ib,1) = vnlph_l(ig,ib,1) + v(ig,1)
       end do
    else
       do ig=1,iba(ik)
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
       if ( .not.allocated(fsr_tmp) ) allocate( fsr_tmp(np_e,nlmta,ista_k:iend_k))
       fsr_tmp = fsr_l
       if ( .not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2) ) then
          if ( .not.allocated(fsi_tmp) ) allocate(fsi_tmp(np_e,nlmta,ista_k:iend_k));
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

end module m_Electronic_Structure
