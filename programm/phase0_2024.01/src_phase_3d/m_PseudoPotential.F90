!!$#define _DEBUG_WRITE_
!!$#define _PAW_DEBUG_WRITE_
#define _PAW_CONTINUE_DATA_PREVIOUS_BEFORE_201403_STYLE_
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MODULE: m_PseudoPotential
!
!  AUTHOR(S): T. Yamasaki,  Y. Morikawa, M. Okamoto, August/20/2003
!      Further modification by T. Yamasaki   Aug/31/2007
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!  patch 10.1 by K. Tagami @adv    2011/06/18
!
!  patch 10.1 :  addition of vec_q_plusG etc for TD-DFT
!=======================================================================
!!
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
#ifdef __TIMER_IODO__
#   define __TIMER_IODO_START(a)   call timer_sta(a)
#   define __TIMER_IODO_STOP(a)    call timer_end(a)
#else
#   define __TIMER_IODO_START(a)
#   define __TIMER_IODO_STOP(a)
#endif
#ifdef __TIMER_IOCOMM__
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)   call timer_barrier(str) ;   call timer_sta(a)
#   define __TIMER_IOCOMM_START(a)       call timer_sta(a)
#   define __TIMER_IOCOMM_STOP(a)        call timer_end(a)
#else
#   define __TIMER_IOCOMM_START_w_BARRIER(str,a)
#   define __TIMER_IOCOMM_START(a)
#   define __TIMER_IOCOMM_STOP(a)
#endif

module m_PseudoPotential
!     (m_PP)
!  $Id: m_PseudoPotential.F90 633 2020-12-01 05:11:03Z jkoga $
!
!  The original subroutine name was "pspot", which had been coded by
!  Y. Morikawa (JRCAT-NAIR) in 1993 or ealier.
!!$c   #1) betar --> betar_l by T.Yamasaki on 3rd Mar. 1995
!!$c   #2) phipdos --> added into the output.
!!$c      for pdos calculation 25th July 1995 Y.M
!!$c   #3) LMTO pseudopotential is introduced.
!!$c       26th Dec. 1995 Y.Morikawa
!!$c   #4) f-non-locality included. 22nd Jan. 1996 Y.Morikawa
!
!  Translated into this module coded using fortran90+mpi by T. Yamasaki in 1999.
!
  use m_Crystal_Structure,  only   : univol,nopr,tau,m_CS_op_in_PUCD,op
  use m_Ionic_System,       only   : ntyp,ivan,iatomn,iatom,natm, iwei, ityp &
       &                           , if_pdos, iproj_group, pos &
       &                           , m_IS_set_ival
  use m_Timing,             only   : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only   : len_xctype, xctype, ipri, iprivloc, istress, icond &
       &                           , ipripp, sw_positron, m_CtrlP_set_xctype &
       &                           , norbital, l_orb, t_orb, rc_orb, k_orb &
       &                           , sw_use_add_proj, sw_berry_phase &
       &                           , sw_orb_popu, initial_chg &
       &                           , sw_hubbard, num_projectors &
       &                           , proj_attribute, projector_type &
       &                           , proj_group, num_proj_elems &
       &                           , intzaj, kimg, sw_wannier, corecharge_cntnbin, sw_epsilon_ele &
       &                           , epsilon_ele, ppprinted &
! === DEBUG by tkato 2011/09/22 ================================================
!      &                           , nspin, af, printable, PAW_switch
       &                           , nspin, af, printable, PAW_switch, sw_fef &
       &                           , omega_exx, sw_rspace_hyb &
! ==============================================================================
       &                           , sw_exclude_orb_tau_neq_1 &
       &                           , sw_exclude_orb_unoccupied &
       &                           , sw_remove_pcc_from_pawpot, sw_wannier90 &
       &                           , sw_optimize_lattice,imdalg, driver   &
       &                           , sw_stress_correction &
#ifndef ENABLE_ESM_PACK
       &                           , m_CtrlP_explicit_cmix
#else
       &                           , m_CtrlP_explicit_cmix &
       &                           , sw_esm
#endif
  use m_Const_Parameters,   only   : DP,PAI,PAI4,NO,YES,SKIP,EXECUT,GGA,LDA,ON,OFF,DELTA &
       &                           , INITIAL, CONTINUATION, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                           , BULK, DEFECT, from_PSEUDOPOTENTIAL_FILE &
       &                           , ATOMIC_ORBITAL, SPHERICAL_HARMONICS &
       &                           , by_pseudo_atomic_orbitals &
       &                           , BUCS, pp_GNCPP1, pp_GNCPP2, pp_GNCPP2_with_AE_WF, pp_PAW,pp_CIAOPP &
       &                           , DELTA10, COORDINATE_CONTINUATION, PT_CONTROL, P_CONTROL &
       &                           , DRIVER_URAMP, DRIVER_SC_DFT, InvHyperFineConst
  use m_Parallelization,    only   : MPI_CommGroup,npes,mype,ierr,ista_kngp,iend_kngp &
       &                           , ista_fs,     iend_fs,     np_fs,     nel_fs &
       &                           , ista_fs_atm, iend_fs_atm, np_fs_atm, nel_fs_atm  &
       &                           , myrank_g, myrank_e, nrank_g, myrank_chg
! ================================== added by K. Tagami ================ 5.0
  use m_Control_Parameters,  only : sw_orbital_cut_smooth, sw_mix_charge_hardpart
! ====================================================================== 5.0
  use m_Control_Parameters,  only : sw_hybrid_functional,sw_screened_exchange


! ====================================== Added by K. Tagami === 10.1
  use m_Control_Parameters,     only : sw_LinearResponse
! ============================================================= 10.1

! ==================================== added by K. Tagami ================ 11.0&13.0U2
  use m_Control_Parameters,   only : noncol, ndim_chgpot, ndim_magmom, &
       &                             sw_potential_mixing, sw_modified_TFW_functional, &
       &                             initial_occmat, long_range_pot_type, &
       &                             sw_add_vlocal_gzero
  use m_Const_Parameters,     only : CMPLDP, CRYSTAL_FIELD_APPROX
! ======================================================================= 11.0&13.0U2

! ==================================== added by K. Tagami ================ 11.0
  use m_Control_Parameters,   only : SpinOrbit_Mode, SpinOrbit_MassCorrection, &
       &                             sw_use_rphi_Hsoc_rphi, &
       &                             sw_use_ival_for_paw_ps_soc, &
       &                             sw_corelevel_spectrum, sw_calc_core_energy
  use m_Const_Parameters,     only : ZeffApprox, ReadFromPP, ByPawPot
! ========================================================================= 11.0

  use m_Control_Parameters,   only : sw_charge_predictor, cntn_bin_paw_format, &
       &                             cntn_bin_paw_format_is_set
  use m_ErrorMessages,        only : INVALID_ATOMIC_NUMBER
! === DEBUG by tkato 2014/03/04 ================================================
  use m_Const_Parameters, only : SmallestPositiveNumber
! ==============================================================================

! ==== KT_add ==== 2014/08/26
  use m_Control_Parameters,   only : sw_use_contracted_psir, orb_popu_method, &
       &                             howto_set_proj_radius, sw_add_corecharge_rspace, &
       &                             eval_corecharge_on_Gspace, &
       &                             sw_read_corecharge_extra_file
! ================ 2014/08/26

! ==== KT_add ==== 2014/09/19 & 13.1XI
  use m_Control_Parameters,   only : sw_calc_ekin_density, use_symm_ekin_density, &
       &                             use_asymm_ekin_density, sw_rspace_ekin_density, &
       &                             sw_excitation, sw_write_rotated_orbitals, &
       &                             use_metagga
! ================ 2014/09/19 & 13.1XI
! == KT_add === 13.2S
  use m_Crystal_Structure,  only : sw_spinorbit_second_variation
! ============= 13.2S
#ifdef __EDA__
  use m_Control_Parameters, only : sw_eda
#endif
  use mpi
  implicit none

  character(len=len("FILEFORMAT")),private,parameter :: tag_fileformat = "FILEFORMAT"
  character(len=len("CIAOPP")),private, parameter    :: tag_ciaopp     = "CIAOPP"
  character(len=len("GNCPP2")),private, parameter    :: tag_gncpp2     = "GNCPP2"
  character(len=len("POTENTIALTYPE")),private, parameter :: tag_potentialtype = "POTENTIALTYPE"
  character(len=len("PAW")),private, parameter       :: tag_paw        = "PAW"
  character(len=len("UltraSoft")),private,parameter  :: tag_ultrasoft  = "ULTRASOFT"

! =============================== added by K. Tagami ==================== 11.0
  character(len=len("SpinOrbit")),private, parameter  :: &
	&                                tag_pot_has_spinorbit = "SpinOrbit"
  character(len=len("ON")),       private, parameter  :: tag_ON = "ON"
! ======================================================================= 11.0

!!$  integer, dimension(ntyp)      :: ivan  ! {0=old | 1=new} format type
  integer, parameter :: MAXIMUM_MPI_SIZE =   100000000
  integer, parameter :: MAXIMUM_BCAST_SIZE = 20000000 !   200000000

  real(kind=DP),allocatable,dimension(:)   :: ival  ! d(ntyp) #valence electrons
  integer,allocatable,dimension(:,:) :: ivanl ! d(nloc,ntyp) {0=Norm-conserv.|1= non-c.}

! ========= KT_mod ============ 13.0U2
!  integer,private,allocatable,dimension(:,:) :: itau  ! d(nloc,ntyp) #reference energies
  integer, allocatable,dimension(:,:) :: itau  ! d(nloc,ntyp) #reference energies
! ============================= 13.0U2

  integer                                    :: ntau  ! maximum of itau
  real(kind=DP),private,allocatable,dimension(:,:,:) :: eps ! d(nloc,ntau,ntyp)
  integer,private,allocatable,dimension(:,:,:) :: irank ! d(nloc,ntau,ntyp)
  integer,allocatable,dimension(:)   :: nmesh ! d(ntyp) log-mesh size
  integer                            :: mmesh ! = maxval(nmesh)
  integer,allocatable,dimension(:)   :: iloc  ! d(ntyp) local orbital
  integer,allocatable,dimension(:)   :: lpsmax! d(ntyp) # valence orbitals
  integer                            :: nloc  ! maxval(lpsmax)
  integer,allocatable,dimension(:)   :: itpcc ! d(ntyp) {1=particla c.c. on | 0= off}
  integer,allocatable,dimension(:)   :: ilmt  ! d(ntyp) #comb. of l, m, and t
  integer                            :: nlmt  ! number of useful combinations of l, m, and tau.
  integer                            :: nlmtt ! summation of nlmt's for species
  integer                            :: nlmta ! number of the mapping Functions, lmta.
  integer,allocatable,dimension(:,:) :: lmta  ! d(nlmt,natm) mapping F. to 'fsr_l,fsi_l'.
  integer,allocatable,dimension(:,:) :: lmtt  ! d(nlmt,ntyp) mapping F. to 'snl' array.
  integer,allocatable,dimension(:,:) :: ltp   ! d(nlmt,ntyp) mapping F. ilmt to l.
  integer,allocatable,dimension(:,:) :: mtp   ! d(nlmt,ntyp) same      ilmt to m.
  integer,allocatable,dimension(:,:) :: taup  ! d(nlmt,ntyp) same      ilmt to ref. enrg.

  integer,allocatable,dimension(:)   :: ilmt_phi  ! d(ntyp) #comb. of l, m, and t
  integer                            :: nlmt_phi  ! number of useful combinations of l, m, and tau.
  integer                            :: nlmtt_phi ! summation of nlmt's for species
  integer                            :: nlmta_phi ! number of the mapping Functions, lmta.
  integer,allocatable,dimension(:,:) :: lmta_phi  ! d(nlmt_phi,natm) mapping F. to 'fsr_l,fsi_l'.
  integer,allocatable,dimension(:,:) :: lmtt_phi  ! d(nlmt_phi,ntyp) mapping F. to 'snl' array.
  integer,allocatable,dimension(:,:) :: ltp_phi   ! d(nlmt_phi,ntyp) mapping F. ilmt to l.
  integer,allocatable,dimension(:,:) :: mtp_phi   ! d(nlmt_phi,ntyp) same      ilmt to m.
  integer,allocatable,dimension(:,:) :: taup_phi  ! d(nlmt_phi,ntyp) same      ilmt to ref. enrg.
  integer,allocatable,dimension(:,:) :: iproj_phi   ! d(nlmt_phi,ntyp) mapping F. ilmt to iproj.
  real(kind=DP),allocatable,dimension(:,:,:) :: k_phi ! d(nloc,ntau,ntyp)
  real(kind=DP),allocatable,dimension(:,:,:) :: rc_phi ! d(nloc,ntau,ntyp)

  integer,allocatable,dimension(:)   :: ilmt_add ! d(ntyp) #comb. of l, and m
  integer                            :: nlmt_add  ! number of useful combinations of l, and m
  integer,save                       :: nlmtt_add = 0 ! summation of nlmt's for species
  integer                            :: nlmta_add ! number of the mapping Functions, lmta.
  integer,allocatable,dimension(:,:) :: lmta_add  ! d(nlmt_add,natm) mapping F. to 'fsr_add_l,fsi_add_l'.
  integer,allocatable,dimension(:,:) :: lmtt_add  ! d(nlmt_add,ntyp) mapping F. to 'snl_add' array.
  integer,allocatable,dimension(:,:) :: ltp_add   ! d(nlmt_add,ntyp) mapping F. ilmt to l.
  integer,allocatable,dimension(:,:) :: mtp_add   ! d(nlmt_add,ntyp) same      ilmt to m.
  real(kind=DP),allocatable,dimension(:)   :: kappa  ! d(ntyp) window function param.
  real(kind=DP),allocatable,dimension(:)   :: rcproj ! d(ntyp) window function param.

  integer,allocatable,dimension(:)   :: ilmt_pao  ! d(ntyp) #comb. of l, m, and t
  integer                            :: nlmt_pao  ! number of useful combinations of l, m, and tau.
  integer,save                       :: nlmtt_pao = 0 ! summation of nlmt's for species
  integer,save                       :: nlmta_pao = 0 ! number of the mapping Functions, lmta.
  integer,allocatable,dimension(:,:) :: lmta_pao  ! d(nlmt_pao,natm) mapping F. to 'fsr_l,fsi_l'.
  integer,allocatable,dimension(:,:) :: lmtt_pao  ! d(nlmt_pao,ntyp) mapping F. to 'snl' array.
  integer,allocatable,dimension(:,:) :: ltp_pao   ! d(nlmt_pao,ntyp) mapping F. ilmt to l.
  integer,allocatable,dimension(:,:) :: mtp_pao   ! d(nlmt_pao,ntyp) same      ilmt to m.
  integer,allocatable,dimension(:,:) :: taup_pao  ! d(nlmt_pao,ntyp) same      ilmt to ref. enrg.
  integer,allocatable,dimension(:) :: ibpao !d(nlmta_pao)

  integer :: num_wfc_pao = 0

  real(kind=DP)                 :: etot1 ! total energy correction from PP
#ifdef __EDA__
! -----  ascat starts modifying  -----
  real(kind=DP),allocatable,dimension(:) :: PP_per_atom
  real(kind=DP),allocatable,dimension(:) :: epc_per_atom
! -----  ascat ceases modifying  -----
#endif
  real(kind=DP)                 :: epc = 0.d0  ! partial core (no overlap) charge xc energy

! ==== KT_add ==== 13.0S
  real(kind=DP)    ::              epc_paw = 0.d0   ! pcc energy when paw
! ================ 13.0S

  real(kind=DP),private         :: schgpc! total partial charge

  real(kind=DP),allocatable,dimension(:,:,:)     :: q    ! d(nlmt,nlmt,ntyp) deficit charge
  integer,allocatable,dimension(:,:,:)           :: if_q_isnotzero !d(nlmt,nlmt,ntyp)
  real(kind=DP),allocatable,dimension(:,:,:)     :: dion ! d(nlmt,nlmt,ntyp)
  integer, allocatable,dimension(:,:,:,:)        :: isph ! d(nlmt,nlmt,6,ntyp)
  integer, allocatable,dimension(:,:,:)          :: il2p ! d(nlmt,nlmt,ntyp)
  real(kind=DP),allocatable,dimension(:,:,:,:)   :: dl2p ! d(nlmt,nlmt,6,ntyp)

  integer                                        :: modnrm, nac
  real(kind=DP), allocatable,dimension(:)        :: fqwei !d(nlmt*ntau*natm)
  integer, allocatable, dimension(:)             :: nlmta1, nlmta2, nac2ia !d(nlmt*ntau*natm)
  integer                                        :: nac_p ! =: nlmt*ntau*nel_fs_atm(myrank_e)
  real(kind=DP), allocatable,dimension(:)        :: fqwei_p !d(nac_p)
  integer, allocatable, dimension(:)             :: nlmta1_p, nlmta2_p !d(nac_p)

  integer                                        :: m_non0_lmtxlmt = -1 !=maxval(n_non0_lmtxlmt)
  integer,public,allocatable,dimension(:,:,:)    :: index_lmt1_lmt2 !d(m_non0_lmtxlmt,ntyp,2)
  integer,private,allocatable,dimension(:,:)     :: w_non0_lmtxlmt  !d(m_non0_lmtxlmt,ntyp)
  integer,public,allocatable,dimension(:)        :: n_non0_lmtxlmt  !d(ntyp)
  real(kind=DP),private,allocatable,dimension(:,:):: q_indp !d(m_non0_lmtxlmt,ntyp)
  real(kind=DP),private,allocatable,dimension(:,:):: dion_indp !d(m_non0_lmtxlmt,ntyp)

  real(kind=DP),allocatable,dimension(:,:) :: psc_l,psc_diff_l    !d(ista_kngp:iend_kngp,ntyp)
  real(kind=DP),allocatable,dimension(:,:) :: qitg_l, qitg_diff_l !d(ista_kngp:iend_kngp,nqitg)
  real(kind=DP),allocatable,dimension(:,:) :: rhpcg_l,rhpcg_diff_l!d(ista_kngp:iend_kngp,ntpcc)

  real(kind=DP),allocatable,dimension(:,:) :: rhvg_l !d(ista_kngp:iend_kngp,ntyp)
  !                                            valence charge density in G-space
  real(kind=DP),allocatable,dimension(:,:) :: rhcg_l,rhceg_l !d(ista_kngp:iend_kngp,ntyp)
  real(kind=DP),allocatable,dimension(:,:) :: rhchg_l        !d(ista_kngp:iend_kngp,ntyp)
  !                                            core charge density in G-space
  real(DP), private :: total_core_charge

  integer, parameter ::               len_ps_xctype = 7
  character(len=len_ps_xctype)  ::    ps_xctype
  integer, allocatable, dimension(:,:,:,:,:,:) :: iqitg  !d(nloc,ntau,nloc,ntau,nloc*2,ntyp)
  integer,save                                 :: ntpcc = 0
  integer                                      :: nqitg
  integer, allocatable, dimension(:)           :: nqitg_sp !d(ntyp)
  integer, private                             :: mmt
  !         mmt is used only in the subroutines of <<m_PP_vanderbilt_type>>

  ! for Berry Phase calc.
  logical, allocatable, dimension(:) :: ae_wavefunctions_are_detected !d(ntyp)
  real(kind=DP) :: dk_BP(1)
  real(kind=DP),allocatable,dimension(:,:) :: qitg_BP !d(1,nqitg)

  ! for Wannier function
  integer :: nwght_wan
  real(kind=DP), allocatable, dimension(:) :: dk_wan ! d(nwght_wan)
  real(kind=DP),allocatable,dimension(:,:) :: qitg_wan !d(nqitg,nwght_wan)

! ============= Added by K. Tagami ======for LinearResponse======= 10.1
  integer :: nmax_q_plus_G_LR, nmax_q_LR
  real(kind=DP), allocatable :: qitg_LR(:,:,:)
  real(kind=DP), allocatable :: norm_q_plus_G_LR(:,:)
  real(kind=DP), allocatable :: vec_q_plus_G_LR(:,:,:)
! ================================================================ 10.1
  ! for Finite Electric Field
  integer :: numef_fef
  real(kind=DP), allocatable, dimension(:) :: dk_fef ! d(numef_fef)
  real(kind=DP),allocatable,dimension(:,:) :: qitg_fef !d(nqitg,numef_fef)

!              ------ temporary arrays ---
  real(DP),private,allocatable,dimension(:)     :: h,chgpc ! d(ntyp)
  real(DP),        allocatable,dimension(:)     :: xh,rmax !(ntyp)
  real(DP),private,allocatable,dimension(:,:)   :: alp,cc  ! d(2,ntyp)
  real(DP),        allocatable,dimension(:,:,:,:):: betar   ! d(mmesh,nloc,ntau,ntyp)
  real(DP),        allocatable,dimension(:,:)   :: betar_add ! d(mmesh,ntyp)
  real(DP),        allocatable,dimension(:,:,:,:):: phirt   ! d(mmesh,nloc,ntau,ntyp)
  real(DP),        allocatable,dimension(:,:,:,:):: paor   ! d(mmesh,nloc,ntau,ntyp)
  real(DP),        allocatable,dimension(:,:)   :: q_phi    ! d(nlmt_phi,ntyp)
  real(DP),        allocatable,dimension(:,:)   :: porb     ! d(nlmta_phi,nspin)
  real(DP),        allocatable,dimension(:)     :: qorb     ! d(nlmta_phi)
  real(DP),private,allocatable,dimension(:)     :: vvv,vlocr,vlocr2,rhvr,vxc,exc ! d(mmesh)
  real(DP),private,allocatable,dimension(:)     :: qrsps,wkx,wky,wkz,rhpcr,rhcr ! d(mmesh)
  !               rhpcr: partial core charge,  rhcr: core charge density in r-space
  real(DP),private,allocatable,dimension(:)     :: qrs ! d(mmesh)
  real(DP),        allocatable,dimension(:)     :: wos,radr
  real(DP),private,allocatable,dimension(:)     :: wkz1,wkz2               ! d(mmesh)
  real(DP),private,allocatable,dimension(:,:,:) :: phir_ae ! d(mmesh,nloc,ntau)
  real(DP),private,allocatable,dimension(:,:,:) :: phir, chir  ! d(mmesh,nloc,ntau)
  real(DP),private,allocatable,dimension(:,:,:) :: qij,qvij    ! d(ntau,ntau,nloc)
  integer, private,parameter                    ::  kord = 20
  real(DP),private,allocatable,dimension(:)     :: copsc       ! d(kord+1)
  real(DP),private,allocatable,dimension(:,:)   :: s, sinv     ! d(ntau,ntau)

  integer, private,parameter                    :: mddiff = 8
  real(DP),private,allocatable, dimension(:,:)  :: fdiff, grdnts,coeff1,coeff2
                      ! d(mmesh,1+mddiff),d(mmesh,3),d(1+mddiff,mddiff)*2
  real(DP),private,allocatable, dimension(:,:)  :: wwk
  real(DP),private,allocatable, dimension(:)    :: rho  ! d(mmesh)

  integer, private,parameter      :: lcmax = 4
  integer, private,parameter      :: kloc_t = 5
  integer, private,parameter      :: ktau_t = 3
  integer, private,parameter      :: klmt_t = 48

  real(DP),private                :: eps_chg  = 1.d-25
  real(DP),private                :: eps_grad = 1.d-40

  integer, private                :: iflag_wos  = OFF
  integer, private                :: iflag_radr = OFF
  integer, private                :: iflag_betar = OFF
  integer, private                :: iflag_ilmt  = OFF
  integer, private, parameter ::     len_str = 80
  character(len=len_str), private :: str

  character(len("core-charge")),private,parameter :: tag_core_charge = "core-charge"

  ! for continuation with the PDOS, Berry-phase, or Wannier calculation
  integer :: iflag_ft = ON

  ! for symmetrization of PDOS
  integer, allocatable :: irorb(:,:,:) !d(2*nloc+1,nlmta_phi,nopr)
  integer, allocatable :: nrorb(:,:) !d(nlmta_phi,nopr)
  real(kind=DP), allocatable :: crorb(:,:,:) !d(2*nloc+1,nlmta_phi,nopr)

  ! ADFT+U: Hubbard model
  real(kind=DP), allocatable :: prodphi(:,:,:) ! d(num_projectors,ntau,ntau)


!              ------ paw related variables ---
  integer,         allocatable,dimension(:) :: ipaw  ! d(ntyp)
  logical                                      :: flg_paw=.false.
  real(DP),private,allocatable,dimension(:) :: chgcr ! d(ntyp)
  real(DP),private,allocatable,dimension(:) :: rhcor ! d(mmesh)
  real(DP),allocatable,dimension(:,:)        :: rhcorpw ! d(mmesh,ntyp)
  real(DP),allocatable,dimension(:,:)        :: rhpcrpw ! d(mmesh,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: psirpw  ! d(mmesh,nloc,ntau,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: phirpw  ! d(mmesh,nloc,ntau,ntyp)
  real(kind=DP),allocatable,dimension(:,:)   :: qrspspw ! d(mmesh,nqitg)
  real(DP),private,allocatable,dimension(:,:,:) &
                                                :: psir ! d(mmesh,nloc,ntau)
! ======== KT_add =================== 13.0U2
  real(DP),allocatable,dimension(:,:,:,:)    :: psir_val  ! d(mmesh,nloc,ntau,ntyp)
! =================================== 13.0U2

! ====== KT_add ==== 2014/08/11
  real(DP),allocatable,dimension(:,:)    :: vlocr_pw  ! d(mmesh,ntyp)
! ================== 2014/08/11

!!$  integer,private,allocatable,dimension(:,:,:) &
  integer,allocatable,dimension(:,:,:) &
                                                :: wf_nrc ! d(nloc,ntau,ntyp)
  integer,allocatable,dimension(:)          :: wf_mnrc ! d(ntyp)
  integer,private                             :: mmpppp
  integer, allocatable, dimension(:,:,:,:,:,:) &
                                            :: ipppp  !d(nltpw,nltpw,nltpw,nltpw,nloc*2,ntyp)
  integer                                  :: npppp
  integer,allocatable,dimension(:)      :: iltpw  ! d(ntyp) #comb. of l and t
  integer                                  :: nltpw  ! number of useful combinations of l and tau.
  integer,allocatable,dimension(:,:)    :: lppw   ! d(nltpw,ntyp) mapping F. ilt to l.
  integer,allocatable,dimension(:,:)    :: tppw   ! d(nltpw,ntyp) same      ilt to ref. enrg.
  real(DP),allocatable,dimension(:)      :: vaeijlm_k     ! d(npppp)
  real(DP),allocatable,dimension(:)      :: vpsijlm_k     ! d(npppp)
  real(DP),allocatable,dimension(:)      :: vqijqlm_k     ! d(npppp)
  real(DP),allocatable,dimension(:)      :: vqijplpm_ks   ! d(npppp)
  real(DP),allocatable,dimension(:,:,:,:)    :: vionaeij      ! d(nloc,ntau,ntau,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: vionpsij      ! d(nloc,ntau,ntau,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: vionpsqij     ! d(nloc,ntau,ntau,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: kin_ae_psij   ! d(nloc,ntau,ntau,ntyp)
  real(DP),allocatable,dimension(:)          :: vloc_scr_ps   ! d(mmesh)
  real(DP),allocatable,dimension(:)          :: vloc_scr_ae   ! d(mmesh)
  integer                                      :: m_clmns_cijkclmk
  integer,allocatable,dimension(:,:,:)      :: n_cijkclmk    ! d(nlmt,nlmt,it)
  integer,allocatable,dimension(:,:,:,:)    :: ilmt3_cijkclmk
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,ntyp)
  integer,allocatable,dimension(:,:,:,:)    :: ilmt4_cijkclmk
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: CijkClmkVVVVijlm_k
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,ntyp)
  real(DP),allocatable,dimension(:,:,:,:)    :: CijkClmkVVVVijlm_k_ae
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,ntyp)
  integer,allocatable,dimension(:,:)          :: index_lmt2lt ! d(nlmt,ntyp)

  real(DP),allocatable,dimension(:,:,:)      :: dion_kin_ion     ! d(nlmt,nlmt,ntyp)
  real(DP),allocatable,dimension(:,:,:)      :: dion_hartree     ! d(nlmt,nlmt,natm)
  real(DP),allocatable,dimension(:,:,:)      :: dion_hartree_now ! d(nlmt,nlmt,natm)
  real(DP),allocatable,dimension(:,:,:,:)    :: dion_vxc         ! d(nlmt,nlmt,nspin,natm)
  real(DP),allocatable,dimension(:,:,:,:)    :: dion_paw         ! d(nlmt,nlmt,nspin,natm)

  real(DP),allocatable,dimension(:,:)        :: radr_paw     ! d(mmesh,ntyp)

  integer                                      :: m_clmns_cijkclmn
  integer,allocatable,dimension(:,:,:,:)    :: n_cijkclmn    ! d(nlmt,nlmt,natm,nopr)
  integer,allocatable,dimension(:,:,:,:,:)  :: ilmt3_cijkclmn
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,natm,nopr)
  integer,allocatable,dimension(:,:,:,:,:)  :: ilmt4_cijkclmn
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,natm,nopr)
  real(DP),allocatable,dimension(:,:,:,:,:)  :: CijkClmnVVVVijlm_kn
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,natm,nopr)
  real(DP),allocatable,dimension(:,:,:,:,:)  :: CijkClmnVVVVijlm_kn_ae
                                                ! d(nlmt,nlmt,m_clmns_cijkclmk,natm,nopr)

  logical:: flg_symmtry=.true.
  integer,pointer,dimension(:,:) :: ia2ia_symmtry_op
  integer,pointer,dimension(:,:) :: ia2ia_symmtry_op_inv
  real(DP),pointer,dimension(:,:,:,:) :: crotylm_paw
  integer,pointer,dimension(:,:,:,:) :: iylm_paw
  integer,pointer,dimension(:,:,:)   :: nylm_paw

  character*80,private :: str_obj

! ============================== added by K. Tagami ============= 11.0
  integer :: nloc_for_l
!
  real(kind=DP), allocatable ::  jpsmax(:) ! d(ntyp) # valence orbitals
  integer, allocatable ::  jtp(:,:) ! d(nlmt,ntyp) mapping F. ilmt to j.
  integer, allocatable ::  jtp_add(:,:)
  integer, allocatable ::  jtp_phi(:,:)
  integer, allocatable ::  jtp_pao(:,:)
!
  integer, allocatable ::  kjtp(:,:) ! d(nlmt,ntyp) mapping F. ilmt to j.
  integer, allocatable ::  kjtp_add(:,:)
  integer, allocatable ::  kjtp_phi(:,:)
  integer, allocatable ::  kjtp_pao(:,:)
!
  complex(kind=CMPLDP), allocatable ::  q_noncl(:,:,:,:)
  complex(kind=CMPLDP), allocatable ::  dion0_noncl(:,:,:,:)

  complex(kind=CMPLDP), allocatable ::  dion_scr_noncl(:,:,:,:)
!
  complex(kind=CMPLDP), allocatable ::  fqwei_noncl(:,:)
  complex(kind=CMPLDP), allocatable ::  fqwei_p_noncl(:,:)
  complex(kind=CMPLDP), allocatable ::  dion_so(:,:,:,:)
!
!  real(kind=DP), allocatable ::  q_so(:,:,:,:)
!
!
  logical :: flg_soc
  logical, allocatable :: pot_has_soc(:)
!
  integer, allocatable :: nums_of_angmom_on_atomtype(:)
!
!
! ============================================================== 11.0

! ============================ added by K. Tagami ========== 11.0
  real(kind=DP), allocatable :: Mat_SOC_Strength_nonpaw(:,:,:,:)
! ========================================================== 11.0

! ======= KT_add ======== 13.0U2
  real(kind=DP), allocatable ::  dion_paw_old(:,:,:,:)
  complex(kind=CMPLDP), allocatable ::  dion_scr_noncl_old(:,:,:,:)
! ======================= 13.0U2

! ==== KT_add === 2014/09/19
  real(kind=DP), allocatable :: kins_qrs(:), kina_qrs(:)
  real(kind=DP), allocatable :: kins_qrsps_mm(:,:), kina_qrsps_mm(:,:)
  real(kind=DP), allocatable :: kins_qitg_l(:,:), kina_qitg_l(:,:)
! =============== 2014/09/19

  type radial_aug_charge
     integer, allocatable :: nmm_il3(:), mm_il3(:,:)
     real(kind=DP), allocatable :: qrsps_mm(:,:)
  end type radial_aug_charge
  type(radial_aug_charge), allocatable, target :: radial_aug_chg(:)

  real(kind=DP), allocatable :: q_phirt_pw(:,:,:,:)  ! (nloc,ntau,ntau,ntyp)

  integer, allocatable, dimension(:)       :: nmm_il3 !d(lcmax+1)
  integer, allocatable, dimension(:,:)     :: mm_il3 !d(nqitg_sp,lcmax+1)
  real(kind=DP), allocatable, dimension(:,:) :: qrsps_mm !d(nmesh(it),nqitg_sp(it))

!  include 'mpif.h'

! subroutines and functions contained in this module
!   1. m_PP_wd_schgpc              <- (PseudoPotential_Construction)
!   2. m_PP_wd_PP_parameters       <- (WriteDownData_onto_Files)
!     { 2.1. wd_kgp_array,     2.2. wd_kgp_array_p}
!   3. m_PP_wd_betar
!   4. m_PP_rd_betar
!   5. m_PP_rd_PP_parameters       <- (PseudoPotential_Construction)
!     { 5.1. rd_kgp_array,     5.2. rd_kgp_array_p,   5.3. bcast_nfcntn_bin}
!   6.*what_is_ps_xctype
!   7. gga_grad_alloc              <- (32.3)
!   8. gga_grad_dealloc            <- (32.3)
!   9. m_PP_alloc0_ps_ntyp         <- (PsuedoPotential_Construction)
!  10. m_PP_alloc_ps_ntyp          <- (PseudoPotential_Construction)
!  11. m_PP_dealloc_ps_ntyp        <- (PseudoPotential_Construction)
!  12. m_PP_alloc_ps0              <- (PseudoPotential_Construction)
!  13. m_PP_alloc_ps               <- (PseudoPotential_Construction)
!  14. m_PP_dealloc_ps             <- (PseudoPotential_Construction)
!  15. m_PP_alloc_NLP              <- (PseudoPotential_Construction)
!  16. m_PP_dealloc_NLP            <- (PseudoPotential_Construction)
!  17. alloc_betar                 <- (4)
!  18.*m_PP_betar_calculated       <- (PseudoPotential_Construction)
!  19. m_PP_alloc_ps_stress        not used
!  20. m_PP_alloc_psc_qitg_rhpcg   <- (PseudoPotential_Construction)
!  21. m_PP_dealloc_psc_qitg_rhpcg <- (Finalization_of_mpi)
!  22. m_PP_alloc_qitg_wan         <- (PseudoPotential_Construction)
!  23. m_PP_set_dk_wan             <- (PseudoPotential_Construction)
!  24. m_PP_partial_core_CD        <- (PseudoPotential_Construction)
!      - epseud_pc, -pcc, -pcc_diff
!  25. m_PP_local_part             <- (PseudoPotential_Construction)
!      - add_localPot_to_vlocr, - localPP_in_Gspace,
!      - localPP_diff_in_Gspace
!      - add_localPotEnergy_to_etot1
!  26. sum_of_vlocr                <- (32.2)
!  27. sum_of_vvv_and_vxc          <- (32.2)
!  28. m_PP_check_gncpp_type       <- (PseudoPotential_Construction)
!  29. ft_valence_charge           <- (32)
!  30. where_does_WFrhvr_damp      <- (32)
!  31. make_index_arrays           <- (5)
!  32. m_PP_vanderbilt_type (nfp,it,nfout,gr_l,paramset) <- (PseudoPotential_Construction)
!      - 32.1  rd_Vloc_psWF_Q_then_chir_qitg
!      - 32.2  rd_itau_etc_then_phir_chir
!      - 32.3  make_index_arrays_nlmt2l_m_t
!      - 32.4  make_index_arrays_nlmt_add2lmt
!      - 32.5  make_index_arrays_nlmt_phi2lmt
!      - 32.6  make_index_arrays_nlmt_pao2lmt
!      - 32.7  cnstrct_of_localWF_Dion_and_q2
!      - 32.8  coulomb_potential_in_Rspace
!      - 32.9  atomic_xc_potential_in_Rspace
!      - 32.10 vlocr_plus_hartree_xc2
!      - 32.11 read_natomn_ival_iloc_itpcc
!      - 32.12 read_ps_xctype()
!      - 32.13 read_alp_cc
!      - 32.14 read_nmesh_xh_rmax
!      - 32.15 determine_lpsmax
!      - 32.16*check_vall
!      - 32.16 smoothing_vlocr()
!    { - rd_Vloc_psWF_Q_then_chir_qitg()  =(32.1)
!            - read_natomn_ival_iloc_itpcc() -(32.11) -> ival(it),iloc(it),itpcc(it) from the file of nfp
!            - read_ps_xctype()              -(32.12)  check of xc-potential type
!            - read_alp_cc()        -(32.13) -> alp(1:2,it),cc(1:2,it)
!            - read_nmesh_xh_rmax() -(32.14) -> nmesh(it),xh(it),rmax(it)
!            - check_vall           -(32.16)
!            - reading vvv          -> vvv(1:nmesh(it))
!            - reading vlocr        -> vloc(1:nmesh(it))
!            - reading rhvr         -> rhvr(1:nmesh(it))
!            - rmeshs() !-(b_P.P.)  -> radr(1:mmesh)
!            - rd_itau_etc_then_phir_chir() -(32.2) -> phir, chir
!                  - read_itau_ivanl        -(b_P.P.) -> itau, ivanl
!                  - read_tau_eps_nrc_mord  -(b_P.P.) -> eps, t1, it
!                  - cnvrtp()        -(b_P.P._f77)
!                  - cnvrtc()        -(b_P.P._f77)
!            - coef_simpson_integeration()      -(b_P.P.)
!            - rd_qrsps_then_iqitg_and_qitgft() -(33)
!                  - read_nrc_mord() -(b_P.P.)
!                  - cnvrtp()        -(b_P.P._f77)
!                  - qitgft()        -(b_P.P.)
!                       - dsjnv()    -(bottom_Subroutines)
!                  - qitgft_diff()   -(b_P.P.)
!                       - dsjnv()    -(bottom_Subroutines)
!                  - qij_qvij_from_qrsps_etc() -(34) -> qij, qvij
!      - make_index_arrays_nlmt2l_m_t()  -(32.3)  -> ltp, mtp, taup, index_lmt1_lmt2, w_non0_lmtxlmt
!      - cnstrct_of_localWF_Dion_and_q() -(32.7)
!            - matrix_inversion()   -(b_P.P.)
!      - where_does_WFrhvr_damp       -(30)
!      - coulomb_potential_in_Rspace  -(32.8) -> vvv
!      - atomic_xc_potential_in_Rspace -(32.9)
!            - gga_grad_alloc         -(7)
!            - get_gradient_of_rho -(b_P.P.)
!                  - getrho
!                  - cpval
!                  - gdiffs
!                  - gtgrad
!                  - cnggrd
!            - gga_xc_potential_atomic -(b_P.P.)
!                  - ggaexch_pw91
!                  - ggacorl_pw91
!                  - ggaexch_pbe
!                  - ggacorl_pbe
!                  - ggabp
!            - gga_grad_dealloc       -(8)
!            - xcpot              -(b_P.P._f77)
!      - vlocr_plus_hartree_xc        -(32.10)
!            - sum_of_vlocr           -(26)
!            - sum_of_vvv_and_vxc     -(27)
!      - smoothing_vlocr              -(32.16)
!     }
!  33. rd_qrsps_then_iqitg_and_qitgft
!      - qitgft_mm,      - wd_nqitg_sp_etc
!  34. qij_qvij_from_qrsps_etc       <- (33)
!  35. wd_index_arrays_etc           <- (32.7)
!  36. m_PP_wd_variables             <- (PseudoPotential_Construction)
!  37. m_PP_tell_lmtt_l_m_tau        <- (m_NonLocal_Potential)
!  38. m_PP_tell_lmtt_l_m_tau_phi    <- (m_NonLocal_Potential)
!  39. m_PP_tell_lmtt_l_m_tau_add    <- (m_NonLocal_Potential)
!  40. m_PP_tell_lmtt_l_m_tau_pao    <- (m_NonLocal_Potential)
!  41. m_PP_make_index_lmtt_2_dl2p   <- (PseudoPotential_Construction)
!  42. m_PP_make_index_lmtt_phi      <- (PseudoPotential_Construction)
!  43. m_PP_make_index_lmtt_add      <- (PseudoPotential_Construction)
!  44. m_PP_make_index_lmtt_pao      <- (PseudoPotential_Construction)
!  45. heap_sort_wrt_rank            <- (44)
!  46. m_PP_find_maximum_l            <- (m_Charge_Density, m_Force,
!                          m_ES_Intgr_VlhxcQlm, m_Stress, m_XC_Potential)
!  47. *m_PP_include_vanderbilt_pot   <-(48, m_Force, m_ES_WF_by_RMM,
!                                         m_Electronic_Structure, m_ES_WF_by_SDorCG,
!                                         m_Stress, m_XC_Potential)
!  48. m_PP_gfqwei()                  <-(Initial_Electronic_Structure)
!  49. *m_PP_max_of_itau                   not used
!  50. m_PP_set_mmesh                      <-(PseudoPotential_Construction)
!  51. m_PP_set_nloc                       <-(PseudoPotential_Construction)
!  52. m_PP_set_m_non0_lmtxlmt             <-(PseudoPotential_Construction)
!  53. m_PP_tell_iorb_ia_l_m_tau           <-(m_ES_dos)
!  54. m_PP_tell_iorb_lmt                  <-(m_ES_dos, m_Electronic_Structure)
!  55. m_PP_make_qorb                      <-(PseudoPotential_Construction)
!  56. make_phi_window_parameter           <-(32.2)
!  57. m_PP_rd_window_param                <-(InputData_Analysis)
!  58. m_PP_cnstrct_crotylm                <-(PseudoPotential_Construction)
!  59. m_PP_ps_xctype                      not used
!  60. m_PP_input_xctype(xctype_is)        <-(Preparation)
!  61. m_PP_set_index_arrays1              <-(m_Charge_Density, m_ES_Intgr_VlhxcQlm, m_Force)
!  62. m_PP_set_index_arrays1              <-(m_Charge_Density, m_ES_Intgr_VlhxcQlm, m_Force)

contains
  subroutine m_PP_wd_schgpc(nfout)
    integer, intent(in) :: nfout
    if(ipripp >= 2) write(nfout,320)  schgpc
320 FORMAT(' !PP TOTAL PARTIAL CHARGE         SCHGPC = ',F20.10)
  end subroutine m_PP_wd_schgpc


  subroutine m_PP_wd_betar(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin
    if(mype==0) then
       if(allocated(betar)) write(nfcntn_bin) betar
    end if
  end subroutine m_PP_wd_betar

  subroutine m_PP_rd_betar(nfcntn_bin)
    integer, intent(in) :: nfcntn_bin
    if(.not.allocated(betar)) call alloc_betar()
    if(mype==0) read(nfcntn_bin) betar
    call mpi_bcast(betar,mmesh*nloc*ntau*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
  end subroutine m_PP_rd_betar


  integer function what_is_ps_xctype()
    if(    ps_xctype == 'GGAPW91' .or. ps_xctype == 'ggapw91' .or.&
         & ps_xctype == 'GGAPBE ' .or. ps_xctype == 'ggapbe ' .or.&
         & ps_xctype == 'GGAPBEX' .or. ps_xctype == 'ggapbex' .or.&
         & ps_xctype == 'GGAPBEY' .or. ps_xctype == 'ggapbey' .or.&
         & ps_xctype == 'VDWDF  ' .or. ps_xctype == 'vdwdf  ' .or.&
! ===> 13.0A
         & ps_xctype == 'REVPBE ' .or. ps_xctype == 'revpbe ' .or.&
         & ps_xctype == 'RPBE   ' .or. ps_xctype == 'rpbe   ' .or.&
         & ps_xctype == 'WC06   ' .or. ps_xctype == 'wc06   ' .or.&
         & ps_xctype == 'HTBS   ' .or. ps_xctype == 'htbs   ' .or.&
         & ps_xctype == 'PBESOL ' .or. ps_xctype == 'pbesol ' .or.&
         & ps_xctype == 'PBEINT ' .or. ps_xctype == 'pbeint ' .or.&
! <=== 13.0A
         & ps_xctype == 'GGAPBEK' .or. ps_xctype == 'ggapbek' .or.&
         & ps_xctype == 'KATOPBE' .or. ps_xctype == 'katopbe' .or.&
         & ps_xctype == 'GGABP  ' .or. ps_xctype == 'ggabp  ' .or.&
         & ps_xctype == 'LDAPW91' .or. ps_xctype == 'ldapw91' .or.&
         & ps_xctype == 'LDAPBE ' .or. ps_xctype == 'ldapbe ' .or.&
         & ps_xctype == 'LDAPBEK' .or. ps_xctype == 'ldapbek' ) then
       what_is_ps_xctype = GGA
    else
       what_is_ps_xctype = LDA
    end if
  end function what_is_ps_xctype

  integer function what_is_input_xctype()
    if(    xctype == 'GGAPW91' .or. xctype == 'ggapw91' .or.&
         & xctype == 'GGAPBE ' .or. xctype == 'ggapbe ' .or.&
         & xctype == 'GGAPBEX' .or. xctype == 'ggapbex' .or.&
         & xctype == 'GGAPBEY' .or. xctype == 'ggapbey' .or.&
         & xctype == 'VDWDF  ' .or. xctype == 'vdwdf  ' .or.&
! ==> 13.0A
         & xctype == 'REVPBE ' .or. xctype == 'revpbe ' .or.&
         & xctype == 'RPBE   ' .or. xctype == 'rpbe   ' .or.&
         & xctype == 'WC06   ' .or. xctype == 'wc06   ' .or.&
         & xctype == 'HTBS   ' .or. xctype == 'htbs   ' .or.&
         & xctype == 'PBESOL ' .or. xctype == 'pbesol ' .or.&
         & xctype == 'PBEINT ' .or. xctype == 'pbeint ' .or.&
         & xctype == 'EV93   ' .or. xctype == 'ev93   ' .or.&
         & xctype == 'EVPW91 ' .or. xctype == 'evpw91 ' .or.&
! <== 13.0A
         & xctype == 'GGAPBEK' .or. xctype == 'ggapbek' .or.&
         & xctype == 'KATOPBE' .or. xctype == 'katopbe' .or.&
         & xctype == 'GGABP  ' .or. xctype == 'ggabp  ' .or.&
         & xctype == 'LDAPW91' .or. xctype == 'ldapw91' .or.&
         & xctype == 'LDAPBE ' .or. xctype == 'ldapbe ' .or.&
         & xctype == 'LDAPBEK' .or. xctype == 'ldapbek' ) then
       what_is_input_xctype = GGA
    else
       what_is_input_xctype = LDA
    end if
  end function what_is_input_xctype

!!$ASASASASAS
!!$  subroutine gga_grad_alloc()
!!$    allocate(fdiff(mmesh,0:mddiff))
!!$    allocate(rho(mmesh))
!!$    allocate(coeff1(0:mddiff,mddiff))
!!$    allocate(coeff2(0:mddiff,2:mddiff))
!!$    allocate(grdnts(mmesh,3))
!!$  end subroutine gga_grad_alloc

  subroutine gga_grad_alloc()
    allocate(fdiff(mmesh,0:mddiff)); fdiff = 0
    allocate(rho(mmesh)); rho = 0
    allocate(coeff1(0:mddiff,mddiff)); coeff1 = 0
    allocate(coeff2(0:mddiff,2:mddiff)); coeff2 = 0
    allocate(grdnts(mmesh,3)); grdnts = 0
  end subroutine gga_grad_alloc
!!$ASASASASAS

  subroutine gga_grad_dealloc()
    deallocate(fdiff)
    deallocate(rho)
    deallocate(coeff1)
    deallocate(coeff2)
    deallocate(grdnts)
  end subroutine gga_grad_dealloc

!!$ASASASASAS
!!$  subroutine gga_grad_alloc2()
!!$    allocate(grdnts(mmesh,0:2))
!!$    allocate(wwk(mmesh,7))
!!$  end subroutine gga_grad_alloc2

  subroutine gga_grad_alloc2()
    allocate(grdnts(mmesh,0:2)); grdnts = 0
    allocate(wwk(mmesh,7)); wwk = 0
  end subroutine gga_grad_alloc2
!!$ASASASASAS

  subroutine gga_grad_dealloc2()
    deallocate(grdnts)
    deallocate(wwk)
  end subroutine gga_grad_dealloc2

!!$ASASASASAS
!!$  subroutine m_PP_alloc0_ps_ntyp
!!$    allocate(ival(ntyp))
!!$    allocate(iloc(ntyp))
!!$    allocate(lpsmax(ntyp))
!!$    allocate(itpcc(ntyp))
!!$    allocate(ilmt(ntyp))
!!$    if(sw_orb_popu == ON) then
!!$       allocate(ilmt_phi(ntyp))
!!$    end if
!!$    if(sw_use_add_proj == ON) then
!!$       allocate(ilmt_add(ntyp))
!!$    end if
!!$    if(intzaj == by_pseudo_atomic_orbitals) then
!!$       allocate(ilmt_pao(ntyp))
!!$    end if
!!$    allocate(nmesh(ntyp))
!!$    allocate(xh(ntyp));  allocate(rmax(ntyp));   allocate(h(ntyp))
!!$    allocate(chgpc(ntyp))
!!$    allocate(alp(2,ntyp)); allocate(cc(2,ntyp))
!!$    allocate(n_non0_lmtxlmt(ntyp))
!!$    allocate(ae_wavefunctions_are_detected(ntyp)); ae_wavefunctions_are_detected = .false.
!!$    allocate(nqitg_sp(ntyp)); nqitg_sp = 0
!!$  end subroutine m_PP_alloc0_ps_ntyp

  subroutine m_PP_alloc0_ps_ntyp
    allocate(ival(ntyp)); ival = 0
    allocate(iloc(ntyp)); iloc = 0
    allocate(lpsmax(ntyp)); lpsmax = 0
    allocate(itpcc(ntyp)); itpcc = 0
    allocate(ilmt(ntyp)); ilmt = 0
    if(sw_orb_popu == ON) then
       allocate(ilmt_phi(ntyp)); ilmt_phi = 0
    end if
    if(sw_use_add_proj == ON) then
       allocate(ilmt_add(ntyp)); ilmt_add = 0
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       allocate(ilmt_pao(ntyp)); ilmt_pao = 0
    end if
    allocate(nmesh(ntyp)); nmesh = 0
    allocate(xh(ntyp));  allocate(rmax(ntyp));   allocate(h(ntyp))
    xh = 0.0; rmax = 0.0; h = 0.0

    allocate(chgpc(ntyp))
    allocate(alp(2,ntyp)); allocate(cc(2,ntyp))
    chgpc = 0; alp = 0; cc = 0

    allocate(n_non0_lmtxlmt(ntyp)); n_non0_lmtxlmt = 0
    allocate(ae_wavefunctions_are_detected(ntyp)); ae_wavefunctions_are_detected = .false.
    allocate(nqitg_sp(ntyp)); nqitg_sp = 0

    allocate(ipaw(ntyp));ipaw=0
    allocate(chgcr(ntyp))
    allocate(iltpw(ntyp));iltpw=0

    if ( sw_hybrid_functional == ON .and. sw_rspace_hyb == OFF ) then
       allocate( radial_aug_chg(ntyp) )
    endif

! ================================= added by K. Tagami ================= 11.0
    if(.not.allocated(pot_has_soc)) allocate( pot_has_soc(ntyp) ); pot_has_soc = .false.
    if(.not.allocated(nums_of_angmom_on_atomtype)) allocate( nums_of_angmom_on_atomtype(ntyp) )
    nums_of_angmom_on_atomtype = 0
! ===================================================================== 11.0

  end subroutine m_PP_alloc0_ps_ntyp
!!$ASASASASAS

  subroutine m_PP_alloc_ps_ntyp(paramset)
    logical, intent(in) :: paramset
    integer :: nloc_t, ntau_t, nlmt_t, nlmt_phi_t, nlmt_add_t, nlmt_pao_t
    integer :: nltpw_t
    if(paramset) then
       nloc_t = kloc_t; ntau_t = ktau_t; nlmt_t = klmt_t
       if(sw_orb_popu == ON) then
          nlmt_phi_t = klmt_t
       end if
       if(sw_use_add_proj == ON) then
          nlmt_add_t = klmt_t
       end if
       if(intzaj == by_pseudo_atomic_orbitals) then
          nlmt_pao_t = klmt_t
       end if
       nltpw_t=klmt_t
       nltpw=0
    else
       nloc_t = nloc; ntau_t = ntau; nlmt_t = nlmt
       if(sw_orb_popu == ON) then
          nlmt_phi_t = nlmt_phi
       end if
       if(sw_use_add_proj == ON) then
          nlmt_add_t = nlmt_add
       end if
       if(intzaj == by_pseudo_atomic_orbitals) then
          nlmt_pao_t = nlmt_pao
       end if
       nltpw_t=nltpw
    end if
    !!$ print '(" -- m_PP_alloc_ps_ntyp -- ")'
    !!$ print '(" nloc_t, ntau_t, nlmt_t = ",3i5)',nloc_t,ntau_t,nlmt_t
    allocate(ivanl(nloc_t,ntyp)); ivanl = 0
    allocate(itau(nloc_t,ntyp)); itau = 0
    allocate(eps(nloc_t,ntau_t,ntyp)); eps = 0.d0
    allocate(lmtt(nlmt_t,ntyp)); lmtt = 0
    allocate(ltp(nlmt_t,ntyp)); ltp = 0
    allocate(mtp(nlmt_t,ntyp)); mtp = 0
    allocate(taup(nlmt_t,ntyp)); taup = 0
    if(sw_orb_popu == ON) then
       allocate(lmtt_phi(nlmt_phi_t,ntyp)); lmtt_phi=0
       allocate(ltp_phi(nlmt_phi_t,ntyp)); ltp_phi=0
       allocate(mtp_phi(nlmt_phi_t,ntyp)); mtp_phi=0
       allocate(taup_phi(nlmt_phi_t,ntyp)); taup_phi=0
       allocate(rc_phi(nloc_t,ntau_t,ntyp)); rc_phi=0.d0
       allocate(k_phi(nloc_t,ntau_t,ntyp)); k_phi=0.d0
       allocate(iproj_phi(nlmt_phi_t,ntyp)); iproj_phi=0
    end if
    if(sw_use_add_proj == ON) then
       allocate(lmtt_add(nlmt_add_t,ntyp)); lmtt_add=0
       allocate(ltp_add(nlmt_add_t,ntyp)); ltp_add=0
       allocate(mtp_add(nlmt_add_t,ntyp)); mtp_add=0
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       allocate(lmtt_pao(nlmt_pao_t,ntyp)); lmtt_pao=0
       allocate(ltp_pao(nlmt_pao_t,ntyp)); ltp_pao=0
       allocate(mtp_pao(nlmt_pao_t,ntyp)); mtp_pao=0
       allocate(taup_pao(nlmt_pao_t,ntyp)); taup_pao=0
       allocate(irank(nloc_t,ntau_t,ntyp)); irank = 0
    end if

! === For continue_bin_paw.data by T.Kato 2011/03/28 ===========================
!   allocate(lppw(nltpw_t,ntyp))
!   allocate(tppw(nltpw_t,ntyp))
    allocate(lppw(nltpw_t,ntyp)); lppw = 0
    allocate(tppw(nltpw_t,ntyp)); tppw = 0
! ==============================================================================

  end subroutine m_PP_alloc_ps_ntyp

  subroutine m_PP_dealloc_ps_ntyp
    deallocate(ivanl);deallocate(itau);deallocate(eps);deallocate(lmtt)
    deallocate(ltp);deallocate(mtp);deallocate(taup)
    deallocate(lppw);deallocate(tppw)
    if(sw_orb_popu == ON) then
       deallocate(lmtt_phi)
       deallocate(ltp_phi);deallocate(mtp_phi);deallocate(taup_phi)
       deallocate(rc_phi);deallocate(k_phi)
       deallocate(iproj_phi)
    end if
    if(sw_use_add_proj == ON) then
       deallocate(lmtt_add)
       deallocate(ltp_add);deallocate(mtp_add)
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       deallocate(lmtt_pao)
       deallocate(ltp_pao)
       deallocate(mtp_pao)
       deallocate(taup_pao)
       deallocate(irank)
    end if
  end subroutine m_PP_dealloc_ps_ntyp

!!$ASASASASAS
!!$  subroutine m_PP_alloc_ps0
!!$    !!$ print '(" -- m_PP_alloc_ps0 --")'
!!$    !!$ print '(" nloc, ntau, ntyp, nlmt, natm = ",5i5)',nloc,ntau,ntyp,nlmt,natm
!!$    allocate(lmta(nlmt,natm)); lmta = 0
!!$    if(sw_orb_popu == ON) then
!!$       allocate(lmta_phi(nlmt_phi,natm)); lmta_phi = 0
!!$       allocate(q_phi(nlmt_phi,ntyp)); q_phi = 0.d0
!!$    end if
!!$    if(sw_use_add_proj == ON) then
!!$       allocate(lmta_add(nlmt_add,natm)); lmta_add = 0
!!$    end if
!!$    if(intzaj == by_pseudo_atomic_orbitals) then
!!$       allocate(lmta_pao(nlmt_pao,natm)); lmta_pao = 0
!!$       allocate(ibpao(nlmta_pao)); ibpao = 0
!!$    end if
!!$!-->
!!$    allocate(q(nlmt,nlmt,ntyp))
!!$    allocate(if_q_isnotzero(nlmt,nlmt,ntyp))
!!$    allocate(dion(nlmt,nlmt,ntyp))
!!$    allocate(isph(nlmt,nlmt,6,ntyp))
!!$    allocate(il2p(nlmt,nlmt,ntyp))
!!$    allocate(dl2p(nlmt,nlmt,6,ntyp))
!!$    allocate(prodphi(num_projectors,ntau,ntau))

!!$    allocate(index_lmt1_lmt2(m_non0_lmtxlmt,ntyp,2))
!!$    allocate(w_non0_lmtxlmt(m_non0_lmtxlmt,ntyp))
!!$    allocate(q_indp(m_non0_lmtxlmt,ntyp))
!!$    allocate(dion_indp(m_non0_lmtxlmt,ntyp))
!!$!<--
!!$    allocate(iqitg(nloc,ntau,nloc,ntau,nloc*2,ntyp)); iqitg = 0
!!$    allocate(fqwei(nlmt*ntau*natm))
!!$    allocate(nlmta1(nlmt*ntau*natm)); allocate(nlmta2(nlmt*ntau*natm)); allocate(nac2ia(nlmt*ntau*natm))
!!$  end subroutine m_PP_alloc_ps0

  subroutine m_PP_alloc_ps0
    !!$ print '(" -- m_PP_alloc_ps0 --")'
    !!$ print '(" nloc, ntau, ntyp, nlmt, natm = ",5i5)',nloc,ntau,ntyp,nlmt,natm
    allocate(lmta(nlmt,natm)); lmta = 0
    if(sw_orb_popu == ON) then
       allocate(lmta_phi(nlmt_phi,natm)); lmta_phi = 0
       allocate(q_phi(nlmt_phi,ntyp)); q_phi = 0.d0
       if( orb_popu_method == 2 ) then
          allocate( q_phirt_pw(nloc,ntau,ntau,ntyp)); q_phirt_pw = 0.d0
       end if
    end if
    if(sw_use_add_proj == ON) then
       allocate(lmta_add(nlmt_add,natm)); lmta_add = 0
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       allocate(lmta_pao(nlmt_pao,natm)); lmta_pao = 0
       allocate(ibpao(nlmta_pao)); ibpao = 0
    end if
!-->
    allocate(q(nlmt,nlmt,ntyp)); q=0
    allocate(if_q_isnotzero(nlmt,nlmt,ntyp)); if_q_isnotzero = 0
    allocate(dion(nlmt,nlmt,ntyp)); dion = 0.0
    allocate(isph(nlmt,nlmt,6,ntyp)); isph = 0
    allocate(il2p(nlmt,nlmt,ntyp)); il2p = 0
    allocate(dl2p(nlmt,nlmt,6,ntyp)); dl2p = 0
    allocate(prodphi(num_projectors,ntau,ntau)); prodphi = 0.0

    allocate(index_lmt1_lmt2(m_non0_lmtxlmt,ntyp,2)); index_lmt1_lmt2=0
    allocate(w_non0_lmtxlmt(m_non0_lmtxlmt,ntyp)); w_non0_lmtxlmt = 0
    allocate(q_indp(m_non0_lmtxlmt,ntyp)); q_indp = 0
    allocate(dion_indp(m_non0_lmtxlmt,ntyp)); dion_indp = 0
!<--
    allocate(iqitg(nloc,ntau,nloc,ntau,nloc*2,ntyp)); iqitg = 0
    allocate(fqwei(nlmt*ntau*natm)); fqwei = 0

    allocate(nlmta1(nlmt*ntau*natm));  nlmta1 = 0
    allocate(nlmta2(nlmt*ntau*natm));  nlmta2 = 0
    allocate(nac2ia(nlmt*ntau*natm)); nac2ia = 0

    allocate(qrspspw(mmesh,0:nqitg));qrspspw=0

! === KT_add === 2014/09/19 &13.1XI &2015/02/23
    if (  .not. flg_paw ) then
       if ( sw_rspace_ekin_density == ON ) then
          allocate(psirpw(mmesh,nloc,ntau,ntyp));psirpw=0.d0
          allocate(phirpw(mmesh,nloc,ntau,ntyp));phirpw=0.d0
          allocate(wf_nrc(nloc,ntau,ntyp));wf_nrc=0
          allocate(wf_mnrc(ntyp));wf_mnrc=0
       else if ( sw_excitation == ON .or. sw_wannier90 == ON ) then
          allocate(psirpw(mmesh,nloc,ntau,ntyp));psirpw=0.d0
          allocate(phirpw(mmesh,nloc,ntau,ntyp));phirpw=0.d0
       else if ( sw_orb_popu == ON ) then
          allocate(wf_nrc(nloc,ntau,ntyp));wf_nrc=0
       else if ( sw_add_corecharge_rspace == ON ) then
          allocate(radr_paw(mmesh,ntyp));radr_paw=0.d0
          allocate(rhcorpw(mmesh,ntyp));rhcorpw=0.d0
       endif
       if ( sw_hubbard == ON .and. initial_occmat == CRYSTAL_FIELD_APPROX ) then
          if (.not.allocated(psirpw)) allocate(psirpw(mmesh,nloc,ntau,ntyp))
          if (.not.allocated(radr_paw)) allocate(radr_paw(mmesh,ntyp))
          psirpw=0.d0; radr_paw=0.d0
       endif
       if ( sw_write_rotated_orbitals == ON ) then
          if (.not.allocated(radr_paw)) allocate(radr_paw(mmesh,ntyp))
          radr_paw=0.d0
       endif
       if ( use_metagga ) then
          allocate(index_lmt2lt(nlmt,ntyp));index_lmt2lt=0
       endif
    endif
! ============== 2014/09/19 &13.1XI &2015/02/23

    if(flg_paw) then
        allocate(psirpw(mmesh,nloc,ntau,ntyp));psirpw=0.d0
        allocate(phirpw(mmesh,nloc,ntau,ntyp));phirpw=0.d0
!        allocate(qrspspw(mmesh,0:nqitg));qrspspw=0
!        allocate(rhcor(mmesh))
!        allocate(psir(mmesh,nloc,ntau))
        allocate(rhcorpw(mmesh,ntyp));rhcorpw=0.d0
        allocate(rhpcrpw(mmesh,ntyp));rhpcrpw=0.d0
        allocate(wf_nrc(nloc,ntau,ntyp));wf_nrc=0
        allocate(wf_mnrc(ntyp));wf_mnrc=0
        if(.not.flg_symmtry) then
            allocate(n_cijkclmk(nlmt,nlmt,ntyp));n_cijkclmk=0
            allocate(CijkClmkVVVVijlm_k(nlmt,nlmt,m_clmns_cijkclmk,ntyp))
            CijkClmkVVVVijlm_k=0.d0
            allocate(CijkClmkVVVVijlm_k_ae(nlmt,nlmt,m_clmns_cijkclmk,ntyp))
            CijkClmkVVVVijlm_k_ae=0.d0
            allocate(ilmt3_cijkclmk(nlmt,nlmt,m_clmns_cijkclmk,ntyp))
            ilmt3_cijkclmk=0
            allocate(ilmt4_cijkclmk(nlmt,nlmt,m_clmns_cijkclmk,ntyp))
            ilmt4_cijkclmk=0
        else
            allocate(n_cijkclmn(nlmt,nlmt,natm,nopr))
            n_cijkclmn=0
            allocate(CijkClmnVVVVijlm_kn &
                                (nlmt,nlmt,m_clmns_cijkclmn,natm,nopr))
            CijkClmnVVVVijlm_kn=0.d0
            allocate(CijkClmnVVVVijlm_kn_ae &
                                (nlmt,nlmt,m_clmns_cijkclmn,natm,nopr))
            CijkClmnVVVVijlm_kn_ae=0.d0
            allocate(ilmt3_cijkclmn(nlmt,nlmt,m_clmns_cijkclmn,natm,nopr))
            ilmt3_cijkclmn=0
            allocate(ilmt4_cijkclmn(nlmt,nlmt,m_clmns_cijkclmn,natm,nopr))
            ilmt4_cijkclmn=0
        end if

        allocate(index_lmt2lt(nlmt,ntyp));index_lmt2lt=0
        allocate(dion_kin_ion(nlmt,nlmt,ntyp));dion_kin_ion=0.d0
        allocate(dion_hartree(nlmt,nlmt,natm));dion_hartree=0.d0
        allocate(dion_hartree_now(nlmt,nlmt,natm));dion_hartree_now=0.d0

! =============================== modified by K. Tagami =================== 11.0
!        allocate(dion_vxc(nlmt,nlmt,nspin,natm));dion_vxc=0.d0
!        allocate(dion_paw(nlmt,nlmt,nspin,natm));dion_paw=0.d0

        if ( noncol ) then
           allocate(dion_vxc(nlmt,nlmt,ndim_magmom,natm));dion_vxc=0.d0
           allocate(dion_paw(nlmt,nlmt,ndim_magmom,natm));dion_paw=0.d0
        else
           allocate(dion_vxc(nlmt,nlmt,nspin,natm));dion_vxc=0.d0
           allocate(dion_paw(nlmt,nlmt,nspin,natm));dion_paw=0.d0
        endif
! ======================================================================== 11.0

! ====== KT_add ====== 13.0U2
        if ( sw_potential_mixing == ON ) then
           if ( noncol ) then
              allocate(dion_paw_old(nlmt,nlmt,ndim_magmom,natm));dion_paw_old=0.d0
           else
              allocate(dion_paw_old(nlmt,nlmt,nspin,natm));dion_paw_old=0.d0
           endif
        endif
! ==================== 13.0U2

        allocate(radr_paw(mmesh,ntyp));radr_paw=0.d0
    end if

! ================================ added by K. Tagami =========== 11.0
!!!!!    call determine_flg_soc

    flg_soc = .false.

! ------------------------------------------------
    if ( flg_soc ) then
       noncol = .true.
                           ! forced to be true irrespective of its val.

!       allocate( dion_so( nlmt,nlmt, ndim_chgpot, ntyp ) ); dion_so = 0.0d0
!       allocate( q_so( nlmt,nlmt, nspin, ntyp ) ); q_so = 0.0d0
    endif

    if ( noncol ) then
       allocate( q_noncl( nlmt,nlmt, ndim_chgpot, ntyp ) )

!       allocate( dion0_noncl( nlmt,nlmt, ndim_chgpot, ntyp ) )
       allocate( dion0_noncl( nlmt,nlmt, ndim_chgpot, natm ) )

       allocate( dion_scr_noncl( nlmt,nlmt, ndim_chgpot, natm ) )
       dion0_noncl = 0.0d0;    dion_scr_noncl = 0.0d0;     q_noncl = 0.0d0

       allocate( fqwei_noncl( nlmt*ntau*natm, ndim_chgpot ));
       fqwei_noncl = 0.0d0

!!!!!!!!!       allocate( fqwei_p_noncl( nac_p, ndim_chgpot ));
!!!!!!!!!!      fqwei_p_noncl = 0.0d0
    endif

! =============================================================== 11.0

! ============================= added by K. Tagami ============== 11.0
    if ( allocated( Mat_SOC_strength_nonpaw ) ) deallocate( Mat_SOC_strength_nonpaw )

!    if ( .not. flg_paw ) then
       if ( SpinOrbit_Mode == ZeffApprox ) then
          allocate( Mat_SOC_strength_nonpaw( nloc,ntau,ntau,ntyp ) )
          Mat_SOC_Strength_nonpaw = 0.0d0
       endif
!    endif
     if ( SpinOrbit_Mode == ReadFromPP ) then
        allocate( Mat_SOC_strength_nonpaw( nloc,ntau,ntau,ntyp ) )
        Mat_SOC_Strength_nonpaw = 0.0d0
     endif
! ================================================================ 11.0

! === KT_add === 2014/08/11 & 13.2S
     if ( noncol .and. SpinOrbit_Mode == ByPawPot ) then
        if ( sw_use_rphi_Hsoc_rphi == ON .and. sw_use_ival_for_paw_ps_soc == OFF ) then
           allocate( vlocr_pw(mmesh,ntyp) );  vlocr_pw = 0.0d0
        endif
     else if ( sw_spinorbit_second_variation==ON .and. SpinOrbit_Mode == ByPawPot ) then
        if ( sw_use_rphi_Hsoc_rphi == ON .and. sw_use_ival_for_paw_ps_soc == OFF ) then
           allocate( vlocr_pw(mmesh,ntyp) );  vlocr_pw = 0.0d0
        endif
     else
        if ( sw_corelevel_spectrum == ON .or. sw_calc_core_energy == ON ) then
           allocate( vlocr_pw(mmesh,ntyp) );  vlocr_pw = 0.0d0
        endif
     endif
! ============== 2014/08/11 & 13.2S

  end subroutine m_PP_alloc_ps0
!!$ASASASASAS

! ============================ addded by K. Tagami ============= 11.0
!  subroutine determine_flg_soc
!    integer :: it

!    flg_soc = .false.
!    Do it=1, ntyp
!       if ( pot_has_soc(it) ) then
!          flg_soc = .true.
!       endif
!    End do
!  end subroutine determine_flg_soc

! ============================================================== 11.0

!!$  subroutine m_PP_alloc_ps
!!$    allocate(qij(ntau,ntau,nloc))
!!$    allocate(qvij(ntau,ntau,nloc))
!!$    allocate(vvv(mmesh))
!!$    allocate(vxc(mmesh))
!!$    allocate(exc(mmesh))
!!$    allocate(vlocr(mmesh))
!!$    allocate(vlocr2(mmesh))
!!$    allocate(rhvr(mmesh))
!!$    allocate(radr(mmesh))
!!$    allocate(qrsps(mmesh))
!!$    if(sw_berry_phase == ON .or. sw_wannier == ON) then
!!$       allocate(qrs(mmesh)); qrs = 0.d0
!!$    end if
!!$    allocate(wkx(mmesh)); wkx = 0.d0
!!$    allocate(wky(mmesh)); wky = 0.d0
!!$    allocate(wkz(mmesh)); wkz = 0.d0
!!$    allocate(rhpcr(mmesh));  rhpcr = 0.d0
!!$    if(sw_positron == BULK .or. sw_positron == DEFECT) then
!!$       allocate(rhcr(mmesh)); rhcr = 0.d0
!!$    end if
!!$    allocate(phir_ae(mmesh,nloc,ntau))
!!$    allocate(phir(mmesh,nloc,ntau))
!!$    allocate(chir(mmesh,nloc,ntau))
!!$    allocate(wos(mmesh))
!!$    allocate(copsc(0:kord)); copsc = 0.d0
!!$    allocate(s(ntau,ntau))
!!$    allocate(sinv(ntau,ntau))
!!$    !!$ print '(" mmesh, nloc,ntau,ntyp = ",4i5)', mmesh,nloc,ntau,ntyp
!!$    if(.not.allocated(betar)) allocate(betar(mmesh,nloc,ntau,ntyp))
!!$    if(sw_use_add_proj == 1) then
!!$       allocate(betar_add(mmesh,ntyp))
!!$       betar_add = 0.d0
!!$    end if
!!$! ----- T. Yamasaki, 3 July 2008 ---
!!$    if(ipripp >= 1)then
!!$       write(6,'(" sw_orb_popu = ",i5)') sw_orb_popu
!!$    end if
!!$    if(sw_orb_popu == 1) then
!!$       allocate(phirt(mmesh,nloc,ntau,ntyp))
!!$! ----- T. Yamasaki, 3 July 2008 ---
!!$       if(ipripp >= 1)then
!!$          write(6,'(" phirt is allocated")')
!!$          write(6,'("  --  mmesh, nloc, ntau, ntyp = ",4i6)') mmesh, nloc, ntau, ntyp
!!$       end if
!!$! ---------------------------------<<
!!$    end if
!!$    if(intzaj == by_pseudo_atomic_orbitals) then
!!$       allocate(paor(mmesh,nloc,ntau,ntyp))
!!$    end if
!!$    if(istress == 1) then
!!$       allocate(wkz1(mmesh)); wkz1 = 0.d0
!!$       allocate(wkz2(mmesh)); wkz2 = 0.d0
!!$    endif
!!$  end subroutine m_PP_alloc_ps

  subroutine m_PP_alloc_ps
    logical,save :: firstcall=.true.
    allocate(qij(ntau,ntau,nloc)); qij = 0
    allocate(qvij(ntau,ntau,nloc)); qvij = 0
    allocate(vvv(mmesh)); vvv = 0
    allocate(vxc(mmesh)); vxc = 0
    allocate(exc(mmesh)); exc = 0
    allocate(vlocr(mmesh)); vlocr = 0
    allocate(vlocr2(mmesh)); vlocr2 = 0
    allocate(rhvr(mmesh)); rhvr = 0
    allocate(radr(mmesh)); radr = 0
    allocate(qrsps(mmesh)); qrsps = 0
    if(sw_berry_phase == ON .or. sw_wannier == ON .or. sw_fef == ON) then
       allocate(qrs(mmesh)); qrs = 0.d0
    end if
! ==== KT_add === 2015/02/23
    if ( sw_wannier90 == ON .or. sw_LinearResponse == ON ) then
       if ( .not.allocated(qrs) ) then
          allocate(qrs(mmesh));  qrs = 0.0d0
       endif
    endif
! ============== 2015/02/23

! === KT_add === 2014/09/19
    if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
       if ( use_asymm_ekin_density ) then
          allocate( kina_qrs(mmesh) ); kina_qrs = 0.0d0
       endif
       if ( use_symm_ekin_density ) then
          allocate( kins_qrs(mmesh) ); kins_qrs = 0.0d0
       endif
    endif
! ============== 2014/09/19

    allocate(wkx(mmesh)); wkx = 0.d0
    allocate(wky(mmesh)); wky = 0.d0
    allocate(wkz(mmesh)); wkz = 0.d0
    allocate(rhpcr(mmesh));  rhpcr = 0.d0
    if(sw_positron == BULK .or. sw_positron == DEFECT) then
       allocate(rhcr(mmesh)); rhcr = 0.d0
    else if ( sw_add_corecharge_rspace == ON ) then
       allocate(rhcr(mmesh)); rhcr = 0.d0
    end if
    allocate(phir_ae(mmesh,nloc,ntau)); phir_ae = 0
    allocate(phir(mmesh,nloc,ntau)); phir = 0
    allocate(chir(mmesh,nloc,ntau)); chir = 0
    allocate(wos(mmesh)); wos = 0
    allocate(copsc(0:kord)); copsc = 0.d0
    allocate(s(ntau,ntau)); s = 0
    allocate(sinv(ntau,ntau)); sinv = 0
    !!$ print '(" mmesh, nloc,ntau,ntyp = ",4i5)', mmesh,nloc,ntau,ntyp
    if(.not.allocated(betar)) then
	allocate(betar(mmesh,nloc,ntau,ntyp)); betar = 0
    endif
    if(sw_use_add_proj == 1) then
       allocate(betar_add(mmesh,ntyp))
       betar_add = 0.d0
    end if
    if(ipripp >= 2) write(6,'(" sw_orb_popu = ",i5)') sw_orb_popu
    if(sw_orb_popu == 1) then
       allocate(phirt(mmesh,nloc,ntau,ntyp));phirt = 0
       if(ipripp >= 2)then
          write(6,'(" phirt is allocated")')
          write(6,'("  --  mmesh, nloc, ntau, ntyp = ",4i6)') mmesh, nloc, ntau, ntyp
       end if
! ---------------------------------<<
    end if

! ======= KT_add =========== 13.0U2
    if ( sw_modified_TFW_functional == ON .or. sw_use_contracted_psir == ON ) then
       allocate( psir_val(mmesh,nloc,ntau,ntyp)); psir_val = 0
    endif
! ========================== 13.0U2

    !if(sw_orb_popu == 1) then
    !   allocate(phirt(mmesh,nloc,ntau,ntyp)); phirt = 0
    !end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       allocate(paor(mmesh,nloc,ntau,ntyp)); paor = 0
    end if
    if(istress == 1) then
       allocate(wkz1(mmesh)); wkz1 = 0.d0
       allocate(wkz2(mmesh)); wkz2 = 0.d0
    endif

    if(flg_paw) then
        allocate(psir(mmesh,nloc,ntau));psir=0.d0
        allocate(ipppp(nltpw,nltpw,nltpw,nltpw,nloc*2,ntyp))
        ipppp=0
        allocate(vaeijlm_k(npppp));vaeijlm_k=0.d0
        allocate(vpsijlm_k(npppp));vpsijlm_k=0.d0
        allocate(vqijqlm_k(npppp));vqijqlm_k=0.d0
        allocate(vqijplpm_ks(npppp));vqijplpm_ks=0.d0
        allocate(vionaeij(nloc,ntau,ntau,ntyp));vionaeij=0.d0
        allocate(vionpsij(nloc,ntau,ntau,ntyp));vionpsij=0.d0
        allocate(vionpsqij(nloc,ntau,ntau,ntyp));vionpsqij=0.d0
        allocate(kin_ae_psij(nloc,ntau,ntau,ntyp));kin_ae_psij=0.d0
        allocate(vloc_scr_ps(mmesh));vloc_scr_ps=0.d0
        allocate(vloc_scr_ae(mmesh));vloc_scr_ae=0.d0
    end if
    if(firstcall.or.sw_optimize_lattice==ON.or.imdalg==P_CONTROL.or.imdalg==PT_CONTROL .or. &
     & driver == DRIVER_URAMP .or. driver == DRIVER_SC_DFT .or.  sw_stress_correction == ON) epc = 0.d0
    firstcall = .false.
  end subroutine m_PP_alloc_ps

  subroutine m_PP_init_epc()
    epc = 0.d0
  end subroutine m_PP_init_epc

#ifdef __EDA__
! -----  ascat starts modifying  -----
  subroutine m_PP_alloc_PP_per_atom_etc
    allocate(PP_per_atom(ntyp))
    allocate(epc_per_atom(ntyp));epc_per_atom = 0.d0
  end subroutine m_PP_alloc_PP_per_atom_etc

  subroutine m_PP_dealloc_PP_per_atom_etc
    deallocate(PP_per_atom)
    deallocate(epc_per_atom)
  end subroutine m_PP_dealloc_PP_per_atom_etc

  subroutine m_PP_set_PP_per_atom_etc_zero
    PP_per_atom = 0.d0
  end subroutine m_PP_set_PP_per_atom_etc_zero
! -----  ascat ceases modifying  -----
#endif

  subroutine m_PP_dealloc_ps
    deallocate(qij)
    deallocate(qvij)
    deallocate(vvv)
    deallocate(vxc)
    deallocate(exc)
    deallocate(vlocr)
    deallocate(vlocr2)
    deallocate(rhvr)
    if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) deallocate(radr)
    deallocate(qrsps)
    if(sw_berry_phase == ON .or. sw_wannier == ON .or. sw_fef == ON) deallocate(qrs)

! ==== KT_add === 2015/02/23
    if ( allocated(qrs) ) deallocate( qrs )
! =============== 2015/02/23

! ==== KT_add === 2014/09/19
    if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
       if ( allocated( kins_qrs ) ) deallocate( kins_qrs )
       if ( allocated( kina_qrs ) ) deallocate( kina_qrs )
    endif
! =============== 2014/09/19

    deallocate(wkx)
    deallocate(wky)
    deallocate(wkz)
!!$    if(.not.flg_paw) deallocate(rhpcr)
    deallocate(rhpcr)
!    if(sw_positron == BULK .or. sw_positron == DEFECT) deallocate(rhcr)
    if ( allocated(rhcr) ) deallocate(rhcr)

    deallocate(phir_ae)
    deallocate(phir)
    deallocate(chir)
    if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) deallocate(wos)
    deallocate(copsc)
    deallocate(s)
    deallocate(sinv)
    if(icond == INITIAL .or. icond == CONTINUATION) deallocate(betar)
    if(icond == INITIAL .or. icond == CONTINUATION) then
       if(sw_use_add_proj == 1) deallocate(betar_add)
    end if
    if ( icond == COORDINATE_CONTINUATION ) then
       if ( allocated(betar) ) deallocate( betar )
       if ( allocated(betar_add) ) deallocate( betar_add )
    endif
! ----- T. Yamasaki, 3 July 2008 ---
!!$    if(sw_orb_popu == 1) deallocate(phirt)
!!$    if(ipripp>=1 .and. sw_orb_popu == 1) write(6,'(" phirt is not allocated")')
! --------------------------------<<
    if(intzaj == by_pseudo_atomic_orbitals) deallocate(paor)
    if(istress == 1) then
       deallocate(wkz1)
       deallocate(wkz2)
    endif

    if(flg_paw) then
        deallocate(psir)
        deallocate(ipppp)
        deallocate(vaeijlm_k)
        deallocate(vpsijlm_k)
        deallocate(vqijqlm_k)
        deallocate(vqijplpm_ks)
        deallocate(vionaeij)
        deallocate(vionpsij)
        deallocate(vionpsqij)
        deallocate(kin_ae_psij)
        deallocate(vloc_scr_ps)
        deallocate(vloc_scr_ae)
    end if

  end subroutine m_PP_dealloc_ps

  subroutine m_PP_alloc_NLP
    if(.not.allocated(wos)) then
       iflag_wos = ON
       allocate(wos(mmesh))
    end if
    if(.not.allocated(radr)) then
       iflag_radr = ON
       allocate(radr(mmesh))
    end if
!!$    if(.not.allocated(betar)) then
!!$       iflag_betar = ON
!!$       allocate(betar(mmesh,nloc,ntau,ntyp))
!!$       write(6,'(" ! betar is allocated (m_PP_alloc_NLP)")')
!!$    end if
!!$    if(.not.allocated(ilmt)) then
!!$       iflag_ilmt = ON
!!$       allocate(ilmt(ntyp))
!!$       write(6,'(" ! ilmt is allocated (m_PP_alloc_NLP)")')
!!$    end if
  end subroutine m_PP_alloc_NLP

  subroutine m_PP_dealloc_NLP
    if(iflag_wos == ON .and. allocated(wos)) then
       if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) deallocate(wos)
       iflag_wos = OFF
    end if
    if(iflag_radr == ON .and. allocated(radr)) then
       if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) deallocate(radr)
       iflag_radr = OFF
    end if
!!$    if(iflag_betar == ON .and. allocated(betar)) then
!!$       deallocate(betar)
!!$       iflag_betar = OFF
!!$    end if
!!$    if(iflag_ilmt  == ON .and. allocated(ilmt)) then
!!$       deallocate(ilmt)
!!$       iflag_ilmt = OFF
!!$    end if
  end subroutine m_PP_dealloc_NLP

  subroutine alloc_betar()
    allocate(betar(mmesh,nloc,ntau,ntyp))
  end subroutine alloc_betar

  logical function m_PP_betar_calculated()
    m_PP_betar_calculated = .false.
    if(allocated(betar)) m_PP_betar_calculated = .true.
  end function m_PP_betar_calculated

  subroutine m_PP_alloc_ps_stress
    allocate(wkz1(mmesh)); wkz1 = 0.d0
    allocate(wkz2(mmesh)); wkz2 = 0.d0
  end subroutine m_PP_alloc_ps_stress

  subroutine m_PP_alloc_psc_qitg_rhpcg


    allocate(psc_l(ista_kngp:iend_kngp,ntyp)); psc_l = 0.d0
    allocate(qitg_l(ista_kngp:iend_kngp,nqitg)); qitg_l = 0.d0
    allocate(rhpcg_l(ista_kngp:iend_kngp,ntpcc)); rhpcg_l = 0.d0

    if(sw_berry_phase == ON) then
       allocate(qitg_BP(1,nqitg)); qitg_BP = 0.d0
    end if
    if(sw_positron == BULK .or. sw_positron == DEFECT) then

       allocate(rhcg_l(ista_kngp:iend_kngp,ntyp)); rhcg_l = 0.d0
       allocate(rhceg_l(ista_kngp:iend_kngp,ntyp)); rhceg_l = 0.d0
       allocate(rhchg_l(ista_kngp:iend_kngp,ntyp)); rhchg_l = 0.d0
    end if

    if ( sw_add_corecharge_rspace == ON .and. eval_corecharge_on_Gspace == ON ) then
       if ( .not.allocated(rhcg_l) ) then
          allocate(rhcg_l(ista_kngp:iend_kngp,ntyp)); rhcg_l = 0.d0
       endif
    endif
    if ( xctype == "tb09" ) then
       if ( .not.allocated(rhcg_l) ) then
          allocate(rhcg_l(ista_kngp:iend_kngp,ntyp)); rhcg_l = 0.d0
       endif
    endif

    if(initial_chg == from_PSEUDOPOTENTIAL_FILE.or.sw_charge_predictor==ON) then

       allocate(rhvg_l(ista_kngp:iend_kngp,ntyp)); rhvg_l = 0.d0
    end if
    if(istress==ON) then


       allocate(psc_diff_l(ista_kngp:iend_kngp,ntyp))
       allocate(qitg_diff_l(ista_kngp:iend_kngp,nqitg))
       allocate(rhpcg_diff_l(ista_kngp:iend_kngp,ntpcc))
       psc_diff_l = 0.d0; qitg_diff_l = 0.d0; rhpcg_diff_l = 0.d0
    endif

! === KT_add ==== 2014/09/19
    if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
       if ( use_asymm_ekin_density ) then
          allocate( kina_qitg_l(ista_kngp:iend_kngp,nqitg) ); kina_qitg_l = 0.d0
       endif
       if ( use_symm_ekin_density ) then
          allocate( kins_qitg_l(ista_kngp:iend_kngp,nqitg) ); kins_qitg_l = 0.d0
       endif
    endif
! =============== 2014/09/19

  end subroutine m_PP_alloc_psc_qitg_rhpcg

  subroutine m_PP_dealloc_psc_qitg_rhpcg
    if(allocated(psc_l)) deallocate(psc_l)
    if(allocated(qitg_l)) deallocate(qitg_l)
    if(allocated(rhpcg_l)) deallocate(rhpcg_l)
    if(sw_berry_phase==ON .and. allocated(qitg_BP)) deallocate(qitg_BP)
    if(sw_positron == BULK .or. sw_positron == DEFECT) deallocate(rhcg_l)
    if(sw_positron == BULK .or. sw_positron == DEFECT) deallocate(rhceg_l)
    if(sw_positron == BULK .or. sw_positron == DEFECT) deallocate(rhchg_l)

    if ( allocated(rhcg_l) ) deallocate( rhcg_l )

    if(initial_chg == from_PSEUDOPOTENTIAL_FILE.or.sw_charge_predictor==ON) then
       if(allocated(rhvg_l)) deallocate(rhvg_l)
    endif
    if(istress==1) then
       if(allocated(psc_diff_l)) deallocate(psc_diff_l)
       if(allocated(qitg_diff_l)) deallocate(qitg_diff_l)
       if(allocated(rhpcg_diff_l)) deallocate(rhpcg_diff_l)
    endif
! ==== KT_add === 2014/09/19
    if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
       if ( allocated( kins_qitg_l ) ) deallocate( kins_qitg_l )
       if ( allocated( kina_qitg_l ) ) deallocate( kina_qitg_l )
    endif
! =============== 2014/09/19
  end subroutine m_PP_dealloc_psc_qitg_rhpcg

  subroutine m_PP_alloc_qitg_wan(nwght)
    integer, intent(in) :: nwght

    nwght_wan = nwght
    ! debug
      write(6,*) 'm_PP_alloc_qitg_wan: nwght_wan=',nwght_wan
    ! end debug
    allocate(qitg_wan(nqitg,nwght_wan)); qitg_wan = 0.d0
    allocate(dk_wan(nwght_wan)); dk_wan = 0.d0
  end subroutine m_PP_alloc_qitg_wan

  subroutine m_PP_set_dk_wan(mp_index,ghat,rltv)
    integer, intent(in) :: mp_index(3)
    integer, intent(in) :: ghat(3,nwght_wan)
    real(kind=DP), intent(in) :: rltv(3,3)

    integer :: i
    real(kind=DP) :: b(3,3),g(3)
    do i=1,3
       b(1:3,i) = rltv(1:3,i)/mp_index(i)
    end do
    do i=1,nwght_wan
       g = matmul(b,ghat(:,i))
       dk_wan(i) = sqrt(dot_product(g,g))
    ! debug
      write(6,*) 'm_PP_set_dk_wan: dk_wan(i)=',dk_wan(i)
    ! end debug
    end do
  end subroutine m_PP_set_dk_wan

  subroutine m_PP_alloc_qitg_fef(numef)
    integer, intent(in) :: numef

    numef_fef = numef
    ! debug
      write(6,*) 'm_PP_alloc_qitg_fef: numef_fef=',numef_fef
    ! end debug
    allocate(qitg_fef(nqitg,numef_fef)); qitg_fef = 0.d0
    allocate(dk_fef(numef_fef)); dk_fef = 0.d0
  end subroutine m_PP_alloc_qitg_fef

  subroutine m_PP_set_dk_fef(mp_index,rltv,elec_id)
    integer, intent(in) :: mp_index(3),elec_id(3)
    real(kind=DP), intent(in) :: rltv(3,3)

    integer :: ig,id
    real(kind=DP) :: g(3)

    do ig=1,3
       id = elec_id(ig)
       if(id == 0) cycle
       g = rltv(1:3,ig)/mp_index(ig)
       dk_fef(id) = sqrt(dot_product(g,g))
    ! debug
      write(6,*) 'm_PP_set_dk_fef: dk_fef(id)=',dk_fef(id)
    ! end debug
    end do
  end subroutine m_PP_set_dk_fef



  subroutine sum_of_vlocr(it,nfout)
    integer, intent(in) :: it,nfout
    real(kind=DP) :: s
    integer       :: i
    s = 0.d0
    do i = 1, nmesh(it)
       s = s + wos(i)*vlocr(i)*radr(i)**2
    end do
    if(ipripp >= 1) then
       write(nfout,'(" !PP it = ",i4," <<sum_of_vlocr>>")') it
       write(nfout,'(" !PP Int_{i=1}^{nmesh} wos(i)vlocr(i)*radr(i)**2 = " &
            & ,d20.8)') s
    end if
  end subroutine sum_of_vlocr

  subroutine sum_of_vvv_and_vxc(it,nfout)
    integer, intent(in) :: it,nfout
    real(kind=DP)   :: s1, s2
    integer         :: i
    s1 = 0.d0; s2 = 0.d0
    do i = 1, nmesh(it)
       s1 = s1 + wos(i)*vvv(i)*radr(i)
       s2 = s2 + wos(i)*vxc(i)*radr(i)
    end do
    if(ipripp>=2) then
       write(nfout,'(" !PP a3 = ", d20.8)') s1
       write(nfout,'(" !PP a2 = ", d20.8)') s2
    end if
  end subroutine sum_of_vvv_and_vxc

  integer function m_PP_check_file_format(nfp,nfout)
    integer, intent(in) :: nfp, nfout
    integer :: fileformat

    if(mype == 0) then
       fileformat = 0
       if(findtag(nfp,nfout,tag_fileformat,str_obj) /= 0) then
          call toupper(str_obj)
          if(index(str_obj,tag_CIAOPP) /= 0) then
             fileformat = pp_CIAOPP
          else if(index(str_obj,tag_GNCPP2) /= 0) then
             fileformat = pp_GNCPP2
          end if
       end if
    end if
    if(npes > 1) &
         & call mpi_bcast(fileformat,1,mpi_integer,0,MPI_CommGroup,ierr)

    if(ipripp>=1) then
       write(nfout,'(" fileformat = ",i8)') fileformat
       write(nfout,'(" ( CIAOPP = ",i8," GNCPP2 = ",i8)') pp_CIAOPP, pp_GNCPP2
    end if

    m_PP_check_file_format = fileformat
  end function m_PP_check_file_format

  subroutine toupper(str)
    character(len=*), intent(inout) :: str
    character(len=1) :: c
    integer :: i, ndif

    ndif = ichar('A') - ichar('a')
    do i = 1, len_trim(str)
       c = str(i:i)
       if(c == '!') exit
       if(ichar(c).ge.ichar('a').and.ichar(c).le.ichar('z')) str(i:i)=char(ichar(c) + ndif)
    end do
  end subroutine toupper

  integer function findtag(nfp,nfout,tagstring,str_obj)
    integer, intent(in)          :: nfp,nfout
    character(len=*), intent(in) :: tagstring
    character(len=*), intent(out):: str_obj

    logical :: comment_statement, tag_is_found
    character(len=1) :: c
    integer :: ic, i, lenc, ic2, ndif

    ndif = ichar('A') - ichar('a')

    findtag = 0

! === DEBUG by tkato 2011/06/29 ================================================
    tag_is_found = .false.
! ==============================================================================
    comment_statement = .true.
    rewind(nfp)
    do while(comment_statement)
       read(nfp,'(a80)') str
       str = adjustl(str)
       lenc = len_trim(str)
       if(lenc == 0) then
          if(ipripp >= 2) &
               & write(nfout,'(" len_trim of ",a,"=",i8)') trim(str),len_trim(str)
          tag_is_found = .false.
          comment_statement = .true.
       else
          if(str(1:1) == '#'.or.str(1:1) == '$' .or. str(1:1) == '!' &
               & .or. str(1:1) == '%' .or. str(1:1) == '*') then
             comment_statement = .true.
             ! to upper
             do i = 2, lenc
                c = str(i:i)
                if(c=='=') exit
                if(ichar(c).ge.ichar('a').and.ichar(c).le.ichar('z')) str(i:i)=char(ichar(c) + ndif)
             end do
             ic = index(str,tagstring)
             if(ic == 0) then
                tag_is_found = .false.
             else
                tag_is_found = .true.
                if(ipripp >= 2) write(nfout,'(" FileFormat is found :",a)') str
             end if
             if(tag_is_found) then
                ic2 = index(str(ic+len("FILEFORMAT"):lenc),'=')
                if(ic2 == 0) then
                   tag_is_found = .false.
                else
                   ic2 = ic2 + ic + len('FILEFORMAT')
                end if
             end if
             if(tag_is_found) then
                str_obj = str(ic2:lenc)
                exit
             end if
          else
             comment_statement = .false.
          end if
       end if
       if(ipripp >= 2 .and. comment_statement) write(nfout,'(a)') str
    end do
    if(tag_is_found) findtag = 1
  end function findtag

  subroutine m_PP_check_gncpp_type(nfp,nfout,is_gncpp)
    ! Original code has been written by K.Mae and/or M.Okamoto, and has
    ! been contained in <PseudoPotential_Construction.F90>.
    ! That subroutine is move here, and mpi_bcast is applied to the
    ! variable of 'is_gncpp' by T. Yamasaki
    !                                   18th Sep. 2003
    integer, intent(in) ::    nfp, nfout
    integer, intent(out) ::   is_gncpp
    integer :: fn_number_of_words, igncpp
    logical :: comment_statement
    integer :: iloc_tmp,itpcc_tmp, igncpp_tmp
    real(kind=DP) :: fval_tmp, natomn_tmp
    integer :: pptype
    character(len=len("PAW (=GNCPP2 with AE wavefunctions)")) :: pptype_char

    is_gncpp = pp_GNCPP1 ! = 1

    if(mype == 0) then
       pptype = 0
       pptype_char = ""
       if(findtag(nfp,nfout,tag_potentialtype,str_obj) /= 0) then
          call toupper(str_obj)
          if(index(str_obj,tag_PAW) /= 0)then
             pptype = pp_PAW
          else if(index(str_obj,tag_ultrasoft) /= 0) then
             pptype = pp_GNCPP2
          end if
       end if
       if(pptype == pp_PAW) then
          is_gncpp = pp_PAW
          pptype_char = "PAW (=GNCPP2 with AE wavefunctions)"
          if(ipripp>=1 .and. .not.ppprinted) then
             write(nfout,'(" !PP PP type --> ",a," , is_gncpp = ",i2)') trim(pptype_char),is_gncpp
!!$          if(ipripp>=1 .and. .not.ppprinted) &
!!$               & write(nfout,'(" !PP PP type PAW (= GNCPP2 with AE wavefunctions)")')
          end if
       else
          comment_statement = .true.
          rewind(nfp)
          if(ipripp >= 1 .and. .not.ppprinted) write(nfout,'(" !PP CHECKING POTENTIAL FILE ",i5)') nfp
          do while(comment_statement)
             read(nfp,'(a80)') str
             if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
                  & .or. str(1:1) == '%' .or. str(1:1) == '*') then
                if(ipripp >= 1 .and. .not.ppprinted) write(nfout,'(a80)') str
             else if(len(trim(str)) == 0) then
                if(ipripp >= 1 .and. .not.ppprinted) write(nfout,'(a80)') str
             else
                comment_statement = .false.
             endif
          enddo

          select case (fn_number_of_words(str))
!          case (9)
          case (4,9)
             is_gncpp = pp_GNCPP1 ! = 1
             pptype_char = "GNCPP1"
!!$             if(ipripp>=1 .and. .not.ppprinted) write(nfout,'(" !PP PP type --> GNCPP1")')
             read(str,*)    natomn_tmp,fval_tmp,iloc_tmp,itpcc_tmp
!!$             if(ipripp>=1 .and. .not.ppprinted) write(nfout,'(" !PP ",2d20.8,2i8)') natomn_tmp,fval_tmp,iloc_tmp,itpcc_tmp
          case (10, 11)
             read(str,*)    natomn_tmp,fval_tmp,iloc_tmp,itpcc_tmp,igncpp_tmp
!!$             if(ipripp>=1 .and. .not.ppprinted) write(nfout,'(" !PP ",2d20.8,3i8)') &
!!$                  & natomn_tmp,fval_tmp,iloc_tmp,itpcc_tmp,igncpp_tmp
             select case (igncpp_tmp)
             case (1)
                is_gncpp = pp_GNCPP1 ! = 1
                pptype_char = "GNCPP1"
!!$                if(ipripp>=1 .and. .not.ppprinted) write(nfout,'(" !PP PP type --> GNCPP1")')
             case (2)
                is_gncpp = pp_GNCPP2 ! = 2
                pptype_char = "GNCPP2"
!!$                if(ipripp>=1 .and. .not.ppprinted) write(nfout,'(" !PP PP type --> GNCPP2")')
             case (-2)
                is_gncpp = pp_PAW ! = -2 = pp_GNCPP2_with_AE_WF
                pptype_char = "GNCPP2 with AE wavefunctions"
!!$                if(ipripp>=1 .and. .not.ppprinted) &
!!$                     & write(nfout,'(" !PP PP type --> GNCPP2 with AE wavefunctions")')
             case default
                if(ipripp>=1 .and. .not.ppprinted) then
                   write(nfout,'(" !PP ### ERROR ### gncpp type is wrong")')
                   write(nfout,'(" !PP   igncpp ... ",i8)') igncpp
                end if
                call phase_error_with_msg(nfout,'### ERROR ### gncpp type is wrong',__LINE__,__FILE__)
             end select
          case default
             if(ipripp >=1 .and. .not.ppprinted) then
                write(nfout,'(" !PP ### ERROR ### number of args of the first line of pp_data")')
                write(nfout,'(" !PP   number of args ...",i8)') fn_number_of_words(str)
             end if
             call phase_error_with_msg(nfout,'### ERROR ### number of args of the first line of pp_data', &
             __LINE__,__FILE__)
          end select
          rewind(nfp)
          if(ipripp>=1 .and. .not.ppprinted) then
             write(nfout,'(" !PP PP type --> ",a," , is_gncpp = ",i2)') trim(pptype_char),is_gncpp
             write(nfout,'(" !PP natomn, fval, iloc, itpcc = ",2f14.6, 2i3)') &
                  & natomn_tmp,fval_tmp,iloc_tmp,itpcc_tmp
          end if
       end if
    end if
    if(npes > 1) then
       if(mype == 0) igncpp = is_gncpp
       call mpi_bcast(igncpp,1,mpi_integer,0,MPI_CommGroup,ierr)
       is_gncpp = igncpp
    end if
    if(ipripp >=2 .and. .not.ppprinted) write(nfout,'(" !PP     is_gncpp = ",i8)') is_gncpp
  end subroutine m_PP_check_gncpp_type

! =============================== added by K. Tagami ======================= 11.0
  subroutine m_PP_check_spinorbit( nfp,nfout, has_spinorbit )
    integer, intent(in) ::    nfp, nfout
    logical, intent(out) ::   has_spinorbit
    integer :: itmp, i_spinorbit

!!$    has_spinorbit =  OFF
    i_spinorbit =  OFF

    if (mype == 0) then
       if ( findtag( nfp,nfout,tag_pot_has_spinorbit,str_obj) /= 0) then
          call toupper(str_obj)
          if ( index(str_obj,tag_ON ) /= 0 ) then
!!$             has_spinorbit = ON
             i_spinorbit = ON
          endif
       end if
    endif
    if (npes > 1) then
!!$       if(mype == 0) itmp = has_spinorbit
       if(mype == 0) itmp = i_spinorbit
       call mpi_bcast( itmp,1,mpi_integer,0,MPI_CommGroup,ierr)
!!$       has_spinorbit = itmp
       i_spinorbit = itmp
    end if

    if (.not.ppprinted.and.((ipripp >=1 .and. i_spinorbit /=OFF).or.ipripp>=2)) &
         & write(nfout,'(" !PP  has_spinorbit = ",i4)') i_spinorbit
    if(i_spinorbit == ON) then
       has_spinorbit = .true.
    else
       has_spinorbit = .false.
    end if

  end  subroutine m_PP_check_spinorbit

! ========================================================================== 11.0


! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q ORG_Parallel
! ==================================================================================================
  subroutine rd_core_charge_only( nfp, it, nfout, mode, paramset, core_charge_found )
    use m_Files, only : m_Files_open_core_charge_dens, m_Files_close_core_charge_dens,&
         &              nf_core_charge_radial
    integer, intent(in)  :: it, nfp, nfout, mode
    integer, intent(out) :: core_charge_found
    logical, intent(in)       :: paramset

    integer :: lun, ifound, i
    real(kind=DP) :: dummy

    if ( mode == 0 ) then     ! from F_POT
       lun = nfp
    else                      ! from F_CORE_CHARGE
       call m_Files_open_core_charge_dens(it)
       lun = nf_core_charge_radial(it)
    endif

    ifound = 0
    if ( mype == 0 ) then
       Do while (.true.)
          read(lun,'(a80)',end=10) str
          ifound = index( str,'CORE-CHARGE' )
          if ( ifound /= 0 ) goto 10
          ifound = index( str,'core-charge' )
          if ( ifound /= 0 ) goto 10
       End do
10     continue

       if ( ifound /= 0 ) then
          read(lun,*);
          if(paramset) then
             read(lun,*) (dummy,i=1,nmesh(it))
          else
             read(lun,*) rhcr(1:nmesh(it))
          endif
       else
          if ( ifound == 0 .and. iatomn(it) /= ival(it) ) then
             write(nfout,*)
             if ( sw_add_corecharge_rspace == ON ) then
                write(nfout,'(A)') " WARNING!!! "
                write(nfout,'(A)') " Although sw_add_corecharge_rspace is set ON,"
             endif
             if ( mode==0 ) then
                write(nfout,'(A,I2,A)') &
                     &    " CORE-CHARGE is not found in F_POT(",it,")"
             else
                write(nfout,'(A,I2,A)') &
                     &    " CORE-CHARGE is not found in F_CORE_CHARGE(",it,")"
             endif
             write(nfout,*)
          endif
       endif
    endif

    core_charge_found = ifound
    if (npes > 1 ) then
       call mpi_bcast(core_charge_found,1,mpi_integer,0,MPI_CommGroup,ierr)
    end if
    if (npes > 1 .and. .not.paramset ) then
       call mpi_bcast(rhcr,mmesh,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(total_core_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    end if

    if ( mode==1 ) call m_Files_close_core_charge_dens(it)

  end subroutine rd_core_charge_only

  subroutine rd_core_charge_then_ft(nfp,it,nfout,gr_l,paramset,is_gncpp,&
       &                            core_charge_found )
    integer, intent(in) :: nfp, it, nfout, is_gncpp
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical, intent(in)       :: paramset
    integer, intent(out)      :: core_charge_found
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      real(kind=DP),             dimension(mmesh) :: wkx, wky
#endif

    integer :: i, n, ifound
    logical :: tf
    real(kind=DP), allocatable, dimension(:) :: enhance
    real(kind=DP) :: gabs,fac,fac2,dummy,valence_charge

    core_charge_found = 0
    if(mype == 0) then
#if 1
       if ( is_gncpp == pp_PAW ) then
          rewind(nfp)
          ifound = 0;  tf = .false.

          Do while (.true.)
             read(nfp,'(a80)',end=10) str
             ifound = index( str,'CORE-CHARGE' )
             if ( ifound /= 0 ) goto 10
          End do
10        continue
          if ( ifound /= 0 ) tf = .true.
       else
          read(nfp,'(a80)',end=1001) str
          call strncmp0(trim(str),tag_core_charge,tf)
       endif
#else
       read(nfp,'(a80)',end=1001) str
       call strncmp0(trim(str),tag_core_charge,tf)
#endif

       if(tf) then
          core_charge_found = 1
          read(nfp,*) total_core_charge
          if(paramset) then
             read(nfp,*) (dummy,i=1,nmesh(it))
          else
             read(nfp,*) (rhcr(i),i=1,nmesh(it))
          end if
       else
          backspace nfp
       end if
1001   continue
    end if
    if(ipripp >= 2) then
       if(core_charge_found == 1) then
          write(nfout,'(" !PP core_charge_found = YES <<rd_core_charge_then_ft>>")')
       else
          write(nfout,'(" !PP core_charge_found = NO  <<rd_core_charge_then_ft>>")')
       end if
    end if
    if(npes > 1 ) then
       call mpi_bcast(core_charge_found,1,mpi_integer,0,MPI_CommGroup,ierr)
    end if

    if(core_charge_found == 1) then
       if(npes > 1 .and. .not.paramset ) then
          call mpi_bcast(rhcr,mmesh,mpi_double_precision,0,MPI_CommGroup,ierr)
          call mpi_bcast(total_core_charge,1,mpi_double_precision,0,MPI_CommGroup,ierr)
       end if
       if(ipripp >=2 ) write(nfout,'(" !PP total_core_charge = ", f8.4)') total_core_charge
       if(sw_positron == ON ) then
          allocate(enhance(mmesh)); enhance = 0.d0
       endif

       if(.not.paramset) then

          valence_charge=0.d0
          do n=1,ntyp
            valence_charge=valence_charge+ival(n)*natm
          end do

          if (sw_positron ==ON) then
             do n = 1, nmesh(it)
                wkz(n)=rhcr(n)/radr(n)/radr(n)/PAI4
!             fac = wkz(n)
                fac = wkz(n)+valence_charge/univol

                call enhance_0(fac,fac2)
                if(sw_epsilon_ele == ON) then
                   call enhance_01(fac,fac2,epsilon_ele)
                end if
                enhance(n) = fac2
             end do
             call epcor_00(1.d-15,wkz,nmesh(it))
!          wkz=wkz*radr(n)*4.d0*3.1415926d0
             wkz=0.d0
          endif

! NEC no check
!!CDIR PARALLEL DO PRIVATE ( gabs, wkx, wky, n, fac )
          do i = ista_kngp, iend_kngp
             gabs = gr_l(i)
             rhcg_l(i,it) = 0.d0
             if ( sw_positron == ON ) then
                rhceg_l(i,it) = 0.d0; rhchg_l(i,it) = 0.d0
             endif
             wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))

             call dsjnv(0,nmesh(it),wkx,wky)
             do n = 1, nmesh(it)
                fac = wos(n)*rhcr(n)/univol
                rhcg_l(i,it) = rhcg_l(i,it) + fac*wky(n)
                if ( sw_positron == ON ) then
                   rhceg_l(i,it) = rhceg_l(i,it) + wos(n)*wkz(n)*wky(n)
                   rhchg_l(i,it) = rhchg_l(i,it) + fac*wky(n)*enhance(n)
                endif
             end do
          end do
       end if
       if(sw_positron == ON) deallocate(enhance)
    end if
  end subroutine rd_core_charge_then_ft
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q END ORG_Parallel
! ==================================================================================================

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q ORG_Parallel
! ==================================================================================================
  subroutine ft_valence_charge(it,nfout,gr_l,ngshell,ngshell_range)
    integer, intent(in) :: it, nfout
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    integer, intent(in) :: ngshell
    integer, intent(in), dimension(ngshell,2) :: ngshell_range

    integer       :: n, igs
    real(kind=DP) :: gabs,fac
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
    real(kind=DP),             dimension(mmesh) :: wkx, wky
#endif
    real(kind=DP) :: rhvg_shell

! NEC no check
!!CDIR PARALLEL DO PRIVATE ( gabs, wkx, wky, n, fac, rhvg_shell )
    do igs = 1, ngshell
       gabs = gr_l(ngshell_range(igs,1))
       wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
       call dsjnv(0,nmesh(it),wkx,wky)
       rhvg_shell = 0.d0
       do n = 1, nmesh(it)
          fac = wos(n)*rhvr(n)
          rhvg_shell = rhvg_shell + fac*wky(n)
       enddo
       rhvg_shell = rhvg_shell/univol
       do n = ngshell_range(igs,1), ngshell_range(igs,2)
          rhvg_l(n,it) = rhvg_shell
       end do
    enddo

!!$    do i = ista_kngp, iend_kngp
!!$       gabs = gr_l(i)
!!$       rhvg_l(i,it) = 0.d0
!!$       wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
!!$       call dsjnv(0,nmesh(it),wkx,wky)
!!$       do n = 1, nmesh(it)
!!$          fac = wos(n)*rhvr(n)/univol
!!$          rhvg_l(i,it) = rhvg_l(i,it) + fac*wky(n)
!!$       end do
!!$    end do
  end subroutine ft_valence_charge
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q_END ORG_Parallel
! ==================================================================================================

  subroutine where_does_WFrhvr_damp(nfout,it,mesh_t)
    integer, intent(in) :: nfout,it
    integer, intent(out) :: mesh_t

    integer     :: i
    real(kind=DP), parameter :: CRDAMP = 1.d0
    real(kind=DP), parameter :: CRDIST = 10.d0
    do i = 10, nmesh(it)-1
       if(rhvr(i)-rhvr(i+1) > CRDAMP .and. radr(i) < CRDIST) then
          mesh_t = i
          if(ipripp >= 2) write(nfout,'(" !PP LMTO pot. r_ws=",i5,f12.6)') i, radr(i)
          return
       end if
    enddo
    mesh_t = nmesh(it)
  end subroutine where_does_WFrhvr_damp

  subroutine make_index_arrays(nfout,it,paramset)
    integer, intent(in) :: nfout,it
    logical, intent(in) :: paramset
    integer mm, il1, im1, il2, im2, ilmt1, ilmt2, ip

    mm = 0
    do ilmt1 = 1, ilmt(it)
       mm = mm + 1
       if(.not.paramset) then
          index_lmt1_lmt2(mm,it,1) = ilmt1
          index_lmt1_lmt2(mm,it,2) = ilmt1
          w_non0_lmtxlmt(mm,it)    = 1
       endif
    enddo

    do ilmt1 = 1, ilmt(it)
       il1 =  ltp(ilmt1,it)
       im1 =  mtp(ilmt1,it)
       do ilmt2 = ilmt1+1, ilmt(it)
          il2 = ltp(ilmt2,it)
          im2 = mtp(ilmt2,it)
          if(il1 == il2 .and. im1 == im2) then
             mm = mm + 1
             if(.not.paramset) then
                index_lmt1_lmt2(mm,it,1) = ilmt1
                index_lmt1_lmt2(mm,it,2) = ilmt2
                w_non0_lmtxlmt(mm,it)    = 2
             endif
          endif
       enddo
    enddo
    n_non0_lmtxlmt(it) = mm

    if(ipripp >= 2) then
       if(paramset) then
          write(nfout,'(" !PP #non0_elements = ",i8)') mm
       else
          write(nfout,'(" !PP --- index_lmt1_lmt2 --")')
          write(nfout,'(" !PP ",10(" (",2i3,")"))') &
               &(index_lmt1_lmt2(ip,it,1),index_lmt1_lmt2(ip,it,2),&
               & ip=1,n_non0_lmtxlmt(it))
       end if
    end if
  end subroutine make_index_arrays

!/////////////////////////////////////////////////////////////////////

!#####################################################################
! GNCPP2 type pseudopotential (M.Okamoto, August, 2003)

!=====================================================================
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q ORG_Parallel
! ==================================================================================================
  subroutine m_PP_vanderbilt_type(nfp,it,nfout,gr_l,ngshell,ngshell_range &
       & ,paramset,is_gncpp,return_after_rd_itau_etc)
!=====================================================================
! This subroutine was modified from m_PP_vanderbilt_type,
! in order to read GNCPP2 type pseudopotential data.
! (M.Okamoto, August 2003)
!
! This subroutine (m_PP_vanderbilt_type_gncp2) was again merged with
! the original subroutine of m_PP_vanderbilt_type and was renamed.
!                       T. Yamasaki, October 2006
! This subroutine was merged with m_PP_vanderbilt_type_gncpm2, which
! is for PAW potential.
!                       T. Yamaskai, May 2010
!---------------------------------------------------------------------
    integer,intent(in) :: nfp, it, nfout,ngshell
    integer,intent(in),dimension(ngshell,2) :: ngshell_range
    real(kind=DP),intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical,intent(in) :: paramset
    integer,intent(in) :: is_gncpp
    logical,intent(in), optional :: return_after_rd_itau_etc
    logical :: ret
    logical            :: vall_is_detected
    integer            :: core_charge_found
    integer            :: mesh_t, iprippex
    integer            :: id_sname = -1

    ret = .false.
    if(present(return_after_rd_itau_etc)) ret = return_after_rd_itau_etc
    iprippex = ipripp
!!$    if(.not.paramset) iprippex = iprippex - 1
    if(paramset) iprippex = iprippex - 1

                                                  __TIMER_SUB_START(1222)
    call tstatc0_begin('m_PP_vanderbilt_type ',id_sname,1)

    if(it == 1) mmt = 0
    if(it == 1) mmpppp = 0
    if(is_gncpp == pp_PAW) then
       if(paramset) then
          if(PAW_switch == ON) then
             ipaw(it) = 1
             flg_paw = .true.
             if(.not.m_CtrlP_explicit_cmix()) then
                sw_mix_charge_hardpart = ON
                if(printable) write(nfout,'(a)') ' !** REMARK : sw_mix_charge_hardpart was set to ON'
             endif
! === KT_add ===== 2014/12/29
#if 0
          else if ( noncol ) then
             if(.not.m_CtrlP_explicit_cmix()) then
                sw_mix_charge_hardpart = ON
                if(printable) write(nfout,'(a)') ' !** REMARK : sw_mix_charge_hardpart is forced to be ON'
             endif
#endif
! ================ 2014/12/29
          end if
       end if
    end if

    call mpi_barrier(MPI_CommGroup,ierr)

! ===============================- added by K. Tagami ================== 11.0
    if ( pot_has_soc(it) ) then
       call case_when_pot_has_soc
    else
       call case_when_pot_has_nosoc
    endif

! ===============================- added by K. Tagami ================== 11.0
    if ( .not. paramset ) then
!!       if ( .not. flg_paw ) then
          if ( .not. pot_has_soc(it) ) then
             if ( SpinOrbit_Mode == ZeffApprox ) then
                call  make_SOC_strength_Zeff_nonpaw
             endif
          endif
!!       endif
    endif
! ====================================================================== 11.0

    call m_IS_set_ival(ntyp,ival)

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1222)

!....................................................................
  contains

! ============================== added by K. Tagami ====================== 11.0
    subroutine case_when_pot_has_nosoc

      call rd_Vloc_psWF_Q_then_chir_qitg      !-(contained here) (potrvd)
      !        -> ps_xctype, etc.

      call make_index_arrays_nlmt2l_m_t       ! (tlmst1)
      if(sw_orb_popu == ON) then
         call make_index_arrays_nlmt_phi2lmt  ! (tlmst1_phi)
      end if

      if(is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) then
         if(sw_use_add_proj == ON) then
            call make_index_arrays_nlmt_add2lmt     ! (tlmst1_add)
         end if
         if(intzaj == by_pseudo_atomic_orbitals) then
            call make_index_arrays_nlmt_pao2lmt     ! (tlmst1_pao)
         end if
      end if

! =========================== added by K. Tagami ======================= 5.0
!!$    if(is_gncpp == pp_GNCPP2.and. sw_mix_charge_hardpart)then
!!$    if(is_gncpp == pp_GNCPP2.and. sw_mix_charge_hardpart == ON)then
      if((is_gncpp == pp_GNCPP2 .or. is_gncpp==pp_PAW .and. PAW_switch == OFF)&
           & .and. ( sw_mix_charge_hardpart == ON .or. af == 1 .or. &
           &         sw_hybrid_functional == ON ) )then

         if(.not.associated(ia2ia_symmtry_op)) then
            allocate(ia2ia_symmtry_op(natm,nopr+af))
            ia2ia_symmtry_op=0
            call set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia_symmtry_op)
         end if
         if(.not.associated(ia2ia_symmtry_op_inv)) then
            allocate(ia2ia_symmtry_op_inv(natm,nopr+af))
            ia2ia_symmtry_op_inv=0
            call set_ia2ia_symmtry_op_inv &
                 (natm,nopr,ia2ia_symmtry_op,ia2ia_symmtry_op_inv)
         end if
!!      call make_ia2ia_symmtry_op_etc
      endif
! ===================================================================== 5.0

      if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
         call make_index_arrays_nlt2l_t
         call make_ia2ia_symmtry_op_etc
!!$       if(iprippex>=2) then
!!$          write(nfout,'(" -- ia2ia_symmtry_op --")')
!!$          do iop = 1, nopr+af
!!$             write(nfout,'(" iop = ",i3)') iop
!!$             write(nfout,'(10i6)') ia2ia_symmtry_op(1:natm,iop)
!!$          end do
!!$          write(nfout,'(" -- ia2ia_symmtry_op_inv --")')
!!$          do iop = 1, nopr+af
!!$             write(nfout,'(" iop = ",i3)') iop
!!$             write(nfout,'(10i6)') ia2ia_symmtry_op_inv(1:natm,iop)
!!$          end do
!!$       end if

         call cnstrct_of_PiPjPlPm_k_Qijk_etc
         if(.not.flg_symmtry) then
            call cnstrct_of_CijkClmkVVVVijlm_k
         else
            call cnstrct_of_CijkClmnVVVVijlm_kn
         end if
      else if ( use_metagga ) then
         call make_index_arrays_nlt2l_t                  ! test !!!
      end if

      if(.not.paramset) then
         call cnstrct_of_localWF_Dion_and_q2    ! betavd
         ! -> betar, dion, q
         call where_does_WFrhvr_damp(nfout,it,mesh_t) ! -> mesh_t
         call coulomb_potential_in_Rspace(5)    ! -> vvv()   (ptclmb)
         call atomic_xc_potential_in_Rspace     ! -> vxc()   (ptxc)
         ! jfp = 0 (non-spin system for atomic charge) y
         call vlocr_plus_hartree_xc2 ! vlocr = vlocr2, vlocr, vvv, vxc
#ifdef _POT_SMOOTHING_
         call smoothing_vlocr()                 ! vlocr smoothing
#endif
         if(is_gncpp == pp_PAW .and. ipaw(it)==ON) then
            call cnstrct_of_VHij_VHpsij_Kinij
            call cnstrct_of_dion_kin_ion
            call init_of_dion_paw
         end if
!!$       else
         if(initial_chg == from_PSEUDOPOTENTIAL_FILE.or.sw_charge_predictor==ON) &
              & call ft_valence_charge(it,nfout,gr_l,ngshell,ngshell_range) ! rhvr -> rhvg_l
        if(sw_positron /= OFF) then
            call rd_core_charge_then_ft(nfp,it,nfout,gr_l,paramset,is_gncpp, &
                 &                      core_charge_found ) ! rhcr -> rhcg_l
         else
            if ( sw_add_corecharge_rspace == ON ) then
               if ( eval_corecharge_on_Gspace == ON ) then
                  call rd_core_charge_then_ft(nfp,it,nfout,gr_l,paramset,is_gncpp,&
                       &                      core_charge_found )
               else
                  if ( is_gncpp /= pp_PAW ) then
                     call rd_core_charge_only( nfp, it, nfout, 0, paramset, &
                          &                    core_charge_found )
                     if ( core_charge_found == 1 ) then
                        rhcorpw(:,it) = rhcr(:)
                     else if ( sw_read_corecharge_extra_file == ON ) then
                        call rd_core_charge_only( nfp, it, nfout, 1, paramset, &
                             &                    core_charge_found )
                        if ( core_charge_found == 1 ) rhcorpw(:,it) = rhcr(:)
                     endif
                  endif
               endif
            endif
         endif
!!$       end if
      end if

! ==== KT_add === 2014/08/11 & 13.2S
      if(.not.paramset) then
         if ( noncol .and. SpinOrbit_Mode == ByPawPot ) then
            if ( sw_use_rphi_Hsoc_rphi == ON .and. &
                 &     sw_use_ival_for_paw_ps_soc == OFF ) then
               vlocr_pw(1:nmesh(it),it) = vlocr2( 1:nmesh(it) )
            endif
         else if ( sw_spinorbit_second_variation==ON .and. SpinOrbit_Mode==ByPawPot ) then
            if ( sw_use_rphi_Hsoc_rphi == ON .and. &
                 &     sw_use_ival_for_paw_ps_soc == OFF ) then
               vlocr_pw(1:nmesh(it),it) = vlocr2( 1:nmesh(it) )
            endif
         else if ( sw_corelevel_spectrum == ON .or. sw_calc_core_energy == ON ) then
            if ( sw_use_rphi_Hsoc_rphi == ON .and. &
                 &     sw_use_ival_for_paw_ps_soc == OFF ) then
               vlocr_pw(1:nmesh(it),it) = vlocr2( 1:nmesh(it) )
            endif
         endif
      endif
! =============== 2014/08/11 & 13.2S

    end subroutine case_when_pot_has_nosoc

    subroutine case_when_pot_has_soc
      call rd_Vloc_psWF_Q_chir_qitg_soc      !-(contained here) (potrvd)
                                             !        -> ps_xctype, etc.
      call mkindx_arrays_nlmt_2_j_l_m_t       ! (tlmst1)
      if (sw_orb_popu == ON) then
         call mkindx_arrays_nlmtphi_2_j_l_m_t   ! (tlmst1_phi)
      end if

      if (is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) then
         if (sw_use_add_proj == ON) then
            call check_l_t_projectors
            call mkindx_arrays_nlmtadd_2_j_l_m_t     ! (tlmst1_add)
         end if
         if(intzaj == by_pseudo_atomic_orbitals) then
            call mkindx_arrays_nlmtpao_2_j_l_m_t     ! (tlmst1_pao)
         end if
      end if

! =========================== added by K. Tagami ======================= 5.0
!!$    if(is_gncpp == pp_GNCPP2.and. sw_mix_charge_hardpart)then
!!$    if(is_gncpp == pp_GNCPP2.and. sw_mix_charge_hardpart == ON)then
      if((is_gncpp == pp_GNCPP2 .or. is_gncpp==pp_PAW .and. PAW_switch == OFF)&
           & .and. sw_mix_charge_hardpart == ON)then
         if(.not.associated(ia2ia_symmtry_op)) then
            allocate(ia2ia_symmtry_op(natm,nopr+af))
            ia2ia_symmtry_op=0
            call set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia_symmtry_op)
         end if
         if(.not.associated(ia2ia_symmtry_op_inv)) then
            allocate(ia2ia_symmtry_op_inv(natm,nopr+af))
            ia2ia_symmtry_op_inv=0
            call set_ia2ia_symmtry_op_inv &
                 (natm,nopr,ia2ia_symmtry_op,ia2ia_symmtry_op_inv)
         end if
!!      call make_ia2ia_symmtry_op_etc
      endif
! ===================================================================== 5.0
      if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
         call phase_error_with_msg(nfout,"kt: Not supported, paw ",__LINE__,__FILE__)

         call make_index_arrays_nlt2l_t
         call make_ia2ia_symmtry_op_etc
         call cnstrct_of_PiPjPlPm_k_Qijk_etc
         if(.not.flg_symmtry) then
            call cnstrct_of_CijkClmkVVVVijlm_k
         else
            call cnstrct_of_CijkClmnVVVVijlm_kn
         end if
      end if

      if(.not.paramset) then
         call cnstrct_LocalWF_Dion_and_q_soc           ! betavd
                                                      ! -> betar, dion, q
         call where_does_WFrhvr_damp(nfout,it,mesh_t) ! -> mesh_t
         call coulomb_potential_in_Rspace(5)    ! -> vvv()   (ptclmb)
         call atomic_xc_potential_in_Rspace     ! -> vxc()   (ptxc)
         ! jfp = 0 (non-spin system for atomic charge) y
         call vlocr_plus_hartree_xc2           ! vlocr = vlocr2, vlocr, vvv, vxc
#ifdef _POT_SMOOTHING_
         call smoothing_vlocr()                 ! vlocr smoothing
#endif
         if(is_gncpp == pp_PAW .and. ipaw(it)==ON) then
            call cnstrct_of_VHij_VHpsij_Kinij
            call cnstrct_of_dion_kin_ion
            call init_of_dion_paw
         end if

         if (initial_chg == from_PSEUDOPOTENTIAL_FILE.or.sw_charge_predictor==ON) &
              & call ft_valence_charge(it,nfout,gr_l,ngshell,ngshell_range)
                                                     ! rhvr -> rhvg_l
         if(sw_positron /= OFF) &
              &   call rd_core_charge_then_ft(nfp,it,nfout,gr_l,paramset,is_gncpp,&
              &                               core_charge_found  )
                                                     ! rhcr -> rhcg_l
      end if

    end subroutine case_when_pot_has_soc
! =========================================================================  11.0

    subroutine skip_commentlines_pp(nf)
      integer, intent(in) :: nf
      integer, parameter  ::    len_str = 80
      character(len=len_str) :: str
      logical :: comment_statement
      comment_statement = .true.
      do while(comment_statement)
         read(nfp,'(a80)',end=1000) str
         if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
              & .or. str(1:1) == '%' .or. str(1:1) == '*') then
            if(iprippex >= 1) write(nfout,'(a80)') str
         else if(len(trim(str)) == 0) then
            if(iprippex >= 1) write(nfout,'(a80)') str
         else
            comment_statement = .false.
         endif
      enddo
      backspace(nfp)
      return
1000  write(nfout,'(" eof is reached")')
      call phase_error_with_msg(nfout,' eof is reached in PP reading',__LINE__,__FILE__)
    end subroutine skip_commentlines_pp

   !================================================
    subroutine rd_Vloc_psWF_Q_then_chir_qitg
   !================================================
      integer       :: i, irc, il, t1
      real(kind=DP) :: dummy

      call read_natomn_ival_iloc_itpcc()

      if ( intzaj == by_pseudo_atomic_orbitals ) then
         if ( is_gncpp == pp_PAW ) call set_irank_minimal_pawpot( nfp, it, 0 )
      endif

      call read_ps_xctype()
      call read_alp_cc()
      call read_nmesh_xh_rmax()

      if(mype == 0) then
         vall_is_detected = check_vall()
         if(vall_is_detected) then
            if(paramset) then
               read(nfp,*) (dummy,i=1,nmesh(it))
            else
               read(nfp,*) (vvv(i),i=1,nmesh(it)) ! ### VAE[scr](r) ### Now vvv is VAE[scr](r)
               if(is_gncpp == pp_PAW .and. ipaw(it) == ON) vloc_scr_ae(1:nmesh(it))=vvv(1:nmesh(it))
            end if
         endif
         if(paramset) then
            read(nfp,*) (dummy,i=1,nmesh(it))
            read(nfp,*) (dummy,i=1,nmesh(it))
            if(is_gncpp==pp_PAW .or. is_gncpp==pp_GNCPP2) read(nfp,*) (dummy,i=1,nmesh(it))
         else
            read(nfp,*) (vlocr(i) ,i=1,nmesh(it)) ! ### Vloc[scr](r) ###
            if(is_gncpp==pp_PAW .and. ipaw(it) == ON) vloc_scr_ps(1:nmesh(it))=vlocr(1:nmesh(it))
            if(is_gncpp==pp_PAW .or. is_gncpp==pp_GNCPP2) &
                 & read(nfp,*) (vlocr2(i),i=1,nmesh(it)) ! ### Vloc[ion](r) ###
            read(nfp,*) (rhvr(i)  ,i=1,nmesh(it)) ! ### 4*pi*r*r*n(r) ###
            call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,h(it)) ! -(b_P.P.) -> radr,h(it)
            if(is_gncpp==pp_PAW ) then
               if ( ipaw(it) == ON .or. sw_add_corecharge_rspace == ON ) then
                  radr_paw(:,it)=radr(:)
               endif
            else
               if ( sw_add_corecharge_rspace == ON ) radr_paw(:,it)=radr(:)
            endif

            if ( sw_write_rotated_orbitals == ON ) radr_paw(:,it)=radr(:)
            if ( sw_hubbard == ON .and. initial_occmat == CRYSTAL_FIELD_APPROX ) &
                 &               radr_paw(:,it)=radr(:)
            ! radr(i) = exp(zmax+(i-nmesh(it))*h(it)) = rmax*exp(h(it)*(i-nmesh(it)))
!!! Epsilon
            !!if(is_gncpp == pp_GNCPP2 .and. sw_use_add_proj == ON) then
            if((is_gncpp==pp_PAW .or. is_gncpp == pp_GNCPP2) .and. sw_use_add_proj == ON) then
            if(rcproj(it) <= 0.d0 ) then
               do i=nmesh(it),1,-1
                  if(dabs(vvv(i)-vlocr(i)) .gt. 1.d-13) then
                     irc = i+1
                     exit
                  end if
               end do
               write(nfout,*) 'Rc(local potential)=',radr(irc)
               write(nfout,*) 'irc=',irc
               rcproj(it) = radr(irc)
               kappa(it)  = 10
            end if
            end if
!!! Epsilon
         end if
      end if
      if(npes > 1 .and. .not.paramset) then
         call mpi_bcast(vvv,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vlocr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         if(is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) &
              & call mpi_bcast(vlocr2,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(rhvr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(h(it),1,        mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(radr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vall_is_detected,1,mpi_logical,0,MPI_CommGroup,ierr)

! ============================ ASMS_DEBUG ========================= 2013/02/07
!         if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
!            call mpi_bcast(vloc_scr_ae,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!            call mpi_bcast(vloc_scr_ps,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!!$      J. Koga 09.06.02
!            call mpi_bcast(radr_paw(:,it),nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!!! Epsilon addproj
!         !!else if(is_gncpp == pp_GNCPP2 .and. sw_use_add_proj == ON) then
!         else if((is_gncpp==pp_PAW .or. is_gncpp == pp_GNCPP2) .and. sw_use_add_proj == ON) then
!            call mpi_bcast(rcproj(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!            call mpi_bcast(kappa(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!         end if
!!!! Epsilon addproj

! -----------
         if( is_gncpp == pp_PAW .and. ipaw(it) == ON ) then
            call mpi_bcast( vloc_scr_ae, nmesh(it), mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
            call mpi_bcast( vloc_scr_ps, nmesh(it), mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
            call mpi_bcast( radr_paw(:,it), nmesh(it),mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
         else if ( sw_write_rotated_orbitals == ON ) then
            call mpi_bcast( radr_paw(:,it), nmesh(it),mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
         else if ( sw_hubbard == ON .and. initial_occmat == CRYSTAL_FIELD_APPROX ) then
            call mpi_bcast( radr_paw(:,it), nmesh(it),mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
         endif
         if( is_gncpp == pp_PAW .and. sw_add_corecharge_rspace == ON ) then
            call mpi_bcast( radr_paw(:,it), nmesh(it),mpi_double_precision, 0, &
                 &          MPI_CommGroup, ierr )
         endif

         if( (is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) &
              &                  .and. sw_use_add_proj == ON ) then

            call mpi_bcast(rcproj(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(kappa(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
         end if
! ============================ ASMS_DEBUG ========================= 2013/02/07

      end if

      call determine_lpsmax ! -> (lpsmax) 3=d or 4=f
      if(iprippex >= 2) write(nfout,'(" !PP lpsmax = ",i8)') lpsmax(it)

! ================================== added by K. Tagami ================ 11.0
      call calc_nums_angmom_on_atomtype
      if(ipri >= 2) then
        write(nfout,*) ' ! nums_of_angmom_on_atomtype = ', &
             &            nums_of_angmom_on_atomtype(it)
      endif
! ====================================================================== 11.0

      call rd_itau_etc_then_phir_chir !->(phir,chir):pWF, chir:(e-T-Vloc)*phir
      if (ret) return

      if(.not.paramset) then
         call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr&
              & ,wos)
         !if(sw_hybrid_functional==ON.and.sw_rspace_hyb==OFF) call alloc_qitg_exx()
      end if
      call rd_qrsps_then_iqitg_and_qitgft(is_gncpp,nfp,it,nfout,gr_l,ngshell,ngshell_range &
           & ,paramset,mmt) ! m_PP

! ===
      if ( is_gncpp == pp_GNCPP2 .or. is_gncpp == pp_PAW ) then
         if (.not.paramset ) then
            if ( intzaj == by_pseudo_atomic_orbitals ) then
               do il=1,lpsmax(it)
                  do t1=1,itau(il,it)
                     paor(1:nmesh(it),il,t1,it) = phir(1:nmesh(it),il,t1)
                  end do
               end do
            end if
         endif

         if ( intzaj == by_pseudo_atomic_orbitals ) then
            if ( is_gncpp == pp_PAW ) then
               call set_irank_minimal_pawpot( nfp, it, 1 )
            else
               call set_irank_default( it )
            endif
         end if
      endif

    end subroutine rd_Vloc_psWF_Q_then_chir_qitg

! ===================================== added by K. Tagami ================ 11.0
    subroutine rd_Vloc_psWF_Q_chir_qitg_soc
      integer       :: i, irc
      real(kind=DP) :: dummy

      call read_natomn_ival_iloc_itpcc()
      call read_ps_xctype()
      call read_alp_cc()
      call read_nmesh_xh_rmax()

      if(mype == 0) then
         vall_is_detected = check_vall()
         if(vall_is_detected) then
            if(paramset) then
               read(nfp,*) (dummy,i=1,nmesh(it))
            else
               read(nfp,*) (vvv(i),i=1,nmesh(it))
                                       ! ### VAE[scr](r) ### Now vvv is VAE[scr](r)
               if(is_gncpp == pp_PAW .and. ipaw(it) == ON) &
                    &            vloc_scr_ae(1:nmesh(it))=vvv(1:nmesh(it))
            end if
         endif
         if(paramset) then
            read(nfp,*) (dummy,i=1,nmesh(it))
            read(nfp,*) (dummy,i=1,nmesh(it))
            if(is_gncpp==pp_PAW .or. is_gncpp==pp_GNCPP2) &
                 &             read(nfp,*) (dummy,i=1,nmesh(it))
         else
            read(nfp,*) (vlocr(i) ,i=1,nmesh(it)) ! ### Vloc[scr](r) ###
            if(is_gncpp==pp_PAW .and. ipaw(it) == ON) &
                 &                  vloc_scr_ps(1:nmesh(it))=vlocr(1:nmesh(it))

            if(is_gncpp==pp_PAW .or. is_gncpp==pp_GNCPP2) &
                 & read(nfp,*) (vlocr2(i),i=1,nmesh(it)) ! ### Vloc[ion](r) ###

            read(nfp,*) (rhvr(i)  ,i=1,nmesh(it)) ! ### 4*pi*r*r*n(r) ###
            call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,h(it))
                                              ! -(b_P.P.) -> radr,h(it)

            if (is_gncpp==pp_PAW .and. ipaw(it) == ON ) radr_paw(:,it)=radr(:)

            ! radr(i) = exp(zmax+(i-nmesh(it))*h(it)) = rmax*exp(h(it)*(i-nmesh(it)))
!!! Epsilon
            if( (is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2 ) &
                 &  .and. sw_use_add_proj == ON ) then

               if(rcproj(it) <= 0.d0 ) then
                  do i=nmesh(it),1,-1
                     if(dabs(vvv(i)-vlocr(i)) .gt. 1.d-13) then
                        irc = i+1
                        exit
                     end if
                  end do
                  write(nfout,*) 'Rc(local potential)=',radr(irc)
                  write(nfout,*) 'irc=',irc
                  rcproj(it) = radr(irc)
                  kappa(it)  = 10
               end if
            end if
!!! Epsilon
         end if
      end if

      if(npes > 1 .and. .not.paramset) then
         call mpi_bcast(vvv,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vlocr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         if(is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) &
              & call mpi_bcast( vlocr2,nmesh(it),mpi_double_precision,0, &
              &                 MPI_CommGroup,ierr )
         call mpi_bcast(rhvr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(h(it),1,        mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(radr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vall_is_detected,1,mpi_logical,0,MPI_CommGroup,ierr)

         if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
            call mpi_bcast(vloc_scr_ae,  nmesh(it),mpi_double_precision,0, &
                 &         MPI_CommGroup,ierr)
            call mpi_bcast(vloc_scr_ps,  nmesh(it),mpi_double_precision,0, &
                 &         MPI_CommGroup,ierr)
            call mpi_bcast(radr_paw(:,it),nmesh(it),mpi_double_precision,0, &
                 &         MPI_CommGroup,ierr)
!!! Epsilon addproj

         else if( (is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2 ) &
              &       .and. sw_use_add_proj == ON ) then
            call mpi_bcast(rcproj(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(kappa(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
         end if
!!! Epsilon addproj

      end if

      call determine_lpsmax ! -> (lpsmax) 3=d or 4=f
      if(iprippex >= 2) write(nfout,'(" !PP lpsmax = ",i8)') lpsmax(it)

      call calc_nums_angmom_on_atomtype
      if(ipri >= 2) then
        write(nfout,*) ' ! nums_of_angmom_on_atomtype = ', nums_of_angmom_on_atomtype(it)
      endif

      call rd_itau_etc_then_phir_chir_soc

      if(.not.paramset) then
         call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr&
              & ,wos)
         !if(sw_hybrid_functional==ON.and.sw_rspace_hyb==OFF) call alloc_qitg_exx()   ! ????????????
      end if

      call rd_qrsps_iqitg_and_qitgft_soc( is_gncpp,nfp,it,nfout,gr_l, &
              &                              ngshell, ngshell_range, &
              &                              paramset, mmt ) ! m_PP

    end subroutine rd_Vloc_psWF_Q_chir_qitg_soc
! ====================================================================== 11.0

   !=============================================
    subroutine rd_itau_etc_then_phir_chir
   !=============================================
      integer       :: il, t1, nrc, i, nm, mord, imax, iatomn_t, nrc0
      real(kind=DP) :: dummy,eps_phi0,eps_phi, drc, dk, fac
      integer, parameter :: PRINTLEVEL = 2

      nm = nmesh(it)
      if(iprippex >= PRINTLEVEL) write(nfout,'(" !PP lpsmax = ",i8)') lpsmax(it)
      loop_L : do il = 1, lpsmax(it)
         if(mype == 0) &
              & call read_itau_ivanl(nfp,nfout,iloc(it),il,itau(il,it)&
              &   ,ivanl(il,it),iprippex) ! -(b_P.P.) ->(itau,ivanl)
         if(npes > 1) then
            call mpi_bcast(itau(il,it), 1,mpi_integer,0,MPI_CommGroup,ierr)
            call mpi_bcast(ivanl(il,it),1,mpi_integer,0,MPI_CommGroup,ierr)
         end if

         loop_tau : do t1 = 1, itau(il,it)
            if(iprippex >= PRINTLEVEL) write(nfout,*) ' t1 = ', t1
            if(mype == 0) then
!               call read_tau_eps_nrc_mord(nfp,nfout,kord &
!                    &                           ,il,t1,eps(il,t1,it),nrc,mord) ! ->(eps,t1,it)
               call read_tau_eps_nrc_mord_nrc0(nfp,nfout,kord &
                    &                           ,il,t1,eps(il,t1,it),nrc,mord,nrc0) ! ->(eps,t1,it)

               if ( is_gncpp==pp_PAW .and. .not. paramset ) then
!!!                  if ( ipaw(it) == ON .or. sw_rspace_ekin_density == ON ) then
                  if ( ipaw(it) == ON .or. sw_rspace_ekin_density == ON &
                       &              .or. sw_orb_popu == ON ) then
                     if(nrc>0) then
                        wf_nrc(il,t1,it)=nrc
                     else if(nrc0>0) then
                        wf_nrc(il,t1,it)=nrc0
                     end if
                  end if
               end if
            end if
            if(is_gncpp == pp_GNCPP2) then
               if(npes > 1) call mpi_bcast(eps(il,t1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
               if(mype == 0 .and. iprippex >= PRINTLEVEL) write(nfout,'(" !PP ",2i4,f12.8,2i6 &
                    & ," : il, ntau, eps, nrc, mord")') il,t1, eps(il,t1,it), nrc,mord
            end if

            if(paramset .and. mype == 0) then
               if(is_gncpp == pp_GNCPP2) then
                  if(ae_wavefunctions_are_detected(it)) read(nfp,*) (dummy,i=1,nm)
               else if(is_gncpp==pp_PAW) then
                  read(nfp,*) (dummy,i=1,nm);
               end if
               if(nrc == 0) then
                  read(nfp,*) (dummy,i=1,nm);      read(nfp,*) (dummy,i=1,nm)
               else if(nrc > 0) then
                  read(nfp,*) (dummy,i=0,mord);    read(nfp,*) (dummy,i=nrc+1,nm)
               end if
            else if(.not.paramset) then
               if(mype == 0) then
                  if(nrc == 0) then
                     if(is_gncpp == pp_GNCPP2) then
                        if(ae_wavefunctions_are_detected(it)) &
                             & read(nfp,*) (phir_ae(i,il,t1),i=1,nm) ! phir_ae = r*phi_ae[n](r)
                     else if(is_gncpp==pp_PAW) then
                        if(ipaw(it) == ON) then
                           read(nfp,*) (psir(i,il,t1),i=1,nm) ! psirpw = r*psi[n](r)
                        else if(ae_wavefunctions_are_detected(it)) then
                           read(nfp,*) (phir_ae(i,il,t1),i=1,nm) ! phir_ae = r*phi_ae[n](r)
                        else
                           read(nfp,*) (dummy,i=1,nm)
                        end if
                     end if
                     read(nfp,*) (phir(i,il,t1),i=1,nm) ! phir = r*phi[n](r)
                     read(nfp,*) (chir(i,il,t1),i=1,nm) ! chir = Vion[n](r)
                     if(is_gncpp == pp_GNCPP1) then
                        chir(1:nm,il,t1) = (chir(1:nm,il,t1)-vlocr(1:nm)) * phir(1:nm,il,t1)
                     else
                        chir(1:nm,il,t1) = (chir(1:nm,il,t1)-vlocr2(1:nm)) * phir(1:nm,il,t1)
                        ! chir = Vion[n](r), vloc2 = Vloc[ion](r) --> chir = r*chi[n](r)
                     end if
                     if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
                        if(nrc0 == 0) then
                           eps_phi0=1.d0
                           !do i=nm,1,-1
                           !   write(0,'(i5,3i3,3f20.15)') i,it,il,t1,phir(i,il,t1),psir(i,il,t1),phir(i,il,t1)-psir(i,il,t1)
                           !enddo
                           !write(0,*)
                           do i=nm,1,-1
                              if(psir(i,il,t1)-phir(i,il,t1).eq.0.d0) then
                                 eps_phi=1.d0
                              else
                                 eps_phi=abs((psir(i,il,t1)-phir(i,il,t1))/psir(i,il,t1))
                              end if
                              ! if(abs(psir(i,il,t1)-phir(i,il,t1)) .gt. eps_phi) then
                              if(eps_phi.gt.5.d0*eps_phi0 ) then
                                 !if(modified_nrc_resolution==OFF .or. abs(psir(i,il,t1)-phir(i,il,t1)).gt.1.d-10) then
                                 if (abs(psir(i,il,t1)-phir(i,il,t1)).gt.1.d-10) then
!!$wf_nrc(il,t1,it)=i+1
                                    if(i>=nm-1) then
                                       wf_nrc(il,t1,it) = i+1
                                    else
                                       wf_nrc(il,t1,it)=i+2
                                    endif
                                    if(iprippex>=2) then
                                       write(nfout,'(a,3i3,a,i5)') 'wf_nrc for ',it,il,t1," : ",i+2
                                    endif
                                    exit
                                 end if
                              end if
                              eps_phi0=eps_phi
                           end do
                        end if
! === KT_add == 2014/09/19
                     else
!                        if (is_gncpp == pp_PAW .and. sw_rspace_ekin_density == ON) then
                        if (is_gncpp == pp_PAW .and.  &
                             &  ( sw_rspace_ekin_density == ON .or. &
                             &    sw_orb_popu == ON ) ) then

                           if(nrc0 == 0) then
                              eps_phi0=1.d0
                              do i=nm,1,-1
                                 if(phir_ae(i,il,t1)-phir(i,il,t1).eq.0.d0) then
                                    eps_phi=1.d0
                                 else
                                    eps_phi = abs((phir_ae(i,il,t1)-phir(i,il,t1)) &
                                         &    /phir_ae(i,il,t1))
                                 end if
                                 if(eps_phi.gt.5.d0*eps_phi0 ) then
                                    if (abs(phir_ae(i,il,t1)-phir(i,il,t1)).gt.1.d-10) then
                                       if(i>=nm-1) then
                                          wf_nrc(il,t1,it) = i+1
                                       else
                                          wf_nrc(il,t1,it)=i+2
                                       endif
                                       if(iprippex>=2) then
                                          write(nfout,'(a,3i3,a,i5)') 'wf_nrc for ', &
                                               &                      it,il,t1," : ",i+2
                                       endif
                                       exit
                                    end if
                                 end if
                                 eps_phi0=eps_phi
                              end do
                           end if
                        endif
! =========== 2014/09/19
                     end if

                  else if(nrc > 0) then
                     if(is_gncpp == pp_GNCPP2) then
                        if(ae_wavefunctions_are_detected(it)) &
                             & read(nfp,*) (phir_ae(i,il,t1),i=1,nm) ! phir_ae = r*phi_ae[n](r)
                     else if(is_gncpp==pp_PAW) then
                        if(ipaw(it) == ON) then
                           read(nfp,*) (psir(i,il,t1),i=1,nm) ! psirpw = r*psi[n](r)
                        else if(ae_wavefunctions_are_detected(it))then
                           read(nfp,*) (phir_ae(i,il,t1),i=1,nm) ! phir_ae = r*phi_ae[n](r)
                        else
                           read(nfp,*) (dummy,i=1,nm)
                        end if
                     end if
                     read(nfp,*) (copsc(i),i=0,mord)
                     read(nfp,*) (phir(i,il,t1),i=nrc+1,nm) ! phir = r*phi[n](r)
                     call cnvrtp(nm,1,nrc,il,mord,copsc,radr,phir(1,il,t1))   !-(b_P.P._f77.f)
                     call cnvrtc(nm,1,nrc,il-1,mord,copsc,radr,eps(il,t1,it)& !-(b_P.P._f77.f)
                          & ,phir(1,il,t1),vlocr,chir(1,il,t1))               ! vlocr = Vloc[scr](r)
                     if(vall_is_detected) then
                        chir(nrc+1:nm,il,t1) = (vvv(nrc+1:nm) &        ! vvv = VAE[scr](r)
                             & - vlocr(nrc+1:nm))*phir(nrc+1:nm,il,t1) ! vlocr = Vloc[scr](r), phir = r*phi[n](r)
                                                                       ! --> chir = r*chi[n](r)
                     else
                        chir(nrc+1:nm,il,t1) = 0.d0
                     end if
                  end if
!                  phirpw(1:nm,il,t1,it)=phir(1:nm,il,t1)
               end if
            end if
!!$            if(is_gncpp == pp_PAW) then
               if(npes > 1) then
                  call mpi_bcast(eps(il,t1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
                  if(.not.paramset ) then
!                     if ( ipaw(it) == ON .or. sw_rspace_ekin_density == ON ) then
                     if ( ipaw(it) == ON .or. sw_rspace_ekin_density == ON &
                          &              .or. sw_orb_popu == ON ) then

                        call mpi_bcast(wf_nrc(il,t1,it),1,mpi_integer,0,MPI_CommGroup,ierr)
                     endif
                  endif
               end if
!!$            end if

         enddo loop_tau
      enddo loop_L
      if(.not.paramset .and. npes > 1) then
         if(is_gncpp == pp_GNCPP2) then
            if(ae_wavefunctions_are_detected(it)) &
                 &  call mpi_bcast(phir_ae,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         else if(is_gncpp == pp_PAW) then
            if(ipaw(it) /= ON .and. ae_wavefunctions_are_detected(it)) &
                 &  call mpi_bcast(phir_ae,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         end if
         call mpi_bcast(phir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(chir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)

         if(is_gncpp == pp_PAW .and. ipaw(it) == ON) &
              & call mpi_bcast(psir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
      end if

!!$      if(is_gncpp /= pp_PAW .and. .not.paramset .and. sw_orb_popu == 1) then
      if(.not.paramset .and. sw_orb_popu == 1) then
         if ( is_gncpp == pp_PAW .and. howto_set_proj_radius == 1 ) then
            call set_projector_radius_if_undef( nfout, it )
         endif

         call make_phi_window_parameter(nfout,it)
         do il=1,lpsmax(it)
            do t1=1,itau(il,it)

! ========================ASMS DEBUG ========= 2013/12/06
!               dk = k_phi(il,t1,it)
               dk = k_phi(il,t1,it) **2
! ============================================ 2013/12/06

               if(iprippex >= PRINTLEVEL) then
                  write(nfout,*) ' !PP il=',il,' t1=',t1,' dk=',dk
               end if
               if(dk > 1.d-5) then
                  do i=2,nm
                     if(radr(i) > rc_phi(il,t1,it)) then
                        imax = i-1
                        drc = radr(imax)
                        exit
                     end if
                  end do
                 if(iprippex >= PRINTLEVEL) write(nfout,*) ' !PP imax=',imax,' drc=',drc
                  fac = 1.d0/(1.d0-exp(-dk*drc**2))
                  phirt(1:imax,il,t1,it) = &
                  & fac*(1.d0-exp(-dk*(radr(1:imax)-drc)**2))*phir(1:imax,il,t1)
                  phirt(imax+1:nm,il,t1,it) = 0.d0
               else
                  if(iprippex >= PRINTLEVEL) write(nfout,*) '!PP without a window function'
                  phirt(1:nm,il,t1,it) = phir(1:nm,il,t1)
               end if
            end do
         end do
      end if

      if(is_gncpp == pp_PAW) then
         if(.not.paramset) then
            if ( ipaw(it) == ON ) then
               psirpw(1:mmesh,1:nloc,1:ntau,it)=psir(1:mmesh,1:nloc,1:ntau)
               phirpw(1:mmesh,1:nloc,1:ntau,it)=phir(1:mmesh,1:nloc,1:ntau)
               wf_mnrc(it)=maxval(wf_nrc(:,:,it))
               if(iprippex >= 2) then
                  write(nfout,'(" wf_mnrc(",i3," ) = ",i8,"<<rd_itau_etc_then_phir_chir>>")') it, wf_mnrc(it)
               end if

            else if ( sw_rspace_ekin_density == ON ) then
               psirpw(1:mmesh,1:nloc,1:ntau,it)=phir_ae(1:mmesh,1:nloc,1:ntau)
               phirpw(1:mmesh,1:nloc,1:ntau,it)=phir(1:mmesh,1:nloc,1:ntau)
               wf_mnrc(it)=maxval(wf_nrc(:,:,it))
               if(iprippex >= 2) then
                  write(nfout,*) "wf_mnrc(", it, ")=", wf_mnrc(it)
               endif
            else if ( sw_excitation == ON .or. sw_wannier90 == ON ) then
               psirpw(1:mmesh,1:nloc,1:ntau,it)=phir_ae(1:mmesh,1:nloc,1:ntau)
               phirpw(1:mmesh,1:nloc,1:ntau,it)=phir(1:mmesh,1:nloc,1:ntau)
            else if ( sw_hubbard == ON .and. initial_occmat==CRYSTAL_FIELD_APPROX ) then
               psirpw(1:mmesh,1:nloc,1:ntau,it)=phir_ae(1:mmesh,1:nloc,1:ntau)
            endif
         end if
      end if

! ======== KT_add ========== 13.0U2
      if ( is_gncpp == pp_PAW .and. (.not.paramset) ) then
         if ( sw_modified_TFW_functional == ON .or. sw_use_contracted_psir == ON ) then
            if ( ipaw(it) == ON ) then
               psir_val(1:mmesh,1:nloc,1:ntau,it) = psir(1:mmesh,1:nloc,1:ntau)
            else
               psir_val(1:mmesh,1:nloc,1:ntau,it) = phir_ae(1:mmesh,1:nloc,1:ntau)
            endif
         endif
      endif
! ========================== 13.0U2

      if(paramset) then
!         if(it == 1) ntau = itau(1,it)
!         do il = 2, lpsmax(it)
!            if(ntau < itau(il,it)) ntau = itau(il,it)
!         end do
         ntau = maxval(itau)
      end if

      if(ipri == PRINTLEVEL) write(nfout,'(" end of <<rd_itau_etc_then_phir_chir>>")')
    end subroutine rd_itau_etc_then_phir_chir

! ===================================== added by K. Tagami ============= 11.0
    subroutine rd_itau_etc_then_phir_chir_soc
      integer       :: il, t1, nrc, i, nm, mord, imax, iatomn_t
      integer :: nrc0

      real(kind=DP) :: dummy,eps_phi0,eps_phi, drc, dk, fac
      integer, parameter :: PRINTLEVEL = 2

      integer :: kj, il1
      integer :: jval

      nm = nmesh(it)
      if (iprippex >= PRINTLEVEL) then
         write(nfout,'(" !PP kjmax = ",i8)') nums_of_angmom_on_atomtype( it )
      endif

!--
      loop_kj : Do kj =1, nums_of_angmom_on_atomtype( it )
         il = get_lp1_in_AngMomList( pot_has_soc(it), kj )
         jval = get_jph_in_AngMomList( pot_has_soc(it), kj )

         if (mype == 0) &
              & call read_itau_ivanl( nfp, nfout, iloc(it), kj, itau(kj,it), &
              &                       ivanl(kj,it), iprippex )
                                              ! -(b_P.P.) ->(itau,ivanl)
         if (npes > 1) then
            call mpi_bcast( itau(kj,it), 1, mpi_integer, 0, MPI_CommGroup, ierr )
            call mpi_bcast( ivanl(kj,it),1, mpi_integer, 0, MPI_CommGroup, ierr )
         end if

         loop_tau : do t1 = 1, itau(kj,it)
            if (iprippex >= PRINTLEVEL) write(nfout,*) ' t1 = ', t1
            if (mype == 0) then
!!!!               call read_tau_eps_nrc_mord( nfp, nfout, kord, &
!!!                    &                      kj, t1, eps(kj,t1,it), nrc, mord )
               call read_tau_eps_nrc_mord_nrc0( nfp, nfout, kord, &
                    &                      kj, t1, eps(kj,t1,it), nrc, mord, nrc0 )
                                             ! ->(eps,t1,it)
               if (is_gncpp==pp_PAW .and. ipaw(it) == ON) then
                  if(.not.paramset) then
                     if ( nrc > 0 ) then
                        wf_nrc(kj,t1,it)=nrc
                     else
                        wf_nrc(kj,t1,it)=nrc0
                     endif
                  end if
               endif
            end if

            if (is_gncpp == pp_GNCPP2) then
               if (npes > 1) then
                  call mpi_bcast( eps(kj,t1,it),1,mpi_double_precision, &
                       &                        0, MPI_CommGroup,ierr )
               endif
               if (mype == 0 .and. iprippex >= PRINTLEVEL) then
                  write(nfout,'(" !PP ",2i4,f12.8,2i6 &
                       & ," : kj, ntau, eps, nrc, mord")') kj,t1, &
                       &      eps(kj,t1,it), nrc,mord
               endif
            end if

            if(paramset .and. mype == 0) then
               if(is_gncpp == pp_GNCPP2) then
                  if(ae_wavefunctions_are_detected(it)) read(nfp,*) (dummy,i=1,nm)
               else if(is_gncpp==pp_PAW) then
                  read(nfp,*) (dummy,i=1,nm);
               end if
               if(nrc == 0) then
                  read(nfp,*) (dummy,i=1,nm);      read(nfp,*) (dummy,i=1,nm)
               else if(nrc > 0) then
                  read(nfp,*) (dummy,i=0,mord);    read(nfp,*) (dummy,i=nrc+1,nm)
               end if
            else if(.not.paramset) then
               if(mype == 0) then
                  if (nrc == 0) then
                     if(is_gncpp == pp_GNCPP2) then
                        if(ae_wavefunctions_are_detected(it)) &
                             & read(nfp,*) (phir_ae(i,kj,t1),i=1,nm)
                                                       ! phir_ae = r*phi_ae[n](r)
                     else if(is_gncpp==pp_PAW) then
                        if(ipaw(it) == ON) then
                           read(nfp,*) (psir(i,kj,t1),i=1,nm)
                                                ! psirpw = r*psi[n](r)
                        else if(ae_wavefunctions_are_detected(it)) then
                           read(nfp,*) (phir_ae(i,kj,t1),i=1,nm)
                                                ! phir_ae = r*phi_ae[n](r)
                        else
                           read(nfp,*) (dummy,i=1,nm)
                        end if
                     end if
                     read(nfp,*) (phir(i,kj,t1),i=1,nm) ! phir = r*phi[n](r)
                     read(nfp,*) (chir(i,kj,t1),i=1,nm) ! chir = Vion[n](r)

                     if(is_gncpp == pp_GNCPP1) then
                        chir(1:nm,kj,t1) = (chir(1:nm,kj,t1)-vlocr(1:nm)) &
                             &              * phir(1:nm,kj,t1)
                     else
                        chir(1:nm,kj,t1) = (chir(1:nm,kj,t1)-vlocr2(1:nm)) &
                             &              * phir(1:nm,kj,t1)
                                       ! chir = Vion[n](r), vloc2 = Vloc[ion](r)
                                       !       --> chir = r*chi[n](r)
                     end if

                     if(is_gncpp == pp_PAW .and. ipaw(it) == ON) then
                        if ( nrc0 == 0 ) then
                           eps_phi0=1.d0
                           do i=nm,1,-1
                              if(psir(i,kj,t1)-phir(i,kj,t1).eq.0.d0) then
                                 eps_phi=1.d0
                              else
                                 eps_phi=abs((psir(i,kj,t1)-phir(i,kj,t1)) &
                                      &  /psir(i,kj,t1))
                              end if

                              if(eps_phi.gt.5.d0*eps_phi0) then
                                 if ( abs(psir(i,kj,t1)-phir(i,kj,t1)).gt.1.0d-10 ) then
                                    if ( i >= nm-1 ) then
                                       wf_nrc(kj,t1,it)=i+1
                                    else
                                       wf_nrc(kj,t1,it)=i+2
                                    endif
                                    exit
                                 end if
                              endif
                              eps_phi0=eps_phi
                           end do
                        end if
                     endif

                  else if(nrc > 0) then
                     if(is_gncpp == pp_GNCPP2) then
                        if(ae_wavefunctions_are_detected(it)) &
                             & read(nfp,*) (phir_ae(i,kj,t1),i=1,nm)
                                                 ! phir_ae = r*phi_ae[n](r)
                     else if(is_gncpp==pp_PAW) then
                        if(ipaw(it) == ON) then
                           read(nfp,*) (psir(i,kj,t1),i=1,nm)
                                                  ! psirpw = r*psi[n](r)
                        else if(ae_wavefunctions_are_detected(it))then
                           read(nfp,*) (phir_ae(i,kj,t1),i=1,nm)
                                                   ! phir_ae = r*phi_ae[n](r)
                        else
                           read(nfp,*) (dummy,i=1,nm)
                        end if
                     end if

                     read(nfp,*) (copsc(i),i=0,mord)
                     read(nfp,*) (phir(i,kj,t1),i=nrc+1,nm) ! phir = r*phi[n](r)

                     call cnvrtp( nm, 1, nrc, kj, mord, copsc, radr, &
                          &       phir(1,kj,t1) )
                                                         !-(b_P.P._f77.f)
                     call cnvrtc( nm, 1, nrc,kj-1, mord, copsc, radr, &
                          &       eps(kj,t1,it), phir(1,kj,t1), &
                          &       vlocr,chir(1,kj,t1))
                                               ! vlocr = Vloc[scr](r)

                     if(vall_is_detected) then
                        chir(nrc+1:nm,kj,t1) = (vvv(nrc+1:nm) - vlocr(nrc+1:nm)) &
                             &                *phir(nrc+1:nm,kj,t1)
                                           ! vvv = VAE[scr](r)
                                           ! vlocr = Vloc[scr](r), phir = r*phi[n](r)
                                           ! --> chir = r*chi[n](r)
                     else
                        chir(nrc+1:nm,kj,t1) = 0.d0
                     end if
                  end if
!                  phirpw(1:nm,kj,t1,it)=phir(1:nm,kj,t1)
               end if
            end if

            if(npes > 1) then
               call mpi_bcast(eps(kj,t1,it),1,mpi_double_precision,0,&
                    &          MPI_CommGroup,ierr )
               if(.not.paramset .and. ipaw(it) == ON) &
                    & call mpi_bcast(wf_nrc(kj,t1,it),1,mpi_integer,0,&
                    &                MPI_CommGroup,ierr)
            end if

         enddo loop_tau
      enddo loop_kj

      if(.not.paramset .and. npes > 1) then
         if(is_gncpp == pp_GNCPP2) then
            if(ae_wavefunctions_are_detected(it)) &
                 &  call mpi_bcast(phir_ae,mmesh*nloc*ntau,mpi_double_precision,0,&
                 &                 MPI_CommGroup,ierr)
         else if(is_gncpp == pp_PAW) then
            if(ipaw(it) /= ON .and. ae_wavefunctions_are_detected(it)) &
                 &  call mpi_bcast(phir_ae,mmesh*nloc*ntau,&
                 &                 mpi_double_precision,0,MPI_CommGroup,ierr)
         end if
         call mpi_bcast(phir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(chir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)

         if(is_gncpp == pp_PAW .and. ipaw(it) == ON) &
              & call mpi_bcast(psir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
      end if


      if(.not.paramset .and. sw_orb_popu == 1) then
         call make_phi_window_parameter_soc(nfout,it)

         Do kj =1, nums_of_angmom_on_atomtype( it )
            il = get_lp1_in_AngMomList( pot_has_soc(it),kj )

            do t1=1,itau(kj,it)
               dk = k_phi(kj,t1,it)**2
               if(iprippex >= PRINTLEVEL) then
                  write(nfout,*) ' !PP kj=',kj,' t1=',t1,' dk=',dk
               end if
               if(dk > 1.d-5) then
                  do i=2,nm
                     if(radr(i) > rc_phi(kj,t1,it)) then
                        imax = i-1
                        drc = radr(imax)
                        exit
                     end if
                  end do
                  if(iprippex >= PRINTLEVEL) write(nfout,*) ' !PP imax=',imax,' drc=',drc

                  fac = 1.d0/(1.d0-exp(-dk*drc**2))
                  phirt(1:imax,kj,t1,it) = &
                  & fac*(1.d0-exp(-dk*(radr(1:imax)-drc)**2))*phir(1:imax,kj,t1)
                  phirt(imax+1:nm,kj,t1,it) = 0.d0
               else
                  if(iprippex >= PRINTLEVEL) write(nfout,*) '!PP without a window function'
                  phirt(1:nm,kj,t1,it) = phir(1:nm,kj,t1)
               end if
            end do
         end do
      end if

      if(is_gncpp == pp_GNCPP2) then
         if(.not.paramset .and. intzaj == by_pseudo_atomic_orbitals) then

            Do kj =1, nums_of_angmom_on_atomtype( it )
               il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj )
               do t1=1,itau(kj,it)
                  paor(1:nm,kj,t1,it) = phir(1:nm,kj,t1)
               end do
            end do
         end if

         if(intzaj == by_pseudo_atomic_orbitals) then
            Do kj =1, nums_of_angmom_on_atomtype( it )
               il = get_lp1_in_AngMomList( pot_has_soc(it),kj )

               iatomn_t = nint(iatomn(it))
               if(iatomn_t<=4) then !1H~4Be
                  if(il/=1) cycle
               else if(iatomn_t<=10) then !5B~10Ne
                  if(il>2) cycle
               else if(iatomn_t<=12) then !11Na~12Mg
                  if(iatomn_t-nint(ival(it))==10) then
                     if(il/=1) cycle
                  else ! sp semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=18) then !13Al~18Ar
                  if(il>2) cycle
               else if(iatomn_t<=20) then !19K~20Ca
                  if(iatomn_t-nint(ival(it))==18) then
                     if(il/=1) cycle
                  else ! sp semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=30) then !21Sc~30Zn
                  if(iatomn_t-nint(ival(it))==18) then
                     if(il/=1.and.il/=3) cycle
                  else ! p semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=36) then !31Ga~36Kr
                  if(iatomn_t-nint(ival(it))==28) then
                     if(il>2) cycle
                  else ! d semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=38) then !37Rb~38Sr
                  if(iatomn_t-nint(ival(it))==36) then
                     if(il/=1) cycle
                  else ! sp semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=48) then !39Y~48Cd
                  if(iatomn_t-nint(ival(it))==36) then
                     if(il/=1.and.il/=3) cycle
                  else ! p semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=54) then !49In~54Xe
                  if(iatomn_t-nint(ival(it))==46) then
                     if(il>2) cycle
                  else ! d semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=56) then !55Cs~56Ba
                  if(iatomn_t-nint(ival(it))==54) then
                     if(il/=1) cycle
                  else ! sp semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=57) then !57La
                  if(iatomn_t-nint(ival(it))==54) then
                     if(il/=1.and.il/=3) cycle
                  else ! p semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=71) then !58Ce~71Lu
                  if(il/=1.or.il/=3.or.il/=4) cycle
               else if(iatomn_t<=80) then !72Hf~80Hg
                  if(iatomn_t-nint(ival(it))==58) then
                     if(il/=1.and.il/=3) cycle
                  else ! f semicore
                     if(il==2.or.il>4) cycle
                  end if
               else if(iatomn_t<=86) then !81Tl~86Rn
                  if(iatomn_t-nint(ival(it))==80) then
                     if(il/=2) cycle
                  else ! s semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=88) then !87Fr~88Ra
                  if(iatomn_t-nint(ival(it))==86) then
                     if(il/=1) cycle
                  else ! sp semicore
                     if(il>2) cycle
                  end if
               else if(iatomn_t<=89) then !89Ac
                  if(iatomn_t-nint(ival(it))==86) then
                     if(il/=1.and.il/=3) cycle
                  else ! p semicore
                     if(il>3) cycle
                  end if
               else if(iatomn_t<=103) then !90Th~103Lr
                  if(il/=1.and.il/=3.and.il/=4) cycle
               else !104Rf~
                  if(iatomn_t-ival(it)==100) then
                     if(il/=1.and.il/=3) cycle
                  else ! f semicore
                     if(il==2.or.il>4) cycle
                  end if
               end if

               do t1 = 1, itau(kj,it)
                  irank(kj,t1,it) = t1
               end do
            end do
         end if
     end if

      if(is_gncpp == pp_PAW) then
         if(.not.paramset .and. ipaw(it) == ON) then
            phirpw(1:mmesh,1:nloc,1:ntau,it)=phir(1:mmesh,1:nloc,1:ntau)
            psirpw(1:mmesh,1:nloc,1:ntau,it)=psir(1:mmesh,1:nloc,1:ntau)
            wf_mnrc(it)=maxval(wf_nrc(:,:,it))
            if(iprippex >= 2) then
               write(nfout,'(" wf_mnrc(",i3," ) = ",i8,"<<rd_itau_etc_then_phir_chir>>")') it, wf_mnrc(it)
            end if
         end if
      end if

      if(paramset) then
!         if(it == 1) ntau = itau(1,it)
!
!         Do kj =2, nums_of_angmom_on_atomtype( it )
!            il = get_lp1_in_AngMomList( pot_has_soc(it),kj )
!
!            if(ntau < itau( kj,it)) ntau = itau(kj,it)
!         end do
         ntau = maxval(itau)
      end if

      if(ipri == PRINTLEVEL) &
           &  write(nfout,'(" end of <<rd_itau_etc_then_phir_chir_soc>>")')
    end subroutine rd_itau_etc_then_phir_chir_soc

! ================================================================================11.0

    subroutine set_irank_default( it )
      integer, intent(in) :: it

      integer :: il, t1, iatomn_t

      do il = 1, lpsmax(it)
         iatomn_t = nint(iatomn(it))
         if(iatomn_t<=4) then !1H~4Be
            if(il/=1) cycle
         else if(iatomn_t<=10) then !5B~10Ne
            if(il>2) cycle
         else if(iatomn_t<=12) then !11Na~12Mg
            if(iatomn_t-nint(ival(it))==10) then
               if(il/=1) cycle
            else ! sp semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=18) then !13Al~18Ar
            if(il>2) cycle
         else if(iatomn_t<=20) then !19K~20Ca
            if(iatomn_t-nint(ival(it))==18) then
               if(il/=1) cycle
            else ! sp semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=30) then !21Sc~30Zn
            if(iatomn_t-nint(ival(it))==18) then
               if(il/=1.and.il/=3) cycle
            else ! p semicore
               if(il>3) cycle
            endif
         else if(iatomn_t<=36) then !31Ga~36Kr
            if(iatomn_t-nint(ival(it))==28) then
               if(il>2) cycle
            else ! d semicore
               if(il>3) cycle
            end if
         else if(iatomn_t<=38) then !37Rb~38Sr
            if(iatomn_t-nint(ival(it))==36) then
               if(il/=1) cycle
            else ! sp semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=48) then !39Y~48Cd
            if(iatomn_t-nint(ival(it))==36) then
               if(il/=1.and.il/=3) cycle
            else ! p semicore
               if(il>3) cycle
            end if
         else if(iatomn_t<=54) then !49In~54Xe
            if(iatomn_t-nint(ival(it))==46) then
               if(il>2) cycle
            else ! d semicore
               if(il>3) cycle
            endif
         else if(iatomn_t<=56) then !55Cs~56Ba
            if(iatomn_t-nint(ival(it))==54) then
               if(il/=1) cycle
            else ! sp semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=57) then !57La
            if(iatomn_t-nint(ival(it))==54) then
               if(il/=1.and.il/=3) cycle
            else ! p semicore
               if(il>3) cycle
            end if
         else if(iatomn_t<=71) then !58Ce~71Lu
            if(il/=1.or.il/=3.or.il/=4) cycle
         else if(iatomn_t<=80) then !72Hf~80Hg
            if(iatomn_t-nint(ival(it))==58) then
               if(il/=1.and.il/=3) cycle
            else ! f semicore
               if(il==2.or.il>4) cycle
            end if
         else if(iatomn_t<=86) then !81Tl~86Rn
            if(iatomn_t-nint(ival(it))==80) then
               if(il/=2) cycle
            else ! s semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=88) then !87Fr~88Ra
            if(iatomn_t-nint(ival(it))==86) then
               if(il/=1) cycle
            else ! sp semicore
               if(il>2) cycle
            end if
         else if(iatomn_t<=89) then !89Ac
            if(iatomn_t-nint(ival(it))==86) then
               if(il/=1.and.il/=3) cycle
            else ! p semicore
               if(il>3) cycle
            end if
         else if(iatomn_t<=103) then !90Th~103Lr
            if(il/=1.and.il/=3.and.il/=4) cycle
         else !104Rf~
            if(iatomn_t-ival(it)==100) then
               if(il/=1.and.il/=3) cycle
            else ! f semicore
               if(il==2.or.il>4) cycle
            end if
         end if
         do t1 = 1, itau(il,it)
            irank(il,t1,it) = t1
            if ( sw_exclude_orb_tau_neq_1 == ON .and. t1 > 1 ) irank(il,t1,it) = 0
         end do

      end do

    end subroutine set_irank_default

    subroutine set_irank_minimal_pawpot( nfp, it, mode )
      integer, intent(in) :: nfp, it, mode

      integer, parameter :: max_orb = 10
      integer, save :: num_orb
      integer, allocatable, save :: flag(:), val_l(:), val_tau(:)
      character*2, allocatable, save :: val_nl(:)

      integer, parameter :: max_orb_ae = 25
      integer :: num_orb_ae
      character*2,  allocatable :: name_orb_ae(:)
      real(kind=DP), allocatable :: occ_orb_ae(:)

      integer :: istart, iend, len1, size1
      integer :: count(4), l1, t1, i, j

      character*132 :: char1, char2
      character*1 :: char3, char4
!
      if ( mode == 0 ) then
         if ( mype == 0 ) then
            if(.not.allocated(flag))    allocate( flag(max_orb) );    flag = 0
            if(.not.allocated(val_nl))  allocate( val_nl(max_orb) )
            if(.not.allocated(val_l))   allocate( val_l(max_orb) );   val_l = 0
            if(.not.allocated(val_tau)) allocate( val_tau(max_orb) ); val_tau = 0

            count = 0
            Do while ( .true. )
               read( nfp,'(A)',end=2) char1

               if ( char1(1:7) == "#   kc" ) then
                  len1 = len( trim(adjustl(char1)) )

                  num_orb = 0
                  loop1: Do while ( .true. )
                     istart = 14 +13*num_orb;    iend = istart +4
                     if ( istart > len1 ) exit loop1
                     char2(1:5) = char1(istart:iend)
                     char3 = char1(istart+3:istart+3)
                     char4 = char1(istart+1:istart+1)
                     num_orb = num_orb +1

                     val_nl(num_orb) = char1(istart:istart+1)

                     if ( char4 == "s" ) then
                        count(1) = count(1) +1
                        val_l(num_orb) = 1;   val_tau(num_orb) = count(1)
                     endif
                     if ( char4 == "p" ) then
                        count(2) = count(2) +1
                        val_l(num_orb) = 2;   val_tau(num_orb) = count(2)
                     endif
                     if ( char4 == "d" ) then
                        count(3) = count(3) +1
                        val_l(num_orb) = 3;   val_tau(num_orb) = count(3)
                     endif
                     if ( char4 == "f" ) then
                        count(4) = count(4) +1
                        val_l(num_orb) = 4;   val_tau(num_orb) = count(4)
                     endif

                     if ( sw_exclude_orb_tau_neq_1 == ON ) then
                        if ( char3 == "1" ) flag(num_orb) = 1
                     else
                        flag(num_orb) = 1
                     endif

                  End Do loop1
               endif
               if ( char1(1:8) == "#  state" ) then
                  len1 = len( trim(adjustl(char1)) )
                  Do i=1, num_orb
                     istart = 13 +13*( i-1 );  iend = istart +7
                     char2(1:8) = char1(istart:iend)
                     if ( char2(1:8) == "unbound" ) flag(i) = 0
                  End Do
                  goto 2
               endif
            End Do
2           continue
         endif

      else if ( mode == 1 ) then
         if ( mype == 0 ) then
            if ( sw_exclude_orb_unoccupied == ON ) then
               allocate( name_orb_ae( max_orb_ae ) )
               allocate( occ_orb_ae( max_orb_ae ) );  occ_orb_ae = 0.0d0
               call read_pp_tail_pawpot( nfp, it, max_orb_ae, num_orb_ae, &
                    &                    name_orb_ae, occ_orb_ae )

               Do i=1, num_orb
                  Do j=1, num_orb_ae
                     if ( val_nl(i) == name_orb_ae(j) ) then
                        if ( occ_orb_ae(j) == 0.0 ) flag(i) = 0
                     endif
                  End Do
               End Do
            endif

            Do i=1, num_orb
               l1 = val_l(i);  t1 = val_tau(i)
               if ( flag(i) > 0 ) irank(l1,t1,it) = t1
            End Do

            deallocate( flag );  deallocate( val_l );   deallocate( val_tau )
            deallocate( val_nl );
            if ( allocated( name_orb_ae) ) deallocate( name_orb_ae )
            if ( allocated( occ_orb_ae) )  deallocate( occ_orb_ae )
         endif

         size1 = sizeof( irank(:,:,it) ) /4
!!         size1 = kloc_t *ktau_t
         call mpi_bcast( irank(:,:,it), size1, mpi_integer, 0, MPI_CommGroup, ierr )

      endif
    end subroutine set_irank_minimal_pawpot


    subroutine read_pp_tail_pawpot( nfp, it, max_orb_ae, norb_ae, &
         &                          name_orb_ae, occ_orb_ae )
      integer, intent(in) :: nfp, it, max_orb_ae
      integer, intent(out) :: norb_ae
      character*2,   intent(out) :: name_orb_ae(max_orb_ae)
      real(kind=DP), intent(out) :: occ_orb_ae(max_orb_ae)

      integer :: num1, num2
      real(kind=DP) :: c1
      character(2) :: char2
      character(256) :: buf, buf2, buf3
      logical :: flag

      flag = .false.
      norb_ae = 0

      Do while ( .true. )
         read(nfp,'(A)',end=22) buf
         buf2 = trim(adjustl(buf))

         if ( buf2(1:16) == 'electron_config' ) then
            flag = .true.
            buf3(1:80) = buf2(17:96)
            read(buf3,*) norb_ae
         endif

         if ( flag ) then
            num2 = 0
            loop1: do while ( .true. )
               read(nfp,'(A)') buf
               if ( buf /= "" ) then
                  read(buf,*,err=2) char2, c1, num1
                  num2 = num2 +1
                  name_orb_ae(num2) = char2
                  occ_orb_ae(num2)  = c1
               endif
               if ( num2 == norb_ae ) then
                  flag = .false. ;  exit loop1
               endif
2              continue
            end do loop1
         endif
      End Do
22    continue
    end subroutine read_pp_tail_pawpot

   !========================================
    subroutine make_index_arrays_nlmt2l_m_t
   !========================================
      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2
      mm = 0
      do il1 = 1, lpsmax(it)
         if(il1 == iloc(it)) cycle
         do tau1 = 1, itau(il1,it)
            do im1 = 1, il1*2-1
               mm = mm + 1
               ltp( mm,it) = il1
               mtp( mm,it) = im1
               taup(mm,it) = tau1
            enddo
         enddo
      enddo
      ilmt(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp(n,it),mtp(n,it),taup(n,it),n=1,ilmt(it))
400      format((' !PP ',4i5))
      end if

      mm = 0
      do ilmt1 = 1, ilmt(it)
         mm = mm + 1
         if(.not.paramset) then
            index_lmt1_lmt2(mm,it,1) = ilmt1
            index_lmt1_lmt2(mm,it,2) = ilmt1
            w_non0_lmtxlmt(mm,it)    = 1
         endif
      enddo

      do ilmt1 = 1, ilmt(it)
         il1 =  ltp(ilmt1,it)
         im1 =  mtp(ilmt1,it)
         do ilmt2 = ilmt1+1, ilmt(it)
            il2 = ltp(ilmt2,it)
            im2 = mtp(ilmt2,it)
            if(il1 == il2 .and. im1 == im2) then
               mm = mm + 1
               if(.not.paramset) then
                  index_lmt1_lmt2(mm,it,1) = ilmt1
                  index_lmt1_lmt2(mm,it,2) = ilmt2
                  w_non0_lmtxlmt(mm,it)    = 2
               endif
            endif
         enddo
      enddo
      n_non0_lmtxlmt(it) = mm
      if(iprippex >= 2) write(nfout,'(" !PP #non0_elements = ",i8)') mm
    end subroutine make_index_arrays_nlmt2l_m_t

! ===================================== added by K. Tagami =================- 11.0
    subroutine mkindx_arrays_nlmt_2_j_l_m_t
      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2
      integer :: kj1
      integer :: ij1, ij2
      integer :: jval1

      mm = 0

      Do kj1 =1, nums_of_angmom_on_atomtype( it )
         il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
         jval1 = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

         if (il1 == iloc(it)) cycle

         do tau1 = 1, itau( kj1,it )
            do im1 = 1, il1*2-1
               mm = mm + 1
               ltp( mm,it ) = il1; mtp( mm,it ) = im1; taup( mm,it ) = tau1
               jtp( mm,it ) = jval1
               kjtp( mm,it ) = kj1

            enddo
         enddo
      End do

      ilmt(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp(n,it),mtp(n,it),taup(n,it),n=1,ilmt(it))
400      format((' !PP ',4i5))
      end if

      mm = 0
      do ilmt1 = 1, ilmt(it)
         mm = mm + 1
         if(.not.paramset) then
            index_lmt1_lmt2(mm,it,1) = ilmt1
            index_lmt1_lmt2(mm,it,2) = ilmt1
            w_non0_lmtxlmt(mm,it)    = 1
         endif
      enddo

      do ilmt1 = 1, ilmt(it)
         il1 = ltp(ilmt1,it)
         im1 = mtp(ilmt1,it)
         ij1 = jtp(ilmt1,it)

         do ilmt2 = ilmt1+1, ilmt(it)
            il2 = ltp(ilmt2,it)
            im2 = mtp(ilmt2,it)
            ij2 = jtp(ilmt2,it)

            if(il1 == il2 .and. im1 == im2 .and. ij1==ij2 ) then
               mm = mm + 1
               if(.not.paramset) then
                  index_lmt1_lmt2(mm,it,1) = ilmt1
                  index_lmt1_lmt2(mm,it,2) = ilmt2
                  w_non0_lmtxlmt(mm,it)    = 2
               endif
            endif
         enddo
      enddo
      n_non0_lmtxlmt(it) = mm

      if(iprippex >= 3) write(nfout,'(" !PP #non0_elements = ",i8)') mm
    end subroutine mkindx_arrays_nlmt_2_j_l_m_t
! ======================================================================== 11.0

   !========================================
    subroutine make_index_arrays_nlmt_add2lmt
   !========================================
      integer mm, il1, im1, n
      mm = 0
      do il1 = 1, lpsmax(it)
         if(il1 /= iloc(it)) cycle
         do im1 = 1, il1*2-1
            mm = mm + 1
            ltp_add( mm,it) = il1
            mtp_add( mm,it) = im1
         enddo
      enddo
      ilmt_add(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP TAUP")')
         write(nfout,400) (n,ltp_add(n,it),mtp_add(n,it),1,n=1,ilmt_add(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_add2lmt

! ===================================== added by K. Tagami ================ 11.0
    subroutine mkindx_arrays_nlmtadd_2_j_l_m_t
      integer mm, il1, im1, n
      integer :: ij1, kj1
      integer :: jval

      mm = 0

      Do kj1 = 1, nums_of_angmom_on_atomtype( it )
         il1  = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
         jval = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

         if (il1 /= iloc(it)) cycle

         do im1 = 1, il1*2-1
            mm = mm + 1
            ltp_add( mm,it ) = il1;  mtp_add( mm,it ) = im1
            jtp_add( mm,it ) = ij1
            kjtp_add( mm,it ) = kj1
         enddo
      End do

      ilmt_add(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP TAUP")')
         write(nfout,400) (n,ltp_add(n,it),mtp_add(n,it),1,n=1,ilmt_add(it))
400      format((' !PP ',4i5))
      end if

    end subroutine mkindx_arrays_nlmtadd_2_j_l_m_t
! ============================================================================ 11.0

    subroutine check_l_t_projectors
      integer :: i, ig, ip, jp, il1, tau1

      do i=1,num_projectors
         if (it/=proj_attribute(i)%ityp) cycle
         ig = proj_attribute(i)%group
         do ip=1,num_proj_elems(ig)
            jp = proj_group(ip,ig)
            il1 = proj_attribute(jp)%l+1
            tau1 = proj_attribute(jp)%t
            if ( il1 < 0 .or. il1 > lpsmax(it) ) then
               if ( mype == 0 ) then
                  write(*,'(A,I2,X,A,I2)') &
                       &  "l of projector no. ", i, "must be 0 <= l <=", lpsmax(it) -1
               endif
               call phase_error_with_msg(nfout,'invalid projector',__LINE__,__FILE__)
            endif
            if ( tau1 < 1 .or. tau1 > itau(il1,it) ) then
               if ( mype == 0 ) then
                  write(*,'(A,I2,X,A,I2)') &
                       &   "t of projector no. ", i, "must be 1 <= t <=", itau(il1,it)
               endif
               call phase_error_with_msg(nfout,'invalid projector',__LINE__,__FILE__)
            endif
         end do
      end do
    end subroutine check_l_t_projectors

   !==========================================
    subroutine make_index_arrays_nlmt_phi2lmt
   !==========================================
      integer mm, il1, tau1, im1, n, i
      integer :: ip_from_l(lpsmax(it)),ip,jp,ig

      do i=1,num_projectors
         if(it/=proj_attribute(i)%ityp) cycle
         ig = proj_attribute(i)%group
         do ip=1,num_proj_elems(ig)
            jp = proj_group(ip,ig)
            il1=proj_attribute(jp)%l+1
            ip_from_l(il1)=ip
         end do
      end do

      mm = 0
      L_loop: do il1 = 1, lpsmax(it)
         T_loop: do tau1 = 1, itau(il1,it)
            !!$do i=1,norbital(it)
            !!$   if(il1 == l_orb(i,it)+1 .and. tau1 == t_orb(i,it)) then
            !!$      do im1 = 1, il1*2-1
            !!$         mm = mm + 1
            !!$         ltp_phi( mm,it) = il1
            !!$         mtp_phi( mm,it) = im1
            !!$         taup_phi(mm,it) = tau1
            !!$      enddo
            !!$   end if
            !!$enddo
            do i=1,num_projectors
               if(it /= proj_attribute(i)%ityp) cycle
               if(il1 == proj_attribute(i)%l+1 .and. &
                & tau1 == proj_attribute(i)%t) then
                  do im1 = 1, il1*2-1
                     mm = mm + 1
                     ltp_phi( mm,it) = il1
                     mtp_phi( mm,it) = im1
                     taup_phi(mm,it) = tau1
                     iproj_phi(mm,it) = ip_from_l(il1)
                  end do
               end if
            end do
         enddo T_loop
      enddo L_loop
      ilmt_phi(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_phi(n,it),mtp_phi(n,it),taup_phi(n,it),n=1,ilmt_phi(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_phi2lmt

! ======================================= added by K. Tagami ============== 11.0
    subroutine mkindx_arrays_nlmtphi_2_j_l_m_t
      integer mm, il1, tau1, im1, n, i
      integer :: ip_from_l(lpsmax(it)),ip,jp,ig

      integer :: kj1
      real(kind=DP) :: jval

      do i=1,num_projectors
         if(it/=proj_attribute(i)%ityp) cycle
         ig = proj_attribute(i)%group
         do ip=1,num_proj_elems(ig)
            jp = proj_group(ip,ig)
            il1=proj_attribute(jp)%l+1
            ip_from_l(il1)=ip
! ------------------------------------
!  kt : under construction
!          ???? ij1 = proj_attribute(jp)%j +0.5
!          ??? ip_from_j( ij1 ) = ip
! --------------------------------

         end do
      end do

      mm = 0

      Do kj1=1, nums_of_angmom_on_atomtype( it )
         il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
         jval = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

         Do tau1 = 1, itau( kj1,it )
            if (it /= proj_attribute(i)%ityp) cycle
            if (il1 == proj_attribute(i)%l+1 .and. &
             & tau1 == proj_attribute(i)%t) then
               do im1 = 1, il1*2-1
                  mm = mm + 1
                  ltp_phi( mm,it) = il1;   mtp_phi( mm,it) = im1
                  taup_phi(mm,it) = tau1;
                  iproj_phi(mm,it) = ip_from_l(il1)
                                                 ! iproj_phi is for l, not for j.
	          jtp_phi( mm,it ) = jval
	          kjtp_phi( mm,it ) = kj1
               end do
            end if
         End do
      End do

      ilmt_phi(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_phi(n,it),mtp_phi(n,it),taup_phi(n,it),n=1,ilmt_phi(it))
400      format((' !PP ',4i5))
      end if

    end subroutine mkindx_arrays_nlmtphi_2_j_l_m_t
! =========================================================================== 11.0

   !========================================
    subroutine make_index_arrays_nlmt_pao2lmt
   !========================================
      integer mm, il1, tau1, im1, n
      mm = 0
      do il1 = 1, lpsmax(it)
         do tau1 = 1, itau(il1,it)
            if(irank(il1,tau1,it)<=0) cycle
            do im1 = 1, il1*2-1
               mm = mm + 1
               ltp_pao( mm,it) = il1
               mtp_pao( mm,it) = im1
               taup_pao(mm,it) = tau1
            enddo
         enddo
      enddo
      ilmt_pao(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_pao(n,it),mtp_pao(n,it),taup_pao(n,it),n=1,ilmt_pao(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_pao2lmt

! ===================================== added by K. Tagami ==================== 11.0
    subroutine mkindx_arrays_nlmtpao_2_j_l_m_t
      integer mm, il1, tau1, im1, n
      integer :: kj1
      real(kind=DP) :: jval

      mm = 0

      Do kj1=1, nums_of_angmom_on_atomtype( it )
         il1  = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
         jval = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

         Do tau1 = 1, itau( kj1,it )
           if(irank( kj1,tau1,it)<=0) cycle
           do im1 = 1, il1*2-1
              mm = mm + 1
              ltp_pao( mm,it) = il1
              mtp_pao( mm,it) = im1
              taup_pao(mm,it) = tau1
              jtp_pao(mm,it) =  jval
              kjtp_pao(mm,it) =  kj1
           enddo
         Enddo
      Enddo

      ilmt_pao(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_pao(n,it),mtp_pao(n,it),taup_pao(n,it),n=1,ilmt_pao(it))
400      format((' !PP ',4i5))
      end if

    end subroutine mkindx_arrays_nlmtpao_2_j_l_m_t
! ============================================================================ 11.0

   !==========================================
    subroutine cnstrct_of_localWF_Dion_and_q2
   !==========================================
      integer       :: il1,t1,t2,i,lmt1,il11,lmt2,im1,il22,im2,ip
      real(kind=DP) :: dion_tmp

      integer :: nrcproj,nrc_hubbard,ih,nrc
      real(kind=DP) :: window,prod,winfunc(mmesh),width, csum

      if(iprippex >= 2) then
         write(nfout,'(" !PP -- <<cnstrct_of_localWF_Dion_and_q2>> --")')
         write(nfout,'(" !PP it = ",i3)') it
         write(nfout,'(" !PP lpsmax ( ",i3,") = ", i3)') it,lpsmax(it)
         write(nfout,'(" !PP iloc   ( ",i3,") = ", i3)') it,iloc(it)
         write(nfout,'(" !PP nmesh  ( ",i3,") = ", i5)') it,nmesh(it)
         write(nfout,'(" !PP wos(100:101) = ",2d16.8)') wos(100),wos(101)
      end if

      dion(:,:,it) = 0.d0; q(:,:,it) = 0.d0; if_q_isnotzero(:,:,it) = NO
      loop_il1 :do il1 = 1, lpsmax(it)
         if(iprippex >= 2) then
            write(nfout,'(" !PP il1 = ",i5)') il1
            write(nfout,'(" !PP itau(",i3,",",i3,") = ",i3)') il1,it,itau(il1,it)
         end if
         if(iloc(it) == il1) cycle
         do t1 = 1, itau(il1,it)
            if(iprippex >= 2) then
               write(nfout,'(" !PP phir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, phir(100,il1,t1),phir(101,il1,t1)
               write(nfout,'(" !PP chir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, chir(100,il1,t1),chir(101,il1,t1)
            end if
            do t2 = 1,itau(il1,it)
               s(t1,t2) = 0.d0
               do i = 1, nmesh(it)
                  s(t1,t2) = s(t1,t2) + wos(i)*phir(i,il1,t1)*chir(i,il1,t2)
               enddo
            enddo
         enddo
         if(iprippex >= 2) then
            write(nfout,'(" !PP")')
            do t1 = 1, itau(il1,it)
               do t2 = 1, itau(il1,it)
                  write(nfout,'(" !PP B[nm] il=",i3," : s (",i3,",",i3,") = ",d16.8)') il1,t1,t2,s(t1,t2)
               end do
            end do
            write(nfout,'(" !PP")')
            do t1 = 1, itau(il1,it)
               do t2 = 1, itau(il1,it)
                  write(nfout,'(" !PP Q[nm] il=",i3," : q (",i3,",",i3,") = ",d16.8)') il1,t1,t2,qij(t1,t2,il1)
               end do
            end do
            write(nfout,'(" !PP")')
            do t1 = 1, itau(il1,it)
               do t2 = 1, itau(il1,it)
                  write(nfout,'(" !PP QV[nm] il=",i3," : qv (",i3,",",i3,") = ",d16.8)') il1,t1,t2,qvij(t1,t2,il1)
               end do
            end do
         end if
         if(iprippex >= 3) write(nfout,'(" !PP -- after t1 loop")')

         do lmt1 = 1, ilmt(it)
            il11 = ltp(lmt1,it)
            im1  = mtp(lmt1,it)
            t1  = taup(lmt1,it)
            do lmt2 = 1, ilmt(it)
               il22 =  ltp(lmt2,it)
               im2  =  mtp(lmt2,it)
               t2   = taup(lmt2,it)
               if(il11==il1 .and. il22==il1 .and. im1==im2) then
                  q   (lmt1,lmt2,it) = qij(t1,t2,il1)
                  dion(lmt1,lmt2,it) = s(t1,t2) - qvij(t1,t2,il1) &
                       & + eps(il1,t2,it)*qij(t1,t2,il1)
               endif
               if(dabs(q(lmt1,lmt2,it)) > DELTA) if_q_isnotzero(lmt1,lmt2,it) = YES
            enddo
         enddo
         if(iprippex >= 3) write(nfout,'(" !PP -- after lmt1 loop")')

         do ip = 1, n_non0_lmtxlmt(it)
            lmt1 = index_lmt1_lmt2(ip,it,1)
            il11 =  ltp(lmt1,it)
            if(il11 == il1) then
               im1  =  mtp(lmt1,it)
               t1   = taup(lmt1,it)
               lmt2 = index_lmt1_lmt2(ip,it,2)
               t2  =  taup(lmt2,it)
               q_indp(ip,it) = qij(t1,t2,il1)
               dion_indp(ip,it) = s(t1,t2) - qvij(t1,t2,il1) &
                    & + eps(il1,t2,it)*qij(t1,t2,il1)
            endif
         enddo

         call matrix_inversion(nfout,ntau,itau(il1,it),s,sinv)  ! s(i,j) = <phir(i)|chir(j)>

         do t1 = 1, itau(il1,it)
            betar(1:nmesh(it),il1,t1,it) = 0.d0
            do t2 = 1, itau(il1,it)
               do i= 1, nmesh(it)
                  betar(i,il1,t1,it) = betar(i,il1,t1,it) &
                       & + sinv(t2,t1)*chir(i,il1,t2)
               enddo
            enddo
         end do
      end do loop_il1

      if(is_gncpp == pp_PAW .or. is_gncpp == pp_GNCPP2) then
      if(sw_use_add_proj == ON) then
         loop_il1_2 :do il1 = 1, lpsmax(it)
            if(iprippex >= 2) write(nfout,'(" !PP il1 = ",i5)') il1
            if(iloc(it) /= il1) cycle
            betar_add(1:nmesh(it),it) = 0.d0
            do i=1,nmesh(it)
               if(radr(i) > rcproj(it)) then
                  nrcproj = i
                  exit
               end if
            end do
            betar_add(1:mmesh,it) = 0.d0
            do i=1,nrcproj
               window = 1.d0 - exp(-kappa(it)*(radr(nrcproj)-radr(i))**2)
               betar_add(i,it) = window*phir(i,il1,1)
            enddo
            prod = 0.d0
            do i=1,nmesh(it)
               prod = prod + wos(i)*betar_add(i,it)*phir(i,il1,1)
            end do
! === DEBUG by tkato 2014/03/04 ================================================
            if(abs(prod) > SmallestPositiveNumber) then
! ==============================================================================
            do i=1,nmesh(it)
               betar_add(i,it) = betar_add(i,it)/prod
            end do
! === DEBUG by tkato 2014/03/04 ================================================
            end if
! ==============================================================================
! debug
           if(iprippex >= 2) then
               write(nfout,*) ' !PP Additioanl projector: it=',it,'il=',il1
               do i=1,min(nrcproj+100,nmesh(it)),40
                  write(nfout,*) ' !PP ',i,radr(i),betar_add(i,it)
               end do
           end if
! end debug
         end do loop_il1_2
      end if
      end if

      do lmt2 = 2, ilmt(it)
         do lmt1 = 1, lmt2-1
            dion_tmp = (dion(lmt1,lmt2,it) + dion(lmt2,lmt1,it))*0.5d0
            dion(lmt1,lmt2,it) = dion_tmp
            dion(lmt2,lmt1,it) = dion_tmp
         enddo
      enddo

!!$      if(is_gncpp /= pp_PAW) then
      if(sw_orb_popu == ON) then
         ! for PDOS
         do lmt1 = 1, ilmt_phi(it)
            t1  = taup_phi(lmt1,it)
            il1 = ltp_phi(lmt1,it)
            q_phi(lmt1,it) = qij(t1,t1,il1)
         end do
      end if
!!$      end if

      if ( is_gncpp == pp_PAW ) then
         if ( sw_orb_popu == ON .and. orb_popu_method == 2 ) then
            do il1=1, lpsmax(it)
               do t1=1,itau(il1,it)
                  do t2=1,itau(il1,it)
                     csum = 0.0d0

                     if ( ipaw(it) == ON ) then
                        Do i=1, nmesh(it)
                           csum = csum +wos(i) *phirt(i,il1,t1,it) &
                                &              *( psir(i,il1,t2) -phir(i,il1,t2) )
                        End Do
                     else
                        Do i=1, nmesh(it)
                           csum = csum +wos(i) *phirt(i,il1,t1,it) &
                                &              *( phir_ae(i,il1,t2) -phir(i,il1,t2) )
                        End Do
                     endif

                     q_phirt_pw(il1,t1,t2,it) = csum
                  End do
               End do
            End do
            q_phi = 0.0d0       ! q_phi is not necessary when orb_popu_method == 2
         endif
      endif

      if( (is_gncpp == pp_GNCPP2 .or. is_gncpp == pp_PAW) .and. num_projectors>0) then
         do ih=1,num_projectors
            if(it /= proj_attribute(ih)%ityp) cycle
            do i=nmesh(it)-1,1,-1
               if(proj_attribute(ih)%radius > radr(i)) then
                  nrc_hubbard=i+1
                  exit
               end if
            end do
            winfunc(1:nmesh(it)) = 0.d0

!  ================================ modified by K. Tagami =================== 5.0
!            if(projector_type == ATOMIC_ORBITAL) then
!               width = proj_attribute(ih)%fwhm*sqrt(log(2.d0))
!               do i=1,nrc_hubbard
!                  winfunc(i) = (1.d0 - exp(-((radr(nrc_hubbard)-radr(i))/width)**2))/(1.d0 - exp(-(radr(nrc_hubbard)/width)**2))
!               end do
!            else if(projector_type == SPHERICAL_HARMONICS ) then
!               do i=1,nrc_hubbard
!                  winfunc(i) = 1.d0
!               end do
!            end if

            if ( sw_orbital_cut_smooth == ON ) then
               width = proj_attribute(ih)%fwhm*sqrt(log(2.d0))
               do i=1,nrc_hubbard
                  winfunc(i) = (1.d0 - exp(-((radr(nrc_hubbard)-radr(i))/width)**2))/(1.d0 - exp(-(radr(nrc_hubbard)/width)**2))
               end do
            else if(projector_type == SPHERICAL_HARMONICS ) then
               do i=1,nrc_hubbard
                  winfunc(i) = 1.d0
               end do
            end if
! ===================================================================== 5.0

            il1 = proj_attribute(ih)%l+1
            if(is_gncpp == pp_PAW .and. flg_paw) then
               do t1=1,itau(il1,it)
                  do t2=1,itau(il1,it)
                     nrc = min(wf_nrc(il1,t1,it),wf_nrc(il1,t2,it))
                     if(proj_attribute(ih)%radius_was_defined) nrc = nrc_hubbard
                     prodphi(ih,t1,t2) = 0.d0
                     !!$do i=1,min(wf_nrc(il1,t1,it),wf_nrc(il1,t2,it))
                     do i=1,nrc
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                             & + wos(i)*psir(i,il1,t1)*psir(i,il1,t2)
                     end do
                  end do
               end do
               if(.not.proj_attribute(ih)%radius_was_defined) nrc_hubbard = maxval(wf_nrc(il1,:,it))
            else if(ae_wavefunctions_are_detected(it)) then
               do t1=1,itau(il1,it)
                  do t2=1,itau(il1,it)
                     prodphi(ih,t1,t2) = 0.d0
                     do i=1,nrc_hubbard
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                     & + wos(i)*winfunc(i)*phir_ae(i,il1,t1)*phir_ae(i,il1,t2)
                     end do
                  end do
               end do
            else
               do t1=1,itau(il1,it)
                  do t2=1,itau(il1,it)
                     prodphi(ih,t1,t2) = qij(t1,t2,il1)
                     do i=1,nrc_hubbard
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                     & + wos(i)*winfunc(i)*phir(i,il1,t1)*phir(i,il1,t2)
                     end do
                  end do
               end do
            end if
            if(iprippex >= 2) then
               if(ae_wavefunctions_are_detected(it)) then
                  write(nfout,'(" !PP ae_wavefunctions_are_detected(it) = .true.")')
                  dion_tmp = 0.d0
! === DEBUG by tkato 2011/09/22 ================================================
                  do t1=1,itau(il1,it)
! ==============================================================================
                  do t2=1,itau(il1,it)
                     do i = 1, nrc_hubbard
                        dion_tmp = dion_tmp + phir_ae(i,il1,t1)*phir_ae(i,il1,t2)
                     end do
                  end do
! === DEBUG by tkato 2011/09/22 ================================================
                  end do
! ==============================================================================
                  write(nfout,'(" sum of phir_ae = ", d20.8)') dion_tmp
               else
                  write(nfout,'(" !PP ae_wavefunctions_are_detected(it) = .false.")')
               end if
               write(nfout,'(" nrc_hubbard = ",i8)') nrc_hubbard
               dion_tmp = 0.d0
               do i = 1, nrc_hubbard
                  dion_tmp = dion_tmp + winfunc(i)
               end do
               write(nfout,'(" sum of winfunc = ", d20.8)') dion_tmp
               write(nfout,'(" !PP Projector: ih = ",i3)') ih
               write(nfout,'(" !PP Products of phi: L=",i3," RC=",f10.5)') &
                  & il1-1, radr(nrc_hubbard)
               do t1=1,itau(il1,it)
                  do t2=1,itau(il1,it)
                     write(nfout,'(" !PP t1,t2,prod=",2i3,f10.5)') t1,t2,prodphi(ih,t1,t2)
                  end do
               end do
            end if
!debug
!!$            do t1=1,itau(il1,it)
!!$               do t2=1,itau(il1,it)
!!$                  prod = 0.d0
!!$                  do i=1,nmesh(it)
!!$                     prod = prod + wos(i)*betar(i,il1,t1,it)*phir(i,il1,t2)
!!$                  end do
!!$                  write(nfout,'(" !PP t1,t2,bp=",2i3,f10.5)') t1,t2,prod
!!$               end do
!!$            end do
!end debug
         end do
      end if

      if(ipripp >= 2) call wd_index_arrays_etc(nfout,it) ! printout dion, q, index_lmt1_lmt2 etc
    end subroutine cnstrct_of_localWF_Dion_and_q2

! =============================== added by K. Tagami =================== 11.0
    subroutine cnstrct_LocalWF_Dion_and_q_soc

      integer       :: il1,t1,t2,i,lmt1,il11,lmt2,im1,il22,im2,ip
      integer :: kj1

      integer :: ij11, ij22
      integer :: kjtmp, iltmp
      real(kind=DP) :: dion_tmp
      integer :: jval1, jvaltmp

      integer :: nrcproj,nrc_hubbard,ih
      real(kind=DP) :: window,prod,winfunc(mmesh),width, csum

      if(iprippex >= 2) then
         write(nfout,'(" !PP -- <<cnstrct_localWF_Dion_and_q_soc>> --")')
         write(nfout,'(" !PP it = ",i3)') it
         write(nfout,'(" !PP lpsmax ( ",i3,") = ", i3)') it,lpsmax(it)
         write(nfout,'(" !PP iloc   ( ",i3,") = ", i3)') it,iloc(it)
         write(nfout,'(" !PP nmesh  ( ",i3,") = ", i5)') it,nmesh(it)
         write(nfout,'(" !PP wos(100:101) = ",2d16.8)') wos(100),wos(101)
      end if

      dion(:,:,it) = 0.d0; q(:,:,it) = 0.d0; if_q_isnotzero(:,:,it) = NO

      loop_kj1 : Do kj1 =1, nums_of_angmom_on_atomtype( it )
         il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
         jval1 = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

         if(iprippex >= 2) then
            write(nfout,'(" !PP il1 = ",i5)') il1
            write(nfout,'(" !PP itau(",i3,",",i3,") = ",i3)') il1,it,itau(il1,it)
         end if
         if(iloc(it) == il1) cycle

         do t1 = 1, itau(kj1,it)
            if(iprippex >= 2) then
               write(nfout,'(" !PP phir(100:101,",i3,",",i3,") = ",2d16.8)') kj1,t1, phir(100,kj1,t1),phir(101,kj1,t1)
               write(nfout,'(" !PP chir(100:101,",i3,",",i3,") = ",2d16.8)') kj1,t1, chir(100,kj1,t1),chir(101,kj1,t1)
            end if
            do t2 = 1,itau(kj1,it)
               s(t1,t2) = 0.d0
               do i = 1, nmesh(it)
                  s(t1,t2) = s(t1,t2) + wos(i)*phir(i,kj1,t1)*chir(i,kj1,t2)
               enddo
            enddo
         enddo
!--------------------------------- ????????????????????? --- sij = dinagonal ?? --
! phir() ha up , down component wo motsunoka ?
! j=l+1/2 ha up,  l-1/2 ha down componet nomi o motsu??
! ----------
         if(iprippex >= 2) then
            write(nfout,'(" !PP")')
            do t1 = 1, itau(kj1,it)
               do t2 = 1, itau(kj1,it)
                  write(nfout,'(" !PP B[nm] kj=",i3," : s (",i3,",",i3,") = ",d16.8)') kj1,t1,t2,s(t1,t2)
               end do
            end do

            write(nfout,'(" !PP")')

            do t1 = 1, itau(kj1,it)
               do t2 = 1, itau(kj1,it)
                  write(nfout,'(" !PP Q[nm] kj=",i3," : q (",i3,",",i3,") = ",d16.8)') kj1,t1,t2,qij(t1,t2,kj1)
               end do
            end do
            write(nfout,'(" !PP")')
            do t1 = 1, itau(kj1,it)
               do t2 = 1, itau(kj1,it)
                  write(nfout,'(" !PP QV[nm] kj=",i3," : qv (",i3,",",i3,") = ",d16.8)') kj1,t1,t2,qvij(t1,t2,kj1)
               end do
            end do
         end if

         if(iprippex >= 3) write(nfout,'(" !PP -- after t1 loop")')

         do lmt1 = 1, ilmt(it)
            il11 = ltp(lmt1,it)
            im1  = mtp(lmt1,it)
            t1  = taup(lmt1,it)
            ij11 = jtp(lmt1,it)

            do lmt2 = 1, ilmt(it)
               il22 =  ltp(lmt2,it)
               im2  =  mtp(lmt2,it)
               t2   = taup(lmt2,it)
               ij22 = jtp(lmt2,it)

               if(il11==il1 .and. il22==il1 .and. im1==im2) then
                  if ( ij11==jval1 .and. ij22==jval1 ) then

                     q   (lmt1,lmt2,it) = qij(t1,t2,kj1)
                     dion(lmt1,lmt2,it) = s(t1,t2) - qvij(t1,t2,kj1) &
                          & + eps(kj1,t2,it)*qij(t1,t2,kj1)

                  endif
               endif
               if(dabs(q(lmt1,lmt2,it)) > DELTA) if_q_isnotzero(lmt1,lmt2,it) = YES
            enddo
         enddo

         if(iprippex >= 3) write(nfout,'(" !PP -- after lmt1 loop")')

         do ip = 1, n_non0_lmtxlmt(it)
            lmt1 = index_lmt1_lmt2(ip,it,1)
            il11 =  ltp(lmt1,it)
            ij11 = jtp(lmt1,it)

            if(il11 == il1 .and. ij11 == jval1 ) then
               im1  =  mtp(lmt1,it)
               t1   = taup(lmt1,it)
               lmt2 = index_lmt1_lmt2(ip,it,2)
               t2  =  taup(lmt2,it)
               q_indp(ip,it) = qij(t1,t2,kj1)
               dion_indp(ip,it) = s(t1,t2) - qvij(t1,t2,kj1) &
                    & + eps(kj1,t2,it)*qij(t1,t2,kj1)
            endif
         enddo

         call matrix_inversion(nfout,ntau,itau(kj1,it),s,sinv)
                                     ! s(i,j) = <phir(i)|chir(j)>

         do t1 = 1, itau(kj1,it)
            betar(1:nmesh(it),kj1,t1,it) = 0.d0
            do t2 = 1, itau(kj1,it)
               do i= 1, nmesh(it)
                  betar(i,kj1,t1,it) = betar(i,kj1,t1,it) &
                       & + sinv(t2,t1)*chir(i,kj1,t2)
               enddo
            enddo
         end do
      end do loop_kj1

!!!      if (is_gncpp == pp_GNCPP2) then
      if ( is_gncpp == pp_GNCPP2 .or. is_gncpp == pp_PAW ) then

         if(sw_use_add_proj == ON) then
            loop_kj1_2 : Do kj1 =1, nums_of_angmom_on_atomtype( it )
               il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
               jval1 = get_jph_in_AngMomList( pot_has_soc(it),kj1 )

               if(iprippex >= 2) write(nfout,'(" !PP kj1 = ",i5)') kj1
               if(iloc(it) /= il1) cycle

               betar_add(1:nmesh(it),it) = 0.d0
               do i=1,nmesh(it)
                  if(radr(i) > rcproj(it)) then
                     nrcproj = i
                     exit
                  end if
               end do

               betar_add(1:mmesh,it) = 0.d0

               do i=1,nrcproj
                  window = 1.d0 - exp(-kappa(it)*(radr(nrcproj)-radr(i))**2)
                  betar_add(i,it) = window*phir(i,kj1,1)
               enddo
               prod = 0.d0
               do i=1,nmesh(it)
                  prod = prod + wos(i)*betar_add(i,it)*phir(i,kj1,1)
               end do
               do i=1,nmesh(it)
                  betar_add(i,it) = betar_add(i,it)/prod
               end do
! debug
               if(iprippex >= 2) then
                  write(nfout,*) ' !PP Additioanl projector: it=',it,'kj=',kj1
                  do i=1,min(nrcproj+100,nmesh(it)),40
                     write(nfout,*) ' !PP ',i,radr(i),betar_add(i,it)
                  end do
               end if
! end debug
            end do loop_kj1_2
         end if
      end if

      do lmt2 = 2, ilmt(it)
         do lmt1 = 1, lmt2-1
            dion_tmp = (dion(lmt1,lmt2,it) + dion(lmt2,lmt1,it))*0.5d0
            dion(lmt1,lmt2,it) = dion_tmp
            dion(lmt2,lmt1,it) = dion_tmp
         enddo
      enddo

!!$      if(is_gncpp /= pp_PAW) then
      if(sw_orb_popu == ON) then
         ! for PDOS
         do lmt1 = 1, ilmt_phi(it)
            t1  = taup_phi(lmt1,it)
            il1 = ltp_phi(lmt1,it)
            kj1 = kjtp_phi(lmt1,it)
            q_phi(lmt1,it) = qij(t1,t1,kj1)
         end do
      end if
!!$      end if

      if( (is_gncpp == pp_GNCPP2 .or. is_gncpp == pp_PAW) .and. num_projectors>0) then
         do ih=1,num_projectors
            if(it /= proj_attribute(ih)%ityp) cycle
            do i=nmesh(it)-1,1,-1
               if(proj_attribute(ih)%radius > radr(i)) then
                  nrc_hubbard=i+1
                  exit
               end if
            end do

            winfunc(1:nmesh(it)) = 0.d0

!  ================================ modified by K. Tagami =================== 5.0
!            if(projector_type == ATOMIC_ORBITAL) then
!               width = proj_attribute(ih)%fwhm*sqrt(log(2.d0))
!               do i=1,nrc_hubbard
!                  winfunc(i) = (1.d0 - exp(-((radr(nrc_hubbard)-radr(i))/width)**2))/(1.d0 - exp(-(radr(nrc_hubbard)/width)**2))
!               end do
!            else if(projector_type == SPHERICAL_HARMONICS ) then
!               do i=1,nrc_hubbard
!                  winfunc(i) = 1.d0
!               end do
!            end if

            if ( sw_orbital_cut_smooth == ON ) then
               width = proj_attribute(ih)%fwhm*sqrt(log(2.d0))
               do i=1,nrc_hubbard
                  winfunc(i) = (1.d0 - exp(-((radr(nrc_hubbard)-radr(i))/width)**2))/(1.d0 - exp(-(radr(nrc_hubbard)/width)**2))
               end do
            else if(projector_type == SPHERICAL_HARMONICS ) then
               do i=1,nrc_hubbard
                  winfunc(i) = 1.d0
               end do
            end if
! ===================================================================== 5.0

            il1 = proj_attribute(ih)%l+1
!!
!!
! ------------------------ Not supported .. ------------------
!!!!!!!!!!!!!!            jval1 = nint( proj_attribute(ih)%j +0.5 )
            call phase_error_with_msg(nfout,"kt: Sorry",__LINE__,__FILE__)
!!
!!
!!
! -----------------------------------------------------------
            kj1 = 0
            Do kjtmp=1, nums_of_angmom_on_atomtype( it )
               iltmp = get_lp1_in_AngMomList( pot_has_soc(it),kjtmp )
               jvaltmp = get_jph_in_AngMomList( pot_has_soc(it),kjtmp )
               if ( il1 == iltmp .and. jval1 == jvaltmp ) then
                  kj1 = kjtmp
               endif
            End do
! -------------------------------

            if(is_gncpp == pp_PAW .and. flg_paw) then
               do t1=1,itau(kj1,it)
                  do t2=1,itau(kj1,it)
                     prodphi(ih,t1,t2) = 0.d0
                     do i=1,min(wf_nrc(kj1,t1,it),wf_nrc(kj1,t2,it))
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                             & + wos(i)*psir(i,kj1,t1)*psir(i,kj1,t2)
                     end do
                  end do
               end do
               nrc_hubbard = maxval(wf_nrc(il1,:,it))

            else if(ae_wavefunctions_are_detected(it)) then
               do t1=1,itau(kj1,it)
                  do t2=1,itau(kj1,it)
                     prodphi(ih,t1,t2) = 0.d0
                     do i=1,nrc_hubbard
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                     & + wos(i)*winfunc(i)*phir_ae(i,kj1,t1)*phir_ae(i,kj1,t2)
                     end do
                  end do
               end do
            else
               do t1=1,itau(kj1,it)
                  do t2=1,itau(kj1,it)
                     prodphi(ih,t1,t2) = qij(t1,t2,kj1)
                     do i=1,nrc_hubbard
                        prodphi(ih,t1,t2) = prodphi(ih,t1,t2) &
                     & + wos(i)*winfunc(i)*phir(i,kj1,t1)*phir(i,kj1,t2)
                     end do
                  end do
               end do
            end if

            if(iprippex >= 2) then
               if(ae_wavefunctions_are_detected(it)) then
                  write(nfout,'(" !PP ae_wavefunctions_are_detected(it) = .true.")')
                  dion_tmp = 0.d0
                  do t2=1,itau(kj1,it)
                     do i = 1, nrc_hubbard
                        dion_tmp = dion_tmp + phir_ae(i,kj1,t1)*phir_ae(i,kj1,t2)
                     end do
                  end do
                  write(nfout,'(" sum of phir_ae = ", d20.8)') dion_tmp
               else
                  write(nfout,'(" !PP ae_wavefunctions_are_detected(it) = .false.")')
               end if
               write(nfout,'(" nrc_hubbard = ",i8)') nrc_hubbard

               dion_tmp = 0.d0
               do i = 1, nrc_hubbard
                  dion_tmp = dion_tmp + winfunc(i)
               end do
               write(nfout,'(" sum of winfunc = ", d20.8)') dion_tmp

               write(nfout,'(" !PP Projector: ih = ",i3)') ih
               write(nfout,'(" !PP Products of phi: L=",i3," RC=",f10.5)') &
                  & il1-1, radr(nrc_hubbard)
               do t1=1,itau(kj1,it)
                  do t2=1,itau(kj1,it)
                     write(nfout,'(" !PP t1,t2,prod=",2i3,f10.5)') t1,t2,prodphi(ih,t1,t2)
                  end do
               end do
            end if

         end do
      end if

      if(ipripp >= 2) call wd_index_arrays_etc(nfout,it)
                            ! printout dion, q, index_lmt1_lmt2 etc
    end subroutine cnstrct_LocalWF_Dion_and_q_soc
! ======================================================================= 11.0

   !==============================================
    subroutine coulomb_potential_in_Rspace(nsize)
   !==============================================
      integer, intent(in) :: nsize
      real(kind=DP),allocatable,dimension(:) :: da, db ! d(nsize)
      real(kind=DP)       :: s2, rhs, rhs1, bm
      integer             :: i
     !+++++++++++++++++++++++++++++++
      allocate(da(nsize)); da = 0.d0
      allocate(db(nsize)); db = 0.d0
     !+++++++++++++++++++++++++++++++
!!$ASASASASAS
      if(rhvr(2).eq.0 .or. rhvr(1).eq.0 .or. h(it).eq.0)then
        deallocate(da)
        deallocate(db)
        return
      endif
!!$ASASASASAS
      s2 = dlog(rhvr(2)/rhvr(1))/h(it)
      rhs = rhvr(1)
      wkx(1)  = rhs*radr(1)/(s2+1)
      wky(1)  = rhs/s2
      db(1)   = h(it)*rhs*3.d0
      da(1)   = db(1)*radr(1)
      rhs1    = rhs
      do i = 2,3
         rhs     = rhvr(i)
         wkx(i)  = wkx(i-1)+h(it)*(rhs*radr(i)+rhs1*radr(i-1))*0.5d0
         wky(i)  = wky(i-1)+h(it)*(rhs        +rhs1          )*0.5d0
         db(i)   = h(it) *rhs*3.d0
         da(i)   = db(i)*radr(i)
         rhs1    = rhs
      enddo
      do i = 4,mesh_t
         rhs    = rhvr(i)
         db(4)  = h(it)*rhs*3.d0
         da(4)  = db(4)*radr(i)
         wkx(i)=(9*wkx(i-1)-wkx(i-3)+da(4)+2.d0*da(3)-da(2))/8.d0
         wky(i)=(9*wky(i-1)-wky(i-3)+db(4)+2.d0*db(3)-db(2))/8.d0
         da(1)  = da(2)
         db(1)  = db(2)
         da(2)  = da(3)
         db(2)  = db(3)
         da(3)  = da(4)
         db(3)  = db(4)
      enddo
      bm               = wky(mesh_t)
     !C--*--COULOMB POTENTIAL RVC
      vvv = 0.d0
      do i = 1,mesh_t
         vvv(i) = wkx(i)+radr(i)*(bm-wky(i))
      enddo
     !+++++++++++++++++++++++++++++++
      deallocate(da); deallocate(db)
     !+++++++++++++++++++++++++++++++
    end subroutine coulomb_potential_in_Rspace

   !=========================================
    subroutine atomic_xc_potential_in_Rspace   ! (ptxc)
   !=========================================
      real(kind=DP) :: chgrsq, vxc1, vxc2, rs, zet, yc
      integer       :: jsm = 1  ! 1: non-mag, 2: magnetic
      real(kind=DP), parameter :: delta40 = 1.d-40
      integer       :: i, ncut

#ifdef GGA_CHECK
      if(iprippex >= 2) write(nfout,'(" !PP -- atomic_xc_potential_in_Rspace -- ")')
#endif
      wkx = 0.d0
      if(itpcc(it) == 0) then
         wkx(1:mesh_t) = rhvr(1:mesh_t)
      else if(itpcc(it) == 1) then
         wkx(1:mesh_t) = rhvr(1:mesh_t) + rhpcr(1:mesh_t)
      endif
      if(what_is_ps_xctype() == GGA) then ! -(m_PseudoPotential)
            call gga_grad_alloc()
            call get_gradient_of_rho(jsm,mddiff,mmesh,mesh_t,wkx,radr&
                 &,h(it),fdiff,coeff1,coeff2,rho,grdnts) ! -(b_PP)  -> rho, grdnts
#ifdef GGA_CHECK
            if(iprippex >= 2) then
               write(nfout,'(" !PP after get_gradient_of_rho ")')
               write(nfout,'(" !PP rad, rho, grdnts(*,1), grdnts(*,2)")')
               do i = 1, mesh_t
                  write(nfout,'(" !PP ( ",i4,") rad, rho, grad, grad2 = ",4d20.8," (get_gradient)")') &
                       &        i, radr(i), rho(i), grdnts(i,1), grdnts(i,2)
               end do
            end if
#endif

            if(         ps_xctype == 'ldapw91' .or. ps_xctype == 'LDAPW91' &
                 & .or. ps_xctype == 'ldapbek' .or. ps_xctype == 'LDAPBEK' ) grdnts = 0.d0

            call gga_xc_potential_atomic(jsm,len_ps_xctype,ps_xctype &
                 &,mmesh,mesh_t,rho,radr,grdnts,fdiff,vxc,eps_chg,eps_grad)
            call gga_grad_dealloc()
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
            do i = mmesh,1,-1
               if(wkx(i)/(PAI4*radr(i)**2) > eps_chg) then
                  ncut = i
                  goto 1001
               end if
            end do
1001        continue
            do i = ncut+1, mmesh
               vxc(i) = 0.d0
            end do
#endif

#ifdef GGA_CHECK
            if(iprippex >= 2) then
               write(nfout,'(" !PP ps_xctype, xctype = ",a7,2x,a7)') ps_xctype, xctype
               do i = 1, nmesh(it)
                  write(nfout,'(" !PP vxc(",i4,") = ",d18.10, " radr(",i4,") = ",d18.10)') i, vxc(i),i, radr(i)
               end do
            end if
#endif
      else
         vxc = 0.d0
         do i = 1, mesh_t
            vxc1 = 0.d0; vxc2 = 0.d0
            chgrsq = wkx(i)
            if(dabs(chgrsq) > delta40) then
               rs = (3*radr(i)*radr(i)/chgrsq)**(1.d0/3.d0)
               zet = 0.d0
          !--*--EXCHANGE-CORRELATION POTENTIALS ARE GIVEN IN RYD. UNITS.
               call xcpot(ps_xctype,vxc1,vxc2,rs,zet,yc)
            endif
            vxc(i) = radr(i)*vxc1*0.5d0
         enddo
      end if
    end subroutine atomic_xc_potential_in_Rspace

   !==================================
    subroutine vlocr_plus_hartree_xc2
   !==================================
      integer :: meshup,i

#ifdef GGA_CHECK
      if(iprippex >=2) write(nfout,*) ' !PP << vlocr_plus_hartree_xc2 >>'
#endif

      if(nmesh(it) > mesh_t) then
         meshup = mesh_t
      else
         meshup = nmesh(it)
      endif
#ifdef _DEBUG_WRITE_
      if(iprippex >= 2) then
         write(nfout,'(" !PP --- vlocr --- ")')
         write(nfout,'(" !PP ",3d26.18)') (vlocr(i),i=1,nmesh(it))
         write(nfout,'(" !PP -- vlocr before adding vcoulomb and xc potential --")')
      end if
#endif
      call sum_of_vlocr(it,nfout)
      call sum_of_vvv_and_vxc(it,nfout) ! Now vvv is the Hartree potential

      if(is_gncpp == pp_GNCPP1) then
         do i = 1, meshup
            vlocr(i) = vlocr(i) - (vvv(i)+vxc(i))/radr(i)
         enddo
      else if(is_gncpp==pp_GNCPP2 .or. is_gncpp==pp_PAW) then
         do i = 1, nmesh(it)
            vlocr(i) = vlocr2(i)
         end do
      end if

#ifdef GGA_CHECK
      if(iprippex >= 2) then
         do i = meshup, meshup-100, -1
            write(nfout,'(" !PP Vlocr, a3, a2: ",i6,3f12.6)') i, vlocr(i),vvv(i),vxc(i)
         end do
      end if
#endif

      if(iprippex >= 3) then
         write(nfout,'(" !PP --- after - hartree - xc --")')
         write(nfout,'(" !PP ",3d30.20)') (vlocr(i),i=1,nmesh(it))
         write(nfout,'(" !PP --- coulomb ---------------")')
         write(nfout,'(" !PP ",6d15.7)') (vvv(i),i=1,nmesh(it))
         write(nfout,'(" !PP --- xc      ---------------")')
         write(nfout,'(" !PP ",6d15.7)') (vxc(i),i=1,nmesh(it))
         write(nfout,'(" !PP --- rhvr    ---------------")')
         write(nfout,'(" !PP ",6d15.7)') (rhvr(i),i=1,nmesh(it))
      end if
    end subroutine vlocr_plus_hartree_xc2

   !=======================================
    subroutine read_natomn_ival_iloc_itpcc
   !=======================================
!!$      integer, parameter  ::    len_str = 80
!!$      character(len=len_str) :: str
      real(kind=DP) :: natomn
!!$      logical :: comment_statement

      if(mype == 0) then
         rewind(nfp)
         if(iprippex >= 2) write(nfout,'(" !PP  READING POTENTIAL FILE ",i5)') nfp
         call skip_commentlines_pp(nfp)
!!$         comment_statement = .true.
!!$         do while(comment_statement)
!!$            read(nfp,'(a80)') str
!!$            if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
!!$                 & .or. str(1:1) == '%') then
!!$               if(iprippex >= 1) write(nfout,'(a80)') str
!!$            else
!!$               comment_statement = .false.
!!$            endif
!!$         enddo
!!$         read(str,*)      natomn,ival(it),iloc(it),itpcc(it)
         read(nfp,*)      natomn,ival(it),iloc(it),itpcc(it)
         if(iprippex >= 2) write(nfout,110) natomn,ival(it),iloc(it),itpcc(it)
110      FORMAT(' !PP ',2f8.4,2I4,'  : NATOMN, IVAL, ILOC, ITPCC ')
! ==== KT_add ==== 2014/07/27
         if ( itpcc(it)==2 ) then
            itpcc(it) = 0
!            if ( iprippex >=2 ) then
               write(nfout,*) "!** itpcc is changed from 2  to  0  for it = ", it
!            endif
         endif
! ================ 2014/07/27
      end if
      if(npes > 1 ) then
         call mpi_bcast(ival(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(iloc(it), 1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(itpcc(it),1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(iatomn(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(natomn,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      end if
      if( abs(iatomn(it) - natomn) > 1.d-8)  then
         if(ipripp >= 2) write(nfout,*) ' !PP iatomn.ne.natomn ',iatomn(it),natomn
         if(iatomn(it)>=1) call phase_execution_error(INVALID_ATOMIC_NUMBER)
         if(ipri >= 1) write(nfout,'(a,i4,a,f7.3)') ' !PP atomic number for element',it,':',natomn
         iatomn(it) = natomn
      end if
      if(iprippex >= 2) then
         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: ival(",i4,") = ",f8.4)') it,ival(it)
         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: iloc(",i4,") = ",i3)') it,iloc(it)
         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: itpcc_(",i3,") = ",i3)') it,itpcc(it)
      end if

    end  subroutine read_natomn_ival_iloc_itpcc

   !============================
    subroutine read_ps_xctype()
   !============================
      integer  ::                     counter = 0
      character(len=len_ps_xctype) :: ps_xctype0
      logical :: tf

! ========================================= modified by K. Tagami ============= 11.0
!      if(iprippex >= 2) write(nfout,'(" !PP read_ps_xctype: MPI_CommGroup = ",i8)') MPI_CommGroup
!
      if (iprippex >= 3) then
         write(nfout,*) " !PP read_ps_xctype: MPI_CommGroup = ", MPI_CommGroup
      endif
! ============================================================================= 11.0

      counter = counter + 1
      if(mype == 0) then   ! MPI
         call skip_commentlines_pp(nfp)
         read(nfp,'(a7)') ps_xctype0
         call nameck(ps_xctype0,len_ps_xctype)
         if(iprippex >=1) write(nfout,120) ps_xctype0
120      FORMAT(' !PP ',A7,'   : NAME ')
         if(counter == 1) then
            ps_xctype = ps_xctype0
         else
            call strncmp0(ps_xctype,ps_xctype0,tf)
            if(.not.tf) then
!!$            if(ps_xctype /= ps_xctype0) then
               if(iprippex >= 1) write(nfout,*) ' !PP XC TYPE OF ', it,'TH ATOM IS NOT CORRECT'
               call phase_error_with_msg(nfout,'invalid xc type',__LINE__,__FILE__)
            endif
            call strncmp0(xctype,ps_xctype0,tf)
            if(.not.tf) then
!!$            if(xctype /= ps_xctype0) then
               if(iprippex >= 1) write(nfout,*) &
                    &' !PP xc-type of ',it,'th atom is not equal to ' &
                    &,'that described in the input file'
               if(iprippex >= 1) write(nfout,'(" !PP xctype = ",a7," ps_xctype0 = ",a7)') xctype,ps_xctype0
            endif
         endif
      end if
      call mpi_bcast(ps_xctype,len_ps_xctype,mpi_character,0,MPI_CommGroup,ierr)
      if(counter == 1 .and. xctype == "nogiven") call m_CtrlP_set_xctype(ps_xctype,nfout)
      if(iprippex >= 1) write(nfout,'(" !PP read_ps_xctype: ps_xctype = ",a7)') ps_xctype
    end subroutine read_ps_xctype

   !=======================
    subroutine read_alp_cc
   !=======================
      if(mype == 0) then          ! MPI
         read(nfp,*) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
         if(iprippex >= 1) write(nfout,130) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
130      FORMAT(' !PP ',4F12.6,'  :   ALP,CC')
      end if
      call mpi_bcast(alp(1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(alp(2,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(cc(1,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(cc(2,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
    end subroutine read_alp_cc

   !==============================
    subroutine read_nmesh_xh_rmax
   !==============================
      if(mype == 0) then  ! MPI
         read(nfp,*) nmesh(it),xh(it),rmax(it)
         if(iprippex >=1) write(nfout,140) nmesh(it),xh(it),rmax(it)
140      format(' !PP ',i4,2f12.6,'  :   nmesh,  xh, rmax')
      end if
      call mpi_bcast(nmesh(it),1,mpi_integer         ,0,MPI_CommGroup,ierr)
      call mpi_bcast(xh(it),   1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(rmax(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
      if(iprippex>=2) then
         write(nfout,'(" !PP read_nmesh_xc_rmax: nmesh(it) = ",i8)') nmesh(it)
         write(nfout,'(" !PP read_nmesh_xc_rmax: xh(it)    = ",f10.6)') xh(it)
         write(nfout,'(" !PP read_nmesh_xc_rmax: rmax(it)  = ",f10.6)') rmax(it)
      end if

    end subroutine read_nmesh_xh_rmax

   !============================
    subroutine determine_lpsmax
   !============================
      character(len=7) :: lpsmaxj
      if(mype == 0) then
         if(iloc(it) == 4) then
            lpsmax(it) = iloc(it)
         else
            lpsmax(it) = 3
         endif
         read(nfp,'(a7)') lpsmaxj
         if(lpsmaxj == 'F-STATE') then
            lpsmax(it) = 4
            if(iprippex >=1) write(6,'(" !PP lpsmaxj == F-STATE")')
         else if(lpsmaxj == 'G-STATE') then
            lpsmax(it) = 5
            if(iprippex >=1) write(nfout,'(" !PP lpsmaxj == G-STATE")')
         else if(lpsmaxj == 'ZEROPOT') then
            lpsmax(it) = 0
            if(iprippex >= 1) write(nfout, '(" !PP lpsmaxj == DUMMY ATOM")')
         else
            backspace nfp
         endif
      end if
      call mpi_bcast(lpsmax(it),1,mpi_integer,0,MPI_CommGroup,ierr)
    end subroutine determine_lpsmax

! ===================================== added by K. Tagami ================= 11.0
    subroutine calc_nums_angmom_on_atomtype
                      ! This gives the number of angular momentum
                      ! l in the case of non-relativisitic pseudopot.
                      ! j in the case of fully relativistic pseudopot.
      integer :: il1, mm

      mm = 0
!
      if ( pot_has_soc(it) ) then
         Do il1=1, lpsmax(it)
            if ( il1 == 1 ) then              ! s-orbital
               mm = mm + 1
            else
               mm = mm + 2
            endif
         End do
      else
         Do il1=1, lpsmax(it)
            mm = mm + 1
         End do
      endif

      nums_of_angmom_on_atomtype(it) = mm

    end subroutine calc_nums_angmom_on_atomtype

! ========================================================================== 11.0

   !==============================
    logical function check_vall()
   !==============================
      character(len=4) :: vall
      read(nfp,'(a4)') vall
      if(vall == 'VALL') then
         check_vall = .true.
      else
         backspace nfp
         check_vall = .false.
      end if
    end function check_vall

#ifdef _POT_SMOOTHING_
   !=============================
    subroutine smoothing_vlocr()
   !=============================
      real(kind=DP), parameter :: delta7 = 1.d-7, delta10 = 1.d-10
      integer       :: icheck, i
      real(kind=DP) :: s
      icheck = 0
      do i = nmesh(it),1,-1
         if((rhvr(i)/pai4/radr(i)**2) < delta7) then
            if(dabs(vlocr(i)+ival(it)/radr(i)) > delta10 ) then
               if(iprippex >= 2) then
                  write(nfout,'(" !PP VLOC:",i6,f12.6,3d16.8)') &
                       & i,radr(i),vlocr(i), -ival(it)/radr(i)&
                       &,          vlocr(i) + ival(it)/radr(i)
               end if
               vlocr(i) = -ival(it)/radr(i)
            endif
         else
            icheck = icheck + 1
         endif
         if(icheck > 30) exit
      enddo
      s = 0.d0
      do i = 1, nmesh(it)
         s = s + wos(i)*vlocr(i)*radr(i)**2
      end do
      if(iprippex >= 2) write(nfout,'(" !PP (smoothing) s = ",d20.8)') s
    end subroutine smoothing_vlocr
#endif
!....................................................................

! =========================== added by K. Tagami ==================== 11.0
    subroutine make_SOC_strength_Zeff_nonpaw

      real(kind=DP) :: HyperFineConst

      integer :: ir, ier
      integer :: nrc, il, it1, it2
      real(kind=DP) :: fac1, fac2,  sum

      real(kind=DP), allocatable :: wght(:)

      HyperFineConst = 1.0d0 / InvHyperFineConst
      fac1 = 0.5d0 * HyperFineConst**2

      nrc = nmesh(it)

      allocate( wght(1:nrc) ); wght = 0.0d0
      call set_weight_exp( ier, 1, nrc, radr, wght )

      Do il=1, lpsmax(it)
         if  ( il == iloc(it) ) cycle

         Do it1=1,itau(il,it)
            Do it2=1,itau(il,it)
               if ( ae_wavefunctions_are_detected(it) ) then

                  sum=0.d0
                  if ( SpinOrbit_MassCorrection == 0 ) then
                     Do ir = 1, nrc
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2)
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum

                  else if ( SpinOrbit_MassCorrection == 1 ) then
                     Do ir = 1, nrc
                        fac2 = 1.0d0 + HyperFineConst**2 &
                             &        *dble(iatomn(it)) /radr(ir)
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2) /fac2
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum

                  else if ( SpinOrbit_MassCorrection == 2 ) then
                     Do ir = 1, nrc
                        fac2 = 1.0d0 + HyperFineConst**2 /2.0d0 &
                             &        *dble(iatomn(it)) /radr(ir)
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2) /fac2
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum
                  endif

               else

                  write(nfout,*) '! *************  *********** '
                  write(nfout,*) '! ** SOC strength cannot be calculated when AE wfns are not found in pseudopotential.'
                  call phase_error_with_msg(nfout,&
                  'SOC strength cannot be calculated when AE wfns are not found in pseudopotential.' &
                  ,__LINE__,__FILE__)

               endif
            End do
         End do
      End do

      deallocate(wght)

    end subroutine make_SOC_strength_Zeff_nonpaw
! ========================================================================= 11.0

   !================================================
    subroutine cnstrct_of_PiPjPlPm_k_Qijk_etc
   !================================================
        integer:: il1,il2,il3,il4,tau1,tau2,tau3,tau4,lks,lkl,mm
        integer:: t2min,lpmx,t3min,t4min,l4min
        integer:: ilk1,ilk2
!        logical:: a,b
        integer:: ilt1,ilt2,ilt3,ilt4,lt4min
!
!        lpmx=lpsmax(it)
!        mm=0
!        do il1=1,lpmx
!            if(iloc(it) == il1) cycle
!            do tau1=1,itau(il1,it)
!                do il2=il1,lpmx
!                    if(iloc(it) == il2) cycle
!                    t2min=1
!                    if(il1 == il2) t2min=tau1
!                    do tau2=t2min,itau(il2,it)
!
!                        do ilk1=abs(il1-il2),il1+il2-2,2
!                            if(ilk1 > lcmax) cycle
!                            do il3=il1,lpmx
!                                if(iloc(it) == il3) cycle
!                                t3min=1
!                                if(il3 == il1) t3min=tau1
!                                do tau3=t3min,itau(il3,it)
!                                    l4min=il3
!                                    if(il1 == il3 .and. tau1 == tau3) l4min=max(il2,il3)
!                                    do il4=l4min,lpmx
!                                        if(iloc(it) == il4) cycle
!                                        a=(il3.eq.il4)
!                                        b=(il1.eq.il3).and.(tau1.eq.tau3).and.(il2.eq.il4)
!                                        t4min=1
!                                        if(a .and. .not.b) then
!                                            t4min=tau3
!                                        else if(.not.a .and. b) then
!                                            t4min=tau2
!                                        else if(a .and. b) then
!                                            t4min=max(tau3,tau2)
!                                        end if
!                                        do tau4=t4min,itau(il4,it)
!
!                                            do ilk2=abs(il3-il4),il3+il4-2,2
!                                                if(ilk2 > lcmax) cycle
!                                                if(ilk1 == ilk2) then
!
!print '(5i4)',il1,il2,il3,il4,ilk1
!                                                    mm=mm+1
!                                                end if
!                                            end do
!
!                                        end do
!                                    end do
!                                end do
!                            end do
!
!                        end do
!
!                    end do
!                end do
!            end do
!        end do
!
!
!10 continue

        mm=0
        do ilt1=1,iltpw(it)
            il1=lppw(ilt1,it)
            do ilt2=ilt1,iltpw(it)
                il2=lppw(ilt2,it)
                do ilk1=abs(il1-il2),il1+il2-2,2
                    if(ilk1 > lcmax) cycle
                    do ilt3=ilt1,iltpw(it)
                        il3=lppw(ilt3,it)
                        lt4min=ilt3
                        if(ilt1==ilt3) lt4min=max(ilt2,ilt3)
                        do ilt4=lt4min,iltpw(it)
                            il4=lppw(ilt4,it)
                            do ilk2=abs(il3-il4),il3+il4-2,2
                                if(ilk2 > lcmax) cycle
                                if(ilk1 == ilk2) then
    !print '(5i4)',il1,il2,il3,il4,ilk1
                                    mm=mm+1
                                    if(.not.paramset) then
                                        ipppp(ilt1,ilt2,ilt3,ilt4,ilk1+1,it)=&
                                                                    mm+mmpppp
                                    ! call makeppppps(ilt1,ilt2,ilt3,ilt4,ilk1,mm+mpppp)
                                        call make_VPiPjPlPm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,0)
!                                        call make_VPiPjPlPm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,0)
                                        call make_VPiPjPlPm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,1)
!                                        call make_VPiPjPlPm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,1)
                                        call make_VQijQlm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!                                        call make_VQijQlm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
                                        call make_VQijPlPm_kSym &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!                                        call make_VQijPlPm_kSym_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!print *,ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do

        mmpppp=mmpppp+mm
        if(paramset) npppp=mmpppp
        if(.not.paramset) then
           if(iprippex >=2) then
              write(nfout,'(" npppp = ",i5)') npppp
           end if
        end if

!if(.not.paramset) then
!print *,'Here'
!    stop
!end if

    end subroutine cnstrct_of_PiPjPlPm_k_Qijk_etc


   !======================================
    subroutine make_index_arrays_nlt2l_t
   !======================================
      integer mm, il1, tau1, n, nn, im1
      mm = 0
      nn = 0
      do il1 = 1, lpsmax(it)
         if(il1 == iloc(it)) cycle
         do tau1 = 1, itau(il1,it)
           mm = mm + 1
           lppw(mm,it) = il1
           tppw(mm,it) = tau1
            if(.not.paramset) then
               do im1=1,il1*2-1
                    nn=nn+1
                    index_lmt2lt(nn,it)=mm
               end do
            end if
         enddo
      enddo
      iltpw(it) = mm
      if(mm.gt.nltpw .and. paramset) nltpw=mm
      if(ipri >= 2) then
         write(nfout,*) ' ILT, LPPW, TPPW  '
         write(nfout,400) (n,lppw(n,it),tppw(n,it),n=1,iltpw(it))
400      format((' ',4i5))
      end if

    end subroutine make_index_arrays_nlt2l_t

   !=============================================================
   subroutine make_VPiPjPlPm_k(ilt1,ilt2,ilt3,ilt4,ilk,mp,ae_or_ps)
   !=============================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp
        integer,intent(in):: ae_or_ps                      ! 0--AE   1--PS

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: pij,plm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        allocate(pij(nnrc))
        allocate(plm(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)* &
                        psirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)* &
                        psirpw(1:nnrc,il4,itau4,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                        phirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                        phirpw(1:nnrc,il4,itau4,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        end if

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*pr3(jr)*pr4(jr)+radr(jr+1)**ilk*pr3(jr+1)*pr4(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*pr3(jr)*pr4(jr)+ &
!                            radr(jr+1)**(-1-ilk)*pr3(jr+1)*pr4(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

!        deallocate(pr1,pr2,pr3,pr4,dsum)
        deallocate(pij,plm,dsum,wght)
        return

   end subroutine make_VPiPjPlPm_k

   !=====================================================================
   subroutine make_VPiPjPlPm_k_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp,ae_or_ps)
   !=====================================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp
        integer,intent(in):: ae_or_ps                      ! 0--AE   1--PS

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: pij,plm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        allocate(pij(nnrc))
        allocate(plm(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)* &
                        psirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc2)=psirpw(1:nnrc2,il3,itau3,it)* &
                         psirpw(1:nnrc2,il4,itau4,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                        phirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc2)=phirpw(1:nnrc2,il3,itau3,it)* &
                         phirpw(1:nnrc2,il4,itau4,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        end if

! ===== Symmetrization ===== !

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)* &
                        psirpw(1:nnrc,il4,itau4,it)
            plm(1:nnrc2)=psirpw(1:nnrc2,il1,itau1,it)* &
                         psirpw(1:nnrc2,il2,itau2,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                        phirpw(1:nnrc,il4,itau4,it)
            plm(1:nnrc2)=phirpw(1:nnrc2,il1,itau1,it)* &
                         phirpw(1:nnrc2,il2,itau2,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=vaeijlm_k(mp)+sum0*PAI4/(2*dble(ilk)+1)
            vaeijlm_k(mp)=vaeijlm_k(mp)*0.5d0
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=vpsijlm_k(mp)+sum0*PAI4/(2*dble(ilk)+1)
            vpsijlm_k(mp)=vpsijlm_k(mp)*0.5d0
        end if


!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*pr3(jr)*pr4(jr)+radr(jr+1)**ilk*pr3(jr+1)*pr4(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*pr3(jr)*pr4(jr)+ &
!                            radr(jr+1)**(-1-ilk)*pr3(jr+1)*pr4(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

!        deallocate(pr1,pr2,pr3,pr4,dsum)
        deallocate(pij,plm,dsum,wght)
        return

   end subroutine make_VPiPjPlPm_k_dbg


   !=====================================================
   subroutine make_VQijQlm_k(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !=====================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,qlm_k,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij,iqlm

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        iqlm=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if
        if(iqlm.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if

        allocate(qij_k(nnrc))
        allocate(qlm_k(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijqlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qlm_k(jr)+radr(jr+1)**ilk*qlm_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qlm_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qlm_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

        deallocate(qij_k,qlm_k,dsum,wght)
        return

   end subroutine make_VQijQlm_k

   !==========================================================
   subroutine make_VQijQlm_k_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==========================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,qlm_k,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij,iqlm

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        iqlm=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if
        if(iqlm.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if

        allocate(qij_k(nnrc2))
        allocate(qlm_k(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        qij_k(1:nnrc2)=qrspspw(1:nnrc2,iqij)
        qlm_k(1:nnrc2)=qrspspw(1:nnrc2,iqlm)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijqlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qlm_k(jr)+radr(jr+1)**ilk*qlm_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qlm_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qlm_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

        deallocate(qij_k,qlm_k,dsum,wght)
        return

   end subroutine make_VQijQlm_k_dbg


   !==========================================================
   subroutine make_VQijPlPm_kSym(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==========================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,plpm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        allocate(qij_k(nnrc))
        allocate(plpm(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            vqijplpm_ks(mp)=0.d0
            goto 10
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                    phirpw(1:nnrc,il4,itau4,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=sum0*PAI4/(2*dble(ilk)+1)

10 continue

        iqij=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            goto 20
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                    phirpw(1:nnrc,il2,itau2,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qij_k(jr)+radr(jr+1)**ilk*qij_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qij_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qij_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

20      continue

        deallocate(qij_k,plpm,dsum,wght)
        return

   end subroutine make_VQijPlPm_kSym

   !==============================================================
   subroutine make_VQijPlPm_kSym_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==============================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,plpm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        allocate(qij_k(nnrc2))
        allocate(plpm(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            vqijplpm_ks(mp)=0.d0
            goto 10
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc2)=phirpw(1:nnrc2,il3,itau3,it)* &
                    phirpw(1:nnrc2,il4,itau4,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=sum0*PAI4/(2*dble(ilk)+1)

! ******* Symmetrization *******
        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)


10 continue

        iqij=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            goto 20
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc2)=phirpw(1:nnrc2,il1,itau1,it)* &
                    phirpw(1:nnrc2,il2,itau2,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

! ***** Symmetrization *****

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

        vqijplpm_ks(mp)=vqijplpm_ks(mp)*0.5d0

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qij_k(jr)+radr(jr+1)**ilk*qij_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qij_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qij_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

20      continue

        deallocate(qij_k,plpm,dsum,wght)
        return

   end subroutine make_VQijPlPm_kSym_dbg

   !=============================================
    subroutine cnstrct_of_VHij_VHpsij_Kinij
   !=============================================
        integer:: il,it1,it2,mp

        if(iprippex>=2)then
        write(nfout,'("--cnstrct_of_VHij_VHpsij_Kinij--")')
        write(nfout,'("    kin_ae_ps     vionaeij     vionpsij     vionpsqij")')
        endif

        do il=1,lpsmax(it)
            if(iloc(it)==il) cycle
            do it1=1,itau(il,it)
                do it2=it1,itau(il,it)
!                    call make_kin_ae_ps(il,it1,it2)
                    call make_kin_ae_ps2(il,it1,it2)
                    call make_vionaeij(il,it1,it2)
                    call make_vionpsij(il,it1,it2)
                    call make_vionpsqij(il,it1,it2)

                    if(iprippex.gt.2)then
                    write(nfout,'(4f15.7)') kin_ae_psij(il,it1,it2,it), &
                                            vionaeij(il,it1,it2,it), &
                                            vionpsij(il,it1,it2,it), &
                                            vionpsqij(il,it1,it2,it)
                    endif
                end do
            end do
        end do
        return
    end subroutine cnstrct_of_VHij_VHpsij_Kinij


   !====================================
    subroutine cnstrct_of_dion_kin_ion
   !====================================
        integer:: lmt1,il1,im1,it1
        integer:: lmt2,il2,im2,it2


        do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it)
            im1=mtp(lmt1,it)
            it1=taup(lmt1,it)
            do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                if(il1==il2 .and. im1==im2) then
                    dion_kin_ion(lmt1,lmt2,it)= &
                                    kin_ae_psij(il1,it1,it2,it)+ &
                                    vionaeij(il1,it1,it2,it)- &
                                    vionpsij(il1,it1,it2,it)- &
                                    vionpsqij(il1,it1,it2,it)
                end if
            end do
        end do

        do lmt1=2,ilmt(it)
            do lmt2=1,lmt1-1
                dion_kin_ion(lmt1,lmt2,it)=dion_kin_ion(lmt2,lmt1,it)
            end do
        end do

        if(iprippex>2)then
        write(nfout,*) ' -- dion_kin_ion ---'
        do lmt1 = 1, ilmt(it)
            write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
              &               ,(dion_kin_ion(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
        enddo
        endif
!do lmt1=1,ilmt(it)
!print '(9f8.5)' ,(dion_kin_ion(lmt1,lmt2,it),lmt2=1,ilmt(it))
!end do
!stop
        return
    end subroutine cnstrct_of_dion_kin_ion


   !=======================================
    subroutine make_kin_ae_ps(il,it1,it2)
   !=======================================
        integer,intent(in):: il,it1,it2

        integer:: nrc,ier,i
        real(DP),pointer,dimension(:):: pai,paj,psi,psj
        real(DP),pointer,dimension(:):: dpai,dpaj,dpsi,dpsj
        real(DP),pointer,dimension(:):: wght
        real(DP):: sum,rr,ll
! real(DP),pointer,dimension(:):: ddpaj,ddpsj
!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(pai(1:nrc))
        allocate(paj(1:nrc))
        allocate(psi(1:nrc))
        allocate(psj(1:nrc))
        allocate(dpai(1:nrc))
        allocate(dpaj(1:nrc))
        allocate(dpsi(1:nrc))
        allocate(dpsj(1:nrc))
        allocate(wght(1:nrc))
!allocate(ddpaj(1:nrc))
!allocate(ddpsj(1:nrc))
        pai(1:nrc)=psirpw(1:nrc,il,it1,it)
        paj(1:nrc)=psirpw(1:nrc,il,it2,it)
        psi(1:nrc)=phirpw(1:nrc,il,it1,it)
        psj(1:nrc)=phirpw(1:nrc,il,it2,it)

        call calc_diff_exp(ier,4,nrc,radr,pai,dpai)
        call calc_diff_exp(ier,4,nrc,radr,paj,dpaj)
        call calc_diff_exp(ier,4,nrc,radr,psi,dpsi)
        call calc_diff_exp(ier,4,nrc,radr,psj,dpsj)
!call calc_ddiff_exp(ier,5,nrc,radr,paj,dpaj,ddpaj)
!call calc_ddiff_exp(ier,5,nrc,radr,psj,dpsj,ddpsj)
        sum=0.d0
        ll=dble(il-1)
        call set_weight_exp(ier,1,nrc,radr,wght)
        do i=1,nrc
            rr=radr(i)
            sum=sum+(dpai(i)*dpaj(i)+ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
                    -dpsi(i)*dpsj(i)-ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wght(i)
        end  do
        kin_ae_psij(il,it1,it2,it)=sum/2.d0

!print *,sum/2.d0
!sum=0.d0
!ll=dble(il-1)
!call set_weight_exp(ier,1,nrc,radr,wght)
!do i=1,nrc
!    rr=radr(i)
!    sum=sum-0.5d0*(pai(i)*ddpaj(i)-ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
!            -psi(i)*ddpsj(i)+ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wght(i)
!end  do
!print *,sum
!stop

        deallocate(pai,paj,psi,psj)
        deallocate(dpai,dpaj,dpsi,dpsj)
        deallocate(wght)
        return
   end subroutine make_kin_ae_ps

   !========================================
    subroutine make_kin_ae_ps2(il,it1,it2)
   !========================================
        integer,intent(in):: il,it1,it2

        integer:: nrc,ier,ir
        real(DP),pointer,dimension(:):: wght1,wght2
        real(DP):: tmp1,tmp2
        real(DP):: vloc_ps,vloc_ae,emhnhm,emsnsm,ene,bmt

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(wght1(1:nrc))
        allocate(wght2(1:nmesh(it)))

        call set_weight_exp(ier,1,nrc,radr,wght1)
        call set_weight_exp(ier,1,nmesh(it),radr,wght2)
        vloc_ps=0.d0
        vloc_ae=0.d0
        emhnhm=0.d0
        emsnsm=0.d0
        do ir=1,nrc
            tmp1=wght1(ir)*phirpw(ir,il,it1,it)*phirpw(ir,il,it2,it)
            tmp2=wght1(ir)*psirpw(ir,il,it1,it)*psirpw(ir,il,it2,it)
            vloc_ps=vloc_ps+tmp1*vloc_scr_ps(ir)
            vloc_ae=vloc_ae+tmp2*vloc_scr_ae(ir)
            emhnhm=emhnhm+tmp1
            emsnsm=emsnsm+tmp2
        end  do
        if(it1.eq.it2) then
            ene=eps(il,it2,it)
        else
            ene=0.5d0*(eps(il,it1,it)+eps(il,it2,it))
        end if
        emhnhm=emhnhm*ene
        emsnsm=emsnsm*ene
        bmt=0.d0
        do ir=1,nmesh(it)
            bmt=bmt+wght2(ir)*phirpw(ir,il,it1,it)*chir(ir,il,it2)
        end do
        if(it1.ne.it2) then
            do ir=1,nmesh(it)
                bmt=bmt+wght2(ir)*phirpw(ir,il,it2,it)*chir(ir,il,it1)
            end do
            bmt=bmt*0.5d0
        end if
        kin_ae_psij(il,it1,it2,it)=emsnsm-vloc_ae-emhnhm+bmt+vloc_ps

!print *,sum/2.d0
!sum=0.d0
!ll=dble(il-1)
!call set_weight_exp(ier,1,nrc,radr,wos)
!do i=1,nrc
!    rr=radr(i)
!    sum=sum-0.5d0*(pai(i)*ddpaj(i)-ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
!            -psi(i)*ddpsj(i)+ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wos(i)
!end  do
!print *,sum
!stop
        deallocate(wght1,wght2)
        return
   end subroutine make_kin_ae_ps2

   !======================================
    subroutine make_vionaeij(il,it1,it2)
   !======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,nrc2,ier,i
        real(DP),pointer,dimension(:):: pipj,dsum,wght
        real(DP):: sum,rr,sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,is

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        nrc2=nmesh(it)
        allocate(pipj(1:nrc))
        allocate(dsum(1:nrc))
        allocate(wght(1:nrc2))
        pipj(1:nrc)=psirpw(1:nrc,il,it1,it)* &
                    psirpw(1:nrc,il,it2,it)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+pipj(i)/radr(i)*wght(i)
        end do
        vionaeij(il,it1,it2,it)=-dble(iatomn(it))*sum

        dsum=0.d0
  !rhcorpw(i,it)
        do ir=1,nrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+rhcorpw(i0+j*is,it)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+rhcorpw(jr,it)*wght(jr)
                end do
            end if
            sum1=sum1*pipj(ir)/radr(ir)
            if(ir==nrc2) then
                sum2=0.d0
            else if((ir<=nrc2-1).and.(ir>=nrc2-4)) then
                do ii=ir,nrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(1.d0/radr(i0+j*is))* &
                                rhcorpw(i0+j*is,it)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nrc2,radr,wght)
                do jr=ir,nrc2
                    sum2=sum2+(1.d0/radr(jr))* &
                                rhcorpw(jr,it)*wght(jr)
                end do
            end if
            sum2=sum2*pipj(ir)
            dsum(ir)=sum1+sum2
        end do

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum0=0.d0
        do ir=1,nrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vionaeij(il,it1,it2,it)=vionaeij(il,it1,it2,it)+sum0

        deallocate(pipj,dsum,wght)
        return
   end subroutine make_vionaeij

   !======================================
    subroutine make_vionpsij(il,it1,it2)
   !======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,ier,i
        real(DP),pointer,dimension(:):: pipj,wght
        real(DP):: sum,rr

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(pipj(1:nrc))
        allocate(wght(1:nrc))
        pipj(1:nrc)=phirpw(1:nrc,il,it1,it)* &
                    phirpw(1:nrc,il,it2,it)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+pipj(i)*vlocr2(i)*wght(i)
        end do
        vionpsij(il,it1,it2,it)=sum

        deallocate(pipj,wght)
        return
   end subroutine make_vionpsij
     !vionpsqij     ! d(nloc,ntau,ntau,ntyp)
   !=======================================
    subroutine make_vionpsqij(il,it1,it2)
   !=======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,ier,i,iq
        real(DP),pointer,dimension(:):: qij,wght
        real(DP):: sum,rr

        iq=iqitg(il,it1,il,it2,1,it)
        if(iq.eq.0) then
            !!$print *,'il1,itau1,il2,itau2,il3 : ',il,it1,il,it2,1
            !!$print *,' is not set in iqitg ! make_vionpsqij '
            vionpsqij(il,it1,it2,it)=0.d0
            return
        end if

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(qij(1:nrc))
        allocate(wght(1:nrc))
        qij(1:nrc)=qrspspw(1:nrc,iq)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+qij(i)*vlocr2(i)*wght(i)
        end do
        vionpsqij(il,it1,it2,it)=sum

        deallocate(qij,wght)
        return
   end subroutine make_vionpsqij

   !==========================================
    subroutine cnstrct_of_CijkClmkVVVVijlm_k
   !==========================================
        integer:: lmt1,il1,it1,im1,ii
        integer:: lmt2,il2,it2,im2,jj
        integer:: lmt3,il3,it3,im3,kk
        integer:: lmt4,il4,it4,im4,ll
        integer:: lmt4min,nn,n,m,mp
        integer:: il12,il34,isph12,isph34
        logical:: matchflg
        logical,save:: initialized=.false.
        real(DP):: cijclm,cijclm_ae

        real(kind=DP), pointer, dimension(:,:,:) :: cr2
        integer, pointer, dimension(:,:,:)       :: isph2
        integer, pointer, dimension(:,:)         :: mmt2

        integer:: lt1,lt2,lt3,lt4,lt5,itmp
        integer,allocatable,dimension(:)        :: ilk

        if(.not.initialized) then
            m_clmns_cijkclmk=0
            initialized=.true.
        end if

        allocate(cr2(16,16,6))
        allocate(isph2(16,16,6))
        allocate(mmt2(16,16))

        call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)

        if(.not.paramset) then
            call m_PP_find_maximum_l(n)    ! n-1: maximum l
            n = (n-1) + (n-1) + 1
            allocate(ilk(n**2)); call substitute_il3(n**2,ilk) ! -(b_Elec..)
        end if

        do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it)
            im1=mtp(lmt1,it)
            it1=taup(lmt1,it)
            ii=(il1-1)**2+im1
            do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2

                nn=0

!                do lmt3=lmt1,ilmt(it)
!                    il3=ltp(lmt3,it)
!                    im3=mtp(lmt3,it)
!                    it3=taup(lmt3,it)
!                    kk=(il3-1)**2+im3
!
!                    lmt4min=lmt3
!                    if(lmt1==lmt3) lmt4min=max(lmt2,lmt3)
!                    do lmt4=lmt4min,ilmt(it)
!                        il4=ltp(lmt4,it)
!                        im4=mtp(lmt4,it)
!                        it4=taup(lmt4,it)
!                        ll=(il4-1)**2+im4

                do lmt3=1,ilmt(it)
                    il3=ltp(lmt3,it)
                    im3=mtp(lmt3,it)
                    it3=taup(lmt3,it)
                    kk=(il3-1)**2+im3
                    do lmt4=lmt3,ilmt(it)
                        il4=ltp(lmt4,it)
                        im4=mtp(lmt4,it)
                        it4=taup(lmt4,it)
                        ll=(il4-1)**2+im4

                        matchflg=.false.
                        cijclm=0.d0
                        cijclm_ae=0.d0
                        do n=1,mmt2(ii,jj)
                            isph12=isph2(ii,jj,n)
                            do m=1,mmt2(kk,ll)
                                isph34=isph2(kk,ll,m)
                                if(isph12==isph34) then
                                    matchflg=.true.
!print *,isph12,ilk(isph12)
                                    if(.not.paramset) then
                                        lt1=index_lmt2lt(lmt1,it)
                                        lt2=index_lmt2lt(lmt2,it)
                                        lt3=index_lmt2lt(lmt3,it)
                                        lt4=index_lmt2lt(lmt4,it)
                                        lt5=ilk(isph12)+1
!print *,'Before ',lt1,lt2,lt3,lt4,lt5
!                                        if(lt1>lt2) then
!                                            itmp=lt1
!                                            lt1=lt2
!                                            lt2=itmp
!                                        end if
!                                        if(lt3>lt4) then
!                                            itmp=lt3
!                                            lt3=lt4
!                                            lt4=itmp
!                                        end if
                                        if(lt1.gt.lt3 .or. (lt1.eq.lt3.and.lt2.gt.lt4)) then
                                            itmp=lt1
                                            lt1=lt3
                                            lt3=itmp
                                            itmp=lt2
                                            lt2=lt4
                                            lt4=itmp
                                        end if
                                        mp=ipppp(lt1,lt2,lt3,lt4,lt5,it)
!print *,'mp=',mp
                                        cijclm=cijclm+ &
                                                (vaeijlm_k(mp) &
                                                -vpsijlm_k(mp) &
                                                -vqijqlm_k(mp) &
                                                -vqijplpm_ks(mp))* &
                                                cr2(ii,jj,n)*cr2(kk,ll,m)
                                        cijclm_ae=cijclm_ae+ &
                                                vaeijlm_k(mp)* &
                                                cr2(ii,jj,n)*cr2(kk,ll,m)
!print '(4e19.6)',vaeijlm_k(mp),vpsijlm_k(mp),&
!vqijqlm_k(mp),vqijplpm_ks(mp)

                                    end if

                                end if
                            end do
                        end do

                        if(matchflg) then
                            nn=nn+1
                            if(.not.paramset) then
                                ilmt3_cijkclmk(lmt1,lmt2,nn,it)=lmt3
                                ilmt4_cijkclmk(lmt1,lmt2,nn,it)=lmt4
                                CijkClmkVVVVijlm_k(lmt1,lmt2,nn,it)=cijclm
                                CijkClmkVVVVijlm_k_ae(lmt1,lmt2,nn,it)=cijclm_ae
!print *,nn,'/',n_cijkclmk(lmt1,lmt2,it)
!print *,lmt1,lmt2,ilmt3_cijkclmk(lmt1,lmt2,nn,it),ilmt4_cijkclmk(lmt1,lmt2,nn,it)
!print '(4i5,e19.6)',lmt1,lmt2,lmt3,lmt4,cijclm
                            end if
!                            print *,nn
!                            print *,ii,jj,kk,ll
!                            print *,lmt1,lmt2,lmt3,lmt4
                        end if

                    end do
                end do
                if(.not.paramset) &
                    n_cijkclmk(lmt1,lmt2,it)=nn
                if(nn.gt.m_clmns_cijkclmk) m_clmns_cijkclmk=nn
            end do
        end do

        deallocate(cr2,isph2,mmt2)
        if(.not.paramset) deallocate(ilk)

!print *,'m_clmns_cijkclmk=',m_clmns_cijkclmk
!if(.not.paramset) stop
        return

   end subroutine cnstrct_of_CijkClmkVVVVijlm_k

   !===============================================
    subroutine cnstrct_of_CijkClmnVVVVijlm_kn
   !===============================================
        integer:: lmt1,il1,it1,im1,ii
        integer:: lmt2,il2,it2,im2,jj
        integer:: lmt3,il3,it3,im3,kk
        integer:: lmt4,il4,it4,im4,ll
        integer:: lmt4min,nn,n,m,mp
        integer:: il12,il34,isph12,isph34
        logical:: matchflg
        logical,save:: initialized=.false.
        real(DP):: cijclm,cijclm_ae

        real(kind=DP), pointer, dimension(:,:,:) :: cr2
        integer, pointer, dimension(:,:,:)       :: isph2
        integer, pointer, dimension(:,:)         :: mmt2

        integer:: lt1,lt2,lt3,lt4,lt5,itmp
        integer,allocatable,dimension(:)        :: ilk

        real(DP),allocatable :: crotylm(:,:,:)
        integer,allocatable :: iylm(:,:,:)
        integer,allocatable :: nylm(:,:)
        real(DP),allocatable :: opr(:,:,:)

        integer:: l1max,mmax,nsph,iopr
        integer:: ia,nm,l

        if(.not.initialized) then
            m_clmns_cijkclmk=0
            m_clmns_cijkclmn = 0
            initialized=.true.
        end if

        allocate(cr2(16,16,6))
        allocate(isph2(16,16,6))
        allocate(mmt2(16,16))

        call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)

!        call m_PP_find_maximum_l(n)    ! n-1: maximum l

        n = 0
        do lmt1 = 1, ilmt(it)
            l = ltp(lmt1,it)
            if(n < l) n = l
        end do

        n = (n-1) + (n-1) + 1
        l1max=n
        mmax=2*l1max-1
        nsph=l1max**2
!        if(.not.paramset) then
        allocate(ilk(n**2)); call substitute_il3(n**2,ilk) ! -(b_Elec..)
!        end if
        allocate(crotylm(mmax,nsph,nopr))
        allocate(iylm(mmax,nsph,nopr))
        allocate(nylm(nsph,nopr))
        allocate(opr(3,3,nopr))

        do ia=1,natm
            if(ityp(ia)/=it) cycle
            do iopr=1,nopr
                if(ia2ia_symmtry_op(ia,iopr).gt.0) then
                    opr(:,:,iopr)=op(:,:,iopr)
                else
                    opr(:,:,iopr)=-op(:,:,iopr)
                end if
            end do
            call get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,opr)

            do iopr=1,nopr

!print '(a,i4)','iopr=',iopr
!print '(3e19.6)',opr(1,:,iopr)
!print '(3e19.6)',opr(2,:,iopr)
!print '(3e19.6)',opr(3,:,iopr)

                do lmt1=1,ilmt(it)
                    il1=ltp(lmt1,it)
                    im1=mtp(lmt1,it)
                    it1=taup(lmt1,it)
                    ii=(il1-1)**2+im1
                    do lmt2=lmt1,ilmt(it)
                        il2=ltp(lmt2,it)
                        im2=mtp(lmt2,it)
                        it2=taup(lmt2,it)
                        jj=(il2-1)**2+im2

                        nn=0

        !                do lmt3=lmt1,ilmt(it)
        !                    il3=ltp(lmt3,it)
        !                    im3=mtp(lmt3,it)
        !                    it3=taup(lmt3,it)
        !                    kk=(il3-1)**2+im3
        !
        !                    lmt4min=lmt3
        !                    if(lmt1==lmt3) lmt4min=max(lmt2,lmt3)
        !                    do lmt4=lmt4min,ilmt(it)
        !                        il4=ltp(lmt4,it)
        !                        im4=mtp(lmt4,it)
        !                        it4=taup(lmt4,it)
        !                        ll=(il4-1)**2+im4

                        do lmt3=1,ilmt(it)
                            il3=ltp(lmt3,it)
                            im3=mtp(lmt3,it)
                            it3=taup(lmt3,it)
                            kk=(il3-1)**2+im3
                            do lmt4=lmt3,ilmt(it)
                                il4=ltp(lmt4,it)
                                im4=mtp(lmt4,it)
                                it4=taup(lmt4,it)
                                ll=(il4-1)**2+im4

                                matchflg=.false.
                                cijclm=0.d0
                                cijclm_ae=0.d0

                                do n=1,mmt2(ii,jj)
                                    isph12=isph2(ii,jj,n)
                                    do m=1,mmt2(kk,ll)
                                        isph34=isph2(kk,ll,m)

                                        if(ilk(isph12).ne.ilk(isph34)) cycle

                                        do nm=1,nylm(isph34,iopr)
                                            if(isph12.eq.iylm(nm,isph34,iopr)) then
                                                matchflg=.true.
!                                           if(isph12==isph34) then
!                                                matchflg=.true.
!print '(2i4,e19.6)',isph12,isph34,crotylm(nm,isph34,iopr)
                                                if(.not.paramset) then
                                                    lt1=index_lmt2lt(lmt1,it)
                                                    lt2=index_lmt2lt(lmt2,it)
                                                    lt3=index_lmt2lt(lmt3,it)
                                                    lt4=index_lmt2lt(lmt4,it)
                                                    lt5=ilk(isph12)+1
            !print *,'Before ',lt1,lt2,lt3,lt4,lt5
            !                                        if(lt1>lt2) then
            !                                            itmp=lt1
            !                                            lt1=lt2
            !                                            lt2=itmp
            !                                        end if
            !                                        if(lt3>lt4) then
            !                                            itmp=lt3
            !                                            lt3=lt4
            !                                            lt4=itmp
            !                                        end if
                                                    if(lt1.gt.lt3 .or. (lt1.eq.lt3.and.lt2.gt.lt4)) then
                                                        itmp=lt1
                                                        lt1=lt3
                                                        lt3=itmp
                                                        itmp=lt2
                                                        lt2=lt4
                                                        lt4=itmp
                                                    end if
                                                    mp=ipppp(lt1,lt2,lt3,lt4,lt5,it)
            !print *,'mp=',mp
                                                    cijclm=cijclm+ &
                                                            (vaeijlm_k(mp) &
                                                            -vpsijlm_k(mp) &
                                                            -vqijqlm_k(mp) &
                                                            -vqijplpm_ks(mp))* &
                                                            cr2(ii,jj,n)*cr2(kk,ll,m)* &
                                                            crotylm(nm,isph34,iopr)
                                                    cijclm_ae=cijclm_ae+ &
                                                            vaeijlm_k(mp)* &
                                                            cr2(ii,jj,n)*cr2(kk,ll,m)* &
                                                            crotylm(nm,isph34,iopr)
            !print '(4e19.6)',vaeijlm_k(mp),vpsijlm_k(mp),&
            !vqijqlm_k(mp),vqijplpm_ks(mp)

                                                end if
                                            end if
                                        end do

                                    end do
                                end do

                                if(matchflg) then
                                    nn=nn+1
                                    if(.not.paramset) then
                                        ilmt3_cijkclmn(lmt1,lmt2,nn,ia,iopr)=lmt3
                                        ilmt4_cijkclmn(lmt1,lmt2,nn,ia,iopr)=lmt4
                                        CijkClmnVVVVijlm_kn(lmt1,lmt2,nn,ia,iopr)=cijclm
                                        CijkClmnVVVVijlm_kn_ae(lmt1,lmt2,nn,ia,iopr)=cijclm_ae
        !print *,nn,'/',n_cijkclmk(lmt1,lmt2,it)
        !print *,lmt1,lmt2,ilmt3_cijkclmk(lmt1,lmt2,nn,it),ilmt4_cijkclmk(lmt1,lmt2,nn,it)
        !print '(4i5,e19.6)',lmt1,lmt2,lmt3,lmt4,cijclm
                                    end if
        !                            print *,nn
        !                            print *,ii,jj,kk,ll
        !                            print *,lmt1,lmt2,lmt3,lmt4
                                end if

                            end do
                        end do
                        if(.not.paramset) &
                            n_cijkclmn(lmt1,lmt2,ia,iopr)=nn
!print *,'lmt1,lmt2,ia,iopr,nn',lmt1,lmt2,ia,iopr,nn
                        if(nn.gt.m_clmns_cijkclmn) m_clmns_cijkclmn=nn
                    end do
                end do

            end do

        end do

        if(ipripp >= 2 .and. .not.paramset) then
           write(nfout,'(" - vaeijlm_k, vpsijlm_k, vqijqlm_k, vqijplpm_ks - <<cnstrct_of_CijkClmnVVVVijlm_kn")')
           do nn = 1, npppp
              write(nfout,'(4d20.8)') vaeijlm_k(nn),vpsijlm_k(nn), vqijqlm_k(nn), vqijplpm_ks(nn)
           end do
           write(nfout,'(" - CijkClmnVVVVijlm_kn -- <<cnstrct_of_CijkClmnVVVVijlm_kn")')
           write(nfout,'(" m_clmns_cijlcln = ", i8)') m_clmns_cijkclmn
           do ia = 1, natm
              do nn = 1, m_clmns_cijkclmn
                 do lmt1 = 1, ilmt(it)
                    do lmt2 = lmt1, ilmt(it)
                       write(nfout,'(" ia, nn, lmt1, lmt2 = ", 4i8)') ia, nn, lmt1, lmt2
                       write(nfout,'(12f13.8)') (cijkClmnVVVVijlm_kn(lmt1,lmt2,nn,ia,iopr),iopr=1,nopr)
                    end do
                 end do
              end do
           end do
        end if
!if(.not.paramset) stop
!print *,m_clmns_cijkclmn
!stop
        deallocate(cr2,isph2,mmt2)
!        if(.not.paramset) deallocate(ilk)
        deallocate(ilk)
        deallocate(crotylm,iylm,nylm,opr)
!        deallocate(ia2ia_symmtry_op)
!print *,'m_clmns_cijkclmk=',m_clmns_cijkclmk
!if(.not.paramset) stop
        return

   end subroutine cnstrct_of_CijkClmnVVVVijlm_kn

   !======================================
    subroutine make_ia2ia_symmtry_op_etc
   !======================================
        if(.not.associated(ia2ia_symmtry_op)) then
            allocate(ia2ia_symmtry_op(natm,nopr+af))
            ia2ia_symmtry_op=0
            call set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia_symmtry_op)
        end if
        if(.not.associated(ia2ia_symmtry_op_inv)) then
            allocate(ia2ia_symmtry_op_inv(natm,nopr+af))
            ia2ia_symmtry_op_inv=0
            call set_ia2ia_symmtry_op_inv &
                        (natm,nopr,ia2ia_symmtry_op,ia2ia_symmtry_op_inv)
        end if
        return
   end subroutine make_ia2ia_symmtry_op_etc

   !==========================================
    subroutine init_of_dion_paw
   !==========================================
      integer:: ia,is,lmt1,lmt2
      integer :: ismax

      if ( noncol ) then
         ismax = 1
      else
         ismax = nspin
      endif

      do ia=1,natm
         if(ityp(ia)/=it) cycle        ! ASMS 2019/12/03
! ================================ modified by K. Tagami ========== 11.0
!         do is=1,nspin
         do is=1, ismax
! ================================================================= 11.0
            do lmt1=1,ilmt(it)
               do lmt2=1,ilmt(it)
                  dion_paw(lmt1,lmt2,is,ia)=dion(lmt1,lmt2,it)
               end do
            end do
         end do
      end do
      return
    end subroutine init_of_dion_paw

!....................................................................
  end subroutine m_PP_vanderbilt_type
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q_END ORG_Parallel
! ==================================================================================================

  subroutine m_PP_qitgft_qmk()
    integer :: it
    if (sw_rspace_hyb == ON) return
    call alloc_qitg_exx()
    do it=1,ntyp
       if(m_PP_include_vanderbilt_pot(it)/=SKIP) then
!          call qitgft_qmk(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h)
          call qitgft_qmk( it, radial_aug_chg(it)%nmm_il3,  &
               &               radial_aug_chg(it)%mm_il3, &
               &               radial_aug_chg(it)%qrsps_mm, lcmax, h )
                                      ! in b_PseudoPotential_EXX
       endif
    enddo
    deallocate(qrsps_mm)
    deallocate(mm_il3,nmm_il3)
    deallocate(radr)
    deallocate(wos)
    deallocate(radial_aug_chg)

  end subroutine m_PP_qitgft_qmk

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q ORG_Parallel
! ==================================================================================================
  subroutine rd_qrsps_then_iqitg_and_qitgft(is_gncpp,nfp,it,nfout,gr_l &
       & ,ngshell,ngshell_range,paramset,mmt)
    integer, intent(in) :: is_gncpp
    integer, intent(in)       :: nfp, it, nfout,ngshell
    integer, intent(in), dimension(ngshell,2) :: ngshell_range
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical, intent(in)       :: paramset
    integer, intent(inout) :: mmt

    integer :: il1, il2, tau1, tau2, il3, tmin, l3s, l3l,  mm, i, igs
    integer :: nrc, mord, iprippex, iw, id
    real(kind=DP) :: dummy
    character(len=32):: cdummy
!    integer, allocatable, dimension(:)       :: nmm_il3 !d(lcmax+1)
!    integer, allocatable, dimension(:,:)     :: mm_il3 !d(nqitg_sp,lcmax+1)
!    real(kind=DP), allocatable, dimension(:,:) :: qrsps_mm !d(nmesh(it),nqitg_sp(it))

! ================= Added by K. Tagami ======== 10.1
    integer   :: nq
! ============================================= 10.1

    mm = 0
    iprippex = ipripp
    if(.not.paramset) then
       qij = 0.d0; qvij = 0.d0
       if(allocated(nmm_il3)) deallocate(nmm_il3)
       allocate(nmm_il3(lcmax+1)); nmm_il3 = 0
       if(allocated(mm_il3)) deallocate(mm_il3)
       allocate(mm_il3(nqitg_sp(it),lcmax+1)); mm_il3 = 0
       if(allocated(qrsps_mm)) deallocate(qrsps_mm)
       allocate(qrsps_mm(mmesh,nqitg_sp(it)))
       iprippex = iprippex-1

! === KT_add === 2014/09/19
       if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
          if ( use_asymm_ekin_density ) then
             allocate(kina_qrsps_mm(mmesh,nqitg_sp(it)))
          endif
          if ( use_symm_ekin_density ) then
             allocate(kins_qrsps_mm(mmesh,nqitg_sp(it)))
          endif
       endif
! ============== 2014/09/19
    end if
    Loop_L1 : do il1 = 1, lpsmax(it)
       if(iloc(it) == il1) cycle
       Loop_tau1 : do tau1 = 1, itau(il1,it)
          Loop_L2 : do il2 = il1, lpsmax(it)
             if((iloc(it)==il2) .or. (ivanl(il1,it)==0.and.ivanl(il2,it)==0)) cycle
             tmin = 1
             if(il1 == il2) tmin = tau1
             Loop_tau2 : do tau2 = tmin, itau(il2,it)
                l3s = abs(il1-il2)
                l3l = abs(il1+il2) - 2
                Loop_L3 : do il3 = l3s+1, l3l+1, 2
                   if(il3-1 > lcmax ) cycle
                   if(ipripp >= 2) write(nfout,'(" !PP il1, tau1, il2, tau2, il3 = ",5i5)') &
                        & il1, tau1, il2, tau2, il3
                   mm = mm + 1
                   if(.not.paramset) iqitg(il1,tau1,il2,tau2,il3,it) = mm + mmt
                   if(mype == 0) call read_nrc_mord(nfp,nfout, kord,&
                        &            mm,il1,tau1,il2,tau2,il3,nrc,mord,iprippex) ! ->(nrc,mord)
                   if(paramset .and. mype == 0) then
                      if(nrc > 0 ) then
                         read(nfp,*) (dummy,i=0,mord)
                         read(nfp,*) (dummy,i=nrc+1,nmesh(it))
                      else if(nrc == 0) then
                         read(nfp,*) (dummy,i=1,nmesh(it))
                      end if
                   else if(.not.paramset) then
                      nmm_il3(il3) = nmm_il3(il3)+1
                      mm_il3(nmm_il3(il3),il3) = mm
                      if(mype == 0) then
                         if(nrc > 0 ) then
                            read(nfp,*) (copsc(i),i=0,mord)
                            read(nfp,*) (qrsps(i),i=nrc+1,nmesh(it))
#ifdef _QRSPS_WITH_L_
#else
                            if(il3 > 1) then
                               do i = nrc+1, nmesh(it)
                                  qrsps(i) = qrsps(i)/(radr(i)**(il3-1))
                               end do
                            end if
#endif
                            call cnvrtp(nmesh(it),1,nrc,il3+1,mord,copsc,radr,qrsps)
!!$20110126
!!$                            if(ipripp>=1) then
!!$#ifdef _QRSPS_WITH_L_
!!$                               write(nfout,'(" _QRSPS_WITH_L_ is defined")')
!!$#else
!!$                               write(nfout,'(" _QRSPS_WITH_L_ is not defined")')
!!$#endif
!!$                               write(nfout,'(" qrsps ")')
!!$                               write(nfout,'(" nrc, mord, nmesh = ",3i8)') nrc, mord, nmesh(it)
!!$                               do i = 1, nmesh(it)
!!$                                  write(nfout,'(i8,2e12.4)') i, radr(i),qrsps(i)
!!$                               end do
!!$                            end if
!!$20110126
                         else if(nrc == 0 ) then
                            read(nfp,*) (qrsps(i),i=1,nmesh(it))
#ifdef _QRSPS_WITH_L_
#else
                            if(il3 > 1) then
                               do i = 1, nmesh(it)
                                  qrsps(i) = qrsps(i)/(radr(i)**(il3-1))
                               end do
                            end if
#endif
                         endif
! -->
!!$                         if(is_gncpp /= pp_PAW) then
! === KT_mod === 2015/02/23
!!!                         if(sw_berry_phase == ON .or. sw_wannier == ON .or. sw_fef == ON) then
                         if(sw_berry_phase == ON .or. sw_wannier == ON &
                              &                  .or. sw_fef == ON &
                              &                  .or. sw_wannier90 == ON ) then
! ============== 2015/02/23
                            if(is_gncpp == pp_GNCPP2 .and. ae_wavefunctions_are_detected(it)) then
                               do i=1,nmesh(it)
                                  qrs(i) = (phir_ae(i,il1,tau1)*phir_ae(i,il2,tau2) - &
                                       & phir(i,il1,tau1)*phir(i,il2,tau2))*radr(i)**(il3-1)
                               end do
                            else
                               qrs = qrsps
                            end if
                         end if
!!$                         end if
! <--
! ================== Added by K. Tagami =================== 10.1
                         if ( sw_LinearResponse == ON ) then
                            if (is_gncpp == 2 .and. ae_wavefunctions_are_detected(it)) then
                               do i=1,nmesh(it)
                                  qrs(i) = (phir_ae(i,il1,tau1)*phir_ae(i,il2,tau2) - &
                                       & phir(i,il1,tau1)*phir(i,il2,tau2))*radr(i)**(il3-1)
                               end do
                            else
                               qrs = qrsps
                            end if
                         endif
! ========================================================= 10.1
                      end if

                      if(npes > 1) &
                           & call mpi_bcast(qrsps,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)

                      qrsps_mm(:,mm) = qrsps(:)

!                      if(is_gncpp == pp_PAW) then
!                         if(.not.paramset .and. ipaw(it) == ON) qrspspw(1:mmesh,mm+mmt)=qrsps(1:mmesh)
!                      end if
                       if(.not.paramset) qrspspw(1:mmesh,mm+mmt)=qrsps(1:mmesh)

! -->
!!$                     if(is_gncpp /= pp_PAW) then
                      if((sw_berry_phase == ON .or. sw_wannier==ON .or. sw_fef==ON) &
                           & .and.npes > 1) &
                           & call mpi_bcast(qrs,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)

                      if ( sw_wannier90 == ON .and.npes > 1 ) &
                           & call mpi_bcast(qrs,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)

! ======================== Added by K. Tagami =========== 10.1
                      if ( sw_LinearResponse == ON .and.npes > 1 ) &
                           & call mpi_bcast(qrs,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
! ======================================================= 10.1

! ==== KT_add == 2014/09/19
                      if ( sw_calc_ekin_density == ON &
                           &          .and. sw_rspace_ekin_density == OFF ) then
                         if ( use_asymm_ekin_density ) then
                            call calc_kin_qrs_asymm( nrc, il1, il2, tau1, tau2, &
                                 &                   kina_qrs )
!                            kina_qrs(1:mmesh) = kina_qrs(1:mmesh) &
!                                 &           *radr(1:mmesh)**(il3-1 )
                           if (npes > 1) then
                               call mpi_bcast( kina_qrs, nmesh(it), &
                                    &          mpi_double_precision, 0, &
                                    &          MPI_CommGroup, ierr )
                            endif
                            kina_qrsps_mm(:,mm) = kina_qrs(:)
                         endif
                         if ( use_symm_ekin_density ) then
                            call calc_kin_qrs_symm( nrc, il1, il2, tau1, tau2, kins_qrs )
!                            kins_qrs(1:mmesh) = kins_qrs(1:mmesh) &
!                                 &           *radr(1:mmesh)**(il3-1 )
                            if (npes > 1) then
                               call mpi_bcast( kins_qrs, nmesh(it), &
                                    &          mpi_double_precision, 0, &
                                    &          MPI_CommGroup, ierr )
                            endif
                            kins_qrsps_mm(:,mm) = kins_qrs(:)
                         endif
                      endif
! ============== 2014/09/19

                      if(sw_berry_phase == ON) then
#if 0
                         if(is_gncpp == pp_GNCPP1) then
                            call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_BP &
                                     &,1.d0,mm+mmt,qitg_BP,wkx,wky)
                         else if(is_gncpp == pp_GNCPP2) then
                            call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_BP &
                                     &,1.d0,mm+mmt,qitg_BP,wkx,wky)
                         end if
#else
                         call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                              &,nfout,il3-1,radr,wos,qrs,dk_BP &
                              &,1.d0,mm+mmt,qitg_BP,wkx,wky)
#endif
                         if(ipripp >= 2) then
                            write(nfout,'(" !PP i=",i3," q_BP=",e25.12)')   mm+mmt,qitg_BP(1,mm+mmt)
                            write(nfout,'(" !PP ",5x," q_BP/pi4=",e25.12)') qitg_BP(1,mm+mmt)/PAI4
                         end if
                      end if

                      if ( sw_wannier == ON .or. sw_wannier90 == ON ) then
                         ! debug
                         if(ipripp >= 2) write(nfout,*) 'rd_qrsps_then_iqitg_and_qitgft: nwght_wan=',nwght_wan
                         ! end debug
                         do iw=1,nwght_wan
#if 0
                            if(is_gncpp == pp_GNCPP1) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_wan(iw) &
                                     &,1.d0,mm+mmt,qitg_wan(1,iw),wkx,wky)
                            else if(is_gncpp == pp_GNCPP2) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_wan(iw) &
                                     &,1.d0,mm+mmt,qitg_wan(1,iw),wkx,wky)
                            end if
#else
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                    &,nfout,il3-1,radr,wos,qrs,dk_wan(iw) &
                                    &,1.d0,mm+mmt,qitg_wan(1,iw),wkx,wky)
#endif
                            if(ipripp >= 2) then
                               write(nfout,'(" !PP i=",i3,"iw=",i3," q_wan=",e25.12)')   mm+mmt,iw,qitg_wan(mm+mmt,iw)
                               write(nfout,'(" !PP ",11x," q_wan/pi4=",e25.12)') qitg_wan(mm+mmt,iw)/PAI4
                            end if
                         end do
                      end if
!!$                     end if
! <--

! ================================= Added by K. Tagami ======== 10.1
                      if ( sw_LinearResponse == ON ) then
                         ! debug
!                         write(6,*) 'rd_qrsps_then_iqitg_and_qitgft: nwght_wan=',nwght_wan                         ! end debug

                         ! ---
                         Do nq=1, nmax_q_LR
#if 0
                            if ( is_gncpp == 1 ) then
                               call qitgft( 1, nmax_q_plus_G_LR, mmesh, &
                                    & nqitg, nmesh(it), &
                                    & ipri, nfout, il3-1, radr, wos, qrsps, &
                                    & norm_q_plus_G_LR(1,nq), &
                                    & 1.d0, mm+mmt, &
                                    & qitg_LR(1,1,nq), &
                                    & wkx,wky )
                            else if ( is_gncpp == 2 ) then
                               call qitgft( 1, nmax_q_plus_G_LR, mmesh, &
                                    & nqitg, nmesh(it), &
                                    & ipri, nfout, il3-1, radr, wos, qrs, &
                                    & norm_q_plus_G_LR(1,nq), &
                                    & 1.d0, mm+mmt, &
                                    & qitg_LR(1,1,nq), &
                                    & wkx,wky )
                            end if
#else
                            call qitgft( 1, nmax_q_plus_G_LR, mmesh, &
                                 & nqitg, nmesh(it), &
                                 & ipri, nfout, il3-1, radr, wos, qrs, &
                                 & norm_q_plus_G_LR(1,nq), &
                                 & 1.d0, mm+mmt, &
                                 & qitg_LR(1,1,nq), &
                                 & wkx,wky )
#endif
                         End do
                      end if
! ================================================================= 10.1
                      if(sw_fef == ON) then
                         ! debug
                           write(6,*) 'rd_qrsps_then_iqitg_and_qitgft: numef_fef=',numef_fef
                         ! end debug
                         do id=1,numef_fef
#if 0
                            if(is_gncpp == 1) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_fef(id) &
                                     &,1.d0,mm+mmt,qitg_fef(1,id),wkx,wky)
                            else if(is_gncpp == 2) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_fef(id) &
                                     &,1.d0,mm+mmt,qitg_fef(1,id),wkx,wky)
                            end if
#else
                            call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                 &,nfout,il3-1,radr,wos,qrs,dk_fef(id) &
                                 &,1.d0,mm+mmt,qitg_fef(1,id),wkx,wky)
#endif
                            if(ipripp >= 2) then
                               write(nfout,'(" !PP i=",i3,"id=",i3," q_fef=",e25.12)')   mm+mmt,id,qitg_fef(mm+mmt,id)
                               write(nfout,'(" !PP ",11x," q_fef/pi4=",e25.12)') qitg_fef(mm+mmt,id)/PAI4
                            end if
                         end do
                      end if

                      if(il3 == 1 .and. il1 == il2) &
                           & call qij_qvij_from_qrsps_etc(it,tau1,tau2,il1)
                   end if
                end do Loop_L3
             end do Loop_tau2
          end do Loop_L2
       end do Loop_tau1
    end do Loop_L1

    mmt = mmt + mm
    if(paramset) then
       nqitg = mmt
       nqitg_sp(it) = mm
    end if

    call wd_nqitg_sp_etc()
!!$    if(.not.paramset) write(6,'(" nqitg = ",i5)') nqitg !K.Mae print->write

    if(.not.paramset .and. iflag_ft == 1) call qitgft_mm() ! contained here

! ==== KT_add ==== 2014/09/19
    if(.not.paramset .and. iflag_ft == 1) then
       if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
          if ( use_asymm_ekin_density ) then
             call kin_qitgft_mm( nqitg_sp(it), kina_qrsps_mm, kina_qitg_l )
          endif
          if ( use_symm_ekin_density ) then
             call kin_qitgft_mm( nqitg_sp(it), kins_qrsps_mm, kins_qitg_l )
          endif
       endif
    endif
! ================ 2014/09/19

!    if(.not.paramset .and. sw_hybrid_functional==ON .and.  m_PP_include_vanderbilt_pot(it)/=SKIP.and.sw_rspace_hyb==OFF) &
!      & call qitgft_qmk(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h) ! in b_PseudoPotential_EXX

    if(.not.paramset .and. sw_hybrid_functional==ON .and. sw_rspace_hyb==OFF) then
       if ( m_PP_include_vanderbilt_pot(it)/=SKIP ) then
          allocate( radial_aug_chg(it)%nmm_il3(lcmax+1) )
          allocate( radial_aug_chg(it)%mm_il3(nqitg_sp(it),lcmax+1) )
          allocate( radial_aug_chg(it)%qrsps_mm(mmesh,nqitg_sp(it)) )
          radial_aug_chg(it)%nmm_il3 = nmm_il3
          radial_aug_chg(it)%mm_il3  = mm_il3
          radial_aug_chg(it)%qrsps_mm = qrsps_mm
       endif
    endif

    if(itpcc(it) == 1) then
       if(mype == 0) then
          call read_chgpc_nrc_mord(iprippex,nfp,nfout,chgpc(it),nrc,mord)
          ! --- > rhpcr,  copsc is a temporary uesd parameter
          if(nrc > 0) then
             if(paramset) then
                read(nfp,*) (dummy,i=0,mord)
                read(nfp,*) (dummy,i=nrc+1,nmesh(it))
             else
                read(nfp,*) (copsc(i),i=0,mord)
                read(nfp,*) (rhpcr(i),i=nrc+1,nmesh(it))
                call cnvrtp(mmesh,1,nrc,2,mord,copsc,radr,rhpcr)
             end if
          else if(nrc == 0) then
             if(paramset) then
                read(nfp,*) (dummy, i=1,nmesh(it))
             else
                read(nfp,*) (rhpcr(i),i=1,nmesh(it))
             end if
          endif
          ! <---- rhpcr
       end if
       if(npes > 1) call mpi_bcast(chgpc(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(.not.paramset) then
          if(npes > 1) call mpi_bcast(rhpcr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
          if(is_gncpp == pp_PAW .and. ipaw(it) == ON) rhpcrpw(:,it)=rhpcr(:)
       end if
    endif

    if(is_gncpp == pp_PAW) then
       if(mype == 0) then
          read(nfp,*) cdummy
          if(ipaw(it) == ON) then
             read(nfp,*) chgcr(it)
             write(nfout,400) chgcr(it)
          else
             read(nfp,*) dummy
          end if
400       format(' ',' chgcr = ',f12.8)
    ! --- > rhpcr,  copsc is a temporary uesd parameter
          if( paramset ) then
             read(nfp,*) (dummy, i=1,nmesh(it))
          else if ( ipaw(it) == OFF .and. sw_add_corecharge_rspace==OFF ) then
             read(nfp,*) (dummy, i=1,nmesh(it))
          else
             read(nfp,*) (rhcorpw(i,it),i=1,nmesh(it))
          end if
        ! <---- rhcr ( rhcorwpw )
       end if

       if(npes > 1 .and. ipaw(it) == ON) then
          call mpi_bcast(chgcr(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
       endif
       if ( .not. paramset ) then
          if ( npes > 1 ) then
             if ( ipaw(it)==ON .or. sw_add_corecharge_rspace ==ON ) then
                call mpi_bcast( rhcorpw(:,it), nmesh(it), mpi_double_precision, &
                     &          0, MPI_CommGroup, ierr )
             endif
          endif
       endif

       if(itpcc(it) == 0 .and. ipaw(it)==ON) then
          if(.not.paramset) rhpcrpw(:,it)=0.d0  !rhcorpw(:,it) !0.d0      !2009.4.18
       end if
    end if

    if(.not.paramset) then
       if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) then
         deallocate(qrsps_mm)
         deallocate(mm_il3,nmm_il3)
       endif
! === KT_add === 2014/09/19
       if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == OFF ) then
          if ( use_asymm_ekin_density  ) then
             deallocate(kina_qrsps_mm)
          endif
          if ( use_symm_ekin_density  ) then
             deallocate(kins_qrsps_mm)
          endif
       endif
! ============== 2014/09/19
    end if

  contains

    subroutine qitgft_mm()
      integer :: mm0, i, il3, mm, mmp, n
      real(kind=DP) :: qitg_sh, qitg_diff_sh, gabs, f2
!!$      integer :: mmp1,mmp2,mmp3,mmp4
!!$      real(kind=DP) :: qitg_sh1, qitg_sh2, qitg_sh3, qitg_sh4
!!$      integer, parameter :: NFOLDING = 4 ! {2|4}
!!$      integer :: mc

      integer, dimension(:), allocatable       :: ip_largerequal_1 ! d(ngshell)

      allocate(ip_largerequal_1(ngshell))
      if(ipripp >= 2) then
         write(nfout,'(" !qitgft_mm  ngshell = ",i8)') ngshell
         write(nfout,'(" !qitgft_mm --- gshell, ip_largerequal_1(=n), radr, G*radr(n-1), G*radr(n)")')
         write(nfout,'(" !PP ngshell_range")')
         write(nfout,'(" !PP ngshell, range1 - range2, nelement, gr")')
         do i = 1, ngshell
            write(nfout,'(" !PP ",i8, i8," - ",i8, i8, f20.8)') i, ngshell_range(i,1)&
                 & ,ngshell_range(i,2),ngshell_range(i,2)-ngshell_range(i,1)+1,gr_l(ngshell_range(i,1))
         end do
      end if

      do i = 1, ngshell
         gabs = gr_l(ngshell_range(i,1))
         if(gabs < DELTA) then
            ip_largerequal_1(i) = nmesh(it)+1
         else
            ip_largerequal_1(i) = ceiling(nmesh(it) - dlog(rmax(it)*gabs)/h(it))
         end if

         if(ipripp >= 2) then
            n = ip_largerequal_1(i)
            if(n <= 1 .or. n > nmesh(it)) then
               if(ipripp >= 3 .or.(ipripp == 2 .and. &
                    & (i <= min(30,ngshell) .or. (i >= max(ngshell-29,1))))) then
                  if(n == 1) then
                     write(nfout,'(" !qitgft_mm : ",i5, i8, "          ",f16.8)') i, n, radr(n)
                  else
                     write(nfout,'(" !qitgft_mm : ",i5, i8)') i, n
                  end if
               end if
            else
               write(nfout,'(" !qitgft_mm : ",i5, i8, 3f16.8)') i,n,radr(n),gabs*radr(n-1),gabs*radr(n)
            end if
         end if
      end do

      mm0 = 0
      do i = 1, it-1
         mm0 = mm0 + nqitg_sp(i)
      end do
      if(ipripp >= 2) write(nfout,'(" !PP mm0 = ", i5)') mm0


      do il3 = 1, lcmax+1
         if(nmm_il3(il3) <= 0) cycle
!!$         mc = nmm_il3(il3)/NFOLDING
         if(istress == 0) then
            if(ipripp >= 2) then
               write(nfout,'(" !PP il3 = ", i5, " nmm = ",i5)') il3, nmm_il3(il3)
               do mm = 1, nmm_il3(il3)
                  write(nfout,'(" !PP mm = ", i5)') mm_il3(mm,il3)
               end do

               write(nfout,'(" qrsps_mm <<qitgft_mm>>")')
               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  write(nfout,'(" il3, mm, mmp = ",3i8)') il3,mm,mmp
                  do i = 1, min(10, nmesh(it))
                     write(nfout,'(" !PP qrsps_mm: ",i6,d16.8)') i, qrsps_mm(i,mmp)
                  end do
                  do i = max(nmesh(it)-9,1),nmesh(it)
                     write(nfout,'(" !PP qrsps_mm: ",i6,d16.8)') i, qrsps_mm(i,mmp)
                  end do
               end do
            end if

            do i = 1, ngshell
               gabs = gr_l(ngshell_range(i,1))
               wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
               call dsjnvn(il3-1,nmesh(it),ip_largerequal_1(i),wkx,wky)

               if(ipripp >= 2) then
                  if(i <= min(5,ngshell)) then
                     write(nfout,'(" !PP ishell = ",i5)') i
                     do mm = 1, min(30, nmesh(it))
                        write(nfout,'(" !PP wkx,wky = ",i5,2d16.8)') mm, wkx(mm),wky(mm)
                     end do
                     do mm = max(nmesh(it)-29,1), nmesh(it)
                        write(nfout,'(" !PP wkx,wky = ",i5,2d16.8)') mm, wkx(mm),wky(mm)
                     end do
                  end if
               end if

               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  qitg_sh = 0.d0
                  do n = 1, nmesh(it)
                     qitg_sh = qitg_sh + wos(n)*qrsps_mm(n,mmp)*wky(n)
                  end do
                  qitg_sh = qitg_sh/univol*PAI4
                  do igs = ngshell_range(i,1), ngshell_range(i,2)
                     qitg_l(igs,mmp+mm0) = qitg_sh
                  end do
               end do

!!$               if(NFOLDING == 2) then
!!$                  do mm = 1, mc*NFOLDING, NFOLDING
!!$                     mmp1 = mm_il3(mm,il3)
!!$                     mmp2 = mm_il3(mm+1,il3)
!!$                     qitg_sh1 = 0.d0
!!$                     qitg_sh2 = 0.d0
!!$                     do n = 1, nmesh(it)
!!$                        f2 = wos(n)*wky(n)
!!$                        qitg_sh1 = qitg_sh1 + f2*qrsps_mm(n,mmp1)
!!$                        qitg_sh2 = qitg_sh2 + f2*qrsps_mm(n,mmp2)
!!$                     end do
!!$                     qitg_sh1 = qitg_sh1/univol*PAI4
!!$                     qitg_sh2 = qitg_sh2/univol*PAI4
!!$                     do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                        qitg_l(igs,mmp1+mm0) = qitg_sh1
!!$                        qitg_l(igs,mmp2+mm0) = qitg_sh2
!!$                     end do
!!$                  end do
!!$               else if(NFOLDING == 4) then
!!$                  do mm = 1, mc*NFOLDING, NFOLDING
!!$                     mmp1 = mm_il3(mm,il3)
!!$                     mmp2 = mm_il3(mm+1,il3)
!!$                     mmp3 = mm_il3(mm+2,il3)
!!$                     mmp4 = mm_il3(mm+3,il3)
!!$                     qitg_sh1 = 0.d0
!!$                     qitg_sh2 = 0.d0
!!$                     qitg_sh3 = 0.d0
!!$                     qitg_sh4 = 0.d0
!!$                     do n = 1, nmesh(it)
!!$                        f2 = wos(n)*wky(n)
!!$                        qitg_sh1 = qitg_sh1 + f2*qrsps_mm(n,mmp1)
!!$                        qitg_sh2 = qitg_sh2 + f2*qrsps_mm(n,mmp2)
!!$                        qitg_sh3 = qitg_sh3 + f2*qrsps_mm(n,mmp3)
!!$                        qitg_sh4 = qitg_sh4 + f2*qrsps_mm(n,mmp4)
!!$                     end do
!!$                     qitg_sh1 = qitg_sh1/univol*PAI4
!!$                     qitg_sh2 = qitg_sh2/univol*PAI4
!!$                     qitg_sh3 = qitg_sh3/univol*PAI4
!!$                     qitg_sh4 = qitg_sh4/univol*PAI4
!!$                     do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                        qitg_l(igs,mmp1+mm0) = qitg_sh1
!!$                        qitg_l(igs,mmp2+mm0) = qitg_sh2
!!$                        qitg_l(igs,mmp3+mm0) = qitg_sh3
!!$                        qitg_l(igs,mmp4+mm0) = qitg_sh4
!!$                     end do
!!$                  end do
!!$               end if
!!$               do mm = mc*NFOLDING+1,nmm_il3(il3)
!!$                  mmp = mm_il3(mm,il3)
!!$                  qitg_sh = 0.d0
!!$                  do n = 1, nmesh(it)
!!$                     qitg_sh = qitg_sh + wos(n)*qrsps_mm(n,mmp)*wky(n)
!!$                  end do
!!$                  qitg_sh = qitg_sh/univol*PAI4
!!$                  do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                     qitg_l(igs,mmp+mm0) = qitg_sh
!!$                  end do
!!$               end do
!!$
            end do
         else if(istress == 1) then
            do i = 1, ngshell
               gabs = gr_l(ngshell_range(i,1))
               wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
               call dsjnvn(il3-1,nmesh(it),ip_largerequal_1(i),wkx,wky)
               if(il3-1 == 0) then
                  wkz1 = 0.d0
               else
                  call dsjnvn(il3-2,nmesh(it),ip_largerequal_1(i),wkx,wkz1)
               end if
               call dsjnvn(il3,nmesh(it),ip_largerequal_1(i),wkx,wkz2)
               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  qitg_sh = 0.d0
                  qitg_diff_sh = 0.d0
                  do n = 1, nmesh(it)
                     f2 = wos(n)*qrsps_mm(n,mmp)
                     qitg_sh = qitg_sh + f2*wky(n)
                     qitg_diff_sh = qitg_diff_sh + f2 * radr(n) &
                          & / (2.d0*(il3-1)+1.d0) * (il3*wkz1(n)-(il3)*wkz2(n))
                  enddo
                  qitg_sh = qitg_sh/univol*PAI4
                  qitg_diff_sh = qitg_diff_sh/univol*PAI4
                  do igs = ngshell_range(i,1), ngshell_range(i,2)
                     qitg_l(igs,mmp+mm0) = qitg_sh
                     qitg_diff_l(igs,mmp+mm0) = qitg_diff_sh
                  end do
               end do
            end do
         end if
      end do

      if(ipripp >= 2) then
         write(nfout,'(" qitg_l <<qitgft_mm>>")')
         do mm = 1, nqitg
            write(nfout,'(" mm = ", i8)') mm
            write(nfout,'(6d12.4)') &
                 & (qitg_l(igs,mm),igs=ista_kngp, min(ista_kngp+60,iend_kngp))
         end do
      end if
      deallocate(ip_largerequal_1)

    end subroutine qitgft_mm

! === KT_add ==== 2014/09/19
    subroutine calc_kin_qrs_asymm( nrc, il1, il2, tau1, tau2, kin_qrs )
      integer, intent(in) :: nrc, il1, il2, tau1, tau2
      real(kind=DP), intent(out) :: kin_qrs( mmesh )

      integer :: ier, ir
      real(kind=DP) :: c1, c2, c3, c4, coeff1, coeff2, ctmp1, ctmp2, r
      real(kind=DP), allocatable ::   psi_i(:), dpsi_i(:), ddpsi_i(:)
      real(kind=DP), allocatable ::   psi_j(:), dpsi_j(:), ddpsi_j(:)

      real(kind=DP), allocatable ::   phi_i(:), dphi_i(:), ddphi_i(:)
      real(kind=DP), allocatable ::   phi_j(:), dphi_j(:), ddphi_j(:)

      allocate( psi_i( nrc )); allocate( dpsi_i( nrc ));   allocate( ddpsi_i( nrc ));
      allocate( psi_j( nrc )); allocate( dpsi_j( nrc ));   allocate( ddpsi_j( nrc ));
      allocate( phi_i( nrc )); allocate( dphi_i( nrc ));   allocate( ddphi_i( nrc ));
      allocate( phi_j( nrc )); allocate( dphi_j( nrc ));   allocate( ddphi_j( nrc ));

      psi_i = 0.0d0;   dpsi_i = 0.0d0;   ddpsi_i = 0.0d0
      psi_j = 0.0d0;   dpsi_j = 0.0d0;   ddpsi_j = 0.0d0
      phi_i = 0.0d0;   dphi_i = 0.0d0;   ddphi_i = 0.0d0
      phi_j = 0.0d0;   dphi_j = 0.0d0;   ddphi_j = 0.0d0

      kin_qrs = 0.0d0

      if (ipaw(it) == ON) then
         psi_i(1:nrc) = psir(1:nrc,il1,tau1)
         psi_j(1:nrc) = psir(1:nrc,il2,tau2)
      else if ( ae_wavefunctions_are_detected(it) ) then
         psi_i(1:nrc) = phir_ae(1:nrc,il1,tau1)
         psi_j(1:nrc) = phir_ae(1:nrc,il2,tau2)
      endif
      phi_i(1:nrc) = phir(1:nrc,il1,tau1)
      phi_j(1:nrc) = phir(1:nrc,il2,tau2)

      call calc_ddiff_exp( ier,5, nrc, radr, psi_i, dpsi_i, ddpsi_i )
      call calc_ddiff_exp( ier,5, nrc, radr, psi_j, dpsi_j, ddpsi_j )

      call calc_ddiff_exp( ier,5, nrc, radr, phi_i, dphi_i, ddphi_i )
      call calc_ddiff_exp( ier,5, nrc, radr, phi_j, dphi_j, ddphi_j )

      coeff1 = dble( il1 )*dble( il1 +1 )
      coeff2 = dble( il2 )*dble( il2 +1 )

      do ir=1,nrc
         r = radr(ir)
         c1 = -psi_j(ir) *ddpsi_i(ir) +coeff1 *psi_i(ir) *psi_j(ir) /r**2
         c2 = -phi_j(ir) *ddphi_i(ir) +coeff1 *phi_i(ir) *phi_j(ir) /r**2

         c3 = -psi_i(ir) *ddpsi_j(ir) +coeff2 *psi_i(ir) *psi_j(ir) /r**2
         c4 = -phi_i(ir) *ddphi_j(ir) +coeff2 *phi_i(ir) *phi_j(ir) /r**2

         ctmp1 = ( c1 -c2 ) /2.0d0
         ctmp2 = ( c3 -c4 ) /2.0d0
         kin_qrs(ir) = ( ctmp1 + ctmp2 )*0.5d0     ! averaging 'cause hsr is symmetric.
      end do

      deallocate( psi_i ); deallocate( dpsi_i ); deallocate( ddpsi_i )
      deallocate( psi_j ); deallocate( dpsi_j ); deallocate( ddpsi_j )

      deallocate( phi_i ); deallocate( dphi_i ); deallocate( ddphi_i )
      deallocate( phi_j ); deallocate( dphi_j ); deallocate( ddphi_j )

    end subroutine calc_kin_qrs_asymm

    subroutine calc_kin_qrs_symm( nrc, il1, il2, tau1, tau2, kin_qrs )
      integer, intent(in) :: nrc, il1, il2, tau1, tau2
      real(kind=DP), intent(out) :: kin_qrs( mmesh )

      integer :: ier, ir
      real(kind=DP) :: c1, c2, coeff, r
      real(kind=DP), allocatable ::   psi_i(:)
      real(kind=DP), allocatable ::  dpsi_i(:)
      real(kind=DP), allocatable ::   psi_j(:)
      real(kind=DP), allocatable ::  dpsi_j(:)

      real(kind=DP), allocatable ::   phi_i(:)
      real(kind=DP), allocatable ::  dphi_i(:)
      real(kind=DP), allocatable ::   phi_j(:)
      real(kind=DP), allocatable ::  dphi_j(:)

      allocate(   psi_i( nrc ));   psi_i = 0.0d0
      allocate(  dpsi_i( nrc ));  dpsi_i = 0.0d0
      allocate(   psi_j( nrc ));   psi_j = 0.0d0
      allocate(  dpsi_j( nrc ));  dpsi_j = 0.0d0

      allocate(   phi_i( nrc ));   phi_i = 0.0d0
      allocate(  dphi_i( nrc ));  dphi_i = 0.0d0
      allocate(   phi_j( nrc ));   phi_j = 0.0d0
      allocate(  dphi_j( nrc ));  dphi_j = 0.0d0

      kin_qrs = 0.0d0

      if (ipaw(it) == ON) then
         psi_i(1:nrc) = psir(1:nrc,il1,tau1)
         psi_j(1:nrc) = psir(1:nrc,il2,tau2)
      else if ( ae_wavefunctions_are_detected(it) ) then
         psi_i(1:nrc) = phir_ae(1:nrc,il1,tau1)
         psi_j(1:nrc) = phir_ae(1:nrc,il2,tau2)
      endif
      phi_i(1:nrc) = phir(1:nrc,il1,tau1)
      phi_j(1:nrc) = phir(1:nrc,il2,tau2)

      call calc_diff_exp( ier,5, nrc, radr, psi_i, dpsi_i )
      call calc_diff_exp( ier,5, nrc, radr, phi_i, dphi_i )
      call calc_diff_exp( ier,5, nrc, radr, psi_j, dpsi_j )
      call calc_diff_exp( ier,5, nrc, radr, phi_j, dphi_j )

      coeff = dble( il2 )*dble( il2 +1 )

      do ir=1,nrc
         r = radr(ir)
         c1 = ( dpsi_i(ir) -psi_i(ir) /r ) *( dpsi_j(ir) -psi_j(ir) /r )
         c1 = c1 + coeff *psi_i(ir) *psi_j(ir) /r**2

         c2 = ( dphi_i(ir) -phi_i(ir) /r ) *( dphi_j(ir) -phi_j(ir) /r )
         c2 = c2 + coeff *phi_i(ir) *phi_j(ir) /r**2
         kin_qrs(ir) = ( c1 -c2 )*0.5d0
      end do

      deallocate( psi_i ); deallocate( dpsi_i )
      deallocate( psi_j ); deallocate( dpsi_j )
      deallocate( phi_i ); deallocate( dphi_i )
      deallocate( phi_j ); deallocate( dphi_j )

    end subroutine calc_kin_qrs_symm

    subroutine kin_qitgft_mm( nqitg_in, kin_qrsps_mm, kin_qitg_l )
      integer, intent(in) :: nqitg_in
      real(kind=DP), intent(in) :: kin_qrsps_mm( mmesh, nqitg_in )
      real(kind=DP), intent(out) :: kin_qitg_l( ista_kngp:iend_kngp, nqitg )

      integer :: mm0, i, il3, mm, mmp, n
      real(kind=DP) :: qitg_sh, qitg_diff_sh, gabs, f2

      integer, dimension(:), allocatable       :: ip_largerequal_1 ! d(ngshell)
!
      allocate(ip_largerequal_1(ngshell))
      do i = 1, ngshell
         gabs = gr_l(ngshell_range(i,1))
         if(gabs < DELTA) then
            ip_largerequal_1(i) = nmesh(it)+1
         else
            ip_largerequal_1(i) = ceiling(nmesh(it) - dlog(rmax(it)*gabs)/h(it))
         end if

      end do

      mm0 = 0
      do i = 1, it-1
         mm0 = mm0 + nqitg_sp(i)
      end do
      if(ipripp >= 2) write(nfout,'(" !PP mm0 = ", i5)') mm0

      do il3 = 1, lcmax+1
         if(nmm_il3(il3) <= 0) cycle

!!         if(istress == 0) then

            do i = 1, ngshell
               gabs = gr_l(ngshell_range(i,1))
               wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
               call dsjnvn(il3-1,nmesh(it),ip_largerequal_1(i),wkx,wky)

               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  qitg_sh = 0.d0
                  do n = 1, nmesh(it)
                     qitg_sh = qitg_sh + wos(n)*kin_qrsps_mm(n,mmp)*wky(n)
                  end do
                  qitg_sh = qitg_sh/univol*PAI4
                  do igs = ngshell_range(i,1), ngshell_range(i,2)
                     kin_qitg_l(igs,mmp+mm0) = qitg_sh
                  end do
               end do
            end do

!!         end if
      end do

      deallocate(ip_largerequal_1)

    end subroutine kin_qitgft_mm
! ==================== 2014/09/19

    subroutine wd_nqitg_sp_etc()
      integer :: il3
      if(ipripp >= 2 .and. .not.paramset) then
         write(nfout,'(" !PP nqitg = ",i5)') nqitg
         write(nfout,'(" !PP nqitg(",i3,") = ",i5)') it, nqitg_sp(it)
         write(nfout,'(" !PP -- il3, nmm_il3 : mm_il3(:) --")')
         do il3 = 1, lcmax+1
            if(nmm_il3(il3) > 0) write(nfout,'(" !PP  ",i3,i5," : ",12i3)') il3,nmm_il3(il3),(mm_il3(i,il3),i=1,nmm_il3(il3))
         end do
      else if(ipripp >= 3 .and. paramset) then
         write(nfout,'(" !PP nqitg = ",i5)') nqitg
         write(nfout,'(" !PP nqitg(",i3,") = ",i5)') it, nqitg_sp(it)
      end if
    end subroutine wd_nqitg_sp_etc

  end subroutine rd_qrsps_then_iqitg_and_qitgft
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q_END ORG_Parallel
! ==================================================================================================

! ================================== added by K. Tagami ==================== 11.0
  subroutine rd_qrsps_iqitg_and_qitgft_soc( is_gncpp,nfp,it,nfout,gr_l, &
             	&                            ngshell,ngshell_range,paramset,mmt )
    integer, intent(in) :: is_gncpp
    integer, intent(in)       :: nfp, it, nfout,ngshell
    integer, intent(in), dimension(ngshell,2) :: ngshell_range
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical, intent(in)       :: paramset
    integer, intent(inout) :: mmt

    integer :: il1, il2, tau1, tau2, il3, tmin, l3s, l3l,  mm, i, igs
    integer :: nrc, mord, iprippex, iw, nq, id

    real(kind=DP) :: dummy
    character(len=32):: cdummy
    integer, allocatable, dimension(:)       :: nmm_il3 !d(lcmax+1)
    integer, allocatable, dimension(:,:)     :: mm_il3 !d(nqitg_sp,lcmax+1)
    real(kind=DP), allocatable, dimension(:,:) :: qrsps_mm !d(nmesh(it),nqitg_sp(it))

    integer :: kj1, kj2
    integer :: jval1, jval2

    mm = 0
    iprippex = ipripp
    if(.not.paramset) then
       qij = 0.d0; qvij = 0.d0
       if(allocated(nmm_il3)) deallocate(nmm_il3)
       allocate(nmm_il3(lcmax+1)); nmm_il3 = 0
       if(allocated(mm_il3)) deallocate(mm_il3)
       allocate(mm_il3(nqitg_sp(it),lcmax+1)); mm_il3 = 0
       if(allocated(qrsps_mm)) deallocate(qrsps_mm)
       allocate(qrsps_mm(mmesh,nqitg_sp(it)))
       iprippex = iprippex-1
    end if

    Loop_kj1 : Do kj1 =1, nums_of_angmom_on_atomtype( it )
       il1 = get_lp1_in_AngMomList( pot_has_soc(it),kj1 )
       jval1 = get_jph_in_AngMomList( pot_has_soc(it),kj1 )
       if(iloc(it) == il1) cycle

       Loop_tau1 : do tau1 = 1, itau(kj1,it)

          Loop_kj2 : Do kj2 =1, nums_of_angmom_on_atomtype( it )
             il2 = get_lp1_in_AngMomList( pot_has_soc(it),kj2 )
             jval2 = get_jph_in_AngMomList( pot_has_soc(it),kj2 )

             if ( (iloc(it)==il2) .or. (ivanl(kj1,it)==0.and.ivanl(kj2,it)==0)) cycle

             tmin = 1
!!!             if( il1 == il2) tmin = tau1
             if( kj1 == kj2 ) tmin = tau1

             Loop_tau2 : do tau2 = tmin, itau(kj2,it)
                l3s = abs(il1-il2)
                l3l = abs(il1+il2) - 2

                Loop_L3 : do il3 = l3s+1, l3l+1, 2
                   if(il3-1 > lcmax ) cycle
                   if(ipripp >= 3) &
                        & write(nfout,'(" !PP il1, tau1, il2, tau2, il3 = ",5i5)') &
                        & il1, tau1, il2, tau2, il3

                   mm = mm + 1
                   if(.not.paramset) iqitg( kj1,tau1,kj2,tau2,il3,it) = mm + mmt

                   if(mype == 0) call read_nrc_mord(nfp,nfout, kord,&
                        &            mm,kj1,tau1,kj2,tau2,il3,nrc,mord,iprippex) ! ->(nrc,mord)
                   if(paramset .and. mype == 0) then
                      if(nrc > 0 ) then
                         read(nfp,*) (dummy,i=0,mord)
                         read(nfp,*) (dummy,i=nrc+1,nmesh(it))
                      else if(nrc == 0) then
                         read(nfp,*) (dummy,i=1,nmesh(it))
                      end if
                   else if(.not.paramset) then
                      nmm_il3(il3) = nmm_il3(il3)+1
                      mm_il3(nmm_il3(il3),il3) = mm
                      if(mype == 0) then
                         if(nrc > 0 ) then
                            read(nfp,*) (copsc(i),i=0,mord)
                            read(nfp,*) (qrsps(i),i=nrc+1,nmesh(it))

! --------------
!                            up, down spin compoment ??????? dou yatte readin ???
! -------

#ifdef _QRSPS_WITH_L_
#else
                            if(il3 > 1) then
                               do i = nrc+1, nmesh(it)
                                  qrsps(i) = qrsps(i)/(radr(i)**(il3-1))
                               end do
                            end if
#endif
                            call cnvrtp(nmesh(it),1,nrc,il3+1,mord,copsc,radr,qrsps)

                         else if(nrc == 0 ) then
                            read(nfp,*) (qrsps(i),i=1,nmesh(it))
#ifdef _QRSPS_WITH_L_
#else
                            if(il3 > 1) then
                               do i = 1, nmesh(it)
                                  qrsps(i) = qrsps(i)/(radr(i)**(il3-1))
                               end do
                            end if
#endif
                         endif


                         if(sw_berry_phase == ON .or. sw_wannier == ON &
                              &       .or. sw_fef==ON ) then
                            if(is_gncpp == pp_GNCPP2 .and. ae_wavefunctions_are_detected(it)) then
                               do i=1,nmesh(it)
                                  qrs(i) = (phir_ae(i,kj1,tau1)*phir_ae(i,kj2,tau2) &
                                  &        -phir(i,kj1,tau1)*phir(i,kj2,tau2) ) &
                                  &         *radr(i)**(il3-1)
                               end do
                            else
                               qrs = qrsps
                            end if
                         end if

                      end if

! ================== Added by K. Tagami =================== 10.1
                         if ( sw_LinearResponse == ON ) then
                            if (is_gncpp == 2 .and. ae_wavefunctions_are_detected(it)) then
                               do i=1,nmesh(it)
                                  qrs(i) = (phir_ae(i,il1,tau1)*phir_ae(i,il2,tau2) - &
                                       & phir(i,il1,tau1)*phir(i,il2,tau2))*radr(i)**(il3-1)
                               end do
                            else
                               qrs = qrsps
                            end if
                         endif
! ========================================================= 10.1

                      if (npes > 1) &
                           & call mpi_bcast(qrsps,nmesh(it),mpi_double_precision,0,&
                           &                MPI_CommGroup,ierr )

                      qrsps_mm(:,mm) = qrsps(:)

!                      if (is_gncpp == pp_PAW) then
!                         if(.not.paramset .and. ipaw(it) == ON) &
!                              &  qrspspw(1:mmesh,mm+mmt)=qrsps(1:mmesh)
!                      end if
                         if(.not.paramset) qrspspw(1:mmesh,mm+mmt)=qrsps(1:mmesh)

                      if( (sw_berry_phase == ON .or. sw_wannier==ON .or. &
                           &  sw_fef == ON) .and.npes > 1) &
                           & call mpi_bcast(qrs,nmesh(it),mpi_double_precision,0,&
                           &                MPI_CommGroup,ierr)

! ======================== Added by K. Tagami =========== 10.1
                      if ( sw_LinearResponse == ON .and.npes > 1 ) &
                           & call mpi_bcast(qrs,nmesh(it),mpi_double_precision,0,&
                           &                MPI_CommGroup,ierr)
! ======================================================= 10.1

                      if(sw_berry_phase == ON) then
                         if(is_gncpp == pp_GNCPP1) then
                            call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_BP &
                                     &,1.d0,mm+mmt,qitg_BP,wkx,wky)
                         else if(is_gncpp == pp_GNCPP2) then
                            call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_BP &
                                     &,1.d0,mm+mmt,qitg_BP,wkx,wky)
                         end if
                         if(ipripp >= 3) then
                            write(nfout,'(" !PP i=",i3," q_BP=",e25.12)') &
                                 &            mm+mmt,qitg_BP(1,mm+mmt)
                            write(nfout,'(" !PP ",5x," q_BP/pi4=",e25.12)') &
                                 &            qitg_BP(1,mm+mmt)/PAI4
                         end if
                      end if

                      if(sw_wannier == ON) then
                         ! debug
                         if(ipripp >= 2) write(nfout,*) &
                              & 'rd_qrsps_then_iqitg_and_qitgft: nwght_wan=',nwght_wan
                         ! end debug

                         do iw=1,nwght_wan
                            if(is_gncpp == pp_GNCPP1) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_wan(iw) &
                                     &,1.d0,mm+mmt,qitg_wan(1,iw),wkx,wky)
                            else if(is_gncpp == pp_GNCPP2) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_wan(iw) &
                                     &,1.d0,mm+mmt,qitg_wan(1,iw),wkx,wky)
                            end if
                            if(ipripp >= 3) then
                               write(nfout,'(" !PP i=",i3,"iw=",i3," q_wan=",e25.12)')&
                                    &          mm+mmt,iw,qitg_wan(mm+mmt,iw)
                               write(nfout,'(" !PP ",11x," q_wan/pi4=",e25.12)') &
                                    &          qitg_wan(mm+mmt,iw)/PAI4
                            end if
                         end do
                      end if

! ================================= Added by K. Tagami ======== 10.1
                      if ( sw_LinearResponse == ON ) then
                         ! debug
!                         write(6,*) 'rd_qrsps_then_iqitg_and_qitgft: nwght_wan=',nwght_wan                         ! end debug

                         ! ---
                         Do nq=1, nmax_q_LR
                            if ( is_gncpp == 1 ) then
                               call qitgft( 1, nmax_q_plus_G_LR, mmesh, &
                                    & nqitg, nmesh(it), &
                                    & ipri, nfout, il3-1, radr, wos, qrsps, &
                                    & norm_q_plus_G_LR(1,nq), &
                                    & 1.d0, mm+mmt, &
                                    & qitg_LR(1,1,nq), &
                                    & wkx,wky )
                            else if ( is_gncpp == 2 ) then
                               call qitgft( 1, nmax_q_plus_G_LR, mmesh, &
                                    & nqitg, nmesh(it), &
                                    & ipri, nfout, il3-1, radr, wos, qrs, &
                                    & norm_q_plus_G_LR(1,nq), &
                                    & 1.d0, mm+mmt, &
                                    & qitg_LR(1,1,nq), &
                                    & wkx,wky )
                            end if
                         End do
                      end if
! ================================================================= 10.1

                      if(sw_fef == ON) then
                         ! debug
                           write(6,*) 'rd_qrsps_then_iqitg_and_qitgft: numef_fef=',numef_fef
                         ! end debug
                         do id=1,numef_fef
                            if(is_gncpp == 1) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrsps,dk_fef(id) &
                                     &,1.d0,mm+mmt,qitg_fef(1,id),wkx,wky)
                            else if(is_gncpp == 2) then
                               call qitgft(1,1,mmesh,nqitg,nmesh(it),ipri &
                                     &,nfout,il3-1,radr,wos,qrs,dk_fef(id) &
                                     &,1.d0,mm+mmt,qitg_fef(1,id),wkx,wky)
                            end if
                            if(ipripp >= 3) then
                               write(nfout,'(" !PP i=",i3,"id=",i3," q_fef=",e25.12)')   mm+mmt,id,qitg_fef(mm+mmt,id)
                               write(nfout,'(" !PP ",11x," q_fef/pi4=",e25.12)') qitg_fef(mm+mmt,id)/PAI4
                            end if
                         end do
                      end if

                      if (il3 == 1 .and. il1 == il2) then
                         call qij_qvij_from_qrsps_etc_soc( it,tau1,tau2,kj1,kj2 )
                      endif

                   end if

                end do Loop_L3
             end do Loop_tau2
          end do Loop_kj2

       end do Loop_tau1
    end do Loop_kj1

    mmt = mmt + mm
    if(paramset) then
       nqitg = mmt
       nqitg_sp(it) = mm
    end if

    call wd_nqitg_sp_etc()
!!$    if(.not.paramset) write(6,'(" nqitg = ",i5)') nqitg !K.Mae print->write

    if(.not.paramset .and. iflag_ft == 1) call qitgft_mm() ! contained here

!    if(.not.paramset .and. sw_hybrid_functional==ON .and.  m_PP_include_vanderbilt_pot(it)/=SKIP.and.sw_rspace_hyb==OFF) &
!      & call qitgft_qmk(it,nmm_il3,mm_il3,qrsps_mm,lcmax,h) ! in b_PseudoPotential_EXX

    if(itpcc(it) == 1) then
       if(mype == 0) then
          call read_chgpc_nrc_mord(iprippex,nfp,nfout,chgpc(it),nrc,mord)
          ! --- > rhpcr,  copsc is a temporary uesd parameter
          if(nrc > 0) then
             if(paramset) then
                read(nfp,*) (dummy,i=0,mord)
                read(nfp,*) (dummy,i=nrc+1,nmesh(it))
             else
                read(nfp,*) (copsc(i),i=0,mord)
                read(nfp,*) (rhpcr(i),i=nrc+1,nmesh(it))
                call cnvrtp(mmesh,1,nrc,2,mord,copsc,radr,rhpcr)
             end if
          else if(nrc == 0) then
             if(paramset) then
                read(nfp,*) (dummy, i=1,nmesh(it))
             else
                read(nfp,*) (rhpcr(i),i=1,nmesh(it))
             end if
          endif
          ! <---- rhpcr
       end if
       if(npes > 1) call mpi_bcast(chgpc(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(.not.paramset) then
          if(npes > 1) call mpi_bcast(rhpcr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
          if(is_gncpp == pp_PAW .and. ipaw(it) == ON) rhpcrpw(:,it)=rhpcr(:)
       end if
    endif

    if(is_gncpp == pp_PAW) then
       if(mype == 0) then
          read(nfp,*) cdummy
          if(ipaw(it) == ON) then
             read(nfp,*) chgcr(it)
             write(nfout,400) chgcr(it)
          else
             read(nfp,*) dummy
          end if
400       format(' ',' chgcr = ',f12.8)
    ! --- > rhpcr,  copsc is a temporary uesd parameter
          if(paramset.or.ipaw(it)/=ON) then
             read(nfp,*) (dummy, i=1,nmesh(it))
          else
             read(nfp,*) (rhcorpw(i,it),i=1,nmesh(it))
          end if
        ! <---- rhpcr
       end if
       if(npes > 1 .and. ipaw(it) == ON) then
          call mpi_bcast(chgcr(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
          if(.not.paramset) call mpi_bcast(rhcorpw(:,it),nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
       end if

       if(itpcc(it) == 0 .and. ipaw(it)==ON) then
          if(.not.paramset) rhpcrpw(:,it)=0.d0  !rhcorpw(:,it) !0.d0      !2009.4.18
       end if
    end if

    if(.not.paramset) then
       if(sw_hybrid_functional == OFF .or. sw_rspace_hyb == ON) then
         deallocate(qrsps_mm)
         deallocate(mm_il3,nmm_il3)
       endif
    end if
  contains

    subroutine qitgft_mm()
      integer :: mm0, i, il3, mm, mmp, n
      real(kind=DP) :: qitg_sh, qitg_diff_sh, gabs, f2
!!$      integer :: mmp1,mmp2,mmp3,mmp4
!!$      real(kind=DP) :: qitg_sh1, qitg_sh2, qitg_sh3, qitg_sh4
!!$      integer, parameter :: NFOLDING = 4 ! {2|4}
!!$      integer :: mc

      integer, dimension(:), allocatable       :: ip_largerequal_1 ! d(ngshell)

      allocate(ip_largerequal_1(ngshell))
      if(ipripp >= 3) then
         write(nfout,'(" !qitgft_mm  ngshell = ",i8)') ngshell
         write(nfout,'(" !qitgft_mm --- gshell, ip_largerequal_1(=n), radr, G*radr(n-1), G*radr(n)")')
         write(nfout,'(" !PP ngshell_range")')
         write(nfout,'(" !PP ngshell, range1 - range2, nelement, gr")')
         do i = 1, ngshell
            write(nfout,'(" !PP ",i8, i8," - ",i8, i8, f20.8)') i, ngshell_range(i,1)&
                 & ,ngshell_range(i,2),ngshell_range(i,2)-ngshell_range(i,1)+1,gr_l(ngshell_range(i,1))
         end do
      end if

      do i = 1, ngshell
         gabs = gr_l(ngshell_range(i,1))
         if(gabs < DELTA) then
            ip_largerequal_1(i) = nmesh(it)+1
         else
            ip_largerequal_1(i) = ceiling(nmesh(it) - dlog(rmax(it)*gabs)/h(it))
         end if

         if(ipripp >= 2) then
            n = ip_largerequal_1(i)
            if(n <= 1 .or. n > nmesh(it)) then
               if(ipripp >= 3 .or.(ipripp == 2 .and. &
                    & (i <= min(30,ngshell) .or. (i >= max(ngshell-29,1))))) then
                  if(n == 1) then
                     write(nfout,'(" !qitgft_mm : ",i5, i8, "          ",f16.8)') i, n, radr(n)
                  else
                     write(nfout,'(" !qitgft_mm : ",i5, i8)') i, n
                  end if
               end if
            else
               write(nfout,'(" !qitgft_mm : ",i5, i8, 3f16.8)') i,n,radr(n),gabs*radr(n-1),gabs*radr(n)
            end if
         end if
      end do

      mm0 = 0
      do i = 1, it-1
         mm0 = mm0 + nqitg_sp(i)
      end do
      if(ipripp >= 2) write(nfout,'(" !PP mm0 = ", i5)') mm0


      do il3 = 1, lcmax+1
         if(nmm_il3(il3) <= 0) cycle
!!$         mc = nmm_il3(il3)/NFOLDING
         if(istress == 0) then
            if(ipripp >= 2) then
               write(nfout,'(" !PP il3 = ", i5, " nmm = ",i5)') il3, nmm_il3(il3)
               do mm = 1, nmm_il3(il3)
                  write(nfout,'(" !PP mm = ", i5)') mm_il3(mm,il3)
               end do

               write(nfout,'(" qrsps_mm <<qitgft_mm>>")')
               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  write(nfout,'(" il3, mm, mmp = ",3i8)') il3,mm,mmp
                  do i = 1, min(10, nmesh(it))
                     write(nfout,'(" !PP qrsps_mm: ",i6,d16.8)') i, qrsps_mm(i,mmp)
                  end do
                  do i = max(nmesh(it)-9,1),nmesh(it)
                     write(nfout,'(" !PP qrsps_mm: ",i6,d16.8)') i, qrsps_mm(i,mmp)
                  end do
               end do
            end if

            do i = 1, ngshell
               gabs = gr_l(ngshell_range(i,1))
               wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
               call dsjnvn(il3-1,nmesh(it),ip_largerequal_1(i),wkx,wky)

               if(ipripp >= 2) then
                  if(i <= min(5,ngshell)) then
                     write(nfout,'(" !PP ishell = ",i5)') i
                     do mm = 1, min(30, nmesh(it))
                        write(nfout,'(" !PP wkx,wky = ",i5,2d16.8)') mm, wkx(mm),wky(mm)
                     end do
                     do mm = max(nmesh(it)-29,1), nmesh(it)
                        write(nfout,'(" !PP wkx,wky = ",i5,2d16.8)') mm, wkx(mm),wky(mm)
                     end do
                  end if
               end if

               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  qitg_sh = 0.d0
                  do n = 1, nmesh(it)
                     qitg_sh = qitg_sh + wos(n)*qrsps_mm(n,mmp)*wky(n)
                  end do
                  qitg_sh = qitg_sh/univol*PAI4
                  do igs = ngshell_range(i,1), ngshell_range(i,2)
                     qitg_l(igs,mmp+mm0) = qitg_sh
                  end do
               end do

!!$               if(NFOLDING == 2) then
!!$                  do mm = 1, mc*NFOLDING, NFOLDING
!!$                     mmp1 = mm_il3(mm,il3)
!!$                     mmp2 = mm_il3(mm+1,il3)
!!$                     qitg_sh1 = 0.d0
!!$                     qitg_sh2 = 0.d0
!!$                     do n = 1, nmesh(it)
!!$                        f2 = wos(n)*wky(n)
!!$                        qitg_sh1 = qitg_sh1 + f2*qrsps_mm(n,mmp1)
!!$                        qitg_sh2 = qitg_sh2 + f2*qrsps_mm(n,mmp2)
!!$                     end do
!!$                     qitg_sh1 = qitg_sh1/univol*PAI4
!!$                     qitg_sh2 = qitg_sh2/univol*PAI4
!!$                     do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                        qitg_l(igs,mmp1+mm0) = qitg_sh1
!!$                        qitg_l(igs,mmp2+mm0) = qitg_sh2
!!$                     end do
!!$                  end do
!!$               else if(NFOLDING == 4) then
!!$                  do mm = 1, mc*NFOLDING, NFOLDING
!!$                     mmp1 = mm_il3(mm,il3)
!!$                     mmp2 = mm_il3(mm+1,il3)
!!$                     mmp3 = mm_il3(mm+2,il3)
!!$                     mmp4 = mm_il3(mm+3,il3)
!!$                     qitg_sh1 = 0.d0
!!$                     qitg_sh2 = 0.d0
!!$                     qitg_sh3 = 0.d0
!!$                     qitg_sh4 = 0.d0
!!$                     do n = 1, nmesh(it)
!!$                        f2 = wos(n)*wky(n)
!!$                        qitg_sh1 = qitg_sh1 + f2*qrsps_mm(n,mmp1)
!!$                        qitg_sh2 = qitg_sh2 + f2*qrsps_mm(n,mmp2)
!!$                        qitg_sh3 = qitg_sh3 + f2*qrsps_mm(n,mmp3)
!!$                        qitg_sh4 = qitg_sh4 + f2*qrsps_mm(n,mmp4)
!!$                     end do
!!$                     qitg_sh1 = qitg_sh1/univol*PAI4
!!$                     qitg_sh2 = qitg_sh2/univol*PAI4
!!$                     qitg_sh3 = qitg_sh3/univol*PAI4
!!$                     qitg_sh4 = qitg_sh4/univol*PAI4
!!$                     do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                        qitg_l(igs,mmp1+mm0) = qitg_sh1
!!$                        qitg_l(igs,mmp2+mm0) = qitg_sh2
!!$                        qitg_l(igs,mmp3+mm0) = qitg_sh3
!!$                        qitg_l(igs,mmp4+mm0) = qitg_sh4
!!$                     end do
!!$                  end do
!!$               end if
!!$               do mm = mc*NFOLDING+1,nmm_il3(il3)
!!$                  mmp = mm_il3(mm,il3)
!!$                  qitg_sh = 0.d0
!!$                  do n = 1, nmesh(it)
!!$                     qitg_sh = qitg_sh + wos(n)*qrsps_mm(n,mmp)*wky(n)
!!$                  end do
!!$                  qitg_sh = qitg_sh/univol*PAI4
!!$                  do igs = ngshell_range(i,1), ngshell_range(i,2)
!!$                     qitg_l(igs,mmp+mm0) = qitg_sh
!!$                  end do
!!$               end do
!!$
            end do
         else if(istress == 1) then
            do i = 1, ngshell
               gabs = gr_l(ngshell_range(i,1))
               wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
               call dsjnvn(il3-1,nmesh(it),ip_largerequal_1(i),wkx,wky)
               if(il3-1 == 0) then
                  wkz1 = 0.d0
               else
                  call dsjnvn(il3-2,nmesh(it),ip_largerequal_1(i),wkx,wkz1)
               end if
               call dsjnvn(il3,nmesh(it),ip_largerequal_1(i),wkx,wkz2)
               do mm = 1, nmm_il3(il3)
                  mmp = mm_il3(mm,il3)
                  qitg_sh = 0.d0
                  qitg_diff_sh = 0.d0
                  do n = 1, nmesh(it)
                     f2 = wos(n)*qrsps_mm(n,mmp)
                     qitg_sh = qitg_sh + f2*wky(n)
                     qitg_diff_sh = qitg_diff_sh + f2 * radr(n) &
                          & / (2.d0*(il3-1)+1.d0) * (il3*wkz1(n)-(il3)*wkz2(n))
                  enddo
                  qitg_sh = qitg_sh/univol*PAI4
                  qitg_diff_sh = qitg_diff_sh/univol*PAI4
                  do igs = ngshell_range(i,1), ngshell_range(i,2)
                     qitg_l(igs,mmp+mm0) = qitg_sh
                     qitg_diff_l(igs,mmp+mm0) = qitg_diff_sh
                  end do
               end do
            end do
         end if
      end do

      if(ipripp >= 2) then
         write(nfout,'(" qitg_l <<qitgft_mm>>")')
         do mm = 1, nqitg
            write(nfout,'(" mm = ", i8)') mm
            write(nfout,'(6d12.4)') &
                 & (qitg_l(igs,mm),igs=ista_kngp, min(ista_kngp+60,iend_kngp))
         end do
      end if
      deallocate(ip_largerequal_1)

    end subroutine qitgft_mm

    subroutine wd_nqitg_sp_etc()
      integer :: il3
      if(ipripp >= 2 .and. .not.paramset) then
         write(nfout,'(" !PP nqitg = ",i5)') nqitg
         write(nfout,'(" !PP nqitg(",i3,") = ",i5)') it, nqitg_sp(it)
         write(nfout,'(" !PP -- il3, nmm_il3 : mm_il3(:) --")')
         do il3 = 1, lcmax+1
            if(nmm_il3(il3) > 0) write(nfout,'(" !PP  ",i3,i5," : ",12i3)') il3,nmm_il3(il3),(mm_il3(i,il3),i=1,nmm_il3(il3))
         end do
      else if(ipripp >= 3 .and. paramset) then
         write(nfout,'(" !PP nqitg = ",i5)') nqitg
         write(nfout,'(" !PP nqitg(",i3,") = ",i5)') it, nqitg_sp(it)
      end if
    end subroutine wd_nqitg_sp_etc

  end subroutine rd_qrsps_iqitg_and_qitgft_soc
!$$#endif
! =========================================================================== 11.0

  subroutine qij_qvij_from_qrsps_etc(it,t1,t2,il)
    integer, intent(in) :: it,t1,t2,il
    integer :: i
    real(kind=DP) :: sij, svij
    sij = 0.d0
    svij = 0.d0
    do i = 1, nmesh(it)
       sij  = sij + wos(i)*qrsps(i)
       svij = svij + wos(i)*qrsps(i)*vlocr(i)
    enddo
    qij( t1,t2,il) = sij
    qvij(t1,t2,il) = svij
    if(t1 /= t2) then
       qij( t2,t1,il) = sij
       qvij(t2,t1,il) = svij
    endif
  end subroutine qij_qvij_from_qrsps_etc
!#####################################################################

! ==================================== added by K. Tagami ================== 11.0
  subroutine qij_qvij_from_qrsps_etc_soc(it,t1,t2,kj1,kj2)
    integer, intent(in) :: it,t1,t2, kj1, kj2
    integer :: i, is

!    real(kind=DP) :: sij, svij
!    sij = 0.d0
!    svij = 0.d0

!    Do is=1, ndim_spinor
!       do i = 1, nmesh(it)
!          sij  = sij + wos(i)*qrsps(i,is)
!          svij = svij + wos(i)*qrsps(i,is)*vlocr(i)
!       enddo
!    End do

!    qij( t1,t2,kj1) = sij
!    qvij(t1,t2,kj1) = svij
!    if (t1 /= t2) then
!       qij( t2,t1,kj1) = sij
!       qvij(t2,t1,kj1) = svij
!    endif

    call phase_error_with_msg(6,'kt: Not supported ',__LINE__,__FILE__)

  end subroutine qij_qvij_from_qrsps_etc_soc
! ============================================================================= 11.0

  subroutine wd_index_arrays_etc(nfout,it)
    integer, intent(in) :: nfout,it
    integer :: lmt1, lmt2, ip

    write(nfout,'(" !PP")')
    if(ipripp>=3) then
       write(nfout,'(" !PP -- dion ---")')
       do lmt1 = 1, ilmt(it)
          write(nfout,'(" !PP ",i3," ",9f8.5,5(/," !PP     ",9f8.5))') lmt1 &
               &               ,(dion(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
       enddo

       write(nfout,'(" !PP --  q   ---")')
       do lmt1 = 1, ilmt(it)
          write(nfout,'(" !PP ",i3," ",9f8.5,5(/," !PP    ",9f8.5))') lmt1 &
               &               ,(q(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
       enddo
       write(nfout,'(" !PP -- index_lmt1_lmt2 --")')
       write(nfout,'(" !PP ",10(" (",2i3,")"))') &
            &(index_lmt1_lmt2(ip,it,1),index_lmt1_lmt2(ip,it,2),&
            & ip=1,n_non0_lmtxlmt(it))
       write(nfout,'(" !PP -- dion_indp --")')
       write(nfout,'(" !PP ",10f8.5)') (dion_indp(ip,it),ip = 1, n_non0_lmtxlmt(it))
       write(nfout,'(" !PP -- q_indp --")')
       write(nfout,'(" !PP ",10f8.5)') (q_indp(ip,it),ip = 1, n_non0_lmtxlmt(it))
       write(nfout,'(" !PP -- weight --")')
       write(nfout,'(" !PP ",10i8)') (w_non0_lmtxlmt(ip,it),ip=1,n_non0_lmtxlmt(it))
    else if(ipripp == 2) then
       write(nfout,'(" !PP -- lmt1, lmt2, dion, q,  weight --")')
       do ip = 1, n_non0_lmtxlmt(it)
          write(nfout,'(" !PP ",2i8,2f10.5, i8)') index_lmt1_lmt2(ip,it,1),index_lmt1_lmt2(ip,it,2) &
               &  , dion_indp(ip,it), q_indp(ip,it), w_non0_lmtxlmt(ip,it)
       end do
    end if
    if(sw_orb_popu == ON) then
       write(nfout,'(" !PP -- q_phi ---")')
       write(nfout,'(" !PP ",10f8.5)') (q_phi(ip,it),ip=1,ilmt_phi(it))
    end if
  end subroutine wd_index_arrays_etc

  subroutine m_PP_wd_variables(nfout)
    integer, intent(in) :: nfout
    if(ipripp>=1) write(nfout,300) etot1
300 FORMAT(' !PP ','TOTAL ENERGY CORRECTION FROM PS     = ',d20.10)
  end subroutine m_PP_wd_variables

  subroutine m_PP_tell_lmtt_l_m_tau(lmt,it,ilmtt,il,im,tau,nspher)
    integer, intent(in)  :: lmt, it
    integer, intent(out) :: ilmtt,il,im,tau,nspher

    ilmtt = lmtt(lmt,it)
    il    = ltp( lmt,it)
    im    = mtp( lmt,it)
    tau   = taup(lmt,it)
    nspher = (il-1)**2 + im
  end subroutine m_PP_tell_lmtt_l_m_tau

  subroutine m_PP_tell_lmtt_l_m_tau_phi(lmt,it,ilmtt,il,im,tau,nspher)
    integer, intent(in)  :: lmt, it
    integer, intent(out) :: ilmtt,il,im,tau,nspher

    ilmtt = lmtt_phi(lmt,it)
    il    = ltp_phi( lmt,it)
    im    = mtp_phi( lmt,it)
    tau   = taup_phi(lmt,it)
    nspher = (il-1)**2 + im
  end subroutine m_PP_tell_lmtt_l_m_tau_phi

  subroutine m_PP_tell_lmtt_l_m_tau_add(lmt,it,ilmtt,il,im,tau,nspher)
    integer, intent(in)  :: lmt, it
    integer, intent(out) :: ilmtt,il,im,tau,nspher

    ilmtt = lmtt_add(lmt,it)
    il    = ltp_add( lmt,it)
    im    = mtp_add( lmt,it)
    tau   = 1
    nspher = (il-1)**2 + im
  end subroutine m_PP_tell_lmtt_l_m_tau_add

  subroutine m_PP_tell_lmtt_l_m_tau_pao(lmt,it,ilmtt,il,im,tau,nspher)
    integer, intent(in)  :: lmt, it
    integer, intent(out) :: ilmtt,il,im,tau,nspher

    ilmtt = lmtt_pao(lmt,it)
    il    = ltp_pao( lmt,it)
    im    = mtp_pao( lmt,it)
    tau   = taup_pao( lmt,it)
    nspher = (il-1)**2 + im
  end subroutine m_PP_tell_lmtt_l_m_tau_pao

  subroutine m_PP_make_index_lmtt_2_dl2p(nfout,paramset)
    logical, intent(in)   :: paramset
    integer, intent(in)   :: nfout

    integer :: mt,ma,it,il1,tau1,im1,lmt1,ii,lmt2,il2,im2,tau2,jj,n,ia,mm,il3
    real(kind=DP), pointer, dimension(:,:,:) :: cr2
    integer, pointer, dimension(:,:,:)       :: isph2
    integer, pointer, dimension(:,:)         :: mmt2

    if(.not.paramset) then
       allocate(cr2(16,16,6))
       allocate(isph2(16,16,6))
       allocate(mmt2(16,16))
!!$ASASASASAS
       cr2 = 0; isph2 = 0; mmt2 = 0
!!$ASASASASAS
       call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)
    endif

    mt = 0;  nlmt = 0
    do it = 1, ntyp
       mm = 0
       do il1 = 1, lpsmax(it)
          if(il1 == iloc(it)) cycle
          do tau1 = 1, itau(il1,it)
             do im1 = 1, il1*2-1
                mm = mm + 1
                lmtt(mm,it) = mm + mt
             end do
          end do
       end do
       mt = mt + mm
       if(nlmt < mm) nlmt = mm
       do lmt1 = 1, mm
          il1  = ltp(lmt1,it)
          im1  = mtp(lmt1,it)
          tau1 = taup(lmt1,it)
          ii   = (il1-1)**2 + im1
          do lmt2 = 1, mm
             il2  = ltp(lmt2,it)
             im2  = mtp(lmt2,it)
             tau2 = taup(lmt2,it)
             jj   = (il2-1)**2 + im2
             if(.not.paramset) then
                il2p(lmt1,lmt2,it) = mmt2(ii,jj)
                do n = 1, il2p(lmt1,lmt2,it)
                   isph(lmt1,lmt2,n,it) = isph2(ii,jj,n)
                   dl2p(lmt1,lmt2,n,it) = cr2  (ii,jj,n)
                enddo
             end if
          end do
       end do
    end do
    nlmtt = mt

    ma = 0
    do ia = 1, natm
       mm = 0
       it = ityp(ia)
       do il1 = 1, lpsmax(it)
          if(il1 == iloc(it)) cycle
          do tau1 = 1, itau(il1,it)
             do im1 = 1, il1*2-1
                mm = mm + 1
                if(.not.paramset) lmta(mm,ia) = mm + ma
             end do
          end do
       end do
       ma = ma + mm
    end do

    nlmta = ma
    if(ipripp >= 2) then
       write(nfout,'(" !PP nlmt  = ",i8)') nlmt
       write(nfout,'(" !PP nlmta = ",i8)') nlmta
       write(nfout,'(" !PP nlmtt = ",i8)') nlmtt
    end if

!--*--*
    IF(ipripp >= 2 .and. .not.paramset) THEN
       WRITE(NFOUT,'(" !PP N, LMTT, LTP, MTP, TAUP")')
       DO IT = 1,NTYP
          WRITE(NFOUT,510) IT
510       FORMAT(' !PP ',' TYPE = ',I4)
          WRITE(NFOUT,520) (N,LMTT(N,IT),LTP(N,IT),MTP(N,IT),TAUP(N,IT),N=1,ILMT(IT))
520       FORMAT((' !PP ',5I5))
       end DO
       WRITE(NFOUT,'(" !PP LMT1,LMT2,N,ISPH,DL2P,ISPH")')
       DO IT = 1,NTYP
          WRITE(NFOUT,540) IT
540       FORMAT(' !PP ',' TYPE = ',I4)
          DO LMT1 = 1,ILMT(IT)
             DO LMT2 = 1,ILMT(IT)
                WRITE(NFOUT,570) LMT1,LMT2
570             FORMAT(' !PP ',2I4)
                WRITE(NFOUT,580) (N,ISPH(LMT1,LMT2,N,IT), &
                     &   DL2P(LMT1,LMT2,N,IT), N=1,IL2P(LMT1,LMT2,IT)   )
580             FORMAT((' !PP ',8X,2I5,F12.8))
             end DO
          end DO
       end DO
       WRITE(NFOUT,'(" !PP LMT1,LMTA")')
       DO IA = 1,NATM
          WRITE(NFOUT,610) IA
610       FORMAT(' !PP ',' ATOM = ',I4)
          IT = ITYP(IA)
          WRITE(NFOUT,620) (LMT1,LMTA(LMT1,IA),LMT1=1,ILMT(IT))
620       FORMAT((' !PP ',5(2I4,3X)))
       end DO
       WRITE(NFOUT,630)
630    FORMAT(' !PP ',' IT, LMT1, LMT2, IL1, IL2, IL3, TAU1, TAU2, IQITG ')
       DO IT=1,NTYP
          IF(IVAN(IT).NE.1) cycle
          DO LMT1=1,ILMT(IT)
             IL1 = LTP(LMT1,IT)
             TAU1= TAUP(LMT1,IT)
             DO LMT2=LMT1,ILMT(IT)
                IL2 = LTP(LMT2,IT)
                TAU2= TAUP(LMT2,IT)
                WRITE(NFOUT,670) (IT,LMT1,LMT2,IL1,IL2,IL3,TAU1,TAU2,&
                     &    IQITG(IL1,TAU1,IL2,TAU2,IL3+1,IT), &
                     &    IL3 = ABS(IL1-IL2),IL1+IL2-2,2)
670             FORMAT(' !PP ',9I5)
             end DO
          end DO
       end DO
    END IF

    if(.not.paramset) then
       deallocate(cr2)
       deallocate(isph2)
       deallocate(mmt2)
    endif
  end subroutine m_PP_make_index_lmtt_2_dl2p

  subroutine m_PP_make_index_lmtt_phi(nfout,paramset)
    logical, intent(in)   :: paramset
    integer, intent(in)   :: nfout

    integer mm, mt, il1, tau1, im1, n, i, ia, ma, it

    mt = 0; ma = 0
    do it = 1, ntyp
       mm = 0
       L_loop2: do il1 = 1, lpsmax(it)
          T_loop2: do tau1 = 1, itau(il1,it)
             !!$do i=1,norbital(it)
             !!$   if(il1 == l_orb(i,it)+1 .and. tau1 == t_orb(i,it)) then
             !!$      do im1 = 1, il1*2-1
             !!$         mm = mm + 1
             !!$         lmtt_phi(mm,it) = mm + mt
             !!$      end do
             !!$   end if
             !!$end do
             do i=1,num_projectors
                if(it /= proj_attribute(i)%ityp) cycle
                if(il1 == proj_attribute(i)%l+1 .and. &
                 & tau1 == proj_attribute(i)%t) then
                   do im1 = 1, il1*2-1
                      mm = mm + 1
                      lmtt_phi(mm,it) = mm + mt
                   end do
                end if
             end do
          end do T_loop2
       end do L_loop2
       mt = mt + mm
       if(nlmt_phi < mm) nlmt_phi = mm
       do ia = 1, natm
          !!$if(if_pdos(ia) == OFF) cycle
!!$ASASASASAS
!!$          if(iproj_group(ia) == 0 ) cycle
!!$ASASASASAS
          mm = 0
          if(ityp(ia) /= it) cycle
          !!$L_loop3: do il1 = 1, lpsmax(it)
          !!$   T_loop3: do tau1 = 1, itau(il1,it)
          !!$      do i=1,norbital(it)
          !!$         if(il1 == l_orb(i,it)+1 .and. tau1 == t_orb(i,it)) then
          !!$            do im1 = 1, il1*2-1
          !!$               mm = mm + 1
          !!$               if(.not.paramset) lmta_phi(mm,ia) = mm + ma
          !!$            end do
          !!$         end if
          !!$      end do
          !!$   end do T_loop3
          !!$end do L_loop3
          L_loop3: do il1 = 1, lpsmax(it)
             T_loop3: do tau1 = 1, itau(il1,it)
                do i=1,num_projectors
                   if(it /= proj_attribute(i)%ityp) cycle
                   if(il1 == proj_attribute(i)%l+1 .and. &
                    & tau1 == proj_attribute(i)%t) then
                      do im1 = 1, il1*2-1
                         mm = mm + 1
                         if(.not.paramset) lmta_phi(mm,ia) = mm + ma
                      end do
                   end if
                end do
             end do T_loop3
          end do L_loop3
          ma = ma + mm
       end do
    end do
    nlmta_phi = ma
    nlmtt_phi = mt
    if(ipripp >= 2) then
       write(nfout,'(" !PP nlmt_phi  = ",i8)') nlmt_phi
       write(nfout,'(" !PP nlmta_phi = ",i8)') nlmta_phi
       write(nfout,'(" !PP nlmtt_phi = ",i8)') nlmtt_phi
    end if

  end subroutine m_PP_make_index_lmtt_phi

  subroutine m_PP_make_index_lmtt_add(nfout,paramset)
    logical, intent(in)   :: paramset
    integer, intent(in)   :: nfout

    integer mm, mt, il1, im1, ia, ma, it

! debug
    if(ipripp >= 2) then
       write(nfout,*) ' !PP debug(m_PP_make_index_lmtt_add)'
       write(nfout,*) ' !PP alloc(lmtt_add) = ',allocated(lmtt_add)
       write(nfout,*) ' !PP alloc(lmta_add) = ',allocated(lmta_add)
    end if
! end debug

    mt = 0; ma = 0
    do it = 1, ntyp
       mm = 0
       L_loop2: do il1 = 1, lpsmax(it)
          if(il1 /= iloc(it)) cycle L_loop2
          do im1 = 1, il1*2-1
             mm = mm + 1
             lmtt_add(mm,it) = mm + mt
          end do
       end do L_loop2
       mt = mt + mm
       if(nlmt_add < mm) nlmt_add = mm
       do ia = 1, natm
          mm = 0
          if(ityp(ia) /= it) cycle
          L_loop3: do il1 = 1, lpsmax(it)
             if(il1 /= iloc(it)) cycle L_loop3
             do im1 = 1, il1*2-1
                mm = mm + 1
                if(.not.paramset) lmta_add(mm,ia) = mm + ma
             end do
          end do L_loop3
          ma = ma + mm
       end do
    end do
    nlmta_add = ma
    nlmtt_add = mt
    if(ipripp >= 1) then
       write(nfout,'(" !PP nlmt_add  = ",i8)') nlmt_add
       write(nfout,'(" !PP nlmta_add = ",i8)') nlmta_add
       write(nfout,'(" !PP nlmtt_add = ",i8)') nlmtt_add
    end if

  end subroutine m_PP_make_index_lmtt_add

  subroutine m_PP_make_index_lmtt_pao(nfout,paramset)
    logical, intent(in)   :: paramset
    integer, intent(in)   :: nfout

    integer :: mt,ma,it,il1,tau1,im1,ia,mm,lmta1
    integer, allocatable :: irank_pao(:)
    integer, allocatable :: ind(:)
    integer :: mm_ma, lmt1, ib1, ib2

! ========================= added by K. Tagami =========================== 11.0
#ifdef forsafe
    integer :: iwei_ia
                                   ! This is introduced in order to overcome
                                   ! a compiler bug ????
#endif
! ======================================================================== 11.0

    if(.not.paramset) then
       allocate(irank_pao(nlmta_pao))
       if(kimg==1) then
          allocate(ind(2*nlmta_pao))
       else
          allocate(ind(nlmta_pao))
       end if
    end if

    mt = 0; ma = 0; nlmt_pao = 0
    do it = 1, ntyp
       mm = 0
       do il1 = 1, lpsmax(it)
          do tau1 = 1, itau(il1,it)
             if(irank(il1,tau1,it)<=0) cycle
             do im1 = 1, il1*2-1
                mm = mm + 1
                lmtt_pao(mm,it) = mm + mt
             end do
          end do
       end do
       mt = mt + mm
       if(nlmt_pao < mm) nlmt_pao = mm
       do ia = 1, natm
          mm = 0
          if(ityp(ia) /= it) cycle

! =================================== added by K. Tagami ============== 11.0
#ifdef forsafe
          iwei_ia = iwei(ia)
#endif
! ===================================================================== 11.0

          do il1 = 1, lpsmax(it)
             do tau1 = 1, itau(il1,it)
                if(irank(il1,tau1,it)<=0) cycle
                do im1 = 1, il1*2-1
                   mm = mm + 1
                   if(.not.paramset) then
                      mm_ma = mm + ma
                      lmta_pao(mm,ia) = mm_ma
                      irank_pao(mm_ma) = irank(il1,tau1,it)

! ============================= modified by K. Tagami ================= 11.0
!                      ind(mm_ma) = mm_ma
!                      if(iwei(ia)/=1) ind(mm_ma) = -ind(mm_ma)
!
#ifdef forsafe
                      if ( iwei_ia == 1 ) then
                         ind(mm+ma) =   mm+ma
                      else
                         ind(mm+ma) = -(mm+ma)
                      end if
#else
                      ind(mm_ma) = mm_ma
                      if(iwei(ia)/=1) ind(mm_ma) = -ind(mm_ma)
#endif
! ===================================================================== 11.0

                   end if
                end do
             end do
          end do
          ma = ma + mm
       end do
    end do
    nlmta_pao = ma
    nlmtt_pao = mt
    if(ipripp >= 1) then
       write(nfout,'(" !PP nlmt_pao  = ",i8)') nlmt_pao
       write(nfout,'(" !PP nlmta_pao = ",i8)') nlmta_pao
       write(nfout,'(" !PP nlmtt_pao = ",i8)') nlmtt_pao
    end if

    if(.not.paramset) then
!! Debug
!      do lmta1=1,nlmta_pao
!         write(nfout,'("lmta1=",i5," ind=",i5)') lmta1, ind(lmta1)
!      end do
!! End Debug
       call heap_sort_wrt_rank(nlmta_pao,ind,irank_pao)

       ma = nlmta_pao
       do lmta1=1,ma
          if(ind(lmta1)<0) then
             do mm=ma,lmta1+1,-1
                ind(mm+1) = ind(mm)
             end do
             ind(lmta1+1) =  ind(lmta1)
             ind(lmta1)   = -ind(lmta1)
             ma = ma+1
          end if
       end do
       do mm=1,ma
          lmta1 = ind(mm)
          if(lmta1>0) then
             ibpao(lmta1) = mm
          end if
       end do

       deallocate(irank_pao);       deallocate(ind)

       num_wfc_pao = 0
       Do ia=1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt_pao(it)
             lmta1 = lmta_pao(lmt1,ia)
             ib1    = ibpao(lmta1)
             ib2   = 0
             if (ib1 < 0) then
                ib1 = abs(ib1);    ib2 = ib1 + 1
             end if
             if (ib1 > 0 ) num_wfc_pao = num_wfc_pao +1
             if (ib2 > 0 ) num_wfc_pao = num_wfc_pao +1
          end do
       end do
!       write(*,*) "ll : num_wfc_pao = ", num_wfc_pao
    end if

  end subroutine m_PP_make_index_lmtt_pao

  subroutine heap_sort_wrt_rank(n,ind,irank)
    integer, intent(in) :: n
    integer, intent(inout) :: ind(n)
    integer, intent(inout) :: irank(n)

    integer :: i,j,k,itmp
    integer :: indt,inds,irant,irans

    !*** initiation: heapfy ind and eps ***
    do i=2,n
       indt = ind(i)
       irant = irank(i)
       j=i
10     itmp=j/2
       if(irank(itmp).ge.irant) go to 20
       ind(j) = ind(itmp)
       irank(j) = irank(itmp)
       j=itmp
       if(j.gt.1) go to 10
20     ind(j) = indt
       irank(j) = irant
    end do
    !*** to be fully sort ***
    do k=n-1,1,-1
       inds = ind(1)
       ind(1) = ind(k+1)
       ind(k+1) = inds
       irans = irank(1)
       irank(1) = irank(k+1)
       irank(k+1) = irans
       indt = ind(1)
       irant = irank(1)
       j=1
       itmp=2
30     if(itmp.gt.k) go to 40
       if(itmp.lt.k) then
          if(irank(itmp+1).gt.irank(itmp)) itmp=itmp+1
       end if
       if(irant.ge.irank(itmp)) go to 40
       ind(j) = ind(itmp)
       irank(j) = irank(itmp)
       j=itmp
       itmp=j*2
       go to 30
40     ind(j) = indt
       irank(j) = irant
    end do

  end subroutine heap_sort_wrt_rank

  subroutine m_PP_find_maximum_l(l_max)
    integer, intent(out) :: l_max
    integer it, lmt1, l
                                                  __TIMER_SUB_START(1058)
    l_max = 0
                                                  __TIMER_DO_START(1080)
    do it = 1, ntyp
!!$       if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do lmt1 = 1, ilmt(it)
          l = ltp(lmt1,it)
          if(l_max < l) l_max = l
       end do
    end do
                                                  __TIMER_DO_STOP(1080)
                                                  __TIMER_SUB_STOP(1058)
  end subroutine m_PP_find_maximum_l

  integer function m_PP_include_vanderbilt_pot(it)
    integer, intent(in) :: it

    integer :: tmodnrm, il
    tmodnrm = NO
    if(ivan(it) /= 0) then
       Loop_iloc :do il = 1, lpsmax(it)
          if(ivanl(il,it) == 1) then
             tmodnrm = YES
             exit
          end if
       end do Loop_iloc
    end if
    m_PP_include_vanderbilt_pot = tmodnrm
  end function m_PP_include_vanderbilt_pot


  integer function m_PP_max_of_itau()
    integer it, n, im
    n = 0
    do it = 1, ntyp
       do im = 1, lpsmax(it)
          if(itau(im,it) > n) n = itau(im,it)
       end do
    end do
    m_PP_max_of_itau = n
  end function m_PP_max_of_itau

  subroutine m_PP_set_mmesh
    mmesh = maxval(nmesh)
    !!$ print '(" ! mmesh = ", i5)',mmesh
  end subroutine m_PP_set_mmesh

  subroutine m_PP_set_nloc(nfout)
    integer, intent(in) :: nfout

! ======================== modified by K. Tagami =================== 11.0
!    nloc = maxval(lpsmax)
!
    nloc = maxval( nums_of_angmom_on_atomtype)
    nloc_for_l = maxval( lpsmax )
! ================================================================== 11.0

    if(ipripp >= 2) write(nfout,'(" !PP nloc = ",i5)') nloc
    !!$ print '(" ! nloc = ",i5)', nloc
  end subroutine m_PP_set_nloc

  subroutine m_PP_set_m_non0_lmtxlmt
    m_non0_lmtxlmt = maxval(n_non0_lmtxlmt)
    !!$ print '(" ! m_non0_lmtxlmt = ",i5)', m_non0_lmtxlmt
  end subroutine m_PP_set_m_non0_lmtxlmt

  subroutine m_PP_tell_iorb_ia_l_m_tau(iorb,ia1,il,im,tau)
    integer, intent(in) :: iorb
    integer, intent(out) :: ia1,il,im,tau

    integer :: ia,it,ilmt

    do ia=1,natm
       !!$if(if_pdos(ia) == OFF) cycle
!!$ASASASASAS
!!$       if(iproj_group(ia) == OFF) cycle
!!$ASASASASAS
       it = ityp(ia)
       do ilmt = 1,ilmt_phi(it)
          if(iorb == lmta_phi(ilmt,ia)) then
             ia1 = ia
             il  = ltp_phi(ilmt,it)
             im  = mtp_phi(ilmt,it)
             tau = taup_phi(ilmt,it)
          end if
       end do
    end do

  end subroutine m_PP_tell_iorb_ia_l_m_tau

  subroutine m_PP_tell_iorb_lmt(iorb,lmt)
    integer, intent(in) :: iorb
    integer, intent(out) :: lmt

    integer :: ia,it,ilmt

    ATOM: do ia=1,natm
       !!$if(if_pdos(ia) == OFF) cycle
!!$ASASASASAS
!!$       if(iproj_group(ia) == 0) cycle
!!$ASASASASAS
       it = ityp(ia)
       do ilmt = 1,ilmt_phi(it)
          if(iorb == lmta_phi(ilmt,ia)) then
             lmt = ilmt
             lmt = lmtt_phi( ilmt,it )
             exit ATOM
          end if
       end do
    end do ATOM

  end subroutine m_PP_tell_iorb_lmt

  subroutine m_PP_tell_iorb_lmtt(iorb,lmtt)
    integer, intent(in) :: iorb
    integer, intent(out) :: lmtt

    integer :: ia,it,ilmt

    ATOM: do ia=1,natm
       it = ityp(ia)
       do ilmt = 1,ilmt_phi(it)
          if(iorb == lmta_phi(ilmt,ia)) then
             lmtt = lmtt_phi(ilmt,it)
             exit ATOM
          end if
       end do
    end do ATOM
  end subroutine m_PP_tell_iorb_lmtt

  subroutine m_PP_make_qorb(nfout)
  integer, intent(in) :: nfout
  integer :: mm,lmt,ia,it,i,ilmta

    allocate(qorb(nlmta_phi))

    do ia=1,natm
       !!$if(if_pdos(ia) == OFF) cycle
!!$ASASASASAS
!!$       if(iproj_group(ia) == 0) cycle
!!$ASASASASAS
       it = ityp(ia)
       do lmt=1,ilmt_phi(it)
          ilmta=lmta_phi(lmt,ia)
          qorb(ilmta) = q_phi(lmt,it)
       end do
    end do

    if(ipripp>=2) then
       write(nfout,*) ' !PP --- qorb --- nlmta_phi = ',nlmta_phi
       do i=1,nlmta_phi
          write(nfout,*) ' !PP i=',i,' qorb = ',qorb(i)
       end do
    end if

  end subroutine m_PP_make_qorb

  subroutine set_projector_radius_if_undef( nfout, it )
    integer, intent(in) :: nfout,it
    integer :: i, il1, tau1, mm

    mm = maxval( wf_nrc(:,:,it) )

    do il1 = 1, lpsmax(it)
       do tau1 = 1,itau(il1,it)
          do i=1, num_projectors
             if ( proj_attribute(i)%radius_was_defined ) cycle
             if (it /= proj_attribute(i)%ityp) cycle
             if (il1 == proj_attribute(i)%l+1 .and. &
                  & tau1 == proj_attribute(i)%t) then

!               mm = wf_nrc( il1, tau1, it )

                proj_attribute(i)%radius = radr( mm )
                proj_attribute(i)%radius_was_defined = .true.
                if ( mype == 0 ) then
                   write(nfout,'(A,I5,A,F15.8)') "** The radius of ", i, &
                        &         "-th projector is set to ", proj_attribute(i)%radius
                endif
             endif
          end do
        end do
     end do
  end subroutine set_projector_radius_if_undef

  subroutine make_phi_window_parameter(nfout,it)
  integer, intent(in) :: nfout,it
  integer :: i, il1, tau1

  do il1 = 1, lpsmax(it)
     do tau1 = 1,itau(il1,it)
        do i=1,num_projectors
           if(it /= proj_attribute(i)%ityp) cycle
           if(il1 == proj_attribute(i)%l+1 .and. &
            & tau1 == proj_attribute(i)%t) then
              rc_phi(il1,tau1,it) = proj_attribute(i)%radius
              k_phi(il1,tau1,it) = 1.d0/(proj_attribute(i)%fwhm*sqrt(log(2.d0)))

! === ASMS === 2014/05/29
              k_phi(il1,tau1,it) = sqrt( k_phi(il1,tau1,it) )
! === ASMS === 2014/05/29

              exit
           end if
        end do
     end do
  end do

  end subroutine make_phi_window_parameter

! ================================ added by K. Tagami ================== 11.0
  subroutine make_phi_window_parameter_soc(nfout,it)
    integer, intent(in) :: nfout,it
    integer :: i, il1, tau1

    call phase_error_with_msg(nfout,"kt: not implemented yet.",__LINE__,__FILE__)

    do il1 = 1, lpsmax(it)
       do tau1 = 1,itau(il1,it)
          do i=1,num_projectors
             if(it /= proj_attribute(i)%ityp) cycle
             if(il1 == proj_attribute(i)%l+1 .and. &
                  & tau1 == proj_attribute(i)%t) then
                rc_phi(il1,tau1,it) = proj_attribute(i)%radius
                k_phi(il1,tau1,it) = 1.d0/(proj_attribute(i)%fwhm*sqrt(log(2.d0)))

! === ASMS === 2014/05/29
                k_phi(il1,tau1,it) = sqrt( k_phi(il1,tau1,it) )
! === ASMS === 2014/05/29

                exit
             end if
          end do
       end do
    end do

  end subroutine make_phi_window_parameter_soc
! ===================================================================== 11.0

  subroutine m_PP_rd_window_param(nfout)
  integer, intent(in) :: nfout

  integer :: i, iret, rint, ip
  integer :: f_selectTop, f_selectBlock, f_selectParentBlock
  integer :: f_selectFirstTableLine, f_selectNextTableLine
  integer :: f_getIntValue, f_getRealValue
  real(kind=DP) :: dret

  ! Tags
  character(len("kappa")),parameter        :: tag_kappa        = "kappa"
  character(len("rcproj")),parameter       :: tag_rcproj       = "rcproj"
  character(len("structure")),parameter    :: tag_structure    = "structure"
  character(len("element_list")),parameter :: tag_element_list = "element_list"
  character(len("id")),parameter           :: tag_id           = "id"
  character(len("no")),parameter           :: tag_no           = "no"

  if(sw_use_add_proj == OFF) return

  if(.not.allocated(kappa)) allocate(kappa(ntyp))
  if(.not.allocated(rcproj)) allocate(rcproj(ntyp)); rcproj=-1.d0

  iret = f_selectTop()
  if( f_selectBlock( tag_structure) == 0) then
     if( f_selectBlock( tag_element_list) == 0) then
        do i = 1, ntyp
           if(ipripp>=2) write(nfout,*) ' !PP debug(window_param): i=',i
           if (i == 1) then
              if(f_selectFirstTableLine() /= 0) then
                 if(ipripp>=2) write(nfout,*) ' !PP debug(window_param): FirstTableLine exit'
                 exit
              end if
           else
              if(f_selectNextTableLine() /= 0) then
                 if(ipripp>=2) write(nfout,*) ' !PP debug(window_param): NextTableLine exit'
                 exit
              end if
           end if
           ip = i
           iret = f_getIntValue(tag_id,rint)
           if(iret /= 0 ) iret = f_getIntValue(tag_no,rint)
           if(iret == 0 .and. rint > 0 .and. rint <= ntyp) ip = rint
           if(f_getRealValue(tag_kappa,dret,'') == 0) kappa(ip) = dret
           if(f_getRealValue(tag_rcproj,dret,'bohr') == 0) rcproj(ip) = dret
        end do
        iret = f_selectParentBlock()
     end if
     iret = f_selectParentBlock()
  end if

  if(ipripp>=2)then
     write(nfout,*) ' !PP **  it   kappa   rcproj **'
     do i=1,ntyp
        write(nfout,'(" !PP ",i3,2(1x,f10.5))') i,kappa(i),rcproj(i)
     end do
  end if

  end subroutine m_PP_rd_window_param

  subroutine m_PP_cnstrct_crotylm(nfout)
    use m_Crystal_Structure, only : op,nopr
    use m_Ionic_System, only : napt
    integer, intent(in) :: nfout

    integer :: l1max,mmax,nsph
    integer :: lmt,ia,ja,it,ilmta,ilmtt,il,im,tau,mm,isph,jsph
    integer :: iopr,lmt1
    real(kind=DP), allocatable :: crotylm(:,:,:)
    integer, allocatable :: iylm(:,:,:)
    integer, allocatable :: nylm(:,:)

! debug
!    write(nfout,*) '=== In m_PP_cnstrct_crotylm ==='
! end debug

    l1max = nloc
    mmax  = 2*l1max-1
    nsph = l1max**2
!   allocate(crotylm(mmax,nsph,nopr))
!   allocate(iylm(mmax,nsph,nopr))
!   allocate(nylm(nsph,nopr))
    allocate(crotylm(mmax,nsph,nopr+af)) !ASMS
    allocate(iylm(mmax,nsph,nopr+af))    !ASMS
    allocate(nylm(nsph,nopr+af))         !ASMS

!   call get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,op)
    call get_crotylm(l1max,mmax,nsph,nopr+af,crotylm,iylm,nylm,op) !ASMS

    !debug
    !  do iopr=1,nopr
    !     do isph=1,nsph
    !        write(nfout,'(2(1x,i3))') iopr,isph
    !        do mm=1,nylm(isph,iopr)
    !           write(nfout,'(i3,"=>",f20.8)') iylm(mm,isph,iopr),crotylm(mm,isph,iopr)
    !        end do
    !     end do
    !  end do
    !end debug

!   allocate(nrorb(nlmta_phi,nopr))
!   allocate(irorb(2*nloc+1,nlmta_phi,nopr))
!   allocate(crorb(2*nloc+1,nlmta_phi,nopr))
    allocate(nrorb(nlmta_phi,nopr+af))          !ASMS
    allocate(irorb(2*nloc+1,nlmta_phi,nopr+af)) !ASMS
    allocate(crorb(2*nloc+1,nlmta_phi,nopr+af)) !ASMS

    do ia=1,natm
       it=ityp(ia)
       do lmt=1,ilmt_phi(it)
          call m_PP_tell_lmtt_l_m_tau_phi(lmt,it,ilmtt,il,im,tau,isph)
          ilmta = lmta_phi(lmt,ia)

!         do iopr=1,nopr
          do iopr=1,nopr+af !ASMS
             ja=napt(ia,iopr)
             if(ja>natm) ja=ja-natm
             nrorb(ilmta,iopr) = nylm(isph,iopr)
             do mm=1,nylm(isph,iopr)
! debug
!                write(nfout,'(4(1x,i3))') ia,lmt,iopr,mm
! end debug
                jsph = iylm(mm,isph,iopr)
                lmt1 = get_lmt(it,jsph,tau)
! debug
!                write(nfout,'(2(1x,i3))') jsph,lmt1
! end debug
                irorb(mm,ilmta,iopr) = lmta_phi(lmt1,ja)
                crorb(mm,ilmta,iopr) = crotylm(mm,isph,iopr)
             end do
          end do
       end do
    end do

    deallocate(crotylm)
    deallocate(iylm)
    deallocate(nylm)

! debug
!    write(nfout,'("== Summary ==")')
!    do iopr=1,nopr
!       write(nfout,*) 'iopr=',iopr
!       do ilmta=1,nlmta_phi
!          write(nfout,*) 'ilmta=',ilmta
!          do mm=1,nrorb(ilmta,iopr)
!             write(nfout,'(i3,"=>",f20.8)') irorb(mm,ilmta,iopr), crorb(mm,ilmta,iopr)
!          end do
!       end do
!    end do
! end debug

    !stop 'debug(m_PP_cnstrct_crotylm)'

  contains

    integer function get_lmt(it,isph,itau)
      integer, intent(in) :: it,isph,itau

      integer :: lmt,jsph,jtau

      do lmt=1,ilmt_phi(it)
         jsph = (ltp_phi(lmt,it)-1)**2+mtp_phi(lmt,it)
         jtau = taup_phi(lmt,it)
         if(jsph .eq. isph .and. jtau .eq. itau) then
            get_lmt = lmt
            exit
         end if
      end do

    end function get_lmt

  end subroutine m_PP_cnstrct_crotylm

  subroutine m_PP_ps_xctype(xctype_is)
    integer, intent(out) :: xctype_is
    if(what_is_ps_xctype() == GGA) then
       xctype_is = GGA
    else
       xctype_is = LDA
    end if
  end subroutine m_PP_ps_xctype

  subroutine m_PP_input_xctype(xctype_is)
    integer, intent(out) :: xctype_is
    if(what_is_input_xctype() == GGA) then
       xctype_is = GGA
    else
       xctype_is = LDA
    end if
  end subroutine m_PP_input_xctype

!!$!BRANCH DEV
!!$  subroutine m_PP_vanderbilt_type0(nfp,it,nfout,gr_l,ngshell,ngshell_range,paramset)
!!$    integer, intent(in)       :: nfp, it, nfout,ngshell
!!$    integer, intent(in), dimension(ngshell,2) :: ngshell_range
!!$    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
!!$    logical, intent(in)       :: paramset
!!$
!!$    logical           :: vall_is_detected
!!$    integer           :: mesh_t
!!$    integer           :: id_sname = -1
!!$    call tstatc0_begin('m_PP_vanderbilt_type0 ',id_sname)
!!$    if(it == 1) mmt = 0
!!$    call mpi_barrier(MPI_CommGroup,ierr)
!!$    call rd_Vloc_psWF_Q_then_chir_qitg        !-(contained here) (potrvd)
!!$    !        -> ps_xctype, etc.
!!$    call make_index_arrays_nlmt2l_m_t         ! (tlmst1)
!!$    if(sw_orb_popu == ON ) then
!!$       call make_index_arrays_nlmt_phi2lmt     ! (tlmst1_phi)
!!$    end if
!!$
!!$    if(.not.paramset) then
!!$       call cnstrct_of_localWF_Dion_and_q     ! betavd
!!$       ! -> betar, dion, q
!!$
!!$       call where_does_WFrhvr_damp(nfout,it,mesh_t)    ! -> mesh_t
!!$       call coulomb_potential_in_Rspace(5)    ! -> vvv()   (ptclmb)
!!$       call atomic_xc_potential_in_Rspace     ! -> vxc()   (ptxc)
!!$       ! jfp = 0 (non-spin system for atomic charge) y
!!$       call vlocr_plus_hartree_xc ! vlocr = vlocr - (vvv + vxc)/r
!!$#ifdef _POT_SMOOTHING_
!!$       call smoothing_vlocr()                 ! vlocr smoothing
!!$#endif
!!$       if(initial_chg == from_PSEUDOPOTENTIAL_FILE) &
!!$       & call ft_valence_charge(it,nfout,gr_l,ngshell,ngshell_range) ! rhvr -> rhvg_l
!!$       if(sw_positron /= OFF) call rd_core_charge_then_ft(nfp,it,nfout,gr_l,paramset) ! rhcr -> rhcg_l
!!$    end if
!!$
!!$    call tstatc0_end(id_sname)
!!$  contains
!!$#ifdef _POT_SMOOTHING_
!!$    subroutine smoothing_vlocr()
!!$      real(kind=DP), parameter :: delta7 = 1.d-7, delta10 = 1.d-10
!!$      integer       :: icheck, i
!!$      real(kind=DP) :: s
!!$      icheck = 0
!!$      do i = nmesh(it),1,-1
!!$         if((rhvr(i)/pai4/radr(i)**2) < delta7) then
!!$            if(dabs(vlocr(i)+ival(it)/radr(i)) > delta10 ) then
!!$               if(ipripp >= 2) then
!!$                  write(nfout,'(" !PP VLOC:",i6,f12.6,3d16.8)') &
!!$                       & i,radr(i),vlocr(i), -ival(it)/radr(i)&
!!$                       &,          vlocr(i) + ival(it)/radr(i)
!!$               end if
!!$               vlocr(i) = -ival(it)/radr(i)
!!$            endif
!!$         else
!!$            icheck = icheck + 1
!!$         endif
!!$         if(icheck > 30) exit
!!$      enddo
!!$      s = 0.d0
!!$      do i = 1, nmesh(it)
!!$         s = s + wos(i)*vlocr(i)*radr(i)**2
!!$      end do
!!$      if(ipripp >= 1) write(nfout,'(" !PP (smoothing) s = ",d20.8)') s
!!$    end subroutine smoothing_vlocr
!!$#endif
!!$
!!$    subroutine vlocr_plus_hartree_xc
!!$      integer :: meshup,i
!!$
!!$#ifdef GGA_CHECK
!!$      if(ipripp >= 1) write(nfout,'(" !PP << vlocr_plus_hartree_xc >>")')
!!$#endif
!!$
!!$      if(nmesh(it) > mesh_t) then
!!$         meshup = mesh_t
!!$      else
!!$         meshup = nmesh(it)
!!$      endif
!!$#ifdef _DEBUG_WRITE_
!!$      if(ipripp >= 1) then
!!$         write(nfout,'(" !PP --- vlocr --- ")')
!!$         write(nfout,'(" !PP ",3d26.18)') (vlocr(i),i=1,nmesh(it))
!!$         write(nfout,'(" !PP -- vlocr before adding vcoulomb and xc potential --")')
!!$      end if
!!$#endif
!!$      call sum_of_vlocr(it,nfout)
!!$      call sum_of_vvv_and_vxc(it,nfout)
!!$
!!$      do i = 1, meshup
!!$         vlocr(i) = vlocr(i) - (vvv(i)+vxc(i))/radr(i)
!!$      enddo
!!$
!!$#ifdef GGA_CHECK
!!$      if(ipripp >= 1) then
!!$         do i = meshup, meshup-100, -1
!!$            write(nfout,'(" !PP Vlocr, a3, a2: ",i6,3f12.6)') i, vlocr(i),vvv(i),vxc(i)
!!$         end do
!!$      end if
!!$#endif
!!$
!!$      if(ipripp >= 2) then
!!$         write(nfout,'(" !PP --- after - hartree - xc --")')
!!$         write(nfout,'(" !PP ",3d30.20)') (vlocr(i),i=1,nmesh(it))
!!$         write(nfout,'(" !PP --- coulomb ---------------")')
!!$         write(nfout,'(" !PP ",6d15.7)') (vvv(i),i=1,nmesh(it))
!!$         write(nfout,'(" !PP --- xc      ---------------")')
!!$         write(nfout,'(" !PP ",6d15.7)') (vxc(i),i=1,nmesh(it))
!!$         write(nfout,'(" !PP --- rhvr    ---------------")')
!!$         write(nfout,'(" !PP ",6d15.7)') (rhvr(i),i=1,nmesh(it))
!!$      end if
!!$
!!$    end subroutine vlocr_plus_hartree_xc
!!$
!!$    subroutine atomic_xc_potential_in_Rspace   ! (ptxc)
!!$      real(kind=DP) :: chgrsq, vxc1, vxc2, rs, zet, yc
!!$      integer       :: jsm = 1  ! 1: non-mag, 2: magnetic
!!$      real(kind=DP), parameter :: delta40 = 1.d-40
!!$
!!$      real(kind=DP) :: exc
!!$      integer       :: i, ncut
!!$#ifdef GGA_CHECK
!!$      if(ipripp >= 2) write(nfout,'(" !PP -- atomic_xc_potential_in_Rspace -- ")')
!!$#endif
!!$      wkx = 0.d0
!!$      if(itpcc(it) == 0) then
!!$         wkx(1:mesh_t) = rhvr(1:mesh_t)
!!$      else if(itpcc(it) == 1) then
!!$         wkx(1:mesh_t) = rhvr(1:mesh_t) + rhpcr(1:mesh_t)
!!$      endif
!!$      if(what_is_ps_xctype() == GGA) then ! -(m_PseudoPotential)
!!$         if(         ps_xctype == "ggapbe " .or. ps_xctype == "GGAPBE " &
!!$              & .or. ps_xctype == "ldapbe " .or. ps_xctype == "LDAPBE " &
!!$              & .or. ps_xctype == "revpbe " .or. ps_xctype == "REVPBE " &
!!$              & .or. ps_xctype == "rpbe   " .or. ps_xctype == "rpbe   ") then
!!$            call gga_grad_alloc2()
!!$            if(ipripp >= 2) &
!!$                 & write(nfout,'(" !PP ps_xctype (just before calling xc_gga_rad = ",a7)') &
!!$                 & ps_xctype
!!$            call xc_gga_rad &
!!$                 & (mmesh,nmesh(it),1,1,nfout,radr,h(it),wos,wkx,ps_xctype &
!!$                 & ,vxc,exc &    ! : output
!!$                 & ,grdnts,wwk)  ! : work
!!$            call gga_grad_dealloc2()
!!$#ifdef GGA_CHECK
!!$            if(ipripp >= 1) write(nfout,'(" !PP --- after xc_gga_rad() ---")')
!!$#endif
!!$            do i = 1, nmesh(it)
!!$#ifdef GGA_CHECK
!!$               if((i < 10 .or. i > nmesh(it)-10+1) .and. ipripp >=1 ) &
!!$                    & write(nfout,'(" !PP ( ",i4," ) radr, vxc = ", 2d18.8)') i, radr(i), vxc(i)
!!$#endif
!!$               vxc(i) = vxc(i)*radr(i)
!!$            end do
!!$         else
!!$            call gga_grad_alloc()
!!$            call get_gradient_of_rho(jsm,mddiff,mmesh,mesh_t,wkx,radr&
!!$                 &,h(it),fdiff,coeff1,coeff2,rho,grdnts) ! -(b_PP)  -> rho, grdnts
!!$#ifdef GGA_CHECK
!!$            if(ipripp>=1) then
!!$               write(nfout,'(" !PP after get_gradient_of_rho ")')
!!$               write(nfout,'(" !PP rad, rho, grdnts(*,1), grdnts(*,2)")')
!!$               do i = 1, mesh_t
!!$                  write(nfout,'(" !PP ( ",i4,") rad, rho, grad, grad2 = ",4d20.8," (get_gradient)")') &
!!$                       &        i, radr(i), rho(i), grdnts(i,1), grdnts(i,2)
!!$               end do
!!$            end if
!!$#endif
!!$
!!$            if(         ps_xctype == 'ldapw91' .or. ps_xctype == 'LDAPW91' &
!!$                 & .or. ps_xctype == 'ldapbek' .or. ps_xctype == 'LDAPBEK' ) grdnts = 0.d0
!!$
!!$            call gga_xc_potential_atomic(jsm,len_ps_xctype,ps_xctype &
!!$                 &,mmesh,mesh_t,rho,radr,grdnts,fdiff,vxc,eps_chg,eps_grad)
!!$            call gga_grad_dealloc()
!!$#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
!!$            do i = mmesh,1,-1
!!$               if(wkx(i)/(PAI4*radr(i)**2) > eps_chg) then
!!$                  ncut = i
!!$                  goto 1001
!!$               end if
!!$            end do
!!$1001        continue
!!$            do i = ncut+1, mmesh
!!$               vxc(i) = 0.d0
!!$            end do
!!$#endif
!!$         endif
!!$
!!$#ifdef GGA_CHECK
!!$         if(ipripp >=1) then
!!$            write(nfout,'(" !PP ps_xctype, xctype = ",a7,2x,a7)') ps_xctype, xctype
!!$            do i = 1, nmesh(it)
!!$               write(nfout,'(" !PP vxc(",i4,") = ",d18.10, " radr(",i4,") = ",d18.10)') i, vxc(i),i, radr(i)
!!$            end do
!!$         end if
!!$#endif
!!$
!!$      else
!!$         vxc = 0.d0
!!$         do i = 1, mesh_t
!!$            vxc1 = 0.d0; vxc2 = 0.d0
!!$            chgrsq = wkx(i)
!!$            if(dabs(chgrsq) > delta40) then
!!$               rs = (3*radr(i)*radr(i)/chgrsq)**(1.d0/3.d0)
!!$               zet = 0.d0
!!$!--*--EXCHANGE-CORRELATION POTENTIALS ARE GIVEN IN RYD. UNITS.
!!$               call xcpot(ps_xctype,vxc1,vxc2,rs,zet,yc)
!!$            endif
!!$            vxc(i) = radr(i)*vxc1*0.5d0
!!$         enddo
!!$      end if
!!$    end subroutine atomic_xc_potential_in_Rspace
!!$
!!$    subroutine coulomb_potential_in_Rspace(nsize)
!!$      integer, intent(in) :: nsize
!!$
!!$      real(kind=DP),allocatable,dimension(:) :: da, db ! d(nsize)
!!$      real(kind=DP)       :: s2, rhs, rhs1, bm
!!$      integer             :: i
!!$
!!$      allocate(da(nsize)); da = 0.d0
!!$      allocate(db(nsize)); db = 0.d0
!!$
!!$      s2 = dlog(rhvr(2)/rhvr(1))/h(it)
!!$      rhs = rhvr(1)
!!$      wkx(1)  = rhs*radr(1)/(s2+1)
!!$      wky(1)  = rhs/s2
!!$      db(1)   = h(it)*rhs*3.d0
!!$      da(1)   = db(1)*radr(1)
!!$      rhs1    = rhs
!!$      do i = 2,3
!!$         rhs     = rhvr(i)
!!$         wkx(i)  = wkx(i-1)+h(it)*(rhs*radr(i)+rhs1*radr(i-1))*0.5d0
!!$         wky(i)  = wky(i-1)+h(it)*(rhs        +rhs1          )*0.5d0
!!$         db(i)   = h(it) *rhs*3.d0
!!$         da(i)   = db(i)*radr(i)
!!$         rhs1    = rhs
!!$      enddo
!!$      do i = 4,mesh_t
!!$         rhs    = rhvr(i)
!!$         db(4)  = h(it)*rhs*3.d0
!!$         da(4)  = db(4)*radr(i)
!!$         wkx(i)=(9*wkx(i-1)-wkx(i-3)+da(4)+2.d0*da(3)-da(2))/8.d0
!!$         wky(i)=(9*wky(i-1)-wky(i-3)+db(4)+2.d0*db(3)-db(2))/8.d0
!!$         da(1)  = da(2)
!!$         db(1)  = db(2)
!!$         da(2)  = da(3)
!!$         db(2)  = db(3)
!!$         da(3)  = da(4)
!!$         db(3)  = db(4)
!!$      enddo
!!$      bm               = wky(mesh_t)
!!$!C--*--COULOMB POTENTIAL RVC
!!$      vvv = 0.d0
!!$      do i = 1,mesh_t
!!$         vvv(i) = wkx(i)+radr(i)*(bm-wky(i))
!!$      enddo
!!$
!!$      deallocate(da); deallocate(db)
!!$    end subroutine coulomb_potential_in_Rspace
!!$
!!$    subroutine cnstrct_of_localWF_Dion_and_q
!!$      integer       :: il1,t1,t2,i,lmt1,il11,lmt2,im1,il22,im2,ip
!!$      real(kind=DP) :: dion_tmp
!!$
!!$      if(ipripp >= 2) then
!!$         write(nfout,'(" !PP -- <<cnstrct_of_localWF_Dion_and_q>> --")')
!!$         write(nfout,'(" !PP it = ",i3)') it
!!$         write(nfout,'(" !PP lpsmax ( ",i3,") = ", i3)') it,lpsmax(it)
!!$         write(nfout,'(" !PP iloc   ( ",i3,") = ", i3)') it,iloc(it)
!!$         write(nfout,'(" !PP nmesh  ( ",i3,") = ", i5)') it,nmesh(it)
!!$         write(nfout,'(" !PP wos(100:101) = ",2d16.8)') wos(100),wos(101)
!!$      end if
!!$
!!$      dion(:,:,it) = 0.d0; q(:,:,it) = 0.d0; if_q_isnotzero(:,:,it) = NO
!!$      loop_il1 :do il1 = 1, lpsmax(it)
!!$         if(ipripp >= 1) then
!!$            write(nfout,'(" !PP il1 = ",i5)') il1
!!$            write(nfout,'(" !PP itau(",i3,",",i3,") = ",i3)') il1,it,itau(il1,it)
!!$         end if
!!$         if(iloc(it) == il1) cycle
!!$         do t1 = 1, itau(il1,it)
!!$            if(ipripp >= 1) then
!!$               write(nfout,'(" !PP phir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, phir(100,il1,t1),phir(101,il1,t1)
!!$               write(nfout,'(" !PP chir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, chir(100,il1,t1),chir(101,il1,t1)
!!$            end if
!!$            do t2 = 1,itau(il1,it)
!!$               s(t1,t2) = 0.d0
!!$               do i = 1, nmesh(it)
!!$                  s(t1,t2) = s(t1,t2) + wos(i)*phir(i,il1,t1)*chir(i,il1,t2)
!!$               enddo
!!$            enddo
!!$         enddo
!!$         if(ipripp >= 1) then
!!$            do t1 = 1, itau(il1,it)
!!$               do t2 = 1, itau(il1,it)
!!$                  write(nfout,'(" !PP s (",i3,",",i3,") = ",d16.8)') t1,t2,s(t1,t2)
!!$               end do
!!$            end do
!!$         end if
!!$         if(ipripp >= 3) write(nfout,'(" !PP -- after t1 loop")')
!!$
!!$         do lmt1 = 1, ilmt(it)
!!$            il11 = ltp(lmt1,it)
!!$            im1  = mtp(lmt1,it)
!!$            t1  = taup(lmt1,it)
!!$            do lmt2 = 1, ilmt(it)
!!$               il22 =  ltp(lmt2,it)
!!$               im2  =  mtp(lmt2,it)
!!$               t2   = taup(lmt2,it)
!!$               if(il11==il1 .and. il22==il1 .and. im1==im2) then
!!$                  q   (lmt1,lmt2,it) = qij(t1,t2,il1)
!!$                  dion(lmt1,lmt2,it) = s(t1,t2) - qvij(t1,t2,il1) &
!!$                       & + eps(il1,t2,it)*qij(t1,t2,il1)
!!$               endif
!!$               if(dabs(q(lmt1,lmt2,it)) > DELTA) if_q_isnotzero(lmt1,lmt2,it) = YES
!!$            enddo
!!$         enddo
!!$         if(ipripp >= 3) write(nfout,'(" !PP -- after lmt1 loop")')
!!$
!!$         do ip = 1, n_non0_lmtxlmt(it)
!!$            lmt1 = index_lmt1_lmt2(ip,it,1)
!!$            il11 =  ltp(lmt1,it)
!!$            if(il11 == il1) then
!!$               im1  =  mtp(lmt1,it)
!!$               t1   = taup(lmt1,it)
!!$               lmt2 = index_lmt1_lmt2(ip,it,2)
!!$               t2  =  taup(lmt2,it)
!!$               q_indp(ip,it) = qij(t1,t2,il1)
!!$               dion_indp(ip,it) = s(t1,t2) - qvij(t1,t2,il1) &
!!$                    & + eps(il1,t2,it)*qij(t1,t2,il1)
!!$            endif
!!$         enddo
!!$
!!$         call matrix_inversion(nfout,ntau,itau(il1,it),s,sinv)  ! s(i,j) = <phir(i)|chir(j)>
!!$
!!$         do t1 = 1, itau(il1,it)
!!$            betar(1:nmesh(it),il1,t1,it) = 0.d0
!!$            do t2 = 1, itau(il1,it)
!!$               do i= 1, nmesh(it)
!!$                  betar(i,il1,t1,it) = betar(i,il1,t1,it) &
!!$                       & + sinv(t2,t1)*chir(i,il1,t2)
!!$               enddo
!!$            enddo
!!$         end do
!!$      end do loop_il1
!!$
!!$      do lmt2 = 2, ilmt(it)
!!$         do lmt1 = 1, lmt2-1
!!$            dion_tmp = (dion(lmt1,lmt2,it) + dion(lmt2,lmt1,it))*0.5d0
!!$            dion(lmt1,lmt2,it) = dion_tmp
!!$            dion(lmt2,lmt1,it) = dion_tmp
!!$         enddo
!!$      enddo
!!$
!!$      if(sw_orb_popu == ON) then
!!$         ! for PDOS
!!$         do lmt1 = 1, ilmt_phi(it)
!!$            t1  = taup_phi(lmt1,it)
!!$            il1 = ltp_phi(lmt1,it)
!!$            q_phi(lmt1,it) = qij(t1,t1,il1)
!!$         end do
!!$      end if
!!$
!!$      if(ipripp >= 1) call wd_index_arrays_etc(nfout,it) ! printout dion, q, index_lmt1_lmt2 etc
!!$    end subroutine cnstrct_of_localWF_Dion_and_q
!!$
!!$    subroutine make_index_arrays_nlmt2l_m_t
!!$      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2
!!$      mm = 0
!!$      do il1 = 1, lpsmax(it)
!!$         if(il1 == iloc(it)) cycle
!!$         do tau1 = 1, itau(il1,it)
!!$            do im1 = 1, il1*2-1
!!$               mm = mm + 1
!!$               ltp( mm,it) = il1
!!$               mtp( mm,it) = im1
!!$               taup(mm,it) = tau1
!!$            enddo
!!$         enddo
!!$      enddo
!!$      ilmt(it) = mm
!!$      if(ipripp >= 1) then
!!$         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
!!$         write(nfout,400) (n,ltp(n,it),mtp(n,it),taup(n,it),n=1,ilmt(it))
!!$400      format((' !PP ',4i5))
!!$      end if
!!$
!!$      call make_index_arrays(nfout,it,paramset)
!!$    end subroutine make_index_arrays_nlmt2l_m_t
!!$
!!$    subroutine make_index_arrays_nlmt_phi2lmt
!!$      integer mm, il1, tau1, im1, n, ilmt1, i
!!$      logical l_exists, t_exists
!!$
!!$      mm = 0
!!$      L_loop: do il1 = 1, lpsmax(it)
!!$         T_loop: do tau1 = 1, itau(il1,it)
!!$            do i=1,norbital(it)
!!$               if(il1 == l_orb(i,it)+1 .and. tau1 == t_orb(i,it)) then
!!$                  do im1 = 1, il1*2-1
!!$                     mm = mm + 1
!!$                     ltp_phi( mm,it) = il1
!!$                     mtp_phi( mm,it) = im1
!!$                     taup_phi(mm,it) = tau1
!!$                  enddo
!!$               end if
!!$            enddo
!!$         enddo T_loop
!!$      enddo L_loop
!!$      ilmt_phi(it) = mm
!!$      if(ipripp >= 1) then
!!$         write(nfout,'(" !PP PDOS: orbital list")')
!!$         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
!!$         write(nfout,400) (n,ltp_phi(n,it),mtp_phi(n,it),taup_phi(n,it),n=1,ilmt_phi(it))
!!$400      format((' !PP ',4i5))
!!$      end if
!!$
!!$    end subroutine make_index_arrays_nlmt_phi2lmt
!!$
!!$    subroutine rd_Vloc_psWF_Q_then_chir_qitg
!!$      integer       :: i
!!$      real(kind=DP) :: dummy
!!$
!!$      call read_natomn_ival_iloc_itpcc()
!!$      call read_ps_xctype()
!!$      call read_alp_cc()
!!$      call read_nmesh_xh_rmax()
!!$
!!$      if(mype == 0) then
!!$         vall_is_detected = check_vall()
!!$         if(vall_is_detected) then
!!$            if(paramset) then
!!$               read(nfp,*) (dummy,i=1,nmesh(it))
!!$            else
!!$               read(nfp,*) (vvv(i),i=1,nmesh(it))
!!$            end if
!!$         endif
!!$         if(paramset) then
!!$            read(nfp,*) (dummy,i=1,nmesh(it))
!!$            read(nfp,*) (dummy,i=1,nmesh(it))
!!$         else
!!$            read(nfp,*) (vlocr(i),i=1,nmesh(it)) ! local potential
!!$            read(nfp,*) (rhvr(i), i=1,nmesh(it)) !
!!$            call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,h(it)) ! -(b_P.P.) -> radr,h(it)
!!$         end if
!!$      end if
!!$      if(npes > 1 .and. .not.paramset) then
!!$         call mpi_bcast(vvv,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(vlocr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(rhvr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(h(it),1,        mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(radr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(vall_is_detected,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$      end if
!!$
!!$
!!$      call decide_lpsmax ! -> (lpsmax)
!!$      if(ipripp >= 2) write(nfout,'(" !PP lpsmax = ",i8)') lpsmax(it)
!!$
!!$      call rd_itau_etc_then_phir_chir   !->(phir,chir):pWF, chir:(e-T-Vloc)*phir
!!$
!!$      if(.not.paramset) then
!!$         call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr&
!!$              & ,wos)
!!$      end if
!!$      call rd_qrsps_then_iqitg_and_qitgft(1,nfp,it,nfout,gr_l,ngshell,ngshell_range,paramset,mmt) ! m_PP
!!$    end subroutine rd_Vloc_psWF_Q_then_chir_qitg
!!$
!!$    subroutine rd_itau_etc_then_phir_chir
!!$      integer       :: il, t1, nrc, i, nm, mord, imax
!!$      real(kind=DP) :: dummy,drc,dk,fac
!!$
!!$      nm = nmesh(it)
!!$      if(ipripp >= 1) write(nfout,'(" !PP lpsmax = ",i8)') lpsmax(it)
!!$      loop_L : do il = 1, lpsmax(it)
!!$         if(mype == 0) &
!!$              & call read_itau_ivanl(nfp,nfout,iloc(it),il,itau(il,it)&
!!$              &   ,ivanl(il,it)) ! -(b_P.P.) ->(itau,ivanl)
!!$         if(npes > 1) then
!!$            call mpi_bcast(itau(il,it), 1,mpi_integer,0,MPI_CommGroup,ierr)
!!$            call mpi_bcast(ivanl(il,it),1,mpi_integer,0,MPI_CommGroup,ierr)
!!$         end if
!!$
!!$         loop_tau : do t1 = 1, itau(il,it)
!!$            if(mype == 0 .and. ipripp >= 2) write(nfout,'(" !PP t1 = ",i8)') t1
!!$            if(mype == 0) call read_tau_eps_nrc_mord(nfp,nfout,kord &
!!$                 &                           ,il,t1,eps(il,t1,it),nrc,mord) ! ->(eps,t1,it)
!!$            if(npes > 1) call mpi_bcast(eps(il,t1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$            if(mype == 0 .and. ipripp >= 2) write(nfout,'(" !PP ",2i4,f12.8,2i6 &
!!$                 & ," : il, ntau, eps, nrc, mord")') il,t1, eps(il,t1,it), nrc,mord
!!$
!!$            if(paramset .and. mype == 0) then
!!$               if(nrc == 0) then
!!$                  read(nfp,*) (dummy,i=1,nm);      read(nfp,*) (dummy,i=1,nm)
!!$               else if(nrc > 0) then
!!$                  read(nfp,*) (dummy,i=0,mord);    read(nfp,*) (dummy,i=nrc+1,nm)
!!$               end if
!!$            else if(.not.paramset) then
!!$               if(mype == 0) then
!!$                  if(nrc == 0) then
!!$                     read(nfp,*) (phir(i,il,t1),i=1,nm)
!!$                     read(nfp,*) (chir(i,il,t1),i=1,nm)
!!$                     chir(1:nm,il,t1) = (chir(1:nm,il,t1)-vlocr(1:nm)) * phir(1:nm,il,t1)
!!$                  else if(nrc > 0) then
!!$                     read(nfp,*) (copsc(i),i=0,mord)
!!$                     read(nfp,*) (phir(i,il,t1),i=nrc+1,nm)
!!$                     call cnvrtp(nm,1,nrc,il,mord,copsc,radr,phir(1,il,t1)) !-(b_P.P._f77.f)
!!$                     call cnvrtc(nm,1,nrc,il-1,mord,copsc,radr,eps(il,t1,it)&
!!$                          & ,phir(1,il,t1),vlocr,chir(1,il,t1))             !-(b_P.P._f77.f)
!!$                     if(vall_is_detected) then
!!$                        chir(nrc+1:nm,il,t1) = (vvv(nrc+1:nm) &
!!$                             & - vlocr(nrc+1:nm))*phir(nrc+1:nm,il,t1)
!!$                     else
!!$                        chir(nrc+1:nm,il,t1) = 0.d0
!!$                     end if
!!$                  end if
!!$               end if
!!$            end if
!!$
!!$         enddo loop_tau
!!$      enddo loop_L
!!$
!!$      if(.not.paramset .and. npes > 1) then
!!$         call mpi_bcast(phir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$         call mpi_bcast(chir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      end if
!!$
!!$      if(.not.paramset .and. sw_orb_popu == 1) then
!!$         call make_phi_window_parameter(nfout,it)
!!$         do il=1,lpsmax(it)
!!$            do t1=1,itau(il,it)
!!$               dk = k_phi(il,t1,it)
!!$               if(dk > 1.d-5) then
!!$                  do i=2,nm
!!$                     if(radr(i) > rc_phi(il,t1,it)) then
!!$                        imax = i-1
!!$                        drc = radr(imax)
!!$                        exit
!!$                     end if
!!$                  end do
!!$                  fac = 1.d0/(1.d0-exp(-dk*drc**2))
!!$                  phirt(1:imax,il,t1,it) = &
!!$                    & fac*(1.d0-exp(-dk*(radr(1:imax)-drc)**2))*phir(1:imax,il,t1)
!!$                  phirt(imax+1:nm,il,t1,it) = 0.d0
!!$               else
!!$                  ! debug
!!$                  if(ipripp >= 2) write(nfout,*) '!PP without a window function'
!!$                  ! end debug
!!$                  phirt(1:nm,il,t1,it) = phir(1:nm,il,t1)
!!$               end if
!!$            end do
!!$         end do
!!$      end if
!!$
!!$      if(paramset) then
!!$         if(it == 1) ntau = itau(1,it)
!!$         do il = 2, lpsmax(it)
!!$            if(ntau < itau(il,it)) ntau = itau(il,it)
!!$         end do
!!$      end if
!!$
!!$    end subroutine rd_itau_etc_then_phir_chir
!!$
!!$    subroutine read_natomn_ival_iloc_itpcc
!!$      integer, parameter  ::    len_str = 80
!!$      character(len=len_str) :: str
!!$
!!$      integer :: natomn
!!$      logical :: comment_statement
!!$
!!$      if(mype == 0) then
!!$         comment_statement = .true.
!!$         rewind(nfp)
!!$         if(ipripp >= 1) write(nfout,'(" !PP READING POTENTIAL FILE ",i8)') nfp
!!$         do while(comment_statement)
!!$            read(nfp,'(a80)') str
!!$            if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
!!$                 & .or. str(1:1) == '%') then
!!$               if(ipripp >= 1) write(nfout,'(a80)') str
!!$            else
!!$               comment_statement = .false.
!!$            endif
!!$         enddo
!!$         read(str,*)      natomn,ival(it),iloc(it),itpcc(it)
!!$         if(ipripp >= 1) write(nfout,110) natomn,ival(it),iloc(it),itpcc(it)
!!$110      FORMAT(' !PP ',I4,f8.4,2I4,'  : NATOMN, IVAL, ILOC, ITPCC ')
!!$         if(iatomn(it) /= natomn) then
!!$            if(ipripp >= 1) write(nfout,'(" !PP iatomn.ne.natomn ",2i8)') iatomn(it),natomn
!!$            stop
!!$         end if
!!$      end if
!!$
!!$      call mpi_bcast(ival(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(iloc(it), 1,mpi_integer,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(itpcc(it),1,mpi_integer,0,MPI_CommGroup,ierr)
!!$      if(ipripp >= 1) then
!!$         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: ival(",i4,") = ",f8.4)') it,ival(it)
!!$         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: iloc(",i4,") = ",i3)') it,iloc(it)
!!$         write(nfout,'(" !PP read_natomn_ival_iloc_itpcc: itpcc(",i3,") = ",i3)') it,itpcc(it)
!!$      end if
!!$
!!$    end  subroutine read_natomn_ival_iloc_itpcc
!!$
!!$    subroutine read_ps_xctype()
!!$      integer  ::                     counter = 0
!!$      character(len=len_ps_xctype) :: ps_xctype0
!!$      logical :: tf
!!$
!!$      counter = counter + 1
!!$      if(mype == 0) then   ! MPI
!!$         read(nfp,'(a7)') ps_xctype0
!!$         call nameck(ps_xctype0,len_ps_xctype)
!!$         if(ipripp >= 1) write(nfout,120) ps_xctype0
!!$120      FORMAT(' !PP ',A7,'   : NAME ')
!!$         if(counter == 1) then
!!$            ps_xctype = ps_xctype0
!!$         else
!!$            call strncmp0(ps_xctype,ps_xctype0,tf)
!!$            if(.not.tf) then
!!$            if(ps_xctype /= ps_xctype0) then
!!$               if(ipripp >= 1) write(nfout,*) ' !PP XC TYPE OF ', it,'TH ATOM IS NOT CORRECT'
!!$               stop
!!$            endif
!!$            call strncmp0(xctype,ps_xctype0,tf)
!!$            if(.not.tf) then
!!$            if(xctype /= ps_xctype0) then
!!$               if(ipripp >= 1) then
!!$                  write(nfout,*) &
!!$                       &' !PP xc-type of ',it,'th atom is not equal to ' &
!!$                       &,'that described in the input file'
!!$                  write(nfout,'(" !PP xctype = ",a7," ps_xctype0 = ",a7)') xctype,ps_xctype0
!!$               end if
!!$            endif
!!$         endif
!!$      end if
!!$      call mpi_bcast(ps_xctype,len_ps_xctype,mpi_character,0,MPI_CommGroup,ierr)
!!$      if(counter == 1 .and. xctype == "nogiven") call m_CtrlP_set_xctype(ps_xctype)
!!$    end subroutine read_ps_xctype
!!$
!!$    subroutine read_alp_cc
!!$      if(mype == 0) then          ! MPI
!!$         read(nfp,*) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
!!$         if(ipripp >= 1) write(nfout,130) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
!!$130      FORMAT(' !PP ',4F12.6,'  :   ALP,CC')
!!$      end if
!!$      call mpi_bcast(alp(1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(alp(2,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(cc(1,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(cc(2,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$    end subroutine read_alp_cc
!!$
!!$    subroutine read_nmesh_xh_rmax
!!$      if(mype == 0) then  ! MPI
!!$         read(nfp,*) nmesh(it),xh(it),rmax(it)
!!$         if(ipripp >= 1) write(nfout,140) nmesh(it),xh(it),rmax(it)
!!$140      format(' !PP ',i4,2f12.6,'  :   nmesh,  xh, rmax')
!!$      end if
!!$      call mpi_bcast(nmesh(it),1,mpi_integer         ,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(xh(it),   1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      call mpi_bcast(rmax(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$    end subroutine read_nmesh_xh_rmax
!!$
!!$    subroutine decide_lpsmax
!!$      character(len=7) :: lpsmaxj
!!$      if(mype == 0) then
!!$         if(iloc(it) == 4) then
!!$            lpsmax(it) = iloc(it)
!!$         else
!!$            lpsmax(it) = 3
!!$         endif
!!$         read(nfp,'(a7)') lpsmaxj
!!$         if(lpsmaxj == 'F-STATE') then
!!$            lpsmax(it) = 4
!!$            if(ipripp >= 1) then
!!$debug T. Y. 2003.02.20
!!$               write(nfout,'(" !PP lpsmaxj == F-STATE")')
!!$debug
!!$            end if
!!$         else
!!$            backspace nfp
!!$         endif
!!$      end if
!!$      call mpi_bcast(lpsmax(it),1,mpi_integer,0,MPI_CommGroup,ierr)
!!$    end subroutine decide_lpsmax
!!$
!!$    logical function check_vall()
!!$      character(len=4) :: vall
!!$      read(nfp,'(a4)') vall
!!$      if(vall == 'VALL') then
!!$         check_vall = .true.
!!$      else
!!$         backspace nfp
!!$         check_vall = .false.
!!$      end if
!!$    end function check_vall
!!$  end subroutine m_PP_vanderbilt_type0
!!$!BRANCH_END DEV

  subroutine m_PP_set_index_arrays1(nfout,mtyp,mqitg,mcritical,mil3,il3 &
       & , maxm,mc &
       & , nqitg_sp, nqitg_sp0,iq2l3,nc,exx_mode)
    ! *****
    ! Coded by T. Yamasaki
    !  This subroutine is translated from the module 'm_ES_Intgr_VlhxcQlm',
    ! at 31 Aug. 2007
    !
    integer, intent(in)                    :: nfout,mtyp,mqitg,mcritical,mil3
    integer, intent(in),  dimension(mil3)  :: il3
    integer, intent(out)                   :: maxm,mc
    integer, intent(out), dimension(ntyp)  :: nqitg_sp,nqitg_sp0
    integer, intent(out), dimension(mqitg) :: iq2l3
    integer, intent(inout), dimension(mcritical,mqitg) :: nc
    logical, intent(in), optional          :: exx_mode

    integer :: it, lmt1, il1, tau1, lmt2, il2, tau2, np, ilm3, l3, iiqitg, m, nqlast
    integer :: lmt2s
    logical :: exx_m
                                                  __TIMER_SUB_START(723)
    exx_m = .false.
    if(present(exx_mode))then
       exx_m = exx_mode
    endif
    maxm = 0
    nqlast = 0
!  -- Revised by T. Yamasaki, 2009/05/28 (pointed out by RIKEN group 2009/02/06)
!!$    nqitg_sp0(1) = 1
    nqitg_sp0 = 1
!  ---------------
                                                  __TIMER_DO_START(838)
    do it = 1, ntyp
       nqitg_sp(it) = 0
       if(m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do lmt1 = 1, ilmt(it)
          il1  = ltp( lmt1, it);   tau1 = taup(lmt1, it)
          lmt2s = lmt1
          if (exx_m) lmt2s=1
          do lmt2 = lmt2s, ilmt(it)
             il2  = ltp( lmt2, it);   tau2 = taup(lmt2, it)
             do np = 1, il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,np,it);   l3 = il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(exx_m)then
                   if (lmt2<lmt1)then
                      iiqitg = iqitg(il2,tau2,il1,tau1,l3+1,it)
                   endif
                endif
                if(iiqitg == 0) cycle
                if(iiqitg > nqitg_sp(it)) nqitg_sp(it)=iiqitg
                iq2l3(iiqitg) = l3
                m = ilm3 - l3*l3
                if(m > maxm) maxm = m
                if(m > mcritical) then
                   if(ipripp >=2 ) write(nfout,'(" !mPP: m, mcritical = ",2i5)') m, mcritical
                end if
                nc(m,iiqitg) = nc(m,iiqitg)+1
             end do
          end do
       end do
       nqitg_sp0(it) = nqlast + 1
       nqlast = nqitg_sp(it)
       if(ipripp >= 2) &
            & write(nfout,'(" !mPP: nqitg_sp0(",i3,") = ",i5," nqitg_sp(",i3,") = ",i5)') &
            & it,nqitg_sp0(it),it, nqitg_sp(it)
    end do
                                                  __TIMER_DO_STOP(838)

    mc = maxval(nc(1:maxm,1:nqitg))

    if(ipripp >= 2) then
       write(nfout,'(" !mPP: nqitg = ",i5)') nqitg
       write(nfout,'(" !mPP: maxm = ",i5)') maxm
       do iiqitg = 1, nqitg
          l3 = iq2l3(iiqitg)
          do m = 1, l3*2+1
             write(nfout,'(" !mPP: nc( ",i3,",",i4," ) = ",i5)') m,iiqitg,nc(m,iiqitg)
          end do
       end do
       write(nfout,'(" !mPP: mc = ", i5)') mc
    end if
                                                  __TIMER_SUB_STOP(723)
  end subroutine m_PP_set_index_arrays1

  subroutine m_PP_set_index_arrays2(nfout,mc,maxm,mqitg,mcritical,mil3,il3,iq2l3&
       & , nc2lmt1,nc2lmt2,nc2n,nc,exx_mode)
    ! *****
    ! Coded by T. Yamasaki
    !  This subroutine is translated from the module 'm_ES_Intgr_VlhxcQlm',
    ! at 31 Aug. 2007
    !
    integer, intent(in)                              :: nfout,mc,maxm,mqitg &
         &                                            , mcritical,mil3
    integer, intent(in),  dimension(mil3)            :: il3
    integer, intent(in),  dimension(mqitg)           :: iq2l3
    integer, intent(out), dimension(mc,maxm,mqitg)   :: nc2lmt1,nc2lmt2,nc2n
    integer, intent(out), dimension(mcritical,mqitg) :: nc
    logical, intent(in), optional                    :: exx_mode

    integer :: it, lmt1, il1, tau1, lmt2, il2, tau2, np, ilm3, l3, iiqitg, m, ip
    integer :: lmt2s
    logical :: exx_m
                                                  __TIMER_SUB_START(724)
    exx_m = .false.
    if(present(exx_mode))then
       exx_m = exx_mode
    endif
    nc = 0
                                                  __TIMER_DO_START(839)
    do it = 1, ntyp
       if( m_PP_include_vanderbilt_pot(it) == SKIP) cycle
       do lmt1 = 1, ilmt(it)
          il1  = ltp( lmt1, it);   tau1 = taup(lmt1, it)
          lmt2s = lmt1
          if (exx_m) lmt2s = 1
          do lmt2 = lmt2s, ilmt(it)
             il2  = ltp( lmt2, it);   tau2 = taup(lmt2, it)
             do np = 1, il2p(lmt1,lmt2,it)
                ilm3 = isph(lmt1,lmt2,np,it);   l3 = il3(ilm3)
                iiqitg = iqitg(il1,tau1,il2,tau2,l3+1,it)
                if(exx_m)then
                   if(lmt2<lmt1)then
                      iiqitg = iqitg(il2,tau2,il1,tau1,l3+1,it)
                   endif
                endif
                if(iiqitg == 0) cycle
                m = ilm3 - l3*l3
                nc(m,iiqitg) = nc(m,iiqitg)+1
                ip = nc(m,iiqitg)
                nc2lmt1(ip,m,iiqitg) = lmt1
                nc2lmt2(ip,m,iiqitg) = lmt2
                nc2n(ip,m,iiqitg) = np
             end do
          end do
       end do
    end do
                                                  __TIMER_DO_STOP(839)
    if(ipripp >= 2) then
       write(nfout,'(" !mPP:  it, lmt1, lmt2, il2p, ilm3, l3, iiqitg, n")')
       np = 0
       do iiqitg = 1, nqitg
          l3 = iq2l3(iiqitg)
          it_loop: do it = 1, ntyp
             if(iiqitg <= nqitg_sp(it)) exit it_loop
          end do it_loop
          if(it > ntyp) it = ntyp
          do m = 1, 2*l3+1
             ilm3 = l3*l3+m
             do ip = 1, nc(m,iiqitg)
                np = np + 1
                write(nfout,'( " !mPP: ",8i5)') &
                     & it, nc2lmt1(ip,m,iiqitg),nc2lmt2(ip,m,iiqitg),nc2n(ip,m,iiqitg),ilm3, l3, iiqitg,np
             end do
          end do
       end do
    end if
                                                  __TIMER_SUB_STOP(724)
  end subroutine m_PP_set_index_arrays2

  subroutine m_PP_dealloc
    if(allocated(ival)) deallocate(ival)
    if(allocated(iloc)) deallocate(iloc)
    if(allocated(lpsmax)) deallocate(lpsmax)
    if(allocated(itpcc)) deallocate(itpcc)
    if(allocated(ilmt)) deallocate(ilmt)
    if(allocated(nmesh)) deallocate(nmesh)
    if(allocated(xh)) deallocate(xh)
    if(allocated(rmax)) deallocate(rmax)
    if(allocated(h)) deallocate(h)
    if(allocated(chgpc)) deallocate(chgpc)
    if(allocated(alp)) deallocate(alp)
    if(allocated(cc)) deallocate(cc)
    if(allocated(n_non0_lmtxlmt)) deallocate(n_non0_lmtxlmt)

    ! deallocate  ivanl, itau, eps, lmtt, ltp, mtp and taup
    if(allocated(ivanl)) call m_PP_dealloc_ps_ntyp

    if(allocated(lmta)) deallocate(lmta)
    if(allocated(q)) deallocate(q)
    if(allocated(if_q_isnotzero)) deallocate(if_q_isnotzero)
    if(allocated(dion)) deallocate(dion)
    if(allocated(isph)) deallocate(isph)
    if(allocated(il2p)) deallocate(il2p)
    if(allocated(dl2p)) deallocate(dl2p)
    if(allocated(prodphi)) deallocate(prodphi)

    if(allocated(index_lmt1_lmt2)) deallocate(index_lmt1_lmt2)
    if(allocated(iqitg)) deallocate(iqitg)
    if(allocated(fqwei)) deallocate(fqwei)
    if(allocated(nlmta1)) deallocate(nlmta1)
    if(allocated(nlmta2)) deallocate(nlmta2)
    if(allocated(nac2ia)) deallocate(nac2ia)
    if(allocated(qrspspw)) deallocate(qrspspw)

    if(allocated(w_non0_lmtxlmt)) deallocate(w_non0_lmtxlmt)
    if(allocated(q_indp)) deallocate(q_indp)
    if(allocated(dion_indp)) deallocate(dion_indp)

    if(allocated(ilmt_add)) deallocate(ilmt_add)
    if(allocated(lmta_add)) deallocate(lmta_add)
    if(allocated(ae_wavefunctions_are_detected)) deallocate(ae_wavefunctions_are_detected)

    if(allocated(nqitg_sp)) deallocate(nqitg_sp)
    if ( allocated(radial_aug_chg) ) deallocate( radial_aug_chg )

    call m_PP_dealloc_psc_qitg_rhpcg()

    !if(allocated(nlmta1_p)) deallocate(nlmta1_p)
    !if(allocated(nlmta2_p)) deallocate(nlmta2_p)
    !if(allocated(fqwei_p)) deallocate(fqwei_p)

    if(allocated(ipaw)) deallocate(ipaw)
    if(allocated(chgcr)) deallocate(chgcr)
    if(allocated(iltpw)) deallocate(iltpw)

    if(allocated(ilmt_phi)) deallocate(ilmt_phi)
    if(allocated(ilmt_pao)) deallocate(ilmt_pao)

    if(allocated(lmta_phi)) deallocate(lmta_phi)
    if(allocated(q_phi)) deallocate(q_phi)

    if(allocated(q_phirt_pw)) deallocate(q_phirt_pw)

    if ( allocated(radr_paw) ) deallocate( radr_paw )
    if ( allocated(psirpw) ) deallocate(psirpw)
    if ( allocated(phirpw) ) deallocate(phirpw)
    if ( allocated(wf_nrc) ) deallocate(wf_nrc)
    if ( allocated(wf_mnrc) ) deallocate(wf_mnrc)
    if ( allocated(rhcorpw) ) deallocate(rhcorpw)

    if(allocated(lmta_pao)) deallocate(lmta_pao)
    if(allocated(ibpao)) deallocate(ibpao)
    if(allocated(rhvg_l)) deallocate(rhvg_l)

    if(allocated(qitg_BP)) deallocate(qitg_BP)

! =========================== modified by K. Tagami ================= 5.0
!    if(flg_paw) call m_PP_dealloc_paw()

    if(flg_paw) then
      call m_PP_dealloc_paw()
!!$    else if ( sw_mix_charge_hardpart ) then
    else if ( sw_mix_charge_hardpart == ON) then
      call m_PP_dealloc_aug_charge_mixing()
    else if ( af == 1 ) then
       if( associated(ia2ia_symmtry_op) )      deallocate( ia2ia_symmtry_op )
       if( associated(ia2ia_symmtry_op_inv) )  deallocate( ia2ia_symmtry_op_inv )
    endif
! ===================================================================== 5.0

! ====================================== added by K. Tagami ============== 11.0
    if ( allocated(dion_scr_noncl) )    deallocate( dion_scr_noncl )
    if ( allocated(dion0_noncl) )       deallocate( dion0_noncl )

    if ( allocated(q_noncl) )           deallocate( q_noncl )
    if ( allocated(fqwei_noncl) )       deallocate( fqwei_noncl )
    if ( allocated(fqwei_p_noncl) )     deallocate( fqwei_p_noncl )
!
!    if ( allocated(dion_so) )     deallocate( dion_so )
!    if ( allocated(q_so) )        deallocate( q_so )
! ======================================================================== 11.0

! ======= KT_add ========== 13.0AS
    if ( allocated(phirt) ) deallocate(phirt)
    if ( allocated(qorb) ) deallocate( qorb )
    if ( allocated(nrorb) ) deallocate( nrorb )
    if ( allocated(irorb) ) deallocate( irorb )
    if ( allocated(crorb) ) deallocate( crorb )
! ========================== 13.0AS

  end subroutine m_PP_dealloc

!#####################################################################

!=====================================================================
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q ORG_Parallel
! ==================================================================================================
!$$#ifndef PARA3D
  subroutine m_PP_vanderbilt_type_gncpm2(nfp,it,nfout,gr_l,paramset)
!=====================================================================
! This subroutine was modified from m_PP_vanderbilt_type_gncp2,
! in order to read GNCPPM2 type pseudopotential data.
! (AAS, August 2009)
!---------------------------------------------------------------------
    integer,intent(in) :: nfp, it, nfout
    real(kind=DP),intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical,intent(in) :: paramset
    logical            :: vall_is_detected
    integer            :: mesh_t, iprippex
    integer            :: id_sname = -1

    iprippex = ipripp
!!$    if(.not.paramset) iprippex = iprippex - 1
    if(paramset) iprippex = iprippex - 1

    call tstatc0_begin('m_PP_vanderbilt_type_gncpm2 ',id_sname,1)
    if(it == 1) mmt = 0
    if(it == 1) mmpppp = 0
    if(paramset) then
        ipaw(it)=1
        flg_paw=.true.
    end if
    call mpi_barrier(MPI_CommGroup,ierr)
    call rd_Vloc_psWF_Q_chir_qitg_gncpm2 !-(contained here) (potrvd)
    !        -> ps_xctype, etc.
    call make_index_arrays_nlmt2l_m_t         ! (tlmst1)



    if(sw_orb_popu == ON) then
       call make_index_arrays_nlmt_phi2lmt  ! (tlmst1_phi)
    end if

    if(sw_use_add_proj == ON) then
       call make_index_arrays_nlmt_add2lmt     ! (tlmst1_add)
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       call make_index_arrays_nlmt_pao2lmt     ! (tlmst1_pao)
    end if



    call make_index_arrays_nlt2l_t
    call make_ia2ia_symmtry_op_etc
    call cnstrct_of_PiPjPlPm_k_Qijk_etc
    if(.not.flg_symmtry) then
        call cnstrct_of_CijkClmkVVVVijlm_k
    else
        call cnstrct_of_CijkClmnVVVVijlm_kn
    end if
    if(.not.paramset) then
       call cnstrct_of_localWF_Dion_and_q2    ! betavd
       ! -> betar, dion, q
       call where_does_WFrhvr_damp0            ! -> mesh_t
       call coulomb_potential_in_Rspace(5)    ! -> vvv()   (ptclmb)
       call atomic_xc_potential_in_Rspace     ! -> vxc()   (ptxc)
       ! jfp = 0 (non-spin system for atomic charge) y
       call vlocr_plus_hartree_xc2 ! vlocr = vlocr2, vlocr, vvv, vxc
#ifdef _POT_SMOOTHING_
       call smoothing_vlocr()                 ! vlocr smoothing
#endif
       call cnstrct_of_VHij_VHpsij_Kinij
       call cnstrct_of_dion_kin_ion
       call init_of_dion_paw
    end if

! ===============================- added by K. Tagami ================== 11.0
    if ( .not. paramset ) then
       if ( .not. flg_paw ) then
          if ( .not. pot_has_soc(it) ) then
             if ( SpinOrbit_Mode == ZeffApprox ) then
                call  make_SOC_strength_Zeff_nonpaw
             endif
          endif
       endif
    endif
! ====================================================================== 11.0

    call tstatc0_end(id_sname)

!....................................................................
  contains

   !================================================
    subroutine rd_Vloc_psWF_Q_chir_qitg_gncpm2
   !================================================
      integer       :: i
      real(kind=DP) :: dummy

      call read_natomn_ival_iloc_itpcc()
      call read_ps_xctype()
      call read_alp_cc()
      call read_nmesh_xh_rmax()

      if(mype == 0) then
         vall_is_detected = check_vall()
         if(vall_is_detected) then
            if(paramset) then
               read(nfp,*) (dummy,i=1,nmesh(it))
            else
               read(nfp,*) (vvv(i),i=1,nmesh(it)) ! ### VAE[scr](r) ### Now vvv is VAE[scr](r)
               vloc_scr_ae(1:nmesh(it))=vvv(1:nmesh(it))
            end if
         endif
         if(paramset) then
            read(nfp,*) (dummy,i=1,nmesh(it))
            read(nfp,*) (dummy,i=1,nmesh(it))
            read(nfp,*) (dummy,i=1,nmesh(it))
         else
            read(nfp,*) (vlocr(i) ,i=1,nmesh(it)) ! ### Vloc[scr](r) ###
            vloc_scr_ps(1:nmesh(it))=vlocr(1:nmesh(it))
            read(nfp,*) (vlocr2(i),i=1,nmesh(it)) ! ### Vloc[ion](r) ###
            read(nfp,*) (rhvr(i)  ,i=1,nmesh(it)) ! ### 4*pi*r*r*n(r) ###
            call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,h(it)) ! -(b_P.P.) -> radr,h(it)
            radr_paw(:,it)=radr(:)
         end if
      end if
      if(npes > 1 .and. .not.paramset) then
         call mpi_bcast(vvv,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vlocr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vlocr2,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(rhvr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(h(it),1,        mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(radr ,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vall_is_detected,1,mpi_logical,0,MPI_CommGroup,ierr)
         call mpi_bcast(vloc_scr_ae,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(vloc_scr_ps,  nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
!!$      J. Koga 09.06.02
         call mpi_bcast(radr_paw(:,it),nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
      end if

      call determine_lpsmax ! -> (lpsmax) 3=d or 4=f
      if(ipri >= 2) write(nfout,*) ' !lpsmax = ', lpsmax(it)

      call rd_itau_etc_phir_chir_gncpm2 !->(phir,chir):pWF, chir:(e-T-Vloc)*phir

      if(.not.paramset) then
         call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr&
              & ,wos)
      end if
      call rd_qrsps_iqitg_and_qitgft !-(contained in m_PP_vanderbilt_type)
    end subroutine rd_Vloc_psWF_Q_chir_qitg_gncpm2

   !=============================================
    subroutine rd_itau_etc_phir_chir_gncpm2
   !=============================================
      integer       :: il, t1, nrc, i, nm, mord
      real(kind=DP) :: dummy,eps_phi0,eps_phi

      nm = nmesh(it)
      write(6,*) ' !! lpsmax = ', lpsmax(it)
      loop_L : do il = 1, lpsmax(it)
         if(mype == 0) &
              & call read_itau_ivanl(nfp,nfout,iloc(it),il,itau(il,it)&
              &   ,ivanl(il,it), iprippex) ! -(b_P.P.) ->(itau,ivanl)
         if(npes > 1) then
            call mpi_bcast(itau(il,it), 1,mpi_integer,0,MPI_CommGroup,ierr)
            call mpi_bcast(ivanl(il,it),1,mpi_integer,0,MPI_CommGroup,ierr)
         end if
         loop_tau : do t1 = 1, itau(il,it)
            if(ipri >= 2) write(nfout,*) ' t1 = ', t1
            if(mype == 0) then
                call read_tau_eps_nrc_mord(nfp,nfout,kord &
                 &                           ,il,t1,eps(il,t1,it),nrc,mord) ! ->(eps,t1,it)
                if(.not.paramset) wf_nrc(il,t1,it)=nrc
            end if
!            if(npes > 1) then
!                call mpi_bcast(eps(il,t1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
!                if(.not.paramset) &
!                call mpi_bcast(wf_nrc(il,t1,it),1,mpi_integer,0,MPI_CommGroup,ierr)
!            end if

            if(paramset .and. mype == 0) then
               if(nrc == 0) then
                  read(nfp,*) (dummy,i=1,nm);
                  read(nfp,*) (dummy,i=1,nm);      read(nfp,*) (dummy,i=1,nm)
               else if(nrc > 0) then
                  read(nfp,*) (dummy,i=1,nm);
                  read(nfp,*) (dummy,i=0,mord);    read(nfp,*) (dummy,i=nrc+1,nm)
               end if
            else if(.not.paramset) then
               if(mype == 0) then
                  if(nrc == 0) then
                     read(nfp,*) (psir(i,il,t1),i=1,nm) ! psirpw = r*psi[n](r)
                     read(nfp,*) (phir(i,il,t1),i=1,nm) ! phir = r*phi[n](r)
                     read(nfp,*) (chir(i,il,t1),i=1,nm) ! chir = Vion[n](r)
                     chir(1:nm,il,t1) = (chir(1:nm,il,t1)-vlocr2(1:nm)) * phir(1:nm,il,t1)
                        ! chir = Vion[n](r), vloc2 = Vloc[ion](r) --> chir = r*chi[n](r)
                     eps_phi0=1.d0
                     do i=nm,1,-1
!print *,i,psir(i,il,t1),phir(i,il,t1),psir(i,il,t1).eq.phir(i,il,t1),abs((psir(i,il,t1)-phir(i,il,t1))/psir(i,il,t1))
                        if(psir(i,il,t1)-phir(i,il,t1).eq.0.d0) then
                            eps_phi=1.d0
                        else
                            eps_phi=abs((psir(i,il,t1)-phir(i,il,t1))/psir(i,il,t1))
                        end if
!                        if(abs(psir(i,il,t1)-phir(i,il,t1)) .gt. eps_phi) then
                        if(eps_phi.gt.5.d0*eps_phi0) then
                            wf_nrc(il,t1,it)=i+1
                            exit
                        end if
                        eps_phi0=eps_phi
                     end do
                  else if(nrc > 0) then
                     read(nfp,*) (psir(i,il,t1),i=1,nm) ! psirpw = r*psi[n](r)
                     read(nfp,*) (copsc(i),i=0,mord)
                     read(nfp,*) (phir(i,il,t1),i=nrc+1,nm) ! phir = r*phi[n](r)
                     call cnvrtp(nm,1,nrc,il,mord,copsc,radr,phir(1,il,t1))   !-(b_P.P._f77.f)
                     call cnvrtc(nm,1,nrc,il-1,mord,copsc,radr,eps(il,t1,it)& !-(b_P.P._f77.f)
                          & ,phir(1,il,t1),vlocr,chir(1,il,t1))               ! vlocr = Vloc[scr](r)
                     if(vall_is_detected) then
                        chir(nrc+1:nm,il,t1) = (vvv(nrc+1:nm) &        ! vvv = VAE[scr](r)
                             & - vlocr(nrc+1:nm))*phir(nrc+1:nm,il,t1) ! vlocr = Vloc[scr](r), phir = r*phi[n](r)
                                                                       ! --> chir = r*chi[n](r)
                     else
                        chir(nrc+1:nm,il,t1) = 0.d0
                     end if
                  end if
!                  phirpw(1:nm,il,t1,it)=phir(1:nm,il,t1)
               end if
            end if

            if(npes > 1) then
                call mpi_bcast(eps(il,t1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
                if(.not.paramset) &
                call mpi_bcast(wf_nrc(il,t1,it),1,mpi_integer,0,MPI_CommGroup,ierr)
            end if

         enddo loop_tau
      enddo loop_L
      if(.not.paramset .and. npes > 1) then
         call mpi_bcast(phir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(chir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(psir,mmesh*nloc*ntau,mpi_double_precision,0,MPI_CommGroup,ierr)
      end if
      if(.not.paramset) then
         phirpw(1:mmesh,1:nloc,1:ntau,it)=phir(1:mmesh,1:nloc,1:ntau)
         psirpw(1:mmesh,1:nloc,1:ntau,it)=psir(1:mmesh,1:nloc,1:ntau)
         wf_mnrc(it)=maxval(wf_nrc(:,:,it))
      end if
      if(paramset) then
!         if(it == 1) ntau = itau(1,it)
!         do il = 2, lpsmax(it)
!            if(ntau < itau(il,it)) ntau = itau(il,it)
!         end do
         ntau = maxval(itau)
      end if

    end subroutine rd_itau_etc_phir_chir_gncpm2

   !=================================================
    subroutine rd_qrsps_iqitg_and_qitgft
   !=================================================
      integer :: il1, il2, tau1, tau2, il3, tmin, l3s, l3l,  mm, i
      integer :: nrc, mord
      real(kind=DP) :: dummy
      character(len=32):: cdummy
!real(8):: fn(mmesh),dfn(mmesh),ddfn(mmesh)
      mm = 0
      if(.not.paramset) then
         qij = 0.d0; qvij = 0.d0
      end if
      Loop_L1 : do il1 = 1, lpsmax(it)
         if(iloc(it) == il1) cycle
         Loop_tau1 : do tau1 = 1, itau(il1,it)
            Loop_L2 : do il2 = il1, lpsmax(it)
               if((iloc(it)==il2) .or. (ivanl(il1,it)==0.and.ivanl(il2,it)==0)) cycle
               tmin = 1
               if(il1 == il2) tmin = tau1
               Loop_tau2 : do tau2 = tmin, itau(il2,it)
                  l3s = abs(il1-il2)
                  l3l = abs(il1+il2) - 2
                  Loop_L3 : do il3 = l3s+1, l3l+1, 2
                    if(iprippex>=2) write(nfout,'(" il3 = ",i5)') il3
                     if(il3-1 > lcmax ) cycle
                     mm = mm + 1
                     if(.not.paramset) iqitg(il1,tau1,il2,tau2,il3,it) = mm + mmt
                     if(mype == 0) call read_nrc_mord(nfp,nfout, kord,&
                          &            mm,il1,tau1,il2,tau2,il3,nrc,mord, iprippex) ! ->(nrc,mord)
                     if(paramset .and. mype == 0) then
                        if(nrc > 0 ) then
                           read(nfp,*) (dummy,i=0,mord)
                           read(nfp,*) (dummy,i=nrc+1,nmesh(it))
                        else if(nrc == 0) then
                           read(nfp,*) (dummy,i=1,nmesh(it))
                        end if
                     else if(.not.paramset) then
                        if(mype == 0) then
                           if(nrc > 0 ) then
                              read(nfp,*) (copsc(i),i=0,mord)
                              read(nfp,*) (qrsps(i),i=nrc+1,nmesh(it))
                              call cnvrtp(nmesh(it),1,nrc,il3+1,mord,copsc,radr,qrsps)
                           else if(nrc == 0 ) then
                              read(nfp,*) (qrsps(i),i=1,nmesh(it))
                           endif
                        end if
                        if(npes > 1) &
                             & call mpi_bcast(qrsps,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)

                        if(.not.paramset) qrspspw(1:mmesh,mm+mmt)=qrsps(1:mmesh)

                        if(istress == 0) then
                           call qitgft(ista_kngp,iend_kngp,mmesh,nqitg,nmesh(it),ipri &
                             &,nfout,il3-1,radr,wos,qrsps,gr_l &
                             &,univol,mm+mmt,qitg_l,wkx,wky)
                        else if(istress == 1) then
                           call qitgft_diff(ista_kngp,iend_kngp,mmesh,nqitg,nmesh(it),ipri &
                             &,nfout,il3-1,radr,wos,qrsps,gr_l &
                             &,univol,mm+mmt,qitg_l,qitg_diff_l &
                             &,wkx,wky,wkz1,wkz2)
                        endif
                        if(il3 == 1 .and. il1 == il2) &
                             & call qij_qvij_from_qrsps_etc0(tau1,tau2,il1)
                     end if
                  end do Loop_L3
               end do Loop_tau2
            end do Loop_L2
         end do Loop_tau1
      end do Loop_L1
      mmt = mmt + mm
      if(paramset) nqitg = mmt
      if(.not.paramset) write(6,'(" nqitg = ",i5)') nqitg !K.Mae print->write
      if(itpcc(it) == 1) then
         if(mype == 0) then
            call read_chgpc_nrc_mord(iprippex,nfp,nfout,chgpc(it),nrc,mord)
            ! --- > rhpcr,  copsc is a temporary uesd parameter
            if(nrc > 0) then
               if(paramset) then
                  read(nfp,*) (dummy,i=0,mord)
                  read(nfp,*) (dummy,i=nrc+1,nmesh(it))
               else
                  read(nfp,*) (copsc(i),i=0,mord)
                  read(nfp,*) (rhpcr(i),i=nrc+1,nmesh(it))
                  call cnvrtp(mmesh,1,nrc,2,mord,copsc,radr,rhpcr)
               end if
            else if(nrc == 0) then
               if(paramset) then
                  read(nfp,*) (dummy, i=1,nmesh(it))
               else
                  read(nfp,*) (rhpcr(i),i=1,nmesh(it))
               end if
            endif
            ! <---- rhpcr
            rhpcrpw(:,it)=rhpcr(:)
         end if
         if(npes > 1) then
            call mpi_bcast(chgpc(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
            if(.not.paramset) then
                call mpi_bcast(rhpcr,nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
                rhpcrpw(:,it)=rhpcr(:)
            end if
         end if
      endif

     if(mype == 0) then
        read(nfp,*) cdummy
        read(nfp,*) chgcr(it)
        write(nfout,400) chgcr(it)
400 format(' ',' chgcr = ',f12.8)
    ! --- > rhpcr,  copsc is a temporary uesd parameter
        if(paramset) then
          read(nfp,*) (dummy, i=1,nmesh(it))
        else
          read(nfp,*) (rhcorpw(i,it),i=1,nmesh(it))
        end if
        ! <---- rhpcr
     end if
     if(npes > 1) then
        call mpi_bcast(chgcr(it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
        if(.not.paramset) call mpi_bcast(rhcorpw(:,it),nmesh(it),mpi_double_precision,0,MPI_CommGroup,ierr)
     end if

     if(itpcc(it) == 0) then
        if(.not.paramset) rhpcrpw(:,it)=0.d0  !rhcorpw(:,it) !0.d0      !2009.4.18
     end if
!if(.not.paramset) then
!do i=1,nmesh(it)
!print '(3e19.6)',radr(i),rhpcrpw(i,it),rhcorpw(i,it)
!end do
!stop
!end if
     ! some test
!if(.not.paramset) then
!fn(1:mmesh)=radr(1:mmesh)*vlocr2(1:mmesh)
!call calc_ddiff_exp(i,5,mmesh,radr,fn,dfn,ddfn)
!do i=2,mmesh-1
!print '(6f19.6)',radr(i),vlocr(i),vlocr2(i),radr(i)*ddfn(i),rhpcrpw(i,it),rhcorpw(i,it)
!end do
!stop
!end if

    end subroutine rd_qrsps_iqitg_and_qitgft

   !========================================
    subroutine make_index_arrays_nlmt2l_m_t
   !========================================
      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2
      mm = 0
      do il1 = 1, lpsmax(it)
         if(il1 == iloc(it)) cycle
         do tau1 = 1, itau(il1,it)
            do im1 = 1, il1*2-1
               mm = mm + 1
               ltp( mm,it) = il1
               mtp( mm,it) = im1
               taup(mm,it) = tau1
            enddo
         enddo
      enddo
      ilmt(it) = mm
      if(ipri >= 2) then
         if(iprippex>=1)write(nfout,*) ' ILMT, LTP, MTP, TAUP  '
         if(iprippex>=1)write(nfout,400) (n,ltp(n,it),mtp(n,it),taup(n,it),n=1,ilmt(it))
400      format((' ',4i5))
      end if

      mm = 0
      do ilmt1 = 1, ilmt(it)
         mm = mm + 1
         if(.not.paramset) then
            index_lmt1_lmt2(mm,it,1) = ilmt1
            index_lmt1_lmt2(mm,it,2) = ilmt1
            w_non0_lmtxlmt(mm,it)    = 1
         endif
      enddo

      do ilmt1 = 1, ilmt(it)
         il1 =  ltp(ilmt1,it)
         im1 =  mtp(ilmt1,it)
         do ilmt2 = ilmt1+1, ilmt(it)
            il2 = ltp(ilmt2,it)
            im2 = mtp(ilmt2,it)
            if(il1 == il2 .and. im1 == im2) then
               mm = mm + 1
               if(.not.paramset) then
                  index_lmt1_lmt2(mm,it,1) = ilmt1
                  index_lmt1_lmt2(mm,it,2) = ilmt2
                  w_non0_lmtxlmt(mm,it)    = 2
               endif
            endif
         enddo
      enddo
      n_non0_lmtxlmt(it) = mm
      if(iprippex>=2)write(nfout,*) ' ! #non0_elements = ',mm
    end subroutine make_index_arrays_nlmt2l_m_t

   !========================================
    subroutine make_index_arrays_nlmt_add2lmt
   !========================================
      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2
      mm = 0
      do il1 = 1, lpsmax(it)
         if(il1 /= iloc(it)) cycle
         do im1 = 1, il1*2-1
            mm = mm + 1
            ltp_add( mm,it) = il1
            mtp_add( mm,it) = im1
         enddo
      enddo
      ilmt_add(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP TAUP")')
         write(nfout,400) (n,ltp_add(n,it),mtp_add(n,it),1,n=1,ilmt_add(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_add2lmt

   !==========================================
    subroutine make_index_arrays_nlmt_phi2lmt
   !==========================================
      integer mm, il1, tau1, im1, n, il2, im2, ilmt1, ilmt2, i
      logical l_exists, t_exists
      integer :: ip_from_l(lpsmax(it)),ip,jp,ig

      do i=1,num_projectors
         if(it/=proj_attribute(i)%ityp) cycle
         ig = proj_attribute(i)%group
         do ip=1,num_proj_elems(ig)
            jp = proj_group(ip,ig)
            il1=proj_attribute(jp)%l+1
            ip_from_l(il1)=ip
         end do
      end do

      mm = 0
      L_loop: do il1 = 1, lpsmax(it)
         T_loop: do tau1 = 1, itau(il1,it)
            !!$do i=1,norbital(it)
            !!$   if(il1 == l_orb(i,it)+1 .and. tau1 == t_orb(i,it)) then
            !!$      do im1 = 1, il1*2-1
            !!$         mm = mm + 1
            !!$         ltp_phi( mm,it) = il1
            !!$         mtp_phi( mm,it) = im1
            !!$         taup_phi(mm,it) = tau1
            !!$      enddo
            !!$   end if
            !!$enddo
            do i=1,num_projectors
               if(it /= proj_attribute(i)%ityp) cycle
               if(il1 == proj_attribute(i)%l+1 .and. &
                & tau1 == proj_attribute(i)%t) then
                  do im1 = 1, il1*2-1
                     mm = mm + 1
                     ltp_phi( mm,it) = il1
                     mtp_phi( mm,it) = im1
                     taup_phi(mm,it) = tau1
                     iproj_phi(mm,it) = ip_from_l(il1)
                  end do
               end if
            end do
         enddo T_loop
      enddo L_loop
      ilmt_phi(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_phi(n,it),mtp_phi(n,it),taup_phi(n,it),n=1,ilmt_phi(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_phi2lmt

   !========================================
    subroutine make_index_arrays_nlmt_pao2lmt
   !========================================
      integer mm, il1, tau1, im1, n, ilmt1
      mm = 0
      do il1 = 1, lpsmax(it)
         do tau1 = 1, itau(il1,it)
            if(irank(il1,tau1,it)<=0) cycle
            do im1 = 1, il1*2-1
               mm = mm + 1
               ltp_pao( mm,it) = il1
               mtp_pao( mm,it) = im1
               taup_pao(mm,it) = tau1
            enddo
         enddo
      enddo
      ilmt_pao(it) = mm
      if(iprippex >= 2) then
         write(nfout,'(" !PP ILMT, LTP, MTP, TAUP")')
         write(nfout,400) (n,ltp_pao(n,it),mtp_pao(n,it),taup_pao(n,it),n=1,ilmt_pao(it))
400      format((' !PP ',4i5))
      end if

    end subroutine make_index_arrays_nlmt_pao2lmt
   !==========================================
    subroutine cnstrct_of_localWF_Dion_and_q2
   !==========================================
      integer       :: il1,t1,t2,i,lmt1,il11,lmt2,im1,il22,im2,ip
      real(kind=DP) :: dion_tmp

      if(iprippex >= 2) then
      write(nfout,'(" -- cnstrct_of_localWF_Dion_and_q --")')
      write(nfout,'(" !! it = ",i3)') it
      write(nfout,'(" !! lpsmax ( ",i3,") = ", i3)') it,lpsmax(it)
      write(nfout,'(" !! iloc   ( ",i3,") = ", i3)') it,iloc(it)
      write(nfout,'(" !! nmesh  ( ",i3,") = ", i5)') it,nmesh(it)
      write(nfout,'(" !! wos(100:101) = ",2d16.8)') wos(100),wos(101)
      endif

      dion(:,:,it) = 0.d0; q(:,:,it) = 0.d0
      loop_il1 :do il1 = 1, lpsmax(it)
         if(iprippex >= 2) then
         write(nfout,'(" !! il1 = ",i5)') il1
         write(nfout,'(" !! itau(",i3,",",i3,") = ",i3)') il1,it,itau(il1,it)
         endif
         if(iloc(it) == il1) cycle
         do t1 = 1, itau(il1,it)
            if(iprippex >= 2) then
            write(nfout,'(" !! phir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, phir(100,il1,t1),phir(101,il1,t1)
            write(nfout,'(" !! chir(100:101,",i3,",",i3,") = ",2d16.8)') il1,t1, chir(100,il1,t1),chir(101,il1,t1)
            endif
            do t2 = 1,itau(il1,it)
               s(t1,t2) = 0.d0
               do i = 1, nmesh(it)
                  s(t1,t2) = s(t1,t2) + wos(i)*phir(i,il1,t1)*chir(i,il1,t2)
               enddo
            enddo
         enddo
         if(iprippex>=2)then
         write(nfout,*)
         do t1 = 1, itau(il1,it)
            do t2 = 1, itau(il1,it)
               write(nfout,'(" ! B[nm] il=",i3," : s (",i3,",",i3,") = ",d16.8)') il1,t1,t2,s(t1,t2)
            end do
         end do
         write(nfout,*)
         do t1 = 1, itau(il1,it)
            do t2 = 1, itau(il1,it)
               write(nfout,'(" ! Q[nm] il=",i3," : q (",i3,",",i3,") = ",d16.8)') il1,t1,t2,qij(t1,t2,il1)
            end do
         end do
         write(nfout,*)
         do t1 = 1, itau(il1,it)
            do t2 = 1, itau(il1,it)
               write(nfout,'(" ! QV[nm] il=",i3," : qv (",i3,",",i3,") = ",d16.8)') il1,t1,t2,qvij(t1,t2,il1)
            end do
         end do

         write(nfout,'(" -- after t1 loop")')
         endif
         do lmt1 = 1, ilmt(it)
            il11 = ltp(lmt1,it)
            im1  = mtp(lmt1,it)
            t1  = taup(lmt1,it)
            do lmt2 = 1, ilmt(it)
               il22 =  ltp(lmt2,it)
               im2  =  mtp(lmt2,it)
               t2   = taup(lmt2,it)
               if(il11==il1 .and. il22==il1 .and. im1==im2) then
                  q   (lmt1,lmt2,it) = qij(t1,t2,il1)
                  dion(lmt1,lmt2,it) = s(t1,t2) - qvij(t1,t2,il1) &
                       & + eps(il1,t2,it)*qij(t1,t2,il1)
               endif
            enddo
         enddo
         if(iprippex>=2) write(nfout,'(" -- after lmt1 loop")')

         do ip = 1, n_non0_lmtxlmt(it)
            lmt1 = index_lmt1_lmt2(ip,it,1)
            il11 =  ltp(lmt1,it)
            if(il11 == il1) then
               im1  =  mtp(lmt1,it)
               t1   = taup(lmt1,it)
               lmt2 = index_lmt1_lmt2(ip,it,2)
               t2  =  taup(lmt2,it)
               q_indp(ip,it) = qij(t1,t2,il1)
               dion_indp(ip,it) = s(t1,t2) - qvij(t1,t2,il1) &
                    & + eps(il1,t2,it)*qij(t1,t2,il1)
            endif
         enddo

         if(itau(il1,it).ne.0) &
         call matrix_inversion(nfout,ntau,itau(il1,it),s,sinv)  ! s(i,j) = <phir(i)|chir(j)>

         do t1 = 1, itau(il1,it)
            betar(1:nmesh(it),il1,t1,it) = 0.d0
            do t2 = 1, itau(il1,it)
               do i= 1, nmesh(it)
                  betar(i,il1,t1,it) = betar(i,il1,t1,it) &
                       & + sinv(t2,t1)*chir(i,il1,t2)
               enddo
            enddo
         end do
      end do loop_il1

      do lmt2 = 2, ilmt(it)
         do lmt1 = 1, lmt2-1
            dion_tmp = (dion(lmt1,lmt2,it) + dion(lmt2,lmt1,it))*0.5d0
            dion(lmt1,lmt2,it) = dion_tmp
            dion(lmt2,lmt1,it) = dion_tmp
         enddo
      enddo

      if(iprippex>=2)then
      write(nfout,*)
      write(nfout,*) ' -- dion ---'
      do lmt1 = 1, ilmt(it)
         write(nfout,'(i3,9f8.5/9f8.5)') lmt1 &
              &               ,(dion(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
      enddo

      write(nfout,*) ' --  q   ---'
      do lmt1 = 1, ilmt(it)
         write(nfout,'(i3,9f8.5/9f8.5)') lmt1 &
              &               ,(q(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
      enddo
      write(nfout,*)
      write(nfout,'(10("(",2i3,")"))') &
           &(index_lmt1_lmt2(ip,it,1),index_lmt1_lmt2(ip,it,1),&
           & ip=1,n_non0_lmtxlmt(it))
      write(nfout,*) ' -- dion_indp --'
      write(nfout,'(10f8.5)') (dion_indp(ip,it),ip = 1, n_non0_lmtxlmt(it))
      write(nfout,*) ' -- q_indp --'
      write(nfout,'(10f8.5)') (q_indp(ip,it),ip = 1, n_non0_lmtxlmt(it))
      write(nfout,*) ' -- weight --'
      write(nfout,'(10i8)') (w_non0_lmtxlmt(ip,it),ip=1,n_non0_lmtxlmt(it))
      endif
    end subroutine cnstrct_of_localWF_Dion_and_q2

   !==================================
    subroutine where_does_WFrhvr_damp0
   !==================================
      integer     :: i
      real(kind=DP), parameter :: CRDAMP = 1.d0
      real(kind=DP), parameter :: CRDIST = 10.d0
      do i = 10, nmesh(it)-1
         if(rhvr(i)-rhvr(i+1) > CRDAMP .and. radr(i) < CRDIST) then
            mesh_t = i
            if(iprippex>=2) write(nfout,'(" LMTO pot. r_ws=",i5,f12.6)') i, radr(i)
            return
         end if
      enddo
      mesh_t = nmesh(it)
    end subroutine where_does_WFrhvr_damp0

   !==============================================
    subroutine coulomb_potential_in_Rspace(nsize)
   !==============================================
      integer, intent(in) :: nsize
      real(kind=DP),allocatable,dimension(:) :: da, db ! d(nsize)
      real(kind=DP)       :: s2, rhs, rhs1, bm
      integer             :: i
     !+++++++++++++++++++++++++++++++
      allocate(da(nsize)); da = 0.d0
      allocate(db(nsize)); db = 0.d0
     !+++++++++++++++++++++++++++++++
      s2 = dlog(rhvr(2)/rhvr(1))/h(it)
      rhs = rhvr(1)
      wkx(1)  = rhs*radr(1)/(s2+1)
      wky(1)  = rhs/s2
      db(1)   = h(it)*rhs*3.d0
      da(1)   = db(1)*radr(1)
      rhs1    = rhs
      do i = 2,3
         rhs     = rhvr(i)
         wkx(i)  = wkx(i-1)+h(it)*(rhs*radr(i)+rhs1*radr(i-1))*0.5d0
         wky(i)  = wky(i-1)+h(it)*(rhs        +rhs1          )*0.5d0
         db(i)   = h(it) *rhs*3.d0
         da(i)   = db(i)*radr(i)
         rhs1    = rhs
      enddo
      do i = 4,mesh_t
         rhs    = rhvr(i)
         db(4)  = h(it)*rhs*3.d0
         da(4)  = db(4)*radr(i)
         wkx(i)=(9*wkx(i-1)-wkx(i-3)+da(4)+2.d0*da(3)-da(2))/8.d0
         wky(i)=(9*wky(i-1)-wky(i-3)+db(4)+2.d0*db(3)-db(2))/8.d0
         da(1)  = da(2)
         db(1)  = db(2)
         da(2)  = da(3)
         db(2)  = db(3)
         da(3)  = da(4)
         db(3)  = db(4)
      enddo
      bm               = wky(mesh_t)
     !C--*--COULOMB POTENTIAL RVC
      vvv = 0.d0
      do i = 1,mesh_t
         vvv(i) = wkx(i)+radr(i)*(bm-wky(i))
      enddo
     !+++++++++++++++++++++++++++++++
      deallocate(da); deallocate(db)
     !+++++++++++++++++++++++++++++++
    end subroutine coulomb_potential_in_Rspace

   !=========================================
    subroutine atomic_xc_potential_in_Rspace   ! (ptxc)
   !=========================================
      real(kind=DP) :: chgrsq, vxc1, vxc2, rs, zet, yc
      integer       :: jsm = 1  ! 1: non-mag, 2: magnetic
      real(kind=DP), parameter :: delta40 = 1.d-40
      real(kind=DP) :: exc
      integer       :: i, ncut
#ifdef GGA_CHECK
      if(iprippex>=2) write(nfout,'(" -- atomic_xc_potential_in_Rspace -- ")')
#endif
      wkx = 0.d0
      if(itpcc(it) == 0) then
         wkx(1:mesh_t) = rhvr(1:mesh_t)
      else if(itpcc(it) == 1) then
         wkx(1:mesh_t) = rhvr(1:mesh_t) + rhpcr(1:mesh_t)
      endif
      if(what_is_ps_xctype() == GGA) then ! -(m_PseudoPotential)
            call gga_grad_alloc()
            call get_gradient_of_rho(jsm,mddiff,mmesh,mesh_t,wkx,radr&
                 &,h(it),fdiff,coeff1,coeff2,rho,grdnts) ! -(b_PP)  -> rho, grdnts
#ifdef GGA_CHECK
            write(6,'(" after get_gradient_of_rho ")')
            write(6,'(" rad, rho, grdnts(*,1), grdnts(*,2)")')
            do i = 1, mesh_t
               write(6,'(" ( ",i4,") rad, rho, grad, grad2 = ",4d20.8," (get_gradient)")') &
                    &        i, radr(i), rho(i), grdnts(i,1), grdnts(i,2)
            end do
#endif

            if(         ps_xctype == 'ldapw91' .or. ps_xctype == 'LDAPW91' &
                 & .or. ps_xctype == 'ldapbek' .or. ps_xctype == 'LDAPBEK' ) grdnts = 0.d0

            call gga_xc_potential_atomic(jsm,len_ps_xctype,ps_xctype &
                 &,mmesh,mesh_t,rho,radr,grdnts,fdiff,vxc,eps_chg,eps_grad)
            call gga_grad_dealloc()
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
            do i = mmesh,1,-1
               if(wkx(i)/(PAI4*radr(i)**2) > eps_chg) then
                  ncut = i
                  goto 1001
               end if
            end do
1001        continue
            do i = ncut+1, mmesh
               vxc(i) = 0.d0
            end do
#endif
!         endif

#ifdef GGA_CHECK
         write(6,'(" ps_xctype, xctype = ",a7,2x,a7)') ps_xctype, xctype
         do i = 1, nmesh(it)
            write(6,'(" vxc(",i4,") = ",d18.10, " radr(",i4,") = ",d18.10)') i, vxc(i),i, radr(i)
         end do
#endif
      else
         vxc = 0.d0
         do i = 1, mesh_t
            vxc1 = 0.d0; vxc2 = 0.d0
            chgrsq = wkx(i)
            if(dabs(chgrsq) > delta40) then
               rs = (3*radr(i)*radr(i)/chgrsq)**(1.d0/3.d0)
               zet = 0.d0
          !--*--EXCHANGE-CORRELATION POTENTIALS ARE GIVEN IN RYD. UNITS.
               call xcpot(ps_xctype,vxc1,vxc2,rs,zet,yc)
            endif
            vxc(i) = radr(i)*vxc1*0.5d0
         enddo
      end if
    end subroutine atomic_xc_potential_in_Rspace

   !==================================
    subroutine vlocr_plus_hartree_xc2
   !==================================
      integer :: meshup,i

#ifdef GGA_CHECK
      if(iprippex>=2) write(nfout,*) ' << vlocr_plus_hartree_xc2 >>'
#endif

      if(nmesh(it) > mesh_t) then
         meshup = mesh_t
      else
         meshup = nmesh(it)
      endif
#ifdef _DEBUG_WRITE_
      write(nfout,'(" --- vlocr --- ")')
      write(nfout,'(3d26.18)') (vlocr(i),i=1,nmesh(it))
      write(nfout,'(" ! -- vlocr before adding vcoulomb and xc potential --")')
#endif
      call sum_of_vlocr(it,nfout)
      call sum_of_vvv_and_vxc(it,nfout) ! Now vvv is the Hartree potential

      !do i = 1, meshup
      !   vlocr(i) = vlocr(i) - (vvv(i)+vxc(i))/radr(i)
      !enddo
      do i = 1, nmesh(it)
         vlocr(i) = vlocr2(i)
      end do

#ifdef GGA_CHECK
      do i = meshup, meshup-100, -1
         if(iprippex>=2) write(nfout,'(" Vlocr, a3, a2: ",i6,3f12.6)') i, vlocr(i),vvv(i),vxc(i)
      end do
#endif

#ifdef _DEBUG_WRITE_
      write(nfout,*) ' --- after - hartree - xc --'
      write(nfout,'(3d30.20)') (vlocr(i),i=1,nmesh(it))
      write(nfout,*) ' --- coulomb ---------------'
      write(nfout,'(6d15.7)') (vvv(i),i=1,nmesh(it))
      write(nfout,*) ' --- xc      ---------------'
      write(nfout,'(6d15.7)') (vxc(i),i=1,nmesh(it))
      write(nfout,*) ' --- rhvr    ---------------'
      write(nfout,'(6d15.7)') (rhvr(i),i=1,nmesh(it))
#endif
    end subroutine vlocr_plus_hartree_xc2

   !=============================================
    subroutine qij_qvij_from_qrsps_etc0(t1,t2,il)
   !=============================================
      integer, intent(in) :: t1,t2,il
      integer :: i
      real(kind=DP) :: sij, svij
      sij = 0.d0
      svij = 0.d0
      do i = 1, nmesh(it)
         sij  = sij + wos(i)*qrsps(i)
         svij = svij + wos(i)*qrsps(i)*vlocr(i)
      enddo
      qij( t1,t2,il) = sij
      qvij(t1,t2,il) = svij
      if(t1 /= t2) then
         qij( t2,t1,il) = sij
         qvij(t2,t1,il) = svij
      endif
    end subroutine qij_qvij_from_qrsps_etc0

   !=======================================
    subroutine read_natomn_ival_iloc_itpcc
   !=======================================
      integer, parameter  ::    len_str = 80
      character(len=len_str) :: str
      integer :: natomn
      logical :: comment_statement

      if(mype == 0) then
         comment_statement = .true.
         rewind(nfp)
         write(nfout,*) ' READING POTENTIAL FILE ', nfp
         do while(comment_statement)
            read(nfp,'(a80)') str
            if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
                 & .or. str(1:1) == '%') then
               write(nfout,'(a80)') str
            else
               comment_statement = .false.
            endif
         enddo
         read(str,*)      natomn,ival(it),iloc(it),itpcc(it)
         write(nfout,110) natomn,ival(it),iloc(it),itpcc(it)
110      FORMAT(' ',I4,f8.4,2I4,'  : NATOMN, IVAL, ILOC, ITPCC ')
         if(iatomn(it) /= natomn) then
            write(nfout,*) ' iatomn.ne.natomn ',iatomn(it),natomn
            call phase_error_with_msg(nfout,'iatomn.ne.natomn',__LINE__,__FILE__)
         end if
      end if
      if(npes > 1 ) then
         call mpi_bcast(ival(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(iloc(it), 1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(itpcc(it),1,mpi_integer,0,MPI_CommGroup,ierr)
      end if
      if(iprippex>=2)then
      write(nfout,'(" !!read_natomn_ival_iloc_itpcc: ival(",i4,") = ",f8.4)') it,ival(it)
      write(nfout,'(" !!read_natomn_ival_iloc_itpcc: iloc(",i4,") = ",i3)') it,iloc(it)
      write(nfout,'(" !!read_natomn_ival_iloc_itpcc: itpcc_(",i3,") = ",i3)') it,itpcc(it)
      endif
    end  subroutine read_natomn_ival_iloc_itpcc

   !============================
    subroutine read_ps_xctype()
   !============================
      integer  ::                     counter = 0
      character(len=len_ps_xctype) :: ps_xctype0

      if(iprippex>=2) write(nfout,'(" !!read_ps_xctype: MPI_CommGroup = ",i8)') MPI_CommGroup

      if(mype == 0) then   ! MPI
         counter = counter + 1
         read(nfp,'(a7)') ps_xctype0
         call nameck(ps_xctype0,len_ps_xctype)
         write(nfout,120) ps_xctype0
120      FORMAT(A7,'   : NAME ')
         if(counter == 1) then
            ps_xctype = ps_xctype0
            if(xctype == "nogiven") call m_CtrlP_set_xctype(ps_xctype,nfout)
         else
            if(ps_xctype /= ps_xctype0) then
               write(nfout,*) ' XC TYPE OF ', it,'TH ATOM IS NOT CORRECT'
               call phase_error_with_msg(nfout,'XC type is not correct',__LINE__,__FILE__)
            endif
            if(xctype /= ps_xctype0) then
               write(nfout,*) &
                    &' !W xc-type of ',it,'th atom is not equal to ' &
                    &,'that described in the input file'
               write(nfout,'(" xctype = ",a7," ps_xctype0 = ",a7)') xctype,ps_xctype0
            endif
         endif
      end if
      call mpi_bcast(ps_xctype,len_ps_xctype,mpi_character,0,MPI_CommGroup,ierr)
      if(iprippex>=2) write(nfout,'(" !!read_ps_xctype: ps_xctype = ",a7)') ps_xctype
    end subroutine read_ps_xctype

   !=======================
    subroutine read_alp_cc
   !=======================
      if(mype == 0) then          ! MPI
         read(nfp,*) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
         write(nfout,130) alp(1,it),alp(2,it),cc(1,it),cc(2,it)
130      FORMAT(' ',4F12.6,'  :   ALP,CC')
      end if
      call mpi_bcast(alp(1,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(alp(2,it),1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(cc(1,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(cc(2,it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
    end subroutine read_alp_cc

   !==============================
    subroutine read_nmesh_xh_rmax
   !==============================
      if(mype == 0) then  ! MPI
         read(nfp,*) nmesh(it),xh(it),rmax(it)
         write(nfout,140) nmesh(it),xh(it),rmax(it)
140      format(' ',i4,2f12.6,'  :   nmesh,  xh, rmax')
      end if
      call mpi_bcast(nmesh(it),1,mpi_integer         ,0,MPI_CommGroup,ierr)
      call mpi_bcast(xh(it),   1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(rmax(it), 1,mpi_double_precision,0,MPI_CommGroup,ierr)
      if(iprippex>=2)then
      write(nfout,'(" !! read_nmesh_xc_rmax: nmesh(it) = ",i8)') nmesh(it)
      write(nfout,'(" !! read_nmesh_xc_rmax: xh(it)    = ",f10.6)') xh(it)
      write(nfout,'(" !! read_nmesh_xc_rmax: rmax(it)  = ",f10.6)') rmax(it)
      endif

    end subroutine read_nmesh_xh_rmax

   !============================
    subroutine determine_lpsmax
   !============================
      character(len=7) :: lpsmaxj
      if(mype == 0) then
         if(iloc(it) == 4) then
            lpsmax(it) = iloc(it)
         else
            lpsmax(it) = 3
         endif
         read(nfp,'(a7)') lpsmaxj
         if(lpsmaxj == 'F-STATE') then
            lpsmax(it) = 4
            write(6,'(" lpsmaxj == F-STATE")')
         else
            backspace nfp
         endif
      end if
      call mpi_bcast(lpsmax(it),1,mpi_integer,0,MPI_CommGroup,ierr)
    end subroutine determine_lpsmax

   !==============================
    logical function check_vall()
   !==============================
      character(len=4) :: vall
      read(nfp,'(a4)') vall
      if(vall == 'VALL') then
         check_vall = .true.
      else
         backspace nfp
         check_vall = .false.
      end if
    end function check_vall

#ifdef _POT_SMOOTHING_
   !=============================
    subroutine smoothing_vlocr()
   !=============================
      real(kind=DP), parameter :: delta7 = 1.d-7, delta10 = 1.d-10
      integer       :: icheck, i
      real(kind=DP) :: s
      icheck = 0
      do i = nmesh(it),1,-1
         if((rhvr(i)/pai4/radr(i)**2) < delta7) then
            if(dabs(vlocr(i)+ival(it)/radr(i)) > delta10 ) then
               if(ipri >= 3) then
                  write(nfout,'("VLOC:",i6,f12.6,3d16.8)') &
                       & i,radr(i),vlocr(i), -ival(it)/radr(i)&
                       &,          vlocr(i) + ival(it)/radr(i)
               end if
               vlocr(i) = -ival(it)/radr(i)
            endif
         else
            icheck = icheck + 1
         endif
         if(icheck > 30) exit
      enddo
      s = 0.d0
      do i = 1, nmesh(it)
         s = s + wos(i)*vlocr(i)*radr(i)**2
      end do
      if(ipri >= 2) write(nfout,*) ' !D (smoothing) s = ', s
    end subroutine smoothing_vlocr
#endif

! =========================== added by K. Tagami ==================== 11.0
    subroutine make_SOC_strength_Zeff_nonpaw

      real(kind=DP) :: HyperFineConst

      integer :: ir, ier
      integer :: nrc, il, it1, it2
      real(kind=DP) :: fac1, fac2, sum

      real(kind=DP), allocatable :: wght(:)

      HyperFineConst = 1.0d0 / InvHyperFineConst
      fac1 = 0.5d0 * HyperFineConst**2

      nrc = nmesh(it)

      allocate( wght(1:nrc) ); wght = 0.0d0
      call set_weight_exp( ier, 1, nrc, radr, wght )

      Do il=1, lpsmax(it)
         if  ( il == iloc(it) ) cycle

         Do it1=1,itau(il,it)
            Do it2=1,itau(il,it)
               if ( ae_wavefunctions_are_detected(it) ) then

                  sum=0.d0
                  if ( SpinOrbit_MassCorrection == 0 ) then
                     Do ir = 1, nrc
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2)
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum

                  else if ( SpinOrbit_MassCorrection == 1 ) then
                     Do ir = 1, nrc
                        fac2 = 1.0d0 + HyperFineConst**2 &
                             &        *dble(iatomn(it)) /radr(ir)
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2) /fac2
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum

                  else if ( SpinOrbit_MassCorrection == 2 ) then
                     Do ir = 1, nrc
                        fac2 = 1.0d0 + HyperFineConst**2 /2.0d0 &
                             &        *dble(iatomn(it)) /radr(ir)
                        sum = sum + wght(ir) *dble(iatomn(it)) /radr(ir)**3 &
                             &               *phir_ae(ir,il,it1) &
                             &               *phir_ae(ir,il,it2) /fac2
                     End do
                     Mat_SOC_strength_nonpaw( il, it1, it2, it ) = fac1 *sum

                  endif

               else

                  write(nfout,*) '! *************  *********** '
                  write(nfout,*) '! ** SOC strength cannot be calculated when AE wfns are not found in pseudopotential.'
                  call phase_error_with_msg(nfout,&
                  'SOC strength cannot be calculated when AE wfns are not found in pseudopotential.'&
                  ,__LINE__,__FILE__)

               endif
            End do
         End do
      End do

      deallocate(wght)

    end subroutine make_SOC_strength_Zeff_nonpaw
! ========================================================================= 11.0

   !================================================
    subroutine cnstrct_of_PiPjPlPm_k_Qijk_etc
   !================================================
        integer:: il1,il2,il3,il4,tau1,tau2,tau3,tau4,lks,lkl,mm
        integer:: t2min,lpmx,t3min,t4min,l4min
        integer:: ilk1,ilk2
!        logical:: a,b
        integer:: ilt1,ilt2,ilt3,ilt4,lt4min
!
!        lpmx=lpsmax(it)
!        mm=0
!        do il1=1,lpmx
!            if(iloc(it) == il1) cycle
!            do tau1=1,itau(il1,it)
!                do il2=il1,lpmx
!                    if(iloc(it) == il2) cycle
!                    t2min=1
!                    if(il1 == il2) t2min=tau1
!                    do tau2=t2min,itau(il2,it)
!
!                        do ilk1=abs(il1-il2),il1+il2-2,2
!                            if(ilk1 > lcmax) cycle
!                            do il3=il1,lpmx
!                                if(iloc(it) == il3) cycle
!                                t3min=1
!                                if(il3 == il1) t3min=tau1
!                                do tau3=t3min,itau(il3,it)
!                                    l4min=il3
!                                    if(il1 == il3 .and. tau1 == tau3) l4min=max(il2,il3)
!                                    do il4=l4min,lpmx
!                                        if(iloc(it) == il4) cycle
!                                        a=(il3.eq.il4)
!                                        b=(il1.eq.il3).and.(tau1.eq.tau3).and.(il2.eq.il4)
!                                        t4min=1
!                                        if(a .and. .not.b) then
!                                            t4min=tau3
!                                        else if(.not.a .and. b) then
!                                            t4min=tau2
!                                        else if(a .and. b) then
!                                            t4min=max(tau3,tau2)
!                                        end if
!                                        do tau4=t4min,itau(il4,it)
!
!                                            do ilk2=abs(il3-il4),il3+il4-2,2
!                                                if(ilk2 > lcmax) cycle
!                                                if(ilk1 == ilk2) then
!
!print '(5i4)',il1,il2,il3,il4,ilk1
!                                                    mm=mm+1
!                                                end if
!                                            end do
!
!                                        end do
!                                    end do
!                                end do
!                            end do
!
!                        end do
!
!                    end do
!                end do
!            end do
!        end do
!
!
!10 continue

        mm=0
        do ilt1=1,iltpw(it)
            il1=lppw(ilt1,it)
            do ilt2=ilt1,iltpw(it)
                il2=lppw(ilt2,it)
                do ilk1=abs(il1-il2),il1+il2-2,2
                    if(ilk1 > lcmax) cycle
                    do ilt3=ilt1,iltpw(it)
                        il3=lppw(ilt3,it)
                        lt4min=ilt3
                        if(ilt1==ilt3) lt4min=max(ilt2,ilt3)
                        do ilt4=lt4min,iltpw(it)
                            il4=lppw(ilt4,it)
                            do ilk2=abs(il3-il4),il3+il4-2,2
                                if(ilk2 > lcmax) cycle
                                if(ilk1 == ilk2) then
    !print '(5i4)',il1,il2,il3,il4,ilk1
                                    mm=mm+1
                                    if(.not.paramset) then
                                        ipppp(ilt1,ilt2,ilt3,ilt4,ilk1+1,it)=&
                                                                    mm+mmpppp
                                    ! call makeppppps(ilt1,ilt2,ilt3,ilt4,ilk1,mm+mpppp)
                                        call make_VPiPjPlPm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,0)
!                                        call make_VPiPjPlPm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,0)
                                        call make_VPiPjPlPm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,1)
!                                        call make_VPiPjPlPm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp,1)
                                        call make_VQijQlm_k &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!                                        call make_VQijQlm_k_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
                                        call make_VQijPlPm_kSym &
                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!                                        call make_VQijPlPm_kSym_dbg &
!                                            (ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp)
!print *,ilt1,ilt2,ilt3,ilt4,ilk1,mm+mmpppp
                                    end if
                                end if
                            end do
                        end do
                    end do
                end do
            end do
        end do

        mmpppp=mmpppp+mm
        if(paramset) npppp=mmpppp
        if(.not.paramset) write(6,'(" npppp = ",i5)')

!if(.not.paramset) then
!print *,'Here'
!    stop
!end if

    end subroutine cnstrct_of_PiPjPlPm_k_Qijk_etc


   !======================================
    subroutine make_index_arrays_nlt2l_t
   !======================================
      integer mm, il1, tau1, n, nn, im1
      mm = 0
      nn = 0
      do il1 = 1, lpsmax(it)
         if(il1 == iloc(it)) cycle
         do tau1 = 1, itau(il1,it)
           mm = mm + 1
           lppw(mm,it) = il1
           tppw(mm,it) = tau1
            if(.not.paramset) then
               do im1=1,il1*2-1
                    nn=nn+1
                    index_lmt2lt(nn,it)=mm
               end do
            end if
         enddo
      enddo
      iltpw(it) = mm
      if(mm.gt.nltpw .and. paramset) nltpw=mm
      if(ipri >= 2) then
         write(nfout,*) ' ILT, LPPW, TPPW  '
         write(nfout,400) (n,lppw(n,it),tppw(n,it),n=1,iltpw(it))
400      format((' ',4i5))
      end if

    end subroutine make_index_arrays_nlt2l_t

   !=============================================================
   subroutine make_VPiPjPlPm_k(ilt1,ilt2,ilt3,ilt4,ilk,mp,ae_or_ps)
   !=============================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp
        integer,intent(in):: ae_or_ps                      ! 0--AE   1--PS

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: pij,plm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        allocate(pij(nnrc))
        allocate(plm(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)* &
                        psirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)* &
                        psirpw(1:nnrc,il4,itau4,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                        phirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                        phirpw(1:nnrc,il4,itau4,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        end if

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*pr3(jr)*pr4(jr)+radr(jr+1)**ilk*pr3(jr+1)*pr4(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*pr3(jr)*pr4(jr)+ &
!                            radr(jr+1)**(-1-ilk)*pr3(jr+1)*pr4(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

!        deallocate(pr1,pr2,pr3,pr4,dsum)
        deallocate(pij,plm,dsum,wght)
        return

   end subroutine make_VPiPjPlPm_k

   !=====================================================================
   subroutine make_VPiPjPlPm_k_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp,ae_or_ps)
   !=====================================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp
        integer,intent(in):: ae_or_ps                      ! 0--AE   1--PS

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: pij,plm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        allocate(pij(nnrc))
        allocate(plm(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)* &
                        psirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc2)=psirpw(1:nnrc2,il3,itau3,it)* &
                         psirpw(1:nnrc2,il4,itau4,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                        phirpw(1:nnrc,il2,itau2,it)
            plm(1:nnrc2)=phirpw(1:nnrc2,il3,itau3,it)* &
                         phirpw(1:nnrc2,il4,itau4,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)
        end if

! ===== Symmetrization ===== !

        select case(ae_or_ps)
        case(0)
!            pr1(1:nnrc)=psirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=psirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=psirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=psirpw(1:nnrc,il3,itau3,it)* &
                        psirpw(1:nnrc,il4,itau4,it)
            plm(1:nnrc2)=psirpw(1:nnrc2,il1,itau1,it)* &
                         psirpw(1:nnrc2,il2,itau2,it)
        case(1)
!            pr1(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)
!            pr2(1:nnrc)=phirpw(1:nnrc,il2,itau2,it)
!            pr3(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)
!            pr4(1:nnrc)=phirpw(1:nnrc,il4,itau4,it)
            pij(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                        phirpw(1:nnrc,il4,itau4,it)
            plm(1:nnrc2)=phirpw(1:nnrc2,il1,itau1,it)* &
                         phirpw(1:nnrc2,il2,itau2,it)
        case default
            write(nfout) 'Error in make_Vijlm_k ! Unknown ae_or_ps :',ae_or_ps
            call phase_error_with_msg(nfout,'Error in make_Vijlm_k ! Unknown ae_or_ps',__LINE__,__FILE__)
        end select

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
            sum1=sum1*radr(ir)**(-1-ilk)*pij(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
!                                pr3(i0+j*is)*pr4(i0+j*is)*wght(i0+j*is)
                                plm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
!                                pr3(jr)*pr4(jr)*wght(jr)
                                plm(jr)*wght(jr)
                end do
            end if
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
            sum2=sum2*(radr(ir)**ilk)*pij(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do
        if(ae_or_ps==0) then
            vaeijlm_k(mp)=vaeijlm_k(mp)+sum0*PAI4/(2*dble(ilk)+1)
            vaeijlm_k(mp)=vaeijlm_k(mp)*0.5d0
        else if(ae_or_ps==1) then
            vpsijlm_k(mp)=vpsijlm_k(mp)+sum0*PAI4/(2*dble(ilk)+1)
            vpsijlm_k(mp)=vpsijlm_k(mp)*0.5d0
        end if


!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*pr3(jr)*pr4(jr)+radr(jr+1)**ilk*pr3(jr+1)*pr4(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*pr1(ir)*pr2(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*pr3(jr)*pr4(jr)+ &
!                            radr(jr+1)**(-1-ilk)*pr3(jr+1)*pr4(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*pr1(ir)*pr2(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

!        deallocate(pr1,pr2,pr3,pr4,dsum)
        deallocate(pij,plm,dsum,wght)
        return

   end subroutine make_VPiPjPlPm_k_dbg


   !=====================================================
   subroutine make_VQijQlm_k(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !=====================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,qlm_k,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij,iqlm

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        iqlm=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if
        if(iqlm.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if

        allocate(qij_k(nnrc))
        allocate(qlm_k(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijqlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qlm_k(jr)+radr(jr+1)**ilk*qlm_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qlm_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qlm_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

        deallocate(qij_k,qlm_k,dsum,wght)
        return

   end subroutine make_VQijQlm_k

   !==========================================================
   subroutine make_VQijQlm_k_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==========================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,qlm_k,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij,iqlm

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        iqlm=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if
        if(iqlm.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijQlm_k'
            vqijqlm_k(mp)=0.d0
            return
        end if

        allocate(qij_k(nnrc2))
        allocate(qlm_k(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        qij_k(1:nnrc2)=qrspspw(1:nnrc2,iqij)
        qlm_k(1:nnrc2)=qrspspw(1:nnrc2,iqlm)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qlm_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qlm_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijqlm_k(mp)=sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qlm_k(jr)+radr(jr+1)**ilk*qlm_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qlm_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qlm_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

        deallocate(qij_k,qlm_k,dsum,wght)
        return

   end subroutine make_VQijQlm_k_dbg


   !==========================================================
   subroutine make_VQijPlPm_kSym(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==========================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,plpm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)

        allocate(qij_k(nnrc))
        allocate(plpm(nnrc))
        allocate(dsum(nnrc))
        allocate(wght(nnrc))

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            vqijplpm_ks(mp)=0.d0
            goto 10
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc)=phirpw(1:nnrc,il3,itau3,it)* &
                    phirpw(1:nnrc,il4,itau4,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=sum0*PAI4/(2*dble(ilk)+1)

10 continue

        iqij=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            goto 20
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc)=phirpw(1:nnrc,il1,itau1,it)* &
                    phirpw(1:nnrc,il2,itau2,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qij_k(jr)+radr(jr+1)**ilk*qij_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qij_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qij_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

20      continue

        deallocate(qij_k,plpm,dsum,wght)
        return

   end subroutine make_VQijPlPm_kSym

   !==============================================================
   subroutine make_VQijPlPm_kSym_dbg(ilt1,ilt2,ilt3,ilt4,ilk,mp)
   !==============================================================
        integer,intent(in):: ilt1,ilt2,ilt3,ilt4,ilk,mp

        integer:: nnrc,nnrc2
        integer:: il1,il2,il3,il4
        integer:: itau1,itau2,itau3,itau4
        real(DP),pointer,dimension(:):: qij_k,plpm,dsum,wght
        real(DP):: sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,ier,is
        integer:: iqij

        il1=lppw(ilt1,it)
        il2=lppw(ilt2,it)
        il3=lppw(ilt3,it)
        il4=lppw(ilt4,it)
        itau1=tppw(ilt1,it)
        itau2=tppw(ilt2,it)
        itau3=tppw(ilt3,it)
        itau4=tppw(ilt4,it)

!        nnrc=max(wf_nrc(il1,itau1,it),wf_nrc(il2,itau2,it), &
!                wf_nrc(il3,itau3,it),wf_nrc(il4,itau4,it))
        nnrc=wf_mnrc(it)
        nnrc2=nmesh(it)

        allocate(qij_k(nnrc2))
        allocate(plpm(nnrc2))
        allocate(dsum(nnrc))
        allocate(wght(nnrc2))

        iqij=iqitg(il1,itau1,il2,itau2,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il1,itau1,il2,itau2,il3 : ',il1,itau1,il2,itau2,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            vqijplpm_ks(mp)=0.d0
            goto 10
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc2)=phirpw(1:nnrc2,il3,itau3,it)* &
                    phirpw(1:nnrc2,il4,itau4,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=sum0*PAI4/(2*dble(ilk)+1)

! ******* Symmetrization *******
        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)


10 continue

        iqij=iqitg(il3,itau3,il4,itau4,ilk+1,it)
        if(iqij.eq.0) then
            !!$print '(a,5i3)','il3,itau3,il4,itau4,il3 : ',il3,itau3,il4,itau4,ilk+1
            !!$print *,'not set in iqitg! make_VQijPlPm_kSym'
            goto 20
        end if

        qij_k(1:nnrc)=qrspspw(1:nnrc,iqij)
!        qlm_k(1:nnrc)=qrspspw(1:nnrc,iqlm)
        plpm(1:nnrc2)=phirpw(1:nnrc2,il1,itau1,it)* &
                    phirpw(1:nnrc2,il2,itau2,it)

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
            if(ir==nnrc) then
                sum2=0.d0
            else if((ir<=nnrc-1).and.(ir>=nnrc-4)) then
                do ii=ir,nnrc-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                qij_k(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc,radr,wght)
                do jr=ir,nnrc
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                qij_k(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

! ***** Symmetrization *****

        dsum=0.d0

        do ir=1,nnrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+(radr(i0+j*is)**ilk)* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+(radr(jr)**ilk)* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum1=sum1*radr(ir)**(-1-ilk)*qij_k(ir)
            if(ir==nnrc2) then
                sum2=0.d0
            else if((ir<=nnrc2-1).and.(ir>=nnrc2-4)) then
                do ii=ir,nnrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(radr(i0+j*is)**(-1-ilk))* &
                                plpm(i0+j*is)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nnrc2,radr,wght)
                do jr=ir,nnrc2
                    sum2=sum2+(radr(jr)**(-1-ilk))* &
                                plpm(jr)*wght(jr)
                end do
            end if
            sum2=sum2*(radr(ir)**ilk)*qij_k(ir)
            dsum(ir)=sum1+sum2
!    print *,radr(ir),dsum(ir)
        end do

        call set_weight_exp(ier,1,nnrc,radr,wght)
        sum0=0.d0
        do ir=1,nnrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vqijplpm_ks(mp)=vqijplpm_ks(mp)+sum0*PAI4/(2*dble(ilk)+1)

        vqijplpm_ks(mp)=vqijplpm_ks(mp)*0.5d0

!        do ir=1,nnrc
!            sum1=0.d0
!            sum2=0.d0
!            do jr=1,ir-1
!                sum1=sum1+(radr(jr)**ilk*qij_k(jr)+radr(jr+1)**ilk*qij_k(jr+1))* &
!                        (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum1=sum1*radr(ir)**(-1-ilk)*plpm(ir)
!            do jr=ir,nnrc-1
!                sum2=sum2+(radr(jr)**(-1-ilk)*qij_k(jr)+ &
!                            radr(jr+1)**(-1-ilk)*qij_k(jr+1))* &
!                            (radr(jr+1)-radr(jr))*0.5d0
!            end do
!            sum2=sum2*(radr(ir)**ilk)*plpm(ir)
!            dsum(ir)=sum1+sum2
!!    print *,radr(ir),dsum(ir)
!        end do
!
!        sum0=0.d0
!        do ir=1,nnrc-1
!            sum0=sum0+(dsum(ir)+dsum(ir+1))*(radr(ir+1)-radr(ir))*0.5d0
!        end do
!        sum0=sum0*PAI4/(2*dble(ilk)+1)

20      continue

        deallocate(qij_k,plpm,dsum,wght)
        return

   end subroutine make_VQijPlPm_kSym_dbg

   !=============================================
    subroutine cnstrct_of_VHij_VHpsij_Kinij
   !=============================================
        integer:: il,it1,it2,mp

        if(iprippex>=2)then
        write(nfout,'("--cnstrct_of_VHij_VHpsij_Kinij--")')
        write(nfout,'("    kin_ae_ps     vionaeij     vionpsij     vionpsqij")')
        endif

        do il=1,lpsmax(it)
            if(iloc(it)==il) cycle
            do it1=1,itau(il,it)
                do it2=it1,itau(il,it)
!                    call make_kin_ae_ps(il,it1,it2)
                    call make_kin_ae_ps2(il,it1,it2)
                    call make_vionaeij(il,it1,it2)
                    call make_vionpsij(il,it1,it2)
                    call make_vionpsqij(il,it1,it2)

                    if(iprippex.gt.2)then
                    write(nfout,'(4f15.7)') kin_ae_psij(il,it1,it2,it), &
                                            vionaeij(il,it1,it2,it), &
                                            vionpsij(il,it1,it2,it), &
                                            vionpsqij(il,it1,it2,it)
                    endif
                end do
            end do
        end do
        return
    end subroutine cnstrct_of_VHij_VHpsij_Kinij


   !====================================
    subroutine cnstrct_of_dion_kin_ion
   !====================================
        integer:: lmt1,il1,im1,it1
        integer:: lmt2,il2,im2,it2


        do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it)
            im1=mtp(lmt1,it)
            it1=taup(lmt1,it)
            do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                if(il1==il2 .and. im1==im2) then
                    dion_kin_ion(lmt1,lmt2,it)= &
                                    kin_ae_psij(il1,it1,it2,it)+ &
                                    vionaeij(il1,it1,it2,it)- &
                                    vionpsij(il1,it1,it2,it)- &
                                    vionpsqij(il1,it1,it2,it)
                end if
            end do
        end do

        do lmt1=2,ilmt(it)
            do lmt2=1,lmt1-1
                dion_kin_ion(lmt1,lmt2,it)=dion_kin_ion(lmt2,lmt1,it)
            end do
        end do

        if(iprippex>2)then
        write(nfout,*) ' -- dion_kin_ion ---'
        do lmt1 = 1, ilmt(it)
            write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
              &               ,(dion_kin_ion(lmt1,lmt2,it),lmt2 = 1, ilmt(it))
        enddo
        endif
!do lmt1=1,ilmt(it)
!print '(9f8.5)' ,(dion_kin_ion(lmt1,lmt2,it),lmt2=1,ilmt(it))
!end do
!stop
        return
    end subroutine cnstrct_of_dion_kin_ion


   !=======================================
    subroutine make_kin_ae_ps(il,it1,it2)
   !=======================================
        integer,intent(in):: il,it1,it2

        integer:: nrc,ier,i
        real(DP),pointer,dimension(:):: pai,paj,psi,psj
        real(DP),pointer,dimension(:):: dpai,dpaj,dpsi,dpsj
        real(DP),pointer,dimension(:):: wght
        real(DP):: sum,rr,ll
! real(DP),pointer,dimension(:):: ddpaj,ddpsj
!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(pai(1:nrc))
        allocate(paj(1:nrc))
        allocate(psi(1:nrc))
        allocate(psj(1:nrc))
        allocate(dpai(1:nrc))
        allocate(dpaj(1:nrc))
        allocate(dpsi(1:nrc))
        allocate(dpsj(1:nrc))
        allocate(wght(1:nrc))
!allocate(ddpaj(1:nrc))
!allocate(ddpsj(1:nrc))
        pai(1:nrc)=psirpw(1:nrc,il,it1,it)
        paj(1:nrc)=psirpw(1:nrc,il,it2,it)
        psi(1:nrc)=phirpw(1:nrc,il,it1,it)
        psj(1:nrc)=phirpw(1:nrc,il,it2,it)

        call calc_diff_exp(ier,4,nrc,radr,pai,dpai)
        call calc_diff_exp(ier,4,nrc,radr,paj,dpaj)
        call calc_diff_exp(ier,4,nrc,radr,psi,dpsi)
        call calc_diff_exp(ier,4,nrc,radr,psj,dpsj)
!call calc_ddiff_exp(ier,5,nrc,radr,paj,dpaj,ddpaj)
!call calc_ddiff_exp(ier,5,nrc,radr,psj,dpsj,ddpsj)
        sum=0.d0
        ll=dble(il-1)
        call set_weight_exp(ier,1,nrc,radr,wght)
        do i=1,nrc
            rr=radr(i)
            sum=sum+(dpai(i)*dpaj(i)+ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
                    -dpsi(i)*dpsj(i)-ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wght(i)
        end  do
        kin_ae_psij(il,it1,it2,it)=sum/2.d0

!print *,sum/2.d0
!sum=0.d0
!ll=dble(il-1)
!call set_weight_exp(ier,1,nrc,radr,wght)
!do i=1,nrc
!    rr=radr(i)
!    sum=sum-0.5d0*(pai(i)*ddpaj(i)-ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
!            -psi(i)*ddpsj(i)+ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wght(i)
!end  do
!print *,sum
!stop

        deallocate(pai,paj,psi,psj)
        deallocate(dpai,dpaj,dpsi,dpsj)
        deallocate(wght)
        return
   end subroutine make_kin_ae_ps

   !========================================
    subroutine make_kin_ae_ps2(il,it1,it2)
   !========================================
        integer,intent(in):: il,it1,it2

        integer:: nrc,ier,ir
        real(DP),pointer,dimension(:):: wght1,wght2
        real(DP):: tmp1,tmp2
        real(DP):: vloc_ps,vloc_ae,emhnhm,emsnsm,ene,bmt

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(wght1(1:nrc))
        allocate(wght2(1:nmesh(it)))

        call set_weight_exp(ier,1,nrc,radr,wght1)
        call set_weight_exp(ier,1,nmesh(it),radr,wght2)
        vloc_ps=0.d0
        vloc_ae=0.d0
        emhnhm=0.d0
        emsnsm=0.d0
        do ir=1,nrc
            tmp1=wght1(ir)*phirpw(ir,il,it1,it)*phirpw(ir,il,it2,it)
            tmp2=wght1(ir)*psirpw(ir,il,it1,it)*psirpw(ir,il,it2,it)
            vloc_ps=vloc_ps+tmp1*vloc_scr_ps(ir)
            vloc_ae=vloc_ae+tmp2*vloc_scr_ae(ir)
            emhnhm=emhnhm+tmp1
            emsnsm=emsnsm+tmp2
        end  do
        if(it1.eq.it2) then
            ene=eps(il,it2,it)
        else
            ene=0.5d0*(eps(il,it1,it)+eps(il,it2,it))
        end if
        emhnhm=emhnhm*ene
        emsnsm=emsnsm*ene
        bmt=0.d0
        do ir=1,nmesh(it)
            bmt=bmt+wght2(ir)*phirpw(ir,il,it1,it)*chir(ir,il,it2)
        end do
        if(it1.ne.it2) then
            do ir=1,nmesh(it)
                bmt=bmt+wght2(ir)*phirpw(ir,il,it2,it)*chir(ir,il,it1)
            end do
            bmt=bmt*0.5d0
        end if
        kin_ae_psij(il,it1,it2,it)=emsnsm-vloc_ae-emhnhm+bmt+vloc_ps

!print *,sum/2.d0
!sum=0.d0
!ll=dble(il-1)
!call set_weight_exp(ier,1,nrc,radr,wos)
!do i=1,nrc
!    rr=radr(i)
!    sum=sum-0.5d0*(pai(i)*ddpaj(i)-ll*(ll+1.d0)*pai(i)*paj(i)/rr/rr &
!            -psi(i)*ddpsj(i)+ll*(ll+1.d0)*psi(i)*psj(i)/rr/rr)*wos(i)
!end  do
!print *,sum
!stop
        deallocate(wght1,wght2)
        return
   end subroutine make_kin_ae_ps2

   !======================================
    subroutine make_vionaeij(il,it1,it2)
   !======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,nrc2,ier,i
        real(DP),pointer,dimension(:):: pipj,dsum,wght
        real(DP):: sum,rr,sum1,sum2,sum0
        integer:: ir,ii,i0,j,jr,is

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        nrc2=nmesh(it)
        allocate(pipj(1:nrc))
        allocate(dsum(1:nrc))
        allocate(wght(1:nrc2))
        pipj(1:nrc)=psirpw(1:nrc,il,it1,it)* &
                    psirpw(1:nrc,il,it2,it)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+pipj(i)/radr(i)*wght(i)
        end do
        vionaeij(il,it1,it2,it)=-dble(iatomn(it))*sum

        dsum=0.d0
  !rhcorpw(i,it)
        do ir=1,nrc
            sum1=0.d0
            sum2=0.d0
            if(ir==1) then
                sum1=0.d0
            else if((ir>=2).and.(ir<=5)) then
                do ii=2,ir
                    i0=ii-1
                    is=1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum1=sum1+rhcorpw(i0+j*is,it)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,1,ir,radr,wght)
                do jr=1,ir
                    sum1=sum1+rhcorpw(jr,it)*wght(jr)
                end do
            end if
            sum1=sum1*pipj(ir)/radr(ir)
            if(ir==nrc2) then
                sum2=0.d0
            else if((ir<=nrc2-1).and.(ir>=nrc2-4)) then
                do ii=ir,nrc2-1
                    i0=ii+1
                    is=-1
                    call set_open_weight_exp(ier,i0,is,radr,wght)
                    do j=1,4
                        sum2=sum2+(1.d0/radr(i0+j*is))* &
                                rhcorpw(i0+j*is,it)*wght(i0+j*is)
                    end do
                end do
            else
                call set_weight_exp(ier,ir,nrc2,radr,wght)
                do jr=ir,nrc2
                    sum2=sum2+(1.d0/radr(jr))* &
                                rhcorpw(jr,it)*wght(jr)
                end do
            end if
            sum2=sum2*pipj(ir)
            dsum(ir)=sum1+sum2
        end do

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum0=0.d0
        do ir=1,nrc
            sum0=sum0+dsum(ir)*wght(ir)
        end do

        vionaeij(il,it1,it2,it)=vionaeij(il,it1,it2,it)+sum0

        deallocate(pipj,dsum,wght)
        return
   end subroutine make_vionaeij

   !======================================
    subroutine make_vionpsij(il,it1,it2)
   !======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,ier,i
        real(DP),pointer,dimension(:):: pipj,wght
        real(DP):: sum,rr

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(pipj(1:nrc))
        allocate(wght(1:nrc))
        pipj(1:nrc)=phirpw(1:nrc,il,it1,it)* &
                    phirpw(1:nrc,il,it2,it)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+pipj(i)*vlocr2(i)*wght(i)
        end do
        vionpsij(il,it1,it2,it)=sum

        deallocate(pipj,wght)
        return
   end subroutine make_vionpsij
     !vionpsqij     ! d(nloc,ntau,ntau,ntyp)
   !=======================================
    subroutine make_vionpsqij(il,it1,it2)
   !=======================================
        integer,intent(in):: il,it1,it2
        integer:: nrc,ier,i,iq
        real(DP),pointer,dimension(:):: qij,wght
        real(DP):: sum,rr

        iq=iqitg(il,it1,il,it2,1,it)
        if(iq.eq.0) then
            !!$print *,'il1,itau1,il2,itau2,il3 : ',il,it1,il,it2,1
            !!$print *,' is not set in iqitg ! make_vionpsqij '
            vionpsqij(il,it1,it2,it)=0.d0
            return
        end if

!        nrc=max(wf_nrc(il,it1,it),wf_nrc(il,it2,it))
        nrc=wf_mnrc(it)
        allocate(qij(1:nrc))
        allocate(wght(1:nrc))
        qij(1:nrc)=qrspspw(1:nrc,iq)

        call set_weight_exp(ier,1,nrc,radr,wght)
        sum=0.d0
        do i=1,nrc
            sum=sum+qij(i)*vlocr2(i)*wght(i)
        end do
        vionpsqij(il,it1,it2,it)=sum

        deallocate(qij,wght)
        return
   end subroutine make_vionpsqij

   !==========================================
    subroutine cnstrct_of_CijkClmkVVVVijlm_k
   !==========================================
        integer:: lmt1,il1,it1,im1,ii
        integer:: lmt2,il2,it2,im2,jj
        integer:: lmt3,il3,it3,im3,kk
        integer:: lmt4,il4,it4,im4,ll
        integer:: lmt4min,nn,n,m,mp
        integer:: il12,il34,isph12,isph34
        logical:: matchflg
        logical,save:: initialized=.false.
        real(DP):: cijclm,cijclm_ae

        real(kind=DP), pointer, dimension(:,:,:) :: cr2
        integer, pointer, dimension(:,:,:)       :: isph2
        integer, pointer, dimension(:,:)         :: mmt2

        integer:: lt1,lt2,lt3,lt4,lt5,itmp
        integer,allocatable,dimension(:)        :: ilk

        if(.not.initialized) then
            m_clmns_cijkclmk=0
            initialized=.true.
        end if

        allocate(cr2(16,16,6))
        allocate(isph2(16,16,6))
        allocate(mmt2(16,16))

        call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)

        if(.not.paramset) then
            call m_PP_find_maximum_l(n)    ! n-1: maximum l
            n = (n-1) + (n-1) + 1
            allocate(ilk(n**2)); call substitute_il3(n**2,ilk) ! -(b_Elec..)
        end if

        do lmt1=1,ilmt(it)
            il1=ltp(lmt1,it)
            im1=mtp(lmt1,it)
            it1=taup(lmt1,it)
            ii=(il1-1)**2+im1
            do lmt2=lmt1,ilmt(it)
                il2=ltp(lmt2,it)
                im2=mtp(lmt2,it)
                it2=taup(lmt2,it)
                jj=(il2-1)**2+im2

                nn=0

!                do lmt3=lmt1,ilmt(it)
!                    il3=ltp(lmt3,it)
!                    im3=mtp(lmt3,it)
!                    it3=taup(lmt3,it)
!                    kk=(il3-1)**2+im3
!
!                    lmt4min=lmt3
!                    if(lmt1==lmt3) lmt4min=max(lmt2,lmt3)
!                    do lmt4=lmt4min,ilmt(it)
!                        il4=ltp(lmt4,it)
!                        im4=mtp(lmt4,it)
!                        it4=taup(lmt4,it)
!                        ll=(il4-1)**2+im4

                do lmt3=1,ilmt(it)
                    il3=ltp(lmt3,it)
                    im3=mtp(lmt3,it)
                    it3=taup(lmt3,it)
                    kk=(il3-1)**2+im3
                    do lmt4=lmt3,ilmt(it)
                        il4=ltp(lmt4,it)
                        im4=mtp(lmt4,it)
                        it4=taup(lmt4,it)
                        ll=(il4-1)**2+im4

                        matchflg=.false.
                        cijclm=0.d0
                        cijclm_ae=0.d0
                        do n=1,mmt2(ii,jj)
                            isph12=isph2(ii,jj,n)
                            do m=1,mmt2(kk,ll)
                                isph34=isph2(kk,ll,m)
                                if(isph12==isph34) then
                                    matchflg=.true.
!print *,isph12,ilk(isph12)
                                    if(.not.paramset) then
                                        lt1=index_lmt2lt(lmt1,it)
                                        lt2=index_lmt2lt(lmt2,it)
                                        lt3=index_lmt2lt(lmt3,it)
                                        lt4=index_lmt2lt(lmt4,it)
                                        lt5=ilk(isph12)+1
!print *,'Before ',lt1,lt2,lt3,lt4,lt5
!                                        if(lt1>lt2) then
!                                            itmp=lt1
!                                            lt1=lt2
!                                            lt2=itmp
!                                        end if
!                                        if(lt3>lt4) then
!                                            itmp=lt3
!                                            lt3=lt4
!                                            lt4=itmp
!                                        end if
                                        if(lt1.gt.lt3 .or. (lt1.eq.lt3.and.lt2.gt.lt4)) then
                                            itmp=lt1
                                            lt1=lt3
                                            lt3=itmp
                                            itmp=lt2
                                            lt2=lt4
                                            lt4=itmp
                                        end if
                                        mp=ipppp(lt1,lt2,lt3,lt4,lt5,it)
!print *,'mp=',mp
                                        cijclm=cijclm+ &
                                                (vaeijlm_k(mp) &
                                                -vpsijlm_k(mp) &
                                                -vqijqlm_k(mp) &
                                                -vqijplpm_ks(mp))* &
                                                cr2(ii,jj,n)*cr2(kk,ll,m)
                                        cijclm_ae=cijclm_ae+ &
                                                vaeijlm_k(mp)* &
                                                cr2(ii,jj,n)*cr2(kk,ll,m)
!print '(4e19.6)',vaeijlm_k(mp),vpsijlm_k(mp),&
!vqijqlm_k(mp),vqijplpm_ks(mp)

                                    end if

                                end if
                            end do
                        end do

                        if(matchflg) then
                            nn=nn+1
                            if(.not.paramset) then
                                ilmt3_cijkclmk(lmt1,lmt2,nn,it)=lmt3
                                ilmt4_cijkclmk(lmt1,lmt2,nn,it)=lmt4
                                CijkClmkVVVVijlm_k(lmt1,lmt2,nn,it)=cijclm
                                CijkClmkVVVVijlm_k_ae(lmt1,lmt2,nn,it)=cijclm_ae
!print *,nn,'/',n_cijkclmk(lmt1,lmt2,it)
!print *,lmt1,lmt2,ilmt3_cijkclmk(lmt1,lmt2,nn,it),ilmt4_cijkclmk(lmt1,lmt2,nn,it)
!print '(4i5,e19.6)',lmt1,lmt2,lmt3,lmt4,cijclm
                            end if
!                            print *,nn
!                            print *,ii,jj,kk,ll
!                            print *,lmt1,lmt2,lmt3,lmt4
                        end if

                    end do
                end do
                if(.not.paramset) &
                    n_cijkclmk(lmt1,lmt2,it)=nn
                if(nn.gt.m_clmns_cijkclmk) m_clmns_cijkclmk=nn
            end do
        end do

        deallocate(cr2,isph2,mmt2)
        if(.not.paramset) deallocate(ilk)

!print *,'m_clmns_cijkclmk=',m_clmns_cijkclmk
!if(.not.paramset) stop
        return

   end subroutine cnstrct_of_CijkClmkVVVVijlm_k

   !===============================================
    subroutine cnstrct_of_CijkClmnVVVVijlm_kn
   !===============================================
        integer:: lmt1,il1,it1,im1,ii
        integer:: lmt2,il2,it2,im2,jj
        integer:: lmt3,il3,it3,im3,kk
        integer:: lmt4,il4,it4,im4,ll
        integer:: lmt4min,nn,n,m,mp
        integer:: il12,il34,isph12,isph34
        logical:: matchflg
        logical,save:: initialized=.false.
        real(DP):: cijclm,cijclm_ae

        real(kind=DP), pointer, dimension(:,:,:) :: cr2
        integer, pointer, dimension(:,:,:)       :: isph2
        integer, pointer, dimension(:,:)         :: mmt2

        integer:: lt1,lt2,lt3,lt4,lt5,itmp
        integer,allocatable,dimension(:)        :: ilk

        real(DP),allocatable :: crotylm(:,:,:)
        integer,allocatable :: iylm(:,:,:)
        integer,allocatable :: nylm(:,:)
        real(DP),allocatable :: opr(:,:,:)

        integer:: l1max,mmax,nsph,iopr
        integer:: ia,nm,l

        if(.not.initialized) then
            m_clmns_cijkclmk=0
            m_clmns_cijkclmn = 0
            initialized=.true.
        end if

        allocate(cr2(16,16,6))
        allocate(isph2(16,16,6))
        allocate(mmt2(16,16))

        call sphset2(nfout,ipri,lcmax,cr2,isph2,mmt2)

!        call m_PP_find_maximum_l(n)    ! n-1: maximum l

        n = 0
        do lmt1 = 1, ilmt(it)
            l = ltp(lmt1,it)
            if(n < l) n = l
        end do

        n = (n-1) + (n-1) + 1
        l1max=n
        mmax=2*l1max-1
        nsph=l1max**2
!        if(.not.paramset) then
        allocate(ilk(n**2)); call substitute_il3(n**2,ilk) ! -(b_Elec..)
!        end if
        allocate(crotylm(mmax,nsph,nopr))
        allocate(iylm(mmax,nsph,nopr))
        allocate(nylm(nsph,nopr))
        allocate(opr(3,3,nopr))

        do ia=1,natm
            if(ityp(ia)/=it) cycle
            do iopr=1,nopr
                if(ia2ia_symmtry_op(ia,iopr).gt.0) then
                    opr(:,:,iopr)=op(:,:,iopr)
                else
                    opr(:,:,iopr)=-op(:,:,iopr)
                end if
            end do
            call get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,opr)

            do iopr=1,nopr

!print '(a,i4)','iopr=',iopr
!print '(3e19.6)',opr(1,:,iopr)
!print '(3e19.6)',opr(2,:,iopr)
!print '(3e19.6)',opr(3,:,iopr)

                do lmt1=1,ilmt(it)
                    il1=ltp(lmt1,it)
                    im1=mtp(lmt1,it)
                    it1=taup(lmt1,it)
                    ii=(il1-1)**2+im1
                    do lmt2=lmt1,ilmt(it)
                        il2=ltp(lmt2,it)
                        im2=mtp(lmt2,it)
                        it2=taup(lmt2,it)
                        jj=(il2-1)**2+im2

                        nn=0

        !                do lmt3=lmt1,ilmt(it)
        !                    il3=ltp(lmt3,it)
        !                    im3=mtp(lmt3,it)
        !                    it3=taup(lmt3,it)
        !                    kk=(il3-1)**2+im3
        !
        !                    lmt4min=lmt3
        !                    if(lmt1==lmt3) lmt4min=max(lmt2,lmt3)
        !                    do lmt4=lmt4min,ilmt(it)
        !                        il4=ltp(lmt4,it)
        !                        im4=mtp(lmt4,it)
        !                        it4=taup(lmt4,it)
        !                        ll=(il4-1)**2+im4

                        do lmt3=1,ilmt(it)
                            il3=ltp(lmt3,it)
                            im3=mtp(lmt3,it)
                            it3=taup(lmt3,it)
                            kk=(il3-1)**2+im3
                            do lmt4=lmt3,ilmt(it)
                                il4=ltp(lmt4,it)
                                im4=mtp(lmt4,it)
                                it4=taup(lmt4,it)
                                ll=(il4-1)**2+im4

                                matchflg=.false.
                                cijclm=0.d0
                                cijclm_ae=0.d0

                                do n=1,mmt2(ii,jj)
                                    isph12=isph2(ii,jj,n)
                                    do m=1,mmt2(kk,ll)
                                        isph34=isph2(kk,ll,m)

                                        if(ilk(isph12).ne.ilk(isph34)) cycle

                                        do nm=1,nylm(isph34,iopr)
                                            if(isph12.eq.iylm(nm,isph34,iopr)) then
                                                matchflg=.true.
!                                           if(isph12==isph34) then
!                                                matchflg=.true.
!print '(2i4,e19.6)',isph12,isph34,crotylm(nm,isph34,iopr)
                                                if(.not.paramset) then
                                                    lt1=index_lmt2lt(lmt1,it)
                                                    lt2=index_lmt2lt(lmt2,it)
                                                    lt3=index_lmt2lt(lmt3,it)
                                                    lt4=index_lmt2lt(lmt4,it)
                                                    lt5=ilk(isph12)+1
            !print *,'Before ',lt1,lt2,lt3,lt4,lt5
            !                                        if(lt1>lt2) then
            !                                            itmp=lt1
            !                                            lt1=lt2
            !                                            lt2=itmp
            !                                        end if
            !                                        if(lt3>lt4) then
            !                                            itmp=lt3
            !                                            lt3=lt4
            !                                            lt4=itmp
            !                                        end if
                                                    if(lt1.gt.lt3 .or. (lt1.eq.lt3.and.lt2.gt.lt4)) then
                                                        itmp=lt1
                                                        lt1=lt3
                                                        lt3=itmp
                                                        itmp=lt2
                                                        lt2=lt4
                                                        lt4=itmp
                                                    end if
                                                    mp=ipppp(lt1,lt2,lt3,lt4,lt5,it)
            !print *,'mp=',mp
                                                    cijclm=cijclm+ &
                                                            (vaeijlm_k(mp) &
                                                            -vpsijlm_k(mp) &
                                                            -vqijqlm_k(mp) &
                                                            -vqijplpm_ks(mp))* &
                                                            cr2(ii,jj,n)*cr2(kk,ll,m)* &
                                                            crotylm(nm,isph34,iopr)
                                                    cijclm_ae=cijclm_ae+ &
                                                            vaeijlm_k(mp)* &
                                                            cr2(ii,jj,n)*cr2(kk,ll,m)* &
                                                            crotylm(nm,isph34,iopr)
            !print '(4e19.6)',vaeijlm_k(mp),vpsijlm_k(mp),&
            !vqijqlm_k(mp),vqijplpm_ks(mp)

                                                end if
                                            end if
                                        end do

                                    end do
                                end do

                                if(matchflg) then
                                    nn=nn+1
                                    if(.not.paramset) then
                                        ilmt3_cijkclmn(lmt1,lmt2,nn,ia,iopr)=lmt3
                                        ilmt4_cijkclmn(lmt1,lmt2,nn,ia,iopr)=lmt4
                                        CijkClmnVVVVijlm_kn(lmt1,lmt2,nn,ia,iopr)=cijclm
                                        CijkClmnVVVVijlm_kn_ae(lmt1,lmt2,nn,ia,iopr)=cijclm_ae
        !print *,nn,'/',n_cijkclmk(lmt1,lmt2,it)
        !print *,lmt1,lmt2,ilmt3_cijkclmk(lmt1,lmt2,nn,it),ilmt4_cijkclmk(lmt1,lmt2,nn,it)
        !print '(4i5,e19.6)',lmt1,lmt2,lmt3,lmt4,cijclm
                                    end if
        !                            print *,nn
        !                            print *,ii,jj,kk,ll
        !                            print *,lmt1,lmt2,lmt3,lmt4
                                end if

                            end do
                        end do
                        if(.not.paramset) &
                            n_cijkclmn(lmt1,lmt2,ia,iopr)=nn
!print *,'lmt1,lmt2,ia,iopr,nn',lmt1,lmt2,ia,iopr,nn
                        if(nn.gt.m_clmns_cijkclmn) m_clmns_cijkclmn=nn
                    end do
                end do

            end do

        end do

!if(.not.paramset) stop
!print *,m_clmns_cijkclmn
!stop
        deallocate(cr2,isph2,mmt2)
!        if(.not.paramset) deallocate(ilk)
        deallocate(ilk)
        deallocate(crotylm,iylm,nylm,opr)
!        deallocate(ia2ia_symmtry_op)
!print *,'m_clmns_cijkclmk=',m_clmns_cijkclmk
!if(.not.paramset) stop
        return

   end subroutine cnstrct_of_CijkClmnVVVVijlm_kn

   !======================================
    subroutine make_ia2ia_symmtry_op_etc
   !======================================
      if(.not.associated(ia2ia_symmtry_op)) then
         allocate(ia2ia_symmtry_op(natm,nopr+af))
         ia2ia_symmtry_op=0
         call set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia_symmtry_op)
      end if
      if(.not.associated(ia2ia_symmtry_op_inv)) then
         allocate(ia2ia_symmtry_op_inv(natm,nopr+af))
         ia2ia_symmtry_op_inv=0
         call set_ia2ia_symmtry_op_inv &
              (natm,nopr,ia2ia_symmtry_op,ia2ia_symmtry_op_inv)
      end if
      return
   end subroutine make_ia2ia_symmtry_op_etc

   !==========================================
   subroutine init_of_dion_paw
   !==========================================
     integer:: ia,is,lmt1,lmt2
     integer :: ismax

     if ( noncol ) then
        ismax = 1
     else
        ismax = nspin
     endif

     do ia=1,natm
        if(ityp(ia)/=it) cycle       ! ASMS 2019/12/03
! ================================ modified by K. Tagami ========== 11.0
!        do is=1,nspin
        do is=1, ismax
! ================================================================= 11.0
           do lmt1=1,ilmt(it)
              do lmt2=1,ilmt(it)
                 dion_paw(lmt1,lmt2,is,ia)=dion(lmt1,lmt2,it)
              end do
           end do
        end do
     end do
     return
   end subroutine init_of_dion_paw

!....................................................................
  end subroutine m_PP_vanderbilt_type_gncpm2
!$$#endif
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
!!BRANCH_Q_END ORG_Parallel
! ==================================================================================================

!#####################################################################

   !====================================
    subroutine m_PP_get_dion_paw(nfout)
   !====================================
        integer,intent(in):: nfout
        integer:: ia,it,is,lmt1,lmt2

        do ia=1,natm
            it=ityp(ia)
            if(ipaw(it)/=1) cycle
            do is=1,nspin
                do lmt1=1,ilmt(it)
                    do lmt2=lmt1,ilmt(it)
                        dion_paw(lmt1,lmt2,is,ia)=  dion_kin_ion(lmt1,lmt2,it)+ &
                                                    dion_hartree(lmt1,lmt2,ia)+ &
                                                    dion_vxc(lmt1,lmt2,is,ia)
                    end do
                end do
                do lmt1=2,ilmt(it)
                    do lmt2=1,lmt1-1
                        dion_paw(lmt1,lmt2,is,ia)=dion_paw(lmt2,lmt1,is,ia)
                    end do
                end do
            end do
        end do

        if(ipripp>=2.and.printable)then
        write(nfout,*)
        write(nfout,*) ' -- dion_paw ---'
        do ia=1,natm
            it=ityp(ia)
            if(ipaw(it)/=1) cycle
            do is=1,nspin
                write(nfout,'(a,i2,a,i2,a)') '(ia,is)=(',ia,',',is,')'
                do lmt1 = 1, ilmt(it)
                    write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
                      &               ,(dion_paw(lmt1,lmt2,is,ia),lmt2 = 1, ilmt(it))
                enddo
            end do
        end do
        endif

        return
    end subroutine m_PP_get_dion_paw

! ======================================== added by K. Tagami ============ 11.0
  subroutine m_PP_get_dion_paw_noncl(nfout)
    integer,intent(in):: nfout
    integer:: ia,it,is,lmt1,lmt2

    do ia=1,natm
       it=ityp(ia)
       if(ipaw(it)/=1) cycle
       do is=1, ndim_magmom
          if ( is == 1 ) then
             do lmt1=1,ilmt(it)
                do lmt2=lmt1,ilmt(it)
                   dion_paw(lmt1,lmt2,is,ia)=  dion_kin_ion(lmt1,lmt2,it)+ &
                        dion_hartree(lmt1,lmt2,ia)+ &
                        dion_vxc(lmt1,lmt2,is,ia)
                end do
             end do
          else
             do lmt1=1,ilmt(it)
                do lmt2=lmt1,ilmt(it)
                   dion_paw(lmt1,lmt2,is,ia) = dion_vxc(lmt1,lmt2,is,ia)
                end do
             end do
          endif
          do lmt1=2,ilmt(it)
             do lmt2=1,lmt1-1
                dion_paw(lmt1,lmt2,is,ia)=dion_paw(lmt2,lmt1,is,ia)
             end do
          end do
       end do
    end do

    if(ipripp>=2.and.printable)then
       write(nfout,*)
       write(nfout,*) ' -- dion_paw ---'
       do ia=1,natm
          it=ityp(ia)
          if(ipaw(it)/=1) cycle
          do is=1,nspin
             write(nfout,'(a,i2,a,i2,a)') '(ia,is)=(',ia,',',is,')'
             do lmt1 = 1, ilmt(it)
                write(nfout,'(i3,15f12.9/15f12.9)') lmt1 &
                     &               ,(dion_paw(lmt1,lmt2,is,ia),lmt2 = 1, ilmt(it))
             enddo
          end do
       end do
    endif

    return
  end subroutine m_PP_get_dion_paw_noncl
! ================================================================= 11.0

  subroutine m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
    integer, intent(in) :: nfout,nfcntn_bin_paw
    integer:: ierr, itmp
    integer             :: id_sname = -1, ityp,itau
    call tstatc0_begin('m_PP_rd_PAW_parameters ',id_sname)

!!$    if(mype==0) then
!!$        if(.not.flg_symmtry) then
!!$           read(nfcntn_bin_paw) &
!!$                &  ipaw,chgcr,iltpw,lppw,tppw,psirpw,phirpw,qrspspw &
!!$                &, rhcorpw,rhpcrpw,wf_nrc,wf_mnrc,n_cijkclmk &
!!$                &, CijkClmkVVVVijlm_k,ilmt3_cijkclmk,ilmt4_cijkclmk &
!!$                &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
!!$                &, dion_vxc,dion_paw,radr_paw,m_clmns_cijkclmk,CijkClmkVVVVijlm_k_ae
!!$        else
!!$           read(nfcntn_bin_paw) &
!!$                &  ipaw,chgcr,iltpw,lppw,tppw,psirpw,phirpw,qrspspw &
!!$                &, rhcorpw,rhpcrpw,wf_nrc,wf_mnrc,n_cijkclmn &
!!$                &, CijkClmnVVVVijlm_kn,ilmt3_cijkclmn,ilmt4_cijkclmn &
!!$                &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
!!$                &, dion_vxc,dion_paw,radr_paw,m_clmns_cijkclmn,CijkClmnVVVVijlm_kn_ae
!!$        end if
!!$    end if

    if(.not.flg_symmtry) then
       if(mype==0) then
          if ( cntn_bin_paw_format_is_set .and. cntn_bin_paw_format == 1 ) then
             read(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw, itmp
          else
#ifdef _PAW_CONTINUE_DATA_PREVIOUS_BEFORE_201403_STYLE_
             read(nfcntn_bin_paw) &
                  &  ipaw,chgcr,iltpw,lppw,tppw &
                  &, wf_nrc,wf_mnrc,n_cijkclmk &
                  &, CijkClmkVVVVijlm_k,ilmt3_cijkclmk,ilmt4_cijkclmk &
                  &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmk,CijkClmkVVVVijlm_k_ae
#else
             read(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
               &, dion_vxc,dion_paw,m_clmns_cijkclmk
#endif
          end if
       endif
    else
       if(mype==0) then
          if ( cntn_bin_paw_format_is_set .and. cntn_bin_paw_format == 1 ) then
             read(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw, itmp
          else
#ifdef _PAW_CONTINUE_DATA_PREVIOUS_BEFORE_201403_STYLE_
             read(nfcntn_bin_paw) &
                  &  ipaw,chgcr,iltpw,lppw,tppw &
                  &, wf_nrc,wf_mnrc,n_cijkclmn &
                  &, CijkClmnVVVVijlm_kn,ilmt3_cijkclmn,ilmt4_cijkclmn &
                  &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn,CijkClmnVVVVijlm_kn_ae
#else
             read(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn
#endif
          endif
       end if
    end if

    call bcast_nfcntn_bin_paw                !-(c.h.)

    do ityp=1,ntyp
       do itau=1,ntau
          call rd_mmesh_array(mmesh,nloc,psirpw(1,1,itau,ityp))
       end do
    end do
    do ityp=1,ntyp
       do itau=1,ntau
          call rd_mmesh_array(mmesh,nloc,phirpw(1,1,itau,ityp))
       end do
    end do
    call rd_mmesh_array(mmesh,nqitg+1,qrspspw(1,0))
    call rd_mmesh_array(mmesh,ntyp,rhcorpw)
    call rd_mmesh_array(mmesh,ntyp,rhpcrpw)
    call rd_mmesh_array(mmesh,ntyp,radr_paw)

    call tstatc0_end(id_sname)

  contains

    subroutine rd_mmesh_array(n1,n2,reduced_array)
      integer,intent(in):: n1,n2
!!$      real(kind=DP),intent(out),dimension(ista:iend,n2) :: reduced_array
      real(kind=DP),intent(out),dimension(1:mmesh,n2) :: reduced_array
      real(DP),allocatable,dimension(:,:) :: a_mpi
      real(DP),allocatable,dimension(:) :: b_mpi
      integer                         :: i, j, ista,iend
      ista=1
      iend=mmesh

      if(n2 >= 1) then
         allocate(a_mpi(n1,n2))
         if(mype == 0) read(nfcntn_bin_paw) a_mpi

         if(n1*n2*8 > MAXIMUM_BCAST_SIZE) then
            if(ipripp >= 2) write(nfout,'(" n1, n2 = ",i10,",",i8, &
                 & "  ( n1*n2*8 > MAXIMUM_BCAST_SIZE = ",i12,") <<rd_mmesh_array>>")') &
                 & n1, n2, MAXIMUM_BCAST_SIZE
            allocate(b_mpi(n1))
            do j = 1, n2
               b_mpi(:) = a_mpi(:,j)
               call mpi_bcast(b_mpi,n1,mpi_double_precision,0,mpi_comm_world,ierr)
               reduced_array(ista:iend,j) = b_mpi(ista:iend)
            end do
            deallocate(b_mpi)
         else
            call mpi_bcast(a_mpi,n1*n2,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
            do j = 1, n2
               reduced_array(ista:iend,j) = a_mpi(ista:iend,j)  ! MPI
            end do
         end if
         deallocate(a_mpi)
      end if
    end subroutine rd_mmesh_array

    subroutine bcast_nfcntn_bin_paw
      real(DP), allocatable,dimension(:,:,:,:) :: b_mpi
      integer,  allocatable,dimension(:,:,:,:) :: bi_mpi
      integer :: nsize, i, j, natm_d, iloop, is, ie

      call mpi_bcast(ipaw,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(chgcr,ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(iltpw,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(lppw,nltpw*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(tppw,nltpw*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
!!$        call mpi_bcast(psirpw,mmesh*nloc*ntau*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$        call mpi_bcast(phirpw,mmesh*nloc*ntau*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$        call mpi_bcast(qrspspw,mmesh*nqitg,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$        call mpi_bcast(rhcorpw,mmesh*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
!!$        call mpi_bcast(rhpcrpw,mmesh*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(wf_nrc,nloc*ntau*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(wf_mnrc,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      if(.not.flg_symmtry) then
         call mpi_bcast(m_clmns_cijkclmk,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(n_cijkclmk,nlmt*nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
         nsize = nlmt*nlmt*m_clmns_cijkclmk*ntyp
         call mpi_bcast(CijkClmkVVVVijlm_k,    nsize, mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(CijkClmkVVVVijlm_k_ae, nsize, mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(ilmt3_cijkclmk,        nsize, mpi_integer,         0,MPI_CommGroup,ierr)
         call mpi_bcast(ilmt4_cijkclmk,        nsize, mpi_integer,         0,MPI_CommGroup,ierr)
      else
         call mpi_bcast(m_clmns_cijkclmn,1,mpi_integer,0,MPI_CommGroup,ierr)
         call mpi_bcast(n_cijkclmn,nlmt*nlmt*natm*nopr,mpi_integer,0,MPI_CommGroup,ierr)
         nsize = nlmt*nlmt*m_clmns_cijkclmn*natm*nopr
         if(nsize*8 > MAXIMUM_BCAST_SIZE) then
            natm_d = max(1,min(MAXIMUM_BCAST_SIZE/(nlmt*nlmt*m_clmns_cijkclmn*nopr*8),natm))
            allocate (b_mpi(nlmt,nlmt,m_clmns_cijkclmn,natm_d))
            allocate (bi_mpi(nlmt,nlmt,m_clmns_cijkclmn,natm_d))
            nsize = nlmt*nlmt*m_clmns_cijkclmn
            do j = 1, nopr
               is = 1; ie = natm_d
               iloop = ceiling(dble(natm)/natm_d)
               do i = 1, iloop
                  if(mype==0) b_mpi(:,:,:,1:ie-is+1) = CijkClmnVVVVijlm_kn(:,:,:,is:ie,j)
                  call mpi_bcast(b_mpi,nsize*(ie-is+1),mpi_double_precision,0,MPI_CommGroup,ierr)
                  if(mype/=0) CijkClmnVVVVijlm_kn(:,:,:,is:ie,j) = b_mpi(:,:,:,1:ie-is+1)

                  if(mype==0) b_mpi(:,:,:,1:ie-is+1) = CijkClmnVVVVijlm_kn_ae(:,:,:,is:ie,j)
                  call mpi_bcast(b_mpi,nsize*(ie-is+1),mpi_double_precision,0,MPI_CommGroup,ierr)
                  if(mype/=0) CijkClmnVVVVijlm_kn_ae(:,:,:,is:ie,j) = b_mpi(:,:,:,1:ie-is+1)

                  if(mype==0) bi_mpi(:,:,:,1:ie-is+1) = ilmt3_cijkclmn(:,:,:,is:ie,j)
                  call mpi_bcast(bi_mpi,nsize*(ie-is+1),mpi_integer,0,MPI_CommGroup,ierr)
                  if(mype/=0) ilmt3_cijkclmn(:,:,:,is:ie,j) = bi_mpi(:,:,:,1:ie-is+1)

                  if(mype==0) bi_mpi(:,:,:,1:ie-is+1) = ilmt4_cijkclmn(:,:,:,is:ie,j)
                  call mpi_bcast(bi_mpi,nsize*(ie-is+1),mpi_integer,0,MPI_CommGroup,ierr)
                  if(mype/=0) ilmt4_cijkclmn(:,:,:,is:ie,j) = bi_mpi(:,:,:,1:ie-is+1)

                  is = ie+1; ie = min(ie+natm_d,natm)
               end do
            end do
            deallocate(bi_mpi)
            deallocate(b_mpi)
         else
            call mpi_bcast(CijkClmnVVVVijlm_kn,   nsize, mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(CijkClmnVVVVijlm_kn_ae,nsize, mpi_double_precision,0,MPI_CommGroup,ierr)
            call mpi_bcast(ilmt3_cijkclmn,        nsize, mpi_integer,         0,MPI_CommGroup,ierr)
            call mpi_bcast(ilmt4_cijkclmn,        nsize, mpi_integer,         0,MPI_CommGroup,ierr)
         end if
      end if
      call mpi_bcast(index_lmt2lt,nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(dion_kin_ion,nlmt*nlmt*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(dion_hartree,nlmt*nlmt*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(dion_hartree_now,nlmt*nlmt*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
! ================================== modified by K. Tagami =========== 11.0
!        call mpi_bcast(dion_vxc,nlmt*nlmt*nspin*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
!        call mpi_bcast(dion_paw,nlmt*nlmt*nspin*natm,mpi_double_precision,0,MPI_CommGroup,ierr)

      call mpi_bcast( dion_vxc, nlmt*nlmt*ndim_magmom*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast( dion_paw, nlmt*nlmt*ndim_magmom*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
! ========================================================================= 11.0

!!$        call mpi_bcast(radr_paw,mmesh*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)

    end subroutine bcast_nfcntn_bin_paw
  end subroutine m_PP_rd_PAW_parameters

  subroutine m_PP_wd_PAW_parameters(nfout,nfcntn_bin_paw)
    integer, intent(in) :: nfout,nfcntn_bin_paw
#ifdef _PAW_DEBUG_WRITE_
    integer :: it, lmt1, ilt1
#endif
    integer :: ityp, itau

!!$    if(mype==0) then
!!$       if(.not.flg_symmtry) then
!!$           write(nfcntn_bin_paw) &
!!$                &  ipaw,chgcr,iltpw,lppw,tppw,psirpw,phirpw,qrspspw &
!!$                &, rhcorpw,rhpcrpw,wf_nrc,wf_mnrc,n_cijkclmk &
!!$                &, CijkClmkVVVVijlm_k,ilmt3_cijkclmk,ilmt4_cijkclmk &
!!$                &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
!!$                &, dion_vxc,dion_paw,radr_paw,m_clmns_cijkclmk,CijkClmkVVVVijlm_k_ae
!!$        else
!!$           write(nfcntn_bin_paw) &
!!$                &  ipaw,chgcr,iltpw,lppw,tppw,psirpw,phirpw,qrspspw &
!!$                &, rhcorpw,rhpcrpw,wf_nrc,wf_mnrc,n_cijkclmn &
!!$                &, CijkClmnVVVVijlm_kn,ilmt3_cijkclmn,ilmt4_cijkclmn &
!!$                &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
!!$                &, dion_vxc,dion_paw,radr_paw,m_clmns_cijkclmn,CijkClmnVVVVijlm_kn_ae
!!$        end if
!!$    end if

#ifdef _PAW_DEBUG_WRITE_
    if(ipripp>=1) then
       write(nfout,'(" **** subroutine <<m_PP_wd_PAW_parameters(nfout,nfcntn_bin_paw)>> ****")')
       write(nfout,'(" m_clmns_cijkclmn  = ",i12)') m_clmns_cijkclmn
       write(nfout,'(" *** sizes of arrays ***")')
       write(nfout,'(" size of dion_hartree           = ",i10)') size(dion_hartree)
       write(nfout,'(" size of dion_hartree_now       = ",i10)') size(dion_hartree_now)
       write(nfout,'(" size of dion_vxc               = ",i10)') size(dion_vxc)
       write(nfout,'(" size of dion_paw               = ",i10)') size(dion_paw)
       call flush(nfout)
    end if
#endif

    if(.not.flg_symmtry) then
       if(mype==0) then
          if ( cntn_bin_paw_format_is_set .and. cntn_bin_paw_format == 1 ) then
             write(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn
          else
#ifdef _PAW_CONTINUE_DATA_PREVIOUS_BEFORE_201403_STYLE_
             write(nfcntn_bin_paw) &
                  &  ipaw,chgcr,iltpw,lppw,tppw &
                  &, wf_nrc,wf_mnrc,n_cijkclmk &
                  &, CijkClmkVVVVijlm_k,ilmt3_cijkclmk,ilmt4_cijkclmk &
                  &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmk,CijkClmkVVVVijlm_k_ae
#else
             write(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn
!!$          write(nfcntn_bin_paw) dion_hartree
!!$          write(nfcntn_bin_paw) dion_hartree_now
!!$          write(nfcntn_bin_paw) dion_vxc
!!$          write(nfcntn_bin_paw) dion_paw
!!$          write(nfcntn_bin_paw) m_clmns_cijkclmk
#endif
          endif
       end if

    else
       if(mype==0) then
          if ( cntn_bin_paw_format_is_set .and. cntn_bin_paw_format == 1 ) then
             write(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn
          else
#ifdef _PAW_CONTINUE_DATA_PREVIOUS_BEFORE_201403_STYLE_
             write(nfcntn_bin_paw) &
                  &  ipaw,chgcr,iltpw,lppw,tppw &
                  &, wf_nrc,wf_mnrc,n_cijkclmn &
                  &, CijkClmnVVVVijlm_kn,ilmt3_cijkclmn,ilmt4_cijkclmn &
                  &, index_lmt2lt,dion_kin_ion,dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn,CijkClmnVVVVijlm_kn_ae
#else
             write(nfcntn_bin_paw) dion_hartree,dion_hartree_now &
                  &, dion_vxc,dion_paw,m_clmns_cijkclmn
!!$          write(nfcntn_bin_paw) dion_hartree
!!$          write(nfcntn_bin_paw) dion_hartree_now
!!$          write(nfcntn_bin_paw) dion_vxc
!!$          write(nfcntn_bin_paw) dion_paw
!!$          write(nfcntn_bin_paw) m_clmns_cijkclmk
#endif
          endif
       end if
    end if

    do ityp=1,ntyp
       do itau=1,ntau
          call wd_mmesh_array(mmesh,nloc,psirpw(1,1,itau,ityp))
       end do
    end do
    do ityp=1,ntyp
       do itau=1,ntau
          call wd_mmesh_array(mmesh,nloc,phirpw(1,1,itau,ityp))
       end do
    end do
    call wd_mmesh_array(mmesh,nqitg+1,qrspspw(1,0))
    call wd_mmesh_array(mmesh,ntyp,rhcorpw)
    call wd_mmesh_array(mmesh,ntyp,rhpcrpw)
    call wd_mmesh_array(mmesh,ntyp,radr_paw)

  contains
    subroutine wd_mmesh_array(n1,n2,reduced_array)
      integer, intent(in)                         :: n1, n2
      real(DP),intent(in),dimension(1:mmesh,n2) :: reduced_array

      real(DP),allocatable,dimension(:,:)         :: a_mpi,b_mpi,a
      real(DP),allocatable,dimension(:)           :: c1,c2
      integer                                     :: i,j,ista,iend

      ista=1
      iend=n1

      if(n2 >= 1) then
         if(n1*n2*8 > MAXIMUM_MPI_SIZE) then
            if(ipripp >= 2) write(nfout,'(" n1, n2 = ",i10,",",i8, &
                 & "  ( n1*n2*8 > MAXIMUM_MPI_SIZE = ",i12,") <<wd_kgp_array>>")') &
                 & n1, n2, MAXIMUM_MPI_SIZE
            allocate(a(n1,n2),c1(n1),c2(n1))
#ifdef _PAW_DEBUG_WRITE_
            if(ipripp>=1) then
               write(nfout,'(" size of a(:,:)                 = ",i10)') size(a)
               call flush(nfout)
            end if
#endif
            a = 0.0d0
            c1 = 0.0d0
            do j = 1, n2
               do i = ista, iend
                  c1(i) = reduced_array(i,j)
               end do
               if(iend-ista+1 .eq. n1) then
                  c2 = c1
               else
                  c2 = 0.0d0
!!$               call mpi_allreduce(c1,c2,n1,mpi_double_precision,mpi_sum &
!!$                    &            , mpi_nrc_world, ierr)
                  call mpi_allreduce(c1,c2,n1,mpi_double_precision,mpi_sum &
                       &            , MPI_CommGroup, ierr)
               end if
               if(mype == 0) a(:,j) = c2(:)
            end do
            if(mype == 0) write(nfcntn_bin_paw) a
            deallocate(c2,c1,a)
         else
            allocate(a_mpi(n1,n2)); a_mpi = 0.d0
            allocate(b_mpi(n1,n2))
#ifdef _PAW_DEBUG_WRITE_
            if(ipripp>=1) then
               write(nfout,'(" size of b(:,:)                 = ",i10)') size(b_mpi)
               call flush(nfout)
            end if
#endif
            do j = 1, n2
               do i = ista, iend                     !for mpi
                  a_mpi(i,j) = reduced_array(i,j)
               enddo
            enddo
            if(iend-ista+1 .eq. n1) then
               b_mpi = a_mpi
            else
!!$            call mpi_allreduce(a_mpi,b_mpi,n1*n2,mpi_double_precision,mpi_sum & ! MPI
!!$                 &                             ,mpi_nrc_world,ierr)  ! MPI
               call mpi_allreduce(a_mpi,b_mpi,n1*n2,mpi_double_precision,mpi_sum & ! MPI
                    &                             ,MPI_CommGroup,ierr)  ! MPI
            end if
            if(mype == 0) write(nfcntn_bin_paw) b_mpi
            deallocate(a_mpi); deallocate(b_mpi)
         end if
      end if
    end subroutine wd_mmesh_array
  end subroutine m_PP_wd_PAW_parameters

  subroutine m_PP_dealloc_paw
    if(allocated(ipaw)) deallocate(ipaw)
    if(allocated(chgcr)) deallocate(chgcr)
    if(allocated(iltpw)) deallocate(iltpw)
    if(allocated(lppw)) deallocate(lppw)
    if(allocated(tppw)) deallocate(tppw)
    if(allocated(psirpw)) deallocate(psirpw)
    if(allocated(phirpw)) deallocate(phirpw)
    if(allocated(qrspspw)) deallocate(qrspspw)
    if(allocated(rhcorpw)) deallocate(rhcorpw)
    if(allocated(rhpcrpw)) deallocate(rhpcrpw)
    if(allocated(wf_nrc)) deallocate(wf_nrc)
    if(allocated(wf_mnrc)) deallocate(wf_mnrc)
    if(.not.flg_symmtry) then
        if(allocated(n_cijkclmk)) deallocate(n_cijkclmk)
        if(allocated(CijkClmkVVVVijlm_k)) deallocate(CijkClmkVVVVijlm_k)
        if(allocated(CijkClmkVVVVijlm_k_ae)) deallocate(CijkClmkVVVVijlm_k_ae)
        if(allocated(ilmt3_cijkclmk)) deallocate(ilmt3_cijkclmk)
        if(allocated(ilmt4_cijkclmk)) deallocate(ilmt4_cijkclmk)
    else
        if(allocated(n_cijkclmn)) deallocate(n_cijkclmn)
        if(allocated(CijkClmnVVVVijlm_kn)) deallocate(CijkClmnVVVVijlm_kn)
        if(allocated(CijkClmnVVVVijlm_kn_ae)) deallocate(CijkClmnVVVVijlm_kn_ae)
        if(allocated(ilmt3_cijkclmn)) deallocate(ilmt3_cijkclmn)
        if(allocated(ilmt4_cijkclmn)) deallocate(ilmt4_cijkclmn)
    end if
    if(allocated(index_lmt2lt)) deallocate(index_lmt2lt)
    if(allocated(dion_kin_ion)) deallocate(dion_kin_ion)
    if(allocated(dion_hartree)) deallocate(dion_hartree)
    if(allocated(dion_hartree_now)) deallocate(dion_hartree_now)
    if(allocated(dion_vxc)) deallocate(dion_vxc)
    if(allocated(dion_paw)) deallocate(dion_paw)
    if(allocated(radr_paw)) deallocate(radr_paw)
    if(associated(ia2ia_symmtry_op)) deallocate(ia2ia_symmtry_op)
    if(associated(ia2ia_symmtry_op_inv)) deallocate(ia2ia_symmtry_op_inv)

! ============= KT_add ======== 13.0U2
    if(allocated(dion_paw_old)) deallocate(dion_paw_old)
! ============================= 13.0U2

! === KT_add === 2014/08/11
    if ( allocated(vlocr_pw) ) deallocate( vlocr_pw )
! ============== 2014/08/11

! === ASMS_DEBUG == 2013/02/07
    if ( .not. associated(crotylm_paw) ) return
                         ! sw_phonon == on .and. sw_calc_force == off
! === ASMS_DEBUG == 2013/02/07

    deallocate(crotylm_paw)
    deallocate(iylm_paw)
    deallocate(nylm_paw)
    return
  end subroutine m_PP_dealloc_paw


! === DEBUG by tkato 2011/09/22 ================================================
#if 0
! ==============================================================================
    subroutine set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia)
        integer,intent(in):: nfout,natm,nopr,af
        integer,intent(out):: ia2ia(natm,nopr+af)
        integer:: ia,it,no,ja,jt
        integer:: i,j,k
        real(DP):: pos0(3),pos1(3),pos2(3)
        real(DP):: distance
        real(kind=DP), allocatable, dimension(:,:) :: pos_t(:,:)
        real(kind=DP), allocatable,dimension(:,:,:)  :: op_pr

        allocate(op_pr(3,3,nopr+af))
        call m_CS_op_in_PUCD(nfout,op_pr,nopr+af)

        allocate(pos_t(natm,3))
        pos_t(1:natm,1:3) = pos(1:natm,1:3)
        do ia = 1, natm
           do i = 1, 3
              pos_t(ia,i) = pos_t(ia,i) - floor(pos_t(ia,i))
           end do
        end do

        do ia=1,natm
           it=ityp(ia)
!            if(ipaw(it)/=1) cycle
           pos0(1:3)=pos_t(ia,1:3)

!!$           write(nfout,*)
!!$           write(nfout,'(" ia = ",i8, "( ",3f12.6, " )")') ia, pos0(1:3)

           OPLoop: do no=1,nopr+af
              pos1(:)=matmul(op_pr(:,:,no),pos0(:))+tau(:,no,BUCS)
              pos1(:) = pos1(:) - floor(pos1(:))
!!$              write(nfout,'("  no = ",i8, "( ",3f12.6, " )")') no, pos1(1:3)

              ATOMLoop: do ja = 1, natm
                 jt=ityp(ja)
                 if(nint(iatomn(it)) /= nint(iatomn(jt))) cycle

                 pos2(1:3)=pos_t(ja,1:3)

                 distance=abs(pos1(1)-pos2(1))+abs(pos1(2)-pos2(2))+abs(pos1(3)-pos2(3))
                 if(distance < 1.d-5) then
                    ia2ia(ia,no)=ja
                    exit ATOMLoop
                 end if

                 if(kimg==1 .and. iwei(ja)==2) then
                    pos2(1:3)=-pos_t(ja,1:3)
                    do i = 1, 3
                       if(pos2(i) < 0.d0) pos2(i) = pos2(i) + 1.d0
                    end do
                    distance=abs(pos1(1)-pos2(1))+abs(pos1(2)-pos2(2))+abs(pos1(3)-pos2(3))
                    if(distance < 1.d-5) then
                       ia2ia(ia,no)=-ja
                       exit ATOMLoop
                    end if
                 end if

                 if(ja==natm) then
                    write(nfout,*) 'Error in set_ia2ia_symmtry_op !'
                    write(nfout,*) 'ia,no =', ia,no
                    call phase_error_with_msg(nfout,'Error in set_ia2ia_symmtry_op !',__LINE__,__FILE__)
                 end if

              end do ATOMLoop
           end do OPLoop

        end do
        deallocate(pos_t)
        deallocate(op_pr)
        return
    end subroutine set_ia2ia_symmtry_op
! === DEBUG by tkato 2011/09/22 ================================================
#else
    subroutine set_ia2ia_symmtry_op(nfout,natm,nopr,af,ia2ia)
        integer,intent(in):: nfout,natm,nopr,af
        integer,intent(out):: ia2ia(natm,nopr+af)
        integer:: ia,it,no,ja,jt, ja2
        integer:: i,j,k
        real(DP):: pos0(3),pos1(3),pos2(3)
        real(DP):: distance, dpos(3)
        real(kind=DP), allocatable, dimension(:,:) :: pos_t(:,:)
        real(kind=DP), allocatable,dimension(:,:,:)  :: op_pr

        allocate(op_pr(3,3,nopr+af))
        call m_CS_op_in_PUCD(nfout,op_pr,nopr+af)

        allocate(pos_t(natm,3))
        pos_t(1:natm,1:3) = pos(1:natm,1:3)
        do ia = 1, natm
           do i = 1, 3
              pos_t(ia,i) = pos_t(ia,i) - floor(pos_t(ia,i))
           end do
        end do

#ifdef _DEBUG_WRITE_
        do ia=1,natm
           write(nfout,'(" ia = ",i8, "( ",3d20.12, " )")') ia, pos_t(ia,1:3)
           write(nfout,'("      ",8x, "( ",3d20.12, " )")')     pos(ia,1:3)
        end do
#endif

        do ia=1,natm
           it=ityp(ia)
!            if(ipaw(it)/=1) cycle
           pos0(1:3)=pos_t(ia,1:3)

#ifdef _DEBUG_WRITE_
           write(nfout,'(" ia = ",i8, "( ",3f12.6, " )")') ia, pos0(1:3)
#endif

           OPLoop: do no=1,nopr+af
              pos1(:)=matmul(op_pr(:,:,no),pos0(:))+tau(:,no,BUCS)
              pos1(:) = pos1(:) - floor(pos1(:))
#ifdef _DEBUG_WRITE_
              write(nfout,'("  no = ",i8, "( ",3f12.6, " )")') no, pos1(1:3)
#endif

              ATOMLoop: do ja = 1, natm
                 jt=ityp(ja)
                 if(nint(iatomn(it)) /= nint(iatomn(jt))) cycle

                 pos2(1:3)=pos_t(ja,1:3)

                 dpos = abs(pos1-pos2)
                 dpos = dpos - floor(dpos + DELTA10)

                 distance = sum(dpos)
#ifdef _DEBUG_WRITE_
                 write(nfout,'("  ja = ",i8, " dpos = ( ",3f12.6, " ), distance =  ",f12.6)') ja, dpos(1:3), distance
                 do i = 1, 3
                    if(dpos(i) < 0.d0) write(nfout,'("     negative  : dpos(", i," ) = ",d20.8)') i, dpos(i)
                    if(1.d0-1.d-14 <  dpos(i) .and. dpos(i) < 1.d0) write(nfout,'("  smaller than 1.d0 : dpos(", i , " ) = ",d20.8)') i, dpos(i)
                 end do
#endif
                 if(distance < 1.d-5) then
                    ia2ia(ia,no)=ja
                    exit ATOMLoop
                 end if

                 if(kimg==1 .and. iwei(ja)==2) then
                    pos2(1:3)=-pos_t(ja,1:3)
                    dpos = abs(pos1-pos2)
                    dpos = dpos - floor(dpos+DELTA10)
                    distance = sum(dpos)
                    if(distance < 1.d-5) then
                       ia2ia(ia,no)=-ja
                       exit ATOMLoop
                    end if
                 end if

                 if(ja==natm) then
                    write(nfout,*) 'Error in set_ia2ia_symmtry_op !'
                    write(nfout,*) 'ia,no =', ia,no
#ifdef _DEBUG_WRITE_
                    ATOMLoop2: do ja2 = 1, natm
                       jt=ityp(ja2)
                       if(nint(iatomn(it)) /= nint(iatomn(jt))) cycle

                       pos2(1:3)=pos_t(ja2,1:3)

                       dpos = abs(pos1-pos2)
                       dpos = dpos - floor(dpos+DELTA10)
                       distance = sum(dpos)
                       write(nfout,'(" ja = ",i8," pos = ",3f12.6, " dist = ",d20.8)') &
                            & ja2, pos2(1:3), distance

                       if(kimg==1 .and. iwei(ja2)==2) then
                          pos2(1:3)=-pos_t(ja2,1:3)
                          dpos = abs(pos1-pos2)
                          dpos = dpos - floor(dpos+DELTA10)
                          distance = sum(dpos)
                          write(nfout,'(" ja = ",i8," pos = ",3f12.6, " dist = ",d20.8)') &
                               & ja2, pos2(1:3), distance
                       end if
                    end do ATOMLoop2
#endif

                    call phase_error_with_msg(nfout,'Error in set_ia2ia_symmtry_op !',__LINE__,__FILE__)
                 end if

              end do ATOMLoop
           end do OPLoop

        end do
        deallocate(pos_t)
        deallocate(op_pr)
        return
    end subroutine set_ia2ia_symmtry_op
#endif
!===============================================================================

    subroutine set_ia2ia_symmtry_op_inv(natm,nopr,ia2ia,ia2ia_inv)
        integer,intent(in) :: natm,nopr
        integer,intent(in) :: ia2ia(natm,nopr+af)
        integer,intent(out):: ia2ia_inv(natm,nopr+af)

        integer:: ia,no,ja

        do ia=1,natm
            do no=1,nopr+af
                ja=ia2ia(ia,no)
                if(ja > 0) then
                    ia2ia_inv(ja,no)=ia
                else
                    ia2ia_inv(-ja,no)=-ia
                end if
            end do
        end do
!do ia=1,natm
!do no=1,nopr
!print *,ia,no,ia2ia_inv(ia,no)
!end do
!end do
!stop
        return
    end subroutine set_ia2ia_symmtry_op_inv

   !=====================================
    subroutine m_PP_cnstrct_crotylm_paw
   !=====================================
        real(DP),allocatable :: crotylm(:,:,:)
        integer,allocatable :: iylm(:,:,:)
        integer,allocatable :: nylm(:,:)
        real(DP),allocatable :: opr(:,:,:)

        integer:: l1max,mmax,nsph,n
        integer:: ia,iopr

        call m_PP_find_maximum_l(n)    ! n-1: maximum l
        n = (n-1) + (n-1) + 1
        l1max=n
        mmax=2*l1max-1
        nsph=l1max**2

!       allocate(crotylm(mmax,nsph,nopr))
!       allocate(iylm(mmax,nsph,nopr))
!       allocate(nylm(nsph,nopr))
!       allocate(opr(3,3,nopr))

!       allocate(crotylm_paw(mmax,nsph,nopr,natm))
!       allocate(iylm_paw(mmax,nsph,nopr,natm))
!       allocate(nylm_paw(nsph,nopr,natm))

        allocate(crotylm(mmax,nsph,nopr+af)) !ASMS
        allocate(iylm(mmax,nsph,nopr+af))    !ASMS
        allocate(nylm(nsph,nopr+af))         !ASMS
        allocate(opr(3,3,nopr+af))           !ASMS

        allocate(crotylm_paw(mmax,nsph,nopr+af,natm)) !ASMS
        allocate(iylm_paw(mmax,nsph,nopr+af,natm))    !ASMS
        allocate(nylm_paw(nsph,nopr+af,natm))         !ASMS

        do ia=1,natm
!           do iopr=1,nopr
            do iopr=1,nopr+af      !ASMS
                if(ia2ia_symmtry_op(ia,iopr).gt.0) then
                    opr(:,:,iopr)=op(:,:,iopr)
                else
                    opr(:,:,iopr)=-op(:,:,iopr)
                end if
            end do
!           call get_crotylm(l1max,mmax,nsph,nopr,crotylm,iylm,nylm,opr)
            call get_crotylm(l1max,mmax,nsph,nopr+af,crotylm,iylm,nylm,opr) !ASMS

            crotylm_paw(:,:,:,ia)=crotylm(:,:,:)
            iylm_paw(:,:,:,ia)=iylm(:,:,:)
            nylm_paw(:,:,ia)=nylm(:,:)
        end do

        deallocate(crotylm,iylm,nylm,opr)

        return

   end subroutine m_PP_cnstrct_crotylm_paw

   !===================================
    subroutine m_PP_alloc_crotylm_paw
   !===================================
        integer:: l1max,mmax,nsph,n

        call m_PP_find_maximum_l(n)    ! n-1: maximum l
        n = (n-1) + (n-1) + 1
        l1max=n
        mmax=2*l1max-1
        nsph=l1max**2

        allocate(crotylm_paw(mmax,nsph,nopr,natm))
        allocate(iylm_paw(mmax,nsph,nopr,natm))
        allocate(nylm_paw(nsph,nopr,natm))

        return
    end subroutine m_PP_alloc_crotylm_paw

!===============================================================================
  subroutine m_PP_gfqwei_3D(nfout)
    ! --------------------------------------------------------
    !  Revised by T. Yamasaki (FUJITSU LABORATORIES LTD.)
    !                      7th Jun. 2003
    !  [ if(fq < DELTA) cycle ] -> [ if(abs(fq) < DELTA) cycle ]
    ! Deficit charges have been assumed not to have negative values.
    ! However, deficit charges may have negative values in some conditions.
    ! So, the decision of [ if(fq <DELTA) ] is not appropriate.
    ! --------------------------------------------------------
    integer, intent(in) :: nfout
    integer       :: i, ia, it, lmt1, lmt2
    real(kind=DP) :: fq
    nac = 0
    modnrm = SKIP
    Loop_ntyp: do i = 1, ntyp
       if(m_PP_include_vanderbilt_pot(i) == YES) then
          modnrm = EXECUT
          exit
       end if
    end do Loop_ntyp

    if(ipripp >= 1) then
       if(modnrm == SKIP) then
          write(nfout,'(" !PP modnrm = ",i5," (SKIP)")') modnrm
       else if(modnrm == EXECUT) then
          write(nfout,'(" !PP modnrm = ",i5," (EXECUT)")') modnrm
       endif
    end if

    if(modnrm == SKIP) then
       fqwei = 0.d0; nlmta1 = 0; nlmta2 = 0; nac2ia = 0
    else if(modnrm == EXECUT) then
       do ia = 1, natm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it)) cycle
                fq = q(lmt1,lmt2,it)*iwei(ia)
                if(abs(fq) < DELTA) cycle
                nac = nac + 1
                nlmta1(nac) = lmta(lmt1,ia)
                nlmta2(nac) = lmta(lmt2,ia)
                fqwei(nac)  = fq
                nac2ia(nac) = ia
             end do
          end do
       end do
!PARA3D
       nac_p = 0
       do ia = ista_fs_atm, iend_fs_atm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it)) cycle
                fq = q(lmt1,lmt2,it)*iwei(ia)
                if(abs(fq) < DELTA) cycle
                nac_p = nac_p + 1
             end do
          end do
       end do
       if(.not.allocated(nlmta1_p)) allocate(nlmta1_p(nac_p))
       if(.not.allocated(nlmta2_p)) allocate(nlmta2_p(nac_p))
       if(.not.allocated(fqwei_p)) allocate(fqwei_p(nac_p))
       i = 0
       do ia = ista_fs_atm, iend_fs_atm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it)) cycle
                fq = q(lmt1,lmt2,it)*iwei(ia)
                if(abs(fq) < DELTA) cycle
                i = i + 1
                nlmta1_p(i) = lmta(lmt1,ia)-ista_fs+1
                nlmta2_p(i) = lmta(lmt2,ia)-ista_fs+1
                fqwei_p(i)  = fq
             end do
          end do
       end do
!PARA3D
       nac_p = 0
       do ia = ista_fs_atm, iend_fs_atm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it)) cycle
                fq = q(lmt1,lmt2,it)*iwei(ia)
                if(abs(fq) < DELTA) cycle
                nac_p = nac_p + 1
             end do
          end do
       end do
       if(.not.allocated(nlmta1_p)) allocate(nlmta1_p(nac_p))
       if(.not.allocated(nlmta2_p)) allocate(nlmta2_p(nac_p))
       if(.not.allocated(fqwei_p)) allocate(fqwei_p(nac_p))
       i = 0
       do ia = ista_fs_atm, iend_fs_atm
          it = ityp(ia)
          do lmt1 = 1, ilmt(it)
             do lmt2 = 1, ilmt(it)
                if(ltp(lmt1,it) /= ltp(lmt2,it) .or. mtp(lmt1,it) /= mtp(lmt2,it)) cycle
                fq = q(lmt1,lmt2,it)*iwei(ia)
                if(abs(fq) < DELTA) cycle
                i = i + 1
                nlmta1_p(i) = lmta(lmt1,ia)-ista_fs+1
                nlmta2_p(i) = lmta(lmt2,ia)-ista_fs+1
                fqwei_p(i)  = fq
             end do
          end do
       end do

    end if
#ifndef FQWEI_CHECK
    if(ipripp >= 2) then
#endif
       write(nfout,'(" !PP nac = ", i5)') nac
       write(nfout,'(" !PP --- iwei ---")')
       write(nfout,'(" !PP ",20i4)') (iwei(ia),ia=1,natm)
       write(nfout,'(" !PP --- m_PP_gfqwei ---")')
       do it = 1, nac
          write(nfout,'(" !PP (nlmta1,nlmta2,nac2ia,fqwei)(",i8,") = ",3i8,d16.8)') &
               & it,nlmta1(it),nlmta2(it),nac2ia(it),fqwei(it)
       end do

       write(nfout,'(" !PP nac_p = ", i5)') nac_p
       write(nfout,'(" !PP ista_fs_atm, iend_fs_atm = ",2i5)') ista_fs_atm, iend_fs_atm
       write(nfout,'(" !PP nlmt, ntau, nel_fs_atm = ",3i5)') nlmt,ntau,np_fs_atm
       write(nfout,'(" !PP --- m_PP_gfqwei ---")')
       do it = 1, nac_p
          write(nfout,'(" !PP (nlmta1_p,nlmta2_p,fqwei_p)(",i8,") = ",2i8,d16.8)') &
               & it,nlmta1_p(it),nlmta2_p(it),fqwei_p(it)
       end do
!PARA3D
       write(nfout,'(" !PP nac_p = ", i5)') nac_p
       write(nfout,'(" !PP ista_fs_atm, iend_fs_atm = ",2i5)') ista_fs_atm, iend_fs_atm
       write(nfout,'(" !PP nlmt, ntau, nel_fs_atm = ",3i5)') nlmt,ntau,np_fs_atm
       write(nfout,'(" !PP --- m_PP_gfqwei ---")')
       do it = 1, nac_p
          write(nfout,'(" !PP (nlmta1_p,nlmta2_p,fqwei_p)(",i8,") = ",2i8,d16.8)') &
               & it,nlmta1_p(it),nlmta2_p(it),fqwei_p(it)
       end do
#ifndef FQWEI_CHECK
    end if
#endif


  end subroutine m_PP_gfqwei_3D

!===============================================================================

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
! ==================================================================================================
!===============================================================================
  subroutine m_PP_partial_core_CD_3D(nfout,it,gr_l,paramset)

    integer, intent(in)                        :: nfout,it
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    logical, intent(in)                        :: paramset

    integer,save              :: mmp
    integer                   :: id_sname = -1
    call tstatc0_begin('m_PP_partial_core_CD ',id_sname,1)

    if(it == 1) mmp = 0
    if(itpcc(it) == 1) then
       if(ipripp >= 2) write(nfout,'(" !PP itpcc(",i5," ) = ",i5)') it,itpcc(it)
       mmp = mmp + 1
       itpcc(it) = mmp
       ntpcc     = mmp
       if(.not.paramset) then
          call epseud_pc()
          if(istress == 0) then
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
             call pcc(gr_l)
#else
             call pcc()
#endif
          elseif(istress == 1) then
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
             call pcc_diff(gr_l)
#else
             call pcc_diff()
#endif
          else
             if(ipripp >= 1) write(nfout,'(" !PP pcc - istress = ",i5)') istress
             call phase_error_with_msg(nfout,' !PP istress /= 0 or 1',__LINE__,__FILE__)
          endif
       end if
    endif

    call tstatc0_end(id_sname)
  contains
    subroutine epseud_pc
      real(kind=DP)            :: zeta=0.d0, chgir,rs,vxc1,vxc2,exc1,eepc
      integer                  :: jsm = 1 ! 1:non-mag, 2: magnetic
      real(kind=DP), parameter :: delta40 = 1.d-40
      integer                  :: i, ncut

      if(it == 1) epc = 0.d0
      eepc = 0.d0
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      epc_per_atom(it) = 0.d0
! -----  ascat ceases modifying  -----
      endif
#endif

      if(ipripp >= 2) write(nfout,'(" !PP -- epseud_pc -- ")')

      if(what_is_ps_xctype() == GGA) then ! -(m_PseudoPotential)
            call gga_grad_alloc()
!                ___________________
            call get_gradient_of_rho(jsm,mddiff,mmesh,nmesh(it),rhpcr,radr&
                 &,h(it),fdiff,coeff1,coeff2,rho,grdnts) ! -(b_PseudoPot..)
#ifdef GGA_CHECK
            if(ipripp >= 2) then
               write(nfout,'(" !PP after get_gradient_of_rho ")')
               write(nfout,'(" !PP rad, rho, grdnts(*,1), grdnts(*,2)")')
               do i = 1, mmesh
                  write(nfout,'(" !PP ( ",i4,") rad, rho, grad, grad2 = ",4d20.8," (get_gradient)")') &
                       &        i, radr(i), rho(i), grdnts(i,1), grdnts(i,2)
               end do
            end if
#endif
            if(         ps_xctype == 'ldapw91' .or. ps_xctype == 'LDAPW91' &
                 & .or. ps_xctype == 'ldapbek' .or. ps_xctype == 'LDAPBEK' ) grdnts = 0.d0
!                _______________________
            call gga_xc_potential_atomic(jsm,len_ps_xctype,ps_xctype &
                 &,mmesh,nmesh(it),rho,radr,grdnts,exc,vxc,eps_chg,eps_grad) ! -(b_PseudoPot..)
#ifdef GGA_ATOMIC_WITH_NEW_GNCPP
!!$ASASASASAS
            ncut = 1
!!$ASASASASAS
            do i = mmesh,1,-1
               if(rhpcr(i)/(PAI4*radr(i)**2) > eps_chg) then
                  ncut = i
                  goto 1001
               end if
            end do
1001        continue
            do i = ncut+1, mmesh
               vxc(i) = 0.d0
            end do
#endif
#ifdef GGA_CHECK
            if(ipripp >= 2) then
               write(nfout,'(" !PP ps_xctype, xctype (epseud_pc) = ",a7,2x,a7)') ps_xctype, xctype
               do i = 1, nmesh(it)
                  write(nfout,'(" !PP exc(",i4,") = ",d18.10, " radr(",i4,") = ",d18.10)') i, exc(i),i, radr(i)
               end do
            end if
#endif
            do i = 1, nmesh(it)
               eepc = eepc + wos(i)*exc(i)*PAI4*radr(i)**2
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               epc_per_atom(it) = epc_per_atom(it) + wos(i)*exc(i)*PAI4*radr(i)**2
! -----  ascat ceases modifying  -----
               endif
#endif
            end do

            call gga_grad_dealloc()
      else
         do i = 1, nmesh(it)
            chgir = max(delta40, rhpcr(i))
            rs    = (radr(i)**2 * 3 / chgir)**(1.d0/3.d0)
            call xcpot(ps_xctype,vxc1,vxc2,rs,zeta,exc1)
!!$            wkx(i) = exc1*0.5d0
            eepc =  eepc + wos(i)*rhpcr(i)*exc1*0.5
#ifdef __EDA__
            if(sw_eda==ON) then
! -----  ascat starts modifying  -----
            epc_per_atom(it) = epc_per_atom(it) + wos(i)*rhpcr(i)*exc1*0.5
! -----  ascat ceases modifying  -----
            endif
#endif
         enddo
      end if
      epc = epc + eepc*iatom(it)
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      epc_per_atom(it) = epc_per_atom(it) * iatom(it)
! -----  ascat ceases modifying  -----
      endif
#endif

      if(ipripp >= 1) write(nfout,400) epc
400   format( ' !PP ',' EPC     ENERGY = ',F20.12,' (HARTREE)'  )
    end subroutine epseud_pc

#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
    subroutine pcc(gr_l)
#else
    subroutine pcc
#endif
      real(kind=DP)            :: s, gabs, fac
      real(kind=DP), parameter :: delta10 = 1.d-5
      integer                  :: n, i
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
      real(kind=DP),             dimension(mmesh) :: wkx, wky
#endif

      if(it == 1) schgpc = 0.d0

      s   = 0.d0
      do n = 1, nmesh(it)
         s   =  s + rhpcr(n)*wos(n)
      enddo
      if(ipripp >= 2) write(nfout,'(" !PP chgpc = ",d20.8)') chgpc(it)
      if(dabs(s-chgpc(it)) >= delta10) then
         if(ipripp >=2) write (nfout,'(" !PP partial core charge  s-chgpc= ",d20.8)') s-chgpc(it)
        call phase_error_with_msg(nfout,' |s-chgpc(it)| >= delta0 <<pcc>>',__LINE__,__FILE__)
      end if

! NEC no check
!!CDIR PARALLEL DO PRIVATE ( gabs, n, wkx, wky, fac )
      do i = ista_kngp, iend_kngp  !for mpi
         gabs = gr_l(i)
         rhpcg_l(i,mmp) = 0.d0
         do n = 1, nmesh(it)
            wkx(n) = gabs * radr(n)
         enddo

         call dsjnv(0,nmesh(it),wkx, wky)

         do  n = 1, nmesh(it)
            fac = wos(n)*rhpcr(n)/univol
            rhpcg_l(i,mmp) = rhpcg_l(i,mmp) + fac*wky(n)
         enddo
      end do

      schgpc = schgpc+chgpc(it)*iatom(it)
    end subroutine pcc

#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
    subroutine pcc_diff(gr_l)
#else
    subroutine pcc_diff
#endif
      real(kind=DP)            :: s, gabs, fac
      real(kind=DP), parameter :: delta10 = 1.d-5
      integer                  :: n, i
      logical                  :: test
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
      real(kind=DP),             dimension(mmesh) :: wkx, wky, wkz1
#endif

      if(it == 1) schgpc = 0.d0

      s   = 0.d0
      do n = 1, nmesh(it)
         s   =  s + rhpcr(n)*wos(n)
      enddo
      if(ipripp >= 2) write(nfout,'(" !PP chgpc = ",d20.8)') chgpc(it)
      if(dabs(s-chgpc(it)) >= delta10) then
         if(ipripp >= 2) write (nfout,'(" !PP partial core charge  s-chgpc= ",d20.8)') s-chgpc(it)
        call phase_error_with_msg(nfout,' |s-chgpc(it)| >= delta0 <<pcc_diff>>',__LINE__,__FILE__)
      end if

! NEC no check
!!CDIR PARALLEL DO PRIVATE ( gabs, n, wkx, wky, fac )
      do i = ista_kngp, iend_kngp  !for mpi
         gabs = gr_l(i)
         rhpcg_l(i,mmp) = 0.d0
         rhpcg_diff_l(i,mmp) = 0.d0
         do n = 1, nmesh(it)
            wkx(n) = gabs * radr(n)
         enddo

         call dsjnv(0,nmesh(it),wkx, wky)
         call dsjnv(1,nmesh(it),wkx, wkz1)

         do  n = 1, nmesh(it)
            fac = wos(n)*rhpcr(n)/univol
            rhpcg_l(i,mmp) = rhpcg_l(i,mmp) + fac*wky(n)
            rhpcg_diff_l(i,mmp) = rhpcg_diff_l(i,mmp) &
              & - fac*radr(n)*wkz1(n)
         enddo
      end do

      schgpc = schgpc+chgpc(it)*iatom(it)
      test = .false.
      if(ipripp >= 2 .and. test .and. mype==0) then
         write(nfout,'(" !PP mmp = ",i5)') mmp
         write(nfout,'(" !PP RHPCG")')
         write(nfout,'(" !PP ",(8f10.5))') (rhpcg_l(i,mmp),i=1,20)
         write(nfout,'(" !PP RHPCG_DIFF")')
         write(nfout,'(" !PP ",(8f10.5))') (rhpcg_diff_l(i,mmp),i=1,20)
      endif
    end subroutine pcc_diff
  end subroutine m_PP_partial_core_CD_3D

  subroutine m_PP_local_part_3D(it,nfout,gr_l,ngshell,ngshell_range)
    integer, intent(in)                        :: it, nfout, ngshell
    integer, intent(in), dimension(ngshell,2)  :: ngshell_range
    real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
    integer          :: id_sname = -1
                                                  __TIMER_SUB_START(1223)

    call tstatc0_begin('m_PP_local_part ',id_sname,1)

    call add_localPot_to_vlocr
    call add_localPotEnergy_to_etot1
    if(istress == 0) then
       call localPP_in_Gspace(gr_l,ngshell,ngshell_range)
    elseif(istress == 1) then
       call localPP_diff_in_Gspace(gr_l,ngshell,ngshell_range)
    else
       if(ipripp>=1) write(nfout,'(" !PP localPP_in_Gspace - istress = ",i8)') istress
       call phase_error_with_msg(nfout,' istress /= 0 or 1 <<m_PP_local_part>>',__LINE__,__FILE__)
    endif

    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1223)
  contains

    subroutine add_localPot_to_vlocr           ! short-range
      integer       :: i
      real(kind=DP) :: alp1, alp2, erf1, erf2
      real(kind=DP) :: derf

      select case (long_range_pot_type)
      case (0)
         alp1 = alp(1,it);    alp2 = alp(2,it)
         do i = 1, nmesh(it)
            erf1 = derf(dsqrt(alp1)*radr(i))
            erf2 = derf(dsqrt(alp2)*radr(i))
            vlocr(i) = vlocr(i)&
                 & + ival(it)*(cc(1,it)*(erf1-erf2) + erf2)/radr(i)
         enddo
         if(ipripp >= 1) write(nfout,'(" !PP alp1, alp2, acoef = ",3f12.8)') alp1,alp2,cc(1,it)
      case (1)
         do i = 1, nmesh(it)
            vlocr(i) = vlocr(i) +ival(it) /radr(i)
         enddo
      end select
    end subroutine add_localPot_to_vlocr

    subroutine localPP_in_Gspace(gr_l,ngshell,ngshell_range)
      real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
      integer, intent(in) :: ngshell
      integer, intent(in), dimension(ngshell,2) :: ngshell_range

      real(kind=DP) :: c1, c2, b1,b2,gabs,gg,fac1, fac
      integer       :: ist !mpi
      integer       :: i, n, igs, igs_start
!!$      integer       :: igt1
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      real(kind=DP),             dimension(mmesh) :: wkx, wky
#endif
      real(kind=DP) :: psc_shell

      fac1  = PAI4/univol

      select case (long_range_pot_type)
      case (0)
         c1 = cc(1,it)*PAI4*ival(it)/univol
         c2 = cc(2,it)*PAI4*ival(it)/univol
         b1 = -0.25d0/alp(1,it)
         b2 = -0.25d0/alp(2,it)
      case (1)
         c1 = PAI4 *ival(it) /univol
      end select

      igs_start = 1
#ifdef ENABLE_ESM_PACK
      if(sw_esm==OFF)then
         if(myrank_chg==0) then
            psc_l(1,it) = 0.d0
            igs_start = 2
         end if
      else
         psc_l(:,it) = 0.d0
      endif
#else
!      if(mype==0) then
      if(myrank_chg==0) then
         psc_l(1,it) = 0.d0
         igs_start = 2
      end if
#endif
!CDIR PARALLEL DO PRIVATE ( gabs, gg, wkx, wky, n, fac, psc_shell )
      do igs = igs_start, ngshell
         gabs = gr_l(ngshell_range(igs,1))
         gg   = gabs*gabs

         select case (long_range_pot_type)
         case (0)
#ifndef ENABLE_ESM_PACK
            psc_shell = -(c1*exp(b1*gg) + c2*exp(b2*gg))/gg
#else
            psc_shell = 0.d0
            if(sw_esm==OFF) psc_shell = -(c1*exp(b1*gg) + c2*exp(b2*gg))/gg
#endif
         case (1)
            psc_shell = -c1 /gg
         end select

         wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
         call dsjnv(0,nmesh(it),wkx,wky)
         do n = 1, nmesh(it)
            fac = wos(n)*vlocr(n)*radr(n)**2*fac1
#ifdef ENABLE_ESM_PACK
            if(sw_esm==OFF) then
               psc_shell = psc_shell + fac*wky(n)
            else
               if(ngshell_range(igs,1).ne.1)then
                  psc_shell = psc_shell + fac*wky(n)
               else
                  psc_shell = psc_shell + fac
               endif
            endif
#else
            psc_shell = psc_shell + fac*wky(n)
#endif
         enddo
         do n = ngshell_range(igs,1), ngshell_range(igs,2)
            psc_l(n,it) = psc_shell
         end do
      enddo

      if ( sw_add_vlocal_gzero == ON ) then
         if (myrank_chg==0) then
            select case (long_range_pot_type)
            case (0)
               psc_shell = -( c1*b1 +c2*b2 )
            case (1)
               psc_shell = 0.0d0
            end select

            do n = 1, nmesh(it)
               fac = wos(n)*vlocr(n)*radr(n)**2*fac1
               psc_shell = psc_shell + fac
            enddo
            igs = 1;  psc_l(igs,it) = psc_shell
         end if
      endif

!!$      if(mype==0) psc_l(1,it) = 0.d0
!!$      if(ista_kngp==1) then
!!$         ist = 2
!!$      else
!!$         ist = ista_kngp
!!$      endif
!!$      do i = ist, iend_kngp  !for mpi
!!$         gabs = gr_l(i)
!!$         gg   = gabs*gabs
!!$         psc_l(i,it) = -(c1*exp(b1*gg) + c2*exp(b2*gg))/gg
!!$         wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
!!$!!$         do n = 1, nmesh(it)
!!$!!$            if(wkx(n) > 1.d0) exit
!!$!!$         end do
!!$!!$         igt1 = n
!!$!!$         if(i < 100) write(6,'(" i, igt1 = ",2i8)') i, igt1
!!$!!$         do n = 1, nmesh(it)
!!$!!$            wkx(n) = gabs * radr(n)
!!$!!$         end do
!!$         call dsjnv(0,nmesh(it),wkx,wky)
!!$         do n = 1, nmesh(it)
!!$            fac = wos(n)*vlocr(n)*radr(n)**2*fac1
!!$            psc_l(i,it) = psc_l(i,it) + fac*wky(n)
!!$         enddo
!!$      enddo
    end subroutine localPP_in_Gspace

    subroutine localPP_diff_in_Gspace(gr_l,ngshell,ngshell_range)
      real(kind=DP), intent(in), dimension(ista_kngp:iend_kngp) :: gr_l
      integer, intent(in) :: ngshell
      integer, intent(in), dimension(ngshell,2) :: ngshell_range

      integer       :: i, n
      real(kind=DP) :: c1, c2, b1,b2,gabs,gg,fac1, fac
      logical       :: test
      integer :: ist, igs_start, igs !mpi

#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      real(kind=DP),             dimension(mmesh) :: wkx, wky, wkz1
#endif
      real(kind=DP)  :: psc_shell, psc_diff_shell

      fac1  = PAI4/univol

      select case (long_range_pot_type)
      case (0)
         c1 = cc(1,it)*PAI4*ival(it)/univol
         c2 = cc(2,it)*PAI4*ival(it)/univol
         b1 = -0.25d0/alp(1,it)
         b2 = -0.25d0/alp(2,it)
      case (1)
         c1 = PAI4 *ival(it) /univol
      end select

      igs_start = 1
!f      if(mype==0) then
      if(myrank_chg==0) then
         psc_l(1,it) = 0.d0
         psc_diff_l(1,it) = 0.d0
         igs_start = 2
      endif

! NEC no check
!!CDIR PARALLEL DO PRIVATE ( gabs, gg, wkx, wky, wkz1, n, fac, psc_shell, psc_diff_shell )
      do igs = igs_start, ngshell
         gabs = gr_l(ngshell_range(igs,1))
         gg   = gabs*gabs

         select case (long_range_pot_type)
         case (0)
            psc_shell = -(c1*exp(b1*gg) + c2*exp(b2*gg))/gg
            psc_diff_shell = (c1*exp(b1*gg) + c2*exp(b2*gg))/(gg*gg) &
                 &         - (c1*exp(b1*gg)*b1 + c2*exp(b2*gg)*b2) / gg
         case (1)
            psc_shell = -c1 /gg
            psc_diff_shell = c1 /gg**2
         end select

         wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
         call dsjnv(0,nmesh(it),wkx,wky)
         call dsjnv(1,nmesh(it),wkx,wkz1)
         do n = 1, nmesh(it)
            fac = wos(n)*vlocr(n)*radr(n)**2*fac1
            psc_shell = psc_shell + fac*wky(n)
            psc_diff_shell = psc_diff_shell &
                 &         - fac*wkz1(n)*radr(n)/gabs/2.d0
         enddo

         do n = ngshell_range(igs,1), ngshell_range(igs,2)
            psc_l(n,it) = psc_shell
            psc_diff_l(n,it) = psc_diff_shell
         end do
      end do

      if ( sw_add_vlocal_gzero == ON ) then
         if (myrank_chg==0) then
            select case (long_range_pot_type)
            case (0)
               psc_shell = -( c1*b1 +c2*b2 )
               psc_diff_shell = -( c1*b1**2 +c2*b2**2 ) /2.0d0
            case (1)
               psc_shell = 0.0d0
               psc_diff_shell = 0.0d0
            end select

            do n = 1, nmesh(it)
               fac = wos(n)*vlocr(n)*radr(n)**2*fac1
               psc_shell = psc_shell + fac
               psc_diff_shell = psc_diff_shell &
                    &         - fac *radr(n)**2 /6.d0
            enddo
            igs = 1;  psc_l(igs,it) = psc_shell
            psc_diff_l(igs,it) = psc_diff_shell
         end if
      endif

!!$      do i = ist, iend_kngp  !for mpi
!!$         gabs = gr_l(i)
!!$         gg   = gabs*gabs
!!$         psc_l(i,it) = -(c1*exp(b1*gg) + c2*exp(b2*gg))/gg
!!$         psc_diff_l(i,it) = (c1*exp(b1*gg) + c2*exp(b2*gg))/(gg*gg) &
!!$                    &     - (c1*exp(b1*gg)*b1 + c2*exp(b2*gg)*b2) / gg
!!$         wkx(1:nmesh(it)) = gabs*radr(1:nmesh(it))
!!$         call dsjnv(0,nmesh(it),wkx,wky)
!!$         call dsjnv(1,nmesh(it),wkx,wkz1)
!!$         do n = 1, nmesh(it)
!!$            fac = wos(n)*vlocr(n)*radr(n)**2*fac1
!!$            psc_l(i,it) = psc_l(i,it) + fac*wky(n)
!!$            psc_diff_l(i,it) = psc_diff_l(i,it) &
!!$                     &       - fac*wkz1(n)*radr(n)/gabs/2.d0
!!$         enddo
!!$      enddo

      test = .false.
      if(ipripp >= 2 .and. test .and. mype==0) then
         write(nfout,'(" !PP it = ",i5)') it
         write(nfout,'(" !PP PSC = ")')
         write(nfout,'(" !PP ",2(8f10.5))') (psc_l(i,it),i=2,20)
         write(nfout,'(" !PP PSC_DIFF = ")')
         write(nfout,'(" !PP ",2(8f10.5))') (psc_diff_l(i,it),i=2,20)
      endif
    end subroutine localPP_diff_in_Gspace

    subroutine add_localPotEnergy_to_etot1
      integer       :: i
      real(kind=DP) :: s

      if(it == 1) etot1 = 0.d0
      s = 0.d0
      do i = 1, nmesh(it)
         s = s + wos(i)*vlocr(i)*radr(i)**2
      enddo
      if(ipripp >= 2) write(nfout,*) ' !PP Int_{i=1}^{nmesh} wos(i)vlocr(i)*radr(i)**2 = ', s

      if ( sw_add_vlocal_gzero == OFF ) then
         select case (long_range_pot_type)
         case (0)
            etot1 = etot1 + iatom(it)*PAI/univol* &
                 & (4*s + ival(it)*(cc(1,it)/alp(1,it)+cc(2,it)/alp(2,it)))
         case (1)
            etot1 = etot1 + iatom(it)*PAI/univol* 4*s
         end select
      end if
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      PP_per_atom(it) = PP_per_atom(it) &
                    & + iatom(it)*PAI/univol*(4*s + ival(it)*(cc(1,it)/alp(1,it)+cc(2,it)/alp(2,it)))
! -----  ascat ceases modifying  -----
      endif
#endif

      if(ipripp >= 1) then
         write(nfout,'(" !PP univol = ",f20.12)') univol
         write(nfout,'(" !PP etot1  = ",f20.12)') etot1
      end if
#ifdef GGA_CHECK
      if(ipripp >= 2) then
         do i = 1, 10
            write(nfout,'(" !PP vlocr(",i5," ) = ",f12.6)') i, vlocr(i)
         end do
         do i = nmesh(it)-9,nmesh(it)
            write(nfout,'(" !PP vlocr(",i5," ) = ",f12.6)') i, vlocr(i)
         end do
      end if
#endif
    end subroutine add_localPotEnergy_to_etot1
  end subroutine m_PP_local_part_3D
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
! ==================================================================================================
  subroutine m_PP_renew_etot1(univolo)
    real(kind=DP),intent(in) :: univolo
    etot1 = etot1*univol/univolo
  end subroutine m_PP_renew_etot1

  subroutine m_PP_wd_PP_parameters_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kgp)
    use m_Parallelization,     only :mpi_ke_world, mpi_chg_world

    integer, intent(in) :: nfout,nfcntn_bin,kgp
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer             :: ista, iend, i
    integer             :: id_sname = -1
                                                  __TIMER_SUB_START(1365)

    call tstatc0_begin('m_PP_wd_PP_parameters ',id_sname)

    ista = ista_kngp; iend = iend_kngp

    if(F_CNTN_BIN_partitioned) then
       write(nfcntn_bin) ival,ivanl,itau,eps,iloc,lpsmax,itpcc,ilmt,nlmtt&
            &, nlmta,lmta,lmtt,ltp,mtp,taup&
            &, etot1,epc,schgpc,q,dion,isph,il2p,dl2p,modnrm,nac,fqwei&
            &, nlmta1,nlmta2,ps_xctype,iqitg
       call wd_kgp_array_p(ntyp,psc_l)      ! -(contained here) write(nfcntn_bin) psc_l
       call wd_kgp_array_p(nqitg,qitg_l)    ! -(contained here) write(nfcntn_bin) qitg_l
       call wd_kgp_array_p(ntpcc,rhpcg_l)   ! -(contained here) write(nfcntn_bin) rhpcg_l
       if(ipripp >= 2) then
          write(nfout,'(" !PP rhpcg_l is written (m_PP_wd_PP_parameters)")')
          write(nfout,'(" !PP size of psc_l   = ", i9)') (iend-ista+1)*ntyp*8
          write(nfout,'(" !PP size of qitg_l  = ", i9)') (iend-ista+1)*nqitg*8
          write(nfout,'(" !PP size of rhpcg_l = ", i9)') (iend-ista+1)*ntpcc*8
       end if
       if(istress==ON) then
          call wd_kgp_array_p(ntyp,psc_diff_l)   ! -(c.h.) write(nfcntn_bin) psc_l
          call wd_kgp_array_p(nqitg,qitg_diff_l) ! -(c.h.) write(nfcntn_bin) qitg_l
          call wd_kgp_array_p(ntpcc,rhpcg_diff_l)! -(c.h.) write(nfcntn_bin) rhpcg_l
       endif
       if(sw_positron /= OFF .and. corecharge_cntnbin == ON) then
          call wd_kgp_array_p(ntyp,rhcg_l)   ! -(c.h.) write(nfcntn_bin) rhcg_l
          call wd_kgp_array_p(ntyp,rhceg_l)  ! -(c.h.) write(nfcntn_bin) rhceg_l
          call wd_kgp_array_p(ntyp,rhchg_l)  ! -(c.h.) write(nfcntn_bin) rhchg_l
       end if
    else
       if(mype==0) then
          write(nfcntn_bin) ival,ivanl,itau,eps,iloc,lpsmax,itpcc,ilmt,nlmtt&
               &, nlmta,lmta,lmtt,ltp,mtp,taup&
               &, etot1,epc,schgpc,q,dion,isph,il2p,dl2p,modnrm,nac,fqwei&
               &, nlmta1,nlmta2,ps_xctype,iqitg
       end if
       call wd_kgp_array(kgp,ntyp,psc_l)      ! -(contained here) write(nfcntn_bin) psc_l
       call wd_kgp_array(kgp,nqitg,qitg_l)    ! -(contained here) write(nfcntn_bin) qitg_l
       call wd_kgp_array(kgp,ntpcc,rhpcg_l)   ! -(contained here) write(nfcntn_bin) rhpcg_l
       if(ipripp >= 2) then
          write(nfout,'(" !PP rhpcg_l is written (m_PP_wd_PP_parameters)")')
          write(nfout,'(" !PP size of psc_l   = ", i9)') iend_kngp*ntyp*8
          write(nfout,'(" !PP size of qitg_l  = ", i9)') iend_kngp*nqitg*8
          write(nfout,'(" !PP size of rhpcg_l = ", i9)') iend_kngp*ntpcc*8
       end if
       if(istress==ON) then
          call wd_kgp_array(kgp,ntyp,psc_diff_l)   ! -(c.h.) write(nfcntn_bin) psc_l
          call wd_kgp_array(kgp,nqitg,qitg_diff_l) ! -(c.h.) write(nfcntn_bin) qitg_l
          call wd_kgp_array(kgp,ntpcc,rhpcg_diff_l)! -(c.h.) write(nfcntn_bin) rhpcg_l
       endif
       if(sw_positron /= OFF .and. corecharge_cntnbin == ON) then
          call wd_kgp_array(kgp,ntyp,rhcg_l)   ! -(c.h.) write(nfcntn_bin) rhcg_l
          call wd_kgp_array(kgp,ntyp,rhceg_l) ! -(c.h.) write(nfcntn_bin) rhceg_l
          call wd_kgp_array(kgp,ntyp,rhchg_l) ! -(c.h.) write(nfcntn_bin) rhchg_l
          if(ipripp >= 2) then
             write(nfout,'(" !PP rhcg_l, rhceg_l, rhchg_l have been written (m_PP_wd_PP_parameters)")')
             write(nfout,'(" !PP size of rhcg_l   = ", i9)') iend_kngp*ntyp*8
             write(nfout,'(" !PP size of rhceg_l  = ", i9)') iend_kngp*ntyp*8
             write(nfout,'(" !PP size of rhchg_l  = ", i9)') iend_kngp*ntyp*8
          end if

       end if
    end if
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1365)
  contains
    subroutine wd_kgp_array(n1,n2,reduced_array)
      integer, intent(in)                         :: n1, n2
      real(DP),intent(in),dimension(ista:iend,n2) :: reduced_array

      real(DP),allocatable,dimension(:,:)         :: a_mpi,b_mpi,a
      real(DP),allocatable,dimension(:)           :: c1,c2
      integer                                     :: i,j
                                                  __TIMER_SUB_START(1366)

      if(n2 >= 1) then
         if(n1*n2*8 > MAXIMUM_MPI_SIZE) then
            if(ipripp >= 2) write(nfout,'(" n1, n2 = ",i10,",",i8, &
                 & "  ( n1*n2*8 > MAXIMUM_MPI_SIZE = ",i12,") <<wd_kgp_array>>")') &
                 & n1, n2, MAXIMUM_MPI_SIZE
            allocate(a(n1,n2),c1(n1),c2(n1))
!!$ASASASASAS
            a = 0.0d0
            c1 = 0.0d0
!!$ASASASASAS
            do j = 1, n2
               do i = ista, iend
                  c1(i) = reduced_array(i,j)
               end do
!!$ASASASASAS
               c2 = 0.0d0
!!$ASASASASAS
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_chg_world,1388)
               call mpi_allreduce(c1,c2,n1,mpi_double_precision,mpi_sum, mpi_chg_world, ierr)
                                                  __TIMER_IOCOMM_STOP(1388)
               if(mype == 0) a(:,j) = c2(:)
            end do
                                                  __TIMER_IODO_START(1389)
            if(mype == 0) write(nfcntn_bin) a
                                                  __TIMER_IODO_STOP(1389)
            deallocate(c2,c1,a)
         else
            allocate(a_mpi(n1,n2)); a_mpi = 0.d0
            allocate(b_mpi(n1,n2))
                                                  __TIMER_IODO_START(1390)
            do j = 1, n2
               do i = ista, iend
                  a_mpi(i,j) = reduced_array(i,j)
               enddo
            enddo
                                                  __TIMER_IODO_STOP(1390)
                                                  __TIMER_IOCOMM_START_w_BARRIER(mpi_chg_world,1391)
            call mpi_allreduce(a_mpi,b_mpi,n1*n2,mpi_double_precision, mpi_sum,mpi_chg_world,ierr)
                                                  __TIMER_IOCOMM_STOP(1391)
                                                  __TIMER_IODO_START(1392)
            if(mype == 0) write(nfcntn_bin) b_mpi
                                                  __TIMER_IODO_STOP(1392)
            deallocate(a_mpi); deallocate(b_mpi)
         end if
      end if
                                                  __TIMER_SUB_STOP(1366)
    end subroutine wd_kgp_array

    subroutine wd_kgp_array_p(n2,reduced_array)
      integer, intent(in)                         :: n2
      real(DP),intent(in),dimension(ista:iend,n2) :: reduced_array

      real(DP),allocatable,dimension(:,:)         :: a_wk
      integer                                     :: i,i2,j

                                                  __TIMER_SUB_START(1367)
      if(n2 >= 1) then
!!$         write(nfcntn_bin) reduced_array(ista:iend,n2)
         allocate(a_wk(iend-ista+1,n2)); a_wk = 0.d0
                                                  __TIMER_IODO_START(1393)
         do j = 1, n2
            do i = 1, iend-ista+1
               i2 = i + ista - 1
               a_wk(i,j) = reduced_array(i2,j)
            enddo
         enddo
                                                  __TIMER_IODO_STOP(1393)
                                                  __TIMER_IODO_START(1394)
         write(nfcntn_bin) a_wk
                                                  __TIMER_IODO_STOP(1394)
         deallocate(a_wk)
      end if
                                                  __TIMER_SUB_STOP(1367)
    end subroutine wd_kgp_array_p

  end subroutine m_PP_wd_PP_parameters_3D

  subroutine m_PP_rd_PP_parameters_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned, kgp)
    use m_Parallelization,     only :mpi_ke_world

    integer, intent(in) :: nfout,nfcntn_bin,kgp
    logical, intent(in) :: F_CNTN_BIN_partitioned
    integer             :: ista, iend, it, i, lmt1, lmt2
    integer             :: id_sname = -1
                                                  __TIMER_SUB_START(1361)

    call tstatc0_begin('m_PP_rd_PP_parameters ',id_sname)

    ista = ista_kngp; iend = iend_kngp
    if(F_CNTN_BIN_partitioned) then
       read(nfcntn_bin) ival,ivanl,itau,eps,iloc,lpsmax,itpcc,ilmt,nlmtt&
            &, nlmta,lmta,lmtt,ltp,mtp,taup&
            &, etot1,epc,schgpc,q,dion,isph,il2p,dl2p,modnrm,nac,fqwei&
            &, nlmta1,nlmta2,ps_xctype,iqitg
       call rd_kgp_array_p(ntyp,psc_l)    !-(c.h.) read(nfcntn_bin) psc_l
       call rd_kgp_array_p(nqitg,qitg_l)  !-(c.h.) read(nfcntn_bin) qitg_l
       call rd_kgp_array_p(ntpcc,rhpcg_l) !-(c.h.) read(nfcntn_bin) rhpcg_l
       if(ipripp >= 2) then
          write(nfout,'(" !PP rhpcg_l is read (m_PP_rd_PP_parameters)")')
          write(nfout,'(" !PP size of psc_l   = ", i9)') (iend-ista+1)*ntyp*8
          write(nfout,'(" !PP size of qitg_l  = ", i9)') (iend-ista+1)*nqitg*8
          write(nfout,'(" !PP size of rhpcg_l = ", i9)') (iend-ista+1)*ntpcc*8
       end if

       if(istress == ON) then
          call rd_kgp_array_p(ntyp,psc_diff_l)    !-(c.h.) read(nfcntn_bin) psc_diff_l
          call rd_kgp_array_p(nqitg,qitg_diff_l)  !-(c.h.) read(nfcntn_bin) qitg_diff_l
          call rd_kgp_array_p(ntpcc,rhpcg_diff_l) !-(c.h.) read(nfcntn_bin) rhpcg_diff_l
       endif

       if(corecharge_cntnbin == ON) then
          if(sw_positron /= OFF) then
             call rd_kgp_array_p(ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             call rd_kgp_array_p(ntyp,rhceg_l)   !-(c.h.) read(nfcntn_bin) rhceg_l
             call rd_kgp_array_p(ntyp,rhchg_l)   !-(c.h.) read(nfcntn_bin) rhchg_l
          else
             allocate(rhcg_l(ista_kngp:iend_kngp,ntyp))
             call rd_kgp_array_p(ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             call rd_kgp_array_p(ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhceg_l
             call rd_kgp_array_p(ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhchg_l
             deallocate(rhcg_l)
          end if
       end if

    else
       if(mype==0) then
          read(nfcntn_bin) ival,ivanl,itau,eps,iloc,lpsmax,itpcc,ilmt,nlmtt&
               &, nlmta,lmta,lmtt,ltp,mtp,taup&
               &, etot1,epc,schgpc,q,dion,isph,il2p,dl2p,modnrm,nac,fqwei&
               &, nlmta1,nlmta2,ps_xctype,iqitg
       end if
       call bcast_nfcntn_bin                  !-(c.h.)

       call rd_kgp_array(kgp,ntyp,psc_l)    !-(c.h.) read(nfcntn_bin) psc_l
       call rd_kgp_array(kgp,nqitg,qitg_l)  !-(c.h.) read(nfcntn_bin) qitg_l
       call rd_kgp_array(kgp,ntpcc,rhpcg_l) !-(c.h.) read(nfcntn_bin) rhpcg_l
       if(ipripp >= 2) then
          write(nfout,'(" !PP psc_l, qitg_l and rhpcg_l have been read (m_PP_rd_PP_parameters)")')
          write(nfout,'(" !PP size of psc_l   = ", i9)') iend_kngp*ntyp*8
          write(nfout,'(" !PP size of qitg_l  = ", i9)') iend_kngp*nqitg*8
          write(nfout,'(" !PP size of rhpcg_l = ", i9)') iend_kngp*ntpcc*8
       end if

       if(istress == ON) then
          call rd_kgp_array(kgp,ntyp,psc_diff_l)    !-(c.h.) read(nfcntn_bin) psc_diff_l
          call rd_kgp_array(kgp,nqitg,qitg_diff_l)  !-(c.h.) read(nfcntn_bin) qitg_diff_l
          call rd_kgp_array(kgp,ntpcc,rhpcg_diff_l) !-(c.h.) read(nfcntn_bin) rhpcg_diff_l
       endif


       if(ipripp >= 2) write(nfout,'(" !PP corecharge_cntnbin = ",i8)') corecharge_cntnbin
       if(corecharge_cntnbin == ON) then
          if(sw_positron /= OFF) then
             call rd_kgp_array(kgp,ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             call rd_kgp_array(kgp,ntyp,rhceg_l)  !-(c.h.) read(nfcntn_bin) rhceg_l
             call rd_kgp_array(kgp,ntyp,rhchg_l)  !-(c.h.) read(nfcntn_bin) rhchg_l
             if(ipripp >= 2) then
                write(nfout,'(" !PP rhcg_l, rhceg_l, and rhchg_l have been read (m_PP_rd_PP_parameters)")')
                write(nfout,'(" !PP size of rhcg_l   = ", i9)') iend_kngp*ntyp*8
                write(nfout,'(" !PP size of rhceg_l  = ", i9)') iend_kngp*ntyp*8
                write(nfout,'(" !PP size of rhchg_l  = ", i9)') iend_kngp*ntyp*8
             end if
          else
             allocate(rhcg_l(ista_kngp:iend_kngp,ntyp))
             call rd_kgp_array(kgp,ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             call rd_kgp_array(kgp,ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             call rd_kgp_array(kgp,ntyp,rhcg_l)    !-(c.h.) read(nfcntn_bin) rhcg_l
             deallocate(rhcg_l)
          end if
       end if
    end if

    do it = 1, ntyp
       if_q_isnotzero(:,:,it) = NO
       call make_index_arrays(nfout,it,.false.)
       do lmt1 = 1, ilmt(it)
          do lmt2 = 1, ilmt(it)
             if(dabs(q(lmt1,lmt2,it)) > DELTA) if_q_isnotzero(lmt1,lmt2,it) = YES
          end do
       end do
    end do
    call tstatc0_end(id_sname)
                                                  __TIMER_SUB_STOP(1361)
  contains
    subroutine rd_kgp_array(n1,n2,reduced_array)
      integer, intent(in)                          :: n1, n2
      real(DP),intent(out),dimension(ista:iend,n2) :: reduced_array
      real(DP),allocatable,dimension(:,:) :: a_mpi
      real(DP),allocatable,dimension(:) :: b_mpi
      integer                         :: i, j

                                                  __TIMER_SUB_START(1362)
      if(n2 >= 1) then
         if(n1*n2*8 > MAXIMUM_BCAST_SIZE) then
            if(ipripp >= 2) write(nfout,'(" n1, n2 = ",i10,",",i8, &
                 & "  ( n1*n2*8 > MAXIMUM_BCAST_SIZE = ",i12,") <<rd_kgp_array>>")') &
                 & n1, n2, MAXIMUM_BCAST_SIZE
            allocate(a_mpi(n1,n2),b_mpi(n1))
                                                  __TIMER_IODO_START(1380)
            if(mype == 0) read(nfcntn_bin) a_mpi
                                                  __TIMER_IODO_STOP(1380)
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1381)
            do j = 1, n2
               b_mpi(:) = a_mpi(:,j)
               call mpi_bcast(b_mpi,n1,mpi_double_precision,0,MPI_CommGroup,ierr)
               reduced_array(ista:iend,j) = b_mpi(ista:iend)
            end do
                                                  __TIMER_IOCOMM_STOP(1381)
            deallocate(b_mpi,a_mpi)
         else
            allocate(a_mpi(n1,n2))
                                                  __TIMER_IODO_START(1382)
            if(mype == 0) read(nfcntn_bin) a_mpi
                                                  __TIMER_IODO_STOP(1382)
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1383)
            call mpi_bcast(a_mpi,n1*n2,mpi_double_precision,0,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1383)
                                                  __TIMER_IODO_START(1384)
            do j = 1, n2
               do i = ista, iend
                  reduced_array(i,j) = a_mpi(i,j)
               end do
            end do
                                                  __TIMER_IODO_STOP(1384)
            deallocate(a_mpi)
         end if
      end if

                                                  __TIMER_SUB_STOP(1362)
    end subroutine rd_kgp_array

    subroutine rd_kgp_array_p(n2,reduced_array)
      integer, intent(in)                          :: n2
      real(DP),intent(out),dimension(ista:iend,n2) :: reduced_array
      real(DP),allocatable,dimension(:,:)          :: a_wk
      integer                                      :: i, j

                                                  __TIMER_SUB_START(1363)
      if(n2 >= 1) then
         allocate(a_wk(iend-ista+1,n2))
                                                  __TIMER_IODO_START(1385)
         read(nfcntn_bin) a_wk
                                                  __TIMER_IODO_STOP(1385)
                                                  __TIMER_IODO_START(1386)
         do j = 1, n2
            do i = ista, iend
               reduced_array(i,j) = a_wk(i-ista+1,j)
            end do
         end do
                                                  __TIMER_IODO_STOP(1386)
         deallocate(a_wk)
      end if
                                                  __TIMER_SUB_STOP(1363)

    end subroutine rd_kgp_array_p

    subroutine bcast_nfcntn_bin
                                                  __TIMER_SUB_START(1364)
                                                  __TIMER_IOCOMM_START_w_BARRIER(MPI_CommGroup,1387)
      call mpi_bcast(ival,ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(ivanl,nloc*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(itau,nloc*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(eps,nloc*ntau*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(iloc,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(lpsmax,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(itpcc,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(ilmt,ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(nlmtt,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(nlmta,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(lmta,nlmt*natm,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(lmtt,nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(ltp,nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(mtp,nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(taup,nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(etot1,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(epc,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(schgpc,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(q,nlmt*nlmt*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(dion,nlmt*nlmt*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(isph,nlmt*nlmt*6*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(il2p,nlmt*nlmt*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(dl2p,nlmt*nlmt*6*ntyp,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(modnrm,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(nac,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(fqwei,nlmt*ntau*natm,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(nlmta1,nlmt*ntau*natm,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(nlmta2,nlmt*ntau*natm,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(ps_xctype,len_ps_xctype,mpi_character,0,MPI_CommGroup,ierr)
      call mpi_bcast(iqitg,nloc*ntau*nloc*ntau*nloc*2*ntyp,mpi_integer,0,MPI_CommGroup,ierr)
                                                  __TIMER_IOCOMM_STOP(1387)
                                                  __TIMER_SUB_STOP(1364)
    end subroutine bcast_nfcntn_bin
  end subroutine m_PP_rd_PP_parameters_3D


!===============================================================================

! ====================== added by K. Tagami ==================== 5.0
  subroutine m_PP_alloc_aug_charge_mixing
    deallocate(ia2ia_symmtry_op)
    deallocate(ia2ia_symmtry_op_inv)
    deallocate(crotylm_paw)
    deallocate(iylm_paw)
    deallocate(nylm_paw)
    return
  end subroutine m_PP_alloc_aug_charge_mixing

  subroutine m_PP_dealloc_aug_charge_mixing
    deallocate(ia2ia_symmtry_op)
    deallocate(ia2ia_symmtry_op_inv)
    deallocate(crotylm_paw)
    deallocate(iylm_paw)
    deallocate(nylm_paw)
    return
  end subroutine m_PP_dealloc_aug_charge_mixing
! ======================================================================= 5.0

! =============================== added by K. Tagami ================= 11.0

    integer function get_lp1_in_AngMomList( has_soc, inum )
      integer, intent(in) :: inum                  ! lp1 means  l plus 1
      logical, intent(in) :: has_soc

      integer :: itmp1
! ---------
!                           in the case of with soc:
!                                inum |   1   2   3   4   5   6  7
!                                  l  |   0   1   1   2   2   3  3
!                                l+1  |   1   2   2   3   3   4  4
! --
      integer :: data_lp1(7)=(/ 1, 2, 2, 3, 3, 4, 4 /)
! --
      if ( has_soc ) then
         get_lp1_in_AngMomList = data_lp1( inum )
      else
         get_lp1_in_AngMomList = inum
      endif

    end function get_lp1_in_AngMomList

    integer function get_jph_in_AngMomList( has_soc, inum )
      integer, intent(in) :: inum                      ! jph = j plus half
      logical, intent(in) :: has_soc

      integer :: itmp1
! ---------
!                           in the case of with soc:
!                                inum |  1    2    3    4    5    6    7
!                                  j  |  1/2  1/2  3/2  3/2  5/2  5/2  7/2
!                              j +1/2 |  1    1    2    2    3    3    4
!
!
      integer :: data_jph(7)=(/ 1, 1, 2, 2, 3, 3, 4 /)
! -----
      if ( has_soc ) then
         get_jph_in_AngMomList = data_jph( inum )
      else
         get_jph_in_AngMomList = inum              ! same value with l
      endif

    end function get_jph_in_AngMomList
! ================================================================= 11.0

    subroutine m_PP_alloc_radr()
       if(.not.allocated(radr)) allocate(radr(mmesh));radr=0.d0
    end subroutine m_PP_alloc_radr

    subroutine m_PP_dealloc_radr()
       if(allocated(radr)) deallocate(radr)
    end subroutine m_PP_dealloc_radr

    subroutine m_PP_rd_ival(nfp,it,nfout,ival)
      integer, intent(in) :: nfp,it,nfout
      real(kind=DP), intent(out) :: ival
      call mpi_barrier(MPI_CommGroup,ierr)
      call read_natomn_ival()

      contains
        subroutine read_natomn_ival()
          integer :: natomn, iloc,itpcc
          if(mype == 0) then
             rewind(nfp)
             if(ipripp >= 2) write(nfout,'(" !PP  READING POTENTIAL FILE ",i5)') nfp
             call skip_commentlines_pp(nfp)
             read(nfp,*)      natomn,ival,iloc,itpcc
             if(ipripp >= 2) write(nfout,110) natomn,ival,iloc,itpcc
110          FORMAT(' !PP ',2f8.4,2I4,'  : NATOMN, IVAL, ILOC, ITPCC ')
             if( abs(iatomn(it) - natomn) > 1.d-8)  then
                if(ipripp >= 2) write(nfout,*) ' !PP iatomn.ne.natomn ',iatomn(it),natomn
                if(iatomn(it)>=1) call phase_execution_error(INVALID_ATOMIC_NUMBER)
                write(nfout,'(a,i4,a,f7.3)') ' !PP atomic number for element',it,':',natomn
             end if
          end if
        end subroutine read_natomn_ival

        subroutine skip_commentlines_pp(nf)
          integer, intent(in) :: nf
          integer, parameter  ::    len_str = 80
          character(len=len_str) :: str
          logical :: comment_statement
          comment_statement = .true.
          do while(comment_statement)
             read(nfp,'(a80)',end=1000) str
             if(str(1:1) == '#'.or. str(1:1) == '$' .or. str(1:1) == '!' &
                  & .or. str(1:1) == '%' .or. str(1:1) == '*') then
                if(ipripp >= 1) write(nfout,'(a80)') str
             else if(len(trim(str)) == 0) then
                if(ipripp >= 1) write(nfout,'(a80)') str
             else
                comment_statement = .false.
             endif
          enddo
          backspace(nfp)
          return
1000      write(nfout,'(" eof is reached")')
          call phase_error_with_msg(nfout,' eof is reached in PP reading',__LINE__,__FILE__)
        end subroutine skip_commentlines_pp

    end subroutine m_PP_rd_ival

end module m_PseudoPotential
