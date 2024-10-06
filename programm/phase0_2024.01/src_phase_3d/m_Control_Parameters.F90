!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev$)
!
!  MODULE: m_Control_Parameters
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004, May/09/2004
!                                    , May 2005,  Jan 2010
!                        J. Koga,  March/01/2010
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!#========================================================================
!
!   patch 0.1 by K. Tagami@adv    2009/05/28
!   patch 0.2 by K. Tagami@adv    2009/10/19
!
!   patch 0.1:  correction for DFT+U by introducing prealloc_kt
!   patch 0.2:  correction for phonon calculation with DFT+U
!
!
!   patch 10.1 by K. Tagami@adv    2011/06/18
!
!   patch 10.1:  introduction of sw_LinearResponse
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
module m_Control_Parameters
!     (m_CtrlP)
! $Id$
!
! This module "m_Control_Parameters" holds parameters that give
! methods and calculational conditions in jobs.
!
! Parameters had been transported from the program coded with
! VPP-fortran.
! Translation from VPP-fortran version to this f90+mpi module was done
! by T. Yamasaki (JRCAT-ATP, FUJITSU LABORATORIES Ltd.) in 1999.
!
  use m_Const_Parameters, only : by_matrix_diagon, by_random_numbers &
       &, by_pseudo_atomic_orbitals &
       &, INITIAL, CONTINUATION, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &, COORDINATE_CONTINUATION &
       &, AUTOMATIC, FIXED_CHARGE_AUTOMATIC, PREPARATION_ONLY, EK_CONVERGED &
! --------------------- T. Hamada 2021.9.23 ------------------- ! UVSOR
       &, NO_SAVE, SAVE, READ, PREP &
! ------------------------------------------------------------- ! UVSOR
       &, SD, SDPRC, MSD, MSDPRC, CG, CGPRC, RMM, RMMPRC, RMM2P, RMM2PPRC &
       &, RMM2, RMM2PRC, RMM3, SDLM, SDLMPRC, MSDLM, MSDLMPRC, bbCG, bbCGPRC &
       &, eazyCG, lmeazyCG, eazyCGPRC, MATRIXDIAGON &
       &, default_sd2cg, Gauss_distrib_func, VERY_NARROW, from_PSEUDOPOTENTIAL_FILE &
       &, from_wave_functions, DP, SP &
       &, lmSD, lmMSD, lmCG, ON, OFF, GRID,YES, NO, MDOLD, MDSMPL, SUBMAT, DAVIDSON, MDDAVIDSON &
       &, MDKOSUGI, UNIT_MATRIX, SPIN_POLARIZED, CRYSTAL_FIELD_APPROX &
       &, SMEARING, FIXED, PARABOLIC, MP, TETRAHEDRON, COLD, LOWEST_AT_EACH_KPT &
       &, irsize, START, END, FINISH &
       &, UMICRO, SIMPLE, BROYD1, BROYD2, DFP, PULAY, ANEW, RENEW&
       &, SMALL, MEDIUM, LARGE, len_tag_isolver, tag_isolver &
       &, BEFORE, AFTER, varLINEAR, varTANH, QUENCHED_MD, GDIIS, VERLET, BFGS, L_BFGS &
       &, T_CONTROL, BLUEMOON, QUENCHED_CONSTRAINT, CG_STROPT, CG_STROPT2, STEEPEST_DESCENT, FIRE &
       &, FILE, SKPS_DIRECT_IN, MESH, GAMMA, MONKHORST_PACK &
       &, HR_ENERGY_UNIT, EV_ENERGY_UNIT &
       &, Hartree, OLD, NEW_, FMAXVALLEN, FMAXUNITLEN, CUBE, VTK, BINARY, DENSITY_ONLY, INTEGRATED, SEPARATE &
       &, PARA, ANTIFERRO, FERRO, NOCONV, LOWER, UPPER, BULK, DEFECT &
       &, Wavefunction, REGULAR_INTERVALS, BY_ATOMIC_POSITIONS, TAG_FORMAT,TAG_LINE,TABLE &
       &, NOPOLAR, POLARIZATION, EFFECTIVE_CHARGE, PIEZOELECTRIC_CONST &
       &, SPHERICAL_HARMONICS, ATOMIC_ORBITAL, WANNIER_FUNCTION &
!!$       &, HOUSEHOLDER, DIVIDEandCONQUER, LIMIT_1DSEARCH &
       &, HOUSEHOLDER, DIVIDEandCONQUER &
       &, TOTAL_ENERGY, MODIFIED_TOTAL_ENERGY, BAND_ENERGY &
       &, PHONON_GAMMA, PHONON_BAND, PHONON_DOS, ONE_BY_ONE, ALL_AT_ONCE, DELTA10 &
       &, DAMPED_MD, VELOCITY_SCALING, DRIVER_GENERAL, NEWTON, LAGRANGE &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!      &, DRIVER_CONSTRAINT, DRIVER_NEB,  DRIVER_MTD
       &, DRIVER_CONSTRAINT, DRIVER_NEB,  DRIVER_MTD, DRIVER_URAMP, DRIVER_SC_DFT, DRIVER_DIMER &
       &, EXPLICIT, LMM_RMM, DAV_RMM, FEF &
       &, VDW_WILLIAMS, VDW_GRIMME, VDW_DFTD3 &
       &, BARE, PE1, PE2 &
       &, MASK_FUNCTION, PREFITTING &
       &, P_CONTROL, PT_CONTROL &
       &, DELTA_ENERGY, DELTA_MOVING_AVERAGE, SLOPE, DELTA_V &
!!$       &, FFT_redundant, FFT_parallel &
       &, WIDE, NARROW, eVunit &
       &, LINEAR_INTERPOLATION, SPLINE_INTERPOLATION &
       &, EXACT, FINE, MODERATE, COARSE, UNDETERMINED, MAX_TIME_REACHED &
       &, ALL, XSF, N2P2, DEEPMD &
       &, BOHR
! ==============================================================================
  use m_Parallelization, only   : MPI_CommGroup,npes,mype,ierr,nrank_e,nrank_k, nrank_g
! ====================================added by K. Tagami =================5.0
  use m_Const_Parameters, only : OccMat_type1, OccMat_Type2, &
       &                         Ueff_From_First, Ueff_Gradually, FLL, AMF, &
       &                         DIAG_CHARGE_DENSITY_MATRIX, &
       &                         DIAG_SPIN_DENSITY_MATRIX, DIAG_LS, &
       &                         LOCAL_POINT_GROUP, LOCAL_DOUBLE_POINT_GROUP
! ======================================================================== 5.0

! ============================ added by K. Tagami ============== 11.0
  use m_Const_Parameters,   only : NONCOLLINEAR,Neglected, BuiltIn, &
       &                           ByProjector, ByPawPot, ZeffApprox, ReadFromPP
! ============================================================== 11.0
  use m_ErrorMessages,        only : INVALID_CHARGE_MIXING

! ====== KT_add ========================================= 13.0E, 13.0U3, positron
  use m_Const_Parameters,     only : FERMI_DIRAC, CONST_kB, STEPWISE, &
       &                             Positron_CONV, Positron_GGGC, Positron_PSN
! ======================================================= 13.0E, 13.0U3, positron

#ifdef LIBXC
  use xc_f03_lib_m
#endif
  use mpi

  implicit none
!  include 'mpif.h'

  integer            :: ekmode=OFF
  integer            :: uvsormode=OFF
  integer            :: multiple_replica_mode = OFF
!!  integer            :: multiple_replica_max_iteration = 100
  integer            :: multiple_replica_max_iteration = -1
  integer,private    :: ntcnvg=0
!!  integer,private,parameter  :: MTIMES_CONVERGENCE = 3
  integer,private,parameter  :: MTIMES_CONVERGENCE = 2
!!!!  integer,private ::    mtimes_convergence_scf = MTIMES_CONVERGENCE
  integer ::    mtimes_convergence_scf = MTIMES_CONVERGENCE
  integer,private ::    mtimes_convergence_ek  = MTIMES_CONVERGENCE
  integer            :: kimg  ! {1|2} 1: when inversion_symmetry == ON,  2:OFF
  real(kind=DP)      :: gmax, gmaxp, gmaxs_given, gmaxp_reduced, gmax_org &
       &              , gmax_buf=0.d0
  integer            :: n_matrix_size = 0
  integer ::            icond = AUTOMATIC
  integer ::            icond_org = AUTOMATIC
  integer ::            precision_WFfile = SP
  integer ::            continuation_using_ppdata = NO
  integer ::            fixed_charge_k_parallel = ALL_AT_ONCE ! {ALL_AT_ONCE|ONE_BY_ONE}
  integer ::            ipriekzaj = 1
  integer ::            iexpl
  integer ::            ipri = 1 &
       &              , ipritiming = -1, ipritiming0 = -1, iprisolver = -1, iprievdff = -1 &
       &              , iprirmm = -1, iprisnl = -1 &
       &              , ipri_spg, ipri_kp, ipripulay = -1, iprimatdiagon = -1 &
       &              , iprivlhxcq = -1, iprigdiis = -1, iprieigenvalue = -1 &
       &              , ipritotalcharge = -1, iprisubmat = -1, ipristrcfctr = -1 &
       &              , ipriinputfile = -1, ipriparallel = -1, iprifftmap = -1 &
       &              , iprivloc = -1, iprimd = -1, iprioccup = -1, iprichargemixing = -1 &
       &              , ipripositron = -1, ipriforce= -1, ipripp=-1, ipridos = -1 &
       &              , iprinegativecharge = -1, iprinegativecharge0 = -1 &
       &              , ipriphig = -1, ipripao = -1, ipriberry = -1, ipriphonon = -1 &
       &              , ipriparadeb = 0, iprijobstatus = -1 &
       &              , ipridavidson = -1, iprihubbard = -1, ipriwf = -1, iprichargedensity = -1 &
       &              , iprivelocity = -1, ipribetar = -1, ipripaw = -1, iprixc = -1 &
       &              , iprifef = -1, iprimddavidson = -1, iprimdkosugi = -1 &
       &              , ipriesm = -1, ipripredictor = -1, ipriunitcell = -1, ipricoefwf = -1 &
       &              , iprirs  = -1, iprirsb = -1, iprisym = -1, iprifcp = -1, iprivdw = -1 &
       &              , iprifire = -1, ipriexx = -1, ipriblsize = -1
  integer ::            sw_timing_2ndlevel = OFF

  integer ::            sw_flatten = OFF, sw_firstlevel_only = ON
  integer ::            sw_details = OFF
  integer ::            measure_count_limit = -1
  integer ::            n_fermi_vicinity=6

  integer, parameter :: Nw_Psicoef = 100
  integer, private,parameter :: lmm_status_store_size = 10
  integer, private,dimension(lmm_status_store_size) :: lmm_status_store
  integer, private       :: lmm_status_pointer = 1
  integer, private       :: lmm_status_stored = 0
  logical                :: in_line_minimization = .false.

!!  ipri_spg,ipri_kp were introduced by mizouchi@adv 2003.2.21  !!
  logical ::            printable = .false.
  logical ::            ppprinted = .false.
  integer, private, parameter    :: N_SUBROUTINES = 10
  integer ::            num_subroutines = N_SUBROUTINES
  integer ::            statistics_in_parallel = 0
  real(DP)::            PCPUDF = 0.03
  integer ::            max_warnings_negativecharge = 10
  integer ::            jobstatus_series = OFF, jobstatus_format = TAG_FORMAT
  integer            :: iconvergence = 0 ! {0|1|2|3}
  integer            :: iconvergence_previous_job = 0
  integer, allocatable, dimension(:) :: iconv_ek_tmp
  integer ::            numk_tmp = 0
  integer ::            numk_zajsaved = 0
  character(len("PCPUDF")), private, parameter :: tag_pcpudf = 'PCPUDF'
  character(len("convergence")),private,parameter :: tag_convergence = "convergence"
  character(len("numk")),private, parameter :: tag_numk = "numk"
  character(len("numk_zajsaved")),private, parameter ::  tag_numk_zajsaved = "numk_zajsaved"
  character(len("convergence_ek")),private, parameter :: tag_convergence_ek = "convergence_ek"

  integer ::            neg = 1       ! number of eigen values for each k-point
  logical ::            neg_is_given = .false.
  integer ::            neg_previous  ! number of eigen values for each k-point in the previous job
  integer ::            neg_fixed = 0 ! number of eigen values that are fixed
  integer ::            num_extra_bands = 0 ! number of extra eigen values for each k-point.
  logical ::            num_extra_bands_was_set = .false.
  !           This is fixed zero when ekmode == OFF, and a default value is 1 when ekmode == ON
  logical ::            neg_is_enlarged = .false.

  integer            :: sw_ekzaj = OFF

  ! Additional projectors
  integer,save :: sw_use_add_proj = OFF

  ! k-point data
! -------------------- T. Hamada 2021.9.23 --------------------
  integer      :: sw_kptdata = NO_SAVE                  ! UVSOR
! -------------------------------------------------------------

  integer, parameter :: len_xctype = 7

  integer, parameter :: LEN_TITLE = 80
  ! ------- Positron start
  integer,save ::       sw_positron = OFF
  integer ::            corecharge_cntnbin = OFF
  real(kind=DP)      :: gmax_positron = 0.d0
  character(len=len_xctype) :: epctype = 'nogiven'
  integer,private ::    positron_ntcnvg = 0
  integer ::            npeg = 1
  integer ::            num_extra_pev = 0 ! number of extra positron eigen values
  real(kind=DP) ::      delta_pev = 1.d-15
  integer,private ::    mtimes_convergence_pev = MTIMES_CONVERGENCE
  integer ::            pev_max_iteration = 300
  integer ::            evaluation_pev_diff = ON
  real(kind=DP) ::      dtim_p = 1.d0
  integer ::            intpzaj = by_random_numbers
  integer ::            sw_submat_p = ON
  integer ::            isolver_p = lmMSD
  integer ::            sw_gga_p = OFF
  integer ::            sw_epsilon_ele = OFF
  real(kind=DP) ::      epsilon_ele
  integer ::            sw_positron_file = ON
  integer ::            positron_filetype = CUBE

  integer ::  positron_method = Positron_CONV

  character(len("positron_method")),private,parameter :: &
       &             tag_positron_method    = "positron_method"
  character(len("conv")),private,parameter ::         tag_conv    = "conv"
  character(len("gggc")),private,parameter ::         tag_gggc    = "gggc"
  character(len("psn")),private,parameter ::          tag_psn     = "psn"

  character(len=LEN_TITLE) ::    positron_title(5)
  data positron_title / &
       &  "positron density", "valence electron density", "e-p pair density" &
       & ,"gradient of valence electron density","e-p pair density (minority)"/
!!$  character(len=LEN_TITLE) ::    positron_title_positr = 'positron density'
!!$  character(len=LEN_TITLE) ::    positron_title_velect = 'valence electron density'
!!$  character(len=LEN_TITLE) ::    positron_title_eppair = 'e-p pair density'
!!$  character(len=LEN_TITLE) ::    positron_title_velgrd = 'gradient of valence electron density'

  character(len("positron_convergence")),private,parameter :: tag_positron_convergence = "positron_convergence"
  character(len("positron")),private,parameter ::     tag_positron    = "positron"
  character(len("bulk")),private,parameter ::         tag_bulk        = "bulk"
  character(len("defect")),private,parameter ::       tag_defect      = "defect"
  character(len("ipripositron")),private,parameter ::   tag_ipripositron      = "ipripositron"
  character(len("dtim")),private,parameter ::         tag_dtim        = "dtim"
  character(len("sw_submat")),private,parameter ::    tag_sw_submat   = "sw_submat"
  character(len("corecharge_cntnbin")),private,parameter ::  tag_corecharge_cntnbin = "corecharge_cntnbin"
  character(len("solver_for_positronWF")),private,parameter :: tag_solver_for_positronWF = "solver_for_positronWF"
  character(len("solver")),private,parameter ::       tag_solver      = "solver"
  character(len("sw_gga")),private,parameter ::       tag_sw_gga = "sw_gga"
  character(len("sw_gga_p")),private,parameter ::     tag_sw_gga_p = "sw_gga_p"
  character(len("epsilon_ele")),private,parameter ::  tag_epsilon_ele = "epsilon_ele"
  character(len("positron_file")),private,parameter ::tag_positron_file = "positron_file"
  character(len("sw_positron_file")),private,parameter:: tag_sw_positron_file="sw_positron_file"
  character(len("title_positron")),private,parameter ::tag_title_positr = "title_positron"
  character(len("title_electron")),private,parameter ::tag_title_velect = "title_electron"
  character(len("title_eppair")),private,parameter ::  tag_title_eppair =   "title_eppair"
  character(len("title_electron_gradient")),private,parameter :: tag_title_velgrd = "title_electron_gradient"
  ! ------- Positron end

  ! ------- LDOS -----
  integer,public ::      dos_write_format = NARROW  ! WIDE, NARROW, or eVunit
  integer, save  ::      set_write_format_count = 0
  integer,public ::      sw_checksum = ON

  integer ::             sw_ldos = OFF, sw_aldos = OFF, sw_layerdos = OFF
  integer,public ::      ldos_method = Gauss_distrib_func
  integer,public ::      sw_save_ldos_weight = ON
  integer,public ::      sw_cal_ldos = ON
  logical,public ::      crtdst_is_given = .false.
  real(kind=DP) ::       crtdst_aldos = 6.0d0
  integer,public ::      naldos_from = 0, naldos_to = 0
  integer ::             slicing_way_winlay = REGULAR_INTERVALS
  integer,public ::      integration_dimension_winlay = 1
  real(kind=DP) ::       deltaz_winlay = 0.5
  integer ::             normal_axis_winlay   = 1
  real(kind=DP) ::       crtdst_winlay = 3.5d0
  integer ::             hardpart_subroutine = 0
  integer ::             sw_rspace_ldos = OFF
  integer ::             sw_ac_mesh = OFF
  integer ::             acmesh_factor = 1

!!$  character(len("ldos_hardpart_fft")),private,parameter ::    tag_ldos_hardpart_fft = "ldos_hardpart_fft"
!!$  character(len("redundant")),private,parameter ::   tag_redundant  = "redundant"
!!$  character(len("parallel")),private,parameter ::    tag_parallel   = "parallel"
  character(len("dos_write_format")),private,parameter :: tag_dos_write_format = "dos_write_format"
  character(len("write_format")),private,parameter :: tag_write_format = "write_format"
  character(len("wide")),private,parameter ::         tag_wide        = "wide"
  character(len("narrow")),private,parameter ::       tag_narrow      = "narrow"
  character(len("evunit")),private,parameter ::       tag_evunit      = "evunit"
  character(len("eV")),private,parameter ::           tag_eV          = "eV"
  character(len("sw_checksum")),private,parameter :: tag_sw_checksum  = "sw_checksum"

  character(len("ldos")),private,parameter ::        tag_ldos         = "ldos"
  character(len("sw_aldos")),private,parameter ::    tag_sw_aldos     = "sw_aldos"
  character(len("sw_atomicdos")),private,parameter :: tag_sw_atomicdos = "sw_atomicdos"
  character(len("sw_layerdos")),private,parameter :: tag_sw_layerdos  = "sw_layerdos"
  character(len("aldos")),private,parameter ::       tag_aldos        = "aldos"
  character(len("atomicdos")),private,parameter ::   tag_atomicdos     = "atomicdos"
  character(len("sw_save_ldos_weight")),private,parameter :: tag_sw_save_ldos_weight = "sw_save_ldos_weight"
  character(len("sw_cal_ldos")),private,parameter :: tag_sw_cal_ldos  = "sw_cal_ldos"
  character(len("crtdst")),private,parameter ::      tag_crtdst       = "crtdst"
  character(len("naldos_from")),private,parameter :: tag_naldos_from  = "naldos_from"
  character(len("naldos_to")),private,parameter ::   tag_naldos_to    = "naldos_to"
  character(len("layerdos")),private,parameter ::    tag_layerdos     = "layerdos"
  character(len("slicing_way")),private,parameter :: tag_slicing_way  = "slicing_way"
  character(len("deltaz")),private,parameter ::      tag_deltaz       = "deltaz"
  character(len("normal_axis")),private,parameter :: tag_normal_axis  = "normal_axis"
  character(len("regular_intervals")),private,parameter :: tag_regular_intervals = "regular_intervals"
  character(len("integration_dimension")),private,parameter:: tag_integration_dimension = "integration_dimension"
!!$  character(len("by_atoms")),private,parameter ::    tag_by_atoms     = "by_atoms"
  character(len("by_atomic_positions")),private,parameter :: tag_by_atomic_positions = "by_atomic_positions"
  character(len("hardpart_subroutine")),private,parameter :: tag_hardpart_subroutine = "hardpart_subroutine"
  character(len("sw_ac_mesh")),private,parameter :: tag_sw_ac_mesh = "sw_ac_mesh"
  character(len("sw_atom_centered_mesh")),private,parameter :: tag_sw_atom_centered_mesh = "sw_atom_centered_mesh"
  character(len("sw_softpart_only")),private,parameter :: tag_sw_softpart_only = "sw_softpart_only"
  character(len("ac_mesh_factor")),private,parameter :: tag_acmesh_factor = "ac_mesh_factor"
  character(len("atom_centered_mesh_factor")),private,parameter :: tag_atom_centered_mesh_factor = "atom_centered_mesh_factor"
  ! ------------------

  ! ------- LBAND -----
  character(len("lband")),private,parameter :: &
       &        tag_lband = "lband"
  character(len("wf_local_decomposition")),private,parameter :: &
       &        tag_wf_local_decomposition = "wf_local_decomposition"
  character(len("atom_decomposition")),private,parameter :: &
       &        tag_atom_decomposition = "atom_decomposition"
  character(len("layer_decomposition")),private,parameter :: &
       &        tag_layer_decomposition = "layer_decomposition"
  character(len("sw_calc_wf_atom_decomposition")),private,parameter :: &
       &        tag_sw_calc_wf_atom_decomp = "sw_calc_wf_atom_decomposition"
  character(len("sw_calc_wf_layer_decomposition")),private,parameter :: &
       &        tag_sw_calc_wf_layer_decomp = "sw_calc_wf_layer_decomposition"
  character(len("sw_wd_only_specified_atoms")),private,parameter :: &
       &        tag_sw_wd_only_specified_atoms = "sw_wd_only_specified_atoms"

  integer ::  sw_calc_wf_atom_decomposition = OFF, &
       &      sw_calc_wf_layer_decomposition = OFF
  integer :: sw_lband = off
  integer :: sw_rspace_lband = off
  integer :: sw_wd_only_specified_atoms = OFF

  ! ------ Dipole ----
  integer :: sw_dipole = OFF
  integer :: sw_layered = OFF
  integer :: sw_dipole_correction = OFF
  integer :: idir_dip = 0 ! 0:all, 1:x, 2:y, 3:z
  integer :: ndiv_dip = 100
  real(kind=DP) :: rvac(3) = (/ 0.d0, 0.d0, 0.d0 /)
  real(kind=DP) :: width_dip = 1.d-2
  real(kind=DP) :: elec_field(3) = (/ 0.d0, 0.d0, 0.d0 /)
!  real(kind=DP) :: amix_dip = 1.d-1
  real(kind=DP) :: amix_dip = 1.d0
  character(len("dipole")),private,parameter ::    tag_dipole     = "dipole"
  character(len("sw_dipole")),private,parameter ::    tag_sw_dipole     = "sw_dipole"
  character(len("sw_layered")),private,parameter ::    tag_sw_layered     = "sw_layered"
  character(len("direction")),private,parameter ::    tag_direction     = "direction"
  character(len("division")),private,parameter ::    tag_division     = "division"
  character(len("width")),private,parameter ::    tag_width     = "width"
  character(len("vacuum")),private,parameter ::    tag_vacuum     = "vacuum"
  !!$character(len("rx")),private,parameter ::    tag_rx     = "rx"
  !!$character(len("ry")),private,parameter ::    tag_ry     = "ry"
  !!$character(len("rz")),private,parameter ::    tag_rz     = "rz"
  character(len("dipole_correction")),private,parameter ::    tag_dipole_correction     = "dipole_correction"
  character(len("sw_dipole_correction")),private,parameter ::    tag_sw_dipole_correction     = "sw_dipole_correction"
  character(len("electric_field")),private,parameter :: tag_electric_field = "electric_field"
  character(len("ex")),private,parameter ::    tag_ex     = "ex"
  character(len("ey")),private,parameter ::    tag_ey     = "ey"
  character(len("ez")),private,parameter ::    tag_ez     = "ez"
  ! ------------------

  ! --- Screening Correction ------
  integer :: sw_screening_correction = OFF
  real(kind=DP) :: screening_alpha = 1.d0
  character(len("sw_screening_correction")),private,parameter :: tag_sw_screening_correction = "sw_screening_correction"
  character(len("screening_correction")),private,parameter :: tag_screening_correction = "screening_correction"
  character(len("alpha")),private,parameter :: tag_screening_alpha = "alpha"
  ! --- External Potential (External Charge) ------
  integer :: sw_external_potential = OFF
  ! ------------------

  ! --- fermi surface ----
  character(len("fermi_surface")),private,parameter :: &
       &     tag_fermi_surface = "fermi_surface"
  character(len("sw_write_bxsf_file")),private,parameter :: &
       &     tag_sw_write_bxsf_file = "sw_write_bxsf_file"
  integer :: sw_write_bxsf_file = OFF
  ! ------------------

  ! ------ Maximally localized Wannier functions ----
  integer :: sw_wannier = OFF
  integer :: sw_random_wannier = OFF
  integer :: sw_potential_wannier = OFF
  integer :: sw_continue_wannier = OFF
  integer :: wannier_opt_method = SD
  integer :: max_iter_wan = 1000
  integer :: wannier_filetype = DENSITY_ONLY
  real(kind=DP) :: eps_wan = 1.d-3
  real(kind=DP) :: dt_wan  = 1.d-4
  character(len("wannier")),private,parameter :: tag_wannier = "wannier"
  character(len("sw_wannier")),private,parameter :: tag_sw_wannier = "sw_wannier"
  character(len("sw_random_wannier")),private,parameter :: tag_sw_random_wannier = "sw_random_wannier"
  character(len("sw_continue")),private,parameter :: tag_sw_continue = "sw_continue"
  character(len("sw_potential")),private,parameter :: tag_sw_potential = "sw_potential"
  character(len("eps_grad")),private,parameter :: tag_eps_grad = "eps_grad"
  ! ------------------

  ! ------ Homogeneous Finite Electric Field Method ----
  integer :: sw_fef = OFF
  integer :: sw_check_polar = OFF
  character(len("fef")),private,parameter :: tag_fef = "fef"
  character(len("sw_fef")),private,parameter :: tag_sw_fef = "sw_fef"
  character(len("sw_check_polar")),private,parameter :: tag_sw_check_polar = "sw_check_polar"
  ! ------------------

  ! --- van der Waals (pair-interaction approximation) ------
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
! integer :: sw_pair_vdw = OFF
! real(kind=DP) :: rcut_vdw = 20.d0 ! vdW cutoff
! character(len("sw_pair_vdw")),private,parameter :: tag_sw_pair_vdw = "sw_pair_vdw"
! character(len("rcut_vdw")),private,parameter :: tag_rcut_vdw = "rcut_vdw"
  integer :: sw_vdw_correction = OFF
  integer :: vdw_method = 0
  real(kind=DP) :: vdw_radius = 20.d0 ! vdW cutoff
  real(kind=DP) :: vdw_scaling_factor = 1.0d0
  real(kind=DP) :: vdw_scaling_factor_r = 1.0d0
  real(kind=DP) :: vdw_damping_factor = 1.0d0
  character(len("sw_vdw_correction")),private,parameter :: tag_sw_vdw_correction = "sw_vdw_correction"
  character(len("vdw_method")),private,parameter :: tag_vdw_method = "vdw_method"
  character(len("williams")),private,parameter :: tag_williams = "williams"
  character(len("grimme")),private,parameter :: tag_grimme = "grimme"
  character(len("dft-d2")),private,parameter :: tag_dft_d2 = "dft-d2"
  character(len("dft-d3")),private,parameter :: tag_dft_d3 = "dft-d3"
  character(len("dft-d3")),private,parameter :: tag_dftd3 = "dftd3"
  character(len("vdw_radius")),private,parameter :: tag_vdw_radius = "vdw_radius"
  character(len("vdw_scaling_factor")),private,parameter :: tag_vdw_scaling_factor = "vdw_scaling_factor"
  character(len("vdw_scaling_factor_r")),private,parameter :: tag_vdw_scaling_factor_r = "vdw_scaling_factor_r"
  character(len("vdw_damping_factor")),private,parameter :: tag_vdw_damping_factor = "vdw_damping_factor"
! ==============================================================================
! --- DFT-D3 -------
  character(len("damping_function")),private,parameter :: &
       &             tag_dftd3_damping_function = "damping_function"
  character(len("zero")),private,parameter :: &
       &             tag_dftd3_damping_zero = "zero"
  character(len("bj")),private,parameter :: &
       &             tag_dftd3_damping_bj = "bj"
  integer :: dftd3_damping_function = 0            ! zero damping

  ! ------ Wannier90 ----
  integer :: sw_wannier90 = OFF
  integer :: nb_wan90 = 0
  integer :: sw_use_hardpart_wan90 = ON
  integer :: sw_write_unk_file = ON
  integer :: spin_component_wan90 = 1

  character(len=LEN_TITLE) :: wan90_seedname = "wan90"
  character(len("sw_wannier90")),private,parameter :: tag_sw_wannier90 = "sw_wannier90"
  character(len("seedname")),private,parameter :: tag_seedname         = "seedname"
  character(len("nb_wan90")),private,parameter :: tag_nb_wan90         = "nb_wan90"
  character(len("sw_use_hardpart_wan90")),private,parameter :: &
       &           tag_sw_use_hardpart_wan90 = "sw_use_hardpart_wan90"
  character(len("sw_write_unk_file")),private,parameter :: &
       &           tag_sw_write_unk_file = "sw_write_unk_file"
  character(len("spin_component_wan90")),private,parameter :: &
       &                   tag_spin_component_wan90      = "spin_component_wan90"

  ! ------------------

  integer            :: af = 0  ! {0|1} antiferro?,  = 0(.not.antiferro), = 1(antiferro)
  logical            :: max_scf_iteration_is_given = .false.
  integer            :: max_scf_iteration = 300    ! max_scf_iteration in one mdstep
  logical ::            max_TS_iteration_is_given = .false.
  logical ::            max_mdstep_is_given   = .false.
!  integer ::            max_total_scf_iteration = 10000  ! max_total_scf_iteration
  integer ::            max_total_scf_iteration = HUGE(1)  ! max_total_scf_iteration
  integer ::            max_mdstep = 10000
  integer            :: max_scdft_iteration = 10
  logical            :: max_scdft_iteration_is_given = .false.
  integer ::            iter_last, ek_max_iteration = 300
  integer            :: ifstop = 1
  integer            :: istop = -1
  integer            :: ldx = 1, ldy = 0, ldz = 0    ! mpifft
  real(kind=DP)      :: forccr = 1.d-3
  real(kind=DP)      :: force_error_check_rangeL = 0.d0
  integer            :: force_error_check_mode = OFF
  real(kind=DP)      :: f_tolerable_norm_error = 1.0
  real(kind=DP)      :: f_tolerable_angle_error = 90
  real(kind=DP)      :: f_tolerable_hyper_angle_error = 90

  real(kind=DP)      :: destm, dtim_1234(4),dtim_initial &
       & ,              dtim_1Dsearch = -1
  real(kind=DP)      :: rmx_1234(4)=0.5, rmx
!!$  real(kind=DP)      :: edelta = 1.d-10
  real(kind=DP)      :: edelta = 1.d-9
!!$  real(kind=DP)      :: edelta_initial = 1.d-10, max_force_edelta_i = 1.d-1, edelta_ontheway = 1.d-10
  real(kind=DP)      :: edelta_initial = 1.d-9, max_force_edelta_i = 1.d-1, edelta_ontheway = 1.d-9
  real(kind=DP)      :: edelta_sampling = 1.d-9
  logical            :: edelta_initial_is_given = .false., max_force_edelta_i_is_given = .false.
  integer, private   :: sub_mtimes_convergence = MTIMES_CONVERGENCE
  integer            :: convergence_criteria = DELTA_ENERGY
  real(kind=DP)      :: sub_delta_factor = 1.0
  logical            :: sub_delta_factor_is_given = .false.
  integer,private    :: sub_ntcnvg = 0
  integer            :: sw_fix = OFF
  real(kind=DP)      :: delta_eigenvalue_conduction = 1.d+15
  logical            :: delta_eigenvalue_cond_is_given = .false.
  integer, private   :: flag_wd_force = 0
  real(kind=DP), parameter :: MAX_FORCE_EDELTA_I_FACTOR = 1.d+4
  real(kind=DP), parameter :: FORCCR_FACTOR_EDELTA_ONTHEWAY = 2.d0
  real(kind=DP), parameter :: Edelta_Critical_value = 1.d-20
!!  real(kind=DP)      :: delta_eigenvalue = 1.d-15
  real(kind=DP)      :: delta_eigenvalue = 1.d-5
  integer            :: way_of_smearing = PARABOLIC
  real(kind=DP)      :: width = 0.001
  real(kind=DP)      :: width_tetra = 1.e-5
  integer            :: idimtetra = 3
  integer            :: sw_correction = OFF
  integer            :: order_mp = 2 ! default order for the Methfessel-Paxton smearing
  real(kind=DP)      :: esearch_factor_mp = 10.d0
  integer            :: esearch = on
  real(kind=DP)      :: cpumax = 86400.d0
  integer            :: waymix, istrbr, nbxmix, hownew, cutoff_mix = LARGE
  integer            :: waymix_b, istrbr_b, nbxmix_b, hownew_b, cutoff_mix_b = LARGE
  logical            :: decomp_into_charge_spin
  logical            :: c_precon
  integer            :: ncachesize_given = -1
  integer            :: sw_use_wfred = OFF
  integer            :: nblocksize_dgemm
  integer            :: nblocksize_mgs
#ifdef SAVE_FFT_TIMES
  integer            :: sw_save_fft = OFF
#endif
  integer            :: recursivesize_mgs
  integer            :: nblocksize_betar_dot_wfs
  integer            :: nblocksize_betar_dot_wfs_nlmta
  integer            :: nblocksize_vnonlocal_w
  integer            :: nblocksize_submat = 0
  integer            :: nblocksize_betar_dot_wfs_npe
#ifndef NO_FORCE_DGEMM
  integer            :: nblocksize_force
#endif
  integer            :: nblocksize_vnonlocal_w_nlmta
  logical            :: nblocksize_vnonlocal_w_nlmta_is_given = .false.

  logical            :: nblocksize_fftw_is_given = .false.
  integer            :: nblocksize_fftw = 32
  integer            :: nblocksize_vnonlocal_w_f
  logical            :: nblocksize_vnonlocal_w_f_is_given = .false.
  integer            :: nblocksize_gather_f
  logical            :: nblocksize_gather_f_is_given = .false.
  integer            :: nblocksize_subspace = 0
  integer            :: fftbox_divide_cube = 0
  integer            :: divide_square = 0
  integer            :: fftbox_3ddiv_1 = 0
  integer            :: fftbox_3ddiv_2 = 0
  integer            :: fftbox_3ddiv_3 = 0
  integer            :: fftbox_div_1 = 0
  integer            :: fftbox_div_2 = 0
  integer            :: sw_fft_xzy   = 0
  integer            :: nblocksize_submat_latter = 0
  logical            :: nblocksize_dgemm_is_given = .false.
  logical            :: nblocksize_mgs_is_given = .false.
  logical            :: recursivesize_mgs_is_given = .false.
  logical            :: nblocksize_betar_is_given = .false.
  logical            :: nblocksize_betar_nlmta_is_given = .false.
  logical            :: nblocksize_vnonlocal_is_given = .false.
  logical            :: nblocksize_submat_is_given = .false.
  logical            :: nblocksize_submat_latter_is_given = .false.
  logical            :: nblocksize_betar_npe_is_given = .false.
#ifndef NO_FORCE_DGEMM
  logical            :: nblocksize_force_is_given = .false.
#ifdef SX
  integer,parameter :: nb_force_default     = 5000
#else
  integer,parameter :: nb_force_default     = 32
#endif
#endif

  integer           :: nblocksize_rspace_betar = 1
  integer           :: nblocksize_rspace_v = 1
  logical            :: submat_uncalled = .true.
#ifdef SX
  integer,parameter :: nb_mgs_default       = 256
  integer,parameter :: nb_betar_default     = 10000
  integer,parameter :: nb_vnonlocal_default = 5000
  integer,parameter :: nb_submat_default    = 5000
#else
  integer,parameter :: nb_mgs_default       = 8
  integer,parameter :: nb_betar_default     = 32
  integer,parameter :: nb_vnonlocal_default = 1000
  integer,parameter :: nb_submat_default    = 32
#endif
! <--

  real(kind=DP)      :: amix = 1.0, bmix = -1.0, amin = -1.0d0
#ifndef _EMPIRICAL_
  integer            :: itercg, imax, indxcg
  integer            :: way_ksample = MONKHORST_PACK
  integer, parameter :: M_STM = 3
  integer            :: n_stm
  integer            :: intzaj = by_random_numbers
  integer            :: imatrix_diagon, imsd, i_2lm, i_sd2another
  integer            :: n_WF_solvers = 1
  integer,allocatable,dimension(:) :: WF_solver, till_n_iteration
  integer, private   :: previous_solver = 0
  integer            :: iwrksz_rmm_phase, evaluation_eko_diff
  real(kind=DP)      :: fftsize_factor_gmaxp = 1.d0
  real(kind=DP)      :: fftsize_factor_gmax  = 1.d0

!!$  integer, public ::   wf_inheritance = OFF
  character(len=len_xctype) ::  xctype = 'nogiven'
  integer            :: ggacmp_parallel = OFF
  integer            :: initial_chg = Gauss_distrib_func
  integer            :: sw_initial_charge_rspace = OFF
  integer            :: initial_charge_filetype = DENSITY_ONLY
  character(len("precalculation")),private,parameter :: tag_precalculation   = "precalculation"
  integer            :: nel_Ylm = 9 ! Number of (l,m) sets for preparation of Ylm
  character(len("nel_Ylm")),private,parameter ::      tag_nel_Ylm = "nel_Ylm"
  logical, public ::    sw_submat_is_on = .false.
  logical, public ::    renew_wf_again_m_CtrlP  = .false.
  logical, public ::    submat_is_done_this_iter = .false.

  integer,parameter  :: len_solvername = 10
  integer, save      :: number_of_solvers_applied  = 0
  integer,parameter  :: m_solvers_applied = 2
  character(len=len_solvername),dimension(m_solvers_applied) :: solver_names_applied
  integer,parameter  :: len_cdmixingname = 10
  integer, save      :: number_of_cdmixing_applied  = 0
  integer,parameter  :: m_cdmixing_applied = 2
  character(len=len_cdmixingname),dimension(m_cdmixing_applied) :: cdmixing_names_applied

  !--- Phonon ---
  integer            :: phonon_method = PHONON_GAMMA
  integer            :: sw_phonon = OFF
  integer            :: sw_calc_force = OFF
  integer            :: sw_calc_force_all = OFF
  logical            :: skip_alloc_phonon = .false.
  integer            :: sw_vibrational_modes = OFF
  integer            :: sw_lo_to_splitting = OFF
  integer            :: sw_lattice_dielectric_tensor = OFF
  integer            :: sw_dielectric_function = OFF
  integer            :: sw_correct_force_constants = ON
  integer            :: sw_polynomial_fit = OFF
  integer            :: norder = 1
  integer            :: sw_int_strain_piezo_tensor = OFF
  integer            :: sw_phonon_oneshot = OFF
  integer            :: num_phonon_calc_mode = 0

  !--- vibrational mode ---
  integer            :: sw_vibrational_mode = OFF
  integer            :: with_mode_effchg = NO

#else
  !--- Phonon ---
  integer            :: sw_phonon = OFF
  integer            :: sw_calc_force = OFF
  integer            :: sw_calc_force_all = OFF
  integer            :: sw_polynomial_fit = OFF
  integer            :: norder = 1
  integer            :: sw_phonon_oneshot = OFF
  integer            :: num_phonon_calc_mode = 0
  !--- vibrational mode ---
!!$  integer            :: sw_vibrational_modes = OFF
  integer            :: sw_vibrational_mode = OFF
  integer            :: with_mode_effchg = NO
#endif

! === KT_add === 13.1R
  integer            :: sw_raman = OFF
  integer            :: sw_phonon_with_epsilon = OFF
  integer            :: sw_calc_dielectric_tensor = OFF
! ============== 13.1R

  integer            :: nspin = 1         ! {1|2}

! ============================== added by K. Tagami ========== 11.0
  ! --- Noncollinear
!
  integer            :: ndim_spinor = 1  ! {1|2}
  integer            :: ndim_chgpot = 1
  integer            :: ndim_magmom = 1
  logical            :: noncol = .false.
! =============================================================== 11.0

  ! --- optimization of the lattice
  integer            :: sw_optimize_lattice = OFF
  integer            :: sw_uniform = OFF
  integer            :: sw_rebuild_pws = ON
  integer            :: sw_optimize_coordinates_once = OFF
  integer            :: lattice_optimization_method = BFGS
  real(kind=DP)      :: lattice_optimization_conv = 1.0d-5
  integer            :: nhistory_stress = 6
  real(kind=DP)      :: delta_stress = 1.d0
  real(kind=DP)      :: max_stress = 1.0d-6
  character(len("lattice")),private,parameter :: tag_lattice = "lattice"
  character(len("sw_optimize_lattice")),private,parameter :: tag_sw_optimize_lattice = "sw_optimize_lattice"
  character(len("sw_optimize_coordinates_once")),private,parameter :: tag_sw_optimize_coordinates_once &
       &      = "sw_optimize_coordinates_once"
  character(len("sw_rebuild_pws")),private,parameter :: tag_sw_rebuild_pws="sw_rebuild_pws"
  character(len("nhistory")),private,parameter :: tag_nhistory="nhistory"
  character(len("delta")),private,parameter :: tag_delta = "delta"
  character(len("max_stress")),private,parameter :: tag_max_stress="max_stress"
  character(len("stress_convergence")),private,parameter :: tag_stress_convergence = "stress_convergence"
  character(len('sw_uniform')),private,parameter :: tag_sw_uniform='sw_uniform'

  character(len("sw_optimize_coords_sametime")),private,parameter :: &
       &       tag_sw_optimize_coords_sametime = "sw_optimize_coords_sametime"
  character(len("lattice_coords_opt_mode")),private,parameter :: &
       &       tag_lattice_coords_opt_mode = "lattice_coords_opt_mode"
  character(len("omega_for_cellopt")),private,parameter :: &
       &       tag_omega_for_cellopt = "omega_for_cellopt"
  character(len("stress_force_mul_factor")),private,parameter :: &
       &       tag_stress_force_mul_factor = "stress_force_mul_factor"
  character(len("threshold_start_cellopt")),private,parameter :: &
       &       tag_threshold_start_cellopt = "threshold_start_cellopt"
  integer :: sw_optimize_coords_sametime = OFF
  integer :: lattice_coords_opt_mode = 2
  real(kind=DP) :: omega_for_cellopt = 1.5d0
  real(kind=DP) :: stress_force_mul_factor = 1.0d0
  real(kind=DP) :: threshold_start_cellopt = 1.0D1

! ========================== JK_add ================== 13.0AS
  character(len('external_stress')),private,parameter :: tag_external_stress='external_stress'
  character(len("s11")),private,parameter :: tag_s11 = "s11"
  character(len("s22")),private,parameter :: tag_s22 = "s22"
  character(len("s33")),private,parameter :: tag_s33 = "s33"
  character(len("s12")),private,parameter :: tag_s12 = "s12"
  character(len("s21")),private,parameter :: tag_s21 = "s21"
  character(len("s23")),private,parameter :: tag_s23 = "s23"
  character(len("s32")),private,parameter :: tag_s32 = "s32"
  character(len("s13")),private,parameter :: tag_s13 = "s13"
  character(len("s31")),private,parameter :: tag_s31 = "s31"
  real(kind=DP),dimension(3,3) :: external_stress
  logical :: estress_has_been_set = .false.
! =================================================== 13.0AS

! ========================== KT_add ===================== 13.0B
! -----------------------
! fixing lattice angles
! -----------------------
  character(len("fix_angle_alpha")),private,parameter :: &
       &                       tag_fix_angle_alpha = "fix_angle_alpha"
  character(len("fix_angle_beta")),private,parameter :: &
       &                       tag_fix_angle_beta = "fix_angle_beta"
  character(len("fix_angle_gamma")),private,parameter :: &
       &                       tag_fix_angle_gamma = "fix_angle_gamma"
!
  integer :: sw_fix_lattice_angles = OFF
  logical :: fix_angle_alpha = .false.
  logical :: fix_angle_beta  = .false.
  logical :: fix_angle_gamma = .false.
!
! -----------------------
! fixing lattice lengths ( added 13.1AS )
! -----------------------
  character(len("fix_length_a")),private,parameter :: &
       &                       tag_fix_length_a = "fix_length_a"
  character(len("fix_length_b")),private,parameter :: &
       &                       tag_fix_length_b = "fix_length_b"
  character(len("fix_length_c")),private,parameter :: &
       &                       tag_fix_length_c = "fix_length_c"
!
  integer :: sw_fix_lattice_lengths = OFF
  logical :: fix_length_a = .false.
  logical :: fix_length_b = .false.
  logical :: fix_length_c = .false.
!
! -----------------------
! fixing lattice shapes of planes
! -----------------------
  character(len("fix_shape_ab_plane")),private,parameter :: &
       &                       tag_fix_shape_ab_plane = "fix_shape_ab_plane"
  character(len("fix_shape_bc_plane")),private,parameter :: &
       &                       tag_fix_shape_bc_plane = "fix_shape_bc_plane"
  character(len("fix_shape_ac_plane")),private,parameter :: &
       &                       tag_fix_shape_ac_plane = "fix_shape_ac_plane"
  integer :: sw_fix_lattice_shapes = OFF
  logical :: fix_shape_ab_plane = .false.
  logical :: fix_shape_bc_plane = .false.
  logical :: fix_shape_ac_plane = .false.
!
! === 2014/11/22
  character(len("sw_neglect_stress_offdiagonal")),private,parameter :: &
       &                 tag_sw_neglect_stress_offdiag = "sw_neglect_stress_offdiagonal"
  integer :: sw_neglect_stress_offdiagonal = OFF
! ==== 2014/11/22

! ==== EXP_CELLOPT === 2015/09/24
! -----------------------
! read nfchgt.data of previous cell
! -----------------------
  character(len("sw_read_nfchgt_prev_cell")),private,parameter :: &
       &                 tag_sw_read_nfchgt_prev_cell = "sw_read_nfchgt_prev_cell"
  character(len("sw_read_nfzaj_prev_cell")),private,parameter :: &
       &                 tag_sw_read_nfzaj_prev_cell = "sw_read_nfzaj_prev_cell"

  integer :: sw_read_nfchgt_prev_cell = OFF
  integer :: sw_read_nfzaj_prev_cell = OFF
! ==================== 2015/09/24

! -------------------
! symmetry during optimization
! ------------------
  character(len("keep_symmetry_strict")),private,parameter :: &
       &                      tag_keep_symmetry_strict = "keep_symmetry_strict"
  integer :: sw_keep_symmetry_strict = OFF
!
! ===================================================== 13.0B

  character(len("sw_interpolate_charge")),private,parameter :: &
       &                      tag_sw_interpolate_charge = "sw_interpolate_charge"
  character(len("sw_interpolate_wfs")),private,parameter :: &
       &                      tag_sw_interpolate_wfs = "sw_interpolate_wfs"
  character(len("interpolation_method_charge")),private,parameter :: &
       &                      tag_interpolation_method_charge = "interpolation_method_charge"
  character(len("interpolation_method_wfs")),private,parameter :: &
       &                      tag_interpolation_method_wfs = "interpolation_method_wfs"
  integer :: sw_interpolate_charge = ON
  integer :: sw_interpolate_wfs = ON
  integer :: interpolation_method_chg = LINEAR_INTERPOLATION
  integer :: interpolation_method_wfs = LINEAR_INTERPOLATION

! -------------------
! charge symmetrization
! ------------------
  character(len("charge_symmetrization")), private, parameter :: &
       &             tag_charge_symmetrization = "charge_symmetrization"
  character(len("charge_symm_mode")), private, parameter :: &
       &             tag_charge_symm_mode = "charge_symm_mode"
  integer :: charge_symm_mode = 0

  ! --- mpi partitioning concerning about nbmx ---
  integer            :: ngnode_nbmx = 8
  integer            :: ncritical_vectorlength_nbmx = 100
  logical            :: flag_mpi_g_dot_r = .false.
  logical            :: flag_mpi_g_dot_r_k = .false.

  ! --- control ---
  character(len("control")),private,parameter ::      tag_control     = "control"
  character(len("condition")),private,parameter ::    tag_condition   = "condition"
  character(len("initial")),private,parameter ::      tag_initial     = "initial"
  character(len("decision_by_file_existence")),private,parameter :: tag_decision_by_file_existence = "decision_by_file_existence"
  character(len("automatic")),private,parameter ::    tag_automatic   = "automatic"
  character(len("fixed_charge_automatic")),private,parameter :: tag_fixed_charge_automatic = "fixed_charge_automatic"
  character(len("preparation")),private,parameter ::    tag_preparation   = "preparation"
  character(len("continuation")),private,parameter :: tag_continuation = "continuation"
  character(len("coordinate_continuation")),private,parameter :: tag_coordinate_continuation = "coordinate_continuation"
  character(len("continuation_using_ppdata")),private,parameter :: tag_continuation_using_ppdata = "continuation_using_ppdata"
  character(len("fixed_charge")),private,parameter :: tag_fixed_charge = "fixed_charge"
  character(len("fixed_charge_continuation")),private,parameter :: &
       &   tag_fixed_charge_continuation = "fixed_charge_continuation"
  character(len("precision_WFfile")),private,parameter:: tag_precision_WFfile = "precision_WFfile"
  character(len("double_precision")),private,parameter:: tag_double_precision = "double_precision"
  character(len("single_precision")),private,parameter:: tag_single_precision = "single_precision"
! --------------------T. Hamada 2021.9.22 -------------------- ! UVSOR
  character(len("null")), private, parameter          :: tag_null = "null"
! ------------------------------------------------------------ ! UVSOR
  character(len("DP")),private,parameter              :: tag_DP = "dp"
  character(len("SP")),private,parameter              :: tag_SP = "sp"
! -------------------- T. Hamada 2021.9.23 ------------------- ! UVSOR
  character(len("NL")),private,parameter              :: tag_NL = "nl"
  character(len("save")),private,parameter            :: tag_save = "save"
  character(len("no_save")),private,parameter         :: tag_no_save = "no_save"
  character(len("read")),private,parameter            :: tag_read = "read"
  character(len("prep")),private,parameter            :: tag_prep =  "prep"
! ------------------------------------------------------------ ! UVSOR
  character(len("fixed_charge_option")),private,parameter :: tag_fixed_charge_option = "fixed_charge_option"
  character(len("kparallel")),private,parameter ::    tag_kparallel  = "kparallel"
  character(len("one_by_one")),private,parameter ::   tag_one_by_one = "one_by_one"

  character(len("cpumax")),private,parameter ::       tag_cpumax      = "cpumax"
  character(len("max_iteration")),private,parameter ::tag_max_iteration = "max_iteration"
  character(len("max_total_scf_iteration")),private,parameter :: tag_max_total_scf_iteration = "max_total_scf_iteration"
  character(len("max_scf_iteration")),private,parameter :: tag_max_scf_iteration = "max_scf_iteration"
  character(len("max_scdft_iteration")),private,parameter ::tag_max_scdft_iteration = "max_scdft_iteration"
  character(len("max_mdstep")),private, parameter ::  tag_max_mdstep = "max_mdstep"
  character(len("nfstopcheck")),private,parameter ::  tag_nfstopcheck = "nfstopcheck"
  character(len("mpifft")),private,parameter ::       tag_mpifft      = "mpifft"
  character(len("ldx")),private,parameter ::          tag_ldx         = "ldx"
  character(len("ldy")),private,parameter ::          tag_ldy         = "ldy"
  character(len("ldz")),private,parameter ::          tag_ldz         = "ldz"
  character(len("cachesize")),private,parameter ::    tag_cachesize   = "cachesize"
#ifdef SAVE_FFT_TIMES
  character(len("sw_save_fft")),private,parameter ::  tag_sw_save_fft = "sw_save_fft"
#endif
  character(len("sw_use_wfred")),private,parameter ::        tag_sw_use_wfred      = "sw_use_wfred"

! --> T. Yamasaki, 26th Aug. 2009
!     modified by T.Kokubo & D.Fukata, Feb. 2010
  character(len("number_of_blocksize")),private,parameter :: tag_number_of_blocksize = "number_of_blocksize"
  character(len("blocksize")),private,parameter ::           tag_blocksize           = "blocksize"
  character(len("nblocksize")),private,parameter ::          tag_nblocksize          = "nblocksize"
  character(len("nb")),private,parameter ::                  tag_NB                  = "nb"
  character(len("nblocksize_dgemm")),private,parameter ::    tag_nblocksize_dgemm    = "nblocksize_dgemm"
  character(len("nblocksize_mgs")),private,parameter ::      tag_nblocksize_mgs      = "nblocksize_mgs"
  character(len("recursivesize_mgs")),private,parameter ::      tag_recursivesize_mgs      = "recursivesize_mgs"
  character(len("nblocksize_betar_dot_wfs")),private,parameter ::tag_nblocksize_betar_dot_wfs = "nblocksize_betar_dot_wfs"
  character(len("nblocksize_betar_dot_wfs_npe")),private,parameter ::tag_nblocksize_betar_dot_wfs_npe = &
       &          "nblocksize_betar_dot_wfs_npe"
  character(len("nblocksize_betar_dot_wfs_nlmta")),private,parameter ::tag_nblocksize_betar_dot_wfs_nlmta = &
       &          "nblocksize_betar_dot_wfs_nlmta"
  character(len("nblocksize_vnonlocal_w")),private,parameter :: tag_nblocksize_vnonlocal_w = "nblocksize_vnonlocal_w"
  character(len("nblocksize_betar")),private,parameter ::    tag_nblocksize_betar    = "nblocksize_betar"
  character(len("nblocksize_vnonlocal")),private,parameter :: tag_nblocksize_vnonlocal = "nblocksize_vnonlocal"
  character(len("nblocksize_submat")),private,parameter    :: tag_nblocksize_submat    = "nblocksize_submat"
  character(len("nblocksize_submat_latter")),private,parameter :: tag_nblocksize_submat_latter = "nblocksize_submat_latter"
! < --
#ifndef NO_FORCE_DGEMM
  character(len("nblocksize_force")),private,parameter     :: tag_nblocksize_force     = "nblocksize_force"
#endif
  character(len("nblocksize_rspace_betar")),private,parameter :: tag_nblocksize_rspace_betar = "nblocksize_rspace_betar"
  character(len("nblocksize_rspace_v")),private,parameter :: tag_nblocksize_rspace_v = "nblocksize_rspace_v"
  character(len("ek")),private,parameter ::           tag_ek          = "ek"
  character(len("use_intermediatefile")),private,parameter :: tag_use_intermediatefile = "use_intermediatefile"
  character(len("sw_ekzaj")),private,parameter ::           tag_sw_ekzaj = "sw_ekzaj"
  character(len("use_additional_projector")),private,parameter :: tag_use_additional_projector = "use_additional_projector"
  character(len("kpoint_data")), private,parameter :: tag_kpoint_data = "kpoint_data"
  character(len("multiple_replica_mode")),private,parameter :: tag_multiple_replica_mode = "multiple_replica_mode"
  character(len("multiple_replica_method")),private,parameter :: tag_multiple_replica_method = "multiple_replica_method"
  character(len("multiple_replica_max_iteration")),private,parameter :: tag_multiple_replica_max_iter = &
 & "multiple_replica_max_iteration"

! === KT_add === 2014/07/20
  character(len("reuse_nfout_for_nfneb")),private,parameter :: &
         tag_reuse_nfout_for_nfneb = "reuse_nfout_for_nfneb"
  integer :: reuse_nfout_for_nfneb = off
! ============== 2014/07/20

  character(len("nblocksize_vnonlocal_w_nlmta")),private,parameter :: tag_nblocksize_vnonlocal_w_nlmta = &
 & "nblocksize_vnonlocal_w_nlmta"

  character(len("nblocksize_fftw")),private,parameter :: tag_nblocksize_fftw = "nblocksize_fftw"
  character(len("nblocksize_vnonlocal_w_f")),private,parameter :: tag_nblocksize_vnonlocal_w_f = "nblocksize_vnonlocal_w_f"
  character(len("nblocksize_gather_f")),private,parameter :: tag_nblocksize_gather_f = "nblocksize_gather_f"
  character(len("fftbox_divide_cube")),private,parameter :: tag_fftbox_divide_cube = "fftbox_divide_cube"
  character(len("divide_square")),private,parameter :: tag_divide_square = "divide_square"
  character(len("fftbox_3ddiv_1")),private,parameter :: tag_fftbox_3ddiv_1 = "fftbox_3ddiv_1"
  character(len("fftbox_3ddiv_2")),private,parameter :: tag_fftbox_3ddiv_2 = "fftbox_3ddiv_2"
  character(len("fftbox_3ddiv_3")),private,parameter :: tag_fftbox_3ddiv_3 = "fftbox_3ddiv_3"
  character(len("fftbox_div_1")),private,parameter :: tag_fftbox_div_1 = "fftbox_div_1"
  character(len("fftbox_div_2")),private,parameter :: tag_fftbox_div_2 = "fftbox_div_2"
  character(len("sw_fft_xzy")),private,parameter :: tag_sw_fft_xzy = "sw_fft_xzy"

  character(len("cntn_bin_paw_format")), private, parameter :: &
       tag_cntn_bin_paw_format = "cntn_bin_paw_format"
  integer :: cntn_bin_paw_format = 0
  logical :: cntn_bin_paw_format_is_set = .false.

  character(len("sw_write_zaj")),private,parameter :: &
         tag_sw_write_zaj = "sw_write_zaj"
  integer :: sw_write_zaj = on
  character(len("sw_write_zaj_socsv")),private,parameter :: &
         tag_sw_write_zaj_socsv = "sw_write_zaj_socsv"
  integer :: sw_write_zaj_socsv = off

  ! --- structure ---
  character(len("structure")),private,parameter ::        tag_structure = "structure"

  ! --- accuracy ---
  character(len("accuracy")),parameter ::                 tag_accuracy = "accuracy"
  character(len("cke_wavefunctions")),private, parameter::tag_cke_wavefunctions = "cke_wavefunctions"
  character(len("cke_chargedensity")),private, parameter::tag_cke_chargedensity = "cke_chargedensity"
  character(len("cke_wf")),private, parameter ::          tag_cke_wf = "cke_wf"
  character(len("cke_cd")),private, parameter ::          tag_cke_cd = "cke_cd"
  character(*),private,parameter :: tag_cke_wf2 = "cutoff_energy_for_wavefunctions"
  character(*),private,parameter :: tag_cke_cd2 = "cutoff_energy_for_chargedensity"
  character(len("cutoff_wf")),private,parameter ::        tag_cke_wf3 = "cutoff_wf"
  character(len("cutoff_cd")),private,parameter ::        tag_cke_cd3 = "cutoff_cd"

  character(len("cke_pwf")),private,parameter ::          tag_cke_pwf    = "cke_pwf"
  character(len("cutoff_pwf")),private,parameter ::       tag_cutoff_pwf = "cutoff_pwf"

  character(len("num_bands")),private,parameter ::        tag_num_bands = "num_bands"
  character(len("num_extra_bands")),private,parameter ::  tag_num_extra_bands = "num_extra_bands"

! ------------------------
! occupations
!
  character(len("occupations")),private,parameter ::   &
       &                tag_occupations = "occupations"
  character(len("fixed")),private,parameter ::   &
       &                tag_occup_fixed = "fixed"
  integer :: occupations = SMEARING

! ------------------------
! smearing
!
  character(len("smearing")),private,parameter ::         tag_smearing  = "smearing"
  character(len("parabolic")),private,parameter ::        tag_parabolic = "parabolic"
  character(len("cold")),private,parameter ::             tag_cold = "cold"

! =============== KT_add =================================== 13.0E
  character(len("fermi_dirac")),private,parameter ::      tag_fermi_dirac = "fermi_dirac"
  character(len("electronic_temp")),private,parameter ::  &
       &            tag_electronic_temp = "electronic_temp"
! ========================================================== 13.0E

!  character(len("mp")),private,parameter ::               tag_mp        = "mp"
  character(len("meth")),private,parameter ::                tag_meth       = "meth"
  character(len("methfessel_paxton")),private,parameter :: tag_methfessel_paxton =  &
  & "methfessel_paxton"
  character(len("esearch_factor")),private,parameter :: tag_esearch_factor = "esearch_factor"
  character(len("esearch")),private,parameter :: tag_esearch = "esearch"
  character(len("tetrahedron")),private,parameter ::      tag_tetrahedron = "tetrahedron"
  character(len("tetrahedral")),private,parameter ::      tag_tetrahedral = "tetrahedral"
  character(len("improved_tetrahedron")),private,parameter :: tag_improved_tetrahedron = "improved_tetrahedron"
  character(len("method")),private,parameter ::           tag_smearing_method = "method"
  character(len("width")),private,parameter ::            tag_smearing_width = "width"
  character(len("dimension")),private,parameter ::        tag_dimension = "dimension"
  character(len("sw_correction")),private,parameter ::    tag_sw_correction = "sw_correction"

  character(len("lowest_at_each_kpt")),private,parameter ::  &
       &                tag_lowest_at_each_kpt = "lowest_at_each_kpt"

! --> T. Yamasaki, 2011/03/01
  character(len("fftsize")),private,parameter ::          tag_fftsize = "fftsize"
  character(len("factor_for_chargedensity")),private,parameter :: tag_factor_for_chargedensity = "factor_for_chargedensity"
  character(len("factor_for_wavefunctions")),private,parameter :: tag_factor_for_wavefunctions = "factor_for_wavefunctions"
  character(len("tag_full")),private,parameter ::         tag_full = "full"
  character(len("tag_large")),private,parameter ::        tag_large = "large"
  character(len("tag_small")),private,parameter ::        tag_small = "small"
  character(len("tag_tiny")),private,parameter ::         tag_tiny = "tiny"
! <--
  character(len("xctype")),private,parameter ::           tag_xctype    = "xctype"
  character(len("ggacmp_parallel")),private,parameter ::  tag_ggacmp_parallel = "ggacmp_parallel"
  character(len("scf_convergence")),private,parameter ::  tag_scf_convergence = "scf_convergence"
  character(len("delta_total_energy")),private,parameter ::tag_delta_total_energy = "delta_total_energy"
  character(len("delta_total_energy_initial")),private,parameter :: &
   & tag_delta_total_energy_initial = "delta_total_energy_initial"

  character(len("delta_total_energy_sampling")),private,parameter :: &
   &            tag_delta_total_energy_sampling = "delta_total_energy_sampling"
! --> T. Yamasaki, 25 July 2008
  character(len("sub_delta_factor")),private,parameter :: tag_sub_delta_factor = "sub_delta_factor"
  character(len("sub_succession")),private,parameter ::   tag_sub_succession = "sub_succession"
  character(len("sw_fix")),private,parameter ::           tag_sw_fix = "sw_fix"
! <--
  character(len("delta_eigenvalue_conduction")),private,parameter &
       &                               :: tag_delta_eigenvalue_conduction = "delta_eigenvalue_conduction"
  character(len("num_fix_bands")),private,parameter ::    tag_num_fix_bands = "num_fix_bands"
  character(len("max_force_for_edelta_initial")),private,parameter :: tag_max_force_edelta_i = "max_force_for_edelta_initial"
  character(len("edelta")),private,parameter ::           tag_edelta         = "edelta"
  character(len("edelta_initial")),private,parameter ::   tag_edelta_initial = "edelta_initial"
  character(len("edelta_ontheway")),private,parameter ::  tag_edelta_ontheway = "edelta_ontheway"
  character(len("neg")),private,parameter ::              tag_neg  = "neg"
  character(len("succession")),private,parameter ::       tag_succession = "succession"
!!$  character(len("sw_eval_eig_diff_extraband")),private,parameter ::tag_sw_eval_eig_diff_scf = "sw_eval_eig_diff_extraband"
  character(len("force_convergence")),private,parameter:: tag_force_convergence = "force_convergence"
  character(len("delta_force")),private,parameter ::      tag_delta_force = "delta_force"
  character(len("max_force")),private,parameter ::        tag_max_force = "max_force"
  character(len("tolerable_force_error")),private,parameter :: tag_tolerable_force_error = "tolerable_force_error"
  character(len("tolerable_error")),private,parameter ::  tag_tolerable_error = "tolerable_error"
  character(len("tolerable_norm_error")),private,parameter :: tag_tolerable_norm_error = "tolerable_norm_error"
  character(len("tolerable_force_norm_error")),private,parameter :: tag_tolerable_force_norm_error &
       &                                                            = "tolerable_force_norm_error"
  character(len("tolerable_angle_error")),private,parameter:: tag_tolerable_angle_error = "tolerable_angle_error"
  character(len("tolerable_hyper_angle_error")),private,parameter:: tag_tolerable_hyper_angle_error &
       &                                                            = "tolerable_hyper_angle_error"
  character(len("force_error_check_rangeL")),private,parameter :: tag_force_error_check_rangeL = "force_error_check_rangeL"
  character(len("error_check_rangeL")),private,parameter::tag_error_check_rangeL = "error_check_rangeL"

  character(len("ek_convergence")),private,parameter ::   tag_ek_convergence = "ek_convergence"
  character(len("num_max_iteration")),private,parameter ::tag_num_max_iteration = "num_max_iteration"
  character(len("sw_eval_eig_diff")),private,parameter :: tag_sw_eval_eig_diff = "sw_eval_eig_diff"
  character(len("delta_eigenvalue")),private,parameter :: tag_delta_eigenvalue = "delta_eigenvalue"
!!$  character(len("wf_inheritance")),private,parameter ::   tag_wf_inheritance = "wf_inheritance"

! --------------------------
! initial wavefunctions
! --------------------------
  character(len("initial_wavefunctions")),private,parameter :: tag_initial_wavefunctions = "initial_wavefunctions"
  character(len("matrix_diagon")),private,parameter ::    tag_matrix_diagon = "matrix_diagon"
  character(len("random_numbers")),private,parameter ::   tag_random_numbers = "random_numbers"

  character(len("pao_setting")),private,parameter :: &
       &                       tag_pao_setting = "pao_setting"
  character(len("sw_exclude_orb_tau_neq_1")),private,parameter :: &
       &                tag_sw_exclude_orb_tau_neq_1 = "sw_exclude_orb_tau_neq_1"
  character(len("sw_exclude_orb_unoccupied")),private,parameter :: &
       &                tag_sw_exclude_orb_unoccupied = "sw_exclude_orb_unoccupied"
  integer :: sw_exclude_orb_tau_neq_1 = yes
  integer :: sw_exclude_orb_unoccupied = yes

  character(len("noise_mode")),private,parameter :: tag_noise_mode = "noise_mode"
  character(len("noise_amplitude")),private,parameter :: &
       &                      tag_noise_amplitude = "noise_amplitude"
  integer :: noise_mode = 1
  real(kind=DP) :: noise_amplitude = 0.2D0
! ===============================

  character(len("cke_initial_matdiagon")),private,parameter::   tag_cke_initial_matdiagon = "cke_initial_matdiagon"
  character(len("msz_initial_matdiagon")),private,parameter ::  tag_msz_initial_matdiagon = "msz_initial_matdiagon"

! --------------------------
! initial charge density
! --------------------------
  character(len("initial_charge_density")),private,parameter :: tag_initial_charge_density = "initial_charge_density"
! --> T. Yamasaki, 28 July 2008
  character(len("initial_charge_density_file")),private,parameter :: tag_initial_charge_density_file = "initial_charge_density_file"
! <--
  character(len("Gauss_distrib_func_over_r")),private,parameter :: tag_Gauss_distrib_func_over_r = "gauss_distrib_func_over_r"
  character(len("from_PseudoPotential_FILE")),private,parameter :: tag_from_PseudoPotential_File = "from_PseudoPotential_FILE"
  character(len("Atomic_charge_density")),private,parameter :: tag_Atomic_charge_density = "Atomic_charge_density"
  character(len("very_narrow")),private,parameter ::      tag_very_narrow = "very_narrow"
  character(len("Gauss")),private,parameter ::                  tag_Gauss = "Gauss"
  character(len("given_by_a_file")),private,parameter ::  tag_given_by_a_file = "given_by_a_file"
  character(len("file")),private,parameter ::             tag_file = "file"

  character(len("off")),private,parameter ::              tag_off = "off"
  character(len("none")),private,parameter ::             tag_none = "none"
  character(len("0")),private,parameter ::              tag_0 = "0"
  character(len("unit_matrix")),private,parameter :: tag_unit_matrix = "unit_matrix"
  character(len("spin_polarized")),private,parameter :: tag_spin_polarized = "spin_polarized"
! === For restart lm+MSD! by tkato 2012/02/16 ==================================
  character(len("dtim_previous")), private, parameter :: tag_dtim_previous = "dtim_previous"
! ==============================================================================

! =============================== added by K. Tagami ================== 11.0
!
!  Ports from collinear to noncollinear
!
  character(len("ports")), private,parameter ::  &
       &        tag_ports = "ports"
  character(len("import_collinear_spindensity")), private,parameter ::  &
       &        tag_import_collinear_spindens = "import_collinear_spindensity"
  character(len("import_collinear_wavefunctions")), private,parameter ::  &
       &        tag_import_collinear_wfns = "import_collinear_wavefunctions"
  character(len("tag_previous_nspin")), private,parameter ::  &
       &        tag_previous_nspin = "previous_nspin"
  character(len("tag_previous_nband")), private,parameter ::  &
       &        tag_previous_nband = "previous_nband"

  integer :: previous_nband_collinear = 0
  integer :: previous_nspin_collinear = 2
  integer :: import_collinear_spindensity = off
  integer :: import_collinear_wavefunctions = off
! ===================================================================== 11.0

! ----
  character(len("read_charge_hardpart")), private, parameter :: &
                 & tag_read_charge_hardpart = "read_charge_hardpart"
  integer :: read_charge_hardpart = YES
! ----

  character(len("sw_add_qex_to_initial_charge")), private, parameter :: &
      & tag_sw_add_qex_to_initial_charge = "sw_add_qex_to_initial_charge"
  integer :: sw_add_qex_to_initial_charge = ON

! --> T. Yamasaki, 21st May 2010
!  integer :: PAW_switch = ON    ! OFF = 0, ON = 1
  integer :: PAW_switch = OFF    ! OFF = 0, ON = 1
  integer :: interpolation_method = LAGRANGE ! NEWTON | LAGRANGE
  integer :: polynomial_order = 4
  type paw_gradient_type
     integer :: interpolation_method = LAGRANGE
     integer :: order = 4
  end type paw_gradient_type
  type(paw_gradient_type),save :: paw_density_gradient

  character(len("PAW")),private,parameter :: tag_PAW = "PAW"
  character(len("PAW_switch")),private,parameter :: tag_PAW_switch = "PAW_switch"
  character(len("paw_density_gradient")),private,parameter :: tag_paw_density_gradient = "paw_density_gradient"
  character(len("interpolation_method")),private,parameter :: tag_interpolation_method = "interpolation_method"
  character(len("newton")),private,parameter :: tag_newton = "newton"
  character(len("lagrange")),private,parameter :: tag_lagrange = "lagrange"
  character(len("order")),private,parameter :: tag_order = "order"
  character(len=3), private, dimension(0:1) :: on_or_off = (/"OFF","ON "/)  ! 0 = OFF, 1 = ON
! <--

  integer,private ::    nitersub = 0

  ! --- Temperature Control ---
!!$  integer            :: nrsv = 1  ! number of heat bath
!!$  logical ::            tag_T_cntrl_is_found
!!$  character(len("Temperature Control")),private :: tag_T_cntrl = "Temperature Control"
!!$  character(len("temperature control")),private :: tag_T_cntrl2 = "temperature control"
!!$  character(len("nrsv")),private ::                tag_nrsv = "nrsv"


  ! --- m_ES_WF_by_MatDiagon ---
  real(kind=DP) ::      eps_solve_Hx_eq_ex = 1.d-15
  character(len("matdiagon")),private,parameter ::  tag_matdiagon        = "matdiagon"
  character(len("eps_solve_Hx_eq_ex")),private,parameter:: tag_eps_solve_Hx_eq_ex = "eps_solve_Hx_eq_ex"

  ! --- MD Algorithm ---
!!  integer ::            imdalg = QUENCHED_MD
  integer ::            imdalg = BFGS
  real(kind=DP) ::      dtio=100.d0
  character(*),private,parameter ::             tag_structure_evolution = "structure_evolution"
  character(len("dt")),private,parameter ::     tag_dt  = "dt"
  character(len("quench")),private,parameter :: tag_quench = "quench"
  character(len("gdiis")),private,parameter ::      tag_gdiis             = "gdiis"
  character(len("velocity_verlet")),private,parameter :: tag_velocity_verlet = "velocity_verlet"
  character(len("verlet")),private,parameter ::  tag_verlet = "verlet"
  character(len("quench_with_constraint")),private,parameter :: tag_quench_with_constraint = "quench_with_constraint"
  character(len("fire")),private,parameter ::  tag_fire = "fire"

  character(len("temperature_control")),private,parameter :: tag_temperature_control = "temperature_control"
  character(len("t_control")),private,parameter :: tag_t_control = "t_control"
  character(len("temp_control")),private,parameter :: tag_temp_control = "temp_control"
  character(len("bluemoon")),private,parameter :: tag_bluemoon = "bluemoon"
  character(len("steepest_descent")),private,parameter :: tag_steepest_descent = "steepest_descent"
  character(len("temperature_pressure_control")), private, parameter :: tag_temperature_pressure_control = &
  & "temperature_pressure_control"
  character(len("temperature_pressure_control")), private, parameter :: tag_pressure_temperature_control = &
  & "pressure_temperature_control"
  character(len("pressure_control")), private, parameter :: tag_pressure_control = &
  & "pressure_control"

!!$ASASASASAS
  character(len("damp")),private,parameter :: tag_damp = "damp"

  character(len("velocity_scaling")),private,parameter :: tag_velocity_scaling = "velocity_scaling"
!!$ASASASASAS
! === Include PAW by tkato =====================================================
  character(len("bfgs")),private, parameter :: tag_bfgs = "bfgs"
  character(len("sw_correct_eigenvalue")),private,parameter :: tag_sw_correct_eigenvalue = "sw_correct_eigenvalue"
  character(len("sw_optimize_alpha")), private, parameter :: tag_sw_optmize_alpha = "sw_optimize_alpha"
  character(len("lbfgs")),private, parameter :: tag_lbfgs = "lbfgs"
!!  integer :: sw_correct_eigenvalue=OFF
  integer :: sw_correct_eigenvalue=ON
  real(DP) :: eigenvalue_threshold=0.01d0
  integer :: sw_optimize_alpha = OFF
! ==============================================================================

  character(len("mu")),private,parameter :: tag_mu = "mu"
  character(len("sw_prec")),private,parameter :: tag_sw_prec = "sw_prec"
  character(len("maxstep")),private,parameter :: tag_maxstep = "maxstep"
  integer :: sw_prec = OFF
  real(kind=DP) :: precon_mu = -1.d0
  real(kind=DP) :: precon_A = 3.d0
  real(kind=DP) :: maxstep = 0.1d0

!=== FIRE optimizer
  real(kind=DP)      :: fire_incre_factor = 1.2d0
  real(kind=DP)      :: fire_decre_factor = 1.d0/1.2d0
  real(kind=DP)      :: fire_decre_factor_alpha = 1.d0/1.2d0
  integer            :: fire_nmin = 3
  real(kind=DP)      :: fire_dtmax = 300.d0
  real(kind=DP)      :: fire_initial_dt = 100.d0
  real(kind=DP)      :: fire_invmass_factor = 2.d-5
  real(kind=DP)      :: fire_alpha_start = 1.d0

  character(len("incre_factor")),private,parameter :: tag_fire_incre_factor = "incre_factor"
  character(len("decre_factor")),private,parameter :: tag_fire_decre_factor = "decre_factor"
  character(len("decre_factor_alpha")),private,parameter :: tag_fire_decre_factor_alpha = "decre_factor_alpha"
  character(len("alpha_start")),private,parameter :: tag_fire_alpha_start="alpha_start"
  character(len("nmin")),private,parameter :: tag_fire_nmin = "nmin"
  character(len("dtmax")),private,parameter :: tag_fire_dtmax = "dtmax"
  character(len("initial_dt")),private,parameter :: tag_fire_initial_dt = "initial_dt"
  character(len("invmass_factor")),private,parameter :: tag_fire_invmass_factor = "invmass_factor"

  integer            :: istress = 0
  integer,private    :: iconstpw = 0
  real(kind=DP)      :: ke_e0 = 12.5d0
  real(kind=DP)      :: ke_sigma = 0.05d0
  real(kind=DP)      :: ke_a = 0.d0
  integer            :: sw_smear_ke = OFF
  character(len("stress")),private,parameter ::   tag_stress = "stress"
  character(len("sw_stress")),private,parameter ::tag_sw_stress = "sw_stress"
  character(len("iconstpw")),private,parameter :: tag_iconstpw = "iconstpw"
  character(len("e0")),private,parameter    :: tag_e0 = "e0"
  character(len("sigma")),private,parameter :: tag_sigma = "sigma"
  character(len("a")),private,parameter :: tag_a = "a"
  character(len("sw_smear_ke")),private,parameter :: tag_sw_smear_ke = "sw_smear_ke"

!!$  character(len("sw_force_file")),private,parameter :: tag_sw_force_file = "sw_force_file"
!!$  integer ::            sw_force_file = 0

  ! --- GDIIS ---
  integer ::            kqnmditer_p = 6 ! size of gdiis box
  integer ::            gdiis_hownew = RENEW
  real(kind=DP) ::      c_forc_prop_region_high = 0.1d0
  real(kind=DP) ::      c_forc_prop_region_low  = 0.0001d0
  real(kind=DP) ::      factor_prop_region     = 0.02d0
!!  real(kind=DP) ::      c_forc2GDIIS           = 0.0025d0
  real(kind=DP) ::      c_forc2GDIIS           = 0.05d0
  integer,save ::       optmode = -999  !  QUENCHED_MD
!!$  integer,save ::       optmode = CG_STROPT !  QUENCHED_MD
  integer ::            c_iteration2GDIIS      = 3
!!  integer ::            initial_method_of_gdiis = QUENCHED_MD
!!  integer ::            initial_method_of_gdiis = CG_STROPT
  integer ::            initial_method_of_gdiis = CG_STROPT2
  character(len("gdiis_box_size")),private,parameter ::  tag_gdiis_box_size    = "gdiis_box_size"
  character(len("initial_method")),private,parameter ::  tag_initial_method    = "initial_method"
  character(*),private,parameter ::              tag_gdiis_hownew      = "gdiis_hownew"
  character(*),private,parameter ::              tag_gdiis_update      = "gdiis_update"
  character(*),private,parameter ::              tag_c_forc_prop_region_high = "c_forc_prop_region_high"
  character(*),private,parameter ::              tag_c_forc_prop_region_low = "c_forc_prop_region_low"
  character(*),private,parameter ::              tag_factor_prop_region = "factor_prop_region"
  character(*),private,parameter ::              tag_c_forc2gdiis      = "c_forc2gdiis"
  character(*),private,parameter ::              tag_c_iteration2GDIIS = "c_iteration2GDIIS"

  ! --- SD ---
  integer, public  ::      mode_fi_coefficient = 0
  real(kind=DP),public ::  fi_coefficient = 1.d0

  character(len("coefficient")),private,parameter :: tag_coefficient = "coefficient"
  character(len("mode_coefficient")),private,parameter :: tag_mode_coefficient = "mode_coefficient"

  ! --- CG ---
  real(kind=DP) ::      etol    = 1.d-3
  character(*),private,parameter ::              tag_etol = "etol"

  ! --- PRINTOUT LEVEL ---
  character(len("printoutlevel")),private,parameter ::  tag_printoutlevel     = "printoutlevel"
  character(len("printlevel")),private,parameter ::     tag_printlevel        = "printlevel"
  character(len("ipribase")),private,parameter ::       tag_ipribase          = "ipribase"
  character(len("ipritiming")),private,parameter ::     tag_ipritiming        = "ipritiming"
  character(len("timing_option")),private,parameter ::  tag_timing_option     = "timing_option"
  character(len("sw_timing_2ndlevel")),private,parameter :: tag_sw_timing_2ndlevel = "sw_timing_2ndlevel"
  character(len("sw_flatten")),private,parameter ::    tag_sw_flatten        = "sw_flatten"
  character(len("sw_firstlevel_only")),private,parameter :: tag_sw_firstlevel_only  = "sw_firstlevel_only"
  character(len("sw_details")),private,parameter ::    tag_sw_details        = "sw_details"
  character(len("measure_count_limit")),private,parameter :: tag_measure_count_limit = "measure_count_limit"

  character(len("num_subroutines")),private,parameter ::tag_num_subroutines   = "num_subroutines"
  character(len("cputime_diff")),private,parameter ::   tag_cputime_diff      = "cputime_diff"
  character(len("statistics_in_parallel")),private,parameter :: tag_statistics_in_parallel = "statistics_in_parallel"
  character(len("iprisolver")),private,parameter ::     tag_iprisolver        = "iprisolver"
  character(len("iprievdff")),private,parameter ::      tag_iprievdff         = "iprievdff"
  character(len("iprirmm")),private,parameter ::        tag_iprirmm           = "iprirmm"
  character(len("ipridavidson")),private,parameter ::   tag_ipridavidson      = "ipridavidson"
  character(len("iprimddavidson")),private,parameter ::   tag_iprimddavidson      = "iprimddavidson"
  character(len("iprimdkosugi")),private,parameter ::   tag_iprimdkosugi      = "iprimdkosugi"
  character(len("ipriesm")),private,parameter ::   tag_ipriesm      = "ipriesm"
  character(len("iprifcp")),private,parameter ::   tag_iprifcp      = "iprifcp"
  character(len("iprivdw")),private,parameter ::   tag_iprivdw      = "iprivdw"
  character(len("iprirs")),private,parameter ::   tag_iprirs      = "iprirs"
  character(len("iprirsb")),private,parameter ::   tag_iprirsb      = "iprirsb"
  character(len("iprihubbard")),private,parameter ::   tag_iprihubbard      = "iprihubbard"
  character(len("iprisnl")),private,parameter ::        tag_iprisnl           = "iprisnl"
  character(len("ipriphig")),private,parameter ::        tag_ipriphig           = "ipriphig"
  character(len("ipripao")),private,parameter ::        tag_ipripao           = "ipripao"
  character(len("ipriberry")),private,parameter ::        tag_ipriberry           = "ipriberry"
! === KT_add === 13.1R
  character(len("iprisym")),private,parameter ::        tag_iprisym = "iprisym"
! ============== 13.1R
  character(len("ipriphonon")),private,parameter ::        tag_ipriphonon           = "ipriphonon"
  character(len("iprifef")),private,parameter ::        tag_iprifef           = "iprifef"
  character(len("iprigdiis")),private,parameter ::      tag_iprigdiis         = "iprigdiis"
  character(len("iprieigenvalue")),private,parameter :: tag_iprieigenvalue    = "iprieigenvalue"
  character(len("iprispg")),private,parameter ::        tag_iprispg           = "iprispg"
  character(len("iprikp")),private,parameter ::         tag_iprikp            = "iprikp"
  character(len("ipripulay")),private,parameter ::      tag_ipripulay         = "ipripulay"
  character(len("iprimatdiagon")),private,parameter ::  tag_iprimatdiagon     = "iprimatdiagon"
  character(len("iprivlhxcq")),private,parameter ::     tag_iprivlhxcq        = "iprivlhxcq"
  character(len("ipritotalcharge")),private,parameter:: tag_ipritotalcharge   = "ipritotalcharge"
  character(len("iprisubmat")),private,parameter ::     tag_iprisubmat        = "iprisubmat"
  character(len("ipribetar")),private,parameter ::      tag_ipribetar         = "ipribetar"
  character(len("ipripaw")),private,parameter ::        tag_ipripaw           = "ipripaw"
  character(len("iprixc")),private,parameter ::         tag_iprixc            = "iprixc"
  character(len("ipristrcfctr")),private,parameter ::   tag_ipristrcfctr      = "ipristrcfctr"
  character(len("ipriinputfile")),private,parameter ::  tag_ipriinputfile     = "ipriinputfile"
  character(len("ipriparallel")),private,parameter ::   tag_ipriparallel      = "ipriparallel"
  character(len("iprifftmap")),private,parameter ::     tag_iprifftmap        = "iprifftmap"
  character(len("iprivloc")),private,parameter ::       tag_iprivloc          = "iprivloc"
  character(len("iprimd")),private,parameter ::         tag_iprimd            = "iprimd"
  character(len("iprioccup")),private,parameter ::      tag_iprioccup         = "iprioccup"
  character(len("iprichargemixing")),private,parameter::tag_iprichargemixing  = "iprichargemixing"
  character(len("ipriekzaj")),private,parameter ::      tag_ipriekzaj         = "ipriekzaj"
  character(len("ipriforce")),private,parameter ::      tag_ipriforce         = "ipriforce"
  character(len("ipridos")),private,parameter ::        tag_ipridos           = "ipridos"
  character(len("iprinegativecharge")),private,parameter :: tag_iprinegativecharge = "iprinegativecharge"
  character(len("negativecharge_option")),private,parameter::tag_negativecharge_option = "negativecharge_option"
  character(len("max_warnings")),private,parameter ::   tag_max_warnings      = "max_warnings"
  character(len("ipripp")),private,parameter ::         tag_ipripp            = "ipripp"
  character(len("ipriwf")),private,parameter ::         tag_ipriwf            = "ipriwf"
  character(len("ipricoefwf")),private,parameter ::     tag_ipricoefwf        = "ipricoefwf"
  character(len("iprichargedensity")),private,parameter ::     tag_iprichargedensity        = "iprichargedensity"
  character(len("iprivelocity")),private,parameter ::   tag_iprivelocity      = "iprivelocity"
  character(len("ipriparadeb")),private,parameter ::    tag_ipriparadeb       = "ipriparadeb"
  character(len("ipriparallel_debug")),private,parameter :: tag_ipriparallel_debug = "ipriparallel_debug"
  character(len("iprijobstatus")),private,parameter ::  tag_iprijobstatus     = "iprijobstatus"
  character(len("ipripredictor")),private,parameter ::  tag_ipripredictor     = "ipripredictor"
  character(len("ipriunitcell")),private,parameter ::   tag_ipriunitcell       = "ipriunitcell"
  character(len("iprifire")),private,parameter ::       tag_iprifire       = "iprifire"
  character(len("ipriexx")),private,parameter ::        tag_ipriexx        = "ipriexx"
  character(len("ipriblsize")),private,parameter ::     tag_ipriblsize     = "ipriblsize"
  character(len("jobstatus_option")),private,parameter :: tag_jobstatus_option = "jobstatus_option"
  character(len("jobstatus_series")),private,parameter :: tag_jobstatus_series = "jobstatus_series"
  character(len("jobstatus_format")),private,parameter :: tag_jobstatus_format = "jobstatus_format"
  character(len("tag_line")),private,parameter ::       tag_tag_line          = "tag_line"
  character(len("tag")),private,parameter ::            tag_tag               = "tag"
  character(len("table")),private,parameter ::          tag_table             = "table"

  character(len("n_fermi_vicinity")),private,parameter :: tag_n_fermi_vicinity = "n_fermi_vicinity"

  ! --- POSTPROCESSING ---
  character(len("Postprocessing")),private,parameter :: tag_postprocessing    = "postprocessing"

  character(len("frequency")),private,parameter :: tag_frequency = "frequency"
  integer, public :: postproc_frequency=-1

#ifndef _EMPIRICAL_
  ! --- RMM ---
  integer,public ::            imGSrmm  = 1
  ! Orthonormalization is done every IMGSRMM-th iteration(in msd_rmm_n_uda).
  real(kind=DP),public ::      rr_Critical_Value = 1.d-15
  integer,public ::            rmm_printout       = 0
  integer,public ::            rmm_precal_phase_matm     = 0
  integer ::            rmm3_bisec_trial_max   = 150
  real(kind=DP) ::      rmm3_bisec_crtcl_value = 1.d-90
!! asms
!!$  real(kind=DP),private ::      edelta_change_to_rmm   = 1.d-7
  real(kind=DP),private ::      edelta_change_to_rmm   = 1.d-3
  real(kind=DP),private ::      edelta_change_to_rmm_md  = 1.d-3

  real(kind=DP) :: delta_residual_occup = 1e-15
  real(kind=DP) :: delta_residual_empty = 1e-15

  logical,private :: edelta_rmm_given = .false.
  integer,public ::            rmm_save_memory_mode = OFF
  character(len("rmm")),private,parameter ::                 tag_rmm               = "rmm"
  character(len("imGSrmm")),private,parameter ::             tag_imGSrmm           = "imgsrmm"
  character(len("rr_critical_value")),private,parameter ::   tag_rr_critical_value = "rr_critical_value"
  character(len("rmm_printout")),private,parameter ::        tag_rmm_printout      = "rmm_printout"
  character(len("rmm_precal_phase_matm")),private,parameter::tag_rmm_precal_phase_matm  = "rmm_precal_phase_matm"
  character(len("rmm3_bisec_trial_max")),private,parameter ::tag_rmm3_bisec_trial_max   = "rmm3_bisec_trial_max"
  character(len("rmm3_bisec_crtcl_value")),private,parameter::tag_rmm3_bisec_crtcl_value = "rmm3_bisec_crtcl_value"
  character(len("edelta_change_to_rmm")),private,parameter:: tag_edelta_change_to_rmm   = "edelta_change_to_rmm"
  character(len("edelta_change_to_rmm_md")),private,parameter:: tag_edelta_change_to_rmm_md   = "edelta_change_to_rmm_md"
  character(len("save_memory_mode")),private,parameter ::    tag_save_memory_mode  = "save_memory_mode"

  character(len("delta_residual_occup")),private,parameter ::   tag_delta_residual_occup = "delta_residual_occup"
  character(len("delta_residual_empty")),private,parameter ::   tag_delta_residual_empty = "delta_residual_empty"

  ! --- Line Minimization ---
  real(kind=DP) ::      delta_lmdenom     = 1.d-12
  real(kind=DP) ::      dt_Lower_CRITICAL = 5.d-2
  real(kind=DP) ::      dt_Upper_CRITICAL = 2.d0
  real(kind=DP),private :: dt_upper_factor = 5.d0
  logical,private ::       dt_upper_factor_is_set = .false.
  real(kind=DP),private :: dt_lower_factor = 0.1d0
  integer,private ::       incre_etot_in_1dsrch = 0
  integer,private ::       incre_etot_in_1dsrch_limit = 10
  real(kind=DP) ::         de_critical_1dsrch = 10.d0

  integer :: num_conduction_bands_lmm = -1
  integer :: energy_evaluation = TOTAL_ENERGY
  integer :: sw_lmm_status_check = NO

  character(len("lineminimization")),private,parameter ::   tag_lineminimization = "lineminimization"
  character(len("line_minimization")),private,parameter ::  tag_line_minimization = "line_minimization"
  character(len("dt_lower_critical")),private,parameter ::   tag_dt_lower_critical = "dt_lower_critical"
  character(len("dt_upper_critical")),private,parameter ::   tag_dt_upper_critical = "dt_upper_critical"
  character(len("dt_upper_factor")),private,parameter ::     tag_dt_upper_factor  = "dt_upper_factor"
  character(len("dt_lower_factor")),private,parameter ::     tag_dt_lower_factor  = "dt_lower_factor"
  character(len("delta_lmdenom")),private,parameter ::       tag_delta_lmdenom     = "delta_lmdenom"

  character(len("energy_evaluation")),private,parameter ::   tag_energy_evaluation = "energy_evaluation"
  character(len("num_conduction_bands")),private,parameter :: tag_num_conduction_bands = "num_conduction_bands"
  character(len("total_energy")),private,parameter      ::   tag_total_energy = "total_energy"
  character(len("modified_total_energy")),private,parameter::tag_modified_total_energy = "modified_total_energy"
  character(len("band_energy")),private,parameter ::         tag_band_energy = "band_energy"

  character(len("sw_lmm_status_check")),private,parameter :: tag_sw_lmm_status_check = "sw_lmm_status_check"
  character(len("lmm_status_check")),private,parameter :: tag_lmm_status_check = "lmm_status_check"

  ! --- Subspace rotation ---
  integer ::            ldiag, meg
  real(kind=DP) ::      damp = 1.d0
  integer ::            submat_period = 1
  real(kind=DP) ::      submat_critical_ratio = 1.d-15
#ifdef _USE_SCALAPACK_
  integer :: sw_scalapack = ON
  integer :: sw_scalapack_md = OFF
#else
  integer :: sw_scalapack = OFF
  integer :: sw_scalapack_md = OFF
#endif
#ifdef USE_EIGENLIB
  integer :: sw_eigen_exa = OFF
#endif
#ifdef _DEFAULT_HOUSEHOLDER_
  integer :: method_scalapack = HOUSEHOLDER
#elif _DEFAULT_DIVIDEandCONQUER_
  integer :: method_scalapack = DIVIDEandCONQUER
#else
  integer :: method_scalapack = HOUSEHOLDER
#endif
!finteger :: block_size = 64
  integer :: block_size = 0
  integer :: nprow = 0, npcol = 0
  integer :: msize_submat = 512 ! MegaByte
!!$  integer :: submat_before_renewal = OFF
  integer :: submat_before_renewal = ON

  integer :: sw_gep = OFF

  character(len("subspace_rotation")),private,parameter ::   tag_subspace_rotation = "subspace_rotation"
  character(len("subspace_matrix_size")),private,parameter:: tag_subspace_matrix_size = "subspace_matrix_size"
  character(len("damping_factor")),private,parameter ::      tag_damping_factor = "damping_factor"
  character(len("one_period")),private,parameter ::          tag_one_period = "one_period"
  character(len("period")),private,parameter ::              tag_period = "period"
  character(len("critical_ratio")),private,parameter ::      tag_critical_ratio = "critical_ratio"
  character(len("scalapack")),private,parameter :: tag_scalapack = "scalapack"
  character(len("sw_scalapack")),private,parameter :: tag_sw_scalapack = "sw_scalapack"
  character(len("block_size")),private,parameter :: tag_block_size = "block_size"
  character(len("nprow")),private,parameter :: tag_nprow = "nprow"
  character(len("npcol")),private,parameter :: tag_npcol = "npcol"
  character(len("householder")),private,parameter :: tag_householder = "householder"
  character(len("divide_and_conquer")),private,parameter :: tag_divide_and_conquer = "divide_and_conquer"
  character(len("memory")),private,parameter :: tag_memory = "memory"
  character(len("sw_gep")),private,parameter :: tag_sw_gep = "sw_gep"

#ifdef USE_EIGENLIB
  character(len("eigen_exa")), private, parameter :: tag_eigen_exa="eigen_exa"
  character(len("sw_eigen_exa")), private, parameter :: tag_sw_eigen_exa="sw_eigen_exa"
#endif

  ! --- CG  ---
  !!$integer :: sw_modified_cg_formula = OFF
  integer :: sw_modified_cg_formula = ON
  character(len("sw_modified_cg_formula")),private,parameter ::   tag_sw_modified_cg_formula = "sw_modified_cg_formula"

  ! --- Davidson ---
! integer :: ndavid = 20
! integer :: max_iter_david = 20
  integer :: ndavid = 5
  integer :: max_iter_david = 3
  integer :: max_subspace_size = 0
  integer :: sw_first_conv_check = ON
  real(kind=DP) :: delta_eig_occup = 1.d-8
  real(kind=DP) :: delta_eig_empty = 1.d-8
! === Default value of eps_david seems to be too large. ==================================
! === This occurs energy divergence when MDDAVIDSON is used. by T.Kato 2013/08/06 ========
! real(kind=DP) :: eps_david = 1.d-4
  real(kind=DP) :: eps_david = 1.d-8
! ========================================================================================
  integer :: sw_dav_scalapack = OFF
  integer :: dav_nprow = 0, dav_npcol = 0
  integer :: dav_divide_square = 0
  integer :: dav_block_size = 32
  integer :: nblock = 8
  character(len("davidson")),private,parameter ::     tag_davidson          = "davidson"
  character(len("ndavid")),private,parameter ::       tag_ndavid = "ndavid"
  character(len("max_iter_david")),private,parameter :: tag_max_iter_david = "max_iter_david"
  character(len("max_subspace_size")),private,parameter :: tag_max_subspace_size = "max_subspace_size"
  character(len("delta_eig_occup")),private,parameter :: tag_delta_eig_occup = "delta_eig_occup"
  character(len("delta_eig_empty")),private,parameter :: tag_delta_eig_empty = "delta_eig_empty"
  character(len("eps_david")),private,parameter :: tag_eps_david = "eps_david"
  character(len("sw_first_conv_check")),private,parameter :: tag_sw_first_conv_check = "sw_first_conv_check"
  character(len("scalapack")),private,parameter :: tag_dav_scalapack = "scalapack"
  character(len("sw_scalapack")),private,parameter :: tag_sw_dav_scalapack = "sw_scalapack"
  character(len("nprow")),private,parameter :: tag_dav_nprow = "nprow"
  character(len("npcol")),private,parameter :: tag_dav_npcol = "npcol"
  character(len("divide_square")),private,parameter :: tag_dav_divide_square = "divide_square"
  character(len("block_size")),private,parameter :: tag_dav_block_size = "block_size"

  ! modified Davidson, modified Kosugi merged
  integer :: sw_MRCV_only = ON
  integer :: npartition_david = 1
  integer :: sw_divide_subspace = ON
  integer :: submat_GE = OFF
  logical :: sw_divide_subspace_changed=.false.
  logical :: sw_npartition_changed=.false.

! md-david, md-kosugi
  character(len("sw_MRCV_only")), private, parameter :: tag_sw_MRCV_only = "sw_MRCV_only"
  character(len("sw_divide_subspace")), private, parameter :: tag_sw_divide_subspace = "sw_divide_subspace"
  character(len("submat_GE")),private,parameter :: tag_submat_GE = "submat_GE"
  character(len("npartition_david")), private, parameter :: tag_npartition_david = "npartition_david"
  character(len("npartition")), private, parameter :: tag_npartition = "npartition"
  character(len("nblock_david")), private, parameter :: tag_nblock_david = "nblock_david"
  character(len("nblock")), private, parameter :: tag_nblock = "nblock"
  character(len("nbands_in_partition")),private,parameter :: tag_nbands_in_partition = "nbands_in_partition"

  ! --- Modified Davidson ---
  integer :: npartition_mddavid = 1
  !!$integer :: max_iter_mddavid = 4
  integer :: max_iter_mddavid = 2
  real(kind=DP) :: delta_eig_occup_md = 1.d-8
  real(kind=DP) :: delta_eig_empty_md = 1.d-8
  real(kind=DP) :: eps_mddavid = 1.d-4
  real(kind=DP) :: eps_residual_mddavid = 1.d-6
  character(len("mddavidson")),private,parameter ::     tag_mddavidson          = "mddavidson"
  character(len("npartition_mddavid")),private,parameter ::tag_npartition_mddavid = &
                                                                             "npartition_mddavid"
  character(len("max_iter_mddavid")),private,parameter :: tag_max_iter_mddavid = "max_iter_mddavid"
!  character(len("delta_eig_occup")),private,parameter :: tag_delta_eig_occup = "delta_eig_occup"
!  character(len("delta_eig_empty")),private,parameter :: tag_delta_eig_empty = "delta_eig_empty"
  character(len("eps_mddavid")),private,parameter :: tag_eps_mddavid = "eps_mddavid"
  character(len("eps_residual_mddavid")),private,parameter :: tag_eps_residual_mddavid = "eps_residual_mddavid"

  character(len("p-davidson")),private,parameter :: tag_pdav1 = "p-davidson"
  character(len("p_davidson")),private,parameter :: tag_pdav2 = "p_davidson"
  character(len("pdavidson")),private,parameter :: tag_pdav3 = "pdavidson"

  character(len("p-kosugi")),private,parameter :: tag_pkos1 = "p-kosugi"
  character(len("p_kosugi")),private,parameter :: tag_pkos2 = "p_kosugi"
  character(len("pkosugi")),private,parameter :: tag_pkos3 = "pkosugi"

  ! --- Modified Kosugi ---
  integer :: npartition_mdkosugi = 1
  !!$integer :: max_iter_mdkosugi = 4
  integer :: max_iter_mdkosugi = 2
  real(kind=DP) :: delta_eig_occup_mdkosugi = 1.d-8
  real(kind=DP) :: delta_eig_empty_mdkosugi = 1.d-8
  real(kind=DP) :: eps_mdkosugi = 1.d-4
  real(kind=DP) :: eps_residual_mdkosugi = 1.d-6
  integer :: sw_apply_gs = OFF
  character(len("mdkosugi")),private,parameter ::     tag_mdkosugi          = "mdkosugi"
  character(len("npartition_mdkosugi")),private,parameter ::tag_npartition_mdkosugi = "npartition_mdkosugi"
  character(len("max_iter_mdkosugi")),private,parameter :: tag_max_iter_mdkosugi = "max_iter_mdkosugi"
  character(len("delta_eig_occup_mdkosugi")),private,parameter :: tag_delta_eig_occup_mdkosugi = "delta_eig_occup_mdkosugi"
  character(len("delta_eig_empty_mdkosugi")),private,parameter :: tag_delta_eig_empty_mdkosugi = "delta_eig_empty_mdkosugi"
  character(len("eps_mdkosugi")),private,parameter :: tag_eps_mdkosugi = "eps_mdkosugi"
  character(len("eps_residual_mdkosugi")),private,parameter :: tag_eps_residual_mdkosugi = "eps_residual_mdkosugi"
  character(len("sw_apply_gs")),private,parameter :: tag_sw_apply_gs = "sw_apply_gs"

  ! --- Control of Sovers for Wave Functions ---
  character(*),private,parameter ::                   tag_wavefunction_solver = "wavefunction_solver"
  character(len("wf_solver")),private,parameter ::    tag_wf_solver         = "wf_solver"
  character(len("for_init_str")),private,parameter :: tag_for_init_str      = "for_init_str"
  character(len("during_str_relax")),private,parameter :: tag_during_str_relax = "during_str_relax"
  character(len("num_solvers")),private,parameter ::  tag_num_solvers       = "num_solvers"
  character(len("solvers")),private,parameter ::      tag_solvers           = "solvers"
  character(len("lmmsd")),private,parameter ::        tag_lmMSD             = "lmmsd"
  character(len("lm+msd")),private,parameter ::       tag_lmMSD2            = "lm+msd"
  character(len("msd")),private,parameter ::          tag_msd               = "msd"
  character(len("lmsd")),private,parameter ::         tag_lmsd              = "lmsd"
  character(len("lm+sd")),private,parameter ::        tag_lmsd2             = "lm+sd"
  character(len("sd")),private,parameter ::           tag_sd                = "sd"
  character(len("rmm2p")),private,parameter ::        tag_rmm2p             = "rmm2p"
  character(len("rmm3")),private,parameter ::         tag_rmm3              = "rmm3"
  character(len("rmm2")),private,parameter ::         tag_rmm2              = "rmm2"
  character(len("matrixdiagon")),private,parameter :: tag_matrixdiagon      = "matrixdiagon"
  character(len("submat")),private,parameter ::       tag_submat            = "submat"
  character(len("subspacerotation")),private,parameter :: tag_subspacerotation = "subspacerotation"
  character(len("cg")),private,parameter ::           tag_cg                = "cg"
  character(len("cg2")),private,parameter ::          tag_cg2               = "cg2"
  character(len("oldcg")),private,parameter ::        tag_oldcg             = "oldcg"
  character(len("lmcg")),private,parameter ::         tag_lmcg              = "lmcg"
  character(len("lm+cg")),private,parameter ::        tag_lmcg2             = "lm+cg"

  character(len("solver_of_WF")),private,parameter :: tag_solver_of_WF      = "solver_of_WF"
  character(len("n_before")),private,parameter ::     tag_n_before          = "n_before"
  character(len("n_after")),private,parameter ::      tag_n_after           = "n_after"
  character(len("sol")),private,parameter ::          tag_sol               = "sol"
  character(len("prec")),private,parameter ::         tag_prec              = "prec"
  character(len("till_n")),private,parameter ::       tag_till_n            = "till_n"
  character(len("dts")),private,parameter ::          tag_dts               = "dts"
  character(len("dte")),private,parameter ::          tag_dte               = "dte"
  character(len("itr")),private,parameter ::          tag_itr               = "itr"
  character(len("var")),private,parameter ::          tag_var               = "var"
  character(len("linear")),private,parameter ::       tag_linear            = "linear"
  character(len("tanh")),private,parameter ::         tag_tanh              = "tanh"
  character(len("cmix")),private,parameter ::         tag_cmix              = "cmix"

  integer, private  :: solver_set =  DAV_RMM

  character(len("before_renewal")),private,parameter ::tag_before_renewal = "before_renewal"

  ! --- band_symmetry_analysis ---
  character(len("sw_band_symmetry_analysis")),private,parameter::tag_sw_band_symmetry_analysis = "sw_band_symmetry_analysis"
  character(len("sw_symmetry_analysis")),     private,parameter::tag_sw_symmetry_analysis      = "sw_symmetry_analysis"
  integer   ::                                                   sw_band_symmetry_analysis


!!!!  integer,private :: n_WF_solvers_before = 1, n_WF_solvers_after = 1
!!!!  integer,private :: n_WF_solvers_all = 2
  integer :: n_WF_solvers_before = 1, n_WF_solvers_after = 1
  integer :: n_WF_solvers_all = 2
  type solver
     integer :: before_or_after_convergence
     integer :: solver
     integer :: till_n_iter
     integer :: precon
     integer :: iter_range
     integer :: variation_way
     integer :: cmix_pointer
     integer :: subspace_rotation
     real(kind=DP) :: dtim_s, dtim_e
!!$     integer :: before_or_after_convergence = BEFORE
!!$     integer :: solver = MSD
!!$     integer :: till_n_iter = -1
!!$     integer :: precon = YES
!!$     integer :: iter_range = 100
!!$     integer :: variation_way = varLINEAR
!!$     integer :: cmix_pointer = 1
!!$     real(kind=DP) :: dtim_s = 0.1, dtim_e=0.1
  end type solver

!!!!  type(solver),private,allocatable,dimension(:) :: w_solver
  type(solver),allocatable,dimension(:) :: w_solver
  logical,public     :: tag_solver_of_WF_is_found = .false.
  integer,private    :: ip_w_solver = 1
  ! --------------

  ! --- Charge Density Mixing ---
  character(len("charge_mixing")),private,parameter:: tag_charge_mixing   = "charge_mixing"
  character(len("n_mixing_way")),private,parameter :: tag_n_mixing_way    = "n_mixing_way"
  character(len("num_mixing_methods")),private,parameter::tag_num_mixing_methods = "num_mixing_methods"
  character(len("mixing_methods")),private,parameter:: tag_mixing_methods = "mixing_methods"
  character(len("method")),private,parameter ::       tag_method          = "method"
  character(len("simple")),private,parameter ::       tag_simple          = "simple"
  character(len("broyden")),private,parameter ::      tag_broyden         = "broyden"
  character(len("broyden2")),private,parameter ::     tag_broyden2        = "broyden2"
  character(len("b2")),private,parameter ::           tag_b2              = "b2"
  character(len("dfp")),private,parameter ::          tag_dfp             = "dfp"
  character(len("pulay")),private,parameter ::        tag_pulay           = "pulay"
  character(len("no")),private,parameter ::           tag_no              = "no"
  character(len("id")),private,parameter ::           tag_id              = "id"
  character(len("hownew")),private,parameter ::       tag_hownew          = "hownew"
  character(len("update")),private,parameter ::       tag_update          = "update"
  character(len("renew")),private,parameter ::        tag_renew           = "renew"
  character(len("anew")),private,parameter ::         tag_anew            = "anew"
  character(len("cutoff")),private,parameter ::       tag_cutoff          = "cutoff"
  character(len("istr")),private,parameter ::         tag_istr            = "istr"
  character(len("nbxmix")),private,parameter ::       tag_nbxmix          = "nbxmix"
  character(len("nbmix")),private,parameter ::        tag_nbmix           = "nbmix"
  character(len("rmxs")),private,parameter ::         tag_rmxs            = "rmxs"
  character(len("rmxe")),private,parameter ::         tag_rmxe            = "rmxe"
  character(len("rmx")),private,parameter ::          tag_rmx             = "rmx"
! --> T. Yamasaki 04 Aug. 2009
  character(len("sw_recomposing")),private,parameter::tag_sw_recomposing  = "sw_recomposing"
  character(len("spin_density_mixfactor")),private,parameter :: tag_spin_density_mixfactor = "spin_density_mixfactor"
!!$  integer       :: sw_recomposing = NO
  integer       :: sw_recomposing = YES
  real(kind=DP) :: spin_density_mixfactor = 1.0
  integer, private, save :: print_sw_recomposing = 0   ! 2024.06.04 T. Yamasaki
! <-

  character(len("occ_matrix_mixfactor")),private,parameter :: tag_occ_matrix_mixfactor="occ_matrix_mixfactor"
  character(len("occ_matrix_amix")),private,parameter :: tag_occ_matrix_amix="occ_matrix_amix"
  character(len("initial_occ_matrix")),private,parameter :: tag_initial_occ_matrix="initial_occ_matrix"
  character(len("initial_es")),private,parameter :: tag_initial_es="initial_es"
  character(len("amix_hsr")), private, parameter :: tag_amix_hsr="amix_hsr"
  real(kind=DP) :: occ_matrix_mixfactor = 1.0
  real(kind=DP) :: occ_matrix_amix = -1
  real(kind=DP) :: amix_hsr = 1
  integer :: initial_occ_matrix=0
  integer :: sw_initial_es=OFF

  ! --
  character(len("charge_preconditioning")),private,parameter :: &
       & tag_charge_preconditioning = "charge_preconditioning"
  character(len("amix")),private,parameter ::   tag_amix = "amix"
  character(len("bmix")),private,parameter ::   tag_bmix = "bmix"

! ==> by J. Koga, 25 Nov 2010
  character(len("metric_ratio")),private, parameter :: tag_metric_ratio = "metric_ratio"
  character(len("amin")),private,parameter  ::   tag_amin = "amin"
  character(len("sw_precon_diff")), private, parameter :: tag_sw_precon_diff = "sw_precon_diff"
  character(len("sw_metric_diff")), private, parameter :: tag_sw_metric_diff = "sw_metric_diff"
  integer :: sw_precon_diff = NO
  integer :: sw_metric_diff = NO
  integer :: metric_ratio = -1

! === KT_add === 13.0U3
  character(len("precon_mode")),private,parameter ::   tag_precon_mode = "precon_mode"
!
  integer :: precon_mode = 1
! ============== 13.0U3

! -- spin density
  character(len("spin_density")), private, parameter :: tag_spin_density="spin_density"
  character(len("sw_apply_precon")), private, parameter :: tag_sw_apply_precon="sw_apply_precon"
  character(len("sw_apply_metric")), private, parameter :: tag_sw_apply_metric="sw_apply_metric"
  character(len("sw_force_simple_mixing")), private, parameter :: &
 & tag_sw_force_simple_mixing="sw_force_simple_mixing"
!!$  integer :: sw_force_simple_mixing=OFF    ! 2024.06.04 by T. Yamasaki
  integer :: sw_force_simple_mixing=OFF

  integer :: cmix_set = PULAY
  logical :: cmix_explicit = .true.
! ============================= added by K. Tagami ================ 5.0
  character(len("eval_energy_before_charge")), private, parameter :: &
                 & tag_eval_energy_before_charge = "eval_energy_before_charge"
  integer :: eval_energy_before_charge  = OFF
  integer :: sw_eval_energy_before_charge = OFF
!
!
  character(len("sw_mix_charge_hardpart")), private, parameter :: &
                         & tag_sw_mix_charge_hardpart = "sw_mix_charge_hardpart"
  character(len("sw_mix_bothspins_sametime")), private, parameter :: &
                         & tag_sw_mix_bothspins_sametime = "sw_mix_bothspins_sametime"
  character(len("sw_mix_occ_matrix")),private,parameter :: &
                         & tag_sw_mix_occ_matrix = "sw_mix_occ_matrix"
!!$  integer :: sw_mix_charge_hardpart = OFF
!!$  integer :: sw_mix_bothspins_sametime = OFF
  integer :: sw_mix_charge_hardpart = OFF
  integer :: sw_mix_bothspins_sametime = ON
!
!
  integer :: sw_mix_occ_matrix = OFF

  character(len("sw_force_simple_mixing_hsr")), private, parameter :: &
                 & tag_sw_force_simplemix_hsr = "sw_force_simple_mixing_hsr"
  character(len("sw_recomposing_hsr")), private, parameter :: &
                 & tag_sw_recomposing_hsr = "sw_recomposing_hsr"

  integer :: sw_force_simple_mixing_hsr = OFF
  integer :: sw_recomposing_hsr = OFF
!
  integer :: sw_mix_bothspins_sametime_hsr = OFF
! ----------
  integer :: sw_update_charge_total = ON
  integer :: sw_update_charge_hsr = ON
! ======================================================================= 5.0

! ============== added by K. Tagami ============= 11.0
  character(len("sw_mix_imaginary_hardpart")), private, parameter :: &
                         & tag_sw_mix_imaginary_hardpart = "sw_mix_imaginary_hardpart"
  integer :: sw_mix_imaginary_hardpart = OFF
! ================================================ 11.0

! ===== KT_add ==== 2014/09/19
  character(len("sw_mix_charge_with_eikndens")), private, parameter :: &
                         & tag_sw_mix_chg_with_ekindens = "sw_mix_charge_with_ekindens"
  integer :: sw_mix_charge_with_ekindens = OFF
! =================2014/09/19

! ==================================== added by K. Tagami ============== 5.0
! ------------------------------------------------------- not used
  integer :: sw_force_simple_mixing_occdiff = OFF
  integer :: sw_recomposing_occmat = OFF
! =======================================================================5.0

  character(len("sw_gradient_simplex")),private,parameter :: tag_sw_gradient_simplex="sw_gradient_simplex"
  integer :: sw_gradient_simplex = OFF
  real(kind=DP) :: alpha_pulay=0.d0
  real(kind=DP) :: alpha_pulay_org=0.d0
  real(kind=DP) :: alpha_pulay_damp=1.0d0
  real(kind=DP) :: alpha_pulay_damp_thres=1.d-4
  character(len("sw_control_stepsize")),private,parameter :: tag_sw_control_stepsize="sw_control_stepsize"
  character(len("max_stepsize")),private,parameter :: tag_max_stepsize="max_stepsize"
  character(len("threshold_alpha")),private,parameter :: tag_threshold_alpha="threshold_alpha"
  integer :: sw_control_stepsize = OFF
  real(kind=DP) :: max_stepsize = 0.4d0

  character(len("ommix_factor")),private,parameter :: tag_ommix_factor = "ommix_factor"
  real(kind=DP) :: ommix_factor = 1.d0

! ========= KT_add ========== 13.0U2
! Potential mixing
!
  character(len("sw_potential_mixing")),private,parameter :: &
       &                        tag_sw_potential_mixing = "sw_potential_mixing"
  integer :: sw_potential_mixing = off
! =========================== 13.0U2

! =============== KT_add ========================== 13.0U2
! Thomas Fermi von Weizsacker functional (for charge/potential mixing)
!
  character(len("kinetic_energy_functional")), private, parameter :: &
       &             tag_kinetic_energy_functional = "kinetic_energy_functional"
  character(len("use_TFW_functional")), private, parameter :: &
       &             tag_use_TFW_functional = "use_TFW_functional"
  character(len("spherical_averaged_nonlocal")), private, parameter :: &
       &             tag_spherical_averaged_nonlocal = "spherical_averaged_nonlocal"
  character(len("weight_ThomasFermi")), private, parameter :: &
       &             tag_weight_ThomasFermi = "weight_ThomasFermi"
  character(len("weight_Weizsacker")), private, parameter :: &
       &             tag_weight_Weizsacker = "weight_Weizsacker"
  character(len("use_averaged_nonlocal")), private, parameter :: &
       &             tag_use_averaged_nonlocal = "use_averaged_nonlocal"
  character(len("use_deltaW")), private, parameter :: &
       &             tag_use_deltaW = "use_deltaW"
  character(len("use_deltaV")), private, parameter :: &
       &             tag_use_deltaV = "use_deltaV"
!
  integer :: sw_modified_TFW_functional = OFF
  integer :: sw_spherical_averaged_nonlocal = OFF
!
  logical :: use_averaged_nonlocal = .true.
  logical :: use_deltaW = .true.
  logical :: use_deltaV = .true.
!
  real(kind=DP) :: weight_TF_functional   = 1.00d0
  real(kind=DP) :: weight_Weiz_functional = 1.00d0

! ------------
! scf iterations
! ------------
  character(len("scf_loop")), private, parameter :: &
       &             tag_scf_loop = "scf_loop"
  character(len("threshold_exit_loop")), private, parameter :: &
       &             tag_threshold_exit_loop = "threshold_exit_loop"
  character(len("threshold_skip_loop")), private, parameter :: &
       &             tag_threshold_skip_loop = "threshold_skip_loop"
  character(len("threshold_enter_loop")), private, parameter :: &
       &             tag_threshold_enter_loop = "threshold_enter_loop"

  character(len("max_iteration_loop")), private, parameter :: &
       &             tag_max_iteration_loop = "max_iteration_loop"
  character(len("use_preconditioning")), private, parameter :: &
       &             tag_use_preconditioning = "use_preconditioning"
!
  real(kind=DP) :: threshold_enter_loop = 1.0D+2
  real(kind=DP) :: threshold_exit_loop = 1.0D-6
  real(kind=DP) :: threshold_skip_loop = 1.0D-6
!
  integer :: max_iteration_loop = 20
  logical :: use_preconditioning = .false.
!
! ------------
! line minimization
! ------------
  character(len("max_iteration_linmin")), private, parameter :: &
       &             tag_max_iteration_linmin = "max_iteration_linmin"
  character(len("threshold_linmin")), private, parameter :: &
       &             tag_threshold_linmin = "threshold_linmin"
!
  integer :: max_iteration_linmin = 10
  real(kind=DP) :: threshold_linmin = 1.0D-7

! ------------
! printing
! ------------
  character(len("ipritfwfunc")),private,parameter :: tag_ipritfwfunc = "ipritfwfunc"
!
  integer :: ipritfwfunc = -1
! ================================================== 13.0U2

! ===== KT_add ========= 13.0U3
! Hard part preconditioning, mixing
!
  character(len("hardpart_mixfactor")),private,parameter :: &
       &                          tag_hardpart_mixfactor = "hardpart_mixfactor"
  real(kind=DP) :: hardpart_mixfactor = 1.0d0
!
!
  character(len("preconditioning_hardpart")),private,parameter :: &
       &                  tag_preconditioning_hardpart = "preconditioning_hardpart"
  character(len("aug_charge_density")),private,parameter :: &
       &                  tag_aug_charge_density = "aug_charge_density"
  character(len("sw_wf_mixing")),private,parameter :: &
       &                  tag_sw_wf_mixing = "sw_wf_mixing"
  character(len("recover_wf_after_mixing")),private,parameter :: &
       &        tag_recover_wf_after_mixing = "recover_wf_after_mixing"
  character(len("edelta_start_wf_mixing")),private,parameter :: &
       &                  tag_edelta_start_wf_mixing = "edelta_start_wf_mixing"

  character(len("amin")),private,parameter ::  tag_amin_wf_precon = "amin"
  character(len("bmix")),private,parameter ::  tag_bmix_wf_precon = "bmix"
  character(len("g0")),  private,parameter ::  tag_g0_wf_precon = "g0"

  logical :: precon_hardpart = .false.
  integer :: sw_wf_mixing = off
  logical :: wf_mixing_is_active = .true.
  logical :: recover_wf_after_mixing = .true.

  real(kind=DP) :: amin_wf_precon = 0.01d0
  real(kind=DP) :: bmix_wf_precon =  1.0d0
  real(kind=DP) :: g0_wf_precon =    1.0d0
  real(kind=DP) :: edelta_start_wf_mixing = 1.0D3
!
! ====================== 13.0U3

  type charge_mixing_
     integer       :: mixing_way
     integer       :: iter_range
!!$     integer       :: variation_way = varLINEAR
!!$     integer       :: precon = YES
!!$     integer       :: hownew = ANEW
!!$     integer       :: cutoff = LARGE
!!$     integer       :: istr   = 1
!!$     integer       :: nbxmix = 0
!!$     real(kind=DP) :: rmxs = 0.10, rmxe = 0.50
     integer       :: variation_way
     integer       :: precon
     integer       :: hownew
     integer       :: cutoff
     integer       :: istr
     integer       :: nbxmix
     real(kind=DP) :: rmxs, rmxe
  end type charge_mixing_

!!!!  type(charge_mixing_),private,allocatable,dimension(:) :: charge_mixing
  type(charge_mixing_),allocatable,dimension(:) :: charge_mixing
!!!!  integer,private    :: n_Charge_Mixing_way = 1
  integer            :: n_Charge_Mixing_way = 1
  logical,public     :: tag_charge_mixing_is_found = .false.
  integer,private    :: ip_charge_mixing = 1

! -- DOS --
  integer,public ::        sw_dos = OFF
  integer,public ::        dos_method = Gauss_distrib_func
  integer,public ::        sw_dos_gaussdistrib = OFF
  integer,public ::        dos_subroutine = 5
  real(kind=DP),public ::  deltaE_dos = 1.d-4
  real(kind=DP),public ::  variance_dos_GaussD = 1.d-6
  real(kind=DP),public ::  dos_smearing_width  = 1.d-3
  integer,public ::        nwd_dos_window_width = 10
  integer,private ::       dos_energy_unit = HR_ENERGY_UNIT
  character(len("dos")),private,parameter ::                  tag_dos        = "dos"
  character(len("sw_dos")),private,parameter ::               tag_sw_dos     = "sw_dos"
  character(len("dos_method")),private,parameter ::           tag_dos_method = "dos_method"
  character(len("dos_subroutine")),private,parameter ::       tag_dos_subroutine = "dos_subroutine"
  character(len("gaussdistrib")),private,parameter ::         tag_gaussdistrib = "gaussdistrib"
  character(len("gaussiandistrib")),private,parameter ::      tag_gaussiandistrib = "gaussiandistrib"
  character(len("gaussian")),private,parameter ::             tag_gaussian   = "gaussian"
  character(len("sw_dos_gaussdistrib")),private,parameter ::  tag_sw_dos_gaussdistrib = "sw_dos_gaussdistrib"
  character(len("deltaE_dos_GaussD")),private,parameter ::    tag_deltaE_dos_GaussD   = "deltaE_dos_GaussD"
  character(len("deltaE_dos")),private,parameter ::           tag_deltaE_dos  = "deltaE_dos"
  character(len("variance_dos_GaussD")),private,parameter ::  tag_variance_dos_GaussD = "variance_dos_GaussD"
  character(len("variance_GaussD")),private,parameter ::      tag_variance_GaussD = "variance_GaussD"
  character(len("variance")),private,parameter ::             tag_variance = "variance"
  character(len("dos_smearing_width")),private,parameter ::  &
       &                      tag_dos_smearing_width = "dos_smearing_width"
  character(len("nwd_dos_window_width")),private,parameter :: tag_nwd_dos_window_width = "nwd_dos_window_width"
  character(len("energy_unit")),private,parameter ::          tag_energy_unit         = "energy_unit"

! ====================== added by K. Tagami ============================= 11.0
  character(len("calc_magmom_contrib")), private,parameter :: &
       &         tag_calc_magmom_contrib = "calc_magmom_contrib"
!
  integer, public ::  calc_dos_magmom_contrib = OFF
! ======================================================================= 11.0

! -- PDOS --
  integer,public ::        sw_orb_popu = OFF
  integer,public ::        sw_pdos = OFF
  integer,public ::        pdos_method = Wavefunction
  integer,allocatable,public :: norbital(:)   ! dim(ntyp)
  integer,public ::        maxorb
  integer,allocatable,public :: l_orb(:,:) ! dim(maxorb,ntyp)
  integer,allocatable,public :: t_orb(:,:) ! dim(maxorb,ntyp)
  real(kind=DP),allocatable,public :: rc_orb(:,:) ! dim(maxorb,ntyp)
  real(kind=DP),allocatable,public :: k_orb(:,:) ! dim(maxorb,ntyp)


! --- Orbital population --
  character(len("orbital_population")),private,parameter ::   &
       &               tag_orbital_population = "orbital_population"
  character(len("eval_population_method")),private,parameter :: &
       &         tag_orb_popu_method  = "eval_population_method"

  character(len("sw_diagonalize_population")),private,parameter :: &
       &         tag_sw_diagonalize_population = "sw_diagonalize_population"
  character(len("sw_write_rotated_orbitals")),private,parameter :: &
       &         tag_sw_write_rotated_orbitals = "sw_write_rotated_orbitals"

  character(len("population_diag_mode")),private,parameter :: &
       &         tag_population_diag_mode = "population_diag_mode"
  character(len("charge_density_matrix")),private,parameter :: &
       &         tag_charge_density_matrix = "charge_density_matrix"
  character(len("spin_density_matrix")),private,parameter :: &
       &         tag_spin_density_matrix = "spin_density_matrix"

  character(len("ls")),private,parameter :: &
       &         tag_ls = "ls"
  character(len("local_point_group")),private,parameter :: &
       &         tag_local_point_group = "local_point_group"
  character(len("local_double_point_group")),private,parameter :: &
       &         tag_local_double_point_group = "local_double_point_group"

  character(len("sw_calc_score_sigma_bond")),private,parameter :: &
       &         tag_sw_calc_score_sigma_bond = "sw_calc_score_sigma_bond"
  character(len("sw_write_orb_dens_mat_file")),private,parameter :: &
       &         tag_sw_write_orb_dens_mat_file = "sw_write_orb_dens_mat_file"
  character(len("sw_read_orb_rot_mat_file")),private,parameter :: &
       &         tag_sw_read_orb_rot_mat_file = "sw_read_orb_rot_mat_file"
  character(len("ipriorb_rot")),private,parameter :: &
       &         tag_ipriorb_rot = "ipriorb_rot"

  integer :: orb_popu_method = 1
  integer :: population_diag_mode = DIAG_CHARGE_DENSITY_MATRIX
  integer :: sw_diagonalize_population = OFF
  integer :: sw_write_rotated_orbitals = OFF

  integer :: sw_calc_score_sigma_bond = ON
  integer :: sw_write_orb_dens_mat_file = OFF
  integer :: sw_read_orb_rot_mat_file = OFF
  integer :: ipriorb_rot = -1

! -- PROCAR --
  character(len("procar")),private,parameter ::   &
       &               tag_procar = "procar"
  character(len("sw_write_procar_file")),private,parameter :: &
       &         tag_sw_write_procar_file = "sw_write_procar_file"
  character(len("sw_procar_full_bz")),private,parameter :: &
       &         tag_sw_procar_full_bz = "sw_procar_full_bz"
  character(len("procar_save_memory_mode")),private,parameter :: &
       &         tag_procar_save_memory_mode = "procar_save_memory_mode"
  character(len("split_procar_file")),private,parameter :: &
       &        tag_split_procar_file  = "split_procar_file"
  character(len("num_procar_files_once")),private,parameter :: &
       &        tag_num_procar_files_once  = "num_procar_files_once"
  character(len("procar_sort_kpt")),private,parameter :: &
       &        tag_procar_sort_kpt = "procar_sort_kpt"

  integer :: sw_write_procar_file = OFF
  integer :: sw_procar_full_bz = OFF
  integer :: procar_save_memory_mode = 1
  integer :: split_procar_file = OFF
  integer :: num_procar_files_once = 8
  integer :: procar_sort_kpt = OFF

! -- Parity --
  character(len("parity")),private,parameter ::   &
       &               tag_parity = "parity"
  character(len("sw_write_parity_file")),private,parameter :: &
       &         tag_sw_write_parity_file = "sw_write_parity_file"
  character(len("eval_parity_on_fftmesh")),private,parameter :: &
       &         tag_eval_parity_on_fftmesh = "eval_parity_on_fftmesh"
  integer :: sw_write_parity_file = OFF
  integer :: eval_parity_on_fftmesh = ON

! -- CHARGE --
  integer,public ::        sw_charge_rspace = OFF
  integer,public ::        sw_partial_charge = OFF
  integer,public ::        charge_filetype = DENSITY_ONLY
  integer,public ::        partial_charge_filetype =  INTEGRATED
  real(kind=DP),public ::  partial_charge_Emin =  -99999.99999 !  -0.4409852
  real(kind=DP),public ::  partial_charge_Emax =   99999.99999 !   0.4409852
  real(kind=DP),public ::  partial_charge_deltaE = 0.036748    !  1.0 eV
  character(len=LEN_TITLE) ::    charge_title = ''
  character(len=LEN_TITLE) ::    charge_title_tmp = ''

  character(len("charge")),private,parameter ::           tag_charge              = "charge"
  character(len("sw_charge_rspace")),private,parameter :: tag_sw_charge_rspace    = "sw_charge_rspace"
  character(len("filetype")),private,parameter   ::       tag_filetype            = "filetype"
  character(len("cube")),private,parameter ::             tag_cube                = "cube"
  character(len("vtk")),private,parameter ::              tag_vtk                 = "vtk"
  character(len("binary")),private,parameter ::           tag_binary              = "binary"
  character(len("density_only")),private,parameter ::     tag_density_only        = "density_only"
  character(len("title")),private,parameter ::            tag_title               = "title"
  character(len("sw_partial_charge")),private,parameter ::tag_sw_partial_charge   = "sw_partial_charge"
  character(len("partial_charge")),private,parameter ::   tag_partial_charge      = "partial_charge"
  character(len("Erange_min")),private,parameter ::       tag_Erange_min          = "Erange_min"
  character(len("Erange_max")),private,parameter ::       tag_Erange_max          = "Erange_max"
  character(len("Erange_delta")),private,parameter ::     tag_Erange_delta        = "Erange_delta"
  character(len("deltaE")),private,parameter ::           tag_deltaE              = "deltaE"
  character(len("deltaE_partial_charge")),private,parameter :: tag_deltaE_partial_charge = "deltaE_partial_charge"
  character(len("partial_charge_filetype")),private,parameter :: tag_partial_charge_filetype = "partial_charge_filetype"
  character(len("outputfiletype")),private,parameter ::   tag_outputfiletype      = "outputfiletype"
  character(len("individual")),private,parameter ::       tag_individual          = "individual"
  character(len("separate")),private,parameter ::         tag_separate            = "separate"
  character(len("integrated")),private,parameter ::       tag_integrated          = "integrated"

  character(len("sw_add_corecharge_rspace")),private,parameter :: &
       &       tag_sw_add_corecharge_rspace    = "sw_add_corecharge_rspace"
  character(len("sw_read_corecharge_extra_file")),private,parameter :: &
       &       tag_sw_read_corecharge_extra   = "sw_read_corecharge_extra_file"
  character(len("eval_corecharge_on_Gspace")),private,parameter :: &
       &       tag_eval_corecharge_on_Gspace = "eval_corecharge_on_Gspace"

  integer,public :: sw_add_corecharge_rspace = OFF
  integer,public :: sw_read_corecharge_extra_file = OFF
  integer,public :: eval_corecharge_on_Gspace = OFF

! ====== KT_add === 2014/06/07
  character(len("sw_spin_magmom_rspace")),private,parameter ::  &
       &                 tag_sw_spin_magmom_rspace = "sw_spin_magmom_rspace"
  integer ::  sw_spin_magmom_rspace = OFF
! ================= 2014/06/07

! ==== parameters for filtering the cube file
  character(len("sw_subset_only")),private,parameter :: tag_sw_subset_only = "sw_subset_only"
  character(len("min_x")),private,parameter :: tag_min_x = "min_x"
  character(len("min_y")),private,parameter :: tag_min_y = "min_y"
  character(len("min_z")),private,parameter :: tag_min_z = "min_z"
  character(len("max_x")),private,parameter :: tag_max_x = "max_x"
  character(len("max_y")),private,parameter :: tag_max_y = "max_y"
  character(len("max_z")),private,parameter :: tag_max_z = "max_z"
  integer :: sw_subset_only = OFF
  real(kind=DP),dimension(3) :: minxyz,maxxyz

! ==============================================================

! --- Wavefunction ---
  integer,public ::        sw_wf_rspace = OFF
  integer,public ::        wf_filetype = DENSITY_ONLY
  character(len=LEN_TITLE) ::    wf_title = ''
  real(kind=DP) :: eigmin_wf=-1.d+10, eigmax_wf=1.d+10
  character(len("wf")),private,parameter ::           tag_wf              = "wf"
  character(len("sw_wf_rspace")),private,parameter :: tag_sw_wf_rspace    = "sw_wf_rspace"
  character(len("eigenvalue")),private,parameter :: tag_eigenvalue = "eigenvalue"
  character(len("eigmin")),private,parameter :: tag_eigmin    = "eigmin"
  character(len("eigmax")),private,parameter :: tag_eigmax    = "eigmax"

! --- Radial wavefunction ---
  integer,public ::        sw_rwf2 = OFF
  integer,public ::        ib_rwf2 = 1
  integer,public ::        ik_rwf2 = 1
  integer,public ::        nr = 100
  real(kind=DP),public ::  rmax = 1.d0
  real(kind=DP),public ::  center(3) = (/0.d0,0.d0,0.d0/)
  character(len("rwf2")),private,parameter ::    tag_rwf2       = "rwf2"
  character(len("sw_rwf2")),private,parameter :: tag_sw_rwf2    = "sw_rwf2"
  character(len("ib")),private,parameter :: tag_ib    = "ib"
  character(len("ik")),private,parameter :: tag_ik    = "ik"
  character(len("nr")),private,parameter :: tag_nr    = "nr"
  character(len("rmax")),private,parameter :: tag_rmax = "rmax"
  character(len("center")),private,parameter :: tag_center = "center"
  character(len("rx")),private,parameter :: tag_rx = "rx"
  character(len("ry")),private,parameter :: tag_ry = "ry"
  character(len("rz")),private,parameter :: tag_rz = "rz"

! ---- WaveFunction Orb-projection
  integer :: sw_calc_wf_orb_projection = OFF
  integer :: wf_orb_proj_print_format = 0
  integer :: use_rotated_compri = OFF
  character(len("wf_orb_projection")),private,parameter ::  &
       &            tag_wf_orb_projection       = "wf_orb_projection"
  character(len("sw_calc_wf_orb_projection")),private,parameter :: &
       &       tag_sw_calc_wf_orb_projection   = "sw_calc_wf_orb_projection"
  character(len("wf_orb_proj_print_format")),private,parameter :: &
       &       tag_wf_orb_proj_print_format = "wf_orb_proj_print_format"
  character(len("use_rotated_compri")),private,parameter :: &
       &       tag_use_rotated_compri = "use_rotated_compri"

! ---- WaveFunction Squared
  integer :: sw_wf_squared_rspace = OFF
  integer,public ::  wf_squared_filetype = CUBE
  integer,public ::  ik_wf_squared = 1
  integer,public ::  ib1_wf_squared = 1
  integer,public ::  ib2_wf_squared = 1

  character(len("wf_squared")),private,parameter ::  &
       &            tag_wf_squared       = "wf_squared"
  character(len("sw_wf_squared_rspace")),private,parameter :: &
       &            tag_sw_wf_squared_rspace    = "sw_wf_squared_rspace"
  character(len("ib1")),private,parameter :: tag_ib1    = "ib1"
  character(len("ib2")),private,parameter :: tag_ib2    = "ib2"

  integer :: sw_wf_integ_moment = OFF
  character(len("sw_wf_integ_moment")),private,parameter :: &
       &            tag_sw_wf_integ_moment   = "sw_wf_integ_moment"

! --- ELF ---
  integer,public :: sw_elf = OFF
  integer,public :: elf_filetype = DENSITY_ONLY
  character(len=LEN_TITLE) ::    elf_title = ''
  character(len("elf")),private,parameter ::           tag_elf              = "elf"
  character(len("sw_elf")),private,parameter ::        tag_sw_elf           = "sw_elf"

  ! --- Berry phase ---
  integer :: sw_berry_phase = 0
  integer :: sw_displace_atom = 0
  integer :: polar_prop = NOPOLAR
  integer :: sw_bp_property = 0
  character(len("polarization")),private,parameter :: tag_polarization = "polarization"
  character(len("sw_bp_property")),private,parameter :: tag_sw_bp_property = "sw_bp_property"
  character(len("property")),private,parameter :: tag_property = "property"
  character(len("effective_charge")),private,parameter :: tag_effective_charge= "effective_charge"
  character(len("piezoelectric_const")),private,parameter :: tag_piezoelectric_const = "piezoelectric_const"

! == KT_add === 2014/06/30
  integer :: sw_check_sumrule_born_charge = YES
  character(len("sw_check_sumrule_born_charge")), private, parameter &
       &  :: tag_sw_chk_sumrule_born_charge = "sw_check_sumrule_born_charge"
! ============= 2014/06/30

  ! --- approximate DFT+U : Hubbard model ---
  integer :: sw_hubbard = OFF
  integer :: sw_constraint = OFF
  integer :: initial_occmat = OFF
  real(kind=DP) :: initial_occmat_factor=1.d0
  integer :: const_site = 0
  real(kind=DP) :: const_alpha = 0.d0
  real(kind=DP) :: critical_ehub = -1.d-6
  real(kind=DP) :: delta_ehub = -1.d-6 ! default value is negative.
  real(kind=DP) :: alpha_hubbard = 1.d0
! === Include PAW by tkato =====================================================
  integer :: sw_force_simple_mixing_hub = OFF
! ==============================================================================
  integer :: occ_matrix_fix_period = -1

  character(len("hubbard")),private,parameter :: tag_hubbard = "hubbard"
  character(len("initial_occmat")),private,parameter ::   tag_initial_occmat = "initial_occmat"
  character(len("initial_occmat_factor")),private,parameter ::   tag_initial_occmat_factor = "initial_occmat_factor"
  character(len("sw_hubbard")),private,parameter ::       tag_sw_hubbard = "sw_hubbard"
  character(len("delta_hubbard_energy")),private,parameter :: tag_delta_hubbard_energy = "delta_hubbard_energy"
  character(len("critical_hubbard_energy")),private,parameter :: tag_critical_hubbard_energy = "critical_hubbard_energy"
  character(len("component")),private,parameter :: tag_component = "component"
  character(len("Ueff")),private,parameter :: tag_Ueff = "Ueff"
  character(len("norbital")),private,parameter :: tag_norbital = "norbital"
  character(len("constraint")),private,parameter :: tag_constraint = "constraint"
  character(len("sw_constraint")),private,parameter :: tag_sw_constraint = "sw_constraint"
  character(len("site")),private,parameter :: tag_site = "site"
  character(len("alpha")),private,parameter :: tag_alpha = "alpha"

  integer :: dftu_type = FLL
  character(len("dftu_type")),private,parameter :: tag_dftu_type = "dftu_type"
  character(len("FLL")),private,parameter :: tag_fll = "fll"
  character(len("AMF")),private,parameter :: tag_amf = "amf"

! U-ramping
  character(len("initialUeff")),private,parameter :: tag_initialUeff = "initialUeff"
  character(len("finalUeff")),private, parameter :: tag_finalUeff = "finalUeff"
  character(len("nUeff")), private, parameter :: tag_nUeff = "nUeff"

  integer :: nUeff = 0

  ! === KT === 2019/12/31
  character(len("crystal_field_approx")),private,parameter :: &
       &         tag_crystal_field_approx = "crystal_field_approx"
  character(len("initial_spin_sum")),private,parameter :: &
       &         tag_initial_spin_sum = "initial_spin_sum"
  character(len("initial_spin_diff")),private,parameter :: &
       &         tag_initial_spin_diff = "initial_spin_diff"
  real(kind=DP) :: initial_spin_sum  = 0.0d0
  real(kind=DP) :: initial_spin_diff = 0.0d0
! === KT === 2019/12/31

! ================================== added by K. Tagami =================== 11.0
  integer :: occmat_diag_mode = 1
!
  character(len("occmat_diag_mode")),private,parameter :: &
       &        tag_occmat_diag_mode = "occmat_diag_mode"
!
  integer :: occmat_file_format = 2
  character(len("occmat_file_format")),private,parameter :: &
       &        tag_occmat_file_format = "occmat_file_format"
! ========================================================================= 11.0

  character(len("occ_matrix_fix_period")),private,parameter :: &
   & tag_occ_matrix_fix_period = "occ_matrix_fix_period"
  ! Projector set for PDOS, Hubbard, and Pseudo-SIC
  integer :: projector_type = ATOMIC_ORBITAL

!--- attribution of initial values is not allowed in the VPP-compiler
!    Initialization values have been removed by T. Yamasaki, 17th Feb. 2006
!!$  type t_projector
!!$     real(kind=DP) :: radius = 0.d0
!!$     real(kind=DP) :: fwhm = 0.1d0
!!$     integer :: l = 0
!!$     integer :: t = 1
!!$     logical, pointer :: component(:) ! dim(2l+1)
!!$     real(kind=DP) :: Ueff = 0.d0
!!$     logical :: frotate = .false.
!!$     real(kind=DP) :: phi=0.d0,theta=0.d0,psi=0.d0
!!$     real(kind=DP), pointer :: crotylm(:,:) ! dim(2l+1,2l+1)
!!$     integer :: ityp = 1
!!$     logical :: strong_correlated = .false.
!!$     integer :: norbital = 0
!!$     integer :: group = 0
!!$     integer :: ielem = 0
!!$  end type
  type t_projector
     real(kind=DP) :: radius
     real(kind=DP) :: fwhm
     integer :: l
     integer :: t
     logical, pointer :: component(:) ! dim(2l+1)
     real(kind=DP) :: Ueff
     logical :: frotate
     real(kind=DP) :: phi,theta,psi
     real(kind=DP), pointer :: crotylm(:,:) ! dim(2l+1,2l+1)
     integer :: ityp
     logical :: strong_correlated
     integer :: norbital
     integer :: group
     integer :: ielem

! ========================= added by K. Tagami ================= 11.0
     real(kind=DP) :: LScoupling0
     real(kind=DP) :: LScoupling_scaling_factor
! ============================================================== 11.0

     logical :: radius_was_defined

     real(kind=DP) :: initialUeff
     real(kind=DP) :: finalUeff
     real(kind=DP) :: deltaUeff
     real(kind=DP) :: initial_spin_sum, initial_spin_diff
  end type t_projector

  type(t_projector), allocatable, target :: proj_attribute(:) ! dim(num_projectors)
  integer :: num_projectors = 0
  integer, allocatable :: proj_group(:,:) ! dim(max_projs,num_proj_group)

! === KT_mod === 2014/06/05
!  integer :: max_projs = 4 ! = s,p,d,f

  integer :: max_projs
  integer :: max_projs_org = 4    ! = s,p,d,f
! ==============  2014/06/05

! === KT_add === 2014/06/05 & 2016/09/24
  integer, allocatable :: num_projs_in_each_group(:)
!
  integer :: sw_allow_maxprojs_gt_4 = OFF
  character(len("sw_allow_maxprojs_gt_4")),private,parameter :: &
       &       tag_sw_allow_maxprojs_gt_4 = "sw_allow_maxprojs_gt_4"

  character(len("howto_set_proj_radius")),private,parameter :: &
       &         tag_howto_set_proj_radius = "howto_set_proj_radius"
  integer :: howto_set_proj_radius = 0
                            ! 0: manual, 1: from-PP (paw pot only), 2: rad-cov
! ============== 2014/06/05 & 2016/09/24

  integer :: num_proj_group = 0
  integer, allocatable :: num_proj_elems(:) ! dim(num_proj_group)
  character(len("projector_list")),private,parameter :: tag_projector_list = "projector_list"
  character(len("projector_type")),private,parameter :: tag_projector_type = "projector_type"
  character(len("spherical_harmonics")),private,parameter :: tag_spherical_harmonics = "spherical_harmonics"
  character(len("atomic_orbital")),private,parameter :: tag_atomic_orbital = "atomic_orbital"
  character(len("wannier_function")),private,parameter :: tag_wannier_function = "wannier_function"

! ==============================-- added by K. Tagami =================== 5.0
  character(len("sw_orbital_cut_smooth")), private,parameter :: &
        &      tag_sw_orbital_cut_smooth = "sw_orbital_cut_smooth"
  character(len("occmat_type")), private,parameter :: &
        &      tag_occmat_type = "occmat_type"
  character(len("type1")), private,parameter :: &
        &      tag_occmat_type1 = "type1"
  character(len("type2")), private,parameter :: &
        &      tag_occmat_type2 = "type2"
!
  integer :: sw_orbital_cut_smooth = ON
!!  integer :: occmat_type = OCCMAT_Type2
  integer :: occmat_type = OCCMAT_Type1

  character(len("Ueff_imposition")), private, parameter :: &
        &      tag_Ueff_imposition = "Ueff_imposition"
  character(len("method")), private, parameter :: &
        &      tag_method_Ueff_imposition = "method"
  character(len("from_first")), private, parameter :: &
        &      tag_Ueff_from_first = "from_first"
  character(len("gradually")), private, parameter :: &
        &      tag_Ueff_gradually = "gradually"
  character(len("transition_period")), private, parameter :: &
        &      tag_Ueff_transition_period = "transition_period"
  character(len("edelta_for_Ueff_starting")), private, parameter :: &
        &      tag_edelta_for_Ueff_starting = "edelta_for_Ueff_starting"
  character(len("iteration_for_Ueff_starting")), private, parameter :: &
        &      tag_iter_for_Ueff_starting = "iteration_for_Ueff_starting"

  integer :: method_Ueff_imposition = Ueff_From_First
  integer :: Ueff_transition_period =  40
  integer :: iteration_for_Ueff_starting = 20
  real(kind=DP) :: edelta_for_Ueff_starting = 1.0D-5
!
  integer :: sw_eval_Ueff_using_iter = OFF
  integer :: sw_eval_Ueff_using_edelta = OFF
! ======================================================================= 5.0

  character(len("projectors")),private,parameter :: tag_projectors = "projectors"
  character(len("radius")),private,parameter :: tag_radius = "radius"
  character(len("fwhm")),private,parameter :: tag_fwhm = "fwhm"
  character(len("l")),private,parameter :: tag_l = "l"
  character(len("t")),private,parameter :: tag_t = "t"
  character(len("phi")),private,parameter :: tag_phi = "phi"
  character(len("theta")),private,parameter :: tag_theta = "theta"
  character(len("psi")),private,parameter :: tag_psi = "psi"
  character(len("group")),private,parameter :: tag_group = "group"

  ! --- fine STM_images ---
  integer,public       :: sw_fine_STM_simulation = OFF
!!  integer,public     :: sw_deficit_charge = OFF
  integer,public       :: sw_deficit_charge = ON
  integer,public       :: z_axis = 3
  real(kind=DP),public :: connect_from = 3.d0
  character(len("STM")),private,parameter ::                    tag_STM = "STM"
  character(len("sw_STM")),private,parameter ::                 tag_sw_stm = "sw_STM"
  character(len("sw_fine_STM")),private,parameter ::            tag_sw_fine_STM = "sw_fine_STM"
  character(len("sw_fine_STM_simulation")),private,parameter :: tag_sw_fine_STM_simulation = "sw_fine_STM_simulation"
  character(len("sw_STM_images")),private,parameter ::          tag_sw_STM_images = "sw_STM_images"
  character(len("sw_deficit_charge")),private,parameter ::      tag_sw_deficit_charge = "sw_deficit_charge"
  character(len("z_axis")), private, parameter ::               tag_z_axis = "z_axis"
  character(len("connect_from")), private, parameter ::         tag_connect_from = "connect_from"

! -- workfunction --
  integer,public     :: sw_add_xc_to_vloc = ON
  integer,public     :: sw_xc_only = OFF
  character(len("sw_add_xc_to_vloc")),private,parameter ::      tag_sw_add_xc_to_vloc = "sw_add_xc_to_vloc"
  character(len("sw_xc_only")),private,parameter        ::      tag_sw_xc_only = "sw_xc_only"
  character(len("workfunc")),private,parameter          ::      tag_workfunc = "workfunc"
  character(len("sw_workfunc")),private,parameter       ::      tag_sw_workfunc = "sw_workfunc"

  ! --- SC_DFT method---
  character(len("sc_dft")),private,parameter    ::   tag_sc_dft = "sc_dft"
  character(len("scdft")),private,parameter    ::   tag_scdft = "scdft"
  character(len("delta_epsilon")), private, parameter :: tag_delta_epsilon = "delta_epsilon"
  real(kind=DP), public :: delta_epsilon = 0.01
  real(kind=DP), public :: epsilon0, epsilon0_previous

  ! --- Hybrid functional method---
  integer, public :: sw_hybrid_functional = OFF
  integer, public :: sw_exchange_only = OFF
  integer, public :: sw_screened_exchange = OFF
  integer, public :: sw_singular_correction = ON
  integer, public :: sw_memory_reduction_exx = OFF
  integer :: reduction_factor_exx(3) = (/1,1,1/)
  real(kind=DP), public :: alpha_exx = 0.25d0
  logical, public :: alpha_exx_is_set_in_SCDFTLoop = .false.
  real(kind=DP), public :: omega_exx = 0.106d0 ! Bohr^-1
  real(kind=DP), public :: omega_exx_pbe = 0.106d0 ! Bohr^-1
  character(len("hybrid_functional")),private,parameter    ::   tag_hybrid_functional = "hybrid_functional"
  character(len("sw_hybrid_functional")),private,parameter ::   tag_sw_hybrid_functional = "sw_hybrid_functional"
  character(len("sw_exchange_only")),private,parameter      ::   tag_sw_exchange_only = "sw_exchange_only"
  character(len("sw_screened_exchange")),private,parameter ::   tag_sw_screened_exchange = "sw_screened_exchange"
  character(len("sw_memory_reduction")),private,parameter ::   tag_sw_memory_reduction = "sw_memory_reduction"
  character(len("sw_singular_correction")),private,parameter  ::   tag_sw_singular_correction = "sw_singular_correction"
  character(len("reduction_factor")),private,parameter      ::   tag_reduction_factor = "reduction_factor"
  character(len("f1")),private,parameter                ::   tag_f1 = "f1"
  character(len("f2")),private,parameter                ::   tag_f2 = "f2"
  character(len("f3")),private,parameter                ::   tag_f3 = "f3"
  !character(len("alpha")),private,parameter                ::   tag_alpha = "alpha"
  character(len("omega")),private,parameter                ::   tag_omega = "omega"
  character(len("omega")),private,parameter                ::   tag_omega_hf  = "omega_hf"
  character(len("omega_pbe")),private,parameter            ::   tag_omega_pbe = "omega_pbe"
  ! HF    : alpha=1.00 & sw_exchange_only=ON
  ! PBE0  : alpha=0.25 with PBE (defalut hybrid functional)
  ! HSE06 : alpha=0.25 & omega=0.106 with PBE

  character(len("functional_type")),private,parameter :: tag_functional_type = "functional_type"
  character(len("pbe0")),private,parameter  :: tag_pbe0  = "pbe0"
  character(len("hse06")),private,parameter :: tag_hse06 = "hse06"
  character(len("hf")),private,parameter    :: tag_hf    = "hf"
  character(len("gaupbe")),private,parameter :: tag_gaupbe = "gaupbe"
  character(len("hse06-hjs")),private,parameter :: tag_hse06_hjs = "hse06-hjs"
  character(len("hsesol-hjs")),private,parameter :: tag_hsesol_hjs = "hsesol-hjs"

  integer, parameter :: len_hybrid_pot_type = 10
  character(len=len_hybrid_pot_type) :: hybrid_functional_type = 'nogiven'
! ===================== KT_Test  ================================== 12.5Exp
  character(len("nmax_G_hyb")), private, parameter :: &
       &                         tag_nmax_G_hyb = "nmax_G_hyb"
  character(len("charge_mesh")), private ,parameter :: &
       &                         tag_charge_mesh = "charge_mesh"
  character(len("exact")),private,parameter :: &
       &                         tag_exact  = "exact"
  character(len("fine")),private,parameter :: &
       &                         tag_fine   = "fine"
  character(len("moderate")),private,parameter :: &
       &                         tag_moderate = "moderate"
  character(len("coarse")),private,parameter :: &
       &                         tag_coarse    = "coarse"
  integer :: nmax_G_hyb = 0
!
  character(len("truncate_vxw_updating")), private, parameter :: &
       &                         tag_truncate_vxw_updating = "truncate_vxw_updating"
  logical :: truncate_vxw_updating = .false.
  integer :: sw_update_vxw = ON
!
  character(len("delta_total_energy_hyb1")), private, parameter :: &
       &                      tag_delta_total_energy_hyb1 = "delta_total_energy_hyb1"
  character(len("delta_total_energy_hyb2")), private, parameter :: &
       &                      tag_delta_total_energy_hyb2 = "delta_total_energy_hyb2"
!
  real(kind=DP) :: edelta_for_hyb_chgfix = 1.0d-5
  real(kind=DP) :: edelta_for_hyb_convgd = 1.0d-7
! ================================================================= 12.5Exp

! ======= KT_add  ===================================== 13.0F
  character(len("gmax_exx_ratio")), private, parameter :: &
       &                         tag_gmax_exx_ratio = "gmax_exx_ratio"
  character(len("gmaxp_exx_ratio")), private, parameter :: &
       &                         tag_gmaxp_exx_ratio = "gmaxp_exx_ratio"
!
  character(len("cutoff_wf_for_exx")), private, parameter :: &
       &                         tag_cutoff_wf_for_exx = "cutoff_wf_for_exx"
  character(len("cutoff_cd_for_exx")), private, parameter :: &
       &                         tag_cutoff_cd_for_exx = "cutoff_cd_for_exx"
!
  character(len("edelta_change_to_hybrid")), private, parameter :: &
       &                      tag_edelta_change_to_hybrid = "edelta_change_to_hybrid"
!
  real(kind=DP) :: gmax_exx = 0.0
  real(kind=DP) :: gmaxp_exx = 0.0
  real(kind=DP) :: edelta_change_to_hybrid = 1.0d-5
!
  logical :: use_hybrid_functional = .false.
  logical :: hybrid_calc_is_active = .false.
!
  logical :: use_fft_exx = .false.
! ==================================================== 13.0F

  integer :: potential_update=1
  character(len("potential_update")), private, parameter :: tag_potential_update="potential_update"
  character(len("always")), private, parameter :: tag_always = "always"
  character(len("minimal")), private, parameter :: tag_minimal = "minimal"

  character(len("sw_change_axis")), private, parameter ::  &
   &         tag_sw_change_axis="sw_change_axis"
!  integer :: sw_change_axis = OFF
  integer :: sw_change_axis = ON

  character(len("sw_output_hybrid_info")), private, parameter :: tag_sw_output_hybrid_info="sw_output_hybrid_info"
  integer :: sw_output_hybrid_info = OFF

  integer :: accuracy_of_exx = 3

! ======= KT_add === 13.0Y
! Partial Core correction ( paw, hybrid,... )
!
  character(len("partial_core_correction")), private, parameter :: &
       &                  tag_partial_core_correction = "partial_core_correction"
  character(len("pcc")), private, parameter :: &
       &                  tag_pcc = "pcc"
  character(len("sw_eval_epc_on_fftmesh")), private, parameter :: &
       &                      tag_sw_eval_epc_on_fftmesh = "sw_eval_epc_on_fftmesh"
  character(len("sw_remove_pcc_from_pawpot")), private, parameter :: &
       &                      tag_sw_remove_pcc_from_pawpot = "sw_remove_pcc_from_pawpot"
!
  integer :: sw_eval_epc_on_fftmesh = off
  integer :: sw_remove_pcc_from_pawpot = off
! ================== 13.0Y

!  integer :: sw_rspace_hyb = OFF
!  integer :: sw_rspace_hyb_dgm = OFF
  integer :: sw_rspace_hyb = ON
  integer :: sw_rspace_hyb_dgm = ON
  character(len("sw_rspace_dgemm")), private, parameter :: tag_sw_rspace_dgm = "sw_rspace_dgemm"
  character(len("sw_eval_vexx")), private, parameter :: tag_sw_eval_vexx = "sw_eval_vexx"
  character(len("sw_retard_eigval_evaluation")), private, parameter :: &
   & tag_sw_retard_eigval_evaluation = "sw_retard_eigval_evaluation"
  character(len("sw_precalculate")), private, parameter :: &
   & tag_sw_precalculate = "sw_precalculate"
  integer :: sw_eval_vexx = OFF
  integer :: sw_retard_eigval_evaluation = ON
  integer :: sw_precalculate = OFF

  ! -- RTTDDFT --
  integer,public :: sw_rttddft = OFF
  integer,public :: time_step_max = 100
  real(kind=DP),public :: time_step_delta = 0.1
  integer,public :: propagator_method = 1
  integer,public :: propagator_order = 4
  integer,public :: ext_ie_elec = 0, ext_ie_hole = 0
  real(kind=DP),public :: ext_pulse_epsilon = 0.0d0
  real(kind=DP),public :: ext_pulse_kx = 0.0d0, ext_pulse_ky = 0.0d0, ext_pulse_kz = 0.0d0
  character(len("rttddft")),private,parameter :: tag_rttddft = "rttddft"
  character(len("sw_rttddft")),private,parameter :: tag_sw_rttddft = "sw_rttddft"
  character(len("sw_time_step_max")),private,parameter :: tag_time_step_max = "time_step_max"
  character(len("sw_time_step_delta")),private,parameter :: tag_time_step_delta = "time_step_delta"
  character(len("sw_propagator_method")),private,parameter :: tag_propagator_method = "propagator_method"
  character(len("sw_propagator_order")),private,parameter :: tag_propagator_order = "propagator_order"
  character(len("sw_ext_ie_elec")),private,parameter :: tag_ext_ie_elec = "ext_ie_elec"
  character(len("sw_ext_ie_hole")),private,parameter :: tag_ext_ie_hole = "ext_ie_hole"
  character(len("sw_ext_pulse_epsilon")),private,parameter :: tag_ext_pulse_epsilon = "ext_pulse_epsilon"
  character(len("sw_ext_pulse_kx")),private,parameter :: tag_ext_pulse_kx = "ext_pulse_kx"
  character(len("sw_ext_pulse_ky")),private,parameter :: tag_ext_pulse_ky = "ext_pulse_ky"
  character(len("sw_ext_pulse_kz")),private,parameter :: tag_ext_pulse_kz = "ext_pulse_kz"

#else
  ! --- Structure evolution ---
  character(len("method")),private,parameter ::       tag_method          = "method"
  character(len("renew")),private,parameter ::        tag_renew           = "renew"
  character(len("anew")),private,parameter ::         tag_anew            = "anew"
  character(len("sd")),private,parameter ::           tag_sd                = "sd"
  character(len("cg")),private,parameter ::           tag_cg                = "cg"
  character(len("cg2")),private,parameter ::           tag_cg2              = "cg2"
  character(len("oldcg")),private,parameter ::         tag_oldcg            = "oldcg"
  ! --- Berry phase ---
  integer :: sw_displace_atom = 0
!!$  ! --- PRINTOUT LEVEL ---
!!$  character(len("printoutlevel")),private,parameter ::  tag_printoutlevel     = "printoutlevel"
!!$  character(len("printlevel")),private,parameter ::     tag_printlevel        = "printlevel"

#endif

! ==================== Added by K. Tagami ======================= 10.1
  integer :: Num_q_Points = 0
  integer :: sw_LinearResponse = 0
! =============================================================== 10.1

! --
  logical,public     :: paramset = .false.

  real(kind=DP), parameter :: cdel_critical = 1.d-5

  integer, parameter :: len_str = 132
  character(len=len_str) :: str
  character(len=FMAXVALLEN),private :: rstr

  real(kind=DP),pointer,private,dimension(:) :: work
  real(kind=DP), save   :: wct_start = 0.d0
#ifdef _CHECK_ELAPSE_
!                 wct_start: remaining time at start(sec).
!                 wct_now  : remaining time now(sec).
  real(kind=DP), parameter :: Critical_Remaining_CPU_TIME = 600  !(sec.)
#else
!                 wct_start: Wall Clock Time Start.
!                 wct_now  : Wall Clock Time Now.
#endif

! ============================== added by K. Tagami ================ 11.0
! Spin-Orbit
! --
  integer :: iprispinorb = -1
  character(len("iprispinorb")),private,parameter :: tag_iprispinorb   = "iprispinorb"

  character(len("spinorbit")),private,parameter :: tag_spinorbit = "spinorbit"
  character(len("mode")),private,parameter ::  tag_spinorbit_mode = "mode"

  character(len("LScoupling")),private,parameter :: tag_LScoupling = "LScoupling"
  character(len("Splitting")),private,parameter :: tag_splitting = "Splitting"
  character(len("scaling_factor")),private,parameter :: tag_scaling_factor &
        &                                  = "scaling_factor"
!
  character(len("projector")), private,parameter :: &
       &                                  tag_spinorbit_by_projector = "projector"
  character(len("neglected")),private,parameter :: tag_spinorbit_neglected = "neglected"
  character(len("builtin")),  private,parameter :: tag_spinorbit_builtin = "builtin"
  character(len("pawpot")),  private,parameter   :: tag_spinorbit_by_pawpot = "pawpot"
  character(len("zeff")),  private,parameter   :: tag_spinorbit_by_zeff = "zeff"
  character(len("read_from_pp")),  private,parameter   :: &
       &                        tag_spinorbit_read_from_pp = "read_from_pp"
!
  character(len("mass_correction")),  private,parameter   :: &
       &                tag_spinorbit_mass_correction = "mass_correction"
!
  character(len("sw_use_rphi_Hsoc_rphi")),  private,parameter   :: &
       &             tag_sw_use_rphi_Hsoc_rphi = "sw_use_rphi_Hsoc_rphi"
  character(len("sw_use_ival_for_paw_ps_soc")),  private,parameter   :: &
       &             tag_sw_use_ival_for_paw_ps_soc = "sw_use_ival_for_paw_ps_soc"
!
  integer :: SpinOrbit_Mode = Neglected
  integer :: SpinOrbit_MassCorrection = 1
!
  integer :: sw_use_rphi_Hsoc_rphi = ON
  integer :: sw_use_ival_for_paw_ps_soc = OFF
! =================================================================== 11.0

! Open core
  character(len("core_electrons")),private,parameter :: &
       &         tag_core_electrons = "core_electrons"
  character(len("sw_opencore")),private,parameter :: &
       &         tag_sw_opencore = "sw_opencore"
  character(len("sw_xc_opencore_ae_only")),private,parameter :: &
       &         tag_sw_xc_opencore_ae_only = "sw_xc_opencore_ae_only"
  character(len("spin_orientation")),private,parameter :: &
       &         tag_spin_orientation = "spin_orientation"
  character(len("parallel")),private,parameter :: &
       &         tag_parallel = "parallel"
  character(len("anti_parallel")),private,parameter :: &
       &         tag_anti_parallel = "anti_parallel"
  character(len("sw_fix_core_spin_pol")),private,parameter :: &
       &         tag_sw_fix_core_spin_pol = "sw_fix_core_spin_pol"
  integer :: sw_opencore = OFF
  integer :: sw_xc_opencore_ae_only = OFF
  real(kind=DP) :: core_spin_pol_factor = 1.0d0       ! parallel
  integer :: sw_fix_core_spin_pol = OFF

! Charge State
  integer :: sw_calc_extfnv_correction = OFF

! Local potential       ! experimental,  2020/09/08
  character(len("local_potential")),private,parameter :: &
       &         tag_local_potential = "local_potential"
  character(len("long_range_pot_type")),private,parameter :: &
       &         tag_long_range_pot_type = "long_range_pot_type"
  character(len("sw_add_vlocal_gzerol")),private,parameter :: &
       &         tag_sw_add_vlocal_gzero = "sw_add_vlocal_gzero"
  integer :: long_range_pot_type = 0
  integer :: sw_add_vlocal_gzero = OFF

! ==== KT_add ==== 2014/08/26
! Orbital decomposition
!
  character(len("orbital_decomposition")),private,parameter :: &
       &         tag_orbital_decomposition = "orbital_decomposition"
  character(len("decomp_mode")),private,parameter :: &
       &         tag_orbital_decomp_mode = "decomp_mode"
  character(len("sw_calc_orbital_moment")),private,parameter :: &
       &          tag_sw_calc_orbital_moment = "sw_calc_orbital_moment"
  character(len("sw_use_contracted_psir")),private,parameter :: &
       &          tag_sw_use_contracted_psir = "sw_use_contracted_psir_val"
! --
  integer :: Orbital_Decomp_Mode = 0
  integer :: sw_calc_orbital_moment = OFF
  integer :: sw_use_contracted_psir = OFF
! ================ 2014/08/26

! ============================== added by K. Tagami ================ 11.0
! Magnetic Moment
! --
  integer :: iprimagmom = -1
  character(len("iprimagmom")),private,parameter :: tag_iprimagmom   = "iprimagmom"
! =================================================================== 11.0

  integer,save :: driver = DRIVER_GENERAL
  character(len("driver")), parameter :: tag_driver ="driver"

#ifdef ENABLE_ESM_PACK
! ====================== ESM ==========================
  character(len("esm")),private,parameter :: tag_esm = "esm"
  character(len("sw_esm")),private,parameter :: tag_sw_esm = "sw_esm"
  character(len("z1")),private,parameter :: tag_z1 = "z1"
  character(len("w")),private,parameter :: tag_w = "w"
  character(len("fix_ef")),private,parameter :: tag_fix_ef = "fix_ef"
  character(len("add_elec")),private,parameter :: tag_add_elec = "add_elec"
  character(len("z_wall")),private,parameter :: tag_z_wall = "z_wall"
  character(len("bar_height")),private,parameter :: tag_bar_height = "bar_height"
  character(len("bar_width")),private,parameter :: tag_bar_width = "bar_width"
  character(len("nosmooth")),private,parameter :: tag_nosmooth = "nosmooth"
  character(len("external_potential")),private,parameter :: tag_external_potential="external_potential"
  character(len("gps")),private,parameter :: tag_gps = "gps"
  character(len("gpe")),private,parameter :: tag_gpe = "gpe"
  character(len("bc")),private,parameter :: tag_bc="bc"
  character(len("bare")),private,parameter :: tag_bare = "bare"
  character(len("pe1")),private,parameter :: tag_pe1 = "pe1"
  character(len("pe2")),private,parameter :: tag_pe2 = "pe2"
  integer :: sw_esm = OFF
  real(kind=DP) :: esm_z1=-1000.d0
  real(kind=DP) :: esm_w=0.d0
  logical :: esm_z1_defined = .false.
  integer :: esm_izwall = 0
  integer :: esm_iexpot = 0
  real(kind=DP) :: esm_z_wall = 0.d0
  real(kind=DP) :: esm_bar_height = 0.d0
  real(kind=DP) :: esm_bar_width = 0.d0
  real(kind=DP) :: esm_e_field = 0.d0
  real(kind=DP) :: esm_fix_ef = 0.d0
  real(kind=DP) :: esm_add_elec=0.d0
  real(kind=DP) :: esm_qbac = 0.d0
  integer :: esm_gps = 1
  integer :: esm_gpe = 1
  integer :: esm_bc = BARE
! ====================== ESM ==========================
#endif

! ====================== FCP ==========================
  character(len("fcp")),private,parameter :: tag_fcp = "fcp"
  character(len("sw_fcp")),private,parameter :: tag_sw_fcp= "sw_fcp"
  character(len("mu")),private,parameter :: tag_fcp_mu = "mu"
  character(len("mass")),private,parameter :: tag_fcp_mass = "mass"
  character(len("temperatue")),private,parameter :: tag_fcp_temperature = "temperature"
  character(len("relax_step")),private,parameter :: tag_fcp_relax_step = "relax_step"
  character(len("relax_crit")),private,parameter :: tag_fcp_relax_crit = "relax_crit"
  character(len("tot_charge_first")),private,parameter :: tag_fcp_tot_charge_first = "tot_charge_first"
  character(len("tot_charge_last")),private,parameter :: tag_fcp_tot_charge_last = "tot_charge_last"
  character(len("qmass")),private,parameter :: tag_fcp_qmass = "qmass"
  integer       :: sw_fcp               = OFF
  real(kind=DP) :: fcp_mu               = 0.0d0
  real(kind=DP) :: fcp_mass             = 10000.0d0
  real(kind=DP) :: fcp_temperature      = 0.0d0
  real(kind=DP) :: fcp_relax_step       = 0.5d0
  real(kind=DP) :: fcp_relax_crit       = 0.001d0
  real(kind=DP) :: fcp_tot_charge_first = 0.0d0
  real(kind=DP) :: fcp_tot_charge_last  = 0.0d0
  real(kind=DP) :: fcp_qmass            = 4000.0d0
! ====================== FCP ==========================

! === variables for charge & wf prediction after the update of atomic coordinates
  character(len("predictor")),private,parameter :: tag_predictor="predictor"
  character(len("sw_charge_predictor")),private,parameter :: tag_sw_charge_predictor = "sw_charge_predictor"
  character(len("sw_wf_predictor")),private,parameter :: tag_sw_wf_predictor = "sw_wf_predictor"
  character(len("sw_extrapolate_charge")),private,parameter :: tag_sw_extrapolate_charge = "sw_extrapolate_charge"
  character(len("rms_threshold")),private,parameter :: tag_rms_threshold="rms_threshold"
  integer :: sw_charge_predictor = ON
  integer :: sw_wf_predictor = OFF
  integer :: sw_extrapolate_charge = ON
  real(kind=DP) :: rms_threshold=0.1d0

  logical,private :: in_initialization = .true.
#ifdef MEMORY_SAVE_ZAJ_OLD
  logical :: RMM2P_is_specified = .false.
#endif

! nonlocal potential in rspace
  character(len("sw_rspace")),private,parameter :: tag_sw_rspace = "sw_rspace"
  character(len("sw_rspace_v")),private,parameter :: tag_sw_rspace_v = "sw_rspace_v"
  character(len("nonlocal_potential")),private,parameter :: tag_nonlocal_potential = "nonlocal_potential"
  character(len("r0_factor")),private,parameter :: tag_r0_factor = "r0_factor"
  character(len("gamma_factor")),private,parameter :: tag_gamma_factor = "gamma_factor"
!  character(len("nq")),private,parameter :: tag_nq = "nq"
  character(len("dq")),private,parameter :: tag_dq = "dq"

  character(len("projector_optimization")),private,parameter :: tag_projector_optimization = "projector_optimization"
  character(len("qr_optimization")),private,parameter :: tag_qr_optimization = "qr_optimization"
  character(len("prefitting")),private,parameter :: tag_prefitting = "prefitting"
  character(len("mask_function")),private,parameter :: tag_mask_function = "mask_function"

  character(len("sw_save_memory")),private,parameter :: tag_sw_save_memory = "sw_save_memory"

  integer :: sw_rspace = OFF
  integer :: sw_rspace_v = OFF
  real(kind=DP) :: r0_factor = 1.9d0
  real(kind=DP) :: r0_factor_q = 1.0d0
  real(kind=DP) :: gamma_factor = 2.d0
  !!$integer :: nq = 1000
  real(kind=DP) :: dq = 0.005d0
  integer :: projector_optimization = MASK_FUNCTION
  integer :: qr_optimization = MASK_FUNCTION
  integer :: sw_save_memory = ON

!-- RSB
  character(len("rsb")),private,parameter :: tag_rsb = "rsb"
  character(len("sw_rsb")),private,parameter :: tag_sw_rsb = "sw_rsb"
  character(len("rsb_test")),private,parameter :: tag_rsb_test="rsb_test"
  character(len("sw_rsb_test")),private,parameter :: tag_sw_rsb_test="sw_rsb_test"
  character(len("bisect_by")), private,parameter :: tag_bisect_by="bisect_by"
  character(len("lmax")), private, parameter :: tag_lmax="lmax"
  character(len("eps_rsb")), private, parameter :: tag_eps_rsb="eps_rsb"
  character(len("sw_valence_electrons_only")), private, parameter :: &
  &  tag_sw_valence_electrons_only = "sw_valence_electrons_only"
  integer :: sw_rsb_test = OFF
  integer :: sw_rsb = OFF
  integer :: sw_valence_electrons_only = OFF
  integer :: bisect_by = 0
  integer :: lmax_rsb = 3
  real(kind=DP) :: eps_rsb = 1.d-2

! --- msb effect
  character(len("msb")), private, parameter :: &
       &                     tag_msb = "msb"
  character(len("sw_calc_contact_density")), private, parameter :: &
       &                     tag_sw_calc_contact_density = "sw_calc_contact_density"
  integer :: sw_calc_contact_density = off

! -- band unfolding
  character(len("band_unfolding")), private, parameter :: &
       &                     tag_band_unfolding = "band_unfolding"
  character(len("sw_band_unfolding")), private, parameter :: &
       &                     tag_sw_band_unfolding = "sw_band_unfolding"
  character(len("tolerance_gvec_matching")), private, parameter :: &
       &                     tag_tol_gvec_matching = "tolerance_gvec_matching"
  integer :: sw_band_unfolding = off
  logical :: band_unfolding_active = .false.
  real(kind=DP) :: tolerance_Gvec_matching = 1.0D-3

  character(len("prepare_masked_compri")), private, parameter :: &
       &           tag_prepare_masked_compri = "prepare_masked_compri"
  integer :: prepare_masked_compri = ON

! ================= KT_add === 13.0S
!-- CoreLevels
!
  character(len("corelevels")), private,parameter :: tag_corelevels="corelevels"
  character(len("sw_calc_core_energy")), private,parameter :: &
       &                        tag_sw_calc_core_energy = "sw_calc_core_energy"
  character(len("corehole")), private,parameter :: &
       &                        tag_corehole = "corehole"
  character(len("atom_id")), private,parameter :: &
       &                        tag_atom_id = "atom_id"
  character(len("orbital")),   parameter  :: tag_orbital  = "orbital"
!
  integer :: sw_calc_core_energy = off
  integer :: atom_with_corehole = 0
  integer :: qnum_n_corehole = 0, qnum_l_corehole = 0
!
! ------
  character(len("eval_core_level_splitting")),   parameter :: &
       &            tag_eval_core_level_splitting = "eval_core_level_splitting"
!
  character(len("pawpot")),  private,parameter :: &
       &    tag_core_level_splitting_pawpot = "pawpot"
  character(len("read_from_pp")),  private,parameter :: &
       &    tag_core_level_splitting_frompp = "read_from_pp"
!
!  integer :: eval_core_level_splitting = ByPawPot
  integer :: eval_core_level_splitting = 0
! ============================ 13.0S

! ===== KT_add ===== 13.0U3
!-- smearing update
!
  character(len("method_change_width")), private, parameter :: &
       &                          tag_method_change_width = "method_change_width"
  character(len("stepwise")),private,parameter ::   &
       &                          tag_stepwise = "stepwise"

  character(len("width_initial")),private,parameter ::   &
       &                          tag_width_initial = "width_initial"
!
  character(len("num_intermid_width")), private, parameter :: &
       &                          tag_num_intermid_width = "num_intermid_width"
  character(len("edelta_change_width_first")), private, parameter :: &
       &                 tag_edelta_change_width_first = "edelta_change_width_first"
  character(len("edelta_change_width_last")), private, parameter :: &
       &                 tag_edelta_change_width_last = "edelta_change_width_last"
!
  integer, parameter :: nmax_intermid_width = 50

  integer :: method_change_smearing_width = 0
  integer :: num_intermid_width = 0
!
  real(kind=DP) :: edelta_change_width_first = 1.0D-6     ! hartree
  real(kind=DP) :: edelta_change_width_last = 1.0D-6     ! hartree
  real(kind=DP) :: width_initial = 0.008
! ================== 13.0U3

! ============== KT_add ====== 13.0XX, 2014/09/19
! -- meta GGA
  character(len("megagga")),   private, parameter :: &
       &                      tag_metagga = "metagga"
  character(len("val_c_mbj09")),   private, parameter :: &
       &                      tag_val_c_mbj09 = "val_c_mbj09"
  character(len("val_c_tb09")),   private, parameter :: &
       &                      tag_val_c_tb09 = "val_c_tb09"
  logical :: use_metagga = .false.
  integer :: sw_calc_val_g_tb09_paw = OFF
  real(kind=DP) :: val_g_tb09_paw = 0.0d0

  integer :: sw_fix_val_c_tb09 = OFF
  real(kind=DP) :: val_c_tb09 = 1.0d0

  logical :: vtau_exists = .false.

#ifdef LIBXC
! --- libxc
  character(len("libxc")), private, parameter :: &
       &                 tag_libxc = "libxc"
  character(len("exch_id")), private, parameter :: &
       &                 tag_exch_id = "exch_id"
  character(len("corr_id")), private, parameter :: &
       &                 tag_corr_id = "corr_id"
  character(len("exch_name")), private, parameter :: &
       &                 tag_exch_name = "exch_name"
  character(len("corr_name")), private, parameter :: &
       &                 tag_corr_name = "corr_name"
  integer :: libxc_exch_id = 0, libxc_corr_id = 0
  integer :: xc_family_exch = 0, xc_family_corr = 0
  character*30 :: libxc_exch_name = "", libxc_corr_name = ""
  character*30 :: xc_name_exch = "", xc_name_corr = ""

  integer, parameter :: num_item_flag = 15
  integer :: xc_flag_exch(0:num_item_flag),  xc_flag_corr(0:num_item_flag)

  TYPE(xc_f03_func_t) :: xc_func_exch,  xc_func_corr
  TYPE(xc_f03_func_info_t) :: xc_info_exch, xc_info_corr
  TYPE(xc_f03_func_reference_t) :: xc_ref_exch, xc_ref_corr
#endif

! -- Kinetic Energy Density
!
  integer :: sw_calc_ekin_density = off
  logical :: use_symm_ekin_density  = .false.
  logical :: use_asymm_ekin_density = .false.
  logical :: ekin_density_is_active = .true.
!
  character(len("ekin_density_type")),   private, parameter :: &
       &                      tag_ekin_density_type = "ekin_density_type"
  integer :: ekin_density_type = 0
  logical :: use_modeled_ekin_density = .false.
!
  character(len("sw_rspace_ekin_density")),   private, parameter :: &
       &                      tag_sw_rspace_ekin_density = "sw_rspace_ekin_density"
  character(len("sw_calc_ekin_density_hardpart")),   private, parameter :: &
       &             tag_sw_calc_ekindens_hardpart = "sw_calc_ekin_density_hardpart"
  integer :: sw_rspace_ekin_density = OFF
  integer :: sw_calc_ekin_density_hardpart = OFF
!
  character(len("sw_add_ekin_hardpart_on_Gspace")),   private, parameter :: &
       &             tag_sw_add_ekin_hard_on_Gspace = "sw_add_ekin_hardpart_on_Gspace"
  integer :: sw_add_ekin_hardpart_on_Gspace = OFF

  character(len("initial_kinetic_density")),private,parameter :: tag_initial_kinetic_density = "initial_ekin_dens"
  integer :: initial_ekin_dens = 0
! =========================== 13.0XX, 2014/09/19

  character(len("sw_cif_output")), private, parameter :: &
     & tag_sw_cif_output = "sw_cif_output"
  integer :: sw_cif_output = OFF

  integer :: sw_output_xc_seperately = OFF

! vdW-DF related parameters

  character(len("vdwdf")),   private, parameter :: tag_vdwdf   = "vdwdf"
  character(len("ndel")),    private, parameter :: tag_ndel    = "ndel"
  character(len("nphiD")),   private, parameter :: tag_nphiD   = "nphiD"
  character(len("nr12")),    private, parameter :: tag_nr12    = "nr12"
  character(len("maxk")),    private, parameter :: tag_maxk    = "maxk"
  character(len("r12max")),  private, parameter :: tag_r12max  = "r12max"
  character(len("lambda")),  private, parameter :: tag_lambda  = "lambda"
  character(len("q0cut")),   private, parameter :: tag_q0cut   = "q0cut"
  character(len("q0min")),   private, parameter :: tag_q0min   = "q0min"
  character(len("ds")),      private, parameter :: tag_ds      = "ds"
  character(len("na_gl")),   private, parameter :: tag_na_gl   = "na_gl"
  character(len("a1")),      private, parameter :: tag_a1      = "a1"
  character(len("a2")),      private, parameter :: tag_a2      = "a2"

  character(len("eval_kernel_by_interpolation")), private, parameter ::  &
  & tag_eval_kernel_by_interpolation = "eval_kernel_by_interpolation"
  character(len("mode")),    private, parameter :: tag_mode    = "mode"
  character(len("oneshot")), private, parameter :: tag_oneshot = "oneshot"
  character(len("scf")),    private, parameter :: tag_scf      = "scf"
  character(len("sw_use_WuGygi_method")),   private, parameter :: &
       &            tag_sw_use_WuGygi_method  = "sw_use_WuGygi_method"

  logical :: eval_kernel_by_interpolation = .true.
  integer :: na_gl=30
  real(kind=DP) :: a1=0.d0,a2=60.d0
  real(kind=DP) :: dq_vdw=0.05d0,lambda=1.03d0,q0cut=3.d0, q0min=0.09d0
  Real(kind=DP) :: ds=0.05d0
  integer :: ndel=200,nphiD=1000
  integer :: nr12=3000,nk=1500
!!$  real(kind=DP) :: maxk=10.d0,r12max=30.0d0  ! r12max=30.0 Bohr
  real(kind=DP) :: maxk=10.d0,r12max=56.7d0   ! r12max=56.7d0 Bohr = 30.0 Angstrom
!                Default r12max is changed from "r12max=30.0" bohr by T. Yamasaki,2017/04/14
  logical :: oneshot = .true.
  logical :: sw_save_memory_vdw = .true.
  integer :: sw_use_WuGygi_method = OFF

! ---- vdwdf2
  character(len("vdwdf_version")),    private, parameter :: &
       &                tag_vdwdf_version      = "vdwdf_version"
  character(len("exchange_pot_type")),    private, parameter :: &
       &                tag_exchange_pot_type    = "exchange_pot_type"
  integer, parameter :: len_exchange_pot_type = 10
!!$  character(len=len_xctype) ::  exchange_pot_type = 'revpbe'
  character(len=len_xctype) ::  exchange_pot_type = 'pbe'

  integer :: vdwdf_version = 1
! ----

  logical :: force_exx_energy1=.false.

! ==== SOI
  character(len("sw_write_soi_on_atoms")),    private, parameter :: &
       &               tag_sw_write_soi_on_atoms = "sw_write_soi_on_atoms"
  integer :: sw_write_soi_on_atoms = OFF

! --- Transntion moment related --
!
! === KT_add === 2014/09/24, 2014/10/14
  integer :: ipriepsilon = 1           !! moved from m_Epslilon.F90
  character(len("ipriepsilon")),private,parameter :: tag_ipriepsilon = "ipriepsilon"
!
  integer :: sw_corelevel_spectrum = OFF   !! moved from m_CoreLevel_Spectrum.F90
  integer :: sw_excitation = OFF           !! moved from m_Excitation.F90

  integer :: sw_local_approx_trans_moment = OFF
  integer :: sw_v2c_xes = OFF
! ============== 2014/09/24, 2014/10/14

! --- Crystal Field (experimental)
  character(len("crystal_field")),private,parameter :: &
       &           tag_crystal_field = "crystal_field"
  character(len("sw_print_crystal_field_param")),private,parameter :: &
       &           tag_sw_print_crys_field_param = "sw_print_crystal_field_param"
  integer :: sw_print_crystal_field_param = OFF

! --- gap
  integer :: iprigap = -1
  character(len("iprigap")),private,parameter :: &
       &                      tag_iprigap = "iprigap"

! -- charged defect analysis
  integer :: iprichgdefect = -1
  character(len("iprichgdefect")),private,parameter :: &
       &                      tag_iprichgdefect = "iprichgdefect"

! --- internal coords. in bravais lattice
  integer :: ipribravpos = -1
  character(len("ipribravpos")),private,parameter :: &
       &                      tag_ipribravpos = "ipribravpos"

  integer, save             :: isolver_now

  integer :: sw_harris_functional = OFF
  character(len("sw_harris_functional")),private,parameter :: tag_sw_harris_functional="sw_harris_functional"

! FFT
  integer :: sw_serial_fft = OFF
  character(len("sw_serial_fft")), private, parameter :: tag_sw_serial_fft = "sw_serial_fft"

! convergence criteria
  character(len("convergence_criteria")),private,parameter :: tag_convergence_criteria = "convergence_criteria"
  character(len("delta_energy")),private,parameter :: tag_delta_energy = "delta_energy"
  character(len("delta_moving_average")),private,parameter :: tag_delta_moving_average = "delta_moving_average"
  character(len("slope")),private,parameter :: tag_slope = "slope"
  character(len("delta_v")),private,parameter :: tag_delta_v = "delta_v"

  integer :: nsamp = 10
  character(len("nsamp")), private, parameter :: tag_nsamp = "nsamp"

! use communicator dedicated for charge or not
  integer :: sw_communicator_for_chg = off
  character(len("sw_communicator_for_chg")), private, parameter :: tag_communicator_for_chg = "sw_communicator_for_chg"

! BoltzTraP output
  integer :: sw_boltztrap = off
  character(len("boltztrap")), private,parameter :: tag_boltztrap = "boltztrap"
  character(len("sw_boltztrap")), private, parameter :: tag_boltztrap_output = "sw_boltztrap"
  character(255) :: boltztrap_prefix
  character(len("prefix")), private, parameter :: tag_boltztrap_prefix = "prefix"
  integer :: boltztrap_version = 1
  character(len("version")), private, parameter :: tag_boltztrap_version = "version"
  character(255) :: boltztrap_header
  character(len("header")), private, parameter :: tag_boltztrap_header = "header"

! stress correction
  integer :: sw_stress_correction = off
  character(len("sw_stress_correction")), private, parameter :: tag_stress_correction="sw_stress_correction"
  real(kind=DP) :: decut_stress_correction = 5.d0
  character(len("delta_ecut")), private, parameter :: tag_delta_ecut="delta_ecut"

! chekpoint file output control
  character(len("checkpoint_file")), private, parameter :: tag_checkpoint_file="checkpoint_file"
  character(len("iteration")), private, parameter :: tag_scf_iteration="iteration"
  character(len("iteration_ionic")), private, parameter :: tag_strevl_iteration="iteration_ionic"
  character(len("iteration_unitcell")), private, parameter :: tag_unitcell_iteration="iteration_unitcell"
  character(len("iteration_neb")), private, parameter :: tag_neb_iteration="iteration_neb"
  character(len("iteration_reac")), private, parameter :: tag_reac_iteration="iteration_reac"
  character(len("cputime")), private, parameter :: tag_time="cputime"
  integer :: cpt_scf_iteration=0,cpt_strevl_iteration=10,cpt_neb_iteration=0,cpt_unitcell_iteration=0 &
          & ,cpt_reac_iteration
  real(kind=DP) :: cpt_time=-1.d0
  integer :: cpt_nhistory=1,cpt_history_count=0

  integer :: terminated_because = UNDETERMINED

  logical :: gmaxp_defined

  character(len("sw_write_pwbs_info")), private, parameter :: tag_sw_write_pwbs_info="sw_write_pwbs_info"
  character(len("sw_read_pwbs_info")), private, parameter :: tag_sw_read_pwbs_info="sw_read_pwbs_info"
  integer :: sw_read_pwbs_info = OFF
  integer :: sw_write_pwbs_info = ON

  character(len("sw_modified_kpoint_increment")), private, parameter :: tag_sw_modified_kpoint_increment= &
              & "sw_modified_kpoint_increment"
  integer :: sw_modified_kpoint_increment = OFF

  character(len("sw_output_kpoint_info")), private, parameter :: tag_sw_output_kpoint_info="sw_output_kpoint_info"
  integer :: sw_output_kpoint_info = OFF
  character(len("nnp_output")),    private, parameter :: tag_nnp_output = "nnp_output"
  character(len("sw_nnp_output")), private, parameter :: tag_sw_nnp_output = &
                                                      & "sw_nnp_output"
  character(len("xsf")),  private, parameter :: tag_xsf  = "xsf"
  character(len("n2p2")), private, parameter :: tag_n2p2 = "n2p2"
  character(len("deepmd")), private, parameter :: tag_deepmd = "deepmd"
  character(len("all")),  private, parameter :: tag_all  = "all"

  integer :: sw_nnp_output = OFF
  integer :: filetype_nnp  = ALL
  integer :: frequency_nnp = 100

  logical :: cutoff_wf_changed = .false.
  interface m_CtrlP_rd_val
    module procedure rd_val_real
    module procedure rd_val_int
    module procedure rd_val_str
  end interface m_CtrlP_rd_val

#ifdef __EDA__
  character(len("eda")), private, parameter :: tag_eda = "eda"
  character(len("sw_eda")), private, parameter :: tag_sw_eda = "sw_eda"
  integer :: sw_eda = OFF
#endif

  character(len("torque_convergence")), private, parameter :: tag_torque_convergence = "torque_convergence"
  character(len("max_torque")), private, parameter :: tag_max_torque = "max_torque"
  character(len("max_translational_force")), private, parameter :: tag_max_translational_force = "max_translational_force"

  real(DP) :: max_force_trans  = 1e-3
  real(DP) :: max_torque       = 1e-3


#ifdef KMATH_FFT3D
  character(len("sw_kmath_fft3d")), private, parameter :: tag_sw_kmath_fft3d = "sw_kmath_fft3d"
  character(len("kmath_fft3d")),    private, parameter :: tag_kmath_fft3d    = "kmath_fft3d"
  character(len("nstage")),         private, parameter :: tag_nstage         = "nstage"
  character(len("nstagex")),        private, parameter :: tag_nstagex        = "nstagex"
  character(len("nstagey")),        private, parameter :: tag_nstagey        = "nstagey"
  character(len("nstagez")),        private, parameter :: tag_nstagez        = "nstagez"
  character(len("nprocx")),         private, parameter :: tag_nprocx         = "nprocx"
  character(len("nprocy")),         private, parameter :: tag_nprocy         = "nprocy"
  character(len("nprocz")),         private, parameter :: tag_nprocz         = "nprocz"
  integer :: sw_kmath_fft3d = OFF
  integer, dimension(3) :: nstage_fft3d, nproc_fft3d
#endif

  character(len("sw_optimize_blocking_parameters")), private, parameter :: tag_sw_optimize_blocking_parameters &
  &                                                  = "sw_optimize_blocking_parameters"
  integer                            :: sw_optimize_blocking_parameters = OFF
  integer                            :: nblsizecand_betar, nblsizecand_mgs, nblsizecand_vnonlocal_w &
  &                                   , nblsizecand_submat
  integer, allocatable, dimension(:) :: blsizecand_betar,  blsizecand_mgs, blsizecand_vnonolocal_w &
  &                                   , blsizecand_submat

#ifdef MPI_FFTW
  integer :: sw_mpi_fftw=OFF
  character(len("sw_mpi_fftw")), private, parameter :: tag_sw_mpi_fftw = "sw_mpi_fftw"
#endif

  integer :: sw_keep_hloc_phi = ON
  character(len("sw_keep_hloc_phi")), private, parameter :: tag_sw_keep_hloc_phi='sw_keep_hloc_phi'
  integer :: sw_betar_dot_wfs_exp = OFF
  integer :: sw_precalculate_phase_vnonlocal = OFF
  character(len("sw_betar_dot_wfs_exp")), private, parameter :: tag_sw_betar_dot_wfs_exp='sw_betar_dot_wfs_exp'
  character(len("sw_precalculate_phase_vnonlocal")), private, parameter :: &
                tag_sw_precalculate_phase_vnonlocal='sw_precalculate_phase_vnonlocal'

  integer :: sw_reduce_fft_for_charge = OFF
  character(len("sw_reduce_fft_for_charge")), private, parameter :: tag_sw_reduce_fft_for_charge='sw_reduce_fft_for_charge'

  integer :: configuration_tag = -1

contains
  ! ---- subroutines
  !    ### reading the inputfile in new format style ###
  ! - alloc_w_solver
  ! - alloc_charge_mixing
  ! m_CtrlP_rd_control
  !   -- set_icond
  ! m_CtrlP_rd_accuracy
  !   -- getgmax, getgmaxp, getgmaxs, set_smearing_method, set_intzaj, set_initial_chg
  ! m_CtrlP_rd_wfsolver
  !   -- set_wfsolvers
  ! - f_readsolver
  ! m_CtrlP_rd_struc_evol
  ! m_CtrlP_rd_chargemix
  !   -- set_cdmixingmethods
  ! - f_readchargemixing
  ! m_CtrlP_rd_printlevel
  !   -- set_ipri
  ! - set_printoutlevel_default
  ! - confirm_printoutlevel
  ! m_CtrlP_rd_postproc
  !
  !    ### reading the inputfile in old format style ###
  ! m_CtrlP_rd_parameters
  !   -- read_charge_mixing_detail, read_solver_detail, down2nextst, skiplines_and_read_neg,
  !      read_gmax_gmaxp_natm_ntyp, read_icond_iconstpw, read_ipri, read_nmd1_max_total_scf_iteration_etc,
  !      read_mixinig_parameters, read_charge_precon, read_dtim_1234, read_dtio_imdalg_iexpl_edelta,
  !      read_width_forccr_istress, read_xctype_nspin, read_destm, read_intzaj_imatrix_diagon,
  !      read_gmaxs_or_n_matrix_size, read_imsd, read_evaluation_eko_diff_submat,
  !      read_solver_numbers, read_fine_STM_simulation, check_of_n_WF_solvers, read_imGSrmm,
  !      read_rr_Critical_Value, read_rmm_printout, read_rmm_precal_phase_matm, read_gdiis_hownew
  !
  ! m_CtrlP_set_ekmod_ON
  ! m_CtrlP_ntcnvg_reset
  ! m_CtrlP_ntcnvg_incre
  ! m_CtrlP_ntcnvg_clear
  ! m_CtrlP_pstrn_ntcnvg_reset
  ! m_CtrlP_pstrn_ntcnvg_incre
  ! m_CtrlP_pstrn_ntcnvg_clear
  ! m_CtrlP_set_rmx
  ! m_CtrlP_set_kimg
  ! - alloc_p_WF_solvers
  ! m_CtrlP_wd_isolver
  ! m_CtrlP_rd_isolver
  ! m_CtrlP_check_inputfilestyle
  ! m_CtrlP_wd_iconvergence
  ! m_CtrlP_rd_iconvergence
!!$  ! m_CtrlP_rd_nrsv
!!$  ! m_CtrlP_wd_nrsv
!!$  ! m_CtrlP_rd_nrsv_stdin
  ! m_CtrlP_decide_dtim_1Dsearch
  ! m_CtrlP_dtim_1Dsearch_now
  ! m_CtrlP_dtim_1Dsearch_is_neg
  ! m_CtrlP_reset_dtim_1Dsearch
  ! m_CtrlP_rd_istop
  ! m_CtrlP_set_wct_start
  ! m_CtrlP_ckcput
  ! m_CtrlP_wd_cpu_total
  ! m_CtrlP_set_paramset_on
  ! m_CtrlP_set_paramset_off
  ! m_CtrlP_way_of_smearing
  ! m_CtrlP_waymix_now
  ! m_CtrlP_set_mix_parameter
  ! m_CtrlP_rmx_now
  ! - set_dx_linear
  ! - set_dx_tanh
  ! m_CtrlP_dtim_now
  ! m_CtrlP_On_or_Off_precon_WFs
  ! - On_or_Off_precon_for_WFs_imsd
  ! m_CtrlP_solver_for_WFs_now
  ! m_CtrlP_set_way_ksample
  ! m_CtrlP_set_nspin_and_af
  ! m_CtrlP_set_af
  ! m_CtrlP_what_is_imdalg
  ! m_CtrlP_set_gdiisoptmode
  ! m_CtrlP_set_iconvergence
  ! set_charge_filetype
  !
  subroutine m_CtrlP_set_printable()
    printable = .false.
    if(mype == 0 .or. ipriparadeb /= 0) printable = .true.
  end subroutine m_CtrlP_set_printable

!!$!BRANCH_P 3D_Parallel
#ifdef _USE_SCALAPACK_
  subroutine m_CtrlP_set_sw_scalapack(printable, nfout)
    logical, intent(in) :: printable
    integer, intent(in) :: nfout
    integer :: nptot, ierr
    call mpi_comm_size(mpi_comm_world,nptot,ierr)
    if(nrank_k >=2 .or. nptot .ne. npes) then
       sw_scalapack_md = OFF
       ! This is a tentative default setting until scalapack parallelization is completed for nrank_k>=2
       if(printable) write(nfout,'(a)') 'the default value for parameter sw_scalapack_md is set to off'
    else
       if(printable) write(nfout,'(a," : sw_scalapack = ",i3)') &
                    & 'the default value for parameter sw_scalapack_md is unchanged',sw_scalapack_md
    end if
  end subroutine m_CtrlP_set_sw_scalapack
#endif
!!$!BRANCH_P_END 3D_Parallel

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_check_matm(nfout,natm)
    integer, intent(in) :: nfout,natm
    if(rmm_precal_phase_matm > natm) then
       rmm_precal_phase_matm = natm
       if(printable) &
            & write(nfout,'(" !** rmm_precal_phase_matm(redefined) = ",i10," <<m_CtrlP_check_matm>>")') rmm_precal_phase_matm
    end if
    if(.not.edelta_rmm_given)then
       edelta_change_to_rmm = 1.d-3/dble(natm)
       edelta_change_to_rmm_md = 1.d-3/dble(natm)
    endif
  end subroutine m_CtrlP_check_matm

  subroutine alloc_w_solver(n)
    integer, intent(in) :: n
    integer :: i
    allocate(w_solver(n))
    do i = 1, n
       w_solver(i)%before_or_after_convergence = BEFORE
       if(ekmode == OFF) then
          w_solver(i)%solver = lmMSD
       else
          w_solver(i)%solver = MSD
       end if
       w_solver(i)%till_n_iter = -1
       w_solver(i)%precon = YES
       w_solver(i)%iter_range = 100
       w_solver(i)%variation_way = varLINEAR
       w_solver(i)%cmix_pointer = 1
!       w_solver(i)%subspace_rotation = OFF
       w_solver(i)%subspace_rotation = ON
       w_solver(i)%dtim_s = 0.2d0
       w_solver(i)%dtim_e = 0.2d0
    end do
  end subroutine alloc_w_solver

  subroutine dealloc_w_solver
    if(allocated(w_solver)) deallocate(w_solver)
  end subroutine dealloc_w_solver

  subroutine m_CtrlP_dealloc()
    if(allocated(w_solver)) deallocate(w_solver)
    if(allocated(WF_solver)) deallocate(WF_solver)
    if(allocated(till_n_iteration)) deallocate(till_n_iteration)
    if(allocated(charge_mixing)) deallocate(charge_mixing)
    if(allocated(norbital)) deallocate(norbital)
    if(allocated(l_orb)) deallocate(l_orb)
    if(allocated(t_orb)) deallocate(t_orb)
    if(allocated(rc_orb)) deallocate(rc_orb)
    if(allocated(k_orb)) deallocate(k_orb)

! ========================== KT_mod ======================== 13.0B
!    if(allocated(proj_attribute)) deallocate(proj_attribute)
!    if(allocated(proj_group)) deallocate(proj_group)
!    if(allocated(num_proj_elems)) deallocate(num_proj_elems)
! ========================================================== 13.0B
  end subroutine m_CtrlP_dealloc

! ========================== KT_add ======================== 13.0B
  subroutine m_CtrlP_dealloc_proj()
    if (allocated(proj_attribute)) deallocate(proj_attribute)
    if (allocated(proj_group)) deallocate(proj_group)
    if (allocated(num_proj_elems)) deallocate(num_proj_elems)
  end subroutine m_CtrlP_dealloc_proj
! ========================================================== 13.0B

  subroutine dealloc_charge_mixing()
    if(allocated(charge_mixing)) deallocate(charge_mixing)
  end subroutine dealloc_charge_mixing

  subroutine alloc_charge_mixing(n)
    integer, intent(in) :: n
    integer :: i
    allocate(charge_mixing(n))
    do i = 1, n
!       charge_mixing(i)%mixing_way    = SIMPLE
       charge_mixing(i)%mixing_way    = PULAY
       charge_mixing(i)%iter_range    = 100
       charge_mixing(i)%variation_way = varLINEAR
       charge_mixing(i)%precon        = YES
       charge_mixing(i)%hownew        = RENEW
       charge_mixing(i)%cutoff        = LARGE
!       charge_mixing(i)%istr          = 1
       charge_mixing(i)%istr          = 3
!       charge_mixing(i)%nbxmix        = 0
       charge_mixing(i)%nbxmix        = 15
       charge_mixing(i)%rmxs          = 0.5d0
       charge_mixing(i)%rmxe          = 0.5d0
    end do
  end subroutine alloc_charge_mixing

  subroutine alloc_proj_attribute(num_projector)
    integer, intent(in) :: num_projector
    integer :: i

    allocate(proj_attribute(num_projector))
    do i = 1, num_projector
       proj_attribute(i)%radius = 0.0d0
       proj_attribute(i)%fwhm   = 0.1d0
       proj_attribute(i)%l      = 0
       proj_attribute(i)%t      = 1
       proj_attribute(i)%Ueff   = 0.d0
       proj_attribute(i)%frotate = .false.
       proj_attribute(i)%phi    = 0.d0
       proj_attribute(i)%theta  = 0.d0
       proj_attribute(i)%psi    = 0.d0
! == KT_mod === 2014/06/06
!       proj_attribute(i)%ityp   = 1
       proj_attribute(i)%ityp   = 0
! ============ 2014/06/06
       proj_attribute(i)%strong_correlated = .false.
       proj_attribute(i)%norbital = 0
       proj_attribute(i)%group  = 0
       proj_attribute(i)%ielem  = 0
       proj_attribute(i)%radius_was_defined = .false.
    end do
  end subroutine alloc_proj_attribute
#endif

  subroutine m_CtrlP_rd_control(nfout,contfile_existence, cont3files_existence)
    ! This subroutine sets following parameters
    !
    !   parameter                  tag_name                      default value
    !   ------------------------------------------------------------------
    !   icond                      "condition"                   AUTOMATIC
    !   fixed_charge_k_parallel    "fixed_charge_option"         ALL_at_ONCE
    !   continuation_using_ppdata  "continuation_using_ppdata"   NO
    !   cpumax                     "cpumax"                      86400.d0
    !   max_total_scf_iteration    "max_total_scf_iteration"     HUGE(1)
    !   max_scf_iteration          "max_scf_iteration"           300
    !   max_mdstep                 "max_mdstep"                  10000
    !   sw_use_wfred               "sw_use_wfred"                OFF
    !   ifstop                     "nfstopcheck"                 1
    !   ipriekzaj                  "use_intermediatefile"        1
    !
    integer, intent(in) :: nfout
    logical, intent(in) :: contfile_existence, cont3files_existence
!!$    character(len=FMAXVALLEN) :: rstr
    integer :: iret, f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock
    real(kind=DP) :: dret
    logical :: contfiles, tf

    if(ipriinputfile >= 3) write(nfout,'(" !** << m_CtrlP_rd_control >>")')
    ! --- Control ---
    if( f_selectBlock( tag_control) == 0 ) then
       if(ipriinputfile >= 2) write(nfout,'(" !** -- tag_control = ",a32," --")') tag_control
       if( f_getStringValue( tag_condition, rstr, LOWER) == 0 ) then
          call set_icond(rstr)  ! -> icond
          if(ipriinputfile >= 1) write(nfout, '(" !** condition = ",a30, " icond = ",i8)') rstr, icond
       end if
       if(icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION .or. icond==FIXED_CHARGE_AUTOMATIC) then
          if(f_selectBlock( tag_fixed_charge_option) == 0) then
             if(ipriinputfile >= 1) write(nfout, '(" !** fixed_charge_option is given")')
             if(f_getStringValue(tag_kparallel,rstr,LOWER) == 0) then
                if(trim(rstr) == trim(tag_one_by_one) ) then
                   fixed_charge_k_parallel = ONE_BY_ONE
                else
                   fixed_charge_k_parallel = ALL_at_ONCE
                end if
                if(ipriinputfile >= 1) write(nfout, '(" !** fixed_charge_k_parallel =",a30)') rstr
             else
                fixed_charge_k_parallel = ALL_at_ONCE
             end if
             if(f_getIntValue(tag_sw_modified_kpoint_increment,iret) == 0) then
               sw_modified_kpoint_increment = iret
               write(nfout,'(" !** modified kpoint increment : ",i3)') sw_modified_kpoint_increment
             endif
             iret = f_selectParentBlock()
          end if
          if(ipriinputfile >= 1) &
               &write(nfout, '(" !** fixed_charge_k_parallel =",i5 &
               &  ," :0=ALL_AT_ONCE, 1=ONE_BY_ONE")') fixed_charge_k_parallel
          if(ekmode == ON) then
             fixed_charge_k_parallel = ONE_by_ONE
             if(ipriinputfile >= 1) then
                write(nfout, '(" !** ekmode = ",i8)') ekmode
                write(nfout, '(" !** fixed_charge_k_parallel =",i5 &
                     &  ," :0=ALL_AT_ONCE, 1=ONE_BY_ONE")') fixed_charge_k_parallel
             end if
          end if

       end if
       if( f_getIntValue( tag_continuation_using_ppdata,iret) == 0) continuation_using_ppdata = iret
       if(ipriinputfile >= 1) write(nfout,'(" !** continuation_using_ppdata = ",i5)') continuation_using_ppdata

       iret = f_getRealValue( tag_cpumax, dret, "sec" )
       if(iret == 0 ) cpumax = dret

!!$       if( f_getIntValue( tag_max_iteration, iret ) == 0 ) max_total_scf_iteration = iret
       tf =  f_getIntValue( tag_max_iteration, iret ) == 0
       if(.not.tf) tf = f_getIntValue( tag_max_total_scf_iteration, iret) == 0
       if(tf) then
          max_total_scf_iteration = iret
          max_TS_iteration_is_given = .true.
       end if

       if( f_getStringValue( tag_precision_WFfile, rstr, LOWER) == 0 ) then
          call set_precision_WFfile(rstr) ! -> precision_WFfile
          if(ipriinputfile >= 1) write(nfout, '(" !** precision_WFfile = ",a30 &
               & , ", precision_WFfile = ",i8, " : ",i4," =DP, ",i4," =SP")') rstr,precision_WFfile,DP,SP
       else
          if(ipriinputfile >= 1) write(nfout, '(" !** precision_WFfile = ",i8, " : ",i4," =DP, ",i4," =SP")') precision_WFfile,DP,SP
       end if

       if( f_getIntValue( tag_max_scf_iteration, iret) == 0) then
          max_scf_iteration = iret
          max_scf_iteration_is_given = .true.   ! TY, 18th Aug. 2009
       end if

       if( f_getIntValue( tag_max_scdft_iteration, iret) == 0) then
          max_scdft_iteration = iret
          max_scdft_iteration_is_given = .true.   ! TY, 18th Aug. 2009
       end if

       if( f_getIntValue( tag_max_mdstep, iret) == 0) then
          max_mdstep = iret
          max_mdstep_is_given = .true.
       end if

       if(max_mdstep_is_given .and. max_TS_iteration_is_given ) then
          if( max_mdstep > max_total_scf_iteration ) max_mdstep = max_total_scf_iteration
       end if

       if(max_TS_iteration_is_given) then
          if( max_scf_iteration > max_total_scf_iteration ) max_scf_iteration = max_total_scf_iteration
       end if

! --> T. Yamasaki, 31 Oct 2008
       if( f_getIntValue( tag_sw_use_wfred,iret)==0) sw_use_wfred = iret
       if(ipriinputfile >= 1) write(nfout,'(" !** sw_use_wfred = ",i3)') sw_use_wfred
! <--
! --> T. Yamasaki, 26th Aug. 2009
!     modified by T.Kokubo & D.Fukata, Feb. 2010
       tf = (f_getIntValue( tag_number_of_blocksize,iret)==0)
       if(.not.tf) tf = (f_getIntValue( tag_blocksize, iret)==0)
       if(.not.tf) tf = (f_getIntValue( tag_nblocksize,iret)==0)
       if(.not.tf) tf = (f_getIntValue( tag_nb,iret) == 0)
       if(.not.tf) tf = (f_getIntValue( tag_nblocksize_dgemm,iret)==0)
       if(tf) then
          nblocksize_dgemm = iret
          nblocksize_dgemm_is_given = .true.

          nblocksize_mgs   = nblocksize_dgemm
          recursivesize_mgs   = nblocksize_dgemm
          nblocksize_betar_dot_wfs = nblocksize_dgemm
          nblocksize_betar_dot_wfs_nlmta = nblocksize_dgemm
          nblocksize_vnonlocal_w   = nblocksize_dgemm
          nblocksize_submat        = nblocksize_dgemm
          nblocksize_mgs_is_given       = .true.
          recursivesize_mgs_is_given    = .false.
          nblocksize_betar_is_given     = .true.
          nblocksize_vnonlocal_is_given = .true.
          nblocksize_submat_is_given    = .true.
          nblocksize_submat_latter_is_given    = .true.
#ifndef NO_FORCE_DGEMM
          nblocksize_force         = nblocksize_dgemm
          nblocksize_force_is_given     = .true.
#endif
       end if
       tf = (f_getIntValue( tag_nblocksize_mgs,iret)==0)
       if(tf) then
          nblocksize_mgs = iret
          nblocksize_mgs_is_given = .true.
       end if
       recursivesize_mgs   = nblocksize_mgs
       tf = (f_getIntValue( tag_recursivesize_mgs,iret)==0)
       if(tf) then
          recursivesize_mgs = iret
          recursivesize_mgs_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_betar,iret) == 0)
       if(.not.tf) tf = (f_getIntValue(tag_nblocksize_betar_dot_wfs,iret) == 0)
       if(tf) then
          nblocksize_betar_dot_wfs = iret
          nblocksize_betar_is_given = .true.
       end if
       tf = (f_getIntValue(tag_nblocksize_betar_dot_wfs_nlmta,iret) == 0)
       if(tf) then
          nblocksize_betar_dot_wfs_nlmta = iret
          nblocksize_betar_nlmta_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_betar_dot_wfs_npe,iret) == 0)
       if(.not.tf) tf = (f_getIntValue(tag_nblocksize_betar_dot_wfs_npe,iret) == 0)
       if(tf) then
          nblocksize_betar_dot_wfs_npe = iret
          nblocksize_betar_npe_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_vnonlocal,iret) == 0)
       if(.not.tf) tf = (f_getIntValue(tag_nblocksize_vnonlocal_w,iret) == 0)
       if(tf) then
          nblocksize_vnonlocal_w = iret
          nblocksize_vnonlocal_is_given = .true.
       end if

       tf = (f_getIntValue( tag_nblocksize_vnonlocal_w_nlmta,iret) == 0)
       if(.not.tf) tf = (f_getIntValue(tag_nblocksize_vnonlocal_w_nlmta,iret) == 0)
       if(tf) then
          nblocksize_vnonlocal_w_nlmta = iret
          nblocksize_vnonlocal_w_nlmta_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_submat,iret) == 0)
       if(tf) then
          nblocksize_submat = iret
          nblocksize_submat_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_submat_latter,iret) == 0)
       if(tf) then
          nblocksize_submat_latter = iret
          nblocksize_submat_latter_is_given = .true.
       end if
#ifndef NO_FORCE_DGEMM
       tf = (f_getIntValue( tag_nblocksize_force,iret) == 0)
       if(tf) then
          nblocksize_force = iret
          nblocksize_force_is_given = .true.
       end if
#endif
       if(f_getIntValue(tag_nblocksize_rspace_betar,iret) == 0) nblocksize_rspace_betar = iret
       if(f_getIntValue(tag_nblocksize_rspace_v,iret) == 0) nblocksize_rspace_v = iret
       tf = (f_getIntValue( tag_nblocksize_vnonlocal_w_nlmta,iret) == 0)
       if(.not.tf) tf = (f_getIntValue(tag_nblocksize_vnonlocal_w_nlmta,iret) == 0)
       if(tf) then
          nblocksize_vnonlocal_w_nlmta = iret
          nblocksize_vnonlocal_w_nlmta_is_given = .true.
       end if
       tf = (f_getIntValue(tag_nblocksize_fftw,iret) == 0)
       if(tf) then
          nblocksize_fftw_is_given = .true.
          nblocksize_fftw = iret
       endif
       tf = (f_getIntValue( tag_nblocksize_vnonlocal_w_f,iret) == 0)
!      if(.not.tf) tf = (f_getIntValue(tag_nblocksize_vnonlocal_w_f,iret) == 0)
       if(tf) then
          nblocksize_vnonlocal_w_f = iret
          nblocksize_vnonlocal_w_f_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_gather_f,iret) == 0)
!      if(.not.tf) tf = (f_getIntValue(tag_nblocksize_gather_f,iret) == 0)
       if(tf) then
          nblocksize_gather_f = iret
          nblocksize_gather_f_is_given = .true.
       end if
       tf = (f_getIntValue( tag_nblocksize_submat,iret) == 0)
       if(tf) then
          nblocksize_submat = iret
       end if
       tf = (f_getIntValue( tag_nblocksize_submat_latter,iret) == 0)
       if(tf) then
          nblocksize_submat_latter = iret
       end if
       tf = (f_getIntValue( tag_fftbox_divide_cube,iret) == 0)
       if(tf) then
          fftbox_divide_cube = iret
       end if
       tf = (f_getIntValue( tag_fftbox_3ddiv_1,iret) == 0)
       if(tf) then
          fftbox_3ddiv_1 = iret
       end if
       tf = (f_getIntValue( tag_fftbox_3ddiv_2,iret) == 0)
       if(tf) then
          fftbox_3ddiv_2 = iret
       end if
       tf = (f_getIntValue( tag_fftbox_3ddiv_3,iret) == 0)
       if(tf) then
          fftbox_3ddiv_3 = iret
       end if
       tf = (f_getIntValue( tag_fftbox_div_1,iret) == 0)
       if(tf) then
          fftbox_div_1 = iret
       end if
       tf = (f_getIntValue( tag_fftbox_div_2,iret) == 0)
       if(tf) then
          fftbox_div_2 = iret
       end if
       tf = (f_getIntValue( tag_sw_fft_xzy,iret) == 0)
       if(tf) then
          sw_fft_xzy = iret
       end if

       if(f_getIntValue(tag_communicator_for_chg,iret) == 0) sw_communicator_for_chg = iret
       if (sw_communicator_for_chg == ON)then
          write(nfout,'(" !** sw_communicator_for_chg = on ")')
       endif

       if(ipriinputfile >= 2) then
#ifndef NO_MGS_DGEMM
          write(nfout,'(" !** MGS_DGEMM is defined")')
#endif
#ifndef NO_NONLOCAL_DGEMM
          write(nfout,'(" !** NONLOCAL_DGEMM is defined")')
#endif
#ifndef NO_NONLOCAL_RMM_DGEMM
          write(nfout,'(" !** NONLOCAL_RMM_DGEMM is defined")')
#endif
#ifndef NO_SUBMAT_DGEMM
          write(nfout,'(" !** SUBMAT_DGEMM is defined")')
#endif
#ifndef NO_FORCE_DGEMM
          write(nfout,'(" !** FORCE_DGEMM is defined")')
#endif
#ifndef NO_MATDIAGON_DGEMM
          write(nfout,'(" !** MATDIAGON_DGEMM is defined")')
#endif
          if(nblocksize_dgemm_is_given) write(nfout,'(" !** nblocksize_dgemm = ",i8)') nblocksize_dgemm
          if(nblocksize_mgs_is_given)   write(nfout,'(" !** nblocksize_mgs   = ",i8)') nblocksize_mgs
          if(recursivesize_mgs_is_given)   write(nfout,'(" !** recursivesize_mgs   = ",i8)') recursivesize_mgs
          if(nblocksize_betar_is_given) &
               & write(nfout,'(" !** nblocksize_betar_dot_wfs = ",i8)') nblocksize_betar_dot_wfs
          if(nblocksize_betar_nlmta_is_given) &
               & write(nfout,'(" !** nblocksize_betar_dot_wfs_nlmta = ",i8)') nblocksize_betar_dot_wfs_nlmta
          if(nblocksize_vnonlocal_is_given) &
               & write(nfout,'(" !** nblocksize_vnonlocal_is_given = ",i8)') nblocksize_vnonlocal_w
          if(nblocksize_betar_npe_is_given) &
               & write(nfout,'(" !** nblocksize_betar_dot_wfs_npe = ",i8)') nblocksize_betar_dot_wfs_npe
          if(nblocksize_submat_is_given) &
               & write(nfout,'(" !** nblocksize_submat_is_given = ",i8)') nblocksize_submat
          if(nblocksize_submat_latter_is_given) &
               & write(nfout,'(" !** nblocksize_submat_latter_is_given = ",i8)') nblocksize_submat_latter
#ifndef NO_FORCE_DGEMM
          if(nblocksize_force_is_given) &
               & write(nfout,'(" !** nblocksize_force_is_given = ",i8)') nblocksize_force
#endif
          if(nblocksize_vnonlocal_w_nlmta_is_given) &
               & write(nfout,'(" !** nblocksize_vnonlocal_w_nlmta = ",i8)') nblocksize_vnonlocal_w_nlmta
          if(nblocksize_fftw_is_given) then
             write(nfout,'(" !** nblocksize_fftw = ",i8)') nblocksize_fftw
          else
             write(nfout,'(" !** nblocksize_fftw is not given, nblocksize_fftw = ",i8)') nblocksize_fftw
          end if
          if(nblocksize_vnonlocal_w_f_is_given) &
        &  write(nfout,'(" !** nblocksize_vnonlocal_w_f = ",i8)') nblocksize_vnonlocal_w_f
          if(nblocksize_gather_f_is_given) &
        &  write(nfout,'(" !** nblocksize_gather_f = ",i8)') nblocksize_gather_f
          if(nblocksize_submat > 0) &
        &  write(nfout,'(" !** nblocksize_submat = ",i8)') nblocksize_submat
          if(nblocksize_submat_latter > 0) &
        &  write(nfout,'(" !** nblocksize_submat_latter = ",i8)') nblocksize_submat_latter
          if(fftbox_divide_cube > 0) &
        &  write(nfout,'(" !** fftbox_divide_cube = ",i8)') fftbox_divide_cube
          if(fftbox_3ddiv_1 > 0) write(nfout,'(" !** fftbox_3ddiv_1 = ",i8)') fftbox_3ddiv_1
          if(fftbox_3ddiv_2 > 0) write(nfout,'(" !** fftbox_3ddiv_2 = ",i8)') fftbox_3ddiv_2
          if(fftbox_3ddiv_3 > 0) write(nfout,'(" !** fftbox_3ddiv_3 = ",i8)') fftbox_3ddiv_3
          if(fftbox_div_1 > 0) write(nfout,'(" !** fftbox_div_1 = ",i8)') fftbox_div_1
          if(fftbox_div_2 > 0) write(nfout,'(" !** fftbox_div_2 = ",i8)') fftbox_div_2
          if(sw_fft_xzy > 0) write(nfout,'(" !** sw_fft_xzy = ",i8)') sw_fft_xzy
       end if
! < --

       if( f_getIntValue( tag_nfstopcheck, iret ) == 0) ifstop = iret
       if( f_selectBlock( tag_mpifft) == 0) then
          if( f_getIntValue( tag_ldx, iret) == 0) ldx = iret
          if( f_getIntValue( tag_ldy, iret) == 0) ldy = iret
          if( f_getIntValue( tag_ldz, iret) == 0) ldz = iret
          if(ipriinputfile >= 1) then
             write(nfout,'(" !* tag_mpifft is not found")')
             write(nfout,'(" !* ldx, ldy, ldz = ",3i8)') ldx, ldy, ldz
          end if
          iret = f_selectParentBlock()
       end if
       if( f_getIntValue( tag_cachesize, iret) == 0) ncachesize_given = iret
#ifdef SAVE_FFT_TIMES
       if( f_getIntValue( tag_sw_save_fft,iret) == 0) sw_save_fft = iret
#endif
       if( f_getIntValue( tag_multiple_replica_mode, iret) == 0) multiple_replica_mode = iret
       if( f_getIntValue( tag_multiple_replica_max_iter, iret) == 0) &
      & multiple_replica_max_iteration = iret

! ===== KT_add ==== 2014/07/20
       if ( f_getIntValue( tag_reuse_nfout_for_nfneb, iret) == 0) then
          reuse_nfout_for_nfneb = iret
          write(nfout,*) '!** reuse_nfout_for_nfneb = ', reuse_nfout_for_nfneb
       endif
! ================= 2014/07/20
       if ( f_getIntValue( tag_cntn_bin_paw_format, iret) == 0) then
          cntn_bin_paw_format = iret
          write(nfout,*) '!** cntn_bin_paw_format = ', cntn_bin_paw_format
          cntn_bin_paw_format_is_set = .true.
       endif
       if ( f_getIntValue( tag_sw_write_zaj, iret) == 0) then
          sw_write_zaj = iret
          write(nfout,*) '!** sw_write_zaj = ', sw_write_zaj
       endif
       if ( f_getIntValue( tag_sw_write_zaj_socsv, iret) == 0) then
          sw_write_zaj_socsv = iret
          write(nfout,*) '*** sw_write_zaj_socsv = ', sw_write_zaj_socsv
       endif

#ifndef _EMPIRICAL_
       if( f_getIntValue( tag_sw_ekzaj, iret ) == 0) sw_ekzaj = iret
       if(ipriinputfile >= 1) write(nfout,'(" !** sw_ekzaj = ",i10)') sw_ekzaj
       if ( f_selectBlock( tag_ek ) == 0 ) then
          if( f_getIntValue( tag_use_intermediatefile, iret) == 0) ipriekzaj = iret
          iret = f_selectParentBlock()
       end if

       if( f_getIntValue( tag_use_additional_projector, iret) == 0) then
          sw_use_add_proj = iret
          if(ipriinputfile >= 1) &
               & write(nfout,'(" !** use_additional_projector = ",i5, " : 1=ON, 0=OFF")') sw_use_add_proj
       end if

!  -------------------- T. Hamada 2021.9.23 --------------------- ! UVSOR
       if( f_getStringValue( tag_kpoint_data, rstr, LOWER) == 0 ) then
           call set_sw_kptdata(rstr)
                     if(ipriinputfile >= 1) &
               & write(nfout,'(" !** kpont_data = ",a30, " sw_kptdata = ",i5)') rstr,sw_kptdata
       end if
! --------------------------------------------------------------- ! UVSOR


       if( f_getStringValue( tag_positron, rstr, LOWER) == 0 ) then
          call set_sw_positron(rstr) ! -> sw_positron
          if(ipriinputfile >= 1) &
               & write(nfout,'(" !** positron = ",a30, " sw_positron = ",i5)') rstr,sw_positron
       end if

#ifdef FFTW3
       if(sw_positron /= OFF) then
          if(ldx /= 0) then
             ldx = 0; if(ipriinputfile >= 1) write(nfout,'(" !* ldx = 0 (ldx is reset from (ldx=1)")')
          end if
          if(ldy /= 0) then
             ldy = 0; if(ipriinputfile >= 1) write(nfout,'(" !* ldy = 0 (ldx is reset from (ldy=1)")')
          end if
          if(ldz /= 0) then
             ldz = 0; if(ipriinputfile >= 1) write(nfout,'(" !* ldz = 0 (ldx is reset from (ldz=1)")')
          end if
       end if

       if( f_getIntValue( tag_sw_fef, iret ) == 0) then
          sw_fef = iret
          if(ipriinputfile >= 1) write(nfout,'(" !** sw_fef = ",i5)') sw_fef
       end if

! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!      if( f_getIntValue( tag_sw_pair_vdw, iret ) == 0) then
!         sw_pair_vdw = iret
!         if(ipriinputfile >= 1) write(nfout,'(" !** sw_pair_vdw = ",i5)') sw_pair_vdw
       if( f_getIntValue( tag_sw_vdw_correction, iret ) == 0) then
          sw_vdw_correction = iret
          if(ipriinputfile >= 1) write(nfout,'(" !** sw_vdw_correction = ",i5)') sw_vdw_correction
! ==============================================================================
       end if

! ========================== KT_mod ========================== 13.0B
!!#endif
#else
       if( f_getIntValue( tag_sw_vdw_correction, iret ) == 0) then
          sw_vdw_correction = iret
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_vdw_correction = ",i5)') sw_vdw_correction
          endif
       end if
#endif
! ============================================================ 13.0B

       if( f_getIntValue( tag_sw_dipole_correction, iret ) == 0) then
          sw_dipole_correction = iret
          if(ipriinputfile >= 1) write(nfout,'(" !** sw_dipole_correction = ",i5)') sw_dipole_correction
       end if

       if( f_getIntValue( tag_sw_screening_correction, iret ) == 0) then
          sw_screening_correction = iret
          if(ipriinputfile >= 1) write(nfout,'(" !** sw_screening_correction = ",i5)') sw_screening_correction
       end if

#endif
       if( f_getIntValue(tag_sw_cif_output,iret)==0) then
          sw_cif_output = iret
          if(ipriinputfile>=1) write(nfout,'(" !** sw_cif_output = ",i5)') sw_cif_output
       endif

       if( f_getIntValue(tag_sw_harris_functional,iret)==0) then
         sw_harris_functional = iret
         if(sw_harris_functional == ON .and. ipriinputfile>=1) then
            write(nfout,'(" !** sw_harris_functional = ON")')
         endif
       endif
       if(f_getIntValue(tag_sw_serial_fft,iret) == 0)then
         sw_serial_fft = iret
         if(sw_harris_functional == ON .and. ipriinputfile>=1) then
            write(nfout,'(" !** sw_serial_fft = ON")')
         endif
       endif

       if(f_selectBlock(tag_checkpoint_file) == 0)then
          if(f_getIntValue(tag_scf_iteration,iret) == 0) cpt_scf_iteration = iret
          if(f_getIntValue(tag_strevl_iteration,iret) == 0) cpt_strevl_iteration = iret
          if(f_getIntValue(tag_unitcell_iteration,iret) == 0) cpt_unitcell_iteration = iret
          if(f_getIntValue(tag_reac_iteration,iret) == 0) cpt_reac_iteration = iret
          if(f_getIntValue(tag_neb_iteration,iret) == 0) cpt_neb_iteration = iret
          if(f_getRealValue(tag_time,dret, "sec") == 0) cpt_time = dret
          if(f_getIntValue(tag_nhistory,iret) == 0) cpt_nhistory = iret
       endif
       cpt_history_count = 0
       if (cpt_nhistory<1) cpt_nhistory = 1
       if(ipriinputfile >= 1)then
           if (cpt_scf_iteration>0 .or. cpt_strevl_iteration>0 .or.  cpt_time>0)then
             write(nfout,'(" !** checkpoint file configuration ")')
           endif
           if(cpt_scf_iteration>0) write(nfout,'(" !** output checkpoint file every ",i5," SCF iterations")') &
                                             &  cpt_scf_iteration
           if(cpt_strevl_iteration>0) write(nfout,'(" !** output checkpoint file every ",i5," strevl steps")') &
                                             &  cpt_strevl_iteration
           if(cpt_time>0) write(nfout,'(" !** output checkpoint file every ",f15.0," s")') &
                                             &  cpt_time
           if (cpt_scf_iteration>0 .or. cpt_strevl_iteration>0 .or.  cpt_time>0)then
             write(nfout,'(" !** store ",i5," checkpoint(s)")') cpt_nhistory
           endif
       endif
       if (f_getIntValue(tag_sw_output_kpoint_info,iret)==0) then
           sw_output_kpoint_info = iret
           if(sw_output_kpoint_info==ON .and. printable) write(nfout,'(a)') ' !** sw_output_kpoint_info == ON'
       endif
#ifdef KMATH_FFT3D
       if (f_getIntValue(tag_sw_kmath_fft3d,iret)==0) then
         sw_kmath_fft3d = iret
       endif
       if(sw_kmath_fft3d == ON) then
         nstage_fft3d = 1
         nproc_fft3d = 1;nproc_fft3d(1) =  nrank_g
         if (f_selectBlock(tag_kmath_fft3d)==0) then
           if(f_getIntValue(tag_nstage,iret)==0) then
             nstage_fft3d = iret
           endif
           if(f_getIntValue(tag_nstagex,iret)==0) then
             nstage_fft3d(1) = iret
           endif
           if(f_getIntValue(tag_nstagey,iret)==0) then
             nstage_fft3d(2) = iret
           endif
           if(f_getIntValue(tag_nstagez,iret)==0) then
             nstage_fft3d(3) = iret
           endif
           if(f_getIntValue(tag_nprocx,iret)==0) then
             nproc_fft3d(1) = iret
           endif
           if(f_getIntValue(tag_nprocy,iret)==0) then
             nproc_fft3d(2) = iret
           endif
           if(f_getIntValue(tag_nprocz,iret)==0) then
             nproc_fft3d(3) = iret
           endif
           if (nproc_fft3d(1)*nproc_fft3d(2)*nproc_fft3d(3) .ne. nrank_g) then
             write(nfout,'(a)') &
             &    '!** the product of each element in nproc must be equal to nrank_g'
             call phase_error_with_msg(nfout,'!** the product of each element in nproc must be equal to nrank_g'&
                                      ,__LINE__,__FILE__)
           endif
           iret = f_selectParentBlock()
         endif
       endif
#endif
       if (f_getIntValue(tag_sw_keep_hloc_phi,iret)==0) then
           sw_keep_hloc_phi = iret
           if(sw_keep_hloc_phi==OFF .and. printable) write(nfout,'(a)') ' !** sw_keep_hloc_phi == OFF'
       endif
       if (f_getIntValue(tag_sw_betar_dot_wfs_exp,iret)==0) then
           sw_betar_dot_wfs_exp = iret
           if(sw_betar_dot_wfs_exp==ON .and. printable) write(nfout,'(a)') ' !** sw_betar_dot_wfs_exp == ON'
       endif
       if (f_getIntValue(tag_sw_precalculate_phase_vnonlocal,iret)==0) then
           sw_precalculate_phase_vnonlocal = iret
           if(sw_precalculate_phase_vnonlocal==ON .and. printable) write(nfout,'(a)') &
           ' !** sw_precalculate_phase_vnonlocal == ON'
       endif
       if (f_getIntValue(tag_sw_reduce_fft_for_charge,iret)==0) then
           sw_reduce_fft_for_charge = iret
           if(sw_reduce_fft_for_charge==ON .and. printable) write(nfout,'(a)') &
           ' !** sw_reduce_fft_for_charge == ON'
       endif
       if(f_getRealValue(tag_pcpudf,dret,'')==0) then
         PCPUDF = dret
       endif

       if(f_getIntValue(tag_sw_optimize_blocking_parameters,iret)==0) sw_optimize_blocking_parameters = iret
       if(sw_optimize_blocking_parameters==ON) call set_blsize_opt_scheme()

#ifdef MPI_FFTW
       if(f_getIntValue(tag_sw_mpi_fftw,iret)==0) sw_mpi_fftw = iret
       if(sw_mpi_fftw==ON .and. ipriinputfile>=1 .and. printable) write(nfout,'(a)') ' !** sw_mpi_fftw = ON'
#endif

       iret = f_selectParentBlock()

    else
       if(ipriinputfile >= 1) then
          write(nfout,'(" !* tag_control is not found")')
          write(nfout,'(" !* default values are applied")')
       end if
    end if

    if(ipriinputfile >= 1) then
       write(nfout,'(" !** cpumax = ",f25.8," (sec)")') cpumax
       write(nfout,'(" !** max_iteration = ",i10)') max_total_scf_iteration
       write(nfout,'(" !** nfstopcheck = ",i10)') ifstop
       write(nfout,'(" !** multiple_replica_mode = ",i10)') multiple_replica_mode
       write(nfout,'(" !** sw_ekzaj = ",i10)') sw_ekzaj
       if(ekmode==ON .or. ekmode == GRID) write(nfout,'(" !** ipriekzaj = ",i10)') ipriekzaj
       if(ncachesize_given >= 0) write(nfout,'(" !** ncachesize_given = ",i10)') ncachesize_given
#ifdef SAVE_FFT_TIMES
       write(nfout,'(" !** sw_save_fft = ",i5)') sw_save_fft
#endif
    end if

    if(icond == AUTOMATIC .or. icond == FIXED_CHARGE_AUTOMATIC) then
       if(icond == AUTOMATIC .and. ekmode == OFF) then
          if(.not.contfile_existence) then
             icond = INITIAL
          else
             icond = CONTINUATION
          end if
       else if((icond==AUTOMATIC .and. ekmode==ON) .or. icond==FIXED_CHARGE_AUTOMATIC) then
          contfiles = contfile_existence
          if(icond == FIXED_CHARGE_AUTOMATIC .and.  continuation_using_ppdata == NO) contfiles = cont3files_existence
          if(.not.contfiles) then
             icond = FIXED_CHARGE
          else
             icond = FIXED_CHARGE_CONTINUATION
          end if
       end if
       if(ipriinputfile >= 1) write(nfout, '(" !** condition is changed: icond = ",i8)') icond
    end if

  contains
    subroutine set_icond(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf

      call strncmp2(rstr, FMAXVALLEN, tag_preparation &
           & , len(tag_preparation),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'-2',2,tf)
      if(tf) then
         icond = PREPARATION_ONLY
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_fixed_charge_automatic &
           & , len(tag_fixed_charge_automatic), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '-3', 2, tf)
      if(tf) then
         icond = FIXED_CHARGE_AUTOMATIC
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_coordinate_continuation &
           & , len(tag_coordinate_continuation), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '-4', 2, tf)
      if(tf) then
         icond = COORDINATE_CONTINUATION
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_decision_by_file_existence &
           & , len(tag_decision_by_file_existence),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_automatic &
           & , len(tag_automatic),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'-1',2,tf)
      if(tf) then
         icond = AUTOMATIC
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_fixed_charge_continuation &
           & , len(tag_fixed_charge_continuation), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '3', 1, tf)
      if(tf) then
         icond = FIXED_CHARGE_CONTINUATION
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_initial, len(tag_initial), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '0', 1, tf)
      if(tf) then
         icond = INITIAL
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_continuation, len(tag_continuation), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '1', 1, tf)
      if(tf) then
         icond = CONTINUATION
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_fixed_charge, len(tag_fixed_charge), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '2', 1, tf)
      if(tf) then
         icond = FIXED_CHARGE
         goto 1001
      end if
1001  continue
    end subroutine set_icond

    subroutine set_precision_WFfile(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf

      call strncmp2(rstr, FMAXVALLEN, tag_double_precision, len(tag_double_precision),tf)
      if(.not.tf) call strncmp2(rstr,FMAXVALLEN, tag_DP, len(tag_DP), tf)
      if(tf) then
         precision_WFfile = DP
         if(ipriinputfile >= 1) write(nfout, '(" !** precision_WFfile =  ",i8)') precision_WFfile
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_single_precision, len(tag_single_precision),tf)
      if(.not.tf) call strncmp2(rstr,FMAXVALLEN, tag_SP, len(tag_SP), tf)
      if(tf) then
         precision_WFfile = SP
         if(ipriinputfile >= 1) write(nfout, '(" !** precision_WFfile =  ",i8)') precision_WFfile
         goto 1001
      end if
! -------------------- THamada 2021.9.22 ------------------- ! UVSOR
      call strncmp2(rstr, FMAXVALLEN, tag_null,len(tag_null),tf)
      if(.not.tf) call strncmp2(rstr,FMAXVALLEN, tag_NL, len(tag_NL), tf)
      if(tf) then
         precision_WFfile = 0
         if(ipriinputfile >= 1) write(nfout, '(" !** precision_WFfile =  ",i8)') precision_WFfile
         goto 1001
      end if
!----------------------------------------------------------- ! UVSOR
1001  continue
    end subroutine set_precision_WFfile

!--------------------- T. Hamada 2021.23 ------------------- ! UVSOR
    subroutine set_sw_kptdata(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp0(trim(rstr), tag_no_save, tf)
      if(tf) then
         sw_kptdata = NO_SAVE
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_save, tf)
      if(tf) then
         sw_kptdata = SAVE
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_read, tf)
      if(tf) then
         sw_kptdata = READ
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_prep, tf)
      if(tf) then
         sw_kptdata = PREP
         goto 1001
      end if
1001  continue
    end subroutine set_sw_kptdata
!------------------------------------------------------------- ! UVSOR

#ifndef _EMPIRICAL_
    subroutine set_sw_positron(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp0(trim(rstr), tag_bulk, tf)
      if(tf) then
         sw_positron = BULK;   positron_method = Positron_CONV
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_defect, tf)
      if(tf) then
         sw_positron = DEFECT;   positron_method = Positron_GGGC
         goto 1001
      end if
1001  continue
    end subroutine set_sw_positron
#endif

    subroutine set_blsize_opt_scheme()
      integer :: i
      nblsizecand_betar       = 8
      nblsizecand_mgs         = 8
      nblsizecand_vnonlocal_w = 8
      nblsizecand_submat      = 8
      if(.not.allocated(blsizecand_betar)) allocate(blsizecand_betar(nblsizecand_betar))
      if(.not.allocated(blsizecand_mgs)) allocate(blsizecand_mgs(nblsizecand_mgs))
      if(.not.allocated(blsizecand_vnonolocal_w)) allocate(blsizecand_vnonolocal_w(nblsizecand_vnonlocal_w))
      if(.not.allocated(blsizecand_submat)) allocate(blsizecand_submat(nblsizecand_submat))

      do i=1,nblsizecand_betar
        blsizecand_betar(i) = 4*(2**i)
      enddo

      do i=1,nblsizecand_mgs
        blsizecand_mgs(i) = 4*(2**i)
      enddo

#ifdef SX
      do i=1,nblsizecand_vnonlocal_w
        blsizecand_vnonolocal_w(i) = 512*(2**i)
      enddo
#else
      do i=1,nblsizecand_vnonlocal_w
        blsizecand_vnonolocal_w(i) = 128*(2**i)
      enddo
#endif

      do i=1,nblsizecand_submat
        !blsizecand_submat(i) = neg/i
        !if(blsizecand_submat(i)==0) blsizecand_submat(i) = 1
        blsizecand_submat(i) = 4*(2**i)
      enddo

    end subroutine set_blsize_opt_scheme

  end subroutine m_CtrlP_rd_control

  subroutine m_CtrlP_set_alpha_exx(nfout,alpha)
    integer , intent(in) :: nfout
    real(kind=DP), intent(in) :: alpha
    alpha_exx = alpha
    alpha_exx_is_set_in_SCDFTLoop = .true.
    write(nfout,'(" !!** alpha_exx = ",f8.4)') alpha_exx
  end subroutine m_CtrlP_set_alpha_exx

  subroutine m_CtrlP_rd_accuracy_paw_switch(nfout)
    integer, intent(in) :: nfout
#ifndef _EMPIRICAL_
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock, f_getIntValue
    integer :: iret
    logical :: tf
    iret = f_selectTop()
    ! --- Accuracy ---
    if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** -- tag_accuracy for PAW_switch --")')
    if( f_selectBlock( tag_accuracy) == 0) then
       tf = f_getIntValue( tag_PAW,iret) == 0
       if(.not.tf) tf = f_getIntValue(tag_PAW_switch,iret)==0
       if(tf) PAW_switch = iret
       if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** PAW_switch = ",a3)') on_or_off(PAW_switch)
       iret = f_selectParentBlock()
    end if
#endif
    return
  end subroutine m_CtrlP_rd_accuracy_paw_switch

  subroutine m_CtrlP_rd_accuracy_hubbard_switch(nfout)
    integer, intent(in) :: nfout
#ifndef _EMPIRICAL_
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock, f_getIntValue
    integer :: iret
    logical :: tf
    iret = f_selectTop()
    ! --- Accuracy ---
    if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** -- tag_accuracy for Hubbard_switch --")')
    if( f_selectBlock( tag_accuracy) == 0) then
       if( f_selectBlock( tag_hubbard) == 0) then
          if( f_getIntValue( tag_sw_hubbard, iret) == 0)  sw_hubbard = iret
          iret = f_selectParentBlock()
       end if
       iret = f_selectParentBlock()
    end if
#endif
    return
  end subroutine m_CtrlP_rd_accuracy_hubbard_switch

  subroutine m_CtrlP_rd_accuracy(nfout)
    !  This subroutine sets following parameters
    !      gmax
    !      gmaxp
    !      gmax_positron
    !      neg
    !      way_of_smearing
    !      width
    !      xctype
    !      edelta
    !      mtimes_convergence_scf
    !      ek_max_iteration
    !      evaluation_eko_diff
    !      mtimes_convergence_ek
    !      intzaj
    !      n_matrix_size
    !      eps_solve_Hx_eq_ex
    !      gmaxs_given
    !      initial_chg
    !      nel_Ylm
    !      forccr
    integer, intent(in) :: nfout
!!$    character(len=FMAXVALLEN) :: rstr
    character(len=FMAXUNITLEN) :: unit_f
    integer :: iret, f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop
    real(kind=DP) :: dret
    logical :: tf, prealloc
!!$   integer :: n_projectors

! ================= Added by K. Tagami ==================== 0.1
    logical, save :: prealloc_kt = .false.
! =========================================================

    iret = f_selectTop()
    ! --- Accuracy ---
    if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_accuracy --")')
    if( f_selectBlock( tag_accuracy) == 0) then
#ifndef _EMPIRICAL_
       call getgmax()  ! -> gmax
!!$       call getgmaxp() ! -> gmaxp
       if(getgmaxp()==0) gmaxp = gmax*2.d0 ! -> gmaxp
       if(gmaxp < gmax*2.d0) then
          if(printable) then
             write(nfout,'(" !* gmaxp (= ",d12.4 &
                  & ," ) is smaller than 2*gmax (= ",d12.4," )")') gmaxp, gmax
             write(nfout,'(" !* gmaxp is enlarged to be 2*gmax")')
          end if
          gmaxp = gmax*2.d0
       end if
       if(sw_positron /= OFF) call getgmax_positron()
       if( f_getIntValue( tag_num_bands, iret ) == 0) then
          neg = iret
          neg_is_given = .true.
       end if

       if( f_getStringValue( tag_occupations, rstr,LOWER) == 0) then
          if ( rstr == tag_smearing ) then
             occupations = SMEARING
          else if ( rstr == tag_occup_fixed ) then
             occupations = FIXED
             way_of_smearing = LOWEST_AT_EACH_KPT
             write(nfout,*) '!** occupations are fixed (lowest energy, zero smearing)'
          else if ( rstr == tag_file ) then
             occupations = FILE
          endif
       endif

       if ( occupations == SMEARING .and. f_selectBlock( tag_smearing) == 0 ) then
          if( f_getStringValue( tag_smearing_method, rstr,LOWER) == 0) then
             call set_smearing_method(rstr) ! way_of_smearing
             write(nfout,*) "!** smearing method is set to", way_of_smearing
          end if
          if(way_of_smearing == MP) then
             width = 0.01
          endif
          if( f_getRealValue( tag_smearing_width, dret, "hartree") == 0) then
             width = dret
             if(way_of_smearing == TETRAHEDRON) width_tetra = width
          endif
          if( f_getRealValue( tag_electronic_temp, dret, "k") == 0) then
             way_of_smearing = FERMI_DIRAC;     width = dret *CONST_kB
             write(nfout,'(A,E20.10,A)') ' ! --- smearing with (Fermi-Dirac): ', &
                  &                          width, " Hartree"
          endif
          if( f_selectBlock( tag_tetrahedron) == 0) then
             if( f_getIntValue( tag_dimension, iret ) == 0) idimtetra = iret
             if( f_getIntValue( tag_sw_correction, iret ) == 0) sw_correction = iret
             iret = f_selectParentBlock()
          end if
          tf = f_selectBlock(tag_methfessel_paxton)==0
!          if(.not.tf) tf = f_selectBlock(tag_mp)==0
          if(.not.tf) tf = f_selectBlock(tag_meth)==0
          if(tf)then
             if( f_getIntValue( tag_order, iret ) == 0) then
                if(iret>=0) then
                   order_mp = iret
                else
                   write(nfout,'(a)') &
                        & ' !** WARNING : order for methfessel-paxton smearing must not be &
                        & negative; using the default value '
                endif
             endif
             if(f_getIntValue(tag_esearch,iret)==0) esearch = iret
             if(f_getRealValue(tag_esearch_factor,dret,'')==0) then
                esearch_factor_mp = dret
             endif
             iret = f_selectParentBlock()
          endif
          if(sw_correction == ON) then
! ==================== added by K. Tagami ===================== 12.0
             write(nfout,*) '=== sw_correction is set ON'
! ============================================================== 12.0
             idimtetra = -idimtetra
          end if

! === KT_add == 13.0U3
          iret = f_getStringValue( tag_method_change_width,rstr,LOWER )
          if ( rstr == tag_stepwise ) then
             method_change_smearing_width = STEPWISE
             write(nfout,*) '** method_change_smearing_width is set to ', &
                  &         method_change_smearing_width
          endif
          if ( method_change_smearing_width > 0 ) then
             if( f_getRealValue( tag_width_initial, dret, "hartree") == 0) then
                width_initial = dret
                write(nfout,*) '** width_initial is set to ', width_initial
             endif
             if( f_getIntValue( tag_num_intermid_width, iret ) == 0) then
                num_intermid_width = iret
                write(nfout,*) '** num_intermid_width is set to ', num_intermid_width
             endif
             if( f_getRealValue( tag_edelta_change_width_first,dret,"hartree")== 0) then
                edelta_change_width_first = dret
                write(nfout,*) '** edelta_change_width_first is set to ', &
                     &             edelta_change_width_first
             endif
             if( f_getRealValue( tag_edelta_change_width_last,dret, "hartree")== 0) then
                edelta_change_width_last = dret
                write(nfout,*) '** edelta_change_width_last is set to ', &
                     &             edelta_change_width_last
             endif
          endif
! ==================== 13.0U3

          iret = f_selectParentBlock()
       end if

! ========== KT_add ================= 13.0E
       if( f_getRealValue( tag_electronic_temp, dret, "k") == 0) then
          way_of_smearing = FERMI_DIRAC
          width = dret *CONST_kB
          write(nfout,'(A,E20.10,A)') ' ! --- smearing with (Fermi-Dirac): ', &
               &                          width, " Hartree"
       endif
! =================================== 13.0E

! --> T. Yamasaki, 2011/03/01
       fftsize_factor_gmaxp = 1.d0
       if( f_selectBlock( tag_fftsize ) == 0) then
          if(f_getStringValue( tag_factor_for_ChargeDensity, rstr, LOWER) == 0) then
!!$             if(printable) write(nfout,'(" rstr = ",a)') trim(rstr)
             if(rstr == tag_full) then
                fftsize_factor_gmaxp = 1.d0
             else if(rstr == tag_large) then
                fftsize_factor_gmaxp = 1.d0
             else if(rstr == tag_small) then
                fftsize_factor_gmaxp = 0.75d0
             else if(rstr == tag_tiny) then
                fftsize_factor_gmaxp = 0.5d0
             else
                if(f_getRealValue(tag_factor_for_ChargeDensity, dret,"") == 0) then
                   fftsize_factor_gmaxp = dret
                else
                   fftsize_factor_gmaxp = 1.d0
                end if
             end if
          end if
          if(1.d0 <fftsize_factor_gmaxp) fftsize_factor_gmaxp = 1.d0
          if(fftsize_factor_gmaxp < 0.5d0) fftsize_factor_gmaxp = 0.5d0
          iret = f_selectParentBlock()
       end if
       gmaxp_reduced = gmaxp * fftsize_factor_gmaxp
       if(printable) write(nfout,'(" !** fftsize_factor_gmaxp = ",f8.4)') fftsize_factor_gmaxp
!!$       if(gmaxp_reduced < gmax*2.d0) then
!!$          gmaxp_reduced = gmax*2.d0
!!$          if(printable) write(nfout,'(" !** fftsize_factor_gmax(redifined) = ",f8.4)') gmaxp_reduced/gmaxp
!!$       end if
       if(printable) write(nfout,'(" !** gmaxp_reduced        = ",f8.4)') gmaxp_reduced

! <--
       if( f_getStringValue( tag_xctype, rstr,LOWER) == 0) xctype = rstr(1:len_xctype)
! <--
       if( f_getStringValue( tag_xctype, rstr,LOWER) == 0) xctype = rstr(1:len_xctype)
! --> T. Yamasaki, 2010/05/21, 06/08

! ---test --
       if( f_getStringValue( tag_xctype, rstr,LOWER) == 0) then
#ifdef LIBXC
          if ( xctype == "libxc" ) then
             if( f_selectBlock( tag_libxc ) == 0 ) then
                if ( f_getIntValue( tag_exch_id, iret ) == 0 ) then
                   libxc_exch_id = iret
                endif
                if ( f_getIntValue( tag_corr_id, iret ) == 0 ) then
                   libxc_corr_id = iret
                endif
                if ( f_getStringValue(tag_exch_name,rstr,LOWER)==0 ) then
                   libxc_exch_name = rstr
                endif
                if ( f_getStringValue(tag_corr_name,rstr,LOWER)==0 ) then
                   libxc_corr_name = rstr
                endif
                iret = f_selectParentBlock()
             end if
          else
             call check_xctype_2016( nfout, rstr, xctype, &
                  &                  exchange_pot_type, vdwdf_version )
          endif
#else
          call check_xctype_2016( nfout, rstr, xctype, exchange_pot_type, vdwdf_version )
#endif
       endif

       tf = f_getIntValue( tag_PAW,iret) == 0
       if(.not.tf) tf = f_getIntValue(tag_PAW_switch,iret)==0
       if(tf) PAW_switch = iret
       if(PAW_switch == ON .and. f_selectBlock(tag_paw_density_gradient) == 0) then
          tf = f_getStringValue(tag_interpolation_method,rstr,LOWER)==0
          if(.not.tf) tf = f_getStringValue(tag_method,rstr,LOWER)==0
          if(tf) then
             if(rstr == tag_newton) then
                paw_density_gradient%interpolation_method = NEWTON
             else if(rstr == tag_lagrange) then
                paw_density_gradient%interpolation_method = LAGRANGE
             end if
          end if
          if(f_getIntValue(tag_order,iret)==0) paw_density_gradient%order = iret
          if(paw_density_gradient%order .le. 0) paw_density_gradient%order = 1
          iret = f_selectParentBlock()
       end if
! <--
! ==============================================================================
       if( f_getIntValue( tag_ggacmp_parallel, iret) == 0) ggacmp_parallel = iret

       if( f_selectBlock( tag_sc_dft) == 0 .or. f_selectBlock( tag_scdft ) == 0) then
          if( f_getRealValue(tag_delta_epsilon, dret,"") == 0) then
             delta_epsilon = dret
             if(printable) then
                write(nfout,'(1x," !** tag block of sw_dft is set")')
                write(nfout,'(1x," !** delta_epsilon = ",f8.4)') delta_epsilon
             end if
          end if
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock( tag_hybrid_functional) == 0) then
          sw_output_hybrid_info = ON
          if( f_getIntValue( tag_sw_hybrid_functional, iret ) == 0) sw_hybrid_functional = iret
          if( f_getIntValue(tag_sw_output_hybrid_info, iret) == 0) sw_output_hybrid_info = iret
          if(sw_hybrid_functional == ON) then
             if( PAW_switch == ON )then
               if(printable) write(nfout,'(a)') &
               ' !** PAW method cannot be used in conjunction with the hybrid functional.'
               call phase_error_with_msg(nfout, &
               ' !** PAW method cannot be used in conjunction with the hybrid functional.' &
               ,__LINE__,__FILE__)
             endif
             if( f_getIntValue( tag_sw_exchange_only, iret ) == 0) sw_exchange_only = iret
             if( f_getIntValue( tag_sw_singular_correction, iret ) == 0) sw_singular_correction = iret
             if( f_getIntValue( tag_sw_screened_exchange, iret ) == 0) sw_screened_exchange = iret
             !!if( f_getIntValue( tag_sw_memory_reduction, iret ) == 0) sw_memory_reduction_exx = iret
             if( f_selectBlock( tag_reduction_factor) == 0) then
                if( f_getIntValue( tag_f1, iret ) == 0) reduction_factor_exx(1) = iret
                if( f_getIntValue( tag_f2, iret ) == 0) reduction_factor_exx(2) = iret
                if( f_getIntValue( tag_f3, iret ) == 0) reduction_factor_exx(3) = iret
                iret = f_selectParentBlock()
             end if
             if(f_getStringValue(tag_functional_type,rstr,LOWER)==0)then
                if (rstr.eq.tag_hf) then
                   write(nfout,'(a)') ' !** functional_type : HF (Hartree-Fock)'
                   hybrid_functional_type = 'hf'
                   alpha_exx = 1.d0
                   sw_screened_exchange = OFF;    sw_exchange_only = ON
                else if(rstr.eq.tag_pbe0) then
                   write(nfout,'(a)') ' !** functional_type : PBE0'
                   hybrid_functional_type = 'pbe0'
                   alpha_exx = 0.25d0
                   sw_screened_exchange = OFF;   sw_exchange_only = OFF
                else if (rstr.eq.tag_hse06_hjs) then
                   write(nfout,'(a)') ' !** functional_type : HSE06 (HJS)'
                   hybrid_functional_type = 'hse06-hjs'
                   alpha_exx = 0.25d0;   omega_exx = 0.106d0
                   sw_screened_exchange = ON;    sw_exchange_only = OFF
                else if (rstr.eq.tag_hse06) then
                   write(nfout,'(a)') ' !** functional_type : HSE06'
                   hybrid_functional_type = 'hse06'
                   alpha_exx = 0.25d0;   omega_exx = 0.106d0
                   sw_screened_exchange = ON;    sw_exchange_only = OFF
                else if (rstr.eq.tag_hsesol_hjs) then
                   write(nfout,'(a)') ' !** functional_type : HSESOL (HJS)'
                   hybrid_functional_type = 'hsesol-hjs';  xctype = 'pbesol'
                   alpha_exx = 0.25d0;   omega_exx = 0.106d0
                   sw_screened_exchange = ON;    sw_exchange_only = OFF
                else if (rstr.eq.tag_gaupbe) then
                   write(nfout,'(a)') ' !** functional_type : Gau-PBE'
                   hybrid_functional_type = 'gaupbe'
                   sw_screened_exchange = ON;     sw_exchange_only = OFF
                else
                   write(nfout,'(a)') ' !** WARNING : invalid functional_type : '//trim(rstr)
                endif
             endif

             if( f_getRealValue( tag_alpha, dret, "") == 0) then
                if(.not.alpha_exx_is_set_in_SCDFTLoop)  alpha_exx = dret
             end if
             if( f_getRealValue( tag_omega, dret, "") == 0) omega_exx = dret
             if( f_getRealValue( tag_omega_hf, dret, "") == 0)  omega_exx = dret
             if( f_getRealValue( tag_omega_pbe, dret, "") == 0) omega_exx_pbe = dret

             if( f_getIntValue(tag_sw_rspace,iret)==0 ) sw_rspace_hyb = iret
             if( f_getIntValue(tag_sw_rspace_dgm,iret)==0 ) sw_rspace_hyb_dgm = iret
             if( f_getIntValue(tag_sw_eval_vexx,iret)==0 ) sw_eval_vexx = iret
             if( f_getIntValue(tag_sw_retard_eigval_evaluation,iret)==0 ) sw_retard_eigval_evaluation = iret
             if( f_getIntValue(tag_sw_precalculate,iret)==0 ) sw_precalculate = iret
             if( f_getRealValue(tag_r0_factor,dret, "")==0 ) r0_factor_q = dret
             if(f_getStringValue(tag_qr_optimization,rstr,LOWER)==0)then
                if(rstr==tag_mask_function)then
                   qr_optimization = MASK_FUNCTION
                elseif(rstr==tag_prefitting)then
                   qr_optimization = PREFITTING
                elseif(rstr==tag_none)then
                   qr_optimization = NO
                else
                   write(nfout,'(a)') 'unsupported projector optimization method : '//trim(rstr)
                   call phase_error_with_msg(nfout,'unsupported projector optimization method : ' &
                   //trim(rstr),__LINE__,__FILE__)
                endif
             endif

             if(f_getIntValue(tag_sw_change_axis,iret)==0)then
                sw_change_axis = iret
             endif

             if(f_getStringValue(tag_accuracy,rstr,LOWER)==0)then
                if(rstr==tag_exact)    accuracy_of_exx = EXACT
                if(rstr==tag_fine)     accuracy_of_exx = FINE
                if(rstr==tag_moderate) accuracy_of_exx = MODERATE
                if(rstr==tag_coarse)   accuracy_of_exx = COARSE
             endif

             select case(accuracy_of_exx)
             case(EXACT)
               nmax_G_hyb = 0
               gmax_exx = gmax
             case(FINE)
               nmax_G_hyb = -1
               gmax_exx = gmax*dsqrt(0.5d0)
             case(MODERATE)
               nmax_G_hyb = -2
               gmax_exx = gmax*dsqrt(0.5d0)
             case(COARSE)
               nmax_G_hyb = -3
               gmax_exx = gmax*dsqrt(0.25d0)
             end select
! ============================= KT_Test ============================ 12.5Exp
             if( f_getIntValue( tag_nmax_G_hyb, iret ) == 0 ) then
                nmax_G_hyb = iret
             else if (f_getStringValue(tag_charge_mesh,rstr,LOWER) == 0) then
                if(rstr==tag_exact)    nmax_G_hyb =  0
                if(rstr==tag_fine)     nmax_G_hyb = -1
                if(rstr==tag_moderate) nmax_G_hyb = -2
                if(rstr==tag_coarse)   nmax_G_hyb = -3
             endif
             if(sw_change_axis==OFF) nmax_G_hyb=0
             if(printable)then
                if(nmax_G_hyb==0)  write(nfout,'(a)') ' !** charge mesh for the hybrid functional calculations : exact'
                if(nmax_G_hyb==-1) write(nfout,'(a)') ' !** charge mesh for the hybrid functional calculations : fine'
                if(nmax_G_hyb==-2) write(nfout,'(a)') ' !** charge mesh for the hybrid functional calculations : moderate'
                if(nmax_G_hyb==-3) write(nfout,'(a)') ' !** charge mesh for the hybrid functional calculations : coarse'
             endif
! ------
             if( f_getIntValue( tag_truncate_vxw_updating, iret ) == 0 ) then
                if(iret==ON) truncate_vxw_updating = .true.
             endif
             if( f_getRealValue( tag_delta_total_energy_hyb1, dret,"hartree" )== 0) &
                  &                 edelta_for_hyb_chgfix = dret
             if( f_getRealValue( tag_delta_total_energy_hyb2, dret,"hartree" )== 0) &
                  &                 edelta_for_hyb_convgd = dret
! =================================================================== 12.5Exp

! ======= KT_add ======================================== 13.0F
             use_hybrid_functional = .true.
             hybrid_calc_is_active = .true.

             !gmax_exx = gmax;   gmaxp_exx = gmaxp
             gmaxp_exx = gmaxp
             if(nmax_G_hyb==-1) gmaxp_exx = gmax*2.d0
             if(nmax_G_hyb==-2) gmaxp_exx = gmax*4.d0**(1.d0/3.d0)
             if(nmax_G_hyb==-3) gmaxp_exx = gmax
! -----
             if ( f_getRealValue( tag_gmax_exx_ratio, dret, "") == 0 ) then
                if( dret < 0.2D0 .or. dret > 1.0D0 ) then
                   write(nfout,*) ' !! Out of range, so, gmax_ratio is set to 1.0'
                   gmax_exx = gmax
                else
                   gmax_exx = gmax *dret
                end if
             end if
             if ( f_getRealValue( tag_cutoff_wf_for_exx, dret, "rydberg") == 0 ) then
                if ( dret > 0.0 ) then
                   gmax_exx = sqrt( dret )
                endif
                if( gmax_exx < gmax*0.2D0 .or. gmax_exx > gmax ) then
                   write(nfout,*) ' !! gmax_exx is too small or large, so is set to gmax'
                   gmax_exx = gmax
                endif
             end if
! -----
             !if ( gmax_exx < gmax ) use_fft_exx = .true.
             if (sw_hybrid_functional == ON) use_fft_exx = .true.
! -----
             if ( f_getRealValue( tag_gmaxp_exx_ratio, dret, "") == 0 ) then
                if ( dret > 0.0 ) then
                   gmaxp_exx = gmaxp *dret
                endif
                if( gmaxp_exx < gmax ) then
                   write(nfout,*) ' !! gmaxp_exx is too small, so is set to gmax'
                   gmaxp_exx = gmax
                else if ( gmax_exx > gmaxp ) then
                   write(nfout,*) ' !! gmaxp_exx is too large, so is set to gmaxp'
                   gmaxp_exx = gmaxp
                endif
             end if
             if ( f_getRealValue( tag_cutoff_cd_for_exx, dret, "rydberg") == 0 ) then
                if ( dret > 0.0 ) then
                   gmaxp_exx = sqrt( dret )
                endif
                if( gmaxp_exx < gmax ) then
                   write(nfout,*) ' !! gmaxp_exx is too small, so is set to gmax'
                   gmaxp_exx = gmax
                else if ( gmax_exx > gmaxp ) then
                   write(nfout,*) ' !! gmaxp_exx is too large, so is set to gmaxp'
                   gmaxp_exx = gmaxp
                endif
             end if
!
             if( f_getRealValue( tag_edelta_change_to_hybrid, dret,"hartree" )== 0) then
                edelta_change_to_hybrid = dret
                if ( edelta_change_to_hybrid > 0.0 ) then
                   hybrid_calc_is_active = .false.
                endif
             endif
! =================================================== 13.0F

             if(f_getStringValue(tag_potential_update,rstr,LOWER)==0) then
                if(rstr.eq.tag_always)   potential_update = 0
                if(rstr.eq.tag_moderate) potential_update = 1
                if(rstr.eq.tag_minimal)  potential_update = 2
             endif

             if(ipriinputfile>=1) then
                write(nfout,'(" <<< Hybrid functional method >>>")')
                if(sw_exchange_only==ON) &
                 & write(nfout,'(" Use of Exchage only. The correlation part will not be included.")')
                write(nfout,'(" Mixing parameter for exact/screened exhange (alpha):",f10.5)') alpha_exx
                if(sw_screened_exchange==ON) &
                 & write(nfout,'(" Screening parameter for screened exhange (omega):",f10.5)') omega_exx
                if(sw_singular_correction==ON) then
                   write(nfout,'(" Singularity of the exhange potential in Fourier space will be corrected.")')
                else
                   write(nfout,'(" G=0 part of the exhange potential in Fourier space will be zero.")')
                end if
                write(nfout,'(" reduction_factor = ",3i5)') reduction_factor_exx(1:3)
                !!write(nfout,'(" sw_memory_reduction_exx = ",i5)') sw_memory_reduction_exx
! ======= KT_add ====================================== 13.0F
                write(nfout,*) ' ------- '
                write(nfout,'(" !** gmax_exx = ",f8.4)') gmax_exx
                write(nfout,'(" !** gmaxp_exx = ",f8.4)') gmaxp_exx

                if ( hybrid_calc_is_active ) then
                   write(nfout,*) 'hybrid_calc is active'
                else
                   write(nfout,*) 'edelta_change_to_hbrid is set to ', dret
                   write(nfout,*) 'hybrid_calc is inactive first'
                endif
                write(nfout,*) ' ------- '
! ===================================================== 13.0F
                if(sw_change_axis == ON.and.ipriinputfile>=1)then
                   write(nfout,'(a)') ' !** parallelization axis for the inner-most loop will be changed &
                    & from G space to B space for the exact-exchange calculation'
                endif
             end if
          end if
          iret = f_selectParentBlock()
       end if

! === KT_add ==== 13.0Y
       if ( f_selectBlock( tag_pcc ) == 0 .or. &
            &  f_selectBlock( tag_partial_core_correction ) == 0 ) then

          if ( f_getIntValue( tag_sw_eval_epc_on_fftmesh, iret ) == 0 ) then
             sw_eval_epc_on_fftmesh = iret
             write(nfout,*) '!** sw_eval_epc_on_fftmesh is set to ', iret
          endif
          if ( f_getIntValue( tag_sw_remove_pcc_from_pawpot, iret ) == 0 ) then
             sw_remove_pcc_from_pawpot = iret
             write(nfout,*) '!** sw_remove_pcc_from_pawpot is set to ', iret
          endif
          iret = f_selectParentBlock()
       endif
! =============== 13.0Y

       if( f_selectBlock( tag_projector_list) == 0) then
          iret = f_getStringValue(tag_projector_type,rstr,LOWER)
          if( rstr == tag_spherical_harmonics) then
             projector_type = SPHERICAL_HARMONICS
          else if( rstr == tag_atomic_orbital) then
             projector_type = ATOMIC_ORBITAL
          else if( rstr == tag_wannier_function) then
             projector_type = WANNIER_FUNCTION
          end if
          if(projector_type == ATOMIC_ORBITAL) sw_orb_popu = ON

! ==== KT_add == 2014/06/05
          if ( f_getIntValue( tag_sw_allow_maxprojs_gt_4, iret ) == 0 ) then
             sw_allow_maxprojs_gt_4 = iret
             write(nfout,*) '!** sw_allow_maxprojs_gt_4 is set to ', iret
          endif
! ============== 2014/06/05
          if ( f_getIntValue( tag_howto_set_proj_radius, iret ) == 0 ) then
             howto_set_proj_radius = iret
             write(nfout,*) '!** howto_set_proj_radius is set to ', iret
          endif

          prealloc = .true.
          call m_CtrlP_set_projectors(prealloc,num_projectors)

! ==================== Added by K. Tagami =============== 0.1
          if ( .not.prealloc_kt ) then
! ==================================================================

!--- allocation and initialization of proj_attribute
!    revised by T. Yamasaki, 17th Feb. 2006
!!$          allocate(proj_attribute(num_projectors))
          call alloc_proj_attribute(num_projectors)
!<---
          prealloc = .false.
          call m_CtrlP_set_projectors(prealloc,num_projectors)

! ======================== Added by K. Tagami =========== 0.1
            prealloc_kt = .true.
         endif
! =================================================================

          iret = f_selectParentBlock()
       end if

       if( f_selectBlock( tag_hubbard) == 0) then

          if( f_getStringValue(tag_initial_occmat, rstr,LOWER) == 0) call set_initial_occmat(rstr)
          if( f_getRealValue(tag_initial_occmat_factor, dret,"") ==0) initial_occmat_factor=dret
! ==================================== added by K. T. ====== 11.0
          if( f_getIntValue( tag_occmat_diag_mode,iret ) == 0 ) then
             occmat_diag_mode = iret              ! only meaningful in noncol system.
          end if
          if( f_getIntValue( tag_occmat_file_format,iret ) == 0 ) then
             occmat_file_format = iret           ! only meaningful in noncol system.
          end if
! ========================================================== 11.0

          if( f_getIntValue( tag_sw_hubbard, iret) == 0)  sw_hubbard = iret
          if(sw_hubbard == ON) then
             if( f_getRealValue( tag_critical_hubbard_energy, dret, "hartree") == 0) critical_ehub = dret
             if( f_getRealValue( tag_delta_hubbard_energy, dret, "hartree") == 0) delta_ehub = dret
             if( f_getRealValue( tag_amix, dret, "hartree") == 0) alpha_hubbard = dret
             if( f_getIntValue( tag_sw_force_simple_mixing,iret)==0) sw_force_simple_mixing_hub=iret
             if( f_getIntValue(tag_nUeff,iret) == 0) nUeff = iret
             call m_CtrlP_set_hubbard_proj(nfout)
             if( f_selectBlock( tag_constraint) == 0) then
                if( f_getIntValue( tag_sw_constraint, iret) == 0)  sw_constraint = iret
                if( f_getIntValue( tag_site, iret) == 0)  const_site = iret
                if( f_getRealValue( tag_alpha, dret, "hartree") == 0)  const_alpha = dret
                iret = f_selectParentBlock()
             end if
             if(f_getIntValue(tag_occ_matrix_fix_period,iret) == 0) occ_matrix_fix_period = iret
! =============================== added by K. Tagami ===========  5.0
             if( f_getIntValue( tag_sw_orbital_cut_smooth, iret) == 0 ) then
               sw_orbital_cut_smooth = iret
             endif
             write(nfout,*) '!** sw_orbital_cut_smooth = ', sw_orbital_cut_smooth
!
             iret = f_getStringValue( tag_occmat_type,rstr,LOWER )
             if ( rstr == tag_occmat_type1 ) then
                occmat_type = OCCMAT_Type1
             else if ( rstr == tag_occmat_type2 ) then
                occmat_type = OCCMAT_Type2
             endif
             write(nfout,*) '!** occmat type = ', occmat_type

!
             iret = f_getStringValue( tag_dftu_type,rstr,LOWER )
             if ( rstr == tag_fll ) then
                dftu_type = FLL
                write(nfout,*) '!** dftu_type = FLL (Full Localization Limit)'
             else if ( rstr == tag_amf ) then
                dftu_type = AMF
                write(nfout,*) '!** dftu_type = AMF (Around Mean Field)'
             endif

             if ( sw_constraint == ON ) dftu_type = 0
             write(nfout,*) '!** initial_occmat = ', initial_occmat

! -------------------------------------------------
             if( f_selectBlock( tag_Ueff_imposition )==0 )then
                write(nfout,*) '!** tag_Ueff_imposition is found'
                if ( f_getStringValue( tag_method_Ueff_imposition, rstr, &
        &                              LOWER ) == 0 ) then
!
                  call strncmp0( tag_Ueff_gradually, trim(rstr), tf )
                  if ( tf ) then
                     method_Ueff_imposition = Ueff_Gradually
                     write(nfout,*) '!** tag_Ueff_Gradualy is found'
                  endif

                  call strncmp0( tag_Ueff_from_first, trim(rstr), tf )
                  if ( tf ) then
                     method_Ueff_imposition = Ueff_From_First
                     write(nfout,*) '!** tag_Ueff_from_first is found'
                  endif
!
                  if ( method_Ueff_imposition == Ueff_Gradually ) then
                     if ( f_getIntValue( tag_Ueff_transition_period, iret ) == 0 ) then
                        Ueff_transition_period = iret
                        write(nfout,*) '!** Ueff_transition_period is set to ', iret
                     else
                        write(nfout,*) '!** Ueff_transition_period is set to default :',&
        &                               Ueff_transition_period
                     endif

                     if ( f_getIntValue( tag_iter_for_Ueff_starting, iret ) == 0 ) then
                        iteration_for_Ueff_starting = iret
                        write(nfout,*) '!**  iteration_for_Ueff_starting is set to ', &
        &                             iteration_for_Ueff_starting

                        sw_eval_Ueff_using_iter = ON
                     else if ( f_getRealValue( tag_edelta_for_Ueff_starting, &
        &                                     dret,'hartree' )==0 ) then
                        edelta_for_Ueff_starting = dret
                        write(nfout,*) '!** edelta_for_Ueff_starting is set to ', dret

                        sw_eval_Ueff_using_edelta = ON
                     else
                        write(nfout,*) '!** edelta_for_Ueff_starting to default :',&
        &               edelta_for_Ueff_starting

                        sw_eval_Ueff_using_edelta = ON
                     endif
                  endif
                endif
                iret = f_selectParentBlock()
             endif
          endif
! ====================================================================== 5.0
          iret = f_selectParentBlock()
       end if

! ================================= added by K. Tagami ================ 11.0
       if ( PAW_switch == ON ) then     ! default
          eval_core_level_splitting = ByPawPot
       else
!          eval_core_level_splitting = ReadFromPP
       endif

       if ( f_selectBlock( tag_spinorbit ) == 0 ) then

          iret = f_getStringValue( tag_spinorbit_mode,rstr,LOWER )
          if ( rstr == tag_spinorbit_by_projector ) then
             SpinOrbit_Mode = ByProjector
             write(nfout,*) &
                  & '!** Spin-orbit interaction strength is set using the projectors'
          else if ( rstr == tag_spinorbit_neglected ) then
             SpinOrbit_Mode = Neglected
             write(nfout,*) '!** Spin-orbit interaction is neglected '
          else if ( rstr == tag_spinorbit_by_pawpot ) then
             SpinOrbit_Mode = ByPawPot
             write(nfout,*) &
                  &   '!** Spin-orbit interaction is calculated from PAW AE potential'
          else if ( rstr == tag_spinorbit_by_zeff ) then
             SpinOrbit_Mode = ZeffApprox
             write(nfout,*) &
                  &   '!** Spin-orbit interaction is calculated with Zeff approxmation'
          else if ( rstr == tag_spinorbit_read_from_pp ) then
             SpinOrbit_Mode = ReadFromPP
             write(nfout,*) &
                  &   '!** Spin-orbit interaction is read from PP'
          else if ( rstr == tag_spinorbit_builtin ) then
             SpinOrbit_Mode = BuiltIn
             write(nfout,*) &
                  &   '!** Spin-orbit interaction is read from the pseudo potential'

          endif

          if ( SpinOrbit_Mode == ByProjector ) then
             call m_CtrlP_set_spinorbit_proj()
          end if
! === KT_add === 2014/08/11
          if ( SpinOrbit_Mode == ByPawPot ) then
             if( f_getIntValue( tag_sw_use_rphi_Hsoc_rphi, iret) == 0) then
                sw_use_rphi_Hsoc_rphi = iret
                write(nfout,*) '!** sw_use_rphi_Hsoc_rphi is ', iret
             endif
             if( f_getIntValue( tag_sw_use_ival_for_paw_ps_soc, iret) == 0) then
                sw_use_ival_for_paw_ps_soc = iret
                write(nfout,*) '!** sw_use_ival_for_paw_ps_soc is ', iret
             endif
          endif
! ============== 2014/08/11

! ==== KT_add === 2015/01/23
          if ( SpinOrbit_mode == ByPawPot ) eval_core_level_splitting = ByPawPot
          if ( SpinOrbit_mode == ReadFromPP ) eval_core_level_splitting = ReadFromPP
! =============== 2015/01/23

          iret = f_selectParentBlock()
       end if
! =========================================================================== 11.0

       if ( f_selectBlock( tag_core_electrons ) == 0 ) then
          if ( f_getIntValue( tag_sw_opencore, iret ) == 0 ) then
             sw_opencore = iret
             write(nfout,*) '!** sw_opencore is ', sw_opencore
          endif
          if ( sw_opencore == ON ) then
             if ( f_getIntValue( tag_sw_xc_opencore_ae_only, iret ) == 0 ) then
                sw_xc_opencore_ae_only = iret
                write(nfout,*) '!** sw_xc_opencore_ae_only is ', sw_xc_opencore_ae_only
             endif
             if ( f_getIntValue( tag_sw_fix_core_spin_pol, iret ) == 0 ) then
                sw_fix_core_spin_pol = iret
                write(nfout,*) '!** sw_fix_core_spin_pol is ', sw_fix_core_spin_pol
             endif
             iret = f_getStringValue( tag_spin_orientation, rstr, LOWER )
             if ( rstr == tag_anti_parallel ) then
                write(nfout,*) '!** (core) spin_orientation is anti_parallel'
                core_spin_pol_factor = -1.0d0
             else if ( rstr == tag_parallel ) then
                write(nfout,*) '!** (core) spin_orientation is parallel'
                core_spin_pol_factor = 1.0d0
             endif
          endif
          iret = f_selectParentBlock()
       endif

! === EXPERIMENTAL ===
       if ( f_selectBlock( tag_local_potential ) == 0 ) then
          if ( f_getIntValue( tag_long_range_pot_type, iret ) == 0 ) then
             long_range_pot_type = iret
             write(nfout,*) '!** long_range_pot_type is ', long_range_pot_type
          endif
          if ( f_getIntValue( tag_sw_add_vlocal_gzero, iret ) == 0 ) then
             sw_add_vlocal_gzero = iret
             write(nfout,*) '!** sw_add_vlocal_gzero is ', &
                  &   sw_add_vlocal_gzero
          endif
          iret = f_selectParentBlock()
       endif

! ======== KT_add ==== 2014/08/26
       if ( f_selectBlock( tag_orbital_decomposition ) == 0 ) then
          if ( f_getIntValue( tag_orbital_decomp_mode, iret ) == 0 ) then
             orbital_decomp_mode = iret
          endif
          if ( f_getIntValue( tag_sw_calc_orbital_moment, iret ) == 0 ) then
             sw_calc_orbital_moment = iret
          endif
          if ( f_getIntValue( tag_sw_use_contracted_psir, iret ) == 0 ) then
             sw_use_contracted_psir = iret
          endif

          if ( sw_calc_orbital_moment == ON ) then
             if ( orbital_decomp_mode == 0 ) orbital_decomp_mode = 1
          endif
          if ( orbital_decomp_mode == 1 ) sw_use_contracted_psir = ON

          write(nfout,*) '!** orbital_decomp_mode is ', orbital_decomp_mode
          write(nfout,*) '!** sw_calc_orbital_moment is ', sw_calc_orbital_moment
          write(nfout,*) '!** sw_use_contracted_psir is ', sw_use_contracted_psir

          iret = f_selectParentBlock()
       endif
! ===================== 2014.08/26

! ==== KT_add === 2014/09/19
       if ( f_selectBlock( tag_metagga ) == 0 ) then
          if ( f_getIntValue( tag_ekin_density_type, iret ) == 0 ) then
             if ( iret >= 0 .and. iret <= 3 ) then
                ekin_density_type = iret
                write(nfout,*) '!** ekin_density_type = ', iret
             endif
             if ( ekin_density_type > 0 ) then
                use_modeled_ekin_density = .true.
                write(nfout,*) '!** use_modeled_ekin_density = ', &
                     &          use_modeled_ekin_density
             endif
          endif
          if ( f_getIntValue( tag_sw_rspace_ekin_density, iret ) == 0 ) then
             sw_rspace_ekin_density = iret
             write(nfout,*) '!** sw_rspace_ekin_density = ', iret
          endif
          if ( f_getIntValue( tag_sw_calc_ekindens_hardpart, iret ) == 0 ) then
             sw_calc_ekin_density_hardpart = iret
             write(nfout,*) '!** sw_calc_ekin_density_hardpart = ', iret
          endif
          if ( f_getIntValue( tag_sw_add_ekin_hard_on_Gspace, iret ) == 0 ) then
             sw_add_ekin_hardpart_on_Gspace = iret
             write(nfout,*) '!** sw_add_ekin_hardpart_on_Gspace = ', iret
          endif
! --- tb09
          if( f_getRealValue(tag_val_c_tb09,dret,"") == 0 &
               & .or. f_getRealValue(tag_val_c_tb09,dret,"") == 0 ) then
             val_c_tb09 = dret
             write(nfout,*) '!** val c of tb09 is set to ',val_c_tb09
             sw_fix_val_c_tb09 = ON
          endif
! ---
          iret = f_selectParentBlock()
       endif
! =============== 2014/09/19

       if( f_selectBlock( tag_scf_convergence) == 0) then
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag_scf_convergence is found --")')
          tf = f_getRealValue( tag_delta_total_energy, dret, "hartree") == 0
          if(.not.tf) tf = f_getRealValue( tag_edelta, dret, "hartree") == 0
          if(tf) edelta = dret

          tf = f_getRealValue( tag_delta_total_energy_initial, dret, "hartree") == 0
          if(.not.tf) tf = f_getRealValue( tag_edelta_initial, dret, "hartree") == 0
          if(tf) then
             edelta_initial = dret
             edelta_initial_is_given = .true.
          end if

          if(.not.edelta_initial_is_given) edelta_initial = edelta
          edelta_ontheway = edelta_initial

          unit_f = 'hartree/bohr'
          tf = f_getRealValue( tag_max_force, dret, unit_f) == 0
          if( f_getRealValue(tag_max_force_edelta_i,dret,unit_f) == 0) then
             max_force_edelta_i = dret
             max_force_edelta_i_is_given = .true.
          end if

          if( f_getIntValue( tag_succession, iret) == 0) then
             if(ipriinputfile>=1) write(nfout,'(" !** tag_succession is found, iret(succession) = ",i5)') iret
             mtimes_convergence_scf = iret
          else
             if(ipriinputfile>=1) write(nfout,'(" !** tag_succession is not found")')
          end if

          tf = f_getRealValue( tag_delta_eigenvalue_conduction, dret, "hartree") == 0
          if(tf) then
             delta_eigenvalue_conduction = dret
             delta_eigenvalue_cond_is_given = .true.
             delta_eigenvalue = delta_eigenvalue_conduction
          end if

          tf = f_getRealValue( tag_sub_delta_factor, dret, "") == 0
          if(tf) then
             sub_delta_factor_is_given = .true.
             sub_delta_factor = dret
             if(ipriinputfile>=1) write(nfout,'(" sub_delta_factor_is_given = .true.")')
             if(ipriinputfile>=1) write(nfout,'(" sub_delta_factor = ",f8.4)') sub_delta_factor
             if( f_getIntValue( tag_sub_succession, iret) == 0) then
                if(ipriinputfile>=1) write(nfout,'(" tag_sub_succession is found, iret(sub_succession) = ",i5)') iret
                sub_mtimes_convergence = iret
             else
                sub_mtimes_convergence = mtimes_convergence_scf*sub_delta_factor
                if(ipriinputfile>=1) write(nfout,'(" tag_sub_succession is not found")')
             end if
             if(ipriinputfile>=1) write(nfout,'(" sub_mtimes_convergence = ",i5)') sub_mtimes_convergence
          end if
          edelta_sampling = edelta_ontheway
          if(f_getRealValue(tag_delta_total_energy_sampling,dret,'hartree')==0) then
             edelta_sampling = dret
          endif
          if(f_getStringValue(tag_convergence_criteria,rstr,LOWER)==0)then
            if(trim(rstr) == trim(tag_delta_energy)) then
               convergence_criteria = DELTA_ENERGY
            else if (trim(rstr) == trim(tag_delta_moving_average)) then
               convergence_criteria = DELTA_MOVING_AVERAGE
            else if (trim(rstr) == trim(tag_slope)) then
               convergence_criteria = SLOPE
            else if (trim(rstr) == trim(tag_delta_v)) then
               convergence_criteria = DELTA_V
            else
               if(ipriinputfile>=1) write(nfout,'(" !** unknown convergence criteria : ",a)') trim(rstr)
            endif
          endif
          if(convergence_criteria == DELTA_MOVING_AVERAGE .or. convergence_criteria == SLOPE) then
             if(f_getIntValue(tag_nsamp,iret) == 0)then
                nsamp = iret
             endif
          endif
          if(printable)then
            if (convergence_criteria == DELTA_MOVING_AVERAGE) then
               write(nfout,'(" !** convergence criteria : delta_moving_average ")')
            else if (convergence_criteria == SLOPE) then
               write(nfout,'(" !** convergence criteria : slope ")')
            else if (convergence_criteria == DELTA_V) then
               write(nfout,'(" !** convergence criteria : delta_v ")')
            else
               write(nfout,'(" !** convergence criteria : delta_energy ")')
            endif
          endif
          iret = f_selectParentBlock()
       end if

       if(ekmode == ON .or. ekmode == GRID .or. &
            & icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION .or. icond==FIXED_CHARGE_AUTOMATIC) then
          !num_extra_bands = 2 ! default value
          evaluation_eko_diff = ON ! default value
          call m_CtrlP_set_def_numextrabands(neg)
          if( f_selectBlock( tag_ek_convergence) == 0) then
             if( f_getIntValue(tag_num_max_iteration, iret) == 0)  ek_max_iteration = iret
             if( f_getIntValue(tag_sw_eval_eig_diff, iret) == 0) evaluation_eko_diff = iret
             if( f_getRealValue(tag_delta_eigenvalue, dret, "hartree") == 0) &
                  & delta_eigenvalue = dret
             if( f_getIntValue(tag_succession, iret) == 0) mtimes_convergence_ek = iret
             if( f_getIntValue(tag_num_extra_bands, iret) == 0 ) then
               num_extra_bands = iret
               num_extra_bands_was_set = .true.
             endif
!!$             if( f_getIntValue(tag_wf_inheritance, iret) == 0) wf_inheritance = iret
             iret = f_selectParentBlock()
             if(ipriinputfile >= 2) write(nfout,'(" !** -- tag_ek_convergence is found --")')
          else
             if(ipriinputfile >= 2) write(nfout,'(" !** -- tag_ek_convergence is not found --")')
          end if
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !** ek_max_iteration = ",i8)') ek_max_iteration
             write(nfout,'(" !** evaluation_eko_diff = ",i3)') evaluation_eko_diff
             write(nfout,'(" !** delta_eigenvalue = ",d20.8)') delta_eigenvalue
             write(nfout,'(" !** mtimes_convergence_ek = ",i4)') mtimes_convergence_ek
             write(nfout,'(" !** num_extra_bands = ",i6)') num_extra_bands
!!$             write(nfout,'(" !** wf_inheritance = ",i4)') wf_inheritance
          end if
       else
          num_extra_bands = 0
       end if

  ! --- charge_symmetrization ---
!       if ( f_getIntValue( tag_charge_symmetrization, iret ) == 0 ) then
!          charge_symm_mode = iret
!       endif
       if ( f_getIntValue( tag_charge_symm_mode, iret ) == 0 ) then
          charge_symm_mode = iret
       endif
       if ( charge_symm_mode == 1 ) then
          write(nfout,*) "! Charge symmetrization mode is set to ", charge_symm_mode
       endif

  ! ------- Positron start
       iret = f_getStringValue(tag_positron_method,rstr,LOWER)
       if( rstr == tag_CONV ) then
          positron_method = Positron_CONV;          sw_positron = BULK
       else if( rstr == tag_GGGC ) then
          positron_method = Positron_GGGC;          sw_positron = DEFECT
       else if( rstr == tag_PSN ) then
          positron_method = Positron_PSN;           sw_positron = DEFECT
          call phase_error_with_msg(nfout,"Positron-PSN : Not implemented",__LINE__,__FILE__)
       end if
       if ( sw_positron /= OFF ) then
          write(nfout,*) "!** Positron_method is set to ", positron_method
          if ( gmax_positron < 0.01 ) call getgmax_positron()
       endif

       if(sw_positron /= OFF) then
          npeg = 1
          num_extra_pev = 0
          if( f_selectBlock( tag_positron_convergence ) == 0) then
             if( f_getIntValue(tag_num_max_iteration, iret) == 0) pev_max_iteration = iret
             if( f_getIntValue(tag_sw_eval_eig_diff, iret) == 0) evaluation_pev_diff = iret
             if( f_getRealValue(tag_delta_eigenvalue, dret, "hartree") == 0) &
                  &  delta_pev = dret
             if( f_getIntValue(tag_succession, iret) == 0) mtimes_convergence_pev = iret
             if( f_getIntValue(tag_num_extra_bands, iret) == 0) num_extra_pev = iret
             if(num_extra_pev < 0 ) num_extra_pev = 0
             if( f_getRealValue(tag_dtim, dret, "hartree") == 0) dtim_p = dret
	     tf = ( f_getIntValue(tag_sw_submat,iret) == 0)
             if(.not.tf) tf = f_getIntValue(tag_submat,iret) == 0
             if(tf) sw_submat_p = iret
             tf = ( f_getstringValue(tag_solver_for_positronWF, rstr, LOWER) == 0)
             if(.not.tf) tf = (f_getstringValue(tag_solver, rstr, LOWER)==0)
             if(tf) then
                if(trim(rstr) == trim(tag_lmMSD)) then
                   isolver_p = lmMSD
                else if(trim(rstr) == trim(tag_lmmsd2)) then
                   isolver_p = lmMSD
                else if(trim(rstr) == trim(tag_msd)) then
                   isolver_p = MSD
                else if(trim(rstr) == trim(tag_lmsd)) then
                   isolver_p = lmSD
                else if(trim(rstr) == trim(tag_lmsd2)) then
                   isolver_p = lmSD
                else if(trim(rstr) == trim(tag_sd)) then
                   isolver_p = SD
                else
                   call phase_error_with_msg(nfout,' ! solver for positronWF is not given properly <<m_CtrlP_rd_accuracy>>'&
                                            ,__LINE__,__FILE__)
                end if
             end if

             tf = (f_getIntValue(tag_sw_gga, iret) == 0)
             if(.not.tf) tf = (f_getIntValue(tag_sw_gga_p,iret) == 0)
             if(tf) sw_gga_p = iret

             if(f_getRealValue(tag_epsilon_ele,dret,"") == 0) then
                epsilon_ele = dret
                sw_epsilon_ele = ON
             end if

             if(f_selectBlock( tag_positron_file)==0) then
                if( f_getIntValue( tag_sw_positron_file, iret) == 0) sw_positron_file = iret
                if( sw_positron_file == ON) then
! --> revised by T. Yamasaki, 28 July 2008
                   if( f_getStringValue( tag_filetype, rstr,LOWER) == 0) call set_charge_filetype(rstr,positron_filetype) !-> positron_filetype
! <--
                   if( f_getStringValue( tag_title,rstr,NOCONV) == 0) then
                      iret = min(len_trim(rstr),LEN_TITLE)
                      if(ipriinputfile>=1) write(nfout,'(" positron_title = ",a80)') rstr
                      if(ipriinputfile>=1) write(nfout,'(" iret = ",i8)') iret
!!$                      write(positron_title(1),'(a80)')
                      positron_title(1) = ''
                      positron_title(1)(1:iret) = rstr(1:iret)
!!$                      positron_title(1)(1:iret) = rstr(1:iret)
                      positron_title(2) = positron_title(1)
                      positron_title(3) = positron_title(1)
                      positron_title(5) = positron_title(1)
                      if(sw_gga_p == ON) positron_title(4) = positron_title(1)
                      if(ipriinputfile>=1) then
                         write(nfout,'(" positron_title(1) = ",a80)') positron_title(1)
                         write(nfout,'(" positron_title(2) = ",a80)') positron_title(2)
                         write(nfout,'(" positron_title(3) = ",a80)') positron_title(3)
                         write(nfout,'(" positron_title(5) = ",a80)') positron_title(5)
                         if(sw_gga_p == ON) write(nfout,'(" positron_title(4) = ",a80)') positron_title(4)
                      end if
                   end if
                   if( f_getStringValue( tag_title_positr,rstr,NOCONV) == 0) &
                        & positron_title(1) = rstr(1:min(len_trim(rstr),LEN_TITLE))
                   if( f_getStringValue( tag_title_velect,rstr,NOCONV) == 0) &
                        & positron_title(2) = rstr(1:min(len_trim(rstr),LEN_TITLE))
                   if( f_getStringValue( tag_title_eppair,rstr,NOCONV) == 0) &
                        & positron_title(3) = rstr(1:min(len_trim(rstr),LEN_TITLE))
                   if( sw_gga_p == ON .and. f_getStringValue( tag_title_velgrd,rstr,NOCONV) == 0) &
                        & positron_title(4) = rstr(1:min(len_trim(rstr),LEN_TITLE))
                end if
                iret = f_selectParentBlock()
             end if

             if(ipriinputfile >= 2 .and. printable) then
                write(nfout,'(" !** -- tag_positron_convergence --")')
                write(nfout,'(" !** pev_max_iteration      = ",i8)') pev_max_iteration
                write(nfout,'(" !** evaluation_pev_diff    = ",i6)') evaluation_pev_diff
                write(nfout,'(" !** delta_pev              = ",d20.8)') delta_pev
                write(nfout,'(" !** mtimes_convergence_p6v = ",i4)') mtimes_convergence_pev
                write(nfout,'(" !** num_extra_pev          = ",i6)') num_extra_pev
                write(nfout,'(" !** dtim_p                 = ",d20.8)') dtim_p
                if(sw_positron == ON) then
                   write(nfout,'(" !** sw_submat              = ",i6)') sw_submat_p
                   write(nfout,'(" !** isolver_p              = ",i6)') isolver_p
                   write(nfout,'(" !** sw_gga_p               = ",i6)') sw_gga_p
                   write(nfout,'(" !** sw_epsilon_ele         = ",i6)') sw_epsilon_ele
                   if(sw_epsilon_ele == ON) &
                        & write(nfout,'(" !** epsilon_ele            = ",d20.8)') epsilon_ele
                   write(nfout,'(" !** sw_positron_file       = ",i6)') sw_positron_file
                   if(sw_positron_file == ON) then
                      write(nfout,'(" !** positron_filetype      = ",i6," :0=DENSITY_ONLY,1=CUBE,2=VTK")') positron_filetype
                      write(nfout,'(" !** positron_title_positr  = ",a80)') positron_title(1)
                      write(nfout,'(" !** positron_title_velect  = ",a80)') positron_title(2)
                      write(nfout,'(" !** positron_title_eppair  = ",a80)') positron_title(3)
                      write(nfout,'(" !** positron_title_eppair2 = ",a80)') positron_title(5)
                      if(sw_gga_p == ON) &
                           & write(nfout,'(" !** positron_title_velgrd  = ",/," !** ",a80)') positron_title(4)
                   end if
                end if
             end if
             iret = f_selectParentBlock()

          end if
          npeg = npeg + num_extra_pev
       end if
  ! ------- Positron end

  ! ------- Dipole correction start
       if(sw_dipole_correction /= OFF) then
          if( f_selectBlock( tag_dipole_correction ) == 0) then
             if(f_getIntValue( tag_direction, iret ) == 0) idir_dip = iret
             if(f_getIntValue( tag_division, iret ) == 0) ndiv_dip = iret
             if( f_selectBlock( tag_vacuum ) == 0) then
                if(f_getRealValue( tag_rx, dret, "") == 0) rvac(1) = dret
                if(f_getRealValue( tag_ry, dret, "") == 0) rvac(2) = dret
                if(f_getRealValue( tag_rz, dret, "") == 0) rvac(3) = dret
                iret = f_selectParentBlock()
             end if
             if( f_selectBlock( tag_electric_field ) == 0) then
                if(f_getRealValue( tag_ex, dret, "hartree/bohr") == 0) elec_field(1) = dret
                if(f_getRealValue( tag_ey, dret, "hartree/bohr") == 0) elec_field(2) = dret
                if(f_getRealValue( tag_ez, dret, "hartree/bohr") == 0) elec_field(3) = dret
                iret = f_selectParentBlock()
             end if
             if(f_getRealValue( tag_amix, dret, "") == 0) amix_dip = dret
             iret = f_selectParentBlock()
          end if
       end if
  ! ------- Dipole correction end

  ! ---- Screening correction --
      if(sw_screening_correction /= OFF) then
         if( f_selectBlock( tag_screening_correction ) == 0) then
            if(f_getRealValue( tag_screening_alpha, dret, '') == 0) screening_alpha = dret
write(nfout,'(" !** sw_screening_correction = ",i5," alpha = ",f8.3)') sw_screening_correction, screening_alpha
         end if
         iret = f_selectParentBlock()
      end if
  ! --------------------

  ! ------- FEF start
       if(sw_fef /= OFF) then
          if( f_selectBlock( tag_fef ) == 0) then
             if( f_getIntValue( tag_sw_check_polar, iret) == 0) sw_check_polar = iret
             if( f_selectBlock( tag_electric_field ) == 0) then
                if(f_getRealValue( tag_ex, dret, "hartree/bohr") == 0) elec_field(1) = dret
                if(f_getRealValue( tag_ey, dret, "hartree/bohr") == 0) elec_field(2) = dret
                if(f_getRealValue( tag_ez, dret, "hartree/bohr") == 0) elec_field(3) = dret
                iret = f_selectParentBlock()
             end if
             iret = f_selectParentBlock()
          end if
       end if
  ! ------- FEF end

  ! ------- vdW start
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!      if(sw_pair_vdw /= OFF) then
!         if(f_getRealValue( tag_rcut_vdw, dret, "bohr") == 0) rcut_vdw = dret
       if(sw_vdw_correction /= OFF) then
          iret = f_getStringValue(tag_vdw_method,rstr,LOWER)
          if( rstr == tag_williams) then
             vdw_method = VDW_WILLIAMS
          else if( rstr == tag_grimme) then
             vdw_method = VDW_GRIMME
          else if( rstr == tag_dft_d2) then
             vdw_method = VDW_GRIMME
          else if( rstr == tag_dft_d3 .or. rstr == tag_dftd3) then
             vdw_method = VDW_DFTD3
          end if

          ! set default vdw parameters
          select case(vdw_method)
          case(VDW_WILLIAMS)
            vdw_radius = 20.0
            vdw_scaling_factor = 0.8095
            vdw_scaling_factor_r = 0.80
            vdw_damping_factor = 3.0
          case(VDW_GRIMME)
            vdw_radius = 30.0
            vdw_scaling_factor = 0.75
            vdw_scaling_factor_r = 1.0
            vdw_damping_factor = 20.0
          end select

          if(f_getRealValue( tag_vdw_radius, dret, "bohr") == 0) vdw_radius = dret
          if(vdw_method == VDW_GRIMME) vdw_radius = vdw_radius / BOHR
          if(f_getRealValue( tag_vdw_scaling_factor, dret, '') == 0) vdw_scaling_factor = dret
          if(f_getRealValue( tag_vdw_scaling_factor_r, dret, '') == 0) vdw_scaling_factor_r = dret
          if(f_getRealValue( tag_vdw_damping_factor, dret, '') == 0) vdw_damping_factor = dret
! ==============================================================================
       end if

       if ( vdw_method == VDW_DFTD3 ) then
          if ( f_selectBlock( tag_dftd3 ) == 0 ) then
             iret = f_getStringValue( tag_dftd3_damping_function,rstr,LOWER )
             if ( rstr == tag_dftd3_damping_zero ) then
                dftd3_damping_function = 0
             else if ( rstr == tag_dftd3_damping_bj ) then
                dftd3_damping_function = 1
             endif
             iret = f_selectParentBlock()
          endif
          write(nfout,*) "!** dftd3_damping function = ", dftd3_damping_function
       endif
  ! ------- vdW end

       if( f_getStringValue(tag_initial_wavefunctions,rstr,LOWER) == 0) call set_intzaj(rstr)
       if(intzaj == by_matrix_diagon) then
          if( f_selectBlock( tag_matrix_diagon) == 0) then
             call getgmaxs()
             if( f_getIntValue( tag_msz_initial_matdiagon, iret) == 0) n_matrix_size = iret
             if( f_getRealValue( tag_eps_solve_Hx_eq_ex, dret, '') == 0) eps_solve_Hx_eq_ex = dret
             if( f_getIntValue( tag_sw_scalapack, iret) == 0) sw_scalapack_md = iret
             if(ipriinputfile >= 1 .and. printable .and. n_matrix_size >= 1) &
                  & write(nfout,'(" !** n_matrix_size = ",i8)') n_matrix_size
             iret = f_selectParentBlock()
          else
             if(ekmode == OFF) then
                gmaxs_given = gmax
             else
                gmaxs_given = gmax*0.5d0
             end if
          end if
       end if
!
       if ( intzaj == by_pseudo_atomic_orbitals ) then
!          if ( noise_mode == 0 ) noise_amplitude = 0.05d0

          if ( f_selectBlock( tag_pao_setting ) == 0 ) then
             if ( f_getIntValue( tag_sw_exclude_orb_tau_neq_1, iret ) == 0 ) then
                sw_exclude_orb_tau_neq_1 = iret         ! non-PAW pot only
                write(nfout,*) '!** sw_exclude_orb_tau_neq_1 =', iret
             endif
             if ( f_getIntValue( tag_sw_exclude_orb_unoccupied, iret ) == 0 ) then
                sw_exclude_orb_unoccupied = iret        ! PAW pot only
                write(nfout,*) '!** sw_exclude_orb_unoccupied =', iret
             endif

             if( f_getIntValue( tag_noise_mode, iret ) == 0 ) then
                noise_mode = iret
                write(nfout,*) "!** noise mode= ", iret
             endif
             if( f_getRealValue( tag_noise_amplitude, dret, '') == 0 ) then
                noise_amplitude = dret
                write(nfout,*) "!** noise amplitude = ", dret
             endif

             iret = f_selectParentBlock()
          endif
       endif
!
       if(intzaj==by_pseudo_atomic_orbitals.and.(icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION))then
          if(printable) write(nfout,&
          & '(" !** WARNING : initial_wavefunctions = atomic_orbitals unsupported for fixed-charge calculations")')
          intzaj = by_random_numbers
       endif

       if( f_getStringValue(tag_initial_charge_density, rstr,LOWER) == 0) call set_initial_chg(rstr) ! -> initial_chg
! --> T. Yamasaki, 28 July 2008
       if(initial_chg == FILE) then
          if( f_selectBlock( tag_initial_charge_density_file) == 0) then
             if( f_getIntValue(tag_sw_charge_rspace,iret) == 0) sw_initial_charge_rspace = iret
             if( f_getStringValue( tag_filetype, rstr, LOWER) == 0) then
                call set_charge_filetype(rstr, initial_charge_filetype)
                if(initial_charge_filetype == VTK) then
                   initial_charge_filetype = CUBE
                   if(ipri>= 1) then
                      write(nfout,'(" !** (initial_charge_filetype = VTK) is invalid")')
                      write(nfout,'(" !** initial_charge_filetype is set CUBE")')
                   end if
                end if
             end if
             iret = f_selectParentBlock()
          end if
       end if
       if(ipri >= 1) then
          write(nfout,'(" !** initial_chg = ",i3," 1:Gauss_distrib_func, 2:Gauss narrow, 3:from_PP_file, 4: file")')
       end if
! <--

       if( f_getStringValue(tag_initial_kinetic_density, rstr,LOWER) == 0) then
          call strncmp0(trim(rstr), tag_given_by_a_file, tf)
          if(.not.tf) call strncmp0(trim(rstr), tag_file, tf)
          if(tf) then
             initial_ekin_dens = FILE
          end if
       endif

! =============================== added by K. Tagami ======================= 11.0
       if( f_selectBlock( tag_ports ) == 0) then
          if( f_getIntValue( tag_import_collinear_spindens, iret) == 0 ) then
             import_collinear_spindensity = iret
             write(nfout,*) "!** import_collinear_spindensity is ", iret
          end if
          if( f_getIntValue( tag_import_collinear_wfns, iret) == 0 ) then
             import_collinear_wavefunctions = iret
             write(nfout,*) "!** import_collinear_wavefunctions is ", iret
          end if
          if( f_getIntValue( tag_previous_nspin, iret) == 0 ) then
             previous_nspin_collinear = iret
             if ( iret > 2 ) previous_nspin_collinear = 2
             if ( iret < 1 ) previous_nspin_collinear = 2
             write(nfout,*) "!** previous_nspin_collinear is ", previous_nspin_collinear
          endif
          previous_nband_collinear = neg /2
          if( f_getIntValue( tag_previous_nband, iret) == 0 ) then
             previous_nband_collinear = iret
             if ( iret < 1 ) previous_nband_collinear = neg /2
          endif
          write(nfout,*) "!** previous_nband_collinear is ", previous_nband_collinear
          iret = f_selectParentBlock()
       endif
! =========================================================================== 11.0

       if( f_getIntValue(tag_sw_add_qex_to_initial_charge,iret)==0) sw_add_qex_to_initial_charge = iret
! -----------
       if( f_getIntValue( tag_read_charge_hardpart, iret) == 0 ) then
          read_charge_hardpart = iret
       endif
! -----------

       if( f_selectBlock( tag_precalculation) == 0) then
          if( f_getIntValue( tag_nel_Ylm, iret) == 0) nel_Ylm = iret
          if( nel_Ylm < 0) nel_Ylm = 0
          iret = f_selectParentBlock()
       end if

       if(ipriinputfile >= 2 .and. printable) then
          write(nfout,'(" !** num_bands = ",i6)') neg
          if(ekmode == ON .or. ekmode == OFF) write(nfout,'(" !** num_extra_bands = ",i6)') num_extra_bands
          write(nfout,'(" !** smearing_method = ",i6)') way_of_smearing
          write(nfout,'(" !** width = ", f8.4)') width
          write(nfout,'(" !** idimtetra = ", i6)') idimtetra
          write(nfout,'(" !** xctype = ",a7)') xctype
          write(nfout,'(" !** PAW_switch = ",a3)') on_or_off(PAW_switch)
          if(PAW_switch == ON) then
             if(paw_density_gradient%interpolation_method==LAGRANGE) then
                write(nfout,'(" !** - interpolation_method = LAGRANGE = ",i2)') paw_density_gradient%interpolation_method
             else if (paw_density_gradient%interpolation_method==NEWTON) then
                write(nfout,'(" !** - interpolation_method = NEWTON   = ",i2)') paw_density_gradient%interpolation_method
             else
                write(nfout,'(" !** - interpolation_method = ",i2)') paw_density_gradient%interpolation_method
             end if
             write(nfout,'(" !** - order                = ",i2)') paw_density_gradient%order
          end if
          write(nfout,'(" !** ggacmp_parallel = ",i6)') ggacmp_parallel
          write(nfout,'(" !** edelta = ",d20.8)') edelta
          if(edelta_initial_is_given)  write(nfout,'(" !** edelta_initial = ",d20.8)') edelta_initial
          if(max_force_edelta_i_is_given) write(nfout,'(" !** max_force_edelta_i = ",d20.8)') max_force_edelta_i
          write(nfout,'(" !** mtimes_convergence_scf = ",i6)') mtimes_convergence_scf
          write(nfout,'(" !** delta_eigenvalue = ",d20.8)') delta_eigenvalue
          write(nfout,'(" !** mtimes_convergence_ek  = ",i6)') mtimes_convergence_ek
          write(nfout,'(" !** initial_chg = ",i6)') initial_chg
! --> T. Yamasaki, 28 July 2008
          if(initial_chg == FILE) then
             write(nfout,'(" !** sw_initial_charge_rspace = ",i6)') sw_initial_charge_rspace
             write(nfout,'(" !** initial_charge_filetype = ",i6)') initial_charge_filetype
          end if
! <--
          write(nfout,'(" !** initial_wavefucntions = ",i6)') intzaj
          write(nfout,'(" !** nel_Ylm = ",i6)') nel_Ylm
          if(sw_hubbard == ON) then
             write(nfout,'(" !** initial_occmat = ",i3, " : 0=OFF, 2=FILE")') initial_occmat
             write(nfout,'(" !** sw_hubbard = ",i2)') sw_hubbard
             write(nfout,'(" !** critical_ehub = ",e12.5)') critical_ehub
             write(nfout,'(" !** delta_ehub = ",e14.5)') delta_ehub
             write(nfout,'(" !** amix = ",e12.5)') alpha_hubbard
             if(sw_constraint == ON) then
                write(nfout,'(" !** sw_constraint = ",i2)') sw_constraint
                write(nfout,'(" !** site = ",i2)') const_site
                write(nfout,'(" !** alpha = ",f10.5)') const_alpha
             end if
          end if
          if(num_projectors > 0) then
             write(nfout,'(" !** projector_type = ",i2)') projector_type
             call m_CtrlP_wd_proj_attr(nfout)
          end if
       end if
       if(ekmode==ON .or. ekmode==GRID .or. icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION ) then
          neg = neg + num_extra_bands
          neg_is_enlarged = .true.
       end if
       if(ipriinputfile >= 1 .and. printable ) then
!!$          if(ekmode == ON .or. ekmode == GRID) &
          if(neg_is_enlarged) write(nfout,'(" !** neg (=num_bands+num_extra_bands) = ",i6)') neg
       end if
       neg_previous = neg

  ! ------- Positron start
       if(ipriinputfile >= 1 .and. sw_positron /= OFF) then
          if(printable) then
             write(nfout,'(" !** --- parameters for positron --")')
             write(nfout,'(" !** pev_max_iteration = ",i6)') pev_max_iteration
             write(nfout,'(" !** num_extra_bands = ",i6)') num_extra_pev
             write(nfout,'(" !** delta_pev = ",d20.8)') delta_pev
             write(nfout,'(" !** mtimes_convergence_pev = ",i6)') mtimes_convergence_pev
             write(nfout,'(" !** pev_max_iteration      = ",i6)') pev_max_iteration
             write(nfout,'(" !** evaluation_pev_diff = ",i6)') evaluation_pev_diff
             write(nfout,'(" !** dtim_p                 = ",d20.8)') dtim_p
             write(nfout,'(" !** sw_submat              = ",i6)') sw_submat_p
             if(isolver_p == lmMSD) then
                write(nfout,'(" !** isolver_p              =  lmMSD")')
             else if(isolver_p == MSD) then
                write(nfout,'(" !** isolver_p              =  MSD")')
             else if(isolver_p == lmSD) then
                write(nfout,'(" !** isolver_p              =  lmSD")')
             else if(isolver_p == SD) then
                write(nfout,'(" !** isolver_p              =  SD")')
             else
                write(nfout,'(" !** isolver_p              =  ",i6)') isolver_p
             end if
             write(nfout,'(" !** sw_gga_p               = ",i6)') sw_gga_p
             write(nfout,'(" !** sw_epsilon_ele         = ",i6)') sw_epsilon_ele
             if(sw_epsilon_ele == ON) &
                  & write(nfout,'(" !** epsilon_ele            = ",d20.8)') epsilon_ele
          end if
       end if
  ! ------- Positron end

  ! ------- Dipole correction start
       if(ipriinputfile >= 1 .and. sw_dipole_correction /= OFF) then
          if(printable) then
             write(nfout,'(" !** --- parameters for dipole correction --")')
             write(nfout,'(" !** direction = ",i6)') idir_dip
             write(nfout,'(" !** division = ",i6)') ndiv_dip
             write(nfout,'(" !** vacuum position (internal)    = ",3f10.5)') rvac(1:3)
             write(nfout,'(" !** electric field (Hartree/Bohr) = ",3f10.5)') elec_field(1:3)
          end if
       end if
  ! ------- Dipole correction end

  ! ------- FEF start
       if(ipriinputfile >= 1 .and. sw_fef /= OFF) then
          if(printable) then
             write(nfout,'(" !** --- parameters for FEF --")')
             write(nfout,'(" !** electric field (Hartree/Bohr) = ",3f10.5)') elec_field(1:3)
          end if
       end if
  ! ------- FEF end

  ! ------- vdW start
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!      if(ipriinputfile >= 1 .and. sw_pair_vdw /= OFF) then
       if(ipriinputfile >= 1 .and. sw_vdw_correction /= OFF) then
! ==============================================================================
          if(printable) then
             write(nfout,'(" !** --- parameters for vdW pair interaction --")')
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!            write(nfout,'(" !** rcut_vdw (Bohr) = ",f10.5)') rcut_vdw
             write(nfout,'(" !** vdw_method = ",i3)') vdw_method
             if ( vdw_method == VDW_WILLIAMS .or. vdw_method == VDW_GRIMME ) then
                write(nfout,'(" !** vdw_radius (Bohr) = ",f10.5)') vdw_radius
                write(nfout,'(" !** vdw_scaling_factor = ",f10.5)') vdw_scaling_factor
                write(nfout,'(" !** vdw_scaling_factor_r = ",f10.5)') &
                     &                                             vdw_scaling_factor_r
                write(nfout,'(" !** vdw_damping_factor = ",f10.5)') vdw_damping_factor
             endif
! ==============================================================================
          end if
       end if
  ! ------- vdW end

#endif

       if( f_selectBlock( tag_force_convergence) == 0) then
          unit_f = 'hartree/bohr'
          tf = f_getRealValue( tag_max_force, dret, unit_f) == 0
          if(.not.tf) tf = f_getRealValue( tag_delta_force, dret, unit_f) == 0
          if(tf) forccr = dret
!!$          if( f_getRealValue( tag_delta_force, dret,unit_f) == 0) forccr = dret
          tf = f_getRealValue( tag_tolerable_norm_error,dret,'') == 0
          if(.not.tf) tf = f_getRealValue( tag_tolerable_force_error,dret,'') == 0
          if(.not.tf) tf = f_getRealValue( tag_tolerable_force_norm_error,dret,'') == 0
          if(.not.tf) tf = f_getRealValue( tag_tolerable_error,dret,'') == 0
          if(tf) then
             f_tolerable_norm_error = dret
             force_error_check_mode = ON
          end if
          tf = f_getRealValue( tag_tolerable_angle_error,dret,'radian') == 0
          if(tf) then
             f_tolerable_angle_error = dret
             force_error_check_mode = ON
          end if
          tf = f_getRealValue( tag_tolerable_hyper_angle_error,dret,'radian') == 0
          if(tf) then
             f_tolerable_hyper_angle_error = dret
             force_error_check_mode = ON
          end if

          tf = f_getRealValue( tag_force_error_check_rangeL, dret,unit_f) == 0
          if(.not.tf) tf = f_getRealValue( tag_error_check_rangeL, dret, unit_f) == 0
          if(tf) then
             force_error_check_rangeL = dret
          else
             force_error_check_rangeL = forccr
          end if

          if(f_getRealValue(tag_max_translational_force,dret,unit_f)==0) then
             max_force_trans = dret
          endif

          iret = f_selectParentBlock()
       end if
       if(f_selectBlock(tag_torque_convergence)==0) then
          if(f_getRealValue(tag_max_torque, dret, 'hartree')==0) then
             max_torque = dret
          endif
          iret = f_selectParentBlock()
       endif
       if(f_selectBlock(tag_stress_convergence)==0)then
          if (f_getRealValue(tag_max_stress,dret,'hartree/bohr3')==0)then
             max_stress = dret
          endif
          iret = f_selectParentBlock()
       endif
       if(ipriinputfile >= 1 .and. printable) then
          write(nfout,'(" !** --- parameters for force_convergence --")')
          write(nfout,'(" !** forccr (=max_force) = ",d20.8)') forccr
          if(force_error_check_mode == ON) then
             write(nfout,'(" !** force_error_check_mode = ON")')
             write(nfout,'(" !** tolerable_norm_error        = ",f10.6)') f_tolerable_norm_error
             write(nfout,'(" !** tolerable_angle_error       = ",f10.6," radian")') f_tolerable_angle_error
             write(nfout,'(" !** tolerable_hyper_angle_error = ",f10.6," radian")') f_tolerable_hyper_angle_error
             write(nfout,'(" !** force_error_check_rangeL    = ",d20.8)') force_error_check_rangeL
          else
             write(nfout,'(" !** force_error_check_mode = OFF")')
          end if
       end if

       if(.not.max_force_edelta_i_is_given) then
          max_force_edelta_i = forccr * MAX_FORCE_EDELTA_I_FACTOR
          if(ipriinputfile >= 1 .and. printable) &
               & write(nfout,'(" !** max_force_edelta_i = ",d20.8)') max_force_edelta_i
       end if

! =================================== added by K. Tagami ============= 5.0
       if( f_getIntValue( tag_eval_energy_before_charge, iret) == 0 ) then
          eval_energy_before_charge = iret
       endif
       if(printable) write(nfout,*) '!** evaluating energy before charge construction is ', &
        &        eval_energy_before_charge
! =================================================================== 5.0

#ifdef ENABLE_ESM_PACK
! ====================== ESM ==========================
       if (f_selectBlock(tag_esm)==0)then
          if(f_getIntValue(tag_sw_esm,iret) == 0 )then
             sw_esm = iret
          endif
          if(sw_esm==ON .and. printable)then
             write(nfout,'(a)') ' !** ESM enabled'
          endif
          if(f_getRealValue(tag_z1,dret,'bohr')==0)then
             esm_z1 = dret
             esm_z1_defined = .true.
          endif
          if(f_getRealValue(tag_w,dret,'bohr')==0)then
             esm_w = dret
          endif
          if(f_getRealValue(tag_fix_ef,dret,'')==0)then
             esm_fix_ef = dret
          endif
          if(f_getRealValue(tag_add_elec,dret,'')==0)then
             esm_add_elec = dret
             esm_qbac = -esm_add_elec
          endif
          if(f_getRealValue(tag_charge,dret,'')==0)then
             esm_qbac = dret
          endif
          if(f_getRealValue(tag_z_wall,dret,'bohr')==0)then
             esm_z_wall = dret
             esm_izwall = 1
          endif
          if(f_getRealValue(tag_bar_height,dret,'hartree')==0)then
             esm_bar_height = dret
          endif
          if(f_getRealValue(tag_bar_width,dret,'bohr')==0)then
             esm_bar_width = dret
          endif
          if(f_getRealValue(tag_electric_field,dret,"hartree/bohr")==0)then
             esm_e_field = dret
          endif
          if(f_getIntValue(tag_nosmooth,iret)==0)then
             if(iret==ON.and.esm_izwall==1)then
                esm_izwall = -1
             endif
          endif
          if(f_getIntValue(tag_external_potential,iret)==0)then
             if(iret==ON) esm_iexpot = ON
          endif
          if(f_getIntValue(tag_gps,iret)==0)then
             esm_gps = iret
          endif
          if(f_getIntValue(tag_gpe,iret)==0)then
             esm_gpe = iret
          endif
          if(f_getStringValue(tag_bc,rstr,LOWER)==0)then
             if(rstr==tag_bare)then
                esm_bc = BARE
             else if (rstr==tag_pe1)then
                esm_bc = PE1
             else if (rstr==tag_pe2)then
                esm_bc = PE2
             else
                write(nfout,'(a)') 'unsupported boundary condition : '//trim(rstr)
                call phase_error_with_msg(nfout,'unsupported boundary condition : '//trim(rstr),__LINE__,__FILE__)
             endif
          endif
          iret = f_selectParentBlock()
       endif
! ====================== ESM ==========================
#endif

! ====================== FCP ==========================
       if (f_selectBlock(tag_fcp)==0)then
          if(f_getIntValue(tag_sw_fcp,iret) == 0 )then
             sw_fcp = iret
          endif
          if(sw_fcp ==ON .and. printable)then
             write(nfout,'(a)') ' !** FCP enabled'
          endif
          if(f_getRealValue(tag_fcp_mu,dret,'hartree')==0)then
             fcp_mu = dret
          endif
          if(f_getRealValue(tag_fcp_mass,dret,'')==0)then
             fcp_mass = dret
          endif
          if(f_getRealValue(tag_fcp_temperature,dret,'k')==0)then
             fcp_temperature = dret
          endif
          if(f_getRealValue(tag_fcp_relax_step,dret,'')==0)then
             fcp_relax_step = dret
          endif
          if(f_getRealValue(tag_fcp_relax_crit,dret,'')==0)then
             fcp_relax_crit = dret
          endif
          if(f_getRealValue(tag_fcp_tot_charge_first,dret,'')==0)then
             fcp_tot_charge_first = dret
          endif
          if(f_getRealValue(tag_fcp_tot_charge_last,dret,'')==0)then
             fcp_tot_charge_last = dret
          endif
          if(f_getRealValue(tag_fcp_qmass,dret,'')==0)then
             fcp_qmass = dret
          endif
          iret = f_selectParentBlock()
       endif
! ====================== FCP ==========================

       if(f_selectBlock(tag_nonlocal_potential)==0)then
          if(f_getIntValue(tag_sw_rspace,iret)==0) sw_rspace = iret
          sw_rspace_v = sw_rspace
          if(f_getIntValue(tag_sw_rspace_v,iret)==0) sw_rspace_v = iret
          if(f_getRealValue(tag_r0_factor,dret,'')==0) r0_factor = dret
          !!if(f_getIntValue(tag_nq,iret)==0) nq = iret
          if(f_getRealValue(tag_dq,dret,'')==0) dq = dret
          if(f_getIntValue(tag_sw_save_memory,iret)==0) sw_save_memory = iret
          if(f_getStringValue(tag_projector_optimization,rstr,LOWER)==0)then
             if(rstr==tag_mask_function)then
                projector_optimization = MASK_FUNCTION
             elseif(rstr==tag_prefitting)then
                projector_optimization = PREFITTING
             elseif(rstr==tag_none)then
                projector_optimization = NO
             else
                write(nfout,'(a)') 'unsupported projector optimization method : '//trim(rstr)
                call phase_error_with_msg(nfout, &
                'unsupported projector optimization method : '//trim(rstr),__LINE__,__FILE__)
             endif
          endif
          if(projector_optimization==PREFITTING) gamma_factor = 3.0d0
          if(f_getRealValue(tag_gamma_factor,dret,'')==0) gamma_factor = dret
          iret = f_selectParentBlock()
       endif
       if(printable.and.ipriinputfile>=1) write(nfout,'(a,i3)') ' !** nonlocal potential in real space : ',sw_rspace
       if(sw_rspace==ON.and.printable.and.ipriinputfile>=1)then
          write(nfout,'(a)')        ' !** parameters for the calclulation of the nonlocal potential in real space'
          write(nfout,'(a,f6.3,a)') ' !** sphere radius = ',r0_factor,   ' x rcut'
          write(nfout,'(a,f6.3,a)') ' !** gamma         = ',gamma_factor,' x Gmax'
          write(nfout,'(a,i3)')     ' !** save memory   = ',sw_save_memory
          if(projector_optimization==MASK_FUNCTION)then
          write(nfout,'(a)')        ' !** projector optimization method : by mask function'
          else if(projector_optimization==PREFITTING)then
          write(nfout,'(a)')        ' !** projector optimization method : by prefitting procedure'
          else if(projector_optimization==NO)then
          write(nfout,'(a)')        ' !** projector optimization method : no optimization (for debugging)'
          endif
       endif

! ====================== RSB ==========================
       if(f_selectBlock(tag_rsb)==0)then
          if(f_getIntValue(tag_sw_rsb,iret)==0) sw_rsb=iret
          if(f_getIntValue(tag_sw_valence_electrons_only,iret)==0) sw_valence_electrons_only = iret
          if(f_getIntValue(tag_bisect_by,iret)==0) bisect_by = iret
          if(f_getIntValue(tag_lmax,iret)==0) lmax_rsb = iret
          if(f_getRealValue(tag_eps_rsb,dret,"")==0) eps_rsb = dret
          iret = f_selectParentBlock()
       endif

! ====================== vdW-DF ==========================
       if(f_selectBlock(tag_vdwdf)==0)then
         if(f_getStringValue(tag_mode,rstr,LOWER)==0) then
            if(rstr==tag_scf)then
               oneshot = .false.
            endif
         endif
         if(f_getIntValue(tag_ndel,iret)==0)       ndel   = iret
         if(f_getIntValue(tag_nphiD,iret)==0)      nphiD  = iret
         if(f_getIntValue(tag_nr12,iret)==0)       nr12   = iret
         if(f_getRealValue(tag_maxk,dret,'')==0)   maxk   = dret
         if(f_getRealValue(tag_r12max,dret,"bohr")==0) r12max = dret
         if(f_getRealValue(tag_dq,dret,'')==0)     dq_vdw = dret
         if(f_getRealValue(tag_lambda,dret,'')==0) lambda = dret
         if(f_getRealValue(tag_q0cut,dret,'')==0)  q0cut  = dret
         if(f_getRealValue(tag_q0min,dret,'')==0)  q0min  = dret
         if(f_getRealValue(tag_ds,dret,'')==0)     ds     = dret
         if(f_getIntValue(tag_na_gl,iret)==0)      na_gl  = iret
         if(f_getRealValue(tag_a1,dret,'')==0)     a1     = dret
         if(f_getRealValue(tag_a2,dret,'')==0)     a2     = dret
         if(f_getIntValue(tag_eval_kernel_by_interpolation,iret)==0) eval_kernel_by_interpolation = iret == ON
         if(f_getIntValue(tag_save_memory_mode,iret)==0) sw_save_memory_vdw = iret == ON
         if(f_getIntValue(tag_sw_use_WuGygi_method,iret)==0) then
            sw_use_WuGygi_method = iret
         endif

         if( f_getIntValue(tag_vdwdf_version,iret)==0 )   vdwdf_version = iret
         if( f_getStringValue( tag_exchange_pot_type, rstr,LOWER ) == 0 ) then
            exchange_pot_type = rstr(1:len_exchange_pot_type)
         endif

         write(nfout,*) '! vdwdf info '
         write(nfout,*) "!** exchange potential type = ", exchange_pot_type
         write(nfout,*) "!** vdwdf version no. = ", vdwdf_version

         iret = f_selectParentBlock()
       endif
       if(f_getIntValue(tag_sw_read_pwbs_info,iret)==0) then
         sw_read_pwbs_info = iret
       endif
       if(f_getIntValue(tag_sw_write_pwbs_info,iret)==0) then
         sw_write_pwbs_info = iret
       endif

       iret = f_selectParentBlock()
    else
       if(printable) write(nfout,'(" !** -- tag_accuracy is not found --")')
       call phase_error_with_msg(nfout, &
       ' tag_accuracy is not given in the inputfile <<CtrlP_rd_accuracy>>', &
       __LINE__,__FILE__)
    end if
#ifndef _EMPIRICAL_
  contains

!!$    subroutine set_positron_filetype(rstr)
!!$      character(len=FMAXVALLEN),intent(in) :: rstr
!!$      logical :: tf
!!$      call strncmp2(rstr, FMAXVALLEN, tag_cube, len(tag_cube),tf)
!!$      if(tf) then
!!$         positron_filetype = CUBE
!!$         goto 1001
!!$      end if
!!$      call strncmp2(rstr, FMAXVALLEN, tag_vtk, len(tag_vtk),tf)
!!$      if(tf) then
!!$         positron_filetype = VTK
!!$         goto 1001
!!$      end if
!!$      call strncmp2(rstr, FMAXVALLEN, tag_density_only, len(tag_density_only),tf)
!!$      if(tf) then
!!$         positron_filetype = DENSITY_ONLY
!!$         goto 1001
!!$      end if
!!$1001  continue
!!$    end subroutine set_positron_filetype

    subroutine set_initial_occmat(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_off, len(tag_off), tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_off, tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_none,tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_0,   tf)
      if(tf) then
         initial_occmat = OFF
         goto 1001
      end if

      call strncmp0(trim(rstr), tag_unit_matrix, tf)
      if(tf) then
         initial_occmat = UNIT_MATRIX
         goto 1001
      endif

      call strncmp0(trim(rstr), tag_spin_polarized, tf)
      if(tf) then
         initial_occmat = SPIN_POLARIZED
         goto 1001
      endif

      call strncmp0(trim(rstr), tag_initial_es, tf)
      if(tf) then
         sw_initial_es = ON
      endif

      call strncmp0(trim(rstr), tag_crystal_field_approx, tf)
      if(tf) then
         initial_occmat = CRYSTAL_FIELD_APPROX
         goto 1001
      endif

      call strncmp0(trim(rstr), tag_given_by_a_file, tf)
      if(.not.tf) call strncmp0(trim(rstr), tag_file, tf)
      if(tf) then
         initial_occmat = FILE
         goto 1001
      end if
1001  continue
    end subroutine set_initial_occmat

    subroutine getgmax()
      logical :: tf
      tf =  f_getRealValue( tag_cke_wavefunctions, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf2, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf3, dret, "rydberg") == 0
      if(tf) then
         if(dret < 0) call phase_error_with_msg(nfout,' !! illegal input value of cutoff_energy_for_wavefunctions'&
                                               ,__LINE__,__FILE__)
         gmax = sqrt(dret)
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** gmax = ",f8.4)') gmax
         gmax_org = gmax
      else
         call phase_error_with_msg(nfout,' !! no tag line for cutoff_energy_for_wavefunctions',__LINE__,__FILE__)
      end if
    end subroutine getgmax

    integer function getgmaxp()
      logical :: tf
      gmaxp_defined = .false.
      tf = f_getRealValue( tag_cke_chargedensity, dret, "rydberg") == 0
      if(.not.tf) tf =  f_getRealValue( tag_cke_cd, dret, "rydberg") == 0
      if(.not.tf) tf =  f_getRealValue( tag_cke_cd2, dret, "rydberg") == 0
      if(.not.tf) tf =  f_getRealValue( tag_cke_cd3, dret, "rydberg") == 0
      if(tf) then
         if(dret < 0) call phase_error_with_msg(nfout,' !! illegal input value of cutoff_energy_for_chargedensity'&
                                               ,__LINE__,__FILE__)
         gmaxp = sqrt(dret)
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** gmaxp = ",f8.4)') gmaxp
         getgmaxp = 1
         gmaxp_defined = .true.
      else
         getgmaxp = 0
!!$         stop ' !! no tag line for cutoff_energy_for_chargedensity'
      end if
    end function getgmaxp

  ! ------- Positron start
    subroutine getgmax_positron()
      logical :: tf
      tf = f_getRealValue( tag_cke_pwf, dret, "rydberg") == 0
      if(.not.tf) tf =  f_getRealValue( tag_cutoff_pwf, dret, "rydberg") == 0
      if(tf) then
         if(dret < 0) call phase_error_with_msg(nfout,' !! illegal input value of cutoff-energy for a positron wave function'&
                                               ,__LINE__,__FILE__)
         gmax_positron = sqrt(dret)
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** gmax_positron = ",f8.4)') gmax_positron
      else
         if(printable) &
              & write(nfout,'(" !! no tag line for cutoff-energy for a positron wave function")')
         gmax_positron = gmax
      end if
      if(gmax_positron > gmaxp*0.5d0) then
         gmax_positron = gmaxp*0.5d0
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** gmax_positron = ",f8.4)') gmax_positron
      end if
    end subroutine getgmax_positron
  ! ------- Positron end

    subroutine getgmaxs()
      logical :: tf
      tf = f_getRealValue( tag_cke_initial_matdiagon, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wavefunctions, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf,dret,"rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf2, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf3, dret, "rydberg") == 0
      if(tf) then
         if(dret < 0) call phase_error_with_msg(nfout,' !! illegal input value of cutoff_energy_for_initial_matdiagon'&
                                               ,__LINE__,__FILE__)
         gmaxs_given = sqrt(dret)
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** gmaxs_given = ",f8.4)') gmaxs_given
      else
!!$         stop ' !! no tag line for cutoff_energy_for_initial_matdiagon'
      end if
    end subroutine getgmaxs

    subroutine set_smearing_method(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_parabolic, len(tag_parabolic), tf)
      if(tf) then
         way_of_smearing = PARABOLIC
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_cold, len(tag_cold), tf)
      if(tf) then
         way_of_smearing = COLD
         goto 1001
      end if
! ================= KT_add ========================= 13.0E
      call strncmp2(rstr, FMAXVALLEN, tag_fermi_dirac, len(tag_fermi_dirac), tf)
      if(tf) then
         way_of_smearing = FERMI_DIRAC;  goto 1001
      end if
! ================================================== 13.0E
!      call strncmp2(rstr, FMAXVALLEN, tag_methfessel_paxton, len(tag_methfessel_paxton), tf)
      call strncmp2(rstr, FMAXVALLEN, tag_meth, len(tag_meth), tf)
      if(tf) then
         way_of_smearing = MP;  goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_improved_tetrahedron, len(tag_improved_tetrahedron), tf)
      if(tf) then
         way_of_smearing = TETRAHEDRON
         sw_correction = ON
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, tag_tetrahedron, len(tag_tetrahedron), tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_tetrahedral,tf)
      if(tf) then
         way_of_smearing = TETRAHEDRON
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_lowest_at_each_kpt, &
           &       len(tag_lowest_at_each_kpt), tf)
      if(tf) then
         way_of_smearing = LOWEST_AT_EACH_KPT
         goto 1001
      end if

1001  continue
    end subroutine set_smearing_method

    subroutine set_intzaj(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp0(tag_matrix_diagon,trim(rstr),tf)
      if(.not.tf .and. rstr(1:1) == '1') tf = .true.
      if(tf) then
         intzaj = by_matrix_diagon
         goto 1001
      end if
      call strncmp0(tag_random_numbers,trim(rstr),tf)
      if(.not.tf .and. rstr(1:1) == '0') tf = .true.
      if(tf) then
         intzaj = by_random_numbers
         goto 1001
      end if
      call strncmp0(tag_given_by_a_file,trim(rstr),tf)
      if(.not.tf) call strncmp0(tag_file,trim(rstr),tf)
      if(tf) then
         intzaj = FILE
         goto 1001
      end if
      call strncmp0(tag_atomic_orbital,trim(rstr),tf)
      if(tf) then
         intzaj = by_pseudo_atomic_orbitals
         goto 1001
      end if
1001  continue
    end subroutine set_intzaj

    subroutine set_initial_chg(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_Gauss, len(tag_Gauss), tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_Gauss,tf)
      if(tf) then
         initial_chg = Gauss_distrib_func
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_very_narrow, tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_Gauss_distrib_func_over_r,tf)
      if(tf) then
         initial_chg = VERY_NARROW
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_from_PseudoPotential_FILE,tf)
      if(.not.tf) call strncmp0(trim(rstr), tag_Atomic_charge_density,tf)
      if(tf) then
         initial_chg = from_PSEUDOPOTENTIAL_FILE
         goto 1001
      end if
      call strncmp0(trim(rstr), tag_given_by_a_file, tf)
      if(.not.tf) call strncmp0(trim(rstr), tag_file, tf)
      if(tf) then
         initial_chg = FILE
         goto 1001
      end if
1001  continue
    end subroutine set_initial_chg

#endif
  end subroutine m_CtrlP_rd_accuracy

  subroutine m_CtrlP_set_default_gmaxp(nfout,uspp)
    integer, intent(in) :: nfout
    logical, intent(in) :: uspp
    if(gmaxp_defined) return
    if(uspp) then
      gmaxp = gmax*3.d0
    else
      gmaxp = gmax*2.d0
    endif
    gmaxp_reduced = gmaxp*fftsize_factor_gmaxp
    if(mype==0) then
      write(nfout,'(a,l)')     ' !** has_uspp: ',uspp
      write(nfout,'(a,f10.5)') ' !** set default gmaxp : ',gmaxp
    endif
  end subroutine m_CtrlP_set_default_gmaxp

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_rd_wfsolver(nfout,natm2)
    integer, intent(in) :: nfout,natm2
    integer :: f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue
    integer :: f_selectParentBlock, f_selectTop
    integer :: iret,i,ba
    real(kind=DP) :: dret
    logical :: flag = .false., prealloc, number_is_given
    logical :: explict_solver = .false.

    iret = f_selectTop()
    ! --- wavefunction_solver ---
    flag = f_selectBlock( tag_wavefunction_solver ) == 0
    if(.not.flag) flag = f_selectBlock( tag_wf_solver ) == 0
    if(flag)then
       if(f_selectBlock(tag_solvers)==0) then
          explict_solver = .true.
!!$        explict_solver = f_selectBlock(tag_solvers)==0
          iret = f_selectParentBlock()
       else if(f_selectBlock(tag_for_init_str) == 0) then
          explict_solver = .true.
          iret = f_selectParentBlock()
       end if
    endif

    ! determine the default value for the davidson-related variables
    if(neg/nrank_e<4) then
       sw_divide_subspace_changed = .true.
       sw_npartition_changed = .true.
    else
       npartition_david = neg/(nblock*nrank_e)
       if (npartition_david<1) npartition_david = 1
       sw_npartition_changed = .true.
       if(printable) write(nfout,'(a,i8)') " !** REMARK: npartition_david was set to : ",npartition_david
    endif

    if(.not.explict_solver)then
       tag_solver_of_WF_is_found = .true.
       if(sw_fef == ON) solver_set = FEF
       call configure_wf_solver(solver_set,natm2)
       meg = neg
    endif

    if(flag) then
       if(explict_solver)then !explicit solver specification
       tag_solver_of_WF_is_found = .true.
       if(ipriinputfile >= 2 .and. printable) &
            & write(nfout,'(" !** -- tag_wavefunction_solver is found --")')
       ! --- count total number of solvers ---
       prealloc = .true.
       if( f_selectBlock( tag_for_init_str) == 0) then
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag_for_init_str is found --")')
          number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
          if(number_is_given) n_WF_solvers_before = iret
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          if(iret < 0) call phase_error_with_msg(nfout,&
          ' solver set (for_init_str) is not given properly <<m_CtrlP_rd_wfsolver>>',__LINE__,__FILE__)
          if(.not.number_is_given) n_WF_solvers_before = iret
          if(number_is_given .and. n_WF_solvers_before>iret) n_WF_solvers_before = iret
          iret = f_selectParentBlock()

          if( f_selectBlock( tag_during_str_relax) == 0) then
             number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
             if(number_is_given) n_WF_solvers_after = iret
             ba = AFTER
             call set_wfsolvers(prealloc,n_WF_solvers_after,n_WF_solvers_before,ba,iret)
             if(iret < 0) call phase_error_with_msg(nfout,&
             ' solver set (during_str_relax) is not given properly <<m_CtrlP_rd_wfsolver>>',__LINE__,__FILE__)
             if(.not.number_is_given) n_WF_solvers_after = iret
             if(number_is_given .and. n_WF_solvers_after > iret) n_WF_solvers_after = iret
             iret = f_selectParentBlock()
          else
             n_WF_solvers_after = 0
          end if
       else
          number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
          if(number_is_given) n_WF_solvers_before = iret
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          if(iret < 0) call phase_error_with_msg(nfout,' solver set is not given properly <<m_CtrlP_rd_wfsolver>>'&
                                                ,__LINE__,__FILE__)
          if(.not.number_is_given) n_WF_solvers_before = iret
          if(number_is_given .and. n_WF_solvers_before > iret) n_WF_solvers_before = iret
          n_WF_solvers_after = 0
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** n_WF_solvers_before, n_WF_solvers_after = ",2i6)') &
               &  n_WF_solvers_before, n_WF_solvers_after

!!$          iret = f_selectParentBlock()
       end if

       !----------------------------------------------------------
       n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
       !----------------------------------------------------------
       if(ipriinputfile >= 2 .and. printable) &
            & write(nfout,'(" !** number of all solvers = ",i6)') n_WF_solvers_all
       ! -- allocation and initialization of "w_solver"
       call alloc_w_solver(n_WF_solvers_all)
       ! --- substitution for w_solver
       prealloc = .false.
       if( f_selectBlock( tag_for_init_str) == 0) then
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          iret = f_selectParentBlock()
          if( f_selectBlock( tag_during_str_relax) == 0) then
             ba = AFTER
             call set_wfsolvers(prealloc,n_WF_solvers_after,n_WF_solvers_before,ba,iret)
             iret = f_selectParentBlock()
          end if
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag_for_init_str is not found again --")')
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_all,0,ba,iret)
       end if

       endif !explicit solver specification

       ! ---- lineminimization ---
       flag = f_selectBlock( tag_lineminimization) == 0
       if(.not.flag) flag = f_selectBlock( tag_line_minimization) == 0
       if(flag) then
          if( f_getRealValue( tag_dt_lower_critical, dret, 'au_time') == 0) then
             dt_Lower_CRITICAL = dret
          end if
          if( f_getRealValue( tag_dt_upper_critical, dret, 'au_time') == 0) then
             if(dret >  dt_Lower_CRITICAL .and. dret > 1.d-3) then
                dt_Upper_CRITICAL = dret
             else
                if(printable) write(nfout &
                     & ,'(" !** the given value of dt_Upper_CRITICAL is invalid: dt_Upper_CRITICAL = ",d15.5)') dret
             end if
          end if
! --> T. Yamasaki, 15th July 2008
          if( f_getRealValue( tag_dt_upper_factor, dret,'') == 0) then
             dt_upper_factor = dret
             if(dt_upper_factor < 1.0) dt_upper_factor = 1.0
             dt_upper_factor_is_set = .true.
          end if
          if( f_getRealValue( tag_dt_lower_factor, dret,'') == 0) then
             dt_lower_factor = dret
             if(dt_lower_factor > 1.0) dt_lower_factor = 1.0
          end if
! <---

          if( f_getRealValue( tag_delta_lmdenom, dret,'') == 0) delta_lmdenom = dret

          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** -- tag for lineminimization is found --")')
             write(nfout,'(" !** dt_lower_critical, dt_upper_critical = ",2f12.8)') dt_lower_critical, dt_upper_critical
             write(nfout,'(" !** delta_lmdenom = ",d20.8)') delta_lmdenom
          end if

! --> T. Yamasaki, 19th June 2009
          if( f_getStringValue( tag_energy_evaluation, rstr, LOWER) == 0) then
             call set_energy_evaluation(rstr) ! --> energy_evaluation
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** -- tag_energy_evaluation is found --")')
                write(nfout,'(" !** energy_evaluation = ",i3," : 1=TOTAL_ENERGY, 2=MODIFIED_TOTAL_ENERGY, 3=BAND_ENERGY")') &
                     & energy_evaluation
             end if
          else
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** -- tag_energy_evaluation is not found --")')
             end if
          end if
          if(energy_evaluation == BAND_ENERGY .or. energy_evaluation == MODIFIED_TOTAL_ENERGY) then
             if( f_getIntValue(tag_num_conduction_bands, iret)== 0) then
                num_conduction_bands_lmm = iret
                if(ipriinputfile >= 1 .and. printable) then
                   write(nfout,'(" !** -- tag_num_conduction_bands is found --")')
                   write(nfout,'(" !**  num_conduction_bands = ",i8)') num_conduction_bands_lmm
                end if
             end if
          end if
! <--
! --> T. Yamasaki, 18th Aug. 2009
          sw_lmm_status_check = -1
          if( f_getIntValue(tag_lmm_status_check,iret) == 0) then
             sw_lmm_status_check = iret
          else if(f_getIntValue(tag_sw_lmm_status_check,iret) == 0) then
             sw_lmm_status_check = iret
          end if
          if(ipriinputfile >= 1 .and. printable) then
             if(sw_lmm_status_check /= -1) &
                  & write(nfout,'(" !** -- tag_sw_lmm_status_check is found --")')
             if(sw_lmm_status_check == -1) sw_lmm_status_check = NO
             write(nfout,'(" !**  sw_lmm_status_check = ",i3)') sw_lmm_status_check
          end if
! <--
          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag for lineminimization is not found --")')
       end if

       ! ---- rmm ---
       if( explict_solver )then
          edelta_change_to_rmm = 1.e-3/dble(natm2)
          edelta_change_to_rmm_md = 1.e-3/dble(natm2)
       endif

       if( f_selectBlock( tag_rmm) == 0) then
          if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_rmm is found --")')
          if( f_getIntValue(tag_imGSrmm, iret) == 0) imGSrmm = iret
          if( f_getRealValue(tag_rr_critical_value,dret,'') == 0) rr_Critical_Value = dret
          if( f_getIntValue(tag_rmm_precal_phase_matm,iret)==0) rmm_precal_phase_matm = iret
!!$          if( f_getIntValue(tag_rmm3_bisec_trial_max, iret)==0) rmm3_bisec_trial_max = iret
!!$          if( f_getRealValue(tag_rmm3_bisec_crtcl_value,dret,'') == 0) rmm3_bisec_crtcl_value = dret
          if( f_getRealValue(tag_edelta_change_to_rmm,dret,'hartree')==0) then
             edelta_change_to_rmm = dret
             edelta_change_to_rmm_md = dret
             edelta_rmm_given = .true.
          endif
          if(f_getRealValue(tag_edelta_change_to_rmm_md,dret,'hartree')==0)then
             edelta_change_to_rmm_md = dret
             edelta_rmm_given = .true.
          endif
          if( f_getRealValue(tag_delta_residual_occup,dret,'') == 0) delta_residual_occup = dret
          if( f_getRealValue(tag_delta_residual_empty,dret,'') == 0) delta_residual_empty = dret

          if( f_getIntValue(tag_save_memory_mode, iret) == 0) rmm_save_memory_mode = iret
          iret = f_selectParentBlock()
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** --- parameters for rmm ---")')
             write(nfout,'(" !** rr_Critical_Value = ",d20.8)') rr_Critical_Value
             write(nfout,'(" !** edelta_change_to_rmm = ",d20.8)') edelta_change_to_rmm
             write(nfout,'(" !** rmm_precal_phase_matm = ",i10)') rmm_precal_phase_matm
             write(nfout,'(" !** save_memory_mode = ",i10)') rmm_save_memory_mode
             write(nfout,'(" !** delta_residual_occup = ",d20.8)') delta_residual_occup
             write(nfout,'(" !** delta_residual_empty = ",d20.8)') delta_residual_empty
          end if
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !* tag_rmm is not found")')
       end if

       ! --- subspace_rotation ---
       flag = f_selectBlock( tag_submat) == 0
       if(.not.flag) flag = f_selectBlock( tag_subspacerotation) == 0
       if(.not.flag) flag = f_selectBlock( tag_subspace_rotation) == 0
       if( flag ) then
          meg = neg
          if(ipriinputfile >= 1 .and. printable) &
               & write(nfout,'(" !** tag_subspace_rotation is found")')
          if( f_getIntValue(tag_subspace_matrix_size,iret) == 0) meg = iret
          if( f_getRealValue(tag_damping_factor,dret,'') == 0) damp = dret
          if( f_getIntValue(tag_before_renewal,iret) == 0) submat_before_renewal = iret
          if( meg > neg .or. meg < 1) then
             if(printable) then
                write(nfout,'(" !** given subspace_matrix_size(=meg) = ",i6)') meg
                write(nfout,'(" !*  subspace_matrix_size (meg) is set to be neg (= ",i6,")")') neg
             end if
             meg = neg
          end if
          if(damp > 1.0 .or. damp < 0.0) then
             if(printable) then
                write(nfout,'(" !** given damping factor of subspace rotation = ",f8.4)') damp
                write(nfout,'(" !*  damping factor (damp) is set to be 1.0")')
             end if
             damp = 1.d0
          end if

          flag =  f_getIntValue(tag_period,iret) == 0
                    if(ipriinputfile >= 2 .and. flag .and. printable) &
                         & write(nfout,'(" !* tag_period is found")')
          if(.not.flag) then
             flag = f_getIntValue(tag_one_period,iret) == 0
                    if(ipriinputfile >= 2 .and. flag .and. printable) &
                         & write(nfout,'(" !* tag_one_period is found")')
          end if
          if(flag) submat_period = iret

          if( f_getRealValue(tag_critical_ratio,dret,'') == 0) then
             submat_critical_ratio = dret
             if(ipriinputfile >= 2 .and. printable) &
                  & write(nfout,'(" !* tag_critical_ratio is found")')
          end if

          if( f_getIntValue(tag_sw_gep, iret) == 0) sw_gep = iret
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** --- parameters for subspace roation ---")')
             write(nfout,'(" !** subspace_matrix_size = ",i6)') meg
             write(nfout,'(" !** damping_factor       = ",f8.4)') damp
             write(nfout,'(" !** submat_period        = ",i6)') submat_period
             write(nfout,'(" !** submat_critical_ratio= ",d20.8)') submat_critical_ratio
             if(submat_before_renewal == ON) then
                write(nfout,'(" !** before_renewal= ON")')
             else
                write(nfout,'(" !** before_renewal= OFF")')
             end if
             write(nfout,'(" !** sw_gep = ",i6)') sw_gep
          end if

!!$!BRANCH_P 3D_Parallel
!!$          if(nrank_k>=2) sw_scalapack = OFF
!!$             ! This is a tentative default setting until scalapack parallelization is completed for nrank_k>=2
!!$!BRANCH_P_END 3D_Parallel
          if( f_selectBlock( tag_scalapack) == 0) then
             if( f_getIntValue(tag_sw_scalapack, iret) == 0) sw_scalapack = iret
#ifndef _USE_SCALAPACK_
             if(sw_scalapack == ON) then
                write(nfout,*) 'ScaLAPACK is unusable.'
                write(nfout,*) 'Use the CPP option -D_USE_SCALAPACK_.'
                call phase_error_with_msg(nfout,'ScaLAPACK is unusable.',__LINE__,__FILE__)
             end if
#endif
             call set_method_scalapack(method_scalapack)
             if( f_getIntValue(tag_block_size, iret) == 0) block_size = iret
             if( f_getIntValue(tag_nprow, iret) == 0) nprow = iret
             if( f_getIntValue(tag_npcol, iret) == 0) npcol = iret
             if( f_getIntValue(tag_divide_square, iret) == 0) divide_square = iret
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** tag_scalapack is found")')
                write(nfout,'(" !** sw_scalapack = ",i6)') sw_scalapack
                write(nfout,'(" !** method = ",i6)') method_scalapack
                write(nfout,'(" !** block_size = ",i6)') block_size
                write(nfout,'(" !** nprow = ",i6)') nprow
                write(nfout,'(" !** npcol = ",i6)') npcol
                write(nfout,'(" !** divide_square = ",i6)') divide_square
             end if
             iret = f_selectParentBlock()
          end if


          if( f_selectBlock( tag_scalapack) == 0) then
             if( f_getIntValue(tag_sw_scalapack, iret) == 0) sw_scalapack = iret
             call set_method_scalapack(method_scalapack)
             if( f_getIntValue(tag_block_size, iret) == 0) block_size = iret
             if( f_getIntValue(tag_nprow, iret) == 0) nprow = iret
             if( f_getIntValue(tag_npcol, iret) == 0) npcol = iret
             if( f_getIntValue(tag_memory, iret) == 0) msize_submat = iret
             if( f_getIntValue(tag_divide_square, iret) == 0) divide_square = iret
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** tag_scalapack is found")')
                write(nfout,'(" !** sw_scalapack = ",i6)') sw_scalapack
                write(nfout,'(" !** method = ",i6)') method_scalapack
                write(nfout,'(" !** block_size = ",i6)') block_size
                write(nfout,'(" !** nprow = ",i6)') nprow
                write(nfout,'(" !** npcol = ",i6)') npcol
                write(nfout,'(" !** memory(MB) = ",i6)') msize_submat
                write(nfout,'(" !** divide_square = ",i6)') divide_square
             end if
             iret = f_selectParentBlock()
          end if

#ifdef USE_EIGENLIB
          if( f_selectBlock(tag_eigen_exa) == 0) then
             if( f_getIntValue(tag_sw_eigen_exa, iret) == 0) sw_eigen_exa = iret
             if(sw_eigen_exa == ON) then
               sw_scalapack = ON
             endif
             write(nfout,'(" !** sw_eigen_exa =", i6)'),sw_eigen_exa
             iret = f_selectParentBlock()
          endif
#endif

          iret = f_selectParentBlock()
       else
          meg = neg
       end if

       ! ---- Davidson ---
       if(ekmode == ON .or. &
         & (icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION)) then
          sw_first_conv_check = OFF
       end if
       if( f_selectBlock( tag_davidson) == 0) then
          if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_davidson is found --")')
          if( f_getIntValue(tag_ndavid,iret)==0) ndavid = iret
          if( f_getIntValue(tag_max_iter_david,iret)==0) max_iter_david = iret
          if( f_getIntValue(tag_max_subspace_size,iret)==0) max_subspace_size = iret
          if( f_getIntValue(tag_sw_first_conv_check,iret)==0) sw_first_conv_check = iret
          if(ndavid<0) ndavid=abs(ndavid)
          if(max_iter_david<0) max_iter_david=abs(max_iter_david)
          if( max_subspace_size < 4*neg) max_subspace_size = 4*neg ! default value
          if( f_getRealValue(tag_delta_eig_occup,dret,"hartree")==0) delta_eig_occup = dret
          if( f_getRealValue(tag_delta_eig_empty,dret,"hartree")==0) delta_eig_empty = dret
          if( f_getRealValue(tag_eps_david,dret,"")==0) eps_david = dret
!! md-david, md-kosugi
          if( f_getIntValue(tag_submat_GE,iret)==0) submat_GE = iret
          if( f_getIntValue(tag_sw_MRCV_only,iret)==0) sw_MRCV_only = iret
          if( f_getIntValue(tag_sw_divide_subspace,iret)==0) then
              sw_divide_subspace = iret
              sw_divide_subspace_changed = .false.
          endif
          if( f_getIntValue(tag_npartition,iret)==0.or.f_getIntValue(tag_npartition_david,iret)==0) then
              npartition_david = iret
              sw_npartition_changed = .false.
          endif
          if( f_getIntValue(tag_nbands_in_partition,iret)==0 .or. &
            & f_getIntValue(tag_nblock_david,iret)==0 .or. &
            & f_getIntValue(tag_nblock,iret)==0) then
             if(iret>0) then
                npartition_david = neg/(iret*nrank_e)
                if (npartition_david<1) npartition_david = 1
                sw_npartition_changed = .false.
             endif
          endif

          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !** ndavid = ",i10)') ndavid
             write(nfout,'(" !** max_iter_david = ",i10)') max_iter_david
             write(nfout,'(" !** max_subspace_size = ",i10)') max_subspace_size
             write(nfout,'(" !** delta_eig_occup = ",e12.5)') delta_eig_occup
             write(nfout,'(" !** delta_eig_empty = ",e12.5)') delta_eig_empty
             write(nfout,'(" !** eps_david = ",e12.5)') eps_david
             write(nfout,'(" !** sw_first_conv_check = ",i10)') sw_first_conv_check
             write(nfout,'(" !** sw_MRCV_only = ",i10)') sw_MRCV_only
             write(nfout,'(" !** submat_GE = ",i10)') submat_GE
             write(nfout,'(" !** sw_divide_subspace = ",i10)') sw_divide_subspace
             write(nfout,'(" !** npartition = ",i10)') npartition_david
          end if
          if( f_selectBlock( tag_dav_scalapack) == 0) then
             if( f_getIntValue(tag_sw_dav_scalapack, iret) == 0) sw_dav_scalapack = iret
#ifndef _USE_SCALAPACK_
             if(sw_dav_scalapack == ON) then
                write(nfout,*) 'ScaLAPACK is unusable.'
                write(nfout,*) 'Use the CPP option -D_USE_SCALAPACK_.'
                call phase_error_with_msg(nfout,'ScaLAPACK is unusable.',__LINE__,__FILE__)
             end if
#endif
             if( f_getIntValue(tag_dav_block_size, iret) == 0) dav_block_size = iret
             if( f_getIntValue(tag_dav_nprow, iret) == 0) dav_nprow = iret
             if( f_getIntValue(tag_dav_npcol, iret) == 0) dav_npcol = iret
             if( f_getIntValue(tag_dav_divide_square, iret) == 0) dav_divide_square = iret
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** tag_dav_scalapack is found")')
                write(nfout,'(" !** sw_dav_scalapack = ",i6)') sw_dav_scalapack
                write(nfout,'(" !** dav_block_size = ",i6)') dav_block_size
                write(nfout,'(" !** dav_nprow = ",i6)') dav_nprow
                write(nfout,'(" !** dav_npcol = ",i6)') dav_npcol
                write(nfout,'(" !** dav_divide_square = ",i6)') dav_divide_square
             end if
             iret = f_selectParentBlock()
          end if
          iret = f_selectParentBlock()
       else
          max_subspace_size = 4*neg ! default value
!!$          if(ipriinputfile >= 2 .and. printable) &
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** tag_davidson is not found")')
             write(nfout,'(" !** max_subspace_size = ",i6)') max_subspace_size
          end if
       end if

       ! ---- Modified Davidson ---
       if( f_selectBlock( tag_mddavidson) == 0) then
          if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_mddavidson is found --")')
          if( f_getIntValue(tag_npartition_mddavid,iret)==0) npartition_mddavid = iret
          if( f_getIntValue(tag_max_iter_mddavid,iret)==0) max_iter_mddavid = iret
          if(npartition_mddavid<0) npartition_mddavid=abs(npartition_mddavid)
          if(max_iter_mddavid<0) max_iter_mddavid=abs(max_iter_mddavid)
          if( f_getRealValue(tag_delta_eig_occup,dret,"hartree")==0) delta_eig_occup_md = dret
          if( f_getRealValue(tag_delta_eig_empty,dret,"hartree")==0) delta_eig_empty_md = dret
          if( f_getRealValue(tag_eps_mddavid,dret,"")==0) eps_mddavid = dret
          if( f_getRealValue(tag_eps_residual_mddavid,dret,"")==0) eps_residual_mddavid = dret
          if( f_getIntValue(tag_sw_apply_gs,iret)==0 ) sw_apply_gs = iret
          max_iter_david = max_iter_mddavid
          delta_eig_occup = delta_eig_occup_md
          delta_eig_empty = delta_eig_empty_md
          eps_david = eps_mddavid
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !** npartition_mddavid = ",i10)') npartition_mddavid
             write(nfout,'(" !** max_iter_mddavid = ",i10)') max_iter_mddavid
             write(nfout,'(" !** delta_eig_occup = ",e12.5)') delta_eig_occup_md
             write(nfout,'(" !** delta_eig_empty = ",e12.5)') delta_eig_empty_md
             write(nfout,'(" !** eps_mddavid = ",e12.5)') eps_mddavid
             write(nfout,'(" !** eps_residual_mddavid = ",e12.5)') eps_residual_mddavid
          end if
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** tag_mddavidson is not found")')
       end if

       ! ---- Modified Kosugi ---
       if( f_selectBlock( tag_mdkosugi) == 0) then
          if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_mdkosugi is found --")')
          if( f_getIntValue(tag_npartition_mdkosugi,iret)==0) npartition_mdkosugi = iret
          if( f_getIntValue(tag_max_iter_mdkosugi,iret)==0) max_iter_mdkosugi = iret
          if(npartition_mdkosugi<0) npartition_mdkosugi=abs(npartition_mdkosugi)
          if(max_iter_mdkosugi<0) max_iter_mdkosugi=abs(max_iter_mdkosugi)
          if( f_getRealValue(tag_delta_eig_occup_mdkosugi,dret,"hartree")==0) delta_eig_occup_mdkosugi = dret
          if( f_getRealValue(tag_delta_eig_empty_mdkosugi,dret,"hartree")==0) delta_eig_empty_mdkosugi = dret
          if( f_getRealValue(tag_eps_mdkosugi,dret,"")==0) eps_mdkosugi = dret
          if( f_getRealValue(tag_eps_residual_mdkosugi,dret,"")==0) eps_residual_mdkosugi = dret
          if( f_getIntValue(tag_sw_apply_gs,iret)==0 ) sw_apply_gs = iret
          max_iter_david = max_iter_mdkosugi
          delta_eig_occup = delta_eig_occup_mdkosugi
          delta_eig_empty = delta_eig_empty_mdkosugi
          eps_david = eps_mdkosugi
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !** npartition_mdkosugi = ",i10)') npartition_mdkosugi
             write(nfout,'(" !** max_iter_mdkosugi = ",i10)') max_iter_mdkosugi
             write(nfout,'(" !** delta_eig_occup_mdkosugi = ",e12.5)') delta_eig_occup_mdkosugi
             write(nfout,'(" !** delta_eig_empty_mdkosugi = ",e12.5)') delta_eig_empty_mdkosugi
             write(nfout,'(" !** eps_mdkosugi = ",e12.5)') eps_mdkosugi
             write(nfout,'(" !** eps_residual_mdkosugi = ",e12.5)') eps_residual_mdkosugi
          end if
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** tag_mdkosugi is not found")')
       end if
       ! ---- CG  ---
       if( f_selectBlock( tag_cg) == 0) then
         if(f_getIntValue(tag_sw_modified_cg_formula,iret)==0) sw_modified_cg_formula=iret
         if(ipriinputfile>=2 .and. printable) &
         & write(nfout,'(a,i3)') 'modified formula for CG : ',sw_modified_cg_formula
       endif

       iret = f_selectParentBlock()
    else
       tag_solver_of_WF_is_found = .true.
       if(ipriinputfile >=2 .and. printable) &
            & write(nfout,'(" !** -- tag_wavefunction_solver is not found --")')
       !n_WF_solvers_before = 1
       !n_WF_solvers_after = 0
       !n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
       !call alloc_w_solver(n_WF_solvers_all)
       !if(ipriinputfile >= 2 .and. printable) &
       !     & write(nfout,'(" !** n_WF_solvers_before, n_WF_solvers_after = ",2i6)') &
       !     &  n_WF_solvers_before, n_WF_solvers_after
       meg = neg
    end if
    if(ipriinputfile >= 1 .and. printable) then
       write(nfout,'(" !** --- id, sol, till_n,  dts, dte, itr, var, prec, cmix, submat ---")')
       do i = 1, n_WF_solvers_all
          if(i == 1) write(nfout,'(" !**      - for_init_str -")')
          if(i == n_WF_solvers_before+1) write(nfout,'(" !**      - during_str_relax -")')
          write(nfout,'(" !** ",i3,i4,i6,f7.3,f7.3,i4,i3,i3,i3,i3)') i, w_solver(i)%solver, w_solver(i)%till_n_iter &
               &              , w_solver(i)%dtim_s, w_solver(i)%dtim_e, w_solver(i)%iter_range &
               &              , w_solver(i)%variation_way, w_solver(i)%precon, w_solver(i)%cmix_pointer &
               &              , w_solver(i)%subspace_rotation
       end do
    end if
  contains
    subroutine set_energy_evaluation(rstr)
      character(len=FMAXVALLEN), intent(in) :: rstr
      logical :: tf
      energy_evaluation = TOTAL_ENERGY
      call strncmp0(tag_total_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = TOTAL_ENERGY
         goto 1001
      end if
      call strncmp0(tag_modified_total_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = MODIFIED_TOTAL_ENERGY
         goto 1001
      end if
      call strncmp0(tag_band_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = BAND_ENERGY
         goto 1001
      end if
      call phase_error_with_msg(nfout,' ! tag for energy_evaluation is invalid <<m_CtrlP_rd_wfsolver.set_energy_evaluation>>'&
                               ,__LINE__,__FILE__)
1001  continue
    end subroutine set_energy_evaluation

    subroutine configure_wf_solver(solver_set,natm2)
       integer, intent(in) :: solver_set,natm2
       integer :: i
       if(solver_set == LMM_RMM)then
          call alloc_w_solver(2)
          w_solver(1)%solver = lmMSD
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = 5
          w_solver(2)%solver = RMM3
          w_solver(2)%till_n_iter = -1
          w_solver(2)%subspace_rotation = ON
          edelta_change_to_rmm = 1.d-4/dble(natm2)
          edelta_change_to_rmm_md = 1.d-4/dble(natm2)
          n_WF_solvers_before = 2
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : lm+msd -> rmm3")')
       else if (solver_set == DAV_RMM)then
          call alloc_w_solver(4)
          if(icond==INITIAL .or. icond==CONTINUATION .or. icond==AUTOMATIC)then
             w_solver(1)%solver = MDDAVIDSON
             w_solver(3)%solver = MDDAVIDSON
          else
             w_solver(1)%solver = MDKOSUGI
             w_solver(3)%solver = MDKOSUGI
          endif
          if(sw_hubbard==ON) then
             w_solver(1)%solver = MDKOSUGI
             w_solver(3)%solver = MDKOSUGI
          endif
! === KT_add === 2015/01/05
          if ( noncol ) then
             w_solver(1)%solver = MDDAVIDSON
             w_solver(3)%solver = MDDAVIDSON
          endif
! ============== 2015/01/05

!          w_solver(1)%solver = MDKOSUGI
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = 5
          w_solver(1)%precon = ON
          w_solver(1)%before_or_after_convergence = BEFORE
          w_solver(2)%solver = RMM3
          w_solver(2)%precon = ON
          w_solver(2)%till_n_iter = -1
          w_solver(2)%subspace_rotation = ON
          w_solver(2)%before_or_after_convergence = BEFORE
          !!w_solver(3)%solver = MDDAVIDSON
          w_solver(3)%subspace_rotation = ON
          w_solver(3)%till_n_iter = 5
          w_solver(3)%precon = ON
          w_solver(3)%before_or_after_convergence = AFTER
          w_solver(4)%solver = RMM3
          w_solver(4)%precon = ON
          w_solver(4)%till_n_iter = -1
          w_solver(4)%subspace_rotation = ON
          w_solver(4)%before_or_after_convergence = AFTER
          edelta_change_to_rmm = 1.d-3/dble(natm2)
          edelta_change_to_rmm_md = 1.d-3/dble(natm2)
          n_WF_solvers_before = 2
          n_WF_solvers_after = 2
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
!          if(neg/nrank_e<4) then
!             sw_divide_subspace=OFF
!             sw_divide_subspace_changed = .true.
!             sw_npartition_changed = .true.
!             if(printable) write(nfout,'(" !** REMARK: sw_divide_subspace is set to OFF ")')
!          else
!             npartition_david = neg/(nblock*nrank_e)
!             if (npartition_david<1) npartition_david = 1
!             sw_npartition_changed = .true.
!             if(printable) write(nfout,'(a,i8)') " !** REMARK: npartition_david was set to : ",npartition_david
!          endif
          !!$if(sw_hubbard==ON.or.nspin>1) sw_divide_subspace=OFF
          !!$if(sw_hubbard==ON) sw_divide_subspace=OFF
          if (printable) write(nfout,'(" !** applied wavefunction solver set : davidson -> rmm3")')
       else if (solver_set == lmMSD) then
          call alloc_w_solver(1)
          w_solver(1)%solver = lmMSD
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : lm+msd")')
       else if (solver_set == DAVIDSON)then
          call alloc_w_solver(1)
          w_solver(1)%solver = DAVIDSON
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%precon = ON
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if(neg/nrank_e<4) then
             sw_divide_subspace=OFF
             sw_divide_subspace_changed = .true.
             write(nfout,'(" !** REMARK: sw_divide_subspace is set to OFF ")')
          endif
          !!$if (sw_hubbard==ON.or.nspin>1) sw_divide_subspace=OFF
          if (printable) write(nfout,'(" !** applied wavefunction solver set : davidson")')
       else if (solver_set == FEF)then
          call alloc_w_solver(1)
          w_solver(1)%solver = CG
          w_solver(1)%subspace_rotation = OFF
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : FEF (cg)")')
       endif
       !if (n_WF_solvers_before>1)then
       !   if(intzaj == by_matrix_diagon)then
       !      do i=1,n_WF_solvers_before-1
       !         w_solver(i)%till_n_iter = i+4
       !      enddo
       !   else
       !      w_solver(1)%till_n_iter = 5
       !   endif
       !endif
    end subroutine configure_wf_solver

    subroutine set_wfsolvers(prealloc,msol,nbase,ba,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: msol,nbase,ba
      integer, intent(out):: iret

      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: i, sol, till_n, itr, var, prec,cmix, no, sw_submat
      integer, dimension(9) :: icf
      real(kind=DP) :: dts,dte
      if( f_selectBlock(tag_solvers) == 0) then
         i = 1
         do while(.true.)
            if(i == 1) then
               if( f_selectFirstTableLine() /= 0 ) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0 ) then
                  exit
               end if
            end if
            if(.not.prealloc) then
               if(i > msol) exit
               no = i
               iret = f_readsolver(icf,no,sol,till_n,dts,dte,itr,var,prec,cmix,sw_submat)

               if(no <= msol) then
                  w_solver(no+nbase)%before_or_after_convergence = ba
                  if(icf(1) == 1) w_solver(no+nbase)%solver = sol
!!$                  if(ekmode == ON) then
!!$                     if(w_solver(no+nbase)%solver == lmMSD) then
!!$                        w_solver(no+nbase)%solver = MSD
!!$                     end if
!!$                  end if
                  if(icf(2) == 1) w_solver(no+nbase)%till_n_iter = till_n
                  if(icf(3) == 1) w_solver(no+nbase)%dtim_s = dts
                  if(icf(4) == 1) then
                     w_solver(no+nbase)%dtim_e = dte
                  else if(icf(3) == 1) then
                     w_solver(no+nbase)%dtim_e = w_solver(no+nbase)%dtim_s
                  end if
                  if(icf(5) == 1) w_solver(no+nbase)%iter_range = itr
                  if(icf(6) == 1) w_solver(no+nbase)%variation_way = var
                  if(icf(7) == 1) w_solver(no+nbase)%precon = prec
                  if(icf(8) == 1) w_solver(no+nbase)%cmix_pointer = cmix
                  if(icf(9) == 1) w_solver(no+nbase)%subspace_rotation = sw_submat
               end if
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      else
         call phase_error_with_msg(nfout,' ! No wf solver is given in the inputfile <<m_CtrlP_rd_wfsolver.set_wfsolvers>>'&
                                  ,__LINE__,__FILE__)
      end if
        iret = i-1
    end subroutine set_wfsolvers

  end subroutine m_CtrlP_rd_wfsolver

  subroutine m_CtrlP_rd_wfsolver2(nfout,natm2, reread)
    integer, intent(in) :: nfout,natm2
    logical, intent(in), optional :: reread
    integer :: f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue
    integer :: f_selectParentBlock, f_selectTop
    integer :: iret,i,ba
    real(kind=DP) :: dret
    logical :: flag = .false., prealloc, number_is_given
    logical :: explict_solver = .false.
    logical :: rr, done_something
    rr = .false.
    if(present(reread)) rr = .true.

    iret = f_selectTop()
    ! --- wavefunction_solver ---
    flag = f_selectBlock( tag_wavefunction_solver ) == 0
    if(.not.flag) flag = f_selectBlock( tag_wf_solver ) == 0
    if(flag)then
       if(f_selectBlock(tag_solvers)==0) then
          explict_solver = .true.
!!$        explict_solver = f_selectBlock(tag_solvers)==0
          iret = f_selectParentBlock()
       else if(f_selectBlock(tag_for_init_str) == 0) then
          explict_solver = .true.
          iret = f_selectParentBlock()
       end if
    endif

    ! determine the default value for the davidson-related variables
    if(neg/nrank_e<4) then
       sw_divide_subspace_changed = .true.
       sw_npartition_changed = .true.
    else
       if(npartition_david /= neg/(nblock*nrank_e)) then
         npartition_david = neg/(nblock*nrank_e)
         if (npartition_david<1) npartition_david = 1
         sw_npartition_changed = .true.
         if(printable .and. .not. rr) write(nfout,'(a,i8)') &
         &  " !** REMARK: npartition_david was set to : ",npartition_david
       endif
    endif

    if(.not.explict_solver)then
       tag_solver_of_WF_is_found = .true.
       if(sw_fef == ON) solver_set = FEF
       if(.not. rr) call configure_wf_solver(solver_set,natm2)
       meg = neg
    endif

    if(flag) then
       if(explict_solver)then !explicit solver specification
       tag_solver_of_WF_is_found = .true.
       if(ipriinputfile >= 2 .and. printable) &
            & write(nfout,'(" !** -- tag_wavefunction_solver is found --")')
       ! --- count total number of solvers ---
       prealloc = .true.
       if( f_selectBlock( tag_for_init_str) == 0) then
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag_for_init_str is found --")')
          number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
          if(number_is_given) n_WF_solvers_before = iret
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          if(iret < 0) call phase_error_with_msg(nfout,&
          ' solver set (for_init_str) is not given properly <<m_CtrlP_rd_wfsolver>>',__LINE__,__FILE__)
          if(.not.number_is_given) n_WF_solvers_before = iret
          if(number_is_given .and. n_WF_solvers_before>iret) n_WF_solvers_before = iret
          iret = f_selectParentBlock()

          if( f_selectBlock( tag_during_str_relax) == 0) then
             number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
             if(number_is_given) n_WF_solvers_after = iret
             ba = AFTER
             call set_wfsolvers(prealloc,n_WF_solvers_after,n_WF_solvers_before,ba,iret)
             if(iret < 0) call phase_error_with_msg(nfout,&
             ' solver set (during_str_relax) is not given properly <<m_CtrlP_rd_wfsolver>>',__LINE__,__FILE__)
             if(.not.number_is_given) n_WF_solvers_after = iret
             if(number_is_given .and. n_WF_solvers_after > iret) n_WF_solvers_after = iret
             iret = f_selectParentBlock()
          else
             n_WF_solvers_after = 0
          end if
       else
          number_is_given = f_getIntValue( tag_num_solvers, iret) == 0
          if(number_is_given) n_WF_solvers_before = iret
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          if(iret < 0) call phase_error_with_msg(nfout,' solver set is not given properly <<m_CtrlP_rd_wfsolver>>'&
                                                ,__LINE__,__FILE__)
          if(.not.number_is_given) n_WF_solvers_before = iret
          if(number_is_given .and. n_WF_solvers_before > iret) n_WF_solvers_before = iret
          n_WF_solvers_after = 0
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** n_WF_solvers_before, n_WF_solvers_after = ",2i6)') &
               &  n_WF_solvers_before, n_WF_solvers_after

!!$          iret = f_selectParentBlock()
       end if
       if (rr .and. &
         n_WF_solvers_all /= (n_WF_solvers_before+n_WF_solvers_after)) then
         call dealloc_w_solver()
         n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
         call alloc_w_solver(n_WF_solvers_all)
       endif
       if (.not. rr) then
         !----------------------------------------------------------
         n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
         !----------------------------------------------------------
         if(ipriinputfile >= 2 .and. printable) &
              & write(nfout,'(" !** number of all solvers = ",i6)') n_WF_solvers_all
         ! -- allocation and initialization of "w_solver"
         call alloc_w_solver(n_WF_solvers_all)
         ! --- substitution for w_solver
       endif
       prealloc = .false.
       if( f_selectBlock( tag_for_init_str) == 0) then
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_before,0,ba,iret)
          iret = f_selectParentBlock()
          if( f_selectBlock( tag_during_str_relax) == 0) then
             ba = AFTER
             call set_wfsolvers(prealloc,n_WF_solvers_after,n_WF_solvers_before,ba,iret)
             iret = f_selectParentBlock()
          end if
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- tag_for_init_str is not found again --")')
          ba = BEFORE
          call set_wfsolvers(prealloc,n_WF_solvers_all,0,ba,iret)
       end if

       endif !explicit solver specification

       ! ---- lineminimization ---
       flag = f_selectBlock( tag_lineminimization) == 0
       if(.not.flag) flag = f_selectBlock( tag_line_minimization) == 0
       if(flag) then
          call m_CtrlP_rd_val(nfout, tag_dt_lower_critical, 'au_time' &
            &, dt_Lower_CRITICAL, rr)
          call m_CtrlP_rd_val(nfout, tag_dt_upper_critical, 'au_time' &
            &, dret, rr)
          if(dret >  dt_Lower_CRITICAL .and. dret > 1.d-3) then
             dt_Upper_CRITICAL = dret
          else
             if(printable) write(nfout &
                  & ,'(" !** the given value of dt_Upper_CRITICAL is invalid: dt_Upper_CRITICAL = ",d15.5)') dret
          end if
! --> T. Yamasaki, 15th July 2008
          call m_CtrlP_rd_val(nfout, tag_dt_upper_factor,'',dt_upper_factor &
          &    , rr, done_something)
          if(done_something) then
             if(dt_upper_factor < 1.0) dt_upper_factor = 1.0
             dt_upper_factor_is_set = .true.
          end if
          call m_CtrlP_rd_val(nfout, tag_dt_lower_factor,'',dt_lower_factor &
          &    , rr, done_something)
          if(done_something) then
             if(dt_lower_factor > 1.0) dt_lower_factor = 1.0
          end if
! <---

          call m_CtrlP_rd_val(nfout, tag_delta_lmdenom, '', delta_lmdenom,rr)

          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             write(nfout,'(" !** -- tag for lineminimization is found --")')
             write(nfout,'(" !** dt_lower_critical, dt_upper_critical = ",2f12.8)') dt_lower_critical, dt_upper_critical
             write(nfout,'(" !** delta_lmdenom = ",d20.8)') delta_lmdenom
          end if

! --> T. Yamasaki, 19th June 2009
          call m_CtrlP_rd_val(nfout, tag_energy_evaluation, LOWER, rstr, rr, done_something)
          if( done_something ) then
             call set_energy_evaluation(rstr) ! --> energy_evaluation
             if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
                write(nfout,'(" !** -- tag_energy_evaluation is found --")')
                write(nfout,'(" !** energy_evaluation = ",i3," : 1=TOTAL_ENERGY, 2=MODIFIED_TOTAL_ENERGY, 3=BAND_ENERGY")') &
                     & energy_evaluation
             end if
          else
             if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
                write(nfout,'(" !** -- tag_energy_evaluation is not found --")')
             end if
          end if
          if(energy_evaluation == BAND_ENERGY .or. energy_evaluation == MODIFIED_TOTAL_ENERGY) then
             call m_CtrlP_rd_val(nfout, tag_num_conduction_bands, num_conduction_bands_lmm, rr, done_something)
             if( done_something ) then
                if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
                   write(nfout,'(" !** -- tag_num_conduction_bands is found --")')
                   write(nfout,'(" !**  num_conduction_bands = ",i8)') num_conduction_bands_lmm
                end if
             end if
          end if
! <--
! --> T. Yamasaki, 18th Aug. 2009
          sw_lmm_status_check = -1
          call m_CtrlP_rd_val(nfout, tag_lmm_status_check, sw_lmm_status_check, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_lmm_status_check, sw_lmm_status_check, rr)
          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             if(sw_lmm_status_check /= -1) &
                  & write(nfout,'(" !** -- tag_sw_lmm_status_check is found --")')
             if(sw_lmm_status_check == -1) sw_lmm_status_check = NO
             write(nfout,'(" !**  sw_lmm_status_check = ",i3)') sw_lmm_status_check
          end if
! <--
          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
               & write(nfout,'(" !** -- tag for lineminimization is not found --")')
       end if

       ! ---- rmm ---
       if( explict_solver )then
          edelta_change_to_rmm = 1.e-3/dble(natm2)
          edelta_change_to_rmm_md = 1.e-3/dble(natm2)
       endif

       if( f_selectBlock( tag_rmm) == 0) then
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) write(nfout,'(" !** -- tag_rmm is found --")')
          call m_CtrlP_rd_val(nfout, tag_imGSrmm, imGSrmm, rr)
          call m_CtrlP_rd_val(nfout, tag_rr_critical_value, '', rr_Critical_Value, rr)
          call m_CtrlP_rd_val(nfout, tag_rmm_precal_phase_matm, &
          &    rmm_precal_phase_matm, rr)
          call m_CtrlP_rd_val(nfout, tag_edelta_change_to_rmm, 'hartree', &
          &    edelta_change_to_rmm, rr, done_something)
          if(done_something) then
             edelta_change_to_rmm_md = dret
             edelta_rmm_given = .true.
          endif
          call m_CtrlP_rd_val(nfout, tag_edelta_change_to_rmm_md, 'hartree', &
          &    edelta_change_to_rmm_md, rr, done_something)
          if(done_something) then
             edelta_rmm_given = .true.
          endif
          call m_CtrlP_rd_val(nfout, tag_save_memory_mode, rmm_save_memory_mode,rr)
          iret = f_selectParentBlock()
          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             write(nfout,'(" !** --- parameters for rmm ---")')
             write(nfout,'(" !** rr_Critical_Value = ",d20.8)') rr_Critical_Value
             write(nfout,'(" !** edelta_change_to_rmm = ",d20.8)') edelta_change_to_rmm
             write(nfout,'(" !** rmm_precal_phase_matm = ",i10)') rmm_precal_phase_matm
             write(nfout,'(" !** save_memory_mode = ",i10)') rmm_save_memory_mode
          end if
       else
          if(ipriinputfile >= 2 .and. printable .and. .not.rr) &
               & write(nfout,'(" !* tag_rmm is not found")')
       end if

       ! --- subspace_rotation ---
       flag = f_selectBlock( tag_submat) == 0
       if(.not.flag) flag = f_selectBlock( tag_subspacerotation) == 0
       if(.not.flag) flag = f_selectBlock( tag_subspace_rotation) == 0
       if( flag ) then
          meg = neg
          if(ipriinputfile >= 1 .and. printable .and. .not.rr) &
               & write(nfout,'(" !** tag_subspace_rotation is found")')
          call m_CtrlP_rd_val(nfout, tag_subspace_matrix_size, meg, rr)
          call m_CtrlP_rd_val(nfout, tag_damping_factor, '', damp, rr)
          call m_CtrlP_rd_val(nfout, tag_before_renewal, &
               &              submat_before_renewal,rr)
          if( meg > neg .or. meg < 1) then
             if(printable) then
                write(nfout,'(" !** given subspace_matrix_size(=meg) = ",i6)') meg
                write(nfout,'(" !*  subspace_matrix_size (meg) is set to be neg (= ",i6,")")') neg
             end if
             meg = neg
          end if
          if(damp > 1.0 .or. damp < 0.0) then
             if(printable) then
                write(nfout,'(" !** given damping factor of subspace rotation = ",f8.4)') damp
                write(nfout,'(" !*  damping factor (damp) is set to be 1.0")')
             end if
             damp = 1.d0
          end if

          call m_CtrlP_rd_val(nfout, tag_period, submat_period, rr, done_something)
          if(ipriinputfile >= 2 .and. done_something .and. printable .and. .not.rr) &
          & write(nfout,'(" !* tag_period is found")')

          if(.not.done_something) then
             call m_CtrlP_rd_val(nfout, tag_one_period, submat_period, rr,&
                  &              done_something)
             if(ipriinputfile >= 2 .and. done_something .and. printable .and. &
             &  .not.rr) write(nfout,'(" !* tag_one_period is found")')
          end if

          call m_CtrlP_rd_val(nfout, tag_critical_ratio, '', submat_critical_ratio &
          & ,  rr, done_something)
          if( done_something ) then
             if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
                  & write(nfout,'(" !* tag_critical_ratio is found")')
          end if

          if( f_getIntValue(tag_sw_gep, iret) == 0) sw_gep = iret
          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             write(nfout,'(" !** --- parameters for subspace roation ---")')
             write(nfout,'(" !** subspace_matrix_size = ",i6)') meg
             write(nfout,'(" !** damping_factor       = ",f8.4)') damp
             write(nfout,'(" !** submat_period        = ",i6)') submat_period
             write(nfout,'(" !** submat_critical_ratio= ",d20.8)') submat_critical_ratio
             if(submat_before_renewal == ON) then
                write(nfout,'(" !** before_renewal= ON")')
             else
                write(nfout,'(" !** before_renewal= OFF")')
             end if
             write(nfout,'(" !** sw_gep = ",i6)') sw_gep
          end if

          if (.not.rr) then
            if( f_selectBlock( tag_scalapack) == 0) then
               if( f_getIntValue(tag_sw_scalapack, iret) == 0) sw_scalapack = iret
#ifndef _USE_SCALAPACK_
               if(sw_scalapack == ON) then
                  write(nfout,*) 'ScaLAPACK is unusable.'
                  write(nfout,*) 'Use the CPP option -D_USE_SCALAPACK_.'
                  call phase_error_with_msg(nfout,'ScaLAPACK is unusable.',__LINE__,__FILE__)
               end if
#endif
               call set_method_scalapack(method_scalapack)
               if( f_getIntValue(tag_block_size, iret) == 0) block_size = iret
               if( f_getIntValue(tag_nprow, iret) == 0) nprow = iret
               if( f_getIntValue(tag_npcol, iret) == 0) npcol = iret
               if( f_getIntValue(tag_divide_square, iret) == 0) divide_square = iret
               if(ipriinputfile >= 1 .and. printable) then
                  write(nfout,'(" !** tag_scalapack is found")')
                  write(nfout,'(" !** sw_scalapack = ",i6)') sw_scalapack
                  write(nfout,'(" !** method = ",i6)') method_scalapack
                  write(nfout,'(" !** block_size = ",i6)') block_size
                  write(nfout,'(" !** nprow = ",i6)') nprow
                  write(nfout,'(" !** npcol = ",i6)') npcol
                  write(nfout,'(" !** divide_square = ",i6)') divide_square
               end if
               iret = f_selectParentBlock()
            end if
            iret = f_selectParentBlock()
          endif
       else
          if(.not. rr) meg = neg
       end if

       ! ---- Davidson ---
       if(ekmode == ON .or. &
         & (icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION)) then
          sw_first_conv_check = OFF
       end if
       if( f_selectBlock( tag_davidson) == 0) then
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
          &  write(nfout,'(" !** -- tag_davidson is found --")')
          call m_CtrlP_rd_val(nfout, tag_ndavid, ndavid, rr)
          call m_CtrlP_rd_val(nfout, tag_max_iter_david, max_iter_david, rr)
          call m_CtrlP_rd_val(nfout, tag_max_subspace_size, &
          &    max_subspace_size, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_first_conv_check, &
          &    sw_first_conv_check, rr)
          if(ndavid<0) ndavid=abs(ndavid)
          if(max_iter_david<0) max_iter_david=abs(max_iter_david)
          if( max_subspace_size < 4*neg) max_subspace_size = 4*neg ! default value

          call m_CtrlP_rd_val(nfout, tag_delta_eig_occup, 'hartree', &
          &    delta_eig_occup, rr)
          call m_CtrlP_rd_val(nfout, tag_delta_eig_empty, 'hartree', &
          &    delta_eig_empty, rr)
          call m_CtrlP_rd_val(nfout,tag_eps_david, '', eps_david, rr)
!! md-david, md-kosugi
          call m_CtrlP_rd_val(nfout, tag_submat_GE, submat_GE, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_MRCV_only, sw_MRCV_only, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_divide_subspace, &
          &    sw_divide_subspace, rr, done_something)
          if( done_something ) then
              sw_divide_subspace_changed = .false.
          endif
          call m_CtrlP_rd_val(nfout, tag_npartition, iret, rr,&
          &    done_something)
          if(.not. done_something) then
            call m_CtrlP_rd_val(nfout, tag_npartition_david, iret, rr,&
            &    done_something)
          endif
          if( done_something ) then
              npartition_david = iret
              sw_npartition_changed = .false.
          endif
          call m_CtrlP_rd_val(nfout, tag_nbands_in_partition, iret, rr, done_something)
          if(.not. done_something) then
            call m_CtrlP_rd_val(nfout, tag_nblock_david, iret, rr, done_something)
          endif
          if(.not. done_something) then
            call m_CtrlP_rd_val(nfout, tag_nblock, iret, rr, done_something)
          endif
          if(done_something) then
             if(iret>0) then
                npartition_david = neg/(iret*nrank_e)
                if (npartition_david<1) npartition_david = 1
                sw_npartition_changed = .false.
             endif
          endif

          if(ipriinputfile >= 2 .and. printable .and. .not. rr) then
             write(nfout,'(" !** ndavid = ",i10)') ndavid
             write(nfout,'(" !** max_iter_david = ",i10)') max_iter_david
             write(nfout,'(" !** max_subspace_size = ",i10)') max_subspace_size
             write(nfout,'(" !** delta_eig_occup = ",e12.5)') delta_eig_occup
             write(nfout,'(" !** delta_eig_empty = ",e12.5)') delta_eig_empty
             write(nfout,'(" !** eps_david = ",e12.5)') eps_david
             write(nfout,'(" !** sw_first_conv_check = ",i10)') sw_first_conv_check
             write(nfout,'(" !** sw_MRCV_only = ",i10)') sw_MRCV_only
             write(nfout,'(" !** submat_GE = ",i10)') submat_GE
             write(nfout,'(" !** sw_divide_subspace = ",i10)') sw_divide_subspace
             write(nfout,'(" !** npartition = ",i10)') npartition_david
          end if
          if( .not. rr) then
            if( f_selectBlock( tag_dav_scalapack) == 0) then
               if( f_getIntValue(tag_sw_dav_scalapack, iret) == 0) sw_dav_scalapack = iret
#ifndef _USE_SCALAPACK_
               if(sw_dav_scalapack == ON) then
                  write(nfout,*) 'ScaLAPACK is unusable.'
                  write(nfout,*) 'Use the CPP option -D_USE_SCALAPACK_.'
                  call phase_error_with_msg(nfout,'ScaLAPACK is unusable.',__LINE__,__FILE__)
               end if
#endif
               if( f_getIntValue(tag_dav_block_size, iret) == 0) dav_block_size = iret
               if( f_getIntValue(tag_dav_nprow, iret) == 0) dav_nprow = iret
               if( f_getIntValue(tag_dav_npcol, iret) == 0) dav_npcol = iret
               if( f_getIntValue(tag_dav_divide_square, iret) == 0) dav_divide_square = iret
               if(ipriinputfile >= 1 .and. printable) then
                  write(nfout,'(" !** tag_dav_scalapack is found")')
                  write(nfout,'(" !** sw_dav_scalapack = ",i6)') sw_dav_scalapack
                  write(nfout,'(" !** dav_block_size = ",i6)') dav_block_size
                  write(nfout,'(" !** dav_nprow = ",i6)') dav_nprow
                  write(nfout,'(" !** dav_npcol = ",i6)') dav_npcol
                  write(nfout,'(" !** dav_divide_square = ",i6)') dav_divide_square
               end if
               iret = f_selectParentBlock()
            end if
          endif
          iret = f_selectParentBlock()
       else
          max_subspace_size = 4*neg ! default value
!!$          if(ipriinputfile >= 2 .and. printable) &
          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             write(nfout,'(" !* tag_davidson is not found")')
             write(nfout,'(" !** max_subspace_size = ",i6)') max_subspace_size
          end if
       end if

       ! ---- Modified Davidson ---
       if( f_selectBlock( tag_mddavidson) == 0) then
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) write(nfout,'(" !** -- tag_mddavidson is found --")')
          call m_CtrlP_rd_val(nfout, tag_npartition_mddavid, npartition_mddavid, rr)
          call m_CtrlP_rd_val(nfout, tag_max_iter_mddavid, max_iter_mddavid, rr)
          if(npartition_mddavid<0) npartition_mddavid=abs(npartition_mddavid)
          if(max_iter_mddavid<0) max_iter_mddavid=abs(max_iter_mddavid)
          call m_CtrlP_rd_val(nfout, tag_delta_eig_occup, 'hartree', &
          &    delta_eig_occup_md, rr)
          call m_CtrlP_rd_val(nfout, tag_delta_eig_empty, 'hartree', &
          &    delta_eig_empty_md, rr)
          call m_CtrlP_rd_val(nfout,tag_eps_mddavid, '', eps_mddavid, rr)
          call m_CtrlP_rd_val(nfout,tag_eps_residual_mddavid, '', eps_residual_mddavid, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_apply_gs, sw_apply_gs, rr)

          max_iter_david = max_iter_mddavid
          delta_eig_occup = delta_eig_occup_md
          delta_eig_empty = delta_eig_empty_md
          eps_david = eps_mddavid
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) then
             write(nfout,'(" !** npartition_mddavid = ",i10)') npartition_mddavid
             write(nfout,'(" !** max_iter_mddavid = ",i10)') max_iter_mddavid
             write(nfout,'(" !** delta_eig_occup = ",e12.5)') delta_eig_occup_md
             write(nfout,'(" !** delta_eig_empty = ",e12.5)') delta_eig_empty_md
             write(nfout,'(" !** eps_mddavid = ",e12.5)') eps_mddavid
             write(nfout,'(" !** eps_residual_mddavid = ",e12.5)') eps_residual_mddavid
          end if
       else
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
               & write(nfout,'(" !* tag_mddavidson is not found")')
       end if

       ! ---- Modified Kosugi ---
       if( f_selectBlock( tag_mdkosugi) == 0) then
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) write(nfout,'(" !** -- tag_mdkosugi is found --")')
          call m_CtrlP_rd_val(nfout, tag_npartition_mdkosugi, npartition_mdkosugi, rr)
          call m_CtrlP_rd_val(nfout, tag_max_iter_mdkosugi, max_iter_mdkosugi, rr)
          if(npartition_mdkosugi<0) npartition_mdkosugi=abs(npartition_mdkosugi)
          if(max_iter_mdkosugi<0) max_iter_mdkosugi=abs(max_iter_mdkosugi)
          call m_CtrlP_rd_val(nfout, tag_delta_eig_occup_mdkosugi, 'hartree', &
          &    delta_eig_occup_mdkosugi, rr)
          call m_CtrlP_rd_val(nfout, tag_delta_eig_empty_mdkosugi, 'hartree', &
          &    delta_eig_empty_mdkosugi, rr)
          call m_CtrlP_rd_val(nfout,tag_eps_mdkosugi, '', eps_mddavid, rr)
          call m_CtrlP_rd_val(nfout,tag_eps_residual_mdkosugi, '', eps_residual_mddavid, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_apply_gs, sw_apply_gs, rr)


          max_iter_david = max_iter_mdkosugi
          delta_eig_occup = delta_eig_occup_mdkosugi
          delta_eig_empty = delta_eig_empty_mdkosugi
          eps_david = eps_mdkosugi
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) then
             write(nfout,'(" !** npartition_mdkosugi = ",i10)') npartition_mdkosugi
             write(nfout,'(" !** max_iter_mdkosugi = ",i10)') max_iter_mdkosugi
             write(nfout,'(" !** delta_eig_occup_mdkosugi = ",e12.5)') delta_eig_occup_mdkosugi
             write(nfout,'(" !** delta_eig_empty_mdkosugi = ",e12.5)') delta_eig_empty_mdkosugi
             write(nfout,'(" !** eps_mdkosugi = ",e12.5)') eps_mdkosugi
             write(nfout,'(" !** eps_residual_mdkosugi = ",e12.5)') eps_residual_mdkosugi
          end if
       else
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
               & write(nfout,'(" !* tag_mdkosugi is not found")')
       end if
       ! ---- CG  ---
       if( f_selectBlock( tag_cg) == 0) then
         call m_CtrlP_rd_val(nfout, tag_sw_modified_cg_formula, sw_modified_cg_formula, rr)
         if(ipriinputfile>=2 .and. printable .and. .not. rr) &
         & write(nfout,'(a,i3)') 'modified formula for CG : ',sw_modified_cg_formula
       endif

       iret = f_selectParentBlock()
    else
       tag_solver_of_WF_is_found = .true.
       if(ipriinputfile >=2 .and. printable .and. .not. rr) &
            & write(nfout,'(" !** -- tag_wavefunction_solver is not found --")')
       meg = neg
    end if
    if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
       write(nfout,'(" !** --- id, sol, till_n,  dts, dte, itr, var, prec, cmix, submat ---")')
       do i = 1, n_WF_solvers_all
          if(i == 1) write(nfout,'(" !**      - for_init_str -")')
          if(i == n_WF_solvers_before+1) write(nfout,'(" !**      - during_str_relax -")')
          write(nfout,'(" !** ",i3,i4,i6,f7.3,f7.3,i4,i3,i3,i3,i3)') i, w_solver(i)%solver, w_solver(i)%till_n_iter &
               &              , w_solver(i)%dtim_s, w_solver(i)%dtim_e, w_solver(i)%iter_range &
               &              , w_solver(i)%variation_way, w_solver(i)%precon, w_solver(i)%cmix_pointer &
               &              , w_solver(i)%subspace_rotation
       end do
    end if
  contains
    subroutine set_energy_evaluation(rstr)
      character(len=FMAXVALLEN), intent(in) :: rstr
      logical :: tf
      energy_evaluation = TOTAL_ENERGY
      call strncmp0(tag_total_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = TOTAL_ENERGY
         goto 1001
      end if
      call strncmp0(tag_modified_total_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = MODIFIED_TOTAL_ENERGY
         goto 1001
      end if
      call strncmp0(tag_band_energy,trim(rstr),tf)
      if(tf) then
         energy_evaluation = BAND_ENERGY
         goto 1001
      end if
      call phase_error_with_msg(nfout,' ! tag for energy_evaluation is invalid <<m_CtrlP_rd_wfsolver.set_energy_evaluation>>'&
                               ,__LINE__,__FILE__)
1001  continue
    end subroutine set_energy_evaluation

    subroutine configure_wf_solver(solver_set,natm2)
       integer, intent(in) :: solver_set,natm2
       integer :: i
       if(solver_set == LMM_RMM)then
          call alloc_w_solver(2)
          w_solver(1)%solver = lmMSD
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = 5
          w_solver(2)%solver = RMM3
          w_solver(2)%till_n_iter = -1
          w_solver(2)%subspace_rotation = ON
          edelta_change_to_rmm = 1.d-4/dble(natm2)
          edelta_change_to_rmm_md = 1.d-4/dble(natm2)
          n_WF_solvers_before = 2
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : lm+msd -> rmm3")')
       else if (solver_set == DAV_RMM)then
          call alloc_w_solver(4)
          if(icond==INITIAL .or. icond==CONTINUATION .or. icond==AUTOMATIC)then
             w_solver(1)%solver = MDDAVIDSON
             w_solver(3)%solver = MDDAVIDSON
          else
             w_solver(1)%solver = MDKOSUGI
             w_solver(3)%solver = MDKOSUGI
          endif
          if(sw_hubbard==ON) then
             w_solver(1)%solver = MDKOSUGI
             w_solver(3)%solver = MDKOSUGI
          endif
! === KT_add === 2015/01/05
          if ( noncol ) then
             w_solver(1)%solver = MDDAVIDSON
             w_solver(3)%solver = MDDAVIDSON
          endif
! ============== 2015/01/05

!          w_solver(1)%solver = MDKOSUGI
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = 5
          w_solver(1)%precon = ON
          w_solver(1)%before_or_after_convergence = BEFORE
          w_solver(2)%solver = RMM3
          w_solver(2)%precon = ON
          w_solver(2)%till_n_iter = -1
          w_solver(2)%subspace_rotation = ON
          w_solver(2)%before_or_after_convergence = BEFORE
          !!w_solver(3)%solver = MDDAVIDSON
          w_solver(3)%subspace_rotation = ON
          w_solver(3)%till_n_iter = 5
          w_solver(3)%precon = ON
          w_solver(3)%before_or_after_convergence = AFTER
          w_solver(4)%solver = RMM3
          w_solver(4)%precon = ON
          w_solver(4)%till_n_iter = -1
          w_solver(4)%subspace_rotation = ON
          w_solver(4)%before_or_after_convergence = AFTER
          edelta_change_to_rmm = 1.d-3/dble(natm2)
          edelta_change_to_rmm_md = 1.d-3/dble(natm2)
          n_WF_solvers_before = 2
          n_WF_solvers_after = 2
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : davidson -> rmm3")')
       else if (solver_set == lmMSD) then
          call alloc_w_solver(1)
          w_solver(1)%solver = lmMSD
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : lm+msd")')
       else if (solver_set == DAVIDSON)then
          call alloc_w_solver(1)
          w_solver(1)%solver = DAVIDSON
          w_solver(1)%subspace_rotation = ON
          w_solver(1)%precon = ON
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if(neg/nrank_e<4) then
             sw_divide_subspace=OFF
             sw_divide_subspace_changed = .true.
             write(nfout,'(" !** REMARK: sw_divide_subspace is set to OFF ")')
          endif
          !!$if (sw_hubbard==ON.or.nspin>1) sw_divide_subspace=OFF
          if (printable) write(nfout,'(" !** applied wavefunction solver set : davidson")')
       else if (solver_set == FEF)then
          call alloc_w_solver(1)
          w_solver(1)%solver = CG
          w_solver(1)%subspace_rotation = OFF
          w_solver(1)%till_n_iter = -1
          n_WF_solvers_before = 1
          n_WF_solvers_after = 0
          n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
          if (printable) write(nfout,'(" !** applied wavefunction solver set : FEF (cg)")')
       endif
    end subroutine configure_wf_solver

    subroutine set_wfsolvers(prealloc,msol,nbase,ba,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: msol,nbase,ba
      integer, intent(out):: iret

      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: i, sol, till_n, itr, var, prec,cmix, no, sw_submat
      integer, dimension(9) :: icf
      real(kind=DP) :: dts,dte
      if( f_selectBlock(tag_solvers) == 0) then
         i = 1
         do while(.true.)
            if(i == 1) then
               if( f_selectFirstTableLine() /= 0 ) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0 ) then
                  exit
               end if
            end if
            if(.not.prealloc) then
               if(i > msol) exit
               no = i
               iret = f_readsolver(icf,no,sol,till_n,dts,dte,itr,var,prec,cmix,sw_submat)

               if(no <= msol) then
                  w_solver(no+nbase)%before_or_after_convergence = ba
                  if(icf(1) == 1 .and. w_solver(no+nbase)%solver /= sol) &
                  &  w_solver(no+nbase)%solver = sol
                  if(icf(2) == 1 .and. w_solver(no+nbase)%till_n_iter /= till_n) &
                  &  w_solver(no+nbase)%till_n_iter = till_n
                  if(icf(3) == 1 .and. w_solver(no+nbase)%dtim_s /= dts) &
                  &  w_solver(no+nbase)%dtim_s = dts
                  if(icf(4) == 1 .and. w_solver(no+nbase)%dtim_e /= dte) then
                     w_solver(no+nbase)%dtim_e = dte
                  else if(icf(3) == 1 .and. w_solver(no+nbase)%dtim_e &
                  &              /= w_solver(no+nbase)%dtim_s) then
                     w_solver(no+nbase)%dtim_e = w_solver(no+nbase)%dtim_s
                  end if
                  if(icf(5) == 1 .and. w_solver(no+nbase)%iter_range /= itr)    &
                  &  w_solver(no+nbase)%iter_range = itr
                  if(icf(6) == 1 .and. w_solver(no+nbase)%variation_way /= var) &
                  &  w_solver(no+nbase)%variation_way = var
                  if(icf(7) == 1 .and. w_solver(no+nbase)%precon /= prec)       &
                  &  w_solver(no+nbase)%precon = prec
                  if(icf(8) == 1 .and. w_solver(no+nbase)%cmix_pointer /= cmix) &
                  &  w_solver(no+nbase)%cmix_pointer = cmix
                  if(icf(9) == 1 .and. w_solver(no+nbase)%subspace_rotation /= sw_submat) &
                  &  w_solver(no+nbase)%subspace_rotation = sw_submat
               end if
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      else
         call phase_error_with_msg(nfout,' ! No wf solver is given in the inputfile <<m_CtrlP_rd_wfsolver.set_wfsolvers>>'&
                                  ,__LINE__,__FILE__)
      end if
        iret = i-1
    end subroutine set_wfsolvers

  end subroutine m_CtrlP_rd_wfsolver2

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_dealloc_wfsolver()
    call dealloc_w_solver()
  end subroutine m_CtrlP_dealloc_wfsolver
#endif

  subroutine m_CtrlP_set_sw_MRCV_only(sw)
    integer, intent(in) :: sw
    sw_MRCV_only = sw
    if(sw_MRCV_only /= ON .and. sw_MRCV_only /= OFF) sw_MRCV_only = OFF
  end subroutine m_CtrlP_set_sw_MRCV_only

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intelligent reading of Wave-function solvers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function f_readsolver(icf,no,sol,till_n,dts,dte,itr,var,prec,cmix,sw_submat)
!  Coded by T. Yamasaki (FUJITSU LABORATORIES LTD), 13th Jun 2003
!  A tag "subspace_rotation" is newly introduced by T. Yamasaki, 27th Aug. 2003
    integer, intent(out), dimension(9) :: icf
    integer, intent(inout) :: no
    integer, intent(out) ::  sol,till_n,itr,var,prec,cmix,sw_submat
    real(kind=DP),intent(inout) :: dts,dte
    integer :: iret
!!$    character(len=FMAXVALLEN) :: rstr
    integer :: f_getstringValue, f_getIntValue, f_getRealValue
    logical :: tf

    icf = 0
    if(f_getIntValue(tag_no, iret) == 0) then
       no = iret
    else if(f_getIntValue(tag_id, iret) == 0) then
       no = iret
    end if

    if( f_getIntValue(tag_till_n,till_n) == 0) icf(2) = 1     ! till_n
    if( f_getRealValue(tag_dts,dts,'au_time') == 0) icf(3) = 1   ! dts
    if( f_getRealValue(tag_dte,dte,'au_time') == 0) icf(4) = 1   ! dte
    if( f_getIntValue(tag_itr,itr) == 0) icf(5) = 1           ! iret

    ! --- var ---
    if( f_getStringValue(tag_var,rstr,LOWER) == 0) then
       call strncmp0(tag_linear,trim(rstr),tf)
       if(tf) then
          icf(6) = 1
          var = varLINEAR
       else
          call strncmp0(tag_tanh,trim(rstr),tf)
          if(tf) then
             icf(6) = 1
             var = varTANH
          end if
       end if
    end if

    if( f_getIntValue(tag_prec,prec) == 0) icf(7) = 1         ! prec
    if( f_getIntValue(tag_cmix,cmix) == 0) icf(8) = 1         ! cmix

    if( f_getIntValue(tag_submat,sw_submat) == 0) icf(9) = 1     ! submat
    if(icf(9) /= 1) then
       if( f_getIntValue(tag_subspacerotation,sw_submat) == 0) icf(9) = 1
    end if
    if(icf(9) /= 1) then
       if( f_getIntValue(tag_subspace_rotation,sw_submat) == 0) icf(9) = 1
    end if

    ! --- sol --
    sol = MSD
    f_readsolver = f_getstringValue(tag_sol,rstr,LOWER)
    if( f_readsolver == 0 ) then
       icf(1) = 1
       if(trim(rstr) == trim(tag_lmMSD)) then
          sol = lmMSD
       else if(trim(rstr) == trim(tag_lmmsd2)) then
          sol = lmMSD
       else if(trim(rstr) == trim(tag_msd)) then
          sol = MSD
       else if(trim(rstr) == trim(tag_lmsd)) then
          sol = lmSD
       else if(trim(rstr) == trim(tag_lmsd2)) then
          sol = lmSD
       else if(trim(rstr) == trim(tag_sd)) then
          sol = SD
       else if(trim(rstr) == trim(tag_rmm2p)) then
          sol = RMM2p
#ifdef MEMORY_SAVE_ZAJ_OLD
          RMM2P_is_specified = .true.
#endif
       else if(trim(rstr) == trim(tag_rmm3)) then
          sol = RMM3
       else if(trim(rstr) == trim(tag_rmm2)) then
          sol = RMM2
       else if(trim(rstr) == trim(tag_matrixdiagon)) then
          sol = MATRIXDIAGON
       else if(trim(rstr) == trim(tag_matdiagon)) then
          sol = MATRIXDIAGON
       else if(trim(rstr) == trim(tag_submat)) then
          sol = SUBMAT
       else if(trim(rstr) == trim(tag_cg)) then
          sol = CG
       else if(trim(rstr) == trim(tag_lmcg)) then
          sol = CG
       else if(trim(rstr) == trim(tag_lmcg2)) then
          sol = CG
       else if(trim(rstr) == trim(tag_davidson)) then
          sol = DAVIDSON
       else if(trim(rstr) == trim(tag_mddavidson)) then
          sol = MDDAVIDSON
       else if(trim(rstr) == trim(tag_pkos1)) then
          sol = MDDAVIDSON
       else if(trim(rstr) == trim(tag_pkos2)) then
          sol = MDDAVIDSON
       else if(trim(rstr) == trim(tag_pkos3)) then
          sol = MDDAVIDSON
       else if(trim(rstr) == trim(tag_mdkosugi)) then
          sol = MDKOSUGI
       else if(trim(rstr) == trim(tag_pdav1)) then
          sol = MDKOSUGI
       else if(trim(rstr) == trim(tag_pdav2)) then
          sol = MDKOSUGI
       else if(trim(rstr) == trim(tag_pdav3)) then
          sol = MDKOSUGI
       else
          call phase_error_with_msg(6,' ! sol is not given properly <<f_readsolver>>',__LINE__,__FILE__)
       end if
    end if
  end function f_readsolver
#endif

  subroutine m_CtrlP_rd_struc_evol(nfout)
    integer, intent(in) :: nfout
    integer :: iret
    logical :: tf
    real(kind=DP) :: dret
!!$    character(len=FMAXVALLEN) :: rstr

    integer :: f_selectTop, f_getIntValue, f_selectParentBlock, f_selectBlock
    integer :: f_getRealValue, f_getStringValue

    iret = f_selectTop()
    if( f_selectBlock( tag_structure_evolution) == 0) then
       if( f_getIntValue(tag_sw_fix,iret)== 0) then
          sw_fix = iret
       end if

       if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
          call strncmp2(trim(rstr),len(trim(rstr)) &
               & ,tag_quench_with_constraint,len(tag_quench_with_constraint),tf)
          if(tf) then
             imdalg = QUENCHED_CONSTRAINT
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_quench,tf)
          if(tf) then
             imdalg = QUENCHED_MD
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_gdiis,tf)
          if(tf) then
             imdalg = GDIIS
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_velocity_verlet,tf)
          if(.not.tf) call strncmp0(trim(rstr),tag_verlet,tf)
          if(tf) then
             imdalg = VERLET
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_temperature_control,tf)
          if(.not.tf) call strncmp0(trim(rstr),tag_t_control,tf)
          if(.not.tf) call strncmp0(trim(rstr),tag_temp_control,tf)
          if(tf) then
             imdalg = T_CONTROL
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_bluemoon,tf)
          if(tf) then
             imdalg = BLUEMOON
             goto 1001
          end if
! === DEBUG by tkato 2014/03/20 ================================================
!         call strncmp2(trim(rstr),FMAXVALLEN,tag_cg2,len(tag_cg2),tf)
          call strncmp2(trim(rstr),len(trim(rstr)),tag_cg2,len(tag_cg2),tf)
! ==============================================================================
          if(tf) then
             imdalg = CG_STROPT2
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_cg,tf)
          if(tf) then
             imdalg = CG_STROPT2
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_oldcg,tf)
          if(tf) then
             imdalg = CG_STROPT
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_fire,tf)
          if(tf) then
             imdalg = FIRE
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_sd,tf)
          if(.not.tf) call strncmp0(trim(rstr),tag_steepest_descent,tf)
          if(tf) then
             imdalg = STEEPEST_DESCENT
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_damp,tf)
          if(tf) then
             imdalg = DAMPED_MD
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_velocity_scaling,tf)
          if(tf) then
             imdalg = VELOCITY_SCALING
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_temperature_pressure_control,tf)
          if(tf) then
             imdalg = PT_CONTROL
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_pressure_temperature_control,tf)
          if(tf) then
             imdalg = PT_CONTROL
             goto 1001
          end if
          call strncmp0(trim(rstr),tag_pressure_control,tf)
          if(tf) then
             imdalg = P_CONTROL
             goto 1001
          end if
! === Include PAW by tkato =====================================================
          call strncmp0(trim(rstr),tag_bfgs,tf)
          if(tf) then
             imdalg = BFGS
             goto 1001
          end if
! ==============================================================================
          call strncmp0(trim(rstr),tag_lbfgs,tf)
          if(tf) then
             imdalg = L_BFGS
             goto 1001
          end if
1001      continue
       end if

       if(imdalg==VELOCITY_SCALING .or. imdalg==T_CONTROL .or. imdalg==VERLET .or. &
      &   imdalg == PT_CONTROL .or. imdalg == P_CONTROL) iprivelocity=2
       if(imdalg == PT_CONTROL .or. imdalg == P_CONTROL) istress = 1

       if( f_getRealValue(tag_dt,dret,'au_time') == 0) dtio = dret
       if(ipriinputfile >= 1 .and. printable) then
          write(nfout,'(" !** --- structure_evolution ---")')
          if(sw_fix == ON) then
             write(nfout,'(" !** sw_fix = ON")')
          end if
          mdalgorithm: select case(imdalg)
          case (QUENCHED_CONSTRAINT)
             write(nfout,'(" !** imdalg = ",i5," (= QUENCHED_CONSTRAINT)")') imdalg
          case (QUENCHED_MD)
             write(nfout,'(" !** imdalg = ",i5," (= QUENCHED_MD)")') imdalg
          case (GDIIS)
             write(nfout,'(" !** imdalg = ",i5," (= GDIIS)")') imdalg
          case (VERLET)
             write(nfout,'(" !** imdalg = ",i5," (= VERLET)")') imdalg
          case (T_CONTROL)
             write(nfout,'(" !** imdalg = ",i5," (= T_CONTROL)")') imdalg
          case (BLUEMOON)
             write(nfout,'(" !** imdalg = ",i5," (= BLUEMOON)")') imdalg
          case (CG_STROPT)
             write(nfout,'(" !** imdalg = ",i5," (= CG_STROPT)")') imdalg
          case (CG_STROPT2)
             write(nfout,'(" !** imdalg = ",i5," (= CG_STROPT2)")') imdalg
          case (STEEPEST_DESCENT)
             write(nfout,'(" !** imdalg = ",i5," (= STEEPEST_DESCENT)")') imdalg
! === Include PAW by tkato =====================================================
          case(BFGS)
             write(nfout,'(" !** imdalg = ",i5," (= BFGS)")') imdalg
! ==============================================================================
          case(L_BFGS)
             write(nfout,'(" !** imdalg = ",i5," (= L_BFGS)")') imdalg
!!$ASASASASAS
          case (DAMPED_MD)
             write(nfout,'(" !** imdalg = ",i5," (= DAMPED_MD)")') imdalg
          case (VELOCITY_SCALING)
             write(nfout,'(" !** imdalg = ",i5," (= VELOCITY_SCALING)")') imdalg
!!$ASASASASAS
          case (PT_CONTROL)
             write(nfout,'(" !** imdalg = ",i5," (= PT_CONTROL)")') imdalg
          case (P_CONTROL)
             write(nfout,'(" !** imdalg = ",i5," (= P_CONTROL)")') imdalg
          case default
             write(nfout,'(" !** imdalg = ",i5," ( is not proper)")') imdalg
          end select mdalgorithm
          write(nfout,'(" !** dtio = ",f10.4)') dtio
       end if

       if(f_getIntValue(tag_sw_prec,iret)==0) then
          sw_prec = iret
       endif
       if(printable) write(nfout,'(a,i3)') ' !** sw_prec = ',sw_prec

!!$       if( f_getIntValue( tag_sw_force_file,iret) == 0) sw_force_file = iret
!!$       if(ipriinputfile >= 1) write(nfout,'(" !** sw_force_file = ",i5)') sw_force_file

       if( f_selectBlock( tag_sd) == 0) then
          tf = f_getIntValue(tag_mode_coefficient,iret)== 0
          if(tf) mode_fi_coefficient = iret
          if(mode_fi_coefficient == ON) then
             tf = f_getRealValue(tag_coefficient,dret,'') == 0
             if(tf) fi_coefficient = dret
          end if
          iret = f_selectParentBlock()
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** mode_fi_coefficient = ",i5)') mode_fi_coefficient
             write(nfout,'(" !**      fi_coefficient = ",f8.4)') fi_coefficient
          end if
       end if

       if( f_selectBlock( tag_gdiis) == 0 .or. f_selectBlock(tag_bfgs)==0 .or. f_selectBlock(tag_lbfgs)==0) then
          if( f_getStringValue( tag_initial_method, rstr, LOWER) == 0) then
             call strncmp0(trim(rstr),tag_quench,tf)
             if(tf) then
                initial_method_of_gdiis = QUENCHED_MD
                goto 1003
             end if
! === DEBUG by T.Kato 2013/07/17 ===============================================
!            call strncmp2(trim(rstr),FMAXVALLEN,tag_cg2,len(tag_cg2),tf)
             call strncmp2(trim(rstr),len(trim(rstr)),tag_cg2,len(tag_cg2),tf)
! ==============================================================================
             if(tf) then
                initial_method_of_gdiis = CG_STROPT2
                goto 1003
             endif
             call strncmp0(trim(rstr),tag_cg,tf)
             if(tf) then
!!$                initial_method_of_gdiis = CG_STROPT
                initial_method_of_gdiis = CG_STROPT2
                goto 1003
             endif
             call strncmp0(trim(rstr),tag_oldcg,tf)
             if(tf) then
                initial_method_of_gdiis = CG_STROPT
                goto 1003
             endif
             call strncmp0(trim(rstr),tag_sd,tf)
             if(.not.tf) call strncmp0(trim(rstr),tag_steepest_descent,tf)
             if(tf) then
                initial_method_of_gdiis = STEEPEST_DESCENT
                goto 1003
             end if
1003         continue
          end if
          if(ipriinputfile >= 2) write(nfout,'(" !** gdiis_box_size = ",i8)') kqnmditer_p
          if( f_getIntValue( tag_gdiis_box_size, iret) == 0) kqnmditer_p = iret
          if( f_getStringValue( tag_gdiis_update,rstr,LOWER) == 0) then
             call strncmp0(trim(rstr),tag_renew,tf)
             if(tf) then
                gdiis_hownew = RENEW
                goto 1002
             end if
             call strncmp0(trim(rstr),tag_anew,tf)
             if(tf) gdiis_hownew = ANEW
1002         continue
          end if
          if( f_getRealValue(tag_c_forc2gdiis,dret,'hartree/bohr') == 0) c_forc2gdiis = dret
          if( f_getIntValue(tag_c_iteration2GDIIS,iret) == 0) c_iteration2GDIIS = iret
! === Include PAW by tkato =====================================================
          if( f_getIntValue(tag_sw_correct_eigenvalue,iret)==0) sw_correct_eigenvalue=iret
! ==============================================================================
          if( f_getIntValue(tag_sw_optmize_alpha,iret)==0) sw_optimize_alpha = iret
          if( f_getRealValue(tag_maxstep,dret,'')==0) maxstep=dret
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2) then
             write(nfout,'(" !** gdiis_box_size = ",i8)') kqnmditer_p
             write(nfout,'(" !** gdiis_hownew   = ",i8)') gdiis_hownew
             write(nfout,'(" !** c_forc2gdiis   = ",f12.6)') c_forc2gdiis
             write(nfout,'(" !** initial_method = ",i8)') initial_method_of_gdiis
             write(nfout,'(" !** c_iteration2GDIIS = ",i8)') c_iteration2GDIIS
          end if
       end if

       if(f_selectBlock(tag_fire) == 0)then
          if(f_getRealValue(tag_fire_incre_factor,dret,'')==0) fire_incre_factor = dret
          if(f_getRealValue(tag_fire_decre_factor,dret,'')==0) fire_decre_factor = dret
          if(f_getRealValue(tag_fire_decre_factor_alpha,dret,'')==0) &
          &  fire_decre_factor_alpha = dret
          if(f_getRealValue(tag_fire_alpha_start,dret,'')==0) fire_alpha_start = dret
          if(f_getIntValue(tag_fire_nmin,iret)==0) fire_nmin = iret
          if(f_getRealValue(tag_fire_dtmax,dret,'au_time')==0) fire_dtmax = dret
          if(f_getRealValue(tag_fire_initial_dt,dret,'au_time')==0) fire_initial_dt = dret
          if(f_getRealValue(tag_fire_invmass_factor,dret,'')==0) fire_invmass_factor = dret
          iret = f_selectParentBlock()
       endif

       if  (f_selectBlock(tag_prec)==0)then
         if(f_getRealValue(tag_mu,dret,'hartree')==0)then
           precon_mu = dret
           if(printable) write(nfout,'(a,f15.10)') ' !** value of mu for the preconditioner : ',precon_mu
         endif
         if(f_getRealValue(tag_a,dret,'')==0) precon_A = dret
         iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_stress) == 0) then
          if( f_getIntValue( tag_sw_stress,iret) == 0) istress = iret
          if( f_getIntValue( tag_iconstpw,iret)  == 0) iconstpw = iret
          if( f_getIntValue( tag_sw_smear_ke,iret) == 0) sw_smear_ke = iret
          if( sw_smear_ke == ON ) then
             ke_sigma = 0.5d0
             ke_e0 = 0.5d0*(gmax*gmax-1.d0)
             ke_a = 0.5d0*(gmax*gmax*0.375d0)
             if( f_getRealValue( tag_e0,dret,'hartree') == 0 ) ke_e0 = dret
             if( f_getRealValue( tag_sigma,dret,'hartree') == 0 ) ke_sigma = dret
             if( f_getRealValue( tag_a,dret,'hartree') == 0 ) ke_a = dret
             if(printable) write(nfout,'(a,3f10.5)') &
             & ' !** smearing parameters e0, sigma and a for the stress tensor calculation : ',&
             & ke_e0,ke_sigma,ke_a
          endif
          if( f_getIntValue( tag_stress_correction, iret) == 0) sw_stress_correction = iret
          if( sw_stress_correction == ON )then
             if( f_getRealValue( tag_delta_ecut, dret, "rydberg") ==0 ) decut_stress_correction = dret
          endif
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock(tag_lattice) == 0) then
          if( f_getIntValue(tag_sw_optimize_lattice,iret)==0) sw_optimize_lattice = iret
          if( f_getIntValue(tag_sw_optimize_coordinates_once,iret)==0) sw_optimize_coordinates_once = iret
          if( f_getIntValue(tag_sw_rebuild_pws,iret)==0) sw_rebuild_pws = iret
          if(sw_optimize_lattice == ON .and. istress == 0) then
              if(printable) then
                 write(nfout,'(a)') ' !** REMARK: sw_stress is set to ON, since sw_optimize_lattice == ON'
              endif
              istress = ON
          endif
          if ( f_getIntValue(tag_sw_optimize_coords_sametime,iret)==0) then
             sw_optimize_coords_sametime = iret
          endif
          if ( f_getIntValue(tag_lattice_coords_opt_mode,iret)==0) then
             lattice_coords_opt_mode = iret
          endif
          if( f_getRealValue(tag_omega_for_cellopt,dret,'')==0 ) &
               &      omega_for_cellopt = dret
          if( f_getRealValue(tag_stress_force_mul_factor,dret,'')==0 ) &
               &      stress_force_mul_factor = dret
          if( f_getRealValue(tag_threshold_start_cellopt,dret,'')==0 ) &
               &      threshold_start_cellopt = dret

          if( f_getIntValue(tag_sw_uniform,iret)==0 ) sw_uniform = iret
          if( f_getIntValue(tag_nhistory,iret)==0 ) nhistory_stress = iret
          if( f_getRealValue(tag_delta,dret,'')==0) delta_stress = dret
          if( f_getRealValue(tag_max_stress,dret,'hartree/bohr3')==0) max_stress = dret
          if( f_getstringValue(tag_method,rstr,LOWER)==0)then
            if(rstr==tag_quench)then
               lattice_optimization_method = QUENCHED_MD
            else if (rstr==tag_steepest_descent) then
               lattice_optimization_method = STEEPEST_DESCENT
            else if (rstr==tag_bfgs) then
               lattice_optimization_method = BFGS
            endif
          endif

! ========================== JK_add ============================ 13.0AS
          if(.not.estress_has_been_set) then
              external_stress = 0.0d0
          else
              external_stress(1,2) = 0.d0
              external_stress(2,1) = 0.d0
              external_stress(1,3) = 0.d0
              external_stress(3,1) = 0.d0
              external_stress(2,3) = 0.d0
              external_stress(3,2) = 0.d0
          endif
          if( f_selectBlock(tag_external_stress)==0) then
            if(.not.estress_has_been_set)then
                if(f_getRealValue(tag_s11,dret,'hartree/bohr3')==0) external_stress(1,1) = dret
            endif
            if(.not.estress_has_been_set)then
                if(f_getRealValue(tag_s22,dret,'hartree/bohr3')==0) external_stress(2,2) = dret
            endif
            if(.not.estress_has_been_set)then
                if(f_getRealValue(tag_s33,dret,'hartree/bohr3')==0) external_stress(3,3) = dret
            endif
            if(f_getRealValue(tag_s12,dret,'hartree/bohr3')==0) external_stress(1,2) = dret
            if(f_getRealValue(tag_s21,dret,'hartree/bohr3')==0) external_stress(2,1) = dret
            if(f_getRealValue(tag_s13,dret,'hartree/bohr3')==0) external_stress(1,3) = dret
            if(f_getRealValue(tag_s31,dret,'hartree/bohr3')==0) external_stress(3,1) = dret
            if(f_getRealValue(tag_s23,dret,'hartree/bohr3')==0) external_stress(2,3) = dret
            if(f_getRealValue(tag_s32,dret,'hartree/bohr3')==0) external_stress(3,2) = dret
            iret = f_selectParentBlock()
          endif
! =========================================================== 13.0AS

          if(printable.and.ipriinputfile>=1.and.sw_optimize_lattice==ON)then
             write(nfout,'(a)') ' !** parameters for lattice optimization'
             write(nfout,'(a,f20.15)') ' !** max. stress : ',max_stress
             write(nfout,'(a,i5)')     ' !** sw_rebuild_pws : ',sw_rebuild_pws
             write(nfout,'(a,i5)')     ' !** sw_uniform  : ',sw_uniform
             write(nfout,'(a,i5)')     ' !** sw_optimize_coordinates_once  : ',sw_optimize_coordinates_once
             write(nfout,'(a,i5)')     ' !** nhistory    : ',nhistory_stress
             write(nfout,'(a,i5)')     ' !** opt. method : ',lattice_optimization_method

             write(nfout,*) "!** sw_optimize_coords_sametime = ", &
                  &                    sw_optimize_coords_sametime
             if ( sw_optimize_coords_sametime == ON ) then
                write(nfout,*) "!** lattice_coords_opt_mode = ", lattice_coords_opt_mode
                write(nfout,*) "!** omega_for_cellopt = ", omega_for_cellopt
                write(nfout,*) "!** stress_force_mul_factor = ", stress_force_mul_factor
                write(nfout,*) "!** threshold start cellopt = ", threshold_start_cellopt
             endif
! ==================================================13.0B
             write(nfout,*) '!** External stress is set as follows. '
             write(nfout,'(a,F12.8,a,F12.8,a,F12.8)') &
                  &     ' !** s11 = ', external_stress(1,1), &
                  &     '     s12 = ', external_stress(1,2), &
                  &     '     s13 = ', external_stress(1,3)
             write(nfout,'(a,F12.8,a,F12.8,a,F12.8)') &
                  &     ' !** s21 = ', external_stress(2,1), &
                  &     '     s22 = ', external_stress(2,2), &
                  &     '     s23 = ', external_stress(2,3)
             write(nfout,'(a,F12.8,a,F12.8,a,F12.8)') &
                  &     ' !** s31 = ', external_stress(3,1), &
                  &     '     s32 = ', external_stress(3,2), &
                  &     '     s33 = ', external_stress(3,3)
! ================================================ 13.0B
          endif

! =================
          if( f_getIntValue(tag_fix_angle_alpha,iret)==0 ) then
             if ( iret == 1 ) fix_angle_alpha = .true.
          endif
          if( f_getIntValue(tag_fix_angle_beta,iret)==0 ) then
             if ( iret == 1 ) fix_angle_beta = .true.
          endif
          if( f_getIntValue(tag_fix_angle_gamma,iret)==0 ) then
             if ( iret == 1 ) fix_angle_gamma = .true.
          endif
!
          if ( fix_angle_alpha .or. fix_angle_beta .or. fix_angle_gamma ) then
             sw_fix_lattice_angles = on
             if ( sw_uniform == on ) then
                sw_uniform = off
                write(nfout,*) '!** sw_uniform is turned off'
             endif
          endif

          if ( printable.and.ipriinputfile>=1.and.sw_fix_lattice_angles==ON ) then
             write(nfout,*) '!** fix_angle_alpha is ', fix_angle_alpha
             write(nfout,*) '!** fix_angle_beta  is ', fix_angle_beta
             write(nfout,*) '!** fix_angle_gamma is ', fix_angle_gamma
          endif

! =================
          if( f_getIntValue(tag_fix_shape_ab_plane,iret)==0 ) then
             if ( iret == 1 ) then
                fix_shape_ab_plane = .true.;     fix_angle_gamma = .true.
             endif
          endif
          if( f_getIntValue(tag_fix_shape_ac_plane,iret)==0 ) then
             if ( iret == 1 ) then
                fix_shape_ac_plane = .true.;     fix_angle_beta = .true.
             endif
          endif
          if( f_getIntValue(tag_fix_shape_bc_plane,iret)==0 ) then
             if ( iret == 1 ) then
                fix_shape_bc_plane = .true.;     fix_angle_alpha = .true.
             endif
          endif
          if ( fix_shape_ab_plane .or. fix_shape_ac_plane .or. fix_shape_bc_plane ) then
             sw_fix_lattice_shapes = on;
             if ( sw_uniform == on ) then
                sw_uniform = off
                write(nfout,*) '!** sw_uniform is turned off'
             endif
          endif
          if ( printable.and.ipriinputfile>=1.and.sw_fix_lattice_shapes==ON ) then
             write(nfout,*) '!** fix_shape_ab_plane is ', fix_shape_ab_plane
             write(nfout,*) '!** fix_shape_ac_plane is ', fix_shape_ac_plane
             write(nfout,*) '!** fix_shape_bc_plane is ', fix_shape_bc_plane
          endif

! =================
          if( f_getIntValue(tag_fix_length_a,iret)==0 ) then
             if ( iret == 1 ) fix_length_a = .true.
          endif
          if( f_getIntValue(tag_fix_length_b,iret)==0 ) then
             if ( iret == 1 ) fix_length_b = .true.
          endif
          if( f_getIntValue(tag_fix_length_c,iret)==0 ) then
             if ( iret == 1 ) fix_length_c = .true.
          endif
!
          if ( fix_length_a .or. fix_length_b .or. fix_length_c ) then
             sw_fix_lattice_lengths = on
             if ( sw_uniform == on ) then
                sw_uniform = off
                write(nfout,*) '!** sw_uniform is turned off'
             endif
          endif
          if ( printable.and.ipriinputfile>=1.and.sw_fix_lattice_lengths==ON ) then
             write(nfout,*) '!** fix_length_a is ', fix_length_a
             write(nfout,*) '!** fix_length_b is ', fix_length_b
             write(nfout,*) '!** fix_length_c is ', fix_length_c
          endif
! === 2014/11/22
          if( f_getIntValue(tag_sw_neglect_stress_offdiag,iret)==0 ) then
             if ( iret == 1 ) sw_neglect_stress_offdiagonal = on
             write(nfout,*) '!** sw_neglect_stress_offdiagonal is  ', iret
          endif
! === 2014/11/22

! === EXP_CELLOPT === 2015/09/24
          if ( f_getIntValue( tag_sw_read_nfchgt_prev_cell, iret ) ==0 ) then
             sw_read_nfchgt_prev_cell = iret
             write(nfout,*) '!** sw_read_nfchgt_prev_cell is  ', iret
          endif
          if ( f_getIntValue( tag_sw_read_nfzaj_prev_cell, iret ) ==0 ) then
             sw_read_nfzaj_prev_cell = iret
             write(nfout,*) '!** sw_read_nfzaj_prev_cell is  ', iret
          endif
! ================== 2015/09/24

          if (f_getIntValue(tag_sw_interpolate_charge,iret)==0)then
             sw_interpolate_charge = iret
          endif
          if (f_getIntValue(tag_sw_interpolate_wfs,iret)==0)then
             sw_interpolate_wfs = iret
          endif
          if (f_getStringValue(tag_interpolation_method_charge, rstr, LOWER)==0) then
            if(rstr == 'spline') interpolation_method_chg = SPLINE_INTERPOLATION
          endif
          if (f_getStringValue(tag_interpolation_method_wfs, rstr, LOWER)==0) then
            if(rstr == 'spline') interpolation_method_wfs = SPLINE_INTERPOLATION
          endif
          if(printable) then
            write(nfout,'(a,i3)') ' !** sw_interpolate_charge : ',sw_interpolate_charge
            if(sw_interpolate_charge == ON)then
              if(interpolation_method_chg == LINEAR_INTERPOLATION)then
              write(nfout,'(a)')  ' !** interpolation method : LINEAR'
              else if (interpolation_method_chg == SPLINE_INTERPOLATION) then
              write(nfout,'(a)')  ' !** interpolation method : SPLINE'
              endif
            endif
          endif
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_cg) == 0) then
          if( f_getRealValue(tag_etol,dret,'hartree') == 0) etol = dret
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,*) '!** == CG for structure optimization =='
             write(nfout,*) '!** etol =', etol
          end if
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock(tag_predictor) == 0 )then
          if(f_getIntValue(tag_sw_charge_predictor,iret)==0) sw_charge_predictor = iret
          if(f_getIntValue(tag_sw_wf_predictor,iret)==0) sw_wf_predictor = iret
          if(f_getIntValue(tag_sw_extrapolate_charge,iret)==0) sw_extrapolate_charge = iret
          if(f_getRealValue(tag_rms_threshold,dret,'bohr')==0) rms_threshold = dret
          iret = f_selectParentBlock()
       endif
       if(printable .and. ipriinputfile>=1)then
          write(nfout,'(a)')    ' !** parameters for charge and wf predictor'
          write(nfout,'(a,i5)') ' !** sw_charge_predictor : ',sw_charge_predictor
          write(nfout,'(a,i5)') ' !** sw_wf_predictor : ',sw_wf_predictor
          write(nfout,'(a,i5)') ' !** sw_extrapolate_charge : ',sw_extrapolate_charge
          write(nfout,'(a,f13.10,a)') ' !** rms_threshold : ',rms_threshold,' bohr'
       endif

! ============================ KT_add ======================== 13.0B
       if( f_getIntValue(tag_keep_symmetry_strict,iret)==0 ) then
          sw_keep_symmetry_strict = iret
       endif

       if ( printable.and.ipriinputfile>=1.and.sw_keep_symmetry_strict ==ON ) then
          write(nfout,*) ' !**  sw_keep_symmetry_strict is set ON'
       endif
! ============================================================ 13.0B

       if(f_selectBlock(tag_nnp_output) == 0) then
          if(f_getIntValue(tag_sw_nnp_output,iret) == 0) then
             sw_nnp_output = iret
          endif
          if(sw_nnp_output == ON) then
             filetype_nnp = ALL
             if(f_getStringValue(tag_filetype,rstr,LOWER) == 0) then
                if(rstr == tag_xsf) then
                   filetype_nnp = XSF
                else if (rstr == tag_n2p2) then
                   filetype_nnp = N2P2
                else if (rstr == tag_deepmd) then
                   filetype_nnp = DEEPMD
                endif
             endif
             if(f_getIntValue(tag_frequency,iret) == 0) then
                frequency_nnp = iret
             endif
             if(printable .and. ipriinputfile>=1) then
                write(nfout,'(a)')        ' !** sw_nnp_output = ON'
                if(filetype_nnp==XSF) then
                write(nfout,'(a)')        ' !** nnp filetype  = XSF'
                endif
                if(filetype_nnp==N2P2) then
                write(nfout,'(a)')        ' !** nnp filetype  = N2P2'
                endif
                if(filetype_nnp==DEEPMD) then
                write(nfout,'(a)')        ' !** nnp filetype  = DEEPMD'
                endif
                if(filetype_nnp==ALL) then
                write(nfout,'(a)')        ' !** nnp filetype  = XSF and N2P2'
                endif
                write(nfout,'(a,i6)')     ' !** sampling frequency  ',frequency_nnp
                write(nfout,'(a,f20.15)') ' !** edelta for sampling ',edelta_sampling,' hartree'
             endif
          endif
          iret = f_selectParentBlock()
       endif

       iret = f_selectParentBlock()
    end if
  end subroutine m_CtrlP_rd_struc_evol

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_rd_chargemix(nfout)
    integer, intent(in) :: nfout
    integer :: iret, i
    logical :: number_is_given, prealloc
    real(kind=DP) :: rret

    integer :: f_selectTop, f_getIntValue, f_selectParentBlock, f_selectBlock
    integer :: f_getRealValue, f_getstringValue
    integer :: mix_method=-1

    iret = f_selectTop()

    if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** << m_CtrlP_rd_chargemix>>")')
    ! ---- count total number of mixing methods ---
    tag_charge_mixing_is_found = .true.
    cmix_set = PULAY
    if(sw_hubbard==ON) cmix_set = BROYD2
    if( f_selectBlock( tag_charge_mixing) == 0) then
       if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_charge_mixing --")')
       if(f_selectBlock(tag_mixing_methods)==0)then
          cmix_set = EXPLICIT
          iret = f_selectParentBlock()
       endif
       prealloc = .true.
       number_is_given = f_getIntValue(tag_num_mixing_methods,iret) == 0
       if(number_is_given) n_Charge_Mixing_way = iret
       if(cmix_set==EXPLICIT) then
          call set_cdmixingmethods(prealloc,n_Charge_Mixing_way,iret)
          if(iret < 0) call phase_error_with_msg(nfout,' charge mixing ways are not given properly <<m_Ctrl_rd_chargemix>>'&
                                                ,__LINE__,__FILE__)
          if(.not.number_is_given) n_Charge_Mixing_way = iret
          if(number_is_given .and. n_Charge_Mixing_way > iret) n_Charge_Mixing_way = iret
       else
          n_Charge_Mixing_way = 1
       endif
! --> T. Yamasaki 04 Aug. 2009
       if(f_getIntValue(tag_sw_recomposing,iret) == 0) sw_recomposing = iret
       if(f_getRealValue(tag_spin_density_mixfactor,rret,'') == 0) spin_density_mixfactor = rret
! <--

! === KT_add ===== 13.0U3
       if (f_getRealValue(tag_hardpart_mixfactor,rret,'') == 0) hardpart_mixfactor = rret
! ================ 13.0U3

       if(f_getRealValue(tag_occ_matrix_mixfactor,rret,'') == 0) occ_matrix_mixfactor = rret
       if(f_getRealValue(tag_occ_matrix_amix,rret,'') == 0 ) occ_matrix_amix = rret
       if(f_getRealValue(tag_amix_hsr,rret,'') == 0 ) amix_hsr = rret
       if(f_getIntValue(tag_initial_occ_matrix,iret)==0) initial_occ_matrix = iret
       iret = f_selectParentBlock()

    else
       n_Charge_Mixing_way = 1
    end if
    if(ipriinputfile >= 1 .and. printable) &
         & write(nfout,'(" !** n_Charge_Mixing_way = ",i5)') n_Charge_Mixing_way

    call alloc_charge_mixing(n_Charge_Mixing_way)
    if(sw_hubbard==ON) sw_mix_charge_hardpart = ON
    if(PAW_switch==ON) sw_mix_charge_hardpart = ON
    if(cmix_set/=EXPLICIT)then
        call configure_charge_mixing(cmix_set)
    endif

    ! --- setting charge_mixing ---
    prealloc = .false.
    if( f_selectBlock( tag_charge_mixing) == 0) then
       if(cmix_set==EXPLICIT) call set_cdmixingmethods(prealloc,n_Charge_Mixing_way,iret)
       if(f_getstringValue(tag_method,rstr,LOWER)==0)then
           if(rstr==tag_simple) mix_method = SIMPLE
           if(rstr==tag_broyden2) mix_method = BROYD2
           if(rstr==tag_pulay) mix_method = PULAY
           if(mix_method>0)then
              do i=1,n_Charge_Mixing_way
                 charge_mixing(i)%mixing_way = mix_method
              enddo
           endif
       endif
       if(f_getRealValue(tag_rmx,rret,'')==0 )then
          do i=1,n_Charge_Mixing_way
             charge_mixing(i)%rmxs = rret
             charge_mixing(i)%rmxe = rret
          enddo
       endif
       if(f_getIntValue(tag_istr,iret)==0)then
          do i=1,n_Charge_Mixing_way
             charge_mixing(i)%istr = iret
          enddo
       endif
       if(f_getIntValue(tag_nbxmix,iret)==0)then
          do i=1,n_Charge_Mixing_way
             charge_mixing(i)%nbxmix = iret
          enddo
       endif
       iret = f_selectParentBlock()
    end if
!!$    if(ipriinputfile >= 1 .and. printable) then
!!$       write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
!!$       write(nfout,'(" !** spin_density_mixfactor = ",f8.4)') spin_density_mixfactor
!!$       write(nfout,'(" !** --- id, method, rmxs, rmxe, itr, var, prec, istr, nbxmix, update ---")')
!!$       do i = 1, n_Charge_Mixing_way
!!$          write(nfout,'(" !** ",4x,i4,2x,i4,2f7.3,6i5)') i, charge_mixing(i)%mixing_way &
!!$               & ,charge_mixing(i)%rmxs, charge_mixing(i)%rmxe, charge_mixing(i)%iter_range&
!!$               & ,charge_mixing(i)%variation_way, charge_mixing(i)%precon &
!!$               & ,charge_mixing(i)%istr, charge_mixing(i)%nbxmix, charge_mixing(i)%hownew
!!$       end do
!!$       print_sw_recomposing = 1
!!$    end if

    if( f_selectBlock( tag_charge_mixing) == 0) then
       ! --- charge preconditioning ---
       if( f_selectBlock( tag_charge_preconditioning) == 0) then
          if( f_getRealValue(tag_amix,rret,'') == 0) amix = rret
          if( f_getRealValue(tag_bmix,rret,'') == 0) bmix = rret
          if( f_getRealValue(tag_amin,rret,'') == 0) amin = rret
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** --- charge preconditioning ---")')
             write(nfout,'(" !**  amix = ",f8.4,"  bmix = ",f8.4,"  amin = ",f8.4)') amix,bmix,amin
          end if
          if (f_getIntValue(tag_sw_precon_diff,iret)==0) sw_precon_diff = iret
          if (f_getIntValue(tag_sw_metric_diff,iret)==0) sw_metric_diff = iret
          if (f_getIntValue(tag_metric_ratio,iret)==0)   metric_ratio = iret

! ====== KT_add ====== 13.0U3
          if (f_getIntValue(tag_precon_mode,iret)==0) then
             precon_mode = iret
             if ( iret == 1 ) write(nfout,*) "!** precon mode is set to default"
             if ( iret == 2 ) write(nfout,*) "!** precon mode is set to Kresse type"
          endif
          if ( precon_mode == 2 ) then
             amin = 0.1d0;     bmix = 1.0d0;         !! amix = 0.40d0
          endif
! ==================== 13.0U3

          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 2 .and. printable) &
               & write(nfout,'(" !** -- charge preconditioning parameters are not given --")')
       end if

! === KT_add === 13.0U3
!
! --- charge preconditioning hardpart ---
       if ( f_selectBlock( tag_aug_charge_density ) == 0 &
            &  .or.  f_selectBlock( tag_preconditioning_hardpart ) == 0 ) then

          precon_hardpart = .true.

          if (f_getIntValue( tag_sw_wf_mixing,iret )==0 ) then
             sw_wf_mixing = iret
          endif

          if ( sw_wf_mixing == ON ) then
             if ( f_getRealValue( tag_amin_wf_precon,rret,'' ) == 0) then
                amin_wf_precon = rret
             endif
             if ( f_getRealValue( tag_g0_wf_precon,rret,'' ) == 0) then
                g0_wf_precon = rret
             endif
             if ( f_getRealValue( tag_bmix_wf_precon,rret,'' ) == 0) then
                !g0_wf_precon = rret
                bmix_wf_precon = rret
             endif

             if ( f_getRealValue( tag_edelta_start_wf_mixing, rret,'hartree' )== 0) then
                edelta_start_wf_mixing = rret
             endif

             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,*) " !** --- charge preconditioning hardpart ---"
                write(nfout,*) ' !** sw_wf_mixing = ', sw_wf_mixing
                write(nfout,'(" !**  amin = ",f8.4,"  g0(=bmix) = ",f8.4)') &
                     &                amin_wf_precon, g0_wf_precon
                write(nfout,*) " !**  edelta_start_wf_mixing is ",edelta_start_wf_mixing
             end if
          endif
          iret = f_selectParentBlock()
       end if
! ============== 13.0U3

! =========================== added by K. Tagami ==================x====== 5.0
       if (f_getIntValue( tag_sw_mix_charge_hardpart,iret ) == 0 )  then
          sw_mix_charge_hardpart = iret
          if(printable) write(nfout,*) '!** sw_mix_charge_hardpart is set to ', &
        &                     sw_mix_charge_hardpart
       else
          if(printable) write(nfout,*) '!** sw_mix_charge_hardpart is set to default, ', &
        &                     sw_mix_charge_hardpart
       endif
! ======================================================================= 5.0

       if(f_getIntValue(tag_sw_mix_occ_matrix,iret)==0) sw_mix_occ_matrix = iret

! ================================ added by K. Tagami ================== 11.0
       if (f_getIntValue( tag_sw_mix_imaginary_hardpart,iret ) == 0 )  then
          sw_mix_imaginary_hardpart = iret
          if(printable) write(nfout,*) '!** sw_mix_imaginary_hardpart is set to ', &
        &                     sw_mix_imaginary_hardpart
       else
          if(printable) write(nfout,*) &
               & '!** sw_mix_imaginary_hardpart is set to default, ', &
               &                     sw_mix_imaginary_hardpart
       endif
! ======================================================================= 11.0

! === KT_add === 2014/09/19
       if (f_getIntValue( tag_sw_mix_chg_with_ekindens,iret ) == 0 )  then
          sw_mix_charge_with_ekindens = iret
          if(printable) write(nfout,*) '!** sw_mix_charge_with_ekindens is set to ', &
               &                     sw_mix_charge_with_ekindens
       else
          if(printable) write(nfout,*) &
               & '!** sw_mix_charge_with_ekindens is set to default, ', &
               &                     sw_mix_charge_with_ekindens
       endif
! ============= 2014/09/19

       if(f_getIntValue(tag_sw_gradient_simplex,iret)==0) sw_gradient_simplex = iret
       if(f_getIntValue(tag_sw_control_stepsize,iret)==0) sw_control_stepsize = iret
       if(f_getRealValue(tag_alpha,rret,'')==0) then
          alpha_pulay = rret
          alpha_pulay_org = alpha_pulay
       endif
       if(f_getRealValue(tag_damping_factor,rret,'')==0) alpha_pulay_damp = rret
       if(f_getRealValue(tag_threshold_alpha,rret,'')==0) alpha_pulay_damp_thres = rret

! ======= KT_add ========= 13.0U2
       if (f_getIntValue(tag_sw_potential_mixing,iret)==0) then
          sw_potential_mixing = iret
       endif
! ======================== 13.0U2

       ! --- spin density ---
       if( f_selectBlock(tag_spin_density)==0)then
          if(f_getRealValue(tag_spin_density_mixfactor,rret,'')==0) spin_density_mixfactor=rret
          if(f_getIntValue(tag_sw_apply_precon,iret)==0) sw_precon_diff=iret
          if(f_getIntValue(tag_sw_apply_metric,iret)==0) sw_metric_diff=iret
          if(f_getIntValue(tag_sw_force_simple_mixing,iret)==0) sw_force_simple_mixing=iret

! ========================== added by K. Tagami ===================== 5.0
          if ( f_getIntValue( tag_sw_mix_bothspins_sametime, iret )==0 ) then
             sw_mix_bothspins_sametime = iret
             write(nfout,*) '!** sw_mix_bothspins_sametime is set to ', &
        &                    sw_mix_bothspins_sametime
          else
             write(nfout,*) '!** sw_mix_bothspins_sametime is set to default, ', &
        &                    sw_mix_bothspins_sametime
          endif
! ====================================================================== 5.0
          iret = f_selectParentBlock()
       endif

! ================ KT_add ======================== 13.0U2
       ! --- kinetic energy functional ---
       if ( f_selectBlock( tag_kinetic_energy_functional ) == 0 ) then
          if ( f_getIntValue( tag_use_TFW_functional, iret ) == 0 ) then
             sw_modified_TFW_functional = iret
             if ( iret == 1 ) then
                write(nfout,*) "!! sw_modified_TFW_functional is set on"
             endif
          endif

          if ( f_getRealValue( tag_weight_ThomasFermi, rret,'' ) == 0 ) then
             if ( rret >= 0.0 ) weight_TF_functional = rret
          endif
          if ( f_getRealValue( tag_weight_Weizsacker, rret,'' ) == 0 ) then
             if ( rret >= 0.0 ) weight_Weiz_functional = rret
          endif

          if ( f_getIntValue( tag_spherical_averaged_nonlocal, iret ) == 0 ) then
             if ( iret == YES ) sw_spherical_averaged_nonlocal = on
             if ( iret == NO )  sw_spherical_averaged_nonlocal = off
          endif

          if ( f_getIntValue( tag_use_averaged_nonlocal, iret ) == 0 ) then
             if ( iret == YES ) use_averaged_nonlocal = .true.
             if ( iret == NO  ) use_averaged_nonlocal = .false.
          endif
          if ( f_getIntValue( tag_use_deltaW, iret ) == 0 ) then
             if ( iret == YES ) use_deltaW = .true.
             if ( iret == NO  ) use_deltaW = .false.
          endif
          if ( f_getIntValue( tag_use_deltaV, iret ) == 0 ) then
             if ( iret == YES ) use_deltaV = .true.
             if ( iret == NO  ) use_deltaV = .false.
          endif

          if ( f_selectBlock( tag_scf_loop ) == 0 ) then
             if ( f_getIntValue( tag_use_preconditioning, iret ) == 0 ) then
                if ( iret == YES ) use_preconditioning = .true.
                if ( iret == NO  ) use_preconditioning = .false.
             endif
             if ( f_getIntValue( tag_max_iteration_loop, iret ) == 0 ) then
                if ( iret >=0 ) max_iteration_loop = iret
             endif

             if ( f_getRealValue( tag_threshold_enter_loop, rret,'' ) == 0 ) then
                if ( rret > 0.0 ) threshold_enter_loop = rret
             endif
             if ( f_getRealValue( tag_threshold_exit_loop, rret,'' ) == 0 ) then
                if ( rret > 0.0 ) threshold_exit_loop = rret
             endif
             if ( f_getRealValue( tag_threshold_skip_loop, rret,'' ) == 0 ) then
                if ( rret > 0.0 ) threshold_skip_loop = rret
             endif
             iret = f_selectParentBlock()
          endif

          if ( f_selectBlock( tag_line_minimization ) == 0 ) then
             if ( f_getIntValue( tag_max_iteration_linmin, iret ) == 0 ) then
                if ( iret > 0 ) max_iteration_linmin = iret
             endif
             if ( f_getRealValue( tag_threshold_linmin, rret,'' ) == 0 ) then
                if ( rret > 0.0 ) threshold_linmin = rret
             endif
             iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       endif

       if ( sw_modified_TFW_functional == ON ) then
          sw_mix_charge_hardpart = on
       endif
! ================================================ 13.0U2

       if(printable .and. ipriinputfile>=1) then
          write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
          write(nfout,'(" !** metric_ratio   = ",i5)') metric_ratio
          print_sw_recomposing = 1
       endif
       if(sw_recomposing==ON .and. printable.and. ipriinputfile>=1)then
          write(nfout,'(" !** --- spin density ---")')
          write(nfout,'(a,i5)') ' !**  force simple mixing to spin density (sw_force_simple_mixing) : ',sw_force_simple_mixing
!234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
          write(nfout,'(a,f8.4)') ' !**  spin_density_mixfactor = ', spin_density_mixfactor
          write(nfout,'(a,i5)') ' !**  apply preconditioning to spin density (sw_precon_diff) :       ',sw_precon_diff
          write(nfout,'(a,i5)') ' !**  apply metric to spin density (sw_metric_diff) :                ',sw_metric_diff
       endif
       if(printable.and.ipriinputfile>=1)then
          write(nfout,'(a   )') ' !**  use the gradient simplex method for storing the history of charge'
          write(nfout,'(a,i5)') ' !**                                   ( sw_gradient_simplex) :      ',sw_gradient_simplex
       endif

       max_stepsize = charge_mixing(1)%rmxs*1.5d0
       if(max_stepsize>0.8d0) max_stepsize=0.8d0
       if(f_getRealValue(tag_max_stepsize,rret,'')==0) max_stepsize=rret

       if(f_getRealValue(tag_ommix_factor,rret,'')==0) ommix_factor = rret
! =========== KT_add ========== 13.0U2
       if (printable .and. ipriinputfile>=1 )then
          write(nfout,'(a,i5)') ' !**  sw_potential_mixing is set to ', &
               &                       sw_potential_mixing
          if ( sw_modified_TFW_functional == ON ) then
             write(nfout,*) "**** TFW functional info ***"
             write(nfout,*) '!! weight_ThomasFermi is ', weight_TF_functional
             write(nfout,*) '!! weight_Weizsacker is ',  weight_Weiz_functional

             write(nfout,*) "!! sw_spherical_averaged_nonlocal is ", &
                  &                       sw_spherical_averaged_nonlocal
             write(nfout,*) "!! use_averaged_nonlocal is ", use_averaged_nonlocal
             write(nfout,*) "!! use_deltaW is ", use_deltaW
             write(nfout,*) "!! use_deltaV is ", use_deltaV

             write(nfout,*) '!! max_iteration_loop is ',  max_iteration_loop
             write(nfout,*) "!! use_precoditioning is ", use_preconditioning

             write(nfout,*) '!! threshold_enter_loop is ', threshold_enter_loop
             write(nfout,*) '!! threshold_exit_loop is ', threshold_exit_loop
             write(nfout,*) '!! threshold_skip_loop is ', threshold_skip_loop

             write(nfout,*) '!! max_iteration_linmin is ',  max_iteration_linmin
             write(nfout,*) '!! threshold_linmin is ', threshold_linmin

             write(nfout,*) "**** ***"
          endif
       endif
! ============================= 13.0U2

       iret = f_selectParentBlock()
    end if

! =============================== added by K. Tagami =============== 5.0
       sw_force_simple_mixing_hsr = sw_force_simple_mixing
       sw_recomposing_hsr = sw_recomposing
       sw_mix_bothspins_sametime_hsr = sw_mix_bothspins_sametime
! ================================================================== 5.0

    if (printable .and. ipriinputfile >=1 .and. print_sw_recomposing==0) then
       write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
       write(nfout,'(" !** metric ratio   = ",i5)') metric_ratio
       if(sw_recomposing==ON) then
          write(nfout,'(" !** --- spin density ---")')
          write(nfout,'(a,i5)') ' !**  force simple mixing to spin density (sw_force_simple_mixing) : ',sw_force_simple_mixing
          write(nfout,'(a,f8.4)') ' !**  spin_density_mixfactor = ', spin_density_mixfactor
          write(nfout,'(a,i5)') ' !**  apply preconditioning to spin density (sw_precon_diff) :       ',sw_precon_diff
          write(nfout,'(a,i5)') ' !**  apply metric to spin density (sw_metric_diff) :                ',sw_metric_diff
       end if
       write(nfout,'(" !** --- id, method, rmxs, rmxe, itr, var, prec, istr, nbxmix, update ---")')
       do i = 1, n_Charge_Mixing_way
          write(nfout,'(" !** ",4x,i4,2x,i4,2f7.3,6i5)') i, charge_mixing(i)%mixing_way &
               & ,charge_mixing(i)%rmxs, charge_mixing(i)%rmxe, charge_mixing(i)%iter_range&
               & ,charge_mixing(i)%variation_way, charge_mixing(i)%precon &
               & ,charge_mixing(i)%istr, charge_mixing(i)%nbxmix, charge_mixing(i)%hownew
       end do
       print_sw_recomposing = 1
    end if

  contains
    subroutine set_cdmixingmethods(prealloc,m,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: m
      integer, intent(out) :: iret
      integer :: i, icf(9)
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: mway, itr, var, prec, hownew, istr, nbxmix, no
      real(kind=DP) :: rmxs, rmxe

      if( f_selectBlock(tag_mixing_methods) == 0) then
         i = 1
         do while(.true.)
            if(i==1) then
               if( f_selectFirstTableLine() /= 0) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0) then
                  exit
               end if
            end if
            if(.not.prealloc) then
               if(i > m) exit
               no = i
               iret = f_readchargemixing(icf,no,mway,rmxs,rmxe,itr,var,prec,istr,nbxmix,hownew)
               if(icf(1)==1) charge_mixing(no)%mixing_way  = mway
               if(icf(2)==1) charge_mixing(no)%rmxs   = rmxs
               if(icf(3)==1) then
                  charge_mixing(no)%rmxe   = rmxe
               else if(icf(2) == 1) then
                  charge_mixing(no)%rmxe   = charge_mixing(no)%rmxs
               end if
               if(icf(4)==1) charge_mixing(no)%iter_range  = itr
               if(icf(5)==1) charge_mixing(no)%variation_way = var
               if(icf(6)==1) charge_mixing(no)%precon = prec
               if(icf(7)==1) charge_mixing(no)%istr   = istr
               if(icf(8)==1) charge_mixing(no)%nbxmix = nbxmix
               if(icf(9)==1) charge_mixing(no)%hownew = hownew
            end if
            i = i + 1
         end do
         iret = f_selectParentBlock()
      else
         call phase_error_with_msg(nfout,&
         ' ! No charge-mixing method is given in the inputfile <<m_CtrlP_rd_chargemix.set_cdmixingmethods>>', &
         __LINE__,__FILE__)
      end if
      iret = i -1
    end subroutine set_cdmixingmethods

    subroutine configure_charge_mixing(cmixset)
        integer, intent(in) :: cmixset
        cmix_explicit = .false.
        if (cmixset == BROYD2 .or. cmixset==PULAY)then
           if(cmixset==BROYD2)then
              charge_mixing(1)%mixing_way = BROYD2
              sw_mix_bothspins_sametime = ON
              if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : broyden2'
           else
              charge_mixing(1)%mixing_way = PULAY
              sw_mix_bothspins_sametime = ON
              if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : pulay'
           endif
           charge_mixing(1)%rmxs          = 0.4d0
           charge_mixing(1)%istr          = 3
           charge_mixing(1)%nbxmix        = 15
           if(nspin>1) then
              charge_mixing(1)%rmxs = 0.1d0
              charge_mixing(1)%istr = 5
              charge_mixing(1)%nbxmix = 5
           endif
           charge_mixing(1)%rmxe = charge_mixing(1)%rmxs
           if(sw_hubbard==ON)then
             charge_mixing(1)%nbxmix  = 5
             charge_mixing(1)%hownew  = ANEW
           endif
        else if (cmix_set == SIMPLE)then
           charge_mixing(1)%mixing_way = SIMPLE
           charge_mixing(1)%rmxs          = 0.4d0
           if(nspin>1) charge_mixing(1)%rmxs = 0.1d0
           charge_mixing(1)%rmxe = charge_mixing(1)%rmxs
           if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : simple'
        end if
        if(printable) write(nfout,'(" !** sw_mix_bothspins_sametime = ",i3," <<configuration_charge_mixing>>")') &
                      sw_mix_bothspins_sametime
    end subroutine configure_charge_mixing

  end subroutine m_CtrlP_rd_chargemix

  subroutine m_CtrlP_rd_chargemix2(nfout, reread)
    integer, intent(in) :: nfout
    logical, intent(in), optional :: reread
    integer :: iret, i
    logical :: number_is_given, prealloc
    real(kind=DP) :: rret

    integer :: f_selectTop, f_getIntValue, f_selectParentBlock, f_selectBlock
    integer :: f_getRealValue, f_getstringValue
    integer :: mix_method=-1
    integer :: n_Charge_Mixing_way_l
    logical :: rr, done_something
    rr = .false.
    if(present(reread)) rr = reread

    iret = f_selectTop()

    if(ipriinputfile >= 2 .and. printable .and. .not. rr) write(nfout,'(" !** << m_CtrlP_rd_chargemix>>")')
    ! ---- count total number of mixing methods ---
    tag_charge_mixing_is_found = .true.
    cmix_set = PULAY
    if(sw_hubbard==ON) cmix_set = BROYD2
    if( f_selectBlock( tag_charge_mixing) == 0) then
       if(ipriinputfile >= 2 .and. printable .and. .not. rr) write(nfout,'(" !** -- tag_charge_mixing --")')
       if(f_selectBlock(tag_mixing_methods)==0)then
          cmix_set = EXPLICIT
          iret = f_selectParentBlock()
       endif
       prealloc = .true.
       number_is_given = f_getIntValue(tag_num_mixing_methods,iret) == 0
       if(number_is_given) n_Charge_Mixing_way_l = iret
       if(cmix_set==EXPLICIT) then
          call set_cdmixingmethods(prealloc,n_Charge_Mixing_way_l,iret)
          if(iret < 0) call phase_error_with_msg(nfout,' charge mixing ways are not given properly <<m_Ctrl_rd_chargemix>>'&
                                                ,__LINE__,__FILE__)
          if(.not.number_is_given) n_Charge_Mixing_way_l = iret
          if(number_is_given .and. n_Charge_Mixing_way_l > iret) n_Charge_Mixing_way_l = iret
       else
          n_Charge_Mixing_way_l = 1
       endif
! --> T. Yamasaki 04 Aug. 2009
       call m_CtrlP_rd_val(nfout, tag_sw_recomposing, sw_recomposing, rr)
       call m_CtrlP_rd_val(nfout, tag_spin_density_mixfactor,'',spin_density_mixfactor,rr)
! <--

! === KT_add ===== 13.0U3
       call m_CtrlP_rd_val(nfout, tag_hardpart_mixfactor, '', hardpart_mixfactor,rr)
! ================ 13.0U3

       call m_CtrlP_rd_val(nfout, tag_occ_matrix_mixfactor, '', occ_matrix_mixfactor, rr)
       call m_CtrlP_rd_val(nfout, tag_occ_matrix_amix, '', occ_matrix_amix, rr)
       call m_CtrlP_rd_val(nfout, tag_amix_hsr, '', amix_hsr, rr)
       call m_CtrlP_rd_val(nfout, tag_initial_occ_matrix, initial_occ_matrix, rr)
       iret = f_selectParentBlock()
    else
       n_Charge_Mixing_way_l = 1
    end if
    if( .not. rr) then
       n_Charge_Mixing_way = n_Charge_Mixing_way_l
    else
       if(n_Charge_Mixing_way_l /= n_Charge_Mixing_way) then
         call dealloc_charge_mixing()
         call alloc_charge_mixing(n_Charge_Mixing_way_l)
       endif
       n_Charge_Mixing_way = n_Charge_Mixing_way_l
    endif
    if(ipriinputfile >= 1 .and. printable .and. .not. rr) &
         & write(nfout,'(" !** n_Charge_Mixing_way = ",i5)') n_Charge_Mixing_way
    if(.not. rr) call alloc_charge_mixing(n_Charge_Mixing_way)

    if(sw_hubbard==ON) sw_mix_charge_hardpart = ON
    if(PAW_switch==ON) sw_mix_charge_hardpart = ON
    if(cmix_set/=EXPLICIT .and. .not. rr)then
        call configure_charge_mixing(cmix_set)
    endif

    ! --- setting charge_mixing ---
    prealloc = .false.
    if( f_selectBlock( tag_charge_mixing) == 0) then
       if(cmix_set==EXPLICIT) call set_cdmixingmethods(prealloc,n_Charge_Mixing_way,iret)
       if(f_getstringValue(tag_method,rstr,LOWER)==0)then
           if(rstr==tag_simple) mix_method = SIMPLE
           if(rstr==tag_broyden2) mix_method = BROYD2
           if(rstr==tag_pulay) mix_method = PULAY
           if(mix_method>0)then
              do i=1,n_Charge_Mixing_way
                 if(charge_mixing(i)%mixing_way /= mix_method) then
                   if(rr) write(nfout,'(a,i8,a,i8)') '! ** F_INP_MOD changed '//&
                   &  trim(tag_method)//' to ',mix_method,' from ',charge_mixing(i)%mixing_way
                   charge_mixing(i)%mixing_way = mix_method
                 endif
              enddo
           endif
       endif
       do i=1,n_Charge_Mixing_way
          call m_CtrlP_rd_val(nfout, tag_rmx, '', &
          &    charge_mixing(i)%rmxs,rr,msg='rmxs')
          call m_CtrlP_rd_val(nfout, tag_rmx, '', &
          &    charge_mixing(i)%rmxe,rr,msg='rmxe')
       enddo
       do i=1,n_Charge_Mixing_way
          call m_CtrlP_rd_val(nfout, tag_istr, charge_mixing(i)%istr,rr)
       enddo
       do i=1,n_Charge_Mixing_way
          call m_CtrlP_rd_val(nfout, tag_nbxmix, charge_mixing(i)%nbxmix, rr)
       enddo
       iret = f_selectParentBlock()
    end if
!!$    if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
!!$       write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
!!$       write(nfout,'(" !** spin_density_mixfactor = ",f8.4)') spin_density_mixfactor
!!$       write(nfout,'(" !** --- id, method, rmxs, rmxe, itr, var, prec, istr, nbxmix, update ---")')
!!$       do i = 1, n_Charge_Mixing_way
!!$          write(nfout,'(" !** ",4x,i4,2x,i4,2f7.3,6i5)') i, charge_mixing(i)%mixing_way &
!!$               & ,charge_mixing(i)%rmxs, charge_mixing(i)%rmxe, charge_mixing(i)%iter_range&
!!$               & ,charge_mixing(i)%variation_way, charge_mixing(i)%precon &
!!$               & ,charge_mixing(i)%istr, charge_mixing(i)%nbxmix, charge_mixing(i)%hownew
!!$       end do
!!$       print_sw_recomposing = 1
!!$    end if

    if( f_selectBlock( tag_charge_mixing) == 0) then
       ! --- charge preconditioning ---
       if( f_selectBlock( tag_charge_preconditioning) == 0) then
          call m_CtrlP_rd_val(nfout, tag_amix, '', amix, rr)
          call m_CtrlP_rd_val(nfout, tag_bmix, '', bmix, rr)
          call m_CtrlP_rd_val(nfout, tag_amin, '', amin, rr)
          if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
             write(nfout,'(" !** --- charge preconditioning ---")')
             write(nfout,'(" !**  amix = ",f8.4,"  bmix = ",f8.4,"  amin = ",f8.4)') amix,bmix,amin
          end if
          call m_CtrlP_rd_val(nfout, tag_sw_precon_diff, sw_precon_diff, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_metric_diff, sw_metric_diff, rr)
          call m_CtrlP_rd_val(nfout, tag_metric_ratio, metric_ratio, rr)

! ====== KT_add ====== 13.0U3
          call m_CtrlP_rd_val(nfout, tag_precon_mode, precon_mode, rr)
          if ( precon_mode == 1 .and. .not. rr) write(nfout,*) "!** precon mode is set to default"
          if ( precon_mode == 2 .and. .not. rr) write(nfout,*) "!** precon mode is set to Kresse type"
          if ( precon_mode == 2 ) then
             amin = 0.1d0;     bmix = 1.0d0;         !! amix = 0.40d0
          endif
! ==================== 13.0U3

          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 2 .and. printable .and. .not. rr) &
               & write(nfout,'(" !** -- charge preconditioning parameters are not given --")')
       end if

! === KT_add === 13.0U3
!
! --- charge preconditioning hardpart ---
       if ( f_selectBlock( tag_aug_charge_density ) == 0 &
            &  .or.  f_selectBlock( tag_preconditioning_hardpart ) == 0 ) then

          precon_hardpart = .true.

          call m_CtrlP_rd_val(nfout, tag_sw_wf_mixing, sw_wf_mixing, rr)

          if ( sw_wf_mixing == ON ) then
             call m_CtrlP_rd_val(nfout, tag_amin_wf_precon, '', amin_wf_precon,rr)
             call m_CtrlP_rd_val(nfout, tag_g0_wf_precon,'',g0_wf_precon,rr)
             call m_CtrlP_rd_val(nfout, tag_bmix_wf_precon, '', bmix_wf_precon,rr)
             call m_CtrlP_rd_val(nfout, tag_edelta_start_wf_mixing, 'hartree', &
             &    edelta_start_wf_mixing, rr)

             if(ipriinputfile >= 1 .and. printable .and. .not. rr) then
                write(nfout,*) " !** --- charge preconditioning hardpart ---"
                write(nfout,*) ' !** sw_wf_mixing = ', sw_wf_mixing
                write(nfout,'(" !**  amin = ",f8.4,"  g0(=bmix) = ",f8.4)') &
                     &                amin_wf_precon, g0_wf_precon
                write(nfout,*) " !**  edelta_start_wf_mixing is ",edelta_start_wf_mixing
             end if
          endif
          iret = f_selectParentBlock()
       end if
! ============== 13.0U3

! =========================== added by K. Tagami ==================x====== 5.0
       call m_CtrlP_rd_val(nfout, tag_sw_mix_charge_hardpart, &
       &    sw_mix_charge_hardpart, rr, done_something)
       if(done_something) then
         if(printable .and. .not. rr) write(nfout,*) '!** sw_mix_charge_hardpart is set to ', &
         &             sw_mix_charge_hardpart
       else
         if(printable .and. .not. rr) write(nfout,*) '!** sw_mix_charge_hardpart is set to default, ', &
         &             sw_mix_charge_hardpart
       endif
! ======================================================================= 5.0

       call m_CtrlP_rd_val(nfout, tag_sw_mix_occ_matrix, sw_mix_occ_matrix, rr)

! ================================ added by K. Tagami ================== 11.0
       call m_CtrlP_rd_val(nfout, tag_sw_mix_imaginary_hardpart, &
       &    sw_mix_imaginary_hardpart, rr, done_something)
       if (done_something) then
          if(printable .and. .not. rr) write(nfout,*) &
               & '!** sw_mix_imaginary_hardpart is set to ', &
               &                     sw_mix_imaginary_hardpart
       else
          if(printable .and. .not. rr) write(nfout,*) &
               & '!** sw_mix_imaginary_hardpart is set to default, ', &
               &                     sw_mix_imaginary_hardpart
       endif
! ======================================================================= 11.0

! === KT_add === 2014/09/19
       call m_CtrlP_rd_val(nfout, tag_sw_mix_chg_with_ekindens, &
       & sw_mix_charge_with_ekindens, rr, done_something)
       if (done_something) then
          if(printable .and. .not. rr) write(nfout,*) &
               & '!** sw_mix_charge_with_ekindens is set to ', &
               &                       sw_mix_charge_with_ekindens
       else
          if(printable .and. .not. rr) write(nfout,*) &
               & '!** sw_mix_charge_with_ekindens is set to default, ', &
               &                     sw_mix_charge_with_ekindens
       endif
! ============= 2014/09/19

       call m_CtrlP_rd_val(nfout, tag_sw_gradient_simplex, &
       &    sw_gradient_simplex,rr)
       call m_CtrlP_rd_val(nfout, tag_sw_control_stepsize, &
       &    sw_control_stepsize,rr)
       call m_CtrlP_rd_val(nfout, tag_alpha,'', alpha_pulay, &
       &    rr, done_something)
       if(done_something) then
          alpha_pulay_org = alpha_pulay
       endif
       call m_CtrlP_rd_val(nfout, tag_damping_factor, '', alpha_pulay_damp, rr)
       call m_CtrlP_rd_val(nfout, tag_threshold_alpha, '', alpha_pulay_damp_thres, rr)

! ======= KT_add ========= 13.0U2
       call m_CtrlP_rd_val(nfout, tag_sw_potential_mixing, sw_potential_mixing,&
       &    rr)
! ======================== 13.0U2

       ! --- spin density ---
       if( f_selectBlock(tag_spin_density)==0)then
          if(ipriinputfile >= 1 .and. printable ) &
               & write(nfout,'(" !** tag_spin_density is found")')
          call m_CtrlP_rd_val(nfout, tag_spin_density_mixfactor, '', &
          &    spin_density_mixfactor, rr)
          call m_CtrlP_rd_val(nfout, tag_sw_apply_precon,sw_precon_diff,rr)
          call m_CtrlP_rd_val(nfout, tag_sw_apply_metric,sw_metric_diff,rr)
          call m_CtrlP_rd_val(nfout, tag_sw_force_simple_mixing, &
          &    sw_force_simple_mixing,rr)

! ========================== added by K. Tagami ===================== 5.0
          call m_CtrlP_rd_val(nfout, tag_sw_mix_bothspins_sametime, &
          &    sw_mix_bothspins_sametime, rr, done_something)
          if ( done_something ) then
             if(.not. rr) write(nfout,*) '!** sw_mix_bothspins_sametime is set to ', &
                  &                              sw_mix_bothspins_sametime
          else
             if(.not.rr) write(nfout,*) '!** sw_mix_bothspins_sametime is set to default, ', &
                  &                             sw_mix_bothspins_sametime
          endif
! ====================================================================== 5.0
          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 1 .and. printable ) &
               & write(nfout,'(" !** tag_spin_density is not found")')
       endif

! ================ KT_add ======================== 13.0U2
       ! --- kinetic energy functional ---
       if ( f_selectBlock( tag_kinetic_energy_functional ) == 0 ) then
          call m_CtrlP_rd_val(nfout, tag_use_TFW_functional, &
          & sw_modified_TFW_functional, rr)
          if ( sw_modified_TFW_functional == 1 ) then
             write(nfout,*) "!! sw_modified_TFW_functional is set on"
          endif

          call m_CtrlP_rd_val(nfout, tag_weight_ThomasFermi,'', &
          &    weight_TF_functional, rr, positive=.true.)
          call m_CtrlP_rd_val(nfout, tag_weight_Weizsacker,'', &
          &    weight_Weiz_functional, rr, positive=.true.)

          call m_CtrlP_rd_val(nfout, tag_spherical_averaged_nonlocal, &
          &    sw_spherical_averaged_nonlocal, rr)

          if ( f_getIntValue( tag_use_averaged_nonlocal, iret ) == 0 ) then
             if ( iret == YES ) use_averaged_nonlocal = .true.
             if ( iret == NO  ) use_averaged_nonlocal = .false.
          endif
          if ( f_getIntValue( tag_use_deltaW, iret ) == 0 ) then
             if ( iret == YES ) use_deltaW = .true.
             if ( iret == NO  ) use_deltaW = .false.
          endif
          if ( f_getIntValue( tag_use_deltaV, iret ) == 0 ) then
             if ( iret == YES ) use_deltaV = .true.
             if ( iret == NO  ) use_deltaV = .false.
          endif

          if ( f_selectBlock( tag_scf_loop ) == 0 ) then
             if ( f_getIntValue( tag_use_preconditioning, iret ) == 0 ) then
                if ( iret == YES ) use_preconditioning = .true.
                if ( iret == NO  ) use_preconditioning = .false.
             endif
             call m_CtrlP_rd_val(nfout, tag_max_iteration_loop, &
             &    max_iteration_loop, rr, positive = .true.)

             call m_CtrlP_rd_val(nfout, tag_threshold_enter_loop, '', &
             &    threshold_enter_loop, rr, positive = .true.)

             call m_CtrlP_rd_val(nfout, tag_threshold_exit_loop, '', &
             &    threshold_exit_loop, rr, positive = .true.)

             call m_CtrlP_rd_val(nfout, tag_threshold_skip_loop, '', &
             &    threshold_skip_loop, rr, positive = .true.)
             iret = f_selectParentBlock()
          endif

          if ( f_selectBlock( tag_line_minimization ) == 0 ) then
             call m_CtrlP_rd_val(nfout, tag_max_iteration_linmin, &
             &    max_iteration_linmin, rr, positive = .true.)
             call m_CtrlP_rd_val(nfout, tag_threshold_linmin,'', &
             &    threshold_linmin, rr, positive = .true.)
             iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       endif

       if ( sw_modified_TFW_functional == ON ) then
          sw_mix_charge_hardpart = on
       endif
! ================================================ 13.0U2

       if(printable .and. ipriinputfile>=1 .and. .not. rr) then
          write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
          write(nfout,'(" !** metric_ratio   = ",i5)') metric_ratio
          print_sw_recomposing = 1
       endif
       if(sw_recomposing==ON .and. printable.and. ipriinputfile>=1 .and. &
       &  .not.  rr)then
          write(nfout,'(" !** --- spin density ---")')
          write(nfout,'(a,i5)') ' !**  force simple mixing to spin density (sw_force_simple_mixing) : ',sw_force_simple_mixing
          write(nfout,'(a,f8.4)') ' !**  spin_density_mixfactor = ',spin_density_mixfactor
          write(nfout,'(a,i5)') ' !**  apply preconditioning to spin density (sw_precon_diff) :       ',sw_precon_diff
          write(nfout,'(a,i5)') ' !**  apply metric to spin density (sw_metric_diff) :                ',sw_metric_diff
       endif
       if(printable.and.ipriinputfile>=1 .and. .not. rr)then
          write(nfout,'(a   )') ' !**  use the gradient simplex method for storing the history of charge'
          write(nfout,'(a,i5)') ' !**                                   ( sw_gradient_simplex) :      ',sw_gradient_simplex
       endif

       max_stepsize = charge_mixing(1)%rmxs*1.5d0
       if(max_stepsize>0.8d0) max_stepsize=0.8d0
       call m_CtrlP_rd_val(nfout, tag_max_stepsize, '', max_stepsize, rr)
       call m_CtrlP_rd_val(nfout, tag_ommix_factor, '', ommix_factor, rr)
! =========== KT_add ========== 13.0U2
       if (printable .and. ipriinputfile>=1 .and. .not. rr)then
          write(nfout,'(a,i5)') ' !**  sw_potential_mixing is set to ', &
               &                       sw_potential_mixing
          if ( sw_modified_TFW_functional == ON ) then
             write(nfout,*) "**** TFW functional info ***"
             write(nfout,*) '!! weight_ThomasFermi is ', weight_TF_functional
             write(nfout,*) '!! weight_Weizsacker is ',  weight_Weiz_functional

             write(nfout,*) "!! sw_spherical_averaged_nonlocal is ", &
                  &                       sw_spherical_averaged_nonlocal
             write(nfout,*) "!! use_averaged_nonlocal is ", use_averaged_nonlocal
             write(nfout,*) "!! use_deltaW is ", use_deltaW
             write(nfout,*) "!! use_deltaV is ", use_deltaV

             write(nfout,*) '!! max_iteration_loop is ',  max_iteration_loop
             write(nfout,*) "!! use_precoditioning is ", use_preconditioning

             write(nfout,*) '!! threshold_enter_loop is ', threshold_enter_loop
             write(nfout,*) '!! threshold_exit_loop is ', threshold_exit_loop
             write(nfout,*) '!! threshold_skip_loop is ', threshold_skip_loop

             write(nfout,*) '!! max_iteration_linmin is ',  max_iteration_linmin
             write(nfout,*) '!! threshold_linmin is ', threshold_linmin

             write(nfout,*) "**** ***"
          endif
       endif
! ============================= 13.0U2

       iret = f_selectParentBlock()
    end if

! =============================== added by K. Tagami =============== 5.0
       sw_force_simple_mixing_hsr = sw_force_simple_mixing
       sw_recomposing_hsr = sw_recomposing
       sw_mix_bothspins_sametime_hsr = sw_mix_bothspins_sametime
! ================================================================== 5.0

    if (printable .and. ipriinputfile >=1 .and. .not. rr .and. print_sw_recomposing==0) then
       write(nfout,'(" !** sw_recomposing = ",i5)') sw_recomposing
       write(nfout,'(" !** metric_ratio   = ",i5)') metric_ratio
       if(sw_recomposing==ON) then
          write(nfout,'(" !** --- spin density ---")')
          write(nfout,'(a,i5)') ' !**  force simple mixing to spin density (sw_force_simple_mixing) : ',sw_force_simple_mixing
          write(nfout,'(a,f8.4)') ' !**  spin_density_mixfactor = ', spin_density_mixfactor
          write(nfout,'(a,i5)') ' !**  apply preconditioning to spin density (sw_precon_diff) :       ',sw_precon_diff
          write(nfout,'(a,i5)') ' !**  apply metric to spin density (sw_metric_diff) :                ',sw_metric_diff
       end if
       write(nfout,'(" !** --- id, method, rmxs, rmxe, itr, var, prec, istr, nbxmix, update ---")')
       do i = 1, n_Charge_Mixing_way
          write(nfout,'(" !** ",4x,i4,2x,i4,2f7.3,6i5)') i, charge_mixing(i)%mixing_way &
               & ,charge_mixing(i)%rmxs, charge_mixing(i)%rmxe, charge_mixing(i)%iter_range&
               & ,charge_mixing(i)%variation_way, charge_mixing(i)%precon &
               & ,charge_mixing(i)%istr, charge_mixing(i)%nbxmix, charge_mixing(i)%hownew
       end do
       print_sw_recomposing = 1
    end if

  contains

    subroutine set_cdmixingmethods(prealloc,m,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: m
      integer, intent(out) :: iret
      integer :: i, icf(9)
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: mway, itr, var, prec, hownew, istr, nbxmix, no
      real(kind=DP) :: rmxs, rmxe

      if( f_selectBlock(tag_mixing_methods) == 0) then
         i = 1
         do while(.true.)
            if(i==1) then
               if( f_selectFirstTableLine() /= 0) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0) then
                  exit
               end if
            end if
            if(.not.prealloc) then
               if(i > m) exit
               no = i
               iret = f_readchargemixing(icf,no,mway,rmxs,rmxe,itr,var,prec,istr,nbxmix,hownew)
               if(icf(1)==1 .and. charge_mixing(no)%mixing_way /= mway) charge_mixing(no)%mixing_way  = mway
               if(icf(2)==1 .and. charge_mixing(no)%rmxs /= rmxs) charge_mixing(no)%rmxs   = rmxs
               if(icf(3)==1) then
                  if(charge_mixing(no)%rmxe /= rmxe ) charge_mixing(no)%rmxe   = rmxe
               else if(icf(2) == 1) then
                  if(charge_mixing(no)%rmxe /= charge_mixing(no)%rmxs) &
                  &  charge_mixing(no)%rmxe   = charge_mixing(no)%rmxs
               end if
               if(icf(4)==1 .and. charge_mixing(no)%iter_range /= itr) charge_mixing(no)%iter_range  = itr
               if(icf(5)==1 .and. charge_mixing(no)%variation_way /= var) &
               & charge_mixing(no)%variation_way = var
               if(icf(6)==1 .and. charge_mixing(no)%precon /= prec) charge_mixing(no)%precon = prec
               if(icf(7)==1 .and. charge_mixing(no)%istr /= istr) charge_mixing(no)%istr   = istr
               if(icf(8)==1 .and. charge_mixing(no)%nbxmix /= nbxmix) charge_mixing(no)%nbxmix = nbxmix
               if(icf(9)==1 .and. charge_mixing(no)%hownew /= hownew) charge_mixing(no)%hownew = hownew
            end if
            i = i + 1
         end do
         iret = f_selectParentBlock()
      else
         call phase_error_with_msg(nfout,&
         ' ! No charge-mixing method is given in the inputfile <<m_CtrlP_rd_chargemix.set_cdmixingmethods>>', &
         __LINE__,__FILE__)
      end if
      iret = i -1
    end subroutine set_cdmixingmethods

    subroutine configure_charge_mixing(cmixset)
        integer, intent(in) :: cmixset
        cmix_explicit = .false.
        if (cmixset == BROYD2 .or. cmixset==PULAY)then
           if(cmixset==BROYD2)then
              charge_mixing(1)%mixing_way = BROYD2
              sw_mix_bothspins_sametime = ON
              if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : broyden2'
           else
              charge_mixing(1)%mixing_way = PULAY
              sw_mix_bothspins_sametime = ON
              if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : pulay'
           endif
           charge_mixing(1)%rmxs          = 0.4d0
           charge_mixing(1)%istr          = 3
           charge_mixing(1)%nbxmix        = 15
           if(nspin>1) then
              charge_mixing(1)%rmxs = 0.1d0
              charge_mixing(1)%istr = 5
              charge_mixing(1)%nbxmix = 5
           endif
           charge_mixing(1)%rmxe = charge_mixing(1)%rmxs
           if(sw_hubbard==ON)then
             charge_mixing(1)%nbxmix  = 5
             charge_mixing(1)%hownew  = ANEW
           endif
        else if (cmix_set == SIMPLE)then
           charge_mixing(1)%mixing_way = SIMPLE
           charge_mixing(1)%rmxs          = 0.4d0
           if(nspin>1) charge_mixing(1)%rmxs = 0.1d0
           charge_mixing(1)%rmxe = charge_mixing(1)%rmxs
           if(printable) write(nfout,'(a)') ' !** applied charge-mixing method : simple'
        end if
    end subroutine configure_charge_mixing

  end subroutine m_CtrlP_rd_chargemix2

  logical function m_CtrlP_explicit_cmix()
     m_CtrlP_explicit_cmix = cmix_explicit
  end function m_CtrlP_explicit_cmix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intelligent reading of charge-density mixing methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function f_readchargemixing(icf,no,mway,rmxs,rmxe,itr,var,prec,istr,nbxmix,hownew)
!  Coded by T. Yamasaki (FUJITSU LABORATORIES LTD), 13th Jun 2003
    integer, intent(out), dimension(9) :: icf
    integer, intent(inout) :: no
    integer, intent(out) ::  mway,itr,var,prec,istr,nbxmix,hownew

    real(kind=DP),intent(inout) :: rmxs,rmxe
    integer :: iret
!!$    character(len=FMAXVALLEN) :: rstr
    integer :: f_getstringValue, f_getIntValue, f_getRealValue
    logical :: tf
!!$    character(len=FMAXUNITLEN) unit_f
!!$    unit_f = 'hartree'

    icf = 0
    if(f_getIntValue(tag_no, iret) == 0) no = iret
    mway = SIMPLE
    if(f_getStringValue(tag_method,rstr,LOWER) == 0) then
       icf(1) = 1
       if(len_trim(rstr) <= len_trim(tag_simple) .and. len_trim(rstr) >= 1) then
          call strncmp0(tag_simple,trim(rstr),tf)
          if(tf) then
             mway = SIMPLE
             goto 1001
          end if
       end if
       if(len_trim(rstr) <= len_trim(tag_dfp) .and. len_trim(rstr) >= 1) then
          call strncmp0(tag_dfp,trim(rstr),tf)
          if(tf) then
             mway = DFP
             goto 1001
          end if
       end if
       if(len_trim(rstr) <= len_trim(tag_pulay) .and. len_trim(rstr) >= 1) then
          call strncmp0(tag_pulay,trim(rstr),tf)
          if(tf) then
             mway = PULAY
             goto 1001
          end if
       end if
       if(len_trim(rstr) == len_trim(tag_broyden) .and. len_trim(rstr) >= 1) then
          call strncmp0(tag_broyden, trim(rstr),tf)
          if(tf) then
             mway = BROYD1
             goto 1001
          end if
       end if
       if(len_trim(rstr) == len_trim(tag_broyden2)) then
          call strncmp0(tag_broyden2, trim(rstr),tf)
          if(tf) then
             mway = BROYD2
             goto 1001
          end if
       end if
       if(len_trim(rstr) == len_trim(tag_b2)) then
          call strncmp0(tag_b2,trim(rstr),tf)
          if(tf) then
             mway = BROYD2
             goto 1001
          end if
       end if
       if( len_trim(rstr) == 1 .and. trim(rstr) == "b") mway = BROYD1
1001   continue
    end if

    if(f_getRealValue(tag_rmxs,rmxs,'') == 0) icf(2) = 1
    if(f_getRealValue(tag_rmxe,rmxe,'') == 0) icf(3) = 1
    if(f_getIntValue(tag_itr,itr) == 0) icf(4) = 1
    if(f_getStringValue(tag_var,rstr,LOWER) == 0) then
       icf(5) = 1
       call strncmp0(tag_linear,trim(rstr),tf)
       if(tf) then
          var = varLINEAR
       else
          call strncmp0(tag_tanh,trim(rstr),tf)
          if(tf) then
             var = varTANH
          end if
       end if
    end if

    if(f_getIntValue(tag_prec,prec) == 0) icf(6) = 1
    if(f_getIntValue(tag_istr,istr) == 0) icf(7) = 1
    if(f_getIntValue(tag_nbmix,nbxmix) == 0) then
       icf(8) = 1
    else
       if(f_getIntValue(tag_nbxmix,nbxmix) == 0) then
          icf(8) = 1
       end if
    end if

    f_readchargemixing = f_getStringValue(tag_update,rstr,LOWER)
    if(f_readchargemixing == 0)  icf(9) = 1
    if(icf(9) /= 1) then
       f_readchargemixing = f_getStringValue(tag_hownew,rstr,LOWER)
       if(f_readchargemixing == 0) icf(9) = 1
    end if
    if(icf(9) == 1) then
       hownew = RENEW
       call strncmp0(tag_renew,trim(rstr),tf)
       if(tf) then
          hownew = RENEW
          goto 1002
       end if
       call strncmp0(tag_anew, trim(rstr),tf)
       if(tf) then
          hownew = ANEW
          goto 1002
       end if
    end if
1002 continue

  end function f_readchargemixing
#endif

  subroutine m_CtrlP_rd_printlevel(nfout)
    integer, intent(in) :: nfout
    integer :: iret, f_selectBlock
    real(kind=DP) :: dret
    integer :: f_selectParentBlock, f_selectTop, f_getIntValue, f_getRealValue
    integer :: f_getStringValue
    logical :: flag = .false.
    logical :: tag_is_found = .false.
    logical :: tf
!!$    character(len=FMAXVALLEN) :: rstr

    call set_printoutlevel_default(ipri)
    iret = f_selectTop()
    ! --- Printoutlevel ---
    if( f_selectBlock( tag_printoutlevel) == 0) flag = .true.
    if(.not.flag .and. f_selectBlock( tag_printlevel) == 0) flag = .true.
    if(flag) then
       call set_ipri(tag_ipriparadeb,ipriparadeb)
       if(.not.tag_is_found) call set_ipri(tag_ipriparallel_debug,ipriparadeb)

       if(mype == 0 .or. ipriparadeb /= 0) then
          printable = .true.
       else
          printable = .false.
       end if

       call set_ipri(tag_ipribase,ipri)
       call set_printoutlevel_default(ipri)
       call set_ipri(tag_ipritiming,ipritiming)
#ifndef _EMPIRICAL_
       call set_ipri(tag_iprisolver,iprisolver)
       call set_ipri(tag_iprievdff,iprievdff)
       call set_ipri(tag_iprirmm,iprirmm); rmm_printout = iprirmm
       call set_ipri(tag_ipridavidson,ipridavidson)
       call set_ipri(tag_iprimddavidson,iprimddavidson)
       call set_ipri(tag_iprimdkosugi,iprimdkosugi)
       call set_ipri(tag_ipriesm,ipriesm)
       call set_ipri(tag_iprifcp,iprifcp)
       call set_ipri(tag_iprivdw,iprivdw)
       call set_ipri(tag_iprirs,iprirs)
       call set_ipri(tag_iprirsb,iprirsb)
       call set_ipri(tag_iprihubbard,iprihubbard)
       call set_ipri(tag_iprisnl,iprisnl)
       call set_ipri(tag_ipriphig,ipriphig)
       call set_ipri(tag_ipripao,ipripao)
       call set_ipri(tag_ipriberry,ipriberry)
! === KT_add === 13.1R
       call set_ipri(tag_iprisym, iprisym )
! ============== 13.1R
       call set_ipri(tag_ipriphonon,ipriphonon)
       call set_ipri(tag_iprifef,iprifef)
       call set_ipri(tag_iprimatdiagon,iprimatdiagon)
       call set_ipri(tag_iprieigenvalue,iprieigenvalue)
       call set_ipri(tag_ipripulay,ipripulay)
       call set_ipri(tag_ipritotalcharge,ipritotalcharge)
       call set_ipri(tag_iprivlhxcq,iprivlhxcq)
       call set_ipri(tag_iprisubmat,iprisubmat)
       call set_ipri(tag_ipribetar,ipribetar)
       call set_ipri(tag_ipripaw,ipripaw)
       call set_ipri(tag_iprixc,iprixc)
       call set_ipri(tag_iprivloc,iprivloc)
       call set_ipri(tag_iprichargemixing,iprichargemixing)
       call set_ipri(tag_ipripositron,ipripositron)
       call set_ipri(tag_iprispg,ipri_spg)
       call set_ipri(tag_iprikp,ipri_kp)
       call set_ipri(tag_iprioccup,iprioccup)
       call set_ipri(tag_ipridos,ipridos)
       call set_ipri(tag_iprinegativecharge,iprinegativecharge)
       call set_ipri(tag_ipripp,ipripp)
       call set_ipri(tag_ipriwf,ipriwf)
       call set_ipri(tag_ipricoefwf,ipricoefwf)
       call set_ipri(tag_iprichargedensity,iprichargedensity)
#endif
       call set_ipri(tag_iprigdiis,iprigdiis)
       call set_ipri(tag_ipristrcfctr,ipristrcfctr)
       call set_ipri(tag_ipriparallel,ipriparallel)
       call set_ipri(tag_iprifftmap,iprifftmap)
       call set_ipri(tag_ipriinputfile,ipriinputfile)
       call set_ipri(tag_iprimd,iprimd)
       call set_ipri(tag_ipriforce,ipriforce)
       call set_ipri(tag_ipriexx,ipriexx)
       call set_ipri(tag_ipriblsize,ipriblsize)
       call set_ipri(tag_iprivelocity,iprivelocity)
       call set_ipri(tag_iprijobstatus,iprijobstatus)
       call set_ipri(tag_ipripredictor,ipripredictor)
       call set_ipri(tag_ipriunitcell,ipriunitcell)
       call set_ipri(tag_iprifire,iprifire)

! ---- for noncol. ---
       call set_ipri(tag_iprimagmom,iprimagmom)
       call set_ipri(tag_iprispinorb,iprispinorb)
! -------------------

! ==== KT_add ==== 13.0U2
! for TFW mixing --
       call set_ipri(tag_ipritfwfunc,ipritfwfunc)
! ================ 13.0U2

! === KT_add ==== 2014/09/24
       call set_ipri( tag_ipriepsilon,ipriepsilon )
! =============== 2014/09/24
       call set_ipri( tag_ipribravpos,ipribravpos )
       call set_ipri( tag_iprigap,iprigap )
       call set_ipri( tag_iprichgdefect,iprichgdefect )
       call set_ipri( tag_ipriorb_rot,ipriorb_rot )

       call multiply_printable_factor()
       if(ipriinputfile >= 2) call confirm_printoutlevel(nfout)

       if(f_getIntValue(tag_n_fermi_vicinity,iret)==0) n_fermi_vicinity = iret

       if( f_selectBlock( tag_timing_option) == 0) then
          if( f_getIntValue(tag_num_subroutines,iret) == 0) num_subroutines = iret
          if( f_getRealValue(tag_cputime_diff,dret,'') == 0) pcpudf = dret
          if(pcpudf < 0.d0) pcpudf = 0.d0
          if(pcpudf > 1.d0) pcpudf = 1.d0
          if( f_getIntValue(tag_sw_timing_2ndlevel,iret)==0) sw_timing_2ndlevel = iret

          if( f_getIntValue(tag_sw_firstlevel_only,iret)==0) sw_firstlevel_only = iret
          if( f_getIntValue(tag_sw_flatten,iret)==0)         sw_flatten = iret
          if( f_getIntValue(tag_sw_details,iret)==0)         sw_details = iret
          if( f_getIntValue(tag_measure_count_limit,iret)==0)     measure_count_limit = iret

          if( f_getIntValue(tag_statistics_in_parallel,iret) == 0) statistics_in_parallel = iret
          iret = f_selectParentBlock()
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !**  (timing_option) num_subroutines = ",i8)') num_subroutines
             write(nfout,'(" !**  (timing option) pcpudf = ",f8.4)') pcpudf
             write(nfout,'(" !**  (timing option) sw_timing_2ndlevel = ",i8)') sw_timing_2ndlevel

             write(nfout,'(" !**  (timing option) sw_flatten         = ",i8)') sw_flatten
             write(nfout,'(" !**  (timing option) sw_firstlevel_only = ",i8)') sw_firstlevel_only
             write(nfout,'(" !**  (timing option) sw_details         = ",i8)') sw_details
             write(nfout,'(" !**  (timing option) measure_count_limit = ",i8)') measure_count_limit
          end if
       end if

       if( f_selectBlock( tag_negativecharge_option) == 0) then
          if( f_getIntValue(tag_max_warnings,iret) == 0) max_warnings_negativecharge = iret
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !**  (negativecharge_option) max_warnings = ",i8)') max_warnings_negativecharge
          end if
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock( tag_jobstatus_option) == 0) then
          if( f_getIntValue(tag_jobstatus_series, iret) == 0) jobstatus_series = iret
          if( f_getStringValue(tag_jobstatus_format,rstr,LOWER) == 0) then
             call strncmp0(tag_tag_line,trim(rstr),tf)
             if(tf) then
                jobstatus_format = TAG_LINE
             endif
             if(.not.tf) then
                call strncmp0(tag_tag,trim(rstr),tf)
                if(tf) then
                   jobstatus_format = TAG_FORMAT
                end if
             end if
             if(.not.tf) then
                call strncmp0(tag_table,trim(rstr),tf)
                if(tf) then
                   jobstatus_format = TABLE
                end if
             end if
          end if
          if(ipriinputfile >= 2 .and. printable) then
             write(nfout,'(" !** (jobstatus_option) jobstatus_series = ",i4)') jobstatus_series
             write(nfout,'(" !** (jobstatus_option) jobstatus_format = ",i4)') jobstatus_format
          end if

          iret = f_selectParentBlock()
       end if

       iret = f_selectParentBlock()
    else
       if(mype == 0) then
          printable = .true.
       else
          printable = .false.
       end if
       call set_printoutlevel_default(ipri)
!!$       ipricoefwf = 0
       call multiply_printable_factor()
    end if
  contains
    subroutine set_ipri(tag,ipri)
      character(*),intent(in) :: tag
      integer, intent(inout) ::    ipri
      integer :: iret, f_getIntValue
      logical :: tf

      tag_is_found = .false.
      if(len_trim(tag) <= 4) return
      tf = f_getIntValue( tag, iret) == 0
      if(.not.tf) tf = f_getIntValue( tag(5:len_trim(tag)),iret) == 0
      if(tf) then
         tag_is_found = .true.
         ipri = iret
      end if
    end subroutine set_ipri
  end subroutine m_CtrlP_rd_printlevel

  subroutine set_printoutlevel_default(ipri)
    integer,intent(in) :: ipri
    ipritiming = ipri
#ifndef _EMPIRICAL_
    iprisolver = ipri
    iprievdff  = ipri
    iprirmm    = ipri
    ipridavidson = ipri
    iprimddavidson = ipri
    iprimdkosugi = ipri
    ipriesm = ipri
    iprifcp = ipri
    iprivdw = ipri
    iprirs = ipri
    iprirsb = ipri
    iprihubbard = ipri
    rmm_printout = ipri
    ipripulay  = ipri
    iprisnl    = ipri
    ipriphig   = ipri
    ipripao    = ipri
    ipriberry  = ipri
! ==== KT_add === 13.1R
    iprisym   = ipri
! ===============13.1R
    ipriphonon = ipri
    iprifef    = ipri
    iprimatdiagon = ipri
    iprivlhxcq = ipri
    iprisubmat = ipri
    ipribetar = ipri
    ipripaw   = ipri
    iprixc    = ipri
    iprieigenvalue = ipri
    ipritotalcharge = ipri
    iprivloc   = ipri
    iprichargemixing = ipri
    ipripositron = ipri
!!$    ipri_spg = ipri
    ipri_spg = 0
    ipri_kp  = ipri
    iprioccup = ipri
    ipridos = ipri
    iprinegativecharge = ipri
    ipripp  = ipri
    ipriwf  = ipri
    ipricoefwf  = ipri
    iprichargedensity = ipri
!!$    ipriekzaj  = ipri
#endif
    iprigdiis  = ipri
    ipristrcfctr = ipri
!!$    ipriparallel = ipri
    ipriparallel = -1
    iprifftmap = ipri
    ipriinputfile = ipri
    iprimd  = ipri
    ipriforce = ipri
    ipriexx = ipri
    ipriblsize = ipri
    iprivelocity = ipri
    iprijobstatus = ipri
    ipripredictor = ipri
    ipriunitcell = ipri
    iprifire = ipri
!!$    ipriparadeb = 0

! ------ for noncol --
    iprimagmom = ipri
    iprispinorb = ipri
! --------------------

! === KT_add === 2013/09/24
    ipriepsilon = ipri
! ============== 2013/09/24
    ipribravpos = ipri
    iprigap = ipri
    iprichgdefect = ipri
    ipriorb_rot = ipri

  end subroutine set_printoutlevel_default

  subroutine multiply_printable_factor()
    integer :: ifactor, ipri_0
    if(printable) then
       ifactor = 1
    else
       ifactor = 0
    end if

    ipri = ipri*ifactor
    ipritiming = ipritiming*ifactor
    ipritiming0 = ipritiming
    if(npes > 1) call mpi_bcast(ipritiming0,1,mpi_integer,0,MPI_CommGroup,ierr)
#ifndef _EMPIRICAL_
    iprisolver = iprisolver*ifactor
    iprievdff  = iprievdff*ifactor
    iprirmm    = iprirmm*ifactor
    ipridavidson = ipridavidson*ifactor
    if ( noncol ) then
       if (npes > 1) call mpi_bcast(ipridavidson,1,mpi_integer,0,MPI_CommGroup,ierr)
    endif
    iprimddavidson = iprimddavidson*ifactor
    iprimdkosugi = iprimdkosugi*ifactor
    ipriesm = ipriesm*ifactor
    iprifcp = iprifcp*ifactor
    iprivdw = iprivdw*ifactor
    iprirs = iprirs*ifactor
    iprirsb = iprirsb*ifactor
    iprihubbard = iprihubbard*ifactor
    rmm_printout = rmm_printout*ifactor
    ipripulay  = ipripulay*ifactor
    iprisnl    = iprisnl*ifactor
    ipriphig   = ipriphig*ifactor
    ipripao    = ipripao*ifactor
    ipriberry  = ipriberry*ifactor
! === KT_add ==== 13.1R
    iprisym = iprisym *ifactor
! =============== 13.1R
    ipriphonon = ipriphonon*ifactor
    iprifef    = iprifef*ifactor
    iprimatdiagon = iprimatdiagon*ifactor
    iprivlhxcq = iprivlhxcq*ifactor
    iprisubmat = iprisubmat*ifactor
    ipribetar = ipribetar*ifactor
    ipripaw   = ipripaw*ifactor
    iprixc    = iprixc*ifactor
    iprieigenvalue = iprieigenvalue*ifactor
    ipritotalcharge = ipritotalcharge*ifactor
    iprivloc   = iprivloc*ifactor
    iprichargemixing = iprichargemixing*ifactor
    ipripositron = ipripositron*ifactor
    ipri_spg = ipri_spg*ifactor
    ipri_kp  = ipri_kp*ifactor
    iprioccup = iprioccup*ifactor
    ipridos = ipridos*ifactor
    iprinegativecharge = iprinegativecharge*ifactor
    ipripp  = ipripp*ifactor
    ipriwf  = ipriwf*ifactor
    ipricoefwf  = ipricoefwf*ifactor
    iprichargedensity = iprichargedensity*ifactor

    iprinegativecharge0 = iprinegativecharge
    if(npes > 1) call mpi_bcast(iprinegativecharge0,1,mpi_integer,0,MPI_CommGroup,ierr)
#endif
    iprigdiis  = iprigdiis*ifactor
    ipristrcfctr = ipristrcfctr*ifactor
    ipriparallel = ipriparallel*ifactor
    iprifftmap = iprifftmap*ifactor
    ipriinputfile = ipriinputfile*ifactor
    iprimd  = iprimd*ifactor
    ipriforce = ipriforce*ifactor
    ipriexx   = ipriexx*ifactor
    ipriblsize = ipriblsize*ifactor
    iprivelocity = iprivelocity*ifactor
    iprijobstatus = iprijobstatus*ifactor
    ipripredictor = ipripredictor*ifactor
    ipriunitcell = ipriunitcell*ifactor
    iprifire = iprifire*ifactor

! --- for noncol --
    iprimagmom = iprimagmom *ifactor
    iprispinorb = iprispinorb *ifactor
! -----------------

! ==== KT_add === 2014/09/24
    ipriepsilon = ipriepsilon *ifactor
! =============== 2014/09/24
    ipribravpos = ipribravpos *ifactor
    iprigap = iprigap *ifactor
    iprichgdefect = iprichgdefect *ifactor
    if(npes > 1) call mpi_bcast( iprigap,1,mpi_integer,0,MPI_CommGroup,ierr )
    if(npes > 1) call mpi_bcast( ipriexx,1,mpi_integer,0,MPI_CommGroup,ierr )

    ipriorb_rot = ipriorb_rot *ifactor

  end subroutine multiply_printable_factor

  subroutine confirm_printoutlevel(nfout)
    integer, intent(in) :: nfout
    write(nfout,'(" !** confirming printoutlevels **")')
    write(nfout,'(" !** ipribase      = ",i4)') ipri
    write(nfout,'(" !** ipritiming    = ",i4)') ipritiming
#ifndef _EMPIRICAL_
    write(nfout,'(" !** iprisolver    = ",i4)') iprisolver
    write(nfout,'(" !** iprievdff     = ",i4)') iprievdff
    write(nfout,'(" !** iprirmm       = ",i4)') iprirmm
    write(nfout,'(" !** ipridavidson  = ",i4)') ipridavidson
    write(nfout,'(" !** iprimddavidson= ",i4)') iprimddavidson
    write(nfout,'(" !** iprimdkosugi  = ",i4)') iprimdkosugi
    write(nfout,'(" !** ipriesm       = ",i4)') ipriesm
    write(nfout,'(" !** iprifcp       = ",i4)') iprifcp
    write(nfout,'(" !** iprivdw       = ",i4)') iprivdw
    write(nfout,'(" !** iprirs        = ",i4)') iprirs
    write(nfout,'(" !** iprirsb       = ",i4)') iprirsb
    write(nfout,'(" !** iprihubbard   = ",i4)') iprihubbard
    write(nfout,'(" !** ipripulay     = ",i4)') ipripulay
    write(nfout,'(" !** iprisnl       = ",i4)') iprisnl
    write(nfout,'(" !** ipriphig      = ",i4)') ipriphig
    write(nfout,'(" !** ipripao       = ",i4)') ipripao
    write(nfout,'(" !** ipriberry     = ",i4)') ipriberry
! === KT_add === 13.1R
    write(nfout,'(" !** iprisym    = ",i4)') iprisym
! ============== 13.1R
    write(nfout,'(" !** ipriphonon    = ",i4)') ipriphonon
    write(nfout,'(" !** iprifef       = ",i4)') iprifef
    write(nfout,'(" !** iprimatdiagon = ",i4)') iprimatdiagon
    write(nfout,'(" !** ipri2vlhxcq   = ",i4)') iprivlhxcq
    write(nfout,'(" !** iprisubmat    = ",i4)') iprisubmat
    write(nfout,'(" !** ipribetar     = ",i4)') ipribetar
    write(nfout,'(" !** ipripaw       = ",i4)') ipripaw
    write(nfout,'(" !** iprixc        = ",i4)') iprixc
    write(nfout,'(" !** iprieigenvalue = ",i3)') iprieigenvalue
    write(nfout,'(" !** ipritotalcharge = ",i2)') ipritotalcharge
    write(nfout,'(" !** iprivloc      = ",i4)') iprivloc
    write(nfout,'(" !** iprichargemixing = ",i4)') iprichargemixing
    write(nfout,'(" !** ipripositron  = ",i4)') ipripositron
    write(nfout,'(" !** ipri_spg      = ",i4)') ipri_spg
    write(nfout,'(" !** ipri_kp       = ",i4)') ipri_kp
    write(nfout,'(" !** iprioccup     = ",i4)') iprioccup
    write(nfout,'(" !** ipridos       = ",i4)') ipridos
    write(nfout,'(" !** iprinegativecharge = ",i4)') iprinegativecharge
    write(nfout,'(" !** ipripp        = ",i4)') ipripp
    write(nfout,'(" !** ipriwf        = ",i4)') ipriwf
    write(nfout,'(" !** ipricoefwf    = ",i4)') ipricoefwf
    write(nfout,'(" !** iprichargedensity = ",i4)') iprichargedensity
!!$    write(nfout,'(" !** ipriekzaj     = ",i4)') ipriekzaj
#endif
    write(nfout,'(" !** iprigdiis     = ",i4)') iprigdiis
    write(nfout,'(" !** ipristrcfctr  = ",i4)') ipristrcfctr
    write(nfout,'(" !** ipriparallel  = ",i4)') ipriparallel
    write(nfout,'(" !** iprifftmap    = ",i4)') iprifftmap
    write(nfout,'(" !** ipriinputfile = ",i4)') ipri_kp
    write(nfout,'(" !** iprimd        = ",i4)') iprimd
    write(nfout,'(" !** ipriforce     = ",i4)') ipriforce
    write(nfout,'(" !** ipriexx       = ",i4)') ipriexx
    write(nfout,'(" !** ipriblsize    = ",i4)') ipriblsize
    write(nfout,'(" !** iprivelocity  = ",i4)') iprivelocity
    write(nfout,'(" !** ipriparadeb   = ",i4)') ipriparadeb
    write(nfout,'(" !** iprijobstatus = ",i4)') iprijobstatus
    write(nfout,'(" !** ipripredictor = ",i4)') ipripredictor
    write(nfout,'(" !** ipriunitcell  = ",i4)') ipriunitcell
    write(nfout,'(" !** iprifire      = ",i4)') iprifire

! ---- for noncol --
    write(nfout,'(" !** iprispinorb = ",i4)') iprispinorb
    write(nfout,'(" !** iprimagmom = ",i4)') iprimagmom
! -----------------

! === KT_add === 2014/09/24
    write(nfout,'(" !** ipriepsilon = ",i4)') ipriepsilon
! ============== 2014/09/24
    write(nfout,'(" !** ipribravpos = ",i4)') ipribravpos
    write(nfout,'(" !** iprigap = ",i4)') iprigap
    write(nfout,'(" !** iprichgdefect = ",i4)') iprichgdefect
    write(nfout,'(" !** ipriorb_rot = ",i4)') ipriorb_rot

  end subroutine confirm_printoutlevel

  subroutine m_CtrlP_rd_postproc(nfout)
    integer, intent(in) :: nfout
!!$    character(len=FMAXVALLEN) :: rstr
    integer :: iret, f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectTop, f_selectParentBlock
    real(kind=DP) :: dret
    character(3) :: str1

#ifndef _EMPIRICAL_
    logical :: tag_is_found, tf
    if(ipriinputfile >= 2) &
         & write(nfout,'(" !** << m_CtrlP_rd_postproc >>")')
    ! --- Postprocessing ---
    iret = f_selectTop()
    if( f_selectBlock( tag_postprocessing) == 0) then
       if(ipriinputfile >= 2) write(nfout,'(" !** -- tag_postprocessing --")')
       if(f_getIntValue(tag_frequency,iret)==0) postproc_frequency=iret
       if(ipriinputfile>=2) write(nfout,'(a,i2)') " !** -- portproc frequency ",postproc_frequency
       ! -- sw_band_symmetry_analysis--
       if( f_getIntValue( tag_sw_band_symmetry_analysis, iret) == 0) sw_band_symmetry_analysis = iret
       if(ipriinputfile>=1) write(nfout,'(a,i2)') " !** -- sw_band_symmetry_analysis ",sw_band_symmetry_analysis

!!$             write(nfout,'(" !** sw_band_symmetry_analysis = ",i6)') sw_band_symmetry_analysis
       ! -- ldos --
       if( f_selectBlock( tag_ldos) == 0) then
          if(printable) write(nfout,'(" !** tag_ldos is found")')
!!$          if( f_getStringValue(tag_ldos_hardpart_fft,rstr,LOWER)==0) call set_ldos_hardpart_fft(rstr) ! ->ldos_hardpart_fft
          tf = f_getStringValue(tag_dos_write_format,rstr,LOWER) == 0
          if(.not.tf) tf = f_getStringValue(tag_write_format,rstr,LOWER)==0
          if(tf) call set_write_format(rstr,dos_write_format) ! -> dos_write_format

          tf = f_getStringValue(tag_sw_checksum,rstr,LOWER) == 0
          if(tf) then
             if(printable) write(nfout,'(" !** tag_sw_checksum is found")')
          else
             if(printable) write(nfout,'(" !** tag_sw_checksum is not found")')
          end if
          if(tf) call set_sw_checksum(rstr,sw_checksum) ! -> sw_checksum

          if(ekmode == OFF) then
             if( f_getStringValue( tag_dos_method, rstr, LOWER) == 0) then
                call set_dos_method(rstr,ldos_method)  ! -> ldos_method
             else if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
                call set_dos_method(rstr,ldos_method)  ! -> ldos_method
             else if( way_of_smearing == TETRAHEDRON) then
                ldos_method = TETRAHEDRON
             else if( way_of_smearing == PARABOLIC) then
                ldos_method = Gauss_distrib_func

! ============================== KT_add ================= 13.0E
             else if( way_of_smearing == Fermi_Dirac ) then
                ldos_method = Fermi_Dirac
! ======================================================= 13.0E
             else
                ldos_method = Gauss_distrib_func
             end if
          else
             ldos_method = Gauss_distrib_func
          end if

          if( f_getIntValue( tag_hardpart_subroutine,iret) == 0) then
             hardpart_subroutine = iret
             if(printable) write(nfout,'(" !** hardart subroutine : ",i3)') hardpart_subroutine
          endif
          if( f_getIntValue(tag_sw_ac_mesh,iret) == 0) then
             sw_ac_mesh = iret
             if(printable) write(nfout,'(" !** sw_ac_mesh : ",i3)') sw_ac_mesh
          endif
          if( f_getIntValue(tag_acmesh_factor,iret) == 0) then
             acmesh_factor = iret
             if(printable) write(nfout,'(" !** ac_mesh_factor : ",i4)') acmesh_factor
          endif
          if( f_getIntValue( tag_sw_save_ldos_weight, iret) == 0) sw_save_ldos_weight = iret
          if(ekmode == ON .and. sw_save_ldos_weight == OFF) then
             sw_save_ldos_weight = ON
             if(printable) &
                  & write(nfout,'(" !** sw_save_ldos_weight has been reset ON from OFF")')
          end if

          if( f_getIntValue( tag_sw_cal_ldos, iret) == 0) sw_cal_ldos = iret

          tag_is_found = f_getIntValue( tag_sw_aldos, iret )==0
          if(.not.tag_is_found) tag_is_found = f_getIntValue( tag_sw_atomicdos,iret)==0
          if(tag_is_found) then
!!$          if ( f_getIntValue( tag_sw_aldos, iret)==0) THEN
             sw_aldos = iret
             if(sw_aldos == ON) then
                tag_is_found = f_selectBlock(tag_aldos)==0
                if(.not.tag_is_found) tag_is_found = f_selectBlock(tag_atomicdos)==0
                if(tag_is_found) then
!!$                if( f_selectBlock( tag_aldos) == 0) then
                   if( f_getRealValue(tag_crtdst, dret, "bohr") == 0) then
                      crtdst_is_given = .true.
                      crtdst_aldos = dret
                   end if
                   if( f_getIntValue( tag_naldos_from, iret ) == 0) naldos_from = iret
                   if( f_getIntValue( tag_naldos_to,   iret ) == 0) naldos_to   = iret

! atom-centered mesh
                   if( f_getIntValue(tag_sw_atom_centered_mesh,iret)==0 .or. &
                     & f_getIntValue(tag_sw_ac_mesh,iret) == 0) then
                      sw_ac_mesh = iret
                      if(printable) write(nfout,'(" !** sw_ac_mesh : ",i3)') sw_ac_mesh
                   endif
                   if( f_getIntValue(tag_atom_centered_mesh_factor,iret)==0 .or. &
                     & f_getIntValue(tag_acmesh_factor,iret) == 0) then
                      acmesh_factor = iret
                      if(printable) write(nfout,'(" !** ac_mesh_factor : ",i4)') acmesh_factor
                   endif
                   iret = f_selectParentBlock()
                end if
             end if
          end if
          if ( f_getIntValue( tag_sw_layerdos, iret) == 0) then
             sw_layerdos = iret
             if(sw_layerdos == ON) then
                if( f_selectBlock( tag_layerdos) == 0) then
                   if( f_getStringValue( tag_slicing_way, rstr, LOWER) == 0) call set_slicing_way(rstr)
                   tf = f_getRealValue(tag_deltaz,dret,"bohr") == 0
                   if(.not.tf) tf = f_getRealValue(tag_delta,dret,"bohr") == 0
                   if(tf) then
                      deltaz_winlay = dret
                   end if
                   if( f_getIntValue(tag_normal_axis,iret) == 0) normal_axis_winlay = iret
                   if( f_getRealValue(tag_crtdst,dret,"bohr") == 0) then
                      crtdst_is_given = .true.
                      crtdst_winlay = dret
                   end if
                   if( f_getIntValue(tag_integration_dimension,iret)==0) integration_dimension_winlay = iret
                   if(integration_dimension_winlay/=1 .and. integration_dimension_winlay/=3) &
                        & integration_dimension_winlay = 1
                   iret = f_selectParentBlock()
                end if
             end if
          end if
          if(f_getIntValue(tag_sw_rspace,iret)==0) sw_rspace_ldos = iret
          if(sw_aldos == ON .or. sw_layerdos == ON) sw_ldos = ON
          iret = f_selectParentBlock()
       end if

       ! --- dos ---
       if( f_selectBlock( tag_dos) == 0) then
          tag_is_found = f_getIntValue( tag_sw_dos, iret) == 0
          if(tag_is_found) then
             sw_dos = iret

             tf = f_getStringValue(tag_dos_write_format,rstr,LOWER) == 0
             if(.not.tf) tf = f_getStringValue(tag_write_format,rstr,LOWER)==0
             if(tf) call set_write_format(rstr,dos_write_format) ! -> dos_write_format

             if(sw_dos == ON) then
                if( f_getStringValue( tag_dos_method, rstr, LOWER) == 0) then
                   if(printable) write(nfout,'(" !** tag_dos_method is found")')
                   call set_dos_method(rstr,dos_method)  ! -> dos_method
                else if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
                   if(printable) write(nfout,'(" !** tag_method is found")')
                   call set_dos_method(rstr,dos_method)  ! -> dos_method
                else if( way_of_smearing == TETRAHEDRON) then
                   dos_method = TETRAHEDRON
                else if( way_of_smearing == PARABOLIC) then
                   dos_method = Gauss_distrib_func

! ======================== KT_add ========================= 13.0E
                else if( way_of_smearing == FERMI_DIRAC ) then
                   dos_method = Fermi_Dirac
! ========================================================= 13.0E
                else
                   dos_method = Gauss_distrib_func
                end if
                if(ipriinputfile >= 1) write(nfout,'(" !** ldos_method         = ",i6)') ldos_method
                if(ipriinputfile >= 2) then
                   write(nfout,'(" !** dos_method          = ",i6," (Gauss_distrib_func=",i2,",TETRAHEDRON=",i2,")")') &
                        & dos_method, Gauss_distrib_func, TETRAHEDRON
                end if
             end if
          else
             if( f_getIntValue( tag_sw_dos_gaussdistrib, iret) == 0) then
                sw_dos_Gaussdistrib = iret
                sw_dos = sw_dos_Gaussdistrib
                dos_method = Gauss_distrib_func
             end if
          end if
          if(f_getIntValue(tag_dos_subroutine,iret)==0) then
             dos_subroutine = iret
             if(dos_subroutine <= 2 .or.dos_subroutine >= 6) then
                dos_subroutine = 3
             end if
          end if
          tag_is_found = f_getRealValue( tag_deltaE_dos,dret,"hartree") == 0
          if(tag_is_found) deltaE_dos = dret
          if(.not.tag_is_found) then
             if( f_getRealValue( tag_deltaE_dos_GaussD, dret, "hartree") == 0) deltaE_dos = dret
          end if
          if( f_getRealValue( tag_variance_dos_GaussD, dret,'') == 0) then
             variance_dos_GaussD = dret
          else if( f_getRealValue( tag_variance_GaussD, dret, '') == 0) then
             variance_dos_GaussD = dret
          else if( f_getRealValue( tag_variance, dret, '') == 0) then
             variance_dos_GaussD = dret
          end if

          if( f_getRealValue( tag_dos_smearing_width, dret, "hartree" ) == 0 ) then
             dos_smearing_width = dret
             variance_dos_GaussD = dos_smearing_width**2
          else
             dos_smearing_width = sqrt( variance_dos_GaussD )
          endif

          if( f_getIntValue( tag_nwd_dos_window_width, iret) == 0) nwd_dos_window_width = iret
          if(ipriinputfile >= 1) then
             if(sw_ldos == ON) then
                write(nfout,'(" !** sw_ldos             = ON")')
                write(nfout,'(" !** ldos_method         = ",i6)') ldos_method
                write(nfout,'(" !** sw_aldos            = ",a3)') On_or_OFF(min(sw_aldos,1))
                write(nfout,'(" !** crtdst_aldos        = ",d12.4)') crtdst_aldos
                if(naldos_from /= 0) write(nfout,'(" !** naldos_from         = ",i6)') naldos_from
                if(naldos_to /= 0)   write(nfout,'(" !** naldos_to           = ",i6)') naldos_to
                write(nfout,'(" !** sw_layerdos         = ",a3)') On_or_OFF(min(sw_layerdos,1))
                write(nfout,'(" !** slicing_way_winlay  = ",i6)') slicing_way_winlay
                write(nfout,'(" !** deltaz_winlay       = ",d12.4)') deltaz_winlay
                write(nfout,'(" !** normal_axis_winlay  = ",i6)') normal_axis_winlay
                write(nfout,'(" !** crtdst_winlay       = ",d12.4)') crtdst_winlay
                write(nfout,'(" !** integration_dimension = ",i5)') integration_dimension_winlay
                write(nfout,'(" !** sw_save_ldos_weight = ",a3)') On_or_Off(min(1,sw_save_ldos_weight))
                write(nfout,'(" !** sw_cal_ldos         = ",i6)') sw_cal_ldos
                write(nfout,'(" !** sw_checksum         = ",i3," = ",a3)') sw_checksum,On_or_OFF(sw_checksum)
!!$                write(nfout,'(" !** ldos_hardpart_fft   = ",i6, " 1:FFT_REDUNDANT, 2;FFT_PARALLEL")') ldos_hardpart_fft
             end if

             if(printable) then
                write(nfout,'(" !** sw_dos              = ",a3)') On_or_Off(min(1,sw_dos))
                write(nfout,'(" !** dos_method          = ",i6," (Gauss_distrib_func=",i2,",TETRAHEDRON=",i2,")")') &
                     & dos_method, Gauss_distrib_func, TETRAHEDRON
                if(.not.tag_is_found) &
                     & write(nfout,'(" !** sw_dos_Gaussdistrib = ",i6)') sw_dos_Gaussdistrib
                write(nfout,'(" !** deltaE_dos          = ",d12.4)') deltaE_dos
                write(nfout,'(" !** variance_dos_GaussD = ",d12.4)') variance_dos_GaussD
                write(nfout,'(" !** dos_smearing_width = ",d12.4)') dos_smearing_width

                write(nfout,'(" !** nwd_dos_window_width = ",i8)') nwd_dos_window_width
                write(nfout,'(" !** dos_subroutine      = ",i5)') dos_subroutine
             end if
             if(set_write_format_count >= 1) then
                if(dos_write_format==WIDE)   write(nfout,'(" !** dos_write_format = wide (default)")')
                if(dos_write_format==NARROW) write(nfout,'(" !** dos_write_format = narrow")')
                if(dos_write_format==eVunit) write(nfout,'(" !** dos_write_format = eV unit")')
             end if
          end if

! =============================== added by K. Tagami ==================== 11.0
          if ( f_getIntValue( tag_calc_magmom_contrib, iret) == 0) then
             calc_dos_magmom_contrib = iret
             if (printable) then
                write(nfout,*) '****************************************** '
                write(nfout,*) '!** calc_dos_magmom_contrib is set to ', &
        &                       calc_dos_magmom_contrib
             endif
          else
             if(printable) then
                write(nfout,*) '!** calc_dos_magmom_contrib is set to default,  ', &
        &                       calc_dos_magmom_contrib
             endif
          endif
! ======================================================================= 11.0

          iret = f_selectParentBlock()
       end if

! -- lband --
       if ( f_selectBlock(tag_wf_local_decomposition) == 0 .or. &
            &  f_selectBlock(tag_lband) == 0 ) then

          tf = f_getStringValue(tag_sw_checksum,rstr,LOWER) == 0
          if(tf) then
             if(printable) write(nfout,'(" !** tag_sw_checksum is found")')
          else
             if(printable) write(nfout,'(" !** tag_sw_checksum is not found")')
          end if
          if(tf) call set_sw_checksum(rstr,sw_checksum) ! -> sw_checksum

          if( f_getIntValue( tag_hardpart_subroutine,iret) == 0) then
             hardpart_subroutine = iret
             if(printable) write(nfout,'(" !** hardart subroutine : ",i3)') hardpart_subroutine
          endif

          if ( f_getIntValue(tag_sw_calc_wf_atom_decomp,iret) == 0) then
             sw_calc_wf_atom_decomposition = iret
          endif
          write(nfout,*) "!** sw_calc_wf_atom_decomposition = ", &
                  &      sw_calc_wf_atom_decomposition

          if ( sw_calc_wf_atom_decomposition == ON ) then
             if ( f_selectBlock(tag_atom_decomposition) == 0 ) then
                if( f_getRealValue(tag_crtdst, dret, "bohr") == 0) then
                   crtdst_is_given = .true.
                   crtdst_aldos = dret
                end if

                if( f_getIntValue(tag_sw_wd_only_specified_atoms,iret)==0 ) then
                   sw_wd_only_specified_atoms = iret
                   write(nfout,*) "!** sw_wd_only_specified_atoms = ", iret
                endif

                ! atom-centered mesh
                if( f_getIntValue(tag_sw_atom_centered_mesh,iret)==0 .or. &
                     & f_getIntValue(tag_sw_ac_mesh,iret) == 0) then
                   sw_ac_mesh = iret
                   if(printable) write(nfout,'(" !** sw_ac_mesh : ",i3)') sw_ac_mesh
                endif
                if( f_getIntValue(tag_atom_centered_mesh_factor,iret)==0 .or. &
                     & f_getIntValue(tag_acmesh_factor,iret) == 0) then
                   acmesh_factor = iret
                   if(printable) write(nfout,'(" !** ac_mesh_factor : ",i4)') acmesh_factor
                endif
                iret = f_selectParentBlock()
             endif
          end if

          if ( f_getIntValue(tag_sw_calc_wf_layer_decomp,iret) == 0 ) then
             sw_calc_wf_layer_decomposition = iret
          endif
          write(nfout,*) "!** sw_calc_wf_layer_decomposition = ", &
               &      sw_calc_wf_layer_decomposition

          if ( sw_calc_wf_layer_decomposition == ON ) then
             if ( f_selectBlock(tag_layer_decomposition) == 0 ) then
                if ( f_getStringValue( tag_slicing_way, rstr, LOWER) == 0) &
                     &           call set_slicing_way(rstr)
                tf = f_getRealValue(tag_deltaz,dret,"bohr") == 0
                if(.not.tf) tf = f_getRealValue(tag_delta,dret,"bohr") == 0
                if(tf) then
                   deltaz_winlay = dret
                   write(nfout,*) "!** deltaz_winlay = ", dret
                end if
                if( f_getIntValue(tag_normal_axis,iret) == 0) then
                   normal_axis_winlay = iret
                   write(nfout,*) "!** normal_axis_winlay = ", iret
                endif
                if( f_getRealValue(tag_crtdst,dret,"bohr") == 0) then
                   crtdst_is_given = .true.
                   crtdst_winlay = dret
                   write(nfout,*) "!** crtdst = ", dret
                end if
                if( f_getIntValue(tag_integration_dimension,iret)==0) &
                     &                integration_dimension_winlay = iret
                if ( integration_dimension_winlay/=1 &
                     &        .and. integration_dimension_winlay/=3 ) &
                     &                integration_dimension_winlay = 1
                iret = f_selectParentBlock()
             end if
          end if

          if (f_getIntValue(tag_sw_rspace,iret)==0) sw_rspace_lband = iret
          if (sw_calc_wf_atom_decomposition == ON .or. &
               &  sw_calc_wf_layer_decomposition == ON ) sw_lband = ON
          iret = f_selectParentBlock()
       end if

! --- orbital population --
       if( f_selectBlock( tag_orbital_population) == 0) then
          if ( f_getIntValue( tag_orb_popu_method, iret ) == 0 ) then
             orb_popu_method = iret
             write(nfout,*) '!** eval_population_method is set to ', iret
          endif
          if ( f_getIntValue( tag_sw_write_orb_dens_mat_file, iret ) == 0 ) then
             sw_write_orb_dens_mat_file = iret
             write(nfout,*) '!** sw_write_orb_dens_mat_file is set to ', iret
          endif

          if ( f_getIntValue( tag_sw_diagonalize_population, iret ) == 0 ) then
             sw_diagonalize_population = iret
             write(nfout,*) '!** sw_diagonalize_population is set to ', iret
          endif
          if ( sw_diagonalize_population == ON ) then
             if ( icond/=FIXED_CHARGE .and. icond/=FIXED_CHARGE_CONTINUATION ) then
                sw_write_orb_dens_mat_file = ON
             endif
             if ( f_getIntValue( tag_sw_write_rotated_orbitals, iret ) == 0 ) then
                sw_write_rotated_orbitals = iret
                write(nfout,*) '!** sw_write_rotated_orbitals is set to ', iret
             endif
             if ( f_getStringValue( tag_population_diag_mode, rstr, LOWER) == 0) then
                if ( trim(rstr) == trim(tag_charge_density_matrix) ) then
                   population_diag_mode = DIAG_CHARGE_DENSITY_MATRIX
                else if ( trim(rstr) == trim(tag_spin_density_matrix) ) then
                   population_diag_mode = DIAG_SPIN_DENSITY_MATRIX
                else if ( trim(rstr) == trim(tag_ls) ) then
                   population_diag_mode = DIAG_LS

                else if ( trim(rstr) == trim(tag_local_point_group) ) then
                   population_diag_mode = LOCAL_POINT_GROUP
                else if ( trim(rstr) == trim(tag_local_double_point_group) ) then
                   population_diag_mode = LOCAL_DOUBLE_POINT_GROUP
                endif
                write(nfout,*) '!** population_diag_mode is set to ', &
                     &          population_diag_mode
             endif
             if ( f_getIntValue( tag_sw_read_orb_rot_mat_file, iret ) == 0 ) then
                sw_read_orb_rot_mat_file = iret
                write(nfout,*) '!** sw_read_orb_rot_mat_file is set to ', iret
             endif
             if ( f_getIntValue( tag_sw_calc_score_sigma_bond, iret ) == 0 ) then
                sw_calc_score_sigma_bond = iret
                write(nfout,*) '!** sw_calc_score_sigma_bond is set to ', iret
             endif
          endif
          iret = f_selectParentBlock()
       end if

! --- procar ---
       if( f_selectBlock( tag_procar ) == 0) then
          if ( f_getIntValue( tag_sw_write_procar_file, iret ) == 0 ) then
             sw_write_procar_file = iret
             write(nfout,*) '!** sw_write_procar_file is set to ', iret
          endif
          if ( sw_write_procar_file == ON ) then
             if( f_getIntValue(tag_sw_procar_full_bz, iret) == 0) then
                sw_procar_full_bz = iret
                write(nfout,*) '!** sw_procar_full_bz is set to ', &
                     &          sw_procar_full_bz
             endif
             if( f_getIntValue(tag_procar_save_memory_mode, iret) == 0) then
                procar_save_memory_mode = iret
                write(nfout,*) '!** procar_save_memory_mode is set to ', &
                     &          procar_save_memory_mode
             endif
             if ( f_getIntValue( tag_split_procar_file, iret ) == 0 ) then
                split_procar_file = iret
                write(nfout,*) '!** split_procat_file is set to ', &
                     &          split_procar_file
             endif
             if ( f_getIntValue( tag_num_procar_files_once, iret ) == 0 ) then
                num_procar_files_once = iret
                write(nfout,*) '!** num_procar_files_once is set to ', &
                     &          num_procar_files_once
             endif
             if ( f_getIntValue( tag_procar_sort_kpt, iret ) == 0 ) then
                procar_sort_kpt = iret
                write(nfout,*) '!** procat_sort_kpt is set to ', &
                     &          procar_sort_kpt
             endif
          endif
          iret = f_selectParentBlock()
       end if

! --- parity ---
       if( f_selectBlock( tag_parity ) == 0) then
          if ( f_getIntValue( tag_sw_write_parity_file, iret ) == 0 ) then
             sw_write_parity_file = iret
             write(nfout,*) '!** sw_write_parity_file is set to ', iret
          endif
          if ( f_getIntValue( tag_eval_parity_on_fftmesh, iret ) == 0 ) then
             eval_parity_on_fftmesh = iret
             write(nfout,*) '!** eval_parity_on_fftmesh set to ', iret
          endif
          iret = f_selectParentBlock()
       end if

! --- charge
       if( f_selectBlock( tag_charge) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_charge")')
          if( f_getIntValue( tag_sw_charge_rspace, iret) == 0) sw_charge_rspace = iret
          if( f_getIntValue( tag_sw_add_corecharge_rspace, iret) == 0) &
               &          sw_add_corecharge_rspace = iret
          if( f_getIntValue( tag_sw_read_corecharge_extra, iret) == 0) &
               &          sw_read_corecharge_extra_file = iret
          if( f_getIntValue( tag_eval_corecharge_on_Gspace, iret) == 0) &
               &          eval_corecharge_on_Gspace = iret
! ===== KT_add === 2014/06/07
          if ( f_getIntValue(tag_sw_spin_magmom_rspace,iret) == 0 ) then
             sw_spin_magmom_rspace = iret
          endif
! ================ 2014/06/07
          if( f_getStringValue( tag_filetype, rstr,LOWER) == 0) call set_charge_filetype(rstr,charge_filetype) !-> charge_filetype
          if( f_getStringValue( tag_title,rstr,NOCONV) == 0) then
             iret = len_trim(rstr)
             if(iret > LEN_TITLE) iret = LEN_TITLE
             charge_title(1:iret) = rstr(1:iret)
          end if
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_charge_rspace   = ",i6)') sw_charge_rspace
             write(nfout,'(" !** charge_filetype    = ",i6)') charge_filetype
             write(nfout,'(" !** charge_title    = ",a80)') charge_title
          end if
          if ( sw_add_corecharge_rspace == ON ) then
             write(nfout,'(" !** sw_add_corecharge_rspace   = ",i6)') &
                  &                  sw_add_corecharge_rspace
             write(nfout,'(" !** sw_read_corecharge_extra_file    = ",i6)') &
                  &                  sw_read_corecharge_extra_file
             write(nfout,'(" !** eval_corecharge_on_Gspace   = ",i6)') &
                  &                  eval_corecharge_on_Gspace
          endif
! ===== KT_add === 2014/06/07
          if(ipriinputfile >= 1) then
             write(nfout,*) '!** sw_spin_magmom_rspace = ', sw_spin_magmom_rspace
          endif
! ================ 2014/06/07

          if(f_selectBlock(tag_cube)==0)then
             call rd_charge_subset()
             iret = f_selectParentBlock()
          endif
          if( f_selectBlock( tag_partial_charge) == 0) then
             if( f_getIntValue( tag_sw_partial_charge, iret) == 0) sw_partial_charge = iret
             if(ipriinputfile >= 2) write(nfout,'(" !* tag_partial_charge")')
             if( f_getRealValue( tag_Erange_min, dret, "hartree") == 0) partial_charge_Emin = dret
             if( f_getRealValue( tag_Erange_max, dret, "hartree") == 0) partial_charge_Emax = dret
             if( f_getRealValue( tag_deltaE, dret, "hartree") == 0) then
                partial_charge_deltaE = dret
             else if( f_getRealValue( tag_deltaE_partial_charge, dret, "hartree") == 0) then
                partial_charge_deltaE = dret
             else if( f_getRealValue( tag_Erange_delta, dret, "hartree") == 0) then
                partial_charge_deltaE = dret
             end if
             if ( f_getStringValue( tag_filetype, rstr, LOWER) == 0) then
                call set_partial_charge_filetype(rstr)
             else if( f_getStringValue( tag_partial_charge_filetype, rstr, LOWER) == 0) then
                call set_partial_charge_filetype(rstr)
             else if( f_getStringValue( tag_outputfiletype, rstr, LOWER) == 0) then
                call set_partial_charge_filetype(rstr)
             end if

             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_partial_charge = ",i3)') sw_partial_charge
                write(nfout,'(" !** Erange_min (partial_charge_Emin) = ",f8.4," (hartree)")') partial_charge_Emin
                write(nfout,'(" !** Erange_max (partial_charge_Emax) = ",f8.4," (hartree)")') partial_charge_Emax
                write(nfout,'(" !** DeltaE (partial_charge_deltaE)   = ",f8.4," (hartree)")') partial_charge_deltaE
                write(nfout,'(" !** outputfiletype = ",i5)') partial_charge_filetype
             end if
             iret = f_selectParentBlock()
          end if
          if(sw_partial_charge == ON .and. sw_charge_rspace == OFF) then
             sw_charge_rspace = ON
             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_charge_rspace   = ",i6, " : This is changed due to sw_partial_charge")') sw_charge_rspace
             end if
          end if
          iret = f_selectParentBlock()
       end if

       if(f_selectBlock(tag_rsb_test)==0)then
          if(f_getIntValue(tag_sw_rsb_test,iret)==0) sw_rsb_test=iret
          if(f_getIntValue(tag_sw_valence_electrons_only,iret)==0) sw_valence_electrons_only = iret
          if(f_getIntValue(tag_bisect_by,iret)==0) bisect_by = iret
          if(f_getIntValue(tag_lmax,iret)==0) lmax_rsb = iret
          if( f_getRealValue( tag_eps_rsb, dret, "") == 0) eps_rsb = dret
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_wf) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_wf")')
          if( f_getIntValue( tag_sw_wf_rspace, iret) == 0) sw_wf_rspace = iret
          if( f_getStringValue( tag_filetype, rstr,LOWER) == 0) call set_wf_filetype(rstr) !-> charge_filetype
          if( f_getStringValue( tag_title,rstr,NOCONV) == 0) then
             iret = len_trim(rstr)
             if(iret > LEN_TITLE) iret = LEN_TITLE
             wf_title(1:iret) = rstr(1:iret)
          end if
          if( f_selectBlock( tag_eigenvalue) == 0) then
             if( f_getRealValue( tag_eigmin, dret, "hartree") == 0) eigmin_wf = dret
             if( f_getRealValue( tag_eigmax, dret, "hartree") == 0) eigmax_wf = dret
             iret = f_selectParentBlock()
          end if
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_wf_rspace   = ",i6)') sw_wf_rspace
             write(nfout,'(" !** wf_filetype    = ",i6)') wf_filetype
             write(nfout,'(" !** wf_title    = ",/a80)') wf_title
! === DEBUG by tkato 2015/03/19 ================================================
!            write(nfout,'(" !** eigmin = ",f8.4," (hartree)")') eigmin_wf
!            write(nfout,'(" !** eigmax = ",f8.4," (hartree)")') eigmax_wf
             write(nfout,'(" !** eigmin = ",e12.4," (hartree)")') eigmin_wf
             write(nfout,'(" !** eigmax = ",e12.4," (hartree)")') eigmax_wf
! ==============================================================================
          end if
          iret = f_selectParentBlock()
       end if
       if( f_selectBlock( tag_rwf2) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_rwf2")')
          if( f_getIntValue( tag_sw_rwf2, iret) == 0) sw_rwf2 = iret
          if( f_getIntValue( tag_ib, iret) == 0) ib_rwf2 = iret
          if( f_getIntValue( tag_ik, iret) == 0) ik_rwf2 = iret
          if( f_getIntValue( tag_nr, iret) == 0) nr = iret
          if( f_getRealValue( tag_rmax, dret, "bohr") == 0) rmax = dret
          if( f_selectBlock( tag_center) == 0) then
             if( f_getRealValue( tag_rx, dret, "bohr") == 0) center(1) = dret
             if( f_getRealValue( tag_ry, dret, "bohr") == 0) center(2) = dret
             if( f_getRealValue( tag_rz, dret, "bohr") == 0) center(3) = dret
             iret = f_selectParentBlock()
          end if
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_rwf2   = ",i6)') sw_rwf2
             write(nfout,'(" !** ib   = ",i6)') ib_rwf2
             write(nfout,'(" !** ik   = ",i6)') ik_rwf2
             write(nfout,'(" !** nr   = ",i6)') nr
             write(nfout,'(" !** rmax = ",f8.4," (bohr)")') rmax
             write(nfout,'(" !** center = ",3f10.4," (bohr)")') center
          end if
          iret = f_selectParentBlock()
       end if

       if ( f_selectBlock( tag_wf_orb_projection ) == 0 ) then
          if ( f_getIntValue( tag_sw_calc_wf_orb_projection, iret ) == 0 ) then
             sw_calc_wf_orb_projection = iret
          endif
          if ( f_getIntValue( tag_use_rotated_compri, iret ) == 0 ) then
              use_rotated_compri = iret
          endif
          if ( f_getIntValue( tag_wf_orb_proj_print_format, iret ) == 0 ) then
             if ( iret < 0 .or. iret > 1 ) iret = 0
             wf_orb_proj_print_format = iret
          endif
          if (ipriinputfile >= 1) then
             write(nfout,'(A,I6)') " !** sw_calc_wf_orb_projection = ", &
                  &                   sw_calc_wf_orb_projection
!             write(nfout,'(A,i6,A)') "!** wf_orb_proj_print_format = ", &
!                  &        wf_orb_proj_print_format, " ( 0: {l m t}, 1: {j l mj t} )"
          endif
!          if ( sw_calc_wf_orb_projection == ON ) then
!             orb_popu_method = 2
!             write(nfout,'(A)') "!** orb_popu_method is forced to be 2"
!          endif
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_wf_squared ) == 0 ) then
          if( f_getIntValue( tag_sw_wf_squared_rspace, iret) == 0) &
               &                             sw_wf_squared_rspace = iret
          if( f_getIntValue( tag_ik, iret) == 0) ik_wf_squared = iret
          if( f_getIntValue( tag_ib, iret) == 0) ib1_wf_squared = iret

          if( f_getIntValue( tag_ib1, iret) == 0) ib1_wf_squared = iret
          if( f_getIntValue( tag_ib2, iret) == 0) then
             ib2_wf_squared = iret
          else
             ib2_wf_squared = ib1_wf_squared
          endif

          if( f_getIntValue( tag_sw_wf_integ_moment, iret) == 0) &
               &                             sw_wf_integ_moment = iret

          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_wf_squared_rspace   = ",i6)') sw_wf_squared_rspace
             if ( sw_wf_squared_rspace == ON ) then
                write(nfout,'(" !** ik   = ",i6)') ik_wf_squared
                write(nfout,'(" !** ib1  = ",i6)') ib1_wf_squared
                write(nfout,'(" !** ib2  = ",i6)') ib2_wf_squared
             endif
             write(nfout,'(" !** sw_wf_integ_moment = ",i6)') sw_wf_integ_moment
          endif

          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_band_unfolding ) == 0) then
          if( f_getIntValue( tag_sw_band_unfolding, iret) == 0) then
             sw_band_unfolding = iret
             write(nfout,*) "!** sw_band_unfolding is ", iret
          endif
          if( f_getIntValue( tag_prepare_masked_compri, iret) == 0) then
             prepare_masked_compri = iret
             write(nfout,*) "!** prepare_masked_compri is ", iret
          endif
          if( f_getRealValue( tag_tol_gvec_matching, dret, "" ) == 0) then
             tolerance_Gvec_matching = dret
             write(nfout,*) "!** tolerance_Gvec_matching is ", dret
          endif
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_msb ) == 0) then
          if( f_getIntValue( tag_sw_calc_contact_density, iret) == 0) then
             sw_calc_contact_density = iret
             write(nfout,*) "!** sw_calc_contact_density is ", iret
          endif
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_elf) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_elf")')
          if( f_getIntValue( tag_sw_elf, iret) == 0) sw_elf = iret
          if( f_getStringValue( tag_filetype, rstr,LOWER) == 0) call set_elf_filetype(rstr) !-> elf_filetype
          if( f_getStringValue( tag_title,rstr,NOCONV) == 0) then
             iret = len_trim(rstr)
             if(iret > LEN_TITLE) iret = LEN_TITLE
             elf_title(1:iret) = rstr(1:iret)
          end if
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_elf   = ",i6)') sw_elf
             write(nfout,'(" !** elf_filetype    = ",i6)') elf_filetype
             write(nfout,'(" !** elf_title    = ",/a80)') elf_title
          end if
          iret = f_selectParentBlock()
       end if
       if( f_selectBlock( tag_polarization) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_polarization")')
          if( f_getIntValue( tag_sw_bp_property, iret) == 0) sw_bp_property = iret
          if( sw_bp_property == ON ) then
             if( f_getStringValue( tag_property, rstr, LOWER) == 0 ) then
                call set_polar_prop(rstr,polar_prop) ! -> polar_prop
             end if
             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_bp_property      = ",i6)') sw_bp_property
                write(nfout,'(" !** polar_prop          = ",i6)') polar_prop
             end if
! === KT_add ==== 2014/06/30
             if( f_getIntValue( tag_sw_chk_sumrule_born_charge, iret) == 0 ) then
                sw_check_sumrule_born_charge = iret
                write(nfout,*) '!** sw_check_sumrule_born_charge =', iret
             endif
! =============== 2014/06/30
          end if
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock( tag_STM ) == 0 .or. f_selectBlock(tag_workfunc) == 0 ) then
          if(ipriinputfile >= 2) write(nfout,'(" !** tag_STM")')
          tag_is_found = f_getIntValue( tag_sw_STM, iret ) == 0
          if(.not.tag_is_found) tag_is_found = f_getIntValue( tag_sw_fine_STM, iret ) == 0
          if(.not.tag_is_found) tag_is_found = f_getIntValue( tag_sw_fine_STM_simulation, iret ) == 0
          if(.not.tag_is_found) tag_is_found = f_getIntValue( tag_sw_STM_images, iret ) == 0
          if(.not.tag_is_found) tag_is_found = f_getIntValue( tag_sw_workfunc, iret ) == 0
          if(tag_is_found) sw_fine_STM_simulation = iret
          tag_is_found = f_getIntValue( tag_sw_deficit_charge, iret ) == 0
          if(tag_is_found) sw_deficit_charge = iret
          if(f_getIntValue(tag_sw_add_xc_to_vloc,iret)==0) sw_add_xc_to_vloc = iret
          if(f_getIntValue(tag_sw_xc_only,iret)==0) sw_xc_only = iret
          if(f_getIntValue(tag_z_axis,iret)==0) z_axis = iret
          if(f_getRealValue(tag_connect_from, dret, "bohr")==0) connect_from = dret
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_STM              = ",i6)') sw_fine_STM_simulation
             if(sw_fine_STM_simulation == 1) then
                write(nfout,'(" !** sw_deficit_charge   = ",i6)') sw_deficit_charge
             endif
          endif
          iret = f_selectParentBlock()
       end if

       if( f_selectBlock( tag_dipole ) == 0 .and. sw_dipole_correction == OFF) then
          if(ipriinputfile >= 2) write(nfout,'(" !** tag_dipole")')
          if(f_getIntValue( tag_sw_dipole, iret ) == 0) sw_dipole = iret
          if(f_getIntValue( tag_sw_layered, iret ) == 0) sw_layered = iret
          if(f_getIntValue( tag_direction, iret ) == 0) idir_dip = iret
          if(f_getIntValue( tag_division, iret ) == 0) ndiv_dip = iret
          if(f_getRealValue( tag_width, dret, "bohr") == 0) width_dip = dret
          if( f_selectBlock( tag_vacuum ) == 0) then
             if(f_getRealValue( tag_rx, dret, "") == 0) rvac(1) = dret
             if(f_getRealValue( tag_ry, dret, "") == 0) rvac(2) = dret
             if(f_getRealValue( tag_rz, dret, "") == 0) rvac(3) = dret
             iret = f_selectParentBlock()
          end if
          iret = f_selectParentBlock()
       end if

       if ( f_selectBlock( tag_fermi_surface ) == 0 ) then
          if ( f_getIntValue( tag_sw_write_bxsf_file, iret ) == 0 ) then
             sw_write_bxsf_file = iret
             write(nfout,*) "!*  sw_write_bxsf_file is ", iret
          endif
          iret = f_selectParentBlock()
       endif

       if ( f_selectBlock( tag_spinorbit ) == 0 ) then
          if ( f_getIntValue( tag_sw_write_soi_on_atoms, iret ) == 0 ) then
             sw_write_soi_on_atoms = iret
             write(nfout,*) "!*  sw_write_soi_on_atoms is ", iret
          endif
          iret = f_selectParentBlock()
       endif

       if( f_selectBlock( tag_wannier ) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !** tag_wannier")')
          if(f_getIntValue( tag_sw_wannier90, iret ) == 0) sw_wannier90 = iret
          if(sw_wannier90 == OFF) then
             if(f_getIntValue( tag_sw_wannier, iret ) == 0) sw_wannier = iret
             if(f_getIntValue( tag_sw_random_wannier, iret ) == 0) sw_random_wannier = iret
             if(f_getIntValue( tag_sw_potential, iret ) == 0) sw_potential_wannier = iret
             if(f_getIntValue( tag_sw_continue, iret ) == 0) sw_continue_wannier = iret
             if(f_getRealValue( tag_eps_grad, dret, "") == 0) eps_wan = dret
             if(f_getRealValue( tag_dt, dret, "") == 0) dt_wan = dret
             if( f_getIntValue( tag_max_iteration, iret ) == 0 ) max_iter_wan = iret
             if( f_getStringValue( tag_filetype, rstr,LOWER) == 0) call set_wannier_filetype(rstr)
             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_wannier          = ",i6)') sw_wannier
                write(nfout,'(" !** sw_random_wannier   = ",i6)') sw_random_wannier
                write(nfout,'(" !** sw_potential_wannier= ",i6)') sw_potential_wannier
                write(nfout,'(" !** sw_continue_wannier = ",i6)') sw_continue_wannier
                write(nfout,'(" !** eps_grad            = ",e12.5)') eps_wan
                write(nfout,'(" !** dt_wan              = ",f10.5)') dt_wan
                write(nfout,'(" !** max_iter_wan        = ",i6)') max_iter_wan
                write(nfout,'(" !** wannier_filetype    = ",i6)') wannier_filetype
             end if
          else
             if( f_getStringValue( tag_seedname,rstr,NOCONV) == 0) then
                iret = len_trim(rstr)
                if(iret > LEN_TITLE) iret = LEN_TITLE
                wan90_seedname = ""
                wan90_seedname(1:iret) = rstr(1:iret)
             end if
             if( f_getIntValue( tag_nb_wan90, iret ) == 0 ) nb_wan90 = iret
             if( f_getIntValue( tag_sw_use_hardpart_wan90, iret ) == 0) then
                sw_use_hardpart_wan90 = iret
             endif
             if( f_getIntValue( tag_sw_write_unk_file, iret ) == 0) then
                sw_write_unk_file = iret
             endif
             if( f_getIntValue( tag_spin_component_wan90, iret ) == 0 ) then
                spin_component_wan90 = iret
             endif

             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_wannier90        = ",i6)') sw_wannier90
                write(nfout,'(" !** seedname            = ",a80)') wan90_seedname
                write(nfout,'(" !** nb_wan90            = ",i6)') nb_wan90
                write(nfout,'(" !** sw_use_hardpart_wan90 = ",i6)') sw_use_hardpart_wan90
                write(nfout,'(" !** sw_write_unk_file = ",i6)') sw_write_unk_file

                if ( .not. noncol .and. nspin == 2 ) then
                   write(nfout,'(" !** spin_component_wan90            = ",i6)') &
                        &              spin_component_wan90
                endif
             end if
          end if
          iret = f_selectParentBlock()
       end if

! ========= experimental ===
       if ( f_selectBlock( tag_crystal_field ) == 0 ) then
          if ( f_getIntValue( tag_sw_print_crys_field_param, iret) == 0 ) then
             sw_print_crystal_field_param = iret
          endif
          iret = f_selectParentBlock()
       endif

! ========= KT_add ============ 13.0S
       if ( f_selectBlock( tag_corelevels ) == 0 ) then
          if ( f_getIntValue( tag_sw_calc_core_energy, iret) == 0 ) then
             sw_calc_core_energy = iret
          endif
! === KT_add === 2014/08/08
          if ( sw_calc_core_energy == ON ) then
             if ( f_selectBlock( tag_corehole ) == 0 ) then
                if ( f_getIntValue( tag_atom_id, iret ) == 0 ) then
                   atom_with_corehole = iret
                endif
                if ( f_getStringValue( tag_orbital, rstr, LOWER ) == 0 ) then
                   str1 = trim(rstr)
                   if ( str1(1:1) == "1" ) qnum_n_corehole = 1
                   if ( str1(1:1) == "2" ) qnum_n_corehole = 2
                   if ( str1(1:1) == "3" ) qnum_n_corehole = 3
                   if ( str1(1:1) == "4" ) qnum_n_corehole = 4
                   if ( str1(1:1) == "5" ) qnum_n_corehole = 5
                   if ( str1(1:1) == "6" ) qnum_n_corehole = 6

                   if ( str1(2:2) == "s" .or. str1(2:2) == "s" ) qnum_l_corehole = 0
                   if ( str1(2:2) == "p" .or. str1(2:2) == "p" ) qnum_l_corehole = 1
                   if ( str1(2:2) == "d" .or. str1(2:2) == "d" ) qnum_l_corehole = 2
                   if ( str1(2:2) == "f" .or. str1(2:2) == "f" ) qnum_l_corehole = 3
                endif
                iret = f_selectParentBlock()
             endif

             iret = f_getStringValue( tag_eval_core_level_splitting,rstr,LOWER )
             if ( rstr == tag_core_level_splitting_pawpot ) then
                eval_core_level_splitting = ByPawPot
                write(nfout,*) '!** CoreLevel splliting is evaluated by paw pot.'
             else if ( rstr == tag_core_level_splitting_frompp ) then
                eval_core_level_splitting = ReadFromPP
                write(nfout,*) '!** CoreLevel splliting is read from PP'
             endif
          endif
! ============= 2014/08/08
          iret = f_selectParentBlock()
       endif
! ============================ 13.0S

       if( f_selectBlock( tag_rttddft) == 0) then
          if(ipriinputfile >= 2) write(nfout,'(" !*  tag_rttddft")')
          if(f_getIntValue( tag_sw_rttddft, iret) == 0) sw_rttddft = iret
          if(f_getIntValue( tag_time_step_max, iret) == 0) time_step_max = iret
          if(f_getRealValue( tag_time_step_delta, dret, "au_time") == 0) time_step_delta = dret
          if(f_getIntValue( tag_propagator_method, iret) == 0) propagator_method = iret
          if(f_getIntValue( tag_propagator_order, iret) == 0) propagator_order = iret
          if(f_getIntValue( tag_ext_ie_elec, iret) == 0) ext_ie_elec = iret
          if(f_getIntValue( tag_ext_ie_hole, iret) == 0) ext_ie_hole = iret
          if(f_getRealValue( tag_ext_pulse_epsilon, dret, '') == 0) ext_pulse_epsilon = dret
          if(f_getRealValue( tag_ext_pulse_kx, dret, '') == 0) ext_pulse_kx = dret
          if(f_getRealValue( tag_ext_pulse_ky, dret, '') == 0) ext_pulse_ky = dret
          if(f_getRealValue( tag_ext_pulse_kz, dret, '') == 0) ext_pulse_kz = dret
          if(ipriinputfile >= 1) then
             write(nfout,'(" !** sw_rttddft          = ",i6)') sw_rttddft
             write(nfout,'(" !** time_step_max       = ",i6)') time_step_max
             write(nfout,'(" !** time_step_delta     = ",e12.5)') time_step_delta
             write(nfout,'(" !** propagator_method   = ",i6)') propagator_method
             write(nfout,'(" !** propagator_order    = ",i6)') propagator_order
             write(nfout,'(" !** ext_ie_elec         = ",i6)') ext_ie_elec
             write(nfout,'(" !** ext_ie_hole         = ",i6)') ext_ie_hole
             write(nfout,'(" !** ext_pulse_epsilon   = ",e12.5)') ext_pulse_epsilon
             write(nfout,'(" !** ext_pulse_kx        = ",e12.5)') ext_pulse_kx
             write(nfout,'(" !** ext_pulse_ky        = ",e12.5)') ext_pulse_ky
             write(nfout,'(" !** ext_pulse_kz        = ",e12.5)') ext_pulse_kz
          endif
          iret = f_selectParentBlock()
       end if
       if(f_selectBlock(tag_boltztrap) ==0) then
          if(f_getIntValue(tag_boltztrap_output,iret) == 0) sw_boltztrap = iret
          if(sw_boltztrap == ON)then
              !!boltztrap_prefix = 'phase0'
              call set_default_btprefix()
              if(f_getStringValue(tag_boltztrap_prefix, rstr, NOCONV) == 0)then
                  boltztrap_prefix = rstr
              endif
              boltztrap_header = boltztrap_prefix
              if(f_getStringValue(tag_boltztrap_header, rstr, NOCONV) == 0)then
                  boltztrap_header = rstr
              endif
              if(f_getIntValue(tag_boltztrap_version,iret) == 0) boltztrap_version = iret
              if(boltztrap_version /= 1 .and. boltztrap_version /= 2) then
                  write(nfout,'(a,i5,a)') ' !** invalid boltztrap version ',boltztrap_version, &
                  & 'using the default value of 1'
                  boltztrap_version = 1
              endif
              write(nfout,'(a)')    " !** sw_boltztrap      = ON"
              write(nfout,'(a)')    " !** prefix            = "//trim(boltztrap_prefix)
              write(nfout,'(a)')    " !** header            = "//trim(boltztrap_header)
              write(nfout,'(a,i5)') ' !** boltztrap version = ',boltztrap_version
          endif
          iret = f_selectParentBlock()
       endif
#ifdef __EDA__
       if(f_selectBlock(tag_eda) == 0) then
          if(f_getIntValue(tag_sw_eda, iret) == 0) sw_eda = iret
          if (sw_eda == ON) write(nfout,'(a)')    " !** sw_eda      = ON"
          iret = f_selectParentBlock()
       endif
#endif

       iret = f_selectParentBlock()
    end if
  contains
    subroutine set_sw_checksum(rstr,sw_checksum)
      character(len=FMAXVALLEN), intent(in) :: rstr
      integer, intent(out) :: sw_checksum
      logical :: tf
      call strncmp2(rstr,FMAXVALLEN,'on',2,tf)
!!$      if(.not.tf) call strncmp0(trim(rstr),tag_on,tf)
      if(tf) then
         sw_checksum = ON
         goto 1001
      end if
      call strncmp2(rstr,FMAXVALLEN,'off',3,tf)
!!$      if(.not.tf) call strncmp0(trim(rstr),tag_off,tf)
      if(tf) then
         sw_checksum = OFF
         goto 1001
      end if
!!$      write(nfout,'("sw_checksum = ",a40, " no set")') rstr(1:40)
1001  continue
      write(nfout,'(" !** sw_checksum = ",a10, " = ",i3)') rstr(1:10),sw_checksum
    end subroutine set_sw_checksum

    subroutine set_write_format(rstr,dos_write_format)
      character(len=FMAXVALLEN),intent(in) :: rstr
      integer, intent(out) :: dos_write_format
      logical :: tf
      call strncmp2(rstr,FMAXVALLEN, tag_wide,len(tag_wide),tf)
      if(tf) then
         dos_write_format = WIDE
         set_write_format_count = set_write_format_count+1
         goto 1001
      end if
      call strncmp2(rstr,FMAXVALLEN, tag_narrow, len(tag_narrow),tf)
      if(tf) then
         dos_write_format = NARROW
         set_write_format_count = set_write_format_count+1
         goto 1001
      end if
      call strncmp2(rstr,FMAXVALLEN, tag_eVunit, len(tag_eVunit),tf)
      if(.not.tf) call strncmp2(rstr,FMAXVALLEN,tag_eV,len(tag_eV),tf)
      if(tf) then
         dos_write_format = eVunit
         set_write_format_count = set_write_format_count+1
         goto 1001
      end if
      if(ipriinputfile>=1) write(nfout,'(" !** dos_write_format = ",a40)') rstr(1:40)
1001  continue
      if(ipriinputfile >=  1) then
         write(nfout,'(" !** dos_write_format = ",i6, " set_time = ",i3)') dos_write_format,set_write_format_count
      end if
    end subroutine set_write_format
!!$    subroutine set_ldos_hardpart_fft(rstr)
!!$      character(len=FMAXVALLEN),intent(in) :: rstr
!!$      logical :: tf
!!$      call strncmp2(rstr,FMAXVALLEN, tag_redundant,len(tag_redundant),tf)
!!$      if(.not.tf) call strncmp2(rstr,FMAXVALLEN,'1',1,tf)
!!$      if(tf) then
!!$         ldos_hardpart_fft = FFT_REDUNDANT
!!$         goto 1001
!!$      end if
!!$      call strncmp2(rstr,FMAXVALLEN, tag_parallel,len(tag_parallel),tf)
!!$      if(.not.tf) call strncmp2(rstr,FMAXVALLEN,'2',1,tf)
!!$      if(tf) then
!!$         ldos_hardpart_fft = FFT_PARALLEL
!!$         goto 1001
!!$      end if
!!$1001  continue
!!$    end subroutine set_ldos_hardpart_fft

    subroutine set_slicing_way(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp0(tag_regular_intervals,trim(rstr),tf)
      if(tf) then
         slicing_way_winlay = REGULAR_INTERVALS
         goto 1001
      end if
      call strncmp0(tag_by_atomic_positions,trim(rstr),tf)
      if(tf) then
         slicing_way_winlay = BY_ATOMIC_POSITIONS
         goto 1001
      end if
1001  continue
    end subroutine set_slicing_way

    subroutine set_dos_method(rstr,dos_method_t)
      character(len=FMAXVALLEN),intent(in) :: rstr
      integer, intent(out) :: dos_method_t
      logical :: tf

      call strncmp0(tag_gaussdistrib, trim(rstr), tf)
      if(.not.tf) call strncmp0(tag_gaussiandistrib, trim(rstr), tf)
      if(.not.tf) call strncmp0(tag_gaussian, trim(rstr), tf)
      if(tf) then
         dos_method_t = Gauss_distrib_func
         goto 1001
      end if
      call strncmp0(tag_tetrahedral, trim(rstr), tf)
      if(.not.tf) &
           & call strncmp0(tag_tetrahedron, trim(rstr), tf)
      if(tf) then
         dos_method_t = TETRAHEDRON
         goto 1001
      end if
1001  continue
!!$      if(dos_method_t == TETRAHEDRON .and. way_of_smearing /= TETRAHEDRON) then
!!$         dos_method_t = Gauss_distrib_func
!!$      end if
    end subroutine set_dos_method

    subroutine set_partial_charge_filetype(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_individual, len(tag_individual),tf)
      if(tf) then
         partial_charge_filetype = SEPARATE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_separate, len(tag_separate),tf)
      if(tf) then
         partial_charge_filetype = SEPARATE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_integrated, len(tag_integrated),tf)
      if(tf) then
         partial_charge_filetype = INTEGRATED
         goto 1001
      end if

1001  continue
    end subroutine set_partial_charge_filetype

    subroutine set_wf_filetype(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_cube, len(tag_cube),tf)
      if(tf) then
         wf_filetype = CUBE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_vtk, len(tag_vtk),tf)
      if(tf) then
         wf_filetype = VTK
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_binary, len(tag_binary),tf)
      if(tf) then
         wf_filetype = BINARY
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_density_only, len(tag_density_only),tf)
      if(tf) then
         wf_filetype = DENSITY_ONLY
         goto 1001
      end if
1001  continue
    end subroutine set_wf_filetype

    subroutine set_wannier_filetype(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_cube, len(tag_cube),tf)
      if(tf) then
         wannier_filetype = CUBE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_vtk, len(tag_vtk),tf)
      if(tf) then
         wannier_filetype = VTK
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_binary, len(tag_binary),tf)
      if(tf) then
         wannier_filetype = BINARY
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_density_only, len(tag_density_only),tf)
      if(tf) then
         wannier_filetype = DENSITY_ONLY
         goto 1001
      end if
1001  continue
    end subroutine set_wannier_filetype

    subroutine set_elf_filetype(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_cube, len(tag_cube),tf)
      if(tf) then
         elf_filetype = CUBE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_density_only, len(tag_density_only),tf)
      if(tf) then
         elf_filetype = DENSITY_ONLY
         goto 1001
      end if
1001  continue
    end subroutine set_elf_filetype

    subroutine rd_charge_subset()
       real(kind=DP) :: tmpmin,tmpmax
       if(f_getIntValue(tag_sw_subset_only,iret)==0) sw_subset_only = iret
       if(sw_subset_only==ON)then
         if(printable.and.ipriinputfile>=1) write(nfout,'(a)') ' !** sw_subset_only = ON'
         minxyz(:) = -1; maxxyz(:) = -1
         if(f_getRealValue(tag_min_x,dret,'')==0) minxyz(1) = dret
         if(f_getRealValue(tag_min_y,dret,'')==0) minxyz(2) = dret
         if(f_getRealValue(tag_min_z,dret,'')==0) minxyz(3) = dret
         if (minxyz(1).gt.0 .and. minxyz(1).gt.1) then
            write(nfout,'(a,f10.5)') ' !** WARNING min_x is invalid : ',minxyz(1)
            minxyz(1) = -1
         endif
         if (minxyz(2).gt.0 .and. minxyz(2).gt.1) then
            write(nfout,'(a,f10.5)') ' !** WARNING min_y is invalid : ',minxyz(2)
            minxyz(2) = -1
         endif
         if (minxyz(3).gt.0 .and. minxyz(3).gt.1) then
            write(nfout,'(a,f10.5)') ' !** WARNING min_z is invalid : ',minxyz(3)
            minxyz(3) = -1
         endif

         if(f_getRealValue(tag_max_x,dret,'')==0) maxxyz(1) = dret
         if(f_getRealValue(tag_max_y,dret,'')==0) maxxyz(2) = dret
         if(f_getRealValue(tag_max_z,dret,'')==0) maxxyz(3) = dret
         if (maxxyz(1).gt.0 .and. (maxxyz(1).gt.1 .or. maxxyz(1).lt.minxyz(1))) then
            write(nfout,'(a,f10.5)') ' !** WARNING max_x is invalid : ',maxxyz(1)
            maxxyz(1) = -1
         endif
         if (maxxyz(2).gt.0 .and. (maxxyz(2).gt.1 .or. maxxyz(2).lt.minxyz(2))) then
            write(nfout,'(a,f10.5)') ' !** WARNING max_y is invalid : ',maxxyz(2)
            maxxyz(2) = -1
         endif
         if (maxxyz(3).gt.0 .and. (maxxyz(3).gt.1 .or. maxxyz(3).lt.minxyz(3))) then
            write(nfout,'(a,f10.5)') ' !** WARNING max_z is invalid : ',maxxyz(3)
            maxxyz(3) = -1
         endif

         if(printable.and.ipriinputfile>=1)then
            write(nfout,'(a)') " !** region which shall be taken into account "
            tmpmin=0;tmpmax=1
            if(minxyz(1)>0) tmpmin=minxyz(1)
            if(maxxyz(1)>0) tmpmax=maxxyz(1)
            write(nfout,'(a,f10.5,a,f10.5)') " !** a-axis : ",tmpmin," to ",tmpmax
            tmpmin=0;tmpmax=1
            if(minxyz(2)>0) tmpmin=minxyz(2)
            if(maxxyz(2)>0) tmpmax=maxxyz(2)
            write(nfout,'(a,f10.5,a,f10.5)') " !** b-axis : ",tmpmin," to ",tmpmax
            tmpmin=0;tmpmax=1
            if(minxyz(3)>0) tmpmin=minxyz(3)
            if(maxxyz(3)>0) tmpmax=maxxyz(3)
            write(nfout,'(a,f10.5,a,f10.5)') " !** c-axis : ",tmpmin," to ",tmpmax
         endif
       endif
    end subroutine rd_charge_subset

#endif
    subroutine set_default_btprefix()
      character(len=255) :: dir
      character(len=255) :: nam
      integer :: i,last_slash,last_index
      call getcwd(dir)
      do i=1,255
        if(dir(i:i)==' ') exit
        if(dir(i:i)=='/') last_slash=i
        last_index = i
      enddo
      if(last_index<i) then
        boltztrap_prefix = trim(dir(last_slash+1:i))
      else
        boltztrap_prefix = 'phase0'
      endif
    end subroutine set_default_btprefix

  end subroutine m_CtrlP_rd_postproc

  subroutine m_CntrlP_set_crtdst(a,b,c)
    real(kind=DP) :: a,b,c
    if(sw_aldos==ON) then
       crtdst_aldos = min(a,b,c)/2.0
       write(6,'(" crtdst_aldos = ",f8.4)') crtdst_aldos
    end if
    if(sw_layerdos==ON) then
       if(normal_axis_winlay == 1) then
          crtdst_winlay = a/2;
       else if(normal_axis_winlay == 2) then
          crtdst_winlay = b/2;
       else if(normal_axis_winlay == 3) then
          crtdst_winlay = c/2;
       end if
       write(6,'(" crtdst_winlay = ",f8.4)') crtdst_winlay
    end if
  end subroutine m_CntrlP_set_crtdst

  ! --------------------------------
  subroutine m_CtrlP_rd_parameters(nfinp,nfout,nlines)
    integer, intent(in) :: nfinp,nfout
    integer, intent(out):: nlines

    integer, parameter :: NWK = 6
    integer :: natm, ntyp, i
!!$    integer :: lentag
    logical :: isterm, skip_to_tagbegin, skip_to_tagbegin2
    if(printable) write(nfout,'("<<< m_CtrlP_rd_parameters >>>")')

    allocate(work(NWK))
    rewind nfinp
    nlines = 0
!       ==========================
    call read_gmax_gmaxp_natm_ntyp(natm, ntyp)
    if(printable) write(nfout,'(" !** natm = ", i6, " ntyp = ", i6)') natm, ntyp
    call skip_lines(nfinp,4+natm+ntyp)
    call read_icond_iconstpw
    call read_ipri
    call read_nmd1_nmd2_etc
    call read_mixing_parameters
    call read_charge_precon(isterm)
    if(isterm) nlines = nlines + 1
    call read_dtim_1234(isterm)
    call read_dtio_imdalg_iexpl_edelta
    call read_width_forccr_istress
#ifndef _EMPIRICAL_
    call read_xctype_nspin
    call read_destm
    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    call skiplines_and_read_neg  ! neg
    call skip_lines(nfinp,1)
!!$    call down2nextst
!!$    call skip_lines(nfinp,1)
!!$    call skip_lines(nfinp,5)
    call read_intzaj_imatrix_diagon
    call read_gmaxs_or_n_matrix_size(isterm)
    call read_imsd(isterm)
    call read_evaluation_eko_diff_submat(neg)

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    rewind nfinp
    call read_solver_numbers  ! -(contained here)

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_precalculation)) then
       call read_itagvalue(nfinp,tag_precalculation,tag_nel_Ylm,nel_Ylm)
    end if

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_matdiagon)) then
       if(printable) write(nfout,'(" !* - tag for matrix diagonalization -")')
       call read_dtagvalue(nfinp,tag_matdiagon,tag_eps_solve_Hx_eq_ex,eps_solve_Hx_eq_ex)
    end if

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
!!$    rewind nfinp
!!$    call read_imGSrmm           ! -(contained here) -> imGSrmm
!!$    call read_rr_Critical_Value ! -(contained here) -> rr_Critical_Value
!!$    call read_rmm_printout       ! -(contained here) -> rmm_printout
!!$    call read_rmm_precal_phase_matm     ! -(contained here) -> rmm_precal_phase_matm

    if(skip_to_tagbegin(nfinp,tag_rmm)) then
       if(printable) write(nfout,'(" !* -- tag for rmm --")')
       call read_itagvalue(nfinp,tag_rmm,tag_imGSrmm, imGSrmm)
       call read_dtagvalue(nfinp,tag_rmm,tag_rr_critical_value,rr_Critical_Value)
       call read_itagvalue(nfinp,tag_rmm,tag_rmm_printout, rmm_printout)
       call read_itagvalue(nfinp,tag_rmm,tag_rmm_precal_phase_matm,rmm_precal_phase_matm)
       call read_itagvalue(nfinp,tag_rmm,tag_rmm3_bisec_trial_max, rmm3_bisec_trial_max)
       call read_dtagvalue(nfinp,tag_rmm,tag_rmm3_bisec_crtcl_value, rmm3_bisec_crtcl_value)
       call read_dtagvalue(nfinp,tag_rmm,tag_edelta_change_to_rmm, edelta_change_to_rmm)

    else
       if(printable) write(nfout,'(" !* -- tag for rmm is not found --")')
    end if

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_lineminimization)) then
       if(printable) write(nfout,'(" !* -- tag for liniminimization --")')
       call read_dtagvalue(nfinp,tag_lineminimization,tag_dt_lower_critical,dt_lower_critical)
       call read_dtagvalue(nfinp,tag_lineminimization,tag_dt_upper_critical,dt_upper_critical)
       call read_dtagvalue(nfinp,tag_lineminimization,tag_delta_lmdenom, delta_lmdenom)
    end if
    call read_fine_STM_simulation ! -(c. h.)        -> sw_fine_STM_simulation

#endif

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_gdiis)) then
       if(printable) write(nfout,'(" !* -- tag for GDIIS --")')
       call read_itagvalue(nfinp,tag_gdiis,tag_gdiis_box_size,kqnmditer_p)
       call read_gdiis_hownew(nfinp,tag_gdiis,tag_gdiis_hownew,gdiis_hownew)
       call read_dtagvalue(nfinp,tag_gdiis,tag_c_forc_prop_region_high,c_forc_prop_region_high)
       call read_dtagvalue(nfinp,tag_gdiis,tag_c_forc_prop_region_low,c_forc_prop_region_low)
       call read_dtagvalue(nfinp,tag_gdiis,tag_factor_prop_region, factor_prop_region)
       call read_dtagvalue(nfinp,tag_gdiis,tag_c_forc2gdiis, c_forc2GDIIS)
       call read_itagvalue(nfinp,tag_gdiis,tag_c_iteration2GDIIS, c_iteration2GDIIS)
    end if

    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    rewind nfinp
    call set_printoutlevel_default(ipri)
    if(skip_to_tagbegin(nfinp,tag_printoutlevel)) then
       if(printable) write(nfout,'(" !* -- tag for printoutlevel --")')
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipritiming, ipritiming)
#ifndef _EMPIRICAL_
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprisolver, iprisolver)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprievdff,  iprievdff)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprirmm,    iprirmm); rmm_printout = iprirmm
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripulay,  ipripulay)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprisnl,    iprisnl)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriphig,   ipriphig)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripao,    ipripao)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriberry,  ipriberry)
! === KT_add == 13.1R
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprisym, iprisym)
! ============= 13.1R
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriphonon, ipriphonon)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprifef,    iprifef)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprimatdiagon,iprimatdiagon)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprivlhxcq, iprivlhxcq)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprisubmat, iprisubmat)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipribetar,  ipribetar)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripaw,  ipripaw)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprixc,   iprixc)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprieigenvalue, iprieigenvalue)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipritotalcharge, ipritotalcharge)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprivloc,   iprivloc)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprichargemixing, iprichargemixing)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripositron, ipripositron)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriekzaj,  ipriekzaj)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprispg,    ipri_spg)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprikp,     ipri_kp)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprioccup, iprioccup)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipridos, ipridos)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprinegativecharge,iprinegativecharge)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripp, ipripp)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriwf, ipriwf)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipricoefwf, ipricoefwf)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprichargedensity, iprichargedensity)
#endif
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprigdiis,  iprigdiis)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipristrcfctr,ipristrcfctr)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriparallel, ipriparallel)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprifftmap, iprifftmap)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriinputfile,ipriinputfile)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprimd, iprimd)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriforce, ipriforce)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriexx,   ipriexx )
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriblsize,   ipriblsize )
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprivelocity, iprivelocity)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriparadeb, ipriparadeb)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprijobstatus, iprijobstatus)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipripredictor, ipripredictor)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriunitcell, ipriunitcell)
       call read_itagvalue(nfinp,tag_printoutlevel, tag_iprifire, iprifire)

! === KT_add === 2014/09/24
       call read_itagvalue(nfinp,tag_printoutlevel, tag_ipriepsilon, ipriepsilon )
! ============== 2014/09/24
    end if

    call multiply_printable_factor()

#ifndef _EMPIRICAL_
    ! ---- wf_solvers ----
    if(skip_to_tagbegin(nfinp,tag_solver_of_WF)) then
       tag_solver_of_WF_is_found = .true.
       call read_itagvalue(nfinp,tag_solver_of_WF, tag_n_before, n_WF_solvers_before)
       call read_itagvalue(nfinp,tag_solver_of_WF, tag_n_after,  n_WF_solvers_after)
       n_WF_solvers_all = n_WF_solvers_before + n_WF_solvers_after
       ! -- allocation and initialization of "w_solver"
       call alloc_w_solver(n_WF_solvers_all)
!!$       allocate(w_solver(n_WF_solvers_all))
!!$       do i = 1, n_WF_solvers_all
!!$          w_solver(i)%before_or_after_convergence = BEFORE
!!$          w_solver(i)%solver = MSD
!!$          w_solver(i)%till_n_iter = -1
!!$          w_solver(i)%precon = YES
!!$          w_solver(i)%iter_range = 100
!!$          w_solver(i)%variation_way = varLINEAR
!!$          w_solver(i)%cmix_pointer = 1
!!$          w_solver(i)%dtim_s = 0.1d0
!!$          w_solver(i)%dtim_e = 0.1d0
!!$       end do
       if(ipriinputfile >= 1) then
          write(nfout,*) " !* tag_solver_of_WF_is_found = ", tag_solver_of_WF_is_found
          write(nfout,'(" !** n_WF_solvers_before = ",i5)') n_WF_solvers_before
          write(nfout,'(" !** n_WF_solvers_after  = ",i5)') n_WF_solvers_after
       end if
       if(skip_to_tagbegin2(nfinp,tag_n_before)) then
          do i = 1, n_WF_solvers_before
             call read_solver_detail(nfinp,i)
             w_solver(i)%before_or_after_convergence = BEFORE
          end do
       else
          if(ipriinputfile >= 1) write(nfout,'(" !* tag_n_before is not found")')
       end if
       if(skip_to_tagbegin2(nfinp,tag_n_after)) then
          do i = n_WF_solvers_before+1, n_WF_solvers_all
             call read_solver_detail(nfinp,i)
             w_solver(i)%before_or_after_convergence = AFTER
          end do
       else
          if(ipriinputfile >= 1) write(nfout,'(" !* tag_n_after is not found")')
       end if
    end if

    ! ---- charge_mixing ----
    if(skip_to_tagbegin(nfinp,tag_charge_mixing)) then
       tag_charge_mixing_is_found = .true.
       call read_itagvalue(nfinp,tag_charge_mixing, tag_n_mixing_way, n_Charge_Mixing_way)
       if(n_Charge_Mixing_way >= 1) then
          call alloc_charge_mixing(n_Charge_Mixing_way)
!!$          allocate(charge_mixing(n_Charge_Mixing_way))
!!$          do i = 1, n_Charge_Mixing_way
!!$             charge_mixing(i)%variation_way = varLINEAR
!!$             charge_mixing(i)%precon = YES
!!$             charge_mixing(i)%hownew = ANEW
!!$             charge_mixing(i)%cutoff = LARGE
!!$             charge_mixing(i)%istr   = 1
!!$             charge_mixing(i)%nbxmix = 0
!!$             charge_mixing(i)%rmxs = 0.1d0
!!$             charge_mixing(i)%rmxe = 0.5d0
!!$          end do
       end if
       do i = 1, n_Charge_Mixing_way
          call read_charge_mixing_detail(nfinp,i)
       end do
    end if

    ! -- dos ---
    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_dos)) then
       call read_itagvalue(nfinp,tag_dos, tag_sw_dos_gaussdistrib, sw_dos_gaussdistrib)
       call read_dtagvalue(nfinp,tag_dos, tag_deltaE_dos_gaussD,   deltaE_dos)
       call read_dtagvalue(nfinp,tag_dos, tag_variance_dos_GaussD, variance_dos_GaussD)
       call read_itagvalue(nfinp,tag_dos, tag_nwd_dos_window_width,nwd_dos_window_width)
       call read_energy_unit(nfinp,tag_dos,tag_energy_unit,dos_energy_unit)
       if(dos_energy_unit == EV_ENERGY_UNIT) then
          deltaE_dos = deltaE_dos/Hartree               ! unit: eV -> Hartree
          variance_dos_GaussD = variance_dos_GaussD/(Hartree*Hartree) ! unit: eV**2 -> Hartree**2
       end if
    end if

    ! -- charge --
    if(npes >= 2) call mpi_barrier(MPI_CommGroup,ierr)
    if(skip_to_tagbegin(nfinp,tag_charge)) then
       call read_itagvalue(nfinp,tag_charge,tag_sw_charge_rspace, sw_charge_rspace)
    end if

#endif
!  endif of "_EMPIRICAL_"
!       ==========================
    deallocate(work)
  contains
#ifndef _EMPIRICAL_
    subroutine read_charge_mixing_detail(nf,i)
      integer, intent(in) :: nf,i
      integer :: ip, ipn
      integer :: strncmp3
      read(nf,'(a132)',end=2, err=2) str

      ! -- mixing_way
      ip = strncmp3(str,":")
      ipn = strncmp3(str,"{")
      if(ipn-1 >= ip+1) then
         if(strncmp3(str(ip+1:ipn-1),"simple")>=1) then
            charge_mixing(i)%mixing_way = SIMPLE
         else if(strncmp3(str(ip+1:ipn-1),"broyden2")>=1) then
            charge_mixing(i)%mixing_way = BROYD2
         else if(strncmp3(str(ip+1:ipn-1),"broyden") >=1) then
            charge_mixing(i)%mixing_way = BROYD1
         else if(strncmp3(str(ip+1:ipn-1),"dfp")>=1) then
            charge_mixing(i)%mixing_way = DFP
         else if(strncmp3(str(ip+1:ipn-1),"pulay") >= 1) then
            charge_mixing(i)%mixing_way = PULAY
         end if
      end if

      ! -- tag_itr
      ip = strncmp3(str,tag_itr)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) charge_mixing(i)%iter_range
      end if
      ! -- tag_var
      ip = strncmp3(str,tag_var)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"linear") >= 1) then
            charge_mixing(i)%variation_way = varLINEAR
         else if(strncmp3(str(ip+ipn:len_str),"tanh") >= 1) then
            charge_mixing(i)%variation_way = varTANH
         end if
      end if
      ! -- tag_prec --
      ip = strncmp3(str,tag_prec)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"yes")>= 1) then
            charge_mixing(i)%precon = YES
         else
            charge_mixing(i)%precon = NO
         end if
      end if
      ! -- tag_hownew --
      ip = strncmp3(str,tag_hownew)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"anew") >= 1) then
            charge_mixing(i)%hownew = ANEW
         else if(strncmp3(str(ip+ipn:len_str),"renew") >= 1) then
            charge_mixing(i)%hownew = RENEW
         end if
      end if
      ! -- tag_cutoff --
      ip = strncmp3(str,tag_cutoff)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"large") >= 1) then
            charge_mixing(i)%cutoff = LARGE
         else if(strncmp3(str(ip+ipn:len_str),"MEDIUM") >= 1) then
            charge_mixing(i)%cutoff = MEDIUM
         else if(strncmp3(str(ip+ipn:len_str),"SMALL") >= 1) then
            charge_mixing(i)%cutoff = SMALL
         end if
      end if
      ! -- tag_istr --
      ip = strncmp3(str,tag_istr)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) charge_mixing(i)%istr
      end if
      ! -- tag_nbxmix --
      ip = strncmp3(str,tag_nbxmix)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) charge_mixing(i)%nbxmix
      end if
      ! -- tag_rmxs --
      ip = strncmp3(str,tag_rmxs)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) charge_mixing(i)%rmxs
      end if
      ! -- tag_rmxe --
      ip = strncmp3(str,tag_rmxe)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) charge_mixing(i)%rmxe
      end if

      if(ipriinputfile >= 1) then
         write(nfout,'(" !** way=",i4," rmxs=",f5.2," rmxe=",f5.2, " itr=", i3 &
              & , " var=", i2, " prec=",i2, " istr=", i2," nbxmix=",i2 &
              & , " cutoff=",i2, " hownew= ",i2)') charge_mixing(i)%mixing_way&
              & , charge_mixing(i)%rmxs, charge_mixing(i)%rmxe, charge_mixing(i)%iter_range &
              & , charge_mixing(i)%variation_way, charge_mixing(i)%precon &
              & , charge_mixing(i)%istr, charge_mixing(i)%nbxmix &
              & , charge_mixing(i)%cutoff, charge_mixing(i)%hownew
      end if
2     continue
    end subroutine read_charge_mixing_detail

    subroutine read_solver_detail(nf,i)
      integer, intent(in) :: nf,i
      integer :: ip, ipn,n
      integer :: strncmp3
      read(nf,'(a132)',end=2,err=2) str

      ! -- tag_sol --
      ip = strncmp3(str,tag_sol)
      ipn = strncmp3(str(ip:len_str),"=")
      if(strncmp3(str(ip+ipn:len_str),"lm+MSD") >= 1) then
         w_solver(i)%solver = lmMSD
      else if(strncmp3(str(ip+ipn:len_str),"MSD") >= 1) then
         w_solver(i)%solver = MSD
      else if(strncmp3(str(ip+ipn:len_str),"lm+SD") >= 1) then
         w_solver(i)%solver = lmSD
      else if(strncmp3(str(ip+ipn:len_str),"SD") >= 1) then
         w_solver(i)%solver = SD
      else if(strncmp3(str(ip+ipn:len_str),"RMM2P") >= 1) then
         w_solver(i)%solver = RMM2P
      else if(strncmp3(str(ip+ipn:len_str),"RMM3") >= 1) then
         w_solver(i)%solver = RMM3
      else if(strncmp3(str(ip+ipn:len_str),"RMM2") >= 1) then
         w_solver(i)%solver = RMM2
      else if(strncmp3(str(ip+ipn:len_str),"MatrixDiagon") >= 1) then
         w_solver(i)%solver = MATRIXDIAGON
      else if(strncmp3(str(ip+ipn:len_str),"SUBMAT") >= 1) then
         w_solver(i)%solver = SUBMAT
      else if(strncmp3(str(ip+ipn:len_str),"CG") >= 1) then
         w_solver(i)%solver = CG
      else
         if(ipriinputfile >= 1) write(nfout,'(" w_solver%solver is not defined")')
      end if

      ! -- tag_prec --
      ip = strncmp3(str,tag_prec)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"yes")>= 1) then
            w_solver(i)%precon = YES
         else
            w_solver(i)%precon = NO
         end if
      end if

      ! -- tag_till_n
      ip = strncmp3(str,tag_till_n)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) n
         w_solver(i)%till_n_iter = n
      end if

      ! -- tag_dts
      ip = strncmp3(str,tag_dts)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) w_solver(i)%dtim_s
      end if
      ! -- tag_dte
      ip = strncmp3(str,tag_dte)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) w_solver(i)%dtim_e
      end if
      ! -- tag_itr
      ip = strncmp3(str,tag_itr)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) w_solver(i)%iter_range
      end if
      ! -- tag_var
      ip = strncmp3(str,tag_var)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         if(strncmp3(str(ip+ipn:len_str),"linear") >= 1) then
            w_solver(i)%variation_way = varLINEAR
         else if(strncmp3(str(ip+ipn:len_str),"tanh") >= 1) then
            w_solver(i)%variation_way = varTANH
         end if
      end if

      ! -- cmix
      ip = strncmp3(str,tag_cmix)
      if(ip >= 1) then
         ipn = strncmp3(str(ip:len_str),"=")
         read(str(ip+ipn:len_str),*) w_solver(i)%cmix_pointer
      end if
      if(ipriinputfile >= 1) then
         write(nfout,'(" !** sol=",i3," till_n=",i4," dts = ",f5.2," dte = ",f5.2 &
              & , " itr=",i4," var=",i4," cmix=",i4, " subspaceroation = ",i4)') w_solver(i)%solver &
              & , w_solver(i)%till_n_iter, w_solver(i)%dtim_s,w_solver(i)%dtim_e &
              & , w_solver(i)%iter_range,  w_solver(i)%variation_way &
              & , w_solver(i)%cmix_pointer, w_solver(i)%subspace_rotation
      end if
2     continue
    end subroutine read_solver_detail
#endif

    subroutine down2nextst
      logical :: tf

1     read(nfinp,'(a132)',end=2,err=2) str
      call strncmp2(str,len_str,'nextst',6,tf)
      if(.not.tf) goto 1
      if(tf) goto 3
2     print *,' reached end of file '
3     continue
    end subroutine down2nextst

    subroutine skiplines_and_read_neg
      logical :: tf

1     read(nfinp,'(a132)',end=2,err=2) str
      call strncmp2(str,len_str,'neg',3,tf)
      if(.not.tf) goto 1
      if(tf) goto 3
2     print *,' reached end of file '
3     continue
      read(str,*) neg
    end subroutine skiplines_and_read_neg

    subroutine read_gmax_gmaxp_natm_ntyp(natm,ntyp)
      integer, intent(out) ::  natm, ntyp
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,4,work,nsize)
      gmax  = work(1)
      gmaxp = work(2)
      ntyp = nint(work(3))
      natm = nint(work(4))
    end subroutine read_gmax_gmaxp_natm_ntyp

    subroutine read_icond_iconstpw
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,2,work,nsize)
      icond = nint(work(1))
      if(nsize < 2) then
         iconstpw = 0
      else
         iconstpw = nint(work(2))
      endif
      if(printable) write(nfout,600) icond, iconstpw
600   format(' !** ',2i4,'   : icond 0-md, 1-cont.md, ', &
           &               '2-wave fn,, 3-charge den.,  iconstpw')
    end subroutine read_icond_iconstpw

    subroutine read_ipri
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,2,work,nsize)
      ipri = nint(work(2))
      if(printable) write(nfout,610) nint(work(1)), ipri
  610 format(' !** ',2i4,'   : ipre, ipri')
!!$!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
!!$      ipri_spg = ipri + 1
!!$      ipri_kp = ipri + 1
!!$!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
    end subroutine read_ipri

    subroutine read_nmd1_nmd2_etc
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,5,work,nsize)
      max_scf_iteration = nint(work(1))
      max_total_scf_iteration = nint(work(2))
      iter_last = nint(work(3))
      cpumax = work(4)
      if(nsize < 5) then
         ifstop = 0
      else
         ifstop = nint(work(5))
         if(ifstop < 0) ifstop = 0
      endif

      if(printable) write(nfout,621) max_scf_iteration,max_total_scf_iteration,iter_last,cpumax, ifstop
621   format(' !** ',3i6,f12.2,i3 &
           &     ,' : nmd1, nmd2, iter_last, cpumax, ifstop')
      if(icond == INITIAL .and. iter_last >= 1) then
         iter_last = 0
         if(printable) write(nfout,'(" !*  last_last is set 0")')
      endif
      if(icond == CONTINUATION .and. (iter_last > max_total_scf_iteration)) then
         if(printable) write(nfout,'(" !*  a combination of max_total_scf_iteration(nmd2) and iter_last  is invalid.")')
         call phase_error_with_msg(nfout,' a combination of nmd2 and iter_last is invalid <<read_nmd1_nmd2_etc>>'&
                                  ,__LINE__,__FILE__)
      endif
    end subroutine read_nmd1_nmd2_etc

    subroutine read_mixing_parameters
      character(len=6) mixnam(PULAY)
      data mixnam/'SIMPLE','BROYD1','BROYD2','DFP','PULAY'/
      character(len=5) newnam(RENEW)
      data newnam/'ANEW ','RENEW'/
      character(len=6) cutoffsize(3)
      data cutoffsize/'SMALL ','MEDIUM','LARGE '/

      integer  :: nsize, i

      read(nfinp,'(a132)') str
      if(printable) write(nfout,'(a132)') str
      call ch4mix(str,len_str,waymix,istrbr,nbxmix,hownew&
           &,cutoff_mix,decomp_into_charge_spin,work,nsize)
      rmx_1234(1:4) = work(1:4)
      if(icond == CONTINUATION .and. nsize >= 5) then
         rmx = work(5)
      else if(icond == CONTINUATION .and. nsize < 5) then
         if(printable) write(nfout,'(" !*  No data for rmx. rmx4 is substituted into rmx")')
         rmx = rmx_1234(4)
      endif
      if(waymix /= SIMPLE) then
         if(printable) write(nfout,60) mixnam(waymix),istrbr,nbxmix,newnam(hownew),work(1),cutoffsize(cutoff_mix)
60       format(' !** ',a6,' istrbr = ',i3,', nbxmix = ',i3,', hownew = ',a5 &
              &   ,', alpha = ',f6.2,' cutoff = ',a6)
         waymix_b = waymix
         istrbr_b = istrbr
         nbxmix_b = nbxmix
         hownew_b = hownew
         cutoff_mix_b = cutoff_mix
      else if(waymix == SIMPLE) then
         if(printable) write(nfout,61) mixnam(waymix), (work(i),i=1,nsize)
61       format(' !** ',a6,' rmx1 = ',f6.2,' rmx2 = ',f6.2,' rmx3 = ',f6.2 &
              &   ,' rmx4 = ',f6.2,' rmx = ',f6.2)
      endif
    end subroutine read_mixing_parameters

    subroutine read_charge_precon(isterm)
      logical, intent(out) :: isterm

      character(len=8):: cprnam(2)
      data cprnam/'cprecond','nocpreco'/
      logical ::         precon_charge_only

      read(nfinp,'(a132)') str
      call ch4prc(str,len_str,c_precon,amix,bmix,precon_charge_only,isterm)
      if(printable) then
         if(c_precon) then
            if(precon_charge_only) then
               write(nfout,655) cprnam(1), amix, bmix
            else
               write(nfout,651) cprnam(1), amix, bmix
            endif
         else
            write(nfout,652) cprnam(2)
         endif
      end if
 651  format(' !* ',a8,' amix = ', f6.2,' bmix = ', f6.2 )
 655  format(' !* ',a8,' amix = ', f6.2,' bmix = ', f6.2, ' precon_charge_only')
 652  format(' !* ',a8)
   end subroutine read_charge_precon

    subroutine read_dtim_1234(isterm)
      logical, intent(in) :: isterm
      integer  :: nsize,mm
      if(isterm) read(nfinp,'(a132)') str
      if(printable) write(nfout,'(a132)') str
      call chnnm(str,len_str,5,work,nsize)
      if(nsize < 4)  call phase_error_with_msg(nfout,' shortage of data for dtims',__LINE__,__FILE__)
      dtim_1234(1:4) = work(1:4)
      if(icond == CONTINUATION .and. nsize >= 5) then
         dtim_initial = work(5)
      else if(icond == CONTINUATION .and. nsize < 5) then
         if(printable) write(nfout,'(" !* shortage of data for dtim_1234(2)")')
         dtim_initial = dtim_1234(4)
      endif
      if(icond == CONTINUATION) then
         if(printable) write(nfout,641) (dtim_1234(mm), mm = 1, 4), dtim_initial
      else
         if(printable) write(nfout,640) (dtim_1234(mm), mm = 1, 4)
      endif
640   format(' !** ',4f6.2,'          : dtim1, dtim2, dtim3, dtim4')
641   format(' !** ',5f6.2,'          : dtim1,dtim2,dtim3,dtim4,dtim')
    end subroutine read_dtim_1234

    subroutine read_dtio_imdalg_iexpl_edelta
      read(nfinp,*) dtio,imdalg,iexpl,edelta
      if(printable) write(nfout,650) dtio,imdalg,iexpl,edelta
650   format(' !** ',f6.2,2i6,d12.2,'    : dtio ,imdalg, iexpl, edelta')
    end subroutine read_dtio_imdalg_iexpl_edelta

    subroutine read_width_forccr_istress
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,3,work,nsize)
      width  = work(1)
      forccr = work(2)
      if(nsize < 3) then
         istress = 0
      else
         istress = nint(work(3))
      endif
      if(printable) write(nfout,660) width,forccr,istress
660   format(' !** ',f8.4,d10.2,i5,'           : width,forccr,istress')
      if(width >= 0.d0) then
         way_of_smearing = PARABOLIC
      else if(width > -10.0) then
         way_of_smearing = MP
      else
         way_of_smearing = TETRAHEDRON
      end if
    end subroutine read_width_forccr_istress

#ifndef _EMPIRICAL_
    subroutine read_xctype_nspin
      read(nfinp,'(a132)') str
      call chchn0(str,len_str,len_xctype,xctype,nspin)
      if(nspin <= 0 .or.nspin >= 3) then
         if(nspin <= 0) then
            nspin = 1
         else if(nspin >= 3) then
            nspin = 2
         endif
         if(printable) write(nfout,'(" nspin is reset to be ",i6)') nspin
      endif
      if(printable) write(nfout,665) xctype, nspin
665   format(' !** ',a7,i10,'                   : xctype,nspin,zeta1')
    end subroutine read_xctype_nspin

    subroutine read_destm
      integer nsize
      read(nfinp,'(a132)') str
      call chnnm(str, len_str, 2, work, nsize)
      if(nsize < 1) then
         call phase_error_with_msg(nfout,' ! no data for destm',__LINE__,__FILE__)
      else
         destm = work(1)
         if(nsize == 1) then
            n_stm = M_STM
         else
            n_stm = work(2) + 0.5
         endif
      endif
      if(nsize == 1) then
         if(printable) write(nfout,666) destm
      else
         if(printable) write(nfout,667) destm, n_stm
      endif
 666  format(' !** ',f6.2,'                             : destm')
 667  format(' !** ',f6.2,' ', i5,'                     : destm, n_stm')
    end subroutine read_destm

    subroutine read_intzaj_imatrix_diagon
      integer :: nsize
      read(nfinp,'(a132)') str
      call chnnm(str,len_str,2,work,nsize)
      if(nsize == 0) then
         intzaj = by_matrix_diagon
      else
         intzaj = nint(work(1))
         if(nsize == 1 .or. intzaj == by_random_numbers) then
            imatrix_diagon = 1
         else
            imatrix_diagon = nint(work(2))
         endif
      endif

      if(intzaj == by_random_numbers) then
         if(printable) write(nfout,730) intzaj
      else if(intzaj == by_matrix_diagon) then
         if(printable) write(nfout,731)intzaj, imatrix_diagon
      else
         if(printable) write(nfout,'(" !** intzaj = ",i5)') intzaj
      endif
730   format(' !** ',i8,'       : 0 = random numbers, 1 = matrix diagon')
731   format(' !** ',i8,i6,' : 0 = random numbers, 1 = matrix diagon(#imd)')
    end subroutine read_intzaj_imatrix_diagon

    subroutine read_gmaxs_or_n_matrix_size(isterm)
      logical, intent(out) :: isterm
      integer :: n
      isterm = .false.
      read(nfinp,'(a132)') str
      n = scan(str,'gG')
      if(n > 0 .and. n < len_str) then
         if(str(n:n+4) == 'gmaxs' .or. str(n:n+4) == 'GMAXS' ) then
            isterm = .true.
            n = scan(str,'=')
            read(str(n+1:len_str),*) gmaxs_given
         endif
      endif
      if(isterm) then
         if(printable) write(nfout,653) gmaxs_given
         return
      endif

      n = scan(str,'nN')
      if(n > 0 .and. n < len_str) then
         if(str(n:n+7) == 'n_matrix' .or. str(n:n+4) == 'N_MATRIX' ) then
            isterm = .true.
            n = scan(str,'=')
            read(str(n+1:len_str),*) n_matrix_size
         endif
      endif
      if(isterm .and. printable) write(nfout,654) n_matrix_size

 653  format(' !**  gmaxs         = ', f8.4)
 654  format(' !**  n_matrix_size = ', i8)
    end subroutine read_gmaxs_or_n_matrix_size

    subroutine read_imsd(isterm)

      logical, intent(in) :: isterm
      integer nsize
      if(isterm) read(nfinp,'(a132)') str
      call chnnm(str,len_str,4,work,nsize)
      imsd = nint(work(1))
      if(printable) write(nfout,'(" !** - imsd = ",i6, " nsize = ", i6)') imsd, nsize
      if(imsd == CG     .or. imsd == CGPRC .or. imsd == RMM .or. &
           &  imsd == RMMPRC .or. imsd == RMM2  .or. imsd == RMM2PRC .or. &
           &  imsd == RMM2P  .or. imsd == RMM2PPRC .or. &
           &  imsd == SDLM   .or. imsd == SDLMPRC .or.&
           &  imsd == eazyCG .or. imsd == eazyCGPRC) then
         if((imsd == RMM .or. imsd == RMMPRC .or. imsd == RMM2 .or. &
              &   imsd == RMM2PRC .or. imsd == RMM .or. &
              &   imsd == RMM2P .or. imsd == RMM2PPRC) .and. nsize >=4) then
            i_2lm            = nint(work(2))
            i_sd2another     = nint(work(3))
            iwrksz_rmm_phase = nint(work(4))
         else if((imsd == SDLM .or. imsd == SDLMPRC) .and. nsize >=2) then
            i_2lm            = nint(work(2))
            i_sd2another     = max_total_scf_iteration + 30
            iwrksz_rmm_phase = 0
         else if((imsd == SDLM .or. imsd == SDLMPRC) .and. nsize ==1) then
            i_2lm            = default_sd2cg
            i_sd2another     = max_total_scf_iteration + 30
            iwrksz_rmm_phase = 0
         else if(nsize >= 3) then
            i_2lm            = nint(work(2))
            i_sd2another     = nint(work(3))
            iwrksz_rmm_phase = 0

         else if(nsize >= 2) then
            i_2lm            = nint(work(2))
            i_sd2another     = default_sd2cg
            iwrksz_rmm_phase = 0
         endif
         iwrksz_rmm_phase = iwrksz_rmm_phase*1024*1024/irsize
      endif
      if(i_sd2another < i_2lm ) i_sd2another = i_2lm
      if((imsd == RMM2P .or. imsd == RMM2PPRC) .and. i_sd2another == i_2lm) &
           & i_sd2another = i_2lm + 1
      if(printable) then
         if(imsd >= RMM) then
            write(nfout,742) imsd, i_2lm, i_sd2another, iwrksz_rmm_phase
         else if(imsd >= CG) then
            write(nfout,741) imsd, i_2lm, i_sd2another
         else
            write(nfout,740) imsd
         endif
      end if
740   format(' !** ',i8,'       : imsd' )
741   format(' !** ',i8,2i7,' : imsd, i_2lm, i_sd2another' )
742   format(' !** ',i8,3i7&
           &     ,'(MB) : imsd, i_2lm, i_sd2another, wksz for phase' )

      if(imsd == CG .or. imsd == CGPRC .or. imsd == RMM .or. &
           & imsd == RMMPRC .or. imsd == RMM2 .or. imsd == RMM2PRC .or. &
           & imsd == RMM2P  .or. imsd == RMM2PPRC .or. &
           & imsd == eazyCG .or. imsd == eazyCGPRC) then
         n_WF_solvers = 3
         call alloc_p_WF_solvers  ! -(m_Control_Parameters)
         WF_solver(1) = MSD   + On_or_Off_precon_for_WFs_imsd(imsd)
         WF_solver(2) = MSDLM + On_or_Off_precon_for_WFs_imsd(imsd)
         WF_solver(3) = imsd
         till_n_iteration(1) = i_2lm        - 1
         till_n_iteration(2) = i_sd2another - 1
         till_n_iteration(3) = max_total_scf_iteration         + 1
      else if(imsd == SDLM .or. imsd == SDLMPRC) then
         n_WF_solvers = 2
         call alloc_p_WF_solvers ! -(m_Control_Parameters)
         WF_solver(1) = SD    - On_or_Off_precon_for_WFs_imsd(imsd)
         WF_solver(2) = imsd
         till_n_iteration(1) = i_2lm        - 1
         till_n_iteration(2) = max_total_scf_iteration         + 1
      else if(imsd == SD .or. imsd == SDPRC .or. imsd == MSD .or. imsd == MSDPRC) then
         n_WF_solvers = 1
         call alloc_p_WF_solvers ! -(m_Control_Parameters)
         WF_solver(1) = imsd
         till_n_iteration(1) = max_total_scf_iteration         + 1
      end if
    end subroutine read_imsd

    subroutine read_evaluation_eko_diff_submat(neg)
      integer, intent(in) :: neg
      read(nfinp,*, end = 1004, err = 1004) evaluation_eko_diff
      read(nfinp,*, end = 1005, err = 1005) ldiag, meg, damp
      if(meg > neg) meg = neg
      goto 1002
1004  evaluation_eko_diff = 0
1005  ldiag  = max_total_scf_iteration + 1
      meg    = neg
      damp   = 1.d0
1002  continue
      if(printable) then
         write(nfout,750) evaluation_eko_diff
         write(nfout,760) ldiag,meg,damp
      end if
750   format(' !** ',i8,'       : evaluation of eko difference.' &
           &     ,'0 = no ,1 = yes')
760   format(' !** ',2i6,f12.6,': ldiag , meg ,damp ')
    end subroutine read_evaluation_eko_diff_submat

    subroutine read_solver_numbers
      logical :: tf
      integer :: nsize, i

1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'solver of wfs',13,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* - tag of solver is found -")')
         call chnnm2(str,len_str,1,work,nsize)
         n_WF_solvers = nint(work(1))
         if(printable) write(nfout,'(" !** n_WF_solvers = ",i5)') n_WF_solvers
         call alloc_p_WF_solvers ! -(m_Control_Parameters)
         if(n_WF_solvers > NWK) then
            deallocate(work)
            allocate(work(n_WF_solvers))
         end if
         read(nfinp,'(a132)') str
         call chnnm2(str,len_str,n_WF_solvers,work,nsize)
         call check_of_n_WF_solvers(nsize) !-(m_Control_Parameters)
         do i = 1, n_WF_solvers
            WF_solver(i) = nint(work(i))
         end do

         read(nfinp,'(a132)') str
         call chnnm2(str,len_str,n_WF_solvers,work,nsize)
         call check_of_n_WF_solvers(nsize) !-(m_Control_Parameters)
         do i = 1, n_WF_solvers
            till_n_iteration(i) = nint(work(i))
         end do
      else
         if(printable) write(nfout,'(" !* no tag of solver")')
      end if

      if(printable) then
         write(nfout,'(" !** n_WF_solvers = ",i6)') n_WF_solvers
         write(nfout,'(" !** solvers    = ",8i6)') (WF_solver(i),i=1,n_WF_solvers)
         write(nfout,'(" !** #iter last = ",8i6)') (till_n_iteration(i),i=1,n_WF_solvers)
      end if
    end subroutine read_solver_numbers

    subroutine read_fine_STM_simulation
      logical :: tf

1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'sw_fine_stm_simulation',22,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* -- tag of sw_fine_stm_simulation is found --")')
         call strncmp2(str,len_str,'ON',2,tf)
         if(tf) sw_fine_stm_simulation = ON
         if(printable) write(nfout,'(" !** sw_fine_STM_simulation = ",i5)') sw_fine_STM_simulation
      end if
    end subroutine read_fine_STM_simulation

    subroutine check_of_n_WF_solvers(nsize)
      integer, intent(in) :: nsize

      if(nsize < n_WF_solvers) then
         if(printable) write(nfout,'(" !** nsize = ",i6," <  n_WF_solvers = ",i6)') &
              & nsize,n_WF_solvers
         call phase_error_with_msg(nfout,'nsize < n_WF_solvers',__LINE__,__FILE__)
      end if
    end subroutine check_of_n_WF_solvers

    subroutine read_imGSrmm
      logical :: tf
      integer :: nsize

      rewind nfinp
1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'imgsrmm',7,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* - tag of imGSrmm is found -")')
         call chnnm2(str,len_str,1,work,nsize)
         imGSrmm = nint(work(1))
         if(printable) write(nfout,'(" !** imGSrmm = ",i5)') imGSrmm
      else
         if(printable) write(nfout,'(" !* -- no tag of imGSrmm -- ")')
      end if
    end subroutine read_imGSrmm

    subroutine read_rr_Critical_Value
      logical :: tf
      integer :: nsize

      rewind nfinp
1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'rr_critical_value',17,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* -- tag of rr_Critical_Value is found --")')
         call chnnm2(str,len_str,1,work,nsize)
         rr_Critical_Value = work(1)
         if(printable) write(nfout,'(" !** rr_Critical_Value = ",d20.8)') rr_Critical_Value
      else
         if(printable) write(nfout,'(" !* -- tag of rr_Critical_Value is not found --")')
      end if
    end subroutine read_rr_Critical_Value

    subroutine read_rmm_printout
      logical :: tf
      integer :: nsize

      rewind nfinp
1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'rmm_printout',12,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* -- tag of rmm_printout is found --")')
         call chnnm2(str,len_str,1,work,nsize)
         rmm_printout = nint(work(1))
         if(printable) write(nfout,'(" !** rmm_printout = ",i5)') rmm_printout
      end if
    end subroutine read_rmm_printout

    subroutine read_rmm_precal_phase_matm
      logical :: tf
      integer :: nsize

1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,len_str,'rmm_precal_phase_matm',21,tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         if(printable) write(nfout,'(" !* - tag of rmm_precal_phase_matm is found -")')
         call chnnm2(str,len_str,1,work,nsize)
         rmm_precal_phase_matm = nint(work(1))
         if(printable) write(nfout,'(" !**  rmm_precal_phase_matm = ",i5)') rmm_precal_phase_matm
      end if
    end subroutine read_rmm_precal_phase_matm

#endif
! endif of "_EMPIRICAL_"
    subroutine read_gdiis_hownew(nfinp,tag_begin,tag,gdiis_hownew)
      integer, intent(in) ::           nfinp
      character(len=*), intent(in) ::  tag_begin,tag
      integer, intent(out) ::          gdiis_hownew

      logical :: tf, tf2, skip_to_tagbegin
      integer :: ip, strncmp3

      if(.not.skip_to_tagbegin(nfinp,tag_begin)) return

      tf = .false.
1     read(nfinp,'(a132)',end = 1001) str
      call strncmp2(str,3,"end",3,tf2)
      if(tf2) goto 1001
      call strncmp2(str,len_str,tag,len(tag),tf)
      if(.not.tf) goto 1
1001  continue

      if(tf) then
         ip = strncmp3(str,"=")
         if(strncmp3(str(ip+1:len_str),"anew") >= 1 ) then
            gdiis_hownew = ANEW
            if(printable) write(nfout,'(" !* gdiis_hownew is set ANEW")')
         else if(strncmp3(str(ip+1:len_str),"renew") >=1 ) then
            gdiis_hownew = RENEW
            if(printable) write(nfout,'(" !* gdiis_hownew is set RENEW")')
         end if
      else
         if(printable) write(nfout,'(" !* -- no tag value -- ",a31,"  default value is set")') tag
      end if

    end subroutine read_gdiis_hownew

  end subroutine m_CtrlP_rd_parameters
  ! -----------------------------------

  subroutine m_CtrlP_set_ekmode_ON()
    ekmode = ON
  end subroutine m_CtrlP_set_ekmode_ON

  subroutine m_CtrlP_set_ekmode_GRID()
    ekmode = GRID
  end subroutine m_CtrlP_set_ekmode_GRID

  subroutine m_CtrlP_set_uvsormode_ON()
    uvsormode = ON
  end subroutine m_CtrlP_set_uvsormode_ON

  subroutine m_CtrlP_ntcnvg_reset()
    ntcnvg = 0
    if(sub_delta_factor_is_given) sub_ntcnvg = 0
  end subroutine m_CtrlP_ntcnvg_reset

  integer function m_CtrlP_ntcnvg_incre()
    ntcnvg = ntcnvg + 1
    m_CtrlP_ntcnvg_incre = ntcnvg
  end function m_CtrlP_ntcnvg_incre

! --> T. Yamasaki, 25 July 2008
  integer function m_CtrlP_sub_ntcnvg_incre()
    sub_ntcnvg = sub_ntcnvg + 1
    m_CtrlP_sub_ntcnvg_incre = sub_ntcnvg
  end function m_CtrlP_sub_ntcnvg_incre
! <--

  logical function m_CtrlP_ntcnvg_clear()
    m_CtrlP_ntcnvg_clear = .false.
!!$    if(ntcnvg >= MTIMES_CONVERGENCE) m_CtrlP_ntcnvg_clear = .true.
    if(ekmode == ON .or. &
         & (icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION)) then
       if(ntcnvg >= mtimes_convergence_ek) m_CtrlP_ntcnvg_clear = .true.
    else
       if(ntcnvg >= mtimes_convergence_scf) m_CtrlP_ntcnvg_clear = .true.
! --> T. Yamasaki, 25 July 2008
       if(sub_delta_factor_is_given) then
          if(sub_ntcnvg >= sub_mtimes_convergence) m_CtrlP_ntcnvg_clear = .true.
       end if
! <--
    end if
  end function m_CtrlP_ntcnvg_clear

  ! ------- Positron start
  subroutine m_CtrlP_pstrn_ntcnvg_reset()
    positron_ntcnvg = 0
  end subroutine m_CtrlP_pstrn_ntcnvg_reset

  integer function m_CtrlP_pstrn_ntcnvg_incre()
    positron_ntcnvg = positron_ntcnvg + 1
    m_CtrlP_pstrn_ntcnvg_incre = positron_ntcnvg
  end function m_CtrlP_pstrn_ntcnvg_incre

  logical function m_CtrlP_pstrn_ntcnvg_clear()
    m_CtrlP_pstrn_ntcnvg_clear = .false.
    if(positron_ntcnvg >= mtimes_convergence_pev) m_CtrlP_pstrn_ntcnvg_clear = .true.
  end function m_CtrlP_pstrn_ntcnvg_clear

  integer function m_CtrlP_solver_for_pWFs_now(iteration_positron_wf)
    integer, intent(in)       :: iteration_positron_wf

    m_CtrlP_solver_for_pWFs_now = isolver_p
  end function m_CtrlP_solver_for_pWFs_now

  integer function m_CtrlP_submat_for_pWFs_now(iteration_positron_wf)
    integer, intent(in)        :: iteration_positron_wf
    m_CtrlP_submat_for_pWFs_now = sw_submat_p
  end function m_CtrlP_submat_for_pWFs_now

  function m_CtrlP_dtim_p_now(iteration_positron_wf)
    real(kind=DP) :: m_CtrlP_dtim_p_now
    integer, intent(in) :: iteration_positron_wf
    m_CtrlP_dtim_p_now = dtim_p
  end function m_CtrlP_dtim_p_now

  ! ------- Positron end

  logical function m_CtrlP_ntcnvg_pre_clear()
    m_CtrlP_ntcnvg_pre_clear = .false.
!!$    if(ntcnvg >= MTIMES_CONVERGENCE) m_CtrlP_ntcnvg_clear = .true.
    if(ekmode == OFF) then
       if(mtimes_convergence_scf >= 2 .and. ntcnvg >= mtimes_convergence_scf-1) &
            & m_CtrlP_ntcnvg_pre_clear = .true.
    else
       if(mtimes_convergence_ek >= 2 .and. ntcnvg >= mtimes_convergence_ek-1) &
            & m_CtrlP_ntcnvg_pre_clear = .true.
    end if
  end function m_CtrlP_ntcnvg_pre_clear

!!$  interface read_tagvalue
!!$     subroutine read_dtagvalue(nf,tag,value)
!!$       integer, intent(in)          :: nf
!!$       character(len=*), intent(in) :: tag
!!$       real(kind=DP), intent(out)   :: value
!!$     end subroutine read_dtagvalue
!!$     subroutine read_itagvalue(nf,tag,value)
!!$       integer, intent(in)          :: nf
!!$       character(len=*), intent(in) :: tag
!!$       integer, intent(out)         :: value
!!$     end subroutine read_itagvalue
!!$  end interface

  subroutine m_CtrlP_set_rmx(x)
    real(kind=DP) :: x
    rmx_1234 = x
  end subroutine m_CtrlP_set_rmx

  subroutine m_CtrlP_set_kimg(inversion_symmetry, nfout)
    integer, intent(in) :: inversion_symmetry, nfout
    if(inversion_symmetry == ON)  kimg = 1
    if(inversion_symmetry == OFF) kimg = 2
! ============================================ modified by K. Tagami ======= 11.0
!    if(printable) write(6,'(" kimg = ",i5, " (m_CtrlP_set_kimg)")') kimg
!
    if ( noncol ) then
       if ( kimg == 1 ) then
          call phase_error_with_msg(nfout,'inversion_symmetry should be OFF in non-collinear systems.', &
          __LINE__,__FILE__)
       endif
       if(printable) write(nfout,'(" kimg = ",i5, " (m_CtrlP_set_kimg)")') kimg
    else
       if(printable) write(nfout,'(" kimg = ",i5, " (m_CtrlP_set_kimg)")') kimg
    endif
! ========================================================================== 11.0

  end subroutine m_CtrlP_set_kimg

#ifndef _EMPIRICAL_
  subroutine alloc_p_WF_solvers
    if(allocated(WF_solver)) deallocate(WF_solver)
    allocate(WF_solver(n_WF_solvers))
    if(allocated(till_n_iteration)) deallocate(till_n_iteration)
    allocate(till_n_iteration(n_WF_solvers))
  end subroutine alloc_p_WF_solvers

  subroutine m_CtrlP_wd_isolver(nfcntn)
    integer, intent(in)       :: nfcntn
    if(printable) then
       write(nfcntn,*) tag_isolver
       write(nfcntn,'(i10)') previous_solver
    end if
  end subroutine m_CtrlP_wd_isolver

  subroutine m_CtrlP_rd_isolver(nfcntn)
    integer, intent(in)       :: nfcntn
    logical    :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len_tag_isolver,tag_isolver &
            &, EOF_reach, tag_is_found,str,len_str)
       if(.not.tag_is_found) call phase_error_with_msg(6,' tag_isolver is not found',__LINE__,__FILE__)
       read(nfcntn,*) previous_solver
    end if
    if(npes > 1) call mpi_bcast(previous_solver,1,mpi_integer,0 &
         & ,MPI_CommGroup,ierr)
    if(printable) write(6,'(i5, " : previous_solver")') previous_solver
  end subroutine m_CtrlP_rd_isolver

#endif

  integer function m_CtrlP_check_inputfilestyle(nfinp)
    integer, intent(in) :: nfinp
    logical :: skip_to_tagbegin
    m_CtrlP_check_inputfilestyle = NEW_
#ifdef OLD_INPUT_SUPPORT
    m_CtrlP_check_inputfilestyle = OLD
    if(skip_to_tagbegin(nfinp,tag_structure)) then
       m_CtrlP_check_inputfilestyle = NEW_
    end if
#endif
!!$    if(skip_to_tagbegin(nfinp,tag_control)) then
!!$       write(6,'(" !*-- tag_control = ",a32)') tag_control
!!$       if(skip_to_tagbegin(nfinp,tag_accuracy)) then
!!$          write(6,'(" !*-- tag_accuracy = ",a32)') tag_accuracy
!!$          m_CtrlP_check_inputfilestyle = NEW_
!!$       end if
!!$    end if
  end function m_CtrlP_check_inputfilestyle

  subroutine m_CtrlP_wd_iconvergence(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,'(a11)') tag_convergence
       write(nfcntn,'(i10)') iconvergence
    end if
  end subroutine m_CtrlP_wd_iconvergence

  subroutine m_CtrlP_rd_iconvergence(nfcntn)
    integer, intent(in) :: nfcntn
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_convergence),tag_convergence &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) call phase_error_with_msg(6,' tag_convergence is not found',__LINE__,__FILE__)
       read(nfcntn,*) iconvergence_previous_job
    end if
    if(npes > 1) call mpi_bcast(iconvergence_previous_job,1,mpi_integer,0 &
         & ,MPI_CommGroup,ierr)
    if(printable) write(6,'(i5, " : iconvergence_previous_job")') iconvergence_previous_job
    if(icond == CONTINUATION .and. neg_previous < neg) then
       iconvergence_previous_job = 0
       if(printable) write(6,'(" iconvergence_previous_job is reset " &
            & ,i2,", because neg_previous < neg")') iconvergence_previous_job
    end if
    !if(sw_optimize_lattice==ON)then
    !   iconvergence_previous_job = 0
    !endif
  end subroutine m_CtrlP_rd_iconvergence

  subroutine m_CtrlP_reset_iconvergence
     iconvergence = 0
     iconvergence_previous_job = 0
  end subroutine m_CtrlP_reset_iconvergence

  subroutine m_CtrlP_wd_numk_zajsaved(nfcntn,nk)
    integer, intent(in) :: nfcntn
    integer, intent(in) :: nk
    if(mype==0) then
       write(nfcntn,'(a13)') tag_numk_zajsaved
       write(nfcntn,'(i10)') nk
    end if
  end subroutine m_CtrlP_wd_numk_zajsaved

  subroutine m_CtrlP_rd_numk_zajsaved(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_numk_zajsaved),tag_numk_zajsaved &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) call phase_error_with_msg(nfout,' tag_numk_zajsaved is not found',__LINE__,__FILE__)
       read(nfcntn,*) numk_zajsaved
    end if
    if(npes > 1) call mpi_bcast(numk_zajsaved,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(printable) write(nfout,'(i5, " : numk_zajsaved")') numk_zajsaved
  end subroutine m_CtrlP_rd_numk_zajsaved

  subroutine m_CtrlP_wd_iconv_ek(kv3_ek,iconv_ek,nfcntn)
    integer, intent(in) :: kv3_ek,nfcntn
    integer, intent(in), dimension(kv3_ek) :: iconv_ek
    integer :: i
    if(mype==0) then
       write(nfcntn,'(a4)') tag_numk
       write(nfcntn,'(i10)') kv3_ek
       write(nfcntn,'(a14)') tag_convergence_ek
       write(nfcntn,'(40i2)') (iconv_ek(i),i=1,kv3_ek)
    end if
  end subroutine m_CtrlP_wd_iconv_ek

  subroutine m_CtrlP_rd_iconv_ek(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    logical :: EOF_reach, tag_is_found
    integer :: i
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_numk),tag_numk, EOF_reach,tag_is_found,str,len_str)
       if(.not.tag_is_found) then
          numk_tmp = 0
       else
          read(nfcntn,*) numk_tmp
       end if
       write(nfout,'(" numk = ",i8)') numk_tmp
    end if
    if(npes > 1) call mpi_bcast(numk_tmp,1,mpi_integer,0,MPI_CommGroup,ierr)

    if(numk_tmp >= 1) then
       allocate(iconv_ek_tmp(numk_tmp)); iconv_ek_tmp = 0
       if(mype == 0) then
          call rewind_to_tag0(nfcntn,len(tag_convergence_ek),tag_convergence_ek &
               & , EOF_reach, tag_is_found, str,len_str)
          if(.not.tag_is_found) then
             write(nfout,'(" tag_convergence_ek is not found")')
          else
             read(nfcntn,*) (iconv_ek_tmp(i),i=1,numk_tmp)
          end if
          do i = 1, numk_tmp
             if(iconv_ek_tmp(i) < 0 ) iconv_ek_tmp(i) = 0
             if(iconv_ek_tmp(i) > EK_CONVERGED) iconv_ek_tmp(i) = EK_CONVERGED
          end do
       end if
    end if
    if(npes > 1) then
       if(numk_tmp >=1) then
          call mpi_bcast(iconv_ek_tmp,numk_tmp,mpi_integer,0,MPI_CommGroup,ierr)

          if(printable) then
             write(nfout,'(" iconvergence_ek_previous_job =")')
             write(nfout,'(40i2)') (iconv_ek_tmp(i),i=1,numk_tmp)
          end if
          if(icond == FIXED_CHARGE_CONTINUATION .and. neg_previous < neg) then
             iconv_ek_tmp = 0
             if(printable) write(nfout,'(" iconvergence_ek_previous_job is reset " &
                  & ,i2,", because neg_previous < neg")') iconv_ek_tmp(1)
          end if
       end if
    end if
  end subroutine m_CtrlP_rd_iconv_ek

  subroutine m_CtrlP_wd_edelta_ontheway(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,'(a15)') tag_edelta_ontheway
       write(nfcntn,'(d24.16)') edelta_ontheway
    end if
  end subroutine m_CtrlP_wd_edelta_ontheway

  subroutine m_CtrlP_rd_edelta_ontheway(nfcntn)
    integer, intent(in) :: nfcntn
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_edelta_ontheway),tag_edelta_ontheway &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          edelta_ontheway = edelta_initial
!!$          stop ' tag_edelta_ontheway is not found'
       else
          read(nfcntn,*) edelta_ontheway
       end if
    end if
    if(npes > 1) call mpi_bcast(edelta_ontheway,1,mpi_double_precision,0 &
         & ,MPI_CommGroup,ierr)
    if(printable) write(6,'(d24.16, " : edelta_ontheway")') edelta_ontheway
  end subroutine m_CtrlP_rd_edelta_ontheway

  subroutine m_CtrlP_set_corecharge_cntnbin(onoroff)
    integer, intent(in) :: onoroff
    corecharge_cntnbin = onoroff
  end subroutine m_CtrlP_set_corecharge_cntnbin

  subroutine m_CtrlP_wd_corecharge_cntnbin(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,'(a18)') tag_corecharge_cntnbin
       write(nfcntn,'(i8)') corecharge_cntnbin
    end if
  end subroutine m_CtrlP_wd_corecharge_cntnbin

  subroutine m_CtrlP_rd_corecharge_cntnbin(nfcntn)
    integer, intent(in) :: nfcntn
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_corecharge_cntnbin),tag_corecharge_cntnbin &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          corecharge_cntnbin = OFF
       else
          read(nfcntn,*) corecharge_cntnbin
       end if
    end if
    if(npes > 1) call mpi_bcast(corecharge_cntnbin,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(printable) write(6,'(i8, " : corecharge_cntnbin")') corecharge_cntnbin
  end subroutine m_CtrlP_rd_corecharge_cntnbin

  subroutine m_CtrlP_wd_neg(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype==0) then
       write(nfcntn,'(a3)') tag_neg
       write(nfcntn,'(i10)') neg
    end if
  end subroutine m_CtrlP_wd_neg

  subroutine m_CtrlP_rd_neg_previous(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_neg),tag_neg &
            & , EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          neg_previous = neg
!!$          stop ' tag_edelta_ontheway is not found'
       else
          read(nfcntn,*) neg_previous
       end if
    end if
    if(npes > 1) call mpi_bcast(neg_previous,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(printable) write(6,'(i10, " : neg_previous")') neg_previous
! ---> T. Yamasaki, 12 July 2008
    if(neg < neg_previous) then
       neg = neg_previous
       meg = neg
       call m_CntrlP_set_meg(neg)
       call m_CntrlP_set_neg(neg,nfout)
       if(ipri >= 1) write(6,'(" ### Warning(1309): Number of bands is enlarged",i12)') neg
       call m_CntrlP_set_davidson_size(neg) ! -> max_subspace_size
    end if
! <---
  end subroutine m_CtrlP_rd_neg_previous

  subroutine m_CtrlP_renew_edelta_ontheway(nfout,etotal,forcmx,iteration_ionic,iteration,iteration_electronic)
    integer, intent(in) ::       iteration_ionic, iteration, iteration_electronic
    integer, intent(in) ::       nfout
    real(kind=DP), intent(in) :: etotal,forcmx
    real(kind=DP) :: fd, ed, edelta_previous

    if(edelta_initial_is_given) edelta_previous = edelta_ontheway
    if(.not.edelta_initial_is_given) then
       edelta_ontheway = edelta
    else if(dabs(forcmx) < forccr*FORCCR_FACTOR_EDELTA_ONTHEWAY ) then
       edelta_ontheway = edelta
    else if( dabs(forcmx) > max_force_edelta_i) then
       edelta_ontheway = edelta_initial
    else
       fd = max_force_edelta_i - FORCCR_FACTOR_EDELTA_ONTHEWAY*forccr
       ed = edelta_initial - edelta
       if(dabs(fd) < Edelta_Critical_Value) then
          edelta_ontheway = edelta
       else if(dabs(ed) < Edelta_Critical_Value) then
          edelta_ontheway = edelta
       else
! if you want to change the way of deciding edelta_ontheway, change this line.
          edelta_ontheway = (ed/fd)*(dabs(forcmx) - FORCCR_FACTOR_EDELTA_ONTHEWAY*forccr) + edelta
!
       end if
    end if
    if(edelta_initial_is_given) then
       if(edelta_ontheway > edelta_previous) edelta_ontheway = edelta_previous
       if(ipri >=1) then
          write(nfout,'(" !edel  ", i5,i8,i5,2f17.10,",   edelta(now, initial, final) = (",d12.4,2d10.2,")")') &
               & iteration_ionic, iteration, iteration_electronic, etotal, forcmx, edelta_ontheway, edelta_initial, edelta
       end if
    end if
  end subroutine m_CtrlP_renew_edelta_ontheway

  subroutine m_CtrlP_edelta_for_sampling(nfout, iteration_ionic)
    integer, intent(in) :: nfout, iteration_ionic
    integer :: iteration1
    real(kind=DP), save :: edelta_buf = -1
    iteration1 = iteration_ionic+1
    if(mod(iteration1,frequency_nnp) == 0) then
       if(edelta_buf<0) edelta_buf = edelta_ontheway
       edelta_ontheway = edelta_sampling
       if(printable) write(nfout,'(a,e12.5)') ' !** changed edelta ' &
                              & , edelta_ontheway
    endif
    if(mod(iteration_ionic,frequency_nnp) == 0) then
       if(edelta_buf<0) edelta_buf = edelta
       edelta_ontheway = edelta_buf
       if(printable) write(nfout,'(a,e12.5)') ' !** retrieved  ' &
                              & , edelta_ontheway
    endif
  end subroutine m_CtrlP_edelta_for_sampling

  logical function m_CtrlP_edelta_final()
    if(edelta_ontheway <= edelta+Edelta_Critical_Value) then
       m_CtrlP_edelta_final = .true.
    else
       m_CtrlP_edelta_final = .false.
    end if
  end function m_CtrlP_edelta_final

  subroutine m_CtrlP_reset_edelta_ontheway()
    edelta_ontheway = edelta_initial
  end subroutine m_CtrlP_reset_edelta_ontheway

  subroutine m_CtrlP_get_edelta(edelta_transfer)
    real(kind=DP), intent(out) :: edelta_transfer
    edelta_transfer = edelta_ontheway
  end subroutine m_CtrlP_get_edelta

!!$  subroutine m_CtrlP_rd_nrsv(nfcntn)
!!$    integer, intent(in)       :: nfcntn
!!$    logical    :: EOF_reach, tag_is_found
!!$
!!$    if(mype==0) then
!!$       call rewind_to_tag0(nfcntn,len_tag_T_cntrl,tag_T_cntrl &
!!$            &, EOF_reach, tag_is_found,str,len_str)
!!$       if(.not.tag_is_found) then
!!$          call rewind_to_tag0(nfcntn,len_tag_T_cntrl,tag_T_cntrl2&
!!$               &, EOF_reach, tag_is_found, str,len_str)
!!$          if(.not.tag_is_found) then
!!$             stop ' tag_T_cntrl is not found'
!!$          end if
!!$       end if
!!$       read(nfcntn,*)
!!$       read(nfcntn,*) nrsv
!!$    endif
!!$    if(npes > 1) call mpi_bcast(nrsv,1,mpi_integer,0,MPI_CommGroup,ierr)
!!$    write(6,'(i5, " : nrsv")') nrsv
!!$  end subroutine m_CtrlP_rd_nrsv
!!$
!!$  subroutine m_CtrlP_wd_nrsv(nfcntn)
!!$    integer, intent(in)       :: nfcntn
!!$    if(mype==0) then
!!$       write(nfcntn,*) tag_T_cntrl
!!$       write(nfcntn,'(" -- nrsv --")')
!!$       write(nfcntn,'(i10)') nrsv
!!$    endif
!!$  end subroutine m_CtrlP_wd_nrsv
!!$
!!$  subroutine m_CtrlP_rd_nrsv_stdin(nfinp)
!!$    integer, intent(in)       :: nfinp
!!$
!!$    logical :: eof_reach, tag_is_found
!!$    real(kind=DP) :: wk(1)
!!$
!!$    if(imdalg /= T_CONTROL .and. imdalg /= BLUEMOON) return
!!$    call rewind_to_tag0(nfinp,len_tag_T_cntrl, tag_T_cntrl &
!!$         &, EOF_reach, tag_T_cntrl_is_found, str, len_str)
!!$    if(.not.tag_T_cntrl_is_found) &
!!$         & call rewind_to_tag0(nfinp,len_tag_T_cntrl, tag_T_cntrl2 &
!!$         &, EOF_reach, tag_T_cntrl_is_found, str, len_str)
!!$    write(6,'(" -- after rewind_to_tag0 --")')
!!$
!!$    if(.not.tag_T_cntrl_is_found) then
!!$       write(6,'(" ! NO TAG of (TEMPERATURE_CONTROL)")')
!!$       return
!!$    end if
!!$    write(6,'(" -tag_T_cntrl_is_found-")')
!!$ !! after rewind_to_tag0 --'
!!$
!!$    nrsv = 1  ! nrsv: number of heat bath
!!$    tag_is_found = .false.
!!$    do while(.not.tag_is_found)
!!$       read(nfinp,'(a132)') str
!!$       call strncmp2(str,len_str,tag_nrsv,len_tag_nrsv,tag_is_found)
!!$       if(tag_is_found) then
!!$          call read_RHS_number(str,len_str,wk,1)
!!$          nrsv = nint(wk(1))
!!$       end if
!!$    end do
!!$    if(nrsv < 1) nrsv = 1
!!$    call strcpy(tag_T_cntrl,len_tag_T_cntrl,str,60)
!!$    print '(a60)', str(1:60)
!!$    print *, ' nrsv = ', nrsv
!!$  end subroutine m_CtrlP_rd_nrsv_stdin

#ifndef _EMPIRICAL_
  function m_CtrlP_decide_dtim_1Dsearch(nfout,etot_trial,dtim1,factor)
!     Original subroutine name was "set_dtim0",
!         which was coded by T. Sanada
!                           @(#)set_dtim.f 9.4 01/11/14 17:19:34
!!$cTS**************************************************
!!$c                                       '96.01.24
!!$c 1.  Time_factors for sub iteration are set.
!!$c 2.  dtim are modified accroding to the previous
!!$c     'optimized-time-step'(dtim_old).
!!$c 3.  dtim_msdv for msdv is set.
!!$cTS**************************************************
!      Translated into this function coded with fortran90
!                                    by T. Yamasaki in 1999.
!
    real(kind=DP)             :: m_CtrlP_decide_dtim_1Dsearch
    integer, intent(in)       :: nfout
    real(kind=DP), intent(in) :: etot_trial(3),dtim1,factor

    real(kind=DP)   :: tf1, tf2, e1, e2, denom, tff, fac1, afactor, dtim_msdv, eforecast
    integer         :: type_ls, i, lmm_status
!!$    real(kind=DP), parameter :: Delta = 1.d-10
    real(kind=DP), parameter :: Rvmxdf = 100.d0
    integer, save :: title_label = -1
!!$    real(kind=DP), parameter :: dt_Lower_CRITICAL = 5.d-2
!!$    real(kind=DP), parameter :: dt_Upper_CRITICAL = 2.d0

    tf1 = 1; tf2 = factor
!!$    do i = 1, 3
!!$       if(etot_trial(i) > LIMIT_1DSEARCH) then
!!$          if(ipri>=1) then
!!$             write(nfout,'(" etot_trial(",i2,") is NaN <<m_CtrlP_decide_dtim_1Dsearch>>")') i
!!$          end if
!!$          stop ' etot_trial is NaN <<m_CtrlP_decide_dtim_1Dsearch>>'
!!$       end if
!!$    end do

    ! assuming an equation of
    !     e = afactor*((t-tff)**2 - tff*tff) = afactor*t*(t-2*tff), --(1)
    ! then
    !     e1 = afactor*tf1*(tf1-2*tff)         --(2)
    !     e2 = afactor*tf2*(tf2-2*tff)         --(3)
    ! afactor and tff are solved from (1) and (2) as
    !     afactor = (tf2*e1-tf1*e2)/(tf1*tf2*(tf1-tf2))             --(4)
    !         tff = 0.5 * (e1*tf2*tf2 - e2*tf1*tf1)/(tf2*e1-tf1*e2) --(5)
    !

    e1  = etot_trial(2) - etot_trial(1)
    e2  = etot_trial(3) - etot_trial(1)
    denom = tf2*e1 - tf1*e2

! --> T. Yamasaki, 22nd July 2008
    if(e1 > de_critical_1dsrch .and. e2 > de_critical_1dsrch) then
! ================ modified by K. Tagami ==== trial === 0.2
!!!!!!       incre_etot_in_1dsrch = incre_etot_in_1dsrch + 1
        if ( sw_calc_force == OFF ) then
          incre_etot_in_1dsrch = incre_etot_in_1dsrch + 1
        endif
! ===================================================== 0.2
    else
       incre_etot_in_1dsrch = 0
    end if
! <-

! =================== added by K. Tagami =============== 11.0
    dtim_msdv = 0.0d0
! ====================================================== 11.0

! ===================== KT_Test ====================== 12.5Exp
    type_ls = 0
    lmm_status = yes
! =================================================== 12.5Exp

    if(dabs(denom) > delta_lmdenom) then
       tff = 0.5 * (e1*tf2*tf2 - e2*tf1*tf1)/denom
    else
       fac1 = 1.d0
       tff = sign(fac1, e1*tf2*tf2 - e2*tf1*tf1)
       tff = tff * Rvmxdf
    end if

    afactor = (tf2*e1-tf1*e2)/(tf1*tf2*(tf1-tf2))


    if(afactor > 0.0d0) then
       if(tff.gt.0.0d0.and.tff.le.factor*2) then

          dtim_msdv = dtim1 * tff
! ------------------- Revised by T. Yamasaki, 22 July 2008 ---
          if(dtim_msdv < dt_Lower_CRITICAL*dt_lower_factor) then
             dtim_msdv = dt_Lower_CRITICAL*dt_lower_factor
             if(dabs(dtim1) < DELTA10) then
                tf1 = DELTA10
             else
                tf1 = dtim_msdv/dtim1
             end if
             eforecast = afactor*tf1*(tf1-2*tff) + etot_trial(1)
             lmm_status = NO    !  T. Yamasaki, 18th Aug. 2009
          else
             eforecast = -tff*tff*afactor + etot_trial(1)
             lmm_status = YES   !  T. Yamasaki, 18th Aug. 2009
          end if
! -----------------------------------------------------------<<
          type_ls = 1

       else if(tff.le.0.0d0) then

#ifdef LMM_PREVIOUS
          dtim_msdv = dtim1 * 0.1d0
#else
! ------------------- Revised by T. Yamasaki, 09 July 2008 ---
! ------------------- Revised by T. Yamasaki, 28 June 2008 ---
!!$          dtim_msdv = dtim1 * 0.1d0
          if(e2 < 0.2) then
             dtim_msdv = dtim1 * tff
          else
             dtim_msdv = dtim1
          end if
! ------------------------------------------------------------<<
! ------------------------------------------------------------<<
! ------------------- Revised by T. Yamasaki, 28 June 2008 ---
! ------------------- Revised by T. Yamasaki, 09 July 2008 ---
!!$          if(dabs(dtim_msdv) > dt_Upper_CRITICAL) dtim_msdv = -dt_Upper_CRITICAL
          if(dtim_msdv >  dt_Upper_CRITICAL) dtim_msdv =  dt_Upper_CRITICAL
          if(dtim_msdv < -dt_Upper_CRITICAL) dtim_msdv = -dt_Upper_CRITICAL
          if(dtim_msdv < 0.d0 .and. dtim_msdv > -dt_Lower_CRITICAL*dt_lower_factor) dtim_msdv = -dt_Lower_CRITICAL*dt_lower_factor
#endif
! --> T. Yamasaki,  16th July 2008
          if(dabs(dtim1) < DELTA10) then
             tf1 = DELTA10
          else
             tf1 = dtim_msdv/dtim1
          end if
          eforecast = afactor*tf1*(tf1-2*tff) + etot_trial(1)
!!$! ------------------------------------------------------------<<
!!$          eforecast = -tff*tff*afactor + etot_trial(1)
!!$! ------------------------------------------------------------<<
          type_ls = 2
          lmm_status = NO    !  T. Yamasaki, 18th Aug. 2009
       else if(tff.gt.factor*2) then
          dtim_msdv = dtim1 * tff
! --> T. Yamasaki,  15th July 2008
          if(dt_upper_factor_is_set) then
             if(dtim_msdv >  dt_Upper_CRITICAL*dt_Upper_factor) dtim_msdv =  dt_Upper_CRITICAL*dt_Upper_factor
             if(dabs(dtim1) < DELTA10) then
                tf1 = DELTA10
             else
                tf1 = dtim_msdv/dtim1
             end if
             eforecast = afactor*tf1*(tf1-2*tff) + etot_trial(1)
          else
             eforecast = -tff*tff*afactor + etot_trial(1)
          end if
! <--
          type_ls = 3
          lmm_status = YES    !  T. Yamasaki, 18th Aug. 2009
       end if
    else

!!$       if( etot_trial(1).le.etot_trial(2).and. &
!!$            &          etot_trial(1).le.etot_trial(3) ) then
       if( etot_trial(1).le.etot_trial(2).and. &
            &          etot_trial(2).le.etot_trial(3) ) then

          dtim_msdv = dtim1 * 0.1d0
! ------------------- Revised by T. Yamasaki, 28 June 2008 ---
!!$          if(dtim_msdv < dt_Lower_CRITICAL) dtim_msdv = dt_Lower_CRITICAL
! -----------------------------------------------------------<<
          type_ls = 4
          lmm_status = NO    !  T. Yamasaki, 18th Aug. 2009
       else if( etot_trial(2).ge.etot_trial(1).and. &
            &               etot_trial(2).ge.etot_trial(3) ) then
#ifdef LMM_PREVIOUS_V700
          dtim_msdv = dtim1 * factor
          if(e2 > 0.0 .and. e1 > 0.0) then
             lmm_status = NO
          else
             lmm_status = YES    !  T. Yamasaki, 18th Aug. 2009
          end if
#else
! ------------------- Revised by T. Yamasaki, 03 July 2008 ---
          if(e2 > 0.0 .and. e1 > 0.0) then
             dtim_msdv = dt_Lower_CRITICAL
             lmm_status = NO    !  T. Yamasaki, 18th Aug. 2009
          else
             dtim_msdv = dtim1 * factor
             lmm_status = YES    !  T. Yamasaki, 18th Aug. 2009
          end if
! -----------------------------------------------------------<<
#endif
          type_ls = 5
       else if( etot_trial(3).le.etot_trial(1).and. &
            &               etot_trial(3).le.etot_trial(2) ) then

          dtim_msdv = dtim1 * factor
          type_ls = 6
          lmm_status = YES    !  T. Yamasaki, 18th Aug. 2009
       end if

! ------------------- Revised by T. Yamasaki, 18 July 2008 ---
       if(dabs(dtim1) < DELTA10) then
          tf1 = DELTA10
       else
          tf1 = dtim_msdv/dtim1
       end if
       eforecast = afactor*tf1*(tf1-2*tff) + etot_trial(1)
! --------------------------------- <<
    end if
!!$    if(dtim_msdv.lt.dt_Lower_CRITICAL) dtim_msdv = dt_Lower_CRITICAL
!!$    if(dtim_msdv.gt.dt_Upper_CRITICAL) dtim_msdv = dt_Upper_CRITICAL
    m_CtrlP_decide_dtim_1Dsearch = dtim_msdv
    dtim_1Dsearch        = dtim_msdv
    if(dtim_1Dsearch < dt_Lower_CRITICAL) dtim_1Dsearch = dt_Lower_CRITICAL
    if(dtim_1Dsearch > dt_Upper_CRITICAL) dtim_1Dsearch = dt_Upper_CRITICAL

    if(type_ls == 1) then
       nitersub = 1
    else if(type_ls == 3) then
       nitersub = 1
    else
       nitersub = 0
    endif

! --> T. Yamasaki 17th Aug. 2009
    if(sw_lmm_status_check == YES) then
       lmm_status_pointer = lmm_status_pointer+1
       lmm_status_stored = lmm_status_stored+1
       if(lmm_status_pointer > lmm_status_store_size) lmm_status_pointer = 1
       lmm_status_store(lmm_status_pointer) = lmm_status
    end if
! <--

!!$cTS**
!!$c  nitersub = 0: abnormal,
!!$c  nitersub = 1: normal ( dtim is determined so that we have the local min.
!!$cTS**

    if(ipri >= 2) then
       write(nfout,*) '***'
       write(nfout,'(" etot (0,1,2)= :",3f20.12)') (etot_trial(i),i=1,3)
       write(nfout,'(" tff         = :",f20.12,",  dt_new = :",f20.12)') tff,dtim_msdv
       write(nfout,'(" type_ls     = :",i5, ", nitersub = :",i5,&
            &", afactor = :",f20.12)') type_ls, nitersub, afactor
       write(nfout,*) '***'
    else if(ipri >= 1) then
!!$       write(nfout,'(" !1Dsrch etot(:)= ",3f16.8," dt_new = ",f8.4)')&
!!$            & (etot_trial(i),i=1,3), dtim_msdv
       if(title_label == -1) then
          write(nfout,'(" !1Dsrch    etot(1)       etot(2)-(1)  etot(3)-(1)     dt     dtn    ls    lmm   eforecast")')
          title_label = 0
       end if
! >>------------------ Revised by T. Yamasaki, 18 July 2008 ---
! ->>----------------- Revised by T. Yamasaki, 28 June 2008 ---
!!$       if(type_ls < 4) then
!!$       success_or_fault = 1
!!$       if(eforecast < etot_trial(1)) then
!!$          success_or_fault = -1
!!$       else
!!$          success_or_fault = 1
!!$       end if
       if(abs(etot_trial(3)-etot_trial(1)) >= 1000 .or. abs(etot_trial(2)-etot_trial(1)) >= 1.000) then
          write(nfout,'(" !1Dsrch ",f15.7,2f13.4,"  ",2f8.3,2i3,f20.12,i4)') &
               & etot_trial(1), e1, e2, dtim1, dtim_msdv, type_ls, lmm_status, eforecast, incre_etot_in_1dsrch
       else
          write(nfout,'(" !1Dsrch ",f15.7,2f13.8,"  ",2f8.3,2i3,f20.12,i4)') &
               & etot_trial(1), e1, e2, dtim1, dtim_msdv, type_ls, lmm_status, eforecast, incre_etot_in_1dsrch
       end if
!!$       else
!!$          if(abs(etot_trial(3)-etot_trial(1)) >= 1000 .or. abs(etot_trial(2)-etot_trial(1)) >= 1.000) then
!!$             write(nfout,'(" !1Dsrch ",f13.7,2f13.4,"   ",2f8.3,i3)') &
!!$                  & etot_trial(1), e1, e2, dtim1, dtim_msdv, type_ls
!!$          else
!!$             write(nfout,'(" !1Dsrch ",f13.7,2f13.8,"   ",2f8.3,i3)') &
!!$                  & etot_trial(1), e1, e2, dtim1, dtim_msdv, type_ls
!!$          end if
!!$       end if
! -----------------------------<<
! ----------------------------- <<
    end if

  end function m_CtrlP_decide_dtim_1Dsearch

  logical function m_CtrlP_etot_1dsrch_divergent(iter_elec)
    integer, intent(in) :: iter_elec
    integer :: incre, i
    if(min(incre_etot_in_1dsrch,iter_elec) >= incre_etot_in_1dsrch_limit) then
       m_CtrlP_etot_1dsrch_divergent = .true.
    else
! --> T. Yamasaki 17th Aug. 2009
       if(sw_lmm_status_check == YES) then
          if(lmm_status_stored >= lmm_status_store_size) then
             m_CtrlP_etot_1dsrch_divergent = .true.
             do i = 1, lmm_status_store_size
                if(lmm_status_store(i) == 1) then
                   m_CtrlP_etot_1dsrch_divergent = .false.
                   exit
                end if
             end do
          else
             m_CtrlP_etot_1dsrch_divergent = .false.
          end if
          if(m_CtrlP_etot_1dsrch_divergent) call m_CtrlP_clear_lmm_status_store()
       else
          m_CtrlP_etot_1dsrch_divergent = .false.
       end if
    end if
! <--
  end function m_CtrlP_etot_1dsrch_divergent

! --> T. Yamasaki, 18th Aug. 2009
  subroutine m_CtrlP_clear_lmm_status_store
    lmm_status_store = 1
    lmm_status_pointer = 1
    lmm_status_stored = 0
  end subroutine m_CtrlP_clear_lmm_status_store
! <--

  function m_CtrlP_dtim_1Dsearch_now(dtim)
    real(kind=DP)             :: m_CtrlP_dtim_1Dsearch_now
    real(kind=DP), intent(in) :: dtim
    if(dtim_1Dsearch < 0)  dtim_1Dsearch = dtim
    m_CtrlP_dtim_1Dsearch_now = dtim_1Dsearch
  end function m_CtrlP_dtim_1Dsearch_now

  logical function m_CtrlP_dtim_1Dsearch_is_neg()
    if(dtim_1Dsearch < 0 ) then
       m_CtrlP_dtim_1Dsearch_is_neg = .true.
    else
       m_CtrlP_dtim_1Dsearch_is_neg = .false.
    end if
  end function m_CtrlP_dtim_1Dsearch_is_neg
#endif

  subroutine m_CtrlP_reset_dtim_1Dsearch
    dtim_1Dsearch = -1
    nitersub = 0
  end subroutine m_CtrlP_reset_dtim_1Dsearch

  subroutine m_CtrlP_rd_istop(nfstop)
    integer, intent(in) :: nfstop
    if(mype == 0) then
       rewind nfstop
       istop = -1
       read(nfstop,*,end = 1, err=  1) istop
1      continue
       if(istop .ne. -1 .and. ipri>=1) write(6,'(" istop = ",i6,"<<m_CtrlP_rd_istop>>")') istop
    endif
    if(npes > 1) call mpi_bcast(istop,1,mpi_integer,0,MPI_CommGroup,ierr)
  end subroutine m_CtrlP_rd_istop

  subroutine m_CtrlP_set_wct_start
#ifdef _CHECK_ELAPSE_
    integer :: iremain
    call chkelaps(iremain)
    wct_start = iremain
#else
    call gettod(wct_start)
#endif
  end subroutine m_CtrlP_set_wct_start

  function m_CtrlP_get_elpsd_time() result(res)
    real(kind=DP) :: res
    real(kind=DP)       :: wct_now, cpu_total
#ifdef _CHECK_ELAPSE_
    integer :: iremain
    call chkelaps(iremain)
    wct_now = iremain
    cpu_total = wct_start - wct_now
#else
    call gettod(wct_now)
    cpu_total = (wct_now - wct_start) * UMICRO
#endif
    if(npes > 1) call mpi_bcast(cpu_total,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    res = cpu_total
    return
  end function m_CtrlP_get_elpsd_time

  logical function m_CtrlP_ckcput()
    integer             :: istats
    real(kind=DP)       :: wct_now, cpu_total
#ifdef _CHECK_ELAPSE_
    integer :: iremain
    call chkelaps(iremain)
    wct_now = iremain
    cpu_total = wct_start - wct_now
#else
    call gettod(wct_now)
    cpu_total = (wct_now - wct_start) * UMICRO
#endif

    if(npes > 1) call mpi_bcast(cpu_total,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    if(cpu_total > cpumax) then
       istats = FINISH
    else
       istats = 0
    end if
#ifdef _CHECK_ELAPSE_
    if(npes > 1) call mpi_bcast(wct_now,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    if(wct_now < Critical_Remaining_CPU_TIME) istats = FINISH
#endif
    if(istats == FINISH) then
       m_CtrlP_ckcput = .true.
    else
       m_CtrlP_ckcput = .false.
    end if
    if(m_CtrlP_ckcput) terminated_because = MAX_TIME_REACHED
  end function m_CtrlP_ckcput

  subroutine m_CtrlP_wd_cpu_total(nfout)
    integer, intent(in), optional :: nfout
    real(kind=DP)       :: wct_now, cpu_total
    integer             :: nfo
#ifdef _CHECK_ELAPSE_
    integer :: iremain
    call chkelaps(iremain)
    wct_now = iremain
    cpu_total = wct_start - wct_now
#else
    call gettod(wct_now)
    cpu_total = (wct_now - wct_start) * UMICRO
#endif
    nfo = 6
    if(present(nfout)) nfo = nfout
    if(npes > 1) call mpi_bcast(cpu_total,1,mpi_double_precision,0,MPI_CommGroup,ierr)
    if(printable) write(nfo,'(" <<Total CPU TIME (elapsed time) = ",f15.5, " (sec.)>>")') cpu_total
  end subroutine m_CtrlP_wd_cpu_total

  subroutine m_CtrlP_set_paramset_on
    paramset = .true.
  end subroutine m_CtrlP_set_paramset_on

  subroutine m_CtrlP_set_paramset_off
    paramset = .false.
  end subroutine m_CtrlP_set_paramset_off

#ifndef _EMPIRICAL_
  integer function m_CtrlP_way_of_smearing()
!!$    if(width > 0.d0) then
!!$       m_CtrlP_way_of_smearing = PARABOLIC
!!$    else if(width > -10.d0) then
!!$       m_CtrlP_way_of_smearing = MP
!!$    else
!!$       m_CtrlP_way_of_smearing = TETRAHEDRON
!!$       if(way_ksample == FILE .or. way_ksample == SKPS_DIRECT_IN &
!!$            & .or. way_ksample == GAMMA) then
!!$          m_CtrlP_way_of_smearing = PARABOLIC
!!$       end if
!!$    end if

    m_CtrlP_way_of_smearing = way_of_smearing
    if(way_ksample == FILE .or. way_ksample == SKPS_DIRECT_IN &
         & .or. way_ksample == MONKHORST_PACK &
         & .or. way_ksample == GAMMA) then
       m_CtrlP_way_of_smearing = PARABOLIC

       if(way_of_smearing == COLD) m_CtrlP_way_of_smearing = COLD
! =============== KT_add ================================= 13.0E
       if(way_of_smearing == FERMI_DIRAC ) m_CtrlP_way_of_smearing = FERMI_DIRAC
! ======================================================== 13.0E
       if(way_of_smearing == MP ) m_CtrlP_way_of_smearing = MP
       if(way_of_smearing == LOWEST_AT_EACH_KPT ) &
            &                     m_CtrlP_way_of_smearing = LOWEST_AT_EACH_KPT
    end if

  end function m_CtrlP_way_of_smearing

  integer function m_CtrlP_waymix_now(iteration_electronic,iteration_ionic &
       &                             ,edeltb_per_atom,mixer_changed)
    integer, intent(in)      :: iteration_electronic, iteration_ionic
    real(kind=DP),intent(in) :: edeltb_per_atom
    logical, intent(out), optional :: mixer_changed
    integer :: ips, ipe, imsd_t, it, ip, i

    if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
       if(iteration_ionic >= 2) then
          if(n_WF_solvers_after == 0) then
             ips = 1
          else
             ips = n_WF_solvers_before + 1
          end if
          ipe = n_WF_solvers_all
       else
          ips = 1
          ipe = n_WF_solvers_before
       end if
       if(ipe < ips) call phase_error_with_msg(6,' n_WF_solvers_(before|after) are illegal (m_CtrlP_solver_for_WFs_now)'&
                                              ,__LINE__,__FILE__)
       do i = ips, ipe
          if(w_solver(i)%till_n_iter < 0) exit
          if(iteration_electronic <= w_solver(i)%till_n_iter) exit
       end do
       it = i
       if(it > ipe) it = ipe
       imsd_t = w_solver(it)%solver
       solver : select case(imsd_t)
          case (RMM, RMM2, RMM2P)
             if(imsd_t /= previous_solver) it = it -1
!!$             if(dabs(edeltb_per_atom) > edelta_change_to_rmm ) it = it - 1
          case default
       end select solver
       if(it < ips) it = ips
       ip_w_solver = it
       ip = w_solver(ip_w_solver)%cmix_pointer
!       if(ip < 0 .or. ip > n_Charge_Mixing_way) stop ' error at m_CtrlP_waymix_now'
       if(ip < 0 .or. ip > n_Charge_Mixing_way) call phase_execution_error(INVALID_CHARGE_MIXING)
       m_CtrlP_waymix_now = charge_mixing(ip)%mixing_way
       if(present(mixer_changed)) then
         mixer_changed = ip_charge_mixing /= ip
       endif
       ip_charge_mixing = ip
    else
       m_CtrlP_waymix_now = waymix
    end if
  end function m_CtrlP_waymix_now

  subroutine m_CtrlP_set_mix_parameter()
    if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
       waymix = charge_mixing(ip_charge_mixing)%mixing_way
       if(charge_mixing(ip_charge_mixing)%precon == YES) then
          c_precon = .true.
       else
          c_precon = .false.
       end if
       istrbr = charge_mixing(ip_charge_mixing)%istr
       nbxmix = charge_mixing(ip_charge_mixing)%nbxmix
       if(waymix == BROYD1 .or. waymix == BROYD2 .or. waymix == DFP .or. waymix == PULAY) then
          if(nbxmix <= 2) then
             if(printable) then
                write(6,'(" !! [nbxmix] given is ",i6, " waymix = ",i6)') nbxmix,waymix
                write(6,'(" !! [nbxmix] should be larger than 2")')
                write(6,'(" !! [nbxmix] is set to be 3")')
             end if
             nbxmix = 3
          end if
       end if
       hownew = charge_mixing(ip_charge_mixing)%hownew
       cutoff_mix = charge_mixing(ip_charge_mixing)%cutoff
    else
       istrbr = istrbr_b
       nbxmix = nbxmix_b
       hownew = hownew_b
       cutoff_mix = cutoff_mix_b
    end if
  end subroutine m_CtrlP_set_mix_parameter

  function m_CtrlP_rmx_now(iteration_electronic, iteration_ionic)
    real(kind=DP) :: m_CtrlP_rmx_now
    integer, intent(in) :: iteration_electronic, iteration_ionic
    real(kind=DP) :: rmxs, rmxe, rmxt
    integer       :: iter_range, itr_p, ipp, ips

    if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
       rmxs = charge_mixing(ip_charge_mixing)%rmxs
       rmxe = charge_mixing(ip_charge_mixing)%rmxe
       iter_range = charge_mixing(ip_charge_mixing)%iter_range
       ipp = ip_w_solver -1
       if(iteration_ionic >= 2) then
          if(n_WF_solvers_after == 0) then
             ips = 1
          else
             ips = n_WF_solvers_before+1
          end if
       else
          ips = 1
       end if
       !!$if(ipp < ips) ipp = ips
       !!$itr_p = w_solver(ipp)%till_n_iter
       if(ipp < ips) then
          itr_p = 1
       else
          itr_p = w_solver(ipp)%till_n_iter
       end if
       if(charge_mixing(ip_charge_mixing)%variation_way == varLINEAR) then
          rmxt = set_dx_linear(rmxs,rmxe,iter_range,iteration_electronic-itr_p)
       else
          rmxt = set_dx_tanh(rmxs,rmxe,iter_range,iteration_electronic-itr_p)
       end if
    else
       if(iexpl > MDOLD) then
          call decide_rmx_case1
       else
          call decide_rmx_case2
       end if
    end if
    m_CtrlP_rmx_now = rmxt
  contains
    subroutine decide_rmx_case1
      if(iteration_ionic > 1) then
         rmxt = rmx_1234(4)
      else
         if(iteration_electronic < max_scf_iteration/4) then
            rmxt = rmx_1234(1)
         else if(iteration_electronic < max_scf_iteration/4*2) then
            rmxt = rmx_1234(2)
         else if(iteration_electronic < max_scf_iteration/4*3) then
            rmxt = rmx_1234(3)
         else
            rmxt = rmx_1234(4)
         end if
      end if
    end subroutine decide_rmx_case1

    subroutine decide_rmx_case2
      integer :: iter_crc
      real(kind=DP) :: rmxs, rmxd, xnmd1, fex, fac
      integer       :: nmd21
      integer, parameter :: COEFF = 10

      if(iteration_ionic == 1) then
         rmxs = rmx_1234(1)
         rmxd = rmx_1234(2) - rmx_1234(1)
         xnmd1 = max_scf_iteration
         fex  = (dble(iteration_electronic)-xnmd1*0.5)/xnmd1
      else
         rmxs = rmx_1234(3)
         rmxd = rmx_1234(4) - rmx_1234(3)
         iter_crc = iteration_electronic
         if(iteration_electronic < max_scf_iteration) &
              &iter_crc = iteration_electronic + max_scf_iteration
         if(max_total_scf_iteration < max_scf_iteration*COEFF) then
            nmd21 = max_total_scf_iteration - max_scf_iteration
         else
            nmd21 = max_scf_iteration*COEFF
         end if
         if(nmd21 == 0) then
            fex = iter_crc
         else
            fex = (dble(iter_crc) - nmd21 * 0.5)/nmd21
         end if
      end if
      fac = dtanh(COEFF*fex)*0.5 + 0.5
      rmxt = rmxs + fac * rmxd
      if(rmxt > 1.d0) rmxt = 1.d0
    end subroutine decide_rmx_case2
  end function m_CtrlP_rmx_now
#endif
! endif of "_EMPIRICAL_"

  real(kind=DP) function set_dx_linear(dxs,dxe,itr_range,itr)
    real(kind=DP), intent(in) :: dxs,dxe
    integer, intent(in)       :: itr_range, itr
    real(kind=DP) :: dx

    if(itr > itr_range) then
       dx = dxe
    else if(itr < 1) then
       dx = dxs
    else
       dx = (dxs*(itr_range-itr) + dxe*itr)/dble(itr_range)
    end if
    set_dx_linear = dx
  end function set_dx_linear

  real(kind=DP) function set_dx_tanh(dxs,dxe,itr_range,itr)
    real(kind=DP), intent(in) :: dxs,dxe
    integer, intent(in)       :: itr_range, itr
    integer, parameter :: COEFF = 10
    real(kind=DP) :: dx,fex,fac

    if(itr > itr_range) then
       dx = dxe
    else if(itr < 1) then
       dx = dxs
    else
       fex    = (itr - itr_range*0.5)/itr_range
       fac    = dtanh(COEFF*fex)*0.5 + 0.5
       dx   = dxs + fac * (dxe-dxs)
    end if
    set_dx_tanh = dx
  end function set_dx_tanh

#ifndef _EMPIRICAL_
  function m_CtrlP_dtim_now(iteration_electronic,iteration_ionic)
    real(kind=DP) :: m_CtrlP_dtim_now
    integer, intent(in) :: iteration_electronic, iteration_ionic
    real(kind=DP) :: dtim, dts, dte
    integer :: i, ips, ipe,is, itr_p, itr_range

    if(tag_solver_of_WF_is_found) then
       if(iteration_ionic >= 2) then
          if(n_WF_solvers_after == 0) then
             ips = 1
          else
             ips = n_WF_solvers_before+1
          end if
          ipe = n_WF_solvers_all
       else
          ips = 1
          ipe = n_WF_solvers_before
       end if
       itr_p = 0
       do i = ips, ipe
          if(w_solver(i)%till_n_iter <= 0) exit
          if(iteration_electronic <= w_solver(i)%till_n_iter) exit
          itr_p = w_solver(i)%till_n_iter
       end do
       if(i > ipe) i = ipe
       is = w_solver(i)%solver
       dts = w_solver(i)%dtim_s
       if(is == lmMSD .or. is == lmSD .or. is == CG .or. is == lmCG &
            & .or. is == RMM2 .or. is == RMM3 .or. is == RMM2P) then
          dtim = dts
       else
          dte = w_solver(i)%dtim_e
          itr_range = w_solver(i)%iter_range
          if(w_solver(i)%variation_way == varLINEAR) then
             dtim = set_dx_linear(dts,dte,itr_range,iteration_electronic-itr_p)
          else
             dtim = set_dx_tanh(dts,dte,itr_range,iteration_electronic-itr_p)
          end if
       end if
    else
       if(iexpl > MDOLD) then
          call decide_dtim_case1
       else if(iexpl == MDSMPL) then
          call decide_dtim_case2
       else
          call decide_dtim_case3
       end if
    end if
    m_CtrlP_dtim_now = dtim
  contains
    subroutine decide_dtim_case1
      if(iteration_ionic > 1) then
         dtim = dtim_1234(4)
      else
         if(iteration_electronic < max_scf_iteration/4) then
            dtim = dtim_1234(1)
         else if(iteration_electronic < max_scf_iteration/4*2) then
            dtim = dtim_1234(2)
         else if(iteration_electronic < max_scf_iteration/4*3) then
            dtim = dtim_1234(3)
         else
            dtim = dtim_1234(4)
         end if
      end if
    end subroutine decide_dtim_case1

    subroutine decide_dtim_case2
      integer :: nmhf, nm34
      real(kind=DP) :: dlast, dtims, dtimxd, fex, fac, xiter
      integer, parameter :: ISTAGE = 15, COEFF = 10

      if(iteration_ionic == 1) then
         nmhf = max_scf_iteration/2
         nm34 = max_scf_iteration*0.75
         dlast = dtim_1234(2)/2
         if(iteration_electronic <= ISTAGE) then
            dtim = dtim_1234(1)
         else if(iteration_electronic <= nmhf) then
            dtim = (dtim_1234(2)-dtim_1234(1))/(nmhf-ISTAGE)
            dtim = dtim*(iteration_electronic - ISTAGE) + dtim_1234(1)
         else if(iteration_electronic <= nm34) then
            dtim = dtim_1234(2)
         else if(iteration_electronic <= max_scf_iteration) then
            dtim = (dlast - dtim_1234(2))/(max_scf_iteration-nm34)
            dtim = dtim*(iteration_electronic - nm34) + dtim_1234(2)
         else
            dtim = dlast
         endif
      else
         dtims  = dtim_1234(3)
         dtimxd = dtim_1234(4) - dtim_1234(3)
         if(iteration_electronic <= ISTAGE) then
            xiter = iteration_electronic
         else
            xiter = 2*ISTAGE - iteration_electronic
         end if
         if(xiter < 1.d0) xiter = 1.d0
         fex    = (xiter - ISTAGE*0.5)/ISTAGE
         fac    = dtanh(COEFF*fex)*0.5 + 0.5
         dtim   = dtims + fac * dtimxd
      end if
    end subroutine decide_dtim_case2

    subroutine decide_dtim_case3
      real(kind=DP) :: dtims, dtimxd, fex, fac, xiter
      integer, parameter :: ISTAGE = 15, COEFF = 10

      if(iteration_ionic == 1) then
         dtims  = dtim_1234(1)
         dtimxd = dtim_1234(2) - dtim_1234(1)
         if(iteration_electronic <= max_scf_iteration) then
            xiter = iteration_electronic
         else
            xiter = 2*max_scf_iteration - iteration_electronic
         end if
         if(xiter < 1.d0) xiter = 1.d0
         fex = (xiter - max_scf_iteration*0.5)/max_scf_iteration
      else
         dtims  = dtim_1234(3)
         dtimxd = dtim_1234(4) - dtim_1234(3)
         if(iteration_electronic <= ISTAGE) then
            xiter = iteration_electronic
         else
            xiter = 2*ISTAGE - iteration_electronic
         end if
         if(xiter < 1.d0) xiter = 1.d0
         fex    = (xiter - ISTAGE*0.5)/ISTAGE
      end if
      fac    = dtanh(COEFF*fex)*0.5 + 0.5
      dtim   = dtims + fac * dtimxd
    end subroutine decide_dtim_case3
  end function m_CtrlP_dtim_now

! === Merge modifications from phase@63 -> phase@64. 2011/10/18 ================
! integer function m_CtrlP_On_or_Off_precon_WFs(iteration_electronic,iteration_ionic)
!   integer, intent(in) :: iteration_electronic,iteration_ionic
!   integer             :: i, imsd_t, ips, ipe
  integer function m_CtrlP_On_or_Off_precon_WFs(iteration_electronic,iteration_ionic,isolver_adopted)
    integer, intent(in) :: iteration_electronic,iteration_ionic,isolver_adopted
    integer             :: i, imsd_t, ips, ipe, isolver, j
! ==============================================================================
    if(tag_solver_of_WF_is_found) then
       if(iteration_ionic >= 2) then
          if(n_WF_solvers_after == 0) then
             ips = 1
          else
             ips = n_WF_solvers_before + 1
          end if
          ipe = n_WF_solvers_all
       else
          ips = 1
          ipe = n_WF_solvers_before
       end if
       if(ipe < ips) call phase_error_with_msg(6,' n_WF_solvers_(before|after) are illegal (m_CtrlP_solver_for_WFs_now)'&
                                              ,__LINE__,__FILE__)
       do i = ips, ipe
          if(w_solver(i)%till_n_iter <= 0) exit
          if(iteration_electronic <= w_solver(i)%till_n_iter) exit
       end do
       if(i > ipe) i = ipe
! === Merge modifications from phase@63 -> phase@64. 2011/10/18 ================
       isolver = w_solver(i)%solver
       if(isolver/= isolver_adopted) then
          do j = i-1, ips, -1
             if(w_solver(j)%solver == isolver_adopted) exit
          end do
          if(j >= ips) i=j
       end if
! ==============================================================================
       m_CtrlP_On_or_Off_precon_WFs = w_solver(i)%precon
    else
       do i = 1, n_WF_solvers
          if(iteration_electronic <= till_n_iteration(i)) exit
       end do
       if( i > n_WF_solvers) i = n_WF_solvers
       imsd_t = WF_solver(i)
! === Merge modifications from phase@63 -> phase@64. 2011/10/18 ================
       if(imsd_t /= isolver_adopted) then
          imsd_t = isolver_adopted
       end if
! ==============================================================================
       m_CtrlP_On_or_Off_precon_WFs = On_or_Off_precon_for_WFs_imsd(imsd_t)

    end if
  end function m_CtrlP_On_or_Off_precon_WFs

  integer function On_or_Off_precon_for_WFs_imsd(imsd)
    integer, intent(in) :: imsd

    solver : select case(imsd)
       case(SD,MSD,CG,RMM,RMM2,RMM2P,SDLM,MSDLM,eazyCG)
         On_or_Off_precon_for_WFs_imsd = OFF
       case(SDPRC,MSDPRC,CGPRC,RMMPRC,RMM2PRC,RMM2PPRC,SDLMPRC,MSDLMPRC,eazyCGPRC)
         On_or_Off_precon_for_WFs_imsd = ON
    end select solver
  end function On_or_Off_precon_for_WFs_imsd

  integer function m_CtrlP_solver_for_WFs_now(iteration_electronic,iteration_ionic &
       &                                     ,intzaj,edeltb_per_atom,sw_submat)
    integer, intent(in)       :: iteration_electronic,iteration_ionic,intzaj
    real(kind=DP), intent(in) :: edeltb_per_atom
    integer, intent(out)      :: sw_submat
    integer                   :: i, imsd_t, ips, ipe
    real(kind=DP)             :: curr_ed_rmm

    sw_submat = OFF
    sw_submat_is_on = .false.
    if(tag_solver_of_WF_is_found) then
       if(iteration_ionic >= 2) then
          if(n_WF_solvers_after == 0) then
             ips = 1
          else
             ips = n_WF_solvers_before + 1
          end if
          ipe = n_WF_solvers_all
       else
          ips = 1
          ipe = n_WF_solvers_before
       end if
       if(ipe < ips) call phase_error_with_msg(6,' n_WF_solvers_(before|after) are illegal (m_CtrlP_solver_for_WFs_now)'&
                                              ,__LINE__,__FILE__)
       do i = ips, ipe
          if(w_solver(i)%till_n_iter < 0) exit
          if(iteration_electronic <= w_solver(i)%till_n_iter) exit
       end do
       if(i > ipe) i = ipe
       imsd_t = w_solver(i)%solver
       sw_submat = w_solver(i)%subspace_rotation
    else
       m_CtrlP_solver_for_WFs_now = -9999
       do i = 1, n_WF_solvers
          if(iteration_electronic <= till_n_iteration(i)) exit
       end do
       if(n_WF_solvers <= 0) call phase_error_with_msg(6,' n_WF_solvers is illegal (m_CtrlP_solver_for_WFs_now)'&
                                                      ,__LINE__,__FILE__)
       if( i > n_WF_solvers) i = n_WF_solvers
       imsd_t = WF_solver(i)
       sw_submat = w_solver(i)%subspace_rotation
    end if
    if(sw_submat == ON) then
       sw_submat_is_on = .true.
       if(.not.m_CtrlP_ntcnvg_clear()) then
          if(submat_period >= 2 .and. iteration_electronic >= 2 .and. &
               & mod(iteration_electronic-1, submat_period) /= 0 ) sw_submat = OFF
       end if
    end if

    curr_ed_rmm = edelta_change_to_rmm
    if(iteration_ionic>1) curr_ed_rmm = edelta_change_to_rmm_md
    solver : select case(imsd_t)
       case (SD,SDPRC)
          m_CtrlP_solver_for_WFs_now = SD
       case (MSD,MSDPRC)
          m_CtrlP_solver_for_WFs_now = MSD
       case (SDLM, SDLMPRC)
          m_CtrlP_solver_for_WFs_now = lmSD
       case (CG, CGPRC)
          m_CtrlP_solver_for_WFs_now = lmCG
       case (RMM,RMMPRC)
          m_CtrlP_solver_for_WFs_now = RMM
!!$          if(dabs(edeltb_per_atom) < edelta_change_to_rmm ) then
          if(dabs(edeltb_per_atom) < curr_ed_rmm ) then
             m_CtrlP_solver_for_WFs_now = RMM
          else
             m_CtrlP_solver_for_WFs_now = previous_solver
             if(previous_solver /= RMM) then
                if(i>=ips+1) sw_submat = w_solver(i-1)%subspace_rotation
                if(ipri>=2) write(6,'(" solver=RMM->previous_solver, i = ",i4)') i
             end if
          end if
       case (RMM2, RMM2PRC)
          m_CtrlP_solver_for_WFs_now = RMM2
!!$          if(dabs(edeltb_per_atom) < edelta_change_to_rmm ) then
          if(dabs(edeltb_per_atom) < curr_ed_rmm ) then
             m_CtrlP_solver_for_WFs_now = RMM2
          else
             m_CtrlP_solver_for_WFs_now = previous_solver
             if(previous_solver /= RMM2) then
                if(i>=ips+1) sw_submat = w_solver(i-1)%subspace_rotation
             end if
          end if
       case (RMM2P, RMM2PPRC)
          m_CtrlP_solver_for_WFs_now = RMM2P
!!$          if(dabs(edeltb_per_atom) < edelta_change_to_rmm ) then
          if(dabs(edeltb_per_atom) < curr_ed_rmm ) then
             m_CtrlP_solver_for_WFs_now = RMM2P
          else
             m_CtrlP_solver_for_WFs_now = previous_solver
             if(previous_solver /= RMM2P) then
                if(i>=ips+1) sw_submat = w_solver(i-1)%subspace_rotation
             end if
          end if
       case (eazyCG, eazyCGPRC)
          m_CtrlP_solver_for_WFs_now = lmeazyCG
       case (MSDLM, MSDLMPRC)
          m_CtrlP_solver_for_WFs_now = MSDLM
!!$       case (CG,CGPRC, RMM,RMMPRC, RMM2,RMM2PRC, RMM2P,RMM2PPRC,eazyCG,eazyCGPRC)
!!$          if(iteration_electronic < i_2lm) then
!!$             m_CtrlP_solver_for_WFs_now = MSD
!!$          else if(iteration_electronic <i_sd2another) then
!!$             m_CtrlP_solver_for_WFs_now = lmMSD
!!$          else
!!$             if(imsd == CG .or. imsd == CGPRC) then
!!$                m_CtrlP_solver_for_WFs_now = lmCG
!!$             else if(imsd == RMM .or. imsd == RMMPRC) then
!!$                m_CtrlP_solver_for_WFs_now = RMM
!!$             else if(imsd == RMM2 .or. imsd == RMM2PRC) then
!!$                m_CtrlP_solver_for_WFs_now = RMM2
!!$             else if(imsd == RMM2P .or. imsd == RMM2PPRC) then
!!$                m_CtrlP_solver_for_WFs_now = RMM2P
!!$             else if(imsd == eazyCG .or. imsd == eazyCGPRC) then
!!$                m_CtrlP_solver_for_WFs_now = lmeazyCG
!!$             endif
!!$          endif
       case (SUBMAT)
          m_CtrlP_solver_for_WFs_now = SUBMAT
       case (MATRIXDIAGON)
          m_CtrlP_solver_for_WFs_now = MATRIXDIAGON
       case (DAVIDSON)
          m_CtrlP_solver_for_WFs_now = DAVIDSON
       case (MDDAVIDSON)
          m_CtrlP_solver_for_WFs_now = MDDAVIDSON
       case (MDKOSUGI)
          m_CtrlP_solver_for_WFs_now = MDKOSUGI
       case default
          if(printable) then
             write(6,'(" ! error in <<m_CtrlP_solver_for_WFs_now>>")')
             write(6,'(" imsd_t = ",i6)') imsd_t
          end if
    end select solver

    if(sw_submat == ON) then
       sw_submat_is_on = .true.
       if(.not.m_CtrlP_ntcnvg_clear()) then
          if(submat_period >= 2 .and. iteration_electronic >= 2 .and. &
               & mod(iteration_electronic-1, submat_period) /= 0 ) sw_submat = OFF
       end if
    end if

    if(intzaj == by_matrix_diagon) then
       if(iteration_electronic == 1 .and. iteration_ionic <= 1) then
          m_CtrlP_solver_for_WFs_now = MATRIXDIAGON
       else if(skip_alloc_phonon.and.iteration_electronic == 1) then
          m_CtrlP_solver_for_WFs_now = MATRIXDIAGON
       end if
    end if

    if(m_CtrlP_solver_for_WFs_now == 0) then
       if(printable) then
          write(6,'(" !!D m_CtrlP_solver_for_WFs_now is happened to be 0.")')
          write(6,'(" !!D m_CtrlP_solver_for_WFs_now is set to be MSD(=1)")')
       end if
       m_CtrlP_solver_for_WFs_now = MSD
    end if

    if(iprisolver >= 2) &
         & write(6,'(" ! previous_solver, m_CtrlP_solver_for_WFs_now = ",2i8)') &
         & previous_solver, m_CtrlP_solver_for_WFs_now
    previous_solver = m_CtrlP_solver_for_WFs_now
    if(m_CtrlP_solver_for_WFs_now <= -9999) then
       call phase_error_with_msg(6,' ! illegal m_CtrlP_solver_for_WFs_now',__LINE__,__FILE__)
    end if
    isolver_now = m_CtrlP_solver_for_WFs_now
  end function m_CtrlP_solver_for_WFs_now

  integer function m_CtrlP_get_isolver_now()
    m_CtrlP_get_isolver_now = isolver_now
  end function m_CtrlP_get_isolver_now

  subroutine m_CtrlP_set_way_ksample(i)
    integer, intent(in) :: i
    way_ksample = i
    if(way_ksample /= GAMMA .and. way_ksample /= FILE &
         & .and. way_ksample /= SKPS_DIRECT_IN  &
         & .and. way_ksample /= MONKHORST_PACK  &
         & .and. way_ksample /= MESH) then
       if(printable) write(6,'(" way_ksample (= ",i5," ) is not proper ")') way_ksample
       call phase_error_with_msg(6,' illegal way_ksample (m_CtrlP_set_way_ksample)',__LINE__,__FILE__)
    end if
  end subroutine m_CtrlP_set_way_ksample
#endif
! endif of "_EMPIRICAL_"

! ==================================== modified by K. Tagami ============= 11.0
!  subroutine m_CtrlP_set_nspin_and_af(nfout,imag)
!    integer,intent(in) :: nfout,imag
!    if(imag == PARA) then
!       nspin = 1
!       af    = 0
!    else if(imag == ANTIFERRO) then
!       nspin = 2
!       af    = 1
!    else if(imag == FERRO) then
!       nspin = 2
!       af    = 0
!!    end if
!    if(ipriinputfile >= 1) write(nfout,'(" !** imag, nspin, af = ",3i6," <<m_CtrlP_set_nspin_and_af>>")') imag,nspin,af
!  end subroutine m_CtrlP_set_nspin_and_af

  subroutine m_CtrlP_set_nspin_and_af( nfout,imag )
    integer,intent(in) :: nfout,imag

    select case (imag)
    case (PARA)
       nspin = 1;   af = 0;   ndim_spinor = 1;
    case (ANTIFERRO)
       nspin = 2;   af = 1;   ndim_spinor = 1;
    case (FERRO)
       nspin = 2;   af = 0;   ndim_spinor = 1
    case (NONCOLLINEAR)
       nspin = 2;   af = 0;   ndim_spinor = 2;     noncol = .true.
    end select
!
    if ( noncol ) then
      ndim_chgpot = ndim_spinor**2;       ndim_magmom = 4
    else
      ndim_chgpot = nspin;                ndim_magmom = nspin
    endif
!
    if (ipriinputfile >= 1) then
       write(nfout,'(" !** imag, nspin, af = ",3i6," <<m_CtrlP_set_nspin_and_af>>")') &
            &           imag,nspin,af
       write(nfout,*) '!** ndim_spinor = ', ndim_spinor
    endif
  end subroutine m_CtrlP_set_nspin_and_af
! ====================================================================== 11.0

  subroutine m_CtrlP_set_af(af_t)
    integer,intent(in) :: af_t
    af = af_t
    if(ipriinputfile >= 2 .and. printable) write(6,'(" -- af = ",i6, " <m_CtrlP_set_af>")') af
  end subroutine m_CtrlP_set_af

  integer function m_CtrlP_what_is_mdalg()
    m_CtrlP_what_is_mdalg = imdalg
  end function m_CtrlP_what_is_mdalg

  subroutine m_CtrlP_reset_optmode()
    optmode = -999
  end subroutine m_CtrlP_reset_optmode

  integer function m_CtrlP_set_gdiisoptmode(iteration_ionic,forcmx)
    integer, intent(in) ::      iteration_ionic
    real(kind=DP),intent(in) :: forcmx
    if(optmode < -100) optmode = initial_method_of_gdiis ! (QUENCHED_MD)
    m_CtrlP_set_gdiisoptmode = optmode
!!$    if(optmode == QUENCHED_MD .and. forcmx <= c_forc2GDIIS) then
    if(optmode == initial_method_of_gdiis .and. forcmx <= c_forc2GDIIS) then
       if(iteration_ionic >= c_iteration2GDIIS) then
          optmode = GDIIS
          m_CtrlP_set_gdiisoptmode = GDIIS
       else
          optmode = initial_method_of_gdiis ! (QUENCHED_MD)
          m_CtrlP_set_gdiisoptmode = optmode
       end if
    end if
!!$    m_CtrlP_set_gdiisoptmode = optmode
  end function m_CtrlP_set_gdiisoptmode

  subroutine m_CtrlP_set_iconvergence(i)
    integer, intent(in) :: i
    iconvergence = i
  end subroutine m_CtrlP_set_iconvergence

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_set_xctype(type,nfout)
    character(len=*), intent(in) :: type
    integer, intent(in) :: nfout
    integer :: ndif, i, icharx
    if(len_trim(type) > len_xctype ) then
       xctype = type(1:len_xctype)
    else
       xctype = type
    end if

    ndif = ichar('A') - ichar('a')
    do i = 1, len_trim(xctype)
       icharx = ichar(xctype(i:i))
       if(icharx >= ichar('A') .and. icharx <= ichar('Z')) then
          icharx = icharx - ndif
          xctype(i:i) = char(icharx)
       end if
    end do
    if(printable) write(nfout,'(" !CtrlP -- xctype is set to be ",a7)') xctype
  end subroutine m_CtrlP_set_xctype

  subroutine m_CtrlP_set_submat(submat_is_done)
    logical, intent(in) :: submat_is_done
    if(npes > 1) then
       call mpi_allreduce(submat_is_done,submat_is_done_this_iter &
                       & ,1,mpi_logical,mpi_lor,MPI_CommGroup,ierr)
    else
       submat_is_done_this_iter = submat_is_done
    end if
  end subroutine m_CtrlP_set_submat

  subroutine m_Ctrlp_set_renew_wf(tf)
    logical, intent(in) :: tf
    renew_wf_again_m_CtrlP = tf
  end subroutine m_Ctrlp_set_renew_wf

  subroutine m_CtrlP_check_naldos_range(nfout,maldos)
    integer, intent(in) :: nfout,maldos

    if(naldos_from == 0 .and. naldos_to == 0) then
       naldos_from = 1
       naldos_to   = maldos
    else
       if(naldos_from < 1)     naldos_from = 1
       if(naldos_from > maldos) naldos_from = maldos
       if(naldos_to < 1)       naldos_to   = 1
       if(naldos_to > maldos)   naldos_to   = maldos
       if(naldos_to < naldos_from) naldos_to = naldos_from
    end if
    if(printable) then
       write(nfout,'(" !!ctrlP naldos_from         = ",i6," <<m_CtrlP_check_naldos_range>>")') naldos_from
       write(nfout,'(" !!ctrlP naldos_to           = ",i6," <<m_CtrlP_check_naldos_range>>")') naldos_to
    end if
  end subroutine m_CtrlP_check_naldos_range

  subroutine set_polar_prop(rstr,polar_prop_t)
    character(len=FMAXVALLEN),intent(in) :: rstr
    integer, intent(out) :: polar_prop_t
    logical :: tf
! debug
!    if(printable) write(*,*) 'debug: polar prop=',rstr
! end debug
    call strncmp0(tag_polarization, trim(rstr), tf)
    if(tf) polar_prop_t = POLARIZATION
    call strncmp0(tag_effective_charge, trim(rstr), tf)
    if(tf) polar_prop_t = EFFECTIVE_CHARGE
    call strncmp0(tag_piezoelectric_const, trim(rstr), tf)
    if(tf) polar_prop_t = PIEZOELECTRIC_CONST
  end subroutine set_polar_prop

  subroutine m_CntrlP_keep_charge_title()
    charge_title_tmp = charge_title
  end subroutine m_CntrlP_keep_charge_title

  subroutine m_CntrlP_retrieve_charge_title()
    charge_title = charge_title_tmp
  end subroutine m_CntrlP_retrieve_charge_title

  subroutine m_CntrlP_set_pcharge_title(nspin,iloop,i,emin,emax)
    integer, intent(in) :: nspin, iloop, i
    real(kind=DP), intent(in) :: emin, emax
    charge_title = ''
    if(nspin == 1) then
       write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," )")') i,emin,emax
    else
       if(iloop == 1) then
          write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) spin = UP")') i,emin,emax
       else if ( iloop == 2 ) then
          write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) spin = DOWN")') i,emin,emax
       end if
! ======= KT_add ====== 2014/06/07
       if ( iloop == -1 ) then
          write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) TOTAL")') &
               &               i, emin, emax
       else if ( iloop == -2 ) then
          write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) MAGMOM")')&
               &               i, emin, emax
       end if
! ===================== 2014/06/07
    end if
  end subroutine m_CntrlP_set_pcharge_title

! ======================  added by K. Tagami ====================== 11.0
  subroutine m_CtrlP_set_pchg_title_noncl(iloop,i,emin,emax)
    integer, intent(in) :: iloop, i
    real(kind=DP), intent(in) :: emin, emax
    charge_title = ''

    select case (iloop)
    case (1)
       write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) tot")') &
                          &  i,emin,emax
    case (2)
       write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) mx")') &
                          &  i,emin,emax
    case (3)
       write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) my")') &
                          &  i,emin,emax
    case (4)
       write(charge_title,'(" nEwindow = ",i5," : ( ",f9.4,",",f9.4," ) mz")') &
                          &  i,emin,emax
    end select
  end subroutine m_CtrlP_set_pchg_title_noncl
! ======================================================================== 11.0

#endif

  integer function m_CtrlP_flag_wd_force()
    m_CtrlP_flag_wd_force = flag_wd_force
  end function m_CtrlP_flag_wd_force

  subroutine m_CtrlP_set_flag_wd_force_1()
    flag_wd_force = 1
  end subroutine m_CtrlP_set_flag_wd_force_1

#ifndef _EMPIRICAL_
  subroutine m_CtrlP_set_projectors(prealloc,m)
    logical, intent(in) :: prealloc
    integer, intent(inout) :: m
! === Include PAW by tkato =====================================================
!   integer :: i, iret, nsize
    integer :: i, iret
! ==============================================================================
    integer :: f_selectBlock, f_selectParentBlock
! === Include PAW by tkato =====================================================
!   integer :: f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_getRealValue, f_getIntValue
! ==============================================================================
    integer :: f_selectFirstTableLine, f_selectNextTableLine
    integer :: no, ig
    real(kind=DP) dret
    integer, allocatable :: igroup(:)

    if(.not.prealloc) allocate(igroup(m))

    if( f_selectBlock(tag_projectors) == 0) then
       i = 1
       do while(.true.)
          if(i==1) then
             if( f_selectFirstTableLine() /= 0) then
                exit
             end if
          else
             if( f_selectNextTableLine() /= 0) then
                exit
             end if
          end if
          if(.not.prealloc) then
             if(i > m) exit
             no = i
             if( f_getIntValue(tag_no,iret) == 0) no=iret
	     proj_attribute(no)%radius=1.0 ! default
             if( f_getRealValue(tag_radius,dret,"bohr") == 0) then
                 proj_attribute(no)%radius=dret
                 proj_attribute(no)%radius_was_defined = .true.
             endif
	     proj_attribute(no)%fwhm=0.1 ! default
             if( f_getRealValue(tag_fwhm,dret,"bohr") == 0) proj_attribute(no)%fwhm=dret
	     proj_attribute(no)%l=0 ! default
             if( f_getIntValue(tag_l,iret) == 0) proj_attribute(no)%l=iret
	     proj_attribute(no)%t=1 ! default
             if( f_getIntValue(tag_t,iret) == 0) proj_attribute(no)%t=iret
	     proj_attribute(no)%phi=0 ! default
             proj_attribute(no)%frotate=.false. ! default
             if( f_getRealValue(tag_phi,dret,"radian") == 0) proj_attribute(no)%phi=dret
             if(abs(dret) > 1.d-10) proj_attribute(no)%frotate=.true.
	     proj_attribute(no)%theta=0 ! default
             if( f_getRealValue(tag_theta,dret,"radian") == 0) proj_attribute(no)%theta=dret
             if(abs(dret) > 1.d-10) proj_attribute(no)%frotate=.true.
	     proj_attribute(no)%psi=0 ! default
             if( f_getRealValue(tag_psi,dret,"radian") == 0) proj_attribute(no)%psi=dret
             if(abs(dret) > 1.d-10) proj_attribute(no)%frotate=.true.
	     proj_attribute(no)%group=0 ! default
             if( f_getIntValue(tag_group,iret) == 0) igroup(no) = iret
             proj_attribute(no)%group=igroup(no)
             proj_attribute(no)%strong_correlated=.false.
             proj_attribute(no)%Ueff=0.d0
             proj_attribute(no)%initialUeff = 0.d0
             proj_attribute(no)%finalUeff = 0.d0
             proj_attribute(no)%deltaUeff = 0.d0

             proj_attribute(no)%initial_spin_sum = 0.d0
             proj_attribute(no)%initial_spin_diff = 0.d0
          end if
          i = i + 1
       end do
       iret = f_selectParentBlock()
    else
       call phase_error_with_msg(6,' ! No projector is given in the inputfile <<m_CtrlP_set_projectors>>', &
       __LINE__,__FILE__)
    end if
    if(prealloc) m = i -1
    if(.not.prealloc) then
       num_proj_group = 0
       do i=1,num_projectors
          if(igroup(i)>num_proj_group) then
             num_proj_group = igroup(i)
          end if
       end do

! === KT_mod ===== 2014/06/05
!       allocate(proj_group(max_projs,num_proj_group))
!
       if ( sw_allow_maxprojs_gt_4 == ON ) then
          allocate( num_projs_in_each_group( num_proj_group ) );
          num_projs_in_each_group = 0
          Do i=1, num_projectors
             ig = igroup(i)
             num_projs_in_each_group(ig) = num_projs_in_each_group(ig) +1
          End do

          max_projs = maxval( num_projs_in_each_group )
          deallocate( num_projs_in_each_group )
       else
          max_projs = max_projs_org
       endif

       allocate( proj_group( max_projs, num_proj_group ) )
! =============== 2014/06/05

       allocate(num_proj_elems(num_proj_group)); num_proj_elems=0
       do i=1,num_projectors
          ig = igroup(i)
          num_proj_elems(ig) = num_proj_elems(ig)+1
          proj_group(num_proj_elems(ig),ig) = i
          proj_attribute(i)%ielem=num_proj_elems(ig)
       end do
       deallocate(igroup)
    end if
  end subroutine m_CtrlP_set_projectors

  subroutine m_CtrlP_set_hubbard_proj(nfout)
    integer, intent(in) :: nfout
    integer :: i, iret, nsize
    integer :: f_selectBlock, f_selectParentBlock
    integer :: f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectFirstTableLine, f_selectNextTableLine
    integer :: no
    real(kind=DP) dret

    if( f_selectBlock(tag_projectors) == 0) then
       i = 1
       do while(.true.)
          if(i==1) then
             if( f_selectFirstTableLine() /= 0) then
                exit
             end if
          else
             if( f_selectNextTableLine() /= 0) then
                exit
             end if
          end if
          if( f_getIntValue(tag_no,iret) == 0) no=iret
          if( no > num_projectors) cycle
          proj_attribute(no)%strong_correlated=.true.
          nsize = 2*proj_attribute(no)%l+1
          allocate(proj_attribute(no)%component(nsize))
          proj_attribute(no)%component = .true.
          if( f_getStringValue(tag_component,rstr,LOWER) == 0) then
             call set_logical_value(rstr,nsize,proj_attribute(no)%component)
          end if
          if( f_getRealValue(tag_Ueff,dret,"hartree") == 0) proj_attribute(no)%Ueff=dret
          proj_attribute(no)%norbital=nsize
          if( f_getIntValue(tag_norbital,iret) == 0) proj_attribute(no)%norbital=iret
          if( f_getRealValue(tag_initialUeff,dret,"hartree") == 0) proj_attribute(no)%initialUeff = dret
          if( f_getRealValue(tag_finalUeff,dret,"hartree") == 0) proj_attribute(no)%finalUeff = dret
          if(nUeff>0 .and. driver == DRIVER_URAMP)then
             proj_attribute(no)%deltaUeff = &
          & (proj_attribute(no)%finalUeff - proj_attribute(no)%initialUeff)/real(nUeff-1)
             proj_attribute(no)%Ueff = proj_attribute(no)%initialUeff
             if(printable) write(nfout,'(a,i5,3f10.5)') ' !** no initialUeff finalUeff deltaUeff', &
             & no,proj_attribute(no)%initialUeff,proj_attribute(no)%finalUeff, proj_attribute(no)%deltaUeff
          endif
          if( f_getRealValue(tag_initial_spin_sum,dret,"") == 0) &
               &           proj_attribute(no)%initial_spin_sum = dret
          if( f_getRealValue(tag_initial_spin_diff,dret,"") == 0) &
               &           proj_attribute(no)%initial_spin_diff = dret
          i = i + 1
       end do
       iret = f_selectParentBlock()
    else
       call phase_error_with_msg(nfout,' ! No projector is given in the inputfile <<m_CtrlP_set_hubbard_proj>>'&
                                ,__LINE__,__FILE__)
    end if
  contains
    subroutine set_logical_value(rstr,nsize,component)
      character(len=FMAXVALLEN), intent(in) :: rstr
      integer, intent(in) :: nsize
      logical, intent(out) :: component(nsize)

      integer :: i,ilen,ival,ist
      character(len=FMAXVALLEN) :: buf

      ilen = len_trim(rstr)
      if(ilen == 0) return

      component = .false.
      ist = 1
      do i=1,ilen
         if(rstr(i:i) == ":".or.i==ilen) then
            if(i/=ilen) then
               buf = rstr(ist:i-1)
            else
               buf = rstr(ist:ilen)
            end if
            read(buf,*) ival
            if(ival <= nsize .and. ival >=1) then
               component(ival) = .true.
            else
               call phase_error_with_msg(nfout,"bad componet value <<m_CtrlP_set_hubbard_proj.set_logical_value>>"&
                                        ,__LINE__,__FILE__)
            end if
            ist = i+1
         end if
      end do
    end subroutine set_logical_value
  end subroutine m_CtrlP_set_hubbard_proj

! =========================== added by K. Tagami ===================  11.0
  subroutine m_CtrlP_set_spinorbit_proj()
    integer :: i, iret, nsize
    integer :: f_selectBlock, f_selectParentBlock
    integer :: f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectFirstTableLine, f_selectNextTableLine
    integer :: no

    real(kind=DP) dret, my_factor, ctmp1

    if( f_selectBlock(tag_projectors) == 0) then
       i = 1
       do while(.true.)
          if(i==1) then
             if( f_selectFirstTableLine() /= 0) then
                exit
             end if
          else
             if( f_selectNextTableLine() /= 0) then
                exit
             end if
          end if
          if( f_getIntValue(tag_no,iret) == 0) no=iret

          if( no > num_projectors ) cycle
!!!!!          proj_attribute(no)%activate_soc=.true.
!!!!!          write(*,*) '** activate_soc     = ',  no, proj_attribute(no)%activate_soc

          nsize = 2*proj_attribute(no)%l+1

          allocate(proj_attribute(no)%component(nsize))
          proj_attribute(no)%component = .true.
          if( f_getStringValue(tag_component,rstr,LOWER) == 0) then
            call set_logical_value(rstr,nsize,proj_attribute(no)%component)
          end if

          if( f_getRealValue( tag_scaling_factor, dret,"") == 0) then
             proj_attribute(no)%LScoupling_scaling_factor = dret
             write(*,*) '** scaling factor       = ', dret
          endif
!
          if( f_getRealValue(tag_LScoupling, dret,"hartree") == 0) then
             proj_attribute(no)%LScoupling0 = dret
             write(*,*) '** LScoupling      = ', dret
          endif
          if( f_getRealValue( tag_Splitting, dret,"hartree") == 0) then
             my_factor = nsize *0.5d0          ! 1/2 *( 2l + 1 )
!
             proj_attribute(no)%LScoupling0 = dret / my_factor
             write(*,*) '** splitting width = ', dret
             write(*,*) '** LScoupling      = ', dret /my_factor
          endif
!
          proj_attribute(no)%norbital=nsize
          if( f_getIntValue(tag_norbital,iret) == 0) proj_attribute(no)%norbital=iret
          i = i + 1
       end do
       iret = f_selectParentBlock()
    else
       call phase_error_with_msg(6,' ! No projector is given in the inputfile <<m_CtrlP_set_hubbard_proj>>'&
                                ,__LINE__,__FILE__)
    end if
  contains
    subroutine set_logical_value(rstr,nsize,component)
      character(len=FMAXVALLEN), intent(in) :: rstr
      integer, intent(in) :: nsize
      logical, intent(out) :: component(nsize)

      integer :: i,ilen,ival,ist
      character(len=FMAXVALLEN) :: buf

      ilen = len_trim(rstr)
      if(ilen == 0) return

      component = .false.
      ist = 1
      do i=1,ilen
         if(rstr(i:i) == ":".or.i==ilen) then
            if(i/=ilen) then
               buf = rstr(ist:i-1)
            else
               buf = rstr(ist:ilen)
            end if
            read(buf,*) ival
            if(ival <= nsize .and. ival >=1) then
               component(ival) = .true.
            else
               call phase_error_with_msg(6,"bad componet value <<m_CtrlP_set_spinorbit_proj.set_logical_value>>"&
                                        ,__LINE__,__FILE__)
            end if
            ist = i+1
         end if
      end do
    end subroutine set_logical_value
  end subroutine m_CtrlP_set_spinorbit_proj
! ================================================================== 11.0

  subroutine m_CtrlP_wd_proj_attr(nfout)
    integer, intent(in) :: nfout
    integer :: i
    type(t_projector), pointer :: p(:)

    p => proj_attribute

    write(nfout,'(" !** Projector attributes **")')
    write(nfout,'(" !** no",2x,"radius",3x,"fwhm",5x,"Ueff",5x,"phi",6x,"theta",4x,"psi",5x,"l",1x,"component")')
    do i=1,num_projectors
       if(p(i)%strong_correlated) then
          write(nfout,'(" !**",i3,6f9.5,i2,7l2)') i &
          & ,p(i)%radius, p(i)%fwhm, p(i)%Ueff &
          & ,p(i)%phi,p(i)%theta,p(i)%psi &
          & ,p(i)%l, p(i)%component
       else
          write(nfout,'(" !**",i3,6f9.5,i2)') i &
          & ,p(i)%radius, p(i)%fwhm, p(i)%Ueff &
          & ,p(i)%phi,p(i)%theta,p(i)%psi &
          & ,p(i)%l
       end if
    end do

  end subroutine m_CtrlP_wd_proj_attr

  subroutine m_CtrlP_set_proj_ityp(ig,ityp)
    integer, intent(in) :: ig, ityp
    integer :: i,no
    if(ig<1) return
    if(.not.allocated(num_proj_elems)) call phase_error_with_msg(6,'num_proj_elems is not allocated.', &
    __LINE__,__FILE__)
    do i=1,num_proj_elems(ig)
       no = proj_group(i,ig)
       proj_attribute(no)%ityp = ityp
    end do
  end subroutine m_CtrlP_set_proj_ityp

  subroutine m_CntrlP_set_neg(i,nfout)
    integer, intent(in) :: i,nfout
    neg = i
    max_subspace_size = 4*neg
    write(nfout,'(" !** REMARK: the number of bands has been redefined")')
    if( sw_divide_subspace_changed )then
       if(neg/nrank_e>=4) then
          sw_divide_subspace = ON
          if(printable) then
             write(nfout,'(" !** REMARK: sw_divide_subspace is reset to ON ")')
          endif
       endif
    endif
    if(sw_npartition_changed) call m_CtrlP_set_npartition_david(nfout)
    if(icond==INITIAL) neg_previous = neg
  end subroutine m_CntrlP_set_neg

  subroutine m_CtrlP_set_npartition_david(nfout)
     integer, intent(in) :: nfout
     if(sw_divide_subspace==ON)then
       npartition_david = neg/(nblock*nrank_e)
       if (npartition_david<1) then
          npartition_david = 1
       endif
       if(printable) write(nfout,'(a,i8)') " !** REMARK: npartition_david is reset to : ",npartition_david
     endif
  end subroutine m_CtrlP_set_npartition_david

  subroutine m_CntrlP_rst_submat_call_stat()
     submat_uncalled = .true.
  end subroutine m_CntrlP_rst_submat_call_stat

  subroutine m_CntrlP_set_meg(i)
    integer, intent(in) :: i
    meg = i
  end subroutine m_CntrlP_set_meg

  subroutine m_Cntrlp_set_davidson_size(i)
    integer, intent(in) :: i
    max_subspace_size = i
  end subroutine m_Cntrlp_set_davidson_size

  subroutine m_CtrlP_set_neg_properly(totch)
    real(kind=DP), intent(in) :: totch
    real(kind=DP), parameter :: p = -0.4d0
    real(kind=DP) :: t
!!$    neg = int((ceiling(totch/2)+1)*(1+0.05/dlog10(max(totch,2.d0))))
!!$    neg = int((ceiling(totch/2)+1)*(1+0.10/dlog10(max(totch,2.d0))))
!!$    if(neg*2 < totch + 4.0) neg = ceiling(totch/2) + 2
    t = totch/2.d0
    if(t<4) then
       neg = 8
    else if(t <= 10) then
       neg = int(t+4.d0)
    else
       neg = int(t*(1+t**p))
    end if
    if(driver == DRIVER_NEB .or. driver == DRIVER_DIMER) neg_previous = neg
    neg_is_given = .true.
    if(neg_previous <= 1) neg_previous = neg
    if(ekmode == ON .or. ekmode == GRID .or. &
    & icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION .or. icond==FIXED_CHARGE_AUTOMATIC) then
      if(.not.num_extra_bands_was_set) call m_CtrlP_set_def_numextrabands(neg)
       neg = neg + num_extra_bands
       neg_is_enlarged = .true.
    endif
  end subroutine m_CtrlP_set_neg_properly

  subroutine m_CtrlP_set_def_numextrabands(neg)
    integer, intent(in) :: neg
    integer :: ntot
    real(kind=DP), parameter :: p = -0.4d0
!    ntot = neg+num_extra_bands
    ntot = neg
    if(ntot<4) then
       ntot = 8
    else if(ntot <= 10) then
       ntot = ntot+4
    else
       ntot = int(ntot*(1+ntot**p))
    end if
    num_extra_bands = ntot-neg
    if(printable) write(6,'(a,i8)') ' !** default value for num_extra_bands ',num_extra_bands
  end subroutine m_CtrlP_set_def_numextrabands

  subroutine m_CtrlP_flag_mpi_G_dot_R(nfout,nbmx)
    integer,intent(in) :: nfout,nbmx
    integer :: i
    flag_mpi_g_dot_r = .false.
    if(npes >= ngnode_nbmx*2) then
       i = ceiling(dble(npes)/ngnode_nbmx)
       if(nbmx/i >= ncritical_vectorlength_nbmx) flag_mpi_g_dot_r = .true.
    end if

    flag_mpi_g_dot_r_k = .false.
    if(nrank_e >= ngnode_nbmx*2) then
       i = ceiling(dble(nrank_e)/ngnode_nbmx)
       if(nbmx/i >= ncritical_vectorlength_nbmx) flag_mpi_g_dot_r_k = .true.
    end if

    if(printable) then
       write(nfout,'(" !CtrlP  nbmx,  ngnode_nbmx          = ",2i12)') nbmx, ngnode_nbmx
       write(nfout,'(" !CtrlP  ncritical_vectorlength_nbmx = ",i12)') ncritical_vectorlength_nbmx
       if(flag_mpi_g_dot_r) then
          write(nfout,'(" !CtrlP flag_mpi_g_dot_r = .true.")')
       else
          write(nfout,'(" !CtrlP flag_mpi_g_dot_r = .false.")')
       end if

       write(nfout,'(" !CtrlP  nrank_e                     = ",2i12)') nrank_e
       if(flag_mpi_g_dot_r_k) then
          write(nfout,'(" !CtrlP flag_mpi_g_dot_r_k = .true.")')
       else
          write(nfout,'(" !CtrlP flag_mpi_g_dot_r_k = .false.")')
       end if
    end if
  end subroutine m_CtrlP_flag_mpi_G_dot_R

  subroutine set_method_scalapack(method)
    integer, intent(inout) :: method

    integer :: f_getStringValue
    logical :: tf
    if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
       call strncmp0(trim(rstr),tag_householder,tf)
       if(tf) then
          method = HOUSEHOLDER
          goto 1001
       end if
       call strncmp0(trim(rstr),tag_divide_and_conquer,tf)
       if(tf) then
          method = DIVIDEandCONQUER
          goto 1001
       end if
1001   continue
    end if
  end subroutine set_method_scalapack

  subroutine set_charge_filetype(rstr,filetype)
    ! Revised by T. Yamasaki, 28 July 2008
    character(len=FMAXVALLEN),intent(in) :: rstr
    integer, intent(out) :: filetype
    logical :: tf
    call strncmp2(rstr, FMAXVALLEN, tag_cube, len(tag_cube),tf)
    if(tf) then
       filetype = CUBE
       goto 1001
    end if
    call strncmp2(rstr, FMAXVALLEN, tag_vtk, len(tag_vtk),tf)
    if(tf) then
       filetype = VTK
       goto 1001
    end if
    call strncmp2(rstr, FMAXVALLEN, tag_density_only, len(tag_density_only),tf)
    if(tf) then
       filetype = DENSITY_ONLY
       goto 1001
    end if
1001 continue
  end subroutine set_charge_filetype

  integer function m_CtrlP_cachesize()
    integer :: ncache
    if(ncachesize_given >= 0) then
       ncache = ncachesize_given
    else
#ifdef VPP
       ncache = 0
#elif SX
       ncache = 0
#elif SX4
       ncache = 0
#elif HIUX
       ncache = 0
#elif _CACHESIZE_FUNCTION_
       ncache = cachesize(3)
#else
       !!$ncache = 1024
       ncache = 256
#endif
    end if
    m_CtrlP_cachesize = ncache
  end function m_CtrlP_cachesize
#endif

  subroutine m_CtrlP_rd_driver(nfout)
    use m_Const_Parameters, only : DRIVER_GENERAL, DRIVER_CONSTRAINT, DRIVER_NEB,  LOWER, ON &
            , DRIVER_DIMER
    integer, intent(in) :: nfout
    integer :: f_selectTop, f_selectBlock, f_getStringValue, f_getIntValue
    integer :: iiret
    character(len=256) :: cret
    driver = DRIVER_GENERAL
    iiret = f_selectTop()
    iiret = f_selectBlock('control')
    if(f_getIntValue(tag_multiple_replica_mode,iiret)==0)then
      if(iiret==ON) driver = DRIVER_NEB
      if(iiret == ON) then
        if(f_getStringValue(tag_multiple_replica_method, cret, LOWER) == 0) then
          if(adjustl(trim(cret)) .eq. 'dimer') driver = DRIVER_DIMER
        endif
      endif
    endif

    if(f_getStringValue(tag_driver,cret,LOWER)==0)then
       if (adjustl(trim(cret)).eq.'constraints')then
         driver = DRIVER_CONSTRAINT
       else if (adjustl(trim(cret)).eq.'meta_dynamics')then
         driver = DRIVER_MTD
       else if (adjustl(trim(cret)).eq.'general')then
         driver = DRIVER_GENERAL
       else if (adjustl(trim(cret)).eq.'neb')then
         driver = DRIVER_NEB
       else if (adjustl(trim(cret)).eq.'dimer')then
         driver = DRIVER_DIMER
       else if (adjustl(trim(cret)).eq.'uramp')then
         driver = DRIVER_URAMP
       else if (adjustl(trim(cret)).eq.tag_sc_dft)then
         driver = DRIVER_SC_DFT
       else
         if(printable) write(nfout,'(a)') 'unknown driver : '//trim(adjustl(cret))
         call phase_error_with_msg(nfout,'unknown driver : '//trim(adjustl(cret)),__LINE__,__FILE__)
       endif
    endif
    iiret = f_selectTop()
    if(ipriinputfile >= 1) then
    if(driver==DRIVER_GENERAL.and.printable)    write(nfout,'(a)') ' !** driver = general'
    if(driver==DRIVER_CONSTRAINT.and.printable) write(nfout,'(a)') ' !** driver = constraints'
    if(driver==DRIVER_MTD.and.printable)        write(nfout,'(a)') ' !** driver = meta_dynamics'
    if(driver==DRIVER_NEB.and.printable)        write(nfout,'(a)') ' !** driver = NEB'
    if(driver==DRIVER_DIMER.and.printable)      write(nfout,'(a)') ' !** driver = DIMER'
    if(driver==DRIVER_URAMP.and.printable)      write(nfout,'(a)') ' !** driver = U_ramping'
    if(driver==DRIVER_SC_DFT.and.printable)     write(nfout,'(a)') ' !** driver = '//trim(tag_sc_dft)
    endif
  end subroutine m_CtrlP_rd_driver

  subroutine m_CtrlP_push_SolverNameApplied(solvername,lensolver)
    integer, intent(in) ::  lensolver
    character(len=lensolver),intent(in) :: solvername
    number_of_solvers_applied = number_of_solvers_applied + 1
    if(number_of_solvers_applied <= m_solvers_applied) then
       solver_names_applied(number_of_solvers_applied) = solvername(1:min(lensolver,len_solvername))
    end if
  end subroutine m_CtrlP_push_SolverNameApplied

  subroutine m_CtrlP_clear_nsolver_applied()
    number_of_solvers_applied = 0
  end subroutine m_CtrlP_clear_nsolver_applied

  subroutine m_CtrlP_push_CDMixingNameApplied(cdmixingname,lencdmixingname)
    integer, intent(in) ::  lencdmixingname
    character(len=lencdmixingname),intent(in) :: cdmixingname
    number_of_cdmixing_applied = number_of_cdmixing_applied + 1
    if(number_of_cdmixing_applied <= m_cdmixing_applied) then
       cdmixing_names_applied(number_of_cdmixing_applied) = cdmixingname(1:min(lencdmixingname,len_cdmixingname))
    end if
  end subroutine m_CtrlP_push_CDMixingNameApplied

  subroutine m_CtrlP_clear_cdmixing_applied()
    number_of_cdmixing_applied = 0
  end subroutine m_CtrlP_clear_cdmixing_applied

! === For restart lm+MSD! by tkato 2012/02/15 ==================================
  subroutine m_CtrlP_wd_dtim_previous(nfcntn)
    integer, intent(in) :: nfcntn
    if(mype == 0) then
       write(nfcntn, '(a13)') tag_dtim_previous
       write(nfcntn, '(d24.16)') dtim_1Dsearch
    end if
  end subroutine m_CtrlP_wd_dtim_previous

  subroutine m_CtrlP_rd_dtim_previous(nfcntn)
    integer, intent(in) :: nfcntn
    logical :: EOF_reach, tag_is_found
    if(mype==0) then
       call rewind_to_tag0(nfcntn, len(tag_dtim_previous), tag_dtim_previous, &
                           EOF_reach, tag_is_found, str, len_str)
       if(.not. tag_is_found) then
          write(0, *) 'NOTE: Restart data of dtim_1Dsearch cannot be read!'
          write(0, *) 'NOTE: (Restart file is of old version?)'
          write(0, *) 'NOTE: dtim_1Dsearch is set as -1 and restart with'
          write(0, *) 'NOTE: lm+MSD should be NOT GOOD!!!'
          dtim_1Dsearch = -1.0d0
       else
          read(nfcntn, *) dtim_1Dsearch
       end if
    end if
    if(npes > 1) call mpi_bcast(dtim_1Dsearch, 1, mpi_double_precision, 0, &
                                MPI_CommGroup, ierr)
    if(printable) write(6, '(d24.16, " : dtim_previous")') dtim_1Dsearch
  end subroutine m_CtrlP_rd_dtim_previous
! ==============================================================================

  subroutine m_CtrlP_set_init_status(logi)
    logical, intent(in) :: logi
    in_initialization = logi
  end subroutine m_CtrlP_set_init_status

  logical function m_CtrlP_in_initialization()
     m_CtrlP_in_initialization = in_initialization
  end function m_CtrlP_in_initialization

! =============================== KT_Test =========================== 12.5Exp
  subroutine m_CtrlP_set_hybrid_parameters
    if (npes>1) then
       call mpi_bcast( nmax_G_hyb, 1, mpi_integer, 0, MPI_CommGroup, ierr )

       call mpi_bcast( truncate_vxw_updating, 1, mpi_logical, 0, &
            &          MPI_CommGroup, ierr )

       call mpi_bcast( edelta_for_hyb_chgfix, 1, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
       call mpi_bcast( edelta_for_hyb_convgd, 1, mpi_double_precision, 0, &
            &          MPI_CommGroup, ierr )
    endif

  end subroutine m_CtrlP_set_hybrid_parameters
! =================================================================== 12.5Exp

  subroutine m_CtrlP_set_ppprinted(ivalue)
    integer, intent(in) :: ivalue
    if(ivalue >= 0) then
       ppprinted = .true.
    end if
  end subroutine m_CtrlP_set_ppprinted

  subroutine check_xctype_2016( nfout, rstr, xctype, exchange_pot_type, vdwdf_version )
    integer, intent(in) :: nfout

    character(len=FMAXVALLEN), intent(in) :: rstr
    character(len=len_xctype), intent(out) :: xctype
    character(len=len_xctype), intent(out) :: exchange_pot_type
    integer, intent(out) :: vdwdf_version

    if ( rstr == "vdwdf" ) then
       xctype = "vdwdf";   exchange_pot_type = "revpbe";   vdwdf_version = 1
    else if ( rstr == "vdwdf2" ) then
       xctype = "vdwdf";   exchange_pot_type = "pw86r";   vdwdf_version = 2
    else if ( rstr == "vdwdf-c09x" ) then
       xctype = "vdwdf";   exchange_pot_type = "c09x";    vdwdf_version = 1
    else if ( rstr == "vdwdf2-c09x" ) then
       xctype = "vdwdf";   exchange_pot_type = "c09x";    vdwdf_version = 2
    else if ( rstr == "vdwdf-optpbe" ) then
       xctype = "vdwdf";   exchange_pot_type = "optpbe";    vdwdf_version = 1
    else if ( rstr == "vdwdf-optb86b" ) then
       xctype = "vdwdf";   exchange_pot_type = "optb86b";    vdwdf_version = 1
    else if ( rstr == "vdwdf2-b86r" ) then
       xctype = "vdwdf";   exchange_pot_type = "b86r";    vdwdf_version = 2
    else if ( rstr == "vdwdf-cx" ) then
       xctype = "vdwdf";   exchange_pot_type = "lvpw86r";    vdwdf_version = 1
    else
       xctype = rstr(1:len_xctype)
       return
    endif
!
    write(nfout,*) ' [ vdwdf info @check_xctype_2016 ]'
    write(nfout,*) "!** exchange potential type = ", exchange_pot_type
    write(nfout,*) "!** vdwdf version no. = ", vdwdf_version

  end subroutine check_xctype_2016

! ===== KT_add ===== 13.0XX
  subroutine m_CtrlP_chkif_metagga( nfout )
    integer, intent(in) :: nfout

    if ( xctype == "tb09" ) then
       use_metagga = .true.
       sw_calc_ekin_density = ON

       write(nfout,*) " !** use_metagga is ", use_metagga
       write(nfout,'(A,I3)') "!** sw_calc_ekin_density is ", sw_calc_ekin_density

       if ( .not. use_modeled_ekin_density  ) then
          use_symm_ekin_density = .true.
!          use_asymm_ekin_density = .true.

          write(nfout,*) " !** use_symm_ekin_density is ",  use_symm_ekin_density
          write(nfout,*) " !** use_asymm_ekin_density is", use_asymm_ekin_density
       endif
    endif

    if ( sw_calc_ekin_density == OFF ) then
       sw_mix_charge_with_ekindens = OFF
       write(nfout,*) "!** sw_mix_charge_ekindens is turned off"
    endif

  end subroutine m_CtrlP_chkif_metagga

#ifdef LIBXC
  subroutine m_CtrlP_init_xctype_libxc( nfout )
    integer, intent(in) :: nfout
    integer :: ispin, i, j, itmp, iflag
    integer :: v1, v2, v3

    call xc_f03_version( v1, v2, v3 )
    write(nfout,*)
    write(nfout,'(A,I2,A,I2,A,I2)') " !** Libxc version ", v1, ".", v2, ".", v3

    if ( libxc_exch_name /= "" ) then
       libxc_exch_id = xc_f03_functional_get_number( libxc_exch_name )
    endif
    if ( libxc_corr_name /= "" ) then
       libxc_corr_id = xc_f03_functional_get_number( libxc_corr_name )
    endif

    if ( libxc_exch_id > 0 ) then
       xc_name_exch = xc_f03_functional_get_name(libxc_exch_id)
    endif
    if ( libxc_corr_id > 0 ) then
       xc_name_corr = xc_f03_functional_get_name(libxc_corr_id)
    endif

    write(nfout,*) " exchange    : ", trim(adjustl(xc_name_exch) )
    write(nfout,*) " correlation : ", trim(adjustl(xc_name_corr) )
    write(nfout,*)

    if ( libxc_exch_id +libxc_corr_id == 0 ) then
       call phase_error_with_msg( nfout, &
            &  'Exchange and/or Correlation should be set ', __LINE__,__FILE__ )
    endif

    if ( ndim_magmom == 1 ) then
       ispin = XC_UNPOLARIZED
    else
       ispin = XC_POLARIZED
    endif

    if ( libxc_exch_id > 0 ) then
       call xc_f03_func_init( xc_func_exch, libxc_exch_id, ispin )
       xc_info_exch = xc_f03_func_get_info(xc_func_exch)
       xc_family_exch = xc_f03_func_info_get_family(xc_info_exch)
    endif
    if ( libxc_corr_id > 0 ) then
       call xc_f03_func_init( xc_func_corr, libxc_corr_id, ispin )
       xc_info_corr = xc_f03_func_get_info(xc_func_corr)
       xc_family_corr = xc_f03_func_info_get_family(xc_info_corr)
    endif

    if ( xc_family_exch == XC_FAMILY_MGGA .or. xc_family_corr == XC_FAMILY_MGGA ) then
       use_metagga = .true.
    endif
    if ( xc_family_exch == XC_FAMILY_HYB_MGGA .or. &
         &   xc_family_corr == XC_FAMILY_HYB_MGGA ) then
       use_metagga = .true.
    endif

    if ( use_metagga ) then
       sw_calc_ekin_density = ON

       write(nfout,*) " !** use_metagga is ", use_metagga
       write(nfout,'(A,I3)') "!** sw_calc_ekin_density is ", sw_calc_ekin_density

       if ( .not. use_modeled_ekin_density  ) then
          use_symm_ekin_density = .true.
!          use_asymm_ekin_density = .true.
          write(nfout,*) " !** use_symm_ekin_density is ",  use_symm_ekin_density
          write(nfout,*) " !** use_asymm_ekin_density is", use_asymm_ekin_density
       endif

       vtau_exists = .true.

    endif
    if ( use_metagga ) then
       if ( libxc_exch_id > 0 ) then
          call xc_f03_func_set_dens_threshold( xc_func_exch, 1.0D-10 )
          call xc_f03_func_set_sigma_threshold( xc_func_exch, 1.0D-10 )
          call xc_f03_func_set_tau_threshold( xc_func_exch, 1.0D-10 )
       endif
       if ( libxc_corr_id > 0 ) then
          call xc_f03_func_set_dens_threshold( xc_func_corr, 1.0D-10 )
          call xc_f03_func_set_sigma_threshold( xc_func_corr, 1.0D-10 )
          call xc_f03_func_set_tau_threshold( xc_func_corr, 1.0D-10 )
       endif
    endif

!    if ( xc_family_exch == XC_FAMILY_HYB_GGA &
!         &          .or. xc_family_exch == XC_FAMILY_HYB_MGGA ) then
!       call set_hyrid_params_exch
!    endif

    if ( xc_family_corr == XC_FAMILY_HYB_GGA &
         &          .or. xc_family_corr == XC_FAMILY_HYB_MGGA ) then
       call set_hyrid_params_corr
    endif

    if ( libxc_exch_id > 0 ) then
       iflag = xc_f03_func_info_get_flags(xc_info_exch)
       xc_flag_exch = 0;  i = iflag
       Do j=num_item_flag, 0, -1
          itmp = i -2**j
          if ( itmp < 0 ) cycle
          i = itmp;   xc_flag_exch(j) = 1
       End Do
       if ( xc_flag_exch(0) == 0 ) write(nfout,*) &
            &        "!** Warning: potential-only exchange may not work properly."
    endif


    if ( libxc_corr_id > 0 ) then
       iflag = xc_f03_func_info_get_flags(xc_info_corr)
       xc_flag_corr = 0;  i = iflag
       Do j=num_item_flag, 0, -1
          itmp = i -2**j
          if ( itmp < 0 ) cycle
          i = itmp;   xc_flag_corr(j) = 1
       End Do
       if ( xc_flag_corr(0) == 0 ) write(nfout,*) &
            &        "!** Warning: potential-only correlation may not work properly."
    endif

    if ( use_metagga ) then
       if ( PAW_switch == ON ) then
          call phase_error_with_msg( nfout, &
               &  'mega-gga is not available when paw=ON',, __LINE__,__FILE__ )
       endif
       if ( istress /= 0 ) then
          call phase_error_with_msg( nfout, &
               &  'mega-gga is not avaiable when sw_stress==ON',, __LINE__,__FILE__ )
       endif
    endif

  end subroutine m_CtrlP_init_xctype_libxc
#endif
! ================== 13.0XX

! ----------- Written by T.Yamasaki, 1 Aug. 2014 ---------->>
!   This is a procedure to check if the solver which is now applied is appropriate to do integration in
!  the real space.
  logical function m_CtrlP_realspace_integ_OK()
    logical :: flag
    !flag = .false.
    !if(trim(solver_names_applied(1))=="RMM3".or.trim(solver_names_applied(2))=="RMM3") flag = .true.
    !if(trim(solver_names_applied(1))=="RMM2".or.trim(solver_names_applied(2))=="RMM2") flag = .true.
    !if(trim(solver_names_applied(1))=="RMM2P".or.trim(solver_names_applied(2))=="RMM2P") flag = .true.
    !if(trim(solver_names_applied(1))=="PKOSUGI".or.trim(solver_names_applied(2))=="PKOSUGI") flag = .true.
    !if(trim(solver_names_applied(1))=="PDAVIDSON".or.trim(solver_names_applied(2))=="PDAVIDSON") flag = .true.
    flag = .true.
    if(trim(solver_names_applied(1))=="MATDIAGON".or.trim(solver_names_applied(2))=="MATDIAGON") flag = .false.
    if(trim(solver_names_applied(1))=="SD"       .or.trim(solver_names_applied(2))=="SD")        flag = .false.
    if(trim(solver_names_applied(1))=="MSD"      .or.trim(solver_names_applied(2))=="MSD")       flag = .false.
    if(trim(solver_names_applied(1))=="lmSD"     .or.trim(solver_names_applied(2))=="lmSD")      flag = .false.
    if(trim(solver_names_applied(1))=="lmMSD"    .or.trim(solver_names_applied(2))=="lmMSD")     flag = .false.
    if(trim(solver_names_applied(1))=="CG"       .or.trim(solver_names_applied(2))=="CG")        flag = .false.
    if(trim(solver_names_applied(1))=="lmeazyCG" .or.trim(solver_names_applied(2))=="lmeazyCG")  flag = .false.
    if(trim(solver_names_applied(1))=="DAVIDSON" .or.trim(solver_names_applied(2))=="DAVIDSON")  flag = .false.
    m_CtrlP_realspace_integ_OK = flag
  end function m_CtrlP_realspace_integ_OK

!   This function is a procedure to check if the solvers applied through the job are all appropriate
!  to do the realspace integration.
  logical function m_CtrlP_rspace_integ_all_OK()
    integer :: i
    logical :: flag
    flag = .true.
    do i = 1, n_WF_solvers_all
       if(w_solver(i)%solver == lmMSD)         flag = .false.
       if(w_solver(i)%solver == MSD)           flag = .false.
       if(w_solver(i)%solver == lmSD)          flag = .false.
       if(w_solver(i)%solver == SD)            flag = .false.
       if(w_solver(i)%solver == MATRIXDIAGON)  flag = .false.
       if(w_solver(i)%solver == CG)            flag = .false.
       if(w_solver(i)%solver == DAVIDSON)      flag = .false.
    end do
    m_CtrlP_rspace_integ_all_OK = flag
  end function m_CtrlP_rspace_integ_all_OK
! <<--------------------------------------------------------

  subroutine m_CtrlP_set_in_line_minimization(onoffswitch)
    integer, intent(in) :: onoffswitch
    if(onoffswitch == ON) then
       in_line_minimization = .true.
    else if(onoffswitch == OFF) then
       in_line_minimization = .false.
    else
    end if
  end subroutine m_CtrlP_set_in_line_minimization

  subroutine m_CtrlP_set_icond(Iset)
    integer, intent(in) :: Iset
    icond = Iset
    if(icond < -2 .or. 2 <icond) icond = AUTOMATIC
  end subroutine m_CtrlP_set_icond

  subroutine m_CtrlP_check_dos_method(nfout)
    !    Coded by T. Yamasaki, 2017/04/21
    integer, intent(in) :: nfout
    if(dos_method == TETRAHEDRON .and. way_ksample /= MESH) then
       dos_method = Gauss_distrib_func
       ldos_method = Gauss_distrib_func
       if(ipri>=1) then
          write(nfout,999) " dos"
          write(nfout,999) "ldos"
       end if
    end if
999 format('!CtrlP ',a4,' method = Gauss_distrib_func, changed from TETRAHEDRON')
  end subroutine m_CtrlP_check_dos_method

  subroutine m_CtrlP_reread_edelta(nfout)
    integer, intent(in) :: nfout
    integer :: iret
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock
    real(kind=DP) :: dret
    iret = f_selectTop()
    if(f_selectBlock(tag_accuracy)==0) then
      if(f_selectBlock(tag_scf_convergence)==0) then
        call m_CtrlP_rd_val(nfout, tag_delta_total_energy, 'hartree' &
                         & , edelta, .true.)
        edelta_ontheway = edelta
        call m_CtrlP_rd_val(nfout, tag_succession, mtimes_convergence_scf, .true.)
        iret = f_selectParentBlock()
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine m_CtrlP_reread_edelta

  subroutine m_CtrlP_reread_max_force(nfout)
    integer, intent(in) :: nfout
    integer :: iret
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock
    real(kind=DP) :: dret
    iret = f_selectTop()
    if(f_selectBlock(tag_accuracy)==0) then
      if(f_selectBlock(tag_force_convergence)==0) then
        call m_CtrlP_rd_val(nfout, tag_max_force, 'hartree/bohr', forccr, .true.)
        iret = f_selectParentBlock()
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine m_CtrlP_reread_max_force

  subroutine m_CtrlP_reread_max_iteration(nfout)
    integer, intent(in) :: nfout
    integer :: iret
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock
    real(kind=DP) :: dret
    logical :: done_something
    iret = f_selectTop()
    if(f_selectBlock(tag_control)==0) then
      call m_CtrlP_rd_val(nfout, tag_max_iteration, max_total_scf_iteration &
      & , .true., done_something)
      if(done_something) max_TS_iteration_is_given = .true.

      call m_CtrlP_rd_val(nfout, tag_max_scf_iteration, max_scf_iteration &
      & , .true., done_something)
      if(done_something) max_scf_iteration_is_given = .true.

      call m_CtrlP_rd_val(nfout, tag_max_mdstep, max_mdstep &
      & , .true., done_something)
      if(done_something) max_mdstep_is_given = .true.

      iret = f_selectParentBlock()
    endif
  end subroutine m_CtrlP_reread_max_iteration

  subroutine m_CtrlP_reread_cutoff_wf(nfout)
    integer, intent(in) :: nfout
    integer :: iret
    integer :: f_selectParentBlock, f_selectTop, f_selectBlock, f_getRealValue
    real(kind=DP) :: dret
    logical :: tf
    if(m_CtrlP_gmax_changed()) return
    iret = f_selectTop()
    if(f_selectBlock(tag_accuracy)==0) then
      tf =  f_getRealValue( tag_cke_wavefunctions, dret,  "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf, dret,  "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf2, dret, "rydberg") == 0
      if(.not.tf) tf = f_getRealValue( tag_cke_wf3, dret, "rydberg") == 0
      if(tf) then
        if(dret < 0) call phase_error_with_msg(nfout,' !! illegal input value of cutoff_energy_for_wavefunctions'&
                                              ,__LINE__,__FILE__)
        gmax_buf = sqrt(dret)
        if(abs(gmax_buf-gmax)>1.e-15) then
          write(nfout,'(a,f10.3,a,f10.3)') '! ** F_INP_MOD changed &
          &cutoff_wf to ',dret,' from ',gmax*gmax
          write(nfout,'(a)') '! ** F_INP_MOD this modification will take effect &
          &after the current SCF iteration has converged'
          cutoff_wf_changed = .true.
        endif
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine m_CtrlP_reread_cutoff_wf

  function m_CtrlP_gmax_changed() result(res)
    logical :: res
    res = cutoff_wf_changed
    return
  end function m_CtrlP_gmax_changed

  subroutine m_CtrlP_reset_gmax_changed()
    cutoff_wf_changed = .false.
  end subroutine m_CtrlP_reset_gmax_changed

  subroutine rd_val_real(nfout, tag, unit_f, val, reread, done_something, &
  &          positive, msg)
    integer, intent(in) :: nfout
    character(len=*), intent(in) :: tag
    character(len=*), intent(in) :: unit_f
    real(kind=DP), intent(inout) :: val
    logical, intent(in),  optional :: reread
    logical, intent(out), optional :: done_something
    logical, intent(in),  optional :: positive
    character(len=*), intent(in), optional :: msg
    real(kind=DP) :: dret
    integer :: f_getRealValue
    logical :: tf, rr, posi
    tf = f_getRealValue(tag,dret,unit_f)==0
    if(present(done_something)) then
      done_something = .false.
    endif
    rr = .false.
    if(present(reread)) rr = reread
    posi = .false.
    if(present(positive)) posi = positive
    if(tf) then
      if(.not.rr) then
         if(.not.(posi .and. dret<0)) then
           val = dret
           if(present(done_something)) then
            done_something = .true.
           endif
         endif
      else
         if(.not.(posi .and. dret<0)) then
           if(abs(val-dret)>1.e-15) then
             if(present(msg)) then
               write(nfout,'(a,e12.5,a,e12.5)') '! ** F_INP_MOD changed ' &
               & //trim(msg)//' to ',dret,' from ',val
             else
               write(nfout,'(a,e12.5,a,e12.5)') '! ** F_INP_MOD changed ' &
               & //trim(tag)//' to ',dret,' from ',val
             endif
             val = dret
             if(present(done_something)) then
               done_something = .true.
             endif
           endif
        endif
      endif
    endif
  end subroutine rd_val_real

  subroutine rd_val_int(nfout, tag, val, reread, done_something, &
  &          positive, msg)
    integer, intent(in) :: nfout
    character(len=*), intent(in) :: tag
    integer, intent(inout) :: val
    logical, intent(in),  optional :: reread
    logical, intent(out), optional :: done_something
    logical, intent(in),  optional :: positive
    character(len=*), intent(in), optional :: msg
    integer :: iret, f_getIntValue
    logical :: tf, rr
    logical :: posi
    tf = f_getIntValue(tag,iret)==0
    if(present(done_something)) then
      done_something = .false.
    endif
    rr = .false.
    if(present(reread)) rr = reread
    posi = .false.
    if(present(positive)) posi = positive
    if(tf) then
      if(.not.rr) then
         if(.not. (posi .and. iret<0)) then
           val = iret
           if(present(done_something)) then
              done_something = .true.
           endif
         endif
      else
         if(.not. (posi .and. iret<0)) then
           if(val /= iret) then
             if(present(msg)) then
               write(nfout,'(a,i8,a,i8)') '! ** F_INP_MOD changed ' &
               & //trim(msg)//' to ',iret,' from ',val
             else
               write(nfout,'(a,i8,a,i8)') '! ** F_INP_MOD changed ' &
               & //trim(tag)//' to ',iret,' from ',val
             endif
             val = iret
             if(present(done_something)) then
                done_something = .true.
             endif
           endif
         endif
      endif
    endif
  end subroutine rd_val_int

  subroutine rd_val_str(nfout, tag, mode, val, reread, done_something,msg)
    integer, intent(in) :: nfout
    character(len=*), intent(in) :: tag
    integer, intent(in) :: mode
    character(len=*), intent(inout) :: val
    logical, intent(in),  optional :: reread
    logical, intent(out), optional :: done_something
    character(len=*), intent(in), optional :: msg
    integer :: f_getStringValue
    character(len=FMAXVALLEN) :: rstr
    logical :: tf, rr
    tf = f_getStringValue(tag, rstr, mode) == 0
    if(present(done_something)) then
      done_something = .false.
    endif
    rr = .false.
    if(present(reread)) rr = reread
    if(tf) then
      if(.not.rr) then
         val = rstr
         if(present(done_something)) then
            done_something = .true.
         endif
      else
         if(val /= rstr) then
           if(present(msg)) then
             write(nfout,'(a,i8,a,i8)') '! ** F_INP_MOD changed ' &
             & //trim(msg)//' to ',trim(rstr),' from ',val
           else
             write(nfout,'(a,i8,a,i8)') '! ** F_INP_MOD changed ' &
             & //trim(tag)//' to ',trim(rstr),' from ',val
           endif
           val = rstr
           if(present(done_something)) then
              done_something = .true.
           endif
         endif
      endif
    endif
  end subroutine rd_val_str

  subroutine m_CtrlP_set_conftag(tag)
    integer, intent(in) :: tag
    configuration_tag = tag
  end subroutine m_CtrlP_set_conftag

  integer function m_CtrlP_get_conftag()
    m_CtrlP_get_conftag = configuration_tag
    return
  end function m_CtrlP_get_conftag

end module m_Control_Parameters
