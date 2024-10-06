!=======================================================================
!
!  PROGRAM  PHASE/0 2023.01
!
!  MODULE: m_Crystal_Structure
!
!  AUTHORS: T. Yamasaki, K. Betsuyaku,    August/20/2003
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
module m_Crystal_Structure
!      (m_CS)
! $Id: m_Crystal_Structure.F90 630 2020-07-29 09:03:26Z ktagami $
!!$  use m_Files,              only : nfout,nfopgr,nfmatbp &
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only : ipri, af, m_CtrlP_set_af, m_CtrlP_set_nspin_and_af &
#ifndef ENABLE_ESM_PACK
       &                         , ipriinputfile, ipri_spg, printable, nspin, iprisym
#else
       &                         , ipriinputfile, ipri_spg, printable, nspin &
       &                         , sw_esm,esm_z1_defined,esm_z1,esm_w, iprisym
#endif
  use m_Const_Parameters,   only : DP,PAI,PAI2,CRDTYP,CARTS,BUCS,GENERAL &
       &                         , GENERAL_LARGER, HCP,SIMPLE_CUBIC, BCC, FCC &
       &                         , DIAMOND, HEXAGONAL, TRIGONAL, FMAXVALLEN &
       &                         , TRIGONAL_lattice, HEXAGONAL_lattice, PRIMITIVE_lattice &
       &                         , FCC_lattice, BCC_lattice, BottomCenteredCubic_lattice &
       &                         , num_d6h, num_oh, d6h_symbol, oh_symbol &
       &                         , PARA, NONMAG, ANTIFERRO, FERRO, NOCONV, LOWER, UPPER, ON, OFF &
       &                         , AUTOMATIC, MANUAL, BRAVAIS, PRIMITIVE, YES, NO &
       &                         , WHOLE, INITIALLY
  use m_Parallelization,    only : mype

! ================================== added by K. Tagami =============== 11.0 & 13.2S
  use m_Control_Parameters,  only : noncol, SpinOrbit_mode
  use m_Const_Parameters,    only : NONCOLLINEAR, Neglected
! ===================================================================== 11.0 & 13.2S

! =============================== KT_add ========================= 13.0U
  use m_Control_Parameters,  only : num_projectors
  use m_Const_Parameters,    only : MAG_MOMENT_VALS_GLOBAL, MAG_MOMENT_DIREC_GLOBAL, &
       &                            MAG_MOMENT_VALS_LOCAL,  MAG_MOMENT_DIREC_LOCAL, &
       &                            MAG_MOMENT_VALS_OCCMAT, &
       &                            MAG_MOMENT_DIREC_HARDPART, &
       &                            ABRUPT, STEPWISE, LINEAR
! ================================================================ 13.0U

! === KT_add ==== 13.1R
  use m_Control_Parameters,  only : sw_keep_symmetry_strict, icond
  use m_Const_Parameters,    only : COORDINATE_CONTINUATION, CONTINUATION, FIXED_CHARGE_CONTINUATION
! =============== 13.1R

  implicit none

!!$  integer :: input_coordinate_system             ! ncord

  integer ::                          inv = 0    ! flag for moving the origin (0:nomove, 1:move)
  integer :: inversion_symmetry = 0              ! ninv , same as inv.  {0|1} 0=noinversion nomove,1=inversion move
!!$  integer :: ntyp = 1
!!$  integer :: natm = 1
  integer :: nbztyp = GENERAL
  integer :: nbztyp_spg = 0      

  !!$integer :: symmetry_method = AUTOMATIC
  integer :: symmetry_method = MANUAL
  integer :: sw_supercell_symmetry = ON
  real(kind=DP) :: misalignment = 0.1D0    ! Bohr

  real(kind=DP) ::                          univol, rvol
  real(kind=DP) ::                          univol_prev
  real(kind=DP) ::                          univol_super, rvol_super
  real(kind=DP) ::                          univol_prim, rvol_prim
  real(kind=DP) ::                          symmetry_check_criterion = 1.d-12
  real(kind=DP), dimension(3,3) ::          altv, rltv
  real(kind=DP), dimension(3,3) ::          altv_super, rltv_super
  real(kind=DP), dimension(3,3) ::          altv_prim, rltv_prim

! ===================================== Added by K. Tagami ==========
  integer :: kt_iuctype
! ===================================================================

  ! ---- nfspg ---
  integer ::                                      il         ! lattice system
  character, private, dimension(-1:4)  ::         Char_lattice_system
  data Char_lattice_system /'R','H','P','F','I','C'/ 
                   ! Trigonal|Hexagonal|Simple|FaceCentered|BodyCentered|BaseCentered(C-basecentered)
                   ! Each of them corresponds to il= -1, 0, 1, 2, 3, and 4, respectively.
  integer ::                                      imag = NONMAG ! para(nonmag), antiferro, ferro
  integer ::                                      sw_fix_total_spin = NO ! effective when imag==FERRO
  integer ::                                      sw_fix_total_spin_in = NO
  integer ::                                      sw_bandgap_constraint = NO
  integer ::                                      spin_fix_period  = WHOLE
  real(kind=DP) ::                                total_spin = 0.d0 ! effecitve when imag=FERRO
  real(kind=DP) ::                                a,b,c      ! lattice parameter (length)
  real(kind=DP) ::                                ca,cb,cc   ! lattice parameter (cos(angle))
  real(kind=DP),dimension(3) ::                   Bravais_lattice_length ! unit:Bohr
  real(kind=DP),dimension(3) ::                   Bravais_lattice_angle  ! unit:degrees
  integer ::                                      ngen = 1   ! # of generators
  integer, allocatable, dimension(:)::            igen       ! d(ngen) rotation of generater
  integer, allocatable, dimension(:,:,:) ::       jgen       ! d(2,3,ngen) nonprimitive translation vector of genertor 
  integer                              ::         iaf=1      ! rotation of generater for AF case
  integer, dimension(2,3)              ::         jaf=reshape((/0,1,0,1,0,1/),(/2,3/)) ! nonprimitive translation vector of generator for AF case
!  integer, dimension(2,3)              ::         jaf ! nonprimitive translation vector of generator for AF case
  
  integer ::                                      nopr ! (=ng1) # of group elements
  integer ::                                      phonon_nopr ! for Phonon code
  real(kind=DP),allocatable,dimension(:,:,:),target :: tau ! (=ta1) d(3,nopr,CRDTYP), nonprimitive translation vector (A system)
  real(kind=DP),allocatable,dimension(:,:,:),target :: op ! (=ra1) d(3,3,nopr),rotation matrix in real space (A system)

  integer ::                                ngen_tl = 0 ! # of additional generators for supercell
!!$  integer, allocatable,dimension(:)::       igen_tl ! d(ngen_tl) rotation of generater
  integer, allocatable,dimension(:,:,:) ::  jgen_tl ! d(2,3,ngen_tl) nonprimitive translation vector of generator
  integer ::                                nopr_tl = 0 ! (=ng1) # of additional group elements for supercell
  real(kind=DP),allocatable,dimension(:,:,:):: tau_tl !(=ta1)d(3,ngen_tl,CRDTYP),nonprimitive translation vector (A system)
  real(kind=DP),allocatable,dimension(:,:,:):: op_tl !(=ra1)d(3,3,ngen_tl),rotation matrix in real space (A system)
  real(kind=DP), allocatable, dimension(:,:,:) :: sa1  ! d(3,3,nopr), rotation matrix in reciprocal space (A system)
  integer, dimension(48)                       :: ig01 ! symmetry operation numbers
  real(kind=DP) :: b2pmat(3,3) = 0.d0 ! transformation matrix (Bravais -> Primitive)
! ==================== Modified by K. Tagami =======
!  real(kind=DP) :: p2bmat(3,3) = 0.d0 ! transformation matrix (Bravais -> Primitive)
  real(kind=DP) :: p2bmat(3,3) = 0.d0 ! transformation matrix (Bravais <- Primitive)
! ===============================================


  character(len=9) :: pg_symbol_system = 'cubic'

 ! strain tensor
  real(kind=DP), dimension(1:3,1:3) :: strain
  integer,save                      :: sw_strained_cell = OFF
  character(len("strain")),private,parameter ::    tag_strain    = "strain"
  character(len("sw_strained_cell")),private,parameter :: tag_sw_strained_cell = "sw_strained_cell"

 ! reference cell for band unfolding
  character(len("reference_cell")), private,parameter :: &
       &                          tag_reference_cell = 'reference_cell'
  integer :: refcell_unit_type = PRIMITIVE
  real(kind=DP) ::  altv_refcell(3,3), rltv_refcell(3,3)
  real(kind=DP) ::  univol_refcell, rvol_refcell
  real(kind=DP) :: tmpmat(3,3)
  real(kind=DP) :: unit_cell_scale = 1.d0

  ! supercell
  integer :: sw_supercell = OFF
  integer :: supercell_unit_type = BRAVAIS
  integer :: n1_sc=1, n2_sc=1, n3_sc=1
  integer :: nlpnt = 1
  real(kind=DP), allocatable :: lpnt(:,:) ! dim(nlpnt,3)
  character(len("supercell")), private,parameter :: tag_supercell = 'supercell'
  character(len("sw_supercell")), private,parameter :: tag_sw_supercell = 'sw_supercell'
  character(len("size")), private,parameter :: tag_size = 'size'
  character(len("n1")), private,parameter :: tag_n1 = 'n1'
  character(len("n2")), private,parameter :: tag_n2 = 'n2'
  character(len("n3")), private,parameter :: tag_n3 = 'n3'

  ! space group
  character(len("method")),private,parameter ::    tag_method    = "method"
  character(len("automatic")),private,parameter ::    tag_automatic    = "automatic"
  character(len("manual")),private,parameter ::    tag_manual    = "manual"
  character(len("misalignment")),private,parameter :: tag_misalignment = "misalignment"
  character(len("tspace")),private,parameter ::    tag_tspace    = "tspace"
  ! sw_supercell_symmetry
  character(len("sw_supercell_symmetry")),private,parameter :: &
       &                                           tag_sw_supercell_symmetry = "sw_supercell_symmetry"
  ! ---> input tag terms
  ! --- Structure ---
  character(len("structure")),public,parameter ::  tag_structure = "structure"
  character(len("unit_cell")),private,parameter :: tag_unit_cell = "unit_cell"
  character(len("unit_cell_type")),private,parameter :: tag_unit_cell_type = "unit_cell_type"
  character(len("a_vector")), private,parameter :: tag_avector  = 'a_vector'
  character(len("b_vector")), private,parameter :: tag_bvector  = 'b_vector'
  character(len("c_vector")), private,parameter :: tag_cvector  = 'c_vector'
  character(len("a")), private,parameter ::        tag_a        = 'a'
  character(len("b")), private,parameter ::        tag_b        = 'b'
  character(len("c")), private,parameter ::        tag_c        = 'c'
  character(len("alpha")),private,parameter ::     tag_alpha    = 'alpha'
  character(len("beta")), private,parameter ::     tag_beta      = 'beta'
  character(len("gamma")), private,parameter ::    tag_gamma     = 'gamma'
  character(len("unit_cell_scale")), private, parameter :: tag_unit_cell_scale = "unit_cell_scale"

  character(len("symmetry")),private,parameter ::  tag_symmetry  = "symmetry"
  character(len("symmetry_check_criterion")),private,parameter :: tag_symmetry_check_criterion = "symmetry_check_criterion"
  character(len("crystal_structure")),private,parameter ::  tag_crystal_structure = "crystal_structure"
  character(len("crystal")),private,parameter ::            tag_crystal           = "crystal"
  character(len("simple_cubic")),private,parameter ::       tag_simple_cubic      = "simple_cubic"
  character(len("diamond")),private,parameter ::            tag_diamond           = "diamond"

  character(len("lattice_system")),private,parameter ::     tag_lattice_system = "lattice_system"
  character(len("system")),private,parameter ::             tag_system    = "system"
  character(len("rhombohedral")),private,parameter ::       tag_rhombohedral = "rhombohedral"
  character(len("trigonal")),private,parameter ::           tag_trigonal  = "trigonal"
  character(len("hexagonal")),private,parameter ::          tag_hexagonal = "hexagonal"
  character(len("hcp")),private,parameter ::                tag_hcp       = "hcp"
  character(len("primitive")),private,parameter ::          tag_primitive = "primitive"
  character(len("simple")),private,parameter ::             tag_simple    = "simple"
  character(len("facecentered_cubic")),private,parameter::  tag_facecentered = "facecentered_cubic"
  character(len("bodycentered_cubic")),private,parameter::  tag_bodycentered = "bodycentered_cubic"
  character(len("bcc")),private,parameter ::                tag_bcc       = "bcc"
  character(len("fcc")),private,parameter ::                tag_fcc       = "fcc"
  character(len("onefacecentered")),private,parameter ::    tag_onefacecentered = "onefacecentered"
  character(len("bottomcentered")),private,parameter ::     tag_bottomcentered = "bottomcentered"
  character(len("basecentered")),private,parameter ::       tag_basecentered = "basecentered"
  character(len("num_generators")),private,parameter ::     tag_num_generators = "num_generators"
  character(len("generators")),private,parameter ::         tag_generators = "generators"
  character(len("additional_translations")),private,parameter::tag_additional_tl = "additional_translations"
  character(len("af_generator")),private,parameter ::       tag_af_generator = "af_generator"
  character(len("rotation")),private,parameter ::           tag_rotation  = "rotation"
  character(len("tx")),private,parameter ::                 tag_tx        = "tx"
  character(len("ty")),private,parameter ::                 tag_ty        = "ty"
  character(len("tz")),private,parameter ::                 tag_tz        = "tz"
  character(len("sw_inversion")),private,parameter ::       tag_sw_inversion = "sw_inversion"
  character(len("bandgap_constraint")),private,parameter :: tag_bandgap_constraint  = "bandgap_constraint"
  character(len("nspin")),private,parameter ::              tag_nspin     = "nspin"
  logical ::                                                tag_nspin_is_found = .false.
  character(len("magnetic_state")),private,parameter ::     tag_magnetic_state  = "magnetic_state"
  logical ::                                                tag_magnetic_state_is_found = .false.
  character(len("para")),private,parameter ::               tag_para      = "para"
  character(len("nonmagnetic")),private,parameter ::        tag_nonmagnetic = "nonmagnetic"
  character(len("nonmag")),private,parameter ::             tag_nonmag    = "nonmag"
  character(len("none")),private,parameter ::               tag_none      = "none"
  character(len("af")),private,parameter ::                 tag_af        = "af"
  character(len("antiferro")),private,parameter ::          tag_antiferro = "antiferro"
  character(len("ferro")),private,parameter ::              tag_ferro     = "ferro"
  character(len("ferromagnetic_state")),private,parameter:: tag_ferromagnetic_state = "ferromagnetic_state"
  character(len("magnetic")),private,parameter ::           tag_magnetic  = "magnetic"
  character(len("mag")),private,parameter ::                tag_mag       = "mag"

  character(len("sw_non_symmorphic")),private,parameter  :: tag_sw_non_symmorphic="sw_non_symmorphic"
  character(len("cpumax")), private, parameter           :: tag_cpumax="cpumax"
  character(len("sw_read_symmetry_from_file")), private, parameter :: tag_sw_read_symmetry_from_file &
                                                         & = "sw_read_symmetry_from_file"

  character(len("sw_fix_total_spin")),private,parameter ::  tag_sw_fix_total_spin = "sw_fix_total_spin"
  character(len("spin_fix_period")),private,parameter   ::  tag_spin_fix_period = "spin_fix_period"
  character(len("whole")),private,parameter ::              tag_whole     = "whole"
  character(len("initially")),private,parameter ::          tag_initially = "initially"
  character(len("total_spin")),private,parameter ::         tag_total_spin = "total_spin" 

  character(len("primitive")),private,parameter ::          val_primitive     = "primitive"
  character(len("bravais")),private,parameter ::          val_bravais  = "bravais"

! --> T. Yamasaki 18th Aug. 2009 
  integer, parameter :: len_str = 132
  character(len=len_str) :: str
  character(len("fix_spin_status")),private,parameter ::    tag_fix_spin_status = "fix_spin_status" 
! <--

! ================================== added by K. Tagami =============== 11.0
!
!!  -- NonCollinear --
!
  character(len("noncollinear")),private,parameter ::   tag_noncollinear = "noncollinear"
  character(len("noncol")),      private,parameter ::   tag_noncol = "noncol"
!
  character(len("noncollinear_state")),private,parameter :: &
       &                         tag_noncollinear_state = "noncollinear_state"
!
!
  character(len("use_magnetic_symmetry")), private, parameter :: &
       &             tag_use_magnetic_symmetry = "use_magnetic_symmetry"
  character(len("use_mag_symm")), private, parameter :: &
       &             tag_use_mag_symm = "use_mag_symm"
  character(len("reduce_sym_by_magmom")), private, parameter :: &
       &             tag_reduce_sym_by_magmom = "reduce_sym_by_magmom"
  character(len("reduce_sym_by_orbital")), private, parameter :: &
       &             tag_reduce_sym_by_orbital = "reduce_sym_by_orbital"

! --- GGA noncl -
  character(len("constrain_on_grad_correction")), private, parameter :: &
       &             tag_constrain_on_grad_corr = "constrain_on_grad_correction"
! ---- PAW noncl --
  character(len("level_of_projection_paw_charge")), private, parameter :: &
       &             tag_level_projection_paw_charge = "level_of_projection_paw_charge"

! --- magnetic moment
  character(len("axis")),private,parameter :: tag_axis = "axis"
  character(len("direction")),private,parameter :: tag_direction = "direction"
  character(len("moment")),private,parameter :: tag_moment = "moment"
  character(len("magnetic_moment")),private,parameter :: tag_magnetic_moment &
       &                                       = "magnetic_moment"

  character(len("mx")),private,parameter :: tag_mx = "mx"
  character(len("my")),private,parameter :: tag_my = "my"
  character(len("mz")),private,parameter :: tag_mz = "mz"
  character(len("mdx")),private,parameter :: tag_mdx = "mdx"
  character(len("mdy")),private,parameter :: tag_mdy = "mdy"
  character(len("mdz")),private,parameter :: tag_mdz = "mdz"

  character(len("norm")),private,parameter :: tag_norm = "norm"
  character(len("theta")),private,parameter :: tag_theta = "theta"
  character(len("phi")),private,parameter :: tag_phi = "phi"
!
  real(kind=DP) :: norm
!
  integer :: sw_use_magnetic_symmetry = ON
  integer :: sw_reduce_sym_by_magmom = ON
  integer :: sw_reduce_sym_by_orbital = OFF
!!
  integer :: sw_constrain_on_grad_correction = OFF
  integer :: level_of_projection_paw_charge = 1
!
  real(kind=DP) :: mag_moment0_global(3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  real(kind=DP) :: mag_direc0_global(3) = (/ 0.0d0, 0.0d0, 1.0d0 /)
! ================================================================== 11.0

! === KT_add === 2014/08/04
!
! -- lowering contanumation in the case of noncol+GGA
!
  character(len("sw_neglect_low_helicity")),private,parameter :: &
       &          tag_sw_neglect_low_helicity = "sw_neglect_low_helicity"
  character(len("sw_neglect_low_vorticity")),private,parameter :: &
       &          tag_sw_neglect_low_vorticity = "sw_neglect_low_vorticity"
  character(len("threshold_helicity")),private,parameter :: &
       &          tag_threshold_helicity = "threshold_helicity"
  character(len("threshold_vorticity")),private,parameter :: &
       &          tag_threshold_vorticity = "threshold_vorticity"
!
  integer :: sw_neglect_low_helicity = OFF
  integer :: sw_neglect_low_vorticity = OFF
!
  integer :: sw_monitor_magnetic_vorticity = OFF
!
  real(kind=DP) :: threshold_helicity = 1.0D-6
  real(kind=DP) :: threshold_vorticity = 1.0D-4
! ============== 2014/08/04

! === KT_add ==== 2014/09/26
  character(len("sw_neglect_magmom")),private,parameter :: &
       &          tag_sw_neglect_magmom = "sw_neglect_magmom"
  integer :: sw_neglect_magmom = OFF
! =============== 2014/09/26

  character(len("sw_allow_mag_sym_inversion")),private,parameter :: &
       &          tag_sw_allow_mag_sym_inversion = "sw_allow_mag_sym_inversion"
  character(len("sw_allow_improper_rotation")),private,parameter :: &
       &          tag_sw_allow_improper_rotation = "sw_allow_improper_rotation"
  character(len("sw_allow_frac_translation")),private,parameter :: &
       &          tag_sw_allow_frac_translation = "sw_allow_frac_translation"
  integer :: sw_allow_mag_sym_inversion = ON
  integer :: sw_allow_improper_rotation = ON
  integer :: sw_allow_frac_translation = ON

! ==== KT_add ==== 2014/08/09 & 13.2S
  character(len("magnetic_restriction")),private,parameter :: &
       &          tag_magnetic_restriction = "magnetic_restriction"
  character(len("sw_fix_global_quantz_axis")),private,parameter :: &
       &          tag_sw_fix_global_quantz_axis = "sw_fix_global_quantz_axis"
  character(len("global_quantz_axis")),private,parameter :: &
       &        tag_global_quantz_axis = "global_quantz_axis"
  character(len("sx")),private,parameter :: tag_sx = "sx"
  character(len("sy")),private,parameter :: tag_sy = "sy"
  character(len("sz")),private,parameter :: tag_sz = "sz"
!
  integer :: sw_fix_global_quantz_axis = OFF
  real(kind=DP) :: Global_Quantz_Axis_Fixed(3) = (/0.0d0, 0.0d0, 1.0d0/)
! ================ 2014/08/09 & 13.2S

! =========================== KT_add ===================== 13.0U
!
! -- constraint on magnetic moment --
!
  character(len("magnetic_constraint")), private, parameter :: &
       &             tag_magnetic_constraint = "magnetic_constraint"
  character(len("constraint_type")), private, parameter :: &
       &             tag_constraint_type = "constraint_type"

  character(len("vals_global")), private, parameter :: &
       &                             tag_moment_vals_global = "vals_global"
  character(len("direc_global")), private, parameter :: &
       &                             tag_moment_direc_global = "direc_global"
  character(len("vals_local")), private, parameter :: &
       &                             tag_moment_vals_local = "vals_local"
  character(len("direc_local")), private, parameter :: &
       &                             tag_moment_direc_local = "direc_local"
  character(len("vals_occmat")), private, parameter :: &
       &                             tag_moment_vals_occmat = "vals_occmat"
  character(len("direc_hardpart")), private, parameter :: &
       &                             tag_moment_direc_hardpart = "direc_hardpart"
!
  character(len("damping_method")), private, parameter :: &
       &                            tag_damping_method = "damping_method"
  character(len("abrupt")), private, parameter :: &
       &                            tag_damp_abrupt = "abrupt"
  character(len("stepwise")), private, parameter :: &
       &                            tag_damp_stepwize = "stepwise"
  character(len("linear")), private, parameter :: &
       &                            tag_damp_linear = "linear"
! --
  character(len("lambda")), private, parameter ::  tag_lambda = "lambda"
  character(len("num_intermid_lambda")), private, parameter :: &
       &                          tag_num_intermid_lambda = "num_intermid_lambda"

  character(len("edelta_change_lambda_first")), private, parameter :: &
       &                 tag_edelta_change_lambda_first = "edelta_change_lambda_first"
  character(len("edelta_change_lambda_last")), private, parameter :: &
       &                 tag_edelta_change_lambda_last = "edelta_change_lambda_last"
!
  character(len("max_iter_elec_mag_constraint")), private, parameter :: &
       &                 tag_max_iter_elec_mag_constr = "max_iter_elec_mag_constraint"
  character(len("max_iter_ion_mag_constraint")), private, parameter :: &
       &                 tag_max_iter_ion_mag_constr = "max_iter_ion_mag_constraint"
  character(len("max_iter_cell_mag_constraint")), private, parameter :: &
       &                 tag_max_iter_cell_mag_constr = "max_iter_cell_mag_constraint"
!
  character(len("sw_fix_charge_after_constraint")), private, parameter :: &
       &              tag_sw_fix_charge_after_constr = "sw_fix_charge_after_constraint"

  integer, parameter :: nmax_intermid_lambda = 100
!
  integer :: sw_magnetic_constraint = OFF
  integer :: mag_constraint_type = 0
  integer :: damping_method_mag_constraint = ABRUPT
  integer :: num_intermid_lambda = 2
!
  integer :: max_iter_elec_mag_constraint = 50
  integer :: max_iter_ion_mag_constraint = 1
!  integer :: max_iter_cell_mag_constraint = 100
  integer :: max_iter_cell_mag_constraint = 1

  integer :: sw_fix_charge_after_constraint = OFF
!
  real(kind=DP) :: mag_constraint_lambda = 0.20d0
  real(kind=DP) :: edelta_change_lambda_first = 1.0D-4     ! hartree
  real(kind=DP) :: edelta_change_lambda_last = 1.0D-4     ! hartree
! ======================================================== 13.0U

! === KT_add === 2014/06/08
! Charged state
!
  character(len("charged_state")),private,parameter :: &
       &               tag_charged_state  = "charged_state"
  character(len("extra_charge")),private,parameter :: &
       &               tag_extra_charge  = "extra_charge"
  character(len("additional_charge")),private,parameter :: &
       &               tag_additional_charge  = "additional_charge"
!
  real(kind=DP) :: additional_charge = 0.0d0
! ============== 2014/06/08

! === KT_add === 13.2S
  character(len("sw_use_magnetic_symmetry")), private, parameter :: &
       &                tag_sw_use_magnetic_symmetry = "sw_use_magnetic_symmetry"
!  integer :: sw_spinorbit_force_theorem = off
  integer :: sw_spinorbit_second_variation = off
! ============== 13.2S

! =============================== added by K. Tagami ============ 12.0A
  logical :: use_altv_rltv = .false.
  character(len("use_altv_rltv")), private, parameter :: &
       &                tag_use_altv_rltv = "use_altv_rltv"
!
  logical :: gen_name_in_carts = .false.
  character(len("gen_name_in_carts")), private, parameter :: &
       &                tag_gen_name_in_carts = "gen_name_in_carts"
! =============================================================== 12.0A

! == KT_add == 13.1R
  integer, save :: nopr_previous
  real(kind=DP), allocatable :: op_previous(:,:,:), tau_previous(:,:,:)
! ============ 13.1R

  logical, private :: symmetry_not_printed  = .true.

  integer :: sw_non_symmorphic = ON
  real(kind=DP) :: cpumax_symsearch = -1.d0
  integer :: sw_read_symmetry_from_file = ON

contains
  subroutine m_CS_put_Bravais_lattice(length,angle)
    real(kind=DP), dimension(3) :: length, angle
    Bravais_lattice_length = length
    Bravais_lattice_angle  = angle
  end subroutine m_CS_put_Bravais_lattice

  subroutine m_CS_set_total_spin(nfout,total_spin_new)
    integer,       intent(in) :: nfout 
    real(kind=DP), intent(in) :: total_spin_new
    total_spin = total_spin_new
    if(ipri>=1) write(nfout,'(" total_spin_new = ", f12.6)') total_spin_new
  end subroutine m_CS_set_total_spin

  subroutine m_CS_free_fixed_total_spin(nfout)
    integer,       intent(in) :: nfout 
    sw_fix_total_spin = NO
    if(ipri >=1) write(nfout,'(" !** (new setting) sw_fix_total_spin = NO <<m_CS_free_fixed_total_spin>>")')
  end subroutine m_CS_free_fixed_total_spin

  subroutine m_CS_reset_fixed_total_spin(nfout)   ! 2024.03.27 T. Yamasaki
    integer,       intent(in) :: nfout 
    sw_fix_total_spin = sw_fix_total_spin_in
    if(ipri >=1) write(nfout,'(" !** (new setting) sw_fix_total_spin = sw_fix_total_spin_in = ",I3," <<m_CS_reset_fixed_total_spin>>")') &
         sw_fix_total_spin
  end subroutine m_CS_reset_fixed_total_spin

  subroutine m_CS_fix_total_spin(nfout)  ! T. Yamasaki 17th Aug. 2009
    integer,       intent(in) :: nfout 
    sw_fix_total_spin = YES
    if(ipri >=1) write(nfout,'(" !** (new setting) sw_fix_total_spin = YES <<m_CS_fix_total_spin>>")')
  end subroutine m_CS_fix_total_spin

  subroutine m_CS_rd_n(nfout)
    ! <m_CS_rd_n> reads the information of crystal system from an input file
    !    formatted in a new style.
    ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), Jun. 2003

    integer, intent(in) :: nfout
    character(len=FMAXVALLEN) :: rstr, dummy
    integer :: i, iret, rint, iuctype=BRAVAIS, ucinptype     !!! iuctype=0:bravais, 1:primitive;  ucinptype=0:length, 1:vector
    integer :: f_selectBlock, f_getStringValue, f_getIntValue
    integer :: f_selectParentBlock, f_getRealValue
    integer :: f_readUnitCell, f_readAfTauVec
    real(kind=DP),allocatable,dimension(:) :: avec,bvec,cvec
    logical :: number_is_given, prealloc, is_hexagonal, hext
    real(kind=DP) :: dret

    il = PRIMITIVE_lattice !!! K.Mae 040121
    if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** << m_CS_rd_n >>")')
    ! --- Structure ---
    if( f_selectBlock( tag_structure) == 0) then
       if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_structure --")')
       ! --- magnetic_state

       if( f_getIntValue( tag_nspin, iret) == 0) then
          call set_imag_from_nspin(iret)
          tag_nspin_is_found = .true.
       else
          tag_nspin_is_found = .false.
       end if
       if( f_getStringValue( tag_magnetic_state, rstr, LOWER) == 0) then ! == modified by K. Tagami == 11.0
          tag_magnetic_state_is_found = .true.
          if(tag_nspin_is_found ) then
             if(ipriinputfile>=1) then
                write(nfout,'(" !** Both tags ",a1,"nspin",a1," and ",a1,"magnetic_state",a1," are defined.")') &
                     & char(34),char(34),char(34),char(34)
                write(nfout,'(" !**   The tag ",a1,"nspin",a1," takes priority")') char(34),char(34)
             end if
          else
             call set_imag(rstr) ! imag = {PARA|ANTIFERRO|FERRO|NONCOLLINEAR}  == modified by K. Tagami == 11.0
          endif
       end if

       if(ipriinputfile >= 1 .and. printable) then
          if(.not.tag_magnetic_state_is_found .and. .not.tag_nspin_is_found) then
             write(nfout,'(" !** Neither ",a1,"nspin",a1," nor ",a1,"magnetic_state",a1 &
                  & ," is defined in the input file")') char(34),char(34),char(34),char(34)
             write(nfout,'(" !** nspin is set 1 (default)")')
          end if
!!$          write(nfout,'(" !** imag = ", i6)') imag
       end if
       call m_CtrlP_set_nspin_and_af(nfout,imag)

       ! --- bandgap_constraint
       if( f_getIntValue( tag_bandgap_constraint, iret) == 0) then
          sw_bandgap_constraint = iret
          if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** bandgap_constraint = ",i6)') &
               & sw_bandgap_constraint
          if(imag .ne. NONMAG .and. sw_bandgap_constraint .ne. NO) then
             sw_bandgap_constraint = OFF
             if(ipriinputfile >=1 .and. printable) write(nfout,'(" !** bandgap_constraint(reset) = ",i6)') &
                  & sw_bandgap_constraint
          end if
       end if

       if(imag == FERRO) then
          if( f_selectBlock(tag_ferromagnetic_state) == 0) then
             if(f_getIntValue( tag_sw_fix_total_spin,iret) == 0) sw_fix_total_spin_in = iret
             if(sw_fix_total_spin_in /= YES .and. sw_fix_total_spin_in /= NO) sw_fix_total_spin_in = NO
             sw_fix_total_spin = sw_fix_total_spin_in
             if(f_getRealValue(tag_total_spin,dret,'') == 0) total_spin = dret
             if( f_getStringValue( tag_spin_fix_period, rstr, LOWER) == 0) &
                  & call set_spin_fix_period(rstr) ! spin_fix_period = {WHOLE|INITIALLY}
             if( f_getIntValue(tag_spin_fix_period,iret)==0) spin_fix_period=iret
             iret = f_selectParentBlock()
             if(ipriinputfile >= 1) then
                write(nfout,'(" !** sw_fix_total_spin_in = ",i3)') sw_fix_total_spin_in
                write(nfout,'(" !** set_spin_fix_period = ",i3," : 0=WHOLE, 1=INITIALLY")') spin_fix_period
                write(nfout,'(" !** total_spin        = ",f15.10)') total_spin
             end if
          end if
       end if

! == KT_add == 13.2S
       if ( .not. noncol ) then
          sw_use_magnetic_symmetry = off
#if 0
          if( f_getIntValue( tag_sw_use_magnetic_symmetry,iret)==0 ) then
             sw_use_magnetic_symmetry = iret
          endif
          if ( sw_use_magnetic_symmetry == ON ) sw_spinorbit_second_variation = on
#endif
          if ( SpinOrbit_mode /= Neglected ) then
             sw_spinorbit_second_variation = on
             if ( imag==FERRO ) sw_use_magnetic_symmetry = on
          endif
       endif
! ============ 13.2S

! =========================== added by K. Tagami ===================== 11.0
       if ( imag == NONCOLLINEAR ) then
          if( f_selectBlock( tag_noncollinear_state ) == 0 ) then

             call set_sw_reduce_sym_by_magmom
             call set_sw_reduce_sym_by_orbital
             call set_sw_magnetic_symmetry
! ============ KT_mod ================== 13.0U
!             if( f_selectBlock( tag_magnetic_constraint ) == 0 ) then
!                call set_magnetic_constraint
!                iret = f_selectParentBlock()
!             end if
! =======================================13.0U
! ==== KT_add === 2014/09/26
             call set_sw_neglect_magmom
! =============== 2014/09/26
             call set_sw_constrain_on_grad_corr
             call set_level_projection_paw_charge

! ==== KT_add === 2014/08/04
             call set_threshold_lowering_contami     ! vorticity, helicity
! =============== 2014/08/04
             iret = f_selectParentBlock()
          end if
       end if
! ===================================================================== 11.0

! ==== KT_add === 2014/08/09
       call set_magnetic_restriction
! =============== 2014/08/09
! ====================== KT_add ======================= 13.0U
       if( f_selectBlock( tag_magnetic_constraint ) == 0 ) then
          call set_magnetic_constraint
          iret = f_selectParentBlock()
       end if
! ===================================================== 13.0U

! == KT_add === 2014/06/08
       if ( f_selectBlock( tag_charged_state ) == 0 ) then
          if (f_getRealValue(tag_extra_charge,dret,'') == 0) &
               &                        additional_charge = dret
          if (f_getRealValue(tag_additional_charge,dret,'') == 0) &
               &                        additional_charge = dret
          if ( additional_charge > 0.0 ) then
             write(nfout,*) "** The system is positively charged : ", additional_charge
          else if ( additional_charge < 0.0 ) then
             write(nfout,*) "** The system is negatively charged : ", -additional_charge
          endif
          iret = f_selectParentBlock()
       endif
! ============= 2014/06/08

       ! --- symmetry ---
       if( f_selectBlock( tag_symmetry) == 0) then
          if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_symmetry --")')
          if( f_getIntValue( tag_sw_inversion, iret) == 0) inversion_symmetry = iret
          if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** inversion_symmetry = ",i6)') inversion_symmetry
          ! --- symmetry_check_criterion
          if(f_getRealValue(tag_symmetry_check_criterion,dret,'') == 0) symmetry_check_criterion = dret
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** symmetry_check_criterion = ",d20.8)') symmetry_check_criterion
          end if
          ! --- method ---
          call set_symmetry_method(symmetry_method)
          if(f_getRealValue(tag_misalignment,dret,'bohr') == 0) misalignment = dret
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** symmetry_method = ",i6)') symmetry_method
             if(symmetry_method == AUTOMATIC) then
                write(nfout,'(" !** misalignment = ",f15.10)') misalignment 
             end if
          end if

          if(symmetry_method == AUTOMATIC) then
             if(f_getIntValue(tag_sw_non_symmorphic,iret) == 0) sw_non_symmorphic = iret  
             if(ipriinputfile>=1 .and. printable) write(nfout,'(a,i6)') ' !** sw_non_symmorphic = ',sw_non_symmorphic
             if(f_getRealValue(tag_cpumax,dret,'sec') == 0) cpumax_symsearch = dret
             if(ipriinputfile>=1 .and. printable .and. cpumax_symsearch>0) write(nfout,'(a,f10.5,a)') &
               & ' !** cpumax (for symmetry search) = ',sw_non_symmorphic, '(s)'
             if(f_getIntValue(tag_sw_read_symmetry_from_file,iret)==0) sw_read_symmetry_from_file = iret
             if(ipriinputfile>=1 .and. printable .and.  sw_read_symmetry_from_file == OFF &
               & .and. (icond == CONTINUATION .or. icond == FIXED_CHARGE_CONTINUATION &
               & .or. icond == COORDINATE_CONTINUATION)) then
               write(nfout,'(a)') ' !** since sw_read_symmetry_from_file is set to OFF, &
               & symmetry will be re-searched '
             endif
          endif

          ! --- sw_supercell_symmetry ---
          if(f_getIntValue( tag_sw_supercell_symmetry, iret) == 0) sw_supercell_symmetry = iret
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** sw_supercell_symmetry = ",i6)') sw_supercell_symmetry
          end if
             
          if( f_getIntValue( tag_sw_allow_mag_sym_inversion,iret)==0 ) then
             sw_allow_mag_sym_inversion = iret
             write(nfout,*) '!** sw_allow_mag_sym_inversion is set to ', iret
          endif
          if( f_getIntValue( tag_sw_allow_improper_rotation,iret)==0 ) then
             sw_allow_improper_rotation = iret
             write(nfout,*) '!** sw_allow_improper_rotation is set to ', iret
          endif
          if( f_getIntValue( tag_sw_allow_frac_translation,iret)==0 ) then
             sw_allow_frac_translation = iret
             write(nfout,*) '!** sw_allow_frac_translation is set to ', iret
          endif
#if 1
          if ( .not. noncol ) then
             sw_use_magnetic_symmetry = off
             if( f_getIntValue( tag_sw_use_magnetic_symmetry,iret)==0 ) then
                sw_use_magnetic_symmetry = iret
                write(nfout,*) '!** sw_use_magnetic_symmetry is set to ', iret
             endif
             if ( imag==FERRO .and. sw_spinorbit_second_variation==ON ) then
                sw_use_magnetic_symmetry = ON
             endif
          endif
#endif
          ! --- tspace ---
          if( f_selectBlock( tag_tspace) == 0) then
             if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !*  tag_tspace is found")')
             call set_lattice_system()  ! -> il
             number_is_given = f_getIntValue( tag_num_generators, iret) == 0
             if(number_is_given) then
                ngen = iret
                if( ngen >= 4 .or. ngen <= 0) then
                   if(ipriinputfile >= 0 .and. printable) &
                        & write(nfout,'(" !* <ngen> given is not proper, so it is redefined to be 1")')
                   ngen = 1
                end if
             end if
             ! --- generators ---
             prealloc = .true.
             !   --- count generators
             call set_igen_jgen(prealloc,iret)
             if(iret < 0) call phase_error_with_msg(nfout, ' symmetry operation generators are not given properly <<m_CS_rd_n>>'&
                                                         , __LINE__, __FILE__)
             if(.not.number_is_given) ngen = iret
             if(number_is_given .and. ngen > iret) ngen = iret
             if(ngen > 3) call phase_error_with_msg(nfout, ' number of operation generators sould be smaller than 4 <<m_CS_rd_n>>'&
                                                         , __LINE__, __FILE__)
                if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** ngen = ",i5," <<m_CS_rd_n>>")') ngen
             !   --- allocate and set generatos for supercell
             call alloc_igen_jgen()
             prealloc = .false.
             call set_igen_jgen(prealloc,iret)  ! -> igen, jgen
             if(imag==ANTIFERRO) call set_iaf_jaf(iret)  ! -> iaf, jaf

             if(sw_supercell_symmetry == ON) then
                ! --- additional translation operators ---
                prealloc = .true.
                !   --- count additional translation vectors ---
                call set_jgen_tl(prealloc,iret)
                if(iret < 0) call phase_error_with_msg(nfout, &
                             ' symmetry operation generators for sc are not given properly <<m_CS_rd_n>>' &
                             , __LINE__, __FILE__)
                ngen_tl = iret
                if(ipriinputfile>=1.and.printable) write(nfout,'(" !** ngen_tl = ",i5," <<m_CS_rd_n>>")') ngen_tl
                !   --- allocate and set additional translations for supercell
                call alloc_jgen_tl()
                prealloc = .false.
                call set_jgen_tl(prealloc,iret)  ! -> jgen_tl
             end if

! =================================== added by K. Tagami ==================== 12.0A
             if (symmetry_method == AUTOMATIC .and. il>=1 ) then
                gen_name_in_carts = .true.
                if(printable)then
                write(nfout,*) '** ** ** ** ** ** **'
                write(nfout,*) '* Generators are assumed to be in Cartesian system'
                write(nfout,*) '** ** ** ** ** ** **'
                endif
             else
                if ( f_getIntValue( tag_gen_name_in_carts, iret) == 0 ) then
                   if ( iret == YES ) then
                      gen_name_in_carts = .true.
                      if(printable)then
                      write(nfout,*) '** ** ** ** ** ** **'
                      write(nfout,*) '* Generators are assumed to be in Cartesian system'
                      write(nfout,*) '** ** ** ** ** ** **'
                      endif
                   else
                      if(printable)then
                      write(nfout,*) '** ** ** ** ** ** **'
                      write(nfout,*) '* Generators are assumed to be in Bravais system'
                      write(nfout,*) '** ** ** ** ** ** **'
                      endif
                   endif
                endif
             endif
!
             if ( f_getIntValue( tag_use_altv_rltv, iret) == 0 ) then
                if ( iret == YES ) then
                   use_altv_rltv = .true.
                   if(printable)then
                   write(nfout,*) '** use_altv_rltv is set to ', use_altv_rltv
                   write(nfout,*) '** ** ** ** ** ** **'
                   endif
                else
                   if(printable)then
                   write(nfout,*) '** use_altv_rltv is set to .false. '
                   write(nfout,*) '** ** ** ** ** ** **'
                   endif
                endif
             endif
! =========================================================================== 12.0A

             iret = f_selectParentBlock()

          else
             call set_crystal_structure()  ! -> nbztyp_spg
             call m_CS_set_default_symm_op(nfout) ! -> igen, jgen
          end if
          iret = f_selectParentBlock()
       else !!! K.Mae 040121
          call set_crystal_structure()
          call m_CS_set_default_symm_op(nfout) 
       end if

       if(ipriinputfile >= 1 .and. printable) then
          write(nfout,'(" !** il = ",i0, 3x, a)') il, Char_lattice_system(il)
       end if

       ! --- strain
       if( f_selectBlock( tag_strain) == 0) then
          call set_strain_tensor()
          if(ipriinputfile >= 1 .and. printable) then
             write(nfout,'(" !** Strain tensor <<m_CS_rd_n>>")')
             write(nfout,'(" !**      [",3(1x,f10.5)," ]")') strain(1,1:3)
             write(nfout,'(" !**  e = [",3(1x,f10.5)," ]")') strain(2,1:3)
             write(nfout,'(" !**      [",3(1x,f10.5)," ]")') strain(3,1:3)
          end if
          if( f_getIntValue( tag_sw_strained_cell, iret) == 0) sw_strained_cell = iret
          iret = f_selectParentBlock()
       end if

       ! supercell
       if( f_selectBlock( tag_supercell) == 0) then
          if( f_getIntValue( tag_sw_supercell, iret) == 0) sw_supercell = iret
          ! --- unit_cell ---
          iret = f_getStringValue(tag_unit_cell_type,rstr, LOWER)
          if( rstr == val_primitive ) then
             supercell_unit_type = PRIMITIVE ! 1
          else
             supercell_unit_type = BRAVAIS   ! 0
          end if
          if( f_selectBlock( tag_size) == 0) then
             if( f_getIntValue( tag_n1, iret) == 0) n1_sc = iret
             if( f_getIntValue( tag_n2, iret) == 0) n2_sc = iret
             if( f_getIntValue( tag_n3, iret) == 0) n3_sc = iret
             iret = f_selectParentBlock()
          end if
          iret = f_selectParentBlock()
       end if

       ! --- unit_cell ---
       iret = f_getStringValue(tag_unit_cell_type,rstr, LOWER)
       if( rstr == val_primitive ) then
          iuctype = PRIMITIVE ! 1
       else
          iuctype = BRAVAIS   ! 0
       end if
! ================================ Added by K. Tagami ===============
       kt_iuctype = iuctype
! ==================================================================
       if(ipriinputfile >= 1 .and. printable) then
          if(iret == 0) then
             write(nfout,'(" !** unit_cell_type = ",i6," [",a,"] <<m_CS_rd_n>>")') iuctype, trim(rstr)
          else
             write(nfout,'(" !** unit_cell_type is not given")')
          end if
       end if
       ! peek lattice_system
       hext = .false.
       if(f_selectBlock(tag_symmetry)==0) then
          if(f_selectBlock(tag_tspace)==0) then
            if(f_getStringValue('lattice_system',rstr,LOWER)==0) then
              if(rstr==tag_hexagonal .or. rstr==tag_hcp .or. &
              &  rstr==tag_trigonal  .or. rstr==tag_rhombohedral) hext = .true.
            endif
            iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       endif
       if(f_getRealValue(tag_unit_cell_scale,dret,'bohr')==0) then
          unit_cell_scale = dret
          if(printable) write(nfout,'(a,f10.5)') ' !** unit_cell_scale : ',unit_cell_scale
       endif
       if( f_selectBlock( tag_unit_cell) == 0) then
          allocate(avec(1:3)); allocate(bvec(1:3)); allocate(cvec(1:3))
          iret = f_readUnitCell(avec,bvec,cvec,a,b,c,ca,cb,cc,ucinptype,'Bohr' &
               & ,tag_avector,tag_bvector,tag_cvector,tag_a,tag_b,tag_c &
               & ,tag_alpha,tag_beta,tag_gamma,unit_cell_scale,hext,printable)
          if(ipriinputfile >= 1 .and. printable) then
             if(iret == 0) then
                write(nfout,'(" !** a,b,c,ca,cb,cc = ",6f12.6," <<m_CS_rd_n>>")') a,b,c,ca,cb,cc
                write(nfout,'(" !** avec = ",3f12.6)') avec(1:3)
                write(nfout,'(" !** bvec = ",3f12.6)') bvec(1:3)
                write(nfout,'(" !** cvec = ",3f12.6)') cvec(1:3)
                write(nfout,'(" !** ucinptype = ",i6)') ucinptype
             else
                write(nfout,'(" !** iret /= 0: iret = ",i6)') iret
             end if
          end if
!!!!
          call m_CS_gnrt_tmatrices(il)

          if ( sw_strained_cell == OFF ) then
             if(is_hexagonal(a,b,ca,cb,cc).and. il==1 .and. symmetry_method == AUTOMATIC) then
                if(printable) then
                   write(nfout,'(" !** lattice_system is converted to hexagonal")')
                endif
                il = 0
             endif
          endif

          if (iuctype == BRAVAIS) then ! Bravais: 0
!             if (ucinptype == 0) then ! length: 0
                call bravais2primitive(nfout,b2pmat,a,b,c,ca,cb,cc,avec,bvec,cvec,il) ! in b_CS

                if ( icond /= COORDINATE_CONTINUATION ) then
                   if(ipriinputfile >= 1 .and. printable) then
                      write(nfout,*) '!** primitive unit cell vectors are given as'
                      write(nfout,'(" !** avec = ",3f12.6)') avec(1:3)
                      write(nfout,'(" !** bvec = ",3f12.6)') bvec(1:3)
                      write(nfout,'(" !** cvec = ",3f12.6)') cvec(1:3)
                   end if
                endif

             else                ! Primitive
                if  (ucinptype == 0) then ! length: 0
                   !if(printable) write(nfout,*) 'You must give a lattice system by the Vector expression for the primitive cell.'
                   call phase_error_with_msg(nfout, &
                   'You must give a lattice system by the Vector expression for the primitive cell.', &
                   __LINE__, __FILE__) 
                end if
                call primitive2bravais(nfout,p2bmat,avec,bvec,cvec,a,b,c,ca,cb,cc,il) ! in b_CS
!             end if
          end if
          if(sw_strained_cell == ON) then
             call set_strained_cell(avec,bvec,cvec)
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** << Strained cell >>")')
                write(nfout,'(" !** avec = ",3f12.6)') avec(1:3)
                write(nfout,'(" !** bvec = ",3f12.6)') bvec(1:3)
                write(nfout,'(" !** cvec = ",3f12.6)') cvec(1:3)
             end if
             call primitive2bravais(nfout,p2bmat,avec,bvec,cvec,a,b,c,ca,cb,cc,il) ! in b_CS
          end if
          if(sw_supercell == ON) then
             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** << supercell >>")')
                write(nfout,'(" !** unit_cell_type =",i5)') supercell_unit_type
                write(nfout,'(" !** n1,n2,n3 = ",3i5)') n1_sc,n2_sc,n3_sc
             end if
          end if
!!!!
          altv(1:3,1) = avec(1:3); altv(1:3,2) = bvec(1:3); altv(1:3,3) = cvec(1:3)
          deallocate(avec); deallocate(bvec); deallocate(cvec)
          call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS
#ifdef ENABLE_ESM_PACK
          if(sw_esm==ON .and. esm_z1_defined)then
             esm_w = esm_z1-altv(3,3)*0.5d0
          endif
#endif
          iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_unit_cell is not found")')
       end if

! === KT 2018/03/07
       if ( symmetry_method == AUTOMATIC .and. il == 1 ) then
          if ( ucinptype == 1 ) then
             use_altv_rltv = .true.
             write(nfout,*) '** use_altv_rltv is turned to ', use_altv_rltv
          endif
       endif
! === KT 2018/03/07

!  ---- refcell for band-unfolding---- 2016/04/14AS
       if( f_selectBlock( tag_reference_cell) == 0) then
          iret = f_getStringValue(tag_unit_cell_type,rstr, LOWER)
          if( rstr == val_bravais ) then
             refcell_unit_type = BRAVAIS   ! 0
          else
             refcell_unit_type = PRIMITIVE ! 1
          end if

          allocate(avec(1:3)); allocate(bvec(1:3)); allocate(cvec(1:3))

!          iret = f_readUnitCell(avec,bvec,cvec,a,b,c,ca,cb,cc,ucinptype,'Bohr' &
!               & ,tag_avector,tag_bvector,tag_cvector,tag_a,tag_b,tag_c &
!               & ,tag_alpha,tag_beta,tag_gamma,printable)
          iret = f_readUnitCell(avec,bvec,cvec,a,b,c,ca,cb,cc,ucinptype,'Bohr' &
               & ,tag_avector,tag_bvector,tag_cvector,tag_a,tag_b,tag_c &
               & ,tag_alpha,tag_beta,tag_gamma,unit_cell_scale,hext,printable)

          if ( refcell_unit_type == BRAVAIS ) then
             tmpmat = 0.0d0
             tmpmat(1,1) = 1.0d0;  tmpmat(2,2) = 1.0d0;  tmpmat(3,3) = 1.0d0
             call bravais2primitive( nfout, tmpmat, a, b, c, ca, cb, cc, &
                  &                  avec, bvec, cvec, 1 )
          endif
          altv_refcell(1:3,1) = avec(1:3);
          altv_refcell(1:3,2) = bvec(1:3);
          altv_refcell(1:3,3) = cvec(1:3)

          deallocate(avec); deallocate(bvec); deallocate(cvec)
          call altv_2_rltv(altv_refcell,rltv_refcell,univol_refcell,rvol_refcell)

          iret = f_selectParentBlock()
       else
          altv_refcell = altv;   rltv_refcell = rltv;
          univol_refcell = univol;  rvol_refcell = rvol
       endif
! --------------------------------- 2016/04/14AS

       iret = f_selectParentBlock()
    else
       if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_structure is not found")')
    end if

    if ( icond /= COORDINATE_CONTINUATION ) then
       if(ipriinputfile >= 1 .and. printable) then
          write(nfout,'(" !** univol = ",f20.12," <<m_CS_rd_n>>")') univol
       end if
    endif

  contains
    subroutine set_symmetry_method(symmetry_method)
      integer, intent(inout) :: symmetry_method
      integer :: iret
      logical :: tf
      iret = f_getStringValue( tag_method, rstr,LOWER)
      if(iret == 0) then
         if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !*  tag_method is found")')
         call strncmp0(trim(rstr),tag_automatic,tf)
         if(tf) then
            symmetry_method = AUTOMATIC
            goto 1001
         end if
         call strncmp0(trim(rstr),tag_manual,tf)
         if(tf) then
            symmetry_method = MANUAL
            goto 1001
         end if
      end if
1001  continue
    end subroutine set_symmetry_method

    subroutine set_crystal_structure()
      integer :: iret
      logical :: tf
      nbztyp_spg = 0
      iret = f_getStringValue( tag_crystal_structure, rstr,LOWER)
      if(iret /= 0) iret = f_getStringValue( tag_crystal, rstr,LOWER)
      if(iret == 0) then
         if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !*  tag_crystal is found")')
         call strncmp0(trim(rstr),tag_trigonal,tf)
         if(tf) then
            nbztyp_spg = TRIGONAL
            goto 1001
         end if
         call strncmp0(trim(rstr),tag_hexagonal,tf)
         if(tf) then
            nbztyp_spg = HEXAGONAL
            goto 1001
         end if
         call strncmp0(trim(rstr),tag_hcp,tf)
         if(tf) then
            nbztyp_spg = HCP
            goto 1001
         end if
         call strncmp0(trim(rstr),tag_simple_cubic,tf)
         if(tf) then
            nbztyp_spg = SIMPLE_CUBIC
            goto 1001
         end if
         call strncmp0(trim(rstr),tag_diamond, tf)
         if(tf) then
            nbztyp_spg = DIAMOND
            goto 1001
         end if

         call strncmp0(trim(rstr),tag_facecentered, tf)
         if(.not.tf) call strncmp0(trim(rstr),tag_fcc, tf)
         if(tf) then
            nbztyp_spg = FCC
            goto 1001
         end if

         call strncmp0(trim(rstr),tag_bodycentered,  tf)
         if(.not.tf) call strncmp0(trim(rstr),tag_bcc, tf)
         if(tf) then
            nbztyp_spg = BCC
            goto 1001
         end if
      end if
1001  continue
    end subroutine set_crystal_structure
         
    subroutine set_imag(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_para, len(tag_para), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_nonmagnetic,len(tag_nonmagnetic),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_nonmag, len(tag_nonmag),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_none,   len(tag_nonmag),tf)
!!$      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '1', 1, tf)
      if(tf) then
         imag = NONMAG
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_af, len(tag_af), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_antiferro ,len(tag_antiferro), tf)
!!$      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '-2', 2, tf)
      if(tf) then
         imag = ANTIFERRO
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_ferro, len(tag_ferro), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_magnetic, len(tag_magnetic),tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_mag, len(tag_mag),tf)
!!$      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '2', 1, tf)
      if(tf) then
         imag = FERRO
         goto 1001
      end if

! =========================== added by K. Tagami ================= 11.0
      call strncmp2(rstr, FMAXVALLEN, tag_noncollinear, len(tag_noncollinear), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, tag_noncol, len(tag_noncol), tf)
!!$      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, '4', 1, tf)
      if(tf) then
         imag = NONCOLLINEAR
         goto 1001
      end if
! ================================================================ 11.0

1001  continue
      if(.not.(imag == NONMAG .or. imag == ANTIFERRO .or. imag == FERRO)) then
         if(printable) then
            write(nfout,'("!ERROR A setting in the inputfile as (( magnetic_state = ",a10," )) is not proper.")') rstr
            write(nfout,'("!ERROR  The value must be nonmag (=para), magnetic (=ferro), or antiferro for phae/0 3D version.")')
         end if
         call phase_error_with_msg(nfout,'(!** nspin is not given properly in the inputfile' &
              & ,__LINE__,__FILE__)
      end if
    end subroutine set_imag

    subroutine set_imag_from_nspin(iret)
      ! Coded by T. Yamasaki, 16th June 2023
      integer,intent(in) :: iret
      logical :: tf
      if(iret==1) then
         imag = NONMAG
      else if(iret==2) then
         imag = FERRO
      else if(iret==-2) then
         imag = ANTIFERRO
      else if(iret==4) then
         imag = NONCOLLINEAR
      else
      end if
      if(.not.(iret==1 .or. iret==2 .or. iret==-2)) then
         if(printable) then
            write(nfout,'("!ERROR A setting in the inputfile as (( nspin = ", i4, " )) is not proper.")') iret
            write(nfout,'("!ERROR  The value must be 1, 2, or -2 for phase/0 3D version")')
         end if
         call phase_error_with_msg(nfout,'(!** nspin is not given properly in the inputfile' &
              & ,__LINE__,__FILE__)
      end if
    end subroutine set_imag_from_nspin

    subroutine set_spin_fix_period(rstr)
      ! Coded by T. Yamasaki, 06th Aug. 2009
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_whole, len(tag_whole), tf)
      if(tf) then
         spin_fix_period = WHOLE
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_initially, len(tag_initially), tf)
      if(tf) then
         spin_fix_period = INITIALLY
         goto 1001
      end if
1001  continue
    end subroutine set_spin_fix_period

! =========================== Added by K. Tagami =============== 11.0
    subroutine set_sw_magnetic_symmetry

      sw_use_magnetic_symmetry = OFF
      if ( sw_reduce_sym_by_magmom == ON .or. sw_reduce_sym_by_orbital == ON ) then
         sw_use_magnetic_symmetry = ON
      endif

    end subroutine set_sw_magnetic_symmetry

    subroutine set_sw_reduce_sym_by_magmom
      logical :: tf

      tf = ( f_getIntValue( tag_reduce_sym_by_magmom, iret ) == 0 )
      if (tf) then
         if ( iret == ON ) then
            sw_reduce_sym_by_magmom = ON
            write(nfout,*) '**** sw_reduce_sym_by_magmom is set ON'
         else
            sw_reduce_sym_by_magmom = OFF
            write(nfout,*) '**** sw_reduce_sym_by_magmom is set OFF'
         endif
      else
         write(nfout,*) '**** sw_reduce_sym_by_magmom is set to default, ON'
      endif

    end subroutine set_sw_reduce_sym_by_magmom

    subroutine set_sw_reduce_sym_by_orbital
      logical :: tf

      tf = ( f_getIntValue( tag_reduce_sym_by_orbital, iret ) == 0 )
      if (tf) then
         if ( iret == ON ) then
            sw_reduce_sym_by_orbital = ON
            write(nfout,*) '**** sw_reduce_sym_by_orbital is set ON'
         else
            sw_reduce_sym_by_orbital = OFF
            write(nfout,*) '**** sw_reduce_sym_by_orbital is set OFF'
         endif
      else
         write(nfout,*) '**** sw_reduce_sym_by_orbital is set to default, OFF'
      endif

    end subroutine set_sw_reduce_sym_by_orbital

! ==== KT_add === 2014/08/09 & 13.2S
    subroutine set_magnetic_restriction
      if ( (.not. noncol) .and. (sw_use_magnetic_symmetry == ON) ) then
         sw_fix_global_quantz_axis = on

         if ( f_selectBlock( tag_magnetic_restriction ) == 0 ) then
            call set_global_quantz_axis
            iret = f_selectParentBlock()
         endif
      else
         if ( f_selectBlock( tag_magnetic_restriction ) == 0 ) then
            if ( f_getIntValue( tag_sw_fix_global_quantz_axis, iret ) == 0 ) then
               sw_fix_global_quantz_axis = iret
            endif
            if ( sw_fix_global_quantz_axis == ON ) call set_global_quantz_axis
            iret = f_selectParentBlock()
         endif
      endif
    end subroutine set_magnetic_restriction

    subroutine set_global_quantz_axis
      integer :: Flag
      Real(kind=DP) :: mdx, mdy, mdz, cnorm
      Real(kind=DP) :: ctmp, theta, phi

      Flag = 0
      mdx = 0.0d0; mdy = 0.0d0; mdz = 0.0d0
      theta = 0.0d0;   phi = 0.0d0

      Global_Quantz_Axis_Fixed(1) = 0.0d0
      Global_Quantz_Axis_Fixed(2) = 0.0d0
      Global_Quantz_Axis_Fixed(3) = 1.0d0

      if ( f_selectBlock( tag_global_quantz_axis ) == 0 ) then
         if( f_getRealValue( tag_sx, dret, '') == 0 ) then
            mdx = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_sy, dret, '') == 0 ) then
            mdy = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_sz, dret, '') == 0 ) then
            mdz = dret;   Flag = 1
         endif
         cnorm = sqrt( mdx**2 + mdy**2 + mdz**2 )

         if ( abs(cnorm) > 1.0E-4 ) then
            Global_Quantz_Axis_Fixed(1) = mdx / cnorm;
            Global_Quantz_Axis_Fixed(2) = mdy / cnorm;
            Global_Quantz_Axis_Fixed(3) = mdz / cnorm;
         endif
         !
         if ( Flag == 0 ) then
            if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
               theta = dret
            end if
            if ( f_getRealValue( tag_phi, dret, "" ) == 0 ) then
               phi = dret
            end if

            theta  = theta / 180.0d0 *PAI
            phi  =  phi / 180.0d0 *PAI
!
            Global_Quantz_Axis_Fixed(1) = sin( theta ) *cos( phi )
            Global_Quantz_Axis_Fixed(2) = sin( theta ) *sin( phi )
            Global_Quantz_Axis_Fixed(3) = cos( theta )
         endif

         iret = f_selectParentBlock()
      endif
!
      write(nfout,*) '!************* Magnetic Restriction ************'
      write(nfout,'(A,F10.5,A,F10.5,A,F10.5)' ) &
           &             '!** global axis : sx = ', Global_Quantz_Axis_Fixed(1), &
           &             ' sy = ', Global_Quantz_Axis_Fixed(2), &
           &             ' sz = ', Global_Quantz_Axis_Fixed(3)

    end subroutine set_global_quantz_axis
! ============== 2014/08/09 & 13.2S

    subroutine set_sw_constrain_on_grad_corr
      logical :: tf

      tf = ( f_getIntValue( tag_constrain_on_grad_corr, iret ) == 0 )
      if (tf) then
         if ( iret == ON ) then
            sw_constrain_on_grad_correction = ON
            write(nfout,*) '**** sw_constrain_on_grad_correction is set ON'
         else
            sw_constrain_on_grad_correction = OFF
            write(nfout,*) '**** sw_constrain_on_grad_correction is set OFF'
         endif
      else
         write(nfout,*) '**** sw_constrain_on_grad_correction is set to default, OFF'
      endif

! === KT_add == 2014/08/09
      if ( sw_fix_global_quantz_axis == ON ) then
         sw_constrain_on_grad_correction = ON
         write(nfout,*) '**** sw_constrain_on_grad_correction is turned ON'
      endif
! ============== 2014/08/09

    end subroutine set_sw_constrain_on_grad_corr

    subroutine set_level_projection_paw_charge
      logical :: tf

      tf = ( f_getIntValue( tag_level_projection_paw_charge, iret ) == 0 )

      if (tf) then
         level_of_projection_paw_charge = iret

         if ( iret > 0 .and. iret < 4 ) then
            write(nfout,*) '**** level_of_projection_paw_charge is set to ', iret
         else
            level_of_projection_paw_charge = 1
            write(nfout,*) '**** level_of_projection_paw_charge is set to ', 1
         endif

      else
         write(nfout,*) '**** level_of_projection_paw_charge is set to default', 1
      endif

    end subroutine set_level_projection_paw_charge
! =============================================================== 11.0

! ==== KT_add ==== 2014/08/04
    subroutine set_threshold_lowering_contami        ! in case of GGA+noncl
      if ( f_getIntValue( tag_sw_neglect_low_vorticity, iret ) == 0 ) then
         sw_neglect_low_vorticity = iret
      endif
      if ( f_getIntValue( tag_sw_neglect_low_helicity, iret ) == 0 ) then
         sw_neglect_low_helicity = iret
      endif

      if ( sw_neglect_low_vorticity == ON .and. sw_neglect_low_helicity == ON ) then
         call phase_error_with_msg(nfout, "Select either sw_neglect_low_vorticity or sw_neglect_low_helicty" &
                                        , __LINE__, __FILE__)
      endif
!
      if ( sw_neglect_low_vorticity == ON ) then
         sw_monitor_magnetic_vorticity = ON
         if ( f_getRealValue( tag_threshold_vorticity, dret, '' ) == 0 ) then
            threshold_vorticity = dret
         endif
         write(nfout,*) '!** sw_neglect_low_vorticity is ', sw_neglect_low_vorticity
         write(nfout,*) '!** threshold_vorticity is ', threshold_vorticity
      endif
      if ( sw_neglect_low_helicity == ON ) then
         sw_monitor_magnetic_vorticity = ON
         if ( f_getRealValue( tag_threshold_helicity, dret, '' ) == 0 ) then
            threshold_helicity = dret
         endif
         write(nfout,*) '!** sw_neglect_low_helicity is ', sw_neglect_low_helicity
         write(nfout,*) '!** threshold_helicity is ', threshold_helicity
      endif
!
    end subroutine set_threshold_lowering_contami
! ================ 2014/08/04

! ==== KT_add === 2014/09/26
    subroutine set_sw_neglect_magmom
      if ( f_getIntValue( tag_sw_neglect_magmom, iret ) == 0 ) then
         sw_neglect_magmom = iret
      endif
      write(nfout,*) '!** sw_neglect_magmom is ', sw_neglect_magmom
    end subroutine set_sw_neglect_magmom
! =============== 2014/09/26

! ================================= KT_add =================== 11.0&13.0U
    subroutine set_magnetic_constraint
      logical :: tf
!
      if ( f_getStringValue( tag_constraint_type,rstr,LOWER ) == 0 ) then

         call strncmp0( tag_moment_vals_global, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_VALS_GLOBAL

         call strncmp0( tag_moment_direc_global, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_DIREC_GLOBAL

         call strncmp0( tag_moment_vals_local, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_VALS_LOCAL

         call strncmp0( tag_moment_direc_local, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_DIREC_LOCAL

         call strncmp0( tag_moment_vals_occmat, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_VALS_OCCMAT

         call strncmp0( tag_moment_direc_hardpart, trim(rstr), tf )
         if ( tf ) mag_constraint_type = MAG_MOMENT_DIREC_HARDPART

      endif

      write(nfout,*) '!**************************************************'

      if ( mag_constraint_type == MAG_MOMENT_VALS_GLOBAL ) then
         call set_mag_constraint_type1;    sw_magnetic_constraint = ON
         write(nfout,*) '!*** Magnetic constraint type : MOMENT_VALS_GLOBAL'

      else if ( mag_constraint_type == MAG_MOMENT_DIREC_GLOBAL ) then
         call set_mag_constraint_type2;    sw_magnetic_constraint = ON
         write(nfout,*) '!*** Magnetic constraint type : MOMENT_DIREC_GLOBAL'

      else if ( mag_constraint_type == MAG_MOMENT_VALS_LOCAL ) then
         sw_magnetic_constraint = ON
         write(nfout,*) '!*** Magnetic constraint type : MOMENT_VALS_LOCAL'

      else if ( mag_constraint_type == MAG_MOMENT_DIREC_LOCAL ) then
         sw_magnetic_constraint = ON
         write(nfout,*) '!*** Magnetic constraint type : MOMENT_DIREC_LOCAL'

      else if ( mag_constraint_type == MAG_MOMENT_DIREC_HARDPART ) then
         sw_magnetic_constraint = ON
         write(nfout,*) '!*** Magnetic constraint type : MOMENT_DIREC_HARDPART'

      else if ( mag_constraint_type == MAG_MOMENT_VALS_OCCMAT ) then
         if ( num_projectors > 0 ) then
            sw_magnetic_constraint = ON
            write(nfout,*) '!*** Magnetic constraint type : MOMENT_VALS_OCCMAT'
         else
            write(nfout,*) '!!!! ---- This constraint requires projectors'
         endif

      endif

! ------------
      if ( sw_magnetic_constraint == OFF ) then
         write(nfout,*) '! mag_constraint is not set '; return
      endif
! ------------

      if ( f_getRealValue( tag_lambda,dret,"hartree" ) == 0 ) then
         mag_constraint_lambda = dret
      end if

      if ( f_getStringValue( tag_damping_method, rstr,LOWER ) == 0 ) then
         call strncmp0( tag_damp_abrupt, trim(rstr), tf )
         if ( tf ) damping_method_mag_constraint = ABRUPT

         call strncmp0( tag_damp_stepwize, trim(rstr), tf )
         if ( tf ) damping_method_mag_constraint = STEPWISE

         call strncmp0( tag_damp_linear, trim(rstr), tf )
         if ( tf ) damping_method_mag_constraint = LINEAR
      endif

      if ( damping_method_mag_constraint == STEPWISE ) then
         if ( f_getIntValue( tag_num_intermid_lambda, iret ) == 0 ) then
            if ( iret < 0 ) num_intermid_lambda = 0
            if ( iret > nmax_intermid_lambda ) num_intermid_lambda = nmax_intermid_lambda
            num_intermid_lambda = iret
         endif

         if ( f_getRealValue( tag_edelta_change_lambda_first,dret,"hartree" )==0 ) then
            edelta_change_lambda_first = dret
         endif
         if ( f_getRealValue( tag_edelta_change_lambda_last,dret,"hartree" )==0 ) then
            edelta_change_lambda_last = dret
         endif
         if ( num_intermid_lambda == 0 ) then
            edelta_change_lambda_last = edelta_change_lambda_first
         endif
      endif

      if ( damping_method_mag_constraint == ABRUPT .or. &
           &  damping_method_mag_constraint == LINEAR ) then

         if ( f_getIntValue( tag_max_iter_elec_mag_constr, iret ) == 0 ) then
            if ( iret < 0 ) max_iter_elec_mag_constraint = 0
            max_iter_elec_mag_constraint = iret
         endif
      endif

      write(nfout,*) '! mag_constraint_lambda is set to ', &
           &            mag_constraint_lambda
      write(nfout,*) '! damping_method is set to ', &
           &            damping_method_mag_constraint

      if ( damping_method_mag_constraint == STEPWISE ) then
         write(nfout,*) '! number_of_intermidiate lambda is ', num_intermid_lambda
         write(nfout,*) '! edelta_change_lambda_first is ', edelta_change_lambda_first
         write(nfout,*) '! edelta_change_lambda_last  is ', edelta_change_lambda_last
      endif

      if ( damping_method_mag_constraint == ABRUPT .or. &
           &  damping_method_mag_constraint == LINEAR ) then
         write(nfout,*) '! max_iter_elec_mag_constraint is ', &
              &            max_iter_elec_mag_constraint
      endif

      if ( f_getIntValue( tag_max_iter_ion_mag_constr, iret ) == 0 ) then
         if ( iret < 0 ) max_iter_ion_mag_constraint = 0
         max_iter_ion_mag_constraint = iret
      endif
      if ( f_getIntValue( tag_max_iter_cell_mag_constr, iret ) == 0 ) then
         if ( iret < 0 ) max_iter_cell_mag_constraint = 0
         max_iter_cell_mag_constraint = iret
      endif

      if ( f_getIntValue( tag_sw_fix_charge_after_constr, iret ) == 0 ) then
         sw_fix_charge_after_constraint = iret
         write(nfout,*) '! sw_fix_charge_after_constraint is ', &
              &             sw_fix_charge_after_constraint
      endif

      write(nfout,*) '! max_iter_ion_mag_constraint is ', max_iter_ion_mag_constraint
      write(nfout,*) '! max_iter_cell_mag_constraint is ', max_iter_cell_mag_constraint

      write(nfout,*) '! *********************************************** '

    end subroutine set_magnetic_constraint

!--
    subroutine set_mag_constraint_type1
      integer :: Flag
      Real(kind=DP) :: mx, my, mz, cnorm
      Real(kind=DP) :: ctmp, theta, phi

      Flag = 0
      mx = 0.0d0; my = 0.0d0; mz = 0.0d0
      theta = 0.0d0;   phi = 0.0d0

      mag_moment0_global = 0.0d0

      if ( f_selectBlock( tag_magnetic_moment ) == 0 ) then

         if ( noncol ) then
            if( f_getRealValue( tag_mx, dret, '') == 0 ) then
               mx = dret;   Flag = 1
            endif
            if( f_getRealValue( tag_my, dret, '') == 0 ) then
               my = dret;   Flag = 1
            endif
            if( f_getRealValue( tag_mz, dret, '') == 0 ) then
               mz = dret;   Flag = 1
            endif

            mag_moment0_global(1) = mx
            mag_moment0_global(2) = my
            mag_moment0_global(3) = mz
            norm = sqrt( mx**2 + my**2 +mz**2 )
         !
            if ( Flag == 0 ) then
               if ( f_getRealValue( tag_norm, dret, "" ) == 0 ) then
                  norm = dret
               end if
               if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
                  norm = dret
               end if
               if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
                  theta = dret
               end if
               if ( f_getRealValue( tag_phi, dret, "" ) == 0 ) then
                  phi = dret
               end if

               theta  = theta / 180.0d0 *PAI
               phi  =  phi / 180.0d0 *PAI
!
               mag_moment0_global(1) = norm *sin( theta ) *cos( phi )
               mag_moment0_global(2) = norm *sin( theta ) *sin( phi )
               mag_moment0_global(3) = norm *cos( theta )
            endif

         else
            if ( f_getRealValue( tag_norm, dret, "" ) == 0 ) then
               norm = dret;   mag_moment0_global(1) = norm
            end if
            if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
               norm = dret;   mag_moment0_global(1) = norm
            end if

         endif

         iret = f_selectParentBlock()
      endif
!
      write(nfout,*) '!************* Magnetic constraint : type 1 ************'
      write(nfout,'(A,F10.5)') '!** magnetic moment is wanted to reach ', norm

      if ( noncol ) then
         write(nfout,'(A,F10.5,A,F10.5,A,F10.5)' ) &
              &             '!** The xyz component is mx =', mag_moment0_global(1), &
              &               ' my = ', mag_moment0_global(2), &
              &               ' mz = ', mag_moment0_global(3)

      else
         call m_CS_set_total_spin( nfout, mag_moment0_global(1) )
      endif

    end subroutine set_mag_constraint_type1

    subroutine set_mag_constraint_type2
      integer :: Flag
      Real(kind=DP) :: mdx, mdy, mdz, cnorm
      Real(kind=DP) :: ctmp, theta, phi

      Flag = 0
      mdx = 0.0d0; mdy = 0.0d0; mdz = 0.0d0
      theta = 0.0d0;   phi = 0.0d0

      mag_direc0_global(1) = 0.0d0
      mag_direc0_global(2) = 0.0d0
      mag_direc0_global(3) = 1.0d0

      if ( f_selectBlock( tag_direction ) == 0 .or. &
           &  f_selectBlock( tag_magnetic_moment ) == 0 ) then

         if( f_getRealValue( tag_mdx, dret, '') == 0 ) then
            mdx = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_mdy, dret, '') == 0 ) then
            mdy = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_mdz, dret, '') == 0 ) then
            mdz = dret;   Flag = 1
         endif

         if( f_getRealValue( tag_mx, dret, '') == 0 ) then
            mdx = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_my, dret, '') == 0 ) then
            mdy = dret;   Flag = 1
         endif
         if( f_getRealValue( tag_mz, dret, '') == 0 ) then
            mdz = dret;   Flag = 1
         endif
         cnorm = sqrt( mdx**2 + mdy**2 + mdz**2 )

         if ( abs(cnorm) > 1.0E-4 ) then
            mag_direc0_global(1) = mdx / cnorm;
            mag_direc0_global(2) = mdy / cnorm;
            mag_direc0_global(3) = mdz / cnorm;
         endif
         !
         if ( Flag == 0 ) then
            if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
               theta = dret
            end if
            if ( f_getRealValue( tag_phi, dret, "" ) == 0 ) then
               phi = dret
            end if

            theta  = theta / 180.0d0 *PAI
            phi  =  phi / 180.0d0 *PAI
!
            mag_direc0_global(1) = sin( theta ) *cos( phi )
            mag_direc0_global(2) = sin( theta ) *sin( phi )
            mag_direc0_global(3) = cos( theta )
         endif

         iret = f_selectParentBlock()
      endif
!
      write(nfout,*) '!************* Magnetic constraint : type 2 ************'
      write(nfout,'(A,F10.5,A,F10.5,A,F10.5)' ) &
           &             '!** The xyz component is mdx =', mag_direc0_global(1), &
           &               ' mdy = ', mag_direc0_global(2), &
           &               ' mdz = ', mag_direc0_global(3)

    end subroutine set_mag_constraint_type2
! ============================================================ 11.0&13.0U

    subroutine set_igen_jgen(prealloc,iret)
      logical, intent(in) ::  prealloc
      integer, intent(out) :: iret
      integer :: ip, linecounter
      integer, parameter :: limit_number_of_generators = 50
!!$      integer :: f_selectParentBlock, f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_readSymmetryGenerator
      integer :: txu,txd,tyu,tyd,tzu,tzd
      if( f_selectBlock( tag_generators) == 0) then
         ip = 1
         linecounter = 0
         Generators: do while(.true.)
            linecounter = linecounter+1
            if(ip == 1) then
               if( f_selectFirstTableLine() /= 0) then
                  exit Generators
               end if
            else
               if( f_selectNextTableLine() /= 0) then
                  exit Generators
               end if
            end if
            if( f_readSymmetryGenerator(il, rint ,txu,txd,tyu,tyd,tzu,tzd, tag_rotation, tag_tx,tag_ty,tag_tz,printable) == 0) then
               if(.not.prealloc) then
                  if(ip > ngen) exit
                  igen(ip) = rint
                  jgen(1,1,ip) = txu; jgen(2,1,ip) = txd
                  jgen(1,2,ip) = tyu; jgen(2,2,ip) = tyd
                  jgen(1,3,ip) = tzu; jgen(2,3,ip) = tzd
               
                  if(ipriinputfile >= 1 .and. printable) then
                     write(nfout,'(" !** igen = ",i6,"  jgen = ",6i5)') igen(ip), jgen(1,1,ip), jgen(2,1,ip) &
                          & , jgen(1,2,ip),jgen(2,2,ip),jgen(1,3,ip),jgen(2,3,ip)
                  end if
               end if
               ip = ip + 1
            end if
            if(linecounter > limit_number_of_generators) then
               write(nfout,'(" linecounter exceeds limit_number_of_generatos: ip = ",i8," at <<set_igen_jgen>>")') ip
               write(nfout,'(" check tag names in <generatos block> of the input file")')
               call phase_error_with_msg(nfout, ' linecounter exceeds limit_number_of_generators at <<set_igen_jgen>>'&
                                              , __LINE__, __FILE__)
            end if
         end do Generators
         iret = f_selectParentBlock()
         iret = ip - 1
      else
         iret = 1
         if(.not.prealloc) call set_igen_jgen_default(1)
      end if
    end subroutine set_igen_jgen

    subroutine set_jgen_tl(prealloc,iret)
      logical, intent(in) ::  prealloc
      integer, intent(out) :: iret
      integer :: ip, linecounter
      integer, parameter :: limit_number_of_generators = 50
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_readSymmetryGenerator, f_readTranslationVectors
      integer :: txu,txd,tyu,tyd,tzu,tzd
      if( f_selectBlock( tag_additional_tl) == 0) then
         ip = 1
         linecounter = 0
         Generators: do while(.true.)
            linecounter = linecounter+1
            if(ip == 1) then
               if( f_selectFirstTableLine() /= 0) then
                  exit Generators
               end if
            else
               if( f_selectNextTableLine() /= 0) then
                  exit Generators
               end if
            end if
            if( f_readTranslationVectors(il, txu,txd,tyu,tyd,tzu,tzd, tag_tx,tag_ty,tag_tz,printable) == 0) then
               if(.not.prealloc) then
                  if(ip > ngen_tl) exit
                  jgen_tl(1,1,ip) = txu; jgen_tl(2,1,ip) = txd
                  jgen_tl(1,2,ip) = tyu; jgen_tl(2,2,ip) = tyd
                  jgen_tl(1,3,ip) = tzu; jgen_tl(2,3,ip) = tzd
               
                  if(ipriinputfile >= 1 .and. printable) then
                     write(nfout,'(" !** jgen_tl = ",6i5)') &
                          &  jgen_tl(1:2,1,ip), jgen_tl(1:2,2,ip),jgen_tl(1:2,3,ip)
                  end if
               end if
               ip = ip + 1
            end if
            if(linecounter > limit_number_of_generators) then
               write(nfout,'(" linecounter exceeds limit_number_of_generatos: ip = ",i8," at <<set_igen_jgen>>")') ip
               write(nfout,'(" check tag names in <generatos block> of the input file")')
               call phase_error_with_msg(nfout, ' linecounter exceeds limit_number_of_generators at <<set_igen_jgen>>' &
                                              , __LINE__, __FILE__)
            end if
         end do Generators
         iret = f_selectParentBlock()
         iret = ip - 1
      else
         iret = 0
!!$         if(.not.prealloc) call set_jgen_default_tl(1)
      end if
    end subroutine set_jgen_tl

    subroutine set_iaf_jaf(iret)
      integer, intent(out) :: iret
      integer :: ip
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_readSymmetryGenerator
      integer :: txu,txd,tyu,tyd,tzu,tzd
      if( f_selectBlock( tag_af_generator) == 0) then
         ip = 1
         Generators: do while(.true.)
            if(ip == 1) then
               if( f_selectFirstTableLine() /= 0) then
                  exit Generators
               end if
            else
               if( f_selectNextTableLine() /= 0) then
                  exit Generators
               end if
            end if
            if( f_readSymmetryGenerator(il, rint ,txu,txd,tyu,tyd,tzu,tzd, tag_rotation, tag_tx,tag_ty,tag_tz,printable) == 0) then
                  if(ip > 1) exit
                  iaf = rint
                  jaf(1,1) = txu; jaf(2,1) = txd
                  jaf(1,2) = tyu; jaf(2,2) = tyd
                  jaf(1,3) = tzu; jaf(2,3) = tzd

                  if(ipriinputfile >= 1 .and. printable) then
                     write(nfout,'(" !** iaf  = ",i6,"  jaf  = ",6i5)') iaf, jaf(1,1), jaf(2,1) &
                          & , jaf(1,2),jaf(2,2),jaf(1,3),jaf(2,3)
                  end if
               ip = ip + 1
            end if
         end do Generators
         iret = f_selectParentBlock()
         iret = ip - 1
      end if
    end subroutine set_iaf_jaf

    subroutine set_igen_jgen_default(i)
      integer, intent(in) :: i
      integer :: j
      igen(i) = 1
      do j = 1, 3
         jgen(1,j,i) = 0
         jgen(2,j,i) = 1
      end do
    end subroutine set_igen_jgen_default

    subroutine set_jgen_default_tl(i)
      integer, intent(in) :: i
      integer :: j
!!$      igen_tl(i) = 1
      do j = 1, 3
         jgen_tl(1,j,i) = 0
         jgen_tl(2,j,i) = 1
      end do
    end subroutine set_jgen_default_tl

    subroutine set_lattice_system()
      integer :: iret
      logical :: proper_lattice_system
      proper_lattice_system = .false.
      iret = f_getstringValue( tag_lattice_system, rstr, LOWER)
      if(iret /= 0) iret = f_getstringValue( tag_system, rstr, LOWER)
      if(iret == 0) then
         proper_lattice_system = set_lattice_system0(rstr) ! -->il
         if(.not.proper_lattice_system) then
            if(ipriinputfile>= 0 .and. printable)  &
                 & write(nfout,'(" !* lattice_system defined in the [tspace] block is not approprieate (",a32,")")') rstr
            il = PRIMITIVE_lattice
            if(ipriinputfile>= 0 .and. printable) &
                 & write(nfout,'(" !* lattice_system is set to be PRIMITIVE lattice at <<m_CS_rd_n.set_lattice_system>>")')
         end if
      else
         il = PRIMITIVE_lattice !!! K. Mae 040121
         if(ipriinputfile >= 0 .and. printable) then
            write(nfout,'(" ! lattice_system is not defined in the [tspace] block")')
            write(nfout,'(" ! lattice_system is set to be PRIMITIVE lattice at <<m_CS_rd_n.set_lattice_system>>")')
         end if
      end if
      if(ipriinputfile >= 1 .and. printable) then
         write(nfout,'(" !** il = ",i5," <<m_CS_rd_n.set_lattice_system>>")') il
         if(il == -1) write(nfout,'(" !** il = ",a32)') tag_trigonal
         if(il ==  0) write(nfout,'(" !** il = ",a32)') tag_hexagonal
         if(il ==  1) write(nfout,'(" !** il = ",a32)') tag_primitive
         if(il ==  2) write(nfout,'(" !** il = ",a32)') tag_facecentered
         if(il ==  3) write(nfout,'(" !** il = ",a32)') tag_bodycentered
         if(il ==  4) write(nfout,'(" !** il = ",a32)') tag_basecentered
      end if
    end subroutine set_lattice_system

    logical function set_lattice_system0(rstr)
      character(len=FMAXVALLEN), intent(in) :: rstr
      logical :: tf
      integer :: lenstr
      set_lattice_system0 = .false.
      lenstr = len(trim(rstr))
      if(lenstr >= 1) call strncmp0(tag_rhombohedral, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 1) call strncmp0(tag_trigonal, trim(rstr), tf)
      if(.not.tf .and. lenstr == 2) call strncmp0(trim(rstr), "-1", tf)
      if(tf) then
         il = TRIGONAL_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if

      if(lenstr >= 1) call strncmp0(tag_hexagonal, trim(rstr), tf)
      if(.not.tf .and. lenstr == 1) call strncmp0(trim(rstr), "0", tf)
      if(tf) then
         il = HEXAGONAL_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if

      if(lenstr >= 1) call strncmp0(tag_primitive, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 1) call strncmp0(tag_simple, trim(rstr), tf)
      if(.not.tf .and. lenstr == 1) call strncmp0(trim(rstr),"1",tf)
      if(tf) then
         il = PRIMITIVE_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if

      if(lenstr >= 1) call strncmp0(tag_facecentered, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 1) call strncmp0(tag_fcc, trim(rstr), tf)
      if(.not.tf .and. lenstr == 1) call strncmp0(trim(rstr),"2",tf)
      if(tf) then
         il = FCC_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if

      if(lenstr >= 3) call strncmp0(tag_bodycentered, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 2) call strncmp0(tag_bcc, trim(rstr), tf)
      if(.not.tf .and. lenstr == 1) call strncmp0(trim(rstr),"3",tf)
      if(tf) then
         il = BCC_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if

      if(lenstr >= 3) call strncmp0(tag_bottomcentered, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 2) call strncmp0(tag_basecentered, trim(rstr), tf)
      if(.not.tf .and. lenstr >= 1) call strncmp0(tag_onefacecentered, trim(rstr),tf)
      if(.not.tf .and. lenstr == 1) call strncmp0(trim(rstr),"4",tf)
      if(tf) then
         il = BottomCenteredCubic_lattice
         set_lattice_system0 = .true.
         goto 1001
      end if
1001  continue
    end function set_lattice_system0

!!$    integer function f_readUnitCell(avec,bvec,cvec,a,b,c,alpha,beta,gamma)
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! intelligent reading of unit cell
!!$! This function has been coded by K. Mae(adv), May-Jun. 2003
!!$! Subtle modifications( (8) -> (DP), PI -> PAI ) were done
!!$!        by  T. Yamasaki(FUJITSU LABORATORIES LTD.),  9th Jun. 2003
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      real(DP), dimension(3), intent(out) :: avec
!!$      real(DP), dimension(3), intent(out) :: bvec
!!$      real(DP), dimension(3), intent(out) :: cvec
!!$      real(DP), intent(out) :: a, b, c
!!$      real(DP), intent(out) :: alpha, beta, gamma
!!$      real(DP) :: s, al, be, ga, cos_al, sin_be, cos_be, sin_ga, cos_ga
!!$      integer :: getRealVectorValue
!!$      integer :: getRealValue
!!$      integer :: iret
!!$
!!$      if( getRealVectorValue( trim(TAG_AVECTOR)//char(0), avec ) == 0 ) then
!!$         iret = getRealVectorValue( trim(TAG_BVECTOR)//char(0), bvec )
!!$         iret = getRealVectorValue( trim(TAG_CVECTOR)//char(0), cvec )
!!$         a = sqrt( sum(avec(:)*avec(:)) )
!!$         b = sqrt( sum(bvec(:)*bvec(:)) )
!!$         c = sqrt( sum(cvec(:)*cvec(:)) )
!!$         alpha = acos( dot_product(bvec,cvec)/(b*c) )*180.d0/PAI
!!$         beta = acos( dot_product(cvec,avec)/(c*a) )*180.d0/PAI
!!$         gamma = acos( dot_product(avec,bvec)/(a*b) )*180.d0/PAI
!!$      else if( getRealValue( trim(TAG_A)//char(0), a ) == 0 ) then
!!$         iret = getRealValue( trim(TAG_B)//char(0), b )
!!$         iret = getRealValue( trim(TAG_C)//char(0), c )
!!$         iret = getRealValue( trim(TAG_ALPHA)//char(0), alpha )
!!$         iret = getRealValue( trim(TAG_BETA)//char(0), beta )
!!$         iret = getRealValue( trim(TAG_GAMMA)//char(0), gamma )
!!$         al = alpha*PAI/180.d0
!!$         be = beta*PAI/180.d0
!!$         ga = gamma*PAI/180.d0
!!$         cos_al = cos(al)
!!$         sin_be = sin(be)
!!$         cos_be = cos(be)
!!$         sin_ga = sin(ga)
!!$         cos_ga = cos(ga)
!!$         s = (cos_ga*cos_be-cos_al)/sin_ga/cos_be
!!$         avec(1) = a
!!$         avec(2) = 0.d0
!!$         avec(3) = 0.d0
!!$         bvec(1) = b*cos_ga
!!$         bvec(2) = b*sin_ga
!!$         bvec(3) = 0.d0
!!$
!!$         cvec(1) = c*cos_be
!!$         cvec(2) = (b*c*cos_al - bvec(2)*cvec(1))/bvec(2)
!!$         cvec(3) = sqrt(c*c - cvec(1)**2 - cvec(2)**2)
!!$		cvec(1) = c*cos_be
!!$		cvec(2) = c*abs(s)*cos_be
!!$		cvec(3) = c*sqrt(sin_be*sin_be-s*s*cos_be*cos_be)
!!$      else
!!$         f_readUnitCell = -1
!!$         return
!!$      end if
!!$	
!!$      f_readUnitCell = 0
!!$      return
!!$    end function f_readUnitCell
  end subroutine m_CS_rd_n

  subroutine m_CS_altv_2_rltv()
    call altv_2_rltv(altv,rltv,univol,rvol)
  end subroutine m_CS_altv_2_rltv

  subroutine m_CS_set_nbztyp_spg()
    nbztyp_spg = nbztyp
  end subroutine m_CS_set_nbztyp_spg

  subroutine m_CS_set_nbztyp(type)
    integer, intent(in) :: type
    nbztyp    = type
  end subroutine m_CS_set_nbztyp

  subroutine m_CS_alloc_op_tau(nfout)
    integer, intent(in) :: nfout
    if(printable) then
       write(nfout,'(" -- m_CS_alloc_op_tau --")')
       write(nfout,'(" !! nopr, af = ",2i5)') nopr,af
    end if
    if(nopr == 0) call phase_error_with_msg(nfout, ' nopr is invalid',__LINE__,__FILE__)
    if(allocated(op)) deallocate(op)
    if(allocated(tau)) deallocate(tau)
    allocate(op(3,3,nopr+af))
    allocate(tau(3,nopr+af,CRDTYP))
  end subroutine m_CS_alloc_op_tau

  subroutine m_CS_alloc_op_tau_tl(nfout)
    integer, intent(in) :: nfout
    if(printable) then
       write(nfout,'(" -- m_CS_alloc_op_tau_tl --")')
       write(nfout,'(" !! ngen_tl = ",i5)') ngen_tl
    end if
    if(ngen_tl > 0) then
       if(allocated(op_tl)) deallocate(op_tl)
       allocate(op_tl(3,3,ngen_tl))
       if(allocated(tau_tl)) deallocate(tau_tl)
       allocate(tau_tl(3,ngen_tl,CRDTYP))
    else
       if(.not.allocated(op_tl)) allocate(op_tl(3,3,1))
       if(.not.allocated(tau_tl)) allocate(tau_tl(3,1,CRDTYP))
    end if
  end subroutine m_CS_alloc_op_tau_tl

  subroutine alloc_igen_jgen
    allocate(igen(ngen)); igen = 0
    allocate(jgen(2,3,ngen)); jgen = 0
  end subroutine alloc_igen_jgen

  subroutine alloc_jgen_tl
!!$    allocate(igen_tl(ngen_tl)); igen_tl = 0
    allocate(jgen_tl(2,3,ngen_tl)); jgen_tl = 0
  end subroutine alloc_jgen_tl

  subroutine dealloc_igen_jgen
    deallocate(igen)
    deallocate(jgen)
  end subroutine dealloc_igen_jgen

  subroutine m_CS_wd_CS_data(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    write(nfcntn) inversion_symmetry
    write(nfcntn) nbztyp
    write(nfcntn) altv,rltv, univol, rvol
  end subroutine m_CS_wd_CS_data

  subroutine m_CS_rd_CS_data(nfcntn,nfout)
    integer, intent(in) :: nfcntn, nfout
    if(printable)  write(nfout,*) ' <<< m_CS_rd_CS_data >>>'
!!$    read(nfcntn) input_coordinate_system
!!$    write(nfout,*) ' input_coordinate_system = ', input_coordinate_system
    read(nfcntn) inversion_symmetry
    read(nfcntn) nbztyp
    read(nfcntn) altv,rltv, univol, rvol
  end subroutine m_CS_rd_CS_data

! --> T. Yamasaki, 18th Aug. 2009
  subroutine m_CS_wd_fix_spin_status(nfcntn)

    integer, intent(in) :: nfcntn

! ================================== added by K. Tagami ====================== 11.0
    if ( noncol ) return
! ============================================================================ 11.0

    if(nspin == 2) then
       if(mype == 0) then
          write(nfcntn,'(a15)') tag_fix_spin_status
          write(nfcntn,'(2i5,d20.8)') sw_fix_total_spin_in, sw_fix_total_spin, total_spin
       end if
    end if
  end subroutine m_CS_wd_fix_spin_status

  subroutine m_CS_rd_fix_spin_status(nfcntn,nfout)
    integer, intent(in) :: nfcntn,nfout
    integer :: i
    logical :: EOF_reach, tag_is_found

! ================================== added by K. Tagami ====================== 11.0
    if ( noncol ) return
! ============================================================================ 11.0

    if(nspin == 2) then
       if(mype == 0) then
          call rewind_to_tag0(nfcntn,len(tag_fix_spin_status),tag_fix_spin_status &
               &, EOF_reach, tag_is_found, str, len_str)
          write(nfout,'(" << m_CS_rd_fix_spin_status>>")')
          if(.not.tag_is_found) then
             write(nfout,'(" sw_fix_total_spin_in, sw_fix_total_spin and total_spin  are not read")')
          else
             read(nfcntn,*) i, sw_fix_total_spin, total_spin
             write(nfout,'(" sw_fix_total_spin_in, sw_fix_total_spin, total_spin = ",2i5,d20.8)') &
                  &         i, sw_fix_total_spin, total_spin
          end if
       end if
    end if
  end subroutine m_CS_rd_fix_spin_status
! <--

! --- subroutine m_CS_read_symmtry_op_etc ---
!  This subroutine is revised version of the subroutine <rdprp>
!                 T. Yamasaki(FUJITSU Lab.) 31 May 2003
  subroutine m_CS_rd_symmetry_op_etc(nf)
    integer, intent(in) :: nf
    integer :: ianti,janti(2,3),i,j,k
    character(len=60) :: cname
    logical :: skip_to_tagbegin

    rewind nf

    if(.not.skip_to_tagbegin(nf,tag_tspace)) then
       if(printable) write(6,'(" ! tag of tspace is not found in the inputfile")')
       call phase_error_with_msg(6, ' stopped at <<m_CS_rd_symmetry_op_etc>>',__LINE__,__FILE__)
    end if

    read(nf,*) i ! jpr

    read(nf,'(a60)') cname
    read(nf,*) i, il, ngen, inv  ! i<-idim

    if(ipri_spg >= 1 .and. printable) then
       write(6,820) cname,i,il,ngen,inv
       write(6,*) 'ngen=',ngen     
    end if
  820 format(' == ',a60,' ==' / &
     &  ' dimension=',i2,'      il=',i2,'   ngen=',i2,'   inv=',i2)

    call alloc_igen_jgen()
    do i = 1, ngen
       read (nf,*)   igen(i),((jgen(j,k,i),j=1,2),k=1,3)
    end do

    read(nf,*) imag
    if(imag == 1) then
       if(printable) write(6,*) ' antiferromagnetic calculation'
       read (nf,*)   ianti,((janti(j,k),j=1,2),k=1,3)
    endif

!!$    read(nf,*) a,b,c,ca,cb,cc

    if(ipri_spg >= 1 .and. printable) then
       write(6,*) ' << m_CS_rd_symmetry_op_etc >>'
       write(6,860) a,b,c,ca,cb,cc
860    format(' a, b, c =',3f12.6/ &
       &       'ca,cb,cc =',3f12.6)
    end if
  end subroutine m_CS_rd_symmetry_op_etc

  subroutine m_CS_set_abccacbcc()
    a=dsqrt(altv(1,1)**2+altv(2,1)**2+altv(3,1)**2)
    b=dsqrt(altv(1,2)**2+altv(2,2)**2+altv(3,2)**2)
    c=dsqrt(altv(1,3)**2+altv(2,3)**2+altv(3,3)**2)
    ca=(altv(1,2)*altv(1,3)+altv(2,2)*altv(2,3)+altv(3,2)*altv(3,3)) /(b*c)
    cb=(altv(1,3)*altv(1,1)+altv(2,3)*altv(2,1)+altv(3,3)*altv(3,1)) /(c*a)
    cc=(altv(1,1)*altv(1,2)+altv(2,1)*altv(2,2)+altv(3,1)*altv(3,2)) /(a*b)
    if(printable) then 
       write(6,*) ' << m_CS_set_abccacbcc >>  '
       write(6,860) a,b,c,ca,cb,cc
    end if
860 format(' a, b, c =',3f12.6/ 'ca,cb,cc =',3f12.6)
  end subroutine m_CS_set_abccacbcc

  ! --- subroutine m_CS_set_default_symm_op
  !  This subroutine is revised from <setspg_default> by T. Yamasaki(FUJITSU Lab.)
  !  The former part of the <setspg_default> is taken.
  !                                                            31 May 2003
  subroutine m_CS_set_default_symm_op(nfout)
    integer, intent(in) :: nfout
    integer ::   ip

    if(ipri_spg >= 1 .and. printable)  write(nfout,'(" !** nbztyp_spg = ",i6)') nbztyp_spg

!!$    if(ipri_spg >= 1 .and. printable) then
!!$       write(nfout,*) ' '
!!$       write(nfout,860) a,b,c,ca,cb,cc
!!$860    format(' !** a, b, c =',3f12.6/ ' !** ca,cb,cc =',3f12.6)
!!$    end if
    if(nbztyp_spg == SIMPLE_CUBIC) then
       il=1
       ngen=3
       inv=1
    else if(nbztyp_spg == BCC) then
       il=3
       ngen=3
       inv=1
    else if(nbztyp_spg == FCC) then
       il=2
       ngen=3
       inv=1
    else if(nbztyp_spg == DIAMOND) then
       il=2
       ngen=3
       inv=0
    else if(nbztyp_spg == HEXAGONAL) then
       il=0
       ngen=1       !!! K.Mae 040315
       !!! ngen=2   !!! K.Mae 040315
       inv=0
    else if(nbztyp_spg == HCP) then
       il=0
       ngen=3
       inv=1
    else if(nbztyp_spg == TRIGONAL) then
       il=-1
       ngen=1
       inv=0
    else
       il=1
       ngen=1
       inv=0
       if(printable) write(nfout,'(" !** nbztyp_spg = ",i6," <m_CS_set_default_symm_op>")') nbztyp_spg
!!$       stop ' nbztyp_spg (m_CS_set_default_symm_op) is illegal'
    end if
    call alloc_igen_jgen()

    if(nbztyp_spg == SIMPLE_CUBIC) then     ! Oh
       igen(1)=5
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
       igen(2)=19
       jgen(1,1,2)=0;         jgen(2,1,2)=1
       jgen(1,2,2)=0;         jgen(2,2,2)=1
       jgen(1,3,2)=0;         jgen(2,3,2)=1
       igen(3)=25
       jgen(1,1,3)=0;         jgen(2,1,3)=1
       jgen(1,2,3)=0;         jgen(2,2,3)=1
       jgen(1,3,3)=0;         jgen(2,3,3)=1
    else if(nbztyp_spg == BCC) then
       igen(1)=5
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
       igen(2)=19
       jgen(1,1,2)=0;         jgen(2,1,2)=1
       jgen(1,2,2)=0;         jgen(2,2,2)=1
       jgen(1,3,2)=0;         jgen(2,3,2)=1
       igen(3)=25
       jgen(1,1,3)=0;         jgen(2,1,3)=1
       jgen(1,2,3)=0;         jgen(2,2,3)=1
       jgen(1,3,3)=0;         jgen(2,3,3)=1
    else if(nbztyp_spg == FCC) then
       igen(1)=5
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
       igen(2)=19
       jgen(1,1,2)=0;         jgen(2,1,2)=1
       jgen(1,2,2)=0;         jgen(2,2,2)=1
       jgen(1,3,2)=0;         jgen(2,3,2)=1
       igen(3)=25
       jgen(1,1,3)=0;         jgen(2,1,3)=1
       jgen(1,2,3)=0;         jgen(2,2,3)=1
       jgen(1,3,3)=0;         jgen(2,3,3)=1
    else if(nbztyp_spg == DIAMOND) then
       igen(1)=5
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
!ccccccc 2nd choice ccccccccccccccccccccc
       igen(2)=19
       jgen(1,1,2)=1;         jgen(2,1,2)=4
       jgen(1,2,2)=1;         jgen(2,2,2)=2
       jgen(1,3,2)=3;         jgen(2,3,2)=4
       igen(3)=25
       jgen(1,1,3)=0;         jgen(2,1,3)=1
       jgen(1,2,3)=0;         jgen(2,2,3)=1
       jgen(1,3,3)=0;         jgen(2,3,3)=1
!ccccccc 2nd choice ccccccccccccccccccccc
    else if(nbztyp_spg == HEXAGONAL) then
      !!! K.Mae 040315, only C3
      igen(1)=3
      jgen(1,1,1)=0;         jgen(2,1,1)=1
      jgen(1,2,1)=0;         jgen(2,2,1)=1
      jgen(1,3,1)=0;         jgen(2,3,1)=1
      !!! commented out by K.Mae 040315, These generators only for SiO2
      ! igen(1)=3
      ! jgen(1,1,1)=0;         jgen(2,1,1)=1
      ! jgen(1,2,1)=0;         jgen(2,2,1)=1
      ! jgen(1,3,1)=2;         jgen(2,3,1)=3
      ! igen(2)=10
      ! jgen(1,1,2)=0;         jgen(2,1,2)=1
      ! jgen(1,2,2)=0;         jgen(2,2,2)=1
      ! jgen(1,3,2)=0;         jgen(2,3,2)=1
      !!! end commented out by K.Mae 040315
    else if(nbztyp_spg == HCP) then
! [nbztyp_spg == HCP] is added by T. Yamasaki(FUJITSU LABORATORIES LTD.) 11th Jun. 2003
       igen(1)=2
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=1;         jgen(2,3,1)=2
       igen(2)=7
       jgen(1,1,2)=0;         jgen(2,1,2)=1
       jgen(1,2,2)=0;         jgen(2,2,2)=1
       jgen(1,3,2)=1;         jgen(2,3,2)=2
       igen(3)=13
       jgen(1,1,3)=0;         jgen(2,1,3)=1
       jgen(1,2,3)=0;         jgen(2,2,3)=1
       jgen(1,3,3)=0;         jgen(2,3,3)=1
    else if(nbztyp_spg == TRIGONAL) then
       igen(1)=3
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
    else
       igen(1)=1
       jgen(1,1,1)=0;         jgen(2,1,1)=1
       jgen(1,2,1)=0;         jgen(2,2,1)=1
       jgen(1,3,1)=0;         jgen(2,3,1)=1
    end if
    if(ipri_spg >= 1 .and. printable) then
       write(nfout,'(" !** il = ",i6, "<<m_CS_set_default_symm_op>>")') il
       do ip = 1, ngen
          write(nfout,'(" !** igen = ",i6,"  jgen = ",6i5)') igen(ip), jgen(1,1,ip), jgen(2,1,ip) &
               & , jgen(1,2,ip),jgen(2,2,ip),jgen(1,3,ip),jgen(2,3,ip)
       end do
    end if
  end subroutine m_CS_set_default_symm_op

  ! ---------- 
  subroutine m_CS_wd_op_and_tau(nfout)
    integer, intent(in) :: nfout
    integer :: no ,i, j, l, k
    real(kind=DP),dimension(3,3) :: op_br
    write(nfout,*) ' --- Symmetry Operations (CARTS, PUCV) --- (in m_CS_wd_op_and_tau)'
    write(nfout,'(" !! nopr+af = ",i6)') nopr+af

    do no = 1, nopr + af
       do i = 1, 3
          do j = 1, 3
             op_br(i,j) = 0.d0
             do l = 1,3
                do k = 1,3
                   op_br(i,j) = op_br(i,j) &
                        &           + altv(l,i)*op(l,k,no)*rltv(k,j)
                enddo
             enddo
             op_br(i,j) = op_br(i,j)/PAI2
          enddo
       enddo
       if( no<= nopr ) then
          if(il<=0 .or. pg_symbol_system == 'hexagonal') then
             write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(ig01(no))
          else
             write(nfout,*) ' #symmetry op. = ', no, oh_symbol(ig01(no))
          end if
       else
          if(il<=0 .or. pg_symbol_system == 'hexagonal') then
             write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(iaf)
          else
             write(nfout,*) ' #symmetry op. = ', no, oh_symbol(iaf)
          end if
       end if
         
       do i = 1, 3
          write(nfout,'(4f8.4,3x,4f8.4)') &
               &           (op   (i,j,no),j=1,3),tau(i,no,CARTS) &
               &          ,(op_br(i,j),j=1,3),tau(i,no,BUCS)
       enddo
    end do
    if(ipri >= 1 .and. symmetry_not_printed) symmetry_not_printed=.false.
  end subroutine m_CS_wd_op_and_tau

  subroutine m_CS_gnrt_symmetry_operations(paramset,nfout)
    logical, intent(in) :: paramset
    integer, intent(in) :: nfout
    integer ::             i, j, iopr, ip
    integer ::             af_t
    integer ::             id_sname = -1
    real(kind=DP),dimension(3,3) :: op_br
#ifdef __TIMER_SUB__
  call timer_sta(1227)
#endif

    call tstatc0_begin('m_CS_gnrt_symmetry_operations ',id_sname)
! ----- iopr, op, tau

! ---> T. Yamasaki(FUJITSU Laboratoreis) 31 May 2003
!!$    if(paramset) then
!!$       if(nbztyp >= 2 .and. &
!!$            & (nbztyp_spg == GENERAL .or. nbztyp_spg == GENERAL_LARGER)) then
!!$          call m_CS_rd_symmetry_op_etc(nfspg)
!!$       else
!!$          call m_CS_set_default_symm_op()  ! <-- nbztyp_spg
!!$       end if
!!$    end if
! <----

    if(printable .and. paramset) then
       write(nfout,'(" ngen = ", i6)') ngen
       do ip = 1, ngen
          write(nfout,'(" igen = ",i6," jgen = ",6i5)') igen(ip) &
               & , jgen(1,1,ip), jgen(2,1,ip), jgen(1,2,ip) &
               & , jgen(2,2,ip), jgen(1,3,ip), jgen(2,3,ip)
       end do
       if(af > 0) then
          write(nfout,'(" iaf  = ",i6," jaf  = ",6i5)') iaf &
                & , jaf(1,1), jaf(2,1), jaf(1,2) &
                & , jaf(2,2), jaf(1,3), jaf(2,3)
       end if
    end if
    call iopr_op_tau(af_t)  ! depends on 'paramset'

    if(ipri >= 3 .and. paramset) write(nfout,'(" !! af_t = ",i10)') af_t
    if(af_t > 1) call phase_error_with_msg(nfout, ' ! illegal af_t value',__LINE__,__FILE__)
    call m_CtrlP_set_af(af_t)
    if(ipri >= 2 .and. paramset) write(nfout,'(" !! af_t (after m_CtrlP_set_af) = " ,i6)') af_t

!!$    if(ipri >= 0 .and. printable) then
    if((.not.paramset .and. ipri >= 3 ) .or. (paramset .and. ipri >= 0 .and. printable)) then
       write(nfout,*) ' << rltv, altv >>'
       write(nfout,'(3f8.4,11x,3f8.4)') &
            & ((rltv(i,j),j=1,3),(altv(i,j),j=1,3),i=1,3)
    endif
    if(paramset) then
       nopr = iopr
       goto 1001
    end if

    call tau_in_PUCV()

!!$    if(ipri>=1) call wd_op_and_tau()

1001 continue

!! --------------- supercell (redundant space-group) -------------
!  transfered by T. Yamasaki from a revised subroutine of fd_symmetry by Usami-san
!                                                          March 2014

!!$    call gnrt_operators_for_supsercell()

!! ---------------------------------------------------------------    

    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
  call timer_end(1227)
#endif
  contains
    subroutine wd_op_and_tau()
      integer :: no ,i, j, l, k
      write(nfout,*) ' --- Symmetry Operations (CARTS, PUCV) --- (in generator)'
      write(nfout,'(" !! nopr+af = ",i6)') nopr+af

      do no = 1, nopr + af
         do i = 1, 3
            do j = 1, 3
               op_br(i,j) = 0.d0
               do l = 1,3
                  do k = 1,3
                     op_br(i,j) = op_br(i,j) &
                          &           + altv(l,i)*op(l,k,no)*rltv(k,j)
                  enddo
               enddo
               op_br(i,j) = op_br(i,j)/PAI2
            enddo
         enddo
         if( no<= nopr ) then
            if(il<=0 .or. pg_symbol_system == 'hexagonal') then
               write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(ig01(no))
            else
               write(nfout,*) ' #symmetry op. = ', no, oh_symbol(ig01(no))
            end if
         else
            if(il<=0 .or. pg_symbol_system == 'hexagonal') then
               write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(iaf)
            else
               write(nfout,*) ' #symmetry op. = ', no, oh_symbol(iaf)
            end if
         end if
         
         do i = 1, 3
            write(nfout,'(4f8.4,3x,4f8.4)') &
                 &           (op   (i,j,no),j=1,3),tau(i,no,CARTS) &
                 &          ,(op_br(i,j),j=1,3),tau(i,no,BUCS)
         enddo
      end do
    end subroutine wd_op_and_tau

    subroutine iopr_op_tau(af_t)  ! --> iopr,af,op,tau 
      integer, intent(out) :: af_t
      integer ::              ipri_spg_t

      af_t = 0

      if(paramset) then
! ====================================== modified by K. Tagami ============ 12.0A
!         call gnrt_op_paramset(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc &
!              &        ,nfout,iopr,af_t,nbztyp_spg,ig01) !-(b_Crystal_Structure)

         call gnrt_op_paramset( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
              &                 a, b, c, ca, cb, cc, &
              &                 nfout, iopr, af_t, nbztyp_spg, ig01, &
              &                 use_altv_rltv, altv, rltv,&
              &                 gen_name_in_carts )
                                                      !-(b_Crystal_Structure)
! ========================================================================== 12.0A
      else
         ipri_spg_t = ipri_spg
         if(nbztyp < 2) then
            iopr = 1
            op = 0.d0
            tau = 0.d0
            op(1,1,iopr) = 1.d0
            op(2,2,iopr) = 1.d0
            op(3,3,iopr) = 1.d0
         end if
!-----  nbztype =100 & 101   - by Tsuyoshi Miyazaki '94.8.16
!      read operation matrix and k-point from the files
!        1.opgr.data  : operation matrix
!        2.matrix.BP  : transformation matrix from
!                       Bravais lattice to Primitive lattice
!        3.kpoint.dat : irreducible k-point
!      these files can be made by using Dr Hamada-san's code
!     (these files must be located at the working directory)

! ====================================== modified by K. Tagami ============ 12.0A
!         call gnrt_op_n(il,ngen,inv,igen,jgen,imag,iaf,jaf,a,b,c,ca,cb,cc &
!              &        ,nopr+af,altv,nfout,iopr,op,tau, af_t, nbztyp_spg,ipri_spg_t,ig01) !-(b_Crystal_Structure)
!
         call gnrt_op_n( il, ngen, inv, igen, jgen, imag, iaf, jaf, &
              &          a, b, c, ca, cb, cc, &
              &          nopr+af, altv, nfout, iopr, op, tau, af_t, &
              &          nbztyp_spg, ipri_spg_t, ig01, &
              &          use_altv_rltv, rltv, &
              &          gen_name_in_carts ) 
                                          !-(b_Crystal_Structure)
! ============================================================================= 12.0A
      endif

    end subroutine iopr_op_tau

    subroutine tau_in_PUCV
      integer i, no
      do i = 1, 3
         do no = 1, iopr + af
            tau(i,no,BUCS) = (rltv(1,i)*tau(1,no,CARTS) &
                 &     +      rltv(2,i)*tau(2,no,CARTS) &
                 &     +      rltv(3,i)*tau(3,no,CARTS))/PAI2
         enddo
      enddo
    end subroutine tau_in_PUCV
  end subroutine m_CS_gnrt_symmetry_operations

  ! ---------- 
  subroutine m_CS_gnrt_symm_operators_tl(paramset,nfout)
    logical, intent(in) :: paramset
    integer, intent(in) :: nfout
    integer ::             i, j, ip
    integer ::             af_t
    integer ::             id_sname = -1

    if(sw_supercell_symmetry == OFF) return

    call tstatc0_begin('m_CS_gnrt_symm_operators_tl ',id_sname)
! -----  op, tau

!!$    if(printable .and. paramset) then
       write(nfout,'(" ngen_tl = ", i6)') ngen_tl
       do ip = 1, ngen_tl
          write(nfout,'(" jgen_tl = ",6i5)') jgen_tl(1:2,1,ip), jgen_tl(1:2,2,ip), jgen_tl(1:2,3,ip)
       end do
!!$    end if
    if(paramset) goto 1001

    call op_tau_tl()  ! depends on 'paramset'
    call tau_in_CARTS()

    if(ipri>=1) call wd_op_and_tau()

1001 continue

    call tstatc0_end(id_sname)
  contains
    subroutine wd_op_and_tau()
      integer :: no ,i, j, l, k
      real(kind=DP),dimension(3,3) :: op_br

      write(nfout,*) ' --- SC Symmetry Operations (CARTS, PUCV) --- (in m_CS_gnrt_symm_operators_tl)'
      write(nfout,'(" !! ngen_tl = ",i6)') ngen_tl

      do no = 1, ngen_tl
         do i = 1, 3
            do j = 1, 3
               op_br(i,j) = 0.d0
               do l = 1,3
                  do k = 1,3
                     op_br(i,j) = op_br(i,j) &
                          &           + altv(l,i)*op_tl(l,k,no)*rltv(k,j)
                  enddo
               enddo
               op_br(i,j) = op_br(i,j)/PAI2
            enddo
         enddo
         if(il<=0 .or. pg_symbol_system == 'hexagonal') then
            write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(ig01(1))
         else
            write(nfout,*) ' #symmetry op. = ', no, oh_symbol(ig01(1))
         end if
         do i = 1, 3
            write(nfout,'(4f8.4,3x,4f8.4)') &
                 &           (op_tl(i,j,no),j=1,3),tau_tl(i,no,CARTS) &
                 &          ,(op_br(i,j),j=1,3),tau_tl(i,no,BUCS)
         enddo
      end do
    end subroutine wd_op_and_tau

    subroutine op_tau_tl()  ! --> op,tau 
      integer ::              i

      if(ngen_tl < 1) return

      op_tl = 0.d0
      do i = 1, ngen_tl
         op_tl(1,1,i) = 1.d0
         op_tl(2,2,i) = 1.d0
         op_tl(3,3,i) = 1.d0
         tau_tl(1,i,BUCS) = dble(jgen_tl(1,1,i))/dble(jgen_tl(2,1,i))
         tau_tl(2,i,BUCS) = dble(jgen_tl(1,2,i))/dble(jgen_tl(2,2,i))
         tau_tl(3,i,BUCS) = dble(jgen_tl(1,3,i))/dble(jgen_tl(2,3,i))
      end do
    end subroutine op_tau_tl

    subroutine tau_in_CARTS
      integer i, no
      if(ngen_tl < 1) return
      do i = 1, 3
         do no = 1, ngen_tl
            tau_tl(i,no,CARTS) = (altv(i,1)*tau_tl(1,no,BUCS) &
                 &     +          altv(i,2)*tau_tl(2,no,BUCS) &
                 &     +          altv(i,3)*tau_tl(3,no,BUCS))
         enddo
      enddo
    end subroutine tau_in_CARTS
  end subroutine m_CS_gnrt_symm_operators_tl

! ================================== modified by K. Tagami ========== 11.0Ex
  subroutine m_CS_op_in_PUCV(nfout,op_br, n, print_flag, comment )
! =================================================================== 11.0Ex
    integer, intent(in) :: nfout, n
    real(kind=DP), intent(out), dimension(3,3,n) :: op_br

! ================================== added by K. Tagami ========== 11.0Ex
    logical, optional, intent(in) :: print_flag
    character(len=*), optional, intent(in) :: comment
!
    logical :: print_go_flag
! =================================================================== 11.0Ex
    integer :: i,j,l,k,no
    op_br = 0.d0
! ================================== added by K. Tagami ========== 11.0Ex
    print_go_flag = .true.
    if ( present(print_flag) ) then
       print_go_flag = print_flag
    endif
! =================================================================== 11.0Ex

    if ( iprisym >=2 ) symmetry_not_printed = .true.
!!    if ( noncol ) symmetry_not_printed = .true.

! ================================== modified by K. Tagami ========== 11.0Ex
    if(ipri >= 1 .and. print_go_flag .and. symmetry_not_printed) then
       if ( present(comment) ) then
          write(nfout,*) comment
       endif
       write(nfout,*) ' --- Symmetry Operations (CARTS, PUCV) --- (m_CS_op_in_PUCV)'
       write(nfout,'(" !! nopr+af, n = ",2i6)') nopr+af,n
    end if
! ==================================================================== 11.0Ex

    do no = 1, nopr + af
       do i = 1, 3
          do j = 1, 3
             op_br(i,j,no) = 0.d0
             do l = 1,3
                do k = 1,3
                   op_br(i,j,no) = op_br(i,j,no) &
                        &           + altv(l,i)*op(l,k,no)*rltv(k,j)
                enddo
             enddo
             op_br(i,j,no) = op_br(i,j,no)/PAI2
          enddo
       enddo

! ================================== modified by K. Tagami ========== 11.0Ex
!       if(ipri >= 1 .and. printable) then
       if(ipri >= 1 .and. printable .and. print_go_flag .and. symmetry_not_printed) then
! ==================================================================== 11.0Ex
! ================================= Added by K. Tagami =============
          if( no<= nopr ) then
! ==================================================================
             if(il<=0 .or. pg_symbol_system == 'hexagonal') then
                write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(ig01(no))
             else
                write(nfout,*) ' #symmetry op. = ', no, oh_symbol(ig01(no))
             end if
          else
             if(il<=0 .or. pg_symbol_system == 'hexagonal') then
                write(nfout,*) ' #symmetry op. = ', no, d6h_symbol(iaf)
             else
                write(nfout,*) ' #symmetry op. = ', no, oh_symbol(iaf)
             end if
          end if
          do i = 1, 3
             write(nfout,'(4f8.4,3x,4f8.4)') &
                  &           (op   (i,j,no),j=1,3),tau(i,no,CARTS) &
                  &          ,(op_br(i,j,no),j=1,3),tau(i,no,BUCS)
          enddo
       endif
    enddo
    if(ipri >= 1.and.printable.and.print_go_flag.and.symmetry_not_printed) symmetry_not_printed = .false.
  end subroutine m_CS_op_in_PUCV


! ================================== modified by K. Tagami ========== 11.0Ex
  subroutine m_CS_op_in_PUCD(nfout,op_pr, n, print_flag, comment )
! =================================================================== 11.0Ex

    integer, intent(in) :: nfout, n
    real(kind=DP), intent(out), dimension(3,3,n) :: op_pr
      
! ================================== added by K. Tagami ========== 11.0Ex
    logical, optional, intent(in) :: print_flag
    character(len=*), optional, intent(in) :: comment
!
    logical :: print_go_flag
! =================================================================== 11.0Ex

    integer :: i,j,l,k,no

    op_pr = 0.d0

! ================================== added by K. Tagami ========== 11.0Ex
    print_go_flag = .true.
    if ( present(print_flag) ) then
       print_go_flag = print_flag
    endif
! =================================================================== 11.0Ex

! ================================== modified by K. Tagami ========== 11.0Ex
    if(ipri >= 1 .and. print_go_flag .and. symmetry_not_printed) then
       if ( present(comment) ) then
           write(nfout,*) comment
       endif
       write(nfout,*) ' --- Symmetry Operations (CARTS, PUCD) ---'
       write(nfout,'(" !! nopr+af, n = ",2i6)') nopr+af,n
    end if
! ==================================================================== 11.0Ex

    do no = 1, nopr + af
       do i = 1, 3
          do j = 1, 3
             op_pr(i,j,no) = 0.d0
             do l = 1,3
                do k = 1,3
                   op_pr(i,j,no) = op_pr(i,j,no) &
                        &           + rltv(l,i)*op(l,k,no)*altv(k,j)
                enddo
             enddo
             op_pr(i,j,no) = op_pr(i,j,no)/PAI2
          enddo
       enddo

! ================================== modified by K. Tagami ========== 11.0Ex
!       if(ipri >= 1) then
       if(ipri >= 1 .and. print_go_flag .and. symmetry_not_printed) then
! ==================================================================== 11.0Ex
          write(nfout,*) ' #symmetry op. = ', no
          do i = 1, 3
             write(nfout,'(4f8.4,3x,4f8.4)') &
                  &           (op   (i,j,no),j=1,3),tau(i,no,CARTS) &
                  &          ,(op_pr(i,j,no),j=1,3),tau(i,no,BUCS)
          enddo
       endif
    enddo
    if(ipri >= 1.and.print_go_flag.and.symmetry_not_printed) symmetry_not_printed=.false.
  end subroutine m_CS_op_in_PUCD

  subroutine m_CS_gnrt_tmatrices(il)
! This subroutine was coded by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003

    integer,       intent(in)    :: il
    real(kind=DP)                :: unimat(3,3) = 0.d0
    real(kind=DP), parameter     :: one=1.d0, half=0.5d0, &
         &                          one3rd = 0.33333333333333333d0, &
         &                          two3rd = 0.66666666666666667d0, &
         &                          epsilon = 1.d-10
    integer                      :: ido, jdo

    b2pmat(1,1) = one ; b2pmat(2,2) = one ; b2pmat(3,3) = one
    p2bmat(1,1) = one ; p2bmat(2,2) = one ; p2bmat(3,3) = one

    if (il == -1) then                  ! Trigonal

       b2pmat(1,1) =  two3rd ; b2pmat(1,2) =  one3rd ; b2pmat(1,3) =  one3rd
       b2pmat(2,1) = -one3rd ; b2pmat(2,2) =  one3rd ; b2pmat(2,3) =  one3rd
       b2pmat(3,1) = -one3rd ; b2pmat(3,2) = -two3rd ; b2pmat(3,3) =  one3rd

       p2bmat(1,1) =  one ; p2bmat(1,2) = -one ; p2bmat(1,3) = 0.d0
       p2bmat(2,1) = 0.d0 ; p2bmat(2,2) =  one ; p2bmat(2,3) = -one
       p2bmat(3,1) =  one ; p2bmat(3,2) =  one ; p2bmat(3,3) =  one

    else if (il == 0) then              ! Primitive (Hexagonal)
    else if (il == 1) then              ! Primitive
    else if (il == 2) then              ! Face-centered

       b2pmat(1,1) =  0.d0  ;   b2pmat(1,2) =  half  ;   b2pmat(1,3) =  half
       b2pmat(2,1) =  half  ;   b2pmat(2,2) =  0.d0  ;   b2pmat(2,3) =  half
       b2pmat(3,1) =  half  ;   b2pmat(3,2) =  half  ;   b2pmat(3,3) =  0.d0

       p2bmat(1,1) = -one ;  p2bmat(1,2) =  one ;  p2bmat(1,3) =  one
       p2bmat(2,1) =  one ;  p2bmat(2,2) = -one ;  p2bmat(2,3) =  one
       p2bmat(3,1) =  one ;  p2bmat(3,2) =  one ;  p2bmat(3,3) = -one

    else if (il == 3) then              ! Body-centered

       b2pmat(1,1) = -half  ;   b2pmat(1,2) =  half  ;   b2pmat(1,3) =  half
       b2pmat(2,1) =  half  ;   b2pmat(2,2) = -half  ;   b2pmat(2,3) =  half
       b2pmat(3,1) =  half  ;   b2pmat(3,2) =  half  ;   b2pmat(3,3) = -half

       p2bmat(1,1) =  0.d0 ; p2bmat(1,2) =  one  ; p2bmat(1,3) =  one
       p2bmat(2,1) =  one  ; p2bmat(2,2) =  0.d0 ; p2bmat(2,3) =  one
       p2bmat(3,1) =  one  ; p2bmat(3,2) =  one  ; p2bmat(3,3) =  0.d0

    else if (il == 4) then              ! Base-centered

       b2pmat(1,1) =  half  ;   b2pmat(1,2) = -half  ;   b2pmat(1,3) =  0.d0
       b2pmat(2,1) =  half  ;   b2pmat(2,2) =  half  ;   b2pmat(2,3) =  0.d0
       b2pmat(3,1) =  0.d0  ;   b2pmat(3,2) =  0.d0  ;   b2pmat(3,3) =  one

       p2bmat(1,1) =  one  ; p2bmat(1,2) =  one  ; p2bmat(1,3) =  0.d0
       p2bmat(2,1) = -one  ; p2bmat(2,2) =  one  ; p2bmat(2,3) =  0.d0
       p2bmat(3,1) =  0.d0 ; p2bmat(3,2) =  0.d0 ; p2bmat(3,3) =  one

    else
       if(printable) write(6,*) 'The value of il is invalid. il =', il
       !stop
       call phase_error_with_msg(6, 'The value of il is invalid. ',__LINE__,__FILE__)
    end if

    unimat = matmul(b2pmat,p2bmat)
    do ido = 1, 3
       do jdo = 1, 3
          if (ido == jdo) then
             if ((unimat(ido,jdo)-one) > epsilon) call phase_error_with_msg(6, 'not unit matrix' &
                                                  , __LINE__, __FILE__)
          else
             if (unimat(ido,jdo) > epsilon) call phase_error_with_msg(6, 'not unit matrix' &
                                                  , __LINE__, __FILE__)
          end if
       end do
    end do

  end subroutine m_CS_gnrt_tmatrices

  subroutine m_CS_phonon_symmetry(on_or_off)
    use m_Files, only : nfout
    integer, intent(in) :: on_or_off
    integer :: i

    if(on_or_off == ON) then
       nopr = phonon_nopr
       if(printable) then
          do i=1,nopr
             write(nfout,'(i2)') i
             write(nfout,'(3(1x,f10.5),2x,f10.5)') op(1,1:3,i),tau(1,i,CARTS)
             write(nfout,'(3(1x,f10.5),2x,f10.5)') op(2,1:3,i),tau(2,i,CARTS)
             write(nfout,'(3(1x,f10.5),2x,f10.5)') op(3,1:3,i),tau(3,i,CARTS)
          end do
       end if
    else
       phonon_nopr = nopr
       nopr = 1
       ngen = 1 
       igen(ngen)=1
       jgen(1,1:3,ngen) = 0
       jgen(2,1:3,ngen) = 1
    end if

  end subroutine m_CS_phonon_symmetry

  subroutine set_strain_tensor()
    integer :: iret
    integer :: f_getRealValue
    real(kind=DP) :: val

    if(f_getRealValue( "e11", val, '') == 0) strain(1,1) = val
    if(f_getRealValue( "e12", val, '') == 0) strain(1,2) = val
    if(f_getRealValue( "e13", val, '') == 0) strain(1,3) = val
    if(f_getRealValue( "e21", val, '') == 0) strain(2,1) = val
    if(f_getRealValue( "e22", val, '') == 0) strain(2,2) = val
    if(f_getRealValue( "e23", val, '') == 0) strain(2,3) = val
    if(f_getRealValue( "e31", val, '') == 0) strain(3,1) = val
    if(f_getRealValue( "e32", val, '') == 0) strain(3,2) = val
    if(f_getRealValue( "e33", val, '') == 0) strain(3,3) = val

  end subroutine set_strain_tensor

  subroutine set_strained_cell(avec,bvec,cvec)
    real(kind=DP), intent(inout), dimension(3) :: avec,bvec,cvec

    integer :: i,j
    real(kind=DP), dimension(3) :: avec0,bvec0,cvec0

    avec0 = avec; bvec0 = bvec; cvec0 = cvec
    do i=1,3
       do j=1,3
          avec(i) = avec(i) + strain(i,j)*avec0(j)
          bvec(i) = bvec(i) + strain(i,j)*bvec0(j)
          cvec(i) = cvec(i) + strain(i,j)*cvec0(j)
       end do
    end do

  end subroutine set_strained_cell

  subroutine m_CS_supercell(nfout)
    integer, intent(in) :: nfout
    integer :: i,j,k,n,l
    real(kind=DP) :: bltv(3,3)
    
    if(sw_supercell == ON) then
       if(supercell_unit_type == BRAVAIS) then
          do j=1,3
             do i=1,3
                bltv(i,j) = sum(altv(i,1:3)*p2bmat(j,1:3))
             end do
          end do
       else ! PRIMITIVE
          bltv = altv
       end if
       altv_super(1:3,1) = n1_sc*bltv(1:3,1)
       altv_super(1:3,2) = n2_sc*bltv(1:3,2)
       altv_super(1:3,3) = n3_sc*bltv(1:3,3)
       call altv_2_rltv(altv_super,rltv_super,univol_super,rvol_super)  ! in b_CS
       if(printable) then
          write(nfout,'("rltv_super",28x,"altv_super")') 
          do i=1,3
             write(nfout,'(3(1x,f10.5),2x,3(1x,f10.5))') rltv_super(i,1:3), altv_super(i,1:3) 
          end do
       end if

       altv_prim = altv
       rltv_prim = rltv
       univol_prim = univol
       rvol_prim = rvol

       altv = altv_super
       rltv = rltv_super
       univol = univol_super
       rvol = rvol_super

       if(il<=1.or.supercell_unit_type==PRIMITIVE) then
          ! trigonal, hexagonal, primitive
          nlpnt = n1_sc*n2_sc*n3_sc
       else if(il==2) then ! facecentered
          nlpnt = n1_sc*n2_sc*n3_sc*4
       else if(il==3.or.il==4) then ! bodycentered, basecentered
          nlpnt = n1_sc*n2_sc*n3_sc*2
       end if
       allocate(lpnt(nlpnt,3))
       n = 0
       do k=0,n3_sc-1
          do j=0,n2_sc-1
             do i=0,n1_sc-1
                if(il<=1.or.supercell_unit_type==PRIMITIVE) then
                ! trigonal, hexagonal, primitive
                   n = n + 1
                   lpnt(n,1:3) = bltv(1:3,1)*i+bltv(1:3,2)*j+bltv(1:3,3)*k
                else if(il==2) then ! facecentered
                   n = n + 4
                   lpnt(n-3,1:3) = bltv(1:3,1)*i+bltv(1:3,2)*j+bltv(1:3,3)*k
                   lpnt(n-2,1:3) = lpnt(n-3,1:3)+0.5d0*(bltv(1:3,2)+bltv(1:3,3))
                   lpnt(n-1,1:3) = lpnt(n-3,1:3)+0.5d0*(bltv(1:3,1)+bltv(1:3,3))
                   lpnt(n  ,1:3) = lpnt(n-3,1:3)+0.5d0*(bltv(1:3,1)+bltv(1:3,2))
                else if(il==3) then ! bodycentered
                   n = n + 2
                   lpnt(n-1,1:3) = bltv(1:3,1)*i+bltv(1:3,2)*j+bltv(1:3,3)*k
                   lpnt(n  ,1:3) = lpnt(n-1,1:3)+0.5d0*(bltv(1:3,1)+bltv(1:3,2)+bltv(1:3,3))
                else if(il==4) then ! basecentered
                   n = n + 2
                   lpnt(n-1,1:3) = bltv(1:3,1)*i+bltv(1:3,2)*j+bltv(1:3,3)*k
                   lpnt(n  ,1:3) = lpnt(n-1,1:3)+0.5d0*(bltv(1:3,1)+bltv(1:3,2))
                end if
             end do
          end do
       end do

    end if
  end subroutine m_CS_supercell

  subroutine m_CS_supercell_on(nfout,uc_type,n1,n2,n3)
    integer, intent(in) :: nfout,uc_type,n1,n2,n3
    sw_supercell = ON
    supercell_unit_type = uc_type
    n1_sc = n1
    n2_sc = n2
    n3_sc = n3
    if(printable) write(nfout,'("supercell_on: uc_type,n1,n2,n3=",4i4)') uc_type,n1,n2,n3
  end subroutine m_CS_supercell_on

  subroutine m_CS_set_altv_prim
    altv = altv_prim
    rltv = rltv_prim
    univol = univol_prim
    rvol = rvol_prim
  end subroutine m_CS_set_altv_prim

  subroutine m_CS_set_altv_super
    altv = altv_super
    rltv = rltv_super
    univol = univol_super
    rvol = rvol_super
  end subroutine m_CS_set_altv_super

  subroutine m_CS_set_inv_sym_off()
    inversion_symmetry = OFF
  end subroutine m_CS_set_inv_sym_off

  subroutine m_CS_dealloc(neb_mode)
    logical, intent(in), optional :: neb_mode
    logical :: neb
    neb = .false.

! === KT_add === 13.1R
    if ( sw_keep_symmetry_strict == ON ) call m_CS_save_op_tau
! ============== 13.1R

    if(present(neb_mode)) neb = neb_mode
    if(allocated(op)) deallocate(op)
    if(allocated(tau)) deallocate(tau)
    if(.not.neb)then
       if(allocated(igen)) deallocate(igen)
       if(allocated(jgen)) deallocate(jgen)
       if(allocated(lpnt)) deallocate(lpnt)
       if(allocated(jgen_tl)) deallocate(jgen_tl)
       if(allocated(op_tl)) deallocate(op_tl)
       if(allocated(tau_tl)) deallocate(tau_tl)
    endif
  end subroutine m_CS_dealloc

! === KT_add == 13.1R
  subroutine m_CS_dealloc_op_tau
    if(allocated(op)) deallocate(op)
    if(allocated(tau)) deallocate(tau)
  end subroutine m_CS_dealloc_op_tau

  subroutine m_CS_save_op_tau
    if ( .not. allocated( op_previous ) )  allocate( op_previous( 3,3,nopr+af) )
    if ( .not. allocated( tau_previous ) ) allocate( tau_previous(3,nopr+af,CRDTYP) )
    nopr_previous = nopr;    op_previous = op;     tau_previous = tau
  end subroutine m_CS_save_op_tau

  subroutine m_CS_read_op_tau_previous
    integer :: iopr, i, j

    if ( sw_keep_symmetry_strict == OFF) return

    nopr = nopr_previous  
    if ( .not. allocated( op ) )  allocate( op( 3,3,nopr+af) )
    if ( .not. allocated( tau ) ) allocate( tau(3,nopr+af,CRDTYP) )

    op = op_previous;     tau = tau_previous

    do iopr=1, nopr +af
       do i=1, 3
          tau(i,iopr,CARTS) = altv(1,i)*tau(1,iopr,BUCS) &
               &     +      altv(2,i)*tau(2,iopr,BUCS) &
               &     +      altv(3,i)*tau(3,iopr,BUCS)
       enddo
    enddo

#if 0
    Do iopr=1, nopr
       write(800+mype,*) 'iopr = ', iopr
       Do i=1, 3
          write(800+mype,*) tau( i,iopr,CARTS), tau_previous( i,iopr, CARTS )
       End do
    End do
#endif

  end subroutine m_CS_read_op_tau_previous
! ============== 13.1R

! ======================== added by K. Tagami ==================== 12.0YAM
  subroutine check_if_sw_inversion_is_valid( nfout )
    integer, intent(in) :: nfout

    integer :: iopr

    real(kind=DP), parameter :: eps = 1.d-4
    real(kind=DP) :: c1

    logical :: Flag_inverse, sw_inversion_is_valid

    if ( inversion_symmetry == OFF ) return

    Flag_inverse = .false.
    sw_inversion_is_valid = .false.

    Do iopr=1, nopr
       if ( abs(op(1,1,iopr)+1.d0) < eps .and. &
            & abs(op(2,2,iopr)+1.d0) < eps .and. &
            & abs(op(3,3,iopr)+1.d0) < eps .and. &
            & abs(op(1,2,iopr)) < eps .and. &
            & abs(op(2,1,iopr)) < eps .and. &
            & abs(op(2,3,iopr)) < eps .and. &
            & abs(op(3,2,iopr)) < eps .and. &
            & abs(op(3,1,iopr)) < eps .and. &
            & abs(op(1,3,iopr)) < eps ) then

          Flag_inverse = .true.
       end if
!
       if ( Flag_inverse ) then
          c1 = abs( tau(1,iopr,BUCS) ) + abs( tau(2,iopr,BUCS) ) &
               &                       + abs( tau(3,iopr,BUCS) )
          if ( c1 < eps ) then
             sw_inversion_is_valid = .true.
          endif
          exit
       endif

    End do
! -------
    if ( sw_inversion_is_valid ) then
       write(nfout,*) '********************************************************'
       write(nfout,*) '******** Confirmed that sw_inversion is avalibale *****'
       write(nfout,*) '********************************************************'
    else
       if ( Flag_inverse ) then
          write(nfout,*) '********************************************************'
          write(nfout,*) '*   sw_inversion is invalid in this system !! '
          write(nfout,*) '*   Translation vector is non-zero for inversion operation'
          write(nfout,*) '********************************************************'
       else
          write(nfout,*) '********************************************************'
          write(nfout,*) '*   sw_inversion is invalid in this system !! '
          write(nfout,*) '*   Inversion operation does not exist'
          write(nfout,*) '********************************************************'
       endif
       call phase_error_with_msg(nfout,'sw_inversion is invalid',__LINE__,__FILE__)
    endif

  end subroutine check_if_sw_inversion_is_valid
! ================================================================ 12.0YAM

  subroutine m_CS_wd_BandSymInput(ipri,nfout,nf)
    integer, intent(in) :: ipri, nfout, nf
    character :: TF_inversion_symmetry ! {'T'|'F'}
    integer :: i, j
    integer :: id_sname = -1
    call tstatc0_begin('m_CS_wd_BandSymInput ',id_sname)
    if(mype == 0) then
       write(nf,'("##CRYSTAL_DATA")')
       write(nf,'("#LATTICE_SYSTEM")')
       if(inversion_symmetry == 0) TF_inversion_symmetry = 'F'
       if(inversion_symmetry == 1) TF_inversion_symmetry = 'T'
       write(nf,'(a,3x,a,3x,"! System{Rhombo|Hexa|Primitive|FCC|BCC|CBaseCentered},InversionSymmetry"/,"#")') &
            & Char_lattice_system(il), TF_inversion_symmetry

       write(nf,'("#GENERATORS")')
       write(nf,'(7(i0,1x))') (igen(i),(jgen(1,j,i),jgen(2,j,i),j=1,3),i=1,ngen);  write(nf,'("#")')

       write(nf,'("#BRAVAIS_LATTICE_PARAMETER"/,3f22.6," ! a,b,c"/,3f22.6," ! alpha,beta,gamma"/"#")') &
            & Bravais_lattice_length, Bravais_lattice_angle

       write(nf,'("#PRIMITIVE_LATTICE_PARAMETER",3(/,3f11.5,4x,3f11.5," ! rltv, altv"),/"#")') &
            &          ((rltv(i,j),j=1,3),(altv(i,j),j=1,3),i=1,3)
       write(nf,'("##"/)')
    end if
    call tstatc0_end(id_sname)
  end subroutine m_CS_wd_BandSymInput

end module m_Crystal_Structure
