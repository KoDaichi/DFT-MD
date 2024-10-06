!#define DEBUG_WRITE
!=======================================================================
!
!  PROGRAM  PHASE/0 2023.01
!
!  MODULE: m_Ionic_System
!
!  AUTHOR(S): T. Yamasaki  K, Betsuyaku, T. Uchiyama, Y. Morikawa,  August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, T. Yamamoto,  January/13/2004, April/10/2007
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
! =========================================
!   patch 0.1 by K. Tagami@adv    2009/10/19
!
!   patch 0.1:  correction for phonon with DFT+U
!=======================================================================
! Modified by T. Yamasaki, May/18/2023
!     A subroutine <<m_IS_md>> is revised to complete RIGID_BODY control functions.
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
! Functions:  [Identifier: 13.3B]
!                 amion is automatically set when it is not set by the user
!
! =============================================================
#ifdef forsafe
#define MASS_AUTOMATIC_SET
#endif

module m_Ionic_System
!     (m_IS)
!  $Id: m_Ionic_System.F90 633 2020-12-01 05:11:03Z jkoga $
!
!  This module is for structure factor, ewald energy,
!  and motions of atoms.
!  Original files are "md","md_ext","md_nebm1","md_nebm2",
!  "md_qn_gdiis", and "mdexpl", which were coded using VPP-fortran.
!  Subroutines contained in the files "md" and "md_ext" had been
!  translated into this module coded using f90+mpi by T. Yamasaki
!  and T. Uchiyama in 1999. Other files will be translated and be
!  contained in this module.
!
! Modified by T. Yamasaki, May/18/2023
!     RIGID_BODY control functions are revised.
!
! ------------------------------
! "md" --
!!$c     #1) imdtyp=10: structure optimization with a constraint
!!$c         on a bond length between two atoms.
!!$c         Y.Morikawa 24th July 1995.
!!$c     #2) imdtyp=8:  In a case that a force direction is parallel to
!!$c              a constraint vector (fcvect), the atom is not
!!$c              constraint the movement.
!!$c         T. Yamasaki 12th Jan. 1996.
!!$c     #3) imdtyp=101: quenched motion of atoms, of which
!!$c                      center-of-mass is constrained.
!!$c         T. Uchiyama, November 9, 1997.
! ------------------------------
! "md_ext"
!!$c This file contains three md subroutines coded by T. Uchiyama
!!$c (JRCAT-ATP) in 1997-1999.
!!$c   1. md_thermo
!!$c   2. md_bluem
!!$c   3. md_cnstr
!!$c******************************************************************
!!$      subroutine md_thermo
!!$c For imdalg= -1 ; Constant temperature molecular dynamics.
!!$c   imdtyp=    0 ; fixed,         imdtyp=    1 ; free,
!!$c   imdtyp= 1001 ; thermostat#1,  imdtyp= 1002 ; thermostat#2,
!!$c    ....
!!$c   imdtyp= 1000+nrsv, where nrsv= No. of thermostats.
!!$c
!!$c  Ref). S. Nose, Mol. Phys. 52 (1984) 255,
!!$c                 J. Chem. Phys. 81 (1984) 511.
!!$c
!!$c     Released.    Dec. 12, 1997,   by T. Uchiyama
!!$c     Bug fixed.   Aug.  6, 1998,   by T. Uchiyama
!!$c******************************************************************
!!$      subroutine md_bluem
!!$c   nrsv=      1 ; No. of thermostat
!!$c   imdtyp=    0 ; fixed
!!$c   imdtyp= 1001 ; thermostat
!!$c
!!$c     Released.      Sep. 21, 1998.   by T. Uchiyama
!!$c     Bug fixed.     Nov.  3, 1998.   by T. Uchiyama
!!$c******************************************************************
!!$      subroutine md_cnstr
!!$c   imdtyp= 0 ; fixed
!!$c   imdtyp= 1 ; quenched
!!$c
!!$c     Released.      Sep. 21, 1998,   by T. Uchiyama
!!$c     Bug fixed.     Nov.  4, 1998,   by T. Uchiyama
!!$c******************************************************************
! "md_nebm1"
!!$C Nudged Elastic Band method after Hannes Jonsson.
!!$c                           29th Feb. 1996 Y.M
!!$c                           @(#)md_nebm1.f 9.2 97/12/08 11:42:46
!!$c  1) open command is commented out by T. Yamasaki. 22nd May. 1996
! "md_nebm2"
!!$C (1) Nudged Elastic Band method after Hannes Jonsson.
!!$c                           29th Feb. 1996 Y.M
!!$c                   revised 14th Oct. 1996 Y.M
!!$c (2)               Aug. 1996 T. Yamasaki
!!$c    cpd --> cpd_l, forc --> forc_l
!!$c
!!$c (3)               revised 14th Oct. 1996 Y.M
!!$*                           @(#)md_nebm2.f 9.2 97/11/21 17:32:00

  use m_Crystal_Structure,  only : nopr, tau, op, altv, rltv, c,univol, nbztyp &
       &                         , inversion_symmetry, p2bmat & ! inverse transformation matrix
       &                         , b2pmat &
       &                         , sw_supercell, n1_sc, n2_sc, n3_sc &
       &                         , nlpnt, lpnt, altv_prim, m_CS_set_inv_sym_off &
       &                         , symmetry_check_criterion, ngen_tl,tau_tl,op_tl &
       &                         , m_CS_altv_2_rltv, sw_supercell_symmetry
  use m_Files,              only : nfout, nfmode, m_Files_open_nfmode &
       &                         , nfimage, m_Files_open_nfimage, nfdynm_cif &
       &                         , m_Files_set_def_fname_pos, m_Files_set_def_fname_attr &
       &                         , m_Files_rd_file_names_data &
       &                         , m_Files_open_nfpos, m_Files_close_nfpos, nfpos &
       &                         , F_INP, F_DYNM, F_POS, F_COORD_ATTR &
       &                         , m_Files_open, m_Files_close, F_IMAGE, nfpath, F_PATH &
       &                         , m_Files_open_xsf, m_Files_close_xsf, nfxsf &
       &                         , m_Files_open_n2p2, m_Files_close_n2p2, nfn2p2 &
       &                         , m_Files_open_file, m_Files_close_file
  use m_Timing,             only : tstatc0_begin, tstatc0_end
  use m_Control_Parameters, only : ipri,ipriinputfile, iprimd, dtio,icond, ipriforce &
       &                         , forccr, istress, kimg, af, iprigdiis,ipristrcfctr &
       &                         , kqnmditer_p, gdiis_hownew, c_forc2GDIIS &
       &                         , c_forc_prop_region_high, c_forc_prop_region_low &
       &                         , factor_prop_region, imdalg &
       &                         , sw_calc_force, sw_displace_atom &
       &                         , sw_calc_force_all, sw_phonon_oneshot &
       &                         , sw_vibrational_mode, with_mode_effchg &
       &                         , sw_vibrational_modes &
       &                         , sw_polynomial_fit, num_phonon_calc_mode, norder &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!      &                         , sw_screening_correction, sw_pair_vdw &
       &                         , sw_screening_correction, sw_vdw_correction &
! ==============================================================================
       &                         , dftd3_damping_function &
       &                         , etol,printable, m_CtrlP_what_is_mdalg &
       &                         , mode_fi_coefficient, fi_coefficient &
       &                         , multiple_replica_mode &
       &                         , sw_correct_eigenvalue, eigenvalue_threshold &
       &                         , sw_extrapolate_charge,rms_threshold, sw_wf_predictor &
       &                         , sw_optimize_alpha, vdw_method, xctype, iprivdw &
       &                         , fire_incre_factor, fire_decre_factor, fire_dtmax &
       &                         , fire_nmin, fire_invmass_factor, fire_initial_dt &
       &                         , fire_decre_factor_alpha, fire_alpha_start, iprifire &
       &                         , sw_stress_correction &
       &                         , sw_calc_wf_atom_decomposition &
       &                         , sw_prec, precon_mu, precon_A, maxstep &
       &                         , m_CtrlP_rd_val &
#ifndef _EMPIRICAL_
       &                         , ipriphonon &
       &                         , m_CtrlP_check_matm
#else
       &                         , ipriphonon
#endif
#ifndef _EMPIRICAL_
  use m_Control_Parameters, only : sw_aldos, sw_pdos, sw_orb_popu, sw_hubbard &
       &                         , proj_attribute, proj_group &
       &                         , num_projectors &
       &                         , num_proj_group, num_proj_elems &
       &                         , m_CtrlP_set_proj_ityp, m_CtrlP_reset_optmode &
#ifdef ENABLE_ESM_PACK
       &                         , driver,sw_optimize_lattice,ipripredictor,sw_esm
#else
       &                         , driver,sw_optimize_lattice,ipripredictor
#endif
#endif
  use m_IterationNumbers,   only : iteration_ionic, m_Iter_mdIterN_increment, iteration_unit_cell
  use m_Const_Parameters,   only : DP, FIX, RELAX, VERLET, QUENCHED_MD,CG_STROPT,SD_MD &
       &                         , FIX_IN_A_PLANE, FIXED_NORMAL_HYPERVECTOR &
       &                         , RIGID_BODY, RIGID_BODY_FIX_IN_A_PLANE &
       &                         , COG_FIX, COG_CNTR, BONDLENGTH_FIX, BONDLENGTH_FIX_1 &
       &                         , COG_and_RIGID_BODY_FIX_L &
       &                         , BONDLENGTH_FIX_2, COG_FIX_L, PAI2, PAI &
       &                         , ORDINA, CNSTRA, FIX_HBATH, RELAX_HBATH &
       &                         , HEXAGONAL,ORTHORHOMBIC,HEX1FOLD, GDIIS, BFGS, L_BFGS, CARTS &
       &                         , INITIAL, DELTA, DELTA10 &
       &                         , BLUEMOON, QUENCHED_CONSTRAINT, FIRE, T_CONTROL &
       &                         , ANEW, RENEW, HEAT_BATH, PUCV, FMAXVALLEN, FMAXUNITLEN &
       &                         , NOCONV, LOWER, UPPER, NOSE, NOSE_HOOVER &
       &                         , SmallestPositiveNumber, ON, OFF, YES, NO &
       &                         , unit_conv_byname, DIRECTIN, FROM_ENDPOINTS, FILE &
       &                         , STEEPEST_DESCENT, PROPORTIONAL, CONST_kB, VELOCITY_SCALING &
       &                         , CONTINUATION, DRIVER_MTD, COORDINATE_CONTINUATION, AUTOMATIC &
       &                         , CG_STROPT2,BOHR, PT_CONTROL, P_CONTROL &
       &                         , GaussLegendre, SphericalHarmonicsExpansion &
       &                         , CYLINDER, BOX &
       &                         , VDW_DFTD3 &
       &                         , PHASE0_INPUT, PHASE0_OUTPUT, NEB_INPUT &
       &                         , HARTREE, BOHR, LANGEVIN &
       &                         , VOLUME, METRIC_TENSOR, LATTICE_VECTOR

  use m_Parallelization,    only : MPI_CommGroup,ista_kngp, iend_kngp, npes,mype,ierr, ista_atm, iend_atm, ista_atm2, iend_atm2 &
       &                         , mype_conf, nrank_conf, conf_para
  use m_Parallelization,    only : mpi_ke_world, mpi_chg_world

! =============================== added by K. Tagami =================== 11.0
  use m_Control_Parameters,     only : noncol, SpinOrbit_mode, sw_calc_core_energy
  use m_Const_Parameters,       only : ByProjector, ByPawPot, ZeffApprox, ReadFromPP
! ====================================================================== 11.0

! === KT_add ==== 2014/07/14
  use m_Files,  only : nfhypervec, m_Files_open_nfhypervec, m_Files_close_nfhypervec
! =============== 2014/07/14
! === KT_add === 13.2S
  use m_Crystal_Structure,  only : sw_use_magnetic_symmetry, &
       &                           sw_spinorbit_second_variation, &
       &                           sw_fix_global_quantz_axis, &
       &                           Global_Quantz_Axis_Fixed
! ============== 13.2S

  use m_Files, only : m_Files_open_dftd3par,m_Files_close_dftd3par,nfdftd3par

  use m_Control_Parameters, only : filetype_nnp
  use m_Const_Parameters,   only : N2P2, XSF, DEEPMD, ALL
#ifdef __EDA__
  use m_Control_Parameters, only : sw_eda
#endif

  use m_Control_Parameters, only : max_torque, max_force_trans
  use mpi

  implicit none
!  include 'mpif.h'

#ifdef PGI
  external ibath
  integer :: ibath
#endif

#ifdef DEBUG_WRITE
  integer, parameter :: DEBUGPRINTLEVEL = 1
#else
  integer, parameter :: DEBUGPRINTLEVEL = 2
#endif

  integer :: input_coordinate_system = PUCV       ! ncord

  integer ::                                      ntyp = 1
  integer ::                                      natm = 1, natm2 = 1
  integer ::                                      natm_prim = 1
  integer ::                                      natm2_prim = 1
  integer ::                                      natm_super = 1
  integer ::                                      natm2_super = 1
  integer ::                                      natmorg = 1
  integer ::                                      natm_mobile = 1

  character(len("Ionic System")),private,parameter :: tag_ionic_system = "Ionic System"
  character(len("Ionic System Attributes")),private,parameter ::  &
      & tag_ionic_system_attributes = "Ionic System Attributes"
  character(len("-- forcp --")),private,parameter  :: tag_forcp        = "-- forcp --"

  character(len("distances_of_planes")),private,parameter :: tag_distances_of_planes = "distances of planes"

  real(kind=DP), allocatable,dimension(:,:)   ::  pos, cps  ! d(natm,3)
  real(kind=DP), allocatable,dimension(:,:)   ::  oldcps    ! d(natm,3)
  real(kind=DP), allocatable,dimension(:,:)   ::  pos_in, cps_in
  ! pos : atomic coordinates in a unit cell vector system
  ! cps : atomic coordinates in the cartesian system
  real(kind=DP), allocatable,dimension(:,:)   ::  pos_end0, pos_end1, cps_end0, cps_end1 ! d(natm,3)
  real(kind=DP), allocatable,dimension(:,:,:) ::  pos_image, cps_image ! d(:,natm,3)
  real(kind=DP), allocatable,dimension(:,:,:) ::  normal_hypervector ! d(natm,3,PUCV:CARTS)
  real(kind=DP), allocatable,dimension(:,:)   ::  pos_prim, cps_prim  ! d(natm_prim,3)
  real(kind=DP), allocatable,dimension(:,:)   ::  pos_super, cps_super  ! d(natm_super,3)
  real(kind=DP), allocatable,dimension(:)     ::  ionic_mass ! d(natm)
  real(kind=DP), allocatable,dimension(:,:)   ::  cpd_l     ! d(natm,3)
  real(kind=DP), allocatable,dimension(:,:,:) ::  cpo_l     ! d(natm,3,3)
  real(kind=DP), allocatable,dimension(:,:,:) ::  cps_history
  integer,private :: ncps_history=0
  integer,private,parameter ::                    LEN_ATOMNAME = 4
  character(len=LEN_ATOMNAME),allocatable, dimension(:) ::  speciesname ! d(ntyp)
  character(len=LEN_ATOMNAME),private,allocatable,dimension(:) :: species_work ! d(natm)
  character(len=LEN_ATOMNAME),private,allocatable,dimension(:) :: species_indp ! d(natm)

  real(kind=DP), allocatable, dimension(:) ::     iatomn, iatom      ! d(ntyp)
#ifdef __EDA__
! -----  ascat starts modifying  -----
  integer,       allocatable, dimension(:) ::     iatomn_EDA         ! d(natm)
  integer, allocatable, dimension(:)       :: icount_EDA             ! d(ntyp)
  real(kind=DP),public,target,allocatable,dimension(:) :: eewald_per_atom    ! d(natm)
  real(kind=DP),private,allocatable,dimension(:) :: eewald_per_atom_Rspace
  real(kind=DP),private,allocatable,dimension(:) :: eewald_per_atom_Gspace
  real(kind=DP), allocatable, dimension(:,:,:) :: zfm3_l_EDA !d(ista_kngp:iend_kngp,ntyp,kimg)
! -----  ascat ceases modifying  -----
#endif
  integer, allocatable, dimension(:) ::           ivan               ! d(ntyp)
  real(kind=DP),allocatable,dimension(:) ::       alfa, amion, zeta1, qex ! d(ntyp)
  integer, allocatable, dimension(:) ::           iwei, imdtyp, ityp ! d(natm)
  integer, allocatable, dimension(:,:) ::         imdtypxyz
#ifndef _EMPIRICAL_
  integer, allocatable, dimension(:) ::           if_pdos, if_aldos ! d(natm)
#endif
  integer, allocatable, dimension(:) ::           numlay             ! d(natm)

  integer ::                                      nfcatm = 0
  integer ::                                      nfcatm_cog = 0, nfcatm_rb = 0

  real(kind=DP) ::                                ekina
  integer ::                                      mdmode = ORDINA

  integer,       allocatable, dimension(:,:)   :: napt   !d(natm,nopr+af)
  integer,       allocatable, dimension(:,:)   :: napt_tl!d(natm,ngen_tl)
  integer,       allocatable, dimension(:,:)   :: napt_prim !d(natm_prim,nopr+af)

  integer                                      :: nopr_supercell
  integer,       allocatable, dimension(:,:)   :: napt_supercell !d(natm,nopr_supercell_red)
  integer,       allocatable, dimension(:)     :: iop_supercell  !d(nopr_supercell)
!!$  real(kind=dP), allocatable, dimension(:,:)   :: tau_supercell  !d(3,nopr_supercell)
!!$  integer                                      :: mnope_supercell !=maxval(nope_supercell(:))
!!$  integer,       allocatable, dimension(:)     :: nope_supercell !d(nopr)
!!$  integer,       allocatable, dimension(:,:)   :: pope_supercell !d(mnope_supercell,nopr)

  real(kind=DP), allocatable, dimension(:,:,:) :: zfm3_l !d(ista_kngp:iend_kngp,ntyp,kimg)
  real(kind=DP)                                :: eewald
  real(kind=DP), allocatable, dimension(:,:)   :: fxyzew_l ! d(natm,3)
  real(kind=DP),private,pointer,dimension(:,:) :: cpd_old  ! d(natm,3)
  integer,      private,allocatable,dimension(:)   :: ipcpd    ! d(natm2)
  real(kind=DP), dimension(3,3)                :: s_ew
  integer                                      :: displaced_atom = 0
  real(kind=DP), dimension(3)                  :: displacement(3)

  character(len=9), public       :: lattice_system_from_m_CS_SG = 'none'

! -- Temperature Control
  integer         ::                                  t_ctrl_method = NOSE_HOOVER

  integer         ::                                  nrsv = 1  ! number of heat bath
  integer         ::                                  set_initial_velocity = ON  ! (by J. Koga)
  integer         ::                                  sw_read_velocities = OFF  ! (by T. Yamamoto)
  integer         ::                                  sw_shift_velocities = OFF
  integer         ::                                  nchain = 1
  real(kind=DP),private ::                            tk_initial = 0.d0
!!$  real(kind=DP),private,allocatable,dimension(:)   :: qmass,tkb,cprv,frsv ! d(nrsv)
  real(kind=DP),allocatable,dimension(:)   :: qmass,tkb,cprv,frsv ! d(nrsv)
  integer, allocatable,dimension(:)        :: natm_per_thermo
  real(kind=DP),private,allocatable,dimension(:,:) :: cpqr                ! d(nrsv,2)
  real(kind=DP),private,allocatable,dimension(:,:) :: forcp               ! d(natm,3)
  real(kind=DP),private,allocatable,dimension(:,:) :: qmass_c,cprv_c,cpqr_c
  real(kind=DP),private                            :: qfactor

  real(kind=DP),private                            :: ekinq,ekbt
  real(kind=DP)                                    :: ega,almda
  real(kind=DP)                                    :: forcmx_constraint_quench = 1.d+2
! -->   T. Yamasaki 18 July 2008
  real(kind=DP),private                            :: forc_norm_hyperplane_vert = 1.d+2
  real(kind=DP)                                    :: forcmx_hyperplane_vert = 1.d+2
! <--
  real(kind=DP),private,allocatable,dimension(:,:) :: gca                 ! d(natm,3)
  integer,      private,allocatable,dimension(:)   :: ia_cnst             ! d(nfcatm) , ia_cnst(ip): pointer to whole atoms
  integer,      private,allocatable,dimension(:)   :: icnst_a             ! d(natm)   , icnst_a(ia): pointer to constrained atoms
  !  (An example)
  !                   ia      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  !               imdtyp    201 201 201 201   1   1   1   1   0   0   0   0   1   1   1   1 201 201 201 201
  !         nfcatm = 8
  !                    i      1   2   3   4                                                   5   6   7   8
  !             ia_cnst(i)    1   2   3   4                                                  17  18  19  20
  !
  !                   ia      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  !             icnst_a(ia)   1   2   3   4   0   0   0   0   0   0   0   0   0   0   0   0   5   6   7   8
  !
  real(kind=DP),private,allocatable,dimension(:,:) :: fcvect              ! d(nfcatm,8)
  real(kind=DP),private,allocatable,dimension(:)   :: rigid_body_vect     ! d(4)
  real(kind=DP)                                    :: max_reach_of_fixed_plane
  real(kind=DP)                                    :: moved_distance_of_fixed_plane = 0.d0

!!$  real(kind=DP),private,allocatable,dimension(:,:) :: fcvect_tmp          ! d(num_planes_atoms_are_fixed,4)
  integer,      private,allocatable,dimension(:)   :: ipfixedplane        ! d(nfcatm)
  integer,      private,allocatable,dimension(:)   :: icount_of_ipfixedplane      ! d(nfcatm)
  integer,      private,allocatable,dimension(:)   :: relax_in_fixedplane     ! d(nfcatm)
  real(kind=DP),private,dimension(4)               :: sgmc
  character(len("forcmx_constraint_quench")),private,parameter :: &
       &                              tag_forcmx_const = "forcmx_constraint_quench"
  character(len("structure_evolution")),private,parameter :: &
       &                              tag_structure_evolution = "structure_evolution"
  character(len("temperature_control")),private,parameter :: &
       &                              tag_temperature_control = "temperature_control"
  character(len("method")),private,parameter ::       tag_method          = "method"
  character(len("num_thermostat")),private,parameter :: tag_num_thermostat = "num_thermostat"
  character(len("num_thermo")),private,parameter :: tag_num_thermo  = "num_thermo"
  character(len("num_chain")),private,parameter :: tag_num_chain  = "num_chain"
  character(len("thermostat")),private,parameter :: tag_thermostat  = "thermostat"
  character(*),private,parameter ::                 tag_weight_thermo = "weight_thermo"
  character(*),private,parameter ::                 tag_weight_thermostat = "weight_thermostat"
  character(len("weight")),private,parameter ::     tag_weight = "weight"
  character(len("qmass")),private,parameter ::      tag_qmass  = "qmass"
  character(len("temperature")),private,parameter :: tag_temperature = "temperature"
  character(len("temp")),private,parameter ::       tag_temp   = "temp"
  character(len("T")),private,parameter ::          tag_T      = "T"
  character(len("NOSE_HOOVER")),private,parameter :: tag_Nose_Hoover = "nose_hoover"
  character(len("HOOVER")),private,parameter ::     tag_Hoover = "hoover"
  character(len("Nose")),private,parameter ::       tag_Nose   = "nose"
  character(len("velocity_scaling")),private,parameter :: tag_velocity_scaling = "velocity_scaling"
  character(len("langevin")), private, parameter       :: tag_langevin = "langevin"
  character(len("iseed")), private, parameter          :: tag_iseed = "iseed"
  character(len("tdamp")),private,parameter :: tag_tdamp = "tdamp"
  character(len("tdamp_baro")),private,parameter :: tag_tdamp_baro = "tdamp_baro"

  character(len("lattice_vector")),private,parameter :: tag_lattice_vector = "lattice_vector"

! === KT_add === 2014/06/10
  character(len("duration")),private,parameter ::      tag_duration  = "duration"
  character(len("temp_s")),private,parameter ::        tag_temp_s  = "temp_s"
  character(len("temp_e")),private,parameter ::        tag_temp_e  = "temp_e"
  character(len("sw_change_temperature_by_step")), private,parameter :: &
       &        tag_sw_change_temp_bystep  = "sw_change_temperature_by_step"

  integer, allocatable :: mdstep_at_start_thermostat(:)
  integer, allocatable :: mdstep_at_end_thermostat(:)
  real(kind=DP), allocatable :: temp_at_start_thermostat(:)
  real(kind=DP), allocatable :: temp_at_end_thermostat(:)

  integer :: sw_change_temperature_by_step = OFF
  integer :: ithermo_now = 1
! ============== 2014/06/10

  logical ::            tag_T_cntrl_is_found
  character(len("Temperature Control")),private,parameter :: tag_T_cntrl = "Temperature Control"
  character(len("temperature_control")),private,parameter :: tag_T_cntrl2 = "temperature_control"
  character(len("nrsv")),private,parameter ::                tag_nrsv = "nrsv"
  character(len("heat bath")),private,parameter ::           tag_heat_bath = "heat bath"
  character(len("atom")),private,parameter ::                tag_atom_velocity = "atom"
  character(len("set_initial_velocity")),private,parameter :: tag_set_initial_velocity = "set_initial_velocity" ! (by J. Koga)
  character(len("sw_read_velocities")),private,parameter :: tag_sw_read_velocities = "sw_read_velocities"
  character(len("initial_temperature")),private,parameter :: tag_initial_temperature = "initial_temperature"
  character(len("sw_shift_velocities")),private,parameter :: tag_sw_shift_velocities = "sw_shift_velocities"

  !     constraint             sigma                       sgmc
  ! 1. BONDLENGTH_FIX_1
  !       (sigma)        $|\vec{r_1} - \vec{r_2}| - b $
  !       (sgmc)              sgmc(1) = b
  ! 2. BONDLENGTH_FIX_2
  !       (sigma)        $|\vec{r_1} - \vec{r_2}|^2 - b $
  !       (sgmc)              sgmc(1) = b
  ! 3. COG_FIX_L
  !       (sigma)        $(\vec{r_c} - \vec{r_o})\cdot\vec{a} - b$
  !                where,
  !                  $\vec{r_c} = (m_1\vec{r_1}+m_2\vec{r_2})/m_g$
  !                  $\vec{r_o} = \sum_{All except fixed}m_i \vec{r_i}/m_r
  !                            m_g = (m_t m_r)/(m_t + m_r)
  !                            m_r = \sum_{All except fixed} m_i
  !                            m_t = m_1 + m_2
  !       (sgmc)            sgmc(1:3) = \vec{a}
  !                           sgmc(4) = b

  integer                                  :: cnst_typ

! -- Work_arrays for GDIIS
!!$  integer, private, parameter :: kqnmditer_p  = 4
  real(kind=DP),private,allocatable,dimension(:,:,:):: u_l ! d(natm,3,kqnmditer_p)
  real(kind=DP),private,allocatable,dimension(:,:,:):: w_l ! d(natm,3,kqnmditer_p)

!! for continuation
  real(kind=DP),private,allocatable,dimension(:,:,:):: u_l_buf ! d(natm,3,kqnmditer_p)
  real(kind=DP),private,allocatable,dimension(:,:,:):: w_l_buf ! d(natm,3,kqnmditer_p)
  integer, private, allocatable, dimension(:) ::       ncrspd_buf ! d(kqnmditer_p)
  logical, private :: diis_continuable = .false.

  real(kind=DP),private,allocatable,dimension(:,:) ::  fc_l ! d(natm,3)
  integer, allocatable, dimension(:)               ::  ncrspd ! d(kqnmditer_p)
  real(kind=DP),private,allocatable,dimension(:,:) ::  f_gdiis ! d(kqnmditer_p,kqnmditer_p)
  real(kind=DP),private,allocatable,dimension(:,:) ::  f_wk ! d(kqnmditer_p**2,2)
  real(kind=DP),private,allocatable,dimension(:)   ::  f_rslv ! d(kqnmditer_p**2)
  real(kind=DP),private,allocatable,dimension(:) ::    g,e_wk,ww1,etot_trial
  real(kind=DP),private,allocatable,dimension(:,:) ::  forc_g ! d(natm,3)
  integer, private, allocatable, dimension(:) ::       ip ! d(kqnmditer_p)
  integer, private, parameter ::                       UNIT = 1
!!$  real(kind=DP),private,parameter ::              c_forc_prop_region_high = 0.1d0
!!$  real(kind=DP),private,parameter ::            c_forc_prop_region_low = 0.0001d0
!!$  real(kind=DP),private,parameter ::                   factor_prop_region = 0.02d0
!!$  real(kind=DP),private,parameter ::            c_forc_QUENCH2GDIIS = 0.0008d0
  integer,private,parameter ::                         ic_E_overshoot = 3
  real(kind=DP),private,save ::                        etot_previous = 1.d+3
  integer,private,save ::                              iincre_at_forc_cal = 0
!  integer,private,save ::                              iter_gdiis = 0
  integer,save ::                                      iter_gdiis = 0
  integer,save ::                                      iter_md = 0
  integer,private,save ::                              if_allocated=0
  real(kind=DP) :: min_alpha = 1.d0,max_alpha = 2.0d0
!  real(kind=DP) :: wolfe_c1=1.d-4,wolfe_c2=0.9
  real(kind=DP) :: wolfe_c1=0.1d0,wolfe_c2=0.9
! -- Work arrays
  real(kind=DP),private,pointer, dimension(:)     :: ekr                 ! d(nrsv)
  integer,      private,pointer, dimension(:)     :: nathm               ! d(nrsv)

  integer,      parameter   :: len_str = 132
  character(len=len_str)       :: str

  ! ---> input tag terms
!  integer, private,save ::                          howtogive_coordinates = DIRECTIN ! {DIRECTIN|FROM_ENDPOINTS}
  integer, private,save ::                          howtogive_coordinates = FROM_ENDPOINTS ! {DIRECTIN|FROM_ENDPOINTS}
  integer, private,save ::                          endpoint_images  = NO            ! {NO|FILE|DIRECTIN}
  ! --- Structure ---
  character(len("structure")),private ::            tag_structure = "structure"
  character(len("atom_list")),private,parameter ::  tag_atomlist = "atom_list"
  character(len("howtogive_coordinates")),private,parameter :: tag_howtogive_coordinates = "howtogive_coordinates"
  character(len("from_endpoint_images")),private,parameter :: tag_from_endpoint_images = "from_endpoint_images"
  character(len("from_endpoints")),private,parameter :: tag_from_endpoints = "from_endpoints"
  character(len("atom_list_end0")),private,parameter :: tag_atom_list_end0 = "atom_list_end0"
  character(len("atom_list_end1")),private,parameter :: tag_atom_list_end1 = "atom_list_end1"
  character(len("endpoint_images")),private,parameter :: tag_endpoint_images = "endpoint_images"
  character(len("directin")),private,parameter ::   tag_directin = "directin"
  character(len("file")),private,parameter ::       tag_file     = "file"
  character(len("nothing")),private,parameter ::    tag_nothing  = "nothing"

  character(len("atom_duplication")),private,parameter ::tag_atomduplication = "atom_duplication"
  character(len("symmetry")),private,parameter ::   tag_ad_symmetry = "symmetry"
  character(len("num_atoms")),private,parameter ::  tag_numatoms     = "num_atoms"
  character(len("coordinate_system")),private,parameter :: tag_coordinate_system = "coordinate_system"
  character(len("cartesian")),private,parameter ::  tag_cartesian = "cartesian"
  character(len("id")),private,parameter ::         tag_id        = "id"
  character(len("no")),private,parameter ::         tag_no        = "no"
  character(len("number")),private,parameter ::     tag_number    = "number"
  character(len("element")),private,parameter ::    tag_element   = "element"
  character(len("XYZ")),private,parameter ::        tag_XYZ       = "xyz"
  character(len("pucv")),private,parameter ::       tag_pucv      = "pucv"
  character(len("internal")),private,parameter ::   tag_internal  = "internal"
  character(len("relative")),private,parameter ::   tag_relative  = "relative"
  character(len("atoms")),private,parameter ::      tag_atoms     = "atoms"
  character(len("rx")),private,parameter ::         tag_rx        = "rx"
  character(len("ry")),private,parameter ::         tag_ry        = "ry"
  character(len("rz")),private,parameter ::         tag_rz        = "rz"
  character(len("vx")),private,parameter ::         tag_vx        = "vx"
  character(len("vy")),private,parameter ::         tag_vy        = "vy"
  character(len("vz")),private,parameter ::         tag_vz        = "vz"
  character(len("fx")),private,parameter ::         tag_fx        = "fx"
  character(len("fy")),private,parameter ::         tag_fy        = "fy"
  character(len("fz")),private,parameter ::         tag_fz        = "fz"
  character(len("normx")),private,parameter ::      tag_normx     = "normx"
  character(len("normy")),private,parameter ::      tag_normy     = "normy"
  character(len("normz")),private,parameter ::      tag_normz     = "normz"
  character(len("increment")),private,parameter ::  tag_increment = "increment"
  character(len("incrvx")),private,parameter ::     tag_incrvx    = "incrvx"
  character(len("incrvy")),private,parameter ::     tag_incrvy    = "incrvy"
  character(len("incrvz")),private,parameter ::     tag_incrvz    = "incrvz"
  character(len("max_reach")),private,parameter::   tag_max_reach = "max_reach"
  character(len("final_value")),private,parameter:: tag_final_value = ",final_value"
  character(len("iplane")),private,parameter ::     tag_iplane    = "iplane"
  character(len("relax_in_fixedplane")),private,parameter :: tag_relax_in_fixedplane = "relax_in_fixedplane"
  character(len("mobile")),private,parameter ::     tag_mobile    = "mobile"
  character(len("mobilex")),private,parameter ::     tag_mobilex    = "mobilex"
  character(len("mobiley")),private,parameter ::     tag_mobiley    = "mobiley"
  character(len("mobilez")),private,parameter ::     tag_mobilez    = "mobilez"
  character(len("weight")),private,parameter ::     tag_a_weight  = "weight"
  character(len("pdos")),private,parameter ::       tag_pdos      = "pdos"
  character(len("aldos")),private,parameter ::      tag_aldos     = "aldos"
  character(len("atom_decomp")),private,parameter ::  tag_atom_decomp = "atom_decomp"
  character(len("thermo_group")),private,parameter :: tag_thermo_group = "thermo_group"
  character(len("thermo_g")),private,parameter ::   tag_thermo_g = "thermo_g"
  character(len("num_layer")),private,parameter ::  tag_num_layer = "num_layer"
  character(len("displacement")),private,parameter :: tag_displacement = "displacement"
  character(len("sw_displace_atom")),private,parameter :: tag_sw_displace_atom = "sw_displace_atom"
  character(len("displaced_atom")),private,parameter :: tag_displaced_atom = "displaced_atom"

  character(len("key")),private,parameter ::  tag_key = "key"
  character(len("atom_key")),private,parameter ::  tag_atom_key = "atom_key"
  integer, allocatable, dimension(:) ::           atom_key             ! d(natm)

! === KT_add === 2015/03/14
  character(len("sw_displacement_in_carts")), private,parameter :: &
       &              tag_sw_displacement_in_carts = "sw_displacement_in_carts"
  integer :: sw_displacement_in_carts = off
! ============== 2015/03/14
  character(len("ux")),private,parameter :: tag_ux = "ux"
  character(len("uy")),private,parameter :: tag_uy = "uy"
  character(len("uz")),private,parameter :: tag_uz = "uz"
  character(len("vibrational_mode")),private,parameter :: tag_vibrational_mode = "vibrational_mode"
  character(len("sw_vibrational_mode")),private,parameter :: tag_sw_vibrational_mode = "sw_vibrational_mode"
  character(len("mode_index")),private,parameter :: tag_mode_index = "mode_index"
  character(len("normal_coordinate")),private,parameter :: tag_normal_coordinate = "normal_coordinate"
  character(len("with_mode_effchg")),private,parameter :: tag_with_mode_effchg = "with_mode_effchg"

  character(len("mobility_by_distance")),private,parameter :: tag_mobility_by_distance = "mobility_by_distance"
  character(len("sw_mobility_by_distance")),private,parameter :: tag_sw_mobility_by_distance &
                                                               &  = "sw_mobility_by_distance"
  character(len("target_atom")),private,parameter :: tag_target_atom="target_atom"
  character(len("target_posx")),private,parameter :: tag_target_posx="target_posx"
  character(len("target_posy")),private,parameter :: tag_target_posy="target_posy"
  character(len("target_posz")),private,parameter :: tag_target_posz="target_posz"
  character(len("distance")),private,parameter :: tag_distance="distance"
  integer :: sw_mobility_by_distance = OFF
  integer :: target_atom = 0
  real(kind=DP), dimension(3) :: target_pos=0.d0
  real(kind=DP) :: distance = 0.d0

  logical,public ::                                 constraints_exist = .false.
  logical,public ::                                 move_constrained_plane = .false.
  integer,public ::                                 iteration_ionic_at_CNSTRA = -1
  integer,private ::                                iteration_ionic_in_constraint = -1
  integer,private ::                                constraint_type   = 0
  integer,private ::                                num_planes_atoms_are_fixed = 1
  integer,private ::                                num_planes_atoms_are_fixed_rb = 0, num_planes_atoms_are_fixed_cog = 0
  integer,private ::                                num_fixed_bonds  = 0
  integer,private,allocatable,dimension(:,:)  ::    bondlength_fix_set  ! (1:2,num_bonds)
!!$  real(kind=DP),private :: fcg_cog(3), fcg_rb(3), fcg_mdfy_cog(3), fcg_mdfy_rb(3)
  real(kind=DP),private,allocatable,dimension(:,:) :: fcg_cog,    fcg_mdfy_cog    ! d(num_planes_atoms_are_fixed_cog,1:3)
  real(kind=DP),private,allocatable,dimension(:,:) :: fcg_rb, fcg_mdfy_rb ! d(num_planes_atoms_are_fixed_rb,1:3)
!!$  real(kind=DP),private ::  tmass = 0.d0, rtmass = 0.d0 ! rtmass = 1/tmass
!!$  real(kind=DP),private ::  tmass_rb = 0.d0, rtmass_rb = 0.d0 ! rtmass_rb = 1/tmass_rb
  real(kind=DP),private,allocatable,dimension(:)  ::  tmass, rtmass         ! d(num_planes_atoms_are_fixed)
  real(kind=DP),private,allocatable,dimension(:)  ::  tmass_cog, rtmass_cog ! rtmass_cog = 1/tmass_cog d(num_planes_atoms_are_fixed_cog)
  real(kind=DP),private,allocatable,dimension(:)  ::  tmass_rb,  rtmass_rb  ! rtmass_rb = 1/tmass_rb d(num_planes_atoms_are_fixed_rb)
  real(kind=DP),private,allocatable,dimension(:)  ::  distance_cog          ! d(num_planes_atoms_are_fixed_cog)
  real(kind=DP),private,allocatable,dimension(:)  ::  distance_rb           ! d(num_planes_atoms_are_fixed_rb)

  character(len("constraint")),private,parameter :: tag_constraint = "constraint"
  character(len("num_fixed_bonds")),private,parameter:: tag_num_fixed_bonds = "num_fixed_bonds"
  character(len("fixed_bond")),private,parameter :: tag_fixed_bond = "fixed_bond"
  character(len("fix_bondlength")),private,parameter :: tag_fix_bondlength = "fix_bondlength"
  character(len("fixed_normal_hypervector")),private,parameter :: &
       &                                            tag_fixed_normal_hypervector = "fixed_normal_hypervector"

! === KT_add === 2014/07/14
  character(len("hypervec_from_file")),private,parameter :: &
       &                      tag_hypervec_from_file = "hypervec_from_file"
  integer :: hypervec_from_file = OFF
! ============== 2014/07/14

!!! nudged_elastic_band_method
  character(len("accuracy")),private,parameter :: tag_accuracy = "accuracy"
  !!!character(len("constraint")),private,parameter :: tag_constraint = "constraint"
  !!!character(len("structure")),private,parameter :: tag_structure = "structure"

  integer :: sw_path_from_dynm = OFF
  character(len("sw_path_from_dynm")),private,parameter :: tag_sw_path_from_dynm="sw_path_from_dynm"
  integer :: neb_max_iteration = 10
  integer :: frame_end0 = -1
  integer :: frame_end1 = -1
  character(len("neb_max_iteration")),private,parameter :: tag_neb_max_iteration = "neb_max_iteration"
  real(kind=DP) :: neb_dt = 20.0d0
  character(len("dt")),private,parameter :: tag_neb_dt = "dt"
  integer :: ci_neb = OFF
  real(kind=DP) :: ci_thres = 1e+30
  real(kind=DP) :: to_dimer_thres = 0.d0
  integer :: ci_index=0
  character(len("end0_energy")),private,parameter :: tag_end0_energy = "end0_energy"
  character(len("end1_energy")),private,parameter :: tag_end1_energy = "end1_energy"
  real(kind=DP) :: end0_energy,end1_energy
  logical :: end0_energy_given = .false.,end1_energy_given = .false.
  character(len("neb_time_integral")),private,parameter :: tag_neb_time_integral = "neb_time_integral"
  integer :: neb_time_integral = -7
  character(len("dimer_time_integral")),private,parameter :: tag_dimer_time_integral = &
                "dimer_time_integral"
  integer :: dimer_time_integral = QUENCHED_MD
  real(kind=DP) :: sd_factor = 50000.d0
  character(len("sd_factor")),private,parameter :: tag_sd_factor = "sd_factor"
  character(len("ci_neb")),private,parameter :: tag_ci_neb = "ci_neb"
  character(len("ci_index")),private,parameter :: tag_ci_index = "ci_index"
  character(len("ci_thres")),private,parameter :: tag_ci_thres = "ci_thres"
!  real(kind=DP) :: sp_k_init = 1.0d0, sp_k_min = 1.0d0, sp_k_max = 1.0d0
  real(kind=DP) :: sp_k_init = 0.03d0, sp_k_min = 0.03d0, sp_k_max = 0.03d0
  character(len("sp_k_init")),private,parameter :: tag_sp_k_init = "sp_k_init"
  character(len("sp_k_min")),private,parameter :: tag_sp_k_min = "sp_k_min"
  character(len("sp_k_max")),private,parameter :: tag_sp_k_max = "sp_k_max"
  integer :: sp_k_variable = OFF
  character(len("sp_k_variable")),private,parameter :: tag_sp_k_variable = "sp_k_variable"
  integer :: penalty_function = OFF
  character(len("penalty_function")),private,parameter :: tag_penalty_function = "penalty_function"
  integer :: neb_convergence_condition = 1
  character(len("neb_convergence_condition")),private,parameter :: tag_neb_convergence_condition = "neb_convergence_condition"
  real(kind=DP) :: neb_convergence_threshold=1.d-5
  character(len("neb_convergence_threshold")),private,parameter :: tag_neb_convergence_threshold = "neb_convergence_threshold"
  character(len("to_dimer_threshold")),private,parameter :: tag_to_dimer_threshold = "to_dimer_threshold"
  character(len("dimer_convergence_threshold")),private,parameter :: &
          tag_dimer_convergence_threshold = "dimer_convergence_threshold"
  character(len("sw_random_unit_vector")), private, parameter :: &
          tag_sw_random_unit_vector="sw_random_unit_vector"
  character(len("coefficient")),private,parameter :: tag_coefficient = "coefficient"
  character(len("mode_coefficient")),private,parameter :: tag_mode_coefficient = "mode_coefficient"
  character(len("bondlength_fix")),private,parameter :: tag_bondlength_fix = "bondlength_fix"
  character(len("cog_fix")),private,parameter ::    tag_cog_fix    = "cog_fix"
  character(len("cog_fix_l")),private,parameter ::  tag_cog_fix_l  = "cog_fix_l"
  character(len("cog_fix_in_a_plane")),private,parameter :: tag_cog_fix_in_a_plane = "cog_fix_in_a_plane"
  character(len("rigid_body_fix_in_a_plane")),private,parameter:: tag_rigid_body_fix_in_a_plane = "rigid_body_fix_in_a_plane"
  character(len("rigid_body_fix_l")),private,parameter:: tag_rigid_body_fix_l = "rigid_body_fix_l"
  character(len("rigid_body")), private, parameter :: tag_rigid_body = "rigid_body"
  character(len("tag_sw_apply_constant_force")),private,parameter :: tag_sw_apply_constant_force = "sw_apply_constant_force"
  character(len("tag_sw_apply_constant_move")),private,parameter  :: tag_sw_apply_constant_move  = "sw_apply_constant_move"
  character(len("tag_sw_rigid_body_cog_movement")),private,parameter :: &
                 tag_sw_rigid_body_cog_movement = "sw_rigid_body_cog_movement"
  character(len("tag_sw_relax_rigid_body_cog")),private,parameter :: tag_sw_relax_rigid_body_cog = "sw_relax_rigid_body_cog"
  character(len("tag_force_apply")),private,parameter :: tag_force_apply = "force_apply"
  integer :: apply_constant_force = 0
  integer :: apply_constant_move = 0
  integer :: relax_rigid_body_cog = 0
  real(kind=DP) :: force_apply = 0.d0
  character(len("type")),private,parameter ::       tag_type       = "type"
  character(len("absolute")),private,parameter ::   tag_absolute   = "absolute"
  character(len("square")),private,parameter ::     tag_square     = "square"
  character(len("atom1")),private,parameter ::       tag_atom1      = "atom1"
  character(len("atom2")),private,parameter ::       tag_atom2      = "atom2"
  character(len("atom3")),private,parameter ::       tag_atom3      = "atom3"
  character(len("atom4")),private,parameter ::       tag_atom4      = "atom4"
  character(len("length")),private,parameter ::     tag_length     = "length"
  character(len("num_planes")),private,parameter :: tag_num_planes = "num_planes"
  character(len("fix_in_a_plane")),private,parameter :: tag_fix_in_a_plane = "fix_in_a_plane"
  character(len("fixed_plane")),private,parameter :: tag_fixed_plane = "fixed_plane"
  character(len("nx")),private,parameter ::         tag_nx         = "nx"
  character(len("ny")),private,parameter ::         tag_ny         = "ny"
  character(len("nz")),private,parameter ::         tag_nz         = "nz"
  character(len("delta")),private,parameter ::      tag_delta      = "delta"

  character(len("element_list")),private,parameter :: tag_element_list = "element_list"
  character(len("atomicnumber")),private,parameter :: tag_atomicnumber = "atomicnumber"
  character(len("mass")),private,parameter ::       tag_mass      = "mass"
  character(len("zeta")),private,parameter ::       tag_zeta      = "zeta"
!!$  character(len("variance")),private,parameter ::   tag_variance  = "variance"
  character(len("deviation")),private,parameter ::  tag_deviation = "deviation"
  character(len("standard_deviation")),private,parameter :: tag_standard_deviation = "standard_deviation"
  character(len("dev")),private,parameter ::        tag_dev       = "dev"
  character(len("qex")),private,parameter ::        tag_qex       = "qex"

  character(len('dimer_method')), private, parameter :: tag_dimer_method = 'dimer_method'
  character(len('dimer_max_iteration')), private, parameter :: tag_dimer_max_iteration = 'dimer_max_iteration'
  character(len('delta_theta')), private, parameter :: tag_delta_theta = 'delta_theta'
  character(len('delta_r')), private, parameter :: tag_delta_r = 'delta_r'
  integer :: sw_dimer_method = OFF
  integer :: dimer_max_iteration = 2000
  real(kind=DP) :: delta_theta = 1.d-3
  real(kind=DP) :: delta_r = 0.01d0
  integer :: sw_random_unit_vector = OFF

  ! --- PHONON FORCE ---
  integer :: num_force_data
  integer :: num_force_calc
  integer :: num_force_calc_mobile
  integer :: istart_phonon = -1
  integer :: iend_phonon   = -1
  integer, allocatable, dimension(:) :: phonon_atom ! dim(num_force_data)
  real(kind=DP), allocatable, dimension(:,:) :: phonon_displacement ! dim(num_force_data,3)
  real(kind=DP) :: u = 0.d0
  integer,       allocatable, dimension(:,:) :: napt_phonon !d(natm_super,num_force_data)
  integer,       allocatable, dimension(:)   :: iequconf !d(num_force_data)
  integer,       allocatable, dimension(:)   :: iopr_equconf !d(num_force_data)
  integer,       allocatable, dimension(:)   :: iconf !d(num_force_calc)

  ! --- Vibrational mode ---
  integer            :: mode_index = 1
  real(kind=DP)      :: normal_coordinate = 0.d0
  real(kind=DP), allocatable, dimension(:,:) :: xi_mode ! dim(natm,3)

  ! --- NEB or Multiple replica ---
  integer :: number_of_replicas = 1
  integer, allocatable, dimension(:) :: replica_howtogive_coordinates ! d(number_of_replica)
  integer, allocatable, dimension(:,:) :: replica_endpoints           ! d(2,number_of_replica)
  character(len("multiple_replica")),private,parameter :: tag_multiple_replica = "multiple_replica"
  character(len("number_of_replicas")),private,parameter :: tag_number_of_replicas = "number_of_replicas"
  character(len("replicas")),private,parameter :: tag_replicas = "replicas"
  character(len("replica_numbers")),private,parameter :: tag_replica_numbers = "replica_numbers"
  character(len("proportional")),private,parameter :: tag_proportional = "proportional"
  character(len("end0")),private,parameter ::         tag_end0 = "end0"
  character(len("end1")),private,parameter ::         tag_end1 = "end1"

  ! --- approximate DFT+U : Hubbard model ---
  integer, allocatable :: ihubbard(:)  ! d(natm)
  character(len("hubbard")),private,parameter :: tag_hubbard = "hubbard"

  ! --- Projector group
  integer, allocatable :: iproj_group(:) ! dim(natm)
  character(len("proj_group")),private,parameter :: tag_proj_group = "proj_group"

! ================================== added by K. Tagami =============== 11.0
!
!!  -- NonCollinear --
!
  character(len("mdx")),private,parameter :: tag_mdx = "mdx"
  character(len("mdy")),private,parameter :: tag_mdy = "mdy"
  character(len("mdz")),private,parameter :: tag_mdz = "mdz"
  character(len("theta")),private,parameter :: tag_theta = "theta"
  character(len("phi")),private,parameter :: tag_phi = "phi"

  character(len("moment")),private,parameter :: tag_moment = "moment"
  character(len("mx")),private,parameter :: tag_mx = "mx"
  character(len("my")),private,parameter :: tag_my = "my"
  character(len("mz")),private,parameter :: tag_mz = "mz"
!
  real(kind=DP), allocatable :: mag_direction0_atomtyp(:,:)
  real(kind=DP), allocatable :: magmom_local_now(:,:)
!
! --
  character(len("lcore_parfil")),private,parameter :: &
       &                         tag_lcore_parfil = "lcore_parfil"
  integer, allocatable :: has_partially_filled_lcore(:)
! ===================================================================== 11.0

! ================ KT_add ===================== 13.0U
!
!! -- Initial Charge --
!
  character(len("ion_charge")),private,parameter :: tag_ion_charge = "ion_charge"
!
  logical :: mag_zeta1_atomtyp_is_defined = .false.
  logical :: mag_moment0_atomtyp_is_defined = .false.
  logical :: mag_moment0_atoms_is_defined = .false.
!
  real(kind=DP), allocatable :: mag_moment0_atomtyp(:,:)
  real(kind=DP), allocatable :: mag_moment0_atoms(:,:)
  real(kind=DP), allocatable :: mag_direction0_atoms(:,:)
  real(kind=DP), allocatable :: ionic_charge_atoms(:)
  real(kind=DP), allocatable :: ionic_charge_atomtyp(:)

  character(len("sw_set_initial_chgmag_by_atom")),private,parameter :: &
       &  tag_sw_initial_chgmag_by_atom = "sw_set_initial_chgmag_by_atom"
  character(len("sw_set_initial_magmom_by_atom")),private,parameter :: &
       &  tag_sw_initial_magmom_by_atom = "sw_set_initial_magmom_by_atom"
  integer :: sw_set_initial_magmom_by_atom = OFF

! ============================================= 13.0U

! ================================ added by K. Tagami ============= 11.0
!
!! --- Spin-Orbit  ---
!
  integer, allocatable :: itab_spinorbit_addition(:)  ! d(natm)
!
  character(len("scaling_so")),private,parameter :: tag_scaling_so = "scaling_so"
  real(kind=DP), allocatable :: scaling_so(:)
! ================================================================= 11.0

  ! --- van der Waals
  integer :: ntyp_vdw
  integer, allocatable :: ityp_vdw(:) !d(natm)
  real(kind=DP) :: evdw ! vdW energy
  real(kind=DP), allocatable, dimension(:,:) :: fxyzvdw_l ! d(natm,3)
  real(kind=DP), allocatable :: cvdw(:,:) !d(ntyp_vdw,ntyp_vdw)
  real(kind=DP), allocatable :: rvdw(:,:) !d(ntyp_vdw,ntyp_vdw)
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
  real(kind=DP), allocatable :: c6vdw(:) !d(ntyp_vdw)
  real(kind=DP), allocatable :: r0vdw(:) !d(ntyp_vdw)
  real(kind=DP), allocatable :: pvdw(:) !d(ntyp_vdw)
! ==============================================================================
  character(len=LEN_ATOMNAME),allocatable, dimension(:) ::  speciesname_vdw ! d(ntyp_vdw)
  character(len=LEN_ATOMNAME),private,allocatable,dimension(:) :: species_vdw_work ! d(natm)
  character(len=LEN_ATOMNAME),private,allocatable,dimension(:) :: species_vdw_indp ! d(natm)
  character(len("vdw")),private,parameter :: tag_vdw = "vdw"
  character(len("vdw_list")),private,parameter :: tag_vdw_list = "vdw_list"
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
! character(len("type1")),private,parameter :: tag_type1 = "type1"
! character(len("type2")),private,parameter :: tag_type2 = "type2"
  !!character(len("type")),private,parameter :: tag_type = "type"
! ==============================================================================
  character(len("c6")),private,parameter :: tag_c6 = "c6"
  character(len("r0")),private,parameter :: tag_r0 = "r0"
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
  character(len("p")),private,parameter :: tag_p = "p"
! ==============================================================================

! ======================================== KT_add ================ 13.0B
  real(kind=DP) :: s_vdw(3,3)
! ================================================================ 13.0B

#ifdef __EDA__
  real(kind=DP), allocatable :: evdw_on_atom(:)
#endif

  real(kind=DP), allocatable, dimension(:) :: ival

  ! --- for storing/retrieving atomic coordinates
  type atomic_configuration_t
    integer :: id
    integer :: group
    character(len=LEN_ATOMNAME) :: element
    character(len=LEN_ATOMNAME) :: element_vdw
    real(kind=DP),dimension(3) :: pos,cps,pos_in,cps_in,cpd_l
    real(kind=DP),dimension(3,3) :: cpo_l
    integer :: iwei,imdtyp,ityp,if_pdos,if_aldos,ihubbard,iproj_group,numlay,exclusion_target
    integer,dimension(3) :: imdtypxyz
    real(kind=DP) :: ionic_mass
    real(kind=DP) :: nvalence
  end type atomic_configuration_t
  type(atomic_configuration_t), allocatable, dimension(:) :: config_buf
  integer :: nconfig_buf

  ! --- for dynamical addition/removal of atoms
  type(atomic_configuration_t), allocatable, dimension(:) :: atom_reservoir
  integer :: natom_reservoir = -1
  integer :: curr_atom_reservoir = 1
  integer :: natom_group = -1
  integer,allocatable,dimension(:) :: natm_per_group
  integer,allocatable,dimension(:,:) :: atomid_in_group
  character(len("reservoir")),private,parameter :: tag_reservoir = "reservoir"
  character(len("frequency")),private,parameter :: tag_frequency = "frequency"
  character(len("sw_rotate_reservoir")),private,parameter :: tag_sw_rotate_reservoir = "sw_rotate_reservoir"
  character(len("atom_addition_criteria")),private,parameter :: &
 & tag_atom_addition_criteria = "atom_addition_criteria"
  character(len("curr_atom_reservoir")), private, parameter :: tag_curr_atom_reservoir = "curr_atom_reservoir"
  character(len("reservoir_group")),private,parameter :: tag_reservoir_group = "reservoir_group"

  character(len("atom_exclusion_criteria")),private,parameter :: &
 & tag_atom_exclusion_criteria = "atom_exclusion_criteria"
  character(len("exclusion_target")),private,parameter :: &
 & tag_exclusion_target="exclusion_target"
  character(len("sw_atom_excludable")),private,parameter :: &
 & tag_sw_atom_excludable="sw_atom_excludable"
  character(len("x_greater_than")),private,parameter :: tag_x_greater_than = "x_greater_than"
  character(len("y_greater_than")),private,parameter :: tag_y_greater_than = "y_greater_than"
  character(len("z_greater_than")),private,parameter :: tag_z_greater_than = "z_greater_than"
  character(len("x_less_than")),private,parameter :: tag_x_less_than = "x_less_than"
  character(len("y_less_than")),private,parameter :: tag_y_less_than = "y_less_than"
  character(len("z_less_than")),private,parameter :: tag_z_less_than = "z_less_than"
  integer :: sw_atom_excludable = OFF
  integer, allocatable, dimension(:) :: exclusion_target
  real(kind=DP), dimension(3) :: exclusion_criteria_min
  real(kind=DP), dimension(3) :: exclusion_criteria_max

  integer :: addition_frequency = -1
  integer :: sw_rotate_reservoir = OFF
  real(kind=DP) :: neg_incre

! === check of atomic coordinates ===
  integer :: sw_check_coordinates = ON
  real(kind=DP) :: bond_length_min = 0.1
  character(len("sw_check_coordinates")),private,parameter :: tag_sw_check_coordinates = "sw_check_coordinates"
  character(len("bond_length_min")),private,parameter :: tag_bond_length_min = "bond_length_min"

  ! for the new CG optimizer
  integer, save :: iter_CG = 1, iter_linmin = 1, iter_CG_max = 1

  real(kind=DP),dimension(3) :: latconst_len,latconst_angle

  real(kind=DP) :: tdamp
  real(kind=DP) :: tdamp_baro

  ! parameters for NPT
  real(kind=DP) :: m_baro = 1.d0
  real(kind=DP) :: target_pressure = 0.d0
  character(len("temperature_pressure_control")), private, parameter :: tag_temperature_pressure_control = &
  & "temperature_pressure_control"
  character(len("temperature_pressure_control")), private, parameter :: tag_pressure_temperature_control = &
  & "pressure_temperature_control"
  character(len("pressure_control")), private, parameter :: tag_pressure_control = &
  & "pressure_control"
  character(len("mass_baro")), private, parameter :: tag_m_baro  = "mass_baro"
  character(len("pressure")), private, parameter :: tag_pressure = "pressure"
  character(len("press")), private, parameter :: tag_press   = "press"
  character(len("pressi")), private, parameter :: tag_pressi = "pressi"
  character(len("pressf")), private, parameter :: tag_pressf = "pressf"
  character(len("m11")), private, parameter :: tag_m11 = 'm11'
  character(len("m22")), private, parameter :: tag_m22 = 'm22'
  character(len("m33")), private, parameter :: tag_m33 = 'm33'
  character(len("m12")), private, parameter :: tag_m12 = 'm12'
  character(len("m13")), private, parameter :: tag_m13 = 'm13'
  character(len("m23")), private, parameter :: tag_m23 = 'm23'
  character(len("control_pressure")), private, parameter :: tag_control_pressure = "control_pressure"
  integer, dimension(3,3) :: mobility_cell
  integer :: control_pressure = ON ! for debugging purposes
  character(len("volume")),private,parameter :: tag_volume="volume"
  character(len("metric_tensor")),private,parameter :: tag_metric_tensor="metric_tensor"
  integer :: p_ctrl_method = METRIC_TENSOR

  character(len("sw_average_diagonal")), private, parameter :: tag_sw_average_diagonal = &
  & "sw_average_diagonal"
  integer :: sw_average_diagonal = OFF

  character(len("barostat")),private,parameter :: tag_barostat  = "barostat"

! paw
  integer, allocatable, dimension(:) :: surface_integral_method_paw
  character(len("surface_integral_method")), private, parameter :: tag_surface_integral_method = &
            & "surface_integral_method"
  character(len("sph")), private, parameter :: tag_sph = "sph"
  character(len("gl")), private, parameter :: tag_gl = "gl"

  character(len("sw_temperature_profile")), private, parameter ::tag_sw_temperature_profile = "sw_temperature_profile"
  integer :: sw_temperature_profile = OFF

  character(len("sw_pressure_profile")), private, parameter ::tag_sw_pressure_profile = "sw_pressure_profile"
  integer :: sw_pressure_profile = OFF

! tags for temperature profile
  character(len("tempi")), private, parameter :: tag_tempi = "tempi"
  character(len("tempf")), private, parameter :: tag_tempf = "tempf"
  character(len("till_n")), private, parameter :: tag_till_n = "till_n"

! type for specifying the temperature profile
  type temperature_profile_t
    integer :: no
    real(kind=DP), pointer, dimension(:) :: tempi,tempf
    real(kind=DP), pointer, dimension(:) :: qmass,tdamp
    integer, pointer, dimension(:) :: till_n
    integer :: nprof, currprof, natm
  end type temperature_profile_t

! type for specifying the pressure profile
  type pressure_profile_t
    integer :: no
    real(kind=DP), pointer, dimension(:) :: pressi,pressf
    real(kind=DP), pointer, dimension(:) :: bmass,tdamp_baro
    integer, pointer, dimension(:) :: till_n
    integer :: nprof, currprof
  end type pressure_profile_t

  type(temperature_profile_t), allocatable, dimension(:) :: temperature_profile
  type(pressure_profile_t) :: pressure_profile

! tags for defining regions
  character(len("region")),      private, parameter :: tag_region="region"
  character(len("cylinder")),    private, parameter :: tag_cylinder="cylinder"
  character(len("box")),         private, parameter :: tag_box="box"
  character(len("xmin")),        private, parameter :: tag_xmin="xmin"
  character(len("xmax")),        private, parameter :: tag_xmax="xmax"
  character(len("ymin")),        private, parameter :: tag_ymin="ymin"
  character(len("ymax")),        private, parameter :: tag_ymax="ymax"
  character(len("zmin")),        private, parameter :: tag_zmin="zmin"
  character(len("zmax")),        private, parameter :: tag_zmax="zmax"
  character(len("radius")),      private, parameter :: tag_radius="radius"
  character(len("orientation")), private, parameter :: tag_orientation="orientation"
  character(len("cylx")),        private, parameter :: tag_cylx="cylx"
  character(len("cyly")),        private, parameter :: tag_cyly="cyly"
  character(len("cylzmin")),     private, parameter :: tag_cylzmin="cylzmin"
  character(len("cylzmax")),     private, parameter :: tag_cylzmax="cylzmax"
  character(len("sigma")),       private, parameter :: tag_sigma="sigma"
  character(len("epsilon")),     private, parameter :: tag_epsilon="epsilon"
  character(len("sw_tally")),       private, parameter :: tag_sw_tally="sw_tally"
  character(len("region_group")),private, parameter :: tag_region_group = "region_group"

  type region_t
    integer :: no
    integer :: region_group
    integer :: region_type
    real(kind=DP) :: xmin,ymin,zmin
    real(kind=DP) :: xmax,ymax,zmax
    real(kind=DP) :: radius
    integer :: orientation
    real(kind=DP) :: cylx,cyly
    real(kind=DP) :: cylzmin,cylzmax
    real(kind=DP) :: sigma,epsi
    real(kind=DP) :: energy
    real(kind=DP), pointer, dimension(:,:) :: forc
    integer :: i1,i2,i3
    logical :: tally
    integer :: ntarget_atoms
    integer, pointer, dimension(:) :: target_atoms
  end type region_t
  type(region_t), allocatable, dimension(:) :: regions
  integer :: num_regions = 0

  character(len("dftd3")), private, parameter :: tag_dftd3="dftd3"
  character(len("k1")), private, parameter :: tag_k1 = "k1"
  character(len("k2")), private, parameter :: tag_k2 = "k2"
  character(len("k3")), private, parameter :: tag_k3 = "k3"
  character(len("s6")), private, parameter :: tag_s6 = "s6"
  character(len("s8")), private, parameter :: tag_s8 = "s8"
  character(len("sr6")), private, parameter :: tag_sr6 = "sr6"
  character(len("sr8")), private, parameter :: tag_sr8 = "sr8"
  character(len("alpha6")), private, parameter :: tag_alpha6 = "alpha6"
  character(len("alpha8")), private, parameter :: tag_alpha8 = "alpha8"
  character(len("rcut")), private, parameter :: tag_rcut = "rcut"
  character(len("rcut_nc")), private, parameter :: tag_rcut_nc = "rcut_nc"
  character(len("a1")), private, parameter :: tag_a1 = "a1"
  character(len("a2")), private, parameter :: tag_a2 = "a2"

  type dftd3par_t
    integer :: maxelem,maxmaxnc,nlines
    integer, pointer, dimension(:) :: maxnc
    real(kind=DP), pointer, dimension(:,:,:,:) :: c6ab,nc1,nc2
    real(kind=DP), pointer, dimension(:,:) :: r0ab
    real(kind=DP), pointer, dimension(:) :: r2r4,covrad
    real(kind=DP) :: k1,k2,k3,s6,sr6,s8,sr8,sr9,alpha6,alpha8,alpha9
    real(kind=DP) :: a1, a2
    real(kind=DP) :: rcut_nc,rcut_vdw
    logical :: read_pars = .false.
  end type dftd3par_t

  type(dftd3par_t) :: dftd3par

  character(len("phase0_input")), private, parameter :: tag_phase0_input = "phase0_input"
  character(len("phase0_output")), private, parameter :: tag_phase0_output = "phase0_output"
  character(len("frame")), private, parameter :: tag_frame = "frame"
  character(len("frame_end0")), private, parameter :: tag_frame_end0 = "frame_end0"
  character(len("frame_end1")), private, parameter :: tag_frame_end1 = "frame_end1"
  character(len("filetype")), private, parameter :: tag_filetype = "filetype"
  integer :: filetype = PHASE0_OUTPUT
  integer :: coord_method = DIRECTIN
  integer :: frame = -1

! partial force data
  integer :: phonon_iteration=1
  logical, allocatable, dimension(:) :: force_was_read

  real(kind=DP), allocatable, dimension(:,:,:) :: xk,fk

  integer :: iseed = -1

  integer, private, allocatable, dimension(:) :: iproj_group_tmp

! rigid body dynamics
  type rigid_body_t
    integer                                :: id, natm
    integer,       pointer, dimension(:)   :: atoms
    real(kind=DP), pointer, dimension(:,:) :: relative_coords,molecular_coords
    real(kind=DP), pointer, dimension(:,:) :: force_per_atm
    real(kind=DP),          dimension(3)   :: COM
    real(kind=DP)                          :: mass
    real(kind=DP),          dimension(4)   :: q, q_h
    real(kind=DP),          dimension(4,4) :: Qmat
    real(kind=DP),          dimension(3)   :: force
    real(kind=DP),          dimension(3)   :: velocity, velocity_old, velocity_h
    real(kind=DP),          dimension(3)   :: angular_momentum, angular_momentum_h
    real(kind=DP),          dimension(3)   :: angular_momentum_old
    real(kind=DP),          dimension(3)   :: torque
    real(kind=DP),          dimension(3,3) :: rotmat, rotmat_old
    real(kind=DP),          dimension(3)   :: inertia
    real(kind=DP)                          :: dt_translation, dt_rotation
    integer, dimension(3)                  :: mobile
    integer                                :: mobilerot
    integer                                :: thermo_group
    logical                                :: explicit_zaxis
    real(kind=DP),          dimension(3)   :: zaxis, xaxis, yaxis
    real(kind=DP)                          :: phi,theta,psi,phi_h,theta_h,psi_h
    integer                                :: fix_in_a_plane

!   plane on which the rigid-body is constrained
    logical                                :: constrain_on_plane, move_plane
    real(kind=DP)                          :: normx,normy,normz
    real(kind=DP)                          :: incx,incy,incz,dinc,incmax
    real(kind=DP)                          :: distance_moved_thusfar
  end type rigid_body_t

  integer :: nrigid_bodies=0
  real(kind=DP) :: torque_max, trans_force_max
  integer :: sw_full_rb_dynamics = OFF
  type(rigid_body_t), allocatable, target, dimension(:) :: rigid_bodies
  logical, allocatable, dimension(:)               :: is_rigid_body
!!$  character(len("rigid_body")), private, parameter :: tag_rigid_body = "rigid_body"
  character(len("dt_translation")), private, parameter :: tag_dt_translation = "dt_translation"
  character(len("dt_rotation")), private, parameter :: tag_dt_rotation = "dt_rotation"
  character(len("mobilerot")), private, parameter :: tag_mobilerot="mobilerot"
  character(len("sw_full_rb_dynamics")), private, parameter :: tag_sw_full_rb_dynamics="sw_full_rb_dynamics"
!!$  character(len("normx")), private, parameter :: tag_normx = "normx"
!!$  character(len("normy")), private, parameter :: tag_normy = "normy"
!!$  character(len("normz")), private, parameter :: tag_normz = "normz"

! fix bonds between atoms of specified elements (typically hydrogen)
  type fix_bond_t
    integer, dimension(2)         :: atoms
    integer, dimension(2)         :: ityp
    real(kind=DP)                 :: bond_length
    real(kind=DP)                 :: lambda, sigma
    real(kind=DP), dimension(2,3) :: dsigma, dsigma_old
  end type fix_bond_t

  real(kind=DP), allocatable, dimension(:,:) :: bond_forc
  integer :: fix_bond_nelements
  integer, allocatable, dimension(:) :: fix_bond_elements
  integer :: sw_fix_bond = OFF
  integer :: nfixed_bonds = 0
  real(kind=DP) :: thres_fix_bond = 1.d-10
  real(kind=DP) :: eps_shake  = 1.d-14
  real(kind=DP) :: eps_rattle = 1.d-14
  type(fix_bond_t), allocatable, target, dimension(:) :: fixed_bond
  real(kind=DP), allocatable, dimension(:)   :: covalent_radii
  real(kind=DP) :: bond_factor = 1.2d0
  integer :: max_iter_fix_bond = 1000
  integer, allocatable, dimension(:) :: nbonds_per_thermo
  character(len("fix_bond")), private, parameter :: tag_fix_bond = "fix_bond"
  character(len("sw_fix_bond")), private, parameter :: tag_sw_fix_bond = "sw_fix_bond"
  character(len("max_iter_fix_bond")), private, parameter :: tag_max_iter_fix_bond = "max_iter_fix_bond"
  character(len("thres_fix_bond")), private, parameter :: tag_thres_fix_bond = "thres_fix_bond"
  character(len("covalent_radius")), private, parameter :: tag_covalent_radius="covalent_radius"
  character(len("bond_factor")), private, parameter :: tag_bond_factor="bond_factor"
  character(len("target_element")) :: tag_target_element="target_element"

! subroutines contained here
!     m_IS_set_iatom
!     m_IS_rd_n
!        -- specify_ityp, wd_atom_list, count_species, set_input_coordinate_system,
!           set_atompos_and_etc, set_element_detail
!     m_IS_alloc_iatomn_etc
!     m_IS_alloc_napt
!     m_IS_alloc_fxyzew
!     m_IS_alloc_zfm3
!     m_IS_gdiis_alloc
!   1.m_IS_rd_pos_and_v
!   2.m_IS_wd_pos_and_v
!     - copy_cpd_l_to_pwork, - copy_cpo_l_to_pwork
!   3.m_IS_alloc_pos_and_v          <-(InputData_Analysis)
!   4.m_IS_cp_cps2cpo
!   5.m_IS_md                            ->(8,7,12,11,10,6,9)
!   7.  check_if_bondlength_fix_exist <-(5)
!
!   8.  md1_alloc
!   9.  md1_dealloc
!  10.  quench_velocities
!  11.  get_ekina
!  12.  evolve_velocities
!  13.m_IS_structure_factor
!     - structure_factor1,   - structure_factor2,  - wd_zfm3
!  14.m_IS_ewald
!     - wd_eewald_and_fxyzew, - ewald_Rspace_summation, - ewald_Gspace_summation
!     - get_zsum, - add_exp_G2_zsum, - ewald_force_Gspace_summation,
!     - set_ewald_parameters, - cpspac, - decide_newldg,  - decide_alf
!  15.m_IS_symm_check_of_pos
!     - symm_check_of_ions_positions_c
!  16.  decide_rxyz_size
!  17.  substitute_rxyz
!  18.m_IS_initialize_mdmode   <-(Initialization)
!  19.m_IS_initialize_cpd_l
!  20.m_IS_cps_to_pos
!  21.m_IS_wd_forc
!
!  --> temparature control
!  22.m_IS_rd_T_parameters
!  23.  T_control_alloc
!  23b. forcp_alloc
!  24.m_IS_rd_forcp_etc
!  25.m_IS_wd_forcp_etc
!  26.  check_imdtyp
!  27.  vlcty_accrd2_vVerlet
!  28.  ekina_ekinq_ekbt_and_ega
!  29.  evolve_crdn_ACCRD2_vVerlet
!  30.  evolve_cprv
!  31.  md2_alloc
!  32.  md2_dealloc
!  33.  heatrsv
!  34.m_IS_md_thermo
!     - check_nrsv
!  35.  rattle_v
!  36.  rattle_r
!     - stop0, evolve_almda, evolve_cps
!  38.m_IS_md_bluem
!     - print_frsv_and_cpqr, - init_md_bluem
!  39.m_IS_md_cnstr
!     - evaluate_forcmx, - quench_velocity_using_ifq, - init_md_cnstr
!  40.m_IS_force_check_md_cnstr
!  41.m_IS_alloc_cnstrvectors_etc
!  42.m_IS_cp_works2fcvect_etc
!  43.m_IS_init_cnstrnt
!  44.m_IS_wd_cpo_and_forc
!     m_IS_gdiis
!       forc_check
!       cps_check
!       cpd_check
!
contains
  subroutine m_IS_put_latconst_len_angle(length,angle)
    real(kind=DP), intent(in), dimension(3) :: length, angle
    latconst_len   = length
    latconst_angle = angle
  end subroutine m_IS_put_latconst_len_angle

  subroutine m_IS_put_lattice_system(lattice_system)
    character(len=9), intent(in) :: lattice_system
    lattice_system_from_m_CS_SG = lattice_system
  end subroutine m_IS_put_lattice_system

  subroutine alloc_normal_hypervector()
    if ( .not. allocated(normal_hypervector) ) then
       allocate(normal_hypervector(natm,3,PUCV:CARTS))
    endif
    normal_hypervector = 0.0d0
  end subroutine alloc_normal_hypervector

  subroutine set_normal_hypervector()
    integer :: i, j
    real(kind=DP) :: r
    do j = 1, 3
       do i = 1, natm
          if(imdtyp(i) == FIX) then
             normal_hypervector(i,j,PUCV) = 0.d0
          else
             r = pos_end1(i,j) - pos_end0(i,j)
             if(r > 0.5d0) then
                normal_hypervector(i,j,PUCV) = r - floor(r+0.5d0)
             else if( r <= -0.5d0) then
                normal_hypervector(i,j,PUCV) = r + floor(-r+0.5d0)
             else
                normal_hypervector(i,j,PUCV) = r
             end if
          end if
       end do
    end do
    call change_of_coordinate_system(altv,normal_hypervector(1,1,PUCV),natm,natm,normal_hypervector(1,1,CARTS))
  end subroutine set_normal_hypervector

! ===== KT_add === 2014/07/14
  subroutine read_hypervec_from_file
    integer :: i, j

    normal_hypervector = 0.0d0
    if ( mype == 0 ) then
       Do i=1, natm
          read(nfhypervec,*) ( normal_hypervector(i,j,CARTS),j=1,3 )
       End do
    endif

    if ( npes > 1 ) then
       call mpi_bcast( normal_hypervector, natm*3*2, &
            &          mpi_double_precision, 0, MPI_CommGroup, ierr )
    endif
  end subroutine read_hypervec_from_file
! ================ 2014/07/14

  subroutine wd_normal_hypervector()
    integer :: i, j
    write(nfout,'(" !** == hypervector ==")')
    write(nfout,'(" !**     no.",10x,"(internal)",25x,"(cartesian)")')
    do i = 1, natm
       write(nfout,'(" !** ",i5," : ", 3f12.8," : ",3f12.8)') i &
            & , (normal_hypervector(i,j,PUCV),j=1,3), (normal_hypervector(i,j,CARTS),j=1,3)
    end do
  end subroutine wd_normal_hypervector

  subroutine alloc_bondlength_fix_set()
    allocate(bondlength_fix_set(2,num_fixed_bonds))
  end subroutine alloc_bondlength_fix_set

  subroutine m_IS_set_iatom(nfout)
    integer, intent(in) :: nfout
    integer :: i, k, mm
    do k = 1, ntyp
       iatom(k) = 0.d0
       do i = 1, natm
          if(ityp(i) /= k) cycle
          iatom(k) = iatom(k) + iwei(i)
       end do
    end do
! -- check of iwei
    mm = nint(sum(iatom(1:ntyp)))
    if(mm /= natm2 .and. printable) write(nfout,340) mm, natm2
340 format(' ',' sum of iwei .ne. natm2 mm,natm2=',2i6)
  end subroutine m_IS_set_iatom

  subroutine m_IS_set_iproj_group()
    integer :: ia
    integer, allocatable :: iproj_tmp(:)
    if(.not.allocated(iproj_group)) return
    if(allocated(iproj_group_tmp)) then
      do ia=1, natm
        if(iproj_group_tmp(ityp(ia)) > 0) then
          iproj_group(ia) = iproj_group_tmp(ityp(ia))
        endif
      enddo
      deallocate(iproj_group_tmp)
    endif
    allocate(iproj_tmp(ntyp));iproj_tmp=0
    do ia=1, natm
      if(iproj_tmp(ityp(ia))==0 .and. iproj_group(ia)>0) then
        iproj_tmp(ityp(ia)) = iproj_group(ia)
        cycle
      endif
      if(iproj_group(ia)>0 .and. iproj_tmp(ityp(ia)) /= iproj_group(ia)) then
        deallocate(iproj_tmp)
        call phase_error_with_msg(nfout,'!** only one proj_group can be assigned to an element',__LINE__,__FILE__)
      endif
    enddo
    deallocate(iproj_tmp)
  end subroutine m_IS_set_iproj_group

  subroutine m_IS_rd_n_pre(nfout)
    integer, intent(in) :: nfout
    character(len=FMAXVALLEN) :: rstr
    integer :: iret
    real(kind=DP) :: dret
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
    logical :: exi
    iret = f_selectTop()

    if( f_selectBlock( tag_structure) == 0) then
       if( f_selectBlock(tag_atomlist) == 0) then
         iret = f_selectParentBlock()
       else
         coord_method = FILE
         filetype = PHASE0_INPUT
       endif
       if( f_getStringValue(tag_method, rstr, LOWER) == 0) then
         if(trim(rstr) == tag_directin) coord_method = DIRECTIN
         if(trim(rstr) == tag_file) coord_method = FILE
       endif
       if(coord_method == FILE)then
         if(f_selectBlock(tag_file) == 0)then
           if(f_getStringValue(tag_filetype, rstr, LOWER) == 0)then
             if(rstr == tag_phase0_input) filetype = PHASE0_INPUT
             if(rstr == tag_phase0_output) filetype = PHASE0_OUTPUT
           endif
           if(f_getIntValue(tag_frame,iret) == 0)then
             frame = iret
           endif
           iret = f_selectParentBlock()
         endif
       endif
       iret = f_selectParentBlock()
    else
       coord_method = FILE
       filetype = PHASE0_INPUT
    endif
    if(printable)then
      if(coord_method == DIRECTIN) then
        write(nfout,'(a)') ' !** coordinate specification method : DIRECTIN'
      else if (coord_method == FILE) then
        write(nfout,'(a)') ' !** coordinate specification method : FILE'
        if(filetype == PHASE0_INPUT) then
        write(nfout,'(a)') ' !** filetype : PHASE0_INPUT'
        else if (filetype == PHASE0_OUTPUT) then
        write(nfout,'(a)') ' !** filetype : PHASE0_OUTPUT'
        write(nfout,'(a,i8)') ' !** frame : ',frame
        endif
      endif
    endif
    if(coord_method == FILE)then
      if(filetype == PHASE0_INPUT) then
        call m_Files_set_def_fname_pos(F_INP)
      else if (filetype == PHASE0_OUTPUT) then
        call m_Files_set_def_fname_pos(F_DYNM)
        call m_Files_set_def_fname_attr()
      endif
      call m_Files_rd_file_names_data()
      inquire(file=F_POS,exist = exi)
      if(.not.exi) then
        write(nfout,'(a)') 'F_POS : '//trim(F_POS)//' does not exist'
        flush(nfout)
        call phase_error_with_msg(nfout,'F_POS : '//trim(F_POS)//' does not exist',__LINE__,__FILE__)
      endif
      if(filetype == PHASE0_OUTPUT)then
        inquire(file=F_COORD_ATTR,exist = exi)
        if(.not.exi)then
          write(nfout,'(a)') 'F_COORD_ATTR : '//trim(F_COORD_ATTR)//' does not exist'
          flush(nfout)
          call phase_error_with_msg(nfout,'F_COORD_ATTR : '//trim(F_COORD_ATTR)//' does not exist',__LINE__,__FILE__)
        endif
      endif
    endif

  end subroutine m_IS_rd_n_pre

  subroutine m_IS_reread_imdtyp(nfout)
    integer, intent(in) :: nfout
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue &
    &        , f_selectParentBlock, f_selectFirstTableLine, f_selectNextTableLine
    integer :: iret, i, ip
    real(kind=DP) :: dret
    logical :: done_something
    iret = f_selectTop()
    if(f_selectBlock(tag_structure) == 0) then
      if(f_selectBlock(tag_atomlist) == 0) then
        if( f_selectBlock(tag_atoms) == 0) then
          i = 1
          do while(.true.)
            if (i == 1) then
              if(f_selectFirstTableLine() /= 0) then
                 exit
              end if
            else
              if(f_selectNextTableLine() /= 0) then
                 exit
              end if
            end if
            if(i > natm) exit
            ip = i
            call m_CtrlP_rd_val(nfout, tag_mobile, imdtyp(ip), .true., done_something)
            if( done_something ) then
              imdtypxyz(ip,1) = imdtyp(ip)
              imdtypxyz(ip,2) = imdtyp(ip)
              imdtypxyz(ip,3) = imdtyp(ip)
            endif
            call m_CtrlP_rd_val(nfout, tag_mobilex, imdtypxyz(ip,1), .true.)
            call m_CtrlP_rd_val(nfout, tag_mobiley, imdtypxyz(ip,2), .true.)
            call m_CtrlP_rd_val(nfout, tag_mobilez, imdtypxyz(ip,3), .true.)
            i = i + 1
          enddo
          iret = f_selectParentBlock()
        endif
        iret = f_selectParentBlock()
      endif
      iret = f_selectParentBlock()
    endif
  end subroutine m_IS_reread_imdtyp

  subroutine m_IS_rd_n(nfout)
    ! <m_IS_rd_n> reads atomic coorindates and information of species from an input file
    !     formatted in a new style.
    ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), Jun. 2003

    integer, intent(in) :: nfout
    character(len=FMAXVALLEN) :: rstr
    integer :: iret, icounted, i
    real(kind=DP) :: dret
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
    logical :: tf, number_is_given, prealloc
    real(kind=DP), allocatable, dimension(:,:) :: work
    real(kind=DP), allocatable, dimension(:,:) :: rltv_t
    real(kind=DP), allocatable, dimension(:,:) :: fcvect_work
    real(kind=DP), allocatable, dimension(:,:) :: cpd_l_t
    integer, allocatable, dimension(:) :: ityp_t
    character(len=LEN_ATOMNAME), allocatable, dimension(:) :: speciesname_t
    real(kind=DP), dimension(3,3) :: altv_t,rltv_tt
    integer, allocatable, dimension(:) ::   imdtyp_set ! d(natm)
    integer, allocatable, dimension(:,:) :: iwork
    integer :: fixed_atoms, input_coordinate_system_t, istat
    integer :: ig
    integer :: seedsize
    integer, allocatable, dimension(:) :: seeds
    logical, save :: read_from_file=.false.

    iret = f_selectTop()

    if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** << m_IS_rd_n >>")')

    ! --- Structure_evolution ---
    if(imdalg == T_CONTROL .or. imdalg == VERLET .or. imdalg == VELOCITY_SCALING .or. &
    &  imdalg == PT_CONTROL .or. imdalg == P_CONTROL) then
       if( f_selectBlock( tag_structure_evolution) == 0) then
          ! --- Temperature_Control ---
          if( f_selectBlock( tag_temperature_control) == 0) then
             if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_temperature_control is found")')
             ! --- Method ---
             if( f_getStringValue( tag_method, rstr, LOWER) == 0) then
                call strncmp0(trim(rstr), tag_nose_hoover,tf)
                if(.not.tf) call strncmp0(trim(rstr), tag_hoover,tf)
                if(tf) then
                   t_ctrl_method = NOSE_HOOVER
                   goto 1001
                end if
                call strncmp0(trim(rstr),tag_nose, tf)
                if(tf) then
                   t_ctrl_method = NOSE
                   goto 1001
                end if
                call strncmp0(trim(rstr),tag_velocity_scaling, tf)
                if(tf) then
                   t_ctrl_method = VELOCITY_SCALING
                   goto 1001
                end if
                call strncmp0(trim(rstr),tag_langevin, tf)
                if(tf) then
                   t_ctrl_method = LANGEVIN
                   goto 1001
                end if
1001            continue
             else
                if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_method is not found")')
                t_ctrl_method = NOSE_HOOVER
             end if
             iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       endif
    endif
    ! --- Structure ---
    if( f_selectBlock( tag_structure) == 0) then
       if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** -- tag_structure --")')

       ! --- check whether the atoms can be excluded during MD ---
       if( f_getIntValue(tag_sw_atom_excludable,iret)==0) sw_atom_excludable=iret


       ! --- atom_list ---
       if( f_selectBlock( tag_atomlist) == 0) then
          ! --- count total number of atoms ---
          prealloc = .true.
          number_is_given = f_getIntValue( tag_numatoms, iret) == 0
          if(number_is_given) natm = iret

          if( f_getIntValue(tag_sw_initial_magmom_by_atom, iret ) == 0) then
             sw_set_initial_magmom_by_atom = iret
             write(nfout,*) "sw_set_initial_magmom_by_atom = ", &
                  &          sw_set_initial_magmom_by_atom
          endif

          if( f_getStringValue( tag_coordinate_system, rstr, LOWER) == 0) then
             call set_input_coordinate_system(rstr) ! -> input_coordinate_system
          end if

          if( f_selectBlock(tag_displacement) == 0) then
             call set_displacement()
             iret = f_selectParentBlock()
          end if
          if( f_selectBlock(tag_vibrational_mode) == 0) then
             call set_vibrational_mode()
             iret = f_selectParentBlock()
          end if
          if( f_selectBlock(tag_atoms) == 0) then
             call set_atompos_and_etc(prealloc, natm, iret)
             if(iret <= 0) call phase_error_with_msg(nfout,' atomic coordinates are not given properly <<m_IS_rd_n>>'&
                                                    ,__LINE__,__FILE__)
             if(number_is_given .and. natm > iret) then
                natm = iret
             else if(.not.number_is_given) then
                natm = iret
             end if
             natmorg = natm
             iret = f_selectParentBlock()
             write(nfout,'(" set_atompos_and_etc(.true.,...) is executed")')
             call flush(nfout)
          else
             call phase_error_with_msg(nfout,' tag_atom is not given <<m_IS_rd_n>>',__LINE__,__FILE__)
          end if

          if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** number of atoms = ",i5)') natm
          allocate(work(natm,3), stat=istat)
          if (istat /= 0) then
             if(printable) write(nfout,*) 'Allocation error in sub.set_atompos_and_etc',natm,istat
             call phase_error_with_msg(nfout,'Allocation error in sub.set_atompos_and_etc',__LINE__,__FILE__)
          end if

          call m_IS_alloc_pos_and_v(nfout) ! d(natm)
          call alloc_species_work()   ! d(natm)
          call alloc_species_vdw_work()   ! d(natm)

!          if( f_getStringValue( tag_coordinate_system, rstr, LOWER) == 0) then
!             call set_input_coordinate_system(rstr) ! -> input_coordinate_system
!          end if

          !if( f_getStringValue( tag_howtogive_coordinates, rstr, LOWER) == 0) then
          !   call set_howtogive_coordinates(rstr)   ! -> howtogive_coordinates
          !   if(ipriinputfile >= 0 .and. printable) &
          !        & write(nfout,'(" !* howtogive_coordinates = ",i5)') howtogive_coordinates
          !end if

          ! --- setting atomic coordinates and etc.
          if( f_selectBlock(tag_atoms) == 0) then
             prealloc = .false.
             call set_atompos_and_etc(prealloc, natm, iret) ! -> pos,imdtyp,element,natm2
             iret = f_selectParentBlock()
          else
             call phase_error_with_msg(nfout,' atom coordinates are not given properly in the inputfile',__LINE__,__FILE__)
          end if
          call count_species()        ! -> ntyp
          call count_species_vdw()    ! -> ntyp_vdw
          deallocate(work)
          iret = f_selectParentBlock()
       else
          if(coord_method == DIRECTIN) call phase_error_with_msg(nfout,' no atom_list <<m_IS_rd_n>>',__LINE__,__FILE__)
       end if

       if(mype_conf>0.and.driver==DRIVER_MTD)then
          if(f_selectBlock(tag_atomlist)==0)then
             write(rstr,*) mype_conf
             if(f_selectBlock(tag_atoms//trim(adjustl(rstr)))==0)then
!!                if(printable) write(nfout,'(a)') 'reading coodinates from : '//tag_atoms//trim(adjustl(rstr))
                call set_atompos_and_etc(prealloc,natm,iret)
                iret=f_selectParentBlock()
             endif
             iret = f_selectParentBlock()
          endif
       endif

       allocate(work(natm,3), stat=istat)
       if (istat /= 0) then
          if(printable) write(nfout,*) 'Allocation error in sub.set_atompos_and_etc',natm,istat
          call phase_error_with_msg(nfout,'Allocation error in sub.set_atompos_and_etc',__LINE__,__FILE__)
       end if
       if(sw_displace_atom == ON) call set_displacement2()
       if(sw_vibrational_mode == ON) call set_vibrational_mode2()
       call set_atompos2()

       deallocate(work)

       constraints_exist = .false.

       ! --- constraint ---
       call initialize_constraint_param     ! use this for cell optimization

       if( f_selectBlock( tag_constraint) == 0) then
          if(ipriinputfile >= 1) write(nfout,'(" !** -- tag_constraint is found --")')
          if(constraint_is_possible(imdalg)/=0) goto 1004

          ! --- bondlength_fix ---
          ! --- count total number of conditions of bondlength fix ---
          number_is_given = f_getIntValue( tag_num_fixed_bonds, iret) == 0
          if(number_is_given) num_fixed_bonds = iret

          tf = f_selectBlock( tag_fixed_bond) == 0
          if(.not.tf) tf = f_selectBlock( tag_fix_bondlength) == 0
          if(.not.tf) tf = f_selectBlock( tag_bondlength_fix) == 0
          if(tf) then
             if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** --tag_bondlength_fix or resembling tag is found --")')
             prealloc = .true.
             call set_fixed_bond_atoms(prealloc,num_fixed_bonds,iret,fixed_atoms)
             if(iret <= 0) call phase_error_with_msg(nfout,' fixed_bond is not given properly <<m_IS_rd_n>>', &
             __LINE__,__FILE__)
             nfcatm = fixed_atoms
             if(number_is_given .and. num_fixed_bonds > iret) then
                num_fixed_bonds = iret
             else if(.not.number_is_given) then
                num_fixed_bonds = iret
             end if
             if(ipriinputfile >= 1 .and. printable) &
                  & write(nfout,'(" !** number of fixed bonds = ",i5)') num_fixed_bonds
             iret = f_selectParentBlock()
          end if
          ! --- setting bondlength fix sets ---
          if (tf) call alloc_bondlength_fix_set()  ! allocate(bondlength_fix_set(2,num_fixed_bonds))
          tf = f_selectBlock( tag_fixed_bond) == 0
          if(.not.tf) tf = f_selectBlock( tag_fix_bondlength) == 0
          if(.not.tf) tf = f_selectBlock( tag_bondlength_fix) == 0
          if(tf) then
             prealloc = .false.
             call set_fixed_bond_atoms(prealloc,num_fixed_bonds,iret)
             constraints_exist = .true.
             constraint_type = BONDLENGTH_FIX

             call m_IS_alloc_cnstrvectors_etc(imdalg)
             call substitute_ia_cnst()
             iret = f_selectParentBlock()

             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !*  bondlength fix")')
                write(nfout,'(" !** number of fixed_atoms = ",i5)') fixed_atoms
                write(nfout,'(" !** nfcatm                = ",i5)') nfcatm
                do i = 1, nfcatm
                   write(nfout,'(" !** ia_cnst(",i5,") = ",i5)') i, ia_cnst(i)
                end do
             end if
          end if

          if(constraints_exist) goto 1003

          ! --- cog_fix_l ---
          ! --- count total number of planes in which cog's are fixed ---
!          number_is_given = f_getIntValue( tag_num_planes, iret) == 0
!          if(number_is_given) num_planes_atoms_are_fixed = iret
!!$          tf = f_selectBlock( tag_cog_fix) == 0
!!$          if(.not.tf) tf = f_selectBlock( tag_cog_fix_l) == 0
!!$          if(.not.tf) tf = f_selectBlock( tag_cog_fix_in_a_plane) == 0
!!$
!!$          if(tf) then
!!$             if(printable) write(nfout,'(" !* tag_cog_fix, tag_cog_fix_l or tag_cog_fix_in_a_plane is found")')
!!$             prealloc = .true.
!!$             call set_fixed_planes(prealloc,COG_FIX_L,num_planes_atoms_are_fixed,iret,nfcatm)
!!$             if(iret <= 0) stop ' cog_fix_planes are not given properly <<m_IS_rd_n>>'
!!$             if(number_is_given .and. num_planes_atoms_are_fixed > iret) then
!!$                num_planes_atoms_are_fixed = iret
!!$             else if(.not.number_is_given) then
!!$                num_planes_atoms_are_fixed = iret
!!$             end if
!!$             if(ipriinputfile >= 1 .and. printable) &
!!$                  & write(nfout,'(" !** number of planes in which COGs are fixed = ",i5)') num_planes_atoms_are_fixed
!!$             iret = f_selectParentBlock()
!!$          end if

! ===================== 2021/04/29   T. Yamasaki ===>>
          ! --- setting constrained atoms under center of gravity control, or as a rigid body,& ---
          ! ---   and normal vectors of the planes on which the center of gravity is fixed or & ---
          ! ---   rigid body center is fixed                                                    ---
          !       --- number of constrained atoms under center of gravity control ---
          constraint_type = 0
          tf = f_selectBlock( tag_cog_fix) == 0
          if(.not.tf) tf = f_selectBlock( tag_cog_fix_l) == 0
          if(.not.tf) tf = f_selectBlock( tag_cog_fix_in_a_plane) == 0
          if(tf) then
             if(ipriinputfile>=2) write(nfout,'(" !* tag_cog_fix_in_a_plane is found")')
             constraint_type = COG_FIX_L
             constraints_exist = .true.
             call set_number_of_constraint_atoms(constraint_type, iret) ! calling set_constraint_atoms --> iret
             nfcatm_cog = iret
             iret = f_selectParentBlock()
          end if
          !       --- number of constrained atoms dealt as a rigid body ---
          tf = f_selectBlock( tag_rigid_body)==0
          if(.not.tf) tf = f_selectBlock( tag_rigid_body_fix_l)==0
          if(.not.tf) tf = f_selectBlock( tag_rigid_body_fix_in_a_plane) == 0
          if(tf) then
             if(ipriinputfile>=DEBUGPRINTLEVEL) write(nfout,'(" !* tag_rigid_body_fix_in_a_plane is found")')
             constraint_type = constraint_type + RIGID_BODY_FIX_IN_A_PLANE
             constraints_exist = .true.
             call set_number_of_constraint_atoms(constraint_type, iret) ! calling set_constraint_atoms --> iret
             nfcatm_rb = iret
             iret = f_selectParentBlock()
          end if

          !       --- allocate arrays for constraint atoms ---
          select case(constraint_type)
          case (COG_FIX_L, RIGID_BODY_FIX_IN_A_PLANE, COG_and_RIGID_BODY_FIX_L)
             select case(constraint_type)
             case (COG_FIX_L)
                nfcatm = nfcatm_cog
             case (RIGID_BODY_FIX_IN_A_PLANE)
                nfcatm = nfcatm_rb
             case (COG_and_RIGID_BODY_FIX_L)
                nfcatm = nfcatm_cog + nfcatm_rb
             end select
             call m_IS_alloc_cnstrvectors_etc(imdalg) ! allocate ia_cnst, fcvect, ipfixedplane
          end select
          if(ipriinputfile>=1) write(nfout,'(" nfcatm_cog, nfcatm_rb, nfcatm = ",3i5)') nfcatm_cog, nfcatm_rb, nfcatm

          allocate(imdtyp_set(natm)); imdtyp_set = 0
          !       --- set constraint atoms under center of gravity control ---
          tf = f_selectBlock( tag_cog_fix) == 0
          if(.not.tf) tf = f_selectBlock( tag_cog_fix_l) == 0
          if(.not.tf) tf = f_selectBlock( tag_cog_fix_in_a_plane) == 0
          if(tf) then
             call set_constraint_atoms_imdtyp_etc(imdalg,nfcatm_cog,0,COG_FIX_L)  ! --> imdtyp, ia_cnst,nfcatm, ipfixedplane
             iret = f_selectParentBlock()
          end if
             !       --- set constraint atoms dealt as a rigid body ---
          tf = f_selectBlock( tag_rigid_body) == 0
          if(.not.tf) tf = f_selectBlock( tag_rigid_body_fix_l) == 0
          if(.not.tf) tf = f_selectBlock( tag_rigid_body_fix_in_a_plane) == 0
          if(tf) then
             call set_constraint_atoms_imdtyp_etc(imdalg,nfcatm_rb,nfcatm_cog,RIGID_BODY_FIX_IN_A_PLANE) ! --> imdtyp, ia_cnst,nfcatm, ipfixedplane
             iret = f_selectParentBlock()
          else
             if(ipriinputfile>=2) write(nfout,'(" !* tag_rigid_body_fix_in_a_plane is not found")')
          end if
          deallocate(imdtyp_set)

!!$          if(tf) iret = f_selectParentBlock()

!!$             call m_IS_alloc_cnstrvectors_etc(imdalg) ! allocate ia_cnst, fcvect
!!$             prealloc = .false.
!!$             if( f_selectBlock(tag_atoms)==0) then
!!$                call set_constraint_atoms(prealloc, nfcatm, 0,constraint_type, iret) ! --> imdtyp, ia_cnst,nfcatm
!!$                iret = f_selectParentBlock()
!!$             end if

          !       --- set fixed_plane for constraint atoms ---
          prealloc = .true.
          if(f_selectBlock( tag_fixed_plane) == 0) then
             if(ipriinputfile>=1) write(nfout,'(" !* tag_fixed_plane is found")')
             call set_constraint_plane(prealloc,nfcatm,iret) ! ia_cnst --> iret
             num_planes_atoms_are_fixed = iret
             if(ipriinputfile>=1) write(nfout,'(" !* num_planes_atoms_are_fixed = ",i8)') num_planes_atoms_are_fixed
             if(num_planes_atoms_are_fixed == 0) then
                stop ' num_planes_which_atoms_are_fixed = 0'
             end if
             iret = f_selectParentBlock()
          end if

          call check_number_of_planes_atoms_are_fixed(nfcatm_rb,nfcatm_cog,num_planes_atoms_are_fixed,constraint_type)
                                       ! -> num_planes_atoms_are_fixed_rb, num_planes_atoms_are_fixed_cog
!!$          call alloc_distances_of_planes() ! num_planes_atoms_are_fixed_rb, num_planes_atoms_are_fixed_cog ->

          prealloc = .false.
          if(f_selectBlock( tag_fixed_plane) == 0) then
             call set_constraint_plane(prealloc,nfcatm,iret) ! ia_cnst --> fcvect
             iret = f_selectParentBlock()
          end if

          if(f_getIntValue(tag_sw_apply_constant_force,iret)==0) then
             apply_constant_force = iret
             if(f_getRealValue(tag_force_apply,dret,'hartree/bohr')==0) then
                force_apply = dret
             end if
             if(abs(force_apply) < SmallestPositiveNumber*1.d5) apply_constant_force = OFF
          else
             apply_constant_force = OFF
          end if

          if(f_getIntValue(tag_sw_apply_constant_move, iret)==0) then
             apply_constant_move = iret
          end if

          if(nfcatm_rb >=1) then
             tf = f_getIntValue(tag_sw_rigid_body_cog_movement,iret)==0
             if(.not.tf) tf = f_getIntValue(tag_sw_relax_rigid_body_cog,iret)==0
             if(tf) then
                relax_rigid_body_cog = iret
             else
                relax_rigid_body_cog = OFF
             end if
          else
             relax_rigid_body_cog = OFF
          end if

          if(nfcatm>=1) call set_pointer_array_to_constrained_atoms(nfcatm) ! --> icnst_a

          move_constrained_plane = .false.
          do i = 1, nfcatm
             if(max_reach_of_fixed_plane<=fcvect(i,5)) max_reach_of_fixed_plane = fcvect(i,5)
             if(abs(fcvect(i,4))>SmallestPositiveNumber*1.d5) move_constrained_plane = .true.
          end do

          if(relax_rigid_body_cog == ON) then
             apply_constant_force = OFF
             apply_constant_move = OFF
             move_constrained_plane = .false.
          end if

          if(apply_constant_move == ON) apply_constant_force = OFF
          if(ipriinputfile >= 1) then
             write(nfout,'(" relax_rigid_body_cog = ",i8)') relax_rigid_body_cog
             write(nfout,'(" apply_constant_force = ",i8,", force_apply = ",f8.4)') apply_constant_force, force_apply
             write(nfout,'(" apply_constant_move  = ",i8)') apply_constant_move
             write(nfout,'(" constraint_type      = ",i8)') constraint_type
             write(nfout,'(" !** nfcatm = ",i5)') nfcatm
             write(nfout, &
             '("  !** ia  ia_cnst  fcvect(ia,1:3) increm  max_reach fcvect(ia,6:8) ipfixedplane icount_of_ipfixedplane")')
             do i = 1, nfcatm
                write(nfout,'(" !** ",2i5, 3f7.3, 2f8.4, 3f7.3, 2i3, " (2)")') &
                     & i,ia_cnst(i), fcvect(i,1:8), ipfixedplane(i), icount_of_ipfixedplane(i)
!!$                write(nfout,'(" !** ia_cnst(",i4,") = ",i5," , fcvect = ",3f7.3, " increm. = ",f6.3 &
!!$                     & , " max_reach = ",f8.4, " fcvect(6:8) = ",3f8.4," ipfixedplane = ",i2," icount_of_ipfixedp. = ",i2, " (2)")') &
!!$                     & i,ia_cnst(i), fcvect(i,1:4), fcvect(i,5), fcvect(i,6:8), ipfixedplane(i), icount_of_ipfixedplane(i)
             end do
             call flush(nfout)
          end if
!!$          stop '<<m_IS_rd_n>>'
!!$             call m_IS_alloc_rigid_body_vectors_etc(imdalg,nfcatm)
!!$                       if(iprimd >= 1) write(nfout,'(" !* nfcatm = ",i8)') nfcatm
!!$             call set_rigid_body_fixed_plane(RIGID_BODY_FIX_IN_A_PLANE)

!!$          iret = f_selectParentBlock()
          if(constraints_exist) goto 1003
! <<===================== 2021/04/29   T. Yamasaki

          ! --- fix_in_a_plane ---
          ! --- count total number of planes in which atoms are fixed ---
          prealloc = .true.
          number_is_given = f_getIntValue( tag_num_planes, iret) == 0
          if(number_is_given) num_planes_atoms_are_fixed = iret
          tf = f_selectBlock( tag_fix_in_a_plane) == 0
          if(.not.tf) tf = f_selectBlock( tag_fixed_plane) == 0
          if(tf) then
             if(ipriinputfile>=1) write(nfout,'(" !* tag_fix_in_a_plane or tag_fixed_plane is found")')
             call set_fixed_planes(prealloc,FIX_IN_A_PLANE,num_planes_atoms_are_fixed,iret, nfcatm)
             if(iret <= 0) call phase_error_with_msg(nfout,' planes are not given properly <<m_IS_rd_n>>',__LINE__,__FILE__)
             if(number_is_given .and. num_planes_atoms_are_fixed > iret) then
                num_planes_atoms_are_fixed = iret
             else if(.not.number_is_given) then
                num_planes_atoms_are_fixed = iret
             end if
             iret = f_selectParentBlock()
             if(ipriinputfile >= 1) &
               & write(nfout,'(" !** number of planes in which atoms are fixed = ",i5)') num_planes_atoms_are_fixed
          end if

          ! --- setting fixed planes and constrained atoms ---
!!$          call alloc_fcvect_tmp()
          tf = f_selectBlock( tag_fix_in_a_plane) == 0
          if(.not.tf) tf = f_selectBlock( tag_fixed_plane) == 0
          if(tf) then
             call m_IS_alloc_cnstrvectors_etc(imdalg)
             prealloc = .false.
             call set_fixed_planes(prealloc, FIX_IN_A_PLANE, num_planes_atoms_are_fixed, iret)
             constraints_exist = .true.
             constraint_type = FIX_IN_A_PLANE
             iret = f_selectParentBlock()

             move_constrained_plane = .false.
             do i = 1, nfcatm
                if(abs(fcvect(i,4)) > SmallestPositiveNumber*1.d5) move_constrained_plane = .true.
             end do

             if(ipriinputfile >= 1 .and. printable) then
                write(nfout,'(" !** nfcatm = ",i5)') nfcatm
                do i = 1, nfcatm
                   write(nfout,'(" !** ia_cnst(",i4,") = ",i5," , &
                        & fcvect = ",4f8.3," max_reach = ",f8.4, " ipfixedplane = ",i5," (4)")') &
                        & i,ia_cnst(i),fcvect(i,1:4), fcvect(i,5),ipfixedplane(i)
                end do
             end if
          end if

          if(constraints_exist) goto 1003

          ! --- setting a fixed normal hypervector ---
          tf = f_selectBlock( tag_fixed_normal_hypervector) == 0
          if(tf) then
             call alloc_normal_hypervector()

! ======== KT_mod ==== 2014/07/14
!             if(endpoint_images /= NO) then
!                call set_normal_hypervector()
!             else
!                stop ' endpoint_images are not given'
!             end if

             if ( f_getIntValue( tag_hypervec_from_file,iret )== 0 ) then
                write(nfout,*) '*** hypervec_from_file is ', iret
                hypervec_from_file = iret
             endif
             if ( hypervec_from_file == ON ) then
                call m_Files_open_nfhypervec()
                call read_hypervec_from_file
                call m_Files_close_nfhypervec()
             else
                if(endpoint_images /= NO) then
                   call set_normal_hypervector()
                else
                   call phase_error_with_msg(nfout,' endpoint_images are not given',__LINE__,__FILE__)
                end if
             endif
! ===================== 2014/07/14

!!$             tf = f_getIntValue(tag_mode_coefficient,iret)== 0
!!$             if(tf) mode_fi_coefficient = iret
!!$             if(mode_fi_coefficient == ON) then
!!$                tf = f_getRealValue(tag_coefficient,dret,'') == 0
!!$                if(tf) fi_coefficient = dret
!!$             end if
             iret = f_selectParentBlock()
             constraints_exist = .true.
             constraint_type = FIXED_NORMAL_HYPERVECTOR

             if(ipriinputfile >= 1) then
                call wd_normal_hypervector()
!!$                write(nfout,'(" !** mode_fi_coefficient = ",i5)') mode_fi_coefficient
!!$                write(nfout,'(" !**      fi_coefficient = ",f8.4)') fi_coefficient
             end if
          end if

1003      continue
          if(ipriinputfile >=1 ) then
             if(move_constrained_plane) then
                write(nfout,'(" !** move_constrained_plane = .true.")')
             else
                write(nfout,'(" !** move_constrained_plane = .false.")')
             end if
          end if

          if(.not.constraints_exist .and. ipriinputfile>=1) &
               & write(nfout,'(" !** constraint details are not described")')
1004      iret = f_selectParentBlock()
       else
          if(ipriinputfile >= 1) write(nfout,'(" !** -- tag_constraint is not found --")')
          if(imdalg == QUENCHED_CONSTRAINT) then
             imdalg = QUENCHED_MD
             if(ipriinputfile >= 1 ) write(nfout,&
                  & '(" !** imdalg is set QUENCHED_MD (from QUENCHED_CONSTRAINT)")')
          end if
       end if

       if(sw_vdw_correction == ON .and. vdw_method == VDW_DFTD3)then
          call m_IS_alloc_vdwdf3()
          call resolve_dftd3_parameters()
          call read_dftd3_parameters()
       endif
       ! ---- vdw_list ----
       if(vdw_method /= VDW_DFTD3 .and. ntyp_vdw>0) then
          call alloc_speciesname_vdw()     ! d(ntyp_vdw)
          call m_IS_alloc_vdw() ! cvdw, rvdw, ityp_vdw, fxyzvdw_l
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!         if( f_selectBlock( tag_vdw_list) == 0) then
!            call set_vdw_parameters()   ! -> cvdw, rvdw
!            iret = f_selectParentBlock()
!         else
!            stop ' the vdW list is not given properly in the inputfile'
!         end if
          call set_vdw_parameters()   ! -> cvdw, rvdw
! ==============================================================================
          call specify_ityp_vdw()  ! -> itpy_vdw
          call dealloc_species_vdw_work()
       end if

       ! --- element_list ---
       call alloc_speciesname()     ! d(ntyp)
       call m_IS_alloc_iatomn_etc() ! iatomn, iatom, ivan, alfa, amion, zeta1, qex
#ifdef forsafe
       goto 100
#endif
       call initial_species_name()
       call set_iatomn_default(LEN_ATOMNAME,ntyp,speciesname,iatomn)
       tf = .true.
       do i=1,ntyp
          if(iatomn(i) < 1 )then
            tf = .false.
            exit
          endif
       enddo
       if(tf) then
         amion = -100.0
         call set_amion_default(ntyp,iatomn,amion)
       endif

100    continue

       if( f_selectBlock( tag_element_list) == 0) then
          call set_element_detail()   ! -> iatomn, alfa, amion, zeta1,qex
          iret = f_selectParentBlock()
       else if (.not. tf) then
          call phase_error_with_msg(nfout,' the element list is not given properly in the inputfile',__LINE__,__FILE__)
       end if
       call specify_ityp()  ! -> itpy
#ifdef __EDA__
       if(sw_eda==ON) then
! -----  ascat starts modifying  -----
       call specify_icount_EDA()
!org   if(ipriinputfile >= 1 .and. printable) call wd_atom_list()
       call wd_atom_list()
! -----  ascat ceases modifying  -----
       endif
#endif
       call m_IS_set_iproj_group()
       if(ipriinputfile >= 1 .and. printable) call wd_atom_list()
       call m_IS_set_iatom(nfout) ! -> iatom
       call dealloc_species_work()

#ifndef _EMPIRICAL_
       !! Projector: set proj_attribute(:)%ityp
       do i=1,natm
          call m_CtrlP_set_proj_ityp(iproj_group(i),ityp(i))
       end do
#endif
       if(sw_atom_excludable==ON)then
           if(sum(exclusion_target)==0)then
               if(printable)then
                   write(nfout,'(a)') ' !** sw_atom_excludable is ON, but exclusion target is undefined.'
               endif
               sw_atom_excludable = OFF
           endif
       endif
       if(sw_atom_excludable==ON)then
           call set_atm_exclusion_criteria()
           if(printable)then
              write(nfout,'(a)') ' !** atoms are excludable during MD simulation '
              write(nfout,'(a)') ' !** atoms targeted for exclusion '
              do i=1,natm
                 if(exclusion_target(i)==ON) write(nfout,'(i8)') i
              enddo
              write(nfout,'(a)') ' !** an atom will be excluded if it is located outside the following box:'
              write(nfout,'(a,e20.10,a,e20.10,a,e20.10)') &
           & ' !** x greater than ',exclusion_criteria_min(1),' and less than ',exclusion_criteria_max(1)
              write(nfout,'(a,e20.10,a,e20.10,a,e20.10)') &
           & ' !** y greater than ',exclusion_criteria_min(2),' and less than ',exclusion_criteria_max(2)
              write(nfout,'(a,e20.10,a,e20.10,a,e20.10)') &
           & ' !** z greater than ',exclusion_criteria_min(3),' and less than ',exclusion_criteria_max(3)
           endif
       endif

       if(f_selectBlock(tag_reservoir)==0.and.icond/=COORDINATE_CONTINUATION) then
          call set_atompos_and_etc_reservoir()
          iret = f_selectParentBlock()
       endif

       ! --- atom addition criteria
       if(f_selectBlock(tag_atom_addition_criteria)==0)then
          if(f_getIntValue(tag_sw_rotate_reservoir,iret)==0) sw_rotate_reservoir = iret
          if(f_getIntValue(tag_frequency,iret)==0) addition_frequency = iret
          iret = f_selectParentBlock()
       endif
       !if(printable .and. addition_frequency>0) write(nfout,'(a,i8)') '!** addition frequency : ',addition_frequency
       if(printable) write(nfout,'(a,i8)') ' !** addition frequency : ',addition_frequency
       pos_in = pos;cps_in=cps

       if(f_getIntValue(tag_sw_check_coordinates,iret)==0) sw_check_coordinates = iret
       if(f_getRealValue(tag_bond_length_min,dret,'bohr')==0) bond_length_min = dret
       !check validity of atomic coordinates

       ! read and allocate 'regions'
       call set_regions(nfout)

       if(sw_check_coordinates==ON) call check_coor()

       if(coord_method == FILE) then
       if(filetype == phase0_output)then
         if(.not.read_from_file)then
           call import_from_dynm_pre(nfpos,F_POS,'F_POS',ntyp,natm)
           if(.not.allocated(cps)) allocate(cps(natm,3))
           if(.not.allocated(pos)) allocate(pos(natm,3))
           if(.not.allocated(cpd_l)) allocate(cpd_l(natm,3));cpd_l=0.d0
           if(.not.allocated(ityp)) allocate(ityp(natm))
           if(.not.allocated(speciesname)) allocate(speciesname(ntyp))
           call import_from_dynm(nfpos,F_POS,'F_POS',frame,natm,ntyp,cps,pos,cpd_l &
           &    ,ityp,speciesname,altv,rltv)
           read_from_file = .true.
         endif
       endif
       endif

       if(f_getIntValue(tag_sw_mobility_by_distance,iret)==0) then
          sw_mobility_by_distance = iret
          if(sw_mobility_by_distance == ON)then
            if(f_selectBlock(tag_mobility_by_distance)==0)then
              if(f_getIntValue(tag_target_atom,iret)==0)then
                target_atom = iret
              endif
              if(f_getRealValue(tag_target_posx,dret,'bohr')==0)then
                target_pos(1) = dret
              endif
              if(f_getRealValue(tag_target_posy,dret,'bohr')==0)then
                target_pos(2) = dret
              endif
              if(f_getRealValue(tag_target_posz,dret,'bohr')==0)then
                target_pos(3) = dret
              endif
              if(f_getRealValue(tag_distance,dret,'bohr')==0)then
                distance = dret
              endif
              if(printable) then
                 if(target_atom>0) then
                   write(nfout,'(a,i8,a,f10.5,a)') ' !** atoms whose distance from atom ' &
               &  ,target_atom,' is less than or equal to ',distance,' bohr will be set as mobile atoms'
                 else
                   write(nfout,'(a,3f10.5,a,f10.5,a)') ' !** atoms whose distance from position ' &
               &  ,target_pos(1:3),' is less than or equal to ',distance,' bohr will be set as mobile atoms'
                 endif
              endif
              iret = f_selectParentBlock()
            endif
            call set_mobility_by_distance(target_atom,target_pos,distance)
          endif
       endif

       ! setup rigid bodies
       if(f_selectBlock(tag_constraint)==0) then
         if(f_getIntValue(tag_sw_full_rb_dynamics,iret)==0) then
           sw_full_rb_dynamics = iret
         endif
         if(sw_full_rb_dynamics==ON) then
           if(f_selectBlock(tag_rigid_body)==0) then
             call read_rigid_body_input2()
             iret = f_selectParentBlock()
           endif
         endif
         if(nrigid_bodies>0) constraint_type = 0
         iret = f_selectParentBlock()
       endif

       if(nrigid_bodies==0) then
         if(f_selectBlock(tag_atomlist)==0) then
           call read_rigid_body_input()
           iret = f_selectParentBlock()
         endif
       endif
       if(nrigid_bodies>0) then
         call initialize_rigid_bodies()
         call read_rigid_body_options()
         call initialize_zaxis()
         call initialize_rotation_matrix()
         call print_rigid_body_status()
       endif

       if(f_selectBlock(tag_fix_bond)==0) then
         if(f_getIntValue(tag_sw_fix_bond,iret)==0) then
           sw_fix_bond = iret
         endif
         if(sw_fix_bond == ON) then
           call read_fix_bond_input()
           call set_covrads()
           call read_fix_bond_options()
           call resolve_bonds_to_be_fixed()
           call print_bonds_to_be_fixed()
!!$           call map_constraint()
         endif

         iret = f_selectParentBlock()
       endif

       iret = f_selectParentBlock()
    end if


    ! --- multiple replica for new format --
    if(multiple_replica_mode == ON) then

      if( f_selectBlock( tag_multiple_replica) == 0) then
      if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !*  tag_multiple_replica")')

        if( f_getStringValue(tag_method, rstr, LOWER) == 0) then
          if(rstr == tag_dimer_method) then
            sw_dimer_method = ON
          endif
        endif
        if( f_selectBlock(tag_dimer_method)==0) then
          if( f_getRealValue(tag_delta_theta,dret,'') == 0) delta_theta = dret
          if( f_getRealValue(tag_delta_r,dret,'bohr') == 0) delta_r = dret
          if( f_getIntValue(tag_sw_random_unit_vector, iret) ==0) &
                  sw_random_unit_vector = iret
          if( f_getRealValue(tag_to_dimer_threshold, dret,'')==0) to_dimer_thres = dret
          if(printable) write(nfout,'(a,3f10.5)') '!** DIMER method delta_theta, delta_r, to_dimer_thres ', &
                                            delta_theta, delta_r, to_dimer_thres
          iret = f_selectParentBlock()
        endif
        if( f_selectBlock( tag_structure) == 0) then
          if (f_getIntValue(tag_sw_path_from_dynm,iret)==0) sw_path_from_dynm = iret
          if (f_getRealValue(tag_end0_energy,dret,'hartree')==0)then
             end0_energy = dret;end0_energy_given = .true.
          endif
          if (f_getRealValue(tag_end1_energy,dret,'hartree')==0)then
             end1_energy = dret;end1_energy_given = .true.
          endif
          if (f_getIntValue(tag_frame_end0,iret)==0) frame_end0 = iret
          if (f_getIntValue(tag_frame_end1,iret)==0) frame_end1 = iret
        ! -- count total number of replicas
          prealloc = .true.
          number_is_given = f_getIntValue( tag_number_of_replicas, iret) == 0
          if(number_is_given) number_of_replicas = iret
          if( f_selectBlock(tag_replicas)  == 0) then
            call set_replica_input_method(prealloc, number_of_replicas, iret)
            if(number_is_given .and. number_of_replicas > iret) then
               number_of_replicas = iret
            else if(.not.number_is_given) then
               number_of_replicas = iret
            end if
            iret = f_selectParentBlock()
!          else
!            stop ' tag_replicas is not given <<m_IS_rd_n>>'
          end if
          if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** number of replicas = ",i5)') number_of_replicas
          ! <--

          if(number_of_replicas > 1) then
            allocate(replica_howtogive_coordinates(number_of_replicas), stat=istat)
            if (istat /= 0) then
              if(printable) write(nfout,*) 'Allocation error of replica_howtogive_coordinates in <<m_IS_rd_n>>'&
                   & ,number_of_replicas,istat
                call phase_error_with_msg(nfout,'Allocation error of replica_howtogive_coordinates in <<m_IS_rd_n>>'&
                                         ,__LINE__,__FILE__)
            end if
            replica_howtogive_coordinates(:) = PROPORTIONAL
            allocate(replica_endpoints(2,number_of_replicas), stat=istat)
            if (istat /= 0) then
              if(printable) write(nfout,*) 'Allocation error of replica_endpoints in <<m_IS_rd_n>>' &
                   & ,number_of_replicas,istat
              call phase_error_with_msg(nfout,'Allocation error of replica_endpoints in <<m_IS_rd_n>>',__LINE__,__FILE__)
            end if
            replica_endpoints(1,:) = 0
            replica_endpoints(2,:) = -1

            ! 0 represents the right-most image, while -1represents the
            ! left-most image
            ! --- setting replica_howtogive_coordinates and endpoints --
            if( f_selectBlock(tag_replicas)  == 0) then
              prealloc = .false.
              call set_replica_input_method(prealloc, number_of_replicas, iret)
                          ! -> replica_howtogive_coordinates, replica_endpoints
              iret = f_selectParentBlock()
            !else
            !  stop ' replicas are not given properly in the inputfile'
            end if
          else
            number_of_replicas = 1
            allocate(replica_howtogive_coordinates(number_of_replicas), stat=istat)
            allocate(replica_endpoints(2,number_of_replicas), stat=istat)
            replica_howtogive_coordinates(1) = DIRECTIN
          end if

          if(ipriinputfile >= 1 .and. printable) then
            if(number_of_replicas > 1) then
              write(nfout,'(" !**   ---  howtogive_coordinates ---")')
              write(nfout,'(" !**         DIRECTIN : ",i4)') DIRECTIN
              write(nfout,'(" !**     PROPORTIONAL : ",i4)') PROPORTIONAL
              write(nfout,'(" !**   FROM_ENDPOINTS : ",i4)') FROM_ENDPOINTS
              write(nfout,'(" !**            FILE  : ",i4)') FILE
            end if
            do i = 1, number_of_replicas
              write(nfout,'(" !** id = ",i8, " howtogive_coordinates = ",i4, " endpoints = ",2i5 )') &
      &          i, replica_howtogive_coordinates(i), replica_endpoints(1:2,i)
            end do
          end if

          if(number_of_replicas > 1) then
            tf = .false.
            do i = 1, number_of_replicas
              if(replica_howtogive_coordinates(i) == FILE) then
                tf = .true.
                exit
              end if
            end do
            if(tf) then
              allocate(pos_image(number_of_replicas,natm,3))
              allocate(cps_image(number_of_replicas,natm,3))
              if(printable) write(nfout,'(" !** pos_image and cps_image are allocated. ")')
              do i = 1, number_of_replicas
                if(replica_howtogive_coordinates(i) == FILE) then
                  input_coordinate_system_t = input_coordinate_system
                  call m_Files_open_nfimage(i)
                  call set_endpoint_atompos_from_file(nfimage,natm, &
                       pos_image(i,:,:),cps_image(i,:,:))
                   close(nfimage)
                  input_coordinate_system = input_coordinate_system_t
                end if
              end do
            end if
          end if

          endpoint_images = DIRECTIN
          if( f_getStringValue( tag_endpoint_images, rstr, LOWER) == 0) &
               & call set_endpoint_images(rstr)  ! -> endpoint_images
          if(ipriinputfile >= 1) write(nfout,'(" !** set_endpoint_images = ",i5)') endpoint_images

          if( f_getStringValue( tag_howtogive_coordinates, rstr, LOWER) == 0) then
            call set_howtogive_coordinates(rstr)   ! -> howtogive_coordinates
            if(ipriinputfile >= 0 .and. printable) &
                & write(nfout,'(" !* howtogive_coordinates = ",i5)') howtogive_coordinates
          end if

          allocate(work(natm,3), stat=istat)
          if (istat /= 0) then
            if(printable) write(nfout,*) 'Allocation error in sub.set_atompos_and_etc',natm,istat
            call phase_error_with_msg(nfout,'Allocation error in sub.set_atompos_and_etc',__LINE__,__FILE__)
          end if

          call alloc_endpoint_pos(nfout)
          if(endpoint_images == DIRECTIN) then
            input_coordinate_system_t = input_coordinate_system
            pos_end0 = pos;cps_end0=cps
            if( f_selectBlock( tag_atom_list_end0) == 0) then
              if( f_getStringValue( tag_coordinate_system, rstr, LOWER) == 0)  &
                call set_input_coordinate_system(rstr) ! -> input_coordinate_system
              if( f_selectBlock(tag_atoms) == 0) then
                prealloc = .false.
                call set_endpoint_atompos(natm,pos_end0, cps_end0) ! -> pos,imdtyp,element,natm2
                iret = f_selectParentBlock()
              else
                call phase_error_with_msg(nfout,' endpoint0 atom coordinates are not given properly in the inputfile'&
                                         ,__LINE__,__FILE__)
              end if
              iret = f_selectParentBlock()
            end if

            pos_end1 = pos;cps_end1=cps
            if( f_selectBlock( tag_atom_list_end1) == 0) then
              if( f_getStringValue( tag_coordinate_system, rstr, LOWER) == 0)  &
                call set_input_coordinate_system(rstr) ! -> input_coordinate_system
              if( f_selectBlock(tag_atoms) == 0) then
                prealloc = .false.
                call set_endpoint_atompos(natm,pos_end1, cps_end1) ! -> pos,imdtyp,element,natm2
                iret = f_selectParentBlock()
              else
                call phase_error_with_msg(nfout,' endpoint1 atom coordinates are not given properly in the inputfile'&
                                         ,__LINE__,__FILE__)
              end if
              iret = f_selectParentBlock()
            end if
            input_coordinate_system = input_coordinate_system_t
          else if(endpoint_images == FILE) then
            allocate(cpd_l_t(natm,3))
            allocate(ityp_t(natm))
            allocate(speciesname_t(ntyp))
            call import_from_dynm(nfimage,F_IMAGE(0),'F_IMAGE(1)',frame_end0, &
            &    natm,ntyp,cps_end0,pos_end0,cpd_l_t,ityp_t,speciesname_t,altv_t,rltv_tt)
            call import_from_dynm(nfimage,F_IMAGE(-1),'F_IMAGE(-1)',frame_end1, &
            &    natm,ntyp,cps_end1,pos_end1,cpd_l_t,ityp_t,speciesname_t,altv_t,rltv_tt)
            deallocate(cpd_l_t)
            deallocate(ityp_t)
            deallocate(speciesname_t)
          end if
              !!$call set_atompos_from_endpoints() ! pos <- pos_end0, pos_end1

          if (sw_path_from_dynm == ON)then
            allocate(cpd_l_t(natm,3))
            allocate(ityp_t(natm))
            allocate(speciesname_t(ntyp))
            if(.not.allocated(pos_end0)) call alloc_endpoint_pos(nfout)
            if(.not.allocated(cps_image)) allocate(cps_image(number_of_replicas,natm,3))
            if(.not.allocated(pos_image)) allocate(pos_image(number_of_replicas,natm,3))
            call import_from_dynm(nfimage,F_PATH,'F_PATH',1,natm,ntyp &
            &    ,cps_end0,pos_end0,cpd_l_t,ityp_t,speciesname_t,altv_t,rltv_tt)
            do i=1,number_of_replicas
            call import_from_dynm(nfpath,F_PATH,'F_PATH',i+1,natm,ntyp &
            &    ,cps_image(i,:,:),pos_image(i,:,:),cpd_l_t,ityp_t,speciesname_t,altv_t,rltv_tt)
            enddo
            call import_from_dynm(nfpath,F_PATH,'F_PATH',number_of_replicas+2,natm,ntyp &
            &    ,cps_end1,pos_end1,cpd_l_t,ityp_t,speciesname_t,altv_t,rltv_tt)
            deallocate(cpd_l_t)
            deallocate(ityp_t)
            deallocate(speciesname_t)
          endif
          deallocate(work)
          iret = f_selectParentBlock()
        end if   ! tag_structure

        if( f_selectBlock( tag_accuracy) == 0) then
           if(f_getRealValue(tag_neb_dt,dret,'au_time') == 0) neb_dt = dret
           if( f_getStringValue(tag_neb_time_integral, rstr, LOWER) == 0 ) then
             call set_neb_time_integral(rstr)
           end if
           if( f_getStringValue(tag_dimer_time_integral, rstr, LOWER) == 0 ) then
             call set_dimer_time_integral(rstr)
           end if
           if(f_getRealValue(tag_sd_factor,dret,'')==0) sd_factor = dret
           if(f_getIntValue(tag_penalty_function, iret) == 0) penalty_function = iret
!!            if(f_getIntValue(tag_neb_convergence_condition,iret) == 0) &
!!            neb_convergence_condition = iret
           if( f_getStringValue(tag_neb_convergence_condition, rstr, LOWER) == 0 ) then
             call set_neb_convergence_condition(rstr)
           end if
           if(f_getRealValue(tag_neb_convergence_threshold,dret,'') == 0) &
              neb_convergence_threshold = dret
           if(f_getRealValue(tag_dimer_convergence_threshold,dret,'') == 0) &
              neb_convergence_threshold = dret
              iret = f_selectParentBlock()
        end if   ! tag_constraint

        if( f_selectBlock( tag_constraint) == 0) then
          if(f_getIntValue(tag_ci_neb, iret) == 0) ci_neb = iret
          if(f_getIntValue(tag_ci_index, iret) == 0) ci_index = iret
          ci_thres = neb_convergence_threshold*2.d0
          if(f_getRealValue(tag_ci_thres, dret,'') == 0) ci_thres = dret
          if(f_getRealValue(tag_sp_k_init,dret,'') == 0) sp_k_init = dret
          if(f_getRealValue(tag_sp_k_min,dret,'') == 0) sp_k_min = dret
          if(f_getRealValue(tag_sp_k_max,dret,'') == 0) sp_k_max = dret
          if(f_getIntValue(tag_sp_k_variable, iret) == 0) sp_k_variable = iret
          iret = f_selectParentBlock()
        end if   ! tag_constraint

        iret = f_selectParentBlock()
      end if  ! tag_multiple_replica
    end if   ! multiple_replica_mode

    if(nfcatm >= 1) then
       cnst_typ = imdtyp(ia_cnst(1))
       imdtypxyz(ia_cnst(1),:) = imdtyp(ia_cnst(1))
       if(cnst_typ > HEAT_BATH) then
          cnst_typ = cnst_typ - HEAT_BATH
       end if
       if(imdalg == QUENCHED_CONSTRAINT) then
          if(cnst_typ /= COG_FIX_L .and. cnst_typ /= BONDLENGTH_FIX_1 &
               & .and. cnst_typ /= BONDLENGTH_FIX_2 .and. cnst_typ /= RIGID_BODY_FIX_IN_A_PLANE) then
             imdalg = QUENCHED_MD
             if(iprimd >= 1) write(nfout,'(" !** imdalg is reset ",i5 &
                  & ," (= QUENCHED_MD), because cnst_typ is not proper")') imdalg
          end if
       end if

       if(iprimd>=1) then
          cnstrainttype: select case (cnst_typ)
          case(COG_FIX)
             write(nfout,'(" !** cnst_typ = ",i5," (= COG_FIX)")') cnst_typ
          case (COG_FIX_L)
             write(nfout,'(" !** cnst_typ = ",i5," (= COG_FIX_L)")') cnst_typ
          case (FIX_IN_A_PLANE)
             write(nfout,'(" !** cnst_typ = ",i5," (= FIX_IN_A_PLANE)")') cnst_typ
          case (BONDLENGTH_FIX)
             write(nfout,'(" !** cnst_typ = ",i5," (= BONDLENGTH_FIX)")') cnst_typ
          case (BONDLENGTH_FIX_1)
             write(nfout,'(" !** cnst_typ = ",i5," (= BONDLENGTH_FIX_1)")') cnst_typ
          case (BONDLENGTH_FIX_2)
             write(nfout,'(" !** cnst_typ = ",i5," (= BONDLENGTH_FIX_2)")') cnst_typ
          case (FIXED_NORMAL_HYPERVECTOR)
             write(nfout,'(" !** cnst_typ = ",i5," (= FIXED_NORMAL_HYPERVECTOR)")') cnst_typ
          case (RIGID_BODY_FIX_IN_A_PLANE)
             write(nfout,'(" !** cnst_typ = ",i5," (= RIGID_BODY_FIX_IN_A_PLANE)")') cnst_typ
          case default
             write(nfout,'(" !** cnst_typ = ",i5)') cnst_typ
          end select cnstrainttype
       end if
       allocate(fcvect_work(natm,8)); fcvect_work = 0.d0
       if(cnst_typ == RIGID_BODY_FIX_IN_A_PLANE) then
          do ig =1, nfcatm
!!$             fcvect_work(ig,1:4) = rigid_body_vect(1:4)
             fcvect_work(ig,1:8) = fcvect(ig,1:8)
          end do
       else
          fcvect_work(1:nfcatm,1:8) = fcvect(1:nfcatm,1:8)
       end if
       call m_IS_init_cnstrnt(natm,fcvect_work) ! -> sgmc
       deallocate(fcvect_work)
    end if

    ! --- Structure_evolution ---
    if(imdalg == T_CONTROL .or. imdalg == VERLET .or. imdalg == VELOCITY_SCALING .or. &
    &  imdalg == PT_CONTROL .or. imdalg == P_CONTROL) then
       mobility_cell(:,:) = ON
       iret = f_selectTop()
       if( f_selectBlock( tag_structure_evolution) == 0) then
          if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_structure_evolution is found")')
          ! --- Temperature_Control ---
          if( f_selectBlock( tag_temperature_control) == 0) then
             if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_temperature_control is found")')
             ! --- Method ---
             !     was read prior to this region
             !
             ! --- count total number of thermostat ---
             ! --- num_thermostat ---

             number_is_given = f_getIntValue( tag_num_thermostat, iret) == 0
             if(.not.number_is_given) &
                  & number_is_given = f_getIntValue( tag_num_thermo,iret) == 0
             if(number_is_given)  nrsv = iret

             if(t_ctrl_method == NOSE) nrsv = 1
             if(printable) write(nfout,'(" !** t_ctrl_method = ",i5)') t_ctrl_method
             if(printable) write(nfout,'(" !** num_Tresrvoir = ",i5)') nrsv

             if(f_getIntValue(tag_num_chain,iret) == 0) nchain = iret
             if(printable) write(nfout,'(" !** num_chain = ",i5)') nchain

             if (f_getIntValue(tag_sw_temperature_profile,iret) == 0) then
                sw_temperature_profile = iret
             endif

             if(sw_fix_bond==ON) call build_nbonds_per_thermo()

             if(sw_temperature_profile == OFF) then
               tdamp = 50.d0*dtio
               if(f_getRealValue(tag_tdamp,dret,'au_time')==0) tdamp = dret
               prealloc = .true.
               call set_thermostat(prealloc,nrsv,icounted)
               if(ipriinputfile >= 1 .and. printable) &
                    & write(nfout,'(" !** icounted = ",i6," at set_thermostat <<m_IS_rd_n>>")') icounted
               if(.not.number_is_given) then
                  if(printable) write(nfout,'(" !** number_is_given is .false.")')
                  nrsv = icounted
               else
                  if(printable) write(nfout,'(" !** number_is_given is .true.")')
               end if
               if(icounted < 0) nrsv = 1
               if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** nrsv = ",i6," <<m_IS_rd_n>>")') nrsv

               ! --- allocation and initialization of qmass, tkb, etc.
                 call T_control_alloc(nrsv)
               if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !** number of T-reserver = ",i6)') nrsv

               ! --- substitution for qmass(es) and temperature(s) of the reserver(s)
               prealloc = .false.
               call set_thermostat(prealloc,nrsv,icounted) ! -> qmass, tkb
               if(printable) write(nfout,'(" !** icounted = ",i6," <<set_thermostat>>")') icounted
               if(icounted < 0) then
                  call set_qmass_and_temp(1)
               end if
               if(ipriinputfile >= 1 .and. printable) then
                  do i = 1, nrsv
                     write(nfout,'(" !** reserver_id = ",i3," qmass = ",f20.6," tkb = ",f12.6)') i, qmass(i), tkb(i)
                  end do
               end if
! new way of defining thermostats
             else
               nrsv = get_nthermo(nfout)
               call alloc_temperature_profile(nrsv)
               call set_temperature_profile_no(nfout,nrsv,temperature_profile)
               call set_temperature_profile(nfout,nrsv,temperature_profile)
               if(sw_temperature_profile==ON) call set_new_temperature(nfout,nrsv,temperature_profile)
             endif

! ==== KT_add == 2014/06/10
             if ( f_getIntValue(tag_sw_change_temp_bystep,iret) == 0 ) then
                 sw_change_temperature_by_step = iret
                 write(nfout,*) '*** sw_change_temperature_by_step is ', iret
             endif
             if ( sw_change_temperature_by_step == ON ) then
                Do i=1, nrsv
                   write(nfout,'(A,i8,A,i8,A,F8.1,A,F8.1)') &
                        &           '!** MD Step (start): ', &
                        &              mdstep_at_start_thermostat(i), &
                        &           ', (end): ', &
                        &               mdstep_at_end_thermostat(i), &
                        &           ", Temp (start): ", &
                        &              temp_at_start_thermostat(i), &
                        &           ", (end): ", &
                        &              temp_at_end_thermostat(i)
                End Do
             endif
! ============== 2014/06/10

             if ( f_getIntValue(tag_iseed, iret)==0 )then
               iseed = iret
               if ( iseed>0) then
                 call random_seed(size=seedsize)
                 allocate(seeds(seedsize))
                 seeds = iseed
                 call random_seed(put=seeds)
                 deallocate(seeds)
               endif
             endif

             ! set initial velocity or not --> by J. Koga 2005
             if ( f_getIntValue(tag_set_initial_velocity,iret) == 0 ) then
                 set_initial_velocity = iret
             endif
             if ( f_getIntValue(tag_sw_read_velocities,iret) == 0 ) then
                 sw_read_velocities = iret
             endif
             if ( f_getIntValue(tag_sw_shift_velocities,iret) == 0 ) then
                 sw_shift_velocities = iret
             endif
             if( f_getRealValue(tag_initial_temperature,dret,"") == 0 ) then
                tk_initial = dret * CONST_kB
                if(ipriinputfile >= 1 .and. printable) then
                   write(nfout,'(" !** initial_temperature = ",f18.10)') dret
                   write(nfout,'(" !** tk_initial          = ",f18.10)') tk_initial
                end if
             end if

             if ( m_CtrlP_what_is_mdalg() == T_CONTROL .and. set_initial_velocity == ON ) then
                call set_initial_velocities(1)
             else if( m_CtrlP_what_is_mdalg() == VELOCITY_SCALING .and. set_initial_velocity == ON ) then
                call set_initial_velocities(1)
             else if( m_CtrlP_what_is_mdalg() == PT_CONTROL .and. set_initial_velocity == ON ) then
                call set_initial_velocities(1)
             else if( m_CtrlP_what_is_mdalg() == P_CONTROL .and. set_initial_velocity == ON ) then
                call set_initial_velocities(1)
             else if( m_CtrlP_what_is_mdalg() == VERLET .and. set_initial_velocity == ON ) then
                call set_initial_velocities(2)
             endif
             ! <--

             iret = f_selectParentBlock()
          end if
          ! --- Pressure control ---
          if( f_selectBlock( tag_pressure_control) == 0 )then
             if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* tag_pressure_control is found")')
             tdamp_baro = 500.d0*dtio
             if(f_getRealValue(tag_tdamp,dret,'au_time') == 0) tdamp_baro = dret
            ! dret = get_average_temperature()/300.d0
            ! if(dret.lt. 0.1d0) dret=0.1d0
             if(f_getStringValue(tag_method,rstr,LOWER)==0) then
                if(rstr == tag_volume)         p_ctrl_method = VOLUME
                if(rstr == tag_metric_tensor)  p_ctrl_method = METRIC_TENSOR
                if(rstr == tag_lattice_vector) p_ctrl_method = LATTICE_VECTOR
             endif
             !m_baro = 5e+12/(univol*tdamp_baro*tdamp_baro)
             !if(p_ctrl_method==VOLUME) m_baro = 2*m_baro
             !if(p_ctrl_method==LATTICE_VECTOR) m_baro = 10*m_baro*univol**(2.d0/3.d0)
             m_baro = get_default_mbaro(tdamp_baro)
             if(f_getRealValue(tag_m_baro,dret,"") == 0) m_baro = dret
             write(nfout,'(" !** mass of barostat ",f20.5)') m_baro
             if(f_getRealValue(tag_pressure,dret,"hartree/bohr3") == 0) target_pressure = dret
             if(f_getIntValue(tag_control_pressure,iret) == 0) control_pressure = iret
             if(f_getIntValue(tag_m11,iret) == 0 ) mobility_cell(1,1) = iret
             if(f_getIntValue(tag_m22,iret) == 0 ) mobility_cell(2,2) = iret
             if(f_getIntValue(tag_m33,iret) == 0 ) mobility_cell(3,3) = iret
             if(f_getIntValue(tag_m12,iret) == 0 ) then
                mobility_cell(1,2) = iret
                mobility_cell(2,1) = iret
             endif
             if(f_getIntValue(tag_m13,iret) == 0 ) then
                mobility_cell(1,3) = iret
                mobility_cell(3,1) = iret
             endif
             if(f_getIntValue(tag_m23,iret) == 0 ) then
                mobility_cell(2,3) = iret
                mobility_cell(3,2) = iret
             endif
             if(f_getIntValue(tag_sw_average_diagonal,iret) == 0) sw_average_diagonal = iret
             if(sw_average_diagonal == ON .and. printable)  write(nfout,'(a)') &
             &  ' !** diagonal elements of the stress tensor will be averaged'

             if (f_getIntValue(tag_sw_pressure_profile,iret) == 0) then
                sw_pressure_profile = iret
             endif
             if(sw_pressure_profile==ON) then
               call set_pressure_profile(nfout,pressure_profile)
               if((imdalg==P_CONTROL .or. imdalg==PT_CONTROL) .and. iteration_unit_cell==1) then
                 call set_new_pressure(nfout,pressure_profile)
               endif
             endif

             iret = f_selectParentBlock()
          endif
          iret = f_selectParentBlock()
       end if
    end if

  contains
    integer function constraint_is_possible(imdalg)
      integer, intent(in) :: imdalg
      character(len=10) :: word1
      select case (imdalg)
      case (CG_STROPT)
         constraint_is_possible = 1;         word1 = "CG"
      case (CG_STROPT2)
         constraint_is_possible = 1;         word1 = "CG2"
      case (T_CONTROL)
         constraint_is_possible = 1;         word1 = "T_CONTROL"
      case default
         constraint_is_possible = 0
      end select
      if(ipriinputfile >=1 .and. constraint_is_possible == 1) then
         write(nfout,'(" !** The constraint tag block is skipped &
           & ,because imdalg == ",a10)') word1
      else if(ipriinputfile >= 1) then
         write(nfout,'(" !** The constraint tag block is executed")')
      end if
      if(ipriinputfile >= 1) call flush(nfout)
    end function constraint_is_possible

    subroutine set_endpoint_images(rstr)
      character(*),intent(in) :: rstr
      call strncmp0(trim(rstr),tag_no,tf)
      if(tf) then
         endpoint_images = NO
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_nothing,tf)
      if(tf) then
         endpoint_images = NO
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_file,tf)
      if(tf) then
         endpoint_images = FILE
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_directin,tf)
      if(tf) then
         endpoint_images = DIRECTIN
         goto 1001
      end if
1001  continue
    end subroutine set_endpoint_images

    subroutine set_fixed_planes(prealloc,constraint_type, num,countedplanes,countedatoms)
      ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), 1st Aug. 2003
      logical, intent(in) ::  prealloc
      integer, intent(in) ::  constraint_type, num
      integer, intent(out) :: countedplanes
      integer, intent(out),optional :: countedatoms

      integer :: i, no, iret, j, icatm
      integer :: f_selectFirstTableLine,f_selectNextTableLine,f_getRealValue,f_getIntValue
      logical :: tf
      real(kind=dP) :: dret, normalvx,normalvy,normalvz,dlt

      character(len=FMAXUNITLEN) :: unit_f

      unit_f = 'bohr'

      icatm = 0
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
         if(prealloc) then
            if(f_getIntValue(tag_atom1,iret) == 0) then
               icatm = icatm + 1
               if(f_getIntValue(tag_atom2,iret) == 0) then
                  icatm = icatm + 1
                  if(f_getIntValue(tag_atom3,iret) == 0) then
                     icatm = icatm + 1
                     if(f_getIntValue(tag_atom4,iret) == 0) then
                        icatm = icatm + 1
                     end if
                  end if
               end if
            end if
         else
!!$            if(printable) write(nfout,'(" !** i = ",i5,"  num = ",i5)') i, num
            if(i > num) exit
            no = i
            if(f_getIntValue(tag_no, iret) == 0) then
               no = iret
            else if(f_getIntValue(tag_id, iret) == 0) then
               no = iret
            end if
            if(no <= num) then
               normalvx = 0.d0; normalvy = 0.d0; normalvz = 0.d0; dlt = 0.d0
               tf = f_getRealValue(tag_nx,dret,unit_f) == 0
               if(.not.tf) tf = f_getRealValue(tag_vx,dret,unit_f) == 0
               if(tf) normalvx = dret

               tf = f_getRealValue(tag_ny,dret,unit_f) == 0
               if(.not.tf) tf = f_getRealValue(tag_vy,dret,unit_f) == 0
               if(tf) normalvy = dret

               tf = f_getRealValue(tag_nz,dret,unit_f) == 0
               if(.not.tf) tf = f_getRealValue(tag_vz,dret,unit_f) == 0
               if(tf) normalvz = dret

               ! -- normalization of a normal vector --
               dret = normalvx**2 + normalvy**2 + normalvz**2
               if(dret < SmallestPositiveNumber*1.d5) then
                  normalvx = 1.d0; normalvy = 0.d0; normalvz = 0.d0
               else
                  dret = dsqrt(dret)
                  normalvx = normalvx/dret
                  normalvy = normalvy/dret
                  normalvz = normalvz/dret
               end if

               if( f_getRealValue(tag_delta,dret,unit_f) == 0)  dlt = dret

               max_reach_of_fixed_plane = c*0.5
               tf = f_getRealValue(tag_max_reach,dret,unit_f) == 0
               if(.not.tf) tf = f_getRealValue(tag_final_value,dret,unit_f) == 0
               if(tf) max_reach_of_fixed_plane = dret

               tf = f_getIntValue(tag_atom1,iret) == 0
               do j = 1, 4
                  if(tf) then
                     icatm = icatm + 1
                     ia_cnst(icatm) = iret
                     fcvect(icatm,1) = normalvx
                     fcvect(icatm,2) = normalvy
                     fcvect(icatm,3) = normalvz
                     fcvect(icatm,4) = dlt
                     ipfixedplane(icatm) = i
                     fcvect(icatm,5) = max_reach_of_fixed_plane
                     imdtyp(iret) = constraint_type
                     imdtypxyz(iret,:) = constraint_type
                     if(j == 1) tf = f_getIntValue(tag_atom2,iret) == 0
                     if(j == 2) tf = f_getIntValue(tag_atom3,iret) == 0
                     if(j == 3) tf = f_getIntValue(tag_atom4,iret) == 0
                     if(printable) write(nfout,'(" !** j = ",i5," normalvectors = ",3f8.4," dlt = ",f8.4)') &
                          & j, normalvx,normalvy,normalvz,dlt
                  else
                     exit
                  end if
               end do
            end if
         end if
         i = i+1
      end do
      countedplanes = i - 1
      if(prealloc) countedatoms = icatm
    end subroutine set_fixed_planes

    subroutine set_number_of_constraint_atoms(constraint_type,nfcatm)
      integer, intent(in) :: constraint_type
      integer, intent(out) :: nfcatm
      integer :: iret
      prealloc = .true.
      if( f_selectBlock(tag_atoms)==0) then
         if(ipriinputfile>=2) write(nfout,'(" !* tag_atoms is found")')
         Call set_constraint_atoms(prealloc, nfcatm, 0, constraint_type, iret) ! --> imdtyp,ia_cnst,nfcatm
         if(iret <= 0) stop ' constraint atoms are not given properly <<m_IS_rd_n>>'
         nfcatm = iret
         iret = f_selectParentBlock()
      else
         if(constraint_type == COG_FIX_L) then
            call phase_error_with_msg(nfout, &
            ' tag_atoms in <cog>, <cog_fix_l> or <cog_fix_in_a_plane> is not given <<m_IS_rd_n>>',__LINE__,__FILE__)
         else if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE .or. constraint_type == COG_and_RIGID_BODY_FIX_L) then
            call phase_error_with_msg(nfout, &
            ' tag_atoms in <rigid_body> or <rigid_body_fix_in_a_plane> is not given <<m_IS_rd_n>>',__LINE__,__FILE__)
         end if
      end if
    end subroutine set_number_of_constraint_atoms

    subroutine set_constraint_atoms_imdtyp_etc(mdalg,nfcatm,nfcatm_ip0,cnstrtype)
      integer, intent(in) :: mdalg,nfcatm, nfcatm_ip0,cnstrtype
      integer :: iret

      if(ipriinputfile>=2) write(nfout,'(" !* mdalg = ",i5, "<<set_constraint_atoms_imdtyp_etc>>")') mdalg
!!$      if(mdalg == GDIIS .or. mdalg == VERLET .or. mdalg == QUENCHED_MD &
!!$           & .or. mdalg == CG_STROPT .or. mdalg==SD_MD .or. mdalg==QUENCHED_CONSTRAINT .or. mdalg==CG_STROPT2 ) then

         prealloc = .false.
         if( f_selectBlock(tag_atoms)==0) then
            if(ipriinputfile>=1) write(nfout,'(" !* tag_atoms is found <<set_constraint_atoms_imdtyp_etc>>")')
            call set_constraint_atoms(prealloc, nfcatm, nfcatm_ip0, cnstrtype,iret) ! --> imdtyp, ia_cnst,nfcatm, ipfixedplane
            iret = f_selectParentBlock()
         else
            if(ipriinputfile>=1) write(nfout,'(" !* tag_atoms is not found <<set_constraint_atoms_imdtyp_etc>>")')
         end if
!!$      else
!!$         if(ipriinputfile>=DEBUGPRINTLEVEL) then
!!$            write(nfout,'(" mdalg = ",i8, " <<set_constraint_atoms_imdtyp_etc>>")')
!!$            write(nfout,'(" set_constraint_atoms_imdtyp_etc is not executed")')
!!$            call flush(nfout)
!!$         end if
!!$      end if
    end subroutine set_constraint_atoms_imdtyp_etc

    subroutine set_constraint_atoms(prealloc, nfcatm, nfcatm0,cnstrtype,iret)
      logical, intent(in)  :: prealloc
      integer, intent(in)  :: nfcatm, nfcatm0,cnstrtype
      integer, intent(out) :: iret
      integer :: i, rint, ip, len_rstrtrimmed, n_cp_name
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
      character(len=LEN_ATOMNAME) :: atomname
      character(len=FMAXVALLEN) :: rstr
      if(ipriinputfile >= 1) then
         write(nfout,'(" nfcatm = ",i5,", nfcatm0 = ",i5," <<set_constraint_atoms>>")') nfcatm, nfcatm0
         call flush(nfout)
      end if

      i = 1
      do while(.true.)
         if( i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
!!$         write(nfout,'(" i = ",i5, " nfcatm = ",i5)') i, nfcatm
!!$         call flush(nfout)
         if(.not.prealloc) then
            if(i > nfcatm) exit
            ip = i
            tf = f_getIntValue(tag_number,rint) == 0
            if(.not.tf) tf = f_getIntValue(tag_no,rint) == 0
            if(tf) then
               if(ipriinputfile >= DEBUGPRINTLEVEL+1) write(nfout,'(" (no.) rint = ",i5)') rint
               if(0<rint .and. rint <=natm) then
                  if(imdtyp_set(rint)==YES) then
                     if(ipriinputfile >= 1) &
                     & write(nfout,'(" constraint atoms number ( = ",i5," ) is doubly signated. <<set_constraint_atoms>>")') rint
                     call phase_error_with_msg(nfout,' constraint atoms number is doubly signated. <<set_constraint_atoms>>' &
                                             , __LINE__, __FILE__)
                  else if(imdtyp(rint) == FIX) then
                     if(ipriinputfile >=1 ) &
                          & write(nfout,'(" imdtyp(",i5,") (= ",i5,") is set FIX previously ")') rint,imdtyp(rint)
                     call phase_error_with_msg(nfout, ' imdtyp(rint) is set FIX previously <<set_constratin_atoms>>' &
                                             , __LINE__, __FILE__)
                  else
                     ia_cnst(nfcatm0+ip) = rint
                  end if
               else
                  if(ipriinputfile >=1)  write(nfout,'(" given atom number is illegal. no = ",i5," <<set_constraint_atoms>>")') &
                                         rint
                  call phase_error_with_msg(nfout, ' Given atom number is illegal. <<set_constraint_atoms>>' &
                                          , __LINE__, __FILE__)
               end if
               if(ipriinputfile>=DEBUGPRINTLEVEL+1) &
               & write(nfout,'(" imdtyp set previouslly = ",i5, ", cnstrtype given = ",i5)') &
               & imdtyp(ia_cnst(nfcatm0+ip)), cnstrtype
               imdtyp_set(rint) = YES
               imdtyp(ia_cnst(nfcatm0+ip)) = cnstrtype
               imdtypxyz(ia_cnst(nfcatm0+ip),1:3) = cnstrtype
            end if

            if( f_getStringValue(tag_element, rstr, NOCONV) == 0) then
               len_rstrtrimmed = len_trim(rstr)
               if(len_rstrtrimmed > LEN_ATOMNAME) then
                  if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* element name is larger than ",i6)') LEN_ATOMNAME
                  n_cp_name = LEN_ATOMNAME
               else
                  n_cp_name = len_rstrtrimmed
               end if
               atomname(1:n_cp_name) = rstr(1:n_cp_name)
               if(trim(species_work(ia_cnst(nfcatm0+ip))) /= trim(rstr(1:n_cp_name))) then
                  if(ipriinputfile >= 1) then
                     write(nfout,'(" !* atomname is not properly given in <atoms> list <<set_constraint_atoms>>")')
                     write(nfout,'(" !* atomname given is ",a4," atom number = ",i5)') rstr(1:n_cp_name), rint
                     call flush(nfout)
                  end if
                  call phase_error_with_msg(nfout, '!* atomname is not properly given in <atoms> list <<set_constraint_atoms>>' &
                                            , __LINE__, __FILE__)
               end if
            end if

            tf = f_getIntValue(tag_iplane,rint)==0
            if(.not.tf) tf = f_getIntValue(tag_id,rint) == 0
            if(tf) then
               if(ipriinputfile >= DEBUGPRINTLEVEL+1) write(nfout,'(" (id) rint = ",i5)') rint
               if(0<rint) then
                  ipfixedplane(nfcatm0+ip) = rint
               else
                  if(ipriinputfile >=1)  &
                  write(nfout,'(" given iplane number is illegal. iplane = ",i5," <<set_constraint_atoms>>")') rint
                  call phase_error_with_msg(nfout,' Given iplane number is illegal. <<set_constraint_atoms>>'  &
                                            , __LINE__, __FILE__)
               end if
            else
               ipfixedplane(nfcatm0+ip) = 1
            end if
            if(ipriinputfile >= DEBUGPRINTLEVEL+1) write(nfout,'(" ipfixedplane(",i5,") = ",i5)') &
            nfcatm0+ip, ipfixedplane(nfcatm0+ip)

         end if
         i = i + 1
         if(i > 1000) exit
      end do
      iret = i - 1
      if(.not.prealloc) then
         if(ipriinputfile >= DEBUGPRINTLEVEL) then
            do i = 1, nfcatm
               write(nfout,'(" !constraint_atoms number = ",i5, " constraint_type = ",i5)') ia_cnst(nfcatm0+i), cnstrtype
            end do
         end if
      end if
    end subroutine set_constraint_atoms

    subroutine set_pointer_array_to_constrained_atoms(nfcatm)
      integer, intent(in) :: nfcatm
      integer :: ia, ifchit, ifc

      if(.not.allocated(icnst_a)) allocate(icnst_a(natm))
      icnst_a = 0
      do ia = 1, natm
         if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L &
              .or. imdtyp(ia) == FIX_IN_A_PLANE .or. imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
            ifchit = 0
            do ifc = 1, nfcatm
               if(ia_cnst(ifc) == ia) then
                  ifchit = ifc
                  exit
               end if
            end do
            if(ifchit>=1) icnst_a(ia) = ifchit
         end if
      end do
    end subroutine set_pointer_array_to_constrained_atoms

    subroutine check_number_of_planes_atoms_are_fixed(nfcatm_rb,nfcatm_cog,num_planes_all,cnstrtype)
      integer, intent(in)  :: nfcatm_rb,nfcatm_cog,num_planes_all,cnstrtype
      integer, allocatable, dimension(:) :: iplane
      integer :: i, j, ia, icount, ip, ir, nfcatm_t, nfcatm0
      logical :: hit


      if(cnstrtype /= COG_FIX_L .and. cnstrtype /= RIGID_BODY_FIX_IN_A_PLANE .and. cnstrtype /= COG_and_RIGID_BODY_FIX_L) then
         num_planes_atoms_are_fixed_cog = 0; num_planes_atoms_are_fixed_rb = 0
         goto 1005
      end if

      if(nfcatm_cog+nfcatm_rb > 0) then
         allocate(iplane(max(nfcatm_cog, nfcatm_rb))); iplane = 0
      else
         num_planes_atoms_are_fixed_cog = 0; num_planes_atoms_are_fixed_rb = 0
         goto 1005
      end if
      nfcatm0 = 0
      do ir = 1, 2
         if((ir==1 .and. nfcatm_cog >=1) .or. (ir==2 .and. nfcatm_rb >=1) ) then
            if(ir==1) num_planes_atoms_are_fixed_cog = 0
            if(ir==2) num_planes_atoms_are_fixed_rb = 0
            icount = 0
            if(ir==1) nfcatm_t = nfcatm_cog
            if(ir==2) nfcatm_t = nfcatm_rb
            do i = 1, nfcatm_t
               ia = ia_cnst(nfcatm0+i)
               if((ir==1.and.imdtyp(ia) == COG_FIX_L).or.(ir==2.and.imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE)) then
                  ip = ipfixedplane(nfcatm0+i)
                  if(icount == 0) then
                     iplane(icount+1) = ip
                     icount = icount + 1
                     icount_of_ipfixedplane(nfcatm0+i) = icount
                  else
                     hit = .false.
                     do j = 1, icount
                        if(ip == iplane(j)) then
                           hit = .true.
                           icount_of_ipfixedplane(nfcatm0+i) = j
                           exit
                        end if
                     end do
                     if(.not.hit) then
                        iplane(icount+1) = ip
                        icount = icount + 1
                        icount_of_ipfixedplane(nfcatm0+i) = icount
                     end if
                  end if
               end if
            end do
            if(ir==1) num_planes_atoms_are_fixed_cog = icount
            if(ir==2) num_planes_atoms_are_fixed_rb = icount
         end if
         nfcatm0 = nfcatm_cog
      end do

      icount = 0
      do i = 1, nfcatm
         ia = ia_cnst(i)
         if(imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
            ip = ipfixedplane(i)
            if(icount==0) then
               iplane(i) = ip
               icount = icount+1
            else
               hit = .false.
               do j = 1, icount
                  if(ip==iplane(j)) then
                     hit = .true.
                     exit
                  end if
               end do
               if(.not.hit) then
                  iplane(icount+1) = ip
                  icount = icount+1
               end if
            end if
         end if
      end do
      if(icount < num_planes_atoms_are_fixed) then
         if(ipriinputfile >=1) write(nfout, &
         '(" num_planes_atoms_are_fixed is reduced to ",i5," from ",i5)') icount, num_planes_atoms_are_fixed
         num_planes_atoms_are_fixed = icount
      end if

1005  continue
      if(allocated(iplane)) deallocate(iplane)
      if(ipriinputfile>=DEBUGPRINTLEVEL) then
         if(nfcatm_cog+nfcatm_rb>0) then
            do i =1, nfcatm_cog+nfcatm_rb
               write(nfout,'(" icount_of_ipfixedplane(",i5," ) = ",i5)') i, icount_of_ipfixedplane(i)
            end do
         end if
      end if
      if(nfcatm_cog>=1) then
         if(ipriinputfile>=DEBUGPRINTLEVEL) write(nfout,'("num_planes_atoms_are_fixed_cog = ",i5 &
              &                  ," <<check_number_of_planes_atoms_are_fixed>>")') num_planes_atoms_are_fixed_cog
         if(num_planes_atoms_are_fixed_cog <=0) then
            if(ipriinputfile>=1) write(nfout,'(" num_planes_atoms_are_fixed_cog = ",i5," <=0, for nfcatm_cog = ",i5 &
                 &               , " <<check_number_of_planes_atoms_are_fixed>>")') num_planes_atoms_are_fixed_cog, nfcatm_cog
            call phase_error_with_msg(nfout &
            , ' num_planes_atoms_are_fixed_cog is smaller than or equal zero <<check_number_of_planes_atoms_are_fixed>>' &
            , __LINE__, __FILE__)
         end if
      end if
      if(nfcatm_rb>=1) then
         if(ipriinputfile>=DEBUGPRINTLEVEL) write(nfout,'("num_planes_atoms_are_fixed_rb  = ",i5 &
              &                  ," <<check_number_of_planes_atoms_are_fixed>>")') num_planes_atoms_are_fixed_rb
         if(num_planes_atoms_are_fixed_rb <=0) then
            if(ipriinputfile>=1) write(nfout,'(" num_planes_atoms_are_fixed_rb = ",i5," <=0, for nfcatm_rb = ",i5 &
                 &               , " <<check_number_of_planes_atoms_are_fixed>>")') num_planes_atoms_are_fixed_rb, nfcatm_rb
            call phase_error_with_msg(nfout &
            , ' num_planes_atoms_are_fixed_rb is smaller than or equal zero <<check_number_of_planes_atoms_are_fixed>>' &
            , __LINE__, __FILE__)
         end if
      end if
      if(num_planes_all < min(num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb)) then
         if(ipriinputfile>=1) then
            write(nfout,'(" num_planes_atoms_are_fixed = ",i5," <<check_number_of_planes_atoms_are_fixed>>")') num_planes_all
            write(nfout,'(" num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb = ",2i5)') &
                 &                                         num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb
            write(nfout, &
            '(" num_planes_atoms_are_fixed is smaller than num_planes_atoms_are_fixed_cog and num_planes_atoms_are_fixed_rb")')
         end if
         call phase_error_with_msg(nfout &
         , ' num_planes_atoms_are_fixed is smaller than num_planes_atoms_are_fixed_cog and num_planes_atoms_are_fixed_rb' &
         , __LINE__, __FILE__)
      end if
      if(num_planes_atoms_are_fixed_cog>=1) then
         allocate(distance_cog(num_planes_atoms_are_fixed_cog)); distance_cog = 0.d0
      end if
      if(num_planes_atoms_are_fixed_rb>=1) then
         allocate(distance_rb(num_planes_atoms_are_fixed_rb)); distance_rb = 0.d0
      end if
    end subroutine check_number_of_planes_atoms_are_fixed

    subroutine set_constraint_plane(prealloc,nfcatm,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: nfcatm
      integer, intent(out):: iret
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
      integer :: i, j, iplane, rint, sw_relax_or_not
      logical :: normalvx_is_given, normalvy_is_given, normalvz_is_given, Is_table_format
      logical :: incrvx_is_given, incrvy_is_given, incrvz_is_given
      real(kind=dP) :: dret, normalvx,normalvy,normalvz,dlt, incrvx, incrvy, incrvz

      character(len=FMAXUNITLEN) :: unit_f

      unit_f = 'bohr'

      i = 1
      do while(.true.)
         if( i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         if(.not.prealloc) then
            tf = f_getIntValue(tag_iplane,rint) == 0
            if(.not.tf) tf = f_getIntValue(tag_id,rint) == 0
            if(tf) then
               iplane = rint
            else
               iplane = 1
            end if

            normalvx = 0.d0; normalvy = 0.d0; normalvz = 0.d0; dlt = 0.d0
            tf = f_getRealValue(tag_nx,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vx,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_normx,dret,"") == 0
            if(tf) normalvx = dret
            tf = f_getRealValue(tag_ny,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vy,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_normy,dret,"") == 0
            if(tf) normalvy = dret
            tf = f_getRealValue(tag_nz,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vz,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_normz,dret,"") == 0
            if(tf) normalvz = dret
            if(ipriinputfile >= DEBUGPRINTLEVEL) &
                 & write(nfout,'(" normalvx,normalvy,normalvz = ",3f8.4, " iplane = ",i8)') normalvx,normalvy,normalvz,iplane

            ! -- normalization of a normal vector --
            dret = normalvx**2 + normalvy**2 + normalvz**2
            if(dret < SmallestPositiveNumber*1.d5) then
               normalvx = 1.d0; normalvy = 0.d0; normalvz = 0.d0
            else
               dret = dsqrt(dret)
               normalvx = normalvx/dret
               normalvy = normalvy/dret
               normalvz = normalvz/dret
            end if
            dlt = dret

            incrvx_is_given = .false.;    incrvy_is_given = .false.;    incrvz_is_given = .false.
            if(f_getRealValue(tag_incrvx,dret,unit_f) == 0) then
               incrvx = dret; incrvx_is_given = .true.
            end if
            if(f_getRealValue(tag_incrvy,dret,unit_f) == 0) then
               incrvy = dret; incrvy_is_given = .true.
            end if
            if(f_getRealValue(tag_incrvz,dret,unit_f) == 0) then
               incrvz = dret; incrvz_is_given = .true.
            end if
            if(.not.incrvx_is_given) incrvx = normalvx
            if(.not.incrvy_is_given) incrvy = normalvy
            if(.not.incrvz_is_given) incrvz = normalvz
            ! -- normalization of a increment vector --
            dret = incrvx**2 + incrvy**2 + incrvz**2
            if(dret < SmallestPositiveNumber*1.d5) then
               incrvx = 1.d0; incrvy = 0.d0; incrvz = 0.d0
            else
               dret = dsqrt(dret)
               incrvx = incrvx/dret
               incrvy = incrvy/dret
               incrvz = incrvz/dret
            end if

!!$            dlt = dret
            tf = f_getRealValue(tag_delta,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_increment,dret,unit_f) == 0
            if(tf) dlt = dret

            max_reach_of_fixed_plane = c*0.5
            tf = f_getRealValue(tag_max_reach,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_final_value,dret,unit_f) == 0
            if(tf) max_reach_of_fixed_plane = dret

            ! -- relax_in_plane --
            sw_relax_or_not = NO
            tf = f_getIntValue(tag_relax_in_fixedplane,iret) == 0
            if(tf) sw_relax_or_not = iret

            do j = 1, nfcatm
               if(ipfixedplane(j)==iplane) then
                  fcvect(j,1) = normalvx
                  fcvect(j,2) = normalvy
                  fcvect(j,3) = normalvz
                  fcvect(j,4) = dlt
                  fcvect(j,5) = max_reach_of_fixed_plane
                  fcvect(j,6) = incrvx
                  fcvect(j,7) = incrvy
                  fcvect(j,8) = incrvz
                  relax_in_fixedplane(j) = sw_relax_or_not
               end if
            end do

         end if
         i = i + 1
      end do
      iret = i - 1
      if(iret >= 1)  Is_table_format = .true.

      if(iret == 0) then
         normalvx_is_given = .false.; normalvy_is_given = .false.; normalvz_is_given = .false.
         normalvx = 0.d0; normalvy = 0.d0; normalvz = 0.d0; dlt = 0.d0
         tf = f_getRealValue(tag_nx,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_vx,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_normx,dret,"") == 0
         if(tf) normalvx_is_given = .true.
         if(.not.prealloc .and. tf) normalvx = dret
         tf = f_getRealValue(tag_ny,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_vy,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_normy,dret,"") == 0
         if(tf) normalvy_is_given = .true.
         if(.not.prealloc .and. tf) normalvy = dret
         tf = f_getRealValue(tag_nz,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_vz,dret,unit_f) == 0
         if(.not.tf) tf = f_getRealValue(tag_normz,dret,"") == 0
         if(tf) normalvz_is_given = .true.
         if(.not. prealloc .and. tf) normalvz = dret

         if(normalvx_is_given .and. normalvy_is_given .and. normalvz_is_given) then
            iret = 1
         else
            iret = 0
         end if
         if(.not.prealloc) then
            if(ipriinputfile >= 1) write(nfout,'(" normalvx,normalvy,normalvz = ",3f8.4)') normalvx,normalvy,normalvz


            ! -- normalization of a normal vector --
            dret = normalvx**2 + normalvy**2 + normalvz**2
            if(dret < SmallestPositiveNumber*1.d5) then
               normalvx = 1.d0; normalvy = 0.d0; normalvz = 0.d0
            else
               dret = dsqrt(dret)
               normalvx = normalvx/dret
               normalvy = normalvy/dret
               normalvz = normalvz/dret
            end if

            dlt = dret

            incrvx_is_given = .false.;    incrvy_is_given = .false.;    incrvz_is_given = .false.
            if(f_getRealValue(tag_incrvx,dret,unit_f)==0) then
               incrvx_is_given = .true. ; incrvx = dret
            end if
            if(f_getRealValue(tag_incrvy,dret,unit_f)==0) then
               incrvy_is_given = .true. ; incrvy = dret
            end if
            if(f_getRealValue(tag_incrvz,dret,unit_f)==0) then
               incrvz_is_given = .true. ; incrvz = dret
            end if
            if(.not.incrvx_is_given) incrvx = normalvx
            if(.not.incrvy_is_given) incrvy = normalvy
            if(.not.incrvz_is_given) incrvz = normalvz

            tf = f_getRealValue(tag_delta,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_increment,dret,unit_f) == 0
            if(tf) dlt = dret

            max_reach_of_fixed_plane = c*0.5
            tf = f_getRealValue(tag_max_reach,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_final_value,dret,unit_f) == 0
            if(tf) max_reach_of_fixed_plane = dret

            ! -- relax_in_plane --
            sw_relax_or_not = NO
            tf = f_getIntValue(tag_relax_in_fixedplane,iret) == 0
            if(tf) sw_relax_or_not = iret

            do j = 1, nfcatm
               fcvect(j,1) = normalvx
               fcvect(j,2) = normalvy
               fcvect(j,3) = normalvz
               fcvect(j,4) = dlt
               fcvect(j,5) = max_reach_of_fixed_plane
               if(ipfixedplane(j)<=0 .or. 1 <ipfixedplane(j)) then
                  if(1 < ipfixedplane(j)) then
                     if(ipriinputfile >=1 ) write(nfout,'(" !* number of fixed_plane is not enough. Ipfixedplane(",i5,") = ",i3)') &
                          &                   j, ipfixedplane(j)
                     call phase_error_with_msg(nfout, ' number of fixed_plane is short <<set_constraint_plane>>' &
                                             , __LINE__, __FILE__)
                  end if
               end if
               fcvect(j,6) = incrvx
               fcvect(j,7) = incrvy
               fcvect(j,8) = incrvz
               ipfixedplane(j) = 1
               relax_in_fixedplane(j) = sw_relax_or_not
            end do
         end if
      end if

      if(iret == 0) then
         if(ipriinputfile >=1 ) write(nfout,'(" !* fixed_plane is not given properly <<set_constraint_plane>>")')
         call phase_error_with_msg(nfout, ' fixed_plane is not given properly <<set_constraint_plane>>', __LINE__, __FILE__)
      end if

      if(Is_table_format) then
         if(ipriinputfile>=2) then
            do j = 1, nfcatm
               write(nfout,'(" !* ipfixedplane(",i5,") = ",i3," <<set_constraint_plane>>")') j, ipfixedplane(j)
            end do
         end if
         do j = 1, nfcatm
            if(ipfixedplane(j)<=0 .and. iret <ipfixedplane(j)) then
               if(ipriinputfile >=1 ) write(nfout,'(" !* ipfixedplane(",i5,") = ",i3, " is smaller than 0 or larger than ",i5)') &
                    & j, ipfixedplane(j), iret
               call phase_error_with_msg(nfout, ' number of fixed_plane is not given properly',__LINE__,__FILE__)
            end if
         end do
      else
         if(iret /= 1) then
            if(ipriinputfile >=1 ) write(nfout,'(" !* number_of_fixed_plane is not 1, iret = ",i5)') iret
            call phase_error_with_msg(nfout, ' !* number_of_fixed_plane is not 1',__LINE__,__FILE__)
         end if
      end if

    end subroutine set_constraint_plane

    subroutine set_rigid_body_fixed_plane(constraint_type)
      ! This subroutine was coded by T. Yamasaki, 27th Apr. 2021
      integer, intent(in) ::  constraint_type

      integer :: i, no, iret, j
      integer :: f_selectFirstTableLine,f_selectNextTableLine,f_getRealValue,f_getIntValue
      logical :: tf
      real(kind=dP) :: dret, normalvx,normalvy,normalvz,dlt

      character(len=FMAXUNITLEN) :: unit_f

      unit_f = 'bohr'

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
         no = i
         if(f_getIntValue(tag_no, iret) == 0) then
            no = iret
         else if(f_getIntValue(tag_id, iret) == 0) then
            no = iret
         end if
         if(no == 1) then
            normalvx = 0.d0; normalvy = 0.d0; normalvz = 0.d0; dlt = 0.d0
            tf = f_getRealValue(tag_nx,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vx,dret,unit_f) == 0
            if(tf) normalvx = dret
            tf = f_getRealValue(tag_ny,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vy,dret,unit_f) == 0
            if(tf) normalvy = dret
            tf = f_getRealValue(tag_nz,dret,unit_f) == 0
            if(.not.tf) tf = f_getRealValue(tag_vz,dret,unit_f) == 0
            if(tf) normalvz = dret

            ! -- normalization of a normal vector --
            dret = normalvx**2 + normalvy**2 + normalvz**2
            if(dret < SmallestPositiveNumber*1.d5) then
               normalvx = 1.d0; normalvy = 0.d0; normalvz = 0.d0
            else
               dret = dsqrt(dret)
               normalvx = normalvx/dret
               normalvy = normalvy/dret
               normalvz = normalvz/dret
            end if

            if( f_getRealValue(tag_delta,dret,unit_f) == 0)  dlt = dret

            rigid_body_vect(1) = normalvx
            rigid_body_vect(2) = normalvy
            rigid_body_vect(3) = normalvz
            rigid_body_vect(4) = dlt
            exit
         end if
      end do
    end subroutine set_rigid_body_fixed_plane

    subroutine set_cog_fix(countedatoms)
      ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), 2nd Aug. 2003
      integer, intent(out) :: countedatoms

      integer :: i, no, iret, j
      integer :: f_selectFirstTableLine,f_selectNextTableLine,f_getIntValue
      logical :: tf

!!$      character(len=FMAXUNITLEN) :: unit_f
!!$      unit_f = 'bohr'

      i = 0
      countedatoms = 0
      if( f_selectBlock(tag_thermostat) == 0) then
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
!!$            if(i > num) exit
            no = i
            if(f_getIntValue(tag_no, iret) == 0) then
               no = iret
            else if(f_getIntValue(tag_id, iret) == 0) then
               no = iret
            end if
            if(no <= natm) then
               tf = f_getIntValue(tag_atom1,iret) == 0
               do j = 1, 4
                  if(tf) then
                     countedatoms = countedatoms + 1
!!$                     ia_cnst(countedatoms) = iret
                     imdtyp(iret) = COG_FIX
                     imdtypxyz(iret,:) = COG_FIX
                     if(j == 1) tf = f_getIntValue(tag_atom2,iret) == 0
                     if(j == 2) tf = f_getIntValue(tag_atom3,iret) == 0
                     if(j == 3) tf = f_getIntValue(tag_atom4,iret) == 0
                  else
                     exit
                  end if
               end do
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      end if
    end subroutine set_cog_fix

    function get_nthermo(nfout) result(res)
      integer, intent(in) :: nfout
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: i,j,ithermo,iret,iiret
      integer, allocatable, dimension(:) :: thermo_ids
      integer :: max_nthermo
      logical :: found
      integer :: res
      max_nthermo = natm
      allocate(thermo_ids(max_nthermo));thermo_ids = 0
      if( f_selectBlock(tag_thermostat) == 0) then
         i = 1
         ithermo = 0
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
            i = i+1
            iiret = f_getIntValue(tag_no,iret)
            if(iiret /= 0) iiret = f_getIntValue(tag_id,iret)
            if(iiret /= 0) iret = 1
            found = .false.
            do j=1,max_nthermo
                if(thermo_ids(j) == iret)then
                call flush(0)
                   found = .true.
                   exit
                endif
            enddo
            if(.not. found) then
                if(ithermo>max_nthermo)then
                   write(nfout,'(a,i5)') ' !** exceeded maximum number of thermostats : ',max_nthermo
                   exit
                endif
                ithermo = ithermo+1
                thermo_ids(ithermo) = iret
            endif
         end do
         iiret = f_selectParentBlock()
      end if
      res = ithermo
      deallocate(thermo_ids)
    end function get_nthermo

    subroutine set_temperature_profile_no(nfout,nrsv,temperature_profile)
      integer, intent(in) :: nfout,nrsv
      type(temperature_profile_t), intent(inout), dimension(nrsv) :: temperature_profile
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: i,j,ithermo,iret,iiret
      logical :: found
      integer, parameter :: max_nthermo = 100
      if( f_selectBlock(tag_thermostat) == 0) then
         i = 1
         ithermo = 0
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
            i = i+1
            iiret = f_getIntValue(tag_no,iret)
            if(iiret /= 0) iiret = f_getIntValue(tag_id,iret)
            if(iiret /= 0) iiret = f_getIntValue(tag_thermo_group,iret)
            if(iiret /= 0) iret = 1
            found = .false.
            do j=1,nrsv
                if(temperature_profile(j)%no == iret)then
                   found = .true.
                   exit
                endif
            enddo
            if(.not. found) then
               ithermo = ithermo+1
               temperature_profile(ithermo)%no = iret
            endif
         end do
         iiret = f_selectParentBlock()
      end if
    end subroutine set_temperature_profile_no

    subroutine set_temperature_profile(nfout,nrsv,temperature_profile)
      integer, intent(in) :: nfout,nrsv
      type(temperature_profile_t), intent(inout), dimension(nrsv) :: temperature_profile
      integer :: i,j,k,iiret,iret
      integer :: nprof,currprof
      real(kind=DP) :: dret
      integer :: f_selectFirstTableLine,f_selectNextTableLine

      ! first count the number of profiles for each thermostat
      if( f_selectBlock(tag_thermostat) == 0) then
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
            i = i+1
            iiret = f_getIntValue(tag_no,iret)
            if(iiret /= 0) iiret = f_getIntValue(tag_id,iret)
            if(iiret /= 0) iret = 1
            do j=1,nrsv
                if(temperature_profile(j)%no == iret)then
                   temperature_profile(j)%nprof = temperature_profile(j)%nprof + 1
                endif
            enddo
         end do
         iiret = f_selectParentBlock()
      end if

      ! allocate
      do j=1,nrsv
         nprof = temperature_profile(j)%nprof
         allocate(temperature_profile(j)%tempi(nprof));temperature_profile(j)%tempi = -1.d0
         allocate(temperature_profile(j)%tempf(nprof));temperature_profile(j)%tempf = -1.d0
         allocate(temperature_profile(j)%till_n(nprof));temperature_profile(j)%till_n = -1
         allocate(temperature_profile(j)%qmass(nprof));temperature_profile(j)%qmass = -1.d0
         allocate(temperature_profile(j)%tdamp(nprof));temperature_profile(j)%tdamp = 50.d0*dtio
         temperature_profile(j)%currprof = 0
      enddo

      if( f_selectBlock(tag_thermostat) == 0) then
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
            i = i+1
            iiret = f_getIntValue(tag_no,iret)
            if(iiret /= 0) iiret = f_getIntValue(tag_id,iret)
            if(iiret /= 0) iret = 1
            do j=1,nrsv
               if(temperature_profile(j)%no == iret)then
                  temperature_profile(j)%currprof = temperature_profile(j)%currprof+1
                  currprof = temperature_profile(j)%currprof
                  if(f_getRealValue(tag_temp,dret,'K') == 0) then
                     temperature_profile(j)%tempi(currprof) = dret * CONST_kB
                     temperature_profile(j)%tempf(currprof) = dret * CONST_kB
                  endif
                  if(f_getRealValue(tag_tempi,dret,'K') == 0) then
                      temperature_profile(j)%tempi(currprof) = dret * CONST_kB
                      temperature_profile(j)%tempf(currprof) = dret * CONST_kB
                  endif
                  if(f_getRealValue(tag_tempf,dret,'K') == 0) &
                  &  temperature_profile(j)%tempf(currprof) = dret * CONST_kB
                  if(f_getIntValue(tag_till_n,iret) == 0) &
                  &  temperature_profile(j)%till_n(currprof) = iret
                  if(f_getRealValue(tag_tdamp,dret,'au_time') == 0)then
                      temperature_profile(j)%tdamp(currprof) = dret
                  endif
                  if(f_getRealValue(tag_qmass,dret,'au_mass') == 0)then
                     temperature_profile(j)%qmass(currprof) = dret
                  endif
               endif
            enddo
         end do
         iiret = f_selectParentBlock()
         if(printable)then
             write(nfout,'(a,i4)')    ' !** number of thermostats ',nrsv
             do j=1,nrsv
                write(nfout,'(a,i4)') ' !** status of thermostat no. ',j
                write(nfout,'(a,i4)') ' !** number of profiles : ',temperature_profile(j)%nprof
                write(nfout,'(a,i8)') ' !** number of atoms    : ',temperature_profile(j)%natm
                do k=1,temperature_profile(j)%nprof
                   write(nfout,'(a,i4)') ' !** profile no. ',k
                   write(nfout,'(a,2f20.3,i8)') ' !** tempi, tempf, till_n : '    &
                   &                           , temperature_profile(j)%tempi(k)/CONST_kB &
                   &                           , temperature_profile(j)%tempf(k)/CONST_kB &
                   &                           , temperature_profile(j)%till_n(k)
                   write(nfout,'(a,2f20.3)') ' !** qmass, tdamp: ', temperature_profile(j)%qmass(k) &
                   &                           , temperature_profile(j)%tdamp(k)
                   if(temperature_profile(j)%tempi(k) <= 0)then
                      call phase_error_with_msg(nfout,'initial temperature must be positive ',__LINE__,__FILE__)
                   endif
                   if(temperature_profile(j)%tempf(k) <= 0)then
                      call phase_error_with_msg(nfout,'final temperature must be positive ',__LINE__,__FILE__)
                   endif
                enddo
             enddo
         endif
      end if
    end subroutine set_temperature_profile

    subroutine set_thermostat(prealloc, num, icounted)
      ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), 1st Aug. 2003
      logical, intent(in) ::  prealloc
      integer, intent(in) ::  num
      integer, intent(out) :: icounted

      integer :: i, no, iret

      integer :: f_selectFirstTableLine, f_selectNextTableLine
      i = 0
      if( f_selectBlock(tag_thermostat) == 0) then
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
               if(i > num) exit
               no = i
               if(f_getIntValue(tag_no, iret) == 0) then
                  no = iret
               else if(f_getIntValue(tag_id, iret) == 0) then
                  no = iret
               end if
               if(no <= num) then
                  call set_qmass_and_temp(no)
               end if
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      end if
      icounted = i - 1
    end subroutine set_thermostat

    subroutine set_qmass_and_temp(no)
      ! This subroutine was coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), 29th Jul. 2003
      integer, intent(in) :: no

      real(kind=DP) :: dret
      integer :: f_getRealValue

      character(len=FMAXUNITLEN) :: unit_f
      real(kind=DP) :: tdam
      integer :: i,ir,j,nfree
      logical :: temp_is_set

      ! --- temperature ---
      temp_is_set = .false.

      unit_f = 'K'
      number_is_given = f_getRealValue(tag_temperature,dret,unit_f) == 0
      if(.not.number_is_given) number_is_given = f_getRealValue(tag_temp,dret,unit_f) == 0
      if(.not.number_is_given) number_is_given = f_getRealValue(tag_T,dret,unit_f) == 0
      if(number_is_given) then
         tkb(no) = dret * CONST_kB
         temp_is_set = .true.
      endif

      ! --- qmass ---
      unit_f = 'au_mass'
      number_is_given = f_getRealValue(tag_weight_thermostat,dret,unit_f) == 0
      if(.not.number_is_given) number_is_given = f_getRealValue(tag_weight_thermo,dret,unit_f) == 0
      if(.not.number_is_given) number_is_given = f_getRealValue(tag_weight,dret,unit_f) == 0
      if(.not.number_is_given) number_is_given = f_getRealValue(tag_qmass,dret,unit_f) == 0
      if(number_is_given) then
         qmass(no) = dret
      else
         tdam = tdamp
         if(f_getRealValue(tag_tdamp,dret,'au_time')==0) tdam = dret
         nfree = 3.d0*natm_per_thermo(no)
         if(sw_fix_bond==ON) nfree = nfree-nbonds_per_thermo(no)
         qmass(no) = (2.d0*nfree*(tdam/PAI2)**2) * tkb(no)
         if(printable) write(nfout,'(a,i3,a,f20.5)') ' !** resolved Qmass for thermostat ',no,' : ',qmass(no)
      endif
      if(nchain>1)then
         qmass_c(1:nrsv,1) = qmass(1:nrsv)
         do i=2,nchain
            qmass_c(1:nrsv,i) = qfactor * qmass(1:nrsv)
         enddo
      endif
      if(printable .and. iprimd>=2) then
        if(nchain==1) then
          do i=1,nrsv
            write(nfout,'(a,i5,a,f20.5)') ' !** qmass for thermostat ',i,':',qmass(i)
          enddo
        else
          do j=1,nchain
            do i=1,nrsv
              write(nfout,'(a,i5,a,i5,a,f20.5)') ' !** qmass for chain ',j,' thermostat ',i,':',qmass_c(i,j)
            enddo
          enddo
        endif
      endif

! ==== KT_add === 2014/06/10
      if (f_getIntValue(tag_duration,iret)==0) then
         if ( no == 1 ) then
            mdstep_at_start_thermostat(no) = 1
            mdstep_at_end_thermostat(no) = iret
         else
            mdstep_at_start_thermostat(no) = mdstep_at_end_thermostat(no-1) +1
            mdstep_at_end_thermostat(no)   = mdstep_at_end_thermostat(no-1) +iret
         endif
      endif

      temp_at_start_thermostat(no) = tkb(no) /const_kB    ! default
      temp_at_end_thermostat(no) = tkb(no) /const_kB

      if ( temp_is_set ) return

      if ( f_getRealValue(tag_temp_s,dret,unit_f) == 0 ) then
         temp_at_start_thermostat(no) = dret
         tkb(no) = dret * CONST_kB
      endif
      if ( f_getRealValue(tag_temp_e,dret,unit_f) == 0 ) then
         temp_at_end_thermostat(no) = dret
         temp_is_set = .true.
      endif
! =============== 2014/06/10
      if ( .not. temp_is_set ) then
         if (printable) write(nfout,'(a)') '!** temperature undefined'
         call phase_error_with_msg(nfout,'temperature undefined',__LINE__,__FILE__)
      endif

    end subroutine set_qmass_and_temp

    subroutine set_pressure_profile(nfout,pressure_profile)
      integer, intent(in) :: nfout
      type(pressure_profile_t), intent(inout) :: pressure_profile
      integer :: i,j,k,iiret,iret
      integer :: nprof,currprof
      real(kind=DP) :: dret
      integer :: f_selectFirstTableLine,f_selectNextTableLine

      ! first count the number of profiles for the barostat
      if( f_selectBlock(tag_barostat) == 0) then
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
            i = i+1
         end do
         iiret = f_selectParentBlock()
      end if
      pressure_profile%nprof = i-1
      ! allocate
      nprof = pressure_profile%nprof
      allocate(pressure_profile%till_n(nprof));pressure_profile%till_n = -1
      allocate(pressure_profile%pressi(nprof));pressure_profile%pressi = 0.d0
      allocate(pressure_profile%pressf(nprof));pressure_profile%pressf = 0.d0
      allocate(pressure_profile%bmass(nprof));pressure_profile%bmass = -1.d0
      allocate(pressure_profile%tdamp_baro(nprof));pressure_profile%tdamp_baro = 500.d0*dtio
      pressure_profile%currprof = 0

      if( f_selectBlock(tag_barostat) == 0) then
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
            i = i+1
            pressure_profile%currprof = pressure_profile%currprof+1
            currprof = pressure_profile%currprof
            if(f_getIntValue(tag_till_n,iret) == 0) &
            &  pressure_profile%till_n(currprof) = iret
            if(f_getRealValue(tag_press,dret,'hartree/bohr3') == 0) then
               pressure_profile%pressi(currprof) = dret
               pressure_profile%pressf(currprof) = dret
            endif
            if(f_getRealValue(tag_pressi,dret,'hartree/bohr3') == 0) then
                pressure_profile%pressi(currprof) = dret
                pressure_profile%pressf(currprof) = dret
            endif
            if(f_getRealValue(tag_pressf,dret,'hartree/bohr3') == 0) &
            &  pressure_profile%pressf(currprof) = dret
            if(f_getRealValue(tag_tdamp_baro,dret,'au_time') == 0)then
                pressure_profile%tdamp_baro(currprof) = dret
            endif
            if(f_getRealValue(tag_m_baro,dret,'') == 0)then
               pressure_profile%bmass(currprof) = dret
            endif
         end do
         iiret = f_selectParentBlock()
         if(printable)then
            write(nfout,'(a)')         ' !** status of the pressure profile'
            write(nfout,'(a,i4)')      ' !** number of profiles : ',pressure_profile%nprof
            do k=1,pressure_profile%nprof
            write(nfout,'(a,2f20.15)') ' !** pressi, pressf  : '    &
            &                           , pressure_profile%pressi(k) &
            &                           , pressure_profile%pressf(k)
            write(nfout,'(a,2f20.3)') ' !** bmass, tdamp_baro: ', pressure_profile%bmass(k) &
            &                           , pressure_profile%tdamp_baro(k)
            enddo
         endif
      end if
    end subroutine set_pressure_profile

    subroutine specify_ityp()
      integer :: i, k, icount
      icount = 0
      do i = 1, natm
         type: do k = 1, ntyp
            if(trim(species_work(i)) == trim(speciesname(k))) then
               icount = icount + 1
               ityp(i) = k
               exit type
            end if
         end do type
      end do
      if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** specified atom number = ",i6)') icount
      if(icount < natm) then
         write(nfout,'(" icount = ", i8, " natm = ",i8)') icount, natm
         do i = 1, natm
            write(nfout,'(" species_work(",i5,")=",a)') i, species_work(i)
         end do
         do k = 1, ntyp
            write(nfout,'(" speciesname(",i5,")=",a)') k, speciesname(k)
         end do

         call phase_error_with_msg(nfout,&
         ' element names in the atom_list and in the element_list are inconsistent <<m_IS_rd_n.specify_ityp>>'&
         ,__LINE__,__FILE__)
      end if

    end subroutine specify_ityp

#ifdef __EDA__
! -----  ascat starts modifying  -----
    subroutine specify_icount_EDA

      integer :: i, k, icount

      allocate(icount_EDA(ntyp))
      do k = 1, ntyp
        icount = 0
        do i = 1, natm
          if(trim(speciesname(k)) == trim(species_work(i))) then
            icount = icount + 1
          endif
        enddo
        icount_EDA(k) = icount
      enddo
    end subroutine specify_icount_EDA
! -----  ascat ceases modifying  -----
#endif

    subroutine specify_ityp_vdw()
      integer :: i, k, icount

      icount = 0
      do i = 1, natm
         type: do k = 1, ntyp_vdw
            if(trim(species_vdw_work(i)) == trim(speciesname_vdw(k))) then
               icount = icount + 1
               ityp_vdw(i) = k
               exit type
            end if
         end do type
      end do
      if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** specified atom number = ",i6)') icount
      if(icount < natm) then
         call phase_error_with_msg(nfout,&
         ' vdW-type names in the atom_list and in the vdw_list are inconsistent <<m_IS_rd_n.specify_ityp_vdw>>'&
         ,__LINE__,__FILE__)
      end if

    end subroutine specify_ityp_vdw

    subroutine wd_atom_list()
      integer :: i
      if(ipriinputfile >= 3) then
         write(nfout,'(" !** === Atomic coordinates expressed in the internal system ===")')
         write(nfout,'(" !** id,  rx,    ry,    rz,    weight,  imdtyp, ityp,  species")')
         do i = 1, natm
            write(nfout,210) i,pos(i,1),pos(i,2),pos(i,3),iwei(i),imdtyp(i),ityp(i),species_work(i)
         end do
210      format(' !** ',i5,3f18.10,i3,i6,i3,3x,a4)
         write(nfout,'(" !** === Atomic coordinates expressed in the cartesian system ===")')
         write(nfout,'(" !** id,  rx,    ry,    rz,    weight,  imdtyp, ityp,  species")')
         do i = 1, natm
            write(nfout,210) i,cps(i,1),cps(i,2),cps(i,3),iwei(i),imdtyp(i),ityp(i),species_work(i)
         end do
      else
         write(nfout,'(" !** === Atomic coordinates ==")')
         write(nfout,'(" !**   id  ( coordinates_in_Intrnal_sys  ) (  coordinates_in_Cartsian_system  ) weight    ityp")')
         write(nfout,'(" !**       (   rx         ry         rz  ) (    rx          ry          rz    )      imdtyp    species")')
!!$                             -----                                 ------------------------------------   ------   ---
         do i = 1, natm
            write(nfout,211) i, pos(i,1:3), cps(i,1:3), iwei(i), imdtyp(i), ityp(i), species_work(i)
         end do
211      format(' !** ',i5,3f11.6,3f12.4,i3,i6,i3,3x,a4)
      end if

      do i = 1, ntyp
         write(nfout,'(" !** i = ",i5," element_name = ",a4)') i, speciesname(i)
      end do
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      if(.not.allocated(iatomn_EDA)) allocate(iatomn_EDA(natm))
      do i = 1, natm
        iatomn_EDA(i) = iatomn(ityp(i))
      enddo
! -----  ascat ceases modifying  -----
      endif
#endif

    end subroutine wd_atom_list

    subroutine count_species()
      integer :: i, k, kcount
      kcount = 0

      species_indp(1) = species_work(1)
      kcount = kcount + 1
      do i = 2, natm
         Registered :do k = 1, kcount
            if(species_indp(k) == species_work(i)) goto 1000
         end do Registered
         kcount = kcount + 1
         species_indp(kcount) = species_work(i)
1000     continue
      end do
      ntyp = kcount
      if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** ntyp = ",i6, " << count_species>>")') ntyp
    end subroutine count_species

    subroutine count_species_vdw()
      integer :: i, k, kcount
      kcount = 0

! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!     if(sw_pair_vdw==OFF) then
      if(sw_vdw_correction==OFF .or. vdw_method == VDW_DFTD3) then
! ==============================================================================
         ntyp_vdw = 0
         return
      end if

      species_vdw_indp(1) = species_vdw_work(1)
      kcount = kcount + 1
      do i = 2, natm
         Registered :do k = 1, kcount
            if(species_vdw_indp(k) == species_vdw_work(i)) goto 1001
         end do Registered
         kcount = kcount + 1
         species_vdw_indp(kcount) = species_vdw_work(i)
1001     continue
      end do
      ntyp_vdw = kcount
      if(ipriinputfile >= 2 .and. printable) write(nfout,'(" !** ntyp_vdw = ",i6, " << count_species>>")') ntyp_vdw
    end subroutine count_species_vdw

    subroutine set_input_coordinate_system(rstr)
      character(*),intent(in) :: rstr
      call strncmp0(trim(rstr),tag_cartesian,tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_XYZ,tf)
      if(tf) then
         input_coordinate_system = CARTS
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_pucv,tf)
      if(tf) then
         input_coordinate_system = PUCV
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_internal,tf)
      if(tf) then
         input_coordinate_system = PUCV
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_relative,tf)
      if(tf) then
         input_coordinate_system = PUCV
         goto 1001
      end if
      if(ipriinputfile >= 0 .and. printable) &
           & write(nfout,'(" !* input_coordinate_system is not defined properly in the input file!")')
      if(ipriinputfile >= 0 .and. printable) &
           & write(nfout,'(" !* So default value of PUCV is assinged")')
1001  continue
    end subroutine set_input_coordinate_system

    subroutine set_howtogive_coordinates(rstr)
      character(*),intent(in) :: rstr
      call strncmp0(trim(rstr),tag_directin,tf)
      if(tf) then
         howtogive_coordinates = DIRECTIN
         goto 1001
      end if
      call strncmp0(trim(rstr),tag_from_endpoint_images,tf)
      if(.not.tf) call strncmp0(trim(rstr),tag_from_endpoints,tf)
      if(tf) then
         howtogive_coordinates = FROM_ENDPOINTS
         goto 1001
      end if
1001  continue
    end subroutine set_howtogive_coordinates

    subroutine set_displacement()
      character(len=FMAXUNITLEN) :: unit_f
      integer :: f_getIntValue, f_getRealValue
      integer :: iret
      real(kind=DP) :: dret

      unit_f = ''
      if(input_coordinate_system == CARTS) unit_f = 'bohr'

      if(f_getIntValue(tag_sw_displace_atom,iret)==0) sw_displace_atom = iret
      displaced_atom = 0         ! defalut value
      displacement(1:3) = 0.d0   ! defalut values

! === KT_add === 2015/03/14
      if ( input_coordinate_system == CARTS ) sw_displacement_in_carts = ON
      if ( f_getIntValue( tag_sw_displacement_in_carts, iret ) == 0 ) then
         sw_displacement_in_carts = iret
      endif
      if ( sw_displacement_in_carts == ON ) then
         unit_f = 'bohr'
      else
         unit_f = ''
      endif
      write(nfout,*) '!** sw_displacement_in_carts = ', sw_displacement_in_carts
! ============= 2015/03/14

      if(sw_displace_atom == ON) then
         if(f_getIntValue(tag_displaced_atom,iret)==0) displaced_atom = iret
         if(f_getRealValue(tag_ux,dret,unit_f) == 0) displacement(1) = dret
         if(f_getRealValue(tag_uy,dret,unit_f) == 0) displacement(2) = dret
         if(f_getRealValue(tag_uz,dret,unit_f) == 0) displacement(3) = dret
      end if

    end subroutine set_displacement

    subroutine set_displacement2()
      integer :: i, istat
      real(kind=DP) :: cvec(3), rltv_t(3,3)

      if(displaced_atom < 1 .or. displaced_atom > natm) then
         if(printable) &
              & write(nfout,'("Because ",i4,"-th atom does not exist, no atom have been displaced.")') displaced_atom
         displaced_atom = 0
         displacement(1:3) = 0.d0
      else
#if 0
         pos(displaced_atom,1:3) =  pos(displaced_atom,1:3) + displacement(1:3)
         if(printable) write(nfout,'(i4,"-th atom was displaced; the displacement =",3(1x,f10.5))') &
              & displaced_atom, displacement(1:3)
         if(input_coordinate_system == PUCV) then
            do i=1,3
               work(1,i) = sum(p2bmat(1:3,i)*displacement(1:3))
            end do
            displacement(1:3) = work(1,1:3)
            do i=1,3
               work(1,i) = sum(altv(i,1:3)*displacement(1:3))
            end do
            displacement(1:3) = work(1,1:3) ! the displacement vector is expressed in the catesian system
         end if
#else
         if(printable) then
            write(nfout,'(i4,"-th atom was displaced; the displacement =",3(1x,f10.5))') &
                 &  displaced_atom, displacement(1:3)
         endif
         if (input_coordinate_system == PUCV ) then
            if ( sw_displacement_in_carts == OFF ) then
               pos(displaced_atom,1:3) = pos(displaced_atom,1:3) +displacement(1:3)
            else
               rltv_t = transpose(rltv) /PAI2
               call change_of_coordinate_system(rltv_t,displacement,1,1,cvec)
               do i=1,3
                  work(1,i) = sum(b2pmat(1:3,i)*cvec(1:3))
               end do
               cvec(1:3) = work(1,1:3)
               pos(displaced_atom,1:3) = pos(displaced_atom,1:3) +cvec(1:3)
            endif

         else if ( input_coordinate_system == CARTS ) then
            if ( sw_displacement_in_carts == ON ) then
               pos(displaced_atom,1:3) = pos(displaced_atom,1:3) +displacement(1:3)
            else
               do i=1,3
                  work(1,i) = sum(p2bmat(1:3,i)*displacement(1:3))
               end do
               displacement(1:3) = work(1,1:3)
               do i=1,3
                  work(1,i) = sum(altv(i,1:3)*displacement(1:3))
               end do
               displacement(1:3) = work(1,1:3)
                      ! the displacement vector is expressed in the catesian system
               pos(displaced_atom,1:3) = pos(displaced_atom,1:3) + displacement(1:3)
            end if
         endif
#endif
      end if

    end subroutine set_displacement2

    subroutine set_vibrational_mode()
      character(len=FMAXUNITLEN) :: unit_f = 'bohr'
      integer :: f_getIntValue, f_getRealValue
      integer :: iret
      real(kind=DP) :: dret

      if(f_getIntValue(tag_sw_vibrational_mode,iret)==0) sw_vibrational_mode = iret
      if(sw_vibrational_mode == ON) then
         if(f_getIntValue(tag_mode_index,iret)==0) mode_index = iret
         if(f_getRealValue(tag_normal_coordinate,dret,unit_f) == 0) normal_coordinate = dret
         if(f_getIntValue(tag_with_mode_effchg,iret)==0) with_mode_effchg = iret
      end if

      if(printable) then
         write(nfout,'(" !** sw_vibrational_mode   = ",i1)') sw_vibrational_mode
         write(nfout,'(" !** mode_index            = ",i4)') mode_index
         write(nfout,'(" !** normal_coordinate (Q) = ",f10.5)') normal_coordinate
         write(nfout,'(" !** with_mode_effchg      = ",i1)') with_mode_effchg
      end if

    end subroutine set_vibrational_mode

    subroutine set_vibrational_mode2()
      real(kind=DP) :: dret
      integer :: i, j

      call m_Files_open_nfmode()
      allocate(xi_mode(natm,3))
      call read_mode_vector(nfmode,natm,xi_mode)
      if(input_coordinate_system == PUCV) then
         do i=1,natm
            work(i,1:3) = (rltv(1,1:3)*xi_mode(i,1) &
                 & + rltv(2,1:3)*xi_mode(i,2) &
                 & + rltv(3,1:3)*xi_mode(i,3) )/PAI2
         end do
         do j = 1, 3
            do i = 1, natm
               work(i,j) = sum(b2pmat(:,j)*work(i,:))
            end do
         end do
      else
         work = xi_mode
      end if
      do i=1,natm
         pos(i,1:3) = pos(i,1:3) + work(i,1:3)
      end do

      if(printable) then
         do i = 1, natm
            write(nfout,'(" !** i = ",i4,"  pos = ",3f8.4, " <<set_vibrational_mode2>>")') &
                 & i, pos(i,1),pos(i,2),pos(i,3)
         end do
      end if
    end subroutine set_vibrational_mode2

    subroutine read_mode_vector(nfmode,natm,xi)
      integer, intent(in) :: nfmode
      integer, intent(in) :: natm
      real(kind=DP), intent(out) :: xi(natm,3)

      integer :: im,ia,idummy
      real(kind=DP) :: xi_tmp(natm,3),mass(natm)

      if(mype == 0) then
      rewind nfmode
      do i=1,6
         read(nfmode,*)
      end do
      do ia=1,natm
         read(nfmode,*) idummy,xi_tmp(ia,1:3),mass(ia)
      end do
      read(nfmode,*)
      read(nfmode,*)
      do im=1,natm*3
         read(nfmode,*)
         read(nfmode,*)
         do ia=1,natm
            read(nfmode,*) idummy,xi_tmp(ia,1:3)
         end do
         if(im == mode_index) xi(1:natm,1:3) = xi_tmp(1:natm,1:3)
         if(with_mode_effchg == YES) then
            read(nfmode,*)
            read(nfmode,*)
         end if
      end do
      end if
      if(npes > 1 ) then
         call mpi_bcast(xi,natm*3 &
                     & ,mpi_double_precision,0,MPI_CommGroup,ierr)
         call mpi_bcast(mass,natm &
                     & ,mpi_double_precision,0,MPI_CommGroup,ierr)
      end if

      if(printable) then
         write(nfout,'(1x,i4,"-th normal mode eigenvector:")') mode_index
         write(nfout,'(3x,"ia",5x,"x",15x,"y",15x,"z")')
         do ia=1,natm
            write(nfout,'(1x,i4,3(1x,f15.10))') ia,xi(ia,1:3)
         end do
      end if

      do ia=1,natm
         xi(ia,1:3) = normal_coordinate*xi(ia,1:3)/sqrt(mass(ia))
      end do

      if(printable) then
         write(nfout,'(1x,"Atomic displacements when Q = ",f10.5,":")') &
              & normal_coordinate
         write(nfout,'(3x,"ia",5x,"x",16x,"y",16x,"z")')
         do ia=1,natm
            write(nfout,'(1x,i4,3(1x,f15.10))') ia,xi(ia,1:3)
         end do
      end if

    end subroutine read_mode_vector

    subroutine set_replica_input_method(prealloc,m,iret)
      logical, intent(in) :: prealloc
      integer, intent(in) :: m
      integer, intent(out) :: iret
      integer :: i, ip, rint
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue

      i = 1
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /=0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         if(.not.prealloc) then
            if(i > m) exit
            ip = i
            iret = f_getIntValue(tag_replica_numbers,rint)
            if(iret == 0) ip = rint
            if( f_getStringValue(tag_howtogive_coordinates,rstr,LOWER) == 0) then
               call strncmp0(trim(rstr),tag_proportional, tf)
               if(tf) then
                  replica_howtogive_coordinates(ip) = PROPORTIONAL
                  goto 1001
               end if
               call strncmp0(trim(rstr),tag_file,tf)
               if(tf) then
                  replica_howtogive_coordinates(ip) = FILE
                  goto 1001
               end if
1001           continue
            end if
            if( f_getIntValue(tag_end0,rint) == 0 ) replica_endpoints(1,ip) = rint
            if( f_getIntValue(tag_end1,rint) == 0 ) replica_endpoints(2,ip) = rint
         end if
         i = i + 1
      end do
      iret = i - 1
    end subroutine set_replica_input_method

    subroutine set_atompos_and_etc(prealloc, m, iret)
!   Partially revised for the transformation from Bravais to Primitive system
!   by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003.
!      use m_Crystal_Structure, only : p2bmat ! inverse transformation matrix
      use m_Files, only : nfmode,m_Files_open_nfmode
      logical, intent(in) :: prealloc
      integer, intent(in) :: m
      integer, intent(out) :: iret
      integer :: i, rint, ip, len_rstrtrimmed, n_cp_name, j, istat = 0
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
      character(len=FMAXVALLEN) :: rstr
      real(kind=DP) :: rx, ry, rz, dret
      real(kind=DP) :: vx, vy, vz
      character(len=FMAXUNITLEN) :: unit_f
      logical :: tf

      real(kind=DP), allocatable, dimension(:,:) :: work
      integer :: ir
      integer :: nthermocount

      unit_f = ''
      if(input_coordinate_system == CARTS) unit_f = 'bohr'

      i = 1
      nthermocount=0
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
!!$         if(printable) write(nfout,'(" !!! i = ",i6," <<set_atompos_and_etc>>")') i
         if(.not.prealloc) then
            if(i > m) exit
            ip = i
            iret = f_getIntValue(tag_id,rint)
            if(iret == 0) ip = rint
            if(imdalg==VERLET .or. imdalg==T_CONTROL.or.imdalg==VELOCITY_SCALING.or.imdalg==PT_CONTROL &
              & .or. sw_calc_force == ON .or. (sw_calc_force == OFF .and. sw_vibrational_modes == ON))then
               imdtyp(ip) = ON
               imdtypxyz(ip,1:3) = ON
            endif
            if( f_getIntValue(tag_mobile, rint) == 0) then
               imdtyp(ip) = rint
               imdtypxyz(ip,1) = rint
               imdtypxyz(ip,2) = rint
               imdtypxyz(ip,3) = rint
            endif
            if( f_getIntValue(tag_mobilex, rint) == 0) imdtypxyz(ip,1) = rint
            if( f_getIntValue(tag_mobiley, rint) == 0) imdtypxyz(ip,2) = rint
            if( f_getIntValue(tag_mobilez, rint) == 0) imdtypxyz(ip,3) = rint
#ifndef _EMPIRICAL_
            if( f_getIntValue(tag_pdos, rint) == 0) if_pdos(ip) = rint
            if( f_getIntValue(tag_aldos, rint) == 0 &
                 &  .or. f_getIntValue(tag_atom_decomp, rint) == 0 ) if_aldos(ip) = rint
            if( f_getIntValue(tag_proj_group, rint) == 0) iproj_group(ip) = rint
#endif
            if( f_getRealValue(tag_mass, dret, 'au_mass') == 0) ionic_mass(ip) = dret
            if( f_getIntValue(tag_a_weight, rint) == 0) then
               if(inversion_symmetry == 0) then
                  if(rint /= 1) then
                     if(printable) then
                        write(nfout,'(" !** iwei(",i5,") = ",i3)') ip,rint
                        write(nfout,'(" !* iwei should be 1 when sw_inversion == OFF")')
                     end if
                  end if
                  iwei(ip) = 1
               else
                  if(rint > 2 .or. rint < 1) then
                     if(printable) then
                        write(nfout,'(" !** iwei(",i5,") = i3")') ip,rint
                        write(nfout,'(" !* iwei should be 1 or 2 ")')
                     end if
                  else
                     iwei(ip) = rint
                  end if
               end if
            end if

            if( f_getIntValue(tag_num_layer, rint) == 0) numlay(ip) = rint ! layer_dos
            if ( f_getIntValue(tag_key, rint) == 0 ) then
               atom_key(ip) = rint
            else if ( f_getIntValue(tag_atom_key, rint) == 0 ) then
               atom_key(ip) = rint
            endif

            if( f_getStringValue(tag_element, rstr, NOCONV) == 0) then
               len_rstrtrimmed = len_trim(rstr)
               if(len_rstrtrimmed > LEN_ATOMNAME) then
                  if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* element name is larger than ",i6)') LEN_ATOMNAME
                  n_cp_name = LEN_ATOMNAME
               else
                  n_cp_name = len_rstrtrimmed
               end if
               species_work(ip)(1:n_cp_name) = rstr(1:n_cp_name)
            end if

            if( f_getStringValue(tag_vdw, rstr, NOCONV) == 0) then
               len_rstrtrimmed = len_trim(rstr)
               if(len_rstrtrimmed > LEN_ATOMNAME) then
                  if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* vdW-type name is larger than ",i6)') LEN_ATOMNAME
                  n_cp_name = LEN_ATOMNAME
               else
                  n_cp_name = len_rstrtrimmed
               end if
               species_vdw_work(ip)(1:n_cp_name) = rstr(1:n_cp_name)
            end if

! === KT_add === 2014/12/29
            if ( sw_set_initial_magmom_by_atom == ON ) then
               call set_initial_ioncharge_atoms( ip )
               call set_initial_magmom_atoms( ip )
            endif
! ============== 2014/12/29

            if( endpoint_images == NO .or. howtogive_coordinates == DIRECTIN) then
               if( f_getRealValue( tag_rx, rx, unit_f) == 0) pos(ip,1) = rx
               if( f_getRealValue( tag_ry, ry, unit_f) == 0) pos(ip,2) = ry
               if( f_getRealValue( tag_rz, rz, unit_f) == 0) pos(ip,3) = rz
               if( f_getRealValue( tag_vx, vx, '') == 0) cpd_l(ip,1) = vx
               if( f_getRealValue( tag_vy, vy, '') == 0) cpd_l(ip,2) = vy
               if( f_getRealValue( tag_vz, vz, '') == 0) cpd_l(ip,3) = vz
            end if

            if(imdalg == T_CONTROL .or. imdalg == PT_CONTROL .or. &
              ((imdalg == VERLET .or. imdalg==P_CONTROL).and. t_ctrl_method == LANGEVIN)) then
!        === modified by T. Y. 2015/01/16 ===
#ifdef THERMOGROUP_DEFAULT_SETTING_FOR_NO_SPECIFIED_ATOMS
               if(imdtyp(ip)/=0)then
                 imdtyp(ip) = NOSE_HOOVER + 1 ! the default is to use the first thermostat !
               endif
#endif
!        ====================================
               tf = f_getIntValue(tag_thermo_group, rint) == 0
               if(.not.tf) tf = f_getIntValue(tag_thermo_g, rint) == 0
               if(tf) then
                  if(imdtyp(ip) /= 0) then
                     imdtyp(ip) = NOSE_HOOVER + rint
!                     imdtypxyz(ip,:) = NOSE_HOOVER + rint
                     nthermocount=nthermocount+1
                  end if
               end if
            end if
            if(sw_atom_excludable==ON)then
                if (f_getIntValue(tag_exclusion_target,iret)==0) then
                    exclusion_target(ip) = iret
                endif
            endif

         end if
         i = i + 1
      end do
      iret = i - 1

!      if(.not.prealloc.and.(imdalg == T_CONTROL .or. imdalg==PT_CONTROL) .and. nthermocount==0) then
!          stop '!** at least one atom must be assigned to a thermostat when using temperature_control'
!      endif
      if(.not.prealloc) then
         natm2 = 0
         do i = 1, natm
            if(iwei(i) <= 0) then
               if(ipriinputfile >= 0 .and. printable) &
                    & write(nfout,'(" !* iwei(",i4,") should be positive value")') i
               call phase_error_with_msg(nfout,' stopped due to illegal value of iwei <<set_atompos_and_etc(m_IS)>>'&
                                        ,__LINE__,__FILE__)
            end if
            natm2 = natm2 + iwei(i)
         end do
         if(printable) write(nfout,'(" !** natm2, natm = ",2i6)') natm2,natm

#ifndef _EMPIRICAL_
         if(sw_pdos == ON .or. sw_orb_popu == ON) then
            iret = 0
            do i=1,natm
               iret = iret + if_pdos(i)
            end do
            if(iret == 0) call phase_error_with_msg(nfout,'if_pdos are 0.',__LINE__,__FILE__)
            if(printable) then
               if(iret == natm) then
                  write(nfout,'(" !** pdos switches for all atoms are ON <<set_atompos_and_etc>>")')
               else
                  do i = 1, natm
                     write(nfout,'(" !** i =",i5," pdos=",i1, " <<set_atompos_and_etc>>")') &
                          & i,if_pdos(i)
                  end do
               end if
            end if
         end if

         if(printable) then
            if(sw_aldos == ON .or. sw_calc_wf_atom_decomposition == ON ) then
               iret = 0
               do i = 1, natm
                  iret = iret + if_aldos(i)
               end do
               if(iret == 0) then
                  write(nfout,'(" !** aldos switches for all atoms are OFF <<set_atompos_and_etc>>")')
               else if(iret == natm) then
                  write(nfout,'(" !** aldos switches for all atoms are ON <<set_atompos_and_etc>>")')
               else
                  write(nfout,'(" !** aldos <<set_atompos_and_etc>>")')
                  write(nfout,'(" !** ",10i8)') (if_aldos(i),i=1,natm)
!!$               do i=1,natm
!!$                  write(nfout,'(" !** i =",i5," aldos=",i1)') i,if_aldos(i)
!!$               end do
               end if
            end if
         end if

         if(sw_hubbard == ON) then
            do i=1,natm
               ig = iproj_group(i)
               if(ig==0) cycle
               do j=1,num_proj_elems(ig)
                  ip=proj_group(j,ig)
                  if(proj_attribute(ip)%strong_correlated) then
                     ihubbard(i) = ip
                  end if
               end do
            end do
         end if

! ============================ added by K. Tagami ================ 11.0
         if ( SpinOrbit_mode == ByProjector ) then
            do i=1,natm
               ig = iproj_group(i)
               if (ig==0) cycle
               do j=1,num_proj_elems(ig)
                  ip=proj_group(j,ig)
!!!!!!!!!                  if ( proj_attribute(ip)%activate_soc ) then

                     itab_spinorbit_addition(i) = ip
!                     write(*,*) 'itab_spinorb = ', i, itab_spinorbit(i)

!!!!!!!!!!!!!!!!!!                  end if
               end do
            end do
         end if
! ================================================================ 11.0

         if(num_projectors > 0.and.printable) then
            iret = 0
            do i=1,natm
               if(ipri>=2 .or. (ipri==1 .and. iproj_group(i) > 0)) then
                  iret = iret+1
                  if(sw_hubbard == ON) then
                     write(nfout,'(" !** i =",i5," iproj_group=",i3," ihubbard=",i3)') i,iproj_group(i),ihubbard(i)
                  else
                     write(nfout,'(" !** i =",i5," iproj_group=",i3)') i,iproj_group(i)
                  end if
               end if
            end do
            if(ipri==1 .and. iret < natm) then
               write(nfout,'(" !** the other iproj_groups being not printed here are all 0")')
            end if

! =================================== added by K. Tagami ========== 11.0
            if ( SpinOrbit_mode == ByProjector ) then
               write(nfout,*) '! ---- Info of projector for Spin-Orbit '
               do i=1,natm
                  write(nfout,'(" !** i =",i5," iproj_group=",i3, &
                       & " itab_spinorbit_addition =",i3)') &
                       &         i,iproj_group(i),itab_spinorbit_addition(i)
               end do
            end if
! ================================================================== 11.0

         end if
#endif

         allocate(work(natm,3), stat=istat)
         if (istat /= 0) then
            if(printable) write(nfout,*) 'Allocation error in sub.set_atompos_and_etc',natm,istat
            call phase_error_with_msg(nfout,'Allocation error in sub.set_atompos_and_etc',__LINE__,__FILE__)
         end if

         natm2 = 0
         do i = 1, natm
            if(iwei(i) <= 0) then
               if(ipriinputfile >= 0 .and. printable) &
                    & write(nfout,'(" !* iwei(",i4,") should be positive value")') i
               call phase_error_with_msg(nfout,' stopped due to illegal value of iwei <<set_atompos_and_etc(m_IS)>>'&
                                        ,__LINE__,__FILE__)
            end if
            natm2 = natm2 + iwei(i)
         end do
         if(printable) write(nfout,'(" !** natm2, natm = ",2i6," <<set_atompos_and_etc>>")') natm2,natm

!!$         if(input_coordinate_system == PUCV) then
!!$            do j = 1, 3
!!$               do i = 1, natm
!!$                  work(i,j) = sum(p2bmat(:,j)*pos(i,:))
!!$               end do
!!$            end do
!!$            pos = work
!!$            call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
!!$!!$            if(printable) then
!!$!!$               do i = 1, natm
!!$!!$                  write(nfout,*) cps(i,:)
!!$!!$               end do
!!$!!$            end if
!!$         else if(input_coordinate_system == CARTS) then
!!$            cps = pos
!!$            allocate(rltv_t(3,3))
!!$            rltv_t = transpose(rltv)/PAI2
!!$            call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
!!$            deallocate(rltv_t)
!!$         end if
         deallocate(work)
      end if
    end subroutine set_atompos_and_etc

    subroutine set_regions(nfout)
      integer, intent(in) :: nfout
      character(len=FMAXVALLEN) :: rstr
      character(len=256) :: str,strstr
      integer :: i,j,ip,iret
      real(kind=DP) :: dret
      integer :: f_selectBlock,f_selectParentBlock, &
   &             f_getIntValue,f_getRealValue,f_getStringValue, &
   &             f_selectFirstTableLine, f_selectNextTableLine
      i=0
      do while(.true.)
        i = i+1
        write(str,*) i
        strstr = trim(adjustl(str))
        if(f_selectBlock(tag_region//strstr)==0)then
          iret = f_selectParentBlock()
        else
          exit
        endif
      enddo
      num_regions = i-1
      if (num_regions < 1) return
      allocate(regions(num_regions))
      do i=1,num_regions
         regions(i)%no = i
         regions(i)%region_group = i
         regions(i)%region_type = BOX
         regions(i)%xmin = -1.d+10
         regions(i)%xmax =  1.d+10
         regions(i)%ymin = -1.d+10
         regions(i)%ymax =  1.d+10
         regions(i)%zmin = -1.d+10
         regions(i)%zmax =  1.d+10
         regions(i)%orientation = 3
         regions(i)%cylx = altv(1,1)
         regions(i)%cyly = altv(2,2)
         regions(i)%cylzmin = -1.d+10
         regions(i)%cylzmax =  1.d+10
         regions(i)%sigma = 1.d0
         regions(i)%epsi = 1.d-3
         regions(i)%energy = 0.d0
         regions(i)%i1 = 1
         regions(i)%i2 = 2
         regions(i)%i3 = 3
         regions(i)%tally = .false.
         allocate(regions(i)%forc(natm,3));regions(i)%forc(natm,3) = 0.d0
         allocate(regions(i)%target_atoms(natm));regions(i)%target_atoms=0
         regions(i)%ntarget_atoms = 0

         write(str,*) i
         strstr = trim(adjustl(str))
         iret = f_selectBlock(tag_region//strstr)
         if(f_getStringValue(tag_type,rstr,LOWER)==0)then
            if(rstr == tag_cylinder) regions(i)%region_type = CYLINDER
         endif
         if(f_getIntValue(tag_region_group,iret)==0)then
            regions(i)%region_group = iret
         endif
         if(f_getIntValue(tag_sw_tally,iret)==0) then
            if(iret == ON) regions(i)%tally = .true.
         endif
         if(regions(i)%region_type == BOX) then
            if(f_getRealValue(tag_xmin,dret,'bohr')==0) regions(i)%xmin = dret
            if(f_getRealValue(tag_xmax,dret,'bohr')==0) regions(i)%xmax = dret
            if(f_getRealValue(tag_ymin,dret,'bohr')==0) regions(i)%ymin = dret
            if(f_getRealValue(tag_ymax,dret,'bohr')==0) regions(i)%ymax = dret
            if(f_getRealValue(tag_zmin,dret,'bohr')==0) regions(i)%zmin = dret
            if(f_getRealValue(tag_zmax,dret,'bohr')==0) regions(i)%zmax = dret
         endif
         if(regions(i)%region_type == CYLINDER) then
            if(f_getRealValue(tag_radius,dret,'bohr')==0) regions(i)%radius = dret
            if(f_getIntValue(tag_orientation,iret)==0) regions(i)%orientation = iret
            if(f_getRealValue(tag_cylx,dret,'bohr')==0) regions(i)%cylx = dret
            if(f_getRealValue(tag_cyly,dret,'bohr')==0) regions(i)%cyly = dret
            if(f_getRealValue(tag_cylzmin,dret,'bohr')==0) regions(i)%cylzmin = dret
            if(f_getRealValue(tag_cylzmax,dret,'bohr')==0) regions(i)%cylzmax = dret
            if(f_getRealValue(tag_sigma,dret,'bohr')==0) regions(i)%sigma = dret
            if(f_getRealValue(tag_epsilon,dret,'hartree')==0) regions(i)%epsi = dret
            if(regions(i)%orientation==1)then
               regions(i)%i1 = 2
               regions(i)%i2 = 3
               regions(i)%i3 = 1
            else if (regions(i)%orientation==2) then
               regions(i)%i1 = 3
               regions(i)%i2 = 1
               regions(i)%i3 = 2
            endif
         endif
         iret = f_selectParentBlock()
      enddo
      i = 1
      iret = f_selectBlock(tag_atomlist)
      iret = f_selectBlock(tag_atoms)
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         ip = i
         if(f_getIntValue(tag_id,iret)==0) ip = iret
         if(f_getIntValue(tag_region_group,iret) == 0)then
           do j=1,num_regions
              if(regions(j)%region_group == iret)then
                 regions(j)%ntarget_atoms = regions(j)%ntarget_atoms+1
                 regions(j)%target_atoms(regions(j)%ntarget_atoms) = ip
              endif
           enddo
         endif
         i = i+1
      enddo
      iret = f_selectParentBlock()
      iret = f_selectParentBlock()
      if(num_regions>0 .and. ipriinputfile>0 .and. printable)then
         write(nfout,'(a)')               ' !** region statistics'
         write(nfout,'(a,i5)')            ' !** num_regions = ',num_regions
         do i=1,num_regions
            write(nfout,'(a,i5)')         ' !** status for region group ',regions(i)%region_group
            if(regions(i)%region_type == BOX) then
               write(nfout,'(a)')         ' !** region type : BOX'
               write(nfout,*) '!** xmin, xmax  : ',regions(i)%xmin,regions(i)%xmax
               write(nfout,*) '!** ymin, ymax  : ',regions(i)%ymin,regions(i)%ymax
               write(nfout,*) '!** zmin, zmax  : ',regions(i)%zmin,regions(i)%zmax
            else if (regions(i)%region_type == CYLINDER)then
               write(nfout,'(a)')         ' !** region type     : CYLINDER'
               write(nfout,'(a,i3,a)')    ' !** orientation     : ',regions(i)%orientation &
                                      & , ' (1->x, 2->y, 3->z)'
               write(nfout,'(a,f20.10)')  ' !** radius          : ',regions(i)%radius
               write(nfout,'(a,2f20.10)') ' !** cylx,cyly       : ',regions(i)%cylx,regions(i)%cyly
               write(nfout,*)             '!** cylzmin,cylzmax : ',regions(i)%cylzmin,regions(i)%cylzmax
            endif
            write(nfout,'(a,2f20.10)')    ' !** sigma, epsilon  : ',regions(i)%sigma,regions(i)%epsi
            write(nfout,'(a,l2)')         ' !** tally           : ',regions(i)%tally
            write(nfout,'(a,i8)')         ' !** n target atoms  : ',regions(i)%ntarget_atoms
            write(nfout,'(8i8)')        (regions(i)%target_atoms(j), j=1,regions(i)%ntarget_atoms)
         enddo
      endif
    end subroutine set_regions

! ===== KT_add ==== 2014/12/29 & 13.2S
    subroutine set_initial_ioncharge_atoms( ip )
      integer, intent(in) :: ip
      real(kind=DP) :: dret

      if ( f_getRealValue( tag_ion_charge, dret, "" ) == 0 ) then
         ionic_charge_atoms(ip) = dret
         mag_moment0_atoms_is_defined = .true.
      endif
    end subroutine set_initial_ioncharge_atoms

    subroutine set_initial_magmom_atoms( ip )
      integer, intent(in) :: ip
      integer :: mode = 0

      real(kind=DP) :: cnorm, mx, my, mz, theta, phi
      real(kind=DP) :: mdx, mdy, mdz, c1

      cnorm = 0.0d0

      if ( noncol ) then

         if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
            cnorm = dret;  mode = 2
         else
            mx = 0.0d0; my = 0.0d0; mz = 0.0d0
            if( f_getRealValue( tag_mx, dret, '') == 0 ) then
               mx = dret;  mode = 1
            endif
            if( f_getRealValue( tag_my, dret, '') == 0 ) then
               my = dret;  mode = 1
            endif
            if( f_getRealValue( tag_mz, dret, '') == 0 ) then
               mz = dret;  mode = 1
            endif
            cnorm = sqrt( mx**2 + my**2 + mz**2 )
         endif

         if ( mode == 2 ) then
            mdx = 0.0d0;  mdy = 0.0d0;  mdz = 0.0d0
            theta = 0.0d0; phi = 0.0d0

            if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
               theta = dret

               if ( f_getRealValue( tag_phi, dret, "" ) == 0 )   phi = dret
               theta = theta /180.0d0 *PAI;  phi = phi /180.0d0 *PAI

               mdx = sin( theta ) *cos( phi )
               mdy = sin( theta ) *sin( phi )
               mdz = cos( theta )

            else
               if ( f_getRealValue( tag_mdx, dret, "" ) == 0 ) mdx = dret
               if ( f_getRealValue( tag_mdy, dret, "" ) == 0 ) mdy = dret
               if ( f_getRealValue( tag_mdz, dret, "" ) == 0 ) mdz = dret

               c1 = sqrt( mdx**2 +mdy**2 + mdz**2 )
               if ( c1 > 1.0D-4 ) then
                  mdx = mdx /c1
                  mdy = mdy /c1
                  mdz = mdz /c1
               else
                  mdx = 0.0d0;  mdy = 0.0d0;  mdz = 1.0d0
               endif
            endif
         endif

         if ( cnorm > 1.0E-8 ) then
            if ( mode == 1 ) then
               mag_moment0_atoms(ip,1) = mx
               mag_moment0_atoms(ip,2) = my
               mag_moment0_atoms(ip,3) = mz
            else if ( mode == 2 ) then
               mag_moment0_atoms(ip,1) = cnorm *mdx
               mag_moment0_atoms(ip,2) = cnorm *mdy
               mag_moment0_atoms(ip,3) = cnorm *mdz
            endif
         endif

         if ( sw_fix_global_quantz_axis == ON ) then
            mx = mag_moment0_atoms(ip,1)
            my = mag_moment0_atoms(ip,2)
            mz = mag_moment0_atoms(ip,3)
            cnorm = sqrt( mx**2 +my**2 +mz**2 )

            c1 = mag_moment0_atoms(ip,1) * Global_Quantz_Axis_Fixed(1) &
                 & +mag_moment0_atoms(ip,2) * Global_Quantz_Axis_Fixed(2) &
                 & +mag_moment0_atoms(ip,3) * Global_Quantz_Axis_Fixed(3)
            if ( c1 > 0.0 ) then
               mag_moment0_atoms(ip,1:3) =  cnorm *Global_Quantz_Axis_Fixed(1:3)
            else
               mag_moment0_atoms(ip,1:3) = -cnorm *Global_Quantz_Axis_Fixed(1:3)
            endif
         endif

      else
         if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
            cnorm = dret;  mode = 1
            if ( abs(cnorm) > 1.0E-8 ) then
               mag_moment0_atoms(ip,1) = dret
            endif
         endif
      endif

      if ( mode /= 0 ) then
         mag_moment0_atoms_is_defined = .true.
         if ( noncol ) then
            mx = mag_moment0_atoms(ip,1)
            my = mag_moment0_atoms(ip,2)
            mz = mag_moment0_atoms(ip,3)
            cnorm = sqrt( mx**2 +my**2 +mz**2 )
            if ( cnorm > 1.0D-4 ) then
               mag_direction0_atoms(ip,1) = mx /cnorm
               mag_direction0_atoms(ip,2) = my /cnorm
               mag_direction0_atoms(ip,3) = mz /cnorm
            else
               mag_direction0_atoms(ip,1) = 0.0d0
               mag_direction0_atoms(ip,2) = 0.0d0
               mag_direction0_atoms(ip,3) = 1.0d0
            endif
         endif
      endif
    end subroutine set_initial_magmom_atoms
! ================= 2014/12/29 & 13.2S

    subroutine set_atompos_and_etc_reservoir()
!   Partially revised for the transformation from Bravais to Primitive system
!   by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003.
!      use m_Crystal_Structure, only : p2bmat ! inverse transformation matrix
      use m_Files, only : nfmode,m_Files_open_nfmode
      integer :: i, rint, ip, len_rstrtrimmed, n_cp_name, j, istat = 0
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue,f_getNumRows
      character(len=FMAXVALLEN) :: rstr
      real(kind=DP) :: rx, ry, rz, dret
      real(kind=DP) :: vx, vy, vz
      character(len=FMAXUNITLEN) :: unit_f
      logical :: tf
      real(kind=DP),allocatable,dimension(:,:) :: cpstt,postt
      real(kind=DP), allocatable, dimension(:,:) :: wk
      real(kind=DP), dimension(3,3) :: rltmp

      i=1
      do
         if(i==1)then
            if(f_selectFirstTableLine()/=0)then
               exit
            endif
         else
            if(f_selectNextTableLine()/=0)then
               exit
            endif
         endif
         i=i+1
      enddo
      natom_reservoir = i-1
      if(allocated(atom_reservoir)) deallocate(atom_reservoir)
      allocate(atom_reservoir(natom_reservoir))
      do i=1,natom_reservoir
         call init_atom(i,atom_reservoir(i))
      enddo

      unit_f = ''
      if(input_coordinate_system == CARTS) unit_f = 'bohr'

      i = 1
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
!!$         if(printable) write(nfout,'(" !!! i = ",i6," <<set_atompos_and_etc>>")') i
         if(i > natom_reservoir) exit
         ip = i
         iret = f_getIntValue(tag_id,rint)
         if(iret == 0) ip = rint
         if(imdalg==VERLET .or. imdalg==T_CONTROL.or.imdalg==VELOCITY_SCALING.or.imdalg==PT_CONTROL)then
            atom_reservoir(ip)%imdtyp = ON
            atom_reservoir(ip)%imdtypxyz(1:3) = ON
         endif
         if( f_getIntValue(tag_mobile, rint) == 0) then
            atom_reservoir(ip)%imdtyp = rint
            atom_reservoir(ip)%imdtypxyz(1:3) = rint
         endif
         if( f_getIntValue(tag_mobilex, rint) == 0) atom_reservoir(ip)%imdtypxyz(1) = rint
         if( f_getIntValue(tag_mobiley, rint) == 0) atom_reservoir(ip)%imdtypxyz(2) = rint
         if( f_getIntValue(tag_mobilez, rint) == 0) atom_reservoir(ip)%imdtypxyz(3) = rint
#ifndef _EMPIRICAL_
         if( f_getIntValue(tag_pdos, rint) == 0) atom_reservoir(ip)%if_pdos = rint
         if( f_getIntValue(tag_aldos, rint) == 0) atom_reservoir(ip)%if_aldos = rint
         if( f_getIntValue(tag_proj_group, rint) == 0) atom_reservoir(ip)%iproj_group = rint
#endif
         if( f_getRealValue(tag_mass, dret, 'au_mass') == 0) atom_reservoir(ip)%ionic_mass = dret
         if( f_getIntValue(tag_a_weight, rint) == 0) then
            if(inversion_symmetry == 0) then
               if(rint /= 1) then
                  if(printable) then
                     write(nfout,'(" !** iwei(",i5,") = ",i3)') ip,rint
                     write(nfout,'(" !* iwei should be 1 when sw_inversion == OFF")')
                  end if
               end if
               atom_reservoir(ip)%iwei = 1
            else
               if(rint > 2 .or. rint < 1) then
                  if(printable) then
                     write(nfout,'(" !** iwei(",i5,") = i3")') ip,rint
                     write(nfout,'(" !* iwei should be 1 or 2 ")')
                  end if
               else
                  atom_reservoir(ip)%iwei = rint
               end if
            end if
         end if

         if( f_getIntValue(tag_num_layer, rint) == 0) atom_reservoir(ip)%numlay = rint ! layer_dos

         if( f_getStringValue(tag_element, rstr, NOCONV) == 0) then
            len_rstrtrimmed = len_trim(rstr)
            if(len_rstrtrimmed > LEN_ATOMNAME) then
               if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* element name is larger than ",i6)') LEN_ATOMNAME
               n_cp_name = LEN_ATOMNAME
            else
               n_cp_name = len_rstrtrimmed
            end if
            atom_reservoir(ip)%element(1:n_cp_name) = rstr(1:n_cp_name)
         end if

         if( f_getStringValue(tag_vdw, rstr, NOCONV) == 0) then
            len_rstrtrimmed = len_trim(rstr)
            if(len_rstrtrimmed > LEN_ATOMNAME) then
               if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* vdW-type name is larger than ",i6)') LEN_ATOMNAME
               n_cp_name = LEN_ATOMNAME
            else
               n_cp_name = len_rstrtrimmed
            end if
            atom_reservoir(ip)%element_vdw(1:n_cp_name) = rstr(1:n_cp_name)
         end if

         if( f_getRealValue( tag_rx, rx, unit_f) == 0) atom_reservoir(ip)%pos(1) = rx
         if( f_getRealValue( tag_ry, ry, unit_f) == 0) atom_reservoir(ip)%pos(2) = ry
         if( f_getRealValue( tag_rz, rz, unit_f) == 0) atom_reservoir(ip)%pos(3) = rz
         if( f_getRealValue( tag_vx, vx, '') == 0) atom_reservoir(ip)%cpd_l(1) = vx
         if( f_getRealValue( tag_vy, vy, '') == 0) atom_reservoir(ip)%cpd_l(2) = vy
         if( f_getRealValue( tag_vz, vz, '') == 0) atom_reservoir(ip)%cpd_l(3) = vz

         if( f_getIntValue(tag_reservoir_group,iret)==0) atom_reservoir(ip)%group = iret
         if(imdalg == T_CONTROL .or. imdalg==PT_CONTROL .or. &
           ((imdalg == VERLET .or. imdalg == P_CONTROL) .and. t_ctrl_method == LANGEVIN)) then
            atom_reservoir(ip)%imdtyp = NOSE_HOOVER + 1
!            atom_reservoir(ip)%imdtypxyz(:) = NOSE_HOOVER + 1
            tf = f_getIntValue(tag_thermo_group, rint) == 0
            if(.not.tf) tf = f_getIntValue(tag_thermo_g, rint) == 0
            if(tf) then
               if(imdtyp(ip) /= 0) then
                  atom_reservoir(ip)%imdtyp = NOSE_HOOVER + rint
!                  atom_reservoir(ip)%imdtypxyz(:) = NOSE_HOOVER + rint
               end if
            end if
         end if
         if(sw_atom_excludable==ON) then
            if (f_getIntValue(tag_exclusion_target,iret)==0) &
          &     atom_reservoir(ip)%exclusion_target = iret
         endif
      i = i + 1
      end do
      iret = i - 1

#ifndef _EMPIRICAL_
      if(sw_hubbard == ON) then
         do i=1,natom_reservoir
            ig = atom_reservoir(i)%iproj_group
            if(ig==0) cycle
            do j=1,num_proj_elems(ig)
               ip=proj_group(j,ig)
               if(proj_attribute(ip)%strong_correlated) then
                  atom_reservoir(i)%ihubbard= ip
               end if
            end do
         end do
      end if
#endif
      do i=1,natom_reservoir
         do j=1,ntyp
            if(trim(atom_reservoir(i)%element)==trim(speciesname(j)))then
               atom_reservoir(i)%ityp = j
               exit
            endif
         enddo
      enddo

      allocate(wk(natom_reservoir,3))
      allocate(cpstt(natom_reservoir,3));cpstt=0.d0
      allocate(postt(natom_reservoir,3));postt=0.d0

      if(input_coordinate_system == PUCV) then
         do j = 1, 3
            do i = 1, natom_reservoir
               wk(i,j) = sum(p2bmat(:,j)*atom_reservoir(i)%pos(:))
            end do
         end do
         do i=1,natom_reservoir
            atom_reservoir(i)%pos(:) = wk(i,:)
         enddo
         do i=1,natom_reservoir
            postt(i,:) = atom_reservoir(i)%pos(:)
         enddo
         call change_of_coordinate_system(altv,postt,natom_reservoir,natom_reservoir,cpstt)!-(b_I.S.) pos -> cps
         do i=1,natom_reservoir
            atom_reservoir(i)%cps(:) = cpstt(i,:)
         enddo
      else if(input_coordinate_system == CARTS) then
         do i=1,natom_reservoir
            atom_reservoir(i)%cps(:) = atom_reservoir(i)%pos(:)
         enddo
         do i=1,natom_reservoir
            cpstt(i,:) = atom_reservoir(i)%cps(:)
         enddo
         rltmp = transpose(rltv)/PAI2
         call change_of_coordinate_system(rltmp,cpstt,natom_reservoir,natom_reservoir,postt) !-(b_I.S.) cps -> pos
      end if

      deallocate(wk)
      deallocate(cpstt)
      deallocate(postt)

      if(printable .and. iprimd>=1)then
         write(nfout,'(a)') ' !** atom reservoir defined in the input file '
         do i=1,natom_reservoir
            call print_atom(atom_reservoir(i))
         enddo
      endif
      call resolve_reservoir_group()
    end subroutine set_atompos_and_etc_reservoir

    subroutine resolve_reservoir_group()
       integer :: i,j,k
       integer, allocatable, dimension(:) :: inds0,inds1,inds2
       integer :: icount
       allocate(inds0(natom_reservoir))
       do i=1,natom_reservoir
          inds0(i) = i
       enddo

       do i=1,natom_reservoir
          if(atom_reservoir(i)%group>0)then
             inds0(i) = atom_reservoir(i)%group
          else
             atom_reservoir(i)%group = i
          endif
       enddo

       natom_group = 0
       allocate(inds1(natom_reservoir));inds1=0
       allocate(inds2(natom_reservoir));inds2=0
       do i=1,natom_reservoir
          if (.not.check_dupli(inds0(i),inds1,natom_reservoir))then
             natom_group = natom_group+1
             inds1(natom_group) = inds0(i)
          endif
          inds2(inds0(i)) = inds2(inds0(i))+1
       enddo

       allocate(natm_per_group(natom_group));natm_per_group = 0
       allocate(atomid_in_group(natom_group,natom_reservoir));atomid_in_group = 0
       do i=1,natom_group
          natm_per_group(i) = inds2(i)
       enddo
       do i=1,natom_group
          icount=0
          loopj:do j=1,natm_per_group(i)
             loopk:do k=1,natom_reservoir
                if(atom_reservoir(k)%group==i .and. &
                & .not.check_dupli(k,atomid_in_group(i,:),natom_reservoir)) then
                   atomid_in_group(i,j) = k
                   icount = icount+1
                   exit loopk
                endif
                if(icount==natm_per_group(i)) exit loopj
             enddo loopk
          enddo loopj
       enddo

       if(printable)then
          write(nfout,'(a,i5)') ' !** number of atom groups : ',natom_group
          do i=1,natom_group
             write(nfout,'(a,i5,a,i5)') ' !** number of atoms in group ',i,' is ',natm_per_group(i)
             do j=1,natm_per_group(i)
                call print_atom(atom_reservoir(atomid_in_group(i,j)))
             enddo
          enddo
       endif

       deallocate(inds0)
       deallocate(inds1)
       deallocate(inds2)
    end subroutine resolve_reservoir_group

    subroutine set_atm_exclusion_criteria()
      real(kind=DP) :: dret
      if(f_selectBlock(tag_atom_exclusion_criteria)==0)then
         if(f_getRealValue(tag_x_greater_than,dret,'bohr')==0)then
             exclusion_criteria_min(1) = dret
         endif
         if(f_getRealValue(tag_y_greater_than,dret,'bohr')==0)then
             exclusion_criteria_min(2) = dret
         endif
         if(f_getRealValue(tag_z_greater_than,dret,'bohr')==0)then
             exclusion_criteria_min(3) = dret
         endif
         if(f_getRealValue(tag_x_less_than,dret,'bohr')==0)then
             exclusion_criteria_max(1) = dret
         endif
         if(f_getRealValue(tag_y_less_than,dret,'bohr')==0)then
             exclusion_criteria_max(2) = dret
         endif
         if(f_getRealValue(tag_z_less_than,dret,'bohr')==0)then
             exclusion_criteria_max(3) = dret
         endif
      endif
    end subroutine set_atm_exclusion_criteria

    ! will return true if inds contain ind.
    logical function check_dupli(ind,inds,n)
      integer, intent(in) :: ind
      integer, dimension(:), intent(in) :: inds
      integer, intent(in) :: n
      integer :: i
      do i=1,n
        if(inds(i).eq.ind)then
          check_dupli = .true.
          return
        endif
      enddo
      check_dupli = .false.
    end function check_dupli

    subroutine set_atompos_from_endpoints()
      integer :: i, j
      do j = 1, 3
         do i = 1, natm
            pos(i,j) = (pos_end0(i,j) + pos_end1(i,j))*0.5d0
            cps(i,j) = (cps_end0(i,j) + cps_end1(i,j))*0.5d0
         end do
      end do

      if(printable) then
         do i = 1, natm
            write(nfout,'(" !** i = ",i4,"  pos = ",3f8.4)') i, pos(i,1),pos(i,2),pos(i,3)
         end do
         write(nfout,*)
         do i = 1, natm
            write(nfout,'(" !** i = ",i4,"  cps = ",3f8.4)') i, cps(i,1),cps(i,2),cps(i,3)
         end do
      end if
    end subroutine set_atompos_from_endpoints

    subroutine set_atompos2()
      integer :: i, j, iret

      if(ipriinputfile >= 3) then
         do i = 1, natm
            write(nfout,'(" !** i = ",i4,"  pos = ",3f10.4, " <<set_atompos2>>")') i, pos(i,1:3)
         end do
      end if

      if(input_coordinate_system == PUCV) then
         do j = 1, 3
            do i = 1, natm
               work(i,j) = sum(p2bmat(:,j)*pos(i,:))
            end do
         end do
         pos = work
         call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
      else if(input_coordinate_system == CARTS) then
         cps = pos
         allocate(rltv_t(3,3))
         rltv_t = transpose(rltv)/PAI2
         call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
         deallocate(rltv_t)
      end if
    end subroutine set_atompos2

    subroutine set_endpoint_atompos(m, pos_end, cps_end)
!   Partially revised for the transformation from Bravais to Primitive system
!   by BETSUYAKU, K. (Fuji Research Institute Co., Ltd.), July 2003.
!      use m_Crystal_Structure, only : p2bmat ! inverse transformation matrix
      integer, intent(in) :: m
      real(kind=DP), intent(out), dimension(m,3) :: pos_end, cps_end
      integer :: i, rint, ip,  j
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
      character(len=FMAXVALLEN) :: rstr
      real(kind=DP) :: rx, ry, rz, dret
      character(len=FMAXUNITLEN) :: unit_f

      unit_f = ''
      if(input_coordinate_system == CARTS) unit_f = 'bohr'

      i = 1
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if

         if(i > m) exit
         ip = i
         iret = f_getIntValue(tag_id,rint)
         if(iret == 0) ip = rint

         if( f_getRealValue( tag_rx, rx, unit_f) == 0) pos_end(ip,1) = rx
         if( f_getRealValue( tag_ry, ry, unit_f) == 0) pos_end(ip,2) = ry
         if( f_getRealValue( tag_rz, rz, unit_f) == 0) pos_end(ip,3) = rz

         i = i + 1
      end do

      if(input_coordinate_system == PUCV) then
         do j = 1, 3
            do i = 1, natm
               work(i,j) = sum(p2bmat(:,j)*pos_end(i,:))
            end do
         end do
         pos_end = work
         call change_of_coordinate_system(altv,pos_end,natm,natm,cps_end) !-(b_I.S.) pos -> cps
      else if(input_coordinate_system == CARTS) then
         cps_end = pos_end
         allocate(rltv_t(3,3))
         rltv_t = transpose(rltv)/PAI2
         call change_of_coordinate_system(rltv_t,cps_end,natm,natm,pos_end) !-(b_I.S.) cps -> pos
         deallocate(rltv_t)
      end if

    end subroutine set_endpoint_atompos

    subroutine set_endpoint_atompos_from_file(nfimage, m, pos_end, cps_end)

      integer, intent(in) :: nfimage, m
      integer :: realConvByUnit
      real(kind=DP), intent(out), dimension(m,3) :: pos_end, cps_end
      character(len=FMAXVALLEN) :: rstr
      integer i, j, id, iret
      character(10) token
      real(kind=DP) :: rx, ry, rz, dret
      character(len=FMAXUNITLEN) :: unit_f, unit_r

      unit_f = ''

      id = 0
      do
        read(nfimage,'(a)',end=1001) rstr
	if(index(rstr,'!') /= 0) rstr = rstr(1:index(rstr,'!')-1)
        if(len_trim(rstr) == 0) cycle

        if(index(rstr,'coordinate_system') /= 0) then
           call set_input_coordinate_system(trim(rstr(index(rstr,'=')+1:)))
	   write(nfout,*) 'coord_stystem: ', input_coordinate_system
	   cycle
        end if
        if(index(rstr,'#units') /= 0) then
	   read(rstr,*) token, unit_f
           write(nfout,*) 'unit: ', unit_f
	   cycle
        end if

	id = id + 1
	read(rstr,*) token, pos_end(id,1),pos_end(id,2),pos_end(id,3)
	write(nfout,'(i5,3f15.10)') id,pos_end(id,1),pos_end(id,2),pos_end(id,3)
      end do
1001  continue

      ! unit translation
      if(input_coordinate_system == PUCV) unit_f = 'bohr'
      unit_r = 'bohr'
      do i = 1, natm
	do j = 1, 3
	  iret = realConvByUnit(pos_end(i,j),pos_end(i,j),unit_f,unit_r)
        end do
      end do

      if(input_coordinate_system == PUCV) then
         allocate(work(natm,3))
         do j = 1, 3
            do i = 1, natm
               work(i,j) = sum(p2bmat(:,j)*pos_end(i,:))
            end do
         end do
         pos_end = work
         deallocate(work)
         call change_of_coordinate_system(altv,pos_end,natm,natm,cps_end) !-(b_I.S.) pos -> cps
      else if(input_coordinate_system == CARTS) then
         cps_end = pos_end
         allocate(rltv_t(3,3))
         rltv_t = transpose(rltv)/PAI2
         call change_of_coordinate_system(rltv_t,cps_end,natm,natm,pos_end) !-(b_I.S.) cps -> pos
         deallocate(rltv_t)
      end if

    end subroutine set_endpoint_atompos_from_file

    subroutine set_neb_convergence_condition(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf

      call strncmp2(rstr, FMAXVALLEN, 'delta_e', len('delta_e'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'1',1,tf)
      if(tf) then
         neb_convergence_condition = 1
         neb_convergence_threshold = 1.d-6
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, 'phase_force', len('phase_force'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'2',1,tf)
      if(tf) then
         neb_convergence_condition = 2
         neb_convergence_threshold = 1.d-3
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, 'neb_force', len('neb_force'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'3',1,tf)
      if(tf) then
         neb_convergence_condition = 3
         neb_convergence_threshold = 1.d-3
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, 'force_at_transition_state', &
	len('force_at_transition_state'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'4',1,tf)
      if(tf) then
         neb_convergence_condition = 4
         neb_convergence_threshold = 1.d-3
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, 'phase_force_normal', &
	len('phase_force_normal'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'5',1,tf)
      if(tf) then
         neb_convergence_condition = 5
         neb_convergence_threshold = 1.d-3
         goto 1001
      end if

1001  continue
    end subroutine set_neb_convergence_condition

    subroutine set_dimer_time_integral(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf

      call strncmp2(rstr, FMAXVALLEN, 'quench', len('quench'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'2',1,tf)
      if(tf) then
         dimer_time_integral = 2
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, 'fire', len('fire'), tf)
      if(tf) then
         dimer_time_integral = FIRE
         goto 1001
      end if
1001  continue
    end subroutine set_dimer_time_integral

    subroutine set_neb_time_integral(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf

      call strncmp2(rstr, FMAXVALLEN, 'steepest_descent', len('steepest_descent'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'12',1,tf)
      if(tf) then
         neb_time_integral = 12
         goto 1001
      end if

      call strncmp2(rstr, FMAXVALLEN, 'quench', len('quench'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'2',1,tf)
      if(tf) then
         neb_time_integral = 2
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, 'lbfgs', len('lbfgs'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'20',1,tf)
      if(tf) then
         neb_time_integral = L_BFGS
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, 'bfgs', len('bfgs'), tf)
      if(tf) then
         neb_time_integral = BFGS
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, 'cg', len('cg'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'11',1,tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN, 'cg2', len('cg2'), tf)
      if(.not.tf) call strncmp2(rstr, FMAXVALLEN,'17',1,tf)
      if(tf) then
         neb_time_integral = CG_STROPT
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, 'fire', len('fire'), tf)
      if(tf) then
         neb_time_integral = FIRE
         goto 1001
      end if

1001  continue
    end subroutine set_neb_time_integral

    subroutine set_fixed_bond_atoms(prealloc,num_fixed_bonds,countedbonds,fixed_atoms)
      ! Coded by T. Yamasaki (FUJITSU LABORATORIES LTD.), 2nd Aug. 2003
      logical,intent(in) :: prealloc
      integer,intent(in) :: num_fixed_bonds
      integer,intent(out) :: countedbonds
      integer,intent(out),optional :: fixed_atoms

      integer :: i, fix_type, rint1, rint2, icatm
      integer :: f_getIntValue, f_getStringValue, f_selectFirstTableLine, f_selectNextTableLine
      logical :: tf

      icatm = 0
      i = 1
      do while(.true.)
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         if(.not.prealloc) then
            if(i > num_fixed_bonds) exit
            if(imdalg == QUENCHED_CONSTRAINT) then
               fix_type = BONDLENGTH_FIX_1
               if( f_getStringValue(tag_type,rstr, LOWER) == 0) then
                  call strncmp0(trim(rstr),tag_absolute,tf)
                  if(tf) then
                     fix_type = BONDLENGTH_FIX_1
                     goto 1001
                  end if
                  call strncmp0(trim(rstr),tag_square,tf)
                  if(tf) then
                     fix_type = BONDLENGTH_FIX_2
                     goto 1001
                  end if
1001              continue
               end if
            else
               fix_type = BONDLENGTH_FIX
            end if
         end if
         if(ipriinputfile >= 1) then
            write(nfout,'(" !** i = ",i5)') i
            if(fix_type == BONDLENGTH_FIX) then
               write(nfout,'(" !**   fix_type = BONDLENGTH_FIX")')
            else if(fix_type == BONDLENGTH_FIX_1) then
               write(nfout,'(" !**   fix_type = BONDLENGTH_FIX_1")')
            else if(fix_type == BONDLENGTH_FIX_2) then
               write(nfout,'(" !**   fix_type = BONDLENGTH_FIX_2")')
            end if
         end if

         rint1 = 0
         rint2 = 0
         if( f_getIntValue(tag_atom1,rint1) == 0)  then
            if(rint1 > 0 .and. rint1 <= natm) then
               if(.not.prealloc) then
                  imdtyp(rint1) = fix_type
                  imdtypxyz(rint1,:) = fix_type
                  bondlength_fix_set(1,i) = rint1
               end if
               icatm = icatm + 1
            else
               rint1 = 0
            end if
         end if

         if( f_getIntValue(tag_atom2,rint2) == 0)  then
            if(rint2 > 0 .and. rint2 <= natm) then
               if(.not.prealloc) then
                  imdtyp(rint2) = fix_type
                  imdtypxyz(rint2,:) = fix_type
                  bondlength_fix_set(2,i) = rint2
               end if
               icatm = icatm + 1
            else
               rint2 = 0
            end if
         end if
         if(prealloc) then
            if(rint1 == rint2) then
               if(printable) write(nfout,'(" !* [ atom1 == atom2 (= ",i5 &
                    & ,")] is illegal <<m_IS_rd_n.set_fixed_bond>>")')
               call phase_error_with_msg(nfout,' [atom1 = atom2, this is illegal <<m_IS_rd_n.set_fixed_bond>>', &
               __LINE__,__FILE__)
            else if(rint1 == 0 .or. rint2 == 0) then
               if(printable) write(nfout,'(" !* input atom number is illegal" &
                    & , " <<m_IS_rd_n.set_fixed_bond>>")')
               call phase_error_with_msg(nfout,'  atom1 or atom2 is illegal <<m_IS_rd_n.set_fixed_bond>>', &
               __LINE__,__FILE__)
            end if
         end if
         i = i + 1
      end do
      countedbonds = i - 1
      if(prealloc) fixed_atoms = icatm
    end subroutine set_fixed_bond_atoms

    subroutine set_element_detail()
      integer :: i, iret, rint, ip, len_rstrtrimmed, n_cp_name, icount, ia
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
!!$      character(len=FMAXVALLEN) :: rstr
      real(kind=DP) :: dret
      integer, allocatable, dimension(:) :: elemmap
      character(len=LEN_ATOMNAME),allocatable, dimension(:) ::  speciesname_wk
      real(kind=DP), allocatable, dimension(:) :: iatomn_wk,amion_wk
      integer :: lens

!!$      character(len=FMAXUNITLEN) unit_f
!!$      unit_f = ''

#ifdef forsafe
      goto 100
#endif
      allocate(elemmap(ntyp));elemmap=0
      allocate(speciesname_wk(ntyp))
      allocate(iatomn_wk(ntyp))
      allocate(amion_wk(ntyp))
      do i = 1, ntyp
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         ip = i
         if( f_getStringValue(tag_element, rstr, NOCONV) == 0) then
            len_rstrtrimmed = len_trim(rstr)
            if(len_rstrtrimmed > LEN_ATOMNAME) then
               if(ipriinputfile >= 1 .and. printable) &
                    & write(nfout,'(" !* element name is larger than ",i6)') LEN_ATOMNAME
               n_cp_name = LEN_ATOMNAME
            else
               n_cp_name = len_rstrtrimmed
            end if
            speciesname_wk(ip)(1:n_cp_name) = rstr(1:n_cp_name)
         end if
      enddo
      call build_elemmap(ntyp,speciesname,speciesname_wk,elemmap)
      speciesname_wk = ''
      do i=1,ntyp
         lens = len(speciesname(i))
         speciesname_wk(elemmap(i))(1:lens) = speciesname(i)(1:lens)
         iatomn_wk(elemmap(i)) = iatomn(i)
         amion_wk(elemmap(i)) = amion(i)
      enddo
      speciesname = speciesname_wk
      iatomn = iatomn_wk
      amion = amion_wk
      deallocate(elemmap)
      deallocate(speciesname_wk)
      deallocate(iatomn_wk)
      deallocate(amion_wk)

100   continue

      if(ipriinputfile >= 1 .and. printable) &
           & write(nfout,'(" !** -- << set_element_detail >> --")')
      icount = 0
      if(.not.allocated(surface_integral_method_paw)) then
         allocate(surface_integral_method_paw(ntyp))
         surface_integral_method_paw = SphericalHarmonicsExpansion
      endif
      do i = 1, ntyp
         if (i == 1) then
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if
         ip = i
         icount = icount + 1
         iret = f_getIntValue(tag_id,rint)
         if(iret /= 0 ) iret = f_getIntValue(tag_no,rint)
         if(iret == 0 .and. rint > 0 .and. rint <= ntyp) ip = rint
         if( f_getRealValue(tag_atomicnumber, dret, '') == 0) iatomn(ip) = dret

         if(printable) write(nfout,'(" !** iatomn(",i6,") = ",f8.4)') ip,iatomn(ip)

         if( f_getRealValue(tag_mass, dret, 'au_mass') == 0) amion(ip) = dret

! ==== KT_mod ==== 13.3B
!         if( amion(ip) < DELTA) then
         if( dret < DELTA ) then
! ================ 13.3B
            if(printable) write(nfout,'(" !** amion(",i6,") = ",d20.8)') ip,amion(ip)
            call phase_error_with_msg(nfout,' amion is too small <<m_IS_rd_n.set_elelment_detail>>',__LINE__,__FILE__)
         end if

#ifndef _EMPIRICAL_
         ! --- zeta ---
         if( f_getRealValue(tag_zeta, dret, '') == 0) then
           zeta1(ip) = dret
           mag_zeta1_atomtyp_is_defined = .true.
         endif
!!$         if( f_getRealValue(tag_variance, dret, 'bohr') ==0 ) alfa(ip) = 1.d0/(dret*dret)
         ! --- parameter for the initial_charge_density ---
         number_is_given = f_getRealValue(tag_deviation, dret, 'bohr') == 0
         if(.not.number_is_given) &
              number_is_given = f_getRealValue(tag_standard_deviation, dret, 'bohr') == 0
         if(.not.number_is_given) number_is_given = f_getRealValue(tag_dev, dret, 'bohr') == 0
         if(number_is_given) then
!!$            alfa(ip) = 1.d0/(dret*dret)
            alfa(ip) = 1.d0/(2*dret*dret)
            if(ipriinputfile >= 1 .and. printable) &
                 & write(nfout,'(" !** deviation(",i2, &
                 & ") of the Gauss. distrib. func. for the initial charge construction = ",f10.5)') ip,dret
         end if
         ! -- excess charge ---
         if( f_getRealValue(tag_qex,dret,'') == 0) qex(ip) = dret

! ====================== KT_add ======================== 11.0&13.0U
!!
!!!       set initial values of ionic_charge and magnetic moment
!
         call set_initial_ioncharge_atomtyp(ip)
         call set_initial_magmom_atomtyp(ip)

         if ( noncol ) then
            if ( .not. mag_moment0_atomtyp_is_defined ) then
               call set_initial_orientation_magmom(ip)
            endif
            call set_elemtype_wrt_lcore_filling(ip)
         endif
! ====================================================== 11.0&13.0U

! ==== KT_add ==== 11.0 & 13.2S
!-----------------------------------------------------------------------
!!
!!!       set scaling factor for Spin-orbit
!!
!----------------------------------------------------------------------
         if ( noncol .or. sw_spinorbit_second_variation == ON &
                     .or. sw_calc_core_energy == ON ) then
            if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
                 &                          .or. SpinOrbit_Mode == ReadFromPP ) then
               if( f_getRealValue( tag_scaling_so, dret, '' ) == 0) then
                  scaling_so(ip) = dret
               endif
            endif
         endif
! =============== 11.0 & 13.2S

#endif
         if( f_getStringValue(tag_element, rstr, NOCONV) == 0) then
            len_rstrtrimmed = len_trim(rstr)
            if(len_rstrtrimmed > LEN_ATOMNAME) then
               if(ipriinputfile >= 1 .and. printable) &
                    & write(nfout,'(" !* element name is larger than ",i6)') LEN_ATOMNAME
               n_cp_name = LEN_ATOMNAME
            else
               n_cp_name = len_rstrtrimmed
            end if
            speciesname(ip)(1:n_cp_name) = rstr(1:n_cp_name)
         end if
         if( f_getStringValue(tag_surface_integral_method, rstr, LOWER) == 0)then
            if(rstr == tag_sph) surface_integral_method_paw(ip) = SphericalHarmonicsExpansion
            if(rstr == tag_gl)  surface_integral_method_paw(ip) = GaussLegendre
         endif
         if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** ityp = ",i6," : iatomn,amion,zeta1,alfa,qex,type = " &
              &              , f8.4,f12.4,3f8.4,3x,a4)') &
              & i,iatomn(ip),amion(ip),zeta1(ip),alfa(ip),qex(ip),speciesname(ip)
         if(f_getIntValue(tag_proj_group, rint)==0) then
           if(rint>0) then
             if(.not.allocated(iproj_group_tmp)) then
               allocate(iproj_group_tmp(ntyp));iproj_group_tmp = 0
             endif
             iproj_group_tmp(i) = rint
           endif
         endif
      end do

! ================== KT_add =============== 11.0&13.0U
      call message_on_chgmagmom_atomtyp

      if ( noncol ) then
         call print_message_on_magmom_direc
         call print_message_on_lcore_filling
         if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
              &                          .or. SpinOrbit_Mode == ReadFromPP ) then
            call print_message_on_scaling_so
         endif
      else
         if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
              &                          .or. SpinOrbit_Mode == ReadFromPP ) then
            call print_message_on_scaling_so
         endif
      endif
! ========================================= 11.0&13.0U

! === KT_add === 13.3B
#ifdef MASS_AUTOMATIC_SET
      call set_amion_default( ntyp, iatomn, amion )
#endif
! ============== 13.3B

    end subroutine set_element_detail

    subroutine build_elemmap(ntyp,speciesname,speciesname_wk,elemmap)
      integer, intent(in) :: ntyp
      character(len=LEN_ATOMNAME), dimension(ntyp),intent(in) :: speciesname,speciesname_wk
      integer, dimension(ntyp), intent(out) :: elemmap
      integer :: i,j
      logical :: logi
      integer :: lens
      do i=1,ntyp
         lens = len(trim(speciesname(i)))
         jloop:do j=1,ntyp
             logi = trim(speciesname(i)(1:lens)) .eq.  trim(speciesname_wk(j)(1:lens))
             if(trim(speciesname(i)(1:lens)) == trim(speciesname_wk(j)(1:lens))) then
                 elemmap(i) = j
                 exit jloop
             endif
         enddo jloop
      enddo
    end subroutine build_elemmap

    subroutine set_vdw_parameters()
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
      use m_Const_Parameters, only : VDW_WILLIAMS, VDW_GRIMME
      use m_Control_Parameters, only : vdw_method, vdw_scaling_factor, vdw_scaling_factor_r
! ==============================================================================
      integer :: i, j, iret, rint, ip, ip1, ip2, len_rstrtrimmed, n_cp_name
      integer :: f_selectFirstTableLine, f_selectNextTableLine, f_getRealValue
      real(kind=DP) :: dret

      if(ipriinputfile >= 1 .and. printable) &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
           & write(nfout,'(" !** -- << set_default_vdw_parameters >> --")')
      if(ipriinputfile >= 1 .and. printable) &
           & write(nfout,'(" !** -- elements --")')

      do i = 1, ntyp_vdw
	select case(vdw_method)
        case(VDW_WILLIAMS)
          select case(speciesname_vdw(i))
	  case('H')
	    c6vdw(i) = 2.831179918
	    r0vdw(i) = 1.17
	    pvdw(i) = 0.387
	  case('F')
	    c6vdw(i) = 3.94987377
	    pvdw(i) = 0.296
	  case('Cl')
	    c6vdw(i) = 3.94987377
	    pvdw(i) = 2.315
	  case('Br')
	    c6vdw(i) = 128.2756865
	    pvdw(i) = 3.013
	  case('I')
	    c6vdw(i) = 309.0603852
	    pvdw(i) = 5.415
	  case('CTE')
	    c6vdw(i) = 22.67403316
	    r0vdw(i) = 1.70
	    pvdw(i) = 1.061
	  case('CTR')
	    c6vdw(i) = 32.61525204
	    r0vdw(i) = 1.70
	    pvdw(i) = 1.352
	  case('CAR')
	    c6vdw(i) = 49.790/vdw_scaling_factor
	    r0vdw(i) = 1.70
	    pvdw(i) = 1.352
	  case('CBR')
	    c6vdw(i) = 54.16430826
	    r0vdw(i) = 1.70
	    pvdw(i) = 1.896
	  case('CDI')
	    c6vdw(i) = 30.15058105
	    r0vdw(i) = 1.70
	    pvdw(i) = 1.283
	  case('NTE')
	    c6vdw(i) = 20.89758657
	    r0vdw(i) = 1.50
	    pvdw(i) = 0.964
	  case('NTR2')
	    c6vdw(i) = 23.08003267
	    r0vdw(i) = 1.50
	    pvdw(i) = 1.030
	  case('NPI2')
	    c6vdw(i) = 25.12582491
	    r0vdw(i) = 1.50
	    pvdw(i) = 1.090
	  case('NDI')
	    c6vdw(i) = 20.63799109
	    r0vdw(i) = 1.50
	    pvdw(i) = 0.956
	  case('OTE')
	    c6vdw(i) = 11.86370812
	    r0vdw(i) = 1.40
	    pvdw(i) = 0.637
	  case('OTR4')
	    c6vdw(i) = 10.01566303
	    r0vdw(i) = 1.40
	    pvdw(i) = 0.569
	  case('OPI2')
	    c6vdw(i) = 3.346856941
	    r0vdw(i) = 1.40
	    pvdw(i) = 0.274
	  case('STE')
	    c6vdw(i) = 121.2531939
	    r0vdw(i) = 1.80
	    pvdw(i) = 3.000
	  case('STR4')
	    c6vdw(i) = 168.0350502
	    r0vdw(i) = 1.80
	    pvdw(i) = 3.729
	  case('SPI2')
	    c6vdw(i) = 103.5277919
	    r0vdw(i) = 1.80
	    pvdw(i) = 2.700
	  case('PTE')
	    c6vdw(i) = 42.11289383
	    r0vdw(i) = 1.80
	    pvdw(i) = 1.538
          end select

	case(VDW_GRIMME)
          select case(speciesname_vdw(i))
	  case('H')
	    c6vdw(i) = 0.14
	    r0vdw(i) = 1.001
	  case('He')
	    c6vdw(i) = 0.08
	    r0vdw(i) = 1.012
	  case('Li')
	    c6vdw(i) = 1.61
	    r0vdw(i) = 0.825
	  case('Be')
	    c6vdw(i) = 1.61
	    r0vdw(i) = 1.408
	  case('B')
	    c6vdw(i) = 3.13
	    r0vdw(i) = 1.485
	  case('C')
	    c6vdw(i) = 1.75
	    r0vdw(i) = 1.452
	  case('N')
	    c6vdw(i) = 1.23
	    r0vdw(i) = 1.397
	  case('O')
	    c6vdw(i) = 0.70
	    r0vdw(i) = 1.342
	  case('F')
	    c6vdw(i) = 0.75
	    r0vdw(i) = 1.287
	  case('Ne')
	    c6vdw(i) = 0.63
	    r0vdw(i) = 1.243
	  case('Na')
	    c6vdw(i) = 5.71
	    r0vdw(i) = 1.144
	  case('Mg')
	    c6vdw(i) = 5.71
	    r0vdw(i) = 1.364
	  case('Al')
	    c6vdw(i) = 10.79
	    r0vdw(i) = 1.716
	  case('Si')
	    c6vdw(i) = 9.23
	    r0vdw(i) = 1.716
	  case('P')
	    c6vdw(i) = 7.84
	    r0vdw(i) = 1.705
	  case('S')
	    c6vdw(i) = 5.57
	    r0vdw(i) = 1.683
	  case('Cl')
	    c6vdw(i) = 5.07
	    r0vdw(i) = 1.639
	  case('Ar')
	    c6vdw(i) = 4.61
	    r0vdw(i) = 1.595
	  case('K')
	    c6vdw(i) = 10.80
	    r0vdw(i) = 1.485
	  case('Ca')
	    c6vdw(i) = 10.80
	    r0vdw(i) = 1.474
	  case('Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn')
	    c6vdw(i) = 10.80
	    r0vdw(i) = 1.562
	  case('Ga')
	    c6vdw(i) = 16.99
	    r0vdw(i) = 1.650
	  case('Ge')
	    c6vdw(i) = 17.10
	    r0vdw(i) = 1.727
	  case('As')
	    c6vdw(i) = 16.37
	    r0vdw(i) = 1.760
	  case('Se')
	    c6vdw(i) = 12.64
	    r0vdw(i) = 1.771
	  case('Br')
	    c6vdw(i) = 12.47
	    r0vdw(i) = 1.749
	  case('Kr')
	    c6vdw(i) = 12.01
	    r0vdw(i) = 1.727
	  case('Rb')
	    c6vdw(i) = 24.67
	    r0vdw(i) = 1.628
	  case('Sr')
	    c6vdw(i) = 24.67
	    r0vdw(i) = 1.606
	  case('Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd')
	    c6vdw(i) = 24.67
	    r0vdw(i) = 1.639
	  case('In')
	    c6vdw(i) = 37.32
	    r0vdw(i) = 1.672
	  case('Sn')
	    c6vdw(i) = 38.71
	    r0vdw(i) = 1.804
	  case('Sb')
	    c6vdw(i) = 38.44
	    r0vdw(i) = 1.881
	  case('Te')
	    c6vdw(i) = 31.74
	    r0vdw(i) = 1.892
	  case('I')
	    c6vdw(i) = 31.50
	    r0vdw(i) = 1.892
	  case('Xe')
	    c6vdw(i) = 29.99
	    r0vdw(i) = 1.881
          end select
	end select

        if(ipriinputfile >= 1 .and. printable) &
              & write(nfout,'(" !** ",i6,1x," : c6vdw,r0vdw = " &
              &              , 3f12.4,3x,a4)') &
              & i,c6vdw(i),r0vdw(i),pvdw(i),speciesname_vdw(i)
      end do

      if(ipriinputfile >= 1 .and. printable) &
! ==============================================================================
           & write(nfout,'(" !** -- << set_vdw_parameters >> --")')
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
      if(ipriinputfile >= 1 .and. printable) &
           & write(nfout,'(" !** -- elements --")')

      if( f_selectBlock( tag_vdw_list) == 0) then

! ==============================================================================
      do i = 1, ntyp_vdw
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!        do j = 1, i
!        if (i == 1 .and. j == 1) then
         if (i == 1) then
! ==============================================================================
            if(f_selectFirstTableLine() /= 0) then
               exit
            end if
         else
            if(f_selectNextTableLine() /= 0) then
               exit
            end if
         end if

         ip1 = 0
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!        if( f_getStringValue(tag_type1, rstr, NOCONV) == 0) then
         if( f_getStringValue(tag_type, rstr, NOCONV) == 0) then
! ==============================================================================
            len_rstrtrimmed = len_trim(rstr)
            if(len_rstrtrimmed > LEN_ATOMNAME) then
               if(ipriinputfile >= 1 .and. printable) &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!                   & write(nfout,'(" !* vdW-type1 name is larger than ",i6)') LEN_ATOMNAME
                    & write(nfout,'(" !* vdW-type name is larger than ",i6)') LEN_ATOMNAME
! ==============================================================================
               n_cp_name = LEN_ATOMNAME
            else
               n_cp_name = len_rstrtrimmed
            end if
            do ip=1,ntyp_vdw
               if(speciesname_vdw(ip)(1:n_cp_name) == rstr(1:n_cp_name)) then
                  ip1 = ip
                  exit
               end if
            end do
            if(ip1 == 0) then
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!              if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* vdW-type1 name did not match.")')
!              stop 'vdw-type1 name did not match.'
!           end if
!        end if

!        ip2 = 0
!        if( f_getStringValue(tag_type2, rstr, NOCONV) == 0) then
!           len_rstrtrimmed = len_trim(rstr)
!           if(len_rstrtrimmed > LEN_ATOMNAME) then
!              if(ipriinputfile >= 1 .and. printable) &
!                   & write(nfout,'(" !* vdW-type2 name is larger than ",i6)') LEN_ATOMNAME
!              n_cp_name = LEN_ATOMNAME
!           else
!              n_cp_name = len_rstrtrimmed
!           end if
!           do ip=1,ntyp_vdw
!              if(speciesname_vdw(ip)(1:n_cp_name) == rstr(1:n_cp_name)) then
!                 ip2 = ip
!                 exit
!              end if
!           end do
!           if(ip2 == 0) then
!              if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* vdW-type2 name did not match.")')
!              stop 'vdw-type2 name did not match.'
               if(ipriinputfile >= 1 .and. printable) write(nfout,'(" !* vdW-type name did not match.")')
               call phase_error_with_msg(nfout,'vdw-type name did not match.',__LINE__,__FILE__)
! ==============================================================================
            end if
         end if

         if( f_getRealValue(tag_c6, dret, '') == 0) then
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!            cvdw(ip1,ip2) = dret
!            cvdw(ip2,ip1) = dret
             c6vdw(ip1) = dret
! ==============================================================================
         end if
         if( f_getRealValue(tag_r0, dret, '') == 0) then
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!            rvdw(ip1,ip2) = dret
!            rvdw(ip2,ip1) = dret
             r0vdw(ip1) = dret
         end if
         if( f_getRealValue(tag_p, dret, '') == 0) then
             pvdw(ip1) = dret
! ==============================================================================
         end if

         if(ipriinputfile >= 1 .and. printable) &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
              & write(nfout,'(" !** ",i6,1x," : c6vdw,r0vdw = " &
              &              , 3f12.4,3x,a4)') &
              & ip1,c6vdw(ip1),r0vdw(ip1),pvdw(ip1),speciesname_vdw(ip1)
      end do

      iret = f_selectParentBlock()
      end if

      if(ipriinputfile >= 1 .and. printable) &
           & write(nfout,'(" !** -- atom pair --")')

      do i = 1, ntyp_vdw
        do j = 1, i
          select case(vdw_method)
          case(VDW_WILLIAMS)
            cvdw(i,j) = -vdw_scaling_factor * (2.0d0*c6vdw(i)*c6vdw(j)*pvdw(i)*pvdw(j)) &
                       / (pvdw(i)**2*c6vdw(i)+pvdw(j)**2*c6vdw(j))
            rvdw(i,j) = vdw_scaling_factor_r * ((2.0d0*r0vdw(i))**3+(2.0d0*r0vdw(j))**3) &
                       / ((2.0d0*r0vdw(i))**2+(2.0d0*r0vdw(j))**2)
          case(VDW_GRIMME)
            cvdw(i,j) = dsqrt(c6vdw(i)*c6vdw(j))
            rvdw(i,j) = r0vdw(i)+r0vdw(j)
            cvdw(i,j) = cvdw(i,j) * 3.8088e-7 * 1.0e6 / BOHR**6  ! J/mol nm6 -> hartree A6
          end select

          rvdw(i,j) = rvdw(i,j) / BOHR  ! ang -> bohr
          cvdw(j,i) = cvdw(i,j)
          rvdw(j,i) = rvdw(i,j)

          if(ipriinputfile >= 1 .and. printable) &
! ==============================================================================
              & write(nfout,'(" !** ",i6,1x,i6," : cvdw,rvdw = " &
              &              , 2f12.4,3x,a4,1x,a4)') &
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!             & ip1,ip2,cvdw(ip1,ip2),rvdw(ip1,ip2),speciesname_vdw(ip1),speciesname_vdw(ip2)
!        end do
              & i,j,cvdw(i,j),rvdw(i,j),speciesname_vdw(i),speciesname_vdw(j)
        end do
! ==============================================================================
      end do

    end subroutine set_vdw_parameters

! ============================= added by K. Tagami ======================== 11.0
#ifndef _EMPIRICAL_
    subroutine set_initial_orientation_magmom( ip )
      integer, intent(in) :: ip
      integer :: mode = 0

      real(kind=DP) :: cnorm, mdx, mdy, mdz
      real(kind=DP) :: theta, phi

      mdx = 0.0d0; mdy = 0.0d0; mdz = 0.0d0

      if( f_getRealValue( tag_mdx, dret, '') == 0 ) then
         mdx = dret;  mode = 1
      endif
      if( f_getRealValue( tag_mdy, dret, '') == 0 ) then
         mdy = dret;  mode = 1
      endif
      if( f_getRealValue( tag_mdz, dret, '') == 0 ) then
         mdz = dret;  mode = 1
      endif

      cnorm = sqrt( mdx**2 + mdy**2 + mdz**2 )

      theta = 0.0d0; phi = 0.0d0
      if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
         theta = dret;  mode = 2
      endif
      if ( f_getRealValue( tag_phi, dret, "" ) == 0 )   then
         phi = dret;  mode = 2
      endif

      theta  = theta / 180.0d0 *PAI
      phi  =  phi / 180.0d0 *PAI

      if ( mode == 1 ) then
         if ( abs(cnorm) > 1.0E-4 ) then
            mag_direction0_atomtyp(ip,1) = mdx / cnorm;
            mag_direction0_atomtyp(ip,2) = mdy / cnorm;
            mag_direction0_atomtyp(ip,3) = mdz / cnorm;
         else
            mag_direction0_atomtyp(ip,1) = 0.0d0
            mag_direction0_atomtyp(ip,2) = 0.0d0
            mag_direction0_atomtyp(ip,3) = 1.0d0
         endif
      else if ( mode == 2 ) then
         mag_direction0_atomtyp(ip,1) = sin( theta ) *cos( phi )
         mag_direction0_atomtyp(ip,2) = sin( theta ) *sin( phi )
         mag_direction0_atomtyp(ip,3) = cos( theta )
      endif

    end subroutine set_initial_orientation_magmom

#endif

    subroutine print_message_on_magmom_direc
      integer :: it

      write(nfout,*) '! ------------ Initial Magnetic Orientation --- '
      write(nfout,*) '     id     name        mdx         mdy         mdz  '
      Do it=1, ntyp
         write(nfout,'(I8,3X,A7,3F12.6)') it, speciesname(it), &
              &                         mag_direction0_atomtyp(it,1), &
              &                         mag_direction0_atomtyp(it,2), &
              &                         mag_direction0_atomtyp(it,3)
      End Do
      write(nfout,*) '! ---------------------------------------------'
    end subroutine print_message_on_magmom_direc

    subroutine set_elemtype_wrt_lcore_filling(ip)
      integer, intent(in) :: ip
      integer :: iret

      if( f_getIntValue( tag_lcore_parfil,iret ) == 0 ) then
         has_partially_filled_lcore(ip) = iret
      endif

    end subroutine set_elemtype_wrt_lcore_filling

    subroutine print_message_on_lcore_filling
      integer :: it

      write(nfout,*) '! ----- Info. of localized core orbitals ( such as 4f )--------'
      write(nfout,*) '     id     name    atom_has_partially_filled_lcore_orb (Yes:1)'
      Do it=1, ntyp
         write(nfout,'(I8,3X,A7,5X,I3)') it, speciesname(it), &
              &                             has_partially_filled_lcore(it)
      End Do
      write(nfout,*) '! ---------------------------------------------'

    end subroutine print_message_on_lcore_filling

    subroutine print_message_on_scaling_so
      integer :: it

      write(nfout,*) '! ----- Info. of scaling param. of spin-orbit strength '
      write(nfout,*) '     id     name    scaling_so'
      Do it=1, ntyp
         write(nfout,'(I8,3X,A7,X,F10.4)') it, speciesname(it), scaling_so(it)
      End Do
      write(nfout,*) '! ---------------------------------------------'

    end subroutine print_message_on_scaling_so
! ========================================================================= 11.0

! ===================== KT_add ================================== 13.0U &13.2S
    subroutine set_initial_magmom_atomtyp( ip )
      integer, intent(in) :: ip
      integer :: mode = 0

      real(kind=DP) :: cnorm, mx, my, mz, theta, phi
      real(kind=DP) :: mdx, mdy, mdz, c1

      cnorm = 0.0d0

      if ( noncol ) then

         if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
            cnorm = dret;  mode = 2
         else
            mx = 0.0d0; my = 0.0d0; mz = 0.0d0
            if( f_getRealValue( tag_mx, dret, '') == 0 ) then
               mx = dret;  mode = 1
            endif
            if( f_getRealValue( tag_my, dret, '') == 0 ) then
               my = dret;  mode = 1
            endif
            if( f_getRealValue( tag_mz, dret, '') == 0 ) then
               mz = dret;  mode = 1
            endif
            cnorm = sqrt( mx**2 + my**2 + mz**2 )
         endif

         if ( mag_zeta1_atomtyp_is_defined ) then
            cnorm = abs( zeta1(ip) );  mode = 2
         endif

         if ( mode == 2 ) then
            mdx = 0.0d0;  mdy = 0.0d0;  mdz = 0.0d0
            theta = 0.0d0; phi = 0.0d0

            if ( f_getRealValue( tag_theta, dret, "" ) == 0 ) then
               theta = dret

               if ( f_getRealValue( tag_phi, dret, "" ) == 0 )   phi = dret
               theta = theta /180.0d0 *PAI;  phi = phi /180.0d0 *PAI

               mdx = sin( theta ) *cos( phi )
               mdy = sin( theta ) *sin( phi )
               mdz = cos( theta )

            else
               if ( f_getRealValue( tag_mdx, dret, "" ) == 0 ) mdx = dret
               if ( f_getRealValue( tag_mdy, dret, "" ) == 0 ) mdy = dret
               if ( f_getRealValue( tag_mdz, dret, "" ) == 0 ) mdz = dret

               c1 = sqrt( mdx**2 +mdy**2 + mdz**2 )
               if ( c1 > 1.0D-4 ) then
                  mdx = mdx /c1
                  mdy = mdy /c1
                  mdz = mdz /c1
               else
                  mdx = 0.0d0;  mdy = 0.0d0;  mdz = 1.0d0
               endif
            endif
         endif

         if ( mode == 1 ) then
            mag_moment0_atomtyp(ip,1) = mx
            mag_moment0_atomtyp(ip,2) = my
            mag_moment0_atomtyp(ip,3) = mz
         else if ( mode == 2 ) then
            mag_moment0_atomtyp(ip,1) = cnorm *mdx
            mag_moment0_atomtyp(ip,2) = cnorm *mdy
            mag_moment0_atomtyp(ip,3) = cnorm *mdz
         endif

         if ( sw_fix_global_quantz_axis == ON ) then
            mx = mag_moment0_atomtyp(ip,1)
            my = mag_moment0_atomtyp(ip,2)
            mz = mag_moment0_atomtyp(ip,3)
            cnorm = sqrt( mx**2 +my**2 +mz**2 )

            c1 = mag_moment0_atomtyp(ip,1) * Global_Quantz_Axis_Fixed(1) &
                 & +mag_moment0_atomtyp(ip,2) * Global_Quantz_Axis_Fixed(2) &
                 & +mag_moment0_atomtyp(ip,3) * Global_Quantz_Axis_Fixed(3)
            if ( c1 > 0.0 ) then
               mag_moment0_atomtyp(ip,1:3) =  cnorm *Global_Quantz_Axis_Fixed(1:3)
            else
               mag_moment0_atomtyp(ip,1:3) = -cnorm *Global_Quantz_Axis_Fixed(1:3)
            endif
         endif

      else
         if ( f_getRealValue( tag_moment, dret, "" ) == 0 ) then
            cnorm = dret;  mode = 1
         endif
         mag_moment0_atomtyp(ip,1) = cnorm
      endif

      if ( mode /= 0 ) then
         mag_moment0_atomtyp_is_defined = .true.
         if ( noncol ) then
            mx = mag_moment0_atomtyp(ip,1)
            my = mag_moment0_atomtyp(ip,2)
            mz = mag_moment0_atomtyp(ip,3)
            cnorm = sqrt( mx**2 +my**2 +mz**2 )
            if ( cnorm > 1.0D-4 ) then
               mag_direction0_atomtyp(ip,1) = mx /cnorm
               mag_direction0_atomtyp(ip,2) = my /cnorm
               mag_direction0_atomtyp(ip,3) = mz /cnorm
            else
               mag_direction0_atomtyp(ip,1) = 0.0d0
               mag_direction0_atomtyp(ip,2) = 0.0d0
               mag_direction0_atomtyp(ip,3) = 1.0d0
            endif
         endif
      endif
    end subroutine set_initial_magmom_atomtyp

    subroutine set_initial_ioncharge_atomtyp( ip )
      integer, intent(in) :: ip
      real(kind=DP) :: dret

      if ( f_getRealValue( tag_ion_charge, dret, "" ) == 0 ) then
         ionic_charge_atomtyp(ip) = dret
      endif

    end subroutine set_initial_ioncharge_atomtyp

    subroutine message_on_chgmagmom_atomtyp
      integer :: it

      if ( noncol ) then
         write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atomtyp) --- '
         write(nfout,*) '     id     name       ion         mx          my          mz'
         Do it=1, ntyp
            write(nfout,'(I8,3X,A7,4F12.6)') it, speciesname(it), &
                 &                         ionic_charge_atomtyp(it), &
                 &                         mag_moment0_atomtyp(it,1), &
                 &                         mag_moment0_atomtyp(it,2), &
                 &                         mag_moment0_atomtyp(it,3)
         End Do
      else if ( sw_fix_global_quantz_axis == ON ) then
         write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atomtyp) --- '
         write(nfout,*) '     id     name       ion         mx          my          mz'
         Do it=1, ntyp
            write(nfout,'(I8,3X,A7,4F12.6)') it, speciesname(it), &
                 &         ionic_charge_atomtyp(it), &
                 &         mag_moment0_atomtyp(it,1)*Global_Quantz_Axis_Fixed(1), &
                 &         mag_moment0_atomtyp(it,1)*Global_Quantz_Axis_Fixed(2), &
                 &         mag_moment0_atomtyp(it,1)*Global_Quantz_Axis_Fixed(3)
         End Do
      else
         write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atomtyp) --- '
         write(nfout,*) '     id     name       ion        moment '
         Do it=1, ntyp
            write(nfout,'(I8,3X,A7,2F12.6)') it, speciesname(it), &
                 &                         ionic_charge_atomtyp(it), &
                 &                         mag_moment0_atomtyp(it,1)
         End Do
      endif
      write(nfout,*) '! ---------------------------------------------'

    end subroutine message_on_chgmagmom_atomtyp
! =============================================================== 13.0U & 13.2S

    subroutine set_initial_velocities(imode)
      integer, intent(in) :: imode
      integer :: i,j,ia,ir, irp, ib
      integer,dimension(natm) :: imdt
      integer :: icnstrnt_typ
      real(kind=DP) :: a,b,p
      real(kind=DP) :: xn
      real(kind=DP),dimension(nrsv,natm,3) :: random
      real(kind=DP)   :: sumt
      real(kind=DP)   :: mcom
      real(kind=DP),dimension(3)   :: pcom
      real(kind=DP),dimension(nrsv)   :: tkin
      integer, dimension(nrsv) :: nir
      integer nrsv_atom, imdalg_t
      integer :: rand_seed = 9
      real(kind=DP) :: sig
      real(kind=DP), parameter :: eps = 1.d-12

      if(conf_para)then
         rand_seed=rand_seed+mype_conf
      endif

      if (iseed>0) then
        rand_seed = iseed
        if(printable.and.ipri>=1) write(nfout,'(a,i8)') ' !** set iseed for initial velocity : ',rand_seed
      endif

      if(imode == 1) then
         do i=1,natm
            if ( imdtyp(i) .le. NOSE_HOOVER ) then
               imdt(i) = NOSE_HOOVER + 1
            else
               imdt(i) = imdtyp(i)
            endif
         enddo
         imdalg_t = T_CONTROL
      else
         do i=1, natm
            imdt(i) = imdtyp(i)
         end do
         imdalg_t = VERLET
      end if

      if(sw_read_velocities == OFF) then

         xn = 0.d0
         a  = 32771.d0
         b  = 1234567891.d0 + dble(rand_seed)*2
         p  = 214783648.d0
         sumt = 0.d0

         ! create normal random numbers.
         do i=1,3
            do ia=1,natm
               ir = icnstrnt_typ(imdt(ia),imdalg_t)
               if ( ir >= 1 ) then
                  do j=1,12
                     xn = dble(mod(xn*a+b,p))
                     sumt = sumt+xn/p
                  enddo
                  random(ir,ia,i) = (sumt-6.0d0)/6.0d0
                  sumt = random(ir,ia,i)
               endif
            enddo
         enddo

         ! randomize velocity
         cpd_l = 0.d0
         do i=1,3
            pcom(i) = 0.d0
            do ia=1,natm
               ir = icnstrnt_typ(imdt(ia),imdalg_t)
               if ((imode == 1 .and. ir >= 1 .and. imdtyp(ia).ne.0).or. &
             &     (imode == 2 .and. ir == 1 .and. imdtyp(ia).ne.0)) then
                  sig = sign(random(ir,ia,i),1.d0)
                  cpd_l(ia,i) = sqrt(abs(random(ir,ia,i))/amion(ityp(ia))) * sig
                  if(iprimd >= 2) write(nfout,'(" !!! ia, ir = ",2i8," cpd_l(ia,",i3,") = " &
                       & ,f12.6, " <<set_initial_velocities>>")') ia, ir, i, cpd_l(ia,i)
                  pcom(i) = pcom(i) + amion(ityp(ia))*cpd_l(ia,i)
               else
                  if(iprimd >= 2) write(nfout,'(" !!! ia, ir = ",2i8," cpd_l(ia,",i3,") = " &
                       & ,f12.6, " * <<set_initial_velocities>>")') ia, ir, i, cpd_l(ia,i)
               endif
            enddo
         enddo

      else ! sw_read_velocities == ON

         ! sum of momenta
         do i=1,3
            pcom(i) = 0.d0
            do ia=1,natm
               ir = icnstrnt_typ(imdt(ia),imdalg_t)
               if (  ir >= 1 ) then
                   pcom(i) = pcom(i) + amion(ityp(ia))*cpd_l(ia,i)
               endif
            enddo
         enddo

      end if

      mcom = 0.d0
      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir >= 1 .and. imdtyp(ia).ne.0 ) then
               mcom = mcom + amion(ityp(ia))
            end if
         end do
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir == 1 .and. imdtyp(ia).ne.0 ) then
               mcom = mcom + amion(ityp(ia))
            end if
         end do
      end if

      ! shift velocity
      if(mcom.gt.eps) pcom(:) = pcom(:)/mcom

      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (   ir >= 1 .and. imdtyp(ia).ne.0 ) then
               cpd_l(ia,:) = cpd_l(ia,:) - pcom(:)
            endif
         enddo
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (   ir == 1 .and. imdtyp(ia).ne.0 ) then
               cpd_l(ia,:) = cpd_l(ia,:) - pcom(:)
            endif
         enddo
      end if

      ! scale velocity
      nir = 0
      tkin = 0.d0
      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir >= 1 .and. imdtyp(ia).ne.0 )  then
               irp = ir
               if(irp > nrsv) irp = 1
               tkin(irp) = tkin(irp)+dot_product(cpd_l(ia,1:3),cpd_l(ia,1:3)) &
                    &               * amion(ityp(ia))*iwei(ia)*0.5d0
               nir(irp) = nir(irp) + 3*iwei(ia)
            endif
         enddo
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir == 1 .and. imdtyp(ia).ne.0 )  then
               tkin(ir) = tkin(ir)+dot_product(cpd_l(ia,1:3),cpd_l(ia,1:3)) &
                    &               * amion(ityp(ia))*iwei(ia)*0.5d0
               nir(ir) = nir(ir) + 3*iwei(ia)
            endif
         enddo
      end if

      if(sw_fix_bond==ON) then
        do ir=1,nrsv
         nir(ir) = nir(ir) - nbonds_per_thermo(ir)
        enddo
!        if(imode == 1) then
!          do ib=1,nfixed_bonds
!            ir = icnstrnt_typ(imdt(fixed_bond(ib)%atoms(1)),imdalg_t)
!            if (  ir >= 1 .and. imdtyp(fixed_bond(ib)%atoms(1)).ne.0 )  then
!               irp = ir
!               if(irp > nrsv) irp = 1
!               nir(irp) = nir(irp) - iwei(fixed_bond(ib)%atoms(1))
!            endif
!          enddo
!        else
!          do ib=1,nfixed_bonds
!            ir = icnstrnt_typ(imdt(fixed_bond(ib)%atoms(1)),imdalg_t)
!            if (  ir == 1 .and. imdtyp(fixed_bond(ib)%atoms(1)).ne.0 )  then
!               irp = ir
!               if(irp > nrsv) irp = 1
!               nir(irp) = nir(irp) - iwei(fixed_bond(ib)%atoms(1))
!            endif
!          enddo
!        endif
      endif

      if(iprimd >= 1) then
         do ir = 1, nrsv
            if ( natm_per_thermo(ir) == 0 ) cycle
            if ( tkin(ir).lt.1e-12) cycle
            write(nfout,'(" !! tkin(",i3,")  = ",d20.8)') ir,tkin(ir)
            write(nfout,'(" !! nir( ",i3,")  = ",i8)') ir,nir(ir)
            write(nfout,'(" !! tkb( ",i3,")  = ",d20.8)') ir,tkb(ir)
            write(nfout,'(" !! dsqrt(0.5d0*nir(ir)*tkb(ir)/tkin(ir)) = ",d20.8)') dsqrt(0.5d0*nir(ir)*tkb(ir)/tkin(ir))
         end do
      end if

      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir >= 1 .and. imdtyp(ia).ne.0 ) then
               irp = ir
               if(irp > nrsv) irp = 1
               if(tkin(irp).gt.1e-12) &
             & cpd_l(ia,:) = cpd_l(ia,:) * dsqrt(0.5d0*nir(irp)*tkb(irp)/tkin(irp))
            endif
         enddo
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir == 1 .and. imdtyp(ia).ne.0 ) then
               if(tkin(ir).gt.1e-12) &
             & cpd_l(ia,:) = cpd_l(ia,:) * dsqrt(0.5d0*nir(ir)*tkb(ir)/tkin(ir))
            endif
         enddo
      end if

      if(iprimd >= 1) then
         write(nfout,*) 'initial velocities (at m_Ionic_System.set_initial_velocities)'
         write(nfout,'("Translational velocity = ",3(1x,e18.9))') pcom(1:3)
         write(nfout,'(" --- initial velocities ---")')
         do ia=1,natm
            write(nfout,'(i5,3(1x,e18.9))') ia,cpd_l(ia,1:3)
         end do
      end if

      if (nrigid_bodies>1) call set_initial_velocities_rb()

    end subroutine set_initial_velocities

    subroutine import_from_dynm_pre(nf,F_IDEN,chara,ntyp,natm)
      integer, intent(in) :: nf
      character(len=*), intent(in) :: F_IDEN
      character(len=*), intent(in) :: chara
      integer, intent(out) :: ntyp,natm
      integer :: natm_tmp,ntyp_tmp,ierr
      call m_Files_open(nf,F_IDEN,chara)
      if(mype == 0)then
        call get_numat_and_ntyp(nf,natm_tmp,ntyp_tmp)
      endif
      call m_Files_close(nf)
      call mpi_bcast(natm_tmp,1,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(ntyp_tmp,1,mpi_integer,0,MPI_CommGroup,ierr)
      natm = natm_tmp;ntyp = ntyp_tmp
    end subroutine import_from_dynm_pre

    subroutine import_from_dynm(nf,F_IDEN,chara,frame,natm,ntyp,cps,pos,cpd_l,ityp,speciesname,altv,rltv)
      integer, intent(in) :: nf
      character(len=*), intent(in) :: F_IDEN
      character(len=*), intent(in) :: chara
      integer, intent(in) :: frame,natm,ntyp
      real(kind=DP), intent(out), dimension(natm,3) :: cps,pos,cpd_l
      integer, intent(out), dimension(natm) :: ityp
      character(len=LEN_ATOMNAME), intent(out), dimension(ntyp) :: speciesname
      real(kind=DP), dimension(3,3), intent(out) :: altv,rltv
      integer :: natm_tmp,ntyp_tmp,ierr
      real(kind=DP), dimension(3) :: avec,bvec,cvec
      real(kind=DP), dimension(3,3) :: rltv_t
      real(kind=DP) :: univol,rvol
      if(printable) then
         write(nfout,'(a)') ' !** importing atomic coordinates from '//trim(F_IDEN)
      endif
      call m_Files_open(nf,F_IDEN,chara)
      if(mype == 0)then
      call get_ac_from_dynm(nf,frame,natm,4,ntyp,avec,bvec,cvec,speciesname,ityp,cps,cpd_l)
      endif
      call m_Files_close(nf)
      call mpi_bcast(avec,3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(bvec,3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call mpi_bcast(cvec,3,mpi_double_precision,0,MPI_CommGroup,ierr)
      altv(1,1)=avec(1);altv(2,1)=avec(2);altv(3,1)=avec(3)
      altv(1,2)=bvec(1);altv(2,2)=bvec(2);altv(3,2)=bvec(3)
      altv(1,3)=cvec(1);altv(2,3)=cvec(2);altv(3,3)=cvec(3)
      call mpi_bcast(speciesname,ntyp,mpi_character,0,MPI_CommGroup,ierr)
      call mpi_bcast(ityp,natm,mpi_integer,0,MPI_CommGroup,ierr)
      call mpi_bcast(cps,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
      call altv_2_rltv(altv,rltv,univol,rvol)
      rltv_t = transpose(rltv)/PAI2
      call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.)
      call mpi_bcast(cpd_l,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
    end subroutine import_from_dynm

    subroutine initial_species_name()
      integer :: ia,j,icount
      logical :: found
      integer :: lene
      icount=0
      do ia=1,natm
         found = .false.
         do j=1,icount
             if(trim(species_work(ia)) == trim(speciesname(j))) then
                 found = .true.
             endif
         enddo
         if(.not.found)then
             icount = icount+1
             lene = len(trim(species_work(ia)))
             speciesname(icount)(1:lene) = species_work(ia)(1:lene)
         endif
      enddo
    end subroutine initial_species_name

  end subroutine m_IS_rd_n

  subroutine set_mobility_by_distance(target_atom,target_pos,distance)
    integer, intent(in) :: target_atom
    real(kind=DP), dimension(3), intent(in) :: target_pos
    real(kind=DP), intent(in) :: distance
    real(kind=DP) :: distance2,r2
    integer       :: neibrd, alen(3)
    real(kind=DP), allocatable :: rxyz(:,:)    ! d(neibrd,3)
    real(kind=DP), allocatable :: rr(:)        ! d(neibrd)
    real(kind=DP), pointer, dimension(:,:) :: cps_fp
    integer,       pointer, dimension(:)   :: ityp_full ! d(natm2)
    integer :: nu,in,ia
    real(kind=DP), dimension(3) :: dr,tpos
    character(len=10) :: smob
    distance2 = distance*distance
    call decide_rxyz_size(distance, alen, neibrd)
    allocate(rxyz(neibrd,3))
    allocate(rr(neibrd))
    call substitute_rxyz(alen, neibrd, rxyz, rr)
    allocate(cps_fp(natm2,3))
    allocate(ityp_full(natm2))
    cps_fp(1:natm,1:3) = pos(1:natm,1:3)
    ityp_full(1:natm)  = ityp(1:natm)
    call rplcps(cps_fp,ityp_full,1,natm2,natm,iwei)
    call cpspac ! -> cps_fp
    imdtyp=OFF;imdtypxyz=OFF
    if(target_atom>0) then
      tpos(:) = cps_fp(target_atom,:)
    else
      tpos = target_pos
    endif
    do in = 1, neibrd
       do ia=1,natm
          dr(1:3) = tpos(1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
          r2 = dot_product(dr,dr)
          if (r2 <= distance2) then
            imdtyp(ia) = ON
            imdtypxyz(ia,1) = ON
            imdtypxyz(ia,2) = ON
            imdtypxyz(ia,3) = ON
          endif
       enddo
    enddo
    if(printable)then
      write(nfout,'(a)') ' !** mobility of atoms'
      write(nfout,'(8i8)') (imdtyp(ia), ia=1,natm)
      if(ipriinputfile>=2)then
        write(nfout,'(i8)') natm
        write(nfout,'(a)')
        do ia=1,natm
          write(smob,'(i1)') imdtyp(ia)
          write(nfout,'(a,3f15.5)') &
          & trim(adjustl(speciesname(ityp_full(ia))))//trim(adjustl(smob)),&
          & BOHR*cps_fp(ia,1),BOHR*cps_fp(ia,2),BOHR*cps_fp(ia,3)
        enddo
      endif
    endif

    deallocate(rxyz)
    deallocate(rr)
    deallocate(cps_fp)
    deallocate(ityp_full)
    contains
    subroutine cpspac
      real(kind=DP), dimension(3) :: catoms(3)
      integer                     :: i
      do i = 1, natm2
         catoms = cps_fp(i,1:3)
         catoms = catoms - nint(catoms)      !Packing
         cps_fp(i,1:3) = matmul(altv,catoms) !Change of coordinate system
      end do
    end subroutine cpspac
  end subroutine set_mobility_by_distance

  subroutine shift_velocities(imode)
      integer, intent(in) :: imode
      integer :: i,j,ia,ir
      real(kind=DP)   :: mcom
      integer,dimension(natm) :: imdt
      integer :: icnstrnt_typ
      integer :: imdalg_t
      real(kind=DP),dimension(3)   :: pcom
      real(kind=DP), parameter :: eps = 1.d-12
      if(imode == 1) then
         do i=1,natm
            if ( imdtyp(i) .le. NOSE_HOOVER ) then
               imdt(i) = NOSE_HOOVER + 1
            else
               imdt(i) = imdtyp(i)
            endif
         enddo
         imdalg_t = T_CONTROL
      else
         do i=1, natm
            imdt(i) = imdtyp(i)
         end do
         imdalg_t = VERLET
      end if

      do i=1,3
         pcom(i) = 0.d0
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if ((imode == 1 .and. ir >= 1 .and. imdtyp(ia).ne.0).or. &
          &     (imode == 2 .and. ir == 1 .and. imdtyp(ia).ne.0)) then
              pcom(i) = pcom(i) + amion(ityp(ia))*cpd_l(ia,i)
            endif
         enddo
      enddo

      mcom = 0.d0
      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir >= 1 .and. imdtyp(ia).ne.0 ) then
               mcom = mcom + amion(ityp(ia))
            end if
         end do
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (  ir == 1 .and. imdtyp(ia).ne.0 ) then
               mcom = mcom + amion(ityp(ia))
            end if
         end do
      end if

      ! shift velocity
      if(mcom.gt.eps) pcom(:) = pcom(:)/mcom

      if(imode == 1) then
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (   ir >= 1 .and. imdtyp(ia).ne.0 ) then
               cpd_l(ia,:) = cpd_l(ia,:) - pcom(:)
            endif
         enddo
      else
         do ia=1,natm
            ir = icnstrnt_typ(imdt(ia),imdalg_t)
            if (   ir == 1 .and. imdtyp(ia).ne.0 ) then
               cpd_l(ia,:) = cpd_l(ia,:) - pcom(:)
            endif
         enddo
      end if
  end subroutine shift_velocities

!!$ 2011.06.06
  subroutine scale_velocity()
     integer :: ia,ir,irp
     integer, dimension(nrsv) :: nir
     real(kind=DP),dimension(nrsv)   :: tkin
     integer,dimension(natm) :: imdt
     integer :: icnstrnt_typ

     if(sw_fix_bond==ON) call fixed_bond_velocities(cpd_l)

     do ia=1,natm
        if ( imdtyp(ia) .le. NOSE_HOOVER ) then
           imdt(ia) = NOSE_HOOVER + 1
        else
           imdt(ia) = imdtyp(ia)
        endif
     enddo
     nir = 0
     tkin = 0.d0
     do ia=1,natm
        ir = icnstrnt_typ(imdt(ia),T_CONTROL)
        if (  ir >= 1 .and. imdtyp(ia).ne.0 )  then
           irp = ir
           if(irp > nrsv) irp = 1
           tkin(irp) = tkin(irp)+dot_product(cpd_l(ia,1:3),cpd_l(ia,1:3)) &
                &               * amion(ityp(ia))*iwei(ia)*0.5d0
           nir(irp) = nir(irp) + 3 * iwei(ia)
        endif
     enddo

     if(sw_fix_bond==ON) then
       do ir=1,nrsv
         nir(ir) = nir(ir) - nbonds_per_thermo(ir)
       enddo
     endif

     do ia=1,natm
        ir = icnstrnt_typ(imdt(ia),T_CONTROL)
        if (  ir >= 1 .and. imdtyp(ia).ne.0 )  then
           irp = ir
           if(irp > nrsv) irp = 1
           if(tkin(irp).gt.1e-12) &
           & cpd_l(ia,:) = cpd_l(ia,:) * dsqrt(0.5d0*nir(irp)*tkb(irp)/tkin(irp))
        endif
     enddo

     if(iprimd >= 1) then
        write(nfout,'(a)') 'scaled the velocities'
        do ir = 1, nrsv
!           if ( natm_per_thermo(ir) == 0 ) cycle
           if(tkin(ir)<1e-14) cycle
           write(nfout,'(" !! tkin(",i3,")  = ",d20.8)') ir,tkin(ir)
           write(nfout,'(" !! nir( ",i3,")  = ",i8)') ir,nir(ir)
           write(nfout,'(" !! tkb( ",i3,")  = ",d20.8)') ir,tkb(ir)
           write(nfout,'(" !! dsqrt(0.5d0*nir(ir)*tkb(ir)/tkin(ir)) = ",d20.8)') dsqrt(0.5d0*(nir(ir))*tkb(ir)/tkin(ir))
        end do
     end if
  end subroutine scale_velocity

  subroutine scale_velocity_rb()
    integer  :: i,ip
    real(DP) :: kin,factor
    integer, dimension(nrsv) :: nir
    real(kind=DP),dimension(nrsv)   :: tkin
    tkin=0.d0;nir=0
    do i=1,nrigid_bodies
      ip = rigid_bodies(i)%thermo_group
      tkin(ip) = tkin(ip) + rb_kinetic_energy(rigid_bodies(i))
      nir(ip) = nir(ip)+1
    enddo
    do i=1,nrigid_bodies
      ip                                    = rigid_bodies(i)%thermo_group
      factor                                = sqrt(3.d0*nir(ip)*tkb(ip)/tkin(ip))
      rigid_bodies(i)%velocity(:)           = rigid_bodies(i)%velocity(:)*factor
      rigid_bodies(i)%angular_momentum(:)   = rigid_bodies(i)%angular_momentum(:)*factor
      rigid_bodies(i)%velocity_old(:)       = rigid_bodies(i)%velocity(:)
      rigid_bodies(i)%angular_momentum_h(:) = rigid_bodies(i)%angular_momentum(:)
    enddo
    tkin=0.d0
    do i=1,nrigid_bodies
      ip       = rigid_bodies(i)%thermo_group
      tkin(ip) = tkin(ip)+rb_kinetic_energy(rigid_bodies(i))
    enddo
  end subroutine scale_velocity_rb

!!$ 2011.06.06

!!$  subroutine alloc_fcvect_tmp()
!!$    allocate(fcvect_tmp(num_planes_atoms_are_fixed,4))
!!$  end subroutine alloc_fcvect_tmp
!!$
!!$  subroutine dealloc_fcvect_tmp()
!!$    deallocate(fcvect_tmp)
!!$  end subroutine dealloc_fcvect_tmp

  subroutine m_IS_alloc_iatomn_etc
    allocate(iatomn(ntyp));iatomn=0.d0
    allocate(iatom(ntyp))
!!$    allocate(iloc_inputf(ntyp))
    allocate(ivan(ntyp)); ivan = 1
    allocate(alfa(ntyp)); alfa = 0.15

! ==== KT_mod ==== 13.3B
!    allocate(amion(ntyp)); amion = 51577.50
!    allocate(amion(ntyp)); amion = 51196.421251715d0 ! mass of Si
!
#ifdef MASS_AUTOMATIC_SET
    allocate(amion(ntyp)); amion = -100.0d0
#else
    allocate(amion(ntyp)); amion = 51196.421251715d0 ! mass of Si
#endif
! ================= 13.3B

    allocate(zeta1(ntyp)); zeta1 = 0.0
    allocate(qex(ntyp)); qex = 0.d0

! ===================================== added by K. Tagami ========== 11.0
    if ( noncol ) then
       allocate( mag_direction0_atomtyp(ntyp,3) )
       mag_direction0_atomtyp(:,1) = 0.0d0
       mag_direction0_atomtyp(:,2) = 0.0d0
       mag_direction0_atomtyp(:,3) = 1.0d0
!
       allocate( has_partially_filled_lcore(ntyp) )
       has_partially_filled_lcore = 0
    endif
! =================================================================== 11.0
! ===================================== added by K. Tagami ========== 11.0 & 13.2S
    if ( noncol .or. sw_spinorbit_second_variation == ON &
         &      .or. sw_calc_core_energy == ON ) then
       if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
            &                          .or. SpinOrbit_Mode == ReadFromPP ) then
          allocate( scaling_so(ntyp) )
          scaling_so(:) = 1.0d0
       endif
    endif
! =================================================================== 11.0 & 13.2S

! ====================== KT_add ================ 13.0U
    if ( noncol ) then
       allocate( mag_moment0_atomtyp(ntyp,3) )
    else
       allocate( mag_moment0_atomtyp(ntyp,1) )
    endif
    mag_moment0_atomtyp = 0.0d0

    allocate( ionic_charge_atomtyp( ntyp) ); ionic_charge_atomtyp = 0.0d0
! ============================================== 13.0U

      if(.not.allocated(surface_integral_method_paw)) then
         allocate(surface_integral_method_paw(ntyp))
         surface_integral_method_paw = SphericalHarmonicsExpansion
      endif

  end subroutine m_IS_alloc_iatomn_etc

  subroutine m_IS_dealloc_iatomn_etc
    if(allocated(iatomn)) deallocate(iatomn)
    if(allocated(iatom)) deallocate(iatom)
    if(allocated(ivan)) deallocate(ivan)
    if(allocated(alfa)) deallocate(alfa)
    if(allocated(amion)) deallocate(amion)
    if(allocated(zeta1)) deallocate(zeta1)
    if(allocated(qex)) deallocate(qex)

! ===== KT_ add =========== 2013/10/30
    if ( noncol ) then
       if (allocated(mag_direction0_atomtyp))     deallocate(mag_direction0_atomtyp)
       if (allocated(has_partially_filled_lcore)) deallocate(has_partially_filled_lcore)
    endif
    if (allocated(scaling_so)) deallocate(scaling_so)
! ========================= 2013/10/30

! ================= KT_add ================ 13.0U
    if ( allocated(mag_moment0_atomtyp) ) deallocate( mag_moment0_atomtyp )
    if ( allocated(ionic_charge_atomtyp) ) deallocate( ionic_charge_atomtyp )
! ========================================= 13.0U

    if (allocated(surface_integral_method_paw)) deallocate(surface_integral_method_paw)

  end subroutine m_IS_dealloc_iatomn_etc

! ===================================== added by K. Tagami ========== 11.0
  subroutine m_IS_alloc_magmom_local
    !    allocate( magmom_local_now(ista_atm:iend_atm,3) )
    if ( allocated(magmom_local_now) ) deallocate( magmom_local_now )
    allocate( magmom_local_now(1:natm,3) )
    magmom_local_now = 0.0d0
  end subroutine m_IS_alloc_magmom_local
! ============================================================ 11.0

! ====== KT_add ==================== 2013/10/31
  subroutine m_IS_dealloc_magmom_local
    deallocate( magmom_local_now )
  end subroutine m_IS_dealloc_magmom_local
! ================================== 2013/10/31

! ======================== added by K. Tagami ===================== 11.0 & 13.2S
  subroutine m_IS_init_magmom_local
    magmom_local_now = 0.0d0

    if ( noncol ) then
       call case_noncol
    else if ( sw_use_magnetic_symmetry == ON ) then
       call case_col
    endif

  contains

    subroutine case_noncol
      integer :: ia, it

      Do ia=1, natm
         it = ityp(ia)
         magmom_local_now(ia,:) = mag_moment0_atomtyp(it,:)
      End Do
! == KT_add == 2014/12/29
      if ( mag_moment0_atoms_is_defined ) then
         Do ia=1, natm
            it = ityp(ia)
            magmom_local_now(ia,:) = mag_moment0_atoms(ia,:)
         End Do
      endif
! ============ 2014/12/29
    end subroutine case_noncol

    subroutine case_col
      integer :: ia, it

      Do ia=1, natm
         it = ityp(ia)
         magmom_local_now(ia,:) = mag_moment0_atomtyp(it,1)*Global_Quantz_Axis_Fixed(:)
      End Do
! == KT_add == 2014/12/29
      if ( mag_moment0_atoms_is_defined ) then
         Do ia=1, natm
            it = ityp(ia)
            magmom_local_now(ia,:) = mag_moment0_atoms(ia,1)*Global_Quantz_Axis_Fixed(:)
         End Do
      endif
! ============ 2014/12/29
    end subroutine case_col

  end subroutine m_IS_init_magmom_local
! ============================================================ 11.0 & 13.2S

  subroutine m_IS_alloc_vdw
    allocate(cvdw(ntyp_vdw,ntyp_vdw))
    allocate(rvdw(ntyp_vdw,ntyp_vdw))
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
    allocate(c6vdw(ntyp_vdw))
    allocate(r0vdw(ntyp_vdw))
    allocate(pvdw(ntyp_vdw))
! ==============================================================================
    if(.not.allocated(ityp_vdw)) allocate(ityp_vdw(natm)); ityp_vdw = 1
    if(.not.allocated(fxyzvdw_l)) allocate(fxyzvdw_l(natm,3))
#ifdef __EDA__
    if ( sw_eda == ON ) then
       if (.not.allocated(evdw_on_atom) ) allocate( evdw_on_atom(natm) )
    endif
#endif
  end subroutine m_IS_alloc_vdw

  subroutine m_IS_alloc_vdwdf3
    allocate(fxyzvdw_l(natm,3))
#ifdef __EDA__
    if ( sw_eda == ON ) then
       allocate( evdw_on_atom(natm) )
    endif
#endif
  end subroutine m_IS_alloc_vdwdf3

! =========================== KT_add ======================= 13.0B
  subroutine m_IS_dealloc_vdw
    if ( allocated(cvdw) ) deallocate(cvdw)
    if ( allocated(rvdw) ) deallocate(rvdw)
    if ( allocated(c6vdw) ) deallocate(c6vdw)
    if ( allocated(r0vdw) ) deallocate(r0vdw)
    if ( allocated(pvdw) )  deallocate(pvdw)
    if ( allocated(ityp_vdw) ) deallocate(ityp_vdw)
    if ( allocated(fxyzvdw_l) ) deallocate(fxyzvdw_l)
#ifdef __EDA__
    if ( allocated(evdw_on_atom) ) deallocate(evdw_on_atom)
#endif
  end subroutine m_IS_dealloc_vdw
! ========================================================= 13.0B

  subroutine m_IS_alloc_napt
    if(printable) then
       write(nfout,'(" -- allocation of napt --")')
       write(nfout,'(" !! natm = ",i5," nopr+af = ",i5)') natm,nopr+af
       if(ngen_tl > 0) write(nfout,'(" !! natm = ",i5," nopr+af+ngen_tl = ",i5)') natm,nopr+af+ngen_tl
    end if
    if(allocated(napt)) deallocate(napt)
    allocate(napt(natm,nopr+af))
    if(nopr == 0) call phase_error_with_msg(nfout,' nopr == 0',__LINE__,__FILE__)
    if(ngen_tl > 0) then
       if(allocated(napt_tl)) deallocate(napt_tl)
       allocate(napt_tl(natm,ngen_tl))
    else
       if(allocated(napt_tl)) deallocate(napt_tl)
       allocate(napt_tl(natm,1))
    end if
  end subroutine m_IS_alloc_napt

  subroutine m_IS_alloc_fxyzew
    if(allocated(fxyzew_l)) deallocate(fxyzew_l)
    allocate(fxyzew_l(natm,3))
  end subroutine m_IS_alloc_fxyzew

  subroutine m_IS_alloc_zfm3_3D
    if(.not.allocated(zfm3_l)) allocate(zfm3_l(ista_kngp:iend_kngp,ntyp,kimg)); zfm3_l = 0.d0
  end subroutine m_IS_alloc_zfm3_3D

  subroutine m_IS_dealloc_zfm3_3D
    if(allocated(zfm3_l)) deallocate(zfm3_l)
  end subroutine m_IS_dealloc_zfm3_3D

  subroutine m_IS_gdiis_alloc(mode_init,if_allocated,nfree)
    integer, intent(in) ::  mode_init
    integer, intent(out) :: if_allocated
    integer, intent(in), optional :: nfree
    integer :: nf
    integer ::             i,j
    real(kind=DP) ::       fac
    nf = natm
    if(present(nfree)) nf = nfree
    if(kqnmditer_p <= 0) call phase_error_with_msg(nfout,' kqnmditer_p is illegal (m_IS_gdiis_alloc)',__LINE__,__FILE__)
    if(printable) write(nfout,'(" !! kqnmditer_p = ",i6," <<m_IS_gdiis_alloc>>")') kqnmditer_p
    if(.not.allocated(u_l)) allocate(u_l(nf,3,kqnmditer_p));u_l=0.d0
    if(.not.allocated(w_l)) allocate(w_l(nf,3,kqnmditer_p));w_l=0.d0
    if(.not.allocated(fc_l)) allocate(fc_l(nf,3))
    do i = 1, nf
       do j=1,3
          if(imdtypxyz(i,j) == 0) then
             fac = 0.d0
          else
             if(mode_init == UNIT) then
                fac = 1.d0
             else
                fac = dtio*dtio/amion(ityp(i))
             end if
          end if
          fc_l(i,j) = -fac
       enddo
    end do

    if(.not.allocated(ncrspd)) allocate(ncrspd(kqnmditer_p))
    ncrspd(:) = (/(i,i=1,kqnmditer_p)/)
    if(.not.allocated(f_gdiis)) allocate(f_gdiis(kqnmditer_p,kqnmditer_p))
    if(.not.allocated(g)) allocate(g(kqnmditer_p))
    if(.not.allocated(f_wk)) allocate(f_wk(kqnmditer_p*kqnmditer_p,2))
    if(.not.allocated(f_rslv)) allocate(f_rslv(kqnmditer_p*kqnmditer_p))
    if(.not.allocated(e_wk)) allocate(e_wk(kqnmditer_p*kqnmditer_p))
    if(.not.allocated(ww1)) allocate(ww1(kqnmditer_p))
    if(.not.allocated(etot_trial)) allocate(etot_trial(0:2))
    if(.not.allocated(forc_g)) allocate(forc_g(nf,3))
    if(.not.allocated(ip)) allocate(ip(kqnmditer_p))
    if_allocated = 1
    if(icond==CONTINUATION .or. icond==AUTOMATIC .and. diis_continuable)then
       if(allocated(u_l_buf)) u_l(:,:,:) = u_l_buf(:,:,:)
       if(allocated(w_l_buf)) w_l(:,:,:) = w_l_buf(:,:,:)
       if(allocated(ncrspd_buf)) ncrspd(:) = ncrspd_buf(:)
    endif
    if(allocated(u_l_buf))    deallocate(u_l_buf)
    if(allocated(w_l_buf))    deallocate(w_l_buf)
    if(allocated(ncrspd_buf)) deallocate(ncrspd_buf)
  end subroutine m_IS_gdiis_alloc

  subroutine m_IS_gdiis_reset()
     iter_gdiis = 0
     if(if_allocated==0) return
     call m_IS_gdiis_dealloc(if_allocated)
  end subroutine m_IS_gdiis_reset

  subroutine m_IS_put_iteration_ionic_in_constraint(iteration_put)
    integer, intent(in) :: iteration_put
    iteration_ionic_in_constraint = iteration_put
  end subroutine m_IS_put_iteration_ionic_in_constraint

  logical function m_IS_iter_ionic_is_over()
    m_IS_iter_ionic_is_over = .false.
    if(iteration_ionic_in_constraint  >= 20) m_IS_iter_ionic_is_over = .true.
    if(iprimd >= DEBUGPRINTLEVEL ) then
       write(nfout,'(" iteration_ionic_in_constraint = ",i8)') iteration_ionic_in_constraint
       call flush(nfout)
    end if
  end function m_IS_iter_ionic_is_over

  subroutine m_IS_freeze()
    integer :: i
    do i=1,natm
       imdtyp(i) = 0
       imdtypxyz(i,1) = 0
       imdtypxyz(i,2) = 0
       imdtypxyz(i,3) = 0
    enddo
  end subroutine m_IS_freeze

  subroutine m_IS_gdiis_dealloc(if_allocated)
    integer, intent(out) :: if_allocated
    deallocate(ip);  deallocate(forc_g);   deallocate(etot_trial)
    deallocate(ww1);   deallocate(e_wk);    deallocate(f_rslv)
    deallocate(f_wk);   deallocate(g);     deallocate(f_gdiis)
    deallocate(ncrspd);   deallocate(fc_l);    deallocate(w_l)
    deallocate(u_l)
    if_allocated = 0
  end subroutine m_IS_gdiis_dealloc

  subroutine m_IS_wd_speciesname_etc(nfdynm)
    use m_IterationNumbers, only : iteration_stress_correction
    integer, intent(in) :: nfdynm
    integer :: i
    if(sw_stress_correction==ON .and. iteration_stress_correction<=2) return
    if(mype == 0) then
       write(nfdynm,'("#")')
       write(nfdynm,'("#   a_vector = ",3f20.10)') altv(1:3,1)
       write(nfdynm,'("#   b_vector = ",3f20.10)') altv(1:3,2)
       write(nfdynm,'("#   c_vector = ",3f20.10)') altv(1:3,3)
       write(nfdynm,'("#   ntyp = ",i8, " natm = ",i8)') ntyp, natm
       write(nfdynm,'("# (natm->type) ",10i5)') (ityp(i),i=1,natm)
       do i = 1, ntyp
          write(nfdynm,'("# (speciesname) ",i5," :   ", a4)') i,speciesname(i)
       end do
       write(nfdynm,'("#")')
    end if
  end subroutine m_IS_wd_speciesname_etc

  subroutine m_IS_rd_moved_distances_of_planes(nfcntn)
    integer, intent(in) :: nfcntn
    logical             :: tag_is_found, EOF_reach
    real(kind=DP) :: d
    integer :: i
    if(mype==0)then
       call rewind_to_tag0(nfcntn,len(tag_distances_of_planes),tag_distances_of_planes, EOF_reach, tag_is_found, str, len_str)
       if(.not.tag_is_found) then
          if(iprimd>=1) write(nfout,'(" tag_distances_of_planes is not found")')
       else
          read(nfcntn,*) num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb
          if(num_planes_atoms_are_fixed_cog>=1) then
             do i = 1, num_planes_atoms_are_fixed_cog
                read(nfcntn,*) d
                distance_cog(i) = d
             end do
          end if
          if(num_planes_atoms_are_fixed_rb>=1) then
             do i = 1, num_planes_atoms_are_fixed_rb
                read(nfcntn,*) d
                distance_rb(i) = d
             end do
          end if
       end if
    end if
    if(npes>1)then
       call mpi_bcast(num_planes_atoms_are_fixed_cog,1,mpi_integer,0, MPI_CommGroup,ierr)
       call mpi_bcast(num_planes_atoms_are_fixed_rb, 1,mpi_integer,0, MPI_CommGroup,ierr)
       if(num_planes_atoms_are_fixed_cog>=1) &
            & call mpi_bcast(distance_cog,num_planes_atoms_are_fixed_cog, mpi_double_precision,0,MPI_CommGroup,ierr)
       if(num_planes_atoms_are_fixed_rb>=1) &
            & call mpi_bcast(distance_rb, num_planes_atoms_are_fixed_rb,  mpi_double_precision,0,MPI_CommGroup,ierr)
    end if
  end subroutine m_IS_rd_moved_distances_of_planes

  subroutine m_IS_wd_moved_distances_of_planes(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i
    if(num_planes_atoms_are_fixed_cog+num_planes_atoms_are_fixed_rb >=1) then
       if(mype == 0) then
          write(nfcntn,*) tag_distances_of_planes
          write(nfcntn,'(2i4," : num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb")') &
               & num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb
          if(num_planes_atoms_are_fixed_cog>=1) then
             do i = 1, num_planes_atoms_are_fixed_cog
                write(nfcntn,'(d16.8," : distance_cog(",i4,")")') distance_cog(i),i
             end do
          end if
          if(num_planes_atoms_are_fixed_rb>=1) then
             do i = 1, num_planes_atoms_are_fixed_rb
                write(nfcntn,'(d16.8," : distance_rb(",i4,")")') distance_rb(i),i
             end do
          end if
       end if
    end if
  end subroutine m_IS_wd_moved_distances_of_planes

  subroutine m_IS_rd_pos_and_v(nfcntn)
    integer, intent(in) :: nfcntn
    integer             :: i, k, natm_t
    logical             :: tag_is_found, EOF_reach
    if(mype==0)then
       call rewind_to_tag0(nfcntn,len(tag_ionic_system),tag_ionic_system &
            &, EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(nfout,' tag_ionic_system is not found',__LINE__,__FILE__)
       else
          read(nfcntn,*)
          read(nfcntn,*) natm_t
       endif
    endif
    if(npes>1)then
       call mpi_bcast(natm_t,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
    endif
    if(natm_t.ne.natm)then
       if(printable)then
          write(nfout,'(a)') ' !** natm_t .ne. natm'
       endif
       natmorg = natm
       natm = natm_t
       call m_IS_dealloc_pos_and_v(nfout)
       call m_IS_alloc_pos_and_v(nfout)
       if(mype==0)then
          call rewind_to_tag0(nfcntn,len(tag_ionic_system_attributes),tag_ionic_system_attributes &
               &, EOF_reach, tag_is_found, str,len_str)
          read(nfcntn,*) &
         & (iwei(i),imdtyp(i),ityp(i),if_pdos(i),if_aldos(i),ihubbard(i), &
         & iproj_group(i),numlay(i),imdtypxyz(i,1:3),i=1,natm)
          natm2 = 0
          do i=1,natm
             natm2 = natm2+iwei(i)
          enddo
       endif
       if(npes>1)then
          call mpi_bcast(natm2,1,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(iwei,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(imdtyp,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(ityp,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(if_pdos,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(if_aldos,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(ihubbard,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(iproj_group,natm,mpi_integer,0,MPI_CommGroup,ierr)
          call mpi_bcast(numlay,natm,mpi_integer,0,MPI_CommGroup,ierr)
       endif
       if(printable)then
          write(nfout,'(a,i8,a,i8)') ' !** natm : ',natm,' natm2 ',natm2
       endif
    endif

    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_ionic_system),tag_ionic_system &
            &, EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(nfout,' tag_ionic_system is not found',__LINE__,__FILE__)
       else
          read(nfcntn,*)
          read(nfcntn,*) natm_t
          read(nfcntn,*)
          read(nfcntn,*) (pos(i,1),pos(i,2),pos(i,3),i=1,natm_t)
          read(nfcntn,*)
          read(nfcntn,*) (cps(i,1),cps(i,2),cps(i,3),i=1,natm_t)
          read(nfcntn,*)
          read(nfcntn,*) (cpd_l(i,1),cpd_l(i,2),cpd_l(i,3),i=1,natm_t)
          do k = 1, 3
             read(nfcntn,*)
             read(nfcntn,*) (cpo_l(i,1,k),cpo_l(i,2,k),cpo_l(i,3,k),i=1,natm_t)
          end do
       end if
       if(imdalg == QUENCHED_CONSTRAINT) then
          call rewind_to_tag0(nfcntn,len(tag_forcmx_const),tag_forcmx_const &
               &,  EOF_reach, tag_is_found, str, len_str)
          if(.not.tag_is_found) then
             forcmx_constraint_quench= forcmx_constraint_quench-1
          else
             read(nfcntn,*) forcmx_constraint_quench
          end if
       end if
       if(sw_optimize_lattice==ON .or. &
          imdalg == PT_CONTROL    .or. &
          imdalg == P_CONTROL) then
          call rewind_to_tag0(nfcntn,len(tag_lattice_vector),tag_lattice_vector &
               &,  EOF_reach, tag_is_found, str, len_str)
          if(tag_is_found)then
             read(nfcntn,*) altv(1,1),altv(2,1),altv(3,1)
             read(nfcntn,*) altv(1,2),altv(2,2),altv(3,2)
             read(nfcntn,*) altv(1,3),altv(2,3),altv(3,3)
          endif
       endif
    end if

    if(npes > 1) then
       call mpi_bcast(pos,natm*3 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(cps,natm*3 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(cpd_l,natm*3 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(cpo_l,natm*3*3 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(forcmx_constraint_quench,1 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(sw_optimize_lattice==ON .or. &
          imdalg == PT_CONTROL    .or. &
          imdalg == P_CONTROL) then
          call mpi_bcast(altv,9,mpi_double_precision,0,MPI_CommGroup,ierr)
       endif
    end if

    if(sw_optimize_lattice==ON .or. &
       imdalg == PT_CONTROL    .or. &
       imdalg == P_CONTROL) then
       call m_CS_altv_2_rltv()
    endif

    call m_IS_set_iatom(nfout)

    if(ipriinputfile >= 1 .and. printable) call wd_atom_list()
    if(nrigid_bodies>0) call m_IS_rb_reinitialize()

  contains
    subroutine wd_atom_list()
      integer :: i
      write(nfout,'(" wd_atom_list from m_IS_rd_pos_and_v(nfcntn)")')
      if(ipriinputfile >= 3) then
         write(nfout,'(" !** === Atomic coordinates expressed in the internal system ===")')
         write(nfout,'(" !** id,  rx,    ry,    rz,    weight,  imdtyp, ityp,  species")')
         do i = 1, natm
            write(nfout,210) i,pos(i,1),pos(i,2),pos(i,3),iwei(i),imdtyp(i),ityp(i),speciesname(ityp(i))
         end do
210      format(' !** ',i5,3f18.10,i3,i6,i3,3x,a4)
         write(nfout,'(" !** === Atomic coordinates expressed in the cartesian system ===")')
         write(nfout,'(" !** id,  rx,    ry,    rz,    weight,  imdtyp, ityp,  species")')
         do i = 1, natm
            write(nfout,210) i,cps(i,1),cps(i,2),cps(i,3),iwei(i),imdtyp(i),ityp(i),speciesname(ityp(i))
         end do
      else
         write(nfout,'(" !** === Atomic coordinates ==")')
         write(nfout,'(" !**   id  ( coordinates_in_Intrnal_sys  ) (  coordinates_in_Cartsian_system  ) weight    ityp")')
         write(nfout,'(" !**       (   rx         ry         rz  ) (    rx          ry          rz    )      imdtyp    species")')
         do i = 1, natm
            write(nfout,211) i, pos(i,1:3), cps(i,1:3), iwei(i), imdtyp(i), ityp(i), speciesname(ityp(i))
         end do
211      format(' !** ',i5,3f11.6,3f12.4,i3,i6,i3,3x,a4)
      end if

!!$      do i = 1, ntyp
!!$         write(nfout,'(" !** i = ",i5," element_name = ",a4)') i, speciesname(i)
!!$      end do

    end subroutine wd_atom_list

  end subroutine m_IS_rd_pos_and_v

  subroutine m_IS_wd_pos_and_v(nfcntn)
    integer, intent(in) :: nfcntn
    integer             :: i, k
    if(iprimd >= 2) write(nfout,'(" tag_ionic_system")')
    if(mype==0) then
       write(nfcntn,*) tag_ionic_system
       write(nfcntn,'("  (natm)",/,i10)') natm
       write(nfcntn,'("  (pos)")')
       write(nfcntn,'(3d24.16)') (pos(i,1),pos(i,2),pos(i,3),i=1,natm)
       write(nfcntn,'("  (cps)")')
       write(nfcntn,'(3d24.16)') (cps(i,1),cps(i,2),cps(i,3),i=1,natm)
       write(nfcntn,'("  (cpd)")')
       write(nfcntn,'(3d24.16)') (cpd_l(i,1),cpd_l(i,2),cpd_l(i,3),i=1,natm)
       do k = 1, 3
          write(nfcntn,'("  (cpo(",i3,"))")') k
          write(nfcntn,'(3d24.16)') &
               & (cpo_l(i,1,k),cpo_l(i,2,k),cpo_l(i,3,k),i=1,natm)
       end do
       write(nfcntn,*) tag_forcmx_const
       write(nfcntn,'(d24.16)') forcmx_constraint_quench
       write(nfcntn,*) tag_ionic_system_attributes
       write(nfcntn,'(11i8)') &
      & (iwei(i),imdtyp(i),ityp(i),if_pdos(i),if_aldos(i),ihubbard(i), &
      & iproj_group(i),numlay(i),imdtypxyz(i,1:3),i=1,natm)
       if(sw_optimize_lattice==ON)then
          write(nfcntn,*) tag_lattice_vector
          write(nfcntn,'(3d24.16)') altv(1,1),altv(2,1),altv(3,1)
          write(nfcntn,'(3d24.16)') altv(1,2),altv(2,2),altv(3,2)
          write(nfcntn,'(3d24.16)') altv(1,3),altv(2,3),altv(3,3)
       endif
    endif

  end subroutine m_IS_wd_pos_and_v

  subroutine alloc_endpoint_pos(nfout)
    integer, intent(in) :: nfout
    allocate(pos_end0(natm,3))
    allocate(pos_end1(natm,3))
    allocate(cps_end0(natm,3))
    allocate(cps_end1(natm,3))
    if(printable) write(nfout,'(" !** pos_end0, pos_end1, cps_end0, and cps_end1 are allocated "&
         & ," <<alloc_endpoint_pos>>")')
  end subroutine alloc_endpoint_pos

  subroutine m_IS_alloc_pos_and_v(nfout)
    integer, intent(in) :: nfout
    allocate(iwei(natm)); iwei = 1
    allocate(imdtyp(natm)); imdtyp = 0
    allocate(imdtypxyz(natm,3));imdtypxyz = 0
    allocate(ityp(natm)); ityp = 1
#ifndef _EMPIRICAL_
    allocate(if_pdos(natm)); if_pdos = 1
    allocate(if_aldos(natm)); if_aldos = 1
    allocate(ihubbard(natm)); ihubbard = 0
    allocate(iproj_group(natm)); iproj_group = 0

! ============================ added by K. Tagami ================== 11.0
    allocate(itab_spinorbit_addition(natm))
    itab_spinorbit_addition = 0
! ================================================================== 11.0
#endif
    allocate(ionic_mass(natm)); ionic_mass = -1.d0
    allocate(pos(natm,3))
    allocate(cps(natm,3))
!    if(.not.allocated(cpd_l)) then
    if(allocated(cpd_l)) deallocate(cpd_l)
    allocate(cpd_l(natm,3));   cpd_l = 0.d0
!    endif
    allocate(cpo_l(natm,3,3)); cpo_l = 0.d0
    allocate(numlay(natm)); numlay = 1
    allocate(atom_key(natm)); atom_key = 0

    allocate(pos_in(natm,3))
    allocate(cps_in(natm,3))
    if(sw_atom_excludable == ON)then
        allocate(exclusion_target(natm));exclusion_target = 0
        exclusion_criteria_min(:) = -1.e+30 ! we don't want to exclude atoms unless explicitly specified
        exclusion_criteria_max(:) = +1.e+30
    endif
    if(sw_extrapolate_charge==ON.or.sw_wf_predictor==ON)then
       allocate(cps_history(natm,3,3));cps_history=0.d0
    endif

! === KT_add === 2014/12/29
    if ( sw_set_initial_magmom_by_atom == ON ) then
       allocate( ionic_charge_atoms(natm) ); ionic_charge_atoms = -100.0d0
       if ( noncol ) then
          allocate( mag_moment0_atoms(natm,3) ); mag_moment0_atoms = -100.0
       else
          allocate( mag_moment0_atoms(natm,1) ); mag_moment0_atoms = -100.0
       endif
       if ( noncol ) then
          allocate( mag_direction0_atoms(natm,3) ); mag_direction0_atoms = 0.0d0
       endif
    endif
! ============== 2014/12/29

!<<===== ASMS === 2021/03/25
    if (vdw_method /= VDW_DFTD3 .and. ntyp_vdw>0) then
       allocate(ityp_vdw(natm)); ityp_vdw = 1
       allocate(fxyzvdw_l(natm,3))
    endif
!======= ASMS === 2021/03/25 >>

    if(ipri >= 3 .and. printable) write(nfout,'(" !** natm = ",i6," <<m_IS_alloc_pos_and_v>>")') natm
    call forcp_alloc()
  end subroutine m_IS_alloc_pos_and_v

  subroutine m_IS_dealloc_pos_and_v(nfout)
    integer, intent(in) :: nfout
    integer :: i
    deallocate(iwei)
    deallocate(imdtyp)
    deallocate(imdtypxyz)
    deallocate(ityp)
    deallocate(if_pdos)
    deallocate(if_aldos)
    deallocate(ihubbard)
    deallocate(iproj_group)

! =========================== added by K. Tagami ================= 11.0
    deallocate( itab_spinorbit_addition )
! ================================================================= 11.0
    deallocate(ionic_mass)
    deallocate(pos)
    deallocate(cps)
    if(sw_extrapolate_charge==ON.or.sw_wf_predictor==ON) deallocate(cps_history)
    deallocate(pos_in)
    deallocate(cps_in)
!    deallocate(cpd_l)
    deallocate(cpo_l)
    deallocate(numlay)
    deallocate(atom_key)
    if(sw_atom_excludable==ON) deallocate(exclusion_target)

!<<===== ASMS === 2021/03/25
    if (vdw_method /= VDW_DFTD3 .and. ntyp_vdw>0) then
       if ( allocated(ityp_vdw) ) deallocate(ityp_vdw)
       if ( allocated(fxyzvdw_l) ) deallocate(fxyzvdw_l)
    endif
!======= ASMS === 2021/03/25 >>

! === KT_add === 2014/12/29
    if ( sw_set_initial_magmom_by_atom == ON ) then
       deallocate( ionic_charge_atoms )
       deallocate( mag_moment0_atoms )
       if ( noncol ) deallocate( mag_direction0_atoms )
    endif
! ============== 2014/12/29
    if(num_regions>0) then
       do i=1,num_regions
          deallocate(regions(i)%forc)
          deallocate(regions(i)%target_atoms)
       enddo
       deallocate(regions)
    endif
  end subroutine m_IS_dealloc_pos_and_v

! === KT_add === 2014/12/29 & 13.2S
  subroutine m_IS_set_ionic_charge_atoms(nfout)
    integer, intent(in) :: nfout

    integer :: ia, it

    if ( .not. mag_moment0_atoms_is_defined ) return

    do ia = 1,natm
       if ( ionic_charge_atoms(ia) < -99.d0 ) then
          it = ityp(ia)
          ionic_charge_atoms(ia) = ionic_charge_atomtyp(it)
       endif
    end do

  end subroutine m_IS_set_ionic_charge_atoms

  subroutine m_IS_set_mag_moment0_atoms(nfout)
    integer, intent(in) :: nfout

    integer :: ia, it

    if ( .not. mag_moment0_atoms_is_defined ) return

    do ia = 1,natm
       if ( mag_moment0_atoms(ia,1) < -99.d0 ) then
          it = ityp(ia)
          if ( noncol ) then
             mag_moment0_atoms(ia,1:3) = mag_moment0_atomtyp(it,1:3)
          else
             mag_moment0_atoms(ia,1) = mag_moment0_atomtyp(it,1)
          endif
       endif
    end do

    if ( noncol ) then
       write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atoms) --- '
       write(nfout,*) '     id     name       ion         mx          my          mz'
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I8,3X,A7,4F12.6)') ia, speciesname(it), &
               &                         ionic_charge_atoms(ia), &
               &                         mag_moment0_atoms(ia,1), &
               &                         mag_moment0_atoms(ia,2), &
               &                         mag_moment0_atoms(ia,3)
       End Do
    else if ( sw_fix_global_quantz_axis == ON ) then
       write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atoms) --- '
       write(nfout,*) '     id     name       ion         mx          my          mz'
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I8,3X,A7,4F12.6)') ia, speciesname(it), &
               &           ionic_charge_atoms(ia), &
               &           mag_moment0_atoms(ia,1) *Global_Quantz_Axis_Fixed(1), &
               &           mag_moment0_atoms(ia,1) *Global_Quantz_Axis_Fixed(2), &
               &           mag_moment0_atoms(ia,1) *Global_Quantz_Axis_Fixed(3)
       End Do
    else
       write(nfout,*) '! ------------ Initial Charge/Magnetic Moment (atoms) --- '
       write(nfout,*) '     id     name       ion        moment '
       Do ia=1, natm
          it = ityp(ia)
          write(nfout,'(I8,3X,A7,2F12.6)') ia, speciesname(it), &
               &                         ionic_charge_atoms(ia), &
               &                         mag_moment0_atoms(ia,1)
       End Do
    endif

    write(nfout,*) '! ---------------------------------------------'

  end subroutine m_IS_set_mag_moment0_atoms
! ================= 2014/12/29

! ====== KT_add === 2014/06/10
  subroutine m_IS_reassgin_thermog
    integer :: i, ia, num
    integer, save :: ithermo_old = 1
    real(kind=DP) :: c1, c2, temp_now

    ithermo_now = 1

    if ( iteration_ionic > mdstep_at_end_thermostat(nrsv) ) then
       ithermo_now = nrsv
    else
       Do i=1, nrsv
          if ( iteration_ionic <= mdstep_at_end_thermostat(i) ) then
             ithermo_now = i
             exit
          endif
       End Do
    endif
!
    if ( ithermo_now /= ithermo_old ) then
       write(nfout,*) '***** Thermo group is changed to ', ithermo_now
    endif
    ithermo_old = ithermo_now

    num = 0
    Do ia=1, natm
       if (imdtyp(ia) > NOSE_HOOVER ) then
          imdtyp(ia) = NOSE_HOOVER + ithermo_now
          num = num +1
       end if
    End Do
    natm_per_thermo = 0
    natm_per_thermo(ithermo_now) = num
!
    c1 = iteration_ionic -mdstep_at_start_thermostat(ithermo_now)
    c2 = mdstep_at_end_thermostat(ithermo_now) -iteration_ionic
!
    temp_now = temp_at_end_thermostat(ithermo_now) *c1 &
         &   + temp_at_start_thermostat(ithermo_now) *c2
    temp_now = temp_now /( c2+c1)
!
    write(nfout,*) "!** temperature is now ", temp_now , " K"

    tkb(ithermo_now) = temp_now *const_KB
    if ( qmass(ithermo_now) < 1.0D-4 ) then
       qmass(ithermo_now) = (2.d0*3.d0*natm_per_thermo(ithermo_now)*(tdamp/PAI2)**2) &
            &             *tkb(ithermo_now)
    endif

  end subroutine m_IS_reassgin_thermog
! =============== 2014/06/10

  subroutine m_IS_recover(n,exclude)
    integer, intent(in) :: n
    logical, dimension(n), intent(in) :: exclude
    integer :: i,ii
    ii=0
    do i=1,n
       if(exclude(i)) cycle
       ii=ii+1
       iwei(ii) = config_buf(i)%iwei
       imdtyp(ii) = config_buf(i)%imdtyp
       imdtypxyz(ii,:) = config_buf(i)%imdtypxyz(:)
       ityp(ii) = config_buf(i)%ityp
       if_pdos(ii) = config_buf(i)%if_pdos
       if_aldos(ii) = config_buf(i)%if_aldos
       ihubbard(ii) = config_buf(i)%ihubbard
       iproj_group(ii) = config_buf(i)%iproj_group
       ionic_mass(ii) = config_buf(i)%ionic_mass
       numlay(ii) = config_buf(i)%numlay
       if(sw_atom_excludable==ON) exclusion_target(ii) = config_buf(i)%exclusion_target
       pos(ii,:)  = config_buf(i)%pos(:)
       cps(ii,:) = config_buf(i)%cps(:)
       pos_in(ii,:) = config_buf(i)%pos_in(:)
       cps_in(ii,:) = config_buf(i)%cps_in(:)
       cpd_l(ii,:) = config_buf(i)%cpd_l(:)
       cpo_l(ii,:,:) = config_buf(i)%cpo_l(:,:)
    enddo
  end subroutine m_IS_recover

  subroutine m_IS_store_current_config()
    integer :: i
    if(allocated(config_buf)) deallocate(config_buf)
    allocate(config_buf(natm))
    nconfig_buf = natm
    do i=1,natm
       !config_buf(i)%element = species_work(i)
       config_buf(i)%iwei = iwei(i)
       config_buf(i)%imdtyp = imdtyp(i)
       config_buf(i)%imdtypxyz(:) = imdtypxyz(i,:)
       config_buf(i)%ityp = ityp(i)
       config_buf(i)%if_pdos = if_pdos(i)
       config_buf(i)%if_aldos = if_aldos(i)
       config_buf(i)%ihubbard = ihubbard(i)
       config_buf(i)%iproj_group = iproj_group(i)
       config_buf(i)%ionic_mass = ionic_mass(i)
       config_buf(i)%numlay = numlay(i)
       config_buf(i)%pos(:) = pos(i,:)
       config_buf(i)%cps(:) = cps(i,:)
       config_buf(i)%pos_in(:) = pos_in(i,:)
       config_buf(i)%cps_in(:) = cps_in(i,:)
       config_buf(i)%cpd_l(:) = cpd_l(i,:)
       config_buf(i)%cpo_l(:,:) = cpo_l(i,:,:)
       !config_buf(i)%nvalence = ival(ityp(i))
    enddo
  end subroutine m_IS_store_current_config

  logical function m_IS_change_natm()
     integer :: i
     m_IS_change_natm = .false.
     if(.not.m_IS_natm_can_change()) return
     natmorg = natm
     if(sw_atom_excludable==ON)then
        do i=1,natm
           if(out_of_bounds(i))then
               if(printable) then
                  write(nfout,'(a,i8,a)') ' !** atom no ',i,'is out of bounds, and will be excluded.'
               endif
               call m_IS_remove_atom(i)
           endif
        enddo
     endif
     if(mod(iteration_ionic,addition_frequency)==0.and.icond/=CONTINUATION)then
        neg_incre = 0
        if(natom_reservoir<curr_atom_reservoir) then
           if(printable .and. iprimd>=1) write(nfout,'(a)') ' !** atom reservoir exhausted.'
        else
           if(printable .and. iprimd>=1) then
              write(nfout,'(a)') ' !** adding new atom(s) from the atom reservoir'
              write(nfout,'(a)') ' !** element fracx fracy fracz cartx carty cartz'
              do i=1,natm_per_group(curr_atom_reservoir)
                 call print_atom(atom_reservoir(atomid_in_group(curr_atom_reservoir,i)))
              enddo
           endif
           do i=1,natm_per_group(curr_atom_reservoir)
              call m_IS_add_atom(atom_reservoir(atomid_in_group(curr_atom_reservoir,i)))
           enddo
           curr_atom_reservoir = curr_atom_reservoir+1
           if(curr_atom_reservoir>natom_group .and. sw_rotate_reservoir==ON) curr_atom_reservoir = 1
           m_IS_change_natm = .true.
        endif
     endif
  end function m_IS_change_natm

  logical function out_of_bounds(iatom)
     integer, intent(in) :: iatom
     integer :: ii
     out_of_bounds = .false.
     do ii=1,3
        if (cps(iatom,ii) .lt. exclusion_criteria_min(ii) ) out_of_bounds = .true.
        if (cps(iatom,ii) .gt. exclusion_criteria_max(ii) ) out_of_bounds = .true.
        if (out_of_bounds) return
     enddo
  end function out_of_bounds

  subroutine m_IS_wd_curr_atom_reservoir(nfcntn)
     integer, intent(in) :: nfcntn
     if(mype==0)then
        write(nfcntn,*) tag_curr_atom_reservoir
        write(nfcntn,'(i8)') curr_atom_reservoir
     endif
  end subroutine m_IS_wd_curr_atom_reservoir

  subroutine m_IS_rd_curr_atom_reservoir(nfcntn)
     integer, intent(in) :: nfcntn
     logical             :: EOF_reach, tag_is_found
     if(mype==0)then
        call rewind_to_tag0(nfcntn,len(tag_curr_atom_reservoir),tag_curr_atom_reservoir &
     &  , EOF_reach, tag_is_found,str,len_str)
        if(.not.tag_is_found) then
           curr_atom_reservoir = 1
        else
           read(nfcntn,*) curr_atom_reservoir
        endif
     endif
     if(npes>1) call mpi_bcast(curr_atom_reservoir,1,mpi_integer,0,MPI_CommGroup,ierr)
  end subroutine m_IS_rd_curr_atom_reservoir

  integer function m_IS_get_neg_incre()
     m_IS_get_neg_incre = neg_incre
     neg_incre=0
  end function m_IS_get_neg_incre

  logical function m_IS_natm_can_change()
    integer :: i
    if(sw_atom_excludable==ON) then
       do i=1,natm
          if(out_of_bounds(i))then
             m_IS_natm_can_change = .true.
             return
          endif
       enddo
    endif
    if(allocated(atom_reservoir) .and. natom_reservoir.gt.0 .and. addition_frequency.gt.1)then
       if(natom_reservoir>=curr_atom_reservoir) then
          m_IS_natm_can_change = .true.
          return
       endif
    endif
    m_IS_natm_can_change = .false.
  end function

  subroutine m_IS_set_ival(nt,iv)
     integer, intent(in) :: nt
     real(kind=DP), dimension(nt),intent(in) :: iv
     integer :: i
     if(.not.allocated(ival)) allocate(ival(nt))
     ival = iv
     if(natom_reservoir>0)then
        do i=1,natom_reservoir
           atom_reservoir(i)%nvalence = ival(atom_reservoir(i)%ityp)
           call print_atom(atom_reservoir(i))
        enddo
     endif
  end subroutine m_IS_set_ival

  subroutine init_atom(i,theatom)
    integer, intent(in) :: i
    type(atomic_configuration_t),intent(inout) :: theatom
    theatom%id = i
    theatom%element = ''
    theatom%group = -1
    theatom%iwei = 1
    theatom%imdtyp = 0
    theatom%imdtypxyz = 0
    theatom%ityp = 1
    theatom%if_pdos = 1
    theatom%if_aldos = 1
    theatom%ihubbard = 0
    theatom%iproj_group = 0
    theatom%exclusion_target = 0
    theatom%ionic_mass = -1.d0
    theatom%numlay = 1
    theatom%pos = 0.0d0
    theatom%cps = 0.0d0
    theatom%pos_in = 0.0d0
    theatom%cps_in = 0.0d0
    theatom%cpd_l = 0.0d0
    theatom%cpo_l = 0.0d0
    theatom%nvalence = 0
  end subroutine init_atom

  subroutine print_atom(theatom)
    type(atomic_configuration_t),intent(in) :: theatom
    integer :: i
    if(printable)then
        write(nfout,'(a,i8,i3,6f10.5,f5.1,i3)') '    '//trim(theatom%element),theatom%id,theatom%ityp,&
       & (theatom%pos(i),i=1,3),(theatom%cps(i),i=1,3),theatom%nvalence,theatom%group
    endif
  end subroutine print_atom

  subroutine resolve_dftd3_parameters()
    integer :: ierr
    integer, parameter     :: len_str = 132
    character(len=len_str) :: str
    logical :: tag_is_found, EOF_reach
    integer :: i,j,ii,jj,kk,ll,ntmp
    real(kind=DP) :: val1,val2,val3

    if(dftd3par%read_pars) return

    dftd3par%k1 = 16.d0
    dftd3par%k2 = 4.d0/3.d0
    dftd3par%k3 = 4.d0
    dftd3par%s6 = 1.d0
    dftd3par%alpha6 = 14.d0
    dftd3par%alpha8 = 16.d0
    dftd3par%sr8 = 1.d0

    if ( dftd3_damping_function == 0 ) then    ! zero damping
!   parameters for PBE
       dftd3par%s8 = 0.722d0
       dftd3par%sr6 = 1.217d0
       dftd3par%sr9 = 4.d0/3.d0
       dftd3par%alpha9 = 16.d0

       if (xctype == 'revpbe') then
          dftd3par%s8 = 1.010d0
          dftd3par%sr6 = 0.923d0
       else if (xctype == 'rpbe  ') then
          dftd3par%s8 = 0.514d0
          dftd3par%sr6 = 0.872d0
       else if (xctype == 'pbe0  ') then
          dftd3par%s8 = 0.928d0
          dftd3par%sr6 = 1.287d0
       else if (xctype == 'hse06 ') then
          dftd3par%s8 = 0.109d0
          dftd3par%sr6 = 1.129d0
       endif

    else if ( dftd3_damping_function == 1 ) then  ! BJ damping
!   parameters for PBE
       dftd3par%s8 =  0.7875d0
       dftd3par%a1 = 0.4289d0
       dftd3par%a2 =  4.4407d0

    endif       

!    dftd3par%rcut_nc = 20.d0
    dftd3par%rcut_nc = 40.d0
    dftd3par%rcut_vdw = 40.d0

    call m_Files_open_dftd3par()
    if(mype == 0)then
      call rewind_to_tag0(nfdftd3par,len('nlines'),'nlines' &
     & , EOF_reach, tag_is_found, str, len_str)
      read(nfdftd3par,*) dftd3par%nlines
      call rewind_to_tag0(nfdftd3par,len('max_elem'),'max_elem' &
     & , EOF_reach, tag_is_found, str, len_str)
      read(nfdftd3par,*) dftd3par%maxelem
      call rewind_to_tag0(nfdftd3par,len('maxmaxnc'),'maxmaxnc' &
     & , EOF_reach, tag_is_found, str, len_str)
      read(nfdftd3par,*) dftd3par%maxmaxnc
    endif
    if(npes>1)then
       call mpi_bcast(dftd3par%nlines,1,mpi_integer,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%maxelem,1,mpi_integer,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%maxmaxnc,1,mpi_integer,0,MPI_CommGroup,ierr)
    endif
    allocate(dftd3par%c6ab(dftd3par%maxelem,dftd3par%maxelem,dftd3par%maxmaxnc,dftd3par%maxmaxnc))
    allocate(dftd3par%nc1 (dftd3par%maxelem,dftd3par%maxelem,dftd3par%maxmaxnc,dftd3par%maxmaxnc))
    allocate(dftd3par%nc2 (dftd3par%maxelem,dftd3par%maxelem,dftd3par%maxmaxnc,dftd3par%maxmaxnc))
    dftd3par%c6ab=0.d0;dftd3par%nc1=0.d0;dftd3par%nc2=0.d0
    allocate(dftd3par%maxnc(dftd3par%maxelem));dftd3par%maxnc=0
    allocate(dftd3par%r2r4(dftd3par%maxelem));dftd3par%r2r4=0.d0
    allocate(dftd3par%covrad(dftd3par%maxelem));dftd3par%covrad=0.d0
    allocate(dftd3par%r0ab(dftd3par%maxelem,dftd3par%maxelem));dftd3par%r0ab=0.d0
    if(mype == 0)then
      call rewind_to_tag0(nfdftd3par,len('maxnc'),'maxnc' &
     & , EOF_reach, tag_is_found, str, len_str)
      do i=1,dftd3par%maxelem
         read(nfdftd3par,*) dftd3par%maxnc(i)
      enddo
      call rewind_to_tag0(nfdftd3par,len('r0ab'),'r0ab' &
     & , EOF_reach, tag_is_found, str, len_str)
      do i=1,dftd3par%maxelem
      do j=1,dftd3par%maxelem
         read(nfdftd3par,*) ii,jj,dftd3par%r0ab(i,j)
      enddo
      enddo
      call rewind_to_tag0(nfdftd3par,len('covrad'),'covrad' &
     & , EOF_reach, tag_is_found, str, len_str)
      do i=1,dftd3par%maxelem
         read(nfdftd3par,*) dftd3par%covrad(i)
      enddo
      call rewind_to_tag0(nfdftd3par,len('r2r4'),'r2r4' &
     & , EOF_reach, tag_is_found, str, len_str)
      do i=1,dftd3par%maxelem
         read(nfdftd3par,*) dftd3par%r2r4(i)
      enddo
    endif
    if(npes>1)then
       call mpi_bcast(dftd3par%maxnc,dftd3par%maxelem,mpi_integer,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%r0ab,dftd3par%maxelem*dftd3par%maxelem,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%covrad,dftd3par%maxelem,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%r2r4,dftd3par%maxelem,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif
    if(mype == 0)then
      call rewind_to_tag0(nfdftd3par,len('c6ab'),'c6ab' &
     & , EOF_reach, tag_is_found, str, len_str)
      do i=1,dftd3par%nlines
         read(nfdftd3par,*) val1,ii,jj,kk,ll,val2,val3
         dftd3par%c6ab(ii,jj,kk,ll) = val1
         dftd3par%c6ab(jj,ii,ll,kk) = val1

         dftd3par%nc1(ii,jj,kk,ll) = val2
         dftd3par%nc1(jj,ii,ll,kk) = val3

         dftd3par%nc2(ii,jj,kk,ll) = val3
         dftd3par%nc2(jj,ii,ll,kk) = val2
      enddo
    endif
    if(npes>1)then
       ntmp = size(dftd3par%c6ab)
       call mpi_bcast(dftd3par%c6ab,ntmp,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%nc1 ,ntmp,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(dftd3par%nc2 ,ntmp,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif
    dftd3par%read_pars = .true.
    call m_Files_close_dftd3par()
  end subroutine resolve_dftd3_parameters

  subroutine read_dftd3_parameters()
    integer :: iret
    real(kind=DP) :: dret
    integer :: f_selectBlock, f_selectParentBlock, f_getRealValue
    if(f_selectBlock('dftd3') == 0) then
       if(f_getRealValue(tag_k1,dret,'') == 0)       dftd3par%k1 = dret
       if(f_getRealValue(tag_k2,dret,'') == 0)       dftd3par%k2 = dret
       if(f_getRealValue(tag_k3,dret,'') == 0)       dftd3par%k3 = dret
       if(f_getRealValue(tag_s6,dret,'') == 0)       dftd3par%s6 = dret
       if(f_getRealValue(tag_s8,dret,'') == 0)       dftd3par%s8 = dret
!       if(f_getRealValue(tag_sr6,dret,'bohr') == 0)  dftd3par%sr6 = dret
!       if(f_getRealValue(tag_sr8,dret,'bohr') == 0)  dftd3par%sr8 = dret
       if(f_getRealValue(tag_sr6,dret,'') == 0)  dftd3par%sr6 = dret
       if(f_getRealValue(tag_sr8,dret,'') == 0)  dftd3par%sr8 = dret
       if(f_getRealValue(tag_alpha6,dret,'') == 0)   dftd3par%alpha6 = dret
       if(f_getRealValue(tag_alpha8,dret,'') == 0)   dftd3par%alpha8 = dret
       if(f_getRealValue(tag_rcut,dret,'bohr') == 0) dftd3par%rcut_vdw = dret
       if(f_getRealValue(tag_rcut_nc,dret,'bohr') == 0) dftd3par%rcut_nc = dret
       if(f_getRealValue(tag_a1,dret,'') == 0)       dftd3par%a1 = dret
       if(f_getRealValue(tag_a2,dret,'bohr') == 0)   dftd3par%a2 = dret
       iret = f_selectParentBlock()
    endif

    if(printable)then
       write(nfout,'(a)')       ' !** dftd3 parameters'
       write(nfout,'(a,d15.5)') ' !** k1     ', dftd3par%k1
       write(nfout,'(a,d15.5)') ' !** k2     ', dftd3par%k2
       write(nfout,'(a,d15.5)') ' !** k3     ', dftd3par%k3
       write(nfout,'(a,d15.5)') ' !** s6     ', dftd3par%s6
       write(nfout,'(a,d15.5)') ' !** s8     ', dftd3par%s8
       write(nfout,'(a,d15.5)') ' !** sr6    ', dftd3par%sr6
       write(nfout,'(a,d15.5)') ' !** sr8    ', dftd3par%sr8
       write(nfout,'(a,d15.5)') ' !** alpha6 ', dftd3par%alpha6
       write(nfout,'(a,d15.5)') ' !** alpha8 ', dftd3par%alpha8
       write(nfout,'(a,d15.5)') ' !** rcut   ', dftd3par%rcut_vdw
       write(nfout,'(a,d15.5)') ' !** rcut_nc', dftd3par%rcut_nc
       if ( dftd3_damping_function == 1 ) then
          write(nfout,'(a,d15.5)') ' !** a1     ', dftd3par%a1
          write(nfout,'(a,d15.5)') ' !** a2     ', dftd3par%a2
       endif
    endif
  end subroutine read_dftd3_parameters

  subroutine alloc_species_work()
    allocate(species_work(natm)); species_work = ""
    allocate(species_indp(natm)); species_indp = ""
  end subroutine alloc_species_work

  subroutine dealloc_species_work()
    if(allocated(species_work)) deallocate(species_work)
    if(allocated(species_indp)) deallocate(species_indp)
  end subroutine dealloc_species_work

  subroutine alloc_speciesname()
    allocate(speciesname(ntyp)); speciesname=''
  end subroutine alloc_speciesname

  subroutine dealloc_speciesname()
    if(allocated(speciesname)) deallocate(speciesname)
  end subroutine dealloc_speciesname

  subroutine alloc_species_vdw_work()
    allocate(species_vdw_work(natm)); species_vdw_work = ""
    allocate(species_vdw_indp(natm)); species_vdw_indp = ""
  end subroutine alloc_species_vdw_work

  subroutine dealloc_species_vdw_work()
    if(allocated(species_vdw_work)) deallocate(species_vdw_work)
    if(allocated(species_vdw_indp)) deallocate(species_vdw_indp)
  end subroutine dealloc_species_vdw_work

  subroutine alloc_speciesname_vdw()
    allocate(speciesname_vdw(ntyp_vdw))
    speciesname_vdw(1:ntyp_vdw) = species_vdw_indp(1:ntyp_vdw)
  end subroutine alloc_speciesname_vdw

! ================================= KT_add ================ 13.0B
  subroutine dealloc_speciesname_vdw()
    if ( allocated(speciesname_vdw) ) deallocate(speciesname_vdw)
  end subroutine dealloc_speciesname_vdw
! ========================================================= 13.0B

  subroutine m_IS_cp_cps2cpo
    cpo_l(1:natm,1:3,1) = cps
  end subroutine m_IS_cp_cps2cpo

  subroutine m_IS_cp_cps2cpo3
    cpo_l(1:natm,1:3,3) = cps
  end subroutine m_IS_cp_cps2cpo3

  subroutine m_IS_fire(forc_l_in)
    real(kind=DP), dimension(natm,3), intent(in) :: forc_l_in
    call m_IS_fire_core(nfout,natm,iteration_ionic,cps,cpd_l,forc_l_in,imdtypxyz)
    if (nrigid_bodies>0) call m_IS_rb_dynamics(forc_l_in, imdalg)
  end subroutine m_IS_fire

  subroutine m_IS_md(mdalg,forc_l_in)
! Modified by T. Yamasaki, May/18/2023
!     Codes related to RIGID_BODY control functions are revised.
    integer, intent(in)     :: mdalg
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP), allocatable, dimension(:,:)   :: forc_l

    integer :: id_sname = -1
    integer :: i,ierr
    call tstatc0_begin('m_IS_md ',id_sname)

    call md1_alloc(mdalg) ! allocate cpd_old, and ipcpd

    allocate(forc_l(natm,3))
    forc_l = forc_l_in

!!$    call check_if_bondlength_fix_exist(ibondlengthfix)
!!$    if(mdmode == ORDINA .and. nfcatm == 0 .and. ibondlengthfix == 0) then

    if(ipriinputfile>=DEBUGPRINTLEVEL) write(nfout,'(" !* sw_fix_bond = ",i3)') sw_fix_bond !!! T. Y. 2023.04.27
    if(sw_fix_bond==ON) then
      allocate(oldcps(natm,3))
      oldcps = cps
    endif
    if(num_planes_atoms_are_fixed>=1) then
       allocate(tmass(num_planes_atoms_are_fixed));tmass= 0.d0
       allocate(rtmass(num_planes_atoms_are_fixed)); rtmass = 0.d0
    end if
    if(num_planes_atoms_are_fixed_cog>=1) then
       allocate(fcg_mdfy_cog(num_planes_atoms_are_fixed_cog,1:3))
       allocate(fcg_cog(num_planes_atoms_are_fixed_cog,1:3))
       allocate(tmass_cog(num_planes_atoms_are_fixed_cog)) ; tmass_cog = 0.d0
       allocate(rtmass_cog(num_planes_atoms_are_fixed_cog)); rtmass_cog = 0.d0

    end if
    if(num_planes_atoms_are_fixed_rb>=1) then
       allocate(fcg_rb(num_planes_atoms_are_fixed_rb,1:3))
       allocate(fcg_mdfy_rb(num_planes_atoms_are_fixed_rb,1:3))
       allocate(tmass_rb(num_planes_atoms_are_fixed_rb)); tmass_rb = 0.d0
       allocate(rtmass_rb(num_planes_atoms_are_fixed_rb)); rtmass_rb = 0.d0
    end if
    if(iprimd >= DEBUGPRINTLEVEL) then
       write(nfout,'(" constraint_type = ",i5)') constraint_type
       write(nfout,'(" mdmode = ",i5)') mdmode
    end if

    if(iprimd >= DEBUGPRINTLEVEL) call w_cps(0)

    if(mdmode == ORDINA) then
                                if(iprimd >= DEBUGPRINTLEVEL+1) write(nfout,'(" mdmode == ORDINA")')
       iter_md = iter_md + 1
! -- Following 3 lines have been revised by Mamoru Usami, Oct. 2004. -->>
       if (mdalg .ne. VERLET) then
          cpd_old = cpd_l
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cpd_l(1,"before evolve_velocities")
          call check_constraint(forc_l)  ! -> fcg_{cog|rb},fcg_mdfy_{cog|rb},tmass_{cog|rb} &
          !                                       : modified group force along in a plane, tmass: group mass
          if(imdalg == STEEPEST_DESCENT .and. mode_fi_coefficient == OFF) then
             fi_coefficient = get_fi_coefficient(sqrt(sum(forc_l(1:natm,1:3)**2)))
                                if(iprimd >= 1) write(nfout,'(a19,f20.10)') ' fi_coefficient =  ', fi_coefficient
          end if
!!$          if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_fi(forc_l)
          if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_hyperplane(forc_l)
          call evolve_velocities(mdalg,forc_l) !-(m_Ionic_System)  cpd_l = cpd_l + dtio/mass * forc_l
       end if

                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cpd_l(2,"after evolve_velocities")
! <<--
       if(mdalg == QUENCHED_MD) call quench_velocities(forc_l) ! -->cpd_l
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cpd_l(2,"after quench_velocities")

!!$       if(constraint_type == COG_FIX_L) call correct_cog_motion()
       if(constraint_type==COG_FIX_L .or. constraint_type==COG_and_RIGID_BODY_FIX_L &
            & .or. constraint_type==RIGID_BODY_FIX_IN_A_PLANE) call correct_cog_andor_rb_motion() ! -->cpd_l

       if(mdalg == QUENCHED_MD .and. apply_constant_force >=1 ) call evolve_velocities_constant_force()

       if(mdalg == VERLET .and. &
            &(nbztyp == HEXAGONAL .or. nbztyp == ORTHORHOMBIC .or. &
            & nbztyp >= HEX1fold)) then
          if(nopr_supercell <= nopr .or. sw_supercell_symmetry /= ON) then
             call fd_symmetrize(natm2,natm,natm,napt,nopr+af,nopr,op,iwei &
                  & , cpd_l, cpd_old, ipcpd)  ! -(b_Ionic_System) -->cpd_l
          else
!!$             call fd_symmetrize(natm2,natm,natm,napt_supercell,nopr_supercell,op,nopr+af,iwei &
!!$                  & , cpd_l, cpd_old, ipcpd, iop_supercell)  ! -(b_Ionic_System) -->cpd_l
             call fd_supercell_symmetrize(natm2,natm,natm,napt_supercell,nopr_supercell,op,nopr+af,iwei &
                  & , cpd_l, cpd_old, ipcpd, iop_supercell)
          end if
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cpd_l(2,"after fd_symmetrize")
       end if
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cps(1)

       if(constraint_type /= BONDLENGTH_FIX) then
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cpd_l(2,"before cps+dtio*cpd_l")
!!$          call evolve_cps()
          cps = cps + dtio*cpd_l
          if(mdalg == QUENCHED_MD .and. apply_constant_move >= 1) then
!!$             call move_atoms_normal_to_plane() ! -> cps
             call evolve_cps_constant_move() ! -> cps
          else
          end if
                                if(iprimd >= DEBUGPRINTLEVEL+1) call w_cps(2)
       else
          cpd_old = cpd_l
          call correct_fixed_bond() ! -> cpd_old, cps = cps+dtio*cpd_old, fcvect,cpd,cps
       end if

    else if(mdmode == CNSTRA) then
!!$       call evolve_cps_constant_move()
       call move_atoms_normal_to_plane()
    else
       call phase_error_with_msg(nfout,' Invalid value of mdmode <<m_IS_md>>',__LINE__,__FILE__)
    end if
    if (sw_fix_bond==ON) then
      call fixed_bond_coords()
      call update_bond_force()
      call update_bond_dsigma_old()
   endif

    if (sw_shift_velocities==ON.and.mdalg==VERLET) call shift_velocities(1)

    if(sw_fix_bond==ON) deallocate(oldcps)
    deallocate(forc_l)
    call md1_dealloc()
    if(num_planes_atoms_are_fixed_cog>=1) deallocate(fcg_mdfy_cog,fcg_cog,rtmass_cog,tmass_cog)
    if(num_planes_atoms_are_fixed_rb>=1)  deallocate(fcg_mdfy_rb, fcg_rb, rtmass_rb, tmass_rb)
    if(num_planes_atoms_are_fixed>=1)  deallocate(rtmass,tmass)
    call tstatc0_end(id_sname)

  contains
    subroutine w_cps(ip)
      integer, intent(in) :: ip
      integer :: i
      do i = 1, natm
         write(nfout,'(" cps(",i5,") = ",3f8.4, " cpd_l = ",3f8.4, " ip = ",i3)') i, cps(i,1:3),cpd_l(i,1:3),ip
      end do
    end subroutine w_cps

    subroutine w_cpd_l(iswitch,name)
      integer, intent(in) :: iswitch
      character(len=*), intent(in) :: name
      integer :: i
      if(iswitch==1) then
         write(nfout,'(" -- cpd_l (",a25,")--")') name
         do i  = 1, natm
            write(nfout,'(i5,3f18.7)') i, cpd_l(i,1:3)
         end do
      else if(iswitch==2) then
         write(nfout,'(" -- cpd_l, cpd_old (",a25,")--")') name
         do i  = 1, natm
            write(nfout,'(i5,6f18.7)') i, cpd_l(i,1:3), cpd_old(i,1:3)
         end do
      end if
    end subroutine w_cpd_l

    subroutine correct_fixed_bond()
      integer :: ib, ia1, ia2, i
      real(kind=DP) :: pdot, fac
      real(kind=DP),allocatable,dimension(:) :: frc3,frc4,frc1,frc2,cps1,cps2,fcv1,fcv2
      allocate(frc1(3));allocate(frc2(3));allocate(frc3(3));allocate(frc4(3))
      allocate(cps1(3));allocate(cps2(3));allocate(fcv1(3));allocate(fcv2(3))

      do ib = 1, num_fixed_bonds
         ia1 = bondlength_fix_set(1,ib); ia2 = bondlength_fix_set(2,ib)
         pdot = dot_product(forc_l(ia1,1:3)-forc_l(ia2,1:3), fcvect(ia1,1:3)) &
              & /fcvect(ia1,4)**2*0.5d0
         frc3=forc_l(ia1,1:3) - pdot*fcvect(ia1,1:3)
         frc4=forc_l(ia2,1:3) + pdot*fcvect(ia1,1:3)
         frc1=frc3+frc4
         frc2=frc3-frc4  ! frc2 = forc_l(ia1,:)-forc_l(ia2,:) - 2*pdot*fcvect(ia1,1:3)
         pdot = dot_product(frc2,fcvect(ia1,1:3))
         if(dabs(pdot) > 1.d-6) then
            if(iprimd >=1)then
               write(nfout,'(" !! frc2(:)       = ",3f12.6," <<correct_fixed_bond>>")') frc2
               write(nfout,'(" !! fcvect(ia1,:) = ",3f12.6," <<correct_fixed_bond>>")') fcvect(ia1,1:3)
            end if
            if(printable) write(nfout,'(" fcvect times frc2(pdot)=",d16.8)') pdot
            call phase_error_with_msg(nfout,' fcvect is not normal to frc2 <<m_IS_md.correct_fixed_bond>>', &
            __LINE__,__FILE__)
         end if
         pdot = dot_product(cpd_old(ia1,:)-cpd_old(ia2,:),fcvect(ia1,1:3))/fcvect(ia1,4)**2*0.5d0
         cpd_old(ia1,:) = cpd_old(ia1,:) - pdot*fcvect(ia1,1:3)
         cpd_old(ia2,:) = cpd_old(ia2,:) + pdot*fcvect(ia1,1:3)
         cps1=cpd_old(ia1,:)+cpd_old(ia2,:)
         cps2=cpd_old(ia1,:)-cpd_old(ia2,:)
         pdot = dot_product(cps2,fcvect(ia1,1:3))
         if(iprimd >= 2) then
            write(nfout,'(" !! cps(    ",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(cps(ia1,1:3))
            write(nfout,'(" !! cps(    ",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia2,(cps(ia2,1:3))
            write(nfout,'(" !! cpd_l(  ",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(cpd_l(ia1,1:3))
            write(nfout,'(" !! cpd_l(  ",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia2,(cpd_l(ia2,1:3))
            write(nfout,'(" !! cpd_old(",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(cpd_old(ia1,1:3))
            write(nfout,'(" !! cpd_old(",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia2,(cpd_old(ia2,1:3))
            write(nfout,'(" !! cps1            ) = ",3f12.6," <<correct_fixed_bond>>")') cps1
            write(nfout,'(" !! cps2            ) = ",3f12.6," <<correct_fixed_bond>>")') cps2
            write(nfout,'(" !! fcvect(",i5,",:)  = ",4f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(fcvect(ia1,1:4))
            write(nfout,'(" !! fcvect(",i5,",:)  = ",4f12.6," <<correct_fixed_bond>>")') &
                 & ia2,(fcvect(ia2,1:4))
            write(nfout,'(" !! ----------")')
         end if

         if(dabs(pdot).gt.1.d-6) then
            if(printable) write(nfout,'(" fcvect times cps2 (pdot)=",d16.8)') pdot
            call phase_error_with_msg(nfout,' fcvect is not normal to cps2 <<m_IS_md.correct_fixed_bond>>', &
            __LINE__,__FILE__)
         end if
         do i = 1, 3
            if(cps1(i)*frc1(i) < -1.d-6) frc1(i) = 0.d0
            if(cps2(i)*frc2(i) < -1.d-6) frc2(i) = 0.d0
         end do
         fac = dtio/amion(ityp(ia1))
         cps1 = cps1 + fac*frc1
         cps2 = cps2 + fac*frc2
         cpd_old(ia1,:) = (cps1+cps2)*0.5d0
         cpd_old(ia2,:) = (cps1-cps2)*0.5d0
      end do

      if(iprimd >= 2) then
         write(nfout,'(" !! ---")')
         do ia1 = 1, natm
            write(nfout,'(" !! cpd_old(",i5," :) = ",3f12.6," << correct_fixed_bond>>")') ia1,cpd_old(ia1,1:3)
         end do
      end if

      cps = cps + dtio*cpd_old

      do ib = 1, num_fixed_bonds
         ia1 = bondlength_fix_set(1,ib); ia2 = bondlength_fix_set(2,ib)
         fcv1(1:3) = (cps(ia1,1:3) - cps(ia2,1:3))*0.5d0        ! = (c1-c2)/2
         pdot = dsqrt(dot_product(fcv1,fcv1))*2.d0              ! = |c1-c2|
         fcv2(1:3) = (cps(ia1,1:3) + cps(ia2,1:3))*0.5d0        ! = (c1+c2)/2

!!$         fcvect(ia1,1:3) = (cps(ia1,1:3) - cps(ia2,1:3))*0.5d0        ! = (c1-c2)/2
!!$         pdot = dsqrt(dot_product(fcvect(ia1,1:3),fcvect(ia1,1:3)))*2.d0  ! = |c1-c2|
!!$         fcvect(ia2,1:3) = (cps(ia1,1:3) + cps(ia2,1:3))*0.5d0        ! = (c1+c2)/2
         if(iprimd >= 2) then
            write(nfout,'(" !! fcv1(1:3) = ",3f12.6," <<correct_fixed_bond>>")') fcv1
            write(nfout,'(" !! fcv2(1:3) = ",3f12.6," <<correct_fixed_bond>>")') fcv2
            write(nfout,'(" !! cps(",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(cps(ia1,1:3))
            write(nfout,'(" !! cps(",i5,",:) = ",3f12.6," <<correct_fixed_bond>>")') &
                 & ia2,(cps(ia2,1:3))
            write(nfout,'(" !! ----------")')
         end if
         cps1 = cps(ia1,:)
         cps2 = cps(ia2,:)
         cps(ia1,:) = fcv2(1:3)+fcv1(1:3)*fcvect(ia1,4)/pdot ! =(c1+c2)/2 + (c1-c2)/2*f4/|c1-c2|
         cps(ia2,:) = fcv2(1:3)-fcv1(1:3)*fcvect(ia1,4)/pdot ! =(c1+c2)/2 - (c1-c2)/2*f4/|c1-c2|
         cpd_old(ia1,:) = cpd_old(ia1,:) + (cps(ia1,:)-cps1)
         cpd_old(ia2,:) = cpd_old(ia2,:) + (cps(ia2,:)-cps2)
         pdot = dot_product(cps(ia1,:)-cps(ia2,:),cps(ia1,:)-cps(ia2,:)) ! =|(c1-c2)*f4/|c1-c2||**2
         if(dabs(dsqrt(pdot)-fcvect(ia1,4)) > 1.d-10) then
            if(printable) then
               write(nfout,'(" !! imdtyp=BONDLENGTH_FIX normalization error")')
               write(nfout,'(" !!   ia1, ia2 = ",2i8)') ia1, ia2
               write(nfout,'(" !!   fcvect(ia1,4) = ",d16.8)') fcvect(ia1,4)
               write(nfout,'(" !!   pdot = ",d16.8)') pdot
               write(nfout,'(" !! fcvect(",i5,",:) = ",4f12.6," <<correct_fixed_bond>>")') &
                    & ia1,(fcvect(ia1,1:4))
               write(nfout,'(" !! fcvect(",i5,",:) = ",4f12.6," <<correct_fixed_bond>>")') &
                    & ia2,(fcvect(ia2,1:4))
            end if
            call phase_error_with_msg(nfout,' normalization error <<m_IS_md.correct_fixed_bond>>',__LINE__,__FILE__)
         end if
         if(iprimd >= 3) &
              & write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<correct_fixed_bond>>")') &
                 & ia1,(fcvect(ia1,i),i=1,4)
      end do

      cpd_l = cpd_old

      deallocate(fcv2);deallocate(fcv1);deallocate(cps2);deallocate(cps1)
      deallocate(frc4);deallocate(frc3);deallocate(frc2);deallocate(frc1)
    end subroutine correct_fixed_bond

!!$    subroutine evolve_cps()
!!$      integer ::             ia, it, itcrspd, iaa
!!$      real(kind=DP),allocatable,dimension(:,:) :: cps_old !d(nfcatm,3)
!!$      real(kind=DP),allocatable,dimension(:,:) :: cps_cog !d(num_planes_atoms_are_fixed,3)
!!$      real(kind=DP),allocatable,dimension(:,:) :: cps_cog_perpend !d(num_planes_atoms_are_fixed,3)
!!$      real(kind=DP) :: denom, pdot
!!$
!!$      if(constraint_type == COG_FIX_L &
!!$           & .or. constraint_type == COG_FIX .or. constraint_type == COG_CNTR) then
!!$         allocate(cps_old(nfcatm,3))
!!$         do it = 1, nfcatm
!!$            ia = ia_cnst(it)
!!$            cps_old(it,1:3) = cps(ia,1:3)
!!$         end do
!!$      end if
!!$
!!$      cps = cps + dtio*cpd_l
!!$
!!$      if(constraint_type == COG_FIX_L &
!!$           & .or. constraint_type == COG_FIX .or. constraint_type == COG_CNTR) then
!!$         allocate(cps_cog(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
!!$         allocate(cps_cog_perpend(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
!!$         do it = 1, nfcatm
!!$            ia = ia_cnst(it)
!!$            iaa = ipfixedplane(it)
!!$            cps_cog(iaa,1:3) = cps_cog(iaa,1:3) + amion(ityp(ia))*(cps(iaa,1:3)-cps_old(it,1:3))
!!$         end do
!!$         do it = 1, num_planes_atoms_are_fixed
!!$            denom = 0.d0
!!$            do itcrspd = 1, nfcatm
!!$               ia = ia_cnst(itcrspd)
!!$               if(ipfixedplane(itcrspd) == it) then
!!$                  denom = denom+amion(ityp(ia))
!!$                  iaa = itcrspd
!!$               end if
!!$            end do
!!$            cps_cog(it,1:3) = cps_cog(it,1:3)/denom
!!$            pdot = dot_product(fcvect(iaa,1:3),cps_cog(it,1:3))
!!$            cps_cog_perpend(it,1:3) = pdot*fcvect(iaa,1:3)   ! perpendicular components
!!$            if(iprimd >= 1) then
!!$               write(6,'(" !!f denom = ",f8.4)') denom
!!$            end if
!!$         end do
!!$         do it = 1, nfcatm
!!$            ia = ia_cnst(it)
!!$            itcrspd = ipfixedplane(it)
!!$            cps(ia,1:3) = cps(ia,1:3) - cps_cog_perpend(itcrspd,1:3)
!!$            if(iprimd >= 1) then
!!$               write(6,'(" !!f cps(        ",i3,") = ",3f8.4)') ia, cps(ia,1:3)
!!$               write(6,'(" !!f cps-cps_old(",i3,") = ",3f8.4)') ia, cps(ia,1:3)-cps_old(it,1:3)
!!$            end if
!!$         end do
!!$         deallocate(cps_old)
!!$         deallocate(cps_cog_perpend)
!!$         deallocate(cps_cog)
!!$      end if
!!$    end subroutine evolve_cps


  end subroutine m_IS_md

  subroutine check_constraint(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l

    integer :: ia, ib,ia1,ia2, i, ip, ifc, ifchit

    real(kind=DP) :: pdot
    real(kind=DP), allocatable, dimension(:,:) :: fcvect_t !d (num_planes_atoms_are_fixed,3)

    if(num_planes_atoms_are_fixed_cog+num_planes_atoms_are_fixed_rb>=1) &
         & allocate(fcvect_t(max(num_planes_atoms_are_fixed_cog, num_planes_atoms_are_fixed_rb),3))

    if(nfcatm>=1) then
       tmass = 0.d0
       do i = 1, nfcatm
          ia = ia_cnst(i)
          ip = ipfixedplane(i)
          tmass(ip) = tmass(ip) + amion(ityp(ia))
       end do
       do ip = 1, num_planes_atoms_are_fixed
          if(tmass(ip) >= (DELTA10-SmallestPositiveNumber) ) then
             rtmass(ip) = 1.d0/tmass(ip)
          else
             rtmass(ip) = 0.d0
          end if
       end do
    end if

    if(nfcatm_cog >= 1) then
       tmass_cog = 0.d0
       fcg_cog = 0.d0
!!$       ic_cog = 0
       do i = 1, nfcatm_cog
!!$          if(relax_in_fixedplane(i) == NO) cycle
          ia = ia_cnst(i)
          ip = icount_of_ipfixedplane(i)
          fcvect_t(ip,1:3) = fcvect(i,1:3)
          tmass_cog(ip) = tmass_cog(ip) + amion(ityp(ia))
          fcg_cog(ip,1:3) = fcg_cog(ip,1:3) + forc_l(ia,1:3)
       end do
       do ip = 1, num_planes_atoms_are_fixed_cog
          pdot = dot_product(fcvect_t(ip,1:3),fcg_cog(ip,1:3))
          fcg_mdfy_cog(ip,1:3) = fcg_cog(ip,1:3) - pdot*fcvect_t(ip,1:3)
          if(tmass_cog(ip) >= (DELTA10-SmallestPositiveNumber) ) then
             rtmass_cog(ip) = 1.d0/tmass_cog(ip)
          else
             rtmass_cog(ip) = 0.d0
          end if
       end do
    end if
    if(nfcatm_rb >=1) then
       tmass_rb = 0.d0
       fcg_rb = 0.d0
!!$       ic_rb = 0
       do i = nfcatm_cog+1, nfcatm_cog+nfcatm_rb
!!$          if(relax_in_fixedplane(i) == NO) cycle
          ia = ia_cnst(i)
          ip = icount_of_ipfixedplane(i)
          fcvect_t(ip,1:3) = fcvect(i,1:3)
          tmass_rb(ip) = tmass_rb(ip) + amion(ityp(ia))
          fcg_rb(ip,1:3) = fcg_rb(ip,1:3) + forc_l(ia,1:3)

          if(iprimd >= DEBUGPRINTLEVEL) then
             write(nfout,'("!check_constraint forc_l(",i5,",1:3) = ",3f15.8," fcg_rb(",i3,",1:3) = ",3f15.8)') &
                  & ia, forc_l(ia,1:3), ip, fcg_rb(ip,1:3)
          end if
       end do
       do ip = 1, num_planes_atoms_are_fixed_rb
          if(relax_rigid_body_cog==0) then
             pdot = dot_product(fcvect_t(ip,1:3),fcg_rb(ip,1:3))
             fcg_mdfy_rb(ip,1:3) = fcg_rb(ip,1:3) - pdot*fcvect_t(ip,1:3)
          else
             fcg_mdfy_rb(ip,1:3) = fcg_rb(ip,1:3)
          end if
          if(tmass_rb(ip) >= (DELTA10-SmallestPositiveNumber) ) then
             rtmass_rb(ip) = 1.d0/tmass_rb(ip)
          else
             rtmass_rb(ip) = 0.d0
          end if
       end do

    end if

    if(iprimd>=1) then
       do ip = 1, num_planes_atoms_are_fixed_cog
          write(nfout, &
          & '("<<check_constraint>> ip = ",i3," fcg_mdfy_cog = ",3f12.8, " tmass_cog = ",f12.4, " rtmass_cog = ",f12.8)') &
          & ip, fcg_mdfy_cog(ip,1:3),tmass_cog(ip), rtmass_cog(ip)
       end do
       do ip = 1, num_planes_atoms_are_fixed_rb
          write(nfout, &
          & '("<<check_constraint>> ip = ",i3," fcg_mdfy_rb  = ",3f12.8, " tmass_rb  = ",f12.4, " rtmass_rb = ",f12.8)') &
          & ip, fcg_mdfy_rb(ip,1:3),tmass_rb(ip), rtmass_rb(ip)
       end do
    end if

    ! --- bondlength check ---
    if(constraint_type == BONDLENGTH_FIX) then
       do ib = 1, num_fixed_bonds
          ia1 = bondlength_fix_set(1,ib)
          ia2 = bondlength_fix_set(2,ib)
          ia = ib*2 - 1
          fcvect(ia,1:3) = cps(ia1,1:3) - cps(ia2,1:3)
          fcvect(ia,4)   = dsqrt(fcvect(ia,1)**2 + fcvect(ia,2)**2 + fcvect(ia,3)**2)
          fcvect(ia+1,1:3) = -fcvect(ia,1:3)
          fcvect(ia+1,4) = fcvect(ia,4)
          if(iprimd >= 3) then
             write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<check_constraint>>")') &
                  & ia,(fcvect(ia,1:4))
             write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<check_constraint>>")') &
                  & ia+1,(fcvect(ia+1,1:4))
          end if
       end do
    end if
    if(num_planes_atoms_are_fixed_cog+num_planes_atoms_are_fixed_rb>=1) deallocate(fcvect_t)
  end subroutine check_constraint

  subroutine check_constraint_cog_or_rb(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l

    integer :: ic_cog, ia, iaa, ib,ia1,ia2, ip, ifchit, ifc, i, j, ic

    real(kind=DP) :: pdot
    real(kind=DP), allocatable, dimension(:,:) :: fcvect_t ! d(num_planes_atoms_are_fixed,3)


    if(num_planes_atoms_are_fixed_cog+num_planes_atoms_are_fixed_rb>=1) &
         & allocate(fcvect_t(max(num_planes_atoms_are_fixed_cog,num_planes_atoms_are_fixed_rb),1:3))

    if(nfcatm_cog >= 1) then
       tmass_cog = 0.d0
       fcg_cog = 0.d0
       ic_cog = 0
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L) then
             ip = icount_of_ipfixedplane(icnst_a(ia))
             ic_cog = ic_cog + 1
             tmass_cog(ip) = tmass_cog(ip) + amion(ityp(ia))
             fcg_cog(ip,1:3) = fcg_cog(ip,1:3) + forc_l(ia,1:3)
             fcvect_t(ip,1:3) = fcvect(icnst_a(ia),1:3)
          end if
       end do
       if(ic_cog <= 1) call phase_error_with_msg(nfout, ' #ic_cog is not enough',__LINE__,__FILE__)

       do ip = 1, num_planes_atoms_are_fixed_cog
          pdot = dot_product(fcvect_t(ip,1:3),fcg_cog(ip,1:3))
          fcg_mdfy_cog(ip,1:3) = fcg_cog(ip,1:3) - pdot*fcvect_t(ip,1:3)
          if(tmass_cog(ip) >= (DELTA10-SmallestPositiveNumber) ) then
             rtmass_cog(ip) = 1.d0/tmass_cog(ip)
          else
             rtmass_cog(ip) = 0.d0
          end if
       end do
    end if
    if(nfcatm_rb >=1) then
       tmass_rb = 0.d0
       fcg_rb = 0.d0
!!$       ic_rb = 0
       do i = nfcatm_cog+1, nfcatm_cog+nfcatm_rb
!!$          if(relax_in_fixedplane(i) == NO) cycle
          ia = ia_cnst(i)
          ip = icount_of_ipfixedplane(i)
          fcvect_t(ip,1:3) = fcvect(i,1:3)
          tmass_rb(ip) = tmass_rb(ip) + amion(ityp(ia))
          fcg_rb(ip,1:3) = fcg_rb(ip,1:3) + forc_l(ia,1:3)
       end do
       do ip = 1, num_planes_atoms_are_fixed_rb
          if(relax_rigid_body_cog==0) then
             pdot = dot_product(fcvect_t(ip,1:3),fcg_rb(ip,1:3))
             fcg_mdfy_rb(ip,1:3) = fcg_rb(ip,1:3) - pdot*fcvect_t(ip,1:3)
          else
             fcg_mdfy_rb(ip,1:3) = fcg_rb(ip,1:3)
          end if
          if(tmass_rb(ip) >= (DELTA10-SmallestPositiveNumber) ) then
             rtmass_rb(ip) = 1.d0/tmass_rb(ip)
          else
             rtmass_rb(ip) = 0.d0
          end if
       end do

    end if
    if(iprimd>=1) then
       do ip = 1, num_planes_atoms_are_fixed_cog
          write(nfout, &
          & '("<<check_constraint>> ip = ",i3," fcg_mdfy_cog = ",3f12.8, " tmass_cog = ",f12.4, " rtmass_cog = ",f12.8)') &
          & ip, fcg_mdfy_cog(ip,1:3),tmass_cog(ip), rtmass_cog(ip)
       end do
       do ip = 1, num_planes_atoms_are_fixed_rb
          write(nfout,'("<<check_constraint>> ip = ",i3," fcg_mdfy_rb  = ",3f12.8, " tmass_rb  = ",f12.4, " rtmass_rb = ",f12.8)') &
               & ip, fcg_mdfy_rb(ip,1:3),tmass_rb(ip), rtmass_rb(ip)
       end do
    end if

    ! --- bondlength check ---
    if(constraint_type == BONDLENGTH_FIX) then
       do ib = 1, num_fixed_bonds
          ia1 = bondlength_fix_set(1,ib)
          ia2 = bondlength_fix_set(2,ib)
          ia = ib*2 - 1
          fcvect(ia,1:3) = cps(ia1,1:3) - cps(ia2,1:3)
          fcvect(ia,4)   = dsqrt(fcvect(ia,1)**2 + fcvect(ia,2)**2 + fcvect(ia,3)**2)
          fcvect(ia+1,1:3) = -fcvect(ia,1:3)
          fcvect(ia+1,4) = fcvect(ia,4)
          if(iprimd >= 3) then
             write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<check_constraint>>")') &
                  & ia,(fcvect(ia,1:4))
             write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<check_constraint>>")') &
                  & ia+1,(fcvect(ia+1,1:4))
          end if
       end do
    end if
    if(num_planes_atoms_are_fixed_cog+num_planes_atoms_are_fixed_rb>=1) deallocate(fcvect_t)

  end subroutine check_constraint_cog_or_rb

  real(kind=DP) function get_fi_coefficient(force_t_norm)
    real(kind=DP), intent(in) :: force_t_norm
    if(force_t_norm < 0.008d0) then
       get_fi_coefficient = 3.5d0
    else if (force_t_norm < 0.02d0) then
       get_fi_coefficient = 2.5d0
    else if (force_t_norm < 0.05d0) then
       get_fi_coefficient = 1.5d0
    else
       get_fi_coefficient = 0.8d0
    end if
  end function get_fi_coefficient

  subroutine move_atoms_normal_to_plane()
    integer :: ia, icatm, i, ip
    icatm = 0
    do ia = 1, natm
       if(imdtyp(ia)==FIX_IN_A_PLANE .or. imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==COG_FIX .or. &
          imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
          do i = 1, nfcatm
             if(ia_cnst(i) == ia) then
                ip = i
                exit
             end if
          end do
          if(ip<=0 .or. nfcatm< ip) then
             if(iprimd >= 1) then
                write(nfout,'(" ip = ",i4," ia = ",i4)') ip,ia
                write(nfout,'(" ip is out of range <<move_atoms_normal_to_plane>>")')
                call flush(nfout)
             end if
             call phase_error_with_msg(nfout, 'ip is out of range <<move_atoms_normal_to_plane>>', &
             __LINE__,__FILE__)
          end if

          if(iprimd >= DEBUGPRINTLEVEL) then
             write(nfout,'(" ia = ",i8, " cps_diff = ",3f12.6)') ia, (fcvect(ip,i)*fcvect(ip,4),i=1,3)
             write(nfout,'(" fcvect = ",4f12.6)') (fcvect(ip,i),i=1,4)
          end if
          if((fcvect(ip,1)**2 + fcvect(ip,2)**2 + fcvect(ip,3)**2) < 1.d-9) then
             if(iprimd >= 1) write(nfout,'(" fcvect(",i3,",1:3) is too small <<move_atoms_normal_to_plane>>")') ip
             call phase_error_with_msg(nfout,' fcvect is too small <<move_atoms_normal_to_plane>>', &
             __LINE__,__FILE__)
          end if
          if(fcvect(ip,4)**2 < 1.d-9) then
             if(iprimd >= 1) write(nfout,'(" fcvect(",i3,",4) increment is too small <<move_atoms_normal_to_plane>>")') ip
             call phase_error_with_msg(nfout, ' fcvect increment is too small <<move_atoms_normal_to_plane>>', &
             __LINE__,__FILE__)
          end if

          cps(ia,1:3) = cps(ia,1:3) + fcvect(ip,1:3)*fcvect(ip,4)
       end if
    end do
    moved_distance_of_fixed_plane = moved_distance_of_fixed_plane + fcvect(1,4)
    cpd_l = 0.d0
  end subroutine move_atoms_normal_to_plane

  subroutine check_if_bondlength_fix_exist(ib)
    integer, intent(out) :: ib
    integer i
    ib = 0
    do i = 1, natm
       if(imdtyp(i) == BONDLENGTH_FIX) ib = ib + 1
    end do
  end subroutine check_if_bondlength_fix_exist

  subroutine md1_alloc(mdalg)
    integer, intent(in) :: mdalg
    allocate(cpd_old(natm,3));cpd_old=0.d0
    if(mdalg == VERLET .and. &
         &(nbztyp == HEXAGONAL .or. nbztyp == ORTHORHOMBIC .or. &
         & nbztyp >= HEX1fold)) then
       allocate(ipcpd(natm2))
    else
       allocate(ipcpd(natm))
    end if
  end subroutine md1_alloc

  subroutine md1_dealloc()
    deallocate(ipcpd)
    deallocate(cpd_old)
  end subroutine md1_dealloc

  subroutine md1_alloc2()
    integer :: ia, ifcatm

    allocate(ipcpd(natm)); ipcpd = 0
    ifcatm = 0
    do ia = 1, natm
       if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L &
            .or. imdtyp(ia) == FIX_IN_A_PLANE) then
          ifcatm = ifcatm + 1
          ipcpd(ia) = ifcatm
       end if
    end do
  end subroutine md1_alloc2

  subroutine md1_dealloc2()
    deallocate(ipcpd)
  end subroutine md1_dealloc2

  subroutine quench_velocities(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer       :: ia, ip
    real(kind=DP) :: cpd_parallel(3), cpd_vertical(3), f, fv, fcpd, v(3), cpd_av(3)
    !  Revised according to an indication by Dr. T. Yamamoto.
    !            T. Yamasaki,   Nov. 2003
    !  Revised according to an indication by Usami-san @adv
    !      An adopted quenching algorithm depends on the nopr value
    !            T. Yamasaki,   Sep. 2006
    !
!!$    real(kind=DP) :: v
!!$    do j = 1, 3
!!$!xocl spread do/ind_natm
!!$       do ia = 1, natm
!!$          if(imdtyp(ia) == BONDLENGTH_FIX) cycle
!!$          v = (cpd_old(ia,j) + cpd_l(ia,j))*0.5
!!$          if(dabs(v) > SmallestPositiveNumber*1.d5 .and. v*forc_l(ia,j) < 0.d0) then
!!$             if(printable) write(nfout,'(" quenched atom = ",i6," --- <<quench_velocities>>")') ia
!!$             cpd_l(ia,j) = 0.d0
!!$          end if
!!$       end do
!!$!xocl end spread
!!$    end do
    if(nopr <= 1) then
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == BONDLENGTH_FIX) cycle
          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) cycle
          v(1:3) = (cpd_old(ia,1:3) + cpd_l(ia,1:3))*0.5
          f = forc_l(ia,1)**2 + forc_l(ia,2)**2 + forc_l(ia,3)**2
!!$          if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE) then
!!$          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE .or. imdtyp(ia) == COG_FIX_L) then
          if(imdtyp(ia) == COG_FIX_L) then
             ip = icount_of_ipfixedplane(icnst_a(ia))
             fv = dot_product(fcg_mdfy_cog(ip,1:3), v(1:3))
!!$          else if(imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
!!$             ip = icount_of_ipfixedplane(icnst_a(ia))
!!$             fv = dot_product(fcg_mdfy_rb(ip,1:3),v)
          else
             fv = dot_product(forc_l(ia,1:3), v)
          end if
          fcpd = dot_product(forc_l(ia,1:3),cpd_l(ia,1:3))
          if(f > SmallestPositiveNumber*1.d5 ) then
             if(imdtyp(ia) == COG_FIX_L) then
                cpd_parallel(1:3) = fcpd/f * fcg_mdfy_cog(ip,1:3)
             else if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                cpd_parallel(1:3) = fcpd/f * fcg_mdfy_rb(ip,1:3)
             else
                cpd_parallel(1:3) = fcpd/f * forc_l(ia,1:3)
             end if
             cpd_vertical(1:3) = cpd_l(ia,1:3) - cpd_parallel(1:3)
             if(fv < 0.d0) then
                cpd_parallel(1:3) = 0.d0
                cpd_l(ia,1:3) = cpd_vertical(1:3)
                if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" quenched atom = ",i6," ---<<quench_velocities>>")') ia
             end if
          else
             cpd_l(ia,1:3) = 0.d0
          end if
       end do
!xocl end spread
    else
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == BONDLENGTH_FIX) cycle
          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) cycle
          v(1:3) = (cpd_old(ia,1:3) + cpd_l(ia,1:3))*0.5
          f = forc_l(ia,1)**2 + forc_l(ia,2)**2 + forc_l(ia,3)**2
          if(imdtyp(ia) == COG_FIX_L) then
             ip = icount_of_ipfixedplane(icnst_a(ia))
             fv = dot_product(fcg_mdfy_cog(ip,1:3), v(1:3))
          else if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
             ip = icount_of_ipfixedplane(icnst_a(ia))
             fv = dot_product(fcg_mdfy_rb(ip,1:3),v)
          else
             fv = dot_product(forc_l(ia,1:3), v)
          end if
          fcpd = dot_product(forc_l(ia,1:3),cpd_l(ia,1:3))
          if(f > SmallestPositiveNumber*1.d5 ) then
             if(imdtyp(ia) == COG_FIX_L) then
                cpd_parallel(1:3) = fcpd/f * fcg_mdfy_cog(ip,1:3)
!!$             else if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
!!$                cpd_parallel(1:3) = fcpd/f * fcg_mdfy_rb(ip,1:3)
             else
                cpd_parallel(1:3) = fcpd/f * forc_l(ia,1:3)
             end if
             cpd_vertical(1:3) = cpd_l(ia,1:3) - cpd_parallel(1:3)
             if(fv < 0.d0) then
                cpd_parallel(1:3) = 0.d0
                cpd_l(ia,1:3) = cpd_vertical(1:3)
                if(iprimd>=1) write(nfout,'(" quenched atom = ",i6," ---<<quench_velocities>>")') ia
             end if
          else
             cpd_l(ia,1:3) = 0.d0
          end if
       end do
!xocl end spread
    end if

  end subroutine quench_velocities

  subroutine quench_velocities_cog(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer       :: ia, ia_t
    real(kind=DP) :: cpd_parallel(3), cpd_vertical(3), f, fv, fcpd, v(3)
    !  Revised according to an indication by Dr. T. Yamamoto.
    !            T. Yamasaki,   Nov. 2003
    !  Revised according to an indication by Usami-san @adv
    !      An adopted quenching algorithm depends on the nopr value
    !            T. Yamasaki,   Sep. 2006
    !
!!$    real(kind=DP) :: v
!!$    do j = 1, 3
!!$!xocl spread do/ind_natm
!!$       do ia = 1, natm
!!$          if(imdtyp(ia) == BONDLENGTH_FIX) cycle
!!$          v = (cpd_old(ia,j) + cpd_l(ia,j))*0.5
!!$          if(dabs(v) > SmallestPositiveNumber*1.d5 .and. v*forc_l(ia,j) < 0.d0) then
!!$             if(printable) write(nfout,'(" quenched atom = ",i6," --- <<quench_velocities>>")') ia
!!$             cpd_l(ia,j) = 0.d0
!!$          end if
!!$       end do
!!$!xocl end spread
!!$    end do
    if(nopr <= 1) then
       ia_t = 0
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L &
               .or. imdtyp(ia) == FIX_IN_A_PLANE) then
             ia_t = ia_t + 1
             v(1:3) = (cpd_old(ia_t,1:3) + cpd_l(ia,1:3))*0.5
             fv = dot_product(forc_l(ia,1:3), v)
             if(fv < 0.d0) then
                if(iprimd>=1) write(nfout,'(" quenched atom = ",i6," ---<<quench_velocities>>")') ia
                cpd_l(ia,1:3) = 0.d0
             end if
          end if
       end do
!xocl end spread
    else
       ia_t = 0
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L &
               .or. imdtyp(ia) == FIX_IN_A_PLANE) then
             ia_t = ia_t + 1
             v(1:3) = (cpd_old(ia_t,1:3) + cpd_l(ia,1:3))*0.5
             f = forc_l(ia,1)**2 + forc_l(ia,2)**2 + forc_l(ia,3)**2
             fv = dot_product(forc_l(ia,1:3), v)
             fcpd = dot_product(forc_l(ia,1:3),cpd_l(ia,1:3))
             if(f > SmallestPositiveNumber*1.d5 ) then
                cpd_parallel(1:3) = fcpd/f * forc_l(ia,1:3)
                cpd_vertical(1:3) = cpd_l(ia,1:3) - cpd_parallel(1:3)
                if(fv < 0.d0) then
                   cpd_parallel(1:3) = 0.d0
                   cpd_l(ia,1:3) = cpd_vertical(1:3)
                   if(iprimd>=1) write(nfout,'(" quenched atom = ",i6," ---<<quench_velocities>>")') ia
                end if
             else
                cpd_l(ia,1:3) = 0.d0
             end if
          end if
       end do
!xocl end spread
    end if

  end subroutine quench_velocities_cog

  subroutine quench_velocities_rb(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer :: ia, ip, imd
    real(kind=DP) :: amass, pdot
    real(kind=DP), allocatable,dimension(:,:) :: dcg, dcg_n, fcvect_t
    real(kind=DP), allocatable,dimension(:)   :: denom

    if(relax_rigid_body_cog==On) return
    if(num_planes_atoms_are_fixed_rb<=0) return
!!$       if(num_planes_atoms_are_fixed_cog >=1) allocate(dcg(num_planes_atoms_are_fixed_cog,3))
    allocate(dcg(num_planes_atoms_are_fixed_rb,3))
    allocate(dcg_n(num_planes_atoms_are_fixed_rb,3))
    allocate(fcvect_t(num_planes_atoms_are_fixed_rb,3))
    allocate(denom(num_planes_atoms_are_fixed_rb))
    dcg = 0.d0
    denom = 0.d0
    do ia = 1, natm
       if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
          amass = amion(ityp(ia))
          ip = icount_of_ipfixedplane(icnst_a(ia))
          fcvect_t(ip,1:3) = fcvect(icnst_a(ia),1:3)
          dcg(ip,:) = dcg(ip,:) + amass*cpd_l(ia,:)
          denom(ip) = denom(ip) + amass
       end if
    end do
    do ip = 1, num_planes_atoms_are_fixed_rb
       dcg(ip,:) = dcg(ip,:)/denom(ip)
       if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" !D dcg(",i3," )   = ",3f12.8,") <<quench_velocites_rb>>")') &
                                         ip,dcg(ip,1:3)
       pdot = dot_product(fcvect_t(ip,1:3),dcg(ip,1:3))
       dcg_n(ip,1:3) = pdot*fcvect_t(ip,1:3)
       if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" !D dcg_n(",i3," ) = ",3f12.8,") <<quench_velocites_rb>>")') &
                                         ip,dcg_n(ip,1:3)
    end do

    do ia = 1, natm
       if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
          ip = icount_of_ipfixedplane(icnst_a(ia))
!!$             pdot = dot_product(fcvect(icnst_a(ia),1:3),cpd_l(ia,1:3))
!!$             dcg_n(ip,:) = pdot*fcvect(icnst_a(ia),1:3)
!!$             cpd_l(ia,:) = cpd_l(ia,:) - dcg_n(ip,:)
          cpd_l(ia,:) = dcg(ip,:) - dcg_n(ip,:)
          if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" !D dcg_n(",i3," ) = ",3f12.8,") <<quench_velocites_rb>>")') &
                                         ip,dcg_n(ip,1:3)
       end if
    end do

    do ia = 1, natm
       imd = imdtyp(ia)
       if(imd == RIGID_BODY_FIX_IN_A_PLANE) then
          ip = icnst_a(ia)
          if(relax_in_fixedplane(ip) == NO) cpd_l(ia,:) = 0.d0
       end if
    end do
    deallocate(denom,fcvect_t,dcg_n,dcg)
  end subroutine quench_velocities_rb

  subroutine correct_cog_andor_rb_motion()
    integer :: ia, ip, imd
    real(kind=DP) :: amass, pdot, f, fv
    real(kind=DP), allocatable,dimension(:,:) :: dcg, dcg_n, fcvect_t
    real(kind=DP), allocatable,dimension(:)   :: denom
    real(kind=DP) :: cpd_parallel(3), cpd_vertical(3), fcpd, v(3), cpd_av(3)

    if(iprimd >= DEBUGPRINTLEVEL) then
       write(nfout,'(" <<correct_cog_andor_rb_motion>>")')
       if(nfcatm>=1) then
          do ia = 1, natm
             if(imdtyp(ia) == COG_FIX_L .or. imdtyp(ia) == COG_FIX .or. imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) &
                  & write(nfout,'(" cpd_l(",i4,") = ",3f13.8," initial, <<correct_cog_andor_rb_motion>>")') ia, cpd_l(ia,:)
          end do
       end if
    end if

    if(num_planes_atoms_are_fixed>=1) then
!!$       if(num_planes_atoms_are_fixed_cog >=1) allocate(dcg(num_planes_atoms_are_fixed_cog,3))
       allocate(dcg(max(num_planes_atoms_are_fixed_cog,num_planes_atoms_are_fixed_rb),3))
       allocate(dcg_n(max(num_planes_atoms_are_fixed_cog,num_planes_atoms_are_fixed_rb),3))
       allocate(fcvect_t(max(num_planes_atoms_are_fixed_cog,num_planes_atoms_are_fixed_rb),3))
       allocate(denom(max(num_planes_atoms_are_fixed_cog,num_planes_atoms_are_fixed_rb)))
       dcg = 0.d0
       denom = 0.d0
    end if
    do ia = 1, natm
       if(imdtyp(ia) == COG_FIX_L .or. imdtyp(ia) == COG_FIX) then
          amass = amion(ityp(ia))
          ip = icount_of_ipfixedplane(icnst_a(ia))
          fcvect_t(ip,1:3) = fcvect(icnst_a(ia),1:3)
          dcg(ip,:) = dcg(ip,:) + amass*cpd_l(ia,:)
          denom(ip) = denom(ip) + amass
          if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" cpd_l(",i4,")=",3f13.6," <<correct_cog_andor_rb_motion>>")') ia,cpd_l(ia,:)
       end if
    end do
    if(num_planes_atoms_are_fixed_cog >= 1) then
       do ip = 1, num_planes_atoms_are_fixed_cog
          dcg(ip,:) = dcg(ip,:)/denom(ip)
          if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" !D dcg(",i3," )   = ",3f12.8,") <<correct_cog_andor_rb_motion>>")') &
                                            ip,dcg(ip,1:3)
          pdot = dot_product(fcvect_t(ip,1:3),dcg(ip,1:3))
          dcg_n(ip,1:3) = pdot*fcvect_t(ip,1:3)
          if(iprimd>=DEBUGPRINTLEVEL) write(nfout,'(" !D dcg_n(",i3," ) = ",3f12.8,") <<correct_cog_andor_rb_motion>>")') &
                                            ip,dcg_n(ip,1:3)
       end do
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX_L .or. imdtyp(ia) == COG_FIX) then
             ip = icount_of_ipfixedplane(icnst_a(ia))
             cpd_l(ia,:) = cpd_l(ia,:) - dcg_n(ip,:)
          end if
       end do
    end if

    if(relax_rigid_body_cog==OFF) then
       if(num_planes_atoms_are_fixed_rb>=1) then
          dcg = 0.d0
          denom = 0.d0
          do ia = 1, natm
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                amass = amion(ityp(ia))
                ip = icount_of_ipfixedplane(icnst_a(ia))
                fcvect_t(ip,1:3) = fcvect(icnst_a(ia),1:3)
                dcg(ip,:) = dcg(ip,:) + amass*cpd_l(ia,:)
                denom(ip) = denom(ip) + amass
             end if
          end do
          do ip = 1, num_planes_atoms_are_fixed_rb
             dcg(ip,:) = dcg(ip,:)/denom(ip)
             if(iprimd>=DEBUGPRINTLEVEL) &
                  & write(nfout,'(" !D dcg(",i3," )   = ",3f12.8,") <<correct_cog_andor_rb_motion>>")') ip,dcg(ip,1:3)
             pdot = dot_product(fcvect_t(ip,1:3),dcg(ip,1:3))
             dcg_n(ip,1:3) = pdot*fcvect_t(ip,1:3)
             if(iprimd>=DEBUGPRINTLEVEL) &
                  & write(nfout,'(" !D dcg_n(",i3," ) = ",3f12.8,") <<correct_cog_andor_rb_motion>>")') ip,dcg_n(ip,1:3)
          end do

          do ia = 1, natm
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                ip = icount_of_ipfixedplane(icnst_a(ia))
!!$             pdot = dot_product(fcvect(icnst_a(ia),1:3),cpd_l(ia,1:3))
!!$             dcg_n(ip,:) = pdot*fcvect(icnst_a(ia),1:3)
!!$             cpd_l(ia,:) = cpd_l(ia,:) - dcg_n(ip,:)
                cpd_l(ia,:) = dcg(ip,:) - dcg_n(ip,:)
                if(iprimd>=DEBUGPRINTLEVEL) &
                  & write(nfout,'(" !D dcg_n(",i3," ) = ",3f12.8,") <<correct_cog_andor_rb_motion>>")') ip,dcg_n(ip,1:3)
            end if
          end do
       end if
    else if(relax_rigid_body_cog==ON) then
       if(num_planes_atoms_are_fixed_rb>=1) then
!!$          do ia = 1, natm
!!$             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
!!$                amass = amion(ityp(ia))
!!$                ip = icount_of_ipfixedplane(icnst_a(ia))
!!$                dcg(ip,:) = dcg(ip,:) + amass*cpd_l(ia,:)
!!$                denom(ip) = denom(ip) + amass
!!$             end if
!!$          end do

          dcg = 0.d0
          denom = 0.d0
          do ia = 1, natm
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                amass = amion(ityp(ia))
                ip = icount_of_ipfixedplane(icnst_a(ia))
                dcg(ip,:) = dcg(ip,:) + amass*cpd_l(ia,:)
                denom(ip) = denom(ip) + amass
             end if
          end do
          do ip = 1, num_planes_atoms_are_fixed_rb
             dcg(ip,:) = dcg(ip,:)/denom(ip)
             if(iprimd>=DEBUGPRINTLEVEL) &
                  & write(nfout,'(" !D dcg and fcg_mdfy_rb(",i3," )  = ",6f12.8,") <<correct_cog_andor_rb_motion>>")') &
                  & ip,dcg(ip,1:3), fcg_mdfy_rb(num_planes_atoms_are_fixed_cog+ip,1:3)
          end do

          do ia = 1, natm
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                ip = icount_of_ipfixedplane(icnst_a(ia))
                f = fcg_mdfy_rb(ip,1)**2 + fcg_mdfy_rb(ip,2)**2 + fcg_mdfy_rb(ip,3)**2
                fv = dot_product(fcg_mdfy_rb(ip,1:3), dcg(ip,1:3))
                if(f > SmallestPositiveNumber*1.d5 ) then
                   cpd_parallel(1:3) = fv/f * fcg_mdfy_rb(ip,1:3)
                   cpd_vertical(1:3) = dcg(ip,1:3) - cpd_parallel(1:3)
                   if(fv < 0.d0) then
                      cpd_parallel(1:3) = 0.d0
                      cpd_l(ia,1:3) = cpd_vertical(1:3)
                   else
                      cpd_l(ia,:) = dcg(ip,:)
                   end if
                   if(iprimd>=DEBUGPRINTLEVEL)then
                      if(fv<0.d0) then
                         write(nfout,'(" quenched atom = ",i6," fv<0.d0  ---<<correct_cog_andor_rb_motion>>")') ia
                      else
                         write(nfout,'(" quenched atom = ",i6," fv>=0.d0 ---<<correct_cog_andor_rb_motion>>")') ia
                      end if
                   end if
                else
                   if(iprimd>=DEBUGPRINTLEVEL) &
                        & write(nfout,'(" quenched atom = ",i6," f=0  ---<<correct_cog_andor_rb_motion>>")') ia
                   cpd_l(ia,1:3) = 0.d0
                end if
             end if
          end do

!!$          do ia = 1, natm
!!$             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
!!$                ip = icount_of_ipfixedplane(icnst_a(ia))
!!$                cpd_l(ia,:) = dcg(ip,:)
!!$             end if
!!$          end do
       end if
    endif

    do ia = 1, natm
       imd = imdtyp(ia)
       if(imd == RIGID_BODY_FIX_IN_A_PLANE .or. imd == COG_FIX_L .or. imd == COG_FIX) then
!!$          ip = icnst_a(ia)
!!$          if(relax_in_fixedplane(ip) == NO) cpd_l(ia,:) = 0.d0
       end if
    end do

    if(nfcatm>=1 .and. iprimd >=DEBUGPRINTLEVEL ) then
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX_L .or. imdtyp(ia) == COG_FIX .or. imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) &
              & write(nfout,'(" cpd_l(",i4,") = ",3f13.8," corrected, <<correct_cog_andor_rb_motion>>")') ia, cpd_l(ia,:)
       end do
    end if

    if(allocated(denom)) deallocate(denom)
    if(allocated(fcvect_t)) deallocate(fcvect_t)
    if(allocated(dcg_n)) deallocate(dcg_n)
    if(allocated(dcg))   deallocate(dcg)
  end subroutine correct_cog_andor_rb_motion

  subroutine get_ekina
    integer ia
    ekina = 0.d0
!xocl spread do/ind_natm
    do ia = 1, natm
       if(imdtyp(ia) == BONDLENGTH_FIX) cycle
       ekina = ekina + ( (cpd_old(ia,1) + cpd_l(ia,1))**2 &
            &           +(cpd_old(ia,2) + cpd_l(ia,2))**2 &
            &           +(cpd_old(ia,3) + cpd_l(ia,3))**2)*0.25 & ! v*v
            &          * amion(ityp(ia))*iwei(ia) * 0.5     ! * mass * 1/2
    end do
!!$!xocl end spread sum(ekina)
    if(nrigid_bodies>0) ekina = ekina + m_IS_rb_kinetic_energies()
  end subroutine get_ekina

  subroutine modify_forc_fi(forc_l)
    real(kind=DP), intent(inout), dimension(natm,3) :: forc_l
    real(kind=DP), allocatable, dimension(:,:) :: force_t_para, force_t_vert
    real(kind=DP) :: f_dot_v, vec_norm2, force_t_norm
    integer :: ia

    do ia = 1, natm
       if(imdtyp(ia) == FIX) then
          forc_l(ia,1:3) = 0.d0
       end if
    end do

    allocate(force_t_para(natm,3), force_t_vert(natm,3))

    f_dot_v = sum( forc_l(1:natm,1:3)*normal_hypervector(1:natm,1:3,CARTS) )
    vec_norm2 = sum( normal_hypervector(1:natm,1:3,CARTS)**2 )
    force_t_para = 0.d0; force_t_vert = 0.d0
    if( vec_norm2 > DELTA10) then
       force_t_para(:,:)=f_dot_v/vec_norm2*normal_hypervector(:,:,CARTS)
       force_t_vert(:,:)=forc_l(:,:)-force_t_para(:,:)
       f_dot_v= sum( force_t_para(:,:)*force_t_vert(:,:) )
       if(iprimd >= 1) then
          write(nfout,*) 'opt_mode is FI'
          write(nfout,'(a19,f20.10)') ' true_force =      ', force_t_norm
          write(nfout,'(a19,f20.10)') ' force_parallel =  ' &
               & , sum(force_t_para(:,:)*normal_hypervector(:,:,CARTS))/sqrt(vec_norm2)
          write(nfout,'(a19,f20.10)') ' force_vertical =  ' &
               & , sqrt(sum(force_t_vert(:,:)**2))
          write(nfout,'(a19,f20.10)') ' f_para x f_vert = ', f_dot_v
          if(iprimd >= 2) then
             write(nfout,'(" forc_l, force_t_vert, force_t_para")')
             do ia = 1, natm
                write(nfout,'(i5,9f10.4)') forc_l(ia,1:3), force_t_vert(ia,1:3),force_t_para(ia,1:3)
             end do
          end if
       end if
       forc_l(:,:)=force_t_vert(:,:)-force_t_para(:,:)
       if(iprimd >= 1) write(nfout,'(a19,f20.10)') 'modified force =  ', sqrt(sum(forc_l(:,:)**2))
    else
       if(iprimd >= 1) then
          write(nfout,*) 'vector is strange for FI or NEB !!'
          write(nfout,*) 'vector =', vec_norm2
          write(nfout,*) 'opt_mode is changed to SD'
          write(nfout,*) ' total force =', force_t_norm
       end if
    end if

    deallocate(force_t_vert, force_t_para)

  end subroutine modify_forc_fi

  subroutine modify_forc_hyperplane(forc_l)
    real(kind=DP), intent(inout), dimension(natm,3) :: forc_l
    real(kind=DP), allocatable, dimension(:,:) :: force_t_para, force_t_vert
    real(kind=DP) :: f_dot_v, vec_norm2, force_t_norm, fa
    integer :: ia

    do ia = 1, natm
       if(imdtyp(ia) == FIX) then
          forc_l(ia,1:3) = 0.d0
       end if
    end do

    allocate(force_t_para(natm,3), force_t_vert(natm,3))

    f_dot_v = sum( forc_l(1:natm,1:3)*normal_hypervector(1:natm,1:3,CARTS) )
    force_t_norm = sqrt(sum( forc_l(1:natm,1:3)**2))
    vec_norm2 = sum( normal_hypervector(1:natm,1:3,CARTS)**2 )
    force_t_para = 0.d0; force_t_vert = 0.d0
    if( vec_norm2 > DELTA10) then
       force_t_para(:,:)=f_dot_v/vec_norm2*normal_hypervector(:,:,CARTS)
       force_t_vert(:,:)=forc_l(:,:)-force_t_para(:,:)
       f_dot_v= sum( force_t_para(:,:)*force_t_vert(:,:) )
       if(iprimd >= 1) then
          write(nfout,*) 'opt_mode is fix_in_a_hyperplane'
          write(nfout,'(a19,f20.10)') ' true_force =      ', force_t_norm
          write(nfout,'(a19,f20.10)') ' force_parallel =  ' &
               & , sum(force_t_para(:,:)*normal_hypervector(:,:,CARTS))/sqrt(vec_norm2)
          write(nfout,'(a19,f20.10)') ' force_vertical =  ' &
               & , sqrt(sum(force_t_vert(:,:)**2))
          write(nfout,'(a19,f20.10)') ' f_para x f_vert = ', f_dot_v
! -->   T. Yamasaki 22 Aug 2008
          if(iprimd >= 2) then
             write(nfout,'(" forc_l, force_t_vert, force_t_para")')
             do ia = 1, natm
                write(nfout,'(i5,9f10.4)') ia, forc_l(ia,1:3), force_t_vert(ia,1:3),force_t_para(ia,1:3)
             end do
          end if
! <--
       end if
!!$       forc_l(:,:)=force_t_vert(:,:)-force_t_para(:,:)
       forc_l(:,:)=force_t_vert(:,:)
! -->   T. Yamasaki 18 July 2008
       forc_norm_hyperplane_vert = sqrt(sum(forc_l(:,:)**2))
       forcmx_hyperplane_vert = 0.d0
       do ia= 1, natm
          fa = dsqrt(forc_l(ia,1)**2 + forc_l(ia,2)**2 + forc_l(ia,3)**2)
          if(fa > forcmx_hyperplane_vert) forcmx_hyperplane_vert = fa
       end do
       if(iprimd >= 1) write(nfout,'(a19,f20.10)') 'forc_norm_hp_v =  ', forc_norm_hyperplane_vert
       if(iprimd >= 1) write(nfout,'(a19,f20.10)') 'forcmx_hp_vert =  ', forcmx_hyperplane_vert
! <--
       if(iprimd >= 1) write(nfout,'(a19,f20.10)') 'modified force =  ', sqrt(sum(forc_l(:,:)**2))
    else
       if(iprimd >= 1) then
          write(nfout,*) 'vector is strange for FI or NEB !!'
          write(nfout,*) 'vector =', vec_norm2
          write(nfout,*) 'opt_mode is changed to SD'
          write(nfout,*) ' total force =', force_t_norm
       end if
    end if

    deallocate(force_t_vert, force_t_para)

  end subroutine modify_forc_hyperplane

! -->   T. Yamasaki 18 July 2008
  logical function m_IS_cnstr_is_fixed_nhp()
    if(constraint_type == FIXED_NORMAL_HYPERVECTOR) then
       m_IS_cnstr_is_fixed_nhp = .true.
    else
       m_IS_cnstr_is_fixed_nhp = .false.
    end if
  end function m_IS_cnstr_is_fixed_nhp

  logical function m_IS_force_check_md_nhp()
    if(forcmx_hyperplane_vert < forccr) then
       m_IS_force_check_md_nhp = .true.
    else
       m_IS_force_check_md_nhp = .false.
    end if
    if(iprimd >= 2) write(nfout,'(" !D forcmx_hyperplane_vert = ",d20.12)') forcmx_hyperplane_vert
  end function m_IS_force_check_md_nhp
! <--

  subroutine initialize_constraint_param
    forcmx_constraint_quench = 1.d+2
    forc_norm_hyperplane_vert = 1.d+2
    forcmx_hyperplane_vert = 1.d+2
  end subroutine initialize_constraint_param

  subroutine evolve_velocities(mdalg,forc_l)
    integer, intent(in) ::                          mdalg
    real(kind=DP), intent(inout), dimension(natm,3) :: forc_l

    integer  :: ia, ifc, j, ip, ic
    real(kind=DP) :: fac(3), frc(3), rm, pdot
    real(kind=DP) :: f_dot_v, vec_norm2
    real(kind=DP), allocatable, dimension(:,:) :: cpd_av, cpd_av_cog
    integer, allocatable, dimension(:) :: ncount, ncount_cog

!xocl spread do/ind_natm
    do ia = 1, natm
       do j=1,3
          if(imdtypxyz(ia,j) == FIX) then
             forc_l(ia,j) = 0.d0
          end if
       enddo
    end do
!xocl end spread

    if(mdalg == STEEPEST_DESCENT) then
       cpd_l(:,:)=forc_l(:,:)/dtio*fi_coefficient                   !! moving force cal
    else
       if(constraint_type==RIGID_BODY_FIX_IN_A_PLANE.or.constraint_type==COG_and_RIGID_BODY_FIX_L &
            & .or.constraint_type==COG_FIX_L) then
          if(num_planes_atoms_are_fixed_cog>=1) then
             allocate(cpd_av_cog(num_planes_atoms_are_fixed_cog,3)); cpd_av_cog = 0.d0
             allocate(ncount_cog(num_planes_atoms_are_fixed_cog))  ; ncount_cog = 0
          end if
          if(num_planes_atoms_are_fixed_rb>=1) then
             allocate(cpd_av(num_planes_atoms_are_fixed_rb,3)); cpd_av = 0.d0
             allocate(ncount(num_planes_atoms_are_fixed_rb))  ; ncount = 0
          end if
          if(nfcatm_cog >=1) then
!xocl spread do/ind_natm
             do ia = 1, nfcatm_cog
                ip = icount_of_ipfixedplane(ia)
                cpd_av_cog(ip,1:3) = cpd_av_cog(ip,1:3) + cpd_l(ia_cnst(ia),1:3)
                ncount_cog(ip) = ncount_cog(ip)+1
             end do
!xocl end spread
             do ip = 1, num_planes_atoms_are_fixed_cog
                cpd_av_cog(ip,1:3) = cpd_av_cog(ip,1:3)/ncount_cog(ip)
                  if(iprimd>=DEBUGPRINTLEVEL) &
                       &  write(nfout,'("!evolve_velocities cpd_av_cog(",i3,") = ",3f11.7, " ncount_cog = ",i5)') &
                       &     ip, cpd_av_cog(ip,1:3),ncount_cog(ip)
             end do
          end if
          if(nfcatm_rb >=1) then
!xocl spread do/ind_natm
             do ia = 1, nfcatm_rb
                ip = icount_of_ipfixedplane(nfcatm_cog+ia)
                cpd_av(ip,1:3) = cpd_av(ip,1:3) + cpd_l(ia_cnst(nfcatm_cog+ia),1:3)
                ncount(ip) = ncount(ip)+1
                if(iprimd>=DEBUGPRINTLEVEL) then
                   write(nfout,'("!evolve_velocities nfcatm_cog+ia = ",i3," cpd_l(",i3,",1:3) = ", 3f12.8)') &
                        & nfcatm_cog+ia, ia_cnst(nfcatm_cog+ia),cpd_l(ia_cnst(nfcatm_cog+ia),1:3)
                end if
             end do
!xocl end spread
             do ip = 1, num_planes_atoms_are_fixed_rb
                cpd_av(ip,1:3) = cpd_av(ip,1:3)/ncount(ip)
                  if(iprimd>=DEBUGPRINTLEVEL) &
                  &   write(nfout,'("!evolve_velocities cpd_av(",i3,")     = ",3f11.7," ncount     = ",i5)') &
                  &   ip, cpd_av(ip,1:3),ncount(ip)
             end do
          end if
       end if

                         if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'("<<evolve_velocities>>")')
!xocl spread do/ind_natm
       do ia = 1, natm
          if(imdtyp(ia) == BONDLENGTH_FIX) cycle

          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
             ifc = icnst_a(ia)
             ic = icount_of_ipfixedplane(ifc)
             if(tmass_rb(ic) > SmallestPositiveNumber*1.d5 ) then
                fac(1:3) = dtio/tmass_rb(ic)
             else
                fac(1:3) = 0.d0
             end if
          else
             do j = 1,3
                if(imdtypxyz(ia,j) == FIX) then
                   fac(j) = 0.d0
                else
                   fac(j) = dtio/amion(ityp(ia))
                end if
             end do
          end if

          if(iprimd >= DEBUGPRINTLEVEL+1) write(nfout,'(i5,3f15.8," fac(1:3)")') ia, fac(1:3)

!!$          do j=1,3
!!$             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
!!$                if(tmass_rb > SmallestPositiveNumber*1.d5 ) then
!!$                   fac(j) = dtio/tmass_rb
!!$                else
!!$                   fac(j) = 0.d0
!!$                end if
!!$             else if(imdtypxyz(ia,j) == FIX) then
!!$                fac(j) = 0.d0
!!$             else
!!$                fac(j) = dtio/amion(ityp(ia))
!!$             end if
!!$          enddo
          if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==FIX_IN_A_PLANE)then
             if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L) then
                ifc = icnst_a(ia)
                ic = icount_of_ipfixedplane(ifc)
                rm = amion(ityp(ia))*rtmass_cog(ic)
                frc = forc_l(ia,1:3) + rm*(fcg_mdfy_cog(ic,1:3) - fcg_cog(ic,1:3))
                !
!!$                frc(3) = frc(3) + fcvect(icnst_a(ia),4)
             else if(imdtyp(ia)==FIX_IN_A_PLANE) then
                ifc = icnst_a(ia)
                pdot = dot_product(fcvect(ifc,1:3),forc_l(ia,1:3))
                frc = forc_l(ia,1:3) - pdot*fcvect(ifc,1:3)
             end if
          else if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
             ifc = icnst_a(ia)
             ic = icount_of_ipfixedplane(ifc)
             frc = fcg_mdfy_rb(ic,1:3)
          else
             frc = forc_l(ia,1:3)
          end if
!!$          if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'("!frc <<evolve_velocities>> ", i5,3f13.6)') ia, frc(1:3)

! -- following 6 lines have been revised by Mamoru Usami, Oct. 2004.--->>
! further modified by T.Yamamoto on Nov. 29, 2005. --->>
          if (mdalg == VERLET .and. iteration_ionic == 1)  then
             !!cpd_l(ia,1:3) = fac*frc(1:3) / 2.d0
             !!cpd_old(ia,1:3) = -cpd_l(ia,1:3)
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                ifc = icnst_a(ia)
                ic = icount_of_ipfixedplane(ifc)
                cpd_l(ia,1:3) = cpd_av(ic,1:3) + 0.5d0*fac(1:3)*frc(1:3)
                if(iprimd >= DEBUGPRINTLEVEL+1) &
                     & write(nfout,'(i5,9f14.8," cpd_l(",i3,",1:3) ,cpd_av(",i3,",1:3),frc <<rigid_body_fix_in_a_plane")') &
                     & ia, cpd_l(ia,1:3), cpd_av(ic,1:3), frc(1:3), ia,ic
             else
                cpd_l(ia,1:3) = cpd_l(ia,1:3) + 0.5d0*fac(1:3)*frc(1:3)
             end if
          else
             if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
                ifc = icnst_a(ia)
                ic = icount_of_ipfixedplane(ifc)
                cpd_l(ia,1:3) = cpd_av(ic,1:3) + fac(1:3)*frc(1:3)
                if(iprimd >= DEBUGPRINTLEVEL+1) write(nfout,'(i5,9f14.8," cpd_l,cpd_av,frc")') &
                     & ia, cpd_l(ia,1:3), cpd_av(ic,1:3),frc(1:3)
             else
                cpd_l(ia,1:3) = cpd_l(ia,1:3) + fac(1:3)*frc(1:3)
                if(iprimd >= DEBUGPRINTLEVEL+1) write(nfout,'(i5,9f14.8," cpd_l,frc,forc_l")') &
                     & ia, cpd_l(ia,1:3), frc(1:3), forc_l(ia,1:3)
             end if
          endif

!  <<-----
       end do
!xocl end spread
!!$       if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE) then
!!$          allocate(cpd_av(3))
!!$          cpd_av = 0.d0
!!$          do ia = 1, nfcatm
!!$             cpd_av = cpd_av + cpd_l(ia_cnst(ia),1:3)
!!$          end do
!!$!!$          cpd_av = cpd_av/nfcatm
!!$          do ia = 1, nfcatm
!!$             cpd_l(ia_cnst(ia),1:3) = cpd_av
!!$          end do
!!$          if(iprimd>=1) then
!!$             write(nfout,'("!evolve_velocities cpd_av = ",3f18.7)') cpd_av(1:3)
!!$          end if
       if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE .or. constraint_type == COG_and_RIGID_BODY_FIX_L &
            & .or. constraint_type==COG_FIX_L) then
          if(allocated(cpd_av)) deallocate(cpd_av)
          if(allocated(ncount)) deallocate(ncount)
          if(allocated(cpd_av_cog)) deallocate(cpd_av_cog)
          if(allocated(ncount_cog)) deallocate(ncount_cog)
       end if
!!$       if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE) then
!!$          deallocate(cpd_av)
!!$       end if
!xocl end spread
    end if

  end subroutine evolve_velocities

  subroutine evolve_cps_constant_move()
    integer :: ia, ip, ic, ifchit, ifc
    real(kind=DP), dimension(3) :: cpsdelta
    integer, allocatable, dimension(:) :: hit_cog, hit_rb ! d(num_planes_atoms_are_fixed_cog), d(num_planes_atoms_are_fixed_rb)

    if(num_planes_atoms_are_fixed_cog>=1) then
       allocate(hit_cog(num_planes_atoms_are_fixed_cog)); hit_cog = 0
    end if
    if(num_planes_atoms_are_fixed_rb>=1) then
       allocate(hit_rb(num_planes_atoms_are_fixed_rb));  hit_rb = 0
    end if

    do ia = 1, natm
       if(imdtyp(ia)==FIX_IN_A_PLANE .or. imdtyp(ia)==COG_FIX .or. imdtyp(ia) == COG_FIX_L &
     .or. imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
          ip = icnst_a(ia)
!!$          cps(ia,1:3) = cps(ia,1:3) + fcvect(ip,1:3)*fcvect(ip,4)
          cps(ia,1:3) = cps(ia,1:3) + fcvect(ip,6:8)*fcvect(ip,4)
          ic = icount_of_ipfixedplane(ip)
          if(iprimd >= DEBUGPRINTLEVEL) then
             write(nfout,'(3i5,7f13.6" : ia, ip, ic, cps(",i4,",1:3), delta-cps, fcvect(",i4,",4")') &
!!$                  ia, ip, ic, cps(ia,1:3), fcvect(ip,1:3)*fcvect(ip,4), ia
                  ia, ip, ic, cps(ia,1:3), fcvect(ip,6:8)*fcvect(ip,4), fcvect(ip,4),ia,ip
             call flush(nfout)
          end if
          if(ip <=nfcatm_cog) then
             if(hit_cog(ic) == 0) then
                distance_cog(ic) = distance_cog(ic) + dabs(fcvect(ip,4))
                hit_cog(ic) = 1
             end if
          end if
          if(nfcatm_cog < ip .and. ip <= nfcatm) then
             if(hit_rb(ic) == 0) then
                distance_rb(ic) = distance_rb(ic) + dabs(fcvect(ip,4))
                hit_rb(ic) = 1
             end if
          end if
       else
          cps(ia,1:3) = cps(ia,1:3) + dtio*cpd_l(ia,1:3)
       end if
    end do

    if(num_planes_atoms_are_fixed_cog>=1) then
       do ic = 1, num_planes_atoms_are_fixed_cog
          if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'(" distance_cog(",i3,") = ",f13.6)') ic, distance_cog(ic)
          if(moved_distance_of_fixed_plane < distance_cog(ic)) moved_distance_of_fixed_plane = distance_cog(ic)
       end do
    end if
    if(num_planes_atoms_are_fixed_rb>=1) then
       do ic = 1, num_planes_atoms_are_fixed_rb
          if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'(" distance_rb (",i3,") = ",f13.6)') ic, distance_rb(ic)
          if(moved_distance_of_fixed_plane < distance_rb(ic)) moved_distance_of_fixed_plane = distance_rb(ic)
       end do
    end if

    if(num_planes_atoms_are_fixed>=1) then
       if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'(" moved_distance_of_fixed_plane = ",f13.6)') moved_distance_of_fixed_plane
       if(allocated(hit_cog)) deallocate(hit_cog)
       if(allocated(hit_rb))  deallocate(hit_rb)
    end if
!!$    moved_distance_of_fixed_plane = moved_distance_of_fixed_plane + dabs(fcvect(1,4))
  end subroutine evolve_cps_constant_move

  subroutine evolve_velocities_constant_force()
    integer  :: ia, ip, ifc
    real(kind=DP) :: frc(3), f, wmass, cpd_old(3), wmass_min
    real(kind=DP), allocatable, dimension(:) ::  delta_moved_distance_cog, delta_moved_distance_rb, f_cog, f_rb


    f = 0.d0
    if(nfcatm>=2 .and. iprimd >= DEBUGPRINTLEVEL) write(nfout,'("<<evolve_velocities_constant_force>>")')
    if(iprimd>=DEBUGPRINTLEVEL) then
       do ia = 1, natm
          if(icnst_a(ia)>0) write(nfout,'(" icnst_a(",i4,") = ",i4)') ia, icnst_a(ia)
       end do
    end if

    if(iprimd>=DEBUGPRINTLEVEL) then
       if(num_planes_atoms_are_fixed_cog>=1) then
          allocate(delta_moved_distance_cog(num_planes_atoms_are_fixed_cog)); delta_moved_distance_cog = 0.d0
       end if
       if(num_planes_atoms_are_fixed_rb>=1) then
          allocate(delta_moved_distance_rb(num_planes_atoms_are_fixed_rb)); delta_moved_distance_rb = 0.d0
       end if
    end if
    if(num_planes_atoms_are_fixed_cog>=1) then
       allocate(f_cog(num_planes_atoms_are_fixed_cog)); f_cog = 0.d0
    end if
    if(num_planes_atoms_are_fixed_rb>=1) then
       allocate(f_rb(num_planes_atoms_are_fixed_cog)); f_rb = 0.d0
    end if

    wmass_min = 1.d10
    do ia = 1, natm
       if(imdtyp(ia) == BONDLENGTH_FIX) cycle
       if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE ) then
!!$          f = force_apply*nfcatm
!!$          frc(1:3) = fcvect(ipcpd(ia),1:3)*f
!!$          cpd_l(ia,1:3) = cpd_l(ia,1:3) + dtio/amion(ityp(ia)) * frc(1:3)
!!$          cpd_old(1:3) = cpd_l(ia,1:3)
          ip = icnst_a(ia)
          ifc = icount_of_ipfixedplane(ip)
          if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L) then
             wmass = tmass_cog(ifc)
             if(apply_constant_force>=1) then
                f_cog(ifc) = f_cog(ifc) + force_apply
             else
                f_cog(ifc) = fcvect(icnst_a(ia),4)
             end if
             if(iprimd>=DEBUGPRINTLEVEL) then
                write(nfout,'(" ia = ",i4," f_cpg(",i3,") = ",f12.8)') ia, ifc, f_cog(ifc)
             end if
          else if(imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
             wmass =  tmass_rb(icount_of_ipfixedplane(ip))
             if(apply_constant_force>=1) then
                f_rb(ifc) = f_rb(ifc) + force_apply
             else
                f_rb(ifc) = fcvect(ip,4)
             end if
             if(iprimd>=DEBUGPRINTLEVEL) then
                write(nfout,'(" ia = ",i4," f_rb(",i3,") = ",f12.8)') ia, ifc, f_rb(ifc)
             end if
          end if
!!$          wmass = tmass(ipfixedplane(ip))
          if(wmass < wmass_min) wmass_min = wmass
!!$          if(wmass >= SmallestPositiveNumber*1.d5) &
!!$               & cpd_l(ia,1:3) = cpd_l(ia,1:3) + dtio/wmass * frc(1:3)
!!$                   if(iprimd >= 1) write(nfout,'(i5,9f13.6," : cpd_old,cpd,frc")') ia, cpd_old(1:3),cpd_l(ia,1:3), frc(1:3)
       end if
    end do
    do ip = 1, nfcatm
       ia = ia_cnst(ip)
       if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
          cpd_old(1:3) = cpd_l(ia,1:3)
          ifc = icount_of_ipfixedplane(ip)
       end if
       if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L) then
          frc(1:3) = fcvect(ip,1:3)*f_cog(ifc)
          cpd_l(ia,1:3) = cpd_l(ia,1:3) + dtio*rtmass_cog(ifc) * frc(1:3)
       else if(imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
          frc(1:3) = fcvect(ip,1:3)*f_rb(ifc)
          cpd_l(ia,1:3) = cpd_l(ia,1:3) + dtio*rtmass_rb(ifc) * frc(1:3)
       end if
         if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'(i5,9f13.6," : cpd_old,cpd,frc")') ia, cpd_old(1:3),cpd_l(ia,1:3), frc(1:3)
    end do


    if(iprimd>=DEBUGPRINTLEVEL) then
       do ip = 1, nfcatm
          ia = ia_cnst(ip)
          ifc = icount_of_ipfixedplane(ip)
          write(nfout,'(" ** ip, ia, ifc = ",3i8)') ip, ia, ifc
          if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L) then
             delta_moved_distance_cog(ifc) = dtio*dtio*rtmass_cog(ifc)*f_cog(ifc)
          else if(imdtyp(ia)==RIGID_BODY_FIX_IN_A_PLANE) then
             delta_moved_distance_rb(ifc) = dtio*dtio*rtmass_rb(ifc)*f_rb(ifc)
          end if
       end do
       moved_distance_of_fixed_plane = 0.d0
       do ip = 1, num_planes_atoms_are_fixed_cog
          if(moved_distance_of_fixed_plane < delta_moved_distance_cog(ip)) &
          moved_distance_of_fixed_plane = delta_moved_distance_cog(ip)
          write(nfout,'(" delta_moved_distance_cog(",i3,") = ",f8.4, " rtmass_cog = ",f12.8, " f_cog = ",f12.8)') &
               & ip, delta_moved_distance_cog(ip),rtmass_cog(ip), f_cog(ip)
       end do
       do ip = 1, num_planes_atoms_are_fixed_rb
          if(moved_distance_of_fixed_plane < delta_moved_distance_rb(ip)) &
          moved_distance_of_fixed_plane = delta_moved_distance_rb(ip)
          write(nfout,'(" delta_moved_distance_rb( ",i3,") = ",f8.4, " rtmass_rb  = ",f12.8, " f_rb  = ",f12.8)') &
               ip, delta_moved_distance_rb(ip), rtmass_rb(ip), f_rb(ip)
       end do
    end if


!!$    if(wmass_min < SmallestPositiveNumber*1.d5) then
!!$       moved_distance_of_fixed_plane = moved_distance_of_fixed_plane
!!$    else
!!$       moved_distance_of_fixed_plane = moved_distance_of_fixed_plane + dtio * dtio/wmass_min * f
!!$    end if
    if(iprimd>=DEBUGPRINTLEVEL) &
    & write(nfout,'(" moved_distance_of_fixed_plane = ",f8.4," max_reach_of_fixed_plane = ",f8.4, " wmass_min = ",d12.4)') &
    & moved_distance_of_fixed_plane, max_reach_of_fixed_plane, wmass_min

    if(iprimd>=DEBUGPRINTLEVEL) then
       if(allocated(delta_moved_distance_cog)) deallocate(delta_moved_distance_cog)
       if(allocated(delta_moved_distance_rb))  deallocate(delta_moved_distance_rb)
    end if

  end subroutine evolve_velocities_constant_force

  subroutine evolve_velocities_cog_or_rb(forc_l)
    real(kind=DP), intent(inout), dimension(natm,3) :: forc_l

    integer  :: ia, ifc, mdalg, ip, j, ic
    real(kind=DP) :: fac(3), frc(3), rm, pdot
    real(kind=DP), allocatable, dimension(:,:) :: cpd_av, cpd_av_cog
    integer, allocatable, dimension(:) :: ncount, ncount_cog

    if(nfcatm_cog <= 0 .and. nfcatm_rb <= 0) return
!xocl spread do/ind_natm
    do ia = 1, natm
       do j=1,3
          if(imdtypxyz(ia,j) == FIX) then
             forc_l(ia,j) = 0.d0
          end if
       enddo
    end do
!xocl end spread

    mdalg = VERLET
    if(constraint_type==RIGID_BODY_FIX_IN_A_PLANE.or.constraint_type==COG_and_RIGID_BODY_FIX_L &
   .or.constraint_type==COG_FIX_L) then
       if(num_planes_atoms_are_fixed_cog>=1) then
          allocate(cpd_av_cog(num_planes_atoms_are_fixed_cog,3)); cpd_av_cog = 0.d0
          allocate(ncount_cog(num_planes_atoms_are_fixed_cog))  ; ncount_cog = 0
       end if
       if(num_planes_atoms_are_fixed_rb>=1) then
          allocate(cpd_av(num_planes_atoms_are_fixed_rb,3)); cpd_av = 0.d0
          allocate(ncount(num_planes_atoms_are_fixed_rb))  ; ncount = 0
       end if
       if(nfcatm_cog >=1) then
!xocl spread do/ind_natm
          do ia = 1, nfcatm_cog
             ip = icount_of_ipfixedplane(ia)
             cpd_av_cog(ip,1:3) = cpd_av_cog(ip,1:3) + cpd_l(ia_cnst(ia),1:3)
             ncount_cog(ip) = ncount_cog(ip)+1
          end do
!xocl end spread
          do ip = 1, num_planes_atoms_are_fixed_cog
             cpd_av_cog(ip,1:3) = cpd_av_cog(ip,1:3)/ncount_cog(ip)
             if(iprimd>=DEBUGPRINTLEVEL) &
                  &  write(nfout,'("!evolve_velocities cpd_av_cog(",i3,") = ",3f11.7, " ncount_cog = ",i5)') &
                  &     ip, cpd_av_cog(ip,1:3),ncount_cog(ip)
          end do
       end if
       if(nfcatm_rb >=1) then
!xocl spread do/ind_natm
          do ia = 1, nfcatm_rb
             ip = icount_of_ipfixedplane(nfcatm_cog+ia)
             cpd_av(ip,1:3) = cpd_av(ip,1:3) + cpd_l(ia_cnst(nfcatm_cog+ia),1:3)
             ncount(ip) = ncount(ip)+1
          end do
!xocl end spread
          do ip = 1, num_planes_atoms_are_fixed_rb
             cpd_av(ip,1:3) = cpd_av(ip,1:3)/ncount(ip)
             if(iprimd>=DEBUGPRINTLEVEL) &
                  &   write(nfout,'("!evolve_velocities cpd_av(",i3,")     = ",3f11.7," ncount     = ",i5)')  &
                  &   ip, cpd_av(ip,1:3),ncount(ip)
          end do
       end if
    end if

                         if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'("<<evolve_velocities>> ia, cpd_l, frc, frc_l")')
!xocl spread do/ind_natm
    do ia = 1, natm
       if(imdtyp(ia) == BONDLENGTH_FIX) cycle

       if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
          ifc = icnst_a(ia)
          ic = icount_of_ipfixedplane(ifc)
          if(tmass_rb(ic) > SmallestPositiveNumber*1.d5 ) then
             fac(1:3) = dtio/tmass_rb(ic)
          else
             fac(1:3) = 0.d0
          end if
       else
          do j = 1,3
             if(imdtypxyz(ia,j) == FIX) then
                fac(j) = 0.d0
             else
                fac(j) = dtio/amion(ityp(ia))
             end if
          end do
       end if


!!$    do ia = 1, natm
!!$       fac = dtio/amion(ityp(ia))
       if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L .or. imdtyp(ia)==FIX_IN_A_PLANE) then
          if(imdtyp(ia)==COG_FIX .or. imdtyp(ia)==COG_FIX_L) then
             ifc = icnst_a(ia)
             ic = icount_of_ipfixedplane(ifc)
             rm = amion(ityp(ia))*rtmass_cog(ic)
             frc = forc_l(ia,1:3) + rm*(fcg_mdfy_cog(ip,1:3) - fcg_cog(ip,1:3))
          else if(imdtyp(ia)==FIX_IN_A_PLANE) then
             ifc = icnst_a(ia)
             pdot = dot_product(fcvect(ifc,1:3),forc_l(ia,1:3))
             frc = forc_l(ia,1:3) - pdot*fcvect(ifc,1:3)
          end if
       else if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
          ifc = icnst_a(ia)
          ic = icount_of_ipfixedplane(ifc)
          frc = fcg_mdfy_rb(ic,1:3)
       else
          frc = forc_l(ia,1:3)
       end if
! -- following 6 lines have been revised by Mamoru Usami, Oct. 2004.--->>
! further modified by T.Yamamoto on Nov. 29, 2005. --->>

       if (mdalg == VERLET .and. iteration_ionic == 1)  then
          !!cpd_l(ia,1:3) = fac*frc(1:3) / 2.d0
          !!cpd_old(ia,1:3) = -cpd_l(ia,1:3)
          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
             ifc = icnst_a(ia)
             ic = icount_of_ipfixedplane(ifc)
             cpd_l(ia,1:3) = cpd_av(ic,1:3) + 0.5d0*fac(1:3)*frc(1:3)
          else
             cpd_l(ia,1:3) = cpd_l(ia,1:3) + 0.5d0*fac(1:3)*frc(1:3)
          end if
       else
          if(imdtyp(ia) == RIGID_BODY_FIX_IN_A_PLANE) then
             ifc = icnst_a(ia)
             ic = icount_of_ipfixedplane(ifc)
             cpd_l(ia,1:3) = cpd_av(ic,1:3) + fac(1:3)*frc(1:3)
          else
             cpd_l(ia,1:3) = cpd_l(ia,1:3) + fac(1:3)*frc(1:3)
          end if
       end if
    end do
    if(constraint_type == RIGID_BODY_FIX_IN_A_PLANE .or. constraint_type == COG_and_RIGID_BODY_FIX_L &
  .or. constraint_type==COG_FIX_L) then
       if(allocated(cpd_av)) deallocate(cpd_av)
       if(allocated(ncount)) deallocate(ncount)
       if(allocated(cpd_av_cog)) deallocate(cpd_av_cog)
       if(allocated(ncount_cog)) deallocate(ncount_cog)
    end if
  end subroutine evolve_velocities_cog_or_rb

  subroutine m_IS_evaluate_v_verlet(mdalg,forc_l_in)
    ! coded by Usami, Oct. 2004
    ! subroutine name is revised by T. Yamasaki
    ! revised by T. Yamasaki, Jun. 2005
    !    'if(mdalg.ne.VERLET) return' is commented out
    integer, intent(in) :: mdalg
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP), allocatable, dimension(:,:)   :: forc_l
!!$    if (mdalg .ne. VERLET) return

    call md1_alloc(mdalg)

    allocate(forc_l(natm,3)); forc_l = forc_l_in

    if(t_ctrl_method == LANGEVIN) then
      if(sw_temperature_profile == ON) then
        call set_new_temperature(nfout,nrsv,temperature_profile)
      endif
      call m_IS_add_friction_force(forc_l)
      call m_IS_add_random_force(forc_l)
    endif

    cpd_old = cpd_l
    call check_constraint(forc_l)
    call evolve_velocities(mdalg,forc_l)
    call get_ekina

    deallocate(forc_l)
    call md1_dealloc
  endsubroutine m_IS_evaluate_v_verlet

  subroutine m_IS_structure_factor_3D(nfout,kgp,ngabc_kngp_l)

    integer, intent(in) :: nfout, kgp
!    integer, intent(in) :: ngabc_kngp_l(kgp,3)
    integer, intent(in) :: ngabc_kngp_l(ista_kngp:iend_kngp,3)

    real(kind=DP), allocatable, dimension(:,:) :: cps_r
    integer             :: i,j
    integer             :: id_sname = -1
#ifdef __TIMER_SUB__
  call timer_sta(1248)
#endif

    call tstatc0_begin('m_IS_structure_factor ',id_sname,1)
#ifdef ENABLE_ESM_PACK
    if(sw_esm==ON) then
       allocate(cps_r(3,natm))
       do i=1,3
          do j=1,natm
             cps_r(i,j) = cps(j,i)
          enddo
       enddo
       call Esm_interface_set_tau(cps_r)
       deallocate(cps_r)
    endif
#endif
    zfm3_l = 0.d0
#ifdef __EDA__
    if(sw_eda==ON) zfm3_l_EDA = 0.d0
#endif
    if(kimg == 1) then
#ifdef NEC_TUNE_SMP
       call structure_factor1(zfm3_l,pos,ngabc_kngp_l)
#else
       call structure_factor1
#endif
    else if(kimg == 2) then
#ifdef NEC_TUNE_SMP
       call structure_factor2(zfm3_l,pos,ngabc_kngp_l)
#else
       call structure_factor2
#endif
    endif

    if(ipristrcfctr >= 2 .and. printable) call wd_zfm3(160)
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
  call timer_end(1248)
#endif
  contains
#ifdef NEC_TUNE_SMP
    subroutine structure_factor1(zfm3_l,pos,ngabc_kngp_l)
#else
    subroutine structure_factor1
#endif
      integer       :: ia, i, it
      real(kind=DP) :: grt
#ifdef NEC_TUNE_SMP
      real(kind=DP) :: zfm3_l(ista_kngp:iend_kngp,ntyp,kimg)
      real(kind=DP) :: pos(natm,3)
!     integer       :: ngabc_kngp_l(kgp,3)
      integer       :: ngabc_kngp_l(ista_kngp:iend_kngp,3)
#endif
#ifdef __TIMER_SUB__
  call timer_sta(1249)
#endif
      do ia = 1, natm
         it = ityp(ia)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
         do i = ista_kngp, iend_kngp  !for mpi
            grt = (pos(ia,1)*ngabc_kngp_l(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
            &    + pos(ia,3)*ngabc_kngp_l(i,3))*PAI2
            zfm3_l(i,it,1) = zfm3_l(i,it,1)+dcos(grt)*iwei(ia)
#ifdef __EDA__
            if(sw_eda==ON) then
! -----  ascat starts modifying  -----
            zfm3_l_EDA(i,ia,1) = zfm3_l_EDA(i,ia,1)+dcos(grt)*iwei(ia)
! -----  ascat ceases modifying  -----
            endif
#endif
         end do
      end do
#ifdef __TIMER_SUB__
  call timer_end(1249)
#endif
    end subroutine structure_factor1

#ifdef NEC_TUNE_SMP
    subroutine structure_factor2(zfm3_l,pos,ngabc_kngp_l)
#else
    subroutine structure_factor2
#endif
      integer       :: ia, i, it
      real(kind=DP) :: grt
#ifdef NEC_TUNE_SMP
      real(kind=DP) :: zfm3_l(ista_kngp:iend_kngp,ntyp,kimg)
      real(kind=DP) :: pos(natm,3)
!     integer       :: ngabc_kngp_l(kgp,3)
      integer       :: ngabc_kngp_l(ista_kngp:iend_kngp,3)
#endif
#ifdef __TIMER_SUB__
  call timer_sta(1250)
#endif
      do ia = 1, natm
         it = ityp(ia)
#ifdef NEC_TUNE_SMP
!CDIR INNER
#endif
         do i = ista_kngp, iend_kngp  !for mpi
            grt = (pos(ia,1)*ngabc_kngp_l(i,1) + pos(ia,2)*ngabc_kngp_l(i,2) &
            &    + pos(ia,3)*ngabc_kngp_l(i,3))*PAI2
            zfm3_l(i,it,1)    = zfm3_l(i,it,1)    + dcos(grt)
            zfm3_l(i,it,kimg) = zfm3_l(i,it,kimg) - dsin(grt)
#ifdef __EDA__
            if(sw_eda==ON) then
! -----  ascat starts modifying  -----
            zfm3_l_EDA(i,ia,1)    = zfm3_l_EDA(i,ia,1)    + dcos(grt)
            zfm3_l_EDA(i,ia,kimg) = zfm3_l_EDA(i,ia,kimg) - dsin(grt)
! -----  ascat ceases modifying  -----
            endif
#endif
         end do
      end do
#ifdef __TIMER_SUB__
  call timer_end(1250)
#endif
    end subroutine structure_factor2

    subroutine wd_zfm3(nelment)
      integer, intent(in) :: nelment

      integer i, it, ij, ik, nnelm
      integer, parameter :: Nwk = 8
      real(kind=DP), pointer, dimension(:) :: zfm3_copy

      nnelm = ((nelment-1)/Nwk + 1)*Nwk
      allocate(zfm3_copy(Nwk))

      write(nfout,'(" === structure factor (first",i5," elements) ===")')&
           &    nnelm
      do it = 1, ntyp
         write(nfout,'(" -- #sp = ",i5," --")') it
         do ik = 1, kimg
            if(kimg == 2) then
               if(ik == 1) then
                  write(nfout,*) ' -- real part      --'
               else
                  write(nfout,*) ' -- imaginary part --'
               endif
            endif
            ij = 1
            do i = 1, nnelm
               if(i >= ista_kngp .and. i <= iend_kngp) then
                  zfm3_copy(ij) = zfm3_l(i,it,ik)
                  if(ij == Nwk) then
                     write(nfout,'(8f10.5)') (zfm3_copy(ij),ij=1,Nwk)
                     ij = 1
                  else
                     ij = ij + 1
                  endif
               end if
            end do
         end do
      end do

      deallocate(zfm3_copy)
    end subroutine wd_zfm3
  end subroutine m_IS_structure_factor_3D

  subroutine m_IS_ewald_3D(nfout,kg,gr_l,kgp,ngabc,ival,ngabc_kngp_l,alf_out)
    integer, intent(in)                     :: nfout, kg,kgp
    real(kind=DP), intent(in)               :: gr_l(ista_kngp:iend_kngp)
    integer, intent(in)                     :: ngabc_kngp_l(ista_kngp:iend_kngp,3)
!    integer, intent(in)                     :: ngabc(kgp,3)
    integer, intent(in)                     :: ngabc(kg,3)
    real(kind=DP), intent(in), dimension(:) :: ival

    real(kind=DP), pointer, dimension(:,:) :: rxyz
    real(kind=DP), pointer, dimension(:)   :: rr
    real(kind=DP), pointer, dimension(:,:) :: cps_fp    ! d(natm2)
    integer,       pointer, dimension(:)   :: ityp_full ! d(natm2)
    real(kind=DP), pointer, dimension(:,:) :: zsum      ! d(newldg)
    real(kind=DP), parameter :: rsphere_radius = 12.5d0
    real(kind=DP), parameter :: phi            =  6.0d0
    integer                  :: neibrd,  alen(3),  newldg
    real(kind=DP) :: alf, alf2, alf24,alfi,alfi2,c1,c2,c3,c4
    real(kind=DP), pointer, dimension(:)   :: ttr
    real(kind=DP),intent(out),optional :: alf_out
#ifdef __EDA__
! -----  ascat starts modifying  -----
    real(kind=DP),allocatable,dimension(:,:,:) :: zsum_per_atom
! -----  ascat ceases modifying  -----
#endif
    integer       :: id_sname = -1
#ifdef __TIMER_SUB__
  call timer_sta(1251)
#endif

    call tstatc0_begin('m_IS_ewald ',id_sname)

    call decide_rxyz_size(rsphere_radius,alen,neibrd)
!!$    if(printable) write(nfout,*) ' ! neibrd = ', neibrd
           allocate(rxyz(neibrd,3))
           allocate(rr(neibrd))
    call substitute_rxyz(alen, neibrd, rxyz, rr)
                   deallocate(rr)
    call decide_alf    ! -> alf
    call decide_newldg ! -> newldg = #Gvectors for summation

           allocate(cps_fp(natm2,3))
           allocate(ityp_full(natm2))
       cps_fp(1:natm,1:3) = pos(1:natm,1:3)
       ityp_full(1:natm)    = ityp(1:natm)
    call rplcps(cps_fp,ityp_full,1,natm2,natm,iwei)
    call cpspac        ! -> cps_fp
    call set_ewald_parameters ! alf2,alf24,alfi,alfi2,c1,c3,c2,c4
#ifdef __EDA__
    if(sw_eda==ON) then
! -----  ascat starts modifying  -----
    allocate(eewald_per_atom_Rspace(natm))
! -----  ascat ceases modifying  -----
    endif
#endif
    call ewald_Rspace_summation
    if(istress==1) call ewald_Rspace_summation_4_stress
                   deallocate(ityp_full)
                   deallocate(cps_fp)
                   deallocate(rxyz)

           allocate(zsum(newldg,kimg))
           allocate(ttr(6))
#ifdef __EDA__
    if(sw_eda==ON) then
! -----  ascat starts modifying  -----
    allocate(zsum_per_atom(newldg,kimg,natm))
    allocate(eewald_per_atom_Gspace(natm))
! -----  ascat ceases modifying  -----
    endif
#endif
#ifdef ENABLE_ESM_PACK
    if(sw_esm==OFF) then
       call ewald_Gspace_summation
       if(istress==1) call ewald_stress_Gspace_summation
    else
        call mpi_allreduce(MPI_IN_PLACE,fxyzew_l,natm*3,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
    endif
#else
    call ewald_Gspace_summation
    if(istress==1) call ewald_stress_Gspace_summation
#endif
                   deallocate(ttr)
                   deallocate(zsum)
    eewald = eewald*0.5d0
#ifdef __EDA__
    if(sw_eda==ON) then
! -----  ascat starts modifying  -----
    eewald_per_atom = eewald_per_atom_Rspace + eewald_per_atom_Gspace
    eewald_per_atom = 0.5d0*eewald_per_atom

    deallocate(zsum_per_atom)
    deallocate(eewald_per_atom_Rspace)
    deallocate(eewald_per_atom_Gspace)
! -----  ascat ceases modifying  -----
    endif
#endif
#ifdef ENABLE_ESM_PACK
    if(printable.and.sw_esm==OFF) call wd_eewald_and_fxyzew
#else
    if(printable) call wd_eewald_and_fxyzew
#endif
    if(present(alf_out)) alf_out = alf
    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
  call timer_end(1251)
#endif
  contains
    subroutine wd_eewald_and_fxyzew
      integer mu

!!$      write(nfout,*) '--- ival --'
!!$      write(nfout,*) ival(1)
      if(printable .and. ipri>=1) write(nfout,'("  Ewald sum = ",d25.12)') eewald
      if(ipri >= 2) then
         do mu = 1, natm
!xocl spread do/ind_natm
!xocl index mu
            write(nfout,710) mu &
            &         ,fxyzew_l(mu,1),fxyzew_l(mu,2),fxyzew_l(mu,2)
!xocl end spread
         end do
710      format(' ',i4,3f25.20)
      end if
    end subroutine wd_eewald_and_fxyzew

    subroutine ewald_Rspace_summation
      real(kind=DP)  :: etr, z, rnms,rnm,rnmc,fc1,fc2,fc3, rnm1, rnm2, rnm3,e0
      integer        :: nu, in, ia
      real(kind=DP)  :: derfc, eewald_rspace, sum_abc1, sum_abc2, sum_abc3
#ifdef __EDA__
! -----  ascat starts modifying  -----
      real(kind=DP)  :: etr2
! -----  ascat ceases modifying  -----
#endif

#ifdef __TIMER_SUB__
  call timer_sta(1258)
#endif

      e0 = -c2 -c3*c4
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      eewald_per_atom_Rspace = 0.d0
      do nu = 1, natm2
        eewald_per_atom_Rspace(nu) = eewald_per_atom_Rspace(nu) - c3*(ival(ityp(nu)))**2
        do ia = 1, natm2
          eewald_per_atom_Rspace(nu) = eewald_per_atom_Rspace(nu) - PAI*alf2/univol*ival(ityp_full(ia))*ival(ityp(nu))
        enddo
      enddo
      do nu = 1, natm
         etr = 0.d0
         do in = 1, neibrd
            do ia = 1, natm
               if(in == 1 .and. ia == nu) cycle
               z = ival(ityp_full(ia))
               rnm1 = cps_fp(nu,1) - cps_fp(ia,1) - rxyz(in,1)
               rnm2 = cps_fp(nu,2) - cps_fp(ia,2) - rxyz(in,2)
               rnm3 = cps_fp(nu,3) - cps_fp(ia,3) - rxyz(in,3)
               rnms = rnm1*rnm1 + rnm2*rnm2 + rnm3*rnm3
               rnm  = dsqrt(rnms)
               fc1  = derfc(alfi*rnm)
               etr2 = fc1/rnm
               eewald_per_atom_Rspace(ia) = eewald_per_atom_Rspace(ia) &
                                 & + 0.5d0*iwei(ia)*ival(ityp_full(ia))*ival(ityp(nu))*etr2
               eewald_per_atom_Rspace(nu) = eewald_per_atom_Rspace(nu) &
                                 & + 0.5d0*iwei(nu)*ival(ityp_full(ia))*ival(ityp(nu))*etr2
            end do
         end do
      enddo
      endif
! -----  ascat ceases modifying  -----
#endif
      if(ipri >= 2 .and. printable) write(6,'(" !! c2+c3*c4 = ",d12.4," <<ewald_Rspace_summation>>")') c2+c3*c4
      z = ival(ityp(1))
      eewald_rspace = 0.d0
      fxyzew_l = 0.d0
#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(etr, z, rxyznm, rnms, rnmc, fc1, fc2, etr, fc3)
#endif
      do nu = ista_atm, iend_atm ! MPI
         etr = 0.d0
         sum_abc1 = 0.d0
         sum_abc2 = 0.d0
         sum_abc3 = 0.d0
         do in = 1, neibrd
            do ia = 1, natm2
               if(in == 1 .and. ia == nu) cycle
               z = ival(ityp_full(ia))
               rnm1 = cps_fp(nu,1) - cps_fp(ia,1) - rxyz(in,1)
               rnm2 = cps_fp(nu,2) - cps_fp(ia,2) - rxyz(in,2)
               rnm3 = cps_fp(nu,3) - cps_fp(ia,3) - rxyz(in,3)
               rnms = rnm1*rnm1 + rnm2*rnm2 + rnm3*rnm3
               rnm  = dsqrt(rnms)
               rnmc = rnm*rnms
               fc1  = derfc(alfi*rnm)
               fc2  = dexp(-alfi2*rnms)*c3
               etr  = etr + z * fc1/rnm
               fc3  = z * (fc1/rnmc+fc2/rnms)
               sum_abc1 = sum_abc1 + rnm1*fc3
               sum_abc2 = sum_abc2 + rnm2*fc3
               sum_abc3 = sum_abc3 + rnm3*fc3
            end do
         end do
#ifdef NEC_TUNE_SMP
!CDIR ATOMIC
#endif
!!$       eewald = eewald + iwei(nu)*ival(ityp(nu))*etr
         eewald_rspace = eewald_rspace + iwei(nu)*ival(ityp(nu))*etr
         fxyzew_l(nu,1) = sum_abc1*ival(ityp(nu))
         fxyzew_l(nu,2) = sum_abc2*ival(ityp(nu))
         fxyzew_l(nu,3) = sum_abc3*ival(ityp(nu))
      end do

      if(npes > 1) then
         call mpi_allreduce(eewald_rspace, eewald, 1, mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
      else
         eewald = eewald_rspace
      end if
#ifdef ENABLE_ESM_PACK
      if(sw_esm==OFF) eewald = eewald + e0
#else
      eewald = eewald + e0
#endif
      if(ipri >= 2 .and. printable) write(6,'(" !! eewald = ",f12.4," <<ewald_Rspace_summation>>")') eewald

#ifdef __TIMER_SUB__
  call timer_end(1258)
#endif
    end subroutine ewald_Rspace_summation

    subroutine ewald_Rspace_summation_4_stress
      real(kind=DP)  :: zu, za, rnms,rnm,rnmc,fc1,fc2,fc3
      integer        :: nu, in, ia, i, j
      real(kind=DP)  :: derfc
      real(kind=DP)  :: rnm1, rnm2, rnm3, s11,s12,s13,s21,s22,s23,s31,s32,s33,f1,f2,f3
      real(kind=DP), pointer, dimension(:,:) :: s_ewt

      s11=0.d0; s12=0.d0; s13=0.d0; s21=0.d0; s22=0.d0; s23=0.d0; s31=0.d0; s32=0.d0; s33=0.d0
      do nu = ista_atm2, iend_atm2 ! MPI
         zu = ival(ityp_full(nu))
         do in = 1, neibrd
            do ia = 1, natm2
               if(in == 1 .and. ia == nu) cycle
               za = ival(ityp_full(ia))
               rnm1 = cps_fp(nu,1) - cps_fp(ia,1) - rxyz(in,1)
               rnm2 = cps_fp(nu,2) - cps_fp(ia,2) - rxyz(in,2)
               rnm3 = cps_fp(nu,3) - cps_fp(ia,3) - rxyz(in,3)
               rnms = rnm1*rnm1 + rnm2*rnm2 + rnm3*rnm3
               rnm  = dsqrt(rnms)
               rnmc = rnm*rnms
               fc1  = derfc(alfi*rnm)
               fc2  = dexp(-alfi2*rnms)*c3
               fc3  = za * zu * (fc1/rnmc+fc2/rnms) * 0.5d0
               f1 = fc3 * (rnm1*rltv(1,1)+rnm2*rltv(2,1)+rnm3*rltv(3,1))/PAI2
               f2 = fc3 * (rnm1*rltv(1,2)+rnm2*rltv(2,2)+rnm3*rltv(3,2))/PAI2
               f3 = fc3 * (rnm1*rltv(1,3)+rnm2*rltv(2,3)+rnm3*rltv(3,3))/PAI2
               s11 = s11 - rnm1 * f1
               s12 = s12 - rnm1 * f2
               s13 = s13 - rnm1 * f3
               s21 = s21 - rnm2 * f1
               s22 = s22 - rnm2 * f2
               s23 = s23 - rnm2 * f3
               s31 = s31 - rnm3 * f1
               s32 = s32 - rnm3 * f2
               s33 = s33 - rnm3 * f3
            end do
         end do
      end do
      s_ew(1,1)=s11; s_ew(1,2)=s12; s_ew(1,3)=s13
      s_ew(2,1)=s21; s_ew(2,2)=s22; s_ew(2,3)=s23
      s_ew(3,1)=s31; s_ew(3,2)=s32; s_ew(3,3)=s33
      if(npes > 1) then
         allocate(s_ewt(3,3))
         call mpi_allreduce(s_ew,s_ewt,9,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
         s_ew = s_ewt
         deallocate(s_ewt)
      end if
      s_ew = s_ew + c2 * 0.5d0 * rltv / PAI2

    end subroutine ewald_Rspace_summation_4_stress

    subroutine ewald_Gspace_summation
#ifdef __TIMER_SUB__
  call timer_sta(1259)
#endif

      call get_zsum
      call add_exp_G2_zsum
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      call ewald_force_Gspace_summation(ngabc)
#else
      call ewald_force_Gspace_summation
#endif
#ifdef __TIMER_SUB__
  call timer_end(1259)
#endif
    end subroutine ewald_Gspace_summation

    subroutine get_zsum
      integer mu, in
      real(kind=DP) :: z, phs, cv1,cv2,cv3
      real(kind=DP),pointer, dimension(:,:) :: zsum_mpi
#ifdef NEC_TUNE4
      real(kind=DP) :: tmp_abc1,tmp_abc2,tmp_abc3
      real(kind=DP) :: tmp_phs0,tmp_phs
#endif

#ifdef __TIMER_SUB__
  call timer_sta(1260)
#endif

      zsum = 0.d0
#ifdef __EDA__
      if(sw_eda==ON) then
      zsum_per_atom = 0.d0
      endif
#endif
      if(kimg == 1) then
         do mu = ista_atm, iend_atm
            z = ival(ityp(mu))*iwei(mu)
            cv1= pos(mu,1)*PAI2
            cv2= pos(mu,2)*PAI2
            cv3= pos(mu,3)*PAI2
#ifdef NEC_TUNE4
!CDIR NODEP
            do in = 1, newldg
               tmp_abc1 = dfloat(ngabc(in,1))
               tmp_abc2 = dfloat(ngabc(in,2))
               tmp_abc3 = dfloat(ngabc(in,3))
               tmp_phs0= 0.0000000000000000e+000
               tmp_phs = tmp_phs0+ cv1*tmp_abc1
               tmp_phs = tmp_phs + cv2*tmp_abc2
               tmp_phs = tmp_phs + cv3*tmp_abc3
               zsum(in,1) = zsum(in,1) + z*dcos(tmp_phs)
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               zsum_per_atom(in,1,mu) = z*dcos(tmp_phs)
! -----  ascat starts modifying  -----
               endif
#endif
            end do
#else
            do in = 1, newldg
               phs = cv1*ngabc(in,1)+cv2*ngabc(in,2)+cv3*ngabc(in,3)
               zsum(in,1) = zsum(in,1) + z*dcos(phs)
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               zsum_per_atom(in,1,mu) = z*dcos(phs)
! -----  ascat starts modifying  -----
               endif
#endif
            end do
#endif
         end do
      else if(kimg == 2) then
         do mu = ista_atm, iend_atm
            z = ival(ityp(mu))*iwei(mu)
            cv1 = pos(mu,1)*PAI2
            cv2 = pos(mu,2)*PAI2
            cv3 = pos(mu,3)*PAI2
#ifdef NEC_TUNE4
!CDIR NODEP
            do in = 1, newldg
               tmp_abc1 = dfloat(ngabc(in,1))
               tmp_abc2 = dfloat(ngabc(in,2))
               tmp_abc3 = dfloat(ngabc(in,3))
               tmp_phs0= 0.0000000000000000e+000
               tmp_phs = tmp_phs0+ cv1*tmp_abc1
               tmp_phs = tmp_phs + cv2*tmp_abc2
               tmp_phs = tmp_phs + cv3*tmp_abc3
               zsum(in,1) = zsum(in,1) + z*dcos(tmp_phs)
               zsum(in,kimg) = zsum(in,kimg) - z*dsin(tmp_phs)
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               zsum_per_atom(in,1,mu)    =  z*dcos(tmp_phs)
               zsum_per_atom(in,kimg,mu) = -z*dsin(tmp_phs)
! -----  ascat starts modifying  -----
               endif
#endif
            end do
#else
            do in = 1, newldg
               phs = cv1*ngabc(in,1)+cv2*ngabc(in,2)+cv3*ngabc(in,3)
               zsum(in,1)    = zsum(in,1)    + z*dcos(phs)
               zsum(in,kimg) = zsum(in,kimg) - z*dsin(phs)
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               zsum_per_atom(in,1,mu)    =  z*dcos(phs)
               zsum_per_atom(in,kimg,mu) = -z*dsin(phs)
! -----  ascat starts modifying  -----
               endif
#endif
            end do
#endif
         end do
      end if
      if(npes > 1) then
         allocate(zsum_mpi(newldg,kimg))
         call mpi_allreduce(zsum, zsum_mpi, newldg*kimg,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
         zsum = zsum_mpi
         deallocate(zsum_mpi)
#ifdef __EDA__
         if(sw_eda==ON) then
         call mpi_allreduce(MPI_IN_PLACE, zsum_per_atom, newldg*kimg*natm,mpi_double_precision, &
         &                  mpi_sum, mpi_ke_world,ierr)
         endif
#endif
      end if

#ifdef __TIMER_SUB__
  call timer_end(1260)
#endif
    end subroutine get_zsum

    subroutine add_exp_G2_zsum
      real(kind=DP)  :: gsc,fc1,fc2,etr
      integer        :: in
      integer        :: ist, iend  !mpi
      real(kind=DP)  :: etr_mpi
#ifdef __EDA__
! -----  ascat starts modifying  -----
      integer        :: iatm, jatm

      if(sw_eda==ON) eewald_per_atom_Gspace = 0.d0
! -----  ascat ceases modifying  -----
#endif
#ifdef __TIMER_SUB__
  call timer_sta(1261)
#endif

      etr = 0.d0
      ist = ista_kngp
      if(ist == 1) ist = 2
      iend = iend_kngp
      if( iend > newldg ) iend = newldg
      if(kimg == 1) then
         if( ist <= iend ) then
            do in = ist, iend  !for mpi
               gsc = gr_l(in)**2
               fc1 = dexp(-gsc*alf24)/gsc
               fc2 = fc1*zsum(in,1)**2
               etr = etr + fc2
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               do iatm = 1, natm
                 do jatm = 1, natm
                   eewald_per_atom_Gspace(iatm) = eewald_per_atom_Gspace(iatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,1,iatm)*zsum_per_atom(in,1,jatm)
                   eewald_per_atom_Gspace(jatm) = eewald_per_atom_Gspace(jatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,1,iatm)*zsum_per_atom(in,1,jatm)
                 enddo
               enddo
! -----  ascat ceases modifying  -----
               endif
#endif
            end do
         endif
      else if(kimg == 2) then
         if( ist <= iend ) then
            do in = ist, iend  !for mpi
               gsc = gr_l(in)**2
               fc1 = dexp(-gsc*alf24)/gsc
               fc2 = fc1*(zsum(in,1)**2 + zsum(in,kimg)**2)
               etr = etr + fc2
#ifdef __EDA__
               if(sw_eda==ON) then
! -----  ascat starts modifying  -----
               do iatm = 1, natm
                 do jatm = 1, natm
                   eewald_per_atom_Gspace(iatm) = eewald_per_atom_Gspace(iatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,1,iatm)*zsum_per_atom(in,1,jatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,kimg,iatm)*zsum_per_atom(in,kimg,jatm)
                   eewald_per_atom_Gspace(jatm) = eewald_per_atom_Gspace(jatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,1,iatm)*zsum_per_atom(in,1,jatm) &
                                              & + 0.5d0*fc1*zsum_per_atom(in,kimg,iatm)*zsum_per_atom(in,kimg,jatm)
                 enddo
               enddo
! -----  ascat ceases modifying  -----
               endif
#endif
            end do
         end if
      end if
      if(npes > 1) then
         call mpi_allreduce(etr,etr_mpi,1,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
         etr = etr_mpi
      end if
      eewald = eewald + etr*c1
#ifdef __EDA__
      if(sw_eda==ON) then
! -----  ascat starts modifying  -----
      eewald_per_atom_Gspace = c1*eewald_per_atom_Gspace
! -----  ascat ceases modifying  -----
      endif
#endif
!!$      if(printable) write(nfout,'(" !! eewald = ",f8.4," etr, c1 = ",2f8.4 " <<add_exp_G2_zsum>>")') eewald,etr,c1
      if(istress==1) s_ew = s_ew - etr * c1 * 0.5d0 * rltv / PAI2
#ifdef __TIMER_SUB__
  call timer_end(1261)
#endif
    end subroutine add_exp_G2_zsum

#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
    subroutine ewald_force_Gspace_summation(ngabc)
#else
    subroutine ewald_force_Gspace_summation
#endif
      integer       :: nu, in
      real(kind=DP) :: fabc(3),gsc,fc1,fc3, phs, g1,g2,g3,p1,p2,p3
      real(kind=DP), pointer, dimension(:,:) :: fxyzew_mpi ! d(natm,3)
      real(kind=DP) :: sum_abc1,sum_abc2,sum_abc3
#ifdef NEC_TUNE4
      real(kind=DP) :: tmp_phs0,tmp_phs
#endif
#if defined(NEC_TUNE2) || defined(NEC_TUNE_SMP)
      integer, intent(in)                     :: ngabc(kgp,3)
#endif

#ifdef __TIMER_SUB__
  call timer_sta(1262)
#endif

      call getttr(rltv,ttr)

#ifdef NEC_TUNE_SMP
!CDIR PARALLEL DO PRIVATE(fabc,cxyz,sum_abc1,sum_abc2,sum_abc3,tmp_abc1,tmp_abc2,tmp_abc3)
!CDIR&PRIVATE (tmp_phs0,tmp_phs,gsc,fc1,fc3,gabc)
#endif
!xocl spread do/ind_natm
      do nu = ista_atm, iend_atm  ! MPI
         fabc = 0.d0
         p1 = pos(nu,1)*PAI2
         p2 = pos(nu,2)*PAI2
         p3 = pos(nu,3)*PAI2
         sum_abc1 = fabc(1)
         sum_abc2 = fabc(2)
         sum_abc3 = fabc(3)
#ifdef NEC_TUNE4
!CDIR NODEP
         do in = 1, newldg - 1
            g1 = dfloat(ngabc(1+in,1))
            g2 = dfloat(ngabc(1+in,2))
            g3 = dfloat(ngabc(1+in,3))
            tmp_phs0= 0.0000000000000000e+000
            tmp_phs = tmp_phs0+ p1*g1
            tmp_phs = tmp_phs + p2*g2
            tmp_phs = tmp_phs + p3*g3
            gsc =    ttr(1)*g1*g1 + ttr(2)*g2*g2 + ttr(3)*g3*g3 + ttr(4)*g1*g2 &
                 & + ttr(5)*g2*g3 + ttr(6)*g3*g1
            fc1 = dexp((-gsc*alf24))/gsc
            fc3 = fc1*dsin(tmp_phs)*zsum(1+in,1)
            if (kimg .eq. 2) then
               fc3 = fc3 + fc1*dcos(tmp_phs)*zsum(1+in,kimg)
            endif
            sum_abc1 = sum_abc1 + g1*fc3
            sum_abc2 = sum_abc2 + g2*fc3
            sum_abc3 = sum_abc3 + g3*fc3
         end do
#else
         do in = 2, newldg
            g1 = ngabc(in,1); g2 = ngabc(in,2); g3 = ngabc(in,3)
            phs = p1*g1 + p2*g2 + p3*g3
            gsc =    ttr(1)*g1*g1 + ttr(2)*g2*g2 + ttr(3)*g3*g3 + ttr(4)*g1*g2 &
                 & + ttr(5)*g2*g3 + ttr(6)*g3*g1
            fc1 = exp(-gsc*alf24)/gsc
            fc3 = fc1*dsin(phs)*zsum(in,1)
            if(kimg .eq. 2) fc3 = fc3 + fc1*dcos(phs)*zsum(in,kimg)
            sum_abc1 = sum_abc1 + g1*fc3
            sum_abc2 = sum_abc2 + g2*fc3
            sum_abc3 = sum_abc3 + g3*fc3
         end do
#endif
         fabc(1) = sum_abc1
         fabc(2) = sum_abc2
         fabc(3) = sum_abc3
         fxyzew_l(nu,1:3) = fxyzew_l(nu,1:3) + matmul(rltv,fabc)*c1*ival(ityp(nu))
      end do
!xocl end spread

      if(npes > 1) then
         allocate(fxyzew_mpi(natm,3))
         call mpi_allreduce(fxyzew_l, fxyzew_mpi, natm*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
         fxyzew_l = fxyzew_mpi
         deallocate(fxyzew_mpi)
      end if

#ifdef __TIMER_SUB__
  call timer_end(1262)
#endif
    end subroutine ewald_force_Gspace_summation

    subroutine ewald_stress_Gspace_summation
      real(kind=DP) :: ga,gb,gc,gsc,fc1,fc2,g1,g2,g3,s11,s12,s13,s21,s22,s23,s31,s32,s33,f1,f2,f3,c0
      integer       :: in,i,j
      integer       :: ist, iend !mpi
      real(kind=DP), pointer, dimension(:,:) :: s_ewt, s_ew_mpi2

      allocate(s_ewt(3,3))
      ist = ista_kngp
      if(ist == 1) ist = 2
      iend = iend_kngp
      if( iend > newldg ) iend = newldg
      s11=0.d0;s12=0.d0;s13=0.d0;s21=0.d0;s22=0.d0;s23=0.d0;s31=0.d0;s32=0.d0;s33=0.d0
      if( ist <= iend ) then
         do in = ist, iend  !for mpi
            ga = ngabc_kngp_l(in,1); gb = ngabc_kngp_l(in,2); gc = ngabc_kngp_l(in,3)
            g1 = rltv(1,1)*ga+rltv(1,2)*gb+rltv(1,3)*gc
            g2 = rltv(2,1)*ga+rltv(2,2)*gb+rltv(2,3)*gc
            g3 = rltv(3,1)*ga+rltv(3,2)*gb+rltv(3,3)*gc
            gsc = gr_l(in)**2
            fc1 = dexp(-gsc*alf24)/gsc
            fc2 = fc1*(zsum(in,1)**2+zsum(in,kimg)**2)/(3.d0-kimg)
            c0 = c1*fc2 /PAI2 * (1.d0/gsc + alf24)
            f1 = c0 * (g1*rltv(1,1)+g2*rltv(2,1)+g3*rltv(3,1))
            f2 = c0 * (g1*rltv(1,2)+g2*rltv(2,2)+g3*rltv(3,2))
            f3 = c0 * (g1*rltv(1,3)+g2*rltv(2,3)+g3*rltv(3,3))

            s11 = s11 + g1*f1
            s12 = s12 + g1*f2
            s13 = s13 + g1*f3
            s21 = s21 + g2*f1
            s22 = s22 + g2*f2
            s23 = s23 + g2*f3
            s31 = s31 + g3*f1
            s32 = s32 + g3*f2
            s33 = s33 + g3*f3
         end do
      end if
      s_ewt(1,1)=s11;s_ewt(1,2)=s12;s_ewt(1,3)=s13
      s_ewt(2,1)=s21;s_ewt(2,2)=s22;s_ewt(2,3)=s23
      s_ewt(3,1)=s31;s_ewt(3,2)=s32;s_ewt(3,3)=s33
      if(npes > 1) then
         allocate(s_ew_mpi2(3,3))
         call mpi_allreduce(s_ewt,s_ew_mpi2,9,mpi_double_precision,mpi_sum,mpi_chg_world,ierr)
         s_ew = s_ew + s_ew_mpi2
         deallocate(s_ew_mpi2)
      else
         s_ew = s_ew + s_ewt
      end if
      deallocate(s_ewt)
    end subroutine ewald_stress_Gspace_summation

    subroutine set_ewald_parameters
      integer       :: ia
      real(kind=DP) :: z
#ifdef __TIMER_SUB__
  call timer_sta(1257)
#endif

      alf2 = alf**2
      alf24= alf2*0.25d0
      alfi = 1.d0/alf
      alfi2= alfi**2
      c1   = 4*PAI/univol
      c3   = 2*alfi/dsqrt(PAI)
      c2   = 0.d0
      c4   = 0.d0
      do ia = 1,natm2
         z  = ival(ityp_full(ia))
         c2 = c2 + z
         c4 = c4 + z**2
      end do
      c2 = PAI*alf2/univol*c2**2
      if(ipri >= 2 .and. printable) then
        write(nfout,'(" natm, natm2 = ",2i5)') natm,natm2
        write(nfout,'(" alf, alf2 = ",2d20.10)') alf, alf2
        write(nfout,160) c1,c2,c3,c4
  160   format(' ',' PI4/UNIVOL = ',F12.6,' Z**2*PI*ALF2/UNIVOL = ',&
     &                         F12.6/,' ',' 2/DSQRT(PI)*ALFI = ',&
     &                        F12.6,' SUM Z(NU)**2 = ',F12.6)
      end if

#ifdef __TIMER_SUB__
  call timer_end(1257)
#endif
    end subroutine set_ewald_parameters

    subroutine cpspac
      real(kind=DP), dimension(3) :: catoms(3)
      integer                     :: i
#ifdef __TIMER_SUB__
  call timer_sta(1256)
#endif

      do i = 1, natm2
         catoms = cps_fp(i,1:3)
         catoms = catoms - nint(catoms)      !Packing
         cps_fp(i,1:3) = matmul(altv,catoms) !Change of coordinate system
      end do
#ifdef __TIMER_SUB__
  call timer_end(1256)
#endif
    end subroutine cpspac

    subroutine decide_newldg
      integer i
      integer  :: iend, newldg_mpi  !mpi

#ifdef __TIMER_SUB__
  call timer_sta(1255)
#endif

      newldg_mpi= 1
      newldg= 1
      if(printable .and. ipri>1) write(nfout,'(" ! kg = ",i9)') kg
      iend = iend_kngp
      if( iend > kg ) iend = kg
      if( ista_kngp <= iend ) then
         do i = ista_kngp, iend  !for mpi
            if(alf*gr_l(i) < phi*2.d0) newldg_mpi = i
         end do
      endif
      if(npes > 1) then
         call mpi_allreduce(newldg_mpi,newldg,1,mpi_integer,mpi_max,mpi_chg_world,ierr)
      else
         newldg = newldg_mpi
      end if
      if(newldg.eq.kg+1) then
         if(printable) write(nfout,'(" **warn alf is too small: alf=",d20.10)') alf
         call phase_error_with_msg(nfout, 'alf is too small',__LINE__,__FILE__)
      endif
      if(printable .and. ipri>1) write(nfout,440) newldg
  440 format(' ',' newldg = ',i8)
#ifdef __TIMER_SUB__
  call timer_end(1255)
#endif
    end subroutine decide_newldg

    subroutine decide_alf
      real(kind=DP)                    :: xalen(3), aamin

#ifdef __TIMER_SUB__
  call timer_sta(1254)
#endif

      xalen = alen*(abs(int(rsphere_radius/alen)) + 1)
      aamin = minval(xalen)
      alf = aamin/phi
      if(printable .and. ipri>1) write(nfout,'("  alf = ",f12.6," aamin = ",f12.6)') alf, aamin
#ifdef __TIMER_SUB__
  call timer_end(1254)
#endif
    end subroutine decide_alf
  end subroutine m_IS_ewald_3D
!===============================================================================

  subroutine m_IS_symm_check_of_pos()

    real(kind=DP), pointer, dimension(:,:) :: cps_full
    integer,       pointer, dimension(:)   :: ityp_full
    real(kind=DP), pointer, dimension(:,:) :: rxyz
    real(kind=DP), pointer, dimension(:)   :: rr
    real(kind=DP), parameter :: rsphere_radius = 12.5d0
!      integer, parameter ::  neibr = 3
!      integer, parameter ::  neibrd =(2*neibr+1)**3
    integer neibrd
    integer,dimension(3) :: alen
    integer              :: id_sname = -1
#ifdef __TIMER_SUB__
  call timer_sta(1228)
#endif

    call tstatc0_begin('m_IS_symm_check_of_pos ',id_sname,1)

    allocate(cps_full(natm2,3))
    allocate(ityp_full(natm2))
    call decide_rxyz_size(rsphere_radius,alen,neibrd)
    allocate(rxyz(neibrd,3))
    allocate(rr(neibrd))
    call substitute_rxyz(alen, neibrd, rxyz, rr)
    deallocate(rr)
    cps_full(1:natm,1:3) = cps(1:natm,1:3)
    ityp_full(1:natm) = ityp(1:natm)
    call rplcps(cps_full, ityp_full, 1, natm2, natm, iwei)
    call symm_check_of_ions_positions_c()

    if( ngen_tl >= 1) call wd_napt_tl()

    deallocate(rxyz)
    deallocate(ityp_full)
    deallocate(cps_full)

    call tstatc0_end(id_sname)
#ifdef __TIMER_SUB__
  call timer_end(1228)
#endif
  contains
    subroutine wd_napt_tl()
      integer :: i,j
      write(nfout,'(" --- napt_tl ---")')
      do i = 1, ngen_tl
         write(nfout,'(i3, " : ",25i4)') i,(napt_tl(j,i),j=1,min(natm,25))
      end do
    end subroutine wd_napt_tl

    subroutine symm_check_of_ions_positions_c()
      integer no, i, j, it, jt
      real(kind=DP) :: f(3),f2(3)
!!$      real(kind=DP), parameter :: ddd = 1.d-12, dde = 1.d-6
      real(kind=DP), parameter :: dde = 1.d-6

      do i = 1, natm
!!$         N_operations :  do no = 1, nopr+af
         N_operations :  do no = 1, nopr+af+ngen_tl
            if(no <= nopr+af) then
               f = matmul(op(1:3,1:3,no),cps_full(i,1:3))+tau(1:3,no,CARTS)
            else
               f = matmul(op_tl(1:3,1:3,no-(nopr+af)),cps_full(i,1:3))+tau_tl(1:3,no-(nopr+af),CARTS)
            end if
            !!$print '(" cps             = ",3f12.8)', cps_full(i,1:3)
            !!$print '(" cps(translated) = ",3f12.8)', f
            AtomSearch: do j = 1,natm2
               f2(1) = abs(cos(sum(rltv(1:3,1)*(f - cps_full(j,1:3))))-1.d0)
               f2(2) = abs(cos(sum(rltv(1:3,2)*(f - cps_full(j,1:3))))-1.d0)
               f2(3) = abs(cos(sum(rltv(1:3,3)*(f - cps_full(j,1:3))))-1.d0)
               !!$print '(" cps             = ",3f12.8)', cps_full(i,1:3)
               !!$print '(" f2              = ",3f12.8)', f2
!!$               if(maxval(f2) <= ddd) then
               if(maxval(f2) <= symmetry_check_criterion) then
                  it=ityp_full(i)
                  jt=ityp_full(j)
                  if(ityp_full(i) /= ityp_full(j) .and. (no <= nopr .or. nopr+af+1<=no) ) then
                     if(printable) write(nfout,9001) i, ityp_full(i), j, ityp_full(j), no
                     call phase_error_with_msg(nfout,'failed symmetry check',__LINE__,__FILE__)
                  endif
                  if( (abs(iatomn(it) - iatomn(jt)) > 1.d-8) .and. (no > nopr .and. no < nopr+af+1)) then
                     if(printable) write(nfout,9002) i,iatomn(it),j,iatomn(jt),no
                     call phase_error_with_msg(nfout,'failed symmetry check',__LINE__,__FILE__)
                  endif
                  if(no<=nopr+af) then
                     napt(i,no) = j
                  else
                     napt_tl(i,no-(nopr+af)) = j
                  end if
                  cycle N_operations
               else if(maxval(f2) <= dde) then
                  if(printable) then
                     write(nfout,'(" -- <<symmetry_check_of_ions_positions_c>> --")')
                     write(nfout,'(" maxval(f2) <= ",d20.8)') dde
                     write(nfout,'(" i = ",i5," no = ",i2," j = ",i5," maxval(f2) = ",d12.5)') i, no, j, maxval(f2)
                     write(nfout,'(" cps(",i5,")        = ",3f20.12)') i,cps_full(i,1:3)
                     write(nfout,'(" op(no)*cps(",i5,") = ",3f20.12)') i,f(1:3)
                     write(nfout,'(" cps(",i5,")        = ",3f20.12)') j,cps_full(j,1:3)
                     write(nfout,'(" f2(1:3)           = ",3f20.12)') f2(1:3)
                  end if
               end if
            end do AtomSearch
            if(printable) then
               write(nfout,*) ' no pair i(atom, no(operation-no.) ', i,no
               write(nfout,'(" cps_full(",i5,")       = ",3f20.12)') i,cps_full(i,1:3)
               write(nfout,'(" op(no)*cps(",i5,")     = ",3f20.12)') i,f(1:3)
            end if
            call phase_error_with_msg(nfout,' no pair of an operated atom <<m_IS_symm_check_of_pos>>',__LINE__,__FILE__)
         enddo N_operations
      enddo
9001  format(i3,'-th site ( atom type = ',i3,' ) is transfered to',/,i3,&
           &     '-th site ( atom type = ',i3,' ) by', i3, '-th operation')
9002  format(i3,'-th site ( atom no.  = ',i3,' ) is transfered to',/,i3,&
           &     '-th site ( atom no.  = ',i3,' ) by', i3, '-th operation')

      if(inversion_symmetry == ON) then
         do i = 1, natm
!!$            N_operations :  do no = 1, nopr+af
            Inv_operations :  do no = 1, 1
               f = -cps_full(i,1:3)
               AtomSearch2: do j = 1,natm2
                  f2(1) = abs(cos(sum(rltv(1:3,1)*(f - cps_full(j,1:3))))-1.d0)
                  f2(2) = abs(cos(sum(rltv(1:3,2)*(f - cps_full(j,1:3))))-1.d0)
                  f2(3) = abs(cos(sum(rltv(1:3,3)*(f - cps_full(j,1:3))))-1.d0)
                  if(maxval(f2) <= symmetry_check_criterion) then
                     it=ityp_full(i)
                     jt=ityp_full(j)
                     if(ityp_full(i) /= ityp_full(j) ) then
                        if(printable) write(nfout,9003) i, ityp_full(i), j, ityp_full(j)
                        call phase_error_with_msg(nfout,'failed symmetry check',__LINE__,__FILE__)
                     endif
                     if( (abs(iatomn(it) - iatomn(jt)) > 1.d-8) ) then
                        if(printable) write(nfout,9004) i,iatomn(it),j,iatomn(jt)
                        call phase_error_with_msg(nfout,'failed symmetry check',__LINE__,__FILE__)
                     endif
                     cycle Inv_operations
                  else if(maxval(f2) <= dde) then
                     if(printable) then
                        write(nfout,'(" -- <<symmetry_check_of_ions_positions_c>> inversion_symmetry--")')
                        write(nfout,'(" maxval(f2) <= ",d20.8)') dde
                        write(nfout,'(" i = ",i5," j = ",i5," maxval(f2) = ",d12.5)') i, j, maxval(f2)
                        write(nfout,'(" cps(",i5,")        = ",3f20.12)') i,cps_full(i,1:3)
                        write(nfout,'(" -cps(",i5,")       = ",3f20.12)') i,f(1:3)
                        write(nfout,'("  cps(",i5,")       = ",3f20.12)') j,cps_full(j,1:3)
                        write(nfout,'(" f2(1:3)            = ",3f20.12)') f2(1:3)
                     end if
                  end if
               end do AtomSearch2
               if(printable) then
                  write(nfout,*) ' no pair i(atom, inversion_symmetry)', i
                  write(nfout,'(" cps_full(",i5,")       = ",3f20.12)') i,cps_full(i,1:3)
                  write(nfout,'(" op(inv)*cps(",i5,")    = ",3f20.12)') i,f(1:3)
               end if
               call phase_error_with_msg(nfout,' "sw_inversion = ON" is invalid <<m_IS_symm_check_of_pos>>', &
               __LINE__,__FILE__)
            enddo Inv_operations
         enddo
9003     format(i3,'-th site ( atom type = ',i3,' ) is transfered to',/,i3,&
              &     '-th site ( atom type = ',i3,' ) by the inversion_symmetry')
9004     format(i3,'-th site ( atom no.  = ',i3,' ) is transfered to',/,i3,&
           &     '-th site ( atom no.  = ',i3,' ) by the inversion_symmetry')
      end if

      if(printable) write(nfout,*) ' -- OK symmetry check of atomic coordinates'
      if(ipri >= 2 .and. printable) then
         do i = 1, natm
            write(nfout,'(" na = ",i5)') i
            write(nfout,'(8i8)') (napt(i,j),j=1,nopr+af)
         enddo
      endif

    end subroutine symm_check_of_ions_positions_c
  end subroutine m_IS_symm_check_of_pos

  subroutine decide_rxyz_size(rsphere_radius, alen, neibrd)
    real(kind=DP), intent(in)  :: rsphere_radius
    integer,       intent(out) :: neibrd

    integer,  intent(out), dimension(3) :: alen
    real(kind=DP)               :: a
    integer i

#ifdef __TIMER_SUB__
  call timer_sta(1252)
#endif

    do i = 1, 3
       a=dabs((rltv(1,i)*altv(1,i)+rltv(2,i)*altv(2,i)+rltv(3,i)*altv(3,i))&
            &/dsqrt(rltv(1,i)*rltv(1,i)+rltv(2,i)*rltv(2,i)+rltv(3,i)*rltv(3,i)))
       alen(i) = abs(int(rsphere_radius/a)) + 1
    enddo

    neibrd = (alen(1)*2+1)*(alen(2)*2+1)*(alen(3)*2+1)
#ifdef __TIMER_SUB__
  call timer_end(1252)
#endif
  end subroutine decide_rxyz_size

  subroutine substitute_rxyz(alen, neibrd, rxyz, rr)
    integer, intent(in), dimension(3) :: alen
    integer, intent(in)               :: neibrd
    real(kind=DP), intent(out), dimension(neibrd,3) :: rxyz
    real(kind=DP), intent(out), dimension(neibrd)   :: rr

    integer i, j, k, mm
    real(kind=DP) :: f(3)

#ifdef __TIMER_SUB__
  call timer_sta(1253)
#endif
#ifdef SX
!CDIR SERIAL
#endif
    mm = 0
#ifdef SX
!CDIR NOVECTOR
#endif
    do i = -alen(1), alen(1)
#ifdef SX
!CDIR NOVECTOR
#endif
       do j = -alen(2), alen(2)
#ifdef SX
!CDIR NOVECTOR
#endif
          do k = -alen(3), alen(3)
             f(1) = i; f(2) = j; f(3) = k
             mm = mm + 1
             rxyz(mm,1:3) = matmul(altv,f)
             rr(mm) = dsqrt(dot_product(rxyz(mm,1:3),rxyz(mm,1:3)))
          enddo
       enddo
    enddo

#ifdef LIBRARY_BUILD
    call hpsort_phase0(neibrd,neibrd,rxyz,rr)
#else
    call hpsort(neibrd,neibrd,rxyz,rr)
#endif
#ifdef SX
!CDIR ENDSERIAL
#endif

#ifdef __TIMER_SUB__
  call timer_end(1253)
#endif
  end subroutine substitute_rxyz

  subroutine m_IS_initialize_mdmode
    mdmode = ORDINA
  end subroutine m_IS_initialize_mdmode

  subroutine m_IS_initialize_cpd_l
    cpd_l = 0.d0
  end subroutine m_IS_initialize_cpd_l

  subroutine m_IS_cps_to_pos
    integer ia
    do ia = 1, natm
       pos(ia,1:3) = matmul(transpose(rltv),cps(ia,1:3))/PAI2
    end do
  end subroutine m_IS_cps_to_pos

  subroutine m_IS_wd_forc(forc_l)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer :: i,j

    if(printable) then
       write(nfout,'(" -- pos, forc_l --")')
       do i = 1, natm
          write(nfout,'(i5,3e13.5,3e12.4)') i,(pos(i,j),j=1,3),(forc_l(i,j),j=1,3)
       end do
    end if
  end subroutine m_IS_wd_forc

  subroutine m_IS_rd_T_parameters(mdalg,nfinp)
    integer, intent(in) :: mdalg,nfinp

    logical             :: EOF_reach, tag_is_found
    real(kind=DP),pointer, dimension(:) :: wk

    if(.not.tag_T_cntrl_is_found) then
       if(mdalg == QUENCHED_CONSTRAINT) call forcp_alloc !-(m_Ionic_System) ->(forcp)
       return
    end if

    allocate(wk(natm*3+1))
!!$    num_Treservoir = nrsv
    call T_control_alloc(nrsv)   ! -(m_Ionic_System) ->(qmass,tkb,cprv,cpqr,forcp)
    call get_qmass_tkb_cprv_cpqr  ! -(contained here) ->(qmass,tkb,cprv,cpqr)
    if(icond == INITIAL) call get_cpd_l   ! -(contained here) ->(cpd_l)
    deallocate(wk)

  contains
    subroutine get_qmass_tkb_cprv_cpqr
      integer :: n

      qmass = 0.d0; tkb = 0.d0

      call rewind_to_tag0(nfinp,len(tag_T_cntrl),tag_T_cntrl &
           &, EOF_reach, tag_is_found, str,len_str)
      if(.not.tag_is_found) &
           call rewind_to_tag0(nfinp,len(tag_T_cntrl2),tag_T_cntrl2&
           &, EOF_reach, tag_is_found, str, len_str)
      do
         read(nfinp,'(a132)',end=1002) str
         call strncmp2(str,len_str,tag_heat_bath,len(tag_heat_bath),tag_is_found)
         if(tag_is_found) then
            call read_LHS_number(str,len_str,wk,4)
            n = nint(wk(1))
            if(icond == INITIAL) then
               call read_RHS_number(str,len_str,wk,4)
            else
               call read_RHS_number(str,len_str,wk,2)
            end if
            if(dabs(qmass(n)) > DELTA .and. printable) &
                 & write(nfout,'(" !D redundant (qmass, tkb) : n = ")') n
            qmass(n) = wk(1);    tkb(n)   = wk(2)
            if(icond == INITIAL) then
               cprv(n)   = wk(3)
               cpqr(n,1) = wk(4);    cpqr(n,2) = 0.d0
            end if
         end if
      end do
1002  continue

      if(printable) then
         do n = 1, nrsv
            write(nfout,'(" heat bath (qmass, tkb ) (",i5,") = (",2d12.4,")")') &
                 & n,qmass(n),tkb(n)
         end do
      end if
    end subroutine get_qmass_tkb_cprv_cpqr

    subroutine get_cpd_l
      integer :: ipnt, i
! initial velocity
      wk = 0.d0
      call rewind_to_tag0(nfinp,len(tag_T_cntrl),tag_T_cntrl,EOF_reach,tag_is_found &
           &, str, len_str)
      if(.not.tag_is_found) &
           call rewind_to_tag0(nfinp,len(tag_T_cntrl2),tag_T_cntrl2 &
           & ,EOF_reach,tag_is_found, str, len_str)

      do while(.true.)
         read(nfinp,'(a132)',end=1003) str
         call strncmp2(str,len_str,tag_atom_velocity,len(tag_atom_velocity),tag_is_found)
         if(tag_is_found) then
            call read_LHS_number(str,len_str,wk(natm*3+1:),1)
            i = nint(wk(natm*3+1))
            if(i <= 0) call phase_error_with_msg(nfout,' ! i < 0 (m_IS_rd_T_parameters)',__LINE__,__FILE__)
            ipnt = 3*(i-1) + 1
            if( dabs(wk(ipnt)) + dabs(wk(ipnt+1)) + dabs(wk(ipnt+2)) > DELTA &
                 & .and. printable) &
                 & write(nfout,'(" !D redundant (velocity) : i = ",i5)') i
            call read_RHS_number(str,len_str,wk(ipnt:),3)
         end if
      end do
1003  continue

!xocl spread do/ind_natm
      do i = 1, natm
         ipnt = 3*(i-1)+1
         cpd_l(i,1) = wk(ipnt)
         cpd_l(i,2) = wk(ipnt+1)
         cpd_l(i,3) = wk(ipnt+2)
      end do
!xocl end spread

      if(printable) then
         write(nfout,'(" !! initial -- velocity --")')
         do i = 1, natm
!xocl spread do/ind_natm
!xocl index i
            write(nfout,'(" velocity atom ",i2," = ", 3d17.9)') i, cpd_l(i,1),cpd_l(i,2),cpd_l(i,3)
!xocl end spread
         enddo
      end if
    end subroutine get_cpd_l
  end subroutine m_IS_rd_T_parameters

  subroutine T_control_dealloc()
    integer :: i
    if(allocated(qmass)) deallocate(qmass)
    if(allocated(tkb)) deallocate(tkb)
    if(allocated(cprv)) deallocate(cprv)
    if(allocated(frsv)) deallocate(frsv)
    if(allocated(cpqr)) deallocate(cpqr)
    if(allocated(forcp)) deallocate(forcp)
    if(allocated(natm_per_thermo)) deallocate(natm_per_thermo)
! ==== KT_add == 2014/06/10
    if (allocated(mdstep_at_start_thermostat) ) deallocate(mdstep_at_start_thermostat)
    if (allocated(mdstep_at_end_thermostat) ) deallocate(mdstep_at_end_thermostat)
    if (allocated(temp_at_start_thermostat) ) deallocate(temp_at_start_thermostat)
    if (allocated(temp_at_end_thermostat) ) deallocate(temp_at_end_thermostat)
! ============== 2014/06/10
    if(nchain>1)then
      deallocate(qmass_c)
      deallocate(cprv_c)
      deallocate(cpqr_c)
    endif
    if(allocated(temperature_profile)) then
       do i=1,nrsv
          deallocate(temperature_profile(i)%tempi)
          deallocate(temperature_profile(i)%tempf)
          deallocate(temperature_profile(i)%qmass)
          deallocate(temperature_profile(i)%tdamp)
          deallocate(temperature_profile(i)%till_n)
       enddo
       deallocate(temperature_profile)
    endif
  end subroutine T_control_dealloc

  subroutine T_control_alloc(n_Treservoir)
    integer, intent(in) :: n_Treservoir
    integer :: ir,iat,icnstrnt_typ,i

    allocate(qmass(n_Treservoir)); qmass = 10.d0
    allocate(tkb(n_Treservoir));   tkb   = 300 * CONST_kB  ! (K)
! ==== KT_add == 2014/06/10
    allocate(mdstep_at_start_thermostat(n_Treservoir))
    allocate(mdstep_at_end_thermostat(n_Treservoir))
    mdstep_at_start_thermostat = 1
    mdstep_at_end_thermostat =   1
    allocate(temp_at_start_thermostat(n_Treservoir))
    allocate(temp_at_end_thermostat(n_Treservoir))
    temp_at_start_thermostat = tkb
    temp_at_start_thermostat = tkb
! ============== 2014/06/10

    allocate(cprv(n_Treservoir));  cprv = 0.d0
    allocate(frsv(n_Treservoir));  frsv = 0.d0
    allocate(cpqr(n_Treservoir,2));  cpqr = 0.d0
    allocate(natm_per_thermo(n_Treservoir)); natm_per_thermo=0
    do iat=1,natm
       ir = icnstrnt_typ(imdtyp(iat),imdalg)
       if (ir>=1) natm_per_thermo(ir) = natm_per_thermo(ir)+1
    enddo
    call forcp_alloc()

    if(nchain>1)then
       allocate(qmass_c(n_Treservoir,nchain));qmass_c = 10.d0
       allocate(cprv_c(n_Treservoir,nchain));cprv_c = 0.d0
       allocate(cpqr_c(n_Treservoir,nchain));cpqr_c = 0.d0
       qmass_c(1:nrsv,1) = qmass(1:nrsv)
       qfactor = 0.d0
       do iat=1,natm
          ir = icnstrnt_typ(imdtyp(iat),imdalg)
          if(ir >= 1) then
              qfactor = qfactor+3.d0*iwei(iat)
          endif
       enddo
       qfactor = 1.d0/dsqrt(qfactor)
       qfactor = 0.1d0*qfactor
       do i=2,nchain
          qmass_c(1:nrsv,i) = qmass(1:nrsv) * qfactor
       enddo
    endif
    !   print *,' -- natm = ',natm
  end subroutine T_control_alloc

  subroutine alloc_temperature_profile(nrsv)
    integer, intent(in) :: nrsv
    integer :: i,iat,ir,icnstrnt_typ
    allocate(temperature_profile(nrsv))
    do i=1,nrsv
       temperature_profile(i)%no = 0
       temperature_profile(i)%natm = 0
       temperature_profile(i)%currprof = 0
       temperature_profile(i)%nprof = 0
    enddo
    do iat=1,natm
       ir = icnstrnt_typ(imdtyp(iat),imdalg)
       if (ir>=1) temperature_profile(ir)%natm = temperature_profile(ir)%natm+1
    enddo
    call T_control_alloc(nrsv)
  end subroutine alloc_temperature_profile

  subroutine forcp_alloc()
    if(allocated(forcp)) deallocate(forcp)
    allocate(forcp(natm,3))
    forcp = 0.d0
    if(ipriinputfile>=2) write(nfout,'(" !! forcp is allocated <<forcp_alloc>>")')
  end subroutine forcp_alloc

  subroutine m_IS_rd_forcp_etc(mdalg,nfcntn)
    integer, intent(in) :: mdalg
    integer, intent(in) :: nfcntn
    integer :: i,j
    logical :: tag_is_found, EOF_reach

    if(iprimd >= 1) write(nfout,'(" -- m_IS_rd_forcp_etc --")')
    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_forcp),tag_forcp, EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(nfout,' tag_forcp is not found <<m_IS_rd_forcp_etc>>',__LINE__,__FILE__)
       else
!!$          read(nfcntn,*)
          if(iprimd >= 1) then
             write(nfout,'(a132)') str
             write(nfout,'(" !** -- forcp -- <<m_IS_rd_forcp_etc>>")')
          end if
          read(nfcntn,*) (forcp(i,1),forcp(i,2),forcp(i,3),i=1,natmorg)
          if(mdalg /= QUENCHED_CONSTRAINT) then
             read(nfcntn,*)
             if(nchain == 1)then
               read(nfcntn,*) (cprv(i),cpqr(i,1),cpqr(i,2),frsv(i),i=1,nrsv)
             else
               read(nfcntn,*) (frsv(i),i=1,nrsv)
               do i=1,nchain
                 read(nfcntn,*) (cprv_c(j,i),cpqr_c(j,i),j=1,nrsv)
               enddo
             endif
          end if
          if(mdalg == BLUEMOON .or. mdalg == QUENCHED_CONSTRAINT) then
             read(nfcntn,*) str
             if(iprimd >= 2) then
                write(nfout,'(a132)') str
                write(nfout,'(" !** -- gca -- <<m_IS_rd_forcp_etc>>")')
             end if
             read(nfcntn,*) (gca(i,1),gca(i,2),gca(i,3),i=1,natmorg)
          end if
       end if
    end if
    if(npes > 1) then
       call mpi_bcast(forcp,natmorg*3 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(nchain==1)then
       call mpi_bcast(cprv,nrsv &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(cpqr,nrsv*2 &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       else
       call mpi_bcast(cprv_c,nrsv*nchain &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       call mpi_bcast(cpqr_c,nrsv*nchain &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       endif
       call mpi_bcast(frsv,nrsv &
            & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       if(mdalg == BLUEMOON .or. mdalg == QUENCHED_CONSTRAINT) then
         call mpi_bcast(gca,natmorg*3 &
              & ,mpi_double_precision,0,MPI_CommGroup,ierr)
       endif
    end if
  end subroutine m_IS_rd_forcp_etc

  subroutine m_IS_wd_forcp_etc(mdalg,nfcntn)
    integer, intent(in) :: mdalg,nfcntn
    integer :: i,j
    if(mype==0) then
       write(nfcntn,'(" -- forcp --")')
       write(nfcntn,'(3d24.16)') (forcp(i,1),forcp(i,2),forcp(i,3),i=1,natm)
       if(mdalg /= QUENCHED_CONSTRAINT) then
          if(nchain == 1)then
            write(nfcntn,'(" -- cprv,cpqr(*,1:2),frsv --")')
            write(nfcntn,'(4d24.16)') &
               & (cprv(i),cpqr(i,1),cpqr(i,2),frsv(i),i=1,nrsv)
          else
            write(nfcntn,'(" -- frsv,cprv_c,cpqr_c --")')
            write(nfcntn,'(4d24.16)') (frsv(i),i=1,nrsv)
            do i=1,nchain
              write(nfcntn,'(4d24.16)') (cprv_c(j,i),cpqr_c(j,i),j=1,nrsv)
            enddo
          endif
       end if
       if(mdalg == BLUEMOON .or. mdalg == QUENCHED_CONSTRAINT) then
          write(nfcntn,'(" -- gca --")')
          write(nfcntn,'(3d24.16)') (gca(i,1),gca(i,2),gca(i,3),i=1,natm)
       end if
    endif
  end subroutine m_IS_wd_forcp_etc

  subroutine m_IS_rd_nrsv(nfcntn)
    integer, intent(in)       :: nfcntn
    logical    :: EOF_reach, tag_is_found

    if(mype==0) then
       call rewind_to_tag0(nfcntn,len(tag_T_cntrl),tag_T_cntrl &
            &, EOF_reach, tag_is_found,str,len_str)
       if(.not.tag_is_found) then
          call rewind_to_tag0(nfcntn,len(tag_T_cntrl2),tag_T_cntrl2&
               &, EOF_reach, tag_is_found, str,len_str)
          if(.not.tag_is_found) then
             call phase_error_with_msg(nfout,' tag_T_cntrl is not found',__LINE__,__FILE__)
          end if
       end if
       read(nfcntn,*)
       read(nfcntn,*) nrsv
    endif
    if(npes > 1) call mpi_bcast(nrsv,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(printable) write(nfout,'(i5, " : nrsv")') nrsv
  end subroutine m_IS_rd_nrsv

  subroutine m_IS_wd_nrsv(nfcntn)
    integer, intent(in)       :: nfcntn
    if(mype==0) then
       write(nfcntn,*) tag_T_cntrl
       write(nfcntn,'(" -- nrsv --")')
       write(nfcntn,'(i10)') nrsv
    endif
  end subroutine m_IS_wd_nrsv

  subroutine m_IS_rd_nrsv_stdin(nfinp)
    integer, intent(in)       :: nfinp

    logical :: eof_reach, tag_is_found
    real(kind=DP) :: wk(1)

    if(imdalg /= T_CONTROL .and. imdalg /= BLUEMOON) return
    call rewind_to_tag0(nfinp,len(tag_T_cntrl), tag_T_cntrl &
         &, EOF_reach, tag_T_cntrl_is_found, str, len_str)
    if(.not.tag_T_cntrl_is_found) &
         & call rewind_to_tag0(nfinp,len(tag_T_cntrl2), tag_T_cntrl2 &
         &, EOF_reach, tag_T_cntrl_is_found, str, len_str)
    if(printable) write(nfout,'(" -- after rewind_to_tag0 --")')

    if(.not.tag_T_cntrl_is_found) then
       if(printable) write(nfout,'(" ! NO TAG of (TEMPERATURE_CONTROL)")')
       return
    end if
    if(printable) write(nfout,'(" -tag_T_cntrl_is_found-")')
!!$ after rewind_to_tag0 --'

    nrsv = 1  ! nrsv: number of heat bath
    tag_is_found = .false.
    do while(.not.tag_is_found)
       read(nfinp,'(a132)') str
       call strncmp2(str,len_str,tag_nrsv,len(tag_nrsv),tag_is_found)
       if(tag_is_found) then
          call read_RHS_number(str,len_str,wk,1)
          nrsv = nint(wk(1))
       end if
    end do
    if(nrsv < 1) nrsv = 1
    call f_strcpy(tag_T_cntrl,len(tag_T_cntrl),str,60)
    if(printable) write(nfout,'(a60)') str(1:60)
    if(printable) write(nfout,'(" nrsv = ",i5)') nrsv
  end subroutine m_IS_rd_nrsv_stdin

  subroutine check_imdtyp(mdalg,external_fix_or_free, ifix,ifree,ifix_P_ifree)
    ! check of the #atoms under constraint and #fixed atoms
    external                external_fix_or_free
    integer, intent(in)  :: mdalg
    integer              :: external_fix_or_free
    integer, intent(out) :: ifix, ifree, ifix_P_ifree
    integer              :: ia,ir

    ifix  = 0
    ifree = 0
    do ia = 1, natm
       ir = external_fix_or_free(imdtyp(ia))
       if(ir > nrsv) then
          if(printable) write(nfout,'(" ia, ir, nrsv = ",3i5)') ia,ir,nrsv
          call phase_error_with_msg(nfout,'invalid thermostat',__LINE__,__FILE__)
       end if
       if(ir == FIX_HBATH)   ifix  = ifix + 1
       if(ir == RELAX_HBATH) ifree = ifree + 1
    end do
    ifix_P_ifree = ifix + ifree
    if(mdalg == BLUEMOON .and. ifree /= 0) then
       if(printable) write(nfout,'(" *** Some atoms are out of control. *** ")')
       call phase_error_with_msg(nfout, ' Some atoms are out of control.',__LINE__,__FILE__)
    end if
    if(ifix == natm .and. nrigid_bodies==0) then
       if(printable) write(nfout,'(" ** All atoms are fixed. ***")')
       if(printable) write(nfout,'(" ! ifix, natm = ",2i5)') ifix, natm
       call phase_error_with_msg(nfout,'All atoms are fixed.',__LINE__,__FILE__)
    else
       if(printable) write(nfout,'(" ! ifix = ",i6)') ifix
       if(printable) write(nfout,'(" ! nrigid_bodies = ",i6)') nrigid_bodies
    end if
  end subroutine check_imdtyp

  subroutine vlcty_accrd2_vVerlet(mdalg,forc_l,ifq,fcg)
    integer              :: icnstrnt_typ
    integer, intent(in) ::                              mdalg
    real(kind=DP), intent(in),dimension(natm,3) ::      forc_l(natm,3)
    integer, intent(out),optional,dimension(natm,3) ::  ifq
    real(kind=DP), optional,      dimension(3) ::       fcg
    integer              :: ia, ir, j
    real(kind=DP)        :: cpdxyz(natm,3),frcxyz(3)
    real(kind=DP)        :: dtml,dtz
    integer, allocatable, dimension(:,:) :: mxyz
    allocate(mxyz(natm,3))
    do ia=1,natm
       mxyz(ia,1:3) = 1
       do j=1,3
          if(imdtypxyz(ia,j) == OFF) mxyz(ia,j) = 0
       enddo
    enddo


    do ia = 1, natm
       ir = icnstrnt_typ(imdtyp(ia),mdalg)         !-(b_Ionic_System)
       if(ir == FIX_HBATH) then
          cpd_l(ia,1:3) = 0.d0
       else
          cpdxyz(ia,1:3) = cpd_l(ia,1:3)
          frcxyz(1:3)    = mxyz(ia,1:3) * (forc_l(ia,1:3) + forcp(ia,1:3))/2.d0
          cpd_l(ia,1:3) = cpdxyz(ia,1:3) + dtio/amion(ityp(ia))*frcxyz(1:3)
       end if
    end do

    if(mdalg == T_CONTROL .or. mdalg == BLUEMOON) then
       do ia = 1, natm
          ir = icnstrnt_typ(imdtyp(ia),mdalg)         !-(b_Ionic_System)
          if(ir >= 1) then
             if(nchain==1)then
             dtz = dtio*cpqr(ir,1)/2.d0
             else
             dtz = dtio*cpqr_c(ir,1)/2.d0
             endif
             cpd_l(ia,1:3) = cpd_l(ia,1:3) - mxyz(ia,1:3) * dtz*cpdxyz(ia,1:3)
          end if
       end do
       if(mdalg == BLUEMOON) then
          do ia = 1, natm
             ir = icnstrnt_typ(imdtyp(ia),mdalg)
             if(ir /= FIX_HBATH) then
                dtml = dtio/amion(ityp(ia))*almda/2
                cpd_l(ia,1:3) = cpd_l(ia,1:3) + mxyz(ia,1:3) * dtml*gca(ia,1:3)
             end if
          end do
       end if
    else if(mdalg == QUENCHED_CONSTRAINT) then
       ifq = 0
       do ia = 1, natm
          ir = icnstrnt_typ(imdtyp(ia),mdalg)         !-(b_Ionic_System)
          if(ir /= FIX_HBATH) then
             dtml = dtio/amion(ityp(ia)) * almda/2
             do j = 1, 3
                cpd_l(ia,j) = cpd_l(ia,j) + mxyz(ia,j) * dtml*gca(ia,j)
                fcg(j)      = forcp(ia,j) + mxyz(ia,j) * almda*gca(ia,j)
                if(cpdxyz(ia,j)*fcg(j) < 0.d0) ifq(ia,j) = 1
             end do
          end if
       end do
    end if
    deallocate(mxyz)
  end subroutine vlcty_accrd2_vVerlet

  subroutine ekina_ekinq_ekbt_and_ega(mdalg,external_irtyp_or_ibath,iw_cnst)
    external external_irtyp_or_ibath
    integer ::                     external_irtyp_or_ibath
    integer,intent(in) ::          mdalg
    integer,intent(in),optional :: iw_cnst
    integer ::                     ia, ir, i, ib
    real(kind=DP) ::               tkin
    nathm = 0      !d(nrsv)
    ekr   = 0.d0   !d(nrsv)
    ekina = 0.d0
    ekbt  = 0.d0
    ekinq = 0.d0

    do ia = 1, natm
       tkin = dot_product(cpd_l(ia,1:3),cpd_l(ia,1:3)) &
            &               * amion(ityp(ia))*iwei(ia)*0.5d0
       ekina = ekina + tkin
       ir = external_irtyp_or_ibath(imdtyp(ia))    !-(b_Ionic_System)
       if(ir >= 1) then
          nathm(ir) = nathm(ir) + 3 * iwei(ia)
          ekr(ir)   = ekr(ir) + tkin
          if(nchain==1)then
             ekbt      = ekbt + 3*iwei(ia)*tkb(ir)*cprv(ir)
          else
             ekbt      = ekbt + 3*iwei(ia)*tkb(ir)*cprv_c(ir,1)
          endif
       end if
    end do
    if(nrigid_bodies>0) ekina = ekina + m_IS_rb_kinetic_energies()
    if(nchain>1)then
       do i=2,nchain
          do ir=1,nrsv
             ekbt = ekbt + tkb(ir) * cprv_c(ir,i)
          enddo
       enddo
    endif

    if(nchain == 1)then
    do ir = 1, nrsv
       ekinq = ekinq + qmass(ir)/2.d0 * cpqr(ir,1)*cpqr(ir,1)
       if(mdalg == T_CONTROL) then
          if(nathm(ir) == 0) cycle
          ekr(ir) = ekr(ir)/nathm(ir)/0.5d0
       end if
    end do
    else
    do i=1,nchain
       do ir=1,nrsv
          ekinq = ekinq + qmass_c(ir,i)/2.d0 * cpqr_c(ir,i) **2
       enddo
    enddo

    if(sw_fix_bond==ON) then
      do ir=1,nrsv
        nathm(ir) = nathm(ir)-nbonds_per_thermo(ir)
      enddo
    endif

    if(mdalg == T_CONTROL) then
      do ir=1,nrsv
          if(nathm(ir) == 0) cycle
          ekr(ir) = ekr(ir)/nathm(ir)/0.5d0
      enddo
    end if
    end if

    if(mdalg == BLUEMOON .or. mdalg == QUENCHED_CONSTRAINT) then
       ir = 1
       ekbt = ekbt - iw_cnst*tkb(ir)*cprv(ir)
       ekr(ir) = ekr(ir)*2/(nathm(ir)-iw_cnst)
    end if

    ega = ekina + ekinq + ekbt

    if(printable) then
       write(nfout,'(" ***** ",i6,"-th time step ; ega= ",1pe16.7 /&
            &," ekina, ekinq, ekbt",3e16.7)') &
            & iteration_ionic, ega, ekina, ekinq, ekbt
       write(nfout,'(" ir, tkb, ekr= ",i3,1p,2e16.7)') (ir,tkb(ir),ekr(ir),ir=1,nrsv)
    end if

    if(ipri >= 2 .and. printable) then
       write(nfout,'(" ! *** ia, cps ***")')
       write(nfout,'(" ",i4,3f20.10)') (ia,cps(ia,1),cps(ia,2),cps(ia,3), ia=1,natm)
       write(nfout,'(" ! *** ir, cprv, cpqr ***")')
       write(nfout,'(" ",i3,2f20.10)') (ir,cprv(ir),cpqr(ir,1), ir=1,nrsv)
    end if
  end subroutine ekina_ekinq_ekbt_and_ega

  subroutine evolve_crdn_ACCRD2_vVerlet(external_irtyp_or_ibath,forc_l)
    external         external_irtyp_or_ibath
    integer       :: external_irtyp_or_ibath
    real(kind=DP),intent(in),dimension(natm,3) :: forc_l
    integer       :: ir, ia, j
    real(kind=DP) :: dtm, dtz
    integer, allocatable, dimension(:,:) :: mxyz
    if(sw_fix_bond==ON) then
      allocate(oldcps(natm,3));oldcps = cps
    endif
    allocate(mxyz(natm,3))
    do ia=1,natm
       mxyz(ia,1:3) = 1
       do j=1,3
          if(imdtypxyz(ia,j) == OFF) mxyz(ia,j) = 0
       enddo
    enddo

    do ia = 1, natm
!---*----*----*----*----*----> Velocity Verlet
       ir = external_irtyp_or_ibath(imdtyp(ia))   !-(b_Ionic_System)
       if(ir == FIX_HBATH) cycle
       dtm = dtio/amion(ityp(ia))/2.d0
       cps(ia,1:3) = cps(ia,1:3)+dtio*(cpd_l(ia,1:3)+dtm*forc_l(ia,1:3)) * mxyz(ia,1:3)
!---*----*----*----*----*----< Velocity Verlet
       if(ir >= 1) then
          if(nchain==1)then
          dtz= dtio*dtio*cpqr(ir,1)/2.d0
          else
          dtz= dtio*dtio*cpqr_c(ir,1)/2.d0
          endif
          cps(ia,1:3)= cps(ia,1:3) - dtz*cpd_l(ia,1:3) * mxyz(ia,1:3)
       end if
!!$       print '(" ia, cps = ",i4,3d16.8)',ia,cps(ia,1),cps(ia,2),cps(ia,3)
    end do
    if(sw_fix_bond==ON) then
      call fixed_bond_coords()
      call update_bond_dsigma_old()
    endif
    deallocate(mxyz)
    if(sw_fix_bond==ON) deallocate(oldcps)
  end subroutine evolve_crdn_ACCRD2_vVerlet

  subroutine evolve_cprv
    integer :: i
      if(nchain == 1)then
       Do i=1, nrsv
       if ( natm_per_thermo(i) == 0 ) cycle
         cprv(i) =  cprv(i) &
              & +dtio*(cpqr(i,1) + 0.5d0*dtio/qmass(i)*frsv(i))
       End Do
       if(printable) write(nfout,'(" cprv = ",5d12.4)') cprv
      else
         cprv_c(1:nrsv,1) = cprv_c(1:nrsv,1)+dtio &
                        & *(cpqr_c(1:nrsv,1)+0.5d0*dtio &
                        & *(frsv(1:nrsv)/qmass_c(1:nrsv,1)-cpqr_c(1:nrsv,1)*cpqr_c(1:nrsv,2)))
         do i=2,nchain-1
         cprv_c(1:nrsv,i) = cprv_c(1:nrsv,i)+dtio &
                        & *(cpqr_c(1:nrsv,i) + 0.5d0*dtio &
                        & *((qmass_c(1:nrsv,i-1)*cpqr_c(1:nrsv,i-1)**2-tkb(1:nrsv))/qmass_c(1:nrsv,i) &
                        & - cpqr_c(1:nrsv,i) * cpqr_c(1:nrsv,i+1)))
         enddo
         cprv_c(1:nrsv,nchain) = cprv_c(1:nrsv,nchain) + dtio &
                        & *(cpqr_c(1:nrsv,nchain) + 0.5d0*dtio &
                        & *(qmass_c(1:nrsv,nchain-1)*cpqr_c(1:nrsv,nchain-1)**2-tkb(1:nrsv))/qmass_c(1:nrsv,nchain))
      endif
  end subroutine evolve_cprv

  subroutine md2_alloc
    allocate(ekr(nrsv))
    allocate(nathm(nrsv))
  end subroutine md2_alloc

  subroutine md2_dealloc
    deallocate(nathm)
    deallocate(ekr)
  end subroutine md2_dealloc

  subroutine heatrsv_chain(mdalg,natm,nrsv,iw_cnst)
    implicit none
    integer, intent(in) ::          mdalg,natm,nrsv
    integer, intent(in),optional :: iw_cnst
    real(kind=DP), allocatable, dimension(:,:) :: cpd_t
    real(kind=DP), allocatable, dimension(:,:) :: cpqr_t,cpqr_prev
    integer :: nmax = 100
    real(kind=DP) :: eps = 1.d-12
    real(kind=DP), allocatable, dimension(:) :: frsv_t
    real(kind=DP)            :: dtz
    integer :: i,j,k,ia,ir,icnstrnt_typ,icount
    logical :: converged

    if(sw_fix_bond==ON) call fixed_bond_velocities(cpd_l)

    allocate(cpd_t(natm,3));cpd_t = cpd_l
    allocate(cpqr_t(nrsv,nchain));cpqr_t = cpqr_c
    allocate(cpqr_prev(nrsv,nchain))
    allocate(frsv_t(nrsv));frsv_t = frsv
    converged = .false.
    do i=1,nmax
       cpqr_prev = cpqr_t
       do ia=1,natm
          ir = icnstrnt_typ(imdtyp(ia),mdalg)
          if(ir>=1)then
             !cpd_t(ia,1:3) = cpd_l(ia,1:3) &
             !            & - 0.5d0*(cpd_t(ia,1:3)*cpqr_t(ir,1)+cpd_l(ia,1:3)*cpqr_c(ir,1)) * dtio
             dtz = 1 + dtio*0.5*cpqr_t(ir,1)
             cpd_t(ia,1:3) = cpd_l(ia,1:3)/dtz
          endif
       enddo
       if(sw_fix_bond==ON) call fixed_bond_velocities(cpd_t)
       call forcrsv(natm,ityp,imdtyp,iwei,amion,cpd_t,tkb,nrsv,mdalg &
            & ,frsv,iw_cnst)         ! -(b_Ionic_System) ->(frsv)
       if(nchain==1) then
       cpqr_t(1:nrsv,1) = cpqr_c(1:nrsv,1) + dtio*((0.5d0/qmass_c(1:nrsv,1))*(frsv(1:nrsv)+frsv_t(1:nrsv)))
       else
       cpqr_t(1:nrsv,1) = cpqr_c(1:nrsv,1) + ((0.5d0/qmass_c(1:nrsv,1)) * (frsv(1:nrsv)+frsv_t(1:nrsv)) &
                      & - (cpqr_t(1:nrsv,1)+cpqr_c(1:nrsv,1))*0.5d0 &
                      & * (cpqr_t(1:nrsv,2)+cpqr_c(1:nrsv,2))*0.5d0) * dtio
       do j=2,nchain-1
       cpqr_t(1:nrsv,j) = cpqr_c(1:nrsv,j) &
                      & + ((1.d0/qmass_c(1:nrsv,j)) &
                      & * (qmass_c(1:nrsv,j-1)*(0.5d0*(cpqr_t(1:nrsv,j-1)+cpqr_c(1:nrsv,j-1)))**2 &
                      & - tkb(1:nrsv)) &
                      & - (cpqr_t(1:nrsv,  j)+cpqr_c(1:nrsv,  j))*0.5d0 &
                      & * (cpqr_t(1:nrsv,j+1)+cpqr_c(1:nrsv,j+1))*0.5d0) * dtio
       enddo
       cpqr_t(1:nrsv,nchain) = cpqr_c(1:nrsv,nchain) &
                      &      + dtio * ((1.d0/qmass_c(1:nrsv,nchain)) &
                      &      * (qmass_c(1:nrsv,nchain-1)*(0.5d0*(cpqr_t(1:nrsv,nchain-1) &
                      &      + cpqr_c(1:nrsv,nchain-1)))**2 - tkb(1:nrsv)))
       endif
       if(printable .and. iprimd>=2)then
       do j=1,nchain
          do k=1,nrsv
          write(nfout,'(a,2i5,3f20.10)') 'chain, thermo, cpqr_t, cpqr_prev', &
            & j,k,cpqr_t(k,j),cpqr_prev(k,j),dabs(cpqr_t(k,j)-cpqr_prev(k,j))
          enddo
       enddo
       endif
       icount = 0
       do j = 1,nchain
          do k=1,nrsv
             if(dabs(cpqr_t(k,j))<eps) then
                icount = icount+1
             else if(dabs(cpqr_prev(k,j)-cpqr_t(k,j))<eps) then
                icount = icount+1
             endif
          enddo
       end do
       if (icount == nchain*nrsv) then
          converged = .true.
          if(printable .and. iprimd>=2) &
       &     write(nfout,'(a,i5,a)') ' !** reservoir iteration converged after ',i,' iterations'
          exit
       endif
    enddo
    if(.not.converged)then
      call phase_error_with_msg(nfout,'!** reservoir iteration did not converge. possible reasons : '//&
      & 'dt/qmass to small or very unstable structure',__LINE__,__FILE__)
    endif
    cpd_l = cpd_t
    cpqr_c = cpqr_t
    deallocate(cpd_t)
    deallocate(cpqr_t)
    deallocate(cpqr_prev)
    deallocate(frsv_t)
  end subroutine heatrsv_chain

  subroutine heatrsv(mdalg,natm,nrsv,iw_cnst)
    ! INPUT  : natm,ityp,imdtyp,iwei,amion,dtio,qmass,tkb,nrsv,imdalg
    ! OUTPUT : cpd_l,cpqr,frsv (,iw_cnst)
    implicit none
    integer, intent(in) ::          mdalg,natm,nrsv
    integer, intent(in),optional :: iw_cnst

    integer :: icnstrnt_typ
    real(kind=DP) :: cpdcor(natm,3),cpqrpre(nrsv),cpqrcor(nrsv),fchk(nrsv)

    integer,       parameter :: mmin = 5, mmax = 100
    real(kind=DP), parameter :: eps = 1.d-12
    integer                  :: ir,ic,ia
    real(kind=DP)            :: dtz

    if(sw_fix_bond==ON) call fixed_bond_velocities(cpd_l)

    cpdcor = cpd_l

    Do ir=1, nrsv
       if ( natm_per_thermo(ir) == 0 ) cycle
       cpqrcor(ir) = cpqr(ir,2) + 2*dtio/qmass(ir)*frsv(ir)
       cpqr(ir,2)  = cpqr(ir,1)
       cpqr(ir,1)  = cpqr(ir,1) + 0.5*dtio/qmass(ir)*frsv(ir)
    End Do

! +++++ predictor-corrector method +++++
    PREDICTOR_CORRECTOR: do ic = 1, mmax
       cpqrpre(1:nrsv)= cpqrcor(1:nrsv)
!xocl spread do/ind_natm
       do ia = 1, natm
          ir = icnstrnt_typ(imdtyp(ia),mdalg)
          if(ir >= 1) then
             dtz = 1 + dtio*0.5*cpqrpre(ir)
             cpdcor(ia,1:3)= cpd_l(ia,1:3)/dtz
          end if
       end do
!xocl end spread

       if(sw_fix_bond==ON) call fixed_bond_velocities(cpdcor)

       call forcrsv(natm,ityp,imdtyp,iwei,amion,cpdcor,tkb,nrsv,mdalg &
            & ,frsv,iw_cnst)         ! -(b_Ionic_System) ->(frsv)

       Do ir=1, nrsv
          if ( natm_per_thermo(ir) == 0 ) cycle
          cpqrcor(ir)= cpqr(ir,1)+frsv(ir)*dtio/qmass(ir)*0.5
       End Do

       if(ic <= mmin) cycle

       FCHECK: do ir = 1, nrsv
          if ( natm_per_thermo(ir) == 0 ) cycle
          fchk(ir) = dabs(cpqrpre(ir)/cpqrcor(ir) - 1.d0)
          if(fchk(ir) >= eps) cycle PREDICTOR_CORRECTOR
       end do FCHECK

       go to 300
    end do PREDICTOR_CORRECTOR

    if(printable) then
       write(nfout,'(" ! Warning in heatrsv ***")')
       write(nfout,'(" cpqrpre, cpqrcor = ",2d20.12)') &
            &                  (cpqrpre(ir),cpqrcor(ir),ir=1,nrsv)
    end if
    call phase_error_with_msg(nfout,'heatrsv unconverged',__LINE__,__FILE__)

300 continue
    if(printable) then
       do ir= 1, nrsv
          write(nfout,'(" ic, cpqrcor, fchk= ",i3,d12.4,d12.4)') ic,cpqrcor(ir),fchk(ir)
       end do
    end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++
    cpd_l = cpdcor
    cpqr(1:nrsv,1) = cpqrcor

  end subroutine heatrsv

  subroutine set_new_temperature(nfout,nrsv,temperature_profile)
    integer, intent(in) :: nfout,nrsv
    type(temperature_profile_t), dimension(nrsv), intent(in) :: temperature_profile
    integer :: i,j,nprof,till_n,ir
    integer :: icount,itemp,nat, iter
    real(kind=DP) :: tempi,tempf,delta_temp,qmasst,tdam
    logical :: temp_was_set
    real(kind=DP) :: qfac
    integer :: nfree
    iter = iteration_ionic
    if(imdalg==P_CONTROL .or. imdalg==PT_CONTROL) iter = iteration_unit_cell
    do i=1,nrsv
       icount = 0
       nprof = temperature_profile(i)%nprof
       temp_was_set = .false.
       do j=1,nprof
          till_n = temperature_profile(i)%till_n(j)
          if(till_n<=2)then
             tkb(i) = temperature_profile(i)%tempf(j)
             temp_was_set = .true.
          else if(iter<=icount+till_n)then
             tempi = temperature_profile(i)%tempi(j)
             tempf = temperature_profile(i)%tempf(j)
             delta_temp = (tempf-tempi)/dble(till_n-1)
             itemp = iter - icount -1
             tkb(i) = tempi + delta_temp * dble(itemp)
             temp_was_set = .true.
          endif
          if(temp_was_set)then
             qmasst = temperature_profile(i)%qmass(j)
             if (qmasst>0) then
                qmass(i) = qmasst
             else
                tdam = temperature_profile(i)%tdamp(j)
                nat = temperature_profile(i)%natm
                nfree = 3*nat
                if(sw_fix_bond==ON) nfree = nfree-nbonds_per_thermo(i)
                qmass(i) = (2.d0*nfree*(tdam/PAI2)**2) * tkb(i)
             endif
             if(printable) write(nfout,'(a,i5,a,i5,a,f15.3,a,f15.3)') &
           &  ' !** thermostat: ',i,' profile: ',j,' temperature: ',tkb(i)/CONST_kB, ' qmass : ',qmass(i)
             exit
          endif
          icount = icount + till_n
       enddo
       if(.not.temp_was_set)then
          tkb(i) = temperature_profile(i)%tempf(nprof)
          qmasst = temperature_profile(i)%qmass(nprof)
          if (qmasst>0) then
             qmass(i) = qmasst
          else
             tdam = temperature_profile(i)%tdamp(nprof)
             nat = temperature_profile(i)%natm
             nfree = 3*nat
             if(sw_fix_bond==ON) nfree = nfree-nbonds_per_thermo(i)
             qmass(i) = (2.d0*nfree*(tdam/PAI2)**2) * tkb(i)
          endif
          if(printable) write(nfout,'(a,i5,a,i5,a,f15.3,a,f15.3)') &
        &  ' !** thermostat: ',i,' profile: ',nprof,' temperature: ',tkb(i)/CONST_kB, ' qmass : ',qmass(i)
       endif
    enddo
    if(nchain>1)then
       qmass_c(1:nrsv,1) = qmass(1:nrsv)
       do i=2,nchain
          qmass_c(1:nrsv,i) = qfactor * qmass(1:nrsv)
       enddo
    endif
  end subroutine set_new_temperature

  subroutine set_new_pressure(nfout,pressure_profile)
    integer, intent(in) :: nfout
    type(pressure_profile_t), intent(in) :: pressure_profile
    integer :: i,j,nprof,till_n,ir
    integer :: icount,ipress,nat
    real(kind=DP) :: pressi,pressf,delta_press,bmasst,tdam
    logical :: press_was_set
    icount = 0
    nprof = pressure_profile%nprof
    press_was_set = .false.
    do j=1,nprof
       till_n = pressure_profile%till_n(j)
       if(till_n<=2)then
          target_pressure = pressure_profile%pressf(j)
          press_was_set = .true.
       else if(iteration_unit_cell<=icount+till_n)then
          pressi = pressure_profile%pressi(j)
          pressf = pressure_profile%pressf(j)
          delta_press = (pressf-pressi)/dble(till_n-1)
          ipress = iteration_unit_cell - icount -1
          target_pressure = pressi + delta_press * dble(ipress)
          press_was_set = .true.
       endif
       if(press_was_set)then
          bmasst = pressure_profile%bmass(j)
          if (bmasst>0) then
             m_baro = bmasst
          else
             tdam = pressure_profile%tdamp_baro(j)
             m_baro = get_default_mbaro(tdam)
          endif
          if(printable) write(nfout,'(a,f20.15,a,f15.3)') ' !** new pressure: ',target_pressure, ' mass_baro : ',m_baro
          exit
       endif
       icount = icount + till_n
    enddo
    if(.not.press_was_set)then
       target_pressure = pressure_profile%pressf(nprof)
       bmasst = pressure_profile%bmass(nprof)
       if (bmasst>0) then
          m_baro = bmasst
       else
          tdam = pressure_profile%tdamp_baro(nprof)
          m_baro = get_default_mbaro(tdam)
       endif
       if(printable) write(nfout,'(a,f20.15,a,f15.3)') ' !** new pressure: ',target_pressure, ' mass_baro : ',m_baro
    endif
  end subroutine set_new_pressure

  subroutine m_IS_md_thermo(forc_l)
    external   irtyp
    integer :: irtyp
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer ::        ifix, ifree
    integer, save  :: ifix_P_ifree = -1
    integer ::        mdalg
    integer ::        id_sname = -1
    integer, parameter :: idummy = -999
    call tstatc0_begin('m_IS_md_thermo ',id_sname)

    if(sw_temperature_profile == ON) then
      call set_new_temperature(nfout,nrsv,temperature_profile)
    endif
    if(imdalg==PT_CONTROL .and. sw_pressure_profile==ON) call set_new_pressure(nfout,pressure_profile)

    mdalg = T_CONTROL
    call md2_alloc                       !-(m_Ionic_System)->(ekr,nathm)

    if(ifix_P_ifree == -1) then
       call check_nrsv                      ! -(contained here)
       call check_imdtyp(mdalg,irtyp,ifix,ifree,ifix_P_ifree) !-(m_Ionic_System)
    end if
    ! --> Velocities at iteration_ionic-th step
    if(iteration_ionic == 1) then
       call forcrsv(natm,ityp,imdtyp,iwei,amion,cpd_l,tkb,nrsv,mdalg,frsv, idummy)
       !  -(b_Ionic_System)(cpd_l,tkb)->(frsv(=force on the thermostat coordinate))
    else
       call vlcty_accrd2_vVerlet(mdalg,forc_l)    !-(m_Ionic_System) ->(cpd_l)
    end if
!!$ 2011.06.06
    if(t_ctrl_method.ne.VELOCITY_SCALING.and.ifix_P_ifree<natm) then
      if(nchain == 1) then
        call heatrsv(mdalg,natm,nrsv)
      else
        call heatrsv_chain(mdalg,natm,nrsv)
      endif
    endif
!!$ 2011.06.06
         ! natm-(ifix_P_ifree):= #atoms in a heat bath
         ! -(m_Ionic_System) (cpd_l,cpqr,frsv,tkb)->(cpd_l,cpqr,frsv)
    forcp = forc_l
    call ekina_ekinq_ekbt_and_ega(mdalg,irtyp)          !-(m_Ionic_System) ->(ekina,..)
    !     <== Kinetic and thermostat energies at iteration_ionic-th step
    call evolve_crdn_ACCRD2_vVerlet(irtyp,forc_l)       !-(m_Ionic_System) ->(cps)
    !     <== Coordinates at (iteration_ionic+1)-th step
    call evolve_cprv                                    !-(m_Ionic_System) ->(cprv)

    if(nrigid_bodies>0) call m_IS_rb_dynamics(forc_l,mdalg)
    if(t_ctrl_method == VELOCITY_SCALING)then
       call scale_velocity()
    endif
    if(nrigid_bodies>0) call scale_velocity_rb()
    if (sw_shift_velocities==ON) call shift_velocities(1)
    call md2_dealloc                                    !-(m_Ionic_System)
    call tstatc0_end(id_sname)
  contains

    subroutine check_nrsv
      if(nrsv > natm) then
         if(printable) then
            write(nfout,'(" *** Too many thermostats. ***")')
            write(nfout,'(" ! nrsv, natm = ")') nrsv,natm
         end if
         call phase_error_with_msg(nfout,'Too many thermostats.',__LINE__,__FILE__)
      end if
    end subroutine check_nrsv

  end subroutine m_IS_md_thermo

  subroutine rattle_v
    real(kind=DP) :: gvu,ggm,dtl,dtlm,almdav
    integer       :: ia, j

    gvu = 0.d0; ggm = 0.d0
    do ia = 1, natm
       gvu = gvu + dot_product(cpd_l(ia,1:3),gca(ia,1:3))
       ggm = ggm + dot_product(gca(ia,1:3),gca(ia,1:3))/amion(ityp(ia))
    end do

    dtl = -gvu/ggm
    do ia = 1, natm
       dtlm  = dtl/amion(ityp(ia))
       cpd_l(ia,1:3)= cpd_l(ia,1:3) + dtlm*gca(ia,1:3)
    end do

    gvu= 0.d0
    do j = 1, 3
       gvu = gvu + dot_product(cpd_l(1:natm,j),gca(1:natm,j))
    end do
    if(ipri >= 2 .and. printable) write(nfout,'(" gvu=",d20.10)') gvu

    almdav= 2.d0* dtl/dtio
    if(ipri >= 2 .and. printable) write(nfout,'(" almdav=",d20.10)') almdav
  end subroutine rattle_v

  subroutine rattle_r(nfcatm)
    integer,       intent(in)          :: nfcatm

    real(kind=DP), allocatable, dimension(:,:) :: gcb
    integer, parameter      :: imax = 50
    real(kind=DP),parameter :: eps  = 1.d-14

    real(kind=DP)  :: dtl,sigma
    integer        :: iter

    allocate(gcb(natm,3))

    dtl = dtio*almda

    iter= 0
!---
    if(cnst_typ == BONDLENGTH_FIX_1 .or. cnst_typ == BONDLENGTH_FIX_2) then
       if(sgmc(1) < DELTA) then
          if(iprimd >= 1) then
             write(nfout,'(" !IS sgmc(1) = ",f12.4)') sgmc(1)
          end if
          call phase_error_with_msg(nfout,  ' sgmc(1) < DELTA <<rattle_r>>',__LINE__,__FILE__)
       end if
    else if(cnst_typ == COG_FIX_L)then
       if(dabs(sgmc(4)) < DELTA) then
          if(iprimd >= 1) then
             write(nfout,'(" !IS sgmc(1:4) = ",4f12.4)') sgmc(1:4)
          end if
          call phase_error_with_msg(nfout,' sgmc(4) < DELTA <<rattle_r>>',__LINE__,__FILE__)
       end if
    end if

    do while(.true.)
       iter = iter+1
       call evolve_cps    !-(contained here), (dtl,gca,dtio,amion)->(cps)
       ! r_{n+1} = r_{n} + (\Delta t)^2/(2 m_i) \lambda g_{i}
       call cnstrnt(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
            &, sgmc,sigma)             !-(b_Ionic_System) ->(sigma)
       if(dabs(sigma) <= eps) exit
       call forc_cnst(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
            & ,sgmc,gcb)              !-(b_Ionic_System) ->(gcb)
       call evolve_almda              !-(contained here)
       !            (ia_cnst,dtio,amion,ityp,gcb,gca,sigma)->(dtl,almda)
       if(iter > imax) then
          if(dabs(sigma) <= eps*10.0) then
             call warning0
             exit
          else
             call stop0     !-(contained here)
          end if
       end if
    end do
    if(iprimd >= 2) then
       write(nfout,'(" !Ionic iter = ",i8," <<rattle_r>>")') iter
       write(nfout,'(" !Ionic   dtl, dtio, almda = ",3d20.8)') dtl,dtio, almda
    end if

    deallocate(gcb)

    if(iprimd > 2 .and. printable) &
         & write(nfout,'(" ! iter, almda, sigma = ",i3,2d12.4)') iter,almda,sigma
  contains
    subroutine warning0
      if(iprimd >= 1) then
         write(nfout,'(" !Ionic *** WARNING *** iter, almda, sigma = ",i5,2d12.4)') iter,almda,sigma
      end if
    end subroutine warning0

    subroutine stop0
      if(printable) then
         write(nfout,'(" *** ERROR in rattle_r ***")')
         write(nfout,'(" ! iter, almda, sigma=",i5,2d12.4)') iter,almda,sigma
      end if
      call phase_error_with_msg(nfout,'<<rattle_r>>',__LINE__,__FILE__)
    end subroutine stop0

    subroutine evolve_almda
      integer       :: ia
      real(kind=DP) :: ggm
      ggm = 0.d0
      do ia = 1, natm
         ggm = ggm + 0.5*dtio/amion(ityp(ia)) &
              &       * dot_product(gcb(ia,1:3),gca(ia,1:3))
      end do

      dtl = -sigma/ggm
      almda = almda + dtl/dtio
      if(iter >= 4 .and. iprimd >= 1) then
         write(nfout,'(" !IS iter = ", i5)') iter
         write(nfout,'(" !IS sigma = ",d20.8)') sigma
         write(nfout,'(" !IS almda = ",d20.8)') almda
         write(nfout,'(" !IS ggm   = ",d20.8)') ggm
         write(nfout,'(" !IS dtio  = ",d20.8)') dtio
         write(nfout,'(" !IS dtl   = ",d20.8)') dtl
         write(nfout,'(" !IS sgmc(1) = ",d20.8)') sgmc(1)
         if(iprimd >= 3) then
            do ia = 1, natm
               write(nfout,'(" !IS ia = ",i5," gcb, gca = ",6f8.4)') ia, gcb(ia,1:3),gca(ia,1:3)
            end do
            do ia = 1, natm
               write(nfout,'(" !IS ia = ",i5," cps = ",3f12.8)') ia, cps(ia,1:3)
            end do
         end if
      end if
    end subroutine evolve_almda

    subroutine evolve_cps
      integer       :: ia
      do ia = 1, natm
         cps(ia,1:3) = cps(ia,1:3)+(dtl*dtio/amion(ityp(ia))/2.d0)*gca(ia,1:3)
      end do
    end subroutine evolve_cps
  end subroutine rattle_r

  subroutine m_IS_md_bluem(forc_l)
  ! === md_bluem ===
  !  -- variables --
  !      tkb  : $k_B T$, temperature
  !     cprv  : $\eta, the thermostat coordinate
  !     cpqr  : $\dot{\eta}, time-derivative of the thermostat coordinate
  !     frsv  : $f_{\eta}$, force on the thermostat coordinate $\eta$
  !      gca  : $\frac{\partial \sigma(\{ \vec{r} \})}{\partial \vec{r_i}},
  !            derivative of constraint
  !    iw_cnst: weight factor for #constraint atom
#ifndef PGI
    external   ibath
    integer :: ibath
#endif
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer :: iw_cnst = -1
    integer :: mdalg
    integer :: id_sname = -1
    call tstatc0_begin('m_IS_md_bluemoon ',id_sname)
    mdalg = BLUEMOON
    call md2_alloc                  !-(m_Ionic_System),->(ekr,nathm)

    nrsv = 1
    if(iw_cnst == -1) call init_md_bluem    !-(contained here) ->(iw_cnst)
    ! --> Velocities at iteration_ionic-th step
    if(iteration_ionic == 1) then
       call forcrsv(natm,ityp,imdtyp,iwei,amion,cpd_l,tkb,nrsv,mdalg&
            &  ,frsv,iw_cnst)       !-(b_Ionic_System) ->(frsv,iw_cnst)
    else
       call vlcty_accrd2_vVerlet(mdalg,forc_l)! -(m_Ionic_System)(forc_l,cpd_l,gca)->(cpd_l)
    end if
    ! --> Rattle
    call forc_cnst(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
         &, sgmc,gca)                  !-(b_Ionic_System) ->(gca)
    if(iteration_ionic /= 1) then
       call rattle_v                   !-(m_Ionic_System)(cpd_l,gca)->(cpd_l)
       call heatrsv(mdalg,natm,nrsv,iw_cnst) !-(m_Ionic_System)
    end if                             !     (tkb)->(cpd_l,cpqr,frsv,iw_cnst)

    forcp = forc_l
    if(printable) call print_frsv_and_cpqr           !-(contained here)
    call ekina_ekinq_ekbt_and_ega(mdalg,ibath,iw_cnst)   !-(m_Ionic_System) ->(ekina,..)
    !     <== Kinetic and thermostat energies at iteration_ionic-th step
    call evolve_crdn_ACCRD2_vVerlet(ibath,forc_l)  !-(m_Ionic_System) ->(cps)
    !     <== Coordinates at (iteration_ionic+1)-th step
    call evolve_cprv            !-(m_Ionic_System)(cprv,cpqr,frsv)->(cprv)
    call rattle_r(nfcatm)           !-(m_Ionic_System) ->(cps,almda)
    if(ipri > 2) print *,' almda = ', almda

    call md2_dealloc
    call tstatc0_end(id_sname)
  contains
    subroutine print_frsv_and_cpqr
      integer :: ir
      do ir = 1, nrsv
         write(nfout,'(" frsv, cpqr= ",2f12.6)') frsv(ir),cpqr(ir,1)
      end do
    end subroutine print_frsv_and_cpqr

    subroutine init_md_bluem
      integer :: ifix, ifree, ifix_P_ifree, m,ia,ir
#ifndef PGI
      external   ibath
      integer :: ibath
#endif

      call check_imdtyp(mdalg,ibath,ifix,ifree,ifix_P_ifree) !-(m_Ionic_System)
      iw_cnst = 1
      do m = 1, nfcatm
         ia = ia_cnst(m)
         ir = ibath(imdtyp(ia))
         if(ir <= RELAX_HBATH) then
            if(printable) then
               write(nfout,'(" *** Inconsistency in imdtyp ***")')
               write(nfout,'(" ! m, ia, ir= ",3i5)') m,ia,ir
            end if
            call phase_error_with_msg(nfout,' ir <= RELAX_HBATH <<m_IS_md_bluem.init_md_bluem>>',__LINE__,__FILE__)
         end if
         if(iwei(ia) == 2) iw_cnst= 2
      end do
      if(printable) write(nfout,'(" ! iw_cnst = ",i5)') iw_cnst
    end subroutine init_md_bluem
  end subroutine m_IS_md_bluem

  subroutine m_IS_md_cnstr(forc_l)
    external   ibath
    integer :: ibath
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l

    integer,       pointer, dimension(:,:)   :: ifq
    real(kind=DP), pointer, dimension(:)     :: fcg
    real(kind=DP), allocatable, dimension(:,:) :: fcvect_work
    integer :: mdalg = QUENCHED_CONSTRAINT
    integer :: id_sname = -1
    call tstatc0_begin('m_IS_md_cnstr ',id_sname)

    if(mdmode == ORDINA) then
    call md2_alloc                  !-(m_Ionic_System)->(ekr,nathm)
    allocate(ifq(natm,3)); allocate(fcg(3))

    call init_md_cnstr                      !-(contained here)

    ! -- Velocities at iteration_ionic-th step
    if(iteration_ionic == 1) then
       almda = 0.d0
    else
       call vlcty_accrd2_vVerlet(mdalg,forc_l,ifq,fcg) !-(m_I.S.)->(cpd_l,ifq)
    end if

    ! --> Rattle
    call forc_cnst(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
         &, sgmc,gca)                  !-(b_Ionic_System) ->(gca)

    call rattle_v                      !-(m_Ionic_System)(cpd_l,gca)->(cpd_l)
    if(iteration_ionic /= 1) then
       call quench_velocity_using_ifq  !-(contained here) (ifq)->(cpd_l)
    end if

    forcp = forc_l
    call evolve_crdn_ACCRD2_vVerlet(ibath,forc_l)       !-(m_Ionic_System) ->(cps)
    !     <== Coordinates at (iteration_ionic+1)-th step
    call rattle_r(nfcatm)              !-(m_Ionic_System) ->(cps,almda)
    call evaluate_forcmx               !-(contained here) ->(forcmx_constraint_quench)
    deallocate(ifq);deallocate(fcg)
    call md2_dealloc

    else if(mdmode == CNSTRA) then
       call move_atoms_normal_to_plane()
       forcmx_constraint_quench = 1.d0

       sgmc = 0.d0
       allocate(fcvect_work(1,8)); fcvect_work = 0.d0
       fcvect_work(1,1:8) = fcvect(1,1:8)
       call m_IS_init_cnstrnt(1,fcvect_work) ! -> sgmc
       deallocate(fcvect_work)
    else
       call phase_error_with_msg(nfout,' Invalid value of mdmode <<m_IS_md>>',__LINE__,__FILE__)
    end if
  contains
    subroutine evaluate_forcmx
      integer       :: ia,ir,icnstrnt_typ
      real(kind=DP) :: fa
      forcmx_constraint_quench = 0.d0
      if(iprimd >= 3) write(nfout,'(" !Ionic(cnstrnt) forcmx_constraint")')
      do ia = 1, natm
         ir = icnstrnt_typ(imdtyp(ia),mdalg)
         if(ir /= FIX_HBATH) then
            fcg(1:3) = forc_l(ia,1:3) + almda*gca(ia,1:3)
            fa = dsqrt(fcg(1)**2 + fcg(2)**2 + fcg(3)**2)
            if(iprimd >= 3) then
               write(nfout,'(" !     fa(",i4,") = ",d20.8)') ia, fa
            end if
            if(fa > forcmx_constraint_quench) forcmx_constraint_quench = fa
         end if
      end do
    end subroutine evaluate_forcmx

    subroutine quench_velocity_using_ifq
      cpd_l = cpd_l * (1-ifq)
    end subroutine quench_velocity_using_ifq

    subroutine init_md_cnstr
      integer :: ifix, ifree, ifix_P_ifree, m,ia,ir
#ifndef PGI
      external   ibath
      integer :: ibath
#endif

      call check_imdtyp(mdalg,ibath,ifix,ifree,ifix_P_ifree) !-(m_Ionic_System)
      do m = 1, nfcatm
         ia = ia_cnst(m)
         ir = ibath(imdtyp(ia))
         if(ir < RELAX_HBATH) then
            if(printable) then
               write(nfout,'(" *** Inconsistency in imdtyp ***")')
               write(nfout,'(" ! m, ia, ir= ",3i5)') m,ia,ir
            end if
            call phase_error_with_msg(nfout,' ir < RELAX_HBATH <<m_IS_md_cnstr.init_md_cnstr>>',__LINE__,__FILE__)
         end if
      end do
    end subroutine init_md_cnstr

  end subroutine m_IS_md_cnstr

  logical function m_IS_moved_distance_of_plane_is_over()
    m_IS_moved_distance_of_plane_is_over = .false.
    if(moved_distance_of_fixed_plane >= max_reach_of_fixed_plane) m_IS_moved_distance_of_plane_is_over = .true.
    if(iprimd >= 1) then
       write(nfout,'(" moved_distance_of_fixed_plane = ",f8.4," max_reach_of_fixed_plane = ",f8.4)')&
            & moved_distance_of_fixed_plane, max_reach_of_fixed_plane
    end if
  end function m_IS_moved_distance_of_plane_is_over

  logical function m_IS_force_check_md_cnstr()
    if(forcmx_constraint_quench < forccr) then
       m_IS_force_check_md_cnstr = .true.
    else
       m_IS_force_check_md_cnstr = .false.
    end if
    if(iprimd >= 2) write(nfout,'(" !D forcmx_constraint_quench = ",d20.12)') forcmx_constraint_quench
  end function m_IS_force_check_md_cnstr

  subroutine m_IS_alloc_cnstrvectors_etc(mdalg)
    integer :: mdalg
    write(nfout,'(" nfcatm = ", i8, "  mdalg = ", i5, " <<m_IS_alloc_cnstrvectors_etc>>")') nfcatm, mdalg
    call flush(nfout)
    if(nfcatm > 0) then
       if(.not.allocated(ia_cnst)) allocate(ia_cnst(nfcatm))
       if(mdalg == GDIIS .or. mdalg == VERLET .or. mdalg == QUENCHED_MD &
            & .or. mdalg == CG_STROPT .or. mdalg==SD_MD .or. mdalg==QUENCHED_CONSTRAINT .or. mdalg==CG_STROPT2 ) then
          if(printable) write(nfout,'(" !!f nfcatm = ",i5," <<m_IS_alloc_cnstrvectors_etc>>")') nfcatm
          ! --- fcvect --
          if(allocated(fcvect)) then
             if(iprimd >= 1) write(nfout,'(" !!f fcvect is already allocated")')
          else
             allocate(fcvect(nfcatm,8)); fcvect = 0.d0
             if(iprimd >= 1) write(nfout,'(" !!f fcvect is allocated and fcvect = 0.d0")')
          end if
          ! --- ipfixedplane --
          if(allocated(ipfixedplane)) then
             if(iprimd >= 1) write(nfout,'(" !!f ipfixedplane is already allocated")')
          else
             allocate(ipfixedplane(nfcatm)); ipfixedplane = 0
             if(iprimd >= 1) write(nfout,'(" !!f ipfixedplane is allocated and ipfixedplane = 0")')
          end if
          ! --- icount_of_ipfixedplane --
          if(allocated(icount_of_ipfixedplane)) then
             if(iprimd >= 1) write(nfout,'(" !!f icount_of_ipfixedplane is already allocated")')
          else
             allocate(icount_of_ipfixedplane(nfcatm)); icount_of_ipfixedplane = 0
             if(iprimd >= 1) write(nfout,'(" !!f icount_of_ipfixedplane is allocated and icount_of_ipfixedplane = 0")')
          end if
          ! -- relax_in_fixedplane --
          if(allocated(relax_in_fixedplane)) then
             if(iprimd >= 1) write(nfout,'(" !!f relax_in_fixedplane is allocated")')
          else
             allocate(relax_in_fixedplane(nfcatm)); relax_in_fixedplane = YES
             if(iprimd >= 1) &
             write(nfout,'(" !!f relax_in_fixedplane is allocated and relax_in_fixedplane is set YES (=",i2,")")') YES
          end if
       end if
       if(mdalg == BLUEMOON .or. mdalg == QUENCHED_CONSTRAINT) then
          if(allocated(gca)) then
             if(iprimd >= 1) write(nfout,'(" !!f gca is already allocated <<m_IS_alloc_cnstrvectors_etc>>")')
          else
             allocate(gca(natm,3));      gca = 0.d0
             if(iprimd >= 1) write(nfout,'(" !!f gca is allocated and gca = 0.d0 <<m_IS_alloc_cnstrvectors_etc>>")')
          end if
       end if
       call forcp_alloc()
    end if
  end subroutine m_IS_alloc_cnstrvectors_etc

!!$  subroutine alloc_distances_of_planes()
!!$    if(num_planes_atoms_are_fixed_cog>=1) then
!!$       allocate(distance_cog(num_planes_atoms_are_fixed_cog)); distance_cog = 0.d0
!!$    end if
!!$    if(num_planes_atoms_are_fixed_rb>=1) then
!!$       allocate(distance_rb(num_planes_atoms_are_fixed_rb)); distance_rb = 0.d0
!!$    end if
!!$  end subroutine alloc_distances_of_planes

  subroutine m_IS_alloc_rigid_body_vectors_etc(mdalg,nfcatm)
    integer,intent(in) :: mdalg
    integer,intent(out) :: nfcatm

    integer :: i, icount
    icount = 0
    do i = 1, natm
       if(imdtyp(i) == RIGID_BODY_FIX_IN_A_PLANE) icount = icount + 1
    end do
    nfcatm = icount

    if(nfcatm > 0) then
       if(.not.allocated(ia_cnst)) allocate(ia_cnst(nfcatm))
       icount = 0
       do i = 1, natm
          if(imdtyp(i) == RIGID_BODY_FIX_IN_A_PLANE) then
             icount = icount + 1
             if(icount <= nfcatm) ia_cnst(icount) = i
          end if
       end do

       if(mdalg == GDIIS .or. mdalg == VERLET .or. mdalg == QUENCHED_MD &
            & .or. mdalg == CG_STROPT .or. mdalg==SD_MD .or. mdalg==QUENCHED_CONSTRAINT .or. mdalg==CG_STROPT2 ) then
          if(printable) write(nfout,'(" !!f nfcatm = ",i5," <<m_IS_alloc_rigid_body_vectors_etc>>")') nfcatm
          ! --- fcvect --
          if(allocated(rigid_body_vect)) then
             if(iprimd >= 1) write(nfout,'(" !!f rigid_body_vect is already allocated")')
          else
             allocate(rigid_body_vect(4)); rigid_body_vect = 0.d0
          end if
          ! --- ipfixedplane --
          if(allocated(ipfixedplane)) then
             if(iprimd >= 1) write(nfout,'(" !!f ipfixedplane is already allocated")')
          else
             allocate(ipfixedplane(nfcatm)); ipfixedplane = 0
             if(iprimd >= 1) write(nfout,'(" !!f ipfixedplane is allocated and ipfixedplane = 1")')
          end if
       end if
       call forcp_alloc()
    end if
  end subroutine m_IS_alloc_rigid_body_vectors_etc

  subroutine substitute_ia_cnst()
    integer :: i, icount
    icount = 0
    do i = 1, natm
       if(imdtyp(i) == BONDLENGTH_FIX .or. imdtyp(i) == BONDLENGTH_FIX_1 &
            & .or. imdtyp(i) == BONDLENGTH_FIX_2 .or. imdtyp(i) == COG_FIX) then
          icount = icount + 1
          if(icount <= nfcatm) ia_cnst(icount) = i
       end if
    end do
  end subroutine substitute_ia_cnst

  subroutine m_IS_cp_works2fcvect_etc(mdalg,n,ia_cnst_work,fcvect_work)
    integer, intent(in) ::                       mdalg,n
    integer, intent(in), dimension(n) ::         ia_cnst_work
    real(kind=DP), intent(in), dimension(n,4) :: fcvect_work
    integer :: i,j

    ia_cnst(1:nfcatm) = ia_cnst_work(1:nfcatm)
    if(mdalg == VERLET .or. mdalg == QUENCHED_MD) then
       fcvect(1:nfcatm,1:4) = fcvect_work(1:nfcatm,1:4)
       if(printable) then
          do i = 1, nfcatm
             write(nfout,'(" fcvect(",i5,",:) = ",4f12.6," <<m_IS_cp_work2fcvect_etc>>")') &
                  & i,(fcvect(i,j),j=1,4)
          end do
       end if
    end if
  end subroutine m_IS_cp_works2fcvect_etc

  subroutine m_IS_init_cnstrnt(n,fcv_wk)
    integer, intent(in)                       :: n
    real(kind=DP), intent(in), dimension(n,4) :: fcv_wk

    integer :: ia
    real(kind=DP) :: sigma

    sgmc = 0.d0;sigma = 0.d0
    if(nfcatm == 0) return

    if(iprimd>=2) then
       do ia = 1, natm
          write(nfout,'("!m_IS_init_cnstrnt:  ia, imdtyp(ia) = ",2i8)') ia, imdtyp(ia)
       end do
       write(nfout,'("!m_IS_init_cnstrnt:  nfcatm = ",i8)') nfcatm
       write(nfout,'("!m_IS_init_cnstrnt:  ia_cnst range = ",2i8)') lbound(ia_cnst,1), ubound(ia_cnst,1)
    end if

    if(cnst_typ == BONDLENGTH_FIX_1 .or. cnst_typ == BONDLENGTH_FIX_2) then
       sgmc(1) = fcv_wk(1,1)
       if(sgmc(1) < DELTA) then
          sgmc(1) = 0.d0
          call cnstrnt(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp&
               &, sgmc,sigma)     !-(b_Ionic_System) ->sigma
          sgmc(1) = sigma
       end if
    else if(cnst_typ == COG_FIX_L) then
       sgmc(1:3) = fcv_wk(1,1:3)
       call cnstrnt(natm,ityp,amion,cps,cnst_typ,nfcatm,ia_cnst,imdtyp &
            &, sgmc, sigma)      !-(b_Ionic_System) ->sigma
       sgmc(4) = sigma
    end if
    if(iprimd>=1) then
       write(nfout,'(" !IS const. of constraint = ",f16.8," <<m_IS_init_cnstrnt>>")') sigma
       write(nfout,'(" !IS sgmc(1:4) = ",4f12.8)') sgmc(1:4)
    end if
  end subroutine m_IS_init_cnstrnt

  subroutine m_IS_wd_cpo_and_forc(nfdynm,forc_l)
    integer, intent(in)                          :: nfdynm
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l
    integer    :: ia
    if(mype == 0) then
       do ia = 1, natm
!xocl spread do/ind_natm
!xocl index ia
          write(nfdynm,'(" ",i4,3f15.9,3f10.5)') ia &
               &, cpo_l(ia,1,1),cpo_l(ia,2,1),cpo_l(ia,3,1) &
               &, forc_l(ia,1), forc_l(ia,2), forc_l(ia,3)
!xocl end spread
       end do
    end if
  end subroutine m_IS_wd_cpo_and_forc

  subroutine m_IS_rd_diis_history(nfcntn,nfree)
    integer, intent(in) :: nfcntn
    integer, intent(in), optional :: nfree
    integer :: nf
    integer :: i,j
    logical             :: tag_is_found, EOF_reach
    integer :: ntmp
    nf = natm
    if(present(nfree)) nf = nfree
    if(.not. (imdalg==GDIIS .or. imdalg==BFGS .or. imdalg==L_BFGS) ) return
    if(printable) write(nfout,'("tag_diis_history")')
    allocate(u_l_buf(nf,3,kqnmditer_p));u_l_buf=0.0d0
    allocate(w_l_buf(nf,3,kqnmditer_p));w_l_buf=0.0d0
    allocate(ncrspd_buf(kqnmditer_p)); ncrspd_buf(:) = (/(i,i=1,kqnmditer_p)/)
    if(mype==0)then
        call rewind_to_tag0(nfcntn,len('(tag_diis_history)'),'(tag_diis_history)' &
    &   ,EOF_reach,tag_is_found,str,len_str)
        if(tag_is_found) then
           read(nfcntn,*)
           read(nfcntn,*) ntmp
           if(ntmp/=kqnmditer_p)then
              if(printable) write(nfout,*) "!!kqnmditer has been changed"
           endif

           if(ntmp==kqnmditer_p)then
              read(nfcntn,*)
              read(nfcntn,*) iter_gdiis

              read(nfcntn,*)
              do i=1,kqnmditer_p
                 do j=1,nf
                    read(nfcntn,*) u_l_buf(j,1,i),u_l_buf(j,2,i),u_l_buf(j,3,i)
                 enddo
              enddo
              read(nfcntn,*)
              do i=1,kqnmditer_p
                 do j=1,nf
                    read(nfcntn,*) w_l_buf(j,1,i),w_l_buf(j,2,i),w_l_buf(j,3,i)
                 enddo
              enddo
              read(nfcntn,*)
              do i=1,kqnmditer_p
                  read(nfcntn,*) ncrspd_buf(i)
              enddo
              diis_continuable = .true.
           endif
        endif
    endif
    if(npes>1) call mpi_bcast(diis_continuable,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(diis_continuable)then
       if(npes>1)then
           call mpi_bcast(u_l_buf,kqnmditer_p*nf*3,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(w_l_buf,kqnmditer_p*nf*3,mpi_double_precision,0,MPI_CommGroup,ierr)
           call mpi_bcast(ncrspd_buf, kqnmditer_p,mpi_integer,0,MPI_CommGroup,ierr)
           call mpi_bcast(iter_gdiis,1,mpi_integer,0,MPI_CommGroup,ierr)
       endif
    endif
  end subroutine m_IS_rd_diis_history

  subroutine m_IS_wd_diis_history(nfcntn,nfree)
    integer, intent(in) :: nfcntn
    integer, intent(in), optional :: nfree
    integer :: nf
    integer :: i,j
    nf = natm
    if(present(nfree)) nf = nfree
    if(.not. (imdalg==GDIIS .or. imdalg==BFGS .or. imdalg==L_BFGS) ) return
    if(.not.allocated(u_l) .or. .not.allocated(w_l) .or. .not.allocated(ncrspd)) then
        if (.not.allocated(u_l_buf) .or. .not.allocated(w_l_buf).or. .not.allocated(ncrspd_buf)) then
            return
        else
            allocate(u_l(nf,3,kqnmditer_p))
            allocate(w_l(nf,3,kqnmditer_p))
            allocate(ncrspd(kqnmditer_p));ncrspd(:) = (/(i,i=1,kqnmditer_p)/)
            u_l = u_l_buf
            w_l = w_l_buf
            ncrspd = ncrspd_buf
        endif
    endif
    if(printable) write(nfout,'(" tag_diis_history")')
    if(mype==0)then
        write(nfcntn,*) '(tag_diis_history)'
        write(nfcntn,*) '(kqnmditer_p)'
        write(nfcntn,*) kqnmditer_p
        write(nfcntn,*) '(iter_gdiis)'
        write(nfcntn,*) iter_gdiis
        write(nfcntn,*) '(u)'
        do i=1,kqnmditer_p
           do j=1,nf
              write(nfcntn,'(3f30.20)') u_l(j,1,i),u_l(j,2,i),u_l(j,3,i)
           enddo
        enddo
        write(nfcntn,*) '(w)'
        do i=1,kqnmditer_p
           do j=1,nf
              write(nfcntn,'(3f30.20)') w_l(j,1,i),w_l(j,2,i),w_l(j,3,i)
           enddo
        enddo
        write(nfcntn,*) '(ncrspd)'
        do i=1,kqnmditer_p
           write(nfcntn,'(i8)') ncrspd(i)
        enddo
    endif
  end subroutine m_IS_wd_diis_history

  subroutine m_IS_wd_cps(nf)
    integer, intent(in) :: nf
    integer :: ia
    if(mype == 0) then
       if(iprimd >= 2 .or. natm <= 100) then
          write(nf,'(" !ion cps and pos at the end of this job")')
          do ia = 1, natm
             write(nf,'(" !ion ",i4,3f14.8,3f12.7)') ia, cps(ia,1),cps(ia,2),cps(ia,3) &
                  & ,pos(ia,1),pos(ia,2),pos(ia,3)
          end do
       end if
    end if
  end subroutine m_IS_wd_cps

  subroutine m_IS_wd_pos_brav(nf)
    use m_Control_Parameters,  only : ipribravpos
    use m_Crystal_Structure,    only : kt_iuctype, a, b, c, ca, cb, cc, il

    integer, intent(in) :: nf
    integer :: ia, j, it
    real(kind=DP) :: work1(3,3), work2(3)

    if ( ipribravpos < 1 ) return
    if ( kt_iuctype /= 0 ) return       ! only for bravias

    if (mype == 0) then
       write(nf,'(a)')
       write(nf,'(" -- Bravais lattice info and internal positions at the end of this job--- ")')
       call primitive2bravais( nf, p2bmat, altv(:,1), altv(:,2), altv(:,3), &
            &                  a, b, c, ca, cb, cc, il )
       write(nf,'(a)')
       write(nf,'(a)') ' --- lattice --- '
       write(nf,'(3(a,f18.12))')  '  a     = ', a, ' , b    = ', b, ' ,  c     = ', c
       write(nf,'(3(a,f18.12))')  '  alpha = ', acos(ca)/PAI*180.0d0, &
            &                           ' , beta = ', acos(cb)/PAI*180.0d0, &
            &                           ' ,  gamma = ', acos(cc)/PAI*180.0d0
       write(nf,'(a)')

       work1(1,:) = altv(:,1);    work1(2,:) = altv(:,2);    work1(3,:) = altv(:,3)
       work1 = matmul( p2bmat, work1 )
       write(nf,'(a,3f20.10)') '    a_vector = ',work1(1,1:3)
       write(nf,'(a,3f20.10)') '    b_vector = ',work1(2,1:3)
       write(nf,'(a,3f20.10)') '    c_vector = ',work1(3,1:3)

       write(nf,*)
       write(nf,'(a)') ' --- ionic positions --- '
       Do ia=1, natm
          it = ityp(ia)
          Do j=1, 3
             work2(j) = sum(b2pmat(:,j)*pos(ia,:))
          end do
          write(nf,'(A4,3f22.18)') speciesname(it), work2(1:3)
       ENd Do
       write(nf,'(a)') ' ---------- '
    end if

  end subroutine m_IS_wd_pos_brav

  subroutine m_IS_rot_ncrspd(nsum)
    integer, intent(out) :: nsum
    integer ::              nbox, istrbr, itemp, it
    nbox = (iter_gdiis-1)/kqnmditer_p
    if(gdiis_hownew == ANEW) then
       istrbr = nbox*kqnmditer_p + 1
    else if (gdiis_hownew == RENEW) then
       if(nbox == 0) then
          istrbr = 1
       else
          istrbr = iter_gdiis - (kqnmditer_p-1)
          itemp = ncrspd(1)
          do it = 1, kqnmditer_p-1
             ncrspd(it) = ncrspd(it+1)
          end do
          ncrspd(kqnmditer_p) = itemp
       end if
    end if
    nsum = iter_gdiis - istrbr + 1
    if(iprigdiis >= 2 .and. printable) then
       write(nfout,'(" -- rot_ncrspd -- ")')
       write(nfout,'("   i  : ",8i8)') (it,it=1,nsum)
       write(nfout,'("ncrspd: ",8i8)') (ncrspd(it),it=1,nsum)
    end if
  end subroutine m_IS_rot_ncrspd

  subroutine m_IS_stor_cps_forc(it,nfree,cps,forc_g)
    integer, intent(in) :: it
    integer, intent(in) :: nfree
    real(kind=DP), intent(in), dimension(nfree,3) :: cps,forc_g
!!$      real(kind=DP) :: xmul
!!$      integer  :: i
    u_l(:,:,it) = cps(:,:)
    w_l(:,:,it) = forc_g(:,:)
  end subroutine m_IS_stor_cps_forc

  subroutine m_IS_do_bfgs(nfout,nsum,nfree,etotal,imdtypxyz,cps,forc_l,factor)
     integer, intent(in) :: nfout
     integer, intent(in) :: nsum
     integer, intent(in) :: nfree
     real(kind=DP), intent(in) :: etotal
     integer, intent(in), dimension(nfree,3) :: imdtypxyz
     real(kind=DP), intent(inout), dimension(nfree,3) :: cps
     real(kind=DP), intent(in), dimension(nfree,3) :: forc_l
     real(kind=DP), intent(in), optional :: factor
     integer :: i,j,k,i1,j1,itr0,itr1,info
     real(DP) :: xgi,gihg
     real(kind=DP), allocatable, dimension(:) :: gdelta,xdelta
     logical :: corrected_eig
     real(DP), allocatable, dimension(:) :: gdotinvh
     real(DP), allocatable, dimension(:,:) :: tmpforc
     real(DP), allocatable, dimension(:,:) :: ihess
     real(DP), allocatable, dimension(:) :: amat
     real(DP), allocatable, dimension(:) :: eigv
     real(DP), allocatable, dimension(:,:) :: eigvec
     real(DP), allocatable, dimension(:) :: workar
     real(DP) :: maxoptforc,tmpmaxoptforc
     real(DP),save :: maxoptforc_old=1000.d0
     real(DP),save :: energy_old=0.0d0
     real(DP),allocatable, dimension(:,:), save :: force_old
     real(DP),save :: alpha_bfgs = 1.0d0
     integer :: nmobile,icounti,icountj
     real(DP) :: c1,c2,mina,maxa,incre
     real(DP) :: tmpsum1,tmpsum2,e0,f0
     integer :: nf
     real(DP) :: fac
     fac = 1.d0
     if(present(factor)) fac = factor
     c1 = wolfe_c1
     c2 = wolfe_c2
     mina = min_alpha
     maxa = max_alpha
     incre = 1.1d0
     nmobile=0
     do i=1,nfree
        do j=1,3
           if(imdtypxyz(i,j)==0)cycle
           nmobile = nmobile+1
        enddo
     enddo
     allocate(xdelta(nmobile));xdelta=0.d0
     allocate(gdelta(nmobile));gdelta=0.d0
     allocate(gdotinvh(nmobile));gdotinvh=0.d0
     allocate(tmpforc(nfree,3));tmpforc=0.d0
     allocate(ihess(nmobile,nmobile));ihess=0.d0
     if(sw_correct_eigenvalue==ON)then
        allocate(amat(nmobile*(nmobile+1)/2));amat=0.d0
        allocate(eigv(nmobile));eigv=0.0d0
        allocate(eigvec(nmobile,nmobile));eigvec=0.0d0
        allocate(workar(3*nmobile));workar=0.0d0
     endif

!      build inverse of the Hessian
     ihess = 0.d0
     do i=1,nmobile
        ihess(i,i) = 1.0d0
     enddo
     do i=2,nsum
        itr1 = ncrspd(i)
        itr0 = ncrspd(i-1)
        icountj=0

        do j=1,nfree
           do k=1,3
              if(imdtypxyz(j,k)==0)cycle
              icountj=icountj+1
              xdelta(icountj) =  u_l(j,k,itr1)-u_l(j,k,itr0)
              gdelta(icountj) = -w_l(j,k,itr1)+w_l(j,k,itr0)
           enddo
        enddo
        xgi = 1.0d0/dot_product(xdelta,gdelta)
        if(xgi<0)then
           if(printable .and. iprigdiis>=2)then
              write(nfout,'(a,i3)') '!** WARNING dx dot dg is negative for history ',itr1
              write(nfout,'(a)') 'skipping this update.'
           endif
           cycle
        endif
        do j=1,nmobile
           gdotinvh(j) = dot_product(ihess(j,:),gdelta(:))
        enddo
        gihg = dot_product(gdelta,gdotinvh)
        do j=1,nmobile
           do k=1,nmobile
              ihess(j,k) = ihess(j,k)+xgi*xgi*(1.0d0/xgi+gihg)*xdelta(j)*xdelta(k) &
        &                - (gdotinvh(j)*xdelta(k)+gdotinvh(k)*xdelta(j))*xgi
           enddo
        enddo
     enddo

!      correct bad eigenvalues present in the Hessian
     if(sw_correct_eigenvalue==ON)then
        corrected_eig=.false.
        do i=1,nmobile
           do j=i,nmobile
              amat(i + (j-1)*j/2) = ihess(i,j)
           enddo
        enddo
        call dspev('V','U',nmobile,amat,eigv,eigvec,nmobile,workar,info)
        if(printable .and. iprigdiis>=2) write(nfout,'(a)') '--- eigenvalues for the approximate Hessian ---'
        do i=1,nmobile
           if(printable.and.iprigdiis>=2) write(nfout,'(i8,f20.10)') i,1.0d0/eigv(i)
           if (1.0d0/eigv(i)<eigenvalue_threshold)then
              eigv(i) = 1.0d0/eigenvalue_threshold
              if(printable.and.iprigdiis>=2) write(nfout,'(a,i8,a,f20.10)') &
              &  'corrected the eigenvalue for the ',i,'-th element to : ',1.0d0/eigv(i)
              corrected_eig=.true.
           endif
        enddo
        if(corrected_eig)then
           ihess=0.d0
           do i=1,nmobile
              do j=1,nmobile
                 do k=1,nmobile
                    ihess(i,j) = ihess(i,j)+eigvec(i,k)*eigvec(j,k)*eigv(k)
                 enddo
              enddo
           enddo
        endif
     endif

!      H^-1 dot g
     maxoptforc=0.d0
     icounti=0
     do i=1,nfree
        do j=1,3
           tmpforc(i,j) = 0.d0
           if(imdtypxyz(i,j)==0)cycle
           icounti=icounti+1
           icountj=0
           do i1=1,nfree
              do j1=1,3
                 if(imdtypxyz(i1,j1)==0)cycle
                 icountj=icountj+1
                 tmpforc(i,j) = tmpforc(i,j) - ihess(icounti,icountj)*forc_l(i1,j1)
              enddo
           enddo
           tmpmaxoptforc=dsqrt(dot_product(tmpforc(i,:),tmpforc(i,:)))
           if(tmpmaxoptforc>maxoptforc)maxoptforc=tmpmaxoptforc
        enddo
     enddo

     if(printable .and. ipri>=1) then
        write(nfout,'(a,f20.10)') 'max. optimal force obtained from the BFGS update : ',maxoptforc
     endif

     if(sw_optimize_alpha == ON) then
     if(.not.allocated(force_old)) then
        allocate(force_old(nmobile,3))
     else
        tmpsum1 = 0.d0
        tmpsum2 = 0.d0
        do i=1,nfree
           do j=1,3
              if (imdtypxyz(i,j)==0) cycle
              tmpsum1 = tmpsum1 + force_old(i,j)*tmpforc(i,j)
              tmpsum2 = tmpsum2 + forc_l(i,j)*tmpforc(i,j)
           enddo
        enddo
        e0 = energy_old + c1*alpha_bfgs*tmpsum1
        f0 = -c2*tmpsum1
        if(e0>etotal .and. f0>abs(tmpsum2))then
           alpha_bfgs = alpha_bfgs*incre
        else
           alpha_bfgs = alpha_bfgs/incre
        endif
        if(alpha_bfgs>maxa) alpha_bfgs = maxa
        if(alpha_bfgs<mina) alpha_bfgs = mina
        if(printable)then
           write(nfout,'(a,2f20.10)') 'Wolfe condition  i) : ',e0,etotal
           write(nfout,'(a,2f20.10)') 'Wolfe condition ii) : ',f0,abs(tmpsum2)
           write(nfout,'(a,f20.10)')  'new alpha           : ',alpha_bfgs
        endif
     endif

     energy_old = etotal
     force_old(:,:)  = forc_l(:,:)
     endif

     if(maxoptforc_old*10<maxoptforc)then
        if (printable) write(nfout,'(a)') 'the estimated force seems to be very large; &
         & update will be done by the steepest-descent method'
        do i=1,nfree
           do j=1,3
              if (imdtypxyz(i,j)==0) cycle
              cps(i,j) = cps(i,j)+forc_l(i,j)
           enddo
        enddo
     else
        tmpsum1=0.d0
        tmpsum2=0.d0
        do i=1,nfree
           do j=1,3
              if (imdtypxyz(i,j)==0) cycle
              cps(i,j) = cps(i,j)-alpha_bfgs*tmpforc(i,j)*fac
           enddo
        enddo
     endif
     maxoptforc_old = maxoptforc

     deallocate(xdelta)
     deallocate(gdelta)
     deallocate(gdotinvh)
     deallocate(tmpforc)
     deallocate(ihess)
     if(sw_correct_eigenvalue==ON)then
        deallocate(amat)
        deallocate(eigv)
        deallocate(eigvec)
        deallocate(workar)
     endif
  end subroutine m_IS_do_bfgs

  subroutine m_IS_do_lbfgs(nfout,nsum,nfree,etotal,imdtypxyz,cps,forc_l)
     integer, intent(in) :: nfout
     integer, intent(in) :: nsum
     integer, intent(in) :: nfree
     real(kind=DP), intent(in) :: etotal
     integer, intent(in), dimension(nfree,3) :: imdtypxyz
     real(kind=DP), intent(inout), dimension(nfree,3) :: cps
     real(kind=DP), intent(in), dimension(nfree,3) :: forc_l
     real(kind=DP), allocatable, dimension(:) :: rhok,alphak,betak
     real(kind=DP), allocatable, dimension(:,:,:) :: sk,yk
     real(kind=DP), allocatable, dimension(:,:) :: q,z,qq
     real(kind=DP), save :: alpha_bfgs = 1.0d0
     real(kind=DP), save :: energy_old
     real(kind=DP), allocatable, dimension(:,:), save :: force_old,cps_old
     real(kind=DP) :: c1
     real(kind=DP), allocatable, dimension(:,:) :: pmat,forcw
     real(kind=DP) :: dotp,H0,gamm,dotp2,fpk,energy_p,maxv,fac,diff,maxoptfor,av0,av1,fsum
     integer, allocatable, dimension(:) :: ipfrc
     integer :: i,j,k,l,inf,nmobile
     c1 = 0.1d0
     maxv = maxstep
     if(nsum>=2)then
       allocate(rhok(2:nsum))
       allocate(alphak(2:nsum))
       allocate(betak(2:nsum))
       allocate(sk(nfree,3,2:nsum))
       allocate(yk(nfree,3,2:nsum))
     endif
     allocate(q(nfree,3));q=-forc_l
     allocate(z(nfree,3))
     if(.not.allocated(force_old)) allocate(force_old(nfree,3))
     if(.not.allocated(cps_old)) allocate(cps_old(nfree,3))

     do i=2,nsum
        yk(:,:,i) = -w_l(:,:,i)+w_l(:,:,i-1)
        sk(:,:,i) =  u_l(:,:,i)-u_l(:,:,i-1)
        dotp = 0.d0
        do j=1,3
           do k=1,nfree
              if(imdtypxyz(k,j)==0)cycle
              dotp = dotp+yk(k,j,i)*sk(k,j,i)
           enddo
        enddo
        rhok(i) = 1.d0/dotp
     enddo
     gamm=1.d0
     if(nsum>=2)then
       dotp = 0.d0
       dotp2 = 0.d0
       do i=1,3
          do j=1,nfree
            if(imdtypxyz(j,i)==0)cycle
            dotp = dotp+sk(j,i,nsum)*yk(j,i,nsum)
            dotp2 = dotp2+yk(j,i,nsum)*yk(j,i,nsum)
          enddo
       enddo
       gamm = dotp/dotp2
       if(sw_prec==ON)then
         allocate(pmat(nfree,nfree));pmat=0.d0
         call m_IS_precon_matrix(nfree,forc_l-force_old,cps-cps_old,precon_A,pmat,gamm)
       endif
     endif
     H0 = gamm
     if(printable .and. ipri>1) write(nfout,*) 'gamm ',gamm
     do i=nsum,2,-1
        dotp = 0.d0
        do j=1,3
           do k=1,nfree
              if(imdtypxyz(k,j)==0)cycle
              dotp = dotp+q(k,j)*sk(k,j,i)
           enddo
        enddo
        alphak(i) = rhok(i)*dotp
        q(:,:) = q(:,:)-alphak(i)*yk(:,:,i)
     enddo
!     if(nsum==1)then
     if(sw_prec==OFF .or. nsum==1)then
        z(:,:) = q(:,:)*H0
     else
        allocate(qq(nfree,3));qq=q
        call dposv('U',nfree,3,pmat,nfree,q,nfree,inf)
        if(inf /= 0) then
          write(nfout,'(a,i8,a)') '!** info from dposv : ',inf, &
          & ' using the default value for the initial estimate of the Hessian'
          z(:,:) = qq(:,:)*H0
        else
!          av0=0.d0;av1=0.d0
!          do j=1,3
!            do i=1,nfree
!                av0 = av0+qq(i,j)*qq(i,j)
!                av1 = av1+q(i,j)*q(i,j)
!            enddo
!          enddo
!          av0=sqrt(av0)
!          av1=sqrt(av1)
!          z(:,:) = H0*(av0/av1)*q(:,:)
          z(:,:) = q(:,:)
        endif
        deallocate(qq)
     endif
     do i=2,nsum
        dotp = 0.d0
        do j=1,3
           do k=1,nfree
              if(imdtypxyz(k,j)==0)cycle
              dotp = dotp+yk(k,j,i)*z(k,j)
           enddo
        enddo
        betak(i) = rhok(i)*dotp
        z(:,:) = z(:,:)+(alphak(i)-betak(i))*sk(:,:,i)
     enddo
     if(nsum>=2)then
       fpk=0.d0
       do i=1,3
          do j=1,nfree
             if(imdtypxyz(j,i)==0)cycle
             fpk = fpk+force_old(j,i)*z(j,i)
          enddo
       enddo
       energy_p = energy_old+alpha_bfgs*c1*fpk
       diff=energy_p-etotal
       if(printable .and. ipri>=1) write(nfout,'(a,2f20.10,l2)') 'energy pred, energy cal : ',energy_p,etotal,diff>=0
       if(diff<0)then
         if(printable .and. ipri>=1) write(nfout,'(a,f20.10)') 'candidate alpha ', &
        & -0.5*alpha_bfgs*fpk/((etotal-energy_old)/alpha_bfgs - fpk)
       endif
     endif
     energy_old = etotal
     force_old = forc_l
     cps_old = cps

     allocate(forcw(natm,3)); allocate(ipfrc(natm2))
     call fd_symmetrize(natm2,natm,natm,napt,nopr+af,nopr,op,iwei &
          &, z, forcw, ipfrc) ! -(bottom_Subroutines_para)
     deallocate(forcw)
     deallocate(ipfrc)
     fac = 1.d0
     maxoptfor = maxval(abs(z))
     if(alpha_bfgs*maxoptfor>maxv)then
       fac = maxv/(alpha_bfgs*maxoptfor)
     endif
     if(printable .and. ipri>=1) write(nfout,'(a,2f20.10)') 'maxval and factor ',maxoptfor,fac
     if(nopr==1) then
       do i=1,3
          fsum = 0.d0
          nmobile = 0
          do j=1,nfree
            if(imdtypxyz(j,i)==0)cycle
            fsum = fsum+z(j,i)
            nmobile = nmobile+1
          enddo
          fsum = fsum/dble(nmobile)
          do j=1,nfree
            if(imdtypxyz(j,i)==0)cycle
            z(j,i) = z(j,i)-fsum
          enddo
       enddo
     endif
     do i=1,3
        do j=1,nfree
           if(imdtypxyz(j,i)==0)cycle
           cps(j,i) = cps(j,i) - alpha_bfgs * z(j,i) * fac
        enddo
     enddo
     if(nsum>=2)then
       deallocate(rhok)
       deallocate(alphak)
       deallocate(betak)
       deallocate(sk)
       deallocate(yk)
       if(sw_prec==ON) deallocate(pmat)
     endif
     deallocate(q)
     deallocate(z)
  end subroutine m_IS_do_lbfgs

  subroutine m_IS_precon_matrix(nfree,difff,diff,Aval,pmat,gamm)
     integer, intent(in) :: nfree
     real(kind=DP), dimension(nfree,3), intent(in) :: difff,diff
     real(kind=DP), intent(in) :: Aval,gamm
     real(kind=DP), dimension(nfree,nfree), intent(out) :: pmat
     real(kind=DP) :: r2,cij
     real(kind=DP), dimension(3) :: dr
     real(kind=DP), save :: rnn,rcut,mu
     real(kind=DP) :: rmin
     real(kind=DP), allocatable, dimension(:) :: rmins
     real(kind=DP), allocatable, dimension(:) :: covrads
     real(kind=DP), allocatable, dimension(:,:) :: tmpar
     integer       :: neibrd, alen(3)
     real(kind=DP), allocatable :: rxyz(:,:)    ! d(neibrd,3)
     real(kind=DP), allocatable :: rr(:)        ! d(neibrd)
     real(kind=DP), pointer, dimension(:,:) :: cps_fp
     integer,       pointer, dimension(:)   :: ityp_full ! d(natm2)
     integer :: i,j,k,ierr,nu,in,ia,ia2
     logical, save :: firstcall=.true.
     real(kind=DP) :: rhs,lhs,diag,sdiag
     real(kind=DP), allocatable, dimension(:,:) :: pmatt
     real(kind=DP), allocatable, dimension(:) :: work
     integer,allocatable, dimension(:) :: ipiv
     integer :: inf

     allocate(cps_fp(natm2,3))
     allocate(ityp_full(natm2))
     cps_fp(1:natm,1:3) = pos(1:natm,1:3)
     ityp_full(1:natm)  = ityp(1:natm)
     call rplcps(cps_fp,ityp_full,1,natm2,natm,iwei)
     call cpspac ! -> cps_fp

     call decide_rxyz_size(12.5d0,alen,neibrd) !-> alen, neibrd
     allocate(rxyz(neibrd,3))
     allocate(rr(neibrd))
     call substitute_rxyz(alen, neibrd, rxyz, rr)

!     if(firstcall)then
       allocate(rmins(ista_atm:iend_atm))
       do nu = ista_atm, iend_atm ! MPI
          rmin = 1.d+30
          do in = 1, neibrd
             do ia=1,natm2
                if(in == 1 .and. ia == nu) cycle
                dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
                r2 = dot_product(dr,dr)
                if(r2<rmin) rmin=r2
             enddo
          enddo
          rmins(nu) = rmin
       enddo
       rnn = maxval(rmins)
       deallocate(rmins)
       call mpi_allreduce(MPI_IN_PLACE,rnn,1,mpi_double_precision,mpi_max, &
        &   mpi_ke_world,ierr)
       if(rnn<1.d0) then
          write(nfout,'(a,f10.3)') '!** WARNING longest nn distance is very short ',sqrt(rnn)
          rnn = 1.d0
       endif
       rnn = sqrt(rnn)
       rcut = 1.1d0*rnn
       if(Aval.gt.0) then
         rcut = 2*rnn
       endif
       if(ipri>1) write(nfout,'(a,2f20.10,2i8)') 'rnn, rcut : ',rnn,rcut

       pmat = 0.d0
       do nu = ista_atm, iend_atm ! MPI
          do in = 1, neibrd
             do ia=1,natm2
                if(in == 1 .and. ia == nu) cycle
                dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
                r2 = dot_product(dr,dr)
                if (r2>(rcut*rcut)) cycle
                ia2 = ia
                if(ia>natm) ia2 = ia-natm
                pmat(nu,ia2) = -get_cij(Aval,sqrt(r2),rnn)
             enddo
          enddo
       enddo
       call mpi_allreduce(MPI_IN_PLACE,pmat,nfree*nfree,mpi_double_precision, &
          & mpi_sum,mpi_ke_world,ierr)

       sdiag = 0.d0
       do ia=1,natm
          diag = 0.d0
          do ia2=1,natm
             if(ia==ia2) cycle
             diag = diag-pmat(ia,ia2)
          enddo
!          pmat(ia,ia) = diag+0.1d0
          pmat(ia,ia) = diag
          sdiag = sdiag+diag
       enddo

!     allocate(pmatt(nfree,nfree));pmatt=pmat
!     allocate(ipiv(nfree));ipiv=0
!     allocate(work(nfree*4));work=0.d0
!     call dgetrf(nfree,nfree,pmatt,nfree,ipiv,inf)
!     if(ipri>1) write(nfout,'(a,i8)') 'matinv info0 ',inf
!     call dgetri(nfree,pmatt,nfree,ipiv,work,nfree*4,inf)
!     if(ipri>1) then
!     write(nfout,'(a,i8)') 'matinv info1 ',inf
!     write(nfout,'(a)') 'precon matrix'
!     write(nfout,'(a,f20.10)') 'diag av ',sdiag/dble(natm)
!     do ia=1,nfree
!        do ia2=1,nfree
!           write(nfout,*) ia,ia2,pmat(ia,ia2)
!        enddo
!     enddo
!     endif
!     deallocate(pmatt)
!     deallocate(ipiv)
!     deallocate(work)
!
!       diag = 0.d0
!       do ia=1,natm
!          diag = diag+pmat(ia,ia)
!       enddo
!       write(nfout,'(a,f20.10)') 'diag ',diag
!       mu = gamm/diag
!       mu = 1.d0/diag
!       write(nfout,'(a,f20.10)') 'mu ',mu

!       allocate(tmpar(natm,3));tmpar=0.d0
!       do i=1,natm
!          tmpar = 0.d0
!          do j=1,natm
!             do k=1,3
!                tmpar(i,k) = tmpar(i,k)+pmat(i,j)*diff(j,k)
!             enddo
!          enddo
!          rhs = 0.d0
!          do j=1,natm
!             do k=1,3
!                rhs = rhs + tmpar(j,k)*diff(j,k)
!             enddo
!          enddo
!       enddo
!       deallocate(tmpar)
!
!       lhs = 0.d0
!       do j=1,natm
!          do k=1,3
!             lhs = lhs+difff(j,k)*diff(j,k)
!          enddo
!       enddo
!       mu = 0.01d0*abs(lhs/rhs)
       !mu = precon_mu
       if(precon_mu>0) then
         mu = precon_mu
       else
!         mu = abs(lhs/rhs/gamm/35.d0)
!         mu = abs(lhs/rhs)
!         if(mu>1) mu=1.d0
!         mu = 1.d0/abs(gamm)
!         mu = 0.16d0/(sdiag/dble(natm))
         mu = abs(1.d0/(sdiag/dble(natm2))/gamm)
         if(printable .and. ipri>0) write(nfout,'(a,f20.10)') ' !** PRECON estimated value for mu :',mu
       endif
!       write(nfout,'(a,4f20.10)') 'rhs,lhs,lhs/rhs,mu',rhs,lhs,abs(lhs/rhs),mu
!       firstcall = .false.
!     endif

     pmat = 0.d0
     do nu = ista_atm, iend_atm ! MPI
        do in = 1, neibrd
           do ia=1,natm2
              if(in == 1 .and. ia == nu) cycle
              dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
              r2 = dot_product(dr,dr)
              if (r2>(rcut*rcut)) cycle
              ia2 = ia
              if(ia>natm) ia2 = ia-natm
              pmat(nu,ia2) = -get_cij(Aval,sqrt(r2),rnn)*mu
           enddo
        enddo
     enddo
     call mpi_allreduce(MPI_IN_PLACE,pmat,nfree*nfree,mpi_double_precision, &
        & mpi_sum,mpi_ke_world,ierr)
     do ia=1,natm
        diag = 0.d0
        do ia2=1,natm
           if(ia==ia2) cycle
           diag = diag-pmat(ia,ia2)
        enddo
        pmat(ia,ia) = diag+mu*0.1d0
     enddo

     deallocate(rr)
     deallocate(rxyz)
     deallocate(cps_fp)
     deallocate(ityp_full)

     contains

     function get_cij(Aval,rij,rnn) result(ret)
       real(kind=DP), intent(in) :: Aval,rij,rnn
       real(kind=DP) :: ret
       ret = exp(-Aval*(rij/rnn-1.d0))
       return
     end function get_cij

     subroutine cpspac
       real(kind=DP), dimension(3) :: catoms(3)
       integer                     :: i
       do i = 1, natm2
          catoms = cps_fp(i,1:3)
          catoms = catoms - nint(catoms)      !Packing
          cps_fp(i,1:3) = matmul(altv,catoms) !Change of coordinate system
       end do
     end subroutine cpspac

  end subroutine m_IS_precon_matrix

  subroutine m_IS_gdiis(forc_l_in,forcmx,etotal)
!!$   Original program was coded by T. Yamasaki(JRCAT-ATP)
!!$                            @(#)md_qn_gdiis.f 9.1 97/05/08 14:48:50
!    Transformed into this using fortran90
!                 by T. Yamasaki (FUJITSU Laboratories Ltd.)
!                                   29th Apr. 2003
!
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP), intent(in)                    :: forcmx, etotal
    real(kind=DP), allocatable, dimension(:,:)   :: forc_l

    integer ::              m, mode_init = 0,sub_optmode,it,nsum,icount,icon, ia, imd_t, ia_t,j
    logical ::              m_IS_md_is_executed
    real(kind=DP) ::        c_dist, forcmx_mdfy
    integer ::              id_sname = -1
    call tstatc0_begin('m_IS_gdiis ',id_sname)

    if(absolute_convergence_of_forc(forc_l_in))then
       write(nfout,'(a)') ' m_IS_gdiis: forces are absolutely converged!! nothing to do...'
       return
    endif

    if(num_planes_atoms_are_fixed>=1) then
       allocate(tmass(num_planes_atoms_are_fixed));tmass= 0.d0
       allocate(rtmass(num_planes_atoms_are_fixed)); rtmass = 0.d0
    end if
    if(num_planes_atoms_are_fixed_cog>=1) then
       allocate(fcg_mdfy_cog(num_planes_atoms_are_fixed_cog,1:3))
       allocate(fcg_cog(num_planes_atoms_are_fixed_cog,1:3))
       allocate(tmass_cog(num_planes_atoms_are_fixed_cog)) ; tmass_cog = 0.d0
       allocate(rtmass_cog(num_planes_atoms_are_fixed_cog)); rtmass_cog = 0.d0

    end if
    if(num_planes_atoms_are_fixed_rb>=1) then
       allocate(fcg_rb(num_planes_atoms_are_fixed_rb,1:3))
       allocate(fcg_mdfy_rb(num_planes_atoms_are_fixed_rb,1:3))
       allocate(tmass_rb(num_planes_atoms_are_fixed_rb)); tmass_rb = 0.d0
       allocate(rtmass_rb(num_planes_atoms_are_fixed_rb)); rtmass_rb = 0.d0
    end if

    if(mdmode == ORDINA) then
       allocate(forc_l(natm,3))
       forc_l = forc_l_in

!!$       if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_fi(forc_l)
       if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_hyperplane(forc_l)
       iter_gdiis = iter_gdiis + 1
       if(if_allocated == 0) call m_IS_gdiis_alloc(mode_init,if_allocated)

       forc_g = forc_l
!!$       call mdfy_forc(forcmx,forcmx_mdfy)  ! -> forc_g
       m = mod(iter_gdiis,kqnmditer_p)
       if(m == 0) m = kqnmditer_p

       call d_fc_l(mode_init,forcmx_mdfy) ! -> fc_l
       call m_IS_rot_ncrspd(nsum)  ! -> ncrspd, nsum
       call d_sub_optmode(etotal,etot_previous,forcmx_mdfy,sub_optmode) ! -> sub_optmode
       it = ncrspd(nsum)
       if(it == 0) it = 1
       call m_IS_stor_cps_forc(it,natm,cps,forc_g) ! cps,forc_g -> u_l, w_l
       if(iprigdiis >= 2 .and. printable) write(nfout,'(" sub_optmode = ", i5)') sub_optmode
       if(sub_optmode == QUENCHED_MD) then
          call m_IS_md(sub_optmode,forc_l)
       else if (imdalg==GDIIS)then
          m_IS_md_is_executed = .false.
          if(nsum > 1) then
             call gdiis_mat(nsum) ! ncrspd,w_l -> f_gdiis
             call getmatinv(nsum,icon) ! -> f_gdiis := f_gdiis^(-1)
             if(icon == 0) then
                call gdiis_g(nsum)   ! -> g
                call opt_forc(nsum)  ! -> forc_g := sum_{it}(g(it)*w_l(,ncrspd(it))
!!$                call mdfy_forc(forcmx,forcmx_mdfy2)
                          if(iprigdiis >= 2 .and. printable) call forc_check(" --- new force ---")
                          if(iprigdiis >= 1) call cmp_fop_fin(it)
                call opt_geom(nsum)  ! -> cps := sum_{it}(g(it)*u_l(,ncrspd(it))
             else   ! when icon /= 0, f_gdiis matrix is singular
                if(printable) write(nfout,'(" !! f_gdiis matrix is singular <<m_IS_gdiis>>")')
                call m_IS_md(QUENCHED_MD,forc_l) ; m_IS_md_is_executed = .true.
             end if
          end if

          if(m_IS_md_is_executed) goto 1002
! --> T. Yamasaki, 19 Mar 2007
          do ia = 1, natm
             do j=1,3
             imd_t = imdtypxyz(ia,j)
             if(imd_t == 0 .or. imd_t == COG_FIX .or. &
                  & imd_t == COG_FIX_L .or. imd_t == FIX_IN_A_PLANE) cycle
             if(imd_t == RIGID_BODY_FIX_IN_A_PLANE) cycle
             cpd_l(ia,j) = -fc_l(ia,j)/dtio*forc_g(ia,j)
             cps(ia,j) = cps(ia,j) + dtio*cpd_l(ia,j)
             enddo
          end do

          select case(constraint_type)
          case (COG_FIX_L, FIX_IN_A_PLANE, RIGID_BODY_FIX_IN_A_PLANE, COG_and_RIGID_BODY_FIX_L)
             if(nfcatm >= 1) allocate(cpd_old(nfcatm,3))
             ia_t = 0
             do ia = 1, natm
                select case(imdtyp(ia))
                case (COG_FIX, COG_FIX_L, FIX_IN_A_PLANE, RIGID_BODY_FIX_IN_A_PLANE)
                   ia_t = ia_t + 1
                   cpd_old(ia_t,:) = cpd_l(ia,:)
                end select
             end do

!!$             select case(constraint_type)
!!$             case (COG_FIX_L, FIX_IN_A_PLANE, COG_and_RIGID_BODY_FIX_L)
             call check_constraint_cog_or_rb(forc_l)

             call evolve_velocities_cog_or_rb(forc_l)

             if(constraint_type==COG_FIX_L .or. constraint_type==COG_and_RIGID_BODY_FIX_L &
                  & .or. constraint_type==RIGID_BODY_FIX_IN_A_PLANE) call correct_cog_andor_rb_motion()

             select case(constraint_type)
             case (COG_FIX_L, FIX_IN_A_PLANE, COG_and_RIGID_BODY_FIX_L)
                call quench_velocities_cog(forc_l)
             case (RIGID_BODY_FIX_IN_A_PLANE)
                call quench_velocities_rb(forc_l)
             end select

             if(nfcatm >= 1) deallocate(cpd_old)
             do ia = 1, natm
                select case(imdtyp(ia))
                case (COG_FIX, COG_FIX_L, FIX_IN_A_PLANE, RIGID_BODY_FIX_IN_A_PLANE )
                   cps(ia,:) = cps(ia,:) + dtio*cpd_l(ia,:)
                   if(iprimd >= DEBUGPRINTLEVEL) write(nfout,'(" cps(",i5,",:) = ",3f8.4)') ia, cps(ia,1:3)
                end select
             end do
          case default
          end select

!!$          if(apply_constant_move >= 1)  call evolve_cps_constant_move() ! -> cps

1002      continue
! <---
          if(constraint_type == COG_FIX_L) call check_cog()
       else if(imdalg==BFGS) then
          call m_IS_do_bfgs(nfout,nsum,natm,etotal,imdtypxyz,cps,forc_l)
       else if(imdalg==L_BFGS) then
!          call do_lbfgs()
          call m_IS_do_lbfgs(nfout,nsum,natm,etotal,imdtypxyz,cps,forc_l)
       end if
       if(imdalg==GDIIS)then
          if(iprigdiis >= DEBUGPRINTLEVEL) call cps_check(" --- new coordinates ---") ! cps

          c_dist = 0.3d0
          call cps_damp(c_dist,it,icount)
          if(icount >= 1 .and. iprigdiis >= 1 .and. printable) then
             call cps_check(" --- new (damped) coordinates ---") ! cps
          endif
       endif

       if(iprigdiis >= 2 .and. printable) call cpd_check(" --- velocity ---") ! cpd_l

       deallocate(forc_l)
    else if(mdmode == CNSTRA) then
       call evolve_cps_constant_move()
!!$       call move_atoms_normal_to_plane()
       iter_gdiis = 0
       if(if_allocated == 1) call m_IS_gdiis_dealloc(if_allocated)
    else
       call phase_error_with_msg(nfout, ' Invalid value of mdmode <<m_IS_gdiis>>',__LINE__,__FILE__)
    end if

    if(num_planes_atoms_are_fixed_cog>=1) deallocate(fcg_mdfy_cog,fcg_cog,rtmass_cog,tmass_cog)
    if(num_planes_atoms_are_fixed_rb>=1)  deallocate(fcg_mdfy_rb, fcg_rb, rtmass_rb, tmass_rb)
    if(num_planes_atoms_are_fixed>=1)  deallocate(rtmass,tmass)

    call tstatc0_end(id_sname)
  contains

    subroutine d_fc_l(mode_init,forcmx)
      integer, intent(in) ::        mode_init
      real(kind=DP), intent(in) :: forcmx
      integer ::             ia, j
      do ia = 1, natm
         do j = 1, 3
            if(imdtypxyz(ia,j) == 0) then
               fc_l(ia,j) = 0.d0
            else if(forcmx >= c_forc_prop_region_high) then
               if(mode_init == UNIT) then
                  fc_l(ia,j) = -1.d0
               else
                  fc_l(ia,j) = -dtio*dtio/amion(ityp(ia))
               end if
            else if(forcmx >= c_forc_prop_region_low) then
               fc_l(ia,j) = - factor_prop_region/forcmx
               if(dabs(fc_l(ia,j)) < 1.d0) then
                  fc_l(ia,j) = -1.d0
               end if
            end if
         end do
      end do
    end subroutine d_fc_l

    subroutine d_sub_optmode(etotal,etot_previous,forcmx_mdfy,sub_optmode)
      real(kind=DP),intent(in) :: etotal, etot_previous
      real(kind=DP),intent(in) :: forcmx_mdfy
      integer, intent(out) ::     sub_optmode
      real(kind=DP) ::            edel
!!$      edel = etotal - etot_previous
!!$      if(edel > 0.d0 .and. forcmx <= c_forc2GDIIS) then
!!$         iincre_at_forc_cal = iincre_at_forc_cal + 1
!!$         if(iincre_at_forc_cal >= ic_E_overshoot) then
!!$            sub_optmode = GDIIS
!!$         end if
!!$      else if(sub_optmode /= GDIIS) then
!!$         sub_optmode = QUENCHED_MD
!!$      end if
      sub_optmode = GDIIS
!!$      if(sub_optmode == GDIIS) then
!!$	if(ipriinputfile >= 1) write(nfout,'(" !! sub_optmode = GDIIS")')
!!$      else if(sub_optmode == QUENCHED_MD) then
!!$         if(ipriinputfile >= 1) write(nfout,'(" !! sub_optmode = QUENCHED_MD")')
!!$      end if
    end subroutine d_sub_optmode

    subroutine gdiis_mat(nsum)
      integer, intent(in) :: nsum

      integer ::             it, jt, itcrspd, jtcrspd,ia,j, imd_t
      real(kind=DP) ::       xmul

      do it = 1, nsum
         itcrspd = ncrspd(it)
         do jt = it, nsum
            jtcrspd = ncrspd(jt)
            xmul = 0.d0
            do ia = 1, natm
               do j = 1, 3
               imd_t = imdtypxyz(ia,j)
! --> T. Yamasaki, 19 Mar 2007
               if(imd_t /= 0 .and. imd_t /= COG_FIX .and. imd_t /= COG_FIX_L &
                    & .and. imd_t /= FIX_IN_A_PLANE .and. imd_t /= RIGID_BODY_FIX_IN_A_PLANE) then
! <---
                     xmul = xmul + w_l(ia,j,jtcrspd)*w_l(ia,j,itcrspd)
               end if
               enddo
            end do
            f_gdiis(it,jt) = xmul
            if(jt /= it) f_gdiis(jt,it) = f_gdiis(it,jt)
         end do
         if(iprigdiis >= 2 .and. printable) then
            if(it == 1) write(nfout,'(" -- f_gdiis -- ")')
            write(nfout,'( 6d12.4)') (f_gdiis(it,jt),jt=1,nsum)
         end if
      end do
    end subroutine gdiis_mat

    subroutine getmatinv(nsum,icon)
      integer ,intent(in) ::  nsum
      integer, intent(out) :: icon

      real(kind=DP) :: div
      integer ::       it, jt, ipfr, ipto
#ifdef _GDIIS_MAT_CHECK_
      integer ::       kt, ik_count, kj_count
#endif

      div = 1.0/f_gdiis(1,1)
      do it = 1, nsum
         do jt = 1, nsum
            ipto = it + (jt-1)*nsum
            f_wk(ipto, 1) = f_gdiis(it,jt)*div
#ifdef _GDIIS_MAT_CHECK_
            f_wk(ipto, 2) = f_gdiis(it,jt)*div
#endif
            if(it == jt) then
               e_wk(ipto) = 1.d0
            else
               e_wk(ipto) = 0.d0
            end if
         end do
      end do

      call rdecomp(nsum,f_wk(1,1),ww1,ip,icon)

      if(icon /= 0) then
         if(printable) then
            write(nfout,'("  [f_wk] after rdecomp <<m_IS_gdiis_getmatinv>>")')
            do it = 1, nsum
               write(nfout,'(" f_wk(:,",i3," ) = ",8d12.4)') it,(f_wk(it+(jt-1)*nsum,1),jt=1,nsum)
            end do
            write(nfout,'(" ip = ",8i6)') (ip(ipto),ipto=1,nsum)
            write(nfout,*) ' LU decomposition is impossible. <<m_IS_gdiis.getmatinv>>'
         end if
         return
!!$         stop ' LU decompoition is impossible <<m_IS_gdiis_getmatinv>>'
      else
         call rsolve(nsum,nsum,f_wk(1,1),e_wk,f_rslv,ip)
      endif
#ifdef _GDIIS_MAT_CHECK_
! ----------- checking of inversion matrix ------------>
      if(printable) then
         write(nfout,*) ' -- below should equal to a unit matrix --'
         do it = 1, nsum
            ww1 = 0.d0
            do jt = 1, nsum
               do kt = 1, nsum
                  ik_count = (kt-1)*nsum + it
                  kj_count = (jt-1)*nsum + kt
                  ww1(jt) = ww1(jt) + f_wk(ik_count,2)*f_rslv(kj_count)
               enddo
            enddo
            write(nfout,9009) it,(ww1(jt),jt=1,nsum)
         enddo
9009     format(i3,8f8.4,/,3x,8f8.4,/,3x,8f8.4,/3x,8f8.4)
      end if
! <--------------------------------------------------
#endif

      do jt = 1, nsum
         do it = 1, nsum
            ipfr = it + (jt-1)*nsum
            f_gdiis(it,jt) = f_rslv(ipfr)
         end do
      end do

      if(iprigdiis >= 2 .and. printable) then
         write(nfout,*) ' **inverse matrix**'
         do jt = 1, nsum
            write(nfout,9008) jt,(f_gdiis(it,jt),it=1,nsum)
         enddo
      end if
 9008 format(i3,(6d12.4))

    end subroutine getmatinv

    subroutine gdiis_g(nsum)
      integer, intent(in) :: nsum
      real(kind=DP) ::       alpha_fac
      integer ::             it,jt
      alpha_fac = 0.d0
      do it = 1, nsum
         do jt = 1, nsum
            alpha_fac = alpha_fac + f_gdiis(it,jt)
         end do
      end do
      alpha_fac = 1.d0/alpha_fac

      g = 0.d0
      do it = 1, nsum
         do jt = 1, nsum
            g(it) = g(it) + alpha_fac*f_gdiis(it,jt)
         end do
      end do

      if(iprigdiis >= 2 .and. printable) then
         write(nfout,*) ' ---- alpha_fac     ----'
         write(nfout,9008) alpha_fac
         write(nfout,*) ' ---- a(1:',nsum,') ----'
         write(nfout,9008) (g(it),it=1,nsum)
      end if
 9008 format(8f20.12)
    end subroutine gdiis_g

    subroutine opt_forc(nsum)
      integer, intent(in) :: nsum
      integer ::             ia, it, itcrspd, imd_t, j

      forc_g = 0.d0
      do ia = 1, natm
         do j=1,3
         imd_t = imdtypxyz(ia,j)
         if(imd_t == 0 .or. imd_t == COG_FIX .or. &
              & imd_t == COG_FIX_L .or. imd_t == FIX_IN_A_PLANE) cycle
         if(imd_t == RIGID_BODY_FIX_IN_A_PLANE) cycle
         do it = 1, nsum
            itcrspd = ncrspd(it)
            forc_g(ia,j) = forc_g(ia,j) + g(it)*w_l(ia,j,itcrspd)
         end do
         enddo
      end do
    end subroutine opt_forc

    subroutine opt_geom(nsum)
      ! Furthre modification: T. Yamasaki, March/15/2007
      integer, intent(in) :: nsum
      integer ::             ia, ifc, itcrspd, iaa, imd_t,j, num_planes, ifc1, iplane, ia_cog, ia_rb
      real(kind=DP) :: pdot
      real(kind=DP),allocatable,dimension(:,:) :: cps_cog !d(num_planes_atoms_are_fixed,3)
      real(kind=DP),allocatable,dimension(:,:) :: cps_cog_n !d(num_planes_atoms_are_fixed,3)
      real(kind=DP),allocatable,dimension(:)   :: denom ! d(num_planes_atoms_are_fixed)
      integer, allocatable, dimension(:) ::       itpcd ! d(num_planes_atoms_are_fixed)
      integer, allocatable ,dimension(:) ::       sw_cog_planes ! d(num_planes_atoms_are_fixed)
      integer, allocatable ,dimension(:) ::       sw_rb_planes ! d(num_planes_atoms_are_fixed)

      do ia = 1, natm
         do j=1,3
         imd_t = imdtypxyz(ia,j)
         if(imd_t == 0 .or. imd_t == COG_FIX .or. &
              & imd_t == COG_FIX_L .or. imd_t == FIX_IN_A_PLANE) cycle
         if(imd_t == RIGID_BODY_FIX_IN_A_PLANE) cycle
         cps(ia,j) = 0.d0
         do ifc = 1, nsum
            itcrspd = ncrspd(ifc)
            cps(ia,j) = cps(ia,j) + g(ifc)*u_l(ia,j,itcrspd)
         end do
         enddo
      end do

      goto 1001

      select case(constraint_type)
      case (COG_FIX, COG_FIX_L, RIGID_BODY_FIX_IN_A_PLANE, COG_and_RIGID_BODY_FIX_L)
         allocate(itpcd(num_planes_atoms_are_fixed))
         allocate(denom(num_planes_atoms_are_fixed)); denom = 0.d0
         allocate(cps_cog(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
         allocate(cps_cog_n(num_planes_atoms_are_fixed,3)); cps_cog_n = 0.d0

         if(nfcatm_cog >=1 .and. num_planes_atoms_are_fixed >=1 ) then
            allocate(sw_cog_planes(num_planes_atoms_are_fixed)); sw_cog_planes = 0

!!$            call count_independent_cog_planes(num_planes) ! --> num_planes
            do ia_cog = 1, nfcatm_cog
               ia = ia_cnst(ia_cog)
               iplane = ipfixedplane(ia_cog) ! ipfixedplane : fixed_plane's id
               sw_cog_planes(iplane) = 1
               itpcd(iplane) = ia_cog
               cps_cog(iplane,1:3) = cps_cog(iplane,1:3) + amion(ityp(ia))*(cps(ia,1:3)-u_l(ia,1:3,it))
               denom(iplane) = denom(iplane) + amion(ityp(ia))
            end do
            do iplane = 1, num_planes_atoms_are_fixed
               if(sw_cog_planes(iplane)==0) cycle

               cps_cog(iplane,1:3) = cps_cog(iplane,1:3)/denom(iplane)
               iaa = itpcd(iplane)
               pdot = dot_product(fcvect(iaa,1:3),cps_cog(iplane,1:3))
               cps_cog_n(iplane,1:3) = pdot*fcvect(iaa,1:3)   ! perpendicular components
!!$               pdot = dot_product(fcvect(1,1:3),cps_cog(iplane,1:3))
!!$               cps_cog_n(iplane,1:3) = pdot*fcvect(1,1:3)   ! perpendicular components
               if(iprigdiis >= 1) then
                  write(nfout,'(" !!f denom = ",f12.4)') denom(iplane)
                  write(nfout,'(" !!f delta cps_cog(    ",i3,") = ",3f8.4)') iplane,cps_cog(iplane,1:3)
                  write(nfout,'(" !!f delta cps_cog_n(  ",i3,") = ",3f8.4)') iplane,cps_cog_n(iplane,1:3)
               end if
            end do
            do ia_cog = 1, nfcatm_cog
               ia = ia_cnst(ia_cog)
               itcrspd = ipfixedplane(ia_cog)
               cps(ia,1:3) = cps(ia,1:3) - cps_cog_n(itcrspd,1:3)
               if(iprigdiis >= 1) then
                  write(nfout,'(" !!f cps(        ",i3,") = ",3f8.4)') ia, cps(ia,1:3)
                  write(nfout,'(" !!f cps-cps_old(",i3,") = ",3f8.4)') ia, cps(ia,1:3)-u_l(ia,1:3,it)
               end if
            end do

            deallocate(sw_cog_planes)
         end if

         if(nfcatm_rb >= 1 .and. num_planes_atoms_are_fixed >= 1) then
            allocate(sw_rb_planes(num_planes_atoms_are_fixed)); sw_rb_planes = 0
!!$            call count_independent_rb_planes(num_planes)
            cps_cog = 0.d0; cps_cog_n = 0.d0; denom = 0.d0
            do ia_rb = 1, nfcatm_rb
               ia = ia_cnst(nfcatm_cog+ia_rb)
               iplane = ipfixedplane(nfcatm_cog+ia_rb)
               sw_rb_planes(iplane) = 1
               itpcd(iplane) = ia_rb
               cps_cog(iplane,1:3) = cps_cog(iplane,1:3) + amion(ityp(ia))*(cps(ia,1:3)-u_l(ia,1:3,it))
               denom(iplane) = denom(iplane) + amion(ityp(ia))
            end do
            do iplane = 1, num_planes_atoms_are_fixed
               if(sw_rb_planes(iplane)==0) cycle
               cps_cog(iplane,1:3) = cps_cog(iplane,1:3)/denom(iplane)
               iaa = itpcd(iplane)
               pdot = dot_product(fcvect(iaa,1:3),cps_cog(iplane,1:3))
            end do
            do ia_rb = 1, nfcatm_rb
               ia = ia_cnst(nfcatm_cog+ia_rb)
               itcrspd = ipfixedplane(nfcatm_cog+ia_rb)
               cps(ia,1:3) = cps(ia,1:3) - cps_cog_n(itcrspd,1:3)
               if(iprigdiis >= 1) then
                  write(nfout,'(" !!f cps(        ",i3,") = ",3f8.4," (rigid_body)")') ia, cps(ia,1:3)
                  write(nfout,'(" !!f cps-cps_old(",i3,") = ",3f8.4," (rigid_body)")') ia, cps(ia,1:3)-u_l(ia,1:3,it)
               end if
            end do
         end if

         deallocate(cps_cog_n)
         deallocate(cps_cog)
         deallocate(denom)
         deallocate(itpcd)

      end select

!!$      if(constraint_type == COG_FIX_L &
!!$           & .or. constraint_type == COG_FIX .or. constraint_type == COG_CNTR) then
!!$         allocate(itpcd(num_planes_atoms_are_fixed))
!!$         allocate(cps_cog(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
!!$         allocate(cps_cog_n(num_planes_atoms_are_fixed,3)); cps_cog_n = 0.d0
!!$         allocate(denom(num_planes_atoms_are_fixed)); denom = 0.d0
!!$
!!$         do ifc = 1, nfcatm
!!$            ia = ia_cnst(ifc)
!!$            iaa = ipfixedplane(ifc)
!!$            itpcd(iaa) = ifc
!!$            cps_cog(iaa,1:3) = cps_cog(iaa,1:3) + amion(ityp(ia))*(cps(ia,1:3)-u_l(ia,1:3,it))
!!$            denom(iaa) = denom(iaa) + amion(ityp(ia))
!!$         end do
!!$         do ifc = 1, num_planes_atoms_are_fixed
!!$            cps_cog(ifc,1:3) = cps_cog(ifc,1:3)/denom(ifc)
!!$            iaa = itpcd(ifc)
!!$!$$!!$            pdot = dot_product(fcvect(iaa,1:3),cps_cog(ifc,1:3))
!!$!$$!!$            cps_cog_n(ifc,1:3) = pdot*fcvect(iaa,1:3)   ! perpendicular components
!!$            pdot = dot_product(fcvect(1,1:3),cps_cog(ifc,1:3))
!!$            cps_cog_n(ifc,1:3) = pdot*fcvect(1,1:3)   ! perpendicular components
!!$            if(iprigdiis >= 1) then
!!$               write(nfout,'(" !!f denom = ",f12.4)') denom(ifc)
!!$               write(nfout,'(" !!f delta cps_cog(    ",i3,") = ",3f8.4)') ifc,cps_cog(ifc,1:3)
!!$               write(nfout,'(" !!f delta cps_cog_n(  ",i3,") = ",3f8.4)') ifc,cps_cog_n(ifc,1:3)
!!$            end if
!!$         end do
!!$         do ifc = 1, nfcatm
!!$            ia = ia_cnst(ifc)
!!$            itcrspd = ipfixedplane(ifc)
!!$! -------------->
!!$            cps(ia,1:3) = cps(ia,1:3) - cps_cog_n(itcrspd,1:3)
!!$! <-------------
!!$            if(iprigdiis >= 1) then
!!$               write(nfout,'(" !!f cps(        ",i3,") = ",3f8.4)') ia, cps(ia,1:3)
!!$               write(nfout,'(" !!f cps-cps_old(",i3,") = ",3f8.4)') ia, cps(ia,1:3)-u_l(ia,1:3,it)
!!$            end if
!!$         end do
!!$
!!$         if(iprigdiis >= 1) then
!!$            cps_cog = 0.d0
!!$            cps_cog_n = 0.d0
!!$            do ifc = 1, nfcatm
!!$               ia = ia_cnst(ifc)
!!$               iaa = ipfixedplane(ifc)
!!$               itpcd(iaa) = ifc
!!$               cps_cog(iaa,1:3) = cps_cog(iaa,1:3) + amion(ityp(ia))*cps(ia,1:3)
!!$               cps_cog_n(iaa,1:3) = cps_cog_n(iaa,1:3) + amion(ityp(ia))*u_l(ia,1:3,it)
!!$            end do
!!$            do ifc = 1, num_planes_atoms_are_fixed
!!$               cps_cog(ifc,1:3) = cps_cog(ifc,1:3)/denom(ifc)
!!$               cps_cog_n(ifc,1:3) = cps_cog_n(ifc,1:3)/denom(ifc)
!!$               write(nfout,'(" !!f real cps_old_cog(    ",i3,") = ",3f8.4)') ifc,cps_cog_n(ifc,1:3)
!!$               write(nfout,'(" !!f real cps_cog(        ",i3,") = ",3f8.4)') ifc,cps_cog(ifc,1:3)
!!$            end do
!!$         end if

!!$         deallocate(denom)
!!$         deallocate(cps_cog_n)
!!$         deallocate(cps_cog)
!!$         deallocate(itpcd)
!!$      end if

1001  continue
    end subroutine opt_geom

    subroutine count_independent_cog_planes(num_planes)
      integer, intent(out) :: num_planes

      integer :: ifc, ifc1
      integer, allocatable, dimension(:) :: i_planes
      allocate(i_planes(num_planes_atoms_are_fixed)); i_planes = 0

      if(nfcatm_cog==1) then
         num_planes = 1
      else if(nfcatm_cog >= 2) then
         i_planes(1:nfcatm_cog) = ipfixedplane(1:nfcatm_cog)
         do ifc = 2, nfcatm_cog
            do ifc1 = 1, ifc-1
               if(i_planes(ifc) == i_planes(ifc1)) then
                  i_planes(ifc) = 0
                  exit
               end if
            end do
         end do
         num_planes = 0
         do ifc = 1, nfcatm_cog
            if(i_planes(ifc) /= 0) num_planes = num_planes+1
         end do
      else
         num_planes = 0
      end if
      deallocate(i_planes)
    end subroutine count_independent_cog_planes

    subroutine count_independent_rb_planes(num_planes)
      integer, intent(out) :: num_planes
      integer :: ifc, ifc1
      integer, allocatable, dimension(:) :: i_planes
      allocate(i_planes(num_planes_atoms_are_fixed)); i_planes = 0

      if(nfcatm_cog==1) then
         num_planes = 1
      else if(nfcatm_cog >= 2) then
         i_planes(1:nfcatm_cog) = ipfixedplane(1:nfcatm_cog)
         do ifc = 2, nfcatm_cog
            do ifc1 = 1, ifc-1
               if(i_planes(ifc) == i_planes(ifc1)) then
                  i_planes(ifc) = 0
                  exit
               end if
            end do
         end do
         num_planes = 0
         do ifc = 1, nfcatm_cog
            if(i_planes(ifc) /= 0) num_planes = num_planes+1
         end do
      else
         num_planes = 0
      end if
      deallocate(i_planes)
    end subroutine count_independent_rb_planes

    subroutine check_cog()
      integer :: ifc, ia, iaa
      real(kind=DP),allocatable,dimension(:,:) :: cps_cog !d(num_planes_atoms_are_fixed,3)
      real(kind=DP),allocatable,dimension(:,:) :: cps_cog_n !d(num_planes_atoms_are_fixed,3)
      real(kind=DP),allocatable,dimension(:)   :: denom ! d(num_planes_atoms_are_fixed)
      integer, allocatable, dimension(:) ::       itpcd ! d(num_planes_atoms_are_fixed)

      if(iprigdiis >= 1) then
         allocate(itpcd(num_planes_atoms_are_fixed))
         allocate(cps_cog(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
         allocate(cps_cog_n(num_planes_atoms_are_fixed,3)); cps_cog_n = 0.d0
         allocate(denom(num_planes_atoms_are_fixed)); denom = 0.d0

         do ifc = 1, nfcatm
            ia = ia_cnst(ifc)
            iaa = ipfixedplane(ifc)
            itpcd(iaa) = ifc
         end do

         cps_cog = 0.d0
         cps_cog_n = 0.d0
         do ifc = 1, nfcatm
            ia = ia_cnst(ifc)
            iaa = ipfixedplane(ifc)
            itpcd(iaa) = ifc
            cps_cog(iaa,1:3) = cps_cog(iaa,1:3) + amion(ityp(ia))*cps(ia,1:3)
            cps_cog_n(iaa,1:3) = cps_cog_n(iaa,1:3) + amion(ityp(ia))*u_l(ia,1:3,it)
            denom(iaa) = denom(iaa) + amion(ityp(ia))
         end do
         do ifc = 1, num_planes_atoms_are_fixed
            cps_cog(ifc,1:3) = cps_cog(ifc,1:3)/denom(ifc)
            cps_cog_n(ifc,1:3) = cps_cog_n(ifc,1:3)/denom(ifc)
            write(nfout,'(" !!f real cps_old_cog(    ",i3,") = ",3f8.4)') ifc,cps_cog_n(ifc,1:3)
            write(nfout,'(" !!f real cps_cog(        ",i3,") = ",3f8.4)') ifc,cps_cog(ifc,1:3)
         end do
         deallocate(denom)
         deallocate(cps_cog_n)
         deallocate(cps_cog)
         deallocate(itpcd)
      end if

    end subroutine check_cog

    subroutine cmp_fop_fin(it)
      integer, intent(in) :: it
      real(kind=DP) :: xmul, xmul_input_f
      integer :: ia,j
      xmul = 0.d0
      xmul_input_f = 0.d0
      do ia = 1, natm
         do j=1,3
         if(imdtypxyz(ia,j) /= 0) then
            xmul = xmul + forc_g(ia,j)**2
            xmul_input_f = xmul_input_f + w_l(ia,j,it)**2
         end if
         enddo
      end do
      if(printable) write(nfout,*) ' norm of optimal force = ', xmul
      if(printable) write(nfout,*) ' norm of input   force = ', xmul_input_f
    end subroutine cmp_fop_fin

    subroutine cps_damp(c_dist,it,icount)
      real(kind=DP), intent(in) :: c_dist
      integer, intent(in) ::      it
      integer ,intent(out) ::     icount
      integer ::   ia,j, iaa, ifc, itcrspd
      real(kind=DP) :: df
!!$      real(kind=DP) :: pdot
      integer, allocatable, dimension(:) :: sw_damp !d(natm)
!!$      real(kind=DP),allocatable,dimension(:,:) :: cps_cog !d(num_planes_atoms_are_fixed,3)
!!$      real(kind=DP),allocatable,dimension(:,:) :: cps_cog_n !d(num_planes_atoms_are_fixed,3)
!!$      real(kind=DP),allocatable,dimension(:)   :: denom ! d(num_planes_atoms_are_fixed)
!!$      integer, allocatable, dimension(:) ::       itpcd ! d(num_planes_atoms_are_fixed)

      allocate(sw_damp(natm)); sw_damp = YES
      if(constraint_type == COG_FIX_L &
           & .or. constraint_type == COG_FIX .or. constraint_type == COG_CNTR) then
!!$         allocate(itpcd(num_planes_atoms_are_fixed))
!!$         allocate(cps_cog(num_planes_atoms_are_fixed,3)); cps_cog = 0.d0
!!$         allocate(cps_cog_n(num_planes_atoms_are_fixed,3)); cps_cog_n = 0.d0
!!$         allocate(denom(num_planes_atoms_are_fixed)); denom = 0.d0
         do j = 1, nfcatm
            ia = ia_cnst(j)
            sw_damp(ia) = NO
!!$            iaa = ipfixedplane(j)
!!$            itpcd(iaa) = j
!!$            cps_cog(iaa,1:3) = cps_cog(iaa,1:3) + amion(ityp(ia))*(cps(ia,1:3)-u_l(ia,1:3,it))
!!$            denom(iaa) = denom(iaa) + amion(ityp(ia))
         end do
!!$         do ifc = 1, num_planes_atoms_are_fixed
!!$            cps_cog(ifc,1:3) = cps_cog(ifc,1:3)/denom(ifc)
!!$            iaa = itpcd(ifc)
!!$            pdot = dot_product(fcvect(1,1:3),cps_cog(ifc,1:3))
!!$            cps_cog_n(ifc,1:3) = pdot*fcvect(1,1:3)   ! perpendicular components
!!$         end do
!!$         do ifc = 1, nfcatm
!!$            ia = ia_cnst(ifc)
!!$            itcrspd = ipfixedplane(ifc)
!!$            cps(ia,1:3) = cps(ia,1:3) - cps_cog_n(itcrspd,1:3)
!!$         end do
!!$         deallocate(denom)
!!$         deallocate(cps_cog_n)
!!$         deallocate(cps_cog)
!!$         deallocate(itpcd)
      end if

      icount = 0
      do ia = 1, natm
         if(sw_damp(ia) == NO) cycle
         do j = 1, 3
            df = cps(ia,j) - u_l(ia,j,it)
            if(df > c_dist) then
               icount = icount + 1
               cps(ia,j) = u_l(ia,j,it) + c_dist
            else if(df < -c_dist) then
               icount = icount + 1
               cps(ia,j) = u_l(ia,j,it) - c_dist
            end if
         end do
      end do
      deallocate(sw_damp)
    end subroutine cps_damp

  end subroutine m_IS_gdiis

  subroutine  mdfy_forc(forcmx,forcmx_mdfy)
    real(kind=DP), intent(in) ::  forcmx
    real(kind=DP), intent(out) :: forcmx_mdfy

    integer :: iaa, ic_cog, ifcatm, ia
    integer, allocatable, dimension(:) :: npfatm  ! d(natm)
    real(kind=DP), dimension(3) :: fcg, fcg_mdfy
    real(kind=DP) :: pdot, fca


    if(.not.constraints_exist) then
       forcmx_mdfy = forcmx
    else

       allocate(npfatm(natm)); npfatm = 0
       ifcatm = 0
       do ia = 1, natm
          if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L &
               .or. imdtyp(ia) == FIX_IN_A_PLANE) then
             ifcatm = ifcatm + 1
             npfatm(ia) = ifcatm
             if(imdtyp(ia) == COG_FIX .or. imdtyp(ia) == COG_FIX_L) iaa = ia
          end if
       end do

       if(constraint_type == COG_CNTR .or. constraint_type == COG_FIX_L ) then
          ic_cog = 0
          do ia = 1, natm
             if(imdtyp(ia) == COG_CNTR .or. imdtyp(ia) == COG_FIX_L) then
                ic_cog = ic_cog + 1
                fcg(1) = fcg(1) + forc_g(ia,1)
                fcg(2) = fcg(2) + forc_g(ia,2)
                fcg(3) = fcg(3) + forc_g(ia,3)
             end if
             ! fcg(1:3) = \sum_{ia} forc_g(ia,1:3)
          end do
          if(ic_cog >= 2) then
             fcg = fcg/ic_cog
             ! fcg(1:3) = \sum_{ia} forc_g(ia,1:3)/(\sum_{ia}) := f
             ifcatm = npfatm(iaa)
             pdot = dot_product(fcvect(ifcatm,1:3),fcg)
             fcg_mdfy = fcg - pdot*fcvect(ifcatm,1:3)   ! parallel components
             ! = f_{parallel}, pdot*fcvect(ifcatm,:) = f_{perpendicular}
          else if(ic_cog == 1) then
             call phase_error_with_msg(nfout,' #ic_cog is not enough <<m_IS_gdiis.mdfy_forc>>',__LINE__,__FILE__)
          else
             fcg_mdfy = 0.d0
          end if
       end if

       do ia = 1, natm
          if(imdtyp(ia) == FIX_IN_A_PLANE .or. imdtyp(ia) == COG_CNTR &
               & .or. imdtyp(ia) == COG_FIX_L) then
             ifcatm = npfatm(ia)
             if(imdtyp(ia) == COG_CNTR .or. imdtyp(ia) == COG_FIX_L) then
!!$                forc_g(ia,1:3) = fcg_mdfy + forc_g(ia,1:3) - fcg
                forc_g(ia,1:3) = forc_g(ia,1:3) - (fcg - fcg_mdfy)
                !  f_{ia} = f_{ia} - f_{perpendicular}
             else
                pdot = dot_product(fcvect(ifcatm,1:3),forc_g(ia,1:3))
                forc_g(ia,1:3) = forc_g(ia,1:3) - pdot*fcvect(ifcatm,1:3)
             end if
           end if
       end do

       forcmx_mdfy = 0.d0
       do ia = 1, natm
          if(imdtyp(ia) /= 0) then
             fca = dsqrt(forc_g(ia,1)**2 + forc_g(ia,2)**2 + forc_g(ia,3)**2)
             if(forcmx_mdfy < fca) forcmx_mdfy = fca
          end if
       end do

       deallocate(npfatm)
    end if
  end subroutine mdfy_forc

  subroutine forc_check(name)
    character(len=*),intent(in) :: name
    integer :: ia,j
    write(nfout,'(a40)') name
    do ia = 1, natm
       write(nfout,'(i4,3f20.12)') ia, (forc_g(ia,j),j=1,3)
    end do
  end subroutine forc_check

  subroutine cps_check(name)
    character(len=*),intent(in) :: name
    integer :: ia,j
    write(nfout,'(a40)') name
    do ia = 1, natm
       write(nfout,'(i4,3f20.12)') ia, (cps(ia,j),j=1,3)
    end do
  end subroutine cps_check

  subroutine cpd_check(name)
    character(len=*),intent(in) :: name
    integer :: ia,j
    write(nfout,'(a40)') name
    do ia = 1, natm
       write(nfout,'(i4,3f20.12)') ia, (cpd_l(ia,j),j=1,3)
    end do
  end subroutine cpd_check

  subroutine m_IS_set_mdmode_cnstra()
    call m_IS_gdiis_reset()
    call m_IS_CG_reset()
    call m_IS_reset_extrpl_status()
    call m_CtrlP_reset_optmode()
    iteration_ionic_at_CNSTRA = iteration_ionic
    mdmode = CNSTRA
  end subroutine m_IS_set_mdmode_cnstra

  subroutine m_IS_set_mdmode_ordina()
    mdmode = ORDINA
  end subroutine m_IS_set_mdmode_ordina

  subroutine m_IS_phonon_force()
    integer :: id,ic,ic0,iat
    id = iteration_ionic + istart_phonon - 1
    ic0 = iconf(id)
    displaced_atom = phonon_atom(ic0)
    do while(.true.)
       if(iteration_ionic>=iend_phonon-istart_phonon+1) then
         return
       endif
       ic = iconf(iteration_ionic+istart_phonon)
       iat = phonon_atom(ic)
       if(imdtyp(iat) /= OFF .and. .not. force_was_read(ic)) exit
       call m_Iter_mdIterN_increment(nfout)
    enddo
    cps(displaced_atom,1:3) = cps(displaced_atom,1:3) - phonon_displacement(ic0,1:3)
    id = iteration_ionic + istart_phonon
    ic = iconf(id)
    displaced_atom = phonon_atom(ic)
    displacement(1:3) = phonon_displacement(ic,1:3)
    cps(displaced_atom,1:3) = cps(displaced_atom,1:3) + displacement(1:3)
    phonon_iteration = phonon_iteration+1
  end subroutine m_IS_phonon_force

  subroutine m_IS_phonon_set_displacement()  ! asms
    integer :: id,ic
    id = iteration_ionic + istart_phonon - 1
    ic = iconf(id)
    displaced_atom    = phonon_atom(ic)
    displacement(1:3) = phonon_displacement(ic,1:3)
  end subroutine m_IS_phonon_set_displacement

  subroutine m_IS_phonon_equilibrium()
    integer :: id,ic
    id = iteration_ionic + istart_phonon - 1
    ic = iconf(id)
    displaced_atom = phonon_atom(ic)
    cps(displaced_atom,1:3) = cps(displaced_atom,1:3) - phonon_displacement(ic,1:3)
  end subroutine m_IS_phonon_equilibrium

  subroutine m_IS_phonon_init_firsthalf()
    use m_Files, only : nfout
    use m_Const_Parameters, only : PUCV,PAI2,PHONON_FORCE
    use m_Crystal_Structure, only : rltv,b2pmat

    real(kind=DP) :: work(3)
    integer :: i,ic

    if(.not.allocated(force_was_read)) then
      allocate(force_was_read(3*natm_prim*norder*2));force_was_read = .false.
    endif

    if(sw_calc_force == ON) then
      imdalg = PHONON_FORCE
      if(printable) write(nfout,*) 'PHONON: imdalg = ',imdalg
    end if

    call set_phonon_displacemets()
    call search_equiv_config()
    if(istart_phonon <= 0 .or. istart_phonon > num_force_calc) &
      & istart_phonon=1
    if(iend_phonon <= 0 .or. iend_phonon > num_force_calc) &
      & iend_phonon=num_force_calc
    if(iend_phonon<istart_phonon) iend_phonon=istart_phonon
    if(printable) write(nfout,*) &
      & 'PHONON: istart, iend = ',istart_phonon,iend_phonon
  contains
    subroutine set_phonon_displacemets()
      ! local variables
      integer :: ia,is,id

      num_force_data = natm_prim*3*norder*2
      if(sw_polynomial_fit == ON) num_force_data = num_force_data + 1
      allocate(phonon_atom(num_force_data))
      allocate(phonon_displacement(num_force_data,3))

      if(sw_polynomial_fit == ON) then
         num_force_data = 1
         phonon_atom(num_force_data) = 1
         phonon_displacement(num_force_data,1:3) = 0.d0
      else
         num_force_data = 0
      end if
      do ia=1,natm_prim
         do id=1,3
            do is=-norder,norder
               if(is == 0) cycle
               num_force_data = num_force_data + 1
               phonon_atom(num_force_data) = ia
               phonon_displacement(num_force_data,1:3) = 0.d0
               phonon_displacement(num_force_data,id) =  u/dble(is)
            end do
         end do
      end do

    end subroutine set_phonon_displacemets

    subroutine search_equiv_config()
      ! local variables
      integer :: i,j,ia,ja,iopr,ia1,ia2,k
      real(kind=DP) :: ui(3),uj(3),t(3),ruj(3),p1(3),p2(3),f(3)
!!$      real(kind=DP), parameter :: ddd = 1.d-12

      allocate(napt_phonon(natm_super,num_force_data)); napt_phonon=0
      allocate(iequconf(num_force_data)); iequconf = 0
      allocate(iopr_equconf(num_force_data)); iopr_equconf = 0

      if(sw_calc_force_all == OFF) then
         do i=2,num_force_data
            ia=phonon_atom(i)
            ui=phonon_displacement(i,1:3)
            Loop_conf: do j=1,i-1
               ja=phonon_atom(j)
               uj=phonon_displacement(j,1:3)
               Loop_op: do iopr=1,nopr
                  if(ia/=napt(ja,iopr)) cycle
                  do k=1,3
                     ruj(k) = dot_product(op(k,1:3,iopr),uj)
                  end do
                  if(sum(abs(ui(1:3)-ruj(1:3)))>1.d-8) cycle
                  do k=1,3
                     t(k) = cps(ia,k) - dot_product(op(k,1:3,iopr),cps(ja,1:3))
                  end do
                  Loop_atom: do ia1=1,natm_super
                     do k=1,3
                        p1(k) = dot_product(op(k,1:3,iopr),cps(ia1,1:3)) + t(k)
                     end do
                     do ia2=1,natm_super
                        if(ityp(ia1) /= ityp(ia2)) cycle
                        p2(1:3) = cps(ia2,1:3)-p1(1:3)
                        f(1) = abs(cos(sum(rltv(1:3,1)*p2))-1.d0)
                        f(2) = abs(cos(sum(rltv(1:3,2)*p2))-1.d0)
                        f(3) = abs(cos(sum(rltv(1:3,3)*p2))-1.d0)
!!$                        if(maxval(f) <= ddd) then
                        if(maxval(f) <= symmetry_check_criterion) then
                           napt_phonon(ia1,i) = ia2
                           cycle Loop_atom
                        end if
                     end do
                     cycle Loop_op
                  end do Loop_atom
                  if(ia==ja .and. sum(abs(ui(1:3)+uj(1:3)))<1.d-8) then
                     iequconf(i) = -j
                  else
                     iequconf(i) = j
                  end if
                  iopr_equconf(i) = iopr
                  exit Loop_conf
               end do Loop_op
            end do Loop_conf
         end do
      end if

      num_force_calc=0
      num_force_calc_mobile = 0
      do i=1,num_force_data
         if(iequconf(i)<=0)then
           num_force_calc=num_force_calc+1
           if(imdtyp(phonon_atom(i))/=OFF) num_force_calc_mobile = num_force_calc_mobile+1
         endif
      end do
      allocate(iconf(num_force_calc))
      j=0
      do i=1,num_force_data
         if(iequconf(i)<=0) then
            j = j+1
            iconf(j)=i
         end if
      end do

      if(ipriphonon>=1) then
         write(nfout,'(" PHONON: equivalent configurations")')
         write(nfout,'(" PHONON: i, iequconf, iopr_equconf")')
         do i=1,num_force_data
            write(nfout,'(" PHONON: ",i5,1x,i5,1x,i5)') i,iequconf(i),iopr_equconf(i)
         end do
         write(nfout,'(" PHONON: number of force calculations = ",i5)') num_force_calc
         write(nfout,'(" PHONON: i, iconf")')
         do i=1,num_force_calc
            write(nfout,'(" PHONON: ",i5,1x,i5)') i,iconf(i)
         end do
         do i=1,num_force_data
            write(nfout,'(" PHONON: napt_phonon i=",i5)') i
            write(nfout,'(" PHONON: ",10i5)') napt_phonon(1:natm_super,i)
         end do
      end if
    end subroutine search_equiv_config

  end subroutine m_IS_phonon_init_firsthalf

  subroutine m_IS_phonon_init_secondhalf()
    use m_Files, only : nfout
    use m_Const_Parameters, only : PUCV,PAI2,PHONON_FORCE
    use m_Crystal_Structure, only : rltv,b2pmat

    real(kind=DP) :: work(3)
    integer :: i,ic

    do while(.true.)
      if(istart_phonon+iteration_ionic-1>num_force_calc) exit
      ic = iconf(istart_phonon+iteration_ionic-1)
      displaced_atom = phonon_atom(ic)
      displacement(1:3) = phonon_displacement(ic,1:3)
      phonon_iteration = 1
      if(imdtyp(phonon_atom(ic)) /= OFF .and. .not. force_was_read(ic)) exit
      call m_Iter_mdIterN_increment(nfout)
    enddo

  end subroutine m_IS_phonon_init_secondhalf

  subroutine m_IS_phonon_initial_disp()
    use m_Crystal_Structure, only : m_CS_phonon_symmetry
    integer :: i, id, ic

! =====ASMS_DEBUG ====================== 2013/02/07
    if ( sw_calc_force == OFF ) return
! =====ASMS_DEBUG ====================== 2013/02/07

    if(sw_phonon_oneshot == ON) then
       id = num_phonon_calc_mode
       if(id > num_force_calc) id = num_force_calc
       ic = iconf(id)
       if(ipri>=1) write(nfout,'(" ic, id = ", 2i8)') ic, id
       displaced_atom = phonon_atom(ic)
       displacement(1:3) = phonon_displacement(ic,1:3)
       cps(displaced_atom,1:3) = cps(displaced_atom,1:3) + displacement(1:3)
    else
       cps(displaced_atom,1:3) = cps(displaced_atom,1:3) + displacement(1:3)
    end if
    call m_IS_cps_to_pos()
    if(sw_calc_force == ON) then
      call m_CS_phonon_symmetry(OFF)
    end if
  end subroutine m_IS_phonon_initial_disp

  subroutine m_IS_set_ionic_mass(nfout)
    integer, intent(in) :: nfout

    integer :: ia, ic, is, ie

    do ia = 1,natm
       if(ionic_mass(ia)<0.d0) ionic_mass(ia) = amion(ityp(ia))
    end do

    if(natom_reservoir>0)then
       do ia=1,natom_reservoir
          if(atom_reservoir(ia)%ionic_mass<0.d0) &
        & atom_reservoir(ia)%ionic_mass = amion(atom_reservoir(ia)%ityp)
       enddo
    endif
    if(printable) then
       if(natm <= 3) then
          write(nfout,'(" !**   ia, ionic_mass ")')
          do ia=1,natm
             write(nfout,'(" !** ",i4,1x,f20.5)') ia,ionic_mass(ia)
          end do
       else
          write(nfout,'(" !**   ia1 - ia2, ityp, ionic_mass ")')
          ic = ityp(1)
          is = 1; ie = 0
          do ia=2,natm
             if(ityp(ia) /= ic) then
                ie = ia-1
                write(nfout,'(" !** ",i4," - ",i4,i4,f20.5)') is,ie,ityp(ie),ionic_mass(ie)
                is = ia
                ic = ityp(ia)
             end if
          end do
          ie = ia-1
          write(nfout,'(" !** ",i4," - ",i4,i4,f20.5)') is,ie,ityp(ie),ionic_mass(ie)
       end if
    end if
  end subroutine m_IS_set_ionic_mass

  subroutine m_IS_pack_all_ions_in_uc(ityp_full_r,cps_full)
    integer, intent(out), dimension(natm2) ::       ityp_full_r
    real(kind=DP),intent(out),dimension(natm2,3) :: cps_full
    ! All packed ion coordinates in the primitive unit cell system
    real(kind=DP),allocatable,dimension(:,:) :: pos_wk
    real(kind=DP),allocatable,dimension(:)   :: r_wk
    integer :: i, m, ucret

    !!!allocate(r_wk(3))

    allocate(pos_wk(natm2,3))
    pos_wk(1:natm,:) = pos(1:natm,:)

    ityp_full_r(1:natm) = ityp(1:natm)
    if(natm2 > natm) then
       call rplcps(pos_wk,ityp_full_r,1,natm2,natm,iwei)  ! -(b_Ionic_System)
       call rbinuc(pos_wk,natm2)  ! -(b_Ionic_System)
       call change_of_coordinate_system(altv,pos_wk,natm2,natm2,cps_full) ! -(b_Ionic_System)
       !! for unit conversion to angstrom
       !!do i = 1, natm2
       !!   !!!! K.Mae 040315
       !!   do m = 1, 3
       !!      ucret = unit_conv_byname( cps_full(i,m), r_wk(m), 'bohr', 'angstrom' )
       !!   end do
       !!   cps_full(i,1:3) = r_wk(1:3)
       !!end do
    else
       call rbinuc(pos_wk,natm)
       call change_of_coordinate_system(altv,pos_wk,natm,natm,cps_full) ! -(b_C.S.)
       !! for unit conversion to angstrom
       !!do i = 1, natm
       !!   !!! K.Mae 040315
       !!   do m = 1, 3
       !!      ucret = unit_conv_byname( cps_full(i,m), r_wk(m), 'bohr', 'angstrom' )
       !!   end do
       !!   cps_full(i,1:3) = r_wk(1:3)
       !!end do
    end if
    deallocate(pos_wk)
    !!!deallocate(r_wk)
  end subroutine m_IS_pack_all_ions_in_uc

  subroutine m_IS_CG_reset()
     iter_CG = 1
     iter_linmin = 1
     iter_CG_max = 1
  end subroutine m_IS_CG_reset

  subroutine m_IS_fire_param_by_args(nfout,ndim,iteration,cps,cpd,forc_l_in,imdtyp, &
   &                                 iteration_from_lq, fire_alpha, fire_dt, forcold)
    integer, intent(in) :: nfout,ndim
    integer, intent(in) :: iteration
    real(kind=DP), intent(inout), dimension(ndim,3) :: cps,cpd
    real(kind=DP), intent(in), dimension(ndim,3) :: forc_l_in
    integer, intent(in), dimension(ndim,3) :: imdtyp
    integer, dimension(ndim) :: iteration_from_lq
    real(kind=DP),  intent(inout), dimension(ndim)  :: fire_alpha
    real(kind=DP),  intent(inout), dimension(ndim)  :: fire_dt
    real(kind=DP),  intent(inout), dimension(ndim,3)  :: forcold
    real(kind=DP), allocatable, dimension(:,:) :: forc
    real(kind=DP), dimension(3) :: r1,unitf,v1,f,v,r,fold
    integer :: i,j,k

    allocate(forc(ndim,3));forc(:,:) = forc_l_in(:,:)
    do i=1,ndim
       do j=1,3
         if (imdtyp(i,j) == FIX) forc(i,j) = 0.d0
       enddo
    enddo
    do i=1,ndim
       r(:) = cps(i,:)
       v(:) = cpd(i,:)
       f(:) = forc(i,:)
       fold(:) = forcold(i,:)
       v1 = v + 0.5d0*fire_dt(i)*(f+fold)*fire_invmass_factor
       v = 0.d0
       unitf=0.d0
       if(dot_product(f,f)/=0) unitf=f/dsqrt(dot_product(f,f))
       if(dot_product(v1,unitf) >= 0.0d0) then
         v = (1.d0-fire_alpha(i))*v1 + fire_alpha(i) * unitf * dot_product(v1,unitf)
         if(iteration-iteration_from_lq(i)>=fire_nmin)then
           fire_dt(i) = min(fire_dt(i)*fire_incre_factor,fire_dtmax)
           fire_alpha(i) = fire_alpha(i)*fire_decre_factor_alpha
           if(iprifire>=2) write(nfout,'(a,i8,2f20.15)') &
           & '!** FIRE new alpha and dt for freedom no. : ',i,fire_alpha(i),fire_dt(i)
         endif
       else
         iteration_from_lq(i) = iteration
         fire_alpha(i) = fire_alpha_start
         fire_dt(i) = fire_dt(i) * fire_decre_factor
         if(iprifire>=2) then
           write(nfout,'(a,i8)') '!** FIRE quenched freedom no. ',i
           write(nfout,'(a,i8,2f20.15)') '!** FIRE new alpha and dt for freedom no. : ' &
          &                              ,i,fire_alpha(i),fire_dt(i)
         endif
       end if
       r1 = r + v*fire_dt(i) + fire_invmass_factor*f*fire_dt(i)**2
       cps(i,:) = r1(:);cpd(i,:) = v(:)
    enddo
    forcold = forc
    deallocate(forc)
  end subroutine m_IS_fire_param_by_args

  subroutine m_IS_fire_core(nfout,ndim,iteration,cps,cpd,forc_l_in,imdtyp)
    integer, intent(in) :: nfout,ndim
    integer, intent(in) :: iteration
    real(kind=DP), intent(inout), dimension(ndim,3) :: cps,cpd
    real(kind=DP), intent(in), dimension(ndim,3) :: forc_l_in
    integer, intent(in), dimension(ndim,3) :: imdtyp
    integer, allocatable, dimension(:), save :: iteration_from_lq
    real(kind=DP), allocatable, dimension(:), save :: fire_alpha
    real(kind=DP), allocatable, dimension(:), save :: fire_dt
    real(kind=DP), allocatable, dimension(:,:), save :: forcold

    if(.not.allocated(fire_alpha))then
       allocate(fire_alpha(ndim));fire_alpha = fire_alpha_start
    endif
    if(.not.allocated(fire_dt)) then
       allocate(fire_dt(ndim));fire_dt = fire_initial_dt
    endif
    if(.not.allocated(iteration_from_lq)) then
       allocate(iteration_from_lq(ndim));iteration_from_lq=0
    endif
    if(.not.allocated(forcold)) then
      allocate(forcold(ndim,3));forcold = forc_l_in
    endif
    call m_IS_fire_param_by_args(nfout,ndim,iteration,cps,cpd,forc_l_in,imdtyp,iteration_from_lq, fire_alpha, fire_dt, forcold)
  end subroutine m_IS_fire_core

!! Added by J.N. 6/Jan/2013
!! Previous version of CG installed by T.Y was deleted.
!! If that verson is needed, check out from svn.
  subroutine m_IS_cg2_core(nfout,ndim,cps,forc_l_in,imdtyp,etotal)
    integer, intent(in) :: nfout,ndim
    real(kind=DP), intent(in), dimension(ndim,3) :: forc_l_in
    real(kind=DP), intent(inout), dimension(ndim,3) :: cps
    integer, intent(in), dimension(ndim,3) :: imdtyp
    real(kind=DP), intent(in)                    :: etotal

    real(kind=DP), save :: alpha, alpha2, gamma
    real(kind=DP), save :: f_para1, f_para2, f_para3
    real(kind=DP), save :: vec_h_norm
    !!$integer, save :: iter_CG = 1, iter_linmin = 1, iter_CG_max = 1
    real(kind=DP), dimension(:,:), allocatable, save :: vec_g, vec_h          ! dim(natm,3)
    real(kind=DP), dimension(:,:), allocatable, save ::       f_total0        ! dim(natm,3)
    real(kind=DP), dimension(:,:), allocatable, save :: cps1, f_total1        ! dim(natm,3)
    real(kind=DP), dimension(:,:), allocatable, save :: cps2, f_total2        ! dim(natm,3)
    integer :: ia,j


    integer ::              id_sname = -1
    call tstatc0_begin('m_IS_cg2',id_sname)

!     if(absolute_convergence_of_forc(forc_l_in))then
!       write(nfout,'(a)') ' CG2: forces are absolutely converged!! nothing to do...'
!       return
!     endif

    if(.not.allocated(vec_g)) allocate(vec_g(ndim,3))
    if(.not.allocated(vec_h)) allocate(vec_h(ndim,3))
    if(.not.allocated(f_total0)) allocate(f_total0(ndim,3))
    if(.not.allocated(cps1)) allocate(cps1(ndim,3))
    if(.not.allocated(f_total1)) allocate(f_total1(ndim,3))
    if(.not.allocated(cps2)) allocate(cps2(ndim,3))
    if(.not.allocated(f_total2)) allocate(f_total2(ndim,3))

    f_total0(:,:) = forc_l_in(:,:)


    if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_hyperplane(f_total0)

     do ia=1,ndim
       do j=1,3
         if (imdtyp(ia,j) == FIX) then
           f_total0(ia,j) = 0.0d0
         endif
       enddo
     end do

     if ( ( iter_CG .eq. 1 ) .and. ( iter_linmin .eq. 1 ) ) then
       call first_cg_step()
     else if ( iter_linmin .eq. 2 ) then
       f_para2 = sum(f_total0(:,:)*vec_h(:,:))/sqrt(sum(vec_h(:,:)**2))
       if(printable) write(nfout,'(a18,2f10.6)') ' CG2: F_parallel =', f_para1, f_para2
       if ( abs(f_para2) .lt. f_para1*0.2d0 ) then
         if(printable) write(nfout,*) 'CG2: F_parallel is small. Then linmin finished.'
         if ( iter_CG .lt. iter_CG_max ) then
           if(printable) write(nfout,*) 'CG2: Go next_cg_step'
           call next_cg_step()
         else
           if(printable) write(nfout,*) 'CG2: Reach to max_CG_step. Go first_cg_step'
           call first_cg_step()
         endif
       else
         call linmin_set_new_cps()
       endif
     else
       f_para3 = sum(f_total0(:,:)*vec_h(:,:))/sqrt(sum(vec_h(:,:)**2))
       if(printable) write(nfout,'(a18,3f10.6)') ' CG2: F_parallel =', f_para1, f_para2, f_para3
       f_para3 = abs(f_para3)
       f_para2 = abs(f_para2)
       if ( ( f_para3 .lt. f_para2 ) .and. ( f_para3 .lt. f_para1 ) ) then
         if ( f_para3 .lt. f_para1*0.2d0 ) then
           if(printable) write(nfout,*) 'CG2: F_parallel is small. Then linmin finished.'
           if ( iter_CG .lt. iter_CG_max ) then
             if(printable) write(nfout,*) 'CG2: Go next_cg_step'
             call next_cg_step()
           else
             if(printable) write(nfout,*) 'CG2: Reach to max_CG_step. Go first_cg_step'
             call first_cg_step()
           endif
         else
           if(printable)then
           write(nfout,*) 'CG2: CG procedure is failed !!'
           write(nfout,*) 'CG2: f_para3 must be less than f_para1*0.2d0.'
           write(nfout,*) 'CG2: Then, new CG procedure starts.'
           endif
           call first_cg_step()
         endif
       else
         if(printable)then
         write(nfout,*) 'CG2: Linmin calc is failed !!'
         write(nfout,*) 'CG2: The 3rd point is not a minimum.'
         write(nfout,*) 'CG2: Then, new CG procedure starts from the 2nd point.'
         endif
         call first_cg_step_from_cps2()
       endif
     endif

  contains

    subroutine first_cg_step()
      f_total1(:,:) = f_total0(:,:)
      vec_g(:,:) = f_total1(:,:)
      vec_h(:,:) = vec_g(:,:)
      cps1(:,:) = cps(:,:)
      vec_h_norm = sqrt( sum( vec_h(1:natm,1:3)**2 ) )
      call set_alpha()
      cps(:,:) = cps1(:,:) + vec_h(:,:) * alpha
      f_para1 = sum(f_total1(:,:)*vec_h(:,:))/sqrt(sum(vec_h(:,:)**2))
      if ( f_para1 .gt. 0.08d0 ) then
        iter_CG_max = 1
      else
        iter_CG_max = 5
      endif
      iter_CG = 1
      iter_linmin = 2
      write(nfout,'(a28,i2,a18,i2)') ' CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine first_cg_step

    subroutine set_alpha()
      if      ( vec_h_norm .lt. 0.004d0 ) then
        alpha = 3.5d0
      else if ( vec_h_norm .lt. 0.007d0 ) then
        alpha = 3.0d0
      else if ( vec_h_norm .lt. 0.010d0 ) then
        alpha = 2.5d0
      else if ( vec_h_norm .lt. 0.015d0 ) then
        alpha = 2.0d0
      else if ( vec_h_norm .lt. 0.020d0 ) then
        alpha = 1.5d0
      else
        alpha = 1.0d0
      endif
    end subroutine set_alpha

    subroutine next_cg_step()
      f_total1(:,:) = f_total0(:,:)
!     gamma = sum(f_total1(:,:)**2) / sum(vec_g(:,:)**2)
      gamma = sum( (f_total1(:,:)-vec_g(:,:))*f_total1(:,:) ) / sum(vec_g(:,:)**2)
      vec_g(:,:) = f_total1(:,:)
      vec_h(:,:) = vec_g(:,:) + gamma * vec_h(:,:)
      cps1(:,:) = cps(:,:)
      vec_h_norm = sqrt( sum( vec_h(1:natm,1:3)**2 ) )
      call set_alpha()
      cps(:,:) = cps1(:,:) + vec_h(:,:) * alpha
      f_para1 = sum(f_total1(:,:)*vec_h(:,:))/sqrt(sum(vec_h(:,:)**2))
      iter_CG = iter_CG + 1
      iter_linmin = 2
     write(nfout,'(a28,i2,a18,i2)') ' CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine next_cg_step

    subroutine linmin_set_new_cps()
      cps2(:,:) = cps(:,:)
      f_total2(:,:) = f_total0(:,:)
      alpha2 = f_para1 * alpha / ( f_para1 - f_para2 )
      write(nfout,*) 'CG2: alpha(force ) =', alpha2
      if ( alpha2 .lt. 0.0d0 ) then
        write(nfout,*) 'CG2: Linmin calc is strange !!'
        write(nfout,*) 'CG2: Alpha must be positive !!'
        write(nfout,'(a13,f15.10)') ' CG2: Alpha =', alpha2
        write(nfout,*) 'CG2: New CG step start from the 2nd point.'
        call first_cg_step()
      else
        if ( alpha2 .gt. 10.0d0*alpha ) then
          write(nfout,*) 'CG2: Linmin calculated.'
          write(nfout,'(a22,f15.10)') ' CG2: Original alpha =', alpha2
          write(nfout,*) 'CG2: Alpha is too large, then, adjusted.'
          alpha2 = 10.0d0*alpha
        endif
        write(nfout,'(a13,f15.10)') ' CG2: Alpha =', alpha2
        cps(:,:) = cps1(:,:) + vec_h(:,:) * alpha2
        iter_linmin = 3
        write(nfout,'(a28,i2,a18,i2)') ' CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
      endif
    end subroutine linmin_set_new_cps

    subroutine first_cg_step_from_cps2()
      f_total1(:,:) = f_total2(:,:)
      vec_g(:,:) = f_total1(:,:)
      vec_h(:,:) = vec_g(:,:)
      cps1(:,:) = cps2(:,:)
      vec_h_norm = sqrt( sum( vec_h(1:natm,1:3)**2 ) )
      call set_alpha()
      cps(:,:) = cps1(:,:) + vec_h(:,:) * alpha
      f_para1 = sum(f_total1(:,:)*vec_h(:,:))/sqrt(sum(vec_h(:,:)**2))
      iter_CG = 1
      iter_linmin = 2
      write(nfout,'(a28,i2,a18,i2)') ' CG2: Next step is iter_CG =', iter_CG, ' and iter_linmin =', iter_linmin
    end subroutine first_cg_step_from_cps2

!! jn_130104

  end subroutine m_IS_cg2_core

  subroutine m_IS_cg2(forc_l_in,etotal)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP), intent(in)                    :: etotal
    call m_IS_cg2_core(nfout,natm,cps,forc_l_in,imdtypxyz,etotal)
  end subroutine m_IS_cg2

  subroutine m_IS_cg(forc_l_in,etotal)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP), intent(in)                    :: etotal

    logical, save :: finit = .true., fconv = .false.
    integer, save :: itime = 1
    integer :: ia,j
    integer, save :: natm_free
    real(kind=DP), save :: e0,e1,e2,e3
    real(kind=DP), save :: dt1,dt2,dt3
    real(kind=DP), save :: de1,de2,de3
    real(kind=DP), dimension(:,:), allocatable, save :: hvec,g0vec,g1vec,cps0 ! dim(natm,3)
    real(kind=DP), allocatable, dimension(:,:)   :: forc_l

    integer ::              id_sname = -1
    call tstatc0_begin('m_IS_cg ',id_sname)

    if(printable) write(nfout,*) 'gCG: iter=',iteration_ionic,' itime=',itime

    allocate(forc_l(natm,3))
    forc_l = forc_l_in

    if(.not.allocated(hvec)) allocate(hvec(natm,3))
    if(.not.allocated(g0vec)) allocate(g0vec(natm,3))
    if(.not.allocated(g1vec)) allocate(g1vec(natm,3))
    if(.not.allocated(cps0)) allocate(cps0(natm,3))

!!$    if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_fi(forc_l)
    if(constraint_type == FIXED_NORMAL_HYPERVECTOR) call modify_forc_hyperplane(forc_l)

    if(finit) then
       hvec = 0.d0
       g0vec = 0.d0
       natm_free = 0
       do ia=1,natm
          do j=1,3
          if(imdtypxyz(ia,j) == 0 ) cycle
          hvec(ia,j) = forc_l(ia,j)
          g0vec(ia,j) = hvec(ia,j)
          natm_free = natm_free + 1
          enddo
       end do
       if(printable) write(nfout,*) 'gCG: # of free atoms=',natm_free
       finit = .false.
    end if

    if(itime == 1) then
       call first()
       itime = 2
    else if(itime == 2) then
       e2 = etotal
       call deriv_energy(de2,hvec)
       if(de2 > 0.d0 .or. e2 > e0 ) then
          call estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
          itime = 3
       else
      !debug
          if(printable) write(nfout,*) 'gCG:== Second (dt2)=='
      !end debug
          e1 = e2
          de1 = de2
          dt1 = dt2
          dt2 = 2*dt1
          cpd_l = dt2 * hvec
          cps = cps0 + cpd_l
      !debug
      ! if(printable)
      ! write(nfout,*) 'gCG:cpd_l:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cpd_l(ia,1:3)
      ! end do
      ! write(nfout,*) 'gCG:cps:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cps(ia,1:3)
      ! end do
      ! end if
      !end debug
       end if
    else
       call estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
       itime = itime + 1
       if(fconv) then
          g1vec = 0.d0
          do ia=1,natm
             do j=1,3
             if(imdtypxyz(ia,j) == 0 ) cycle
             g1vec(ia,j) = forc_l(ia,j)
             enddo
          end do
          call conjugate_grad(hvec,g0vec,g1vec)
          call first()
          itime = 2
          fconv = .false.
       end if
    end if


    deallocate(forc_l)

    call tstatc0_end(id_sname)
  contains

    subroutine first()
      !debug
           if(printable) write(nfout,*) 'gCG:== First (dt1)=='
      !end debug
      cps0 = cps
      e0 = etotal
      dt1 = 0.d0
      e1 = e0
      call deriv_energy(de1,hvec)
      dt2 = 1.d0
      cpd_l = dt2 * hvec

      cps = cps0 + cpd_l
      !debug
      ! if(printable)
      ! write(nfout,*) 'gCG:cpd_l:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cpd_l(ia,1:3)
      ! end do
      ! write(nfout,*) 'gCG:cps:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cps(ia,1:3)
      ! end do
      ! end if
      !end debug
    end subroutine first

    subroutine deriv_energy(de,hvec)
      real(kind=DP), intent(out) :: de
      real(kind=DP), intent(in) :: hvec(natm,3)

      integer :: ia,j

      de = 0.d0
      do ia=1,natm
         do j=1,3
            if(imdtypxyz(ia,j) == 0 ) cycle
            de = de + forc_l(ia,j)*hvec(ia,j)
         enddo
      end do
      de = -de

    end subroutine deriv_energy

    subroutine estimate_min(dt1,dt2,e1,e2,de1,de2,fconv)
      real(kind=DP), intent(inout) :: dt1,dt2,e1,e2,de1,de2
      logical, intent(out) :: fconv

      real(kind=DP), save :: dt3,e3,de3
      real(kind=DP) :: c1,c2,c3,c4
      real(kind=DP) :: ddt,sdt,pdt2,pdt3,sq,dtm,dtp
      real(kind=DP) :: dfp,dfm,fp,fm
      logical, save :: finit = .true.

      if(.not.finit) then
         e3 = etotal
         call deriv_energy(de3,hvec)
         !debug
         if(printable) then
           write(nfout,*) 'gCG:== Estimate Emin (old dt1,dt2) =='
           write(nfout,*) 'gCG:dt1=',dt1,' e1=',e1
           write(nfout,*) 'gCG:de1=',de1
           write(nfout,*) 'gCG:dt2=',dt2,' e2=',e2
           write(nfout,*) 'gCG:de2=',de2
           write(nfout,*) 'gCG:dt3=',dt3,' e3=',e3
           write(nfout,*) 'gCG:de3=',de3
           write(nfout,*) 'gCG:de3/natm_free=',de3/dble(natm_free)
           write(nfout,*) 'gCG:etol=',etol
         end if
         !end debug
         if(abs(de3)/dble(natm_free) .le. etol) then
            fconv = .true.
            finit = .true.
            return
         end if
         if(de3 > 0.d0 .or. e3 > e1) then
            e2 = e3
            de2 = de3
            dt2 = dt3
         else
            e1 = e3
            de1 = de3
            dt1 = dt3
         end if
      end if
      finit = .false.
      !debug
      if(printable) then
        write(nfout,*) 'gCG:== Estimate Emin (new dt1,dt2) =='
        write(nfout,*) 'gCG:dt1=',dt1,' dt2=',dt2
        write(nfout,*) 'gCG: e1=',e1,' e2=',e2
        write(nfout,*) 'gCG:de1=',de1,' de2=',de2
      end if
      !end debug

      ddt = dt1 - dt2
      sdt = dt1 + dt2
      pdt2 = dt1*dt1
      pdt3 = dt1*pdt2
      c1 = ((de1+de2)*ddt-2.d0*(e1-e2))/ddt**3
      c2 = 0.5d0*(de1-de2-3.d0*c1*ddt*sdt)/ddt
      c3 = de1 - 3.d0*c1*pdt2 - 2.d0*c2*dt1
      c4 = e1 - c1*pdt3 - c2*pdt2 - c3*dt1

!!$ASASASASASAS
!!$      sq = c2+sqrt(c2*c2-3.d0*c1*c3)
      sq = c2*c2-3.d0*c1*c3
      if (sq > 0.d0) then
         sq = sqrt(sq)
      else
         write(nfout,*) 'gCG: Warning! Cannot estimate the minimum.'
         sq = 0.d0
      endif
      sq = sq + c2
!!$ASASASASASAS
      dtp = -c3/sq
      dtm = -sq/(3.d0*c1)

      dfp = 3.d0*c1*dtp**2+2.d0*c2*dtp+c3
      dfm = 3.d0*c1*dtm**2+2.d0*c2*dtm+c3

      fp = c1*dtp**3+c2*dtp**2+c3*dtp+c4
      fm = c1*dtm**3+c2*dtm**2+c3*dtm+c4

      !debug
      if(printable) then
        write(nfout,*) 'gCG:== Estimate Emin (c1,c2,c3) =='
        write(nfout,*) 'gCG: c1=',c1,' c2=',c2
        write(nfout,*) 'gCG: c3=',c3,' c4=',c4
        write(nfout,*) 'gCG:dtp=',dtp,' dtm=',dtm
        write(nfout,*) 'gCG:dfp=',dfp,' dfm=',dfm
        write(nfout,*) 'gCG: fp=',fp,' fm=',fm
      end if
      !end debug

      if(dtp >= dt1 .and. dtp <= dt2) then
         dt3 = dtp
      else if(dtm >= dt1 .and. dtm <= dt2) then
         dt3 = dtm
      else
        !!stop 'I cant estimated dt3 at which the total energy is minimum'
        write(nfout,*) 'gCG: Warning! line-minimization failed'
        dt3 = 1.d0
      end if

      cpd_l = dt3 * hvec
      cps = cps0 + cpd_l

      !debug
      if(printable) then
       write(nfout,*) 'gCG:== Estimate Emin (dt3) =='
       write(nfout,*) 'gCG:dt1=',dt1,' dt2=',dt2
       write(nfout,*) 'gCG:dtp=',dtp
       write(nfout,*) 'gCG:dtm=',dtm
       write(nfout,*) 'gCG:dt3=',dt3
      ! write(nfout,*) 'gCG:cpd_l:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cpd_l(ia,1:3)
      ! end do
      ! write(nfout,*) 'gCG:cps:'
      ! do ia=1,natm
      !    write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,cps(ia,1:3)
      ! end do
      end if
      !end debug

    end subroutine estimate_min

    subroutine conjugate_grad(hvec,g0vec,g1vec)
      real(kind=DP), dimension(natm,3), intent(inout) :: hvec,g0vec
      real(kind=DP), dimension(natm,3), intent(in) :: g1vec

      integer :: i,ia
      real(kind=DP) :: gg,gam

      gg = 0.d0
      gam = 0.d0
      do i=1,3
         do ia=1,natm
            gg = gg + g0vec(ia,i)**2
            gam = gam + (g1vec(ia,i)-g0vec(ia,i))*g1vec(ia,i)
         end do
      end do
      if(gg > 1.d-10*gam) then
         gam = gam/gg
      else
         gam = 1.d0
      end if

      hvec  = g1vec + gam * hvec
      g0vec = g1vec

      !debug
      if(printable) then
       write(nfout,*) 'gCG:== CG =='
       write(nfout,*) 'gCG:gam=',gam
       write(nfout,*) 'gCG:gg=',gg
       write(nfout,*) 'gCG: H (conjugate_grad)'
       do ia=1,natm
          write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,hvec(ia,1:3)
       end do
       write(nfout,*) 'gCG: G'
       do ia=1,natm
          write(nfout,'("gCG:",i3,3(1x,f17.9))') ia,g0vec(ia,1:3)
       end do
      end if
      !end debug

    end subroutine conjugate_grad

  end subroutine m_IS_cg

  subroutine m_IS_remove_atom(irem)
    integer, intent(in) :: irem
    logical, allocatable,dimension(:) :: excl
    call m_IS_store_current_config()
    allocate(excl(nconfig_buf));excl=.false.
    call m_IS_dealloc_pos_and_v(nfout)
    excl(irem) = .true.
    natm = natm-1
    natm2 = natm2-config_buf(irem)%iwei
    call m_IS_alloc_pos_and_v(nfout)
    call m_IS_recover(nconfig_buf,excl)
    call m_CtrlP_check_matm(nfout,natm)
    neg_incre = neg_incre - int(ceiling(config_buf(irem)%nvalence*0.5d0))
    call m_IS_set_iatom(nfout)
    deallocate(excl)
  end subroutine m_IS_remove_atom

  subroutine m_IS_add_atom(newatom)
    type(atomic_configuration_t), intent(in) :: newatom
    logical, allocatable,dimension(:) :: excl
    call m_IS_store_current_config()
    allocate(excl(nconfig_buf));excl=.false.
    call m_IS_dealloc_pos_and_v(nfout)
    natm = natm+1
    natm2 = natm2+newatom%iwei
    call m_IS_alloc_pos_and_v(nfout)
    call m_IS_recover(nconfig_buf,excl)
    ! add the new atom
    !!species_work(natm) = newatom%element
    iwei(natm) = newatom%iwei
    imdtyp(natm) = newatom%imdtyp
    imdtypxyz(natm,:) = newatom%imdtypxyz(:)
    ityp(natm) = newatom%ityp
    if_pdos(natm) = newatom%if_pdos
    if_aldos(natm) = newatom%if_aldos
    ihubbard(natm) = newatom%ihubbard
    iproj_group(natm) = newatom%iproj_group
    ionic_mass(natm) = newatom%ionic_mass
    numlay(natm) = newatom%numlay
    pos(natm,:)  = newatom%pos(:)
    cps(natm,:) = newatom%cps(:)
    pos_in(natm,:) = newatom%pos_in(:)
    cps_in(natm,:) = newatom%cps_in(:)
    cpd_l(natm,:) = newatom%cpd_l(:)
    cpo_l(natm,:,:) = newatom%cpo_l(:,:)
    neg_incre = neg_incre + int(ceiling(newatom%nvalence*0.5d0))
    call m_CtrlP_check_matm(nfout,natm)
    call m_IS_set_iatom(nfout)
    deallocate(excl)
  end subroutine m_IS_add_atom

  subroutine m_IS_supercell(nfout)
    integer, intent(in) :: nfout
    integer :: i
    if(sw_supercell==ON) then
       natm_mobile = 0
       do i=1,natm
          if(imdtyp(i)/=OFF) natm_mobile = natm_mobile+1
       enddo
       natm_prim = natm
       natm2_prim = natm2
       natm_super = nlpnt*natm
       natm2_super = nlpnt*natm2
       natm = natm_super
       natm2 = natm2_super
       call spread_atoms_on_supercell(natm_prim,natm)
       if(printable) call print_atoms_supercell(nfout)
#ifndef _EMPIRICAL_
       call m_CtrlP_check_matm(nfout,natm)
#endif
    end if
  contains
    subroutine spread_atoms_on_supercell(natm_prim,natm)
      integer, intent(in) :: natm_prim
      integer, intent(in) :: natm

      integer :: i,j,k,ia,ja
      integer, dimension(natm_prim) :: iwei_wk, imdtyp_wk, ityp_wk, if_pdos_wk, if_aldos_wk
      real(kind=DP), dimension(natm_prim) :: ionic_mass_wk
      real(kind=DP) :: rltv_t(3,3)

! =============== Added by K. Tagami  =============== 0.1
      integer, dimension(natm_prim) :: ihubbard_wk, iproj_group_wk
! ===================================================

! =============== Added by K. Tagami  =============== 11.0
      integer, dimension(natm_prim) :: itab_spinorbit_addition_wk
! =================================================== 11.0

!<<===== ASMS === 2021/03/25
      integer, allocatable :: ityp_vdw_wk(:)
!======= ASMS === 2021/03/25 >>

      integer, allocatable :: atom_key_wk(:)

!<<===== ASMS === 2021/03/25
      if (vdw_method /= VDW_DFTD3 .and. ntyp_vdw>0) then
         allocate(ityp_vdw_wk(natm_prim) );  ityp_vdw_wk = ityp_vdw
      endif
!======= ASMS === 2021/03/25 >>

      allocate( atom_key_wk(natm_prim) );  atom_key_wk = atom_key

      iwei_wk = iwei
      imdtyp_wk = imdtyp
      ityp_wk = ityp
      if_pdos_wk = if_pdos
      if_aldos_wk = if_aldos
      ionic_mass_wk = ionic_mass

! =========== Added by K. Tagami ================= 0.1
      if ( sw_hubbard == ON ) then
         iproj_group_wk = iproj_group
         ihubbard_wk = ihubbard
      endif
! ================================================ 0.1

! =============================== Added by K. Tagami ================= 11.0
      if ( SpinOrbit_Mode == ByProjector ) then
         iproj_group_wk = iproj_group
         itab_spinorbit_addition_wk = itab_spinorbit_addition
      endif
! ================================================================= 11.0

      if(nlpnt > 1) then
         do ia=1,natm_prim
            call mod1(pos(ia,1))
            call mod1(pos(ia,2))
            call mod1(pos(ia,3))
         end do
      end if
      call change_of_coordinate_system(altv_prim,pos,natm_prim,natm_prim,cps) !-(b_I.S.) pos -> cps
      allocate(pos_prim(natm_prim,3))
      allocate(cps_prim(natm_prim,3))
      pos_prim = pos
      cps_prim = cps

      call m_IS_dealloc_pos_and_v(nfout)
      call m_IS_alloc_pos_and_v(nfout)

      do i=1,nlpnt
         do ia=1,natm_prim
            ja = ia + natm_prim*(i-1)
            cps(ja,1:3) = cps_prim(ia,1:3) + lpnt(i,1:3)
            iwei(ja) = iwei_wk(ia)
            imdtyp(ja) = imdtyp_wk(ia)
            ityp(ja) = ityp_wk(ia)
            if_pdos(ja) = if_pdos_wk(ia)
            if_aldos(ja) = if_aldos_wk(ia)
            ionic_mass(ja) = ionic_mass_wk(ia)

! ==================== Added by K. Tagami ====== trial ======= 0.1
            if ( sw_hubbard == ON ) then
               ihubbard(ja) = ihubbard_wk(ia)
               iproj_group(ja) = iproj_group_wk(ia)
           endif
! ============================================================ 0.1

! =================== Added by K. Tagami ====== trial ==================== 11.0
            if ( SpinOrbit_mode == ByProjector ) then
               itab_spinorbit_addition(ja) = itab_spinorbit_addition_wk(ia)
               iproj_group(ja) = iproj_group_wk(ia)
           endif
! ======================================================================= 11.0

!<<===== ASMS === 2021/03/25
           if (vdw_method /= VDW_DFTD3 .and. ntyp_vdw>0) then
              ityp_vdw(ja) = ityp_vdw_wk(ia)
           endif
!======= ASMS === 2021/03/25 >>

           atom_key(ja) = atom_key_wk(ia)
         end do
      end do
      rltv_t = transpose(rltv)/PAI2
      call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
      do ia=1,natm
         call mod1(pos(ia,1))
         call mod1(pos(ia,2))
         call mod1(pos(ia,3))
      end do

      allocate(pos_super(natm,3))
      allocate(cps_super(natm,3))
      pos_super = pos
      cps_super = cps

      call m_IS_set_iatom(nfout) ! -> iatom
!<<===== ASMS === 2021/02/25
      if ( allocated(ityp_vdw_wk) ) deallocate( ityp_vdw_wk )
!======= ASMS === 2021/02/25 >>
      if ( allocated(atom_key_wk) ) deallocate( atom_key_wk )

    end subroutine spread_atoms_on_supercell

    subroutine mod1(t)
      real(kind=DP), intent(inout) :: t
      real(kind=DP), parameter :: eps = 1.d-6
      t = mod(t,1.d0)
      if(t < -eps) t = t + 1.d0
      if(t > 1.d0 - eps) t = t - 1.d0
    end subroutine mod1

    subroutine print_atoms_supercell(nfout)
      integer, intent(in) :: nfout
      integer :: ia
      write(nfout,'("natm_super,natm2_super=",2i4)') natm,natm2
      write(nfout,'("ia,cps(3),pos(3),ityp")')
      do ia=1,natm
         write(nfout,'(i4,6(1x,f10.5),1x,i1)') ia,cps(ia,1:3),pos(ia,1:3),ityp(ia)
      end do
    end subroutine print_atoms_supercell

  end subroutine m_IS_supercell

  subroutine m_IS_set_natm_prim
    natm = natm_prim
    natm2 = natm2_prim
    cps(1:natm,1:3) = cps_prim(1:natm,1:3)
    pos(1:natm,1:3) = pos_prim(1:natm,1:3)
    call m_IS_set_iatom(nfout) ! -> iatom
#ifndef _EMPIRICAL_
    call m_CtrlP_check_matm(nfout,natm)
#endif
  end subroutine m_IS_set_natm_prim

  subroutine m_IS_set_natm_super
    natm = natm_super
    natm2 = natm2_super
    cps(1:natm,1:3) = cps_super(1:natm,1:3)
    pos(1:natm,1:3) = pos_super(1:natm,1:3)
    call m_IS_set_iatom(nfout) ! -> iatom
#ifndef _EMPIRICAL_
    call m_CtrlP_check_matm(nfout,natm)
#endif
  end subroutine m_IS_set_natm_super

  subroutine m_IS_set_napt_prim
    if(.not.allocated(napt_prim)) return
! ========================= modified by K. Tagami ========= 0.1
!    deallocate(napt)
    if ( allocated( napt ) )  deallocate(napt)
! ======================================================== 0.1
    allocate(napt(natm_prim,nopr+af))
    napt = napt_prim
  end subroutine m_IS_set_napt_prim

  subroutine m_IS_set_napt_super
    integer :: ia
    if(.not.allocated(napt_prim)) then
       allocate(napt_prim(natm_prim,nopr+af))
       napt_prim=napt
    end if
    deallocate(napt)
    allocate(napt(natm_super,1))
    do ia=1,natm_super
       napt(ia,1) = ia
    end do
  end subroutine m_IS_set_napt_super

  subroutine m_IS_inv_sym_off(nfout)
    integer, intent(in) :: nfout

    integer :: i,n
    integer, dimension(natm2) :: imdtyp_wk, ityp_wk, if_pdos_wk, if_aldos_wk
    integer, dimension(natm2) :: ihubbard_wk, iproj_group_wk

! ========================== added by K. Tagami ================= 11.0
    integer, dimension(natm2) :: itab_spinorbit_addition_wk
! ================================================================ 11.0

    real(kind=DP), dimension(natm2,3) :: cps_wk,pos_wk
    real(kind=DP), dimension(natm2) :: ionic_mass_wk
    if(inversion_symmetry == ON) then
      if(printable) write(nfout,*) ' Inversion symmety will be OFF. (kimg=2)'
      call m_CS_set_inv_sym_off() ! -> inversion_symmetry = OFF
      cps_wk = cps(1:natm,1:3)
      pos_wk = pos(1:natm,1:3)
      imdtyp_wk = imdtyp(1:natm)
      ityp_wk = ityp(1:natm)
      if_pdos_wk = if_pdos(1:natm)
      if_aldos_wk = if_aldos(1:natm)
      ihubbard_wk = ihubbard(1:natm)

! ============================== added by K. Tagami ==================== 11.0
      itab_spinorbit_addition_wk = itab_spinorbit_addition(1:natm)
! ====================================================================== 11.0

      iproj_group_wk = iproj_group(1:natm)
      ionic_mass_wk = ionic_mass(1:natm)
      n = natm
      do i=1,natm
         if(iwei(i)==1) cycle
         n = n + 1
         cps_wk(n,1:3) = -cps(i,1:3)
         pos_wk(n,1:3) = -pos(i,1:3)
         imdtyp_wk(n) = imdtyp(i)
         ityp_wk(n) = ityp(i)
         if_pdos_wk(n) = if_pdos(i)
         if_aldos_wk(n) = if_aldos(i)
         ihubbard_wk(n) = ihubbard(i)

! =============================== added by K. Tagami =================== 11.0
         itab_spinorbit_addition_wk(n) = itab_spinorbit_addition(i)
! ====================================================================== 11.0

         iproj_group_wk(n) = iproj_group(i)
         ionic_mass_wk(n) = ionic_mass(i)
      end do
      deallocate(cps,pos,imdtyp,ityp,if_pdos,if_aldos,ihubbard,iproj_group,ionic_mass,iwei)
! ============================= added by K. Tagami ============ 11.0
      deallocate( itab_spinorbit_addition )
! ============================================================= 11.0

      natm = natm2
      allocate(cps(natm,3)); cps = cps_wk
      allocate(pos(natm,3)); pos = pos_wk
      allocate(imdtyp(natm)); imdtyp = imdtyp_wk
      allocate(ityp(natm)); ityp = ityp_wk
      allocate(if_pdos(natm)); if_pdos = if_pdos_wk
      allocate(if_aldos(natm)); if_aldos = if_aldos_wk
      allocate(ihubbard(natm)); ihubbard = ihubbard_wk

! =================================== added by K. Tagami ========== 11.0
      allocate(itab_spinorbit_addition(natm))
      itab_spinorbit_addition = itab_spinorbit_addition_wk
! ================================================================= 11.0

      allocate(iproj_group(natm)); iproj_group = iproj_group_wk
      allocate(ionic_mass(natm)); ionic_mass = ionic_mass_wk
      allocate(iwei(natm)); iwei = 1
    end if
  end subroutine m_IS_inv_sym_off

  subroutine m_IS_symmetrize_atom_pos(nfout)
#ifdef SX
!CDIR BEGIN NOVECTOR
#endif
    integer, intent(in) :: nfout

    real(kind=DP), dimension(natm2,3) :: cps_wk,cps_wk2
    real(kind=DP), dimension(natm,3) :: cpso,poso
    real(kind=DP), dimension(3,3) :: rltv_t
    real(kind=DP), dimension(3) :: p,di,dimin
    real(kind=DP) :: df,dfmin
    integer, dimension(natm2) :: ityp_wk
    integer :: i,n,ia,ja,iia
    cps_wk(1:natm,1:3)  = cps(1:natm,1:3)
    ityp_wk(1:natm) = ityp(1:natm)
    n = natm
    do i=1,natm
       if(iwei(i)==1) cycle
       n = n + 1
       cps_wk(n,1:3) = -cps(i,1:3)
       ityp_wk(n) = ityp(i)
    end do
    cps_wk2 = 0.d0
    do n=1,nopr
       do ia=1,natm2
          p = matmul(op(1:3,1:3,n),cps_wk(ia,1:3)) + tau(1:3,n,CARTS)
          iia = 0
          dfmin = 1.d10
          do ja=1,natm2
             if(ityp_wk(ia) /= ityp_wk(ja)) cycle
             di = matmul(transpose(rltv),(p - cps_wk(ja,1:3)))
             df = sum(abs(cos(di(1:3))-1.d0))
             if(df < dfmin) then
                iia = ja
                dfmin = df
                dimin = di/PAI2
             end if
          end do
          if(iia == 0) call phase_error_with_msg(nfout,'m_IS_symmetrize_atom_pos: error iia=0',__LINE__,__FILE__)
          p = p - matmul(altv,nint(dimin))
          cps_wk2(iia,1:3) = cps_wk2(iia,1:3) + p(1:3)
       end do
    end do
    cps_wk2 = cps_wk2/nopr

    cpso = cps
    poso = pos
    cps = cps_wk2(1:natm,1:3)
    rltv_t = transpose(rltv)/PAI2
    call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
    if(printable) then
       write(nfout,*) 'Atomic coordinates were symmetrized.'
       !!$write(nfout,'(20x,"Inputted Cartesian coordinate -> symmetrized Cartesian coordinate")')
       !!$do ia=1,natm
       !!$   write(nfout,'(i4,3f15.8," -> ",3f15.8)') ia,cpso(ia,1:3),cps(ia,1:3)
       !!$end do
       !!$write(nfout,'(20x,"Inputted internal coordinate  -> symmetrized internal coordinate")')
       !!$do ia=1,natm
       !!$   write(nfout,'(i4,3f15.8," -> ",3f15.8)') ia,poso(ia,1:3),pos(ia,1:3)
       !!$end do
       write(nfout,'(" === Symmetrized Cartesian coordinates and errors===")')
       do ia=1,natm
          write(nfout,'(i4,4f18.9)') ia,cps(ia,1:3),sqrt(sum((cps(ia,1:3)-cpso(ia,1:3))**2))
       end do
       write(nfout,'(" === Symmetrized internal coordinates ===")')
       do ia=1,natm
          write(nfout,'(i4,3f18.9)') ia,pos(ia,1:3)
       end do
    end if
#ifdef SX
!CDIR END
#endif
  end subroutine m_IS_symmetrize_atom_pos

  subroutine m_IS_force_af_symmetry(nfout)
#ifdef SX
!CDIR BEGIN NOVECTOR
#endif
    integer, intent(in) :: nfout

    real(kind=DP), dimension(natm2,3) :: cps_wk,cps_wk2
    real(kind=DP), dimension(natm,3) :: cpso,poso
    real(kind=DP), dimension(3,3) :: rltv_t
    real(kind=DP), dimension(3) :: p,di,dimin
    real(kind=DP) :: df,dfmin
    integer, dimension(natm2) :: ityp_wk,ityp_af
    integer :: i,n,ia,ja,iia
    cps_wk(1:natm,1:3)  = cps(1:natm,1:3)
    do i=1,natm
       ityp_af(i) = nint(iatomn(ityp(i)))
    end do
    n = natm
    do i=1,natm
       if(iwei(i)==1) cycle
       n = n + 1
       cps_wk(n,1:3) = -cps(i,1:3)
       ityp_af(n) = nint(iatomn(ityp(i)))
    end do
    cps_wk2 = 0.d0
    do ia=1,natm2
       p = matmul(op(1:3,1:3,nopr+af),cps_wk(ia,1:3)) + tau(1:3,nopr+af,CARTS)
       iia = 0
       dfmin = 1.d10
       do ja=1,natm2
          if(ityp_af(ia) /= ityp_af(ja)) cycle
          di = matmul(transpose(rltv),(p - cps_wk(ja,1:3)))
          df = sum(abs(cos(di(1:3))-1.d0))
          if(df < dfmin) then
             iia = ja
             dfmin = df
             dimin = di/PAI2
          end if
       end do
       if(iia == 0) call phase_error_with_msg(nfout,'m_IS_symmetrize_atom_pos: error iia=0',__LINE__,__FILE__)
       p = p - matmul(altv,nint(dimin))
       cps_wk2(iia,1:3) = cps_wk2(iia,1:3) + p(1:3)
       cps_wk(iia,1:3) = cps_wk2(iia,1:3)
    end do

    cpso = cps
    poso = pos
    cps = cps_wk2(1:natm,1:3)
    rltv_t = transpose(rltv)/PAI2
    call change_of_coordinate_system(rltv_t,cps,natm,natm,pos) !-(b_I.S.) cps -> pos
    if(printable.and.ipri>1) then
       write(nfout,*) 'Atomic coordinates were symmetrized.'
       !!$write(nfout,'(20x,"Inputted Cartesian coordinate -> symmetrized Cartesian coordinate")')
       !!$do ia=1,natm
       !!$   write(nfout,'(i4,3f15.8," -> ",3f15.8)') ia,cpso(ia,1:3),cps(ia,1:3)
       !!$end do
       !!$write(nfout,'(20x,"Inputted internal coordinate  -> symmetrized internal coordinate")')
       !!$do ia=1,natm
       !!$   write(nfout,'(i4,3f15.8," -> ",3f15.8)') ia,poso(ia,1:3),pos(ia,1:3)
       !!$end do
       write(nfout,'(" === Symmetrized Cartesian coordinates and errors===")')
       do ia=1,natm
          write(nfout,'(i4,7f18.9)') ia,cpso(ia,1:3),cps(ia,1:3),sqrt(sum((cps(ia,1:3)-cpso(ia,1:3))**2))
       end do
       write(nfout,'(" === Symmetrized internal coordinates ===")')
       do ia=1,natm
          write(nfout,'(i4,6f18.9)') ia,poso(ia,1:3),pos(ia,1:3)
       end do
    end if
#ifdef SX
!CDIR END
#endif
  end subroutine m_IS_force_af_symmetry

  subroutine m_IS_dealloc(neb_mode)
    logical, intent(in), optional :: neb_mode
    logical :: neb
    neb = .false.
    if(present(neb_mode)) neb = neb_mode
    if(allocated(napt)) deallocate(napt)
    if(allocated(napt_tl)) deallocate(napt_tl)
    if(allocated(fxyzew_l)) deallocate(fxyzew_l)

    if(allocated(zfm3_l)) deallocate(zfm3_l)

    if(.not.neb)then
       call m_IS_dealloc_pos_and_v(nfout)
       call dealloc_species_vdw_work()
       call dealloc_speciesname()
       call m_IS_dealloc_iatomn_etc()

! == KT_add ============================== 2013/10/31
       if ( noncol ) call m_IS_dealloc_magmom_local
! ======================================= 2013/10/31

       call T_control_dealloc()

! ===================================== KT_add ================ 13.0B
       call dealloc_speciesname_vdw()
       call m_IS_dealloc_vdw()
! ============================================================= 13.0B

    endif
    call dealloc_supercell_symmetry()
  end subroutine m_IS_dealloc

  subroutine m_IS_vdwdf3(nfout)
    integer, intent(in) :: nfout
    real(kind=DP) :: esum
    integer       :: id_sname = -1
    integer       :: neibrd, alen(3)
    real(kind=DP), allocatable :: rxyz(:,:)    ! d(neibrd,3)
    real(kind=DP), allocatable :: rr(:)        ! d(neibrd)
    real(kind=DP), pointer, dimension(:,:) :: cps_fp
    integer,       pointer, dimension(:)   :: ityp_full ! d(natm2)
    integer :: nu,u,v
    real(kind=DP) :: rc2,r2
    integer, allocatable, dimension(:) :: nnei
    real(kind=DP), allocatable, dimension(:,:) :: rlist
    integer, allocatable, dimension(:,:) :: indlist
    real(kind=DP), allocatable, dimension(:,:,:) :: drs
    real(kind=DP), allocatable, dimension(:) :: ncs
    real(kind=DP), allocatable, dimension(:,:) :: dncs
    real(kind=DP), allocatable, dimension(:) :: dncsdr
    real(kind=DP), allocatable, dimension(:,:) :: forc
    real(kind=DP) :: c6,c8,fdamp6,fdamp8,r,r6,r8,dcdr6,dcdr8,dfdamp6,dfdamp8,ddr
    real(kind=DP), dimension(3) :: dc6,dc8,ftmp
    real(kind=DP), dimension(3) :: dr,tmpar1,tmpar2
    integer :: in,ielem,jelem,kelem,ia,iai,iaj,icount,neimax,neimaxatm
    integer :: ia1,ia2,ia3,i3
    real(kind=DP) :: rab,rac,rbc,c9,cosa,cosb,cosc, ctmp
    real(kind=DP),dimension(3) :: rabvec,racvec,rbcvec

    real(kind=DP) :: dncsdr_wk, dncsdr1_wk, dncsdr2_wk, vec1(3), rcut_max
    integer :: iaj1, iaj2, jelem1, jelem2

    call tstatc0_begin('m_IS_vdwdf3 ',id_sname)
!    call decide_rxyz_size(dftd3par%rcut_vdw,alen,neibrd) !-> alen, neibrd

    rcut_max = max( dftd3par%rcut_vdw, dftd3par%rcut_nc )
    call decide_rxyz_size( rcut_max, alen,neibrd) !-> alen, neibrd

    allocate(rxyz(neibrd,3))
    allocate(rr(neibrd))
    call substitute_rxyz(alen, neibrd, rxyz, rr)
    deallocate(rr)
    allocate(cps_fp(natm2,3))
    allocate(ityp_full(natm2))
    allocate(nnei(ista_atm:iend_atm));nnei=0
    cps_fp(1:natm,1:3) = pos(1:natm,1:3)
    ityp_full(1:natm)  = ityp(1:natm)
    call rplcps(cps_fp,ityp_full,1,natm2,natm,iwei)
    call cpspac ! -> cps_fp

!!    rc2 = dftd3par%rcut_vdw**2
    rc2 = rcut_max **2

    neimax = 0
    do nu = ista_atm, iend_atm ! MPI
       neimaxatm = 0
       do in = 1, neibrd
          do ia=1,natm2
             if(in == 1 .and. ia == nu) cycle
             dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
             r2 = dot_product(dr,dr)
             if(r2 > rc2) cycle
             neimaxatm = neimaxatm+1
          enddo
       enddo
       if(neimaxatm>neimax) neimax = neimaxatm
    enddo
    allocate(rlist(neimax,ista_atm:iend_atm))
    allocate(indlist(neimax,ista_atm:iend_atm))
    allocate(drs(neimax,ista_atm:iend_atm,3));drs=0.d0
    do nu = ista_atm, iend_atm ! MPI
       do in = 1, neibrd
          do ia=1,natm2
             if(in == 1 .and. ia == nu) cycle
             dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
             r2 = dot_product(dr,dr)
             if(r2 > rc2) cycle
             nnei(nu) = nnei(nu)+1
             rlist(nnei(nu),nu) = sqrt(r2)
             indlist(nnei(nu),nu) = ia
             drs(nnei(nu),nu,1:3) = dr(1:3)
          enddo
       enddo
    enddo
    allocate(ncs(natm));ncs=0
    allocate(dncs(natm,3));dncs=0.d0
    allocate(dncsdr(natm));dncsdr=0.d0
    allocate(forc(natm,3));forc=0.d0
    do ia=ista_atm,iend_atm
       call eval_nc_and_der(ia,neimax,nnei,rlist,indlist,drs,ncs(ia),tmpar2,dncsdr(ia))
       dncs(ia,1:3) = tmpar2(1:3)
    enddo
    if(npes>1) then
       call mpi_allreduce(MPI_IN_PLACE,ncs,natm,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,dncs,natm*3,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,dncsdr,natm,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
    endif
    if(iprivdw>1) then
       do ia=1,natm
          write(nfout,'(a,i8,a,d15.5)') '!** coordination number for atom ',ia,' : ',ncs(ia)
       enddo
    endif
    esum = 0.d0
    forc = 0.d0
    s_vdw = 0.0d0

#ifdef __EDA__
    if ( sw_eda == ON ) evdw_on_atom = 0.0d0
#endif

    do iai = ista_atm, iend_atm ! MPI
       ielem = iatomn(ityp(iai))
       do ia=1,nnei(iai)
          iaj = indlist(ia,iai)
          jelem = iatomn(ityp(iaj))

          tmpar1(1:3) = dncs(iai,1:3)
          call eval_deriv_nc_off_diagonal( iai, ia, neimax, rlist, indlist, &
         &                                 drs, tmpar2, dncsdr_wk )
!          call eval_c6_and_der( ielem, jelem, ncs(iai), ncs(iaj), tmpar1, tmpar2, &
!               &                dncsdr(iai), dncsdr_wk, c6, dc6, dcdr6 )
          call eval_c6_and_der_mod( ielem, jelem, ncs(iai), ncs(iaj), tmpar1, tmpar2, &
               &                    dncsdr(iai), dncsdr_wk, c6, dc6, dcdr6 )
          call eval_c8_and_der( ielem,jelem,c6,dc6,dcdr6,c8,dc8,dcdr8 )

          r = rlist(ia,iai)
          if ( r > dftd3par%rcut_vdw ) cycle

          r6 = r**6
          r8 = r6*r*r
          call eval_fdamp_and_der(6,r,ielem,jelem,fdamp6,dfdamp6)
          call eval_fdamp_and_der(8,r,ielem,jelem,fdamp8,dfdamp8)
          !!fdamp6=1.d0;dfdamp6=0.d0;fdamp8=1.d0;dfdamp8=0.d0

          ctmp = (fdamp6 * dftd3par%s6*c6/r6 &
          &           +  fdamp8 * dftd3par%s8*c8/r8)*iwei(iai)
          esum = esum +ctmp

          ftmp(1:3) = dftd3par%s6*((-6.d0/(r6*r))*(drs(ia,iai,1:3)/r)*c6*fdamp6 &
          &         + (1.d0/r6)*dc6(1:3)*fdamp6 &
          &         + (c6/r6)*dfdamp6*drs(ia,iai,1:3)/r)     &
          &         + dftd3par%s8*((-8.d0/(r8*r))*(drs(ia,iai,1:3)/r)*c8*fdamp8 &
          &         + (1.d0/r8)*dc8(1:3)*fdamp8 &
          &         + (c8/r8)*dfdamp8*drs(ia,iai,1:3)/r)

          ddr       = dftd3par%s6*((-6.d0/(r6*r))*c6*fdamp6 &
          &         + (1.d0/r6)*dcdr6*fdamp6 &
          &         + (c6/r6)*dfdamp6)  &
          &         + dftd3par%s8*((-8.d0/(r8*r))*c8*fdamp8 &
          &         + (1.d0/r8)*dcdr8*fdamp8 &
          &         + (c8/r8)*dfdamp8)

          forc(iai,1:3) = forc(iai,1:3) + ftmp(1:3)
!!          forc(iaj,1:3) = forc(iaj,1:3) - ftmp(1:3)
          do u=1,3
             do v=1,3
! ---> ASMS 2024/04/03 test
!!                s_vdw(u,v) = s_vdw(u,v) + ddr *drs(ia,iai,u) *drs(ia,iai,v)/r
                s_vdw(u,v) = s_vdw(u,v) + ftmp(u) *drs(ia,iai,v)
! <--- ASMS 2024/04/03 test
             end do
          end do

#ifdef __EDA__
          if ( sw_eda == ON ) evdw_on_atom(iai) = evdw_on_atom(iai) +ctmp
#endif
       enddo
    enddo

    do iai = ista_atm, iend_atm ! MPI
       ielem = iatomn(ityp(iai))

       do ia1=1,nnei(iai)
          iaj1 = indlist(ia1,iai)
          jelem1 = iatomn(ityp(iaj1))

          do ia2=1,nnei(iai)
             if ( ia1 == ia2 ) cycle
!!!!!!!!!             if ( ia2 < ia1 ) cycle

             iaj2 = indlist(ia2,iai)
             jelem2 = iatomn(ityp(iaj2))

             call eval_deriv_nc_off_diagonal( iai, ia1, neimax, rlist, indlist, &
                  &                           drs, tmpar1, dncsdr1_wk )
             call eval_deriv_nc_off_diagonal( iai, ia2, neimax, rlist, indlist, &
                  &                           drs, tmpar2, dncsdr2_wk )

!            call eval_c6_and_der( jelem1, jelem2, ncs(iaj1), ncs(iaj2), &
!                  &                tmpar1, tmpar2, dncsdr1_wk, dncsdr2_wk, &
!                  &                c6, dc6, dcdr6 )
             call eval_c6_and_der_mod( jelem1, jelem2, ncs(iaj1), ncs(iaj2), &
                  &                    tmpar1, tmpar2, dncsdr1_wk, dncsdr2_wk, &
                 &                    c6, dc6, dcdr6 )
             call eval_c8_and_der( jelem1, jelem2,c6,dc6,dcdr6,c8,dc8,dcdr8 )

             vec1(:) = drs(ia2,iai,:) -drs(ia1,iai,:)
             r = sqrt( dot_product( vec1, vec1 ) )
             if ( r > dftd3par%rcut_vdw ) cycle

             r6 = r**6
             r8 = r6*r*r
             call eval_fdamp_and_der(6,r,jelem1,jelem2,fdamp6,dfdamp6)
             call eval_fdamp_and_der(8,r,jelem1,jelem2,fdamp8,dfdamp8)
             !!fdamp6=1.d0;dfdamp6=0.d0;fdamp8=1.d0;dfdamp8=0.d0

             ftmp(1:3) = dftd3par%s6*( (1.d0/r6)*dc6(1:3)*fdamp6 ) &
                  &     + dftd3par%s8*( (1.d0/r8)*dc8(1:3)*fdamp8 )

             ddr       = dftd3par%s6*((1.d0/r6)*dcdr6*fdamp6 ) &
                  &    + dftd3par%s8*((1.d0/r8)*dcdr8*fdamp8 )

             forc(iai,1:3) = forc(iai,1:3) + ftmp(1:3) /2.0

             do u=1,3
                do v=1,3
! ---> ASMS 2024/04/03 test
!                   s_vdw(u,v) = s_vdw(u,v) + ddr *drs(ia1,iai,u) &
!                        &                   *drs(ia2,iai,v)/r
                   s_vdw(u,v) = s_vdw(u,v) + ftmp(u) *drs(ia2,iai,v)
! <--- ASMS 2024/04/03 test
                end do
             end do
          end do
       enddo
    enddo

1000 continue

    if(npes>1) then
       call mpi_allreduce(MPI_IN_PLACE,esum,1,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,forc,3*natm,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
       call mpi_allreduce(MPI_IN_PLACE,s_vdw,9,mpi_double_precision,mpi_sum,mpi_ke_world,ierr)
#ifdef __EDA__
       if ( sw_eda == ON ) then
          call mpi_allreduce( MPI_IN_PLACE, evdw_on_atom, natm, mpi_double_precision, &
               &              mpi_sum, mpi_ke_world, ierr )
       endif
#endif
    endif

#ifdef __EDA__
    if ( sw_eda == ON ) evdw_on_atom = -evdw_on_atom *0.5d0
#endif

    evdw = -esum*0.5d0
!!    fxyzvdw_l = forc*0.5d0
    fxyzvdw_l = forc
    s_vdw = s_vdw*0.5d0/ univol

    if(iprivdw>=2)then
       write(nfout,'(a,d15.5)') ' !** dftd3 energy ',evdw
       write(nfout,'(a)')       ' !** dftd3 forces'
       do ia=1,natm
          write(nfout,'(i8,3d15.5)') ia,fxyzvdw_l(ia,1:3)
       enddo
       write(nfout,'(a)')       ' !** dftd3 stress tensor'
       do u=1,3
          write(nfout,'(3d15.5)') s_vdw(u,1:3)
       enddo
    endif
    call flush(nfout)
    deallocate(rxyz)
    deallocate(drs)
    deallocate(cps_fp)
    deallocate(ityp_full)
    deallocate(nnei)
    deallocate(rlist)
    deallocate(indlist)
    deallocate(ncs)
    deallocate(dncs)
    deallocate(forc)
    deallocate(dncsdr)
    call tstatc0_end(id_sname)

  contains

    subroutine eval_c6_and_der(ielem,jelem,c1,c2,dc1,dc2,dcr1,dcr2,resc6,dresc6,dresdrc6)
       integer, intent(in) :: ielem,jelem
       real(kind=DP),intent(in) :: c1,c2
       real(kind=DP),intent(in),dimension(3) :: dc1,dc2
       real(kind=DP),intent(in) :: dcr1,dcr2
       real(kind=DP),intent(out) :: resc6
       real(kind=DP),intent(out),dimension(3) :: dresc6
       real(kind=DP),intent(out) :: dresdrc6
       integer :: i,j
       real(kind=DP) :: w,z,c6,n1,n2,dtmp,dtmpr,dwr,dzr,wi,resc6wi
       real(kind=DP), dimension(3) :: ddtmp,dw,dz
       !!$real(kind=DP) :: eps=1.d-10
       real(kind=DP) :: eps=0.d0
       z = 0.d0;w=0.d0;dw=0.d0;dz=0.d0;dwr=0.d0;dzr=0.d0
       do i=1,dftd3par%maxnc(ielem)
          do j=1,dftd3par%maxnc(jelem)
             c6 = dftd3par%c6ab(ielem,jelem,i,j)
             n1 = dftd3par%nc1(ielem,jelem,i,j)
             n2 = dftd3par%nc2(ielem,jelem,i,j)
             dtmp = exp(-dftd3par%k3*((c1-n1)**2+(c2-n2)**2))
             dtmpr =      dtmp*(-dftd3par%k3)*(2.d0*(c1-n1)*dcr1    +2.d0*(c2-n2)*dcr2)
             ddtmp(1:3) = dtmp*(-dftd3par%k3)*(2.d0*(c1-n1)*dc1(1:3)+2.d0*(c2-n2)*dc2(1:3))
             w = w+dtmp
             z = z+c6*dtmp
             dw = dw+ddtmp
             dz = dz+c6*ddtmp
             dwr = dwr+dtmpr
             dzr = dzr+c6*dtmpr
          enddo
       enddo
       if(w.gt.eps) then
          wi = 1.d0/w
          resc6 = z*wi
          resc6wi  =  resc6*wi
          dresc6   = -resc6wi*dw+wi*dz
          dresdrc6 = -resc6wi*dwr+wi*dzr
       endif
    end subroutine eval_c6_and_der

    subroutine eval_c6_and_der_mod( ielem, jelem, c1, c2, dc1, dc2, dcr1, dcr2, &
         &                          resc6, dresc6, dresdrc6 )
      integer, intent(in) :: ielem,jelem
      real(kind=DP),intent(in) :: c1,c2
      real(kind=DP),intent(in),dimension(3) :: dc1,dc2
      real(kind=DP),intent(in) :: dcr1,dcr2
      real(kind=DP),intent(out) :: resc6
      real(kind=DP),intent(out),dimension(3) :: dresc6
      real(kind=DP),intent(out) :: dresdrc6
      integer :: i,j
      real(kind=DP) :: w,z,c6,n1,n2,dtmp,dtmpr,dwr,dzr,wi,resc6wi
      real(kind=DP) :: exponent_wk, exponent_min, cwk
      real(kind=DP), dimension(3) :: ddtmp,dw,dz
       !!$real(kind=DP) :: eps=1.d-10
!       real(kind=DP) :: eps=0.d0

      z = 0.d0;w=0.d0;dw=0.d0;dz=0.d0;dwr=0.d0;dzr=0.d0
      
      exponent_min = 1.0D90
      do i=1,dftd3par%maxnc(ielem)
         do j=1,dftd3par%maxnc(jelem)
            n1 = dftd3par%nc1(ielem,jelem,i,j)
            n2 = dftd3par%nc2(ielem,jelem,i,j)
            exponent_wk = (c1-n1)**2 +(c2-n2)**2
            exponent_min = min( exponent_min, exponent_wk )
         end do
      end do
      
      do i=1,dftd3par%maxnc(ielem)
         do j=1,dftd3par%maxnc(jelem)
            c6 = dftd3par%c6ab(ielem,jelem,i,j)
            n1 = dftd3par%nc1(ielem,jelem,i,j)
            n2 = dftd3par%nc2(ielem,jelem,i,j)
            
            cwk = (c1-n1)**2+(c2-n2)**2 -exponent_min
            
            dtmp = exp( -dftd3par%k3 *cwk )
            dtmpr =      dtmp*(-dftd3par%k3)*(2.d0*(c1-n1)*dcr1 +2.d0*(c2-n2)*dcr2)
            ddtmp(1:3) = dtmp*(-dftd3par%k3)*(2.d0*(c1-n1)*dc1(1:3) &
                 &                           +2.d0*(c2-n2)*dc2(1:3))
            w   = w   +dtmp;      z   = z   +c6*dtmp
            dw  = dw  +ddtmp;     dz  = dz  +c6*ddtmp
            dwr = dwr +dtmpr;     dzr = dzr +c6*dtmpr
         enddo
      enddo

      wi = 1.d0/w
      resc6 = z *wi
      resc6wi  =  resc6 *wi
      dresc6   = -resc6wi *dw +wi *dz
      dresdrc6 = -resc6wi *dwr +wi *dzr

    end subroutine eval_c6_and_der_mod

    subroutine eval_c8_and_der(ielem,jelem,c6,dc6,dcdr6,c8,dc8,dcdr8)
       integer, intent(in) :: ielem,jelem
       real(kind=DP), intent(in) :: c6
       real(kind=DP), intent(in), dimension(3) :: dc6
       real(kind=DP), intent(in) :: dcdr6
       real(kind=DP), intent(out) :: c8
       real(kind=DP), intent(out), dimension(3) :: dc8
       real(kind=DP), intent(out) :: dcdr8
       real(kind=DP) :: res
       c8 = 3.d0*c6*dftd3par%r2r4(ielem)*dftd3par%r2r4(jelem)
       dc8 = 3.d0*dc6*dftd3par%r2r4(ielem)*dftd3par%r2r4(jelem)
       dcdr8 = 3.d0*dcdr6*dftd3par%r2r4(ielem)*dftd3par%r2r4(jelem)
    end subroutine eval_c8_and_der

    subroutine eval_nc_and_der(ia,neimax,nnei,rlist,indlist,drs,nca,dnca,dncadr)
       integer, intent(in) :: ia,neimax
       integer, dimension(ista_atm:iend_atm) :: nnei
       real(kind=DP), intent(in), dimension(neimax,ista_atm:iend_atm) :: rlist
       integer, intent(in), dimension(neimax,ista_atm:iend_atm) :: indlist
       real(kind=DP), intent(in), dimension(neimax,ista_atm:iend_atm,3) :: drs
       real(kind=DP), intent(out) :: nca
       real(kind=DP), intent(out), dimension(3) :: dnca
       real(kind=DP), intent(out) :: dncadr
       real(kind=DP) :: rab,racov,rbcov,dtmp,etmp,rtmp,rinv
       integer :: i,ielem,iaelem,itmp

       real(kind=DP) :: rr, scale_tmp
       scale_tmp = dftd3par%k2 / (4.0d0/3.0d0 )

       nca = 0.d0
       dnca = 0.d0
       dncadr = 0.d0
       iaelem = nint(iatomn(ityp(ia)))
       racov = dftd3par%covrad(iaelem) *scale_tmp

       do i=1,nnei(ia)
          itmp = indlist(i,ia)
          ielem = nint(iatomn(ityp_full(itmp)))
          rbcov = dftd3par%covrad(ielem) *scale_tmp

          rr = rlist(i,ia)
          if ( rr > dftd3par%rcut_nc ) cycle

          rinv = 1.d0 /rr

          rtmp = (racov+rbcov)*rinv
          etmp = exp(-dftd3par%k1*(rtmp-1.d0))
          dtmp = 1.d0+etmp
          nca = nca+1.d0/dtmp
          dtmp = -(1.d0/dtmp**2)*etmp*dftd3par%k1*rtmp*rinv
          dnca(1:3) = dnca(1:3)+dtmp*drs(i,ia,1:3)*rinv
          dncadr = dncadr + dtmp
       enddo
    end subroutine eval_nc_and_der

    subroutine eval_deriv_nc_off_diagonal( ia, i, neimax, rlist, indlist, &
         &                                 drs, dnca, dncadr )
      integer, intent(in) :: ia, neimax, i
      real(kind=DP), intent(in), dimension(neimax,ista_atm:iend_atm) :: rlist
      integer, intent(in), dimension(neimax,ista_atm:iend_atm) :: indlist
      real(kind=DP), intent(in), dimension(neimax,ista_atm:iend_atm,3) :: drs
      real(kind=DP), intent(out), dimension(3) :: dnca
      real(kind=DP), intent(out) :: dncadr
      real(kind=DP) :: rab,racov,rbcov,dtmp,etmp,rtmp,rinv
      integer :: ielem,iaelem,itmp

      real(kind=DP) :: rr, scale_tmp
      scale_tmp = dftd3par%k2 / (4.0d0/3.0d0 )

      dnca = 0.d0
      dncadr = 0.d0
      iaelem = nint(iatomn(ityp(ia)))
      racov = dftd3par%covrad(iaelem) *scale_tmp
      
      itmp = indlist(i,ia)
      ielem = nint(iatomn(ityp_full(itmp)))
      rbcov = dftd3par%covrad(ielem) *scale_tmp
      
      rr = rlist(i,ia)
      if ( rr > dftd3par%rcut_nc ) return
      
      rinv = 1.d0 /rr

      rtmp = (racov+rbcov)*rinv
      etmp = exp(-dftd3par%k1*(rtmp-1.d0))
      dtmp = 1.d0+etmp
      
      dtmp = -(1.d0/dtmp**2)*etmp*dftd3par%k1*rtmp*rinv
      dnca(1:3) = dtmp*drs(i,ia,1:3)*rinv
      dncadr = dtmp
      
    end subroutine eval_deriv_nc_off_diagonal

    subroutine eval_fdamp_and_der(norder,r,ielem,jelem,fdamp,dfdamp)
       integer, intent(in) :: norder
       real(kind=DP), intent(in) :: r
       integer, intent(in) :: ielem,jelem
       real(kind=DP), intent(out) :: fdamp
       real(kind=DP), intent(out) :: dfdamp
       real(kind=DP) :: res,sr,alpha,srr0ab
       real(kind=DP) :: c1, c2, rad_bj

       select case ( dftd3_damping_function )
       case (0)
          sr = dftd3par%sr6;   alpha = dftd3par%alpha6
          if (norder == 8) then
             sr = dftd3par%sr8;    alpha = dftd3par%alpha8
          endif
          srr0ab = sr*dftd3par%r0ab(ielem,jelem)
          fdamp  = 1.d0/(1.d0+6.d0*(r/srr0ab)**(-alpha))
          dfdamp = fdamp*fdamp*(6.d0*alpha/srr0ab**(-alpha))*r**(-alpha-1)

       case (1)
          rad_bj = sqrt( 3.d0 *dftd3par%r2r4(ielem) *dftd3par%r2r4(jelem) )
          c1 = dftd3par%a1 *rad_bj +dftd3par%a2
          c2 = ( c1 /r )**norder
          fdamp = 1.d0 /( 1.0d0 +c2 )
          dfdamp = fdamp**2 *norder *c2 /r

       end select

    end subroutine eval_fdamp_and_der

    subroutine cpspac
      real(kind=DP), dimension(3) :: catoms(3)
      integer                     :: i
      do i = 1, natm2
         catoms = cps_fp(i,1:3)
         catoms = catoms - nint(catoms)      !Packing
         cps_fp(i,1:3) = matmul(altv,catoms) !Change of coordinate system
      end do
    end subroutine cpspac

  end subroutine m_IS_vdwdf3

  subroutine m_IS_vdw(nfout)
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!   use m_Control_Parameters, only : rcut_vdw
    use m_Control_Parameters, only : vdw_method, vdw_radius, &
                                     vdw_scaling_factor, vdw_scaling_factor_r, vdw_damping_factor
    use m_Const_Parameters, only : VDW_WILLIAMS, VDW_GRIMME
! ==============================================================================
    integer, intent(in) :: nfout

    integer       :: nu,ia,in, u, v
    integer       :: it1,it2
    integer       :: neibrd,  alen(3)
    real(kind=DP) :: dr(3),r,r2,r6,r8,x,x7
    real(kind=DP) :: fx, gx, exp1, exp2, exp3
    real(kind=DP) :: c6, r0, fac, frc(3), esum, rc2
    integer,       allocatable :: ityp_vdw_full(:) ! d(natm2)
    real(kind=DP), allocatable :: rxyz(:,:)    ! d(neibrd,3)
    real(kind=DP), allocatable :: rr(:)        ! d(neibrd)
    real(kind=DP), allocatable :: cps_fp(:,:)  ! d(natm2,3)
    real(kind=DP), allocatable :: fxyzvdw_mpi(:,:) ! d(natm,3)

! ================================ KT_add ===================== 13.0B
    real(kind=DP), allocatable :: s_vdw_mpi(:,:) ! d(3,3)
! ============================================================= 13.0B

    integer       :: id_sname = -1
    call tstatc0_begin('m_IS_vdw ',id_sname)

! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!   rc2 = rcut_vdw**2
    rc2 = vdw_radius**2
! ==============================================================================

! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!   call decide_rxyz_size(rcut_vdw,alen,neibrd) !-> alen, neibrd
    call decide_rxyz_size(vdw_radius,alen,neibrd) !-> alen, neibrd
! ==============================================================================
    allocate(rxyz(neibrd,3))
    allocate(rr(neibrd))
    call substitute_rxyz(alen, neibrd, rxyz, rr)
    deallocate(rr)
    allocate(cps_fp(natm2,3))
    allocate(ityp_vdw_full(natm2))
    cps_fp(1:natm,1:3) = pos(1:natm,1:3)
    ityp_vdw_full(1:natm)  = ityp_vdw(1:natm)
    call rplcps(cps_fp,ityp_vdw_full,1,natm2,natm,iwei)
    call cpspac ! -> cps_fp

    evdw = 0.d0
    fxyzvdw_l = 0.d0
! ================================= KT_add ================== 13.0B
    s_vdw = 0.0d0
! =========================================================== 13.0B

#ifdef __EDA__
    if ( sw_eda == ON ) evdw_on_atom = 0.0d0
#endif

    do nu = 1,natm
       it1 = ityp_vdw(nu)
       esum = 0.d0
       do in = 1, neibrd
          do ia=1,natm2
             if(in == 1 .and. ia == nu) cycle
             it2 = ityp_vdw_full(ia)
             dr(1:3) = cps_fp(nu,1:3) - cps_fp(ia,1:3) - rxyz(in,1:3)
             r2 = dot_product(dr,dr)
             if(r2 > rc2) cycle
! === Apply modifications for vdW function. by tkato 2012/06/14 ================
!            r  = sqrt(r2)
!            r6 = r2*r2*r2
!            r8 = r6*r2
!            r0 = rvdw(it1,it2)
!            x = r/r0
!            x7 = x**7.d0
!            exp1 = exp(-3.d0*x7)
!            exp2 = 1.d0-exp1
!            exp3 = exp2*exp2
!            fx = exp3*exp3
!            gx = 4.d0*7.d0*3.d0*x7*exp1/exp2
!            c6 = cvdw(it1,it2)
!            esum = esum + fx * c6 / r6
!            fac = -c6*(gx-6.d0)*fx / r8
!            fxyzvdw_l(nu,1:3) = fxyzvdw_l(nu,1:3) + fac * dr(1:3)
! !!          write(nfout,*) 'nu,ia,it1,it2=',nu,ia,it1,it2
! !!          write(nfout,*) 'c6,r0=',c6,r0
! !!          write(nfout,*) 'x,fx,gx=',x,fx,gx
! !!          write(nfout,*) 'frc=', fac * dr(1:3)
! !!          write(nfout,*) 'dr=',dr(1:3)

             select case(vdw_method)
             case(VDW_WILLIAMS)
               r  = sqrt(r2)
               r6 = r2*r2*r2
               r8 = r6*r2
               r0 = rvdw(it1,it2)
               x = r/r0
               x7 = x**7.d0
               !exp1 = exp(-3.d0*x7)
               exp1 = exp(-vdw_damping_factor*x7)
               exp2 = 1.d0-exp1
               exp3 = exp2*exp2
               fx = exp3*exp3
               !gx = 4.d0*7.d0*3.d0*x7*exp1/exp2
               gx = 4.d0*7.d0*vdw_damping_factor*x7*exp1/exp2
               c6 = cvdw(it1,it2)
               esum = esum + fx * c6 / r6
               fac = -c6*(gx-6.d0)*fx / r8
               fxyzvdw_l(nu,1:3) = fxyzvdw_l(nu,1:3) + fac * dr(1:3)
  !!           write(nfout,*) 'nu,ia,it1,it2=',nu,ia,it1,it2
  !!           write(nfout,*) 'c6,r0=',c6,r0
  !!           write(nfout,*) 'x,fx,gx=',x,fx,gx
  !!           write(nfout,*) 'frc=', fac * dr(1:3)
  !!           write(nfout,*) 'dr=',dr(1:3)

             case(VDW_GRIMME)
               r  = sqrt(r2)
               r6 = r2*r2*r2
               r0 = rvdw(it1,it2)
               x = r/r0
               exp1 = exp(-vdw_damping_factor*(x-1.0d0))
               fx = 1.0d0 / (1.0d0 + exp1)
               gx = vdw_damping_factor*exp1*fx**2 / r0
               c6 = cvdw(it1,it2)
               esum = esum + (-vdw_scaling_factor) * fx * c6 / r6

! ============================== KT_mod ================================= 13.0B
!!               fac = +c6 * vdw_scaling_factor * (gx - 6.0d0*fx/r) / r6
!!               fxyzvdw_l(nu,1:3) = fxyzvdw_l(nu,1:3) + fac * dr(1:3)/r

               fac = +c6 * vdw_scaling_factor * (gx - 6.0d0*fx/r) / r6 /r
               fxyzvdw_l(nu,1:3) = fxyzvdw_l(nu,1:3) + fac * dr(1:3)
! ======================================================================= 13.0B

! =========================== KT_add ======================== 13.0B
               Do u=1, 3
                  Do v=1, 3
                     s_vdw(u,v) = s_vdw(u,v) +fac *dr(u) *dr(v)
                  End do
               End do
! ============================================================ 13.0B

             end select
! ==============================================================================
          end do
       end do
       evdw = evdw + iwei(nu) * esum
#ifdef __EDA__
       if ( sw_eda == ON ) evdw_on_atom(nu) = evdw_on_atom(nu) +iwei(nu) *esum
#endif
    end do
    evdw = 0.5d0*evdw

#ifdef __EDA__
    if ( sw_eda == ON ) evdw_on_atom = evdw_on_atom /2.0d0
#endif

! ================================ KT_add ===================== 13.0B
!    s_vdw = s_vdw / 2.0d0 / univol
    s_vdw = s_vdw / 2.0d0 / univol            ! ASMS 2015/03/24
! ============================================================= 13.0B

    deallocate(rxyz)
    deallocate(cps_fp)
    deallocate(ityp_vdw_full)

    !debug
    !write(nfout,'(" debug -- evdw --")')
    !write(nfout,'("evdw=",f16.8)') evdw
    !write(nfout,'(" debug -- fxyzvdw_l --")')
    !do i = 1, natm
    !   write(nfout,'(i5,3f16.8)') i,(fxyzvdw_l(i,j),j=1,3)
    !end do
    !!stop 'debug vdw'
    !end debug

    call tstatc0_end(id_sname)
  contains
    subroutine cpspac
      real(kind=DP), dimension(3) :: catoms(3)
      integer                     :: i
      do i = 1, natm2
         catoms = cps_fp(i,1:3)
         catoms = catoms - nint(catoms)      !Packing
         cps_fp(i,1:3) = matmul(altv,catoms) !Change of coordinate system
      end do
    end subroutine cpspac
  end subroutine m_IS_vdw

  subroutine check_region(nfout)
    integer, intent(in) :: nfout
    integer :: i,iat,iatt
    real(kind=DP) :: ratio = 1.5d0
    real(kind=DP) :: sigma,tmpr1,tmpr2
    logical :: too_close
    do i=1,num_regions
       sigma = regions(i)%sigma
       do iat=1,regions(i)%ntarget_atoms
          too_close=.false.
          iatt = regions(i)%target_atoms(iat)
          if(regions(i)%region_type == BOX)then
             if( sigma/dabs(regions(i)%xmin-cps(iatt,1))>ratio .or. &
               & sigma/dabs(regions(i)%xmax-cps(iatt,1))>ratio .or. &
               & sigma/dabs(regions(i)%ymin-cps(iatt,2))>ratio .or. &
               & sigma/dabs(regions(i)%ymax-cps(iatt,2))>ratio .or. &
               & sigma/dabs(regions(i)%zmin-cps(iatt,3))>ratio .or. &
               & sigma/dabs(regions(i)%zmax-cps(iatt,3))>ratio ) too_close = .true.
          else if (regions(i)%region_type == CYLINDER)then
             tmpr1 = dsqrt((cps(iatt,regions(i)%i1)-regions(i)%cylx)**2+(cps(iatt,regions(i)%i2)-regions(i)%cyly)**2)
             tmpr2 = sigma/(regions(i)%radius-tmpr1)
             if(tmpr2>ratio) too_close = .true.
          endif
          if(too_close.and.printable) &
         & write(nfout,'(a,i8,a,i5)') ' !** WARNING atom ',iat,' might be too close to region ',i
       enddo
    enddo
  end subroutine check_region

  subroutine m_IS_regions(nfout)
    integer, intent(in) :: nfout
    integer :: i,j,iat,iatt
    real(kind=DP) :: tmpf,tmpr1,tmpr2,tmpr6,tmpeps,tmps,tmpf1,tmpf2
    integer, allocatable, dimension(:) :: target_atms
    integer :: n_target_atms
    allocate(target_atms(natm));target_atms=0
    call check_region(nfout)
    do i=1,num_regions
       regions(i)%forc = 0.d0
       regions(i)%energy = 0.d0
       tmpeps = regions(i)%epsi
       tmps   = regions(i)%sigma**12
       target_atms=0
       n_target_atms=0
       do iatt=1,regions(i)%ntarget_atoms
          iat = regions(i)%target_atoms(iatt)
          if(regions(i)%region_type == BOX .and. &
          &  regions(i)%xmin<cps(iat,1) .and. regions(i)%xmax>cps(iat,1) .and. &
          &  regions(i)%ymin<cps(iat,2) .and. regions(i)%ymax>cps(iat,2) .and. &
          &  regions(i)%zmin<cps(iat,3) .and. regions(i)%zmax>cps(iat,3)) then
             n_target_atms = n_target_atms+1
             target_atms(n_target_atms) = iat
          endif
          if(regions(i)%region_type == CYLINDER .and. &
          &  regions(i)%cylzmin<cps(iat,regions(i)%i3) .and. regions(i)%cylzmax>cps(iat,regions(i)%i3) .and. &
          &  in_cylinder(cps(iat,:),regions(i))) then
             n_target_atms = n_target_atms+1
             target_atms(n_target_atms) = iat
          endif
       enddo
       write(nfout,'(a,i5)') ' !** number of target atoms : ',n_target_atms
       do j=1,n_target_atms
          iat = target_atms(j)
          if(regions(i)%region_type == BOX) then
             if(cps(iat,1)>regions(i)%xmin .and. cps(iat,1)<regions(i)%xmax) then
                 tmpr1 = cps(iat,1)-regions(i)%xmin
                 tmpr2 = regions(i)%xmax-cps(iat,1)
                 tmpf = -12.d0*tmpeps*tmps*(-(1.d0/tmpr1)**13 + (1.d0/tmpr2)**13)
                 regions(i)%forc(iat,1) = regions(i)%forc(iat,1)+tmpf
                 regions(i)%energy = regions(i)%energy + tmpeps*tmps*(1.d0/tmpr1**12+1.d0/tmpr2**12)
             endif
             if(cps(iat,2)>regions(i)%ymin .and. cps(iat,2)<regions(i)%ymax) then
                 tmpr1 = cps(iat,2)-regions(i)%ymin
                 tmpr2 = regions(i)%ymax-cps(iat,2)
                 tmpf = -12.d0*tmpeps*tmps*(-(1.d0/tmpr1)**13 + (1.d0/tmpr2)**13)
                 regions(i)%forc(iat,2) = regions(i)%forc(iat,2)+tmpf
                 regions(i)%energy = regions(i)%energy + tmpeps*tmps*(1.d0/tmpr1**12+1.d0/tmpr2**12)
             endif
             if(cps(iat,3)>regions(i)%zmin .and. cps(iat,3)<regions(i)%zmax) then
                 tmpr1 = cps(iat,3)-regions(i)%zmin
                 tmpr2 = regions(i)%zmax-cps(iat,3)
                 tmpf = -12.d0*tmpeps*tmps*(-(1.d0/tmpr1)**13 + (1.d0/tmpr2)**13)
                 regions(i)%forc(iat,3) = regions(i)%forc(iat,3)+tmpf
                 regions(i)%energy = regions(i)%energy + tmpeps*tmps*(1.d0/tmpr1**12+1.d0/tmpr2**12)
             endif
          endif
          if(regions(i)%region_type == CYLINDER) then
             tmpr1 = dsqrt((cps(iat,regions(i)%i1)-regions(i)%cylx)**2+(cps(iat,regions(i)%i2)-regions(i)%cyly)**2)
             tmpr2 = regions(i)%radius-tmpr1
             tmpf1 = -(12.d0*tmpeps*tmps/tmpr2**13) * (cps(iat,regions(i)%i1)-regions(i)%cylx)/tmpr1
             tmpf2 = -(12.d0*tmpeps*tmps/tmpr2**13) * (cps(iat,regions(i)%i2)-regions(i)%cyly)/tmpr1
             regions(i)%forc(iat,regions(i)%i1) = regions(i)%forc(iat,regions(i)%i1) + tmpf1
             regions(i)%forc(iat,regions(i)%i2) = regions(i)%forc(iat,regions(i)%i2) + tmpf2
             regions(i)%energy = regions(i)%energy + tmpeps*tmps*(1.d0/tmpr2**12)
          endif
       enddo
    enddo
    if(printable.and.ipriforce>=2)then
       write(nfout,'(" -- force from regions --")')
       do i=1,num_regions
          do iatt = 1, regions(i)%ntarget_atoms
             iat = regions(i)%target_atoms(iatt)
             write(nfout,'(i5,3f16.8)') iat,(regions(i)%forc(iat,j),j=1,3)
          end do
       enddo
    endif
    if(printable)then
       do i=1,num_regions
          write(nfout,'(a,i5,a,f20.10)') ' energy from region ',i,' : ',regions(i)%energy
       enddo
    endif
  end subroutine m_IS_regions

  logical function in_cylinder(cps,region)
    real(kind=DP), dimension(3), intent(in) :: cps
    type(region_t), intent(in) :: region
    real(kind=DP) :: r2,rad2
    r2 = (region%cylx-cps(region%i1))**2+(region%cyly-cps(region%i2))**2
    rad2 = region%radius**2
    in_cylinder = r2<rad2
  end function in_cylinder


  subroutine m_IS_update_cps_history()
    if(sw_extrapolate_charge==ON .or. sw_wf_predictor==ON) then
      cps_history(:,:,3) = cps_history(:,:,2)
      cps_history(:,:,2) = cps_history(:,:,1)
      cps_history(:,:,1) = cps(:,:)
      ncps_history = ncps_history+1
    endif
  end subroutine m_IS_update_cps_history

  subroutine m_IS_reset_extrpl_status()
    ncps_history = 0
    if(allocated(cps_history)) cps_history = 0.d0
  end subroutine m_IS_reset_extrpl_status

  logical function m_IS_is_extrpl_ready()
    m_IS_is_extrpl_ready = ncps_history>=3
  end function

  subroutine m_IS_get_extpl_factor(alpha,beta,rms,nextpl)
    real(kind=DP), intent(out) :: alpha,beta,rms
    integer, intent(out) :: nextpl
    integer :: iatm,ic,i
    real(kind=DP) :: p1,p2,p3,p4,p5,p6
    real(kind=DP),allocatable,dimension(:,:) :: cps_predicted
    alpha = 0.d0;beta = 0.d0;rms=-1.d0
    nextpl = ncps_history
    if(ncps_history<2) return
    allocate(cps_predicted(natm,3));cps_predicted=0.d0
    p1=0.d0;p2=0.d0;p3=0.d0;p4=0.d0;p5=0.d0
    do ic=1,3
       do iatm=1,natm
          p1 = p1 + (cps_history(iatm,ic,1)-cps_history(iatm,ic,2))**2
          p2 = p2 + (cps_history(iatm,ic,1)-cps_history(iatm,ic,2))*(cps_history(iatm,ic,2)-cps_history(iatm,ic,3))
          p3 = p3 + (cps_history(iatm,ic,2)-cps_history(iatm,ic,3))**2
          p4 = p4 + (cps_history(iatm,ic,1)-cps(iatm,ic))*(cps_history(iatm,ic,1)-cps_history(iatm,ic,2))
          p5 = p5 + (cps_history(iatm,ic,1)-cps(iatm,ic))*(cps_history(iatm,ic,2)-cps_history(iatm,ic,3))
       enddo
    enddo
    p6 = p2*p2-p1*p3
    if(dabs(p6)>1.d-12.and.ncps_history>2)then
       alpha = (p4*p3-p5*p2)/p6
       beta  = (p5*p1-p4*p2)/p6
    else if (dabs(p1)>1.d-12) then
       alpha = p4/p1
       beta = 0.d0
    endif

    if(printable.and.ipripredictor>=2) write(nfout,'(a,f15.10,a,f15.10)') &
      &  ' extrapolation factor : alpha = ',alpha,' beta = ',beta
    do ic=1,3
       do iatm=1,natm
          cps_predicted(iatm,ic) = cps_history(iatm,ic,1) &
       & + alpha*(cps_history(iatm,ic,1)-cps_history(iatm,ic,2)) &
       & + beta* (cps_history(iatm,ic,2)-cps_history(iatm,ic,3))
       enddo
    enddo

    rms = 0.d0
    do iatm=1,natm
       do ic=1,3
          rms = rms+(cps(iatm,ic)-cps_predicted(iatm,ic))**2
       enddo
       if(printable.and.ipripredictor>=2)then
          write(nfout,'(a,6f15.10)') 'cps, cps_predicted ',cps(iatm,1:3),cps_predicted(iatm,1:3)
       endif
    enddo
    rms = dsqrt(rms)/dble(natm)
    if(printable.and.ipripredictor>=2) write(nfout,'(a,f15.10)') 'average RMS : ',rms
    if(rms>rms_threshold) then
       if(printable.and.ipripredictor>=2) &
     & write(nfout,'(a,f13.10)') &
     & '!** WARN rms of the predicted coordinates is greater than the threshold : ',rms_threshold
    endif
    deallocate(cps_predicted)
  end subroutine m_IS_get_extpl_factor

  subroutine alloc_supercell_symmetry()
    if(.not.allocated(napt_supercell)) allocate(napt_supercell(natm,nopr_supercell))
    if(.not.allocated(iop_supercell)) allocate(iop_supercell(nopr_supercell))
!!$    allocate(tau_supercell(3,nopr_supercell))
!!$    allocate(nope_supercell(nopr))
!!$    allocate(pope_supercell(mnope_supercell,nopr))
  end subroutine alloc_supercell_symmetry

  subroutine dealloc_supercell_symmetry()
    if(allocated(napt_supercell)) deallocate(napt_supercell)
    if(allocated(iop_supercell)) deallocate(iop_supercell)
!!$    deallocate(tau_supercell)
!!$    deallocate(nope_supercell)
!!$    deallocate(pope_supercell)
  end subroutine dealloc_supercell_symmetry

  subroutine m_IS_gnrt_supercell_symmetry(paramset,nfout)
! =========== coded by T. Yamasaki after a provided code by Usami-san on Nov. 2013. April 2014 ===
    logical, intent(in) :: paramset
    integer, intent(in) :: nfout

    integer, allocatable, dimension(:,:) :: napt_local
    integer, allocatable, dimension(:)   :: iop_local
    real(kind=DP),allocatable,dimension(:,:) :: tau_local
    integer :: dim_supercell,dim2e
    integer, allocatable, dimension(:) ::   nope_local
    integer, allocatable, dimension(:,:) :: pope_local
    integer, parameter :: modeTAU = 1
    integer, parameter :: modeAPT = 2
    integer :: ii,ia
    integer ::             id_sname = -1

    if(sw_supercell_symmetry == OFF) return
    call tstatc0_begin('m_IS_gnrt_supercell_symmetry ',id_sname)

    dim_supercell = natm*nopr*2
    dim2e = natm*2
    allocate(napt_local(natm,dim_supercell))
    allocate(iop_local(dim_supercell))
!    allocate(tau_local(3,dim_supercell))
!    allocate(nope_local(nopr))
!    allocate(pope_local(dim2e,nopr))
    !call gnrt_supercell_symm_operations(natm2,natm,natm,napt,nopr,nopr,op &
    !     & ,nopr+af,tau,ngen_tl,tau_tl,napt_tl &
    !     & ,lattice_system_from_m_CS_SG,modeAPT &
    !     & ,dim_supercell, napt_local, iop_local, tau_local, nopr_supercell &
    !     & ,dim2e, nope_local, pope_local,iwei)
    if(mype==0) &
    & call gnrt_supercell_symm_operations_apt(natm2,natm,natm,napt,nopr,nopr,op &
         & ,nopr+af,tau,ngen_tl,tau_tl,napt_tl &
         & ,lattice_system_from_m_CS_SG,modeAPT &
         & ,dim_supercell, napt_local, iop_local, nopr_supercell &
         & ,dim2e, iwei)
    !                    -(b_Ionic_System) --> nopr_supercell
    call mpi_bcast(nopr_supercell,1,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(napt_local,natm*nopr_supercell,mpi_integer,0,MPI_CommGroup,ierr)
    call mpi_bcast(iop_local,nopr_supercell,mpi_integer,0,MPI_CommGroup,ierr)

    if(.not.paramset) then
!!$       mnope_supercell = maxval(nope_local(1:nopr))
       call alloc_supercell_symmetry()
       napt_supercell(:,1:nopr_supercell) = napt_local(:,1:nopr_supercell)
       iop_supercell(1:nopr_supercell) = iop_local(1:nopr_supercell)
!       if(mype==0) then
!         do ii=1,nopr_supercell
!            write(0,*) 'iop ',iop_supercell(ii)
!            do ia=1,natm
!               write(0,'(a,3i5)') 'napt ',ii,ia,napt_supercell(ia,ii)
!            enddo
!         enddo
!       endif
!!$       tau_supercell(:,1:nopr_supercell) = tau_local(:,1:nopr_supercell)
!!$       nope_supercell(1:nopr) = nope_local(1:nopr)
!!$       pope_supercell(1:mnope_supercell,1:nopr) = pope_local(1:mnope_supercell,1:nopr)
    end if

    if(ipriinputfile>=2) then
       if(nopr_supercell > nopr) then
          call wd_nopr_supercell_etc1()
          call wd_nopr_supercell_etc2()
       end if
    end if

!    deallocate(tau_local)
    deallocate(iop_local)
    deallocate(napt_local)
!    deallocate(nope_local)
!    deallocate(pope_local)


    call tstatc0_end(id_sname)

    contains
      subroutine wd_nopr_supercell_etc1()
        integer :: iop, i, ic, j, k, iop_temp, ia

        write(nfout,*) ' --- Supercell Symmetry Operations ---'
        write(nfout,'(" !! nopr_supercell = ",i8)') nopr_supercell
       if(.not.paramset) then
        do iop = 1, nopr
!!$           write(nfout,'(" #symmetry op. = ",i8," #elements in this op. = ",i8)')  iop, nope_local(iop)
           do i = 1, 3
              write(nfout,'(3f8.4)') (op(i,k,iop),k=1,3)
           end do
           ic = 0
           do j = 1, nopr_supercell
              if(iop_local(j) == iop) then
                 ic = ic+1
              end if
           end do
           write(nfout,'(" # operation",i3," = ",i8)') iop,ic
        end do
        write(nfout,'(" --- iop_supercell ---")')
        write(nfout,'(24i3)') (iop_supercell(j),j=1,nopr_supercell)
       end if

!!$        if(ipri>=2) then
!!$         if(.not.paramset) then
!!$           do i = 1, nopr
!!$              write(nfout,'(" -- tau_local for #op = ",i8)') i
!!$              do iop = 1, nope_local(i)
!!$                 iop_temp = pope_local(iop,i)
!!$                 write(nfout,'(i3,x,": ",3f10.5)') iop,(tau_local(1:3,iop_temp))
!!$              end do
!!$           end do
!!$        end if
!!$           write(nfout,'(" nopr_supercell = ",i8)') nopr_supercell
!!$           write(nfout,'(" -- iop_supercell, tau_local -- ")')
!!$           do i = 1, nopr_supercell
!!$              write(nfout,'(i8,2x,i8,2x,3f8.4)') i, iop_supercell(i),tau_local(1:3,i)
!!$           end do
!!$        end if
      end subroutine wd_nopr_supercell_etc1
!-----------------------------------------------------------------------------------
      subroutine wd_nopr_supercell_etc2()
        integer :: iop, i, ic, j, k, iop_temp, ia

        write(nfout,'(" -- napt_supercell --")')
        do iop = 1, nopr_supercell
           write(nfout,'(i4,x,": ",64i3)') iop,(napt_supercell(i,iop),i=1,natm)
        end do

        if(ipri>=2) then
           do i = 1, nopr
              do iop = 1, nope_local(i)
                 iop_temp = pope_local(iop,i)
                 write(nfout,'(i3," : ",32i3)') iop, (napt_supercell(ia,iop_temp),ia=1,min(32,natm))
              end do
           end do
        end if
      end subroutine wd_nopr_supercell_etc2
  end subroutine m_IS_gnrt_supercell_symmetry

!***************************** Note ************************************
! This program checks the distances of the atoms.
! The distances between atoms in the neighbor cells and those in the
! main unit cell will also be checked. When the program finds the pair
! which distance are unnaturally small, ID Numbers of those atoms will
! be indicated.
! The density of the atoms in the unit cell will also be checked.
! The unit cell will be divided into several small blocks and the
! density will be check in each of them.
!
!                                   Written by Youky Ono in 2013/Jan
!***********************************************************************
  subroutine check_coor

  Implicit none
!+++++++++++++++++++++++++++++ VARIABLES +++++++++++++++++++++++++++++++
  Integer i,j,k,l,m,num,exnum
  Real(8) aM(3,3),exaM(3,3),exbM(3,3),blaM(3,3),blbM(3,3),detr,dist,min,min1
!  Parameter (min=0.1)
  Real(8), Allocatable  ::  &
&     abc(:,:),xyz(:,:),exxyz(:,:),blxyz(:,:,:,:,:)
  Integer, Allocatable :: bl_atom_num(:,:,:),exID(:),bl_atom_ID(:,:,:,:)

  Real(8) w,a,b,c,x,y,z
  Parameter (w=4.0)

  Real(8) abs_a,abs_b,abs_c,abs_cross,p,pa,pb,pc
  Real(8) dnM(3),d0M(3),daM(3),dbM(3),dcM(3),dM(3)

  Integer nbl_a,nbl_b,nbl_c,cbl_a,cbl_b,cbl_c
  Real(8) bw_a,bw_b,bw_c,bp,vol_bl,dens_bl,dens_max,pi,au,bond
  Parameter (bp=10.0,bond=0.9)

! 'min' is a minimum limit of bond length. [in Bohr]
! 'w' is a width of the margin of the unit cell. [in Bohr]
! 'bp' is a width of the small cell devided from the unit cell. [in Bohr]
!+++++++++++++++++++++++++++ end VARIABLES +++++++++++++++++++++++++++++



!------------------------- Calculation Start ---------------------------
  au = BOHR
  pi = PAI
  min = bond_length_min
  min1 = 1.0d0

!  i = 0
!  open(16,FILE='coor',status='old')
!  Do
!    Read(16,'()',END=999)
!    i = i + 1
!  Enddo
!  999 close(16)
!  num = i - 3


!  Open(111,FILE='coor')
!  Allocate(abc(3,num))
!  Allocate(xyz(3,num))
!  Read(111,*) (aM(1,i),i=1,3)
!  Read(111,*) (aM(2,i),i=1,3)
!  Read(111,*) (aM(3,i),i=1,3)
!  Do i=1,num
!     Read(111,*) (abc(j,i),j=1,3)
!  End do
!  Close(111)

  num = natm

  Allocate(abc(3,num))
  Allocate(xyz(3,num))

  aM(1,:) = altv(:,1)
  aM(2,:) = altv(:,2)
  aM(3,:) = altv(:,3)

  do i=1,natm
     do j=1,3
        abc(j,i) = pos(i,j)
     enddo
  enddo

! STEP 1  :  Roll up all abc into a single unit cell.
  Do i=1,num
    Do j=1,3
      x = abc(j,i)
      abc(j,i) = (x-AINT(x)+1.0)-AINT(x-AINT(x)+1.0)
    Enddo
    a = abc(1,i)
    b = abc(2,i)
    c = abc(3,i)
    xyz(1,i) = (a*aM(1,1) + b*aM(2,1) + c*aM(3,1))
    xyz(2,i) = (a*aM(1,2) + b*aM(2,2) + c*aM(3,2))
    xyz(3,i) = (a*aM(1,3) + b*aM(2,3) + c*aM(3,3))
  Enddo


! STEP 2  :  Make the lattice vectors of the expanded unit cell.
  abs_a = SQRT(aM(1,1)**2.0 + aM(1,2)**2.0 + aM(1,3)**2.0)
  abs_b = SQRT(aM(2,1)**2.0 + aM(2,2)**2.0 + aM(2,3)**2.0)
  abs_c = SQRT(aM(3,1)**2.0 + aM(3,2)**2.0 + aM(3,3)**2.0)

  dnM(1) = aM(1,1)/abs_a + aM(2,1)/abs_b + aM(3,1)/abs_c
  dnM(2) = aM(1,2)/abs_a + aM(2,2)/abs_b + aM(3,2)/abs_c
  dnM(3) = aM(1,3)/abs_a + aM(2,3)/abs_b + aM(3,3)/abs_c

  abs_cross = SQRT((aM(1,2)*dnM(3) - aM(1,3)*dnM(2))**2.0 +  &
&                  (aM(1,3)*dnM(1) - aM(1,1)*dnM(3))**2.0 +  &
&                  (aM(1,1)*dnM(2) - aM(1,2)*dnM(1))**2.0 )
  pa = w*abs_a/abs_cross

  abs_cross = SQRT((aM(2,2)*dnM(3) - aM(2,3)*dnM(2))**2.0 +  &
&                  (aM(2,3)*dnM(1) - aM(2,1)*dnM(3))**2.0 +  &
&                  (aM(2,1)*dnM(2) - aM(2,2)*dnM(1))**2.0 )
  pb = w*abs_b/abs_cross

  abs_cross = SQRT((aM(3,2)*dnM(3) - aM(3,3)*dnM(2))**2.0 +  &
&                  (aM(3,3)*dnM(1) - aM(3,1)*dnM(3))**2.0 +  &
&                  (aM(3,1)*dnM(2) - aM(3,2)*dnM(1))**2.0 )
  pc = w*abs_c/abs_cross

  p = MAX(pa,pb,pc)

  d0M(1) = p*dnM(1)
  d0M(2) = p*dnM(2)
  d0M(3) = p*dnM(3)

  dnM(1) =  aM(1,1)/abs_a - aM(2,1)/abs_b - aM(3,1)/abs_c
  dnM(2) =  aM(1,2)/abs_a - aM(2,2)/abs_b - aM(3,2)/abs_c
  dnM(3) =  aM(1,3)/abs_a - aM(2,3)/abs_b - aM(3,3)/abs_c
  daM(1) =  p*dnM(1)
  daM(2) =  p*dnM(2)
  daM(3) =  p*dnM(3)

  dnM(1) = -aM(1,1)/abs_a + aM(2,1)/abs_b - aM(3,1)/abs_c
  dnM(2) = -aM(1,2)/abs_a + aM(2,2)/abs_b - aM(3,2)/abs_c
  dnM(3) = -aM(1,3)/abs_a + aM(2,3)/abs_b - aM(3,3)/abs_c
  dbM(1) =  p*dnM(1)
  dbM(2) =  p*dnM(2)
  dbM(3) =  p*dnM(3)

  dnM(1) = -aM(1,1)/abs_a - aM(2,1)/abs_b + aM(3,1)/abs_c
  dnM(2) = -aM(1,2)/abs_a - aM(2,2)/abs_b + aM(3,2)/abs_c
  dnM(3) = -aM(1,3)/abs_a - aM(2,3)/abs_b + aM(3,3)/abs_c
  dcM(1) =  p*dnM(1)
  dcM(2) =  p*dnM(2)
  dcM(3) =  p*dnM(3)

  exaM(1,1) =  d0M(1) + aM(1,1) + daM(1)
  exaM(1,2) =  d0M(2) + aM(1,2) + daM(2)
  exaM(1,3) =  d0M(3) + aM(1,3) + daM(3)
  exaM(2,1) =  d0M(1) + aM(2,1) + dbM(1)
  exaM(2,2) =  d0M(2) + aM(2,2) + dbM(2)
  exaM(2,3) =  d0M(3) + aM(2,3) + dbM(3)
  exaM(3,1) =  d0M(1) + aM(3,1) + dcM(1)
  exaM(3,2) =  d0M(2) + aM(3,2) + dcM(2)
  exaM(3,3) =  d0M(3) + aM(3,3) + dcM(3)


! STEP 3 : Expand the unit cell and paste extra atoms in the expanded unit cell.
  detr = (exaM(1,1)*exaM(2,2)*exaM(3,3)+exaM(1,2)*exaM(2,3)*exaM(3,1)+exaM(1,3)*exaM(2,1)*exaM(3,2)) &
&      - (exaM(1,1)*exaM(2,3)*exaM(3,2)+exaM(1,2)*exaM(2,1)*exaM(3,3)+exaM(1,3)*exaM(2,2)*exaM(3,1))

  exbM(1,1) =  (exaM(2,2)*exaM(3,3)-exaM(2,3)*exaM(3,2))/detr
  exbM(2,1) = -(exaM(2,1)*exaM(3,3)-exaM(2,3)*exaM(3,1))/detr
  exbM(3,1) =  (exaM(2,1)*exaM(3,2)-exaM(2,2)*exaM(3,1))/detr
  exbM(1,2) = -(exaM(1,2)*exaM(3,3)-exaM(1,3)*exaM(3,2))/detr
  exbM(2,2) =  (exaM(1,1)*exaM(3,3)-exaM(1,3)*exaM(3,1))/detr
  exbM(3,2) = -(exaM(1,1)*exaM(3,2)-exaM(1,2)*exaM(3,1))/detr
  exbM(1,3) =  (exaM(1,2)*exaM(2,3)-exaM(1,3)*exaM(2,2))/detr
  exbM(2,3) = -(exaM(1,1)*exaM(2,3)-exaM(1,3)*exaM(2,1))/detr
  exbM(3,3) =  (exaM(1,1)*exaM(2,2)-exaM(1,2)*exaM(2,1))/detr

  m = 0
  Do i = 1,num
    Do j = -1,1
    Do k = -1,1
    Do l = -1,1
      a = abc(1,i) + Real(j)
      b = abc(2,i) + Real(k)
      c = abc(3,i) + Real(l)
      x = d0M(1) + (a*aM(1,1) + b*aM(2,1) + c*aM(3,1))
      y = d0M(2) + (a*aM(1,2) + b*aM(2,2) + c*aM(3,2))
      z = d0M(3) + (a*aM(1,3) + b*aM(2,3) + c*aM(3,3))

      a = x*exbM(1,1) + y*exbM(2,1) + z*exbM(3,1)
      b = x*exbM(1,2) + y*exbM(2,2) + z*exbM(3,2)
      c = x*exbM(1,3) + y*exbM(2,3) + z*exbM(3,3)

      If(a.GE.0.AND.a.LE.1.AND.      &
&        b.GE.0.AND.b.LE.1.AND.      &
&        c.GE.0.AND.c.LE.1) Then
        m = m + 1
      Endif

    Enddo
    Enddo
    Enddo
  Enddo
  exnum = m

  Allocate(exxyz(3,exnum))
  Allocate(exID(exnum))

  m = 0
  Do i = 1,num
    Do j = -1,1
    Do k = -1,1
    Do l = -1,1
      a = abc(1,i) + Real(j)
      b = abc(2,i) + Real(k)
      c = abc(3,i) + Real(l)
      x = d0M(1) + (a*aM(1,1) + b*aM(2,1) + c*aM(3,1))
      y = d0M(2) + (a*aM(1,2) + b*aM(2,2) + c*aM(3,2))
      z = d0M(3) + (a*aM(1,3) + b*aM(2,3) + c*aM(3,3))

      a = x*exbM(1,1) + y*exbM(2,1) + z*exbM(3,1)
      b = x*exbM(1,2) + y*exbM(2,2) + z*exbM(3,2)
      c = x*exbM(1,3) + y*exbM(2,3) + z*exbM(3,3)

      If(a.GE.0.AND.a.LE.1.AND.      &
&        b.GE.0.AND.b.LE.1.AND.      &
&        c.GE.0.AND.c.LE.1) Then
        m = m + 1

        exxyz(1,m) = x - d0M(1)
        exxyz(2,m) = y - d0M(2)
        exxyz(3,m) = z - d0M(3)
        exID(m) = i
      Endif

    Enddo
    Enddo
    Enddo
  Enddo


! STEP 4 : Divide the extended unit cell into small blocks.
  nbl_a = 1 + AINT(abs_a/bp)
   bw_a = abs_a/Real(nbl_a)
  nbl_b = 1 + AINT(abs_b/bp)
   bw_b = abs_b/Real(nbl_b)
  nbl_c = 1 + AINT(abs_c/bp)
   bw_c = abs_c/Real(nbl_c)

  Allocate(blxyz(nbl_a,nbl_b,nbl_c,3,exnum))
  Allocate(bl_atom_num(nbl_a,nbl_b,nbl_c))
  Allocate(bl_atom_ID(nbl_a,nbl_b,nbl_c,exnum))
  blxyz = 0

  blaM(1,1) = d0M(1) + aM(1,1)/Real(nbl_a) + daM(1)
  blaM(1,2) = d0M(2) + aM(1,2)/Real(nbl_a) + daM(2)
  blaM(1,3) = d0M(3) + aM(1,3)/Real(nbl_a) + daM(3)
  blaM(2,1) = d0M(1) + aM(2,1)/Real(nbl_b) + dbM(1)
  blaM(2,2) = d0M(2) + aM(2,2)/Real(nbl_b) + dbM(2)
  blaM(2,3) = d0M(3) + aM(2,3)/Real(nbl_b) + dbM(3)
  blaM(3,1) = d0M(1) + aM(3,1)/Real(nbl_c) + dcM(1)
  blaM(3,2) = d0M(2) + aM(3,2)/Real(nbl_c) + dcM(2)
  blaM(3,3) = d0M(3) + aM(3,3)/Real(nbl_c) + dcM(3)

  dens_max = 1.0 / ((4.0/3.0) * pi * (bond/au)**3)
  vol_bl = blaM(1,1)*(blaM(2,2)*blaM(3,3)-blaM(2,3)*blaM(3,2)) &
&        + blaM(1,2)*(blaM(2,3)*blaM(3,1)-blaM(2,1)*blaM(3,3)) &
&        + blaM(1,3)*(blaM(2,1)*blaM(3,2)-blaM(2,2)*blaM(3,1))

  detr = (blaM(1,1)*blaM(2,2)*blaM(3,3)+blaM(1,2)*blaM(2,3)*blaM(3,1)+blaM(1,3)*blaM(2,1)*blaM(3,2)) &
&      - (blaM(1,1)*blaM(2,3)*blaM(3,2)+blaM(1,2)*blaM(2,1)*blaM(3,3)+blaM(1,3)*blaM(2,2)*blaM(3,1))

  blbM(1,1) =  (blaM(2,2)*blaM(3,3)-blaM(2,3)*blaM(3,2))/detr
  blbM(2,1) = -(blaM(2,1)*blaM(3,3)-blaM(2,3)*blaM(3,1))/detr
  blbM(3,1) =  (blaM(2,1)*blaM(3,2)-blaM(2,2)*blaM(3,1))/detr
  blbM(1,2) = -(blaM(1,2)*blaM(3,3)-blaM(1,3)*blaM(3,2))/detr
  blbM(2,2) =  (blaM(1,1)*blaM(3,3)-blaM(1,3)*blaM(3,1))/detr
  blbM(3,2) = -(blaM(1,1)*blaM(3,2)-blaM(1,2)*blaM(3,1))/detr
  blbM(1,3) =  (blaM(1,2)*blaM(2,3)-blaM(1,3)*blaM(2,2))/detr
  blbM(2,3) = -(blaM(1,1)*blaM(2,3)-blaM(1,3)*blaM(2,1))/detr
  blbM(3,3) =  (blaM(1,1)*blaM(2,2)-blaM(1,2)*blaM(2,1))/detr

  Do cbl_a = 1,nbl_a
  Do cbl_b = 1,nbl_b
  Do cbl_c = 1,nbl_c

    dM(1) = aM(1,1)*Real(cbl_a - 1)/Real(nbl_a) + &
&           aM(2,1)*Real(cbl_b - 1)/Real(nbl_b) + &
&           aM(3,1)*Real(cbl_c - 1)/Real(nbl_c) - d0M(1)
    dM(2) = aM(1,2)*Real(cbl_a - 1)/Real(nbl_a) + &
&           aM(2,2)*Real(cbl_b - 1)/Real(nbl_b) + &
&           aM(3,2)*Real(cbl_c - 1)/Real(nbl_c) - d0M(2)
    dM(3) = aM(1,3)*Real(cbl_a - 1)/Real(nbl_a) + &
&           aM(2,3)*Real(cbl_b - 1)/Real(nbl_b) + &
&           aM(3,3)*Real(cbl_c - 1)/Real(nbl_c) - d0M(3)

    m = 0
    Do i = 1,exnum
      x = exxyz(1,i) - dM(1)
      y = exxyz(2,i) - dM(2)
      z = exxyz(3,i) - dM(3)
      a = x*blbM(1,1) + y*blbM(2,1) + z*blbM(3,1)
      b = x*blbM(1,2) + y*blbM(2,2) + z*blbM(3,2)
      c = x*blbM(1,3) + y*blbM(2,3) + z*blbM(3,3)
      If(a.GE.0.AND.a.LE.1.AND.      &
&        b.GE.0.AND.b.LE.1.AND.      &
&        c.GE.0.AND.c.LE.1) Then
        m = m + 1
        blxyz(cbl_a,cbl_b,cbl_c,1,m) = x + dM(1)
        blxyz(cbl_a,cbl_b,cbl_c,2,m) = y + dM(2)
        blxyz(cbl_a,cbl_b,cbl_c,3,m) = z + dM(3)
        bl_atom_ID(cbl_a,cbl_b,cbl_c,m) = exID(i)
      Endif
    Enddo
    bl_atom_num(cbl_a,cbl_b,cbl_c) = m

    dens_bl = m / vol_bl

    If(dens_bl.GT.dens_max) Then
      if(printable) write(nfout,'(a,f20.10)') '### Warning(6101): Density is too large : ',dens_bl
    Endif

  Enddo
  Enddo
  Enddo


! STEP 5  : Check each distance between xyz(:,:) and blxyz(:,:,:,:,:)
  Do i = 1,num
    j = AINT(abc(1,i)*Real(nbl_a)) + 1
    k = AINT(abc(2,i)*Real(nbl_b)) + 1
    l = AINT(abc(3,i)*Real(nbl_c)) + 1

    Do m = 1,bl_atom_num(j,k,l)
       dist = SQRT(                                   &
&             (xyz(1,i) - (blxyz(j,k,l,1,m)) )**2.0 + &
&             (xyz(2,i) - (blxyz(j,k,l,2,m)) )**2.0 + &
&             (xyz(3,i) - (blxyz(j,k,l,3,m)) )**2.0   &
&                 )

      If(i.LT.bl_atom_ID(j,k,l,m).AND.dist.LT.min) Then
         if(printable) write(nfout,12) i,bl_atom_ID(j,k,l,m),dist
         call phase_error_with_msg(nfout, 'Distance between atoms is too small',__LINE__,__FILE__)
      Else if (i.LT.bl_atom_ID(j,k,l,m).and.dist.LT.min1) Then
         if(printable) write(nfout,13) i,bl_atom_ID(j,k,l,m),dist
         call phase_error_with_msg(nfout, 'Distance between atoms is too small',__LINE__,__FILE__)
      Endif
    Enddo
  Enddo

  deallocate(abc,xyz,exxyz,exID,blxyz,bl_atom_num,bl_atom_ID)

 12 Format('### ERROR(6102): Distance between the atom No.' ,I0 , ' and No.' ,I0 , ' is too small : ',f15.10,' bohr')
 13 Format('### WARNING(6102): Distance between the atom No.' ,I0 , ' and No.' ,I0 , ' is too small : ',f15.10,' bohr')


  end subroutine check_coor

  subroutine m_IS_dump_cif(nfdynm_cif,app)
     integer, intent(in) :: nfdynm_cif
     character(len=*), intent(in) :: app
     integer :: iat
     real(kind=DP),allocatable,dimension(:,:) :: cps_full,pos_full
     integer, allocatable,dimension(:) :: ityp_full
     real(kind=DP), dimension(3,3) :: rltv_t
     write(nfdynm_cif,'(a)')        "data_"//trim(adjustl(app))
     write(nfdynm_cif,'(a)')        "_symmetry_cell_setting triclinic"
     write(nfdynm_cif,'(a)')        "_symmetry_space_group_name_H-M  'P1'"
     write(nfdynm_cif,'(a)')        "_symmetry_Int_Tables_number  1"
     write(nfdynm_cif,'(a,f20.10)') "_cell_length_a   ",latconst_len(1)*BOHR
     write(nfdynm_cif,'(a,f20.10)') "_cell_length_b   ",latconst_len(2)*BOHR
     write(nfdynm_cif,'(a,f20.10)') "_cell_length_c   ",latconst_len(3)*BOHR
     write(nfdynm_cif,'(a,f20.10)') "_cell_angle_alpha",latconst_angle(1)
     write(nfdynm_cif,'(a,f20.10)') "_cell_angle_beta ",latconst_angle(2)
     write(nfdynm_cif,'(a,f20.10)') "_cell_angle_gamma",latconst_angle(3)
     write(nfdynm_cif,'(a,f20.10)') "_cell_volume     ",univol*BOHR*BOHR*BOHR
     write(nfdynm_cif,'(a)')        "loop_"
     write(nfdynm_cif,'(a)')        "    _atom_site_label"
     write(nfdynm_cif,'(a)')        "    _atom_site_fract_x"
     write(nfdynm_cif,'(a)')        "    _atom_site_fract_y"
     write(nfdynm_cif,'(a)')        "    _atom_site_fract_z"
     write(nfdynm_cif,'(a)')        "    _atom_site_type_symbol"
     write(nfdynm_cif,'(a)')        "    _atom_site_occupancy"

     allocate(pos_full(natm2,3))
     allocate(ityp_full(natm2))
     if(natm2>natm)then
        allocate(cps_full(natm2,3))
        call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
        rltv_t = transpose(rltv)/PAI2
        call change_of_coordinate_system(rltv_t,cps_full,natm2,natm2,pos_full)
     else
        pos_full(1:natm,1:3) = pos(1:natm,1:3)
        ityp_full(1:natm)  = ityp(1:natm)
     endif
     do iat=1,natm2
        write(nfdynm_cif,'(a5,3f20.10,a5,i3)') &
       & speciesname(ityp_full(iat)),pos_full(iat,1),pos_full(iat,2),pos_full(iat,3),speciesname(ityp_full(iat)),1
     enddo
     if(natm2>natm)then
        deallocate(cps_full)
     endif
     deallocate(ityp_full)
     deallocate(pos_full)
  end subroutine m_IS_dump_cif

  subroutine m_IS_dump_nnpfiles(energy, forc_l)
     real(kind=DP), intent(in) :: energy
     real(kind=DP), dimension(natm,3), intent(in) :: forc_l
     if(filetype_nnp == N2P2 .or. filetype_nnp == ALL) then
        call m_IS_dump_n2p2(nfn2p2, energy, forc_l)
     endif
     if(filetype_nnp == XSF .or. filetype_nnp == ALL) then
        call m_IS_dump_xsf(nfxsf, energy, forc_l)
     endif
     if(filetype_nnp == DEEPMD .or. filetype_nnp == ALL) then
        call m_IS_dump_deepmd(energy, forc_l)
     endif
  end subroutine m_IS_dump_nnpfiles

  subroutine m_IS_dump_xsf(nfxsf,energy,forc_l)
     integer, intent(in) :: nfxsf
     real(kind=DP), intent(in) :: energy
     real(kind=DP), dimension(natm,3), intent(in) :: forc_l
     real(kind=DP) :: factor
     integer :: ia
     call m_Files_open_xsf()
     if(mype == 0) then
       write(nfout,'(a)') ' !** output coordinates in XSF format '
       factor = HARTREE/BOHR
       write(nfxsf,'(a,f25.10,a)') '# total energy = ',energy*HARTREE,' eV'
       write(nfxsf,*)
       write(nfxsf,'(a)')          'CRYSTAL'
       write(nfxsf,'(a)')          'PRIMVEC'
       write(nfxsf,'(3f25.10)')     altv(1:3,1)*BOHR
       write(nfxsf,'(3f25.10)')     altv(1:3,2)*BOHR
       write(nfxsf,'(3f25.10)')     altv(1:3,3)*BOHR
       write(nfxsf,'(a)')          'PRIMCOORD'
       write(nfxsf,'(2i8)')        natm2,1
       do ia=1, natm
          write(nfxsf,'(a5,6f20.10)') &
         & speciesname(ityp(ia)),cps(ia,1)*BOHR,cps(ia,2)*BOHR,cps(ia,3)*BOHR &
         & , forc_l(ia,1)*factor, forc_l(ia,2)*factor, forc_l(ia,3)*factor
       enddo
     endif
     call m_Files_close_xsf
  end subroutine m_IS_dump_xsf

  subroutine m_IS_dump_n2p2(nfn2p2,energy,forc_l)
     integer, intent(in) :: nfn2p2
     real(kind=DP), intent(in) :: energy
     real(kind=DP), dimension(natm,3), intent(in) :: forc_l
     real(kind=DP) :: factor
     integer :: ia
     call m_Files_open_n2p2()
     factor = HARTREE/BOHR
     if(mype == 0) then
       write(nfout,'(a)') ' !** output coordinates in N2P2 format '
       write(nfn2p2,'(a)')         'begin'
       write(nfn2p2,'(a,i8)')      'comment iteration ',iteration_ionic
       write(nfn2p2,'(a,3f25.10)') 'lattice', altv(1:3,1)*BOHR
       write(nfn2p2,'(a,3f25.10)') 'lattice', altv(1:3,2)*BOHR
       write(nfn2p2,'(a,3f25.10)') 'lattice', altv(1:3,3)*BOHR
       do ia=1, natm
          write(nfn2p2,'(a,3f20.10,a,5f20.10)') &
         &                        'atom ',cps(ia,1)*BOHR,cps(ia,2)*BOHR,cps(ia,3)*BOHR  &
         &                               ,' '//speciesname(ityp(ia))//' ', 0.0, 0.0 &
         &                               ,forc_l(ia,1)*factor &
         &                               ,forc_l(ia,2)*factor &
         &                               ,forc_l(ia,3)*factor
       enddo
       write(nfn2p2,'(a,f20.10)')  'energy ',energy*HARTREE
       write(nfn2p2,'(a)')         'charge 0.0'
       write(nfn2p2,'(a)')         'end'
     endif
     call m_Files_close_n2p2()
  end subroutine m_IS_dump_n2p2

  subroutine m_IS_dump_deepmd(energy, forc_l)
    real(kind=DP), intent(in) :: energy
    real(kind=DP), dimension(natm,3), intent(in) :: forc_l
    integer :: nfdeep, i, j
    logical :: ext
    real(kind=DP), allocatable, dimension(:) :: ar
    character(len=256) :: ch

    if(mype==0) then
      write(nfout,'(a)') ' !** output coordinates in deepmd format '
      nfdeep = get_unused_unitnumber()
      inquire(file='set.000',exist=ext)
      if(.not. ext) call system('mkdir -p set.000')
      inquire(file='type.raw',exist=ext)
      if(.not.ext .or. icond==INITIAL) then
        open(nfdeep,file='type.raw')
        do i=1,natm
          write(nfdeep,*) ityp(i)-1
        enddo
        close(nfdeep)
      endif
      inquire(file='type_map.raw',exist=ext)
      if(.not.ext .or. icond==INITIAL) then
        open(nfdeep,file='type_map.raw')
        do i=1,ntyp
          write(nfdeep,*) trim(speciesname(i))
        enddo
        close(nfdeep)
      endif
      call m_Files_open_file(nfdeep, 'set.000/box.raw ', 'F_DEEP_BOX')
      write(nfdeep,'(9f20.10)') altv(1,1)*BOHR,altv(2,1)*BOHR,altv(3,1)*BOHR, &
                                altv(1,2)*BOHR,altv(2,2)*BOHR,altv(3,2)*BOHR, &
                                altv(1,3)*BOHR,altv(2,3)*BOHR,altv(3,3)*BOHR
      call m_Files_close_file(nfdeep)

      call m_Files_open_file(nfdeep, 'set.000/energy.raw ', 'F_DEEP_ENERGY')
      write(nfdeep,'(f25.15)') energy*HARTREE
      call m_Files_close_file(nfdeep)

      call m_Files_open_file(nfdeep, 'set.000/coord.raw ', 'F_DEEP_COORD')
      allocate(ar(3*natm))
      do i=1,natm
        do j=1,3
          ar((i-1)*3+j) = cps(i,j)*BOHR
        enddo
      enddo
      write(ch,*) natm*3
      write(nfdeep,'('//trim(adjustl(ch))//'f25.15)') ar(1:natm*3)
      call m_Files_close_file(nfdeep)

      call m_Files_open_file(nfdeep, 'set.000/force.raw ', 'F_DEEP_FORCE')
      do i=1,natm
        do j=1,3
          ar((i-1)*3+j) = forc_l(i,j)*HARTREE/BOHR
        enddo
      enddo
      write(ch,*) natm*3
      write(nfdeep, '('//trim(adjustl(ch))//'f25.15)') ar(1:natm*3)
      call m_Files_close_file(nfdeep)
      deallocate(ar)
    endif

    contains
    integer function get_unused_unitnumber()
    integer :: i
    logical :: op
    do i=500,9999
      inquire(unit=i,opened=op)
      if( .not.op ) then
        get_unused_unitnumber = i
        return
        endif
    enddo
    call phase_error_with_msg(nfout,'  insufficient file handle',__LINE__,__FILE__)
  end function get_unused_unitnumber

  end subroutine m_IS_dump_deepmd

  logical function absolute_convergence_of_forc(forc_l_in)
    real(kind=DP), intent(in), dimension(natm,3) :: forc_l_in
    real(kind=DP),allocatable,dimension(:,:) :: forct
    integer :: ia
    allocate(forct(natm,3));forct=forc_l_in
    absolute_convergence_of_forc = .false.
    do ia=1,natm
       if(imdtyp(ia)==0) forct(ia,:)=0.d0
    enddo
    if(sum(forct(1:natm,1:3)**2).lt.1e-15)then
        absolute_convergence_of_forc = .true.
    endif
    deallocate(forct)
  end function absolute_convergence_of_forc

  subroutine m_IS_get_tauinv_from_qmass(tauinv)
    real(kind=DP), dimension(nrsv) :: tauinv
    integer :: ia, ir, i, ibath
    call md2_alloc()
    nathm = 0
    do ia = 1, natm
      ir = ibath(imdtyp(ia))    !-(b_Ionic_System)
      if(ir>=1) then
        nathm(ir) = nathm(ir) + 3 * iwei(ia)
      endif
    enddo
    do i = 1, nrsv
      tauinv(i) = 1.d0/(PAI2*sqrt(qmass(i)/dble(nathm(i)*2.d0*tkb(i))))
    enddo
    call md2_dealloc()
    return

  end subroutine m_IS_get_tauinv_from_qmass

  subroutine m_IS_add_friction_force(forc_l)
    real(kind=DP), dimension(natm,3), intent(inout) :: forc_l
    real(kind=DP), allocatable, dimension(:) :: tauinv
    integer :: ia, j, ir, ibath
    allocate(tauinv(nrsv))
    call m_IS_get_tauinv_from_qmass(tauinv)
    do ia=1,natm
      ir = ibath(imdtyp(ia))    !-(b_Ionic_System)
      if(ir>=1) then
        do j=1,3
          forc_l(ia,j) = forc_l(ia,j) - tauinv(ir) * amion(ityp(ia)) * cpd_l(ia,j)
        enddo
      endif
    enddo
    deallocate(tauinv)
  end subroutine m_IS_add_friction_force

  subroutine m_IS_add_random_force(forc_l)
    real(kind=DP), dimension(natm,3), intent(inout) :: forc_l
    integer :: i, ia, j, ir, ibath
    real(kind=DP) :: idamp, idt, fac, z1, z2
    real(kind=DP), allocatable, dimension(:) :: tauinv
    real(kind=DP), allocatable, dimension(:) :: rand

    allocate(tauinv(nrsv))
    call m_IS_get_tauinv_from_qmass(tauinv)

    allocate(rand(natm*3))

    do ia = 1, natm*3, 2
      call normal_random_number(z1, z2)
      rand(ia)   = z1
      if(ia+1<=natm*3) rand(ia+1) = z2
    enddo

    idt = 1.d0/dtio
    do ia = 1, natm
      ir = ibath(imdtyp(ia))    !-(b_Ionic_System)
      if(ir>=1) then
        fac = sqrt(2.d0 * amion(ityp(ia)) * idt * tauinv(ir) * tkb(ir))
        do j = 1, 3
          forc_l(ia,j) = forc_l(ia,j) + fac * rand((ia-1)*3+j)
        enddo
      endif
    enddo

    deallocate(tauinv)
    deallocate(rand)
  end subroutine m_IS_add_random_force

#ifdef __EDA__
! -----  ascat starts modifying  -----
  subroutine m_IS_alloc_eewald_per_atom(natm)

    integer,intent(in) :: natm

    allocate(eewald_per_atom(natm))
  end subroutine m_IS_alloc_eewald_per_atom

  subroutine m_IS_dealloc_eewald_per_atom
    deallocate(eewald_per_atom)
  end subroutine m_IS_dealloc_eewald_per_atom

  subroutine m_IS_alloc_zfm3_EDA(natm)

    integer,intent(in) :: natm

    allocate(zfm3_l_EDA(ista_kngp:iend_kngp,natm,kimg)); zfm3_l_EDA = 0.d0
  end subroutine m_IS_alloc_zfm3_EDA

  subroutine m_IS_dealloc_zfm3_EDA
    deallocate(zfm3_l_EDA)
  end subroutine m_IS_dealloc_zfm3_EDA
! -----  ascat ceases modifying  -----
#endif

  logical function m_IS_rigid_body_exists()
    m_IS_rigid_body_exists = nrigid_bodies>0
  end function m_IS_rigid_body_exists

  subroutine m_IS_rigid_body_map_force(forc_l)
    real(kind=DP), dimension(natm,3), intent(in) :: forc_l
    integer :: i,j,ia,na
    do i=1,nrigid_bodies
      na = rigid_bodies(i)%natm
      rigid_bodies(i)%force = 0.d0
      do j=1,na
        ia = rigid_bodies(i)%atoms(j)
        rigid_bodies(i)%force_per_atm(j,:) = forc_l(ia,:)
        rigid_bodies(i)%force(:) = rigid_bodies(i)%force(:) + forc_l(ia,:)
      enddo
      if(printable .and. iprimd>=2) write(nfout,'(a,i5,3f10.5)') '!** translational force for rb ',i,rigid_bodies(i)%force(1:3)
    enddo
    call update_rb_torque()
    return
  end subroutine m_IS_rigid_body_map_force

  subroutine update_rb_torque()
    integer :: i, j, na
    real(kind=DP), dimension(3) :: cps, forc, vout
    do i=1,nrigid_bodies
      rigid_bodies(i)%torque = 0.d0
      na = rigid_bodies(i)%natm
      do j=1,na
        cps  = rigid_bodies(i)%relative_coords(j,:)
        forc = rigid_bodies(i)%force_per_atm(j,:)
        call cross_product(cps,forc,vout)
        rigid_bodies(i)%torque(:) = rigid_bodies(i)%torque(:) + vout(:)
      enddo
      if(printable .and. iprimd>=2) write(nfout,'(a,i5,3f10.5)') '!** torque for rb ',i,rigid_bodies(i)%torque(1:3)
    enddo
    return
  end subroutine update_rb_torque

  subroutine cross_product(a,b,c)
    implicit none
    real(DP), dimension(3), intent(in) :: a,b
    real(DP), dimension(3), intent(out)::  c
    c=0.d0
    c(1) = a(2)*b(3)-a(3)*b(2)
    c(2) = a(3)*b(1)-a(1)*b(3)
    c(3) = a(1)*b(2)-a(2)*b(1)
  end subroutine cross_product

  subroutine update_Qmat(rigid_body,qhalfway)
    type(rigid_body_t), target, intent(inout) :: rigid_body
    logical, intent(in), optional             :: qhalfway
    real(kind=DP), pointer, dimension(:)      :: q
    real(kind=DP), pointer, dimension(:,:)    :: Qmat
    logical                                   :: qh
    qh = .false.
    if(present(qhalfway)) qh = qhalfway
    if(.not.qh) then
      q  => rigid_body%q
    else
      q  => rigid_body%q_h
    endif
    Qmat => rigid_body%Qmat
    Qmat(1,1) = -q(3)
    Qmat(1,2) = -q(4)
    Qmat(1,3) =  q(2)
    Qmat(1,4) =  q(1)
    Qmat(2,1) =  q(4)
    Qmat(2,2) = -q(3)
    Qmat(2,3) = -q(1)
    Qmat(2,4) =  q(2)
    Qmat(3,1) =  q(1)
    Qmat(3,2) =  q(2)
    Qmat(3,3) =  q(4)
    Qmat(3,4) =  q(3)
    Qmat(4,1) = -q(2)
    Qmat(4,2) =  q(1)
    Qmat(4,3) = -q(3)
    Qmat(4,4) =  q(4)

    Qmat      =  0.5d0*Qmat
    return
  end subroutine update_Qmat

  subroutine update_rotation_matrix(rigid_body, qhalfway)
    type(rigid_body_t), target, intent(inout) :: rigid_body
    logical, intent(in), optional             :: qhalfway
    real(kind=DP), pointer, dimension(:)      :: q
    real(kind=DP), pointer, dimension(:,:)    :: rotmat
    logical                                   :: qh
    real(kind=DP)                             :: phi,theta,psi
    real(kind=DP)                             :: cphi,sphi,cpsi,spsi,ctheta,stheta
    integer                                   :: i
    qh = .false.
    if(present(qhalfway)) qh = qhalfway
    if(.not.qh) then
      q    => rigid_body%q
    else
      q    => rigid_body%q_h
    endif
    rotmat => rigid_body%rotmat

    if(rigid_body%fix_in_a_plane == OFF) then
      rotmat(1,1) = -q(1)**2 + q(2)**2 - q(3)**2 + q(4)**2
      rotmat(2,2) =  q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2
      rotmat(3,3) = -q(1)**2 - q(2)**2 + q(3)**2 + q(4)**2
      rotmat(1,2) =  2.d0*(q(3)*q(4)-q(1)*q(2))
      rotmat(1,3) =  2.d0*(q(2)*q(3)+q(1)*q(4))
      rotmat(2,1) = -2.d0*(q(1)*q(2)+q(3)*q(4))
      rotmat(2,3) =  2.d0*(q(2)*q(4)-q(1)*q(3))
      rotmat(3,1) =  2.d0*(q(2)*q(3)-q(1)*q(4))
      rotmat(3,2) = -2.d0*(q(1)*q(3)+q(2)*q(4))
    else
      if(.not.qh) then
        phi    = rigid_body%phi
        psi    = rigid_body%psi
        theta  = rigid_body%theta
      else
        phi    = rigid_body%phi_h
        psi    = rigid_body%psi_h
        theta  = rigid_body%theta_h
      endif
      cphi   = cos(phi)
      sphi   = sin(phi)
      cpsi   = cos(psi)
      spsi   = sin(psi)
      ctheta = cos(theta)
      stheta = sin(theta)

      rotmat(1,1) = cpsi*cphi-ctheta*sphi*spsi
      rotmat(2,1) = -spsi*cphi-ctheta*sphi*cpsi
      rotmat(3,1) = stheta*sphi
      rotmat(1,2) = cpsi*sphi+ctheta*cphi*spsi
      rotmat(2,2) = -spsi*sphi+ctheta*cphi*cpsi
      rotmat(3,2) = -stheta*cphi
      rotmat(1,3) = spsi*stheta
      rotmat(2,3) = cpsi*stheta
      rotmat(3,3) = ctheta
      rotmat = transpose(rotmat)
      write(nfout,'(a,3f10.5)') 'euler angles',180.d0*phi/PAI,180.d0*theta/PAI,180.d0*psi/PAI
      write(nfout,'(a)')        'rotation matrix'
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(1,i),i=1,3)
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(2,i),i=1,3)
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(3,i),i=1,3)
    endif

    return
  end subroutine update_rotation_matrix

  subroutine normalize_quaternion(rigid_body, qhalfway)
    type(rigid_body_t), target, intent(inout) :: rigid_body
    logical, intent(in), optional             :: qhalfway
    real(kind=DP), pointer, dimension(:)      :: q
    real(kind=DP) :: sm
    integer :: i
    logical :: qh
    qh = .false.
    if(present(qhalfway)) qh = qhalfway
    if(.not.qh) then
      q  => rigid_body%q
    else
      q  => rigid_body%q_h
    endif
    sm =  0.d0
    do i=1,4
      sm = sm+q(i)*q(i)
    enddo
    sm = sqrt(sm)
    q = q/sm
  end subroutine normalize_quaternion

  subroutine read_rigid_body_options()
    integer :: i, iret
    real(kind=DP) :: dret
    character(len=256) :: idstr
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock

    if(f_selectBlock(tag_rigid_body)==0) then
      if(f_getRealValue(tag_dt_rotation,dret,'au_time')==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%dt_rotation = dret
        enddo
      endif
      if(f_getRealValue(tag_dt_translation,dret,'au_time')==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%dt_translation = dret
        enddo
      endif
      if(f_getIntValue(tag_mobile,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%mobile  = iret
        enddo
      endif
      if(f_getIntValue(tag_mobilex,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%mobile(1) = iret
        enddo
      endif
      if(f_getIntValue(tag_mobiley,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%mobile(2) = iret
        enddo
      endif
      if(f_getIntValue(tag_mobilez,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%mobile(3) = iret
        enddo
      endif
      if(f_getIntValue(tag_mobilerot,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%mobilerot = iret
        enddo
      endif
      if(f_getRealValue(tag_normx,dret,'')==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(1) = dret
        enddo
      endif
      if(f_getRealValue(tag_normy,dret,'')==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(2) = dret
        enddo
      endif
      if(f_getRealValue(tag_normz,dret,'')==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(3) = dret
        enddo
      endif
      if(f_getIntValue(tag_fix_in_a_plane,iret)==0) then
        do i=1,nrigid_bodies
          rigid_bodies(i)%fix_in_a_plane = iret
        enddo
      endif
      iret = f_selectParentBlock()
    endif
    do i=1,nrigid_bodies
      write(idstr,*) rigid_bodies(i)%id
      if(f_selectBlock(tag_rigid_body//trim(adjustl(idstr)))==0)then
        if(f_getRealValue(tag_dt_rotation,dret,'au_time')==0) then
          rigid_bodies(i)%dt_rotation = dret
        endif
        if(f_getRealValue(tag_dt_translation,dret,'au_time')==0) then
          rigid_bodies(i)%dt_translation = dret
        endif
        if(f_getIntValue(tag_mobile,iret)==0) then
          rigid_bodies(i)%mobile = iret
        endif
        if(f_getIntValue(tag_mobilex,iret)==0) then
          rigid_bodies(i)%mobile(1) = iret
        endif
        if(f_getIntValue(tag_mobiley,iret)==0) then
          rigid_bodies(i)%mobile(2) = iret
        endif
        if(f_getIntValue(tag_mobilez,iret)==0) then
          rigid_bodies(i)%mobile(3) = iret
        endif
        if(f_getIntValue(tag_mobilerot,iret)==0) then
          rigid_bodies(i)%mobilerot = iret
        endif
        if(f_getIntValue(tag_thermo_group,iret)==0) then
          rigid_bodies(i)%thermo_group = iret
        endif
        if(f_getRealValue(tag_normx,dret,'')==0) then
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(1) = dret
        endif
        if(f_getRealValue(tag_normy,dret,'')==0) then
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(2) = dret
        endif
        if(f_getRealValue(tag_normz,dret,'')==0) then
          rigid_bodies(i)%explicit_zaxis = .true.
          rigid_bodies(i)%zaxis(3) = dret
        endif
        if(f_getIntValue(tag_fix_in_a_plane,iret)==0) then
          rigid_bodies(i)%fix_in_a_plane = iret
        endif
        iret = f_selectParentBlock()
      endif
    enddo
  end subroutine read_rigid_body_options

  subroutine initialize_zaxis()
    integer :: i
    do i=1,nrigid_bodies
      if (.not. rigid_bodies(i)%explicit_zaxis) cycle
      call resolve_xy_axis(rigid_bodies(i))
    enddo

    contains

    subroutine resolve_xy_axis(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP), dimension(3) :: xaxis, yaxis, zaxis, vec
      integer :: j
      vec(:) = rigid_body%zaxis(:)
      if(abs(maxval(vec)) < 1.d-10) then
        if(printable) write(nfout,'(a)') '!** the components of the specified zaxis are too small'
        rigid_body%explicit_zaxis = .false.
        return
      endif
      call normalize(vec,zaxis)
      rigid_body%zaxis = zaxis
      vec(:) = zaxis(:)
      vec(1) = vec(1)+0.123
      vec(2) = vec(2)+0.348
      vec(3) = vec(3)-0.5
      call normalize(vec,xaxis)
      vec(:) = zaxis(:)
      vec(1) = vec(1)-0.5
      vec(2) = vec(2)-0.112
      vec(3) = vec(3)+0.34
      call normalize(vec,yaxis)

      ! GS
      vec(:) = xaxis(:) - dot_product(xaxis, zaxis) * zaxis(:)
      call normalize(vec, xaxis)

      rigid_body%xaxis = xaxis

      vec(:) = yaxis(:) - dot_product(yaxis, zaxis) * zaxis(:) - dot_product(yaxis, xaxis) * xaxis(:)
      call normalize(vec, yaxis)
      rigid_body%yaxis = yaxis

      if(printable) then
        write(nfout,'(a,i5)')     '!** x-axis, y-axis and z-axis for the rotation axis of rigid body ',rigid_body%id
        write(nfout,'(a,3f10.5)') '!** x-axis ',(rigid_body%xaxis(j),j=1,3)
        write(nfout,'(a,3f10.5)') '!** y-axis ',(rigid_body%yaxis(j),j=1,3)
        write(nfout,'(a,3f10.5)') '!** z-axis ',(rigid_body%zaxis(j),j=1,3)
        write(nfout,'(a,f10.5)')  '!** x dot y',dot_product(rigid_body%xaxis,rigid_body%yaxis)
        write(nfout,'(a,f10.5)')  '!** x dot z',dot_product(rigid_body%xaxis,rigid_body%zaxis)
        write(nfout,'(a,f10.5)')  '!** y dot z',dot_product(rigid_body%yaxis,rigid_body%zaxis)
      endif
    end subroutine resolve_xy_axis

    subroutine normalize(a,b)
      real(kind=DP), dimension(3), intent(in)  :: a
      real(kind=DP), dimension(3), intent(out) :: b
      real(kind=DP)                            :: sm
      sm   = sqrt(dot_product(a,a))
      b(:) = a(:)/sm
      return
    end subroutine normalize

  end subroutine initialize_zaxis

  subroutine read_rigid_body_input2()
    integer :: iret
    real(kind=DP) :: dret
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
    if( f_selectBlock(tag_atoms) == 0) then
      call set_rigid_bodies()
      iret = f_selectParentBlock()
    endif

    return

    contains

    subroutine set_rigid_bodies()
      integer :: nrbcount
      integer :: i,j,k,ind
      integer, allocatable, dimension(:) :: rbbuf, icount, jcount
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
      integer :: iicount

      i = 1
      nrbcount = 0
      allocate(rbbuf(natm));rbbuf=0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_id, iret)==0) then
          if(iret /= 0) then
            if(new_rigid_body(iret, rbbuf)) then
              nrbcount = nrbcount+1
              rbbuf(i) = iret
            endif
          endif
        endif
        i = i+1
      enddo
      nrigid_bodies = nrbcount
      if(nrbcount==0) then
        return
      endif
      allocate(rigid_bodies(nrigid_bodies))

      do j=1, nrigid_bodies
        rigid_bodies(j)%id   = 0
        rigid_bodies(j)%natm = 0
      enddo
      i = 1
      nrbcount = 0
      iicount = 0
      allocate(jcount(natm));jcount=0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_id,iret)==0) then
          iicount = iicount+1
          jcount(iicount) = iret
          if(iret /= 0) then
            call find_rigid_body(iret, nrbcount, ind)
            if (ind>0) then
              rigid_bodies(ind)%natm      = rigid_bodies(ind)%natm+1
            else
              nrbcount                    = nrbcount+1
              rigid_bodies(nrbcount)%id   = iret
              rigid_bodies(nrbcount)%natm = rigid_bodies(nrbcount)%natm+1
            endif
          endif
        endif
        i = i+1
      enddo

      do j=1, nrigid_bodies
        allocate(rigid_bodies(j)%atoms(rigid_bodies(j)%natm))
      enddo

      i = 1
      nrbcount = 0
      iicount = 0
      allocate(icount(nrigid_bodies));icount=0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_no,iret)==0) then
          if(iret /= 0) then
            iicount = iicount+1
            !call find_rigid_body(iret, nrigid_bodies, ind)
            ind = jcount(iicount)
            icount(ind) = icount(ind)+1
            rigid_bodies(ind)%atoms(icount(ind)) = iret
          endif
        endif
        i = i+1
      enddo

      allocate(is_rigid_body(natm))
      is_rigid_body = .false.
      do i=1,natm
        jloop:do j=1,nrigid_bodies
          do k=1,rigid_bodies(j)%natm
            if(rigid_bodies(j)%atoms(k)==i) then
              is_rigid_body(i) = .true.
              imdtyp(i) = OFF
              imdtypxyz(i,:) = OFF
              exit jloop
            endif
          enddo
        enddo jloop
      enddo
      deallocate(rbbuf)
      deallocate(icount)
      return
    end subroutine set_rigid_bodies

    logical function new_rigid_body(icand, rbbuf)
      integer, intent(in)      :: icand
      integer, dimension(natm) :: rbbuf
      integer :: i
      do i=1,natm
        if(icand==rbbuf(i))then
          new_rigid_body = .false.
          return
        endif
      enddo
      new_rigid_body = .true.
      return
    end function new_rigid_body

    subroutine find_rigid_body(id,nrbcount,ind)
      integer, intent(in)  :: id
      integer, intent(in)  :: nrbcount
      integer, intent(out) :: ind
      integer :: i
      do i=1,nrbcount
        if(rigid_bodies(i)%id==id) then
          ind = i
          return
        endif
      enddo
      ind = -1
      return
    end subroutine find_rigid_body

  end subroutine read_rigid_body_input2

  subroutine read_rigid_body_input()
    integer :: iret
    real(kind=DP) :: dret
    integer :: f_selectTop, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
    if( f_selectBlock(tag_atoms) == 0) then
      call set_rigid_bodies()
!      if(nrigid_bodies>0) then
!        call initialize_rigid_bodies()
!        call print_rigid_body_status()
!      endif
      iret = f_selectParentBlock()
    endif

    return

    contains

    subroutine set_rigid_bodies()
      integer :: nrbcount
      integer :: i,j,k,ind
      integer, allocatable, dimension(:) :: rbbuf, icount
      integer :: f_selectFirstTableLine, f_selectNextTableLine
      integer :: f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue, f_selectParentBlock
      i = 1
      nrbcount = 0
      allocate(rbbuf(natm));rbbuf=0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_rigid_body, iret)==0) then
          if(iret /= 0) then
            if(new_rigid_body(iret, rbbuf)) then
              nrbcount = nrbcount+1
              rbbuf(i) = iret
            endif
          endif
        endif
        i = i+1
      enddo
      nrigid_bodies = nrbcount
      if(nrbcount==0) then
        return
      endif
      allocate(rigid_bodies(nrigid_bodies))

      do j=1, nrigid_bodies
        rigid_bodies(j)%id   = 0
        rigid_bodies(j)%natm = 0
      enddo
      i = 1
      nrbcount = 0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_rigid_body,iret)==0) then
          if(iret /= 0) then
            call find_rigid_body(iret, nrbcount, ind)
            if (ind>0) then
              rigid_bodies(ind)%natm      = rigid_bodies(ind)%natm+1
            else
              nrbcount                    = nrbcount+1
              rigid_bodies(nrbcount)%id   = iret
              rigid_bodies(nrbcount)%natm = rigid_bodies(nrbcount)%natm+1
            endif
          endif
        endif
        i = i+1
      enddo

      do j=1, nrigid_bodies
        allocate(rigid_bodies(j)%atoms(rigid_bodies(j)%natm))
      enddo

      i = 1
      nrbcount = 0
      allocate(icount(nrigid_bodies));icount=0
      do while(.true.)
        if (i == 1) then
          if(f_selectFirstTableLine() /= 0) then
            exit
          end if
        else
          if(f_selectNextTableLine() /= 0) then
            exit
          end if
        end if
        if(f_getIntValue(tag_rigid_body,iret)==0) then
          if(iret /= 0) then
            call find_rigid_body(iret, nrigid_bodies, ind)
            icount(ind) = icount(ind)+1
            rigid_bodies(ind)%atoms(icount(ind)) = i
          endif
        endif
        i = i+1
      enddo

      allocate(is_rigid_body(natm))
      is_rigid_body = .false.
      do i=1,natm
        jloop:do j=1,nrigid_bodies
          do k=1,rigid_bodies(j)%natm
            if(rigid_bodies(j)%atoms(k)==i) then
              is_rigid_body(i) = .true.
              imdtyp(i) = OFF
              imdtypxyz(i,:) = OFF
              exit jloop
            endif
          enddo
        enddo jloop
      enddo
      deallocate(rbbuf)
      deallocate(icount)
      return
    end subroutine set_rigid_bodies

    logical function new_rigid_body(icand, rbbuf)
      integer, intent(in)      :: icand
      integer, dimension(natm) :: rbbuf
      integer :: i
      do i=1,natm
        if(icand==rbbuf(i))then
          new_rigid_body = .false.
          return
        endif
      enddo
      new_rigid_body = .true.
      return
    end function new_rigid_body

    subroutine find_rigid_body(id,nrbcount,ind)
      integer, intent(in)  :: id
      integer, intent(in)  :: nrbcount
      integer, intent(out) :: ind
      integer :: i
      do i=1,nrbcount
        if(rigid_bodies(i)%id==id) then
          ind = i
          return
        endif
      enddo
      ind = -1
      return
    end subroutine find_rigid_body

  end subroutine read_rigid_body_input

  subroutine initialize_rigid_bodies()
    integer :: i,j, it, ia, na
    real(kind=DP),dimension(3) :: sm
    do i=1,nrigid_bodies
      rigid_bodies(i)%mass = 0.d0
      rigid_bodies(i)%COM  = 0.d0
      rigid_bodies(i)%angular_momentum   = 0.d0
      rigid_bodies(i)%angular_momentum_h = 0.d0
      sm = 0.d0
      do j=1,rigid_bodies(i)%natm
        sm(:) = sm(:)+cpd_l(rigid_bodies(i)%atoms(j),:)
      enddo
      rigid_bodies(i)%velocity(1:3)  = sm(1:3)
      rigid_bodies(i)%velocity_h(1:3)= 0.d0
      rigid_bodies(i)%dt_translation = dtio
      rigid_bodies(i)%dt_rotation    = dtio
      rigid_bodies(i)%mobile         = ON
      rigid_bodies(i)%mobilerot      = ON
      rigid_bodies(i)%thermo_group   = 1
      na = rigid_bodies(i)%natm
      do j=1,na
        ia = rigid_bodies(i)%atoms(j)
        it = ityp(ia)
        rigid_bodies(i)%mass   = rigid_bodies(i)%mass   + amion(it)
        rigid_bodies(i)%COM(:) = rigid_bodies(i)%COM(:) + amion(it)*cps(ia,:)
      enddo
      rigid_bodies(i)%COM = rigid_bodies(i)%COM/rigid_bodies(i)%mass
      rigid_bodies(i)%explicit_zaxis = .false.
      rigid_bodies(i)%zaxis = 0.d0
      rigid_bodies(i)%fix_in_a_plane = OFF
      rigid_bodies(i)%constrain_on_plane = .false.
      rigid_bodies(i)%move_plane = .false.
      allocate(rigid_bodies(i)%relative_coords(na,3))
      allocate(rigid_bodies(i)%molecular_coords(na,3))
      allocate(rigid_bodies(i)%force_per_atm(na,3))
      call resolve_plane(rigid_bodies(i))
    enddo
    return

    contains

    subroutine resolve_plane(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      integer :: ic
      if (allocated(icnst_a)) ic = icnst_a(rigid_body%atoms(1))
      if (allocated(icount_of_ipfixedplane)) then
        ic = icnst_a(rigid_body%atoms(1))
        rigid_body%constrain_on_plane = icount_of_ipfixedplane(ic)>0
        rigid_body%move_plane         = rigid_body%constrain_on_plane .and. move_constrained_plane
        if(rigid_body%constrain_on_plane) then
          rigid_body%normx = fcvect(ic,1)
          rigid_body%normy = fcvect(ic,2)
          rigid_body%normz = fcvect(ic,3)
        endif
        if(rigid_body%move_plane) then
           rigid_body%dinc   = fcvect(ic,4)
           rigid_body%incmax = fcvect(ic,5)
           rigid_body%incx   = fcvect(ic,6)
           rigid_body%incy   = fcvect(ic,7)
           rigid_body%incz   = fcvect(ic,8)
           rigid_body%mobile = relax_in_fixedplane(ic)
           rigid_body%distance_moved_thusfar = 0.d0
        endif
      endif
    end subroutine resolve_plane

  end subroutine initialize_rigid_bodies

  subroutine initialize_rotation_matrix()
    call principle_axis_of_inertia()
    call build_quaternion_from_rotmat()
    call map_cps_to_molecular_coords()
    call map_cps_to_relative_coords()
    return
  end subroutine initialize_rotation_matrix

  subroutine print_rigid_body_status()
    use m_Files, only : nfout
    use m_Control_Parameters, only : printable
    real(kind=DP) :: f
    integer :: i,j
    if(printable)then
      write(nfout,'(a,i5)')        ' !** number of rigid bodies defined : ',nrigid_bodies
      f = 180.d0/PAI
      do i=1,nrigid_bodies
        write(nfout,'(a,i5,a,i5)') ' !** rigid body no. ',i,' ID ',rigid_bodies(i)%id
        write(nfout,'(a,i5)')      ' !** number of atoms defined in this rigid body ',rigid_bodies(i)%natm
        write(nfout,'(8i8)')       (rigid_bodies(i)%atoms(j), j=1,rigid_bodies(i)%natm)
        write(nfout,'(a,f20.5)')   ' !** mass         ', rigid_bodies(i)%mass
        write(nfout,'(a,3f20.5)')  ' !** COM          ', (rigid_bodies(i)%COM(j),j=1,3)
        write(nfout,'(a,3f20.5)')  ' !** inertia      ',(rigid_bodies(i)%inertia(j),j=1,3)
        write(nfout,'(a,4f20.5)')  ' !** quaternion   ',(rigid_bodies(i)%q(j),j=1,4)
        write(nfout,'(a,3f20.5)')  ' !** euler angles ',rigid_bodies(i)%phi*f,rigid_bodies(i)%theta*f,rigid_bodies(i)%psi*f
        write(nfout,'(a)')         ' !** initial rotation matrix'
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(1,j),j=1,3)
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(2,j),j=1,3)
        write(nfout,'(3f10.5)')    (rigid_bodies(i)%rotmat(3,j),j=1,3)
        if(rigid_bodies(i)%constrain_on_plane) then
        write(nfout,'(a,3f10.5)')  ' !** norm of the plane on which this rigid body shall be constrained '&
                                   , rigid_bodies(i)%normx,rigid_bodies(i)%normy,rigid_bodies(i)%normz
        if(.not. rigid_bodies(i)%move_plane) then
        write(nfout,'(a)')         ' !** the plane is imobile'
        else
        write(nfout,'(a,3f10.5)')  ' !** the plane will move along direction                             ' &
                                   ,rigid_bodies(i)%incx,rigid_bodies(i)%incy,rigid_bodies(i)%incz
        endif
        endif
      enddo
      write(nfout,'(a)')           ' !** atom belongs to some rigid body '
      write(nfout,'(8l3)')         (is_rigid_body(j),j=1,natm)
    endif
    return
  end subroutine print_rigid_body_status

  subroutine build_quaternion_from_rotmat()
    real(kind=DP), pointer, dimension(:,:) :: rotmat
    integer :: i,j
    real(kind=DP) :: theta, phi, psi, sinphi, cosphi, cospsi, sinpsi, b, tmp, costheta, sintheta
    real(kind=DP) :: x1,x2,x3,x4,a,p1,p2,p3,p4
    real(kind=DP) :: very_small = 1.d-12
    do i=1,nrigid_bodies
      rotmat => rigid_bodies(i)%rotmat
      theta  =  acos(rotmat(3,3))
      costheta = rotmat(3,3)
      sintheta = sin(theta)
      phi    =  0.d0
      psi    =  0.d0
!      if(abs(theta).gt.very_small)then
!        tmp = rotmat(3,2)/dsin(theta)
!        if(tmp.gt.1)tmp=1
!        if(tmp.lt.-1)tmp=-1
!        psi = acos(tmp)
!!        tmp = rotmat(1,3)/sin(theta)
!        if(tmp.gt.1)tmp=1
!        if(tmp.lt.-1)tmp=-1
!        phi = asin(tmp)
!      else
!        theta = 0.d0
!        psi   = 0.d0
!        phi   = dacos(rotmat(1,1))
!      endif

      if (abs(sintheta)>very_small) then
        sinphi =  rotmat(3,1)/sintheta
        cosphi = -rotmat(3,2)/sintheta
!        if (sinphi>0) then
        phi    =  acos(cosphi)
!        else
!        phi    =  PAI2-acos(cosphi)
!        endif
        sinpsi =  rotmat(1,3)/sintheta
        cospsi =  rotmat(2,3)/sintheta
!        if (sinpsi>0) then
        psi    =  acos(cospsi)
!        else
!        psi    =  PAI2-acos(cospsi)
!        endif
        !cospsi=(cosphi*rotmat(1,1)/(costheta*sinphi)-rotmat(2,1))/(cosphi**2/(costheta*sinphi)+costheta*sinphi)
        !sinphi=(cosphi*cospsi-rotmat(1,1))/(costheta*sinphi)
        !if (sinpsi>0) then
        !psi    =  acos(cospsi)
        !else
        !psi    =  PAI2-acos(cospsi)
        !endif
      endif
!      phi   = atan(rotmat(2,3)/rotmat(1,3))
!      theta = atan(sqrt(1.d0-rotmat(3,3)**2)/rotmat(3,3))
!      psi   = atan(rotmat(3,2)/-rotmat(3,1))
!      p1 = rotmat(1,1)
!      p2 = rotmat(2,1)
!      p3 = rotmat(1,2)
!      p4 = rotmat(2,2)
!      a  = cos(theta)
!      x2 = (p4-a*p1)/(1+a*a)
!      x1 = p1+a*x2
!      x4 = (p3+a*p2)/(1-a*a)
!      x3 = -p2-a*x4
!      phi = atan(x3/x1)
!      psi = atan(x2/x4)
!      write(nfout,'(a,2f10.5)') 'phi, psi from a different path ',180.d0*atan(x3/x1)/PAI,180.30*atan(x2/x4)/PAI
      !psi = 0.d0
      !if (abs(cos(theta))>very_small .and. abs(sinphi)>very_small) then
      !  b = cos(theta)*sin(phi)
      !  cospsi = ((cos(phi)*rotmat(1,1))/b-rotmat(2,1))/(cos(phi)**2/b+b)
      !  sinpsi =  (cos(phi)*cospsi-rotmat(1,1))/b
      !  psi    =  sign(1.d0,sinpsi)*acos(cospsi)
      !endif
      !write(nfout,'(a,3f10.5)') 'EULER ',theta,psi,phi
      rigid_bodies(i)%q(1)  = sin(0.5d0*theta)*sin(0.5d0*(psi-phi))
      rigid_bodies(i)%q(2)  = sin(0.5d0*theta)*cos(0.5d0*(psi-phi))
      rigid_bodies(i)%q(3)  = cos(0.5d0*theta)*sin(0.5d0*(psi+phi))
      rigid_bodies(i)%q(4)  = cos(0.5d0*theta)*cos(0.5d0*(psi+phi))
      rigid_bodies(i)%phi   = phi
      rigid_bodies(i)%theta = theta
      rigid_bodies(i)%psi   = psi
      call normalize_quaternion(rigid_bodies(i))
      call update_rotation_matrix(rigid_bodies(i))
    enddo
    return
  end subroutine build_quaternion_from_rotmat

  subroutine set_initial_velocities_rb()
    integer                                  :: i,j
    real(kind=DP)                            :: rn1,rn2
    real(kind=DP), dimension(3)              :: pcom,a
    real(kind=DP), allocatable, dimension(:) :: rand
    real(kind=DP)                            :: mcom, factor
    real(kind=DP),dimension(nrsv)            :: tkin
    integer,      dimension(nrsv)            :: nir
    integer                                  :: irp

    if(nrigid_bodies<2) then
      if(printable) write(nfout,'(a)') '!** initial velocities/angular momentum for rigid bodies can only be set when&
                                       & nrigid_bodies>=2'
      return
    endif

    ! normal random numbers
    allocate(rand(2*nrigid_bodies*3))
    do i=1,2*nrigid_bodies*3, 2
      call normal_random_number(rn1, rn2)
      rand(i)   = rn1
      if(i+1<=2*nrigid_bodies*3) rand(i+1) = rn2
    enddo

    ! assign random velocity/angular momentum
    mcom = 0.d0
    pcom = 0.d0
    a    = 0.d0
    do i=1,nrigid_bodies
      mcom = mcom+rigid_bodies(i)%mass
      do j=1,3
        rigid_bodies(i)%velocity(j)         = rand((i-1)*3+j)
        rigid_bodies(i)%angular_momentum(j) = rand(nrigid_bodies*3+(i-1)*3+j)
        pcom(j) = pcom(j) + rigid_bodies(i)%velocity(j)*rigid_bodies(i)%mass
        a(j)    = a(j)    + rigid_bodies(i)%angular_momentum(j)
      enddo
    enddo

    deallocate(rand)

    if(mcom.gt.1e-12) then
      pcom(:) = pcom(:)/mcom
    endif

    ! shift velocity
    do i=1,nrigid_bodies
      rigid_bodies(i)%velocity(:)           = rigid_bodies(i)%velocity(:)         - pcom(:)
      rigid_bodies(i)%angular_momentum(:)   = rigid_bodies(i)%angular_momentum(:) - a(:)
      rigid_bodies(i)%velocity_old(:)       = rigid_bodies(i)%velocity(:)
      rigid_bodies(i)%angular_momentum_h(:) = rigid_bodies(i)%angular_momentum(:)
    enddo

    ! scale
    call scale_velocity_rb()
    return
  end subroutine set_initial_velocities_rb

  subroutine map_cps_to_molecular_coords()
    integer :: i, j, na, ia
    do i=1,nrigid_bodies
      na = rigid_bodies(i)%natm
      do j=1,na
        ia = rigid_bodies(i)%atoms(j)
        rigid_bodies(i)%molecular_coords(j,:) = matmul(rigid_bodies(i)%rotmat(:,:),cps(ia,:)-rigid_bodies(i)%COM(:))
      enddo
    enddo
    return
  end subroutine map_cps_to_molecular_coords

  subroutine map_cps_to_relative_coords()
    integer :: i, j, na, ia
    do i=1,nrigid_bodies
      na = rigid_bodies(i)%natm
      do j=1,na
        ia = rigid_bodies(i)%atoms(j)
        rigid_bodies(i)%relative_coords(j,:) = cps(ia,:)-rigid_bodies(i)%COM(:)
      enddo
    enddo
    return
  end subroutine map_cps_to_relative_coords

  subroutine map_molecular_coords_to_cps()
    integer :: i,j,ia,na
    real(kind=DP), dimension(3,3) :: rotinv
    do i=1,nrigid_bodies
      na = rigid_bodies(i)%natm
      rotinv = transpose(rigid_bodies(i)%rotmat)
      do j=1,na
        ia        = rigid_bodies(i)%atoms(j)
        cps(ia,:) = rigid_bodies(i)%COM(:)+matmul(rotinv(:,:), rigid_bodies(i)%molecular_coords(j,:))
      enddo
    enddo

    call m_IS_cps_to_pos()
    return
  end subroutine map_molecular_coords_to_cps

  logical function is_molecular_dynamics()
    integer :: mdalg
    mdalg = m_CtrlP_what_is_mdalg()
    is_molecular_dynamics = mdalg == VERLET     .or. mdalg == T_CONTROL &
    &                  .or. mdalg == PT_CONTROL .or. mdalg == P_CONTROL
    return
  end function is_molecular_dynamics

  subroutine cal_max_force_and_torque()
    torque_max = get_rb_max_torque()
    trans_force_max = get_rb_max_trans()
  end subroutine cal_max_force_and_torque

  logical function m_IS_rigid_body_converged()
    real(kind=DP) :: tor, trans
    if(nrigid_bodies==0) then
      m_IS_rigid_body_converged = .true.
      return
    endif
    if (is_molecular_dynamics()) then
      m_IS_rigid_body_converged = .false.
      return
    endif
    torque_max      = get_rb_max_torque()
    trans_force_max = get_rb_max_trans()
    if (torque_max.gt.max_torque .or. trans_force_max.gt.max_force_trans) then
      m_IS_rigid_body_converged = .false.
    else
      m_IS_rigid_body_converged = .true.
    endif
    if(printable) write(nfout,'(a,2f10.5)') '!** max trans. force and torque ',trans_force_max,torque_max
    return
  end function m_IS_rigid_body_converged

  real(kind=DP) function get_rb_max_trans()
    integer :: i,j
    real(kind=DP) :: maxv,v
    maxv = 0.d0
    do i=1,nrigid_bodies
      v = 0.d0
      do j=1,3
        if(rigid_bodies(i)%mobile(j)==ON) v = v+rigid_bodies(i)%force(j)**2
      enddo
      v = sqrt(v)
      if (v>maxv) maxv = v
    enddo
    get_rb_max_trans = maxv
    return
  end function get_rb_max_trans

  real(kind=DP) function get_rb_max_torque()
    integer :: i
    real(kind=DP) :: maxv,v
    maxv = 0.d0
    do i=1,nrigid_bodies

      if(rigid_bodies(i)%fix_in_a_plane==OFF) then
        v = sqrt(dot_product(rigid_bodies(i)%torque,rigid_bodies(i)%torque))
      else
        v = sqrt((rigid_bodies(i)%torque(1)**2+rigid_bodies(i)%torque(2)**2))
      endif
      if(rigid_bodies(i)%mobilerot==OFF) v=0.d0
      if (v>maxv) maxv = v
    enddo
    get_rb_max_torque = maxv
    return
  end function get_rb_max_torque

  subroutine m_IS_rb_dealloc()
    integer :: i
    do i=1,nrigid_bodies
      deallocate(rigid_bodies(i)%atoms)
      deallocate(rigid_bodies(i)%relative_coords)
      deallocate(rigid_bodies(i)%force_per_atm)
    enddo
    deallocate(rigid_bodies)
    return
  end subroutine m_IS_rb_dealloc

  subroutine m_IS_rb_dynamics(forc_l, mdalg)
    real(kind=DP), dimension(natm,3), intent(in) :: forc_l
    integer, intent(in) :: mdalg
    integer :: i
    call md1_alloc(VERLET)
    if(mdalg==FIRE) then
      call rb_trans_fire(forc_l)
    endif
    do i=1, nrigid_bodies
      if (is_molecular_dynamics() .or. mdalg==QUENCHED_MD) then
        call evolve_translational_velocity(rigid_bodies(i))
        call evolve_COM(rigid_bodies(i))
      endif
      if(rigid_bodies(i)%mobilerot==ON) then
        call evolve_angular_momentum_pre(rigid_bodies(i))
        if(rigid_bodies(i)%fix_in_a_plane==OFF) then
          call evolve_quaternion_pre(rigid_bodies(i))
        else
          call evolve_euler_pre(rigid_bodies(i))
        endif
        call evolve_angular_momentum(rigid_bodies(i))

        if(rigid_bodies(i)%fix_in_a_plane==OFF) then
          call evolve_quaternion(rigid_bodies(i))
        else
          call evolve_euler(rigid_bodies(i))
        endif
      endif
    enddo
    call map_molecular_coords_to_cps()
    call map_cps_to_relative_coords()
    call get_ekina()
    call md1_dealloc
    call cal_max_force_and_torque()
    return

    contains

    subroutine evolve_translational_velocity(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      integer :: i
      real(kind=DP) :: factor,denom,amass,pdot
      real(kind=DP), dimension(3) :: dcg, dcg_n, dcgh, dcgh_n, fcvect_t
      rigid_body%velocity_old = rigid_body%velocity
      factor = 1.d0
      if(iteration_ionic==1) factor = 0.5d0
      do i=1,3
        if(rigid_body%mobile(i)/=ON) cycle
        rigid_body%velocity(i)   = rigid_body%velocity(i)    &
      &                          + rigid_body%dt_translation &
      &                          * factor * rigid_body%force(i)/rigid_body%mass
        rigid_body%velocity_h(i) = rigid_body%velocity(i)    &
      &                          + 0.5d0*rigid_body%dt_translation &
      &                          * factor * rigid_body%force(i)/rigid_body%mass
      enddo
      if(printable .and. iprimd>=2) write(nfout,'(a,i5,3f15.10)') '!** translational velocity for rb ', &
      &                             rigid_body%id,rigid_body%velocity_h(1:3)

      if(rigid_body%constrain_on_plane) then
        denom  = 0.d0
        dcg    = 0.d0
        dcg_n  = 0.d0
        dcgh   = 0.d0
        dcgh_n = 0.d0
        fcvect_t(1) = rigid_body%normx
        fcvect_t(2) = rigid_body%normy
        fcvect_t(3) = rigid_body%normz
        do i=1,rigid_body%natm
          amass     = amion(ityp(rigid_body%atoms(i)))
          denom     = denom+amass
          dcg(1:3)  = dcg(1:3)+amass*rigid_body%velocity(1:3)
          dcgh(1:3) = dcgh(1:3)+amass*rigid_body%velocity_h(1:3)
        enddo
        dcg                 = dcg/denom
        pdot                = dot_product(fcvect_t,dcg)
        dcg_n               = pdot*fcvect_t
        rigid_body%velocity = dcg - dcg_n

        dcgh                  = dcgh/denom
        pdot                  = dot_product(fcvect_t,dcgh)
        dcgh_n                = pdot*fcvect_t
        rigid_body%velocity_h = dcgh - dcgh_n
      endif

      if(mdalg==QUENCHED_MD .and. dot_product(rigid_body%velocity,rigid_body%force)<0.d0) then
        rigid_body%velocity = 0.d0
        if(printable) write(nfout,'(a,i5)') '!** quenched translational velocity for rigid body ',rigid_body%id
      endif
    end subroutine evolve_translational_velocity

    subroutine evolve_COM(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP) :: sx,sy,sz
      integer :: i
      logical, save :: reached_final_value=.false.
      do i=1,3
        if(rigid_body%mobile(i)/=ON) cycle
        rigid_body%COM(i) = rigid_body%COM(i) + rigid_body%dt_translation * rigid_body%velocity(i)
      enddo
      if(rigid_body%move_plane) then
        if(rigid_body%distance_moved_thusfar < rigid_body%incmax) then
          sx = rigid_body%incx*rigid_body%dinc
          sy = rigid_body%incy*rigid_body%dinc
          sz = rigid_body%incz*rigid_body%dinc
          rigid_body%COM(1) = rigid_body%COM(1) + sx
          rigid_body%COM(2) = rigid_body%COM(2) + sy
          rigid_body%COM(3) = rigid_body%COM(3) + sz
          rigid_body%distance_moved_thusfar = rigid_body%distance_moved_thusfar+sqrt(sx*sx+sy*sy+sz*sz)
          if(printable .and. iprimd>=2) then
            write(nfout,'(a,i5,a,f20.5,a)') '!** the plane associated to rigid body ',rigid_body%id,' moved a total of '&
                                      , rigid_body%distance_moved_thusfar,' bohr thus far'
          endif
        else
            if(.not.reached_final_value) &
            write(nfout,'(a)')              '!** reached final_value, so the plane will not move any more'
            reached_final_value = .true.
        endif
      endif
    end subroutine evolve_COM

    subroutine evolve_angular_momentum_pre(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP)                     :: dotp
      rigid_body%rotmat_old = rigid_body%rotmat
      rigid_body%angular_momentum_old = rigid_body%angular_momentum
      rigid_body%angular_momentum_h = rigid_body%angular_momentum &
      &                           + 0.5d0*rigid_body%dt_rotation*rigid_body%torque
      if(printable .and. iprimd>=2) write(nfout,'(a,i5,3f25.10)') '!** angular momentum for rb       ', &
      &                             rigid_body%id,rigid_body%angular_momentum_h(1:3)
      if(rigid_body%fix_in_a_plane==OFF) then
        dotp = dot_product(rigid_body%angular_momentum, rigid_body%torque)
        if(mdalg==QUENCHED_MD .and. dotp<0.d0) then
          rigid_body%angular_momentum_h = 0.d0
          if(printable) write(nfout,'(a,i5)') '!** quenched angular momentum for rigid body ',rigid_body%id
        endif
      else
        if(rigid_body%angular_momentum_h(1)*rigid_body%torque(1)<0.d0) rigid_body%angular_momentum_h(1) = 0.d0
        if(rigid_body%angular_momentum_h(2)*rigid_body%torque(2)<0.d0) rigid_body%angular_momentum_h(2) = 0.d0
      endif
    end subroutine evolve_angular_momentum_pre

    subroutine evolve_quaternion_pre(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP), dimension(3) :: romega
      real(kind=DP), dimension(4) :: omega,qomega
      integer :: i,j
      call update_Qmat(rigid_body)
      romega = matmul(rigid_body%rotmat,rigid_body%angular_momentum_h)
      do i=1,3
        omega(i) = romega(i)/rigid_body%inertia(i)
      enddo
      omega(4) = 0.d0
      qomega = matmul(rigid_body%Qmat,omega)
      rigid_body%q_h(:) = rigid_body%q(:) + 0.5d0*rigid_body%dt_rotation*qomega(:)
      call normalize_quaternion(rigid_body,qhalfway=.true.)
      call update_rotation_matrix(rigid_body,qhalfway=.true.)
      call update_Qmat(rigid_body,qhalfway=.true.)
    end subroutine evolve_quaternion_pre

    subroutine evolve_euler_pre(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP), dimension(3) :: der
      integer :: i,j
      call time_derivative_of_euler_angles(rigid_body, .true., der)
      rigid_body%phi_h   = rigid_body%phi   - 0.5d0*rigid_body%dt_rotation*der(1)
      if(rigid_body%fix_in_a_plane==OFF) then
        rigid_body%theta_h = rigid_body%theta - 0.5d0*rigid_body%dt_rotation*der(2)
        rigid_body%psi_h   = rigid_body%psi   - 0.5d0*rigid_body%dt_rotation*der(3)
      endif
      call update_rotation_matrix(rigid_body,.true.)
    end subroutine evolve_euler_pre

    subroutine evolve_angular_momentum(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP)                     :: dotp
      rigid_body%angular_momentum = rigid_body%angular_momentum &
                                  + rigid_body%dt_rotation * rigid_body%torque
      if(rigid_body%fix_in_a_plane==OFF) then
        dotp = dot_product(rigid_body%angular_momentum, rigid_body%torque)
        if(mdalg==QUENCHED_MD .and. dotp<0.d0) then
          rigid_body%angular_momentum = 0.d0
          if(printable) write(nfout,'(a,i5)') '!** quenched angular momentum for rigid body ',rigid_body%id
        endif
      else
        if(rigid_body%angular_momentum(1)*rigid_body%torque(1)<0.d0) rigid_body%angular_momentum(1) = 0.d0
        if(rigid_body%angular_momentum(2)*rigid_body%torque(2)<0.d0) rigid_body%angular_momentum(2) = 0.d0
      endif
    end subroutine evolve_angular_momentum

    subroutine evolve_quaternion(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP), dimension(3) :: romega
      real(kind=DP), dimension(4) :: omega,qomega
      integer :: i,j
      romega = matmul(rigid_body%rotmat,rigid_body%angular_momentum)
      do i=1,3
        omega(i) = romega(i)/rigid_body%inertia(i)
      enddo
      omega(4) = 0.d0
      qomega = matmul(rigid_body%Qmat,omega)
      rigid_body%q(:) = rigid_body%q(:) + rigid_body%dt_rotation*qomega(:)
      call normalize_quaternion(rigid_body)
      call update_rotation_matrix(rigid_body)
    end subroutine evolve_quaternion

    subroutine evolve_euler(rigid_body)
      type(rigid_body_t), intent(inout) :: rigid_body
      real(kind=DP), dimension(3) :: der
      integer :: i,j
      call time_derivative_of_euler_angles(rigid_body, .false., der)
      rigid_body%phi   = rigid_body%phi   - rigid_body%dt_rotation*der(1)
      if(rigid_body%fix_in_a_plane==OFF) then
        rigid_body%theta = rigid_body%theta - rigid_body%dt_rotation*der(2)
        rigid_body%psi   = rigid_body%psi   - rigid_body%dt_rotation*der(3)
      endif
      if(printable) write(nfout,'(a,i5,3f15.10)') '!** der ',rigid_body%id,der(1),der(2),der(3)
      call update_rotation_matrix(rigid_body,.false.)
    end subroutine evolve_euler

    subroutine time_derivative_of_euler_angles(rigid_body, halfway, der)
      type(rigid_body_t),          intent(in)  :: rigid_body
      logical,                     intent(in)  :: halfway
      real(kind=DP), dimension(3), intent(out) :: der
      real(kind=DP), dimension(3) :: omega
      real(kind=DP) :: phi,theta,psi,dphi,dtheta,dpsi
      real(kind=DP) :: cphi,sphi,ctheta,stheta,cpsi,spsi,ttheta
!      if(halfway) then
!        omega  = matmul(rigid_body%rotmat,rigid_body%angular_momentum_h)/rigid_body%inertia
!      else
!        omega  = matmul(rigid_body%rotmat,rigid_body%angular_momentum)/rigid_body%inertia
!      endif
      if(halfway) then
        omega  = rigid_body%angular_momentum_h/rigid_body%inertia
      else
        omega  = rigid_body%angular_momentum/rigid_body%inertia
      endif
      phi    = rigid_body%phi
      theta  = rigid_body%theta
      psi    = rigid_body%psi
      ctheta = cos(theta)
      stheta = sin(theta)
      cpsi   = cos(psi)
      spsi   = sin(psi)
      der(1) = (1.d0/stheta)   * (omega(1)*spsi+omega(2)*cpsi)
      der(2) = omega(1)*cpsi-omega(2)*spsi
      der(3) = omega(3)-ctheta/stheta * (omega(1)*spsi+omega(2)*cpsi)
    end subroutine time_derivative_of_euler_angles


    subroutine rb_trans_fire(forc_l)
      real(kind=DP), dimension(natm,3), intent(in) :: forc_l
      integer, allocatable, dimension(:), save :: iteration_from_lq
      real(kind=DP), allocatable, dimension(:,:) :: forc_l_in
      real(kind=DP), allocatable, dimension(:), save :: fire_alpha
      real(kind=DP), allocatable, dimension(:), save :: fire_dt
      real(kind=DP), allocatable, dimension(:,:), save :: forcold
      real(kind=DP), allocatable, dimension(:,:) :: cpscom,vel
      integer, allocatable, dimension(:,:) :: imdtyp
      if(.not.allocated(fire_alpha))then
         allocate(fire_alpha(nrigid_bodies));fire_alpha = fire_alpha_start
      endif
      if(.not.allocated(fire_dt)) then
         allocate(fire_dt(nrigid_bodies));fire_dt = fire_initial_dt
      endif
      if(.not.allocated(iteration_from_lq)) then
         allocate(iteration_from_lq(nrigid_bodies));iteration_from_lq=0
      endif
      if(.not.allocated(forcold)) then
        allocate(forcold(nrigid_bodies,3))
        do i=1,nrigid_bodies
          forcold(i,:) = rigid_bodies(i)%force(:)
        enddo
      endif
      allocate(forc_l_in(nrigid_bodies,3))
      do i=1,nrigid_bodies
        forc_l_in(i,:) = rigid_bodies(i)%force(:)
      enddo
      allocate(cpscom(nrigid_bodies,3))
      do i=1,nrigid_bodies
        cpscom(i,:) = rigid_bodies(i)%COM(:)
      enddo
      allocate(vel(nrigid_bodies,3))
      do i=1,nrigid_bodies
        vel(i,:) = rigid_bodies(i)%velocity(:)
      enddo
      allocate(imdtyp(nrigid_bodies,3))
      do i=1,nrigid_bodies
        imdtyp(i,:) = rigid_bodies(i)%mobile(:)
      enddo
      call m_IS_fire_param_by_args(nfout,nrigid_bodies,iteration_ionic,cpscom,vel,forc_l_in,imdtyp &
      &                           ,iteration_from_lq, fire_alpha, fire_dt, forcold)
      do i=1,nrigid_bodies
        rigid_bodies(i)%COM(:) = cpscom(i,:)
      enddo
      deallocate(imdtyp)
      deallocate(cpscom)
      deallocate(forc_l_in)
    end subroutine rb_trans_fire

  end subroutine m_IS_rb_dynamics

  real(kind=DP) function m_IS_rb_kinetic_energies()
    integer :: i
    real(kind=DP) :: ret
    ret = 0.d0
    do i=1,nrigid_bodies
      ret = ret + rb_kinetic_energy(rigid_bodies(i))
    enddo
    m_IS_rb_kinetic_energies = ret
    return
  end function m_IS_rb_kinetic_energies

  real(kind=DP) function rb_kinetic_energy(rigid_body)
    type(rigid_body_t), intent(in) :: rigid_body
    integer :: j
    real(kind=DP) :: ret
    real(kind=DP),dimension(3) :: omega !,omegao
    omega  = matmul(rigid_body%rotmat,rigid_body%angular_momentum_h)
    omega  = omega/rigid_body%inertia
!    omegao = matmul(rigid_body%rotmat_old,rigid_body%angular_momentum_old)
!    omegao = omegao/rigid_body%inertia
    ret = 0.d0
    do j=1,3
      ret = ret + rigid_body%inertia(j) * omega(j)**2 &
    &           + rigid_body%mass       * (0.5d0*(rigid_body%velocity_old(j)+rigid_body%velocity(j)))**2
    enddo
    rb_kinetic_energy = 0.5d0*ret
    return
  end function rb_kinetic_energy

  subroutine m_IS_rd_rigid_body(nfcntn)
    integer, intent(in)        :: nfcntn
    integer                    :: nrt,i,j,natmt,ierr
    logical                    :: tag_is_found, EOF_reach
    type(rigid_body_t),pointer :: rigid_body

    if(nrigid_bodies==0) return

    if(mype==0)then
       call rewind_to_tag0(nfcntn,len(tag_rigid_body),tag_rigid_body &
            &, EOF_reach, tag_is_found, str,len_str)
       if(.not.tag_is_found) then
          call phase_error_with_msg(nfout,' tag_rigid_body is not found',__LINE__,__FILE__)
       else
          read(nfcntn,*)
          read(nfcntn,*) nrt
       endif
    endif
    if(npes>1)then
       call mpi_bcast(nrt,1 &
            & ,mpi_integer,0,MPI_CommGroup,ierr)
    endif

    if(nrt /= nrigid_bodies) then
      if(printable) write(nfout,'(a)') 'nrigid_bodies from nfcntn /= nrigid_bodies from nfinp'
      return
    endif

    if (mype==0) then
      do i=1,nrigid_bodies
        rigid_body => rigid_bodies(i)
        read(nfcntn,*)
        read(nfcntn,*)
        read(nfcntn,*)
        read(nfcntn,*) rigid_body%velocity(1:3)
        read(nfcntn,*)
        read(nfcntn,*) rigid_body%angular_momentum(1:3)
        if(rigid_body%move_plane) then
        read(nfcntn,*) rigid_body%distance_moved_thusfar
        endif
      enddo
    endif

    if(npes>1) then
        do i=1,nrigid_bodies
          call mpi_bcast(rigid_bodies(i)%velocity,3,mpi_double_precision,0,MPI_CommGroup,ierr)
          call mpi_bcast(rigid_bodies(i)%angular_momentum,3,mpi_double_precision,0,MPI_CommGroup,ierr)
          if(rigid_bodies(i)%move_plane) then
          call mpi_bcast(rigid_bodies(i)%distance_moved_thusfar,1,mpi_double_precision,0,MPI_CommGroup,ierr)
          endif
        enddo
    endif
    return
  end subroutine m_IS_rd_rigid_body

  subroutine m_IS_rb_reinitialize()
      integer :: i,j, it, ia, na
      do i=1,nrigid_bodies
        na = rigid_bodies(i)%natm
        rigid_bodies(i)%mass = 0.d0
        rigid_bodies(i)%COM  = 0.d0
        do j=1,na
          ia = rigid_bodies(i)%atoms(j)
          it = ityp(ia)
          rigid_bodies(i)%mass   = rigid_bodies(i)%mass   + amion(it)
          rigid_bodies(i)%COM(:) = rigid_bodies(i)%COM(:) + amion(it)*cps(ia,:)
        enddo
        rigid_bodies(i)%COM = rigid_bodies(i)%COM/rigid_bodies(i)%mass
      enddo
      call principle_axis_of_inertia()
      call build_quaternion_from_rotmat()
      call map_cps_to_molecular_coords()
      call map_cps_to_relative_coords()
  end subroutine m_IS_rb_reinitialize

  subroutine principle_axis_of_inertia()
    integer :: na, ia, i, j, it, k
    real(kind=DP) :: mass
    real(kind=DP), dimension(3,3)            :: inertia
    real(kind=DP), dimension(3)              :: tcps, eigenvalue
    integer, parameter                       :: lwork  = 20
    real(kind=DP), allocatable, dimension(:) :: work
    integer, parameter                       :: liwork = 20
    integer, dimension(liwork)               :: iwork
    integer                                  :: inf
    real(kind=DP), dimension(3,3)            :: mat33
    allocate(work(lwork))
    do i=1, nrigid_bodies
      inertia = 0.d0
      na = rigid_bodies(i)%natm
      do j=1,na
        ia           = rigid_bodies(i)%atoms(j)
        it           = ityp(ia)
        tcps(:)      = cps(ia,:) - rigid_bodies(i)%COM(:)
        inertia(1,1) = inertia(1,1) + amion(it)*(tcps(2)**2+tcps(3)**2)
        inertia(2,2) = inertia(2,2) + amion(it)*(tcps(3)**2+tcps(1)**2)
        inertia(3,3) = inertia(3,3) + amion(it)*(tcps(1)**2+tcps(2)**2)
        inertia(1,2) = inertia(1,2) - amion(it)*(tcps(1)*tcps(2))
        inertia(2,3) = inertia(2,3) - amion(it)*(tcps(2)*tcps(3))
        inertia(1,3) = inertia(1,3) - amion(it)*(tcps(3)*tcps(1))
      enddo

      call DSYEV('V','U',3,inertia,3,eigenvalue,work,lwork,inf)

      rigid_bodies(i)%inertia(:) = eigenvalue(:)
      if(rigid_bodies(i)%explicit_zaxis) then
        mat33(:,1) = rigid_bodies(i)%xaxis(:)
        mat33(:,2) = rigid_bodies(i)%yaxis(:)
        mat33(:,3) = rigid_bodies(i)%zaxis(:)
        rigid_bodies(i)%rotmat     = mat33
      else
        rigid_bodies(i)%rotmat     = transpose(inertia)
      endif
    enddo
    deallocate(work)
    return
  end subroutine principle_axis_of_inertia

  subroutine m_IS_wd_rigid_body(nfcntn)
    integer, intent(in) :: nfcntn
    integer :: i
    if(mype/=0 .or. nrigid_bodies==0) return
    write(nfcntn,'(a)') tag_rigid_body
    write(nfcntn,'(a)') ' (nrigid_bodies)'
    write(nfcntn,'(i8)') nrigid_bodies
    do i=1,nrigid_bodies
      call wd_rigid_body(rigid_bodies(i))
    enddo

    contains

    subroutine wd_rigid_body(rigid_body)
      type(rigid_body_t), intent(in) :: rigid_body
      integer :: ia
      write(nfcntn,'(a)') ' (ID)'
      write(nfcntn,'(i8)') rigid_body%id
!      write(nfcntn,'(a)') ' (natm)'
!      write(nfcntn,'(i8)') rigid_body%natm
!      write(nfcntn,'(a)') ' (atoms)'
!      do ia=1,rigid_body%natm
!        write(nfcntn,'(i8,9d24.16)') rigid_body%atoms(ia), rigid_body%relative_coords(ia,1:3), &
!                                   & rigid_body%molecular_coords(ia,1:3),rigid_body%force_per_atm(ia,1:3)
!      enddo
!      write(nfcntn,'(a)') ' (COM)'
!      write(nfcntn,'(3d24.16)') rigid_body%COM(1:3)
!      write(nfcntn,'(a)') ' (mass)'
!      write(nfcntn,'(d24.16)') rigid_body%mass
!      write(nfcntn,'(a)') ' (quaternion)'
!      write(nfcntn,'(4d24.16)') rigid_body%q(1:4)
!      write(nfcntn,'(a)') ' (force)'
!      write(nfcntn,'(3d24.16)') rigid_body%force(1:3)
      write(nfcntn,'(a)') ' (velocity)'
      write(nfcntn,'(6d24.16)') rigid_body%velocity(1:3)
      write(nfcntn,'(a)') ' (angular_momentum)'
      write(nfcntn,'(3d24.16)') rigid_body%angular_momentum(1:3)
      if(rigid_body%move_plane) then
      write(nfcntn,'(d24.16)') rigid_body%distance_moved_thusfar
      endif
!      write(nfcntn,'(a)') ' (torque)'
!      write(nfcntn,'(3d24.16)') rigid_body%torque(1:3)
!      write(nfcntn,'(a)') ' (inertia)'
!      write(nfcntn,'(3d24.16)') rigid_body%inertia(1:3)
!      write(nfcntn,'(a)') ' (rotmat)'
!      write(nfcntn,'(3d24.16)') rigid_body%rotmat(1,1:3)
!      write(nfcntn,'(3d24.16)') rigid_body%rotmat(2,1:3)
!      write(nfcntn,'(3d24.16)') rigid_body%rotmat(3,1:3)
    end subroutine wd_rigid_body

  end subroutine m_IS_wd_rigid_body

  subroutine read_fix_bond_input()
    integer :: i,j,n,iret
    character(len=FMAXVALLEN) :: rstr
    integer ::  f_getIntValue, f_getRealValue, f_getStringValue
    logical :: found
    character(len=256) :: idstr
    n = 0

    if(f_getStringValue(tag_target_element,rstr, LOWER)==0) then
      found = .false.
      do j=1,ntyp
        if(rstr == speciesname(j)) then
          n = n+1
          found = .true.
        endif
      enddo
      if (.not.found) then
        if(printable) write(nfout,'(a)') '!** WARNING undefined element '//trim(adjustl(rstr))//'; will be ignored'
      endif
    endif
    i = 0
    do
      i = i+1
      write(idstr,*) i
      if(f_getStringValue(tag_target_element//trim(adjustl(idstr)),rstr,LOWER)==0) then
        n = n+1
      else
        exit
      endif
    enddo

    if(n==0) then
      fix_bond_nelements = 1
      allocate(fix_bond_elements(1))
      fix_bond_elements(1) = 1
      return
    endif

    fix_bond_nelements = n
    n = 0
    allocate(fix_bond_elements(fix_bond_nelements))
    if(f_getStringValue(tag_target_element,rstr,LOWER)==0) then
      do j=1,ntyp
        if(rstr == speciesname(j)) then
          n = n+1
          fix_bond_elements(n) = int(iatomn(j))
        endif
      enddo
    endif

    return
  end subroutine read_fix_bond_input

  subroutine set_covrads()
    allocate(covalent_radii(ntyp))
    call set_covrad_default(ntyp,iatomn,covalent_radii)
  end subroutine set_covrads

  subroutine read_fix_bond_options()
    integer :: f_getRealValue, f_getIntValue
    integer :: iret
    real(kind=DP) :: dret
    if(f_getRealValue(tag_bond_factor,dret,'')==0) then
      bond_factor = dret
    endif
    if(f_getIntValue(tag_max_iter_fix_bond,iret)==0) then
      max_iter_fix_bond = iret
    endif
    if(f_getRealValue(tag_thres_fix_bond,dret,'')==0) then
      thres_fix_bond = dret
    endif
  end subroutine read_fix_bond_options

  subroutine resolve_bonds_to_be_fixed()
    integer :: ia,ja,ib,jb,ie, nbonds
    real(kind=DP) :: blen2max, dist2
    integer, allocatable, dimension(:) :: ityp_work
    integer, allocatable, dimension(:,:) :: points_work
    logical :: found
    allocate(ityp_work(natm))
    do ia=1,natm
      ityp_work(ia) = int(iatomn(ityp(ia)))
    enddo
    nbonds = 0
    allocate(points_work(natm*natm,2));points_work=0
    do ib=1,fix_bond_nelements
      ie = fix_bond_elements(ib)
      do ia=1,natm
        if(ityp_work(ia)==ie) then
          do ja=1,natm
            if(ia==ja) cycle
            blen2max = (bond_factor*covalent_radii(ityp(ia))+covalent_radii(ityp(ja)))**2
            dist2 = dot_product(cps(ia,:)-cps(ja,:),cps(ia,:)-cps(ja,:))
            if (dist2<blen2max) then
              found = .false.
              do jb=1,nbonds
                if((points_work(jb,1)==ia .and. points_work(jb,2)==ja)  .or. &
                &  (points_work(jb,1)==ja .and. points_work(jb,2)==ia)) then
                  found = .true.
                  exit
                endif
              enddo
              if(.not.found) then
                nbonds = nbonds+1
                points_work(nbonds,1) = ia
                points_work(nbonds,2) = ja
              endif
            endif
          enddo
        endif
      enddo
    enddo
    nfixed_bonds = nbonds
    nbonds = 0
    points_work = 0
    allocate(fixed_bond(nfixed_bonds))
    allocate(bond_forc(natm,3));bond_forc=0.d0
    do ib=1,fix_bond_nelements
      ie = fix_bond_elements(ib)
      do ia=1,natm
        if(ityp_work(ia)==ie) then
          do ja=1,natm
            if(ia==ja) cycle
            blen2max = (bond_factor*covalent_radii(ityp(ia))+covalent_radii(ityp(ja)))**2
            dist2 = dot_product(cps(ia,:)-cps(ja,:),cps(ia,:)-cps(ja,:))
            if (dist2<blen2max) then
              found = .false.
              do jb=1,nbonds
                if((points_work(jb,1)==ia .and. points_work(jb,2)==ja)  .or. &
                &  (points_work(jb,1)==ja .and. points_work(jb,2)==ia)) then
                  found = .true.
                  exit
                endif
              enddo
              if(.not.found) then
                nbonds = nbonds+1
                points_work(nbonds,1) = ia
                points_work(nbonds,2) = ja
                fixed_bond(nbonds)%atoms(1) = ia
                fixed_bond(nbonds)%atoms(2) = ja
                fixed_bond(nbonds)%ityp(1) = ityp(ia)
                fixed_bond(nbonds)%ityp(2) = ityp(ja)
                fixed_bond(nbonds)%bond_length = sqrt(dist2)
                call update_bond_sigma_dsigma(fixed_bond(nbonds))
              endif
            endif
          enddo
        endif
      enddo
    enddo
    call update_bond_dsigma_old()
    deallocate(points_work)
    deallocate(ityp_work)
    deallocate(fix_bond_elements)
  end subroutine resolve_bonds_to_be_fixed

  subroutine print_bonds_to_be_fixed()
    integer :: ib
    if(printable) then
      write(nfout,'(a,i8)') '!** number of bonds to be fixed during MD',nfixed_bonds
      do ib=1,nfixed_bonds
        write(nfout,'(a,i8,a,2i8,a,f10.5)') &
        &                  '!** bond no ',ib,' assc. atoms ',fixed_bond(ib)%atoms(1),fixed_bond(ib)%atoms(2), &
        &                  ' bond length ', fixed_bond(ib)%bond_length
      enddo
!      write(nfout,'(a)')       '!** fix_bond options'
!      write(nfout,'(a,f8.2)')  '!** bond factor            ',bond_factor
!      write(nfout,'(a,i8)')    '!** SHAKE/RATTLE max iter  ',max_iter_fix_bond
!      write(nfout,'(a,e15.5)') '!** SHAKE/RATTLE threshold ',thres_fix_bond
    endif
  end subroutine print_bonds_to_be_fixed

  subroutine update_bond_sigma_dsigma(fbond)
    type(fix_bond_t), intent(inout) :: fbond
    real(kind=DP), dimension(3) :: cps1, cps2, diffvec
    real(kind=DP) :: distance2
    cps1(1:3)    = cps(fbond%atoms(1),1:3)
    cps2(1:3)    = cps(fbond%atoms(2),1:3)
    diffvec(1:3) = cps1(1:3)-cps2(1:3)
    distance2    = dot_product(diffvec,diffvec)
    fbond%sigma  =  distance2 - fbond%bond_length**2
    fbond%dsigma(1,:) =  2.d0*diffvec(:)
    fbond%dsigma(2,:) = -2.d0*diffvec(:)
  end subroutine update_bond_sigma_dsigma

  subroutine update_bond_dsigma_old()
    integer :: i
    do i=1,nfixed_bonds
      fixed_bond(i)%dsigma_old = fixed_bond(i)%dsigma
    enddo
  end subroutine update_bond_dsigma_old

  subroutine fixed_bond_coords()
    real(kind=DP) :: currlambda,dt2,dt2inv,factor,denom
    integer :: i,j,k,i1,ierr,ia1,ia2
    logical :: all_bonds_conv
    real(kind=DP) :: maxsig,muij,a,b,c
    real(kind=DP),dimension(2,3) :: cps_org
    if(mype==0) then
      dt2 = 0.5d0*dtio*dtio
      do i=1, nfixed_bonds
        fixed_bond(i)%lambda = 0.d0
      enddo
      dt2inv = 1.d0/dtio/dtio
      do i=1,max_iter_fix_bond
        maxsig = 0.d0
        do j=1,nfixed_bonds
          ia1 = fixed_bond(j)%atoms(1)
          ia2 = fixed_bond(j)%atoms(2)
          muij = 1.d0/amion(ityp(ia1))+1.d0/amion(ityp(ia2))
          a = muij*muij*dot_product(oldcps(ia1,:)-oldcps(ia2,:),oldcps(ia1,:)-oldcps(ia2,:))
          b = 2.d0*muij*dot_product(oldcps(ia1,:)-oldcps(ia2,:),cps(ia1,:)-cps(ia2,:))
          c = dot_product(cps(ia1,:)-cps(ia2,:),cps(ia1,:)-cps(ia2,:))-fixed_bond(j)%bond_length**2
          currlambda = (-b+sqrt(b**2-4.d0*a*c))/(2.d0*a)/(2.d0*dtio*dtio)
          fixed_bond(j)%lambda = fixed_bond(j)%lambda+currlambda
          do k=1,2
            i1 = fixed_bond(j)%atoms(k)
            !factor = dt2/amion(ityp(i1))
            factor = dtio*dtio/amion(ityp(i1))
            !cps(i1,:) = cps(i1,:) + factor * currlambda * fixed_bond(j)%dsigma_old(k,:)
            cps(i1,:) = cps(i1,:) + factor * currlambda * fixed_bond(j)%dsigma_old(k,:)
          enddo
          call update_bond_sigma_dsigma(fixed_bond(j))
        enddo
        all_bonds_conv = .true.
        do j=1,nfixed_bonds
          if(abs(fixed_bond(j)%sigma) .gt. thres_fix_bond) then
            all_bonds_conv = .false.
            exit
          endif
        enddo
        if(all_bonds_conv) then
          if(iprimd>=2) write(nfout,'(a,i8,a)') '!** fixed_bond : SHAKE converged after ',i,' iterations'
          exit
        endif
        if(i.eq.max_iter_fix_bond) then
          call phase_error_with_msg(nfout,'!** fixed_bond : SHAKE iteration unconverged',__LINE__,__FILE__)
        endif
      enddo
    endif
    if(npes>1) then
      do j=1,nfixed_bonds
        call mpi_bcast(fixed_bond(j)%lambda,1,mpi_double_precision,0,MPI_CommGroup,ierr)
      enddo
      call mpi_bcast(cps,natm*3,mpi_double_precision,0,MPI_CommGroup,ierr)
    endif
  end subroutine fixed_bond_coords

  subroutine fixed_bond_coords_by_opt()
    real(kind=DP) :: currlambda,dt2,dt2inv,factor,denom
    integer :: i,j,k,i1
    logical :: all_bonds_conv
    real(kind=DP) :: maxsig
    dt2 = 0.5d0*dtio*dtio
    do i=1, nfixed_bonds
      fixed_bond(j)%lambda = 0.d0
    enddo
    dt2inv = 1.d0/dtio/dtio
    do i=1,max_iter_fix_bond
      maxsig = 0.d0
      do j=1,nfixed_bonds
        call update_bond_sigma_dsigma(fixed_bond(j))
        if(abs(fixed_bond(j)%sigma)>maxsig) maxsig = abs(fixed_bond(j)%sigma)
        denom = 0.d0
        do k=1,2
          i1 = fixed_bond(j)%atoms(k)
          factor = 0.5d0/amion(ityp(i1))
          denom  = denom + factor * dot_product(fixed_bond(j)%dsigma_old(k,:),fixed_bond(j)%dsigma(k,:))
        enddo
        if(denom.lt.eps_shake) then
          call phase_error_with_msg(nfout,'!** fixed_bond : SHAKE denom is too small',__LINE__,__FILE__)
        endif
        currlambda           = dt2inv * fixed_bond(j)%sigma/denom
        fixed_bond(j)%lambda = fixed_bond(j)%lambda+currlambda
        do k=1,2
          i1 = fixed_bond(j)%atoms(k)
          factor = dt2/amion(ityp(i1))
          cps(i1,:) = cps(i1,:) - factor * currlambda * fixed_bond(j)%dsigma_old(k,:)
        enddo
      enddo
      all_bonds_conv = .true.
      do j=1,nfixed_bonds
        if(abs(fixed_bond(j)%sigma) .gt. thres_fix_bond) then
          all_bonds_conv = .false.
          exit
        endif
      enddo
      if(all_bonds_conv) then
        if(iprimd>=2) write(nfout,'(a,i8,a)') '!** fixed_bond : SHAKE converged after ',i,' iterations'
        exit
      endif
      if(i.eq.max_iter_fix_bond) then
        call phase_error_with_msg(nfout,'!** fixed_bond : SHAKE iteration unconverged',__LINE__,__FILE__)
      endif
    enddo
  end subroutine fixed_bond_coords_by_opt

  subroutine fixed_bond_velocities(veloc)
    real(kind=DP), dimension(natm,3), intent(inout) :: veloc
    integer :: i,j,k,i1
    real(kind=DP) :: numera, denom, tempeta, massinv, tmpsum, tmpeta
    logical :: all_bonds_conv
    do i=1,max_iter_fix_bond
      do j=1,nfixed_bonds
        numera = 0.d0
        denom  = 0.d0
        do k=1,2
          i1 = fixed_bond(j)%atoms(k)
          numera = numera + dot_product(veloc(i1,:),fixed_bond(j)%dsigma(k,:))
          denom  = denom  + dot_product(fixed_bond(j)%dsigma(k,:),fixed_bond(j)%dsigma(k,:))/amion(ityp(i1))
        enddo
        if(denom.lt.eps_rattle) then
          call phase_error_with_msg(nfout,'!** fixed_bond : RATTLE denom is too small',__LINE__,__FILE__)
        endif
        tmpeta = numera/denom
        do k=1,2
          i1 = fixed_bond(j)%atoms(k)
          massinv = 1.d0/amion(ityp(i1))
          veloc(i1,:) = veloc(i1,:) - massinv*tmpeta*fixed_bond(j)%dsigma(k,:)
        enddo
      enddo
      all_bonds_conv = .true.
      do j=1, nfixed_bonds
        tmpsum = 0.d0
        do k=1,2
          i1 = fixed_bond(j)%atoms(k)
          tmpsum = tmpsum + dot_product(veloc(i1,:),fixed_bond(j)%dsigma(k,:))
        enddo
        if(tmpsum.gt.thres_fix_bond) then
          all_bonds_conv = .false.
          exit
        endif
      enddo
      if(all_bonds_conv) then
        if(iprimd>=2) write(nfout,'(a,i8,a)') '!** fixed_bond : RATTLE converged after ',i,' iterations'
        exit
      endif
      if(i.eq.max_iter_fix_bond) then
        call phase_error_with_msg(nfout,'!** fixed_bond : RATTLE iteration unconverged',__LINE__,__FILE__)
      endif
    enddo
  end subroutine fixed_bond_velocities

  subroutine build_nbonds_per_thermo()
    external   ibath
    integer :: ibath
    integer :: ib, ir
    if(.not.allocated(nbonds_per_thermo)) allocate(nbonds_per_thermo(nrsv))
    nbonds_per_thermo = 0
    do ib=1,nfixed_bonds
       ir = ibath(imdtyp(fixed_bond%atoms(1)))
       if(ir >= 1) then
         nbonds_per_thermo(ir) = nbonds_per_thermo(ir)+1
       endif
    enddo
  end subroutine build_nbonds_per_thermo

  subroutine update_bond_force()
    integer :: ib,j
    bond_forc=0.d0
    do ib=1,nfixed_bonds
      do j=1,2
        bond_forc(fixed_bond(ib)%atoms(j),:) = bond_forc(fixed_bond(ib)%atoms(j),:) &
        &                                    - fixed_bond(ib)%lambda*fixed_bond(ib)%dsigma(j,:)
      enddo
    enddo
  end subroutine update_bond_force

  subroutine map_constraint()
    integer :: ib,icount,fix_type
    icount = 0
    nfcatm = 2*nfixed_bonds
    call m_IS_alloc_cnstrvectors_etc(imdalg)
    num_fixed_bonds = nfixed_bonds
    call alloc_bondlength_fix_set()
    fix_type = BONDLENGTH_FIX
    if(imdalg == QUENCHED_CONSTRAINT) fix_type = BONDLENGTH_FIX_1
    do ib=1,nfixed_bonds
      icount = icount + 1
      ia_cnst(icount) = fixed_bond(ib)%atoms(1)
      icount = icount + 1
      ia_cnst(icount) = fixed_bond(ib)%atoms(2)
      bondlength_fix_set(1,ib) = fixed_bond(ib)%atoms(1)
      bondlength_fix_set(2,ib) = fixed_bond(ib)%atoms(2)
      imdtyp(fixed_bond(ib)%atoms(1))      = fix_type
      imdtypxyz(fixed_bond(ib)%atoms(1),:) = fix_type
      imdtyp(fixed_bond(ib)%atoms(2))      = fix_type
      imdtypxyz(fixed_bond(ib)%atoms(2),:) = fix_type
    enddo
    constraint_type = fix_type
    constraints_exist = .true.
  end subroutine map_constraint

  function get_average_temperature() result(ret)
    integer :: it
    integer :: ntot
    real(kind=DP) :: ret
    ret = 0.d0
    ntot = 0
    do it=1,nrsv
      ntot = ntot+natm_per_thermo(it)
    enddo
    do it=1,nrsv
      ret = ret+natm_per_thermo(it)*tkb(it)/CONST_kB/dble(ntot)
    enddo
  end function get_average_temperature

  function get_default_mbaro(tdamp) result(res)
    real(kind=DP), intent(in) :: tdamp
    real(kind=DP) :: res
    res = (1.d0/univol)*CONST_kB*get_average_temperature()*(tdamp/PAI2)**2

    if (p_ctrl_method==LATTICE_VECTOR) then
      res = 2.d0*res*univol**(2.d0/3.d0)
    endif
  end function get_default_mbaro

end module m_Ionic_System

