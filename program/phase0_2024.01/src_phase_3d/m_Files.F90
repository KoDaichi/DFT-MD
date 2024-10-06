!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 633 $)
!
!  MODULE: m_Files
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!  patch 10.1 by K. Tagami @adv    2011/06/18
!
!  patch 10.1 : addition of output file for LinearResponse Spectrum
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
module m_Files
!  $Id: m_Files.F90 633 2020-12-01 05:11:03Z jkoga $
!
!  Operations concerning to files as "open", and "close",
!  should be done in this module.
!
  use m_Control_Parameters,only : paramset, ekmode, printable, iprijobstatus,jobstatus_series, ipriparadeb &
       &                        , multiple_replica_mode, icond, ppprinted, fixed_charge_k_parallel &
       &                        , sw_cif_output, imdalg
  use m_Const_Parameters,  only : OLD, NEW_, UNKNOWN, FORMATTED &
       &                        , UNFORMATTED, check_file_name_on &
       &                        , check_file_name_off, DP, ON, GENERAL, FILE &
       &                        , YES, NO, UNIFIED, PARTITIONED, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                        , ON, OFF, INITIAL, ALL_AT_ONCE, ONE_BY_ONE &
       &                        , P_CONTROL, PT_CONTROL
  use m_ErrorMessages,     only : F_POT_FILE_NOT_EXIST, FILE_NOT_EXIST, ERROR_IN_INPUTFILE_OPENING,FILENAMES_FORMAT_ERROR &
       &                        , FILENAMES_FORMAT_ERROR_NEB, FILENAMES_NOT_EXIST,m_EMsg_Warning
  use m_Parallelization, only   : MPI_CommGroup,mype,npes,ierr,workdir,sw_wdir &
       &                        , mype_conf,conf_para, nrank_conf
  use mpi

  implicit none
!  include 'mpif.h'
  integer istatus(mpi_status_size)

  logical, public ::  file_existence_contfiles = .false.
  logical, public ::  file_existence_contfiles_when_paw_on = .false.
  logical, public ::  file_existence_3contfiles = .false.
  logical, private :: file_existence_nfcntn = .false.
  logical, private :: file_existence_nfcntn_bin = .false.
  logical, private :: file_existence_nfcntn_bin_paw = .false.  ! T.Y. 2020/02/28
  logical, private :: file_existence_nfzaj =  .false.
  logical, private :: file_existence_nfchgt = .false.
  logical, private :: file_existence_nfoccmat = .false.        ! T.Y. 2020/02/28

  integer, parameter :: min_ext_length_mype = 3
  integer, parameter ::  MAXNSP = 16
  integer  nfinp,nfpot(MAXNSP)&  ! (1 + MAXNSP) files
       & ,nfpkb,nfpd,nfppc,nfstop,nfopgr,nfmatbp,nfkpoint,nfkindex & !8 files
       & ,nfvibrate,nfzak,nf_neord,nfcnst,nfotp,nfdynm,nfenf,nfgpt & !8 files
       & ,nfchgt,nfchgo,nfchgu,nfchgs,nfout,nfcps,nffor,nfcntn &     !8 files
       & ,nfcntn_bin,nfzaj,nfrsp,nfwf1,nfwf2,nfwf3,nfeng,nfvlc &     !8 files
       & ,nfcntn_bin_stm, nfstatus &                                 !2 files
       & ,nfspg,nfkpgn  &                                            !2 files
       & ,nfdos, nfchr, nfelf, nfegrid, nfldos &                     !5 files
       & ,nfberry,nfeffchg,nfforce,nfmode,nfepsilon &                !5 files
       & ,nfepsout,nfphdos,nfoccmat,nfwfk,nfstrfrc &                 !5 files
       & ,nfwannier,nfcntn_wannier,nfpot_wannier,nfefermi &          !4 files
! --------------------T. Hamada 2021.9. -------------------- ! UVSOR
!      & ,nfnlo,nfmagopt,nfepscont   &                               !3 files UVSOR
       & ,nfnlo,nfmagopt,nfepscont,nfkptdata,nfkptepsdata &                    !5 files
! ---------------------------------------------------------- ! UVSOR
       & ,nfpstrn, nfvelec, nfeppair, nfeppair2, nfvelec_grad &      !5 files Positron
       & ,nfpsicoef, nfbandsyminput &                                !2 file
! ============================= Added by K. Tagami =========== 10.1
!       & ,nfcntn_bin_paw                                            !1 file PAW
       & ,nfcntn_bin_paw              &                              !1 file PAW
       &, nf_LR_spectra               &                              !1 file LR
! ============================================================ 10.1
       &, nfcntn_berry &
       &, nf_ekindens             &                          ! 1 file 360 metagga
       &, nfqpoint                &                          ! 1 file q-point phonon
       &, nf_core_charge_radial(MAXNSP) &                      ! MAXNSP files 500-515

! ====================== KT_add ======= 13.1R
       &, nfeps_ph &                                         ! 1 file EPS_Phonon
       &, nfoptical_coeff &                                  ! 1 file optical coeff
       &, nframan_spectra &                                  ! 1 file Raman spectra
! ===================================== 13.1R
! ===== KT_add ==== 13.1XI
       &, nf_excitation_spectra &                            ! 1 file Excitation
! ================= 13.1XI
       &, nf_mag_core &                                      ! 1 file Opencore
       &, nf_elecpot_bin, nf_elecpot_bin_ref &
                                                             ! 2 files ChargedDefect
! ====================== KT_add ======= 13.0S
       &, nfcore_energy_out, nfcore_energy_initial, nfcore_energy_final &
                                                             ! 3 files CoreLevel
! ===================================== 13.0S
       &, nfdynm_cif & ! CIF output
       &, nfwfk_sq   &   ! squared wf
       &, nfwfk_integ_mom  &  ! moment of wf ( integerated over space )
       &, nfwfk_orb_proj      &      ! orbital-projection of wf
       &, nfwfk_local_decomp  &      ! local-decomposition of wf (lband)
       &, nfband_spectr_wght  &      ! spectral weight for band unfolding
       &, nfporb_dens_mat     &      ! orbital population density matrix
       &, nfporb_rot_mat      &      ! rotation matrix of orbital population
       &, nflatconst,nfmetric & ! for pressurecontrol simulations
       &, nfdftd3par & ! parameters for DFT-D3
       &, nfchkpoint & ! checkpoint file directory
       &, nfhybridinfo & ! hybrid functional information
       &, nfscfzaj & ! WF from the preceding SCF calculation
       &, nfpos, nfcoordattr &! file from which the coordinates are read when structure.method = file
       &, nfdefaultpp &
       &, nfpwbs &
       &, nfxsf, nfn2p2, nfinp_mod &
       &, nfzaj_kall, nfkpoint_bin

  data  &
       &   nfinp,nfpot &                                               ! 31,(37,38,39,40,45,46,11-19,36)
       &  ,nfpkb,nfpd,nfppc,nfstop,nfopgr,nfmatbp,nfkpoint,nfkindex &  ! 33,34,35,49,21,22,23,24
       &  ,nfvibrate,nfzak,nf_neord,nfcnst,nfotp,nfdynm,nfenf,nfgpt &  ! 47,48,50,61,51,52,53,54
       &  ,nfchgt,nfchgo,nfchgu,nfchgs,nfout,nfcps,nffor,nfcntn &      ! 55,56,57,58, 6, 6, 6,42
       &  ,nfcntn_bin,nfzaj,nfrsp,nfwf1,nfwf2,nfwf3,nfeng,nfvlc &      ! 43,44,90,95,96,97,59,60
       &  ,nfcntn_bin_stm, nfstatus &                                  ! 62,75
       &  ,nfspg,nfkpgn &                                              ! 63,64
       &  ,nfdos, nfchr, nfelf, nfegrid, nfldos &                      ! 5 files 65,66,67,68,69
       &  ,nfberry,nfeffchg,nfforce,nfmode,nfepsilon &                 ! 5 files 70,71,72,73,74
       &  ,nfepsout,nfphdos,nfoccmat,nfwfk,nfstrfrc &                  ! 5 files 20,41,76,77,78
       &  ,nfwannier,nfcntn_wannier,nfpot_wannier,nfefermi &           ! 4 files 79,80,81,82
!-------------------- T. Hamada 2021.9.28 --------------------! UVSOR
!       &  ,nfnlo,nfmagopt,nfepsconta  &                               ! 3 files 10,30,83 !UVSOR
       &  ,nfnlo,nfmagopt,nfepscont,nfkptdata,nfkptepsdata &           ! 5 files 10,30,83,84,85  UVSOR
!------------------------------------------------------------=! UVSOR

       &  ,nfpstrn, nfvelec, nfeppair, nfeppair2, nfvelec_grad &       ! 5 files 300,310,320,322,330 Positron
       &  ,nfpsicoef, nfbandsyminput &                                 ! 2 file  340, 341
       &  ,nfcntn_bin_paw &                                            ! 1 file  321 PAW
! ===================== Added by K. Tagami ================= 10.1
       &  ,nf_LR_spectra           &                     ! 1 file  400 Linear Response
! ========================================================== 10.1
       &  ,nfcntn_berry           &
       &  ,nf_ekindens             &                     ! 1 file 360 metagga
       &  ,nfqpoint                &                     ! 1 file 325 q-point for phonon
       &  ,nf_core_charge_radial     &             ! MAXNSP files 500-515

! ======= KT_add === 13.1R
       &   , nfeps_ph &                                  ! 1 file 407 Eps_Phonon
       &   , nfoptical_coeff &                           ! 1 file 408 optical coeff
       &   , nframan_spectra  &                          ! 1 file 409 Raman spectra
! ================== 13.1R
! ===== KT_add ==== 13.1XI
       &   , nf_excitation_spectra  &                    ! 1 file 420 Excitation spectra
! ================= 13.1XI
       &   , nf_mag_core &                               ! 1 files 378
       &   , nf_elecpot_bin, nf_elecpot_bin_ref &        ! 2 files 380, 381 defect
! ============================== KT_add ========= 13.0S
       &  ,nfcore_energy_out, nfcore_energy_initial, nfcore_energy_final &
! =============================================== 13.0S
       &  ,nfdynm_cif &
       &  ,nfwfk_sq &
       &  ,nfwfk_integ_mom &
       &  ,nfwfk_orb_proj &
       &  ,nfwfk_local_decomp &
       &  ,nfband_spectr_wght &
       &  ,nfporb_dens_mat, nfporb_rot_mat &
! NPT, NPH-MD
       &  ,nflatconst &
       &  ,nfmetric &
       &  ,nfdftd3par &
       &  ,nfchkpoint &
       &  ,nfhybridinfo &
       &  ,nfscfzaj &
       &  ,nfpos, nfcoordattr &
       &  ,nfdefaultpp &
       &  ,nfpwbs &
       &  ,nfxsf, nfn2p2 &
       &  , nfinp_mod, nfzaj_kall &
       &  , nfkpoint_bin  &

       &    /31,37,38,39,40,45,46 &
       & ,11,12,13,14,15,16,17,18,19,36 &
       & ,33,34,35,49,21,22,23,24 &
       & ,47,48,50,61,51,52,53,54 &
       & ,55,56,57,58, 6, 6, 6,42 &
       & ,43,44,90,95,96,97,59,60 &
       & ,62,75,63,64,65,66,67,68,69 &
       & ,70,71,72,73,74,20,41,76,77,78 &
       & ,79,80,81,82 &
! -------------------- T. Hamada 2021.9.23 --------------------
!       & ,10,30,83, 300,310,320,322,330,340,341 &
       & ,10,30,83,84, 85, 300,310,320,322,330,340,341 &
! -------------------------------------------------------------
! ================================ Added by K. Tagami ========== 10.1
       & ,321 &
       & ,400 &
! ============================================================== 10.1
       & ,410 &
       & ,360 &
       & ,325 &
       & ,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515 &

! =============== KT_add =========== 13.1R
       & ,407 &
       & ,408 &
       & ,409 &
! ================================== 13.1R
! =============== KT_add ====== 13.1XI
       & ,420 &
! ============================= 13.1XI
       &, 378 &           ! mag-core
       &, 380, 381 &
! =============== KT_add =========== 13.0S
       &, 370, 371, 372  &
! ================================== 13.0S
       &, 349 &
       &, 350 &
       &, 351 &
       &, 352 &
       &, 345 &
       &, 346 &
       &, 355, 356 &
! NPT, NPH-MD
       &, 100 &
       &, 101 &
       &, 102 &
       &, 103 &
       &, 104 &
       &, 105 &
       &, 106,107 & ! dummy unit number
       &, 108 &
       &, 109 &
       &, 111 &
       &, 112 &
       &, 113 &
       &, 114 &
       &, 115 /

  integer,private,parameter :: number_of_all_files = 104 + 2*MAXNSP

  integer,private, dimension(number_of_all_files) :: n_file
  data n_file &
       &    /31,37,38,39,40,45,46 &
       & ,11,12,13,14,15,16,17,18,19,36 &
       & ,33,34,35,49,21,22,23,24 &
       & ,47,48,50,61,51,52,53,54 &
       & ,55,56,57,58, 6, 6, 6,42 &
       & ,43,44,90,95,96,97,59,60 &
       & ,62,75,63,64,65,66,67,68,69 &
       & ,70,71,72,73,74,20,41,76,77,78 &
       & ,79,80,81,82 &
! -------------------- T. Hamada 2021.9.28 --------------------
!      & ,10,30, 300,310,320,330,340 &
       & ,10,30,83,84,85,300,310,320,330,340 &
! -------------------------------------------------------------
! ===================================== Added by K. Tagami ======= 10.1
!       & ,321/
       & ,321 &
       & ,400 &
! ================================================================ 10.1
       & ,410 &
       & ,360 &
       & ,325 &
       & ,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515 &

! =============== KT_add =========== 13.1R
       & ,407 &
       & ,408 &
       & ,409 &
! ================================== 13.1R
! =============== KT_add ====== 13.1XI
       & ,420 &
! ============================= 13.1XI
       & ,378 &          ! mag-core
       & ,380, 381 &
! ================= KT_add ======== 13.0S
       & ,370, 371, 372 &
! ================================= 13.0S
       & ,349, 350, 351, 352, 345, 346, 355, 356 &
! NPT, NPH-MD
       & ,100, 101 &
       & ,102 &
       & ,103 &
       & ,104 &
       & ,105 &
       & ,106, 107 &
       & ,108 &
       & ,109 &
       & ,111 &
       & ,112 &
       & ,113 &
       & ,114 &
       & ,115 /

  integer,private,parameter :: stringlength_of_filenames = 261
  character(len=stringlength_of_filenames) ::  &
       &  F_INP,F_POT(MAXNSP),F_PKB,F_PD,F_PPC,F_STOP,F_OPGR &         ! 31,(37-40,45,46,11-19,36),33,34,35,49,21
       & ,F_MATBP,F_KPOINT,F_KINDEX, F_VIBRATE, F_ZAK,F_NEORD,F_CNST & ! 22,23,24,47,48,50,61
       & ,F_OTP, F_DYNM,F_ENF, F_GPT, F_CHGT, F_CHGO, F_CHGU, F_CHGS & ! 51,52,53,54,55,56,57,58
       & ,F_OUT, F_CNTN, F_CNTN_BIN, F_ZAJ, F_RSP, F_WF(3), F_ENERG  & ! 6, 42,43,44,90,95,96,97,59
       & ,F_VLC, F_CNTN_BIN_STM, F_STATUS, F_file_names_data &         ! 60,62,75
       & ,F_SPG, F_KPGN, F_DOS, F_CHR, F_ELF, F_EGRID, F_LDOS &        ! 63,64,65,66,67,68,69
       & ,F_BERRY,F_EFFCHG,F_FORCE,F_MODE,F_EPSILON &                  ! 70,71,72,73,74
       & ,F_EPSOUT,F_PHDOS,F_OCCMAT,F_WFk,F_STRFRC &                   ! 20,41,76,77,78
       & ,F_WANNIER,F_CNTN_WAN,F_POT_WAN,F_EFERMI &                    ! 79,80,81,82
! ^^^^^^^^^^^^^^^^^^^^^^^T. Hamada 2021.9.28 ---------------------------! UVSOR
!       & ,F_NLO,F_MAGOPT,F_EPSCONT                                    ! 10,30,83 UVSOR
       & ,F_NLO,F_MAGOPT,F_EPSCONT,F_KPT_DATA,F_KPTEPS_DATA &          ! 10,30,83,84,85 UVSOR
! ----------------------------------------------------------------------! UVSOR
       & ,F_PSTRN, F_VELEC, F_EPPAIR, F_EPPAIR2, F_VELEC_GRAD  &       ! 300,310,320,322,330 Positron
       & ,F_PSICOEF, F_BAND_SYM_INPUT       &                          ! 340,341 Psi_coef, Psi_coef2
! ================================== Added by K. Tagami ======= 10.1
!       & ,F_CNTN_BIN_PAW                                               ! 321 PAW
       & ,F_CNTN_BIN_PAW &                                              ! 321 PAW
       & ,F_LR_SPECTRA   &                                         ! 400 LinearResponse
! ============================================================= 10.1
       & ,F_CNTN_BERRY &
       & ,F_EKINDENS &                                             ! 360 metagga
       & ,F_QPOINT   &                                ! 325 Q-POINT
       &, F_CORE_CHARGE(MAXNSP)   &              ! 500-515
! ================= KT_add ====== 13.1R
       &, F_EPS_PHONON &                              ! 407 EPS_Phonon
       &, F_OPTICAL_COEFF &                           ! 408 Optical Coeff
       &, F_RAMAN_SPECTRA &                           ! 409 Raman spectra
! =============================== 13.1R

! ================= KT_add ====== 13.1XI
       &, F_EXCITATION_SPECTRA &                       ! 420 Excitation spectra
! =============================== 13.1XI
       &, F_MAG_CORE &                                 ! 378 Mag OpenCore
       &, F_ELECPOT_BIN, F_ELECPOT_BIN_REF &           ! 380, 381 Defect
! ================= KT_add ====== 13.0S
       &, F_CORE_ENERGY_OUT, F_CORE_ENERGY_INITIAL,  F_CORE_ENERGY_FINAL &
! =============================== 13.0S
       &, F_DYNM_CIF, F_WFk_Squared, F_WFk_IntegMoment &
       &, F_WFK_ORB_PROJ, F_WFK_LOCAL_DECOMP &
       &, F_BAND_SPECTR_WGHT, F_PORB_DENS_MAT, F_PORB_ROT_MAT &
! NPT, NPH-MD
       &, F_LATCONST, F_METRIC &
       &, F_DFTD3PAR &
       &, F_CHKPNT &
       &, F_HYBRIDINFO &
       &, F_SCF_ZAJ &
       &, F_POS, F_COORD_ATTR, F_DEFAULTPP, F_PWBS &
       &, F_XSF, F_N2P2, F_INP_MOD, F_ZAJ_KALL, F_KPOINT_BIN &
       &, F_POT_ORG(MAXNSP)

  namelist/fnames/ &
       &  F_INP, F_POT, F_PKB, F_PD, F_PPC, F_STOP, F_OPGR &
       & ,F_MATBP,F_KPOINT,F_KINDEX, F_VIBRATE, F_ZAK,F_NEORD,F_CNST &
       & ,F_OTP, F_DYNM,F_ENF, F_GPT, F_CHGT, F_CHGO, F_CHGU, F_CHGS &
       & ,F_OUT, F_CNTN, F_CNTN_BIN, F_ZAJ, F_RSP,  F_WF,   F_ENERG  &
       & ,F_VLC, F_CNTN_BIN_STM, F_STATUS &
       & ,F_SPG, F_KPGN, F_DOS, F_CHR, F_ELF, F_EGRID, F_LDOS &
       & ,F_BERRY,F_EFFCHG,F_FORCE,F_MODE,F_EPSILON &
       & ,F_EPSOUT,F_PHDOS,F_OCCMAT,F_WFk,F_STRFRC &
       & ,F_WANNIER,F_CNTN_WAN,F_POT_WAN,F_EFERMI &
!-------------------- T. Hamada 2021.9.28 -------------------- ! UVSOR
!      & ,F_NLO,F_MAGOPT,F_EPSCONTA  &
       & ,F_NLO,F_MAGOPT,F_EPSCONT,F_KPT_DATA,F_KPTEPS_DATA  &
!------------------------------------------------------------- ! UVSOR
       & ,F_PSTRN, F_VELEC, F_EPPAIR, F_EPPAIR2, F_VELEC_GRAD &        ! Positron
       & ,F_PSICOEF, F_BAND_SYM_INPUT &
! ================================== Added by K. Tagami ======= 10.1
       & ,F_CNTN_BIN_PAW &
       & ,F_LR_SPECTRA   &
! ============================================================= 10.1
       & ,F_CNTN_BERRY &
       & ,F_EKINDENS &
       & ,F_QPOINT &
       &, F_CORE_CHARGE &
! ================= KT_add ====== 13.1R
       &, F_EPS_PHONON &
       &, F_OPTICAL_COEFF &
       &, F_RAMAN_SPECTRA &
! =============================== 13.1R

! ================= KT_add === 13.1XI
       &, F_EXCITATION_SPECTRA &
! ============================ 13.1XI
       &, F_MAG_CORE &
       &, F_ELECPOT_BIN, F_ELECPOT_BIN_REF &
! ================= KT_add ====== 13.0S
       &, F_CORE_ENERGY_OUT, F_CORE_ENERGY_INITIAL,  F_CORE_ENERGY_FINAL &
! =============================== 13.0S
       &, F_DYNM_CIF, F_WFk_Squared, F_Wfk_IntegMoment &
       &, F_WFK_ORB_PROJ, F_WFK_LOCAL_DECOMP &
       &, F_BAND_SPECTR_WGHT, F_PORB_DENS_MAT, F_PORB_ROT_MAT &
       &, F_LATCONST, F_METRIC &
       &, F_DFTD3PAR &
       &, F_CHKPNT &
       &, F_HYBRIDINFO &
       &, F_SCF_ZAJ &
       &, F_POS, F_COORD_ATTR &
       &, F_DEFAULTPP &
       &, F_PWBS &
       &, F_XSF, F_N2P2, F_INP_MOD &
       &, F_ZAJ_KALL, F_KPOINT_BIN


  logical ::             F_ZAJ_partitioned      = .false.
  logical ::             F_CHGT_partitioned     = .false.
  logical ::             F_CNTN_BIN_partitioned = .false.
  logical ::             F_CNTN_partitioned     = .false.
  logical ::             F_ZAJ_in_partitioned      = .false.
  logical ::             F_CHGT_in_partitioned     = .false.
  logical ::             F_CNTN_BIN_in_partitioned = .false.
  logical ::             F_CNTN_in_partitioned     = .false.

! ======================== KT_add ================= 13.0D
  logical ::             F_CNTN_BIN_PAW_partitioned = .false.
  logical ::             F_CNTN_BIN_PAW_in_partitioned  = .false.
! ================================================= 13.0D

  character(len=stringlength_of_filenames) &
       &  ::             F_ZAJ_filetype,      F_CHGT_filetype &
       &              ,  F_CNTN_BIN_filetype,    F_CNTN_filetype
  character(len=stringlength_of_filenames) &
       &  ::             F_ZAJ_in_filetype,   F_CHGT_in_filetype &
       &               , F_CNTN_BIN_in_filetype, F_CNTN_in_filetype

! ======================== KT_add ================= 13.0D
  character(len=stringlength_of_filenames) &
       &                 F_CNTN_BIN_PAW_filetype, F_CNTN_BIN_PAW_in_filetype
! ================================================= 13.0D

  namelist/mpifiletypes/ F_ZAJ_filetype,         F_CHGT_filetype &
       &               , F_CNTN_BIN_filetype,    F_CNTN_filetype &
       &               , F_ZAJ_in_filetype,      F_CHGT_in_filetype &
       &               , F_CNTN_BIN_in_filetype, F_CNTN_in_filetype &
! ======================== KT_add ================= 13.0D
       &               , F_CNTN_BIN_PAW_filetype, F_CNTN_BIN_PAW_in_filetype
! ================================================= 13.0D


  character(len=stringlength_of_filenames) &
       &  ::             F_ZAJ_in, F_CHGT_in, F_CNTN_BIN_in, F_CNTN_in
  character(len=stringlength_of_filenames) &
       &  ::             F_ZAJ_bak, F_CHGT_bak, F_CNTN_BIN_bak, F_CNTN_bak

! ======================== KT_add ================= 13.0D
  character(len=stringlength_of_filenames) :: F_CNTN_BIN_PAW_in
! ================================================= 13.0D

  logical, private :: fout_is_given = .false.
  character(len=stringlength_of_filenames+3) :: F_STATUS_ext
  integer, dimension(-1:1) :: nfstm

  integer nfprm
  data nfprm/93/
  character(len=stringlength_of_filenames) ::  F_PRM
  namelist/f_param_name/F_PRM

  character, private :: name_jobstep*3

  integer  nfimage, nfnebstop, nfneb, nfnebcntn, nfnebenf, nfnebdynm, nfpath
  data  nfimage, nfnebstop, nfneb, nfnebcntn, nfnebenf, nfnebdynm, nfpath &
	/201,202,203,204,205,206,207/
  character(len=stringlength_of_filenames) :: F_IMAGE(-1:99), F_NEB_STOP, F_NEB_OUT, F_NEB_CNTN &
         &                                  , F_NEB_ENF, F_NEB_DYNM, F_PATH
  namelist/nebfiles/ F_IMAGE, F_NEB_STOP, F_NEB_OUT, F_NEB_CNTN, F_NEB_ENF, F_NEB_DYNM, F_PATH

  character(len=5) :: prefix='_conf'

! === KT_add === 2014/07/20
  character(len=100) :: F_OUT_BASE, F_CONF_EXTENSION, F_PARA_EXTENSION
! ============== 2014/07/20

! ==== KT_add ==== 2014/07/14
  integer :: nfhypervec
  data nfhypervec /210/
  character(len=100) :: F_HYPERVEC
!
  namelist /constraints/ F_HYPERVEC

  integer :: conftag=-1


! ================ 2014/07/14
!   Following interface block is moved from subroutines of 'm_Files_open_ps_files' and 'm_Files_open_ps_file",
!  according to a report from ASMS Co.ltd, 10 March 2016.
    interface
       subroutine phase_error(ierror,nfout,nfpp,filename,line,file)
         integer,         intent(in)          :: ierror, nfout, nfpp
         character(len=*),intent(in)          :: filename
         integer,         intent(in),optional :: line
         character(len=*),intent(in),optional :: file
       end subroutine phase_error
    end interface

contains

  subroutine m_Files_set_conftag(tag)
    integer, intent(in) :: tag
    conftag = tag
  end subroutine m_Files_set_conftag

  subroutine m_Files_set_stdout(nfo)
    integer, intent(in) :: nfo
    nfout = nfo
  end subroutine m_Files_set_stdout

  subroutine m_Files_open_standardout()
    character :: name_mype*12
    integer :: name_length
    integer ::   i, js, ls, imax
    logical ::   existence

    character(len=256) :: cid,cketa, f_tmp
    integer :: iketa=2
    integer :: itmp

    if(npes > 1) then
       name_length = 12
       imax = int(log10(dble(npes-1)))+1
       if(imax < min_ext_length_mype) imax = min_ext_length_mype
       if(imax > name_length) then
          call phase_error_with_msg(nfout,' number of pes (=npes) is larger than 1,000,000,000,000',__LINE__,__FILE__)
       end if
       write(name_mype,'(i12)') mype
       do i = 1+name_length-imax, name_length
          if(name_mype(i:i) == ' ') name_mype(i:i)='0'
       enddo
    end if

    ls = len(trim(F_OUT))

! === KT_add === 2014/07/20
    F_CONF_EXTENSION = "";    F_PARA_EXTENSION = ""

    if ( mype /= 0 ) then
       F_PARA_EXTENSION = '_'//name_mype(name_length-imax+1:name_length)
    endif
! ============== 2014/07/20

    if(.not.conf_para .or. conf_para.and.mype_conf==0) then

    cid = ''
    if(conftag>=0) then
      itmp=int(log10(real(conftag)))+1
      if(itmp>2)iketa=itmp
      write(cketa,*) iketa
      write(cid,'(i'//trim(adjustl(cketa))//'.'//trim(adjustl(cketa))//')') conftag
      cid = '.'//cid
    endif
    if(ls < 1) then
       fout_is_given = .false.
       if(mype==0) then
          do js = 0, 999
             write(name_jobstep,'(i3)') js
             do i = 1, 3
                if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
             enddo
             F_OUT = "output"//name_jobstep//trim(cid)
             if(sw_wdir == ON) then
                inquire(file=trim(workdir)//F_OUT, exist=existence)
             else
                inquire(file=F_OUT, exist=existence)
             end if
!!$             inquire(file=F_OUT, exist=existence)

! ==== KT_add ====== 2014/07/20
#ifdef NEB_NEW_FILENAMES
             if ( .not. existence ) then
                f_tmp = trim(F_OUT) //'_conf0'
                inquire(file=f_tmp, exist=existence)
             endif
             if ( .not. existence ) then
                f_tmp = trim(F_OUT) // '_conf00'
                inquire(file=f_tmp, exist=existence)
             endif
             if ( .not. existence ) then
                f_tmp = trim(F_OUT) // '_conf000'
                inquire(file=f_tmp, exist=existence)
             endif
#endif
! ================== 2014/07/20

             if(.not.existence) goto 1001
          end do
1001      if(js >=1000) js = 999
       end if
       if(npes > 1) call mpi_bcast(js,1,mpi_integer,0,MPI_CommGroup,ierr)
       write(name_jobstep,'(i3)') js
       do i = 1, 3
          if(name_jobstep(i:i) == ' ') name_jobstep(i:i)='0'
       enddo
! =========== KT_mod ==== 2014/07/20
!       if(mype==0) then
!          F_OUT = "output"//name_jobstep
!       else
!          F_OUT = "output"//name_jobstep//'_'//name_mype(name_length-imax+1:name_length)
!       end if
!
       F_OUT_BASE = "output"//name_jobstep
       if(mype==0) then
          F_OUT = trim(F_OUT_BASE) // trim(cid)
       else
          F_OUT = trim(F_OUT_BASE) // trim(F_PARA_EXTENSION) // trim(cid)
       end if
! ======================= 2014/07/20
    else
       fout_is_given = .true.
       if(mype > 0) then
! =========== KT_mod ==== 2014/07/20
!          if(ls+imax+1 > stringlength_of_filenames) then
!             F_OUT = F_OUT(1:stringlength_of_filenames-(imax+1))//'_'//name_mype(name_length-imax+1:name_length)
!          else
!             F_OUT = F_OUT(1:ls)//'_'//name_mype(name_length-imax+1:name_length)
!          end if
!
          if(ls+imax+1 > stringlength_of_filenames) then
             F_OUT_BASE = F_OUT(1:stringlength_of_filenames-(imax+1))
          else
             F_OUT_BASE = F_OUT(1:ls)
          end if
          F_OUT = trim(F_OUT_BASE) // F_PARA_EXTENSION
! ======================= 2014/07/20
       end if
    end if

    endif

    if(conf_para)then
      call mpi_bcast(F_OUT,len(F_OUT),mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(name_jobstep,len(name_jobstep),mpi_character,0,mpi_comm_world,ierr)

! === KT_add === 2014/07/20
      call mpi_bcast(F_OUT_BASE,len(F_OUT_BASE),mpi_character,0,mpi_comm_world,ierr)
! ============== 2014/07/20

!      if(mype_conf/=0.and.mype==0) then
#ifndef NEB_NEW_FILENAMES
      if(mype_conf/=0) then
#endif
        itmp=int(log10(real(nrank_conf)))+1
        if(itmp>2)iketa=itmp
        write(cketa,*) iketa
        write(cid,'(i'//trim(adjustl(cketa))//'.'//trim(adjustl(cketa))//')') mype_conf

! ==== KT_mod === 2014/07/20
!        F_OUT = trim(F_OUT)//trim(adjustl(prefix))//cid
!
        F_CONF_EXTENSION = trim(adjustl(prefix))//cid
        F_OUT = trim(F_OUT) // F_CONF_EXTENSION
! ============== 2014/07/20

#ifndef NEB_NEW_FILENAMES
      endif
#endif

! === KT_mod === 2014/07/20
!     if(mype/=0) F_OUT = trim(F_OUT)//'_'//name_mype(name_length-imax+1:name_length)
!
      if (mype/=0) F_OUT = trim(F_OUT) // F_PARA_EXTENSION
! ============== 2014/07/20

    endif

    if(sw_wdir == ON) then
!       open(6, file=trim(workdir)//F_OUT, status='unknown', form='formatted')
       open(nfout, file=trim(workdir)//F_OUT, status='unknown', form='formatted')
    else
!       open(6, file=F_OUT, status='unknown', form='formatted')
       open(nfout, file=F_OUT, status='unknown', form='formatted')
    end if

!!$    if(printable) write(nfout,'(" length of F_OUT = ",i8)') ls
!!$    if(printable) write(nfout,'(" F_OUT = ",a32)') F_OUT

  end subroutine m_Files_open_standardout

  subroutine m_Files_open_nfvlc
    if(mype == 0) &
         & call open0(nfvlc,F_VLC,'F_VLC     ',unknown,unformatted&
         &     , check_file_name_on)
  end subroutine m_Files_open_nfvlc

  subroutine m_Files_open_nfcntn
    logical open
    inquire(unit = nfcntn, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfcntn,F_CNTN,'F_CNTN    ',unknown,formatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_nfcntn

  subroutine m_Files_open_nfcntn_bin
    logical open
    inquire(unit = nfcntn_bin, opened = open)
    if(open .and. printable) write(nfout,'("! nfcntn_bin is already open")')
    if(printable) &
         & write(nfout,'(" F_CNTN_BIN_in = ",a40," <<m_Files_open_nfcntn_bin>>")') F_CNTN_BIN_in

    if(F_CNTN_BIN_in_partitioned) then
       if(.not.open .and. mype == 0) then
          call open0(nfcntn_bin,F_CNTN_BIN_in,'F_CNTN_BIN',unknown,unformatted&
               &     ,check_file_name_on)
       else if(.not.open .and. mype /= 0) then
          call open0(nfcntn_bin,F_CNTN_BIN_in,'F_CNTN_BIN',unknown,unformatted&
               &     ,check_file_name_off)
       end if
    else
       if(.not.open .and. mype==0) &
            & call open0(nfcntn_bin,F_CNTN_BIN_in,'F_CNTN_BIN',unknown,unformatted&
            &     ,check_file_name_on)
    end if
  end subroutine m_Files_open_nfcntn_bin

  Subroutine m_Files_reopen_nfcntn
    logical open
    inquire(unit = nfcntn, opened = open)
    if(open .and. mype == 0) close(nfcntn)
    if(mype==0) &
         & call open0(nfcntn,F_CNTN,'F_CNTN    ',unknown,formatted&
         &     ,check_file_name_on)
  end subroutine m_Files_reopen_nfcntn

  subroutine m_Files_reopen_nfcntn_bin
    logical open
    inquire(unit = nfcntn_bin, opened = open)
!!$    if(open .and. printable) write(nfout,'("! nfcntn_bin is already open")')
    if(open) close(nfcntn_bin)
    if(F_CNTN_BIN_partitioned) then
       if(mype == 0) then
          call open0(nfcntn_bin, F_CNTN_BIN, 'F_CNTN_BIN',unknown, unformatted,check_file_name_on)
       else
          call open0(nfcntn_bin, F_CNTN_BIN, 'F_CNTN_BIN',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype == 0) call open0(nfcntn_bin,F_CNTN_BIN,'F_CNTN_BIN',unknown &
            &                         ,unformatted,check_file_name_on)
    end if
  end subroutine m_Files_reopen_nfcntn_bin

  subroutine m_Files_open_nfcntn_bin_stm
    logical open
    inquire(unit = nfcntn_bin_stm, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfcntn_bin_stm,F_CNTN_BIN_STM,'F_CNTN_STM',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_nfcntn_bin_stm

  subroutine m_Files_open_nfcntn_berry()
    logical open
    inquire(unit = nfcntn_berry, opened = open)
    if(open) write(6,'("! nfcntn_berry is already opened")')
    if(.not.open .and. mype==0) &
    & call open0(nfcntn_berry,F_CNTN_BERRY,'F_CNTN_BERRY',unknown,unformatted &
    &      ,check_file_name_off)
  end subroutine m_Files_open_nfcntn_berry

  subroutine m_Files_close_nfcntn_berry()
    if(mype==0) close(nfcntn_berry)
  end subroutine m_Files_close_nfcntn_berry

  logical function m_Files_nfcntn_bin_paw_exists()
    logical :: ex
    inquire(file=F_CNTN_BIN_PAW,exist=ex)
    m_Files_nfcntn_bin_paw_exists = ex
  end function m_Files_nfcntn_bin_paw_exists

  logical function m_Files_nfoccmat_exists()
    logical :: ex
    inquire(file=F_OCCMAT,exist=ex)
    m_Files_nfoccmat_exists = ex
  end function m_Files_nfoccmat_exists

  subroutine m_Files_open_nfcntn_bin_paw
    logical open
    inquire(unit = nfcntn_bin_paw, opened = open)
    if(open) write(6,'("! nfcntn_bin_paw is already open")')
    if(open) close(nfcntn_bin_paw)
    if(open) open=.false.
    if(.not.open .and. mype==0) &
         & call open0(nfcntn_bin_paw,F_CNTN_BIN_PAW,'F_CNTN_BIN_PAW',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_nfcntn_bin_paw

! ======================== KT_add ================= 13.0D
  subroutine m_Files_close_nfcntn_bin_paw()
    if(mype==0) close(nfcntn_bin_paw)
  end subroutine m_Files_close_nfcntn_bin_paw

  subroutine m_Files_reopen_nfcntn_bin_paw
    logical open
    inquire(unit = nfcntn_bin_paw, opened = open)
    if(open) close(nfcntn_bin_paw)

    if(.not.open .and. mype==0) &
         & call open0(nfcntn_bin_paw,F_CNTN_BIN_PAW,'F_CNTN_BIN_PAW',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine m_Files_reopen_nfcntn_bin_paw
! ================================================= 13.0D

  subroutine m_Files_open_nfstop
    if(mype == 0) call open0(nfstop,F_STOP,'F_STOP    '&
         &                  ,unknown,formatted,check_file_name_off)       ! MPI
  end subroutine m_Files_open_nfstop

  subroutine m_Files_close_nfstop
    if(mype == 0) close(nfstop,status='keep')  ! MPI
  end subroutine m_Files_close_nfstop

  subroutine m_Files_open_nfstatus
    integer :: ls
    integer :: itmp,iketa
    character(len=256) :: suf,cketa
    integer :: ierr
    if(mype == 0) then
!!$       write(nfout,'(" !! name_jobstep = ",a3)') name_jobstep
!!$       write(nfout,*) ' !! fout_is_given = ', fout_is_given
       if(fout_is_given) then
       else
!!$          if(.not.conf_para .or. conf_para.and.mype_conf==0) then
          ls = len(trim(F_STATUS))
          if(ls+3 > stringlength_of_filenames) then
             F_STATUS = F_STATUS(1:stringlength_of_filenames-3-iketa)//name_jobstep
          else
             F_STATUS = F_STATUS(1:ls)//name_jobstep
          end if
!!$          endif
          fout_is_given = .true.
!!$ '09.11.11
          if(conf_para)then
!!$            call mpi_bcast(F_STATUS,len(F_STATUS),mpi_character,0,mpi_comm_world,ierr)

#ifndef NEB_NEW_FILENAMES
             if(mype_conf/=0) then
#endif
                iketa=2
                itmp = int(log10(real(nrank_conf)))+1
                if(itmp>2)iketa=itmp
                write(cketa,*) iketa
                write(suf,'(i'//trim(adjustl(cketa))//'.'//trim(adjustl(cketa))//')') mype_conf
                F_STATUS = trim(F_STATUS)//trim(adjustl(prefix))//trim(adjustl(suf))
#ifndef NEB_NEW_FILENAMES
             endif
#endif
          endif
!!$ '09.11.11
       end if
!!$       open(nfstatus, file=F_STATUS, status='unknown',form='formatted')
       if(sw_wdir == ON) then
          open(nfstatus, file=trim(workdir)//F_STATUS, status='unknown',form='formatted')
       else
          open(nfstatus, file=F_STATUS, status='unknown',form='formatted')
       end if
    end if
  end subroutine m_Files_open_nfstatus

  subroutine m_Files_close_nfstatus
    if(mype == 0) close(nfstatus,status='keep')  ! MPI
  end subroutine m_Files_close_nfstatus

  subroutine m_Files_open_qpoint_files(way_ksample,nbztyp)
    integer, intent(in) :: way_ksample,nbztyp
    if(way_ksample == FILE .or. nbztyp < GENERAL) then
       call open0(nfqpoint,F_QPOINT,'F_QPOINT  ', old, formatted,check_file_name_on)
    end if
  end subroutine m_Files_open_qpoint_files

  subroutine m_Files_open_kpoint_files(way_ksample,nbztyp)
    integer, intent(in) :: way_ksample,nbztyp
!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
!!    call open0(nfmatbp, F_MATBP, 'F_MATBP   ',unknown,formatted,check_file_name_on)
!!    call open0(nfkpoint,F_KPOINT,'F_KPOINT  ', old,   formatted,check_file_name_on)
    if(way_ksample == FILE .or. nbztyp < GENERAL) then
       call open0(nfkpoint,F_KPOINT,'F_KPOINT  ', old,   formatted,check_file_name_on)
!  <-- modified by T. Yamasaki 4th May 2003
    end if
!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
!!$    if(nbztyp >= GENERAL) then
!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
!!       call open0(nfopgr, F_OPGR, 'F_OPGR    ', old,   formatted,check_file_name_on)
!!$       call open0(nfspg,  F_SPG,  'F_SPG     ', old,   formatted,check_file_name_on)
!!$       call open0(nfkpgn, F_KPGN, 'F_KPGN    ', old,   formatted,check_file_name_on)
!!!!!!!! modified by mizouchi@adv 2003.2.21 !!!!!
!!$    end if
  end subroutine m_Files_open_kpoint_files

  subroutine m_Files_open_kpoint_bin()
#ifdef FIX_ENDIANNESS
    if(mype==0) call open0(nfkpoint_bin, F_KPOINT_BIN, 'F_KPOINT_BIN     ',unknown, &
    &  unformatted,check_file_name_on,fix_endianness=.true.)
#else
    if(mype==0) call open0(nfkpoint_bin, F_KPOINT_BIN, 'F_KPOINT_BIN     ',unknown, &
    &  unformatted,check_file_name_on)
#endif
  end subroutine m_Files_open_kpoint_bin

  subroutine m_Files_close_kpoint_bin()
    if(mype==0) close(nfkpoint_bin)
  end subroutine m_Files_close_kpoint_bin

  subroutine m_Files_open_nfspg()
    call open0(nfspg,  F_SPG,  'F_SPG     ', old,   formatted,check_file_name_on)
  end subroutine m_Files_open_nfspg

  subroutine m_Files_open_nfkindex
    call open0(nfkindex, F_KINDEX, 'F_KINDEX  ',old,formatted,check_file_name_on)
  end subroutine m_Files_open_nfkindex

  subroutine set_filenumbers_of_NFSTM
    nfstm(-1) = nfchgu
    nfstm( 0) = nfchgs
    nfstm( 1) = nfchgo
  end subroutine set_filenumbers_of_NFSTM

  subroutine checkfilenumbers
    if(printable) then
       write(nfout,*) ' nfstm(-1) = ', nfstm(-1)
       write(nfout,*) ' nfstm( 0) = ', nfstm( 0)
       write(nfout,*) ' nfstm( 1) = ', nfstm( 1)
    end if
  end subroutine checkfilenumbers

!!$  subroutine m_Files_set_default_filenames_paramset
!!$    F_PRM       = './m_ArraySize_Parameters.f90'
!!$  end subroutine m_Files_set_default_filenames_paramset

  subroutine m_Files_set_default_filenames
    integer :: i
    character(len=stringlength_of_filenames) :: pppath
    F_file_names_data   = "./file_names.data"
!c---------- Input files.
    call set_filenumbers_of_NFSTM
!!$    if(.not.paramset) call checkfilenumbers
    F_INP       = "./nfinp.data"
    F_POT(1)    = "./pot.01"   !  fort.37'
    F_POT(2)    = "./pot.02"   !  fort.38'
    F_POT(3)    = "./pot.03"   !  fort.39'
    F_POT(4)    = "./pot.04"   !  fort.40'
    F_POT(5)    = "./pot.05"   !  fort.45'
    F_POT(6)    = "./pot.06"   !  fort.46'
    F_POT(7)    = "./pot.07"   !  fort.11'
    F_POT(8)    = "./pot.08"   !  fort.12'
    F_POT(9)    = "./pot.09"   !  fort.13'
    F_POT(10)   = "./pot.10"   !  fort.14'
    F_POT(11)   = "./pot.11"   !  fort.15'
    F_POT(12)   = "./pot.12"   !  fort.16'
    F_POT(13)   = "./pot.13"   !  fort.17'
    F_POT(14)   = "./pot.14"   !  fort.18'
    F_POT(15)   = "./pot.15"   !  fort.19'
    F_POT(16)   = "./pot.16"   !  fort.36'
    do i=1,16
    F_POT_ORG(i) = F_POT(i)
    enddo
    F_PKB       = "./vkb.data"
    F_PD        = "./vd.data"
    F_PPC       = "./vpc.data"
    F_STOP      = "./nfstop.data"
    F_CNST      = "./nfcnst.data"
!c---------- nbztype >= 100
    F_OPGR      = "./opgr.data"
    F_MATBP     = "./matrix.BP"
    F_KPOINT    = "./kpoint.data"
    F_QPOINT    = "./qpoint.data"
    F_KINDEX    = "./f.kp0"
!!!!!!!! added by mizouchi@adv 2003.2.21 !!!!!
    F_SPG       = "./bnprp4.i5"
    F_KPGN      = "./bnkpgn.i5"
!!!!!!!! added by mizouchi@adv 2003.2.21 !!!!!
!c---------- Output files.
    F_OTP       = "./nfotp.data"
    F_DYNM      = "./nfdynm.data"
    F_DYNM_CIF  = "./nfdynm.cif"
!!$    F_FORCE     = "./nfforce.data"
    F_EGRID     = "./nfegrid.data"
    F_LDOS      = "./nfldos.data"
    F_ENF       = "./nfefn.data"
!!$    F_GPT       = "./nfgpt.data"
    F_CHGT      = "./nfchgt.data"
    F_CHGO      = "./nfchgo.data"
    F_CHGU      = "./nfchgu.data"
!#ifdef _STMIMG_EVERY_TIME_
    F_CHGS      = "./nfchgs.data"
!#endif
    F_ENERG     = "./nfenergy.data"
    F_OUT       = ""
    F_STATUS    = "./jobstatus"
!c----------- Information files
!c      F_IENG     = "./nfieng.data"
!c----------- Continue file
    F_CNTN      = "./continue.data"
    F_CNTN_BIN  = "./continue_bin.data"
    F_ZAJ       = "./zaj.data"
    F_CNTN_BIN_PAW  = "./continue_bin_paw.data"
    F_VLC       = "./nfvlc.data"
    F_CNTN_BIN_STM = "./continue_bin_stm.data"

!c----------- Work files.
    F_WF(1)     = "./ftn95.data"
    F_WF(2)     = "./ftn96.data"
    F_WF(3)     = "./ftn97.data"
!c---------------------
    if(paramset) F_PRM       = "./m_ArraySize_Parameters.f90"
    F_DOS       = "./dos.data"
    F_CHR       = "./nfchr.data"
    F_WFk       = "./nfwfk.data"
    F_WANNIER   = "./nfwannier.data"
    F_CNTN_WAN  = "./nfcontinue_wannier.data"
    F_POT_WAN   = "./nfpotential_wannier.data"
    F_ELF       = "./nfelf.data"
    F_BERRY     = "./berry.data"
    F_EFFCHG    = "./effchg.data"
    F_FORCE     = "./force.data"
    F_MODE      = "./mode.data"
    F_EPSILON   = "./epsilon.data"
    F_EPSOUT    = "./eps.data"
    F_NLO       = "./nlo.data"         ! UVSOR
    F_MAGOPT    = "./magopt.data"      ! UVSOR
    F_EPSCONT   = "./eps_continue.data"! UVSOR
!-------------------- T. Hamada 2021.9.23 -------------------- ! UVSOR
    F_KPT_DATA  = "./kpt.data"
    F_KPTEPS_DATA = "./kpteps.data"
!------------------------------------------------------------- ! UVSOR
    F_PHDOS     = "./phdos.data"
    F_PHDOS     = "./phdos.data"
    F_OCCMAT    = "./occmat.data"
    F_STRFRC    = "./strfrc.data"
    F_PSTRN     = "./positron.cube"
    F_VELEC     = "./electron.cube"
    F_EPPAIR    = "./ep_pair.cube"
    F_EPPAIR2   = "./ep_pair2.cube"
    F_VELEC_GRAD= "./electron_grad.cube"
    F_PSICOEF   = "./psicoef.data"
    F_BAND_SYM_INPUT = "./band_sym_input.data"
    F_EFERMI    = "./nfefermi.data"
! ---------- mpifiletypes
    F_ZAJ_filetype = "unified"
    F_CHGT_filetype = "unified"
    F_CNTN_BIN_filetype = "unified"
    F_CNTN_filetype = "unified"
    F_CNTN_BIN_PAW_filetype = "unified"

    F_ZAJ_in      = F_ZAJ
    F_CNTN_in     = F_CNTN
    F_CNTN_BIN_in = F_CNTN_BIN
    F_CHGT_in     = F_CHGT

    F_ZAJ_bak      = F_ZAJ
    F_CNTN_bak     = F_CNTN
    F_CNTN_BIN_bak = F_CNTN_BIN
    F_CHGT_bak     = F_CHGT

    F_ZAJ_in_filetype = "nogiven"
    F_CHGT_in_filetype = "nogiven"
    F_CNTN_BIN_in_filetype = "nogiven"
    F_CNTN_in_filetype = "nogiven"

! ---------- nebfiles
    F_NEB_OUT  = "./output_neb"
    F_IMAGE(0)  = "./endpoint0.data"
    F_IMAGE(-1) = "./endpoint1.data"
    F_NEB_STOP = "./nfnebstop.data"
    F_NEB_CNTN = "./neb_continue.data"
    F_NEB_ENF = "./nfnebenf.data"
    F_NEB_DYNM = "./nfnebdynm.data"
    F_PATH = "./nfdynm.data"

! ---------- PAW continuation file
    F_CNTN_BIN_PAW          = "./continue_bin_paw.data"

! ======================== KT_add ================= 13.0D
    F_CNTN_BIN_PAW_in = F_CNTN_BIN_PAW
    F_CNTN_BIN_PAW_in_filetype = 'nogiven'
! ================================================= 13.0D

! ============================== Added by K. Tagami ======== 10.1
    F_LR_SPECTRA = "./spectrum.data"
! ========================================================== 10.1

    F_CNTN_BERRY = './continue_bin_berry.data'
    F_EKINDENS   = './ekindens_bin.data'

! ======== KT_add =========== 13.0S
    F_CORE_ENERGY_OUT = './core_energy.data'
    F_CORE_ENERGY_INITIAL = './core_energy.initial'
    F_CORE_ENERGY_FINAL   = './core_energy.final'
! =========================== 13.0S

! ====== KT_add === 13.1R
    F_EPS_PHONON = './nfeps_phonon.data'
    F_OPTICAL_COEFF = './nfoptical_coeff.data'
    F_RAMAN_SPECTRA = './raman_spectra.data'
! ================= 13.1R

! ====== KT_add ====== 13.1XI
    F_EXCITATION_SPECTRA = './excitation_spectra.data'
! ================= 13.1XI

! ====== KT_add === 2014/07/14/
    F_HYPERVEC = './hypervector.data'
! ================= 2014/07/14/

    F_WFk_Squared  = "./wfnsq.cube"
    F_WFk_IntegMoment   = "./wfn_integ_moment.data"
    F_WFK_ORB_PROJ    = "./wfn_orb_proj.data"
    F_WFK_LOCAL_DECOMP    = "./wfn_local_decomp.data"

    F_BAND_SPECTR_WGHT = "./nfband_spectr_wght.data"
    F_PORB_DENS_MAT = "./porb_density_matrix.data"
    F_PORB_ROT_MAT = "./porb_rot_matrix.data"

    F_MAG_CORE = './nfmag_core.data'
    F_ELECPOT_BIN = './elecpot_bin.data'
    F_ELECPOT_BIN_REF = './elecpot_bin.data.ref'

    F_LATCONST = "./nflatconst.data"
    F_METRIC   = "./nfmetric.data"

    F_DFTD3PAR = "./dftd3par.data"

    F_CHKPNT   = "chkpnt"

    F_HYBRIDINFO   = "./nfhybridinfo.data"
    F_SCF_ZAJ   = "./zaj.data"

    F_POS = "./nfinp.data"
    F_COORD_ATTR = "./nfinp.data"

    call getenv('PHASE_PP_PATH',pppath)
    if(len(pppath)==0) then
    F_DEFAULTPP = './defaultpp.data'
    else
    F_DEFAULTPP = trim(pppath)//'/defaultpp.data'
    endif

    F_PWBS = './nfpwbs.data'

    F_XSF  = './nfdynm.xsf'
    F_N2P2 = './nfdynm.n2p2'

    F_INP_MOD = './nfinp_mod.data'

    F_ZAJ_KALL = './zaj_kall.data'

    F_KPOINT_BIN = './kpoint_bin.data'

  end subroutine m_Files_set_default_filenames

  subroutine m_Files_set_def_fname_pos(fname)
    character(len=*), intent(in) :: fname
    F_POS = fname
  end subroutine m_Files_set_def_fname_pos

  subroutine m_Files_set_def_fname_attr()
    F_COORD_ATTR = F_INP
  end subroutine m_Files_set_def_fname_attr

  subroutine m_Files_rd_file_names_data
    integer, parameter ::  nffile = 99
    logical :: existence
    integer :: iwarning

    if(sw_wdir==ON)then
       inquire(file=trim(workdir)//F_file_names_data, exist = existence)
    else
       inquire(file=F_file_names_data, exist = existence)
    endif
    if(existence) then
       call open0(nffile, F_file_names_data, "FFILENAMES", old, formatted,check_file_name_off)
       !      open(nffile,file="./file_names.data",status='unknown')
       rewind nffile
       read(nffile,NML = fnames, err = 1007, end = 1004)
1004   continue
       rewind nffile
       read(nffile,NML = nebfiles, err = 1008, end = 1005)
1005   continue
       close(nffile,status='keep')
    else
       iwarning = FILENAMES_NOT_EXIST
       if(mype==0) call m_EMsg_Warning(iwarning, nfout)
    end if

    F_ZAJ_bak      = F_ZAJ
    F_CHGT_bak     = F_CHGT
    F_CNTN_BIN_bak = F_CNTN_BIN
    F_CNTN_bak     = F_CNTN

    return
!1007 stop ' file_names.data fnames format error'
!1008 stop ' file_names.data nebfiles format error'
1007 call phase_execution_error(FILENAMES_FORMAT_ERROR)
1008 call phase_execution_error(FILENAMES_FORMAT_ERROR_NEB)
  end subroutine m_Files_rd_file_names_data

  logical function m_Files_check_nfcntn_existence()
    logical :: existence
    existence = .false.
!!$    if((.not.F_CNTN_in_partitioned .and. mype == 0) .or. F_CNTN_in_partitioned) then
    if(mype == 0) then
       if(sw_wdir == ON) then
          inquire(file=trim(workdir)//F_CNTN_in,exist = existence)
       else
          inquire(file=F_CNTN_in,exist = existence)
       end if
    end if
!!$    if(.not.F_CNTN_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
    if(npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
    m_Files_check_nfcntn_existence = existence
  end function m_Files_check_nfcntn_existence

  subroutine m_Files_check_file_existence
    integer, parameter ::  nffile = 99
    logical :: existence, existence_all
    integer, allocatable, dimension(:,:) :: file_existence !d(npes,cont files)
    integer, allocatable, dimension(:) :: file_existence_local !d(cont files)
    integer :: i, j, npartitioned_files
    character :: name_mype*12
!!!    logical :: existence
    integer :: name_length, imax
    logical, save :: firstcall_to_add_name_mype=.true.
    if(sw_wdir==ON)then
       inquire(file=trim(workdir)//F_file_names_data, exist = existence)
    else
       inquire(file=F_file_names_data, exist = existence)
    endif
    if(existence) then
       call open0(nffile, F_file_names_data, "FFILENAMES", old, formatted,check_file_name_on)
       rewind nffile
       read(nffile,NML = mpifiletypes, err = 1005, end = 1005)
    end if

    if(F_ZAJ_in_filetype == "nogiven") F_ZAJ_in_filetype = F_ZAJ_filetype
    if(F_CHGT_in_filetype == "nogiven") F_CHGT_in_filetype = F_CHGT_filetype
    if(F_CNTN_BIN_in_filetype == "nogiven") F_CNTN_BIN_in_filetype = F_CNTN_BIN_filetype
    if(F_CNTN_in_filetype == "nogiven") F_CNTN_in_filetype = F_CNTN_filetype

! ======================== KT_add ================= 13.0D
    if(F_CNTN_BIN_PAW_in_filetype == 'nogiven') &
         &            F_CNTN_BIN_PAW_in_filetype = F_CNTN_BIN_PAW_filetype
! ================================================= 13.0D

    if(printable) then
       write(nfout,'(" !! mpifiletypes is read")')
       write(nfout,'(" !! F_ZAJ_filetype         = ",a12)') F_ZAJ_filetype
       write(nfout,'(" !! F_CHGT_filetype        = ",a12)') F_CHGT_filetype
       write(nfout,'(" !! F_CNTN_BIN_filetype    = ",a12)') F_CNTN_BIN_filetype
       write(nfout,'(" !! F_CNTN_filetype        = ",a12)') F_CNTN_filetype
       write(nfout,'(" !! F_ZAJ_in_filetype      = ",a12)') F_ZAJ_in_filetype
       write(nfout,'(" !! F_CHGT_in_filetype     = ",a12)') F_CHGT_in_filetype
       write(nfout,'(" !! F_CNTN_BIN_in_filetype = ",a12)') F_CNTN_BIN_in_filetype
       write(nfout,'(" !! F_CNTN_in_filetype     = ",a12)') F_CNTN_in_filetype

! ======================== KT_add ================= 13.0D
       write(nfout,'(" !! F_CNTN_BIN_PAW_in_filetype     = ",a12)') &
            &                            F_CNTN_BIN_PAW_in_filetype
! ================================================= 13.0D
    end if

    call strncmp0('partitioned',trim(F_ZAJ_filetype),     F_ZAJ_partitioned)
    call strncmp0('partitioned',trim(F_CHGT_filetype),    F_CHGT_partitioned)
    call strncmp0('partitioned',trim(F_CNTN_BIN_filetype),F_CNTN_BIN_partitioned)
    call strncmp0('partitioned',trim(F_CNTN_filetype),    F_CNTN_partitioned)
    call strncmp0('partitioned',trim(F_ZAJ_in_filetype),     F_ZAJ_in_partitioned)
    call strncmp0('partitioned',trim(F_CHGT_in_filetype),    F_CHGT_in_partitioned)
    call strncmp0('partitioned',trim(F_CNTN_BIN_in_filetype),F_CNTN_BIN_in_partitioned)
    call strncmp0('partitioned',trim(F_CNTN_in_filetype),    F_CNTN_in_partitioned)

! ======================== KT_add ================= 13.0D
    call strncmp0('partitioned',trim(F_CNTN_BIN_PAW_in_filetype),&
         &                                  F_CNTN_BIN_PAW_in_partitioned )
! ================================================= 13.0D

    if(npes > 1) then
       if(printable) then
          write(nfout,*) '!! F_ZAJ_partitioned         = ',F_ZAJ_partitioned
          write(nfout,*) '!! F_CHGT_partitioned        = ',F_CHGT_partitioned
          write(nfout,*) '!! F_CNTN_BIN_partitioned    = ',F_CNTN_BIN_partitioned
          write(nfout,*) '!! F_CNTN_partitioned        = ',F_CNTN_partitioned
          write(nfout,*) '!! F_ZAJ_in_partitioned      = ',F_ZAJ_in_partitioned
          write(nfout,*) '!! F_CHGT_in_partitioned     = ',F_CHGT_in_partitioned
          write(nfout,*) '!! F_CNTN_BIN_in_partitioned = ',F_CNTN_BIN_in_partitioned
          write(nfout,*) '!! F_CNTN_in_partitioned     = ',F_CNTN_in_partitioned
       end if
    end if

    goto 1006
1005 continue
    if(printable) write(nfout,'(" !! mpifiletypes is not read")')
1006 continue

    F_ZAJ_in      = F_ZAJ
    F_CHGT_in     = F_CHGT
    F_CNTN_BIN_in = F_CNTN_BIN
    F_CNTN_in     = F_CNTN

! ======================== KT_add ================= 13.0D
    F_CNTN_BIN_PAW_in = F_CNTN_BIN_PAW
! ================================================= 13.0D

    if(npes > 1) then
       name_length = 12
       imax = int(log10(dble(npes-1)))+1
       if(imax < min_ext_length_mype) imax = min_ext_length_mype
       if(imax > name_length) then
          call phase_error_with_msg(nfout,' number of pes (=npes) is larger than 1,000,000,000,000',__LINE__,__FILE__)
       end if
       write(name_mype,'(i12)') mype
       do i = 1+name_length-imax, name_length
          if(name_mype(i:i) == ' ') name_mype(i:i)='0'
       enddo
!!$       do i = 1, 4
!!$          if(name_mype(i:i) == ' ') name_mype(i:i) = '0'
!!$       end do
       if(firstcall_to_add_name_mype) then
         if(F_ZAJ_partitioned)      call add_name_mype(F_ZAJ)
         if(F_CHGT_partitioned)     call add_name_mype(F_CHGT)
         if(F_CNTN_BIN_partitioned) call add_name_mype(F_CNTN_BIN)
         if(F_CNTN_partitioned)     call add_name_mype(F_CNTN)
         if(F_ZAJ_in_partitioned)      call add_name_mype(F_ZAJ_in)
         if(F_CHGT_in_partitioned)     call add_name_mype(F_CHGT_in)
         if(F_CNTN_BIN_in_partitioned) call add_name_mype(F_CNTN_BIN_in)
         if(F_CNTN_in_partitioned)     call add_name_mype(F_CNTN_in)

! ======================== KT_add ================= 13.0D
         if (F_CNTN_BIN_PAW_in_partitioned) call add_name_mype(F_CNTN_BIN_PAW_in)
! ================================================= 13.0D
         firstcall_to_add_name_mype = .false.
       endif

    end if

    if(printable) write(nfout,'(" F_POT(1) = ",a60)') F_POT(1)
    !!$ print *, ' paramset = ',paramset

    if(paramset) call read_file_names_data_paramset

!!$    file_existence_contfiles = .true.
!!$    file_existence_3contfiles = .true.

!!$    inquire(file=F_CNTN_in,exist = existence)
    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_CNTN_in,exist = existence)
    else
       inquire(file=F_CNTN_in,exist = existence)
    end if

    if(.not.F_CNTN_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
!!345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!!$    if(.not.existence) file_existence_contfiles = .false.
!!$    if(.not.existence) file_existence_3contfiles = .false.
    file_existence_nfcntn = existence

!!$    inquire(file=F_CNTN_BIN_in, exist = existence)
    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_CNTN_BIN_in, exist = existence)
    else
       inquire(file=F_CNTN_BIN_in, exist = existence)
    end if
    if(.not.F_CNTN_BIN_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$    if(.not.existence) file_existence_contfiles = .false.
    file_existence_nfcntn_bin = existence

!   -------------------------------  ! T. Y. 2020/02/28
    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_CNTN_BIN_PAW_in,exist = existence)
    else
       inquire(file=F_CNTN_BIN_PAW_in,exist = existence)
    end if

    if(.not.F_CNTN_BIN_PAW_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
    file_existence_nfcntn_bin_paw = existence
!   -------------------------------

    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_ZAJ, exist = existence)
    else
       inquire(file=F_ZAJ_in, exist = existence)
    end if
    if(.not.F_ZAJ_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$    if(.not.existence) file_existence_contfiles = .false.
!!$    if(.not.existence) file_existence_3contfiles = .false.
    file_existence_nfzaj = existence

!!$    inquire(file=F_CHGT_in, exist = existence)
    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_CHGT, exist = existence)
    else
       inquire(file=F_CHGT_in, exist = existence)
    end if
    if(.not.F_CHGT_in_partitioned .and. npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
!!$    if(.not.existence) file_existence_contfiles = .false.
!!$    if(.not.existence) file_existence_3contfiles = .false.
    file_existence_nfchgt = existence

!   -------------------------------
    if(sw_wdir == ON) then
       inquire(file=trim(workdir)//F_OCCMAT,exist = existence)
    else
       inquire(file=F_OCCMAT,exist = existence)
    end if

    if(npes >1) call mpi_bcast(existence,1,mpi_logical,0,MPI_CommGroup,ierr)
    file_existence_nfoccmat = existence
!   -------------------------------

    if(npes > 1) then
       npartitioned_files = 0
       if(F_CNTN_in_partitioned)         npartitioned_files = npartitioned_files+1
       if(F_CNTN_BIN_in_partitioned)     npartitioned_files = npartitioned_files+1
       if(F_CNTN_BIN_PAW_in_partitioned) npartitioned_files = npartitioned_files+1
                                                             ! T. Y. 2020/02/28
       if(F_ZAJ_in_partitioned)          npartitioned_files = npartitioned_files+1
       if(F_CHGT_in_partitioned)         npartitioned_files = npartitioned_files+1
       if(printable) write(nfout,'(" !! npartitioned_files = ",i5)') npartitioned_files
       if(npartitioned_files >=1) then
          allocate(file_existence(0:npes-1,npartitioned_files))
          allocate(file_existence_local(npartitioned_files))
          file_existence_local = 0
          j = 0
          if(F_CNTN_in_partitioned ) then
             j = j+1; if(file_existence_nfcntn)     file_existence_local(j) = 1
          end if
          if(F_CNTN_BIN_in_partitioned) then
             j = j+1; if(file_existence_nfcntn_bin) file_existence_local(j) = 1
          end if
          if(F_CNTN_BIN_PAW_in_partitioned) then
             j = j+1; if(file_existence_nfcntn_bin_paw) file_existence_local(j) = 1
          end if
          if(F_ZAJ_in_partitioned) then
             j = j+1; if(file_existence_nfzaj)      file_existence_local(j) = 1
          end if
          if(F_CHGT_in_partitioned) then
             j = j+1; if(file_existence_nfchgt)     file_existence_local(j) = 1
          end if

#ifdef FFT_ALLTOALL
          call mpi_gather(file_existence_local, npartitioned_files, mpi_integer, &
                          file_existence, npartitioned_files, mpi_integer, &
                          0, MPI_CommGroup, ierr)
#else
          do i = 0, npes-1
             if(i >= 1) then
                if(mype == i ) then
                   call mpi_send(file_existence_local,npartitioned_files,mpi_integer,0,1,MPI_CommGroup,ierr)
                else if(mype == 0) then
                   call mpi_recv(file_existence_local,npartitioned_files,mpi_integer,i,1,MPI_CommGroup,istatus,ierr)
                end if
             end if
             if(mype == 0) file_existence(i,1:npartitioned_files) = file_existence_local(1:npartitioned_files)
         end do
#endif
       end if

       if(mype == 0) then
          j = 0
          if(F_CNTN_in_partitioned) then
             j = j+1; call check_file_existence_all(j,existence_all); file_existence_nfcntn = existence_all
          end if
          if(F_CNTN_BIN_in_partitioned) then
             j = j+1; call check_file_existence_all(j,existence_all); file_existence_nfcntn_bin = existence_all
          end if
          if(F_CNTN_BIN_PAW_in_partitioned) then
             j = j+1; call check_file_existence_all(j,existence_all); file_existence_nfcntn_bin_paw = existence_all
          end if
          if(F_ZAJ_in_partitioned) then
             j = j+1; call check_file_existence_all(j,existence_all); file_existence_nfzaj = existence_all
          end if
          if(F_CHGT_in_partitioned) then
             j = j+1; call check_file_existence_all(j,existence_all); file_existence_nfchgt = existence_all
          end if
       end if

       if(npartitioned_files >= 1) then
          deallocate(file_existence_local)
          deallocate(file_existence)
       end if
    end if

    file_existence_contfiles = .true.
    file_existence_contfiles_when_paw_on = .true.
    if(.not.file_existence_nfcntn)     file_existence_contfiles = .false.
    if(printable) write(nfout,'(" file_existence_nfcntn     = ",L3," file_existence_contfiles = ",L3)') &
         & file_existence_nfcntn, file_existence_contfiles
    if(.not.file_existence_nfcntn_bin) file_existence_contfiles = .false.
    if(printable) write(nfout,'(" file_existence_nfcntn_bin = ",L3," file_existence_contfiles = ",L3)') &
         & file_existence_nfcntn_bin, file_existence_contfiles
    if(.not.file_existence_nfzaj)      file_existence_contfiles = .false.
    if(printable) write(nfout,'(" file_existence_nfzaj      = ",L3," file_existence_contfiles = ",L3)') &
         & file_existence_nfzaj, file_existence_contfiles
    if(.not.file_existence_nfchgt)     file_existence_contfiles = .false.
    if(printable) write(nfout,'(" file_existence_nfchgt     = ",L3," file_existence_contfiles = ",L3)') &
         & file_existence_nfchgt, file_existence_contfiles

    if(.not.file_existence_contfiles .or. .not.file_existence_nfcntn_bin_paw) &
         &   file_existence_contfiles_when_paw_on = .false.    ! T. Y. 2020/02/28
    if(printable)  write(nfout,&
         &  '(" file_existence_nfcntn_bin_paw = ",L3," file_existence_contfiles_when_paw_on = ",L3)') &
         &      file_existence_nfcntn_bin_paw,         file_existence_contfiles_when_paw_on  ! T. Y. 2020/02/28

    file_existence_3contfiles = .true.
    if(.not.file_existence_nfcntn)     file_existence_3contfiles = .false.
    if(.not.file_existence_nfzaj)      file_existence_3contfiles = .false.
    if(.not.file_existence_nfchgt)     file_existence_3contfiles = .false.

    if(printable) then
       write(nfout,'(" --- existence check of continue files ---")')
       call check_file_existence('F_CNTN          ',F_CNTN_in,file_existence_nfcntn)
       call check_file_existence('F_CNTN_BIN      ',F_CNTN_BIN_in,file_existence_nfcntn_bin)
       call check_file_existence('F_CNTN_BIN_PAW  ',F_CNTN_BIN_PAW_in,file_existence_nfcntn_bin_paw) ! T. Y. 2020/02/28
       call check_file_existence('F_ZAJ           ',F_ZAJ_in,file_existence_nfzaj)
       call check_file_existence('F_CHGT          ',F_CHGT_in,file_existence_nfchgt)
       call check_file_existence('contfiles       ','continue files',file_existence_contfiles)
       call check_file_existence('contfiles paw on','continue files when paw on' & ! T. Y. 2020/02/28
            &                    ,file_existence_contfiles_when_paw_on)
       call check_file_existence('3contfiles      ','3continue files',file_existence_3contfiles)
       write(nfout,'(" -----------------------------------------")')
    end if
    if(existence) close(nffile,status='keep')
  contains
    subroutine add_name_mype(F)
      character(len=*), intent(inout) :: F
      integer :: ls
      ls = len(trim(F))
!!$      if(ls+5 > stringlength_of_filenames) then
      if(ls+imax+1 > stringlength_of_filenames) then
         F = F(1:stringlength_of_filenames-(imax+1))//'_'//name_mype(name_length-imax+1:name_length)
      else
         F = F(1:ls)//'_'//name_mype(name_length-imax+1:name_length)
      end if
    end subroutine add_name_mype

    subroutine check_file_existence_all(npfile,tf)
      integer, intent(in)  :: npfile
      logical, intent(out) :: tf
      integer :: i
      tf = .true.
      do i = 0, npes-1
         if(file_existence(i,npfile) == 0) then
            tf = .false.
            exit
         end if
      end do
    end subroutine check_file_existence_all

    subroutine check_file_existence(tagname,filename,torf)
      character(len=16),intent(in) :: tagname
      character(len=*), intent(in) :: filename
      logical, intent(in) ::          torf
      integer :: ls
      character(len=12) :: exist_or_not

      if(torf) then
         exist_or_not = 'existing'
      else
         exist_or_not = 'not existing'
      end if

      ls = len(trim(filename))
      write(nfout,*) tagname," (= ",trim(filename)," ) ", trim(exist_or_not)
    end subroutine check_file_existence
  end subroutine m_Files_check_file_existence

  function m_Files_check_nfzaj_existence()
    logical m_Files_check_nfzaj_existence
    m_Files_check_nfzaj_existence = file_existence_nfzaj
  end function m_Files_check_nfzaj_existence

  function m_Files_check_nfchgt_existence()
    logical m_Files_check_nfchgt_existence
    m_Files_check_nfchgt_existence = file_existence_nfchgt
  end function m_Files_check_nfchgt_existence
  subroutine m_Files_check_file_names()
    if(printable) then
       write(nfout,'(" --- check of file_names ---")')
       call repeat_filename('F_INP     ',F_INP)
       call repeat_filename('F_POT(1)  ',F_POT(1))
       call repeat_filename('F_POT(2)  ',F_POT(2))
       call repeat_filename('F_POT(3)  ',F_POT(3))
!!$       write(nfout,'(" F_INP      = ", a60)') F_INP
!!$       write(nfout,'(" F_POT(1)   = ", a60)') F_POT(1)
!!$       write(nfout,'(" F_POT(2)   = ", a60)') F_POT(2)
!!$       write(nfout,'(" F_POT(3)   = ", a60)') F_POT(3)
!!$       call repeat_filename('F_CNST    ',F_CNST)
       call repeat_filename('F_KPOINT  ',F_KPOINT)
       call repeat_filename('F_CHGT    ',F_CHGT_in)
       call repeat_filename('F_CNTN    ',F_CNTN_in)
       call repeat_filename('F_CNTN_BIN',F_CNTN_BIN_in)
       call repeat_filename('F_ZAJ     ',F_ZAJ_in)
       call repeat_filename('F_STOP    ',F_STOP)
!!$       write(nfout,'(" F_CNST     = ", a60)') F_CNST
!!$       write(nfout,'(" F_CHGT     = ", a60)') F_CHGT
!!$!!$    write(nfout,'(" F_MATBP    = ", a60)') F_MATBP
!!$       write(nfout,'(" F_KPOINT   = ", a60)') F_KPOINT
!!$!!$       write(nfout,'(" F_KINDEX   = ", a60)') F_KINDEX
!!$       write(nfout,'(" F_CNTN     = ", a60)') F_CNTN
!!$       write(nfout,'(" F_CNTN_BIN = ", a60)') F_CNTN_BIN
!!$       write(nfout,'(" F_ZAJ      = ", a60)') F_ZAJ

!!$       if(ekmode==ON) write(nfout,'(" F_ENERG    = ", a60)') F_ENERG
       if(ekmode==ON) call repeat_filename('F_ENERG   ',F_ENERG)

       write(nfout,'(" F_NEB_OUT   = ", a60)') F_NEB_OUT
       write(nfout,'(" F_IMAGE(0)  = ", a60)') F_IMAGE(0)
       write(nfout,'(" F_IMAGE(-1) = ", a60)') F_IMAGE(-1)

! ======================== KT_mod ================= 13.0D
!!       write(nfout,'(" F_CNTN_BIN_PAW = ", a60)') F_CNTN_BIN_PAW
       call repeat_filename('F_CNTN_BIN_PAW',F_CNTN_BIN_PAW_in)
! ================================================= 13.0D

    end if
  contains
    subroutine repeat_filename(tagname,filename)
      character(len=*), intent(in) :: tagname
      character(len=*),  intent(in) :: filename
      integer :: ls
      character(len=10) :: tagname2
      ls = len(trim(tagname))
      if(ls < 10) then
         tagname2 = "          "
         tagname2 = trim(tagname)
         write(nfout,'(1x,a10,a3,a)') tagname2, " = ", trim(filename)
      else
         write(nfout,*) trim(tagname)," = ", trim(filename)
      end if
    end subroutine repeat_filename
  end subroutine m_Files_check_file_names

  subroutine read_file_names_data_paramset
    logical ::   existence
    integer, parameter :: nffile = 5
    if(sw_wdir==ON)then
      inquire(file=trim(workdir)//F_file_names_data, exist = existence)
    else
      inquire(file=F_file_names_data, exist = existence)
    endif
    if(existence) then
       call open0(nffile,F_file_names_data,'FFILENAMES',old,formatted,check_file_name_off)
       rewind nffile
       read(nffile,NML = f_param_name, err = 1004, end = 1004)
1004   continue
       close(nffile,status='keep')
    end if

  end subroutine read_file_names_data_paramset

  subroutine m_Files_echo_nfinp()
     character(len=1024) :: buf
     if(printable)then
        write(nfout,'(a)')
        write(nfout,'(a)') ' !** contents of the input parameter file : '//trim(F_INP)
        if(sw_wdir==ON)then
           open(nfinp,file=trim(workdir)//trim(F_INP),status='old')
        else
           open(nfinp,file=trim(F_INP),status='old')
        endif
        do
           read(nfinp,'(a)',end=1000,err=1000) buf
           write(nfout,'(a)') trim(buf)
        enddo
1000    continue
        write(nfout,'(a)') ' !**'
        close(nfinp)
     endif
  end subroutine m_Files_echo_nfinp

  subroutine m_Files_swap_F_INP(mode,unum,fname)
    integer, intent(in) :: mode
    integer, intent(in) :: unum
    character(len=stringlength_of_filenames),intent(in) :: fname
    integer :: iret, f_closeInputFile, f_openInputFile, ierror
    logical :: open
    if(mode==0)then
      inquire(unit = nfinp, opened = open)
      if(open) close(nfinp,status='keep')
      iret = f_closeInputFile()
      call m_Files_open_finp(1,unum,fname)
    else if (mode==1)then
      inquire(unit = unum, opened = open)
      if(open) close(unum,status='keep')
      iret = f_closeInputFile()
      call m_Files_open_finp(1,nfinp,F_INP)
    endif
  end subroutine m_Files_swap_F_INP

  subroutine m_Files_open_finp(mode,unum,fname, verbose)
    integer, intent(in) :: mode
    integer, intent(in) :: unum
    character(len=stringlength_of_filenames),intent(in) :: fname
    logical, intent(in), optional :: verbose
    logical :: open
    integer :: iret, f_closeInputFile, f_openInputFile, ierror
    logical :: verb
    integer :: check_fname
    verb = .true.
    if(present(verbose)) verb = verbose
    if(mode == 1) then
       inquire(unit = unum, opened = open)
       if(open) close(unum,status='keep')
       if(printable.and.verb) write(nfout,'(" !!  F_INP = ",a32)') fname
!!$       iret = f_openInputFile('./testinp.data')
!!$       iret = f_openInputFile(F_INP)
       if(sw_wdir == ON) then
          iret = f_openInputFile(trim(workdir)//fname)
       else
          iret = f_openInputFile(fname)
       end if
       if( iret < 0) then
          if(printable) write(nfout,'(" !!! Error in opening of inputfile: ",a32)') fname
!!$          stop ' stop at << m_Files_reopen_nfinp >>'
       else if(iret > 0) then
          if(printable) write(nfout,'(" !!! There is something wrong in the input file: ",a32)') fname
       end if
       if(iret/=0) then
          ierror = ERROR_IN_INPUTFILE_OPENING
          if(sw_wdir == ON) then
             call phase_error(ierror,nfout,unum,trim(workdir)//fname)
          else
             call phase_error(ierror,nfout,unum,fname)
          end if
       end if
    else
       iret = f_closeInputFile()
       check_fname = check_file_name_on
       if(.not.verb) check_fname = 0
       call open0(unum, fname, 'F_INP     ',    old,   formatted,check_fname)
    end if
  end subroutine m_Files_open_finp

  subroutine m_Files_reopen_nfinp(mode)
    integer, intent(in) :: mode
    logical :: open
    integer :: iret, f_closeInputFile, f_openInputFile, ierror
    call m_Files_open_finp(mode,nfinp,F_INP)
  end subroutine m_Files_reopen_nfinp

  subroutine m_Files_reopen_nfinp_mod(mode)
    integer, intent(in) :: mode
    logical :: open
    integer :: iret, f_closeInputFile, f_openInputFile, ierror
    call m_Files_open_finp(mode, nfinp_mod, F_INP_MOD, verbose=.false.)
  end subroutine m_Files_reopen_nfinp_mod

  subroutine m_Files_close_nfinp_mod()
    logical :: op
    if(mype == 0) then
       inquire(unit = nfinp_mod, opened = op)
       if(op) close(nfinp_mod)
    end if
  end subroutine m_Files_close_nfinp_mod

#ifndef _EMPIRICAL_
  subroutine m_Files_reopen_nfzaj
    logical :: open
    inquire(unit=nfzaj,opened=open)
    if(open) close(nfzaj,status='keep')

    if(F_ZAJ_partitioned) then
       if(mype == 0) then
          call open0(nfzaj, F_ZAJ, 'F_ZAJ     ',unknown, unformatted,check_file_name_on)
       else
          call open0(nfzaj, F_ZAJ, 'F_ZAJ     ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) call open0(nfzaj, F_ZAJ, 'F_ZAJ     ',unknown, unformatted,check_file_name_on)
    end if
  end subroutine m_Files_reopen_nfzaj

  subroutine m_Files_close_nfzaj()
    if(F_ZAJ_partitioned) then
       close(nfzaj)
    else
       if(mype==0) close(nfzaj)
    end if
  end subroutine m_Files_close_nfzaj

  subroutine m_Files_reopen_nfzaj_kall
    logical :: open
    inquire(unit=nfzaj_kall,opened=open)
    if(open) close(nfzaj_kall,status='keep')
    if(mype==0) call open0(nfzaj_kall, F_ZAJ_KALL, 'F_ZAJ_KALL',unknown, unformatted,check_file_name_on)
  end subroutine m_Files_reopen_nfzaj_kall

  subroutine m_Files_close_nfzaj_kall()
    if(mype==0) close(nfzaj_kall)
  end subroutine m_Files_close_nfzaj_kall

  subroutine m_Files_reopen_nfchgt()
    logical :: open
    inquire(unit=nfchgt,opened=open)
    if(open) close(nfchgt,status='keep')

    if(F_CHGT_partitioned) then
       if(mype == 0) then
          call open0(nfchgt, F_CHGT, 'F_CHGT    ',unknown, unformatted,check_file_name_on)
       else
          call open0(nfchgt, F_CHGT, 'F_CHGT    ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) call open0(nfchgt, F_CHGT, 'F_CHGT    ',unknown, unformatted,check_file_name_on)
    end if
  end subroutine m_Files_reopen_nfchgt

#endif

  subroutine m_Files_open_nfzaj()
    if(F_ZAJ_in_partitioned) then
       if(mype == 0) then
          call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
       else
          call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfzaj

  subroutine m_Files_open_nfzaj_kall()
    if(mype==0) call open0(nfzaj_kall, F_ZAJ_KALL, 'F_ZAJ_KALL  ',unknown, unformatted,check_file_name_on)
  end subroutine m_Files_open_nfzaj_kall

  subroutine m_Files_open_nfzaj_append_kall()
    if(mype==0) call open1(nfzaj_kall, F_ZAJ_KALL, 'F_ZAJ_KALL  ',unknown, unformatted,check_file_name_on)
  end subroutine m_Files_open_nfzaj_append_kall

  subroutine m_Files_open_nfzaj_with_check()
    logical :: open
    if(F_ZAJ_in_partitioned) then
       if(mype == 0) then
          inquire(unit=nfzaj,opened=open)
          if(.not.open) call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
       else
          inquire(unit=nfzaj,opened=open)
          if(.not.open) call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) then
          inquire(unit=nfzaj,opened=open)
          if(.not.open) call open0(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
       end if
    end if
! ================================ modified by K. Tagami =========== 11.0
!    if(open .and.mype==0) then
    if ( mype==0 .and. open ) then
! =================================================================== 11.0
       write(nfout,'(" nfzaj is opened")')
    end if
  end subroutine m_Files_open_nfzaj_with_check

  subroutine m_Files_open_nfzaj_append()
    if(F_ZAJ_in_partitioned) then
       if(mype == 0) then
          call open1(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
       else
          call open1(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) call open1(nfzaj, F_ZAJ_in, 'F_ZAJ_in  ',unknown, unformatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfzaj_append

  subroutine m_Files_open_nfchgt()
    if(F_CHGT_in_partitioned) then
       if(mype == 0) then
          call open0(nfchgt,F_CHGT_in,'F_CHGT_in ',unknown, unformatted,check_file_name_on)
       else
          call open0(nfchgt,F_CHGT_in,'F_CHGT_in ',unknown, unformatted,check_file_name_off)
       end if
    else
       if(mype==0) call open0(nfchgt,F_CHGT_in,'F_CHGT_in ',unknown, unformatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfchgt

  subroutine m_Files_open_files_initially
    call open0(nfinp, F_INP, 'F_INP     ',    old,   formatted,check_file_name_on)

    if(ekmode==ON) then
       if(mype==0) call open0(nfeng, F_ENERG,'F_ENERG   ',unknown,formatted,check_file_name_on)
    else
       if(file_existence_contfiles) then
          if(mype==0) call open1(nfdynm,F_DYNM,'F_DYNM    ',unknown,   formatted,check_file_name_on)
          if(mype==0) call open1(nfenf, F_ENF, 'F_ENF     ',unknown,   formatted,check_file_name_on)
       else
          if(mype==0) call open0(nfdynm,F_DYNM,'F_DYNM    ',unknown,   formatted,check_file_name_on)
          if(mype==0) call open0(nfenf, F_ENF, 'F_ENF     ',unknown,   formatted,check_file_name_on)
       end if
    end if
!!$    if(mype == 0) then
!!$       if(iprijobstatus >= 1) then
!!$          F_STATUS = 'jobstatus'//name_jobstep
!!$          open(nfstatus, file=F_STATUS, status='unknown',form='formatted')
!!$       end if
!!$    end if
  end subroutine m_Files_open_files_initially

  subroutine m_Files_open_nfdynm_cif_initially
    if(.not.file_existence_contfiles) call m_Files_check_file_existence()
    if(file_existence_contfiles) then
       if(mype==0.and.sw_cif_output==ON) &
           &       call open1(nfdynm_cif,F_DYNM_CIF,'F_DYNM_CIF',unknown,formatted,check_file_name_on)
    else
       if(mype==0.and.sw_cif_output==ON) &
           &       call open0(nfdynm_cif,F_DYNM_CIF,'F_DYNM_CIF',unknown,formatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfdynm_cif_initially

  subroutine m_Files_open_dftd3par()
    logical open
    inquire(unit = nfdftd3par, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfdftd3par,F_DFTD3PAR,'F_DFTD3PAR',unknown,formatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_dftd3par

  subroutine m_Files_close_dftd3par()
    logical open
    inquire(unit = nfdftd3par, opened = open)
    if(open .and. mype == 0) close(nfdftd3par)
  end subroutine m_Files_close_dftd3par

  subroutine m_Files_open_pcont_files_ini()
    if(.not.file_existence_contfiles) call m_Files_check_file_existence()
    if(file_existence_contfiles) then
       if(imdalg == P_CONTROL .or. imdalg == PT_CONTROL .and. mype == 0) then
          call open1(nflatconst,F_LATCONST,'F_LATCONST',unknown,   formatted,check_file_name_on)
          call open1(nfmetric, F_METRIC,   'F_METRIC  ',unknown,   formatted,check_file_name_on)
       endif
    else
       if(imdalg == P_CONTROL .or. imdalg == PT_CONTROL .and. mype == 0)then
          call open0(nflatconst,F_LATCONST,'F_LATCONST',unknown,   formatted,check_file_name_on)
          call open0(nfmetric, F_METRIC,   'F_METRIC  ',unknown,   formatted,check_file_name_on)
          write(nflatconst,*) '#iteration a b c alpha beta gamma vol'
          write(nfmetric,*) '#iteration'
          write(nfmetric,*) '#m11 m12 m13'
          write(nfmetric,*) '#m21 m22 m23'
          write(nfmetric,*) '#m31 m32 m33'
       endif
    endif
  end subroutine m_Files_open_pcont_files_ini

  subroutine m_Files_flush_nfdynm()
    if(mype==0) then
       close(nfdynm,status='keep')
       call open1(nfdynm,F_DYNM,'F_DYNM    ',unknown,   formatted,check_file_name_on)
       if(sw_cif_output==ON)then
          close(nfdynm_cif,status='keep')
          call open1(nfdynm_cif,F_DYNM_CIF,'F_DYNM_CIF',unknown,   formatted,check_file_name_on)
       endif
    end if
  end subroutine m_Files_flush_nfdynm

  subroutine m_Files_flush_pcontrol_files
    if(mype==0)then
       close(nflatconst,status='keep')
       call open1(nflatconst,F_LATCONST,'F_LATCONST',unknown,   formatted,check_file_name_on)
       close(nfmetric,status='keep')
       call open1(nfmetric,F_METRIC,  'F_METRIC  ',unknown,   formatted,check_file_name_on)
    endif
  end subroutine m_Files_flush_pcontrol_files

  subroutine m_Files_flush_nfenf()
    if(mype==0) then
       close(nfenf,status='keep')
       call open1(nfenf, F_ENF,'F_ENF     ',unknown,   formatted,check_file_name_on)
    end if
  end subroutine m_Files_flush_nfenf

  subroutine m_Files_close_files_initial0
    logical :: open
    if(ekmode==ON) then
       if(mype==0) then
          inquire(unit=nfeng, opened=open)
          if(open) close(nfeng,status='keep')
       end if
    else
       if(mype==0) then
          inquire(unit=nfdynm,opened=open)
          if(open) close(nfdynm,status='keep')
          inquire(unit=nfenf,opened=open)
          if(open) close(nfenf,status='keep')
          if(sw_cif_output==ON)then
            inquire(unit=nfdynm_cif,opened=open)
            if(open) close(nfdynm_cif,status='keep')
          endif
          if(imdalg == P_CONTROL .or. imdalg == PT_CONTROL)then
            inquire(unit=nflatconst,opened=open)
            if(open) close(nflatconst,status='keep')
            inquire(unit=nfmetric,opened=open)
            if(open) close(nfmetric,status='keep')
          endif
       end if
    end if
  end subroutine m_Files_close_files_initial0

  subroutine m_Files_open_nfeng(icond)
    integer, intent(in) :: icond
    logical :: open
    if(mype==0) then
       inquire(unit = nfeng, opened = open)
       if(open) close(nfeng)
       if(icond==FIXED_CHARGE_CONTINUATION .and. fixed_charge_k_parallel==ONE_BY_ONE) then
          call open1(nfeng, F_ENERG,'F_ENERG   ',unknown,formatted,check_file_name_on)
       else
          call open0(nfeng, F_ENERG,'F_ENERG   ',unknown,formatted,check_file_name_on)
       end if
    end if
  end subroutine m_Files_open_nfeng

  subroutine m_Files_close_nfeng()
    logical :: open
    if(mype==0) then
       inquire(unit=nfeng, opened = open)
       if(open) close(nfeng)
    end if
  end subroutine m_Files_close_nfeng

! ===
  subroutine m_Files_open_nfwfk_orb_proj(icond)
    integer, intent(in) :: icond
    logical :: open

    if (mype==0) then
       inquire(unit = nfwfk_orb_proj, opened = open)

       if (open) close(nfwfk_orb_proj)

       if (icond == FIXED_CHARGE .or. &
            & (icond==FIXED_CHARGE_CONTINUATION &
            &   .and. fixed_charge_k_parallel==ALL_AT_ONCE)) then
          call open0( nfwfk_orb_proj, F_WFK_ORB_PROJ, 'F_WFK_ORB_PROJ',&
               &      unknown, formatted, check_file_name_on )
       else if (icond==FIXED_CHARGE_CONTINUATION &
            &        .and. fixed_charge_k_parallel==ONE_BY_ONE) then
          call open1( nfwfk_orb_proj, F_WFK_ORB_PROJ, 'F_WFK_ORB_PROJ', &
               &      unknown,formatted, check_file_name_on )
       else
          call open1( nfwfk_orb_proj, F_WFK_ORB_PROJ, 'F_WFK_ORB_PROJ', &
               &      unknown,formatted, check_file_name_on )
       end if
    end if
  end subroutine m_Files_open_nfwfk_orb_proj

  subroutine m_Files_close_nfwfk_orb_proj()
    logical :: open

    if (mype==0) then
       inquire( unit=nfwfk_orb_proj, opened = open )
       if(open) close( nfwfk_orb_proj )
    end if
  end subroutine m_Files_close_nfwfk_orb_proj

  subroutine m_Files_open_nfwfk_lband(icond)
    integer, intent(in) :: icond
    logical :: open

    if (mype==0) then
       inquire(unit = nfwfk_local_decomp, opened = open)

       if (open) close(nfwfk_local_decomp)

       if (icond == FIXED_CHARGE .or. &
            & (icond==FIXED_CHARGE_CONTINUATION &
            &   .and. fixed_charge_k_parallel==ALL_AT_ONCE)) then
          call open0( nfwfk_local_decomp, F_WFK_LOCAL_DECOMP, 'F_WFK_LOCAL_DECOMP',&
               &      unknown, formatted, check_file_name_on )
       else if (icond==FIXED_CHARGE_CONTINUATION &
            &        .and. fixed_charge_k_parallel==ONE_BY_ONE) then
          call open1( nfwfk_local_decomp, F_WFK_LOCAL_DECOMP, 'F_WFK_LOCAL_DECOMP', &
               &      unknown,formatted, check_file_name_on )
       else
          call open1( nfwfk_local_decomp, F_WFK_LOCAL_DECOMP, 'F_WFK_LOCAL_DECOMP', &
               &      unknown,formatted, check_file_name_on )
       end if
    end if
  end subroutine m_Files_open_nfwfk_lband

  subroutine m_Files_close_nfwfk_lband()
    logical :: open

    if (mype==0) then
       inquire( unit=nfwfk_local_decomp, opened = open )
       if(open) close( nfwfk_local_decomp )
    end if
  end subroutine m_Files_close_nfwfk_lband

  subroutine m_Files_open_nfband_spwt(icond)
    integer, intent(in) :: icond
    logical :: open

    if (mype==0) then
       inquire(unit = nfband_spectr_wght, opened = open)

       if (open) close(nfband_spectr_wght)

       if (icond == FIXED_CHARGE .or. &
            & (icond==FIXED_CHARGE_CONTINUATION &
            &   .and. fixed_charge_k_parallel==ALL_AT_ONCE)) then
          call open0( nfband_spectr_wght, F_BAND_SPECTR_WGHT, 'F_BAND_SPECTR_WGHT',&
               &      unknown, formatted, check_file_name_on )
       else if (icond==FIXED_CHARGE_CONTINUATION &
            &        .and. fixed_charge_k_parallel==ONE_BY_ONE) then
          call open1( nfband_spectr_wght, F_BAND_SPECTR_WGHT, 'F_BAND_SPECTR_WGHT ', &
               &      unknown,formatted, check_file_name_on )
       else
          call open1( nfband_spectr_wght, F_BAND_SPECTR_WGHT, 'F_BAND_SPECTR_WGHT', &
               &      unknown,formatted, check_file_name_on )
       end if
    end if
  end subroutine m_Files_open_nfband_spwt

  subroutine m_Files_close_nfband_spwt()
    logical :: open

    if (mype==0) then
       inquire( unit=nfband_spectr_wght, opened = open )
       if(open) close( nfband_spectr_wght )
    end if
  end subroutine m_Files_close_nfband_spwt

  subroutine m_Files_open_nfporb_dens_mat(mode)
    integer, intent(in) :: mode
    logical :: open

    if (mype==0) then
       inquire(unit = nfporb_dens_mat, opened = open)
       if (open) close(nfporb_dens_mat)
       if ( mode == 0 ) then               ! write
          call open0( nfporb_dens_mat, F_PORB_DENS_MAT, 'F_PORB_DENS_MAT',&
               &      unknown, unformatted, check_file_name_on )
       else                                ! read
          call open0( nfporb_dens_mat, F_PORB_DENS_MAT, 'F_PORB_DENS_MAT',&
               &      old, unformatted, check_file_name_on )
       endif
    end if
  end subroutine m_Files_open_nfporb_dens_mat

  subroutine m_Files_close_nfporb_dens_mat()
    logical :: open

    if (mype==0) then
       inquire( unit=nfporb_dens_mat, opened = open )
       if(open) close( nfporb_dens_mat )
    end if
  end subroutine m_Files_close_nfporb_dens_mat

  subroutine m_Files_open_nfporb_rot_mat(mode)
    integer, intent(in) :: mode
    logical :: open

    if (mype==0) then
       inquire(unit = nfporb_rot_mat, opened = open)
       if (open) close(nfporb_rot_mat)
       if ( mode == 0 ) then               ! write
          call open0( nfporb_rot_mat, F_PORB_ROT_MAT, 'F_PORB_ROT_MAT',&
               &      unknown, unformatted, check_file_name_on )
       else                                ! read
          call open0( nfporb_rot_mat, F_PORB_ROT_MAT, 'F_PORB_ROT_MAT',&
               &      old, unformatted, check_file_name_on )
       endif
    end if
  end subroutine m_Files_open_nfporb_rot_mat

  subroutine m_Files_close_nfporb_rot_mat()
    logical :: open

    if (mype==0) then
       inquire( unit=nfporb_rot_mat, opened = open )
       if(open) close( nfporb_rot_mat )
    end if
  end subroutine m_Files_close_nfporb_rot_mat

! ====
  logical function m_Files_nfinp_is_opened()
    logical :: op
    integer :: ierr
    op = .false.
    if(mype == 0) then
       inquire(unit = nfinp, opened = op)
    end if
    call mpi_bcast(op,1,mpi_logical,0,MPI_CommGroup,ierr)
    m_Files_nfinp_is_opened = op
  end function m_Files_nfinp_is_opened

  subroutine m_Files_close_nfinp()
    logical :: op
    if(mype == 0) then
       inquire(unit = nfinp, opened = op)
       if(op) close(nfinp)
    end if
  end subroutine m_Files_close_nfinp

  subroutine m_Files_close_nfcntn()
    logical :: op
    if(mype == 0) then
       inquire(unit = nfcntn, opened = op)
       if(op) close(nfcntn)
    end if
  end subroutine m_Files_close_nfcntn

  subroutine m_Files_open_nfdos(iter,reacid)
    integer,intent(in),optional :: iter
    integer, intent(in), optional :: reacid
    character(len=256) :: striter
    character(len=105) :: retstr,retstr0
    logical :: op
    if(mype==0)then
       inquire(unit = nfdos, opened = op)
       if(op) close(nfdos)
       if(.not.present(iter))then
           call open0(nfdos,F_DOS,'F_DOS     ',unknown,formatted,check_file_name_on)
       else
           write(striter,*) iter
           call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_DOS,retstr)
           if(present(reacid))then
              write(striter,*) reacid
              call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
              retstr = retstr0
           endif
           call open0(nfdos,retstr,'F_DOS     ',unknown,formatted,check_file_name_on)
       endif
    endif
  end subroutine m_Files_open_nfdos

  subroutine m_Files_open_nfegrid()
    if(mype==0) call open0(nfegrid,F_EGRID,'F_EGRID   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfegrid

  subroutine m_Files_open_nfldos()
    if(mype==0) then
       if(icond==FIXED_CHARGE_CONTINUATION .and. ekmode==ON)then
         call open1(nfldos,F_LDOS,'F_LDOS    ',unknown,formatted,check_file_name_on)
       else
         call open0(nfldos,F_LDOS,'F_LDOS    ',unknown,formatted,check_file_name_on)
       endif
    endif
  end subroutine m_Files_open_nfldos

  subroutine m_Files_close_nfldos()
    logical :: open
    if(mype==0) then
       inquire(unit=nfldos,opened=open)
       if(open) close(nfldos,status='keep')
    end if
  end subroutine m_Files_close_nfldos


!!$  subroutine m_Files_close_nfdos()
!!$    logical :: open
!!$    if(mype==0) then
!!$       inquire(unit=nfdos,opened=open)
!!$       if(open) close(nfdos,status='keep')
!!$    end if
!!$  end subroutine m_Files_close_nfdos

  subroutine m_Files_skiptoend(nf)
    integer, intent(in) :: nf
    if(mype == 0) then
1      continue
       read(nf,*, err=2, end=2)
       goto 1
2      continue
       backspace nf
    end if
  end subroutine m_Files_skiptoend

  subroutine m_Files_open_files_initially_p
    call open0(nfinp, F_INP, 'F_INP     ',    old,   formatted,check_file_name_on)
    call open0(nfprm, F_PRM, 'F_PRM     ',unknown,   formatted,check_file_name_on)
  end subroutine m_Files_open_files_initially_p


  subroutine m_Files_close_all_paramset
    close(nfinp,status='keep')
    close(nfprm,status='keep')
  end subroutine m_Files_close_all_paramset

  subroutine open0(nfile,filename,tagname,status, form, write_or_no,fix_endianness)
    integer, intent(in) :: nfile
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: tagname
    integer, intent(in) :: status, form, write_or_no
    logical, intent(in), optional :: fix_endianness
    logical             :: open
    character(len=5) ::    torf
    character(len=14) ::   filestatus
    integer ::             ls,ierror
    logical :: existence
    logical :: endi

    character(len=50) :: wd
    endi = .false.
    if(present(fix_endianness)) endi = fix_endianness

    inquire(unit = nfile, opened = open)

    wd = workdir
!!$    write(6,'(" workdir = ",a)') workdir
    ls = min(max(index(filename,' ')-1,1),stringlength_of_filenames)

    if(.not.open) then
       inquire(file=trim(wd)//filename(1:ls), exist=existence)
       if(status==old .and. .not.existence) then
          ierror = FILE_NOT_EXIST
#ifdef DEBUG_ERRORS
          call phase_error(ierror,nfout,nfile,trim(wd)//filename(1:ls), &
               & __LINE__,__FILE__)
#else
          call phase_error(ierror,nfout,nfile,trim(wd)//filename(1:ls))
#endif
       else
          if(status==old .and. form==formatted) then
             open(nfile,file=trim(wd)//filename(1:ls), status='old', form='formatted')
          else if(status==old .and. form==unformatted) then
#ifdef FIX_ENDIANNESS
          if (endi) then
             open(nfile,file=trim(wd)//filename(1:ls), status='old', form='unformatted', &
          &  convert='little_endian')
          else
             open(nfile,file=trim(wd)//filename(1:ls), status='old', form='unformatted')
          endif
#else
             open(nfile,file=trim(wd)//filename(1:ls), status='old', form='unformatted')
#endif
          else if(status==unknown .and. form==formatted) then
             open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='formatted')
          else if(status==unknown .and. form==unformatted) then
#ifdef FIX_ENDIANNESS
          if (endi) then
             open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='unformatted', &
          &  convert='little_endian')
          else
             open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='unformatted')
          endif
#else
             open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='unformatted')
#endif
          endif
       end if
       torf = 'false'
       filestatus = 'newly opened'
    else
       torf = 'true '
       filestatus = 'already opened'
    end if

!!$    ls = len(trim(wd)//trim(filename))
    if(write_or_no == check_file_name_on .and. printable) then
       if(len(tagname) <= 10) then
          write(nfout,'(1x,a10,a3,a,a3,a)') tagname, " = ", trim(filename(1:ls))," , ",filestatus
       else
          write(nfout,'(1x,a,a3,a,a3,a)') trim(tagname), " = ", trim(filename(1:ls))," , ",filestatus
       end if
    endif
  end subroutine open0

  subroutine open1(nfile,filename,tagname,status, form, write_or_no)
    integer, intent(in) :: nfile
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: tagname
    integer, intent(in) :: status, form, write_or_no
    logical             :: open
    character(len=5) ::    torf
    character(len=14) ::   filestatus
    integer ::             ls

    character(len=50) :: wd

    inquire(unit = nfile, opened = open)

    wd = workdir

!!$    ls = index(filename,' ')-1
!!$    ls = min(max(ls,1),stringlength_of_filenames)
    ls = min(max(index(filename,' ')-1,1),stringlength_of_filenames)

    if(.not.open) then
       if(status==old .and. form==formatted) then
          open(nfile,file=trim(wd)//filename(1:ls), status='old', form='formatted', position='append')
       else if(status==old .and. form==unformatted) then
          open(nfile,file=trim(wd)//filename(1:ls), status='old', form='unformatted', position='append')
       else if(status==unknown .and. form==formatted) then
          open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='formatted', position='append')
       else if(status==unknown .and. form==unformatted) then
          open(nfile,file=trim(wd)//filename(1:ls), status='unknown', form='unformatted', position='append')
       endif
       torf = 'false'
       filestatus = 'newly opened'
    else
       torf = 'true '
       filestatus = 'already opened'
    end if

!!$    ls = len(trim(wd)//trim(filename))
    if(write_or_no == check_file_name_on .and. printable) then
       if(len(tagname)<=10) then
          write(nfout,'(1x,a10,a3,a,a3,a)') tagname," = ",trim(filename(1:ls))," , ",filestatus
       else
          write(nfout,'(1x,a,a3,a,a3,a)') trim(tagname)," = ",trim(filename(1:ls))," , ",filestatus
       end if
    endif
  end subroutine open1

  subroutine m_Files_close_and_clear_nfstop()
    logical :: open
    inquire(unit = nfstop, opened = open)
    if(open) then
       if(sw_wdir==ON) then
          open(nfstop, file=trim(workdir)//F_STOP, status='replace', form='formatted')
       else
          open(nfstop, file=F_STOP, status='replace', form='formatted')
       end if
       close(nfstop, status='keep')
    endif
  end subroutine m_Files_close_and_clear_nfstop

  subroutine m_Files_close_all()
    integer :: i
    logical :: open
    do i = 1, number_of_all_files
       if(n_file(i) == 6) cycle
       inquire(unit = n_file(i), opened = open)
       if(open) then
          close(n_file(i),status='keep')
          if(printable) write(nfout,*) ' closed filenumber = ', n_file(i)
	  if (multiple_replica_mode == OFF) then
            if (n_file(i) == nfstop) then
!!$               open(n_file(i), file=F_STOP, status='replace', form='formatted')
               if(sw_wdir==ON) then
                  open(n_file(i), file=trim(workdir)//F_STOP, status='replace', form='formatted')
               else
                  open(n_file(i), file=F_STOP, status='replace', form='formatted')
               end if
               close(n_file(i), status='keep')
            end if
          end if
       else
	  if (multiple_replica_mode == OFF) then
            if (n_file(i) == nfstop .and. mype == 0) then
               if(sw_wdir==ON) then
                  open(n_file(i), file=trim(workdir)//F_STOP, status='replace', form='formatted')
               else
                  open(n_file(i), file=F_STOP, status='replace', form='formatted')
               end if
               close(n_file(i), status='keep')
            end if
          end if
       endif
    enddo

!!    if(final_c)then
!!$    if(.not.conf_para .or. (conf_para .and. mype_conf==0))then
!!    if(mype == 0) then
!!       close(6, status='keep')
!!    else if(ipriparadeb == 0) then
!!       close(6, status='delete')
!!    else
!!       close(6, status='keep')
!!    end if

!!$    endif
!!    endif

  end subroutine m_Files_close_all

  subroutine m_Files_close_logfile()
    if(mype == 0) then
       close(nfout, status='keep')
#ifndef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
    else if(ipriparadeb == 0) then
       close(nfout, status='delete')
#endif
    else
       close(nfout, status='keep')
    end if
  end subroutine m_Files_close_logfile

  subroutine open_newps_files(ntyp,ierror)
    integer, intent(in)                  :: ntyp
    integer, intent(out)                 :: ierror

    integer :: it, nfp
    logical :: existence

    ierror = 0
    do it = 1, ntyp
       nfp = nfpot(it)
       if(printable) write(nfout,*) ' F_POT(it) = ', F_POT(it)
       if(sw_wdir==ON) then
          inquire(file=workdir//trim(F_POT(it)), exist=existence)
       else
          inquire(file=trim(F_POT(it)), exist=existence)
       end if
       if(.not.existence) then
          ierror = F_POT_FILE_NOT_EXIST
#ifdef DEBUG_ERRORS
          call phase_error(ierror,nfout,nfp,F_POT(it),__LINE__, __FILE__)
#else
          call phase_error(ierror,nfout,nfp,F_POT(it))
#endif
!!$          if(printable) write(nfout,'(" file F_POT(",i3,") does not exist")') it
!!$          if(printable) write(nfout,'(a100)') F_POT(it)
          goto 1001
       end if
       call open0(nfp,F_POT(it),'F_POT     ',OLD,formatted,check_file_name_on)
    enddo
1001 return
  end subroutine open_newps_files

#ifndef _EMPIRICAL_
  subroutine m_Files_open_ps_files(ivan,iatomn,ntyp,ierror)
    integer, intent(in)                  :: ntyp
    integer, intent(in), dimension(ntyp) :: ivan
    real(kind=DP), intent(in), dimension(ntyp) :: iatomn
    integer, intent(out)                 :: ierror

    integer     :: nfpp, nfp, it, n_non_vanderbilt, n_pcc, itpcc, check_file_name
    real(kind=DP) :: ival
    real(kind=DP), dimension(2)   :: alp, cc
    logical :: existence

!!$    if(ppprinted) then
!!$       check_file_name = check_file_name_off
!!$    else
!!$       check_file_name = check_file_name_on
!!$    end if
    checK_file_name = check_file_name_off

    ierror = 0
    if(mype /= 0) return
    n_non_Vanderbilt = 0
    n_pcc            = 0
    nfpp             = 0
    do it = 1, ntyp
       if(ivan(it) /= OLD) then
          nfpp = nfpp + 1
          nfp = nfpot(nfpp)
          if(printable .and. .not.ppprinted ) write(nfout,*) ' F_POT(nfpp) = ', trim(F_POT(nfpp))
          if(sw_wdir==ON)then
             inquire(file=trim(workdir)//F_POT(nfpp), exist=existence)
          else
             inquire(file=F_POT(nfpp), exist=existence)
          endif
          if(.not.existence) then
             ierror = F_POT_FILE_NOT_EXIST
#ifdef DEBUG_ERRORS
             call phase_error(ierror,nfout,nfp,F_POT(nfpp),__LINE__, __FILE__)
#else
             call phase_error(ierror,nfout,nfp,(F_POT(nfpp)))
#endif
!!$             if(printable) write(nfout,'(" file F_POT(",i3,") does not exist")') nfpp
!!$             if(printable) write(nfout,'(a100)') F_POT(nfpp)
             goto 1001
          end if
          call open0(nfp,F_POT(nfpp),'F_POT     ',OLD,formatted,check_file_name)
       else if(ivan(it) == OLD) then
          n_non_Vanderbilt = n_non_Vanderbilt + 1
          call psbhs0(nfout,nint(iatomn(it)),ival,itpcc,alp,cc)  ! -(b_PseudoPotential)
          if(itpcc == ON) n_pcc = n_pcc + 1
       endif
    enddo
    if(n_non_Vanderbilt >= 1 ) then
       call open0(nfpkb,F_PKB,'F_PKB     ',OLD,formatted,check_file_name)
       call open0(nfpd, F_PD, 'F_PD      ',OLD,formatted,check_file_name)
       rewind nfpkb
       rewind nfpd
       if(n_pcc >= 1) then
          call open0(nfppc,F_PPC,'F_PPC     ',OLD,formatted,check_file_name)
          rewind nfppc
       endif
    endif
1001 return
  end subroutine m_Files_open_ps_files

  subroutine m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
    integer, intent(in)                  :: ntyp,it
    integer, intent(in), dimension(ntyp) :: ivan
    real(kind=DP), intent(in), dimension(ntyp) :: iatomn
    integer, intent(out)                 :: ierror

    integer     ::  nfp,  n_non_vanderbilt, n_pcc, itpcc, check_file_name
    integer, save :: nfpp=0
    real(kind=DP) :: ival
    real(kind=DP), dimension(2)   :: alp, cc
    logical :: existence

!!$    if(ppprinted) then
!!$       check_file_name = check_file_name_off
!!$    else
!!$       check_file_name = check_file_name_on
!!$    end if
    checK_file_name = check_file_name_off

    ierror = 0
    if(mype /= 0) return
    n_non_Vanderbilt = 0
    n_pcc            = 0
    if(ivan(it) /= OLD) then
       nfpp = nfpp + 1
       nfp = nfpot(nfpp)
       if(printable .and. .not.ppprinted ) write(nfout,*) ' F_POT(nfpp) = ', trim(F_POT(nfpp))
       if(sw_wdir==ON)then
          inquire(file=trim(workdir)//F_POT(nfpp), exist=existence)
       else
          inquire(file=F_POT(nfpp), exist=existence)
       endif
       if(.not.existence) then
          ierror = F_POT_FILE_NOT_EXIST
#ifdef DEBUG_ERRORS
          call phase_error(ierror,nfout,nfp,F_POT(nfpp),__LINE__, __FILE__)
#else
          call phase_error(ierror,nfout,nfp,(F_POT(nfpp)))
#endif
!!$             if(printable) write(nfout,'(" file F_POT(",i3,") does not exist")') nfpp
!!$             if(printable) write(nfout,'(a100)') F_POT(nfpp)
          goto 1001
       end if
       call open0(nfp,F_POT(nfpp),'F_POT     ',OLD,formatted,check_file_name)
    else if(ivan(it) == OLD) then
       n_non_Vanderbilt = n_non_Vanderbilt + 1
       call psbhs0(nfout,nint(iatomn(it)),ival,itpcc,alp,cc)  ! -(b_PseudoPotential)
       if(itpcc == ON) n_pcc = n_pcc + 1
    endif

    if(n_non_Vanderbilt >= 1 ) then
       call open0(nfpkb,F_PKB,'F_PKB     ',OLD,formatted,check_file_name)
       call open0(nfpd, F_PD, 'F_PD      ',OLD,formatted,check_file_name)
       rewind nfpkb
       rewind nfpd
       if(n_pcc >= 1) then
          call open0(nfppc,F_PPC,'F_PPC     ',OLD,formatted,check_file_name)
          rewind nfppc
       endif
    endif

    if(it==ntyp) nfpp = 0
1001 return
  end subroutine m_Files_open_ps_file

  subroutine m_Files_close_ps_files
    integer :: i
    logical :: open
    if(mype /= 0) return
    do i = 1, MAXNSP
       inquire(unit = nfpot(i), opened = open)
       if(open) close(nfpot(i),status='keep')
    enddo
    inquire(unit = nfpkb, opened = open)
    if(open) close(nfpkb,status='keep')
    inquire(unit = nfpd,  opened = open)
    if(open) close(nfpd, status='keep')
    inquire(unit = nfppc, opened = open)
    if(open) close(nfppc,status='keep')
  end subroutine m_Files_close_ps_files

  subroutine m_Files_close_ps_file(it)
    integer, intent(in) :: it
    logical :: open
    inquire(unit = nfpot(it), opened = open)
    if(open) close(nfpot(it),status='keep')
    inquire(unit = nfpkb, opened = open)
    if(open) close(nfpkb,status='keep')
    inquire(unit = nfpd,  opened = open)
    if(open) close(nfpd, status='keep')
    inquire(unit = nfppc, opened = open)
    if(open) close(nfppc,status='keep')
  end subroutine m_Files_close_ps_file
#endif


  subroutine insert_str_before_dot(pref,origstr,retstr)
    character(len=*),intent(in) :: pref
    character(len=*),intent(in) :: origstr
    character(len=105),intent(out) :: retstr
    integer :: lenstr, ip, i
    lenstr = len_trim(origstr)
    ip = lenstr+1
    findingdot: do i = lenstr, 1, -1
      if(origstr(i:i) == '/') then
         exit findingdot
      else if(origstr(i:i) == '.') then
         ip = i-1
         exit findingdot
      end if
    end do findingdot

    if(ip >= 1 .and. ip <= lenstr) then
       if(ip+1 <= lenstr) then
          retstr = origstr(1:ip)//trim(adjustl(pref))//origstr(ip+1:lenstr)
       else
          retstr = origstr(1:ip)//trim(adjustl(pref))
       end if
    else
       retstr = origstr(1:lenstr)//trim(adjustl(pref))
    end if
  end subroutine insert_str_before_dot

  subroutine m_Files_open_nfchr( nspin, ispin, iter, reacid, add_core )
    integer, intent(in) :: nspin, ispin
    integer, intent(in), optional :: iter
    integer, intent(in), optional :: reacid
    logical, intent(in), optional :: add_core

    integer :: lenstr, ip, i
    character(len=105) :: F_CHRUD
    character(len=256) :: striter
    character(len=105) :: retstr,retstr0
    logical :: open
    logical :: flag_add_core

    if ( present(add_core) ) then
       flag_add_core = add_core
    else
       flag_add_core = .false.
    endif

    if(mype == 0) then
       if(nspin == 1) then
          inquire(unit=nfchr,opened=open)
          if(open) close(nfchr,status='keep')

          if(.not.present(iter)) then
             if ( .not. flag_add_core ) then
                call open0(nfchr,F_CHR,'F_CHR     ',unknown,formatted,check_file_name_on)
             else
                call insert_str_before_dot( '_ae', F_CHR, retstr )
                call open0(nfchr,retstr,'F_CHR     ',unknown,formatted,check_file_name_on)
             endif
          else
             write(striter,*) iter
             call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHR,retstr)
             if(present(reacid))then
                write(striter,*) reacid
                call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                retstr = retstr0
             endif
             call open0(nfchr,retstr,'F_CHR     ',unknown,formatted,check_file_name_on)
          endif

       else if(nspin == 2) then
          lenstr = len_trim(F_CHR)
          ip = lenstr+1
          findingdot: do i = lenstr, 1, -1
             if(F_CHR(i:i) == '/') then
                exit findingdot
             else if(F_CHR(i:i) == '.') then
                ip = i
                exit findingdot
             end if
          end do findingdot

          if(ispin == 1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRUD = F_CHR(1:ip)//'up.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRUD = F_CHR(1:ip)//'up'
                end if
             else
                F_CHRUD = F_CHR(1:lenstr)//'.up'
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRUD (up) ---'
                write(nfout,*) F_CHRUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')

             if(.not.present(iter)) then
                if ( .not. flag_add_core ) then
                   call open0( nfchr, F_CHRUD, 'F_CHR.UP  ', unknown, &
                        &      formatted, check_file_name_on )
                else
                   call insert_str_before_dot( '_ae', F_CHRUD, retstr )
                   call open0( nfchr,retstr, 'F_CHR.UP  ', unknown, &
                        &      formatted,check_file_name_on )
                endif

             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRUD,retstr)
                call open0(nfchr,retstr,'F_CHR.UP  ',unknown,formatted,check_file_name_on)
             endif
          else if(ispin == 2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRUD = F_CHR(1:ip)//'down.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRUD = F_CHR(1:ip)//'down'
                end if
             else
                F_CHRUD = F_CHR(1:lenstr)//'.down'
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRUD (down) ---'
                write(nfout,*) F_CHRUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')

             if(.not.present(iter)) then
               if ( .not. flag_add_core ) then
                   call open0( nfchr, F_CHRUD, 'F_CHR.DOWN', unknown, &
                        &      formatted, check_file_name_on )
                else
                   call insert_str_before_dot( '_ae', F_CHRUD, retstr )
                   call open0( nfchr, retstr, 'F_CHR.DOWN', unknown, &
                        &      formatted, check_file_name_on )
                endif

             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRUD,retstr)
                call open0(nfchr,retstr,'F_CHR.DOWN',unknown,formatted,check_file_name_on)
             endif
          end if

! ======= KT_add ========== 2014/06/07
          if (ispin == -1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRUD = F_CHR(1:ip)//'tot.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRUD = F_CHR(1:ip)//'tot'
                end if
             else
                F_CHRUD = F_CHR(1:lenstr)//'.tot'
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRUD (tot) ---'
                write(nfout,*) F_CHRUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')

             if(.not.present(iter)) then
                if ( .not. flag_add_core ) then
                   call open0( nfchr, F_CHRUD, 'F_CHR.TOT ', &
                        &      unknown, formatted, check_file_name_on )
                else
                   call insert_str_before_dot( '_ae', F_CHRUD, retstr )
                   call open0( nfchr, retstr, 'F_CHR.TOT ',unknown, &
                        &      formatted, check_file_name_on )
                endif
             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRUD,retstr)
                call open0(nfchr,retstr,'F_CHR.TOT ',unknown,formatted,check_file_name_on)
             endif
          else if(ispin == -2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRUD = F_CHR(1:ip)//'mag.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRUD = F_CHR(1:ip)//'mag'
                end if
             else
                F_CHRUD = F_CHR(1:lenstr)//'.mag'
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRUD (mag) ---'
                write(nfout,*) F_CHRUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')

             if(.not.present(iter)) then
                call open0(nfchr,F_CHRUD,'F_CHR.MAG ',unknown,formatted,check_file_name_on)
             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRUD,retstr)
                call open0(nfchr,retstr,'F_CHR.MAG ',unknown,formatted,check_file_name_on)
             endif
          endif
! ========================= 2014/06/07

       end if
    end if
  end subroutine m_Files_open_nfchr

!====================================== added by K. Tagami =========== 11.0
  subroutine m_Files_open_nfchr_noncl( iloop, iter, add_core )
    integer, intent(in) :: iloop
    integer, intent(in), optional :: iter
    logical, intent(in), optional :: add_core

    integer :: lenstr, ip, i
    character(len=105) :: F_CHRMAG
    character(len=105) :: retstr,retstr0
    logical :: open
    logical :: flag_add_core

    if ( present(add_core) ) then
       flag_add_core = add_core
    else
       flag_add_core = .false.
    endif

    if(mype == 0) then

       lenstr = len_trim(F_CHR)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
         if(F_CHR(i:i) == '/') then
            exit findingdot
         else if(F_CHR(i:i) == '.') then
            ip = i
            exit findingdot
         end if
       end do findingdot

       select case (iloop)
       case(1)
           if(ip >= 1 .and. ip <= lenstr) then
              if(ip+1 <= lenstr) then
                 F_CHRMAG = F_CHR(1:ip)//'tot.'//F_CHR(ip+1:lenstr)
              else
                 F_CHRMAG = F_CHR(1:ip)//'tot'
              end if
           else
              F_CHRMAG = F_CHR(1:lenstr)//'.tot'
           end if
           if(printable) then
              write(nfout,*) ' --- F_CHRMAG (tot) ---'
              write(nfout,*) F_CHRMAG
           end if
           inquire(unit=nfchr,opened=open)
           if(open) close(nfchr,status='keep')

           if ( .not. flag_add_core ) then
              call open0( nfchr, F_CHRMAG, 'F_CHR.TOT ',unknown, &
                   &           formatted, check_file_name_on )
           else
              call insert_str_before_dot( '_ae', F_CHRMAG, retstr )
              call open0( nfchr, retstr, 'F_CHR.TOT ',unknown, &
                   &      formatted, check_file_name_on )
           endif

       case(2)
           if(ip >= 1 .and. ip <= lenstr) then
              if(ip+1 <= lenstr) then
                 F_CHRMAG = F_CHR(1:ip)//'mx.'//F_CHR(ip+1:lenstr)
              else
                 F_CHRMAG = F_CHR(1:ip)//'mx'
              end if
           else
              F_CHRMAG = F_CHR(1:lenstr)//'.mx'
           end if
           if(printable) then
              write(nfout,*) ' --- F_CHRMAG (mx) ---'
              write(nfout,*) F_CHRMAG
           end if
           inquire(unit=nfchr,opened=open)
           if(open) close(nfchr,status='keep')
           call open0( nfchr, F_CHRMAG, 'F_CHR.MX  ',unknown, &
                &           formatted, check_file_name_on )

       case(3)
           if(ip >= 1 .and. ip <= lenstr) then
              if(ip+1 <= lenstr) then
                 F_CHRMAG = F_CHR(1:ip)//'my.'//F_CHR(ip+1:lenstr)
              else
                 F_CHRMAG = F_CHR(1:ip)//'my'
              end if
           else
              F_CHRMAG = F_CHR(1:lenstr)//'.my'
           end if
           if(printable) then
              write(nfout,*) ' --- F_CHRMAG (my) ---'
              write(nfout,*) F_CHRMAG
           end if
           inquire(unit=nfchr,opened=open)
           if(open) close(nfchr,status='keep')
           call open0( nfchr, F_CHRMAG, 'F_CHR.MY  ',unknown, &
                &           formatted, check_file_name_on )

       case(4)
           if(ip >= 1 .and. ip <= lenstr) then
              if(ip+1 <= lenstr) then
                 F_CHRMAG = F_CHR(1:ip)//'mz.'//F_CHR(ip+1:lenstr)
              else
                 F_CHRMAG = F_CHR(1:ip)//'mz'
              end if
           else
              F_CHRMAG = F_CHR(1:lenstr)//'.mz'
           end if
           if(printable) then
              write(nfout,*) ' --- F_CHRMAG (mz) ---'
              write(nfout,*) F_CHRMAG
           end if
           inquire(unit=nfchr,opened=open)
           if(open) close(nfchr,status='keep')
           call open0( nfchr, F_CHRMAG, 'F_CHR.MZ  ',unknown, &
                &           formatted, check_file_name_on )
       end select
    end if
  end subroutine m_Files_open_nfchr_noncl
! ================================================================== 11.0

  subroutine m_Files_open_nfchr_pc(nspin,ispin,iw,iter,reacid)
    integer, intent(in) :: nspin, ispin, iw
    integer, intent(in),optional :: iter
    integer, intent(in),optional :: reacid
    integer :: lenstr, ip, i
    character(len=105) :: F_CHRpcUD, F_CHRpc
    character(len=256) :: striter
    character(len=105) :: retstr,retstr0
    integer, parameter :: len_nnumber = 4
    character(len=len_nnumber) :: nnumber
    logical :: open

    if(mype == 0) then
       write(nnumber,'(i4)') iw
       do i = 1, len_nnumber
          if(nnumber(i:i) == ' ') nnumber(i:i) = '0'
       end do

       lenstr = len_trim(F_CHR)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
          if(F_CHR(i:i) == '/') then
             exit findingdot
          else if(F_CHR(i:i) == '.') then
             ip = i
             exit findingdot
          end if
       end do findingdot

       if(nspin == 1) then
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_CHRpc = F_CHR(1:ip)//nnumber//'.'//F_CHR(ip+1:lenstr)
             else
                F_CHRpc = F_CHR(1:ip)//nnumber
             end if
          else
             F_CHRpc = F_CHR(1:lenstr)//nnumber
          end if
          inquire(unit=nfchr,opened=open)
          if(open) close(nfchr,status='keep')
          if(.not.present(iter)) then
             call open0(nfchr,F_CHRpc,'F_CHRpc   ',unknown,formatted,check_file_name_on)
          else
             write(striter,*) iter
             call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRpc,retstr)
             if(present(reacid))then
                write(striter,*) reacid
                call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                retstr = retstr0
             endif
             call open0(nfchr,retstr,'F_CHRpc   ',unknown,formatted,check_file_name_on)
          endif

       else if(nspin == 2) then
          if(ispin == 1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcUD = F_CHR(1:ip)//'up.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcUD = F_CHR(1:ip)//'up.'//nnumber
                end if
             else
                F_CHRpcUD = F_CHR(1:lenstr)//'.up.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcUD (up) ---'
                write(nfout,*) F_CHRpcUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             if(.not.present(iter)) then
                call open0(nfchr,F_CHRpcUD,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRpcUD,retstr)
                if(present(reacid))then
                  write(striter,*) reacid
                  call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                  retstr = retstr0
                endif
                call open0(nfchr,retstr,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             endif
          else if(ispin == 2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcUD = F_CHR(1:ip)//'down.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcUD = F_CHR(1:ip)//'down.'//nnumber
                end if
             else
                F_CHRpcUD = F_CHR(1:lenstr)//'.down.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcUD (down) ---'
                write(nfout,*) F_CHRpcUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             if(.not.present(iter)) then
                 call open0(nfchr,F_CHRpcUD,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             else
                 write(striter,*) iter
                 call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRpcUD,retstr)
                 if(present(reacid))then
                   write(striter,*) reacid
                   call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                   retstr = retstr0
                 endif
                 call open0(nfchr,retstr,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             endif
          end if
! =========== KT_add ======== 2014/06/07
          if(ispin == -1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcUD = F_CHR(1:ip)//'tot.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcUD = F_CHR(1:ip)//'tot.'//nnumber
                end if
             else
                F_CHRpcUD = F_CHR(1:lenstr)//'.tot.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcUD (tot) ---'
                write(nfout,*) F_CHRpcUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             if(.not.present(iter)) then
                call open0(nfchr,F_CHRpcUD,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             else
                write(striter,*) iter
                call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRpcUD,retstr)
                if(present(reacid))then
                  write(striter,*) reacid
                  call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                  retstr = retstr0
                endif
                call open0(nfchr,retstr,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             endif
          else if(ispin == -2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcUD = F_CHR(1:ip)//'mag.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcUD = F_CHR(1:ip)//'mag.'//nnumber
                end if
             else
                F_CHRpcUD = F_CHR(1:lenstr)//'.mag.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcUD (mag) ---'
                write(nfout,*) F_CHRpcUD
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             if(.not.present(iter)) then
                 call open0(nfchr,F_CHRpcUD,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             else
                 write(striter,*) iter
                 call insert_str_before_dot('_iter'//trim(adjustl(striter)),F_CHRpcUD,retstr)
                 if(present(reacid))then
                   write(striter,*) reacid
                   call insert_str_before_dot('_reac'//trim(adjustl(striter)),retstr,retstr0)
                   retstr = retstr0
                 endif
                 call open0(nfchr,retstr,'F_CHRpcUD ',unknown,formatted,check_file_name_on)
             endif
          end if
! =========================== 2014/06/07
       end if
    end if
  end subroutine m_Files_open_nfchr_pc

! ================================ added by K. Tagami ================== 11.0
  subroutine m_Files_open_nfchr_pc_noncl( iloop, iw, iter )
    integer, intent(in) :: iloop, iw
    integer, intent(in), optional :: iter

    integer :: lenstr, ip, i
    character(len=105) :: F_CHRpcMAG, F_CHRpc
    integer, parameter :: len_nnumber = 4
    character(len=len_nnumber) :: nnumber
    logical :: open

    if(mype == 0) then
       write(nnumber,'(i4)') iw
       do i = 1, len_nnumber
          if(nnumber(i:i) == ' ') nnumber(i:i) = '0'
       end do

       lenstr = len_trim(F_CHR)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
          if(F_CHR(i:i) == '/') then
             exit findingdot
          else if(F_CHR(i:i) == '.') then
             ip = i
             exit findingdot
          end if
       end do findingdot

       select case (iloop)
       case(1)
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcMAG = F_CHR(1:ip)//'tot.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcMAG = F_CHR(1:ip)//'tot.'//nnumber
                end if
             else
                F_CHRpcMAG = F_CHR(1:lenstr)//'.tot.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcMAG (tot) ---'
                write(nfout,*) F_CHRpcMAG
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             call open0( nfchr, F_CHRpcMAG, 'F_CHRpcMAG ', unknown, &
	&                formatted, check_file_name_on )

       case(2)
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcMAG = F_CHR(1:ip)//'mx.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcMAG = F_CHR(1:ip)//'mx.'//nnumber
                end if
             else
                F_CHRpcMAG = F_CHR(1:lenstr)//'.mx.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcMAG (mx) ---'
                write(nfout,*) F_CHRpcMAG
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             call open0( nfchr, F_CHRpcMAG, 'F_CHRpcMAGx ', unknown, &
	&                formatted, check_file_name_on )

       case(3)
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcMAG = F_CHR(1:ip)//'my.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcMAG = F_CHR(1:ip)//'my.'//nnumber
                end if
             else
                F_CHRpcMAG = F_CHR(1:lenstr)//'.my.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcMAG (my) ---'
                write(nfout,*) F_CHRpcMAG
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             call open0( nfchr, F_CHRpcMAG, 'F_CHRpcMAGy ', unknown, &
	&                formatted, check_file_name_on )

       case(4)
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_CHRpcMAG = F_CHR(1:ip)//'mz.'//nnumber//'.'//F_CHR(ip+1:lenstr)
                else
                   F_CHRpcMAG = F_CHR(1:ip)//'mz.'//nnumber
                end if
             else
                F_CHRpcMAG = F_CHR(1:lenstr)//'.mz.'//nnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_CHRpcMAG (tot) ---'
                write(nfout,*) F_CHRpcMAG
             end if
             inquire(unit=nfchr,opened=open)
             if(open) close(nfchr,status='keep')
             call open0( nfchr, F_CHRpcMAG, 'F_CHRpcMAGz ', unknown, &
	&                formatted, check_file_name_on )
       end select

    end if
  end subroutine m_Files_open_nfchr_pc_noncl
! ===================================================================== 11.0

  subroutine m_Files_open_nfwfk(nspin,ik,ib)
    integer, intent(in) :: nspin, ik, ib
    integer :: lenstr, ip, i, ispin, nk
    character(len=105) :: F_WFnUD, F_WFn
    integer, parameter :: len_nnumber = 4
    integer, parameter :: len_mnumber = 5
    character(len=len_nnumber) :: nnumber
    character(len=len_mnumber) :: mnumber
    logical :: open

    ispin = mod(ik-1,nspin)+1
    nk = (ik-1)/nspin+1

    if(mype == 0) then
       write(nnumber,'(i4)') nk
       do i = 1, len_nnumber
          if(nnumber(i:i) == ' ') nnumber(i:i) = '0'
       end do
       write(mnumber,'(i5)') ib
       do i = 1, len_mnumber
          if(mnumber(i:i) == ' ') mnumber(i:i) = '0'
       end do

       lenstr = len_trim(F_WFk)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
          if(F_WFk(i:i) == '/') then
             exit findingdot
          else if(F_WFk(i:i) == '.') then
             ip = i
             exit findingdot
          end if
       end do findingdot

       if(nspin == 1) then
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn = F_WFk(1:ip)//"k"//nnumber//"n"//mnumber//'.'//F_WFk(ip+1:lenstr)
             else
                F_WFn = F_WFk(1:ip)//"k"//nnumber//"n"//mnumber
             end if
          else
             F_WFn = F_WFk(1:lenstr)//".k"//nnumber//"n"//mnumber
          end if
          inquire(unit=nfwfk,opened=open)
          if(open) close(nfwfk,status='keep')
          call open0(nfwfk,F_WFn,'F_WFn     ',unknown,formatted,check_file_name_on)

       else if(nspin == 2) then
          if(ispin == 1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_WFnUD = F_WFk(1:ip)//'up.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk(ip+1:lenstr)
                else
                   F_WFnUD = F_WFk(1:ip)//'up.'//"k"//nnumber//"n"//mnumber
                end if
             else
                F_WFnUD = F_WFk(1:lenstr)//'.up.'//"k"//nnumber//"n"//mnumber
             end if
              if(printable) then
                write(nfout,*) ' --- F_WFnUD (up) ---'
                write(nfout,*) F_WFnUD
             end if
             inquire(unit=nfwfk,opened=open)
             if(open) close(nfwfk,status='keep')
             call open0(nfwfk,F_WFnUD,'F_WFnUD   ',unknown,formatted,check_file_name_on)
          else if(ispin == 2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_WFnUD = F_WFk(1:ip)//'down.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk(ip+1:lenstr)
                else
                   F_WFnUD = F_WFk(1:ip)//'down.'//"k"//nnumber//"n"//mnumber
                end if
             else
                F_WFnUD = F_WFk(1:lenstr)//'.down.'//"k"//nnumber//"n"//mnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_WFnUD (down) ---'
                write(nfout,*) F_WFnUD
             end if
             inquire(unit=nfwfk,opened=open)
             if(open) close(nfwfk,status='keep')
             call open0(nfwfk,F_WFnUD,'F_WFnUD   ',unknown,formatted,check_file_name_on)
          end if
       end if
    end if
  end subroutine m_Files_open_nfwfk

! ====
  subroutine m_Files_open_nfwfksq_noncl(iloop,ik,ib)
    integer, intent(in) :: ik, ib, iloop
    integer :: lenstr, ip, i
    character(len=105) :: F_WFn_MAG
    integer, parameter :: len_nnumber = 4
    integer, parameter :: len_mnumber = 5
    character(len=len_nnumber) :: nnumber
    character(len=len_mnumber) :: mnumber
    logical :: open

    if(mype == 0) then
       write(nnumber,'(i4)') ik
       do i = 1, len_nnumber
          if(nnumber(i:i) == ' ') nnumber(i:i) = '0'
       end do
       write(mnumber,'(i5)') ib
       do i = 1, len_mnumber
          if(mnumber(i:i) == ' ') mnumber(i:i) = '0'
       end do

       lenstr = len_trim(F_WFk_squared)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
          if(F_WFk_squared(i:i) == '/') then
             exit findingdot
          else if(F_WFk_squared(i:i) == '.') then
             ip = i
             exit findingdot
          end if
       end do findingdot

       select case (iloop)
       case (1)
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn_Mag = F_WFk_squared(1:ip)//'tot.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk_squared(ip+1:lenstr)
             else
                F_WFn_Mag = F_WFk_squared(1:ip)//'tot.'//"k"//nnumber//"n"//mnumber
             end if
          else
             F_WFn_Mag = F_WFk_squared(1:lenstr)//'.tot.'//"k"//nnumber//"n"//mnumber
          end if
          if(printable) then
             write(nfout,*) ' --- F_WFn_Mag (tot) ---'
             write(nfout,*) F_WFn_Mag
          end if

          inquire(unit=nfwfk_sq,opened=open)
          if(open) close(nfwfk_sq,status='keep')

          call open0( nfwfk_sq, F_WFn_Mag, 'F_WFn_Mag  ', unknown, formatted, &
               &      check_file_name_on )

       case (2)
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn_Mag = F_WFk_squared(1:ip)//'mx.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk_squared(ip+1:lenstr)
             else
                F_WFn_Mag = F_WFk_squared(1:ip)//'mx.'//"k"//nnumber//"n"//mnumber
             end if
          else
             F_WFn_Mag = F_WFk_squared(1:lenstr)//'.mx.'//"k"//nnumber//"n"//mnumber
          end if
          if(printable) then
             write(nfout,*) ' --- F_WFn_Mag (mx) ---'
             write(nfout,*) F_WFn_Mag
          end if

          inquire(unit=nfwfk_sq,opened=open)
          if(open) close(nfwfk_sq,status='keep')

          call open0(nfwfk_sq,F_WFn_Mag,'F_WFn_Mag  ',unknown,formatted,check_file_name_on)

       case (3)
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn_Mag = F_WFk_squared(1:ip)//'my.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk_squared(ip+1:lenstr)
             else
                F_WFn_Mag = F_WFk_squared(1:ip)//'my.'//"k"//nnumber//"n"//mnumber
             end if
          else
             F_WFn_Mag = F_WFk_squared(1:lenstr)//'.my.'//"k"//nnumber//"n"//mnumber
          end if
          if(printable) then
             write(nfout,*) ' --- F_WFn_Mag (my) ---'
             write(nfout,*) F_WFn_Mag
          end if

          inquire(unit=nfwfk_sq,opened=open)
          if(open) close(nfwfk_sq,status='keep')

          call open0(nfwfk_sq,F_WFn_Mag,'F_WFn_Mag  ',unknown,formatted,check_file_name_on)

       case (4)
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn_Mag = F_WFk_squared(1:ip)//'mz.'//"k"//nnumber//"n"//mnumber//'.'//F_WFk_squared(ip+1:lenstr)
             else
                F_WFn_Mag = F_WFk_squared(1:ip)//'mz.'//"k"//nnumber//"n"//mnumber
             end if
          else
             F_WFn_Mag = F_WFk_squared(1:lenstr)//'.mz.'//"k"//nnumber//"n"//mnumber
          end if
          if(printable) then
             write(nfout,*) ' --- F_WFn_Mag (mz) ---'
             write(nfout,*) F_WFn_Mag
          end if

          inquire(unit=nfwfk_sq,opened=open)
          if(open) close(nfwfk_sq,status='keep')

          call open0(nfwfk_sq,F_WFn_Mag,'F_WFn_Mag  ',unknown,formatted,check_file_name_on)
       end select
    end if
  end subroutine m_Files_open_nfwfksq_noncl

  subroutine m_Files_open_nfwfk_integ_mom(icond)
    integer, intent(in) :: icond
    logical :: open
    if(mype==0) then
       inquire(unit = nfwfk_integ_mom, opened = open)
       if(open) close(nfwfk_integ_mom)
       if(icond == FIXED_CHARGE .or. &
            & ( icond==FIXED_CHARGE_CONTINUATION &
            &       .and. fixed_charge_k_parallel==ALL_AT_ONCE ) ) then
          call open0( nfwfk_integ_mom, F_WFk_IntegMoment, 'F_WFk_IntegMoment', &
               &     unknown, formatted, check_file_name_on )
       else if(icond==FIXED_CHARGE_CONTINUATION &
            &       .and. fixed_charge_k_parallel==ONE_BY_ONE) then
          call open1( nfwfk_integ_mom, F_WFk_IntegMoment, 'F_WFk_IntegMoment', &
               &      unknown,formatted, check_file_name_on )

       else if (icond == INITIAL) then
          call open0( nfwfk_integ_mom, F_WFk_IntegMoment, 'F_WFk_IntegMoment', &
               &     unknown,formatted, check_file_name_on )
       else
          call open1( nfwfk_integ_mom,  F_WFk_IntegMoment,'F_WFk_IntegMoment', &
               &     unknown,formatted, check_file_name_on )
       end if
    end if
  end subroutine m_Files_open_nfwfk_integ_mom

  subroutine m_Files_close_nfwfk_integ_mom()
    logical :: open
    if(mype==0) then
       inquire(unit= nfwfk_integ_mom, opened = open)
       if(open) close( nfwfk_integ_mom )
    end if
  end subroutine m_Files_close_nfwfk_integ_mom
! ====

  subroutine m_Files_open_nfwannier(nspin,ispin,ib)
    integer, intent(in) :: nspin, ispin, ib
    integer :: lenstr, ip, i
    character(len=105) :: F_WFnUD, F_WFn
    integer, parameter :: len_mnumber = 5
    character(len=len_mnumber) :: mnumber
    logical :: open

    if(mype == 0) then
       write(mnumber,'(i5)') ib
       do i = 1, len_mnumber
          if(mnumber(i:i) == ' ') mnumber(i:i) = '0'
       end do

       lenstr = len_trim(F_WANNIER)
       ip = lenstr+1
       findingdot: do i = lenstr, 1, -1
          if(F_WFk(i:i) == '/') then
             exit findingdot
          else if(F_WANNIER(i:i) == '.') then
             ip = i
             exit findingdot
          end if
       end do findingdot

       if(nspin == 1) then
          if(ip >= 1 .and. ip <= lenstr) then
             if(ip+1 <= lenstr) then
                F_WFn = F_WANNIER(1:ip)//mnumber//'.'//F_WANNIER(ip+1:lenstr)
             else
                F_WFn = F_WANNIER(1:ip)//mnumber
             end if
          else
             F_WFn = F_WANNIER(1:lenstr)//"."//mnumber
          end if
          inquire(unit=nfwannier,opened=open)
          if(open) close(nfwannier,status='keep')
          call open0(nfwannier,F_WFn,'F_WFn     ',unknown,formatted,check_file_name_on)

       else if(nspin == 2) then
          if(ispin == 1) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_WFnUD = F_WANNIER(1:ip)//'up.'//mnumber//'.'//F_WANNIER(ip+1:lenstr)
                else
                   F_WFnUD = F_WANNIER(1:ip)//'up.'//mnumber
                end if
             else
                F_WFnUD = F_WANNIER(1:lenstr)//'.up.'//mnumber
             end if
              if(printable) then
                write(nfout,*) ' --- F_WFnUD (up) ---'
                write(nfout,*) F_WFnUD
             end if
             inquire(unit=nfwannier,opened=open)
             if(open) close(nfwannier,status='keep')
             call open0(nfwannier,F_WFnUD,'F_WFnUD   ',unknown,formatted,check_file_name_on)
          else if(ispin == 2) then
             if(ip >= 1 .and. ip <= lenstr) then
                if(ip+1 <= lenstr) then
                   F_WFnUD = F_WANNIER(1:ip)//'down.'//mnumber//'.'//F_WANNIER(ip+1:lenstr)
                else
                   F_WFnUD = F_WANNIER(1:ip)//'down.'//mnumber
                end if
             else
                F_WFnUD = F_WANNIER(1:lenstr)//'.down.'//mnumber
             end if
             if(printable) then
                write(nfout,*) ' --- F_WFnUD (down) ---'
                write(nfout,*) F_WFnUD
             end if
             inquire(unit=nfwannier,opened=open)
             if(open) close(nfwannier,status='keep')
             call open0(nfwannier,F_WFnUD,'F_WFnUD   ',unknown,formatted,check_file_name_on)
          end if
       end if
    end if
  end subroutine m_Files_open_nfwannier

  subroutine m_Files_open_nfelf()
    if(mype==0) call open0(nfelf,F_ELF,'F_ELF     ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfelf

  subroutine m_Files_open_nfberry()
    if(mype==0) call open0(nfberry,F_BERRY,'F_BERRY   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfberry

  subroutine m_Files_open_nfeffchg()
    if(mype==0) call open0(nfeffchg,F_EFFCHG,'F_EFFCHG  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfeffchg

!  subroutine m_Files_open_nfforce()
!    if(mype==0) call open0(nfforce,F_FORCE,'F_FORCE   ',unknown,formatted,check_file_name_on)
!  end subroutine m_Files_open_nfforce
  subroutine m_Files_open_nfforce(append)
    logical, intent(in) :: append
    logical :: opd
    if(mype/=0) return

    if(append) then
     call open1(nfforce,F_FORCE,'F_FORCE   ',unknown,formatted,check_file_name_on)
    else
     inquire(file=F_FORCE, opened=opd)
     if(opd) close(nfforce)
     call open0(nfforce,F_FORCE,'F_FORCE   ',unknown,formatted,check_file_name_on)
    endif
  end subroutine m_Files_open_nfforce

! ==== KT_add === 13.1R
  subroutine m_Files_open_nfeps_phonon(append)
    logical, intent(in) :: append
    logical :: opd

    if(mype/=0) return

    if(append) then
       call open1(nfeps_ph,F_EPS_PHONON,'F_EPS_PHONON',unknown,formatted, &
            &     check_file_name_on)
    else
       inquire(file=F_EPS_PHONON, opened=opd)
       if(opd) close(nfeps_ph)
       call open0(nfeps_ph,F_EPS_PHONON,'F_EPS_PHONON',unknown,formatted,&
            &     check_file_name_on)
    endif
  end subroutine m_Files_open_nfeps_phonon

  subroutine m_Files_open_nfoptical_coeff(append)
    logical, intent(in) :: append
    logical :: opd

    if(mype/=0) return
    if(append) then
       call open1(nfoptical_coeff,F_OPTICAL_COEFF,'F_OPTICAL_Coeff',unknown,formatted,&
            &     check_file_name_on)
    else
       inquire(file=F_OPTICAL_COEFF, opened=opd)
       if(opd) close(nfoptical_coeff)
       call open0(nfoptical_coeff,F_OPTICAL_COEFF,'F_OPTICAL_Coeff',unknown,formatted,&
            &     check_file_name_on)
    endif
  end subroutine m_Files_open_nfoptical_coeff

  subroutine m_Files_open_nframan_spectra()
    if (mype==0) call open0( nframan_spectra, F_RAMAN_SPECTRA, 'F_RAMAN_SPECTRA', &
         &                   unknown, formatted, check_file_name_on )
  end subroutine m_Files_open_nframan_spectra
! =============== 13.1R

! ==== KT_add ==== 13.1XI
  subroutine m_Files_open_nf_xi_spectra()
    if (mype==0) call open0( nf_excitation_spectra, F_EXCITATION_SPECTRA, &
         &                   'F_EXCITATION_SPECTRA', &
         &                   unknown, formatted, check_file_name_on )
  end subroutine m_Files_open_nf_xi_spectra

  subroutine m_Files_close_nf_xi_spectra
    close(nf_excitation_spectra, status = 'keep')
  end subroutine m_Files_close_nf_xi_spectra
! ================ 13.1XI

  subroutine m_Files_open_nfmode()
    if(mype==0) call open0(nfmode,F_MODE,'F_MODE    ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfmode

  subroutine m_Files_open_nfepsilon()
    if(mype==0) call open0(nfepsilon,F_EPSILON,'F_EPSILON ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfepsilon

  subroutine m_Files_open_nfoccmat()
    if(mype==0) call open0(nfoccmat,F_OCCMAT,'F_OCCMAT  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfoccmat

  subroutine m_Files_open_nfeps( name_add )
    character*(*), optional, intent(in) :: name_add
    logical open
    character(len=105) :: retstr, retstr0

    if ( .not.present(name_add) ) then
       retstr = F_EPSOUT
    else
       retstr0 = "_" // trim(name_add)
       call insert_str_before_dot( retstr0, F_EPSOUT, retstr )
    endif

    inquire(unit = nfepsout, opened = open)
    if ( .not.open .and. mype==0 ) then
       call open0( nfepsout, retstr, 'F_EPSOUT  ', unknown, &
            &      formatted, check_file_name_on )
    endif
  end subroutine m_Files_open_nfeps

  subroutine m_Files_close_nfeps
    close(nfepsout, status = 'keep')
  end subroutine m_Files_close_nfeps

  subroutine m_Files_open_nfnlo
    logical open
    inquire(unit = nfnlo, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfnlo,F_NLO,'F_NLO     ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfnlo

  subroutine m_Files_close_nfnlo
    close(nfnlo, status = 'keep')
  end subroutine m_Files_close_nfnlo

  subroutine m_Files_open_nfmagopt
    logical open
    inquire(unit = nfmagopt, opened = open)
    if(.not.open.and.mype==0) &
         & call open0(nfmagopt,F_MAGOPT,'F_MAGOPT  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfmagopt

  subroutine m_Files_close_nfmagopt
    close(nfmagopt, status = 'keep')
  end subroutine m_Files_close_nfmagopt

  subroutine m_Files_open_nfepscont
    logical open
    inquire(unit = nfepscont, opened = open)
    if(.not.open.and.mype==0) &
         & call open0(nfepscont,F_EPSCONT,'F_EPSCONT     ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfepscont

  subroutine m_Files_close_nfepscont
    close(nfepscont, status = 'keep')
  end subroutine m_Files_close_nfepscont

!-------------------- T. Hamada 2021.9.28 -------------------- ! UVSOR
  subroutine m_Files_open_nfkptdata
    logical open
    inquire(unit = nfkptdata, opened = open)
    if(.not.open .and. mype==0) &
   & call open0(nfkptdata,F_KPT_DATA,'F_KPT_DATA',unknown,unformatted,check_file_name_on)
  end subroutine m_Files_open_nfkptdata

  subroutine m_Files_close_nfkptdata
    close(nfkptdata, status = 'keep')
  end subroutine m_Files_close_nfkptdata

    subroutine m_Files_open_nfkptepsdata
    logical open
    inquire(unit = nfkptepsdata, opened = open)
    if(.not.open .and. mype==0) &
   & call open0(nfkptepsdata,F_KPTEPS_DATA,'F_KPTEPS_DATA',unknown,unformatted,check_file_name_on)
  end subroutine m_Files_open_nfkptepsdata

  subroutine m_Files_close_nfkptepsdata
    close(nfkptepsdata, status = 'keep')
  end subroutine m_Files_close_nfkptepsdata
! ------------------------------------------------------------ ! UVSOR

  subroutine m_Files_open_nfphdos()
    if(mype==0) call open0(nfphdos,F_PHDOS,'F_PHDOS   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfphdos

  subroutine m_Files_open_nfstrfrc()
    if(mype==0) call open0(nfstrfrc,F_STRFRC,'F_STRFRC  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfstrfrc

  subroutine m_Files_open_nfcntn_wannier
    logical open
    inquire(unit = nfcntn_wannier, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfcntn_wannier,F_CNTN_WAN,'F_CNTN_WAN',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_nfcntn_wannier

  subroutine m_Files_open_nfpot_wannier
    logical open
    inquire(unit = nfpot_wannier, opened = open)
    if(.not.open .and. mype==0) &
         & call open0(nfpot_wannier,F_POT_WAN,'F_POT_WAN ',unknown,unformatted&
         &     ,check_file_name_on)
  end subroutine m_Files_open_nfpot_wannier

  subroutine m_Files_open_nfpstrn
    if(mype == 0) &
         & call open0(nfpstrn,F_PSTRN,'F_PSTRN   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfpstrn

  subroutine m_Files_open_nfPsicoef()
    if(ekmode == ON .or. (ekmode==OFF.and..not.file_existence_contfiles) .or. icond==INITIAL) then
       if(mype==0) call open0(nfpsicoef,F_PsiCoef,'F_PSICOEF ',unknown,formatted,check_file_name_on)
    else if(ekmode==OFF.and.file_existence_contfiles) then
       if(mype==0) call open1(nfpsicoef,F_PsiCoef,'F_PSICOEF ',unknown,formatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfPsicoef

  subroutine m_Files_open_nfBandSymInput()
    if((ekmode==ON .or. (ekmode==OFF.and..not.file_existence_contfiles).or.icond==INITIAL) .and. mype==0) then
       call open0(nfbandsyminput,F_Band_Sym_Input,'F_BAND_SYM_INPUT ',unknown,formatted,check_file_name_on)
    else if((ekmode==OFF.and.file_existence_contfiles) .and. mype==0) then
       call open1(nfbandsyminput,F_Band_Sym_Input,'F_BAND_SYM_INPUT ',unknown,formatted,check_file_name_on)
    end if
  end subroutine m_Files_open_nfBandSymInput

  subroutine m_Files_close_nfPsicoef()
    logical :: op
    if(mype == 0) then
       inquire(unit = nfpsicoef, opened = op)
       if(op) close(nfpsicoef)
    end if
  end subroutine m_Files_close_nfPsicoef

  subroutine m_Files_close_nfBandSymInput()
    logical :: op
    if(mype == 0) then
       inquire(unit = nfbandsyminput, opened = op)
       if(op) close(nfbandsyminput)
    end if
  end subroutine m_Files_close_nfBandSymInput

  subroutine m_Files_open_nfvelec
    if(mype == 0) &
         & call open0(nfvelec,F_VELEC,'F_VELEC   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfvelec

  subroutine m_Files_open_nfeppair
    if(mype == 0) &
         & call open0(nfeppair,F_EPPAIR,'F_EPPAIR  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfeppair

  subroutine m_Files_open_nfeppair2
    if(mype == 0) &
         & call open0(nfeppair2,F_EPPAIR2,'F_EPPAIR2  ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfeppair2

  subroutine m_Files_open_nfvelec_grad
    if(mype == 0) &
         & call open0(nfvelec_grad,F_VELEC_GRAD,'F_VELECGRA',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfvelec_grad

  subroutine m_Files_close_nfpstrn
    logical open
    inquire(unit=nfpstrn,opened=open)
    if(open) close(nfpstrn,status='keep')
  end subroutine m_Files_close_nfpstrn

  subroutine m_Files_close_nfvelec
    logical open
    inquire(unit=nfvelec,opened=open)
    if(open) close(nfvelec,status='keep')
  end subroutine m_Files_close_nfvelec

  subroutine m_Files_close_nfeppair
    logical open
    inquire(unit=nfeppair,opened=open)
    if(open) close(nfeppair,status='keep')
  end subroutine m_Files_close_nfeppair

  subroutine m_Files_close_nfeppair2
    logical open
    inquire(unit=nfeppair2,opened=open)
    if(open) close(nfeppair2,status='keep')
  end subroutine m_Files_close_nfeppair2

  subroutine m_Files_close_nfvelec_grad
    logical open
    inquire(unit=nfvelec_grad,opened=open)
    if(open) close(nfvelec_grad,status='keep')
  end subroutine m_Files_close_nfvelec_grad

  subroutine m_Files_open_nfefermi()
    if(mype == 0) &
         & call open0(nfefermi,F_EFERMI,'F_FERMI   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfefermi

  subroutine m_Files_close_nfefermi()
    logical open
    inquire(unit=nfefermi,opened=open)
    if(open) close(nfefermi,status='keep')
  end subroutine m_Files_close_nfefermi

  subroutine m_Files_open_nfimage(i)
    integer  i
    !!$if(mype==0) call open0(nfimage,F_IMAGE(i),'F_IMAGE   ',unknown,formatted,check_file_name_on)
    call open0(nfimage,F_IMAGE(i),'F_IMAGE   ',unknown,formatted,check_file_name_on)
  end subroutine m_Files_open_nfimage


! ======================== Added by K. Tagami ================= 10.1
  subroutine m_Files_open_nf_LR()
    logical open
    inquire(unit = nf_LR_spectra, opened = open)
    if(.not.open .and. mype==0) &
      & call open0( nf_LR_spectra,F_LR_SPECTRA,'F_LR_SPECTRA   ', &
      &             unknown,formatted,check_file_name_on )
  end subroutine m_Files_open_nf_LR

  subroutine m_Files_close_nf_LR
    close(nf_LR_spectra, status = 'keep')
  end subroutine m_Files_close_nf_LR
! ============================================================= 10.1

! ====== KT_add ========== 13.0S
  subroutine m_Files_open_core_energy_file()
    logical open
    inquire(unit = nfcore_energy_out, opened = open)

    if(.not.open .and. mype==0) &
      & call open0( nfcore_energy_out, F_CORE_ENERGY_OUT, 'F_CORE_ENERGY_OUT', &
      &             unknown,formatted,check_file_name_on )
  end subroutine m_Files_open_core_energy_file

  subroutine m_Files_close_core_energy_file
    close(nfcore_energy_out, status = 'keep')
  end subroutine m_Files_close_core_energy_file

  subroutine m_Files_open_core_ene_initial()
    logical open
    inquire(unit = nfcore_energy_initial, opened = open)

    if(.not.open .and. mype==0) &
      & call open0( nfcore_energy_initial, F_CORE_ENERGY_INITIAL, &
      &            'F_CORE_ENERGY_INITIAL', &
      &             unknown,formatted,check_file_name_on )
  end subroutine m_Files_open_core_ene_initial

  subroutine m_Files_close_core_ene_initial
    close(nfcore_energy_initial, status = 'keep')
  end subroutine m_Files_close_core_ene_initial

  subroutine m_Files_open_core_ene_final()
    logical open
    inquire(unit = nfcore_energy_final, opened = open)

    if(.not.open .and. mype==0) &
      & call open0( nfcore_energy_final, F_CORE_ENERGY_FINAL, &
      &            'F_CORE_ENERGY_FINAL', &
      &             unknown,formatted,check_file_name_on )
  end subroutine m_Files_open_core_ene_final

  subroutine m_Files_close_core_ene_final
    close(nfcore_energy_final, status = 'keep')
  end subroutine m_Files_close_core_ene_final
! ======================== 13.0S

  subroutine m_Files_open_mag_opencore()
    logical open
    inquire(unit = nf_mag_core, opened = open)
    if(.not.open .and. mype==0) &
      & call open0( nf_mag_core, F_MAG_CORE, 'F_MAG_CORE', &
      &             unknown, formatted, check_file_name_on )
  end subroutine m_Files_open_mag_opencore

  subroutine m_Files_close_mag_opencore
    close(nf_mag_core, status = 'keep')
  end subroutine m_Files_close_mag_opencore

! === defect
  subroutine m_Files_open_electrosta_pot()
    logical open
    inquire(unit = nf_elecpot_bin, opened = open)
    if(.not.open .and. mype==0) &
      & call open0( nf_elecpot_bin, F_ELECPOT_BIN, 'F_ELECPOT_BIN', &
      &             unknown, unformatted, check_file_name_on )
  end subroutine m_Files_open_electrosta_pot

  subroutine m_Files_close_electrosta_pot
    close(nf_elecpot_bin, status = 'keep')
  end subroutine m_Files_close_electrosta_pot

  subroutine m_Files_open_electrosta_pot_ref()
    logical open
    inquire(unit = nf_elecpot_bin_ref, opened = open)
    if(.not.open .and. mype==0) &
      & call open0( nf_elecpot_bin_ref, F_ELECPOT_BIN_REF, 'F_ELECPOT_BIN_REF', &
      &             unknown, unformatted, check_file_name_on )
  end subroutine m_Files_open_electrosta_pot_ref

  subroutine m_Files_close_electrosta_pot_ref
    close(nf_elecpot_bin_ref, status = 'keep')
  end subroutine m_Files_close_electrosta_pot_ref
! === defect


! ==== KT_add === 2014/07/14
  subroutine m_Files_open_nfhypervec()
    if (mype==0) then
       call open0( nfhypervec, F_HYPERVEC,'F_HYPERVEC',unknown,formatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open_nfhypervec

  subroutine m_Files_close_nfhypervec()
    logical open
    inquire(unit=nfhypervec,opened=open)
    if(open) close(nfhypervec,status='keep')
  end subroutine m_Files_close_nfhypervec
! =============== 2014/07/14

  subroutine m_Files_open_nf_ekindens()
    logical open
    inquire(unit = nf_ekindens, opened = open)
    if (.not.open .and. mype==0) &
         & call open0( nf_ekindens, F_EKINDENS, 'F_EKINDENS',  unknown, unformatted &
         &      ,check_file_name_off )
  end subroutine m_Files_open_nf_ekindens

  subroutine m_Files_close_nf_ekindens()
    if(mype==0) close(nf_ekindens)
  end subroutine m_Files_close_nf_ekindens

  subroutine m_Files_open(nf,F_IDEN,chara)
    integer, intent(in) :: nf
    character(len=stringlength_of_filenames), intent(in) :: F_IDEN
    character(len=*), intent(in):: chara
    if (mype==0)then
       call open0( nf, F_IDEN,chara,unknown,formatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open

  subroutine m_Files_close(nf)
    integer, intent(in) :: nf
    logical open
    inquire(unit=nf,opened=open)
    if(open) close(nf,status='keep')
  end subroutine m_Files_close

  subroutine m_Files_open_nfpos()
    if (mype==0) then
       call open0( nfpos, F_POS,'F_POS',unknown,formatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open_nfpos

  subroutine m_Files_close_nfpos()
    logical open
    inquire(unit=nfpos,opened=open)
    if(open) close(nfpos,status='keep')
  end subroutine m_Files_close_nfpos

  subroutine m_Files_open_nfcoordattr()
    if (mype==0) then
       call open0( nfcoordattr, F_COORD_ATTR,'F_COORD_ATTR',unknown,formatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open_nfcoordattr

  subroutine m_Files_close_nfcoordattr()
    logical open
    inquire(unit=nfcoordattr,opened=open)
    if(open) close(nfcoordattr,status='keep')
  end subroutine m_Files_close_nfcoordattr

  subroutine m_Files_open_nfhybridinfo()
    if (mype==0) then
       call open0( nfhybridinfo, F_HYBRIDINFO,'F_HYBRIDINFO',unknown,unformatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open_nfhybridinfo

  subroutine m_Files_close_nfhybridinfo()
    logical open
    inquire(unit=nfhybridinfo,opened=open)
    if(open) close(nfhybridinfo,status='keep')
  end subroutine m_Files_close_nfhybridinfo

  subroutine m_Files_open_nfscfzaj()
    if (mype==0) then
       call open0( nfscfzaj, F_SCF_ZAJ,'F_SCF_ZAJ ',unknown,unformatted,&
            &      check_file_name_on )
    endif
  end subroutine m_Files_open_nfscfzaj

  subroutine m_Files_close_nfscfzaj()
    logical open
    inquire(unit=nfscfzaj,opened=open)
    if(open) close(nfscfzaj,status='keep')
  end subroutine m_Files_close_nfscfzaj

  subroutine m_Files_checkpoint_dir(update_counter)
    use m_Control_Parameters, only : cpt_nhistory,cpt_history_count
    use m_IterationNumbers, only : iteration_ionic
    logical, intent(in) :: update_counter
    character(len=10) :: cid1,cidmax,cid2
    integer :: i,nn,ierr
    if (cpt_nhistory<2) return
    if (cpt_history_count==0)then
        if(update_counter) cpt_history_count = cpt_history_count + 1
        return
    endif

    if (cpt_history_count<cpt_nhistory)then
      cid1 = get_count_char(cpt_nhistory,cpt_history_count)
      call system('mkdir -p '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1)))
      cidmax = cid1
    else
      cidmax = get_count_char(cpt_nhistory,cpt_nhistory-1)
      nn = cpt_nhistory-1
      if(mype==0 .and. nn>=2)then
        do i=2,nn
          cid1 = get_count_char(cpt_nhistory,i)
          cid2 = get_count_char(cpt_nhistory,i-1)
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_CHGT//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_ZAJ//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_CNTN//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_CNTN_BIN//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          if (m_Files_nfcntn_bin_paw_exists()) then
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_CNTN_BIN_PAW//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          endif
          if (m_Files_nfoccmat_exists()) then
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_OCCMAT//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
          endif
        enddo
      endif
    endif

    if(mype==0)then
    call system('cp -f '//F_CHGT//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    call system('cp -f '//F_ZAJ//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    call system('cp -f '//F_CNTN//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    call system('cp -f '//F_CNTN_BIN//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    if (m_Files_nfcntn_bin_paw_exists()) then
      call system('cp -f '//F_CNTN_BIN_PAW//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    endif
    if (m_Files_nfoccmat_exists()) then
      call system('cp -f '//F_OCCMAT//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')
    endif
    endif

    if(update_counter) cpt_history_count = cpt_history_count+1
    call mpi_barrier(MPI_CommGroup,ierr)

  end subroutine m_Files_checkpoint_dir

  subroutine m_Files_checkpoint_dir_neb()
    use m_Control_Parameters, only : cpt_nhistory,cpt_history_count
    use m_IterationNumbers, only : iteration_ionic
    character(len=10) :: cid1,cidmax,cid2
    integer :: i,nn,ierr
    if (cpt_nhistory<2) return
    if (cpt_history_count==0)then
        cpt_history_count = cpt_history_count + 1
        return
    endif

    if (cpt_history_count<cpt_nhistory)then
      cid1 = get_count_char(cpt_nhistory,cpt_history_count)
      call system('mkdir -p '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1)))
      cidmax = cid1
    else
      cidmax = get_count_char(cpt_nhistory,cpt_nhistory-1)
      nn = cpt_nhistory-1
      if(mype==0 .and. nn>=2)then
        do i=2,nn
          cid1 = get_count_char(cpt_nhistory,i)
          cid2 = get_count_char(cpt_nhistory,i-1)
          call system('cp -f '//trim(adjustl(F_CHKPNT))//trim(adjustl(cid1))//'/'//F_NEB_CNTN//' '// &
                        &  trim(adjustl(F_CHKPNT))//trim(adjustl(cid2))//'/')
        enddo
      endif
    endif

    if(mype==0) call system('cp -f '//F_NEB_CNTN//' '//trim(adjustl(F_CHKPNT))//trim(adjustl(cidmax))//'/')

    cpt_history_count = cpt_history_count+1
    call mpi_barrier(MPI_CommGroup,ierr)

  end subroutine m_Files_checkpoint_dir_neb

  character(len=10) function get_count_char(maxcount,currcount)
    integer, intent(in) :: maxcount,currcount
    integer :: iketa=2
    integer :: itmp
    character(len=256) :: ch
    character(len=256) :: string
    itmp=int(log10(real(maxcount)))+1
    if(itmp>2)iketa=itmp
    write(ch,*) iketa
    write(get_count_char,'(i'//trim(adjustl(ch))//'.'//trim(adjustl(ch))//')') currcount
  end function get_count_char

  subroutine m_Files_resolve_default_pp(ntyp,iatomn,xctype)
    implicit none
    integer, intent(in) :: ntyp
    real(kind=DP), dimension(ntyp), intent(in) :: iatomn
    character(len=*), intent(in) :: xctype
    integer, parameter :: max_elements = 120
    integer :: ityp
    logical :: exi
    integer :: icount
    character(len=stringlength_of_filenames) :: pppath
    character(len=stringlength_of_filenames), dimension(max_elements) :: ggapbe,ldapw91
    namelist/defaultpp/ggapbe,ldapw91
    icount = 0
    do ityp=1,ntyp
       inquire(file=F_POT(ityp),exist=exi)
       if (exi) icount = icount+1
    enddo
    if(icount==ntyp) return

    icount = 0
    do ityp=1,ntyp
       if (F_POT(ityp) /= F_POT_ORG(ityp)) icount = icount+1
    enddo
    if(icount==ntyp) return

    call getenv('PHASE_PP_PATH',pppath)
    if(len(pppath)==0) return

    include 'defaultppfiles'
    inquire(file=trim(F_DEFAULTPP),exist=exi)
    if(exi)then
      open(nfdefaultpp,file=trim(F_DEFAULTPP))
      read(nfdefaultpp,NML = defaultpp, err=30, end=10)
10    continue
      close(nfdefaultpp)
    endif

    do ityp=1,ntyp
       inquire(file=trim(F_POT(ityp)),exist=exi)
       if(.not.exi) then
          F_POT(ityp) = trim(pppath)//'/'//trim(ggapbe(int(iatomn(ityp))))
          if (xctype == 'ldapw91')then
             F_POT(ityp) = trim(pppath)//'/'//trim(ldapw91(int(iatomn(ityp))))
          endif
          write(nfout,'(a,i3,a)') ' !** default PP file adopted for element no. '&
         &                        ,ityp,' : '//trim(F_POT(ityp))
       endif
    enddo

    return

30  call phase_execution_error(FILENAMES_FORMAT_ERROR)
  end subroutine m_Files_resolve_default_pp

  subroutine m_Files_open_pwbs()
    logical open
    inquire(unit = nfpwbs, opened = open)
    if(.not.open .and. mype == 0) &
         & call open0(nfpwbs,F_PWBS,'F_PWBS    ',unknown,unformatted&
         &     , check_file_name_on)
  end subroutine m_Files_open_pwbs

  subroutine m_Files_close_pwbs()
    if(mype==0) close(nfpwbs)
  end subroutine m_Files_close_pwbs

  logical function m_Files_file_exists(fname)
    character(len=*) :: fname
    logical :: exi
    inquire(file=trim(fname),exist=exi)
    m_Files_file_exists = exi
  end function m_Files_file_exists

  subroutine m_Files_open_core_charge_dens(it)
    integer, intent(in) :: it
    integer :: n1

    logical open
    inquire(unit = nf_core_charge_radial(it), opened = open)
    if( .not.open) &
         & call open0( nf_core_charge_radial(it), F_CORE_CHARGE(it), &
         &             'F_CORE_CHARGE ', unknown, formatted,  &
         &             check_file_name_on )

  end subroutine m_Files_open_core_charge_dens

  subroutine m_Files_close_core_charge_dens(it)
    integer, intent(in) :: it

    close( nf_core_charge_radial(it) )
  end subroutine m_Files_close_core_charge_dens

  subroutine m_Files_open_xsf()
    use m_IterationNumbers, only : iteration_ionic
    character(len=105) :: F_XSF_id
    integer, parameter :: len_nnumber = 6
    character(len=len_nnumber) :: nnumber
    integer :: i,ip,lenstr,id
    integer :: max_nfiles = 999999
    logical open,exi,found_new_file

    if(mype==0)then
      lenstr = len_trim(F_XSF)
      ip = lenstr+1
      findingdot: do i = lenstr, 1, -1
        if(F_CHR(i:i) == '/') then
          exit findingdot
        else if(F_CHR(i:i) == '.') then
          ip = i
          exit findingdot
        end if
      end do findingdot
      write(nnumber,'(i6)') iteration_ionic
      do i = 1, len_nnumber
         if(nnumber(i:i) == ' ') nnumber(i:i) = '0'
      end do
      if(ip >= 1 .and. ip <= lenstr) then
        if(ip+1 <= lenstr) then
          F_XSF_id = F_XSF(1:ip)//nnumber//F_XSF(ip+1:lenstr)
        else
          F_XSF_id = F_XSF(1:ip)//nnumber
        end if
      else
        F_XSF_id = F_XSF(1:lenstr)//nnumber
      end if
      inquire(unit = nfxsf, opened = open)
      if(.not.open .and. mype == 0) &
           & call open0(nfxsf,F_XSF_id,'F_XSF',unknown,formatted&
           &     , check_file_name_on)
    endif
  end subroutine m_Files_open_xsf

  subroutine m_Files_close_xsf()
    if(mype==0) close(nfxsf)
  end subroutine m_Files_close_xsf

  subroutine m_Files_open_n2p2()
    logical open
    inquire(unit = nfn2p2, opened = open)
    if(.not.open .and. mype == 0) &
         & call open1(nfn2p2,F_N2P2,'F_N2P2',unknown,formatted&
         &     , check_file_name_on)
  end subroutine m_Files_open_n2p2

  subroutine m_Files_close_n2p2()
    if(mype==0) close(nfn2p2)
  end subroutine m_Files_close_n2p2

  subroutine m_Files_open_file(nfile, fname, tagname)
    integer, intent(in) :: nfile
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: tagname
    logical open
    inquire(unit = nfile, opened = open)
    if(.not.open .and. mype == 0) &
         & call open1(nfile,fname,tagname,unknown,formatted&
         &     , check_file_name_on)
  end subroutine m_Files_open_file

  subroutine m_Files_close_file(nfile)
    integer, intent(in) :: nfile
    if(mype==0) close(nfile)
  end subroutine m_Files_close_file

end module m_Files
