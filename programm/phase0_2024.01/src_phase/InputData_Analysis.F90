!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: cnstr_fcvect_work_alloc, cnstr_fcvect_work_dealloc,
!             get_CS_and_ionic_system_data, read_ntyp_natm_natm2, read_altv,
!             read_ncord_ninv, read_atomic_coordinates,
!             read_iatomn_alfa_amion_etc, read_nbztyp, InputData_Analysis
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
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
subroutine InputData_Analysis()
! $Id: InputData_Analysis.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Files, only :             nfinp, nfout, nfcntn, file_existence_contfiles, file_existence_3contfiles &
       &                        ,file_existence_contfiles_when_paw_on &
       &                        , m_Files_open_nfcntn &
       &                        , m_Files_check_file_names &
       &                        , m_Files_reopen_nfinp &
       &                        , m_Files_reopen_nfcntn &
       &                        , F_ZAJ_in, F_CNTN_in, F_CNTN_BIN_in &
       &                        , F_CHGT_in, F_ZAJ, F_CNTN, F_CNTN_BIN &
       &                        , F_CHGT, F_CNTN_bak &
       &                        , F_CHGT_bak, F_CNTN_BIN_bak &
       &                        , F_ZAJ_bak &
       &                        , m_Files_check_file_existence &
       &                        , m_Files_check_nfcntn_existence, m_Files_echo_nfinp &
       &                        , m_Files_open_nfdynm_cif_initially &
       &                        , m_Files_open_pcont_files_ini &
       &                        , m_Files_swap_F_INP, nfpos,nfcoordattr,F_COORD_ATTR,F_POS &
       &                        , m_Files_resolve_default_pp
  use m_Crystal_Structure, only : altv,rltv,univol,rvol&
       &                        , inversion_symmetry&
       &                        , nbztyp,nbztyp_spg &
       &                        , m_CS_set_nbztyp_spg, m_CS_set_nbztyp &
       &                        , m_CS_rd_symmetry_op_etc &
       &                        , m_CS_set_default_symm_op &
       &                        , m_CS_set_abccacbcc &
       &                        , m_CS_rd_n, m_CS_rd_fix_spin_status
  use m_Ionic_System, only :      input_coordinate_system &
       &                        , ntyp, natm, natm2, nfcatm, iwei, imdtyp  &
       &                        , ityp,cnst_typ, iatomn,iatom,ivan &
       &                        , pos, cps, alfa, amion,zeta1 &
       &                        , m_IS_alloc_pos_and_v, m_IS_rd_pos_and_v &
       &                        , m_IS_alloc_iatomn_etc &
       &                        , m_IS_alloc_cnstrvectors_etc &
       &                        , m_IS_init_cnstrnt &
       &                        , m_IS_cp_works2fcvect_etc &
       &                        , m_IS_rd_n, m_IS_set_iatom &
       &                        , m_IS_rd_nrsv_stdin, m_IS_rd_T_parameters &
       &                        , m_IS_set_ionic_mass &
       &                        , m_IS_rd_diis_history &
       &                        , m_IS_natm_can_change &
       &                        , m_IS_rd_curr_atom_reservoir, m_IS_rd_n_pre, m_IS_rd_moved_distances_of_planes &
       &                        , coord_method, filetype, iatomn, ntyp &
       &                        , m_IS_rd_rigid_body
  use m_Control_Parameters, only: imdalg, paramset,icond,icond_org,nspin &
       &                        , iconvergence_previous_job, ekmode, printable &
       &                        , sw_phonon, sw_positron &
       &                        , driver, xctype &
       &                        , sw_phonon, sw_wannier90, PAW_switch &
       &                        , fixed_charge_k_parallel, neg_is_given &
       &                        , m_CtrlP_rd_parameters, m_CtrlP_rd_iconvergence &
       &                        , m_CtrlP_rd_iconv_ek &
       &                        , m_CtrlP_rd_numk_zajsaved &
       &                        , m_CtrlP_check_inputfilestyle &
       &                        , m_CtrlP_rd_control &
       &                        , m_CtrlP_rd_accuracy &
       &                        , m_CtrlP_rd_accuracy_paw_switch &     ! T. Y. 2020/02/28
       &                        , m_CtrlP_rd_accuracy_hubbard_switch & ! T. Y. 2020/02/28
       &                        , m_CtrlP_rd_driver &
#ifndef _EMPIRICAL_
       &                        , m_CtrlP_dealloc_wfsolver &
       &                        , m_CtrlP_rd_wfsolver &
       &                        , m_CtrlP_rd_chargemix &
       &                        , m_CtrlP_check_matm &
       &                        , m_CtrlP_rd_edelta_ontheway, m_CtrlP_rd_neg_previous &
#endif
       &                        , m_CtrlP_rd_struc_evol &
       &                        , m_CtrlP_rd_printlevel &
       &                        , m_CtrlP_rd_postproc &
       &                        , m_CtrlP_rd_corecharge_cntnbin, m_CtrlP_check_dos_method &
       &                        , sw_optimize_lattice ,sw_rebuild_pws,ipriinputfile &
       &                        , m_CtrlP_set_default_gmaxp, gmaxp_defined
#ifndef _EMPIRICAL_
  use m_Kpoints,            only: m_Kp_rd_way_ksample_etc, m_Kp_rd_n, m_Kp_rd_kv3
  use m_ES_dos,             only: m_ESdos_rd_pdos_param, m_ESdos_rd_pcohp_param
#endif
  use m_Const_Parameters, only :DP, PUCV, CARTS, PAI2 &
       &, FIX, RELAX, BONDLENGTH_FIX, BONDLENGTH_FIX_1, BONDLENGTH_FIX_2 &
       &, COG_CNTR, COG_FIX_L, HEAT_BATH &
       &, FIX_IN_A_PLANE, CONTINUATION, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION,DELTA &
       &, GENERAL, GENERAL_LARGER &
       &, OLD, NEW_, GRID, ON, OFF &
       &, GDIIS, BFGS, L_BFGS &
       &, INITIAL, PREPARATION_ONLY, DRIVER_CONSTRAINT, DRIVER_NEB, DRIVER_MTD &
       &, COORDINATE_CONTINUATION, ONE_BY_ONE, ALL_AT_ONCE &
       &, P_CONTROL, PT_CONTROL &
       &, FILE, PHASE0_INPUT, PHASE0_OUTPUT
  use m_IterationNumbers, only : m_Iter_rd_iteration_numbers,m_Iter_rd_iters_and_nk,iteration_unit_cell
#ifndef _EMPIRICAL_
  use m_Phonon,     only : m_Phonon_rd_param, m_Phonon_set_supercell
  use m_BerryPhase, only : m_BP_rd_parameters
  use m_PseudoPotential,  only : m_PP_rd_window_param

  use m_Parallelization, only : mype_conf, nrank_conf, mype, npes, MPI_CommGroup, conf_para

  use m_PAW_ChargeDensity,  only : m_PAWCD_rd_ntheta_nphi
  use m_Wannier90, only : m_Wan90_init
  use m_UnitCell, only : m_UC_rd_cntn_data, m_UC_rd_md_cntn_data

!!$  use m_Parallel,   only : m_Parallel_rd_npes_etc
#endif
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
  use m_Control_Parameters, only: m_CtrlP_rd_dtim_previous
! ==============================================================================

! ================================ added by K. Tagami ============== 11.0&13.0XX
  use m_Control_Parameters,    only : noncol, ndim_spinor, m_CtrlP_chkif_metagga
! =================================================================== 11.0&13.0XX

  use m_ES_occup, only : m_ES_rd_occ_ext
  use m_ErrorMessages

! ==== KT_add === 13.1XI
  use m_Excitation,  only : m_XI_read_input
! =============== 13.1XI

! === KT_add === 2014/12/29
  use m_Ionic_System, only :  m_IS_set_mag_moment0_atoms, m_IS_set_ionic_charge_atoms, &
       &                      sw_set_initial_magmom_by_atom
! ============== 2014/12/29

  use m_ES_ChargeState,  only : m_ESCS_read_param
  use m_Potential_Average,  only : m_PotAvg_read_param

  use m_Epsilon_ek, only : m_CtrlP_rd_epsilon

#ifdef LIBXC
  use m_Control_Parameters,    only : m_CtrlP_init_xctype_libxc
#endif

  implicit none

  real(kind=DP), allocatable, dimension(:)       :: work
  integer, parameter                         :: NWK = 13
  integer                                    :: nlines
  integer,       allocatable, dimension(:)   :: ia_cnst_work    ! d(natm)
  real(kind=DP), allocatable, dimension(:,:) :: fcvect_work     ! d(natm,4)
  integer :: ierror

  integer,      parameter   :: len_str = 132
  character(len=len_str)       :: str
  character(len=260)           :: zajbuf,cntnbuf,cntn_binbuf,chgtbuf
  character(len=256) :: to_string

  integer :: iret, f_closeInputFile
  logical :: ex
  logical,save :: first_call = .true.
  logical :: initialization_required, inputfilestyle_is_new = .false.,isSCDFT
  logical :: file_existence_contfiles_t, has_uspp

  if(.not.initialization_required()) return

  if(m_CtrlP_check_inputfilestyle(nfinp) == NEW_) then
     if(printable) write(nfout,'(" !*--- input-file style = NEW")')
     inputfilestyle_is_new = .true.
     call m_Files_reopen_nfinp(1)
     if(first_call)then
        call m_CtrlP_rd_printlevel(nfout)
        if(ipriinputfile>0) call m_Files_echo_nfinp()
        call m_CtrlP_rd_driver(nfout)
        call m_CtrlP_rd_accuracy_paw_switch(nfout)     ! T. Y. 2020/02/28
        call m_CtrlP_rd_accuracy_hubbard_switch(nfout) ! set sw_hubbard T. Y. 2020/02/28
        file_existence_contfiles_t = file_existence_contfiles
        if(PAW_switch == ON) file_existence_contfiles_t = file_existence_contfiles_when_PAW_on
        call m_CtrlP_rd_control(nfout,file_existence_contfiles_t,file_existence_3contfiles)  ! -> icond, cpumax, nmd2, ifstop, ipriekzaj
     endif
     call m_Files_open_nfdynm_cif_initially()
     if(first_call) call m_CtrlP_rd_accuracy(nfout)   ! -> neg_is_given
     call m_Phonon_rd_param(nfout)
     call m_IS_rd_n_pre(nfout)
     if(coord_method == FILE .and. filetype ==PHASE0_INPUT)then
       call m_Files_swap_F_INP(0,nfpos,F_POS)
     else if (coord_method == FILE .and. filetype == PHASE0_OUTPUT) then
       call m_Files_swap_F_INP(0,nfcoordattr,F_COORD_ATTR)
     endif
     call m_CS_rd_n(nfout)
     if(coord_method == FILE .and. filetype ==PHASE0_INPUT)then
       call m_Files_swap_F_INP(1,nfpos,F_POS)
     else if (coord_method == FILE .and. filetype == PHASE0_OUTPUT) then
       call m_Files_swap_F_INP(1,nfcoordattr,F_COORD_ATTR)
     endif
     call m_CtrlP_rd_struc_evol(nfout)
     call m_CtrlP_rd_postproc(nfout)
     call m_CtrlP_rd_epsilon(nfout)
#ifndef _EMPIRICAL_
     if(sw_wannier90==ON) call m_Wan90_init(nfout)
#endif
!     call m_CS_rd_n(nfout)
     if(coord_method == FILE .and. filetype ==PHASE0_INPUT)then
       call m_Files_swap_F_INP(0,nfpos,F_POS)
     else if (coord_method == FILE .and. filetype == PHASE0_OUTPUT) then
       call m_Files_swap_F_INP(0,nfcoordattr,F_COORD_ATTR)
     endif
     call m_IS_rd_n(nfout)
     if(coord_method == FILE .and. filetype ==PHASE0_INPUT)then
       call m_Files_swap_F_INP(1,nfpos,F_POS)
     else if (coord_method == FILE .and. filetype == PHASE0_OUTPUT) then
       call m_Files_swap_F_INP(1,nfcoordattr,F_COORD_ATTR)
     endif
     call m_IS_set_ionic_mass(nfout)
     call m_Files_resolve_default_pp(ntyp,iatomn,xctype)
     if(.not.gmaxp_defined) call m_CtrlP_set_default_gmaxp(nfout,has_uspp())

     call m_Files_open_pcont_files_ini()

#ifndef _EMPIRICAL_
     if(neg_is_given) then
        call m_CtrlP_rd_wfsolver(nfout,natm2)
        call m_CtrlP_rd_chargemix(nfout)
     end if
#endif

! === KT_add ==== 2014/12/29
     if ( sw_set_initial_magmom_by_atom == ON ) then
        call m_IS_set_ionic_charge_atoms(nfout)
        call m_IS_set_mag_moment0_atoms(nfout)
     endif
! =============== 2014/12/29

#ifndef _EMPIRICAL_
     call m_PP_rd_window_param(nfout)
     call m_ESdos_rd_pdos_param(nfout)
     call m_ESdos_rd_pcohp_param(nfout)
     call m_CtrlP_check_matm(nfout,natm) ! rmm_precal_phase_matm

! === KT_add ===== 13.0XX
#ifdef LIBXC
     if ( xctype == "libxc") then
        call m_CtrlP_init_xctype_libxc( nfout )
     else
        call m_CtrlP_chkif_metagga(nfout)
     endif
#else
     call m_CtrlP_chkif_metagga(nfout)
#endif
! ================ 13.0XX

! === KT_add ==== 13.1XI
     call m_XI_read_input
! =============== 13.1XI

     if(ekmode /= GRID) then
        call m_Kp_rd_n(nfout)
        call m_BP_rd_parameters
     end if
#endif

     call m_ESCS_read_param
     call m_PotAvg_read_param

     call m_CtrlP_check_dos_method(nfout)

     call m_PAWCD_rd_ntheta_nphi

     call m_ES_rd_occ_ext(nfout)

     ! setup supercell

#ifndef _EMPIRICAL_
     if(sw_phonon == ON) call m_Phonon_set_supercell
#endif

     icond_org=icond
     if(driver == DRIVER_CONSTRAINT .and. conf_para .and. .not.(icond==PREPARATION_ONLY)) then
       icond=INITIAL
     endif

     if(driver == DRIVER_MTD .and.mype_conf>0)then
        zajbuf=F_ZAJ;cntnbuf=F_CNTN;cntn_binbuf=F_CNTN_BIN;chgtbuf=F_CHGT
        F_ZAJ=trim(zajbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CNTN=trim(cntnbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CNTN_BIN=trim(cntn_binbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CHGT=trim(chgtbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))

        F_ZAJ_in=trim(zajbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CNTN_in=trim(cntnbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CNTN_BIN_in=trim(cntn_binbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
        F_CHGT_in=trim(chgtbuf)//'_conf'//trim(adjustl(to_string(mype_conf,2)))
     endif

!!$     call parse_args()
!     call m_Files_reopen_nfinp(2)
!!$     stop '  !! stop after [call m_Files_reopen_nfinp(2) <<InputData_Analysis>>'
  else
     if(printable) write(nfout,'(" !! --- input-file style = OLD")')
     allocate(work(NWK))
     call m_CtrlP_rd_parameters(nfinp,nfout,nlines)
     call get_CS_and_ionic_system_data(nlines) !-(contained here) CS: Crystal Structure, contained here

#ifndef _EMPIRICAL_
     if(ekmode /= GRID) call m_Kp_rd_way_ksample_etc(nfinp,nlines) ! <-(nfinp-file)
#endif
     call m_CS_set_abccacbcc()
     if(nbztyp >= 2 .and. &
          & (nbztyp_spg == GENERAL .or. nbztyp_spg == GENERAL_LARGER)) then
!!$        call m_Files_open_nfspg()
!!$        call m_CS_rd_symmetry_op_etc(nfspg)
        call m_CS_rd_symmetry_op_etc(nfinp)
     else
        call m_CS_set_default_symm_op(nfout)  ! <-- nbztyp_spg
     end if
!!$     call m_CS_rd_symmetry_operations(nfinp)

     call m_IS_rd_nrsv_stdin(nfinp)   ! ->(nrsv)
     call m_IS_rd_T_parameters(imdalg,nfinp)
     !                                            ->(qmass,tkb,cprv,cpqr,forcp)
     deallocate(work)
  end if

!!$  call m_Files_check_file_existence()
  ex = file_existence_contfiles.and.file_existence_3contfiles

  if(driver==DRIVER_MTD.and.mype_conf>0 .and. .not.ex.and.icond==CONTINUATION)then
     F_CNTN=F_CNTN_bak
     F_CNTN_BIN=F_CNTN_BIN_bak
     F_CHGT=F_CHGT_bak
     F_ZAJ=F_ZAJ_bak
     F_CNTN_in=F_CNTN_bak
     F_CNTN_BIN_in=F_CNTN_BIN_bak
     F_CHGT_in=F_CHGT_bak
     F_ZAJ_in=F_ZAJ_bak
     if(printable) then
        write(nfout,'(a)') 'temporarily renamed continuation files'
        write(nfout,'(a)') '  F_CNTN     : '//trim(F_CNTN)
        write(nfout,'(a)') '  F_CNTN_BIN : '//trim(F_CNTN_BIN)
        write(nfout,'(a)') '  F_CHGT     : '//trim(F_CHGT)
        write(nfout,'(a)') '  F_ZAJ      : '//trim(F_ZAJ)
     endif
  endif

  if(driver .ne. DRIVER_NEB) then
  if(icond == CONTINUATION .or.icond==COORDINATE_CONTINUATION .or. &
       & (icond == FIXED_CHARGE_CONTINUATION .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
     if(sw_phonon == OFF) then
        if(m_Files_check_nfcntn_existence()) then
           call m_Files_open_nfcntn()
        else
           ierror = FILE_NOT_EXIST
           call phase_error(ierror, nfout, nfcntn, F_CNTN, __LINE__, __FILE__)
        end if
     endif
  endif
#ifndef DISABLE_CONSTRAINTS
  if(driver==DRIVER_CONSTRAINT) call constrained_dynamics_init(nfout)
#endif
  if(icond == CONTINUATION .or.icond==COORDINATE_CONTINUATION .or. &
       & (icond == FIXED_CHARGE_CONTINUATION .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
   if(sw_phonon == OFF) then
     if(m_Files_check_nfcntn_existence()) then
        call m_Files_open_nfcntn()
     else
        ierror = FILE_NOT_EXIST
        call phase_error(ierror, nfout, nfcntn, F_CNTN, __LINE__, __FILE__)
     end if
     call m_Iter_rd_iteration_numbers(nfcntn,icond)
     call m_IS_rd_pos_and_v(nfcntn)
     call m_IS_rd_rigid_body(nfcntn)

     call m_IS_rd_curr_atom_reservoir(nfcntn)
     if ((imdalg==GDIIS .or. imdalg==BFGS .or. imdalg==L_BFGS).and.icond.ne.COORDINATE_CONTINUATION) then
        call m_IS_rd_diis_history(nfcntn)
     endif
#ifndef _EMPIRICAL_
     if(icond==CONTINUATION) call m_CtrlP_rd_neg_previous(nfcntn,nfout)
     call m_CtrlP_rd_edelta_ontheway(nfcntn)
#endif
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
   call m_CtrlP_rd_dtim_previous(nfcntn)
! ==============================================================================
     call m_CtrlP_rd_iconvergence(nfcntn)
     call m_CtrlP_rd_corecharge_cntnbin(nfcntn) ! -> status_cntnbin_positron
     call m_CS_rd_fix_spin_status(nfcntn,nfout)  ! <-- T. Yamasaki, 18th Aug. 2009
     if(icond==CONTINUATION.and.sw_optimize_lattice==ON.and.sw_rebuild_pws==ON) &
      & call m_UC_rd_cntn_data(nfcntn)
     if(icond==CONTINUATION.and.(imdalg==P_CONTROL .or. imdalg==PT_CONTROL).and.sw_rebuild_pws==ON) &
      & call m_UC_rd_md_cntn_data(nfcntn)
     if(isSCDFT()) call rd_eps0(nfcntn)
     call m_IS_rd_moved_distances_of_planes(nfcntn)
   end if
  else if(icond == FIXED_CHARGE .and. fixed_charge_k_parallel == ALL_AT_ONCE) then
     !!if( m_Files_check_nfcntn_existence()) then
     !!   call m_Files_open_nfcntn()
     !!   call m_IS_rd_pos_and_v(nfcntn)
     !!endif
  else if(icond == FIXED_CHARGE_CONTINUATION.and.fixed_charge_k_parallel == ONE_BY_ONE) then
     if(m_Files_check_nfcntn_existence()) then
        call m_Files_open_nfcntn()
     else
        ierror = FILE_NOT_EXIST
        call phase_error(ierror, nfout, nfcntn, F_CNTN, __LINE__, __FILE__)
     end if
     call m_CtrlP_rd_iconvergence(nfcntn)

! ================================================ modified by K. Tagami ======= 11.0
!     call m_Iter_rd_iters_and_nk(nfcntn,nspin)
!!$     if(fixed_charge_k_parallel == ALL_AT_ONCE) then
!!$        call m_Iter_rd_iteration_numbers(nfcntn,icond)
!!$     else if(fixed_charge_k_parallel == ONE_BY_ONE) then
     if ( noncol ) then
        call m_Iter_rd_iters_and_nk(nfcntn,ndim_spinor,nfout)
     else
        call m_Iter_rd_iters_and_nk(nfcntn,nspin,nfout)
     endif
! ============================================================================== 11.0

!!$     call m_Parallel_rd_npes_etc(nfcntn)
     call m_Kp_rd_kv3(nfcntn,nfout)
     if(sw_positron /= OFF) call m_CtrlP_rd_corecharge_cntnbin(nfcntn) ! -> status_cntnbin_positron
     call m_CtrlP_rd_iconv_ek(nfcntn,nfout)  ! -> iconv_ek_tmp
     call m_CtrlP_rd_numk_zajsaved(nfcntn,nfout) ! -> numk_zajsaved
  end if
  end if

  !!$call m_Files_check_file_names  ! for debugging

  iret = f_closeInputFile()
  first_call = .false.

  if(.not.neg_is_given) then
     call Check_of_Pseudopotential()  !-(PseudoPotential_Construction.F90), set neg properly

     if(inputfilestyle_is_new) then
        call m_Files_reopen_nfinp(1)
#ifndef _EMPIRICAL_
!!$        call m_CtrlP_dealloc_wfsolver()
        call m_CtrlP_rd_wfsolver(nfout,natm2)
        call m_CtrlP_rd_chargemix(nfout)
#endif
        iret = f_closeInputFile()
     else
        if(printable) then
           write(nfout,'(" !** inputfilestyle = OLD, and you cannot omit number_of_bands in the inputfile")')
        end if
        call phase_error_with_msg(nfout, 'inputfilestyle = OLD, and you cannot omit number_of_bands in the inputfile',&
       & __LINE__,__FILE__)
        !stop  'inputfilestyle = OLD, and you cannot omit number_of_bands in the inputfile'
     end if
  end if
contains
  subroutine cnstr_fcvect_work_alloc
    if(printable) write(nfout,'(" ! natm  = ",i6," << cnstr_fcvect_work_alloc >>")') natm
    allocate(ia_cnst_work(natm))
    allocate(fcvect_work(natm,4));fcvect_work = 0.d0
  end subroutine cnstr_fcvect_work_alloc

  subroutine cnstr_fcvect_work_dealloc
    if(allocated(fcvect_work)) deallocate(fcvect_work)
    if(allocated(ia_cnst_work)) deallocate(ia_cnst_work)
  end subroutine cnstr_fcvect_work_dealloc

  subroutine get_CS_and_ionic_system_data(nlines)
    integer, intent(in) :: nlines
!        =================
    rewind nfinp
    call read_ntyp_natm_natm2(nfinp)       ! gmax, gmaxp, ntyp, natm, natm2
    call read_altv(nfinp)                  ! altv, rltv, univol, rvol
    call read_ncord_ninv(nfinp)            ! ncord, ninv

    ! --> atomic coordinates, iwei, imdtyp, ityp, fcvect etc.
    call m_IS_alloc_pos_and_v(nfout)       ! pos, cps, cpd_l, cpo_l, iwei, imdtyp, ityp
    call cnstr_fcvect_work_alloc()          !-(InputData_Analysis) (*) ->ia_cnst_work,fcvect_work
    call read_atomic_coordinates(nfinp)      ! -> nfcatm
!!$    call m_IS_rd_atomic_coordinates(nfinp,nfout,work,ia_cnst_work,fcvect_work)      ! -> nfcatm
    call m_IS_alloc_iatomn_etc()           ! iatomn,iatom,ivan,alfa,amion,zeta1
    call m_IS_alloc_cnstrvectors_etc(imdalg) ! (fcvect,ia_cnst,gca)
    call m_IS_cp_works2fcvect_etc(imdalg,nfcatm,ia_cnst_work,fcvect_work)
    !   (fcvect_work,ia_cnst_work) ->(fcvect,ia_cnst)

    call m_IS_set_iatom(nfout)             ! -> iatom
    call read_iatomn_alfa_amion_etc(nfinp) ! iatomn,alfa,amion
    call m_IS_init_cnstrnt(natm,fcvect_work)
    call cnstr_fcvect_work_dealloc         !-(InputData_Analysis) (*)

    call skip_lines(nfinp,nlines+9)  ! nbztyp
    call read_nbztyp(nfinp)
!!$    call skip_lines(nfinp,2)  ! neg
!!$    call read_neg(nfinp)
!        =================
  end subroutine get_CS_and_ionic_system_data

  subroutine read_ntyp_natm_natm2(nfinp)
    integer, intent(in) :: nfinp

    real(kind=DP) :: wk1, wk2
    read(nfinp,*) wk1, wk2, ntyp, natm, natm2
    if(printable) write(nfout,100) wk1, wk2, ntyp, natm, natm2
100 format(' ',2f8.4,3i4,'   : gmax, gmaxp, ntyp, natm, natm2')
  end subroutine read_ntyp_natm_natm2

  subroutine read_altv(nfinp)
    integer, intent(in) :: nfinp
    integer      :: i,j
    do i = 1, 3
       read(nfinp,*) altv(1,i), altv(2,i), altv(3,i)
       if(printable) write(nfout,'(2H  ,3f20.10)') altv(1,i), altv(2,i), altv(3,i)
    enddo

    if(printable) write(nfout,'("      <<<RECIPROCAL LATTICE VECTOR>>>")')
    call altv_2_rltv(altv,rltv,univol,rvol)
    if(printable) write(nfout,100) ((rltv(i,j),i=1,3),j=1,3)
100 format((' ',3f18.10))
    if(printable) write(nfout,'("VOLUME OF A UNIT CELL=",2f20.10)') univol,rvol

  end subroutine read_altv

  subroutine read_ncord_ninv(nfinp)
    integer, intent(in) :: nfinp

    read(nfinp,*) input_coordinate_system, inversion_symmetry
    if(printable) write(nfout,120) input_coordinate_system, inversion_symmetry
120 format(' ',2i4,'   : ncord, ninv,   : iwei, imdtyp, ityp')
  end subroutine read_ncord_ninv

  subroutine read_atomic_coordinates(nfinp)
    integer, intent(in)  :: nfinp

    integer          :: i, imdt, incunt
    real(kind=DP)    :: v

    if(printable) write(nfout,'(" !! << read_atomic_coordinates >>")')
    nfcatm = 0
    do i = 1, natm
       read(nfinp,*) pos(i,1),pos(i,2),pos(i,3),iwei(i),imdtyp(i),ityp(i)
       if(printable) write(nfout,210) pos(i,1),pos(i,2),pos(i,3),iwei(i),imdtyp(i),ityp(i)
       imdt = imdtyp(i)
       if((imdt<HEAT_BATH .and. imdt/=FIX .and. imdt/=RELAX .and.imdt/=BONDLENGTH_FIX) &
            & .or.&
            & (imdt == HEAT_BATH+BONDLENGTH_FIX_1 .or. imdt == HEAT_BATH+BONDLENGTH_FIX_2 &
            &  .or. imdt == HEAT_BATH+COG_FIX_L)) then
          nfcatm = nfcatm + 1
          if(printable) write(nfout,*)' nfcatm = ', nfcatm
          ia_cnst_work(nfcatm) = i
          if(nfcatm == 1) then
             if(imdt > HEAT_BATH) then
                cnst_typ = imdt - HEAT_BATH
             else
                cnst_typ = imdt
             end if
          end if
          backspace nfinp
          read(nfinp,'(a132)') str
          call chnnm(str,len_str,10,work,incunt)
          if(incunt == 7) then
             fcvect_work(nfcatm,1) = work(7)
          else if(incunt >= 9) then
             if(input_coordinate_system == PUCV) then
                work(11:13) = matmul(altv,work(7:9))
             else if(input_coordinate_system == CARTS) then
                work(11:13) = work(7:9)
             end if
    !        ---> Normalization of directional vectors fcvect(*,1:3).
             v = dsqrt(dot_product(work(11:13),work(11:13)))
             fcvect_work(nfcatm,1:3) = work(11:13)/v

             if(incunt == 10) fcvect_work(nfcatm,4) = work(10)
             if(printable) write(nfout,9001) nfcatm &
                  &     ,fcvect_work(nfcatm,1),fcvect_work(nfcatm,2) &
                  &     ,fcvect_work(nfcatm,3),fcvect_work(nfcatm,4)
9001         format(' Fcvect(',i4,') = ',3f8.4, d20.8)
          endif
       endif
    enddo
210 format(' ',3f18.10,i3,i6,i3)

    if(input_coordinate_system == PUCV) then
       call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_C.S.)
    else if(input_coordinate_system == CARTS) then
       cps = pos
       rltv = transpose(rltv)/PAI2
       call change_of_coordinate_system(rltv,cps,natm,natm,pos) !-(b_C.S.)
       call altv_2_rltv(altv,rltv,univol,rvol) !-(b_Crystal_Structure)
       !                   (altv)->(rltv,univol,rvol)
    end if

    if(printable) then
       write(nfout,*) ' ! nfcatm = ', nfcatm
       write(nfout,*) ' ! cnst_typ = ', cnst_typ
    end if

!!$    call m_IS_init_cnstrnt(natm,fcvect_work)   ! ->sgmc

    if(printable) write(nfout,*) ' ! univol = ', univol
  end subroutine read_atomic_coordinates

  subroutine read_iatomn_alfa_amion_etc(nfinp)
    integer, intent(in) :: nfinp
    integer nsize
    integer i
    do i = 1, ntyp
       read(nfinp,'(a132)') str
       call chnnm(str,len_str,NWK,work,nsize)
       iatomn(i) = work(1)
       alfa(i)   = work(2)
       amion(i)  = work(3)
!!$       iloc_inputf(i)   = nint(work(4))
       ivan(i)   = nint(work(5))
       zeta1(i)  = work(6)
       if(printable) then
          if(nsize >=6 .and. dabs(zeta1(i)) > 1.d-20) then
             write(nfout,511) iatomn(i),alfa(i),amion(i),nint(work(4)),ivan(i),zeta1(i),i
          else
             write(nfout,510) iatomn(i),alfa(i),amion(i),nint(work(4)),ivan(i),i
          endif
       end if
    enddo
510 format(' ',2f8.4,f10.2, 2i2 &
         & ,' :type',i2,'iatomn,alfa,amion,iloc,ivan')
511 format(' ',2f8.4,f10.2, 2i2, f10.2 &
         & ,' :type',i2, 'iatomn,alfa,amion,iloc,ivan,zeta1')
  end subroutine read_iatomn_alfa_amion_etc

  subroutine read_nbztyp(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) nbztyp
    if(printable) write(6,9000) nbztyp
9000 format(' ',i3,'  : nbztyp 0-sf, 1-bk, 2-sc, 3-bcc, 4-fcc, ', &
          &                      '5-dia, 6-hex')
    call m_CS_set_nbztyp_spg() ! -> nbztyp_spg <- nbztyp
!!$    nbztyp1 = nbztyp
    if(nbztyp /= 0 .and. nbztyp /= 1 .and. nbztyp /= GENERAL .and. nbztyp /= GENERAL_LARGER) then
       call m_CS_set_nbztyp(GENERAL) ! nbztyp <- GENERAL
!!$        nbztyp = GENERAL
    end if
  end subroutine read_nbztyp

#ifndef SX
  subroutine parse_args()
    integer :: i,narg,ieq
    character(len=32) :: arg,argval
    narg = command_argument_count()
    do i=1,narg
       call get_command_argument(i,arg)
       argval = arg
       ieq = index(arg,'=')
       if (ieq>0) then
          if(.not.(ieq.eq.len_trim(arg)-1)) then
          call phase_error_with_msg(nfout, 'invalid argument : '//trim(arg),&
         & __LINE__,__FILE__)
          endif
          argval = arg(1:ieq-1)
       endif
       select case (argval)
       case('-p','--preparation')
          icond = PREPARATION_ONLY
       end select
    enddo
  end subroutine parse_args
#endif
end subroutine InputData_Analysis

subroutine InputData_Analysis_neb()
! $Id: InputData_Analysis.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Files, only :             nfout, nfcntn &
       &                        , m_Files_reopen_nfcntn
  use m_Ionic_System, only :      m_IS_rd_pos_and_v
  use m_Control_Parameters, only: icond,nspin &
       &                        , iconvergence_previous_job &
#ifndef _EMPIRICAL_
       &                        , m_CtrlP_rd_edelta_ontheway &
#endif
       &                        , m_CtrlP_rd_iconvergence
  use m_Const_Parameters, only :  CONTINUATION, FIXED_CHARGE_CONTINUATION
  use m_IterationNumbers, only :  m_Iter_rd_iteration_numbers,m_Iter_rd_iters_and_nk

! ============================= added by K. Tagami ==================== 11.0
  use m_Control_Parameters,    only : noncol, ndim_spinor
! ===================================================================== 11.0

  implicit none

  if(icond == CONTINUATION) then
     call m_Files_reopen_nfcntn()
     call m_Iter_rd_iteration_numbers(nfcntn,icond)
!!!!!     call m_IS_rd_pos_and_v(nfcntn)
     call m_CtrlP_rd_iconvergence(nfcntn)
#ifndef _EMPIRICAL_
     call m_CtrlP_rd_edelta_ontheway(nfcntn)
#endif
  else if(icond == FIXED_CHARGE_CONTINUATION) then
     call m_Files_reopen_nfcntn()
     call m_CtrlP_rd_iconvergence(nfcntn)
!!$     call m_Iter_rd_iters_and_nk(iconvergence_previous_job,nfcntn,nspin)

! ========================================= modified by K. Tagami ========= 11.0
!     call m_Iter_rd_iters_and_nk(nfcntn,nspin)

     if ( noncol ) then
       call m_Iter_rd_iters_and_nk(nfcntn,ndim_spinor,nfout)
     else
       call m_Iter_rd_iters_and_nk(nfcntn,nspin,nfout)
     endif
! ======================================================================== 11.0

  end if

end subroutine InputData_Analysis_neb
