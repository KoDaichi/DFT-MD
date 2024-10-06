!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE:  fft_box_finding_way, Preparation, Preparation_ek
!
!  AUTHOR(S): T. Yamasaki and H. Mizouchi   August/20/2003
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
#ifdef JRCATFFT_WS
#define _INCLUDE_EXX_
#elif FFTW3
#define _INCLUDE_EXX_
#endif

#ifdef FJ_TIMER
#   define __TIMER_FJ_START_w_BARRIER(str,a)   call mpi_barrier(str,ierr) ;   call timer_sta(a)
#   define __TIMER_FJ_START(a)   call timer_sta(a)
#   define __TIMER_FJ_STOP(a)    call timer_end(a)
#else
#   define __TIMER_FJ_START(a)
#   define __TIMER_FJ_STOP(a)
#endif

subroutine Preparation(imode)
! $Id: Preparation.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters,  only:DP,FILE,GENERAL,OUTER,INNER,SIMPLE_CUBIC &
       &                       ,HEXAGONAL,TETRAHEDRON &
       &                       ,INITIAL,CONTINUATION,FIXED_CHARGE &
       &                       ,FIXED_CHARGE_CONTINUATION, ON, OFF, GRID &
       &                       ,AUTOMATIC, MANUAL, PREPARATION_ONLY, LDA &
       &                       ,ALL_AT_ONCE, ONE_BY_ONE, DRIVER_NEB &
       &                       ,COORDINATE_CONTINUATION,by_matrix_diagon, MESH &
       &                       ,P_CONTROL, PT_CONTROL &
       &                       ,DRIVER_URAMP,DRIVER_SC_DFT,RETURN_AFTER_SYMMCHECK,SKIP_SYMMCHECK, DRIVER_DIMER
  use m_Parallelization,   only:MPI_CommGroup &
       &                       ,m_Parallel_init_mpi_kngp &
       &                       ,m_Parallel_init_mpi_cdfft &
       &                       ,m_Parallel_end_mpi,   ierr &
       &                       ,ista_kngp,iend_kngp,mype &
       &                       ,m_Parallel_resolve_decomp
  use m_Control_Parameters,only:paramset, ipriparallel, ipri, icond, ekmode &
       &                       ,gmaxs_given,n_matrix_size, sw_positron &
       &                       ,sw_ldos, fixed_charge_k_parallel &
!!$       &                       ,m_CtrlP_set_wct_start &
       &                       ,m_CtrlP_set_kimg &
       &                       ,m_CtrlP_way_of_smearing &
       &                       ,sw_pdos, sw_berry_phase, sw_phonon &
       &                       ,num_projectors, ggacmp_parallel, sw_fef &
       &                       ,sw_wannier, numk_tmp, iconv_ek_tmp, driver &
       &                       ,m_CtrlP_rd_neg_previous &
       &                       ,m_CtrlP_rd_edelta_ontheway &
       &                       ,m_CtrlP_rd_iconvergence &
       &                       ,m_CtrlP_rd_corecharge_cntnbin &
       &                       ,sw_hybrid_functional, m_CntrlP_rst_submat_call_stat &
       &                       ,sw_optimize_lattice,nhistory_stress &
       &                       ,sw_optimize_coordinates_once &
       &                       ,m_CtrlP_reset_iconvergence &
       &                       ,gmaxp, sw_rebuild_pws, way_of_smearing &
       &                       ,m_CtrlP_in_initialization &
       &                       ,m_CtrlP_reset_optmode &
       &                       ,m_CtrlP_set_npartition_david &
       &                       ,intzaj, imdalg &
       &                       ,sw_communicator_for_chg &
       &                       ,gmax,gmax_org, printable, nspin, m_CtrlP_set_npartition_david &
       &                       ,sw_charge_rspace, sw_partial_charge
#ifdef _USE_SCALAPACK_
  use m_Control_Parameters,only: m_CtrlP_set_sw_scalapack
#endif
  use m_Files,             only:nfinp,nfout,nfkpgn,nfcntn &
       &                       ,m_Files_open_kpoint_files &
       &                       ,m_Files_check_file_names &
       &                       ,m_Files_open_nfldos &
       &                       ,m_Files_close_files_initial0 &
       &                       ,m_Files_open_nfeng &
       &                       ,m_Files_open_nfcntn &
       &                       ,m_Files_reopen_nfzaj_kall &
       &                       ,nfdynm
  use m_FFT,               only:m_FFT_set_box_sizes, m_FFT_setup,fft_box_size_CD &
       &                       ,m_FFT_set_box_size_prev, m_FFT_set_box_size_cd_prev
  use m_PlaneWaveBasisSet, only:n_rGv,n_rGpv,n_rGpv_reduced, n_rGv_pstrn,kgp &
       &                       ,n_rGv_prev, n_rGpv_prev &
       &                       ,m_pwBS_decide_cutoff_mix &
       &                       ,m_pwBS_assume_G_rhombohedron &
       &                       ,m_pwBS_for_each_WF &
       &                       ,m_pwBS_generate_G_vectors &
       &                       ,m_pwBS_G_trans_functions  &
       &                       ,m_pwBS_alloc_ngpt_igfp_gr &
       &                       ,m_pwBS_calc_length_of_G &
       &                       ,m_pwBS_setup_FFTmapfunctions &
       &                       ,m_pwBS_rebuild_gr_l &
       &                       ,m_pwBS_positronWF &
       &                       ,m_pwBS_cp_iba_to_iba_ek &
       &                       ,m_pwBS_rebuild_gr_l, ngabc
  use m_Crystal_Structure, only:nbztyp, inversion_symmetry &
       &                       ,nbztyp_spg &
       &                       ,symmetry_method &
       &                       ,univol,altv,rltv, &
       &                           check_if_sw_inversion_is_valid &
       &                       ,m_CS_gnrt_symmetry_operations &
       &                       ,m_CS_gnrt_symm_operators_tl &
       &                       ,m_CS_alloc_op_tau &
       &                       ,m_CS_alloc_op_tau_tl &
       &                       ,symmetry_method &
       &                       ,m_CS_set_altv_prim, m_CS_set_altv_super &
       &                       ,m_CS_supercell &
       &                       ,m_CS_rd_fix_spin_status,altv &
       &                       ,m_CS_wd_op_and_tau,tau

  use m_CS_SpaceGroup,     only:m_CS_SG_auto_gnrt_sym_op &
       &                       ,m_CS_SG_print_space_group_name &
       &                       ,m_CS_remove_frac_translation
  use m_Ionic_System,      only:m_IS_alloc_napt,m_IS_symm_check_of_pos &
       &                       ,m_IS_phonon_initial_disp &
       &                       ,m_IS_phonon_init_firsthalf &
       &                       ,m_IS_phonon_init_secondhalf &
       &                       ,m_IS_phonon_set_displacement &
       &                       ,m_IS_set_napt_super, m_IS_supercell &
       &                       ,m_IS_inv_sym_off &
       &                       ,m_IS_symmetrize_atom_pos &
       &                       ,m_IS_rd_pos_and_v &
       &                       ,m_IS_natm_can_change &
       &                       ,m_IS_change_natm &
       &                       ,m_IS_wd_speciesname_etc &
       &                       ,m_IS_get_neg_incre,natm,cpd_l &
       &                       ,m_IS_gdiis_reset &
       &                       ,m_IS_CG_reset &
       &                       ,m_IS_freeze &
       &                       ,m_IS_gnrt_supercell_symmetry &
       &                       ,natm,ityp,iatomn,ntyp &
       &                       ,natm,ityp,iatomn,ntyp &
       &                       ,m_IS_reset_extrpl_status
  use m_PseudoPotential,   only:m_PP_input_xctype,ival,m_PP_renew_etot1
  use m_Kpoints,           only:way_ksample,kv3 &
       &                       ,m_Kp_gnrt_or_rd_k_points &
       &                       ,m_Kp_alloc_kpoints &
       &                       ,m_Kp_cr_kpoints_table &
       &                       ,m_Kp_alloc_kpoints_ek &
       &                       ,m_Kp_cp_vkxyz_to_vkxyz_ek &
       &                       ,m_Kp_realloc_kpoints &
       &                       ,m_Kp_set_mesh_super &
       &                       ,m_Kp_set_ek_group, m_Kp_realloc_kpoints2
! === DEBUG for icond=2,3 & one_by_one by tkato 2014/01/24 =====================
  use m_Kpoints,           only:m_Kp_cp_vkxyz_ek_to_vkxyz, m_Kp_set_star_of_k, &
       &                        sw_force_kpt_inside_bz
! ==============================================================================

  use m_Ldos,              only:m_Ldos_preparation
  use m_BerryPhase,        only:m_BP_gen_Kpoints
  use m_Orbital_Population,only: m_OP_alloc
  use m_Phonon,            only:m_Phonon_alloc_qvec, m_Phonon_set_qvec, m_Phonon_read_forces, sw_read_forces_pre
  use m_Electronic_Structure,only:m_ES_set_num_bands_super &
       &                       ,m_ES_cp_iconv, m_ES_add_neg
  use m_Dipole,            only:m_Dipole_calc
  use m_Wannier,           only:m_Wan_gen_weight
  use m_FiniteElectricField, only:m_FEF_init
  use m_IterationNumbers, only : m_Iter_rd_iteration_numbers,iteration_unit_cell,m_Iter_cmix_reset&
       &                       ,iteration_uramp,iteration_scdft, iteration_stress_correction, iteration_dimer
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
  use m_Control_Parameters, only: m_CtrlP_rd_dtim_previous
! ==============================================================================

  use m_UnitCell, only : m_UC_doit,m_UC_get_univol_old,m_UnitCell_md
  use m_XC_Potential, only : m_XC_rst_npsend
  use m_FFT, only : m_FFT_reset_CD_setup_stat
  use m_CD_mixing, only : m_CD_force_dealloc
  use m_ES_WF_by_RMM, only : m_ESrmm_dealloc_r_norm_flag
  use m_ES_WF_by_submat, only : m_ESsubmat_dealloc

! ================================= added by K. Tagami =================== 11.0&13.0U
  use m_Control_Parameters,   only : noncol
  use m_Ionic_System,         only : m_IS_alloc_magmom_local, m_IS_init_magmom_local
  use m_CS_Magnetic,     only : m_CS_set_Magnetic_Sym, m_CS_set_inverse_operation, &
       &                        m_CS_chk_determinant_op
  use m_Crystal_Structure,   only : sw_use_magnetic_symmetry, sw_allow_frac_translation

  use m_CD_Mag_Moment,      only : m_CD_alloc_rad_cov, m_CD_set_rad_cov_default, &
       &                           m_CD_set_rad_cov_now, sw_monitor_atomcharge, &
       &                           m_CD_set_sw_monitor_atomcharge, &
       &                           m_CD_set_proj_radius
! ======================================================================== 11.0&13.0U


  use m_Charge_Density,       only : m_CD_cp_chgq_to_chgqo

  use m_UnitCell,             only : m_UC_rd_cntn_data

#ifdef _INCLUDE_EXX_
! ============== KT_add ======================================= 13.0F
  use m_Control_Parameters,  only : use_fft_exx
  use m_PlaneWaveBasisSet, only : n_rGv_exx, n_rGpv_exx, m_pwBS_exxWF, m_pwBS_exxCD, m_pwBS_store_prev_kg1_kgp
  use m_FFT, only : m_FFT_set_box_size_exx, m_FFT_set_box_size_cd_exx, m_FFT_setup_exx
! ============================================================= 13.0F
#endif

! ==== KT_add ==== 13.1R
  use m_Control_Parameters,  only : sw_keep_symmetry_strict, sw_optimize_coords_sametime
  use m_Ionic_System, only : m_IS_symmetrize_atom_pos
  use m_Crystal_Structure, only : m_CS_read_op_tau_previous, m_CS_dealloc_op_tau
! ================ 13.1R
! ===== KT_add === 13.1XI
  use m_Control_Parameters,  only : sw_phonon_with_epsilon, sw_excitation
  use m_Excitation,   only : m_XI_set_value_nmax_G
! ================ 13.1XI

  use m_ES_WF_by_MatDiagon, only : m_ESmat_set_reduced_basis_mode

! ==== EXP_CELLOPT === 2015/09/24
  use m_PlaneWaveBasisSet, only : m_pwBS_store_prev_kg1_kgp
! ==================== 2015/09/24
  use m_Control_Parameters,   only : sw_band_unfolding, howto_set_proj_radius, sw_lband
  use m_Band_Unfolding,       only : m_BU_set_GVec_flag_refcell
  use m_Files,                only : m_Files_open_nfwfk_lband

  use m_Force,                only : forc_l

  use m_Control_Parameters,   only : sw_stress_correction
  use m_Stress,               only : m_Stress_in_correction, m_Stress_decut_for_correction &
   &                               , m_Stress_correction, m_Stress_rd_cdata_4_correction, m_Stress_correction
  use mpi

  implicit none
  integer, intent(in) :: imode
  integer :: imd
  integer :: outer_or_inner
!  include 'mpif.h'
  integer :: ggacmp_parallel_rev, xctype_is
  integer :: ne,i

  real(kind=DP), dimension(3,3) :: stress_tensor
  real(kind=DP) :: gmax_factor
  logical :: ini,initialization_required
  logical :: logi

  ini = initialization_required()
  if(ini) then
  call m_Files_check_file_names()
!!$  call m_CtrlP_set_wct_start       ! (ckcput)

  if(ekmode /= GRID .and. sw_phonon == OFF) then
     ! setup supercell
     call m_CS_supercell(nfout)
     call m_IS_supercell(nfout)
  end if
  if (m_IS_natm_can_change())then
     if(m_IS_change_natm()) then
        call m_IS_wd_speciesname_etc(nfdynm)
        ne = m_IS_get_neg_incre()
        call m_ES_add_neg(ne)
     endif
  endif
  endif

  if(icond==CONTINUATION .and.  m_Stress_in_correction())then
     call m_Files_open_nfcntn()
     call m_Stress_rd_cdata_4_correction(nfout,nfcntn)
     call m_Stress_correction(nfout,.true.)
  endif

  if(m_Stress_in_correction(3)) then
     gmax_factor = gmaxp/gmax
     gmax = dsqrt((gmax_org*gmax_org) + m_Stress_decut_for_correction())
     gmaxp = gmax * gmax_factor
     if(printable) then
        write(nfout,'(a,f10.5)') ' !** new cutoff : ',gmax_org*gmax_org+m_Stress_decut_for_correction()
        write(nfout,'(a,2f10.5)') ' !** new gmax and gmaxp : ',gmax,gmaxp
     endif
  endif

  if ((((sw_optimize_lattice==ON .or. imdalg == P_CONTROL .or. imdalg == PT_CONTROL) &
       & .and. iteration_unit_cell>1) .or.  &
       & ((driver == DRIVER_URAMP)  .and. (iteration_uramp>1)) .or.  &
       & ((driver == DRIVER_SC_DFT) .and. (iteration_scdft>1)) .or.  &
       & ((iteration_stress_correction<=3) .and. (iteration_stress_correction>1)))  &
       & .and. &
       & .not.(icond==CONTINUATION.and.m_CtrlP_in_initialization()))then
     if(ini)then
        call m_XC_rst_npsend()
        call m_FFT_reset_CD_setup_stat()
        call m_CD_force_dealloc()
        call m_ESrmm_dealloc_r_norm_flag()
        call m_ESsubmat_dealloc()
        call m_pwBS_store_prev_kg1_kgp()
     endif
     call m_IS_gdiis_reset()
     call m_IS_CG_reset()
     if ( sw_optimize_coords_sametime == OFF ) then
        call m_IS_reset_extrpl_status()
     endif
     call m_CtrlP_reset_optmode()
     if(sw_optimize_lattice==ON .and. sw_optimize_coordinates_once==ON)then
        call m_IS_freeze()
     endif
     if(sw_optimize_lattice == ON .and. .not. m_Stress_in_correction(3)) call m_UC_doit()
     if(imdalg == PT_CONTROL .or. imdalg == P_CONTROL) call m_UnitCell_md(forc_l)
     if(driver /= DRIVER_DIMER) call m_IS_wd_speciesname_etc(nfdynm)
  endif
!  call m_CntrlP_rst_submat_call_stat()

  if(intzaj==by_matrix_diagon) call m_ESmat_set_reduced_basis_mode(.true.)
!  if( way_ksample == MESH) call rd_spg_tetra_vars()

  if(.not.ini) then
     call m_pwBS_rebuild_gr_l()
     call m_Kp_gnrt_or_rd_k_points(nfinp,preallocation=.true.) !(kstep) -> kv3
!!$     if( m_CtrlP_way_of_smearing() == TETRAHEDRON ) call m_Kp_cr_kpoints_table()
     if( way_ksample == MESH ) call m_Kp_cr_kpoints_table()
     if(ipri>=1) write(nfout,'(" !** smearing_method = ",i6)') way_of_smearing

     if(sw_berry_phase == OFF) then
        call m_Kp_alloc_kpoints    ! <- kv3  -> vkxyz, qwgt
     else
        call m_BP_gen_Kpoints(preallocation=.true.) ! -> kv3
     end if
     call m_Kp_alloc_kpoints    ! with using the given value of kv3, allocate vkxyz and qwgt
     if(sw_berry_phase == ON) then
        call m_BP_gen_Kpoints(preallocation=.false.) ! -> vkxyz, qwgt
     end if
     call m_Kp_gnrt_or_rd_k_points(nfinp,preallocation=.false.) !(kstep) -> kv3
     if(imode /= SKIP_SYMMCHECK)then
       if(symmetry_method == AUTOMATIC) then
          call m_CS_SG_auto_gnrt_sym_op(paramset,nfout) ! paramset == .false.
          call m_IS_symmetrize_atom_pos(nfout) ! -> cps,pos
       else
          call m_CS_gnrt_symmetry_operations(paramset,nfout) ! paramset == .false.
          call m_CS_gnrt_symm_operators_tl(paramset,nfout) ! -(m_Crystal_Structure) -> op_tl, tau_tl
       end if

       if ( sw_phonon == OFF ) call m_IS_gnrt_supercell_symmetry(paramset,nfout)
     endif

     return
  endif

! === KT_add ==== 13.1R
  if (sw_optimize_lattice==ON.and.iteration_unit_cell>1.and. &
  &    .not.(icond==CONTINUATION.and.m_CtrlP_in_initialization()))then
     if (  sw_keep_symmetry_strict == ON ) then
        call m_CS_read_op_tau_previous
        call m_IS_symmetrize_atom_pos(nfout)
        call m_CS_dealloc_op_tau
     endif
  endif
! =============== 13.1R

! ============================== added by K. Tagami ===================== 11.0&13.0U
  if ( noncol .or. sw_use_magnetic_symmetry == ON ) then
     call m_IS_alloc_magmom_local
     call m_IS_init_magmom_local
  endif
!
  call m_CD_set_sw_monitor_atomcharge
  if ( sw_monitor_atomcharge == ON ) then
     call m_CD_alloc_rad_cov
     call m_CD_set_rad_cov_default
     call m_CD_set_rad_cov_now
  endif
  if ( howto_set_proj_radius == 2 ) then
     call m_CD_set_rad_cov_default
     call m_CD_set_proj_radius
  endif
! ======================================================================== 11.0&13.0U

  if(imode /= SKIP_SYMMCHECK)then
    if(symmetry_method == AUTOMATIC) then
       call m_CS_SG_auto_gnrt_sym_op(.true.,nfout) ! -(m_CS_SpaceGroup) -> nopr,af
    else
                                                     __TIMER_FJ_START(31)
     call m_CS_gnrt_symmetry_operations(.true.,nfout) ! -(m_Crystal_Structure) -> nopr,af
!!$     call m_CS_gnrt_symm_operators_tl(.true.,nfout) ! -(m_Crystal_Structure) -> nopr,af
                                                     __TIMER_FJ_STOP(31)
  end if
  call m_CS_alloc_op_tau(nfout)
  call m_CS_alloc_op_tau_tl(nfout)
  if(symmetry_method == AUTOMATIC) then
     call m_CS_SG_auto_gnrt_sym_op(paramset,nfout) ! paramset == .false.
!     if(.not. m_Stress_in_correction()) call m_IS_symmetrize_atom_pos(nfout) ! -> cps,pos
     call m_IS_symmetrize_atom_pos(nfout) ! -> cps,pos
  else
                                                     __TIMER_FJ_START(31)
     call m_CS_gnrt_symmetry_operations(paramset,nfout) ! paramset == .false.
     call m_CS_gnrt_symm_operators_tl(paramset,nfout) ! -(m_Crystal_Structure) -> nopr,af
                                                     __TIMER_FJ_STOP(31)
  end if

  call m_CS_SG_print_space_group_name(nfout)
  endif

! ============================== added by K. Tagami ===================== 11.0Ex
!  if ( noncol .and. symmetry_method /= AUTOMATIC &
!       &      .and. sw_use_magnetic_symmetry == ON ) then
!     call m_CS_set_Magnetic_Sym
!  endif
  if ( noncol ) then
     if ( sw_use_magnetic_symmetry == ON ) then
        call m_CS_set_Magnetic_Sym
     endif
  else
     call m_CS_set_Magnetic_Sym
  endif
  call m_CS_set_inverse_operation
  call m_CS_chk_determinant_op
  if ( sw_allow_frac_translation == OFF ) call m_CS_remove_frac_translation
! ======================================================================== 11.0Ex

! ============================== added by K. Tagami ===================== 12.0YAM
!  commented by T. Yamasaki 2012/12/06
!!$  call check_if_sw_inversion_is_valid( nfout )
! ======================================================================= 12.0YAM

  if(ekmode /= GRID .and. sw_phonon == ON) then
     call m_IS_inv_sym_off(nfout) ! -> inversion_symmetry
  end if
  call m_IS_alloc_napt()
                                                     __TIMER_FJ_START(31)
  if(imode /= SKIP_SYMMCHECK)then
    call m_CS_wd_op_and_tau(nfout)
  endif
  call m_IS_symm_check_of_pos()
                                                     __TIMER_FJ_STOP(31)
  if(ekmode /= GRID .and. sw_phonon == ON) then
     call m_Phonon_alloc_qvec()
     call m_Phonon_set_qvec(nfout)
     call m_CS_supercell(nfout)
     call m_IS_supercell(nfout)
     call m_IS_phonon_init_firsthalf()
     if(sw_read_forces_pre==ON) call m_Phonon_read_forces()
     call m_IS_phonon_init_secondhalf()
     call m_IS_set_napt_super
     call m_IS_phonon_initial_disp()
     call m_ES_set_num_bands_super
     call m_Kp_set_mesh_super
     if(driver .ne. DRIVER_NEB .and. icond == CONTINUATION) then
        call m_Files_open_nfcntn()
        call m_Iter_rd_iteration_numbers(nfcntn,icond)
	call m_IS_phonon_set_displacement()  ! asms
        call m_IS_rd_pos_and_v(nfcntn)
#ifndef _EMPIRICAL_
        call m_CtrlP_rd_neg_previous(nfcntn,nfout)
        call m_CtrlP_rd_edelta_ontheway(nfcntn)
#endif
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
        call m_CtrlP_rd_dtim_previous(nfcntn)
! ==============================================================================
        call m_CtrlP_rd_iconvergence(nfcntn)
        call m_CtrlP_rd_corecharge_cntnbin(nfcntn) ! -> status_cntnbin_positron
        call m_CS_rd_fix_spin_status(nfcntn,nfout)  ! <-- T. Yamasaki, 18th Aug. 2009
     end if
  end if

  call m_CtrlP_set_kimg(inversion_symmetry,nfout) ! ->kimg

  call m_pwBS_assume_G_rhombohedron() ! ->n_rGv,n_rGpv, n_rGv_pstrn
  call fft_box_finding_way(outer_or_inner)     ! -(contained here)
!!$  call m_FFT_set_box_sizes(n_rGv,n_rGpv,n_rGv_pstrn,outer_or_inner) ! ->fft_box_size_WF,CD
  call m_FFT_set_box_sizes(n_rGv,n_rGpv_reduced,n_rGv_pstrn,outer_or_inner) ! ->fft_box_size_WF,CD

#ifdef _INCLUDE_EXX_
! ================ KT_add ==================================== 13.0F
  if ( sw_hybrid_functional == ON .and. use_fft_exx ) then
     call m_FFT_set_box_size_exx( n_rGv_exx,outer_or_inner) ! ->fft_box_size_exx
  endif
  if(sw_hybrid_functional == ON ) call m_FFT_set_box_size_cd_exx( n_rGpv_exx, outer_or_inner)
! ============================================================ 13.0F
#endif

  if(ekmode /= GRID) call m_Kp_gnrt_or_rd_k_points(nfinp,preallocation=.true.) !(kstep) -> kv3
  if(sw_berry_phase == OFF) then
    call m_Parallel_resolve_decomp(nfout,kv3/nspin,logi)
    if (logi) then
#ifdef _USE_SCALAPACK_
      call m_CtrlP_set_sw_scalapack(printable)
#endif
      call m_CtrlP_set_npartition_david(nfout)
    endif
  endif
  call m_pwBS_generate_G_vectors()    ! ->n_rGv,n_rGpv ->kgp
  call m_Parallel_init_mpi_kngp(nfout,ipriparallel,kgp)  ! -(m_Parallelization) ->ista_kngp,iend_kngp

!!$  if(ipri>=1) write(nfout,'(" nbztyp (nbztyp_spg) = ", i3)') nbztyp_spg
!  if(ekmode /= GRID .and. (nbztyp_spg >= GENERAL .or. way_ksample == FILE)) &
!       & call m_Files_open_kpoint_files(way_ksample,nbztyp_spg)

  call m_pwBS_alloc_ngpt_igfp_gr()
  call m_pwBS_calc_length_of_G()         ! -> gr_l
  call m_pwBS_G_trans_functions()   ! -> ngpt_l: Set of G-vectors translated according to symmetry operations

  ggacmp_parallel_rev = ggacmp_parallel
  call m_PP_input_xctype(xctype_is) ; if(xctype_is == LDA) ggacmp_parallel_rev = OFF
!!$  if(ipri >= 1) write(nfout,'(" ggacmp_parallel_rev = ",i3)') ggacmp_parallel_rev
  call m_Parallel_init_mpi_cdfft(nfout,ipriparallel,ggacmp_parallel_rev)
                          ! -> nrank_ggacmp, npes_cdfft, and nrest_cdfft
  call m_FFT_setup(inversion_symmetry,paramset) ! paramset == .false.

#ifdef _INCLUDE_EXX_
#ifdef FFTW3
! ================ KT_add =================================== 13.0F
  if ( sw_hybrid_functional == ON .and. use_fft_exx ) then
     call m_FFT_setup_exx( inversion_symmetry,paramset ) ! paramset == .false.
  endif
! =========================================================== 13.0F
#endif
#endif


  call m_pwBS_setup_FFTmapfunctions()
  if(num_projectors>0) call m_OP_alloc()

  if(ekmode /= GRID) then
!!$     if( m_CtrlP_way_of_smearing() == TETRAHEDRON ) call m_Kp_cr_kpoints_table()
     if( way_ksample == MESH ) call m_Kp_cr_kpoints_table()
     if(sw_berry_phase == OFF) then
        call m_Kp_alloc_kpoints    ! <- kv3  -> vkxyz, qwgt
     else
        call m_BP_gen_Kpoints(preallocation=.true.) ! -> kv3
     end if
     call m_Kp_alloc_kpoints    ! with using the given value of kv3, allocate vkxyz and qwgt
     if(sw_berry_phase == ON) then
        call m_BP_gen_Kpoints(preallocation=.false.) ! -> vkxyz, qwgt
     end if
     call m_Kp_gnrt_or_rd_k_points(nfinp,preallocation=.false.) ! -> vkxyz, qwgt (,tk_initial)

     if(icond == PREPARATION_ONLY .or. icond == INITIAL .or. icond == CONTINUATION .or. &
        & icond==COORDINATE_CONTINUATION) then
                                                     __TIMER_FJ_START(28)
        call m_pwBS_for_each_WF(preallocation=paramset) ! -> kg1, nbase,iba (when paramset==.false.)
                                                     __TIMER_FJ_STOP(28)
     else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
!!$        if(ekmode == OFF)&
!!$             & stop ' ! combination of ekmode and icond is illegal (Preparation)'
        call m_Files_close_files_initial0()
        call m_Files_open_nfeng(icond)
        if(ekmode == OFF .and. fixed_charge_k_parallel == ALL_AT_ONCE) then
                                                     __TIMER_FJ_START(28)
           call m_pwBS_for_each_WF(preallocation=paramset) ! -> kg1, nbase,iba (when paramset==.false.)
                                                     __TIMER_FJ_STOP(28)
           if(star_of_k()) call m_Kp_set_star_of_k
        else
                                                     __TIMER_FJ_START(28)
           call m_pwBS_for_each_WF(preallocation=.true.) ! -> kg1, iba
                                                     __TIMER_FJ_STOP(28)
           call m_Kp_alloc_kpoints_ek()  ! -> kv3_ek (=kv3), allocate(vkxyz_ek,qwgt_ek)
           call m_Kp_cp_vkxyz_to_vkxyz_ek()
           if(star_of_k()) call m_Kp_set_star_of_k

           if(fixed_charge_k_parallel == ONE_BY_ONE) then
              call m_Kp_set_ek_group()
              call m_Kp_realloc_kpoints2()
! === DEBUG for icond=2,3 & one_by_one by tkato 2014/01/24 =====================
!              call m_Kp_cp_vkxyz_ek_to_vkxyz(1)
! ==============================================================================
           else
              call m_Kp_realloc_kpoints()
           end if
           call m_pwBS_cp_iba_to_iba_ek()
           if(icond == FIXED_CHARGE_CONTINUATION) &
                & call m_ES_cp_iconv(numk_tmp,iconv_ek_tmp)
                                                     __TIMER_FJ_START(28)
           call m_pwBS_for_each_WF(preallocation=.false.) ! -> kg1, iba
                                                     __TIMER_FJ_STOP(28)
           if( icond == FIXED_CHARGE .and. &
           &  (sw_ldos==ON .or. sw_pdos==ON .or. sw_partial_charge==ON .or. sw_charge_rspace==ON)) &
           &   call m_Files_reopen_nfzaj_kall()
        end if
     else
        !stop ' icond is illegal (Preparation)'
        call phase_error_with_msg(nfout,'icond is illegal (Preparation) ',__LINE__,__FILE__)
     end if
     call m_pwBS_decide_cutoff_mix()    ! ->kgpm

     if(sw_positron /= OFF) call m_pwBS_positronWF()

#ifdef _INCLUDE_EXX_
! ============= KT_add ============================= 13.0F
     if (sw_hybrid_functional == ON .and. use_fft_exx ) then
        call m_pwBS_exxWF()
        call m_pwBS_exxCD()
     endif
! ========== ======================================= 13.0F
     if(ipri>=1) write(nfout,'(" _INCLUDE_EXX_ is defined")')
#else
     if(ipri>=1) write(nfout,'(" _INCLUDE_EXX_ is not defined")')
#endif

     if(sw_ldos == ON .or. sw_lband==ON ) then
        call m_Ldos_preparation()
        if(sw_ldos == ON) call m_Files_open_nfldos()
     end if

     if(sw_wannier == ON) then
        call m_Wan_gen_weight(nfout)
     end if

     if(sw_fef == ON) call m_FEF_init(nfout)

  end if

!  if( way_ksample == MESH) call wd_spg_tetra_vars()

  if(imode == RETURN_AFTER_SYMMCHECK)then
    return
  endif

!!$  call m_Parallel_resolve_decomp(nfout,kv3/nspin)

  if(icond == PREPARATION_ONLY) then
     call mpi_barrier(MPI_CommGroup,ierr)
     call m_Parallel_end_mpi()
     stop 'The preparation has been done.'
  end if

  if(sw_optimize_lattice==ON.and.icond==CONTINUATION.and.sw_rebuild_pws==OFF) then
     call m_Files_open_nfcntn()
     call m_UC_rd_cntn_data(nfcntn)
  endif

! ===== KT_add ==== 13.1XI
!  if ( sw_phonon_with_epsilon == ON ) then
!     sw_excitation = ON
!     write(nfout,*) "!** XI: sw_excitation is turned on"
!  endif
  if ( sw_excitation == ON ) call m_XI_set_value_nmax_G
! ================= 13.1XI

  if ( sw_phonon == OFF ) call m_IS_gnrt_supercell_symmetry(paramset,nfout)

  if ( sw_band_unfolding == ON ) then
     if ( sw_force_kpt_inside_bz == OFF ) call m_BU_set_GVec_flag_refcell
  endif

  call m_Iter_cmix_reset()

contains
  subroutine fft_box_finding_way(outer_or_inner)
    integer,intent(out)::  outer_or_inner
    if( nbztyp == GENERAL &
         & .or.(nbztyp >= 30 .and. nbztyp <= 32) &
         & .or.(nbztyp >= SIMPLE_CUBIC  .and. nbztyp <= HEXAGONAL)) then
       outer_or_inner = OUTER
    else
       outer_or_inner = INNER
    endif
  end subroutine fft_box_finding_way

  function star_of_k() result(res)
    use m_Control_Parameters, only : charge_symm_mode, sw_excitation, &
  & sw_local_approx_trans_moment, sw_write_procar_file
    use m_Const_Parameters, only : ON, chg_symm_level2
    logical :: res
    res = charge_symm_mode == chg_symm_level2 .or. &
  & sw_excitation == ON .or. sw_local_approx_trans_moment == ON .or. &
  & sw_write_procar_file == ON
  end function star_of_k

end subroutine Preparation

subroutine Preparation_ek()
  use m_IterationNumbers,   only : nk_in_the_process
  use m_PlaneWaveBasisSet,  only : m_pwBS_for_each_WF, m_pwBS_cp_iba_ek_to_iba
  use m_Kpoints,            only : m_Kp_cp_vkxyz_ek_to_vkxyz
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : sw_dos, sw_ldos
!  use m_Ldos,               only : m_Ldos_alloc_weiwsc_etc

  call m_Kp_cp_vkxyz_ek_to_vkxyz(nk_in_the_process) !(kreset)
  call m_pwBS_cp_iba_ek_to_iba(nk_in_the_process)
  call m_pwBS_for_each_WF(preallocation=.false.)     !(basnum) ->(nbase)
!  if(sw_dos == ON .and. sw_ldos == ON) then
!    call m_Ldos_alloc_weiwsc_etc()
!  endif
end subroutine Preparation_ek

subroutine Preparation_grid(nk,kxyz)
  use m_Files, only           : nfout
  use m_Const_Parameters,only :   DP, GRID
  use m_Control_Parameters,only : ekmode
  use m_Kpoints, only :           m_Kp_set_kv3 &
       &                        , m_Kp_alloc_kpoints &
       &                        , m_Kp_realloc_kpoints &
       &                        , m_Kp_cp_kxyz_to_vkxyz
  use m_PlaneWaveBasisSet, only : m_pwBS_for_each_WF &
       &                        , m_pwBS_increase_kg1
  implicit none
  integer, intent(in) :: nk
  real(kind=DP), dimension(3), intent(in) :: kxyz

  if(ekmode /= GRID)&
  & call phase_error_with_msg(nfout, ' combination of ekmode and icond is illegal (Preparation_grid)', __LINE__,__FILE__)

  call m_Kp_set_kv3(nk)
  call m_Kp_alloc_kpoints()
  call m_Kp_cp_kxyz_to_vkxyz(nk,kxyz) ! kxyz -> vkxyz

  call m_pwBS_for_each_WF(preallocation=.true.) ! -> kg1
  call m_pwBS_increase_kg1(30)

  call m_Kp_realloc_kpoints()

end subroutine Preparation_grid

subroutine Preparation_ek_grid2(kxyz)
  use m_Const_Parameters, only : DP
  use m_PlaneWaveBasisSet,only : m_pwBS_for_each_WF
  use m_Kpoints,          only : m_Kp_cp_kxyz_to_vkxyz
  real(kind=DP),intent(in),dimension(3):: kxyz

  call m_Kp_cp_kxyz_to_vkxyz(1,kxyz)
  call m_pwBS_for_each_WF(preallocation=.false.)
end subroutine Preparation_ek_grid2

