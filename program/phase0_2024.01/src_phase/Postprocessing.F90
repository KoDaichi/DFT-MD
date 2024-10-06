!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 1200)
!
!  SUBROUINE: Postprocessing
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
! ====================================
!  patch 0.1 by K. Tagami @adv    2008/08/13
!  patch 0.2 by J. Koga   @adv    2008/08/20
!  patch 0.3 by K. Tagami @adv    2008/08/26
!  patch 0.4 by K. Tagami @adv    2009/04/21
!  patch 0.5 by K. Tagami @adv    2009/05/09
!  patch 0.6 by K. Tagami @adv    2009/06/02
!
!     patch 0.1 :  using new subroutine m_ESIO_wd_ElectroStaticPot
!     patch 0.2 :  correction in wd_Header_For_Cube
!     patch 0.3 :  using new subroutine m_ESIO_wd_HartreePot
!                                       m_ESIO_wd_LocalPPot
!     patch 0.4 :  calculation of orbital_population only in case of ekmode=off
!                ( This is introduced for band calculation with DFT+U )
!     patch 0.5 :  calculation of DOS in case of DFT+U
!                ( ekcal , only trial )
!     patch 0.6 :  treatment in the case of sw_calc_force = ON
!                   ( trial )
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
!$$#ifndef PARA3D
subroutine Postprocessing(ignore_convergence)
! $Id: Postprocessing.F90 634 2020-12-22 12:17:00Z ktagami $
  use m_Const_Parameters, only :   DP, ON, OFF, FORCE_CONVERGED, INITIAL, CONTINUATION &
       &                         , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                         , Gauss_distrib_func, EK, SCF &
       &                         , MESH, TETRAHEDRON, BULK, DEFECT, TOTAL &
       &                         , Valence_plus_PC_Charge &
       &                         , VXC_AND_EXC,CARTS,BUCS &
       &                         , unit_conv_byname, YES, SEPARATE, INTEGRATED &
       &                         , POLARIZATION, EFFECTIVE_CHARGE &
       &                         , PIEZOELECTRIC_CONST, PHONON_GAMMA, PHONON_DOS &
       &                         , PARABOLIC, TETRAHEDRON, COLD, DRIVER_CONSTRAINT, PAI2
  use m_Timing,             only : tstatc_wd, tstatc_init, tstatc_iter, tstatc_wd0
  use m_Control_Parameters, only : nspin,icond,sw_dos,sw_orb_popu,dos_method, ldos_method &
       &                         , iconvergence,iconvergence_previous_job &
       &                         , neg, kimg, ekmode, ipri, ipridos,partial_charge_filetype, ipritiming0 &
       &                         , sw_charge_rspace, sw_positron, sw_elf, sw_ldos, sw_cal_ldos &
       &                         , sw_berry_phase, polar_prop, sw_vibrational_modes &
       &                         , sw_lattice_dielectric_tensor, sw_dielectric_function &
       &                         , sw_lo_to_splitting, phonon_method &
       &                         , sw_int_strain_piezo_tensor &
       &                         , sw_fine_STM_simulation, sw_deficit_charge, sw_add_xc_to_vloc  &
       &                         , sw_partial_charge, sw_phonon_oneshot &
       &                         , m_CntrlP_keep_charge_title &
       &                         , m_CntrlP_retrieve_charge_title &
       &                         , m_CntrlP_set_pcharge_title &
!!$       &                         , m_CtrlP_way_of_smearing &
       &                         , m_CntrlP_set_crtdst &
       &                         , sw_dipole,sw_layered,sw_dipole_correction &
       &                         , sw_wf_rspace, sw_rwf2, nr, rmax, center &
       &                         , ik_rwf2, ib_rwf2, sw_wannier, sw_wannier90 &
       &                         , printable, sw_xc_only, driver,sw_rsb_test &
       &                         , sw_spin_magmom_rspace, sw_write_unk_file &
       &                         , sw_diagonalize_population, crtdst_is_given &
       &                         , sw_boltztrap, boltztrap_prefix, boltztrap_header &
#ifdef __EDA__
       &                         , sw_eda &
#endif
       &                         , boltztrap_version,num_extra_bands, ekmode &
       &                         , z_axis, connect_from
  use m_IterationNumbers,   only : iteration, first_iteration_of_this_job, iteration_electronic, iteration_ionic
  use m_Crystal_Structure,  only : altv,m_CS_set_altv_prim, a,b,c, op, nopr, rltv
  use m_Ionic_System,       only : natm, natm2, ityp, iwei, iatomn, pos &
       &                         , m_IS_pack_all_ions_in_uc &
       &                         , m_IS_set_natm_prim,m_IS_set_napt_prim,speciesname, cps
  use m_PseudoPotential,    only : ival, flg_paw
  use m_Files, only :              nfout,nfdos,nfchr,nfwfk,nfldos,nfelf &
       &                         , Nfvlc, nfcntn_bin_stm &
       &                         , m_Files_open_nfdos, m_Files_open_nfchr &
       &                         , m_Files_open_nfldos, m_Files_close_nfldos &
       &                         , m_Files_open_nfelf &
       &                         , m_Files_open_nfvlc &
       &                         , m_Files_open_nfcntn_bin_stm &
       &                         , m_Files_open_nfchr_pc &
       &                         , m_Files_skiptoend &
       &                         , m_Files_open_nfwfk
  use m_Kpoints, only :            way_ksample, kv3, kv3_ek, vkxyz, &
       &                           sw_force_kpt_inside_bz, vkxyz_ek, qwgt
  use m_Parallelization, only :    mype, ista_k, iend_k, map_e, map_k, map_z, MPI_CommGroup, myrank_e, myrank_k
  use m_PlaneWaveBasisSet,   only: kg,kgp,kg1,m_pwBS_wd_ngabc_etc
  use m_FFT,                 only: fft_box_size_WF, fft_box_size_CD &
       &                         , m_FFT_wd_box_sizes &
       &                         , m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work
  use m_ES_dos, only :             m_ESdos_gaussdistrib, m_ESdos_gaussdistrib_ek &
       &                         , m_ESdos_tetrahedral &
       &                         , m_ESdos_put_dos_weight &
       &                         , m_ESdos_write_dos_header &
       &                         , m_ESdos_alloc_dos_weight
  use m_Electronic_Structure,only: efermi, vbm, metalic_system, totch &
       &                         , check_if_metalic_flag &
       &                         , neordr, nrvf_ordr &
       &                         , m_ES_wd_zaj_small_portion &
       &                         , m_ES_orbital_population, m_ES_sym_comp, vlhxc_l &
       &                         , eko_ek
  use m_ES_nonlocal,         only: m_ES_phir_dot_WFs

  use m_ES_occup, only :           m_ESoc_check_if_metalic, m_ESoc_set_nEwindows_pc &
       &                         , m_ESoc_keep_occup,       m_ESoc_if_elec_state &
       &                         , m_ESoc_retrieve_occup,   m_ESoc_free_nEwindows &
       &                         , m_ESoc_substitute_occup, m_ESoc_check_if_metalic
!!$       &                         , m_ESoc_fermi_parabolic, m_ESoc_fermi_tetrahedron &
!!$       &                         , m_ESoc_fermi_ColdSmearing
  use m_ES_LHXC, only :            m_ESlhxc_potential
  use m_ES_IO,               only: m_ESIO_wd_WFn
  use m_ES_IO,               only: m_ESIO_wd_EigenValues_etc &
       &                         , m_ESIO_wd_vlhxc, m_ESIO_check_energy
  use m_Charge_Density, only :     chgsoft, chgq_l, chgq_l_pc, nEwindow_ek &
       &                         , m_CD_alloc_rspace_charge &
       &                         , m_CD_dealloc_rspace_charge &
       &                         , m_CD_rspace_charge &
       &                         , m_CD_keep_chgq_l &
       &                         , m_CD_rspace_put_headermark, m_CD_rspace_put_endmark &
       &                         , m_CD_retrieve_chgq, m_CD_keep_retrieve_hsr &
       &                         , m_CD_softpart, m_CD_hardpart &
       &                         , m_CD_softpart_rwf2 &
       &                         , m_CD_partial_charge_ek_finalize
  use m_XC_Potential,        only: vxc_l, m_XC_cal_potential
  use m_Ldos, only :               weiwsc,weilay, sw_save_ldos_weight &
       &                         , m_Ldos_cal &
       &                         , m_Ldos_cal_ek &
       &                         , m_Ldos_alloc_weiwsc_etc &
       &                         , m_Ldos_dealloc_weiwsc_etc &
       &                         , m_Ldos_what_is_n_total_ldos &
       &                         , m_Ldos_get_ldos_index  &
       &                         , m_Ldos_get_dos_weight  &
       &                         , m_Ldos_load_dos_weight
  use m_BerryPhase,  only        : m_BP_Polarization
  use m_BP_Properties, only      : m_BP_read_Berry_phase &
       &                         , m_BP_calc_Polarization &
       &                         , m_BP_calc_Effective_charge &
       &                         , m_BP_read_Effective_charge &
       &                         , m_BP_calc_piezoelectric_const
  use m_Phonon, only             : m_Phonon_read_forces &
       &                         , m_Phonon_calc_dynamical_matrix &
       &                         , m_Phonon_calc_vib_modes &
       &                         , m_Phonon_det_irr_rep &
       &                         , m_Phonon_calc_static_dielectric &
       &                         , m_Phonon_write_vib_modes &
       &                         , m_Phonon_write_epsilon &
       &                         , m_Phonon_calc_dos &
       &                         , m_Phonon_read_strain_forces &
       &                         , m_Phonon_calc_int_strain_piezo
  use m_ELF, only                : m_ELF_cnstrct_elf,m_ELF_wd_elf,m_ELF_alloc,m_ELF_dealloc
  use m_Dipole, only             : m_Dipole_calc, m_Dipole_layer, m_Dipole_wd_vlhxc
  use m_Wannier, only            : m_Wan_opt_Omega

  use m_Wannier90, only          : m_Wan90_gen_amat, m_Wan90_gen_mmat, m_Wan90_wd_eig, &
       &                           m_Wan90_write_unk, &
       &                           m_Wan90_gen_mmat_mpi, &
       &                           m_Wan90_wd_eig_mpi, m_Wan90_write_unk_mpi, &
       &                           m_Wan90_gen_mat_spn

!!$  use m_Files, only :              m_Files_open_nfdos, m_Files_close_nfdos

! ======================= Added by K. Tagami ================ 0.6
  use m_Phonon, only             : sw_calc_force
! ===========================================================

! =============================== added by K. Tagami =================== 11.0
  use m_Control_Parameters,     only : noncol, ndim_magmom, &
       &                               m_CtrlP_set_pchg_title_noncl
  use m_Charge_Density,         only :  m_CD_rspace_charge_noncl
  use m_Files,                  only : m_Files_open_nfchr_noncl, &
       &                               m_Files_open_nfchr_pc_noncl
  use m_Electronic_Structure,    only : m_ES_Orbital_population_noncl
  use m_ES_dos,                 only :  m_ESdos_gaussdistrib_noncl, &
       &                                m_ESdos_gaussdistrib_ek_noncl, &
       &                                m_ESdos_tetrahedral_noncl, &
       &                                m_ESdos_set_dos_weight_noncl, &
       &                                m_ESdos_put_dos_weight_noncl, &
       &                                m_ESdos_alloc_dos_wght_noncl
  use m_Ldos,                   only : m_Ldos_cal_noncl, &
       &                               m_Ldos_get_dos_weight_noncl
! ====================================================================== 11.0

! =============================== added by K. Tagami =================== 11.0
  use m_Charge_Density,         only :  m_CD_softpart_noncl, &
       &                                m_CD_hardpart_noncl
  use m_ES_LHXC,                only :  m_ESlhxc_potential_noncl
! ====================================================================== 11.0

! =============================== added by K. Tagami =================== 12.0Exp
  use m_Const_Parameters,    only  : ONE_BY_ONE, EK_CONVERGED, ALL_AT_ONCE
  use m_Control_Parameters,  only : fixed_charge_k_parallel
! ====================================================================== 12.0Exp

#ifndef DISABLE_CONSTRAINTS
  use m_constraints,         only : m_cnstr_reac_coords_variable,m_cnstr_get_id
  use m_velocity_verlet,     only : m_vv_get_curr_md_step
#endif

  use m_ES_RSB,              only : m_ES_RSB_doit,m_ES_RSB_alloc,m_ES_RSB_dealloc

! ================= KT_add ================== 13.0E
  use m_Control_Parameters, only : m_CtrlP_way_of_smearing
  use m_Const_Parameters,  only : Fermi_Dirac
  use m_ES_dos,            only : m_ESdos_FDiracDistrib, m_ESdos_FDiracDistrib_ek
  use m_ES_occup,          only : m_ESoc_count_charge_belowEF, &
       &                          m_ESoc_count_charge_belowEF_ek
! =========================================== 13.0E

! ======== KT_add ======= 13.0S
  use m_Control_Parameters,  only : sw_calc_core_energy, sw_write_soi_on_atoms
  use m_CoreLevel_Spectrum,  only : m_CLS_calc_core_energy
  use m_SpinOrbit_Potential,  only : m_SO_print_SOC_on_atoms
! ======================= 13.0S

! ==== KT_add === 13.1R
  use m_Control_Parameters,  only : sw_phonon_with_epsilon, sw_calc_dielectric_tensor
  use m_Raman,  only : m_Raman_read_dielec_tensors, &
       &               m_Raman_calc_Mat_dEps_dR
  use m_Phonon,  only : m_Phonon_print_Raman_susc
! =============== 13.1R

! ==== KT_add === 2014/09/01
  use m_Control_Parameters,  only : Orbital_decomp_mode, sw_calc_orbital_moment
  use m_CD_Mag_Moment,    only : m_CD_orb_decomposed_chg_magmom
  use m_OP_Moment,        only : m_OP_orb_decomposed_orb_magmom
! =============== 2014/09/01
! === KT_add === 13.2S
  use m_Crystal_Structure,  only : sw_spinorbit_second_variation
  use m_SpinOrbit_SecondVariation,  only : m_SO_calc_band_energy_socsv, &
       &                                   m_SO_init_second_variation, &
       &                                   m_SO_finalize_second_variation
! ============== 13.2S

  use m_Control_Parameters, only : sw_wf_squared_rspace, charge_filetype, &
       &                           sw_wf_integ_moment,  sw_calc_contact_density, &
       &                           sw_add_corecharge_rspace, sw_write_procar_file, &
       &                           sw_write_parity_file, &
       &                           sw_print_crystal_field_param, &
       &                           sw_calc_extfnv_correction, &
       &                           sw_band_unfolding, prepare_masked_compri, &
       &                           sw_calc_wf_orb_projection, &
       &                           sw_write_rotated_orbitals, sw_write_orb_dens_mat_file
  use m_ES_IO,              only : m_ESIO_wd_Wfn_squared_noncl, &
       &                           m_ESIO_wd_Wfn_integ_magmom, m_ESIO_wd_wfn_parity
  use m_OP_decomp_band,     only : m_OP_wd_PROCAR
  use m_PAW_ChargeDensity,  only : m_PAWCD_calc_contact_density
  use m_Crystal_Field,      only : m_CF_calc_CF_param

  use m_Band_Unfolding,     only : m_BU_calc_spectral_weight, &
       &                           m_BU_set_GVec_flag_refcell, &
       &                           m_BU_phir_dot_WFs
  use m_OP_rotation,        only : m_OP_rotation_col, m_OP_rotation_noncl,  &
       &                           m_OP_wd_phirt2_rotated,  m_OP_wd_WFn_orb_proj, &
       &                           m_OP_calc_dm_col, m_OP_calc_dm_noncl

  use m_Electronic_Structure, only : eko_l
  use m_routines, only : get_unused_unitnumber

  use m_ES_ChargeState,   only : sw_write_electrostatic_pot, &
       &                         m_ESCS_wd_electrostatic_pot, &
       &                         m_ESCS_wd_extfnv_correction
  use m_Potential_Average,  only :  sw_calc_pot_avg, m_PotAvg_calc_pot_on_atoms

#ifdef __EDA__
  use m_PlaneWaveEDA, only : execute_PlaneWaveEDA, m_PW_EDA_allocation, m_PW_EDA_deallocation
#endif

  use m_Control_Parameters,   only : sw_lband, sw_calc_wf_atom_decomposition, &
       &                             sw_calc_wf_layer_decomposition
  use m_Files,                only : nfwfk_local_decomp, &
       &                             m_Files_open_nfwfk_lband, &
       &                             m_Files_close_nfwfk_lband
  use mpi

  implicit none
!  include 'mpif.h'

  logical, intent(in) :: ignore_convergence
  logical :: Already_Converged
  integer :: iloop, dos_method_act, i, aldos_or_layerdos, k_total, n_total_ldos, icomponent, nEwindows, it
  integer :: ik,ib
  real(kind=DP), allocatable, dimension(:,:) :: dos_weight ! d(neg,kv3|kv3_ek)

! ========================== added by K. Tagami ================ 11.0
  real(kind=DP), allocatable, dimension(:,:,:) :: dos_weight_noncl
! ============================================================== 11.0

  real(kind=DP) :: emin, emax
  real(kind=DP), allocatable :: radius(:)
  integer :: nfrwf2

  if(ignore_convergence)then
     if(printable)then
#ifndef DISABLE_CONSTRAINTS
        if(driver==DRIVER_CONSTRAINT)then
          write(nfout,'(a,i8)') ' !** running the post-processor for ionic iteration : ',m_vv_get_curr_md_step()
        else
          write(nfout,'(a,i8)') ' !** running the post-processor for ionic iteration : ',iteration_ionic
        endif
#else
        write(nfout,'(a,i8)') ' !** running the post-processor for ionic iteration : ',iteration_ionic
#endif
     endif
  endif


  call tstatc_iter(iteration, first_iteration_of_this_job)
  it = iteration
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) it = iteration_electronic
  if(ipritiming0 >= 1 .and. iteration > first_iteration_of_this_job) call tstatc_wd(it)
  call tstatc_init

  if(iconvergence >= FORCE_CONVERGED .or. iconvergence_previous_job >= FORCE_CONVERGED .or. ignore_convergence) then
!!$     if(sw_positron == BULK .or. sw_positron == DEFECT) then
!!$        write(nfout,'(" --- m_positron_lifetime --")')
!!$        call m_positron_lifetime()
!!$     end if

     if((icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION ) &
          & .and. fixed_charge_k_parallel == ALL_AT_ONCE ) then
        if(.not.check_if_metalic_flag) call m_ESoc_check_if_metalic(nfout)
     end if

     if(ipri >= 1) write(nfout,'(" ! totch = ",f16.8)') totch
!!$     if(efermi < -1.d+10) call FermiEnergyLevel()

     if(ipri >= 1) then
        if(vbm > -9.9d10) then
           write(nfout,'(" !  totch, efermi, vbm = ",f16.8,2f18.10,"  <<Postprocessing>>")') &
                & totch, efermi, vbm
        else
           write(nfout,'(" !  totch, efermi      = ",f16.8,f18.10,"  <<Postprocessing>>")') &
                & totch, efermi
        end if
        if(check_if_metalic_flag) then
           if(metalic_system) then
              write(nfout,'(" !  metalic_system = true")')
           else
              write(nfout,'(" !  metalic_system = false")')
           end if
        end if
     end if

     if ( m_CtrlP_way_of_smearing() == Fermi_Dirac ) then
        if ( ekmode == ON ) then
           call m_ESoc_count_charge_belowEF_ek( nfout )
        else
           call m_ESoc_count_charge_belowEF( nfout )
        endif
     endif

! -------- Population ----
! ==== KT_add ===== 2014/09/01
     if ( Orbital_decomp_mode == 2 ) then
        call m_CD_orb_decomposed_chg_magmom
        if ( noncol .and. sw_calc_orbital_moment == ON ) then
           call m_OP_orb_decomposed_orb_magmom
        endif
     endif
! ================ 2014/09/01

     if(sw_orb_popu == ON) then
!        if ( ekmode == OFF .and. sw_calc_force == OFF .and. &
!             &  (icond /=FIXED_CHARGE .and. icond/=FIXED_CHARGE_CONTINUATION) ) then
        if ( ekmode == OFF .and. sw_calc_force == OFF) then
           call m_ES_phir_dot_WFs(nfout)
           call m_ES_sym_comp(nfout)
           if ( noncol ) then
              call m_ES_orbital_population_noncl( nfout )
              if ( sw_write_orb_dens_mat_file == ON ) call m_OP_calc_dm_noncl
              if ( sw_diagonalize_population == ON ) then
                 call m_OP_rotation_noncl(nfout)
                 if ( sw_write_rotated_orbitals == ON ) call m_OP_wd_phirt2_rotated
              endif
           else
              call m_ES_orbital_population(nfout)
              if ( sw_write_orb_dens_mat_file == ON ) call m_OP_calc_dm_col
              if ( sw_diagonalize_population == ON ) then
                 call m_OP_rotation_col(nfout)
                 if ( sw_write_rotated_orbitals == ON ) call m_OP_wd_phirt2_rotated
              endif
           endif
           if ( sw_write_procar_file == ON ) then
              call m_OP_wd_PROCAR
           endif
        endif

     end if

! --------- DOS ------
     if (sw_dos == ON) then
        if ( noncol ) then
           call m_ESdos_alloc_dos_wght_noncl()
           call m_ESdos_set_dos_weight_noncl()
           call calc_totaldos_noncl()
        else
           call m_ESdos_alloc_dos_weight()
           call calc_totaldos()
        endif
     end if

! ---- LDOS ----------
!
     if ( sw_dos == ON .and. sw_ldos == ON ) then
        if(.not.crtdst_is_given) then
           call m_CntrlP_set_crtdst(a,b,c) ! -> crtdst_aldos,crtdst_winlay
        end if
        if ( noncol ) then
           call calc_localdos_noncl()
        else
           call calc_localdos()
        endif
     end if
!
! ------ Charge Distrib ----
!
     if (sw_charge_rspace == ON) then
        if(Already_Converged()) then
           call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l, VXC_AND_EXC)  ! -> vxc_l, afft
        end if

        if ( noncol ) then
           call calc_spatial_chg_distrib_noncl(.false.)
           if ( sw_add_corecharge_rspace == ON ) then
              call calc_spatial_chg_distrib_noncl(.true.)
           endif
        else
           call calc_spatial_chg_distrib(.false.)
           if ( sw_add_corecharge_rspace == ON ) then
              call calc_spatial_chg_distrib(.true.)
           endif
        endif
     endif

! ------ Parital Charge ---
!
     if ( sw_charge_rspace == ON .and. sw_partial_charge == ON .and. &
       &  (icond==INITIAL .or. icond==CONTINUATION)) then
        if ( noncol ) then
           call calc_partial_charge_noncl()
        else
           call calc_partial_charge()
        endif
     end if

     if ( sw_calc_pot_avg == ON ) call m_PotAvg_calc_pot_on_atoms
     if ( sw_charge_rspace == ON .and. icond == FIXED_CHARGE .and. ekmode == ON ) then
       call calc_spatial_chg_distrib_ek()
     endif
     if ( sw_charge_rspace == ON .and. sw_partial_charge == ON .and. &
       &  icond == FIXED_CHARGE .and. ekmode == ON ) then
           call calc_partial_charge_ek()
     endif

! =============================== added by K. Tagami ============= 11.0
     if ( noncol ) goto 550
! ================================================================ 11.0

!! test test test
     if(sw_rsb_test==ON) then
        call m_ES_RSB_alloc()
        call m_ES_RSB_doit(.true.,fermi=.false.)
        call m_ES_RSB_dealloc()
     endif
!! test test test

     if(sw_wf_rspace == ON) then
        call m_FFT_alloc_WF_work
        do ik=1,kv3
           do ib=1,neg
              if(m_ESIO_check_energy(ik,ib)) then
                 call m_Files_open_nfwfk(nspin,ik,ib)
                 call m_ESIO_wd_WFn(nfout,nfwfk,ik,ib)
              end if
           end do
        end do
        call m_FFT_dealloc_WF_work
     end if

     if(sw_rwf2 == ON) then
        allocate(radius(nr))
        call m_FFT_alloc_WF_work
        do i=1,nr
           radius(i) = (i-1)*rmax/(nr-1)
        end do
        nfrwf2 = 2000
        call m_CD_softpart_rwf2(nfout,nfrwf2,kv3,ik_rwf2,ib_rwf2,center,nr,radius)
        deallocate(radius)
        call m_FFT_dealloc_WF_work
     end if

     if(sw_elf == ON) then
        if(icond == INITIAL .or. icond == CONTINUATION) then
           call m_ELF_alloc()
           call m_ELF_cnstrct_elf(nfout)
           call m_Files_open_nfelf
           call m_ELF_wd_elf(nfout,nfelf)
           call m_ELF_dealloc()
        end if
     end if

     if(sw_berry_phase == ON) then
        if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
           call m_BP_Polarization
        end if
     end if

! =================================== added by K. Tagami ================= 11.0
550  continue
! ======================================================================== 11.0

     if ( sw_wf_squared_rspace == ON ) then
        if ( noncol ) call m_ESIO_wd_Wfn_squared_noncl(nfout)
     endif

! ------ STM or Work function ---
!
     if (sw_fine_STM_simulation == ON) then
        if ( noncol ) then
           call write_potential_for_STM_noncl
        else
           call write_potential_for_STM
        endif
        write(nfout, *) '!!STM:   template input'
        call write_template_nfstminput(nfout)
     endif


! --- Dipole ?
!
     if (sw_dipole_correction == ON) then
        call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgq_l, VXC_AND_EXC )
                                                              ! chgsoft -> vxc_l
        call m_Dipole_wd_vlhxc( nfout, vxc_l )
     end if

! =============================== added by K. Tagami ============= 11.0
     if ( noncol ) goto 620
! ================================================================ 11.0

     if (sw_wannier == ON) then
        call m_Wan_opt_Omega(nfout)
     end if

620  continue

     if (sw_wannier90 == ON) then
        call m_Wan90_gen_amat(nfout)
#ifndef USE_MPI_WANNIER90
        call m_Wan90_gen_mmat(nfout)
        call m_Wan90_wd_eig(nfout)
        if ( sw_write_unk_file == ON ) call m_Wan90_write_unk(nfout)
#else
        call m_Wan90_gen_mmat_mpi(nfout)
        call m_Wan90_wd_eig_mpi(nfout)
        if ( sw_write_unk_file == ON ) call m_Wan90_write_unk_mpi(nfout)
        if ( noncol ) call m_Wan90_gen_mat_spn
#endif
     end if

! =================================== added by K. Tagami ================= 11.0
650  continue
! ======================================================================== 11.0

  end if

! =============================== added by K. Tagami ============= 11.0
!  if ( noncol ) goto 750
! ================================================================ 11.0

  if( polar_prop == POLARIZATION .or. &
    & polar_prop == EFFECTIVE_CHARGE .or. &
    & polar_prop == PIEZOELECTRIC_CONST) then
     call m_BP_read_Berry_phase(.false.)
     call m_BP_read_Berry_phase(.true.)
  end if

  if( polar_prop == POLARIZATION ) then
     call m_BP_calc_Polarization
  else if( polar_prop == EFFECTIVE_CHARGE .or. sw_lo_to_splitting == ON &
    & .or. sw_lattice_dielectric_tensor == ON &
    & .or. sw_int_strain_piezo_tensor == ON) then
     if(sw_vibrational_modes == ON) then
        call m_CS_set_altv_prim()
        call m_IS_set_natm_prim()
        call m_IS_set_napt_prim()
! == KT_add === 2014/06/20
     else if ( polar_prop == EFFECTIVE_CHARGE ) then
        call m_IS_set_napt_prim()
! ============= 2014/06/20
     endif
     call m_BP_calc_Effective_charge
  else if( polar_prop == PIEZOELECTRIC_CONST) then
     call m_BP_calc_piezoelectric_const
  end if

  if(sw_vibrational_modes == ON .and. sw_phonon_oneshot == OFF &
    & .and. (iconvergence >= FORCE_CONVERGED .or. sw_calc_force == OFF)) then  ! asms
     call m_Phonon_read_forces()
     call m_CS_set_altv_prim()
     call m_IS_set_natm_prim()
     call m_IS_set_napt_prim()
     call m_Phonon_calc_dynamical_matrix()
     call m_Phonon_calc_vib_modes()
     call m_Phonon_det_irr_rep()
     if(phonon_method == PHONON_GAMMA .and. sw_lattice_dielectric_tensor == ON) &
       & call m_Phonon_calc_static_dielectric()
     call m_Phonon_write_vib_modes()
     if(phonon_method == PHONON_GAMMA .and. sw_dielectric_function == ON) &
       & call m_Phonon_write_epsilon()
     if(phonon_method == PHONON_DOS) call m_Phonon_calc_dos()
     if(phonon_method == PHONON_GAMMA .and. sw_int_strain_piezo_tensor == ON) then
	call m_Phonon_read_strain_forces()
        call m_Phonon_calc_int_strain_piezo()
     end if

! === KT_add === 13.1R
!     if ( phonon_method == PHONON_GAMMA .and. sw_phonon_with_epsilon == ON &
!          &                             .and. sw_calc_dielectric_tensor == OFF ) then
     if ( phonon_method == PHONON_GAMMA  &
          &        .and. sw_phonon_with_epsilon == ON &
          &        .and. ( iconvergence >= FORCE_CONVERGED &
          &                .or. sw_calc_dielectric_tensor == OFF ) ) then
        call m_Raman_read_dielec_tensors
        call m_Raman_calc_Mat_dEps_dR
        call m_Phonon_print_Raman_susc
     endif
! ============== 13.1R

  end if

! =================================== added by K. Tagami ================= 11.0
750 continue
! ======================================================================== 11.0
  if ( sw_print_crystal_field_param == ON ) call m_CF_calc_CF_param

  if (sw_dipole == ON) then
     if(sw_layered == ON) then
        call m_Dipole_layer(nfout)
     else
        call m_Dipole_calc(nfout)
     end if
  end if
  if ( sw_write_soi_on_atoms == ON ) then
     call m_SO_print_SOC_on_atoms
  endif
  if ( ekmode == OFF .and. sw_write_parity_file == ON ) call m_ESIO_wd_wfn_parity

  if ( sw_write_electrostatic_pot == ON ) then
     call m_ESCS_wd_electrostatic_pot
  endif
  if ( sw_calc_extfnv_correction == ON ) then
     call m_ESCS_wd_extfnv_correction
  endif

! ========= KT_add === 13.0S
  if ( icond < 2 .and. sw_calc_core_energy == ON ) call m_CLS_calc_core_energy
! ==================== 13.0S
! === KT_add === 13.2S
  if ( .not. noncol .and. sw_spinorbit_second_variation == ON ) then
     if ( ekmode == OFF ) then
        call m_SO_init_second_variation
        call m_SO_calc_band_energy_socsv
        call m_SO_finalize_second_variation
     endif
  endif
! ============== 13.2S

  if ( noncol ) then
     if ( ekmode == OFF .and. sw_wf_integ_moment == ON ) then
        call m_ESIO_wd_Wfn_integ_magmom
     endif
  endif
  if ( sw_calc_contact_density == ON ) then
     if ( flg_paw ) call m_PAWCD_calc_contact_density
  endif

  if ( ekmode == OFF ) then
     if ( sw_band_unfolding == ON ) then
        if ( sw_force_kpt_inside_bz == ON ) call m_BU_set_GVec_flag_refcell
        call m_BU_calc_spectral_weight
     endif
     if ( sw_orb_popu == ON .and. sw_calc_wf_orb_projection == ON ) then
        if ( sw_band_unfolding == ON .and. prepare_masked_compri == ON ) then
           call m_BU_phir_dot_WFs
        else
           call m_ES_phir_dot_WFs( nfout )
        endif
        call m_OP_wd_WFn_orb_proj
     end if
     if ( sw_lband == ON ) then
        if(.not.crtdst_is_given) then
           call m_CntrlP_set_crtdst(a,b,c) ! -> crtdst_aldos,crtdst_winlay
        end if
        call calc_lband
     endif
  endif

  if(sw_boltztrap == ON)then
    call output_data_for_bt(nfout,boltztrap_prefix,boltztrap_header,boltztrap_version)
  endif

#ifdef __EDA__
  if(sw_eda == ON) then
    call m_PW_EDA_allocation()
    call execute_PlaneWaveEDA()
    call m_PW_EDA_deallocation()
  endif
#endif

!!$  call tstatc_wd0()

contains
!!$  subroutine FermiEnergyLevel()
!!$    integer :: way_of_smearing
!!$    way_of_smearing = m_CtrlP_way_of_smearing()
!!$    if(way_of_smearing == PARABOLIC) then
!!$       if(ipri>=1) write(nfout,'(" way_of_smearing = PARABOLIC <<Postprocessing>>")')
!!$       call m_ESoc_fermi_parabolic(nfout)
!!$!!$  else if(way_of_smearing == MP) then
!!$!!$     call fermi_mesfessel_paxton(nfout)
!!$    else if(way_of_smearing == TETRAHEDRON) then
!!$       if(ipri>=1) write(nfout,'(" way_of_smearing = TETRAHEDRON <<Postprocessing>>")')
!!$        call m_ESoc_fermi_tetrahedron(nfout)
!!$    else if(way_of_smearing == COLD) then
!!$       if(ipri>=1) write(nfout,'(" way_of_smearing = COLD <<Postprocessing>>")')
!!$        call m_ESoc_fermi_ColdSmearing(nfout)
!!$    end if
!!$  end subroutine FermiEnergyLevel

  subroutine calc_totaldos()

    if(.not.check_if_metalic_flag.or.ignore_convergence) call m_ESoc_check_if_metalic(nfout)

    if(.not.ignore_convergence)then
        call m_Files_open_nfdos()
    else
#ifndef DISABLE_CONSTRAINTS
        if(driver==DRIVER_CONSTRAINT)then
          if(m_cnstr_reac_coords_variable())then
             call m_Files_open_nfdos(iter=m_vv_get_curr_md_step(),reacid=m_cnstr_get_id())
          else
             call m_Files_open_nfdos(iter=m_vv_get_curr_md_step())
          endif
        else
#endif
          call m_Files_open_nfdos(iter=iteration_ionic)
#ifndef DISABLE_CONSTRAINTS
        endif
#endif
    endif

    dos_method_act = dos_method
    if(dos_method_act == TETRAHEDRON .and. way_ksample /= MESH) then
       dos_method_act = Gauss_distrib_func
    end if

    if(ipridos>=1) then
       if(dos_method_act == TETRAHEDRON) then
          write(nfout,'(" dos_method_act = TETRAHEDRON")')

! ======================== KT_add ============= 13.0E
       else if ( dos_method_act == Fermi_Dirac ) then
          write(nfout,'(" dos_method_act = FERMI_DIRAC")')
! ============================================== 13.0E
       else
          write(nfout,'(" dos_method_act = Gauss_distrib_funcN")')
       end if
    end if

    icomponent = TOTAL
    if(dos_method_act == Gauss_distrib_func) then
       if(icond == INITIAL .or. icond == CONTINUATION ) then
          call m_ESdos_gaussdistrib(nfdos,icomponent)
       else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
!  ================================ Modified by K. Tagami ====== 0.5
!          if(sw_orb_popu /= ON) &
!          & call m_ESdos_gaussdistrib_ek(nfdos,icomponent)
       if( fixed_charge_k_parallel == ONE_BY_ONE ) then
           call m_ESdos_gaussdistrib_ek(nfdos,icomponent)
       else
           call m_ESdos_gaussdistrib(nfdos,icomponent)
       endif
! ===============================================================
       end if

! ================== KT_add  ================== 13.0E
    else if ( dos_method_act == Fermi_Dirac ) then
       if(icond == INITIAL .or. icond == CONTINUATION ) then
          call m_ESdos_FDiracDistrib(nfdos,icomponent)
       else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
           call m_ESdos_FDiracDistrib_ek(nfdos,icomponent)
       end if
! ============================================= 13.0E

    else  ! dos_method_act == Tetrahedral
       if(icond == INITIAL .or. icond == CONTINUATION ) then
          call m_ESdos_tetrahedral(nfdos,icomponent,mode=SCF)
       else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
!  ============== Modified by K. Tagami ================== 0.5
!!          if(sw_orb_popu /= ON) &
!!!!          & call m_ESdos_tetrahedral(nfdos,icomponent,mode=EK)

! =================================== KT ======================== 12.0Exp
!             call m_ESdos_tetrahedral(nfdos,icomponent,mode=EK)
!
          if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
             call m_ESdos_tetrahedral(nfdos,icomponent,mode=EK)
          else
             call m_ESdos_tetrahedral(nfdos,icomponent,mode=SCF)
          endif
! ==================================================================== 12.0Exp

! ======================================================= 0.5
       end if
    end if
!!$     call m_Files_close_nfdos()
  end subroutine calc_totaldos

  subroutine calc_totaldos_noncl()

    if(.not.check_if_metalic_flag.or.ignore_convergence) call m_ESoc_check_if_metalic(nfout)

    call m_Files_open_nfdos()

    dos_method_act = dos_method
    if(dos_method_act == TETRAHEDRON .and. way_ksample /= MESH) then
       dos_method_act = Gauss_distrib_func
    end if

    if(ipridos>=1) then
       if(dos_method_act == TETRAHEDRON) then
          write(nfout,'(" dos_method_act = TETRAHEDRON")')
       else
          write(nfout,'(" dos_method_act = Gauss_distrib_funcN")')
       end if
    end if

    icomponent = TOTAL

    if (dos_method_act == Gauss_distrib_func) then

       if (icond == INITIAL .or. icond == CONTINUATION ) then
          call m_ESdos_gaussdistrib_noncl(nfdos,icomponent)
       else if (icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
           call m_ESdos_gaussdistrib_ek_noncl(nfdos,icomponent)
       end if

    else                      ! dos_method_act == Tetrahedral
       if (icond == INITIAL .or. icond == CONTINUATION ) then
          call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=SCF )
       else if (icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then

          if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
             call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=EK )
          else
             call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=SCF )
          endif

       end if
    end if

  end subroutine calc_totaldos_noncl

  subroutine calc_localdos()
    character(len=40) :: tagwords

    if (ekmode == OFF) then
       call m_Ldos_alloc_weiwsc_etc()
       call m_Ldos_cal( nfldos, .true. ) !-> weiwsc, weilay
    else
       call m_Ldos_alloc_weiwsc_etc()
       call m_Ldos_cal_ek( nfldos ) !-> weiwsc, weilay
    end if
    call m_Files_close_nfldos()

    if (sw_cal_ldos == ON) then
       if(sw_dos /= ON) call calc_totaldos()

       call m_Files_open_nfldos()
       dos_method_act = ldos_method
       if(dos_method_act == TETRAHEDRON .and. way_ksample /= MESH) then
          dos_method_act = Gauss_distrib_func
       end if

       n_total_ldos = m_Ldos_what_is_n_total_ldos()

       if(n_total_ldos > 0 ) then
          if(ekmode == OFF) then
             k_total = kv3
          else if(ekmode == ON) then
             k_total = kv3_ek
          end if
          allocate(dos_weight(neg,k_total))
       end if

       if(ekmode==ON) call m_Ldos_load_dos_weight()
       do i=1, n_total_ldos
!!$              write(nfout,'(" i = ",i8, " <<Postprocessing>>")') i
          call m_Ldos_get_ldos_index(i,aldos_or_layerdos,icomponent,tagwords) ! -> aldos_or_layerdos,icomponent
          call m_Ldos_get_dos_weight( aldos_or_layerdos, icomponent, nfldos,&
               &                      neg, k_total, dos_weight )     !  weiwsc or weilay -> dos_weight
          call m_ESdos_put_dos_weight( neg, k_total, dos_weight )    ! -> dos_weight (in m_ES_dos)
          call m_ESdos_write_dos_header( nfdos, aldos_or_layerdos, icomponent, tagwords )

          if (dos_method_act == Gauss_distrib_func) then
             if(icond == INITIAL .or. icond == CONTINUATION ) then
                call m_ESdos_gaussdistrib(nfdos,icomponent)
             else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
                call m_ESdos_gaussdistrib_ek(nfdos,icomponent)
             end if

! ============================= KT_add ================= 13.0E
          else if ( dos_method_act == Fermi_Dirac ) then
             if(icond == INITIAL .or. icond == CONTINUATION ) then
                call m_ESdos_FDiracDistrib(nfdos,icomponent)
             else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
                call m_ESdos_FDiracDistrib_ek(nfdos,icomponent)
             end if
! ====================================================== 13.0Ea

          else  ! dos_method_act = Tetrahedral
             if (icond == INITIAL .or. icond == CONTINUATION ) then
                call m_ESdos_tetrahedral(nfdos,icomponent,mode=SCF)
             else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then

! =============================================== KT ================= 12.0Exp
!                call m_ESdos_tetrahedral(nfdos,icomponent,mode=EK)
!
               if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
                   call m_ESdos_tetrahedral(nfdos,icomponent,mode=EK)
                else
                   call m_ESdos_tetrahedral(nfdos,icomponent,mode=SCF)
                endif
! ===================================================================== 12.0Exp

             end if
          end if
       end do
    end if

    if (n_total_ldos > 0 ) deallocate(dos_weight)
    call m_Ldos_dealloc_weiwsc_etc()

  end subroutine calc_localdos

  subroutine calc_localdos_noncl()
    character(len=40) :: tagwords

    if (ekmode == OFF) then
       call m_Ldos_alloc_weiwsc_etc()
       call m_Ldos_cal_noncl( nfldos, .true. )
    endif
    call m_Files_close_nfldos()

    if (sw_cal_ldos == ON) then

       if (sw_dos /= ON) then
          call m_ESdos_alloc_dos_wght_noncl()
          call m_ESdos_set_dos_weight_noncl()
          call calc_totaldos_noncl()
       endif

       call m_Files_open_nfldos()
       dos_method_act = ldos_method
       if(dos_method_act == TETRAHEDRON .and. way_ksample /= MESH) then
          dos_method_act = Gauss_distrib_func
       end if

       n_total_ldos = m_Ldos_what_is_n_total_ldos()

       if(n_total_ldos > 0 ) then
          if(ekmode == OFF) then
             k_total = kv3
          else if(ekmode == ON) then
             k_total = kv3_ek
          end if
          allocate(dos_weight_noncl(neg,k_total,ndim_magmom))
       end if

       do i = 1, n_total_ldos
          call m_Ldos_get_ldos_index( i, aldos_or_layerdos, icomponent, tagwords )
                               ! -> aldos_or_layerdos, icomponent

          call m_Ldos_get_dos_weight_noncl( aldos_or_layerdos, icomponent, &
               &                            nfldos, neg, k_total, &
               &                            dos_weight_noncl )
          call m_ESdos_put_dos_weight_noncl( neg, k_total, dos_weight_noncl )

          call m_ESdos_write_dos_header( nfdos, aldos_or_layerdos, icomponent, tagwords )

          if (dos_method_act == Gauss_distrib_func) then
             if(icond == INITIAL .or. icond == CONTINUATION ) then
                call m_ESdos_gaussdistrib_noncl( nfdos, icomponent )
             else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
                call m_ESdos_gaussdistrib_ek_noncl( nfdos, icomponent )
             end if
          else                     ! dos_method_act = Tetrahedral
             if (icond == INITIAL .or. icond == CONTINUATION ) then
                call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=SCF )
             else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then

                if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
                   call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=EK )
                else
                   call m_ESdos_tetrahedral_noncl( nfdos, icomponent, mode=SCF )
                endif

             end if
          end if
       end do
    end if

    if ( n_total_ldos > 0 ) deallocate(dos_weight_noncl)
    call m_Ldos_dealloc_weiwsc_etc()

  end subroutine calc_localdos_noncl

  subroutine calc_lband()
    call m_Ldos_alloc_weiwsc_etc()
    call m_Files_open_nfwfk_lband(icond)
    if ( noncol ) then
!       stop "Does not work"
       call m_Ldos_cal_noncl( nfwfk_local_decomp, .true. ) !-> weiwsc, weilay
    else
       call m_Ldos_cal( nfwfk_local_decomp, .true. ) !-> weiwsc, weilay
    endif
    call m_Files_close_nfwfk_lband()
    call m_Ldos_dealloc_weiwsc_etc()
  end subroutine calc_lband

  subroutine calc_spatial_chg_distrib( flg_add_core )
    logical, intent(in) :: flg_add_core
    integer :: iloop, iloop2

    if (icond == INITIAL .or. icond == CONTINUATION) then
       call m_CD_alloc_rspace_charge()
       do iloop = 1, nspin
! ======================= KT_add === 2014/06/07
          if ( nspin == 2 .and. sw_spin_magmom_rspace == ON ) then
             iloop2 = -iloop
          else
             iloop2 = iloop
          endif
! =================================== 2014/06/07
          if ( flg_add_core .and. iloop2 == -2 ) cycle

          if(.not.ignore_convergence)then
             call m_Files_open_nfchr( nspin, iloop2, add_core=flg_add_core )
          else
#ifndef DISABLE_CONSTRAINTS
             if(driver==DRIVER_CONSTRAINT)then
                if(m_cnstr_reac_coords_variable())then
                  call m_Files_open_nfchr(nspin,iloop2,iter=m_vv_get_curr_md_step(),reacid=m_cnstr_get_id())
                else
                  call m_Files_open_nfchr(nspin,iloop2,iter=m_vv_get_curr_md_step())
                endif
             else
#endif
                call m_Files_open_nfchr(nspin,iloop2,iter=iteration_ionic)
#ifndef DISABLE_CONSTRAINTS
             endif
#endif
          endif
          call m_CD_rspace_charge(nspin,iloop2,nfchr,nfout,add_core=flg_add_core)
       end do
       call m_CD_dealloc_rspace_charge()
    endif

  end subroutine calc_spatial_chg_distrib

  subroutine calc_spatial_chg_distrib_ek()
    use m_IterationNumbers, only : m_Iter_reset_nk_in_the_process
    use m_Charge_Density, only : m_CD_softpart, m_CD_hardpart, chgq_l
    use m_ES_IO, only : m_ESIO_rd_next_wfs_ek
    use m_ES_nonlocal, only : m_ES_betar_dot_WFs_4_each_k
    use m_Parallelization, only : ista_kngp, iend_kngp
    use m_Files, only : m_Files_open_nfzaj_kall, nfzaj_kall
    use m_Control_Parameters, only : af
    use m_Parallelization, only : ista_kngp, iend_kngp
    use m_Electronic_Structure, only : vbm, efermi
    use m_ES_occup, only : m_ESoc_occup_under_ef
    integer :: iloop, iloop2
    integer :: nk, is
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_l_sum
    call m_CD_alloc_rspace_charge()
    allocate(chgq_l_sum(ista_kngp:iend_kngp,kimg,nspin));chgq_l_sum = 0.d0
    call m_Files_open_nfzaj_kall()
    chgq_l_sum = 0.d0
    rewind(nfzaj_kall)
    call m_Iter_reset_nk_in_the_process()
    do nk=1, kv3_ek, kv3
      call KpointNumber_Setting2()
      call Preparation_ek()
      call Preparation_for_mpi_ek()
      call PseudoPotential_ek()
      call m_ESoc_occup_under_ef(nfout, nk)
      do is=1, nspin, af+1
        do ik=is, kv3+is-nspin, nspin
          call m_ESIO_rd_next_wfs_ek(ik, nfout, nfzaj_kall)
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_betar_dot_WFs_4_each_k(nfout, ik)
        enddo
      enddo
      call m_CD_softpart(nfout, kv3)
      call m_CD_hardpart(nfout, kv3)
      chgq_l_sum =  chgq_l_sum + chgq_l
    enddo
    chgq_l = chgq_l_sum
    do is=1, nspin, af+1
      call m_Files_open_nfchr(nspin, is)
      call m_CD_rspace_charge(nspin,is,nfchr,nfout)
    enddo
    call m_CD_dealloc_rspace_charge()
    deallocate(chgq_l_sum)
  end subroutine calc_spatial_chg_distrib_ek

  subroutine calc_partial_charge_ek()
    use m_Control_Parameters, only : af
    use m_ES_occup, only : m_ESoc_set_nEwindows_pc_ek
    use m_Files, only : m_Files_open_nfzaj_kall, nfzaj_kall
    use m_IterationNumbers, only : m_Iter_reset_nk_in_the_process
    use m_ES_IO, only : m_ESIO_rd_next_wfs_ek
    use m_ES_nonlocal, only : m_ES_betar_dot_WFs_4_each_k
    use m_Charge_Density, only : m_CD_softpart, m_CD_hardpart, chgq_l
    use m_Parallelization, only : ista_kngp, iend_kngp
    integer :: i, iloop, iloop2
    integer :: nk, is
    real(kind=DP), allocatable, dimension(:,:,:) :: chgq_l_pc

    call m_ESoc_set_nEwindows_pc_ek(nfout,nEwindows)
    call m_ESoc_keep_occup()
    call m_CD_keep_chgq_l()
    call m_CD_keep_retrieve_hsr(.true.)
    call m_CntrlP_keep_charge_title()
    call m_Files_open_nfzaj_kall()

    allocate(chgq_l_pc(ista_kngp:iend_kngp,kimg,nspin));chgq_l_pc = 0.d0
    do i = 1, nEwindows
       if(m_ESoc_if_elec_state(nfout,i,emin,emax) /= YES) cycle
       chgq_l_pc = 0.d0
       rewind(nfzaj_kall)
       call m_Iter_reset_nk_in_the_process()
       do nk=1, kv3_ek, kv3
         call KpointNumber_Setting2()
         call Preparation_ek()
         call Preparation_for_mpi_ek()
         call PseudoPotential_ek()
         call m_ESoc_substitute_occup(nfout,i,nk)

         do is=1, nspin, af+1
           do ik=is, kv3+is-nspin, nspin
             call m_ESIO_rd_next_wfs_ek(ik, nfout, nfzaj_kall)
             if(map_k(ik) /= myrank_k) cycle
             call m_ES_betar_dot_WFs_4_each_k(nfout, ik)
           enddo
         enddo
         call m_CD_softpart(nfout, kv3)
         call m_CD_hardpart(nfout, kv3)
         chgq_l_pc =  chgq_l_pc + chgq_l
       enddo
       chgq_l = chgq_l_pc

       call m_CD_alloc_rspace_charge()

       do iloop = 1, nspin
! =========================== KT_add === 2014/06/07
          if ( nspin == 2 .and. sw_spin_magmom_rspace == ON ) then
             iloop2 = -iloop
          else
             iloop2 = iloop
          endif
! ===================================== 2014/06/07
          call m_CntrlP_set_pcharge_title(nspin,iloop2,i,emin,emax)
          if(partial_charge_filetype == SEPARATE) then
             call m_Files_open_nfchr_pc(nspin,iloop2,i)
          else
             call m_Files_open_nfchr(nspin,iloop2)
             call m_Files_skiptoend(nfchr)
             call m_CD_rspace_put_headermark(nfchr,nspin,iloop2,i)
          end if
          call m_CD_rspace_charge(nspin,iloop2,nfchr,nfout)
          if(partial_charge_filetype == INTEGRATED) & ! /= SEPARATE
               & call m_CD_rspace_put_endmark(nfchr)
       end do
       call m_CD_dealloc_rspace_charge()
    end do
    call m_ESoc_retrieve_occup()
    call m_CD_retrieve_chgq()
    call m_CD_keep_retrieve_hsr(.false.)
    call m_CntrlP_retrieve_charge_title()
    call m_ESoc_free_nEwindows()
    deallocate(chgq_l_pc)
  end subroutine calc_partial_charge_ek

  subroutine calc_spatial_chg_distrib_noncl( flg_add_core )
    logical, intent(in) :: flg_add_core
    integer :: iloop

    if (icond == INITIAL .or. icond == CONTINUATION) then
       call m_CD_alloc_rspace_charge()
       do iloop = 1, ndim_magmom
          if ( flg_add_core .and. iloop > 1 ) cycle

          if (.not.ignore_convergence)then
             call m_Files_open_nfchr_noncl(iloop,add_core=flg_add_core)
          else
             call m_Files_open_nfchr_noncl(iloop,iter=iteration_ionic)
          endif
          call m_CD_rspace_charge_noncl( iloop, nfchr, nfout, charge_filetype, &
               &                         add_core=flg_add_core )
       end do
       call m_CD_dealloc_rspace_charge()
    endif

  end subroutine calc_spatial_chg_distrib_noncl

  subroutine calc_partial_charge
    integer :: i, iloop, iloop2

    if(icond == INITIAL .or. icond == CONTINUATION) then

       if(.not.check_if_metalic_flag.or.ignore_convergence) then
          call m_ESoc_check_if_metalic(nfout)
          if(ipri >= 1) then
!!$                    write(nfout,'(" !  totch, efermi, vbm = ",2f20.12,d20.8,"  <<Postprocessing>>")') &
!!$                         & totch, efermi, vbm
             if(metalic_system) then
                write(nfout,'(" !  metalic_system = true")')
             else
                write(nfout,'(" !  metalic_system = false")')
             end if
          end if
       end if

       call m_ESoc_set_nEwindows_pc(nfout,nEwindows)
       call m_ESoc_keep_occup()
       call m_CD_keep_chgq_l()
       call m_CD_keep_retrieve_hsr(.true.)
       call m_CntrlP_keep_charge_title()

       do i = 1, nEwindows
          if(m_ESoc_if_elec_state(nfout,i,emin,emax) /= YES) cycle
          call m_ESoc_substitute_occup(nfout,i)

          call m_CD_softpart(nfout,kv3)
          call m_CD_hardpart(nfout,kv3)

          call m_CD_alloc_rspace_charge()

          do iloop = 1, nspin
! =========================== KT_add === 2014/06/07
             if ( nspin == 2 .and. sw_spin_magmom_rspace == ON ) then
                iloop2 = -iloop
             else
                iloop2 = iloop
             endif
! ===================================== 2014/06/07
             call m_CntrlP_set_pcharge_title(nspin,iloop2,i,emin,emax)
             if(partial_charge_filetype == SEPARATE) then
                if(.not.ignore_convergence)then
                   call m_Files_open_nfchr_pc(nspin,iloop2,i)
                else
#ifndef DISABLE_CONSTRAINTS
                   if(driver==DRIVER_CONSTRAINT)then
                      if(m_cnstr_reac_coords_variable())then
                         call m_Files_open_nfchr_pc(nspin,iloop2,i,iter=m_vv_get_curr_md_step(),reacid=m_cnstr_get_id())
                      else
                         call m_Files_open_nfchr_pc(nspin,iloop2,i,iter=m_vv_get_curr_md_step())
                      endif
                   else
#endif
                      call m_Files_open_nfchr_pc(nspin,iloop2,i,iter=iteration_ionic)
#ifndef DISABLE_CONSTRAINTS
                   endif
#endif
                endif
             else
                if(.not.ignore_convergence)then
                   call m_Files_open_nfchr(nspin,iloop2)
                else
#ifndef DISABLE_CONSTRAINTS
                   if(driver==DRIVER_CONSTRAINT)then
                      if(m_cnstr_reac_coords_variable())then
                         call m_Files_open_nfchr(nspin,iloop2,iter=m_vv_get_curr_md_step(),reacid=m_cnstr_get_id())
                      else
                         call m_Files_open_nfchr(nspin,iloop2,iter=m_vv_get_curr_md_step())
                      endif
                   else
#endif
                      call m_Files_open_nfchr(nspin,iloop2,iter=iteration_ionic)
#ifndef DISABLE_CONSTRAINTS
                   endif
#endif
                endif
                call m_Files_skiptoend(nfchr)
                call m_CD_rspace_put_headermark(nfchr,nspin,iloop2,i)
             end if
             call m_CD_rspace_charge(nspin,iloop2,nfchr,nfout)
             if(partial_charge_filetype == INTEGRATED) & ! /= SEPARATE
                  & call m_CD_rspace_put_endmark(nfchr)
          end do

          call m_CD_dealloc_rspace_charge()

       end do
       call m_ESoc_retrieve_occup()
       call m_CD_retrieve_chgq()
       call m_CD_keep_retrieve_hsr(.false.)
       call m_CntrlP_retrieve_charge_title()
       call m_ESoc_free_nEwindows()
    end if

  end subroutine calc_partial_charge

  subroutine calc_partial_charge_noncl
    integer :: i, iloop

    if(icond == INITIAL .or. icond == CONTINUATION) then

       if(.not.check_if_metalic_flag.or.ignore_convergence) then
          call m_ESoc_check_if_metalic(nfout)
          if(ipri >= 1) then
!!$                    write(nfout,'(" !  totch, efermi, vbm = ",2f20.12,d20.8,"  <<Postprocessing>>")') &
!!$                         & totch, efermi, vbm
             if(metalic_system) then
                write(nfout,'(" !  metalic_system = true")')
             else
                write(nfout,'(" !  metalic_system = false")')
             end if
          end if
       end if

       call m_ESoc_set_nEwindows_pc(nfout,nEwindows)
       call m_ESoc_keep_occup()
       call m_CD_keep_chgq_l()
       call m_CD_keep_retrieve_hsr(.true.)
       call m_CntrlP_keep_charge_title()

       do i = 1, nEwindows
          if(m_ESoc_if_elec_state(nfout,i,emin,emax) /= YES) cycle
          call m_ESoc_substitute_occup(nfout,i)

          call m_CD_softpart_noncl(nfout,kv3)
          call m_CD_hardpart_noncl(nfout,kv3)

          call m_CD_alloc_rspace_charge()

          do iloop = 1, ndim_magmom
             call m_CtrlP_set_pchg_title_noncl(iloop,i,emin,emax)
             if(partial_charge_filetype == SEPARATE) then
                if(.not.ignore_convergence)then
                   call m_Files_open_nfchr_pc_noncl(iloop,i)
                else
                   call m_Files_open_nfchr_pc_noncl(iloop,i,iteration_ionic)
                endif
             else
                if(.not.ignore_convergence)then
                   call m_Files_open_nfchr_noncl(iloop)
                else
                   call m_Files_open_nfchr_noncl(iloop,iteration_ionic)
                endif
                call m_Files_skiptoend(nfchr)

                call m_CD_rspace_put_headermark(nfchr,nspin,iloop,i)
             end if
             call m_CD_rspace_charge_noncl(iloop,nfchr,nfout,charge_filetype)
             if(partial_charge_filetype == INTEGRATED) & ! /= SEPARATE
                  & call m_CD_rspace_put_endmark(nfchr)
          end do

          call m_CD_dealloc_rspace_charge()

       end do
       call m_ESoc_retrieve_occup()
       call m_CD_retrieve_chgq()
       call m_CD_keep_retrieve_hsr(.false.)
       call m_CntrlP_retrieve_charge_title()
       call m_ESoc_free_nEwindows()
    end if

  end subroutine calc_partial_charge_noncl


  subroutine write_potential_for_STM
    integer :: ismax

    ismax = nspin

    if(.not.check_if_metalic_flag.or.ignore_convergence) then
       call m_ESoc_check_if_metalic(nfout)
       if(ipri >= 1) then
!!$                 write(nfout,'(" !  totch, efermi, vbm = ",2f20.12,d20.8,"  <<Postprocessing>>")') &
!!$                      & totch, efermi, vbm
          if(metalic_system) then
             write(nfout,'(" !  metalic_system = true")')
          else
             write(nfout,'(" !  metalic_system = false")')
          end if
       end if
    end if

    if (ipri >= 1) call wd_fine_STM_parameters( ismax )
    call m_Files_open_nfvlc()
    if (mype==0) rewind nfvlc

    if (sw_deficit_charge == 1) then
!       call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgq_l, VXC_AND_EXC )
!       call m_ESlhxc_potential(nfout,chgq_l,vxc_l)      ! chq_l, vxc_l -> vlhxc
           if(sw_add_xc_to_vloc==ON)then
              call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC) ! chgsoft -> vxc_l
           else
              vxc_l = 0.0d0
           endif
           if(sw_xc_only==ON)then
               vlhxc_l(:,:,:) = vxc_l(:,:,:)
           else
               call m_ESlhxc_potential(nfout,chgq_l,vxc_l) ! chgsoft, vxc_l -> vlhxc
           endif

    else
           ! sw_deficit_charge==0
!       call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgsoft, VXC_AND_EXC )
!       call m_ESlhxc_potential(nfout,chgsoft,vxc_l)     ! chgsoft, vxc_l -> vlhxc
           if(sw_add_xc_to_vloc==ON)then
              call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgsoft,VXC_AND_EXC) ! chgsoft -> vxc_l
           else
              vxc_l = 0.d0
           endif
           if(sw_xc_only==ON)then
               vlhxc_l(:,:,:) = vxc_l(:,:,:)
           else
               call m_ESlhxc_potential(nfout,chgsoft,vxc_l) ! chgsoft, vxc_l -> vlhxc
           endif

    endif
    call m_ESIO_wd_vlhxc( nfvlc, ismax )

    if(ipri >= 3) call check_neordr_nrvf_ordr()
    call m_Files_open_nfcntn_bin_stm()
    if(mype==0) rewind nfcntn_bin_stm

    call wd_ArraySize_Parameters_For_STM( nfcntn_bin_stm, ismax )
    call m_ESIO_wd_EigenValues_etc(nfcntn_bin_stm,.false.,totch_flag=OFF)
    call m_pwBS_wd_ngabc_etc(nfcntn_bin_stm)
    call m_FFT_wd_box_sizes(nfcntn_bin_stm)
    call wd_Header_For_Cube(nfcntn_bin_stm)

  end subroutine write_potential_for_STM

  subroutine write_potential_for_STM_noncl
    integer :: ismax

    ismax = ndim_magmom

    if (ipri >= 1) call wd_fine_STM_parameters( ismax )
    call m_Files_open_nfvlc()
    if (mype==0) rewind nfvlc

    if (sw_deficit_charge == 1) then
       call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgq_l, VXC_AND_EXC )
                                                        ! chg_l -> vxc_l
       call m_ESlhxc_potential_noncl( nfout, chgq_l,vxc_l )
    else
       call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgsoft, VXC_AND_EXC )
                                                        ! chgsoft -> vxc_l
       call m_ESlhxc_potential_noncl( nfout, chgsoft,vxc_l )
    endif

    call m_ESIO_wd_vlhxc( nfvlc, ismax )

    if (ipri >= 3) call check_neordr_nrvf_ordr()
    call m_Files_open_nfcntn_bin_stm()
    if (mype==0) rewind nfcntn_bin_stm

    call wd_ArraySize_Parameters_For_STM( nfcntn_bin_stm, ismax )
    call m_ESIO_wd_EigenValues_etc(nfcntn_bin_stm,.false.,totch_flag=OFF)
    call m_pwBS_wd_ngabc_etc(nfcntn_bin_stm)
    call m_FFT_wd_box_sizes(nfcntn_bin_stm)
    call wd_Header_For_Cube(nfcntn_bin_stm)

  end subroutine write_potential_for_STM_noncl

  subroutine write_template_nfstminput(nfout)
    integer, intent(in) :: nfout
    real(kind=DP),allocatable,dimension(:,:) :: cps_full
    integer, allocatable,dimension(:) :: ityp_full
    integer :: iat, iz=3
    real(kind=DP) :: zh, posz, cpoint, minz, maxz

    if(mype == 0) then
      iz = z_axis
      write(nfout, '(i8,i2)') kv3, nspin 
      do ik = 1, kv3
         write(nfout,'(i3,3f8.4,3x,3f8.4,f8.4)') ik &
              &   , (vkxyz(ik,i,CARTS),i=1,3),(vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
      enddo
      write(nfout, '(i8, f20.10)') neg, altv(iz,iz)
      write(nfout, '(a)') "-1.0 0.0"
      if(kimg==2) then
        write(nfout, '(i8)') fft_box_size_CD(iz,1)
      else
        write(nfout, '(i8)') fft_box_size_CD(iz,1)/2
      endif

      allocate(cps_full(natm2,3))
      allocate(ityp_full(natm2))
      cps_full(1:natm,1:3) = cps(1:natm,1:3)
      ityp_full(1:natm) = ityp(1:natm)
      call rplcps(cps_full, ityp_full, 1, natm2, natm, iwei)
      maxz = 0.d0
      minz = altv(iz, iz)
      do iat=1, natm2
        posz = cps_full(iat, iz)
        if(posz>altv(iz,iz)*0.5d0) then
          posz = posz-altv(iz, iz)
        endif
        if(posz>maxz) maxz = posz
        if(posz<minz) minz = posz
      enddo
      cpoint = maxz+connect_from
      zh = 0.5d0*(minz+altv(iz, iz)-maxz)+maxz

      write(nfout, '(i8, i8)') nint((cpoint/altv(iz,iz))*fft_box_size_CD(iz, 1)), &
     &                         nint((zh/altv(iz,iz))*dble(fft_box_size_CD(iz, 1)))
      write(nfout, '(i2)') iz
      write(nfout, '(a)') '0.10 0.30 100'
      write(nfout, '(a)') '0.01'
    
      deallocate(ityp_full,cps_full)
    endif
  end subroutine write_template_nfstminput

  subroutine wd_fine_STM_parameters( ismax )
    integer, intent(in) :: ismax

    integer :: ik
    write(nfout,'(" !!STM:    kg(kng)   = ",i8)') kg
    write(nfout,'(" !!STM:    kgp(kngp) = ",i8)') kgp
    write(nfout,'(" !!STM:    kg1(kng1) = ",i8)') kg1
    write(nfout,'(" !!STM:    neg(keg)  = ",i8)') neg
    write(nfout,'(" !!STM:    kimg      = ",i8)') kimg
    write(nfout,'(" !!STM:    fft_box_size_WF(1,1)(knl)  = ",i8)') fft_box_size_WF(1,1)
    write(nfout,'(" !!STM:    fft_box_size_WF(2,1)(knm)  = ",i8)') fft_box_size_WF(2,1)
    write(nfout,'(" !!STM:    fft_box_size_WF(3,1)(knn)  = ",i8)') fft_box_size_WF(3,1)
    write(nfout,'(" !!STM:    fft_box_size_WF(1,0)(kid)  = ",i8)') fft_box_size_WF(1,0)
    write(nfout,'(" !!STM:    fft_box_size_CD(1,1)(knlp) = ",i8)') fft_box_size_CD(1,1)
    write(nfout,'(" !!STM:    fft_box_size_CD(2,1)(knmp) = ",i8)') fft_box_size_CD(2,1)
    write(nfout,'(" !!STM:    fft_box_size_CD(3,1)(knnp) = ",i8)') fft_box_size_CD(3,1)
    write(nfout,'(" !!STM:    fft_box_size_CD(1,0)(kidp) = ",i8)') fft_box_size_CD(1,0)
    write(nfout,'(" !!STM:    kv3(knv3) = ",i8)') kv3

! ====================== modiifed by K. Tagami ======================= 11.0
!    write(nfout,'(" !!STM:    nspin(kspin) = ",i8)') nspin
    write(nfout,'(" !!STM:    ndim of magmom (kspin) = ",i8)') ismax
! =================================================================== 11.0

    write(nfout,'(" !!STM:  == k-points ==")')
    write(nfout,'(" !!STM: ik",8x,"CARTS",22x,"PUCS",16x,"QITG")')
    do ik = 1, kv3
       write(nfout,'(" !!STM: ",i3,3f8.4,3x,3f8.4,f8.4)') ik &
            &   , (vkxyz(ik,i,CARTS),i=1,3),(vkxyz(ik,i,BUCS) ,i=1,3),qwgt(ik)
    enddo
  end subroutine wd_fine_STM_parameters

  subroutine check_neordr_nrvf_ordr()
    integer :: i
    do i = ista_k, iend_k
       write(nfout,'(" !!Postprocessing -- neordr --, ik = ",i8)') i
       write(nfout,'(" !!Postprocessing ",10i6)') (neordr(iloop,i),iloop=1,neg)
    end do
    do i = ista_k, iend_k
       write(nfout,'(" !!Postprocessing -- nrvf_ordr --, ik = ",i8)') i
       write(nfout,'(" !!Postprocessing ",10i6)') (nrvf_ordr(iloop,i),iloop=1,neg)
    end do
  end subroutine check_neordr_nrvf_ordr

  subroutine wd_ArraySize_Parameters_For_STM( nf_bin, ismax )
    integer, intent(in) :: ismax

    integer, intent(in) :: nf_bin
    if(mype == 0) then
       write(nf_bin) kg
       write(nf_bin) kgp
       write(nf_bin) kg1
       write(nf_bin) neg
       write(nf_bin) kimg
       write(nf_bin) fft_box_size_WF(1,1)
       write(nf_bin) fft_box_size_WF(2,1)
       write(nf_bin) fft_box_size_WF(3,1)
       write(nf_bin) fft_box_size_WF(1,0)
       write(nf_bin) fft_box_size_CD(1,1)
       write(nf_bin) fft_box_size_CD(2,1)
       write(nf_bin) fft_box_size_CD(3,1)
       write(nf_bin) fft_box_size_CD(1,0)
       write(nf_bin) kv3
! ====================== modiifed by K. Tagami ======================= 11.0
!       write(nf_bin) nspin
       write(nf_bin) ismax
! ==================================================================== 11.0
    endif
  end subroutine wd_ArraySize_Parameters_For_STM

  subroutine wd_Header_For_Cube(nf_bin)
    integer, intent(in) :: nf_bin
    integer :: i, m
    real(kind=DP) :: x,y,z
    real(kind=DP),allocatable,dimension(:,:) :: cps_full
    integer, allocatable,dimension(:) :: ityp_full
    !!$real(kind=DP), dimension(3) :: r_wk    !!! K.Mae 040315

    if(mype == 0) then

      x = 0.d0; y = 0.d0; z = 0.d0
      write(nf_bin) natm2, x,y,z, natm

      do i = 1, 3
         !!$do m = 1, 3
         !!$   ucret = unit_conv_byname( altv(m,i), r_wk(m), 'bohr', 'angstrom' )
         !!$end do
         !!$write(nf_bin) fft_box_size_CD(i,1), r_wk(1:3)/dble(fft_box_size_CD(i,1))
         write(nf_bin) fft_box_size_CD(i,1), altv(1:3,i)/dble(fft_box_size_CD(i,1))
      end do

      allocate(cps_full(natm2,3))
      allocate(ityp_full(natm2))
      call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
! -
!	write(*,*) 'natom 2 = ',natm2
!	stop
      do i = 1, natm2
         m = ityp_full(i)
! =========================== Modified by J, Koga ============== 0.2
!         write(nf_bin) iatomn(m), ival(m), cps_full(i,1:3)
          write(nf_bin) nint(iatomn(m)), ival(m), cps_full(i,1:3)
!          write(*,*) ' I = ',i, nint(iatomn(m)), ival(m), cps_full(i,3)
! ===========================================================
      end do
      deallocate(ityp_full,cps_full)

    endif

  end subroutine wd_Header_For_Cube

  subroutine output_data_for_bt(nfout,prefix,header,version)
    integer, intent(in) :: nfout
    character(len=*), intent(in) :: prefix,header
    integer, intent(in) :: version
    integer :: un
    integer :: i,j,k,l,ik,ie,ierr,no
    integer :: mpi_comm1,mpi_comm2
    real(DP),allocatable, dimension(:,:) :: e_wk
    real(DP), dimension(3,3) :: op_br
    real(kind=DP),allocatable,dimension(:,:) :: cps_full
    integer, allocatable, dimension(:) :: ityp_full
    integer,allocatable,dimension(:) :: neordr_ek
    integer :: ib, jb, ibo, jbo
    real(kind=DP), parameter :: delta = 1.d-12
    integer :: numk
    real(kind=DP), pointer :: vkxyz_p(:,:,:)

    if(ekmode == OFF) then
      numk = kv3
      vkxyz_p => vkxyz
    else if (ekmode == ON) then
      numk = kv3_ek
      vkxyz_p => vkxyz_ek
    endif
    allocate(e_wk(neg,numk));e_wk=0.d0
    if(ekmode == OFF) then
      do ik = 1, kv3                                     ! MPI
         if(map_k(ik) /= myrank_k) cycle                 ! MPI
         do ie = 1, neg                                  ! MPI
            if(map_e(ie) /= myrank_e) cycle              ! MPI
            e_wk(ie,ik) = eko_l(map_z(ie),ik)            ! MPI
         end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,e_wk,neg*kv3,mpi_double_precision,mpi_sum,MPI_CommGroup,ierr)
    else if(ekmode == ON) then
      allocate(neordr_ek(neg))
      do ik = 1, kv3_ek
         neordr_ek = (/(ib,ib=1,neg)/)
         do ib = 1,neg-1
            do jb = ib+1, neg
               ibo = neordr_ek(ib)
               jbo = neordr_ek(jb)
               if(eko_ek(jbo,ik) < eko_ek(ibo,ik)-delta) then
                  neordr_ek(jb) = ibo
                  neordr_ek(ib) = jbo
               end if
            end do
         end do
! Substitution eko_ek for eko in the order of thier values
         do ie = 1, neg
            e_wk(ie,ik) = eko_ek(neordr_ek(ie),ik)
         end do
      end do
      deallocate(neordr_ek)
    endif

    if(mype == 0)then
    un = get_unused_unitnumber()
    write(nfout,'(a,i5)') ' !** output data for BoltzTraP; prefix: '//trim(prefix)//', version: ',version
    if (version == 1) then
       write(nfout,'(a)') ' !** template BoltzTraP input written to: '//trim(prefix)//'.intrans'
       open(un,file=trim(prefix)//'.intrans',status='replace')
       write(un,'(a)') 'GENE                     # Format of DOS'
       write(un,'(3i4,f10.5)') 0,0,0,0.0
       write(un,'(3f10.5,i8,a)') efermi*2.0d0,0.0005,0.4,int(totch) &
          &, ' # Fermilevel (Ry), energygrid, energy span around Fermilevel, number of electrons'
       write(un,'(a)') 'CALC                     # CALC (calculate expansion coeff), NOCALC read from file'
       write(un,'(a)') '5                        # lpfac, number of latt-points per k-point'
       write(un,'(a)') 'BOLTZ                    # run mode (only BOLTZ is supported)'
       write(un,'(a)') '.15                      # (efcut) energy range of chemical potential'
       write(un,'(a)') '800. 50.                 # Tmax, temperature grid'
       write(un,'(a)') '-1.                      '// &
        &'# energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)'
       write(un,'(a)') 'HISTO'
       write(un,'(a)') '0 0 0 0 0                # For scattering model undocumented'
       write(un,'(a)') '2                        # number of'
       write(un,'(a)') '1E20 -1E20               #  fixed doping levels in cm-3'
       close(un)

       write(nfout,'(a)') ' !** eigenvalue data is written to: '//trim(prefix)//'.energy'
       open(un,file=trim(prefix)//'.energy', status='replace')
       write(un,'(a)') trim(header)
       write(un,*) numk
       do i=1,numk
           write(un,'(3f10.5,i5)') vkxyz_p(i,1:3,BUCS),neg-num_extra_bands
           do j=1,neg-num_extra_bands
              write(un,'(f10.5)') e_wk(j,i)*2.d0 ! to Rydberg units
           enddo
       enddo
       close(un)
       write(nfout,'(a)') ' !** structure  data is written to: '//trim(prefix)//'.struct'
       open(un,file=trim(prefix)//'.struct', status='replace')
       write(un,'(a)') trim(header)
       write(un,'(3f15.5)') altv(1:3,1)
       write(un,'(3f15.5)') altv(1:3,2)
       write(un,'(3f15.5)') altv(1:3,3)
       write(un,'(i5)') nopr
       do no=1,nopr
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
          write(un,'(9i3)') ((nint(op_br(j,i)),j=1,3),i=1,3)
       enddo
       close(un)
    else if (version == 2)then
       !version 2 format : energy in Ryd, length in Bohr
       write(nfout,'(a)') ' !** eigenvalue data is written to: '//trim(prefix)//'.energy'
       open(un,file=trim(prefix)//'.energy',status='replace')
       write(un,'(a)') trim(header)
       write(un,'(i8,i3,f10.5)') numk,nspin,efermi*2.d0
       do i=1,numk
           write(un,'(3f10.5,i5)') vkxyz_p(i,1:3,BUCS),neg-num_extra_bands
           do j=1,neg-num_extra_bands
              write(un,'(f10.5)') e_wk(j,i)*2.d0
           enddo
       enddo
       close(un)
       write(nfout,'(a)') ' !** structure  data is written to: '//trim(prefix)//'.struct'
       open(un,file=trim(prefix)//'.structure',status='replace')
       write(un,'(a)') trim(header)
       write(un,'(3f15.5)') altv(1:3,1)
       write(un,'(3f15.5)') altv(1:3,2)
       write(un,'(3f15.5)') altv(1:3,3)
       write(un,'(i8)') natm2
       allocate(cps_full(natm2,3))
       allocate(ityp_full(natm2))
       call m_IS_pack_all_ions_in_uc(ityp_full,cps_full)
       do i=1,natm2
          write(un,'(a,3f15.5)') speciesname(ityp_full(i)),cps_full(i,1:3)
       enddo
       close(un)
    endif
    endif
    deallocate(e_wk)
  end subroutine output_data_for_bt

end subroutine Postprocessing
!$$#endif

subroutine Postprocessing_k
!$$#ifndef PARA3D
  !  T. Yamasaki,  Feb. 2004
  use m_Const_Parameters, only :   ON, OFF
  use m_Control_Parameters, only : ekmode, sw_ldos, noncol, sw_lband
  use m_IterationNumbers, only :   nk_in_the_process

  use m_Ldos, only :               m_Ldos_cal &
       &                         , m_Ldos_wd_natm2_and_totch &
       &                         , m_Ldos_alloc_weiwsc_etc &
       &                         , m_Ldos_dealloc_weiwsc_etc
  use m_Ldos,                 only :   m_Ldos_cal_noncl

! === KT_add === 13.2S
  use m_Crystal_Structure,  only : sw_spinorbit_second_variation
  use m_SpinOrbit_SecondVariation,  only : m_SO_calc_band_energy_socsv
! ============== 13.2S

  use m_Control_Parameters,    only :  sw_wf_squared_rspace, sw_wf_integ_moment, &
       &                               sw_orb_popu, sw_band_unfolding, &
       &                               prepare_masked_compri, &
       &                               sw_calc_wf_orb_projection, sw_write_parity_file, &
       &                               sw_pdos, sw_partial_charge, sw_charge_rspace
  use m_Files, only : nfout, nfzaj_kall, m_Files_open_nfzaj_append_kall, &
       &              nfwfk_local_decomp
  use m_ES_IO, only : m_ESIO_wd_Wfn_squared_noncl, m_ESIO_wd_Wfn_integ_magmom, &
       &              m_ESIO_wd_wfn_parity, m_ESIO_wd_WFs

  use m_OP_rotation,  only : m_OP_wd_WFn_orb_proj, m_OP_wd_phirt2_rotated

  use m_ES_nonlocal,  only : m_ES_phir_dot_WFs, m_ES_add_betar_dot_WFs

  use m_Kpoints,  only : kv3, vkxyz, sw_force_kpt_inside_bz

  use m_Band_Unfolding,       only : m_BU_calc_spectral_weight, &
       &                             m_BU_set_GVec_flag_refcell, &
       &                             m_BU_phir_dot_WFs, &
       &                             m_BU_betar_dot_WFs2, m_BU_set_unfolding_active
  use m_Electronic_Structure,  only : m_ES_keep_retrieve_fsr

  use m_Control_Parameters,  only : sw_use_add_proj, sw_lband
  use m_Excitation,   only : m_XI_set_ilmt_val, m_XI_set_ngpt_XI, m_XI_wd_WAVEDERF, &
       &                     sw_write_wavederf_file, m_XI_set_mat_dipole_corr_vv
  use m_ValenceBand_Spectrum, only : m_VBS_set_data_ppc_from_pp
  use m_ES_occup,             only : m_ES_occ_with_fermi_fixed

  use m_Control_Parameters, only : m_CntrlP_set_crtdst, crtdst_is_given
  use m_Crystal_Structure,  only : a, b, c
  use m_Files,              only : nfwfk_local_decomp, &
       &                           m_Files_open_nfwfk_lband, m_Files_close_nfwfk_lband

  use m_Parallelization,  only :   mype

  implicit none

  logical :: First = .true.

!  if(ekmode == ON .and. sw_ldos == ON) then
!     if(nk_in_the_process == 1)  call m_Ldos_wd_natm2_and_totch()
!     call m_Ldos_alloc_weiwsc_etc()
!
!     if ( noncol ) then
!        call m_Ldos_cal_noncl()
!     else
!        call m_Ldos_cal()
!     endif
!
!     call m_Ldos_dealloc_weiwsc_etc()
!  end if

! === KT_add === 13.2S
  if ( .not. noncol .and. sw_spinorbit_second_variation == ON ) then
     call m_SO_calc_band_energy_socsv
  endif
! ============== 13.2S

  if ( ekmode == ON .and. sw_write_parity_file == ON ) call m_ESIO_wd_wfn_parity
  if ( sw_band_unfolding == ON ) then
     if ( sw_force_kpt_inside_bz == ON ) call m_BU_set_GVec_flag_refcell
     call m_BU_calc_spectral_weight
  endif

  if ( noncol ) then
     if ( sw_wf_squared_rspace == ON ) call m_ESIO_wd_Wfn_squared_noncl(nfout)
     if ( sw_wf_integ_moment == ON ) call m_ESIO_wd_Wfn_integ_magmom
  endif

  if ( sw_orb_popu == ON .and. sw_calc_wf_orb_projection == ON ) then
     if ( sw_band_unfolding == ON .and. prepare_masked_compri == ON ) then
        call m_BU_phir_dot_WFs
     else
        call m_ES_phir_dot_WFs( nfout )
     endif
     call m_OP_wd_WFn_orb_proj
  endif

  if ( sw_write_wavederf_file == ON ) then
     call m_ES_occ_with_fermi_fixed
     if ( First ) then
        call m_XI_set_ilmt_val
        call m_XI_set_ngpt_XI
        call m_VBS_set_data_ppc_from_pp
        call m_XI_set_mat_dipole_corr_vv
        First = .false.
     endif
     if ( sw_use_add_proj == ON ) call m_ES_add_betar_dot_WFs(nfout)
     call m_XI_wd_WAVEDERF
  end if

  if ( sw_lband == ON ) then
     if ( sw_band_unfolding == ON ) then
        call phase_error_with_msg(6,'sw_lband does not work with sw_band_unfolding', &
             &                    __LINE__,__FILE__)
     endif

     if(.not.crtdst_is_given) then
        call m_CntrlP_set_crtdst(a,b,c) ! -> crtdst_aldos,crtdst_winlay
     end if
     call m_Ldos_alloc_weiwsc_etc()
     if ( mype == 0 ) then
        if ( ekmode == ON ) then
           if ( nk_in_the_process == 1 ) then
              call m_Files_open_nfwfk_lband(2)
           else
              call m_Files_open_nfwfk_lband(3)
           endif
        else
           call m_Files_open_nfwfk_lband(2)
        endif
    endif

     if ( sw_band_unfolding == ON ) then
        call m_ES_keep_retrieve_fsr( .true. )
        call m_BU_set_unfolding_active( .true. )
        call m_BU_betar_dot_WFs2(nfout)
     endif
     if ( noncol ) then
        call m_Ldos_cal_noncl( nfwfk_local_decomp, .true. ) !-> weiwsc, weilay
     else
        call m_Ldos_cal( nfwfk_local_decomp, .true. ) !-> weiwsc, weilay
     endif

     if ( sw_band_unfolding == ON ) then
        call m_ES_keep_retrieve_fsr( .false. )
        call m_BU_set_unfolding_active( .false. )
     endif

     call m_Files_close_nfwfk_lband()
     call m_Ldos_dealloc_weiwsc_etc()
  endif

  if(sw_ldos == ON .or. sw_pdos == ON &
  &                .or. sw_partial_charge == ON &
  &                .or. sw_charge_rspace  == ON ) then
    call m_Files_open_nfzaj_append_kall()
    call m_ESIO_wd_WFs(nfout, nfzaj_kall, .false., rew=.false.)
  endif

!$$#endif
end subroutine Postprocessing_k
