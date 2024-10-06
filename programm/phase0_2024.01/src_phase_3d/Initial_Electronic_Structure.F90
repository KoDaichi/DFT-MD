#ifdef JRCATFFT_WS
#define _INCLUDE_EXX_
#elif FFTW3
#define _INCLUDE_EXX_
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: Initial_Electronic_Structure, Initial_WaveFunctions_ek
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
subroutine Initial_Electronic_Structure
! $Id: Initial_Electronic_Structure.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Const_Parameters,   only : Gauss_distrib_func, from_wave_functions&
       &                         , INITIAL, CONTINUATION, FIXED_CHARGE &
       &                         , FIXED_CHARGE_CONTINUATION, ON, OFF, EXECUT &
       &                         , COORDINATE_CONTINUATION &
       &                         , by_random_numbers, by_matrix_diagon, FILE &
       &                         , by_pseudo_atomic_orbitals &
       &                         , Valence_plus_PC_Charge, VXC_AND_EXC, EXC_ONLY &
       &                         , ONE_BY_ONE, ALL_AT_ONCE, DP, YES &
       &                         , Partial_Core_Charge, PT_CONTROL
  use m_IterationNumbers,   only : iteration, nk_in_the_process
  use m_Files,              only : nfout,nfchgt,nfzaj,nfcntn_bin,nfeng,nfefermi &
       &                         , F_ZAJ_in_partitioned, F_CHGT_in_partitioned &
       &                         , F_CNTN_BIN_in_partitioned &
       &                         , nfchr &
       &                         , m_Files_open_nfchr &
       &                         , m_Files_open_nfcntn_bin &
       &                         , m_Files_skiptoend, m_Files_open_nfzaj, m_Files_open_nfchgt &
       &                         , m_Files_open_nfefermi, m_Files_close_nfefermi &
       &                         , nfcntn_bin_paw &
       &                         , file_existence_contfiles &
       &                         , m_Files_check_nfzaj_existence, m_Files_check_nfchgt_existence &
       &                         , m_Files_check_file_existence &
       &                         , m_Files_open_nfcntn_bin_paw,m_Files_close_nfcntn_bin_paw &
       &                         , m_Files_nfcntn_bin_paw_exists
  use m_Control_Parameters, only : ipri, iprichargedensity, initial_chg, icond, nspin, intzaj &
       &                         , evaluation_eko_diff &
       &                         , skip_alloc_phonon, sw_phonon, sw_calc_force, neg, neg_previous &
       &                         , delta_eigenvalue_cond_is_given, sw_dipole_correction, sw_hubbard &
       &                         , ekmode, fixed_charge_k_parallel, ipriparallel, printable &
       &                         , sw_initial_charge_rspace &
       &                         , m_CtrlP_set_neg_properly, m_CntrlP_set_neg, m_CntrlP_set_meg &
       &                         , m_Cntrlp_set_davidson_size &
       &                         , sw_hybrid_functional, sw_screening_correction &
       &                         , sw_external_potential, sw_fef, initial_occmat,kimg &
       &                         , sw_berry_phase, sw_rsb, sw_eval_epc_on_fftmesh &
       &                         , sw_initial_es, imdalg, sw_interpolate_charge, sw_interpolate_wfs &
       &                         , sw_pdos, sw_optimize_blocking_parameters &
       &                         , sw_precalculate_phase_vnonlocal
  use m_Kpoints,            only : kv3
  use m_PlaneWaveBasisSet,  only : kg1, m_pwBS_alloc_ylm_l,kgp,ngabc, kg1_prev, kgp_prev, m_pwBS_rd_prev_pws
  use m_Total_Energy,      only  : m_TE_set_etotal_old
  use m_Electronic_Structure,only: totch, m_ES_gtotch &
       &                         , m_ES_energy_eigen_values_ext_3D &
       &                         , m_ES_alloc_vlhxc, m_ES_alloc_vlhxcQ &
       &                         , m_ES_alloc_zaj_etc &
       &                         , m_ES_alloc_eko_ek, m_ES_alloc_eko1 &
       &                         , m_ES_cpeko &
       &                         , m_ES_alloc_Dhub,vlhxc_l &
#ifdef FFTW3
       &                         , m_ES_cp_zaj_prev_to_zaj &
       &                         , m_ES_gen_zaj_from_prev_zaj
#else
       &                         , m_ES_cp_zaj_prev_to_zaj
#endif
  use m_ES_ortho,           only : m_ES_modified_gram_schmidt
  use m_ES_wf_extrpl,       only : m_ES_wf_extrpl_alloc
  use m_ES_nonlocal,        only : m_ES_betar_dot_WFs_3D
  use m_ES_IO,              only : m_ESIO_rd_Efermi, m_ESIO_rd_EigenValues_etc, m_ESIO_rd_WFs
  use m_ES_initialWF,       only : m_ESIW_by_randomnumbers_3D  &
       &                         , m_ESIW_by_randomnumbers0_3D
  use m_ES_initialWF,       only : m_ESIW_by_atomic_orbitals
  use m_ES_WF_by_SDorCG,    only : m_ESsd_alloc_dzajn2, m_ESsd_alloc_zaj_old
  use m_Charge_Density,     only : chgq_l &
       &                         , m_CD_alloc_hsr &
       &                         , m_CD_rd_hsr &
       &                         , m_CD_rd_chgq  &
       &                         , m_CD_initial_CD_by_file_rspace &
       &                         , m_CD_adjust_spindensity  &
#ifdef FFTW3
       &                         , m_CD_cp_chgq_prev_to_chgq &
       &                         , m_CD_gen_chgq_from_prev_chgq
#else
       &                         , m_CD_cp_chgq_prev_to_chgq
#endif
  use m_Crystal_Structure, only  : sw_magnetic_constraint,  univol, &
       &                           sw_spinorbit_second_variation
! <--
  use m_XC_Potential,      only  : vxc_l, exc

  use m_Orbital_Population,only  : m_OP_rd_occ_mat, m_OP_mix_om, m_OP_occ_mat_init, m_OP_cp_ommix_to_omold
  use m_FiniteElectricField,only : m_FEF_Constract_of_ftq
  use m_Parallelization,      only : m_Parallel_dealloc_mpi_elec, m_Parallel_init_mpi_elec_3D &
                       &           , np_g1k_prev, ista_kngp,iend_kngp,nrank_e,mype
  use m_Charge_Density,       only : m_CD_wd_chgq_l_small_portion   &
 &                                 , m_CD_cp_chgq_to_chgqo
  use m_Charge_Density,       only : m_CD_initial_CD_by_Gauss_func, m_CD_initial_CD_by_Gauss_kt
  use m_Charge_Density,       only : m_CD_softpart_3D, m_CD_hardpart
  use m_Dipole,               only : m_Dipole_vdip_alloc_3D
  use m_Electronic_Structure, only : vlhxc_l             &
 &                                 , zaj_l,   zaj_ball    &
 &                                 , fsr_l                 &
 &                                 , fsi_l                 &
 &                                 , eko_l              &
 &                                 , occup_l            &
 &                                 , neordr,  nrvf_ordr ,vlhxcQ
  use m_ES_LHXC,              only : m_ESlhxc_potential_3D
  use m_ES_Intgr_VlhxcQlm,    only : m_ESiVQ_integrate_VlhxcQlm_3D
  use m_Ionic_System,         only : zfm3_l
  use m_Parallelization,      only : myrank_k, map_k, init_zaj_para, myrank_e
  use m_PlaneWaveBasisSet,    only : ylm_l        &
  &                                , gr_l          &
  &                                , m_pwBS_sphrp_l_3D
  use m_PseudoPotential,      only : qitg_l ,  nqitg     &
 &                                 , psc_l               &
 &                                 , rhvg_l              &
 &                                 , modnrm, m_PP_gfqwei_3D, flg_paw, epc, m_PP_rd_PAW_parameters
  use m_XC_Potential,         only : m_XC_cal_potential_3D
  use m_Control_Parameters,   only : nspin, af, kimg, from_PSEUDOPOTENTIAL_FILE  &
 &                                 , istress,sw_fine_STM_simulation
  use m_Const_Parameters,     only : GAMMA
  use m_Kpoints,              only : kv3,vkxyz,k_symmetry
!F---
  use m_PlaneWaveBasisSet   ,only : iba
  use m_Control_Parameters  ,only : nblocksize_mgs,  nblocksize_mgs_is_given
  use m_Electronic_Structure,only : nblocksize_mgs_default,m_ES_alloc_zaj_l_prev
  use m_Parallelization     ,only : make_index_band_3D, make_index_band_for_Gdiv_3D, make_ball_buff
!F---
!!$#ifdef FJ_TIMER
!!$  use m_Parallelization,      only : MPI_CommGroup
!!$#endif
  use m_ES_ExactExchange,  only  : m_ES_EXX_gather_valence_states, m_ES_EXX_kernel &
       &                         , m_ES_EXX_occup,  m_ES_EXX_ngpt &
       &                         , m_ES_EXX_init &
       &                         , m_ES_EXX_update, m_ES_EXX_crotylm, m_ES_EXX_ylm
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================
  use m_Ionic_System,      only  : zeta1,m_IS_natm_can_change


   use m_Control_Parameters,     only : m_CtrlP_set_init_status, read_charge_hardpart, &
        &                               use_metagga, sw_read_pwbs_info, &
        &                               sw_calc_wf_orb_projection, sw_orb_popu, &
        &                               initial_ekin_dens
   use m_OP_rotation,            only : m_OP_init_Wfn_orb_proj
#ifdef MEMORY_SAVE_ZAJ_OLD
   use m_Control_Parameters,     only : RMM2P_is_specified
#endif
   use m_ErrorMessages
! ===================== KT_Test ========================= 12.5Exp
   use m_Control_Parameters,    only : m_CtrlP_set_hybrid_parameters, xctype, &
        &                              sw_opencore
   use m_ES_ExactExChange,      only : m_ES_EXX_set_nmax_G_hyb
   use m_PS_opencore,         only : m_PS_set_rmag_opencore, m_PS_set_pcc_opencore, &
        &                            m_PS_init_mag_opencore_pol, &
        &                            m_PS_rd_mag_opencore_pol
! ======================================================= 12.5Exp


! ================================== added by K. Tagami =============== 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor, ndim_magmom, &
       &                             import_collinear_spindensity, &
       &                             import_collinear_wavefunctions

!  use m_Charge_Density,  only : m_CD_initCD_by_file_rsp_noncl, &
!       &                        m_CD_adjust_spindensity_noncl, &
!       &                        m_CD_initial_CD_by_Gauss_kt, &
!       &                        m_CD_wd_chgq_l_portion_noncl
!  use m_ES_LHXC,         only : m_ESlhxc_potential_noncl

!  use m_Electronic_Structure,   only : m_ES_energy_eigen_vals_noncl, &
!       &                               m_ES_energy_eigenvals_ext_noncl, &
!       &                               m_ES_alloc_zaj_etc_noncl
!  use m_ES_Ortho,               only : m_ES_modified_gramschmidt_noncl

  use m_Files,                  only : m_Files_open_nfchr_noncl
!!  use m_PseudoPotential,        only : m_PP_set_fqwei_noncl
!!  use m_ES_NonCollinear,        only :  m_ES_init_Mat_dion0_noncl, &
!!       &                                m_ES_set_Mat_q_noncl
!
!!  use m_ES_IO,                   only : m_ESIO_rd_WFs_import_frm_collin
!!  use m_Charge_Density,          only : m_CD_rd_chgq_import_frm_collin, &
!!       &                                m_CD_rd_chgq_noncl
  use m_CD_Mag_Moment,           only : m_CD_alloc_RhoMag_on_atom, &
       &                                m_CD_calc_ChgMagMom_in_sphere, &
       &                                sw_monitor_atomcharge
  use m_ES_Mag_Constraint,       only : m_ES_add_MagConstraintPot_chgql, &
       &                                m_ES_add_MagConstraintPot_hsr
  use m_Orbital_Population,      only  : m_OP_rd_occ_mat_noncl
! ===================================================================== 11.0

! ================================== added by K. Tagami =============== 11.0
  use m_Control_Parameters,      only : SpinOrbit_mode
  use m_Const_Parameters,        only : Neglected, BuiltIn, ByProjector
  use m_SpinOrbit_Potential,     only : m_SO_alloc_dsoc, &
       &                                m_SO_set_Dsoc_potential1, &
       &                                m_SO_calc_MatLS_orb_s_to_f, &
       &                                m_SO_diagonalize_MatLS, &
       &                                m_SO_set_MatU_ylm_RC
  use m_ES_NonCollinear,         only : m_ES_alloc_factor_fss, &
       &                                 m_ES_set_factor_fss

! ==================================================================== 11.0

! ================================== added by K. Tagami =============== 11.0
  use m_Const_Parameters,        only : ByPawPot, ZeffApprox
  use m_SpinOrbit_Potential,    only : m_SO_alloc_Mat_SOC_strenth, &
       &                                m_SO_set_Dsoc_potential2
  use m_SpinOrbit_RadInt,       only : m_SO_calc_SOC_strength_pawpot, &
       &                               m_SO_calc_SOC_strength_zeff, &
       &                               m_SO_check_mode_Builtin, &
       &                               m_SO_check_mode_Pawpot, &
       &                               m_SO_check_mode_Zeff
! ======================================================================= 11.0
  use m_BerryPhase,             only : m_BP_rd_cntn_data

! ============= KT_add ==================== 13.0D
  use m_Files,  only : m_Files_reopen_nfcntn_bin_paw
! ========================================= 13.0D


  use m_Const_Parameters,    only : CRYSTAL_FIELD_APPROX
  use m_Orbital_Population,  only : m_OP_occ_mat_gen_kt

! ==== EXP_CELLOPT === 2015/09/24
  use m_IterationNumbers,     only : iteration_unit_cell,iteration_stress_correction
  use m_Control_Parameters,   only :  sw_read_nfchgt_prev_cell, sw_read_nfzaj_prev_cell
! ==================== 2015/09/24

  use m_SpinOrbit_SecondVariation, only : m_SO_init_second_variation
  use m_ES_dos,               only : m_ES_dos_alloc_compr_compi_ek

  use m_ES_nonlocal,            only : m_ES_AtaulmnaG
  use m_NonLocal_Potential,   only : m_NLP_alloc_taulmaG

  implicit none
  integer :: iloop
  logical, save :: is_charge_density_read = .false.
! <--
  integer :: neg_incre

  integer :: ierror
  logical :: read_pwbs

  integer :: ispin, ik, iksnl,ib
!!$#ifdef FJ_TIMER
!!$  integer :: ierr
!!$#endif

! -----> T. Yamasaki, 11 July 2008 ----
  real(kind=DP) :: zeta_sum
  logical :: initialization_required
  if(.not.initialization_required()) then
!     call m_CD_softpart_3D(nfout,kv3)
!     call m_CD_hardpart(nfout,kv3)
!     call m_CD_cp_chgq_to_chgqo()
     call m_PP_gfqwei_3D(nfout)  ! -> modnrm, fqwei, nlmta1, nlmta2
     call m_CtrlP_set_init_status(.false.)
     call UpdateUeff()
     if(sw_hybrid_functional == ON) then
       call m_ES_EXX_init()
       call EXX()
     endif
     return
  endif

  call UpdateUeff()

  read_pwbs = .false.
  if(sw_read_pwbs_info==ON)then
    call m_pwBS_rd_prev_pws(read_pwbs)
  endif

! ====================== KT_Test ========================= 12.5Exp
!  if(sw_hybrid_functional == ON) then
!     call m_CtrlP_set_hybrid_parameters
!     call m_ES_EXX_set_nmax_G_hyb
!  endif
! ======================================================== 12.5Exp

  if(sw_hybrid_functional == ON) call m_ES_EXX_update(ON)

  if(icond == INITIAL .or. icond==COORDINATE_CONTINUATION) then
     call m_ES_gtotch(nfout)
     if(icond == INITIAL) then
        if(totch >= neg*2.0) then
           if(ipri >= 1) then
              write(nfout,'(" ### Warning(1309): Number of bands(neg) is insufficient:")')
              write(nfout,'("                totch = ",f10.3," >= neg*2.0 = ",f10.3)') totch, neg*2.0
           end if
           if(dabs(sum(zeta1)) > 0.d0) then
              call m_CtrlP_set_neg_properly(1.3*totch) ! -> neg
           else
              call m_CtrlP_set_neg_properly(totch) ! -> neg
           end if
           call m_CntrlP_set_meg(neg)
           call m_CntrlP_set_neg(neg,nfout)
!           if(ipri >= 1) write(nfout,'(" neg (=num_bands) is enlarged ",i12)') neg
           if(ipri >= 1) write(nfout,'(" ### Warning(1309): Number of bands is enlarged",i12)') neg
           call m_CntrlP_set_davidson_size(neg) ! -> max_subspace_size
           call m_Parallel_dealloc_mpi_elec()


           call m_Parallel_init_mpi_elec_3D(nfout,ipriparallel,printable,neg,kv3,nspin,kg1,iba)
           call make_index_band_3D(nfout,ipriparallel,printable,kv3,neg &
                & , nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
           call make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
           call make_ball_buff
        end if
     end if
  else if(icond == CONTINUATION) then
     if(neg_previous > neg) then
        call m_CntrlP_set_neg(neg_previous,nfout)
        call m_CntrlP_set_meg(neg)
        call m_CntrlP_set_davidson_size(neg) ! -> max_subspace_size
        call m_Parallel_dealloc_mpi_elec()


        call m_Parallel_init_mpi_elec_3D(nfout,ipriparallel,printable,neg,kv3,nspin,kg1,iba)
        call make_index_band_3D(nfout,ipriparallel,printable,kv3,neg &
             &   , nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
        call make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
        call make_ball_buff
        if(ipri >= 1) write(nfout,'(" neg is enlarged ",f12.5)') neg
     end if
  end if
! <------

#ifdef _INCLUDE_EXX_
  if(sw_hybrid_functional == ON) call m_ES_EXX_init()
#endif

  if(nrank_e > neg) call phase_execution_error(PARALLELIZATION_INVALID_NE)

  if(.not.skip_alloc_phonon) then
     call m_CD_alloc_hsr()
     call m_ES_alloc_vlhxc()
     call m_ES_alloc_vlhxcQ()
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(30)
#endif

     call m_SO_set_MatU_ylm_RC

     if ( sw_monitor_atomcharge == ON ) then
        call m_CD_alloc_RhoMag_on_atom
     endif

     call m_ES_alloc_zaj_etc()

! === Use extrapolation! =======================================================
     call m_ES_wf_extrpl_alloc()
! ==============================================================================

#ifdef FJ_TIMER
                    call timer_end(30)
#endif
#ifdef MEMORY_SAVE_ZAJ_OLD
     if(RMM2P_is_specified) then
#endif
     call m_ESsd_alloc_zaj_old()
#ifdef MEMORY_SAVE_ZAJ_OLD
     end if
#endif
     call m_ESsd_alloc_dzajn2()

     if(sw_optimize_blocking_parameters==ON) call Optimize_Blocking_Parameters()

     if(sw_hubbard == ON) call m_ES_alloc_Dhub()
     if(sw_hybrid_functional == ON) call m_ES_EXX_ngpt()
     if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
        if(ekmode == ON) then
           call m_ES_alloc_eko_ek()
        else if(ekmode == OFF .and. fixed_charge_k_parallel == ONE_BY_ONE) then
           call m_ES_alloc_eko_ek()
        end if

        if(evaluation_eko_diff == ON) then
           call m_ES_alloc_eko1()
        else if(delta_eigenvalue_cond_is_given) then
           call m_ES_alloc_eko1()
           call m_ES_cpeko()
        end if
     else
        if(delta_eigenvalue_cond_is_given) then
           call m_ES_alloc_eko1()
           call m_ES_cpeko()
        end if
     end if
     if(sw_dipole_correction == ON) then
        call m_Dipole_vdip_alloc_3D()
     end if
  end if


  if(sw_phonon == ON .and. sw_calc_force == OFF) return

!!$  call m_PP_gfqwei(nfout)  ! -> modnrm, fqwei, nlmta1, nlmta2
!!$  if(modnrm == EXECUT .and. (icond == INITIAL .or. icond == CONTINUATION)) then
  if(.not.skip_alloc_phonon) then

     if ( sw_opencore == ON .and. nspin == 2 ) then
        call m_PS_set_rmag_opencore
        call m_PS_set_pcc_opencore
     endif


     call m_PP_gfqwei_3D(nfout)  ! -> modnrm, fqwei, nlmta1, nlmta2

     call m_NLP_alloc_taulmaG()
     if(modnrm == EXECUT ) then
        call m_pwBS_alloc_ylm_l()
        call m_pwBS_sphrp_l_3D()                           ! -> ylm_l
     end if
  end if

  if(sw_fef==ON) call m_FEF_Constract_of_ftq()

! ==== ASMS ====
  if ( sw_hybrid_functional == ON ) sw_eval_epc_on_fftmesh = ON
  if ( xctype == 'vdwdf' ) sw_eval_epc_on_fftmesh = ON
  if ( sw_opencore == ON .and. nspin == 2 ) sw_eval_epc_on_fftmesh = ON
  if ( flg_paw ) sw_eval_epc_on_fftmesh = off
  if ( sw_eval_epc_on_fftmesh == ON .and. sw_opencore == OFF ) then
     call m_XC_cal_potential_3D( nfout, Partial_Core_Charge, chgq_l, VXC_AND_EXC )
     epc = exc
  endif
! ==== ASMS ====

  if(iteration == 0) call m_TE_set_etotal_old()

  if(sw_screening_correction == ON) then
     call screening_potential
  end if
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
! ==============================================================================

  if(icond == INITIAL .or. icond == COORDINATE_CONTINUATION .or. &
       & (icond==FIXED_CHARGE .and. ekmode==OFF .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
     if(icond /= INITIAL) call m_ES_gtotch(nfout)
     if(icond == INITIAL) call check_neg()

     !---- set a charge density ----
     if(.not.is_charge_density_read.and.(initial_chg == FILE .or. (icond==FIXED_CHARGE.and.ekmode==OFF))) then
        call read_charge_density(condition=1)
        call read_efermi()
        if ( use_metagga ) then
           if ( initial_ekin_dens == FILE ) then
              call read_ekindens()
           else
!              call m_KE_init_ekindens
           endif
        endif
        if(icond==FIXED_CHARGE .and. ekmode == OFF) call copy_chgq_to_chgqo()
     else
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(30)
#endif
        if (iteration_unit_cell > 1 .and. sw_read_nfchgt_prev_cell == ON ) then
           call read_charge_density(condition=-4)
!        call m_CD_initial_CD_by_Gauss_func(nfout)   ! (intchg) -> chgq_l
#ifdef FFTW3
        else if ((iteration_unit_cell > 1 .or. iteration_stress_correction>1) .and. &
        &  kgp_prev .eq. kgp .and. sw_interpolate_charge == ON) then
           call m_CD_cp_chgq_prev_to_chgq(nfout)
        else if ((iteration_unit_cell > 1 .or. iteration_stress_correction>1) .and. &
        &  kgp_prev .ne. kgp .and. sw_interpolate_charge == ON) then
           call m_CD_gen_chgq_from_prev_chgq(nfout)
#endif
        else
          call m_CD_initial_CD_by_Gauss_kt(nfout)   ! (intchg) -> chgq_l
        endif
#ifdef FJ_TIMER
                    call timer_end(30)
#endif
!        if ( use_metagga ) call m_KE_init_ekindens
     end if

     if ( sw_opencore == ON ) then
        call m_PS_init_mag_opencore_pol
        if ( sw_eval_epc_on_fftmesh == ON ) then
           call m_XC_cal_potential_3D( nfout, Partial_Core_Charge, chgq_l, EXC_ONLY )
           epc = exc
        endif
     endif

     !---- set wave functions ----
     if ( ((iteration_unit_cell > 1 .or. iteration_stress_correction>1) .and.  imdalg .ne. PT_CONTROL) .and. &
!     &  (sw_read_nfzaj_prev_cell == ON .or. is_charge_density_read)) then
     &  sw_read_nfzaj_prev_cell == ON ) then
        call read_zaj( condition =-4 )
     else if((iteration_unit_cell > 1 .or. iteration_stress_correction>1) .and. &
     &  kg1_prev .eq. kg1 .and. sw_interpolate_wfs == ON) then
        call m_ES_cp_zaj_prev_to_zaj()
     else if((iteration_unit_cell > 1 .or. iteration_stress_correction>1) .and. &
     &  kg1_prev .ne. kg1 .and. sw_interpolate_wfs == ON) then
#ifdef FFTW3
        call m_ES_gen_zaj_from_prev_zaj(nfout)
#else
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(40)
#endif
        if ( init_zaj_para ) then
           call m_ESIW_by_randomnumbers_3D(nfout,kv3,1,neg)      ! (rndzaj) -> zaj_l
        else
           call m_ESIW_by_randomnumbers0_3D(nfout,kv3,1,neg)      ! (rndzaj) -> zaj_l
        endif
#ifdef FJ_TIMER
                    call timer_end(40)
#endif
#endif
     else if(intzaj == by_random_numbers) then
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(40)
#endif
        if ( init_zaj_para ) then
           call m_ESIW_by_randomnumbers_3D(nfout,kv3,1,neg)      ! (rndzaj) -> zaj_l
        else
           call m_ESIW_by_randomnumbers0_3D(nfout,kv3,1,neg)      ! (rndzaj) -> zaj_l
        endif
#ifdef FJ_TIMER
                    call timer_end(40)
#endif
     else if(intzaj == by_pseudo_atomic_orbitals) then
        call m_ESIW_by_atomic_orbitals(nfout,kv3,1,neg)    ! (paozaj) -> zaj_l
     else if(intzaj == by_matrix_diagon) then

     else if(intzaj == FILE) then
        call read_zaj( condition = 1 )
! === ik is not defined here!!! by tkato 2013/02/12 ============================
!!$       call m_ES_betar_dot_WFs_3D(nfout,ik)         ! (fsrfsi)
       do ik = 1, kv3, af+1
          if(map_k(ik) /= myrank_k) cycle
          call m_ES_betar_dot_WFs_3D(nfout,ik)         ! (fsrfsi)
       end do

        if(initial_chg == FILE .or. (icond==FIXED_CHARGE.and.ekmode==OFF)) then
        else
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(34)
#endif
           call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge,chgq_l, VXC_AND_EXC)  ! -> vxc_l
#ifdef FJ_TIMER
                    call timer_end(34)
#endif
           call Local_Hartree_XC_potential()
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(36)
#endif
           call m_ESiVQ_integrate_VlhxcQlm_3D(nfout) ! (lclchh) -> vlhxcQ
           if(sw_precalculate_phase_vnonlocal==ON) then
             call m_ES_AtaulmnaG(hardpart=.true.)
           endif
#ifdef FJ_TIMER
                    call timer_end(36)
#endif
           call energy_eigen_values()
           call ChargeDensity_Construction(1)
        end if
     endif

     !---- set a potential ----
     call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge,chgq_l, VXC_AND_EXC) ! -> vxc_l
     call Local_Hartree_XC_potential()
     if(intzaj == by_random_numbers .or. intzaj == by_pseudo_atomic_orbitals .or. intzaj == FILE) then
        if ( sw_monitor_atomcharge == ON ) then
           call m_CD_calc_ChgMagMom_in_sphere
        endif
        if ( sw_magnetic_constraint == ON ) then
           call m_ES_add_MagConstraintPot_chgql
           call m_ES_add_MagConstraintPot_hsr
        endif
     end if

! === DEBUG by tkato 2014/08/14 ================================================
!    call modified_gramschmidt()
     if(intzaj /= by_matrix_diagon) call modified_gramschmidt()
! ==============================================================================
!     if(intzaj /= by_matrix_diagon) then
      if(sw_hybrid_functional == ON) call EXX()
!     end if
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(36)
#endif
     call m_ESiVQ_integrate_VlhxcQlm_3D(nfout)   ! (lclchh) -> vlhxcQ
     if(sw_precalculate_phase_vnonlocal==ON) then
       call m_ES_AtaulmnaG(hardpart=.true.)
     endif
#ifdef FJ_TIMER
                    call timer_end(36)
#endif
     if(intzaj /= by_matrix_diagon) call energy_eigen_values()

     if(iprichargedensity >= 2) call m_CD_wd_chgq_l_small_portion(nfout)
     if(icond==COORDINATE_CONTINUATION) icond = INITIAL

     if(sw_hubbard == ON) then
        if(intzaj==FILE)then
           if(initial_occmat == FILE .or. icond==FIXED_CHARGE .or.  iteration_unit_cell > 1) then
              call read_occ_mat() ! -> om
           else
              call ChargeDensity_Construction(1)
              call Renewal_of_OccMat(.false.,ON, .false.)
           endif
        else
           if(initial_occmat == FILE .or. icond==FIXED_CHARGE .or.  iteration_unit_cell > 1) then
              call read_occ_mat() ! -> om
           else if (sw_initial_es==ON)then
              call ChargeDensity_Construction(1)
              call Renewal_of_OccMat(.false.,ON, .false.)
           else
              if ( initial_occmat == CRYSTAL_FIELD_APPROX ) then
                 call m_OP_occ_mat_gen_kt(nfout) ! -> om
              else
                 call m_OP_occ_mat_init(nfout) ! -> om
              endif
           end if
        endif
        call m_OP_mix_om(1.d0) ! om -> ommix
        call Renewal_of_Hubbard_Potential() ! ommix -> dhub
     end if

  else if(icond == CONTINUATION .or. &
       & (icond==FIXED_CHARGE_CONTINUATION .and. ekmode==OFF &
       & .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
     call m_Files_open_nfcntn_bin
     call m_ESIO_rd_EigenValues_etc(nfout,nfcntn_bin,F_CNTN_BIN_in_partitioned)
!!$     call m_ESIO_wd_EigenValues(nfout,2,nooccupation=NO)
     call m_Files_open_nfzaj()
     if(read_pwbs)then
       call m_ES_alloc_zaj_l_prev(np_g1k_prev)
     endif
     call m_ESIO_rd_WFs(nfout,nfzaj,F_ZAJ_in_partitioned,read_pwbs)
#ifdef FFTW3
     if(read_pwbs)then
       call m_ES_gen_zaj_from_prev_zaj(nfout)
       call modified_gramschmidt()
     endif
#endif
     if(neg_previous < neg) then
        if ( init_zaj_para ) then
           call m_ESIW_by_randomnumbers_3D(nfout,kv3,neg_previous+1,neg)      ! (rndzaj) -> zaj_l
        else
           call m_ESIW_by_randomnumbers0_3D(nfout,kv3,neg_previous+1,neg)      ! (rndzaj) -> zaj_l
        endif
     end if
     call read_charge_density(condition = 2)
     if(flg_paw) call m_CD_rd_hsr(nfcntn_bin_paw)
     if ( use_metagga ) call read_ekindens()
     if ( sw_opencore == ON ) then
        call m_PS_rd_mag_opencore_pol
        if ( sw_eval_epc_on_fftmesh == ON ) then
           call m_XC_cal_potential_3D( nfout, Partial_Core_Charge, chgq_l, EXC_ONLY )
           epc = exc
        endif
     endif

     if ( sw_monitor_atomcharge == ON ) then
        call m_CD_calc_ChgMagMom_in_sphere
     endif

     if(icond /= CONTINUATION) call read_efermi()

     do ik = 1, kv3, af+1
        if(map_k(ik) /= myrank_k) cycle         ! MPI
        call m_ES_betar_dot_WFs_3D(nfout,ik)            ! (fsrfsi) ->fsr_l
     end do
     if(sw_hybrid_functional == ON) call EXX()

     if(neg_previous < neg) then
        call modified_gramschmidt()
        call m_ES_energy_eigen_values_ext_3D(nfout)
     end if
     call Renewal_of_Potential()
     if(sw_hubbard == ON) then
        call read_occ_mat()    ! -> om
        call m_OP_mix_om(1.d0) ! om -> ommix
        call Renewal_of_Hubbard_Potential() ! ommix -> dhub
     end if
     if(icond==FIXED_CHARGE_CONTINUATION .and. ekmode == OFF) then
        call copy_chgq_to_chgqo()
        if(sw_hubbard == ON) call m_OP_cp_ommix_to_omold()
     end if

  else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
     call m_ES_gtotch(nfout)

     call read_charge_density(condition = 2)
     if(flg_paw.and.ekmode/=OFF) call m_CD_rd_hsr(nfcntn_bin_paw)
     if ( use_metagga ) call read_ekindens()
     if ( sw_opencore == ON ) then
        call m_PS_rd_mag_opencore_pol
        if ( sw_eval_epc_on_fftmesh == ON ) then
           call m_XC_cal_potential_3D( nfout, Partial_Core_Charge, chgq_l, EXC_ONLY )
           epc = exc
        endif
     endif

     if ( sw_monitor_atomcharge == ON ) then
        call m_CD_calc_ChgMagMom_in_sphere
     endif

     call read_efermi()

     if(sw_hybrid_functional == ON) call EXX()
     call Renewal_of_Potential()
     if(sw_hubbard==ON) then
        call read_occ_mat() ! -> om
        call m_OP_mix_om(1.d0) ! om -> ommix
        call Renewal_of_Hubbard_Potential() ! ommix -> dhub
     end if
     if(icond == FIXED_CHARGE_CONTINUATION .and. nk_in_the_process >= 2) call m_Files_skiptoend(nfeng)
     if(icond == FIXED_CHARGE_CONTINUATION .and. sw_berry_phase == ON ) call m_BP_rd_cntn_data()
     if(sw_pdos == ON) then
        call m_ES_dos_alloc_compr_compi_ek()
     endif
  endif

  if (icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
     if ( sw_orb_popu == ON .and. sw_calc_wf_orb_projection == ON ) then
        call m_OP_init_Wfn_orb_proj(nfout)
     endif
     if ( .not. noncol .and. sw_spinorbit_second_variation == ON ) then
        if ( ekmode == ON ) call m_SO_init_second_variation
     endif
  endif

  call m_CtrlP_set_init_status(.false.)


contains
  subroutine read_charge_density(condition)
    use m_Files,              only : m_Files_open_nfchr
    use m_Charge_Density,     only : m_CD_rd_chgq &
         &                         , m_CD_initial_CD_by_file_rspace &
         &                         , m_CD_adjust_spindensity
    integer, intent(in) :: condition
    integer :: ispin
    integer :: iloop

    call m_Files_open_nfchgt()
!!$    write(6,'( " icond = ",i8, " condition = ", i8,"  (read_charge_density)")') icond,condition
    is_charge_density_read = .true.
    if(condition == 1) then
       if(sw_initial_charge_rspace == ON) then
             do iloop = 1, nspin
                call m_Files_open_nfchr(nspin,iloop)
                call m_CD_initial_CD_by_file_rspace(nspin,iloop,nfout,nfchr)
             end do
       else
             call m_CD_rd_chgq(nfout,nfchgt,F_CHGT_in_partitioned,read_pwbs)
          if ( flg_paw .and. read_charge_hardpart == YES .and. m_Files_nfcntn_bin_paw_exists()) then
             call m_Files_open_nfcntn_bin_paw()
             call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
             call m_CD_rd_hsr(nfcntn_bin_paw)
             call m_Files_close_nfcntn_bin_paw()
          endif
       end if

          call m_CD_adjust_spindensity(nfout)

    else if(condition == 2) then
          call m_CD_rd_chgq(nfout,nfchgt,F_CHGT_in_partitioned,read_pwbs)
!!$       if ( flg_paw .and. ) then
!!$          call m_CD_rd_hsr(nfcntn_bin_paw)
!!$       endif

    else if(condition == -4) then      ! coordinate-continuation
!
    end if

  end subroutine read_charge_density

  subroutine read_ekindens
    use m_Files,  only :  m_Files_open_nf_ekindens, m_Files_close_nf_ekindens
    use m_KineticEnergy_Density, only : m_KE_rd_ekindens

    call m_Files_open_nf_ekindens
    call m_KE_rd_ekindens()
    call m_Files_close_nf_ekindens

  end subroutine read_ekindens

  subroutine read_efermi()
    use m_Files,              only : nfout,m_Files_open_nfefermi, m_Files_close_nfefermi
    use m_ES_IO,              only : m_ESIO_rd_Efermi

    call m_Files_open_nfefermi()
    call m_ESIO_rd_Efermi(nfout,nfefermi)
    call m_Files_close_nfefermi()
  end subroutine read_efermi

  subroutine copy_chgq_to_chgqo()
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
        call timer_sta(30)
#endif
    call m_CD_cp_chgq_to_chgqo()
#ifdef FJ_TIMER
        call timer_end(30)
#endif
  end subroutine copy_chgq_to_chgqo

  subroutine Local_Hartree_XC_potential()
    use m_ES_LHXC,            only : m_ESlhxc_potential_3D
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(35)
#endif
    call m_ESlhxc_potential_3D(nfout,chgq_l,vxc_l) ! (stlhxc) ->vlhxc_l
#ifdef FJ_TIMER
                    call timer_end(35)
#endif
  end subroutine Local_Hartree_XC_potential

  subroutine energy_eigen_values()
    use m_Electronic_Structure,   only : m_ES_energy_eigen_values_3D
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(33)
#endif
    call m_ES_energy_eigen_values_3D(nfout)   ! (eigen0) -> eko_l,neordr
#ifdef FJ_TIMER
                    call timer_end(33)
#endif
  end subroutine energy_eigen_values

  subroutine modified_gramschmidt()
    use m_ES_ortho,           only : m_ES_modified_gram_schmidt
#ifdef FJ_TIMER
!                    call mpi_barrier(MPI_CommGroup, ierr)
                    call timer_sta(32)
#endif
       call m_ES_modified_gram_schmidt(nfout) ! ->zaj_l, fsr_l,fsi_l
#ifdef FJ_TIMER
                    call timer_end(32)
#endif
  end subroutine modified_gramschmidt

  subroutine EXX()
    use m_ES_ExactExchange,  only  : m_ES_EXX_gather_valence_states, m_ES_EXX_kernel &
         &                         , m_ES_EXX_occup  &
         &                         , sw_rspace_hyb, m_ES_EXX_init0 &
         &                         , m_ES_EXX_ylm, m_ES_EXX_crotylm
    call m_ES_EXX_init0()
    if(modnrm == EXECUT ) then
       if(sw_rspace_hyb==OFF) call m_ES_EXX_ylm()
       call m_ES_EXX_crotylm()
    end if
    call m_ES_EXX_occup(nfout)
    call m_ES_EXX_gather_valence_states(nfout)
    call m_ES_EXX_kernel(nfout)
  end subroutine EXX

  subroutine read_occ_mat()
    use m_Orbital_Population,      only  : m_OP_rd_occ_mat
       call m_OP_rd_occ_mat(nfout) ! -> om
  end subroutine read_occ_mat

  subroutine read_zaj( condition )
    use m_ES_IO,                   only : m_ESIO_import_WFs_prev_cell
    use m_Parallelization,         only : np_g1k_prev
    use m_Electronic_Structure,    only : m_ES_alloc_zaj_l_prev
    use m_PlaneWaveBasisSet,       only : kg1_prev
    integer, intent(in) :: condition

    call m_Files_open_nfzaj()

    if ( condition == 1 ) then
          if(read_pwbs)then
            call m_ES_alloc_zaj_l_prev(np_g1k_prev)
          endif
          call m_ESIO_rd_WFs(nfout,nfzaj,F_ZAJ_in_partitioned,read_pwbs)
#ifdef FFTW3
          if(read_pwbs)then
            call m_ES_gen_zaj_from_prev_zaj(nfout)
            call modified_gramschmidt()
          endif
#endif

    else if ( condition == -4 ) then       ! coordinate-continuation
       call m_ESIO_import_WFs_prev_cell(nfout,nfzaj,F_ZAJ_in_partitioned)
    endif

  end subroutine read_zaj

  subroutine check_neg()
    if(totch >= neg*2.0) then
       if(ipri >= 1) write(nfout,'(" This job stops because totch = ",f12.5, " > neg*2.0 = ",f12.5)') &
            & totch, neg*2.0
       call phase_error_with_msg(nfout,' totch > neg*2.0',__LINE__,__FILE__)
    end if
  end subroutine check_neg
end subroutine Initial_Electronic_Structure

subroutine Initial_WaveFunctions_ek
!                           @(#)Initial_Electronic_Structure.f90 1.4 03/02/19 00:43:18
  use m_Const_Parameters,   only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                         , by_random_numbers, by_matrix_diagon, YES, ON, OFF &
       &                         , by_pseudo_atomic_orbitals
  use m_IterationNumbers,   only : nk_in_the_process &
       &                         , first_kpoint_in_this_job &
       &                         , iteration_electronic &
       &                         , m_Iter_reset_iter_electronic
  use m_Control_Parameters, only : ekmode,icond,intzaj,evaluation_eko_diff &
       &                         , ipri,ipriekzaj, sw_ekzaj,neg, printable, numk_zajsaved &
       &                         , sw_hybrid_functional, kimg
  use m_Files,              only : nfout,nfzaj, m_Files_open_nfzaj, m_Files_open_nfzaj_append
  use m_Kpoints,            only : kv3, kv3_ek, vkxyz, vkxyz_ek
  use m_Electronic_Structure,only: m_ES_energy_eigen_values_3D &
       &                         , m_ES_cpeko &
       &                         , m_ES_wd_zaj_small_portion0
  use m_ES_nonlocal,        only : m_ES_betar_dot_WFs_3D
  use m_ES_ortho,           only : m_ES_modified_gram_schmidt
  use m_ES_initialWF,       only : m_ESIW_by_randomnumbers_3D &
       &                         , m_ESIW_by_atomic_orbitals
  use m_ES_IO,              only : m_ESIO_rd_WFs_and_EVs_ek
  use m_ES_ExactExchange,  only  : m_ES_EXX_gather_valence_states, m_ES_EXX_kernel &
       &                         , m_ES_EXX_occup

! =========================================== added by K. Tagami ========== 11.0
  use m_Control_Parameters,   only : noncol, ndim_magmom
  use m_Control_Parameters,   only : af
  use m_Parallelization,      only : myrank_k, map_k
! ========================================================================= 11.0

  implicit none
  integer :: ik

  if(printable) write(nfout,'("!! icond, nk_in_the_process, first_kpoint_in_this_job = ",3i8)') &
       & icond, nk_in_the_process, first_kpoint_in_this_job
  if(ekmode == ON) then
     if(icond == FIXED_CHARGE &
          & .or. (icond == FIXED_CHARGE_CONTINUATION &
          &         .and. nk_in_the_process /= first_kpoint_in_this_job)) then
        if(sw_ekzaj == OFF) then
           if(intzaj == by_random_numbers .or. intzaj == by_pseudo_atomic_orbitals) then
              if(nk_in_the_process==1)then
                 if(intzaj == by_random_numbers) then
                    call m_ESIW_by_randomnumbers_3D(nfout,kv3,1,neg)! (rndzaj) -> zaj_l
                 else
                    call m_ESIW_by_atomic_orbitals(nfout,kv3,1,neg)! (paozaj) -> zaj_l
                 end if
              endif
              do ik = 1, kv3, af+1
                 if(map_k(ik) /= myrank_k) cycle
                 call m_ES_betar_dot_WFs_3D(nfout,ik)     ! (fsrfsi,sumset)
              end do

! ======================================= modified by K. Tagami ========== 11.0
!              call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
!
              call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
! ======================================================================== 11.0

              if(sw_hybrid_functional == ON) then
                 call m_ES_EXX_occup(nfout)
                 call m_ES_EXX_gather_valence_states(nfout)
                 call m_ES_EXX_kernel(nfout)
              end if

! ======================================= modified by K. Tagami ========== 11.0
!              call m_ES_energy_eigen_values(nfout)   ! (eigen0) -> eko_l,neordr
!
              call m_ES_energy_eigen_values_3D(nfout)   ! (eigen0) -> eko_l,neordr
! ======================================================================= 11.0

              if(intzaj == by_matrix_diagon) then

              end if
           end if
           if(icond == FIXED_CHARGE_CONTINUATION.and.iteration_electronic == 0) &
                & call m_Files_open_nfzaj_append()
        else
           call m_Files_open_nfzaj()
           call m_ESIO_rd_WFs_and_EVs_ek(nfout,nfzaj) ! -> zaj_l,eko_l,neordr
           do ik = 1, kv3, af+1
              if(map_k(ik) /= myrank_k) cycle
              call m_ES_betar_dot_WFs_3D(nfout,ik)     ! (fsrfsi,sumset)
           end do
        end if
     else if (icond == FIXED_CHARGE_CONTINUATION &
          & .and. nk_in_the_process == first_kpoint_in_this_job) then
        call m_Files_open_nfzaj()
        call m_ESIO_rd_WFs_and_EVs_ek(nfout,nfzaj)
        do ik = 1, kv3, af+1
           if(map_k(ik) /= myrank_k) cycle
           call m_ES_betar_dot_WFs_3D(nfout,ik)
        end do
     else
        call phase_error_with_msg(nfout,' icond is illegal (Initial_Electronic_Structure_ek)', &
        __LINE__,__FILE__)
     end if
  else if(ekmode == OFF) then
     if(ipri>=1) write(nfout,'(" kv3 = ",i8)') kv3
     if(icond == FIXED_CHARGE_CONTINUATION .and. first_kpoint_in_this_job == nk_in_the_process) then
        if(ipri>=1) write(nfout,'("!! icond = ",i8)') icond
        call m_Files_open_nfzaj()
        call m_ESIO_rd_WFs_and_EVs_ek(nfout,nfzaj)
        if(nk_in_the_process+kv3-1 > numk_zajsaved) then
           call m_ESIW_by_randomnumbers_3D(nfout,kv3,1,neg)
           do ik = 1, kv3, af+1
              if(map_k(ik) /= myrank_k) cycle
              call m_ES_betar_dot_WFs_3D(nfout,ik)     ! (fsrfsi,sumset)
           end do

! ================================ modified by K. Tagami ==================== 11.0
!          call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
           call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
! ========================================================================== 11.0

           if(sw_hybrid_functional == ON) then
              call m_ES_EXX_occup(nfout)
              call m_ES_EXX_gather_valence_states(nfout)
              call m_ES_EXX_kernel(nfout)
           end if

! ================================= modified by K. Tagami =============== 11.0
!           call m_ES_energy_eigen_values(nfout)   ! (eigen0) -> eko_l,neordr
          call m_ES_energy_eigen_values_3D(nfout)   ! (eigen0) -> eko_l,neordr
! ====================================================================== 11.0

        else
           do ik = 1, kv3, af+1
              if(map_k(ik) /= myrank_k) cycle
              call m_ES_betar_dot_WFs_3D(nfout,ik)
           end do
        end if
        call m_ES_wd_zaj_small_portion0(" -- Initial_WF_ek --",20)
     else
        if(icond==FIXED_CHARGE) then
           call m_Files_open_nfzaj()
        else
           call m_Files_open_nfzaj_append()
        end if
        if(intzaj == by_random_numbers .or. intzaj == by_pseudo_atomic_orbitals) then
           if(intzaj == by_random_numbers) then
              call m_ESIW_by_randomnumbers_3D(nfout,kv3,1,neg)! (rndzaj) -> zaj_l
           else
              call m_ESIW_by_atomic_orbitals(nfout,kv3,1,neg)! (paozaj) -> zaj_l
           end if
           do ik = 1, kv3, af+1
              if(map_k(ik) /= myrank_k) cycle
              call m_ES_betar_dot_WFs_3D(nfout,ik)     ! (fsrfsi,sumset)
           end do

! =========================================== modified by K. Tagami =========== 11.0
!           call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
           call m_ES_modified_gram_schmidt(nfout) ! (grmsmd)
! ============================================================================== 11.0

           if(sw_hybrid_functional == ON) then
              call m_ES_EXX_occup(nfout)
              call m_ES_EXX_gather_valence_states(nfout)
              call m_ES_EXX_kernel(nfout)
           end if

! ================================ modifie by K. Tagami ==================== 11.0
!           call m_ES_energy_eigen_values(nfout)   ! (eigen0) -> eko_l,neordr
           call m_ES_energy_eigen_values_3D(nfout)   ! (eigen0) -> eko_l,neordr
! ======================================================================== 11.0

        end if
        if(intzaj == by_matrix_diagon .or. ipriekzaj == 0) &
             & call m_Iter_reset_iter_electronic(nfout)
     end if
  end if
  if(evaluation_eko_diff == YES) call m_ES_cpeko()
end subroutine Initial_WaveFunctions_ek


