!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: ChargeDensity_Mixing
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
subroutine WriteDownData_onto_Files(final_call)
! $Id: WriteDownData_onto_Files.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Files,               only : m_Files_open_nfcntn,    m_Files_close_all &
       &                          , m_Files_open_nfstatus,  m_Files_skiptoend &
       &                          , m_Files_close_nfinp &
       &                          , m_Files_reopen_nfcntn,  m_Files_reopen_nfcntn_bin &
       &                          , m_Files_close_nfcntn &
       &                          , m_Files_open_nfvlc  &
       &                          , m_Files_open_nfcntn_bin_paw, m_Files_close_nfcntn_bin_paw &
       &                          , m_Files_reopen_nfzaj,   m_Files_reopen_nfchgt &
       &                          , m_Files_open_nfefermi,  m_Files_close_nfefermi &
       &                          , m_Files_open_nfPsicoef, m_Files_close_nfPsicoef &
       &                          , m_Files_open_nfBandSymInput, m_Files_close_nfBandSymInput &
       &                          , m_Files_open_nfeng,     m_Files_close_nfeng &
       &                          , m_Files_checkpoint_dir, m_Files_close_logfile &
       &                          , nfstatus, F_ZAJ_partitioned, F_CHGT_partitioned &
       &                          , F_ZAJ_in_partitioned, F_CHGT_in_partitioned &
       &                          , F_CNTN_BIN_partitioned, F_CNTN_BIN_in_partitioned &
#ifdef _EMPIRICAL_
       &                     ,nfcntn,nfout
#else
       &                     ,nfcntn,nfchgt,nfzaj,nfpsicoef,nfbandsyminput,nfeng &
       &                     ,nfcntn_bin,nfout,nfvlc,nfefermi  &
       &                     ,nfcntn_bin_paw, F_ZAJ
#endif
  use m_Const_Parameters,    only: T_CONTROL, BLUEMOON, QUENCHED_CONSTRAINT &
       &                         , ON,OFF,YES,NO,FORCE_CONVERGED, Valence_plus_PC_Charge &
       &                         , VXC_AND_EXC, FINISH, INITIAL, EK_CONVERGED &
       &                         , CONTINUATION &
       &                         , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, ONE_BY_ONE, EK &
       &                         , DRIVER_CONSTRAINT, BFGS, GDIIS, ALL_AT_ONCE, GRID, P_CONTROL, PT_CONTROL &
       &                         , DRIVER_NEB
  use m_IterationNumbers,    only: m_Iter_wd_iteration_numbers &
       &                         , iteration, iteration_ionic, iteration_electronic
  use m_Timing,              only: tstatc_wd0, m_Timing_wd_status, m_Timing_wd_timenow
  use m_Ionic_System,        only: m_IS_wd_pos_and_v, m_IS_wd_forcp_etc &
       &                         , m_IS_wd_nrsv,      m_IS_wd_cps, m_IS_wd_diis_history &
       &                         , m_IS_natm_can_change, m_IS_wd_curr_atom_reservoir &
       &                         , m_IS_wd_pos_brav &
       &                         , m_IS_wd_rigid_body &
       &                         ,  m_IS_wd_moved_distances_of_planes
  use m_Total_Energy,        only: m_TE_wd_total_energy
  use m_Control_Parameters, only : imdalg,iconvergence_previous_job &
       &                         , continuation_using_ppdata &
       &                         , driver &
       &                         , m_CtrlP_wd_iconvergence &
       &                         , m_CtrlP_wd_edelta_ontheway &
       &                         , m_CtrlP_set_corecharge_cntnbin &
       &                         , m_CtrlP_wd_corecharge_cntnbin &
#ifndef _EMPIRICAL_
       &                         , sw_ekzaj, sw_hubbard, sw_positron &
       &                         , m_CtrlP_wd_neg &
       &                         , sw_rttddft,iconvergence &
#endif
       &                         , iprijobstatus, jobstatus_series, jobstatus_format &
       &                         , icond, sw_write_zaj &
       &                         , sw_optimize_lattice, ipricoefwf &
       &                         , sw_band_symmetry_analysis, iprimagmom &
       &                         , ipritiming &
       &                         , sw_output_hybrid_info
#ifndef _EMPIRICAL_
  use m_IterationNumbers,    only: m_Iter_wd_iters_and_nk
  use m_Control_Parameters,  only: sw_fine_STM_simulation, fixed_charge_k_parallel &
       &                         , m_CtrlP_wd_isolver &
       &                         , m_CtrlP_wd_iconv_ek
  use m_Kpoints,             only: kv3, kv3_ek, m_Kp_wd_kv3, m_Kp_wd_BandSymInput
  use m_PseudoPotential,     only: flg_paw,m_PP_wd_PAW_parameters
  use m_NonLocal_Potential,  only: m_NLP_wd_snl_3D
  use m_PseudoPotential,     only: m_PP_wd_PP_parameters_3D
  use m_ES_IO,               only: m_ESIO_wd_EigenValues, m_ESIO_wd_EigenValues_etc &
       &                         , m_ESIO_wd_Efermi
  use m_ES_IO,               only: m_ESIO_wd_WFs_3D
!!$       &                         , m_ESIO_wd_WFs_dp_3D &  ! write file DP test
  use m_Electronic_Structure,only: iconv_ek, vbm, efermi, metalic_system, check_if_metalic_flag
  use m_Charge_Density,      only: m_CD_wd_hsr
  use m_Charge_Density,      only: m_CD_wd_chgq
  use m_Orbital_Population,  only: m_OP_wd_occ_mat
  use m_PlaneWaveBasisSet,   only: kgp, m_pwBS_wd_curr_pws
  use m_PAW_ChargeDensity,   only: m_PAWCD_wd_cd
  use m_Crystal_Structure,   only: m_CS_wd_fix_spin_status, m_CS_wd_BandSymInput
#endif
  use m_Parallelization,     only: mype, npes, conf_para,  MPI_CommGroup
  use m_Parallelization,     only: m_Parallel_wd_npes_etc_3D
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
  use m_Control_Parameters,  only: m_CtrlP_wd_dtim_previous
! ==============================================================================

  use m_Control_Parameters,   only : noncol, ndim_spinor, num_extra_bands, neg, nspin, ekmode &
       &                           , neg_is_enlarged, sw_write_pwbs_info
  use m_UnitCell,            only: m_UC_wd_cntn_data, m_UC_wd_md_cntn_data

  use m_ES_occup,             only: m_ESoc_check_num_bands, m_ESoc_check_if_metalic
! === Restart with phase ek-mode is supported. by tkato 2014/01/23 =============
  use m_Control_Parameters,  only : m_CtrlP_wd_numk_zajsaved
  use m_IterationNumbers,    only : nk_in_the_process
! ==============================================================================
! ======================================= added by K. Tagami ============ 11.0
  use m_Control_Parameters,   only :  noncol
!  use m_Charge_Density,       only : m_CD_wd_chgq_noncl
  use m_Orbital_Population,   only: m_OP_wd_occ_mat_noncl
! ======================================================================== 11.0

  use m_ES_ExactExchange,     only : m_ES_EXX_wd_hybridinfo
  use m_Control_Parameters,  only : use_metagga, sw_opencore
  use m_Files,  only :  m_Files_open_nf_ekindens, m_Files_close_nf_ekindens
  use m_KineticEnergy_Density, only : m_KE_wd_ekindens
  use m_PS_opencore,    only : m_PS_wd_mag_opencore_pol

  use m_Control_Parameters, only : sw_stress_correction
  use m_Stress, only : m_Stress_wd_cdata_4_correction
  use m_CS_SpaceGroup, only : m_CS_SG_wd_cntn
  use mpi

  implicit none
  logical,intent(in),optional :: final_call
  logical :: final_c = .true.
  integer :: status_wdmode
  logical :: flg_wd_zaj = .true.

  integer :: mpi_err
  integer :: sw_wd_Psi_coef = OFF
! === Restart with phase ek-mode is supported. by tkato 2014/01/23 =============
  integer :: numk_zajsaved
! ==============================================================================

  logical :: isSCDFT

#ifndef _EMPIRICAL_
  logical :: ChargeDensity_is_Converged
#endif

!  include 'mpif.h'                                      ! MPI

! === An IF BLOCK STRUCTURE ===
!=  if((icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) &
!=       & .and. fixed_charge_k_parallel == ONE_BY_ONE) then
!=     if(iconvergence_previous_job < EK_CONVERGED) then
!=        ....
!=     else
!=         write(nfout,'(" --- iconvergence_previous >= EK_CONVERGED ---")')
!=     end if
!=  else
!=     if(iconvergence_previous_job < FORCE_CONVERGED .or. sw_rttddft == ON) then
!=        ....
!=     else
!=        if((F_ZAJ_partitioned .neqv. F_ZAJ_in_partitioned) .or. icond == INITIAL) then
!=           ....
!=        end if
!=        if((F_CHGT_partitioned .neqv. F_CHGT_in_partitioned) .or. icond == INITIAL) then
!=           ....
!=        end if
!=        ....
!=     end if
!=  end if

  flg_wd_zaj = .true.
  if ( F_ZAJ == "/dev/null" .or. sw_write_zaj == OFF ) flg_wd_zaj = .false.

  if (present(final_call)) final_c = final_call
  call m_Files_checkpoint_dir(driver /= DRIVER_NEB)
#ifndef _EMPIRICAL_
  if((icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) &
       & .and. fixed_charge_k_parallel == ONE_BY_ONE) then
     if(iconvergence_previous_job < EK_CONVERGED) then
        call m_Files_open_nfcntn
        if(mype == 0) rewind nfcntn
        call m_Iter_wd_iters_and_nk(nfcntn)
        call m_Parallel_wd_npes_etc_3D(nfcntn)
        call m_Kp_wd_kv3(nfcntn)
        call m_CtrlP_wd_iconvergence(nfcntn)
        call m_CtrlP_wd_iconv_ek(kv3_ek,iconv_ek,nfcntn)
     else
        write(nfout,'(" --- iconvergence_previous >= EK_CONVERGED ---")')
     end if
  else
#endif
!  if(iconvergence_previous_job < FORCE_CONVERGED) then
  if(iconvergence_previous_job < FORCE_CONVERGED .or. sw_rttddft == ON) then
     if(sw_rttddft == ON .and. iconvergence_previous_job>=FORCE_CONVERGED) iconvergence=FORCE_CONVERGED
     if (.not.conf_para) then
       call m_Files_open_nfcntn                        ! -(m_Files)
     else
       call m_Files_reopen_nfcntn
     endif
     if(mype==0) rewind nfcntn
     call m_Iter_wd_iteration_numbers(nfcntn)
     call m_IS_wd_pos_and_v(nfcntn)
     call m_IS_wd_rigid_body(nfcntn)
     call m_IS_wd_cps(nfout)
     call m_TE_wd_total_energy(nfcntn)
#ifndef _EMPIRICAL_
     call m_CtrlP_wd_isolver(nfcntn)
#endif
     call m_IS_wd_pos_brav(nfout)

     if(imdalg == T_CONTROL .or. imdalg == BLUEMOON .or. &
          & imdalg == QUENCHED_CONSTRAINT) then
        if(imdalg /= QUENCHED_CONSTRAINT) &
             & call m_IS_wd_nrsv(nfcntn)
        call m_IS_wd_forcp_etc(imdalg,nfcntn)
     end if
     if(imdalg == BFGS .or. imdalg == GDIIS)then
        call m_IS_wd_diis_history(nfcntn)
     endif
     call m_CtrlP_wd_iconvergence(nfcntn)
     call m_CtrlP_wd_edelta_ontheway(nfcntn)
! === For restart lm+MSD! by tkato 2012/02/15 ==================================
     call m_CtrlP_wd_dtim_previous(nfcntn)
! ==============================================================================
     if(isSCDFT()) call wd_eps0(nfcntn)
#ifndef _EMPIRICAL_
     call m_ESIO_wd_EigenValues(nfout,2,nooccupation=NO)
     call m_Files_open_nfeng(icond)
     if((icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) &
          & .and. fixed_charge_k_parallel == ALL_AT_ONCE) then
!!$        if ( noncol ) then
        if(mype == 0) write(nfeng,'(" num_kpoints = ",i6)') kv3 /ndim_spinor
        if(neg_is_enlarged) then
           if(mype == 0) write(nfeng,'(" num_bands   = ",i6)') neg-num_extra_bands
        else
           if(mype == 0) write(nfeng,'(" num_bands   = ",i6)') neg
        end if
        if(mype == 0) write(nfeng,'(" nspin       = ",i6)') nspin / ndim_spinor
        if(m_ESoc_check_num_bands()) then
           call m_ESoc_check_if_metalic(nfout)
!!$           if(metalic_system .or. .not.check_if_metalic_flag) then
           if(metalic_system) then
              if(mype == 0) write(nfeng,'(" Fermi energy level = ",f10.6)') efermi
           else
!!$           write(nfeng,'(" The Highest occupied band energy = ",f10.6/)') vbm
              if(mype == 0) write(nfeng,'(" Valence band max   = ",f10.6)') vbm
           end if
        else
           if(mype == 0) write(nfeng,'(" Fermi energy level = unknown")')
        end if
        call m_ESIO_wd_EigenValues(nfeng,2,nooccupation=NO)
     else
        call m_ESIO_wd_EigenValues(nfeng,2,nooccupation=NO)
     end if
     call m_Files_close_nfeng()

     if(sw_positron == ON) then
        call m_CtrlP_set_corecharge_cntnbin(ON)
     else
        call m_CtrlP_set_corecharge_cntnbin(OFF)
     end if
     call m_CtrlP_wd_corecharge_cntnbin(nfcntn)
     call m_CtrlP_wd_neg(nfcntn)
     call m_CS_wd_fix_spin_status(nfcntn)
     call m_IS_wd_curr_atom_reservoir(nfcntn)
     if(sw_stress_correction==ON)then
       call m_Stress_wd_cdata_4_correction(nfcntn)
     endif
     call m_IS_wd_moved_distances_of_planes(nfcntn)

     call m_CS_SG_wd_cntn(nfcntn)

     call m_Files_reopen_nfcntn_bin()
!!$     call m_Files_open_nfcntn_bin
!!$     if(mype==0) rewind nfcntn_bin
     if(continuation_using_ppdata == YES) then
        call m_PP_wd_PP_parameters_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kgp)
        call m_NLP_wd_snl_3D(nfout,nfcntn_bin,F_CNTN_BIN_partitioned,kv3)
     end if
     call m_ESIO_wd_EigenValues_etc(nfcntn_bin,F_CNTN_BIN_partitioned,totch_flag=ON)

     if(sw_optimize_lattice==ON) call m_UC_wd_cntn_data(nfcntn)
     if(imdalg == P_CONTROL .or. imdalg == PT_CONTROL) call m_UC_wd_md_cntn_data(nfcntn)

     if((F_ZAJ_partitioned .neqv. F_ZAJ_in_partitioned) .or. icond == INITIAL .or. conf_para) &
          & call m_Files_reopen_nfzaj()
     if((F_CHGT_partitioned .neqv. F_CHGT_in_partitioned) .or. icond == INITIAL .or. conf_para) &
          & call m_Files_reopen_nfchgt()

! ===================== added by K. Tagami ======================= 12.0Exp
     if ( (icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) .and. &
          &        fixed_charge_k_parallel /= ONE_BY_ONE ) then
        call m_Files_reopen_nfzaj()
     endif
! ================================================================ 12.0Exp

     if ( flg_wd_zaj ) then
        if(sw_ekzaj == OFF) then
           call m_ESIO_wd_WFs_3D(nfout,nfzaj,F_ZAJ_partitioned)
        else
        end if
     endif

     if(sw_write_pwbs_info==ON) call m_pwBS_wd_curr_pws()

! ==== DEBUG TEST ==== 2020/01/28
     if ( icond /= FIXED_CHARGE .and. icond /= FIXED_CHARGE_CONTINUATION ) then
        if(flg_paw) then
           call m_Files_open_nfcntn_bin_paw
           if(mype==0) rewind nfcntn_bin_paw
           call m_PP_wd_PAW_parameters(nfout,nfcntn_bin_paw)
           call m_CD_wd_hsr(nfcntn_bin_paw)
!!$call m_PAWCD_wd_cd(nfout)
           call m_Files_close_nfcntn_bin_paw
        end if

        call m_CD_wd_chgq(nfchgt,F_CHGT_partitioned) ! write chgq_l


! ================================= modified by K. Tagami ============= 11.0
!     if(sw_hubbard == ON) call m_OP_wd_occ_mat(nfout) ! write occupation matrix
!
        if(sw_hubbard == ON) then
           if ( noncol ) then
              call m_OP_wd_occ_mat_noncl(nfout)
           else
              call m_OP_wd_occ_mat(nfout) ! write occupation matrix
           endif
        endif
! ===================================================================== 11.0

        if ( sw_opencore == ON ) call m_PS_wd_mag_opencore_pol
        if ( use_metagga ) then
           call m_Files_open_nf_ekindens()
           call m_KE_wd_ekindens()
           call m_Files_close_nf_ekindens()
        endif
     endif
! ==== DEBUG TEST ==== 2020/01/28

     call m_Files_open_nfefermi()
     call m_ESIO_wd_Efermi(nfout,nfefermi)
     call m_Files_close_nfefermi()

     if(sw_output_hybrid_info == ON .and.  (icond==INITIAL.or.icond==CONTINUATION)) &
     & call m_ES_EXX_wd_hybridinfo()

  else
     if((F_ZAJ_partitioned .neqv. F_ZAJ_in_partitioned) .or. icond == INITIAL) then
        if ( flg_wd_zaj ) then
           if(sw_ekzaj == OFF) then
              call m_Files_reopen_nfzaj()
              call m_ESIO_wd_WFs_3D(nfout,nfzaj,F_ZAJ_partitioned)
           end if
        end if
     endif


     if((F_CHGT_partitioned .neqv. F_CHGT_in_partitioned) .or. icond == INITIAL) then
        call m_Files_reopen_nfchgt()

        call m_CD_wd_chgq(nfchgt,F_CHGT_partitioned) ! write chgq_l

     end if


! ==== DEBUG TEST ==== 2020/01/28
     if(flg_paw) then
        if ( icond /= FIXED_CHARGE .and. icond /= FIXED_CHARGE_CONTINUATION ) then
           call m_Files_open_nfcntn_bin_paw
           if(mype==0) rewind nfcntn_bin_paw
           call m_PP_wd_PAW_parameters(nfout,nfcntn_bin_paw)
           call m_CD_wd_hsr(nfcntn_bin_paw)
           call m_Files_close_nfcntn_bin_paw
!!$call m_PAWCD_wd_cd(nfout)
        end if
     endif
! ==== DEBUG TEST ==== 2020/01/28
#endif
  end if
#ifndef _EMPIRICAL_
  end if
#endif

! ============================= KT_Test ================ 12.5Exp
#ifdef USE_ZAJ_HISTORY
  call m_ES_wf_wd_zaj_history
#endif
! ====================================================== 12.5Exp

  if(iprijobstatus >= 1) then
!!$     write(nfout,'(" WDD onto Files")')
     call m_Files_open_nfstatus()
     if(jobstatus_series == ON) then
!!$        switch_header = OFF
        call m_Files_skiptoend(nfstatus)
     else
!!$        switch_header = ON
     end if
     status_wdmode = FINISH
     call m_Timing_wd_status(nfstatus,jobstatus_format,jobstatus_series,status_wdmode &
          &                , iteration,iteration_ionic,iteration_electronic)
  end if
  if(ipritiming>=2) call tstatc_wd0

  if (conf_para.and.driver==DRIVER_CONSTRAINT) call mpi_barrier(mpi_comm_world,mpi_err)

  if(final_c) then
     call m_Files_close_all()                    ! -(m_Files)
     call PrintStatus()
     call m_Files_close_logfile()
  else
     call m_Files_close_nfinp()
     call m_Files_close_nfcntn()
!!$     if(mype==0)then
!!$        close(nfcntn)
!!$        close(nfinp)
!!$     endif
  endif
contains
  integer function bcast_ipri(ipri_in)
    integer, intent(in)  :: ipri_in
    integer :: ipri_out, ierr
    if(npes > 1) then
       if(mype == 0) ipri_out = ipri_in
       call mpi_bcast(ipri_out,1,mpi_integer,0,MPI_CommGroup,ierr)
    else
       ipri_out = ipri_in
    end if
!!$    if(ipri_out /= OFF) ipri_out = ON
    if(ipri_out /= ON) ipri_out = OFF
    bcast_ipri = ipri_out
  end function bcast_ipri
end subroutine WriteDownData_onto_Files

!!$subroutine WriteCheckPointData
!!$  use m_IterationNumbers,   only : m_Iter_wd_iteration_numbers_b, iteration_electronic
!!$  use m_Files,              only : m_Files_open_nfcheckpoint, m_Files_close_nfcheckpoint &
!!$       &                         , nfout, nfcheckpoint
!!$  use m_ES_IO,              only : m_ESIO_wd_EigenValues_etc, m_ESIO_wd_WFs
!!$  use m_Charge_Density,     only : m_CD_wd_chgq
!!$  implicit none
!!$
!!$  if((icond == CONTINUATION or. icond == INITIAL).and. iteration_electronic == 10) then
!!$     call m_Files_open_nfcheckpoint
!!$     call m_Iter_wd_iteration_numbers_b(nfcheckpoint,F_CHKPNT_partitioned)
!!$     call m_ESIO_wd_EigenValues_etc(nfcheckpoint,F_CHKPNT_partitioned,totch_flag=ON)
!!$     call m_ESIO_wd_WFs(nfout,nfcheckpoint,F_CHKPNT_partitioned)
!!$     call m_CD_wd_chgq(nfcheckpoint, F_CHKPNT_partitioned)
!!$     call m_Files_close_nfcheckpoint
!!$  end if
!!$end subroutine WriteCheckPointData
!!$
!!$subroutine ReadCheckPointData_if_needed
!!$  use m_IterationNumbers, only : iteration_ionic, iteration_electronic
!!$  implicit none
!!$  if((icond == CONTINUATION or. icond == INITIAL).and. iteration_electronic == 15) then
!!$     call m_Files_open_nfcheckpoint
!!$     call m_Iter_rd_iteration_numbers_b(nfcheckpoint,F_CHKPNT_partitioned)
!!$     call m_ESIO_rd_EigenValues_etc(
!!$
!!$
!!$end subroutine ReadCheckPointData_if_needed

subroutine PrintStatus()
  use m_Files, only : nfout, m_Files_open_nfstatus, nfstatus, m_Files_skiptoend
  use m_Const_Parameters, only : MAX_SCF_ITERATION_REACHED, MAX_TIME_REACHED, FSTOP &
     &                         , FORCE_CONVERGENCE_REACHED, STRESS_CONVERGENCE_REACHED, MAX_MDSTEPS_REACHED &
     &                         , CHG_CONVERGENCE_REACHED &
     &                         , FORCE_CONVERGED, WF_CONVERGENCE_REACHED, MAX_PHSTEPS_REACHED &
     &                         , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, FINISH, ON
  use m_Control_Parameters,  only : terminated_because, iconvergence, iconvergence_previous_job, icond
  use m_IterationNumbers, only : iteration_ionic, iteration_unit_cell, iteration_stress_correction

  implicit none
  integer, dimension(8) :: ipresent_time
  character(len('reason: ')), parameter :: tag_reason = 'reason: '
  integer :: status_wdmode

  call date_and_time(values=ipresent_time)
  write(nfout,*)
write(nfout,'(a)') '------------------------------------------------------------------------'
  write(nfout,'("PHASE/0 terminated at ",i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4," reason: ")') &
        & ipresent_time(5:7),ipresent_time(3),ipresent_time(2),ipresent_time(1)
!write(nfout,'(a)') "RRRRRR   EEEEEEE     A      SSSSS   OOOOOOO  N     N  ::: "
!write(nfout,'(a)') "R     R  E          A A    S     S  O     O  NN    N  ::: "
!write(nfout,'(a)') "R     R  E         A   A   S        O     O  N N   N  ::: "
!write(nfout,'(a)') "RRRRRR   EEEEE    A     A   SSSSS   O     O  N  N  N       "
!write(nfout,'(a)') "R   R    E        A######        S  O     O  N   N N  ::: "
!write(nfout,'(a)') "R    R   E        A     A  S     S  O     O  N    NN  ::: "
!write(nfout,'(a)') "R     R  EEEEEEE  A     A   SSSSS   OOOOOOO  N     N  ::: "
  select case (terminated_because)
  case (MAX_SCF_ITERATION_REACHED)
!    write(nfout,'(a)') tag_reason//'scf iteration reached max_iteration'
write(nfout,'(a)') "M     M     A     X     X  III  TTTTTTT  EEEEEEE  RRRRRR  "
write(nfout,'(a)') "MM   MM    A A     X   X    I      T     E        R     R "
write(nfout,'(a)') "M M M M   A   A     X X     I      T     E        R     R "
write(nfout,'(a)') "M  M  M  A     A     X      I      T     EEEEE    RRRRRR  "
write(nfout,'(a)') "M     M  AAAAAAA    X X     I      T     E        R   R   "
write(nfout,'(a)') "M     M  A     A   X   X    I      T     E        R    R  "
write(nfout,'(a)') "M     M  A     A  X     X  III     T     EEEEEEE  R     R "
  case (MAX_TIME_REACHED)
!    write(nfout,'(a)') tag_reason//'max. time reached cpumax'
write(nfout,'(a)') " CCCCC   PPPPPP   U     U  M     M     A     X     X  "
write(nfout,'(a)') "C     C  P     P  U     U  MM   MM    A A     X   X   "
write(nfout,'(a)') "C        P     P  U     U  M M M M   A   A     X X    "
write(nfout,'(a)') "C        PPPPPP   U     U  M  M  M  A     A     X     "
write(nfout,'(a)') "C        P        U     U  M     M  AAAAAAA    X X    "
write(nfout,'(a)') "C     C  P        U     U  M     M  A     A   X   X   "
write(nfout,'(a)') " CCCCC   P         UUUUU   M     M  A     A  X     X  "
  case (MAX_MDSTEPS_REACHED)
!    write(nfout,'(a)') tag_reason//'mdsteps reached max_mdstep'
write(nfout,'(a)') "M     M     A     X     X  M     M  DDDDDD    SSSSS   TTTTTTT  PPPPPP   "
write(nfout,'(a)') "MM   MM    A A     X   X   MM   MM  D     D  S     S     T     P     P  "
write(nfout,'(a)') "M M M M   A   A     X X    M M M M  D     D  S           T     P     P  "
write(nfout,'(a)') "M  M  M  A     A     X     M  M  M  D     D   SSSSS      T     PPPPPP   "
write(nfout,'(a)') "M     M  AAAAAAA    X X    M     M  D     D        S     T     P        "
write(nfout,'(a)') "M     M  A     A   X   X   M     M  D     D  S     S     T     P        "
write(nfout,'(a)') "M     M  A     A  X     X  M     M  DDDDDD    SSSSS      T     P        "
  case (FSTOP)
!    write(nfout,'(a)') tag_reason//'integer smaller than the current iteration was detected from the F_STOP file'
write(nfout,'(a)') "FFFFFFF       SSSSS   TTTTTTT  OOOOOOO  PPPPPP   "
write(nfout,'(a)') "F            S     S     T     O     O  P     P  "
write(nfout,'(a)') "F            S           T     O     O  P     P  "
write(nfout,'(a)') "FFFFF         SSSSS      T     O     O  PPPPPP   "
write(nfout,'(a)') "F                  S     T     O     O  P        "
write(nfout,'(a)') "F            S     S     T     O     O  P        "
write(nfout,'(a)') "F     _____   SSSSS      T     OOOOOOO  P        "
  case (FORCE_CONVERGENCE_REACHED)
!    write(nfout,'(a)') tag_reason//'force converged'
write(nfout,'(a)') "FFFFFFF  OOOOOOO  RRRRRR    CCCCC     CCCCC   OOOOOOO  N     N  V     V "
write(nfout,'(a)') "F        O     O  R     R  C     C   C     C  O     O  NN    N  V     V "
write(nfout,'(a)') "F        O     O  R     R  C         C        O     O  N N   N  V     V "
write(nfout,'(a)') "FFFFF    O     O  RRRRRR   C         C        O     O  N  N  N  V     V "
write(nfout,'(a)') "F        O     O  R   R    C         C        O     O  N   N N   V   V  "
write(nfout,'(a)') "F        O     O  R    R   C     C   C     C  O     O  N    NN    V V   "
write(nfout,'(a)') "F        OOOOOOO  R     R   CCCCC     CCCCC   OOOOOOO  N     N     V    "
  case (STRESS_CONVERGENCE_REACHED)
!   write(nfout,'(a)') tag_reason//'stress converged'
write(nfout,'(a)') " SSSSS   TTTTTTT  RRRRRR    SSSSS     CCCCC   OOOOOOO  N     N  V     V "
write(nfout,'(a)') "S     S     T     R     R  S     S   C     C  O     O  NN    N  V     V "
write(nfout,'(a)') "S           T     R     R  S         C        O     O  N N   N  V     V "
write(nfout,'(a)') " SSSSS      T     RRRRRR    SSSSS    C        O     O  N  N  N  V     V "
write(nfout,'(a)') "      S     T     R   R          S   C        O     O  N   N N   V   V  "
write(nfout,'(a)') "S     S     T     R    R   S     S   C     C  O     O  N    NN    V V   "
write(nfout,'(a)') " SSSSS      T     R     R   SSSSS     CCCCC   OOOOOOO  N     N     V    "
  case(CHG_CONVERGENCE_REACHED)
write(nfout,'(a)') " CCCCC   H     H    GGGGG      CCCCC   OOOOOOO  N     N  V     V"
write(nfout,'(a)') "C     C  H     H   G     G    C     C  O     O  NN    N  V     V"
write(nfout,'(a)') "C        H     H   G          C        O     O  N N   N  V     V"
write(nfout,'(a)') "C        HHHHHHH   G  GGGG    C        O     O  N  N  N  V     V"
write(nfout,'(a)') "C        H     H   G     G    C        O     O  N   N N   V   V "
write(nfout,'(a)') "C     C  H     H   G     G    C     C  O     O  N    NN    V V  "
write(nfout,'(a)') " CCCCC   H     H    GGGGG      CCCCC   OOOOOOO  N     N     V   "
  case(WF_CONVERGENCE_REACHED)
write(nfout,'(a)') "W     W  FFFFFFF    CCCCC   OOOOOOO  N     N  V     V "
write(nfout,'(a)') "W  W  W  F         C     C  O     O  NN    N  V     V "
write(nfout,'(a)') "W  W  W  F         C        O     O  N N   N  V     V "
write(nfout,'(a)') "W  W  W  FFFFF     C        O     O  N  N  N  V     V "
write(nfout,'(a)') "W  W  W  F         C        O     O  N   N N   V   V  "
write(nfout,'(a)') "W  W  W  F         C     C  O     O  N    NN    V V   "
write(nfout,'(a)') " WW WW   F          CCCCC   OOOOOOO  N     N     V    "
  case(MAX_PHSTEPS_REACHED)
write(nfout,'(a)') "M     M     A     X     X  PPPPPP   H     H   SSSSS   TTTTTTT  PPPPPP   "
write(nfout,'(a)') "MM   MM    A A     X   X   P     P  H     H  S     S     T     P     P  "
write(nfout,'(a)') "M M M M   A   A     X X    P     P  H     H  S           T     P     P  "
write(nfout,'(a)') "M  M  M  A     A     X     PPPPPP   HHHHHHH   SSSSS      T     PPPPPP   "
write(nfout,'(a)') "M     M  AAAAAAA    X X    P        H     H        S     T     P        "
write(nfout,'(a)') "M     M  A     A   X   X   P        H     H  S     S     T     P        "
write(nfout,'(a)') "M     M  A     A  X     X  P        H     H   SSSSS      T     P        "
  case default
    write(nfout,'(a)') 'unknown'
  end select
  if(terminated_because==FORCE_CONVERGENCE_REACHED  .or. &
     terminated_because==STRESS_CONVERGENCE_REACHED .or. &
     terminated_because==MAX_MDSTEPS_REACHED        .or.&
     terminated_because==CHG_CONVERGENCE_REACHED    .or. &
     terminated_because==WF_CONVERGENCE_REACHED     .or. &
     terminated_because==MAX_PHSTEPS_REACHED) then
write(nfout,'(a)') '------------------------------------------------------------------------'
     return
  endif
  if(.not. (iconvergence_previous_job>=FORCE_CONVERGED) .and. &
   & .not. (icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) &
   & .and.  iconvergence == 0 .and. iteration_ionic==1 .and.  iteration_unit_cell==1 &
   & .and. iteration_stress_correction==1 ) then
!     write(nfout,'(a)') 'WARNING: charge unconverged!!'
     write(nfout,'(a)')
write(nfout,'(a)') "W     W     A     RRRRRR   N     N  III  N     N   GGGGG  "
write(nfout,'(a)') "W  W  W    A A    R     R  NN    N   I   NN    N  G     G "
write(nfout,'(a)') "W  W  W   A   A   R     R  N N   N   I   N N   N  G       "
write(nfout,'(a)') "W  W  W  A     A  RRRRRR   N  N  N   I   N  N  N  G  GGGG "
write(nfout,'(a)') "W  W  W  AAAAAAA  R   R    N   N N   I   N   N N  G     G "
write(nfout,'(a)') "W  W  W  A     A  R    R   N    NN   I   N    NN  G     G "
write(nfout,'(a)') " WW WW   A     A  R     R  N     N  III  N     N   GGGGG  "
write(nfout,'(a)') " CCCCC   H     H    GGGGG    U     U   CCCCC   OOOOOOO  N     N  V     V"
write(nfout,'(a)') "C     C  H     H   G     G   U     U  C     C  O     O  NN    N  V     V"
write(nfout,'(a)') "C        H     H   G         U     U  C        O     O  N N   N  V     V"
write(nfout,'(a)') "C        HHHHHHH   G  GGGG   U     U  C        O     O  N  N  N  V     V"
write(nfout,'(a)') "C        H     H   G     G   U     U  C        O     O  N   N N   V   V "
write(nfout,'(a)') "C     C  H     H   G     G   U     U  C     C  O     O  N    NN    V V  "
write(nfout,'(a)') " CCCCC   H     H    GGGGG     UUUUU    CCCCC   OOOOOOO  N     N     V   "
  endif
write(nfout,'(a)') '------------------------------------------------------------------------'
end subroutine PrintStatus

