!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  MAIN PROGRAM: PHASE
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
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
!  $Id: mdmain0.F90 633 2020-12-01 05:11:03Z jkoga $
!
program PHASE
#ifdef NEC_TIMER
  use nec_timer
#endif
  use m_Const_Parameters, only : DRIVER_CONSTRAINT, DRIVER_NEB, DRIVER_MTD, DRIVER_SC_DFT, DRIVER_DIMER
  use m_Parallelization,      only : mype
  implicit none
  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Already_Converged, Already_Converged2
  logical  :: Positron_scf, Positron_nonscf
  logical  :: Hubbard_model
  logical  :: Forces_are_Converged, Ending_Time, Force_errors_are_tolerable,UnitCell_Converged &
  &         , Forces_are_Converged2
!!$  logical  :: ChargeDensity_is_Fixed
#ifdef NEC_ITER_REG
  integer  :: count_for_ftrace
#endif

  integer :: initmpi
  logical :: confpara
  integer :: driver
  logical :: ending_t,force_conv,mpi_initialized,uconv
  logical :: initialization_required
  
#ifdef NEC_ITER_REG
     count_for_ftrace = 0
     call FTRACE_REGION_BEGIN("INITIAL")
#endif

  confpara = Resolve_Config_Parallel()
  initmpi=1
  if(confpara) initmpi=0
  mpi_initialized = .false.

  do

  ending_t = .false.
  if (mpi_initialized) initmpi=0
  call Initialization(initmpi)
  mpi_initialized = .true.
  call InputData_Analysis
  driver   = Resolve_Driver()

  if(ekmode()) call Initialization_set_ekmode_ON()

  if(driver == DRIVER_NEB) then
    call do_neb()
  else if (driver==DRIVER_DIMER) then
    call do_dimer_method()
  else if(driver==DRIVER_CONSTRAINT) then
#ifndef DISABLE_CONSTRAINTS
    call constrained_dynamics()
#endif
  else if(driver==DRIVER_MTD) then
#ifndef DISABLE_CONSTRAINTS
    call meta_dynamics()
#endif
  else

  outer_loop:do
  if(driver==DRIVER_SC_DFT .and. .not.Initial_SCDFTLoop()) then
     call InputData_Analysis()
     if(ekmode()) call Initialization_set_ekmode_ON()
  end if

  call Preparation(0)                ! Basis set, symmetry check etc.
  if(initialization_required())then
     call Preparation_for_mpi(1)     ! mpi
  endif
  call PseudoPotential_Construction
#ifdef ENABLE_ESM_PACK
  if(initialization_required())then
     call Preparation_for_ESM
  endif
#endif

  call Ewald_and_Structure_Factor
  if ( Positron_scf() ) call Initial_pWaveFunctions()

  call Initial_Electronic_Structure

  if(ChargeDensity_is_Fixed() .and. One_by_one_in_each_rank_k()) then ! icond=2, 3
     call ekcal()  ! contained here
     exit outer_loop
  else
     call Initial_MD_Condition
#ifdef NEC_ITER_REG
     call FTRACE_REGION_END("INITIAL")
#endif
     force_conv = Already_Converged()
     if(.not.force_conv) then
#ifdef NEC_ITER_REG
        call FTRACE_REGION_BEGIN("SOLVE-FIRST")
#endif
        StressLoop: do
           AtomicConfiguration: do
              force_conv = Forces_are_Converged2()
              if(force_conv) then
                exit AtomicConfiguration
              endif
              ChargeDensity:    do
                 force_conv=.false.
#ifdef NEC_ITER_REG
                 count_for_ftrace = count_for_ftrace + 1
                 if(count_for_ftrace .eq. 2) then
                    call FTRACE_REGION_END("SOLVE-FIRST")
                    call FTRACE_REGION_BEGIN("SOLVE-CORE")
                 end if
#endif
                 call IterationNumber_Setting

! ============================ added by K. Tagami ============- 5.0
                 call Renewal_of_Chg_Ctrl_Param
! ============================================================ 5.0

                 call Renewal_of_WaveFunctions
                 if ( Positron_scf() ) call Renewal_of_pWaveFunctions                 

                 call ChargeDensity_Construction(1)
                 call Potential_Construction

                 if ( PotentialMix() ) then
                    call Potential_Mixing
                 else
                    call ChargeDensity_Mixing
                 endif

                 ending_t = Ending_Time()
                 if(ending_t)                      exit StressLoop
                 if(TotalEnergy_is_Divergent())    exit StressLoop

                 if ( PotentialMix() ) then
                 else
                    call Renewal_of_Potential
                    if ( Positron_scf() ) call Renewal_of_pPotential
                    if (Hubbard_model() ) then
                       call Renewal_of_Hubbard_Parameters
                       call Renewal_of_Hubbard_Potential
                    end if
                 endif

                 if(ChargeDensity_is_Converged())  exit ChargeDensity
              enddo ChargeDensity
#ifdef LIBRARY_BUILD
              call Forces_phase0
#else
              call Forces
#endif
              if ( epsilon_post_scf() ) call Epsilon_PostScf

              force_conv = Forces_are_Converged()
              if ( Enforce_Update_Cell() )  exit AtomicConfiguration

              if(force_conv) exit AtomicConfiguration
              if(Force_errors_are_tolerable()) then
                 call Postprocessing_during_MD()
                 call Move_Ions
                 call MDIterationNumber_Setting
                 call Ewald_and_Structure_Factor
                 if ( Hubbard_model() ) then
                    call Renewal_of_Hubbard_Parameters
                    call Renewal_of_Hubbard_Potential
                 end if
!!$                 call MDIterationNumber_Setting
              end if
!              if ( Enforce_Update_Cell() )  exit AtomicConfiguration

              if(BreakMD(force_conv))then
                 exit AtomicConfiguration
              endif
              call Postproc_after_SCF_convergence()
           enddo AtomicConfiguration
           if(driver==DRIVER_SC_DFT) exit StressLoop
           call Stress
           exit StressLoop
        end do StressLoop
#ifdef NEC_ITER_REG
        if(count_for_ftrace .eq. 1) then
           call FTRACE_REGION_END("SOLVE-FIRST")
        else
           call FTRACE_REGION_END("SOLVE-CORE")
        end if
#endif
     end if

     if (driver == DRIVER_SC_DFT .and. .not.ending_t) then
! * T.Hamada 2016.1.8 + T. Yamasaki 2016.02.21
        write(6,'(" driver = DRIVER_SC_DFT")')
        call flush(6)
        call Epsilon_Paramset
        call Epsilon_Postscf
        call setalpha_exx()
! * T.Hamada 2015.1.8
        if(Break_SC_DFT()) exit outer_loop
!!$! * T.Hamada 2015.11.4
        call SCDFTIterationNumber_Setting
        call Set_Icond
        call WriteDownData_onto_Files(.false.)
        if(.not.ending_t) call Array_Deallocate()
     else if ( Already_Converged2() ) then
        if ( Positron_nonscf() ) then
#ifdef NEC_ITER_REG
           call FTRACE_REGION_BEGIN("POSITRON")
#endif
           call Initial_pWaveFunctions()
           call Renewal_of_pPotential()
           call Solve_pWaveFunctions()
#ifdef NEC_ITER_REG
           call FTRACE_REGION_END("POSITRON")
#endif
        else if ( Positron_scf() ) then
           call Write_Positron_LifeTime
        end if
        exit outer_loop
     else
        exit outer_loop
     endif

#ifdef NEC_ITER_REG
     call FTRACE_REGION_BEGIN("FINAL")
#endif
  end if

  enddo outer_loop
  uconv = force_conv
  if(.not. ending_t) &
  uconv = UnitCell_Converged(force_conv)
  if(.not.uconv.and..not.ending_t)then
     call MDIterationNumber_Setting2()
     call WriteDownData_onto_Files(.false.)
  else
     if(ChargeDensity_is_Fixed() .and. One_by_one_in_each_rank_k()) then ! icond=2, 3
        !call Postprocessing(.false.)
        !call WriteDownData_onto_Files_ek()
     else
        call WriteDownData_onto_Files(.false.)
        if(uconv) call Postprocessing(.false.)
        if ( uconv .and. epsilon_post_proc() ) call Epsilon_PostScf
        call rttddft_main
!        call WriteDownData_onto_Files(ending_t.or.uconv)
        if(ending_t .or. uconv) call final_output()
        if(uconv) exit
     end if
  endif
  if(.not.ending_t) call MDiterationNumber_Setting_pre()

  if(.not.ending_t) then
     call Array_Deallocate()
     call Continuation_Mode()
  endif

  endif

  if (ending_t.or.OneShot()) exit

  enddo
#ifdef NEC_TIMER
  call print_timer()
#endif
  if(driver/=DRIVER_NEB) call Finalization_of_mpi           ! mpi
#ifdef NEC_ITER_REG
  call FTRACE_REGION_END("FINAL")
#endif
contains
  subroutine final_output()
    use m_Files, only : m_Files_close_all, m_Files_close_logfile
    use m_Control_Parameters, only : m_CtrlP_wd_cpu_total
    implicit none
    call m_CtrlP_wd_cpu_total()
    call m_Files_close_all()                    ! -(m_Files)
    call PrintStatus()
    call m_Files_close_logfile()
  end subroutine final_output

  subroutine ekcal()
    use m_Epsilon_ek, only : sw_epsilon
    use m_Const_Parameters, only : ON
    logical :: AllKpoints_are_Calculated2,AllKpoints_are_Calculated, Already_Converged_for_Kgroup
    logical :: EigenValues_are_Converged, AllKpoints_are_Converged, Already_Converged
    integer :: nk
    logical :: all_conv = .true.
! --------------------T. Hamada 2021.9.28 --------------------
    logical  :: EPS_Prep_Check
!-------------------------------------------------------------
    nk = 0
    if(sw_epsilon == ON) then
       call Initialization_Epsilon        ! Epsilon
       call Shift_Kpoint
    endif
    KPOINT_GROUP: do
       call KpointNumber_Setting2()
       if(AllKpoints_are_Converged()) exit KPOINT_GROUP
       call Preparation_ek()
       call Preparation_for_mpi_ek
       if(sw_epsilon==ON)then
          call PseudoPotential_ek_Epsilon
       else
          call PseudoPotential_ek
       endif
       call Initial_WaveFunctions_ek
       if(.not.Already_Converged_for_kgroup()) then
          all_conv = .false.
          SolveWaveFunctions: do
!-------------------- T. Hamada 2021.9.28 --------------------
             if(EPS_Prep_Check())              exit KPOINT_GROUP
! ------------------------------------------------------------
             if(Ending_Time())                 exit KPOINT_GROUP
             call IterationNumber_Setting()
             call Renewal_of_WaveFunctions()
             if(EigenValues_are_Converged()) then
                if(sw_epsilon==ON)then
                   call Transition_moment_Epsilon                           ! Epsilon
                   call Dealloc_Radr_and_Wos_Epsilon                        ! Epsilon
                endif
                exit SolveWaveFunctions
             endif
          enddo SolveWaveFunctions
          call Postprocessing_k()
          if(AllKpoints_are_Calculated2(nk)) then
             all_conv = .true.
             exit KPOINT_GROUP
          endif
       else
          exit KPOINT_GROUP
       end if
    enddo KPOINT_GROUP
    call Postprocessing(.false.)
!-------------------- T. Hamada 2021.9.28 --------------------
    if((all_conv.and.sw_epsilon==ON) .or.EPS_Prep_Check()) then
!-------------------------------------------------------------
      call Reset_Kpoint                                                   ! Epsilon
      call Prep_for_Calc_Epsilon                                          ! Epsilon
      call Calc_Epsilon                                                   ! Epsilon
      call Calc_Nonlinear_optics                                          ! Epsilon
      call WriteDownData_onto_Files_Eps                                   ! Epsilon
    end if
    call WriteDownData_onto_Files_ek()
  end subroutine ekcal

  subroutine Continuation_Mode()
    use m_Const_Parameters,   only : COORDINATE_CONTINUATION,ON,OFF
    use m_Control_Parameters, only : icond, sw_optimize_lattice,sw_rebuild_pws
    implicit none
    logical :: unitcell_can_change
    if(unitcell_can_change() .and. sw_rebuild_pws==OFF) return
    icond = COORDINATE_CONTINUATION
  end subroutine continuation_mode

  subroutine Postprocessing_scf()
    use m_Const_Parameters,     only : ON
    use m_Control_Parameters,   only : sw_rsb,neg,damp,noncol
    use m_Files,                only : nfout
    use m_ES_WF_by_submat,      only : m_ESsubmat_alloc,m_ESsubmat_renew_WF,m_ESsubmat_renew_WF_noncl
    use m_Electronic_Structure, only : m_ES_energy_eigen_values
    use m_ES_occup,             only : m_ESoc_fermi_parabolic
    use m_ES_IO,                only : m_ESIO_wd_EigenValues
    implicit none
    if(sw_rsb==ON)then
       call m_ESsubmat_alloc()
#ifndef DISABLE_NONCL
! ============================================ modified by K. Tagami ======== 11.0
       if ( noncol ) then
          call m_ESsubmat_renew_WF_noncl(nfout,neg,damp)
       else
#endif
          call m_ESsubmat_renew_WF(nfout,neg,damp)
#ifndef DISABLE_NONCL
       endif
! =========================================================================== 11.0
#endif
       call ChargeDensity_Construction(1)
       call Renewal_of_Potential()
       call m_ES_energy_eigen_values(nfout,.true.)
       call m_ESoc_fermi_parabolic(nfout)
       call m_ESIO_wd_EigenValues(0,2,0)
    endif
  end subroutine Postprocessing_scf

  subroutine Postprocessing_during_MD()
    use m_IterationNumbers, only : iteration_ionic
    use m_Control_Parameters, only : postproc_frequency
    if(postproc_frequency<=0) return
    if(mod(iteration_ionic,postproc_frequency)==0)then
       call Postprocessing(.true.)
    endif
  end subroutine Postprocessing_during_MD

  logical function Enforce_Update_Cell()
    use m_Control_Parameters,   only : sw_optimize_coords_sametime
    use m_Const_Parameters,     only : ON

    Enforce_update_Cell = .false.
    if ( sw_optimize_coords_sametime == ON ) then
       Enforce_update_Cell = .true.
    endif
  end function Enforce_Update_Cell
  
  logical function BreakMD(conv)
    use m_IterationNumbers, only : iteration_ionic
    use m_Ionic_System, only : addition_frequency,m_IS_natm_can_change
    logical, intent(in) :: conv
    if(.not.m_IS_natm_can_change())then
       BreakMD = conv
       return
    endif
    BreakMD = mod(iteration_ionic,addition_frequency)==0
  end function BreakMD

  logical function Break_SC_DFT()
    use m_IterationNumbers, only : iteration, first_iteration_of_this_job, iteration_scdft
    use m_Control_Parameters, only : delta_epsilon, max_scdft_iteration, epsilon0, epsilon0_previous
    if(iteration-first_iteration_of_this_job <= 0 .or. iteration_scdft<=1) then
       Break_SC_DFT = .false.
    else
       if(dabs(epsilon0 - epsilon0_previous) < max(0.d0, delta_epsilon)) then
          Break_SC_DFT = .true.
       else
          if(iteration_scdft > max_scdft_iteration) then
             Break_SC_DFT = .true.
          else
             Break_SC_DFT = .false.
          end if
       end if
    end if
    write(6,'(" Break_SC_DFT = ",L3)') Break_SC_DFT
    call flush(6)
  end function Break_SC_DFT

  logical function Initial_SCDFTLoop()
    use m_IterationNumbers, only : iteration_scdft, iteration_scdft_initial
    Initial_SCDFTLoop = .true.
    if(Iteration_scdft > Iteration_scdft_initial) Initial_SCDFTLoop = .false.
!!$    if(mype==0) then
       write(6,*) ' ! Iteration_scdft, iteration_scdft_initial = ', iteration_scdft,iteration_scdft_initial
       write(6,*) ' I Initial_SCDFTLoop = ',Initial_SCDFTLoop
       call flush(6)
!!$    end if
  end function Initial_SCDFTLoop

  subroutine setalpha_exx()
    use m_Control_Parameters, only : m_CtrlP_set_alpha_exx, epsilon0
    use m_Files,              only : nfout
    use m_Const_Parameters,   only : DP
    real(kind=DP) :: alpha
    alpha = 1.d0/epsilon0
    call m_CtrlP_set_alpha_exx(nfout,alpha)
  end subroutine setalpha_exx

  logical function OneShot()
    use m_Ionic_System, only : m_IS_natm_can_change
    use m_Control_Parameters, only : sw_optimize_lattice
    use m_Const_Parameters, only : OFF
    logical :: unitcell_can_change
    OneShot = .not.m_IS_natm_can_change()
    if(OneShot) then
        OneShot = .not.unitcell_can_change()
    endif
  end function OneShot

! ======= KT_add ============ 13.0U2
  logical function PotentialMix()
    use m_Control_Parameters, only : sw_potential_mixing
    use m_Const_Parameters,   only : ON

    if ( sw_potential_mixing == ON ) then
       PotentialMix = .true.
    else
       PotentialMix = .false.
    endif
  end function PotentialMix
! =========================== 13.0U2

! === KT_add ==== 13.1R
  logical function epsilon_post_scf()
    use m_Const_Parameters,   only : ON, OFF
    use m_Epsilon_ek,  only : sw_epsilon, m_Eps_chkif_sw_epsilon
    use m_Control_Parameters, only : sw_phonon, sw_phonon_with_epsilon, &
         &                           sw_calc_dielectric_tensor, sw_excitation

    logical, save :: First = .true.

    epsilon_post_scf = .false.

!    if ( First ) then
!       call m_Eps_chkif_sw_epsilon;        First = .false.
!    endif
    if ( sw_epsilon == OFF .and. sw_excitation == OFF ) return

    if ( sw_phonon == ON ) then
       if ( sw_phonon_with_epsilon == ON .and. sw_calc_dielectric_tensor == ON ) then
          epsilon_post_scf = .true.
       endif
    else
!       epsilon_post_scf = .true.
       epsilon_post_scf = .false.
    endif

  end function epsilon_post_scf

  logical function epsilon_post_proc()
    use m_Const_Parameters,   only : ON, OFF
    use m_Epsilon_ek,  only : sw_epsilon, m_Eps_chkif_sw_epsilon
    use m_Control_Parameters, only : sw_phonon, sw_phonon_with_epsilon, &
         &                           sw_calc_dielectric_tensor, sw_excitation

    logical, save :: First = .true.

    epsilon_post_proc = .false.

!    if ( First ) then
!       call m_Eps_chkif_sw_epsilon;        First = .false.
!    endif
    if ( sw_epsilon == OFF .and. sw_excitation == OFF ) return

    if ( sw_phonon == ON ) then
       epsilon_post_proc = .false.
    else
       epsilon_post_proc = .true.
    endif

  end function epsilon_post_proc
! =============== 13.1R

  logical function ChargeDensity_is_Fixed()
    use m_Control_Parameters, only : icond
    use m_Const_Parameters,   only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION
    if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
       ChargeDensity_is_Fixed = .true.
    else
       ChargeDensity_is_Fixed = .false.
    end if
  end function ChargeDensity_is_Fixed

  logical function One_by_one_in_each_rank_k()
    use m_Control_Parameters, only : fixed_charge_k_parallel
    use m_Const_Parameters,   only : ONE_BY_ONE
    if(fixed_charge_k_parallel == ONE_BY_ONE) then
       One_by_one_in_each_rank_k = .true.
    else
       One_by_one_in_each_rank_k = .false.
    end if
  end function One_by_one_in_each_rank_k

  integer function Resolve_Driver()
    use m_Control_Parameters, only : driver
    Resolve_Driver = driver
  end function Resolve_Driver

  logical function ekmode()
    use m_Control_Parameters, only : icond,fixed_charge_k_parallel
    use m_Const_Parameters, only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION,ONE_BY_ONE
    implicit none
    ekmode = .false.
    if((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION) .and. &
    &  fixed_charge_k_parallel .eq. ONE_BY_ONE) then
       ekmode = .true.
    endif
  end function ekmode

  logical function Resolve_Config_Parallel()
    use m_Parallelization, only : m_Parallel_resolve_conf_para
    implicit none
    Resolve_Config_Parallel = m_Parallel_resolve_conf_para()
  end function Resolve_Config_Parallel

   subroutine Set_Icond
     use m_Const_Parameters,   only : INITIAL
     use m_Control_Parameters, only : m_CtrlP_set_icond
     call m_CtrlP_set_icond(INITIAL)
   end subroutine Set_Icond

  subroutine Array_Deallocate()
     use m_Const_Parameters, only : OFF, ON, DRIVER_URAMP, DRIVER_SC_DFT, DRIVER_NEB
     use m_Control_Parameters, only : m_CtrlP_dealloc,sw_rebuild_pws,m_CtrlP_set_init_status, driver
     use m_Crystal_Structure, only : m_CS_dealloc
     use m_Ionic_System, only : m_IS_dealloc
     use m_PlaneWaveBasisSet, only : m_pwBS_dealloc
     use m_Parallelization, only : m_Parallel_dealloc,m_Parallel_dealloc_mpi_nlmta,m_Parallel_dealloc_mpi_exx &
                        &        , m_Parallel_cp_g1k
     use m_Kpoints, only : m_Kp_dealloc, kv3
     use m_Force, only : m_Force_dealloc
     use m_Charge_Density, only : m_CD_dealloc
     use m_XC_Potential, only : m_XC_dealloc_vxc
     use m_PseudoPotential, only : m_PP_dealloc, flg_paw
     use m_NonLocal_Potential, only : m_NLP_dealloc
     use m_Electronic_Structure, only : m_ES_dealloc
     use m_ES_WF_by_SDorCG, only : m_ESsd_dealloc
     use m_PAW_ChargeDensity, only : m_PAW_dealloc
     use m_PAW_XC_Potential, only : m_PAW_XC_dealloc_vxc
     use m_ES_wf_extrpl, only : m_ES_wf_extrpl_dealloc

! ===== KT_add ============== 13.0AS
     use m_Control_Parameters,  only : num_projectors, sw_hubbard
     use m_Orbital_Population,only: m_OP_dealloc
     use m_Electronic_Structure,  only : m_ES_dealloc_Dhub
! ============================13.0AS

! ======== KT_add ============= 2013/10/31
     use m_Control_Parameters, only: noncol, SpinOrbit_Mode
     use m_Const_Parameters,   only: Neglected, BuiltIn, ByPawPot, ZeffApprox, &
          &                          ByProjector, ReadFromPP
     use m_SpinOrbit_Potential,only: m_SO_dealloc_Dsoc, m_SO_dealloc_Mat_SOC_Strenth
! ============================= 2013/10/31

! ========= KT_add ========= 13.0U2
  use m_Control_Parameters,   only : sw_modified_TFW_functional
  use m_ThomasFermiW_Potential,  only : m_TFW_dealloc_ChgDensityBasisFn
! ========================== 13.0U2

! ====== KT_add === 2014/08/01
  use m_Orbital_QuantumNum, only : m_OP_Qnum_dealloc_array
! ================= 2014/08/01

  use m_Ldos, only : m_Ldos_dealloc

#ifdef FFTW3
  use m_FFT, only : m_FFT_finalize
#endif

     implicit none

     if (sw_rebuild_pws==OFF .or. driver == DRIVER_URAMP .or. driver == DRIVER_SC_DFT ) then
        call m_PP_dealloc()
        call m_NLP_dealloc()
        call m_Parallel_dealloc_mpi_nlmta()
        call m_Parallel_dealloc_mpi_exx()
        call m_Kp_dealloc
!        call m_ES_dealloc
!        call m_ESsd_dealloc
        return
     endif

     call m_CtrlP_dealloc
     call m_ES_dealloc
     call m_Parallel_cp_g1k(kv3)
     call m_Parallel_dealloc
     call m_CS_dealloc
     call m_IS_dealloc
  !call m_Parallel_dealloc_mpi_elec
     call m_Kp_dealloc
     call m_PP_dealloc
     call m_CD_dealloc
     call m_Force_dealloc
     call m_NLP_dealloc
     call m_ESsd_dealloc
     call m_XC_dealloc_vxc
     call m_pwBS_dealloc
     call m_PAW_dealloc
     if(flg_paw) then
         call m_PAW_XC_dealloc_vxc
     endif

! ==== KT_add ============== 2013/10/31
     if ( noncol ) then
        if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
             &                          .or. SpinOrbit_Mode == ReadFromPP ) then
           call m_SO_dealloc_Mat_SOC_strenth
           call m_SO_dealloc_Dsoc
        endif
        if ( SpinOrbit_Mode == ByProjector ) then
           call m_SO_dealloc_Dsoc
        endif
     endif
! ====================== 2013/10/31

! ====== KT_add ======== 13.0U2
     if ( sw_modified_TFW_functional /= OFF ) then
        call m_TFW_dealloc_ChgDensityBasisFn
     end if
! ====================== 13.0U2

! ====== KT_add =========== 13.0AS
     if (num_projectors>0)  call m_OP_dealloc
     if (sw_hubbard == ON) call m_ES_dealloc_Dhub
! ========================= 13.0AS

! ==== KT_add ==== 2014/08/01
     call m_OP_Qnum_dealloc_array
! ================ 2014/08/01

     call m_ES_wf_extrpl_dealloc()
     call m_Ldos_dealloc()
     call m_ES_wf_extrpl_dealloc()
#ifdef FFTW3
     call m_FFT_finalize()
#endif
     call m_CtrlP_set_init_status(.true.)
   end subroutine Array_Deallocate

   subroutine Write_Positron_LifeTime
     use m_Positron_Wave_Functions, only : m_pWF_wlifetime

     call m_pWF_wlifetime()
   end subroutine Write_Positron_LifeTime

end program PHASE
