#define POST3D
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
!  $Id: mdmain.F90 633 2020-12-01 05:11:03Z jkoga $
!
#ifdef FJ_TIMER
#  define __TIMER_FJ_START_w_BARRIER(str,a)  call mpi_barrier(str,ierr); call timer_sta(a)
#  define __TIMER_FJ_START(a)                call timer_sta(a)
#  define __TIMER_FJ_STOP(a)                 call timer_end(a)
#else
#  define __TIMER_FJ_START_w_BARRIER(str,a)
#  define __TIMER_FJ_START(a)
#  define __TIMER_FJ_STOP(a)
#endif
#ifdef __TIMER__
#  define __TIMER_START_w_BARRIER(str,a)  call mpi_barrier(str,ierr) ;   call timer_sta(a)
#  define __TIMER_START(a)                call timer_sta(a)
#  define __TIMER_STOP(a)                 call timer_end(a)
#else
#  define __TIMER_START_w_BARRIER(str,a)
#  define __TIMER_START(a)
#  define __TIMER_STOP(a)
#endif
!
program PHASE

#ifdef PARA3D
!!$  use m_Parallelization,     only : myrank_k_3D  , map_k_3D, np_e   &
!!$ &                                , ista_k       , iend_k           &
!!$ &                                , ista_kngp    , iend_kngp        &
!!$ &                                , ista_k_3D    , iend_k_3D        &
!!$ &                                , ista_kngp_3D , iend_kngp_3D     &
!!$ &                                , ista_snl, iend_snl &
!!$ &                                , np_fs_3D     , np_e_3D , np_g1k_3D
!!$  use m_Electronic_Structure,only : zaj_l, zaj_l_3D       &
!!$ &                                , eko_l, eko_l_3D       &
!!$ &                                , fsr_l, fsr_l_3D       &
!!$ &                                , fsi_l, fsi_l_3D       &
!!$ &                                , occup_l, occup_l_3D    &
!!$ &                                , vlhxc_l, vlhxc_l_3D   &
!!$ &                                , nrvf_ordr,neordr      &
!!$ &                                , vlhxcQ
!!$  use z_interface_3D, only        : decomp_eko_l_3D       &
!!$       &                          , decomp_eko_l_3D_2     &
!!$       &                          , decomp_eko_l_r_3D     &
!!$       &                          , decomp_eko_l_r_3D_2   &
!!$       &                          , decomp_fsr_l_3D       &
!!$       &                          , decomp_fsr_l_3D_ik    &
!!$       &                          , decomp_snl_l_3D       &
!!$       &                          , decomp_snl_l_r_3D     &
!!$       &                          , decomp_snl_l_3D_2     &
!!$       &                          , decomp_zaj_l_3D       &
!!$       &                          , decomp_zaj_l_3D_ik    &
!!$       &                          , decomp_fsr_l_r_3D     &
!!$       &                          , decomp_fsr_l_r_3D_ik  &
!!$       &                          , decomp_vnlph_l_r_3D   &
!!$       &                          , decomp_vnlph_l_3D     &
!!$       &                          , decomp_zaj_l_r_3D     &
!!$       &                          , decomp_zaj_l_r_3D_ik  &
!!$       &                          , decomp_vlhxc_l_3D     &
!!$       &                          , decomp_vlhxc_l_r_3D   &
!!$       &                          , decomp_wfsd_l_3D      &
!!$       &                          , decomp_wfsd_l_r_3D    &
!!$       &                          , decomp_zfm3_l_3D      &
!!$       &                          , decomp_zfm3_l_r_3D    &
!!$       &                          , decomp_psc_l_3D       &
!!$       &                          , decomp_psc_l_r_3D     &
!!$       &                          , decomp_qitg_l_3D      &
!!$       &                          , decomp_qitg_l_r_3D    &
!!$       &                          , decomp_gr_l_3D        &
!!$       &                          , decomp_gr_l_r_3D      &
!!$       &                          , decomp_occup_l_3D     &
!!$       &                          , decomp_occup_l_r_3D   &
!!$       &                          , decomp_ngpt_l_3D      &
!!$       &                          , decomp_ngpt_l_r_3D    &
!!$       &                          , decomp_igfp_l_3D      &
!!$       &                          , decomp_igfp_l_r_3D    &
!!$       &                          , replacement_zaj_ball_eigenvalue  &
!!$       &                          , replacement_zaj_ball_sequence
!!$  use m_Parallelization,    only  : mype
!!$  use m_IterationNumbers,   only  : iteration
!!$  use m_Control_Parameters,  only : nspin,kimg,neg,af,istress,sw_fine_STM_simulation, ON &
!!$       &                          , initial_chg, sw_positron
!!$  use m_Crystal_Structure,  only : nopr
!!$  use m_Const_Parameters,    only : GAMMA, from_PSEUDOPOTENTIAL_FILE, DEFECT, BULK
!!$  use m_PlaneWaveBasisSet,   only : kg1,kg,kgp           &
!!$  &                               , gr_l   , gr_l_3D     &
!!$  &                               , ngpt_l , ngpt_l_3D   &
!!$  &                               , igfp_l , igfp_l_3D
!!$  use m_Kpoints,             only : kv3, k_symmetry
!!$  use m_PseudoPotential,     only : nlmta   , nlmtt            &
!!$ &                                , qitg_l  , qitg_l_3D  , nqitg   &
!!$ &                                , psc_l   , psc_l_3D             &
!!$ &                                , rhpcg_l , rhpcg_l_3D , ntpcc   &
!!$ &                                , qitg_diff_l  , qitg_diff_l_3D  &
!!$ &                                , psc_diff_l   , psc_diff_l_3D   &
!!$ &                                , rhpcg_diff_l , rhpcg_diff_l_3D &
!!$ &                                , rhcg_l  , rhcg_l_3D   &
!!$ &                                , rhceg_l , rhceg_l_3D  &
!!$ &                                , rhchg_l , rhchg_l_3D  &
!!$ &                                , rhvg_l  , rhvg_l_3D
!!$  use m_NonLocal_Potential,  only : snl, snl_l_3D
!!$  use m_Charge_Density,      only : chgq_l , chgq_l_3D              &
!!$       &                          , chgqo_l, chgqo_l_3D             &
!!$       &                          , chgsoft, chgsoft_3D
!!$  use m_XC_Potential,        only : vxc_l  , vxc_l_3D               &
!!$       &                          , vxcpc_l, vxcpc_l_3D
!!$  use m_Ionic_System,        only : zfm3_l, zfm3_l_3D, ntyp
#endif
  use m_Files,                only : nfout
  use m_Ionic_System,         only : natm
  use m_PseudoPotential,      only : mmesh
  use m_Parallelization,      only : m_Parallel_init_mpi_paw_3D
  use m_Control_Parameters,   only : ipriparallel

  use m_Const_Parameters,    only : DRIVER_NEB,DRIVER_CONSTRAINT,DRIVER_MTD,DRIVER_SC_DFT, DRIVER_DIMER
  use m_Parallelization,     only : MPI_CommGroup

! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
use mod_timer
use m_Parallelization, only  : mype
#endif
! === TIMERTIMERTIMER ==========================================================

  implicit none
  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Already_Converged, Already_Converged2
  logical  :: Positron_bulk, Positron_defect, Structure_is_fixed
  logical  :: Hubbard_model
  logical  :: Forces_are_Converged, Ending_Time, Force_errors_are_tolerable,UnitCell_Converged &
  &         , Forces_are_Converged2
  logical  :: from_Initialize
  logical  :: InUeffRamping
!!$  logical  :: ChargeDensity_is_Fixed
   interface
     subroutine Initialization(init_mpi)
       integer, intent(in), optional :: init_mpi
     end subroutine Initialization
   end interface
#ifdef NEC_ITER_REG
  integer  :: count_for_ftrace
#endif

  integer :: ispin, ik, iksnl, ierR

  integer :: initmpi
  logical :: confpara
  integer :: driver
  logical :: uconv,ending_t,force_conv
  logical :: initialization_required

#ifdef NEC_ITER_REG
     count_for_ftrace = 0
     call FTRACE_REGION_BEGIN("INITIAL")
#endif

  confpara = Resolve_Config_Parallel()
  initmpi=1
  if(confpara) initmpi=0

  do

  ending_t = .false.
                                                  __TIMER_START(17)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,21)
!  call Initialization(init_mpi=1)
  call Initialization(initmpi)
                                                  __TIMER_FJ_STOP(21)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,22)
  call InputData_Analysis

  driver   = Resolve_Driver()

  if(ekmode()) then
    call Initialization_set_ekmode_ON()
  endif

  if (driver==DRIVER_NEB)then
#ifndef DISABLE_CONSTRAINTS
      call do_neb()
  else if (driver==DRIVER_DIMER) then
      call do_dimer_method()
  else if (driver==DRIVER_CONSTRAINT) then
      call constrained_dynamics()
  else if (driver==DRIVER_MTD) then
      call meta_dynamics()
#endif
  else
  outer_loop:do

!  if(driver==DRIVER_SC_DFT .and. .not.Initial_SCDFTLoop()) then
!     call InputData_Analysis()
!     if(ekmode()) call Initialization_set_ekmode_ON()
!  end if
                                                  __TIMER_FJ_STOP(22)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,23)
  call Preparation(0)                  ! Basis set, symmetry check etc.
                                                  __TIMER_FJ_STOP(23)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,24)
  if(ekmode()) call Preparation_for_mpi_ek()
  if(initialization_required()) then
    call Preparation_for_mpi(1)        ! mpi
  endif
                                                  __TIMER_FJ_STOP(24)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,25)
  call PseudoPotential_Construction
                                                  __TIMER_FJ_STOP(25)
  call  m_Parallel_init_mpi_paw_3D(nfout,ipriparallel,natm,mmesh)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,26)
#ifdef ENABLE_ESM_PACK
  if(initialization_required())then
    call Preparation_for_ESM
  endif
#endif

  call Ewald_and_Structure_Factor
                                                  __TIMER_FJ_STOP(26)
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,27)
  call Initial_Electronic_Structure
                                                  __TIMER_FJ_STOP(27)
                                                  __TIMER_STOP(17)

!!  from_Initialize = .false.
!$$  from_Initialize = .true.
!$$  if(from_Initialize) call Para3d_Entry()

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
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,16)
!!$                 call ReadCheckPointData_if_needed
                 call Renewal_of_WaveFunctions
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,11)
                 call ChargeDensity_Construction(1)
                 call ChargeDensity_Mixing
                                                  __TIMER_STOP(11)
                 ending_t = Ending_Time()
                 if(ending_t)                      exit outer_loop
                 if(TotalEnergy_is_Divergent())    exit outer_loop
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,11)
                 call Renewal_of_Potential
                                                  __TIMER_STOP(11)
                 if(Hubbard_model()) then
                    call Renewal_of_Hubbard_Potential
                 end if
!!$                 call WriteCheckPointData  ! if necessary
                 if(ChargeDensity_is_Converged()) then
                                                  __TIMER_STOP(16)
                    exit ChargeDensity
                 end if
                                                  __TIMER_STOP(16)
              enddo ChargeDensity

!$$!              call PARA_3Dto2D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!

              if(Structure_is_fixed()) then
                  force_conv = .true.
                  exit StressLoop
              end if
!$$!              call PARA_2Dto3D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!
#ifdef LIBRARY_BUILD
              call Forces_phase0
#else
              call Forces
#endif
!$$!              call PARA_3Dto2D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!
              force_conv = Forces_are_Converged()
              if ( Enforce_Update_Cell() )  exit AtomicConfiguration

              if(force_conv) then
                 if ( epsilon_post_scf() ) call Epsilon_PostScf
                 exit AtomicConfiguration
              endif
              if(Force_errors_are_tolerable()) then
                 call Postprocessing_during_MD()
                 call Move_Ions
!$$!              call PARA_2Dto3D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!
!$$!              call PARA_3Dto2D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!
                 call MDIterationNumber_Setting
                 call Ewald_and_Structure_Factor
              end if

              if(BreakMD(force_conv))then
                 exit AtomicConfiguration
              endif
              call Postproc_after_SCF_convergence()
           enddo AtomicConfiguration
           if(driver==DRIVER_SC_DFT) exit StressLoop
           call Stress
           exit StressLoop
        end do StressLoop

!$$!        call PARA_3Dto2D()     !!!!!!!!!! 3D → 2D !!!!!!!!!!


#ifdef NEC_ITER_REG
        if(count_for_ftrace .eq. 1) then
           call FTRACE_REGION_END("SOLVE-FIRST")
        else
           call FTRACE_REGION_END("SOLVE-CORE")
        end if
#endif
     end if
     if (InUeffRamping().and. .not.ending_t)then
        call MDiterationNumber_Setting3()
        if(.not.ending_t) call Array_Deallocate()
     else if (driver == DRIVER_SC_DFT.and. .not.ending_t) then
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
     else
!!$#ifdef POST3D
!!$        if(.not.InUeffRamping()) call Postprocessing(.false.)
!!$#endif
        exit outer_loop
     endif

  end if

  enddo outer_loop
  uconv = force_conv
  if (.not. ending_t) &
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
#ifdef POST3D
       if(uconv) call Postprocessing(.false.)
#endif
!!$        if(driver/=DRIVER_SC_DFT) then
!!$!BRANCH_P ORG_Parallel
!!$! * T.Hamada 2016.1.8
!!$           call Epsilon_Paramset
!!$           call Epsilon_Postscf
!!$        end if
!!$! * T.Hamada 2015.1.8
!!$!BRANCH_P_END ORG_Parallel
        if ( uconv .and. epsilon_post_proc() ) call Epsilon_PostScf
        call rttddft_main
!        call WriteDownData_onto_Files(ending_t.or.uconv)
        if(ending_t .or. uconv) call final_output()
        if(uconv) exit
     end if
  endif
  if(.not.ending_t) call MDiterationNumber_Setting_pre()

#ifdef NEC_ITER_REG
     call FTRACE_REGION_BEGIN("FINAL")
#endif

!$$!        call PARA_2Dto3D()
                                                  __TIMER_FJ_START_w_BARRIER(MPI_CommGroup,39)
  if(.not.ending_t) then
     call Array_Deallocate()
     call Continuation_Mode()
  endif

  endif

  if (ending_t .or. OneShot()) exit

  enddo

  if(driver/=DRIVER_NEB) call Finalization_of_mpi           ! mpi
#ifdef NEC_ITER_REG
  call FTRACE_REGION_END("FINAL")
#endif
! === TIMERTIMERTIMER ==========================================================
#ifdef __TIMER_FFT__
call print_timer(mype)
#endif
! === TIMERTIMERTIMER ==========================================================
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
      call Initialization_Epsilon(0)        ! Epsilon
      call Shift_Kpoint
    endif
    KPOINT_GROUP: do
       call KpointNumber_Setting2()
       if(AllKpoints_are_Converged()) exit KPOINT_GROUP
       call Preparation_ek()
       call Preparation_for_mpi_ek
       call epsmain_reallocate()
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
       call Finalization_ek()
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

  subroutine Finalization_ek()
    use m_Const_Parameters,   only : ON
    use m_Control_Parameters, only : sw_dos, sw_ldos
    use m_Ldos,               only : m_Ldos_dealloc_weiwsc_etc
    implicit none
!    if(sw_dos == ON .and. sw_ldos == ON) then
!      call m_Ldos_dealloc_weiwsc_etc()
!    endif
  end subroutine Finalization_ek

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

#ifdef PARA3D
!!$  subroutine Para3d_Entry()
!!$    if (allocated(vlhxc_l))    deallocate(vlhxc_l)
!!$    if (allocated(chgq_l))     deallocate(chgq_l)
!!$    if (allocated(chgqo_l))    deallocate(chgqo_l)
!!$    if(istress == ON .or. sw_fine_STM_simulation == ON) then
!!$       if (allocated(chgsoft)) deallocate(chgsoft)
!!$    end if
!!$    if (allocated(vxc_l))      deallocate(vxc_l)
!!$    if (allocated(vxcpc_l))    deallocate(vxcpc_l)
!!$    if (allocated(snl))        deallocate(snl)
!!$    if (allocated(eko_l))      deallocate(eko_l)
!!$    if (allocated(zaj_l))      deallocate(zaj_l)
!!$    if (allocated(fsr_l))      deallocate(fsr_l)
!!$    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!!$       if (allocated(fsi_l))   deallocate(fsi_l)
!!$    end if
!!$    if (allocated(occup_l))    deallocate(occup_l)
!!$    if (allocated(zfm3_l))     deallocate(zfm3_l)
!!$    if (allocated(gr_l))       deallocate(gr_l)
!!$    if (allocated(psc_l))      deallocate(psc_l)
!!$    if (allocated(qitg_l))     deallocate(qitg_l)
!!$    if (allocated(rhpcg_l))    deallocate(rhpcg_l)
!!$    if(istress==ON) then
!!$       if (allocated(psc_diff_l))      deallocate(psc_diff_l)
!!$       if (allocated(qitg_diff_l))     deallocate(qitg_diff_l)
!!$       if (allocated(rhpcg_diff_l))    deallocate(rhpcg_diff_l)
!!$    endif
!!$    if (allocated(ngpt_l))     deallocate(ngpt_l)
!!$    if (allocated(igfp_l))    deallocate(igfp_l)
!!$  end subroutine Para3d_Entry
!!$
!!$  subroutine PARA_2Dto3D()
!!$    !                timer_fj sta 37
!!$
!!$    allocate(vlhxc_l_3D(ista_kngp_3D:iend_kngp_3D,kimg,nspin));             vlhxc_l_3D = 0.0d0
!!$    do ispin = 1, nspin, (af+1)
!!$       call decomp_vlhxc_l_3D(vlhxc_l,vlhxc_l_3D,ispin)
!!$    end do
!!$    deallocate(vlhxc_l)
!!$
!!$    allocate(chgq_l_3D(ista_kngp_3D:iend_kngp_3D,kimg,nspin));              chgq_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       call decomp_vlhxc_l_3D(chgq_l,chgq_l_3D,ispin)
!!$    end do
!!$    deallocate(chgq_l)
!!$
!!$    allocate(chgqo_l_3D(ista_kngp_3D:iend_kngp_3D,kimg,nspin));             chgqo_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       call decomp_vlhxc_l_3D(chgqo_l,chgqo_l_3D,ispin)
!!$    end do
!!$    deallocate(chgqo_l)
!!$
!!$    if(istress == ON .or. sw_fine_STM_simulation == ON) then
!!$       allocate(chgsoft_3D(ista_kngp_3D:iend_kngp_3D,kimg,nspin));          chgsoft_3D = 0.0d0
!!$       do ispin = 1, nspin, (af+1)
!!$          call decomp_vlhxc_l_3D(chgsoft,chgsoft_3D,ispin)
!!$       end do
!!$       deallocate(chgsoft)
!!$    end if
!!$
!!$    allocate(vxc_l_3D(ista_kngp_3D:iend_kngp_3D,kimg,nspin));               vxc_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       call decomp_vlhxc_l_3D(vxc_l,vxc_l_3D,ispin)
!!$    end do
!!$    deallocate(vxc_l)
!!$
!!$    allocate(vxcpc_l_3D(ista_kngp_3D:iend_kngp_3D,kimg));                   vxcpc_l_3D = 0.d0
!!$    deallocate(vxcpc_l)
!!$
!!$    allocate(snl_l_3D(maxval(np_g1k_3D),nlmtt,ista_snl:iend_snl));          snl_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$              call decomp_snl_l_3D_2(snl,snl_l_3D,ik)
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(snl)
!!$
!!$    allocate(eko_l_3D(np_e_3D,ista_k_3D:iend_k_3D));                        eko_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$              call decomp_eko_l_3D_2(eko_l,eko_l_3D,ik,nrvf_ordr,'sort')
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(eko_l)
!!$
!!$    allocate(zaj_l_3D(maxval(np_g1k_3D),np_e_3D,ista_k_3D:iend_k_3D,kimg)); zaj_l_3D = 0.0d0
!!$    do ispin = 1, nspin, (af+1)
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$              call decomp_zaj_l_3D_ik(zaj_l,zaj_l_3D,ik,nrvf_ordr,"sort")
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(zaj_l)
!!$
!!$    allocate(fsr_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D));               fsr_l_3D = 0.0d0
!!$    do ispin = 1, nspin, (af+1)
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$              call decomp_fsr_l_3D_ik(fsr_l,fsr_l_3D,ik,nrvf_ordr,'sort')
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(fsr_l)
!!$
!!$    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!!$       allocate(fsi_l_3D(np_e_3D,np_fs_3D,ista_k_3D:iend_k_3D));            fsi_l_3D = 0.d0
!!$       do ispin = 1, nspin, (af+1)
!!$          do ik = ispin, kv3-nspin+ispin, nspin
!!$             if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$                call decomp_fsr_l_3D_ik(fsi_l,fsi_l_3D,ik,nrvf_ordr,'sort')
!!$             end if
!!$          end do
!!$       end do
!!$       deallocate(fsi_l)
!!$    end if
!!$
!!$    allocate(occup_l_3D(np_e_3D,ista_k_3D:iend_k_3D));                      occup_l_3D = 0.d0
!!$    do ispin = 1, nspin, (af+1)
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$              call decomp_occup_l_3D(occup_l,occup_l_3D,ik,nrvf_ordr,'sort')
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(occup_l)
!!$
!!$    allocate(zfm3_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp,kimg));               zfm3_l_3D = 0.d0
!!$    call decomp_zfm3_l_3D(zfm3_l,zfm3_l_3D)
!!$    deallocate(zfm3_l)
!!$
!!$    allocate(gr_l_3D(ista_kngp_3D:iend_kngp_3D));                           gr_l_3D =0.0d0
!!$    call decomp_gr_l_3D(gr_l,gr_l_3D)
!!$    deallocate(gr_l)
!!$
!!$    allocate(psc_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));                     psc_l_3D =0.0d0
!!$    call decomp_psc_l_3D(psc_l,psc_l_3D)
!!$    deallocate(psc_l)
!!$
!!$    allocate(qitg_l_3D(ista_kngp_3D:iend_kngp_3D,nqitg));                   qitg_l_3D =0.0d0
!!$    call decomp_qitg_l_3D(qitg_l,qitg_l_3D,nqitg)
!!$    deallocate(qitg_l)
!!$
!!$    allocate(rhpcg_l_3D(ista_kngp_3D:iend_kngp_3D,ntpcc));                  rhpcg_l_3D =0.0d0
!!$    call decomp_qitg_l_3D(rhpcg_l,rhpcg_l_3D,ntpcc)
!!$    deallocate(rhpcg_l)
!!$
!!$    if(istress==ON) then
!!$       allocate(psc_diff_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));              psc_diff_l_3D =0.0d0
!!$       call decomp_psc_l_3D(psc_diff_l,psc_diff_l_3D)
!!$       deallocate(psc_diff_l)
!!$
!!$       allocate(qitg_diff_l_3D(ista_kngp_3D:iend_kngp_3D,nqitg));            qitg_diff_l_3D =0.0d0
!!$       call decomp_qitg_l_3D(qitg_diff_l,qitg_diff_l_3D,nqitg)
!!$       deallocate(qitg_diff_l)
!!$
!!$       allocate(rhpcg_diff_l_3D(ista_kngp_3D:iend_kngp_3D,ntpcc));           rhpcg_diff_l_3D =0.0d0
!!$       call decomp_qitg_l_3D(rhpcg_diff_l,rhpcg_diff_l_3D,ntpcc)
!!$       deallocate(rhpcg_diff_l)
!!$    endif
!!$
!!$    allocate(ngpt_l_3D(ista_kngp_3D:iend_kngp_3D,nopr+af));         ngpt_l_3D =0.0d0
!!$    call decomp_ngpt_l_3D(ngpt_l,ngpt_l_3D,nopr+af)
!!$    deallocate(ngpt_l)
!!$
!!$    allocate(igfp_l_3D(ista_kngp_3D:iend_kngp_3D));                 igfp_l_3D =0.0d0
!!$    call decomp_igfp_l_3D(igfp_l,igfp_l_3D)
!!$    deallocate(igfp_l)
!!$
!!$   if(initial_chg == from_PSEUDOPOTENTIAL_FILE) then
!!$      allocate(rhvg_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));                     rhvg_l_3D =0.0d0
!!$      call decomp_psc_l_3D(rhvg_l,rhvg_l_3D)
!!$      deallocate(rhvg_l)
!!$   endif
!!$
!!$   if(sw_positron == BULK .or. sw_positron == DEFECT) then
!!$      allocate(rhcg_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));                     rhcg_l_3D =0.0d0
!!$      allocate(rhceg_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));                    rhceg_l_3D =0.0d0
!!$      allocate(rhchg_l_3D(ista_kngp_3D:iend_kngp_3D,ntyp));                    rhchg_l_3D =0.0d0
!!$      call decomp_psc_l_3D(rhcg_l,rhcg_l_3D)
!!$      call decomp_psc_l_3D(rhceg_l,rhceg_l_3D)
!!$      call decomp_psc_l_3D(rhchg_l,rhchg_l_3D)
!!$      deallocate(rhcg_l)
!!$      deallocate(rhceg_l)
!!$      deallocate(rhchg_l)
!!$   endif
!!$
!!$   !                    timer_fj end 37
!!$  end subroutine PARA_2Dto3D
!!$
!!$  subroutine PARA_3Dto2D()
!!$    !                   timer_fj sta 38
!!$
!!$    allocate(vlhxc_l(ista_kngp:iend_kngp,kimg,nspin)); vlhxc_l = 0.0d0
!!$    do ispin = 1, nspin, af + 1
!!$       call decomp_vlhxc_l_r_3D(vlhxc_l,vlhxc_l_3D,ispin)
!!$    end do
!!$    deallocate(vlhxc_l_3D)
!!$
!!$    allocate(chgq_l(ista_kngp:iend_kngp,kimg,nspin)); chgq_l = 0.d0
!!$    do ispin = 1, nspin, af + 1
!!$       call decomp_vlhxc_l_r_3D(chgq_l,chgq_l_3D,ispin)
!!$    end do
!!$    deallocate(chgq_l_3D)
!!$
!!$    allocate(chgqo_l(ista_kngp:iend_kngp,kimg,nspin)); chgqo_l = 0.d0
!!$    do ispin = 1, nspin, af + 1
!!$       call decomp_vlhxc_l_r_3D(chgqo_l,chgqo_l_3D,ispin)
!!$    end do
!!$    deallocate(chgqo_l_3D)
!!$
!!$    if(istress == ON .or. sw_fine_STM_simulation == ON) then
!!$       allocate(chgsoft(ista_kngp:iend_kngp,kimg,nspin)) ; chgsoft = 0.0d0
!!$       do ispin = 1, nspin, af + 1
!!$          call decomp_vlhxc_l_r_3D(chgsoft,chgsoft_3D,ispin)
!!$       end do
!!$       deallocate(chgsoft_3D)
!!$    end if
!!$
!!$    allocate(vxc_l(ista_kngp:iend_kngp,kimg,nspin)); vxc_l = 0.d0
!!$    do ispin = 1, nspin, af + 1
!!$       call decomp_vlhxc_l_r_3D(vxc_l,vxc_l_3D,ispin)
!!$    end do
!!$    deallocate(vxc_l_3D)
!!$
!!$    allocate(vxcpc_l(ista_kngp:iend_kngp,kimg)); vxcpc_l = 0.d0
!!$    deallocate(vxcpc_l_3D)
!!$
!!$    allocate(snl(kg1,nlmtt,ista_snl:iend_snl)); snl = 0.d0
!!$    call decomp_snl_l_r_3D(snl_l_3D,snl)
!!$    deallocate(snl_l_3D)
!!$
!!$    allocate(eko_l(np_e,ista_k:iend_k));                   eko_l = 0.d0
!!$    do ispin = 1, nspin, af + 1
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$             call decomp_eko_l_r_3D_2(eko_l,eko_l_3D,ik,neordr,'sort')
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(eko_l_3D)
!!$
!!$    allocate(zaj_l(kg1,np_e,ista_k:iend_k,kimg));          zaj_l = 0.0d0
!!$    do ispin = 1, nspin, af + 1
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$             call decomp_zaj_l_r_3D_ik(zaj_l,zaj_l_3D,ik,neordr,"sort")
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(zaj_l_3D)
!!$
!!$    allocate(fsr_l(np_e,nlmta,ista_k:iend_k));             fsr_l = 0.0d0
!!$    do ispin = 1, nspin, af + 1
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$             call decomp_fsr_l_r_3D_ik(fsr_l,fsr_l_3D,ik,neordr,"sort",0)
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(fsr_l_3D)
!!$
!!$    if(.not.(kv3/nspin == 1 .and. k_symmetry(1) == GAMMA .and. kimg == 2)) then
!!$       allocate(fsi_l(np_e,nlmta,ista_k:iend_k));          fsi_l = 0.d0
!!$       do ispin = 1, nspin, af + 1
!!$          do ik = ispin, kv3-nspin+ispin, nspin
!!$             if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$                call decomp_fsr_l_r_3D_ik(fsi_l,fsi_l_3D,ik,neordr,"sort",0)
!!$             end if
!!$          end do
!!$       end do
!!$       deallocate(fsi_l_3D)
!!$    end if
!!$
!!$    allocate(occup_l(np_e,ista_k:iend_k));                 occup_l = 0.d0
!!$    do ispin = 1, nspin, af + 1
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$             iksnl = (ik-1)/nspin + 1
!!$             call decomp_occup_l_r_3D(occup_l,occup_l_3D,ik,neordr,'sort')
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(occup_l_3D)
!!$
!!$    allocate(zfm3_l(ista_kngp:iend_kngp,ntyp,kimg));       zfm3_l = 0.d0
!!$    call decomp_zfm3_l_r_3D(zfm3_l,zfm3_l_3D)
!!$    deallocate(zfm3_l_3D)
!!$
!!$    allocate(gr_l(ista_kngp:iend_kngp));                   gr_l =0.0d0
!!$    call decomp_gr_l_r_3D(gr_l,gr_l_3D)
!!$    deallocate(gr_l_3D)
!!$
!!$    allocate(psc_l(ista_kngp:iend_kngp,ntyp));             psc_l =0.0d0
!!$    call decomp_psc_l_r_3D(psc_l,psc_l_3D)
!!$    deallocate(psc_l_3D)
!!$
!!$    allocate(qitg_l(ista_kngp:iend_kngp,nqitg));           qitg_l =0.0d0
!!$    call decomp_qitg_l_r_3D(qitg_l,qitg_l_3D,nqitg)
!!$    deallocate(qitg_l_3D)
!!$
!!$    allocate(rhpcg_l(ista_kngp:iend_kngp,ntpcc));          rhpcg_l =0.0d0
!!$    call decomp_qitg_l_r_3D(rhpcg_l,rhpcg_l_3D,ntpcc)
!!$    deallocate(rhpcg_l_3D)
!!$
!!$    if(istress==ON) then
!!$       allocate(psc_diff_l(ista_kngp:iend_kngp,ntyp));             psc_diff_l =0.0d0
!!$       call decomp_psc_l_r_3D(psc_diff_l,psc_diff_l_3D)
!!$       deallocate(psc_diff_l_3D)
!!$
!!$       allocate(qitg_diff_l(ista_kngp:iend_kngp,nqitg));           qitg_diff_l =0.0d0
!!$       call decomp_qitg_l_r_3D(qitg_diff_l,qitg_diff_l_3D,nqitg)
!!$       deallocate(qitg_diff_l_3D)
!!$
!!$       allocate(rhpcg_diff_l(ista_kngp:iend_kngp,ntpcc));          rhpcg_diff_l =0.0d0
!!$       call decomp_qitg_l_r_3D(rhpcg_diff_l,rhpcg_diff_l_3D,ntpcc)
!!$       deallocate(rhpcg_diff_l_3D)
!!$    endif
!!$
!!$    allocate(ngpt_l(ista_kngp:iend_kngp,nopr+af));         ngpt_l =0.0d0
!!$    call decomp_ngpt_l_r_3D(ngpt_l,ngpt_l_3D,nopr+af)
!!$    deallocate(ngpt_l_3D)
!!$
!!$    allocate(igfp_l(ista_kngp:iend_kngp));                 igfp_l =0.0d0
!!$    call decomp_igfp_l_r_3D(igfp_l,igfp_l_3D)
!!$    deallocate(igfp_l_3D)
!!$
!!$   if(initial_chg == from_PSEUDOPOTENTIAL_FILE) then
!!$      allocate(rhvg_l(ista_kngp:iend_kngp,ntyp));                     rhvg_l =0.0d0
!!$      call decomp_psc_l_r_3D(rhvg_l,rhvg_l_3D)
!!$      deallocate(rhvg_l_3D)
!!$   endif
!!$
!!$   if(sw_positron == BULK .or. sw_positron == DEFECT) then
!!$      allocate(rhcg_l(ista_kngp:iend_kngp,ntyp));                     rhcg_l =0.0d0
!!$      allocate(rhceg_l(ista_kngp:iend_kngp,ntyp));                    rhceg_l =0.0d0
!!$      allocate(rhchg_l(ista_kngp:iend_kngp,ntyp));                    rhchg_l =0.0d0
!!$      call decomp_psc_l_r_3D(rhcg_l,rhcg_l_3D)
!!$      call decomp_psc_l_r_3D(rhceg_l,rhceg_l_3D)
!!$      call decomp_psc_l_r_3D(rhchg_l,rhchg_l_3D)
!!$      deallocate(rhcg_l_3D)
!!$      deallocate(rhceg_l_3D)
!!$      deallocate(rhchg_l_3D)
!!$   endif
!!$
!!$   !                timer_fj end 38
!!$  end subroutine PARA_3Dto2D
!!$
!!$  subroutine output_zaj()
!!$   integer::ispin, ik
!!$
!!$    if (.not. allocated(zaj_l)) then
!!$       allocate(zaj_l(kg1,np_e,ista_k:iend_k,kimg))
!!$    end if
!!$    zaj_l = 0.0d0
!!$
!!$    do ispin = 1, nspin, af + 1
!!$       do ik = ispin, kv3-nspin+ispin, nspin
!!$          if(map_k_3D(ik) == myrank_k_3D) then           ! MPI
!!$             call decomp_zaj_l_r_3D_ik(zaj_l,zaj_l_3D,ik,neordr,"sort")
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    call dump_zaj(zaj_l)
!!$
!!$    deallocate(zaj_l)
!!$
!!$  endsubroutine output_zaj
#endif

  subroutine Continuation_Mode()
    use m_Const_Parameters,   only : COORDINATE_CONTINUATION,ON,OFF
    use m_Control_Parameters, only : icond, sw_optimize_lattice,sw_rebuild_pws
    implicit none
    logical :: unitcell_can_change
    if(unitcell_can_change() .and. sw_rebuild_pws==OFF) return
    icond = COORDINATE_CONTINUATION
  end subroutine continuation_mode

  logical function Enforce_Update_Cell()
    use m_IterationNumbers,     only : iteration_unit_cell
    use m_Control_Parameters,   only : sw_optimize_coords_sametime, &
         &                             threshold_start_cellopt
    use m_Const_Parameters,     only : ON
    use m_Force,                only : forcmx

    Enforce_update_Cell = .false.
    if ( sw_optimize_coords_sametime == ON ) then
       if ( iteration_unit_cell > 1 ) then
          Enforce_update_Cell = .true.
       else
          if ( forcmx < threshold_start_cellopt ) then
             Enforce_update_Cell = .true.
          endif
       endif
    endif

  end function Enforce_Update_Cell

  logical function BreakMD(conv)
    use m_IterationNumbers, only : iteration_ionic
    use m_Ionic_System, only : addition_frequency,m_IS_natm_can_change
    implicit none
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
    implicit none
    logical :: unitcell_can_change
    OneShot = .not.m_IS_natm_can_change()
    if(OneShot) then
        OneShot = .not.unitcell_can_change()
    endif
  end function OneShot

! === KT_add ==== 13.1R
  logical function epsilon_post_scf()
    use m_Const_Parameters,   only : ON, OFF
    use m_Epsilon_ek,  only : sw_epsilon, m_Eps_chkif_sw_epsilon
    use m_Control_Parameters, only : sw_phonon, sw_phonon_with_epsilon, &
         &                           sw_calc_dielectric_tensor, sw_excitation

    implicit none
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

    implicit none
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

  subroutine Postprocessing_during_MD()
    use m_IterationNumbers, only : iteration_ionic
    use m_Control_Parameters, only : postproc_frequency
    if(postproc_frequency<=0) return
    if(mod(iteration_ionic,postproc_frequency)==0)then
       call Postprocessing(.true.)
    endif
  end subroutine Postprocessing_during_MD

   subroutine Set_Icond
     use m_Const_Parameters,   only : INITIAL
     use m_Control_Parameters, only : m_CtrlP_set_icond
     call m_CtrlP_set_icond(INITIAL)
   end subroutine Set_Icond

end program PHASE

  subroutine Array_Deallocate()
     use m_Const_Parameters, only : OFF, ON, DRIVER_URAMP, DRIVER_SC_DFT, DRIVER_NEB
     use m_Control_Parameters, only : m_CtrlP_dealloc,sw_rebuild_pws,m_CtrlP_set_init_status
     use m_Crystal_Structure, only : m_CS_dealloc
     use m_Ionic_System, only : m_IS_dealloc
     use m_PlaneWaveBasisSet, only : m_pwBS_dealloc
     use m_Parallelization, only : m_Parallel_dealloc,m_Parallel_dealloc_mpi_nlmta,m_Parallel_dealloc_mpi_elec &
    &                            , m_Parallel_dealloc_mpi_nval,m_Parallel_dealloc_mpi_fft_box &
    &                            , m_Parallel_dealloc_mpi_kngp_B,m_Parallel_fft_onto_wf_dealloc_3D &
    &                            , m_Parallel_dealloc_mpi_paw_3D, m_Parallel_dealloc_mpi_exx &
    &                            , m_Parallel_cp_g1k
     use m_Kpoints, only : m_Kp_dealloc, kv3
     use m_Force, only : m_Force_dealloc
     use m_Charge_Density, only : m_CD_dealloc
     use m_XC_Potential, only : m_XC_dealloc_vxc_3D
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
     use m_Control_Parameters, only: noncol, SpinOrbit_Mode, driver
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
        call m_Parallel_dealloc_mpi_paw_3D()
!        call m_ES_dealloc
!        call m_ESsd_dealloc
        return
     endif

     call m_CtrlP_dealloc
     call m_ES_dealloc
     call m_pwBS_dealloc
     call m_Parallel_cp_g1k(kv3)
     call m_Parallel_dealloc
     call m_Parallel_dealloc_mpi_elec
     call m_Parallel_dealloc_mpi_fft_box
     call m_Parallel_dealloc_mpi_nval
     call m_Parallel_fft_onto_wf_dealloc_3D
     call m_Parallel_dealloc_mpi_kngp_B
     call m_CS_dealloc
     call m_IS_dealloc
  !call m_Parallel_dealloc_mpi_elec
     call m_Kp_dealloc
     call m_PP_dealloc
     call m_CD_dealloc
     call m_Force_dealloc
     call m_NLP_dealloc
     call m_ESsd_dealloc
     call m_XC_dealloc_vxc_3D
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

