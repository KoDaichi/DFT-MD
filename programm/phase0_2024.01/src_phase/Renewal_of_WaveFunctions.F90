!#define _DEBUG_WRITE_
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 603 $)
!
!  SUBROUINE: Renewal_of_WaveFunctions, Renew_WF_by_lmSDorlmCG,
!             check_new_solver
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
!                        J. Koga,  March/01/2010
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
#   define __TIMER_START_w_BARRIER(str,a)
#   define __TIMER_START(a)
#   define __TIMER_STOP(a)
subroutine Renewal_of_WaveFunctions()
! $Id: Renewal_of_WaveFunctions.F90 603 2020-04-07 03:25:30Z jkoga $
  use m_Control_Parameters,   only : m_CtrlP_solver_for_WFs_now &
       &                           , m_CtrlP_On_or_Off_precon_WFs &
       &                           , m_CtrlP_dtim_now &
       &                           , m_CtrlP_clear_nsolver_applied &
       &                           , m_CtrlP_push_SolverNameApplied &
       &                           , m_CtrlP_set_sw_MRCV_only &
       &                           , meg,damp,icond,iprievdff &
       &                           , intzaj,ipri,ipriwf,renew_wf_again_m_CtrlP &
       &                           , sw_ekzaj,ekmode,sw_hubbard &
       &                           , submat_before_renewal, sw_rsb &
       &                           , sw_hybrid_functional
  use m_Const_Parameters,     only : DP,SD,MSD,CG &
       &                           , lmSD, lmMSD, lmCG, lmeazyCG & ! , eazyCG
       &                           , RMM,RMM2,RMM2P,SUBMAT,MATRIXDIAGON &
       &                           , DAVIDSON, MDDAVIDSON, MDKOSUGI &
       &                           , YES, ON, OFF, NO &
       &                           , INITIAL, CONTINUATION, FIXED_CHARGE &
       &                           , FIXED_CHARGE_CONTINUATION &
       &                           , COORDINATE_CONTINUATION,INVERSE,DIRECT
  use m_IterationNumbers,     only : iteration, iteration_electronic, iteration_ionic &
       &                           , m_Iter_set_rmm_start, nk_in_the_process
  use m_Files,                only : nfout
  use m_Ionic_System,         only : natm2
  use m_ES_WF_by_SDorCG,      only : m_ESsd_alloc_wfsd, m_ESsd_dealloc_wfsd
  use m_ES_WF_by_RMM,         only : m_ESrmm_reset_ng_maxmin &
       &                           , m_ESrmm_reset_r_norm_flag
  use m_Electronic_Structure, only : m_ES_what_is_evdff_now,zaj_l
  use m_ES_IO,                only : m_ESIO_wd_WFs_standardout

  use m_Total_Energy,         only : m_TE_what_is_edeltb_now
  use m_Charge_Density,       only : m_CD_cp_chgq_to_chgqo
  use m_Orbital_Population,   only : m_OP_cp_ommix_to_omold
  use m_ES_WF_by_submat,      only : m_ESsubmat_Renew_WF, m_ESsubmat_alloc, m_ESsubmat_utransform_wf, disable_utransform
#ifdef __TIMER__
  use m_Parallelization,      only : MPI_CommGroup
#endif

#ifndef DISABLE_NONCL
! ================================ added by K. Tagami ================ 11.0
  use m_Control_Parameters,    only : noncol
  use m_ES_WF_by_submat,          only : m_ESsubmat_Renew_WF_noncl
! =================================================================== 11.0
#endif
! === KT_add ====== 13.0U3
  use m_Control_Parameters,   only : precon_hardpart, sw_wf_mixing, wf_mixing_is_active
  use m_ES_WF_mixing,         only : m_ES_WF_cp_zaj_to_zaj_prev
! ================= 13.0U3

  use m_ES_RSB, only : m_ES_RSB_doit,m_ES_RSB_alloc,m_ES_RSB_dealloc
  use m_Kpoints, only : kv3
  use m_ES_ExactExchange, only : m_ES_EXX_gather_valence_states
  use m_Electronic_Structure, only : m_ES_sort_eigen_values

  implicit none
  integer       :: what_is_the_core_solver
  integer       :: isolver, precon, previous_isolver = -100, sw_submat, isolver_core
  integer, save :: is_wfsd_allocated = NO
  real(kind=DP) :: dtim, edeltb_now
#ifdef __TIMER__
  integer :: ierr
#endif
  integer :: ik
  logical :: push_solver=.true.
  logical :: hybfunc=.false.

  !! For EXX test
  !! call m_ES_EXX_update(OFF)

  call m_CtrlP_clear_nsolver_applied()

  if(icond == INITIAL .or. icond == CONTINUATION .or. icond == COORDINATE_CONTINUATION .or. &
       & ((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION)&
       &                             .and. ekmode == OFF)   ) then
     edeltb_now = m_TE_what_is_edeltb_now()/natm2
!!$     write(nfout,'(" edeltb_per_atom = ", d20.12)') edeltb_per_atom
  else
     edeltb_now = m_ES_what_is_evdff_now()
  end if

  if(iprievdff >= 2) write(nfout,'(" edeltb_now = ",d20.10)') edeltb_now

  isolver = m_CtrlP_solver_for_WFs_now(iteration_electronic,iteration_ionic,intzaj &
       &                              ,edeltb_now, sw_submat)
  call check_new_solver  ! -(contained here) comparing isolver with previous_isolver
  !                        Especially, important for the RMM operation.
!!$  previous_isolver = isolver
! === Merge modifications from phase@63 -> phase@64. 2011/10/18 ================
  precon  = m_CtrlP_On_or_Off_precon_WFs(iteration_electronic,iteration_ionic,isolver)
! ==============================================================================
  dtim = m_CtrlP_dtim_now(iteration_electronic,iteration_ionic)
  if(iprievdff >= 2) write(nfout,'(" iter_elec, dtim = ",i6,f10.4)') &
       & iteration_electronic, dtim

! === Merge modifications from phase@63 -> phase@64. 2011/10/18 ================
#ifdef _DEBUG_WRITE_
  if(ipri>=1) write(nfout,'(" !isolver, precon = ",2i8)') isolver, precon
#endif
! ==============================================================================

  isolver_core = what_is_the_core_solver(isolver) ! -(in this file)
  if(isolver_core == CG .and. is_wfsd_allocated == NO) then
     call m_ESsd_alloc_wfsd()
     is_wfsd_allocated = YES
  else if(isolver_core /= CG .and. is_wfsd_allocated == YES) then
     call m_ESsd_dealloc_wfsd()
     is_wfsd_allocated = NO
  end if

  !! for calc. of Berry phase at the Gamma point only
  !! providing that well converged wavefunctions are supplied from nfzaj file
  if(sw_ekzaj == ON .and. ekmode == ON) then
     if(ipri >= 1) write(nfout,*) &
	& '** Returning from Renewal_of_WaveFunctions for Berry phase calculation with samping the Gamma point only **'
     return
  end if

!!$   select case(isolver)
!!$     case (lmSD, lmMSD, lmCG, lmeazyCG)
!!$        if(iteration_electronic >= 10) then
!!$           if(incre_etot_in_1dsrch >= incre_etot_in_1dsrch_limit) then
!!$              call ReadCheckPointData()
!!$           else if (flag_wdcheckpointdata >= 0)
!!$           end if
!!$        end if
!!$     end select

!!!!BRANCH_P ORG_Parallel
!!!  if(sw_rsb==ON)then
!!!    call m_ES_RSB_alloc()
!!!    call m_ES_RSB_doit(verbose=.false.)
!!!  endif
!!!!BRANCH_P_END ORG_Parallel

  if(sw_rsb==ON)then
    if(iteration_electronic==1) then
       push_solver=.false.
       call exec_submat_before_or_after()  ! contained here
       push_solver=.true.
    endif
    do ik=1,kv3
       call UnitaryTransform_WF(INVERSE,ik,zaj_l,.false.,.false.)
    enddo
    call m_ES_RSB_alloc()
    call m_ES_RSB_doit(verbose=.false.)
    call m_ES_EXX_gather_valence_states(nfout,.false.)
  endif


!!$ if(submat_before_renewal==ON.and.previous_isolver/=MATRIXDIAGON)then
  if(submat_before_renewal==ON)then
     if((sw_submat == ON .and. &
          & (isolver /= DAVIDSON .and. isolver /= SUBMAT .and.  &
          &  isolver /= MATRIXDIAGON .and. isolver/=MDDAVIDSON .and. isolver /= MDKOSUGI)) &
          & .or.renew_wf_again_m_CtrlP ) then
!!$        if(ipri>=1)write(nfout,'(" !solver doing subspace rotation before WF renewal")')
        call exec_submat_before_or_after()  ! contained here
     end if
  end if

! === KT_add === 13.0U3
  if ( precon_hardpart ) then
     if ( sw_wf_mixing == ON .and. wf_mixing_is_active ) then
        call m_ES_WF_cp_zaj_to_zaj_prev
     endif
  endif
! ============== 13.0U3

  if(hybfunc) sw_hybrid_functional=ON
  solver: select case(isolver)

     case (MATRIXDIAGON)
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,15)
        if(sw_hybrid_functional==ON.and.iteration==1) then
           hybfunc=.true.
           sw_hybrid_functional=OFF
        endif
        if(.not.(ekmode==ON .and. nk_in_the_process>1)) then
          call exec_matdiagon(isolver,iteration,iteration_electronic)
          call cp_chgq_to_chgqo()
        endif
                                                  __TIMER_STOP(15)
     case (SD,MSD)
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,12)
      call exec_sd(isolver,iteration,iteration_electronic,precon,dtim)
      call cp_chgq_to_chgqo()
                                                  __TIMER_STOP(12)
     case (lmSD, lmMSD, lmCG, lmeazyCG)
        call exec_lmm(isolver,previous_isolver,precon,dtim)

     case (RMM, RMM2, RMM2P)
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,18)
        call exec_rmm(isolver,previous_isolver,precon,dtim)
        call cp_chgq_to_chgqo()
                                                  __TIMER_STOP(18)

     case (SUBMAT)
        call exec_submat(isolver)
        call cp_chgq_to_chgqo()

     case (DAVIDSON)
        call exec_davidson(isolver,precon)
        call cp_chgq_to_chgqo()

     case (MDDAVIDSON)
        call m_CtrlP_set_sw_MRCV_only(ON)      ! sw_MRCV_only = ON
        call exec_mddavidson(isolver,precon)
        call cp_chgq_to_chgqo()

     case (MDKOSUGI)
        call m_CtrlP_set_sw_MRCV_only(OFF)     ! sw_MRCV_only = OFF
        call exec_mddavidson(isolver,precon)
        call cp_chgq_to_chgqo()

     case default
        if(ipri >= 1) write(6,'(" !! isolver = ",i5)') isolver
        call error
        !stop ' error at (Renewal_of_WaveFunctions)'
        call phase_error_with_msg(nfout, 'error at (Renewal_of_WaveFunctions)',__LINE__,__FILE__)
  end select solver

  if(submat_before_renewal/=ON) then
     if((sw_submat == ON .and.  &
          & (isolver /= DAVIDSON .and. isolver /= SUBMAT .and.  &
          &  isolver /= MATRIXDIAGON .and. isolver/=MDDAVIDSON .and. isolver /= MDKOSUGI)) &
          & .or.renew_wf_again_m_CtrlP) then
        call exec_submat_before_or_after()  ! contained here
     end if
  end if

  previous_isolver = isolver

  call m_ESIO_wd_WFs_standardout(nfout,ipriwf)

contains
  subroutine exec_submat_before_or_after()
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,14)
    if(sw_rsb==ON) disable_utransform = .true.
    if(push_solver) call m_CtrlP_push_SolverNameApplied("SUBMAT",6)
    call m_ESsubmat_alloc()

#ifndef DISABLE_NONCL
! ====================================== modified by K. Tagami =============== 11.0
    if ( noncol ) then
       if ( isolver /= MDKOSUGI .and. isolver /= MDDAVIDSON ) then
          call m_ESsubmat_renew_WF_noncl(nfout,meg,damp)
       endif
    else
#endif
       call m_ESsubmat_renew_WF(nfout,meg,damp)
#ifndef DISABLE_NONCL
    endif
! ============================================================================ 11.0
#endif

    if(sw_rsb==ON) disable_utransform = .false.
                                                  __TIMER_STOP(14)
  end subroutine exec_submat_before_or_after

  subroutine cp_chgq_to_chgqo()
    if(icond == INITIAL .or. icond == CONTINUATION) then
       call m_CD_cp_chgq_to_chgqo()
       if(sw_hubbard == ON) call m_OP_cp_ommix_to_omold()
    end if
  end subroutine cp_chgq_to_chgqo

  subroutine check_new_solver
    ! If "isolver" is not "previous_isolver", the value of "isolver" is printed.
    ! And if "isolver" is RMM?*, "iteration_rmm_start" is set in the module m_IterationNumbers
    !
    if(isolver /= previous_isolver) then
       if(ipri>=2) write(nfout,*) ' isolver = ', isolver
       if(isolver == RMM .or. isolver == RMM2 .or. isolver == RMM2P) then
          call m_Iter_set_rmm_start     ! iteration_rmm_start = iteration_electronic
          if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) &
               &  call m_ESrmm_reset_ng_maxmin()
          if(ipri >= 2) write(nfout,'(" !! iteration_rmm_start is set")')
          call m_ESrmm_reset_r_norm_flag ! ->rr_is_over_or_under(= UNDER)
       end if
    end if
  end subroutine check_new_solver

end subroutine Renewal_of_WaveFunctions

function EigenValue(ilevel)
  use m_Const_Parameters, only     : DP, Hartree
  use m_Electronic_Structure, only : m_ES_get_energy
  real(kind=DP) ::       EigenValue
  integer, intent(in) :: ilevel

  EigenValue = m_ES_get_energy(ilevel)*Hartree  ! energy unit is eV
end function EigenValue

! === EXEC_MATDIAGON ==
subroutine exec_matdiagon(isolver,iteration,iteration_electronic)
  use m_Files,                only : nfout
#ifndef DISABLE_NONCL
  use m_Control_Parameters,   only : noncol
  use m_ES_WF_by_MatDiagon,   only : m_ESmat_solve_Hx_eq_eSx_noncl
#endif
  use m_ES_WF_by_MatDiagon,   only : m_ESmat_solve_Hx_eq_eSx
  implicit none
  integer, intent(in) :: isolver,iteration,iteration_electronic

  call push_solver_name_applied(isolver)

#ifndef DISABLE_NONCL
! ======================================= modified by K. Tagami ============== 11.0
  if ( noncol ) then
     call m_ESmat_solve_Hx_eq_eSx_noncl(nfout,iteration,iteration_electronic)
  else
#endif
     call m_ESmat_solve_Hx_eq_eSx(nfout,iteration,iteration_electronic)
#ifndef DISABLE_NONCL
  endif
#endif
! ============================================================================= 11.0
end subroutine exec_matdiagon

! === EXEC_SD ===
subroutine exec_sd(isolver,iteration,iteration_electronic,precon,dtim)
  use m_Const_Parameters,     only : DP,ON
  use m_Files,                only : nfout
#ifndef DISABLE_NONCL
  use m_Control_Parameters,   only : noncol,sw_rsb
  use m_ES_WF_by_SDorCG,      only : m_ESsd_renew_WF_by_SDorCG_noncl
  use m_Electronic_Structure, only : m_ES_sort_eigen_vals_noncl
#endif
  use m_Electronic_Structure, only : m_ES_sort_eigen_values
  use m_ES_WF_by_SDorCG,      only : m_ESsd_renew_WF_by_SDorCG
  implicit none
  integer, intent(in) :: isolver,precon,iteration,iteration_electronic
  real(kind=DP), intent(in)  :: dtim

  call push_solver_name_applied(isolver)
#ifndef DISABLE_NONCL
! ========================================= modified by K. Tagami ============== 11.0
!         call m_ESsd_renew_WF_by_SDorCG(nfout,isolver,precon,dtim)
!
  if ( noncol ) then
     call m_ESsd_renew_WF_by_SDorCG_noncl(nfout,isolver,precon,dtim)
  else
#endif
     call m_ESsd_renew_WF_by_SDorCG(nfout,isolver,precon,dtim)
#ifndef DISABLE_NONCL
  endif
! =============================================================================== 11.0
#endif

! ------------------- Revised by T. Yamasaki,  21 July 2009 --->>
#ifndef DISABLE_NONCL
! ================================= modified by K. Tagami ============= 11.0
  if ( noncol ) then
     if(sw_rsb/=ON) call m_ES_sort_eigen_vals_noncl()
  else
#endif
     if(sw_rsb/=ON) call m_ES_sort_eigen_values() ! -> neordr, nrvf_ordr
#ifndef DISABLE_NONCL
  endif
! ===================================================================== 11.0
#endif
! <----
end subroutine exec_sd

! === EXEC_LMM ===
subroutine exec_lmm(isolver,previous_isolver,precon,dtim)
  use m_Const_Parameters,     only : DP, DAVIDSON, ON
  use m_Files,                only : nfout
!!$  use m_Control_Parameters,   only : m_CtrlP_dtim_1dsearch_now
  use m_ES_WF_by_SDorCG,      only : m_ESsd_renew_WF_by_SDorCG, m_ESsd_evolve_WFs_again
#ifndef DISABLE_NONCL
  use m_Control_Parameters,   only : noncol, sw_rsb
  use m_ES_WF_by_SDorCG,      only : m_ESsd_renew_WF_by_SDorCG_noncl
  use m_Electronic_Structure, only : m_ES_sort_eigen_vals_noncl
  use m_ES_ortho,             only : m_ES_modified_gramschmidt_noncl
#endif
  use m_Electronic_Structure, only : m_ES_sort_eigen_values
#ifdef MEMORY_SAVE_ZAJ_OLD
  use m_Control_Parameters, only : RMM2P_is_specified
  use m_ES_WF_by_SDorCG,    only : m_ESsd_alloc_zaj_old, m_ESsd_dealloc_zaj_old
#endif
  implicit none
  integer, intent(in) :: isolver,previous_isolver,precon
  real(kind=DP), intent(in)  :: dtim

  call push_solver_name_applied(isolver)
!!$         if(icond == INITIAL .or. icond == CONTINUATION) then
#ifdef MEMORY_SAVE_ZAJ_OLD
  if(.not. RMM2P_is_specified) call m_ESsd_alloc_zaj_old()
#endif

#ifndef DISABLE_NONCL
! ========================================= modified by K. Tagami ============== 11.0
!        call Renew_WF_by_lmSDorlmCG(isolver,precon,dtim)       ! -(contained here)
!
  if ( noncol ) then
     if ( previous_isolver == DAVIDSON ) then
        call m_ES_modified_gramschmidt_noncl(nfout)
     endif
     call Renew_WF_by_lmSDorlmCG_noncl( isolver,precon,dtim )
  else
#endif
     call Renew_WF_by_lmSDorlmCG(isolver,precon,dtim)       ! -(contained here)
#ifndef DISABLE_NONCL
  endif
! ============================================================================== 11.0
#endif
#ifdef MEMORY_SAVE_ZAJ_OLD
  if(.not. RMM2P_is_specified) call m_ESsd_dealloc_zaj_old()
#endif

#ifndef DISABLE_NONCL
! ================================= modified by K. Tagami ============= 11.0
  if ( noncol ) then
     if(sw_rsb/=ON) call m_ES_sort_eigen_vals_noncl()
  else
#endif
     if(sw_rsb/=ON) call m_ES_sort_eigen_values() ! -> neordr, nrvf_ordr
#ifndef DISABLE_NONCL
  endif
! ===================================================================== 11.0
#endif
! <----
                                                  __TIMER_STOP(13)

end subroutine exec_lmm

#ifdef CG_PREVIOUS
subroutine Renew_WF_by_lmSDorlmCG(isolver,precon,dtim)
  use m_Const_Parameters,     only : DP, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, BAND_ENERGY,CG &
       &                           , ORTHONORMALIZATION, CONTINUATION &
       &                           , TOTAL_ENERGY, MODIFIED_TOTAL_ENERGY, ON, INITIAL, OFF
  use m_Control_Parameters,   only : energy_evaluation, icond, sw_use_wfred, ipri &
       &                           , m_CtrlP_decide_dtim_1Dsearch, m_CtrlP_dtim_1Dsearch_now &
       &                           , sw_hybrid_functional,force_exx_energy1,sw_retard_eigval_evaluation &
       &                           , potential_update, m_CtrlP_set_in_line_minimization
  use m_Files,                only : nfout
  use m_Kpoints,              only : kv3
  use m_ES_WF_by_SDorCG,      only : m_ESsd_evolve_WFs_again &
       &                           , m_ESsd_alloc_wfred, m_ESsd_dealloc_wfred &
       &                           , m_ESsd_copy_zaj_to_zaj_old, m_ESsd_decide_CG_direction &
       &                           , m_ESsd_renew_WF_by_SDorCG
  use m_Charge_Density,       only : m_CD_cp_chgq_to_chgqo
  use m_Total_Energy,         only : m_TE_tell_total_energy &
       &                           , m_TE_tell_total_energy0 &
       &                           , m_TE_tell_band_energy &
       &                           , m_TE_tell_extended_band_energy
  use m_ES_ExactExchange,     only : m_ES_EXX_store_wfv, m_ES_EXX_restore_wfv
  implicit none
  integer, intent(in)         :: isolver,precon
  real(kind=DP),intent(in)    :: dtim

  real(kind=DP), dimension(3) :: etotal
  real(kind=DP), parameter    :: factor = 2
  real(kind=DP)               :: dtim_new = 0.d0, dtim_msdv
  integer                     :: isolver_core, mode

  integer                     :: energy_evaluation_lmm
  integer                     :: what_is_the_core_solver
  integer                     :: ik

  energy_evaluation_lmm = energy_evaluation
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) energy_evaluation_lmm = BAND_ENERGY
! <--

  isolver_core = what_is_the_core_solver(isolver) ! -(in this file)

  mode = ORTHONORMALIZATION

! ---------------- Added by T. Yamasaki, 28 June 2008, 31 Oct. 2008 ---
  if(sw_use_wfred == ON) call m_ESsd_alloc_wfred()
! ---------------------------------<<
  call m_ESsd_copy_zaj_to_zaj_old()

  call m_CtrlP_set_in_line_minimization(ON)
!!$  in_line_minimization = .true.

  if(energy_evaluation_lmm == TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy()
  else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
          &                                + m_TE_tell_extended_band_energy(nfout,kv3)
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(1) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

! ==============================================================================
! === WARN!!! by tkato 2012/02/11                                            ===
! === Right after re-started, the value of "dtim_1Dsearch" is -1             ===
! === and "dtim_msdv" is over-written with "dtim" here.                      ===
! === "dtim_1Dsearch" should be stored into restart file???                  ===
! === -> Resloved on phaseUnif@97 by tkato 2012/02/18                        ===
! === Subroutine m_CtrlP_wd_dtim_previous and m_CtrlP_rd_dtim_previous       ===
! === have been added into m_Control_Parameters.F90.                         ===
! === They write/read "dtim_1Dsearch" to/from "dtim_previous" entry          ===
! === in continue.data.                                                      ===
! ==============================================================================
  dtim_msdv = m_CtrlP_dtim_1Dsearch_now(dtim)
! ==============================================================================
  if(ipri >= 2) write(nfout,'(" !! dtim_msdv = ",f8.4)') dtim_msdv

  if(isolver_core == CG) then
     call m_ESsd_decide_CG_direction(precon)  !  -> betacg
     !    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv,iupdate=0)
     !    ~~~~~~~~~~~~~~~~~~~~~~~~~
  else
     call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv,iupdate=0)
     !    ~~~~~~~~~~~~~~~~~~~~~~~~~
  end if

  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     call m_CD_cp_chgq_to_chgqo()
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,7)
     if(sw_hybrid_functional==ON) call m_ES_EXX_store_wfv()
     call ChargeDensity_Construction(1)
                                                  __TIMER_STOP(7)
     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy0() -  m_TE_tell_band_energy(nfout,kv3) &
             &                                +  m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     if(icond == INITIAL .or. icond == CONTINUATION) call m_CD_cp_chgq_to_chgqo()
     etotal(2) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,dtim_msdv,factor*dtim_msdv,iupdate=0)
  ! (msdv_grad)

  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     if(sw_hybrid_functional==ON) call m_ES_EXX_restore_wfv()
                                                  __TIMER_START_w_BARRIER(MPI_CommGroup,7)
     call ChargeDensity_Construction(1)
                                                  __TIMER_STOP(7)
     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
             &                                + m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(3) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  if(iteration_electronic==1)then
     etotal(1) = etotal(2);etotal(2) = etotal(3)
     call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,dtim_msdv,factor*dtim_msdv,iupdate=0)
     call ChargeDensity_Construction(1)
     if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        call ChargeDensity_Construction(1)
        if(energy_evaluation_lmm == TOTAL_ENERGY) then
           etotal(3) = m_TE_tell_total_energy()
        else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
           etotal(3) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
                &                                + m_TE_tell_extended_band_energy(nfout,kv3)
        end if
     else if(energy_evaluation_lmm == BAND_ENERGY) then
        etotal(3) = m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  endif

  dtim_new = m_CtrlP_decide_dtim_1Dsearch(nfout,etotal,dtim_msdv,factor)
  if(ipri >= 2) write(nfout,'(" etotal(1:3) = ",3f15.6," dtim_msdv, dtim_new = ",2f15.6)') etotal(1:3), dtim_msdv, dtim_new
  mode = ORTHONORMALIZATION
  if(sw_retard_eigval_evaluation==ON.and.potential_update>0) force_exx_energy1=.true.
! ------------------- Revised by T. Yamasaki,  03 July 2008, 31 Oct 2008 --->>
  if(sw_use_wfred == ON) then
     call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,dtim_msdv,dtim_new,iupdate=0)
  else
     call m_ESsd_evolve_WFs_again(nfout,isolver_core,mode,factor*dtim_msdv,dtim_new,iupdate=0)
  end if

  call m_CtrlP_set_in_line_minimization(OFF)
!!$  in_line_minimization = .false.

! ------------------------------------------------------------<<
! ---------------- Added by T. Yamasaki, 28 June 2008, 31 Oct. 2008 ---
  if(sw_use_wfred == ON) call m_ESsd_dealloc_wfred()
! ------------------------------------------------------<<
end subroutine Renew_WF_by_lmSDorlmCG
#else

subroutine Renew_WF_by_lmSDorlmCG(isolver,precon,dtim)
  use m_Const_Parameters,     only : DP, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, BAND_ENERGY &
       &                           , ORTHONORMALIZATION, INITIAL, CONTINUATION &
       &                           , TOTAL_ENERGY, MODIFIED_TOTAL_ENERGY, CG, MSD, ON, OFF

  use m_Control_Parameters,   only : m_CtrlP_decide_dtim_1Dsearch, m_CtrlP_dtim_1Dsearch_now &
       &                           , sw_hybrid_functional, ipri &
       &                           , energy_evaluation, icond, sw_use_wfred, sw_retard_eigval_evaluation &
       &                           , potential_update, force_exx_energy1, m_CtrlP_set_in_line_minimization
  use m_Files,                only : nfout
  use m_Kpoints,              only : kv3
  use m_ES_WF_by_SDorCG,      only : m_ESsd_evolve_WFs_again,  m_ESsd_decide_CG_direction &
       &                           , m_ESsd_alloc_wfred, m_ESsd_dealloc_wfred &
       &                           , m_ESsd_copy_zaj_to_zaj_old, m_ESsd_decide_CG_direction &
       &                           , m_ESsd_renew_WF_by_SDorCG
  use m_Charge_Density,       only : m_CD_cp_chgq_to_chgqo
  use m_Total_Energy,         only : m_TE_tell_total_energy &
       &                           , m_TE_tell_total_energy0, m_TE_tell_band_energy &
       &                           , m_TE_tell_extended_band_energy
  use m_ES_ExactExchange,     only : m_ES_EXX_store_wfv, m_ES_EXX_restore_wfv
  use m_IterationNumbers,     only : iteration_electronic
  implicit none
  integer, intent(in)         :: isolver,precon
  real(kind=DP),intent(in)    :: dtim

  real(kind=DP), dimension(3) :: etotal
  real(kind=DP), parameter    :: factor = 2
  real(kind=DP)               :: dtim_new = 0.d0, dtim_msdv
  integer                     :: isolver_core, mode, isol
  integer                     :: energy_evaluation_lmm
  integer                     :: what_is_the_core_solver
  integer                     :: iupdate
  integer                     :: ispin,ik

  energy_evaluation_lmm = energy_evaluation
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) energy_evaluation_lmm = BAND_ENERGY

  isolver_core = what_is_the_core_solver(isolver) ! -(in this file)
  isol = isolver_core
  if(isolver_core==CG) isol=MSD

  mode = ORTHONORMALIZATION
  dtim_msdv = m_CtrlP_dtim_1Dsearch_now(dtim)

  if(sw_use_wfred == ON) call m_ESsd_alloc_wfred()
  call m_ESsd_copy_zaj_to_zaj_old()

  call m_CtrlP_set_in_line_minimization(ON)
!!$  in_line_minimization = .true.

  if(energy_evaluation_lmm == TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy()
  else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
          &                                + m_TE_tell_extended_band_energy(nfout,kv3)
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(1) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  if(ipri >= 2) write(nfout,'(" !! dtim_msdv = ",f8.4)') dtim_msdv
  call m_ESsd_renew_WF_by_SDorCG(nfout,isol,precon,dtim_msdv,iupdate=0)

  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     call m_CD_cp_chgq_to_chgqo()
     if(sw_hybrid_functional==ON) call m_ES_EXX_store_wfv()
     call ChargeDensity_Construction(1)
     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy0() -  m_TE_tell_band_energy(nfout,kv3) &
             &                                +  m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     if(icond == INITIAL .or. icond == CONTINUATION) call m_CD_cp_chgq_to_chgqo()
     etotal(2) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  call m_ESsd_evolve_WFs_again(nfout,isol,mode,dtim_msdv,factor*dtim_msdv,iupdate=0)
	 ! (msdv_grad)
  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     if(sw_hybrid_functional==ON.and.iteration_electronic/=1) call m_ES_EXX_restore_wfv()
     call ChargeDensity_Construction(1)
     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
             &                                + m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(3) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  if(iteration_electronic==1)then
     etotal(1) = etotal(2);etotal(2) = etotal(3)
     call m_ESsd_evolve_WFs_again(nfout,isol,mode,dtim_msdv,factor*dtim_msdv,iupdate=0)

! === ASMS DEBUG === 2015/04/13
!     call ChargeDensity_Construction(1)
     if ( icond /= FIXED_CHARGE .and. icond /= FIXED_CHARGE_CONTINUATION ) then
        call ChargeDensity_Construction(1)
     endif
! ================== 2015/04/13

     if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        if(sw_hybrid_functional==ON) call m_ES_EXX_restore_wfv()
        call ChargeDensity_Construction(1)
        if(energy_evaluation_lmm == TOTAL_ENERGY) then
           etotal(3) = m_TE_tell_total_energy()
        else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
           etotal(3) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
                &                                + m_TE_tell_extended_band_energy(nfout,kv3)
        end if
     else if(energy_evaluation_lmm == BAND_ENERGY) then
        etotal(3) = m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  endif

  dtim_new = m_CtrlP_decide_dtim_1Dsearch(nfout,etotal,dtim_msdv,factor)
  if(ipri >= 2) write(nfout,'(" etotal(1:3) = ",3f15.6," dtim_msdv, dtim_new = ",2f15.6)') etotal(1:3), dtim_msdv, dtim_new

  mode = ORTHONORMALIZATION
  if(sw_retard_eigval_evaluation==ON.and.potential_update>0) force_exx_energy1=.true.
  if(sw_use_wfred == ON) then
     call m_ESsd_evolve_WFs_again(nfout,isol,mode,dtim_msdv,dtim_new,iupdate=0)
  else
     call m_ESsd_evolve_WFs_again(nfout,isol,mode,factor*dtim_msdv,dtim_new,iupdate=0)
  end if
  if(isolver_core == CG) then
     call m_ESsd_decide_CG_direction(precon)  !  -> betacg
     !    ~~~~~~~~~~~~~~~~~~~~~~~~~~
       !call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv)
     call m_ESsd_renew_WF_by_SDorCG(nfout,isolver_core,precon,1.0d0)
     !    ~~~~~~~~~~~~~~~~~~~~~~~~~
  endif
  if(sw_use_wfred == ON) call m_ESsd_dealloc_wfred()
  call m_CtrlP_set_in_line_minimization(OFF)
!!$  in_line_minimization = .false.
end subroutine Renew_WF_by_lmSDorlmCG
#endif

#ifndef DISABLE_NONCL
! ===================================== added by K. Tagami =============== 11.0
subroutine Renew_WF_by_lmSDorlmCG_noncl(isolver,precon,dtim)
  use m_Const_Parameters,     only : DP, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, BAND_ENERGY &
       &                           , ORTHONORMALIZATION &
       &                           , TOTAL_ENERGY, MODIFIED_TOTAL_ENERGY, ON, CG &
       &                           , INITIAL, CONTINUATION, OFF

  use m_Control_Parameters,   only : energy_evaluation, icond, sw_use_wfred, ipri&
       &                           , sw_hybrid_functional &
       &                           , m_CtrlP_decide_dtim_1Dsearch &
       &                           , m_CtrlP_dtim_1Dsearch_now, m_CtrlP_set_in_line_minimization
  use m_Files,                only : nfout
  use m_Kpoints,              only : kv3
  use m_ES_WF_by_SDorCG,      only : m_ESsd_evolve_WFs_again,  m_ESsd_decide_CG_direction &
       &                           , m_ESsd_alloc_wfred, m_ESsd_dealloc_wfred &
       &                           , m_ESsd_copy_zaj_to_zaj_old, m_ESsd_decide_CG_direction &
       &                           , m_ESsd_renew_WF_by_SDorCG
  use m_ES_WF_by_SDorCG,      only : m_ESsd_renew_WF_by_SDorCG_noncl &
       &                           , m_ESsd_evolve_WFs_again_noncl &
       &                           , m_ESsd_decide_CG_direc_noncl
  use m_Charge_Density,       only : m_CD_cp_chgq_to_chgqo
  use m_Total_Energy,         only : m_TE_tell_total_energy &
       &                           , m_TE_tell_total_energy0, m_TE_tell_band_energy &
       &                           , m_TE_tell_extended_band_energy
  use m_ES_ExactExchange,     only : m_ES_EXX_store_wfv, m_ES_EXX_restore_wfv
  implicit none
  integer, intent(in)         :: isolver,precon
  real(kind=DP),intent(in)    :: dtim

  real(kind=DP), dimension(3) :: etotal
  real(kind=DP), parameter    :: factor = 2
  real(kind=DP)               :: dtim_new = 0.d0, dtim_msdv
  integer                     :: isolver_core, mode

  integer                     :: energy_evaluation_lmm
  integer                     :: what_is_the_core_solver

  energy_evaluation_lmm = energy_evaluation
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) energy_evaluation_lmm = BAND_ENERGY

  isolver_core = what_is_the_core_solver(isolver) ! -(in this file)

  mode = ORTHONORMALIZATION

  if(sw_use_wfred == ON) call m_ESsd_alloc_wfred()
  call m_ESsd_copy_zaj_to_zaj_old()
  call m_CtrlP_set_in_line_minimization(ON)
!!$  in_line_minimization = .true.
  if(energy_evaluation_lmm == TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy()
  else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     etotal(1) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
          &                                + m_TE_tell_extended_band_energy(nfout,kv3)
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(1) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  dtim_msdv = m_CtrlP_dtim_1Dsearch_now(dtim)
  if(ipri >= 2) write(nfout,'(" !! dtim_msdv = ",f8.4)') dtim_msdv

  if(isolver_core == CG) then
     call m_ESsd_decide_CG_direc_noncl(precon)  !  -> betacg
     call m_ESsd_renew_WF_by_SDorCG_noncl(nfout,isolver_core,precon,dtim_msdv)
  else
     call m_ESsd_renew_WF_by_SDorCG_noncl(nfout,isolver_core,precon,dtim_msdv)
  end if

  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     call m_CD_cp_chgq_to_chgqo()
     if(sw_hybrid_functional==ON) call m_ES_EXX_store_wfv()
     call ChargeDensity_Construction(1)

     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(2) = m_TE_tell_total_energy0() -  m_TE_tell_band_energy(nfout,kv3) &
             &                                +  m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     if(icond == INITIAL .or. icond == CONTINUATION) call m_CD_cp_chgq_to_chgqo()
     etotal(2) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  call m_ESsd_evolve_WFs_again_noncl( nfout, isolver_core, mode, &
       &                               dtim_msdv, factor*dtim_msdv )
                                         ! (msdv_grad)

  if(energy_evaluation_lmm == TOTAL_ENERGY .or. energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
     if(sw_hybrid_functional==ON) call m_ES_EXX_restore_wfv()
     call ChargeDensity_Construction(1)
     if(energy_evaluation_lmm == TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy()
     else if(energy_evaluation_lmm == MODIFIED_TOTAL_ENERGY) then
        etotal(3) = m_TE_tell_total_energy0() - m_TE_tell_band_energy(nfout,kv3) &
             &    + m_TE_tell_extended_band_energy(nfout,kv3)
     end if
  else if(energy_evaluation_lmm == BAND_ENERGY) then
     etotal(3) = m_TE_tell_extended_band_energy(nfout,kv3)
  end if

  dtim_new = m_CtrlP_decide_dtim_1Dsearch(nfout,etotal,dtim_msdv,factor)

  if(ipri >= 2) write(nfout,'(" etotal(1:3) = ",3f15.6," dtim_msdv, dtim_new = ",2f15.6)') etotal(1:3), dtim_msdv, dtim_new

  mode = ORTHONORMALIZATION

  if(sw_use_wfred == ON) then
     call m_ESsd_evolve_WFs_again_noncl( nfout, isolver_core, mode, &
          &                                  dtim_msdv, dtim_new )
  else
     call m_ESsd_evolve_WFs_again_noncl( nfout, isolver_core, mode, &
          &                                  factor*dtim_msdv, dtim_new )
  end if

  if(sw_use_wfred == ON) call m_ESsd_dealloc_wfred()
  call m_CtrlP_set_in_line_minimization(OFF)
!!$  in_line_minimization = .false.

end subroutine Renew_WF_by_lmSDorlmCG_noncl
! ===================================================================== 11.0
#endif

integer function what_is_the_core_solver(isolver)
  use m_Const_Parameters,     only : SD,MSD,CG &
       &                          , lmSD, lmMSD, lmCG, lmeazyCG, eazyCG
  implicit none
  integer, intent(in) :: isolver

  if(isolver == lmSD) then
     what_is_the_core_solver = SD
  else if(isolver == lmMSD) then
     what_is_the_core_solver = MSD
  else if(isolver == lmCG) then
     what_is_the_core_solver = CG
  else if(isolver == lmeazyCG) then
     what_is_the_core_solver = eazyCG
  else
     what_is_the_core_solver = MSD
  end if
end function what_is_the_core_solver

! === EXEC_RMM ===
subroutine exec_rmm(isolver,previous_isolver,precon,dtim)
  use m_Const_Parameters,     only : DP, DAVIDSON
  use m_Files,                only : nfout
#ifndef DISABLE_NONCL
  use m_Control_Parameters,   only : noncol
  use m_ES_ortho,             only : m_ES_modified_gramschmidt_noncl
  use m_ES_WF_by_RMM,          only : m_ESrmm_renew_WF_noncl
#endif
  use m_ES_WF_by_RMM,         only : m_ESrmm_renew_WF

#ifdef MEMORY_SAVE_ZAJ_OLD
  use m_Control_Parameters, only : RMM2P_is_specified
#endif
  implicit none
  integer, intent(in) :: isolver,previous_isolver,precon
  real(kind=DP), intent(in)  :: dtim

  call push_solver_name_applied(isolver)
#ifndef DISABLE_NONCL
! ================================= modified by K. Tagami ============= 11.0
!         call m_ESrmm_renew_WF(nfout,isolver,precon,dtim)
!
  if ( noncol ) then
     if ( previous_isolver == DAVIDSON ) then
        call m_ES_modified_gramschmidt_noncl(nfout)
     endif
     call m_ESrmm_renew_WF_noncl(nfout,isolver,precon,dtim)
  else
#endif
     call m_ESrmm_renew_WF(nfout,isolver,precon,dtim)

#ifndef DISABLE_NONCL
  endif
! ====================================================================== 11.0
#endif
end subroutine exec_rmm

! === EXEC_SUBMAT ===
subroutine exec_submat(isolver)
  use m_Files,                only : nfout
#ifndef DISABLE_NONCL
  use m_Control_Parameters,    only : noncol
  use m_ES_WF_by_submat,       only : m_ESsubmat_Renew_WF_noncl
#endif
  use m_Control_Parameters,    only : meg, damp
  use m_ES_WF_by_submat,       only : m_ESsubmat_Renew_WF, m_ESsubmat_alloc
  implicit none
  integer, intent(in) :: isolver

  call push_solver_name_applied(isolver)
  call m_ESsubmat_alloc()
#ifndef DISABLE_NONCL
! ============================================ modified by K. Tagami ======== 11.0
  if ( noncol ) then
     call m_ESsubmat_renew_WF_noncl(nfout,meg,damp)
  else
#endif
     call m_ESsubmat_renew_WF(nfout,meg,damp)
#ifndef DISABLE_NONCL
  endif
! =========================================================================== 11.0
#endif
end subroutine exec_submat

! === EXEC_DAVIDSON ===
subroutine exec_davidson(isolver,precon)
  use m_Files,                only : nfout
  use m_IterationNumbers,     only : iteration_electronic, iteration_ionic
#ifndef DISABLE_NONCL
  use m_Control_Parameters,    only : noncol
  use m_ES_WF_by_Davidson,     only : m_ESdavidson_Renew_WF_noncl
  use m_ES_ortho,              only : m_ES_modified_gramschmidt_noncl
#endif
  use m_ES_WF_by_Davidson,    only : m_ESdavidson_Renew_WF
  use m_ES_ortho,             only : m_ES_modified_gram_schmidt
  implicit none
  integer, intent(in) :: isolver, precon

  call push_solver_name_applied(isolver)
#ifndef DISABLE_NONCL
! ============================== modified by K. Tagami ============== 11.0
!    if(iteration_electronic==1 .and. iteration_ionic>1) call m_ES_modified_gram_schmidt(nfout)
!     call m_ESdavidson_Renew_WF(nfout,precon)
  if ( noncol ) then
     if ( iteration_electronic==1 .and. iteration_ionic>1 ) then
        call m_ES_modified_gramschmidt_noncl(nfout)
     endif
     call m_ESdavidson_Renew_WF_noncl(nfout,precon)
  else
#endif

     if ( iteration_electronic==1 .and. iteration_ionic>1 ) then
        call m_ES_modified_gram_schmidt(nfout)
     endif
     call m_ESdavidson_Renew_WF(nfout,precon)
#ifndef DISABLE_NONCL
  endif
! ==================================================================== 11.0
#endif
end subroutine exec_davidson

subroutine exec_mddavidson(isolver,precon)
  use m_Const_Parameters,     only : ON, OFF
  use m_Control_Parameters,   only : m_CtrlP_push_SolverNameApplied &
       &                           , submat_before_renewal,meg,damp, submat_GE
  use m_Files,                only : nfout
  use m_IterationNumbers,     only : iteration_electronic, iteration_ionic
#ifndef DISABLE_NONCL
  use m_Control_Parameters,    only : noncol
  use m_ES_WF_by_ModifiedDavidson,only : m_ESmddavid_Subspace_Rot_noncl &
       &                            ,    m_ESmddavid_Renew_WF_noncl
  use m_ES_ortho,             only : m_ES_modified_gramschmidt_noncl
  use m_ES_WF_by_submat,      only : m_ESsubmat_Renew_WF_noncl
#endif
  use m_ES_WF_by_ModifiedDavidson,only : m_ESmddavid_Subspace_Rotation &
       &                            ,    m_ESmddavid_Renew_WF
  use m_ES_ortho,             only : m_ES_modified_gram_schmidt
  use m_ES_WF_by_submat,      only : m_ESsubmat_Renew_WF, m_ESsubmat_alloc
  implicit none
  integer, intent(in) :: isolver,precon

#ifndef DISABLE_NONCL
! ========================= modified by K. Tagami ================ 11.0
!     if(submat_before_renewal==ON) call m_ESmddavid_Subspace_Rotation(nfout)
!     if(iteration_electronic==1 .and. iteration_ionic>1) call m_ES_modified_gram_schmidt(nfout)
!     call m_ESmdkosugi_Renew_WF(nfout,precon)
!     if(submat_before_renewal/=ON) call m_ESmddavid_Subspace_Rotation(nfout)

  if ( noncol ) then
     if (submat_before_renewal==ON) then
        if ( submat_GE == OFF ) then
           call m_ES_modified_gramschmidt_noncl(nfout)
           call m_CtrlP_push_SolverNameApplied("SUBMAT",6)
           call m_ESsubmat_alloc()
           call m_ESsubmat_renew_WF_noncl(nfout,meg,damp)
        else
           call m_ESmddavid_Subspace_Rot_noncl(nfout)
        endif
     endif
     if(submat_GE /= OFF) then
        if (iteration_electronic==1 .and. iteration_ionic>1) then
           call m_ES_modified_gramschmidt_noncl(nfout)
        endif
     endif

     call push_solver_name_applied(isolver)
     call m_ESmddavid_Renew_WF_noncl(nfout,precon)

     if (submat_before_renewal/=ON) then
        if(submat_GE == OFF) then
           call m_ES_modified_gramschmidt_noncl(nfout)
           call m_CtrlP_push_SolverNameApplied("SUBMAT",6)
           call m_ESsubmat_alloc()
           call m_ESsubmat_renew_WF_noncl(nfout,meg,damp)
        else
           call m_ESmddavid_Subspace_Rot_noncl(nfout)
        end if
     endif
  else
#endif
     if(submat_before_renewal == ON) then
        if(submat_GE == OFF) then
           call m_ES_modified_gram_schmidt(nfout)
           call m_CtrlP_push_SolverNameApplied("SUBMAT",6)
           call m_ESsubmat_alloc()
           call m_ESsubmat_renew_WF(nfout,meg,damp)
        else
           call m_ESmddavid_Subspace_Rotation(nfout)
        end if
     end if
     if(submat_GE /= OFF) then
! === DEBUG by tkato 2012/11/05 ================================================
        if((iteration_electronic == 1 .and. iteration_ionic > 1)) then
           call m_ES_modified_gram_schmidt(nfout)
        end if
! ==============================================================================
     end if

     call push_solver_name_applied(isolver)
     call m_ESmddavid_Renew_WF(nfout,precon)
     if(submat_before_renewal /= ON) then
        if(submat_GE == OFF) then
           call m_ES_modified_gram_schmidt(nfout)
           call m_CtrlP_push_SolverNameApplied("SUBMAT",6)
           call m_ESsubmat_alloc()
           call m_ESsubmat_renew_WF(nfout,meg,damp)
        else
           call m_ESmddavid_Subspace_Rotation(nfout)
        end if
     end if
#ifndef DISABLE_NONCL
  endif
#endif
! ==================================================================== 11.0
end subroutine exec_mddavidson

subroutine push_solver_name_applied(isolver)
  use m_Const_Parameters,     only : MATRIXDIAGON, SD, MSD, lmSD, lmMSD, lmCG, lmeazyCG, RMM, RMM2 &
       &                           , RMM2P, SUBMAT, DAVIDSON, MDDAVIDSON, MDKOSUGI
  use m_Control_Parameters,   only : m_CtrlP_push_SolverNameApplied, len_solvername
  implicit none
  integer, intent(in) :: isolver
  character(len=len_solvername) :: solvername
  character*20 :: name
  integer :: len_char
  name = ""
  solver: select case(isolver)
  case (MATRIXDIAGON)
     name = "MATDIAGON"
  case (SD)
     name = "SD"
  case (MSD)
     name = "MSD"
  case (lmSD)
     name = "lmSD"
  case (lmMSD)
     name = "lmMSD"
  case (lmCG)
     name = "CG"
  case (lmeazyCG)
     name = "lmeazyCG"
  case (RMM)
     name = "RMM3"
  case (RMM2)
     name = "RMM2"
  case (RMM2P)
     name = "RMM2P"
  case (SUBMAT)
     name = "SUBMAT"
  case (DAVIDSON)
     name = "DAVIDSON"
  case (MDDAVIDSON)
!!$     name = "MDDAVIDSON"
     name = "PKOSUGI"
  case (MDKOSUGI)
!!$     name = "MDKOSUGI"
     name = "PDAVIDSON"
  case default
     name = "notDefined"
  end select solver
  len_char = len_trim(name)
  solvername = ""
  solvername(1:min(len_char,len_solvername)) = name(1:min(len_char,len_solvername))
  call m_CtrlP_push_SolverNameApplied(solvername,min(len_char,len_solvername))
end subroutine push_solver_name_applied

subroutine wd_wfn()
  use m_FFT, only : m_FFT_alloc_WF_work,m_FFT_dealloc_WF_work
  use m_Control_Parameters, only : neg,nspin
  use m_Kpoints, only : kv3
  use m_Files, only : m_Files_open_nfwfk,nfwfk,nfout
  use m_ES_IO, only : m_ESIO_wd_WFn
  implicit none
  integer :: ik,ib
  call m_FFT_alloc_WF_work
  do ik=1,kv3
     do ib=1,neg
        write(nfout,*) 'foobar'
        call flush(nfout)
        call m_Files_open_nfwfk(nspin,ik,ib)
        write(nfout,*) 'hoge'
        call flush(nfout)
        call m_ESIO_wd_WFn(nfout,nfwfk,ik,ib)
        write(nfout,*) 'bibi'
        call flush(nfout)
     end do
  end do
  call m_FFT_dealloc_WF_work
end subroutine wd_wfn

