!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 626 $)
!
!  FUNCTION:  Ending_Time(), ckiter(), ChargeDensity_is_Converged(),
!             TotalEnergy_is_Divergent(), Forces_are_Converged(),
!             Already_Converged(), EigenValues_are_Converged(),
!              AllKpoints_are_Calculated()
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, December/01/2003, May/16/2005,
!                                     June/04/2005, April/10/2007
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
logical function Ending_Time()
! $Id: Convergence_Check.F90 626 2020-06-03 02:18:38Z jkoga $
  use m_Const_Parameters,   only : INITIAL, CONTINUATION, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                         , MAX_SCF_ITERATION_REACHED, MAX_MDSTEPS_REACHED, FSTOP
  use m_Files,              only : nfstop,nfout,m_Files_open_nfstop,m_Files_close_nfstop
  use m_IterationNumbers,   only : iteration, iteration_ionic, iteration_electronic
  use m_Control_Parameters, only : icond, ifstop, istop, max_total_scf_iteration, max_mdstep &
       &                         , max_TS_iteration_is_given, max_mdstep_is_given, printable &
       &                         , m_CtrlP_ckcput, m_CtrlP_rd_istop, terminated_because
  use m_Parallelization,    only : mype,MPI_CommGroup
  use mpi

  implicit none
!  include 'mpif.h'                                      ! MPI

  logical  :: torf
  integer :: ierr

  torf = m_CtrlP_ckcput()
  Ending_Time = torf

  if(.not.torf) &
       & Ending_Time = ckiter()           ! -(contained here)
!  call mpi_barrier(mpi_comm_world, ierr)
  if( .not. Ending_Time) then
    call Modify_Input_Settings()
  end if
  call mpi_barrier(MPI_CommGroup, ierr)
contains
  logical function ckiter()
    integer :: imax, icaseflag
! --> T. Yamasaki, 14 July 2008
    integer, parameter :: FILE_F_INP = 1, FILE_F_INP_MDSTEP = 2
! <--
    character*8 :: file_cname
    character*9 :: c_iteration
    icaseflag = 0
    ckiter = .false.
    if(icond == INITIAL .or. icond == CONTINUATION) then
       imax = max(iteration, iteration_ionic)
#ifndef _EMPIRICAL_
       if(iteration == 0) imax = 0
#endif
!                                                                                 (max_TS_iteration_is_given, max_mdstep_is_given)
       if(max_TS_iteration_is_given .and. max_mdstep_is_given) then                   ! (.true., .true.)
! --> T. Yamasaki, 07 Aug. 2009
!!$          if(max_total_scf_iteration <= imax .or. max_mdstep <= iteration_ionic) then
          if(max_total_scf_iteration <= imax .or. max_mdstep <= iteration_ionic-1) then
! <--
             ckiter = .true.
             if(max_total_scf_iteration <= imax) then
                icaseflag = FILE_F_INP
                terminated_because = MAX_SCF_ITERATION_REACHED
             else
                icaseflag = FILE_F_INP_MDSTEP
                terminated_because = MAX_MDSTEPS_REACHED
             end if
          else if(ifstop >= 1) then
             ckiter = ckiter_ifstop(imax)
             if(ckiter) terminated_because = FSTOP
          end if
       else if(.not.max_mdstep_is_given) then                                         ! (.true., .false) or (.false., .false)
          if(max_total_scf_iteration <= imax) then
             ckiter = .true.
             icaseflag = FILE_F_INP
             terminated_because = MAX_SCF_ITERATION_REACHED
          else if(ifstop >= 1) then
             ckiter = ckiter_ifstop(imax)
             if(ckiter) terminated_because = FSTOP
          end if
       else if(.not.max_TS_iteration_is_given .and. max_mdstep_is_given) then         ! (.false., .true.)
! --> T. Yamasaki, 06 Aug. 2009
!!$          if(max_mdstep <= iteration_ionic) then
          if(max_mdstep <= iteration_ionic-1) then
! <--
             ckiter = .true.
             icaseflag = FILE_F_INP_MDSTEP
             terminated_because = MAX_MDSTEPS_REACHED
          else if(ifstop >= 1 ) then
             ckiter = ckiter_ifstop(imax)
             if(ckiter) terminated_because = FSTOP
          end if
       end if


!!$       if(max_total_scf_iteration <= imax .or. max_mdstep <= iteration_ionic) then
!!$          ckiter = .true.
!!$          if(max_total_scf_iteration <= imax) then
!!$             icaseflag = FILE_F_INP
!!$          else
!!$             icaseflag = FILE_F_INP_MDSTEP
!!$          end if
!!$       else if(ifstop >= 1) then
!!$          if(mod(imax,ifstop) == 0) then
!!$             call m_Files_open_nfstop
!!$             call m_CtrlP_rd_istop(nfstop)        ! ->istop
!!$             if(istop >= 0 .and. istop <= imax) ckiter = .true.
!!$             call m_Files_close_nfstop
!!$          else
!!$             if(istop >= 0 .and. istop <= imax) ckiter = .true.
!!$          end if
!!$       end if
    else if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then

       if(iteration > max_total_scf_iteration) then
          ckiter = .true.
          icaseflag = FILE_F_INP
          terminated_because = MAX_SCF_ITERATION_REACHED
       else
          imax = max(iteration, iteration_electronic)
!!$          imax = iteration_electronic
!!$       if(imax > ek_max_iteration) then
!!$          ckiter = .true.
!!$       else
!!$          write(nfout,'(" iteration_electronic = ",i6)') iteration_electronic
!!$          write(nfout,'("  ifstop = ",i6)') ifstop
          if(ifstop >= 1) then
             ckiter = ckiter_ifstop(imax)
             if(ckiter) terminated_because = FSTOP
          end if
       end if
    end if
    if(ckiter) then
       if(icaseflag .eq. FILE_F_INP_MDSTEP) then
          file_cname = ' <F_INP>'
          imax = max_mdstep
          c_iteration = ' mdstep  '
       else if(icaseflag .eq. FILE_F_INP) then
          file_cname = ' <F_INP>'
          imax = max_total_scf_iteration
          c_iteration = 'iteration'
       else
          file_cname = '<F_STOP>'
          imax = istop
          c_iteration = 'iteration'
       end if
       if(printable) write(nfout,'(" --- This job stops because the ",a9," has exceeded " &
                              & ,i6," set in the ",a8," file ---")') c_iteration, imax,file_cname
    end if
  end function ckiter

  logical function ckiter_ifstop(imax)
    integer, intent(in) :: imax
    ckiter_ifstop = .false.
    if(mod(imax,ifstop) == 0) then
       call m_Files_open_nfstop
       call m_CtrlP_rd_istop(nfstop)        ! ->istop
       if(istop >= 0 .and. istop <= imax) ckiter_ifstop = .true.
       call m_Files_close_nfstop
    else
       if(istop >= 0 .and. istop <= imax) ckiter_ifstop = .true.
    end if

  end function ckiter_ifstop

end function Ending_Time

#ifndef _EMPIRICAL_
logical function ChargeDensity_is_Converged()
  use m_Files,          only : nfout
  use m_Const_Parameters, only : ON, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, WHOLE, INITIALLY, DP, YES &
       &                   ,     CHARGE_CONVERGED
  use m_Control_Parameters, only : printable, sw_submat_is_on, submat_is_done_this_iter, icond &
       &                     , m_CtrlP_set_renew_wf, sw_hubbard, nspin, m_CtrlP_get_edelta &
       &                     , m_CtrlP_clear_lmm_status_store &
       &                     , max_scf_iteration, max_scf_iteration_is_given &
       &                     , m_CtrlP_set_iconvergence, terminated_because, istress
  use m_Charge_Density, only : m_CD_is_converged_core
  use m_Total_Energy,   only : m_TE_is_converged &
      &                      , m_TE_Converged_Hubbard_Energy, m_TE_edeltb
  use m_ES_WF_by_submat,only : m_ESsubmat_dealloc, m_ESsubmat_reset_non_diagon
  use m_Force,          only : m_Force_get_fec_count
  use m_ES_occup,       only : m_ESoc_efermi_diff, m_ESoc_new_total_spin, m_ESoc_total_spin0
  use m_Crystal_Structure,only : sw_fix_total_spin_in, sw_fix_total_spin,spin_fix_period &
       &                     , m_CS_set_total_spin, m_CS_free_fixed_total_spin, m_CS_fix_total_spin
  use m_Ionic_System,   only : natm2
  use m_IterationNumbers, only : iteration_electronic, iteration_ionic

! ====================== added by K. Tagami ======================= 12.0Exp
  use m_Const_Parameters,    only  : ALL_AT_ONCE, EK_CONVERGED
  use m_Control_Parameters,  only : fixed_charge_k_parallel, &
       &                            m_CtrlP_set_iconvergence,xctype
! ================================================================ 12.0Exp

! ====================== KT_Test ========================================= 12.5Exp
  use m_Const_Parameters,    only  : OFF
  use m_Control_Parameters,   only : sw_hybrid_functional, &
       &                             truncate_vxw_updating, sw_update_vxw, oneshot
! ======================================================================== 12.5Exp

! === Postitron SCF === 2015/11/28
  use m_Control_Parameters,  only : sw_positron, positron_method
  use m_Const_Parameters,   only : positron_GGGC
  use m_Positron_Wave_Functions,  only : m_pWF_update_lifetime
! ===================== 2015/11/28

  use m_Const_Parameters, only : CHG_CONVERGENCE_REACHED

  implicit none
  logical, save :: renew_wf_again = .false.
  logical :: EigenValues_are_Converged
  integer ::       fec_counter
!!$  real(kind=DP) :: total_spin_new, m_ESoc_new_total_spin, m_ESoc_efermi_diff &
!!$       &          ,Delta_Efermi
  real(kind=DP) :: total_spin_new, Delta_Efermi
  logical :: totalenergy_is_converged
  real(kind=DP) :: edeltb, edelta
  real(kind=DP),parameter :: delta_efermi_lower = 0.002
  real(kind=DP),parameter :: factor_initial = 10.0
  real(kind=DP),save :: factor = factor_initial

  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
     if(EigenValues_are_Converged()) then
        ChargeDensity_is_Converged = .true.
! ==================================== added by K. Tagami ============ 12.0Exp
        if ( fixed_charge_k_parallel == ALL_AT_ONCE ) then
           call m_CtrlP_set_iconvergence(EK_CONVERGED)
        endif
! ==================================================================== 12.0Exp
     else
        ChargeDensity_is_Converged = .false.
     end if
     return
  end if

  call m_Force_get_fec_count(fec_counter)
  if(fec_counter /= 0) then
     ChargeDensity_is_Converged = .true.
     if(printable) &
          & write(nfout,'(" !Convergence_Check: fec_counter = ",i8," is not zero")') fec_counter
  else
     if(renew_wf_again) then
        ChargeDensity_is_Converged = .true.
     else
! --> T. Yamasaki 07th Aug. 2009
        totalenergy_is_converged = m_TE_is_converged(nfout)
!!$        if(m_TE_is_converged(nfout)) then
        if(totalenergy_is_converged) then
           if(nspin == 2 .and. sw_fix_total_spin == YES .and. spin_fix_period == INITIALLY) then
              ChargeDensity_is_Converged = .false.
           else
              ChargeDensity_is_Converged = m_CD_is_converged_core(nfout)
           end if
! <--
        else
           ChargeDensity_is_Converged = .false.
        endif
! --> T. Yamasaki 06th Aug. 2009
        if(nspin == 2 .and. sw_fix_total_spin == YES .and. spin_fix_period >= INITIALLY) then
! -->       17th Aug. 2009
           if(spin_fix_period.gt.1 .and. iteration_electronic>=spin_fix_period)then
              call m_CS_free_fixed_total_spin(nfout)
           else if(iteration_electronic>=20 .and. spin_fix_period==INITIALLY) then
!!$                 Delta_efermi = m_ESoc_efermi_diff()
!!$              if(totalenergy_is_converged .and. dabs(Delta_efermi) < delta_efermi_lower ) then
!!$              if(totalenergy_is_converged .and. dabs(Delta_efermi) < delta_efermi_lower ) then
!!$                 if(printable) write(nfout,'(" !** dabs(Delta_efermi) = dabs(",f16.8 &
!!$                      &                        ,") < Delta_efermi_lower = ",f16.8, " and total energy is converged")') &
!!$                      &               delta_efermi, delta_efermi_lower
!!$                 call m_CS_free_fix_total_spin(nfout) ! sw_fix_total_spin = NO
              edeltb = m_TE_edeltb()
              call m_CtrlP_get_edelta(edelta)
              if(dabs(edeltb) < edelta*natm2*factor) then
                 factor = factor/2.0
                 if(factor < 1.0) factor = 1.0
                 if(printable) write(nfout,'(" dabs(edeltb) = dabs(",f16.12,") < edelta*natm2*factor = ",f16.12)') &
                      & edeltb, edelta*natm2*factor
! -->
                 if(iteration_ionic >= 10) then   ! This is a test case
! <--
                    Delta_efermi = m_ESoc_efermi_diff()
                    if( dabs(Delta_efermi) < delta_efermi_lower ) then
                       if(printable) write(nfout,'(" !** dabs(Delta_efermi) = dabs(",f16.8 &
                            &                        ,") < Delta_efermi_lower = ",f16.8, " and total energy is converged")') &
                            &               delta_efermi, delta_efermi_lower
                       call m_CS_free_fixed_total_spin(nfout) ! sw_fix_total_spin = NO
                    else
                       if(totalenergy_is_converged) ChargeDensity_is_Converged = m_CD_is_converged_core(nfout)
                       if(.not.ChargeDensity_is_Converged) then
                          total_spin_new = m_ESoc_new_total_spin()
                          call m_CS_set_total_spin(nfout,total_spin_new)
                       end if
                    end if
                 else
                    if(totalenergy_is_converged) ChargeDensity_is_Converged = m_CD_is_converged_core(nfout)
                    if(.not.ChargeDensity_is_Converged) then
                       total_spin_new = m_ESoc_new_total_spin()
                       call m_CS_set_total_spin(nfout,total_spin_new)
                    end if
                 end if
              end if
           end if
!!$           if(iteration_electronic>=30 .and. mod(iteration_electronic,10)==1) then
!!$              Delta_efermi = m_ESoc_efermi_diff()
!!$              if(totalenergy_is_converged .and. dabs(Delta_efermi) < delta_efermi_lower ) then
!!$                 if(printable) write(nfout,'(" !** dabs(Delta_efermi) = dabs(",f16.8 &
!!$                      &                        ,") < Delta_efermi_lower = ",f16.8, " and total energy is converged")') &
!!$                      &               delta_efermi, delta_efermi_lower
!!$                 call m_CS_free_fix_total_spin(nfout) ! sw_fix_total_spin = NO
!!$              else
!!$                 total_spin_new = m_ESoc_new_total_spin()
!!$                 call m_CS_set_total_spin(nfout,total_spin_new)
!!$              end if
!!$           end if
! <--
        end if
! <--
     end if

     if ( sw_positron /= OFF ) then
        if ( positron_method == positron_GGGC ) call m_pWF_update_lifetime
     endif

! =========================== KT_Test ================= 12.5Exp
     if ( sw_hybrid_functional == ON ) then
        if ( truncate_vxw_updating .and. sw_update_vxw == OFF ) then
           ChargeDensity_is_Converged = .true.
        end if
     end if
! ====================================================== 12.5Exp

     if(sw_hubbard == ON) then
        if(m_TE_Converged_Hubbard_Energy(nfout)) then
           ChargeDensity_is_Converged = .true.
        end if
     end if

! --> T. Yamasaki, 18th Aug. 2009
     if(max_scf_iteration_is_given) then
        if(iteration_electronic >= max_scf_iteration) then
           if( .not. ChargeDensity_is_Converged .and. istress==ON) then
             stop 'max_scf_iteration cannot be used in conjunction with stress calculations'
           endif
           ChargeDensity_is_Converged = .true.
        end if
     end if
! <--
     if(ChargeDensity_is_Converged) then
        call m_CtrlP_set_iconvergence(CHARGE_CONVERGED)
        if(renew_wf_again) then
           renew_wf_again = .false.
           call m_ESsubmat_dealloc()
        else if(sw_submat_is_on .and. .not.renew_wf_again) then
           if(.not.submat_is_done_this_iter) then
              if(printable) write(nfout,'(" !! sw_submat_is_on = .true. and renew_wf_again = .false.")')
              ChargeDensity_is_Converged = .false.
              renew_wf_again = .true.
              call m_ESsubmat_reset_non_diagon()
           end if
        else ! if(.not.sw_submat_is_on .and. .not.renew_wf_again) then
           call m_ESsubmat_dealloc()
        end if
! --> T. Yamasaki, 17th Aug. 2009
        if(nspin == 2 .and. sw_fix_total_spin_in == YES .and. spin_fix_period == INITIALLY) then
           call m_CS_fix_total_spin(nfout)
           total_spin_new = m_ESoc_total_spin0()
           call m_CS_set_total_spin(nfout,total_spin_new)
           factor = factor_initial
        end if
! <--
        call m_CtrlP_clear_lmm_status_store() ! --> T. Yamasaki, 18th Aug. 2009
#ifndef DISABLE_VDWDF
        if(xctype == 'vdwdf'.and.oneshot) call vdW_oneshot()
#endif
     end if
     call m_CtrlP_set_renew_wf(renew_wf_again)
  end if
  if(ChargeDensity_is_Converged) then
    terminated_because = CHG_CONVERGENCE_REACHED
  endif
end function ChargeDensity_is_Converged

logical function TotalEnergy_is_Divergent()
  use m_Const_Parameters, only : SCF_ITERATION
  use m_Total_Energy, only : m_TE_is_Divergent_core
  use m_Files,        only : nfout
! --> T. Yamasaki, 22nd July 2008
  use m_Control_Parameters,only : ipri,m_CtrlP_etot_1dsrch_divergent
  use m_IterationNumbers,  only : iteration_ionic, iteration_electronic
  implicit none
  call Checkpoint_File(SCF_ITERATION)
  TotalEnergy_is_Divergent = m_TE_is_Divergent_core(nfout)
  if(iteration_ionic > 1 .and. .not.TotalEnergy_is_Divergent) then
     TotalEnergy_is_Divergent = m_CtrlP_etot_1dsrch_divergent(iteration_electronic)
     if(TotalEnergy_is_Divergent) then
        if(ipri >= 1) write(nfout,'(" m_CtrlP_etot_1dsrch_divergent is .true.")')
     end if
  end if
! <--
end function TotalEnergy_is_Divergent

#endif

logical function AllForces_are_Converged()
  AllForces_are_Converged = .true.
end function AllForces_are_Converged

logical function Forces_are_Converged2()
  use m_Control_Parameters, only : sw_phonon, m_CtrlP_set_iconvergence, terminated_because
  use m_Const_Parameters, only : ON, OFF, FORCE_CONVERGED, FORCE_CONVERGENCE_REACHED &
  &                            , MAX_PHSTEPS_REACHED
  use m_Ionic_System, only : iend_phonon, istart_phonon
  use m_IterationNumbers, only : iteration_ionic
  use m_Crystal_Structure, only : m_CS_phonon_symmetry
  implicit none
  Forces_are_Converged2 = .false.
  if (sw_phonon == OFF) return
  Forces_are_Converged2 = iteration_ionic>iend_phonon-istart_phonon+1
  if(Forces_are_Converged2) then
    call m_CtrlP_set_iconvergence(FORCE_CONVERGED)
    !terminated_because = FORCE_CONVERGENCE_REACHED
    terminated_because = MAX_PHSTEPS_REACHED
  endif
  if(Forces_are_Converged2) call m_CS_phonon_symmetry(ON)
end function Forces_are_Converged2


logical function Forces_are_Converged()
  use m_Control_Parameters, only : imdalg, ipriforce, icond, sw_phonon_oneshot &
       &                         , m_CtrlP_set_iconvergence &
       &                         , m_CtrlP_what_is_mdalg,sw_optimize_lattice, sw_stress_correction &
       &                         , sw_optimize_coords_sametime &
       &                         , sw_displace_atom, terminated_because &
#ifndef _EMPIRICAL_
       &                         , m_CtrlP_renew_edelta_ontheway, m_CtrlP_edelta_final &
       &                         , m_CtrlP_reset_edelta_ontheway
#else
       &                         , m_CtrlP_renew_edelta_ontheway, m_CtrlP_edelta_final
#endif
  use m_Const_Parameters,   only : DP,QUENCHED_CONSTRAINT,BLUEMOON &
       &                         , TEMPERATURE_CONTROL, FORCE_CONVERGED &
       &                         , VERLET, CNSTRA, ORDINA, WITHOUTTAG &
       &                         , WITHTAG, PHONON_FORCE &
       &                         , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, ON, INITIAL, OFF &
       &                         , PT_CONTROL, P_CONTROL, STREVL_ITERATION, FORCE_CONVERGENCE_REACHED &
       &                         , MAX_PHSTEPS_REACHED, QUENCHED_MD
  use m_Files,              only : nfout
  use m_IterationNumbers,   only : iteration_ionic,iteration,iteration_electronic
  use m_Ionic_System,       only : constraints_exist, move_constrained_plane, apply_constant_force &
       &                         , ekina,ega,almda,mdmode &
       &                         , forcmx_constraint_quench, cps, natm &
       &                         , m_IS_force_check_md_cnstr &
       &                         , m_IS_set_mdmode_cnstra &
       &                         , m_IS_set_mdmode_ordina &
       &                         , m_IS_evaluate_v_verlet &
       &                         , m_IS_moved_distance_of_plane_is_over & ! T. Yamaksai, 29 Apr. 2021
       &                         , m_IS_iter_ionic_is_over &              ! T. Yamasaki,  9 Jun  2021
! --> T. Yamasaki, 18 July 2008
       &                         , m_IS_cnstr_is_fixed_nhp &
       &                         , m_IS_force_check_md_nhp &
       &                         , m_IS_rigid_body_converged &
       &                         , m_IS_rigid_body_exists
! <--
#ifndef _EMPIRICAL_
  use m_Phonon,             only : m_Phonon_Check_iteration
#endif
  use m_Force,              only : forcmx, forc_l &
       &                         , m_Forces_are_Converged_core
!!$       &                         , m_Forces_clear_fec_counter
  use m_Total_Energy,       only : etotal
  use m_Stress,             only : m_Stress_in_correction

  implicit none
  integer :: mdalg, iret
  logical :: unitcell_can_change

  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
     Forces_are_Converged = .true.
     return
  end if

!  if(m_Stress_in_correction(2)) then
  if(m_Stress_in_correction(2)) then
    Forces_are_Converged = .true.
    return
  endif

!  call Checkpoint_File(STREVL_ITERATION)
  if(sw_displace_atom == ON) then
     Forces_are_Converged = .true.
     !return
  endif

  if(imdalg == QUENCHED_CONSTRAINT) then
     if(constraints_exist .and. move_constrained_plane) then
        if(m_IS_force_check_md_cnstr() .and. m_CtrlP_edelta_final()) then
           call m_IS_set_mdmode_cnstra()
#ifndef _EMPIRICAL_
           call m_CtrlP_reset_edelta_ontheway()
#endif
        else
           call m_IS_set_mdmode_ordina()
        end if
        Forces_are_Converged = .false.
        if( m_IS_moved_distance_of_plane_is_over()) Forces_are_Converged = .true.
     else
        call m_IS_set_mdmode_ordina()
        Forces_are_Converged = m_IS_force_check_md_cnstr()
        if(.not.m_CtrlP_edelta_final()) Forces_are_Converged = .false.
        if(ipriforce >= 3)  then
           if(Forces_are_Converged) then
              write(nfout,'(" Forces_are_Converged(2) = .true.")')
           else
              write(nfout,'(" Forces_are_Converged(2) = .false.")')
           end if
        end if
     end if
  else if(imdalg == PT_CONTROL .or. imdalg == P_CONTROL) then
     Forces_are_Converged = .true.
  else if(imdalg == TEMPERATURE_CONTROL .or. imdalg == VERLET) then
     Forces_are_Converged = .false.
#ifndef _EMPIRICAL_
  else if(imdalg == PHONON_FORCE) then
     if(sw_phonon_oneshot == ON) then
        Forces_are_Converged = .true.
     else
        Forces_are_Converged = m_Phonon_Check_iteration(iteration_ionic)
	icond = INITIAL  ! asms
     end if
#endif
  else if(imdalg /= BLUEMOON) then
     if(constraints_exist .and. move_constrained_plane) then
        if(apply_constant_force <= 0) then
           if((m_Forces_are_Converged_core() .and. m_CtrlP_edelta_final()) &
                &                        .or. m_IS_iter_ionic_is_over()) then ! This is a different
           !                                                                  line to the case of
           !                                                                 (imdalg ==
           !                                                                   QUENCHED_CONSTRAINT).
           call m_IS_set_mdmode_cnstra()
#ifndef _EMPIRICAL_
           call m_CtrlP_reset_edelta_ontheway()
#endif
           else
              call m_IS_set_mdmode_ordina()
           end if
        else
           call m_IS_set_mdmode_ordina()
        end if
        Forces_are_Converged = .false.
        if( m_IS_moved_distance_of_plane_is_over()) Forces_are_Converged = .true. ! T. Yamasaki, 29 Apr. 2021
!  --> T. Yamasaki, 18 July 2008
     else if(constraints_exist .and. m_IS_cnstr_is_fixed_nhp()) then
        if(m_IS_force_check_md_nhp()) then
           Forces_are_Converged = .true.
        else
           Forces_are_Converged = .false.
        end if
! <--
     else
        call m_IS_set_mdmode_ordina()
        Forces_are_Converged = m_Forces_are_Converged_core()
        if(.not.m_CtrlP_edelta_final()) Forces_are_Converged = .false.
     end if
  end if
  if(Forces_are_Converged .and. m_IS_rigid_body_exists()) then
     Forces_are_Converged = m_IS_rigid_body_converged()
  endif
  if(Forces_are_Converged) call m_CtrlP_set_iconvergence(FORCE_CONVERGED)
  mdalg = m_CtrlP_what_is_mdalg()

  if(Forces_are_Converged) then
!!$     call m_Forces_clear_fec_counter()
     if(mdalg == VERLET) call m_IS_evaluate_v_verlet(mdalg,forc_l)

     if ( sw_optimize_coords_sametime == OFF ) then
        call wd_forces_cps_etotal_and_etc(mdalg,.true.) ! nfdynm
     endif

     if(.not.unitcell_can_change()) then
       call wd_forces_cps_etotal_and_etc(mdalg,.false.) ! nfenf
     else if (sw_optimize_lattice==OFF .and. sw_stress_correction==ON) then
       call wd_forces_cps_etotal_and_etc(mdalg,.false.) ! nfenf
     endif
#ifndef _EMPIRICAL_
     call m_CtrlP_reset_edelta_ontheway()
#endif
     if(imdalg == PHONON_FORCE)then
       terminated_because = MAX_PHSTEPS_REACHED
     else
       terminated_because = FORCE_CONVERGENCE_REACHED
     endif
  else
     if(imdalg /= BLUEMOON .and. imdalg /= TEMPERATURE_CONTROL &
          & .and. imdalg /= VERLET .and. imdalg /= PHONON_FORCE) then
        if(imdalg /= QUENCHED_CONSTRAINT) then
           call m_CtrlP_renew_edelta_ontheway(nfout,etotal,forcmx,iteration_ionic,iteration,iteration_electronic)
        else if(imdalg == QUENCHED_CONSTRAINT) then
           call m_CtrlP_renew_edelta_ontheway(nfout,etotal,forcmx_constraint_quench,iteration_ionic,iteration,iteration_electronic)
        end if
     end if
  end if
  if(.not. Forces_are_Converged) call Checkpoint_File(STREVL_ITERATION)
  write(nfout,'(" <<Forces_are_Converged>>, Forces_are_Converged = ",L4)') Forces_are_Converged
  call flush(nfout)

end function Forces_are_Converged

logical function UnitCell_Converged(force_converged)
   use m_Control_Parameters, only : sw_optimize_lattice,imdalg,m_CtrlP_reset_iconvergence,icond &
        &                           , terminated_because,sw_stress_correction &
        &                           , sw_optimize_coords_sametime
   use m_Const_Parameters, only : ON, OFF,FIXED_CHARGE,FIXED_CHARGE_CONTINUATION,UNITCELL_ITERATION &
      &                         , STRESS_CONVERGENCE_REACHED, PT_CONTROL, P_CONTROL, MAX_MDSTEPS_REACHED
   use m_UnitCell, only : m_UC_converged
   use m_Stress,   only : m_Stress_in_correction
   implicit none
   logical, intent(in) :: force_converged
   logical :: unitcell_can_change,Rightafter_stress_correction
   if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
      UnitCell_Converged = .true.
      return
   endif
   if(m_Stress_in_correction()) then
      UnitCell_Converged = .false.
      return
   endif
   if(.not.unitcell_can_change())then
      UnitCell_Converged = force_converged
      return
   endif
   if(sw_optimize_lattice==OFF .and. sw_stress_correction==ON .and. .not.m_Stress_in_correction()) then
      UnitCell_Converged = force_converged
      return
   endif
   UnitCell_Converged = m_UC_converged() .and. force_converged
   if ( sw_optimize_coords_sametime == ON ) then
      call wd_forces_cps_etotal_and_etc(imdalg,.false.) ! nfenf
      call wd_forces_cps_etotal_and_etc(imdalg,.true.) ! nfdynm
   else
      if(force_converged) call wd_forces_cps_etotal_and_etc(imdalg,.false.) ! nfenf
   endif

   if(.not.UnitCell_Converged) then
      call m_CtrlP_reset_iconvergence()
      call Checkpoint_File(UNITCELL_ITERATION)
   else
      if(imdalg == PT_CONTROL .or. imdalg == P_CONTROL) then
        terminated_because = MAX_MDSTEPS_REACHED
      else
        terminated_because = STRESS_CONVERGENCE_REACHED
      endif
   endif
end function UnitCell_Converged

subroutine wd_forces_cps_etotal_and_etc(mdalg,flag_dynm)
  use m_Const_Parameters,   only : DP,QUENCHED_CONSTRAINT, BLUEMOON &
       &                         , TEMPERATURE_CONTROL &
       &                         , VERLET, CNSTRA, WITHOUTTAG &
       &                         , WITHTAG, OFF, ON, P_CONTROL, PT_CONTROL, DRIVER_URAMP, DRIVER_SC_DFT
  use m_Files,              only : nfout, nfenf, nfdynm, m_Files_flush_nfdynm, m_Files_flush_nfenf, nfdynm_cif
  use m_IterationNumbers,   only : iteration_ionic,iteration,iteration_electronic,iteration_unit_cell &
       &                         , iteration_uramp,iteration_scdft
  use m_Parallelization,    only : mype
  use m_Control_Parameters, only : ipriforce, iprivelocity &
       &                         , m_CtrlP_flag_wd_force &
       &                         , m_CtrlP_set_flag_wd_force_1 &
       &                         , sw_optimize_lattice &
       &                         , sw_cif_output, xctype, printable, driver &
       &                         , sw_nnp_output, m_CtrlP_edelta_for_sampling, frequency_nnp &
       &                         , filetype_nnp
  use m_Force,              only : forcmx, forc_l &
       &                         , m_Force_wd_force_and_cps &
       &                         , m_Force_wd_force_cps_cpd &
!!$       &                         , m_Force_cal_absforc &
       &                         , m_Force_wd_force_and_cps_lim &
       &                         , m_Force_wd_force_cps_cpd_lim
  use m_Ionic_System,       only : ekina,ega,almda,mdmode &
       &                         , forcmx_constraint_quench, cps,cpd_l,natm &
       &                         , constraints_exist, forcmx_hyperplane_vert &
       &                         , m_IS_cnstr_is_fixed_nhp &
       &                         , m_IS_wd_speciesname_etc &
       &                         , m_IS_dump_cif &
       &                         , m_IS_dump_nnpfiles &
       &                         , m_IS_rigid_body_exists, torque_max, trans_force_max
  use m_Total_Energy,       only : etotal
  use m_Stress,             only : m_Stress_get_stressmx
  use m_XC_Potential,       only : ecor
  use m_UnitCell,           only : m_UC_get_hamiltonian,m_UC_get_curr_pressure,m_UC_get_curr_temperature

  implicit none
  integer, intent(in) :: mdalg
  logical, intent(in) :: flag_dynm

  integer, parameter :: N_FORC = 10

  integer :: flag_wd_force

  character :: cion*256,ctot*256
  integer, parameter :: big_ion=100000, big_tot=big_ion*1000
  integer :: iketa

  if(iteration_ionic<big_ion) then
     cion='i5'
  else
     iketa = int(floor(log10(dble(iteration_ionic))))+2
     write(cion,*) iketa
     cion = 'i'//trim(adjustl(cion))
  endif

  if(iteration<big_tot) then
     ctot='i8'
  else
     iketa = int(floor(log10(dble(iteration))))+2
     write(ctot,*) iketa
     ctot = 'i'//trim(adjustl(ctot))
  endif

  if(flag_dynm) then
     flag_wd_force = m_CtrlP_flag_wd_force()
     call wd_forces_and_cps(flag_wd_force)
     if(flag_wd_force == 0) call m_CtrlP_set_flag_wd_force_1()
  else
     call wd_etotal_etc()
  end if

contains
    subroutine wd_etotal_etc
    integer        :: wd_flag = 0
    real(kind=DP)  :: econst
    real(kind=DP)  :: stressmx,temp
    real(kind=DP)  :: big=10000000.d0
    if(wd_flag == 0) then
       if( m_IS_rigid_body_exists()) then

       if(mdalg == VERLET .or. mdalg == TEMPERATURE_CONTROL) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter_total, etotal, ekina, econst, forcmx, transmx, torqmx'
       else if(mdalg == BLUEMOON) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter,etotal,ekina,econst,almda,forcmx'
       else if(mdalg == P_CONTROL .or. mdalg == PT_CONTROL)then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter_total, etotal, ekina, econst, forcmx, pressure, transmx, torqmx'
       else if(driver == DRIVER_URAMP) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_uramp, iter_ion, iter_total, etotal, forcmx, transmx, torqmx'
       else if(driver == DRIVER_SC_DFT) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_scdft, iter_ion, iter_total, etotal, forcmx, transmx, torqmx'
       else
          if(sw_optimize_lattice==OFF)then
          if(mype == 0) write(nfenf,*) ' iter_ion, iter_total, etotal, forcmx, transmx, torqmx'
          else
          if(mype == 0) write(nfenf,*) ' iter_unitcell, iter_ion, iter_total, etotal, forcmx, stressmx, transmx, torqmx'
          endif
       end if

       else

       if(mdalg == VERLET .or. mdalg == TEMPERATURE_CONTROL) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter_total, etotal, ekina, econst, forcmx'
       else if(mdalg == BLUEMOON) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter,etotal,ekina,econst,almda,forcmx'
       else if(mdalg == P_CONTROL .or. mdalg == PT_CONTROL)then
          if(mype == 0) write(nfenf,*) &
               & ' iter_ion, iter_total, etotal, ekina, econst, forcmx, pressure'
       else if(driver == DRIVER_URAMP) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_uramp, iter_ion, iter_total, etotal, forcmx'
       else if(driver == DRIVER_SC_DFT) then
          if(mype == 0) write(nfenf,*) &
               & ' iter_scdft, iter_ion, iter_total, etotal, forcmx'
       else
          if(sw_optimize_lattice==OFF)then
          if(mype == 0) write(nfenf,*) ' iter_ion, iter_total, etotal, forcmx'
          else
          if(mype == 0) write(nfenf,*) ' iter_unitcell, iter_ion, iter_total, etotal, forcmx, stressmx'
          endif
       end if

       endif
       wd_flag = 1
    end if
    if(xctype=='ggapbey')then
       if(printable) then
          write(nfout,'(a)') ' !** subtracting the correlation energy from the total energy'
          write(nfout,'(a)') ' !** since xctype is ggapbey '
       endif
       etotal = etotal-ecor
    endif

    if(m_IS_rigid_body_exists()) then

    if(mdalg == VERLET .or. mdalg == TEMPERATURE_CONTROL) then
       if(mdalg == TEMPERATURE_CONTROL) then
          econst = etotal + ega
       else
          econst = etotal + ekina
       end if
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',6f20.10,"  converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, ekina, econst,  forcmx, trans_force_max, torque_max
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                  & '," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10"  converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, ekina, econst,  forcmx, trans_force_max, torque_max
             endif
          else
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',6f20.10)') iteration_ionic, iteration &
                  &, etotal, ekina, econst,  forcmx, trans_force_max, torque_max
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                  & '," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') iteration_ionic, iteration &
                  &, etotal, ekina, econst,  forcmx, trans_force_max, torque_max
             endif
          end if
       end if
    else if(mdalg == BLUEMOON) then
       econst = etotal + ega
       if(mype == 0) then
            if(abs(etotal)<big)then
            write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',7f20.8)') iteration_ionic, iteration &
            &, etotal, ekina, econst, almda, forcmx, trans_force_max, torque_max
            else
            write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
            & '," ",f0.8," ",f0.8," ",f0.8," ",f0.8," ",f0.8," ",f0.8," ",f0.8)') iteration_ionic, iteration &
            &, etotal, ekina, econst, almda, forcmx, trans_force_max, torque_max
            endif
       endif
    else if(mdalg == QUENCHED_CONSTRAINT) then
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10," converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, forcmx_constraint_quench, trans_force_max, torque_max
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                  &'," ",f0.10," ",f0.10," ",f0.10," ",f0.10," converged_in_a_plane")') iteration_ionic, iteration &
                  &, etotal, forcmx_constraint_quench, trans_force_max, torque_max
             endif
          else
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') iteration_ionic, iteration &
                  &, etotal, forcmx_constraint_quench, trans_force_max, torque_max
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
                  & iteration_ionic, iteration, etotal, forcmx_constraint_quench, trans_force_max, torque_max
             endif
          end if
       end if
    else if (mdalg == PT_CONTROL .or. mdalg == P_CONTROL) then
        econst = m_UC_get_hamiltonian()
        stressmx = m_UC_get_curr_pressure()
        temp = m_UC_get_curr_temperature()
        if(mype == 0)then
        if(abs(etotal)<big)then
        write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',7f20.10)') iteration_unit_cell, iteration &
             &, etotal, temp, econst, forcmx, stressmx, trans_force_max, torque_max
        else
        write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
        &     '," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
        &     iteration_unit_cell, iteration, etotal, temp, econst, forcmx, stressmx, trans_force_max, torque_max
        endif
        endif
    else if (driver == DRIVER_URAMP) then
       if(mype == 0)then
          if(abs(etotal)<big)then
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') &
               & iteration_uramp,iteration_ionic, iteration &
               &, etotal, forcmx, trans_force_max, torque_max
          else
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
               & '," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') iteration_uramp,iteration_ionic, iteration &
               &, etotal, forcmx, trans_force_max, torque_max
          endif
       endif
    else if (driver == DRIVER_SC_DFT) then
       if(mype == 0)then
          if(abs(etotal)<big)then
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') &
               & iteration_scdft,iteration_ionic, iteration, etotal, forcmx, trans_force_max, torque_max
          else
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
               & '," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') iteration_scdft,iteration_ionic, iteration &
               &, etotal, forcmx, trans_force_max, torque_max
          endif
       endif
    else
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10," converged_in_a_plane")') &
                  &  iteration_ionic, iteration &
                  &, etotal, forcmx, trans_force_max, torque_max
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//&
                  &  '," ",f0.10," ",f0.10," ",f0.10," ",f0.10," converged_in_a_plane")') iteration_ionic, iteration &
                  &, etotal, forcmx, trans_force_max, torque_max
             endif
          else
             if(sw_optimize_lattice==OFF)then
!  --> T. Yamasaki, 18 July 2008
             if(constraints_exist .and. m_IS_cnstr_is_fixed_nhp() &
                  & .and. forcmx_hyperplane_vert < forcmx ) then
                if(abs(etotal)<big)then
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') iteration_ionic, iteration &
                     &, etotal, forcmx_hyperplane_vert, trans_force_max, torque_max
                else
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
                     &  iteration_ionic, iteration, etotal, forcmx_hyperplane_vert, trans_force_max, torque_max
                endif
             else
                if(abs(etotal)<big)then
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') iteration_ionic, iteration &
                     &, etotal, forcmx, trans_force_max, torque_max
                else
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
                     &  iteration_ionic, iteration, etotal, forcmx, trans_force_max, torque_max
                endif
             end if
!  <--
             else
                stressmx = m_Stress_get_stressmx()
                if(stressmx>=0) then
                  if(abs(etotal)<big)then
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',5f20.10)') &
                       & iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx,stressmx, trans_force_max, torque_max
                  else
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                       &  '," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
                       &  iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx,stressmx, trans_force_max, torque_max
                  endif
                else
                  if(abs(etotal)<big)then
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)')  &
                       &  iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx, trans_force_max, torque_max
                  else
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                       &  '," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') iteration_unit_cell, &
                       &  iteration_ionic, iteration, etotal, forcmx, trans_force_max, torque_max
                  endif
                endif
             endif
          end if
       end if
    end if

    else

    if(mdalg == VERLET .or. mdalg == TEMPERATURE_CONTROL) then
       if(mdalg == TEMPERATURE_CONTROL) then
          econst = etotal + ega
       else
          econst = etotal + ekina
       end if
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10,"  converged_in_a_plane")') &
                  &  iteration_ionic, iteration &
                  &, etotal, ekina, econst,  forcmx
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//&
                  &  '," ",f0.10," ",f0.10," ",f0.10," ",f0.10,"  converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, ekina, econst,  forcmx
             endif
          else
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',4f20.10)') iteration_ionic, iteration &
                  &, etotal, ekina, econst,  forcmx
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
                  &  iteration_ionic, iteration, etotal, ekina, econst,  forcmx
             endif
          end if
       end if
    else if(mdalg == BLUEMOON) then
       econst = etotal + ega
       if(mype == 0) then
            if(abs(etotal)<big)then
            write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',5f20.8)') iteration_ionic, iteration &
            &, etotal, ekina, econst, almda, forcmx
            else
            write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.8," ",f0.8," ",f0.8," ",f0.8," ",f0.8)') &
            &  iteration_ionic, iteration, etotal, ekina, econst, almda, forcmx
            endif
       endif
    else if(mdalg == QUENCHED_CONSTRAINT) then
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10," converged_in_a_plane")') &
                  &  iteration_ionic, iteration, etotal, forcmx_constraint_quench
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," converged_in_a_plane")') &
                  &  iteration_ionic, iteration &
                  &, etotal, forcmx_constraint_quench
             endif
          else
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') iteration_ionic, iteration &
                  &, etotal, forcmx_constraint_quench
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10)') &
                  &  iteration_ionic, iteration, etotal, forcmx_constraint_quench
             endif
          end if
       end if
    else if (mdalg == PT_CONTROL .or. mdalg == P_CONTROL) then
        econst = m_UC_get_hamiltonian()
        stressmx = m_UC_get_curr_pressure()
        temp = m_UC_get_curr_temperature()
        if(mype == 0)then
        if(abs(etotal)<big)then
        write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',5f20.10)') iteration_unit_cell, iteration &
             &, etotal, temp, econst, forcmx, stressmx
        else
        write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," ",f0.10," ",f0.10," ",f0.10)') &
        &     iteration_unit_cell, iteration, etotal, temp, econst, forcmx, stressmx
        endif
        endif
    else if (driver == DRIVER_URAMP) then
       if(mype == 0)then
          if(abs(etotal)<big)then
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') &
               &  iteration_uramp,iteration_ionic, iteration, etotal, forcmx
          else
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10)') &
               & iteration_uramp,iteration_ionic, iteration, etotal, forcmx
          endif
       endif
    else if (driver == DRIVER_SC_DFT) then
       if(mype == 0)then
          if(abs(etotal)<big)then
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') &
               & iteration_scdft,iteration_ionic, iteration, etotal, forcmx
          else
          write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10)') &
               & iteration_scdft,iteration_ionic, iteration, etotal, forcmx
          endif
       endif
    else
       if(mype == 0) then
          if(mdmode == CNSTRA) then
             if(abs(etotal)<big)then
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10," converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, forcmx
             else
             write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10," converged_in_a_plane")') &
                  & iteration_ionic, iteration, etotal, forcmx
             endif
          else
             if(sw_optimize_lattice==OFF)then
!  --> T. Yamasaki, 18 July 2008
             if(constraints_exist .and. m_IS_cnstr_is_fixed_nhp() &
                  & .and. forcmx_hyperplane_vert < forcmx ) then
                if(abs(etotal)<big)then
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') iteration_ionic, iteration &
                     &, etotal, forcmx_hyperplane_vert
                else
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10)') &
                     & iteration_ionic, iteration, etotal, forcmx_hyperplane_vert
                endif
             else
                if(abs(etotal)<big)then
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') iteration_ionic, iteration &
                     &, etotal, forcmx
                else
                write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ",f0.10," ",f0.10)') &
                     & iteration_ionic, iteration, etotal, forcmx
                endif
             end if
!  <--
             else
                stressmx = m_Stress_get_stressmx()
                if(stressmx>=0) then
                  if(abs(etotal)<big)then
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',3f20.10)') &
                       & iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx,stressmx
                  else
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                       & '," ",f0.10," ",f0.10," ",f0.10)') iteration_unit_cell, iteration_ionic, iteration &
                       &, etotal, forcmx,stressmx
                  endif
                else
                  if(abs(etotal)<big)then
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))//',2f20.10)') &
                       & iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx
                  else
                  write(nfenf,'(" ",'//trim(adjustl(cion))//','//trim(adjustl(cion))//','//trim(adjustl(ctot))// &
                       & '," ",f0.10," ",f0.10)') iteration_unit_cell, iteration_ionic, iteration, etotal, forcmx
                  endif
                endif
             endif
          end if
       end if
    end if

    endif
    call m_Files_flush_nfenf()

    if(printable) write(nfout,'(a,f20.10,a)') '!** total energy for the current iteration : ',etotal,' hartree'

  end subroutine wd_etotal_etc

  subroutine wd_forces_and_cps(flag_wd_force)
    use m_Control_Parameters, only : sw_stress_correction, imdalg
!
! (iprivelocity)
!    Modified by T. Yamasaki, April/10/2007
!
    integer, intent(in) :: flag_wd_force
    character(len=256) :: ch
!!$    real(kind=DP) :: forc_lower

    if(flag_wd_force == 0 .and. sw_stress_correction /= ON) then
       call m_IS_wd_speciesname_etc(nfdynm)
    end if

    if(mype == 0 .and. iprivelocity <= 1) then
       if(mdmode == CNSTRA) then
          write(nfdynm,'(" cps and forc at (iter_ion, iter_total = " &
               & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ) converged_in_a_plane")') iteration_ionic, iteration
       else
          write(nfdynm,'(" cps and forc at (iter_ion, iter_total = " &
               & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," )")') iteration_ionic, iteration
       end if
    else if(mype == 0 .and. iprivelocity >= 2) then
       if(mdmode == CNSTRA) then
          write(nfdynm,'(" cps cpd and forc at (iter_ion, iter_total = " &
               & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," ) converged_in_a_plane")') iteration_ionic, iteration
       else
          if(.not.(imdalg==PT_CONTROL .or. imdalg==P_CONTROL)) then
            write(nfdynm,'(" cps cpd and forc at (iter_ion, iter_total = " &
                 & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," )")') iteration_ionic, iteration
          else
            write(nfdynm,'(" cps cpd and forc at (iter_ion, iter_total = " &
                 & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," )")') iteration_unit_cell, iteration
          endif
       end if
    end if

    if(iprivelocity <= 1.and.mype==0) then
       call m_Force_wd_force_and_cps(nfdynm, WITHOUTTAG, cps, natm)
    else if(iprivelocity >= 2.and.mype==0) then
       call m_Force_wd_force_cps_cpd(nfdynm, WITHOUTTAG,cps,cpd_l,natm)
    end if

    if(mype==0.and.sw_cif_output==ON)then
       write(ch,*) iteration_ionic
       call m_IS_dump_cif(nfdynm_cif,'iteration_ionic_'//trim(adjustl(ch)))
    endif

    if(sw_nnp_output == ON) then
       call m_CtrlP_edelta_for_sampling(nfout, iteration_ionic)
       if(mod(iteration_ionic,frequency_nnp)==0) then
          call m_IS_dump_nnpfiles(etotal, forc_l)
       endif
    endif

    call m_Files_flush_nfdynm()

    if(ipriforce >= 1) then
       if(mype == 0 .and. iprivelocity <= 1) then
          write(nfout,'(" !forc cps and forc at (iter_ion, iter_total = " &
               & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," )")') iteration_ionic, iteration
       else if(mype == 0 .and. iprivelocity >= 2) then
          write(nfout,'(" !forc cps cpd and forc at (iter_ion, iter_total = " &
               & ,'//trim(adjustl(cion))//','//trim(adjustl(ctot))//'," )")') iteration_ionic, iteration
       end if
    end if
#ifdef _DEBUG_WRITE_
    if(ipriforce >= 1) then
#else
    if(ipriforce == 1 .and. natm > N_FORC) then
!!$       call m_Force_cal_absforc(N_FORC, forc_lower)
!!$       write(nfout,'(" forc_lower = ",e12.4)') forc_lower
       if(iprivelocity <= 1) then
          call m_Force_wd_force_and_cps_lim(nfout, WITHTAG, cps, natm, N_FORC)
       else if(iprivelocity >= 2) then
          call m_Force_wd_force_cps_cpd_lim(nfout, WITHTAG, cps,cpd_l,natm, N_FORC)
       end if
    else if(ipriforce >= 2 .or. (ipriforce == 1 .and. natm <= N_FORC) ) then
#endif
       if(iprivelocity <= 1) then
          call m_Force_wd_force_and_cps(nfout, WITHTAG, cps, natm)
       else if(iprivelocity >= 2) then
          call m_Force_wd_force_cps_cpd(nfout, WITHTAG, cps,cpd_l,natm)
       end if
    end if

  end subroutine wd_forces_and_cps

end subroutine wd_forces_cps_etotal_and_etc

#ifndef _EMPIRICAL_
logical function Force_errors_are_tolerable()
  use m_Const_Parameters,   only : OFF, ON
  use m_Files,              only : nfout
  use m_Control_Parameters, only : force_error_check_mode &
       &                         , force_error_check_rangeL, sw_phonon
!!$       &                         , m_CtrlP_what_is_mdalg
  use m_Force,              only : forcmx, m_Force_error_tolerable &
       &                         , m_Force_clear_fec_counter

  implicit none
!!$  integer :: mdalg

!!$  if(force_error_check_mode == OFF .or. forcmx < force_error_check_rangeL) then
  if(sw_phonon == ON .or. force_error_check_mode == OFF .or. forcmx < force_error_check_rangeL) then
     Force_errors_are_tolerable = .true.
     ! --> T. Yamasaki, August 8 2008
     call m_Force_clear_fec_counter()
     ! <--
  else
     if(m_Force_error_tolerable(nfout)) then
        Force_errors_are_tolerable = .true.
     else
        Force_errors_are_tolerable = .false.
     end if
  end if

!!$  mdalg = m_CtrlP_what_is_mdalg()
!!$  if(Force_errors_are_tolerable) call wd_forces_cps_etotal_and_etc(mdalg)
end function Force_errors_are_tolerable
#endif

logical function Already_Converged()
  use m_Control_Parameters, only : iconvergence_previous_job &
      &                          , sw_phonon, sw_calc_force, max_total_scf_iteration &
      &                          , terminated_because, icond
  use m_Const_Parameters, only :   FORCE_CONVERGED, FORCE_CONVERGENCE_REACHED &
      &                          , MAX_PHSTEPS_REACHED, WF_CONVERGENCE_REACHED &
      &                          , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION
  implicit none
  if(iconvergence_previous_job >= FORCE_CONVERGED) then
     Already_Converged = .true.
  else
     Already_Converged = .false.
  end if
  if(sw_phonon == 1 .and. sw_calc_force == 0) Already_Converged = .true.
  if(max_total_scf_iteration <= 0) Already_Converged = .true.
  if(Already_Converged) then
    if(sw_phonon==1) then
      terminated_because = MAX_PHSTEPS_REACHED
    else if (icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) then
      terminated_because = WF_CONVERGENCE_REACHED
    else
      terminated_because = FORCE_CONVERGENCE_REACHED
    endif
  endif

end function Already_Converged

logical function MultiReplicaMode()
  use m_Control_Parameters, only  : multiple_replica_mode
  use m_Ionic_System,       only  : number_of_replicas
  use m_Const_Parameters, only    : ON
  if(multiple_replica_mode == ON .and. number_of_replicas >= 2) then
     MultiReplicaMode = .true.
  else
     MultiReplicaMode = .false.
  end if
end function MultiReplicaMode

#ifndef _EMPIRICAL_
logical function EigenValues_are_Converged()
  use m_Const_Parameters,     only : EK, SCF,GRID, ON, YES, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                           , ONE_BY_ONE, OFF, WF_CONVERGENCE_REACHED
  use m_Control_Parameters,   only : ekmode,nspin,neg,max_total_scf_iteration,ipri,ipriekzaj,iprievdff &
       &                           , ek_max_iteration, sw_berry_phase, sw_ekzaj, printable, icond &
       &                           , fixed_charge_k_parallel, terminated_because, sw_dos, sw_ldos &
       &                           , sw_charge_rspace, sw_partial_charge, sw_pdos
  use m_IterationNumbers,     only : iteration_electronic, iteration, nk_converged &
       &                           , m_Iter_set_converged_nk, nkgroup, nk_in_the_process
  use m_Files,                only : nfout, nfeng, nfzaj &
       &                           , m_Files_open_nfzaj_with_check
  use m_Kpoints,              only : kv3_ek, kv3
  use m_Parallelization,      only : nrank_k, npes, ierr, mype, MPI_CommGroup
  use m_Electronic_Structure, only : m_ES_eekdif, m_ES_eekdif2 &
       &                           , m_ES_cp_eko_l_to_eko_ek &
       &                           , m_ES_cp_eko_l_to_eko_ek2 &
       &                           , m_ES_wd_eko &
       &                           , m_ES_sym_comp
  use m_ES_nonlocal,          only : m_ES_phir_dot_WFs

  use m_ES_IO,                only : m_ESIO_wd_WFs_and_EVs_ek
  use m_ES_IO,                only : m_ESIO_wd_EigenValues
  use m_BerryPhase,           only : m_BP_write_wfbp, m_BP_write_wfbp_gshift
  use mpi

  implicit none
!  include 'mpif.h'                                      ! MPI
  integer :: ipri0, kv3_t

  if(ekmode == OFF .and. &
       & (icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) &
       & .and. fixed_charge_k_parallel == ONE_BY_ONE) then
     EigenValues_are_Converged = m_ES_eekdif2()
!!$     if(nrank_k >= 2) then
!!$        call m_ESIO_wd_EigenValues(nfout,iprievdff,nooccupation=YES)
!!$     end if
  else
     EigenValues_are_Converged = m_ES_eekdif()
  end if

  if(sw_ekzaj==ON .and. ekmode==ON ) EigenValues_are_Converged = .true.

  if(EigenValues_are_Converged) then
!!$     if(nk_in_the_process >= kv3_ek) stop ' nk_in_the_process >= kv3_ek <<Eigenvalues_are_converged>>'
     if(printable) write(nfout,'(" Wave Functions have converged at ",i7," -th iteration")') iteration_electronic

! ===================== modified by K. Tagami ========================= 12.0Exp
!     if(nk_converged == 0 .and. mype==0) then
!!!$     write(nfeng,'(" nk_converged = ",i6)') nk_converged
!        write(nfeng,'(" num_kpoints = ",i6)') kv3_ek
!        write(nfeng,'(" num_bands   = ",i6)') neg
!        write(nfeng,'(" nspin       = ",i6)') nspin
!        write(nfeng,'(" Fermi energy level = "/)')
!     end if
!
     terminated_because = WF_CONVERGENCE_REACHED
     if ( fixed_charge_k_parallel == ONE_BY_ONE ) then
        if (nk_converged == 0 .and. mype==0) then
           write(nfeng,'(" num_kpoints = ",i6)') kv3_ek
           write(nfeng,'(" num_bands   = ",i6)') neg
           write(nfeng,'(" nspin       = ",i6)') nspin
           write(nfeng,'(" Fermi energy level = "/)')
        end if
     else
        if (nk_converged == 0 .and. mype==0) then
           write(nfeng,'(" num_kpoints = ",i6)') kv3
           write(nfeng,'(" num_bands   = ",i6)') neg
           write(nfeng,'(" nspin       = ",i6)') nspin
           write(nfeng,'(" Fermi energy level = "/)')
        end if
     endif
! ==================================================================== 12.0Exp

     if(ekmode == ON) then
! === DEBUG by tkato 2013/10/18 ================================================
!       call m_ES_cp_eko_l_to_eko_ek()
        call m_ES_cp_eko_l_to_eko_ek2()
! ==============================================================================
     else if((icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) .and. &
           fixed_charge_k_parallel == ONE_BY_ONE) then
        call m_ES_cp_eko_l_to_eko_ek2()
     end if

     if(nrank_k >= 2) then
        call m_ESIO_wd_EigenValues(nfout,2,nooccupation=YES)
        call m_ESIO_wd_EigenValues(nfeng,2,nooccupation=YES)
     else
        call m_ES_wd_eko(nfout,mode=EK)
        call m_ES_wd_eko(nfeng,mode=EK)
     end if
!!$     if(nk_in_the_process >= kv3_ek) stop ' nk_in_the_process >= kv3_ek (2)<<Eigenvalues_are_converged>>'
  else if(iteration_electronic >= ek_max_iteration) then
!!$       & .or. iteration >= max_total_scf_iteration ) then
     if(ekmode == ON) then
! === DEBUG by tkato 2013/10/18 ================================================
!       call m_ES_cp_eko_l_to_eko_ek()
!       write(nfout,'(" m_ES_cp_eko_l_to_eko_ek")')
        call m_ES_cp_eko_l_to_eko_ek2()
! ==============================================================================
     else if((icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) .and. &
           fixed_charge_k_parallel == ONE_BY_ONE) then
        call m_ES_cp_eko_l_to_eko_ek2()
!!$        write(nfout,'(" m_ES_cp_eko_l_to_eko_ek2")')
     end if
     if(mype == 0) then
        write(nfout,'(" Wave Functions have not converged")')
        write(nfout,'(" iteration_electronic reached to ek_max_iteration (= ",i5," )")') ek_max_iteration
     end if
     if(nrank_k >= 2) then
        call m_ESIO_wd_EigenValues(nfout,2,nooccupation=YES)
        call m_ESIO_wd_EigenValues(nfeng,2,nooccupation=YES)
     else
        call m_ES_wd_eko(nfout,mode=EK)
        call m_ES_wd_eko(nfeng,mode=EK)
     end if
     EigenValues_are_Converged = .true.
!!$     if(iteration >= max_total_scf_iteration) write(nfout,'(" iteration has reached to max_total_scf_iteration")')
     if(iteration_electronic >= ek_max_iteration .and. mype==0) &
      & write(nfout,'(" iteration_electronic has reached to ek_max_iteration")')
  end if

! ============================== added by K. Tagami =============== 12.0Exp
  if ( fixed_charge_k_parallel /= ONE_BY_ONE ) return
! ================================================================= 12.0Exp
  if(EigenValues_are_Converged) then
     if(ipri >= 2) write(nfout,'(" ipriekzaj = ",i8," <<EigenValues_are_Converged>>")') ipriekzaj
     if(ipriekzaj >= 1 .and. ekmode /= GRID .and. sw_ekzaj/=ON) then
        call m_Files_open_nfzaj_with_check()
        call m_ESIO_wd_WFs_and_EVs_ek(nfout,nfzaj)
     end if
     if(sw_berry_phase == 1) then
        call m_BP_write_wfbp_gshift ! for the initial k-point of a strig
        call m_BP_write_wfbp
     end if
!!$     call m_Iter_set_converged_nk(nspin)
     kv3_t = min(kv3, kv3_ek-kv3*(nkgroup-1))
     call m_Iter_set_converged_nk(kv3_t)
  else
     ipri0 = iprievdff
!!$     if(nrank_k>=2) call mpi_bcast(ipri0,1,mpi_integer,0,MPI_CommGroup,ierr)
     if(npes>=2)    call mpi_bcast(ipri0,1,mpi_integer,0,MPI_CommGroup,ierr)
     if((ipri0 >= 1 .and. mod(iteration_electronic,10)==1) &
          .or. ipri0 >= 2 ) then
        if(printable) write(nfout,'(" iteration_electronic = ",i6)') iteration_electronic
        if(nrank_k >= 2) then
           call m_ESIO_wd_EigenValues(nfout,2,nooccupation=YES)
        else
           call m_ES_wd_eko(nfout,mode=EK)
        end if
     end if
  end if
end function EigenValues_are_Converged

logical function AllKpoints_are_Calculated()
  use m_Files,              only : nfout
  use m_Const_Parameters,   only : EK_CONVERGED
  use m_IterationNumbers,   only : nk_converged, m_Iter_nk_incre2
  use m_Control_Parameters, only : nspin, ipri, m_CtrlP_set_iconvergence
  use m_Kpoints,            only : kv3_ek

  AllKpoints_are_Calculated = .false.
  if(ipri >= 2) write(nfout,'(" !!! nk_converged = ",i8)') nk_converged
  if(nk_converged+nspin-1 >= kv3_ek) then
     AllKpoints_are_Calculated = .true.
     if(ipri >= 1) write(nfout,'(" All kpoints have been calculated ")')
  end if
  if(AllKpoints_are_Calculated) then
     call m_CtrlP_set_iconvergence(EK_CONVERGED)
     call m_Iter_nk_incre2(nspin,kv3_ek)
  end if
end function AllKpoints_are_Calculated

logical function AllKpoints_are_Converged()
  use m_Files,              only : nfout
  use m_Const_Parameters,   only : EK_CONVERGED
  use m_IterationNumbers,   only : nk_in_the_process
  use m_Control_Parameters, only : ipri, m_CtrlP_set_iconvergence
  use m_Kpoints,            only : kv3_ek

  AllKpoints_are_Converged = .false.
  if(nk_in_the_process > kv3_ek) then
     AllKpoints_are_Converged = .true.
     if(ipri >= 1) write(nfout,'(" All Kpoints have been converged")')
  end if
  if(AllKpoints_are_Converged) call m_CtrlP_set_iconvergence(EK_CONVERGED)
end function AllKpoints_are_Converged

logical function AllKpoints_are_Calculated2(nk)
  use m_Files,              only : nfout
  use m_Const_Parameters,   only : EK_CONVERGED
  use m_IterationNumbers,   only : nk_converged, nk_in_the_process, m_Iter_nk_incre2
  use m_Control_Parameters, only : ipri, nspin, m_CtrlP_set_iconvergence
  use m_Kpoints,            only : kv3_ek, kv3

  integer, intent(inout) :: nk

  nk = nk + kv3
  if(ipri >= 1) write(nfout,'(" nk, nk_in_the_process = ",2i8)') nk,nk_in_the_process
!!$  if(nk >= kv3_ek) then
!!$  if(nk_in_the_process >= kv3_ek) then
!!$  if(nk_in_the_process+kv3 >= kv3_ek) then
  if(nk_converged >= kv3_ek) then
     AllKpoints_are_Calculated2 = .true.
     call m_Iter_nk_incre2(nspin,kv3_ek) !-> nk_in_the_process
  else
     AllKpoints_are_Calculated2 = .false.
  end if
  if(AllKpoints_are_Calculated2) call m_CtrlP_set_iconvergence(EK_CONVERGED)
end function AllKpoints_are_Calculated2

subroutine Set_converged_nk()
  use m_IterationNumbers,   only : m_Iter_set_converged_nk
  use m_Control_Parameters, only : nspin
  call m_Iter_set_converged_nk(nspin)
end subroutine Set_converged_nk

logical function Already_Converged2()
  ! coded by T. Yamasaki, 31 May 2007
  use m_Const_Parameters, only : FORCE_CONVERGED
  use m_Control_Parameters, only : iconvergence, iconvergence_previous_job
  if(iconvergence >= FORCE_CONVERGED .or. iconvergence_previous_job >= FORCE_CONVERGED) then
     Already_Converged2 = .true.
  else
     Already_Converged2 = .false.
  end if
end function Already_Converged2

logical function Already_Converged_for_kgroup()
  use m_Kpoints, only           : kv3_ek, kv3
  use m_IterationNumbers, only  : nkgroup, nk_in_the_process
  use m_Control_Parameters,only : ipri,nspin
  use m_Const_Parameters, only  : EK_CONVERGED
  use m_Electronic_Structure, only : iconv_ek
  use m_Files, only             : nfout
  integer :: i
  Already_Converged_for_kgroup = .true.
  do i = nk_in_the_process, min(nk_in_the_process + kv3, kv3_ek)
     if(iconv_ek(i) < EK_CONVERGED) then
        Already_Converged_for_kgroup = .false.
        exit
     end if
  end do
  if(ipri >= 2 .and. .not.Already_Converged_for_kgroup) then
     write(nfout,'(" -- Convergence status for ek --")')
     write(nfout,'(40i2)') (iconv_ek(i),i=nk_in_the_process, min(nk_in_the_process+kv3,kv3_ek))
  end if
!!$  if(Already_Converged_for_kgroup) call m_Iter_nkgroup_incre() !-> nk_in_the_process
end function Already_Converged_for_kgroup

logical function Positron_Bulk()
  ! Coded by T. Yamasaki, 29 Oct. 2003
  use m_Const_Parameters, only :   BULK
  use m_Control_Parameters, only : sw_positron
  if(sw_positron == BULK) then
     Positron_Bulk = .true.
  else
     Positron_Bulk = .false.
  end if
end function Positron_Bulk

logical function Positron_Defect()
  ! Coded by T. Yamasaki, 11 Dec. 2005
  use m_Const_Parameters, only :   DEFECT
  use m_Control_Parameters, only : sw_positron
  if(sw_positron == DEFECT) then
     Positron_Defect = .true.
  else
     Positron_Defect = .false.
  end if
end function Positron_Defect

! ==== Positron SCF === 2015/11/28
logical function Positron_scf()
  use m_Const_Parameters, only :   OFF, Positron_CONV
  use m_Control_Parameters, only : sw_positron, positron_method

  if ( sw_positron /= OFF ) then
     if ( positron_method == Positron_CONV ) then
        Positron_scf = .false.
     else
        Positron_scf = .true.
     end if
  else
     Positron_scf = .false.
  endif
end function Positron_scf

logical function Positron_nonscf()
  use m_Const_Parameters, only :   OFF, Positron_CONV
  use m_Control_Parameters, only : sw_positron, positron_method

  if ( sw_positron /= OFF ) then
     if ( positron_method == Positron_CONV ) then
        Positron_nonscf = .true.
     else
        Positron_nonscf = .false.
     end if
  else
     Positron_nonscf = .false.
  endif
end function Positron_nonscf
! ========== 2015/11/28

logical function Structure_is_fixed()
  ! Coded by T. Yamasaki, 25 Jul. 2008
  use m_Const_Parameters, only :   ON, VERLET
  use m_Control_Parameters, only : sw_fix &
       &                         , m_CtrlP_set_iconvergence &
       &                         , m_CtrlP_what_is_mdalg
  use m_Const_Parameters,   only : FORCE_CONVERGED
  use m_Ionic_System,       only : m_IS_evaluate_v_verlet
  use m_Force,              only : forc_l
  integer :: mdalg

  if(sw_fix == ON) then
     Structure_is_fixed = .true.
  else
     Structure_is_fixed = .false.
  end if
  if(Structure_is_fixed) call m_CtrlP_set_iconvergence(FORCE_CONVERGED)
  if(Structure_is_fixed) then
     mdalg = m_CtrlP_what_is_mdalg()
     if(mdalg == VERLET) call m_IS_evaluate_v_verlet(mdalg,forc_l)
     call wd_forces_cps_etotal_and_etc(mdalg,.true.) ! nfdynm
     call wd_forces_cps_etotal_and_etc(mdalg,.false.) ! nfenf
  end if
end function Structure_is_fixed

logical function Hubbard_model()
  use m_Const_Parameters, only : ON
  use m_Control_Parameters, only : sw_hubbard
  use m_Total_Energy,   only : m_TE_Converged_Hubbard_Energy
  implicit none
  Hubbard_model=.false.
  if(sw_hubbard == ON) then
     Hubbard_model=.true.
     !!if(m_TE_Converged_Hubbard_Energy(nfout)) Hubbard_model=.false.
  end if
  !!$write(6,*) "Hubbard_model=",Hubbard_model
end function Hubbard_model
#endif

subroutine Checkpoint_File(mode)
  use m_Const_Parameters, only : DP,STREVL_ITERATION,SCF_ITERATION,UNITCELL_ITERATION,REAC_ITERATION, DRIVER_CONSTRAINT
  use m_Control_Parameters, only :  cpt_scf_iteration, cpt_nhistory, cpt_strevl_iteration, cpt_reac_iteration &
                                 &, cpt_time, cpt_unitcell_iteration, m_CtrlP_get_elpsd_time, printable, driver
  use m_IterationNumbers, only : iteration, iteration_ionic, iteration_unit_cell
  use m_Files, only : nfout
  use m_Parallelization, only : mype,MPI_CommGroup
  use m_constraints, only : m_cnstr_get_id
  implicit none
  integer, intent(in) :: mode
  real(kind=DP) :: curr_time
  real(kind=DP), save :: last_time=0.d0
  logical :: write_chkpnt,lscf,lstrevl,ltime,lunit,lreac
  integer :: ierr
  write_chkpnt = .false.;lscf=.false.;lstrevl=.false.;ltime=.false.;lunit=.false.;lreac=.false.
  if (mode==SCF_ITERATION .and. cpt_scf_iteration>0) then
  if (mod(iteration,cpt_scf_iteration) == 0) then
     write_chkpnt = .true.
     lscf = .true.
  endif
  endif
  if (mode==STREVL_ITERATION .and. cpt_strevl_iteration>0) then
  if (mod(iteration_ionic,cpt_strevl_iteration) == 0) then
     write_chkpnt = .true.
     lstrevl = .true.
  endif
  endif
  if (mode==UNITCELL_ITERATION .and. cpt_unitcell_iteration>0) then
  if (mod(iteration_unit_cell,cpt_unitcell_iteration) == 0) then
     write_chkpnt = .true.
     lunit = .true.
  endif
  endif
  if (mode==REAC_ITERATION .and. cpt_reac_iteration>0) then
  if (mod(m_cnstr_get_id(),cpt_reac_iteration) == 0) then
     write_chkpnt = .true.
     lreac = .true.
  endif
  endif
  curr_time = m_CtrlP_get_elpsd_time()
  if (cpt_time>0 .and. (curr_time-last_time)>cpt_time) then
     write_chkpnt = .true.
     ltime = .true.
     last_time = curr_time
     call mpi_barrier(MPI_CommGroup,ierr)
  endif
  if(write_chkpnt) then
     call WriteDownData_onto_Files(.false.)
     if(driver == DRIVER_CONSTRAINT) call constrained_dynamics_dump()
     if(printable) then
        write(nfout,'(a)')              '!** dumped checkpoint files because '
        if (lscf)    write(nfout,'(a)') '!** iteration'
        if (lstrevl) write(nfout,'(a)') '!** iteration_ionic'
        if (lunit)   write(nfout,'(a)') '!** iteration_unitcell'
        if (lreac)   write(nfout,'(a)') '!** iteration_reac'
        if (ltime)   write(nfout,'(a)') '!** cputime'
        write(nfout,'(a)')              '!** met the criterion'
     endif
  endif
end subroutine Checkpoint_File

subroutine Modify_Input_Settings()
  use m_Files, only : m_Files_reopen_nfinp_mod, m_Files_close_nfinp_mod &
     &       , m_Files_file_exists, F_INP_MOD, nfout
  use m_Ionic_System, only : m_IS_reread_imdtyp, natm2
  use m_Control_Parameters, only : m_CtrlP_reread_edelta &
     &       , m_CtrlP_reread_max_force, m_CtrlP_reread_max_iteration &
     &       , m_CtrlP_rd_wfsolver2, m_CtrlP_rd_chargemix2 &
     &       , m_CtrlP_reread_cutoff_wf
  implicit none
  integer :: iret, f_closeInputFile
  if(.not. m_Files_file_exists(F_INP_MOD)) then
    return
  endif
  call m_Files_reopen_nfinp_mod(1)
  call m_IS_reread_imdtyp(nfout)
  call m_CtrlP_reread_edelta(nfout)
  call m_CtrlP_reread_max_force(nfout)
  call m_CtrlP_reread_max_iteration(nfout)
  call m_CtrlP_rd_wfsolver2(nfout, natm2, .true.)
  call m_CtrlP_rd_chargemix2(nfout, .true.)
  call m_CtrlP_reread_cutoff_wf(nfout)
  iret = f_closeInputFile()
  call m_Files_close_nfinp_mod()
end subroutine Modify_Input_Settings

subroutine Postproc_after_SCF_convergence(nfout)
  use m_Control_Parameters, only : m_CtrlP_gmax_changed
#ifndef FFTW3
  use m_Parallelization, only : mype
#endif
  implicit none
  integer, intent(in) :: nfout

  if(m_CtrlP_gmax_changed()) then
#ifdef FFTW3
    call change_gmax()
#else
    if (mype==0) then
        write(nfout,'(a)') &
        & '!** WARN: changing gmax during execution is supported only if phase was built with FFTW3'
    endif
#endif
  endif

#ifdef FFTW3
  contains

  subroutine change_gmax()
    use m_Const_Parameters,   only : ON, OFF, LDA, EXECUT
    use m_Control_Parameters, only : gmax, gmax_buf, ipriparallel &
    & , sw_communicator_for_chg, paramset, printable, ggacmp_parallel, neg &
    & , m_CtrlP_reset_gmax_changed, noncol
    use m_PlaneWaveBasisSet, only : m_pwBS_wd_curr_pws, m_pwBS_dealloc &
    & , m_pwBS_store_prev_kg1_kgp, m_pwBS_assume_G_rhombohedron &
    & , m_pwBS_generate_G_vectors, m_pwBS_alloc_ngpt_igfp_gr, m_pwBS_calc_length_of_G &
    & , m_pwBS_G_trans_functions, m_pwBS_setup_FFTmapfunctions &
    & , n_rGv,n_rGpv,n_rGpv_reduced, n_rGv_pstrn,kgp, n_rGv_prev &
    & , n_rGpv_prev, kgp, m_pwBS_rd_prev_pws, m_pwBS_for_each_WF &
    & , m_pwBS_alloc_ylm_l, m_pwBS_sphrp_l
    use m_Parallelization,    only : m_Parallel_init_mpi_kngp_prev &
    & ,  np_g1k, m_Parallel_init_mpi_kngp, m_Parallel_init_mpi_cdfft
    use m_PseudoPotential,    only : m_PP_input_xctype, modnrm, m_PP_dealloc &
    & , m_PP_gfqwei, m_PP_set_fqwei_noncl, m_PP_init_epc
    use m_FFT,                only : m_FFT_set_box_sizes, m_FFT_setup &
    & , m_FFT_reset_CD_setup_stat
    use m_Electronic_Structure, only : m_ES_alloc_zaj_l_prev &
    & , m_ES_dealloc_zaj_l_prev, m_ES_cp_zaj_to_zaj_prev  &
#ifdef FFTW3
    & , m_ES_gen_zaj_from_prev_zaj &
#endif
    & , m_ES_dealloc, m_ES_alloc_zaj_etc, m_ES_realloc_zaj &
    & , m_ES_alloc_vlhxc, m_ES_alloc_vlhxcQ
    use m_Kpoints, only : kv3
    use m_Crystal_Structure, only : inversion_symmetry
#ifdef FFTW3
    use m_Charge_Density, only : m_CD_gen_chgq_from_prev_chgq, m_CD_dealloc_chgq, m_CD_alloc_chgq
#else
    use m_Charge_Density, only : m_CD_dealloc_chgq, m_CD_alloc_chgq
#endif
    use m_ES_initialWF, only : m_ESIW_by_randomnumbers
    use m_ES_WF_by_SDorCG, only : m_ESsd_dealloc_dzajn2, m_ESsd_dealloc_zaj_old &
    & , m_ESsd_alloc_zaj_old, m_ESsd_alloc_dzajn2
    use m_XC_Potential, only : m_XC_rst_npsend
    use m_CD_mixing, only : m_CD_force_dealloc
    use m_ES_WF_by_submat, only : m_ESsubmat_dealloc
    use m_ES_WF_by_RMM, only : m_ESrmm_dealloc_r_norm_flag, m_ESrmm_reset_r_norm_flag, m_ESrmm_reset_ng_maxmin
    use m_NonLocal_Potential, only : m_NLP_dealloc

    implicit none

    logical :: rd
    integer :: outer_or_inner, xctype_is, ggacmp_parallel_rev
    call m_pwBS_wd_curr_pws()
    call m_pwBS_rd_prev_pws(rd)
    if(.not. rd) then
      if(printable) write(nfout,'(a)') '!** failed to load pwbs info. gmax &
      &  will not be changed.'
      return
    endif

    call m_XC_rst_npsend()
    call m_FFT_reset_CD_setup_stat()
    call m_CD_force_dealloc()
    call m_ESrmm_reset_r_norm_flag()
    call m_ESrmm_reset_ng_maxmin()
    call m_ESsubmat_dealloc()
    call m_pwBS_store_prev_kg1_kgp()
    call m_pwBS_dealloc(dealloc_all = .true.)
    gmax = gmax_buf
    call m_pwBS_assume_G_rhombohedron()
    call fft_box_finding_way(outer_or_inner)
    call m_FFT_set_box_sizes(n_rGv,n_rGpv_reduced,n_rGv_pstrn,outer_or_inner) ! ->fft_box_size_WF,CD
    call m_pwBS_generate_G_vectors()    ! ->n_rGv,n_rGpv ->kgp
    call m_Parallel_init_mpi_kngp(nfout,ipriparallel,kgp)  ! -(m_Parallelization) ->ista_kngp,iend_kngp
    call m_pwBS_alloc_ngpt_igfp_gr()
    call m_pwBS_calc_length_of_G()         ! -> gr_l
    call m_pwBS_G_trans_functions()   ! -> ngpt_l: Set of G-vectors translated according to symmetry operations
    ggacmp_parallel_rev = ggacmp_parallel
    call m_Parallel_init_mpi_cdfft(nfout,ipriparallel,ggacmp_parallel_rev)
    call m_FFT_setup(inversion_symmetry,paramset) ! paramset == .false.
  call m_pwBS_setup_FFTmapfunctions()
    call m_pwBS_for_each_WF(paramset) ! -> kg1, nbase,iba (when paramset==.false.)
    call m_ES_dealloc(store_prev = .true.)
    call m_ESsd_dealloc_dzajn2()
    call m_ESsd_dealloc_zaj_old()
    call pwbs_mpi_initialization()
    call m_ES_alloc_vlhxc()
    call m_ES_alloc_vlhxcQ()
    call m_ES_alloc_zaj_etc()
    call m_ESsd_alloc_zaj_old()
    call m_ESsd_alloc_dzajn2()
    call m_ESIW_by_randomnumbers(nfout,kv3,1,neg)      ! (rndzaj) -> zaj_l
    if(modnrm == EXECUT ) then
       call m_pwBS_alloc_ylm_l()
       call m_pwBS_sphrp_l()                           ! -> ylm_l
    end if
#ifdef FFTW3
    call m_ES_gen_zaj_from_prev_zaj(nfout)
#endif
    call m_CD_dealloc_chgq(store_prev = .true.)
    call m_CD_alloc_chgq()
#ifdef FFTW3
    call m_CD_gen_chgq_from_prev_chgq(nfout)
#endif

    call m_PP_dealloc()
    call m_NLP_dealloc()
    call PseudoPotential_Construction()
    call m_PP_init_epc()
! ==================================== modified by K. Tagami ================ 11.0
!     call m_PP_gfqwei(nfout)  ! -> modnrm, fqwei, nlmta1, nlmta2
    if ( noncol ) then
       call m_PP_set_fqwei_noncl( nfout )
    else
       call m_PP_gfqwei(nfout)  ! -> modnrm, fqwei, nlmta1, nlmta2
    endif
! ======================================================================== 11.0
    call m_CtrlP_reset_gmax_changed()
  end subroutine change_gmax

  subroutine fft_box_finding_way(outer_or_inner)
    use m_Const_Parameters,  only : GENERAL, SIMPLE_CUBIC, HEXAGONAL &
      , OUTER, INNER
    use m_Crystal_Structure, only : nbztyp
    implicit none
    integer,intent(out)::  outer_or_inner
    if( nbztyp == GENERAL &
         & .or.(nbztyp >= 30 .and. nbztyp <= 32) &
         & .or.(nbztyp >= SIMPLE_CUBIC  .and. nbztyp <= HEXAGONAL)) then
       outer_or_inner = OUTER
    else
       outer_or_inner = INNER
    endif
  end subroutine fft_box_finding_way

  subroutine pwbs_mpi_initialization()
    use m_Files, only : nfout
    use m_Control_Parameters, only : noncol, printable, ipriparallel, neg, nspin, ndim_spinor
    use m_Kpoints, only : kv3
    use m_PlaneWaveBasisSet, only : kg1, iba
    use m_Parallelization, only : m_Parallel_dealloc_mpi_elec, m_Parallel_init_mpi_elec, m_Parallel_init_mpi_iba
    call m_Parallel_dealloc_mpi_elec()
    if ( noncol ) then
        call m_Parallel_init_mpi_elec( nfout, ipriparallel, printable, neg, &
             &                                   kv3, ndim_spinor, kg1 )
     else
        call m_Parallel_init_mpi_elec( nfout, ipriparallel, printable, neg, &
             &                                   kv3, nspin, kg1 )
     endif
     call m_Parallel_init_mpi_iba(nfout,ipriparallel,printable,kv3,iba) ! -> np_g1k, mp_g1k
  end subroutine pwbs_mpi_initialization
#endif

end subroutine Postproc_after_SCF_convergence

