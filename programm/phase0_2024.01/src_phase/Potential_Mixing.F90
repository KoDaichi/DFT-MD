!=======================================================================
subroutine Potential_Mixing
! $Id: Potential_Mixing.F90 593 2019-06-20 03:47:31Z jkoga $

  use m_Const_Parameters,    only : DP,SIMPLE,BROYD1,BROYD2,DFP,PULAY,RMM2P,ON &
       &                          , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, CONTINUATION, SKIP

  use m_Charge_Density,      only : m_CD_check

  use m_Potential_mixing,     only : m_Pot_simple_mixing, m_Pot_prepare_precon, &
       &                             m_Pot_mix_broyden2, m_Pot_mix_pulay, &
       &                             m_Pot_simple_mixing_vhsr, &
       &                             m_Pot_set_array_vnonlocal, &
       &                             m_Pot_cp_vnonlocal_to_old
  use m_Ionic_System,        only : natm
  use m_Total_Energy,        only : m_TE_what_is_edeltb_now, m_TE_wd_total_energy_with_solvers
  use m_Control_Parameters,  only : c_precon,waymix,intzaj,icond &
       &                          , tag_solver_of_WF_is_found &
       &                          , tag_charge_mixing_is_found &
       &                          , iprichargemixing, nspin, af &
       &                          , m_CtrlP_rmx_now &
       &                          , m_CtrlP_solver_for_WFs_now &
       &                          , m_CtrlP_set_rmx &
       &                          , m_CtrlP_waymix_now &
       &                          , m_CtrlP_set_mix_parameter &
       &                          , m_CtrlP_clear_cdmixing_applied &
       &                          , sw_hubbard, sw_modified_TFW_functional
  use m_Files,               only : nfout
  use m_IterationNumbers,    only : iteration_electronic, iteration_ionic &
       &                          , m_Iter_cmix_reset
  use m_PseudoPotential,     only : flg_paw, modnrm

  use m_Const_Parameters,    only : OFF
  use m_Control_Parameters,  only : sw_mix_charge_hardpart, ekmode
  use m_Potential_mixing,   only : m_Pot_mix_broyden2_with_vhsr, &
       &                             m_Pot_mix_pulay_with_vhsr,  &
       &                             m_Pot_simple_mixing_hard, &
       &                             m_Pot_add_vnonlocal_to_vlhxcQ, &
       &                             vnonlocal_i, vnonlocal_r
  use m_IterationNumbers,     only : iteration, iteration_for_cmix, iteration_electronic
  use m_Orbital_Population,    only :  m_OP_cp_om_to_ommix

  use m_Control_Parameters,   only : noncol, SpinOrbit_mode
  use m_Charge_Density,       only : m_CD_check_noncl, m_CD_check
  use m_CD_Mag_Moment,        only : m_CD_calc_ChgMagMom_in_sphere, &
       &                             m_CD_print_ChgMagmom_on_atom, &
       &                             m_CD_estim_magmom_local, &
       &                             m_CD_print_magmom_local, &
       &                             m_CD_set_rad_cov_now, &
       &                             m_CD_set_rad_cov_default, &
       &                             sw_monitor_atomcharge
  use m_Const_Parameters,      only : Neglected
#ifndef WO_OP_MOMENT
  use m_OP_Moment,             only : m_OP_calc_orbmom_from_OCC
#endif

  use m_Electronic_Structure,   only : m_ES_cp_vlhxcQ_to_old, m_ES_cp_vlhxc_to_old, &
       &                               vlhxcQ, dhub
  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_integrate_VlhxcQlm
  use m_ES_NonCollinear,   only : m_ES_update_Mat_dion0_noncl, &
       &                          m_ES_set_Mat_dion_scr_noncl

  implicit none

  real(kind=DP) :: rmxt
  real(kind=DP) :: edeltb_per_atom
  integer,save  :: waymix_at_Potmix_previous = 0
  integer       :: waymix_at_Potmix
  logical       :: ini = .false.
  real(kind=DP) :: rmxt_tot, rmxt_hard
#ifdef __TIMER_SUB__
  call timer_sta(1101)
#endif

  call m_CtrlP_clear_cdmixing_applied()

  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) return

! ------------------------------------------------------------
! --- copy the last potential to old ones, and calc the latest potential 
!

#if 0
  if ( iteration <= 1 .or. (icond==CONTINUATION .and. .not.ini))then
     call m_Pot_set_array_vnonlocal          ! vnonlocal = dhub + ....
     ini = .true.
  endif

  if ( sw_modified_TFW_functional == OFF ) then
!!!     call m_ES_cp_vlhxc_to_old
     if ( sw_mix_charge_hardpart == ON ) call m_Pot_cp_vnonlocal_to_old
  endif

#endif

  call m_Pot_set_array_vnonlocal
! ------------------------------------------------------------

  edeltb_per_atom = m_TE_what_is_edeltb_now()/natm 
  waymix_at_Potmix = m_CtrlP_waymix_now(iteration_electronic, iteration_ionic &
       &                             , edeltb_per_atom)
  if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
     if(waymix_at_Potmix_previous /= waymix_at_Potmix) call m_Iter_cmix_reset()
  end if

  rmxt = m_CtrlP_rmx_now(iteration_electronic, iteration_ionic)
  call m_CtrlP_set_mix_parameter()

  if(c_precon) call m_Pot_prepare_precon(nfout,rmxt)
  if(iprichargemixing >= 2) &
       & write(6,'(" waymix_at_Potmix = ",i4," rmxt = ",f8.4)') waymix_at_Potmix, rmxt


  call push_pot_mixing_name_applied(waymix_at_Potmix)
! ------------------------------------------------------------
! --- main part of potential mixnig 
! 
  rmxt_tot = rmxt;       rmxt_hard = rmxt

#if 0
  if ( sw_mix_charge_hardpart == ON ) then
     call mix_potential_with_hardpart()      ! does not work properly.
  else
     call mix_potential()
  endif

#else
  call mix_potential()
  if ( sw_mix_charge_hardpart == ON ) call mix_potential_only_hardpart         
#endif
! ------------------------------------------------------------

! ------------------------------------------------------------
! --- update vlhxcQ and related variables
!
  call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ
  call m_Pot_add_vnonlocal_to_vlhxcQ

  if ( noncol ) then
     if ( flg_paw ) then
        if ( ekmode==ON .or. iteration_electronic >= 1 ) then
           call m_ES_update_Mat_dion0_noncl
        endif
     endif
     call m_ES_set_Mat_dion_scr_noncl( VlhxcQ, vnonlocal_i )
  endif

! ================================== modified by K. Tagami ========== 11.0
  if ( noncol ) then
     call m_CD_check_noncl(nfout)

     if ( .not. flg_paw ) call m_CD_estim_magmom_local(nfout)
     call m_CD_print_magmom_local(nfout)

     call m_CD_calc_ChgMagMom_in_sphere
     call m_CD_print_ChgMagmom_on_atom(nfout)

     if ( SpinOrbit_mode /= Neglected ) then
#ifndef WO_OP_MOMENT
        if ( sw_hubbard == ON ) call m_OP_calc_orbmom_from_OCC
#endif
     endif

  else
     call m_CD_check(nfout)
! ============================ KT_add ========== 13.0U
     if ( sw_monitor_atomcharge == ON ) then
        call m_CD_calc_ChgMagMom_in_sphere
        call m_CD_print_ChgMagmom_on_atom(nfout)
     endif
! ============================================== 13.0U

  endif
! =================================================================== 11.0

  if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
     waymix_at_Potmix_previous = waymix_at_Potmix
  end if

#ifdef __TIMER_SUB__
    call timer_end(1101)
#endif

! ========== KT_add ========== 13.0U2
  if ( sw_modified_TFW_functional == ON ) call ThomasFermiWeiz_loop
! ============================ 13.0U2
  
  call m_TE_wd_total_energy_with_solvers(nfout)

contains

  subroutine mix_potential()
    mixing_way: select case(waymix_at_Potmix)
    case (SIMPLE)
       call m_Pot_simple_mixing(nfout,rmxt_tot)
    case (BROYD2)
       call m_Pot_mix_broyden2(rmxt_tot)
    case (PULAY)
       call m_Pot_mix_pulay(nfout,rmxt_tot)
    case default
       !stop ' ! waymix is invalid'
       call phase_error_with_msg(nfout,' waymix is invalid',__LINE__,__FILE__)
    end select mixing_way

  end subroutine mix_potential

  subroutine mix_potential_only_hardpart
     call m_Pot_simple_mixing_vhsr(nfout,rmxt)
   end subroutine mix_potential_only_hardpart

  subroutine mix_potential_with_hardpart()

    mixing_way: select case(waymix_at_Potmix)

    case (SIMPLE)
       call m_Pot_simple_mixing( nfout,rmxt_tot )
       call m_Pot_simple_mixing_hard( nfout, rmxt_hard )

    case (BROYD2)
       call m_Pot_mix_broyden2_with_vhsr(rmxt_tot)

    case (PULAY)
       call m_Pot_mix_pulay_with_vhsr(nfout,rmxt_tot)

    case default
       !stop ' ! waymix is invalid'
       call phase_error_with_msg(nfout,' waymix is invalid',__LINE__,__FILE__)
    end select mixing_way
  end subroutine mix_potential_with_hardpart

end subroutine Potential_Mixing

subroutine push_pot_mixing_name_applied(waymix)
  use m_Const_Parameters,   only :  SIMPLE , BROYD2, PULAY
  use m_Control_Parameters, only :  m_CtrlP_push_CDMixingNameApplied, len_cdmixingname
  implicit none
  integer, intent(in) :: waymix
  character(len=len_cdmixingname) :: cdmixname
  character*20 :: name
  integer :: len_char
  name = ""
  solver: select case(waymix)
  case (SIMPLE)
     name = "SIMPLE"
  case (BROYD2)
     name = "BROYD2"
  case (PULAY)
     name = "PULAY"
  case default
     name = "notDefined"
  end select solver
  len_char = len_trim(name)
  cdmixname = ""
  cdmixname(1:min(len_char,len_cdmixingname)) = name(1:min(len_char,len_cdmixingname))
  call m_CtrlP_push_CDMixingNameApplied(cdmixname,len_cdmixingname)
end subroutine push_pot_mixing_name_applied
