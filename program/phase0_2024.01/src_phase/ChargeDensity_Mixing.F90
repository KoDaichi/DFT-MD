!=======================================================================
!
!  PROGRAM  PHASE/0 2023.01
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
subroutine ChargeDensity_Mixing
! $Id: ChargeDensity_Mixing.F90 615 2020-05-08 13:58:30Z ktagami $
  use m_Const_Parameters,    only : DP,SIMPLE,BROYD1,BROYD2,DFP,PULAY,RMM2P,ON &
       &                          , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, CONTINUATION, SKIP
  use m_Charge_Density,      only : m_CD_check
  use m_CD_mixing,           only : m_CD_simple_mixing, m_CD_prepare_precon &
       &                           ,m_CD_mix_broyden1,m_CD_mix_broyden2,m_CD_mix_DFP &
       &                           ,m_CD_mix_pulay, m_CD_simple_mixing_hsr, m_CD_mixing_write_DEFINITION
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
       &                          , sw_hubbard, sw_mix_occ_matrix, sw_opencore
  use m_Files,               only : nfout
  use m_IterationNumbers,    only : iteration_electronic, iteration_ionic &
       &                          , m_Iter_cmix_reset
  use m_Orbital_Population,  only : m_OP_simple_mixing, m_OP_mix_broyden1 &
       &                          , m_OP_mix_broyden2, m_OP_mix_DFP, m_OP_mix_pulay
  use m_PseudoPotential,     only : flg_paw, modnrm
  use m_PS_opencore,          only : m_PS_print_core_Magmom_on_atom
! ================================= added by K. Tagami =============== 5.0
  use m_Const_Parameters,    only : OFF
  use m_Control_Parameters,  only : sw_update_charge_total, &
       &                            sw_mix_charge_hardpart, occ_matrix_mixfactor
  use m_Charge_Density,      only : m_CD_cp_hsr_to_hsro
  use m_CD_mixing,           only : m_CD_mix_broyden2_with_hsr, m_CD_mix_pulay_with_hsr &
                                  , m_CD_simple_mixing_hard
  use m_IterationNumbers,     only : iteration, iteration_for_cmix
  use m_Orbital_Population,    only :  m_OP_cp_om_to_ommix
! ==================================================================== 5.0

! ================================= added by K. Tagami ================== 11.0&13.0U
  use m_Control_Parameters,   only : noncol, SpinOrbit_mode
  use m_Charge_Density,       only : m_CD_check_noncl
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
! ======================================================================== 11.0&13.0U


! === KT_add ========== 13.0U2-3, 13.0XX
  use m_Control_Parameters,  only : sw_modified_TFW_functional, hardpart_mixfactor, &
       &                            precon_hardpart, sw_wf_mixing, &
       &                            wf_mixing_is_active, &
       &                            recover_wf_after_mixing
  use m_Kpoints,         only : kv3
  use m_ES_WF_mixing,     only : m_ES_WF_precon_hsr, m_ES_WF_kerker_mixing
  use m_Charge_Density,     only : m_CD_hardpart_hsr, m_CD_hardpart_hsr_noncl
! ===================== 13.0U2-3, 13.0XX

! ==== KT_add === 2014/09/01
  use m_Control_Parameters, only : sw_calc_orbital_moment, sw_use_add_proj, &
       &                          Orbital_decomp_mode
  use m_Charge_Density,    only : m_CD_hardpart_hsr_add, m_CD_hardpart_hsr_add_noncl
  use m_CD_Mag_Moment,   only : m_CD_calc_SpinMagMom_method2
  use m_OP_Moment,  only :        m_OP_calc_OrbMagMom_method2
! =============== 2014/09/01

! === KT_add ==== 2014/09/19
  use m_Control_Parameters,  only : sw_calc_ekin_density, ekin_density_type, &
       &                            sw_mix_charge_with_ekindens
  use m_CD_mixing,  only : m_CD_mix_broyden2_intg, m_CD_mix_pulay_intg, &
       &                   m_CD_simple_mixing_intg, m_CD_force_dealloc
  use m_KE_mixing,  only : m_KE_simple_mixing, m_KE_mix_broyden2, m_KE_mix_pulay, &
       &                   m_KE_prepare_precon, &
       &                   m_KE_simple_mixing_00, m_KE_kerker_mixing, &
       &                   m_KE_set_pointer_ekinq_00, m_KE_set_pointer_ekinq
! =============== 2014/09/19

  implicit none
  real(kind=DP) :: rmxt
  real(kind=DP) :: edeltb_per_atom
  integer,save  :: waymix_at_CDmix_previous = 0
  integer       :: waymix_at_CDmix
  logical       :: ini = .false., mixer_changed
  real(kind=DP) :: rmxt_tot, rmxt_hard, rmxt_occmat
  integer,save :: definition_check_in_m_CD_mixing = 0
  interface
     subroutine Renewal_of_OccMat( flag, mode, skip )
        logical flag, skip
        integer mode
        optional :: skip
     end subroutine
  end interface
#ifdef __TIMER_SUB__
    call timer_sta(1101)
#endif
  call m_CtrlP_clear_cdmixing_applied()

  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) then
!      call procedure_charge_fixed;  return
      return
  endif

  if(iprichargemixing>=1) then
     if(definition_check_in_m_CD_mixing == 0) then
        call m_CD_mixing_write_DEFINITION(nfout)
        definition_check_in_m_CD_mixing = 1
     end if
  end if

  edeltb_per_atom = m_TE_what_is_edeltb_now()/natm 
  waymix_at_CDmix = m_CtrlP_waymix_now(iteration_electronic, iteration_ionic &
       &                             , edeltb_per_atom,mixer_changed)
  if(mixer_changed) call m_CD_force_dealloc()
  if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
     if(waymix_at_CDmix_previous /= waymix_at_CDmix) call m_Iter_cmix_reset()
  end if

  rmxt = m_CtrlP_rmx_now(iteration_electronic, iteration_ionic)
  call m_CtrlP_set_mix_parameter()
  if(c_precon) call m_CD_prepare_precon(nfout,rmxt)
  if(iprichargemixing >= 2) &
       & write(6,'(" waymix_at_CDmix = ",i4," rmxt = ",f8.4)') waymix_at_CDmix, rmxt

  call push_cd_mixing_name_applied(waymix_at_CDmix)

! ===================== added by K. Tagami========================= 5.0
  if ( sw_update_charge_total == OFF ) then
     iteration_for_cmix = iteration_for_cmix -1
  endif
  if ( iteration <= 1 .or. (icond==CONTINUATION .and. .not.ini))then
     call m_CD_cp_hsr_to_hsro
     ini = .true.
  endif
!
  rmxt_tot = rmxt
  rmxt_hard = rmxt *hardpart_mixfactor
!
  rmxt_occmat = rmxt *occ_matrix_mixfactor
!
! === KT_add === 2014/09/19
  if ( sw_calc_ekin_density == ON .and. ekin_density_type == 0 ) then
     if ( sw_mix_charge_with_ekindens == ON ) then
        call m_KE_set_pointer_ekinq
     else
        call m_KE_set_pointer_ekinq_00
        if ( c_precon ) call m_KE_prepare_precon(nfout,rmxt)
     endif
  endif
! ============= 2014/09/19

  if ( sw_calc_ekin_density == ON .and. ekin_density_type == 0 ) then
     if ( sw_mix_charge_hardpart == OFF ) then
        call mix_charge_total_intg()
        if ( flg_paw ) call mix_charge_only_hsr          !! only exception
        if ( sw_hubbard == ON ) call mix_charge_only_occmat
     else
        if ( precon_hardpart ) call precon_hardpart_kt
        call mix_charge_total_intg()
     endif
     if ( sw_mix_charge_with_ekindens == OFF ) call mix_kinetic_energy_density()
  else
     if ( sw_update_charge_total == ON ) then
        if ( modnrm == SKIP .or. sw_mix_charge_hardpart == OFF ) then
           call mix_charge_total()
           if ( flg_paw ) call mix_charge_only_hsr          !! only exception
           if ( sw_hubbard == ON ) call mix_charge_only_occmat
        else
           if ( precon_hardpart ) call precon_hardpart_kt
           call mix_charge_total_with_hsr()
        endif
     endif
  endif
! ================================================================= 5.0

! ================================== modified by K. Tagami ========== 11.0
  if ( noncol ) then
     call m_CD_check_noncl(nfout)

! ============================ KT_mod ========== 13.0U
!     if ( iteration_ionic >1 .and. iteration_electronic ==1 ) then
!        call m_CD_set_rad_cov_default
!        call m_CD_set_rad_cov_now
!     endif
! ============================================== 13.0U

     call m_CD_calc_ChgMagMom_in_sphere
     call m_CD_print_ChgMagmom_on_atom(nfout)
     if ( sw_opencore == ON ) call m_PS_print_core_Magmom_on_atom( nfout )

! ==== KT_add === 2014/09/01
     if ( Orbital_decomp_mode == 2 ) then
        call m_CD_calc_SpinMagMom_method2
     else
        if ( .not. flg_paw ) call m_CD_estim_magmom_local(nfout)
        call m_CD_print_magmom_local(nfout)
     endif

     if ( sw_use_add_proj == ON .and. Orbital_decomp_mode == 2 ) then
        call m_CD_hardpart_hsr_add_noncl( nfout,kv3 )
     endif
     if ( sw_calc_orbital_moment == ON ) then
        if ( Orbital_decomp_mode == 2 ) call m_OP_calc_OrbMagMom_method2
     endif
! =============== 2014/09/01

     if ( SpinOrbit_mode /= Neglected ) then
#ifndef WO_OP_MOMENT
        if ( sw_hubbard == ON ) call m_OP_calc_orbmom_from_OCC
#endif
     endif

  else
     call m_CD_check(nfout)
     if ( sw_monitor_atomcharge == ON ) then
        call m_CD_calc_ChgMagMom_in_sphere
        call m_CD_print_ChgMagmom_on_atom(nfout)
        if ( sw_opencore == ON ) call m_PS_print_core_Magmom_on_atom( nfout )
     endif

! ==== KT_add === 2014/08/28
     if ( sw_use_add_proj == ON .and. Orbital_decomp_mode == 2 ) then
        call m_CD_hardpart_hsr_add( nfout,kv3 )
     endif
! =============== 2014/08/28
  endif
! =================================================================== 11.0


  if(tag_solver_of_WF_is_found .and. tag_charge_mixing_is_found) then
     waymix_at_CDmix_previous = waymix_at_CDmix
  end if

#ifdef __TIMER_SUB__
    call timer_end(1101)
#endif

! ========== KT_add ========== 13.0U2
  if ( sw_modified_TFW_functional == ON ) call ThomasFermiWeiz_loop
! ============================ 13.0U2

  call m_TE_wd_total_energy_with_solvers(nfout)

contains

  subroutine mix_charge_total()
    mixing_way: select case(waymix_at_CDmix)
    case (SIMPLE)
       call m_CD_simple_mixing(nfout,rmxt_tot)
    case (BROYD1)
       call m_CD_mix_broyden1(rmxt_tot)
    case (BROYD2)
       call m_CD_mix_broyden2(nfout,rmxt_tot)
    case (DFP)
       call m_CD_mix_DFP(rmxt_tot)
    case (PULAY)
       call m_CD_mix_pulay(nfout,rmxt_tot)
    case default
       call phase_error_with_msg(nfout,' ! waymix is invalid',__LINE__,__FILE__)
    end select mixing_way

  end subroutine mix_charge_total

  subroutine mix_charge_only_hsr
     call m_CD_simple_mixing_hsr(nfout,rmxt)
  end subroutine mix_charge_only_hsr

  subroutine mix_charge_only_occmat
    mixing_way: select case(waymix_at_CDmix)
    case (SIMPLE)
        call m_OP_simple_mixing(nfout,rmxt_occmat)
    case (BROYD1)
        call m_OP_mix_broyden1(rmxt_occmat)
    case (BROYD2)
        call m_OP_mix_broyden2(rmxt_occmat)
    case (DFP)
        call m_OP_mix_DFP(rmxt_occmat)
    case (PULAY)
        call m_OP_mix_Pulay(rmxt_occmat)
    case default
       call phase_error_with_msg(nfout,' ! waymix is invalid',__LINE__,__FILE__)
    end select mixing_way

  end subroutine mix_charge_only_occmat

  logical function fix_occ_matrix()
    use m_Control_Parameters, only : occ_matrix_fix_period
    use m_IterationNumbers, only : iteration_electronic
    fix_occ_matrix = iteration_electronic <= occ_matrix_fix_period 
    return
  end function fix_occ_matrix

  subroutine mix_charge_total_with_hsr()
    mixing_way: select case(waymix_at_CDmix)

    case (SIMPLE)
       call m_CD_simple_mixing( nfout,rmxt_tot )
       call m_CD_simple_mixing_hard( nfout, rmxt_hard )
       if ( sw_hubbard == ON .and.sw_mix_occ_matrix==OFF ) then
          call Renewal_of_OccMat( .false., ON, fix_occ_matrix() )           ! hsr --> om 
       endif
       if (sw_hubbard==ON) call m_OP_cp_om_to_ommix( nfout, rmxt_hard )       ! om --> ommix

    case (BROYD1)
!       call m_CD_mix_broyden1_with_hsr(rmxt_tot)
       write(*,*) 'Not supported '
       call phase_error_with_msg(nfout,' broyd1 unsupported',__LINE__,__FILE__)

    case (BROYD2)
       call m_CD_mix_broyden2_with_hsr(nfout,rmxt_tot,sw_mix_occ_matrix==ON)
       if ( sw_hubbard == ON .and. sw_mix_occ_matrix==OFF ) then
          call Renewal_of_OccMat(.false., ON, fix_occ_matrix() )           ! hsr --> om 
       endif
       if(sw_hubbard==ON) call m_OP_cp_om_to_ommix( nfout, rmxt_hard )      ! om --> ommix
    case (DFP)
!       call m_CD_mix_DFP(rmxt_tot)
       !write(*,*) 'Not supported '
!       stop
       call phase_error_with_msg(nfout,' dfp unsupported',__LINE__,__FILE__)
    case (PULAY)
       call m_CD_mix_pulay_with_hsr(nfout,rmxt_tot,sw_mix_occ_matrix==ON)
       if ( sw_hubbard == ON .and. sw_mix_occ_matrix==OFF ) then
          call Renewal_of_OccMat(.false., ON, fix_occ_matrix() )           ! hsr --> om 
       endif
       if(sw_hubbard==ON) call m_OP_cp_om_to_ommix( nfout, rmxt_hard )      ! om --> ommix
    case default
       call phase_error_with_msg(nfout,' ! waymix is invalid',__LINE__,__FILE__)
    end select mixing_way
  end subroutine mix_charge_total_with_hsr

! === KT_add == 2014/09/19
  subroutine mix_charge_total_intg()
    mixing_way: select case(waymix_at_CDmix)

    case (SIMPLE)
       call m_CD_simple_mixing_intg( nfout,rmxt_tot )
       if ( sw_hubbard == ON ) then
          call Renewal_of_OccMat( .false., ON )           ! hsr --> om
          call m_OP_cp_om_to_ommix( nfout, rmxt_hard )       ! om --> ommix
       endif

    case (BROYD2)
       call m_CD_mix_broyden2_intg(nfout,rmxt_tot,sw_mix_occ_matrix==ON)
       if ( sw_hubbard == ON ) then
          call Renewal_of_OccMat(.false., ON )           ! hsr --> om
          call m_OP_cp_om_to_ommix( nfout, rmxt_hard )       ! om --> ommix
       endif
    case (PULAY)
       call m_CD_mix_pulay_intg(nfout,rmxt_tot,sw_mix_occ_matrix==ON)
       if ( sw_hubbard == ON ) then
          call Renewal_of_OccMat(.false., ON )           ! hsr --> om
          call m_OP_cp_om_to_ommix( nfout, rmxt_hard )       ! om --> ommix
       endif

    case default
       call phase_error_with_msg(nfout,' ! waymix is invalid',__LINE__,__FILE__)
    end select mixing_way

  end subroutine mix_charge_total_intg

  subroutine mix_kinetic_energy_density()
    mixing_way: select case(waymix_at_CDmix)

    case (SIMPLE)
       call m_KE_simple_mixing(nfout,rmxt_tot)
    case (BROYD2)
       call m_KE_mix_broyden2(nfout,rmxt_tot)
    case (PULAY)
       call m_KE_mix_pulay(nfout,rmxt_tot)
    case default
       call phase_error_with_msg(nfout,' ! waymix is invalid',__LINE__,__FILE__)
    end select mixing_way

  end subroutine mix_kinetic_energy_density

  subroutine procedure_charge_fixed
    if ( sw_calc_ekin_density == ON ) then
       call m_KE_simple_mixing_00( rmxt_tot )
!     call m_KE_Kerker_mixing( rmxt_tot )
    endif
  end subroutine procedure_charge_fixed
! ======= 2014/09/19

! === KT_add ======= 13.0U3
  subroutine precon_hardpart_kt
    real(kind=DP) :: threshold = 2.5D-2

    if ( sw_wf_mixing == OFF ) return
    if ( .not. wf_mixing_is_active ) return

#if 0
    if ( abs(edeltb_per_atom) > threshold ) return
#endif

    if ( sw_wf_mixing == ON ) then
!       call m_ES_WF_precon_hsr
       call m_ES_WF_kerker_mixing( .true. )
    endif

  end subroutine precon_hardpart_kt
! ================== 13.0U3

end subroutine ChargeDensity_Mixing

subroutine push_cd_mixing_name_applied(waymix)
  use m_Const_Parameters,   only :  SIMPLE , BROYD1, BROYD2, DFP, PULAY
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
  case (BROYD1)
     name = "BROYD1"
  case (BROYD2)
     name = "BROYD2"
  case (DFP)
     name = "DFP"
  case (PULAY)
     name = "PULAY"
  case default
     name = "notDefined"
  end select solver
  len_char = len_trim(name)
  cdmixname = ""
  cdmixname(1:min(len_char,len_cdmixingname)) = name(1:min(len_char,len_cdmixingname))
  call m_CtrlP_push_CDMixingNameApplied(cdmixname,len_cdmixingname)
end subroutine push_cd_mixing_name_applied
