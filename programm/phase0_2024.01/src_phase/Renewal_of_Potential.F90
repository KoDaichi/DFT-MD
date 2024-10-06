!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 604 $)
!
!  SUBROUINE: Renewal_of_Potential
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
subroutine Renewal_of_Potential()
! $Id: Renewal_of_Potential.F90 604 2020-04-09 05:52:05Z jkoga $
  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_integrate_VlhxcQlm
  use m_ES_LHXC,             only : m_ESlhxc_potential
  use m_Charge_Density,      only : chgq_l
  use m_XC_Potential,        only : m_XC_cal_potential, vxc_l
  use m_Control_Parameters,  only : icond
  use m_Const_Parameters,    only : Valence_plus_PC_Charge, VXC_AND_EXC &
       &                          , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                          , ON
  use m_Files,               only : nfout
  use m_PAW_Hartree,         only : m_PAWH_get_dion_hartree
  use m_PAW_XC_Potential,    only : m_PAW_XC_cal_potential &
                                    ,m_PAW_XC_cal_potential_sphex2 &
                                    ,m_PAW_XC_get_dion_vxc
  use m_PAW_ChargeDensity,   only : calcGaussLegendreIntegration &
                                    ,calcSphericalHarmonicsExpansion
  use m_PseudoPotential,     only : m_PP_get_dion_paw,flg_symmtry, flg_paw, modnrm

! ================================= added by K. Tagami =============== 11.0&13.0U
  use m_Control_Parameters,   only : noncol
  use m_Electronic_Structure,  only : vlhxcQ,vlhxc_l
  use m_ES_nonlocal,           only : m_ES_AtaulmnaG
  use m_Parallelization,       only : ista_kngp,iend_kngp
  use m_ES_NonCollinear,      only : m_ES_set_Mat_dion_scr_noncl, &
       &                             m_ES_update_Mat_dion0_noncl
  use m_ES_LHXC,              only : m_ESlhxc_potential_noncl
!
  use m_ES_Mag_Constraint,    only : m_ES_add_MagConstraintPot_chgql, &
       &                             m_ES_add_MagConstraintPot_hsr
  use m_PAW_XC_Potential,    only : m_PAW_XC_get_dion_vxc_noncl, &
       &                             m_PAW_XC_get_dion_vxc_noncl2, &
       &                             m_PAW_XC_get_dion_vxc_noncl3
  use m_PseudoPotential,     only : m_PP_get_dion_paw_noncl
  use m_IterationNumbers,    only : iteration_electronic
  use m_Crystal_Structure,   only : level_of_projection_paw_charge,&
       &                            sw_magnetic_constraint
! ==================================================================== 11.0&13.0U

! ================================= added by K. Tagami =============== 11.0
  use m_Control_Parameters,   only : SpinOrbit_mode, ekmode, sw_hubbard
  use m_SpinOrbit_Potential,  only : m_SO_set_Dsoc_potential2
  use m_Const_Parameters,    only : ByPawPot, ON, OFF
  use m_SpinOrbit_RadInt,    only : m_SO_calc_SOC_strength_pawpot
! ==================================================================== 11.0

  use m_Control_Parameters,  only : sw_opencore, nspin
  use m_PS_Opencore,  only : m_PS_set_mag_opencore_pol

  use m_Control_Parameters,  only : sw_precalculate_phase_vnonlocal

  implicit none

!!$  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) return

  if ( sw_opencore == ON .and. nspin == 2 ) then
     call m_PS_set_mag_opencore_pol
  endif

  call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC)
                                              ! (xcfft) -> vxc_l
! ====================== added by K. Tagami ================= 11.0
  if ( noncol ) then
     if ( ekmode==OFF .and. iteration_electronic < 1 ) goto 100
  endif
! ========================================================== 11.0

  if(calcSphericalHarmonicsExpansion) then
    call m_PAW_XC_cal_potential_sphex2(nfout,VXC_AND_EXC)
#ifdef _MEMORY_CONSUMPTION_CHECK_
     call memsize(nfout)
#endif
!    call m_PAW_XC_cal_potential_sphex3(nfout,VXC_AND_EXC)
  end if
  if(calcGaussLegendreIntegration)then
!!!$  if(flg_paw) then
     call m_PAW_XC_cal_potential(nfout,VXC_AND_EXC,flg_symmtry)
!!!$  if(flg_symmtry) then
!!!$    call m_PAW_XC_cal_potential_sym(nfout,VXC_AND_EXC)
!!!$  else
!!!$    call m_PAW_XC_cal_potential(nfout,VXC_AND_EXC)
!!!$  end if
  end if

! ============================= modified by K. Tagami ============= 11.0
!  if(flg_paw) call m_PAW_XC_get_dion_vxc(nfout)
  if (flg_paw) then

     if ( noncol ) then
        select case(level_of_projection_paw_charge)
        case (1)
           call m_PAW_XC_get_dion_vxc_noncl(nfout)
        case (2)
           call m_PAW_XC_get_dion_vxc_noncl2(nfout)
        case (3)
           call m_PAW_XC_get_dion_vxc_noncl3(nfout)
        end select

     else
        call m_PAW_XC_get_dion_vxc(nfout)
     endif
  endif
! ================================================================= 11.0

! =================== added by K. Tagami ================== 11.0
100 continue
! ========================================================= 11.0

#ifdef _MEMORY_CONSUMPTION_CHECK_
     call memsize(nfout)
#endif

!! if(flg_paw) call m_PAW_XC_get_dion_vxc(nfout)
!! call m_PAW_XC_get_dion_vxc_dbg(nfout)

! ================================ modified by K. Tagami ============== 11.0
!!  call m_ESlhxc_potential(nfout,chgq_l,vxc_l) ! (stlhxc) -> vlhxc_l

  if ( noncol ) then
     call m_ESlhxc_potential_noncl(nfout,chgq_l,vxc_l)
  else
     call m_ESlhxc_potential(nfout,chgq_l,vxc_l) ! (stlhxc) -> vlhxc_l
  endif
! ====================================================================== 11.0

! ================================ added by K. Tagami ============== 11.0&13.0U
!
  if ( sw_magnetic_constraint == ON .and. ekmode == OFF ) then
     call m_ES_add_MagConstraintPot_chgql
  endif
!
! ================================================================= 11.0&13.0U

  call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ

  if(sw_precalculate_phase_vnonlocal==ON .and. modnrm==ON) then
    call m_ES_AtaulmnaG()
  endif

! ======================== modified by K. Tagami =============== 11.0
!  if(flg_paw)then
!     call m_PAWH_get_dion_hartree(nfout)
!     call m_PP_get_dion_paw(nfout)
!  endif

  if (flg_paw) then
     if ( noncol ) then
        if ( iteration_electronic >= 1 .or. ekmode == ON ) then
           call m_PAWH_get_dion_hartree(nfout)
           call m_PP_get_dion_paw_noncl(nfout)

           if ( SpinOrbit_Mode == ByPawPot ) then
              call m_SO_calc_SOC_strength_pawpot
              call m_SO_set_Dsoc_potential2
           endif
        endif

     else
        call m_PAWH_get_dion_hartree(nfout)
        call m_PP_get_dion_paw(nfout)
     endif
  endif
! ============================================================ 11.0

! ================================ added by K. Tagami ============== 11.0&13.0U
  if ( sw_magnetic_constraint == ON .and. ekmode == OFF ) then
     call m_ES_add_MagConstraintPot_hsr
  endif
! ================================================================= 11.0&13.0U

! ================================= added by K. Tagami ================ 11.0
  if ( noncol ) then
     if ( flg_paw ) then
        if ( ekmode==ON .or. iteration_electronic >= 1 ) then 
           call m_ES_update_Mat_dion0_noncl
        endif
     endif
     if ( sw_hubbard == OFF ) call m_ES_set_Mat_dion_scr_noncl( VlhxcQ )
  endif
! ===================================================================== 11.0

end subroutine Renewal_of_Potential



