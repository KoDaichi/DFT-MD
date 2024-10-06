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
  use m_Charge_Density,      only : chgq_l 
  use m_Const_Parameters,    only : Valence_plus_PC_Charge, VXC_AND_EXC &
       &                          , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, ON, OFF
  use m_Control_Parameters,  only : nspin, af, kimg, icond, ekmode
  use m_Electronic_Structure,only : vlhxc_l,  vlhxcQ
  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_integrate_VlhxcQlm_3D
  use m_ES_LHXC,             only : m_ESlhxc_potential_3D
  use m_ES_Mag_Constraint,    only : m_ES_add_MagConstraintPot_chgql, &
       &                             m_ES_add_MagConstraintPot_hsr
  use m_Files,               only : nfout
  use m_PAW_Hartree,         only : m_PAWH_get_dion_hartree
  use m_PAW_XC_Potential,    only : m_PAW_XC_cal_potential &
                                  , m_PAW_XC_get_dion_vxc &
                                  , m_PAW_XC_allocation, m_PAW_XC_deallocation
! === DEBUG to make 3D_Parallel. by tkato 2011/07/18 ===========================
!                                 , m_PAW_XC_cal_potential_sym &
! ==============================================================================
  use m_PseudoPotential,     only : m_PP_get_dion_paw,flg_symmtry, flg_paw, modnrm
  use m_XC_Potential,        only : vxc_l         &
 &                                , m_XC_cal_potential_3D
! === DEBUG by tkato 2011/10/05 ================================================
  use m_PAW_XC_Potential,    only : m_PAW_XC_cal_potential_sphex2
  use m_PAW_ChargeDensity,   only : calcGaussLegendreIntegration &
                                   ,calcSphericalHarmonicsExpansion
!===============================================================================
  use m_Crystal_Structure,   only : sw_magnetic_constraint

  use m_Control_Parameters,  only : sw_opencore, nspin
  use m_PS_Opencore,  only : m_PS_set_mag_opencore_pol

  use m_Control_Parameters,  only : sw_precalculate_phase_vnonlocal

  use m_ES_nonlocal,         only : m_ES_AtaulmnaG

  use m_XC_Potential_2D,        only : m_XC_cal_potential

  implicit none

  if ( sw_opencore == ON .and. nspin == 2 ) then
     call m_PS_set_mag_opencore_pol
  endif

#if 1
  call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC)
#else
  call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC)
#endif

! === DEBUG by tkato 2011/10/05 ================================================
! if(flg_paw) then
!===============================================================================

! === DEBUG by tkato 2011/10/05 ================================================
!    call m_PAW_XC_cal_potential(nfout,VXC_AND_EXC,flg_symmtry)

  call m_PAW_XC_allocation(nfout)
  if(calcSphericalHarmonicsExpansion) then
     call m_PAW_XC_cal_potential_sphex2(nfout,VXC_AND_EXC)
  end if

  if(calcGaussLegendreIntegration)then
     call m_PAW_XC_cal_potential(nfout,VXC_AND_EXC,flg_symmtry)
  end if
#ifdef _MEMORY_CONSUMPTION_CHECK_
     call memsize(nfout)
#endif
! ==============================================================================

! === DEBUG by tkato 2011/10/05 ================================================
  if(flg_paw) then
!===============================================================================
     call m_PAW_XC_get_dion_vxc(nfout)
     !call m_PAW_XC_get_dion_vxc_dbg(nfout)

  endif

  call m_PAW_XC_deallocation(nfout)

  call m_ESlhxc_potential_3D(nfout,chgq_l,vxc_l) ! (stlhxc) -> vlhxc_l

! ================================ added by K. Tagami ============== 11.0&13.0U
!
  if ( sw_magnetic_constraint == ON .and. ekmode == OFF ) then
     call m_ES_add_MagConstraintPot_chgql
  endif
!
! ================================================================= 11.0&13.0U

  call m_ESiVQ_integrate_VlhxcQlm_3D(nfout) ! (lclchh) -> vlhxcQ
  if(sw_precalculate_phase_vnonlocal==ON .and. modnrm==ON) then
    call m_ES_AtaulmnaG()
  endif

  if(flg_paw)then
     call m_PAWH_get_dion_hartree(nfout)
     call m_PP_get_dion_paw(nfout)
  endif

! ================================ added by K. Tagami ============== 11.0&13.0U
  if ( sw_magnetic_constraint == ON .and. ekmode == OFF ) then
     call m_ES_add_MagConstraintPot_hsr
  endif
! ================================================================= 11.0&13.0U
end subroutine Renewal_of_Potential

