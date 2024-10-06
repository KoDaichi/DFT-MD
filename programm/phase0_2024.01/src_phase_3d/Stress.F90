!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 577 $)
!
!  SUBROUINE: Stress
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
subroutine Stress
  use m_Stress,           only : hsrd &
         &                     , m_Stress_alloc,m_Stress_ke, m_Stress_dealloc, m_Stress_hsr_diff &
         &                     , m_Stress_alloc_grinv_g, m_Stress_lclchg_diff, m_Stress_dealloc_grinv_g &
         &                     , m_Stress_xc, m_Stress_loc, m_Stress_h, m_Stress_nl, m_Stress_ew &
         &                     , m_Stress_or, m_Stress_sum, m_Stress_Hub, m_Stress_vdw &
#ifndef DISABLE_VDWDF
         &                     , m_Stress_check_g12, m_Stress_vdwdf
#else
         &                     , m_Stress_check_g12
#endif
  use m_Const_Parameters, only : Valence_plus_PC_Charge,STRESS_, CRDTYP, ON
  use m_Kpoints,          only : kv3, vkxyz
  use m_Files,            only : nfout
  use m_Control_Parameters, only: istress, sw_vdw_correction, xctype, sw_hubbard
  use m_XC_Potential,     only : m_XC_alloc_s_gga, m_XC_dealloc_s_gga, m_XC_cal_potential_3D
  use m_Charge_Density,   only : chgsoft, chgq_l, hsr

  implicit none

  if (istress == 0) return

  call m_Stress_alloc()
  call m_XC_alloc_s_gga()
  call m_Stress_ke(nfout,kv3,CRDTYP,vkxyz)
  call m_Stress_hsr_diff(nfout,kv3)
  call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge,chgq_l,STRESS_ &
       &   ,chgsoft,hsr,hsrd)
  call m_Stress_check_g12(nfout)
  call m_Stress_alloc_grinv_g()
  call m_Stress_lclchg_diff(nfout)
  call m_Stress_xc(nfout)
  call m_Stress_loc(nfout)
  call m_Stress_h(nfout)
  call m_Stress_dealloc_grinv_g()
  call m_Stress_nl(nfout)
  call m_Stress_ew(nfout)
  call m_Stress_or(nfout)

! =============================== Added by K. Tagami ============ p1100
  if ( sw_hubbard == ON ) call m_Stress_Hub(nfout)
! =============================================================== p1100

! =================================== KT_add ============== 13.0B
  if ( sw_vdw_correction == ON ) call m_Stress_vdw(nfout)
! ========================================================= 13.0B

#ifndef DISABLE_VDWDF
  if ( xctype == "vdwdf" ) call m_Stress_vdwdf(nfout)
#endif

  call m_Stress_sum(nfout)
  call m_Stress_dealloc()
  call m_XC_dealloc_s_gga()
end subroutine Stress
