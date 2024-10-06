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
!$$#ifndef PARA3D
! $Id: Stress.F90 577 2017-12-15 05:26:08Z jkoga $
  use m_Stress,           only : hsrd &
       &                       , m_Stress_alloc          , m_Stress_dealloc &
       &                       , m_Stress_alloc_grinv_g  , m_Stress_dealloc_grinv_g &
       &                       , m_Stress_ke             , m_Stress_hsr_diff &
       &                       , m_Stress_lclchg_diff    , m_Stress_xc &
       &                       , m_Stress_loc            , m_Stress_h &
       &                       , m_Stress_nl             , m_Stress_ew &
       &                       , m_Stress_or             , m_Stress_sum &
       &                       , m_Stress_check_g12
  use m_XC_Potential,     only : m_XC_alloc_s_gga,m_XC_dealloc_s_gga ,m_XC_cal_potential
  use m_Charge_Density,   only : chgsoft, chgq_l, hsr
  use m_Kpoints,          only : kv3, vkxyz
  use m_Files,            only : nfout
  use m_Control_Parameters,only: istress
  use m_Const_Parameters, only : Valence_plus_PC_Charge,STRESS_, CRDTYP

! ================================ Added by K. Tagami ============ p1100
  use m_Stress,           only : m_Stress_Hub
! ================================================================ p1100

! ================= added by K. Tagami ==================== 11.0
  use m_Const_Parameters,       only : ON
  use m_Control_Parameters,    only : sw_hubbard, noncol
  use m_Stress,                only : m_Stress_ke_noncl, &
       &                              m_Stress_hsr_diff_noncl, &
       &                              m_Stress_lclchg_diff_noncl, &
       &                              m_Stress_xc_noncl, &
       &                              m_Stress_nl_nonclA, &
       &                              m_Stress_or_nonclA, &
       &                              m_Stress_nl_nonclB, &
       &                              m_Stress_or_nonclB, &
       &                              m_Stress_Hub_nonclA, &
       &                              m_Stress_Hub_nonclB
! ========================================================= 11.0

! =================================== KT_add ============== 13.0B
  use m_Control_Parameters,    only : sw_vdw_correction, xctype
#ifndef DISABLE_VDWDF
  use m_Stress,                only : m_Stress_vdw,m_Stress_vdwdf
#else
  use m_Stress,                only : m_Stress_vdw
#endif
 
! ========================================================= 13.0B

  implicit none

  if (istress == 0) return

  if ( noncol ) then
     call case_noncollinear
  else
     call case_collinear
  endif


! ============================ added by K. Tagami ============= 11.0
contains

  subroutine case_collinear
    call m_Stress_alloc()
    call m_XC_alloc_s_gga
    call m_Stress_ke(nfout,kv3,CRDTYP,vkxyz)
    call m_Stress_hsr_diff(nfout,kv3)
!!$  call m_M_Stress_xcfft_diff(nfout)

    call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,STRESS_ &
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

  end subroutine case_collinear

  subroutine case_noncollinear
    call m_Stress_alloc()
    call m_XC_alloc_s_gga

    call m_Stress_ke_noncl( nfout, kv3, CRDTYP, vkxyz )
    call m_Stress_hsr_diff_noncl( nfout, kv3 )

    call m_XC_cal_potential( nfout, Valence_plus_PC_Charge, chgq_l, STRESS_,  &
         &                   chgsoft, hsr, hsrd )
    call m_Stress_check_g12(nfout)

    call m_Stress_alloc_grinv_g()

    call m_Stress_lclchg_diff_noncl(nfout)
    call m_Stress_xc_noncl(nfout)

    call m_Stress_loc(nfout)
    call m_Stress_h(nfout)

    call m_Stress_dealloc_grinv_g()

!!    call m_Stress_nl_nonclA(nfout)
    call m_Stress_nl_nonclB(nfout)

    call m_Stress_ew(nfout)

!!    call m_Stress_or_nonclA(nfout)
    call m_Stress_or_nonclB(nfout)

    if ( sw_hubbard == ON ) call m_Stress_Hub_nonclB(nfout)

! =================================== KT_add ============== 13.0B
    if ( sw_vdw_correction == ON ) call m_Stress_vdw(nfout)
! ========================================================= 13.0B

    call m_Stress_sum(nfout)
    call m_Stress_dealloc()
    call m_XC_dealloc_s_gga()

  end subroutine case_noncollinear

!$$#endif
end subroutine Stress
