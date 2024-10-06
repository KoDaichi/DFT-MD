!=======================================================================
!
!  SOFTWARE NAME : PHASE (ver. 7.01)
!
!  SUBROUINE: Renewal_of_Hubbard_Potential
!
!  AUTHOR(S): T. Yamamoto   May/08/2005
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
! ==========================
!     patch 0.1 by K. Tagami @adv    2009/04/21
!
!     patch 0.1 : commented out an if sentence for band calculation with DFT+U
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
subroutine Renewal_of_Hubbard_Potential()
! $Id: Renewal_of_Hubbard_Potential.f90 376 2014-06-17 07:48:31Z jkoga $
  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_add_dhub_to_vlhxcQ
  use m_Files,               only : nfout
  use m_Hubbard,             only : m_Hubbard_Potential
  use m_Control_Parameters,  only : icond
  use m_Const_Parameters,    only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION

! ================================== K. Tagami ================== 5.0
  use m_Hubbard,             only : m_Hubbard_Potential2
! =============================================================== 5.0


! ============================ added by K. Tagami ================ 11.0
  use m_Control_Parameters,   only : noncol
  use m_Electronic_Structure,  only : vlhxcQ
! ===================================================================== 11.0

  implicit none
 
! =========================== Modified by K. Tagami ============= 0.1
!!!  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) return
! ==============================================================

! =================================== modifed by K. Tagami ============- 5.0
!  call m_Hubbard_Potential(nfout)       ! -> dhub
  call m_Hubbard_Potential2(nfout)       ! -> dhub
! ====================================================================== 5.0
  

  call m_ESiVQ_add_dhub_to_vlhxcQ(nfout)! vlhxcQ = vlhxcQ + dhub
  

end subroutine Renewal_of_Hubbard_Potential
