!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Renewal_of_pPotential
!
!  AUTHOR(S): M. Saito, T. Yamasaki   November/04/2003
!  
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
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
!
subroutine Renewal_of_pPotential()
! $Id: Renewal_of_pPotential.f90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Files,               only : nfout
!!  use m_ES_Intgr_VlhxcQlm,   only : m_ESiVQ_integrate_VlhxcQlm
  use m_Charge_Density,      only : chgq_l
  use m_epc_potential,       only : m_epc_alloc, m_epc_cal_potential, &
       &                            m_epc_ESlhxc_potential, vepc_l, &
       &                            m_epc_ESlhxc_potential_mod
  use m_Control_Parameters,  only : positron_method
  use m_Const_Parameters,    only : Positron_CONV

  implicit none

  call m_epc_alloc()
  call m_epc_cal_potential(nfout,chgq_l) ! chgq_l -> vepc_l, tchgq_l
                                              ! (xcfft) -> vxc_l

  if ( positron_method == Positron_CONV ) then
     call m_epc_ESlhxc_potential(nfout)     ! vepc_l, tchgq_l ->  vlhxc_l
  else
     call m_epc_ESlhxc_potential_mod(nfout)
  endif

!!  call m_ESiVQ_integrate_VlhxcQlm(nfout) ! (lclchh) -> vlhxcQ

end subroutine Renewal_of_pPotential
