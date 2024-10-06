!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : PseudoPotential_ek_Epsilon                         #
!#                                                                       #
!#                                Written by T. Hamada 2007/5/8          #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo CISS              #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: PseudoPotential_ek_Epsilon.f90 376 2014-06-17 07:48:31Z jkoga $
subroutine PseudoPotential_ek_Epsilon
!
! Pseudopotential_ek for ek_Epsilon
! The original program is Pseudopotential_ek in PseudoPotential_Construction
!
  use m_Const_Parameters,     only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, ON
  use m_Control_Parameters,   only : icond, sw_pdos, sw_use_add_proj
  use m_Kpoints,              only : kv3, vkxyz, kv3_previous
  use m_NonLocal_Potential,   only : m_NLP_rd_snl_3D, m_NLP_betar_dot_PWs_3D &
       &                           , m_NLP_betar_dot_PWs_diff_3D &
       &                           , m_NLP_phir_dot_PWs_3D &
       &                           , m_NLP_add_betar_dot_PWs_3D
  use m_PseudoPotential,      only : m_PP_alloc_NLP, m_PP_dealloc_NLP, m_PP_betar_calculated
  use m_IterationNumbers,     only : iteration_electronic,nk_in_the_process
  use m_Files,                only : nfout, nfcntn_bin, m_Files_open_nfcntn_bin &
       &                           , F_CNTN_BIN_in_partitioned
  implicit none

  if(icond == FIXED_CHARGE &
       & .or. (icond == FIXED_CHARGE_CONTINUATION &
       &           .and.( iteration_electronic == 0 .or. kv3_previous /= kv3))) then
     call m_PP_alloc_NLP()
     call m_NLP_betar_dot_PWs_3D(nfout,kv3,vkxyz) !(kbint) --> snl
     if(sw_pdos == ON) then
        call m_NLP_phir_dot_PWs_3D(nfout,kv3,vkxyz) !(kbint) --> phig
     end if
     if(sw_use_add_proj == ON) then
        call m_NLP_add_betar_dot_PWs_3D(nfout,kv3,vkxyz) !(kbint) --> snl_add
     end if
!    call m_PP_dealloc_NLP()
!!$     if(nk_in_the_process == 3) stop 'PseudoPotential_ek'
  else if(icond == FIXED_CHARGE_CONTINUATION) then
!!$     call m_Files_open_nfcntn_bin()
     call m_NLP_rd_snl_3D(nfout,nfcntn_bin,F_CNTN_BIN_in_partitioned,kv3)
  end if
end subroutine PseudoPotential_ek_Epsilon
