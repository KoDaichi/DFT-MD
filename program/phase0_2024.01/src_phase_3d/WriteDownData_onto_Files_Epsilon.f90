!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.01                                       #
!#                                                                       #
!#      Sub Routine : WriteDownData_onto_Files_Eps                       #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: WriteDownData_onto_Files_Epsilon.f90 593 2019-06-20 03:47:31Z jkoga $
subroutine WriteDownData_onto_Files_Eps
  use m_Epsilon_ek,           only : nlo, wd_eps, dealloc_m_Epsilon, &
       &                             kv3_in_the_ek_process
  use m_Files,                only : m_Files_close_nfeps, m_Files_close_nfepscont, &
       &                             m_Files_close_nfnlo
  use m_Kpoints,              only : kv3

  use m_Control_Parameters,   only : sw_corelevel_spectrum, sw_local_approx_trans_moment, ipritiming
  use m_Const_Parameters,     only : ON
  use m_Timing,               only : tstatc_wd0

  implicit none

  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
     call dealloc_m_Epsilon
  else
     call wd_eps
     call dealloc_m_Epsilon
     call m_Files_close_nfeps
     call m_Files_close_nfepscont
  endif
  if(nlo/=0) call m_Files_close_nfnlo
  if(ipritiming>=2) call tstatc_wd0

end subroutine WriteDownData_onto_Files_Eps
