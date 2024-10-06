!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Transition_moment_Epsilon                          #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
!
! $Id: Transition_moment_Epsilon.f90 593 2019-06-20 03:47:31Z jkoga $
subroutine Transition_moment_Epsilon
  use m_Epsilon_ek,            only : eigen_value_ordering_eps_ek &
       &                              ,calc_transition_moment_eps_ek
  use m_Control_Parameters,  only : sw_corelevel_spectrum, sw_local_approx_trans_moment
  use m_Const_Parameters,      only : ON
  use m_CLS_dipquad,           only : m_CLS_calc_transmom_ek
  
  implicit none

  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
     call m_CLS_calc_transmom_ek
  else
     call eigen_value_ordering_eps_ek
     call calc_transition_moment_eps_ek
  endif
end subroutine Transition_moment_Epsilon
