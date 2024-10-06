!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Calc_Epsilon                                       #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: Calc_Epsilon.f90 593 2019-06-20 03:47:31Z jkoga $
subroutine Calc_Epsilon
  use m_Epsilon_ek,    only : BZintegration_eps, calc_drude_eps, kkt_eps, optics_eps

! ========= KT_add ======== 13.0S
  use m_Control_Parameters,  only : sw_corelevel_spectrum, sw_local_approx_trans_moment
  use m_Const_Parameters, only : ON
  use m_Epsilon_ek,       only : corelevel_eps
! ========================= 13.0S

  implicit none

  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
  else
     call BZintegration_eps
     call kkt_eps
  endif

  if ( sw_corelevel_spectrum == ON ) then
     call corelevel_eps
  else
     call calc_drude_eps
     call optics_eps
  endif

end subroutine Calc_Epsilon
