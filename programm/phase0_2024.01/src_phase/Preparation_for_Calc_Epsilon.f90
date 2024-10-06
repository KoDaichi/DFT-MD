!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Prep_for_Calc_Epsilon                              #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: Preparation_for_Calc_Epsilon.f90 593 2019-06-20 03:47:31Z jkoga $
subroutine Prep_for_Calc_Epsilon
  use m_Epsilon_ek,      only : gen_k_points_eps_ek, calc_tm_square_eps_ek

  implicit none

  call gen_k_points_eps_ek
  call calc_tm_square_eps_ek
end subroutine Prep_for_Calc_Epsilon
