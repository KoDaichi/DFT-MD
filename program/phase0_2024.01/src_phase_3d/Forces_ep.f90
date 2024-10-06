!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Forces
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
subroutine Forces
! $Id: Forces_ep.f90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Control_Parameters, only : m_CtrlP_what_is_mdalg
  use m_Force,              only : forc_l, m_Force_stlweb
  use m_Ionic_System,       only : m_IS_evaluate_v_verlet
  implicit none
  integer :: mdalg

  call m_Force_stlweb()
  mdalg = m_CtrlP_what_is_mdalg()
  call m_IS_evaluate_v_verlet(mdalg,forc_l)
end subroutine Forces
