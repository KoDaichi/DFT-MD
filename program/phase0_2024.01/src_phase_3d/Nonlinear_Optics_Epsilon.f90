!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Calc_Nonlinear_optics                              #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: Nonlinear_Optics_Epsilon.f90 376 2014-06-17 07:48:31Z jkoga $
   subroutine Calc_Nonlinear_optics
     use m_Files,                 only : m_Files_open_nfnlo
     use m_Epsilon_ek,            only : nlo, calc_nlo
     if(nlo/=0) then
        call m_Files_open_nfnlo
        call calc_nlo
     end if
   end subroutine Calc_Nonlinear_optics
