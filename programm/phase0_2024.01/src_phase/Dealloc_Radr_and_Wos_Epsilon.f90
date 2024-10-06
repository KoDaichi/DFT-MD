!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Dealloc_Radr_and_Wos_Epsilon                       #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: Dealloc_Radr_and_Wos_Epsilon.f90 376 2014-06-17 07:48:31Z jkoga $
   subroutine Dealloc_Radr_and_Wos_Epsilon
    use m_PseudoPotential,             only : m_PP_dealloc_NLP
     call m_PP_dealloc_NLP
   end subroutine Dealloc_Radr_and_Wos_Epsilon 
