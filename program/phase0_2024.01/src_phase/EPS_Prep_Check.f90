!#========================================================================
!#                                                                       #
!# Software Name : PHASE/0 2019.01                                       #
!#                                                                       #
!#      Function    : EPS_prep_check                                     #
!#                                                                       #
!#                                Written by T. Hamada 2021/9/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
   logical function EPS_Prep_Check()
     use m_Epsilon_ek,            only : kpt_file_mode
     use m_Parallelization,    only : mype,MPI_CommGroup
     implicit none
     EPS_Prep_Check = .false.
     if(kpt_file_mode == 3) EPS_Prep_Check = .true.
   end function EPS_Prep_Check
