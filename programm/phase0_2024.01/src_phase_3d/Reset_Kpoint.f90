!#========================================================================
!#                                                                       #
!# Software Name : UVSOR ver. 3.00                                       #
!#                                                                       #
!#      Sub Routine : Reset_Kpoint                                       #
!#                                                                       #
!#                                Written by T. Hamada 2006/4/28         #
!#                                                                       #
!#      Contact address :  IIS,The University of Tokyo RSS21 project     #
!#                                                                       #
!#"Multiscale Simulation System for Functional Analysis of Nanomaterials"#
!#                                                                       #
!#========================================================================
!
! $Id: Reset_Kpoint.f90 376 2014-06-17 07:48:31Z jkoga $
   subroutine Reset_Kpoint
     use m_Files,                 only : nfout
     use m_Kpoints,               only : kv3_ek, vkxyz_ek
     use m_Epsilon_ek,            only : sw_mass, ikshift, vkxyz_ek_org
     use m_Control_Parameters,    only : printable
     implicit none
     if(sw_mass==0.or.ikshift==0.0d0) return
     vkxyz_ek(1:kv3_ek,1:3,1:2) = vkxyz_ek_org(1:kv3_ek,1:3,1:2)
     if(printable) write(nfout,'(1x,"!* k-points are reseted")')
   end subroutine Reset_Kpoint
