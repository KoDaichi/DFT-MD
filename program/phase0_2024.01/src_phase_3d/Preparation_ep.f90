!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Preparation_ep
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!     The original version of this set of the computer programs "PHASE"
!  was developed by the members of the Theory Group of Joint Research
!  Center for Atom Technology (JRCAT), based in Tsukuba, in the period
!  1993-2001.
!
!     Since 2002, this set has been tuned and new functions have been
!  added to it as a part of the national project "Frontier Simulation 
!  Software for Industrial Science (FSIS)",  which is supported by
!  the IT program of the Ministry of Education, Culture, Sports,
!  Science and Technology (MEXT) of Japan. 
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.
!   Since 2013, this program set has been further developed centering on PHASE System
!  Consortium.
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
subroutine Preparation_ep
! $Id: Preparation_ep.f90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Control_Parameters,only: m_CtrlP_set_wct_start
  use m_Files,             only: nfout
  use m_Crystal_Structure, only: nbztyp_spg, m_CS_alloc_op_tau &
       &                       , m_CS_gnrt_symmetry_operations
  use m_Ionic_System,      only: m_IS_alloc_napt &
       &                       , m_IS_symm_check_of_pos
  use m_Force,             only: m_Force_alloc

  implicit none

  call m_CtrlP_set_wct_start                   !  (ckcput)
!!$  if(nbztyp_spg >= GENERAL) call m_Files_open_
!!$  call m_CS_gnrt_symmetry_operations(.true.,nfout,nfspg)
  call m_CS_gnrt_symmetry_operations(.true.,nfout)
!!$  if(.not.paramset) then
     call m_CS_alloc_op_tau(nfout)
!!$     call m_CS_gnrt_symmetry_operations(.false.,nfout,nfspg)
     call m_CS_gnrt_symmetry_operations(.false.,nfout)
     call m_IS_alloc_napt()
     call m_IS_symm_check_of_pos()
!!$  end if
  call m_Force_alloc()
end subroutine Preparation_ep
