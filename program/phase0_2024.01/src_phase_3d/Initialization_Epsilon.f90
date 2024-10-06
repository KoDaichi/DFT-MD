!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 593 $)
!
!  SUBROUINE: Initialization_Epsilon
!
!  AUTHOR(S): T. Hamada   November/24/2004
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
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
!
!
! $Id: Initialization_Epsilon.f90 593 2019-06-20 03:47:31Z jkoga $
subroutine Initialization_Epsilon(mode)
  use m_Epsilon_ek,           only : initialization_eps_ek
  use m_Control_Parameters,   only : m_CtrlP_set_uvsormode_ON
  use m_Files,                only : m_Files_open_nfeps, m_Files_open_nfepscont

  use m_Control_Parameters,   only : sw_corelevel_spectrum, sw_local_approx_trans_moment
  use m_Const_Parameters,     only : ON

  implicit none
  integer, intent(in) :: mode

  call initialization_eps_ek(mode)

  if ( sw_corelevel_spectrum == ON .and. sw_local_approx_trans_moment == ON ) then
  else
     call m_Files_open_nfeps()
     call m_Files_open_nfepscont
  endif

  call m_CtrlP_set_uvsormode_ON()

end subroutine Initialization_Epsilon
