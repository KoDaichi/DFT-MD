!================================================
!  Software name : STM
!  Subroutine(s) : Preparation
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
!
!  Contact address :  IIS,The University of Tokyo RSS21 project
!  
!  "Multiscale Simulation System for Function Analysis of Nanomaterials"  
!
!================================================
!
!     The original version of this program "STM" was developed in the
!  period 1998-2001 at JRCAT by Koichi Kato (Toshiba Corporate
!  Research and Development Center) and Takahiro Yamasaki (Fujitsu 
!  Laboratories Ltd.) through a joint research between JRCAT and
!  Toshiba Corporate Research and Development Center and another joint
!  research between JRCAT and FUJITSU Laboratories Ltd.
!     Since 2003, this program has been tuned and revised as it
!  matches to the first-principles molecular dynamics program "PHASE"
!  as a part of the national project "Frontier Simulation Software for
!  Industrial Science (FSIS)", which is supported by the IT program of
!  the Ministry of Education, Culture, Sports, Science and Technology
!  (MEXT) of Japan.
!     Since 2006, this program set has been developed as a part of the
!  national project "Revolutionary Simulation Software (RSS21)", which
!  is supported by the next-generation IT program of MEXT of Japan.

subroutine Preparation
! $Id: Preparation.f90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
  use m_Const_Parameters,   only : ON
  use m_Files,              only : nfout
  use m_Control_Parameters, only : m_CtrlP_set_wct_start
  use m_FFT,  only :               fft_box_size_WF, fft_box_size_CD, fft_box_size_fine &
       &                         , m_FFT_query_inversion_symmetry &
       &                         , m_FFT_setup, m_FFT_set_fft_box_size_fine
  use m_PlaneWaveBasisSet, only  : m_pwBS_set_FFT_mapfunctions

  implicit none
  integer :: inversion_symmetry = ON

  call m_CtrlP_set_wct_start()

  call m_FFT_query_inversion_symmetry(inversion_symmetry) !-(m_FFT)
  call m_FFT_set_fft_box_size_fine
  call m_FFT_setup(nfout,inversion_symmetry)
  call m_pwBS_set_FFT_mapfunctions&
       & (fft_box_size_WF, fft_box_size_CD, fft_box_size_fine)

end subroutine Preparation
