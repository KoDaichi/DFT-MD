!================================================
!  Software name : STM
!  Subroutine(s) : InputData_Analysis
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
!
!  FURTHER MODIFICATION: Junichiro Koga (June 24, 2004)
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

subroutine InputData_Analysis
! $Id: InputData_Analysis.F90,v 1.1.1.1 2004/06/26 11:21:36 yamasaki Exp $
  use m_Files, only                : nfinp,nfcntn_bin
  use m_Control_Parameters, only   : m_CtrlP_rd_energy_range &
       &                           , m_CtrlP_rd_z_range &
       &	                   , m_CtrlP_rd_mixing_ratio &
       &                           , m_CtrlP_rd_error_limit &
       &                           , m_CtrlP_rd_nc_z_axis
  use m_Charge_File,          only : m_CF_rd_Header_For_Cube
  use m_ArraySize_Parameters, only : m_ArraySize_Parameters_rd
  use m_Electronic_Structure, only : m_ES_rd_EigenValues_etc
  use m_PlaneWaveBasisSet, only    : m_pwBS_rd_data
  use m_FFT, only                  : m_FFT_rd_fft_box_size &
	&                          , m_FFT_rd_fine_fft_box_size
  use m_Kpoints, only              : m_Kp_rd_kv3, m_Kp_alloc, m_Kp_rd_kpoints
  use m_Crystal_Structure, only    : m_CS_rd_neg_and_zl

  implicit none

  call m_Kp_rd_kv3(nfinp)
  call m_Kp_alloc()
  call m_Kp_rd_kpoints(nfinp)
  call m_CS_rd_neg_and_zl(nfinp)
  call m_CtrlP_rd_energy_range(nfinp)
  call m_FFT_rd_fine_fft_box_size(nfinp)
  call m_CtrlP_rd_z_range(nfinp)
  call m_CtrlP_rd_nc_z_axis(nfinp)
  call m_CtrlP_rd_mixing_ratio(nfinp)
  call m_CtrlP_rd_error_limit(nfinp)
  call m_ArraySize_Parameters_rd(nfcntn_bin)
  call m_ES_rd_EigenValues_etc(nfcntn_bin)
  call m_pwBS_rd_data(nfcntn_bin)
  call m_FFT_rd_fft_box_size(nfcntn_bin)
  call m_CF_rd_Header_For_Cube(nfcntn_bin)

end subroutine InputData_Analysis
