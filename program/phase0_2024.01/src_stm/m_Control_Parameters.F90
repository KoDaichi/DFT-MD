!================================================
!  Software name : STM
!  Module : m_Control_Parameters
!  Subroutine(s) : m_CtrlP_rd_nc_z_axis, m_CtrlP_rd_energy_range,
!                  m_CtrlP_rd_z_range, m_CtrlP_set_wct_start,
!                  m_CtrlP_rd_mixing_ratio, m_CtrlP_rd_error_limit
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
!
!  The license of the code and contact address :
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)
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

module m_Control_Parameters
! $Id: m_Control_Parameters.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_Const_Parameters, only : DP
  implicit none

  integer       :: ipri
  real(DP)      :: e1, e2      ! energy range, Efermi+e1 - Efermi+e2
  integer       :: izi, izf    ! z range top and bottom
  real(kind=DP) :: rini, rfin  ! mixing ratio in iteration
  integer       :: nfin        ! rmix = rini + (rfin - rini)*(iter-1)/nfin
  integer       :: nc_z        ! column position of the z-component,
                               ! which is normaly 1, sometimes 3
  integer, dimension(3) :: n_fc,n_fc_i
  real(kind=DP) :: erlmt ! error limit
  real(kind=DP) :: wct_start = 0.d0
#ifdef VPP
!                 wct_start: remaining time at start(sec).
!                 wct_now  : remaining time now(sec).
  real(kind=DP), parameter :: Critical_Remaining_CPU_TIME = 600  !(sec.)
#else
!                 wct_start: Wall Clock Time Start.
!                 wct_now  : Wall Clock Time Now.
#endif
contains
  subroutine m_CtrlP_rd_nc_z_axis(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) nc_z
    print '(" nc_z = ",i5)',nc_z
    if(nc_z .eq. 1) then
       n_fc(1) = 1; n_fc(2) = 2; n_fc(3) = 3
       n_fc_i(1) = 1; n_fc_i(2) = 2; n_fc_i(3) = 3
    else if(nc_z .eq. 2) then
       n_fc(1) = 2; n_fc(2) = 3; n_fc(3) = 1
       n_fc_i(1) = 3; n_fc_i(2) = 1; n_fc_i(3) = 2
    else if(nc_z .eq. 3) then
       n_fc(1) = 3; n_fc(2) = 1; n_fc(3) = 2
       n_fc_i(1) = 2; n_fc_i(2) = 3; n_fc_i(3) = 1
    end if
  end subroutine m_CtrlP_rd_nc_z_axis

  subroutine m_CtrlP_rd_energy_range(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) e1, e2
    print '(" Bias Voltage (V) = ", 2f10.5)', e1, e2
    if(e1 > e2 ) stop ' e1 > e2'
  end subroutine m_CtrlP_rd_energy_range

  subroutine m_CtrlP_rd_z_range(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) izi, izf
    print '(" izi = ", i5, " izf = ",i5)',izi,izf
    if(izi > izf ) stop ' !! izi > izf (m_CtrlP_rd_z_range)'
  end subroutine m_CtrlP_rd_z_range

  subroutine m_CtrlP_set_wct_start
#ifdef VPP
    integer :: iremain
    call chkelaps(iremain)
    wct_start = iremain
#else
    call gettod(wct_start)
#endif
  end subroutine m_CtrlP_set_wct_start

  subroutine m_CtrlP_rd_mixing_ratio(nfinp)
    integer, intent(in) :: nfinp
    read(nfinp,*) rini, rfin, nfin
    print '(" r_{ini}, r_{fin}, n_{fin} = ",2f10.5,i5)', rini,rfin,nfin
  end subroutine m_CtrlP_rd_mixing_ratio

  subroutine m_CtrlP_rd_error_limit(nfinp)
    integer, intent(in) :: nfinp
    read (nfinp,*) erlmt
    print '(" erlmt(error limit) = ", d20.8)', erlmt
  end subroutine m_CtrlP_rd_error_limit
    
end module m_Control_Parameters
