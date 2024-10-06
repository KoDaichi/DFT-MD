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
  use m_Files, only : nfinp_exists
  use m_Const_Parameters, only : DP
  use m_Kpoints, only : m_Kp_set_nspin
  implicit none

  integer       :: ipri
  real(DP)      :: e1, e2      ! energy range, Efermi+e1 - Efermi+e2
  integer       :: izi, izf    ! z range top and bottom
  real(kind=DP) :: rini, rfin  ! mixing ratio in iteration
  integer       :: nfin        ! rmix = rini + (rfin - rini)*(iter-1)/nfin
  integer       :: nc_z=3      ! column position of the z-component,
                               ! which is normaly 3
  integer, dimension(3) :: n_fc,n_fc_i
  real(kind=DP) :: erlmt ! error limit
  real(kind=DP) :: wct_start = 0.d0
contains
  subroutine print_help()
    write(6,'(a)') 'work function postprocessing utility for PHASE'
    write(6,'(a)') 'usage : workfunc options'
    write(6,'(a)') 'options :' 
!    write(6,'(a)') '-s  or --spin enable this option when spin is taken into account.'
    write(6,'(a)') '-z=ZAXIS  or --zaxis=ZAXIS use this option in order to specify &
              & the index for the z-axis (one of 1, 2 or 3). defaults to 3.'
    write(6,'(a)')
  end subroutine

  subroutine m_CtrlP_rd_nc_z_axis(nfinp)
    integer, intent(in) :: nfinp
    integer :: nspin
    integer :: i,narg,ieq
    character(len=32) :: arg,argval
    if(nfinp_exists)then
       read(nfinp,*) nc_z, nspin
       call m_Kp_set_nspin(nspin)
    endif
    narg = command_argument_count()
    do i=1,narg
       call get_command_argument(i,arg)
       argval = arg 
       ieq = index(arg,'=')
       if (ieq>0) then
          if(.not.(ieq.eq.len_trim(arg)-1)) then
             write(0,*) 'invalid argument : '//trim(arg)
             stop
          endif
          argval = arg(1:ieq-1)
       endif
       select case (argval)
       case('-s','--spin')
          nspin = 2
          call m_Kp_set_nspin(nspin)
       case ('-h','--help')
          call print_help()
          stop
       case ('-z','--zaxis')
          if(ieq>0) then
              read(arg(ieq+1:ieq+1),*) nc_z
          else
              if (i==narg) then
                 stop 'the -z option takes one of : 1, 2 or 3'
              endif
              call get_command_argument(i+1,arg)
              read(arg(1:1),*) nc_z
          endif
          if(nc_z<1 .or. nc_z>3) then
             stop 'the -z option takes one of : 1, 2 or 3'
          endif
       end select
    enddo
!!    print '(" nc_z = ",i5)',nc_z
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

end module m_Control_Parameters
