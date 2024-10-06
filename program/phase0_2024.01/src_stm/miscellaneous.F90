!================================================
!  Software name : STM
!  Subroutine(s) :
!  Author(s)     : Takahiro Yamasaki (June 7, 2004)
!
!  FURTHER MODIFICATION: Takahiro Yamasaki (June 24, 2004)
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
!                           @(#)miscellaneous.F 1.1 01/03/13 00:35:39
!
! $Id: miscellaneous.F90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
#ifndef VPP
#ifdef _USE_DATE_AND_TIME_
subroutine gettod(x)
  use m_Const_Parameters, only :DP
  real(kind=DP), parameter   :: milisec     = 1.d3
  real(kind=DP), parameter   :: t_mil2micro = 1.d3
  real(kind=DP), intent(out) :: x

  integer ipresent_time(8)
  real(kind=DP)       :: tstart
  real(kind=DP), save ::tstart0
  data istart/0/
  data tstart/0.0/
  ! Month   1   2   3   4   5   6   7   8   9  10  11  12
  ! days   31  28  31  30  31  30  31  31  30  31  30  31
  !sumdays 31  59  90 120 151 181 212 243 273 304 334 365
  integer, dimension(12) :: sdays
  integer, parameter     :: days_to_1970 = 719460
  data sdays/0, 31, 59, 90,120,151,181,212,243,273,304,334/
  integer  :: days_from_January_1st, month, year, day, days_from_1970_to_this_year

  call date_and_time(values=ipresent_time)
  year = ipresent_time(1); month = ipresent_time(2); day = ipresent_time(3)
!!$  days_from_1970_to_this_year = year*365 + year/4 - year/100 + year/400 - days_to_1970
  days_from_January_1st = sdays(month) + day - 1
  if(month >= 3 .and. mod(year,4) == 0) &
       & days_from_January_1st = days_from_January_1st + 1
!!$  tstart =(((days_from_1970_to_this_year + days_from_January_1st*24 & !hours
  tstart =(((days_from_January_1st*24 & !hours
       & + ipresent_time(5))*60 &             !minutes
       & + ipresent_time(6))*60 &             !seconds
       & + ipresent_time(7))*milisec  &       !mili-seconds
       & + ipresent_time(8)
!!$  if(istart <= 3) then
!!$     write(6,'(" total_days_from_January_1st = ",i10)') total_days_from_January_1st
!!$     write(6,'(" ipresent_time(1-3) = ",3i10)') ipresent_time(1),ipresent_time(2),ipresent_time(3)
!!$     write(6,'(" ipresent_time(4-6) = ",3i10)') ipresent_time(4),ipresent_time(5),ipresent_time(6)
!!$     write(6,'(" ipresent_time(7-8) = ",3i10)') ipresent_time(7),ipresent_time(8)
!!$  end if

  if(istart.eq.0) tstart0 = tstart

  istart = istart + 1
  x = (tstart - tstart0)*t_mil2micro
  return
end subroutine gettod
#else
subroutine gettod(x)
  use m_Const_Parameters, only :DP
  real(kind=DP), intent(out) :: x
  integer :: count, count_rate, count_max
  integer, save :: count_prev
  integer, save :: initialized = 0
  integer, save :: istart = 0
  real (kind=DP)       :: tstart
  real (kind=DP), save :: tstart0 
  call system_clock(count, count_rate, count_max)
  if(count_prev > count) initialized = initialized + 1
  tstart = (count+initialized*count_max)*1.d0/count_rate
  if(istart <= 2) then
     if(istart <= 1) write(6,'(" !@ -- system_clock check --")')
     write(6,'(" !@ count, count_rate, count_max = ",3i12)') count, count_rate, count_max
     write(6,'(" !@ tstart = ",d20.12)') tstart
  end if
  if(tstart == 0) tstart0 = tstart
  istart = istart + 1
  x  = (tstart-tstart0)*1.d6
  count_prev = count
  return
end subroutine gettod
#endif
#else
subroutine dummy
end subroutine dummy
#endif

#ifdef CRAY
function derf(x)
  implicit real*8(a-h,o-z)
  derf = erf(x)
  return
end function derf

function derfc(x)
  implicit real*8(a-h,o-z)
  derfc = erfc(x)
  return
end function derfc
#endif
