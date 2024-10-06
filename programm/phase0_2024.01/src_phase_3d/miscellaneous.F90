!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: gettod(x), dummy
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
! $Id: miscellaneous.F90 570 2017-04-21 20:34:50Z yamasaki $
#ifndef VPP
#ifdef _USE_DATE_AND_TIME_
subroutine gettod(x)
  use m_Const_Parameters, only :DP
  real(kind=DP), parameter   :: milisec     = 1.d3
  real(kind=DP), parameter   :: t_mil2micro = 1.d3
  real(kind=DP), intent(out) :: x

  integer ipresent_time(8)
  real(kind=DP)       :: tstart
  real(kind=DP), save ::tstart0 = 0.d0
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
  integer, save :: count_prev = -1
  integer, save :: initialized = 0
  integer, save :: istart = 0
  real (kind=DP)       :: tstart
  real (kind=DP), save :: tstart0  = 0.d0
!------------------------------------------------------------------------
  character(len= 8) :: cdate
  character(len=10) :: ctime
  integer :: idate
  real(kind=DP)       :: rtime
  real(kind=DP), save :: rtime_prev = 0
  call date_and_time(date = cdate, time = ctime)
  read(cdate, '(i8)') idate
  read(ctime, '(f10.3)') rtime
  rtime = (idate - 20000000) * 1.d6 + rtime
!------------------------------------------------------------------------
  call system_clock(count, count_rate, count_max)
  if(count_prev > count) then
     if(rtime > rtime_prev) then
        initialized = initialized + 1
!        write(6,*) 'DEBUG:', count
     else
        write(6,*) 'Caution!  : Clock tamperring is detected'
        write(6,*) 'DEBUG:', count
     endif
  endif
  tstart = (count+initialized*real(count_max, kind = DP))/count_rate
!!$  if(istart <= 2) then
!!$     if(istart <= 1) write(6,'(" !@ -- system_clock check --")')
!!$     write(6,'(" !@ count, count_rate, count_max = ",3i12)') count, count_rate, count_max
!!$     write(6,'(" !@ tstart = ",d20.12)') tstart
!!$  end if
  if(istart == 0) tstart0 = tstart
  istart = istart + 1
  x  = (tstart-tstart0)*1.d6
  count_prev = count
  rtime_prev = rtime
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

#elif Linux

#ifdef _RPL_ERF_

function derf(x) result(y)
use m_Const_Parameters, only : PAI
implicit none

real(kind=8) :: y
real(kind=8), intent(in) :: x

! local variables
real(kind=8), parameter :: pi   = PAI                      ! PI
real(kind=8), parameter :: pi4  = 4.d0*pi                  ! 4*PI                    
real(kind=8), parameter :: ipi  = 0.31830988618379067153d0 ! 1/PI
real(kind=8), parameter :: sipi = 0.56418958354775628694d0 ! 1/sqrt(PI)
real(kind=8), parameter :: sipi2 = 2.d0*sipi               ! 2/sqrt(PI)
real(kind=8) :: z

z = abs(x)
if(z <= 0.1d0) then
  y = derf1(z)
else if(z < 6) then
  y = derf2(z)
else if(z <= 100) then
  y = derf3(z)
else
!  y = derf4(z)
  y = 1.d0
end if

if(x < 0.d0) y = -y

return

contains
  ! for 0 <= x <= 0.1 
  function derf1(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: c1 = -1.d0/3.d0
  real(kind=8), parameter :: c2 =  0.1d0
  real(kind=8), parameter :: c3 = -1.d0/42.d0
  real(kind=8), parameter :: c4 =  1.d0/216.d0
  real(kind=8), parameter :: c5 = -1.d0/1320.d0
  real(kind=8) :: x2
  real(kind=8) :: z1,z2,z3,z4,z5
 
  x2 = x*x

  z1 = x2*x
  z2 = z1*x2
  z3 = z2*x2
  z4 = z3*x2
  z5 = z4*x2

  y = sipi2*( x + c1*z1 + c2*z2 + c3*z3 + c4*z4 + c5*z5 )

  return
  end function derf1

  ! for 0.1 <= x <= 6
  function derf2(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  y = derf3(x) + 2.d0/(exp(pi4*x)-1.d0)

  return
  end function derf2

  ! for 6 < x <= 100
  function derf3(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: h(13) = (/ 0.25d0, 1.d0, 2.25d0, 4.d0, 6.25d0, 9.d0, 12.25d0, 16.d0, &
                           & 20.25d0, 25.d0, 30.25d0, 36.d0, 42.25d0 /)
  real(kind=8), parameter :: e(13) = (/ 0.7788007830714049d0, &
                           & 0.3678794411714423d0, &
                           & 0.1053992245618643d0, &
                           & 0.1831563888873418d-1, &
                           & 0.1930454136227709d-2, &
                           & 0.1234098040866796d-3, &
                           & 0.4785117392129009d-5, &
                           & 0.1125351747192591d-6, &
                           & 0.1605228055185612d-8, &
                           & 0.1388794386496402d-10, &
                           & 0.7287724095819693d-13, &
                           & 0.2319522830243569d-15, & 
                           & 0.4477732441718302d-18 /)
  real(kind=8) :: x2

  x2 = x*x
  
  y = sum(e(1:13)/(h(1:13)+x2))
  y = ipi*exp(-x2)*(0.5d0/x+x*y)
  y = 1-y

  return
  end function derf3

  ! for 100 < x
  function derf4(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: c1 = -0.5d0
  real(kind=8), parameter :: c2 =  3.d0/4.d0
  real(kind=8), parameter :: c3 = -15.d0/8.d0
  real(kind=8), parameter :: c4 =  105.d0/16.d0
  real(kind=8), parameter :: c5 = -945.d0/32.d0
  real(kind=8) :: x2,xi
  real(kind=8) :: z1,z2,z3,z4,z5

  x2 = x*x
  xi = 1.d0/x
  
  z1 = xi*xi
  z2 = z1*z1
  z3 = z2*z1
  z4 = z3*z1
  z5 = z4*z1
  
  y = sipi*exp(-x2)*xi*( 1.d0 + c1*z1 + c2*z2 + c3*z3 + c4*z4 + c5*z5 )
  y = 1-y

  return
  end function derf4

end function derf
function derfc(x) result(y)
use m_Const_Parameters, only : PAI
implicit none

real(kind=8) :: y
real(kind=8), intent(in) :: x

! local variables
real(kind=8), parameter :: pi   = PAI                      ! PI
real(kind=8), parameter :: pi4  = 4.d0*pi                  ! 4*PI                    
real(kind=8), parameter :: ipi  = 0.31830988618379067153d0 ! 1/PI
real(kind=8), parameter :: sipi = 0.56418958354775628694d0 ! 1/sqrt(PI)
real(kind=8), parameter :: sipi2 = 2.d0*sipi               ! 2/sqrt(PI)
real(kind=8) :: z

z = abs(x)
if(z <= 0.1d0) then
  y = derfc1(z)
else if(z < 6) then
  y = derfc2(z)
else if(z <= 100) then
  y = derfc3(z)
else
  !y = derfc4(z)
  y = 0.d0
end if

if(x < 0.d0) y = 2.d0-y

return

contains
  ! for 0 <= x <= 0.1 
  function derfc1(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: c1 = -1.d0/3.d0
  real(kind=8), parameter :: c2 =  0.1d0
  real(kind=8), parameter :: c3 = -1.d0/42.d0
  real(kind=8), parameter :: c4 =  1.d0/216.d0
  real(kind=8), parameter :: c5 = -1.d0/1320.d0
  real(kind=8) :: x2
  real(kind=8) :: z1,z2,z3,z4,z5
 
  x2 = x*x

  z1 = x2*x
  z2 = z1*x2
  z3 = z2*x2
  z4 = z3*x2
  z5 = z4*x2

  y = sipi2*( x + c1*z1 + c2*z2 + c3*z3 + c4*z4 + c5*z5 )
  y = 1-y

  return
  end function derfc1

  ! for 0.1 <= x <= 6
  function derfc2(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  y = derfc3(x) - 2.d0/(exp(pi4*x)-1.d0)

  return
  end function derfc2

  ! for 6 < x <= 100
  function derfc3(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: h(13) = (/ 0.25d0, 1.d0, 2.25d0, 4.d0, 6.25d0, 9.d0, 12.25d0, 16.d0, &
                           & 20.25d0, 25.d0, 30.25d0, 36.d0, 42.25d0 /)
  real(kind=8), parameter :: e(13) = (/ 0.7788007830714049d0, &
                           & 0.3678794411714423d0, &
                           & 0.1053992245618643d0, &
                           & 0.1831563888873418d-1, &
                           & 0.1930454136227709d-2, &
                           & 0.1234098040866796d-3, &
                           & 0.4785117392129009d-5, &
                           & 0.1125351747192591d-6, &
                           & 0.1605228055185612d-8, &
                           & 0.1388794386496402d-10, &
                           & 0.7287724095819693d-13, &
                           & 0.2319522830243569d-15, & 
                           & 0.4477732441718302d-18 /)
  real(kind=8) :: x2

  x2 = x*x
  
  y = sum(e(1:13)/(h(1:13)+x2))
  y = ipi*exp(-x2)*(0.5d0/x+x*y)

  return
  end function derfc3

  ! for 100 < x
  function derfc4(x) result(y)
  implicit none

  real(kind=8) :: y
  real(kind=8), intent(in) :: x

  ! local variables
  real(kind=8), parameter :: c1 = -0.5d0
  real(kind=8), parameter :: c2 =  3.d0/4.d0
  real(kind=8), parameter :: c3 = -15.d0/8.d0
  real(kind=8), parameter :: c4 =  105.d0/16.d0
  real(kind=8), parameter :: c5 = -945.d0/32.d0
  real(kind=8) :: x2,xi
  real(kind=8) :: z1,z2,z3,z4,z5

  x2 = x*x
  xi = 1.d0/x
  
  z1 = xi*xi
  z2 = z1*z1
  z3 = z2*z1
  z4 = z3*z1
  z5 = z4*z1
  
  y = sipi*exp(-x2)*xi*( 1.d0 + c1*z1 + c2*z2 + c3*z3 + c4*z4 + c5*z5 )

  return
  end function derfc4

end function derfc

#endif

#endif
