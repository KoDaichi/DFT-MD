!================================================
!  Software name : STM
!  Module m_Timing
!  Subroutine(s) : tstatc0_begin, wd_the_subroutine_name, tstac0_end
!                  tstatc_init, tstatc_wd0, tstatc_wd, tstatc_iter
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

module m_Timing
! $Id: m_Timing.f90,v 1.1.1.1 2004/06/26 11:21:37 yamasaki Exp $
  use m_Const_Parameters,   only : DP
  use m_Control_Parameters, only : ipri

  implicit none

  integer, parameter :: MSBRNM         = 100
  integer, parameter :: long_name_size = 32
  character*(long_name_size), dimension(0:MSBRNM) :: subroutine_names
  character*(long_name_size)                 :: temp_str
  real(DP), dimension(0:MSBRNM)              :: tstart=0, ecpu=0
  integer, dimension(MSBRNM)  :: iorder
  integer                     :: n_sub_names = 0

  real(DP)               :: ecpu_previous = 0.d0
  real(DP), parameter    :: UMICRO = 1.d-6
  real(DP), parameter    :: PCPUDF = 0.03
  integer, parameter     :: N_SUBROUTINES = 10

contains
  subroutine tstatc0_begin(a_sub_name,id)
    character*(*), intent(in) :: a_sub_name
    integer, intent(inout)    :: id

    real(kind=DP)             :: t_start
    integer                   :: i, len

    len = len_trim(a_sub_name)

! -- finding the pointer of subroutine_names --'    
    if(id <= 0) then
       do i = 1, n_sub_names
          if(a_sub_name(1:len) == subroutine_names(i)(1:len)) then
             id = i
             exit
          end if
       end do
    end if
    if(id <= 0) then
       if( n_sub_names < MSBRNM) then
          n_sub_names = n_sub_names + 1
          id          = n_sub_names
          write(subroutine_names(id),'(a32)') a_sub_name
       else
          print *,' !! Size of an array of subroutine_names should be enlarged (->m_Timing)'
       end if
    end if

    if(ipri >= 2) call wd_the_subroutine_name
!!$   if(ipri >= 2) print '(" <<< ",a32," >>>")',subroutine_names(id)

    if(id >= 1 .and. id <= MSBRNM) then
       call gettod(t_start)
       tstart(id) = t_start
    end if
  contains
    subroutine wd_the_subroutine_name
      integer :: i, ip
      ip = long_name_size - len
      ip = ip / 2
      do i = 1, ip
         temp_str(i:i) = '-'
      end do
      do i = 1, len
         temp_str(i+ip:i+ip) = a_sub_name(i:i)
      end do
      do i = ip + len + 1, long_name_size
         temp_str(i:i) = '-'
      end do
      print '(" <<<",a32,">>>")',temp_str
    end subroutine wd_the_subroutine_name
  end subroutine tstatc0_begin

  subroutine tstatc0_end(id)
    integer, intent(in) :: id
    real(kind=DP) :: t_end, t_used

    call gettod(t_end)
    t_used = t_end - tstart(id)
    ecpu(id) = ecpu(id) + t_used * UMICRO
  end subroutine tstatc0_end

  subroutine tstatc_init
    tstart(1:MSBRNM) = 0.d0; ecpu(1:MSBRNM) = 0.d0
  end subroutine tstatc_init

  subroutine tstatc_wd0
    integer i, ip, nsub
    print '(" n_sub_names = ",i5)', n_sub_names
    print '(" << cpu time statistics >>")'

    call tsort(ecpu(1),MSBRNM,iorder)

    nsub = n_sub_names

    do i = 1, nsub
       ip = iorder(i)
       if(ecpu(ip) < UMICRO) exit
       print '(i5,2x,a32,f11.5,"(sec.)")',i,subroutine_names(ip),ecpu(ip)
    end do
  end subroutine tstatc_wd0

  subroutine tstatc_wd(iteration)
    integer, intent(in) :: iteration
    integer i, ip, nsub
    real(kind=DP) :: cpudif, pecpu

    cpudif = dabs(ecpu(0) - ecpu_previous)
    if(cpudif < ecpu_previous * PCPUDF .or. ecpu(0) < UMICRO) then
!!$       print '(" ! cpudif  = ",f20.10," ecpu_prev = ",f20.10)', cpudif, ecpu_previous
!!$       print '(" ! ecpu(0) = ",f20.10)',ecpu(0)
       goto 1001
    end if

    print '(" << CPU Time Consumption -- TOP",i4," Subroutines (",i5,") >>")' &
         & , N_SUBROUTINES, iteration
    call tsort(ecpu(1),MSBRNM,iorder)

    nsub = min(n_sub_names, N_SUBROUTINES)
    do i = 1, nsub
       ip = iorder(i)
       if(ecpu(ip) < UMICRO) exit
       pecpu = ecpu(ip)/ecpu(0) * 100
       print '(2i4,2x,a32,f11.5,"(sec.)",f6.2,"(%)")' &
            & ,i,ip, subroutine_names(ip),ecpu(ip),pecpu
    end do
    print '(6x,"Total cpu time of this iteration",4x,f11.5,"(sec.)")',ecpu(0)

1001 ecpu_previous = ecpu(0)
  end subroutine tstatc_wd

  subroutine tstatc_iter(iteration, first_iteration_of_this_job)
    integer, intent(in) :: iteration, first_iteration_of_this_job
    real(kind=DP) :: t_end, t_used
    call gettod(t_end)
    if(iteration > first_iteration_of_this_job) then
       t_used  = t_end - tstart(0)
       ecpu(0) = t_used * UMICRO
    end if
    tstart(0) = t_end
  end subroutine tstatc_iter

end module m_Timing
