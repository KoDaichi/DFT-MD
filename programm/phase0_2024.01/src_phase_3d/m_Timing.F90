!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  MODULE:  m_Timing
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki, January/13/2004
!  FURTHER MODIFICATION: T. Yamasaki, September 2009
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
!
!tex--\documentstyle{jarticle}
!tex--\title{module m\_Timing}
!tex--\author{T. Yamasaki}
!tex--\textwidth 480pt
!tex--\oddsidemargin -10pt
!tex--\evensidemargin 0pt
!tex--\begin{document}
!tex--\maketitle
!tex--\newenvironment{mydescription}{%
!tex--	\list{$\bullet$}{%
!tex--		\topsep 0pt \itemsep 0pt }}%
!tex--{\endlist}
!tex--
module m_Timing
!
!
! ***************************************************************
!  $Id: m_Timing.F90 570 2017-04-21 20:34:50Z yamasaki $
!
!  The first version of this module was coded by T. Yamasaki
! (JRCAT-ATP, FUJITSU LABORATORIES Ltd.) in 1999.
!
!tex--\section{機能概要}
!tex-- サブルーチン毎、あるいは１-iteraionの時間を測定し、表示する。
!tex--\section{Using modules}
  use m_Control_Parameters, only : ipri, ipritiming, printable,num_subroutines, pcpudf &
       &                         , statistics_in_parallel, sw_flatten, sw_firstlevel_only, sw_details &
       &                         , measure_count_limit
  use m_Const_Parameters,   only : DP, ON, OFF, START,ITERATIVE,FINISH,TAG_FORMAT,TAG_LINE,TABLE,ERROR
  use m_Parallelization,    only : MPI_CommGroup,npes, mype, ierr
  use m_Files,              only : nfout
  use mpi

  implicit none

!tex--\section{Parameters,variables,arrays}
!tex--\subsection{Public}
!tex--   none
!tex--\subsection{Private}
!tex--\begin{itemize}
!tex--      \item {\tt [MSBRNM]:} 時間測定をするsubroutineの最大数
!tex--      \item {\tt [long\_name\_size]:} 配列{\tt subroutine\_names}の各要素の制限文字数
!tex--      \item {\tt [subroutine\_names]:} suburoutine名を蓄える配列
!tex--      \item {\tt [tstart] (0:MSBRNM):} 測定開始時刻を蓄える配列
!tex--      \item {\tt [ecpu] (0:MSBRNM):} 測定した経過時間を蓄える配列。{\tt ecpu(0)}に
!tex--            はSCF１回分の経過時間が蓄えられる。
!tex--      \item {\tt [iorder] (MMSBRM):} 測定時間の大きさの順序を与える配列
!tex--      \item {\tt [n\_sub\_names]:} 測定したsubroutineの数
!tex--      \item {\tt [PCPUDF = 0.03]:} SCF-iterationの経過時間の前回の経過時間
!tex-- ({\tt ecpu\_previous})との差が{\tt PCPUDF}$\times${\tt ecpu\_previous}より大きいときに
!tex-- 測定結果を書き出す。{\tt subroutine tstatc\_wd}内で用いる。
!tex--      \item {\tt [N\_SUBROUTINES] = 10:} {\tt subroutine tstatc\_wd}で書き出すsubroutineの数の上限
!tex--\end{itemize}
!  integer, private, parameter                            :: MSBRNM         = 100
  integer, private                                       :: MSBRNM         = 100
  integer, private, parameter                            :: MSBRNM_EXT     = 10
  integer, private, parameter                            :: INCRE_MSBRNM   = 10
  integer, private, parameter                            :: long_name_size = 32
!  character(len=long_name_size),private,dimension(0:MSBRNM) :: subroutine_names
!  integer, private, dimension(0:MSBRNM)                  :: subroutine_name_length
  character(len=long_name_size),private                     :: temp_str
!  real(DP),private, dimension(0:MSBRNM)                  :: tstart=0, ecpu=0
!  integer, private, dimension(MSBRNM)                    :: iorder, counter=0
!  integer, private, dimension(MSBRNM)                    :: ilevel=-1
  character(len=long_name_size),private,allocatable      :: subroutine_names(:)
  integer, private, allocatable                          :: subroutine_name_length(:)
  real(DP),private, allocatable                          :: tstart(:), ecpu(:)
  integer, private, allocatable                          :: iorder(:), counter(:)
  integer, private, allocatable                          :: ilevel(:)
  integer, private                                       :: levelset = OFF
  integer, private                                       :: n_sub_names = 0
  real(DP),private               :: ecpu_previous = 0.d0
  real(DP),private, parameter    :: UMICRO = 1.d-6
  integer, private :: switch_of_timing = OFF
  real(DP),private :: wct_start = 0.d0
! --> T. Yamasaki, 2th Sep. 2009
  integer, private :: id_ext, id_ext_max = 0, n_allocated_subhir = 0
  integer, private, dimension(MSBRNM_EXT) :: id_ext_stack=0
  integer, private :: id_ext_pointer = 0
  type t_subhir
     integer :: id
     integer :: level
     integer :: id1
     integer :: id2
     integer :: id3
     integer :: flag
  end type t_subhir
  type(t_subhir), allocatable,  dimension(:) :: subhir
  integer,private :: level_previous,level_now=0, id123_now(1:3)=0
  real(kind=DP),private,allocatable,dimension(:) :: ecpu_ext
  integer, private, allocatable, dimension(:) :: counter_ext
  integer, private, allocatable, dimension(:) :: ip_1stlevel_subroutines
  integer :: n_1stlevel_subroutines
  integer, private, parameter :: MODE_EVALUATION = 1, MODE_WRITE = 2, MODE_SINGLE=1, MODE_PARALLEL = 2
  real(kind=DP), private :: tmp_timer_start_time=0.d0
! <--

!!$  real(DP),private, parameter    :: PCPUDF = 0.03
!!$  integer, private, parameter    :: N_SUBROUTINES = 20

!tex--
!tex--\section{含まれるsubroutines}
!tex-- {\tt a)tstatc0\_begin, b)tstatc0\_end, c)tstatc\_init, d)tstatc\_wd0, e)tstatc\_wd, f)tstatc\_iter}
!tex--\section{使い方}
!tex--測定したいsubroutineでは、最初に次の例のように{\tt tstatc0_begin}を呼び出し
!tex--\begin{verbatim}
!tex--    integer :: id_sname = -1
!tex--    call tstatc0_begin('FFT_WF ',id_sname)
!tex--\end{verbatim}
!tex--最後に
!tex--\begin{verbatim}
!tex--    call tstatc0_end(id_sname)
!tex--\end{verbatim}
!tex--と{\tt tstatc0_end}を呼び出す。これにより、{\tt FFT_WF}に対応する
!tex--配列に{\tt tstatc0_begin}と{\tt tstatc0_end}の間の経過時間が加えられる。
!tex--SCFの１回の経過時間は、{\tt tstatc_iter}を呼び出すことで測られる。
!tex--SCF-iterationの終りに{\tt tstatc_wd}を呼び出すと、前回の
!tex--SCF１回分の経過時間との差が3％({\tt PCPUDF})以上あれば、経過時間
!tex--の長かった方から10個(={\tt N_SUBROUTINES})の
!tex--subroutineの名前と経過時間、またその回の経過時間を出力する。
!tex--具体的には、SCF-iterationの終りに次のようにする。
!tex--\begin{verbatim}
!tex--  call tstatc_iter(iteration, first_iteration_of_this_job)
!tex--  if(iteration == first_iteration_of_this_job) then
!tex--     call tstatc_wd0
!tex--  else
!tex--     call tstatc_wd(iteration)
!tex--  end if
!tex--  call tstatc_init
!tex--\end{verbatim}
!tex--{\tt call tstatc_wd(iteration)}により次のような出力を得る。
!tex--\begin{verbatim}
!tex-- << CPU Time Consumption -- TOP  10 Subroutines (    2) >>
!tex--   1  26                           FFT_WF     1.79200(sec.) 52.58(%)
!tex--   2  31      evolve_WFs_in_MSD_direction     1.04297(sec.) 30.60(%)
!tex--   3  17                     xc_potential     0.94531(sec.) 27.74(%)
!tex--   4  19               ggaxcp(in xc_pot.)     0.89258(sec.) 26.19(%)
!tex--   5  28              energy_eigen_values     0.57715(sec.) 16.93(%)
!tex--   6  33                      CD_softpart     0.49121(sec.) 14.41(%)
!tex--   7  36   Total_Energy(including xc_pot)     0.47168(sec.) 13.84(%)
!tex--   8  18                   inverse_FFT_CD     0.23242(sec.)  6.82(%)
!tex--   9  20                    direct_FFT_CD     0.14356(sec.)  4.21(%)
!tex--  10  34                      CD_hardpart     0.13574(sec.)  3.98(%)
!tex--      Total cpu time of (        2 )-th iteration        3.40820(sec.)
!tex--\end{verbatim}
!tex--\end{document}
!  include 'mpif.h'
contains

  subroutine alloc_subroutine_names_etc(n,nprev)
    integer, intent(in) :: n, nprev
    character(len=long_name_size),allocatable :: subroutine_namest(:)
    integer,  allocatable                     :: subroutine_name_lengtht(:)
    real(DP), allocatable                     :: tstartt(:), ecput(:)
    integer,  allocatable                     :: iordert(:), countert(:)
    integer,  allocatable                     :: ilevelt(:)
    if(nprev<n) then
      allocate(subroutine_namest(0:n))
      allocate(subroutine_name_lengtht(0:n))
      allocate(tstartt(0:n));tstartt=0
      allocate(ecput(0:n));ecput=0
      allocate(iordert(n))
      allocate(countert(n));countert=0
      allocate(ilevelt(n));ilevelt=-1
      subroutine_namest(0:nprev)       = subroutine_names(0:nprev)
      subroutine_name_lengtht(0:nprev) = subroutine_name_length(0:nprev)
      tstartt(0:nprev)                 = tstart(0:nprev)
      ecput(0:nprev)                   = ecpu(0:nprev)
      iordert(1:nprev)                 = iorder(1:nprev)
      countert(1:nprev)                = counter(1:nprev)
      ilevelt(1:nprev)                 = ilevel(1:nprev)
    endif
    if(allocated(subroutine_names))       deallocate(subroutine_names)
    if(allocated(subroutine_name_length)) deallocate(subroutine_name_length)
    if(allocated(tstart))                 deallocate(tstart)
    if(allocated(ecpu))                   deallocate(ecpu)
    if(allocated(iorder))                 deallocate(iorder)
    if(allocated(counter))                deallocate(counter)
    if(allocated(ilevel))                 deallocate(ilevel)
    allocate(subroutine_names(0:n))
    allocate(subroutine_name_length(0:n))
    allocate(tstart(0:n));tstart=0
    allocate(ecpu(0:n));ecpu=0
    allocate(iorder(n))
    allocate(counter(n));counter=0
    allocate(ilevel(n));ilevel=-1
    if(nprev<n) then
      subroutine_names(0:nprev)       = subroutine_namest(0:nprev)
      subroutine_name_length(0:nprev) = subroutine_name_lengtht(0:nprev)
      tstart(0:nprev)                 = tstartt(0:nprev)
      ecpu(0:nprev)                   = ecput(0:nprev)
      iorder(1:nprev)                 = iordert(1:nprev)
      counter(1:nprev)                = countert(1:nprev)
      ilevel(1:nprev)                 = ilevelt(1:nprev)
      deallocate(subroutine_namest)
      deallocate(subroutine_name_lengtht)
      deallocate(tstartt)
      deallocate(ecput)
      deallocate(iordert)
      deallocate(countert)
      deallocate(ilevelt)
    endif
  end subroutine alloc_subroutine_names_etc

  subroutine tstatc0_begin(a_sub_name,id,level)
    character(len=*), intent(in) :: a_sub_name
    integer, intent(inout) ::       id
    integer, intent(in), optional :: level

    real(kind=DP)             :: t_start
    integer                   :: i, len_str
    if(.not.allocated(subroutine_names)) then
      call alloc_subroutine_names_etc(MSBRNM,MSBRNM)
    endif

    if(switch_of_timing == OFF) then
       call gettod(wct_start)
       switch_of_timing = ON
    end if

    if(sw_firstlevel_only == ON .and. sw_details == OFF) then
         if(.not.present(level)) then
           return
         else
           if(level/=1) return
         end if
    end if
!      & (.not.present(level).or.(present(level).and.level/=1))) return

    len_str = len_trim(a_sub_name)
    if(len_str > long_name_size) len_str = long_name_size

! -- finding the pointer in subroutine_names --'
    if(id <= 0) then
       do i = 1, n_sub_names
!!$          if(a_sub_name(1:len_str) == subroutine_names(i)(1:len_str)) then
          if(  len_str == subroutine_name_length(i) .and. &
               & a_sub_name(1:len_str) == subroutine_names(i)(1:len_str)) then
             id = i
              exit
          end if
       end do
    end if
    if(id <= 0) then
       if( n_sub_names>=MSBRNM) then
         call alloc_subroutine_names_etc(MSBRNM+INCRE_MSBRNM,MSBRNM)
         if(ipritiming>=2 .and. printable) write(nfout, &
         '(" !! array size for the subroutine_names has been enlarged from",i6," to ",i6,"(->m_Timing)")') &
         MSBRNM,MSBRNM+INCRE_MSBRNM
         MSBRNM = MSBRNM+INCRE_MSBRNM
       endif
       if( n_sub_names < MSBRNM) then
          n_sub_names = n_sub_names + 1
          id          = n_sub_names
          write(subroutine_names(id),'(a32)') a_sub_name
          subroutine_name_length(id) = len_str
       else
          if(printable) write(nfout,'(" !! Size of an array of subroutine_names should be enlarged (->m_Timing)")')
       end if
    end if


    if(id >= 1 .and. id <= MSBRNM) then
       if(measure_count_limit > 0 .and. counter(id) > measure_count_limit) then
       else
          call gettod(t_start)
          tstart(id) = t_start
       end if
       if(sw_firstlevel_only == ON .or. (sw_firstlevel_only==OFF .and. sw_flatten==OFF)) then
          if(present(level) .and. ilevel(id)==-1) then
             ilevel(id) = level
             if(level/=-1 .and. levelset==OFF) levelset = ON
          end if
       end if
    end if

    if(sw_details == ON) then
       call alloc_subhir()
       call setlevels() ! --> level_now, level_previous, id123_now
       call search_id_ext() ! --> id_ext
       if(id_ext > id_ext_max) then
          subhir(id_ext)%id    = id
          subhir(id_ext)%level = level_now
          if(level_now==1) subhir(id_ext)%id1 = 0
          if(level_now>=2) subhir(id_ext)%id1 = id123_now(1)
          if(level_now>=3) subhir(id_ext)%id2 = id123_now(2)
          if(level_now>=4) subhir(id_ext)%id3 = id123_now(3)
          id_ext_max = id_ext
       end if
!!$       id_ext_pointer = min(id_ext_pointer + 1, MSBRNM_EXT)
!!$       id_ext_pointer = id_ext_pointer + 1
       id_ext_pointer = level_now
       if(id_ext_pointer <= MSBRNM_EXT) id_ext_stack(id_ext_pointer) = id_ext
    end if

    if(ipritiming >= 2 .and. printable) call wd_the_subroutine_name
!!$   if(ipri >= 2) print '(" <<< ",a32," >>>")',subroutine_names(id)
  contains
    subroutine search_id_ext()
      integer :: i
      id_ext = 0
      if(id_ext_max == 0) then
         if(id>=1) then
            id_ext = 1
         end if
      else
         if(level_now <= 1) then
            do i = 1, id_ext_max
               if(subhir(i)%id == id) then
                  id_ext = i
                  exit
               end if
            end do
         else
            do i = 1, id_ext_max
               if(subhir(i)%id == id .and. subhir(i)%level == level_now) then
!!$                  if((subhir(i)%level == level_now).and.(subhir(i)%id1 == id123_now(1))) then
                  if(subhir(i)%id1 == id123_now(1)) then
                     id_ext = i
                     exit
                  end if
               end if
            end do
         end if
         if(id_ext == 0) id_ext = id_ext_max+1
      end if
    end subroutine search_id_ext

    subroutine setlevels()
      if(present(level)) then
         level_previous = level_now
         level_now = level
         if(level_now >= 1 .and. level_now <= 3) id123_now(level_now) = id
      else
         if(level_now >= 0) then
            level_previous = level_now
            level_now = level_now+1
            if(level_now >= 1 .and. level_now <= 3) id123_now(level_now) = id
         end if
      end if
    end subroutine setlevels

    subroutine alloc_subhir()
      type(t_subhir), allocatable,  dimension(:) :: subhir2
      real(kind=DP),allocatable, dimension(:) :: ecpu2
      integer, allocatable, dimension(:) :: counter2
      if(.not.allocated(subhir)) then
         allocate(subhir(MSBRNM_EXT))
         do i = 1, MSBRNM_EXT
            subhir(i)%id = 0
            subhir(i)%level = 0
            subhir(i)%id1 = 0
            subhir(i)%id2 = 0
            subhir(i)%id3 = 0
            subhir(i)%flag = 0
         end do
         allocate(ecpu_ext(MSBRNM_EXT)); ecpu_ext = 0.d0
         allocate(counter_ext(MSBRNM_EXT)); counter_ext = 0
         n_allocated_subhir = MSBRNM_EXT
      else if(id_ext_max >= n_allocated_subhir) then
         allocate(subhir2(id_ext_max))
         allocate(ecpu2(id_ext_max))
         allocate(counter2(id_ext_max))
         if(allocated(subhir)) then
            subhir2 = subhir
            deallocate(subhir)
            n_allocated_subhir = n_allocated_subhir+INCRE_MSBRNM
            allocate(subhir(n_allocated_subhir))
            subhir(1:min(id_ext_max,n_allocated_subhir)) = subhir2(1:min(id_ext_max,n_allocated_subhir))
            do i = min(id_ext_max, n_allocated_subhir)+1, n_allocated_subhir
               subhir(i)%id = 0
               subhir(i)%level = 0
               subhir(i)%id1 = 0
               subhir(i)%id2 = 0
               subhir(i)%id3 = 0
               subhir(i)%flag = 0
            end do
         end if
         if(allocated(ecpu_ext)) then
            ecpu2 = ecpu_ext
            deallocate(ecpu_ext)
            allocate(ecpu_ext(n_allocated_subhir)); ecpu_ext = 0.d0
            ecpu_ext(1:min(id_ext_max,n_allocated_subhir)) = ecpu2(1:min(id_ext_max,n_allocated_subhir))
         end if
         if(allocated(counter_ext)) then
            counter2 = counter_ext
            deallocate(counter_ext)
            allocate(counter_ext(n_allocated_subhir)); counter_ext = 0
            counter_ext(1:min(id_ext_max,n_allocated_subhir)) = counter2(1:min(id_ext_max,n_allocated_subhir))
         end if
         deallocate(counter2)
         deallocate(ecpu2)
         deallocate(subhir2)
      end if
    end subroutine alloc_subhir

    subroutine wd_the_subroutine_name
      integer :: i, ip
      ip = long_name_size - len_str
      ip = ip / 2
      do i = 1, ip
         temp_str(i:i) = '-'
      end do
      do i = 1, len_str
         temp_str(i+ip:i+ip) = a_sub_name(i:i)
      end do
      do i = ip + len_str + 1, long_name_size
         temp_str(i:i) = '-'
      end do
      write(nfout,'(" <<<",a32,">>> (",i3,"), level_now = ",i3)') temp_str,ilevel(id), level_now
    end subroutine wd_the_subroutine_name
  end subroutine tstatc0_begin

  subroutine tstatc0_end(id)
    integer, intent(in) :: id
    real(kind=DP) :: t_end, t_used

    if(id==-1)  return

    if(measure_count_limit > 0 .and. counter(id) >measure_count_limit) then
       t_used = 0.d0
    else
       call gettod(t_end)
       t_used = t_end - tstart(id)
       ecpu(id) = ecpu(id) + t_used * UMICRO
    end if
    counter(id) = counter(id) + 1

    if(sw_details == ON) then
       if(level_now >= 1) then
          level_previous = level_now
          level_now = level_now-1
       end if
       if(id_ext_pointer >=1 .and. id_ext_pointer <= MSBRNM_ext) then
          id_ext = id_ext_stack(id_ext_pointer)
          if(id_ext >=1 .and. id_ext <= n_allocated_subhir) then
             ecpu_ext(id_ext) = ecpu_ext(id_ext) + t_used *UMICRO
             counter_ext(id_ext) = counter_ext(id_ext) + 1
          end if
       else
          id_ext = 0
       end if
!!$       id_ext_pointer = id_ext_pointer-1
       id_ext_pointer = level_now
       if(ipritiming >= 2 .and. printable) call wd_the_subroutine_name2
    end if

    if(ipritiming >= 2 .and. printable) call wd_the_subroutine_name
  contains
    subroutine wd_the_subroutine_name
      write(nfout,'(" >>>",a32,"<<<")') subroutine_names(id)
    end subroutine wd_the_subroutine_name
    subroutine wd_the_subroutine_name2
      write(nfout,'(" >>>",a32,"<<< id_ext, id_ext_pointer = ",i5,i5)') subroutine_names(id), id_ext, id_ext_pointer
    end subroutine wd_the_subroutine_name2

  end subroutine tstatc0_end

  subroutine tstatc_init
    tstart(1:MSBRNM) = 0.d0; ecpu(1:MSBRNM) = 0.d0; counter(1:MSBRNM) = 0
    if(sw_details == ON) then
       if(allocated(ecpu_ext)) ecpu_ext = 0.d0
       if(allocated(counter_ext)) counter_ext = 0
       id_ext_max = 0
    end if
  end subroutine tstatc_init

  subroutine tstatc_wd0
    integer :: i, ip, nsub, in, nsub0, npmax, npmin
    real(kind=DP) :: cpu_used, wct_now &
         & , cputmax, cputmin, cputmean, cputvariance, cputstddevia
    real(kind=DP), allocatable, dimension(:,:) :: ecpu_wk, ecpu_sorted
    integer, allocatable, dimension(:) :: iorder_p
    integer :: ipri0, lev, levmax, inu, j, num_sub_wd, cput_mode

    ipri0 = ipritiming
    call mpi_bcast(ipri0,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(ipri0<=0) return

    if(measure_count_limit > 0) call correct_ecpu() ! counter,measure_count_limit,ecpu-> ecpu
    call tsort(ecpu(1),MSBRNM,iorder)

    nsub = n_sub_names

    if(levelset == ON) then
       levmax = 2
       if(sw_firstlevel_only==ON) levmax = 1
    else
       levmax = 1
    end if

    if(statistics_in_parallel == 1 .and. npes >= 2) then
       cput_mode = MODE_PARALLEL
       if(mype == 0) then
          nsub0 = 0
          do i = 1, nsub
             if(ecpu(iorder(i)) < UMICRO) exit
             nsub0 = i
          end do
       end if
       call mpi_bcast(nsub0,1,mpi_integer,0,MPI_CommGroup,ierr)
       allocate(iorder_p(nsub0))
       if(mype == 0) iorder_p(1:nsub0) = iorder(1:nsub0)
       call mpi_bcast(iorder_p,nsub0,mpi_integer,0,MPI_CommGroup,ierr)

       allocate(ecpu_wk(0:npes-1,nsub0)); ecpu_wk = 0.d0
       allocate(ecpu_sorted(0:npes-1,nsub0))
       do i = 1, nsub0
          ip = iorder_p(i)
          ecpu_wk(mype,i) = ecpu(ip)
       end do
       call mpi_allreduce(ecpu_wk, ecpu_sorted,nsub0*npes, mpi_double_precision &
            &            ,mpi_sum, MPI_CommGroup,ierr)

    else
       cput_mode = MODE_SINGLE
       nsub0 = 0
       do i = 1, nsub
          if(ecpu(iorder(i)) < UMICRO) exit
          nsub0 = i
       end do
    end if

    do j = MODE_EVALUATION, MODE_WRITE
       if(j==MODE_EVALUATION) num_sub_wd = 0
       if(j==MODE_WRITE) call wd_headerlines()

       do lev = 1, levmax
          if(j==MODE_WRITE) call wd_pauseline()
          inu = 0
          do i = 1, nsub0
             if(cput_mode == MODE_PARALLEL) then
                ip = iorder_p(i)
             else if(cput_mode == MODE_SINGLE) then
                ip = iorder(i)
             end if

             if(levelset == ON) then
                if(levmax >= 2) then
                   if(.not.((lev==levmax .and. ilevel(ip)==-1) .or. (lev<levmax .and. ilevel(ip)==lev))) cycle
                else if(levmax == 1) then
                   if(ilevel(ip)/=lev) cycle
                end if
             end if
             inu = inu + 1

             if(j==MODE_WRITE) then
                if(cput_mode == MODE_PARALLEL) then
                   call cpu_statistics_in_parallel(ecpu_sorted(0,i) &
                        & ,cputmin,cputmax,cputmean,cputvariance,cputstddevia,npmax,npmin)
                   if(printable) write(nfout,'(2i5,2x,a32,f11.5,i8,4f10.4,3i5)') &
                        & i,ip,subroutine_names(ip),ecpu(ip),counter(ip) &
                        & ,cputmax,cputmin,cputmean,cputstddevia,npmax,npmin,inu
                else
                   if(printable) write(nfout,'(2i5,2x,a32,f11.5,2i8)') &
                        & i,ip,subroutine_names(ip),ecpu(ip),counter(ip), inu
                end if
             end if
          end do
          if(j==MODE_EVALUATION) num_sub_wd = num_sub_wd+inu
       end do
    end do
    if(cput_mode == MODE_PARALLEL) deallocate(iorder_p,ecpu_sorted,ecpu_wk)

!!$    else
!!$       if(printable) then
!!$          nsub0 = 0
!!$          do i = 1, nsub
!!$             if(ecpu(iorder(i)) < UMICRO) exit
!!$             nsub0 = i
!!$          end do
!!$
!!$          do j = MODE_EVALUATION, MODE_WRITE
!!$             if(j==MODE_EVALUATION) num_sub_wd = 0
!!$             if(j==MODE_WRITE) call wd_headerlines()
!!$
!!$             do lev = 1, levmax
!!$                if(j==MODE_WRITE) call wd_pauseline()
!!$                inu = 0
!!$                do i = 1, nsub0
!!$                   ip = iorder(i)
!!$                   if(levelset == ON) then
!!$                      if(levmax >= 2) then
!!$                         if(.not.((lev==levmax .and. ilevel(ip)==-1) &
!!$                              & .or. (lev<levmax .and. ilevel(ip)==lev))) cycle
!!$                      else if(levmax == 1) then
!!$                         if(ilevel(ip)/=lev) cycle
!!$                      end if
!!$                   end if
!!$                   inu = inu+1
!!$                   if(j==MODE_WRITE .and. printable) then
!!$                      write(nfout,'(2i5,2x,a32,f11.5,2i8)') &
!!$                           & i,ip,subroutine_names(ip),ecpu(ip),counter(ip), inu
!!$                   end if
!!$                end do
!!$                if(j==MODE_EVALUATION) num_sub_wd = num_sub_wd+inu
!!$             end do
!!$          end do
!!$!!$       if(sw_details == ON) then
!!!$          if(printable) then
!!!$             write(nfout,'(" --- detail flat---")')
!!!$             write(nfout,'(" id_ext_max = ", i8)') id_ext_max
!!!$             if(printable) write(nfout,'("   no",2x," id ","        subroutine name         ","  time(sec) "&
!!!$                  &               ,"  count","   level","   id1")')
!!!$             do i = 1, min(id_ext_max,100)
!!!$                inu = subhir(i)%id
!!!$                if(inu >= 1) write(nfout,'(i5,2x,i4,a32,f11.5,3i8)') i, inu,subroutine_names(inu) &
!!!$                     & , ecpu_ext(i), counter_ext(i), subhir(i)%level, subhir(i)%id1
!!!$             end do
!!!$          end if
!!!$       end if
!!$       end if
!!$    end if
    call gettod(wct_now)
    cpu_used = (wct_now - wct_start) * UMICRO
    if(printable) write(nfout,'(" <<Total elapsed CPU Time until now =",f12.5," (sec.)>>")') cpu_used
    if(sw_details == ON .and. printable) call wd_ecpu_details(0)

  contains
    subroutine wd_pauseline()
      if(levelset == ON .and. levmax>1 .and. printable) then
         if(lev < levmax) then
            write(nfout,'("    --- level ",i3," ---")') lev
         else
            write(nfout,'("    --- ---")')
         end if
      end if
    end subroutine wd_pauseline
!!$    subroutine wd_headerlines_parallel()
!!$      if(printable) then
!!$         write(nfout,'(" n_sub_names = ",i5," num_subroutines_statistics = ",i5)') n_sub_names, num_sub_wd
!!$         write(nfout,'(" << cpu time statistics >>")')
!!$      end if
!!$    end subroutine wd_headerlines_parallel
    subroutine wd_headerlines()
      if(printable) then
         write(nfout,'(" n_sub_names = ",i5," num_subroutines_statistics = ",i5)') n_sub_names, num_sub_wd
         if(num_sub_wd >=1) then
            write(nfout,'(" << cpu time statistics >>")')
            if(cput_mode==MODE_SINGLE) then
               write(nfout,'("  no "," id ",2x,"        subroutine name         ","  time(sec) "&
                    &               ,"  count","   no(2)")')
            else if(cput_mode == MODE_PARALLEL) then
               write(nfout,'("   no","   ip",2x,"        subroutine name         ","  time(sec) "&
                    &               ,"  count  ", " max_time    min       mean     stddevia" &
                    &  ," nmax nmin  no(2)")')
            end if
         end if
      end if
    end subroutine wd_headerlines

  end subroutine tstatc_wd0

  subroutine cpu_statistics_in_parallel(ecpu_sorted &
       & ,cputmin,cputmax,cputmean,cputvariance,cputstddevia,npmax,npmin)
    real(kind=DP),intent(in),dimension(0:npes-1) :: ecpu_sorted
    real(kind=DP),intent(out)     :: cputmin,cputmax,cputmean,cputvariance,cputstddevia
    integer, intent(out) :: npmax,npmin
    integer :: in

    cputmax = ecpu_sorted(0); npmax = 0
    cputmin = ecpu_sorted(0); npmin = 0
    do in = 1, npes-1
       if(cputmax < ecpu_sorted(in)) then
          cputmax = ecpu_sorted(in);   npmax = in
       end if
       if(cputmin > ecpu_sorted(in)) then
          cputmin = ecpu_sorted(in);   npmin = in
       end if
    end do
!!$          cputmax = maxval(ecpu_sorted(0:npes-1,ic))
!!$          cputmin = minval(ecpu_sorted(0:npes-1,ic))
    cputmean = 0.d0
    do in = 0, npes-1
       cputmean = cputmean + ecpu_sorted(in)
    end do
    cputmean = cputmean/npes
    cputvariance = 0.d0
    do in =0, npes-1
       cputvariance = cputvariance + (ecpu_sorted(in)-cputmean)*(ecpu_sorted(in)-cputmean)
    end do
    cputvariance = cputvariance/npes
    cputstddevia = dsqrt(cputvariance)

  end subroutine cpu_statistics_in_parallel

  subroutine wd_ecpu_details(iteration)
    integer, intent(in) :: iteration
    integer :: i, ip, ic, inu, ips, id, icflag

    n_1stlevel_subroutines = 0
    do i = 1, n_sub_names
       if(ilevel(i) == 1) n_1stlevel_subroutines = n_1stlevel_subroutines+1
    end do
    if(n_1stlevel_subroutines ==0) return
    allocate(ip_1stlevel_subroutines(n_1stlevel_subroutines))
    n_1stlevel_subroutines = 0
    do i = 1, n_sub_names
       ip = iorder(i)
       if(ilevel(ip) == 1) then
          n_1stlevel_subroutines = n_1stlevel_subroutines+1
          ip_1stlevel_subroutines(n_1stlevel_subroutines) = ip
       end if
    end do

    write(nfout,'("--- cputiming details BEGIN ( iteration = ",i8," ) ----")') iteration
    write(nfout,'("   no "," id ",2x,"        subroutine name         ","  time(sec) "&
         &               ,"  count","   level","   id1")')
    subhir(1:id_ext_max)%flag = 0

    icflag = 0
    do i = 1, n_1stlevel_subroutines
       ip = ip_1stlevel_subroutines(i)
       ic = 0
       do inu = 1, id_ext_max
          ips = subhir(inu)%id1
          if(subhir(inu)%level>=2 .and. ips == ip) then
             if(ecpu_ext(inu) > UMICRO) ic = ic+1
             ic = ic+1
          end if
       end do
       if(ic >=1) then
          write(nfout,'(" ",60("="))')
          write(nfout,'("  *  ",i5,2x,a32,f11.5,2i8,"    *")') ip,subroutine_names(ip), ecpu(ip),counter(ip),ilevel(ip)
          write(nfout,'(" ",60("-"))')
!!$          write(nfout,'("  ----------------------------------------")')
          ic = 0
          do inu = 1, id_ext_max
             ips = subhir(inu)%id1
             if(subhir(inu)%level>=2 .and. ips == ip) then
                if(ecpu_ext(inu) > UMICRO) then
                   ic = ic+1
                   id = subhir(inu)%id
                   if(id >=1) then
                      write(nfout,'(2i5,2x,a32,f11.5,3i8)') ic,id,subroutine_names(id) &
                           & , ecpu_ext(inu),counter_ext(inu),subhir(inu)%level, subhir(inu)%id1
                   else
                      write(nfout,'("id = ",i8," <=0")') id
                   end if
                   subhir(inu)%flag = 1
                   icflag = icflag + 1
                end if
             end if
          end do
       end if
    end do
    if(icflag < id_ext_max) then
       ic = 0
       do i = 1, id_ext_max
          if(subhir(i)%flag == 0 .and. subhir(i)%level >= 2) then
             if(ecpu_ext(i) > UMICRO) then
                ic = ic + 1
             end if
          else if(subhir(i)%flag == 0 .and. subhir(i)%level==1) then
             id = subhir(i)%id
             if(1<=id .and. id <=MSBRNM .and. ilevel(id) /= 1) then
                if(ecpu_ext(i) > UMICRO) then
                   ic = ic + 1
                end if
             end if
          end if
       end do
       if(ic>=1) then
          write(nfout,'(" === other subroutines ",40("="))')
          ic = 0
          do i = 1, id_ext_max
             icflag = 0
             if(subhir(i)%flag == 0 .and. subhir(i)%level >= 2) then
                if(ecpu_ext(i) > UMICRO) then
                   ic = ic + 1
                   id = subhir(i)%id
                   icflag = 1
                end if
             else if(subhir(i)%flag == 0 .and. subhir(i)%level<=1) then
                id = subhir(i)%id
                if(1<=id .and. id <=MSBRNM .and. ilevel(id) /= 1) then
                   if(ecpu_ext(i) > UMICRO) then
                      ic = ic + 1
                      icflag = 1
                   end if
                end if
             end if
             if(icflag == 1)  write(nfout,'(2i5,2x,a32,f11.5,3i8)') ic, id,subroutine_names(id) &
                  & , ecpu_ext(i), counter_ext(i), subhir(i)%level, subhir(i)%id1
          end do
       end if
    end if
    write(nfout,'("--- cputiming details END ( iteration = ",i8," ) ----")') iteration

    deallocate(ip_1stlevel_subroutines)
  end subroutine wd_ecpu_details

  subroutine tstatc_wd(iteration)
    integer, intent(in) :: iteration
!!$    integer i, ip, nsub, in, ic, nsub0, npmax, npmin
    integer i, ip, nsub, ic, nsub0, npmax, npmin
    real(kind=DP) :: cpudif, pecpu, wct_now, cputmax, cputmin, cputmean, cputvariance, cputstddevia
    real(kind=DP), allocatable, dimension(:,:) :: ecpu_wk, ecpu_sorted
    integer, allocatable, dimension(:) :: iorder_p
    integer :: ipri0, lev, levmax, inu, j, num_sub_wd,  cput_mode

    ipri0 = ipritiming
    call mpi_bcast(ipri0,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(ipri0<=0) return

    ipri0 = 0
    cpudif = dabs(ecpu(0) - ecpu_previous)
    if(cpudif < ecpu_previous * PCPUDF .or. ecpu(0) < UMICRO) ipri0=1
    call mpi_bcast(ipri0,1,mpi_integer,0,MPI_CommGroup,ierr)
    if(ipri0>0) then
       goto 1001
    end if

    if(measure_count_limit > 0) call correct_ecpu() ! counter,measure_count_limit,ecpu-> ecpu
    call tsort(ecpu(1),MSBRNM,iorder)

    if(levelset == ON) then
       levmax = 2
       if(sw_firstlevel_only==ON) levmax = 1
    else
       levmax = 1
    end if
    nsub = 0
    if(sw_firstlevel_only == ON .and. levmax == 1) then
      ic = 0
       do i = 1, n_sub_names
          ip = iorder(i)
          if(ecpu(ip) < UMICRO) exit
          if(ilevel(ip) == 1) then
             ic = ic + 1
             if(ic > num_subroutines) exit
             nsub = i
          end if
       end do
    else
       nsub = min(n_sub_names, num_subroutines)
    end if

!!$    if(printable) write(nfout,'(" nsub = ",i8)') nsub

    if(statistics_in_parallel == 1 .and. npes >= 2) then
       cput_mode = MODE_PARALLEL
       call mpi_barrier(MPI_CommGroup,ierr)
       if(mype == 0) then
          nsub0 = 0
          do i = 1, nsub
             if(ecpu(iorder(i)) < UMICRO) exit
             nsub0 = i
          end do
       end if
       call mpi_bcast(nsub0,1,mpi_integer,0,MPI_CommGroup,ierr)
       allocate(iorder_p(nsub0))
       if(mype == 0) iorder_p(1:nsub0) = iorder(1:nsub0)
       call mpi_bcast(iorder_p,nsub0,mpi_integer,0,MPI_CommGroup,ierr)

       allocate(ecpu_wk(0:npes-1,0:nsub0)); ecpu_wk = 0.d0
       allocate(ecpu_sorted(0:npes-1,0:nsub0))
       ecpu_wk(mype,0) = ecpu(0)
       do i = 1, nsub0
          ip = iorder_p(i)
          ecpu_wk(mype,i) = ecpu(ip)
       end do
       call mpi_allreduce(ecpu_wk, ecpu_sorted,(nsub0+1)*npes, mpi_double_precision &
            &            ,mpi_sum, MPI_CommGroup,ierr)
       subroutine_names(0) = " total cpu of this iteration "
       subroutine_name_length(0) = 30
       nsub0 = nsub0 + 1
    else
       cput_mode = MODE_SINGLE
       nsub0 = 0
       do i = 1, nsub
          if(ecpu(iorder(i)) < UMICRO) exit
          nsub0 = i
       end do
    end if

    do j = MODE_EVALUATION, MODE_WRITE
       if(j==MODE_EVALUATION) num_sub_wd = 0
       if(j==MODE_WRITE) call wd_headerlines()

       do lev = 1, levmax
          if(j==MODE_WRITE) call wd_pauseline()
          inu = 0
          do i = 1, nsub0
             if(cput_mode == MODE_PARALLEL) then
                if(i <= nsub0-1) then
                   ip = iorder_p(i)
                   ic = i
                else if(i == nsub0) then
                   ip = 0
                   ic = 0
                end if
             else if(cput_mode == MODE_SINGLE) then
                ip = iorder(i)
             end if

!!$             if(levelset == ON .and. ip/=0) then
             if(levelset == ON .and. ip>0) then
                if(levmax >= 2) then
                   if(.not.((lev==levmax .and. ilevel(ip)==-1) .or. (lev<levmax .and. ilevel(ip)==lev))) cycle
                else if(levmax == 1) then
                   if(ilevel(ip)/=lev) cycle
                end if
!!$             else if(levelset == ON .and. ip == 0) then
             else if(levelset == ON .and. ip <= 0) then
                if(lev < levmax) cycle
             end if
             inu = inu+1

             if(j==MODE_WRITE) then
                if(cput_mode == MODE_PARALLEL) then
                   call cpu_statistics_in_parallel(ecpu_sorted(0,ic) &
                        & ,cputmin,cputmax,cputmean,cputvariance,cputstddevia,npmax,npmin)
                end if
                pecpu = ecpu(ip)/ecpu(0) * 100

                if(printable) then
                   if(cput_mode == MODE_PARALLEL) then
                      if(ip /= 0) then
                         write(nfout,'(2i4,2x,a32,f10.4,f6.2,i9,4f10.4,3i5)') &
                              & ic,ip, subroutine_names(ip),ecpu(ip),pecpu,counter(ip) &
                              & ,cputmax,cputmin,cputmean,cputstddevia,npmax,npmin,inu
                      else if(ip == 0) then
                         write(nfout,'(2i4,2x,a32,f10.4,f6.1,i9,4f10.4,3i5)') &
                              & ic,ip, subroutine_names(ip),ecpu(ip),pecpu,1 &
                              & ,cputmax,cputmin,cputmean,cputstddevia,npmax,npmin,inu
                      end if
                   else
                      write(nfout,'(2i4,2x,a32,f11.5,f6.2,2i8)') &
                           & i,ip, subroutine_names(ip),ecpu(ip),pecpu,counter(ip),inu
                   end if
                end if
             end if
          end do
          if(j==MODE_EVALUATION) num_sub_wd = num_sub_wd+inu
       end do
    end do
    if(cput_mode == MODE_PARALLEL)  deallocate(iorder_p,ecpu_sorted,ecpu_wk)
!!$    else
!!$       nsub0 = 0
!!$       do i = 1, nsub
!!$          if(ecpu(iorder(i)) < UMICRO) exit
!!$          nsub0 = i
!!$       end do
!!$
!!$       do j = MODE_EVALUATION, MODE_WRITE
!!$          if(j==MODE_EVALUATION) num_sub_wd = 0
!!$          if(j==MODE_WRITE) call wd_header()
!!$          do lev = 1, levmax
!!$             if(j==MODE_WRITE) call wd_pauseline()
!!$
!!$             inu = 0
!!$             do i = 1, nsub0
!!$                ip = iorder(i)
!!$                if(levelset == ON) then
!!$                   if(levmax >= 2) then
!!$                      if(.not.((lev==levmax .and. ilevel(ip)==-1) .or. (lev<levmax .and. ilevel(ip)==lev))) cycle
!!$                   else if(levmax == 1) then
!!$                      if(ilevel(ip)/=lev) cycle
!!$                   end if
!!$                end if
!!$                inu = inu+1
!!$                if(j==MODE_WRITE) then
!!$                   pecpu = ecpu(ip)/ecpu(0) * 100
!!$
!!$                   if(printable) write(nfout,'(2i4,2x,a32,f11.5,f6.2,2i8)') &
!!$                        & i,ip, subroutine_names(ip),ecpu(ip),pecpu,counter(ip),inu
!!$                end if
!!$             end do
!!$             if(j==MODE_EVALUATION) num_sub_wd = num_sub_wd+inu
!!$          end do
!!$       end do
!!$!!!$       if(sw_details == ON .and. printable) then
!!$!!!$          write(nfout,'(" --- detail ---")')
!!$!!!$          write(nfout,'(" id_ext_max = ", i8)') id_ext_max
!!$!!!$          if(printable) write(nfout,'("   no",2x," id ","        subroutine name         ","  time(sec) "&
!!$!!!$               &               ,"  count","   level","   id1")')
!!$!!!$          do i = 1, min(id_ext_max,100)
!!$!!!$             inu = subhir(i)%id
!!$!!!$             if(inu >= 1) write(nfout,'(i5,2x,i4,a32,f11.5,3i8)') i, inu,subroutine_names(inu) &
!!$!!!$                  & , ecpu_ext(i), counter_ext(i), subhir(i)%level, subhir(i)%id1
!!$!!!$          end do
!!$!!!$       end if
!!$    end if
    call gettod(wct_now)
    cpudif = (wct_now - wct_start)*UMICRO
!!$    write(nfout,'(6x,"Total cpu time of ( ",i8," )-th iteration",4x,f11.5 &
!!$         & ,"(sec.)")') iteration, ecpu(0)
    if(printable) then
       if(iteration < 10000) then
          if(cpudif < 10**6) then
             write(nfout,'(1x,"Total cputime of ( ",i4," )-th iteration",3x,f11.5," /",f10.3 &
                  & ," (sec.)")') iteration, ecpu(0), cpudif
          else
             write(nfout,'(1x,"Total cputime of ( ",i4," )-th iteration",3x,f11.5," /",f12.3 &
                  & ," (sec.)")') iteration, ecpu(0), cpudif
          end if
       else
          if(cpudif < 10**6) then
             write(nfout,'(1x,"Total cputime of ( ",i8," )-th iteration",1x,f11.5," /",f10.3 &
                  & ," (sec.)")') iteration, ecpu(0), cpudif
          else
             write(nfout,'(1x,"Total cputime of ( ",i8," )-th iteration",1x,f11.5," /",f12.3 &
                  & ," (sec.)")') iteration, ecpu(0), cpudif
          end if
       end if
    end if

    if(sw_details == ON .and. printable) call wd_ecpu_details(iteration)
1001 ecpu_previous = ecpu(0)

  contains
    subroutine wd_pauseline()
      if(levelset == ON .and. levmax>1 .and. printable) then
         if(lev < levmax) then
            write(nfout,'("    --- level ",i3," ---")') lev
         else
            write(nfout,'("    --- ---")')
         end if
      end if
    end subroutine wd_pauseline

    subroutine wd_headerlines()
      integer :: num_sub
      if(printable) then
         if(cput_mode == MODE_PARALLEL) then
            num_sub = num_sub_wd-1
         else
            num_sub = num_sub_wd
         end if
         write(nfout,'(" << CPU Time Consumption -- TOP",i4," Subroutines (",i5,") >>")') &
                 &     num_sub, iteration
         if(cput_mode == MODE_PARALLEL) then
            write(nfout,'("  no "," id ",2x,"        subroutine name         ","time(sec) "&
                 &               ," r(%) ","   count  ", "max_time    min       mean     stddevia"," nmax nmin  no(2)")')
         else
            write(nfout,'("  no "," id ",2x,"        subroutine name         "," time(sec) "&
                 &               ," r(%) ","   count","   no(2)")')
         end if
      end if
    end subroutine wd_headerlines

!!$    subroutine wd_headerlines()
!!$      if(printable) then
!!$         write(nfout,'(" << CPU Time Consumption -- TOP",i4," Subroutines (",i5,") >>")') &
!!$              &     num_sub_wd, iteration
!!$         write(nfout,'("  no "," id ",2x,"        subroutine name         "," time(sec) "&
!!$              &               ," r(%) ","   count","   no(2)")')
!!$      end if
!!$    end subroutine wd_headerlines
  end subroutine tstatc_wd

  subroutine tstatc_iter(iteration, first_iteration_of_this_job)
    integer, intent(in) :: iteration, first_iteration_of_this_job
    real(kind=DP) :: t_end, t_used
#ifdef NEC_ITER_REG
    character(len=7) :: str_reg
#endif
    call gettod(t_end)
    if(iteration > first_iteration_of_this_job) then
#ifdef NEC_ITER_REG
       write(str_reg, '(A3,I4.4)') 'IT_',iteration
       call FTRACE_REGION_END(str_reg)
#endif
       t_used  = t_end - tstart(0)
       ecpu(0) = t_used * UMICRO
    end if
    tstart(0) = t_end
#ifdef NEC_ITER_REG
       write(str_reg, '(A3,I4.4)') 'IT_',iteration+1
       call FTRACE_REGION_BEGIN(str_reg)
#endif
  end subroutine tstatc_iter

  subroutine m_Timing_wd_timenow(chars)
    character(len=*), intent(in) :: chars

    integer, dimension(8) :: ipresent_time
    integer :: i, j, ilen
    integer, parameter :: max_chars = 80, max_chars0 = 55
    character(len=80) :: aline
    character(len=20) :: dateandtime = ''

    call date_and_time(values=ipresent_time)

    ilen = len_trim(chars)
    if(ilen > max_chars0) ilen = max_chars0

!!$    aline = ''
!!$    do i = 1, ilen
!!$       aline(i:i) = chars(i:i)
!!$    end do
!!$    write(aline(ilen+6:max_chars),'(",i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4)') &
!!$         & ipresent_time(5:7),ipresent_time(3),ipresent_time(2),ipresent_time(1)
    write(dateandtime,'(i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4)') &
         & ipresent_time(5:7),ipresent_time(3),ipresent_time(2),ipresent_time(1)
    do i = 1,20
       if(i >= 9 .and. i <= 10) cycle
       if(dateandtime(i:i) == ' ' ) dateandtime(i:i) = '0'
    end do

    aline = chars(1:ilen)//'     '//dateandtime

    write(nfout,'(a80)') aline

  end subroutine m_Timing_wd_timenow

  subroutine m_Timing_init_timer()
#ifdef HIUX
    real(kind=DP) :: time0
#endif
    call gettod(wct_start)
    switch_of_timing = ON
#ifdef HIUX
    call xclock(time0,7)
    call xclock(time0,3)
#endif
  end subroutine m_Timing_init_timer

  subroutine m_Timing_wd_status(nfstatus, jobstatus_format,jobstatus_series,status_wdmode &
       &       , iteration, iteration_ionic, iteration_electronic)
    integer, intent(in) :: nfstatus, jobstatus_format, jobstatus_series, status_wdmode &
         &               , iteration, iteration_ionic, iteration_electronic
    integer :: switch_header
    character(len=9) :: modename
    real(kind=DP) :: elapsed_time

#ifdef HIUX
    real(kind=DP) :: cpu_time
    call xclock(elapsed_time,8)
    call xclock(cpu_time,5)
#else
    real(kind=DP) :: wct_now
    call gettod(wct_now)
    elapsed_time = (wct_now - wct_start) * UMICRO
#endif
    if(status_wdmode == START) then
       modename = 'START'
    else if(status_wdmode == ITERATIVE) then
       modename = 'ITERATIVE'
    else if(status_wdmode == FINISH) then
       modename = 'FINISHED'
    else if(status_wdmode == ERROR) then
       modename = 'ERROR'
    end if
    switch_header = ON
    if(jobstatus_format == TAG_FORMAT .or. jobstatus_format == TAG_LINE) then
       switch_header = OFF
    else if(jobstatus_format == TABLE) then
       if(jobstatus_series == ON .and. status_wdmode /= START) switch_header = OFF
    else
       if(printable) write(nfout,'(" jobstatus_format (= ",i5,") is illegal")') jobstatus_format
    end if
    if(jobstatus_format == TABLE) then
#ifdef HIUX
       if(switch_header == ON) &
            & write(nfstatus,'("status     iteration  iter_ionic  iter_elec  elapsed_time  cpu_time")')
       write(nfstatus,'(a9,i11,i11,i11,f14.4,f14.4)') modename, iteration &
            & , iteration_ionic, iteration_electronic, elapsed_time, cpu_time
#else
       if(switch_header == ON) &
            & write(nfstatus,'("status     iteration  iter_ionic  iter_elec  elapsed_time")')
       write(nfstatus,'(a9,i11,i11,i11,f14.4)') modename, iteration &
            & , iteration_ionic, iteration_electronic, elapsed_time
#endif
    else if(jobstatus_format == TAG_FORMAT) then
       write(nfstatus,'(" status       = ",a14)') modename
       write(nfstatus,'(" iteration    = ",i14)') iteration
       write(nfstatus,'(" iter_ionic   = ",i14)') iteration_ionic
       write(nfstatus,'(" iter_elec    = ",i14)') iteration_electronic
       write(nfstatus,'(" elapsed_time = ",f14.4)') elapsed_time
#ifdef HIUX
       write(nfstatus,'(" cpu_time     = ",f14.4)') cpu_time
#endif
    else if(jobstatus_format == TAG_LINE) then
#ifdef HIUX
       write(nfstatus,'(" status = ",a9,", iteration = ",i10,", iter_ionic = ",i10 &
            & ,", iter_elec = ",i10,", elapsed_time = ",f14.4,", cpu_time = ",f14.4)') &
            & modename,iteration,iteration_ionic,iteration_electronic &
            &,elapsed_time, cpu_time
#else
       write(nfstatus,'(" status = ",a9,", iteration = ",i10,", iter_ionic = ",i10 &
            & ,", iter_elec = ",i10,", elapsed_time = ",f14.4)') &
            & modename,iteration,iteration_ionic,iteration_electronic &
            &,elapsed_time
#endif
    end if

  end subroutine m_Timing_wd_status

  subroutine correct_ecpu()
    !   input  : ecpu, counter, measure_count_limit
    !      condition : measure_count_limit > 0
    !   output : ecpu, ecpu_ext
    integer :: i
    real(kind=DP) :: d, ecpu0
    if(measure_count_limit <= 0) return
    d = 1.d0/measure_count_limit
    do i = 1, MSBRNM
       if(counter(i) > measure_count_limit) then
          ecpu0 = ecpu(i)
          ecpu(i) = ecpu(i)*counter(i)/measure_count_limit
          if(ipritiming >= 2) write(nfout,'(" id = ",i8, " ecpu new, old = ",2f8.4)') i, ecpu(i),ecpu0
       end if
    end do
    if(sw_details == ON) then
       do i = 1, n_allocated_subhir
          if(counter_ext(i) > measure_count_limit) then
             ecpu0 = ecpu_ext(i)
             ecpu_ext(i) = ecpu_ext(i)*counter_ext(i)/measure_count_limit
             if(ipritiming >= 2) write(nfout,'(" id = ",i8, " ecpu_ext new, old = ",2f8.4)') i, ecpu_ext(i),ecpu0
          end if
       end do
    end if
  end subroutine correct_ecpu

  subroutine m_Timing_initialize_tmp_timer()
     call gettod(tmp_timer_start_time)
  end subroutine m_Timing_initialize_tmp_timer

  function m_Timing_get_elapsed_time_tmp() result(res)
     real(kind=DP) :: res
     real(kind=DP) :: ct
     call gettod(ct)
     res = (ct-tmp_timer_start_time)*UMICRO
  end function m_Timing_get_elapsed_time_tmp

  subroutine m_Timing_finalize_tmp_timer()
     tmp_timer_start_time = 0.d0
  end subroutine m_Timing_finalize_tmp_timer

end module m_Timing
