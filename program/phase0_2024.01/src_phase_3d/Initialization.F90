!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 622 $)
!
!  SUBROUINE: Initialization_set_ekmode_ON, Initialization, aavers
!
!  AUTHORS: T. Yamasaki, K. Betsuyaku,   August/20/2003
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
subroutine Initialization_set_ekmode_ON
! $Id: Initialization.F90 622 2020-05-26 05:22:07Z jkoga $
  use m_Control_Parameters, only : m_CtrlP_set_ekmode_ON
  implicit none

  call m_CtrlP_set_ekmode_ON()
end subroutine Initialization_set_ekmode_ON

subroutine Initialization_set_ekmode_GRID
  use m_Control_Parameters, only : m_CtrlP_set_ekmode_GRID
  implicit none

  call m_CtrlP_set_ekmode_GRID()
end subroutine Initialization_set_ekmode_GRID

logical function initialization_required()
  use m_Const_Parameters, only : OFF,ON,CONTINUATION,FIXED_CHARGE_CONTINUATION,DRIVER_URAMP,DRIVER_SC_DFT,DRIVER_NEB
  use m_Control_Parameters, only : sw_optimize_lattice,sw_rebuild_pws,icond &
   & , m_CtrlP_in_initialization, driver
  use m_IterationNumbers, only : iteration_unit_cell, iteration_uramp, iteration_scdft
  use m_Ionic_System, only : m_IS_natm_can_change
  implicit none
  logical :: unitcell_can_change
  if((icond==CONTINUATION.or.icond==FIXED_CHARGE_CONTINUATION).and.m_CtrlP_in_initialization())then
     initialization_required = .true.
     return
  endif
  if( driver == DRIVER_URAMP .and. iteration_uramp > 1)then
     initialization_required = .false.
    return
  endif
  if( driver == DRIVER_SC_DFT .and. iteration_scdft > 1)then
     initialization_required = .false.
    return
  endif
  if (.not.unitcell_can_change().or.iteration_unit_cell==1) then
    initialization_required = .true.
    return
  endif
  if(m_IS_natm_can_change())then
    initialization_required = .true.
    return
  endif
  if(sw_rebuild_pws==ON) then
     initialization_required = .true.
     return
  endif

  initialization_required = .false.
end function initialization_required

subroutine Initialization_lib(mpi_comm_root, tag, ne, nk, ng, stdout)
! $Id: Initialization.F90 622 2020-05-26 05:22:07Z jkoga $
  use m_Parallelization, only : m_Parallel_init_comm_group, npes
  use m_Parallelization, only : nrank_e, nrank_k, nrank_g, read_from_args
  use m_Parallelization, only : m_Parallel_get_nproc_from_arg_3D

#ifdef _USE_SCALAPACK_
  use m_Control_Parameters,only:  m_CtrlP_set_sw_scalapack
#endif
  use m_Timing,          only : m_Timing_wd_timenow, m_Timing_init_timer
  use m_Files,           only : nfout &
       &                      , m_Files_open_standardout &
       &                      , m_Files_set_default_filenames &
       &                      , m_Files_rd_file_names_data &
       &                      , m_Files_open_files_initially &
       &                      , m_Files_check_file_existence &
       &                      , m_Files_set_conftag &
       &                      , m_Files_set_stdout
  use m_Ionic_System,    only : m_IS_initialize_mdmode
  use m_Control_Parameters,only: m_CtrlP_set_printable, printable &
       &                      , m_CtrlP_set_wct_start, m_CtrlP_set_sw_scalapack &
       &                      , m_CtrlP_set_conftag
  implicit none

  integer, intent(in) :: mpi_comm_root, tag, ne, nk, ng, stdout
  logical,save :: first_call=.true.
  logical :: initialization_required
  integer :: ierr

  if(.not.initialization_required()) return
!    -----------
  if(first_call)then
#ifdef _USE_SCALAPACK_
!!  call m_CtrlP_set_sw_scalapack(printable,nfout)
#endif
  call m_CtrlP_set_conftag(tag)
  call m_Files_set_conftag(tag)
  call m_Files_set_stdout(stdout)
  call m_Parallel_init_comm_group(mpi_comm_root)
  call m_CtrlP_set_printable()             !-> printable
  call m_Timing_init_timer()
  call m_Files_set_default_filenames()
  call m_Files_rd_file_names_data()
  call m_Files_open_standardout()
  if(printable) call print_title()
  if(printable) call m_Timing_wd_timenow("program start")
  if ((ne.ge.1) .and. (nk.ge.1) .and. (ng.ge.1) .and. (npes .eq. (ne*nk*ng))) then
    call m_Parallel_get_nproc_from_arg_3D(printable,nfout,ne,nk,ng)
    if(printable) write(nfout,'(a, 4i8)') '!** np ne nk ng ',npes,nrank_e,nrank_k,nrank_g
  else
    call m_Parallel_get_nproc_from_arg_3D(printable,nfout)
  endif
#ifdef _USE_SCALAPACK_
!!  call m_CtrlP_set_sw_scalapack(printable, nfout)
#endif
  call aavers                      ! -(here)
  call m_Files_check_file_existence
  endif
  call m_Files_open_files_initially
  call m_IS_initialize_mdmode()
  if(first_call) call m_CtrlP_set_wct_start                   !  (ckcput)
!  call m_CtrlP_set_wct_start                   !  (ckcput)
!!$  call m_CtrlP_set_paramset_off
  first_call = .false.

contains
  subroutine print_title()
    integer :: i
    write(nfout,'(1x,79a)') ("*",i=1,79)
!       10        20                  40        50        60        70        80        90       100       110       120       130      140
!2345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    !                        123456789a      123456789a      123456789a      123456789a      123456789a      1234567      12345678
!!$    write(nfout,'(1x,"*",3x,"PPPPPPP   ",1x,"HH     HH ",1x,"    AA    ",1x,"   SSSSS  ",1x," EEEEEEEE ",1x,"       ",1x,"        ",3x,"*")') ! a
!!$    write(nfout,'(1x,"*",3x,"PPPPPPPP  ",1x,"HH     HH ",1x,"   AAAA   ",1x,"  SSSSSSS ",1x," EEEEEEEE ",1x,"      /",1x,"        ",3x,"*")') ! b
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"   AAAA   ",1x,"  SS   SS ",1x," EE       ",1x,"     //",1x,"  0000  ",3x,"*")') ! c
!!$!   write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SS       ",1x," EE       ",1x,"     / ",1x," 000000 ",3x,"*")') ! d
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SS       ",1x," EE       ",1x,"    // ",1x,"000  000",3x,"*")') ! e
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SSS      ",1x," EE       ",1x,"    /  ",1x,"00    00",3x,"*")') ! f
!!$    write(nfout,'(1x,"*",3x,"PPPPPPPP  ",1x,"HHHHHHHHH ",1x,"  AA  AA  ",1x,"  SSSSS   ",1x," EEEEEEEE ",1x,"   //  ",1x,"00    00",3x,"*")') ! g
!!$    write(nfout,'(1x,"*",3x,"PPPPPPP   ",1x,"HHHHHHHHH ",1x," AA    AA ",1x,"    SSSS  ",1x," EEEEEEEE ",1x,"   /   ",1x,"00    00",3x,"*")') ! h
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AA    AA ",1x,"      SSS ",1x," EE       ",1x,"   /   ",1x,"00    00",3x,"*")') ! i
!!$!   write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AAAAAAAA ",1x,"       SS ",1x," EE       ",1x,"  //   ",1x,"00    00",3x,"*")') ! j
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AAAAAAAA ",1x,"       SS ",1x," EE       ",1x,"  /    ",1x,"00    00",3x,"*")') ! k
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AAAAAAAAAA",1x," SS   SSS ",1x," EE       ",1x," //    ",1x,"00    00",3x,"*")') ! l
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AA      AA",1x," SSSSSSS  ",1x," EEEEEEEE ",1x," /     ",1x," 000000 ",3x,"*")') ! m
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AA      AA",1x,"  SSSSS   ",1x," EEEEEEEE ",1x,"//     ",1x,"  0000  ",3x,"*")') ! n
!!$    write(nfout,'(1x,79a)') ("*",i=1,79)
    write(nfout,2011) "PPPPPPP   ","HH     HH ","    AA    ","   SSSSS  "," EEEEEEEE ","       ","        "  ! a
    write(nfout,2011) "PPPPPPPP  ","HH     HH ","   AAAA   ","  SSSSSSS "," EEEEEEEE ","      /","        "  ! b
    write(nfout,2011) "PP     PP ","HH     HH ","   AAAA   ","  SS   SS "," EE       ","     //","  0000  "  ! c
!   write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SS       "," EE       ","     / "," 000000 "  ! d
    write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SS       "," EE       ","    // ","000  000"  ! e
    write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SSS      "," EE       ","    /  ","00    00"  ! f
    write(nfout,2011) "PPPPPPPP  ","HHHHHHHHH ","  AA  AA  ","  SSSSS   "," EEEEEEEE ","   //  ","00    00"  ! g
    write(nfout,2011) "PPPPPPP   ","HHHHHHHHH "," AA    AA ","    SSSS  "," EEEEEEEE ","   /   ","00    00"  ! h
    write(nfout,2011) "PP        ","HH     HH "," AA    AA ","      SSS "," EE       ","   /   ","00    00"  ! i
!   write(nfout,2011) "PP        ","HH     HH "," AAAAAAAA ","       SS "," EE       ","  //   ","00    00"  ! j
    write(nfout,2011) "PP        ","HH     HH "," AAAAAAAA ","       SS "," EE       ","  /    ","00    00"  ! k
    write(nfout,2011) "PP        ","HH     HH ","AAAAAAAAAA"," SS   SSS "," EE       "," //    ","00    00"  ! l
    write(nfout,2011) "PP        ","HH     HH ","AA      AA"," SSSSSSS  "," EEEEEEEE "," /     "," 000000 "  ! m
    write(nfout,2011) "PP        ","HH     HH ","AA      AA","  SSSSS   "," EEEEEEEE ","//     ","  0000  "  ! n
    write(nfout,'(1x,79a)') ("*",i=1,79)
2011 format(1x,'*',3x,a10,1x,a10,1x,a10,1x,a10,1x,a10,1x,a7,1x,a8,3x,'*')
!!$    write(nfout,'(1x,"****",10a,1x,10a,1x,10a,1x,10a,1x,10a,3x,"*")') ("**********",i=1,5)
  end subroutine print_title

  subroutine aavers
    include 'version.h' ! commit ID
    character(len=80) :: vers, system, codename
    !write(vers,'("Revision:",i5," -- 3D_Parallel --")') commit_id
!    write(vers,'("phase/0 2021.01 Revision:",i5," -- 3D_Parallel --")') commit_id
    vers = "phase/0 2024.01 Revision:"//commit_id//" -- 3D_Parallel --"
    codename = 'phaseUnif'
    system = ''

#ifdef VPP
    system = '@(#)system=vpp'
#elif DEC
    system = '@(#)system=dec'
#elif HP
    system = '@(#)system=hp'
#elif SUN
    system = '@(#)system=sun'
#elif ONYX
    system = '@(#)system=onyx'
#elif IRIX64
    system = '@(#)system=irix64'
#elif CRAY
    system = '@(#)system=crayxmp'
#elif HIUX
    system = '@(#)system=hi-ux'
#elif SX
    system = '@(#)system=sx'
#elif SP2
    system = '@(#)system=aix'
#elif Linux
#ifdef PGI
    system = '@(#)system=linux_pgi'
#else
    system = '@(#)system=linux'
#endif
#endif
    if(printable) then
!!$       write(nfout,*) ' ********************************'
!!$       write(nfout,*) ' *  PHASE/0 ver.2013.11         *'
!!$       write(nfout,*) ' *  (phaseUnif_noncol_r340rev+) *'
!!$       write(nfout,*) ' ********************************'
       write(nfout,'(a80)') trim(vers)
       write(nfout,*) trim(system)
       write(nfout,*) trim(codename)
       write(nfout,'(" --- << CPP options defined in the makefile >> --")')
#ifdef _NO_MPI_
#ifdef _OPENMP
       write(nfout,'(" Parallization: OpenMP")')
#else
       write(nfout,'(" Parallization: Serial")')
#endif
#else
#ifdef _OPENMP
       write(nfout,'(" Parallization: OpenMP and MPI")')
#else
       write(nfout,'(" Parallization: MPI")')
#endif
#endif
    end if
    if(printable) then
#ifdef TRANSPOSE
       write(nfout,'(" MGS  = TRANSPOSE")')
#elif CYCLIC
       write(nfout,'(" MGS  = CYCLIC")')
#endif
#ifdef _EMPIRICAL_
       write(nfout,'(" _EMPIRICAL_")')
#endif
!!$#ifdef _SIMPLE_SORT_
!!$       write(nfout,'(" SORTING = _SIMPLE_SORT_")')
!!$#elif  _HEAP_SORT_
!!$       write(nfout,'(" SORTING = _HEAP_SORT_")')
!!$#endif
#ifdef DECFFT
       write(nfout,'(" FFT WF  = DECFFT")')
#elif  MKLFFT
       write(nfout,'(" FFT WF  = MKLFFT")')
#elif  FFTW3
       write(nfout,'(" FFT WF  = FFTW3")')
#elif  ACMLFFT
       write(nfout,'(" FFT WF  = ACMLFFT")')
#elif  SCSLFFT
       write(nfout,'(" FFT WF  = SCSLFFT")')
#elif  WF_SRFFT
       write(nfout,'(" FFT WF  = WF_SRFFT")')
#elif  WF_JRCATFFT
       write(nfout,'(" FFT WF  = WF_JRCATFFT")')
#elif  WF_JRCATFFT_WS
       write(nfout,'(" FFT WF  = WF_JRCATFFT_WS")')
#endif
    end if
#ifdef CD_SRFFT
    if(printable) write(nfout,'(" FFT CD  = CD_SRFFT")')
#elif CD_JRCATFFT
    if(printable) write(nfout,'(" FFT CD  = CD_JRCATFFT")')
#ifdef WF_JRCATFFT_WS
    if(printable) write(nfout,'(" Do NOT use WF_JRCATFFT_WS and CD_JRCATFFT at the same time")')
    call phase_error_with_msg(nfout,' Do NOT use WF_JRCATFFT_WS and CD_JRCATFFT at the same time',__LINE__,__FILE__)
#endif
#elif  CD_JRCATFFT_WS
    if(printable) write(nfout,'(" FFT CD  = CD_JRCATFFT_WS")')
#ifdef WF_JRCATFFT
    if(printable) write(nfout,'(" Do NOT use WF_JRCATFFT and CD_JRCATFFT_WS at the same time")')
    call phase_error_with_msg(nfout,'Do NOT use WF_JRCATFFT and CD_JRCATFFT_WS at the same time',__LINE__,__FILE__)
#endif
#elif  FFTW3
    if(printable) write(nfout,'(" FFT CD  = FFTW3")')
#endif
#ifdef NO_MGS_DGEMM
    if(printable) write(nfout,'(" NO_MGS_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_MGS_DGEMM is not defined")')
#endif
#ifdef NO_NONLOCAL_DGEMM
    if(printable) write(nfout,'(" NO_NONLOCAL_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_NONLOCAL_DGEMM is not defined")')
#endif
#ifdef NO_NONLOCAL_RMM_DGEMM
    if(printable) write(nfout,'(" NO_NONLOCAL_RMM_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_NONLOCAL_RMM_DGEMM is not defined")')
#endif
#ifdef NO_SUBMAT_DGEMM
    if(printable) write(nfout,'(" NO_SUBMAT_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_SUBMAT_DGEMM is not defined")')
#endif
#ifdef NO_FORCE_DGEMM
    if(printable) write(nfout,'(" NO_FORCE_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_FORCE_DGEMM is not defined")')
#endif
#ifdef NO_MATDIAGON_DGEMM
    if(printable) write(nfout,'(" NO_MATDIAGON_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_MATDIAGON_DGEMM is not defined")')
#endif


#ifdef LMM_PREVIOUS
    if(printable) write(nfout,'(" LMM_PREVIOUS is defined")')
#else
    if(printable) write(nfout,'(" LMM_PREVIOUS is not defined")')
#endif
    if(printable) write(nfout,'(" ----------------------------------------------")')
  end subroutine aavers
  
end subroutine Initialization_lib

subroutine Initialization(init_mpi)
! $Id: Initialization.F90 622 2020-05-26 05:22:07Z jkoga $
  use m_Parallelization, only : m_Parallel_init_comm_world
  use m_Parallelization, only : m_Parallel_get_nproc_from_arg_3D

#ifdef _USE_SCALAPACK_
  use m_Control_Parameters,only:  m_CtrlP_set_sw_scalapack
#endif
  use m_Timing,          only : m_Timing_wd_timenow, m_Timing_init_timer
  use m_Files,           only : nfout &
       &                      , m_Files_open_standardout &
       &                      , m_Files_set_default_filenames &
       &                      , m_Files_rd_file_names_data &
       &                      , m_Files_open_files_initially &
       &                      , m_Files_check_file_existence
  use m_Ionic_System,    only : m_IS_initialize_mdmode
  use m_Control_Parameters,only: m_CtrlP_set_printable, printable &
       &                      , m_CtrlP_set_wct_start
  implicit none

  integer, optional, intent(in) :: init_mpi
  logical,save :: first_call=.true.
  logical :: initialization_required

  if(.not.initialization_required()) return

!    -----------
  if(first_call)then
  call m_Parallel_init_comm_world(init_mpi)
  call m_CtrlP_set_printable()             !-> printable
  call m_Timing_init_timer()
  call m_Files_set_default_filenames()
  call m_Files_rd_file_names_data()
  call m_Files_open_standardout()
  if(printable) call print_title()
  if(printable) call m_Timing_wd_timenow("program start")
  call m_Parallel_get_nproc_from_arg_3D(printable,6)
#ifdef _USE_SCALAPACK_
!!  call m_CtrlP_set_sw_scalapack(printable, nfout)
#endif
  call aavers                      ! -(here)
  call m_Files_check_file_existence
  endif
  call m_Files_open_files_initially
  call m_IS_initialize_mdmode()
  if(first_call) call m_CtrlP_set_wct_start                   !  (ckcput)
!  call m_CtrlP_set_wct_start                   !  (ckcput)
!!$  call m_CtrlP_set_paramset_off
  first_call = .false.
contains
  subroutine print_title()
    integer :: i
    write(nfout,'(1x,79a)') ("*",i=1,79)
!       10        20                  40        50        60        70        80        90       100       110       120       130      140
!2345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    !                        123456789a      123456789a      123456789a      123456789a      123456789a      1234567      12345678
!!$    write(nfout,'(1x,"*",3x,"PPPPPPP   ",1x,"HH     HH ",1x,"    AA    ",1x,"   SSSSS  ",1x," EEEEEEEE ",1x,"       ",1x,"        ",3x,"*")') ! a
!!$    write(nfout,'(1x,"*",3x,"PPPPPPPP  ",1x,"HH     HH ",1x,"   AAAA   ",1x,"  SSSSSSS ",1x," EEEEEEEE ",1x,"      /",1x,"        ",3x,"*")') ! b
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"   AAAA   ",1x,"  SS   SS ",1x," EE       ",1x,"     //",1x,"  0000  ",3x,"*")') ! c
!!$!   write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SS       ",1x," EE       ",1x,"     / ",1x," 000000 ",3x,"*")') ! d
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SS       ",1x," EE       ",1x,"    // ",1x,"000  000",3x,"*")') ! e
!!$    write(nfout,'(1x,"*",3x,"PP     PP ",1x,"HH     HH ",1x,"  AA  AA  ",1x," SSS      ",1x," EE       ",1x,"    /  ",1x,"00    00",3x,"*")') ! f
!!$    write(nfout,'(1x,"*",3x,"PPPPPPPP  ",1x,"HHHHHHHHH ",1x,"  AA  AA  ",1x,"  SSSSS   ",1x," EEEEEEEE ",1x,"   //  ",1x,"00    00",3x,"*")') ! g
!!$    write(nfout,'(1x,"*",3x,"PPPPPPP   ",1x,"HHHHHHHHH ",1x," AA    AA ",1x,"    SSSS  ",1x," EEEEEEEE ",1x,"   /   ",1x,"00    00",3x,"*")') ! h
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AA    AA ",1x,"      SSS ",1x," EE       ",1x,"   /   ",1x,"00    00",3x,"*")') ! i
!!$!   write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AAAAAAAA ",1x,"       SS ",1x," EE       ",1x,"  //   ",1x,"00    00",3x,"*")') ! j
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x," AAAAAAAA ",1x,"       SS ",1x," EE       ",1x,"  /    ",1x,"00    00",3x,"*")') ! k
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AAAAAAAAAA",1x," SS   SSS ",1x," EE       ",1x," //    ",1x,"00    00",3x,"*")') ! l
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AA      AA",1x," SSSSSSS  ",1x," EEEEEEEE ",1x," /     ",1x," 000000 ",3x,"*")') ! m
!!$    write(nfout,'(1x,"*",3x,"PP        ",1x,"HH     HH ",1x,"AA      AA",1x,"  SSSSS   ",1x," EEEEEEEE ",1x,"//     ",1x,"  0000  ",3x,"*")') ! n
!!$    write(nfout,'(1x,79a)') ("*",i=1,79)
    write(nfout,2011) "PPPPPPP   ","HH     HH ","    AA    ","   SSSSS  "," EEEEEEEE ","       ","        "  ! a
    write(nfout,2011) "PPPPPPPP  ","HH     HH ","   AAAA   ","  SSSSSSS "," EEEEEEEE ","      /","        "  ! b
    write(nfout,2011) "PP     PP ","HH     HH ","   AAAA   ","  SS   SS "," EE       ","     //","  0000  "  ! c
!   write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SS       "," EE       ","     / "," 000000 "  ! d
    write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SS       "," EE       ","    // ","000  000"  ! e
    write(nfout,2011) "PP     PP ","HH     HH ","  AA  AA  "," SSS      "," EE       ","    /  ","00    00"  ! f
    write(nfout,2011) "PPPPPPPP  ","HHHHHHHHH ","  AA  AA  ","  SSSSS   "," EEEEEEEE ","   //  ","00    00"  ! g
    write(nfout,2011) "PPPPPPP   ","HHHHHHHHH "," AA    AA ","    SSSS  "," EEEEEEEE ","   /   ","00    00"  ! h
    write(nfout,2011) "PP        ","HH     HH "," AA    AA ","      SSS "," EE       ","   /   ","00    00"  ! i
!   write(nfout,2011) "PP        ","HH     HH "," AAAAAAAA ","       SS "," EE       ","  //   ","00    00"  ! j
    write(nfout,2011) "PP        ","HH     HH "," AAAAAAAA ","       SS "," EE       ","  /    ","00    00"  ! k
    write(nfout,2011) "PP        ","HH     HH ","AAAAAAAAAA"," SS   SSS "," EE       "," //    ","00    00"  ! l
    write(nfout,2011) "PP        ","HH     HH ","AA      AA"," SSSSSSS  "," EEEEEEEE "," /     "," 000000 "  ! m
    write(nfout,2011) "PP        ","HH     HH ","AA      AA","  SSSSS   "," EEEEEEEE ","//     ","  0000  "  ! n
    write(nfout,'(1x,79a)') ("*",i=1,79)
2011 format(1x,'*',3x,a10,1x,a10,1x,a10,1x,a10,1x,a10,1x,a7,1x,a8,3x,'*')
!!$    write(nfout,'(1x,"****",10a,1x,10a,1x,10a,1x,10a,1x,10a,3x,"*")') ("**********",i=1,5)
  end subroutine print_title

  subroutine aavers
    include 'version.h' ! commit_id
    character(len=80) :: vers, system, codename
!    write(vers,'("Revision:",i5, " --- 3D_Parallel --")') commit_id
!    write(vers,'("phase/0 2021.02 Revision:",i5, " --- 3D_Parallel --")') commit_id
    vers = "phase/0 2023.01.01 Revision:"//commit_id//" -- 3D_Parallel --"
    codename = 'phaseUnif'
    system = ''

#ifdef VPP
    system = '@(#)system=vpp'
#elif DEC
    system = '@(#)system=dec'
#elif HP
    system = '@(#)system=hp'
#elif SUN
    system = '@(#)system=sun'
#elif ONYX
    system = '@(#)system=onyx'
#elif IRIX64
    system = '@(#)system=irix64'
#elif CRAY
    system = '@(#)system=crayxmp'
#elif HIUX
    system = '@(#)system=hi-ux'
#elif SX
    system = '@(#)system=sx'
#elif SP2
    system = '@(#)system=aix'
#elif Linux
#ifdef PGI
    system = '@(#)system=linux_pgi'
#else
    system = '@(#)system=linux'
#endif
#endif
    if(printable) then
!!$       write(nfout,*) ' ********************************'
!!$       write(nfout,*) ' *  PHASE/0 ver.2013.11         *'
!!$       write(nfout,*) ' *  (phaseUnif_noncol_r340rev+) *'
!!$       write(nfout,*) ' ********************************'
       write(nfout,'(a80)') trim(vers)
       write(nfout,*) trim(system)
       write(nfout,*) trim(codename)
       write(nfout,'(" --- << CPP options defined in the makefile >> --")')
#ifdef _NO_MPI_
#ifdef _OPENMP
       write(nfout,'(" Parallization: OpenMP")')
#else
       write(nfout,'(" Parallization: Serial")')
#endif
#else
#ifdef _OPENMP
       write(nfout,'(" Parallization: OpenMP and MPI")')
#else
       write(nfout,'(" Parallization: MPI")')
#endif
#endif
    end if
    if(printable) then
#ifdef TRANSPOSE
       write(nfout,'(" MGS  = TRANSPOSE")')
#elif CYCLIC
       write(nfout,'(" MGS  = CYCLIC")')
#endif
#ifdef _EMPIRICAL_
       write(nfout,'(" _EMPIRICAL_")')
#endif
#ifdef _SIMPLE_SORT_
       write(nfout,'(" SORTING = _SIMPLE_SORT_")')
#elif  _HEAP_SORT_
       write(nfout,'(" SORTING = _HEAP_SORT_")')
#endif
#ifdef DECFFT
       write(nfout,'(" FFT WF  = DECFFT")')
#elif  MKLFFT
       write(nfout,'(" FFT WF  = MKLFFT")')
#elif  FFTW3
       write(nfout,'(" FFT WF  = FFTW3")')
#elif  ACMLFFT
       write(nfout,'(" FFT WF  = ACMLFFT")')
#elif  SCSLFFT
       write(nfout,'(" FFT WF  = SCSLFFT")')
#elif  WF_SRFFT
       write(nfout,'(" FFT WF  = WF_SRFFT")')
#elif  WF_JRCATFFT
       write(nfout,'(" FFT WF  = WF_JRCATFFT")')
#elif  WF_JRCATFFT_WS
       write(nfout,'(" FFT WF  = WF_JRCATFFT_WS")')
#endif
    end if
#ifdef CD_SRFFT
    if(printable) write(nfout,'(" FFT CD  = CD_SRFFT")')
#elif CD_JRCATFFT
    if(printable) write(nfout,'(" FFT CD  = CD_JRCATFFT")')
#ifdef WF_JRCATFFT_WS
    if(printable) write(nfout,'(" Do NOT use WF_JRCATFFT_WS and CD_JRCATFFT at the same time")')
    call phase_error_with_msg(nfout,'Do NOT use WF_JRCATFFT_WS and CD_JRCATFFT at the same time',__LINE__,__FILE__)
#endif
#elif  CD_JRCATFFT_WS
    if(printable) write(nfout,'(" FFT CD  = CD_JRCATFFT_WS")')
#ifdef WF_JRCATFFT
    if(printable) write(nfout,'(" Do NOT use WF_JRCATFFT and CD_JRCATFFT_WS at the same time")')
    call phase_error_with_msg(nfout,'Do NOT use WF_JRCATFFT and CD_JRCATFFT_WS at the same time',__LINE__,__FILE__)
#endif
#elif  FFTW3
    if(printable) write(nfout,'(" FFT CD  = FFTW3")')
#endif
#ifdef NO_MGS_DGEMM
    if(printable) write(nfout,'(" NO_MGS_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_MGS_DGEMM is not defined")')
#endif
#ifdef NO_NONLOCAL_DGEMM
    if(printable) write(nfout,'(" NO_NONLOCAL_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_NONLOCAL_DGEMM is not defined")')
#endif
#ifdef NO_NONLOCAL_RMM_DGEMM
    if(printable) write(nfout,'(" NO_NONLOCAL_RMM_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_NONLOCAL_RMM_DGEMM is not defined")')
#endif
#ifdef NO_SUBMAT_DGEMM
    if(printable) write(nfout,'(" NO_SUBMAT_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_SUBMAT_DGEMM is not defined")')
#endif
#ifdef NO_FORCE_DGEMM
    if(printable) write(nfout,'(" NO_FORCE_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_FORCE_DGEMM is not defined")')
#endif
#ifdef NO_MATDIAGON_DGEMM
    if(printable) write(nfout,'(" NO_MATDIAGON_DGEMM is defined")')
#else
    if(printable) write(nfout,'(" NO_MATDIAGON_DGEMM is not defined")')
#endif


#ifdef LMM_PREVIOUS
    if(printable) write(nfout,'(" LMM_PREVIOUS is defined")')
#else
    if(printable) write(nfout,'(" LMM_PREVIOUS is not defined")')
#endif
    if(printable) write(nfout,'(" ----------------------------------------------")')
  end subroutine aavers
  
end subroutine Initialization
