!#define _DEBUG_WRITE_DFTU_MPI_PROCESSES_
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE:  MDiterationNumber_Setting, MDiterationNumber_Setting_ep,
!             IterationNumber_Setting, IterationNumber_Setting_g,
!             KpointNumber_Setting, KpointNumber_Setting2
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!      Further modification by T. Yamasaki   May 2004
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
! $Id: IterationNumbers_Setting.f90 614 2020-05-07 03:24:24Z jkoga $
subroutine MDiterationNumber_Setting
  use m_IterationNumbers,   only : m_Iter_mdIterN_increment&
       &, m_Iter_electronic_reset
  use m_Control_Parameters, only : m_CtrlP_reset_dtim_1Dsearch
  use m_Files,              only : nfout
  implicit none

  call m_CtrlP_reset_dtim_1Dsearch()
  call m_Iter_mdIterN_increment(nfout)
  call m_Iter_electronic_reset
end subroutine MDiterationNumber_Setting

subroutine SCDFTiterationNumber_Setting
  use m_IterationNumbers,   only : m_Iter_scdftIterN_increment &
       &, m_Iter_electronic_reset, m_Iter_reset_iter_ionic
  use m_Control_Parameters, only : m_CtrlP_reset_dtim_1Dsearch, m_CtrlP_reset_iconvergence
  use m_Files,              only : nfout
  implicit none

  call m_CtrlP_reset_dtim_1Dsearch()
  call m_Iter_scdftIterN_increment(nfout)
  call m_Iter_electronic_reset
  call m_Iter_reset_iter_ionic()
  call m_CtrlP_reset_iconvergence()
end subroutine SCDFTiterationNumber_Setting

subroutine MDiterationNumber_Setting2
  use m_IterationNumbers,   only : m_Iter_mdIterN_increment&
       &, m_Iter_electronic_reset,iteration_electronic,m_Iter_reset_iter_ionic &
       &, m_Iter_unitcell_increment &
       &, m_Iter_stress_correction_incre
  use m_Control_Parameters, only : m_CtrlP_reset_dtim_1Dsearch, sw_optimize_lattice &
       &, m_CtrlP_reset_iconvergence,imdalg, sw_stress_correction
  use m_Const_Parameters, only : ON, PT_CONTROL, P_CONTROL
  use m_Stress, only : m_Stress_in_correction, m_Stress_correction
  use m_Files, only : nfout
  implicit none
  logical :: Rightafter_stress_correction
  if(sw_optimize_lattice==ON .or. imdalg == PT_CONTROL .or. imdalg == P_CONTROL .or. m_Stress_in_correction(2))then
     call m_Iter_electronic_reset
     if(sw_optimize_lattice == ON) call m_Iter_reset_iter_ionic()
     if(.not.m_Stress_in_correction(2)) call m_Iter_unitcell_increment()
     call m_CtrlP_reset_iconvergence()
  endif
end subroutine MDiterationNumber_Setting2

logical function Rightafter_stress_correction()
  use m_IterationNumbers, only  : iteration_stress_correction
  implicit none
  Rightafter_stress_correction = iteration_stress_correction == 4
end function

subroutine MDiterationNumber_Setting_pre
  use m_IterationNumbers,   only : m_Iter_electronic_reset, m_Iter_stress_correction_incre
  use m_Control_Parameters, only : m_CtrlP_reset_iconvergence, sw_stress_correction
  use m_Const_Parameters, only : ON
  use m_Stress, only : m_Stress_in_correction, m_Stress_correction
  use m_Files, only : nfout
  implicit none
  if(sw_stress_correction == ON .and. m_Stress_in_correction(4)) then
     call m_Stress_correction(nfout)
     call m_Iter_electronic_reset
     call m_Iter_stress_correction_incre()
     call m_CtrlP_reset_iconvergence()
  endif
end subroutine MDiterationNumber_Setting_pre

subroutine MDiterationNumber_Setting3
  use m_Const_Parameters,   only : DRIVER_URAMP
  use m_Control_Parameters, only : driver
  use m_IterationNumbers,   only : m_Iter_electronic_reset,m_Iter_reset_iter_ionic &
       &, m_Iter_uramp_increment
  use m_Control_Parameters, only : driver, m_CtrlP_reset_iconvergence
  implicit none
  if(driver == DRIVER_URAMP)then
     call m_Iter_electronic_reset
     call m_Iter_reset_iter_ionic()
     call m_Iter_uramp_increment()
     call m_CtrlP_reset_iconvergence()
  endif
end subroutine MDiterationNumber_Setting3

subroutine MDiterationNumber_Setting_ep
  use m_IterationNumbers,   only : m_Iter_mdIterN_increment&
       &                         , m_Iter_total_increment
  use m_Files,              only : nfout
  implicit none
  call m_Iter_mdIterN_increment(nfout)
  call m_Iter_total_increment()
end subroutine MDiterationNumber_Setting_ep

subroutine IterationNumber_Setting
  use m_Control_Parameters,only: icond, ipritiming0&
       &                       , iprijobstatus, jobstatus_series, jobstatus_format, driver
  use m_IterationNumbers, only : m_Iter_electronic_incre &
       &                       , m_Iter_total_increment &
       &                       , iteration, first_iteration_of_this_job &
       &                       , iteration_electronic, iteration_ionic &
       &                       , iteration_unit_cell, iteration_uramp, iteration_scdft
  use m_Timing,           only : tstatc_wd, tstatc_wd0, tstatc_init, tstatc_iter &
       &                       , m_Timing_wd_status
  use m_Const_Parameters, only : INITIAL, CONTINUATION &
       &                       , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION &
       &                       , START, ITERATIVE, ON, OFF, DRIVER_URAMP
  use m_Files,            only : nfstatus,nfout &
       &                       , m_Files_open_nfstatus &
       &                       , m_Files_close_nfstatus &
       &                       , m_Files_skiptoend
  implicit none
  integer :: it, status_wdmode
  logical :: unitcell_can_change,Uramping,isSCDFT

  call tstatc_iter(iteration, first_iteration_of_this_job)
  it = iteration
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) it = iteration_electronic
  if(iteration == first_iteration_of_this_job) then
     call tstatc_wd0
     call flush(nfout)
  else
     if(ipritiming0 >= 1) call tstatc_wd(it)
  end if
  call tstatc_init

  if(iprijobstatus >=1 ) then
     call m_Files_open_nfstatus()
     if(iteration == first_iteration_of_this_job) then
        status_wdmode = START
     else
        status_wdmode = ITERATIVE
        if(jobstatus_series == ON) then
           call m_Files_skiptoend(nfstatus)
        else
        end if
     end if
     call m_Timing_wd_status(nfstatus,jobstatus_format,jobstatus_series,status_wdmode &
          & ,iteration,iteration_ionic,iteration_electronic)
     call m_Files_close_nfstatus()
  end if

  call m_Iter_electronic_incre()
  call m_Iter_total_increment()

#ifndef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
  if(iprijobstatus >=1 ) then
#endif
     if(.not.unitcell_can_change()) then
       if(Uramping()) then
          write(nfout,'(" ---- iteration(total, uramp, ionic, electronic) = ",4i8," ----")') &
          & iteration, iteration_uramp, iteration_ionic, iteration_electronic
       else if (isSCDFT()) then
          write(nfout,'(" ---- iteration(total, scdft, ionic, electronic) = ",4i8," ----")') &
          & iteration, iteration_scdft, iteration_ionic, iteration_electronic
       else
          write(nfout,'(" ---- iteration(total, ionic, electronic) = ",3i8," ----")') &
          & iteration, iteration_ionic, iteration_electronic
       endif
     else
       write(nfout,'(" ---- iteration(total, unitcell, ionic, electronic) = ",4i8," ----")') &
          & iteration, iteration_unit_cell,iteration_ionic, iteration_electronic
     endif
#ifndef _DEBUG_WRITE_DFTU_MPI_PROCESSES_
  end if
#endif

end subroutine IterationNumber_Setting

subroutine IterationNumber_Setting_g()
  use m_Control_Parameters,only: icond,ipritiming0
  use m_IterationNumbers, only : iteration, first_iteration_of_this_job &
       &                       , iteration_electronic &
       &                       , m_Iter_total_increment &
       &                       , m_Iter_electronic_incre
  use m_Timing          , only : tstatc_wd, tstatc_wd0, tstatc_init, tstatc_iter
  use m_Const_Parameters, only : INITIAL, CONTINUATION, FIXED_CHARGE, FIXED_CHARGE_CONTINUATION
  implicit none
  integer :: it

  call tstatc_iter(iteration, first_iteration_of_this_job)
  it = iteration
  if(icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION) it = iteration_electronic
  if(iteration_electronic == 0) then
     call tstatc_wd0()
  else
     if(ipritiming0 >= 1) call tstatc_wd(it)
  end if
  call tstatc_init()

  call m_Iter_electronic_incre()
  call m_Iter_total_increment()

end subroutine IterationNumber_Setting_g

subroutine IterationNumber_reset()
  use m_Timing          , only : tstatc_init
  use m_IterationNumbers, only : m_Iter_electronic_reset

  call m_Iter_electronic_reset()
  call tstatc_init()
end subroutine IterationNumber_reset

subroutine KpointNumber_Setting()
  use m_IterationNumbers, only : nk_in_the_process, first_kpoint_in_this_job &
       &                       , iteration_electronic, first_iteration_electronic &
       &                       , m_Iter_nk_incre &
       &                       , m_Iter_wd_electronic &
       &                       , m_Iter_wd_nk &
       &                       , m_Iter_electronic_reset &
       &                       , m_Iter_electronic_set
  use m_Control_Parameters,only: nspin, ipriekzaj, m_CtrlP_ntcnvg_reset
  use m_Files,only             : nfout

! ===================== added by K. Tagami ================== 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor

  implicit none
! =========================================================== 11.0

  call m_CtrlP_ntcnvg_reset()
! =================================== modified by K. Tagami ======== 11.0
!!  call m_Iter_nk_incre(nspin)
  if ( noncol ) then
    call m_Iter_nk_incre(ndim_spinor)
  else
    call m_Iter_nk_incre(nspin)
  endif
! ================================================================== 11.0

  if(nk_in_the_process == first_kpoint_in_this_job) then
     call m_Iter_electronic_set()
  else
     call m_Iter_electronic_reset()
  end if
  if(ipriekzaj <= 0) call m_Iter_electronic_reset()
  call m_Iter_wd_nk(nfout)
  call m_Iter_wd_electronic(nfout)
  
end subroutine KpointNumber_Setting

subroutine KpointNumber_Setting2()
  use m_IterationNumbers, only : nk_in_the_process, first_kpoint_in_this_job &
       &                       , iteration_electronic, first_iteration_electronic &
       &                       , m_Iter_nk_incre2 &
       &                       , m_Iter_wd_electronic &
       &                       , m_Iter_wd_nk2 &
       &                       , m_Iter_electronic_reset &
       &                       , m_Iter_electronic_set
!!$       &                       , m_Iter_nkgroup_set
  use m_Control_Parameters,only: nspin, ipriekzaj, m_CtrlP_ntcnvg_reset
  use m_Kpoints, only          : kv3, kv3_ek
  use m_Files,only             : nfout

! ===================== added by K. Tagami ================== 11.0
  use m_Control_Parameters,   only : noncol, ndim_spinor

  implicit none
! =========================================================== 11.0

  call m_CtrlP_ntcnvg_reset()

! ======================================= added by K. Tagami ============ 11.0
!  call m_Iter_nk_incre2(nspin,kv3_ek)
  if ( noncol ) then
    call m_Iter_nk_incre2( ndim_spinor,kv3_ek )
  else
    call m_Iter_nk_incre2( nspin,kv3_ek )
  endif
! ====================================================================== 11.0

  if(nk_in_the_process == first_kpoint_in_this_job) then
     call m_Iter_electronic_set()
!!$     call m_Iter_nkgroup_set(kv3)
  else
     call m_Iter_electronic_reset()
  end if
  if(ipriekzaj <= 0) call m_Iter_electronic_reset()
  call m_Iter_wd_nk2(nfout,kv3)
  call m_Iter_wd_electronic(nfout)
  
end subroutine KpointNumber_Setting2

!!$subroutine pIterationNumber_Setting()
!!$  use m_IterationNumbers, only : m_Iter_positron_set()
!!$  implicit none
!!$  call m_Iter_positron_set()
!!$end subroutine pIterationNumber_Setting
