#ifdef NEC_TIMER
#  define START_TIMER(a) call start_timer(a)
#  define STOP_TIMER(a)  call stop_timer(a)
#else
#  define START_TIMER(a)
#  define STOP_TIMER(a)
#endif
!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Forces
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
#ifdef LIBRARY_BUILD
subroutine Forces_phase0
#else
subroutine Forces
#endif
! $Id: Forces.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters,  only : WITHOUTTAG, WITHTAG, ON
  use m_Parallelization,   only : mype
  use m_Control_Parameters,only : ipriforce, sw_calc_force, sw_dipole_correction &
       &                        , sw_fef, sw_hybrid_functional
  use m_Files,             only : nfout, nfdynm
  use m_IterationNumbers,  only : iteration_ionic, iteration
  use m_Kpoints,           only : kv3
  use m_Ionic_System, only : pos, napt, fxyzew_l, fxyzvdw_l &
       &                   , m_IS_rigid_body_exists &
       &                   , m_IS_rigid_body_map_force
  use m_Force,        only : forc_l, fexx_l &
!!$       &                   , m_Force_alloc_zfsin_zfcos &
       &                   , m_Force_initialize &
       &                   , m_Force_sumup_and_symmetrize &
!!$       &                   , m_Force_dealloc_zfsin_zfcos &
       &                   , m_Force_cal_forcmx &
       &                   , m_Force_term_dipole &
       &                   , m_Force_term_fef &
!fj$$#ifndef PARA3D
       &                   , m_Force_term_Elocal_and_Epc &
       &                   , m_Force_term_drv_of_VlhxcQ &
       &                   , m_Force_term_drv_of_flmt
!fj$$#endif
  use m_Phonon,       only : m_Phonon_write_forces
#ifdef __TIMER__
  use m_Parallelization,     only : MPI_CommGroup
#endif
#ifdef NEC_TIMER
  use nec_timer
#endif
  use m_ES_ExactExchange, only : m_ES_EXX_Force

! ================================ added by K. Tagami ============= 11.0
  use m_Control_Parameters,      only : noncol
  use m_Force,                   only : m_Force_term_drv_of_flmt_noncl
! ================================================================ 11.0


  implicit none
  integer, save :: flag_wd_force = 0

#ifdef __TIMER__
  integer :: ierr
  call mpi_barrier(MPI_CommGroup, ierr)
  call timer_sta(19)
#endif
  if(ipriforce >= 2) write(nfout,'(" -- Forces -- ")')
!!$  call m_Force_alloc_zfsin_zfcos()
  call m_Force_initialize()
!fj$$#ifndef PARA3D
START_TIMER('Force')

! ==================== modified by K. Tagami ============== 11.0
!  call m_Force_term_drv_of_flmt(kv3,pos,napt)      ! (vnlsum) ->fnlxyz_l
!
    if ( noncol ) then
      call m_Force_term_drv_of_flmt_noncl(kv3,pos,napt)
    else
      call m_Force_term_drv_of_flmt(kv3,pos,napt)      ! (vnlsum) ->fnlxyz_l
    endif
! =========================================================== 11.0

STOP_TIMER('Force')
  call m_Force_term_drv_of_VlhxcQ(nfout,pos,napt)  ! (lclchg) ->flhxcq_l
  call m_Force_term_Elocal_and_Epc(nfout,pos)      ! (force)  ->forc_l
!fj$$#endif
  if(sw_dipole_correction == ON) then
     call m_Force_term_dipole()                       ! -> fdip_l,fext_l
  end if
  if(sw_fef==ON) call m_Force_term_fef()           ! -> ffef_l
  if(sw_hybrid_functional==ON) call m_ES_EXX_Force(fexx_l)
  call m_Force_sumup_and_symmetrize(nfout,fxyzew_l,napt)
  !                         forc_l,fnlxyz_l,flhxcq_l,fxyzew_l  --> forc_l
!!$  call m_Force_dealloc_zfsin_zfcos()
  call m_Force_cal_forcmx()                        ! -> forcmx
  if(sw_calc_force == ON) call m_Phonon_write_forces()

  if(m_IS_rigid_body_exists()) then
    call m_IS_rigid_body_map_force(forc_l)
  endif
!!$  call m_IS_evaluate_v_verlet(forc_l)
!!$

#ifdef __TIMER__
  call timer_end(19)
#endif
#ifdef LIBRARY_BUILD
end subroutine Forces_phase0
#else
end subroutine Forces
#endif
