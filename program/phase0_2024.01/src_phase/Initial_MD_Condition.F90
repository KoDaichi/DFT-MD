!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Initial_MD_Condition
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
subroutine Initial_MD_Condition
! $Id: Initial_MD_Condition.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters,   only : CONTINUATION, T_CONTROL, BLUEMOON, QUENCHED_CONSTRAINT &
 &                               , COORDINATE_CONTINUATION, INITIAL, FIXED_CHARGE_CONTINUATION
  use m_Files,              only : nfinp, nfcntn, nfout
#ifdef _EMPIRICAL_
  use m_Control_Parameters, only : icond, imdalg
  use m_IterationNumbers,   only : m_Iter_total_increment
#else
  use m_Control_Parameters, only : icond, imdalg, iprimd &
       &                         , m_CtrlP_rd_isolver
#endif
  use m_Total_Energy,       only : m_TE_rd_total_energy
  use m_Ionic_System,       only : m_IS_rd_forcp_etc &
       &                         , m_IS_rd_nrsv, m_IS_natm_can_change
!!$  use m_Ionic_System,       only : m_IS_rd_T_parameters, m_IS_rd_forcp_etc &
!!$       &                         , m_IS_rd_nrsv, m_IS_rd_nrsv_stdin

  implicit none
  if(icond == CONTINUATION.or.icond==COORDINATE_CONTINUATION.or.icond == FIXED_CHARGE_CONTINUATION) then
     call m_TE_rd_total_energy(nfcntn)
     if(iprimd >= 3 ) write(nfout,'(" imdalg = ",i8)') imdalg
     if(imdalg==T_CONTROL.or.imdalg==BLUEMOON.or.imdalg==QUENCHED_CONSTRAINT) then
        if(imdalg /= QUENCHED_CONSTRAINT) call m_IS_rd_nrsv(nfcntn)
        call m_IS_rd_forcp_etc(imdalg,nfcntn)
     end if
#ifndef _EMPIRICAL_
     call m_CtrlP_rd_isolver(nfcntn)
#endif
     if(icond==CONTINUATION.and.m_IS_natm_can_change())then
        icond = INITIAL
     endif
  end if

#ifdef _EMPIRICAL_
  call m_Iter_total_increment()
#endif

end subroutine Initial_MD_Condition
