! ================================= added by K. Tagami  ==================== 5.0
subroutine Renewal_of_Hubbard_Parameters()
  use m_Files,               only : nfout
  use m_Ionic_System,        only : natm
  use m_Const_Parameters,     only : OFF, ON, DP
  use m_Control_Parameters,      only :  sw_eval_Ueff_using_edelta, &
       &                                 sw_eval_Ueff_using_iter, &
       &                                 iteration_for_Ueff_starting, &
       &                                 edelta_for_Ueff_starting, &
       &                                 Ueff_transition_period
  use m_IterationNumbers,     only : iteration_electronic
  use m_Total_Energy,         only : m_TE_what_is_edeltb_now
  use m_Hubbard,             only : Ueff_prefactor

  use m_Control_Parameters,   only : printable

  implicit none

  real(kind=DP) :: edeltab_per_atom
!
! -------------------- 
  edeltab_per_atom = m_TE_what_is_edeltb_now() / natm
  edeltab_per_atom = abs(edeltab_per_atom )

  if ( sw_eval_Ueff_using_edelta==ON .or. sw_eval_Ueff_using_iter==ON ) then
    call set_new_Ueff_prefactor
  endif

contains

  subroutine set_new_Ueff_prefactor

    integer,save :: iteration0_save = 0
    integer, save :: Start_Ueff_Transition = OFF

    Ueff_prefactor = 1.0D0
    if ( sw_eval_Ueff_using_iter == ON ) then
       Ueff_prefactor = dble( iteration_electronic - iteration_for_Ueff_starting ) &
            &          / dble(Ueff_transition_period)
       
       if ( Ueff_prefactor < 0.0 ) Ueff_prefactor = 0.0D0
       if ( Ueff_prefactor > 1.0 ) Ueff_prefactor = 1.0D0
       
    else if ( sw_eval_Ueff_using_edelta == ON ) then
       
       if ( Start_Ueff_Transition == OFF ) then
          Ueff_prefactor = 0.0D0
          if ( edeltab_per_atom < edelta_for_Ueff_starting ) then
             iteration0_save = iteration_electronic
             Start_Ueff_Transition = ON
          endif
       endif
       
       if ( Start_Ueff_Transition == ON ) then
          Ueff_prefactor = dble( iteration_electronic - iteration0_save ) &
               &          / dble(Ueff_transition_period)
          
          if ( Ueff_prefactor < 0.0 ) Ueff_prefactor = 0.0D0
          if ( Ueff_prefactor > 1.0 ) Ueff_prefactor = 1.0D0
       endif
    endif
    if ( printable ) then
       write(nfout,*)  '!== * Ueff_prefactor_now = ', &
            &        iteration_electronic, Ueff_prefactor
    endif

  end subroutine set_new_Ueff_prefactor

end subroutine Renewal_of_Hubbard_Parameters

