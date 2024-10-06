!
! =========== Contributions ===================================
!
! Through the courtesy of contributors, the following functions are added.
!
! Company:  ASMS Co.,Ltd.
! Functions:  [Identifier: 13.1R]
!            Dielectric functions can be calculated when condition = INITIAL 
!
! Functions:  [Identifier: 13.1XI]
!            Dielectric functions can be calculated in the EXCITATION module
!
! =============================================================
subroutine Epsilon_Postscf

  use m_Const_Parameters, only : ON
  use m_Control_Parameters,  only : sw_phonon_with_epsilon, sw_excitation
  use m_Epsilon_ek,  only : sw_epsilon, auto_mode

  implicit none

     if(auto_mode/=0) sw_epsilon = on
!  if ( sw_phonon_with_epsilon == ON ) then
  !   if ( sw_epsilon == ON ) then
        call using_epsilon_ek
        write(6,'(" using_epsilon_ek")')
  !   else if ( sw_excitation == ON ) then
  !      call using_excitation
  !      write(6,'(" using_excitation")')
  !   endif
!  else
!     if ( sw_excitation == ON ) call using_excitation
!     if ( sw_epsilon == ON )    call using_epsilon_ek
!  endif


end subroutine Epsilon_Postscf

subroutine using_epsilon_ek

  use m_Control_Parameters, only : sw_phonon_with_epsilon
  use m_Const_Parameters,  only : on

  use m_KPoints, only : kv3, kv3_ek, m_Kp_cp_vkxyz_to_vkxyz_ek, &
       &                m_Kp_alloc_kpoints_ek
  use m_IterationNumbers, only : nk_in_the_process, nk_converged
!
  use m_Electronic_Structure, only : m_ES_cp_eko_l_to_eko_ek2, &
       &                             m_ES_alloc_eko_ek

  use m_Epsilon_ek,  only : set_dielectric_tensor
  use m_Raman,  only : m_Raman_write_dielec_tensors

  implicit none

  logical, save :: FirstFlag = .true.

  nk_in_the_process = 1
  nk_converged = kv3

  call m_Kp_alloc_kpoints_ek
  call m_Kp_cp_vkxyz_to_vkxyz_ek
  call m_ES_alloc_eko_ek
  call m_ES_cp_eko_l_to_eko_ek2
!
  if ( FirstFlag ) then
     call Initialization_Epsilon(0)
     FirstFlag = .false.
  else
     call Initialization_Epsilon(1)
  endif

  call Transition_moment_Epsilon
  call Prep_for_Calc_Epsilon
  call Calc_Epsilon
  call Calc_Nonlinear_optics

  if ( sw_phonon_with_epsilon  == ON ) then
     call set_dielectric_tensor
     call m_Raman_write_dielec_tensors
  else
     call WriteDownData_onto_Files_Eps
  endif

end subroutine Using_epsilon_ek

