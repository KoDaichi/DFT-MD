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
     if ( sw_epsilon == ON ) then
        call using_epsilon_ek
     else if ( sw_excitation == ON ) then
        call using_excitation
     endif
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

subroutine using_excitation
  use m_Const_Parameters,  only : ON, OFF, CMPLDP
  use m_Control_Parameters,  only : sw_phonon_with_epsilon, sw_use_add_proj, &
       &                            sw_corelevel_spectrum

  use m_Excitation,  only : m_XI_set_ilmt_val, m_XI_init_spectrum, &
       &                    m_XI_set_ngpt_XI, &
       &                    m_XI_set_qitg_XI_vv, &
       &                    m_XI_set_dl2p_XI_vv, &
       &                    m_XI_set_mat_dipole_corr_vv, &
       &                    m_XI_calc_transition_moment_vv, &
       &                    m_XI_set_RhoTilde_vv, &
       &                    m_XI_calc_spectr_fn_vv, &
       &                    m_XI_calc_spectr_tensor_vv, &
       &                    SpectrFn_vv, SpectrTensor_vv, &
       &                    m_XI_set_qitg_XI_vc, &
       &                    m_XI_set_dl2p_XI_vc, &
       &                    m_XI_set_ilmt_core, &
       &                    m_XI_set_mat_dipole_corr_vc, &
       &                    m_XI_calc_transition_moment_vc,  &
       &                    m_XI_set_RhoTilde_vc, &
       &                    m_XI_calc_spectr_fn_vc, &
       &                    m_XI_calc_transition_moment_cc,  &
       &                    m_XI_set_RhoTilde_cc, &
       &                    m_XI_dealloc_arrays, &
       &                    SpectrFn_vc, &
       &                    nstep, e_step, sw_tddft, sw_bse, &
       &                    Kernel_Coulomb, m_XI_set_Kernel_Coulomb, epsinv_omega0, &
       &                    sw_spin_decomposition, m_XI_calc_spectr_fn_vv_spindec, &
       &                    sw_write_wavederf_file, m_XI_wd_WAVEDERF

  use m_ValenceBand_Spectrum,  only : m_VBS_set_data_ppc_from_pp
  use m_CoreLevel_Spectrum,   only : m_CLS_set_data_core2val_from_pp, &
       &                             m_CLS_alloc_wfn_core_states, &
       &                             m_CLS_set_wfn_core_states, &
       &                             m_CLS_set_corelevels_virtually, &
       &                             m_CLS_set_ene_core_states

  use m_ES_nonlocal,    only : m_ES_add_betar_dot_WFs
  use m_Files,  only : nfout
  use m_Raman,  only : dielectric, m_Raman_write_dielec_tensors
  use m_Kpoints,  only : m_Kp_set_star_of_k

  implicit none

  if ( sw_phonon_with_epsilon == ON ) then        ! nmax_G == 1
     call case_phonon_with_epsilon
  else
     call case_epsilon_fn
  endif

contains

  subroutine case_phonon_with_epsilon
    logical :: First = .true.

    if ( First ) then
       call m_Kp_set_star_of_k
       call m_XI_set_ilmt_val
       call m_VBS_set_data_ppc_from_pp
       call m_XI_set_mat_dipole_corr_vv
       call m_XI_init_spectrum
       First = .false.
    endif

    if ( sw_use_add_proj == ON ) call m_ES_add_betar_dot_WFs(nfout)

    call m_XI_calc_transition_moment_vv
    call m_XI_calc_spectr_tensor_vv

#ifdef USE_ASMS_EXCITATION
    call ASMS_XI_set_dielectric_tensor( SpectrTensor_vv, dielectric )
#endif

    call m_Raman_write_dielec_tensors

  end subroutine case_phonon_with_epsilon

  subroutine case_epsilon_fn
    integer :: i

    logical :: First = .true.

    if ( First ) then
       if ( sw_write_wavederf_file == OFF ) call m_Kp_set_star_of_k
       call m_XI_set_ilmt_val
       call m_XI_set_ngpt_XI
!
       if ( sw_corelevel_spectrum == ON ) then
          call m_CLS_set_data_core2val_from_pp
          call m_CLS_alloc_wfn_core_states( .true. )
          call m_CLS_set_wfn_core_states( .true. )
!          call m_CLS_set_corelevels_virtually
          call m_CLS_set_ene_core_states
          call m_XI_set_ilmt_core
          call m_XI_set_mat_dipole_corr_vc
!
          call m_XI_set_qitg_XI_vc
          call m_XI_set_dl2p_XI_vc

          if ( sw_bse == ON ) then
             call m_VBS_set_data_ppc_from_pp
             call m_XI_set_mat_dipole_corr_vv
          endif
       else
          call m_XI_set_qitg_XI_vv
          call m_XI_set_dl2p_XI_vv
          call m_VBS_set_data_ppc_from_pp
          call m_XI_set_mat_dipole_corr_vv
       endif

       if ( sw_write_wavederf_file == OFF ) call m_XI_init_spectrum
       First = .false.
    endif

    if ( sw_use_add_proj == ON ) call m_ES_add_betar_dot_WFs(nfout)

    if ( sw_corelevel_spectrum == ON ) then
       if ( sw_bse == ON ) then
          call m_XI_calc_transition_moment_vv
          call m_XI_calc_transition_moment_vc
          call m_XI_calc_transition_moment_cc
          call m_XI_set_RhoTilde_vv
          call m_XI_set_RhoTilde_vc
          call m_XI_set_RhoTilde_cc

          call m_XI_calc_spectr_fn_vv
          call m_XI_set_Kernel_Coulomb
#ifdef USE_ASMS_EXCITATION
          call ASMS_XI_epsinv_omg0_fr_spctrfn( SpectrFn_vv, Kernel_Coulomb, &
               &                               epsinv_omega0 )
#endif

       else if ( sw_tddft == ON ) then
          call m_XI_calc_transition_moment_vc
          call m_XI_set_RhoTilde_vc
          call m_XI_calc_spectr_fn_vc
          call m_XI_set_Kernel_Coulomb
#ifdef USE_ASMS_EXCITATION
          call ASMS_XI_calc_tddft_spectrum( SpectrFn_vc, Kernel_Coulomb )
#endif
       else
          call m_XI_calc_transition_moment_vc
          call m_XI_calc_spectr_fn_vc
#ifdef USE_ASMS_EXCITATION
          call ASMS_XI_calc_indep_spectrum( SpectrFn_vc )
#endif
       endif

    else
       if ( sw_bse == ON ) then
          call m_XI_calc_transition_moment_vv
          call m_XI_set_RhoTilde_vv
          call m_XI_calc_spectr_fn_vv
          call m_XI_set_Kernel_Coulomb
#ifdef USE_ASMS_EXCITATION
          call ASMS_XI_epsinv_omg0_fr_spctrfn( SpectrFn_vv, Kernel_Coulomb, &
               &                               epsinv_omega0 )
#endif
          stop "PPP"
       else if ( sw_tddft == ON ) then
          call m_XI_calc_transition_moment_vv
          call m_XI_set_RhoTilde_vv

!!          call m_XI_calc_spectr_fn_vv
          if ( sw_spin_decomposition == ON ) then
             call m_XI_calc_spectr_fn_vv_spindec
          else
             call m_XI_calc_spectr_fn_vv
          endif

          call m_XI_set_Kernel_Coulomb
#ifdef USE_ASMS_EXCITATION
          call ASMS_XI_calc_tddft_spectrum( SpectrFn_vv, Kernel_Coulomb )
#endif
       else
          if ( sw_write_wavederf_file == ON ) then
             call m_XI_wd_WAVEDERF
          else
             call m_XI_calc_transition_moment_vv
!!            call m_XI_calc_spectr_fn_vv
             if ( sw_spin_decomposition == ON ) then
                call m_XI_calc_spectr_fn_vv_spindec
             else
                call m_XI_calc_spectr_fn_vv
             endif
#ifdef USE_ASMS_EXCITATION
             call ASMS_XI_calc_indep_spectrum( SpectrFn_vv )
#endif
          endif
       endif
    endif

  end subroutine case_epsilon_fn

end subroutine Using_excitation
