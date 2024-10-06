!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 570 $)
!
!  SUBROUINE: Renewal_of_pWaveFunctions
!
!  AUTHOR(S): M. Saito, T. Yamasaki   November/04/2003
!  
!  FURTHER MODIFICATION: M. Saito, March/26/2007
!
!  The license of the code and contact address:
!  See the files, COPYRIGHT and LICENSE (or LICENSE_J.pdf)                  
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
!
subroutine Renewal_of_pWaveFunctions()
! $Id: Renewal_of_pWaveFunctions.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Const_Parameters, only     :  ON, OFF, DP, ORTHONORMALIZATION &
       &                            , CG, SD, MSD, SUBMAT, MATRIXDIAGON, lmMSD, lmSD, lmCG 
  use m_IterationNumbers,     only :  iteration, iteration_positron_wf
  use m_Control_Parameters,   only :  dtim_p, ipripositron &
       &                            , m_CtrlP_dtim_p_now &
       &                            , m_CtrlP_dtim_1Dsearch_now &
       &                            , m_CtrlP_submat_for_pWFs_now &
       &                            , m_CtrlP_decide_dtim_1Dsearch &
       &                            , m_CtrlP_solver_for_pWFs_now
!!$       &                           , m_CtrlP_On_or_Off_precon_pWFs
  use m_Files, only :                 nfout
  use m_Positron_Wave_Functions,only: m_pWF_renew_WF_by_SDorCG &
       &                            , m_pWF_tell_band_energy &
       &                            , m_pWF_copy_pzaj_to_pzaj_old &
       &                            , m_pWF_wd_pzaj &
       &                            , m_pWF_submat &
       &                            , m_pWF_evolve_wfs_again

! === POSITRON SCF ===
  use m_Const_Parameters,     only : Positron_CONV, Positron_GGGC, BAND_Energy, TOTAL_Energy
  use m_Control_Parameters,   only : positron_method
  use m_Total_Energy,  only : m_TE_tell_total_energy
  use m_Positron_Wave_Functions, only : pchg_l, pchgo_l
! =====

  implicit none 
  integer :: isolver, precon, sw_submat
  real(kind=DP) :: dtim

  isolver = m_CtrlP_solver_for_pWFs_now(iteration_positron_wf)
!!$  isolver = MSD
  sw_submat = m_CtrlP_submat_for_pWFs_now(iteration_positron_wf)

  dtim = m_CtrlP_dtim_p_now(iteration_positron_wf)
  if(ipripositron >= 3) write(nfout,'(" !! dtim = ",f8.4 &
       & ," <<Renewal_of_pWaveFunctions>>")') dtim

!!$  precon = m_CtrlP_On_or_Off_precon_pWFs(iteration_electronic,iteration_ionic)
  precon = OFF
  solver: select case(isolver)
  case (SD, MSD)
     call m_pWF_renew_WF_by_SDorCG(nfout,isolver,precon,dtim)
  case (lmSD, lmMSD)
     call Renew_pWF_by_lmSDorlmCG(isolver,precon,dtim)
!!$     case (RMM, RMM2, RMM2P)
!!$        call m_pWF_renew_WF(nfout,isolver,precon,dtim)
  case default
     write(nfout,'(" !! isolver = ",i5," Renewal_of_pWaveFunctions")') isolver
     call error
     stop ' error at (Renewal_of_WaveFunctions)'
  end select solver

   if(sw_submat == ON .and. (isolver /= SUBMAT .and. isolver /= MATRIXDIAGON)) then
#ifdef TRANSPOSE
      call m_pWF_submat(nfout)
#else
     if(ipripositron >= 1) write(nfout, &
          & '(" -DTRANSPOSE is required in the makefile if you want to set SUBMAT option in the p_solvers")')
#endif
   end if

   pchgo_l = pchg_l

contains
  subroutine Renew_pWF_by_lmSDorlmCG(isolver,precon,dtim)
    integer, intent(in)         :: isolver,precon
    real(kind=DP)               :: dtim
    real(kind=DP), dimension(3) :: etotal
    real(kind=DP), parameter    :: factor = 2
    real(kind=DP)               :: dtim_new = 0.d0, dtim_msdv
    integer                     :: isolver_core, mode, what_is_the_pcore_solver

    integer ::   energy_evaluation_lmm

    if ( positron_method == Positron_CONV ) then
       energy_evaluation_lmm = BAND_ENERGY
    else if ( positron_method == Positron_GGGC ) then
       energy_evaluation_lmm = TOTAL_ENERGY
       energy_evaluation_lmm = BAND_ENERGY
    endif

    isolver_core = what_is_the_pcore_solver(isolver) ! -(in this file)

    mode = ORTHONORMALIZATION

    call m_pWF_copy_pzaj_to_pzaj_old()

    if ( energy_evaluation_lmm == TOTAL_ENERGY ) then
       etotal(1) = m_TE_tell_total_energy()
    else if ( energy_evaluation_lmm == BAND_ENERGY ) then
       etotal(1) = m_pWF_tell_band_energy()
    endif

    if(ipripositron >= 2) write(nfout,'(" !! etotal(1) = ",f10.6)') etotal(1)
    dtim_msdv = m_CtrlP_dtim_1Dsearch_now(dtim)
    if(ipripositron >= 2) write(nfout,'(" !! dtim_msdv = ",f8.4)') dtim_msdv
!!$    if(isolver_core == CG) then
!!$       call m_pWF_decide_CG_direction(precon)  !  -> betacg
!!$       !    ~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$       call m_pWF_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv)
!!$       !    ~~~~~~~~~~~~~~~~~~~~~~~~~
!!$    else
    call m_pWF_renew_WF_by_SDorCG(nfout,isolver_core,precon,dtim_msdv)
       !    ~~~~~~~~~~~~~~~~~~~~~~~~~
!!$    end if

    if ( energy_evaluation_lmm == TOTAL_ENERGY ) then
       etotal(2) = m_TE_tell_total_energy()
    else if(energy_evaluation_lmm == BAND_ENERGY) then
       etotal(2) = m_pWF_tell_band_energy()
    endif

    if(ipripositron >= 2) write(nfout,'(" !! etotal(2) = ",f10.6)') etotal(2)

!!$    if(isolver_core == CG .or. isolver_core == eazyCG) mode = NORMALIZATION
    call m_pWF_evolve_WFs_again(nfout,mode,dtim_msdv,factor*dtim_msdv)

    if ( energy_evaluation_lmm == TOTAL_ENERGY ) then
       etotal(3) = m_TE_tell_total_energy()
    else if(energy_evaluation_lmm == BAND_ENERGY) then
       etotal(3) = m_pWF_tell_band_energy()
    endif

    if(ipripositron >= 2) write(nfout,'(" !! etotal(3) = ",f10.6)') etotal(3)
    dtim_new = m_CtrlP_decide_dtim_1Dsearch(nfout,etotal,dtim_msdv,factor)
    mode = ORTHONORMALIZATION
    call m_pWF_evolve_WFs_again(nfout,mode,factor*dtim_msdv,dtim_new)

  end subroutine Renew_pWF_by_lmSDorlmCG

end subroutine Renewal_of_pWaveFunctions

subroutine Solve_pWaveFunctions()
! $Id: Renewal_of_pWaveFunctions.F90 570 2017-04-21 20:34:50Z yamasaki $
  use m_Control_parameters, only :      pev_max_iteration, ipripositron
  use m_IterationNumbers, only :        iteration_positron_wf, m_Iter_positron_incre
  use m_Files, only :                   nfout
  use m_Positron_Wave_Functions, only : m_pWF_construct_pcharge &
       &                              , m_pWF_charge_rspace, m_pWF_wd_pev &
       &                              , m_pWF_wlifetime
!!  use m_Electronic_Structure, only :    m_ES_wd_vlhxcQ
  implicit none
  logical :: pEigenValues_are_Converged

!!  call Initial_pWaveFunctions()
!!  if(ipripositron >= 2) call m_ES_wd_vlhxcQ()

  Solve_pWaveFunction: do
     call m_Iter_positron_incre()
     call Renewal_of_pWaveFunctions()
     if(pEigenValues_are_Converged()) exit Solve_pWaveFunction
     call m_pWF_construct_pcharge()
     call m_pWF_charge_rspace()
  end do Solve_pWaveFunction

  call m_pWF_construct_pcharge()
  call m_pWF_charge_rspace()
  call m_pWF_wlifetime()

contains
  logical function pIterationNumber_limit()
    if(pev_max_iteration <= iteration_positron_wf) then
       pIterationNumber_limit = .true.
    else
       pIterationNumber_limit = .false.
    end if
  end function pIterationNumber_limit
end subroutine Solve_pWaveFunctions

logical function pEigenValues_are_Converged()
  use m_Files, only :                   nfout
  use m_IterationNumbers, only :        iteration_positron_wf
  use m_Control_Parameters, only :      ipripositron, pev_max_iteration
  use m_Positron_Wave_Functions, only : m_pWF_pevdff, m_pWF_wd_pev
  implicit none
  pEigenValues_are_Converged = m_pWF_pevdff()

  if(pEigenValues_are_Converged) then
     if(ipripositron >= 1) &
          & write(nfout,'(" Positron Wave Functions have converged at ",i7," -th iteration")') &
          & iteration_positron_wf
     call m_pWF_wd_pev(nfout)
  else if(iteration_positron_wf >= pev_max_iteration) then
     if(ipripositron >= 1) then
        write(nfout,'(" Positron Wave Functions have not converged")')
        write(nfout,'(" iteration_positrion_wf have reached to pev_max_iteration (= ",i5," )")') &
             & pev_max_iteration
     end if
     call m_pWF_wd_pev(nfout)
     pEigenValues_are_Converged = .true.
  end if

end function pEigenValues_are_Converged

subroutine Initial_pWaveFunctions()
  use m_Const_Parameters, only   : by_random_numbers, by_matrix_diagon
  use m_Control_Parameters, only : intpzaj, ipripositron
  use m_Files, only :              nfout
  use m_positron_Wave_Functions, only : m_pWF_IW_by_randomnumbers &
       &                              , m_pWF_modified_gram_schmidt &
       &                              , m_pWF_energy_eigen_values &
       &                              , m_pWF_allocate_pzaj_etc &

       &                              , m_pWF_wd_pev &
       &                              , m_pWF_alloc_afft_etc &
       &                              , m_pWF_dealloc_afft_etc &
       &                              , m_pWF_wd_pzaj

! == POSITRON SCF === 2015/11/28
  use m_epc_Potential,  only : m_epc_alloc, m_epc_alloc_vlhxc_p
! =================== 2015/11/28

  implicit none

  call m_epc_alloc
  call m_epc_alloc_vlhxc_p

  call m_pWF_allocate_pzaj_etc()

  if(intpzaj == by_random_numbers) then
     call m_pWF_IW_by_randomnumbers()
     call m_pWF_wd_pzaj(nfout,'after by_randum_numbers',23)
     call m_pWF_modified_gram_schmidt()
     call m_pWF_wd_pzaj(nfout,'after modified_gram_schmidt',27)
     call m_pWF_alloc_afft_etc()
     call m_pWF_energy_eigen_values()
     call m_pWF_wd_pzaj(nfout,'after energy_eigen_values',25)
     call m_pWF_dealloc_afft_etc()
  else if(intpzaj == by_matrix_diagon) then
  end if
  if(ipripositron >= 1) write(nfout,'(" --- initial positron energy eigen values ---")')
  call m_pWF_wd_pev(nfout)

end subroutine Initial_pWaveFunctions

integer function what_is_the_pcore_solver(isolver)
  use m_Const_Parameters,  only : lmSD, SD, lmMSD, MSD, lmCG, CG
  integer, intent(in) :: isolver

  if(isolver == lmSD) then
     what_is_the_pcore_solver = SD
  else if(isolver == lmMSD) then
     what_is_the_pcore_solver = MSD
!!$  else if(isolver == lmCG) then
!!$     what_is_the_core_solver = CG
  else
     what_is_the_pcore_solver = MSD
  end if
end function what_is_the_pcore_solver
