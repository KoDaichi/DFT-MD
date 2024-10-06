!! Last updated: 01 March 2012
!! real-time electron dynamics
!! based on the time-dependent density functional theory
!!
!! mdmain.F90      -->  rttddft_main.f90
!! rttddft_main.f90 -->  m_Control_Parameters.F90 (input parameters)
!!                  -->  m_rttddft.f90
!!
!! some unstable parts commented out, 12 Sep 2012 at isir, osaka univ.
!!

subroutine rttddft_main

use m_Const_Parameters,only : &
    DP,ON,OFF,FORCE_CONVERGED,INITIAL,CONTINUATION,FIXED_CHARGE,FIXED_CHARGE_CONTINUATION &
   ,Valence_plus_PC_Charge,VXC_AND_EXC, AU_TIME
use m_Control_Parameters,only : &
    iconvergence,iconvergence_previous_job &
   ,sw_rttddft,time_step_max,time_step_delta,propagator_method,propagator_order
use m_IterationNumbers,only : &
    iteration, iteration_electronic
use m_Files,only : &
    nfout
use m_Kpoints,only : &
    kv3
use m_Electronic_Structure,only : &
    m_ES_energy_eigen_values
use m_ES_LHXC,only : &
    m_ESlhxc_potential
use m_Charge_Density,only : &
    chgq_l,m_CD_softpart,m_CD_hardpart,m_CD_cp_chgq_to_chgqo
use m_XC_Potential,only : &
    vxc_l, m_XC_cal_potential
use m_Total_Energy,only : &
    m_TE_total_energy
use m_ES_Intgr_VlhxcQlm,only : &
    m_ESiVQ_integrate_VlhxcQlm
use m_ES_WF_by_MatDiagon,only : &
    m_ESmat_solve_Hx_eq_eSx
use m_rttddft,only : &
    m_rttddft_init_occup_control,m_rttddft_init_impulse_field &
   ,m_rttddft_print_eko_occup,m_rttddft_wd_charge &
   ,m_rttddft_propagate_wf_by_taylor,m_rttddft_propagate_wf_by_split &
   ,m_rttddft_check_wf,m_rttddft_dipole,m_rttddft_current

implicit none
integer :: i_time,it
logical :: TotalEnergy_is_Divergent

if(sw_rttddft == OFF) return
if(iconvergence < FORCE_CONVERGED .and. iconvergence_previous_job < FORCE_CONVERGED) return


write(nfout,'(/,"#### RT-TDDFT real-time electron dynamics: start ####")')

!!
!!-- Initial Electronic Structure
!!
write(nfout,'("# time_step=     0")')

call m_rttddft_init_occup_control  !<<
call m_rttddft_init_impulse_field  !<<

call m_CD_softpart(nfout,kv3)
call m_CD_hardpart(nfout,kv3)
call m_CD_cp_chgq_to_chgqo

call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC)
call m_ESlhxc_potential(nfout,chgq_l,vxc_l)
call m_ESiVQ_integrate_VlhxcQlm(nfout)

call m_ES_energy_eigen_values(nfout)
call m_TE_total_energy(nfout,.true.,kv3)
if(TotalEnergy_is_Divergent()) write(nfout,'("** WARNING: initial")')

! call m_rttddft_wd_charge      !<<
! call m_rttddft_check_wf
call m_rttddft_dipole
call m_rttddft_current

time_step: do i_time=1,time_step_max

  write(nfout,'(/,"# time_step=",i6,2x,"time=",e12.4," au =",e12.4," fs")') &
    i_time,time_step_delta*i_time,time_step_delta*i_time* AU_TIME * 1.0d15  ! 0.0241888

  call IterationNumber_Setting

!!
!!-- Time-Propagation of Wavefunction
!!
  if(propagator_method==1) then
    call m_rttddft_propagate_wf_by_taylor   !<<
  elseif(propagator_method==2) then
    stop '**** ERROR(RT-TDDFT): propagator_method=2 not recommended, pls set 1 '
    call m_rttddft_propagate_wf_by_split(i_time)    !<<
  else
    stop '**** ERROR(RT-TDDFT): invalid input, propagator_method '
  endif

!!
!!-- Charge Density Update
!!
  call m_CD_cp_chgq_to_chgqo
  call m_CD_softpart(nfout,kv3)
  call m_CD_hardpart(nfout,kv3)

!!
!!-- Potential Update
!!
  call m_XC_cal_potential(nfout,Valence_plus_PC_Charge,chgq_l,VXC_AND_EXC)
  call m_ESlhxc_potential(nfout,chgq_l,vxc_l)
  call m_ESiVQ_integrate_VlhxcQlm(nfout)

!!
!!-- Energy Eigen Values (Diagonal Part)
!!
  call m_ES_energy_eigen_values(nfout)

!!
!!-- Total Energy and Ion Kinetic Energy
!!
  call m_TE_total_energy(nfout,.true.,kv3)
  if(TotalEnergy_is_Divergent()) exit time_step

!!
!!-- Print for Accuracy Check
!!
  if(i_time==1.or.mod(i_time,1000)==0.or.i_time==time_step_max) then
    call m_rttddft_print_eko_occup
  ! call m_rttddft_check_wf
  endif

!!
!!-- Calculations of Physical Quantities
!!
  call m_rttddft_dipole
  call m_rttddft_current

!!
!!-- Check Off-Diagonal Matrix Element
!!
!!  call m_ESmat_solve_Hx_eq_eSx(nfout,iteration,iteration_electronic)

!!
!!-- Calculate Force and Move Ions (velocity verlet should be selected)
!!
! call Forces
! call Move_Ions
! call Ewald_and_Structure_Factor
  call MDIterationNumber_Setting

enddo time_step

write(nfout,'("#### RT-TDDFT real-time electron dynamics: end ####",/)')

end subroutine rttddft_main
