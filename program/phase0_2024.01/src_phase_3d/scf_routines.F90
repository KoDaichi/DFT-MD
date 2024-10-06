!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 599 $)
!
!  SUBROUINE: scf_initialize, scf_do_scf_and_force, scf_finalize
!
!  AUTHOR(S): J. Koga March/24/2009 
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
!***************************************************************
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
!!!!!BRANCH_P ORG_Parallel
! ==============================================================================
subroutine scf_initialize()
  implicit none
  call Initialization()
  call InputData_Analysis()
  call Preparation(0)                 ! Basis set, symmetry check etc.
  call Preparation_for_mpi(1)         ! mpi
  call PseudoPotential_Construction()
  call Ewald_and_Structure_Factor()
  call Initial_Electronic_Structure()
  call Initial_MD_Condition()
end subroutine scf_initialize

subroutine scf_do_scf_and_force(exi,skip_ewald)
  implicit none
  logical, intent(out) :: exi
  logical, intent(in)  :: skip_ewald
  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Already_Converged, Already_Converged2
  logical  :: Positron_bulk, Positron_defect
  logical  :: Hubbard_model
  logical  :: Ending_Time
  logical  :: tor
  exi = .false.
  if(.not.skip_ewald) call Ewald_and_Structure_Factor()
  ChargeDensity:    do
    call IterationNumber_Setting()
    call Renewal_of_WaveFunctions()
    call ChargeDensity_Construction(1)
    call ChargeDensity_Mixing()
    if(Ending_Time()) then
      exi = .true.
      return
    endif
    if(TotalEnergy_is_Divergent()) then
      exi = .true.
      return
    endif
    call Renewal_of_Potential()
    if(Hubbard_model()) then
      call Renewal_of_Hubbard_Potential()
    end if
    if(ChargeDensity_is_Converged())  then
#ifdef LIBRARY_BUILD
      call Forces_phase0()
#else
      call Forces()
#endif
      call post_force()
      call MDIterationNumber_Setting()
      return
    endif
  enddo ChargeDensity
end subroutine scf_do_scf_and_force

subroutine post_force()
end subroutine post_force

subroutine scf_finalize()
  implicit none
  call Postprocessing(.false.)
  call WriteDownData_onto_Files(.true.)
  !call Finalization_of_mpi()          ! mpi
end subroutine scf_finalize

subroutine scf_rd_wf_and_chg(logi)
  use m_Files, only : nfout,nfzaj,nfchgt,F_ZAJ_in_partitioned &
 &                  , F_CHGT_in_partitioned, m_Files_reopen_nfchgt, m_Files_reopen_nfzaj &
 &                  , nfefermi, m_Files_open_nfefermi, m_Files_close_nfefermi,F_ZAJ,F_CHGT
  use m_ES_IO, only : m_ESIO_rd_WFs, m_ESIO_rd_Efermi
! === DEBUG by tkato 2011/11/09 ================================================
! use m_Electronic_Structure, only : m_ES_betar_dot_WFs, m_ES_energy_eigen_values
  use m_Electronic_Structure, only : m_ES_energy_eigen_values_3D
!  use m_ES_nonlocal,          only : m_ES_betar_dot_WFs_3D
! ==============================================================================
  use m_Charge_Density, only : m_CD_rd_chgq, m_CD_wd_chgq_l_small_portion, chgq_l
  use m_XC_Potential, only: m_XC_cal_potential_3D, vxc_l
  use m_Orbital_Population, only : m_OP_rd_occ_mat, m_OP_mix_om
  use m_Const_Parameters, only : ON, Valence_plus_PC_Charge, VXC_AND_EXC
  use m_Control_Parameters, only : printable, sw_hubbard
  use m_ES_LHXC, only : m_ESlhxc_potential_3D
  use m_ES_Intgr_VlhxcQlm,  only : m_ESiVQ_integrate_VlhxcQlm_3D
  use m_Kpoints, only : kv3
  use m_Control_Parameters, only : af
  use m_Parallelization, only : map_k, myrank_k
  use m_ES_nonlocal,     only : m_ES_betar_dot_WFs_4_each_k_3D
  implicit none

  logical, intent(out) :: logi
  logical :: exi
  integer :: ik
  logi = .true.
  inquire(file=trim(F_CHGT),exist=exi)
  if(.not.exi) then
    logi=.false.
    return
  endif
  inquire(file=trim(F_ZAJ),exist=exi)
  if(.not.exi)then
    logi=.false.
    return
  endif
  
  call m_Files_reopen_nfzaj()
  call m_ESIO_rd_WFs(nfout,nfzaj,F_ZAJ_in_partitioned)
  call m_Files_reopen_nfchgt()
  call m_CD_rd_chgq(nfout,nfchgt,F_CHGT_in_partitioned)

  call m_Files_open_nfefermi()
  call m_ESIO_rd_Efermi(nfout,nfefermi)
  call m_Files_close_nfefermi()

  do ik = 1, kv3, af+1
     if(map_k(ik) /= myrank_k) cycle         ! MPI
     call m_ES_betar_dot_WFs_4_each_k_3D(nfout,ik)
  end do

  call m_CD_wd_chgq_l_small_portion(nfout)
  call Renewal_of_Potential()
  !call m_XC_cal_potential_3D(nfout,Valence_plus_PC_Charge,chgq_l, VXC_AND_EXC) ! -> vxc_l
  !call m_ESlhxc_potential_3D(nfout,chgq_l,vxc_l) ! (stlhxc) ->vlhxc_l
  !call m_ESiVQ_integrate_VlhxcQlm_3D(nfout) ! (lclchh) -> vlhxcQ
  !call m_ES_energy_eigen_values_3D(nfout)   ! (eigen0) -> eko_l,neordr

  if(sw_hubbard == ON) then
     !!$call Renewal_of_OccMat(.true.) ! -> om
     call m_OP_rd_occ_mat(nfout) ! -> om
     call m_OP_mix_om(1.d0) ! om -> ommix
     call Renewal_of_Hubbard_Potential() ! ommix -> dhub
  end if

end subroutine scf_rd_wf_and_chg
! === "phase" repository merge: To make 3D_Parallel by tkato 2011/11/10 ========
!!!!BRANCH_P_END ORG_Parallel
! ==============================================================================

