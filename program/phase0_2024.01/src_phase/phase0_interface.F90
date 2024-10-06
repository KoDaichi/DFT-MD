
subroutine phase0_initialize(comm_world, tag, ne, nk, ng, stdout)
  use m_Const_Parameters, only : ON
  implicit none
  integer, intent(in) :: comm_world 
  integer, intent(in) :: tag
  integer, intent(in) :: ne, nk, ng
  integer, intent(in) :: stdout
  logical :: initialization_required
  call Initialization_lib(comm_world,tag,ne,nk,ng,stdout)
  call InputData_Analysis()
  call phase0_initialize_ppetc()
end subroutine phase0_initialize

subroutine phase0_initialize_ppetc()
  call Preparation(0)
  if(initialization_required())then
    call Preparation_for_mpi(1)
  endif
  call PseudoPotential_Construction
#ifdef ENABLE_ESM_PACK
  if(initialization_required())then
     call Preparation_for_ESM
  endif
#endif
  call Ewald_and_Structure_Factor
  call Initial_Electronic_Structure

!  call Initial_MD_Condition()
end subroutine phase0_initialize_ppetc

subroutine phase0_set_coordinate_xyz(coord_x, coord_y, coord_z)
  use m_Const_Parameters,   only : DP
  use m_Ionic_System,       only : natm, cps, m_IS_cps_to_pos
  real(kind=DP), intent(in),  dimension(natm) :: coord_x, coord_y, coord_z
  cps(1:natm,1) = coord_x(1:natm)
  cps(1:natm,2) = coord_y(1:natm)
  cps(1:natm,3) = coord_z(1:natm)
  call m_IS_cps_to_pos()
end subroutine phase0_set_coordinate_xyz

subroutine phase0_set_coordinate(cpsin)
  use m_Const_Parameters,   only : DP
  use m_Ionic_System,       only : m_IS_cps_to_pos, cps, pos
  use m_Ionic_System,       only : natm
  implicit none
  real(kind=DP), intent(in),  dimension(natm,3) :: cpsin
  cps  = cpsin
  call m_IS_cps_to_pos()
end subroutine phase0_set_coordinate

subroutine phase0_set_unitcell(altvin)
  use m_Const_Parameters,   only : DP, ON, DRIVER_SC_DFT
  use m_Crystal_Structure,  only : m_CS_altv_2_rltv,p2bmat,a,b,c,ca,cb,cc,il,kt_iuctype,rltv &
  &                              , altv, univol,rvol, b2pmat
  use m_Control_Parameters, only : sw_optimize_lattice, driver
  use m_IterationNumbers,   only : iteration_unit_cell, iteration_scdft
  use m_Ionic_System,       only : natm, cps, pos
  use m_Files,              only : nfout,m_Files_set_default_filenames, m_Files_rd_file_names_data
  implicit none
  real(kind=DP), intent(in),  dimension(3,3)    :: altvin
  altv = altvin
  driver = DRIVER_SC_DFT
  iteration_scdft = 2
  call altv_2_rltv(altv,rltv,univol,rvol)  ! in b_CS
  call change_of_coordinate_system(altv,pos,natm,natm,cps) !-(b_I.S.) pos -> cps
  call primitive2bravais(nfout,p2bmat,altv(:,1),altv(:,2),altv(:,3),a,b,c,ca,cb,cc,il) ! in b_CS
  call Array_Deallocate()
  call Continuation_Mode()
  call phase0_initialize_ppetc()
end subroutine phase0_set_unitcell

logical function phase0_stress_enabled()
  use m_Const_Parameters,   only : ON
  use m_Control_Parameters, only : istress
  implicit none
  phase0_stress_enabled = istress == ON
end function phase0_stress_enabled

subroutine phase0_scf(energy, force_x, force_y, force_z, stress_tensor)
  use m_Control_Parameters, only : istress
  use m_Const_Parameters,   only : DP, ON
  use m_Ionic_System,       only : natm,cps
  use m_Force,              only : forc_l
  use m_Stress,             only : m_Stress_get_curr_stress
  use m_Total_Energy,       only : etotal
  implicit none
  real(kind=DP), intent(out)                  :: energy
  real(kind=DP), intent(out), dimension(natm) :: force_x, force_y, force_z
  real(kind=DP), intent(out), dimension(3,3)  :: stress_tensor

  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Hubbard_model
  logical  :: Ending_Time
  stress_tensor=0.d0
  call Ewald_and_Structure_Factor
  ChargeDensity:    do
    call IterationNumber_Setting()
    call Renewal_of_WaveFunctions()
    call ChargeDensity_Construction(1)
    call Potential_Construction()
    call ChargeDensity_Mixing
    if(Ending_Time()) then
      return
    endif
    if(TotalEnergy_is_Divergent()) then
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
      call MDIterationNumber_Setting()
      energy = etotal
      force_x(1:natm)  = forc_l(1:natm,1)
      force_y(1:natm)  = forc_l(1:natm,2)
      force_z(1:natm)  = forc_l(1:natm,3)
      if(istress==ON) then
        call Stress()
        stress_tensor = m_Stress_get_curr_stress()
      endif
      exit ChargeDensity
    endif
  enddo ChargeDensity
  call MDIterationNumber_Setting2()
  call WriteDownData_onto_Files(.false.)
end subroutine phase0_scf

subroutine phase0_energy_force_stress(energy, force, stress_tensor)
  use m_Control_Parameters, only : istress
  use m_Const_Parameters,   only : DP, ON
  use m_Ionic_System,       only : natm,cps
  use m_Force,              only : forc_l
  use m_Stress,             only : m_Stress_get_curr_stress
  use m_Total_Energy,       only : etotal
  implicit none
  real(kind=DP), intent(out)                    :: energy
  real(kind=DP), intent(out), dimension(natm,3) :: force
  real(kind=DP), intent(out), dimension(3,3)    :: stress_tensor

  logical  :: ChargeDensity_is_Converged, TotalEnergy_is_Divergent
  logical  :: Hubbard_model
  logical  :: Ending_Time
  energy = 0.d0; force = 0.d0; stress_tensor=0.d0
  call Ewald_and_Structure_Factor
  ChargeDensity:    do
    call IterationNumber_Setting()
    call Renewal_of_WaveFunctions()
    call ChargeDensity_Construction(1)
    call Potential_Construction()
    call ChargeDensity_Mixing
    if(Ending_Time()) then
      return
    endif
    if(TotalEnergy_is_Divergent()) then
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
      energy = etotal
      force  = forc_l
      if(istress==ON) then
        call Stress()
        stress_tensor = m_Stress_get_curr_stress()
      endif
      exit ChargeDensity
    endif
  enddo ChargeDensity
  call MDIterationNumber_Setting2()
  call WriteDownData_onto_Files(.false.)
end subroutine phase0_energy_force_stress

subroutine phase0_finalize()
  use m_Files, only : m_Files_close_all, m_Files_close_logfile
  use m_Control_Parameters, only : m_CtrlP_wd_cpu_total
  implicit none
  call m_CtrlP_wd_cpu_total()
  call m_Files_close_all()                    ! -(m_Files)
  call PrintStatus()
  call m_Files_close_logfile()
end subroutine phase0_finalize

  subroutine Array_Deallocate()
     use m_Const_Parameters, only : OFF, ON, DRIVER_URAMP, DRIVER_SC_DFT, DRIVER_NEB
     use m_Control_Parameters, only : m_CtrlP_dealloc,sw_rebuild_pws,m_CtrlP_set_init_status, driver
     use m_Crystal_Structure, only : m_CS_dealloc
     use m_Ionic_System, only : m_IS_dealloc
     use m_PlaneWaveBasisSet, only : m_pwBS_dealloc
     use m_Parallelization, only : m_Parallel_dealloc,m_Parallel_dealloc_mpi_nlmta,m_Parallel_dealloc_mpi_exx &
                        &        , m_Parallel_cp_g1k
     use m_Kpoints, only : m_Kp_dealloc, kv3
     use m_Force, only : m_Force_dealloc
     use m_Charge_Density, only : m_CD_dealloc
     use m_XC_Potential, only : m_XC_dealloc_vxc
     use m_PseudoPotential, only : m_PP_dealloc, flg_paw
     use m_NonLocal_Potential, only : m_NLP_dealloc
     use m_Electronic_Structure, only : m_ES_dealloc
     use m_ES_WF_by_SDorCG, only : m_ESsd_dealloc
     use m_PAW_ChargeDensity, only : m_PAW_dealloc
     use m_PAW_XC_Potential, only : m_PAW_XC_dealloc_vxc
     use m_ES_wf_extrpl, only : m_ES_wf_extrpl_dealloc

! ===== KT_add ============== 13.0AS
     use m_Control_Parameters,  only : num_projectors, sw_hubbard
     use m_Orbital_Population,only: m_OP_dealloc
     use m_Electronic_Structure,  only : m_ES_dealloc_Dhub
! ============================13.0AS

! ======== KT_add ============= 2013/10/31
     use m_Control_Parameters, only: noncol, SpinOrbit_Mode
     use m_Const_Parameters,   only: Neglected, BuiltIn, ByPawPot, ZeffApprox, &
          &                          ByProjector, ReadFromPP
     use m_SpinOrbit_Potential,only: m_SO_dealloc_Dsoc, m_SO_dealloc_Mat_SOC_Strenth
! ============================= 2013/10/31

! ========= KT_add ========= 13.0U2
  use m_Control_Parameters,   only : sw_modified_TFW_functional
  use m_ThomasFermiW_Potential,  only : m_TFW_dealloc_ChgDensityBasisFn
! ========================== 13.0U2

! ====== KT_add === 2014/08/01
  use m_Orbital_QuantumNum, only : m_OP_Qnum_dealloc_array
! ================= 2014/08/01

  use m_Ldos, only : m_Ldos_dealloc

#ifdef FFTW3
  use m_FFT, only : m_FFT_finalize
#endif

     implicit none

     if (sw_rebuild_pws==OFF .or. driver == DRIVER_URAMP .or. driver == DRIVER_SC_DFT ) then
        call m_PP_dealloc()
        call m_NLP_dealloc()
        call m_Parallel_dealloc_mpi_nlmta()
        call m_Parallel_dealloc_mpi_exx()
        call m_Kp_dealloc
!        call m_ES_dealloc
!        call m_ESsd_dealloc
        return
     endif

     call m_CtrlP_dealloc
     call m_ES_dealloc
     call m_Parallel_cp_g1k(kv3)
     call m_Parallel_dealloc
     call m_CS_dealloc
     call m_IS_dealloc
  !call m_Parallel_dealloc_mpi_elec
     call m_Kp_dealloc
     call m_PP_dealloc
     call m_CD_dealloc
     call m_Force_dealloc
     call m_NLP_dealloc
     call m_ESsd_dealloc
     call m_XC_dealloc_vxc
     call m_pwBS_dealloc
     call m_PAW_dealloc
     if(flg_paw) then
         call m_PAW_XC_dealloc_vxc
     endif

! ==== KT_add ============== 2013/10/31
     if ( noncol ) then
        if ( SpinOrbit_Mode == ByPawPot .or. SpinOrbit_Mode == ZeffApprox &
             &                          .or. SpinOrbit_Mode == ReadFromPP ) then
           call m_SO_dealloc_Mat_SOC_strenth
           call m_SO_dealloc_Dsoc
        endif
        if ( SpinOrbit_Mode == ByProjector ) then
           call m_SO_dealloc_Dsoc
        endif
     endif
! ====================== 2013/10/31

! ====== KT_add ======== 13.0U2
     if ( sw_modified_TFW_functional /= OFF ) then
        call m_TFW_dealloc_ChgDensityBasisFn
     end if
! ====================== 13.0U2

! ====== KT_add =========== 13.0AS
     if (num_projectors>0)  call m_OP_dealloc
     if (sw_hubbard == ON) call m_ES_dealloc_Dhub
! ========================= 13.0AS

! ==== KT_add ==== 2014/08/01
     call m_OP_Qnum_dealloc_array
! ================ 2014/08/01

     call m_ES_wf_extrpl_dealloc()
     call m_Ldos_dealloc()
     call m_ES_wf_extrpl_dealloc()
#ifdef FFTW3
     call m_FFT_finalize()
#endif
     call m_CtrlP_set_init_status(.true.)
   end subroutine Array_Deallocate

  subroutine Continuation_Mode()
    use m_Const_Parameters,   only : COORDINATE_CONTINUATION,ON,OFF
    use m_Control_Parameters, only : icond, sw_optimize_lattice,sw_rebuild_pws
    implicit none
    logical :: unitcell_can_change
!    if(unitcell_can_change() .and. sw_rebuild_pws==OFF) return
    icond = COORDINATE_CONTINUATION
  end subroutine continuation_mode

