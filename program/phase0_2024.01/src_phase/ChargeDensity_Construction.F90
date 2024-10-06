!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 633 $)
!
!  SUBROUINE: ChargeDensity_Construction, FermiEnergyLevel,
!             CD_Softpart_plus_Hardpart
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
subroutine ChargeDensity_Construction(ic)
! $Id: ChargeDensity_Construction.F90 633 2020-12-01 05:11:03Z jkoga $
  use m_Total_Energy,         only : m_TE_total_energy, ehartr
  use m_Charge_Density,       only : m_CD_convergence_check &
       &                           , m_CD_softpart, m_CD_hardpart &
       &                           , chgq_l &
       &                           , m_CD_hardpart_hsr
  use m_Electronic_Structure, only : metalic_system
  use m_ES_occup,             only : m_ESoc_fermi_parabolic, m_ESoc_fermi_tetrahedron &
       &                           , m_ESoc_fermi_ColdSmearing, m_ESoc_check_num_bands &
       &                           , m_ESoc_occup_fix &
       &                           , check_if_metalic, check_if_metalic_noncl
  use m_ES_IO,                only : m_ESIO_wd_EigenValues
  use m_Kpoints,              only : kv3
  use m_Files,                only : nfout
  use m_Crystal_Structure,    only : sw_bandgap_constraint, imag
  use m_Control_Parameters,   only : iprieigenvalue, projector_type, num_projectors &
       &                           , sw_hubbard, alpha_hubbard, icond, ekmode, occ_matrix_fix_period &
#ifdef ENABLE_ESM_PACK
       &                           , m_CtrlP_way_of_smearing,sw_esm, kimg, nspin
#else
       &                           , m_CtrlP_way_of_smearing
#endif
  use m_Const_Parameters,     only : PARABOLIC, MP, TETRAHEDRON, COLD, YES, OFF, ON &
       &                           , SPHERICAL_HARMONICS, ATOMIC_ORBITAL &
       &                           , INITIAL,CONTINUATION,FIXED_CHARGE,FIXED_CHARGE_CONTINUATION,CMPLDP
  use m_Orbital_Population,   only : m_OP_occ_mat_is_not_read

! ============================= added by K. Tagami =================== 5.0
  use m_Control_Parameters,  only : sw_eval_energy_before_charge, &
       &                            sw_update_charge_total, iprigap
! ==================================================================== 5.0

  use m_Parallelization,      only : ista_kngp,iend_kngp,MPI_CommGroup,mype

! =========================================== added by K. Tagami ========== 11.0
  use m_Control_Parameters,      only : noncol
  use m_Charge_Density,        only : m_CD_hardpart_hsr_noncl, &
       &                              m_CD_softpart_noncl, &
       &                              m_CD_hardpart_noncl
  use m_Total_Energy,          only : m_TE_total_energy_noncl
! ========================================================================= 11.0


  use m_PlaneWaveBasisSet,       only : kgp
  use m_FFT,                     only : fft_box_size_CD

  use m_IterationNumbers, only : iteration_electronic


! ============================= KT_add ================ 13.0E
  use m_Const_Parameters,       only : FERMI_Dirac
  use m_ES_Occup,               only : m_ESoc_fermi_Dirac
! ===================================================== 13.0E
  use m_Const_Parameters,       only : MP, FERMI_DIRAC, LOWEST_AT_EACH_KPT
  use m_ES_Occup,               only : m_ESoc_fermi_Dirac, m_ESoc_methfessel_paxton

! ============================= KT_add ================ 13.0U
  use m_CD_Mag_Moment,        only : m_CD_calc_ChgMagMom_in_sphere, &
       &                             m_CD_print_ChgMagmom_on_atom, &
       &                             sw_monitor_atomcharge
! ===================================================== 13.0U

! ======= KT_add ==== 13.0XX
  use m_Control_Parameters,     only : sw_calc_ekin_density, use_symm_ekin_density, &
       &                               use_asymm_ekin_density, ekin_density_is_active

  use m_KineticEnergy_Density,  only : m_KE_calc_symm_ekin_density, &
       &                               m_KE_calc_asymm_ekin_density, &
       &                               m_KE_cp_ekin_density_to_old
! =================== 13.0XX

! === KT_add === 2014/09/26
  use m_Crystal_Structure,  only : sw_fix_global_quantz_axis, sw_neglect_magmom
  use m_ES_Mag_Constraint,  only : m_ES_proj_magmom_G_quantz_axis, m_ES_set_magmom_zero
! ============== 2014/09/26

! === KT_add === 2014/08/26
  use m_Control_Parameters, only : sw_calc_orbital_moment, sw_use_add_proj, &
       &                          Orbital_decomp_mode
  use m_Charge_Density,    only : m_CD_hardpart_hsr_add, m_CD_hardpart_hsr_add_noncl
  use m_OP_Moment,  only :        m_OP_calc_OrbMagMom_in_sphere, &
       &                          m_OP_calc_OrbMagMom_method1, &
       &                          m_OP_calc_OrbMagMom_method2, &
       &                          m_OP_print_OrbMagMom_on_atom
! ============== 2014/08/26

! ==== Positron SCF === 2015/11/28
  use m_Control_Parameters,  only : sw_positron, positron_method
  use m_Const_Parameters,    only : positron_CONV
  use m_Positron_Wave_Functions, only : m_pWF_construct_pcharge, m_pWF_charge_rspace
! ===================== 2015/11/28
  use mpi

  implicit none
!  include 'mpif.h'

  integer, intent(in) :: ic
  logical             :: display_on, enough_bands
  complex(kind=CMPLDP),allocatable,dimension(:) :: vhar
  complex(kind=CMPLDP),allocatable, dimension(:,:) :: chgc
  integer :: ig,is,nfftcd,ierr, iflag
  integer :: nspin_m

#ifdef __TIMER_SUB__
  call timer_sta(701)
#endif

  if(ic == 0) display_on = .true.
  if(ic /= 0) display_on = .false.


  enough_bands = m_ESoc_check_num_bands()
  if(enough_bands) then
     call FermiEnergyLevel()            ! -(contained here)
     if ( iprigap > 1 ) then
        if ( noncol ) then
           call check_if_metalic_noncl( nfout, metalic_system )
        else
           call check_if_metalic( nfout, metalic_system )
        endif
     endif
  else
     if(iprieigenvalue >= 1) write(nfout,'(" the number of bands is not enough")')
     stop ' the number of bands is not enough'
  end if

! =================== added by K. Tagami ============ 5.0
  if ( sw_eval_energy_before_charge == ON ) then
! ================================= modified by K.Tagami ====== 11.0
!     call m_TE_total_energy(nfout,display_on,kv3)
!
     if ( noncol ) then
        call m_TE_total_energy_noncl(nfout,display_on,kv3)
     else
        call m_TE_total_energy(nfout,display_on,kv3)
     endif
! ============================================================= 11.0
  endif
! =================================================== 5.0

! ======================== modified by K. Tagami ======= 5.0
!  if(icond == INITIAL .or. icond == CONTINUATION) then
!     call CD_Softpart_plus_Hardpart()   ! -(contained here)
!     if(ic == 0) call m_CD_conversion_check(nfout)
!  else if ((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION).and.ekmode==OFF) then
!     call m_CD_hardpart_hsr(nfout,kv3)  ! fsr_l, fsi_l -> hsr
!  end if
  if ( sw_update_charge_total == ON ) then
     if ( icond == INITIAL .or. icond == CONTINUATION ) then
! ======================== modifided by K. Tagami ========================= 11.0
!       call CD_Softpart_plus_Hardpart()   ! -(contained here)
!
        if ( noncol ) then
           call CD_Softpart_plus_Hardpart_noncl()   ! -(contained here)
        else
           call CD_Softpart_plus_Hardpart()   ! -(contained here)
        endif
! ========================================================================= 11.0

! === KT_add == 2014/09/26
        if ( noncol ) then
           if ( sw_fix_global_quantz_axis == ON ) then
              call m_ES_proj_magmom_G_quantz_axis
           endif
           if ( sw_neglect_magmom == ON ) then
              call m_ES_set_magmom_zero
           endif
        endif
! ============= 2014/09/26

        if ( ic == 0 ) call m_CD_convergence_check(nfout)
     else if ( ( icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION ) &
          &         .and. ekmode==OFF ) then
! ======================== modifided by K. Tagami ========================= 11.0
!!     call m_CD_hardpart_hsr(nfout,kv3)             ! fsr_l, fsi_l -> hsr
        if ( noncol ) then
           call m_CD_hardpart_hsr_noncl(nfout,kv3)
        else
           call m_CD_hardpart_hsr(nfout,kv3)             ! fsr_l, fsi_l -> hsr
        endif
! ========================================================================== 11.0
     end if
     if( icond == INITIAL .or. icond == CONTINUATION .or. &
          & ((icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION ) .and. &
          & m_OP_occ_mat_is_not_read())) then
        call Renewal_of_OccMat( .false., OFF, iteration_electronic<=occ_matrix_fix_period ) ! evaluated with new om, hsr
     end if
  endif
! ========================================================= 5.0


#ifdef ENABLE_ESM_PACK
  if ( noncol ) then
     nspin_m = 1
  else
     nspin_m = nspin
  endif
  if(sw_esm==ON)then
     nfftcd = fft_box_size_CD(1,0)*fft_box_size_CD(2,0)*fft_box_size_CD(3,0)
     allocate(vhar(nfftcd));vhar=(0.d0,0.d0)
!     allocate(chgc(iend_kngp-ista_kngp+1,nspin));chgc=(0.d0,0.d0)
     allocate(chgc(1:kgp,nspin_m));chgc=(0.d0,0.d0)
     if(kimg==1)then
        do ig=ista_kngp,iend_kngp
!           chgc(ig-ista_kngp+1,1:nspin) = dcmplx(chgq_l(ig,1,1:nspin),0.d0)
           chgc(ig,1:nspin_m) = dcmplx(chgq_l(ig,1,1:nspin_m),0.d0)
        enddo
     else
        do ig=ista_kngp,iend_kngp
!           chgc(ig-ista_kngp+1,1:nspin) = dcmplx(chgq_l(ig,1,1:nspin),chgq_l(ig,2,1:nspin))
           chgc(ig,1:nspin_m) = dcmplx(chgq_l(ig,1,1:nspin_m),chgq_l(ig,2,1:nspin_m))
        enddo
     endif
     call mpi_allreduce(mpi_in_place,chgc,kgp*nspin_m,mpi_double_complex,mpi_sum,MPI_CommGroup,ierr)
     call esm_hartree(chgc,ehartr,vhar)
     ehartr  = 0.5d0*ehartr  !Ry -> Ha
     deallocate(chgc)
     deallocate(vhar)
  endif
#endif

! ======= KT_add ===== 13.0XX
!  if ( sw_calc_ekin_density == ON .and. ekin_density_is_active ) then
  if ( sw_calc_ekin_density == ON ) then
     if ( icond == INITIAL .or. icond == CONTINUATION ) then
        call m_KE_cp_ekin_density_to_old
        if ( use_symm_ekin_density )  call m_KE_calc_symm_ekin_density
        if ( use_asymm_ekin_density ) call m_KE_calc_asymm_ekin_density
     endif
  endif
! ==================== 13.0XX

! ============================= KT_add ================ 13.0U
  if ( sw_monitor_atomcharge == ON ) then
     call m_CD_calc_ChgMagMom_in_sphere
  endif
! ===================================================== 13.0U

! === KT_add === 2014/08/26
  if ( noncol .and. sw_calc_orbital_moment == ON ) then
     call m_OP_calc_OrbMagMom_in_sphere
     call m_OP_print_OrbMagMom_on_atom(nfout)

     if ( Orbital_decomp_mode == 1 ) then
        call m_OP_calc_OrbMagMom_method1
#if 0
     else if ( Orbital_decomp_mode == 2 ) then
        call m_OP_calc_OrbMagMom_method2
#endif
     endif
  endif
! ============== 2014/08/26

  if ( sw_positron /= OFF ) then
     if ( positron_method /= Positron_CONV ) then
        call m_pWF_construct_pcharge
        call m_pWF_charge_rspace
     endif
  endif

  if ( sw_eval_energy_before_charge == OFF ) then
! ================================== modified by K. Tagami ================ 11.0
     if ( noncol ) then
        call m_TE_total_energy_noncl(nfout,display_on,kv3)
     else
        call m_TE_total_energy(nfout,display_on,kv3)
     endif
! ========================================================================== 11.0
  endif
! ======================================================== 5.0

  if(ic == 0) call m_ESIO_wd_EigenValues(nfout,iprieigenvalue,nooccupation=YES)

#ifdef __TIMER_SUB__
  call timer_end(701)
#endif
contains
  subroutine FermiEnergyLevel()
    integer :: way_of_smearing
#ifdef __TIMER_SUB__
  call timer_sta(702)
#endif
    way_of_smearing = m_CtrlP_way_of_smearing()
    if(way_of_smearing == PARABOLIC) then
       call m_ESoc_fermi_parabolic(nfout)
!!$  else if(way_of_smearing == MP) then
!!$     call fermi_mesfessel_paxton(nfout)
    else if(way_of_smearing == TETRAHEDRON) then
        call m_ESoc_fermi_tetrahedron(nfout)
    else if(way_of_smearing == COLD) then
        call m_ESoc_fermi_ColdSmearing(nfout)

! ================ KT_add ========================= 13.0E
    else if(way_of_smearing == FERMI_DIRAC) then
        call m_ESoc_fermi_Dirac(nfout)
! ================================================= 13.0E
    else if(way_of_smearing == MP) then
        call m_ESoc_methfessel_paxton(nfout)
    else if(way_of_smearing == LOWEST_AT_EACH_KPT) then
       call m_ESoc_occup_fix(nfout)
    end if
    if(sw_bandgap_constraint == ON) then
       call m_ESoc_occup_fix(nfout)
    end if
#ifdef __TIMER_SUB__
  call timer_end(702)
#endif
  end subroutine FermiEnergyLevel

  subroutine CD_Softpart_plus_Hardpart
! $Id: ChargeDensity_Construction.F90 633 2020-12-01 05:11:03Z jkoga $
!fj#ifdef __TIMER_SUB__
!fj  call timer_sta(716)
!fj#endif
    call m_CD_softpart(nfout,kv3)
    call m_CD_hardpart(nfout,kv3)
!fj#ifdef __TIMER_SUB__
!fj  call timer_end(716)
!fj#endif
  end subroutine CD_Softpart_plus_Hardpart

! ============================ added by K. Tagami ============== 11.0
  subroutine CD_Softpart_plus_Hardpart_noncl
    call m_CD_softpart_noncl(nfout,kv3)
    call m_CD_hardpart_noncl(nfout,kv3)
  end subroutine CD_Softpart_plus_Hardpart_noncl
! ============================================================= 11.0

end subroutine ChargeDensity_Construction
