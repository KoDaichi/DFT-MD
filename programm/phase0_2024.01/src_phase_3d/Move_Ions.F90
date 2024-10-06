!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 599 $)
!
!  SUBROUINE:  Move_Ions, wd_cps_and_forces
!
!  AUTHOR(S): T. Yamasaki   August/20/2003
!
!  FURTHER MODIFICATION: T. Yamasaki   January/13/2004
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
subroutine Move_Ions
! $Id: Move_Ions.F90 599 2019-12-03 05:25:28Z ktagami $
  use m_Control_Parameters, only : iprimd, c_iteration2GDIIS &
       &                         , m_CtrlP_what_is_mdalg &
       &                         , m_CtrlP_set_gdiisoptmode &
       &                         , sw_charge_predictor,sw_wf_predictor,sw_rspace,af &
       &                         , sw_fcp, sw_betar_dot_wfs_exp, sw_precalculate_phase_vnonlocal &
       &                         , sw_optimize_coords_sametime &
       &                         , threshold_start_cellopt
  use m_Const_Parameters, only   : DP, TEMPERATURE_CONTROL, VERLET &
       &, BLUEMOON, QUENCHED_CONSTRAINT, QUENCHED_MD, NORMAL_MODE_ANALYSIS &
       &, HYPERPLANE_ADAPTIVE_COORDINATE, HAC, T_CONTROL, GDIIS &
       &, PHONON_FORCE, CG_STROPT, STEEPEST_DESCENT, BFGS, ON, CG_STROPT2 &
       &, PT_CONTROL, P_CONTROL, L_BFGS, FIRE
  use m_Files,        only : nfenf, nfdynm,nfout
  use m_Total_Energy, only : etotal
  use m_Ionic_System, only : natm,ekina,ega,m_IS_wd_cpo_and_forc, iteration_ionic_at_CNSTRA &
       &                   , m_IS_put_iteration_ionic_in_constraint &
       &                   , m_IS_md_thermo, m_IS_md_bluem, m_IS_md_cnstr &
       &                   , m_IS_gdiis, m_IS_wd_forc, m_IS_md, m_IS_cps_to_pos &
       &                   , m_IS_cp_cps2cpo,m_IS_wd_pos_and_v &
       &                   , m_IS_phonon_force, m_IS_cg, m_IS_cg2 &
       &                   , m_IS_evaluate_v_verlet &
       &                   , m_IS_update_cps_history &
       &                   , m_IS_force_af_symmetry &
       &                   , m_IS_fire
!!$       &                   , forcmx_constraint_quench, almda, mdmode &
  use m_Force,        only : forc_l, forcmx
  use m_IterationNumbers, only : iteration_ionic,iteration, iteration_unit_cell
  use m_Parallelization,    only : mype
#ifdef __TIMER__
  use m_Parallelization,    only : MPI_CommGroup
#endif

  use m_Realspace,          only : m_RS_resolve_mesh_soft
  use m_NonLocal_Potential, only : m_NLP_build_snl_in_rspace
  use m_PseudoPotential,    only : m_PP_alloc_radr,m_PP_dealloc_radr

! ================================ KT_add =================== 13.0B
  use m_Ionic_System,       only : m_IS_symmetrize_atom_pos
  use m_Control_Parameters,  only : sw_keep_symmetry_strict
! =========================================================== 13.0B

! === KT_add === 2014/06/10
  use m_Ionic_System,  only : sw_change_temperature_by_step, m_IS_reassgin_thermog
! ============== 2014/06/10

  use m_Fcp, only : m_Fcp_md_thermo, m_Fcp_cg, m_Fcp_cg2, m_Fcp_md, m_Fcp_gdiis &
  &               , m_Fcp_print_status

  use m_UnitCell, only : m_UnitCell_md
  use m_Control_Parameters, only :  sw_hybrid_functional, sw_rspace_hyb
  use m_RealSpace,    only : m_RS_resolve_mesh_hard, m_RS_build_qr_clm_ylm
  use m_FFT, only : fft_box_size_CD_exx,fft_box_size_CD_exx_nonpara

  use m_Ionic_System, only : m_IS_rigid_body_exists, m_IS_rb_dynamics

  use m_NonLocal_Potential, only : m_NLP_cal_i_l_exp_snl
  use m_ES_nonlocal,        only : m_ES_AtaulmnaG

  implicit none

  integer :: mdalg
  integer, save ::                optmode = QUENCHED_MD
  real(kind=DP),parameter :: c_forc2GDIIS = 0.0050d0
  integer, save ::                 iter_f = 0
  integer :: iteration_i

#ifdef __TIMER__
  integer :: ierr
  call mpi_barrier(MPI_CommGroup, ierr)
  call timer_sta(20)
#endif

  if(sw_charge_predictor==ON.or.sw_wf_predictor==ON) call m_IS_update_cps_history()

  mdalg = m_CtrlP_what_is_mdalg()

  if(mdalg == VERLET) call m_IS_evaluate_v_verlet(mdalg,forc_l)
  if(m_IS_rigid_body_exists() .and. (mdalg==VERLET.or.mdalg==QUENCHED_MD)) call m_IS_rb_dynamics(forc_l,mdalg)
  call wd_forces_cps_etotal_and_etc(mdalg,.true.) ! nfdynm  ! (Convergence_Check)
  if(mdalg /= T_CONTROL .and. mdalg /= BLUEMOON) then
     call wd_forces_cps_etotal_and_etc(mdalg,.false.) ! nfefn  ! (Convergence_Check)
  end if

!  if ( sw_optimize_coords_sametime == ON ) return
  if ( sw_optimize_coords_sametime == ON ) then
     if ( iteration_unit_cell > 1 ) then
        return
     else if ( forcmx < threshold_start_cellopt ) then
        return
     endif
  endif

! === KT_add === 2014/06/10
  if ( sw_change_temperature_by_step == ON ) then
     call m_IS_reassgin_thermog
  endif
! ============== 2014/06/10

  call m_IS_cp_cps2cpo()                 ! -> cpo_l
  mdalgorithm: select case(mdalg)
     case (GDIIS,BFGS,L_BFGS)
        if(iteration_ionic_at_CNSTRA<=0) then
           iteration_i = iteration_ionic
        else
           iteration_i = iteration_ionic - iteration_ionic_at_CNSTRA
        end if
        call m_IS_put_iteration_ionic_in_constraint(iteration_i)
!!$        optmode = m_CtrlP_set_gdiisoptmode(iteration_ionic,forcmx)
        optmode = m_CtrlP_set_gdiisoptmode(iteration_i,forcmx)
        if(optmode == QUENCHED_MD .or. optmode == STEEPEST_DESCENT) then
           call m_IS_md(optmode,forc_l)
           if(sw_fcp == ON ) call m_Fcp_md(optmode)
        else if(optmode == CG_STROPT) then
           call m_IS_cg(forc_l,etotal)
           if(sw_fcp == ON ) call m_Fcp_cg()
        else if(optmode == CG_STROPT2) then
           call m_IS_cg2(forc_l,etotal)
           if(sw_fcp == ON ) call m_Fcp_cg2()
        else
           call m_IS_gdiis(forc_l,forcmx,etotal)
           if(sw_fcp == ON ) call m_Fcp_gdiis()
!!$           if(mdmode == CNSTRA) call m_CtrlP_reset_optmode()
!!$     else if(mdalg == HAC) then
!!$        call me_nebm2(nfout,forc_l)
        end if
     case (CG_STROPT)
        call m_IS_cg(forc_l,etotal)
        if(sw_fcp == ON ) call m_Fcp_cg()
     case (CG_STROPT2)
        call m_IS_cg2(forc_l,etotal)
        if(sw_fcp == ON ) call m_Fcp_cg2()
     case (VERLET, QUENCHED_MD,STEEPEST_DESCENT)
        call m_IS_md(mdalg,forc_l)
        if(sw_fcp == ON ) call m_Fcp_md(mdalg)
     case (FIRE)
        call m_IS_fire(forc_l)
     case (T_CONTROL)
        call m_IS_md_thermo(forc_l)
        if(sw_fcp == ON ) call m_Fcp_md_thermo()
     case (BLUEMOON)
        call m_IS_md_bluem(forc_l)
     case (QUENCHED_CONSTRAINT)
        call m_IS_md_cnstr(forc_l)
     case (PHONON_FORCE)
        call m_IS_wd_forc(forc_l)
        call m_IS_phonon_force()
     case (PT_CONTROL,P_CONTROL)
        continue
     case default
        stop ' mdalg error at (Move_Ions)'
  end select mdalgorithm

  if(sw_fcp == ON) call m_Fcp_print_status()

  if(mdalg == T_CONTROL .or. mdalg == BLUEMOON) then
     call wd_forces_cps_etotal_and_etc(mdalg,.false.) ! nfefn  ! (Convergence_Check)
  end if

!!$  call wd_etotal_etc                           ! -(contained here)
!!$  call wd_cps_and_forces()                        ! -(contained here)
  call m_IS_cps_to_pos()

! ================================== KT_add  =================== 13.0B
  if ( sw_keep_symmetry_strict == ON ) then
     call m_IS_symmetrize_atom_pos(nfout)        ! -> cps,pos
  endif
! ============================================================== 13.0B

  if (af/=0) call m_IS_force_af_symmetry(nfout)

  if(iprimd >= 2) call m_IS_wd_forc(forc_l)

  if(sw_rspace==ON)then
     call m_PP_alloc_radr()
     call m_RS_resolve_mesh_soft(nfout)
     call m_NLP_build_snl_in_rspace(nfout)
     call m_PP_dealloc_radr()
  endif
  if (sw_hybrid_functional==ON .and. sw_rspace_hyb==ON ) then
     call m_PP_alloc_radr()
     call m_RS_resolve_mesh_hard(nfout,fft_box_size_CD_exx,fft_box_size_CD_exx_nonpara)
     call m_RS_build_qr_clm_ylm(deriv=.true.,box_size=fft_box_size_CD_exx)
     call m_PP_dealloc_radr()
  endif
  if(sw_betar_dot_wfs_exp==ON) then
    call m_NLP_cal_i_l_exp_snl()
  endif
  if(sw_precalculate_phase_vnonlocal==ON) then
    call m_ES_AtaulmnaG(hardpart=.true.)
  endif

#ifdef __TIMER__
  call timer_end(20)
#endif
contains
!!$  subroutine wd_cps_and_forces()
!!$
!!$    if(flag_wd_nfdynm == 0) then
!!$       call m_IS_wd_speciesname_etc(nfdynm)
!!$       flag_wd_nfdynm = 1
!!$    end if
!!$    if(mype == 0) then
!!$       if(mdmode == CNSTRA) then
!!$          write(nfdynm,'(" cps and forc at (iter_ion, iter_total = " &
!!$               & ,i5,i8," ) converged_in_a_plane")') iteration_ionic, iteration
!!$       else
!!$          write(nfdynm,'(" cps and forc at (iter_ion, iter_total = " &
!!$               & ,i5,i8," )")') iteration_ionic, iteration
!!$       end if
!!$    end if
!!$
!!$    call m_IS_wd_cpo_and_forc(nfdynm,forc_l)
!!$
!!$  end subroutine wd_cps_and_forces

  integer function set_optmode()
    set_optmode = optmode
    if(optmode == QUENCHED_MD .and. forcmx <= c_forc2GDIIS) then
       if(iter_f >= 3) then
          set_optmode = GDIIS
       else
          set_optmode = QUENCHED_MD
       end if
    end if
    iter_f = iter_f + 1
  end function set_optmode

end subroutine Move_Ions

