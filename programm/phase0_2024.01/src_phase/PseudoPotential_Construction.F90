!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 614 $)
!
!  SUBROUINE:  check_gncpp_type, PP_construction_paramset
!             PseudoPotential_Construction, PseudoPotential_ek
!
!  AUTHOR(S): T. Yamasaki and M. Okamoto   August/20/2003
!
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!
!
!
!=======================================================================
!
! =====================================================================
!   patch 10.1 by K. Tagami@adv    2011/06/18
!
!   patch 10.1:  allocation of vec_q_plus_G_LR etc for TDDFT
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
#   define FJ_TIMER_START(a)
#   define FJ_TIMER_STOP(a)
subroutine PseudoPotential_Construction
! $Id: PseudoPotential_Construction.F90 614 2020-05-07 03:24:24Z jkoga $
  use m_PseudoPotential,      only : m_PP_alloc0_ps_ntyp,   m_PP_alloc_ps_ntyp &
       &                           , m_PP_dealloc_ps_ntyp  &
       &                           , m_PP_set_mmesh,        m_PP_set_nloc &
#ifdef __EDA__
! -----  ascat starts modifying  -----
       &                           , m_PP_set_PP_per_atom_etc_zero &
! -----  ascat ceases modifying  -----
#endif
       &                           , m_PP_set_m_non0_lmtxlmt &
       &                           , m_PP_alloc_ps0,        m_PP_alloc_psc_qitg_rhpcg &
       &                           , m_PP_alloc_ps,         m_PP_check_gncpp_type &
! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
       &                           , m_PP_vanderbilt_type &
       &                           , m_PP_local_part,       m_PP_partial_core_CD &
! ==================================================================================================
       &                           , m_PP_wd_schgpc,        m_PP_make_index_lmtt_2_dl2p &
       &                           , m_PP_make_index_lmtt_phi &
       &                           , m_PP_make_index_lmtt_add &
       &                           , m_PP_make_index_lmtt_pao &
       &                           , m_PP_wd_variables,     m_PP_dealloc_ps &
       &                           , m_PP_rd_betar &
       &                           , m_PP_rd_PP_parameters &
       &                           , m_PP_make_qorb, iflag_ft &
       &                           , m_PP_cnstrct_crotylm &
       &                           , flg_paw, m_PP_rd_PAW_parameters &
       &                           , m_PP_cnstrct_crotylm_paw &
       &                           , m_PP_check_file_format
  use m_PseudoPotential,      only : nlmta,nlmt,lmta,ilmt
  use m_Parallelization,      only : m_Parallel_init_mpi_nlmta
  use m_PseudoPotential, only : ae_wavefunctions_are_detected
  use m_PseudoPotential, only : m_PP_alloc_qitg_wan, m_PP_set_dk_wan &
                              , m_PP_alloc_qitg_fef, m_PP_set_dk_fef
  use m_Wannier,         only : nwght,ghat
  use m_Crystal_Structure,  only :  rltv

  use m_Const_Parameters,     only : INITIAL,CONTINUATION &
       &                           , FIXED_CHARGE,FIXED_CHARGE_CONTINUATION &
       &                           , ON, OFF, by_pseudo_atomic_orbitals, ALL_AT_ONCE &
       &                           , pp_CIAOPP, pp_PAW
  use m_Kpoints,              only : kv3, vkxyz, mp_index
  use m_NonLocal_Potential,   only : m_NLP_alloc_snl  &
                                   , m_NLP_betar_dot_PWs_diff,  m_NLP_betar_dot_PWs &
       &                           , m_NLP_rd_snl &
       &                           , m_NLP_alloc_phig, m_NLP_phir_dot_PWs &
       &                           , m_NLP_alloc_snl_add, m_NLP_add_betar_dot_PWs &
       &                           , m_NLP_alloc_paog, m_NLP_paor_dot_PWs &
       &                           , m_NLP_build_snl_in_rspace, m_NLP_alloc_betar_optimized &
       &                           , m_NLP_cal_i_l_exp_snl
  use m_Ionic_System,         only : ityp, natm, m_IS_scale_magmom0_atomtyp
  use m_Ionic_System,         only : ntyp, iatomn, ivan
  use m_Control_Parameters,   only : paramset,icond,istress,evaluation_eko_diff &
       &                           , sw_orb_popu, sw_use_add_proj &
       &                           , sw_phonon, sw_calc_force, ipripp, ipriparallel &
       &                           , num_projectors, projector_type, intzaj &
       &                           , sw_wannier, sw_berry_phase, corecharge_cntnbin, sw_fef &
       &                           , ekmode, fixed_charge_k_parallel, continuation_using_ppdata &
       &                           , m_CtrlP_set_ppprinted,sw_rspace, sw_rspace_hyb, sw_hybrid_functional &
       &                           , sw_ldos, sw_rspace_ldos, m_CtrlP_rspace_integ_all_OK,nmax_G_hyb &
       &                           , sw_betar_dot_wfs_exp, sw_lband, sw_rspace_lband
  use m_PlaneWaveBasisSet,    only : kgp,gr_l,ngshell,ngshell_range
  use m_Const_Parameters,     only : OLD,ON, SPHERICAL_HARMONICS, NO
  use m_Files,                only : m_Files_open_ps_files,m_Files_close_ps_files,nfcntn_bin &
       &                           , m_Files_open_nfcntn_bin,nfpot,nfout &
       &                           , F_CNTN_BIN_in_partitioned &
       &                           , m_Files_open_nfcntn_bin_paw, nfcntn_bin_paw &
       &                           , m_Files_open_ps_file, m_Files_close_ps_file

  use m_Orbital_Population,   only : m_OP_ilmt_yy,m_OP_crotylm
  use m_FiniteElectricField,  only : numef,elec_id,m_FEF_Constract_of_ftq

  use m_Charge_Density,       only : m_CD_alloc_chgsoft
  use m_PAW_ChargeDensity,    only : m_PAWCD_init_surf_integration
  use m_PAW_XC_Potential,     only : m_PAW_XC_alloc_vxc
#ifdef FJ_TIMER
  use m_Parallelization,      only : MPI_CommGroup
#endif
! =================================== addded by K. Tagami =============== 5.0
  use m_Control_Parameters,   only : sw_mix_charge_hardpart
! ======================================================================= 5.0
  use m_Control_Parameters,   only : sw_hybrid_functional,m_CtrlP_set_hybrid_parameters
! ======================== Added by K. Tagami ================ 10.1
  use m_Control_Parameters,        only : sw_LinearResponse
  use m_PseudoPotential,           only : norm_q_plus_G_LR, &
       &                                  vec_q_plus_G_LR
!
  use m_LinearResponse_Qpt, only : m_LR_alloc_Vecs_q_plus_G, &
       &                           m_LR_set_Vecs_q_plus_G, &
       &                           m_LR_alloc_Qitg
! ============================================================ 10.1

! ============================== added by K. Tagami ============= 11.0
  use m_Control_Parameters,    only : noncol
  use m_PseudoPotential,    only : pot_has_soc, m_PP_check_spinorbit
! =============================================================== 11.0

! === KT_add === 2014/08/01
  use m_Orbital_QuantumNum,  only : m_OP_Qnum_init, m_OP_Qnum_read_orbital_index &
       &                          , m_OP_Qnum_orb_ind_data_num,m_OP_Qnum_read_orbital_index_it  &
       &                          , m_OP_Qnum_alloc_orb_ind_data
! ============== 2014/08/01

  use m_Realspace, only : m_RS_resolve_mesh_soft,m_RS_resolve_mesh_hard, m_RS_build_qr_clm_ylm

  use m_ES_ExactExchange, only : m_ES_EXX_set_nmax_G_hyb

  use m_FFT, only : fft_box_size_CD_exx,fft_box_size_CD_exx_nonpara

! ==== KT_add === 2014/09/19
  use m_Control_Parameters,  only : sw_calc_ekin_density, sw_rspace_ekin_density
! =============== 2014/09/19

! ==== ASMS_add === 2015/02/23
  use m_Control_Parameters,   only : sw_wannier90
  use m_PseudoPotential,      only : m_PP_set_dk_wan
  use m_Wannier90,           only : nntot, m_Wan90_set_dk_wan
! ================= 2015/02/23

#ifdef __EDA__
  use m_Control_Parameters, only : sw_eda
#endif

  implicit none

  integer       :: it, nfpp, is_gncpp , ierror
  integer  :: iprippd
#ifdef DEBUG_LDOS
  iprippd = 2
#else
  iprippd = ipripp
#endif

! ====================== KT_Test ========================= 12.5Exp
  if(sw_hybrid_functional == ON) then
     call m_CtrlP_set_hybrid_parameters
     call m_ES_EXX_set_nmax_G_hyb
  endif
! ======================================================== 12.5Exp

  call m_PP_alloc0_ps_ntyp()
  call m_PP_alloc_ps_ntyp(.true.)
  call PP_construction_paramset                !-(contained here)
  if(paramset) return
  call m_PP_dealloc_ps_ntyp
  call m_PP_alloc_ps_ntyp(.false.)

  call m_PP_alloc_ps0()
  call m_PP_alloc_psc_qitg_rhpcg()
  if(sw_wannier == ON) then
     call m_PP_alloc_qitg_wan(nwght)
     call m_PP_set_dk_wan(mp_index,ghat,rltv)
  end if
  if(sw_fef == ON) then
     call m_PP_alloc_qitg_fef(numef)
     call m_PP_set_dk_fef(mp_index,rltv,elec_id)
  end if
! ==== KT_add === 2015/02/23
  if ( sw_wannier90 == ON ) then
     call m_PP_alloc_qitg_wan( nntot )
     call m_Wan90_set_dk_wan
  endif
! =============== 2015/02/23

! ================================= Added by K. Tagami ======= 10.1
  if ( sw_LinearResponse == ON ) then
     call m_LR_alloc_Vecs_q_plus_G
     call m_LR_set_Vecs_q_plus_G
!     if ( modnrm == EXECUT ) then
        call m_LR_alloc_Qitg
!     endif
  endif
! ============================================================ 10.1

                                                  FJ_TIMER_START(29)
  if(sw_rspace/=ON .or. .not. m_CtrlP_rspace_integ_all_OK()) call m_NLP_alloc_snl()
                                                  FJ_TIMER_STOP(29)
  call m_NLP_alloc_snl_add()
  call m_NLP_alloc_phig()
  call m_NLP_alloc_paog()

  if(flg_paw .and. ekmode==ON .and. icond==FIXED_CHARGE) then
     call m_Files_open_nfcntn_bin_paw()
     call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
  end if

  if(sw_phonon == ON .and. sw_calc_force == OFF) return

  if(icond == INITIAL .or. icond == FIXED_CHARGE .or. continuation_using_ppdata == NO) then
!     call m_Files_open_ps_files(ivan,iatomn,ntyp,ierror)
!     if(ierror /= 0) call mpi_stop(nfout)
     call m_PP_alloc_ps()
     nfpp  = 0
#ifdef __EDA__
     if(sw_eda==ON) then
! -----  ascat starts modifying  -----
     call m_PP_set_PP_per_atom_etc_zero
! -----  ascat ceases modifying  -----
     endif
#endif
     do it = 1, ntyp
        call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
        if(ierror /= 0) call mpi_stop(nfout)
        if(iprippd >= 2) write(nfout,'(" !! PP_Construction: it = ",i3)') it
        if(ivan(it) /= OLD) then
           nfpp = nfpp + 1
           call m_PP_check_gncpp_type(nfpot(nfpp),nfout,is_gncpp)
           if(is_gncpp == pp_PAW) then
              ae_wavefunctions_are_detected(nfpp) = .true.
              if(iprippd >= 2) write(nfout,'(" !! AE wavefunction was detected.")')
           end if

! ============================== added by K. Tagami =================== 11.0
           call m_PP_check_spinorbit( nfpot(nfpp), nfout, pot_has_soc(it) )
! ===================================================================== 11.0

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
                                                  FJ_TIMER_START(29)
           call m_PP_vanderbilt_type(nfpot(nfpp),it,nfout,gr_l,ngshell,ngshell_range,.false.,is_gncpp)
                                                  FJ_TIMER_STOP(29)
! ==================================================================================================
        else if(ivan(it) == OLD) then
!!$        call pseudo_potential_norm_conserving_type(it,nfpkb,nfpd,nfpcc &
!!$            & ,iatomn(it), univol,nfout,ipri)
        endif
        call m_PP_local_part(it,nfout,gr_l,ngshell,ngshell_range)   ! (psloc), ->(vlocr,psc_l,etot1)
        call m_PP_partial_core_CD(nfout,it,gr_l,paramset=.false.)
        !                                      -(m_P.P.) epseud_pc, pcc
        call m_Files_close_ps_file(it)
     enddo

     call m_PP_wd_schgpc(nfout)  ! write down schgpc
     call m_PP_make_index_lmtt_2_dl2p(nfout,paramset)
     !           -(m_PseudoPotential) ->lmtt,il2p,lmta,nlmta,nlmtt,dl2p

     if(sw_rspace==ON)then
        !!call m_NLP_resolve_mesh_for_rspace(nfout)
        call m_RS_resolve_mesh_soft(nfout)
     endif
!     if(sw_ldos==ON .and. sw_rspace_ldos==ON)then
     if ( ( sw_ldos==ON .and. sw_rspace_ldos==ON ) &
          &   .or. ( sw_lband==ON .and. sw_rspace_lband==ON ) ) then
        call m_RS_resolve_mesh_hard(nfout)
        call m_RS_build_qr_clm_ylm(deriv=.false.)
     else if (sw_hybrid_functional==ON .and. sw_rspace_hyb==ON ) then
        call m_RS_resolve_mesh_hard(nfout,fft_box_size_CD_exx,fft_box_size_CD_exx_nonpara)
        call m_RS_build_qr_clm_ylm(deriv=.true.,box_size=fft_box_size_CD_exx)
! == KT_add == 2014/09/19
     else if ( sw_calc_ekin_density == ON .and. sw_rspace_ekin_density == ON ) then
        call m_RS_resolve_mesh_hard(nfout)
! =========== 2014/09/19
     endif

     if(sw_orb_popu == ON) then
        call m_PP_make_index_lmtt_phi(nfout,paramset)
     end if
     !           -(m_PseudoPotential) ->lmtt_phi,lmta_phi,nlmta_phi,nlmtt_phi
     if(sw_use_add_proj == ON) then
        call m_PP_make_index_lmtt_add(nfout,paramset)
     end if
     !           -(m_PseudoPotential) ->lmtt_add,lmta_add,nlmta_add,nlmtt_add
     if(intzaj == by_pseudo_atomic_orbitals) then
        call m_PP_make_index_lmtt_pao(nfout,paramset)
     end if
     !           -(m_PseudoPotential) ->lmtt_pao,lmta_pao,nlmta_pao,nlmtt_pao
     call m_PP_wd_variables(nfout)

! ========================= ASMS_Test =========================== 12.1
!     if(icond == INITIAL .or. (icond == FIXED_CHARGE .and. ekmode==OFF &
!!!$          & .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
!          & .and. fixed_charge_k_parallel == ALL_AT_ONCE) .or. &
!          & (continuation_using_ppdata == NO .and. .not.(icond == FIXED_CHARGE .and. ekmode==ON))) then
!
     if ( icond == INITIAL &
          & .or. ( ( icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION ) &
          &                   .and. ekmode==OFF &
          &                   .and. fixed_charge_k_parallel == ALL_AT_ONCE )&
          & .or. ( .not. ( (icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) &
          &                   .and. ekmode==ON ) &
          &                   .and. continuation_using_ppdata == NO ) ) then
! =============================================================== 12.1
        if(iprippd >= 2) then
           write(nfout,'(" icond  = ",i8)') icond
           write(nfout,'(" ekmode = ",i8)') ekmode
           write(nfout,'(" fixed_charge_k_parallel = ",i8)') fixed_charge_k_parallel
           write(nfout,'(" continuation_using_ppdata = ",i8)') continuation_using_ppdata
        end if
        if(istress == OFF) then
                                                  FJ_TIMER_START(29)
           if(sw_rspace == ON) then
              call m_NLP_alloc_betar_optimized()
              call m_NLP_build_snl_in_rspace(nfout)
           end if

           if(sw_rspace/=ON .or. .not. m_CtrlP_rspace_integ_all_OK()) then
              call m_NLP_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl
           end if
                                                  FJ_TIMER_STOP(29)
        else if(istress == ON) then
           call m_NLP_betar_dot_PWs_diff(nfout,kv3,vkxyz) !(kbint_diff) --> snl,snld
        else
           if(iprippd >= 1) write(nfout,'("  istress =",i4," <<PseudoPotential_Construction>>")') istress
           call phase_error_with_msg(nfout,'istress is not 0 or 1',__LINE__,__FILE__)
        endif
        if(iprippd >= 2) write(nfout,'(" sw_orb_popu = ",i8)') sw_orb_popu
        if(sw_orb_popu == ON) then
           if(iprippd >= 2) write(nfout,'(" m_PP_make_qorb")')
           call m_PP_make_qorb(nfout)
           call m_NLP_phir_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> phig
           call m_PP_cnstrct_crotylm(nfout) !-> crotylm
        end if
        if(sw_use_add_proj == ON) then
           call m_NLP_add_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl_add
        end if
        if(intzaj == by_pseudo_atomic_orbitals) then
           call m_NLP_paor_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> paog
        end if
        if(sw_betar_dot_wfs_exp == ON) then
           call m_NLP_cal_i_l_exp_snl()
        endif
        if(num_projectors > 0) then
           call m_OP_ilmt_yy
           call m_OP_crotylm(nfout)
        end if
     end if
!     if(icond==CONTINUATION.and.flg_paw) then
!         call m_Files_open_nfcntn_bin_paw()
!         call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
!     end if
! ============= ASMS_DEBUG =========================== 12.1
!     if(icond==FIXED_CHARGE_CONTINUATION.and.flg_paw) then
!         call m_Files_open_nfcntn_bin_paw()
!         call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
!     end if
! ============= ASMS_DEBUG =========================== 12.1
!
! ====== KT_mod ==== 13.0S
!     if( icond == FIXED_CHARGE .and. ekmode==ON ) then
!
     if( ( icond == FIXED_CHARGE .or. icond == FIXED_CHARGE_CONTINUATION ) &
          &  .and. ekmode==ON ) then
! =================13.0S

        call m_NLP_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl
        if(sw_orb_popu == ON) then
           call m_PP_make_qorb(nfout)
           call m_NLP_phir_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> phig
           call m_PP_cnstrct_crotylm(nfout) !-> crotylm
        end if

        if(sw_use_add_proj == ON) then
           call m_NLP_add_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl_add
        end if
        if(intzaj == by_pseudo_atomic_orbitals) then
           call m_NLP_paor_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> paog
        end if

        if(num_projectors > 0) then
           call m_OP_ilmt_yy
           call m_OP_crotylm(nfout)
        end if
     endif
!
     call m_PP_dealloc_ps
!     call m_Files_close_ps_files
  else if(icond == CONTINUATION .or. icond == FIXED_CHARGE_CONTINUATION) then
!!$  else if(icond == CONTINUATION .or. (icond == FIXED_CHARGE_CONTINUATION .and.ekmode==OFF &
!!$       & .and. fixed_charge_k_parallel == ALL_AT_ONCE)  &
!!$       & .or. (icond == FIXED_CHARGE_CONTINUATION .and. ekmode == ON)) then
     call m_Files_open_nfcntn_bin()
!!$     if(istress == ON) call m_PP_alloc_ps_stress()
     call m_PP_rd_PP_parameters(nfout,nfcntn_bin,F_CNTN_BIN_in_partitioned, kgp)
     if(iprippd >= 2) write(nfout,'(" ! PP_parameters were read!")')
     if(icond == CONTINUATION) &
          & call m_NLP_rd_snl(nfout,nfcntn_bin,F_CNTN_BIN_in_partitioned,kv3)
     if(icond == FIXED_CHARGE_CONTINUATION ) call m_PP_rd_betar(nfcntn_bin) !->betar

!     if(flg_paw) then
!         call m_Files_open_nfcntn_bin_paw()
!         call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
!     end if

     if(sw_orb_popu == ON .or. sw_use_add_proj == ON .or. num_projectors>0 &
          & .or. sw_berry_phase == ON .or. sw_wannier == ON &
          & .or. sw_fef == ON &
          & .or. corecharge_cntnbin == OFF) then
        iflag_ft = OFF
        if(iprippd >= 2) then
           write(nfout,'(" !pp sw_orb_popu     = ",i8)') sw_orb_popu
           write(nfout,'(" !pp sw_use_add_proj = ",i8)') sw_use_add_proj
           write(nfout,'(" !pp num_projectors  = ",i8)') num_projectors
           write(nfout,'(" !pp sw_berry_phase  = ",i8)') sw_berry_phase
           write(nfout,'(" !pp sw_wannier      = ",i8)') sw_wannier
           write(nfout,'(" !pp sw_fef          = ",i8)') sw_fef
           write(nfout,'(" !pp corecharge_cntnbin = ",i8)') corecharge_cntnbin
        end if
!        call m_Files_open_ps_files(ivan,iatomn,ntyp,ierror)
!        if(ierror/=0) call mpi_stop(nfout)
        call m_PP_alloc_ps()
        nfpp  = 0
        do it = 1, ntyp
           call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
           if(ierror/=0) call mpi_stop(nfout)
           if(iprippd >= 2) write(nfout,'(" !! PP_Construction for PDOS: it = ",i3)') it
           if(ivan(it) /= OLD) then
              nfpp = nfpp + 1
              call m_PP_check_gncpp_type(nfpot(nfpp),nfout,is_gncpp)
              if(is_gncpp == -2) then
                 !!$is_gncpp = 2
                 ae_wavefunctions_are_detected(nfpp) = .true.
!!$                 if(ipripp >= 2 ) write(nfout,'(" !! AE wavefunction was detected.")')
              end if

! ============================== added by K. Tagami =================== 11.0
              call m_PP_check_spinorbit( nfpot(nfpp), nfout, pot_has_soc(it) )
! ===================================================================== 11.0

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
              call m_PP_vanderbilt_type(nfpot(nfpp),it,nfout,gr_l,ngshell,ngshell_range,.false.,is_gncpp)
! ==================================================================================================
           else
           end if
!!$ASASASASAS
           call m_PP_partial_core_CD(nfout,it,gr_l,.true.)
!!$ASASASASAS
           call m_Files_close_ps_file(it)
        enddo

        if(sw_orb_popu == ON) then
           call m_PP_make_index_lmtt_phi(nfout,paramset)
           call m_PP_make_qorb(nfout)
           call m_NLP_phir_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> phig
           call m_PP_cnstrct_crotylm(nfout) !-> crotylm
        end if
        if(sw_use_add_proj == ON) then
           call m_PP_make_index_lmtt_add(nfout,paramset)
           call m_NLP_add_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl_add
        end if
        if(num_projectors > 0) then
           call m_OP_ilmt_yy
           call m_OP_crotylm(nfout)
        end if
        call m_PP_dealloc_ps
!        call m_Files_close_ps_files
     end if
  end if

! ==== DEBUG TEST ==== 2020/01/28
  if ( flg_paw ) then
     if ( icond == CONTINUATION .or. icond == FIXED_CHARGE &
          &       .or. icond == FIXED_CHARGE_CONTINUATION ) then
        call m_Files_open_nfcntn_bin_paw()
        call m_PP_rd_PAW_parameters(nfout,nfcntn_bin_paw)
     endif
  end if
! ==== DEBUG TEST ==== 2020/01/28

  if(flg_paw) then
     call m_PAWCD_init_surf_integration
     call m_PAW_XC_alloc_vxc

     call m_CD_alloc_chgsoft
     call m_PP_cnstrct_crotylm_paw
! =============================== modified by K. Tagami ================ 5.0
!  end if

  else if ( sw_mix_charge_hardpart == ON ) then
     call m_PAWCD_init_surf_integration
     call m_PP_cnstrct_crotylm_paw
  end if
! ====================================================================== 5.0

  call m_Parallel_init_mpi_nlmta(nfout,ipriparallel &
       & ,nlmta,nlmt,natm,lmta,ntyp,ilmt,ityp) ! -> np_fsri
!!$#endif

! ==== KT_add ==== 2014/08/01
!  call m_Files_open_ps_files(ivan,iatomn,ntyp,ierror)
!  if(ierror /= 0) call mpi_stop(nfout)
!  call m_OP_Qnum_init
!  call m_OP_Qnum_read_orbital_index
!  call m_Files_close_ps_files

  call m_OP_Qnum_init
  do it=1,ntyp
     call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
     call m_OP_Qnum_orb_ind_data_num(it)
     call m_Files_close_ps_file(it)
  enddo
  call m_OP_Qnum_alloc_orb_ind_data()
  do it=1,ntyp
     call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
     call m_OP_Qnum_orb_ind_data_num(it)
     call m_OP_Qnum_read_orbital_index_it(it)
     call m_Files_close_ps_file(it)
  enddo
! ================ 2014/08/01

  call m_IS_scale_magmom0_atomtyp

contains
  subroutine PP_construction_paramset
!    call m_Files_open_ps_files(ivan, iatomn, ntyp,ierror)
!    if(ierror/=0) call mpi_stop(nfout)
    nfpp = 0
    do it = 1, ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
       if(ierror/=0) call mpi_stop(nfout)
       if(ipripp >= 2) write(nfout,'(" !!PP_construction_paramset: it = ",i3)') it
       if(ivan(it) /= OLD) then
          nfpp = nfpp + 1
           call m_PP_check_gncpp_type(nfpot(nfpp),nfout,is_gncpp)
           if(is_gncpp == -2) then
             !!$is_gncpp = 2
             ae_wavefunctions_are_detected(nfpp) = .true.
             if(ipripp >= 2) write(nfout,'(" !! AE wavefunction was detected.")')
           end if

! ============================== added by K. Tagami =================== 11.0
              call m_PP_check_spinorbit( nfpot(nfpp), nfout, pot_has_soc(it) )
! ===================================================================== 11.0

! === m_PP_vanderbilt_type and m_PP_vanderbilt_type_3D look the same!!! by tkato 2011/08/04 ========
                                                  FJ_TIMER_START(29)
           call m_PP_vanderbilt_type(nfpot(nfpp),it,nfout,gr_l,ngshell,ngshell_range,.true.,is_gncpp)
                                                  FJ_TIMER_STOP(29)
! ==================================================================================================
       else if(ivan(it) == OLD) then
!!$        call pseudo_potential_norm_conserving_type(it,nfpkb,nfpd,nfpcc &
!!$             & ,iatomn(it), univol,nfout,ipri)
       endif
           call m_PP_partial_core_CD(nfout,it,gr_l,.true.)
       ! epseud_pc, pcc
       call m_Files_close_ps_file(it)
    enddo
    call m_PP_make_index_lmtt_2_dl2p(nfout,.true.) ! ->nlmt,nlmta,nlmtt
    if(sw_orb_popu == ON) then
       call m_PP_make_index_lmtt_phi(nfout,.true.)    ! ->nlmt_phi,nlmta_phi,nlmtt_phi
    end if
    if(sw_use_add_proj == ON) then
       call m_PP_make_index_lmtt_add(nfout,.true.) ! ->lmtt_add,lmta_add,nlmta_add,nlmtt_add
    end if
    if(intzaj == by_pseudo_atomic_orbitals) then
       call m_PP_make_index_lmtt_pao(nfout,.true.) ! ->lmtt_pao,lmta_pao,nlmta_pao,nlmtt_pao
    end if

    call m_PP_set_mmesh
    call m_PP_set_nloc(nfout)
    call m_PP_set_m_non0_lmtxlmt
!    call m_Files_close_ps_files
    call m_CtrlP_set_ppprinted(ON)
  end subroutine PP_construction_paramset

end subroutine PseudoPotential_Construction

subroutine Check_of_Pseudopotential
  use m_Const_Parameters,     only : DP, OLD, ON
  use m_Files,                only : nfout,nfpot &
       &                           , m_Files_open_ps_file, m_Files_close_ps_file
#ifndef ENABLE_ESM_PACK
  use m_Control_Parameters,   only : ipripp &
#else
  use m_Control_Parameters,   only : ipripp, esm_qbac &
#endif
      & , neg,m_CtrlP_set_neg_properly, m_CntrlP_set_neg, m_CntrlP_set_meg
  use m_Crystal_Structure,    only : additional_charge
  use m_Ionic_System,         only : ntyp, iatom, iatomn, ivan, zeta1, qex
  use m_PseudoPotential,      only : m_PP_rd_ival
  use m_Parallelization,      only : MPI_CommGroup
  use m_Ionic_System,         only : natm,ionic_charge_atomtyp, ionic_charge_atoms &
       &                           , mag_moment0_atoms_is_defined, ityp, iwei &
       &                           , sw_set_initial_magmom_by_atom
  use mpi
  implicit none
!  include 'mpif.h'                                      ! MPI
  integer :: it, ia, ierror, nfpp
  real(kind=DP),allocatable,dimension(:)   :: ivalt  ! d(ntyp) #valence electrons
  real(kind=DP) :: totch_t, ival
  allocate(ivalt(ntyp))
  nfpp = 0
  do it = 1, ntyp
     call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
     if(ierror /= 0) call mpi_stop(nfout)
     if(ipripp >= 2) write(nfout,'(" !! PP_Construction: it = ",i3)') it
     if(ivan(it) /= OLD) then
        nfpp = nfpp + 1
        call m_PP_rd_ival(nfpot(nfpp),it,nfout,ival) ! -> ival
     else if(ivan(it) == OLD) then
     endif

     call mpi_bcast(ival,1,mpi_double_precision,0,MPI_CommGroup,ierror)
     ivalt(it) = ival
     call m_Files_close_ps_file(it)
  end do
  totch_t = 0.d0
  do it = 1, ntyp
     totch_t = totch_t + ivalt(it)*iatom(it) + qex(it)
     write(nfout,'(" ## it = ",i5,", ivalt = ",f8.4, " iatom = ",f8.4," qex = ",f8.4)') it, ivalt(it),iatom(it),qex(it)
  end do
  call mpi_bcast(totch_t,1,mpi_double_precision,0,MPI_CommGroup,ierror)
#ifdef ENABLE_ESM_PACK
  totch_t = totch_t - esm_qbac
#endif
! ===== KT_add === 2014/06/08
  totch_t = totch_t - additional_charge      ! totch is num. of electrons
! ================ 2014/06/08

  if ( sw_set_initial_magmom_by_atom == ON ) then
     if ( mag_moment0_atoms_is_defined ) then
        do ia=1, natm
           totch_t = totch_t -ionic_charge_atoms(ia) *iwei(ia)
        end do
     else
        Do ia=1, natm
           it = ityp(ia)
           totch_t = totch_t -ionic_charge_atomtyp(it) *iwei(ia)
        end do
     endif
  else
     Do ia=1, natm
        it = ityp(ia)
        totch_t = totch_t -ionic_charge_atomtyp(it) *iwei(ia)
     end do
  endif
  
  if(totch_t >= neg*2.0) then
     if(ipripp >= 0) then
        write(nfout,'(" ### Warning(1309): Number of bands(neg) is insufficient:")')
        write(nfout,'("                totch = ",f10.3," >= neg*2.0 = ",f10.3)') totch_t, neg*2.0
     end if
     if(dabs(sum(zeta1)) > 0.d0) then
        call m_CtrlP_set_neg_properly(1.3*totch_t) ! -> neg
     else
        call m_CtrlP_set_neg_properly(totch_t) ! -> neg
     end if
     call m_CntrlP_set_meg(neg)
!!$     call m_CntrlP_set_neg(neg)
     if(ipripp >= 0) then
        write(nfout,'(" ### Reset value of neg = ",i8)') neg
        write(nfout,'("                totch = ",f10.3," <= neg*2.0 = ",f10.3)') totch_t, neg*2.0
     end if
  else
     if(ipripp >= 0) then
        write(nfout,'(" ###  Number of bands(neg) is sufficient:")')
        write(nfout,'("                totch = ",f10.3," < neg*2.0 = ",f10.3)') totch_t, neg*2.0
     end if
  end if

  deallocate(ivalt)
end subroutine Check_of_Pseudopotential

subroutine PseudoPotential_ek
  use m_Const_Parameters,     only : FIXED_CHARGE, FIXED_CHARGE_CONTINUATION, ON, OFF, &
       &                             by_pseudo_atomic_orbitals
  use m_Control_Parameters,   only : icond, sw_orb_popu, sw_pdos, sw_use_add_proj, uvsormode,sw_rspace &
       &                            ,m_CtrlP_rspace_integ_all_OK, intzaj

  use m_Kpoints,              only : kv3, vkxyz, kv3_previous
  use m_NonLocal_Potential,   only : m_NLP_rd_snl, m_NLP_betar_dot_PWs &
       &                           , m_NLP_betar_dot_PWs_diff &
       &                           , m_NLP_phir_dot_PWs &
       &                           , m_NLP_add_betar_dot_PWs &
       &                           , m_NLP_build_snl_in_rspace, m_NLP_paor_dot_PWs
  use m_PseudoPotential,      only : m_PP_alloc_NLP, m_PP_dealloc_NLP, m_PP_betar_calculated
  use m_IterationNumbers,     only : iteration_electronic,nk_in_the_process
  use m_Files,                only : nfout, nfcntn_bin, m_Files_open_nfcntn_bin &
       &                           , F_CNTN_BIN_in_partitioned
! ======================================= Added by K. Tagami ========= 0.2
!!$  use m_Control_Parameters,   only : num_projectors, projector_type, paramset &
!!$       &                           , intzaj, by_pseudo_atomic_orbitals
!!$  use m_Orbital_Population,   only : m_OP_ilmt_yy,m_OP_crotylm
!!$  use m_PseudoPotential,      only : m_PP_make_qorb, m_PP_cnstrct_crotylm &
!!$       &                           , m_PP_make_index_lmtt_phi &
!!$       &                           , m_PP_make_index_lmtt_add &
!!$       &                           , m_PP_make_index_lmtt_pao &
!!$       &                           , m_PP_alloc_ps, m_PP_dealloc_ps
! =====================================================================
  implicit none

  if(icond == FIXED_CHARGE &
       & .or. (icond == FIXED_CHARGE_CONTINUATION &
       &           .and.( iteration_electronic == 0 .or. kv3_previous /= kv3))) then
     call m_PP_alloc_NLP()

     if ( intzaj == by_pseudo_atomic_orbitals ) then
        call m_NLP_paor_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> paog
     end if
     if(sw_rspace == ON) then
        call m_NLP_build_snl_in_rspace(nfout)
! ----------- Revised by T.Yamasaki, 1 Aug. 2014 ---------->>
     end if
     if(sw_rspace/=ON .or. .not. m_CtrlP_rspace_integ_all_OK()) then
        call m_NLP_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl
! <<--------------------------------------------------------
     endif
     if (sw_orb_popu == ON ) then
        call m_NLP_phir_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> phig
     end if
     if(sw_use_add_proj == ON) then
        call m_NLP_add_betar_dot_PWs(nfout,kv3,vkxyz) !(kbint) --> snl_add
     end if
     if(uvsormode == OFF) call m_PP_dealloc_NLP()
!!$     if(nk_in_the_process == 3) stop 'PseudoPotential_ek'
  else if(icond == FIXED_CHARGE_CONTINUATION) then
!!$     call m_Files_open_nfcntn_bin()
     call m_NLP_rd_snl(nfout,nfcntn_bin,F_CNTN_BIN_in_partitioned,kv3)
  end if
end subroutine PseudoPotential_ek

logical function has_uspp()
  use m_Files,                only : nfout, nfpot, m_Files_open_ps_file, m_Files_close_ps_file
  use m_Ionic_System,         only : ntyp, iatomn, ivan
  use m_PseudoPotential,      only : m_PP_alloc0_ps_ntyp, m_PP_alloc_ps_ntyp, m_PP_dealloc &
     &                             , m_PP_vanderbilt_type,m_PP_include_vanderbilt_pot      &
     &                             , m_PP_check_gncpp_type, ae_wavefunctions_are_detected
  use m_Parallelization,      only : ista_kngp, iend_kngp
  use m_Const_Parameters,     only : YES, DP
  use mpi
  implicit none
!  include 'mpif.h'                                      ! MPI
  integer :: nfpp,it,iret,istabuf,iendbuf,ierror,is_gncpp
  integer, allocatable, dimension(:,:) :: ngsh
  real(kind=DP), allocatable, dimension(:) :: gr
  call m_PP_alloc0_ps_ntyp()
  call m_PP_alloc_ps_ntyp(.true.)
  has_uspp = .false.
  nfpp = 0
  istabuf = ista_kngp;iendbuf = iend_kngp
  ista_kngp = 1;iend_kngp=1
  allocate(gr(1))
  allocate(ngsh(1,1))
  do it = 1, ntyp
     call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierror)
     if(ierror/=0) call mpi_stop(nfout)
     nfpp = nfpp + 1
     call m_PP_check_gncpp_type(nfpot(nfpp),nfout,is_gncpp)
     if(is_gncpp == -2) then
        ae_wavefunctions_are_detected(nfpp) = .true.
     end if
     call m_PP_vanderbilt_type(nfpot(nfpp),it,nfout,gr,1,ngsh,.true.,is_gncpp,.true.)
     call m_Files_close_ps_file(it)
     iret = m_PP_include_vanderbilt_pot(it)
     if(iret==YES) then
       has_uspp = .true.
     endif
  enddo
  deallocate(gr)
  deallocate(ngsh)
  ista_kngp = istabuf
  iend_kngp = iendbuf
  call m_PP_dealloc()
end function has_uspp

