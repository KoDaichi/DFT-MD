!=======================================================================
!
!  PROGRAM  PHASE/0 2016.01 ($Rev: 603 $)
!
!  SUBROUINE:  Preparation_for_mpi
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
subroutine Preparation_for_mpi(prepare_communicators)
! $Id: Preparation_for_mpi.F90 603 2020-04-07 03:25:30Z jkoga $
!                           @(#)Preparation_for_mpi.F90 1.10 03/02/19 00:49:14
  use m_Const_Parameters,     only : INITIAL, CONTINUATION &
       &                           , FIXED_CHARGE, FIXED_CHARGE_CONTINUATION,OFF &
       &                           , ALL_AT_ONCE, ON, DP, COORDINATE_CONTINUATION &
       &                           , ONE_BY_ONE
  use m_Kpoints,              only : kv3, kv3_ek
  use m_Files,                only : nfout
  use m_Control_Parameters,   only : ipriparallel,nspin,neg,printable,ngnode_nbmx &
       &                           , flag_mpi_g_dot_r,flag_mpi_g_dot_r_k &
       &                           , icond, ekmode, fixed_charge_k_parallel, sw_rsb, neg_is_given &
       &                           , m_CtrlP_flag_mpi_G_dot_R, fftbox_divide_cube &
       &                           , fftbox_3ddiv_1, fftbox_3ddiv_2, fftbox_3ddiv_3 &
       &                           , fftbox_div_1, fftbox_div_2, sw_fft_xzy &
       &                           , sw_communicator_for_chg
#ifdef KMATH_FFT3D
  use m_Control_Parameters,   only : nproc_fft3d
#endif

  use m_PlaneWaveBasisSet,    only : kg1, kgpm, nbmx, iba, kgp, kgp_reduced
#ifndef PARAMSET
  use m_Kpoints,              only : k_symmetry
  use m_Parallelization,      only : m_Parallel_init_mpi_elec_3D, m_Parallel_init_mpi_iba_3D &
       &                           , make_index_band_3D, m_Parallel_mpi_fft_box &
       &                           , make_index_band_for_Gdiv_3D, m_Parallel_mpi_fft_box_cd &
       &                           , m_Parallel_mpi_fft_box_3div, m_Parallel_mpi_fft_box_cd_3div &
       &                           , m_Parallel_mpi_fft_box_xyz , m_Parallel_mpi_fft_box_cd_xyz  &
       &                           , m_Parallel_init_mpi_kngp_3D, m_Parallel_init_mpi_atm &
       &                           , m_Parallel_init_mpi_atm2 &
       &                           , m_Parallel_init_mpi_atm_f &
       &                           , m_Parallel_init_mpi_kngp_B_3D, m_Parallel_init_mpi_atm_B_3D &
       &                           , m_Parallel_init_mpi_atm_ke &
       &                           , m_Parallel_chgq_onto_fftcd_3D, m_Parallel_fftcd_onto_chgq_3D &
       &                           , m_Parallel_fft_onto_chgq_3D ,m_Parallel_init_mpi_mix &
       &                           , m_Parallel_wf_onto_fft_3D, m_Parallel_fft_onto_wf_3D &
       &                           , nel_fft_x, nel_fft_y, nel_fft_z &
       &                           , m_Parallel_init_mpi_gga     &
       &                           , m_Parallel_init_mpi_nbmx    &
       &                           , m_Parallel_init_mpi_snl_3D  &
       &                           , m_Parallel_init_mpi_ffth &
       &                           , m_Parallel_init_mpi_kv3_ek
#ifdef MPI_FFTW
  use m_Parallelization,      only : m_Parallel_wf_onto_fft_mpifftw, m_Parallel_fft_onto_wf_mpifftw &
       &                           , m_Parallel_fft_onto_chgq_mpifftw
#endif
  use m_Control_Parameters, only :   nblocksize_mgs          &
       &                           , nblocksize_mgs_is_given , kimg, GAMMA
  use m_Electronic_Structure,only:  nblocksize_mgs_default
  use m_FFT,                  only : nfft &
       &                           , fft_box_size_WF               &
       &                           , fft_box_size_CD_3D            &
       &                           , m_FFT_Direct_3D,  m_FFT_Direct_XYZ_3D  &
       &                           , m_FFT_Inverse_3D, m_FFT_Inverse_XYZ_3D  &
       &                           , m_FFT_Inverse_XYZ_3D_oo_place &
#ifdef MPI_FFTW
       &                           , nfftp, nfftps, m_FFT_init_mpifftw
#else
       &                           , nfftp, nfftps
#endif
  use m_PlaneWaveBasisSet,    only : m_pwBS_alloc_ngpt_igfp_gr_3D  &
       &                           , m_pwBS_calc_length_of_G_3D    &
       &                           , m_pwBS_G_trans_functions_3D   &
       &                           , m_pwBS_setup_FFTmapfunctions_3D &
       &                           , m_pwBS_set_ngabc_B_3D       &
       &                           , m_pwBS_set_ngabc_3D         &
       &                           , igfp_l                       &
       &                           , igf, nbase, nbase_gamma, kg, kg_gamma
  use m_Charge_Density,       only : m_CD_alloc_chgq
  use m_XC_Potential,         only : m_XC_alloc_vxc_3D
  use m_Ionic_System,         only : m_IS_alloc_zfm3_3D
  use m_Ionic_System,         only : natm, natm2, m_IS_alloc_fxyzew
  use m_Force,                only : m_Force_alloc
#ifdef __EDA__
  use m_XC_Potential_2D,      only : m_XC_alloc_exc_on_a_grid
  use m_Ionic_System,         only : m_IS_alloc_eewald_per_atom, m_IS_alloc_zfm3_EDA
  use m_PseudoPotential,      only : m_PP_alloc_PP_per_atom_etc
#endif
#endif

! ============================== added by K. Tagami ================== 11.0&13.0XX
  use m_Control_Parameters,   only : noncol, ndim_spinor, sw_calc_ekin_density
  use m_KineticEnergy_Density,  only : m_KE_alloc_ekin_density
! ==================================================================== 11.0&13.0XX

#ifdef __EDA__
  use m_Control_Parameters,   only : sw_eda
#endif

#ifdef KMATH_FFT3D
  use m_Parallelization, only : m_Parallel_kmath3d_init
  use m_Control_Parameters, only : sw_kmath_fft3d, nstage_fft3d
#endif

#ifdef MPI_FFTW
  use m_Control_Parameters,only: sw_mpi_fftw
#endif

  implicit none

  integer, intent(in) :: prepare_communicators
  integer             :: lsize
  real(kind=DP), allocatable, dimension(:,:) :: dfft_l

#ifndef PARAMSET
  call m_IS_alloc_fxyzew()
  call m_Force_alloc()

! === KT_add ==== 13.0XX
  if ( sw_calc_ekin_density == ON ) call m_KE_alloc_ekin_density
! =============== 13.0XX

  if(prepare_communicators==ON)then

!!$     if(neg_is_given) then
        call m_Parallel_init_mpi_elec_3D(nfout,ipriparallel,printable,neg,kv3,nspin,kg1,iba)
        call make_index_band_3D(nfout,ipriparallel,printable,kv3,neg &
             & , nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
        call make_index_band_for_Gdiv_3D(neg, nblocksize_mgs,nblocksize_mgs_is_given,nblocksize_mgs_default)
!!$     end if
#ifdef FFT_3D_DIVISION
     call m_Parallel_mpi_fft_box_3div(nfout,ipriparallel,printable,fft_box_size_WF,kimg, &
    &                            fftbox_3ddiv_1, fftbox_3ddiv_2,fftbox_3ddiv_3)
#else
     if (sw_fft_xzy > 0) then
        call m_Parallel_mpi_fft_box(nfout,ipriparallel,printable,fft_box_size_WF,kimg,fftbox_divide_cube)
     else
        call m_Parallel_mpi_fft_box_xyz(nfout,ipriparallel,printable,fft_box_size_WF,kimg, &
    &                            fftbox_div_1, fftbox_div_2)
#ifdef KMATH_FFT3D
        if(sw_kmath_fft3d==ON) call m_Parallel_kmath3d_init(fft_box_size_WF,nstage_fft3d)
#endif
     end if
#endif
#ifdef FFT_3D_DIVISION_CD
     call m_Parallel_mpi_fft_box_cd_3div(nfout,ipriparallel,printable,fft_box_size_CD_3D,kimg, &
    &                               fftbox_3ddiv_1,fftbox_3ddiv_2,fftbox_3ddiv_3)
#else
     if (sw_fft_xzy > 0) then
        call m_Parallel_mpi_fft_box_cd(nfout,ipriparallel,printable,fft_box_size_CD_3D,kimg,fftbox_divide_cube)
     else
        call m_Parallel_mpi_fft_box_cd_xyz(nfout,ipriparallel,printable,fft_box_size_CD_3D,kimg, &
    &                               fftbox_div_1,fftbox_div_2)
     endif
#endif
!     call m_Parallel_init_mpi_kngp_3D(nfout,ipriparallel,kgp)  ! -(m_Parallelization) ->ista_kngp,iend_kngp
     call m_Parallel_init_mpi_kngp_B_3D(nfout,ipriparallel,kgp)  ! -(m_Parallelization) ->ista_kngp_B,iend_kngp_B
     call m_pwBS_set_ngabc_B_3D
     call m_pwBS_set_ngabc_3D
!FNS     call m_pwBS_dealloc_ngabc
!!!  call m_CtrlP_flag_mpi_G_dot_R(nfout,nbmx) ! -> flag_mpi_g_dot_r
     call m_Parallel_init_mpi_nbmx(nfout,ipriparallel,printable,nbmx,kg1,ngnode_nbmx,flag_mpi_g_dot_r,flag_mpi_g_dot_r_k)
     call m_Parallel_init_mpi_gga(nfout,ipriparallel,printable,nfftp,nfftps)
     if(sw_rsb==ON) call m_Parallel_init_mpi_ffth(nfout,ipriparallel,printable,nfft)
     call m_Parallel_init_mpi_snl_3D(nfout,ipriparallel,printable,nspin)
     !call m_Parallel_init_mpi_atm_f(nfout,ipriparallel,printable,natm)
     call m_Parallel_init_mpi_atm_ke(nfout,ipriparallel,printable,natm)
     call m_Parallel_init_mpi_atm(nfout,ipriparallel,printable,natm)
     call m_Parallel_init_mpi_atm2(nfout,ipriparallel,printable,natm2)
     call m_Parallel_init_mpi_atm_B_3D(nfout,ipriparallel,printable,natm,sw_communicator_for_chg == ON)
     call m_Parallel_init_mpi_mix(nfout,ipriparallel,printable,kgpm,sw_communicator_for_chg == ON)
#ifdef __EDA__
     if (sw_eda == ON) then
! -----  ascat starts modifying  -----
     call m_XC_alloc_exc_on_a_grid()
     call m_IS_alloc_zfm3_EDA(natm)
     call m_IS_alloc_eewald_per_atom(natm)
     call m_PP_alloc_PP_per_atom_etc
! -----  ascat ceases modifying  -----
  endif
#endif

  endif

  if((icond == INITIAL .or. icond == CONTINUATION .or.  icond==COORDINATE_CONTINUATION) &
       & .or.((icond==FIXED_CHARGE.or.icond==FIXED_CHARGE_CONTINUATION).and.ekmode==OFF &
       & .and. fixed_charge_k_parallel == ALL_AT_ONCE)) then
     call m_Parallel_init_mpi_iba_3D(nfout,ipriparallel,printable,kv3,iba) ! -> np_g1k, mp_g1k
  end if
!!!!$!BRANCH_P ORG_Parallel
  if((icond==FIXED_CHARGE .or. icond==FIXED_CHARGE_CONTINUATION) .and. &
    & fixed_charge_k_parallel==ONE_BY_ONE) &
    & call m_Parallel_init_mpi_kv3_ek(nfout,ipriparallel,printable,kv3_ek,nspin)
!!!!$!BRANCH_P_END ORG_Parallel


  call m_IS_alloc_zfm3_3D()
  call m_CD_alloc_chgq()
  call m_XC_alloc_vxc_3D()
!  call m_pwBS_alloc_ngpt_igfp_gr_3D()
!  call m_pwBS_calc_length_of_G_3D()
!  call m_pwBS_G_trans_functions_3D()
!  call m_pwBS_setup_FFTmapfunctions_3D()
! === For epsmain by tkato 2013/11/14 ==========================================
! call m_Parallel_wf_onto_fft_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma,k_symmetry,GAMMA,kg,kg_gamma,kv3,1)
! call m_Parallel_fft_onto_wf_3D(nfout,fft_box_size_WF,igf,nbase,kg,kv3,nfft,1)
  if(ekmode == OFF) then
     call m_Parallel_wf_onto_fft_3D(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
                                    k_symmetry,GAMMA,kg,kg_gamma,kv3,1)
     call m_Parallel_fft_onto_wf_3D(nfout,fft_box_size_WF,igf,nbase,kg,kv3,nfft,1)
#ifdef MPI_FFTW
     if(sw_mpi_fftw==ON) then
       call m_Parallel_wf_onto_fft_mpifftw(nfout,fft_box_size_WF,igf,nbase,nbase_gamma, &
                                      k_symmetry,GAMMA,kg,kg_gamma,kv3,1)
       call m_Parallel_fft_onto_wf_mpifftw(nfout,fft_box_size_WF,igf,nbase,kg,kv3,nfft,1)
       call m_FFT_init_mpifftw()
     endif
#endif
  end if
! ==============================================================================
  call m_Parallel_fft_onto_chgq_3D(nfout,fft_box_size_WF,igf,kg,nfft)
#ifdef MPI_FFTW
  if(sw_mpi_fftw==ON) then
    call m_Parallel_fft_onto_chgq_mpifftw(nfout,fft_box_size_WF,igf,kg,nfft)
  endif
#endif
  call m_Parallel_chgq_onto_fftcd_3D(nfout,kgp,fft_box_size_CD_3D,igfp_l)
  call m_Parallel_fftcd_onto_chgq_3D(nfout,kgp_reduced,fft_box_size_CD_3D,igfp_l)
#ifndef FFT_USE_SSL2
#ifndef FFT_3D_DIVISION
  lsize = max(maxval(nel_fft_x(:)),maxval(nel_fft_y(:)),maxval(nel_fft_z(:)))
  if (kimg==1) lsize=lsize*2
  allocate(dfft_l(lsize*kimg,1))
  if (sw_fft_xzy > 0) then
     call m_FFT_Direct_3D(nfout,dfft_l,lsize,1)
     call m_FFT_Inverse_3D(nfout,dfft_l,lsize,1)
  else
     call m_FFT_Direct_XYZ_3D(nfout,dfft_l,lsize,1)
     call m_FFT_Inverse_XYZ_3D(nfout,dfft_l,lsize,1)
  end if
  deallocate(dfft_l)
#endif
#endif

#endif

end subroutine Preparation_for_mpi

subroutine Preparation_for_mpi_ek
! $Id: Preparation_for_mpi.F90 603 2020-04-07 03:25:30Z jkoga $
!                           @(#)Preparation_for_mpi.F90 1.10 03/02/19 00:49:14
  use m_Kpoints,              only : kv3
  use m_Files,                only : nfout
  use m_Control_Parameters,   only : ipriparallel, printable
  use m_PlaneWaveBasisSet,    only : iba
  use m_Parallelization,      only : m_Parallel_init_mpi_iba_3D

  call m_Parallel_init_mpi_iba_3D(nfout,ipriparallel,printable,kv3,iba) !  -> np_g1k, mp_g1k

end subroutine Preparation_for_mpi_ek

subroutine Preparation_for_mpi_PAW()
  use m_Files,                only : nfout
  use m_Ionic_System,         only : natm
  use m_PseudoPotential,      only : mmesh, flg_paw
  use m_Parallelization,      only : m_Parallel_init_mpi_paw_3D
  use m_Control_Parameters,   only : ipriparallel
  if (flg_paw) then
    call m_Parallel_init_mpi_paw_3D(nfout,ipriparallel,natm,mmesh)
  endif
end subroutine Preparation_for_mpi_PAW
