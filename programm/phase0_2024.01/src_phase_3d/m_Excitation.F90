module m_Excitation
!
  use m_Const_Parameters,   only :  ON, OFF, DP, CMPLDP, LOWER, FMAXVALLEN, zi, &
       &                            ELECTRON, DIRECT, BUCS, PAI2, PAI4, Hartree, YES, Hartree
  use m_Control_Parameters,  only : sw_use_add_proj, kimg, neg, ndim_spinor, nspin, &
       &                            sw_excitation, noncol, ipri
  use m_Files,  only : nfout, nf_excitation_spectra, m_Files_open_nf_xi_spectra

  use m_PseudoPotential,    only :  nlmt, ilmt, nlmt_add, ilmt_add, itau, iloc, &
       &                            radr, wos, mmesh, nmesh, lpsmax, nloc, ntau, &
       &                            psirpw, phirpw, ltp, mtp, taup, ltp_add, mtp_add, &
       &                            nlmta, nlmta_add, lmta, lmta_add
  use m_NonLocal_Potential, only : new_radr_and_wos

  use m_Ionic_System,       only : ntyp, ityp, natm, pos
  use m_Crystal_Structure,   only : nopr, rltv, univol, op, tau
  use m_Kpoints,       only : kv3, vkxyz, qwgt, kv3_fbz, num_star_of_k, star_of_k, &
       &                      iopr_k_fbz_to_ibz, trev_k_fbz_to_ibz

  use m_PlaneWaveBasisSet,   only : kg1, kgp, ttr, ngabc, ngpt_l, iba, nbase, igf
  use m_Parallelization,     only : npes, ista_kngp, iend_kngp, MPI_CommGroup, ierr, &
       &                            ista_k, iend_k, np_e, myrank_e, myrank_k, &
       &                            map_e, map_k, map_z, ista_e, iend_e, istep_e, &
       &                            mpi_k_world, mype, map_ek, mpi_e_world, &
       &                            nrank_e, nrank_k

  use m_Electronic_Structure,  only : eko_l, zaj_l, fsr_l, fsi_l, fsr_add_l, fsi_add_l, &
       &                              m_ES_WF_in_Rspace1, m_ES_WF_in_Rspace_kt, &
       &                              efermi, occup_l

  use m_FFT,                    only : nfft, fft_box_size_WF, m_FFT_WF, &
       &                               m_FFT_alloc_WF_work, m_FFT_dealloc_WF_work

  use m_ValenceBand_Spectrum,  only : dipole_dxyz_us, m_VBS_find_ptrans_index_ek

  use m_Corelevel_Spectrum,   only : sw_corelevel_spectrum, atom_to_probe, &
       &                             e_low_cls => e_low, e_high_cls => e_high, &
       &                             e_step_cls => e_step, &
       &                             qnum_n_to_probe, qnum_l_to_probe, &
       &                             psir_core_ae_wfns, &
       &                             dipole_dxyz_core2val, &
       &                             m_CLS_chk_sw_corelevel_spectrum, &
       &                             m_CLS_rd_n_main, &
       &                             m_CLS_find_orb_index_to_probe, &
       &                             m_CLS_find_ptrans_indx_core2val, &
       &                             num_core_states, ene_core_states, &
       &                             psig_core_states, &
       &                             fsr_core_states, fsi_core_states, &
       &                             eshift_corelevel_virtual
  use mpi

  implicit none
!  include 'mpif.h'

! -- General --
  character(len("excitation")), private,parameter :: tag_excitation = "excitation"
  character(len("sw_excitation")), private,parameter :: &
       &                              tag_sw_excitation = "sw_excitation"
!
!!!  integer :: sw_excitation = OFF

! --- G expansion --
  character(len("G_expansion")), private,parameter :: &
       &                    tag_G_expansion = "G_expansion"
  character(len("nmax_G")),private,parameter :: tag_nmax_G_given = "nmax_G"
  character(len("sw_spin_decomposition")), parameter ::  &
       &            tag_sw_spin_decomposition = "sw_spin_decomposition"
  character(len("sw_spin_decomp")), parameter ::  &
       &            tag_sw_spin_decomp = "sw_spin_decomp"

  integer :: nmax_G_given = 100
  integer :: sw_spin_decomposition = OFF
!
  integer :: nspin_m
  integer :: nmax_G, nmax_G_spin

! --- Momentum transfer --
  character(len("sw_longwavelimit")),private,parameter :: &
       &           tag_sw_longwavelimit = "sw_longwavelimit"
  integer :: sw_LongWaveLimit = ON
  real(kind=DP) :: qvec(3) = (/0.0d0, 0.0d0, 1.0d-3/)

! -- band gap correction ---
  character(len("band_gap_correction")), parameter :: &
       &                     tag_band_gap_correction = "band_gap_correction"
  character(len("scissor_operator")), parameter ::  tag_scissor  = "scissor_operator"
  character(len("sw_scissor_renormalization")), parameter         :: &
       &           tag_sw_scissor_renormalization   = "sw_scissor_renormalization"
!
  real(kind=DP) :: scissor = 0.0d0
  integer :: sw_scissor_renormalization = OFF

! -- transtion moment ( G=0, q-> 0 )
  real(kind=DP) :: delta_omega = 1.0D-14

! --- Brillouin zone integral
  character(len("BZ_integration")), parameter ::  &
       &                  tag_BZ_integration  = "BZ_integration"
  character(len("method")), parameter ::  tag_method = "method"
  character(len("gaussian")), parameter ::  tag_gaussian = "gaussian"
  character(len("lorentzian")), parameter ::  tag_lorentzian = "lorentzian"

  character(len("g")), parameter ::  tag_g = "g"
  character(len("l")), parameter ::  tag_l = "l"

  character(len("width")), parameter ::  tag_width = "width"
  character(len("extent")), parameter ::  tag_extent = "extent"

  integer, parameter :: GAUSSIAN_B = 1, LORENTZIAN_B = 2
  integer :: way_BZ_integral = GAUSSIAN_B

  real(kind=DP) :: width = 0.5d0 / Hartree  ! 0.0183745D0              ! 0.5 eV
!  real(kind=DP) :: gaussian_extent = 3.893d0
  real(kind=DP) :: gaussian_extent = 4.00d0

! -- Energy range --
  character(len("energy_range")), parameter ::  tag_erange = "energy_range"
  character(len("e_low")), parameter ::  tag_elow = "e_low"
  character(len("e_high")), parameter ::  tag_ehigh = "e_high"
  character(len("e_step")), parameter ::  tag_estep = "e_step"

  integer :: nstep = 500
  real(kind=DP) :: e_low = 0.0d0, e_high = 1.0d0, e_step = 0.002d0

  real(kind=DP), allocatable :: e(:)

! -- Spectrum --
  character(len("spectrum")), parameter ::  tag_spectrum = "spectrum"
  character(len("type")), parameter ::  tag_type = "type"
  character(len("optics")), parameter ::  tag_optics = "optics"
  character(len("eels")), parameter ::  tag_eels = "eels"

  integer, parameter :: OPTICS = 1, EELS = 2

  integer :: spectrum_type = OPTICS
  integer :: sw_calc_epsilon_indep_only = ON

! --- WAVEDERF
  character(len("wavederf")), parameter :: &
       &                     tag_wavederf = "wavederf"
  character(len("sw_write_wavederf_file")), parameter :: &
       &                     tag_sw_write_wavederf_file = "sw_write_wavederf_file"
  character(len("sw_wavederf_full_bz")),private,parameter :: &
       &         tag_sw_wavederf_full_bz = "sw_wavederf_full_bz"
  character(len("split_wavederf_file")),private,parameter :: &
       &        tag_split_wavederf_file  = "split_wavederf_file"
  character(len("num_wavederf_files_once")),private,parameter :: &
       &        tag_num_wavederf_files_once  = "num_wavederf_files_once"
  character(len("wavederf_save_memory_mode")),private,parameter ::  &
       &        tag_wavederf_save_memory_mode  = "wavederf_save_memory_mode"
  character(len("wavederf_sort_kpt")),private,parameter ::  &
       &        tag_wavederf_sort_kpt  = "wavederf_sort_kpt"
  integer :: sw_write_wavederf_file = OFF
  integer :: sw_wavederf_full_bz = OFF
  integer :: split_wavederf_file = OFF
  integer :: num_wavederf_files_once = 8
  integer :: wavederf_save_memory_mode = 0
  integer :: wavederf_sort_kpt = OFF
!
  character(len("wavederf_window_width")),private,parameter ::  &
       &        tag_wavederf_window_width = "wavederf_window_width"
  real(kind=DP) :: wavederf_window_width = -99.0d0

! -- Interaction --
  complex(kind=CMPLDP), allocatable :: Kernel_Coulomb(:)

! -- LR-TDDFT --
  character(len("tddft")), parameter ::  tag_tddft = "tddft"
  character(len("sw_tddft")), parameter ::  tag_sw_tddft = "sw_tddft"
  character(len("sw_nlf")), parameter ::  tag_sw_nlf = "sw_nlf"

  character(len("xc_kernel")), parameter ::  tag_xc_kernel = "xc_kernel"
  character(len("rpa")), parameter ::  tag_rpa = "rpa"
  character(len("lrc")), parameter ::  tag_lrc = "lrc"
  character(len("lrc_scf")), parameter ::  tag_lrc_scf = "lrc_scf"
  character(len("alpha_lrc")), parameter ::  tag_alpha_lrc = "alpha_lrc"
  character(len("bootstrap")), parameter ::  tag_bootstrap = "bootstrap"

  integer :: sw_tddft = OFF
  integer :: sw_nlf = OFF

  integer, parameter :: RPA = 1, LRC = 5, LRC_scf = 6, BOOTSTRAP = 10
  integer :: kernel_xc_type = RPA
  real(kind=DP) :: alpha_LRC = 0.2d0

! -- MBPT Bethe-Salpeter --
  character(len("bse")), parameter ::  tag_bse = "bse"
  character(len("sw_bse")), parameter ::  tag_sw_bse = "sw_bse"

  integer :: sw_bse = OFF

  integer :: num_occ_bands, num_unocc_bands, matsize_bse, ista_matbse, iend_matbse
  complex(kind=CMPLDP), allocatable :: epsinv_omega0(:,:)
  complex(kind=CMPLDP), allocatable :: MatH_bse(:,:)

! -- variables --
! --
! valence <-> valence
!
  integer :: nqitg_XI_vv_sum
  integer, allocatable :: iqitg_XI_vv(:,:,:,:,:,:)
  integer, allocatable :: nqitg_XI_vv(:)
  real(kind=DP), allocatable :: qitg_XI_vv(:,:)
!
  integer, allocatable :: isph_XI_vv(:,:,:,:)
  integer, allocatable :: il2p_XI_vv(:,:,:)
  real(kind=DP), allocatable :: dl2p_XI_vv(:,:,:,:)
!
  real(kind=DP), allocatable :: Mat_dipole_corr_vv(:,:,:,:)
!
  complex(kind=CMPLDP), allocatable :: trm_vv(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: RhoTilde_vv(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: SpectrFn_vv(:,:,:)
  complex(kind=CMPLDP), allocatable :: SpectrTensor_vv(:,:,:,:)

! --
! core <-> valence
!
  integer :: nqitg_XI_vc
  integer, allocatable :: iqitg_XI_vc(:,:,:)
  real(kind=DP), allocatable :: qitg_XI_vc(:,:)
!
  integer, allocatable :: isph_XI_vc(:,:,:)
  integer, allocatable :: il2p_XI_vc(:,:)
  real(kind=DP), allocatable :: dl2p_XI_vc(:,:,:)
!
  real(kind=DP), allocatable :: Mat_dipole_corr_vc(:,:,:)

  complex(kind=CMPLDP), allocatable :: trm_vc(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: RhoTilde_vc(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: SpectrFn_vc(:,:,:)

! --
! core <-> core
!
  integer :: nqitg_XI_cc
  integer, allocatable :: iqitg_XI_cc(:,:,:)
  real(kind=DP), allocatable :: qitg_XI_cc(:,:)
!
  integer, allocatable :: isph_XI_cc(:,:,:)
  integer, allocatable :: il2p_XI_cc(:,:)
  real(kind=DP), allocatable :: dl2p_XI_cc(:,:,:)
!
  complex(kind=CMPLDP), allocatable :: trm_cc(:,:,:,:)
  complex(kind=CMPLDP), allocatable :: RhoTilde_cc(:,:,:,:)

! symmetry
  integer, allocatable :: ngpt_XI(:,:)

! others
  integer :: nlmt_val, nlmt_core
  integer, allocatable :: ilmt_val(:)

  integer, parameter :: lcmax = 4
  integer :: ipriexcitation = 2
  integer :: imple_method = 1       ! 1 : like gpaw, 2 :like yambo
!
contains

! --------------------------------------------------------
!    Parameter setting
! --------------------------------------------------------

  subroutine m_XI_read_input
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop

    integer :: iret
    real(kind=DP) :: dret
    character(len=FMAXVALLEN) :: rstr

    iret = f_selectTop()
    if( f_selectBlock( tag_excitation ) == 0 ) then

       if ( f_getIntValue( tag_sw_excitation, iret ) == 0 ) then
          sw_excitation = iret
          write(nfout,*)
          write(nfout,*) "!** XI: sw_excitation is set to ", iret
       endif
       if ( sw_excitation == OFF ) return
       !
       call m_CLS_chk_sw_corelevel_spectrum
       if ( sw_corelevel_spectrum == ON ) call m_CLS_rd_n_main

       ! ---
       call set_lr_tddft
       call set_mbpt_bse
       call set_nmax_G_for_expansion
       call set_scissor_operator
       call set_BZ_integration_method
       call set_energy_range
       call set_spectrum_type

       if ( f_selectBlock( tag_wavederf ) == 0 ) then
          if ( f_getIntValue( tag_sw_write_wavederf_file, iret ) == 0 ) then
             sw_write_wavederf_file = iret
             write(nfout,*) '!** sw_write_wavederf_file is set to ', iret
          endif
          if ( sw_write_wavederf_file == ON ) then
             if( f_getIntValue(tag_sw_wavederf_full_bz, iret) == 0) then
                sw_wavederf_full_bz = iret
                write(nfout,*) '!** sw_wavederf_full_bz is set to ', &
                     &          sw_wavederf_full_bz
             endif
             if( f_getIntValue(tag_wavederf_save_memory_mode, iret) == 0) then
                wavederf_save_memory_mode = iret
                write(nfout,*) '!** wavederf_save_memory_mode is set to ', &
                     &          wavederf_save_memory_mode
             endif
             if ( f_getIntValue( tag_split_wavederf_file, iret ) == 0 ) then
                split_wavederf_file = iret
                write(nfout,*) '!** split_wavederf_file is set to ', &
                     &          split_wavederf_file
             endif
             if ( f_getIntValue( tag_num_wavederf_files_once, iret ) == 0 ) then
                num_wavederf_files_once = iret
                write(nfout,*) '!** num_wavederf_files_once is set to ', &
                     &          num_wavederf_files_once
             endif
             if( f_getIntValue(tag_wavederf_sort_kpt, iret) == 0) then
                wavederf_sort_kpt = iret
                write(nfout,*) '!** wavederf_sort_kpt is set to ', &
                     &          wavederf_sort_kpt
             endif
             if( f_getRealValue(tag_wavederf_window_width, dret, "hartree" ) == 0) then
                wavederf_window_width = dret
                write(nfout,*) '!** wavederf_window_width is set to ', &
                     &          wavederf_window_width
             endif
          endif
!          write(nfout,*) "!** sw_use_add_proj is turned on"
!          sw_use_add_proj = ON
      endif
      iret = f_selectParentBlock()
   end if

  contains

    subroutine set_nmax_G_for_expansion
      integer :: f_getIntValue

      if ( f_selectBlock( tag_G_expansion ) == 0 ) then
         if ( f_getIntValue( tag_nmax_G_given, iret ) == 0 ) then
            if ( iret > 0 .and. iret < kg1 ) nmax_G_given = iret
         endif
         if ( f_getIntValue( tag_sw_spin_decomposition, iret ) == 0 ) then
            sw_spin_decomposition = iret
         endif
         if ( f_getIntValue( tag_sw_spin_decomp, iret ) == 0 ) then
            sw_spin_decomposition = iret
         endif

         iret = f_selectParentBlock()
      endif
      write(nfout,*) "!** XI: sw_spin_decomp is set to ", sw_spin_decomposition
    end subroutine set_nmax_G_for_expansion

    subroutine set_scissor_operator
      integer :: f_getRealValue

      if ( f_selectBlock( tag_band_gap_correction ) == 0 ) then
         if ( f_getRealValue( tag_scissor, dret, 'hartree' ) == 0 ) then
            scissor = dret
         endif
         if ( f_getIntValue( tag_sw_scissor_renormalization, iret ) == 0 ) then
            sw_scissor_renormalization = iret
         endif
         iret = f_selectParentBlock()
      endif
      write(nfout,*) "!** XI: scissor operator is set to ", scissor
      write(nfout,*) "!** XI: sw_scissor_renormalization is set to ", sw_scissor_renormalization
    end subroutine set_scissor_operator

    subroutine set_BZ_integration_method
      integer :: f_getRealValue, f_getStringValue
      logical :: tf

      if ( f_selectBlock(tag_BZ_integration) == 0 ) then
         if ( f_getStringValue( tag_method, rstr, LOWER ) == 0 ) then
            call strncmp0( tag_gaussian, trim(rstr), tf)
            if (tf) way_BZ_integral = Gaussian_B
            call strncmp0( tag_g, trim(rstr), tf)
            if (tf) way_BZ_integral = Gaussian_B

            call strncmp0( tag_lorentzian, trim(rstr), tf)
            if (tf) way_BZ_integral = Lorentzian_B
            call strncmp0( tag_l, trim(rstr), tf)
            if (tf) way_BZ_integral = Lorentzian_B
         endif

         if( f_getRealValue( tag_width, dret, 'hartree' ) == 0) then
            if ( dret > 0.0 ) width = dret
         endif
         if( f_getRealValue( tag_extent, dret, '' ) == 0) then
            if ( dret > 2.0 ) gaussian_extent = dret
         endif
         iret = f_selectParentBlock()
      end if

      write(nfout,*) "!** XI: BZ integration method is set to ", way_BZ_integral
      write(nfout,*) "!** XI: width is set to ", width
      if ( way_BZ_integral == Gaussian_B ) then
         write(nfout,*) "!** XI: gaussian_extent is set to ", gaussian_extent
      endif

    end subroutine set_BZ_integration_method

    subroutine set_energy_range
      integer :: f_getRealValue

       if ( sw_corelevel_spectrum == ON ) then
          e_low  = e_low_CLS
          e_high = e_high_CLS
          e_step = e_step_CLS
       endif

      if ( f_selectBlock( tag_erange ) == 0 ) then
         if( f_getRealValue( tag_elow, dret, 'hartree' ) == 0) then
            if ( dret > 0.0 ) e_low = dret
         endif
         if( f_getRealValue( tag_ehigh, dret,'hartree' ) == 0) then
            if ( dret > 0.0 .and. dret > e_low ) e_high = dret
         endif
         if( f_getRealValue( tag_estep, dret,'hartree' ) == 0) then
            if ( dret > 0.0 ) e_step = dret
         endif
         iret = f_selectParentBlock()
       end if

       write(nfout,*) "!** XI: e_low  is ", e_low
       write(nfout,*) "!** XI: e_high is ", e_high
       write(nfout,*) "!** XI: e_step is ", e_step

       nstep = nint( ( e_high -e_low ) /e_step )

       write(nfout,*) "!** XI: number of energy steps is ", nstep

    end subroutine set_energy_range

    subroutine set_lr_tddft
      integer :: f_getIntValue, f_getStringValue, f_getRealValue
      logical :: tf

      if ( f_selectBlock( tag_tddft ) == 0 ) then
         if ( f_getIntValue( tag_sw_tddft, iret ) == 0 ) then
            sw_tddft = iret
         endif
         if ( f_getIntValue( tag_sw_nlf, iret ) == 0 ) then
            sw_nlf = iret
         endif
         if ( f_getStringValue( tag_xc_kernel, rstr, LOWER ) == 0 ) then
            call strncmp0( tag_rpa, trim(rstr), tf)
            if (tf) kernel_xc_type = RPA

            if (len_trim(rstr) == len_trim(tag_lrc_scf) ) then
               call strncmp0( tag_lrc_scf, trim(rstr), tf)
               if (tf) kernel_xc_type = LRC_scf
            endif

            if (len_trim(rstr) == len_trim(tag_lrc) ) then
               call strncmp0( tag_lrc, trim(rstr), tf)
               if (tf) kernel_xc_type = LRC
            endif

            call strncmp0( tag_bootstrap, trim(rstr), tf)
            if (tf) kernel_xc_type = BOOTSTRAP
         endif

100      continue

         if( f_getRealValue( tag_alpha_lrc, dret,'' ) == 0) then
            alpha_lrc = dret
         endif

         write(nfout,*) "!** XI: sw_tddft is ", sw_tddft
         if ( sw_tddft == ON ) then
            write(nfout,*) "!** XI: sw_nlf   is ", sw_nlf
            write(nfout,*) "!** XI: kernel_xc_type is ", kernel_xc_type
            if ( kernel_xc_type == LRC ) then
               write(nfout,*) "!** XI: alpha_LRC is ", alpha_LRC
            endif
         endif

         if ( sw_tddft == ON ) sw_calc_epsilon_indep_only = OFF

         iret = f_selectParentBlock()
      endif
    end subroutine set_lr_tddft

    subroutine set_mbpt_bse
      integer :: f_getIntValue

      if ( f_selectBlock( tag_bse ) == 0 ) then
         if ( f_getIntValue( tag_sw_bse, iret ) == 0 ) then
            sw_bse = iret
         endif

         write(nfout,*) "!** XI: sw_bse is ", sw_bse
         if ( sw_bse == ON ) sw_calc_epsilon_indep_only = OFF

         iret = f_selectParentBlock()
      endif
    end subroutine set_mbpt_bse

    subroutine set_spectrum_type
      integer :: f_getIntValue, f_getStringValue
      logical :: tf

      if ( f_selectBlock( tag_spectrum ) == 0 ) then
         if ( f_getStringValue( tag_type, rstr, LOWER ) == 0 ) then
            call strncmp0( tag_optics, trim(rstr), tf)
            if (tf) spectrum_type = OPTICS

            call strncmp0( tag_eels, trim(rstr), tf)
            if (tf) spectrum_type = EELS
         endif
         iret = f_selectParentBlock()
      end if

      write(nfout,*) "!** XI: Spectrum type is set to ", spectrum_type

    end subroutine set_spectrum_type

  end subroutine m_XI_read_input

!  subroutine m_XI_set_value_nmax_G( nmax_G_given, nmax_G )
  subroutine m_XI_set_value_nmax_G
!    integer, intent(in) :: nmax_G_given
!    integer, intent(out) :: nmax_G

    integer :: i, count, ispin
    real(kind=DP) :: g1, g1_prev, dg

    if ( sw_spin_decomposition == ON ) then
       nspin_m = nspin
    else
       nspin_m = 1
    endif

    if ( sw_calc_epsilon_indep_only == ON ) then
       nmax_G = 1;  nmax_G_spin = nspin_m;  return
    endif

    count = 0
    do i = 1, kgp
       g1 = dsqrt( ttr(1)*ngabc(i,1)*ngabc(i,1) &
            &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
            &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
            &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
            &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
            &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
       dg = abs( g1 -g1_prev )
       if ( dg < 1.0D-5 ) then
          count = count +1
       else
          if ( count > nmax_G_given ) exit
          count = count +1
       end if
       g1_prev = g1
    enddo
    nmax_G = count
    nmax_G_spin = nmax_G *nspin_m

    write(nfout,*)
    write(nfout,*) "!** XI: nmax_G is set to ", nmax_G
    write(nfout,*)

  end subroutine m_XI_set_value_nmax_G

! --------------------------------------------------------
!    Preparation of arrays for treating wfns
! --------------------------------------------------------

  subroutine m_XI_set_ilmt_val
    integer :: it

    allocate( ilmt_val(ntyp) )

    if ( sw_use_add_proj == ON ) then
       nlmt_val = nlmt +nlmt_add
       Do it=1, ntyp
          ilmt_val(it) = ilmt(it) +ilmt_add(it)
       End do
    else
       nlmt_val = nlmt
       Do it=1, ntyp
          ilmt_val(it) = ilmt(it)
       End do
    endif

  end subroutine m_XI_set_ilmt_val

  subroutine m_XI_set_ilmt_core
    nlmt_core = 2 *qnum_l_to_probe +1
  end subroutine m_XI_set_ilmt_core

  subroutine m_XI_set_ngpt_XI
    integer, allocatable :: ngpt_wk(:,:), ngpt_wk2(:,:)

    integer :: i

    allocate( ngpt_XI( nmax_G, nopr) ); ngpt_XI = 0

    if ( npes > 1 ) then
       allocate( ngpt_wk( kgp,nopr) ); ngpt_wk  = 0
       allocate( ngpt_wk2(kgp,nopr) ); ngpt_wk2 = 0
       ngpt_wk(ista_kngp:iend_kngp,:) = ngpt_l(ista_kngp:iend_kngp,:)

       call mpi_allreduce( ngpt_wk, ngpt_wk2, kgp*nopr, mpi_integer, mpi_sum, &
            &              MPI_CommGroup, ierr )
       ngpt_XI( 1:nmax_G, 1:nopr ) = ngpt_wk2( 1:nmax_G, 1:nopr )
       deallocate( ngpt_wk ); deallocate( ngpt_wk2 )
    else
       ngpt_XI( 1:nmax_G, 1:nopr ) = ngpt_l( 1:nmax_G, 1:nopr )
    endif

  end subroutine m_XI_set_ngpt_XI

! --------------------------------------------------------
!    Preparation of arrays for transtion moments (val-val)
! --------------------------------------------------------

  subroutine m_XI_set_qitg_XI_vv
    allocate( nqitg_XI_vv(ntyp) ); nqitg_XI_vv = 0
    call set_size_qitg_XI_vv

    if ( ipriexcitation >=2 ) then
       write(nfout,*) '!** XI: nqitg_XI_vv_sum  = ',nqitg_XI_vv_sum
    endif

    allocate( iqitg_XI_vv( nloc, ntau, nloc, ntau, nloc**2, ntyp ) )
    allocate( qitg_XI_vv( nmax_G, nqitg_XI_vv_sum ) )
    iqitg_XI_vv = 0;  qitg_XI_vv = 0.0d0

    call set_array_qitg_XI_vv

  contains

    subroutine set_size_qitg_XI_vv
      integer :: it, il1, il2, tau1, tau2, tmin
      integer :: l3s, l3l, il3, mm_sum, mm

      mm_sum = 0

      Do it=1, ntyp
         mm = 0

         Loop_L1 : do il1 = 1, lpsmax(it)
            Loop_tau1 : do tau1 = 1, itau(il1,it)

               Loop_L2 : do il2 = il1, lpsmax(it)
                  tmin = 1
                  if(il1 == il2) tmin = tau1

                  Loop_tau2 : do tau2 = tmin, itau(il2,it)
                     l3s = abs(il1-il2); l3l = abs(il1+il2) - 2

                     Loop_L3 : do il3 = l3s+1, l3l+1, 2
                        if(il3-1 > lcmax ) cycle
                        mm = mm +1
                     End do Loop_L3

                  End do Loop_tau2
               End do Loop_L2

            End do Loop_tau1
         End do Loop_L1

         nqitg_XI_vv(it) = mm
         mm_sum = mm_sum +mm
      end Do

      nqitg_XI_vv_sum = mm_sum

    end subroutine set_size_qitg_XI_vv

    subroutine set_array_qitg_XI_vv
      integer :: it, il1, il2, tau1, tau2, tmin
      integer :: l3s, l3l, il3, mm_sum, mm
      integer :: i, ir

      real(kind=DP), allocatable :: qrs(:), gr_XI(:), wkx(:), wky(:)

      allocate( qrs( mmesh ) ); qrs = 0.d0
      allocate( wkx( mmesh ) ); wkx = 0.0d0
      allocate( wky( mmesh ) ); wky = 0.0d0
      allocate( gr_XI(nmax_G) ); gr_XI = 0.0d0

      allocate( radr(mmesh) );  allocate( wos(mmesh) )

      do i = 1, nmax_G
         gr_XI(i) =      dsqrt(ttr(1)*ngabc(i,1)*ngabc(i,1) &
              &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
              &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
              &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
              &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
              &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
      enddo

      mm = 0
      Do it=1, ntyp
         call new_radr_and_wos(ista_k,it)

         Loop_L1 : do il1 = 1, lpsmax(it)
            Loop_tau1 : do tau1 = 1, itau(il1,it)

               Loop_L2 : do il2 = il1, lpsmax(it)
                  tmin = 1
                  if (il1 == il2) tmin = tau1

                  Loop_tau2 : do tau2 = tmin, itau(il2,it)
                     l3s = abs(il1-il2);  l3l = abs(il1+il2) - 2

                     Loop_L3 : do il3 = l3s+1, l3l+1, 2
                        if(il3-1 > lcmax ) cycle
                        mm = mm +1
                        iqitg_XI_vv( il1, tau1, il2, tau2, il3, it ) = mm
                        iqitg_XI_vv( il2, tau2, il1, tau1, il3, it ) = mm  !! ??

                        qrs = 0.0d0
                        Do ir=1, nmesh(it)
                           qrs( ir ) &
                                &   = ( psirpw(ir,il1,tau1,it) *psirpw(ir,il2,tau2,it) &
                                &     - phirpw(ir,il1,tau1,it) *phirpw(ir,il2,tau2,it) )
                        End do

                        call qitgft( 1, nmax_G, mmesh, nqitg_XI_vv_sum, nmesh(it), &
                             &       ipriexcitation, nfout, il3-1, radr, wos, &
                             &       qrs, gr_XI, &
                             &       1.d0, mm, qitg_XI_vv, wkx, wky )

                     End do Loop_L3

                  End do Loop_tau2
               End do Loop_L2

            End do Loop_tau1
         End do Loop_L1
      End Do

      deallocate( radr ); deallocate( wos );  deallocate( wkx ); deallocate( wky )
      deallocate( qrs ); deallocate( gr_XI )

    end subroutine set_array_qitg_XI_vv

  end subroutine m_XI_set_qitg_XI_vv

  subroutine m_XI_set_dl2p_XI_vv
    integer :: it, lmt1, lmt2, il1, im1, il2, im2, itmp
    integer :: nspher1, nspher2, n

    integer, allocatable :: isph2(:,:,:), mmt2(:,:)
    real(kind=DP), allocatable :: cr2(:,:,:)

    allocate( isph_XI_vv( nlmt_val, nlmt_val, 6, ntyp )); isph_XI_vv = 0
    allocate( il2p_XI_vv( nlmt_val, nlmt_val, ntyp ) );   il2p_XI_vv = 0
    allocate( dl2p_XI_vv( nlmt_val, nlmt_val, 6, ntyp));  dl2p_XI_vv = 0

    allocate( cr2(  16,16,6) ); cr2 = 0.0d0
    allocate( isph2(16,16,6) ); isph2 = 0.0d0
    allocate( mmt2( 16,16)   ); mmt2 = 0

    call sphset2( nfout, ipri, lcmax, cr2, isph2, mmt2 )

    Do it=1, ntyp
       Do lmt1=1, ilmt_val(it)
          if ( lmt1 <= ilmt(it) ) then
             il1 = ltp(lmt1,it);  im1 = mtp(lmt1,it)
          else
             itmp = lmt1 -ilmt(it)
             il1 = ltp_add(itmp,it);  im1 = mtp_add(itmp,it)
          endif
          nspher1 = ( il1 -1 )**2 +im1

          Do lmt2=1, ilmt_val(it)
             if ( lmt2 <= ilmt(it) ) then
                il2 = ltp(lmt2,it);  im2 = mtp(lmt2,it)
             else
                itmp = lmt2 -ilmt(it)
                il2 = ltp_add(itmp,it);  im2 = mtp_add(itmp,it)
             endif
             nspher2 = ( il2 -1 )**2 +im2

             il2p_XI_vv( lmt1, lmt2, it ) = mmt2( nspher1, nspher2 )
             Do n=1, il2p_XI_vv(lmt1,lmt2,it)
                isph_XI_vv( lmt1, lmt2, n, it ) = isph2( nspher1, nspher2, n )
                dl2p_XI_vv( lmt1, lmt2, n, it ) = cr2  ( nspher1, nspher2, n )
             End do
          End Do
       End Do
    End Do

  end subroutine m_XI_set_dl2p_XI_vv

  subroutine m_XI_set_mat_dipole_corr_vv
    integer :: it, lmt1, lmt2, il1, il2, im1, im2, tau1, tau2
    integer :: nspher1, nspher2, itmp, ifact, index

    allocate( Mat_Dipole_corr_vv( ntyp, nlmt_val, nlmt_val, 3 ) )
    Mat_Dipole_corr_vv = 0.0d0

    Do it=1, ntyp
       Do lmt1=1, ilmt_val(it)
          if ( lmt1 <= ilmt(it) ) then
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it); tau1 = taup(lmt1,it)
          else
             itmp = lmt1 -ilmt(it)
             il1 = ltp_add(itmp,it); im1 = mtp_add(itmp,it);  tau1 = 1
          endif
          nspher1 = ( il1 -1 )**2 +im1

          Do lmt2=1, ilmt_val(it)
             if ( lmt2 <= ilmt(it) ) then
                il2 = ltp(lmt2,it); im2 = mtp(lmt2,it); tau2 = taup(lmt2,it)
             else
                itmp = lmt2 -ilmt(it)
                il2 = ltp_add(itmp,it); im2 = mtp_add(itmp,it);  tau2 = 1
             endif
             nspher2 = ( il2 -1 )**2 +im2

             call m_VBS_find_ptrans_index_ek( it, nspher1, nspher2, tau1, tau2, &
                  &                           index, ifact )

             if ( index /= 0 ) then
                Mat_Dipole_Corr_vv( it, lmt1, lmt2, 1:3 ) &
                     &      = dipole_dxyz_us( it, index, 1:3 ) *ifact
             endif
          End do
       End Do
    End Do

  end subroutine m_XI_set_mat_dipole_corr_vv

! --------------------------------------------------------
!    Preparation of arrays for transtion moments (val-core)
! --------------------------------------------------------

  subroutine m_XI_set_qitg_XI_vc
    nqitg_XI_vc = 0
    call set_size_qitg_XI_vc

    if ( ipriexcitation >= 2 ) then
       write(nfout,*) '** nqitg_XI_vc_sum  = ',nqitg_XI_vc
    endif

    allocate( iqitg_XI_vc( nloc, ntau, nloc**2 ) )
    allocate( qitg_XI_vc( nmax_G, nqitg_XI_vc ) )

    iqitg_XI_vc = 0;  qitg_XI_vc = 0.0d0

    call set_array_qitg_XI_vc

  contains

    subroutine set_size_qitg_XI_vc
      integer :: mm_sum, mm
      integer :: it, il1, il2, tau1, tau2, tmin
      integer :: l3s, l3l, il3, atomtype_to_probe

      mm_sum = 0

      atomtype_to_probe = ityp( atom_to_probe )

      Do it = atomtype_to_probe, atomtype_to_probe
         mm = 0

         Loop_L1 : do il1 = 1, lpsmax(it)
            Loop_tau1 : do tau1 = 1, itau(il1,it)

               Loop_L2 : do il2 = qnum_l_to_probe+1, qnum_l_to_probe+1
                  tmin = 1
                  if(il1 == il2) tmin = tau1

                  Loop_tau2 : do tau2 = 1, 1
                     l3s = abs(il1-il2); l3l = abs(il1+il2) - 2

                     Loop_L3 : do il3 = l3s+1, l3l+1, 2
                        if(il3-1 > lcmax ) cycle
                        mm = mm +1
                     End do Loop_L3

                  End do Loop_tau2
               End do Loop_L2

            End do Loop_tau1
         End do Loop_L1

         nqitg_XI_vc = mm
         mm_sum = mm_sum +mm
      end Do

!      nqitg_XI_vc_sum = mm_sum

    end subroutine set_size_qitg_XI_vc

    subroutine set_array_qitg_XI_vc
      integer :: mm_sum, mm
      integer :: it, il1, il2, tau1, tau2, tmin
      integer :: l3s, l3l, il3
      integer :: i, ir, orb_index, atomtype_to_probe

      real(kind=DP), allocatable :: qrs(:), gr_XI(:), wkx(:), wky(:)

      allocate( qrs( mmesh ) ); qrs = 0.d0
      allocate( wkx( mmesh ) ); wkx = 0.0d0
      allocate( wky( mmesh ) ); wky = 0.0d0
      allocate( gr_XI(nmax_G) ); gr_XI = 0.0d0

      allocate( radr(mmesh) );  allocate( wos(mmesh) )

      do i = 1, nmax_G
         gr_XI(i) =      dsqrt(ttr(1)*ngabc(i,1)*ngabc(i,1) &
              &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
              &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
              &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
              &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
              &             + ttr(6)*ngabc(i,3)*ngabc(i,1))
      enddo

!      write(*,*) "psircore", allocated( psir_core_ae_wfns )

      orb_index = m_CLS_find_orb_index_to_probe()

      mm = 0

      atomtype_to_probe = ityp( atom_to_probe )

      Do it = atomtype_to_probe, atomtype_to_probe
         call new_radr_and_wos(ista_k,it)

         Loop_L1 : do il1 = 1, lpsmax(it)
            Loop_tau1 : do tau1 = 1, itau(il1,it)

               Loop_L2 : do il2 = qnum_l_to_probe+1, qnum_l_to_probe+1
                  tmin = 1
                  if (il1 == il2) tmin = tau1

                  Loop_tau2 : do tau2 = 1, 1
                     l3s = abs(il1-il2);  l3l = abs(il1+il2) - 2

                     Loop_L3 : do il3 = l3s+1, l3l+1, 2
                        if(il3-1 > lcmax ) cycle
                        mm = mm +1
                        iqitg_XI_vc( il1, tau1, il3 ) = mm

                        qrs = 0.0d0
                        Do ir=1, nmesh(it)
                           qrs( ir ) &
                                &   =  (  psirpw(ir,il1,tau1,it) &
                                &        -phirpw(ir,il1,tau1,it) ) &
                                &          *psir_core_ae_wfns(ir,orb_index)
                        End do

                        call qitgft( 1, nmax_G, mmesh, nqitg_XI_vc, nmesh(it), &
                             &       ipriexcitation, nfout, il3-1, radr, wos, &
                             &       qrs, gr_XI, &
                             &       1.d0, mm, qitg_XI_vc, wkx, wky )

                     End do Loop_L3

                  End do Loop_tau2
               End do Loop_L2

            End do Loop_tau1
         End do Loop_L1
      End Do

      deallocate( radr ); deallocate( wos );  deallocate( wkx ); deallocate( wky )
      deallocate( qrs ); deallocate( gr_XI )

    end subroutine set_array_qitg_XI_vc

  end subroutine m_XI_set_qitg_XI_vc

  subroutine m_XI_set_dl2p_XI_vc
    integer :: it, lmt1, lmt2, il1, im1, il2, im2, itmp
    integer :: nspher1, nspher2, n, atomtype_to_probe

    integer, allocatable :: isph2(:,:,:), mmt2(:,:)
    real(kind=DP), allocatable :: cr2(:,:,:)

    allocate( isph_XI_vc( nlmt_val, nlmt_core, 6 )); isph_XI_vc = 0
    allocate( il2p_XI_vc( nlmt_val, nlmt_core ) );   il2p_XI_vc = 0
    allocate( dl2p_XI_vc( nlmt_val, nlmt_core, 6 )); dl2p_XI_vc = 0

    allocate( cr2(  16,16,6) ); cr2 = 0.0d0
    allocate( isph2(16,16,6) ); isph2 = 0.0d0
    allocate( mmt2( 16,16)   ); mmt2 = 0

    call sphset2( nfout, ipri, lcmax, cr2, isph2, mmt2 )

    atomtype_to_probe = ityp( atom_to_probe )

    Do it = atomtype_to_probe, atomtype_to_probe
       Do lmt1=1, ilmt_val(it)
          if ( lmt1 <= ilmt(it) ) then
             il1 = ltp(lmt1,it);  im1 = mtp(lmt1,it)
          else
             itmp = lmt1 -ilmt(it)
             il1 = ltp_add(itmp,it);  im1 = mtp_add(itmp,it)
          endif
          nspher1 = ( il1 -1 )**2 +im1

          Do lmt2=1, nlmt_core
             il2 = qnum_l_to_probe +1;  im2 = lmt2
             nspher2 = ( il2 -1 )**2 +im2

             il2p_XI_vc( lmt1, lmt2 ) = mmt2( nspher1, nspher2 )
             Do n=1, il2p_XI_vc(lmt1,lmt2)
                isph_XI_vc( lmt1, lmt2, n ) = isph2( nspher1, nspher2, n )
                dl2p_XI_vc( lmt1, lmt2, n ) = cr2  ( nspher1, nspher2, n )
             End do
          End Do
       End Do
    End Do

  end subroutine m_XI_set_dl2p_XI_vc

  subroutine m_XI_set_mat_dipole_corr_vc
    integer :: it, lmt1, lmt2, il1, il2, im1, im2, tau1, tau2
    integer :: nspher1, nspher2, itmp, ifact, index, atomtype_to_probe

    allocate( Mat_Dipole_corr_vc( nlmt_val, nlmt_core, 3 ) )
    Mat_Dipole_corr_vc = 0.0d0

    atomtype_to_probe = ityp( atom_to_probe )

    Do it = atomtype_to_probe, atomtype_to_probe
       Do lmt1=1, ilmt_val(it)
          if ( lmt1 <= ilmt(it) ) then
             il1 = ltp(lmt1,it); im1 = mtp(lmt1,it); tau1 = taup(lmt1,it)
          else
             itmp = lmt1 -ilmt(it)
             il1 = ltp_add(itmp,it); im1 = mtp_add(itmp,it);  tau1 = 1
          endif
          nspher1 = ( il1 -1 )**2 +im1

          Do lmt2=1, nlmt_core
             il2 = qnum_l_to_probe +1;  im2 = lmt2;  tau2 = 1
             nspher2 = ( il2 -1 )**2 +im2

             call m_CLS_find_ptrans_indx_core2val( qnum_n_to_probe, nspher2, nspher1, &
                  &                                tau1, index )
             ifact = -1

             if ( index /= 0 ) then
                Mat_Dipole_Corr_vc( lmt1, lmt2, 1:3 ) &
                     &      = dipole_dxyz_core2val( index, 1:3 ) *dble(ifact)
             endif
          End do
       End Do
    End Do

  end subroutine m_XI_set_mat_dipole_corr_vc

! --------------------------------------------------------
!    Calculation of transtion moment/matrix (val-val)
! --------------------------------------------------------

  subroutine m_XI_wd_WAVEDERF
    use m_Kpoints, only : m_Kp_set_star_of_k
    use m_CS_Magnetic,   only : m_CS_set_inverse_operation, invop, magmom_dir_inversion_opr_flag

    integer :: ib_min, ib_max, ib_num

    if ( allocated( trm_vv ) ) deallocate( trm_vv )
    allocate( trm_vv( np_e, neg, ista_k:iend_k, 3 ) );   trm_vv = 0.0d0

    call set_ib_min_max

    call calc_soft_part
    call calc_hard_part

    if ( sw_wavederf_full_bz == ON ) call m_Kp_set_star_of_k
    if ( noncol ) then
       if ( split_wavederf_file == YES ) then
          call print_wavederf_noncl_split
       else
          select case (wavederf_save_memory_mode)
          case (0)
             if ( sw_wavederf_full_BZ == ON .and. wavederf_sort_kpt == ON  ) then
                call print_wavederf_noncl_fbz_sort
             else
                call print_wavederf_noncl
             endif
          case (1)
             call print_wavederf_noncl2
          end select
       endif
    else
       if ( split_wavederf_file == YES ) then
          call print_wavederf_col_split
       else
          select case (wavederf_save_memory_mode)
          case (0)
             if ( sw_wavederf_full_BZ == ON .and. wavederf_sort_kpt == ON  ) then
                call print_wavederf_col_fbz_sort
             else
                call print_wavederf_col
             endif
          case (1)
             call print_wavederf_col2
          end select
       endif
    endif
    deallocate( trm_vv )

  contains

    subroutine set_ib_min_max
      integer :: ik, ib, ierr
      real(kind=DP) :: emin, emax, ene

      ib_min = 1;   ib_max = neg
      ib_num = ib_max -ib_min +1

      if ( wavederf_window_width < 0.0 ) return

      ib_min = neg;   ib_max = 1
      emin = efermi -wavederf_window_width /2.0d0
      emax = efermi +wavederf_window_width /2.0d0

      Do ik=ista_k, iend_k
         Do ib=ista_e, iend_e, istep_e
            ene = eko_l( map_z(ib),ik )
            if ( ene >= emin .and. ene <=emax ) then
               ib_min = min( ib, ib_min )
               ib_max = max( ib, ib_max )
            endif
         End Do
      End Do
      call mpi_allreduce( mpi_in_place, ib_min, 1, mpi_integer, &
           &              mpi_min, MPI_CommGroup, ierr )
      call mpi_allreduce( mpi_in_place, ib_max, 1, mpi_integer, &
           &              mpi_max, MPI_CommGroup, ierr )

      ib_num = ib_max -ib_min +1

    end subroutine set_ib_min_max

    subroutine print_wavederf_col_split
      integer :: ik, is, ib1, ib2
      integer :: lun
      integer :: trev, iopr, i, j, nr, ii, jk, num
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk2(:,:,:)
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      character*4 char1
      character*72 file0, file1

      integer istatus(mpi_status_size)

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) );   zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif

      lun = 10000 +myrank_k
      write(char1,'(I4.4)') myrank_k

      nr = nrank_k /num_wavederf_files_once +1
      num = 0

      iloop: Do i=1, nr
         Do j=1, num_wavederf_files_once
            num = num +1
            if ( num > npes )    goto 200
            if ( num /= myrank_k +1 ) cycle
            if ( myrank_e /= 0 )   cycle

            Do is=1, nspin
               if ( nspin == 1 ) then
                  file0 = "WAVEDERF_phase"
               else
                  if ( is == 1 ) then
                     file0 = "WAVEDERF_phase.up"
                  else
                     file0 = "WAVEDERF_phase.dn"
                  endif
               endif

               file1 = trim(adjustl(file0)) // '.' // char1
               if ( myrank_e == 0 ) then
                  open( lun, file=file1, status="unknown", form="formatted" )
               endif

               if ( mype == 0 .and. is == 1 ) then
                  if ( sw_wavederf_full_bz == ON ) then
                     write(lun,*) nspin, kv3_fbz/nspin, ib_num
                  else
                     write(lun,*) nspin, kv3/nspin, ib_num
                  endif
               endif

               ikloop1: Do ik=ista_k+is-1, iend_k, nspin

                  if ( sw_wavederf_full_bz == ON ) then
                     weight = kv3 *qwgt(ik) /dble(ndim_spinor)

                     Do ii=1, num_star_of_k(ik)
                        jk = star_of_k(ik,ii)
                        iopr = iopr_k_fbz_to_ibz(jk)
                        trev = trev_k_fbz_to_ibz(jk)

                        Do ib1=ib_min, ib_max
                           Do ib2=ib_min, ib_max
                              occ1 = occ_wk(ib1,ik) / weight
                              occ2 = occ_wk(ib2,ik) / weight
                              zvec(1) = zwk(ib1,ib2,1,ik)
                              zvec(2) = zwk(ib1,ib2,2,ik)
                              zvec(3) = zwk(ib1,ib2,3,ik)

                              zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                              if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
#if 1
                              if ( allocated(magmom_dir_inversion_opr_flag) ) then
                                 if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                                    zvec_rot = -conjg(zvec_rot)
                                 endif
                              endif
#endif
                              write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                                   &    ib1, ( eko_wk(ib1,ik)-efermi )*Hartree, occ1, &
                                   &    ib2, ( eko_wk(ib2,ik)-efermi )*Hartree, occ2, &
                                   &    real(zvec_rot(1)), aimag(zvec_rot(1)), &
                                   &    real(zvec_rot(2)), aimag(zvec_rot(2)), &
                                   &    real(zvec_rot(3)), aimag(zvec_rot(3))
                           End Do
                        End Do
                     End Do
                  else
                     weight = kv3 *qwgt(ik) /dble(ndim_spinor)
                     Do ib1=ib_min, ib_max
                        Do ib2=ib_min, ib_max
                           occ1 = occ_wk(ib1,ik) / weight
                           occ2 = occ_wk(ib2,ik) / weight
                           zvec(1) = zwk(ib1,ib2,1,ik)
                           zvec(2) = zwk(ib1,ib2,2,ik)
                           zvec(3) = zwk(ib1,ib2,3,ik)
                           write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                                &    ib1, ( eko_wk(ib1,ik)-efermi )*Hartree, occ1, &
                                &    ib2, ( eko_wk(ib2,ik)-efermi )*Hartree, occ2, &
                                &    real(zvec(1)), aimag(zvec(1)), &
                                &    real(zvec(2)), aimag(zvec(2)), &
                                &    real(zvec(3)), aimag(zvec(3))
                        End Do
                     End Do
                  end if
               End Do ikloop1
               if ( myrank_e == 0 ) close(lun)
            End Do
         End Do
200      continue
         call mpi_barrier( MPI_CommGroup, ierr )
      End Do iloop

      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )
    end subroutine print_wavederf_col_split

    subroutine print_wavederf_col_fbz_sort
      integer :: ik, is, ib1, ib2
      integer :: lun
      integer :: trev, iopr, ii, jk, iktmp, jksnl, iksnl
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      real(kind=DP), allocatable :: eko_wk2(:), occ_wk2(:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk2(:,:,:)

      integer istatus(mpi_status_size)

      lun = 1000
      if ( mype == 0 ) then
         open( lun, file="WAVEDERF_phase", status="unknown", form="formatted" )
         if ( sw_wavederf_full_bz == ON ) then
            write(lun,*) nspin, kv3_fbz/nspin, ib_num
         else
            write(lun,*) nspin, kv3/nspin, ib_num
         endif
      endif

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) );   zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif
!
      if ( mype == 0 ) then
         allocate( eko_wk2(neg) );   allocate( occ_wk2(neg) )
         allocate( zwk2(neg,neg,3) )
      endif

      Do is=1, nspin
         jkloop1: Do jk=is, kv3_fbz, nspin
            if ( myrank_e /= 0 ) cycle

            jksnl = (jk-1)/nspin + 1

            iopr = iopr_k_fbz_to_ibz(jk)
            trev = trev_k_fbz_to_ibz(jk)
            ikloop: Do ik=1, kv3, ndim_spinor
               Do ii=1, num_star_of_k(ik)
                  if ( jk == star_of_k(ik,ii) ) then
                     iktmp = ik;  exit ikloop
                  endif
               End Do
            end Do ikloop
            ik = iktmp
            iksnl = (ik-1)/nspin + 1

            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( eko_wk(:,ik), neg, mpi_double_precision, &
                    &         0, 1, &
                    &         MPI_CommGroup, ierr )
               call mpi_send( occ_wk(:,ik), neg, mpi_double_precision, &
                    &         0, 2, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  eko_wk2(:) = eko_wk(:,ik);     occ_wk2(:) = occ_wk(:,ik)
               else
                  call mpi_recv( eko_wk2(:), neg, mpi_double_precision, &
                       &         map_ek(1,ik), 1, &
                       &         MPI_CommGroup, istatus, ierr )
                  call mpi_recv( occ_wk2(:), neg, mpi_double_precision, &
                       &         map_ek(1,ik), 2, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif

            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( zwk(:,:,:,ik), 2*neg*neg*3, mpi_double_precision, &
                    &         0, 3, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  zwk2(:,:,:) = zwk(:,:,:,ik)
               else
                  call mpi_recv( zwk2(:,:,:), 2*neg*neg*3, mpi_double_precision, &
                       &         map_ek(1,ik), 3, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif

            if ( mype == 0 ) then
               weight = kv3 *qwgt(ik) /dble(ndim_spinor)

               Do ib1=ib_min, ib_max
                  Do ib2=ib_min, ib_max
                     occ1 = occ_wk2(ib1) / weight
                     occ2 = occ_wk2(ib2) / weight
                     zvec(1) = zwk2(ib1,ib2,1)
                     zvec(2) = zwk2(ib1,ib2,2)
                     zvec(3) = zwk2(ib1,ib2,3)

                     zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                     if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
#if 1
                     if ( allocated(magmom_dir_inversion_opr_flag) ) then
                        if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                           zvec_rot = -conjg(zvec_rot)
                        endif
                     endif
#endif
                     write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                          &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                          &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                          &       real(zvec_rot(1)), aimag(zvec_rot(1)), &
                          &       real(zvec_rot(2)), aimag(zvec_rot(2)), &
                          &       real(zvec_rot(3)), aimag(zvec_rot(3))
                  End Do
               End Do
            endif
         End Do jkloop1
      End Do

      if ( mype == 0 ) then
         deallocate( zwk2 );   deallocate( eko_wk2 );  deallocate( occ_wk2 )
      endif
      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_col_fbz_sort

    subroutine print_wavederf_col
      integer :: ik, is, ib1, ib2
      integer :: lun
      integer :: trev, iopr, ii, jk
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      real(kind=DP), allocatable :: eko_wk2(:), occ_wk2(:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk2(:,:,:)

      integer istatus(mpi_status_size)

      lun = 1000
      if ( mype == 0 ) then
         open( lun, file="WAVEDERF_phase", status="unknown", form="formatted" )
         if ( sw_wavederf_full_bz == ON ) then
            write(lun,*) nspin, kv3_fbz/nspin, ib_num
         else
            write(lun,*) nspin, kv3/nspin, ib_num
         endif
      endif

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) );   zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif
!
      if ( mype == 0 ) then
         allocate( eko_wk2(neg) );   allocate( occ_wk2(neg) )
         allocate( zwk2(neg,neg,3) )
      endif

      Do is=1, nspin
         ikloop1: Do ik=is, kv3, nspin
            if ( myrank_e /= 0 ) cycle

            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( eko_wk(:,ik), neg, mpi_double_precision, &
                    &         0, 1, &
                    &         MPI_CommGroup, ierr )
               call mpi_send( occ_wk(:,ik), neg, mpi_double_precision, &
                    &         0, 2, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  eko_wk2(:) = eko_wk(:,ik);     occ_wk2(:) = occ_wk(:,ik)
               else
                  call mpi_recv( eko_wk2(:), neg, mpi_double_precision, &
                       &         map_ek(1,ik), 1, &
                       &         MPI_CommGroup, istatus, ierr )
                  call mpi_recv( occ_wk2(:), neg, mpi_double_precision, &
                       &         map_ek(1,ik), 2, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif

            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( zwk(:,:,:,ik), 2*neg*neg*3, mpi_double_precision, &
                    &         0, 3, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  zwk2(:,:,:) = zwk(:,:,:,ik)
               else
                  call mpi_recv( zwk2(:,:,:), 2*neg*neg*3, mpi_double_precision, &
                       &         map_ek(1,ik), 3, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif

            if ( mype == 0 ) then
               if ( sw_wavederf_full_bz == ON ) then
                  weight = kv3 *qwgt(ik) /dble(ndim_spinor)

                  Do ii=1, num_star_of_k(ik)
                     jk = star_of_k(ik,ii)
                     iopr = iopr_k_fbz_to_ibz(jk)
                     trev = trev_k_fbz_to_ibz(jk)

                     Do ib1=ib_min, ib_max
                        Do ib2=ib_min, ib_max
                           occ1 = occ_wk2(ib1) / weight
                           occ2 = occ_wk2(ib2) / weight
                           zvec(1) = zwk2(ib1,ib2,1)
                           zvec(2) = zwk2(ib1,ib2,2)
                           zvec(3) = zwk2(ib1,ib2,3)

                           zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                           if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
#if 1
                           if ( allocated(magmom_dir_inversion_opr_flag) ) then
                              if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                                 zvec_rot = -conjg(zvec_rot)
                              endif
                           endif
#endif
                           write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                                &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                                &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                                &       real(zvec_rot(1)), aimag(zvec_rot(1)), &
                                &       real(zvec_rot(2)), aimag(zvec_rot(2)), &
                                &       real(zvec_rot(3)), aimag(zvec_rot(3))
                        End Do
                     End Do
                  End Do
               else
                  weight = kv3 *qwgt(ik) /dble(ndim_spinor)
                  Do ib1=ib_min, ib_max
                     Do ib2=ib_min, ib_max
                        occ1 = occ_wk2(ib1) / weight
                        occ2 = occ_wk2(ib2) / weight
                        zvec(1) = zwk2(ib1,ib2,1)
                        zvec(2) = zwk2(ib1,ib2,2)
                        zvec(3) = zwk2(ib1,ib2,3)
                        write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                             &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                             &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                             &       real(zvec(1)), aimag(zvec(1)), &
                             &       real(zvec(2)), aimag(zvec(2)), &
                             &       real(zvec(3)), aimag(zvec(3))
                     End Do
                  End Do
               end if
            endif
         End Do ikloop1
      End Do

      if ( mype == 0 ) then
         deallocate( zwk2 );   deallocate( eko_wk2 );  deallocate( occ_wk2 )
      endif
      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_col

    subroutine print_wavederf_col2    ! save_memory_mode
      integer :: ik, is, ib1, ib2
      integer :: lun
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:), occ_wk(:), eko_mpi(:), occ_mpi(:)
      complex(kind=CMPLDP) :: z1, z2, z3
      complex(kind=CMPLDP), allocatable :: zwk(:,:)

      integer istatus(mpi_status_size)

      lun = 1000
      if ( mype == 0 ) then
         open( lun, file="WAVEDERF_phase", status="unknown", form="formatted" )
         write(lun,*) nspin, kv3/nspin, ib_num
      endif

      allocate( eko_wk( neg ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg ) ); occ_wk = 0.0d0

      Do is=1, nspin
         ikloop1: Do ik=is, kv3, nspin
            eko_wk = 0.0d0;     occ_wk = 0.0d0

            if ( map_k(ik) == myrank_k ) then
               Do ib1=ista_e, iend_e, istep_e
                  eko_wk(ib1) = eko_l( map_z(ib1),ik )
                  occ_wk(ib1) = occup_l( map_z(ib1),ik )
               End do
               if ( npes > 1 ) then
                  allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
                  call mpi_allreduce( eko_wk, eko_mpi, neg, &
                       &              mpi_double_precision, mpi_sum, &
                       &              mpi_k_world(myrank_k), ierr )
                  eko_wk = eko_mpi
                  deallocate( eko_mpi )
                  allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
                  call mpi_allreduce( occ_wk, occ_mpi, neg, &
                       &              mpi_double_precision, mpi_sum, &
                       &              mpi_k_world(myrank_k), ierr )
                  occ_wk = occ_mpi
                  deallocate( occ_mpi )
               endif
            endif

            if ( map_ek(1,ik) == mype ) then
               if ( mype /= 0 ) then
                  call mpi_send( eko_wk, neg, mpi_double_precision, &
                       &         0, 1, MPI_CommGroup, ierr )
                  call mpi_send( occ_wk, neg, mpi_double_precision, &
                       &         0, 2, MPI_CommGroup, ierr )
               endif
            else if ( mype == 0 .and. map_ek(1,ik) /= 0 ) then
               call mpi_recv( eko_wk, neg, mpi_double_precision, &
                    &         map_ek(1,ik), 1, &
                    &         MPI_CommGroup, istatus, ierr )
               call mpi_recv( occ_wk, neg, mpi_double_precision, &
                    &         map_ek(1,ik), 2, &
                    &         MPI_CommGroup, istatus, ierr )
            endif

            allocate( zwk( neg, 3 ) )
            weight = kv3 *qwgt(ik) /dble(ndim_spinor)

            Do ib1=ib_min, ib_max
               if ( map_ek(ib1,ik) == mype ) then
                  zwk( 1:neg, 1:3 ) = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
                  if ( mype /= 0 ) then
                     call mpi_send( zwk(:,:), 2*neg*3, mpi_double_precision, &
                          &         0, 3, &
                          &         MPI_CommGroup, ierr )
                  endif
               else if ( mype == 0 .and. map_ek(ib1,ik) /= 0 ) then
                  call mpi_recv( zwk(:,:), 2*neg*3, mpi_double_precision, &
                       &         map_ek(ib1,ik), 3, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
               if ( mype == 0 ) then
                  Do ib2=ib_min, ib_max
                     occ1 = occ_wk(ib1) / weight
                     occ2 = occ_wk(ib2) / weight
                     z1 = zwk(ib2,1)
                     z2 = zwk(ib2,2)
                     z3 = zwk(ib2,3)
                     write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                          &       ib1, ( eko_wk(ib1)-efermi )*Hartree, occ1, &
                          &       ib2, ( eko_wk(ib2)-efermi )*Hartree, occ2, &
                          &       real(z1), aimag(z1), real(z2), aimag(z2), &
                          &       real(z3), aimag(z3)
                  End Do
               end if
            End Do
            deallocate( zwk )
         End Do ikloop1
      End Do

      if ( mype == 0 ) close(lun)
      deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_col2

    subroutine print_wavederf_noncl
      integer :: ik, is, ib1, ib2, iktmp
      integer :: lun, num, j, nr
      integer :: trev, iopr, ii, jk
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      real(kind=DP), allocatable :: eko_wk2(:), occ_wk2(:)
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk2(:,:,:,:)

      character*72 file0
      integer istatus(mpi_status_size)

      lun = 1000
      file0 = "WAVEDERF_phase"
      if ( mype == 0 ) then
         open( lun, file=file0, status="unknown", form="formatted" )
         if ( sw_wavederf_full_bz == ON ) then
            write(lun,*) nspin, kv3_fbz/nspin, ib_num
         else
            write(lun,*) nspin, kv3/nspin, ib_num
         endif
      endif

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) );   zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif
!
      if ( mype == 0 ) then
         allocate( eko_wk2(neg) );   allocate( occ_wk2(neg) )
         allocate( zwk2(neg,neg,3,ndim_spinor) )
      endif

      Do ik=1, kv3, ndim_spinor
         if ( myrank_e /= 0 ) cycle

         if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
            call mpi_send( eko_wk(:,ik), neg, mpi_double_precision, &
                 &         0, 1, &
                 &         MPI_CommGroup, ierr )
            call mpi_send( occ_wk(:,ik), neg, mpi_double_precision, &
                 &         0, 2, &
                 &         MPI_CommGroup, ierr )
         else if ( mype == 0 ) then
            if ( map_k(ik) == myrank_k ) then
               eko_wk2(:) = eko_wk(:,ik);     occ_wk2(:) = occ_wk(:,ik)
            else
               call mpi_recv( eko_wk2(:), neg, mpi_double_precision, &
                    &         map_ek(1,ik), 1, &
                    &         MPI_CommGroup, istatus, ierr )
               call mpi_recv( occ_wk2(:), neg, mpi_double_precision, &
                    &         map_ek(1,ik), 2, &
                    &         MPI_CommGroup, istatus, ierr )
            endif
         endif
         Do is=1, ndim_spinor
            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( zwk(:,:,:,ik+is-1), 2*neg*neg*3, mpi_double_precision, &
                    &         0, 3 +is, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  zwk2(:,:,:,is) = zwk(:,:,:,ik+is-1)
               else
                  call mpi_recv( zwk2(:,:,:,is), 2*neg*neg*3, mpi_double_precision, &
                       &         map_ek(1,ik), 3+is, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif
         End Do
!
         if ( mype == 0 ) then
            if ( sw_wavederf_full_bz == ON ) then
               weight = kv3 *qwgt(ik) /dble(ndim_spinor)

               Do ii=1, num_star_of_k(ik)
                  jk = star_of_k(ik,ii)
                  iopr = iopr_k_fbz_to_ibz(jk)
                  trev = trev_k_fbz_to_ibz(jk)

                  Do ib1=ib_min, ib_max
                     Do ib2=ib_min, ib_max
                        occ1 = occ_wk2(ib1) / weight
                        occ2 = occ_wk2(ib2) / weight
                        zvec(1) = zwk2(ib1,ib2,1,1) +zwk2(ib1,ib2,1,2)
                        zvec(2) = zwk2(ib1,ib2,2,1) +zwk2(ib1,ib2,2,2)
                        zvec(3) = zwk2(ib1,ib2,3,1) +zwk2(ib1,ib2,3,2)

                        zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                        if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
                        if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                           zvec_rot = -conjg(zvec_rot)
                        endif

                        write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                             &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                             &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                             &       real(zvec_rot(1)), aimag(zvec_rot(1)), &
                             &       real(zvec_rot(2)), aimag(zvec_rot(2)), &
                             &       real(zvec_rot(3)), aimag(zvec_rot(3))
                     End Do
                  End Do
               End Do
            else
               weight = kv3 *qwgt(ik) /dble(ndim_spinor)
               Do ib1=ib_min, ib_max
                  Do ib2=ib_min, ib_max
                     occ1 = occ_wk2(ib1) / weight
                     occ2 = occ_wk2(ib2) / weight
                     zvec(1) = zwk2(ib1,ib2,1,1) +zwk2(ib1,ib2,1,2)
                     zvec(2) = zwk2(ib1,ib2,2,1) +zwk2(ib1,ib2,2,2)
                     zvec(3) = zwk2(ib1,ib2,3,1) +zwk2(ib1,ib2,3,2)
                     write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                          &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                          &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                          &       real(zvec(1)), aimag(zvec(1)), &
                          &       real(zvec(2)), aimag(zvec(2)), &
                          &       real(zvec(3)), aimag(zvec(3))
                  End Do
               End Do
            endif
         end if
      End Do
      if ( mype == 0 ) then
         deallocate( zwk2 );   deallocate( eko_wk2 );  deallocate( occ_wk2 )
      endif
      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_noncl

    subroutine print_wavederf_noncl_fbz_sort
      integer :: ik, is, ib1, ib2, iktmp
      integer :: lun, num, j, nr
      integer :: trev, iopr, ii, jk, jksnl, iksnl
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      real(kind=DP), allocatable :: eko_wk2(:), occ_wk2(:)
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)
      complex(kind=CMPLDP), allocatable :: zwk2(:,:,:,:)

      character*72 file0
      integer istatus(mpi_status_size)

      lun = 1000
      file0 = "WAVEDERF_phase"
      if ( mype == 0 ) then
         open( lun, file=file0, status="unknown", form="formatted" )
         if ( sw_wavederf_full_bz == ON ) then
            write(lun,*) nspin, kv3_fbz/nspin, ib_num
         else
            write(lun,*) nspin, kv3/nspin, ib_num
         endif
      endif

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) );   zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif
!
      if ( mype == 0 ) then
         allocate( eko_wk2(neg) );   allocate( occ_wk2(neg) )
         allocate( zwk2(neg,neg,3,ndim_spinor) )
      endif

      Do jk=1, kv3_fbz, ndim_spinor
         if ( myrank_e /= 0 ) cycle

         jksnl = (jk-1)/nspin + 1

         iopr = iopr_k_fbz_to_ibz(jk)
         trev = trev_k_fbz_to_ibz(jk)
         ikloop: Do ik=1, kv3, ndim_spinor
            Do ii=1, num_star_of_k(ik)
               if ( jk == star_of_k(ik,ii) ) then
                  iktmp = ik;  exit ikloop
               endif
            End Do
         end Do ikloop
         ik = iktmp
         iksnl = (ik-1)/nspin + 1

         if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
            call mpi_send( eko_wk(:,ik), neg, mpi_double_precision, &
                 &         0, 1, &
                 &         MPI_CommGroup, ierr )
            call mpi_send( occ_wk(:,ik), neg, mpi_double_precision, &
                 &         0, 2, &
                 &         MPI_CommGroup, ierr )
         else if ( mype == 0 ) then
            if ( map_k(ik) == myrank_k ) then
               eko_wk2(:) = eko_wk(:,ik);     occ_wk2(:) = occ_wk(:,ik)
            else
               call mpi_recv( eko_wk2(:), neg, mpi_double_precision, &
                    &         map_ek(1,ik), 1, &
                    &         MPI_CommGroup, istatus, ierr )
               call mpi_recv( occ_wk2(:), neg, mpi_double_precision, &
                    &         map_ek(1,ik), 2, &
                    &         MPI_CommGroup, istatus, ierr )
            endif
         endif
         Do is=1, ndim_spinor
            if ( mype /= 0 .and. map_k(ik) == myrank_k ) then
               call mpi_send( zwk(:,:,:,ik+is-1), 2*neg*neg*3, mpi_double_precision, &
                    &         0, 3 +is, &
                    &         MPI_CommGroup, ierr )
            else if ( mype == 0 ) then
               if ( map_k(ik) == myrank_k ) then
                  zwk2(:,:,:,is) = zwk(:,:,:,ik+is-1)
               else
                  call mpi_recv( zwk2(:,:,:,is), 2*neg*neg*3, mpi_double_precision, &
                       &         map_ek(1,ik), 3+is, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            endif
         End Do
!
         if ( mype == 0 ) then
            weight = kv3 *qwgt(ik) /dble(ndim_spinor)
            Do ib1=ib_min, ib_max
               Do ib2=ib_min, ib_max
                  occ1 = occ_wk2(ib1) / weight
                  occ2 = occ_wk2(ib2) / weight
                  zvec(1) = zwk2(ib1,ib2,1,1) +zwk2(ib1,ib2,1,2)
                  zvec(2) = zwk2(ib1,ib2,2,1) +zwk2(ib1,ib2,2,2)
                  zvec(3) = zwk2(ib1,ib2,3,1) +zwk2(ib1,ib2,3,2)

                  zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                  if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
                  if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                     zvec_rot = -conjg(zvec_rot)
                  endif

                  write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                       &       ib1, ( eko_wk2(ib1)-efermi )*Hartree, occ1, &
                       &       ib2, ( eko_wk2(ib2)-efermi )*Hartree, occ2, &
                       &       real(zvec_rot(1)), aimag(zvec_rot(1)), &
                       &       real(zvec_rot(2)), aimag(zvec_rot(2)), &
                       &       real(zvec_rot(3)), aimag(zvec_rot(3))
               End Do
            End Do
         end if
      End Do
      if ( mype == 0 ) then
         deallocate( zwk2 );   deallocate( eko_wk2 );  deallocate( occ_wk2 )
      endif
      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_noncl_fbz_sort

    subroutine print_wavederf_noncl2      ! save_memory_mode: 1
      integer :: ik, is, ib1, ib2, iktmp
      integer :: lun
      real(kind=DP) :: occ1, occ2, weight
      real(kind=DP), allocatable :: eko_wk(:), occ_wk(:), eko_mpi(:), occ_mpi(:)
      complex(kind=CMPLDP) :: z1, z2, z3
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:)

      integer istatus(mpi_status_size)

      lun = 1000
      if ( mype == 0 ) then
         open( lun, file="WAVEDERF_phase", status="unknown", form="formatted" )
         write(lun,*) 1, kv3/nspin, ib_num
      endif

      allocate( eko_wk( neg ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg ) ); occ_wk = 0.0d0

      ikloop: Do ik=1, kv3, ndim_spinor
         eko_wk = 0.0d0;     occ_wk = 0.0d0

         if ( map_k(ik) == myrank_k ) then
            Do ib1=ista_e, iend_e, istep_e
               eko_wk(ib1) = eko_l( map_z(ib1),ik )
               occ_wk(ib1) = occup_l( map_z(ib1),ik )
            End do
            if ( npes > 1 ) then
               allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
               call mpi_allreduce( eko_wk, eko_mpi, neg, mpi_double_precision, mpi_sum, &
                    &              mpi_k_world(myrank_k), ierr )
               eko_wk = eko_mpi
               deallocate( eko_mpi )
               allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
               call mpi_allreduce( occ_wk, occ_mpi, neg, mpi_double_precision, mpi_sum, &
                    &              mpi_k_world(myrank_k), ierr )
               occ_wk = occ_mpi
               deallocate( occ_mpi )
            endif
         endif

         if ( map_ek(1,ik) == mype ) then
            if ( mype /= 0 ) then
               call mpi_send( eko_wk, neg, mpi_double_precision, &
                    &         0, 1, MPI_CommGroup, ierr )
               call mpi_send( occ_wk, neg, mpi_double_precision, &
                    &         0, 2, MPI_CommGroup, ierr )
            endif
         else if ( mype == 0 .and. map_ek(1,ik) /= 0 ) then
            call mpi_recv( eko_wk, neg, mpi_double_precision, &
                 &         map_ek(1,ik), 1, &
                 &         MPI_CommGroup, istatus, ierr )
            call mpi_recv( occ_wk, neg, mpi_double_precision, &
                 &         map_ek(1,ik), 2, &
                 &         MPI_CommGroup, istatus, ierr )
         endif

         weight = kv3 *qwgt(ik) /dble(ndim_spinor)
         allocate( zwk( neg, 3, ndim_spinor ) )

         Do ib1=ib_min, ib_max
            zwk = 0.0d0
            Do is=1, ndim_spinor
               iktmp = ik +is -1
               if ( map_ek(ib1,ik) == mype ) then
                  zwk( 1:neg, 1:3, is ) = trm_vv( map_z(ib1), 1:neg, iktmp, 1:3 ) *(-zi)
                  if ( mype /= 0 ) then
                     call mpi_send( zwk(:,:,is), 2*neg*3, mpi_double_precision, &
                          &         0, 3+is, &
                          &         MPI_CommGroup, ierr )
                  endif
               else if ( mype == 0 .and. map_ek(ib1,ik) /= 0 ) then
                  call mpi_recv( zwk(:,:,is), 2*neg*3, mpi_double_precision, &
                       &         map_ek(ib1,iktmp), 3+is, &
                       &         MPI_CommGroup, istatus, ierr )
               endif
            End Do
            if ( mype == 0 ) then
               Do ib2=ib_min, ib_max
                  occ1 = occ_wk(ib1) / weight
                  occ2 = occ_wk(ib2) / weight
                  z1 = zwk(ib2,1,1) +zwk(ib2,1,2)
                  z2 = zwk(ib2,2,1) +zwk(ib2,2,2)
                  z3 = zwk(ib2,3,1) +zwk(ib2,3,2)
                  write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                       &       ib1, ( eko_wk(ib1)-efermi )*Hartree, occ1, &
                       &       ib2, ( eko_wk(ib2)-efermi )*Hartree, occ2, &
                       &       real(z1), aimag(z1), real(z2), aimag(z2), &
                       &       real(z3), aimag(z3)
               End Do
            End if
         End Do
         deallocate( zwk )
      End Do ikloop

      if ( mype == 0 ) close(lun)
      deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_noncl2

    subroutine print_wavederf_noncl_split
      integer :: ik, is, ib1, ib2, iktmp
      integer :: lun, num, i, j, nr, ii, jk, iopr, trev
      real(kind=DP) :: occ1, occ2, weight
      complex(kind=CMPLDP) :: zvec(3), zvec_rot(3)
      real(kind=DP), allocatable :: eko_wk(:,:), occ_wk(:,:), eko_mpi(:,:), occ_mpi(:,:)
      complex(kind=CMPLDP), allocatable :: zwk(:,:,:,:), zwk_mpi(:,:,:,:)

      character*4 char1
      character*72 file0, file1
      integer istatus(mpi_status_size)

      allocate( eko_wk( neg,ista_k:iend_k ) ); eko_wk = 0.0d0
      allocate( occ_wk( neg,ista_k:iend_k ) ); occ_wk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1,ik) = eko_l( map_z(ib1),ik )
            occ_wk(ib1,ik) = occup_l( map_z(ib1),ik )
         End do
      End Do
      if ( nrank_e > 1 ) then
         allocate( eko_mpi( neg,ista_k:iend_k ) ); eko_mpi = 0.0d0
         call mpi_allreduce( eko_wk, eko_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         eko_wk = eko_mpi
         deallocate( eko_mpi )
         allocate( occ_mpi( neg,ista_k:iend_k ) ); occ_mpi = 0.0d0
         call mpi_allreduce( occ_wk, occ_mpi, neg*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         occ_wk = occ_mpi
         deallocate( occ_mpi )
      endif

      allocate( zwk( neg, neg, 3, ista_k:iend_k ) )
      zwk = 0.0d0

      Do ik=ista_k, iend_k
         Do ib1=ista_e, iend_e, istep_e
            zwk( ib1, 1:neg, 1:3, ik ) &
                 &  = trm_vv( map_z(ib1), 1:neg, ik, 1:3 ) *(-zi)
         End Do
      End Do
      if ( nrank_e > 1 ) then
         allocate( zwk_mpi( neg, neg, 3, ista_k:iend_k ) ); zwk_mpi = 0.0d0
         call mpi_allreduce( zwk, zwk_mpi, neg*neg*3*2*(iend_k-ista_k+1), &
              &              mpi_double_precision, mpi_sum, &
              &              mpi_k_world(myrank_k), ierr )
         zwk = zwk_mpi
         deallocate( zwk_mpi )
      endif

      lun = 10000 +myrank_k
      write(char1,'(I4.4)') myrank_k

      nr = nrank_k /num_wavederf_files_once +1
      num = 0

      iloop: Do i=1, nr
         Do j=1, num_wavederf_files_once
            num = num +1
            if ( num > npes )    goto 200
            if ( num /= myrank_k +1 ) cycle
            if ( myrank_e /= 0 )   cycle

            file0 = "WAVEDERF_phase"
            file1 = trim(adjustl(file0)) // '.' // char1

            if ( myrank_e == 0 ) then
               open( lun, file=file1, status="unknown", form="formatted" )
            endif
            if ( mype == 0 ) then
               if ( sw_wavederf_full_bz == ON ) then
                  write(lun,*) nspin, kv3_fbz/nspin, ib_num
               else
                  write(lun,*) nspin, kv3/nspin, ib_num
               endif
            endif

            Do ik=ista_k, iend_k, ndim_spinor

               if ( sw_wavederf_full_bz == ON ) then
                  weight = kv3 *qwgt(ik) /dble(ndim_spinor)
                  Do ii=1, num_star_of_k(ik)
                     jk = star_of_k(ik,ii)
                     iopr = iopr_k_fbz_to_ibz(jk)
                     trev = trev_k_fbz_to_ibz(jk)

                     Do ib1=ib_min, ib_max
                        Do ib2=ib_min, ib_max
                           occ1 = occ_wk(ib1,ik) / weight
                           occ2 = occ_wk(ib2,ik) / weight
                           zvec(1) = zwk(ib1,ib2,1,ik) +zwk(ib1,ib2,1,ik+1)
                           zvec(2) = zwk(ib1,ib2,2,ik) +zwk(ib1,ib2,2,ik+1)
                           zvec(3) = zwk(ib1,ib2,3,ik) +zwk(ib1,ib2,3,ik+1)

                           zvec_rot = matmul( op(:,:,invop(iopr)), zvec(:) )

                           if ( trev == 1 ) zvec_rot = conjg( zvec_rot )
                           if ( magmom_dir_inversion_opr_flag(iopr) == -1 ) then
                              zvec_rot = -conjg(zvec_rot)
                           endif

                           write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                                &    ib1, ( eko_wk(ib1,ik)-efermi )*Hartree, occ1, &
                                &    ib2, ( eko_wk(ib2,ik)-efermi )*Hartree, occ2, &
                                &    real(zvec_rot(1)), aimag(zvec_rot(1)), &
                                &    real(zvec_rot(2)), aimag(zvec_rot(2)), &
                                &    real(zvec_rot(3)), aimag(zvec_rot(3))
                        End Do
                     End Do
                  End do

               else
                  weight = kv3 *qwgt(ik) /dble(ndim_spinor)
                  Do ib1=ib_min, ib_max
                     Do ib2=ib_min, ib_max
                        occ1 = occ_wk(ib1,ik) / weight
                        occ2 = occ_wk(ib2,ik) / weight
                        zvec(1) = zwk(ib1,ib2,1,ik) +zwk(ib1,ib2,1,ik+1)
                        zvec(2) = zwk(ib1,ib2,2,ik) +zwk(ib1,ib2,2,ik+1)
                        zvec(3) = zwk(ib1,ib2,3,ik) +zwk(ib1,ib2,3,ik+1)
                        write(lun,'(I5,2F15.8,I5,2F15.8,6F20.14)') &
                             &    ib1, ( eko_wk(ib1,ik)-efermi )*Hartree, occ1, &
                             &    ib2, ( eko_wk(ib2,ik)-efermi )*Hartree, occ2, &
                             &    real(zvec(1)), aimag(zvec(1)), &
                             &    real(zvec(2)), aimag(zvec(2)), &
                             &    real(zvec(3)), aimag(zvec(3))
                     End Do
                  End Do
               endif
            End Do

            if ( myrank_e == 0 ) close(lun)
         End Do
200      continue
         call mpi_barrier( MPI_CommGroup, ierr )
      End Do iloop

      deallocate( zwk );    deallocate( eko_wk );   deallocate( occ_wk )

    end subroutine print_wavederf_noncl_split

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, ig
      real(kind=DP), allocatable :: qx(:), qy(:), qz(:)
      real(kind=DP), allocatable :: wk_zaj(:,:)
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2, zsum(3)

      allocate( wk_zaj( kg1, kimg ) ); wk_zaj = 0.0d0

      allocate( qx(kg1) );  allocate( qy(kg1) ); allocate( qz(kg1) );
      qx = 0.0d0;  qy = 0.0d0;  qz = 0.0d0

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, vkxyz, ngabc, &
              &                   rltv, qx, qy, qz )

         Do ib2=1, neg
            wk_zaj = 0.0d0
            if ( map_e(ib2) == myrank_e ) then
               wk_zaj(1:iba(ik),1:kimg) = zaj_l(1:iba(ik),map_z(ib2),ik,1:kimg)
            endif
            call mpi_bcast( wk_zaj, kg1*kimg, mpi_double_precision, map_e(ib2), &
                 &          mpi_k_world(myrank_k), ierr )

            Do ib1=ista_e, iend_e, istep_e
               zsum = 0.0d0
               if ( kimg == 1 ) then
                  Do ig=1, kg1
                     c1 = zaj_l(ig,map_z(ib1),ik,1)
                     c2 = c1 *wk_zaj(ig,1)
                     zsum(1) = zsum(1) +c2 *qx(ig)
                     zsum(2) = zsum(2) +c2 *qy(ig)
                     zsum(3) = zsum(3) +c2 *qz(ig)
                  End do
               else
                  Do ig=1, kg1
                     z1 = dcmplx( zaj_l(ig,map_z(ib1),ik,1), zaj_l(ig,map_z(ib1),ik,2) )
                     z2 = conjg(z1) *dcmplx( wk_zaj(ig,1), wk_zaj(ig,2) )
                     zsum(1) = zsum(1) +z2 *qx(ig)
                     zsum(2) = zsum(2) +z2 *qy(ig)
                     zsum(3) = zsum(3) +z2 *qz(ig)
                  End do
               endif
               trm_vv( map_z(ib1), ib2, ik, 1:3 ) = zsum(1:3)
            End do
         End do
      End Do

      trm_vv = trm_vv *zi

      deallocate( wk_zaj )
      deallocate( qx ); deallocate( qy ); deallocate( qz )

    end subroutine calc_soft_part

    subroutine calc_hard_part
      integer :: ia, it, ik, ib1, ib2
      integer :: lmt1, lmt2, il1, il2, im1, im2, it1, it2, itmp
      integer :: lmta1, lmta2
      real(kind=DP) :: fac
      complex(kind=CMPLDP) :: z1, z2, wf1, wf2
      complex(kind=CMPLDP), allocatable :: wk_fsri(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri_add(:)

      allocate( wk_fsri(nlmta) ); wk_fsri = 0.0d0
      if ( sw_use_add_proj == ON ) then
         allocate( wk_fsri_add(nlmta_add) ); wk_fsri_add = 0.0d0
      endif

      Do ia=1, natm
         it = ityp(ia)

         Do ik=1, kv3
            if ( map_k(ik) /= myrank_k ) cycle

            Do ib2=1, neg
               wk_fsri = 0.0d0
               if ( map_e(ib2) == myrank_e ) then
                  wk_fsri(:) = dcmplx( fsr_l( map_z(ib2),:,ik ), &
                       &               fsi_l( map_z(ib2),:,ik ) )
               endif
               call mpi_bcast( wk_fsri, 2*nlmta, mpi_double_precision, map_e(ib2), &
                    &          mpi_k_world(myrank_k), ierr )

               if ( sw_use_add_proj == ON ) then
                  wk_fsri_add = 0.0d0
                  if ( map_e(ib2) == myrank_e ) then
                     wk_fsri_add(:) = dcmplx( fsr_add_l( map_z(ib2),:,ik ), &
                          &                   fsi_add_l( map_z(ib2),:,ik ) )
                  endif
                  call mpi_bcast( wk_fsri_add, 2*nlmta_add, mpi_double_precision, &
                       &          map_e(ib2), mpi_k_world(myrank_k), ierr )
               endif

               DO ib1=ista_e, iend_e, istep_e

                  Do lmt1=1, ilmt_val(it)
                     if ( lmt1 <= ilmt(it) ) then
                        il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                        lmta1 = lmta( lmt1,ia )
                        wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_l( map_z(ib1), lmta1, ik ) )
                     else
                        itmp = lmt1 -ilmt(it)

                        il1 = ltp_add(itmp,it); it1 = 1
                        lmta1 = lmta_add( itmp,ia )
                        wf1 = dcmplx( fsr_add_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_add_l( map_z(ib1), lmta1, ik ) )
                     endif

                     Do lmt2=1, ilmt_val(it)
                        fac = 1.0d0

                        if ( lmt2 <= ilmt(it) ) then
                           il2 = ltp(lmt2,it); it2 = taup(lmt2,it)
                           lmta2 = lmta( lmt2,ia )
                           wf2 = wk_fsri(lmta2)
                        else
                           itmp = lmt2 -ilmt(it)

                           il2 = ltp_add( itmp,it); it2 = 1
                           lmta2 = lmta_add( itmp,ia )
                           wf2 = wk_fsri_add(lmta2)
                        endif

                        z1 = conjg(wf1) *wf2

                        trm_vv( map_z(ib1), ib2, ik, 1:3 ) &
                             &  = trm_vv( map_z(ib1), ib2, ik, 1:3 ) &
                             &    + z1 *Mat_Dipole_Corr_vv( it,lmt1,lmt2,1:3 )
                     End do
                  End Do
               End Do
            End DO
         End Do
      End Do

      deallocate( wk_fsri )
      if ( allocated( wk_fsri_add ) ) deallocate( wk_fsri_add )

    end subroutine calc_hard_part

  end subroutine m_XI_wd_WAVEDERF

  subroutine m_XI_calc_transition_moment_vv
    if ( allocated( trm_vv ) ) deallocate( trm_vv )
    allocate( trm_vv( np_e, neg, ista_k:iend_k, 3 ) );   trm_vv = 0.0d0

    call calc_soft_part
    call calc_hard_part
    call multiply_factor

  contains

    subroutine multiply_factor
      integer :: ik, ib1, ib2, is
      real(kind=DP) :: ediff
      real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)

      allocate( eko_wk( neg ) ); eko_wk = 0.0d0

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

         eko_wk = 0.0d0
         Do ib1=ista_e, iend_e, istep_e
            eko_wk(ib1) = eko_l( map_z(ib1),ik )
         End do
         if ( npes > 1 ) then
            allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
            call mpi_allreduce( eko_wk, eko_mpi, neg, mpi_double_precision, mpi_sum, &
                 &              mpi_k_world(myrank_k), ierr )
            eko_wk = eko_mpi
            deallocate( eko_mpi )
         endif

         Do ib2=1, neg
            Do ib1=ista_e, iend_e, istep_e
               ediff = eko_wk(ib2) -eko_l(map_z(ib1),ik)
!
               if ( ib2 == ib1 ) then
                  ediff = delta_omega
               else if ( abs(ediff) < delta_omega ) then
                  if ( ediff >= 0.0 ) then
                     ediff = delta_omega
                  else
                     ediff = -delta_omega
                  endif
               else
                  if ( sw_scissor_renormalization == ON ) then
                     if ( eko_wk(ib2) > efermi .and. eko_l(map_z(ib1),ik ) <= efermi ) then
                        ediff = ediff +scissor
                     endif
                     if ( eko_wk(ib2) <= efermi .and. eko_l(map_z(ib1),ik ) > efermi ) then
                        ediff = ediff -scissor
                     endif
                  endif
               endif

               Do is=1, ndim_spinor
                  trm_vv( map_z(ib1),ib2,ik+is-1,: ) &
                       &   = trm_vv( map_z(ib1),ib2,ik+is-1,: ) /ediff
               End do

            End Do
        End Do
      End Do

      deallocate( eko_wk )

    end subroutine multiply_factor

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, ig
      real(kind=DP), allocatable :: qx(:), qy(:), qz(:)
      real(kind=DP), allocatable :: wk_zaj(:,:)
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2, zsum(3)

      allocate( wk_zaj( kg1, kimg ) ); wk_zaj = 0.0d0

      allocate( qx(kg1) );  allocate( qy(kg1) ); allocate( qz(kg1) );
      qx = 0.0d0;  qy = 0.0d0;  qz = 0.0d0

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, vkxyz, ngabc, &
              &                   rltv, qx, qy, qz )

         Do ib2=1, neg
            wk_zaj = 0.0d0
            if ( map_e(ib2) == myrank_e ) then
               wk_zaj(1:iba(ik),1:kimg) = zaj_l(1:iba(ik),map_z(ib2),ik,1:kimg)
            endif
            call mpi_bcast( wk_zaj, kg1*kimg, mpi_double_precision, map_e(ib2), &
                 &          mpi_k_world(myrank_k), ierr )

            Do ib1=ista_e, iend_e, istep_e
               zsum = 0.0d0
               if ( kimg == 1 ) then
                  Do ig=1, kg1
                     c1 = zaj_l(ig,map_z(ib1),ik,1)
                     c2 = c1 *wk_zaj(ig,1)
                     zsum(1) = zsum(1) +c2 *qx(ig)
                     zsum(2) = zsum(2) +c2 *qy(ig)
                     zsum(3) = zsum(3) +c2 *qz(ig)
                  End do
               else
                  Do ig=1, kg1
                     z1 = dcmplx( zaj_l(ig,map_z(ib1),ik,1), zaj_l(ig,map_z(ib1),ik,2) )
                     z2 = conjg(z1) *dcmplx( wk_zaj(ig,1), wk_zaj(ig,2) )
                     zsum(1) = zsum(1) +z2 *qx(ig)
                     zsum(2) = zsum(2) +z2 *qy(ig)
                     zsum(3) = zsum(3) +z2 *qz(ig)
                  End do
               endif
               trm_vv( map_z(ib1), ib2, ik, 1:3 ) = zsum(1:3)
            End do
         End do
      End Do

      trm_vv = trm_vv *zi

      deallocate( wk_zaj )
      deallocate( qx ); deallocate( qy ); deallocate( qz )

    end subroutine calc_soft_part

    subroutine calc_hard_part
      integer :: ia, it, ik, ib1, ib2
      integer :: lmt1, lmt2, il1, il2, im1, im2, it1, it2, itmp
      integer :: lmta1, lmta2
      real(kind=DP) :: fac
      complex(kind=CMPLDP) :: z1, z2, wf1, wf2
      complex(kind=CMPLDP), allocatable :: wk_fsri(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri_add(:)

      allocate( wk_fsri(nlmta) ); wk_fsri = 0.0d0
      if ( sw_use_add_proj == ON ) then
         allocate( wk_fsri_add(nlmta_add) ); wk_fsri_add = 0.0d0
      endif

      Do ia=1, natm
         it = ityp(ia)

         Do ik=1, kv3
            if ( map_k(ik) /= myrank_k ) cycle

            Do ib2=1, neg
               wk_fsri = 0.0d0
               if ( map_e(ib2) == myrank_e ) then
                  wk_fsri(:) = dcmplx( fsr_l( map_z(ib2),:,ik ), &
                       &               fsi_l( map_z(ib2),:,ik ) )
               endif
               call mpi_bcast( wk_fsri, 2*nlmta, mpi_double_precision, map_e(ib2), &
                    &          mpi_k_world(myrank_k), ierr )

               if ( sw_use_add_proj == ON ) then
                  wk_fsri_add = 0.0d0
                  if ( map_e(ib2) == myrank_e ) then
                     wk_fsri_add(:) = dcmplx( fsr_add_l( map_z(ib2),:,ik ), &
                          &                   fsi_add_l( map_z(ib2),:,ik ) )
                  endif
                  call mpi_bcast( wk_fsri_add, 2*nlmta_add, mpi_double_precision, &
                       &          map_e(ib2), mpi_k_world(myrank_k), ierr )
               endif

               DO ib1=ista_e, iend_e, istep_e

                  Do lmt1=1, ilmt_val(it)
                     if ( lmt1 <= ilmt(it) ) then
                        il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                        lmta1 = lmta( lmt1,ia )
                        wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_l( map_z(ib1), lmta1, ik ) )
                     else
                        itmp = lmt1 -ilmt(it)

                        il1 = ltp_add(itmp,it); it1 = 1
                        lmta1 = lmta_add( itmp,ia )
                        wf1 = dcmplx( fsr_add_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_add_l( map_z(ib1), lmta1, ik ) )
                     endif

                     Do lmt2=1, ilmt_val(it)
                        fac = 1.0d0

                        if ( lmt2 <= ilmt(it) ) then
                           il2 = ltp(lmt2,it); it2 = taup(lmt2,it)
                           lmta2 = lmta( lmt2,ia )
                           wf2 = wk_fsri(lmta2)
                        else
                           itmp = lmt2 -ilmt(it)

                           il2 = ltp_add( itmp,it); it2 = 1
                           lmta2 = lmta_add( itmp,ia )
                           wf2 = wk_fsri_add(lmta2)
                        endif

                        z1 = conjg(wf1) *wf2

                        trm_vv( map_z(ib1), ib2, ik, 1:3 ) &
                             &  = trm_vv( map_z(ib1), ib2, ik, 1:3 ) &
                             &    + z1 *Mat_Dipole_Corr_vv( it,lmt1,lmt2,1:3 )
                     End do
                  End Do
               End Do
            End DO
         End Do
      End Do

      deallocate( wk_fsri )
      if ( allocated( wk_fsri_add ) ) deallocate( wk_fsri_add )

    end subroutine calc_hard_part

  end subroutine m_XI_calc_transition_moment_vv

  subroutine m_XI_set_RhoTilde_vv

    allocate( RhoTilde_vv( nmax_G, np_e, neg, ista_k:iend_k ) )
    RhoTilde_vv = 0.0d0

    call calc_soft_part
    call calc_hard_part

  contains

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, i, i1, ngrid
      real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)
      real(kind=DP), allocatable :: psi_l(:,:,:,:)

      allocate( afft(nfft) ); afft = 0.0d0
      allocate( bfft(nfft) ); bfft = 0.0d0
      allocate( cfft(nfft) ); cfft = 0.0d0

      call m_FFT_alloc_WF_work()
! --
      ngrid = product(fft_box_size_WF(1:3,1))

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         allocate( psi_l( kg1, np_e, ik:ik, kimg ) ); psi_l = 0.0d0

         Do ib2=1, neg
            bfft = 0.0d0
            if ( map_e(ib2) == myrank_e ) then
               psi_l(1:iba(ik),map_z(ib2),ik,:) = zaj_l(1:iba(ik),map_z(ib2),ik,:)
               call m_ES_WF_in_Rspace1( ik, ik, ik, ib2, psi_l, bfft )
            endif
            call mpi_bcast( bfft, nfft, mpi_double_precision, map_e(ib2), &
                 &          mpi_k_world(myrank_k), ierr )

            Do ib1=ista_e, iend_e, istep_e
               psi_l(1:iba(ik),map_z(ib1),ik,:) = zaj_l(1:iba(ik),map_z(ib1),ik,:)
               call m_ES_WF_in_Rspace1( ik, ik, ik, ib1, psi_l, afft )

               cfft = 0.0d0
               if ( imple_method == 1 ) then
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) =-afft(i)*bfft(i+1) +afft(i+1)*bfft(i)
                  End Do
               else
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) = afft(i)*bfft(i+1) -afft(i+1)*bfft(i)
                  End do
               endif

               call m_FFT_WF( ELECTRON, nfout, cfft, DIRECT, OFF )    ! R-->G
               cfft = cfft /dble(ngrid)

               if ( kimg == 1 ) then
                  Do i=1, nmax_G
                     i1 = igf( nbase(i,ik) )
                     RhoTilde_vv( i, map_z(ib1), ib2, ik ) = cfft(i1)
                  End do
               else
                  if ( imple_method == 2  ) then
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_vv( i, map_z(ib1), ib2, ik ) &
                             &          = dcmplx( cfft(i1), -cfft(i1+1) )
                     End do
                  else
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_vv( i, map_z(ib1), ib2, ik ) &
                             &          = dcmplx( cfft(i1), cfft(i1+1) )
                     End do
                  endif
               endif
            End Do
         End Do

         deallocate( psi_l )
      End Do

      deallocate( afft ); deallocate( bfft ); deallocate( cfft )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_soft_part

    subroutine calc_hard_part
      integer :: ik, ib1, ib2, ia, it, ig
      integer :: lmt1, lmt2, it1, it2, il1, il2, n
      integer :: lmta1, lmta2, ilm3, l3, iiqitg, mdvdb
      real(kind=DP) :: fac
      complex(kind=CMPLDP) :: wf1, wf2, z1, z2, zph

      integer, allocatable :: il3(:)
      real(kind=DP), allocatable :: zfsin(:), zfcos(:), qx(:), qy(:), qz(:), ylm(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri(:), wk_fsri_add(:)

      integer :: j, itmp

! -- init
      allocate( zfcos( nmax_G ) ); zfcos = 0.0d0
      allocate( zfsin( nmax_G ) ); zfsin = 0.0d0

      n = nloc
      n=(n-1)+(n-1)+1

      allocate(il3(n**2));call substitute_il3(n**2,il3)

      allocate( qx(nmax_G) ); allocate( qy(nmax_G) ); allocate( qz(nmax_G) )
      allocate( ylm(nmax_G) )

      Do ig=1, nmax_G
         qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
         qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
         qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
      End do

      allocate( wk_fsri(nlmta) ); wk_fsri = 0.0d0
      if ( sw_use_add_proj == ON ) then
         allocate( wk_fsri_add(nlmta_add) ); wk_fsri_add = 0.0d0
      endif

! -- start
      Do ia=1, natm
         it = ityp(ia)

         call calc_phase2(natm,pos,ia,kgp,ngabc,1,nmax_G,zfcos,zfsin)

         Do ik=1, kv3
            if ( map_k(ik) /= myrank_k ) cycle

            Do ib2=1, neg
               wk_fsri = 0.0d0
               if ( map_e(ib2) == myrank_e ) then
                  wk_fsri(:) = dcmplx( fsr_l( map_z(ib2),:,ik ), &
                       &               fsi_l( map_z(ib2),:,ik ) )
               endif
               call mpi_bcast( wk_fsri, 2*nlmta, mpi_double_precision, map_e(ib2), &
                    &          mpi_k_world(myrank_k), ierr )

               if ( sw_use_add_proj == ON ) then
                  wk_fsri_add = 0.0d0
                  if ( map_e(ib2) == myrank_e ) then
                     wk_fsri_add(:) = dcmplx( fsr_add_l( map_z(ib2),:,ik ), &
                          &                   fsi_add_l( map_z(ib2),:,ik ) )
                  endif
                  call mpi_bcast( wk_fsri_add, 2*nlmta_add, mpi_double_precision, &
                       &          map_e(ib2), mpi_k_world(myrank_k), ierr )
               endif

               DO ib1=ista_e, iend_e, istep_e

                  Do lmt1=1, ilmt_val(it)
                     if ( lmt1 <= ilmt(it) ) then
                        il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                        lmta1 = lmta( lmt1,ia )
                        wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_l( map_z(ib1), lmta1, ik ) )
                     else
                        itmp = lmt1 -ilmt(it)

                        il1 = ltp_add(itmp,it); it1 = 1
                        lmta1 = lmta_add( itmp,ia )
                        wf1 = dcmplx( fsr_add_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_add_l( map_z(ib1), lmta1, ik ) )
                     endif

!#if 0
!                     Do lmt2=lmt1, ilmt(it)
!#else
                     Do lmt2=1, ilmt_val(it)
!#endif
                        if ( lmt2 <= ilmt(it) ) then
                           il2 = ltp(lmt2,it); it2 = taup(lmt2,it)
                           lmta2 = lmta( lmt2,ia )
                           wf2 = wk_fsri(lmta2)
                        else
                           itmp = lmt2 -ilmt(it)

                           il2 = ltp_add( itmp,it); it2 = 1
                           lmta2 = lmta_add( itmp,ia )
                           wf2 = wk_fsri_add(lmta2)
                        endif

!#if 0
!                        if ( lmt1 == lmt2 ) then
!                           fac = 1.0d0
!                        else
!                           fac = 2.0d0
!                        endif
!#else
                        fac = 1.0d0
!#endif
                        do n=1,il2p_XI_vv(lmt1,lmt2,it)

                           ilm3 = isph_XI_vv(lmt1,lmt2,n,it); l3=il3(ilm3)
                           iiqitg = iqitg_XI_vv(il1,it1,il2,it2,l3+1,it)

                           if(iiqitg == 0) cycle

                           call sphr( nmax_G, ilm3, qx, qy, qz, ylm )

!#if 0
!                           z1 = conjg(wf1) *wf2
!                           z2 = wf1 *conjg(wf2)
!                           z1 = ( z1+z2) /2.0d0 *zi**l3 *dl2p(lmt1,lmt2,n,it)
!#else
                           z1 = conjg(wf1) *wf2 *zi**l3 *dl2p_XI_vv(lmt1,lmt2,n,it)
!#endif
                           if ( imple_method == 1 ) then
                              Do ig=1, nmax_G
                                 z2 = qitg_XI_vv( ig,iiqitg ) *ylm(ig)
                                 zph = dcmplx( zfcos(ig), -zfsin(ig) )

                                 RhoTilde_vv( ig, map_z(ib1), ib2, ik ) &
                                      & = RhoTilde_vv( ig, map_z(ib1), ib2, ik ) &
                                      &   + fac *z1 *z2 *zph
                              End do
                           else
                              Do ig=1, nmax_G
                                 z2 = qitg_XI_vv( ig,iiqitg ) *ylm(ig)
                                 zph = dcmplx( zfcos(ig), zfsin(ig) )

                                 RhoTilde_vv( ig, map_z(ib1), ib2, ik ) &
                                      & = RhoTilde_vv( ig, map_z(ib1), ib2, ik ) &
                                      &   + fac *z1 *z2 *zph
                              End do
                           endif

                        End do
                     End do
                 End do
               End do
            End Do
         End DO
      End Do

      deallocate( il3 ); deallocate( ylm )
      deallocate( qx ); deallocate( qy ); deallocate( qz );
      deallocate( zfcos ); deallocate( zfsin )

      deallocate( wk_fsri )
      if ( allocated( wk_fsri_add ) ) deallocate( wk_fsri_add )

    end subroutine calc_hard_part

  end subroutine m_XI_set_RhoTilde_vv

! --------------------------------------------------------
!    Calculation of transtion moment/matrix (val-core)
! --------------------------------------------------------

  subroutine m_XI_calc_transition_moment_vc

    if ( allocated( trm_vc ) ) deallocate( trm_vc )
    allocate( trm_vc( np_e, num_core_states, ista_k:iend_k, 3 ) ); trm_vc = 0.0d0

    call calc_soft_part
    call calc_hard_part
    call multiply_factor

  contains

    subroutine multiply_factor
      integer :: ik, ib1, ib2, is
      real(kind=DP) :: ediff
      real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)

      allocate( eko_wk( num_core_states ) ); eko_wk = 0.0d0

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

         eko_wk = ene_core_states

         Do ib2=1, num_core_states
            Do ib1=ista_e, iend_e, istep_e
               ediff = eko_wk(ib2) -eko_l(map_z(ib1),ik)

               if ( abs(ediff) < delta_omega ) then
                  if ( ediff >= 0.0 ) then
                     ediff = delta_omega
                  else
                     ediff = -delta_omega
                  endif
               else
                  if ( sw_scissor_renormalization == ON ) then
                     if ( eko_l(map_z(ib1),ik ) > efermi ) then
                        ediff = ediff -scissor
                     endif
                  endif
               endif

               Do is=1, ndim_spinor
                  trm_vc( map_z(ib1),ib2,ik+is-1,: ) &
                       &   = trm_vc( map_z(ib1),ib2,ik+is-1,: ) /ediff
               End do

            End Do
         End Do
      End Do

      deallocate( eko_wk )

    end subroutine multiply_factor

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, ig
      real(kind=DP), allocatable :: qx(:), qy(:), qz(:)
      real(kind=DP), allocatable :: wk_zaj(:,:)
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2, zsum(3)

      allocate( wk_zaj( kg1, kimg ) ); wk_zaj = 0.0d0

      allocate( qx(kg1) );  allocate( qy(kg1) ); allocate( qz(kg1) );
      qx = 0.0d0;  qy = 0.0d0;  qz = 0.0d0

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, vkxyz, ngabc, &
              &                   rltv, qx, qy, qz )

         Do ib2=1, num_core_states
            wk_zaj(1:iba(ik),1:kimg) = psig_core_states(1:iba(ik),ib2,ik,1:kimg)

            Do ib1=ista_e, iend_e, istep_e
               zsum = 0.0d0
               if ( kimg == 1 ) then
                  Do ig=1, kg1
                     c1 = zaj_l(ig,map_z(ib1),ik,1)
                     c2 = c1 *wk_zaj(ig,1)
                     zsum(1) = zsum(1) +c2 *qx(ig)
                     zsum(2) = zsum(2) +c2 *qy(ig)
                     zsum(3) = zsum(3) +c2 *qz(ig)
                  End do
               else
                  Do ig=1, kg1
                     z1 = dcmplx( zaj_l(ig,map_z(ib1),ik,1), zaj_l(ig,map_z(ib1),ik,2) )
                     z2 = conjg(z1) *dcmplx( wk_zaj(ig,1), wk_zaj(ig,2) )
                     zsum(1) = zsum(1) +z2 *qx(ig)
                     zsum(2) = zsum(2) +z2 *qy(ig)
                     zsum(3) = zsum(3) +z2 *qz(ig)
                  End do
               endif
               trm_vc( map_z(ib1), ib2, ik, 1:3 ) = zsum(1:3)
            End do
         End do
      End Do

      trm_vc = trm_vc *zi

      deallocate( wk_zaj )
      deallocate( qx ); deallocate( qy ); deallocate( qz )

    end subroutine calc_soft_part

    subroutine calc_hard_part
      integer :: ia, it, ik, ib1, ib2
      integer :: lmt1, lmt2, il1, il2, im1, im2, it1, it2, itmp
      integer :: lmta1, lmta2, atomtype_to_probe
      real(kind=DP) :: fac
      complex(kind=CMPLDP) :: z1, z2, wf1, wf2
      complex(kind=CMPLDP), allocatable :: wk_fsri(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri_add(:)

      allocate( wk_fsri(nlmt_core) ); wk_fsri = 0.0d0

      atomtype_to_probe = ityp( atom_to_probe )
      it = atomtype_to_probe

      Do ia = atom_to_probe, atom_to_probe
         it = ityp(ia)

         Do ik=1, kv3
            if ( map_k(ik) /= myrank_k ) cycle

            Do ib2=1, num_core_states
               wk_fsri(:) = dcmplx( fsr_core_states( ib2,:,ik ), &
                    &               fsi_core_states( ib2,:,ik ) )

               DO ib1=ista_e, iend_e, istep_e

                  Do lmt1=1, ilmt_val(it)
                     if ( lmt1 <= ilmt(it) ) then
                        il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                        lmta1 = lmta( lmt1,ia )
                        wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_l( map_z(ib1), lmta1, ik ) )
                     else
                        itmp = lmt1 -ilmt(it)

                        il1 = ltp_add(itmp,it); it1 = 1
                        lmta1 = lmta_add( itmp,ia )
                        wf1 = dcmplx( fsr_add_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_add_l( map_z(ib1), lmta1, ik ) )
                     endif

                     Do lmt2=1, nlmt_core
                        fac = 1.0d0

                        wf2 = wk_fsri(lmt2)
                        z1 = conjg(wf1) *wf2

                        trm_vc( map_z(ib1), ib2, ik, 1:3 ) &
                             &  = trm_vc( map_z(ib1), ib2, ik, 1:3 ) &
                             &    + z1 *Mat_Dipole_Corr_vc( lmt1,lmt2,1:3 )
                     End do
                  End Do
               End Do
            End DO
         End Do
      End Do

      deallocate( wk_fsri )

    end subroutine calc_hard_part

  end subroutine m_XI_calc_transition_moment_vc

  subroutine m_XI_set_RhoTilde_vc

    allocate( RhoTilde_vc( nmax_G, np_e, num_core_states, ista_k:iend_k ) )
    RhoTilde_vc = 0.0d0

    call calc_soft_part
    call calc_hard_part

  contains

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, i, i1, ngrid
      real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)
      real(kind=DP), allocatable :: psi_l(:,:,:,:)

      allocate( afft(nfft) ); afft = 0.0d0
      allocate( bfft(nfft) ); bfft = 0.0d0
      allocate( cfft(nfft) ); cfft = 0.0d0

      call m_FFT_alloc_WF_work()
! --
      ngrid = product(fft_box_size_WF(1:3,1))

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         allocate( psi_l( kg1, 1, ik:ik, kimg ) ); psi_l = 0.0d0

         Do ib2=1, num_core_states
            bfft = 0.0d0
            psi_l(1:iba(ik),1,ik,:) = psig_core_states(1:iba(ik),ib2,ik,:)
            call m_ES_WF_in_Rspace_kt( ik, ik, ik, psi_l, bfft )

            Do ib1=ista_e, iend_e, istep_e
               psi_l(1:iba(ik),1,ik,:) = zaj_l(1:iba(ik),map_z(ib1),ik,:)
               call m_ES_WF_in_Rspace_kt( ik, ik, ik, psi_l, afft )

               cfft = 0.0d0
               if ( imple_method == 1 ) then
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) =-afft(i)*bfft(i+1) +afft(i+1)*bfft(i)
                  End Do
               else
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) = afft(i)*bfft(i+1) -afft(i+1)*bfft(i)
                  End do
               endif

               call m_FFT_WF( ELECTRON, nfout, cfft, DIRECT, OFF )    ! R-->G
               cfft = cfft /dble(ngrid)

               if ( kimg == 1 ) then
                  Do i=1, nmax_G
                     i1 = igf( nbase(i,ik) )
                     RhoTilde_vc( i, map_z(ib1), ib2, ik ) = cfft(i1)
                  End do
               else
                  if ( imple_method == 1 ) then
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_vc( i, map_z(ib1), ib2, ik ) &
                             &          = dcmplx( cfft(i1), -cfft(i1+1) )  ! gpaw
                     End do
                  else
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_vc( i, map_z(ib1), ib2, ik ) &
                             &          = dcmplx( cfft(i1), cfft(i1+1) )
                     End do
                  endif
               endif
            End Do
         End Do

         deallocate( psi_l )
      End Do

      deallocate( afft ); deallocate( bfft ); deallocate( cfft )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_soft_part

    subroutine calc_hard_part
      integer :: ik, ib1, ib2, ia, it, ig
      integer :: lmt1, lmt2, it1, it2, il1, il2, n
      integer :: lmta1, lmta2, ilm3, l3, iiqitg, mdvdb
      real(kind=DP) :: fac
      complex(kind=CMPLDP) :: wf1, wf2, z1, z2, zph

      integer, allocatable :: il3(:)
      real(kind=DP), allocatable :: zfsin(:), zfcos(:), qx(:), qy(:), qz(:), ylm(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri(:)
      complex(kind=CMPLDP), allocatable :: wk_fsri_add(:)

      integer :: j, itmp

! -- init
      allocate( zfcos( nmax_G ) ); zfcos = 0.0d0
      allocate( zfsin( nmax_G ) ); zfsin = 0.0d0

      n = nloc
      n=(n-1)+(n-1)+1

      allocate(il3(n**2));call substitute_il3(n**2,il3)

      allocate( qx(nmax_G) ); allocate( qy(nmax_G) ); allocate( qz(nmax_G) )
      allocate( ylm(nmax_G) )

      Do ig=1, nmax_G
         qx(ig) = rltv(1,1)*ngabc(ig,1) +rltv(1,2)*ngabc(ig,2) +rltv(1,3)*ngabc(ig,3)
         qy(ig) = rltv(2,1)*ngabc(ig,1) +rltv(2,2)*ngabc(ig,2) +rltv(2,3)*ngabc(ig,3)
         qz(ig) = rltv(3,1)*ngabc(ig,1) +rltv(3,2)*ngabc(ig,2) +rltv(3,3)*ngabc(ig,3)
      End do

      allocate( wk_fsri(nlmt_core) ); wk_fsri = 0.0d0

! -- start
      Do ia=1, natm
         it = ityp(ia)

         call calc_phase2(natm,pos,ia,kgp,ngabc,1,nmax_G,zfcos,zfsin)

         Do ik=1, kv3
            if ( map_k(ik) /= myrank_k ) cycle

            Do ib2=1, num_core_states
               wk_fsri = 0.0d0
               wk_fsri(:) = dcmplx( fsr_core_states( ib2,:,ik ), &
                       &            fsi_core_states( ib2,:,ik ) )

               DO ib1=ista_e, iend_e, istep_e

                  Do lmt1=1, ilmt_val(it)
                     if ( lmt1 <= ilmt(it) ) then
                        il1 = ltp(lmt1,it); it1 = taup(lmt1,it)
                        lmta1 = lmta( lmt1,ia )
                        wf1 = dcmplx( fsr_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_l( map_z(ib1), lmta1, ik ) )
                     else
                        itmp = lmt1 -ilmt(it)

                        il1 = ltp_add(itmp,it); it1 = 1
                        lmta1 = lmta_add( itmp,ia )
                        wf1 = dcmplx( fsr_add_l( map_z(ib1), lmta1, ik ), &
                             &        fsi_add_l( map_z(ib1), lmta1, ik ) )
                     endif

                     Do lmt2=1, nlmt_core
                        il2 = qnum_l_to_probe +1; it2 = 1
                        wf2 = wk_fsri(lmt2)

                        fac = 1.0d0

                        do n=1,il2p_XI_vc(lmt1,lmt2)

                           ilm3 = isph_XI_vc(lmt1,lmt2,n); l3=il3(ilm3)
                           iiqitg = iqitg_XI_vc(il1,it1,l3+1)

                           if(iiqitg == 0) cycle

                           call sphr( nmax_G, ilm3, qx, qy, qz, ylm )

!#if 0
!                           z1 = conjg(wf1) *wf2
!                           z2 = wf1 *conjg(wf2)
!                           z1 = ( z1+z2) /2.0d0 *zi**l3 *dl2p(lmt1,lmt2,n,it)
!#else
                           z1 = conjg(wf1) *wf2 *zi**l3 *dl2p_XI_vc(lmt1,lmt2,n)
!#endif
                           if ( imple_method == 1 ) then
                              Do ig=1, nmax_G
                                 z2 = qitg_XI_vc( ig,iiqitg ) *ylm(ig)
                                 zph = dcmplx( zfcos(ig), -zfsin(ig) )

                                 RhoTilde_vc( ig, map_z(ib1), ib2, ik ) &
                                      & = RhoTilde_vc( ig, map_z(ib1), ib2, ik ) &
                                      &   + fac *z1 *z2 *zph
                              End do
                           else
                              Do ig=1, nmax_G
                                 z2 = qitg_XI_vc( ig,iiqitg ) *ylm(ig)
                                 zph = dcmplx( zfcos(ig), zfsin(ig) )

                                 RhoTilde_vc( ig, map_z(ib1), ib2, ik ) &
                                      & = RhoTilde_vc( ig, map_z(ib1), ib2, ik ) &
                                      &   + fac *z1 *z2 *zph
                              End do
                           endif

                        End do
                     End do
                  End do
               End do
            End Do
         End DO
      End Do

      deallocate( il3 );  deallocate( ylm )
      deallocate( qx ); deallocate( qy ); deallocate( qz )
      deallocate( zfcos ); deallocate( zfsin )

      deallocate( wk_fsri )

    end subroutine calc_hard_part

  end subroutine m_XI_set_RhoTilde_vc

! --------------------------------------------------------
!    Calculation of transtion moment/matrix (core-core)
! --------------------------------------------------------
  subroutine m_XI_calc_transition_moment_cc
!
    if ( allocated( trm_cc ) ) deallocate( trm_cc )
    allocate( trm_cc( num_core_states, num_core_states, ista_k:iend_k, 3 ) );
    trm_cc = 0.0d0

    call calc_soft_part
    call multiply_factor

  contains

    subroutine multiply_factor
      integer :: ik, ib1, ib2, is
      real(kind=DP) :: ediff
      real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)

      Do ik=1, kv3, ndim_spinor
         if ( map_k(ik) /= myrank_k ) cycle

         Do ib2=1, num_core_states
            Do ib1=1, num_core_states

               ediff = ene_core_states(ib2) -ene_core_states(ib1)

               if ( ib2 == ib1 ) then
                  ediff = delta_omega
               else if ( abs(ediff) < delta_omega ) then
                  if ( ediff >= 0.0 ) then
                     ediff = delta_omega
                  else
                     ediff = -delta_omega
                  endif
               endif

               Do is=1, ndim_spinor
                  trm_cc( ib1,ib2,ik+is-1,: ) &
                       &   = trm_cc( ib1,ib2,ik+is-1,: ) /ediff
               End do

            End Do
         End Do
      End Do

    end subroutine multiply_factor

    subroutine calc_soft_part
      integer :: ik, ib1, ib2, ig
      real(kind=DP), allocatable :: qx(:), qy(:), qz(:)
      real(kind=DP), allocatable :: wk_zaj(:,:)
      real(kind=DP) :: c1, c2
      complex(kind=CMPLDP) :: z1, z2, zsum(3)

      allocate( qx(kg1) );  allocate( qy(kg1) ); allocate( qz(kg1) );
      qx = 0.0d0;  qy = 0.0d0;  qz = 0.0d0

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, vkxyz, ngabc, &
              &                   rltv, qx, qy, qz )

         Do ib2=1, num_core_states

            Do ib1=1, num_core_states
               zsum = 0.0d0
               if ( kimg == 1 ) then
                  Do ig=1, kg1
                     c1 = psig_core_states(ig,ib1,ik,1)
                     c2 = c1 *psig_core_states(ig,ib2,ik,1)

                     zsum(1) = zsum(1) +c2 *qx(ig)
                     zsum(2) = zsum(2) +c2 *qy(ig)
                     zsum(3) = zsum(3) +c2 *qz(ig)
                  End do
               else
                  Do ig=1, kg1
                     z1 = dcmplx( psig_core_states(ig,ib1,ik,1), &
                          &       psig_core_states(ig,ib1,ik,2) )
                     z2 = conjg(z1) *dcmplx( psig_core_states(ig,ib2,ik,1), &
                          &                  psig_core_states(ig,ib2,ik,2) )

                     zsum(1) = zsum(1) +z2 *qx(ig)
                     zsum(2) = zsum(2) +z2 *qy(ig)
                     zsum(3) = zsum(3) +z2 *qz(ig)
                  End do
               endif
               trm_cc( ib1, ib2, ik, 1:3 ) = zsum(1:3)
            End do
         End do
      End Do
      trm_cc = trm_cc *zi

      deallocate( qx ); deallocate( qy ); deallocate( qz )

    end subroutine calc_soft_part

  end subroutine m_XI_calc_transition_moment_cc

  subroutine m_XI_set_RhoTilde_cc

    allocate( RhoTilde_cc( nmax_G, num_core_states, num_core_states, ista_k:iend_k ) )
    RhoTilde_cc = 0.0d0

    call calc_soft_part

  contains
    subroutine calc_soft_part
      integer :: ik, ib1, ib2, i, i1, ngrid
      real(kind=DP), allocatable :: afft(:), bfft(:), cfft(:)
      real(kind=DP), allocatable :: psi_l(:,:,:,:)

      allocate( afft(nfft) ); afft = 0.0d0
      allocate( bfft(nfft) ); bfft = 0.0d0
      allocate( cfft(nfft) ); cfft = 0.0d0

      call m_FFT_alloc_WF_work()
! --
      ngrid = product(fft_box_size_WF(1:3,1))

      Do ik=1, kv3
         if ( map_k(ik) /= myrank_k ) cycle

         allocate( psi_l( kg1, 1, ik:ik, kimg ) ); psi_l = 0.0d0

         Do ib2=1, num_core_states
            bfft = 0.0d0
            psi_l(1:iba(ik),1,ik,:) = psig_core_states(1:iba(ik),ib2,ik,:)
            call m_ES_WF_in_Rspace_kt( ik, ik, ik, psi_l, bfft )

            Do ib1=1, num_core_states
               psi_l(1:iba(ik),1,ik,:) = psig_core_states(1:iba(ik),ib1,ik,:)
               call m_ES_WF_in_Rspace_kt( ik, ik, ik, psi_l, afft )

               cfft = 0.0d0
               if ( imple_method == 1 ) then
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) =-afft(i)*bfft(i+1) +afft(i+1)*bfft(i)
                  End Do
               else
                  Do i=1, nfft, 2
                     cfft(i)   = afft(i)*bfft(i)   +afft(i+1)*bfft(i+1)
                     cfft(i+1) = afft(i)*bfft(i+1) -afft(i+1)*bfft(i)
                  End do
               endif

               call m_FFT_WF( ELECTRON, nfout, cfft, DIRECT, OFF )    ! R-->G
               cfft = cfft /dble(ngrid)

               if ( kimg == 1 ) then
                  Do i=1, nmax_G
                     i1 = igf( nbase(i,ik) )
                     RhoTilde_cc( i, ib1, ib2, ik ) = cfft(i1)
                  End do
               else
                  if ( imple_method == 1 ) then
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_cc( i, ib1, ib2, ik ) &
                             &          = dcmplx( cfft(i1), -cfft(i1+1) )  ! gpaw
                     End do
                  else
                     Do i=1, nmax_G
                        i1 =2*igf( nbase(i,ik) ) -1
                        RhoTilde_cc( i, ib1, ib2, ik ) &
                             &          = dcmplx( cfft(i1), cfft(i1+1) )
                     End do
                  endif

               endif
            End Do
         End Do

         deallocate( psi_l )
      End Do

      deallocate( afft ); deallocate( bfft ); deallocate( cfft )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_soft_part

  end subroutine m_XI_set_RhoTilde_cc

! --------------------------------------------------------
!    Finalization
! --------------------------------------------------------

  subroutine m_XI_dealloc_arrays
    if ( allocated(ilmt_val) ) deallocate( ilmt_val )
    if ( allocated(ngpt_XI) ) deallocate( ngpt_XI )

! val-val
    if ( allocated( iqitg_XI_vv ) ) deallocate( iqitg_XI_vv )
    if ( allocated( nqitg_XI_vv ) ) deallocate( nqitg_XI_vv )
    if ( allocated(  qitg_XI_vv ) ) deallocate(  qitg_XI_vv )
    if ( allocated(  isph_XI_vv ) ) deallocate(  isph_XI_vv )
    if ( allocated(  il2p_XI_vv ) ) deallocate(  il2p_XI_vv )
    if ( allocated(  dl2p_XI_vv ) ) deallocate(  dl2p_XI_vv )
!
    if ( allocated( Mat_dipole_corr_vv ) ) deallocate( Mat_dipole_corr_vv )

    if ( allocated( trm_vv ) )  deallocate( trm_vv )
    if ( allocated( RhoTilde_vv ) ) deallocate( RhoTilde_vv )
    if ( allocated( SpectrFn_vv ) ) deallocate( SpectrFn_vv )
    if ( allocated( SpectrTensor_vv ) ) deallocate( SpectrTensor_vv )

! val-core
    if ( allocated( iqitg_XI_vc ) ) deallocate( iqitg_XI_vc )
    if ( allocated(  qitg_XI_vc ) ) deallocate(  qitg_XI_vc )
    if ( allocated(  isph_XI_vc ) ) deallocate(  isph_XI_vc )
    if ( allocated(  il2p_XI_vc ) ) deallocate(  il2p_XI_vc )
    if ( allocated(  dl2p_XI_vc ) ) deallocate(  dl2p_XI_vc )
!
    if ( allocated( Mat_dipole_corr_vc ) ) deallocate( Mat_dipole_corr_vc )

    if ( allocated( trm_vc ) )  deallocate( trm_vc )
    if ( allocated( RhoTilde_vc ) ) deallocate( RhoTilde_vc )
    if ( allocated( SpectrFn_vc ) ) deallocate( SpectrFn_vc )

! core-core
    if ( allocated( trm_cc ) )  deallocate( trm_cc )
    if ( allocated( RhoTilde_cc ) ) deallocate( RhoTilde_cc )

! Kernel
    if ( allocated( Kernel_Coulomb ) ) deallocate( Kernel_Coulomb )

  end subroutine m_XI_dealloc_arrays

! --------------------------------------------------------
!    Spectrum
! --------------------------------------------------------

  subroutine m_XI_init_spectrum
    real(kind=DP) :: emin

    if ( sw_bse == ON ) then
       allocate( epsinv_omega0( nmax_G, nmax_G ) ) ; epsinv_omega0 = 0.0d0
    endif
    call m_Files_open_nf_xi_spectra()

    if ( sw_corelevel_spectrum == ON ) then
       emin = efermi - maxval( ene_core_states ) +scissor
       e_low  = e_low +emin
       e_high = e_high +emin
    endif

#ifdef USE_ASMS_EXCITATION
    call ASMS_XI_setup( nspin, kv3_fbz, univol, nmax_G, qvec, &
         &              way_BZ_integral, width, gaussian_extent, &
         &              nstep, e_low, e_high, e_step, &
         &              nfout, nf_excitation_spectra, &
         &              npes, mype, MPI_CommGroup, &
         &              sw_nlf, kernel_xc_type, alpha_LRC, spectrum_type, &
         &              sw_spin_decomposition )
#endif

  end subroutine m_XI_init_spectrum

! --------------------------------------------------------
!    Calculation of SpectrFn (val-val)
! --------------------------------------------------------

  subroutine m_XI_calc_spectr_fn_vv_spindec
    integer :: i, ik, jk, iopr, trev
    integer :: ig1, ig2, itmp, ib1, ib2, is1, is2, ispin
    complex(kind=CMPLDP), allocatable :: work(:,:,:)
    real(kind=DP) :: qvec_rotated(3), weight

    real(kind=DP) :: ediff, odiff, occ1, occ2, qv2, txyz(3), fp
    complex(kind=CMPLDP) :: zrho1, zrho2, ztmp, z1, zph

    real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)
    real(kind=DP), allocatable :: occ_mpi(:), occ_wk(:)
    complex(kind=CMPLDP), allocatable :: ztrm2(:,:)

    allocate( SpectrFn_vv( nmax_G_spin, nmax_G_spin, nstep ) );  SpectrFn_vv = 0.0d0
    allocate( ztrm2( nmax_G_spin, nmax_G_spin ) ); ztrm2 = 0.0d0

    allocate( eko_wk( neg ) ); eko_wk = 0.0d0
    allocate( occ_wk( neg ) ); occ_wk = 0.0d0

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k ) cycle

       eko_wk = 0.0d0
       Do ib1=ista_e, iend_e, istep_e
          eko_wk(ib1) = eko_l( map_z(ib1),ik )
       End do
       if ( npes > 1 ) then
          allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
          call mpi_allreduce( eko_wk, eko_mpi, neg, mpi_double_precision, mpi_sum, &
               &              mpi_k_world(myrank_k), ierr )
          eko_wk = eko_mpi
          deallocate( eko_mpi )
       endif

       occ_wk = 0.0d0
       Do ib1=ista_e, iend_e, istep_e
          occ_wk(ib1) = occup_l( map_z(ib1),ik )
       End do
       if ( npes > 1 ) then
          allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
          call mpi_allreduce( occ_wk, occ_mpi, neg, mpi_double_precision, mpi_sum, &
               &              mpi_k_world(myrank_k), ierr )
          occ_wk = occ_mpi
          deallocate( occ_mpi )
       endif

       weight = kv3 *qwgt(ik) /dble(ndim_spinor)

       ispin = 1
       if ( (.not. noncol) .and. nspin == 2 ) then
          ispin = mod( ik-1, 2 )+1
       endif

       Do i=1, num_star_of_k(ik)
          jk = star_of_k(ik,i)
          iopr = iopr_k_fbz_to_ibz(jk)
          trev = trev_k_fbz_to_ibz(jk)

          qvec_rotated = matmul( op(:,:,iopr), qvec(:) )
          if ( trev == 1 ) qvec_rotated = -qvec_rotated

          txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

          Do ib1=ista_e, iend_e, istep_e
             if ( eko_l( map_z(ib1),ik ) > efermi ) cycle

             Do ib2=1, neg
                if ( eko_wk(ib2) < efermi ) cycle

                ediff = eko_wk(ib2)   -eko_l(map_z(ib1),ik) +scissor
                occ1 = occup_l(map_z(ib1),ik) / weight
                occ2 = occ_wk(ib2) / weight

                if ( imple_method == 1 ) then
                   odiff = occ1 -occ2
                else
                   odiff = occ1 *( 1.0d0 -occ2 )
                endif

                if ( ediff < 0.0d0 ) cycle

                Do is1=1, ndim_spinor
                   Do is2=1, ndim_spinor

                      Do ig1=1, nmax_G
                         if ( ig1 == 1 ) then
                            zrho1 = trm_vv( map_z(ib1),ib2,ik+is1-1,1 )*qvec_rotated(1) &
                                 & +trm_vv( map_z(ib1),ib2,ik+is1-1,2 )*qvec_rotated(2) &
                                 &+ trm_vv( map_z(ib1),ib2,ik+is1-1,3 )*qvec_rotated(3)
                            if ( imple_method == 1 ) then
                               zrho1 = -zrho1 *zi
                            else
                               zrho1 = zrho1 *zi
                            endif

                         else
                            itmp = ngpt_XI( ig1, iopr )
                            zrho1 = 0.0d0

                            fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                                 &                     + ngabc(itmp,3)*txyz(3)
                            if ( imple_method == 1 ) then
                               zph = dcmplx( cos(fp), -sin(fp) )
                            else
                               zph = dcmplx( cos(fp), sin(fp) )
                            endif

                            z1 = RhoTilde_vv( itmp, map_z(ib1), ib2, ik+is1-1 )*zph
                            if ( trev == 1 ) then
                               zrho1 = zrho1 + conjg(z1)
                            else
                               zrho1 = zrho1 + z1
                            endif
                         end if

                         Do ig2=1, nmax_G
                            if ( ig2 == 1 ) then
                               zrho2 = trm_vv( map_z(ib1),ib2,ik+is2-1,1 ) &
                                    &                *qvec_rotated(1) &
                                    & +trm_vv( map_z(ib1),ib2,ik+is2-1,2 ) &
                                    &                *qvec_rotated(2) &
                                    & + trm_vv( map_z(ib1),ib2,ik+is2-1,3 ) &
                                    &                *qvec_rotated(3)
                               if ( imple_method == 1 ) then
                                  zrho2 = -zrho2 *zi
                               else
                                  zrho2 = zrho2 *zi
                               endif

                            else
                               itmp = ngpt_XI( ig2, iopr )
                               zrho2 = 0.0d0

                               fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                                    &                     + ngabc(itmp,3)*txyz(3)
                               if ( imple_method == 1 ) then
                                  zph = dcmplx( cos(fp), -sin(fp) )
                               else
                                  zph = dcmplx( cos(fp), sin(fp) )
                               endif

                               z1 = RhoTilde_vv( itmp, map_z(ib1), ib2, ik+is2-1 )*zph
                               if ( trev == 1 ) then
                                  zrho2 = zrho2 + conjg(z1)
                               else
                                  zrho2 = zrho2 + z1
                               endif
                            endif
!
                            if ( noncol ) then
                               ztrm2( nmax_G*(is1-1)+ig1, nmax_G*(is2-1)+ig2 ) &
                                    &           = zrho1 *conjg( zrho2 )
                            else
                               ztrm2( nmax_G*(ispin-1)+ig1, nmax_G*(ispin-1)+ig2 ) &
                                    &           = zrho1 *conjg( zrho2 )
                            endif

                         End do  ! ig2
                      End Do  ! ig1
                   End do
                End Do
#ifdef USE_ASMS_EXCITATION
                call ASMS_XI_add_spectral_function( ediff, odiff, &
                     &                              ztrm2, SpectrFn_vv )
#endif
             End do
          End Do
       End do
    End Do

    if ( npes > 1 ) then
       allocate( work( nmax_G_spin, nmax_G_spin, nstep ) ); work = 0.0d0
       call mpi_allreduce( SpectrFn_vv, work, nmax_G_spin**2 *2 *nstep, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       SpectrFn_vv = work
       deallocate( work )
    endif

    SpectrFn_vv = SpectrFn_vv /dble(kv3_fbz) /univol
    if ( nspin == 1 ) SpectrFn_vv = SpectrFn_vv *2.0d0

#if 0
    Do i=1, nstep
       ztmp = 0.0d0
       Do is1=1, nspin_m
          Do is2=1, nspin_m
             ztmp = ztmp +SpectrFn_vv( nmax_G*(is1-1)+1, nmax_G*(is2-1)+1, i )
          End do
       End Do
       write(2550+mype,*) i*e_step*Hartree, real(ztmp), aimag(ztmp)
    End Do
#endif

  end subroutine m_XI_calc_spectr_fn_vv_spindec

  subroutine m_XI_calc_spectr_fn_vv
    integer :: i, ik, jk, iopr, trev
    integer :: ig1, ig2, itmp, ib1, ib2, is
    complex(kind=CMPLDP), allocatable :: work(:,:,:)
    real(kind=DP) :: qvec_rotated(3), weight

    real(kind=DP) :: ediff, odiff, occ1, occ2, qv2, txyz(3), fp
    complex(kind=CMPLDP) :: zrho1, zrho2, ztmp, z1, zph

    real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)
    real(kind=DP), allocatable :: occ_mpi(:), occ_wk(:)
    complex(kind=CMPLDP), allocatable :: ztrm2(:,:), zrho_work(:)

    allocate( SpectrFn_vv( nmax_G, nmax_G, nstep ) );  SpectrFn_vv = 0.0d0
    allocate( ztrm2(nmax_G,nmax_G) ); ztrm2 = 0.0d0

    allocate( eko_wk( neg ) ); eko_wk = 0.0d0
    allocate( occ_wk( neg ) ); occ_wk = 0.0d0
    allocate( zrho_work( nmax_G ) ); zrho_work = 0.0d0

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k ) cycle

       eko_wk = 0.0d0
       Do ib1=ista_e, iend_e, istep_e
          eko_wk(ib1) = eko_l( map_z(ib1),ik )
       End do
       if ( npes > 1 ) then
          allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
          call mpi_allreduce( eko_wk, eko_mpi, neg, mpi_double_precision, mpi_sum, &
               &              mpi_k_world(myrank_k), ierr )
          eko_wk = eko_mpi
          deallocate( eko_mpi )
       endif

       occ_wk = 0.0d0
       Do ib1=ista_e, iend_e, istep_e
          occ_wk(ib1) = occup_l( map_z(ib1),ik )
       End do
       if ( npes > 1 ) then
          allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
          call mpi_allreduce( occ_wk, occ_mpi, neg, mpi_double_precision, mpi_sum, &
               &              mpi_k_world(myrank_k), ierr )
          occ_wk = occ_mpi
          deallocate( occ_mpi )
       endif

       weight = kv3 *qwgt(ik) /dble(ndim_spinor)

       Do i=1, num_star_of_k(ik)
          jk = star_of_k(ik,i)
          iopr = iopr_k_fbz_to_ibz(jk)
          trev = trev_k_fbz_to_ibz(jk)

          qvec_rotated = matmul( op(:,:,iopr), qvec(:) )
          if ( trev == 1 ) qvec_rotated = -qvec_rotated

          txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

          Do ib1=ista_e, iend_e, istep_e
             if ( eko_l( map_z(ib1),ik ) > efermi ) cycle

             Do ib2=1, neg
                if ( eko_wk(ib2) < efermi ) cycle

                ediff = eko_wk(ib2)   -eko_l(map_z(ib1),ik) +scissor
                occ1 = occup_l(map_z(ib1),ik) / weight
                occ2 = occ_wk(ib2) / weight

                if ( imple_method == 1 ) then
                   odiff = occ1 -occ2
                else
                   odiff = occ1 *( 1.0d0 -occ2 )
                endif

                if ( ediff < 0.0d0 ) cycle
! -----
                Do ig1=1, nmax_G
                   zrho1 = 0.0d0
                   if ( ig1 == 1 ) then
                      Do is=1, ndim_spinor
                         zrho1 = zrho1 &
                              & + trm_vv( map_z(ib1),ib2,ik+is-1,1 )*qvec_rotated(1) &
                              & + trm_vv( map_z(ib1),ib2,ik+is-1,2 )*qvec_rotated(2) &
                              & + trm_vv( map_z(ib1),ib2,ik+is-1,3 )*qvec_rotated(3)
                      End do
                      if ( imple_method == 1 ) then
                         zrho1 = -zrho1 *zi
                      else
                         zrho1 = zrho1 *zi
                      endif
                   else
                      itmp = ngpt_XI( ig1, iopr )
                      fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                           &                     + ngabc(itmp,3)*txyz(3)
                      if ( imple_method == 1 ) then
                         zph = dcmplx( cos(fp), -sin(fp) )
                      else
                         zph = dcmplx( cos(fp), sin(fp) )
                      endif

                      Do is=1, ndim_spinor
                         z1 = RhoTilde_vv( itmp, map_z(ib1), ib2, ik+is-1 )*zph
                         if ( trev == 1 ) then
                            zrho1 = zrho1 + conjg(z1)
                         else
                            zrho1 = zrho1 + z1
                         endif
                      End do
                   end if
                   zrho_work(ig1) = zrho1
                End Do
! -----
                Do ig1=1, nmax_G
                   zrho1 = zrho_work(ig1)
                   Do ig2=1, nmax_G
                      zrho2 = zrho_work(ig2)
                      ztrm2( ig1, ig2 ) = zrho1 *conjg( zrho2 )
                   End do  ! ig2
                End Do  ! ig1

#ifdef USE_ASMS_EXCITATION
                call ASMS_XI_add_spectral_function( ediff, odiff, ztrm2, &
                     &                              SpectrFn_vv )
#endif
             End do
          End Do
       End do
    End Do

    if ( npes > 1 ) then
       allocate( work(nmax_G, nmax_G, nstep ) ); work = 0.0d0
       call mpi_allreduce( SpectrFn_vv, work, nmax_G**2 *2 *nstep, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       SpectrFn_vv = work
       deallocate( work )
    endif

    deallocate( zrho_work ); deallocate( ztrm2 )
    deallocate( eko_wk ); deallocate( occ_wk )

    SpectrFn_vv = SpectrFn_vv /dble(kv3_fbz/nspin) /univol
    if ( nspin == 1 ) SpectrFn_vv = SpectrFn_vv *2.0d0

#if 0
    Do i=1, nstep
       ztmp = SpectrFn_vv(1,1,i)
       write(2450+mype,*) i*e_step*Hartree, real(ztmp), aimag(ztmp)
    End Do
#endif

  end subroutine m_XI_calc_spectr_fn_vv

  subroutine m_XI_calc_spectr_tensor_vv
    integer :: i, ik, jk, iopr, trev
    integer :: ig1, ig2, itmp, ib1, ib2, is, itns

    real(kind=DP) :: qv, qvec1(3), qvec2(3)
    real(kind=DP) :: qvec1_rotated(3), qvec2_rotated(3), weight
    real(kind=DP) :: ediff, odiff, occ1, occ2, txyz(3), fp
    complex(kind=CMPLDP) :: zrho1, zrho2, ztmp, z1, zph

    real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)
    real(kind=DP), allocatable :: occ_mpi(:), occ_wk(:)
    complex(kind=CMPLDP), allocatable :: work(:,:,:)
    complex(kind=CMPLDP), allocatable :: ztrm2(:,:)

    if ( allocated( SpectrTensor_vv ) ) deallocate( SpectrTensor_vv )
    allocate( SpectrTensor_vv( nmax_G, nmax_G, nstep, 6 ) );  SpectrTensor_vv = 0.0d0

    allocate( ztrm2(nmax_G,nmax_G) ); ztrm2 = 0.0d0
    allocate( eko_wk( neg ) ); eko_wk = 0.0d0
    allocate( occ_wk( neg ) ); occ_wk = 0.0d0
!
    qv = sqrt( qvec(1)**2 +qvec(2)**2 + qvec(3)**2 )

    Do itns=1, 6
       qvec1 = 0.d0;   qvec2 = 0.0d0
       select case (itns)
       case (1)
          qvec1(1) = qv;  qvec2(1) = qv
       case (2)
          qvec1(2) = qv;  qvec2(2) = qv
       case (3)
          qvec1(3) = qv;  qvec2(3) = qv
       case (4)
          qvec1(1) = qv;  qvec2(2) = qv
       case (5)
          qvec1(1) = qv;  qvec2(3) = qv
       case (6)
          qvec1(2) = qv;  qvec2(3) = qv
       end select

       Do ik=1, kv3, ndim_spinor
          if ( map_k(ik) /= myrank_k ) cycle

          eko_wk = 0.0d0
          Do ib1=ista_e, iend_e, istep_e
             eko_wk(ib1) = eko_l( map_z(ib1),ik )
          End do
          if ( npes > 1 ) then
             allocate( eko_mpi( neg ) ); eko_mpi = 0.0d0
             call mpi_allreduce( eko_wk, eko_mpi, neg, mpi_double_precision, mpi_sum, &
                  &              mpi_k_world(myrank_k), ierr )
             eko_wk = eko_mpi
             deallocate( eko_mpi )
          endif

          occ_wk = 0.0d0
          Do ib1=ista_e, iend_e, istep_e
             occ_wk(ib1) = occup_l( map_z(ib1),ik )
          End do
          if ( npes > 1 ) then
             allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
             call mpi_allreduce( occ_wk, occ_mpi, neg, mpi_double_precision, mpi_sum, &
                  &              mpi_k_world(myrank_k), ierr )
             occ_wk = occ_mpi
             deallocate( occ_mpi )
          endif

          weight = kv3 *qwgt(ik) /dble(ndim_spinor)

          Do i=1, num_star_of_k(ik)
             jk = star_of_k(ik,i)
             iopr = iopr_k_fbz_to_ibz(jk)
             trev = trev_k_fbz_to_ibz(jk)

             qvec1_rotated = matmul( op(:,:,iopr), qvec1(:) )
             qvec2_rotated = matmul( op(:,:,iopr), qvec2(:) )

             if ( trev == 1 ) qvec1_rotated = -qvec1_rotated
             if ( trev == 1 ) qvec2_rotated = -qvec2_rotated

             txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

             Do ib1=ista_e, iend_e, istep_e
                if ( eko_l( map_z(ib1),ik ) > efermi ) cycle

                Do ib2=1, neg
                   if ( eko_wk(ib2) < efermi ) cycle

                   ediff = eko_wk(ib2)   -eko_l(map_z(ib1),ik) +scissor
                   occ1 = occup_l(map_z(ib1),ik) / weight
                   occ2 = occ_wk(ib2) / weight

                   if ( imple_method == 1 ) then
                      odiff = occ1 -occ2
                   else
                      odiff = occ1 *( 1.0d0 -occ2 )
                   endif

                   if ( ediff < 0.0d0 ) cycle

                   Do ig1=1, nmax_G
                      if ( ig1 == 1 ) then
                         zrho1 = 0.0d0
                         Do is=1, ndim_spinor
                            zrho1 = zrho1 &
                                 & + trm_vv(map_z(ib1),ib2,ik+is-1,1)*qvec1_rotated(1) &
                                 & + trm_vv(map_z(ib1),ib2,ik+is-1,2)*qvec1_rotated(2) &
                                 & + trm_vv(map_z(ib1),ib2,ik+is-1,3)*qvec1_rotated(3)
                         End do
                         if ( imple_method == 1 ) then
                            zrho1 = -zrho1 *zi
                         else
                            zrho1 = zrho1 *zi
                         endif

                      else
                         itmp = ngpt_XI( ig1, iopr )
                         zrho1 = 0.0d0

                         fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                              &                     + ngabc(itmp,3)*txyz(3)
                         if ( imple_method == 1 ) then
                            zph = dcmplx( cos(fp), -sin(fp) )
                         else
                            zph = dcmplx( cos(fp), sin(fp) )
                         endif

                         Do is=1, ndim_spinor
                            z1 = RhoTilde_vv( itmp, map_z(ib1), ib2, ik+is-1 )*zph
                            if ( trev == 1 ) then
                               zrho1 = zrho1 + conjg(z1)
                            else
                               zrho1 = zrho1 + z1
                            endif
                         End do

                      endif

                      Do ig2=1, nmax_G
                         if ( ig2 == 1 ) then
                            zrho2 = 0.0d0
                            Do is=1, ndim_spinor
                               zrho2 = zrho2 &
                                    & + trm_vv(map_z(ib1),ib2,ik+is-1,1) &
                                    &                   *qvec2_rotated(1) &
                                    & + trm_vv(map_z(ib1),ib2,ik+is-1,2) &
                                    &                   *qvec2_rotated(2) &
                                    & + trm_vv(map_z(ib1),ib2,ik+is-1,3) &
                                    &                   *qvec2_rotated(3)
                            End do
                            if ( imple_method == 1 ) then
                               zrho2 = -zrho2 *zi
                            else
                               zrho2 = zrho2 *zi
                            endif

                         else
                            itmp = ngpt_XI( ig2, iopr )
                            zrho2 = 0.0d0

                            fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                                 &                     + ngabc(itmp,3)*txyz(3)
                            if ( imple_method == 1 ) then
                               zph = dcmplx( cos(fp), -sin(fp) )
                            else
                               zph = dcmplx( cos(fp), sin(fp) )
                            endif

                            Do is=1, ndim_spinor
                               z1 = RhoTilde_vv( itmp, map_z(ib1), ib2, ik+is-1 )*zph
                               if ( trev == 1 ) then
                                  zrho2 = zrho2 + conjg(z1)
                               else
                                  zrho2 = zrho2 + z1
                               endif
                            End do

                         endif

                         ztrm2( ig1, ig2 ) = zrho1 *conjg( zrho2 )

                      End do  ! ig2
                   End Do  ! ig1

#ifdef USE_ASMS_EXCITATION
                   call ASMS_XI_add_spectral_function( ediff, odiff, ztrm2, &
                        &                              SpectrTensor_vv(:,:,:,itns) )
#endif
                End do
             End Do
          End do
       End Do

       if ( npes > 1 ) then
          allocate( work(nmax_G, nmax_G, nstep ) ); work = 0.0d0
          call mpi_allreduce( SpectrTensor_vv(:,:,:,itns), work, nmax_G**2 *2 *nstep, &
               &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
          SpectrTensor_vv(:,:,:,itns) = work(:,:,:)
          deallocate( work )
       endif
    End Do

    SpectrTensor_vv = SpectrTensor_vv /dble(kv3_fbz/nspin) /univol
    if ( nspin == 1 ) SpectrTensor_vv = SpectrTensor_vv *2.0d0

#if 0
    Do i=1, nstep
       ztmp = SpectrTensor_vv(1,1,i,1)
       write(2450+mype,*) i*e_step*Hartree, real(ztmp), aimag(ztmp)
    End Do
#endif

    deallocate( ztrm2 ); deallocate( eko_wk ); deallocate( occ_wk )

  end subroutine m_XI_calc_spectr_tensor_vv

! --------------------------------------------------------
!    Calculation of SpectrFn (val-core)
! --------------------------------------------------------

  subroutine m_XI_calc_spectr_fn_vc
    integer :: i, ik, jk, iopr, trev
    integer :: ig1, ig2, itmp, ib1, ib2, is
    complex(kind=CMPLDP), allocatable :: work(:,:,:)
    real(kind=DP) :: qvec_rotated(3), weight

    real(kind=DP) :: ediff, odiff, occ1, occ2, qv2, txyz(3), fp
    complex(kind=CMPLDP) :: zrho1, zrho2, ztmp, z1, zph
    real(kind=DP) :: Delta03 = 1.0D-3

    real(kind=DP), allocatable :: eko_mpi(:), eko_wk(:)
    real(kind=DP), allocatable :: occ_mpi(:), occ_wk(:)
    complex(kind=CMPLDP), allocatable :: ztrm2(:,:)

    write(*,*) "nstep = ", nstep

    allocate( SpectrFn_vc( nmax_G, nmax_G, nstep ) );  SpectrFn_vc = 0.0d0
    allocate( ztrm2(nmax_G,nmax_G) ); ztrm2 = 0.0d0

    allocate( eko_wk( num_core_states ) ); eko_wk = 0.0d0

    Do ik=1, kv3, ndim_spinor
       if ( map_k(ik) /= myrank_k ) cycle

       eko_wk = 0.0d0
       eko_wk(:) = ene_core_states(:)

       weight = kv3 *qwgt(ik) /dble(ndim_spinor)

       Do i=1, num_star_of_k(ik)
          jk = star_of_k(ik,i)
          iopr = iopr_k_fbz_to_ibz(jk)
          trev = trev_k_fbz_to_ibz(jk)

          qvec_rotated = matmul( op(:,:,iopr), qvec(:) )
          if ( trev == 1 ) qvec_rotated = -qvec_rotated

          txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

          Do ib1=ista_e, iend_e, istep_e
             occ1 = occup_l(map_z(ib1),ik) / weight
             if ( 1.0D0 -occ1 < Delta03 ) cycle

             Do ib2=1, num_core_states

                ediff = eko_wk(ib2)   -eko_l(map_z(ib1),ik) -scissor
                occ2 = 1.0d0

                if ( imple_method == 1 ) then
                   odiff = occ1 -occ2
                else
                   odiff = occ1 *( 1.0d0 -occ2 )
                endif
                ediff = -ediff;    odiff = -odiff

                Do ig1=1, nmax_G
                   if ( ig1 == 1 ) then
                      zrho1 = 0.0d0
                      Do is=1, ndim_spinor
                         zrho1 = zrho1 &
                              & + trm_vc( map_z(ib1),ib2,ik+is-1,1 )*qvec_rotated(1) &
                              & + trm_vc( map_z(ib1),ib2,ik+is-1,2 )*qvec_rotated(2) &
                              & + trm_vc( map_z(ib1),ib2,ik+is-1,3 )*qvec_rotated(3)
                      End do
                      if ( imple_method == 1 ) then
                         zrho1 = -zrho1 *zi
                      else
                         zrho1 = zrho1 *zi
                      endif

                   else
                      itmp = ngpt_XI( ig1, iopr )
                      zrho1 = 0.0d0
                      fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                           &                     + ngabc(itmp,3)*txyz(3)
                      if ( imple_method == 1 ) then
                         zph = dcmplx( cos(fp), -sin(fp) )
                      else
                         zph = dcmplx( cos(fp), sin(fp) )
                      endif
                      Do is=1, ndim_spinor
                         z1 = RhoTilde_vc( itmp, map_z(ib1), ib2, ik+is-1 ) *zph
                         if ( trev == 1 ) then
                            zrho1 = zrho1 + conjg(z1)
                         else
                            zrho1 = zrho1 + z1
                         endif
                      End do
                   endif

                   Do ig2=1, nmax_G
                      if ( ig2 == 1 ) then
                         zrho2 = 0.0d0
                         Do is=1, ndim_spinor
                            zrho2 = zrho2 &
                                 & + trm_vc( map_z(ib1),ib2,ik+is-1,1 )*qvec_rotated(1) &
                                 & + trm_vc( map_z(ib1),ib2,ik+is-1,2 )*qvec_rotated(2) &
                                 & + trm_vc( map_z(ib1),ib2,ik+is-1,3 )*qvec_rotated(3)
                         End do
                         if ( imple_method == 1 ) then
                            zrho2 = -zrho2 *zi
                         else
                            zrho2 = zrho2 *zi
                         endif

                      else
                         itmp = ngpt_XI( ig2, iopr )
                         zrho2 = 0.0d0

                         fp = ngabc(itmp,1)*txyz(1) + ngabc(itmp,2)*txyz(2) &
                              &                     + ngabc(itmp,3)*txyz(3)
                         if ( imple_method == 1 ) then
                            zph = dcmplx( cos(fp), -sin(fp) )
                         else
                            zph = dcmplx( cos(fp), sin(fp) )
                         endif

                         Do is=1, ndim_spinor
                            z1 = RhoTilde_vc( itmp, map_z(ib1), ib2, ik+is-1 )*zph
                            if ( trev == 1 ) then
                               zrho2 = zrho2 + conjg(z1)
                            else
                               zrho2 = zrho2 + z1
                            endif

                         End do
                      endif

                      ztrm2( ig1, ig2 ) = zrho1 *conjg( zrho2 )

                   End do  ! ig2
                End Do  ! ig1

#ifdef USE_ASMS_EXCITATION
                call ASMS_XI_add_spectral_function( ediff, odiff, ztrm2, &
                     &                              SpectrFn_vc )
#endif
             End do
          End Do
       End do
    End Do

    if ( npes > 1 ) then
       allocate( work(nmax_G, nmax_G, nstep ) ); work = 0.0d0
       call mpi_allreduce( SpectrFn_vc, work, nmax_G**2 *2 *nstep, &
            &              mpi_double_precision, mpi_sum, MPI_CommGroup, ierr )
       SpectrFn_vc = work
       deallocate( work )
    endif

    SpectrFn_vc = SpectrFn_vc /dble(kv3_fbz/nspin) /univol
    if ( nspin == 1 ) SpectrFn_vc = SpectrFn_vc *2.0d0

#if 0
    Do i=1, nstep
       ztmp = SpectrFn_vc(1,1,i)
       write(2400+mype,*) e(i)*Hartree, real(ztmp), aimag(ztmp)
    End do
#endif

  end subroutine m_XI_calc_spectr_fn_vc

! --------------------------------------------------------
! Coulomb Kernel
! --------------------------------------------------------
  subroutine m_XI_set_Kernel_Coulomb
    integer :: i
    real(kind=DP) :: g2, qv2

    allocate( Kernel_Coulomb(nmax_G) ); Kernel_Coulomb = 0.0d0

    Do i=2, nmax_G
       g2 = ttr(1)*ngabc(i,1)*ngabc(i,1) &
            &             + ttr(2)*ngabc(i,2)*ngabc(i,2) &
            &             + ttr(3)*ngabc(i,3)*ngabc(i,3) &
            &             + ttr(4)*ngabc(i,1)*ngabc(i,2) &
            &             + ttr(5)*ngabc(i,2)*ngabc(i,3) &
            &             + ttr(6)*ngabc(i,3)*ngabc(i,1)
       Kernel_Coulomb(i) = PAI4 / g2
    End do
    qv2 = qvec(1)**2 +qvec(2)**2 + qvec(3)**2
    Kernel_Coulomb(1) = PAI4 / qv2

  end subroutine m_XI_set_Kernel_Coulomb

  subroutine m_XI_calc_num_occupied_bands
    integer :: ik, ib1
    real(kind=DP) :: weight, c1
    real(kind=DP), allocatable :: occ_mpi(:), occ_wk(:)

    allocate( occ_wk( neg ) ); occ_wk = 0.0d0

    Do ik=1, kv3
       if ( map_k(ik) /= myrank_k ) cycle

       occ_wk = 0.0d0
       Do ib1=ista_e, iend_e, istep_e
          occ_wk(ib1) = occup_l( map_z(ib1),ik )
       End do
       if ( npes > 1 ) then
          allocate( occ_mpi( neg ) ); occ_mpi = 0.0d0
          call mpi_allreduce( occ_wk, occ_mpi, neg, mpi_double_precision, mpi_sum, &
               &              mpi_k_world(myrank_k), ierr )
          occ_wk = occ_mpi
          deallocate( occ_mpi )
       endif

       weight = kv3 *qwgt(ik) /dble(ndim_spinor)

       num_occ_bands = 0
       Do ib1=1, neg
          c1 = occ_wk(ib1) /weight
          if ( c1 > 0.5 ) num_occ_bands = num_occ_bands +1
       End do
    End Do

    num_unocc_bands = neg -num_occ_bands

    deallocate( occ_wk )

  end subroutine m_XI_calc_num_occupied_bands

  subroutine m_XI_set_matsize_bse
    integer :: i, ista, iend, n1, n2
    integer, allocatable :: matsize_bse_on_pe(:)

    if ( sw_corelevel_spectrum == ON ) then
       matsize_bse = num_core_states *num_unocc_bands *kv3_fbz
    else
       matsize_bse = num_occ_bands *num_unocc_bands *kv3_fbz
    endif
    matsize_bse = matsize_bse /dble(ndim_spinor)

    allocate( matsize_bse_on_pe( npes ) ); matsize_bse_on_pe = 0

    n1 = int( matsize_bse /npes );  n2 = matsize_bse -n1 *npes
    Do i=1, npes
       matsize_bse_on_pe(i) = n1
    End do
    Do i=1, n2
       matsize_bse_on_pe(i) = matsize_bse_on_pe(i) +1
    End do
!
    ista = 0;  iend = 0

    Do i=1, npes
       ista = iend +1
       iend = iend +matsize_bse_on_pe(i)
       if ( mype == i-1 ) then
          ista_matbse = ista;  iend_matbse = iend
       endif
    End Do
    deallocate( matsize_bse_on_pe )

  end subroutine m_XI_set_matsize_bse

  subroutine m_XI_alloc_Hamiltonian_bse
  end subroutine m_XI_alloc_Hamiltonian_bse

  subroutine m_XI_set_Hamiltonian_bse
    integer :: ista_band, jsta_band

    allocate( MatH_bse( ista_matbse:iend_matbse, matsize_bse ) )
    MatH_bse = 0.0d0

    ista_band = num_occ_bands +1
    jsta_band = ista_band

    if ( sw_corelevel_spectrum == ON ) then
       call set_diagonal_part
       call set_exchange_part_core2val
       call set_screened_coulomb_core2val
    else
    endif
  contains

    subroutine set_exchange_part_core2val       ! electron-hole exchange
      integer :: ik, jk, ib1, ib2, jb1, jb2
      integer :: i, j, is, ig
      integer :: ik_fbz, jk_fbz, iopr, jopr, trev_i, trev_j, ipos, jpos
      integer :: itmp1, itmp2, itmp3, jtmp1, jtmp2, jtmp3
      real(kind=DP) :: fp_i, fp_j, txyz_i(3), txyz_j(3)
      complex(kind=CMPLDP) :: zsum, zph_i, zph_j, z1, z2

      complex(kind=CMPLDP), allocatable :: work_i(:,:), work_j(:,:)

      allocate( work_i(nmax_G,num_core_states) ); work_i = 0.0d0
      allocate( work_j(nmax_G,num_core_states) ); work_j = 0.0d0

      Do ik=1, kv3, ndim_spinor
         Do ib1=ista_band, neg
            work_i = 0.0d0
            if ( map_ek(ib1,ik) == mype ) then
               Do ib2=1, num_core_states
                  Do is=1, ndim_spinor
                     work_i(:,ib2) = work_i(:,ib2) &
                          &        + RhoTilde_vc(:,map_z(ib1),ib2,ik+is-1)
                  End Do
               End Do
            endif
            call mpi_bcast( work_i, nmax_G *num_core_states, &
                 &          mpi_double_precision, &
                 &          map_ek(ib1,ik), MPI_CommGroup, ierr )

            Do jk=1, kv3, ndim_spinor
               Do jb1=jsta_band, neg
                  work_j = 0.0d0
                  if ( map_ek(jb1,jk) == mype ) then
                     Do jb2=1, num_core_states
                        Do is=1, ndim_spinor
                           work_j(:,jb2) = work_j(:,jb2) &
                                &        + RhoTilde_vc(:,map_z(jb1),jb2,jk+is-1)
                        End Do
                     End Do
                  endif
                  call mpi_bcast( work_j, nmax_G *num_core_states, &
                       &          mpi_double_precision, &
                       &          map_ek(jb1,jk), MPI_CommGroup, ierr )

                  Do i=1, num_star_of_k(ik)
                     ik_fbz = star_of_k(ik,i)
                     iopr = iopr_k_fbz_to_ibz(ik_fbz)
                     trev_i = trev_k_fbz_to_ibz(ik_fbz)
                     txyz_i(1:3) = tau(1:3,iopr,BUCS) *PAI2

                     Do j=1, num_star_of_k(jk)
                        jk_fbz = star_of_k(jk,j)
                        jopr = iopr_k_fbz_to_ibz(jk_fbz)
                        trev_j = trev_k_fbz_to_ibz(jk_fbz)
                        txyz_j(1:3) = tau(1:3,jopr,BUCS) *PAI2

                        Do ib2=1, num_core_states
                           Do jb2=1, num_core_states
                              itmp1 = ( neg -ista_band +1 )*num_core_states
                              itmp2 = ( ib1 -ista_band +1 )*num_core_states
                              ipos = ( ik_fbz -1 )*itmp1 +itmp2 +ib2

                              jtmp1 = ( neg -jsta_band +1 )*num_core_states
                              jtmp2 = ( jb1 -jsta_band +1 )*num_core_states
                              jpos = ( jk_fbz -1 )*jtmp1 +jtmp2 +jb2
                              !
                              if ( ipos >= ista_matbse .and. ipos <= iend_matbse ) then
                                 zsum = 0.0d0

                                 Do ig=2, nmax_G         ! assuming sw_longwavelimit
                                    itmp3 = ngpt_XI( ig,iopr )
                                    jtmp3 = ngpt_XI( ig,jopr )
                                    fp_i = ngabc(itmp3,1)*txyz_i(1) &
                                         &   + ngabc(itmp3,2)*txyz_i(2) &
                                         &   + ngabc(itmp3,3)*txyz_i(3)
                                    fp_j = ngabc(jtmp3,1)*txyz_j(1) &
                                         &   + ngabc(jtmp3,2)*txyz_j(2) &
                                         &   + ngabc(jtmp3,3)*txyz_j(3)

                                    if ( imple_method == 1 ) then
                                       zph_i = dcmplx( cos(fp_i), -sin(fp_i) )
                                       zph_j = dcmplx( cos(fp_j), -sin(fp_j) )
                                    else
                                       zph_i = dcmplx( cos(fp_i), sin(fp_i) )
                                       zph_j = dcmplx( cos(fp_j), sin(fp_j) )
                                    endif

                                    Do is=1, ndim_spinor
                                       if ( trev_i == 1 ) then
                                          z1 = conjg( work_i(ig,is) *zph_i )
                                       else
                                          z1 = work_i(ig,is) *zph_i
                                       endif
                                       if ( trev_j == 1 ) then
                                          z2 = conjg( work_j(ig,is) *zph_j )
                                       else
                                          z2 = work_j(ig,is) *zph_j
                                       endif
                                       zsum = zsum +z1 *conjg(z2) *kernel_Coulomb(ig)
                                    End Do
                                 End Do
                                 MatH_bse( ipos, jpos ) = MatH_bse( ipos, jpos ) &
                                      &                 + zsum /univol
                              endif
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

      deallocate( work_i ); deallocate( work_j )

    end subroutine set_exchange_part_core2val

    subroutine set_screened_coulomb_core2val
      integer :: ik, ib1, ib2, jb1, jb2
      integer :: i, is, ik_fbz, iopr, trev
      integer :: itmp1, itmp2, jtmp1, jtmp2, ipos, jpos
      integer :: ig1, ig2, itmp3
      real(kind=DP) :: qvec_rotated(3), txyz(3), fp
      complex(kind=CMPLDP) :: zsum, zrho1, zrho2, zph, z1, z2

      complex(kind=CMPLDP), allocatable :: work1_rho_vv(:,:), work1_trm_vv(:,:)
      complex(kind=CMPLDP), allocatable :: work2_rho_cc(:,:), work2_trm_cc(:,:)

      allocate( work1_rho_vv(nmax_G,neg) ); work1_rho_vv = 0.0d0
      allocate( work1_trm_vv(neg,3) );      work1_trm_vv = 0.0d0
      allocate( work2_rho_cc(nmax_G,num_core_states) ); work2_rho_cc = 0.0d0
      allocate( work2_trm_cc(num_core_states,3) );      work2_trm_cc = 0.0d0

      Do ik=1, kv3, ndim_spinor

         Do ib1=ista_band, neg
            work1_rho_vv = 0.0d0;  work1_trm_vv = 0.0d0
            if ( map_ek(ib1,ik) == mype ) then
               Do jb1=1, neg
                  Do is=1, ndim_spinor
                     work1_rho_vv(:,jb1) = work1_rho_vv(:,jb1) &
                          &              + RhoTilde_vv(:,map_z(ib1),jb1,ik+is-1)
                  End Do
               End Do
               Do jb1=1, neg
                  Do is=1, ndim_spinor
                     work1_trm_vv(jb1,1:3) = work1_trm_vv(jb1,1:3) &
                          &                + trm_vv(map_z(ib1),jb1,ik+is-1,1:3)
                  End Do
               End Do
            endif
            call mpi_bcast( work1_rho_vv, nmax_G *neg, mpi_double_precision, &
                 &          map_ek(ib1,ik), MPI_CommGroup, ierr )
            call mpi_bcast( work1_trm_vv, neg*3, mpi_double_precision, &
                 &          map_ek(ib1,ik), MPI_CommGroup, ierr )

            Do ib2=1, num_core_states
               work2_rho_cc = 0.0d0;    work2_trm_cc = 0.0d0
               if ( map_k(ik) == myrank_k ) then
                  Do jb2=1, num_core_states
                     Do is=1, ndim_spinor
                        work2_rho_cc(:,jb2) = work2_rho_cc(:,ib2) &
                             &                 + RhoTilde_cc(:,ib2,jb2,ik+is-1)
                     End Do
                  End Do
                  Do jb2=1, num_core_states
                     Do is=1, ndim_spinor
                        work2_trm_cc(jb2,1:3) = work2_trm_cc(jb2,1:3) &
                             &                + trm_cc(ib2,jb2,ik+is-1,1:3)
                     End Do
                  End Do
               endif
               call mpi_bcast( work2_rho_cc, nmax_G *num_core_states, &
                    &          mpi_double_precision, &
                    &          map_k(ik), mpi_e_world(myrank_e), ierr )
               call mpi_bcast( work2_trm_cc, num_core_states*3, &
                    &          mpi_double_precision, &
                    &          map_k(ik), mpi_e_world(myrank_e), ierr )

               Do i=1, num_star_of_k(ik)
                  ik_fbz = star_of_k(ik,i)
                  iopr = iopr_k_fbz_to_ibz(ik_fbz)
                  trev = trev_k_fbz_to_ibz(ik_fbz)

                  qvec_rotated = matmul( op(:,:,iopr), qvec(:) )
                  if ( trev == 1 ) qvec_rotated = -qvec_rotated

                  txyz(1:3) = tau(1:3,iopr,BUCS)*PAI2

                  Do jb1=ista_band, neg
                     Do jb2=1, num_core_states
                        itmp1 = ( neg -ista_band +1 )*num_core_states
                        itmp2 = ( ib1 -ista_band +1 )*num_core_states
                        ipos = ( ik_fbz -1 )*itmp1 +itmp2 +ib2
                              !
                        jtmp1 = ( neg -jsta_band +1 )*num_core_states
                        jtmp2 = ( jb1 -jsta_band +1 )*num_core_states
                        jpos = ( ik_fbz -1 )*jtmp1 +jtmp2 +jb2
                              !
                        if ( ipos >= ista_matbse .and. ipos <= iend_matbse ) then
                           zsum = 0.0d0

                           Do ig1=1, nmax_G
                              if ( ig1 == 1 ) then
                                 zrho1 = 0.0d0
                                 Do is=1, ndim_spinor
                                    zrho1 = zrho1 &
                                         & + work1_trm_vv(jb1,1)*qvec_rotated(1) &
                                         & + work1_trm_vv(jb1,2)*qvec_rotated(2) &
                                         & + work1_trm_vv(jb2,3)*qvec_rotated(3)
                                 End do
                                 if ( imple_method == 1 ) then
                                    zrho1 = -zrho1 *zi
                                 else
                                    zrho1 = zrho1 *zi
                                 endif
                              else
                                 itmp3 = ngpt_XI( ig1, iopr )
                                 zrho1 = 0.0d0

                                 fp = ngabc(itmp3,1)*txyz(1) + ngabc(itmp3,2)*txyz(2) &
                                      &                      + ngabc(itmp3,3)*txyz(3)
                                 if ( imple_method == 1 ) then
                                    zph = dcmplx( cos(fp), -sin(fp) )
                                 else
                                    zph = dcmplx( cos(fp), sin(fp) )
                                 endif

                                 z1 = work1_rho_vv(itmp3,jb1) *zph
                                 if ( trev == 1 ) then
                                    zrho1 = zrho1 + conjg(z1)
                                 else
                                    zrho1 = zrho1 + z1
                                 endif
                              endif

                              Do ig2=1, nmax_G
                                 if ( ig2 == 1 ) then
                                    zrho2 = 0.0d0
                                    Do is=1, ndim_spinor
                                       zrho2 = zrho2 &
                                            & + work2_trm_cc(jb2,1)*qvec_rotated(1) &
                                            & + work2_trm_cc(jb2,2)*qvec_rotated(2) &
                                            & + work2_trm_cc(jb2,3)*qvec_rotated(3)
                                    End do
                                    if ( imple_method == 1 ) then
                                       zrho2 = -zrho2 *zi
                                    else
                                       zrho2 = zrho2 *zi
                                    endif
                                 else
                                    itmp3 = ngpt_XI( ig2, iopr )
                                    zrho2 = 0.0d0
                                    fp = ngabc(itmp3,1)*txyz(1) +ngabc(itmp3,2)*txyz(2) &
                                         &                      +ngabc(itmp3,3)*txyz(3)
                                    if ( imple_method == 1 ) then
                                       zph = dcmplx( cos(fp), -sin(fp) )
                                    else
                                       zph = dcmplx( cos(fp), sin(fp) )
                                    endif

                                    z1 = work2_rho_cc(itmp3,jb2) *zph
                                    if ( trev == 1 ) then
                                       zrho2 = zrho2 + conjg(z1)
                                    else
                                       zrho2 = zrho2 + z1
                                    endif
                                 endif
!
                                 zsum = zsum +zrho1 *conjg(zrho2) *kernel_Coulomb(ig1) &
                                      &             *epsinv_omega0(ig1,ig2)
                              End do
                           End Do
                           MatH_bse(ipos,jpos) = MatH_bse(ipos,jpos) &
                                &               - zsum /2.0d0 /univol
                        endif
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

      deallocate( work1_trm_vv ); deallocate( work1_rho_vv )
      deallocate( work2_trm_cc ); deallocate( work2_rho_cc )

    end subroutine set_screened_coulomb_core2val

    subroutine set_diagonal_part
      integer :: ik, ib1, ib2
      integer :: i, iopr, ik_fbz, trev
      integer :: itmp1, itmp2, ipos
      real(kind=DP) :: work1_eko_v, work1_eko_c, c1

      Do ik=1, kv3, ndim_spinor
         Do ib1=ista_band, neg
            work1_eko_v = 0.0d0
            if ( map_ek(ib1,ik) == mype ) then
               work1_eko_v = eko_l(map_z(ib1),ik)
            endif
            call mpi_bcast( work1_eko_v, 1, mpi_double_precision, &
                 &          map_ek(ib1,ik), MPI_CommGroup, ierr )

            Do ib2=1, num_core_states
               work1_eko_c = ene_core_states(ib2)

               Do i=1, num_star_of_k(ik)
                  ik_fbz = star_of_k(ik,i)
                  iopr = iopr_k_fbz_to_ibz(ik_fbz)
                  trev = trev_k_fbz_to_ibz(ik_fbz)

                  itmp1 = ( neg -ista_band +1 )*num_core_states
                  itmp2 = ( ib1 -ista_band +1 )*num_core_states
                  ipos = ( ik_fbz -1 )*itmp1 +itmp2 +ib2

                  if ( ipos >= ista_matbse .and. ipos <= iend_matbse ) then
                     c1 = work1_eko_v -work1_eko_c
                     MatH_bse(ipos,ipos) = MatH_bse(ipos,ipos) +c1
                  endif
               End Do
            End Do
         End Do
      End Do

    end subroutine set_diagonal_part

  end subroutine m_XI_set_Hamiltonian_bse

end module m_Excitation

