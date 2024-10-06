module m_CoreLevel_Spectrum

  use m_Control_Parameters,  only : kimg, nspin, printable, ipriinputfile, noncol, &
       &                            ndim_spinor, sw_corelevel_spectrum, &
       &                            sw_local_approx_trans_moment, &
       &                            sw_v2c_xes, ndim_chgpot
  use m_Const_Parameters, only : DP, CMPLDP, PAI2, PAI4, zi, BUCS, FMAXVALLEN, LOWER, &
       &                         OFF, ON, Hartree, CARTS, GAMMA, ELECTRON, DIRECT
  use m_Kpoints,    only : kv3, vkxyz, qwgt, k_symmetry

  use m_Files,  only : nfout, nfpot, nfcore_energy_out, &
       &               nfcore_energy_initial, nfcore_energy_final, &
       &               m_Files_open_ps_files, m_Files_close_ps_files, &
       &               m_Files_open_core_energy_file, m_Files_close_core_energy_file, &
       &               m_Files_open_core_ene_initial, m_Files_close_core_ene_initial, &
       &               m_Files_open_core_ene_final, m_Files_close_core_ene_final, &
       &               m_Files_open_ps_file, m_Files_close_ps_file


  use m_Parallelization,  only : mype, ierr, myrank_k, map_k, npes, ista_k, iend_k, &
       &                         np_e, ista_kngp, iend_kngp, MPI_CommGroup, ierr, &
       &                         mpi_k_world, nrank_e, np_g1k, ista_g1k

  use m_PlaneWaveBasisSet,  only : kg1, nbase, iba, ngabc, kgp, kg, igf
  use m_Ionic_System,  only : ityp, ivan, iatomn, ntyp, pos, natm, cps

  use m_PseudoPotential,  only : nmesh, radr, wos, xh, rmax, ltp, mtp, taup, &
       &                         ilmt, lmta, psirpw, wf_mnrc, radr_paw, mmesh, flg_paw

  use m_Crystal_Structure,  only : univol, rltv
  use m_NonLocal_Potential,  only : new_radr_and_wos
  use m_Electronic_Structure,     only : efermi, eko_l, occup_l, zaj_l, fsr_l, fsi_l, &
       &                                 vlhxc_l

#if 0
  use m_SPinOrbit_RadInt,  only : m_SO_calc_contrib_corelevels, &
       &                          m_SO_calc_core_energy_atoms, &
       &                          m_SO_calc_corelevel_splitting
#else
  use m_SPinOrbit_RadInt,  only : m_SO_calc_corelevel_splitting, &
       &                          m_SO_calc_contrib_corelevels
#endif

  use m_SpinOrbit_Potential, only :  EigenWfns_MatLS_L1,  EigenWfns_MatLS_L2, &
       &                             EigenWfns_MatLS_L3, &
       &                             m_SO_set_MatU_ylm_RC, &
       &                             m_SO_calc_MatLS_orb_s_to_f, &
       &                             m_SO_diagonalize_MatLS
  use m_Total_Energy,       only : etotal

! === KT_add ==== 2014/08/08
  use m_Const_Parameters,   only : ByPawPot, ReadFromPP
  use m_Control_Parameters,  only :  eval_core_level_splitting
  use m_Orbital_QuantumNum, only : num_orb_index_data, &
       &                           qnum_n_orb_index, qnum_l_orb_index, &
       &                           qnum_t_orb_index, qnum_tau_orb_index
! =============== 2014/08/08

  use m_FFT,            only :  nfft, fft_box_size_WF, m_FFT_alloc_WF_work, &
       &                        m_FFT_dealloc_WF_work, m_FFT_WF

! === KT_add === 2015/01/23
  use m_Control_Parameters,  only : SpinOrbit_mode
! ============== 2015/01/12
  use mpi

  implicit none
!  include 'mpif.h'
!
! --------------
!  General
!
  character(len("corelevel_spectrum")),   parameter  &
       &            ::   tag_corelevel_spectrum   = "corelevel_spectrum"
  character(len("sw_corelevel_spectrum")),   parameter  &
       &            ::   tag_sw_corelevel_spectrum   = "sw_corelevel_spectrum"
  character(len("sw_local_approx_trans_moment")),   parameter  &
       &            ::   tag_sw_local_approx_trans_mom = "sw_local_approx_trans_moment"
  character(len("sw_v2c_xes")),   parameter  &
       &            ::   tag_sw_v2c_xes = "sw_v2c_xes"
  character(len("on")), parameter          ::    tag_on =  "on"
  character(len("off")), parameter          ::    tag_off =  "off"
!
!!!  integer :: sw_corelevel_spectrum = OFF
!
! --------------
!  Core states to probe
!
  character(len("probe")),   parameter  ::   tag_probe   = "probe"
  character(len("atom_id")),   parameter  :: tag_atomid  = "atom_id"
  character(len("spin")),   parameter  :: tag_spin = "spin"
  character(len("orbital")),   parameter  :: tag_orbital  = "orbital"

  integer :: atom_to_probe =0
  integer :: spin_to_probe = 0
  integer :: qnum_n_to_probe =1, qnum_l_to_probe =0     ! K-edge default
  integer :: num_core_states
!
  integer :: ndim_spinor_core_states = 1
!
  real(kind=DP), allocatable :: psig_core_states(:,:,:,:)
  real(kind=DP), allocatable :: fsr_core_states(:,:,:), fsi_core_states(:,:,:)
!
  real(kind=DP), allocatable :: ene_core_states(:)
!
! --------------
! initial/final states
!
  character(len("initial_state_level")),   parameter :: &
       &                        tag_initial_state_level  = "initial_state_level"
  character(len("initial_state_splitting")),   parameter :: &
       &                        tag_initial_state_splitting  = "initial_state_splitting"

  character(len("initial_state_energy")),   parameter :: &
       &                        tag_initial_state_energy  = "initial_state_energy"
  character(len("final_state_energy")),   parameter :: &
       &                        tag_final_state_energy  = "final_state_energy"
  character(len("read_core_energy_from_file")),   parameter :: &
       &                  tag_read_core_energy_from_file = "read_core_energy_from_file"
  character(len("mimic_soc_split_spectrum")),   parameter :: &
       &                  tag_mimic_soc_split_spectrum = "mimic_soc_split_spectrum"

  real(kind=DP) :: ene_initial_state_level = 0.0d0
  real(kind=DP) :: ene_initial_state_splitting = 0.0d0
  real(kind=DP) :: etot_initial_state = 0.0d0
  real(kind=DP) :: etot_final_state = 0.0d0
!
  real(kind=DP) :: etot_change_initial_to_final = 0.0d0
!
  logical :: read_core_energy_from_file = .true.
  logical :: mimic_soc_split_spectrum = .true.
  logical :: initial_state_splitting_is_set = .false.
!
! ---
  real(kind=DP) :: eshift_corelevel_virtual = -0.1d0

! --------------
!  Momentum transfer
!
  character(len("momentum")),   parameter  ::  tag_momentum   = "momentum"
  character(len("qx")),   parameter  :: tag_qx  = "qx"
  character(len("qy")),   parameter  :: tag_qy  = "qy"
  character(len("qz")),   parameter  :: tag_qz  = "qz"
!
  real(kind=DP) :: vec_q(3) = (/0.0D0, 0.0D0, 1.0D0/)
!
  character(len("polarization")),   parameter  ::  tag_polarization   = "polarization"
  character(len("ux")),   parameter  :: tag_ux  = "ux"
  character(len("uy")),   parameter  :: tag_uy  = "uy"
  character(len("uz")),   parameter  :: tag_uz  = "uz"
!
! --------------
!  Energy range
!
  character(len("energy")),   parameter  ::  tag_energy   = "energy"
  character(len("low")),      parameter  ::  tag_low      = "low"
  character(len("high")),     parameter  ::  tag_high     = "high"
  character(len("step")),     parameter  ::  tag_step     = "step"
!
  real(kind=DP) :: e_low = -0.2D0,  e_high = 2.0D0,  e_step = 0.002d0   ! default
  integer :: ndiv_erange
!
! --------------
! Density Response Function
!
  character(len("smearing_width")),   parameter  ::  tag_eta = "smearing_width"

! real(kind=DP) :: eta = 0.01837451d0               ! 0.5 eV
  complex(kind=CMPLDP), allocatable :: Chi0_00(:)
!
! --------------
! AE wave functions of core orbs
!
  integer :: num_core_ae_wfns
  integer, allocatable :: qnum_n_core_ae_wfns(:)
  integer, allocatable :: qnum_l_core_ae_wfns(:)
  real(kind=DP), allocatable :: psir_core_ae_wfns(:,:)
  real(kind=DP), allocatable :: enelevel_core_ae_wfns(:)
  real(kind=DP), allocatable :: focc_core_ae_wfns(:)
  real(kind=DP), allocatable :: level_splitting_core_ae_wfns(:)
!
! --------------
!  Dipole correction term
!
  integer :: num_dipole_dxyz_core2val
  integer, allocatable :: dipole_dxyz_core2val_n1(:),   dipole_dxyz_core2val_t2(:)
  integer, allocatable :: dipole_dxyz_core2val_ylm1(:), dipole_dxyz_core2val_ylm2(:)
  real(kind=DP), allocatable :: dipole_dxyz_core2val(:,:)

contains

  subroutine m_CLS_initialize
    vec_q(1:3) = (/0.0D0, 0.0D0, 1.0D0/)
  end subroutine m_CLS_initialize

  subroutine m_CLS_chk_sw_corelevel_spectrum
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop, iret

    real(kind=DP) :: norm, dret
    character(3) :: str1
    character(len=FMAXVALLEN) :: rstr
    logical :: tf

    if( f_getStringValue( tag_sw_corelevel_spectrum, rstr, LOWER) == 0) then
       call strncmp0( tag_on, trim(rstr), tf)
       if(tf) sw_corelevel_spectrum = on
    end if
    if( f_getStringValue( tag_sw_v2c_xes, rstr, LOWER) == 0) then
       call strncmp0( tag_on, trim(rstr), tf)
       if(tf) sw_v2c_xes = on
    end if
    write(nfout,*) '!** CLS: sw_corelevel_spectrum is ', sw_corelevel_spectrum
    write(nfout,*) '!** CLS: sw_v2c_xes is ', sw_v2c_xes
  end subroutine m_CLS_chk_sw_corelevel_spectrum

  subroutine m_CLS_chk_sw_local_approx_trans
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop, iret

    real(kind=DP) :: norm, dret
    character(3) :: str1
    character(len=FMAXVALLEN) :: rstr
    logical :: tf

    if( f_getStringValue( tag_sw_local_approx_trans_mom, rstr, LOWER) == 0) then
       call strncmp0( tag_on, trim(rstr), tf)
       if(tf) sw_local_approx_trans_moment = on
    end if
    write(nfout,*) '!** CLS: sw_local_approx_trans_moment is ', &
         &                   sw_local_approx_trans_moment
  end subroutine m_CLS_chk_sw_local_approx_trans

  subroutine m_CLS_rd_n
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop, iret

    real(kind=DP) :: norm, dret
    character(3) :: str1
    character(len=FMAXVALLEN) :: rstr

    iret = f_selectTop()
    if ( ipriinputfile >= 2 .and. printable ) then
       write(nfout,'(" !*  tag_corelevel_spectrum")')
    endif

    if ( f_selectBlock( tag_corelevel_spectrum ) == 0 ) then
       sw_corelevel_spectrum = ON
    endif
    if ( sw_corelevel_spectrum == ON ) call m_CLS_rd_n_main

  end subroutine m_CLS_rd_n

  subroutine m_CLS_rd_n_main
    integer :: f_selectBlock, f_getStringValue, f_getRealValue, f_getIntValue
    integer :: f_selectParentBlock, f_selectTop, iret

    real(kind=DP) :: norm, dret
    character(3) :: str1
    character(len=FMAXVALLEN) :: rstr
    logical :: tf

! --
    initial_state_splitting_is_set = .false.

! ----- probe ---
    if ( f_selectBlock( tag_probe ) == 0 ) then
       if ( f_getIntValue( tag_atomid, iret ) == 0 ) then
!          if ( iret < 1 .or. iret > natm ) then
!             write(nfout,*) '!! AtomID is not set'
!             call phase_error_with_msg(nfout,'!! AtomID is not set',__LINE__,__FILE__)
!          else
             atom_to_probe = iret
!          endif
       endif
       if ( f_getIntValue( tag_spin, iret ) == 0 ) then
          spin_to_probe = iret
       endif
       if ( f_getStringValue( tag_orbital, rstr, LOWER ) == 0 ) then
          str1 = trim(rstr)

          if ( str1(1:1) == "1" ) qnum_n_to_probe = 1
          if ( str1(1:1) == "2" ) qnum_n_to_probe = 2
          if ( str1(1:1) == "3" ) qnum_n_to_probe = 3
          if ( str1(1:1) == "4" ) qnum_n_to_probe = 4
          if ( str1(1:1) == "5" ) qnum_n_to_probe = 5
          if ( str1(1:1) == "6" ) qnum_n_to_probe = 6

          if ( str1(2:2) == "s" .or. str1(2:2) == "s" ) qnum_l_to_probe = 0
          if ( str1(2:2) == "p" .or. str1(2:2) == "p" ) qnum_l_to_probe = 1
          if ( str1(2:2) == "d" .or. str1(2:2) == "d" ) qnum_l_to_probe = 2
          if ( str1(2:2) == "f" .or. str1(2:2) == "f" ) qnum_l_to_probe = 3

          if( f_getRealValue( tag_initial_state_level, dret, 'hartree') == 0) then
             ene_initial_state_level  = dret
          endif
          if( f_getRealValue( tag_initial_state_splitting, dret, 'hartree') == 0) then
             ene_initial_state_splitting  = dret
          endif

          if( f_getRealValue( tag_initial_state_energy, dret, 'hartree') == 0) then
             etot_initial_state  = dret
          endif
          if( f_getRealValue( tag_final_state_energy, dret, 'hartree') == 0) then
             etot_final_state  = dret
          endif

          if( f_getRealValue( tag_initial_state_splitting, dret, 'hartree') == 0) then
             ene_initial_state_splitting  = dret
             initial_state_splitting_is_set = .true.
          endif

          if ( f_getStringValue( tag_read_core_energy_from_file, rstr, LOWER )==0 ) then
             call strncmp0( tag_on, trim(rstr), tf )
             if ( tf ) read_core_energy_from_file = .true.
             call strncmp0( tag_off, trim(rstr), tf )
             if ( tf ) read_core_energy_from_file = .false.
          endif
          if ( f_getStringValue( tag_mimic_soc_split_spectrum, rstr, LOWER )==0 ) then
             call strncmp0( tag_on, trim(rstr), tf )
             if ( tf ) mimic_soc_split_spectrum = .true.
             call strncmp0( tag_off, trim(rstr), tf )
             if ( tf ) mimic_soc_split_spectrum = .false.
          endif

       endif
! ----- energy range ----
       if ( f_selectBlock( tag_energy ) == 0 ) then
          if( f_getRealValue( tag_low, dret, 'hartree') == 0) e_low  = dret
          if( f_getRealValue( tag_high, dret,'hartree') == 0) e_high = dret
          if( f_getRealValue( tag_step, dret,'hartree') == 0) e_step = dret
          iret = f_selectParentBlock()
       endif
! ----- momentum transfer
       if ( f_selectBlock( tag_momentum ) == 0 ) then
          if ( f_getRealValue( tag_qx, dret," " ) == 0 ) vec_q(1) = dret
          if ( f_getRealValue( tag_qy, dret," " ) == 0 ) vec_q(2) = dret
          if ( f_getRealValue( tag_qz, dret," " ) == 0 ) vec_q(3) = dret
          iret = f_selectParentBlock()
       endif
       if ( f_selectBlock( tag_polarization ) == 0 ) then
          if ( f_getRealValue( tag_ux, dret," " ) == 0 ) vec_q(1) = dret
          if ( f_getRealValue( tag_uy, dret," " ) == 0 ) vec_q(2) = dret
          if ( f_getRealValue( tag_uz, dret," " ) == 0 ) vec_q(3) = dret
          iret = f_selectParentBlock()
       endif

       iret = f_selectParentBlock()
    endif

    norm = sqrt( vec_q(1)**2 +vec_q(2)**2 +vec_q(3)**2 )
    if ( norm > 1.0D-8 ) then
       vec_q = vec_q /norm
    else
       vec_q = 0.0;  vec_q(3) = 1.0D0
    endif

    if ( etot_initial_state /= 0.0d0 .and. etot_final_state /= 0.0d0 ) then
       etot_change_initial_to_final = etot_final_state -etot_initial_state
    endif

    if ( qnum_l_to_probe > 0 ) then
       if ( nspin == 1 ) mimic_soc_split_spectrum = .true.
    endif

! === kT_add === 2014/09/22
    if ( ndim_spinor == 2 ) mimic_soc_split_spectrum = .false.
! ============== 2014/09/22

    num_core_states = 2 *qnum_l_to_probe +1

    if ( noncol ) then
       num_core_states = 2 *num_core_states
       ndim_spinor_core_states = 2
    else
       if ( qnum_l_to_probe > 0 ) then
          if ( .not. mimic_soc_split_spectrum ) then
             num_core_states = 2 *num_core_states
             ndim_spinor_core_states = 2
          endif
       endif
    endif

    if ( printable ) then
       write(nfout,*)
       write(nfout,*) '!** CLS: ****** Configuration of CoreLevel Spectrum *****'
       write(nfout,*) '!** CLS: Atom ID =', atom_to_probe
       write(nfout,*) '!** CLS: spin =', spin_to_probe
       write(nfout,*) '!** CLS: qnum_n, qnum_l =', qnum_n_to_probe, qnum_l_to_probe
       write(nfout,'(A,3F15.8)') ' !** CLS: Momentum transfer direction = ', &
            &                    vec_q(1), vec_q(2), vec_q(3)
       write(nfout,'(3(A,F15.8))') &
            &          ' !** CLS: energy range =  low: ', e_low, " high: ", e_high, &
            &          " step: ", e_step
       if ( ene_initial_state_level /= 0.0d0 ) then
          write(nfout,'(A,F15.8)') " !** CLS: initial_state_level = ", &
               &              ene_initial_state_level
       endif
       if ( ene_initial_state_splitting /= 0.0d0 ) then
          write(nfout,'(A,F15.8)') " !** CLS: initial_state_splitting = ", &
               &                        ene_initial_state_splitting
       endif

       if ( etot_initial_state /= 0.0d0 ) then
          write(nfout,'(A,F15.8)') " !! initial state total energy= ",etot_initial_state
       endif
       if ( etot_final_state /= 0.0d0 ) then
          write(nfout,'(A,F15.8)') " !! final   state total energy= ",etot_final_state
       endif
       if ( etot_change_initial_to_final /= 0.0d0 ) then
          write(nfout,'(A,F15.8)') " !! etot_change_initial_to_final = ", &
               &                    etot_change_initial_to_final
       endif

       write(nfout,*) '!** CLS: read_core_energy_from_file is ', &
            &                read_core_energy_from_file
       write(nfout,*) '!** CLS: mimic_soc_split_spectrum is ', mimic_soc_split_spectrum

       write(nfout,*) '!** CLS: num_core_states = ', num_core_states
       write(nfout,*) '!** CLS: ************************************************ '
       write(nfout,*)
    endif

  end subroutine m_CLS_rd_n_main

! === KT_add === 2014/10/02
  integer function m_CLS_find_orb_index_to_probe()
    integer :: index, i

    index = 0
    Do i=1, num_core_ae_wfns
       if ( qnum_n_core_ae_wfns(i) == qnum_n_to_probe .and. &
            qnum_l_core_ae_wfns(i) == qnum_l_to_probe ) then
          index = i
          exit
       endif
    End Do
    m_CLS_find_orb_index_to_probe = index

  end function m_CLS_find_orb_index_to_probe
! ============= 2014/10/02

  subroutine m_CLS_read_core_energy
    real(kind=DP) :: c1

    call m_Files_open_core_ene_initial
    call m_Files_open_core_ene_final

    if ( mype == 0 ) then
       call set_splitting_and_etotal( nfcore_energy_initial, c1, etot_initial_state )
       ene_initial_state_splitting = c1

       call set_splitting_and_etotal( nfcore_energy_final,   c1, etot_final_state   )
       etot_change_initial_to_final = etot_final_state -etot_initial_state

       write(nfout,*) '**********************************************'
       write(nfout,*) '** etot_initial_state is read and set to ', &
            &             etot_initial_state
       write(nfout,*) '** etot_final_state is read and set to ', &
            &             etot_final_state
       write(nfout,*) '** ene_initial_state_splitting is read and set to ', &
            &             ene_initial_state_splitting
       write(nfout,*) '**********************************************'
    endif

    if ( npes > 1 ) then
       call mpi_bcast( ene_initial_state_splitting, 1, mpi_double_precision, &
            &          0, mpi_comm_world, ierr )
       call mpi_bcast( etot_change_initial_to_final, 1, mpi_double_precision, &
            &          0, mpi_comm_world, ierr )
    endif

    call m_Files_close_core_ene_initial
    call m_Files_close_core_ene_final

  contains

    subroutine set_splitting_and_etotal( lun, ene_splitting, ene_total )
      integer, intent(in) :: lun
      real(kind=DP), intent(out) :: ene_splitting, ene_total

      integer :: n1, n2, n3, n4, ia, it
      real(kind=DP) :: c1
      logical :: found

      ene_splitting = 0.0d0;  ene_total = 0.0d0
      ia = atom_to_probe;    it = ityp(ia)

      found = .false.

      read(lun,*);   read(lun,*)
      Do while ( .true. )
         read(lun,*) n1, n2, n3, n4, c1
#if 0
         if ( n1 == it .and. n2 == ia .and. n3 == qnum_n_to_probe &
              &                       .and. n4 == qnum_l_to_probe ) then
#else
         if ( n2 == ia .and. n3 == qnum_n_to_probe &
              &        .and. n4 == qnum_l_to_probe ) then
#endif
            ene_splitting = c1
            found = .true.
         endif
         if ( n1 == 0 ) exit
      End do

      read(lun,*);  read(lun,*);   read(lun,*) ene_total

      if ( .not. found ) then
         if ( qnum_l_to_probe > 0 ) call phase_error_with_msg(nfout,"ene_splitting is not found",__LINE__,__FILE__)
      endif

    end subroutine set_splitting_and_etotal

  end subroutine m_CLS_read_core_energy

! ==== ONLY TEST ===
  subroutine m_CLS_set_corelevels_virtually
    integer :: orb_index, ia, it, i
    real(kind=DP) :: ene_lower, ene_higher
    real(kind=DP), allocatable :: level_splitting(:)

    ia = atom_to_probe;  it = ityp(ia)

    allocate( level_splitting( num_core_ae_wfns ) )
    level_splitting = 0.0d0

    if ( flg_paw .and. eval_core_level_splitting == ByPawPot ) then
       call m_SO_calc_corelevel_splitting( ia, it, &
            &                              num_core_ae_wfns, &
            &                              qnum_l_core_ae_wfns, &
            &                              nmesh(it), psir_core_ae_wfns, &
            &                              level_splitting )
    endif
!
    orb_index = m_CLS_find_orb_index_to_probe()
    ene_higher = efermi +eshift_corelevel_virtual
    ene_lower  = ene_higher -level_splitting( orb_index )
!
    if ( noncol ) then
       Do i=1, num_core_states
          if ( i <= 2*qnum_l_to_probe ) then
             ene_core_states(i) = ene_lower
          else
             ene_core_states(i) = ene_higher
          endif
       End Do
    else
       ene_core_states = ene_higher
    endif

    deallocate( level_splitting )

  end subroutine m_CLS_set_corelevels_virtually

#if 0
  subroutine m_CLS_estimate_corelevels
    integer :: ia, it, j, orb_index
    real(kind=DP) :: e_level_hard, e_level_soc

    real(kind=DP), allocatable :: e_levels_soft(:)

    ia = atom_to_probe;    it = ityp(ia)
    orb_index = m_CLS_find_orb_index_to_probe()

    allocate( e_levels_soft( num_core_states ) )
    call calc_softpart_contrib

!    write(*,*) "orb_index = ", orb_index

    call m_SO_calc_contrib_corelevels( ia, it, num_core_ae_wfns, &
         &                             orb_index, qnum_l_to_probe, &
         &                             nmesh(it), psir_core_ae_wfns, &
         &                             e_level_hard, e_level_soc )
    if ( mype == 0 ) then
       Do j=1, num_core_states
          write(9000+mype,*) j, e_levels_soft(j), e_level_hard, e_level_soc
       End do
    endif

    deallocate( e_levels_soft )

  contains

    subroutine calc_softpart_contrib
      integer :: ib, ik, ispin, i, ri, iend, i1, ig
      real(kind=DP) :: fac, epot, ekin, ekin_sum, c1
      real(kind=DP), allocatable :: afft(:), bfft(:), qxyz(:,:)
      real(kind=DP), allocatable :: psi_l(:,:,:,:), tmp_chgq_l(:,:,:)

!      fac = 2.d0 / (univol*kv3*product(fft_box_size_WF(1:3,1)))
      fac = 2.d0 / (univol*product(fft_box_size_WF(1:3,1)))

      allocate( afft(nfft) ); afft = 0.0d0
      allocate( bfft(nfft) ); bfft = 0.0d0

      allocate( qxyz(kg1,3 ) ); qxyz = 0.0d0
      allocate( tmp_chgq_l( ista_kngp:iend_kngp,kimg,nspin) ); tmp_chgq_l = 0.0d0

      call m_FFT_alloc_WF_work()

      Do ib=1, num_core_states
         afft = 0.0d0

         ekin = 0.0d0
         Do ispin=1, nspin
            Do ik = ispin, kv3+ispin-nspin, nspin
!            Do ik = 1, 1
               if ( map_k(ik) /= myrank_k ) cycle

               call k_plus_G_vectors_m( ik, kgp, kg1, kv3, iba, nbase, vkxyz, ngabc, &
                    &                   rltv, qxyz(:,1), qxyz(:,2), qxyz(:,3) )

               Do ig=1, iba(ik)
                  c1 = qxyz(ig,1)**2 + qxyz(ig,2)**2 + qxyz(ig,3)**2
                  c1 = c1 *qwgt(ik) /2.0d0
                  Do ri=1, kimg
                     ekin = ekin +psig_core_states(ig,ib,ik,ri)**2 *c1
                  End do
               End do

               allocate( psi_l( kg1, 1, ik:ik, kimg ) ); psi_l = 0.0d0
               psi_l(1:iba(ik),1,ik,:) = psig_core_states(1:iba(ik),ib,ik,:)

               call m_ES_WF_in_Rspace_kt( ik, ik, ik, psi_l, bfft )

               do i = 1, nfft-1, 2
                  afft(i) = afft(i) + ( bfft(i)**2 +bfft(i+1)**2 ) *qwgt(ik)
!                  afft(i) = afft(i) + ( bfft(i)**2 +bfft(i+1)**2 )
               end do
               deallocate( psi_l )
            End do
         End Do

         if ( npes > 1 ) then
            call mpi_allreduce( afft, bfft, nfft, mpi_double_precision, mpi_sum, &
                 &              MPI_CommGroup, ierr )
            afft = bfft /dble(nrank_e)
         endif

         call m_FFT_WF( ELECTRON, nfout, afft, DIRECT, OFF )

         Do ispin=1, nspin
            do ri = 1, kimg
               iend = iend_kngp
               if( iend_kngp > kg ) iend = kg
               if( ista_kngp <= iend ) then
                  do i = ista_kngp, iend  !for mpi
                     i1 = kimg*igf(i) + (ri - kimg)
                     tmp_chgq_l(i,ri,ispin) = afft(i1) *fac

                     write(9020+mype,*) i, ri, ispin, tmp_chgq_l(i,ri,ispin)*univol, vlhxc_l(i,ri,ispin)
                  end do
               endif
            end do
         End Do

         epot = 0.0d0
         Do ispin=1, nspin
            Do ri=1, kimg
               Do i=ista_kngp, iend_kngp
                  epot = epot +tmp_chgq_l(i,ri,ispin)*vlhxc_l(i,ri,ispin)
               End do
            End Do
         End do

         if ( npes > 1 ) then
            call mpi_allreduce( ekin, ekin_sum, 1, mpi_double_precision, mpi_sum, &
                 &              MPI_CommGroup, ierr )
            ekin = ekin_sum /dble(nrank_e)
         endif

         epot = epot /2.0d0 *univol

         e_levels_soft(ib) = epot +ekin
!
         write(9000+mype,*) "ib csum = ", ib, epot, ekin
      End Do
      call mpi_barrier( MPI_CommGroup, ierr )

      deallocate( afft ); deallocate( bfft ); deallocate( tmp_chgq_l )
      call m_FFT_dealloc_WF_work()

    end subroutine calc_softpart_contrib

  end subroutine m_CLS_estimate_corelevels
#endif

  subroutine m_CLS_set_data_core2val_from_pp
    integer :: ia, it
    integer :: nfp

!    if ( mype == 0 ) call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
!    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
!    if (ierr/=0) call mpi_stop(nfout)

    ia = atom_to_probe;  it = ityp(ia)
    write(nfout,*) 'Reading core info from pp, ia = ', ia

! == KT_DEBUG === 2015/06/15 ==
!    call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
!                                      ! this does not work properly by unknown reason.
    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
! =============== 2015/06/15

    if (ierr/=0) call mpi_stop(nfout)

    nfp = nfpot(it)

    if ( mype == 0 ) call read_num_core_ae_wfns( nfp, num_core_ae_wfns, it )
    if ( npes > 1 ) then
       call mpi_bcast(  num_core_ae_wfns, 1, mpi_integer, 0, mpi_comm_world, ierr )
    endif

    if ( num_core_ae_wfns == 0 ) return

    allocate( qnum_n_core_ae_wfns( num_core_ae_wfns ) )
    allocate( qnum_l_core_ae_wfns( num_core_ae_wfns ) )
    allocate( psir_core_ae_wfns( nmesh(it),  num_core_ae_wfns ) )
    allocate( enelevel_core_ae_wfns( num_core_ae_wfns ) )
    allocate( focc_core_ae_wfns( num_core_ae_wfns ) )
    allocate( level_splitting_core_ae_wfns( num_core_ae_wfns ) )

    if ( mype == 0 ) then
       call read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh(it), &
            &                       qnum_n_core_ae_wfns, qnum_l_core_ae_wfns, &
            &                       psir_core_ae_wfns, &
            &                       enelevel_core_ae_wfns, focc_core_ae_wfns )
    endif

    if ( npes > 1 ) then
       call mpi_bcast( qnum_n_core_ae_wfns, num_core_ae_wfns, &
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( qnum_l_core_ae_wfns, num_core_ae_wfns, &
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( enelevel_core_ae_wfns, num_core_ae_wfns, &
            &          mpi_double_precision, 0, mpi_comm_world, ierr )
       call mpi_bcast( focc_core_ae_wfns, num_core_ae_wfns, &
            &          mpi_double_precision, 0, mpi_comm_world, ierr )
       call mpi_bcast( psir_core_ae_wfns, nmesh(it)*num_core_ae_wfns, &
            &          mpi_double_precision, 0, mpi_comm_world, ierr )
    endif

! === KT_add ===== 2015/01/23
    if ( eval_core_level_splitting == ReadFromPP ) then
       if ( mype == 0 ) then
          call read_data_core_level_splitting( nfp, it, num_core_ae_wfns, &
               &                               qnum_n_core_ae_wfns, &
               &                               qnum_l_core_ae_wfns, &
               &                               level_splitting_core_ae_wfns )
       endif
       if ( npes > 1 ) then
          call mpi_bcast( level_splitting_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
       endif
    endif
! =============== 2015/01/23

    if ( mype == 0 ) then
       call read_num_dipole_dxyz_core2val( nfp, num_dipole_dxyz_core2val,it )
    endif
    if ( npes > 1 ) then
       call mpi_bcast(  num_dipole_dxyz_core2val, 1, &
            &           mpi_integer, 0, mpi_comm_world, ierr )
    endif

    allocate( dipole_dxyz_core2val_n1( num_dipole_dxyz_core2val ) )
    allocate( dipole_dxyz_core2val_t2( num_dipole_dxyz_core2val ) )
    allocate( dipole_dxyz_core2val_ylm1( num_dipole_dxyz_core2val ) )
    allocate( dipole_dxyz_core2val_ylm2( num_dipole_dxyz_core2val ) )
    allocate( dipole_dxyz_core2val( num_dipole_dxyz_core2val,3 ) )

    if ( mype == 0 ) then
       if ( num_orb_index_data(it) > 0 ) then
          call read_dipole_dxyz_core2val_A( nfp, num_dipole_dxyz_core2val, &
               &                            dipole_dxyz_core2val_n1, &
               &                            dipole_dxyz_core2val_t2, &
               &                            dipole_dxyz_core2val_ylm1, &
               &                            dipole_dxyz_core2val_ylm2, &
               &                            dipole_dxyz_core2val, it )
       else
          call read_dipole_dxyz_core2val_B( nfp, num_dipole_dxyz_core2val, &
               &                            dipole_dxyz_core2val_n1, &
               &                            dipole_dxyz_core2val_t2, &
               &                            dipole_dxyz_core2val_ylm1, &
               &                            dipole_dxyz_core2val_ylm2, &
               &                            dipole_dxyz_core2val )
       endif
    endif
!
    if ( npes > 1 ) then
       call mpi_bcast( dipole_dxyz_core2val_n1, num_dipole_dxyz_core2val,&
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( dipole_dxyz_core2val_t2, num_dipole_dxyz_core2val,&
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( dipole_dxyz_core2val_ylm1, num_dipole_dxyz_core2val,&
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( dipole_dxyz_core2val_ylm2, num_dipole_dxyz_core2val,&
            &          mpi_integer, 0, mpi_comm_world, ierr )
       call mpi_bcast( dipole_dxyz_core2val, num_dipole_dxyz_core2val*3,&
            &          mpi_double_precision, 0, mpi_comm_world, ierr )
    endif

    call m_Files_close_ps_file(it)

  end subroutine m_CLS_set_data_core2val_from_pp

  subroutine m_CLS_alloc_ene_core_states
    if ( .not. allocated(ene_core_states) ) then
       allocate( ene_core_states( num_core_states ) )
       ene_core_states = 0.0d0
    endif
  end subroutine m_CLS_alloc_ene_core_states

  subroutine m_CLS_alloc_wfn_core_states( psig_required )
    logical, intent(in) :: psig_required

    if ( allocated( psig_core_states ) ) deallocate( psig_core_states )
    if ( allocated( fsr_core_states ) ) deallocate( fsr_core_states )
    if ( allocated( fsi_core_states ) ) deallocate( fsi_core_states )

    if ( psig_required ) then
       allocate( psig_core_states( maxval(np_g1k), num_core_states, ista_k:iend_k, 2 ) )
       psig_core_states = 0.0d0
    endif
!
    allocate( fsr_core_states(num_core_states, 2*qnum_l_to_probe+1, ista_k:iend_k) )
    allocate( fsi_core_states(num_core_states, 2*qnum_l_to_probe+1, ista_k:iend_k ) )
    fsr_core_states = 0.0d0;  fsi_core_states = 0.0d0
!
  end subroutine m_CLS_alloc_wfn_core_states

  subroutine m_CLS_dealloc_ene_core_states
    if ( allocated( ene_core_states ) )  deallocate( ene_core_states )
  end subroutine m_CLS_dealloc_ene_core_states

  subroutine m_CLS_dealloc_wfn_core_states
    if ( allocated( psig_core_states ) ) deallocate( psig_core_states )
    if ( allocated( fsr_core_states ) )  deallocate( fsr_core_states )
    if ( allocated( fsi_core_states ) )  deallocate( fsi_core_states )
  end subroutine m_CLS_dealloc_wfn_core_states

  subroutine m_CLS_dealloc_core_ae_wfns
    if ( allocated( psir_core_ae_wfns ) ) deallocate( psir_core_ae_wfns )
    if ( allocated( enelevel_core_ae_wfns ) ) deallocate( enelevel_core_ae_wfns )
    if ( allocated( focc_core_ae_wfns ) )   deallocate( focc_core_ae_wfns )
    if ( allocated( qnum_n_core_ae_wfns ) ) deallocate( qnum_n_core_ae_wfns )
    if ( allocated( qnum_l_core_ae_wfns ) ) deallocate( qnum_l_core_ae_wfns )
    if ( allocated( level_splitting_core_ae_wfns ) ) then
       deallocate( level_splitting_core_ae_wfns )
    endif
  end subroutine m_CLS_dealloc_core_ae_wfns

  subroutine m_CLS_dealloc_dipole_core2val
    deallocate(dipole_dxyz_core2val_n1)
    deallocate(dipole_dxyz_core2val_t2)
    deallocate(dipole_dxyz_core2val_ylm1)
    deallocate(dipole_dxyz_core2val_ylm2)
    deallocate(dipole_dxyz_core2val)
  end subroutine m_CLS_dealloc_dipole_core2val

  subroutine m_CLS_set_ene_core_states
    integer :: i, ia, it, il1, index
    real(kind=DP) :: c1, c2, c0
    real(kind=DP), allocatable :: level_splitting(:)

    ia = atom_to_probe;   it = ityp(ia)
    il1 = qnum_l_to_probe +1

    index = 0
    Do i=1, num_core_ae_wfns
       if ( qnum_n_core_ae_wfns(i) == qnum_n_to_probe .and. &
            qnum_l_core_ae_wfns(i) == qnum_l_to_probe ) then
          index = i
          exit
       endif
    End Do
    if ( index == 0 ) call phase_error_with_msg(nfout,"index of core states is not found",__LINE__,__FILE__)

    ene_core_states = enelevel_core_ae_wfns(index)

! -------------
! set edge energy
! -------------
    if ( flg_paw .and. read_core_energy_from_file ) then        ! set automatically
       call m_CLS_read_core_energy
       if ( sw_v2c_xes == ON ) then
          ene_core_states = efermi +etot_change_initial_to_final     !! ????
       else
          ene_core_states = efermi -etot_change_initial_to_final     !! ????
       endif
!!!       ene_core_states = -etot_change_initial_to_final

    else
       if ( ene_initial_state_level /= 0.0 ) then
          ene_core_states = ene_initial_state_level
       endif
       if ( etot_change_initial_to_final /= 0.0d0 ) then
          ene_core_states = efermi -etot_change_initial_to_final
       endif
    endif

! ----------
! set spinorbit strength
! ----------
    if ( flg_paw .and. read_core_energy_from_file ) then
    else
       if ( .not. initial_state_splitting_is_set ) then
          if ( SpinOrbit_mode == ReadFromPP ) then
             ene_initial_state_splitting = level_splitting_core_ae_wfns( index )

          else if ( SpinOrbit_mode == ByPawPot ) then
             allocate( level_splitting( num_core_ae_wfns ) ); level_splitting = 0.0d0
             call m_SO_calc_corelevel_splitting( ia, it, &
                  &                              num_core_ae_wfns, &
                  &                              qnum_l_core_ae_wfns, &
                  &                              nmesh(it), psir_core_ae_wfns, &
                  &                              level_splitting )
             ene_initial_state_splitting = level_splitting( index )
             deallocate( level_splitting )
          endif
          write(nfout,*) "!** CLS: ene_initial_state_splitting = ", &
               &          ene_initial_state_splitting
       endif
    endif

    if ( qnum_l_to_probe > 0 ) then
       if ( mimic_soc_split_spectrum ) then
          c0 = ene_initial_state_splitting / (2 *qnum_l_to_probe +1)
          c1 = -c0 *( qnum_l_to_probe +1 )
          c2 =  c0 *qnum_l_to_probe
          ene_core_states = ene_core_states +c2

       else
          c0 = ene_initial_state_splitting / (2 *qnum_l_to_probe +1)
          c1 = -c0 *( qnum_l_to_probe +1 )
          c2 =  c0 *qnum_l_to_probe
!
          Do i=1, num_core_states
             if ( i<= qnum_l_to_probe*2 ) then
                ene_core_states(i) = ene_core_states(i) +c1
             else
                ene_core_states(i) = ene_core_states(i) +c2
             endif
          End Do
       endif
    endif

  end subroutine m_CLS_set_ene_core_states

  subroutine m_CLS_set_wfn_core_states( psig_required )         ! psig_core_stats
    logical, intent(in) :: psig_required

    integer :: i, ia, it, il1
    integer :: index

    real(kind=DP), allocatable :: psig_core_orb(:,:,:)

    ia = atom_to_probe;   it = ityp(ia)
    il1 = qnum_l_to_probe +1

    index = 0
    Do i=1, num_core_ae_wfns
       if ( qnum_n_core_ae_wfns(i) == qnum_n_to_probe .and. &
            qnum_l_core_ae_wfns(i) == qnum_l_to_probe ) then
          index = i
          exit
       endif
    End Do
    if ( index == 0 ) call phase_error_with_msg(nfout,"index of core states is not found",__LINE__,__FILE__)

    if ( psig_required ) then
       allocate( psig_core_orb( maxval(np_g1k), num_core_states/ndim_spinor_core_states, kv3/nspin ) )
       psig_core_orb = 0.0d0;

       call fourier_trans_of_core_ae_wfns
    endif

    if ( noncol ) then
       call set_core_states_with_phase_A
    else
       if ( ndim_spinor_core_states == 1 ) then
          call set_core_states_with_phase
       else
          call set_core_states_with_phase_A
       endif
    end if

    if ( allocated(psig_core_orb) ) deallocate( psig_core_orb );

  contains

    subroutine fourier_trans_of_core_ae_wfns
      integer :: ik, iksnl, im1, ig, nspher, n
      real(kind=DP), allocatable :: snl2(:), wka(:), wkb(:), ylm(:), vlength(:)
      real(kind=DP), allocatable :: qx(:), qy(:), qz(:)
      real(kind=DP) :: fac, facr, norm, hn

      allocate( snl2(kg1) ); allocate( vlength(kg1) );
      allocate( wka(kg1) ); allocate( wkb(kg1) );  allocate( ylm(kg1) )
      allocate( qx(kg1) ); allocate( qy(kg1) ); allocate( qz(kg1) );
      if ( allocated(radr) ) deallocate (radr)
      if ( allocated(wos) ) deallocate (wos)
      allocate( radr(nmesh(it) ) ); allocate( wos(nmesh(it) ) )

      fac = PAI4/dsqrt(univol)

! === KT_add == 2015/06/16
      call rmeshs(nmesh(it),nmesh(it),xh(it),rmax(it),radr,hn) ! -(b_PP)
      call coef_simpson_integration(nmesh(it),nmesh(it),xh(it),radr,wos) ! -(b_PP)
! ============= 2015/06/16

      Do ik=1, kv3, nspin
         if(map_k(ik) /= myrank_k) cycle
         iksnl =(ik-1) /nspin +1

         call k_plus_G_vectors_3D( ik, kg, kg1, kv3, iba, nbase, vkxyz, ngabc, rltv,&
              &                 qx, qy, qz, vlength )
!         call new_radr_and_wos(ik,it)

         Do im1=1, 2*qnum_l_to_probe +1
            nspher = ( il1 -1 )**2 +im1

            call sphr(np_g1k(ik),nspher,qx,qy,qz,ylm)        ! -(bottom_Subr.)
            snl2 = 0.0d0

            do n=1, nmesh(it)
               facr = fac*wos(n)*radr(n)*psir_core_ae_wfns( n, index )
               do ig = 1, np_g1k(ik)
                  wka(ig) = vlength(ig)*radr(n)
               end do
               call dsjnv(il1-1,np_g1k(ik),wka,wkb)     ! -(bottom_Subr.)
               do ig = 1, np_g1k(ik)
                  snl2(ig) = snl2(ig) + facr *wkb(ig) *ylm(ig)
               end do
            end do
#if 0
! --- if normailzation is required ---
            norm = sum( snl2(1:kg1)*snl2(1:kg1) )
            do ig = 1, iba(ik)
               psig_core_orb(ig,im1,iksnl) = snl2(ig) /sqrt(norm)
            end do
#else
            do ig = 1, np_g1k(ik)
               psig_core_orb(ig,im1,iksnl) = snl2(ig)
            end do
#endif
         End do
      End do

      deallocate( snl2 ); deallocate( vlength );
      deallocate( wka ); deallocate( wkb );  deallocate( ylm )
      deallocate( qx ); deallocate( qy ); deallocate( qz );
      deallocate( radr ); deallocate( wos )

    end subroutine fourier_trans_of_core_ae_wfns

    subroutine set_core_states_with_phase
      integer::ik, iksnl,  mil, im1, ig, ii, ib
      real(kind=DP) :: ph, f1, f2, f3, ga, gb, gc
      real(kind=DP), allocatable :: zfcos(:), zfsin(:), tmp_wfn(:,:,:)

      if ( .not. psig_required ) goto 1000

      mil = mod( il1, 4 )

      allocate( zfcos(maxval(np_g1k)) );  allocate( zfsin(maxval(np_g1k)) )

      Do ik=1, kv3
         if(map_k(ik) /= myrank_k) cycle
         iksnl = ( ik-1 ) /nspin +1

         f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2

         Do ig = 1, np_g1k(ik)
            ii = nbase( ig+ista_g1k(ik)-1, ik )
!            if ( ii <= 0 ) cycle
!            ga = vkxyz(ik,1,BUCS) + real(ngabc(ii,1),kind=DP)
!            gb = vkxyz(ik,2,BUCS) + real(ngabc(ii,2),kind=DP)
!            gc = vkxyz(ik,3,BUCS) + real(ngabc(ii,3),kind=DP)
            ga = real(ngabc(ii,1),kind=DP)
            gb = real(ngabc(ii,2),kind=DP)
            gc = real(ngabc(ii,3),kind=DP)
            ph = ga *f1 +gb*f2 + gc*f3
            zfcos(ig) = dcos(ph);    zfsin(ig) = -dsin(ph)
         End do

         Do ib=1, num_core_states
            im1 = ib
            Do ig=1, np_g1k(ik)
               if ( kimg == 1 ) then
!                  psig_core_states(i,ib,ik,1) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                  if ( mil == 1 .or. mil ==3 ) then
                     psig_core_states(ig,ib,ik,1) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                     psig_core_states(ig,ib,ik,2) = zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                  else
                     psig_core_states(ig,ib,ik,1) =-zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                     psig_core_states(ig,ib,ik,2) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                  endif
               else
                  if ( mil == 1 .or. mil ==3 ) then
                     psig_core_states(ig,ib,ik,1) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                     psig_core_states(ig,ib,ik,2) = zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                  else
                     psig_core_states(ig,ib,ik,1) =-zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                     psig_core_states(ig,ib,ik,2) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                  endif
               endif
            End do

         End Do

      End Do

      deallocate( zfcos, zfsin )

1000  continue
      Do ik=1, kv3
         if(map_k(ik) /= myrank_k) cycle
         Do ib=1, num_core_states
            im1 = ib
            fsr_core_states(ib,im1,ik) = 1.0d0
            fsi_core_states(ib,im1,ik) = 0.0d0
         End do
      End do

    end subroutine set_core_states_with_phase

    subroutine set_core_states_with_phase_A             ! nspin=2 only
      integer::ik, iksnl,  mil, im1, ig, ii, ib, jb, is
      integer :: ndim_spinor_bkup
      integer :: ndim_chgpot_bkup

      real(kind=DP) :: ph, f1, f2, f3, ga, gb, gc
      complex(kind=CMPLDP) :: z1, z2

      real(kind=DP), allocatable :: zfcos(:), zfsin(:), tmp_wfn(:,:,:)
!
      if ( .not. noncol ) then
         ndim_spinor_bkup = ndim_spinor;    ndim_spinor = 2
         ndim_chgpot_bkup = ndim_chgpot;    ndim_chgpot = 4
      endif
!
      call m_SO_set_MatU_ylm_RC
      call m_SO_calc_MatLS_orb_s_to_f
      call m_SO_diagonalize_MatLS
!
      if ( .not. noncol ) then
         ndim_spinor = ndim_spinor_bkup
         ndim_chgpot = ndim_chgpot_bkup
      endif
!
      if ( .not. psig_required ) goto 1000

      mil = mod( il1, 4 )
      allocate( zfcos( maxval(np_g1k) ) );
      allocate( zfsin( maxval(np_g1k) ) )
      allocate( tmp_wfn( maxval(np_g1k),num_core_states/ndim_spinor_core_states,2  ) )

      psig_core_states = 0.0d0

      Do ik=1, kv3, nspin
         if(map_k(ik) /= myrank_k) cycle
         iksnl = ( ik-1 ) /nspin +1

         f1 = pos(ia,1)*PAI2; f2 = pos(ia,2)*PAI2; f3 = pos(ia,3)*PAI2
         Do ig = 1, np_g1k(ik)
            ii = nbase( ig+ista_g1k(ik)-1, ik )
            ga = real(ngabc(ii,1),kind=DP)
            gb = real(ngabc(ii,2),kind=DP)
            gc = real(ngabc(ii,3),kind=DP)
            ph = ga *f1 +gb*f2 + gc*f3
            zfcos(ig) = dcos(ph);    zfsin(ig) = -dsin(ph)
         End do

         tmp_wfn = 0.0d0

         Do ib=1, num_core_states /ndim_spinor_core_states
            im1 = ib
            Do ig=1, np_g1k(ik)
               if ( kimg == 1 ) then
                  if ( mil == 1 .or. mil ==3 ) then
                     tmp_wfn(ig,ib,1) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                     tmp_wfn(ig,ib,2) = zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                  else
                     tmp_wfn(ig,ib,1) =-zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                     tmp_wfn(ig,ib,2) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                  endif
               else
                  if ( mil == 1 .or. mil ==3 ) then
                     tmp_wfn(ig,ib,1) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                     tmp_wfn(ig,ib,2) = zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                  else
                     tmp_wfn(ig,ib,1) =-zfsin(ig)*psig_core_orb(ig,im1,iksnl)
                     tmp_wfn(ig,ib,2) = zfcos(ig)*psig_core_orb(ig,im1,iksnl)
                  endif
               endif
            End do
         End Do
!
         Do ib=1, num_core_states
            Do jb=1, num_core_states
               im1 = mod( jb-1, 2*qnum_l_to_probe+1 ) +1
               is = int( (jb-1)/(2*qnum_l_to_probe+1) ) +1

               z1 = 0.0d0
               if ( qnum_l_to_probe == 1 ) then
                  z1 = EigenWfns_MatLS_L1(jb,ib)
               else if ( qnum_l_to_probe == 2 ) then
                  z1 = EigenWfns_MatLS_L2(jb,ib)
               else if ( qnum_l_to_probe == 3 ) then
                  z1 = EigenWfns_MatLS_L3(jb,ib)
               else
                  if ( ib == jb ) z1 = 1.0d0
               endif
!!!!!!!!!!               z1 = conjg(z1)

               Do ig=1, np_g1k(ik)
                  z2 = z1 *dcmplx(tmp_wfn(ig,im1,1), tmp_wfn(ig,im1,2) )
                  psig_core_states(ig,ib,ik+is-1,1) &
                       &      = psig_core_states(ig,ib,ik+is-1,1) +real(z2)
                  psig_core_states(ig,ib,ik+is-1,2) &
                       &      = psig_core_states(ig,ib,ik+is-1,2) +aimag(z2)
               End Do
            End Do
         End Do

      End Do

      deallocate( zfcos, zfsin );    deallocate( tmp_wfn )

1000  continue

      Do ik=1, kv3, nspin
         if(map_k(ik) /= myrank_k) cycle

         Do ib=1, num_core_states
            Do jb=1, num_core_states
               im1 = mod( jb-1, 2*qnum_l_to_probe+1 ) +1
               is = int( (jb-1)/(2*qnum_l_to_probe+1) ) +1

               z1 = 0.0d0
               if ( qnum_l_to_probe == 1 ) then
                  z1 = EigenWfns_MatLS_L1(jb,ib)
               else if ( qnum_l_to_probe == 2 ) then
                  z1 = EigenWfns_MatLS_L2(jb,ib)
               else if ( qnum_l_to_probe == 3 ) then
                  z1 = EigenWfns_MatLS_L3(jb,ib)
               else
                  if ( ib == jb ) z1 = 1.0d0
               endif

               fsr_core_states(ib,im1,ik+is-1) &
                    &     = fsr_core_states(ib,im1,ik+is-1) +real(z1)
               fsi_core_states(ib,im1,ik+is-1) &
                    &     = fsi_core_states(ib,im1,ik+is-1) +aimag(z1)
!!!!!!                    &     = fsi_core_states(ib,im1,ik+is-1) -aimag(z1)
            End do
         End Do
      End do

    end subroutine set_core_states_with_phase_A

  end subroutine m_CLS_set_wfn_core_states

  subroutine read_num_core_ae_wfns( nfp, nums, it )
    integer, intent(in) :: nfp, it
    integer, intent(out) :: nums

    integer :: length, ierr
    character(30) :: search_key

    nums = 0

    search_key = "CORE STATES";  length = len(search_key)
    call read_size_of_array_from_pp(nfp, nums, length, search_key, ierr)

    if ( ierr /= 0 ) then
       write(nfout,*) '----------------------'
       write(nfout,'(A,I2)') '!!! Keyword CORE STATES is not found in the PP', it
       write(nfout,*) '----------------------'
    endif

  end subroutine read_num_core_ae_wfns

  subroutine read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh, &
       &                             qnum_n, qnum_l, psir_core, ene_level, focc )
    implicit none

    integer, intent(in) :: num_core_ae_wfns
    integer, intent(in) :: nmesh, nfp
    integer, intent(out) :: qnum_l(num_core_ae_wfns)
    integer, intent(out) :: qnum_n(num_core_ae_wfns)
    real(kind=8), intent(out) :: psir_core(nmesh, num_core_ae_wfns), &
         &                       ene_level(num_core_ae_wfns), focc(num_core_ae_wfns)

    integer :: i, n1, l1, k
    real(kind=DP) :: ene1, f1

    Do i=1, num_core_ae_wfns
       read(nfp,*) n1, l1, ene1, f1
       read(nfp,*) (psir_core(k,i),k=1,nmesh)
       !
       qnum_n(i) = n1;  qnum_l(i) = l1;  ene_level(i) = ene1;   focc(i) = f1
    End Do
  end subroutine read_data_core_ae_wfns

  subroutine read_data_core_energy_contrib( nfp, it, ekin_core, eion_core, ehart_core )
    implicit none

    integer, intent(in) :: nfp, it
    real(kind=DP), intent(out) :: ekin_core, eion_core, ehart_core

    integer :: ifound, ier
    character(30) :: search_tag
    character(30) :: line1

    ierr = 0;
    ekin_core = 0.0d0;  eion_core = 0.0d0;  ehart_core = 0.0d0

    search_tag = "CORE ENERGY CONTRIB"

    Do while (.true.)
       read(nfp,'(a30)',end=10) line1
       ifound = index( line1, search_tag )
       if ( ifound /= 0 ) goto 20
    End do

10  ierr = 1;
    write(nfout,*) '----------------------'
    write(nfout,'(A,I2)') '!!! Keyword CORE ENERGY CONTRIB is not found in the PP', it
    write(nfout,*) '----------------------'
    return

20  continue
    read(nfp,*) ekin_core,  line1
    read(nfp,*) eion_core,  line1
    read(nfp,*) ehart_core, line1

  end subroutine read_data_core_energy_contrib

! == KT_add === 2014/08/08
  subroutine read_data_core_level_splitting( nfp, it, num_core_ae_wfns, &
       &                                     qnum_n, qnum_l, level_splitting )
    implicit none

    integer, intent(in) :: nfp, it, num_core_ae_wfns
    integer, intent(in) :: qnum_n(num_core_ae_wfns)
    integer, intent(in) :: qnum_l(num_core_ae_wfns)
    real(kind=DP), intent(out) :: level_splitting(num_core_ae_wfns)

    integer :: ifound, ier
    integer :: count, n1, l1, i, j
    real(kind=DP) :: c1
    character(30) :: search_tag
    character(30) :: line1

    ierr = 0;
    level_splitting = 0.0d0

    search_tag = "SOC-CORE"

    Do while (.true.)
       read(nfp,'(a30)',end=10) line1
       ifound = index( line1, trim(search_tag) )
       if ( ifound /= 0 ) goto 20
    End do

10  ierr = 1;
    write(nfout,*) '----------------------'
    write(nfout,'(A,I2)') '!!! Keyword SOC-CORE is not found in the PP', it
    write(nfout,*) '----------------------'
    return

20  continue

    read(nfp,*) count
    Do i=1, count
       read(nfp,*) n1, l1, c1
       Do j=1, num_core_ae_wfns
          if ( qnum_n(j) == n1 .and. qnum_l(j) == l1 ) then
             level_splitting(j) = c1 *( 2.0d0 *l1 +1.0d0 ) /2.0d0
             exit
          endif
       End Do
    End Do
  end subroutine read_data_core_level_splitting
! ========= 2014/08/08

  subroutine read_num_dipole_dxyz_core2val( nfp, nums, it )
    integer, intent(in) :: nfp, it
    integer, intent(out) :: nums

    integer :: length, ierr
    character(30) :: search_key

    nums = 0

    search_key = "DIPOLE-CORE-TO-VALENCE";  length = len(search_key)

    call read_size_of_array_from_pp(nfp, nums, length, search_key, ierr)
    if ( ierr /= 0 ) then
       write(nfout,*) '----------------------'
       write(nfout,'(A,I2)') '!!! Keyword DIPOLE-CORE-TO-VALENCE not found in PP', it
       write(nfout,*) '----------------------'
    endif

  end subroutine read_num_dipole_dxyz_core2val

! === KT_add === 2014/08/11
  subroutine read_dipole_dxyz_core2val_A( nfp, ndata, dipole_n1, dipole_t2, &
       &                                  phase_ylm1, phase_ylm2, data_dipole, it )
    integer, intent(in) :: nfp, ndata, it
    integer, intent(out) :: dipole_n1(ndata)
    integer, intent(out) :: dipole_t2(ndata)
    integer, intent(out) :: phase_ylm1(ndata), phase_ylm2(ndata)
    real(kind=DP), intent(out) :: data_dipole(:,:)

    integer :: i, j1
    integer :: n1, l1, t1, m1, n2, l2, t2, m2

    do i=1, ndata
       read(nfp,53) n1, l1, t1, m1, n2, l2, t2, m2, &
            &       data_dipole(i,1), data_dipole(i,2), data_dipole(i,3), &
            &       phase_ylm1(i), phase_ylm2(i)

       dipole_n1(i) = n1;   dipole_t2(i) = t2

       Do j1=1, num_orb_index_data(it)
          if ( n2 == qnum_n_orb_index(it,j1) &
               &    .and. l2 == qnum_l_orb_index(it,j1) &
               &    .and. t2 == qnum_t_orb_index(it,j1) ) then
             dipole_t2(i) = qnum_tau_orb_index(it,j1)
             exit
          end if
       End Do
    end do

53  format(1x,8i3,3e18.10,2i3)

  end subroutine read_dipole_dxyz_core2val_A
! ==================== 2014/08/11

  subroutine read_dipole_dxyz_core2val_B( nfp, ndata, dipole_n1, dipole_t2, &
       &                                  phase_ylm1, phase_ylm2, data_dipole )
    integer, intent(in) :: nfp, ndata
    integer, intent(out) :: dipole_n1(ndata)
    integer, intent(out) :: dipole_t2(ndata)
    integer, intent(out) :: phase_ylm1(ndata), phase_ylm2(ndata)
    real(kind=DP), intent(out) :: data_dipole(:,:)

    integer, allocatable :: dipole_n2(:)
    integer, allocatable :: dipole_l1(:), dipole_l2(:)
    integer, allocatable :: dipole_m1(:), dipole_m2(:)
    integer, allocatable :: dipole_t1(:)

    integer, allocatable :: igeta(:,:)

    integer :: lmax, nmax, i, lval, count
    integer :: n2, l2, t2, old_n2, old_t2

    lmax = 3         ! s:0, p:1,  d:2,  f:3
    nmax = 8         ! 8s, 8p, 8d .. etc.

    allocate( dipole_n2(ndata) )
    allocate( dipole_l1(ndata) );  allocate( dipole_l2(ndata) )
    allocate( dipole_m1(ndata) );  allocate( dipole_m2(ndata) )
    allocate( dipole_t1(ndata) )

    allocate( igeta(nmax,0:lmax) );  igeta = -100

    do i=1, ndata
       read(nfp,53) dipole_n1(i), dipole_l1(i), dipole_t1(i), dipole_m1(i), &
            &       dipole_n2(i), dipole_l2(i), dipole_t2(i), dipole_m2(i), &
            &       data_dipole(i,1), data_dipole(i,2), data_dipole(i,3), &
            &       phase_ylm1(i), phase_ylm2(i)
    End do

    Do lval=0, lmax
       count = 0

       Do i=1, ndata
          if ( dipole_l2(i) == lval ) then
             n2 = dipole_n2(i);  t2 = dipole_t2(i)

             if ( count==0 ) then
                count = count +1
                if ( igeta(n2,lval) == -100 ) igeta(n2,lval) = 0
             else
                if ( old_n2 /=n2 .and. igeta(n2,lval)== -100 ) igeta(n2,lval)=old_t2
             endif
             old_n2 = n2;  old_t2 = t2
          endif
       End do
    End do

    Do i=1, ndata
       n2 = dipole_n2(i); l2 = dipole_l2(i);  t2 = dipole_t2(i)

!       dipole_t1(i) = dipole_t1(i) +igeta(n1,l1)
       dipole_t2(i) = dipole_t2(i) +igeta(n2,l2)
    end do

53  format(1x,8i3,3e18.10,2i3)

    deallocate( dipole_n2 )
    deallocate( dipole_l1 );  deallocate( dipole_l2 )
    deallocate( dipole_m1 );  deallocate( dipole_m2 )
    deallocate( dipole_t1 )

    deallocate( igeta )

  end subroutine read_dipole_dxyz_core2val_B

  subroutine m_CLS_find_ptrans_indx_core2val( n1, nspher1, nspher2, tau2, index )
!   find core-repair term with dipole_n1 = n1 (quantum number n of core states),
!                              dipole_tau2 = tau2 (valence orb)
!                              phase_ylm = nspher (combination of (l,m) )

    integer, intent(in)  :: n1, nspher1, nspher2, tau2
    integer, intent(out) :: index

    integer              :: i

    index = 0

    do i=1, num_dipole_dxyz_core2val
       if( dipole_dxyz_core2val_ylm1(i) == nspher1 &
            & .and. dipole_dxyz_core2val_ylm2(i) == nspher2 ) then
          if ( dipole_dxyz_core2val_n1(i) == n1 &
               &  .and. dipole_dxyz_core2val_t2(i) == tau2 ) then
             index = i
             exit
          end if
       end if
    end do

  end subroutine m_CLS_find_ptrans_indx_core2val

  subroutine read_size_of_array_from_pp( nfp, size_of_array, length, search_tag, ierr )
    implicit none

    integer, intent(in) :: nfp, length
    integer, intent(out) :: size_of_array, ierr
    character(length), intent(in) :: search_tag

    integer :: ifound
    character(30) :: line1

    size_of_array = 0;  ierr = 0

    Do while (.true.)
       read(nfp,'(a30)',end=10) line1
       ifound = index( line1, search_tag )
       if ( ifound /= 0 ) goto 20
    End do

10  ierr = 1; return

20  read(nfp,*) size_of_array

  end subroutine read_size_of_array_from_pp

! ----------------------------------------------
!!
!!! The following routine is used only in the postprocess of PHASE
!!!           ( before calculating spectrum )
!!
! ----------------------------------------------
  subroutine m_CLS_calc_core_energy
    use m_Control_Parameters, only : atom_with_corehole, qnum_n_corehole, &
         &                           qnum_l_corehole
    use m_PAW_XC_Potential,   only : m_PAW_XC_allocation, m_PAW_XC_deallocation, &
         &                            m_PAW_XC_cal_potential_sphex2, &
         &                            m_PAW_XC_cal_potential
    use m_PAW_ChargeDensity,    only : calcSphericalHarmonicsExpansion, &
         &                             calcGaussLegendreIntegration
    use m_Const_Parameters, only : VXC_AND_EXC
    use m_PseudoPotential,  only : flg_symmtry

    integer :: it, ia, i, nfp
    real(kind=DP) :: etot_core_sum
    real(kind=DP) :: ekin_core, ehart_core, eion_core, level_splitting_corehole
    real(kind=DP), allocatable :: etot_core(:), level_splitting(:)

     if ( flg_paw .and. eval_core_level_splitting == ByPawPot ) then
        call m_PAW_XC_allocation(nfout)
        if(calcSphericalHarmonicsExpansion) then
           call m_PAW_XC_cal_potential_sphex2(nfout,VXC_AND_EXC)
        end if
        if(calcGaussLegendreIntegration)then
           call m_PAW_XC_cal_potential(nfout,VXC_AND_EXC,flg_symmtry)
        endif
     endif

!    call m_Files_open_ps_files(ivan,iatomn,ntyp,ierr)
!    if (ierr/=0) call mpi_stop(nfout)

    call m_Files_open_core_energy_file

    allocate( etot_core(ntyp) ); etot_core = 0.0d0;
    etot_core_sum = 0.0d0
    level_splitting_corehole = 0.0d0

    if ( mype == 0 ) call print_header

    Do it=1, ntyp
       call m_Files_open_ps_file(ivan,iatomn,ntyp,it,ierr)
       if (ierr/=0) call mpi_stop(nfout)
       nfp = nfpot(it)

       if ( mype == 0 ) call read_num_core_ae_wfns( nfp, num_core_ae_wfns, it )
       if ( npes > 1 ) then
          call mpi_bcast(  num_core_ae_wfns, 1, mpi_integer, 0, mpi_comm_world, ierr )
       endif

       if ( num_core_ae_wfns == 0 ) cycle

       allocate( qnum_n_core_ae_wfns( num_core_ae_wfns ) )
       allocate( qnum_l_core_ae_wfns( num_core_ae_wfns ) )
       allocate( psir_core_ae_wfns( nmesh(it),  num_core_ae_wfns ) )
       allocate( enelevel_core_ae_wfns( num_core_ae_wfns ) )
       allocate( focc_core_ae_wfns( num_core_ae_wfns ) )
       allocate( level_splitting( num_core_ae_wfns ) )
       level_splitting = 0.0d0

       if ( mype == 0 ) then
          call read_data_core_ae_wfns( nfp, num_core_ae_wfns, nmesh(it), &
               &                       qnum_n_core_ae_wfns, qnum_l_core_ae_wfns, &
               &                       psir_core_ae_wfns, &
               &                       enelevel_core_ae_wfns, focc_core_ae_wfns )
          call read_data_core_energy_contrib( nfp, it, ekin_core, eion_core, ehart_core )
          etot_core(it) = ekin_core +eion_core +ehart_core
       endif

       if ( npes > 1 ) then
          call mpi_bcast( qnum_n_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_integer, 0, mpi_comm_world, ierr )
          call mpi_bcast( qnum_l_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_integer, 0, mpi_comm_world, ierr )
          call mpi_bcast( enelevel_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( focc_core_ae_wfns, num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( psir_core_ae_wfns, nmesh(it)*num_core_ae_wfns, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
          call mpi_bcast( etot_core(it), 1, &
               &          mpi_double_precision, 0, mpi_comm_world, ierr )
       endif

       write(nfout,'(A,I5,F25.15)') " ! ", it, etot_core(it)

! === KT_add ===== 2014/08/08
       if ( eval_core_level_splitting == ReadFromPP ) then
          if ( mype == 0 ) then
             call read_data_core_level_splitting( nfp, it, num_core_ae_wfns, &
                  &                               qnum_n_core_ae_wfns, &
                  &                               qnum_l_core_ae_wfns, &
                  &                               level_splitting )
          endif
          if ( npes > 1 ) then
             call mpi_bcast( level_splitting, num_core_ae_wfns, &
                  &          mpi_double_precision, 0, mpi_comm_world, ierr )
          endif
       endif
! =============== 2014/08/08

       Do ia=1, natm
          if ( ityp(ia) /= it ) cycle

          etot_core_sum = etot_core_sum +etot_core(it)

          if ( flg_paw .and. eval_core_level_splitting == ByPawPot ) then
             call m_SO_calc_corelevel_splitting( ia, it, &
                  &                              num_core_ae_wfns, &
                  &                              qnum_l_core_ae_wfns, &
                  &                              nmesh(it), psir_core_ae_wfns, &
                  &                              level_splitting )
          endif

          if ( ia == atom_with_corehole ) then
             Do i=1, num_core_ae_wfns
                if ( qnum_n_core_ae_wfns(i) == qnum_n_corehole &
                     &  .and. qnum_l_core_ae_wfns(i) == qnum_l_corehole ) then
                   if ( focc_core_ae_wfns(i) > 2*qnum_l_corehole ) then
                      level_splitting_corehole = level_splitting(i)
                   endif
                endif
             End Do
          endif

          if ( mype == 0 ) call print_body
       End do

       deallocate( qnum_n_core_ae_wfns ); deallocate( qnum_l_core_ae_wfns )
       deallocate( psir_core_ae_wfns );   deallocate( enelevel_core_ae_wfns )
       deallocate( focc_core_ae_wfns )
       deallocate( level_splitting )
       call m_Files_close_ps_file(it)
    End Do

!    call m_Files_close_ps_files
    deallocate( etot_core )

    if ( mype == 0 ) then
       call print_tail
       call add_soc_corehole
    endif

    call m_Files_close_core_energy_file


    if ( flg_paw .and. eval_core_level_splitting == ByPawPot ) then
       call m_PAW_XC_deallocation(nfout)
    endif

  contains

    subroutine print_header
      write(nfout,*) '**********************************************'
      write(nfout,*) "! ** Total Energy of Core electrons for atom types "
      write(nfout,*) "!    it        energy "

      write(nfcore_energy_out,'(A)') "# Spin Orbit splitting of core orbitals"
      write(nfcore_energy_out,'(A)') " it    ia     n     l            value"

    end subroutine print_header

    subroutine print_body
      integer :: i

      Do i=1, num_core_ae_wfns
         if ( qnum_l_core_ae_wfns(i) == 0 ) cycle
         write(nfcore_energy_out,'(I3,I6,I6,I6,F25.15)') &
              &   it, ia, qnum_n_core_ae_wfns(i), qnum_l_core_ae_wfns(i), &
              &           level_splitting(i)
      End Do
    end subroutine print_body

    subroutine print_tail
      write(nfout,*) '!'
      write(nfout,*) '!  Etotal (Core+Valence) = ', etot_core_sum +etotal
      write(nfout,*) '**********************************************'

      write(nfcore_energy_out,'(I3,I6,I6,I6,F25.15)') 0, 0, 0, 0, 0.0d0
      write(nfcore_energy_out,*)
      write(nfcore_energy_out,'(A)') "# Etotal (Core+Valence)"
      write(nfcore_energy_out,*) etot_core_sum +etotal

    end subroutine print_tail

    subroutine add_soc_corehole
      real(kind=DP) :: c1, c2, e1, e2

      if ( level_splitting_corehole <= 0.0 ) return

      c1 = etot_core_sum +etotal
      c2 = level_splitting_corehole /( 2*qnum_l_corehole +1 )
      e1 = c2 *qnum_l_corehole
      e2 = c2 *( qnum_l_corehole +1 )

      write(nfcore_energy_out,*)
      write(nfcore_energy_out,'(A)') "# Etotal (Core+Valence+Soc_corehole)"
      write(nfcore_energy_out,'(A,I1,A,10X,F20.10)') &
           &            "# J = ", 2*qnum_l_corehole +1, "/2:", c1 -e1
      write(nfcore_energy_out,'(A,I1,A,10X,F20.10)') &
           &            "# J = ", 2*qnum_l_corehole -1, "/2:", c1 +e2
    end subroutine add_soc_corehole

  end subroutine m_CLS_calc_core_energy

end module m_CoreLevel_Spectrum
