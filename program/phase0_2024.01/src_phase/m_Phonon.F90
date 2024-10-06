!=======================================================================
!
!  SOFTWARE NAME : PHASE ($Revision: 603 $)
!
!  MODULE: m_Phonon
!
!  AUTHOR(S): T. Yamamoto   May/01/2004
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
!
module m_Phonon
  ! $Id: m_Phonon.F90 603 2020-04-07 03:25:30Z jkoga $
  use m_Const_parameters,   only : DP, ON, OFF, NOCONV, FMAXVALLEN, PAI, PAI2, PAI4 &
       &     , PHONON_GAMMA, PHONON_BAND, PHONON_DOS &
       &     , LOWER, CARTS, UNIT_PIEZO_CONST &
       &     , CRDTYP, FILE, FERMI_DIRAC, Hartree, AMU, CONST_EV, PLANCK, speed_of_light
  use m_Control_Parameters, only : sw_phonon, sw_calc_force, sw_calc_force_all &
       &     , sw_vibrational_modes, sw_lo_to_splitting &
       &     , sw_lattice_dielectric_tensor &
       &     , sw_dielectric_function &
       &     , sw_correct_force_constants &
       &     , sw_polynomial_fit, norder &
       &     , sw_int_strain_piezo_tensor, sw_phonon_oneshot &
       &     , phonon_method, num_phonon_calc_mode &
       &     , ipriphonon, printable, way_of_smearing, width, width_tetra
  use m_Crystal_Structure,  only : altv,rltv,altv_super,rltv_super,nlpnt,lpnt,nbztyp_spg
  use m_Ionic_System,       only : num_force_data, phonon_atom &
       &, phonon_displacement &
       &,displaced_atom, displacement &
       &,input_coordinate_system,u &
       &,natm,ionic_mass,speciesname,ityp,ntyp &
       &,istart_phonon,iend_phonon &
       &,natm_prim,natm_super &
       &,pos,cps,phonon_iteration, natm_mobile, force_was_read, atom_key
  use m_BP_Properties,      only : zeff
  use m_Files,              only : nfout,nfkpoint,nfmatbp,m_Files_open_kpoint_files, &
       &                           nfqpoint, m_Files_open_qpoint_files
  use m_Parallelization,    only : mype,npes,ierr, MPI_CommGroup

  ! === KT_add === 13.1R
  use m_Control_Parameters,  only : sw_phonon_with_epsilon
  use m_Raman,               only : m_Raman_read_nfinp, m_Raman_initialize, &
       &                            m_Raman_print_param, m_Raman_calc_susceptibility, &
       &                            dielectric_read_from_file => dielectric
  ! ============== 13.1R

  use m_Kpoints,             only : itrs
  use m_Crystal_Structure,  only : rltv_refcell, altv_refcell
  use mpi

  implicit none
!  include 'mpif.h'

  real(kind=DP),allocatable,dimension(:,:,:) :: force_data ! dim(natm_super,3,num_force_data)
  real(kind=DP),allocatable,dimension(:,:,:) :: forces ! dim(natm_super*3,natm_prim*3,norder*2)
  real(kind=DP),allocatable,dimension(:) :: force0 ! dim(natm_super*3)
  real(kind=DP) :: dielectric(3,3) = 0.d0
  real(kind=DP) :: kvec_dir(3) = 0.d0
  character(len=3), allocatable, dimension(:) :: point_group ! dim(nqvec)
  character(len=3) :: pg_gamma = ''
  integer :: nlat(3) ! lattice size
  real(kind=DP),allocatable,dimension(:,:) :: qvec ! dim(nqvec,3)
  real(kind=DP),allocatable,dimension(:,:) :: qvin ! dim(nqvec,3)
  real(kind=DP),allocatable,dimension(:)   :: wght ! dim(nqvec)
  real(kind=DP),allocatable,dimension(:,:) :: qdir ! dim(nqvec,3)
  integer :: nsp = 1
  real(kind=DP),allocatable,dimension(:,:) :: spcpnt ! dim(nsp,3)
  character(len=3), allocatable, dimension(:)  :: pg_spcpnt ! dim(nsp)
  integer :: nsl = 1
  integer,allocatable,dimension(:,:) :: symln  ! dim(nsl,2)
  integer,allocatable,dimension(:) :: ndiv     ! dim(nsl)
  character(len=3), allocatable, dimension(:)  :: pg_symln ! dim(nsl)
  logical, allocatable, dimension(:,:) :: fksym ! dim(48,nqvec)
  real(kind=DP) :: rho_loto = 1.4d0
  ! strain force coupling const.
  real(kind=DP),allocatable         :: strfrc(:,:,:)   !d(3,natm,6)
  real(kind=DP) :: piezoelectric(3,6) = 0.d0

  ! energy range for dielectric function
  real(kind=DP), private :: max_energy = 0.01d0
  real(kind=DP), private :: min_energy = 0.d0
  integer, private :: num_division = 100
  real(kind=DP), private :: damping_factor = 1.d-2

  ! Phonon DOS
  integer, private :: phonon_dos_mesh(3)
  integer                           :: np0,np1,np2
  integer,allocatable, dimension(:) :: ip20,ip10,ip02,ip12
  integer,allocatable, dimension(:) :: ip01,ip21,iu21,iv21,iwt,ip2cub
  integer, dimension(3)             :: nxyz_tetra
  real(kind=DP)                     :: DeltaE_dos = 1.d-5


  ! calculated
  integer :: nmodes ! = 3*natm
  integer :: nqvec = 1
  integer :: qvec_from_file = ON
  real(kind=DP), allocatable,dimension(:,:,:) :: dynmat  ! dim(3*natm,3*natm,nlpnt)
  complex(kind=DP), allocatable,dimension(:,:,:) :: dynmatq ! dim(3*natm,3*natm,nqvec)
  real(kind=DP), allocatable,dimension(:,:)   :: omega   ! dim(3*natm,nqvec)
  complex(kind=DP), allocatable,dimension(:,:,:) :: modes   ! dim(3*natm,3*natm,nqvec)
  character(len=3), allocatable,dimension(:,:) :: mode_irr_rep ! dim(3*natm,nqvec)
  character(len=4), allocatable,dimension(:,:) :: mode_active  ! dim(3*natm,nqvec)
  real(kind=DP), allocatable,dimension(:)   :: msqinv     ! dim(3*natm)
  real(kind=DP), allocatable,dimension(:,:) :: zeff_mode  ! dim(3,nmodes)
  real(kind=DP) :: static_dielectric(3,3) ! static dielectric tensor
  real(kind=DP) :: lattice_dielectric(3,3) ! lattice dielectric tensor
  real(kind=DP), allocatable,dimension(:,:) :: strphcpl  ! dim(6,nmodes)

! real(DP), parameter :: ev2cminv = 1.d0/1.23984185d-4 ! 1 cm^-1 = 1.233984185e-4 eV
  real(DP), parameter :: ev2cminv = CONST_EV / PLANCK / speed_of_light * 1.0d-2

  ! Ewald sum
  real(kind=DP) :: deteps ! Det(epsilon)
  real(kind=DP), dimension(3,3) :: epsinv ! (epsilon)^-1

  integer :: sw_read_forces_pre = ON

  integer, allocatable, dimension(:) :: mode2atom
  integer, allocatable, dimension(:) :: atm2mobileatm
  ! Tags
  character(len("phonon")),private,parameter :: tag_phonon = "phonon"
  character(len("sw_phonon")),private,parameter :: tag_sw_phonon = "sw_phonon"
  character(len("sw_calc_force")),private,parameter :: tag_sw_calc_force = "sw_calc_force"
  character(len("sw_calc_force_all")),private,parameter :: tag_sw_calc_force_all = "sw_calc_force_all"
  character(len("sw_polynomial_fit")),private,parameter :: tag_sw_polynomial_fit = "sw_polynomial_fit"
  character(len("sw_phonon_oneshot")),private,parameter :: tag_sw_phonon_oneshot = "sw_phonon_oneshot"
  character(len("num_phonon_calc_mode")),private,parameter :: tag_num_phonon_calc_mode = "num_phonon_calc_mode"
  character(len("norder")),private,parameter :: tag_norder = "norder"
  character(len("displacement")),private,parameter :: tag_displacement = "displacement"
  character(len("force_calc")),private,parameter :: tag_force_calc = "force_calc"
  character(len("start")),private,parameter :: tag_start = "start"
  character(len("end")),private,parameter :: tag_end = "end"

  character(len("sw_vibrational_modes")),private,parameter :: tag_sw_vibrational_modes = "sw_vibrational_modes"
  character(len("sw_correct_force_constants")),private,parameter :: tag_sw_correct_force_constants = "sw_correct_force_constants"
  character(len("point_group")),private,parameter :: tag_point_group = "point_group"
  character(len("sw_lo_to_splitting")),private,parameter :: tag_sw_lo_to_splitting = "sw_lo_to_splitting"
  character(len("electronic_dielectric_constant")),private,parameter :: &
       & tag_electronic_dielectric_cnst = "electronic_dielectric_constant"
  character(len("exx")),private,parameter :: tag_exx = "exx"
  character(len("eyy")),private,parameter :: tag_eyy = "eyy"
  character(len("ezz")),private,parameter :: tag_ezz = "ezz"
  character(len("exy")),private,parameter :: tag_exy = "exy"
  character(len("eyz")),private,parameter :: tag_eyz = "eyz"
  character(len("ezx")),private,parameter :: tag_ezx = "ezx"
  character(len("k_vector")),private,parameter :: tag_k_vector = "k_vector"
  character(len("kx")),private,parameter :: tag_kx = "kx"
  character(len("ky")),private,parameter :: tag_ky = "ky"
  character(len("kz")),private,parameter :: tag_kz = "kz"

  character(len("sw_lattice_dielectric_tensor")),private,parameter :: &
       & tag_sw_lattice_dielec_tensor = "sw_lattice_dielectric_tensor"
  character(len("sw_internal_strain_piezoelectric_tensor")),private,parameter :: &
       & tag_sw_int_strain_piezo_tensor = "sw_internal_strain_piezoelectric_tensor"

  character(len("sw_dielectric_function")),private,parameter :: tag_sw_dielectric_function = "sw_dielectric_function"
  character(len("energy_range")),private,parameter :: tag_energy_range = "energy_range"
  character(len("min_energy")),private,parameter :: tag_min_energy = "min_energy"
  character(len("max_energy")),private,parameter :: tag_max_energy = "max_energy"
  character(len("division_number")),private,parameter :: tag_division_number = "division_number"
  character(len("damping_factor")),private,parameter :: tag_damping_factor = "damping_factor"


  character(len("lattice")),private,parameter :: tag_lattice = "lattice"
  character(len("l1")),private,parameter :: tag_l1 = "l1"
  character(len("l2")),private,parameter :: tag_l2 = "l2"
  character(len("l3")),private,parameter :: tag_l3 = "l3"

  character(len("no")),private,parameter :: tag_no = "no"
  character(len("id")),private,parameter :: tag_id = "id"
  character(len("special_points")),private,parameter :: tag_special_points = "special_points"
  character(len("symmetry_lines")),private,parameter :: tag_symmetry_lines = "symmetry_lines"
  character(len("k1")),private,parameter :: tag_k1 = "k1"
  character(len("k2")),private,parameter :: tag_k2 = "k2"
  character(len("k3")),private,parameter :: tag_k3 = "k3"
  character(len("num_division")),private,parameter :: tag_num_division = "num_division"

  character(len("method")),private,parameter :: tag_method = "method"
  character(len("calc_type")),private,parameter :: tag_calc_type = "calc_type"
  character(len("gamma")),private,parameter  :: tag_gamma  = "gamma"
  character(len("zone_center")),private,parameter :: tag_zone_center = "zone_center"
  character(len("band")),private,parameter   :: tag_band   = "band"
  character(len("dispersion")),private,parameter :: tag_dispersion = "dispersion"
  character(len("dos")),private,parameter    :: tag_dos    = "dos"

  character(len("mesh")),private,parameter   :: tag_mesh   = "mesh"
  character(len("nx")),private,parameter     :: tag_nx     = "nx"
  character(len("ny")),private,parameter     :: tag_ny     = "ny"
  character(len("nz")),private,parameter     :: tag_nz     = "nz"
  character(len("deltae")),private,parameter :: tag_deltae = "deltae"

  character(len("rho")),private,parameter :: tag_rho = "rho"

  character(len("qvec_from_file")),private, parameter :: tag_qvec_from_file = "qvec_from_file"
  character(len("sw_read_forces_pre")),private, parameter :: tag_sw_read_forces_pre = "sw_read_forces_pre"

  character(len("use_qpoint_data_file")),private, parameter &
       &   :: tag_use_qpoint_data_file = "use_qpoint_data_file"
  integer :: use_qpoint_data_file = OFF

! -------------------------
! phonon band unfolding
!
  character(len("phonon_band_unfolding")),private, parameter :: &
       &      tag_phonon_band_unfolding = "phonon_band_unfolding"
  character(len("sw_phonon_band_unfolding")),private, parameter :: &
       &      tag_sw_phonon_band_unfolding = "sw_phonon_band_unfolding"

  character(len("ngx")),private,parameter     :: tag_ngx     = "ngx"
  character(len("ngy")),private,parameter     :: tag_ngy     = "ngy"
  character(len("ngz")),private,parameter     :: tag_ngz     = "ngz"

  character(len("tolerance_gvec_matching")), private, parameter :: &
       &                     tag_tol_ph_gvec_matching = "tolerance_gvec_matching"

  integer :: sw_phonon_band_unfolding = OFF
  integer :: nmax_G_phonon(3) = (/4,4,4/)
  real(kind=DP) :: tolerance_ph_Gvec_matching = 1.0D-3

  integer, allocatable :: phonon_GVec_on_refcell(:,:,:)
  real(kind=DP),allocatable :: unfolding_weight(:,:)
  real(kind=DP),allocatable :: qvin_refcell(:,:)     ! dim(nqvec,3)

contains

  subroutine m_Phonon_rd_param(nfout)

    integer, intent(in) :: nfout

    integer :: iret, f_selectBlock, f_getIntValue, f_getRealValue, f_getStringValue
    integer :: f_selectParentBlock, f_selectTop
    real(kind=DP) :: dret
    character(len=FMAXVALLEN) :: cret
    real(kind=DP) :: dnorm

    iret = f_selectTop()
    if( f_selectBlock( tag_phonon) == 0) then
       if( f_getIntValue( tag_sw_phonon, iret) == 0) sw_phonon = iret
       if(sw_phonon == ON) then
          sw_calc_force = ON
          if( f_getStringValue( tag_method, cret,LOWER) == 0) then
             call set_phonon_method(cret) ! -> phonon_method
          else if( f_getStringValue( tag_calc_type, cret,LOWER) == 0) then
             call set_phonon_method(cret) ! -> phonon_method
          end if

          if( f_getIntValue( tag_sw_calc_force, iret) == 0) sw_calc_force = iret
          if( f_getIntValue( tag_sw_calc_force_all, iret) == 0) sw_calc_force_all = iret
          if( f_getIntValue( tag_sw_polynomial_fit, iret) == 0) sw_polynomial_fit = iret
          if( f_getIntValue( tag_norder, iret) == 0) norder = iret
          if( f_getRealValue( tag_displacement, dret, 'bohr') == 0) then
             u = dret
          else
             u = 0.1d0 ! default value
          end if
          if( f_getRealValue( tag_damping_factor, dret, '') == 0) damping_factor = dret
          damping_factor = abs(damping_factor)
          if( f_getIntValue( tag_sw_vibrational_modes, iret) == 0) sw_vibrational_modes = iret
          if( f_getIntValue( tag_sw_correct_force_constants, iret) == 0) sw_correct_force_constants = iret
          if( f_getStringValue( tag_point_group, cret, NOCONV) == 0) then
             pg_gamma = cret(1:3)
          else
             pg_gamma = '' ! default value
          end if

          if(f_getIntValue(tag_sw_read_forces_pre,iret)==0) then
             sw_read_forces_pre = iret
          endif

          if( f_selectBlock( tag_phonon_band_unfolding ) == 0) then
             if( f_getIntValue( tag_sw_phonon_band_unfolding, iret)==0) &
                  &      sw_phonon_band_unfolding = iret

             if( f_getRealValue( tag_tol_ph_gvec_matching, dret, "" ) == 0) then
                tolerance_ph_Gvec_matching = dret
             endif

             if( f_getIntValue( tag_ngx, iret) == 0) nmax_G_phonon(1) = iret
             if( f_getIntValue( tag_ngy, iret) == 0) nmax_G_phonon(2) = iret
             if( f_getIntValue( tag_ngz, iret) == 0) nmax_G_phonon(3) = iret
             iret = f_selectParentBlock()
          end if

          if( f_getIntValue( tag_sw_phonon_oneshot, iret)==0) sw_phonon_oneshot = iret
          if( f_getIntValue( tag_num_phonon_calc_mode,iret)==0) num_phonon_calc_mode = iret

          ! ==== KT_add ==== 13.1R
          call m_Raman_read_nfinp
          ! =============== 13.1R

          if( f_getIntValue( tag_sw_lattice_dielec_tensor, iret) == 0) &
               &  sw_lattice_dielectric_tensor = iret
          if( f_getIntValue( tag_sw_int_strain_piezo_tensor, iret) == 0) &
               &  sw_int_strain_piezo_tensor = iret

          if( f_getIntValue( tag_sw_dielectric_function, iret) == 0) sw_dielectric_function = iret

          if( f_getIntValue( tag_sw_lo_to_splitting, iret) == 0) sw_lo_to_splitting = iret
          if(sw_lo_to_splitting == ON) then
             if( f_getRealValue( tag_rho, dret, '') == 0) rho_loto = dret
          end if
          if( f_selectBlock( tag_force_calc) == 0) then
             if( f_getIntValue( tag_start, iret) == 0) istart_phonon = iret
             if( f_getIntValue( tag_end, iret) == 0) iend_phonon = iret
             iret = f_selectParentBlock()
          end if
          if( f_selectBlock( tag_electronic_dielectric_cnst) == 0) then
             if( f_getRealValue( tag_exx, dret, '') == 0) dielectric(1,1) = dret
             if( f_getRealValue( tag_eyy, dret, '') == 0) dielectric(2,2) = dret
             if( f_getRealValue( tag_ezz, dret, '') == 0) dielectric(3,3) = dret
             if( f_getRealValue( tag_exy, dret, '') == 0) dielectric(1,2) = dret
             dielectric(2,1) = dielectric(1,2)
             if( f_getRealValue( tag_eyz, dret, '') == 0) dielectric(2,3) = dret
             dielectric(3,2) = dielectric(2,3)
             if( f_getRealValue( tag_ezx, dret, '') == 0) dielectric(3,1) = dret
             dielectric(1,3) = dielectric(3,1)
             iret = f_selectParentBlock()
          end if
          if( f_selectBlock( tag_energy_range ) == 0) then
             if( f_getRealValue( tag_min_energy, dret, 'hartree') == 0) min_energy = dret
             if( f_getRealValue( tag_max_energy, dret, 'hartree') == 0) max_energy = dret
             if( f_getIntValue( tag_division_number, iret) == 0) num_division = iret
             iret = f_selectParentBlock()
          end if

          if( f_getIntValue(tag_qvec_from_file,iret)==0) qvec_from_file = iret

          ! ==== KT_add ====== 13.1R
          if( phonon_method == PHONON_GAMMA ) qvec_from_file = OFF
          ! ================== 13.1R

          if( phonon_method == PHONON_BAND .and. qvec_from_file== ON ) then
             if( f_getIntValue(tag_use_qpoint_data_file,iret)==0 ) then
                 use_qpoint_data_file = iret
              endif
          endif

          if(qvec_from_file==OFF)then
             if( f_selectBlock( tag_k_vector) == 0) then
                if( f_getRealValue( tag_kx, dret, '') == 0) kvec_dir(1) = dret
                if( f_getRealValue( tag_ky, dret, '') == 0) kvec_dir(2) = dret
                if( f_getRealValue( tag_kz, dret, '') == 0) kvec_dir(3) = dret
                iret = f_selectParentBlock()
                dnorm = sqrt(sum(kvec_dir(1:3)**2))
                if(abs(dnorm) < 1.d-10) then
                   kvec_dir(1:2)=0.d0
                   kvec_dir(3)=1.d0
                   dnorm = 1.d0
                end if
                kvec_dir(1:3) = kvec_dir(1:3)/dnorm
                iret = f_selectParentBlock()
             else
                kvec_dir(3) = 1.d0  ! default value
             end if
          endif

          ! === KT_add === 13.1R
          if ( sw_phonon_with_epsilon == ON ) then
             call m_Raman_initialize( kvec_dir )
             dielectric = dielectric_read_from_file
          endif
          ! ================ 13.1R

          if( phonon_method /= PHONON_GAMMA) then
             if( f_selectBlock( tag_lattice) == 0) then
                if( f_getIntValue( tag_l1, iret) == 0) nlat(1) = iret
                if( f_getIntValue( tag_l2, iret) == 0) nlat(2) = iret
                if( f_getIntValue( tag_l3, iret) == 0) nlat(3) = iret
                iret = f_selectParentBlock()
             else
                nlat(1:3) = 2
             end if
          else
             nlat(1:3) = 0
             if( f_selectBlock( tag_lattice) == 0) then
                if( f_getIntValue( tag_l1, iret) == 0) nlat(1) = iret
                if( f_getIntValue( tag_l2, iret) == 0) nlat(2) = iret
                if( f_getIntValue( tag_l3, iret) == 0) nlat(3) = iret
                iret = f_selectParentBlock()
             endif
          end if

          if( phonon_method == PHONON_BAND) then
             call rd_spcpnt_and_symln()
          else if( phonon_method == PHONON_DOS) then
             phonon_dos_mesh(1:3) = 4 ! default values
             if( f_selectBlock( tag_dos) == 0) then
                if( f_selectBlock( tag_mesh) == 0) then
                   if( f_getIntValue( tag_nx, iret) == 0) phonon_dos_mesh(1) = iret
                   if( f_getIntValue( tag_ny, iret) == 0) phonon_dos_mesh(2) = iret
                   if( f_getIntValue( tag_nz, iret) == 0) phonon_dos_mesh(3) = iret
                   iret = f_selectParentBlock()
                end if
                if( f_getRealValue( tag_deltae, dret, '') == 0) DeltaE_dos = dret
                iret = f_selectParentBlock()
             end if
          end if
          if(f_getIntValue(tag_sw_read_forces_pre,iret)==0) sw_read_forces_pre = iret

          iret = f_selectParentBlock()
       end if

       iret = f_selectParentBlock()
    end if

    if(sw_lattice_dielectric_tensor == ON) then
       phonon_method = PHONON_GAMMA
       sw_lo_to_splitting = OFF
    end if

    if(ipriphonon >= 1) then
       write(nfout,*) '!** sw_phonon = ',sw_phonon
       if(sw_phonon == ON) then
          write(nfout,*) '!** phonon_method = ',phonon_method
          write(nfout,*) '!** sw_calc_force = ',sw_calc_force
          write(nfout,*) '!** sw_calc_force_all = ',sw_calc_force_all
          write(nfout,*) '!** istart_phonon = ',istart_phonon
          write(nfout,*) '!** iend_phonon   = ',iend_phonon
          write(nfout,*) '!** sw_read_forces_pre   = ',sw_read_forces_pre
          write(nfout,*) '!** displacement (Phonon) = ',u
          write(nfout,*) '!** sw_vibrational_modes = ',sw_vibrational_modes
          write(nfout,*) '!** sw_polynomial_fit = ',sw_polynomial_fit
          write(nfout,*) '!** norder = ',norder
          write(nfout,*) '!** nlat = ',nlat
          write(nfout,*) '!** sw_correct_force_constants = ',sw_correct_force_constants
          write(nfout,*) '!** point_group = ',pg_gamma
          write(nfout,*) '!** sw_lo_to_splitting = ',sw_lo_to_splitting
          write(nfout,'(" !**           [",3(1x,f10.5),1x,"]")') dielectric(1:3,1)
          write(nfout,'(" !** epsilon = [",3(1x,f10.5),1x,"]")') dielectric(1:3,2)
          write(nfout,'(" !**           [",3(1x,f10.5),1x,"]")') dielectric(1:3,3)
          write(nfout,'(" !** k_vector =",3(1x,f10.5))') kvec_dir(1:3)
          write(nfout,*) '!** sw_lattice_dielectric_tensor = ',sw_lattice_dielectric_tensor
          write(nfout,*) '!** sw_dielectric_function = ',sw_dielectric_function
          write(nfout,*) '!** energy_range = ',min_energy,max_energy,num_division
          write(nfout,*) '!** damping_factor = ',damping_factor
          write(nfout,*) '!** sw_internal_strain_piezoelectric_tensor = ',sw_int_strain_piezo_tensor
          if( phonon_method == PHONON_DOS) then
             write(nfout,*) '!** phonon_dos_mesh = ',phonon_dos_mesh(1:3)
          end if
          if( phonon_method == PHONON_BAND) then
             write(nfout,*) '!** use_qpoint_data_file = ', use_qpoint_data_file
             if ( sw_phonon_band_unfolding == ON ) then
                write(nfout,*) '!** sw_phonon_band_unfolding = ', &
                     &          sw_phonon_band_unfolding
                write(nfout,*) "!** tolerance_phonon_Gvec_matching is ", &
                     &              tolerance_ph_Gvec_matching
                write(nfout,*) '!** nmax_G_phonon(1) = ', nmax_G_phonon(1)
                write(nfout,*) '!** nmax_G_phonon(2) = ', nmax_G_phonon(2)
                write(nfout,*) '!** nmax_G_phonon(3) = ', nmax_G_phonon(3)
             endif
          endif

          write(nfout,*) '!** sw_phonon_oneshot = ', sw_phonon_oneshot
          if(sw_phonon_oneshot == ON) then
             write(nfout,*) '!** num_phonon_calc_mode = ',num_phonon_calc_mode
          end if
          ! ====== KT_add ==== 13.1R
          if ( sw_phonon_with_epsilon == ON ) call m_Raman_print_param
          ! ================== 13.1R
       end if
    end if

  contains
    subroutine rd_spcpnt_and_symln()
      integer :: num

      call read_spcpnt(.true.,nsp,num)  ! paramset = .true.  -> nsp
      nsp = num
      allocate(spcpnt(nsp,3))
      allocate(pg_spcpnt(nsp))
      call read_spcpnt(.false.,nsp,num) ! paramset = .false. -> spcpnt,pg_spcpnt

      call read_symln(.true.,nsl,num)  ! paramset = .true.  -> nsl
      nsl = num
      allocate(symln(nsl,3))
      allocate(ndiv(nsl))
      allocate(pg_symln(nsl))
      call read_symln(.false.,nsl,num) ! paramset = .false. -> symln,ndiv,pg_symln
    end subroutine rd_spcpnt_and_symln

    subroutine read_spcpnt(paramset,num,icounted)
      logical, intent(in) :: paramset
      integer, intent(in) ::  num
      integer, intent(out) :: icounted

      integer :: i, no, iret
      integer :: f_selectFirstTableLine,f_selectNextTableLine,f_getRealValue,f_getIntValue
      real(kind=dP) :: dret
      character(len=FMAXVALLEN) :: cret

      i = 1
      if( f_selectBlock( tag_special_points) == 0) then
         do
            if(i == 1) then
               if( f_selectFirstTableLine() /= 0 ) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0 ) then
                  exit
               end if
            end if
            if(.not.paramset) then
               if(i > num) exit
               no = i
               if(f_getIntValue(tag_no, iret) == 0) then
                  no = iret
               else if(f_getIntValue(tag_id, iret) == 0) then
                  no = iret
               end if
               if(no <= num) then
                  if(f_getRealValue(tag_k1,dret,'') == 0) spcpnt(no,1) = dret
                  if(f_getRealValue(tag_k2,dret,'') == 0) spcpnt(no,2) = dret
                  if(f_getRealValue(tag_k3,dret,'') == 0) spcpnt(no,3) = dret
                  if( f_getStringValue( tag_point_group, cret, NOCONV) == 0) then
                     pg_spcpnt(no) = cret(1:3)
                  else
                     pg_spcpnt(no) = '' ! default value
                  end if
               end if
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      end if
      icounted = i - 1

    end subroutine read_spcpnt

    subroutine read_symln(paramset,num,icounted)
      logical, intent(in) :: paramset
      integer, intent(in) ::  num
      integer, intent(out) :: icounted

      integer :: i, no, iret
      integer :: f_selectFirstTableLine,f_selectNextTableLine,f_getRealValue,f_getIntValue
      character(len=FMAXVALLEN) :: cret

      i = 1
      if( f_selectBlock( tag_symmetry_lines) == 0) then
         do
            if(i == 1) then
               if( f_selectFirstTableLine() /= 0 ) then
                  exit
               end if
            else
               if( f_selectNextTableLine() /= 0 ) then
                  exit
               end if
            end if
            if(.not.paramset) then
               if(i > num) exit
               no = i
               if(f_getIntValue(tag_no, iret) == 0) then
                  no = iret
               else if(f_getIntValue(tag_id, iret) == 0) then
                  no = iret
               end if
               if(no <= num) then
                  if(f_getIntValue(tag_k1,iret) == 0) symln(no,1) = iret
                  if(f_getIntValue(tag_k2,iret) == 0) symln(no,2) = iret
                  if( f_getStringValue( tag_point_group, cret, NOCONV) == 0) then
                     pg_symln(no) = cret(1:3)
                  else
                     pg_symln(no) = '' ! default value
                  end if
                  if(f_getIntValue(tag_num_division,iret) == 0) ndiv(no) = iret
                  if(iret <= 0) ndiv(no) = 1
               end if
            end if
            i = i+1
         end do
         iret = f_selectParentBlock()
      end if
      icounted = i - 1

    end subroutine read_symln

    subroutine set_phonon_method(rstr)
      character(len=FMAXVALLEN),intent(in) :: rstr
      logical :: tf
      call strncmp2(rstr, FMAXVALLEN, tag_gamma, len(tag_gamma), tf)
      if(.not.tf) &
           & call strncmp2(rstr, FMAXVALLEN, tag_zone_center, len(tag_zone_center), tf)
      if(tf) then
         phonon_method = PHONON_GAMMA
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_band, len(tag_band), tf)
      if(.not.tf) &
           & call strncmp2(rstr, FMAXVALLEN, tag_dispersion, len(tag_dispersion), tf)
      if(tf) then
         phonon_method = PHONON_BAND
         goto 1001
      end if
      call strncmp2(rstr, FMAXVALLEN, tag_dos, len(tag_dos), tf)
      if(tf) then
         phonon_method = PHONON_DOS
         goto 1001
      end if
1001  continue
    end subroutine set_phonon_method

  end subroutine m_Phonon_rd_param

  subroutine m_Phonon_alloc_qvec()
    integer :: isl,ishf
    real(kind=DP),allocatable,dimension(:,:,:) :: vkxyz_t
    real(kind=DP),allocatable,dimension(:) :: qwgt_t
    integer :: qnv3
    if(phonon_method == PHONON_GAMMA) then
       nqvec = 1
    else if(phonon_method == PHONON_BAND) then
       if(qvec_from_file==OFF)then
          nqvec = 0
          do isl=1,nsl
             ishf = 0
             if(isl/=nsl.and.symln(isl,2)==symln(isl+1,1)) ishf = 1
             nqvec = nqvec + ndiv(isl) + 1 - ishf
          end do
       else
          allocate(vkxyz_t(1,1,1))
          allocate(qwgt_t(1))
          qnv3=0

          if ( use_qpoint_data_file == OFF ) then
             call m_Files_open_kpoint_files(FILE,nbztyp_spg)
             call readk0(.true.,nfout,ipriphonon,qnv3,rltv,nfkpoint,nfmatbp &
                  & ,nqvec,vkxyz_t,qwgt_t)
          else
             call m_Files_open_qpoint_files(FILE,nbztyp_spg)
             call readk0(.true.,nfout,ipriphonon,qnv3,rltv,nfqpoint,nfmatbp &
                  & ,nqvec,vkxyz_t,qwgt_t)
          endif
          write(nfout,'(a,i5)') &
               & ' !** number of q-points for the phonon-band calculation : ',nqvec

          deallocate(vkxyz_t)
          deallocate(qwgt_t)
       endif
    else if(phonon_method == PHONON_DOS) then
       call make_qpoint_set(0,nqvec,.true.) ! ipri=0, paramset = .true.
       allocate(wght(nqvec))
    end if

    allocate(qvec(nqvec,3))
    allocate(qvin(nqvec,3))
    allocate(qdir(nqvec,3))
    allocate(point_group(nqvec)); point_group = ''
    allocate(fksym(48,nqvec)); fksym = .false.

    if ( sw_phonon_band_unfolding == ON ) allocate(qvin_refcell(nqvec,3))

  end subroutine m_Phonon_alloc_qvec

  subroutine m_Phonon_set_qvec(nfout)
    integer, intent(in) :: nfout

    call set_qvec(nfout)
    call set_qvec_point_group(nfout)
    call wd_qvec(nfout)
  contains
    subroutine set_qvec(nfout)
      integer,intent(in) :: nfout

      integer :: k1,k2
      integer :: i,isl,ishf,iqvec
      real(kind=DP) :: q1(3),q2(3),dq(3)
      real(kind=DP), allocatable, dimension(:,:,:) :: qvec_t, qvec_t_refcell
      real(kind=DP), allocatable, dimension(:) :: qw_t

      if(phonon_method == PHONON_GAMMA) then
         nqvec = 1
         qvec(nqvec,1:3) = 0.d0
         qvin(nqvec,1:3) = 0.d0
         qdir(nqvec,1:3) = kvec_dir(1:3)
         point_group(nqvec) = pg_gamma
      else if(phonon_method == PHONON_BAND) then
         if(qvec_from_file==OFF)then
            nqvec = 0
            do isl=1,nsl
               k1 = symln(isl,1)
               k2 = symln(isl,2)
               q1(1:3) = rltv(1:3,1)*spcpnt(k1,1) &
                    & + rltv(1:3,2)*spcpnt(k1,2) &
                    & + rltv(1:3,3)*spcpnt(k1,3)
               q2(1:3) = rltv(1:3,1)*spcpnt(k2,1) &
                    & + rltv(1:3,2)*spcpnt(k2,2) &
                    & + rltv(1:3,3)*spcpnt(k2,3)
               dq(1:3) = (q2(1:3)-q1(1:3))/ndiv(isl)
               ishf = 0
               if(isl/=nsl.and.symln(isl,2)==symln(isl+1,1)) ishf = 1
               iqvec = nqvec + 1
               do i=0,ndiv(isl)-ishf
                  nqvec = nqvec + 1
                  qvec(nqvec,1:3) = i*dq(1:3)+q1(1:3)
                  qdir(nqvec,1:3) = dq(1:3)
                  point_group(nqvec) = pg_symln(isl)
               end do
               point_group(iqvec) = pg_spcpnt(k1)
               if(ishf==0) point_group(nqvec) = pg_spcpnt(k2)
            end do
         else
            allocate(qvec_t(nqvec,3,CRDTYP));qvec_t=0.d0
            allocate(qw_t(nqvec));qw_t=1.d0

            if ( sw_phonon_band_unfolding == ON ) then
               allocate(qvec_t_refcell(nqvec,3,CRDTYP));qvec_t_refcell=0.d0
            endif

            if ( sw_phonon_band_unfolding == OFF ) then
               if ( use_qpoint_data_file == OFF ) then
                  call readk0(.false.,nfout,ipriphonon,nqvec,rltv,nfkpoint,nfmatbp &
                       & , nqvec, qvec_t, qw_t)
               else
                  call readk0(.false.,nfout,ipriphonon,nqvec,rltv,nfqpoint,nfmatbp &
                       & , nqvec, qvec_t, qw_t)
               endif
            else
               if ( use_qpoint_data_file == OFF ) then
                  call readk0( .false., nfout, ipriphonon, nqvec, rltv_refcell, &
                       &       nfkpoint, nfmatbp, nqvec, qvec_t_refcell, qw_t )
!                  call readk0_for_band_unfolding( .false., nfout, ipriphonon, &
!                       &       nqvec, rltv_refcell, altv, &
!                       &       nfkpoint, nfmatbp, nqvec, qvec_t, qw_t )
               else
                  call readk0( .false., nfout, ipriphonon, nqvec, rltv_refcell, &
                       &       nfqpoint, nfmatbp, nqvec, qvec_t_refcell, qw_t )
!                  call readk0_for_band_unfolding( .false., nfout, ipriphonon, &
!                       &       nqvec, rltv_refcell, altv, &
!                       &       nfqpoint, nfmatbp, nqvec, qvec_t, qw_t )
               endif
            endif

            qdir=0.d0
            if ( sw_phonon_band_unfolding == OFF ) then
               do i=1,nqvec
                  qvec(i,1:3) = qvec_t(i,1:3,CARTS)
               enddo
            else
               do i=1,nqvec
                  qvec(i,1:3) = qvec_t_refcell(i,1:3,CARTS)
               enddo
            endif
            do i=1,nqvec-1
               qdir(i,1) = qvec(i+1,1)-qvec(i,1)
               qdir(i,2) = qvec(i+1,2)-qvec(i,2)
               qdir(i,3) = qvec(i+1,3)-qvec(i,3)
            enddo
            if(nqvec>1) qdir(nqvec,:) = qdir(nqvec-1,:)
         endif
         do iqvec=1,nqvec
            do i=1,3
               qvin(iqvec,i) = sum(altv(1:3,i)*qvec(iqvec,1:3))/PAI2
            end do
         end do
         if ( sw_phonon_band_unfolding == ON ) then
            do iqvec=1,nqvec
               do i=1,3
                  qvin_refcell(iqvec,i) &
                       &  = sum(altv_refcell(1:3,i)*qvec(iqvec,1:3))/PAI2
               end do
            end do
         endif

      else if(phonon_method == PHONON_DOS) then
         call make_qpoint_set(1,nqvec,.false.) ! ipri=1, paramset = .false.
         do iqvec=1,nqvec
            do i=1,3
               qvin(iqvec,i) = sum(altv(1:3,i)*qvec(iqvec,1:3))/PAI2
            end do
            qdir(iqvec,1:3) = rltv(1:3,1)+rltv(1:3,2)+rltv(1:3,3)
         end do
      end if
    end subroutine set_qvec

    subroutine set_qvec_point_group(nfout)
      use m_Const_Parameters,  only: AUTOMATIC
      use m_Crystal_Structure, only: nopr,op,ig01,il,symmetry_method
      use m_CS_SpaceGroup,     only: lattice_system, &
           & get_kpoint_symmetry,get_point_group_name
      integer,intent(in) :: nfout

      integer :: i,j
      integer :: nopro,ig01o(48)
      integer :: nopro2,ig01o2(48)
      real(kind=DP) :: kvec(3)
      character(len=3) :: pg_name
      character(len=5) :: pg_name_i
      character(len=9) :: system

      if(symmetry_method==AUTOMATIC) then
         system = lattice_system
      else
         if(il>0) then
            system = 'cubic'
         else
            system = 'hexagonal'
         end if
      end if

      do i=1,nqvec
         if(point_group(i) == '') then
            kvec = qvec(i,1:3)
            call get_kpoint_symmetry(nopr,ig01,op,kvec,nopro,ig01o,nopro2,ig01o2)
            do j=1,nopro2
               fksym(ig01o2(j),i) = .true.
            end do
            call get_point_group_name(nfout,system,nopro,ig01o,pg_name,pg_name_i)
            point_group(i) = pg_name
         else
            fksym(1:nopr,i) = .true.
         end if
      end do
    end subroutine set_qvec_point_group

    subroutine wd_qvec(nfout)
      integer,intent(in) :: nfout

      integer :: i
      if(ipriphonon > 0) then
         write(nfout,'("=== q-vectors ===")')
         write(nfout,'("number of q-vectors = ",i5)') nqvec
         write(nfout,'(4x,"i",5x,"qvec(3)",26x,"qvin(3)",22x,"point_group")')
         do i=1,nqvec
            write(nfout,'(i5,3(1x,f10.5),3(1x,f10.5),1x,a3)') i,qvec(i,1:3),qvin(i,1:3),point_group(i)
         end do
      end if
    end subroutine wd_qvec

  end subroutine m_Phonon_set_qvec

  subroutine make_qpoint_set(ipri,nqvec,paramset)
    use m_Crystal_Structure, only: il,inv,ngen,igen,jgen,a,b,c,ca,cb,cc
    integer, intent(in) :: ipri
    integer, intent(inout) :: nqvec
    logical, intent(in) :: paramset

    integer :: nx1,ny1,nz1
    integer :: nx,ny,nz
    integer :: nxx,nyy,nzz
    integer :: nd, lmnp0, lmnp1, lmnp2
    real(kind=DP), pointer, dimension(:,:) :: pa0,pb0,pb
    integer,       pointer, dimension(:,:) :: ka0,ka2
    integer,       pointer, dimension(:,:) :: ip2cub_wk
    integer,       pointer, dimension(:)   :: nstar2
    integer                                :: i
    real(kind=DP) :: nv(3),trmat(3,3)

    nx = phonon_dos_mesh(1)
    ny = phonon_dos_mesh(2)
    nz = phonon_dos_mesh(3)
    call nskma0(il,nx,ny,nz,nxx,nyy,nzz,nx1,ny1,nz1,nd)

    lmnp0=(nxx+1)*(nyy+1)*(nzz+1)
    lmnp1=lmnp0
    lmnp2=lmnp0

    nxyz_tetra = 0
    allocate(ip10(lmnp0))  ; ip10 = 0
    if(paramset) allocate(ip20(lmnp0))
    ip20 = 0
    allocate(ip01(lmnp1))    ; ip01 = 0
    allocate(ip02(lmnp2))    ; ip02 = 0
    allocate(ip21(lmnp1))    ; ip21 = 0
    allocate(ip12(lmnp2))    ; ip12 = 0
    allocate(iu21(lmnp1))    ; iu21 = 0
    allocate(iv21(lmnp1))    ; iv21 = 0
    allocate(nstar2(lmnp2))  ; nstar2 = 0
    allocate(pa0(3,lmnp0))  ; pa0 = 0
    allocate(pb0(3,lmnp0))  ; pb0 = 0
    allocate(pb(3,lmnp2))  ; pb = 0
    allocate(ka0(4,lmnp0))  ; ka0 = 0
!!$  allocate(ka2(4,lmnp2))  ; ka2 = 0

    call setkp0_default_n(il,ngen,inv,igen,jgen,a,b,c,ca,cb,cc &
         &               ,nx,ny,nz &
         &               ,np2,np1,np0,lmnp0,lmnp1,lmnp2 &
         &               ,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3)&
         &               ,ip10,ip20,ip01,ip02,ip21,ip12,iu21,iv21 &
         &               ,nstar2,pa0,pb0,pb,ka0 &
         &               ,ipri,itrs)

    if(.not.paramset) then
       call get_trmat1  !-(contained here) ->(trmat)
       do i=1,nqvec
          nv(1:3) = dble(ka0(1:3,ip02(i)))/dble(ka0(4,ip02(i)))
          qvec(i,1:3) = matmul(trmat,nv)
          wght(i)     = dble(nstar2(i))/dble(np1)
       end do
    else
       nqvec = np2
    end if

    deallocate(ip10)
    deallocate(ip01)
    deallocate(ip02)
    deallocate(ip21)
    deallocate(ip12)
    deallocate(iu21)
    deallocate(iv21)
    deallocate(nstar2)
    deallocate(pa0)
    deallocate(pb0)
    deallocate(pb)
    deallocate(ka0)
!!$  deallocate(ka2)

    if(paramset) return

    allocate(iwt(np2)) ; iwt = 0
    allocate(ip2cub(np1)) ; ip2cub = 0

    allocate(ip2cub_wk(9,nxyz_tetra(1)*nxyz_tetra(2)*nxyz_tetra(3)))
    ip2cub_wk = 0.d0
    call wtetra &
         &  (nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3),np0,np2,ip20 &
         &  ,iwt,ip2cub &
         &  ,ip2cub_wk)
    deallocate(ip2cub_wk)
  contains
    subroutine get_trmat1
      real(kind=DP), dimension(3,3) :: trbp,trpb,mat1,mat2
      !    make translation matrix  trpb (P -> B)

      call getspgtab(trbp)  ! spg+tetra

      call inver3n(3,trbp,trpb)

      mat1 = transpose(trpb)
      call inver3n(3,mat1,mat2)
      call matpr3(rltv,mat2,trmat)

    end subroutine get_trmat1
  end subroutine make_qpoint_set

  subroutine m_Phonon_write_forces()
    use m_IterationNumbers, only : iteration_ionic, iteration_electronic
    use m_Files, only : nfforce, m_Files_open_nfforce
    use m_Ionic_system, only : natm,num_force_calc,iconf
    use m_Control_Parameters, only : num_phonon_calc_mode
    use m_Force, only : forc_l
#ifdef WINDOWS
    use ifport
#endif

    ! local variables
    integer :: i,ic
    logical :: ini

    if(.not.allocated(force_data)) allocate(force_data(natm_super,3,num_force_data))
    ic = iconf(iteration_ionic)
    force_data(1:natm_super,1:3,ic) = forc_l(1:natm_super,1:3)

    if(mype /= 0) return

    ini = .false.
    do i=1,3*natm_prim*norder*2
       if(force_was_read(i)) then
         ini = .true.
         exit
       endif
    enddo
!    if(iteration_ionic == 1) then
    if(phonon_iteration == 1 .and. .not. ini) then
       call m_Files_open_nfforce(.false.)  ! new file (asms)
       if(istart_phonon == 1) then
          if(sw_phonon_oneshot == ON) then
             write(nfforce,'(4(1x,i7),"  : num_force_calc, norder, sw_polynomial_fit, num_phonon_calc_mode")') &
                  & num_force_calc,  norder, sw_polynomial_fit, num_phonon_calc_mode
          else
             write(nfforce,'(3(1x,i7),"  : num_force_calc, norder, sw_polynomial_fit")') &
                  & num_force_calc, norder, sw_polynomial_fit
          end if
       end if
    else
       call m_Files_open_nfforce(.true.)  ! append (asms)
    end if
    write(nfforce,'(i4,3(1x,e25.12),2i6,"  : displaced_atom, displacement(1:3), iteration_ionic,iteration_electronic")') &
         & displaced_atom,displacement(1:3),iteration_ionic,iteration_electronic
    do i=1,natm_super
       write(nfforce,'(i4,3(1x,e25.12))') i,forc_l(i,1:3)
    end do
    call flush(nfforce)

  end subroutine m_Phonon_write_forces

  ! ==== KT_add == 13.1R
  subroutine m_Phonon_print_Raman_susc
    call m_Raman_calc_susceptibility( nmodes, omega, modes, &
         &                            mode_active, mode_irr_rep )
  end subroutine m_Phonon_print_Raman_susc
  ! =============== 13.1R

  function m_Phonon_Check_iteration(iteration_ionic)
    use m_Control_Parameters, only : skip_alloc_phonon
    use m_Crystal_Structure, only : m_CS_phonon_symmetry
    use m_Ionic_System, only : m_IS_phonon_equilibrium, m_IS_cps_to_pos &
         & , istart_phonon, iend_phonon, num_force_calc &
         & , iconf
    logical :: m_Phonon_Check_iteration
    integer, intent(in) :: iteration_ionic
    integer :: id

    if(ipriphonon>=1) then
       id = iteration_ionic+istart_phonon-1
       write(nfout,*) 'PHONON: num_force_calc  = ',num_force_calc
       write(nfout,*) 'PHONON: istart, iend = ',istart_phonon, iend_phonon
       write(nfout,*) 'PHONON: id = ',id
       write(nfout,*) 'PHONON: iconf = ',iconf(id)
       write(nfout,*) 'PHONON: iteration_ionic = ',iteration_ionic
    end if
    if(iteration_ionic >= iend_phonon-istart_phonon+1 ) then
       m_Phonon_Check_iteration = .true.
       call m_IS_phonon_equilibrium()
       call m_IS_cps_to_pos()
       call m_CS_phonon_symmetry(ON)
    else
       m_Phonon_Check_iteration = .false.
    end if
    skip_alloc_phonon = .true.

  end function m_Phonon_Check_iteration

  subroutine m_Phonon_read_forces()
    use m_Files,        only : m_Files_open_nfforce, nfforce
    use m_Ionic_System, only : num_force_calc, napt_phonon, iequconf &
         &                    , iopr_equconf, iconf, phonon_atom,imdtyp
    use m_Crystal_Structure, only : op
#ifdef WINDOWS
    use ifport
#endif
    implicit none

    ! local variables
    integer :: n,i,ia,ja,ind,j,ic,iopr
    integer :: istart,idummy,iter
    character(len=5) :: num
    integer :: icount,jcount
    logical :: exi

    if(ipriphonon>=1) write(nfout,*) '<< m_Phonon_read_forces >>: Reading forces'
    call flush(nfout)
!    nmodes = natm_prim*3
    nmodes = natm_mobile*3

    !if(sw_calc_force == OFF) then
    !   if(ipriphonon>=2) write(nfout,'(" sw_calc_force = OFF")')

    !num_force_data = nmodes*norder*2
    num_force_data = natm_prim*3*norder*2
    if(sw_polynomial_fit == ON) num_force_data = num_force_data + 1

    if(sw_polynomial_fit == ON) then
      if(.not.allocated(forces)) allocate(forces(3*natm_mobile*nlpnt,nmodes,norder*2+1))
    else
      if(.not.allocated(forces)) allocate(forces(3*natm_mobile*nlpnt,nmodes,norder*2))
    endif
    if(.not.allocated(force0)) allocate(force0(3*natm_mobile*nlpnt))
    !if(.not.allocated(force_data)) allocate(force_data(natm_super,3,num_force_data))
    if(.not.allocated(force_data)) allocate(force_data(natm_super,3,num_force_data))
    if(.not.allocated(force_was_read)) then
      allocate(force_was_read(num_force_data));force_was_read = .false.
    endif

    if(mype == 0) then
       call m_Files_open_nfforce(.false.)
       inquire(unit=nfforce,exist=exi)
       if(exi)then
         rewind nfforce
         read(nfforce,*,end=100,err=100) num_force_calc, norder, sw_polynomial_fit
         do while(.true.)
            read(nfforce,*,end=100,err=100) displaced_atom,displacement(1:3),iter,idummy
            ic = iconf(iter)
            do ia=1,natm_super
               read(nfforce,*,end=100,err=100) idummy,force_data(ia,1:3,ic)
            end do
            force_was_read(ic) = .true.
         end do
100      continue
       endif
    end if
    !else
    !   if(ipriphonon>=2) write(nfout,'(" sw_calc_force = ON")')
    !end if

    if(ipriphonon>=2) write(nfout,'(" num_force_data = ",i8)') num_force_data
    do i=1,num_force_data
       ic = iequconf(i)
       if(ipriphonon>=2) write(nfout,'(" i, ic = ",2i8)') i, ic
       if(imdtyp(phonon_atom(i))==OFF) cycle
       if(ic >= 0) cycle
       ic = abs(ic)
       iopr = iopr_equconf(i)
       if(ipriphonon>=2) write(nfout,'(" i, ic, iopr = ",3i8)') i, ic,iopr
       do ia=1,natm_super
          ja = napt_phonon(ia,i)
          do j=1,3
             force_data(ia,j,ic) = (force_data(ia,j,ic) + dot_product(op(1:3,j,iopr),force_data(ja,1:3,i)))*0.5d0
          end do
       end do
    end do

    if(ipriphonon>=2) write(nfout,'(" num_force_data = ",i8)') num_force_data
    do i=1,num_force_data
       ic = iequconf(i)
       if(ic == 0) cycle
       ic = abs(ic)
       iopr = iopr_equconf(i)
       if(ipriphonon>=2) write(nfout,'(" i, ic, iopr = ",3i8)') i, ic,iopr
       do ia=1,natm_super
          ja = napt_phonon(ia,i)
          do j=1,3
             force_data(ja,j,i) =  dot_product(op(j,1:3,iopr),force_data(ia,1:3,ic))
          end do
       end do
    end do

    if(sw_polynomial_fit == ON) then
       ind=1
       icount = 0
       do ia=1,natm_super
          if(imdtyp(ia)==OFF) cycle
          icount = icount+1
          !istart = 3*(ia-1)+1
          istart = 3*(icount-1)+1
          force0(istart:istart+2) = force_data(ia,1:3,ind)
       end do
    else
       ind=0
    end if
!    do i=1,nmodes
    jcount = 0
    do i=1,3*natm_prim
       if (imdtyp((i-1)/3+1) /= OFF) jcount = jcount+1
       do n=1,norder*2
          ind = ind+1
          if (imdtyp((i-1)/3+1) == OFF) cycle
          icount = 0
          do ia=1,natm_super
             if(imdtyp(ia)==OFF) cycle
             icount = icount+1
             !istart = 3*(ia-1)+1
             istart = 3*(icount-1)+1
             !forces(istart:istart+2,i,n) = force_data(ia,1:3,ind)
             forces(istart:istart+2,jcount,n) = force_data(ia,1:3,ind)
          end do
       end do
    end do

    if(npes>1) then
       call mpi_bcast(forces,3*natm_mobile*nlpnt*nmodes*norder*2,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(force0,3*natm_mobile*nlpnt,mpi_double_precision,0,MPI_CommGroup,ierr) ! MPI
       call mpi_bcast(force_was_read,num_force_data,mpi_logical,0,MPI_CommGroup,ierr)
    end if
    if(ipriphonon>1)then
       do ia=1,3*natm_prim*norder*2
          write(nfout,'(a,i8,l2)') ' !** force was read ',ia,force_was_read(ia)
       enddo
    endif

    if(ipriphonon>=1) write(nfout,*) '<< m_Phonon_read_forces >>: Forces are read.'
    call flush(nfout)

    ! check sum of force
    !!if(ipriphonon>=0) then
    if(ipriphonon>=2) then
       write(nfout,*) "check sum of force"
       do i=1,num_force_data
          ic = iconf(i)
          if(ic>num_force_data) cycle
          if(imdtyp(phonon_atom(ic))==OFF) cycle
          do j=1,3
             write(nfout,*) i,j,sum(force_data(1:natm_super,j,i))
          end do
          !do ia=1,natm_super
          !   write(nfout,'(a,2i5,3f20.10)') 'ia, i, f ',ia,i,force_data(ia,1:3,i)
          !enddo
       end do
    end if
    return
  end subroutine m_Phonon_read_forces


  subroutine m_Phonon_calc_dynamical_matrix()
    use m_Crystal_Structure, only : univol, altv
    use m_Ionic_System, only : imdtyp
    implicit none

    ! local variables
    integer :: i,j,k
    integer :: istart
    real(kind=DP) :: mi,ui

    real(kind=DP) :: coefna
    complex(kind=DP), allocatable :: nadynmat(:,:,:) ! dim(nmodes,nmodes,nqvec)
    complex(kind=DP), allocatable :: nafd(:,:,:) ! dim(3,3,natm)
    real(kind=DP) :: sum_rule(nmodes,3)
    integer :: na,ia,ja,ir,jr,ira,jra,n
    real(kind=DP) :: ae,ai,amelem
    integer :: icell,nre,jas
    real(kind=DP) :: r(3),re(3,8)
    real(kind=DP) :: ph,w
    complex(kind=DP) :: expqr
    integer :: ii,jj,iu,ju
    integer :: iq
    real(kind=DP) :: s(3)
    integer :: ng,nr
    real(kind=DP), allocatable :: gv(:,:),rv(:,:),cpstmp(:,:)
    real(kind=DP) :: g2bz

    real(kind=DP), parameter :: eps = 1.d-6
    integer :: icount,jcount

    allocate(mode2atom(nmodes))
    icount=0
    do i=1,natm_prim
       if(imdtyp(i) /= OFF)then
         icount = icount+1
         mode2atom(3*(icount-1)+1) = i
         mode2atom(3*(icount-1)+2) = i
         mode2atom(3*(icount-1)+3) = i
       endif
    enddo
    allocate(dynmat(nmodes,nmodes,nlpnt))
    allocate(dynmatq(nmodes,nmodes,nqvec))
    allocate(msqinv(nmodes))

    call m_Phonon_force_constant_matrix(dynmat)

    deallocate(forces)
    deallocate(force0)
    deallocate(force_was_read)

    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'calculated force constants:'
       call write_dynmat(dynmat)
    end if

    allocate(cpstmp(natm_mobile,3))
    icount=0
    do ia=1,natm
      if(imdtyp(ia) ==  OFF) cycle
      icount = icount+1
      cpstmp(icount,1:3) = cps(ia,1:3)
    enddo
    ! symmetrized force constants in Real space
    call symmetrize_force_const_Real(dynmat,natm_mobile,cpstmp)
    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'symmetrized force constants:'
       call write_dynmat(dynmat)
    end if

    ! check sum rule of force constants
    if(ipriphonon >=2) &
         & write(nfout,*) 'j ir sum_of_force_constants'

    icount = 0
    do ia=1,natm_prim
       if(imdtyp(ia)==OFF) cycle
       icount = icount+1
       do ir=1,3
          !i = (ia-1)*3+ir
          i = (icount-1)*3+ir
          do jr=1,3
             sum_rule(i,jr)=0
             do n=1,nlpnt
                jcount = 0
                do ja=1,natm_prim
                   if(imdtyp(ja)==OFF) cycle
                   jcount = jcount+1
                   !j=(ja-1)*3+jr
                   j=(jcount-1)*3+jr
                   !if(n==1.and.ia==ja) cycle
                   if(n==1.and.icount==jcount) cycle
                   sum_rule(i,jr) = sum_rule(i,jr) + dynmat(i,j,n)
                end do
             end do
             if(ipriphonon >=2) then
                !j=(ia-1)*3+jr
                j=(icount-1)*3+jr
                write(nfout,*) i,jr,sum_rule(i,jr)+dynmat(i,j,1)
             end if
          end do
       end do
    end do
    ! correct force constants
    if(sw_correct_force_constants == ON) then
       do i=1,nmodes
          ia = (i+2)/3
          do jr=1,3
             j=(ia-1)*3+jr
             dynmat(i,j,1) = - sum_rule(i,jr)
          end do
       end do
    end if

    ! output
    if(sw_correct_force_constants == ON) then
       if(printable) write(nfout,*) '*** Force constants were corrected. ***'
       if(ipriphonon >=2) then
          write(nfout,*) 'Corrected force constants:'
          call write_dynmat(dynmat)
       end if
    else
       if(printable) write(nfout,*) '*** Force constants were NOT corrected. ***'
    end if

    ! Fourier transform of force constants
    dynmatq(:,:,:) = dcmplx(0.d0,0.d0)
    do n=1,nqvec
!!$write(nfout,*) 'qvin=',qvin(n,1:3)
!!$write(nfout,*) 'qvec=',qvec(n,1:3)
       do icell=1,nlpnt
          do j=1,nmodes
             !ja = (j-1)/3 + 1
             ja = mode2atom(j)
!!$ju = mod(j-1,3)
             jas = ja + (icell-1)*natm
             do i=1,nmodes
                !ia = (i-1)/3 + 1
                ia = mode2atom(i)
!!$iu = mod(i-1,3)
!!$r(1:3) = cps(jas,1:3)-cps(ia,1:3)
                r(1:3) = cps(ia,1:3)-cps(jas,1:3)
                call get_pos_in_extcell(r,nre,re)
                expqr = dcmplx(0.d0,0.d0)
                do k=1,nre
                   ph = sum(qvec(n,1:3)*re(1:3,k))
                   expqr = expqr + dcmplx(cos(ph),sin(ph))
                end do
                expqr = expqr/dble(nre)
!!$write(nfout,*) 'ia,ja,i,j=',ia,ja,i,j
!!$write(nfout,'("r =",3f10.5)') r(1:3)
!!$do k=1,nre
!!$   write(nfout,'("k,re =",i3,3f10.5)') k,re(1:3,k)
!!$end do
!!$write(nfout,*) 'ph,expqr=',ph,expqr

                dynmatq(i,j,n) = dynmatq(i,j,n)+dynmat(i,j,icell)*expqr
!!$ii = (ja-1)*3 + iu + 1
!!$jj = (ia-1)*3 + ju + 1
!!$dynmatq(ii,jj,n) = dynmatq(ii,jj,n)+dynmat(i,j,icell)*expqr
             end do
          end do
       end do
    end do
    deallocate(dynmat)

    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'Force constant matrix in Fourier space:'
       call write_dynmatq(dynmatq)
    end if

    ! symmetrized force constants in Fourier space
    call symmetrize_force_const_Fourier(dynmatq,natm_mobile,cpstmp)

    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'Symmetrized force constants Fourier space:'
       call write_dynmatq(dynmatq)
    end if


    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'Dynamical matrix in Fourier space:'
       call write_dynmatq(dynmatq)
    end if

    ! non-analytic part
    if(sw_lo_to_splitting == ON) then
       allocate(nadynmat(nmodes,nmodes,nqvec))
       call get_g2bz(g2bz)
       do iq=1,nqvec
          s = qvin(iq,1:3)
          call mod1(s)
          if(abs(sum(s(1:3)**2)) < eps) then
             ! Gamma point
             call calc_napart_gamma(iq,nadynmat(1,1,iq))
          else
             call calc_napart_parlinski(iq,nadynmat(1,1,iq),g2bz)
          end if
       end do

       if(ipriphonon >=2) then
          ! output
          write(nfout,*) 'Non-analytic part of force constant matrix in Fourier space:'
          call write_dynmatq(nadynmat)
       end if

       do iq=1,nqvec
          do j=1,nmodes
             do i=1,nmodes
                dynmatq(i,j,iq) = dynmatq(i,j,iq) + nadynmat(i,j,iq)
             end do
          end do
       end do

       if(ipriphonon >=2) then
          ! output
          write(nfout,*) 'Total force constant matrix in Fourier space:'
          call write_dynmatq(dynmatq)
       end if

       deallocate(nadynmat)

    end if

    ! calculate 1/sqrt(mass)
    icount=0
    do i=1,natm_prim
       if(imdtyp(i) == OFF) cycle
       !istart=3*(i-1)+1
       icount = icount+1
       istart=3*(icount-1)+1
       mi = 1.d0/sqrt(ionic_mass(i))
       msqinv(istart) = mi
       msqinv(istart+1) = mi
       msqinv(istart+2) = mi
    end do
    if(ipriphonon > 1) then
       do i=1,nmodes
          write(nfout,*) 'debug: i= ',i,' msqinv = ',msqinv(i)
       end do
    end if

    ! calculate dynamical matrix in Fourier space
    do n=1,nqvec
       do j=1,nmodes
          do i=1,nmodes
             dynmatq(i,j,n) = dynmatq(i,j,n)*msqinv(i)*msqinv(j)
          end do
       end do
    end do


    do iq=1,nqvec
       !search the absolute maximum element
       amelem = abs(dynmatq(1,1,iq))
       do i=1,nmodes
          do j=i,nmodes
             ae = abs(dynmatq(i,j,iq))
             if(ae > amelem) amelem = ae
          end do
       end do
       if(ipriphonon >=2) then
          write(nfout,*) 'amelem = ',amelem
       end if
       do i=1,nmodes
          do j=1,nmodes
             ae = dble(dynmatq(i,j,iq))
             ai = dimag(dynmatq(i,j,iq))
             if(abs(ae) < amelem*1.d-14) ae = 0.d0
             if(abs(ai) < amelem*1.d-14) ai = 0.d0
             dynmatq(i,j,iq) = dcmplx(ae,ai)
          end do
       end do
    end do
    if(ipriphonon >=2) then
       ! output
       write(nfout,*) 'Dynamical matrix after removing small elements:'
       call write_dynmatq(dynmatq)
    end if

    deallocate(cpstmp)
    deallocate(mode2atom)
    return

  contains
    subroutine get_pos_in_extcell(r,nre,re)
      real(kind=DP), intent(in) :: r(3)
      integer, intent(out) :: nre
      real(kind=DP), intent(out) :: re(3,8)

      integer :: i,j,k,n(3),l
      real(kind=DP) :: s(3),p(3)
      real(kind=DP) :: eps = 1.d-6

      do i=1,3
         s(i) = sum(rltv_super(1:3,i)*r(1:3))/PAI2
      end do
      call mod1(s)
      do i=1,3
         n(i) = 1
         if(abs(s(i)-0.5d0)<eps) n(i) = 2
      end do
      do i=1,3
         if(s(i)>0.5d0-eps) s(i) = s(i) - 1.d0
      end do
      nre = 0
      do k=1,n(3)
         do j=1,n(2)
            do i=1,n(1)
               nre = nre + 1
               p(1) = s(1)+i-1
               p(2) = s(2)+j-1
               p(3) = s(3)+k-1
               do l=1,3
                  re(l,nre) = sum(altv_super(l,1:3)*p(1:3))
               end do
            end do
         end do
      end do
    end subroutine get_pos_in_extcell

    subroutine mod1(t)
      real(kind=DP), intent(inout) :: t(3)
      real(kind=DP), parameter :: eps = 1.d-6
      integer :: k
      t(1:3) = mod(t(1:3),1.d0)
      do k=1,3
         if(t(k) < -eps) t(k) = t(k) + 1.d0
         if(t(k) > 1.d0 - eps) t(k) = t(k) - 1.d0
      end do
    end subroutine mod1

    subroutine symmetrize_force_const_Real(fc,natm,cps)
      use m_Crystal_Structure, only : nopr,op,tau
      real(kind=DP), intent(inout) :: fc(nmodes,nmodes,nlpnt)
      integer, intent(in) :: natm
      real(kind=DP), intent(in) :: cps(natm,3)
      ! local variables
      integer :: i,j,icell,jcell,kcell,irc,jrc
      integer :: nr(natm,nlpnt),l(3,nlpnt),lx(3),ll(3)
      integer :: lr(natm,nlpnt),llr(nlpnt,nlpnt)
!!$  real(kind=DP) :: nonsym_fc(nmodes,nmodes,nlpnt)
      real(kind=DP),allocatable,dimension(:,:,:) :: nonsym_fc
      real(kind=DP) :: tensor(3,3)
      real(kind=DP) :: x(3),y(3)
      real(kind=DP) :: aa(3,3)
      real(kind=DP) :: eps = 1.d-6
      allocate(nonsym_fc(nmodes,nmodes,nlpnt));nonsym_fc=0.d0
      nonsym_fc(1:nmodes,1:nmodes,1:nlpnt) = fc(1:nmodes,1:nmodes,1:nlpnt)
      fc(1:nmodes,1:nmodes,1:nlpnt) = 0.d0

      do icell=1,nlpnt
         l(1:3,icell) = nint(matmul(transpose(rltv),lpnt(icell,1:3))/PAI2)
      end do

      aa = matmul(transpose(rltv_super),altv)/PAI2
      do icell=1,nlpnt
         do jcell=1,nlpnt
            ll = l(1:3,jcell)-l(1:3,icell)
            do kcell=1,nlpnt
               x = matmul(aa,ll-l(1:3,kcell))
               if(sum(abs(x-nint(x)))<eps) then
                  llr(icell,jcell) = kcell
                  exit
               end if
            end do
         end do
      end do

      do n=1,nopr
         do icell=1,nlpnt
            ATOM: do i=1,natm
               x(1:3) = cps(i,1:3) + lpnt(icell,1:3)
               y = matmul(op(1:3,1:3,n),x) + tau(1:3,n,CARTS)
               do j=1,natm
                  x = y-cps(j,1:3)
                  x = matmul(transpose(rltv),x)/PAI2
                  lx = nint(x(1:3))
                  if(sum(abs(x-lx))<eps) then
                     nr(i,icell) = j
                     do jcell=1,nlpnt
                        ll = lx-l(1:3,jcell)
                        x = matmul(aa,ll)
                        if(sum(abs(x-nint(x)))<eps) then
                           lr(i,icell) = jcell
                           exit
                        end if
                     end do
                     cycle ATOM
                  end if
               end do
               write(nfout,'("No equivalent atom for the atom ",i5," in the cell ",i5)') i,icell
               call phase_error_with_msg(nfout,'No equivalent atom',__LINE__,__FILE__)
            end do ATOM
         end do
         icell = 1
         do jcell=1,nlpnt
            do j=1,natm
               ja=3*(j-1)+1
               jr=nr(j,jcell)
               jra=3*(jr-1)+1
               jrc = lr(j,jcell)
               do i=1,natm
                  ia=3*(i-1)+1
                  ir=nr(i,icell)
                  ira=3*(ir-1)+1
                  irc = lr(i,icell)
                  kcell = llr(irc,jrc)
                  tensor(1:3,1:3)=nonsym_fc(ia:ia+2,ja:ja+2,jcell)
                  call rotate_tensor(tensor,n)
                  fc(ira:ira+2,jra:jra+2,kcell)= &
                       & fc(ira:ira+2,jra:jra+2,kcell) + tensor(1:3,1:3)
               end do
            end do
         end do
      end do
      fc(1:nmodes,1:nmodes,1:nlpnt) = fc(1:nmodes,1:nmodes,1:nlpnt)/nopr
      deallocate(nonsym_fc)
    end subroutine symmetrize_force_const_Real

    subroutine rotate_tensor(tensor,n)
      use m_Crystal_Structure, only : op
      real(kind=DP), intent(inout) :: tensor(3,3)
      integer, intent(in) :: n ! index of rotation

      ! local variables
      integer :: i,j,k,l
      real(kind=DP) :: tnew(3,3)

      tnew(1:3,1:3)=0.d0

      do l=1,3
         do i=1,3
            do k=1,3
               do j=1,3
!!$tnew(i,l)=tnew(i,l)+op(j,i,n)*tensor(j,k)*op(k,l,n)
                  tnew(i,l)=tnew(i,l)+op(i,j,n)*tensor(j,k)*op(l,k,n)
               end do
            end do
         end do
      end do

      tensor(1:3,1:3)=tnew(1:3,1:3)
    end subroutine rotate_tensor

    subroutine symmetrize_force_const_Fourier(fcq,natm,cps)
      use m_Ionic_System, only : napt
      use m_Crystal_Structure, only : nopr,op
      complex(kind=DP), intent(inout) :: fcq(nmodes,nmodes,nqvec)
      integer, intent(in) :: natm
      real(kind=DP), dimension(natm,3) :: cps
      ! local variables
      integer :: i,j,nop,iq
!!$  complex(kind=DP) :: nonsym_fcq(nmodes,nmodes,nqvec)
      complex(kind=DP),allocatable,dimension(:,:,:) :: nonsym_fcq
      complex(kind=DP) :: tensor(3,3)
      complex(kind=DP) :: expgr
      real(kind=DP) :: ph,g(3)

      allocate(nonsym_fcq(nmodes,nmodes,nqvec));nonsym_fcq=0.d0
      nonsym_fcq(1:nmodes,1:nmodes,1:nqvec) = fcq(1:nmodes,1:nmodes,1:nqvec)
      fcq(1:nmodes,1:nmodes,1:nqvec) = 0.d0

      do iq=1,nqvec
         nop = 0
         do n=1,nopr
            if(.not.fksym(n,iq)) cycle
            nop = nop + 1
            g = matmul(op(1:3,1:3,n),qvec(iq,1:3)) - qvec(iq,1:3)
            do j=1,natm
               ja=3*(j-1)+1
               jr=napt(j,n)
               jra=3*(jr-1)+1
               do i=1,natm
                  ia=3*(i-1)+1
                  ir=napt(i,n)
                  ira=3*(ir-1)+1
                  ph = dot_product(g,cps(j,1:3)-cps(i,1:3))
                  expgr = dcmplx(cos(ph),sin(ph))
                  tensor(1:3,1:3)=nonsym_fcq(ira:ira+2,jra:jra+2,iq)
                  call rotate_complex_tensor(tensor,n)
                  fcq(ia:ia+2,ja:ja+2,iq)= &
                       & fcq(ia:ia+2,ja:ja+2,iq) + tensor(1:3,1:3)*expgr
               end do
            end do
         end do
         fcq(1:nmodes,1:nmodes,iq) = fcq(1:nmodes,1:nmodes,iq)/nop
      end do
      deallocate(nonsym_fcq)
      do iq=1,nqvec
         do j=1,nmodes
            do i=1,j
               fcq(i,j,iq)=(fcq(i,j,iq)+dconjg(fcq(j,i,iq)))*0.5d0
               fcq(j,i,iq)=dconjg(fcq(i,j,iq))
            end do
         end do
      end do

    end subroutine symmetrize_force_const_Fourier


    subroutine rotate_complex_tensor(tensor,n)
      use m_Crystal_Structure, only : op
      complex(kind=DP), intent(inout) :: tensor(3,3)
      integer, intent(in) :: n ! index of rotation

      ! local variables
      integer :: i,j,k,l
      complex(kind=DP) :: tnew(3,3)

      tnew(1:3,1:3)=0.d0

      do l=1,3
         do i=1,3
            do k=1,3
               do j=1,3
                  tnew(i,l)=tnew(i,l)+op(j,i,n)*tensor(j,k)*op(k,l,n)
               end do
            end do
         end do
      end do

      tensor(1:3,1:3)=tnew(1:3,1:3)
    end subroutine rotate_complex_tensor

    subroutine write_dynmat(dynmat)
      real(kind=DP), intent(in) :: dynmat(nmodes,nmodes,nlpnt)
      integer :: icell,l1,l2,l3,n
      do icell=1,nlpnt
         l1 = nint(sum(rltv(1:3,1)*lpnt(icell,1:3))/PAI2)
         l2 = nint(sum(rltv(1:3,2)*lpnt(icell,1:3))/PAI2)
         l3 = nint(sum(rltv(1:3,3)*lpnt(icell,1:3))/PAI2)
         write(nfout,'("icell=",i4," l1=",i3," l2=",i3," l3=",i3," lp=",3f10.5)') icell,l1,l2,l3,lpnt(icell,1:3)
         do n=1,nmodes
            write(nfout,'(50(1x,e12.5))') dynmat(n,1:nmodes,icell)
         end do
      end do
    end subroutine write_dynmat

    subroutine write_dynmatq(dynmatq)
      complex(kind=DP), intent(in) :: dynmatq(nmodes,nmodes,nqvec)
      integer :: iq
      do iq=1,nqvec
         write(nfout,'(i4," q=(",f8.3,",",f8.3,",",f8.3,")")') iq,qvin(iq,1:3)
         write(nfout,'("=== Real part ===")')
         do n=1,nmodes
            write(nfout,'(50(1x,e12.5))') dble(dynmatq(n,1:nmodes,iq))
         end do
         write(nfout,'("=== Imag part ===")')
         do n=1,nmodes
            write(nfout,'(50(1x,e12.5))') dimag(dynmatq(n,1:nmodes,iq))
         end do
      end do
    end subroutine write_dynmatq

    subroutine calc_napart_gamma(iq,nafc)
      integer, intent(in) :: iq
      complex(kind=DP), intent(out) :: nafc(nmodes,nmodes)

      real(kind=DP) :: coefna,kzeff(nmodes)
      integer :: i,j,na,ia

      coefna = 0.d0
      do j=1,3
         do i=1,3
            coefna = coefna + qdir(iq,i)*dielectric(i,j)*qdir(iq,j)
         end do
      end do
      coefna = PAI4/(univol*coefna)
      i=0
      do na=1,natm
         if (imdtyp(na) == OFF) cycle
         do ia=1,3
            i=i+1
            kzeff(i) = sum(qdir(iq,1:3)*zeff(1:3,ia,na))
            !kzeff(i) = sum(qdir(iq,1:3)*zeff(1:3,mode2atom(i),na))
         end do
      end do
      do j=1,nmodes
         do i=1,nmodes
            nafc(i,j) = coefna*kzeff(i)*kzeff(j)
         end do
      end do

    end subroutine calc_napart_gamma

    subroutine calc_napart_parlinski(iq,nafc,g2bz)
      integer, intent(in) :: iq
      complex(kind=DP), intent(out) :: nafc(nmodes,nmodes)
      real(kind=DP), intent(in) :: g2bz

      real(kind=DP) :: coefna,kzeff(nmodes)
      real(kind=DP) :: expk,rho2
      integer :: i,j,na,ia

      rho2 = rho_loto**2*g2bz*0.25d0

      coefna = 0.d0
      do j=1,3
         do i=1,3
            coefna = coefna + qvec(iq,i)*dielectric(i,j)*qvec(iq,j)
         end do
      end do
      coefna = PAI4/(univol*coefna)
      i=0
      do na=1,natm
         if(imdtyp(na)==OFF)cycle
         do ia=1,3
            i=i+1
            kzeff(i) = sum(qvec(iq,1:3)*zeff(1:3,ia,na))
            !kzeff(i) = sum(qvec(iq,1:3)*zeff(1:3,mode2atom(i),na))
         end do
      end do
      expk = exp(-PAI**2*sum(qvec(iq,1:3)**2)/rho2)
      do j=1,nmodes
         do i=1,nmodes
            nafc(i,j) = coefna*kzeff(i)*kzeff(j)*expk
         end do
      end do

    end subroutine calc_napart_parlinski

    subroutine get_g2bz(g2)
      real(kind=DP), intent(out) :: g2

      integer :: i,j,k
      real(kind=DP) :: g(3),gg

      g2 = dot_product(rltv(1:3,1),rltv(1:3,1))
      do i=-2,2
         do j=-2,2
            do k=-2,2
               g(1:3) = i*rltv(1:3,1)+j*rltv(1:3,2)+k*rltv(1:3,3)
               gg = dot_product(g,g)
               if(abs(gg) > eps .and. gg < g2) then
                  g2 = gg
               end if
            end do
         end do
      end do
    end subroutine get_g2bz

  end subroutine m_Phonon_calc_dynamical_matrix

  subroutine m_Phonon_force_constant_matrix(force_const_mat)
    implicit none

    real(kind=DP), intent(out) :: force_const_mat(nmodes,nmodes,nlpnt)

    ! local variables
    integer :: is,nc,ind,i,j,jj,ic,icell
    real(kind=DP) :: ui,x
    real(kind=DP), allocatable :: cmat(:,:) !dim(nc,nc)
    real(kind=DP), allocatable :: fvec(:) !dim(nc)

    if(sw_polynomial_fit == ON) then
       nc = norder*2+1
       allocate(cmat(nc,nc)); cmat = 0.d0
       allocate(fvec(nc))
       cmat(1:nc,1) = 1.d0
       ind = 1
       do is = -norder,norder
          if(is==0) cycle
          ind=ind+1
          x = u/dble(is)
          do i=2,nc
             cmat(ind,i) = x**(i-1)
          end do
       end do
       do ic=1,nlpnt
          do j=1,nmodes
             jj = nmodes*(ic-1)+j
             do i=1,nmodes
                fvec(1) = force0(jj)
                fvec(2:nc) = forces(jj,i,1:nc-1)
                force_const_mat(i,j,ic) = -polyfit_diff(nc,cmat,fvec)
             end do
          end do
       end do
    else
       ui=dble(norder)/u
       if(printable .and. ipriphonon>=2) then
          write(nfout,*) 'u=',u
          write(nfout,*) 'ui= norder/u =',ui
       endif
       if(norder == 1) then
          ui = -0.5d0*ui
          do ic=1,nlpnt
             do j=1,nmodes
                jj = nmodes*(ic-1)+j
                do i=1,nmodes
                   force_const_mat(i,j,ic) = (forces(jj,i,2)-forces(jj,i,1))*ui
                end do
             end do
          end do
       else if(norder == 2) then
          ui = -ui/12.d0
          do ic=1,nlpnt
             do j=1,nmodes
                jj = nmodes*(ic-1)+j
                do i=1,nmodes
                   force_const_mat(i,j,ic) = (-8.d0*forces(jj,i,1)+forces(jj,i,2) &
                        &  -forces(jj,i,3)+8.d0*forces(jj,i,4))*ui
                end do
             end do
          end do
       else
          call phase_error_with_msg(nfout,'Your inputed norder is not supported for the differencial approximation.'&
                                   ,__LINE__,__FILE__)
       end if
    end if

  contains

    function polyfit_diff(nd,mat,vec)
      implicit none

      real(kind=DP) :: polyfit_diff

      integer, intent(in) :: nd
      real(kind=DP), intent(in) :: mat(nd,nd)
      real(kind=DP), intent(in) :: vec(nd)

      ! local variables
      integer :: ipiv(nd),info,i
      real(kind=DP) :: a(nd,nd),b(nd)
      real(kind=DP) :: x0,df
      real(kind=DP) :: xold
      real(kind=DP) :: eps = 1.d-10

      real(kind=DP) :: ux

      ! debug
      !  write(nfout,*) '<<polyfit_diff>>'
      !  call flush(nfout)
      ! end debug

      a = mat
      b = vec

      !debug
      !do i=1,nd
      !  write(nfout,*) 'a=',a(i,1:nd)
      !end do
      !do i=1,nd
      !  write(nfout,*) 'b=',b(i)
      !end do
      !  call flush(nfout)
      !end debug

      call dgesv(nd,1,a,nd,ipiv,b,nd,info)
      if(info .ne. 0 ) then
         if(printable) write(nfout,*) 'LAPACK routine DGETRS failure: info=',info
         call phase_error_with_msg(nfout,'LAPACK routine DGETRS failure',__LINE__,__FILE__)
      end if
      ! debug
      !  do i=1,nd
      !    print *,'i=',i,' b=',b(i)
      !  end do
      !  write(nfout,'("ux="f10.5," f=",2(1x,f15.10))') ux,vec(1),polyfunc(nd,b,0.d0)
      !  i=0
      !  do is=-norder,norder
      !    if(is==0) cycle
      !    i=i+1
      !    ux = u/dble(is)
      !    write(nfout,'("ux="f10.5," f=",2(1x,f15.10))') ux,vec(i+1),polyfunc(nd,b,ux)
      !  end do
      !  call flush(nfout)
      ! end debug

      x0 = 0.d0
      df = diff_polyfunc(nd,b,x0)
      ! debug
      !print *,'df=',df
      ! end debug
      polyfit_diff = df

    end function polyfit_diff

    function polyfunc(nd,b,x)
      implicit none

      real(kind=DP) :: polyfunc

      integer, intent(in) :: nd
      real(kind=DP), intent(in) :: b(nd),x

      integer :: i
      real(kind=DP) :: f

      f = b(1)
      do i=2,nd
         f = f + b(i)*x**(i-1)
      end do
      polyfunc = f

    end function polyfunc

    function diff_polyfunc(nd,b,x)
      implicit none

      real(kind=DP) :: diff_polyfunc

      integer, intent(in) :: nd
      real(kind=DP), intent(in) :: b(nd),x

      integer :: i
      real(kind=DP) :: df

      df = b(2)
      do i=3,nd
         df = df + b(i)*dble(i-1)*x**(i-2)
      end do
      diff_polyfunc = df

    end function diff_polyfunc

  end subroutine m_Phonon_force_constant_matrix


  subroutine m_Phonon_calc_vib_modes()
    implicit none

    ! local variables
    integer :: i,j
!!$real(DP) :: eps=1.d-8
    real(DP) :: eps=1.d-16
    real(DP) :: norm

    allocate(omega(nmodes,nqvec),modes(nmodes,nmodes,nqvec))

    write(nfout,*) 'calc_vib_modes is started.'

    do j=1,nqvec
       write(nfout,'(i4," q=(",f8.3,",",f8.3,",",f8.3,")")') j,qvin(j,1:3)
       call solve_eigenproblem(dynmatq(1,1,j) ,omega(1,j),modes(1,1,j))
       do i=1,nmodes
          if(ipriphonon>=2) &
               & write(nfout,*) 'omega**2 = ',omega(i,j)
          if(abs(omega(i,j)) < eps) omega(i,j)=0.d0
          if(omega(i,j) > 0.d0) then
             omega(i,j)=sqrt(omega(i,j))
          else
             omega(i,j)=-sqrt(abs(omega(i,j))) ! for the soft mode
          end if
          if(ipriphonon >=2) &
               & write(nfout,*) '==> omega = ',omega(i,j)
       end do
    end do

    deallocate(dynmatq)

  contains
    subroutine solve_eigenproblem(dynmat,omega,modes)
      complex(kind=DP), intent(inout) :: dynmat(nmodes,nmodes)
      real(kind=DP), intent(out)      :: omega(nmodes)
      complex(kind=DP), intent(out)   :: modes(nmodes,nmodes)

      real(kind=DP) :: abstol
      integer :: lcwork,lrwork,liwork
      complex(kind=DP), allocatable :: cwork(:) ! dim(lcwork)
      real(kind=DP), allocatable :: rwork(:) ! dim(lrwork)
      integer, allocatable :: iwork(:) ! dim(liwork)
      integer, allocatable :: ifail(:) ! dim(nmodes)
      integer :: info
      integer :: m
      integer :: il,iu
      real(DP) :: vl,vu

      real(DP), external :: dlamch

      abstol=2*dlamch('S')
      !abstol=4.45014771701440D-308
      !print *,'abstol=',abstol
!!$lwork=8*nmodes
      lcwork=4*nmodes
      lrwork=7*nmodes
      liwork=5*nmodes
!!$allocate(work(lwork))
      allocate(cwork(lcwork))
      allocate(rwork(lrwork))
      allocate(iwork(liwork),ifail(nmodes))
      lcwork=-1
      call zheevx('V','A','U',nmodes,dynmat,nmodes,vl,vu,il,iu,abstol, &
           &      m,omega,modes,nmodes,cwork,lcwork,rwork,iwork,ifail,info)
!!$call dsyevx('V','A','U',nmodes,dynmat,nmodes,vl,vu,il,iu, &
!!$ &      abstol,m,omega,modes,nmodes,work,lwork,iwork, &
!!$ &      ifail,info)
      lcwork = cwork(1)
      !print *,'lwork=',lwork
      !print *,'info=',info
      deallocate(cwork)
      allocate(cwork(lcwork))
      call zheevx('V','A','U',nmodes,dynmat,nmodes,vl,vu,il,iu,abstol, &
           &      m,omega,modes,nmodes,cwork,lcwork,rwork,iwork,ifail,info)
!!$call dsyevx('V','A','U',nmodes,dynmat,nmodes,vl,vu,il,iu, &
!!$ &      abstol,m,omega,modes,nmodes,work,lwork,iwork, &
!!$ &      ifail,info)
      !print *,'info=',info
      !print *,'ifail=',ifail(1:nmodes)
      deallocate(cwork,rwork,iwork,ifail)
    end subroutine solve_eigenproblem
  end subroutine m_Phonon_calc_vib_modes

  subroutine m_Phonon_det_irr_rep()
    use m_Representation, only : get_irr_rep
    use m_Crystal_Structure, only : nopr
    use m_Ionic_System, only : imdtyp
    implicit none

    ! local variabes
    integer :: i,j,ia,ja,ist,n
    real(DP) :: norm
    real(DP) :: eir(3),ejr(3),xir(1:nmodes)
    real(DP) :: eii(3),eji(3),xii(1:nmodes)
    integer :: nrep,nsym
    integer :: rep(48,12)
    character(len=3) :: symbol_irrep(12)
    character(len=4) :: active_rep(12)
    real(DP) :: eps = 1.d-1
    integer :: iq
    integer :: icount,jcount

    allocate(atm2mobileatm(natm_prim));atm2mobileatm=0
    icount=0
    do i=1,natm_prim
       if(imdtyp(i) /= OFF)then
         icount = icount+1
         atm2mobileatm(i) = icount
       endif
    enddo

    write(nfout,*) '** Representation **'

    allocate(mode_irr_rep(nmodes,nqvec))
    allocate(mode_active(nmodes,nqvec))
    do iq=1,nqvec
       write(nfout,'(i4," q=(",f8.3,",",f8.3,",",f8.3,")")') iq,qvin(iq,1:3)
       call get_irr_rep(point_group(iq),nsym,nrep,rep,symbol_irrep,active_rep)

       if(printable)  then
          write(nfout,'("Character table:",a3)') point_group(iq)
          do j=1,nrep
             write(nfout,'(a3,1x,48i3)') symbol_irrep(j), rep(1:nsym,j)
          end do
       end if
       ! debug
       !write(nfout,'("Product table:",a3)') point_group
       !do j=1,nrep
       !  write(nfout,'(a3,1x,48i3)') symbol_irrep(j),(sum(rep(1:nopr,j)*rep(1:nopr,i)),i=1,nrep)
       !end do
       ! end debug

       if ( printable ) then
          write(nfout,*) '** mode  repr. active **'
       endif

       do i=1,nmodes
          ! debug
          !        write(nfout,*) 'index of the mode =',i
          !        write(nfout,*) 'norm =',sum(dble(modes(1:nmodes,i,iq))**2+dimag(modes(1:nmodes,i,iq))**2)
          ! end debug
          mode_irr_rep(i,iq)='NON'
          Representation: do j=1,nrep
             xir(1:nmodes)=0.d0
             xii(1:nmodes)=0.d0
             do n=1,nopr
                icount = 0
                do ia=1,natm
                   if(imdtyp(ia)==OFF)cycle
                   !ist = 3*(ia-1)+1
                   icount = icount+1
                   ist = 3*(icount-1)+1
                   eir(1:3) = dble(modes(ist:ist+2,i,iq))
                   eii(1:3) = dimag(modes(ist:ist+2,i,iq))
                   call rotate_eigenvector(n,ia,ja,eir,ejr)
                   call rotate_eigenvector(n,ia,ja,eii,eji)
                   if ( atm2mobileatm(ja) ==0 ) cycle
                   ist = 3*(atm2mobileatm(ja)-1)+1
                   xir(ist:ist+2)=xir(ist:ist+2)+rep(n,j)*ejr(1:3)
                   xii(ist:ist+2)=xii(ist:ist+2)+rep(n,j)*eji(1:3)
                   ! debug
                   !        write(nfout,*) 'j =',j
                   !        write(nfout,*) 'n =',n
                   !        write(nfout,*) 'ia=',ia
                   !        write(nfout,*) 'eir=',eir(1:3)
                   !        write(nfout,*) 'eii=',eii(1:3)
                   !        write(nfout,*) 'ja=',ja
                   !        write(nfout,*) 'ejr=',ejr(1:3)
                   !        write(nfout,*) 'eji=',eji(1:3)
                   ! end debug
                end do
             end do
             !xi(1:nmodes) = xi(1:nmodes)/nopr
             norm = sum(xir(1:nmodes)**2+xii(1:nmodes)**2)
             ! debug
             !write(nfout,*) 'j=',j,' norm=',norm
             ! end debug
             if(norm > eps) then
                mode_irr_rep(i,iq)= symbol_irrep(j)
                mode_active(i,iq) = active_rep(j)
                ! exit Representation
             end if
          end do Representation

          if(printable) then
             write(nfout,'(1x,i5,5x,a3,4x,a4)') i,mode_irr_rep(i,iq),mode_active(i,iq)
          endif
       end do
    end do ! q loop

    deallocate(atm2mobileatm)
  contains

    subroutine rotate_eigenvector(n,ia,ja,ei,ej)
      use m_Ionic_System, only : napt
      use m_Crystal_Structure, only : op
      integer, intent(in) :: n, ia
      integer, intent(out) :: ja
      real(kind=DP), intent(in) :: ei(3)
      real(kind=DP), intent(out) :: ej(3)

      ! local variables
      integer :: i

      ja = napt(ia,n)

      do i=1,3
         ej(i) = sum(op(i,1:3,n)*ei(1:3))
      end do

    end subroutine rotate_eigenvector

  end subroutine m_Phonon_det_irr_rep


  subroutine m_Phonon_calc_static_dielectric()
    use m_Crystal_Structure, only : univol
    implicit none

    real(DP), parameter :: ha2ev = Hartree ! 1 Ha = 27.211383411 eV

    ! local variables
    integer :: im,i,na,ia,j
    real(DP) :: zsum
    real(DP) :: coef

    if(sw_lo_to_splitting == ON) return
    if(nqvec/=1) return

    allocate(zeff_mode(3,nmodes))
    do im=1,nmodes
       do i=1,3
          zsum = 0.d0
          j=0
          do na=1,natm
             do ia=1,3
                j=j+1
                zsum = zsum + zeff(i,ia,na)*dble(modes(j,im,nqvec))*msqinv(j)
             end do
          end do
          zeff_mode(i,im) = zsum
       end do
    end do

    coef = 16.d0*atan(1.d0)/univol
    lattice_dielectric(1:3,1:3) = 0.d0

    if (printable) then
       write(nfout,*) "*** Contribution to diagonal sum of dielectric tensor ***"
       write(nfout,'(A,15X,A)') '  mode   Repr. Active            Freq [cm-1]', &
            &                  ' Contrib.'
    endif

    do im=1,nmodes
       if(omega(im,nqvec)>1.d-10) then
          do j=1,3
             do i=1,3
                lattice_dielectric(i,j) = lattice_dielectric(i,j) +  &
                     & zeff_mode(i,im)*zeff_mode(j,im)/omega(im,nqvec)**2
             end do
          end do
          if(printable) write(nfout,'("im=",i6,1x,a3,1x,a4,2(1x,f25.10))')  &
               &  im,mode_irr_rep(im,nqvec),mode_active(im,nqvec) &
               & ,omega(im,nqvec)*ha2ev *ev2cminv &
               & ,coef*sum(zeff_mode(1:3,im)**2)/omega(im,nqvec)**2
       end if
    end do
    lattice_dielectric(1:3,1:3) = coef*lattice_dielectric(1:3,1:3)

    static_dielectric(1:3,1:3) =  &
         & lattice_dielectric(1:3,1:3) + dielectric(1:3,1:3)

  end subroutine m_Phonon_calc_static_dielectric


  subroutine m_Phonon_write_vib_modes()
    use m_Files, only : m_Files_open_nfmode, nfmode
    use m_Ionic_System, only : cps,imdtyp
    use m_Crystal_Structure, only : altv
    implicit none

    ! local variabes
!    real(DP), parameter :: amu = AMU ! 1 amu = 1822.88847974 Me
    real(DP), parameter :: ha2ev = Hartree ! 1 Ha = 27.211383411 eV
!    real(DP), parameter :: ev2cminv = 1.d0/1.23984185d-4 ! 1 cm^-1 = 1.233984185e-4 eV
    integer :: i,j,istart,iq
    real(DP) :: omega_ev,omega_cm
    real(DP) :: zamu(3),znorm
    real(DP) :: x,hbarW,free_e
    integer :: icount,jcount

    call  m_Files_open_nfmode()

    if(mype /= 0) return

    if ( sw_phonon_band_unfolding == ON ) then
       call m_Phonon_Unfold_spectr_weight
       write(nfmode,'(1x,"--- supercell lattice vectors ---")')
       do i=1,3
          write(nfmode,'(1x,3(f14.10,1x))') altv(1:3,i)
       end do
       write(nfmode,'(1x,"--- reference cell lattice vectors ---")')
       do i=1,3
          write(nfmode,'(1x,3(f14.10,1x))') altv_refcell(1:3,i)
       end do
    else
       write(nfmode,'(1x,"--- primitive lattice vectors ---")')
       do i=1,3
          write(nfmode,'(1x,3(f14.10,1x))') altv(1:3,i)
       end do
    endif

    write(nfmode,'(1x,"--- Equilibrium position and mass of each atom---")')
    write(nfmode,'(1x,"Natom=",i5)') natm
    do i=1,natm
       write(nfmode,'(1x,i5,3(1x,f14.10),1x,f14.5,1x,a4,i3)') &
            &                      i,cps(i,1:3),ionic_mass(i),speciesname(ityp(i)), &
            &                      atom_key(i)
    end do
    write(nfmode,'(1x,"--- Vibrational modes ---")')
    if(phonon_method == PHONON_GAMMA) then
       write(nfmode,'(1x,"Nmode=",i5," Natom=",i5)') nmodes,natm
    else
       write(nfmode,'(1x,"Nmode=",i5," Natom=",i5," Nqvec",i5)') nmodes,natm,nqvec
    end if

    if(phonon_method==PHONON_DOS .and. way_of_smearing==FERMI_DIRAC) free_e=0.d0

    do iq=1,nqvec
       if(phonon_method == PHONON_DOS) then
          write(nfmode,'(1x,"iq=",i5," q=(",f15.10,",",f15.10,",",f15.10,") (",f15.10,",",f15.10,",",f15.10,")",1x,f15.10)') &
               & iq,qvin(iq,1:3),qvec(iq,1:3),wght(iq)

!       else if(phonon_method /= PHONON_GAMMA) then
       else if(phonon_method == PHONON_BAND ) then
          if ( sw_phonon_band_unfolding == OFF ) then
             write(nfmode,'(1x,"iq=",i5," q=(",f15.10,",",f15.10,",",f15.10,") (",f15.10,",",f15.10,",",f15.10,")")') &
                  & iq,qvin(iq,1:3),qvec(iq,1:3)
          else
             write(nfmode,'(1x,"iq=",i5," q=(",f15.10,",",f15.10,",",f15.10,") (",f15.10,",",f15.10,",",f15.10,")")') &
                  & iq, qvin_refcell(iq,1:3), qvec(iq,1:3)
          endif
       end if

       do i=1,nmodes
          omega_ev = omega(i,iq)*ha2ev
          omega_cm = omega_ev*ev2cminv
          if ( sw_phonon_band_unfolding == ON ) then
             write(nfmode,'(1x,"n= ",i5, 1x, a3, 1x, a4, 5x,"weight=", e15.8, &
                  &                  /5x,"hbarW= ",e15.8," Ha = ",e15.8, &
                  &                " eV; nu= ",e15.8," cm^-1")') &
                  &                  i, mode_irr_rep(i,iq), mode_active(i,iq), &
                  &                  unfolding_weight(i,iq), &
                  &                  omega(i,iq), omega_ev, omega_cm
          else
             write(nfmode,'(1x,"n= ",i5, 1x, a3, 1x, a4, &
                  &                  /5x,"hbarW= ",e15.8," Ha = ",e15.8,&
                  &                  " eV; nu= ",e15.8," cm^-1")') &
                  &                  i, mode_irr_rep(i,iq), mode_active(i,iq), &
                  &                  omega(i,iq), omega_ev, omega_cm
          endif

          jcount=0
          do j=1,natm
             if(imdtyp(j)==OFF)cycle
             jcount = jcount+1
             istart = 3*(jcount-1)+1
             write(nfmode,'(1x,i5,1x,3f14.10)') j,dble(modes(istart:istart+2,i,iq))
          end do
          if(nqvec>1) then
             jcount=0
             do j=1,natm
                if(imdtyp(j)==OFF)cycle
                jcount = jcount+1
                istart = 3*(jcount-1)+1
                write(nfmode,'(1x,i5,1x,3f14.10)') j,dimag(modes(istart:istart+2,i,iq))
             end do
          end if
          if(phonon_method == PHONON_GAMMA.and.sw_lattice_dielectric_tensor == ON) then
             write(nfmode,'(3x,"Mode effective charge and its norm:")')
             zamu(1:3)=zeff_mode(1:3,i)*sqrt(amu)
             znorm = sqrt(sum(zamu(1:3)**2))
             write(nfmode,'(3x,"Z=",3(f14.10,1x),"Norm=",f14.10)') zamu(1:3),znorm
          end if
          if(phonon_method == PHONON_DOS.and.way_of_smearing==FERMI_DIRAC)then
             hbarW = omega(i,iq)
             x = hbarW/width
             if(x.gt.1.d-10) free_e = free_e + wght(iq)*(0.5d0*hbarW+width*log(1.d0-exp(-x)))
          endif
       end do
    end do
    if(phonon_method == PHONON_GAMMA.and.sw_lattice_dielectric_tensor == ON) then
       write(nfmode,'(1x,"--- Lattice and static dielectric tensors ---")')
       do i=1,3
          write(nfmode,'(1x,"[",3(1x,f10.4)," ]",2x,"[",3(1x,f10.4)," ]")') &
               & lattice_dielectric(i,1:3), static_dielectric(i,1:3)
       end do
    end if
    if(way_of_smearing==FERMI_DIRAC.and.phonon_method==PHONON_DOS.and.printable)then
       write(nfout,'(a,f15.10,a)') '!** Free energy due to lattice vibration : ',free_e,' hartree/cell'
    endif

    return
  end subroutine m_Phonon_write_vib_modes

  subroutine m_Phonon_write_epsilon()
    use m_Crystal_Structure, only : univol
    use m_Files, only : m_Files_open_nfepsilon, nfepsilon
    implicit none

    ! local variables
    real(DP), parameter :: ha2ev = Hartree ! 1 Ha = 27.211383411 eV
    integer :: n,im,i,j
    real(kind=DP) :: e,de,e2,e3,e4
    real(kind=DP) :: coef
    real(kind=DP) :: epsilon(3,3,2,0:num_division)

    if(nqvec /= 1) return

    coef = 16.d0*atan(1.d0)/univol
    de = (max_energy-min_energy)/num_division

    do n=0,num_division
       e = min_energy + de * n
       epsilon(1:3,1:3,1:2,n) = 0.d0
       do im=1,nmodes
          if(omega(im,nqvec)>1.d-5) then
             e2 = omega(im,nqvec)**2-e**2
             e3 = e*omega(im,nqvec)*damping_factor
             e4 = e2**2+e3**2
             if(abs(damping_factor) < 1.d-10) then
                if(abs(omega(im,nqvec)-e) < 1.d-10) then
                   if(e2 > 0.d0) then
                      e2 = 1.d-10
                   else
                      e2 = -1.d-10
                   end if
                end if
                do j=1,3
                   do i=1,j
                      epsilon(i,j,1,n) = epsilon(i,j,1,n) + zeff_mode(i,im)*zeff_mode(j,im)/e2
                      epsilon(i,j,2,n) = 0.d0
                   end do
                end do
             else
                do j=1,3
                   do i=1,j
                      epsilon(i,j,1,n) = epsilon(i,j,1,n) + zeff_mode(i,im)*zeff_mode(j,im)*e2/e4
                      epsilon(i,j,2,n) = epsilon(i,j,2,n) + zeff_mode(i,im)*zeff_mode(j,im)*e3/e4
                   end do
                end do
             end if
          end if
       end do
       epsilon(1:3,1:3,1,n)=dielectric(1:3,1:3)+coef*epsilon(1:3,1:3,1,n)
       epsilon(1:3,1:3,2,n)=coef*epsilon(1:3,1:3,2,n)
    end do

    call m_Files_open_nfepsilon()
    if(mype /= 0) return
    write(nfepsilon,'("Energy(eV) E1xx E1yy E1zz E1yz E1zx E1xy E2xx E2yy E2zz E2yz E2zx E2xy")')
    do n=0,num_division
       e = (min_energy + de * n) * ha2ev
       write(nfepsilon,'(13(1x,e25.10))') &
            & e, ((epsilon(i,i,j,n),i=1,3),epsilon(2,3,j,n),epsilon(1,3,j,n),epsilon(1,2,j,n),j=1,2)
    end do

  end subroutine m_Phonon_write_epsilon

  subroutine m_Phonon_set_supercell
    use m_Crystal_Structure, only: m_CS_supercell_on
    integer :: n1,n2,n3
    n1 = nlat(1)
    n2 = nlat(2)
    n3 = nlat(3)
    if(n1==0.and.n2==0.and.n3==0) then
       n1 = 1; n2 = 1; n3 = 1
       call m_CS_supercell_on(nfout,1,n1,n2,n3) ! 1: PRIMITIVE
    else
       call m_CS_supercell_on(nfout,0,n1,n2,n3) ! 0: BRAVAIS
    end if
  end subroutine m_Phonon_set_supercell

  subroutine m_Phonon_calc_dos()
    integer :: idim = 3
    integer :: iq,im,ie
    integer :: newin
    real(kind=DP) :: emin,emax
    !  real(kind=DP), dimension(np2,nmodes) :: omega2
    !  real(kind=DP), dimension(np0) :: eawk, cdwk, cswk
    real(kind=DP), allocatable, dimension(:,:) :: omega2 ! dim(np2,nmodes)      <ASMS>
    real(kind=DP), allocatable, dimension(:)   :: eawk, cdwk, cswk ! dim(np0)   <ASMS>
    real(kind=DP), allocatable, dimension(:,:,:) :: cdos,cind ! dim(np2,nmodes,0:newin)
    real(kind=DP), allocatable, dimension(:) :: e,dos,sumdos ! dim(0:newin)

    allocate(omega2(np2,nmodes))                     ! <ASMS>
    allocate(eawk(np0),cdwk(np0),cswk(np0))          ! <ASMS>
    eawk=0.d0
    cdwk=0.d0
    cswk=0.d0

    do iq=1,np2
       do im=1,nmodes
          omega2(iq,im) = omega(im,iq) * 1.0d3
       end do
    end do

    emin = minval(omega2) - 0.0005 ! (mhartree)
    emax = maxval(omega2) + 0.0005 ! (mhartree)
    newin = (emax - emin)/(DeltaE_dos*1.0d3) + 1

    allocate(cdos(np2,nmodes,0:newin)); cdos = 0.d0
    allocate(cind(np2,nmodes,0:newin)); cind = 0.d0

    allocate(e(0:newin))
    allocate(dos(0:newin)); dos=0.d0
    allocate(sumdos(0:newin)); sumdos=0.d0
    e(0:newin) = (/(emin + (DeltaE_dos*1.0d3)*ie,ie=0,newin)/)

    call nstt3i(idim,newin,e,nxyz_tetra(1),nxyz_tetra(2),nxyz_tetra(3) &
         &  ,np2,np2,nmodes,omega2 &
         &  ,ip20,np0,eawk,cdwk,cswk,np2,nmodes,cdos,cind,width_tetra )

    do iq=1,np2
       do im=1,nmodes
          do ie=0,newin
             dos(ie)    = dos(ie)    + cdos(iq,im,ie)
             sumdos(ie) = sumdos(ie) + cind(iq,im,ie)
          end do
       end do
    end do
    deallocate(cdos)
    deallocate(cind)

    call write_dos(newin,e,dos,sumdos)

    deallocate(e)
    deallocate(dos)
    deallocate(sumdos)
  contains

    subroutine write_dos(newin,e,dos,sumdos)
      use m_Files, only : m_Files_open_nfphdos, nfphdos
      integer, intent(in) :: newin
      real(kind=DP), intent(in) :: e(0:newin), dos(0:newin), sumdos(0:newin)
      integer :: ie
      real(DP), parameter :: ha2ev = Hartree
      real(DP) :: eev, ecminv, dosev, doscminv

      if(mype /= 0) return
      call m_Files_open_nfphdos()

      write(nfphdos,'("# Index Omega(mHa) Omega(eV) Omega(cm-1) DOS(States/mHa) DOS(States/eV) DOS(States/cm-1) IntDOS(States)")')
      do ie = 0, newin
         eev = e(ie)*ha2ev/1.0d3  ! Hartree
         ecminv = eev*ev2cminv
         dosev = dos(ie)/ha2ev*1.0d3 ! Hartree
         doscminv = dosev/ev2cminv
         write(nfphdos,'(i7,7(1x,f15.8))') ie,e(ie), eev, ecminv, dos(ie), dosev, doscminv, sumdos(ie)
      end do
    end subroutine write_dos
  end subroutine m_Phonon_calc_dos

  subroutine m_Phonon_read_strain_forces
    use m_Files, only : nfout, nfstrfrc, m_Files_open_nfstrfrc
    integer :: i, ia, num_force, istr, ipm
    real(kind=DP) :: str, force(3,natm,6,2), strain(6,2), eps(3,3)
    logical :: present(6)

    strain = 0.d0
    present = .false.

    allocate(strfrc(3,natm,6))
    strfrc = 0.d0

    call m_Files_open_nfstrfrc()
    read(nfstrfrc,*) num_force
    do i=1,num_force
       read(nfstrfrc,*) istr, str
       if(str>0.d0) then
          ipm = 1
       else
          ipm = 2
       end if
       strain(istr,ipm) = str
       do ia=1,natm
          read(nfstrfrc,*) force(1:3,ia,istr,ipm)
       end do
    end do

    do istr=1,6
       if(strain(istr,1) > 1.d-8 .and. strain(istr,2) < -1.d-8 &
            & .and. abs(abs(strain(istr,1))-abs(strain(istr,2))) < 1.d-10) then
          present(istr) = .true.
       end if
    end do

    write(nfout,*) 'present strain : ',present

    do istr=1,6
       if(.not.present(istr)) cycle
       str = strain(istr,1)-strain(istr,2)
       eps = 0.d0
       if(istr==1) then
          eps(1,1) = 0.5d0
       else if(istr==2) then
          eps(2,2) = 0.5d0
       else if(istr==3) then
          eps(3,3) = 0.5d0
       else if(istr==4) then
          eps(2,3) = 0.25d0
          eps(3,2) = 0.25d0
       else if(istr==5) then
          eps(3,1) = 0.25d0
          eps(1,3) = 0.25d0
       else if(istr==6) then
          eps(1,2) = 0.25d0
          eps(2,1) = 0.25d0
       end if
       do ia=1,natm
          strfrc(1:3,ia,istr) = (force(1:3,ia,istr,1)-force(1:3,ia,istr,2))/str &
               & + matmul(eps,force(1:3,ia,istr,1)+force(1:3,ia,istr,2))
       end do
    end do

    write(nfout,*) '=== Strain-force coupling constants ==='
    do istr=1,6
       if(.not.present(istr)) cycle
       write(nfout,'("Strain : ",i1)') istr
       do ia=1,natm
          Write(nfout,'(i5,3(1x,f15.5))') ia,strfrc(1:3,ia,istr)
       end do
    end do

  end subroutine m_Phonon_read_strain_forces


  subroutine m_Phonon_calc_int_strain_piezo
    use m_Crystal_Structure, only : univol
    implicit none

    ! local variables
    integer :: im,i,na,ia,j
    real(DP) :: zsum,csum
    real(DP) :: coef

    if(sw_lo_to_splitting == ON) return
    if(nqvec/=1) return

    allocate(zeff_mode(3,nmodes))
    allocate(strphcpl(6,nmodes))
    do im=1,nmodes
       do i=1,3
          zsum = 0.d0
          j=0
          do na=1,natm
             do ia=1,3
                j=j+1
                zsum = zsum + zeff(i,ia,na)*dble(modes(j,im,nqvec))*msqinv(j)
             end do
          end do
          zeff_mode(i,im) = zsum
       end do
       do i=1,6
          csum = 0.d0
          j=0
          do na=1,natm
             do ia=1,3
                j=j+1
                csum = csum + strfrc(ia,na,i)*dble(modes(j,im,nqvec))*msqinv(j)
             end do
          end do
          strphcpl(i,im) = csum
       end do
    end do

    coef = 1.d0/univol
    piezoelectric(1:3,1:6) = 0.d0
    do j=1,6
       do i=1,3
          write(nfout,'("Direction:",i1," Strain:",i1)') i,j
          do im=1,nmodes
             if(omega(im,nqvec)>1.d-10) then
                piezoelectric(i,j) = piezoelectric(i,j) +  &
                     & zeff_mode(i,im)*strphcpl(j,im)/omega(im,nqvec)**2
                if(printable) write(nfout,'("im=",i4,1x,a3,1x,a4,3(1x,f25.10))')  &
                     &  im,mode_irr_rep(im,nqvec),mode_active(im,nqvec) &
!                     & ,omega(im,nqvec)*Hartree/1.23984185d-4 &
                     & ,omega(im,nqvec)*Hartree * ev2cminv &
                     & ,coef*zeff_mode(i,im)*strphcpl(j,im)/omega(im,nqvec)**2
             end if
          end do
       end do
    end do
    piezoelectric(1:3,1:6) = coef*piezoelectric(1:3,1:6)

    if(printable) then
       write(nfout,*) '=== Internal-strain piezoelectric tensor (a.u.) ==='
       do j=1,6
          write(nfout,'(i5,3(1x,f18.10))') j, piezoelectric(1:3,j)
       end do
       write(nfout,*) '=== Internal-strain piezoelectric tensor (C/m^2) ==='
       do j=1,6
          write(nfout,'(i5,3(1x,f18.10))') j, piezoelectric(1:3,j)*UNIT_PIEZO_CONST
       end do
    end if

  end subroutine m_Phonon_calc_int_strain_piezo

  subroutine m_Phonon_Unfold_spectr_weight
    use m_Ionic_System,          only : imdtyp
    use m_Const_Parameters,      only : CMPLDP, zi

    integer :: iq, i, ixyz, ig, ia
    integer :: jcount, istart
    integer :: nx, ny, nz
    real(kind=DP) :: csum1, csum2, c1, qq(3)
    complex(kind=CMPLDP) :: zsum1, zsum2, zfac

    integer :: nxmax, nymax, nzmax

    call m_Phonon_set_GVec_Flag_refcell
    allocate( unfolding_weight(nmodes,nqvec) )

    nxmax = nmax_G_phonon(1)
    nymax = nmax_G_phonon(2)
    nzmax = nmax_G_phonon(3)

    Do iq=1, nqvec
       Do i=1, nmodes
          csum1 = 0.0d0;    csum2 = 0.0d0
          Do ixyz=1, 3
             Do nx=-nxmax, nxmax
                Do ny=-nymax, nymax
                   Do nz=-nzmax, nzmax
!                      qq(1:3) = nx *rltv(1:3,1) +ny *rltv(1:3,2) +nz *rltv(1:3,3)
                      qq(1) = nx;     qq(2) = ny;     qq(3) = nz

                      jcount = 0;     zsum1 = 0.0d0

                      Do ia=1, natm
                         if (imdtyp(ia)==OFF) cycle
                         jcount = jcount+1
                         istart = 3*(jcount-1)+1

                         c1 = qq(1) *pos(ia,1) +qq(2) *pos(ia,2) +qq(3) *pos(ia,3)
!                         c1 = qq(1) *cps(ia,1) +qq(2) *cps(ia,2) +qq(3) *cps(ia,3)
                         zfac = exp( -zi *c1 *PAI2 )

                         zsum1 = zsum1 +zfac *modes( istart+ixyz-1, i, iq )
                      End Do
                      csum1 = csum1 +conjg(zsum1) *zsum1
                      if ( phonon_GVec_on_refcell(nx,ny,nz) == 1 ) then
                         csum2 = csum2 +conjg(zsum1) *zsum1
                      endif
                   End Do
                End Do
             End Do
          End Do
          unfolding_weight(i,iq) = csum2 /csum1
       End Do
    End Do

  end subroutine m_Phonon_Unfold_spectr_weight

  subroutine m_Phonon_set_GVec_Flag_refcell
    integer :: k, n1, n2, n3, nx, ny, nz
    real(kind=DP) :: vec1(3), cx, cy, cz, delta

    n1 = nmax_G_phonon(1)
    n2 = nmax_G_phonon(2)
    n3 = nmax_G_phonon(3)

    allocate( Phonon_GVec_on_refcell(-n1:n1,-n2:n2,-n3:n3) )
    Phonon_Gvec_on_refcell = 0

    delta = tolerance_ph_Gvec_matching

    Do nx=-n1, n1
       Do ny=-n2, n2
          Do nz=-n3, n3
             vec1(1) = nx *rltv(1,1) +ny *rltv(1,2) +nz *rltv(1,3)
             vec1(2) = nx *rltv(2,1) +ny *rltv(2,2) +nz *rltv(2,3)
             vec1(3) = nx *rltv(3,1) +ny *rltv(3,2) +nz *rltv(3,3)

             cx = 0.0d0; cy = 0.0d0; cz = 0.0d0
             Do k=1, 3
                cx = cx +altv_refcell(k,1) *vec1(k)
                cy = cy +altv_refcell(k,2) *vec1(k)
                cz = cz +altv_refcell(k,3) *vec1(k)
             End Do

             cx = cx /PAI2;     cy = cy /PAI2;      cz = cz /PAI2
             cx = abs( cx -nint(cx) );
             cy = abs( cy -nint(cy) );
             cz = abs( cz -nint(cz) );

             if ( cx < delta .and. cy < delta .and. cz < delta ) then
                phonon_GVec_on_refcell(nx,ny,nz) = 1
             endif
          ENd Do
       End Do
    ENd Do
  end subroutine m_Phonon_set_GVec_Flag_refcell

end module m_Phonon
