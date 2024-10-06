!=======================================================================
!
!  PROGRAM  PHASE/0 2014.01 (rev.375)
!
!  "First-principles Electronic Structure Calculation Program"
!
!  MODULE: m_Const_Parameters
!
!  AUTHOR(S): T. Yamasaki and K. Mae   August/20/2003
!  
!  Contact address :  Phase System Consortium
!                     E-mail: phase_system@nims.go.jp URL https://azuma.nims.go.jp
!  
!
!
!=======================================================================
!
!   The original version of this set of the computer programs "PHASE" was developed by 
!  the members of the Theory Group of Joint Research Center for Atom Technology 
!  (JRCAT), based in Tsukuba, in the period 1993-2001.  
!   Since 2002, this program set had been intensively developed as a part of the following 
!  national projects supported by the Ministry of Education, Culture, Sports, Science and 
!  Technology (MEXT) of Japan; "Frontier Simulation Software for Industrial Science 
!  (FSIS)" from 2002 to 2005, "Revolutionary Simulation Software (RSS21)" from 2006 to 
!  2008. "Research and Development of Innovative Simulation Software (RISS)" from 2008 
!  to 2013. These projects is lead by the Center for Research on Innovative Simulation 
!  Software (CISS), the Institute of Industrial Science (IIS), the University of Tokyo.
!   Since 2013, this program set has been further developed centering on PHASE System 
!  Consortium. 
!   The activity of development of this program set has been supervised by Takahisa Ohno.
!
!***************************************************************
module m_Const_Parameters
! $Id: m_Const_Parameters.f90 817 2009-05-28 05:21:30Z yamasaki $
  implicit none

  integer, parameter :: DRIVER_GENERAL=0, DRIVER_CONSTRAINT=1, DRIVER_NEB=2, DRIVER_PHONON=3
  integer, parameter :: SCF    = 0, EK   = -1
  integer, parameter :: NODATA = 0, PUCV = 0, CARTS = 1, BUCS = 2
  integer, parameter :: CRDTYP = BUCS
  integer, parameter :: NONAME = 0, MESH = 1, SKPS_DIRECT_IN = 4 &
       &     , GAMMA = 3, FILE = 2, MONKHORST_PACK = 5
  integer, parameter :: GAMMA_base_symmetrization = 6
  integer, parameter :: DIRECTIN = 2, FROM_ENDPOINTS = 3, PROPORTIONAL = 4
  integer, parameter :: BEGIN = 1, END = 2
  integer, parameter :: ITER01 = 0, START = 0, FINISH = 3 &
       &     , DFINIS = 4, FLUSH = 5, ITERATIVE = 1
  integer, parameter :: TAG_FORMAT = 0, TAG_LINE = 1, TABLE = 2
  integer, parameter :: DP = kind(1.d0)
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: CMPLDP = kind((1.d0, 0.d0))
  integer, parameter :: CMPLSP = kind((1.0 , 0.0 ))
  integer, parameter :: irsize = 8
  real(kind=DP), parameter :: PAI  = 3.141592653589793238462D0
  real(kind=DP), parameter :: PAI2 = PAI*2.d0, PAI4 = PAI*4.d0
  integer, parameter :: SD  = -1, SDPRC = 0, MSD = 1, MSDPRC = 2 &
       &        , CG  = 3, CGPRC = 4, RMM = 5, RMMPRC = 6, RMM2P = 17,RMM2PPRC = 18 &
       &	, RMM2 = 7, RMM2PRC = 8, SDLM = 9, SDLMPRC = 10 &
       &        , MSDLM = 11, MSDLMPRC = 12, bbCG = 13, bbCGPRC = 14 &
       &        , eazyCG = 15, eazyCGPRC = 16 &
       &        , default_sd2cg = 30&
       &        , SUBMAT = 51, MATRIXDIAGON = 53, DAVIDSON = 55
  integer, parameter :: lmSD = SDLM, lmMSD = MSDLM, lmCG = CG &
       &        , lmRMM2 = RMM2, lmRMM = RMM, lmeazyCG = eazyCG&
       &        , RMM3 = RMM, RMM3PRC = RMMPRC
  integer, parameter :: SURFACE = 0, WHOLE_BZ = 1, SIMPLE_CUBIC = 2 &
       &        , BCC = 3,  FCC = 4, DIAMOND = 5, HEXAGONAL = 6, ORTHORHOMBIC = 7 &
       &        , HCP = 33,RUTILE = 8, C2v_SURFACE = 20, Csy_SURFACE = 21 &
       &        , GENERAL = 100, GENERAL_LARGER = 101 &
       &        , HEX1fold = 30, HEX2fold = 31, HEX3fold = 32, TRIGONAL = -11
  integer, parameter :: REGULAR_INTERVALS = 1, EAZY = 1, BY_ATOMIC_POSITIONS = 2, BY_ATOMS = 2
  integer, parameter :: TRIGONAL_lattice = -1, HEXAGONAL_lattice = 0, PRIMITIVE_lattice = 1 &
       &              , FCC_lattice = 2, BCC_lattice = 3, BottomCenteredCubic_lattice = 4
  integer, parameter :: BRAVAIS = 0, PRIMITIVE = 1
  integer, parameter :: SMALL = 1, MEDIUM = 2, LARGE = 3
  integer, parameter :: OLD = 0, NEW_ = 1, UNKNOWN = 2, NEXT = 1
  integer, parameter :: FORMATTED = 0, UNFORMATTED = 1
  integer, parameter :: ON = 1, OFF = 0
  integer, parameter :: GRID = 2
  integer, parameter :: YES = 1, NO = 0
  integer, parameter :: PARA = 0, ANTIFERRO = 1, FERRO = 2, NONMAG = 0
  integer, parameter :: check_file_name_on  = 1
  integer, parameter :: check_file_name_off = 0
  integer, parameter :: by_matrix_diagon    = 1
  integer, parameter :: by_random_numbers   = 0
  integer, parameter :: by_pseudo_atomic_orbitals = 3
  integer, parameter :: Gauss_distrib_func  = 1
  integer, parameter :: PROJECTOR = 1, WAVEFUNCTION = 2, MULLIKEN = 3
  integer, parameter :: SPHERICAL_HARMONICS = 1, ATOMIC_ORBITAL = 2, WANNIER_FUNCTION = 3
  integer, parameter :: VERY_NARROW         = 3
  integer, parameter :: from_PSEUDOPOTENTIAL_FILE = 4
  integer, parameter :: from_wave_functions = 2
  integer, parameter :: INITIAL = 0, CONTINUATION = 1, FIXED_CHARGE = 2 &
       & , FIXED_CHARGE_CONTINUATION = 3, AUTOMATIC = -1, PREPARATION_ONLY = -2 &
       & , FIXED_CHARGE_AUTOMATIC    = -3 
  integer, parameter :: ALL_AT_ONCE = 0, ONE_BY_ONE = 1
  integer, parameter :: MANUAL = -2
  integer, parameter :: Partial_Core_Charge = 0, Valence_plus_PC_Charge = 1
  integer, parameter :: Core_Charge = 2
  integer, parameter :: GGA = 1, LDA  = 0
  integer, parameter :: UP  = 1, DOWN = 2
  integer, parameter :: SKIP = 0, EXECUT = 1
  integer, parameter :: SOFTPART = 1, HARDPART = 2
  integer, parameter :: WITHOUTTAG = 0, WITHTAG = 1
  integer, parameter :: INVERSE = 1, DIRECT = 2
  integer, parameter :: ORTHONORMALIZATION = 0, ORTHOGONALIZATION = 1, NORMALIZATION = 2
  integer, parameter :: OVER = 0, UNDER = 1
  integer, parameter :: OUTER = 1, INNER = -1
  integer, parameter :: ALL_BANDS = 0, OTHER_BANDS = 1, SAME_BAND = 2
  complex(kind=CMPLDP), parameter :: zi = (0.d0,1.d0)
  integer, parameter :: MDOLD = -1, MDSMPL = -2
  integer, parameter :: PARABOLIC = 0, MP = 1, TETRAHEDRON = 2, COLD = 3
  integer, parameter :: TEMPERATURE_CONTROL = -1, VERLET = 1&
       &, QUENCHED_MD    = 2,       NORMAL_MODE_ANALYSIS = 3 &
       &, HYPERPLANE_ADAPTIVE_COORDINATE = 4 &
       &, GDIIS  = 5 &
       &, BLUEMOON = -5, QUENCHED_CONSTRAINT = -6 &
       &, HAC = HYPERPLANE_ADAPTIVE_COORDINATE &
       &, T_CONTROL = TEMPERATURE_CONTROL &
       &, PHONON_FORCE = 10, CG_STROPT = 11 &
       &, SD_MD = 12, STEEPEST_DESCENT = SD_MD, MC = 13 &
       &, DAMPED_MD = 14, VELOCITY_SCALING = 15
  integer, parameter :: ORDINA = 0, CNSTRA = 1
  integer, parameter :: ALDOS = 1, LAYERDOS = 2
  integer, parameter :: FIX    = 0, RELAX  = 1, BONDLENGTH_FIX = 10 &
       &, COG_FIX = 101, COG_CNTR = 101, HEAT_BATH = 1000, FIX_IN_A_PLANE = 5 &
       &, BONDLENGTH_FIX_1 = 11, BONDLENGTH_FIX_2 = 12, COG_FIX_L = 102 &
       &, FIX_HBATH = FIX - 1, RELAX_HBATH = RELAX - 1 &
       &, FIXED_NORMAL_HYPERVECTOR = 13, NUDGED_ELASTIC_BAND_METHOD = 14
  integer, parameter :: NOSE = 2000, NOSE_HOOVER = 1000
  integer, parameter :: IINCRE_CRITICAL = 1000
  real(kind=DP), parameter :: UMICRO = 1.d-6
  real(kind=DP), parameter :: LIMIT_1DSEARCH = 1.d+10
  integer, parameter :: EXC_ONLY = 0, VXC_AND_EXC = 1, STRESS_ = 2
  integer, parameter :: SIMPLE = 1, BROYD1 = 2, BROYD2 = 3, DFP = 4 &
       &, PULAY = 5, ANEW = 1, RENEW = 2
  integer, parameter :: NORMCONSERVATION = 0, VANDERBILT_TYPE = 1, VDB = VANDERBILT_TYPE
  real(kind=DP),parameter :: BOHR = 0.529177d0
  real(kind=DP),parameter :: Hartree = 27.211396311
  integer, parameter ::      HR_ENERGY_UNIT =1, EV_ENERGY_UNIT = 2
  real(kind=DP),parameter :: DELTA = 1.d-60, DELTA10 = 1.d-10, DELTA07 = 1.d-07
  real(kind=DP),parameter :: DELTAevdff = 1.d-10
  real(kind=DP),parameter :: DELTA_FermiSearchRange = 1.d0
  integer,parameter ::       ELECTRON = 0, POSITRON = 1
  integer,parameter ::       TOTAL = 0
  integer,parameter ::       MAPPED = 1, NOTMAPPED = 0
  integer,parameter ::       WHOLE = 0,  INITIALLY = 1

  integer,parameter :: NOPOLAR = 0, POLARIZATION = 1, EFFECTIVE_CHARGE = 2, PIEZOELECTRIC_CONST = 3
  integer,parameter :: UNIFIED = 0, PARTITIONED = 1
  integer,parameter :: HOUSEHOLDER = 0, DIVIDEandCONQUER = 1

  integer,parameter :: PHONON_GAMMA = 0, PHONON_BAND = 1, PHONON_DOS = 2

  integer, parameter :: len_tag_ionic_system = 12
  character(len=len_tag_ionic_system) :: tag_ionic_system
  data tag_ionic_system/'Ionic System'/

  integer, parameter :: len_tag_iteration = 48
  character(len=len_tag_iteration) :: tag_iteration
  data tag_iteration/'iteration, iteration_ionic, iteration_electronic'/

  integer, parameter :: len_tag_iters_and_nk = 50
  character(len=len_tag_iters_and_nk) :: tag_iters_and_nk
  data tag_iters_and_nk/'iteration, iteration_electronic, nk_in_the_process'/

  integer, parameter :: len_tag_npes_etc = 22
  character(len=len_tag_npes_etc) :: tag_npes_etc
  data tag_npes_etc /'npes, nrank_e, nrank_k'/

  integer, parameter :: len_tag_total_energy = 12
  character(len=len_tag_total_energy) :: tag_total_energy
  data tag_total_energy/'Total Energy'/

  ! -- isolver
  integer,parameter  :: len_tag_isolver = len("isolver")
  character(len=len_tag_isolver) :: tag_isolver = "isolver"

  real(kind=DP) :: x
  real(kind=DP),parameter :: SmallestPositiveNumber = tiny(x)

  integer, parameter :: TOTAL_ENERGY = 1, MODIFIED_TOTAL_ENERGY = 2, BAND_ENERGY = 3

  integer, parameter :: BEFORE = 0, AFTER = 1
  integer, parameter :: varLINEAR = 0, varTANH = 1

  integer, parameter :: FORCE_CONVERGED = 2, CHARGE_CONVERGED = 1 &
       & , STRESS_CONVERGED = 3, EK_CONVERGED = 2

  integer, parameter :: BINARY = 3, VTK = 2, CUBE = 1, DENSITY_ONLY = 0, INTEGRATED = 0, SEPARATE = 1

  integer, parameter :: FMAXTAGLEN = 256
  integer, parameter :: FMAXVALLEN = 256
  integer, parameter :: FMAXUNITLEN = 20
  integer, parameter :: FMAXDEPTH = 16

  integer, parameter :: num_d6h = 24
  integer, parameter :: num_oh = 48
  character(5), dimension(num_d6h ), parameter :: d6h_symbol = (/ &
                                                & 'e    ', &	! 1
						& 'c6+  ', &	! 2
						& 'c3+  ', &	! 3
						& 'c2   ', &	! 4
						& 'c3-  ', &	! 5
						& 'c6-  ', &	! 6
						& 'c211 ', &	! 7
						& 'c221 ', &	! 8
						& 'c231 ', &	! 9
						& 'c212 ', &	! 10
						& 'c222 ', &	! 11
						& 'c232 ', &	! 12
						& 'ie   ', &	! 13
						& 'ic6+ ', &	! 14
						& 'ic3+ ', &	! 15
						& 'ic2  ', &	! 16
						& 'ic3- ', &	! 17
						& 'ic6- ', &	! 18
						& 'ic211', &	! 19
						& 'ic221', &	! 20
						& 'ic231', &	! 21
						& 'ic212', &	! 22
						& 'ic222', &	! 23
						& 'ic232'  /)	! 24


  character(5), dimension(num_oh), parameter :: oh_symbol = (/ &
        & 'e    ', 'c2x  ', 'c2y  ', 'c2z  ', &	! 1 2 3 4
	& 'c31+ ', 'c32+ ', 'c33+ ', 'c34+ ', &	! 5 6 7 8
	& 'c31- ', 'c32- ', 'c33- ', 'c34- ', &	! 9 10 11 12
	& 'c2a  ', 'c2b  ', 'c2c  ', 'c2d  ', &	! 13 14 15 16
	& 'c2e  ', 'c2f  ', 'c4x+ ', 'c4y+ ', &	! 17 18 19 20
	& 'c4z+ ', 'c4x- ', 'c4y- ', 'c4z- ', &	! 21 22 23 24
	& 'ie   ', 'ic2x ', 'ic2y ', 'ic2z ', &	! 25 26 27 28
        & 'ic31+', 'ic32+', 'ic33+', 'ic34+', &	! 29 30 21 32
	& 'ic31-', 'ic32-', 'ic33-', 'ic34-', &	! 33 34 35 36
	& 'ic2a ', 'ic2b ', 'ic2c ', 'ic2d ', &	! 37 38 39 40
        & 'ic2e ', 'ic2f ', 'ic4x+', 'ic4y+', &	! 41 42 43 44
        & 'ic4z+', 'ic4x-', 'ic4y-', 'ic4z-' /)	! 45 46 47 48
  ! --> m_input_interface
  integer, parameter :: TL_LINENO = 0
  integer, parameter :: TL_FROMCUR = 1
  integer, parameter :: TL_NEXT = 1
  integer, parameter :: TL_BACK = -1
  
  integer, parameter :: NOCONV = 0
  integer, parameter :: LOWER = 1
  integer, parameter :: UPPER = 2

  integer, parameter :: numdefaultunits = 9
  character(FMAXUNITLEN), dimension(numdefaultunits), parameter :: defaultunits = (/ &
       & 'hartree             ', 'bohr                ', 'au_time             ', 'au_mass             ', 'degree              ', &
       & 'k                   ', 'bohr/au_time        ', 'hartree/bohr        ', 'hartree/bohr3       '                     /)

  type unitsystem
      integer :: energy
      logical :: inpfg_energy
      integer :: length
      logical :: inpfg_length
      integer :: time
      logical :: inpfg_time
      integer :: mass
      logical :: inpfg_mass
      integer :: angle
      logical :: inpfg_angle
      integer :: temperature
      logical :: inpfg_temperature
      integer :: velocity
      logical :: inpfg_velocity
      integer :: force
      logical :: inpfg_force
      integer :: pressure
      logical :: inpfg_pressure
  end type unitsystem
  type(unitsystem), dimension(FMAXDEPTH) :: usys_file
  integer :: fblkdepth = 1
  ! <-- m_input_interface

  ! --> m_unit_conv 
real(DP), parameter :: CONST_EV = 1.6021773d-19
real(DP), parameter :: CONST_NA = 6.022137d+23
real(DP), parameter :: CONST_kB = 3.16676087d-6 ! = 1/11604.6996795594d0/27.21139615d0 K -> Hartree

real(DP), parameter :: UNIT_PIEZO_CONST = 57.214762537d0 ! C/m^2 = 1 a.u.

!!$real(DP), parameter :: CONST_PI = 3.14159265358979323846d0
real(DP), parameter :: CONST_PI = PAI

integer, parameter :: TYPE_ENERGY = 1
integer, parameter :: TYPE_LENGTH = 2
integer, parameter :: TYPE_TIME = 3
integer, parameter :: TYPE_MASS = 4
integer, parameter :: TYPE_ANGLE = 5
integer, parameter :: TYPE_TEMPERATURE = 6
integer, parameter :: TYPE_VELOCITY = 7
integer, parameter :: TYPE_FORCE = 8
integer, parameter :: TYPE_PRESSURE = 9

integer, parameter :: UNITNAMELEN = 16
type unitlist
    character(UNITNAMELEN) :: name
    real(kind=DP) :: factor
    integer ::       type
end type unitlist

integer, parameter :: BULK   = 1
integer, parameter :: DEFECT = 2

integer, parameter :: unit_list_size = 50

type(unitlist),parameter :: unit_1  = unitlist( 'ev                  ',    1.d0,                      TYPE_ENERGY )
type(unitlist),parameter :: unit_2  = unitlist( 'rydberg             ',   13.605698075d0,             TYPE_ENERGY )
type(unitlist),parameter :: unit_3  = unitlist( 'hartree             ',   27.21139615d0,              TYPE_ENERGY )
type(unitlist),parameter :: unit_4  = unitlist( 'j/mol               ',    1.d0/CONST_EV/CONST_NA,    TYPE_ENERGY )
type(unitlist),parameter :: unit_5  = unitlist( 'cal/mol             ',    4.184d0/CONST_EV/CONST_NA, TYPE_ENERGY )
type(unitlist),parameter :: unit_6  = unitlist( 'angstrom            ',    1.d0,                      TYPE_LENGTH )
type(unitlist),parameter :: unit_7  = unitlist( 'bohr                ',    0.5291772480d0,            TYPE_LENGTH )
type(unitlist),parameter :: unit_8  = unitlist( 'nm                  ',    1.d+1,                     TYPE_LENGTH )
type(unitlist),parameter :: unit_9  = unitlist( 's                   ',    1.d0,                      TYPE_TIME   )
type(unitlist),parameter :: unit_10 = unitlist( 'sec                 ',    1.d0,                      TYPE_TIME   )
type(unitlist),parameter :: unit_11 = unitlist( 'au_time             ',    2.418884327d-17,           TYPE_TIME   )
type(unitlist),parameter :: unit_12 = unitlist( 'fs                  ',    1.d-15,                    TYPE_TIME   )
type(unitlist),parameter :: unit_13 = unitlist( 'ps                  ',    1.d-12,                    TYPE_TIME   )
type(unitlist),parameter :: unit_14 = unitlist( 'ns                  ',    1.d-9,                     TYPE_TIME   )
type(unitlist),parameter :: unit_15 = unitlist( 'min                 ',   60.d0,                      TYPE_TIME   )
type(unitlist),parameter :: unit_16 = unitlist( 'hour                ', 3600.d0,                      TYPE_TIME   )
type(unitlist),parameter :: unit_17 = unitlist( 'day                 ',   24.d0*3600.d0,              TYPE_TIME   )
type(unitlist),parameter :: unit_18 = unitlist( 'kg                  ',    1.d0,                      TYPE_MASS   )
type(unitlist),parameter :: unit_19 = unitlist( 'g                   ',    1.d-3,                     TYPE_MASS   )
type(unitlist),parameter :: unit_20 = unitlist( 'au_mass             ',    9.1093897d-31,             TYPE_MASS   )
type(unitlist),parameter :: unit_21 = unitlist( 'atomic_mass         ',    1.66053d-27,               TYPE_MASS   )
type(unitlist),parameter :: unit_22 = unitlist( 'degree              ',    1.d0,                      TYPE_ANGLE  )
type(unitlist),parameter :: unit_23 = unitlist( 'radian              ',  180.d0/CONST_PI,             TYPE_ANGLE  )
type(unitlist),parameter :: unit_24 = unitlist( 'k                   ',    0.d0,                 TYPE_TEMPERATURE )
type(unitlist),parameter :: unit_25 = unitlist( 'centigrade          ',  273.15d0,               TYPE_TEMPERATURE )
type(unitlist),parameter :: unit_26 = unitlist( 'angstrom/s          ',    1.d0,                    TYPE_VELOCITY )
type(unitlist),parameter :: unit_27 = unitlist( 'angstrom/fs         ',    1.d0/1.d-15,             TYPE_VELOCITY )
type(unitlist),parameter :: unit_28 = unitlist( 'angstrom/au_time    ',    1.d0/2.418884327d-17,    TYPE_VELOCITY )
type(unitlist),parameter :: unit_29 = unitlist( 'bohr/fs             ',    0.5291772480d0/1.d-15,   TYPE_VELOCITY )
type(unitlist),parameter :: unit_30 = unitlist( 'bohr/au_time        ',    0.5291772480d0/2.418884327d-17, TYPE_VELOCITY )
type(unitlist),parameter :: unit_31 = unitlist( 'nm/fs               ',    1.d+1/1.d-15,                   TYPE_VELOCITY )
type(unitlist),parameter :: unit_32 = unitlist( 'nm/au_time          ',    1.d+1/2.418884327d-17,          TYPE_VELOCITY )
type(unitlist),parameter :: unit_33 = unitlist( 'ev/angstrom         ',    1.d0,                           TYPE_FORCE )
type(unitlist),parameter :: unit_34 = unitlist( 'ev/bohr             ',    1.d0/0.5291772480d0,            TYPE_FORCE )
type(unitlist),parameter :: unit_35 = unitlist( 'ev/nm               ',    1.d0/1.d+1,                     TYPE_FORCE )
type(unitlist),parameter :: unit_36 = unitlist( 'rydberg/angstrom    ',    13.605698075d0,                 TYPE_FORCE )
type(unitlist),parameter :: unit_37 = unitlist( 'rydberg/bohr        ',    13.605698075d0/0.5291772480d0,  TYPE_FORCE )
type(unitlist),parameter :: unit_38 = unitlist( 'rydberg/nm          ',    13.605698075d0/1.d+1,           TYPE_FORCE )
type(unitlist),parameter :: unit_39 = unitlist( 'hartree/angstrom    ',    27.21139615d0,                  TYPE_FORCE )
type(unitlist),parameter :: unit_40 = unitlist( 'hartree/bohr        ',    27.21139615d0/0.5291772480d0,   TYPE_FORCE )
type(unitlist),parameter :: unit_41 = unitlist( 'hartree/nm          ',    27.21139615d0/1.d+1,            TYPE_FORCE )
type(unitlist),parameter :: unit_42 = unitlist( 'ev/angstrom3        ',    1.d0,                                TYPE_PRESSURE )
type(unitlist),parameter :: unit_43 = unitlist( 'ev/bohr3            ',    1.d0/(0.5291772480d0**3),            TYPE_PRESSURE )
type(unitlist),parameter :: unit_44 = unitlist( 'ev/nm3              ',    1.d0/(1.d+1**3),                     TYPE_PRESSURE )
type(unitlist),parameter :: unit_45 = unitlist( 'rydberg/angstrom3   ',    13.605698075d0,                      TYPE_PRESSURE )
type(unitlist),parameter :: unit_46 = unitlist( 'rydberg/bohr3       ',    13.605698075d0/(0.5291772480d0**3),  TYPE_PRESSURE )
type(unitlist),parameter :: unit_47 = unitlist( 'rydberg/nm3         ',    13.605698075d0/(1.d+1**3),           TYPE_PRESSURE )
type(unitlist),parameter :: unit_48 = unitlist( 'hartree/angstrom3   ',    27.21139615d0,                       TYPE_PRESSURE )
type(unitlist),parameter :: unit_49 = unitlist( 'hartree/bohr3       ',    27.21139615d0/(0.5291772480d0**3),   TYPE_PRESSURE )
type(unitlist),parameter :: unit_50 = unitlist( 'hartree/nm3         ',    27.21139615d0/(1.d+1**3),            TYPE_PRESSURE )

type(unitlist), dimension(unit_list_size), parameter ::  unit_list = (/ &
     & unit_1,  unit_2,  unit_3,  unit_4,  unit_5,  unit_6,  unit_7,  unit_8,  unit_9,  unit_10, &
     & unit_11, unit_12, unit_13, unit_14, unit_15, unit_16, unit_17, unit_18, unit_19, unit_20, &
     & unit_21, unit_22, unit_23, unit_24, unit_25, unit_26, unit_27, unit_28, unit_29, unit_30, &
     & unit_31, unit_32, unit_33, unit_34, unit_35, unit_36, unit_37, unit_38, unit_39, unit_40, &
     & unit_41, unit_42, unit_43, unit_44, unit_45, unit_46, unit_47, unit_48, unit_49, unit_50  &
     & /)

contains

integer function get_unit_id( unit, idno, unittype )
    implicit none
    character(*), intent(in) :: unit
    integer, intent(out) :: idno
    integer, intent(out) :: unittype
    integer i
!!$    character(UNITNAMELEN) :: wkunit

!!$    wkunit = unit
    do i = 1, unit_list_size
        if( unit == unit_list(i)%name ) then
            idno = i
            unittype = unit_list(i)%type
            get_unit_id = 0
            return
        end if
    end do 

    if( i > unit_list_size ) then
        get_unit_id = -1
        return
    end if

end function get_unit_id

integer function unit_conv_byname( valin, valout, unitin, unitout )
    implicit none
    real(kind=DP), intent(in) :: valin
    real(kind=DP), intent(out) :: valout
    character(*), intent(in) :: unitin
    character(*), intent(in) :: unitout
    integer i, fgin, fgout, typein, typeout
    real(kind=DP) factorin, factorout

    fgin = 0
    fgout = 0
    do i = 1, unit_list_size
        if( unitin == unit_list(i)%name ) then
            factorin = unit_list(i)%factor
            typein = unit_list(i)%type
            fgin = 1
        else if( unitout == unit_list(i)%name ) then
            factorout = unit_list(i)%factor
            typeout = unit_list(i)%type
            fgout = 1
        end if
        if( fgin == 1 .and. fgout == 1 ) then
            exit
        end if
    end do 

    if( i > unit_list_size ) then
        if( fgin == 0 ) then
            unit_conv_byname = -1
        else
            unit_conv_byname = -2
        end if
        return
    end if

    if( typein /= typeout ) then
        unit_conv_byname = -3
        return
    end if

    if( typein == TYPE_TEMPERATURE ) then
        valout = valin + factorin - factorout
    else
        valout = valin*factorin/factorout
    end if

    unit_conv_byname = 0
    return
end function unit_conv_byname

integer function unit_conv_byid( valin, valout, idin, idout )
    implicit none
    real(DP), intent(in) :: valin
    real(DP), intent(out) :: valout
    integer, intent(in) :: idin
    integer, intent(in) :: idout
    integer i, fgin, fgout, typein, typeout
    real(DP) factorin, factorout

    factorin = unit_list(idin)%factor
    typein = unit_list(idin)%type

    factorout = unit_list(idout)%factor
    typeout = unit_list(idout)%type

    if( typein /= typeout ) then
        unit_conv_byid = -3
        return
    end if

    if( typein == TYPE_TEMPERATURE ) then
        valout = valin + factorin - factorout
    else
        valout = valin*factorin/factorout
    end if

    unit_conv_byid = 0
    return
end function unit_conv_byid


end module m_Const_Parameters
